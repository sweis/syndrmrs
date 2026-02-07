//! HQC-KEM orchestration.
//!
//! Implements KEM.Keygen / Encaps / Decaps per HQC spec.
//! - Uses FO-style transform as specified (hashes + reencrypt-and-compare).
//! - Uses constant-time comparison and selection from the `subtle` crate.
//! - Treats malformed ciphertexts carefully: no panics, constant-time behavior.
//!
//! Spec: HQC 2025-08-22, §3.6 (HQC-KEM), algorithms on pp. 26–28.
//! Ref: https://gitlab.com/pqc-hqc/hqc/-/blob/d622142a50f3ce6b6e1f5b15a5119d96c67194e0/src/common/kem.c (crypto_kem_keypair/enc/dec).

use hybrid_array::{Array, typenum::U32};
use subtle::{ConditionallySelectable, ConstantTimeEq};

use crate::{
    ParameterSet,
    param::{HqcSeed, KemCtBytes, KemDkBytes, KemEkBytes, Sigma},
    pke::{HqcPke, PkeCiphertext, PkeDecryptionKey, PkeEncryptionKey},
    vect::{TruncatedVect, Vect},
    xof::{XOF_DOMAIN_SEP, Xof, hash_g, hash_h, hash_j},
};

/// Message type for KEM encapsulation (k bytes, where k is the security parameter)
pub use crate::param::Msg;

/// Salt size in bytes (16 bytes = 128 bits)
const SALT_BYTES: usize = 16;

/// KEM ciphertext: (u, v, salt)
#[derive(Clone, Debug, PartialEq)]
pub struct KemCiphertext<P>
where
    P: ParameterSet,
{
    pub(crate) c_pke: PkeCiphertext<P>,
    pub(crate) salt: [u8; SALT_BYTES],
}

/// KEM encapsulation key (public key)
#[derive(Clone, Debug, PartialEq)]
pub struct KemEncapsulationKey<P>
where
    P: ParameterSet,
{
    pub(crate) ek_pke: PkeEncryptionKey<P>,
}

/// KEM decapsulation key (secret key)
///
/// Contains:
/// - ek_pke: The PKE encryption key (for re-encryption check)
/// - dk_pke: The PKE decryption key
/// - sigma: Random value for implicit rejection (k bytes)
/// - seed_kem: The original KEM seed (for serialization)
#[derive(Clone, Debug)]
pub struct KemDecapsulationKey<P>
where
    P: ParameterSet,
{
    pub(crate) ek_pke: PkeEncryptionKey<P>,
    pub(crate) dk_pke: PkeDecryptionKey<P>,
    pub(crate) sigma: Sigma<P>,
    pub(crate) seed_kem: HqcSeed,
}

/// KEM seed (32 bytes)
pub type KemSeed = HqcSeed;

/// KEM shared secret (32 bytes)
pub type KemSharedSecret = Array<u8, U32>;

/// HQC-KEM as defined in §3.6
pub struct HqcKem;

impl HqcKem {
    /// KEM.KeyGen: Generate encapsulation and decapsulation keys from seed_kem.
    ///
    /// # Algorithm
    /// 1. (seed_pke, sigma) = XOF(seed_kem)
    /// 2. (ek_pke, dk_pke) = PKE.KeyGen(seed_pke)
    /// 3. ek_kem = ek_pke
    /// 4. dk_kem = (ek_pke, dk_pke, sigma, seed_kem)
    #[must_use]
    pub fn keygen<P: ParameterSet>(
        seed_kem: &KemSeed,
    ) -> (KemEncapsulationKey<P>, KemDecapsulationKey<P>) {
        // Step 1: Derive seed_pke and sigma from seed_kem
        let mut xof = Xof::init(seed_kem.as_slice(), XOF_DOMAIN_SEP);
        let mut seed_pke = HqcSeed::default();
        let mut sigma = Sigma::<P>::default();
        xof.squeeze(seed_pke.as_mut_slice());
        xof.squeeze(sigma.as_mut_slice());

        // Step 2: Generate PKE keypair
        let (ek_pke, dk_pke) = HqcPke::keygen::<P>(&seed_pke);

        // Step 3 & 4: Construct KEM keys
        let ek = KemEncapsulationKey {
            ek_pke: ek_pke.clone(),
        };
        let dk = KemDecapsulationKey {
            ek_pke,
            dk_pke,
            sigma,
            seed_kem: *seed_kem,
        };

        (ek, dk)
    }

    /// KEM.Encaps: Encapsulate a shared secret using the encapsulation key.
    ///
    /// # Algorithm
    /// 1. Sample random m and salt
    /// 2. hash_ek = H(ek_kem)
    /// 3. (K, theta) = G(hash_ek, m, salt)
    /// 4. c_pke = PKE.Encrypt(ek_pke, m, theta)
    /// 5. c_kem = (c_pke, salt)
    /// 6. Return (K, c_kem)
    #[must_use]
    pub fn encaps_deterministic<P: ParameterSet>(
        ek: &KemEncapsulationKey<P>,
        m: &Msg<P>,
        salt: &[u8; SALT_BYTES],
    ) -> (KemSharedSecret, KemCiphertext<P>) {
        // Step 2: Compute hash of encryption key
        let ek_bytes = ek.to_bytes();
        let hash_ek = hash_h(&ek_bytes);

        // Step 3: Derive K and theta
        let (k, theta) = hash_g(&hash_ek, m.as_slice(), salt);

        // Step 4: Encrypt message
        let mut theta_seed = HqcSeed::default();
        theta_seed.copy_from_slice(&theta);
        let c_pke = HqcPke::encrypt::<P>(&ek.ek_pke, m, &theta_seed);

        // Step 5: Construct ciphertext
        let c_kem = KemCiphertext { c_pke, salt: *salt };

        (k, c_kem)
    }

    /// KEM.Decaps: Decapsulate a shared secret using the decapsulation key.
    ///
    /// # Algorithm (FO transform with implicit rejection)
    /// 1. m' = PKE.Decrypt(dk_pke, c_pke)
    /// 2. hash_ek = H(ek_pke)
    /// 3. (K', theta') = G(hash_ek, m', salt)
    /// 4. c_pke' = PKE.Encrypt(ek_pke, m', theta')
    /// 5. K_bar = J(hash_ek, sigma, c_kem)
    /// 6. If c_pke == c_pke' then return K' else return K_bar (constant-time)
    #[must_use]
    pub fn decaps<P: ParameterSet>(
        dk: &KemDecapsulationKey<P>,
        c_kem: &KemCiphertext<P>,
    ) -> KemSharedSecret {
        // Step 1: Decrypt to get m'
        let m_prime = HqcPke::decrypt::<P>(&dk.dk_pke, &c_kem.c_pke);

        // Step 2: Compute hash of encryption key
        let ek = KemEncapsulationKey {
            ek_pke: dk.ek_pke.clone(),
        };
        let ek_bytes = ek.to_bytes();
        let hash_ek = hash_h(&ek_bytes);

        // Step 3: Derive K' and theta'
        let (k_prime, theta_prime) = hash_g(&hash_ek, m_prime.as_slice(), &c_kem.salt);

        // Step 4: Re-encrypt to get c_pke'
        let mut theta_seed = HqcSeed::default();
        theta_seed.copy_from_slice(&theta_prime);
        let c_pke_prime = HqcPke::encrypt::<P>(&dk.ek_pke, &m_prime, &theta_seed);

        // Step 5: Compute rejection key K_bar
        let c_kem_bytes = c_kem.to_bytes();
        let k_bar = hash_j(&hash_ek, &dk.sigma, &c_kem_bytes);

        // Step 6: Constant-time comparison and selection
        let u_eq = c_kem.c_pke.u.as_slice().ct_eq(c_pke_prime.u.as_slice());
        let v_eq = c_kem.c_pke.v.as_slice().ct_eq(c_pke_prime.v.as_slice());
        let valid = u_eq & v_eq;

        // Constant-time select: if valid then K' else K_bar
        let mut result = Array::<u8, U32>::default();
        for i in 0..32 {
            result[i] = u8::conditional_select(&k_bar[i], &k_prime[i], valid);
        }

        result
    }
}

impl<P: ParameterSet> KemEncapsulationKey<P> {
    /// Serialize the encapsulation key to bytes.
    /// Format: seed_ek || s
    pub fn to_bytes(&self) -> KemEkBytes<P> {
        let mut bytes = KemEkBytes::<P>::default();
        bytes[..32].copy_from_slice(self.ek_pke.seed_ek.as_slice());
        bytes[32..].copy_from_slice(self.ek_pke.s.as_slice());
        bytes
    }
}

impl<P: ParameterSet> KemDecapsulationKey<P> {
    /// Serialize the decapsulation key to bytes.
    /// Format: ek || seed_dk || sigma || seed_kem
    pub fn to_bytes(&self) -> KemDkBytes<P> {
        let mut bytes = KemDkBytes::<P>::default();
        let mut offset = 0;

        // Append KEM encaps key (ek)
        let ek_bytes = self.ek_pke.to_bytes();
        let ek_len = ek_bytes.len();
        bytes[offset..offset + ek_len].copy_from_slice(ek_bytes.as_slice());
        offset += ek_len;

        // Append PKE decryption key (seed_dk)
        bytes[offset..offset + 32].copy_from_slice(self.dk_pke.seed_dk.as_slice());
        offset += 32;

        // Append sigma
        let sigma_len = self.sigma.len();
        bytes[offset..offset + sigma_len].copy_from_slice(self.sigma.as_slice());
        offset += sigma_len;

        // Append seed_kem
        bytes[offset..offset + 32].copy_from_slice(self.seed_kem.as_slice());

        bytes
    }
}

impl<P: ParameterSet> KemCiphertext<P> {
    /// Serialize the ciphertext to bytes.
    /// Format: u || v || salt
    pub fn to_bytes(&self) -> KemCtBytes<P> {
        let mut bytes = KemCtBytes::<P>::default();
        let u_len = self.c_pke.u.len();
        let v_len = self.c_pke.v.len();
        bytes[..u_len].copy_from_slice(self.c_pke.u.as_slice());
        bytes[u_len..u_len + v_len].copy_from_slice(self.c_pke.v.as_slice());
        bytes[u_len + v_len..u_len + v_len + SALT_BYTES].copy_from_slice(&self.salt);
        bytes
    }

    /// Deserialize the ciphertext from bytes.
    /// Format: u || v || salt
    #[must_use]
    pub fn from_bytes(bytes: &[u8]) -> Self {
        let mut u = Vect::<P>::default();
        let mut v = TruncatedVect::<P>::default();
        let mut salt = [0u8; SALT_BYTES];

        let u_len = u.len();
        let v_len = v.len();

        u.as_mut_slice().copy_from_slice(&bytes[..u_len]);
        v.as_mut_slice()
            .copy_from_slice(&bytes[u_len..u_len + v_len]);
        salt.copy_from_slice(&bytes[u_len + v_len..u_len + v_len + SALT_BYTES]);

        KemCiphertext {
            c_pke: PkeCiphertext { u, v },
            salt,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hqc1::Hqc1Params;
    use crate::hqc3::Hqc3Params;
    use crate::hqc5::Hqc5Params;
    use crate::test_util::{extract_hex_field, inter_kats};

    extern crate std;

    #[test]
    fn kem_keygen_matches_intermediates_hqc1() {
        assert_kem_keygen::<Hqc1Params>("hqc-1");
    }

    #[test]
    fn kem_keygen_matches_intermediates_hqc3() {
        assert_kem_keygen::<Hqc3Params>("hqc-3");
    }

    #[test]
    fn kem_keygen_matches_intermediates_hqc5() {
        assert_kem_keygen::<Hqc5Params>("hqc-5");
    }

    fn assert_kem_keygen<P: ParameterSet>(variant: &str) {
        let contents = inter_kats(variant);

        let seed_kem = extract_hex_field(contents, "seed_kem");
        let seed_pke_expected = extract_hex_field(contents, "seed_pke");
        let sigma_expected = extract_hex_field(contents, "sigma");
        let ek_kem_expected = extract_hex_field(contents, "ek_kem");

        let mut seed = KemSeed::default();
        seed.copy_from_slice(&seed_kem);

        let (ek, dk) = HqcKem::keygen::<P>(&seed);

        // Verify sigma
        assert_eq!(
            sigma_expected.as_slice(),
            dk.sigma.as_slice(),
            "sigma mismatch"
        );

        // Verify seed_pke by checking XOF output
        let mut xof = Xof::init(seed.as_slice(), XOF_DOMAIN_SEP);
        let mut seed_pke = HqcSeed::default();
        xof.squeeze(seed_pke.as_mut_slice());
        assert_eq!(
            seed_pke_expected.as_slice(),
            seed_pke.as_slice(),
            "seed_pke mismatch"
        );

        // Verify ek_kem
        let ek_bytes = ek.to_bytes();
        assert_eq!(
            ek_kem_expected.as_slice(),
            ek_bytes.as_slice(),
            "ek_kem mismatch"
        );
    }

    #[test]
    fn kem_encaps_matches_intermediates_hqc1() {
        assert_kem_encaps::<Hqc1Params>("hqc-1");
    }

    #[test]
    fn kem_encaps_matches_intermediates_hqc3() {
        assert_kem_encaps::<Hqc3Params>("hqc-3");
    }

    #[test]
    fn kem_encaps_matches_intermediates_hqc5() {
        assert_kem_encaps::<Hqc5Params>("hqc-5");
    }

    fn assert_kem_encaps<P: ParameterSet>(variant: &str) {
        let contents = inter_kats(variant);

        // Get inputs
        let ek_kem_bytes = extract_hex_field(contents, "ek_kem");
        let m_bytes = extract_hex_field(contents, "m");
        let salt_bytes = extract_hex_field(contents, "salt");

        // Get expected outputs
        let k_expected = extract_hex_field(contents, "secret1");
        let c_kem_expected = extract_hex_field(contents, "c_kem");

        // Construct encapsulation key from bytes
        let mut seed_ek = HqcSeed::default();
        seed_ek.copy_from_slice(&ek_kem_bytes[..32]);
        let mut s = Vect::<P>::default();
        s.as_mut_slice().copy_from_slice(&ek_kem_bytes[32..]);

        let ek = KemEncapsulationKey {
            ek_pke: PkeEncryptionKey { seed_ek, s },
        };

        // Construct message and salt
        let mut m = Msg::<P>::default();
        m.as_mut_slice().copy_from_slice(&m_bytes);

        let mut salt = [0u8; SALT_BYTES];
        salt.copy_from_slice(&salt_bytes);

        // Encapsulate
        let (k, c_kem) = HqcKem::encaps_deterministic::<P>(&ek, &m, &salt);

        assert_eq!(k_expected.as_slice(), k.as_slice(), "K mismatch");

        let c_kem_bytes = c_kem.to_bytes();
        assert_eq!(
            c_kem_expected.as_slice(),
            c_kem_bytes.as_slice(),
            "c_kem mismatch"
        );
    }

    #[test]
    fn kem_decaps_matches_intermediates_hqc1() {
        assert_kem_decaps::<Hqc1Params>("hqc-1");
    }

    #[test]
    fn kem_decaps_matches_intermediates_hqc3() {
        assert_kem_decaps::<Hqc3Params>("hqc-3");
    }

    #[test]
    fn kem_decaps_matches_intermediates_hqc5() {
        assert_kem_decaps::<Hqc5Params>("hqc-5");
    }

    fn assert_kem_decaps<P: ParameterSet>(variant: &str) {
        let contents = inter_kats(variant);

        // First generate the keys from seed_kem
        let seed_kem = extract_hex_field(contents, "seed_kem");
        let mut seed = KemSeed::default();
        seed.copy_from_slice(&seed_kem);
        let (_, dk) = HqcKem::keygen::<P>(&seed);

        // Get ciphertext
        let c_kem_bytes = extract_hex_field(contents, "c_kem");

        // Parse ciphertext
        let u_len = Vect::<P>::default().len();
        let v_len = TruncatedVect::<P>::default().len();

        let mut u = Vect::<P>::default();
        let mut v = TruncatedVect::<P>::default();
        let mut salt = [0u8; SALT_BYTES];

        u.as_mut_slice().copy_from_slice(&c_kem_bytes[..u_len]);
        v.as_mut_slice()
            .copy_from_slice(&c_kem_bytes[u_len..u_len + v_len]);
        salt.copy_from_slice(&c_kem_bytes[u_len + v_len..]);

        let c_kem = KemCiphertext {
            c_pke: PkeCiphertext { u, v },
            salt,
        };

        // Decapsulate
        let k = HqcKem::decaps::<P>(&dk, &c_kem);

        // Get expected K (should match since ciphertext is valid)
        let k_expected = extract_hex_field(contents, "secret1");
        assert_eq!(k_expected.as_slice(), k.as_slice(), "K mismatch");
    }

    #[test]
    fn kem_roundtrip_hqc1() {
        roundtrip_test::<Hqc1Params>();
    }

    #[test]
    fn kem_roundtrip_hqc3() {
        roundtrip_test::<Hqc3Params>();
    }

    #[test]
    fn kem_roundtrip_hqc5() {
        roundtrip_test::<Hqc5Params>();
    }

    fn roundtrip_test<P: ParameterSet>() {
        // Create test seed
        let mut seed_kem = KemSeed::default();
        for (i, byte) in seed_kem.as_mut_slice().iter_mut().enumerate() {
            *byte = (i * 7 + 13) as u8;
        }

        // Generate keys
        let (ek, dk) = HqcKem::keygen::<P>(&seed_kem);

        // Create test message and salt
        let mut m = Msg::<P>::default();
        for (i, byte) in m.as_mut_slice().iter_mut().enumerate() {
            *byte = (i * 17 + 42) as u8;
        }

        let mut salt = [0u8; SALT_BYTES];
        for (i, byte) in salt.iter_mut().enumerate() {
            *byte = (i * 23 + 99) as u8;
        }

        // Encapsulate
        let (k_encaps, c_kem) = HqcKem::encaps_deterministic::<P>(&ek, &m, &salt);

        // Decapsulate
        let k_decaps = HqcKem::decaps::<P>(&dk, &c_kem);

        assert_eq!(k_encaps.as_slice(), k_decaps.as_slice());
    }

    #[test]
    fn kem_implicit_rejection_hqc1() {
        implicit_rejection_test::<Hqc1Params>();
    }

    fn implicit_rejection_test<P: ParameterSet>() {
        // Create test seed
        let mut seed_kem = KemSeed::default();
        for (i, byte) in seed_kem.as_mut_slice().iter_mut().enumerate() {
            *byte = (i * 7 + 13) as u8;
        }

        // Generate keys
        let (ek, dk) = HqcKem::keygen::<P>(&seed_kem);

        // Create test message and salt
        let mut m = Msg::<P>::default();
        for (i, byte) in m.as_mut_slice().iter_mut().enumerate() {
            *byte = (i * 17 + 42) as u8;
        }

        let mut salt = [0u8; SALT_BYTES];
        for (i, byte) in salt.iter_mut().enumerate() {
            *byte = (i * 23 + 99) as u8;
        }

        // Encapsulate
        let (k_encaps, mut c_kem) = HqcKem::encaps_deterministic::<P>(&ek, &m, &salt);

        // Corrupt the ciphertext
        c_kem.c_pke.u[0] ^= 0x01;

        // Decapsulate - should get K_bar, not K
        let k_decaps = HqcKem::decaps::<P>(&dk, &c_kem);

        // The decapsulated key should NOT match the encapsulated key
        assert_ne!(k_encaps.as_slice(), k_decaps.as_slice());
    }
}
