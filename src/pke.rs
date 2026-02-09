//! HQC-PKE (public-key encryption) core.
//!
//! Implements PKE.Keygen / Encrypt / Decrypt per spec.
//! - Keygen outputs (ek_pke, dk_pke)
//! - Encrypt uses fixed-weight sampling + convolution + concatenated codeword
//! - Decrypt uses syndrome-ish reconstruction + concatenated decoding
//!
//! Spec: HQC 2025-08-22, §3.5 (HQC-PKE), algorithms on pp. 22–25.
//! Ref: https://gitlab.com/pqc-hqc/hqc/-/blob/d622142a50f3ce6b6e1f5b15a5119d96c67194e0/src/ref/hqc.c (hqc_pke_keygen/encrypt/decrypt).

use core::marker::PhantomData;

use crate::{
    ParameterSet,
    concat::Concat,
    param::{DupRmEncoded, HqcSeed, KemEkBytes, Msg},
    vect::{
        TruncatedVect, Vect, sample_fixed_weight_vect_biased, sample_fixed_weight_vect_uniform,
        sample_vect, vect_add, vect_mul,
    },
    xof::{XOF_DOMAIN_SEP, Xof, hash_i},
};

/// PKE encryption key: seed_ek || s
///
/// The encryption key consists of:
/// - seed_ek (32 bytes): Used to derive h deterministically
/// - s (n bytes): Public vector s = h·y + x
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PkeEncryptionKey<P>
where
    P: ParameterSet,
{
    pub(crate) seed_ek: HqcSeed,
    pub(crate) s: Vect<P>,
}

impl<P: ParameterSet> PkeEncryptionKey<P> {
    /// Serialize the encryption key to bytes.
    /// Format: seed_ek || s
    pub fn to_bytes(&self) -> KemEkBytes<P> {
        let mut bytes = KemEkBytes::<P>::default();
        bytes[..32].copy_from_slice(self.seed_ek.as_slice());
        bytes[32..].copy_from_slice(self.s.as_slice());
        bytes
    }
}

/// PKE decryption key: seed_dk
///
/// The decryption key is just the seed used to derive y.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PkeDecryptionKey<P>
where
    P: ParameterSet,
{
    pub(crate) seed_dk: HqcSeed,
    _marker: PhantomData<P>,
}

/// PKE ciphertext: (u, v) where u is n bits and v is n1*n2 bits
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PkeCiphertext<P>
where
    P: ParameterSet,
{
    pub(crate) u: Vect<P>,
    pub(crate) v: TruncatedVect<P>,
}

pub type PkeMessage<P> = Msg<P>;
pub type PkeRandomness = HqcSeed;
pub type PkeSeed = HqcSeed;

/// HQC-PKE as defined in §3.5
pub struct HqcPke;

impl HqcPke {
    /// PKE.KeyGen: Generate encryption and decryption keys from seed_pke.
    ///
    /// # Algorithm
    /// 1. (seed_dk, seed_ek) = I(seed_pke)
    /// 2. Sample y, x from XOF(seed_dk) with fixed weight ω
    /// 3. Sample h from XOF(seed_ek) uniformly
    /// 4. Compute s = h·y + x
    /// 5. ek_pke = (seed_ek, s), dk_pke = seed_dk
    pub fn keygen<P: ParameterSet>(
        seed_pke: &PkeSeed,
    ) -> (PkeEncryptionKey<P>, PkeDecryptionKey<P>) {
        // Step 1: Derive seed_dk and seed_ek using hash function I
        let (seed_dk, seed_ek) = hash_i(seed_pke.as_slice());

        // Step 2: Sample y and x from XOF(seed_dk) with fixed weight ω
        let mut dk_xof = Xof::init(seed_dk.as_slice(), XOF_DOMAIN_SEP);
        let y = sample_fixed_weight_vect_uniform::<P>(&mut dk_xof);
        let x = sample_fixed_weight_vect_uniform::<P>(&mut dk_xof);

        // Step 3: Sample h from XOF(seed_ek) uniformly
        let mut ek_xof = Xof::init(seed_ek.as_slice(), XOF_DOMAIN_SEP);
        let h = sample_vect::<P>(&mut ek_xof);

        // Step 4: Compute s = h·y + x
        let hy = vect_mul::<P>(&h, &y);
        let s = vect_add::<P>(&hy, &x);

        // Step 5: Return keys
        let ek = PkeEncryptionKey { seed_ek, s };
        let dk = PkeDecryptionKey {
            seed_dk,
            _marker: PhantomData,
        };

        (ek, dk)
    }

    /// PKE.Encrypt: Encrypt a message using the encryption key.
    ///
    /// # Algorithm
    /// 1. Sample r2, e, r1 from XOF(theta) with fixed weight ω_r
    /// 2. Derive h from seed_ek
    /// 3. Compute u = h·r2 + r1
    /// 4. Compute v = Truncate(s·r2 + e) + Encode(m)
    pub fn encrypt<P: ParameterSet>(
        ek: &PkeEncryptionKey<P>,
        m: &PkeMessage<P>,
        theta: &PkeRandomness,
    ) -> PkeCiphertext<P> {
        // Step 1: Sample r2, e, r1 from XOF(theta) using biased (fast) sampling
        let mut xof = Xof::init(theta.as_slice(), XOF_DOMAIN_SEP);
        let r2 = sample_fixed_weight_vect_biased::<P>(&mut xof);
        let e = sample_fixed_weight_vect_biased::<P>(&mut xof);
        let r1 = sample_fixed_weight_vect_biased::<P>(&mut xof);

        // Step 2: Derive h from seed_ek
        let mut ek_xof = Xof::init(ek.seed_ek.as_slice(), XOF_DOMAIN_SEP);
        let h = sample_vect::<P>(&mut ek_xof);

        // Step 3: Compute u = h·r2 + r1
        let hr2 = vect_mul::<P>(&h, &r2);
        let u = vect_add::<P>(&hr2, &r1);

        // Step 4: Compute v
        // First: s·r2 + e
        let sr2 = vect_mul::<P>(&ek.s, &r2);
        let sr2_plus_e = vect_add::<P>(&sr2, &e);

        // Truncate to n1*n2 bits
        let truncated = truncate::<P>(&sr2_plus_e);

        // Encode message using concatenated code
        let m_encoded = Concat::encode::<P>(m);
        let m_encoded_bytes = dup_rm_encoded_to_bytes::<P>(&m_encoded);

        // v = truncated + m_encoded
        let mut v = TruncatedVect::<P>::default();
        for i in 0..v.len() {
            v[i] = truncated[i] ^ m_encoded_bytes[i];
        }

        PkeCiphertext { u, v }
    }

    /// PKE.Decrypt: Decrypt a ciphertext using the decryption key.
    ///
    /// # Algorithm
    /// 1. Derive y from seed_dk
    /// 2. Compute v - Truncate(u·y)
    /// 3. Decode using concatenated code to recover message
    pub fn decrypt<P: ParameterSet>(
        dk: &PkeDecryptionKey<P>,
        ct: &PkeCiphertext<P>,
    ) -> PkeMessage<P> {
        // Step 1: Derive y from seed_dk
        let mut xof = Xof::init(dk.seed_dk.as_slice(), XOF_DOMAIN_SEP);
        let y = sample_fixed_weight_vect_uniform::<P>(&mut xof);

        // Step 2: Compute u·y and truncate
        let uy = vect_mul::<P>(&ct.u, &y);
        let truncated = truncate::<P>(&uy);

        // Step 3: Compute v - truncated (XOR in GF(2))
        let mut v_minus = TruncatedVect::<P>::default();
        for i in 0..v_minus.len() {
            v_minus[i] = ct.v[i] ^ truncated[i];
        }

        // Step 4: Decode using concatenated code
        let encoded = bytes_to_dup_rm_encoded::<P>(&v_minus);
        Concat::decode::<P>(&encoded)
    }
}

/// Truncate a vector from n bits to n1*n2 bits.
fn truncate<P: ParameterSet>(v: &Vect<P>) -> TruncatedVect<P> {
    let mut out = TruncatedVect::<P>::default();
    let n1n2_bytes = out.len();

    // Copy the first n1*n2 bits
    out.as_mut_slice()[..n1n2_bytes].copy_from_slice(&v.as_slice()[..n1n2_bytes]);

    // Mask the last byte if n1*n2 is not a multiple of 8
    // For HQC, n1*n2 is always a multiple of 8 (384*46 = 17664, etc.)
    // but let's be safe
    let n1n2_bits = n1n2_bytes * 8; // Assuming byte-aligned
    let rem = n1n2_bits % 8;
    if rem != 0
        && let Some(last) = out.as_mut_slice().last_mut()
    {
        *last &= (1u8 << rem) - 1;
    }

    out
}

/// Convert DupRmEncoded to a flat byte array.
fn dup_rm_encoded_to_bytes<P: ParameterSet>(encoded: &DupRmEncoded<P>) -> TruncatedVect<P> {
    let mut out = TruncatedVect::<P>::default();
    let mut offset = 0;

    for dup_codeword in encoded.as_slice() {
        for rm_codeword in dup_codeword.as_slice() {
            out.as_mut_slice()[offset..offset + 16].copy_from_slice(rm_codeword.as_slice());
            offset += 16;
        }
    }

    out
}

/// Convert a flat byte array to DupRmEncoded.
fn bytes_to_dup_rm_encoded<P: ParameterSet>(bytes: &TruncatedVect<P>) -> DupRmEncoded<P> {
    let mut encoded = DupRmEncoded::<P>::default();
    let mut offset = 0;

    for dup_codeword in encoded.as_mut_slice() {
        for rm_codeword in dup_codeword.as_mut_slice() {
            rm_codeword
                .as_mut_slice()
                .copy_from_slice(&bytes.as_slice()[offset..offset + 16]);
            offset += 16;
        }
    }

    encoded
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
    fn pke_keygen_matches_intermediates_hqc1() {
        assert_pke_keygen::<Hqc1Params>("hqc-1");
    }

    #[test]
    fn pke_keygen_matches_intermediates_hqc3() {
        assert_pke_keygen::<Hqc3Params>("hqc-3");
    }

    #[test]
    fn pke_keygen_matches_intermediates_hqc5() {
        assert_pke_keygen::<Hqc5Params>("hqc-5");
    }

    fn assert_pke_keygen<P: ParameterSet>(variant: &str) {
        let contents = inter_kats(variant);

        let seed_pke = extract_hex_field(contents, "seed_pke");
        let seed_dk_expected = extract_hex_field(contents, "seed_dk");
        let seed_ek_expected = extract_hex_field(contents, "seed_ek");
        let x_expected = extract_hex_field(contents, "x");
        let y_expected = extract_hex_field(contents, "y");
        let s_expected = extract_hex_field(contents, "s");

        let mut seed = PkeSeed::default();
        seed.copy_from_slice(&seed_pke);

        let (ek, dk) = HqcPke::keygen::<P>(&seed);

        assert_eq!(
            seed_dk_expected.as_slice(),
            dk.seed_dk.as_slice(),
            "seed_dk mismatch"
        );
        assert_eq!(
            seed_ek_expected.as_slice(),
            ek.seed_ek.as_slice(),
            "seed_ek mismatch"
        );
        assert_eq!(s_expected.as_slice(), ek.s.as_slice(), "s mismatch");

        // Verify y and x by re-deriving them
        let mut dk_xof = Xof::init(dk.seed_dk.as_slice(), XOF_DOMAIN_SEP);
        let y = sample_fixed_weight_vect_uniform::<P>(&mut dk_xof);
        let x = sample_fixed_weight_vect_uniform::<P>(&mut dk_xof);

        assert_eq!(y_expected.as_slice(), y.as_slice(), "y mismatch");
        assert_eq!(x_expected.as_slice(), x.as_slice(), "x mismatch");
    }

    #[test]
    fn pke_encrypt_matches_intermediates_hqc1() {
        assert_pke_encrypt::<Hqc1Params>("hqc-1");
    }

    #[test]
    fn pke_encrypt_matches_intermediates_hqc3() {
        assert_pke_encrypt::<Hqc3Params>("hqc-3");
    }

    #[test]
    fn pke_encrypt_matches_intermediates_hqc5() {
        assert_pke_encrypt::<Hqc5Params>("hqc-5");
    }

    fn assert_pke_encrypt<P: ParameterSet>(variant: &str) {
        let contents = inter_kats(variant);

        // Get encryption inputs
        let seed_ek = extract_hex_field(contents, "seed_ek");
        let s = extract_hex_field(contents, "s");
        let theta = extract_hex_field(contents, "theta");
        let m = extract_hex_field(contents, "m");

        // Get expected outputs
        let u_expected = extract_hex_field(contents, "c_pke->u");
        let v_expected = extract_hex_field(contents, "c_pke->v");

        // Construct encryption key
        let mut ek = PkeEncryptionKey::<P> {
            seed_ek: HqcSeed::default(),
            s: Vect::<P>::default(),
        };
        ek.seed_ek.copy_from_slice(&seed_ek);
        ek.s.as_mut_slice().copy_from_slice(&s);

        // Construct message and theta
        let mut msg = PkeMessage::<P>::default();
        msg.as_mut_slice().copy_from_slice(&m);

        let mut theta_arr = PkeRandomness::default();
        theta_arr.copy_from_slice(&theta);

        // Encrypt
        let ct = HqcPke::encrypt::<P>(&ek, &msg, &theta_arr);

        assert_eq!(u_expected.as_slice(), ct.u.as_slice(), "u mismatch");
        assert_eq!(v_expected.as_slice(), ct.v.as_slice(), "v mismatch");
    }

    #[test]
    fn pke_decrypt_matches_intermediates_hqc1() {
        assert_pke_decrypt::<Hqc1Params>("hqc-1");
    }

    #[test]
    fn pke_decrypt_matches_intermediates_hqc3() {
        assert_pke_decrypt::<Hqc3Params>("hqc-3");
    }

    #[test]
    fn pke_decrypt_matches_intermediates_hqc5() {
        assert_pke_decrypt::<Hqc5Params>("hqc-5");
    }

    fn assert_pke_decrypt<P: ParameterSet>(variant: &str) {
        let contents = inter_kats(variant);

        // Get decryption inputs
        let dk_pke = extract_hex_field(contents, "dk_pke");
        let u = extract_hex_field(contents, "c_pke.u");
        let v = extract_hex_field(contents, "c_pke.v");

        // Get expected output
        let m_prime_expected = extract_hex_field(contents, "m_prime");

        // Construct decryption key
        let mut dk = PkeDecryptionKey::<P> {
            seed_dk: HqcSeed::default(),
            _marker: PhantomData,
        };
        dk.seed_dk.copy_from_slice(&dk_pke);

        // Construct ciphertext
        let mut ct = PkeCiphertext::<P> {
            u: Vect::<P>::default(),
            v: TruncatedVect::<P>::default(),
        };
        ct.u.as_mut_slice().copy_from_slice(&u);
        ct.v.as_mut_slice().copy_from_slice(&v);

        // Decrypt
        let m_prime = HqcPke::decrypt::<P>(&dk, &ct);

        assert_eq!(
            m_prime_expected.as_slice(),
            m_prime.as_slice(),
            "m_prime mismatch"
        );
    }

    #[test]
    fn pke_roundtrip_hqc1() {
        roundtrip_test::<Hqc1Params>();
    }

    #[test]
    fn pke_roundtrip_hqc3() {
        roundtrip_test::<Hqc3Params>();
    }

    #[test]
    fn pke_roundtrip_hqc5() {
        roundtrip_test::<Hqc5Params>();
    }

    fn roundtrip_test<P: ParameterSet>() {
        // Create a test seed for keygen
        let mut seed_pke = PkeSeed::default();
        for (i, byte) in seed_pke.as_mut_slice().iter_mut().enumerate() {
            *byte = (i * 7 + 13) as u8;
        }

        // Generate keys
        let (ek, dk) = HqcPke::keygen::<P>(&seed_pke);

        // Create a test message
        let mut msg = PkeMessage::<P>::default();
        for (i, byte) in msg.as_mut_slice().iter_mut().enumerate() {
            *byte = (i * 17 + 42) as u8;
        }

        // Create test randomness
        let mut theta = PkeRandomness::default();
        for (i, byte) in theta.as_mut_slice().iter_mut().enumerate() {
            *byte = (i * 23 + 99) as u8;
        }

        // Encrypt and decrypt
        let ct = HqcPke::encrypt::<P>(&ek, &msg, &theta);
        let decrypted = HqcPke::decrypt::<P>(&dk, &ct);

        assert_eq!(msg.as_slice(), decrypted.as_slice());
    }

    // ---- Serialization tests ----

    #[test]
    fn pke_ek_to_bytes_hqc1() {
        ek_to_bytes_test::<Hqc1Params>();
    }

    #[test]
    fn pke_ek_to_bytes_hqc3() {
        ek_to_bytes_test::<Hqc3Params>();
    }

    #[test]
    fn pke_ek_to_bytes_hqc5() {
        ek_to_bytes_test::<Hqc5Params>();
    }

    fn ek_to_bytes_test<P: ParameterSet>() {
        let mut seed_pke = PkeSeed::default();
        for (i, byte) in seed_pke.as_mut_slice().iter_mut().enumerate() {
            *byte = (i * 7 + 13) as u8;
        }

        let (ek, _dk) = HqcPke::keygen::<P>(&seed_pke);
        let bytes = ek.to_bytes();

        // First 32 bytes should be seed_ek
        assert_eq!(&bytes[..32], ek.seed_ek.as_slice());
        // Remaining bytes should be s
        assert_eq!(&bytes[32..], ek.s.as_slice());
    }
}
