//! Hash and XOF utilities
//!
//! Spec: HQC 2025-08-22, §3.1 (XOF and Hash functions), Table 1 (domain separation).
//! Ref: https://gitlab.com/pqc-hqc/hqc/-/blob/d622142a50f3ce6b6e1f5b15a5119d96c67194e0/src/common/symmetric.c,
//!      https://gitlab.com/pqc-hqc/hqc/-/blob/d622142a50f3ce6b6e1f5b15a5119d96c67194e0/src/common/symmetric.h.

use hybrid_array::{Array, typenum::U32};
use sha3::digest::{Digest, ExtendableOutput, XofReader};
use sha3::{Sha3_256, Sha3_512, Shake256, Shake256Reader};

pub(crate) const XOF_DOMAIN_SEP: u8 = 1;
pub(crate) const G_DOMAIN_SEP: u8 = 0;
pub(crate) const H_DOMAIN_SEP: u8 = 1;
pub(crate) const I_DOMAIN_SEP: u8 = 2;
pub(crate) const J_DOMAIN_SEP: u8 = 3;

pub struct Xof {
    reader: Shake256Reader,
}

impl Xof {
    pub fn init(seed: &[u8], domain_sep: u8) -> Self {
        let mut hasher = Shake256::default();
        sha3::digest::Update::update(&mut hasher, seed);
        sha3::digest::Update::update(&mut hasher, &[domain_sep]);
        let reader = hasher.finalize_xof();
        Self { reader }
    }

    /// Squeeze bytes from the XOF without alignment constraints.
    ///
    /// **WARNING**: Use `squeeze_aligned()` instead when implementing HQC vector sampling!
    /// The reference implementation reads XOF bytes in 8-byte aligned chunks, and using
    /// this method directly will cause KAT test failures due to XOF state desynchronization.
    pub fn squeeze(&mut self, out: &mut [u8]) {
        self.reader.read(out);
    }

    /// Squeeze bytes with 8-byte chunk alignment for compatibility with reference implementation.
    ///
    /// **CRITICAL**: This method MUST be used for HQC vector sampling to match the reference
    /// implementation's behavior. The reference code reads XOF output in 8-byte aligned chunks,
    /// and any non-aligned reads will cause the XOF state to diverge, resulting in different
    /// sampled vectors and KAT test failures.
    ///
    /// # Why this matters
    ///
    /// When sampling fixed-weight vectors, we need buffers of size `3 * HQC_OMEGA` bytes
    /// (e.g., 198 bytes for HQC-1). Since 198 is not a multiple of 8, naive byte-by-byte
    /// or full-buffer reads will consume bytes differently than the reference implementation,
    /// causing the XOF internal state to desynchronize.
    ///
    /// The reference implementation (symmetric.c) always reads in 8-byte chunks:
    /// - For the main portion (floor(len/8) * 8 bytes), read directly
    /// - For the remainder (len % 8 bytes), read 8 bytes into a temp buffer and copy only what's needed
    ///
    /// This ensures that the XOF state advances in consistent 8-byte increments, maintaining
    /// compatibility across implementations.
    ///
    /// # Reference
    ///
    /// See `symmetric.c` in the HQC reference implementation:
    /// ```c
    /// void shake_prng(uint8_t *output, size_t outlen, uint8_t *seed, size_t seedlen) {
    ///     // ... initialization ...
    ///     shake256_inc_squeeze(output, outlen - (outlen % 8), &state);
    ///     if (outlen % 8) {
    ///         uint8_t tmp[8];
    ///         shake256_inc_squeeze(tmp, 8, &state);
    ///         memcpy(output + outlen - (outlen % 8), tmp, outlen % 8);
    ///     }
    /// }
    /// ```
    pub fn squeeze_aligned(&mut self, out: &mut [u8]) {
        const CHUNK: usize = 8;
        let rem = out.len() % CHUNK;
        let main_len = out.len() - rem;

        if main_len > 0 {
            self.reader.read(&mut out[..main_len]);
        }

        if rem != 0 {
            let mut tmp = [0u8; CHUNK];
            self.reader.read(&mut tmp);
            out[main_len..].copy_from_slice(&tmp[..rem]);
        }
    }
}

pub fn hash_g(hash_ek: &[u8], m: &[u8], salt: &[u8]) -> (Array<u8, U32>, Array<u8, U32>) {
    let mut hasher = Sha3_512::new();
    hasher.update(hash_ek);
    hasher.update(m);
    hasher.update(salt);
    hasher.update([G_DOMAIN_SEP]);
    let out = hasher.finalize();
    let mut a = Array::<u8, U32>::default();
    let mut b = Array::<u8, U32>::default();

    a.copy_from_slice(&out[..32]);
    b.copy_from_slice(&out[32..]);
    (a, b)
}

pub fn hash_h(ek: &[u8]) -> Array<u8, U32> {
    let mut hasher = Sha3_256::new();
    hasher.update(ek);
    hasher.update([H_DOMAIN_SEP]);
    let out = hasher.finalize();
    let mut a = Array::<u8, U32>::default();
    a.copy_from_slice(&out[..32]);
    a
}

pub fn hash_i(seed: &[u8]) -> (Array<u8, U32>, Array<u8, U32>) {
    let mut hasher = Sha3_512::new();
    hasher.update(seed);
    hasher.update([I_DOMAIN_SEP]);
    let out = hasher.finalize();
    let mut a = Array::<u8, U32>::default();
    let mut b = Array::<u8, U32>::default();

    a.copy_from_slice(&out[..32]);
    b.copy_from_slice(&out[32..]);
    (a, b)
}

pub fn hash_j(hash_ek: &[u8], sigma: &[u8], c_kem: &[u8]) -> Array<u8, U32> {
    let mut hasher = Sha3_256::new();
    hasher.update(hash_ek);
    hasher.update(sigma);
    hasher.update(c_kem);
    hasher.update([J_DOMAIN_SEP]);
    let out = hasher.finalize();
    let mut a = Array::<u8, U32>::default();
    a.copy_from_slice(&out[..32]);
    a
}

#[cfg(test)]
mod tests {
    use super::{XOF_DOMAIN_SEP, Xof, hash_g, hash_h};
    use crate::test_util::{extract_hex_field, inter_kats};
    use std::vec;

    fn assert_xof_seed_outputs(variant: &str) {
        let contents = inter_kats(variant);

        let seed_kem = extract_hex_field(contents, "seed_kem");
        let seed_pke = extract_hex_field(contents, "seed_pke");
        let sigma = extract_hex_field(contents, "sigma");

        let mut xof = Xof::init(&seed_kem, XOF_DOMAIN_SEP);
        let mut out = vec![0u8; seed_pke.len() + sigma.len()];
        xof.squeeze(&mut out);

        assert_eq!(&out[..seed_pke.len()], seed_pke.as_slice());
        assert_eq!(&out[seed_pke.len()..], sigma.as_slice());
    }

    #[test]
    fn xof_derives_seed_pke_sigma_hqc1() {
        assert_xof_seed_outputs("hqc-1");
    }

    #[test]
    fn xof_derives_seed_pke_sigma_hqc3() {
        assert_xof_seed_outputs("hqc-3");
    }

    #[test]
    fn xof_derives_seed_pke_sigma_hqc5() {
        assert_xof_seed_outputs("hqc-5");
    }

    fn assert_hashes_from_intermediates(variant: &str) {
        let contents = inter_kats(variant);

        let ek_kem = extract_hex_field(contents, "ek_kem");
        let m = extract_hex_field(contents, "m");
        let salt = extract_hex_field(contents, "salt");
        let h_expected = extract_hex_field(contents, "H(ek_kem)");
        let k_expected = extract_hex_field(contents, "secret1");
        let theta_expected = extract_hex_field(contents, "theta");

        let h = hash_h(&ek_kem);
        assert_eq!(h_expected, h.as_slice());

        let (k, theta) = hash_g(&h, &m, &salt);
        assert_eq!(k_expected, k.as_slice());
        assert_eq!(theta_expected, theta.as_slice());
    }

    #[test]
    fn hashes_match_intermediates_hqc1() {
        assert_hashes_from_intermediates("hqc-1");
    }

    #[test]
    fn hashes_match_intermediates_hqc3() {
        assert_hashes_from_intermediates("hqc-3");
    }

    #[test]
    fn hashes_match_intermediates_hqc5() {
        assert_hashes_from_intermediates("hqc-5");
    }
}
