//! Duplicated RM(1,7) inner code.
//!
//! Spec: HQC 2025-08-22, §3.4.3 (Duplicated Reed-Muller codes), Table 4.
//! Ref: https://gitlab.com/pqc-hqc/hqc/-/blob/d622142a50f3ce6b6e1f5b15a5119d96c67194e0/src/ref/reed_muller.c (encode/hadamard/expand+sum/find_peaks).

use crate::{
    ParameterSet,
    param::{DupRmCodeword, DupRmEncoded, RmCodeword, RsCodeword},
};

pub struct DupRM;
impl DupRM {
    pub fn encode<P: ParameterSet>(msg: &RsCodeword<P>) -> DupRmEncoded<P> {
        let mut out = DupRmEncoded::<P>::default();
        let msg_bytes = msg.as_slice().iter();
        let codewords = out.as_mut_slice();
        for (m, c) in msg_bytes.zip(codewords) {
            *c = Self::encode_byte::<P>(*m);
        }
        out
    }

    pub fn decode<P: ParameterSet>(encoded: &DupRmEncoded<P>) -> RsCodeword<P> {
        let mut msg = RsCodeword::<P>::default();
        let msg_bytes = msg.as_mut_slice().iter_mut();
        let codewords = encoded.as_slice();
        for (m, c) in msg_bytes.zip(codewords) {
            *m = Self::decode_codeword::<P>(c);
        }
        msg
    }

    fn encode_byte<P: ParameterSet>(m: u8) -> DupRmCodeword<P> {
        let mut out = DupRmCodeword::<P>::default();
        let codewords = out.as_mut_slice();
        for c in codewords {
            *c = RM::encode_byte(m);
        }
        out
    }

    /// Decodes the duplicated codeword at once.
    ///
    /// This follows the usual RM(1,7) decoding recipe used in HQC reference code:
    /// 1) "expand+sum": map bits to ±1 and sum duplicates per position
    /// 2) Walsh–Hadamard transform (FHT) to get correlations with all affine functions
    /// 3) pick argmax |corr| with a stable tie-break
    /// 4) recover message:
    ///      - bits 0..6 from the index
    ///      - bit 7 from the sign (negative => add the all-ones vector => constant term 1)
    fn decode_codeword<P: ParameterSet>(codeword: &DupRmCodeword<P>) -> u8 {
        // Each DupRmCodeword<P> is `[RmCodeword; dup]` (as implied by encode_byte),
        // so `codeword.as_slice()` yields the duplicated RM blocks.
        let blocks: &[RmCodeword] = codeword.as_slice();

        // ---- Step 1: expand+sum (0 -> +1, 1 -> -1), sum over duplicates ----
        // scores[i] corresponds to position i in the length-128 RM word.
        // Range: [-dup, +dup], dup is 3 or 5 in HQC.
        let mut scores = [0i16; 128];
        for (i, score) in scores.iter_mut().enumerate() {
            let mut s: i16 = 0;
            for cw in blocks {
                let bit = Self::rm_bit_at(cw, i);
                // bit=0 => +1, bit=1 => -1
                s += if bit == 0 { 1 } else { -1 };
            }
            *score = s;
        }

        // ---- Step 2: Walsh–Hadamard transform in-place ----
        // After this, scores[k] is the correlation with the affine function indexed by k.
        // For RM(1,7), k ranges over 7-bit linear coefficients.
        Self::hadamard_128(&mut scores);

        // ---- Step 3: find peak (max abs), tie-break by smallest index ----
        let mut best_idx: usize = 0;
        let mut best_val: i16 = scores[0];

        for (i, &v) in scores.iter().enumerate().skip(1) {
            let av = v.abs();
            let ab = best_val.abs();
            if av > ab || (av == ab && i < best_idx) {
                best_idx = i;
                best_val = v;
            }
        }

        // ---- Step 4: recover message byte ----
        // Convention: negative correlation means the best match is the complement,
        // i.e., constant term (bit 7) = 1.
        let linear = best_idx as u8; // 7-bit index, bits 0..6
        let const_bit = u8::from(best_val < 0);

        linear | (const_bit << 7)
    }

    /// Extract bit i (0..127) from an RM codeword produced by `encode_byte`.
    #[inline]
    fn rm_bit_at(cw: &RmCodeword, i: usize) -> u8 {
        debug_assert!(i < 128);
        let byte = cw[i >> 3];
        (byte >> (i & 7)) & 1
    }

    /// In-place Walsh–Hadamard transform on 128 signed scores.
    ///
    /// This is the standard FHT:
    /// for len=1,2,4,...:
    ///   (a,b) -> (a+b, a-b)
    #[inline]
    fn hadamard_128(v: &mut [i16; 128]) {
        let mut len = 1usize;
        while len < 128 {
            let step = len << 1;
            let mut i = 0usize;
            while i < 128 {
                for j in 0..len {
                    let a = v[i + j];
                    let b = v[i + j + len];
                    v[i + j] = a + b;
                    v[i + j + len] = a - b;
                }
                i += step;
            }
            len = step;
        }
    }
}

struct RM;
impl RM {
    fn bit0mask(x: u32) -> u32 {
        (0u32).wrapping_sub(x & 1)
    }

    /// Encodes one byte into one RM[1,7] codeword
    pub fn encode_byte(msg: u8) -> RmCodeword {
        let mut first_word: u32 = Self::bit0mask(u32::from(msg >> 7));
        first_word ^= Self::bit0mask(u32::from(msg)) & 0xaaaaaaaa;
        first_word ^= Self::bit0mask(u32::from(msg >> 1)) & 0xcccccccc;
        first_word ^= Self::bit0mask(u32::from(msg >> 2)) & 0xf0f0f0f0;
        first_word ^= Self::bit0mask(u32::from(msg >> 3)) & 0xff00ff00;
        first_word ^= Self::bit0mask(u32::from(msg >> 4)) & 0xffff0000;

        let mut u32s = [0u32; 4];
        u32s[0] = first_word;
        first_word ^= Self::bit0mask(u32::from(msg >> 5));
        u32s[1] = first_word;
        first_word ^= Self::bit0mask(u32::from(msg >> 6));
        u32s[3] = first_word;
        first_word ^= Self::bit0mask(u32::from(msg >> 5));
        u32s[2] = first_word;

        let mut out = RmCodeword::default();
        out[..4].copy_from_slice(&u32s[0].to_le_bytes());
        out[4..8].copy_from_slice(&u32s[1].to_le_bytes());
        out[8..12].copy_from_slice(&u32s[2].to_le_bytes());
        out[12..16].copy_from_slice(&u32s[3].to_le_bytes());
        out
    }

    /// Decodes one RM[1,7] codeword into one byte
    ///
    /// This is a simplified version of DupRM::decode_codeword that works on a single
    /// non-duplicated codeword. It uses the same Walsh-Hadamard decoding approach.
    #[allow(dead_code)] // May be used in future
    pub fn decode_byte(codeword: &RmCodeword) -> u8 {
        // Step 1: expand (0 -> +1, 1 -> -1)
        let mut scores = [0i16; 128];
        for i in 0..128 {
            let byte = codeword[i >> 3];
            let bit = (byte >> (i & 7)) & 1;
            // bit=0 => +1, bit=1 => -1
            scores[i] = if bit == 0 { 1 } else { -1 };
        }

        // Step 2: Walsh–Hadamard transform in-place
        DupRM::hadamard_128(&mut scores);

        // Step 3: find peak (max abs), tie-break by smallest index
        let mut best_idx: usize = 0;
        let mut best_val: i16 = scores[0];

        for (i, &v) in scores.iter().enumerate().skip(1) {
            let av = v.abs();
            let ab = best_val.abs();
            if av > ab || (av == ab && i < best_idx) {
                best_idx = i;
                best_val = v;
            }
        }

        // Step 4: recover message byte
        let linear = best_idx as u8;
        let const_bit = u8::from(best_val < 0);

        linear | (const_bit << 7)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ParameterSet;
    use crate::param::DupRmEncoded;
    use crate::test_util::TestRng;
    use crate::test_util::{extract_hex_field, inter_kats};
    use std::vec;

    #[test]
    fn rm_encoder_is_deterministic() {
        for b in 0u16..=255 {
            let b = b as u8;
            let c1 = RM::encode_byte(b);
            let c2 = RM::encode_byte(b);
            assert_eq!(c1.as_slice(), c2.as_slice(), "nondeterministic encode?");
        }
    }

    fn assert_rm_encode_from_intermediates<P: ParameterSet>(variant: &str) {
        let contents = inter_kats(variant);

        let rs_codeword = extract_hex_field(contents, "Reed-Solomon code word");
        let concatenated = extract_hex_field(contents, "Concatenated code word");

        // Convert the rs_codeword bytes to RsCodeword<P>
        let mut rs = RsCodeword::<P>::default();
        rs.as_mut_slice().copy_from_slice(&rs_codeword);

        // Encode using our RM implementation
        let encoded = DupRM::encode::<P>(&rs);

        // Flatten the nested structure to compare with expected bytes
        let mut encoded_bytes = vec![0u8; concatenated.len()];
        let mut offset = 0;
        for dup_codeword in encoded.as_slice() {
            for rm_codeword in dup_codeword.as_slice() {
                encoded_bytes[offset..offset + 16].copy_from_slice(rm_codeword.as_slice());
                offset += 16;
            }
        }

        assert_eq!(encoded_bytes.as_slice(), concatenated.as_slice());
    }

    fn assert_rm_decode_from_intermediates<P: ParameterSet>(variant: &str) {
        let contents = inter_kats(variant);

        let rs_codeword_expected = extract_hex_field(contents, "Reed-Solomon code word");
        let concatenated = extract_hex_field(contents, "Concatenated code word");

        // Convert the concatenated bytes to DupRmEncoded<P>
        // DupRmEncoded<P> is Array<DupRmCodeword<P>, LittleK>
        // where DupRmCodeword<P> is Array<RmCodeword, RMMultiplicity>
        // and RmCodeword is Array<u8, U16>
        let mut encoded = DupRmEncoded::<P>::default();
        let mut offset = 0;
        for dup_codeword in encoded.as_mut_slice() {
            for rm_codeword in dup_codeword.as_mut_slice() {
                rm_codeword
                    .as_mut_slice()
                    .copy_from_slice(&concatenated[offset..offset + 16]);
                offset += 16;
            }
        }

        // Decode using our RM implementation
        let decoded = DupRM::decode::<P>(&encoded);

        assert_eq!(decoded.as_slice(), rs_codeword_expected.as_slice());
    }

    #[test]
    fn rm_encode_matches_intermediates_hqc1() {
        assert_rm_encode_from_intermediates::<crate::hqc1::Hqc1Params>("hqc-1");
    }

    #[test]
    fn rm_encode_matches_intermediates_hqc3() {
        assert_rm_encode_from_intermediates::<crate::hqc3::Hqc3Params>("hqc-3");
    }

    #[test]
    fn rm_encode_matches_intermediates_hqc5() {
        assert_rm_encode_from_intermediates::<crate::hqc5::Hqc5Params>("hqc-5");
    }

    #[test]
    fn rm_decode_matches_intermediates_hqc1() {
        assert_rm_decode_from_intermediates::<crate::hqc1::Hqc1Params>("hqc-1");
    }

    #[test]
    fn rm_decode_matches_intermediates_hqc3() {
        assert_rm_decode_from_intermediates::<crate::hqc3::Hqc3Params>("hqc-3");
    }

    #[test]
    fn rm_decode_matches_intermediates_hqc5() {
        assert_rm_decode_from_intermediates::<crate::hqc5::Hqc5Params>("hqc-5");
    }

    #[test]
    fn rm_encode_decode_byte_roundtrip() {
        // Test all 256 possible byte values
        for b in 0u16..=255 {
            let b = b as u8;
            let codeword = RM::encode_byte(b);
            let decoded = RM::decode_byte(&codeword);
            assert_eq!(decoded, b, "Failed to roundtrip byte 0x{b:02x}");
        }
    }

    #[test]
    fn dup_rm_encode_decode_roundtrip_hqc1() {
        dup_rm_roundtrip_test::<crate::hqc1::Hqc1Params>();
    }

    #[test]
    fn dup_rm_encode_decode_roundtrip_hqc3() {
        dup_rm_roundtrip_test::<crate::hqc3::Hqc3Params>();
    }

    #[test]
    fn dup_rm_encode_decode_roundtrip_hqc5() {
        dup_rm_roundtrip_test::<crate::hqc5::Hqc5Params>();
    }

    fn dup_rm_roundtrip_test<P: ParameterSet>() {
        let mut rng = TestRng::new();

        // Test with 100 random RS codewords
        for _ in 0..100 {
            let mut rs_codeword = RsCodeword::<P>::default();
            // Fill with random bytes
            for byte in rs_codeword.as_mut_slice() {
                *byte = (rng.next_u32() & 0xff) as u8;
            }

            // Encode and decode
            let encoded = DupRM::encode::<P>(&rs_codeword);
            let decoded = DupRM::decode::<P>(&encoded);

            assert_eq!(
                decoded.as_slice(),
                rs_codeword.as_slice(),
                "Failed to roundtrip RS codeword"
            );
        }
    }

    #[test]
    fn dup_rm_decode_with_errors_hqc1() {
        dup_rm_error_correction_test::<crate::hqc1::Hqc1Params>();
    }

    #[test]
    fn dup_rm_decode_with_errors_hqc3() {
        dup_rm_error_correction_test::<crate::hqc3::Hqc3Params>();
    }

    #[test]
    fn dup_rm_decode_with_errors_hqc5() {
        dup_rm_error_correction_test::<crate::hqc5::Hqc5Params>();
    }

    fn dup_rm_error_correction_test<P: ParameterSet>() {
        let mut rng = TestRng::new();

        // Test that the duplication provides error correction
        // RM(1,7) can correct errors when duplicated, as the decoder uses majority voting
        for _ in 0..20 {
            let mut rs_codeword = RsCodeword::<P>::default();
            for byte in rs_codeword.as_mut_slice() {
                *byte = (rng.next_u32() & 0xff) as u8;
            }

            let mut encoded = DupRM::encode::<P>(&rs_codeword);

            // Introduce a small number of bit flips to the encoded data
            // We flip bits in only one of the duplicates so majority voting should still work
            for dup_codeword in encoded.as_mut_slice().iter_mut().take(5) {
                // Flip a single bit in the first duplicate only
                let rm_codeword = &mut dup_codeword.as_mut_slice()[0];
                rm_codeword[0] ^= 1; // Flip the LSB of first byte
            }

            let decoded = DupRM::decode::<P>(&encoded);

            assert_eq!(
                decoded.as_slice(),
                rs_codeword.as_slice(),
                "Failed to decode with minor errors"
            );
        }
    }
}
