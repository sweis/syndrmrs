//! Concatenated code combining Reed-Solomon and Duplicated Reed-Muller codes.
//!
//! Spec: HQC 2025-08-22, §3.4 (Concatenated code).
//!
//! The concatenated code works as follows:
//! - Encode: message → RS encode → DupRM encode → concatenated codeword
//! - Decode: concatenated codeword → DupRM decode → RS decode → message
//!
//! The encoding produces n1 * n2 bits (stored as DupRmEncoded), where:
//! - n1 is the RS block length (46, 56, or 90 for HQC-1/3/5)
//! - n2 is the RM block length (384 or 640 bits = 3*128 or 5*128)

use crate::{
    ParameterSet,
    param::{DupRmEncoded, Msg},
    rm::DupRM,
    rs::RS,
};

pub struct Concat;

impl Concat {
    /// Encodes a message using the concatenated RS + DupRM code.
    ///
    /// # Flow
    /// 1. RS encode: k-byte message → n1-byte RS codeword
    /// 2. DupRM encode: n1-byte RS codeword → n1 * multiplicity * 16 bytes
    ///
    /// # Reference
    /// §3.4: "The encoding is the composition of the RS encoding and the DupRM encoding."
    pub fn encode<P: ParameterSet>(msg: &Msg<P>) -> DupRmEncoded<P> {
        // Step 1: RS encode the message
        let rs_codeword = RS::encode::<P>(msg);

        // Step 2: DupRM encode the RS codeword
        DupRM::encode::<P>(&rs_codeword)
    }

    /// Decodes a concatenated codeword back to the original message.
    ///
    /// # Flow
    /// 1. DupRM decode: DupRmEncoded → n1-byte RS codeword (with error correction)
    /// 2. RS decode: n1-byte RS codeword → k-byte message (with error correction)
    ///
    /// # Reference
    /// §3.4: "The decoding is the composition of the DupRM decoding and the RS decoding."
    pub fn decode<P: ParameterSet>(encoded: &DupRmEncoded<P>) -> Msg<P> {
        // Step 1: DupRM decode to get RS codeword
        let rs_codeword = DupRM::decode::<P>(encoded);

        // Step 2: RS decode to get original message
        RS::decode::<P>(&rs_codeword)
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
    fn concat_encode_decode_roundtrip_hqc1() {
        roundtrip_test::<Hqc1Params>();
    }

    #[test]
    fn concat_encode_decode_roundtrip_hqc3() {
        roundtrip_test::<Hqc3Params>();
    }

    #[test]
    fn concat_encode_decode_roundtrip_hqc5() {
        roundtrip_test::<Hqc5Params>();
    }

    fn roundtrip_test<P: ParameterSet>() {
        // Create a test message
        let mut msg = Msg::<P>::default();
        for (i, byte) in msg.as_mut_slice().iter_mut().enumerate() {
            *byte = (i * 17 + 42) as u8; // Deterministic pattern
        }

        // Encode and decode
        let encoded = Concat::encode::<P>(&msg);
        let decoded = Concat::decode::<P>(&encoded);

        assert_eq!(msg, decoded);
    }

    #[test]
    fn concat_encode_matches_intermediates_hqc1() {
        encode_matches_intermediates::<Hqc1Params>("hqc-1");
    }

    #[test]
    fn concat_encode_matches_intermediates_hqc3() {
        encode_matches_intermediates::<Hqc3Params>("hqc-3");
    }

    #[test]
    fn concat_encode_matches_intermediates_hqc5() {
        encode_matches_intermediates::<Hqc5Params>("hqc-5");
    }

    fn encode_matches_intermediates<P: ParameterSet>(variant: &str) {
        let contents = inter_kats(variant);

        // Get the message (m) and expected concatenated codeword
        let msg_bytes = extract_hex_field(contents, "m");
        let expected_concat = extract_hex_field(contents, "Concatenated code word");

        // Convert message bytes to Msg<P>
        let mut msg = Msg::<P>::default();
        msg.as_mut_slice().copy_from_slice(&msg_bytes);

        // Encode using our implementation
        let encoded = Concat::encode::<P>(&msg);

        // Convert encoded to flat bytes for comparison
        let mut encoded_bytes = std::vec![0u8; expected_concat.len()];
        let mut offset = 0;
        for dup_codeword in encoded.as_slice() {
            for rm_codeword in dup_codeword.as_slice() {
                encoded_bytes[offset..offset + 16].copy_from_slice(rm_codeword.as_slice());
                offset += 16;
            }
        }

        assert_eq!(encoded_bytes.as_slice(), expected_concat.as_slice());
    }
}
