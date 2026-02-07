#![no_std]
#![doc = include_str!("../README.md")]
#![warn(missing_docs)]

//!
//! # Security Warning
//!
//! **DO NOT USE THIS LIBRARY IN PRODUCTION.**
//!
//! This is an educational implementation for learning and experimentation.
//! It has not been audited, may contain timing side-channels, and provides
//! no security guarantees.
//!
//! # Usage
//!
//! ```
//! use syndrmrs::kem::{HqcKem, KemCiphertext, KemDecapsulationKey, KemEncapsulationKey, KemSeed, Msg};
//! use syndrmrs::{hqc3::Hqc3Params, ParameterSet};
//!
//! // Obviously don't do this irl
//! let mut seed = KemSeed::from([0x42u8; 32]);
//! let mut msg = Msg::<Hqc3Params>::from([0x42u8; 24]);
//! let salt = [0x42u8; 16];
//!
//! let (ek, dk) = HqcKem::keygen::<Hqc3Params>(&seed);
//! let (ss_sendt, ct) = HqcKem::encaps_deterministic::<Hqc3Params>(&ek, &msg, &salt);
//! let ss_recvt = HqcKem::decaps::<Hqc3Params>(&dk, &ct);
//! assert_eq!(ss_sendt, ss_recvt);
//! ```

#[cfg(test)]
mod test_util;

mod precomputed;

/// §3.1
mod xof;

/// §3.4.3
mod rm;

/// p.18
mod gf256;

/// p.20
mod fft;

/// §3.4.2
mod rs;

/// §3.2 and §3.3
mod vect;

/// §3.4
mod concat;

/// §3.5
mod pke;

/// §3.6
pub mod kem;

/// §4
mod param;

#[cfg(test)]
extern crate std;

pub use param::ParameterSet;

/// HQC-1 parameter set (NIST Security Level 1)
pub mod hqc1 {
    use super::{ParameterSet, param};
    use hybrid_array::Array;
    // Import all needed types explicitly from sizes (our custom types)
    use hybrid_array::sizes::{
        U3, U15, U16, U46, U66, U75, U277, U384, U554, U2208, U2209, U2241, U2321, U4433,
    };

    /// Reed-Solomon generator polynomial coefficients for HQC-1
    /// RS-S1: n=46, k=16, delta=15, g(x) has degree 2*delta=30
    static RS_POLY: param::RSPoly<Hqc1Params> = Array([
        89, 69, 153, 116, 176, 117, 111, 75, 73, 233, 242, 233, 65, 210, 21, 139, 103, 173, 67,
        118, 105, 210, 174, 110, 74, 69, 228, 82, 255, 181, 1,
    ]);

    /// HQC-1 parameter set implementation
    #[derive(Default, Clone, Debug, PartialEq)]
    pub struct Hqc1Params;

    impl ParameterSet for Hqc1Params {
        const HQC_N: u32 = 17669;
        const HQC_OMEGA: u32 = 66;
        const HQC_OMEGA_R: u32 = 75;

        type LittleN = U2209; // U17669;
        type LittleNWords = U277; // ceil(17669 / 64)
        type LittleNWords2 = U554; // 2 * LittleNWords
        type LittleN1 = U46;
        type LittleN2 = U384;
        type LittleN1N2 = U2208; // U17664;
        type LittleK = U16;
        type LittleOmega = U66;
        type LittleOmegaR = U75;

        type RSLittleK = U16;
        type RSDelta = U15;
        type RMMultiplicity = U3;

        type KemEkSize = U2241; // 32 + 2209
        type KemDkSize = U2321; // 2241 + 32 + 16 + 32
        type KemCtSize = U4433; // 2209 + 2208 + 16

        fn rs_poly() -> &'static param::RSPoly<Self> {
            &RS_POLY
        }
        fn alpha_ij_pow() -> &'static param::AlphaIJPow<Self> {
            &crate::precomputed::HQC_1_ALPHA_IJ_POW
        }
    }
}

/// HQC-3 parameter set (NIST Security Level 3)
pub mod hqc3 {
    use super::{ParameterSet, param};
    use hybrid_array::Array;
    // Import all needed types explicitly from sizes (our custom types)
    use hybrid_array::sizes::{
        U5, U16, U24, U56, U100, U114, U561, U640, U1122, U4480, U4482, U4514, U4602, U8978,
    };

    /// Reed-Solomon generator polynomial coefficients for HQC-3
    /// RS-S2: n=56, k=24, delta=16, g(x) has degree 2*delta=32
    static RS_POLY: param::RSPoly<Hqc3Params> = Array([
        45, 216, 239, 24, 253, 104, 27, 40, 107, 50, 163, 210, 227, 134, 224, 158, 119, 13, 158, 1,
        238, 164, 82, 43, 15, 232, 246, 142, 50, 189, 29, 232, 1,
    ]);

    /// HQC-3 parameter set implementation
    #[derive(Default, Clone, Debug, PartialEq)]
    pub struct Hqc3Params;
    impl ParameterSet for Hqc3Params {
        const HQC_N: u32 = 35851;
        const HQC_OMEGA: u32 = 100;
        const HQC_OMEGA_R: u32 = 114;

        type LittleN = U4482; // U35851;
        type LittleNWords = U561; // ceil(35851 / 64)
        type LittleNWords2 = U1122; // 2 * LittleNWords
        type LittleN1 = U56;
        type LittleN2 = U640;
        type LittleN1N2 = U4480; // U35840;
        type LittleK = U24;
        type LittleOmega = U100;
        type LittleOmegaR = U114;

        type RSLittleK = U24;
        type RSDelta = U16;
        type RMMultiplicity = U5;

        type KemEkSize = U4514; // 32 + 4482
        type KemDkSize = U4602; // 4514 + 32 + 24 + 32
        type KemCtSize = U8978; // 4482 + 4480 + 16

        fn rs_poly() -> &'static param::RSPoly<Self> {
            &RS_POLY
        }
        fn alpha_ij_pow() -> &'static param::AlphaIJPow<Self> {
            &crate::precomputed::HQC_3_ALPHA_IJ_POW
        }
    }
}

/// HQC-5 parameter set (NIST Security Level 5)
pub mod hqc5 {
    use super::{ParameterSet, param};
    use hybrid_array::Array;
    // Import all needed types explicitly from sizes (our custom types)
    use hybrid_array::sizes::{
        U5, U29, U32, U90, U131, U149, U640, U901, U1802, U7200, U7205, U7237, U7333, U14421,
    };

    /// Reed-Solomon generator polynomial coefficients for HQC-5
    /// RS-S3: n=90, k=32, delta=29, g(x) has degree 2*delta=58
    static RS_POLY: param::RSPoly<Hqc5Params> = Array([
        49, 167, 49, 39, 200, 121, 124, 91, 240, 63, 148, 71, 150, 123, 87, 101, 32, 215, 159, 71,
        201, 115, 97, 210, 186, 183, 141, 217, 123, 12, 31, 243, 180, 219, 152, 239, 99, 141, 4,
        246, 191, 144, 8, 232, 47, 27, 141, 178, 130, 64, 124, 47, 39, 188, 216, 48, 199, 187, 1,
    ]);

    /// HQC-5 parameter set implementation
    #[derive(Default, Clone, Debug, PartialEq)]
    pub struct Hqc5Params;
    impl ParameterSet for Hqc5Params {
        const HQC_N: u32 = 57637;
        const HQC_OMEGA: u32 = 131;
        const HQC_OMEGA_R: u32 = 149;

        type LittleN = U7205; // U57637;
        type LittleNWords = U901; // ceil(57637 / 64)
        type LittleNWords2 = U1802; // 2 * LittleNWords
        type LittleN1 = U90;
        type LittleN2 = U640;
        type LittleN1N2 = U7200; // U57600;
        type LittleK = U32;
        type LittleOmega = U131;
        type LittleOmegaR = U149;

        type RSLittleK = U32;
        type RSDelta = U29;
        type RMMultiplicity = U5;

        type KemEkSize = U7237; // 32 + 7205
        type KemDkSize = U7333; // 7237 + 32 + 32 + 32
        type KemCtSize = U14421; // 7205 + 7200 + 16

        fn rs_poly() -> &'static param::RSPoly<Self> {
            &RS_POLY
        }

        fn alpha_ij_pow() -> &'static param::AlphaIJPow<Self> {
            &crate::precomputed::HQC_5_ALPHA_IJ_POW
        }
    }
}
