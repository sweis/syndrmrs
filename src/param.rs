use core::fmt::Debug;
use core::ops::{Add, Mul, Sub};

use hybrid_array::typenum::operator_aliases::{Diff, Sum};
use hybrid_array::{
    Array, ArraySize,
    typenum::{Prod, U1, U3, U16, U32, U256},
};

/// HQC parameter set trait defining security level and code parameters
///
/// This trait defines all the compile-time parameters for an HQC variant.
/// Implementations exist for HQC-1, HQC-3, and HQC-5.
pub trait ParameterSet: Default + Clone + Debug + PartialEq + Eq
where
    Self::RSDelta: Add<U1> + Add<Self::RSDelta>,
    <Self::RSDelta as Add<U1>>::Output: ArraySize,
    <Self::RSDelta as Add<Self::RSDelta>>::Output: ArraySize + Add<U1>,
    <<Self::RSDelta as Add<Self::RSDelta>>::Output as Add<U1>>::Output: ArraySize,
    Self::LittleN1: Sub<U1>,
    <Self::LittleN1 as Sub<U1>>::Output: ArraySize,
    Self::LittleOmega: Mul<U3>,
    <Self::LittleOmega as Mul<U3>>::Output: ArraySize,
{
    /// HQC_N: length of the ambient space as a u32, denoted `n` in Table 5
    /// Used for Barrett reduction and sampling
    const HQC_N: u32;

    /// HQC_OMEGA: weight of vectors x and y, denoted `ω` in Table 5
    const HQC_OMEGA: u32;

    /// HQC_OMEGA_R: weight of vectors r1, r2, e, denoted `ω_r` in Table 5
    const HQC_OMEGA_R: u32;

    /// length of the ambient space, denoted `n` in Table 5
    type LittleN: ArraySize;

    /// number of u64 words needed to store n bits: ceil(n / 64)
    type LittleNWords: ArraySize;

    /// number of u64 words for polynomial multiplication result: 2 * LittleNWords
    type LittleNWords2: ArraySize;

    /// block length of shortened Reed-Solomon code, denoted `n1` in Table 5, and denoted `n` in Table 3
    type LittleN1: ArraySize;

    /// length of the Reed-Muller code, denoted `n2` in Table 5
    type LittleN2: ArraySize;

    /// length of the concatenated code, product of `n1` and `n2` in Table 5
    type LittleN1N2: ArraySize;

    /// dimension of the concat code, denoted `k` in §4.1
    type LittleK: ArraySize;

    /// weight of the vectors x and y, denoted `\omega` in Table 5
    type LittleOmega: ArraySize;

    /// weight of the vectors r1, r2, and e, denoted `\omega_r` in Table 5
    type LittleOmegaR: ArraySize;

    /// the k in the RS[n,k,d] code, denoted `k` in Table 3
    /// NOT the same as the `k` in Table 5 (LittleK)
    type RSLittleK: ArraySize;

    /// correcting capacity of Reed-Solomon code, see Table 3
    type RSDelta: ArraySize;

    /// multiplicity of duplicated Reed-Muller code, see Table 4
    type RMMultiplicity: ArraySize;

    /// Reed-Solomon generator polynomial coefficients
    /// Length is 2*delta + 1 (see Table 3)
    fn rs_poly() -> &'static RSPoly<Self>;

    /// Reed-Solomon galois field powers
    /// Returns a 2D array of precomputed powers [2*delta][n1-1]
    fn alpha_ij_pow() -> &'static AlphaIJPow<Self>;

    /// KEM encapsulation key size in bytes: 32 (seed) + n (s vector)
    type KemEkSize: ArraySize;

    /// KEM decapsulation key size in bytes: KemEKSize + 32 (dk_pke.seed_dk) + k + 32 (seed_kem)
    type KemDkSize: ArraySize;

    /// KEM ciphertext size in bytes: n (u) + n1*n2 (v) + 16 (salt)
    type KemCtSize: ArraySize;
}

// Shared params

/// |seed| in page 29
pub type HqcSeed = Array<u8, U32>;

// Derived type aliases computed from base parameters using typenum arithmetic

/// 2*delta, see Table 3
pub type RS2Delta<P> = Sum<<P as ParameterSet>::RSDelta, <P as ParameterSet>::RSDelta>;

/// delta + 1, used for error evaluator polynomial z(x)
pub type RSDeltaPlus1<P> = Sum<<P as ParameterSet>::RSDelta, U1>;

/// n1 - 1 where n1 is the RS block length, see Table 3
pub type RSN1Minus1<P> = Diff<<P as ParameterSet>::LittleN1, U1>;

/// Generator polynomial degree = 2*delta + 1, see Table 3
pub type RSPolyLen<P> = Sum<RS2Delta<P>, U1>;

/// the random message to encrypt, denoted `m` in HQC-KEM.Encaps
pub type Msg<P> = Array<u8, <P as ParameterSet>::LittleK>;

/// Reed-Solomon code word
pub type RsCodeword<P> = Array<u8, <P as ParameterSet>::LittleN1>;

/// Reed-Solomon generator polynomial
/// Length is 2*delta + 1
pub type RSPoly<P> = Array<u8, RSPolyLen<P>>;

/// Reed-Solomon galois field powers
/// 2D array of dimensions [2*delta][n1-1]
pub type AlphaIJPow<P> = Array<Array<u16, RSN1Minus1<P>>, RS2Delta<P>>;

/// Reed-Solomon syndromes
/// Array of 2*delta syndrome values for error detection
pub type RsSyndromes<P> = Array<u16, RS2Delta<P>>;

/// Reed-Solomon error locator polynomial σ(x)
/// For GF(256), we need up to 256 coefficients for FFT operations
/// The actual degree will be ≤ delta, but FFT requires power-of-2 sizing
pub type RsErrorLocatorPoly = Array<u16, U256>;

/// Reed-Solomon error evaluator polynomial z(x)
/// Length is delta + 1
pub type RsErrorEvaluatorPoly<P> = Array<u16, RSDeltaPlus1<P>>;

/// Reed-Solomon error locations
/// Binary array indicating which positions have errors (size 256 for GF(256))
pub type RsErrorLocations = Array<u8, U256>;

/// Reed-Solomon error values
/// The magnitude of errors at each codeword position
pub type RsErrorValues<P> = Array<u16, <P as ParameterSet>::LittleN1>;

/// all RM codewords in HQC are 128 bits, §3.4.1
pub type RmCodeword = Array<u8, U16>;

/// multiplicity is 3 for HQC-1 and 5 for HQC-3 & HQC-5, see Table 4
pub type DupRmCodeword<P> = Array<RmCodeword, <P as ParameterSet>::RMMultiplicity>;

/// message encoded with the duplicated Reed-Muller code
/// This encodes N1 bytes (the RS codeword), not K bytes (the original message)
pub type DupRmEncoded<P> = Array<DupRmCodeword<P>, <P as ParameterSet>::LittleN1>;

/// 3 * omega bytes, used for seed expansion during vector sampling
pub type Omega3<P> = Array<u8, Prod<<P as ParameterSet>::LittleOmega, U3>>;

/// Array of u64 words to store n bits
pub type VectNWords<P> = Array<u64, <P as ParameterSet>::LittleNWords>;

/// Array of u64 words for polynomial multiplication (2 * n_words)
pub type VectNWords2<P> = Array<u64, <P as ParameterSet>::LittleNWords2>;

/// KEM encapsulation key as bytes: seed_ek || s
pub type KemEkBytes<P> = Array<u8, <P as ParameterSet>::KemEkSize>;

/// KEM decapsulation key as bytes: ek || seed_dk || sigma || seed_kem
/// This matches the NIST API where the encapsulation key is appended to the decapsulation key.
pub type KemDkBytes<P> = Array<u8, <P as ParameterSet>::KemDkSize>;

/// KEM ciphertext as bytes: u || v || salt
pub type KemCtBytes<P> = Array<u8, <P as ParameterSet>::KemCtSize>;

/// Sigma value for implicit rejection: k bytes where k is the dimension of the concatenated code
pub type Sigma<P> = Array<u8, <P as ParameterSet>::LittleK>;
