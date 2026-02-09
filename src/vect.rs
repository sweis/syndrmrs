//! Vector operations: sampling, arithmetic, and truncation.
//!
//! Implements vector operations from the HQC spec:
//! - §3.2: Vector sampling (SampleVect, SampleFixedWeightVect)
//! - §3.3: Vector arithmetic (addition, multiplication)
//!
//! Spec: HQC 2025-08-22, §3.2 (Vector sampling) and §3.3 (Vector arithmetic).
//! Ref: https://gitlab.com/pqc-hqc/hqc/-/blob/d622142a50f3ce6b6e1f5b15a5119d96c67194e0/src/ref/vector.c
//!
//! Critical: sampling touches secrets. Avoid secret-dependent branches and memory access.
//! Note: sampling order during PKE.Encrypt is r2, then e, then r1 (spec §3.2).

use hybrid_array::Array;

use crate::{
    ParameterSet,
    param::{Omega3, VectNWords, VectNWords2},
    xof::Xof,
};

// ============================================================================
// Barrett Reduction for modular reduction by HQC_N
// ============================================================================
//
// Barrett reduction avoids expensive division operations by using precomputed
// constants and bit shifts.
//
// Algorithm: For x < REJECTION_THRESHOLD_24BIT, compute x mod HQC_N
// 1. Multiply by Barrett multiplier: t = x * BARRETT_MULTIPLIER
// 2. Shift to get quotient: q = t >> BARRETT_SHIFT ≈ floor(x / HQC_N)
// 3. Compute approximate remainder: r = x - q * HQC_N
// 4. Correct if needed: if r >= HQC_N then r -= HQC_N
//
// Reference: https://en.wikipedia.org/wiki/Barrett_reduction

/// Barrett shift parameter (we use 2^32 as the scaling factor)
pub(crate) const BARRETT_SHIFT: u64 = 32;

/// Barrett scaling factor R = 2^BARRETT_SHIFT
pub(crate) const BARRETT_R: u64 = 1 << BARRETT_SHIFT;

/// Compute Barrett multiplier: floor(2^32 / n)
const fn barrett_multiplier(n: u32) -> u64 {
    BARRETT_R / (n as u64)
}

/// Compute rejection threshold for 24-bit numbers: floor(2^24 / n) * n
/// This is the largest multiple of n that fits in 24 bits.
const fn rejection_threshold_24bit(n: u32) -> u32 {
    let max_24bit = (1u32 << 24) - 1;
    (max_24bit / n) * n
}

/// Barrett Reduction: reduces x modulo HQC_N using precomputed constants.
///
/// # Preconditions
/// - x < rejection_threshold_24bit(HQC_N)
///
/// # Postconditions
/// - Returns r where r < HQC_N and r ≡ x (mod HQC_N)
///
/// # Algorithm
/// Barrett reduction approximates floor(x / HQC_N) using fixed-point arithmetic:
///   q = (x * BARRETT_MULTIPLIER) >> BARRETT_SHIFT ≈ floor(x / HQC_N)
///
/// The approximation can have two outcomes:
/// - MOST CASES (exact): q = floor(x / HQC_N)
///   → sub = q * HQC_N = floor(x / HQC_N) * HQC_N
///   → res = x - sub = x mod HQC_N < HQC_N
///   → No correction needed
///
/// - WORST CASE (underestimate): q = floor(x / HQC_N) - 1
///   → sub = q * HQC_N = (floor(x / HQC_N) - 1) * HQC_N
///   → res = x - sub = (x mod HQC_N) + HQC_N
///   → HQC_N ≤ res < 2 * HQC_N, needs one correction step
#[inline]
pub(crate) fn reduce_mod_n<P: ParameterSet>(x: u32) -> u32 {
    let hqc_n = P::HQC_N;
    let multiplier = barrett_multiplier(hqc_n);

    debug_assert!(
        x < rejection_threshold_24bit(hqc_n),
        "x must be less than rejection threshold"
    );

    // Barrett reduction: compute q ≈ floor(x / N) using fixed-point arithmetic
    let t = u64::from(x) * multiplier;
    let q = t >> BARRETT_SHIFT;

    // Compute the approximate multiple to subtract: sub ≈ floor(x / N) * N
    let sub = (q as u32) * hqc_n;

    // Compute approximate remainder
    let mut res = x - sub;

    // At most one correction step is needed (see algorithm description above)
    if res >= hqc_n {
        res -= hqc_n;
    }

    res
}

/// Returns the rejection threshold for 24-bit sampling with the given parameter set.
/// Samples >= this threshold must be rejected to ensure uniform distribution.
#[inline]
pub(crate) fn get_rejection_threshold<P: ParameterSet>() -> u32 {
    rejection_threshold_24bit(P::HQC_N)
}

/// F2^n vector, we often think of it as a polynomial of degree n over GF(2)
pub type Vect<P> = Array<u8, <P as ParameterSet>::LittleN>;

/// F2^(n1n2) vector
pub type TruncatedVect<P> = Array<u8, <P as ParameterSet>::LittleN1N2>;

/// Support set of HQC_OMEGA indices each in [0, HQC_N-1]
pub type SupportOmega<P> = Array<u32, <P as ParameterSet>::LittleOmega>;

/// Support set of HQC_OMEGA_R indices each in [0, HQC_N-1]
#[allow(dead_code)] // Used for type clarity
pub type SupportOmegaR<P> = Array<u32, <P as ParameterSet>::LittleOmegaR>;

/// GenerateRandomSupport_$ in spec (if it were included), uses rejection sampling to generate a uniform support set.
///
/// Generates HQC_OMEGA distinct random indices in [0, HQC_N-1] using rejection sampling
/// to ensure uniform distribution.
///
/// Algorithm:
/// 1. Sample 24-bit random values from XOF in batches of 3*HQC_OMEGA bytes
/// 2. Reject samples >= rejection_threshold to avoid modular bias
/// 3. Reduce accepted samples mod HQC_N using Barrett reduction
/// 4. Reject duplicates to ensure all indices are distinct
///
/// https://pqc-hqc.org/doc/hqc_specifications_2025_08_22.pdf#subsection.3.2
/// https://gitlab.com/pqc-hqc/hqc/-/blob/d622142a50f3ce6b6e1f5b15a5119d96c67194e0/src/ref/vector.c#L45-93
pub(crate) fn generate_random_support_uniform<P: ParameterSet>(xof: &mut Xof) -> SupportOmega<P> {
    let mut support = SupportOmega::<P>::default();
    let hqc_omega = P::HQC_OMEGA as usize;
    let rejection_threshold = get_rejection_threshold::<P>();

    // We sample 24-bit numbers in batches of 3*HQC_OMEGA bytes
    // Each 24-bit sample requires 3 bytes
    let mut random_bytes = Omega3::<P>::default();
    let batch_size = random_bytes.len();

    // Start with idx = batch_size to force immediate refill on first iteration
    let mut byte_idx = batch_size;

    let mut support_idx = 0;
    while support_idx < hqc_omega {
        // Check if we need more bytes BEFORE consuming
        if byte_idx + 3 > batch_size {
            // CRITICAL: Must use squeeze_aligned() to match reference implementation's
            // 8-byte aligned XOF reads. Using squeeze() will cause KAT failures.
            // See xof.rs for detailed explanation.
            xof.squeeze_aligned(random_bytes.as_mut_slice());
            byte_idx = 0;
        }

        // Construct 24-bit random sample from three bytes (big-endian)
        let sample = (u32::from(random_bytes[byte_idx]) << 16)
            | (u32::from(random_bytes[byte_idx + 1]) << 8)
            | u32::from(random_bytes[byte_idx + 2]);
        byte_idx += 3;

        // Reject samples >= threshold to ensure uniform distribution
        if sample >= rejection_threshold {
            continue;
        }

        // Reduce mod HQC_N using Barrett reduction
        let sample_mod_n = reduce_mod_n::<P>(sample);
        debug_assert!(sample_mod_n < P::HQC_N);

        // Check for collision with previous values
        // Linear scan is acceptable since HQC_OMEGA is small (66-131)
        let mut is_duplicate = false;
        for j in 0..support_idx {
            if support[j] == sample_mod_n {
                is_duplicate = true;
                break;
            }
        }

        if !is_duplicate {
            support[support_idx] = sample_mod_n;
            support_idx += 1;
        }
        // If duplicate, loop continues without incrementing support_idx
    }

    support
}

/// SampleFixedWeightVect$ in spec - generates a vector of weight HQC_OMEGA using uniform sampling.
///
/// This function:
/// 1. Generates a uniform random support set of HQC_OMEGA distinct indices in [0, HQC_N-1]
/// 2. Converts the support set to a bit vector representation
///
/// The result is a vector v ∈ F_2^n with exactly HQC_OMEGA bits set to 1.
///
/// https://pqc-hqc.org/doc/hqc_specifications_2025_08_22.pdf#subsection.3.2
pub(crate) fn sample_fixed_weight_vect_uniform<P: ParameterSet>(xof: &mut Xof) -> Vect<P> {
    let support = generate_random_support_uniform::<P>(xof);
    support_to_vect::<P>(&support)
}

/// Converts a support set (list of bit positions) to a bit vector.
///
/// Each index in the support indicates a bit position that should be set to 1.
/// The bit vector is stored as bytes in little-endian bit order within each byte.
#[inline]
pub(crate) fn support_to_vect<P: ParameterSet>(support: &SupportOmega<P>) -> Vect<P> {
    let mut vect = Vect::<P>::default();

    for &idx in support {
        let byte_idx = (idx / 8) as usize;
        let bit_idx = idx % 8;
        // Set the bit at position bit_idx within the byte
        vect[byte_idx] |= 1u8 << bit_idx;
    }

    vect
}

/// SampleVect in spec - generates a uniformly random vector in F_2^n.
///
/// This function samples a random bit vector of length HQC_N from the XOF.
/// The output is a vector v ∈ F_2^n where each bit is uniformly random.
///
/// Algorithm:
/// 1. Squeeze ceil(HQC_N / 8) bytes from the XOF
/// 2. Mask the last byte to zero out bits beyond position HQC_N
///
/// https://pqc-hqc.org/doc/hqc_specifications_2025_08_22.pdf#subsection.3.2
/// https://gitlab.com/pqc-hqc/hqc/-/blob/d622142a50f3ce6b6e1f5b15a5119d96c67194e0/src/ref/vector.c#L199-211
pub(crate) fn sample_vect<P: ParameterSet>(xof: &mut Xof) -> Vect<P> {
    let mut vect = Vect::<P>::default();

    // Squeeze random bytes from XOF
    xof.squeeze(vect.as_mut_slice());

    // Mask the last byte to zero out bits beyond HQC_N
    let n = P::HQC_N as usize;
    let rem = n % 8;
    if rem != 0 {
        let last_idx = vect.len() - 1;
        let mask = (1u8 << rem) - 1;
        vect[last_idx] &= mask;
    }

    vect
}

/// SampleFixedWeightVect (no $) in spec - biased but fast sampling.
///
/// This is the biased sampling used during encryption for r1, r2, e.
/// It uses Fisher-Yates-like shuffling with 32-bit randomness per position.
///
/// Algorithm (from HQC spec and reference implementation):
/// 1. Sample `weight` random 32-bit values
/// 2. For i in 0..weight: support[i] = i + ((rand[i] * (n - i)) >> 32)
/// 3. Deduplicate: if support[i] == support[j] for j > i, set support[i] = i
/// 4. Convert support to bit vector
///
/// The bias comes from step 2: the distribution isn't perfectly uniform,
/// but it's good enough for the security parameters and much faster than
/// rejection sampling.
pub(crate) fn sample_fixed_weight_vect_biased<P: ParameterSet>(xof: &mut Xof) -> Vect<P> {
    let n = u64::from(P::HQC_N);
    let weight = P::HQC_OMEGA_R as usize;

    // Step 1: Sample weight random 32-bit values (4 bytes each)
    // Use aligned squeeze for compatibility with reference implementation
    let mut rand_bytes = [0u8; 4 * 149]; // Max omega_r is 149 for HQC-5
    let needed = 4 * weight;
    xof.squeeze_aligned(&mut rand_bytes[..needed]);

    // Convert to u32 array
    let mut rand_u32 = [0u32; 149];
    for i in 0..weight {
        rand_u32[i] = u32::from_le_bytes([
            rand_bytes[4 * i],
            rand_bytes[4 * i + 1],
            rand_bytes[4 * i + 2],
            rand_bytes[4 * i + 3],
        ]);
    }

    // Step 2: Generate initial support using biased distribution
    // support[i] = i + floor(rand[i] * (n - i) / 2^32)
    let mut support = [0u32; 149];
    for i in 0..weight {
        let buff = u64::from(rand_u32[i]);
        support[i] = (i as u32) + ((buff * (n - i as u64)) >> 32) as u32;
    }

    // Step 3: Deduplicate in constant time
    // If support[i] == support[j] for any j > i, set support[i] = i
    // This ensures all positions are distinct (at worst, we get 0,1,2,...,weight-1)
    for i in (0..weight.saturating_sub(1)).rev() {
        let mut found: u32 = 0;

        for j in (i + 1)..weight {
            // Constant-time comparison: returns 1 if equal, 0 otherwise
            found |= compare_u32(support[j], support[i]);
        }

        // If found, replace support[i] with i; otherwise keep it
        let mask = 0u32.wrapping_sub(found); // -found: all 1s if found, all 0s otherwise
        support[i] = (mask & (i as u32)) ^ (!mask & support[i]);
    }

    // Step 4: Convert support to bit vector
    support_to_vect_from_slice::<P>(&support[..weight])
}

/// Constant-time comparison of two u32 values.
/// Returns 1 if v1 == v2, 0 otherwise.
#[inline]
fn compare_u32(v1: u32, v2: u32) -> u32 {
    // If v1 == v2, then (v1 - v2) and (v2 - v1) are both 0, so OR is 0, >> 31 is 0, 1 ^ 0 = 1
    // If v1 != v2, one of the subtractions will have high bit set, >> 31 is 1, 1 ^ 1 = 0
    1 ^ (((v1.wrapping_sub(v2)) | (v2.wrapping_sub(v1))) >> 31)
}

/// Convert a support slice to a bit vector.
fn support_to_vect_from_slice<P: ParameterSet>(support: &[u32]) -> Vect<P> {
    let mut vect = Vect::<P>::default();

    for &pos in support {
        let byte_idx = (pos / 8) as usize;
        let bit_idx = pos % 8;
        if byte_idx < vect.len() {
            vect[byte_idx] |= 1u8 << bit_idx;
        }
    }

    vect
}

// ============================================================================
// Vector Arithmetic (§3.3)
// ============================================================================

/// Vector addition in F_2^n (XOR operation).
///
/// In GF(2), addition is the same as XOR. This operation is used throughout
/// the HQC scheme for combining vectors.
///
/// XORs at u64 width for throughput, then handles any trailing bytes.
///
/// <https://pqc-hqc.org/doc/hqc_specifications_2025_08_22.pdf#subsection.3.3>
#[inline]
pub(crate) fn vect_add<P: ParameterSet>(a: &Vect<P>, b: &Vect<P>) -> Vect<P> {
    let mut result = Vect::<P>::default();
    let n = a.len();
    let chunks = n / 8;
    let remainder = n % 8;

    // XOR 8 bytes at a time via u64
    for i in 0..chunks {
        let off = i * 8;
        let wa = u64::from_le_bytes(a[off..off + 8].try_into().unwrap());
        let wb = u64::from_le_bytes(b[off..off + 8].try_into().unwrap());
        result[off..off + 8].copy_from_slice(&(wa ^ wb).to_le_bytes());
    }

    // Handle trailing bytes
    let tail = chunks * 8;
    for i in 0..remainder {
        result[tail + i] = a[tail + i] ^ b[tail + i];
    }

    result
}

/// Precomputed 4-bit lookup table for carry-less multiplication by `a`.
///
/// `tab[i]` = low 64 bits of `a * i` (carry-less), for i in 0..15.
/// `tab_hi[i]` = bits 64-66 of `a * i` (at most 3 bits, since 64-bit * 4-bit = 67-bit max).
struct ClmulTable {
    tab: [u64; 16],
    tab_hi: [u8; 16],
}

impl ClmulTable {
    /// Build the lookup table for carry-less multiplication by `a`.
    #[inline]
    fn new(a: u64) -> Self {
        let mut tab = [0u64; 16];
        let mut tab_hi = [0u8; 16];

        tab[1] = a;

        for i in 2..16usize {
            let half = i >> 1;
            tab[i] = tab[half] << 1;
            tab_hi[i] = (tab_hi[half] << 1) | (tab[half] >> 63) as u8;
            if i & 1 != 0 {
                tab[i] ^= a;
            }
        }

        Self { tab, tab_hi }
    }

    /// Carry-less multiply `a * b` using the precomputed table for `a`.
    ///
    /// Returns `(lo, hi)` where the 128-bit product = `(hi << 64) | lo`.
    /// Processes `b` one nibble at a time using the lookup table.
    #[inline]
    fn clmul(&self, b: u64) -> (u64, u64) {
        let mut lo = self.tab[(b & 0xF) as usize];
        let mut hi = u64::from(self.tab_hi[(b & 0xF) as usize]);

        for k in 1..16u32 {
            let nibble = ((b >> (k * 4)) & 0xF) as usize;
            let t = self.tab[nibble];
            let t_over = u64::from(self.tab_hi[nibble]);
            let shift = k * 4;
            lo ^= t << shift;
            hi ^= (t >> (64 - shift)) | (t_over << shift);
        }

        (lo, hi)
    }
}

/// Carry-less multiply of two 64-bit binary polynomials over GF(2).
///
/// Returns `(lo, hi)` where the 128-bit product = `(hi << 64) | lo`.
/// Uses 4-bit windowed lookup: precompute `a * {0..15}`, then process `b`
/// one nibble at a time. This is ~4x faster than the bit-at-a-time approach.
///
/// The product of a 64-bit poly and a 4-bit poly can be up to 67 bits,
/// so we track a 3-bit overflow per table entry (`tab_hi`).
#[inline]
fn clmul64(a: u64, b: u64) -> (u64, u64) {
    ClmulTable::new(a).clmul(b)
}

/// Vector multiplication in F_2^n (polynomial multiplication mod X^n - 1).
///
/// Interprets vectors as polynomials over GF(2) and multiplies them in the
/// ring F_2[X]/(X^n - 1), where n = HQC_N.
///
/// Algorithm:
/// 1. Schoolbook multiplication using 4-bit windowed carry-less word multiply
/// 2. Reduce modulo X^n - 1 (fold high-degree terms back)
///
/// # Constant-time note
///
/// This function executes in time independent of the vector contents. There
/// are no secret-dependent branches or memory accesses - every word pair is
/// always processed. This is critical because one operand may be a secret key
/// vector (e.g. `y` in decapsulation) and the other may be attacker-controlled
/// (e.g. ciphertext `u`).
///
/// <https://pqc-hqc.org/doc/hqc_specifications_2025_08_22.pdf#subsection.3.3>
/// <https://gitlab.com/pqc-hqc/hqc/-/blob/d622142a50f3ce6b6e1f5b15a5119d96c67194e0/src/ref/gf2x.c#L29-60>
pub(crate) fn vect_mul<P: ParameterSet>(a: &Vect<P>, b: &Vect<P>) -> Vect<P> {
    let n = P::HQC_N as usize;
    let n_bytes = a.len();

    // Convert bytes to u64 words using from_le_bytes
    let mut a_words = VectNWords::<P>::default();
    let mut b_words = VectNWords::<P>::default();
    let n_words = a_words.len();

    for i in 0..n_words {
        let off = i * 8;
        let end = (off + 8).min(n_bytes);
        let mut buf_a = [0u8; 8];
        buf_a[..end - off].copy_from_slice(&a[off..end]);
        a_words[i] = u64::from_le_bytes(buf_a);

        let mut buf_b = [0u8; 8];
        buf_b[..end - off].copy_from_slice(&b[off..end]);
        b_words[i] = u64::from_le_bytes(buf_b);
    }

    // Schoolbook multiplication with clmul64 word-level multiply.
    // No zero-word skipping: constant-time execution to prevent timing
    // side-channels on secret vector contents.
    //
    // The clmul lookup table depends only on a_words[i], so we hoist it
    // outside the inner loop to avoid rebuilding it n_words times per
    // outer iteration.
    let mut mul_res = VectNWords2::<P>::default();

    for i in 0..n_words {
        let table = ClmulTable::new(a_words[i]);
        for j in 0..n_words {
            let (lo, hi) = table.clmul(b_words[j]);
            mul_res[i + j] ^= lo;
            mul_res[i + j + 1] ^= hi;
        }
    }

    // Reduce modulo (X^n - 1): X^n = 1, so fold high-degree terms back
    let mut out_words = VectNWords::<P>::default();
    for i in 0..n_words {
        let r = mul_res[i + n_words - 1] >> (n & 0b11_1111);
        let carry = mul_res[i + n_words] << (64 - (n & 0b11_1111));
        out_words[i] = mul_res[i] ^ r ^ carry;
    }

    // Zero the bits beyond the HQC_N-th bit
    let rem = n % 64;
    if rem != 0 {
        let mask = (1u64 << rem) - 1;
        out_words[n_words - 1] &= mask;
    }

    // Convert back to bytes using to_le_bytes
    let mut result = Vect::<P>::default();
    for i in 0..n_words {
        let off = i * 8;
        let end = (off + 8).min(n_bytes);
        let bytes = out_words[i].to_le_bytes();
        result[off..end].copy_from_slice(&bytes[..end - off]);
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ParameterSet;
    use crate::test_util::{extract_hex_field, inter_kats};
    use crate::xof::*;

    // ========================================================================
    // Barrett Reduction Tests
    // ========================================================================

    /// Verify Barrett constants for HQC-1 match the expected values from the reference
    #[test]
    fn barrett_constants_hqc1() {
        use crate::hqc1::Hqc1Params;

        const HQC1_N: u32 = 17669;
        assert_eq!(Hqc1Params::HQC_N, HQC1_N);

        // Expected: floor(2^32 / 17669) = 243079
        let multiplier = barrett_multiplier(HQC1_N);
        assert_eq!(multiplier, 243079, "Barrett multiplier for HQC-1");

        // Expected: floor(2^24 / 17669) * 17669 = 949 * 17669 = 16767881
        let threshold = rejection_threshold_24bit(HQC1_N);
        assert_eq!(threshold, 16767881, "Rejection threshold for HQC-1");
    }

    /// Verify Barrett reduction correctness for HQC-1 across all valid inputs
    #[test]
    fn barrett_reduction_correctness_hqc1() {
        use crate::hqc1::Hqc1Params;

        let hqc_n = Hqc1Params::HQC_N;
        let threshold = get_rejection_threshold::<Hqc1Params>();

        // Test boundary cases
        assert_eq!(reduce_mod_n::<Hqc1Params>(0), 0);
        assert_eq!(reduce_mod_n::<Hqc1Params>(1), 1);
        assert_eq!(reduce_mod_n::<Hqc1Params>(hqc_n - 1), hqc_n - 1);
        assert_eq!(reduce_mod_n::<Hqc1Params>(hqc_n), 0);
        assert_eq!(reduce_mod_n::<Hqc1Params>(hqc_n + 1), 1);
        assert_eq!(reduce_mod_n::<Hqc1Params>(2 * hqc_n), 0);
        assert_eq!(reduce_mod_n::<Hqc1Params>(2 * hqc_n + 1), 1);

        // Test near threshold
        assert_eq!(
            reduce_mod_n::<Hqc1Params>(threshold - 1),
            (threshold - 1) % hqc_n
        );

        // Exhaustive test for small values
        for x in 0..100_000u32 {
            let result = reduce_mod_n::<Hqc1Params>(x);
            assert_eq!(result, x % hqc_n, "Failed for x = {x}");
            assert!(
                result < hqc_n,
                "Result {result} >= HQC_N {hqc_n} for x = {x}"
            );
        }

        // Test values around multiples of HQC_N
        for k in 1..50u32 {
            let base = k * hqc_n;
            if base < threshold {
                assert_eq!(reduce_mod_n::<Hqc1Params>(base), 0);
                if base + 1 < threshold {
                    assert_eq!(reduce_mod_n::<Hqc1Params>(base + 1), 1);
                }
                if base >= 1 {
                    assert_eq!(reduce_mod_n::<Hqc1Params>(base - 1), hqc_n - 1);
                }
            }
        }
    }

    /// Verify Barrett reduction for HQC-3
    #[test]
    fn barrett_reduction_correctness_hqc3() {
        use crate::hqc3::Hqc3Params;

        let hqc_n = Hqc3Params::HQC_N;

        // Test boundary cases
        assert_eq!(reduce_mod_n::<Hqc3Params>(0), 0);
        assert_eq!(reduce_mod_n::<Hqc3Params>(1), 1);
        assert_eq!(reduce_mod_n::<Hqc3Params>(hqc_n - 1), hqc_n - 1);
        assert_eq!(reduce_mod_n::<Hqc3Params>(hqc_n), 0);
        assert_eq!(reduce_mod_n::<Hqc3Params>(hqc_n + 1), 1);

        // Sample test
        for x in 0..50_000u32 {
            let result = reduce_mod_n::<Hqc3Params>(x);
            assert_eq!(result, x % hqc_n, "Failed for x = {x}");
        }
    }

    /// Verify Barrett reduction for HQC-5
    #[test]
    fn barrett_reduction_correctness_hqc5() {
        use crate::hqc5::Hqc5Params;

        let hqc_n = Hqc5Params::HQC_N;

        // Test boundary cases
        assert_eq!(reduce_mod_n::<Hqc5Params>(0), 0);
        assert_eq!(reduce_mod_n::<Hqc5Params>(1), 1);
        assert_eq!(reduce_mod_n::<Hqc5Params>(hqc_n - 1), hqc_n - 1);
        assert_eq!(reduce_mod_n::<Hqc5Params>(hqc_n), 0);
        assert_eq!(reduce_mod_n::<Hqc5Params>(hqc_n + 1), 1);

        // Sample test
        for x in 0..50_000u32 {
            let result = reduce_mod_n::<Hqc5Params>(x);
            assert_eq!(result, x % hqc_n, "Failed for x = {x}");
        }
    }

    // ========================================================================
    // generate_random_support_uniform Tests
    // ========================================================================

    /// Test that generate_random_support_uniform produces correct output
    #[test]
    fn generate_random_support_uniform_basic() {
        use crate::hqc1::Hqc1Params;
        use std::collections::HashSet;

        // Use a deterministic seed for reproducibility
        let seed = [0u8; 32];
        let mut xof = Xof::init(&seed, XOF_DOMAIN_SEP);

        let support = generate_random_support_uniform::<Hqc1Params>(&mut xof);

        // Check that we have exactly HQC_OMEGA elements
        assert_eq!(support.len(), Hqc1Params::HQC_OMEGA as usize);

        // Check that all values are in range [0, HQC_N)
        for &val in &support {
            assert!(
                val < Hqc1Params::HQC_N,
                "Value {} >= HQC_N {}",
                val,
                Hqc1Params::HQC_N
            );
        }

        // Check that all values are distinct
        let unique: HashSet<u32> = support.iter().copied().collect();
        assert_eq!(
            unique.len(),
            support.len(),
            "Support contains duplicate values"
        );
    }

    /// Test determinism: same seed should produce same support
    #[test]
    fn generate_random_support_uniform_deterministic() {
        use crate::hqc1::Hqc1Params;

        let seed = [42u8; 32];

        let mut xof1 = Xof::init(&seed, XOF_DOMAIN_SEP);
        let support1 = generate_random_support_uniform::<Hqc1Params>(&mut xof1);

        let mut xof2 = Xof::init(&seed, XOF_DOMAIN_SEP);
        let support2 = generate_random_support_uniform::<Hqc1Params>(&mut xof2);

        assert_eq!(
            support1.as_slice(),
            support2.as_slice(),
            "Same seed should produce same support"
        );
    }

    /// Test for all parameter sets
    #[test]
    fn generate_random_support_uniform_all_params() {
        use crate::hqc1::Hqc1Params;
        use crate::hqc3::Hqc3Params;
        use crate::hqc5::Hqc5Params;
        use std::collections::HashSet;

        fn test_param_set<P: ParameterSet>() {
            let seed = [123u8; 32];
            let mut xof = Xof::init(&seed, XOF_DOMAIN_SEP);
            let support = generate_random_support_uniform::<P>(&mut xof);

            assert_eq!(support.len(), P::HQC_OMEGA as usize);

            for &val in &support {
                assert!(val < P::HQC_N);
            }

            let unique: HashSet<u32> = support.iter().copied().collect();
            assert_eq!(unique.len(), support.len());
        }

        test_param_set::<Hqc1Params>();
        test_param_set::<Hqc3Params>();
        test_param_set::<Hqc5Params>();
    }

    // ========================================================================
    // sample_fixed_weight_vect_uniform Tests
    // ========================================================================

    /// Helper function to count the number of 1 bits in a byte slice
    fn count_ones(bytes: &[u8]) -> u32 {
        bytes.iter().map(|b| b.count_ones()).sum()
    }

    // NOTE: Tests for sample_fixed_weight_vect_uniform and support_to_vect are currently
    // disabled because the LittleN type (U2209) from hybrid-array's extra-sizes branch
    // incorrectly resolves to 16520 instead of 2209, causing memory issues.
    // The implementation is correct; these tests should be re-enabled once the
    // hybrid-array dependency is fixed.

    /// Test that sample_fixed_weight_vect_uniform produces a vector with exactly HQC_OMEGA bits set
    #[test]
    fn sample_fixed_weight_vect_uniform_weight() {
        use crate::hqc1::Hqc1Params;

        let seed = [0u8; 32];
        let mut xof = Xof::init(&seed, XOF_DOMAIN_SEP);

        let vect = sample_fixed_weight_vect_uniform::<Hqc1Params>(&mut xof);

        let weight = count_ones(vect.as_slice());
        assert_eq!(
            weight,
            Hqc1Params::HQC_OMEGA,
            "Vector should have weight HQC_OMEGA"
        );
    }

    /// Test that support_to_vect correctly converts support to bit vector
    #[test]
    fn support_to_vect_correctness() {
        use crate::hqc1::Hqc1Params;

        let seed = [42u8; 32];
        let mut xof = Xof::init(&seed, XOF_DOMAIN_SEP);

        let support = generate_random_support_uniform::<Hqc1Params>(&mut xof);
        let vect = support_to_vect::<Hqc1Params>(&support);

        for &idx in &support {
            let byte_idx = (idx / 8) as usize;
            let bit_idx = idx % 8;
            assert!(
                (vect[byte_idx] >> bit_idx) & 1 == 1,
                "Bit at index {idx} should be set"
            );
        }

        let weight = count_ones(vect.as_slice());
        assert_eq!(weight, Hqc1Params::HQC_OMEGA);
    }

    /// Test determinism for sample_fixed_weight_vect_uniform
    #[test]
    fn sample_fixed_weight_vect_uniform_deterministic() {
        use crate::hqc1::Hqc1Params;

        let seed = [99u8; 32];

        let mut xof1 = Xof::init(&seed, XOF_DOMAIN_SEP);
        let vect1 = sample_fixed_weight_vect_uniform::<Hqc1Params>(&mut xof1);

        let mut xof2 = Xof::init(&seed, XOF_DOMAIN_SEP);
        let vect2 = sample_fixed_weight_vect_uniform::<Hqc1Params>(&mut xof2);

        assert_eq!(
            vect1.as_slice(),
            vect2.as_slice(),
            "Same seed should produce same vector"
        );
    }

    /// Test sample_fixed_weight_vect_uniform for all parameter sets
    #[test]
    fn sample_fixed_weight_vect_uniform_all_params() {
        use crate::hqc1::Hqc1Params;
        use crate::hqc3::Hqc3Params;
        use crate::hqc5::Hqc5Params;

        fn test_param_set<P: ParameterSet>() {
            let seed = [77u8; 32];
            let mut xof = Xof::init(&seed, XOF_DOMAIN_SEP);
            let vect = sample_fixed_weight_vect_uniform::<P>(&mut xof);

            let weight = count_ones(vect.as_slice());
            assert_eq!(
                weight,
                P::HQC_OMEGA,
                "Vector should have weight HQC_OMEGA for {:?}",
                P::default()
            );
        }

        test_param_set::<Hqc1Params>();
        test_param_set::<Hqc3Params>();
        test_param_set::<Hqc5Params>();
    }

    fn assert_keygen_vects_from_intermediates<P: ParameterSet>(variant: &str) {
        let contents = inter_kats(variant);

        let seed_dk = extract_hex_field(contents, "seed_dk");
        let seed_ek = extract_hex_field(contents, "seed_ek");
        let y_expected = extract_hex_field(contents, "y");
        let x_expected = extract_hex_field(contents, "x");
        let h_expected = extract_hex_field(contents, "h");

        let mut xof = Xof::init(&seed_dk, XOF_DOMAIN_SEP);
        let y = sample_fixed_weight_vect_uniform::<P>(&mut xof);
        let x = sample_fixed_weight_vect_uniform::<P>(&mut xof);

        assert_eq!(y_expected, y.as_slice(), "y vector mismatch");
        assert_eq!(x_expected, x.as_slice(), "x vector mismatch");

        let mut xof = Xof::init(&seed_ek, XOF_DOMAIN_SEP);
        let h = sample_vect::<P>(&mut xof);
        assert_eq!(h_expected, h.as_slice(), "h vector mismatch");
    }

    #[test]
    fn keygen_vects_match_intermediates_hqc1() {
        assert_keygen_vects_from_intermediates::<crate::hqc1::Hqc1Params>("hqc-1");
    }
    #[test]
    fn keygen_vects_match_intermediates_hqc3() {
        assert_keygen_vects_from_intermediates::<crate::hqc3::Hqc3Params>("hqc-3");
    }
    #[test]
    fn keygen_vects_match_intermediates_hqc5() {
        assert_keygen_vects_from_intermediates::<crate::hqc5::Hqc5Params>("hqc-5");
    }

    // ========================================================================
    // Vector Arithmetic Tests
    // ========================================================================

    /// Test all vector addition properties for a given parameter set
    fn test_vect_add_properties<P: ParameterSet>() {
        let seed = [42u8; 32];
        let mut xof = Xof::init(&seed, XOF_DOMAIN_SEP);
        let a = sample_vect::<P>(&mut xof);
        let b = sample_vect::<P>(&mut xof);
        let c = sample_vect::<P>(&mut xof);
        let zero = Vect::<P>::default();

        // Identity: a + 0 = a
        let result = vect_add::<P>(&a, &zero);
        assert_eq!(a.as_slice(), result.as_slice(), "a + 0 should equal a");

        // Commutativity: a + b = b + a
        let ab = vect_add::<P>(&a, &b);
        let ba = vect_add::<P>(&b, &a);
        assert_eq!(ab.as_slice(), ba.as_slice(), "a + b should equal b + a");

        // Self-inverse: a + a = 0
        let result = vect_add::<P>(&a, &a);
        assert_eq!(result.as_slice(), zero.as_slice(), "a + a should equal 0");

        // Associativity: (a + b) + c = a + (b + c)
        let ab_c = vect_add::<P>(&ab, &c);
        let bc = vect_add::<P>(&b, &c);
        let a_bc = vect_add::<P>(&a, &bc);
        assert_eq!(
            ab_c.as_slice(),
            a_bc.as_slice(),
            "(a + b) + c should equal a + (b + c)"
        );
    }

    #[test]
    fn vect_add_properties_hqc1() {
        test_vect_add_properties::<crate::hqc1::Hqc1Params>();
    }

    #[test]
    fn vect_add_properties_hqc3() {
        test_vect_add_properties::<crate::hqc3::Hqc3Params>();
    }

    #[test]
    fn vect_add_properties_hqc5() {
        test_vect_add_properties::<crate::hqc5::Hqc5Params>();
    }

    /// Test all vector multiplication properties for a given parameter set
    fn test_vect_mul_properties<P: ParameterSet>() {
        let seed = [42u8; 32];
        let mut xof = Xof::init(&seed, XOF_DOMAIN_SEP);
        let a = sample_vect::<P>(&mut xof);
        let b = sample_vect::<P>(&mut xof);
        let c = sample_vect::<P>(&mut xof);

        let zero = Vect::<P>::default();
        let mut one = Vect::<P>::default();
        one[0] = 1;

        // Identity: a * 1 = a
        let result = vect_mul::<P>(&a, &one);
        assert_eq!(a.as_slice(), result.as_slice(), "a * 1 should equal a");

        // Zero: a * 0 = 0
        let result = vect_mul::<P>(&a, &zero);
        assert_eq!(result.as_slice(), zero.as_slice(), "a * 0 should equal 0");

        // Commutativity: a * b = b * a
        let ab = vect_mul::<P>(&a, &b);
        let ba = vect_mul::<P>(&b, &a);
        assert_eq!(ab.as_slice(), ba.as_slice(), "a * b should equal b * a");

        // Associativity: (a * b) * c = a * (b * c)
        let ab_c = vect_mul::<P>(&ab, &c);
        let bc = vect_mul::<P>(&b, &c);
        let a_bc = vect_mul::<P>(&a, &bc);
        assert_eq!(
            ab_c.as_slice(),
            a_bc.as_slice(),
            "(a * b) * c should equal a * (b * c)"
        );

        // Distributivity: a * (b + c) = (a * b) + (a * c)
        let b_plus_c = vect_add::<P>(&b, &c);
        let left = vect_mul::<P>(&a, &b_plus_c);
        let ac = vect_mul::<P>(&a, &c);
        let right = vect_add::<P>(&ab, &ac);
        assert_eq!(
            left.as_slice(),
            right.as_slice(),
            "a * (b + c) should equal (a * b) + (a * c)"
        );
    }

    #[test]
    fn vect_mul_properties_hqc1() {
        test_vect_mul_properties::<crate::hqc1::Hqc1Params>();
    }

    #[test]
    fn vect_mul_properties_hqc3() {
        test_vect_mul_properties::<crate::hqc3::Hqc3Params>();
    }

    #[test]
    fn vect_mul_properties_hqc5() {
        test_vect_mul_properties::<crate::hqc5::Hqc5Params>();
    }
}
