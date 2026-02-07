//! Shortened Reed–Solomon over GF(256).
//!
//! Spec: HQC 2025-08-22, §3.4.2 (Reed–Solomon decoding).
//! Ref: https://gitlab.com/pqc-hqc/hqc/-/blob/d622142a50f3ce6b6e1f5b15a5119d96c67194e0/src/ref/reed_solomon.c (encode, compute_syndromes, Berlekamp, FFT roots).
//!
//! Implements:
//! - systematic encoding using generator polynomial from spec
//! - decoding:
//!   1) syndromes
//!   2) Berlekamp (or Berlekamp–Massey) -> error locator polynomial
//!   3) root finding (via fft::find_roots_additive) with bruteforce oracle in tests
//!   4) compute magnitudes, correct
//!

use crate::param::RS2Delta;
use crate::{
    ParameterSet, fft, gf256,
    param::{
        Msg, RSPoly, RsCodeword, RsErrorEvaluatorPoly, RsErrorLocations, RsErrorLocatorPoly,
        RsErrorValues, RsSyndromes,
    },
};
use hybrid_array::typenum::Unsigned;
use subtle::{ConditionallySelectable, ConstantTimeEq, ConstantTimeGreater};

pub struct RS;
impl RS {
    /// Encodes a message of K bytes to a Reed-Solomon codeword of N1 bytes.
    ///
    /// Following the reference implementation, we perform systematic encoding using
    /// a linear (N1 - K)-stage shift register with feedback connections based on
    /// the generator polynomial.
    ///
    /// The systematic format means the message appears in the last K positions,
    /// and the parity symbols are in the first (N1 - K) positions.
    ///
    /// # High-Level Encoding
    ///
    /// ```text
    /// Input message (K bytes):
    /// ┌────────────────────────────┐
    /// │ m[0] │ m[1] │ ... │ m[K-1] │
    /// └────────────────────────────┘
    ///          │
    ///          │ (process through shift register)
    ///          ↓
    /// Output codeword (N1 bytes):
    /// ┌──────────────────────────┬────────────────────────────┐
    /// │ Parity (2*delta bytes)   │ Message (K bytes)          │
    /// │ p[0] │ p[1] │ ... │ p[P] │ m[0] │ m[1] │ ... │ m[K-1] │
    /// └──────────────────────────┴────────────────────────────┘
    ///   ↑ computed via LFSR        ↑ original message (systematic)
    /// ```
    ///
    /// # Shift Register Visualization
    ///
    /// ```text
    /// Input: msg[K-1-i] ──┬──> XOR ──> gate_value
    ///                     │            │
    ///                     │            ├──> × g[0] ──> codeword[0]
    ///                     │            │
    ///                     │            ├──> × g[1] ──> XOR ──> codeword[1]
    ///                     │            │                │
    ///                     │            │           codeword[0] (shifted)
    ///                     │            │
    ///                     │            ├──> × g[2] ──> XOR ──> codeword[2]
    ///                     │            │                │
    ///                     │            │           codeword[1] (shifted)
    ///                     │            │
    ///                     │           ...              ...
    ///                     │            │
    ///                     │            └──> × g[parity_len-1] ──> XOR ──> codeword[parity_len-1]
    ///                     │                                         │
    ///                     └─────────────────────────────────────────┘
    ///                                                         codeword[parity_len-2] (shifted)
    /// ```
    ///
    /// After processing all K message bytes, the parity symbols are in codeword[0..parity_len],
    /// and we copy the message to codeword[parity_len..N1] to form the systematic codeword.
    pub fn encode<P: ParameterSet>(msg: &Msg<P>) -> RsCodeword<P> {
        let g: &RSPoly<P> = P::rs_poly();
        let k = msg.len();
        let n1 = k + g.len() - 1;
        let parity_len = n1 - k;

        let mut codeword = RsCodeword::<P>::default();

        // Perform systematic encoding using shift register
        // Process message bytes from right to left (MSB first)
        for i in 0..k {
            let gate_value = msg[k - 1 - i] ^ codeword[parity_len - 1];

            // Shift right and XOR with g[j] * gate_value
            // We iterate in reverse order (parity_len-1 down to 1), reading from codeword[j-1]
            // and writing to codeword[j]. This is safe without a temporary array because:
            // - At iteration j, we read codeword[j-1] (which hasn't been modified yet)
            // - Then we write to codeword[j] (which we won't read again in this loop)
            for j in (1..parity_len).rev() {
                codeword[j] =
                    codeword[j - 1] ^ (gf256::mul(u16::from(gate_value), u16::from(g[j])) as u8);
            }
            codeword[0] = gf256::mul(u16::from(gate_value), u16::from(g[0])) as u8;
        }

        // Copy message to the end of the codeword (systematic encoding)
        codeword[parity_len..n1].copy_from_slice(msg.as_slice());

        codeword
    }

    /// Decodes a received Reed-Solomon codeword to recover the original message.
    ///
    /// This function implements the standard RS decoding algorithm with error correction
    /// up to delta errors, where delta is the error correction capacity from Table 3.
    ///
    /// # High-Level Decoding Pipeline
    ///
    /// ```text
    /// Input: Received codeword (N1 bytes, possibly corrupted)
    /// ┌──────────────────────────┬────────────────────────────┐
    /// │ Parity (2*delta bytes)   │ Message (K bytes)          │
    /// │ p[0] │ p[1] │ ... │ p[P] │ m[0] │ m[1] │ ... │ m[K-1] │
    /// └──────────────────────────┴────────────────────────────┘
    ///          │
    ///          │ Step 1: Compute Syndromes (2*delta syndromes)
    ///          ↓
    /// ┌────────────────────────────────┐
    /// │ S[0], S[1], ..., S[2*delta-1]  │  ← All zeros if no errors
    /// └────────────────────────────────┘
    ///          │
    ///          │ Step 2: Berlekamp-Massey Algorithm
    ///          ↓
    /// ┌────────────────────────────────┐
    /// │ Error Locator Polynomial σ(x)  │  ← Degree ≤ delta
    /// │ σ(x) = σ[0] + σ[1]x + ... + σ[d]x^d
    /// └────────────────────────────────┘
    ///          │
    ///          │ Step 3: Find Roots via Additive FFT
    ///          ↓
    /// ┌────────────────────────────────┐
    /// │ Error Locations                │  ← Positions of errors
    /// │ e[i] = 1 if position i has error
    /// └────────────────────────────────┘
    ///          │
    ///          │ Step 4: Compute Error Evaluator Polynomial z(x)
    ///          ↓
    /// ┌────────────────────────────────┐
    /// │ z(x) = z[0] + z[1]x + ...      │
    /// └────────────────────────────────┘
    ///          │
    ///          │ Step 5: Compute Error Values (Forney Algorithm)
    ///          ↓
    /// ┌────────────────────────────────┐
    /// │ Error Values e_val[i]          │  ← Magnitude at each error position
    /// └────────────────────────────────┘
    ///          │
    ///          │ Step 6: Correct Errors (XOR error values)
    ///          ↓
    /// ┌─────────────────────────────┐
    /// │ Corrected Message (K bytes) │
    /// │ m[0] │ m[1] │ ... │ m[K-1]  │
    /// └─────────────────────────────┘
    /// ```
    ///
    /// # Syndrome Computation
    ///
    /// For each syndrome index i ∈ [0, 2*delta):
    /// ```text
    /// S[i] = Σ(j=0 to N1-1) cdw[j] · α^(i·j)
    ///      = cdw[0] ⊕ cdw[1]·α^i ⊕ cdw[2]·α^(2i) ⊕ ... ⊕ cdw[N1-1]·α^(i·(N1-1))
    ///
    /// where α is the primitive element of GF(256) and operations are in GF(256).
    /// We use the precomputed powers α^(i·j) from alpha_ij_pow[i][j-1].
    /// ```
    ///
    /// # Error Correction Guarantee
    ///
    /// The RS code can correct up to `delta` symbol errors, where:
    /// - HQC-1: delta = 15 (can correct up to 15 byte errors)
    /// - HQC-3: delta = 16 (can correct up to 16 byte errors)
    /// - HQC-5: delta = 29 (can correct up to 29 byte errors)
    ///
    /// If more than delta errors occur, decoding may fail or produce incorrect results.
    pub fn decode<P: ParameterSet>(encoded: &RsCodeword<P>) -> Msg<P> {
        // Step 1: Compute syndromes
        let syndromes = Self::compute_syndromes::<P>(encoded);

        // Step 2: Compute error locator polynomial using Berlekamp-Massey
        let (sigma, degree) = Self::compute_elp::<P>(&syndromes);

        // Step 3: Find error locations via additive FFT
        // Note: fft::fft() does not take parameter P because it's a parameter-independent
        // operation on GF(256). It always evaluates at all 256 field elements regardless
        // of which HQC variant is being used.
        let error_locations = fft::fft(&sigma, (degree + 1) as usize);

        // Step 4: Compute error evaluator polynomial z(x)
        let z = Self::compute_z_poly::<P>(&sigma, degree, &syndromes);

        // Step 5: Compute error values using Forney algorithm
        let error_values = Self::compute_error_values::<P>(&z, &error_locations);

        // Step 6: Correct errors
        let corrected = Self::correct_errors::<P>(encoded, &error_values);

        // Extract message from systematic codeword (last K bytes)
        let parity_len = crate::param::RSPolyLen::<P>::USIZE - 1; // 2*delta
        let mut msg = Msg::<P>::default();
        msg.copy_from_slice(&corrected[parity_len..]);

        msg
    }

    /// Computes the 2*delta syndromes from a received codeword.
    ///
    /// Syndromes are used to detect errors. All syndromes are zero if no errors occurred.
    ///
    /// # Algorithm
    /// For each syndrome index i ∈ [0, 2*delta):
    /// ```text
    /// S[i] = cdw[0] ⊕ Σ(j=1 to N1-1) cdw[j] · α^(i·j)
    /// ```
    /// where operations are in GF(256) and α^(i·j) is precomputed in `alpha_ij_pow[i][j-1]`.
    ///
    /// # Reference
    /// See `compute_syndromes` in `src/ref/reed_solomon.c:120-127`
    fn compute_syndromes<P: ParameterSet>(codeword: &RsCodeword<P>) -> RsSyndromes<P> {
        let alpha_ij_pow = P::alpha_ij_pow();
        let mut syndromes = RsSyndromes::<P>::default();

        // For each syndrome index i ∈ [0, 2*delta)
        for i in 0..syndromes.len() {
            // S[i] = Σ(j=1 to N1-1) cdw[j] · α^(i·j)
            for j in 1..codeword.len() {
                syndromes[i] ^= gf256::mul(u16::from(codeword[j]), alpha_ij_pow[i][j - 1]);
            }
            // S[i] ^= cdw[0]
            syndromes[i] ^= u16::from(codeword[0]);
        }

        syndromes
    }

    /// Computes the error locator polynomial σ(x) using Berlekamp-Massey algorithm.
    ///
    /// The error locator polynomial is a key component in Reed-Solomon decoding. Its roots
    /// correspond to the positions of errors in the received codeword. For a codeword with
    /// `e` errors at positions `i₁, i₂, ..., iₑ`, the error locator polynomial is:
    ///
    /// ```text
    /// σ(x) = (1 - α^(i₁)·x)(1 - α^(i₂)·x)...(1 - α^(iₑ)·x)
    /// ```
    ///
    /// where α is the primitive element of GF(256).
    ///
    /// # Algorithm Pseudocode
    ///
    /// ```text
    /// Input: Syndromes S[0], S[1], ..., S[2δ-1]
    ///
    /// 1. Initialize:
    ///    σ(x) = 1                          // Error locator polynomial
    ///    σₚ(x) = x                         // Auxiliary polynomial
    ///    deg(σ) = 0, deg(σₚ) = 0           // Polynomial degrees
    ///    ρ = -1                            // Last iteration where degree increased
    ///    dₚ = 1                            // Previous discrepancy
    ///    d = S[0]                          // Current discrepancy
    ///
    /// 2. For μ = 0 to 2δ-1:
    ///
    ///    a) Save current state:
    ///       σ_copy = σ(x)
    ///       deg_copy = deg(σ)
    ///
    ///    b) Update σ(x) to reduce discrepancy:
    ///       σ(x) ← σ(x) - (d/dₚ)·x·σₚ(x)
    ///
    ///    c) Check if degree should increase:
    ///       deg_x = μ - ρ                  // Degree of x^(μ-ρ)
    ///       deg_x_σₚ = deg_x + deg(σₚ)     // Degree of x^(μ-ρ)·σₚ(x)
    ///
    ///       If d ≠ 0 AND deg_x_σₚ > deg(σ):
    ///         deg(σ) ← deg_x_σₚ            // Increase degree
    ///         ρ ← μ                        // Record iteration
    ///         dₚ ← d                       // Save discrepancy
    ///         σₚ(x) ← σ_copy               // Update auxiliary polynomial
    ///         deg(σₚ) ← deg_copy
    ///       Else:
    ///         σₚ(x) ← x·σₚ(x)              // Just shift auxiliary polynomial
    ///
    ///    d) Compute next discrepancy:
    ///       d ← S[μ+1] + Σᵢ₌₁^(μ+1) σᵢ·S[μ+1-i]
    ///
    /// 3. Return: σ(x) and deg(σ)
    ///    - deg(σ) equals the number of errors detected
    /// ```
    ///
    /// # Berlekamp-Massey Intuition
    ///
    /// The algorithm iteratively builds σ(x) by:
    /// 1. Computing a "discrepancy" d that measures how well the current σ(x) fits the syndromes
    /// 2. If d ≠ 0, adjusting σ(x) using a previously saved polynomial σₚ(x)
    /// 3. The adjustment is scaled by d/dₚ to minimize the discrepancy
    ///
    /// The key insight: when the degree increases, we save the current σ(x) as σₚ(x) for
    /// future corrections. This creates a feedback mechanism that efficiently finds the
    /// minimal polynomial that satisfies all syndrome equations.
    ///
    /// # Constant-Time Implementation
    ///
    /// This implementation uses the `subtle` crate for constant-time operations:
    ///
    /// ```text
    /// // Check both conditions without branching
    /// let should_update = !d.ct_eq(&0u16) & deg_x_σₚ.ct_gt(&deg_σ);
    /// // Evaluates to true if: d ≠ 0 AND deg(x·σₚ) > deg(σ)
    ///
    /// // Conditional selection without branching
    /// value = u16::conditional_select(&old_value, &new_value, should_update);
    /// // Returns new_value if should_update is true, old_value otherwise
    /// ```
    ///
    /// This ensures constant execution time regardless of the error pattern, which is
    /// important for cryptographic applications to prevent timing attacks.
    ///
    /// # Example: Single Error
    ///
    /// For a codeword with 1 error at position j:
    /// ```text
    /// Syndromes: S[i] = e·α^(i·j) for i ∈ [0, 2δ)
    /// where e is the error magnitude
    ///
    /// After iteration μ = 0:
    /// σ(x) = 1 - (S[1]/S[0])·x
    ///      = 1 - α^j·x
    ///
    /// The root of σ(x) is α^(-j), which identifies position j
    /// ```
    ///
    /// # Returns
    ///
    /// `(sigma, degree)` where:
    /// - `sigma`: Error locator polynomial coefficients (256 elements for FFT compatibility)
    ///   - Only coefficients `sigma[0..=degree]` are meaningful
    ///   - `sigma[0]` is always 1 (constant term)
    ///   - Remaining elements are zero-padded for FFT input
    /// - `degree`: Actual degree of σ(x), which equals the number of errors detected
    ///   - `degree = 0` means no errors (σ(x) = 1)
    ///   - `degree ≤ delta` (error correction capacity)
    ///   - If `degree > delta`, decoding may fail
    ///
    /// # Reference
    ///
    /// See `compute_elp` in `src/ref/reed_solomon.c:143-205`
    ///
    /// # Implementation Notes
    ///
    /// - Uses 256-element array for FFT compatibility (requires power-of-2 sizing)
    /// - Intermediate arrays (sigma_copy, x_sigma_p) sized for max delta = 29
    /// - All GF(256) operations use precomputed log/exp tables
    /// - Wrapping arithmetic handles the ρ = -1 initialization correctly
    fn compute_elp<P: ParameterSet>(syndromes: &RsSyndromes<P>) -> (RsErrorLocatorPoly, u16) {
        let delta = P::RSDelta::USIZE as u16;
        let two_delta = RS2Delta::<P>::USIZE as u16;

        let mut sigma = RsErrorLocatorPoly::default();
        let mut sigma_copy = [0u16; 64]; // Max delta is 29, so 64 is enough
        let mut x_sigma_p = [0u16; 64];

        // Initialize: σ(x) = 1, σₚ(x) = x
        sigma[0] = 1;
        x_sigma_p[1] = 1;

        let mut deg_sigma = 0u16;
        let mut deg_sigma_p = 0u16;
        let mut pp = u16::MAX; // 2ρ, initialized to -1 (wraps to max)
        let mut d_p = 1u16;
        let mut d = syndromes[0];

        // Berlekamp-Massey iteration: build σ(x) from syndromes
        for mu in 0..two_delta {
            // Save current state before modification
            let delta_usize = delta as usize;
            sigma_copy[..=delta_usize].copy_from_slice(&sigma[..=delta_usize]);
            let deg_sigma_copy = deg_sigma;

            // Update σ(x) to reduce discrepancy: σ(x) ← σ(x) - (d/dₚ)·x·σₚ(x)
            let dd = gf256::mul(d, gf256::inv(d_p));
            let update_limit = (mu + 1).min(delta) as usize;

            for i in 1..=update_limit {
                sigma[i] ^= gf256::mul(dd, x_sigma_p[i]);
            }

            // Check if degree should increase: d ≠ 0 AND deg(x·σₚ) > deg(σ)
            let deg_x = mu.wrapping_sub(pp);
            let deg_x_sigma_p = deg_x.wrapping_add(deg_sigma_p);
            let should_update = !d.ct_eq(&0u16) & deg_x_sigma_p.ct_gt(&deg_sigma);

            // Conditionally update degree
            deg_sigma = u16::conditional_select(&deg_sigma, &deg_x_sigma_p, should_update);

            // Last iteration: no need to prepare for next
            if mu == two_delta - 1 {
                break;
            }

            // Conditionally update ρ and dₚ
            pp = u16::conditional_select(&pp, &mu, should_update);
            d_p = u16::conditional_select(&d_p, &d, should_update);

            // Conditionally update σₚ(x): either σₚ ← σ_copy or σₚ ← x·σₚ
            for i in (1..=delta_usize).rev() {
                x_sigma_p[i] = u16::conditional_select(
                    &x_sigma_p[i - 1],  // if !should_update: shift x·σₚ
                    &sigma_copy[i - 1], // if should_update: use σ_copy
                    should_update,
                );
            }
            x_sigma_p[0] = 0; // x·σₚ(x) has no constant term

            deg_sigma_p = u16::conditional_select(&deg_sigma_p, &deg_sigma_copy, should_update);

            // Compute next discrepancy: d ← S[μ+1] + Σᵢ σᵢ·S[μ+1-i]
            d = syndromes[(mu + 1) as usize];
            for i in 1..=update_limit {
                d ^= gf256::mul(sigma[i], syndromes[(mu + 1) as usize - i]);
            }
        }

        (sigma, deg_sigma)
    }

    /// Computes the error evaluator polynomial z(x).
    ///
    /// # Algorithm
    /// ```text
    /// z[0] = 1
    /// For i in 1..=delta:
    ///     z[i] = sigma[i] (if i <= degree)
    ///     z[i] ^= syndromes[i-1]
    ///     For j in 1..i:
    ///         z[i] ^= gf256::mul(sigma[j], syndromes[i-j-1])
    /// ```
    fn compute_z_poly<P: ParameterSet>(
        sigma: &RsErrorLocatorPoly,
        degree: u16,
        syndromes: &RsSyndromes<P>,
    ) -> RsErrorEvaluatorPoly<P> {
        let delta = RS2Delta::<P>::USIZE / 2;
        let mut z = RsErrorEvaluatorPoly::<P>::default();

        z[0] = 1;

        for i in 1..=delta {
            // Use constant-time mask: mask = 0xffff if i <= degree, else 0
            // C code: mask = -((uint16_t)(i - degree - 1) >> 15)
            // This is 0xffff when (i - degree - 1) is negative (i.e., i <= degree)
            let mask = ((i as u16).wrapping_sub(degree).wrapping_sub(1) >> 15).wrapping_neg();
            z[i] = mask & sigma[i];
        }

        z[1] ^= syndromes[0];

        for i in 2..=delta {
            let mask = ((i as u16).wrapping_sub(degree).wrapping_sub(1) >> 15).wrapping_neg();
            z[i] ^= mask & syndromes[i - 1];

            for j in 1..i {
                z[i] ^= mask & gf256::mul(sigma[j], syndromes[i - j - 1]);
            }
        }

        z
    }

    /// Computes error values using Forney's algorithm.
    ///
    /// Computes the magnitude of errors at each error position.
    ///
    /// # Algorithm
    /// - Computes β_j values (error locator numbers) from error positions
    /// - Evaluates z(x) and σ'(x) at error locations
    /// - Error magnitude = z(β_j^(-1)) / σ'(β_j^(-1))
    fn compute_error_values<P: ParameterSet>(
        z: &RsErrorEvaluatorPoly<P>,
        error_locations: &RsErrorLocations,
    ) -> RsErrorValues<P> {
        let delta = RS2Delta::<P>::USIZE / 2;
        let n1 = P::LittleN1::USIZE;

        let mut error_values = RsErrorValues::<P>::default();
        let mut beta_j = [0u16; 128]; // Max PARAM_DELTA for any parameter set
        let mut e_j = [0u16; 128];

        // Step 1: Compute beta_j values (error locator numbers)
        // beta_j[i] = gf_exp[position of i-th error]
        let mut delta_counter = 0usize;
        for i in 0..n1 {
            let mut found = 0u16;
            // mask1 = 0xffff if error[i] != 0, else 0
            let mask1 = (i32::from(error_locations[i]).wrapping_neg() >> 31) as u16;

            for (j, beta_val) in beta_j.iter_mut().enumerate().take(delta) {
                // mask2 = 0xffff if j == delta_counter, else 0
                let mask2 = !((((j ^ delta_counter) as i32).wrapping_neg() >> 31) as u16);
                *beta_val = beta_val.wrapping_add(mask1 & mask2 & gf256::GF_EXP[i]);
                found = found.wrapping_add(mask1 & mask2 & 1);
            }
            delta_counter = delta_counter.wrapping_add(found as usize);
        }
        let delta_real_value = delta_counter;

        // Step 2: Compute e_j values using Forney's algorithm
        for i in 0..delta {
            let mut tmp1 = 1u16;
            let mut tmp2 = 1u16;
            let inverse = gf256::inv(beta_j[i]);
            let mut inverse_power_j = 1u16;

            // Evaluate z(inverse) = z[0] + z[1]*inverse + ... + z[delta]*inverse^delta
            for j in 1..=delta {
                inverse_power_j = gf256::mul(inverse_power_j, inverse);
                tmp1 ^= gf256::mul(inverse_power_j, z[j]);
            }

            // Compute derivative of error locator polynomial at inverse
            // σ'(inverse) = product of (1 + inverse * beta_j[k]) for k != i
            for k in 1..delta {
                tmp2 = gf256::mul(tmp2, 1 ^ gf256::mul(inverse, beta_j[(i + k) % delta]));
            }

            // mask1 = 0xffff if i < delta_real_value, else 0
            let mask1 = (((i as i16).wrapping_sub(delta_real_value as i16)) >> 15) as u16;
            e_j[i] = mask1 & gf256::mul(tmp1, gf256::inv(tmp2));
        }

        // Step 3: Place error values at correct positions
        delta_counter = 0;
        for i in 0..n1 {
            let mut found = 0u16;
            let mask1 = (i32::from(error_locations[i]).wrapping_neg() >> 31) as u16;

            for (j, &e_val) in e_j.iter().enumerate().take(delta) {
                let mask2 = !((((j ^ delta_counter) as i32).wrapping_neg() >> 31) as u16);
                error_values[i] = error_values[i].wrapping_add(mask1 & mask2 & e_val);
                found = found.wrapping_add(mask1 & mask2 & 1);
            }
            delta_counter = delta_counter.wrapping_add(found as usize);
        }

        error_values
    }

    /// Corrects errors in the codeword by XORing with error values.
    ///
    /// # Arguments
    /// - `codeword`: The received codeword (possibly corrupted)
    /// - `error_values`: The error magnitudes at each position
    ///
    /// # Returns
    /// A corrected codeword with errors fixed
    ///
    /// # Algorithm
    /// ```text
    /// For i in 0..N1:
    ///     corrected[i] = codeword[i] ^ error_values[i]
    /// ```
    fn correct_errors<P: ParameterSet>(
        codeword: &RsCodeword<P>,
        error_values: &RsErrorValues<P>,
    ) -> RsCodeword<P> {
        let n1 = P::LittleN1::USIZE;
        let mut corrected = codeword.clone();
        for i in 0..n1 {
            corrected[i] ^= error_values[i] as u8;
        }
        corrected
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ParameterSet;
    use crate::param::{Msg, RsCodeword};
    use crate::test_util::TestRng;
    use crate::test_util::{extract_hex_field, inter_kats};

    fn assert_rs_encode_from_intermediates<P: ParameterSet>(variant: &str) {
        let contents = inter_kats(variant);

        let msg_bytes = extract_hex_field(contents, "m");
        let rs_codeword_expected = extract_hex_field(contents, "Reed-Solomon code word");

        // Convert the message bytes to Msg<P>
        let mut msg = Msg::<P>::default();
        msg.as_mut_slice().copy_from_slice(&msg_bytes);

        // Encode using our RS implementation
        let encoded = RS::encode::<P>(&msg);

        assert_eq!(encoded.as_slice(), rs_codeword_expected.as_slice());
    }

    fn assert_rs_decode_from_intermediates<P: ParameterSet>(variant: &str) {
        let contents = inter_kats(variant);

        let msg_expected = extract_hex_field(contents, "m");
        let rs_codeword_bytes = extract_hex_field(contents, "Reed-Solomon code word");

        // Convert the RS codeword bytes to RsCodeword<P>
        let mut rs_codeword = RsCodeword::<P>::default();
        rs_codeword
            .as_mut_slice()
            .copy_from_slice(&rs_codeword_bytes);

        // Decode using our RS implementation
        let decoded = RS::decode::<P>(&rs_codeword);

        assert_eq!(decoded.as_slice(), msg_expected.as_slice());
    }

    #[test]
    fn rs_encode_matches_intermediates_hqc1() {
        assert_rs_encode_from_intermediates::<crate::hqc1::Hqc1Params>("hqc-1");
    }

    #[test]
    fn rs_encode_matches_intermediates_hqc3() {
        assert_rs_encode_from_intermediates::<crate::hqc3::Hqc3Params>("hqc-3");
    }

    #[test]
    fn rs_encode_matches_intermediates_hqc5() {
        assert_rs_encode_from_intermediates::<crate::hqc5::Hqc5Params>("hqc-5");
    }

    #[test]
    fn rs_decode_matches_intermediates_hqc1() {
        assert_rs_decode_from_intermediates::<crate::hqc1::Hqc1Params>("hqc-1");
    }

    #[test]
    fn rs_decode_matches_intermediates_hqc3() {
        assert_rs_decode_from_intermediates::<crate::hqc3::Hqc3Params>("hqc-3");
    }

    #[test]
    fn rs_decode_matches_intermediates_hqc5() {
        assert_rs_decode_from_intermediates::<crate::hqc5::Hqc5Params>("hqc-5");
    }

    #[test]
    fn rs_encode_decode_roundtrip_hqc1() {
        rs_roundtrip_test::<crate::hqc1::Hqc1Params>();
    }

    #[test]
    fn rs_encode_decode_roundtrip_hqc3() {
        rs_roundtrip_test::<crate::hqc3::Hqc3Params>();
    }

    #[test]
    fn rs_encode_decode_roundtrip_hqc5() {
        rs_roundtrip_test::<crate::hqc5::Hqc5Params>();
    }

    fn rs_roundtrip_test<P: ParameterSet>() {
        let mut rng = TestRng::new();

        // Test with 100 random messages
        for _ in 0..100 {
            let mut msg = Msg::<P>::default();
            // Fill with random bytes
            for byte in msg.as_mut_slice() {
                *byte = (rng.next_u32() & 0xff) as u8;
            }

            // Encode and decode
            let encoded = RS::encode::<P>(&msg);
            let decoded = RS::decode::<P>(&encoded);

            assert_eq!(
                decoded.as_slice(),
                msg.as_slice(),
                "Failed to roundtrip message"
            );
        }
    }

    #[test]
    fn rs_decode_with_errors_hqc1() {
        rs_error_correction_test::<crate::hqc1::Hqc1Params>();
    }

    #[test]
    fn rs_decode_with_errors_hqc3() {
        rs_error_correction_test::<crate::hqc3::Hqc3Params>();
    }

    #[test]
    fn rs_decode_with_errors_hqc5() {
        rs_error_correction_test::<crate::hqc5::Hqc5Params>();
    }

    fn rs_error_correction_test<P: ParameterSet>() {
        let mut rng = TestRng::new();

        // Test that RS can correct up to delta errors
        // For hqc-1: delta = 15, for hqc-3: delta = 16, for hqc-5: delta = 29
        // We'll introduce delta/2 errors to be safe
        let delta = std::any::type_name::<P>();
        let num_errors = if delta.contains("Hqc1") {
            7 // delta=15, test with 7 errors
        } else if delta.contains("Hqc3") {
            8 // delta=16, test with 8 errors
        } else {
            14 // delta=29, test with 14 errors
        };

        for _ in 0..20 {
            let mut msg = Msg::<P>::default();
            for byte in msg.as_mut_slice() {
                *byte = (rng.next_u32() & 0xff) as u8;
            }

            let mut encoded = RS::encode::<P>(&msg);

            // Introduce random errors in random positions
            let mut error_positions = std::vec::Vec::new();
            while error_positions.len() < num_errors {
                let pos = rng.gen_usize(encoded.as_slice().len());
                if !error_positions.contains(&pos) {
                    error_positions.push(pos);
                }
            }

            for &pos in &error_positions {
                // Flip a random byte value (make sure it actually changes)
                let error = ((rng.next_u32() & 0xfe) as u8) + 1; // Range: 1-255
                encoded.as_mut_slice()[pos] ^= error;
            }

            let decoded = RS::decode::<P>(&encoded);

            assert_eq!(
                decoded.as_slice(),
                msg.as_slice(),
                "Failed to decode with {num_errors} errors"
            );
        }
    }

    // Tests for compute_syndromes

    #[test]
    fn syndromes_zero_for_valid_codeword_hqc1() {
        syndromes_zero_for_valid_codeword::<crate::hqc1::Hqc1Params>();
    }

    #[test]
    fn syndromes_zero_for_valid_codeword_hqc3() {
        syndromes_zero_for_valid_codeword::<crate::hqc3::Hqc3Params>();
    }

    #[test]
    fn syndromes_zero_for_valid_codeword_hqc5() {
        syndromes_zero_for_valid_codeword::<crate::hqc5::Hqc5Params>();
    }

    fn syndromes_zero_for_valid_codeword<P: ParameterSet>() {
        let mut rng = TestRng::new();

        // Test that syndromes are all zero for valid codewords
        for _ in 0..50 {
            let mut msg = Msg::<P>::default();
            for byte in msg.as_mut_slice() {
                *byte = (rng.next_u32() & 0xff) as u8;
            }

            let encoded = RS::encode::<P>(&msg);
            let syndromes = RS::compute_syndromes::<P>(&encoded);

            // All syndromes should be zero for a valid codeword
            for (i, &syndrome) in syndromes.as_slice().iter().enumerate() {
                assert_eq!(
                    syndrome, 0,
                    "Syndrome {i} should be zero for valid codeword"
                );
            }
        }
    }

    #[test]
    fn syndromes_nonzero_for_corrupted_codeword_hqc1() {
        syndromes_nonzero_for_corrupted_codeword::<crate::hqc1::Hqc1Params>();
    }

    #[test]
    fn syndromes_nonzero_for_corrupted_codeword_hqc3() {
        syndromes_nonzero_for_corrupted_codeword::<crate::hqc3::Hqc3Params>();
    }

    #[test]
    fn syndromes_nonzero_for_corrupted_codeword_hqc5() {
        syndromes_nonzero_for_corrupted_codeword::<crate::hqc5::Hqc5Params>();
    }

    fn syndromes_nonzero_for_corrupted_codeword<P: ParameterSet>() {
        let mut rng = TestRng::new();

        // Test that syndromes are non-zero when errors are introduced
        for _ in 0..50 {
            let mut msg = Msg::<P>::default();
            for byte in msg.as_mut_slice() {
                *byte = (rng.next_u32() & 0xff) as u8;
            }

            let mut encoded = RS::encode::<P>(&msg);

            // Introduce a single error
            let error_pos = rng.gen_usize(encoded.as_slice().len());
            let error_value = ((rng.next_u32() & 0xff) as u8) + 1; // Non-zero error
            encoded.as_mut_slice()[error_pos] ^= error_value;

            let syndromes = RS::compute_syndromes::<P>(&encoded);

            // At least one syndrome should be non-zero
            let has_nonzero = syndromes.as_slice().iter().any(|&s| s != 0);
            assert!(
                has_nonzero,
                "At least one syndrome should be non-zero for corrupted codeword"
            );
        }
    }

    #[test]
    fn syndromes_match_kat_hqc1() {
        syndromes_match_kat::<crate::hqc1::Hqc1Params>("hqc-1");
    }

    #[test]
    fn syndromes_match_kat_hqc3() {
        syndromes_match_kat::<crate::hqc3::Hqc3Params>("hqc-3");
    }

    #[test]
    fn syndromes_match_kat_hqc5() {
        syndromes_match_kat::<crate::hqc5::Hqc5Params>("hqc-5");
    }

    fn syndromes_match_kat<P: ParameterSet>(variant: &str) {
        let contents = inter_kats(variant);

        // Extract RS codeword from KAT
        let rs_codeword_bytes = extract_hex_field(contents, "Reed-Solomon code word");
        let mut rs_codeword = RsCodeword::<P>::default();
        rs_codeword
            .as_mut_slice()
            .copy_from_slice(&rs_codeword_bytes);

        // Compute syndromes
        let syndromes = RS::compute_syndromes::<P>(&rs_codeword);

        // Extract expected syndromes from KAT
        // Format: "The syndromes: 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
        let syndrome_line = contents
            .lines()
            .find(|line| line.starts_with("The syndromes:"))
            .expect("Could not find syndrome line in KAT");

        let syndrome_values: std::vec::Vec<u16> = syndrome_line
            .strip_prefix("The syndromes:")
            .unwrap()
            .split_whitespace()
            .map(|s| s.parse::<u16>().expect("Failed to parse syndrome value"))
            .collect();

        // Verify syndromes match
        assert_eq!(
            syndromes.as_slice().len(),
            syndrome_values.len(),
            "Syndrome count mismatch"
        );

        for (i, (&computed, &expected)) in syndromes
            .as_slice()
            .iter()
            .zip(syndrome_values.iter())
            .enumerate()
        {
            assert_eq!(
                computed, expected,
                "Syndrome {i} mismatch: computed={computed}, expected={expected}"
            );
        }
    }

    // Tests for compute_elp (Error Locator Polynomial)

    #[test]
    fn elp_is_one_for_valid_codeword_hqc1() {
        elp_is_one_for_valid_codeword::<crate::hqc1::Hqc1Params>();
    }

    #[test]
    fn elp_is_one_for_valid_codeword_hqc3() {
        elp_is_one_for_valid_codeword::<crate::hqc3::Hqc3Params>();
    }

    #[test]
    fn elp_is_one_for_valid_codeword_hqc5() {
        elp_is_one_for_valid_codeword::<crate::hqc5::Hqc5Params>();
    }

    fn elp_is_one_for_valid_codeword<P: ParameterSet>() {
        let mut rng = TestRng::new();

        // Test that ELP is σ(x) = 1 (degree 0) for valid codewords
        for _ in 0..50 {
            let mut msg = Msg::<P>::default();
            for byte in msg.as_mut_slice() {
                *byte = (rng.next_u32() & 0xff) as u8;
            }

            let encoded = RS::encode::<P>(&msg);
            let syndromes = RS::compute_syndromes::<P>(&encoded);
            let (sigma, degree) = RS::compute_elp::<P>(&syndromes);

            // For valid codeword: σ(x) = 1, so sigma[0] = 1 and degree = 0
            assert_eq!(degree, 0, "Degree should be 0 for valid codeword");
            assert_eq!(sigma[0], 1, "sigma[0] should be 1 for valid codeword");

            // All other coefficients should be 0
            for i in 1..10 {
                assert_eq!(sigma[i], 0, "sigma[{i}] should be 0 for valid codeword");
            }
        }
    }

    #[test]
    fn elp_degree_matches_error_count_hqc1() {
        elp_degree_matches_error_count::<crate::hqc1::Hqc1Params>();
    }

    #[test]
    fn elp_degree_matches_error_count_hqc3() {
        elp_degree_matches_error_count::<crate::hqc3::Hqc3Params>();
    }

    #[test]
    fn elp_degree_matches_error_count_hqc5() {
        elp_degree_matches_error_count::<crate::hqc5::Hqc5Params>();
    }

    fn elp_degree_matches_error_count<P: ParameterSet>() {
        let mut rng = TestRng::new();

        // Test that ELP degree matches the number of errors introduced
        // We'll test with 1, 2, and 3 errors
        for num_errors in 1..=3 {
            for _ in 0..20 {
                let mut msg = Msg::<P>::default();
                for byte in msg.as_mut_slice() {
                    *byte = (rng.next_u32() & 0xff) as u8;
                }

                let mut encoded = RS::encode::<P>(&msg);

                // Introduce exactly num_errors errors
                let mut error_positions = std::vec::Vec::new();
                while error_positions.len() < num_errors {
                    let pos = rng.gen_usize(encoded.as_slice().len());
                    if !error_positions.contains(&pos) {
                        error_positions.push(pos);
                    }
                }

                for &pos in &error_positions {
                    // Generate a non-zero error value (1-255)
                    let error = ((rng.next_u32() % 255) as u8) + 1;
                    encoded.as_mut_slice()[pos] ^= error;
                }

                let syndromes = RS::compute_syndromes::<P>(&encoded);
                let (_sigma, degree) = RS::compute_elp::<P>(&syndromes);

                assert_eq!(
                    degree as usize, num_errors,
                    "ELP degree should match number of errors"
                );
            }
        }
    }

    #[test]
    fn elp_matches_kat_hqc1() {
        elp_matches_kat::<crate::hqc1::Hqc1Params>("hqc-1");
    }

    #[test]
    fn elp_matches_kat_hqc3() {
        elp_matches_kat::<crate::hqc3::Hqc3Params>("hqc-3");
    }

    #[test]
    fn elp_matches_kat_hqc5() {
        elp_matches_kat::<crate::hqc5::Hqc5Params>("hqc-5");
    }

    fn elp_matches_kat<P: ParameterSet>(variant: &str) {
        let contents = inter_kats(variant);

        // Extract RS codeword from KAT
        let rs_codeword_bytes = extract_hex_field(contents, "Reed-Solomon code word");
        let mut rs_codeword = RsCodeword::<P>::default();
        rs_codeword
            .as_mut_slice()
            .copy_from_slice(&rs_codeword_bytes);

        // Compute syndromes and ELP
        let syndromes = RS::compute_syndromes::<P>(&rs_codeword);
        let (sigma, degree) = RS::compute_elp::<P>(&syndromes);

        // Extract expected ELP from KAT
        // Format: "The error locator polynomial: sigma(x) = 1" (for no errors)
        // or "The error locator polynomial: sigma(x) = coef0 + coef1*x + ..."
        let elp_line = contents
            .lines()
            .find(|line| line.starts_with("The error locator polynomial:"))
            .expect("Could not find ELP line in KAT");

        // For valid codeword, it should be "sigma(x) = 1"
        if elp_line.contains("sigma(x) = 1") && !elp_line.contains('+') {
            assert_eq!(degree, 0, "Degree should be 0 for sigma(x) = 1");
            assert_eq!(sigma[0], 1, "sigma[0] should be 1");
        }
        // If there are errors, the format would be different
        // For now, we only test the no-error case from the KAT
    }
}
