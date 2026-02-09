//! GF(256) arithmetic over F2^8 with the primitive polynomial specified by HQC.
//!
//! Spec: HQC 2025-08-22, §3.4 (Coding theory preliminaries for GF(256)).
//! Ref: https://gitlab.com/pqc-hqc/hqc/-/blob/d622142a50f3ce6b6e1f5b15a5119d96c67194e0/src/ref/hqc-1/parameters.h (PARAM_GF_POLY = 0x11D),
//!      https://gitlab.com/pqc-hqc/hqc/-/blob/d622142a50f3ce6b6e1f5b15a5119d96c67194e0/src/ref/gf.h (gf_exp/gf_log tables).

use subtle::{ConditionallySelectable, ConstantTimeEq};

#[inline]
#[must_use]
#[allow(dead_code)] // XOR is obvious, but kept for API completeness
pub const fn add(a: u16, b: u16) -> u16 {
    a ^ b
}

/// Multiplies two elements in GF(256).
///
/// Takes u16 arguments (matching the reference implementation) but assumes values are < 256.
/// Returns u16 for consistency with syndrome/polynomial computations.
///
/// Uses the `subtle` crate to select between 0 and the table result without branching.
#[must_use]
pub fn mul(a: u16, b: u16) -> u16 {
    debug_assert!(a < 256 && b < 256, "GF(256) elements must be < 256");

    let log_a = GF_LOG[a as usize] as usize;
    let log_b = GF_LOG[b as usize] as usize;
    let idx = (log_a + log_b) % 255;
    let result = GF_EXP[idx];

    let either_zero = a.ct_eq(&0) | b.ct_eq(&0);
    u16::conditional_select(&result, &0, either_zero)
}

/// Computes the multiplicative inverse of an element in GF(256).
///
/// Uses the `subtle` crate to select between 0 and the table result without branching.
#[must_use]
pub fn inv(a: u16) -> u16 {
    let log_a = GF_LOG[a as usize] as usize;
    let idx = (255 - log_a) % 255;
    let result = GF_EXP[idx];

    let is_zero = a.ct_eq(&0);
    u16::conditional_select(&result, &0, is_zero)
}

/// Squares an element in GF(256).
///
/// In GF(2^m), squaring is a linear operation.
/// We compute a² = a * a using the multiplication function.
#[inline]
#[must_use]
pub fn square(a: u16) -> u16 {
    mul(a, a)
}

pub const GF_EXP: [u16; 258] = [
    1, 2, 4, 8, 16, 32, 64, 128, 29, 58, 116, 232, 205, 135, 19, 38, 76, 152, 45, 90, 180, 117,
    234, 201, 143, 3, 6, 12, 24, 48, 96, 192, 157, 39, 78, 156, 37, 74, 148, 53, 106, 212, 181,
    119, 238, 193, 159, 35, 70, 140, 5, 10, 20, 40, 80, 160, 93, 186, 105, 210, 185, 111, 222, 161,
    95, 190, 97, 194, 153, 47, 94, 188, 101, 202, 137, 15, 30, 60, 120, 240, 253, 231, 211, 187,
    107, 214, 177, 127, 254, 225, 223, 163, 91, 182, 113, 226, 217, 175, 67, 134, 17, 34, 68, 136,
    13, 26, 52, 104, 208, 189, 103, 206, 129, 31, 62, 124, 248, 237, 199, 147, 59, 118, 236, 197,
    151, 51, 102, 204, 133, 23, 46, 92, 184, 109, 218, 169, 79, 158, 33, 66, 132, 21, 42, 84, 168,
    77, 154, 41, 82, 164, 85, 170, 73, 146, 57, 114, 228, 213, 183, 115, 230, 209, 191, 99, 198,
    145, 63, 126, 252, 229, 215, 179, 123, 246, 241, 255, 227, 219, 171, 75, 150, 49, 98, 196, 149,
    55, 110, 220, 165, 87, 174, 65, 130, 25, 50, 100, 200, 141, 7, 14, 28, 56, 112, 224, 221, 167,
    83, 166, 81, 162, 89, 178, 121, 242, 249, 239, 195, 155, 43, 86, 172, 69, 138, 9, 18, 36, 72,
    144, 61, 122, 244, 245, 247, 243, 251, 235, 203, 139, 11, 22, 44, 88, 176, 125, 250, 233, 207,
    131, 27, 54, 108, 216, 173, 71, 142, 1, 2, 4,
];

pub const GF_LOG: [u16; 256] = [
    0, 0, 1, 25, 2, 50, 26, 198, 3, 223, 51, 238, 27, 104, 199, 75, 4, 100, 224, 14, 52, 141, 239,
    129, 28, 193, 105, 248, 200, 8, 76, 113, 5, 138, 101, 47, 225, 36, 15, 33, 53, 147, 142, 218,
    240, 18, 130, 69, 29, 181, 194, 125, 106, 39, 249, 185, 201, 154, 9, 120, 77, 228, 114, 166, 6,
    191, 139, 98, 102, 221, 48, 253, 226, 152, 37, 179, 16, 145, 34, 136, 54, 208, 148, 206, 143,
    150, 219, 189, 241, 210, 19, 92, 131, 56, 70, 64, 30, 66, 182, 163, 195, 72, 126, 110, 107, 58,
    40, 84, 250, 133, 186, 61, 202, 94, 155, 159, 10, 21, 121, 43, 78, 212, 229, 172, 115, 243,
    167, 87, 7, 112, 192, 247, 140, 128, 99, 13, 103, 74, 222, 237, 49, 197, 254, 24, 227, 165,
    153, 119, 38, 184, 180, 124, 17, 68, 146, 217, 35, 32, 137, 46, 55, 63, 209, 91, 149, 188, 207,
    205, 144, 135, 151, 178, 220, 252, 190, 97, 242, 86, 211, 171, 20, 42, 93, 158, 132, 60, 57,
    83, 71, 109, 65, 162, 31, 45, 67, 216, 183, 123, 164, 118, 196, 23, 73, 236, 127, 12, 111, 246,
    108, 161, 59, 82, 41, 157, 85, 170, 251, 96, 134, 177, 187, 204, 62, 90, 203, 89, 95, 176, 156,
    169, 160, 81, 11, 245, 22, 235, 122, 117, 44, 215, 79, 174, 213, 233, 230, 231, 173, 232, 116,
    214, 244, 234, 168, 80, 88, 175,
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_properties() {
        // Addition in GF(256) is XOR
        assert_eq!(add(0, 0), 0);
        assert_eq!(add(1, 1), 0);
        assert_eq!(add(0xFF, 0xFF), 0);
        assert_eq!(add(0x12, 0x34), 0x12 ^ 0x34);
        assert_eq!(add(100, 200), 100 ^ 200);

        for a in 0..=255u16 {
            // Identity: 0 is the additive identity
            assert_eq!(add(a, 0), a);
            assert_eq!(add(0, a), a);

            // Inverse: every element is its own additive inverse
            assert_eq!(add(a, a), 0);

            for b in 0..=255u16 {
                // Commutative
                assert_eq!(add(a, b), add(b, a));
            }
        }

        // Associative (test on subset for performance)
        for a in (0..=255).step_by(17) {
            for b in (0..=255).step_by(17) {
                for c in (0..=255).step_by(17) {
                    assert_eq!(add(add(a, b), c), add(a, add(b, c)));
                }
            }
        }
    }

    #[test]
    fn mul_properties() {
        for a in 0..=255u16 {
            // Zero: multiplication by zero
            assert_eq!(mul(a, 0), 0);
            assert_eq!(mul(0, a), 0);

            // Identity: 1 is the multiplicative identity
            assert_eq!(mul(a, 1), a);
            assert_eq!(mul(1, a), a);

            for b in 0..=255u16 {
                // Commutative
                assert_eq!(mul(a, b), mul(b, a));
            }
        }

        // Associative (test on subset for performance)
        for a in (0..=255).step_by(17) {
            for b in (0..=255).step_by(17) {
                for c in (0..=255).step_by(17) {
                    assert_eq!(mul(mul(a, b), c), mul(a, mul(b, c)));
                }
            }
        }
    }

    #[test]
    fn distributive_law() {
        // a * (b + c) = a * b + a * c
        for a in (0..=255).step_by(17) {
            for b in (0..=255).step_by(17) {
                for c in (0..=255).step_by(17) {
                    let left = mul(a, add(b, c));
                    let right = add(mul(a, b), mul(a, c));
                    assert_eq!(
                        left, right,
                        "Distributive law failed for a={a}, b={b}, c={c}"
                    );
                }
            }
        }
    }

    #[test]
    fn inv_properties() {
        // By convention, inv(0) = 0
        assert_eq!(inv(0), 0);

        // inv(1) = 1
        assert_eq!(inv(1), 1);

        // For all non-zero a: a * inv(a) = 1 and inv(inv(a)) = a
        for a in 1..=255 {
            let inv_a = inv(a);

            // Multiplicative inverse
            assert_eq!(mul(a, inv_a), 1, "inv({a}) = {inv_a} failed");
            assert_eq!(mul(inv_a, a), 1);

            // Involutive: inv(inv(a)) = a
            assert_eq!(inv(inv_a), a);
        }
    }

    #[test]
    fn generator_is_2() {
        // The generator element is 2 (primitive element)
        // 2^255 should equal 1 (order of multiplicative group is 255)
        let mut result = 1u16;
        for _ in 0..255 {
            result = mul(result, 2);
        }
        assert_eq!(result, 1, "Generator 2 should have order 255");
    }

    #[test]
    fn exp_log_properties() {
        // Verify that exp and log tables are consistent
        // For all non-zero a: exp[log[a]] = a
        for a in 1..=255 {
            let log_a = GF_LOG[a as usize];
            let exp_log_a = GF_EXP[log_a as usize];
            assert_eq!(exp_log_a, a, "exp[log[{a}]] = {exp_log_a} != {a}");
        }

        // The exp table extends beyond 255 to allow wrapping
        // exp[i] = exp[i + 255] for i in 0..3 (the wraparound entries)
        for i in 0..3 {
            assert_eq!(GF_EXP[i], GF_EXP[i + 255]);
        }
    }

    #[test]
    fn known_multiplication_values() {
        // Test some known multiplication results
        // These can be computed by hand or verified against reference implementation
        assert_eq!(mul(2, 2), 4);
        assert_eq!(mul(2, 3), 6);
        assert_eq!(mul(3, 3), 5); // 3 * 3 = 9 = 0b1001, reduced by polynomial
        assert_eq!(mul(16, 16), 29); // 16^2 = 256 = 0x100, reduced
    }

    #[test]
    fn known_inverse_values() {
        // Test some known inverses
        assert_eq!(inv(2), 142); // Can be verified: mul(2, 142) should = 1
        assert_eq!(mul(2, 142), 1);

        assert_eq!(inv(3), 244);
        assert_eq!(mul(3, 244), 1);
    }

    #[test]
    fn primitive_polynomial_verification() {
        // The primitive polynomial is 0x11D = x^8 + x^4 + x^3 + x^2 + 1
        // Verify that 2^8 (which would be 256) reduces correctly
        // 2^8 = 256 = 0x100, which should reduce to 0x1D (29) under modulo 0x11D
        let mut result = 1u16;
        for _ in 0..8 {
            result = mul(result, 2);
        }
        assert_eq!(result, 29, "2^8 should reduce to 29 under polynomial 0x11D");
    }
}
