# syndrmrs

[![Rust](https://github.com/initsecret/syndrmrs/actions/workflows/rust.yml/badge.svg)](https://github.com/initsecret/syndrmrs/actions/workflows/rust.yml)

Educational Rust implementation of Hamming Quasi-Cyclic (HQC) KEM as described in the [2025-08-22 draft](https://pqc-hqc.org/doc/hqc_specifications_2025_08_22.pdf).

## ⚠️ DO NOT USE IN PRODUCTION ⚠️

**This library is for educational purposes only.**

- This is an **experimental** implementation of a **draft** standard
- It has **not** been audited
- Timing side-channels may still exist in some code paths
- The API is unstable and will change without notice
- We make **no security guarantees whatsoever**

If you need a post-quantum KEM in production, use something else.

## Status

HQC-1, HQC-3, and HQC-5 all pass reference KATs.

**Note:** This crate depends on a [fork of `hybrid-array`](https://github.com/initsecret/hybrid-array/tree/add-hqc-sizes) with HQC sizes.

## Testing

```bash
# Run all tests except KATs
cargo test

# Run all KATs (use release mode since KATs are slow in debug)
cargo test --release -- --ignored
```

## Performance

Run benchmarks with:
```bash
# Software carry-less multiply (any platform)
cargo bench --bench kem

# Hardware PCLMULQDQ (x86_64 with SSE2 + PCLMULQDQ)
RUSTFLAGS="-C target-feature=+pclmulqdq" cargo bench --bench kem
```

#### Informal Results

I thought it would be humbling to benchmark this against ML-KEM from aws-lc-rs 1.15.4 at equivalent security levels. I was right. Note that these are entirely different algorithms, so the comparison isn't really fair — but neither are the numbers, so it fits.

**KeyGen**
| Security | syndrmrs | aws-lc-rs | slowdown |
|----------|----------|-----------|----------|
| 128-bit | 74 µs | 22 µs | 3.4x |
| 192-bit | 265 µs | 33 µs | 8.0x |
| 256-bit | 627 µs | 49 µs | 12.8x |

**Encaps**
| Security | syndrmrs | aws-lc-rs | slowdown |
|----------|----------|-----------|----------|
| 128-bit | 146 µs | 25 µs | 5.8x |
| 192-bit | 527 µs | 39 µs | 13.5x |
| 256-bit | 1.24 ms | 53 µs | 23.4x |

**Decaps**
| Security | syndrmrs | aws-lc-rs | slowdown |
|----------|----------|-----------|----------|
| 128-bit | 230 µs | 15 µs | 15.3x |
| 192-bit | 844 µs | 23 µs | 36.7x |
| 256-bit | 2.01 ms | 34 µs | 59.1x |

~55x improvement from constant-time Karatsuba + hardware PCLMULQDQ carry-less multiply 🏎️

---

#### Acknowledgments

The structure of this crate was influenced by [RustCrypto's ml-kem](https://github.com/RustCrypto/KEMs/tree/525c6307021b215f9a3dc4a6f0f63d9dd07a7374/ml-kem). This project is not affiliated with RustCrypto.

#### Disclaimer

This is a personal project and is not affiliated with or endorsed by my employer.
