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
cargo bench --bench kem
```

#### Informal Results

On my MacBook Pro with Apple M2 Pro, I thought it would be humbling to benchmark this against ML-KEM from aws-lc-rs 1.15.4 at equivalent security levels. I was right. Note that these are entirely different algorithms, so the comparison isn't really fair — but neither are the numbers, so it fits.

**KeyGen**
| Security | syndrmrs | aws-lc-rs | slowdown |
|----------|----------|-----------|----------|
| 128-bit | 723 µs | 8.5 µs | 85x |
| 192-bit | 2.76 ms | 12.8 µs | 216x |
| 256-bit | 7.41 ms | 16.8 µs | 441x |

**Encaps**
| Security | syndrmrs | aws-lc-rs | slowdown |
|----------|----------|-----------|----------|
| 128-bit | 1.38 ms | 10.0 µs | 138x |
| 192-bit | 5.61 ms | 14.0 µs | 401x |
| 256-bit | 14.6 ms | 19.9 µs | 734x |

**Decaps**
| Security | syndrmrs | aws-lc-rs | slowdown |
|----------|----------|-----------|----------|
| 128-bit | 2.46 ms | 8.3 µs | 296x |
| 192-bit | 8.13 ms | 12.9 µs | 630x |
| 256-bit | 22.5 ms | 19.6 µs | 1148x |

~6x improvement from constant-time Karatsuba + clmul64 carry-less multiply 🏎️

---

#### Acknowledgments

The structure of this crate was influenced by [RustCrypto's ml-kem](https://github.com/RustCrypto/KEMs/tree/525c6307021b215f9a3dc4a6f0f63d9dd07a7374/ml-kem). This project is not affiliated with RustCrypto.

#### Disclaimer

This is a personal project and is not affiliated with or endorsed by my employer.
