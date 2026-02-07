# syndrmrs

[![Rust](https://github.com/initsecret/syndrmrs/actions/workflows/rust.yml/badge.svg)](https://github.com/initsecret/syndrmrs/actions/workflows/rust.yml)

Educational Rust implementation of Hamming Quasi-Cyclic (HQC) KEM as described in the [2025-08-22 draft](https://pqc-hqc.org/doc/hqc_specifications_2025_08_22.pdf).

## ⚠️ DO NOT USE IN PRODUCTION ⚠️

**This library is for educational purposes only.**

- This is an **experimental** implementation of a **draft** standard
- It has **not** been audited
- It is **not** constant-time (timing side-channels may exist)
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
| 128-bit | 4.0 ms | 8.5 µs | 470x |
| 192-bit | 16.1 ms | 12.8 µs | 1260x |
| 256-bit | 41.4 ms | 16.8 µs | 2460x |

**Encaps**
| Security | syndrmrs | aws-lc-rs | slowdown |
|----------|----------|-----------|----------|
| 128-bit | 7.9 ms | 10.0 µs | 790x |
| 192-bit | 31.8 ms | 14.0 µs | 2270x |
| 256-bit | 82.5 ms | 19.9 µs | 4140x |

**Decaps**
| Security | syndrmrs | aws-lc-rs | slowdown |
|----------|----------|-----------|----------|
| 128-bit | 12.1 ms | 8.3 µs | 1460x |
| 192-bit | 48.3 ms | 12.9 µs | 3740x |
| 256-bit | 123.8 ms | 19.6 µs | 6320x |

syndrmrs is not optimized (yet!) 🤷

---

#### Acknowledgments

The structure of this crate was influenced by [RustCrypto's ml-kem](https://github.com/RustCrypto/KEMs/tree/525c6307021b215f9a3dc4a6f0f63d9dd07a7374/ml-kem). This project is not affiliated with RustCrypto.

#### Disclaimer

This is a personal project and is not affiliated with or endorsed by my employer.
