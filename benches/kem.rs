use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};

use aws_lc_rs::kem::{Algorithm, Ciphertext, DecapsulationKey, EncapsulationKey};
use syndrmrs::kem::{HqcKem, KemCiphertext, KemDecapsulationKey, KemEncapsulationKey, KemSeed};
use syndrmrs::{ParameterSet, hqc1, hqc3, hqc5};

// ============================================================================
// HQC Benchmarks
// ============================================================================

/// Generate a deterministic seed for benchmarking
fn bench_seed() -> KemSeed {
    let mut seed = KemSeed::default();
    seed.as_mut_slice().copy_from_slice(&[0x42u8; 32]);
    seed
}

/// Generate deterministic encapsulation inputs for benchmarking
fn bench_encaps_inputs<P: ParameterSet>() -> (syndrmrs::kem::Msg<P>, [u8; 16]) {
    let mut m = syndrmrs::kem::Msg::<P>::default();
    for (i, byte) in m.as_mut_slice().iter_mut().enumerate() {
        *byte = (i & 0xff) as u8;
    }
    let salt = [0x42u8; 16];
    (m, salt)
}

/// Benchmark HQC KEM.KeyGen for a specific parameter set
fn bench_hqc_keygen<P: ParameterSet>(c: &mut Criterion, name: &str) {
    let seed = bench_seed();

    c.bench_with_input(BenchmarkId::new("keygen", name), &seed, |b, seed| {
        b.iter(|| {
            let (_ek, _dk): (KemEncapsulationKey<P>, KemDecapsulationKey<P>) =
                HqcKem::keygen::<P>(seed);
        });
    });
}

/// Benchmark HQC KEM.Encaps for a specific parameter set
fn bench_hqc_encaps<P: ParameterSet>(c: &mut Criterion, name: &str) {
    let seed = bench_seed();
    let (ek, _dk): (KemEncapsulationKey<P>, KemDecapsulationKey<P>) = HqcKem::keygen::<P>(&seed);
    let (m, salt) = bench_encaps_inputs::<P>();

    c.bench_with_input(
        BenchmarkId::new("encaps", name),
        &(&ek, &m, &salt),
        |b, (ek, m, salt)| {
            b.iter(|| {
                let (_ss, _ct): (_, KemCiphertext<P>) =
                    HqcKem::encaps_deterministic::<P>(ek, m, salt);
            });
        },
    );
}

/// Benchmark HQC KEM.Decaps for a specific parameter set
fn bench_hqc_decaps<P: ParameterSet>(c: &mut Criterion, name: &str) {
    let seed = bench_seed();
    let (ek, dk): (KemEncapsulationKey<P>, KemDecapsulationKey<P>) = HqcKem::keygen::<P>(&seed);
    let (m, salt) = bench_encaps_inputs::<P>();
    let (_ss, ct): (_, KemCiphertext<P>) = HqcKem::encaps_deterministic::<P>(&ek, &m, &salt);

    c.bench_with_input(
        BenchmarkId::new("decaps", name),
        &(&dk, &ct),
        |b, (dk, ct)| {
            b.iter(|| {
                let _ss = HqcKem::decaps::<P>(dk, ct);
            });
        },
    );
}

// ============================================================================
// ML-KEM Benchmarks (aws-lc-rs)
// ============================================================================

/// Benchmark ML-KEM KeyGen
fn bench_mlkem_keygen(c: &mut Criterion, alg: &'static Algorithm, name: &str) {
    c.bench_function(&format!("keygen/{name}"), |b| {
        b.iter(|| {
            let _dk = DecapsulationKey::generate(alg).unwrap();
        });
    });
}

/// Benchmark ML-KEM Encaps
fn bench_mlkem_encaps(c: &mut Criterion, alg: &'static Algorithm, name: &str) {
    let dk = DecapsulationKey::generate(alg).unwrap();
    let ek_bytes = dk.encapsulation_key().unwrap().key_bytes().unwrap();
    let ek = EncapsulationKey::new(alg, ek_bytes.as_ref()).unwrap();

    c.bench_with_input(BenchmarkId::new("encaps", name), &ek, |b, ek| {
        b.iter(|| {
            let (_ss, _ct) = ek.encapsulate().unwrap();
        });
    });
}

/// Benchmark ML-KEM Decaps
fn bench_mlkem_decaps(c: &mut Criterion, alg: &'static Algorithm, name: &str) {
    let dk = DecapsulationKey::generate(alg).unwrap();
    let ek_bytes = dk.encapsulation_key().unwrap().key_bytes().unwrap();
    let ek = EncapsulationKey::new(alg, ek_bytes.as_ref()).unwrap();
    let (ct, _ss) = ek.encapsulate().unwrap();
    let ct_bytes: Vec<u8> = ct.as_ref().to_vec();

    c.bench_with_input(
        BenchmarkId::new("decaps", name),
        &(&dk, ct_bytes),
        |b, (dk, ct_bytes): &(&DecapsulationKey, Vec<u8>)| {
            b.iter(|| {
                let ct: Ciphertext = ct_bytes.as_slice().into();
                let _ss = dk.decapsulate(ct).unwrap();
            });
        },
    );
}

// ============================================================================
// Benchmark Groups
// ============================================================================

fn kem_benchmarks(c: &mut Criterion) {
    // HQC-1 (128-bit security)
    bench_hqc_keygen::<hqc1::Hqc1Params>(c, "hqc1");
    bench_hqc_encaps::<hqc1::Hqc1Params>(c, "hqc1");
    bench_hqc_decaps::<hqc1::Hqc1Params>(c, "hqc1");

    // ML-KEM-512 (128-bit security, comparable to HQC-1)
    bench_mlkem_keygen(c, &aws_lc_rs::kem::ML_KEM_512, "ml-kem-512");
    bench_mlkem_encaps(c, &aws_lc_rs::kem::ML_KEM_512, "ml-kem-512");
    bench_mlkem_decaps(c, &aws_lc_rs::kem::ML_KEM_512, "ml-kem-512");

    // HQC-3 (192-bit security)
    bench_hqc_keygen::<hqc3::Hqc3Params>(c, "hqc3");
    bench_hqc_encaps::<hqc3::Hqc3Params>(c, "hqc3");
    bench_hqc_decaps::<hqc3::Hqc3Params>(c, "hqc3");

    // ML-KEM-768 (192-bit security, comparable to HQC-3)
    bench_mlkem_keygen(c, &aws_lc_rs::kem::ML_KEM_768, "ml-kem-768");
    bench_mlkem_encaps(c, &aws_lc_rs::kem::ML_KEM_768, "ml-kem-768");
    bench_mlkem_decaps(c, &aws_lc_rs::kem::ML_KEM_768, "ml-kem-768");

    // HQC-5 (256-bit security)
    bench_hqc_keygen::<hqc5::Hqc5Params>(c, "hqc5");
    bench_hqc_encaps::<hqc5::Hqc5Params>(c, "hqc5");
    bench_hqc_decaps::<hqc5::Hqc5Params>(c, "hqc5");

    // ML-KEM-1024 (256-bit security, comparable to HQC-5)
    bench_mlkem_keygen(c, &aws_lc_rs::kem::ML_KEM_1024, "ml-kem-1024");
    bench_mlkem_encaps(c, &aws_lc_rs::kem::ML_KEM_1024, "ml-kem-1024");
    bench_mlkem_decaps(c, &aws_lc_rs::kem::ML_KEM_1024, "ml-kem-1024");
}

criterion_group!(benches, kem_benchmarks);
criterion_main!(benches);
