//! Known Answer Tests (KATs) from PQCkemKAT_*.rsp files
//!
//! These tests validate the full KEM operations against reference test vectors.

use sha3::{
    Shake256,
    digest::{ExtendableOutput, Update, XofReader},
};
use std::collections::HashMap;
use syndrmrs::kem::HqcKem;
use syndrmrs::kem::{KemDecapsulationKey, KemEncapsulationKey};

/// PRNG for KAT testing - matches reference implementation
struct KatPrng {
    reader: sha3::Shake256Reader,
}

impl KatPrng {
    fn new(seed: &[u8]) -> Self {
        const PRNG_DOMAIN: u8 = 0;
        let mut hasher = Shake256::default();
        hasher.update(seed);
        hasher.update(&[PRNG_DOMAIN]);
        Self {
            reader: hasher.finalize_xof(),
        }
    }

    fn get_bytes(&mut self, out: &mut [u8]) {
        self.reader.read(out);
    }
}

/// Parse a PQCkemKAT_*.rsp file into test vectors
fn parse_kat_file(contents: &str) -> Vec<HashMap<String, Vec<u8>>> {
    let mut vectors = Vec::new();
    let mut current = HashMap::new();

    for line in contents.lines() {
        let line = line.trim();

        // Skip comments and empty lines
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        // Parse key = value
        if let Some((key, value)) = line.split_once('=') {
            let key = key.trim();
            let value = value.trim();

            // Decode hex value
            if !value.is_empty()
                && let Ok(bytes) = hex::decode(value)
            {
                current.insert(key.to_string(), bytes);

                // After 'ss' field, we have a complete test vector
                if key == "ss" {
                    vectors.push(current.clone());
                    current.clear();
                }
            }
        }
    }

    vectors
}

/// Run a single KAT test vector
macro_rules! run_kat_vector {
    ($i:expr, $vector:expr, $Params:ty, $m_len:expr) => {{
        use hybrid_array::Array;

        let prng_seed = &$vector["seed"];
        let pk_expected = &$vector["pk"];
        let sk_expected = &$vector["sk"];
        let ct_expected = &$vector["ct"];
        let ss_expected = &$vector["ss"];

        // Initialize PRNG with 48-byte seed
        let mut prng = KatPrng::new(prng_seed);

        // Generate keygen seed (32 bytes)
        let mut seed_kem_bytes = [0u8; 32];
        prng.get_bytes(&mut seed_kem_bytes);
        let seed_kem = Array::from(seed_kem_bytes);

        // Test keygen
        let (ek, dk): (KemEncapsulationKey<$Params>, KemDecapsulationKey<$Params>) =
            HqcKem::keygen(&seed_kem);

        let pk_bytes = ek.to_bytes();
        let sk_bytes = dk.to_bytes();

        assert_eq!(
            pk_expected.as_slice(),
            pk_bytes.as_slice(),
            "Vector {}: pk mismatch",
            $i
        );
        assert_eq!(
            sk_expected.as_slice(),
            sk_bytes.as_slice(),
            "Vector {}: sk mismatch",
            $i
        );

        // Generate message m
        let mut m_bytes = [0u8; $m_len];
        prng.get_bytes(&mut m_bytes);
        let m = Array::from(m_bytes);

        // Generate salt (16 bytes)
        let mut salt = [0u8; 16];
        prng.get_bytes(&mut salt);

        // Test encaps with deterministic values
        let (ss, ct) = HqcKem::encaps_deterministic::<$Params>(&ek, &m, &salt);
        let ct_bytes = ct.to_bytes();

        assert_eq!(
            ct_expected.as_slice(),
            ct_bytes.as_slice(),
            "Vector {}: ct mismatch",
            $i
        );
        assert_eq!(
            ss_expected.as_slice(),
            ss.as_slice(),
            "Vector {}: ss mismatch",
            $i
        );

        // Test decaps
        let ss_decaps = HqcKem::decaps::<$Params>(&dk, &ct);
        assert_eq!(
            ss.as_slice(),
            ss_decaps.as_slice(),
            "Vector {}: decaps ss mismatch",
            $i
        );
    }};
}

#[test]
#[ignore = "very slow on debug builds"]
fn test_hqc1_kats() {
    use syndrmrs::hqc1::Hqc1Params;

    let contents = include_str!("kats/ref/hqc-1/PQCkemKAT_2321.rsp");
    let vectors = parse_kat_file(contents);

    assert_eq!(vectors.len(), 100, "Expected 100 test vectors for HQC-1");

    for (i, vector) in vectors.iter().enumerate() {
        run_kat_vector!(i, vector, Hqc1Params, 16);
    }
}

#[test]
#[ignore = "very slow on debug builds"]
fn test_hqc3_kats() {
    use syndrmrs::hqc3::Hqc3Params;

    let contents = include_str!("kats/ref/hqc-3/PQCkemKAT_4602.rsp");
    let vectors = parse_kat_file(contents);

    assert_eq!(vectors.len(), 100, "Expected 100 test vectors for HQC-3");

    for (i, vector) in vectors.iter().enumerate() {
        run_kat_vector!(i, vector, Hqc3Params, 24);
    }
}

#[test]
#[ignore = "very slow on debug builds"]
fn test_hqc5_kats() {
    use syndrmrs::hqc5::Hqc5Params;

    let contents = include_str!("kats/ref/hqc-5/PQCkemKAT_7333.rsp");
    let vectors = parse_kat_file(contents);

    assert_eq!(vectors.len(), 100, "Expected 100 test vectors for HQC-5");

    for (i, vector) in vectors.iter().enumerate() {
        run_kat_vector!(i, vector, Hqc5Params, 32);
    }
}
