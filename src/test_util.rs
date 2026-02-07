use rand_chacha::rand_core::{RngCore, SeedableRng};

/// Deterministic tiny RNG for repeatable tests
pub struct TestRng {
    rng: rand_chacha::ChaCha8Rng,
}
impl TestRng {
    pub fn new() -> Self {
        let rng = rand_chacha::ChaCha8Rng::seed_from_u64(0xdead_beef);
        Self { rng }
    }
    pub fn next_u32(&mut self) -> u32 {
        self.rng.next_u32()
    }
    pub fn gen_usize(&mut self, upper: usize) -> usize {
        (self.next_u32() as usize) % upper
    }
}

impl Default for TestRng {
    fn default() -> Self {
        Self::new()
    }
}

use std::vec::Vec;

/// Get contents of intermediate KAT file for given variant
pub fn inter_kats(variant: &str) -> &'static str {
    match variant {
        "hqc-1" => include_str!("../tests/kats/ref/hqc-1/intermediates_values"),
        "hqc-3" => include_str!("../tests/kats/ref/hqc-3/intermediates_values"),
        "hqc-5" => include_str!("../tests/kats/ref/hqc-5/intermediates_values"),
        _ => panic!("unknown variant: {variant}"),
    }
}

/// Extract a hex field from an intermediate KAT file
///
/// Searches for a line with format "key: value" and returns the hex-decoded value.
/// Panics if the key is not found.
pub fn extract_hex_field(contents: &str, key: &str) -> Vec<u8> {
    for line in contents.lines() {
        let line = line.trim();
        let Some((label, value)) = line.split_once(':') else {
            continue;
        };
        if label.trim() != key {
            continue;
        }
        let value = value.trim();
        if value.is_empty() {
            continue;
        }
        if !value.chars().all(|c| c.is_ascii_hexdigit()) {
            continue;
        }
        assert!(
            value.len() % 2 == 0,
            "hex string length must be even for {key}"
        );
        return hex::decode(value).expect("valid hex string");
    }
    panic!("missing hex field: {key}");
}
