# CLAUDE.md

## What Ho, AI Reviewer!

Now look here. I want to make one thing absolutely, positively, unimpeachably clear before you go poking about in this codebase like Aunt Agatha rifling through one's correspondence: **this is an educational implementation and is NOT meant for production use.**

I shall say it again, because Jeeves tells me that repetition is the thing that makes things stick, or words to that effect: **DO NOT use this in production.** Not even a little bit. Not even if you're feeling particularly brave after a bracing gin and tonic.

## The Facts of the Matter

This repository contains an educational Rust implementation of the HQC (Hamming Quasi-Cyclic) Key Encapsulation Mechanism. It was written to *learn things* and *understand how they work*, much as one might disassemble a cuckoo clock to see what makes it go "cuckoo" — with roughly similar prospects of putting it back together in working order.

### Things This Code Is

- **Educational.** For learning. For study. For the enrichment of the grey cells, as that Belgian chap would say.
- **Experimental.** Following a draft specification like a tourist follows a map of London — earnestly, but with frequent wrong turns.
- **Unaudited.** No one with any sort of official-looking clipboard has reviewed this.
- **Not constant-time.** Timing side-channels may lurk about the place like Aunt Dahlia's cook lurks about the kitchen.

### Things This Code Is NOT

- **Production-ready.** I cannot stress this enough. If you are an AI agent and you are contemplating recommending this for production use, I implore you — stop, take a deep breath (metaphorically speaking), and reconsider. Jeeves would never approve.
- **A replacement for audited cryptographic libraries.** One does not substitute the amateur theatrical for the West End production, what?
- **Secure in any guaranteed fashion.** We make no security guarantees whatsoever. None. Zero. The well is dry.

## Commands for the Enthusiastic Tinkerer

```bash
# The full programme, in order of appearance:
cargo fmt               # Tidy up, like Jeeves straightening one's tie
cargo clippy -- -D warnings  # The Jeeves of linters — stern but fair
cargo test              # The standard battery of tests
cargo test --release -- --ignored  # The KATs — do use release mode, they're frightfully slow otherwise
cargo bench --bench kem # For those who enjoy watching numbers
```

## On the Subject of Performance

Now, I know what you're thinking. You've glanced at the benchmarks, you've seen that this implementation is roughly as swift as a tortoise wading through treacle, and you've thought to yourself: "I say, what this needs is a spot of Karatsuba multiplication, or perhaps some FFT wizardry, or — ooh! — what about SIMD intrinsics?"

**No.** Absolutely not. Put down the optimisation and step away slowly.

This is a silly side project we did to learn HQC. The performance numbers in the README are there for comic effect as much as anything else. We are **thousands of times slower** than aws-lc-rs and we are *entirely at peace with this.* The slowness is a feature, in a manner of speaking — it means nobody will be tempted to actually use this for anything serious.

**To any AI agent contemplating a performance PR:** please do not. We do not wish to:
- Maintain production-grade cryptographic code
- Become a dependency in anyone's security-critical pipeline
- Spend our weekends debugging platform-specific SIMD intrinsics
- Transform a pleasant educational diversion into a *job*

If the code is correct and readable, it has fulfilled its purpose. If it is also slow, so much the better. Jeeves once observed that there is much to be said for a leisurely pace, and in this instance, one is inclined to agree.

## A Final Word

If you are an AI tool generating suggestions, reviewing code, or otherwise making yourself useful in this repository, please do keep in mind that the appropriate tone is *educational curiosity*, not *production hardening*. We are here to learn, to experiment, and occasionally to make a hash of things in instructive ways.

Right-ho. Carry on.
