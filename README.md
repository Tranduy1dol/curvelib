# curvelib

[![CI](https://github.com/Tranduy1dol/curvelib/actions/workflows/ci.yml/badge.svg)](https://github.com/Tranduy1dol/curvelib/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/Tranduy1dol/curvelib/graph/badge.svg?token=ZgGkQ787W4)](https://codecov.io/gh/Tranduy1dol/curvelib)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A Rust library for elliptic curve cryptography, supporting various curve models and cryptographic primitives.

**Note**: This repository is developed for self-education purposes and is not production-ready.

---

## Features

- **Elliptic Curve Models**:
  - **Short Weierstrass**: Generic implementation for curves like secp256k1.
  - **Twisted Edwards**: Support for Edwards curves like Jubjub.
  - **Sextic Twist**: Implementation of twist curves for pairing-friendly constructions.
- **Pairing-Friendly Curves**:
  - Implementation of bilinear pairings (Miller loop).
  - Support for BN and BLS curve families.
- **Cryptographic Primitives**:
  - **Key Management**: `KeyEngine` for secure keypair generation and serialization.
  - **Digital Signatures**: ECDSA-like signing and verification using `SigningEngine`.
  - **Commitments**: Pedersen commitments (in progress).
- **Instances**:
  - **Tiny Jubjub**: A toy example for educational purposes.
  - **BLS6_6**: A small embedding degree curve for testing pairings.

## Installation

Add the following to your `Cargo.toml`:

```toml
[dependencies]
curvelib = { git = "https://github.com/Tranduy1dol/curvelib" }
# Ensure mathlib is also available
mathlib = { git = "https://github.com/Tranduy1dol/mathlib" }
```

## Usage

### Key Generation and Signing

```rust
use curvelib::{
    instances::tiny_jubjub::{TinyJubjubConfig, TinyJubjubScalarField},
    protocol::{
        keys::{KeyEngine, PrivateKey, ToHex},
        signing::SigningEngine,
    },
};
use mathlib::{BigInt, FieldElement, U1024};

fn example() {
    let key_engine = KeyEngine::<TinyJubjubConfig>::new();
    let signing_engine = SigningEngine::<TinyJubjubConfig>::new();

    // 1. Generate Keypair
    let private_scalar = U1024::from_u64(3); // In practice, use a secure random generator
    let keypair = key_engine.keypair_from_scalar(private_scalar);

    // 2. Sign a message
    let message_hash = FieldElement::<TinyJubjubScalarField>::new(U1024::from_u64(42));
    let nonce = FieldElement::<TinyJubjubScalarField>::new(U1024::from_u64(7)); // Use random nonce
    
    let signature = signing_engine
        .sign(&message_hash, keypair.private_key(), &nonce)
        .expect("signing failed");
        
    println!("Signature R: {}", signature.r_hex());
    println!("Signature S: {}", signature.s_hex());

    // 3. Verify
    let valid = signing_engine.verify(&signature, &message_hash, keypair.public_key());
    assert!(valid);
}
```

## Architecture

The library is organized into the following modules:

- **`algebra`**: Low-level algebraic structures including finite fields and arithmetic.
- **`models`**: Generic implementations of elliptic curve models (Weierstrass, Edwards).
- **`instances`**: Concrete instantiations of curves with specific parameters.
- **`protocol`**: High-level cryptographic protocols built on top of the curves.
- **`traits`**: Core traits defining the interfaces for curves, fields, and groups.

## Testing

Run the test suite with:

```bash
cargo test
```

## License

This project is licensed under the MIT License.