//! Cryptographic protocol implementations.
//!
//! This module provides various cryptographic protocols built on elliptic curves:
//!
//! - **commitment**: Polynomial commitment schemes (KZG, etc.)
//! - **keys**: Key pair management (private/public keys with ScalarField)
//! - **signing**: ECDSA signature generation and verification
//! - **pairing**: Bilinear pairing operations

pub mod commitment;
pub mod keys;
pub mod pairing;
pub mod signing;

// Re-export commonly used items
pub use keys::{HexError, KeyEngine, KeyPair, PrivateKey, PublicKey};
pub use signing::{Signature, SignatureError, SigningEngine};

// Re-export commitment traits and KZG types for convenience
pub use commitment::kzg::{G1Point, G2Point, Kzg, KzgParams};
pub use commitment::{BatchCommitment, PolynomialCommitment};
