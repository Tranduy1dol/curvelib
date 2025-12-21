//! Cryptographic protocol implementations.
//!
//! This module provides various cryptographic protocols built on elliptic curves:
//!
//! - **commitment**: Polynomial commitment schemes (KZG, etc.)
//! - **keys**: Key pair management (private/public keys with ScalarField)
//! - **signing**: ECDSA signature generation and verification
//! - **pairing**: Bilinear pairing operations
//! - **kzg** (legacy): Direct KZG access for backward compatibility

pub mod commitment;
pub mod keys;
pub mod kzg; // Legacy module - use commitment::kzg instead
pub mod pairing;
pub mod signing;

// Re-export commonly used items
pub use keys::{HexError, KeyEngine, KeyPair, PrivateKey, PublicKey};
pub use signing::{Signature, SignatureError, SigningEngine};

// Re-export commitment traits
pub use commitment::{BatchCommitment, PolynomialCommitment};
