//! Key management for elliptic curve cryptography.
//!
//! This module provides types and functionality for managing key pairs,
//! including private keys, public keys, and hex encoding/decoding.

mod hex;
mod keypair;

pub use hex::{FromHex, HexError, ToHex};
pub use keypair::{KeyPair, PrivateKey, PublicKey};

use mathlib::U1024;

use crate::traits::{CurveConfig, ProjectivePoint};

/// Key engine for generating and managing key pairs.
///
/// Generic over `C: CurveConfig` which defines the curve.
pub struct KeyEngine<C: CurveConfig> {
    _marker: std::marker::PhantomData<C>,
}

impl<C: CurveConfig> KeyEngine<C> {
    /// Create a new key engine.
    pub fn new() -> Self {
        Self {
            _marker: std::marker::PhantomData,
        }
    }

    /// Generate a key pair from a private scalar.
    ///
    /// # Parameters
    ///
    /// * `scalar` - The private key scalar (will be reduced mod curve order)
    ///
    /// # Returns
    ///
    /// A `KeyPair` containing the private and public keys.
    pub fn keypair_from_scalar(&self, scalar: U1024) -> KeyPair<C> {
        let private_key = PrivateKey::<C>::new(scalar);
        let generator = C::generator();
        // Use the reduced scalar for multiplication
        let reduced_scalar = private_key.to_u1024();
        let public_point = generator.mul(&reduced_scalar);
        let public_key = PublicKey::new(public_point);

        KeyPair::new(private_key, public_key)
    }

    /// Parse a private key from hex string.
    pub fn private_from_hex(&self, hex: &str) -> Result<PrivateKey<C>, HexError> {
        PrivateKey::<C>::from_hex(hex)
    }

    /// Encode a public key to hex string.
    pub fn public_to_hex(&self, public_key: &PublicKey<C>, compressed: bool) -> String {
        public_key.to_hex(compressed)
    }
}

impl<C: CurveConfig> Default for KeyEngine<C> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instances::bls6_6::Bls6_6G1Config;

    #[test]
    fn test_keypair_generation() {
        let engine = KeyEngine::<Bls6_6G1Config>::new();
        let scalar = U1024::from_u64(5);
        let keypair = engine.keypair_from_scalar(scalar);

        // Scalar should be reduced mod 13 (5 mod 13 = 5)
        assert_eq!(keypair.private_key().to_u1024(), U1024::from_u64(5));
        assert!(!keypair.public_key().point().is_identity());
    }

    #[test]
    fn test_keypair_with_large_scalar() {
        let engine = KeyEngine::<Bls6_6G1Config>::new();
        // Use a scalar larger than the order (13)
        let large_scalar = U1024::from_u64(20);
        let keypair = engine.keypair_from_scalar(large_scalar);

        // Should be reduced: 20 mod 13 = 7
        assert_eq!(keypair.private_key().to_u1024(), U1024::from_u64(20 % 13));
    }
}
