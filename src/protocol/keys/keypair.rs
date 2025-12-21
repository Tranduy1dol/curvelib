//! Key pair types for elliptic curve cryptography.

use std::marker::PhantomData;

use mathlib::{FieldElement, U1024};

use crate::traits::{CurveConfig, ProjectivePoint};

/// A private key (scalar in the curve's scalar field).
///
/// The scalar is stored as a `FieldElement<C::ScalarField>` to ensure
/// automatic modular reduction and valid range.
#[derive(Clone, Debug)]
pub struct PrivateKey<C: CurveConfig> {
    scalar: FieldElement<C::ScalarField>,
    _marker: PhantomData<C>,
}

impl<C: CurveConfig> PrivateKey<C> {
    /// Create a new private key from a U1024 scalar.
    ///
    /// The scalar is automatically reduced modulo the curve order.
    pub fn new(scalar: U1024) -> Self {
        Self {
            scalar: FieldElement::<C::ScalarField>::new(scalar),
            _marker: PhantomData,
        }
    }

    /// Create a new private key from a FieldElement.
    pub fn from_field_element(scalar: FieldElement<C::ScalarField>) -> Self {
        Self {
            scalar,
            _marker: PhantomData,
        }
    }

    /// Get the scalar as a FieldElement.
    pub fn scalar(&self) -> &FieldElement<C::ScalarField> {
        &self.scalar
    }

    /// Get the scalar as U1024.
    pub fn to_u1024(&self) -> U1024 {
        self.scalar.to_u1024()
    }

    /// Convert to hex string.
    pub fn to_hex_string(&self) -> String {
        let val = self.scalar.to_u1024();
        let mut hex = String::new();
        let mut started = false;
        for limb in val.0.iter().rev() {
            if *limb != 0 || started {
                if started {
                    hex.push_str(&format!("{:016x}", limb));
                } else {
                    hex.push_str(&format!("{:x}", limb));
                    started = true;
                }
            }
        }
        if hex.is_empty() {
            hex.push('0');
        }
        hex
    }

    /// Create from hex string.
    pub fn from_hex_string(hex: &str) -> Option<Self> {
        let hex = hex.strip_prefix("0x").unwrap_or(hex);
        let hex = hex.strip_prefix("0X").unwrap_or(hex);

        // Validate hex characters
        for c in hex.chars() {
            if !c.is_ascii_hexdigit() {
                return None;
            }
        }

        let scalar = U1024::from_hex(hex);
        Some(Self::new(scalar))
    }

    /// Check if the private key is zero (invalid).
    pub fn is_zero(&self) -> bool {
        self.scalar.is_zero()
    }
}

impl<C: CurveConfig> PartialEq for PrivateKey<C> {
    fn eq(&self, other: &Self) -> bool {
        self.scalar == other.scalar
    }
}

impl<C: CurveConfig> Eq for PrivateKey<C> {}

/// A public key (curve point).
///
/// Wraps a projective point on the curve.
#[derive(Clone, Debug)]
pub struct PublicKey<C: CurveConfig> {
    point: C::Projective,
}

impl<C: CurveConfig> PublicKey<C> {
    /// Create a new public key from a curve point.
    pub fn new(point: C::Projective) -> Self {
        Self { point }
    }

    /// Get the curve point.
    pub fn point(&self) -> &C::Projective {
        &self.point
    }

    /// Check if this is the identity point (invalid public key).
    pub fn is_identity(&self) -> bool {
        self.point.is_identity()
    }
}

impl<C: CurveConfig> PartialEq for PublicKey<C> {
    fn eq(&self, other: &Self) -> bool {
        self.point == other.point
    }
}

impl<C: CurveConfig> Eq for PublicKey<C> {}

/// A key pair containing both private and public keys.
#[derive(Clone, Debug)]
pub struct KeyPair<C: CurveConfig> {
    private_key: PrivateKey<C>,
    public_key: PublicKey<C>,
}

impl<C: CurveConfig> KeyPair<C> {
    /// Create a new key pair.
    pub fn new(private_key: PrivateKey<C>, public_key: PublicKey<C>) -> Self {
        Self {
            private_key,
            public_key,
        }
    }

    /// Get the private key.
    pub fn private_key(&self) -> &PrivateKey<C> {
        &self.private_key
    }

    /// Get the public key.
    pub fn public_key(&self) -> &PublicKey<C> {
        &self.public_key
    }

    /// Consume the key pair and return individual keys.
    pub fn into_parts(self) -> (PrivateKey<C>, PublicKey<C>) {
        (self.private_key, self.public_key)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instances::bls6_6::Bls6_6G1Config;

    #[test]
    fn test_private_key_hex_roundtrip() {
        let scalar = U1024::from_u64(42);
        let key = PrivateKey::<Bls6_6G1Config>::new(scalar);
        let hex = key.to_hex_string();
        let recovered = PrivateKey::<Bls6_6G1Config>::from_hex_string(&hex).unwrap();
        assert_eq!(key, recovered);
    }

    #[test]
    fn test_private_key_modular_reduction() {
        // Test that large scalars are reduced modulo curve order
        let large = U1024::from_u64(100); // larger than order 13
        let key = PrivateKey::<Bls6_6G1Config>::new(large);
        // Should be reduced: 100 mod 13 = 9
        assert_eq!(key.to_u1024(), U1024::from_u64(100 % 13));
    }

    #[test]
    fn test_public_key_from_generator() {
        let g = <Bls6_6G1Config as CurveConfig>::generator();
        let pk = PublicKey::<Bls6_6G1Config>::new(g.clone());
        assert!(!pk.is_identity());
        assert_eq!(pk.point(), &g);
    }
}
