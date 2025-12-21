//! Curve trait for elliptic curve definitions.

use std::fmt::Debug;

use mathlib::U1024;

use crate::traits::point::ProjectivePoint;

/// Trait representing a U1024 conversion.
pub trait ToU1024 {
    fn to_u1024(&self) -> U1024;
}

/// Trait for elliptic curves.
///
/// This trait defines the properties and operations available on an elliptic curve,
/// including access to the identity point, generator, and curve parameters.
pub trait Curve: Clone + Debug {
    /// The type of points on this curve.
    type Point: ProjectivePoint;

    /// The base field type for coordinates.
    type BaseField;

    /// Returns the identity point (point at infinity).
    fn identity(&self) -> Self::Point;

    /// Checks if the affine point (x, y) lies on this curve.
    fn is_on_curve(&self, x: &Self::BaseField, y: &Self::BaseField) -> bool;

    /// Returns the curve generator point.
    fn generator(&self) -> Self::Point;

    /// Returns the curve cofactor as little-endian u64 limbs.
    fn cofactor(&self) -> &'static [u64];

    /// Computes the public key from a private scalar.
    fn generate_keypair(&self, private_key: &U1024) -> Self::Point {
        self.generator().mul(private_key)
    }
}
