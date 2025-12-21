//! Projective point trait for elliptic curve operations.

use std::fmt::Debug;

use mathlib::U1024;

use crate::traits::field::Field;

/// Trait for points on an elliptic curve in projective coordinates.
///
/// This trait defines the operations available on curve points, including
/// addition, doubling, scalar multiplication, and negation.
pub trait ProjectivePoint: Sized + Clone + Debug + PartialEq + Eq {
    /// The field type for coordinates.
    type Field: Field;

    /// Returns true if this is the identity (point at infinity).
    fn is_identity(&self) -> bool;

    /// Adds this point to another point.
    fn add(&self, rhs: &Self) -> Self;

    /// Doubles this point (computes 2P).
    fn double(&self) -> Self;

    /// Converts from projective to affine coordinates (x, y).
    fn to_affine(&self) -> (Self::Field, Self::Field);

    /// Scalar multiplication: computes \[`scalar`\] * P.
    fn mul(&self, scalar: &U1024) -> Self;

    /// Negates this point.
    fn neg(&self) -> Self;
}
