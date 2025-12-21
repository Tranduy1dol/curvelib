//! Configuration traits for elliptic curves.
//!
//! This module defines traits that bundle together compile-time constants
//! and types for curve instances. Curves are treated as groups of points,
//! not as objects with cryptographic methods.
//!
//! # Design Note
//!
//! With mathlib 1.0.0, field elements use type-level configuration via
//! `FieldConfig`, enabling const construction and clean generic code.
//!
//! # Example
//!
//! ```rust,ignore
//! use curvelib::traits::{CurveConfig, ShortWeierstrassConfig};
//! use curvelib::instances::bls6_6::Bls6_6G1Config;
//!
//! // Access curve parameters via the config traits
//! let a = Bls6_6G1Config::coeff_a();
//! let b = Bls6_6G1Config::coeff_b();
//! let gx = Bls6_6G1Config::generator_x();
//!
//! // Use helper methods for optimized arithmetic
//! let elem = SomeFieldElement::new(...);
//! let ax = Bls6_6G1Config::mul_by_a(elem);
//! ```

use std::fmt::Debug;

use mathlib::{BigInt, FieldConfig, FieldElement, U1024};

use crate::traits::point::ProjectivePoint;

/// Configuration trait for elliptic curves.
///
/// This trait bundles together the associated types and static parameters
/// that define a specific curve instance. Following the arkworks pattern,
/// curves are treated as algebraic groups rather than runtime objects.
///
/// # Associated Types
///
/// - `BaseField`: The FieldConfig for the base field F_p
/// - `ScalarField`: The FieldConfig for the scalar field F_r
pub trait CurveConfig: Clone + Debug + Sized + Send + Sync + 'static {
    /// The base field configuration.
    type BaseField: FieldConfig;

    /// The scalar field configuration.
    type ScalarField: FieldConfig;

    /// The projective point type for this curve.
    type Projective: ProjectivePoint<Field = FieldElement<Self::BaseField>>;

    /// Returns the generator point of the curve group.
    fn generator() -> Self::Projective;

    /// The cofactor of the curve as little-endian u64 limbs.
    /// The cofactor h = #E(F_p) / r, where r is the prime subgroup order.
    const COFACTOR: &'static [u64];

    /// Check if the cofactor is one (optimization for prime-order curves).
    #[inline]
    fn cofactor_is_one() -> bool {
        Self::COFACTOR.len() == 1 && Self::COFACTOR[0] == 1
    }

    /// Get the cofactor as a U1024 for scalar multiplication.
    fn cofactor_as_u1024() -> U1024 {
        let mut result = U1024::zero();
        for (i, &limb) in Self::COFACTOR.iter().enumerate() {
            if i < 16 {
                result.0[i] = limb;
            }
        }
        result
    }
}

/// Configuration for Short Weierstrass curves: y² = x³ + ax + b
///
/// This trait extends `CurveConfig` with accessors and helpers specific
/// to the Short Weierstrass curve model.
pub trait ShortWeierstrassConfig: CurveConfig {
    /// Get coefficient `a` of the curve equation y² = x³ + ax + b
    fn coeff_a() -> FieldElement<Self::BaseField>;

    /// Get coefficient `b` of the curve equation y² = x³ + ax + b
    fn coeff_b() -> FieldElement<Self::BaseField>;

    /// Get generator point x-coordinate
    fn generator_x() -> FieldElement<Self::BaseField>;

    /// Get generator point y-coordinate
    fn generator_y() -> FieldElement<Self::BaseField>;

    /// Check if the curve coefficient a is zero (optimization).
    fn a_is_zero() -> bool;

    /// Compute `elem * a` efficiently.
    ///
    /// Override for special cases like a = 0 or a = -3.
    #[inline]
    fn mul_by_a(elem: FieldElement<Self::BaseField>) -> FieldElement<Self::BaseField> {
        if Self::a_is_zero() {
            FieldElement::<Self::BaseField>::zero()
        } else {
            elem * Self::coeff_a()
        }
    }

    /// Compute `elem + b` efficiently.
    #[inline]
    fn add_b(elem: FieldElement<Self::BaseField>) -> FieldElement<Self::BaseField> {
        elem + Self::coeff_b()
    }
}

/// Configuration for Twisted Edwards curves: ax² + y² = 1 + dx²y²
///
/// This trait extends `CurveConfig` with the additional parameter `d`
/// and accessors specific to the Twisted Edwards model.
pub trait TwistedEdwardsConfig: CurveConfig {
    /// Get the Edwards parameter `a` in ax² + y² = 1 + dx²y²
    fn coeff_a() -> FieldElement<Self::BaseField>;

    /// Get the Edwards parameter `d` in ax² + y² = 1 + dx²y²
    fn coeff_d() -> FieldElement<Self::BaseField>;

    /// Get generator point x-coordinate
    fn generator_x() -> FieldElement<Self::BaseField>;

    /// Get generator point y-coordinate
    fn generator_y() -> FieldElement<Self::BaseField>;

    /// Check if `a` is -1 (common optimization for Ed curves)
    fn a_is_minus_one() -> bool;

    /// Compute `elem * a` efficiently.
    ///
    /// Override for special cases like a = -1.
    #[inline]
    fn mul_by_a(elem: FieldElement<Self::BaseField>) -> FieldElement<Self::BaseField> {
        if Self::a_is_minus_one() {
            -elem
        } else {
            elem * Self::coeff_a()
        }
    }
}
