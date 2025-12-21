//! Tiny Jubjub curve instance.
//!
//! This module defines a small Twisted Edwards curve for testing:
//! - Base field: F_13 (prime p = 13)
//! - Scalar field: F_5 (group order r = 5)
//! - Curve: 3x² + y² = 1 + 8x²y²
//! - Generator: (6, 9)
//!
//! # Example
//! ```rust,ignore
//! use curvelib::instances::tiny_jubjub::TinyJubjubConfig;
//! use curvelib::traits::TwistedEdwardsConfig;
//!
//! let a = TinyJubjubConfig::coeff_a();
//! let d = TinyJubjubConfig::coeff_d();
//! ```

use mathlib::{FieldElement, fp};

use crate::models::EdwardsCurve;
use crate::traits::{CurveConfig, TwistedEdwardsConfig};

/// Base field configuration for Tiny Jubjub (modulus p = 13)
#[derive(mathlib::FieldConfig, Clone, Copy, Debug, Default, PartialEq, Eq)]
#[modulus = "0x0D"] // 13 in decimal
pub struct TinyJubjubBaseField;

/// Scalar field configuration for Tiny Jubjub (group order r = 5)
#[derive(mathlib::FieldConfig, Clone, Copy, Debug, Default, PartialEq, Eq)]
#[modulus = "0x05"] // 5 in decimal
pub struct TinyJubjubScalarField;

// Type aliases for convenience
pub type Fp13 = FieldElement<TinyJubjubBaseField>;
pub type Fr5 = FieldElement<TinyJubjubScalarField>;

/// Curve configuration for Tiny Jubjub.
///
/// This implements the arkworks-style `CurveConfig` pattern where curve
/// parameters are defined at the type level.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct TinyJubjubConfig;

impl CurveConfig for TinyJubjubConfig {
    type BaseField = TinyJubjubBaseField;
    type ScalarField = TinyJubjubScalarField;
    type Projective = crate::models::te::Projective<Self>;

    /// Cofactor h = 1 (prime order subgroup)
    const COFACTOR: &'static [u64] = &[1];

    fn generator() -> Self::Projective {
        crate::models::te::Projective::<Self>::generator()
    }
}

impl TwistedEdwardsConfig for TinyJubjubConfig {
    /// Coefficient a = 3 for 3x² + y² = 1 + 8x²y²
    fn coeff_a() -> FieldElement<Self::BaseField> {
        fp!(3u64, TinyJubjubBaseField)
    }

    /// Coefficient d = 8 for 3x² + y² = 1 + 8x²y²
    fn coeff_d() -> FieldElement<Self::BaseField> {
        fp!(8u64, TinyJubjubBaseField)
    }

    /// Generator x-coordinate = 6
    fn generator_x() -> FieldElement<Self::BaseField> {
        fp!(6u64, TinyJubjubBaseField)
    }

    /// Generator y-coordinate = 9
    fn generator_y() -> FieldElement<Self::BaseField> {
        fp!(9u64, TinyJubjubBaseField)
    }

    /// a = 3, not -1
    fn a_is_minus_one() -> bool {
        false
    }
}

/// Generator point x-coordinate (6)
pub const GENERATOR_X: u64 = 6;

/// Generator point y-coordinate (9)
pub const GENERATOR_Y: u64 = 9;

/// Get the Tiny Jubjub curve: 3x² + y² = 1 + 8x²y² over F_13
///
/// Uses the `TinyJubjubConfig` to construct the curve from type-level parameters.
pub fn get_curve() -> EdwardsCurve<TinyJubjubBaseField> {
    EdwardsCurve::new(
        TinyJubjubConfig::coeff_a(),
        TinyJubjubConfig::coeff_d(),
        TinyJubjubConfig::generator_x(),
        TinyJubjubConfig::generator_y(),
    )
}

/// Create a base field element (F_13)
#[inline]
pub fn base_field(value: u64) -> Fp13 {
    fp!(value, TinyJubjubBaseField)
}

/// Create a scalar field element (F_5)
#[inline]
pub fn scalar_field(value: u64) -> Fr5 {
    fp!(value, TinyJubjubScalarField)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::traits::Curve;

    #[test]
    fn test_field_config() {
        let zero = Fp13::zero();
        let one = Fp13::one();
        assert_ne!(zero, one);

        // Test arithmetic in F_13
        let a = fp!(10u64, TinyJubjubBaseField);
        let b = fp!(5u64, TinyJubjubBaseField);
        let sum = a + b;
        // 10 + 5 = 15 ≡ 2 (mod 13)
        assert_eq!(sum, fp!(2u64, TinyJubjubBaseField));
    }

    #[test]
    fn test_generator_on_curve() {
        let curve = get_curve();
        let g = curve.generator();
        let (x, y) = crate::traits::ProjectivePoint::to_affine(&g);
        assert!(curve.is_on_curve(&x, &y));
    }

    #[test]
    fn test_curve_config_trait() {
        // Test that the CurveConfig trait is properly implemented
        assert!(TinyJubjubConfig::cofactor_is_one());

        // Test TwistedEdwardsConfig
        let a = TinyJubjubConfig::coeff_a();
        let d = TinyJubjubConfig::coeff_d();
        assert!(!a.is_zero());
        assert!(!d.is_zero());
        assert!(!TinyJubjubConfig::a_is_minus_one());
    }
}
