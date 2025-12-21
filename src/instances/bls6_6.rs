//! BLS6_6 curve instance.
//!
//! This module defines the BLS6_6 pairing-friendly curve with:
//! - Base field: F_43 (prime p = 43)
//! - Scalar field: F_13 (group order r = 13)
//! - Embedding degree: 6
//!
//! # Example
//! ```rust,ignore
//! use curvelib::instances::bls6_6::{Bls6_6G1Config, get_g1_curve};
//! use curvelib::traits::ShortWeierstrassConfig;
//!
//! // Access curve parameters via the config trait
//! let a = Bls6_6G1Config::coeff_a();
//! let b = Bls6_6G1Config::coeff_b();
//! ```

use mathlib::{FieldElement, fp};

use crate::algebra::fields::Fp2;
use crate::models::{SexticTwist, TwistPoint, WeierstrassCurve, WeierstrassPoint};
use crate::traits::{Curve, CurveConfig, ShortWeierstrassConfig};

/// Base field configuration for BLS6_6 (modulus p = 43)
#[derive(mathlib::FieldConfig, Clone, Copy, Debug, Default, PartialEq, Eq)]
#[modulus = "0x2B"] // 43 in decimal
pub struct Bls6_6BaseField;

/// Scalar field configuration for BLS6_6 (group order r = 13)
#[derive(mathlib::FieldConfig, Clone, Copy, Debug, Default, PartialEq, Eq)]
#[modulus = "0x0D"] // 13 in decimal
pub struct Bls6_6ScalarField;

// Type aliases for convenience
pub type Fp43 = FieldElement<Bls6_6BaseField>;
pub type Fr13 = FieldElement<Bls6_6ScalarField>;

// Point type aliases
pub type G1Point = WeierstrassPoint<Bls6_6BaseField>;
pub type G2Point = TwistPoint<Bls6_6BaseField>;

/// G1 curve configuration for BLS6_6.
///
/// This implements the arkworks-style `CurveConfig` pattern where curve
/// parameters are defined at the type level.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Bls6_6G1Config;

impl CurveConfig for Bls6_6G1Config {
    type BaseField = Bls6_6BaseField;
    type ScalarField = Bls6_6ScalarField;
    type Projective = crate::models::sw::Projective<Self>;

    /// Cofactor h = 1 (BLS6_6 G1 is prime order)
    const COFACTOR: &'static [u64] = &[1];

    fn generator() -> Self::Projective {
        crate::models::sw::Projective::<Self>::generator()
    }
}

impl ShortWeierstrassConfig for Bls6_6G1Config {
    /// Coefficient a = 0 for y² = x³ + 6
    fn coeff_a() -> FieldElement<Self::BaseField> {
        FieldElement::<Bls6_6BaseField>::zero()
    }

    /// Coefficient b = 6 for y² = x³ + 6
    fn coeff_b() -> FieldElement<Self::BaseField> {
        fp!(6u64, Bls6_6BaseField)
    }

    /// Generator x-coordinate = 13
    fn generator_x() -> FieldElement<Self::BaseField> {
        fp!(13u64, Bls6_6BaseField)
    }

    /// Generator y-coordinate = 15
    fn generator_y() -> FieldElement<Self::BaseField> {
        fp!(15u64, Bls6_6BaseField)
    }

    /// a = 0, so this is true
    fn a_is_zero() -> bool {
        true
    }
}

/// Final Exponentiation Power: (43^6 - 1) / 13 = 486258696
pub const FINAL_EXPONENT: u64 = 486_258_696;

/// Get the G1 curve: y² = x³ + 6 over F_43
///
/// Uses the `Bls6_6G1Config` to construct the curve from type-level parameters.
pub fn get_g1_curve() -> WeierstrassCurve<Bls6_6BaseField> {
    WeierstrassCurve::new(
        Bls6_6G1Config::coeff_a(),
        Bls6_6G1Config::coeff_b(),
        Bls6_6G1Config::generator_x(),
        Bls6_6G1Config::generator_y(),
    )
}

/// Get the G1 generator point (13, 15)
pub fn get_g1_generator() -> G1Point {
    let curve = get_g1_curve();
    curve.generator()
}

/// Get the G2 curve: y² = x³ + 6 over F_43²
pub fn get_g2_curve() -> SexticTwist<Bls6_6BaseField> {
    let zero = Fp2::<Bls6_6BaseField>::zero();
    let six = Fp2::new(
        fp!(6u64, Bls6_6BaseField),
        FieldElement::<Bls6_6BaseField>::zero(),
    );
    let gx = Fp2::new(
        fp!(13u64, Bls6_6BaseField),
        FieldElement::<Bls6_6BaseField>::zero(),
    );
    let gy = Fp2::new(
        fp!(15u64, Bls6_6BaseField),
        FieldElement::<Bls6_6BaseField>::zero(),
    );

    SexticTwist::new(zero, six, gx, gy)
}

/// Get the G2 generator point
pub fn get_g2_generator() -> G2Point {
    let curve = get_g2_curve();
    curve.generator()
}

/// Create a base field element (F_43)
#[inline]
pub fn base_field(value: u64) -> Fp43 {
    fp!(value, Bls6_6BaseField)
}

/// Create a scalar field element (F_13)
#[inline]
pub fn scalar_field(value: u64) -> Fr13 {
    fp!(value, Bls6_6ScalarField)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::traits::ProjectivePoint;

    #[test]
    fn test_field_config_base() {
        let zero = Fp43::zero();
        let one = Fp43::one();
        assert_ne!(zero, one);

        // Test arithmetic in F_43
        let a = fp!(40u64, Bls6_6BaseField);
        let b = fp!(5u64, Bls6_6BaseField);
        let sum = a + b;
        // 40 + 5 = 45 ≡ 2 (mod 43)
        assert_eq!(sum, fp!(2u64, Bls6_6BaseField));
    }

    #[test]
    fn test_field_config_scalar() {
        let a = fp!(10u64, Bls6_6ScalarField);
        let b = fp!(5u64, Bls6_6ScalarField);
        let sum = a + b;
        // 10 + 5 = 15 ≡ 2 (mod 13)
        assert_eq!(sum, fp!(2u64, Bls6_6ScalarField));
    }

    #[test]
    fn test_g1_generator_on_curve() {
        let curve = get_g1_curve();
        let g = curve.generator();
        let (x, y) = g.to_affine();

        // Verify: 13³ + 6 = 2197 + 6 = 2203 ≡ 10 (mod 43)
        // And 15² = 225 ≡ 10 (mod 43)
        assert!(curve.is_on_curve(&x, &y));
    }

    #[test]
    fn test_g1_identity() {
        let curve = get_g1_curve();
        let id = curve.identity();
        assert!(id.is_identity());
    }

    #[test]
    fn test_curve_config_trait() {
        // Test that the CurveConfig trait is properly implemented
        assert!(Bls6_6G1Config::cofactor_is_one());

        // Test ShortWeierstrassConfig
        let a = Bls6_6G1Config::coeff_a();
        let b = Bls6_6G1Config::coeff_b();
        assert!(a.is_zero());
        assert!(!b.is_zero());
        assert!(Bls6_6G1Config::a_is_zero());
    }
}
