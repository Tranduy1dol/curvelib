//! Short Weierstrass curve implementation: y² = x³ + ax + b.
//!
//! This module provides projective and affine point types for
//! working with elliptic curves in Short Weierstrass form.
//!
//! # Design
//!
//! Following the arkworks pattern, curve parameters are defined at the type level
//! via `ShortWeierstrassConfig`. Points are generic over the config type, not
//! over runtime curve instances.

use std::marker::PhantomData;
use std::ops::Neg;

use mathlib::{FieldConfig, FieldElement, U1024};

use crate::traits::{Curve, Field, ShortWeierstrassConfig};

/// A point on a Short Weierstrass curve in affine coordinates (x, y).
///
/// Generic over `P: ShortWeierstrassConfig` which defines the curve parameters.
#[derive(Clone, Copy, Debug)]
pub struct Affine<P: ShortWeierstrassConfig> {
    /// X coordinate
    pub x: FieldElement<P::BaseField>,
    /// Y coordinate
    pub y: FieldElement<P::BaseField>,
    /// Whether this is the point at infinity
    pub infinity: bool,
    _marker: PhantomData<P>,
}

impl<P: ShortWeierstrassConfig> PartialEq for Affine<P> {
    fn eq(&self, other: &Self) -> bool {
        if self.infinity && other.infinity {
            true
        } else if self.infinity || other.infinity {
            false
        } else {
            self.x == other.x && self.y == other.y
        }
    }
}

impl<P: ShortWeierstrassConfig> Eq for Affine<P> {}

impl<P: ShortWeierstrassConfig> Affine<P> {
    /// Create the identity (point at infinity).
    pub fn identity() -> Self {
        Self {
            x: FieldElement::<P::BaseField>::zero(),
            y: FieldElement::<P::BaseField>::one(),
            infinity: true,
            _marker: PhantomData,
        }
    }

    /// Create the generator point from the curve configuration.
    pub fn generator() -> Self {
        Self {
            x: P::generator_x(),
            y: P::generator_y(),
            infinity: false,
            _marker: PhantomData,
        }
    }

    /// Create a new affine point from coordinates.
    pub fn new(x: FieldElement<P::BaseField>, y: FieldElement<P::BaseField>) -> Self {
        Self {
            x,
            y,
            infinity: false,
            _marker: PhantomData,
        }
    }

    /// Check if this is the identity point.
    pub fn is_identity(&self) -> bool {
        self.infinity
    }

    /// Check if this point lies on the curve.
    pub fn is_on_curve(&self) -> bool {
        if self.infinity {
            return true;
        }
        // y² = x³ + ax + b
        let y2 = self.y * self.y;
        let x2 = self.x * self.x;
        let x3 = x2 * self.x;
        let ax = P::mul_by_a(self.x);
        let rhs = P::add_b(x3 + ax);
        y2 == rhs
    }

    /// Convert to projective coordinates.
    pub fn into_projective(self) -> Projective<P> {
        if self.infinity {
            Projective::identity()
        } else {
            Projective {
                x: self.x,
                y: self.y,
                z: FieldElement::<P::BaseField>::one(),
                _marker: PhantomData,
            }
        }
    }
}

impl<P: ShortWeierstrassConfig> Neg for Affine<P> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        if self.infinity {
            self
        } else {
            Self {
                x: self.x,
                y: -self.y,
                infinity: false,
                _marker: PhantomData,
            }
        }
    }
}

/// A point on a Short Weierstrass curve in projective coordinates (X : Y : Z).
///
/// Generic over `P: ShortWeierstrassConfig` which defines the curve parameters.
#[derive(Clone, Copy, Debug)]
pub struct Projective<P: ShortWeierstrassConfig> {
    /// X coordinate
    pub x: FieldElement<P::BaseField>,
    /// Y coordinate
    pub y: FieldElement<P::BaseField>,
    /// Z coordinate
    pub z: FieldElement<P::BaseField>,
    _marker: PhantomData<P>,
}

impl<P: ShortWeierstrassConfig> PartialEq for Projective<P> {
    fn eq(&self, other: &Self) -> bool {
        if self.is_identity() && other.is_identity() {
            return true;
        }
        if self.is_identity() || other.is_identity() {
            return false;
        }
        // Compare: X1*Z2 == X2*Z1 and Y1*Z2 == Y2*Z1
        let lhs_x = self.x * other.z;
        let rhs_x = other.x * self.z;
        let lhs_y = self.y * other.z;
        let rhs_y = other.y * self.z;
        lhs_x == rhs_x && lhs_y == rhs_y
    }
}

impl<P: ShortWeierstrassConfig> Eq for Projective<P> {}

impl<P: ShortWeierstrassConfig> Projective<P> {
    /// Create the identity (point at infinity).
    pub fn identity() -> Self {
        Self {
            x: FieldElement::<P::BaseField>::one(),
            y: FieldElement::<P::BaseField>::one(),
            z: FieldElement::<P::BaseField>::zero(),
            _marker: PhantomData,
        }
    }

    /// Create the generator point from the curve configuration.
    pub fn generator() -> Self {
        Self {
            x: P::generator_x(),
            y: P::generator_y(),
            z: FieldElement::<P::BaseField>::one(),
            _marker: PhantomData,
        }
    }

    /// Check if this is the identity point.
    pub fn is_identity(&self) -> bool {
        self.z.is_zero()
    }

    /// Convert to affine coordinates.
    pub fn to_affine(&self) -> Affine<P> {
        if self.is_identity() {
            Affine::identity()
        } else {
            let z_inv = Field::inv(&self.z).unwrap();
            Affine::new(self.x * z_inv, self.y * z_inv)
        }
    }

    /// Helper to create small field constants.
    fn from_u64(val: u64) -> FieldElement<P::BaseField> {
        FieldElement::<P::BaseField>::new(U1024::from_u64(val))
    }

    /// Add two projective points.
    pub fn add(&self, rhs: &Self) -> Self {
        if self.is_identity() {
            return rhs.clone();
        }
        if rhs.is_identity() {
            return self.clone();
        }

        let u1 = self.x * rhs.z;
        let u2 = rhs.x * self.z;
        let s1 = self.y * rhs.z;
        let s2 = rhs.y * self.z;

        if u1 == u2 {
            return if s1 == s2 {
                self.double()
            } else {
                Self::identity()
            };
        }

        let h = u2 - u1;
        let r = s2 - s1;
        let hh = h * h;
        let hhh = hh * h;
        let v = u1 * hh;
        let two = Self::from_u64(2);

        let x3 = (r * r) - hhh - (two * v);
        let y3 = (r * (v - x3)) - (s1 * hhh);
        let z3 = self.z * rhs.z * h;

        Self {
            x: x3,
            y: y3,
            z: z3,
            _marker: PhantomData,
        }
    }

    /// Double this point.
    pub fn double(&self) -> Self {
        if self.is_identity() {
            return self.clone();
        }

        let xx = self.x * self.x;
        let yy = self.y * self.y;
        let yyyy = yy * yy;
        let zz = self.z * self.z;

        let two = Self::from_u64(2);
        let three = Self::from_u64(3);
        let eight = Self::from_u64(8);

        let s = two * ((self.x * yy) * two);
        let zzzz = zz * zz;
        let m = (three * xx) + P::mul_by_a(zzzz);

        let x_new = (m * m) - (s * two);
        let z_new = (self.y * self.z) * two;
        let t = eight * yyyy;
        let y_new = (m * (s - x_new)) - t;

        Self {
            x: x_new,
            y: y_new,
            z: z_new,
            _marker: PhantomData,
        }
    }

    /// Scalar multiplication using double-and-add.
    pub fn mul(&self, scalar: &U1024) -> Self {
        let mut result = Self::identity();
        let mut base = self.clone();

        for i in 0..1024 {
            let limb_idx = i / 64;
            let bit_idx = i % 64;
            if (scalar.0[limb_idx] >> bit_idx) & 1 == 1 {
                result = result.add(&base);
            }
            base = base.double();
        }
        result
    }

    /// Negate this point.
    pub fn neg(&self) -> Self {
        Self {
            x: self.x,
            y: -self.y,
            z: self.z,
            _marker: PhantomData,
        }
    }
}

use crate::traits::ProjectivePoint;

impl<P: ShortWeierstrassConfig> ProjectivePoint for Projective<P> {
    type Field = FieldElement<P::BaseField>;

    fn is_identity(&self) -> bool {
        self.is_identity()
    }

    fn add(&self, rhs: &Self) -> Self {
        self.add(rhs)
    }

    fn double(&self) -> Self {
        self.double()
    }

    fn to_affine(&self) -> (Self::Field, Self::Field) {
        let affine = self.to_affine();
        (affine.x, affine.y)
    }

    fn mul(&self, scalar: &U1024) -> Self {
        self.mul(scalar)
    }

    fn neg(&self) -> Self {
        self.neg()
    }
}

impl<P: ShortWeierstrassConfig> Neg for Projective<P> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            x: self.x,
            y: -self.y,
            z: self.z,
            _marker: PhantomData,
        }
    }
}

/// Legacy Weierstrass curve struct (for backward compatibility).
///
/// **Deprecated**: Prefer using `Projective<P>` with a `ShortWeierstrassConfig`.
#[derive(Clone, Debug)]
pub struct WeierstrassCurve<C: FieldConfig> {
    pub a: FieldElement<C>,
    pub b: FieldElement<C>,
    pub generator_x: FieldElement<C>,
    pub generator_y: FieldElement<C>,
    _marker: PhantomData<C>,
}

impl<C: FieldConfig> WeierstrassCurve<C> {
    pub fn new(
        a: FieldElement<C>,
        b: FieldElement<C>,
        generator_x: FieldElement<C>,
        generator_y: FieldElement<C>,
    ) -> Self {
        Self {
            a,
            b,
            generator_x,
            generator_y,
            _marker: PhantomData,
        }
    }
}

/// Legacy point type (for backward compatibility).
#[derive(Clone, Debug)]
pub struct WeierstrassPoint<C: FieldConfig> {
    pub x: FieldElement<C>,
    pub y: FieldElement<C>,
    pub z: FieldElement<C>,
    pub curve: WeierstrassCurve<C>,
}

impl<C: FieldConfig> PartialEq for WeierstrassPoint<C> {
    fn eq(&self, other: &Self) -> bool {
        if self.z.is_zero() && other.z.is_zero() {
            return true;
        }
        if self.z.is_zero() || other.z.is_zero() {
            return false;
        }
        let lhs_x = self.x * other.z;
        let rhs_x = other.x * self.z;
        let lhs_y = self.y * other.z;
        let rhs_y = other.y * self.z;
        lhs_x == rhs_x && lhs_y == rhs_y
    }
}

impl<C: FieldConfig> Eq for WeierstrassPoint<C> {}

impl<C: FieldConfig> Neg for WeierstrassPoint<C> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        WeierstrassPoint {
            x: self.x,
            y: -self.y,
            z: self.z,
            curve: self.curve,
        }
    }
}

impl<C: FieldConfig> WeierstrassPoint<C> {
    /// Create a new point from affine coordinates.
    ///
    /// # Panics
    ///
    /// Panics if the coordinates do not lie on the curve.
    /// For the identity point, use `WeierstrassPoint::identity()` instead.
    pub fn new_affine(x: FieldElement<C>, y: FieldElement<C>, curve: WeierstrassCurve<C>) -> Self {
        // Validate that the point lies on the curve
        assert!(
            curve.is_on_curve(&x, &y),
            "Point ({:?}, {:?}) does not lie on the curve",
            x.to_u1024(),
            y.to_u1024()
        );

        WeierstrassPoint {
            x,
            y,
            z: FieldElement::<C>::one(),
            curve,
        }
    }

    /// Create the identity (point at infinity).
    pub fn identity(curve: WeierstrassCurve<C>) -> Self {
        WeierstrassPoint {
            x: FieldElement::<C>::one(),
            y: FieldElement::<C>::one(),
            z: FieldElement::<C>::zero(),
            curve,
        }
    }

    fn from_u64(val: u64) -> FieldElement<C> {
        FieldElement::<C>::new(U1024::from_u64(val))
    }
}

impl<C: FieldConfig> Curve for WeierstrassCurve<C> {
    type Point = WeierstrassPoint<C>;
    type BaseField = FieldElement<C>;

    fn identity(&self) -> Self::Point {
        WeierstrassPoint {
            x: FieldElement::<C>::one(),
            y: FieldElement::<C>::one(),
            z: FieldElement::<C>::zero(),
            curve: self.clone(),
        }
    }

    fn is_on_curve(&self, x: &Self::BaseField, y: &Self::BaseField) -> bool {
        let y2 = *y * *y;
        let x2 = *x * *x;
        let x3 = x2 * *x;
        let ax = self.a * *x;
        let rhs = x3 + ax + self.b;
        y2 == rhs
    }

    fn generator(&self) -> Self::Point {
        WeierstrassPoint {
            x: self.generator_x,
            y: self.generator_y,
            z: FieldElement::<C>::one(),
            curve: self.clone(),
        }
    }

    fn cofactor(&self) -> &'static [u64] {
        &[1]
    }
}

impl<C: FieldConfig> crate::traits::ProjectivePoint for WeierstrassPoint<C> {
    type Field = FieldElement<C>;

    fn is_identity(&self) -> bool {
        self.z.is_zero()
    }

    fn add(&self, rhs: &Self) -> Self {
        if self.is_identity() {
            return rhs.clone();
        }
        if rhs.is_identity() {
            return self.clone();
        }

        let u1 = self.x * rhs.z;
        let u2 = rhs.x * self.z;
        let s1 = self.y * rhs.z;
        let s2 = rhs.y * self.z;

        if u1 == u2 {
            return if s1 == s2 {
                self.double()
            } else {
                self.curve.identity()
            };
        }

        let h = u2 - u1;
        let r = s2 - s1;
        let hh = h * h;
        let hhh = hh * h;
        let v = u1 * hh;
        let two = Self::from_u64(2);

        let x3 = (r * r) - hhh - (two * v);
        let y3 = (r * (v - x3)) - (s1 * hhh);
        let z3 = self.z * rhs.z * h;

        Self {
            x: x3,
            y: y3,
            z: z3,
            curve: self.curve.clone(),
        }
    }

    fn double(&self) -> Self {
        if self.is_identity() {
            return self.clone();
        }

        let xx = self.x * self.x;
        let yy = self.y * self.y;
        let yyyy = yy * yy;
        let zz = self.z * self.z;

        let two = Self::from_u64(2);
        let three = Self::from_u64(3);
        let eight = Self::from_u64(8);

        let s = two * ((self.x * yy) * two);
        let zzzz = zz * zz;
        let m = (three * xx) + (self.curve.a * zzzz);

        let x_new = (m * m) - (s * two);
        let z_new = (self.y * self.z) * two;
        let t = eight * yyyy;
        let y_new = (m * (s - x_new)) - t;

        Self {
            x: x_new,
            y: y_new,
            z: z_new,
            curve: self.curve.clone(),
        }
    }

    fn to_affine(&self) -> (Self::Field, Self::Field) {
        if self.is_identity() {
            return (FieldElement::<C>::zero(), FieldElement::<C>::zero());
        }
        let z_inv = Field::inv(&self.z).unwrap();
        (self.x * z_inv, self.y * z_inv)
    }

    fn mul(&self, scalar: &U1024) -> Self {
        let mut result = self.curve.identity();
        let mut base = self.clone();
        for i in 0..1024 {
            let limb_idx = i / 64;
            let bit_idx = i % 64;
            if (scalar.0[limb_idx] >> bit_idx) & 1 == 1 {
                result = result.add(&base);
            }
            base = base.double();
        }
        result
    }

    fn neg(&self) -> Self {
        WeierstrassPoint {
            x: self.x,
            y: -self.y,
            z: self.z,
            curve: self.curve.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instances::bls6_6::{Bls6_6BaseField, Bls6_6G1Config};

    #[test]
    fn test_identity() {
        let id = Projective::<Bls6_6G1Config>::identity();
        assert!(id.is_identity());
    }

    #[test]
    fn test_generator() {
        let g = Projective::<Bls6_6G1Config>::generator();
        assert!(!g.is_identity());

        let g_affine = g.to_affine();
        assert!(g_affine.is_on_curve());
    }

    #[test]
    fn test_generator_on_curve() {
        let g = Affine::<Bls6_6G1Config>::generator();
        assert!(g.is_on_curve());
    }

    #[test]
    fn test_legacy_curve() {
        // Test backward compatibility with legacy types
        use crate::instances::bls6_6::get_g1_curve;
        let curve = get_g1_curve();
        let g = curve.generator();
        let (x, y) = g.to_affine();
        assert!(curve.is_on_curve(&x, &y));
    }

    #[test]
    fn test_legacy_new_affine_validation() {
        use crate::instances::bls6_6::get_g1_curve;
        let curve = get_g1_curve();
        let g = curve.generator();
        let (x, y) = g.to_affine();

        // Should succeed for valid point on curve
        let _p = WeierstrassPoint::new_affine(x, y, curve.clone());
    }

    #[test]
    #[should_panic(expected = "does not lie on the curve")]
    fn test_legacy_new_affine_rejects_invalid_point() {
        use crate::instances::bls6_6::get_g1_curve;
        let curve = get_g1_curve();

        // (0, 0) is not on the curve and should panic
        let zero = FieldElement::<Bls6_6BaseField>::zero();
        let _p = WeierstrassPoint::new_affine(zero, zero, curve);
    }

    #[test]
    fn test_legacy_identity_constructor() {
        use crate::instances::bls6_6::get_g1_curve;
        let curve = get_g1_curve();

        // Use explicit identity constructor
        let id = WeierstrassPoint::identity(curve);
        assert!(id.is_identity());
    }
}
