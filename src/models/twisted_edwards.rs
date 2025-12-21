//! Twisted Edwards curve implementation: ax² + y² = 1 + dx²y²
//!
//! This module provides projective and affine point types for
//! working with elliptic curves in Twisted Edwards form.
//!
//! # Design
//!
//! Following the arkworks pattern, curve parameters are defined at the type level
//! via `TwistedEdwardsConfig`. Points are generic over the config type, not
//! over runtime curve instances.

use std::marker::PhantomData;
use std::ops::Neg;

use mathlib::{FieldConfig, FieldElement, U1024};

use crate::traits::{Curve, Field, TwistedEdwardsConfig};

/// A point on a Twisted Edwards curve in affine coordinates (x, y).
///
/// Generic over `P: TwistedEdwardsConfig` which defines the curve parameters.
#[derive(Clone, Copy, Debug)]
pub struct Affine<P: TwistedEdwardsConfig> {
    /// X coordinate
    pub x: FieldElement<P::BaseField>,
    /// Y coordinate
    pub y: FieldElement<P::BaseField>,
    _marker: PhantomData<P>,
}

impl<P: TwistedEdwardsConfig> PartialEq for Affine<P> {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl<P: TwistedEdwardsConfig> Eq for Affine<P> {}

impl<P: TwistedEdwardsConfig> Affine<P> {
    /// Create the identity point (0, 1) on Edwards curves.
    pub fn identity() -> Self {
        Self {
            x: FieldElement::<P::BaseField>::zero(),
            y: FieldElement::<P::BaseField>::one(),
            _marker: PhantomData,
        }
    }

    /// Create the generator point from the curve configuration.
    pub fn generator() -> Self {
        Self {
            x: P::generator_x(),
            y: P::generator_y(),
            _marker: PhantomData,
        }
    }

    /// Create a new affine point from coordinates.
    pub fn new(x: FieldElement<P::BaseField>, y: FieldElement<P::BaseField>) -> Self {
        Self {
            x,
            y,
            _marker: PhantomData,
        }
    }

    /// Check if this is the identity point.
    pub fn is_identity(&self) -> bool {
        self.x.is_zero() && self.y == FieldElement::<P::BaseField>::one()
    }

    /// Check if this point lies on the curve.
    pub fn is_on_curve(&self) -> bool {
        // ax² + y² = 1 + dx²y²
        let x2 = self.x * self.x;
        let y2 = self.y * self.y;
        let lhs = P::mul_by_a(x2) + y2;
        let rhs = FieldElement::<P::BaseField>::one() + P::coeff_d() * x2 * y2;
        lhs == rhs
    }

    /// Convert to extended projective coordinates.
    pub fn into_projective(self) -> Projective<P> {
        Projective {
            x: self.x,
            y: self.y,
            z: FieldElement::<P::BaseField>::one(),
            t: self.x * self.y,
            _marker: PhantomData,
        }
    }
}

impl<P: TwistedEdwardsConfig> Neg for Affine<P> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: self.y,
            _marker: PhantomData,
        }
    }
}

/// A point on a Twisted Edwards curve in extended projective coordinates (X : Y : Z : T).
///
/// The extended coordinate T = X * Y / Z satisfies the relation.
/// Generic over `P: TwistedEdwardsConfig` which defines the curve parameters.
#[derive(Clone, Copy, Debug)]
pub struct Projective<P: TwistedEdwardsConfig> {
    /// X coordinate
    pub x: FieldElement<P::BaseField>,
    /// Y coordinate
    pub y: FieldElement<P::BaseField>,
    /// Z coordinate
    pub z: FieldElement<P::BaseField>,
    /// Extended coordinate T = X * Y / Z
    pub t: FieldElement<P::BaseField>,
    _marker: PhantomData<P>,
}

impl<P: TwistedEdwardsConfig> PartialEq for Projective<P> {
    fn eq(&self, other: &Self) -> bool {
        // Compare in projective coordinates
        let lhs_x = self.x * other.z;
        let rhs_x = other.x * self.z;
        let lhs_y = self.y * other.z;
        let rhs_y = other.y * self.z;
        lhs_x == rhs_x && lhs_y == rhs_y
    }
}

impl<P: TwistedEdwardsConfig> Eq for Projective<P> {}

impl<P: TwistedEdwardsConfig> Projective<P> {
    /// Create the identity point (0, 1, 1, 0) on Edwards curves.
    pub fn identity() -> Self {
        Self {
            x: FieldElement::<P::BaseField>::zero(),
            y: FieldElement::<P::BaseField>::one(),
            z: FieldElement::<P::BaseField>::one(),
            t: FieldElement::<P::BaseField>::zero(),
            _marker: PhantomData,
        }
    }

    /// Create the generator point from the curve configuration.
    pub fn generator() -> Self {
        let x = P::generator_x();
        let y = P::generator_y();
        Self {
            x,
            y,
            z: FieldElement::<P::BaseField>::one(),
            t: x * y,
            _marker: PhantomData,
        }
    }

    /// Check if this is the identity point.
    pub fn is_identity(&self) -> bool {
        self.x.is_zero() && self.y == self.z
    }

    /// Convert to affine coordinates.
    pub fn to_affine(&self) -> Affine<P> {
        if self.z.is_zero() {
            return Affine::identity();
        }
        let z_inv = Field::inv(&self.z).unwrap();
        Affine::new(self.x * z_inv, self.y * z_inv)
    }

    /// Helper to create small field constants.
    fn from_u64(val: u64) -> FieldElement<P::BaseField> {
        FieldElement::<P::BaseField>::new(U1024::from_u64(val))
    }

    /// Add two extended projective points.
    pub fn add(&self, rhs: &Self) -> Self {
        // Extended Edwards addition formula
        let a = self.x * rhs.x;
        let b = self.y * rhs.y;
        let c = self.t * P::coeff_d() * rhs.t;
        let d = self.z * rhs.z;
        let e = (self.x + self.y) * (rhs.x + rhs.y) - a - b;
        let f = d - c;
        let g = d + c;
        let h = b - P::mul_by_a(a);

        let x3 = e * f;
        let y3 = g * h;
        let t3 = e * h;
        let z3 = f * g;

        Self {
            x: x3,
            y: y3,
            z: z3,
            t: t3,
            _marker: PhantomData,
        }
    }

    /// Double this point.
    pub fn double(&self) -> Self {
        // Extended Edwards doubling formula
        let a = self.x * self.x;
        let b = self.y * self.y;
        let two = Self::from_u64(2);
        let c = two * self.z * self.z;
        let d = P::mul_by_a(a);
        let e = (self.x + self.y) * (self.x + self.y) - a - b;
        let g = d + b;
        let f = g - c;
        let h = d - b;

        let x3 = e * f;
        let y3 = g * h;
        let t3 = e * h;
        let z3 = f * g;

        Self {
            x: x3,
            y: y3,
            z: z3,
            t: t3,
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
            x: -self.x,
            y: self.y,
            z: self.z,
            t: -self.t,
            _marker: PhantomData,
        }
    }
}

use crate::traits::ProjectivePoint;

impl<P: TwistedEdwardsConfig> ProjectivePoint for Projective<P> {
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

impl<P: TwistedEdwardsConfig> Neg for Projective<P> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: self.y,
            z: self.z,
            t: -self.t,
            _marker: PhantomData,
        }
    }
}

/// Legacy Edwards curve struct (for backward compatibility).
#[derive(Clone, Debug)]
pub struct EdwardsCurve<C: FieldConfig> {
    pub a: FieldElement<C>,
    pub d: FieldElement<C>,
    pub generator_x: FieldElement<C>,
    pub generator_y: FieldElement<C>,
    _marker: PhantomData<C>,
}

impl<C: FieldConfig> EdwardsCurve<C> {
    pub fn new(
        a: FieldElement<C>,
        d: FieldElement<C>,
        generator_x: FieldElement<C>,
        generator_y: FieldElement<C>,
    ) -> Self {
        Self {
            a,
            d,
            generator_x,
            generator_y,
            _marker: PhantomData,
        }
    }
}

/// Legacy point type (for backward compatibility).
#[derive(Clone, Debug)]
pub struct EdwardsPoint<C: FieldConfig> {
    pub x: FieldElement<C>,
    pub y: FieldElement<C>,
    pub z: FieldElement<C>,
    pub t: FieldElement<C>,
    pub curve: EdwardsCurve<C>,
}

impl<C: FieldConfig> PartialEq for EdwardsPoint<C> {
    fn eq(&self, other: &Self) -> bool {
        let lhs_x = self.x * other.z;
        let rhs_x = other.x * self.z;
        let lhs_y = self.y * other.z;
        let rhs_y = other.y * self.z;
        lhs_x == rhs_x && lhs_y == rhs_y
    }
}

impl<C: FieldConfig> Eq for EdwardsPoint<C> {}

impl<C: FieldConfig> Neg for EdwardsPoint<C> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        EdwardsPoint {
            x: -self.x,
            y: self.y,
            z: self.z,
            t: -self.t,
            curve: self.curve,
        }
    }
}

impl<C: FieldConfig> EdwardsPoint<C> {
    pub fn new_affine(x: FieldElement<C>, y: FieldElement<C>, curve: EdwardsCurve<C>) -> Self {
        EdwardsPoint {
            x,
            y,
            z: FieldElement::<C>::one(),
            t: x * y,
            curve,
        }
    }

    fn from_u64(val: u64) -> FieldElement<C> {
        FieldElement::<C>::new(U1024::from_u64(val))
    }
}

impl<C: FieldConfig> Curve for EdwardsCurve<C> {
    type Point = EdwardsPoint<C>;
    type BaseField = FieldElement<C>;

    fn identity(&self) -> Self::Point {
        EdwardsPoint {
            x: FieldElement::<C>::zero(),
            y: FieldElement::<C>::one(),
            z: FieldElement::<C>::one(),
            t: FieldElement::<C>::zero(),
            curve: self.clone(),
        }
    }

    fn is_on_curve(&self, x: &Self::BaseField, y: &Self::BaseField) -> bool {
        let x2 = *x * *x;
        let y2 = *y * *y;
        let lhs = self.a * x2 + y2;
        let rhs = FieldElement::<C>::one() + self.d * x2 * y2;
        lhs == rhs
    }

    fn generator(&self) -> Self::Point {
        EdwardsPoint {
            x: self.generator_x,
            y: self.generator_y,
            z: FieldElement::<C>::one(),
            t: self.generator_x * self.generator_y,
            curve: self.clone(),
        }
    }

    fn cofactor(&self) -> &'static [u64] {
        &[1]
    }
}

impl<C: FieldConfig> crate::traits::ProjectivePoint for EdwardsPoint<C> {
    type Field = FieldElement<C>;

    fn is_identity(&self) -> bool {
        self.x.is_zero() && self.y == self.z
    }

    fn add(&self, rhs: &Self) -> Self {
        let a = self.x * rhs.x;
        let b = self.y * rhs.y;
        let c = self.t * self.curve.d * rhs.t;
        let d = self.z * rhs.z;
        let e = (self.x + self.y) * (rhs.x + rhs.y) - a - b;
        let f = d - c;
        let g = d + c;
        let h = b - self.curve.a * a;

        EdwardsPoint {
            x: e * f,
            y: g * h,
            z: f * g,
            t: e * h,
            curve: self.curve.clone(),
        }
    }

    fn double(&self) -> Self {
        let a = self.x * self.x;
        let b = self.y * self.y;
        let two = Self::from_u64(2);
        let c = two * self.z * self.z;
        let d = self.curve.a * a;
        let e = (self.x + self.y) * (self.x + self.y) - a - b;
        let g = d + b;
        let f = g - c;
        let h = d - b;

        EdwardsPoint {
            x: e * f,
            y: g * h,
            z: f * g,
            t: e * h,
            curve: self.curve.clone(),
        }
    }

    fn to_affine(&self) -> (Self::Field, Self::Field) {
        if self.z.is_zero() {
            return (FieldElement::<C>::zero(), FieldElement::<C>::one());
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
        EdwardsPoint {
            x: -self.x,
            y: self.y,
            z: self.z,
            t: -self.t,
            curve: self.curve.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instances::tiny_jubjub::TinyJubjubConfig;

    #[test]
    fn test_identity() {
        let id = Projective::<TinyJubjubConfig>::identity();
        assert!(id.is_identity());
    }

    #[test]
    fn test_generator() {
        let g = Projective::<TinyJubjubConfig>::generator();
        assert!(!g.is_identity());

        let g_affine = g.to_affine();
        assert!(g_affine.is_on_curve());
    }

    #[test]
    fn test_generator_on_curve() {
        let g = Affine::<TinyJubjubConfig>::generator();
        assert!(g.is_on_curve());
    }
}
