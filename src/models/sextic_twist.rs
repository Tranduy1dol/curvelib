//! Sextic Twist curve implementation for G2 in pairing-friendly curves.
//!
//! This module provides the `SexticTwist` and `TwistPoint` types for
//! working with curves over Fp2 (quadratic extension).

use std::marker::PhantomData;
use std::ops::Neg;

use mathlib::{FieldConfig, FieldElement, U1024};

use crate::algebra::fields::Fp2;
use crate::traits::{Curve, ProjectivePoint};

/// A Sextic Twist curve over Fp2: y² = x³ + ax + b
///
/// Generic over `C: FieldConfig` which defines the base prime field.
/// The curve is defined over Fp2, the quadratic extension of Fp.
#[derive(Clone, Debug)]
pub struct SexticTwist<C: FieldConfig> {
    /// Coefficient a (in Fp2)
    pub a: Fp2<C>,
    /// Coefficient b (in Fp2)
    pub b: Fp2<C>,
    /// Generator x-coordinate (in Fp2)
    pub generator_x: Fp2<C>,
    /// Generator y-coordinate (in Fp2)
    pub generator_y: Fp2<C>,
    _marker: PhantomData<C>,
}

impl<C: FieldConfig> SexticTwist<C> {
    /// Create a new Sextic Twist curve.
    pub fn new(a: Fp2<C>, b: Fp2<C>, generator_x: Fp2<C>, generator_y: Fp2<C>) -> Self {
        Self {
            a,
            b,
            generator_x,
            generator_y,
            _marker: PhantomData,
        }
    }
}

impl<C: FieldConfig> Curve for SexticTwist<C> {
    type Point = TwistPoint<C>;
    type BaseField = Fp2<C>;

    fn identity(&self) -> Self::Point {
        TwistPoint {
            x: Fp2::one(),
            y: Fp2::one(),
            z: Fp2::zero(),
            curve: self.clone(),
        }
    }

    fn is_on_curve(&self, x: &Self::BaseField, y: &Self::BaseField) -> bool {
        // y² = x³ + ax + b (all in Fp2)
        let y2 = *y * *y;
        let x2 = *x * *x;
        let x3 = x2 * *x;
        let ax = self.a * *x;
        let rhs = x3 + ax + self.b;
        y2 == rhs
    }

    fn generator(&self) -> Self::Point {
        TwistPoint {
            x: self.generator_x,
            y: self.generator_y,
            z: Fp2::one(),
            curve: self.clone(),
        }
    }

    fn cofactor(&self) -> &'static [u64] {
        &[1]
    }
}

/// A point on a Sextic Twist curve in projective coordinates over Fp2.
#[derive(Clone, Debug)]
pub struct TwistPoint<C: FieldConfig> {
    /// X coordinate (in Fp2)
    pub x: Fp2<C>,
    /// Y coordinate (in Fp2)
    pub y: Fp2<C>,
    /// Z coordinate (in Fp2)
    pub z: Fp2<C>,
    /// The curve this point belongs to
    pub curve: SexticTwist<C>,
}

impl<C: FieldConfig> PartialEq for TwistPoint<C> {
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

impl<C: FieldConfig> Eq for TwistPoint<C> {}

impl<C: FieldConfig> Neg for TwistPoint<C> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        TwistPoint {
            x: self.x,
            y: -self.y,
            z: self.z,
            curve: self.curve,
        }
    }
}

impl<C: FieldConfig> TwistPoint<C> {
    /// Create a point from affine Fp2 coordinates.
    pub fn new_affine(x: Fp2<C>, y: Fp2<C>, curve: SexticTwist<C>) -> Self {
        if x.is_zero() && y.is_zero() {
            curve.identity()
        } else {
            TwistPoint {
                x,
                y,
                z: Fp2::one(),
                curve,
            }
        }
    }

    /// Helper to create small Fp2 constants from a base field value.
    fn fp2_from_u64(val: u64) -> Fp2<C> {
        let fp_val = FieldElement::<C>::new(U1024::from_u64(val));
        Fp2::new(fp_val, FieldElement::<C>::zero())
    }
}

impl<C: FieldConfig> ProjectivePoint for TwistPoint<C> {
    type Field = Fp2<C>;

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
        let two = Self::fp2_from_u64(2);

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

        let two = Self::fp2_from_u64(2);
        let three = Self::fp2_from_u64(3);
        let eight = Self::fp2_from_u64(8);

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
            return (Fp2::zero(), Fp2::zero());
        }
        let z_inv = self.z.inv().unwrap_or(Fp2::zero());
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
        TwistPoint {
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
    use crate::instances::bls6_6::Bls6_6BaseField;

    #[test]
    fn test_identity() {
        let a = Fp2::<Bls6_6BaseField>::zero();
        let b = Fp2::<Bls6_6BaseField>::one();
        let gx = Fp2::<Bls6_6BaseField>::zero();
        let gy = Fp2::<Bls6_6BaseField>::one();
        let curve = SexticTwist::new(a, b, gx, gy);

        let id = curve.identity();
        assert!(id.is_identity());
    }
}
