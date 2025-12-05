use mathlib::field::{element::FieldElement, montgomery::MontgomeryParams};
use mathlib::{BigInt, U1024};

use crate::traits::{Curve, ProjectivePoint};

#[derive(Clone, Debug)]
pub struct WeierstrassCurve<'a> {
    pub a: FieldElement<'a>,
    pub b: FieldElement<'a>,
    pub params: &'a MontgomeryParams,
}

impl<'a> WeierstrassCurve<'a> {
    pub fn new(a: FieldElement<'a>, b: FieldElement<'a>, params: &'a MontgomeryParams) -> Self {
        Self { a, b, params }
    }
}

impl<'a> Curve<'a> for WeierstrassCurve<'a> {
    type Point = SWPoint<'a>;

    fn identity(&self) -> Self::Point {
        let curve = self.clone();
        let one = FieldElement::one(curve.params);
        let zero = FieldElement::zero(curve.params);
        SWPoint {
            x: one,
            y: one,
            z: zero,
            curve,
        }
    }

    fn is_on_curve(&self, x: &FieldElement, y: &FieldElement) -> bool {
        let y2 = *y * *y;
        let x2 = *x * *x;
        let x3 = x2 * *x;
        let ax = self.a * *x;
        let rhs = x3 + ax + self.b;
        y2 == rhs
    }
}

#[derive(Clone, Debug)]
pub struct SWPoint<'a> {
    pub x: FieldElement<'a>,
    pub y: FieldElement<'a>,
    pub z: FieldElement<'a>,
    pub curve: WeierstrassCurve<'a>,
}

impl<'a> PartialEq for SWPoint<'a> {
    fn eq(&self, other: &Self) -> bool {
        if self.is_identity() {
            return other.is_identity();
        }
        if other.is_identity() {
            return self.is_identity();
        }

        let (x1, y1) = self.to_affine();
        let (x2, y2) = other.to_affine();
        x1 == x2 && y1 == y2
    }
}

impl<'a> Eq for SWPoint<'a> {}

impl<'a> SWPoint<'a> {
    pub fn new_affine(
        x: FieldElement<'a>,
        y: FieldElement<'a>,
        curve: &'a WeierstrassCurve<'a>,
    ) -> Self {
        if x.value == U1024::zero() && y.value == U1024::zero() {
            return curve.identity();
        }
        let z = FieldElement::one(curve.params);
        Self {
            x,
            y,
            z,
            curve: curve.clone(),
        }
    }
}

impl<'a> ProjectivePoint for SWPoint<'a> {
    fn is_identity(&self) -> bool {
        self.z.value == U1024::zero()
    }

    fn add(&self, rhs: &Self) -> Self {
        if self.is_identity() {
            return rhs.clone();
        }
        if rhs.is_identity() {
            return self.clone();
        }

        let z1z1 = self.z * self.z;
        let z2z2 = rhs.z * rhs.z;

        let u1 = self.x * z2z2;
        let u2 = rhs.x * z1z1;

        let s1 = self.y * (rhs.z * z2z2);
        let s2 = rhs.y * (self.z * z1z1);

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
        let two = FieldElement::new(U1024::from_u64(2), self.curve.params);

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

        let two = FieldElement::new(U1024::from_u64(2), self.curve.params);
        let s = two * ((self.x * yy) * two);

        let three = FieldElement::new(U1024::from_u64(3), self.curve.params);
        let zzzz = zz * zz;
        let m = (three * xx) + (self.curve.a * zzzz);

        let x_new = (m * m) - (s * two);

        let z_new = (self.y * self.z) * two;

        let eight = FieldElement::new(U1024::from_u64(8), self.curve.params);
        let t = eight * yyyy;
        let y_new = (m * (s - x_new)) - t;

        Self {
            x: x_new,
            y: y_new,
            z: z_new,
            curve: self.curve.clone(),
        }
    }

    fn to_affine(&self) -> (FieldElement<'a>, FieldElement<'a>) {
        if self.is_identity() {
            let zero = FieldElement::zero(self.curve.params);
            return (zero, zero);
        }

        let z_inv = self.z.inv();
        let z2_inv = z_inv * z_inv;
        let z3_inv = z2_inv * z_inv;

        let x_aff = self.x * z2_inv;
        let y_aff = self.y * z3_inv;

        (x_aff, y_aff)
    }

    fn mul(&self, scalar: &U1024) -> Self {
        let mut res = self.curve.identity();

        let num_bits = scalar.bits();
        for i in (0..num_bits).rev() {
            res = res.double();
            if scalar.bit(i) {
                res = res.add(self);
            }
        }
        res
    }
}
