use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{BigInt, FieldElement};

use crate::algebra::fields::{Fp, Fp2};
use crate::def_weierstrass_curve;
use crate::traits::{Curve, Field, ProjectivePoint};

def_weierstrass_curve!(SexticTwist, Fp2<'a>);

impl<'a> Curve<'a> for SexticTwist<'a> {
    type Point = STPoint<'a>;

    fn identity(&self) -> Self::Point {
        let zero_fp = FieldElement::zero(self.params);
        let one_fp = FieldElement::one(self.params);

        let zero_fp2 = Fp2::new(Fp::from(zero_fp), Fp::from(zero_fp));
        let one_fp2 = Fp2::new(Fp::from(one_fp), Fp::from(zero_fp)); // 1 + 0u

        STPoint {
            x: one_fp2,
            y: one_fp2,
            z: zero_fp2,
            curve: self.clone(),
        }
    }

    fn is_on_curve(
        &self,
        x: &<Self::Point as ProjectivePoint<'a>>::Field,
        y: &<Self::Point as ProjectivePoint<'a>>::Field,
    ) -> bool {
        // Check affine: Y^2 = X^3 + aX + b (calculated in Fp2)
        let y2 = *y * *y;
        let x2 = *x * *x;
        let x3 = x2 * *x;
        let ax = self.a * *x;

        let rhs = x3 + ax + self.b;
        y2 == rhs
    }

    fn scalar_params(&self) -> &'a MontgomeryParams {
        self.scalar_params
    }

    fn generator(&self) -> Self::Point {
        let x = self.generator_x;
        let y = self.generator_y;
        let _zero_fp = FieldElement::zero(self.params);
        let _zero_fp2 = Fp2::new(Fp::zero(self.params), Fp::zero(self.params));
        let z = Fp2::new(
            Fp::from(FieldElement::one(self.params)),
            Fp::from(FieldElement::zero(self.params)),
        );

        STPoint {
            x,
            y,
            z,
            curve: self.clone(),
        }
    }
}

// --- 2. G2 POINT DEFINITION (Jacobian Coordinates) ---
#[derive(Clone, Debug)]
pub struct STPoint<'a> {
    pub x: Fp2<'a>,
    pub y: Fp2<'a>,
    pub z: Fp2<'a>,
    pub curve: SexticTwist<'a>,
}

impl<'a> PartialEq for STPoint<'a> {
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
impl<'a> Eq for STPoint<'a> {}

impl<'a> STPoint<'a> {
    pub fn new(x: Fp2<'a>, y: Fp2<'a>, z: Fp2<'a>, curve: SexticTwist<'a>) -> Self {
        Self { x, y, z, curve }
    }

    pub fn neg(&self) -> Self {
        if self.is_identity() {
            return self.clone();
        }
        // -P = (X, -Y, Z)
        let neg_y = -self.y; // Fp2 implements Neg
        Self {
            x: self.x,
            y: neg_y,
            z: self.z,
            curve: self.curve.clone(),
        }
    }
}

impl<'a> ProjectivePoint<'a> for STPoint<'a> {
    type Field = Fp2<'a>;

    fn is_identity(&self) -> bool {
        // Z == 0 in Fp2
        self.z.c0.value == mathlib::U1024::zero() && self.z.c1.value == mathlib::U1024::zero()
    }

    fn add(&self, rhs: &Self) -> Self {
        if self.is_identity() {
            return rhs.clone();
        }
        if rhs.is_identity() {
            return self.clone();
        }

        let z1z1 = self.z.square();
        let z2z2 = rhs.z.square();

        let u1 = self.x * z2z2;
        let u2 = rhs.x * z1z1;

        let s1 = self.y * (rhs.z * z2z2); // Y1 * Z2^3
        let s2 = rhs.y * (self.z * z1z1); // Y2 * Z1^3

        if u1 == u2 {
            return if s1 == s2 {
                self.double()
            } else {
                self.curve.identity()
            };
        }

        let h = u2 - u1;
        let r = s2 - s1;
        let hh = h.square();
        let hhh = hh * h;
        let v = u1 * hh;

        let two_val = mathlib::U1024::from_u64(2);
        let zero_fp = FieldElement::zero(self.curve.params);
        let two_fp = FieldElement::new(two_val, self.curve.params);
        let two_fp2 = Fp2::new(Fp::from(two_fp), Fp::from(zero_fp));

        let x3 = r.square() - hhh - (two_fp2 * v);
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

        // Jacobian doubling formulas (same as SWPoint but on Fp2)
        let xx = self.x.square();
        let yy = self.y.square();
        let yyyy = yy.square();
        let zz = self.z.square();

        // 2
        let zero_fp = FieldElement::zero(self.curve.params);
        let two_val = mathlib::U1024::from_u64(2);
        let two_fp = FieldElement::new(two_val, self.curve.params);
        let two_fp2 = Fp2::new(Fp::from(two_fp), Fp::from(zero_fp));

        // S = 2 * ((X * YY) * 2) = 4XY^2
        let s = two_fp2 * ((self.x * yy) * two_fp2);

        // M = 3*X^2 + a*Z^4
        let three_val = mathlib::U1024::from_u64(3);
        let three_fp = FieldElement::new(three_val, self.curve.params);
        let three_fp2 = Fp2::new(Fp::from(three_fp), Fp::from(zero_fp));

        let zzzz = zz.square();
        let m = (three_fp2 * xx) + (self.curve.a * zzzz);

        // X' = M^2 - 2S
        let x_new = m.square() - (s * two_fp2);

        // Z' = (Y + Z)^2 - YY - ZZ = 2YZ
        // Or Z' = 2 * Y * Z
        let z_new = (self.y * self.z) * two_fp2;

        // Y' = M(S - X') - 8YYYY
        let eight_val = mathlib::U1024::from_u64(8);
        let eight_fp = FieldElement::new(eight_val, self.curve.params);
        let eight_fp2 = Fp2::new(Fp::from(eight_fp), Fp::from(zero_fp));

        let t = eight_fp2 * yyyy;
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
            let zero = FieldElement::zero(self.curve.params);
            let zero_fp2 = Fp2::new(Fp::from(zero), Fp::from(zero));
            return (zero_fp2, zero_fp2);
        }

        let z_inv = self.z.inv().unwrap(); // Fp2 invert
        let z2_inv = z_inv.square();
        let z3_inv = z2_inv * z_inv;

        let x_aff = self.x * z2_inv;
        let y_aff = self.y * z3_inv;
        (x_aff, y_aff)
    }

    fn mul(&self, scalar: &mathlib::U1024) -> Self {
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
