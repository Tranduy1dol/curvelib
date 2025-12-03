use mathlib::field::element::FieldElement;
use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{BigInt, U1024};

use crate::traits::ProjectivePoint;

#[derive(Clone, Debug)]
pub struct EdwardsCurve<'a> {
    pub a: FieldElement<'a>,
    pub d: FieldElement<'a>,
    pub params: &'a MontgomeryParams,
}

impl<'a> EdwardsCurve<'a> {
    pub fn new(a: FieldElement<'a>, d: FieldElement<'a>, params: &'a MontgomeryParams) -> Self {
        Self { a, d, params }
    }
}

#[derive(Clone, Debug)]
pub struct TePoint<'a> {
    pub x: FieldElement<'a>,
    pub y: FieldElement<'a>,
    pub z: FieldElement<'a>,
    pub t: FieldElement<'a>,
    pub curve: &'a EdwardsCurve<'a>,
}

impl<'a> PartialEq for TePoint<'a> {
    fn eq(&self, other: &Self) -> bool {
        if self.is_identity() {
            return other.is_identity();
        }
        if other.is_identity() {
            return self.is_identity();
        }

        let x1z2 = self.x * other.z;
        let x2z1 = other.x * self.z;

        let y1z2 = self.y * other.z;
        let y2z1 = other.y * self.z;

        x1z2 == x2z1 && y1z2 == y2z1
    }
}
impl<'a> Eq for TePoint<'a> {}

impl<'a> TePoint<'a> {
    pub fn new_affine(
        x: FieldElement<'a>,
        y: FieldElement<'a>,
        curve: &'a EdwardsCurve<'a>,
    ) -> Self {
        let z = FieldElement::one(curve.params);
        let t = x * y;
        Self { x, y, z, t, curve }
    }

    pub fn identity(curve: &'a EdwardsCurve<'a>) -> Self {
        let zero = FieldElement::zero(curve.params);
        let one = FieldElement::one(curve.params);
        Self {
            x: zero,
            y: one,
            z: one,
            t: zero,
            curve,
        }
    }
}

impl<'a> ProjectivePoint for TePoint<'a> {
    fn identity() -> Self {
        panic!("Use TePoint::identity(curve)");
    }

    fn is_identity(&self) -> bool {
        let zero = FieldElement::zero(self.curve.params);
        self.x == zero && self.y == self.z
    }

    fn to_affine(&self) -> (FieldElement<'a>, FieldElement<'a>) {
        if self.z.value == U1024::zero() {
            let zero = FieldElement::zero(self.curve.params);
            let one = FieldElement::one(self.curve.params);
            return (zero, one);
        }

        let z_inv = self.z.inv();
        let x_aff = self.x * z_inv;
        let y_aff = self.y * z_inv;
        (x_aff, y_aff)
    }

    fn double(&self) -> Self {
        let a = self.x * self.x;
        let b = self.y * self.y;
        let two = FieldElement::new(U1024::from_u64(2), self.curve.params);
        let c = two * (self.z * self.z);

        let d = self.curve.a * a;

        let x_plus_y = self.x + self.y;
        let e = (x_plus_y * x_plus_y) - a - b;

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
            curve: self.curve,
        }
    }

    fn add(&self, rhs: &Self) -> Self {
        let a = self.x * rhs.x;
        let b = self.y * rhs.y;

        let t1t2 = self.t * rhs.t;
        let c = self.curve.d * t1t2;

        let d = self.z * rhs.z;

        let x1_plus_y1 = self.x + self.y;
        let x2_plus_y2 = rhs.x + rhs.y;
        let e = (x1_plus_y1 * x2_plus_y2) - a - b;

        let f = d - c;
        let g = d + c;

        let h = b - (self.curve.a * a);

        let x3 = e * f;
        let y3 = g * h;
        let t3 = e * h;
        let z3 = f * g;

        Self {
            x: x3,
            y: y3,
            z: z3,
            t: t3,
            curve: self.curve,
        }
    }
}
