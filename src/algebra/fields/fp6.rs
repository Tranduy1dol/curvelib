use std::ops::Mul;

use mathlib::field::montgomery::MontgomeryParams;

use crate::algebra::fields::{Fp, Fp2};
use crate::def_fp6;
use crate::traits::Field;

def_fp6!(Fp6, Fp2<'a>);

// Implement Add, Sub (Vector addition)
// Add, Sub, Neg implemented by macro

impl<'a> Mul for Fp6<'a> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        // Cubic extension field multiplication: (a + bv + cv²)(d + ev + fv²)
        // Using Karatsuba-like optimization for 6 -> 5 multiplications
        // Result = (ad) + (ae+bd)v + (af+be+cd)v² + (bf+ce)v³ + (cf)v⁴
        // With v³ = ξ (non-residue), we need: v³ = ξ, v⁴ = ξv, v⁵ = ξv²

        // Generic implementation assuming ξ = 1 (identity for the non-residue)
        // For specific curves, this would need to be parameterized

        let a = self.c0;
        let b = self.c1;
        let c = self.c2;
        let d = rhs.c0;
        let e = rhs.c1;
        let f = rhs.c2;

        // Schoolbook multiplication approach - simpler and generic
        // (a + bv + cv²)(d + ev + fv²) = ad + (ae+bd)v + (af+be+cd)v² + (bf+ce)v³ + (cf)v⁴
        // Assuming v³ = ξ = 1: v³ -> 1, v⁴ -> v, v⁵ -> v²

        let params = self.c0.c0.params;
        let zero = Fp::zero(params);
        let one = Fp::one(params);
        // xi = u = 0 + 1u
        let xi = Fp2::new(zero, one);

        let ad = a * d;
        let ae = a * e;
        let af = a * f;
        let bd = b * d;
        let be = b * e;
        let bf = b * f;
        let cd = c * d;
        let ce = c * e;
        let cf = c * f;

        Self {
            c0: ad + xi * (bf + ce), // Constant term: ad + ξ(bf+ce)
            c1: ae + bd + xi * cf,   // v coefficient: ae + bd + ξ(cf)
            c2: af + be + cd,        // v² coefficient: af + be + cd
        }
    }
}

impl<'a> Field<'a> for Fp6<'a> {
    fn zero(params: &'a MontgomeryParams) -> Self {
        let z = Fp2::zero(params);
        Self {
            c0: z,
            c1: z,
            c2: z,
        }
    }
    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero() && self.c2.is_zero()
    }
    fn one(params: &'a MontgomeryParams) -> Self {
        let z = Fp2::zero(params);
        let o = Fp2::one(params);
        Self {
            c0: o,
            c1: z,
            c2: z,
        }
    }

    fn inv(&self) -> Option<Self> {
        None
    }
    fn double(&self) -> Self {
        *self + *self
    }
    fn mul(&self, rhs: &Self) -> Self {
        *self * *rhs
    }
    fn add(&self, rhs: &Self) -> Self {
        *self + *rhs
    }

    fn square(&self) -> Self {
        *self * *self
    }
}
