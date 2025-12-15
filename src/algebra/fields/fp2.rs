use std::ops::Mul;

use mathlib::field::montgomery::MontgomeryParams;

use crate::algebra::fields::Fp;
use crate::def_fp2;
use crate::traits::Field;

def_fp2!(Fp2, Fp<'a>);

// Macro implements add, Sub, Neg

impl<'a> Mul for Fp2<'a> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        // Quadratic extension field multiplication: (a + bu)(c + du)
        // Using Karatsuba: need only 3 multiplications instead of 4
        // v0 = a*c, v1 = b*d, v2 = (a+b)(c+d)
        // Result: (v0 + β*v1) + (v2 - v0 - v1)u
        // where β is the quadratic non-residue

        let v0 = self.c0 * rhs.c0; // a * c
        let v1 = self.c1 * rhs.c1; // b * d  
        let v2 = (self.c0 + self.c1) * (rhs.c0 + rhs.c1); // (a+b)(c+d)

        Self {
            c0: v0 - v1,      // Real part: ac + β·bd where β = -1
            c1: v2 - v0 - v1, // Imaginary part
        }
    }
}

// Neg implemented by macro

// Trait Field
impl<'a> Field<'a> for Fp2<'a> {
    fn zero(params: &'a MontgomeryParams) -> Self {
        let z = Fp::zero(params);
        Self { c0: z, c1: z }
    }
    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    fn one(params: &'a MontgomeryParams) -> Self {
        let z = Fp::zero(params);
        let o = Fp::one(params);
        Self { c0: o, c1: z }
    }

    fn inv(&self) -> Option<Self> {
        // For quadratic extension (a + bu), inverse is computed as:
        // (a + bu)^-1 = (a' - bu') / norm
        // where norm must match our multiplication law
        // For β=1: norm = a^2 - b^2

        let a = self.c0;
        let b = self.c1;

        let a_sq = a.square();
        let b_sq = b.square();
        let norm = a_sq + b_sq; // For β = -1: norm = a² - β·b² = a² + b²

        let inv_norm = Field::inv(&norm)?; // Return None if norm = 0

        let z = Fp::zero(b.params);
        let neg_b = z - b;

        Some(Self {
            c0: a * inv_norm,
            c1: neg_b * inv_norm,
        })
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
