//! Sextic extension field Fp6 = Fp2\[v\] / (v³ - ξ).
//!
//! This module provides the sextic extension field built as a cubic extension
//! over Fp2, where `v`³ = ξ = u (the imaginary unit in Fp2).

use std::ops::{Add, Mul, Neg, Sub};

use mathlib::FieldConfig;

use crate::algebra::fields::Fp2;
use crate::traits::Field;

/// Sextic extension field element: c0 + c1·`v` + c2·`v`² where `v`³ = ξ.
///
/// `Fp6<C>` represents elements of the sextic extension field over the base
/// prime field defined by `C: FieldConfig`. The extension is constructed as
/// Fp2\[v\]/(v³ - ξ) where ξ = u is the imaginary unit in Fp2.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fp6<C: FieldConfig> {
    /// Coefficient of 1
    pub c0: Fp2<C>,
    /// Coefficient of v
    pub c1: Fp2<C>,
    /// Coefficient of v²
    pub c2: Fp2<C>,
}

impl<C: FieldConfig> Fp6<C> {
    /// Create a new Fp6 element from three Fp2 coefficients.
    pub fn new(c0: Fp2<C>, c1: Fp2<C>, c2: Fp2<C>) -> Self {
        Self { c0, c1, c2 }
    }

    /// Returns the additive identity (zero).
    pub fn zero() -> Self {
        Self {
            c0: Fp2::zero(),
            c1: Fp2::zero(),
            c2: Fp2::zero(),
        }
    }

    /// Returns the multiplicative identity (one).
    pub fn one() -> Self {
        Self {
            c0: Fp2::one(),
            c1: Fp2::zero(),
            c2: Fp2::zero(),
        }
    }

    /// Returns true if this element is zero.
    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero() && self.c2.is_zero()
    }

    /// Computes the square of this element.
    pub fn square(&self) -> Self {
        *self * *self
    }

    /// Computes 2 * self.
    pub fn double(&self) -> Self {
        *self + *self
    }
}

impl<C: FieldConfig> Add for Fp6<C> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            c0: self.c0 + rhs.c0,
            c1: self.c1 + rhs.c1,
            c2: self.c2 + rhs.c2,
        }
    }
}

impl<C: FieldConfig> Sub for Fp6<C> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            c0: self.c0 - rhs.c0,
            c1: self.c1 - rhs.c1,
            c2: self.c2 - rhs.c2,
        }
    }
}

impl<C: FieldConfig> Neg for Fp6<C> {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            c0: -self.c0,
            c1: -self.c1,
            c2: -self.c2,
        }
    }
}

impl<C: FieldConfig> Mul for Fp6<C> {
    type Output = Self;

    /// Multiply two Fp6 field elements.
    ///
    /// The product is computed in the cubic extension Fp6 = Fp2\[v\]/(v³ - ξ)
    /// where `v`³ = ξ = u (the imaginary unit in Fp2).
    fn mul(self, rhs: Self) -> Self {
        let a = self.c0;
        let b = self.c1;
        let c = self.c2;
        let d = rhs.c0;
        let e = rhs.c1;
        let f = rhs.c2;

        // ξ = u = Fp2(0, 1)
        let xi = Fp2::u();

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

/// Implement Field trait for Fp6
impl<C: FieldConfig> Field for Fp6<C> {
    fn zero() -> Self {
        Fp6::zero()
    }

    fn is_zero(&self) -> bool {
        Fp6::is_zero(self)
    }

    fn one() -> Self {
        Fp6::one()
    }

    fn inv(&self) -> Option<Self> {
        // Fp6 inversion is complex; not implemented yet
        None
    }

    fn double(&self) -> Self {
        Fp6::double(self)
    }

    fn mul(&self, rhs: &Self) -> Self {
        *self * *rhs
    }

    fn add(&self, rhs: &Self) -> Self {
        *self + *rhs
    }

    fn square(&self) -> Self {
        Fp6::square(self)
    }
}
