//! Quadratic extension field Fp2 = Fp\[u\] / (u² + 1).
//!
//! This module provides the quadratic extension field where `u`² = -1.

use std::ops::{Add, Mul, Neg, Sub};

use mathlib::{FieldConfig, FieldElement};

use crate::traits::Field;

/// Quadratic extension field element: a + b·u where u² = -1.
///
/// `Fp2<C>` represents elements of the quadratic extension of the prime field
/// defined by `C: FieldConfig`. The extension is formed by adjoining a root
/// of x² + 1 = 0 to the base field.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fp2<C: FieldConfig> {
    /// Real part (coefficient of 1)
    pub c0: FieldElement<C>,
    /// Imaginary part (coefficient of u)
    pub c1: FieldElement<C>,
}

impl<C: FieldConfig> Fp2<C> {
    /// Create a new Fp2 element from two base field elements.
    pub fn new(c0: FieldElement<C>, c1: FieldElement<C>) -> Self {
        Self { c0, c1 }
    }

    /// Returns the additive identity (zero).
    pub fn zero() -> Self {
        Self {
            c0: FieldElement::<C>::zero(),
            c1: FieldElement::<C>::zero(),
        }
    }

    /// Returns the multiplicative identity (one).
    pub fn one() -> Self {
        Self {
            c0: FieldElement::<C>::one(),
            c1: FieldElement::<C>::zero(),
        }
    }

    /// Returns the imaginary unit u (where u² = -1).
    pub fn u() -> Self {
        Self {
            c0: FieldElement::<C>::zero(),
            c1: FieldElement::<C>::one(),
        }
    }

    /// Returns true if this element is zero.
    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    /// Computes the square of this element.
    pub fn square(&self) -> Self {
        *self * *self
    }

    /// Computes the multiplicative inverse, if it exists.
    pub fn inv(&self) -> Option<Self> {
        // For (a + bu), inverse is (a - bu) / (a² + b²)
        // since (a + bu)(a - bu) = a² - b²u² = a² + b² (as u² = -1)
        let a_sq = self.c0.square();
        let b_sq = self.c1.square();
        let norm = a_sq + b_sq;

        let inv_norm = if norm.is_zero() {
            return None;
        } else {
            norm.inv()
        };

        Some(Self {
            c0: self.c0 * inv_norm,
            c1: -self.c1 * inv_norm,
        })
    }

    /// Computes 2 * self.
    pub fn double(&self) -> Self {
        *self + *self
    }
}

impl<C: FieldConfig> Add for Fp2<C> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            c0: self.c0 + rhs.c0,
            c1: self.c1 + rhs.c1,
        }
    }
}

impl<C: FieldConfig> Sub for Fp2<C> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            c0: self.c0 - rhs.c0,
            c1: self.c1 - rhs.c1,
        }
    }
}

impl<C: FieldConfig> Neg for Fp2<C> {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            c0: -self.c0,
            c1: -self.c1,
        }
    }
}

impl<C: FieldConfig> Mul for Fp2<C> {
    type Output = Self;

    /// Multiplies two elements using Karatsuba-like method.
    ///
    /// (a + bu)(c + du) = (ac - bd) + (ad + bc)u
    /// where u² = -1
    fn mul(self, rhs: Self) -> Self {
        // Karatsuba: 3 multiplications instead of 4
        let v0 = self.c0 * rhs.c0; // a * c
        let v1 = self.c1 * rhs.c1; // b * d
        let v2 = (self.c0 + self.c1) * (rhs.c0 + rhs.c1); // (a+b)(c+d)

        Self {
            c0: v0 - v1,      // ac - bd (since u² = -1)
            c1: v2 - v0 - v1, // ad + bc
        }
    }
}

/// Implement Field trait for Fp2
impl<C: FieldConfig> Field for Fp2<C> {
    fn zero() -> Self {
        Fp2::zero()
    }

    fn is_zero(&self) -> bool {
        Fp2::is_zero(self)
    }

    fn one() -> Self {
        Fp2::one()
    }

    fn inv(&self) -> Option<Self> {
        Fp2::inv(self)
    }

    fn double(&self) -> Self {
        Fp2::double(self)
    }

    fn mul(&self, rhs: &Self) -> Self {
        *self * *rhs
    }

    fn add(&self, rhs: &Self) -> Self {
        *self + *rhs
    }

    fn square(&self) -> Self {
        Fp2::square(self)
    }
}
