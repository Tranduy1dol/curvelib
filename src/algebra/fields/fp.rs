use std::ops::{Add, Deref, Mul, Neg, Sub};

use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{FieldElement, U1024};

use crate::traits::curve::ToU1024;
use crate::traits::field::Field;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fp<'a>(pub FieldElement<'a>);

impl<'a> Fp<'a> {
    pub fn new(value: U1024, params: &'a MontgomeryParams) -> Self {
        Self(FieldElement::new(value, params))
    }
}

impl<'a> From<FieldElement<'a>> for Fp<'a> {
    fn from(f: FieldElement<'a>) -> Self {
        Self(f)
    }
}

// Deref to a field element for convenience (accessing .params, .value)
impl<'a> Deref for Fp<'a> {
    type Target = FieldElement<'a>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<'a> Add for Fp<'a> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl<'a> Sub for Fp<'a> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0)
    }
}

impl<'a> Mul for Fp<'a> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self(self.0 * rhs.0)
    }
}

impl<'a> Neg for Fp<'a> {
    type Output = Self;
    fn neg(self) -> Self {
        let zero = FieldElement::zero(self.0.params);
        Self(zero - self.0)
    }
}

impl<'a> Field<'a> for Fp<'a> {
    fn zero(params: &'a MontgomeryParams) -> Self {
        Self(FieldElement::zero(params))
    }

    fn is_zero(&self) -> bool {
        self.0 == FieldElement::zero(self.0.params)
    }

    fn one(params: &'a MontgomeryParams) -> Self {
        Self(FieldElement::one(params))
    }

    fn inv(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(Self(self.0.inv()))
        }
    }

    fn double(&self) -> Self {
        Self(self.0 + self.0)
    }

    fn mul(&self, rhs: &Self) -> Self {
        Self(self.0 * rhs.0)
    }

    fn add(&self, rhs: &Self) -> Self {
        Self(self.0 + rhs.0)
    }

    fn square(&self) -> Self {
        Self(self.0 * self.0)
    }
}

impl<'a> ToU1024 for Fp<'a> {
    fn to_u1024_val(&self) -> U1024 {
        self.0.to_u1024()
    }
}
