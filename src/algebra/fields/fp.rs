use mathlib::{FieldElement, MontgomeryParams, U1024};

use crate::traits::{Field, ToU1024};

/// Type alias for field elements in the base prime field.
///
/// `Fp` is an alias for `FieldElement` from mathlib, providing a convenient shorthand
/// for working with elements in the prime field defined by Montgomery parameters.
pub type Fp<'a> = FieldElement<'a>;

// Implement curvelib's Field trait for FieldElement
impl<'a> Field<'a> for FieldElement<'a> {
    fn zero(params: &'a MontgomeryParams) -> Self {
        FieldElement::zero(params)
    }

    fn is_zero(&self) -> bool {
        *self == FieldElement::zero(self.params)
    }

    fn one(params: &'a MontgomeryParams) -> Self {
        FieldElement::one(params)
    }

    fn inv(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(FieldElement::inv(self))
        }
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

impl<'a> ToU1024 for Fp<'a> {
    /// Convert this field element into its canonical `U1024` integer representation.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = get_tiny_params();
    /// let fp = Fp::new(U1024::from_u64(11), params);
    /// let int_val = fp.to_u1024();
    /// assert_eq!(int_val, U1024::from_u64(11));
    /// ```
    fn to_u1024(&self) -> U1024 {
        self.to_u1024()
    }
}
