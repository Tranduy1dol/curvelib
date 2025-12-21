//! Field trait for finite field arithmetic.

use mathlib::FieldConfig;

/// Trait for finite field elements.
///
/// This trait defines the basic operations needed for field arithmetic.
/// It is now generic over a `FieldConfig` type parameter instead of using lifetimes.
pub trait Field: Sized + Clone + Copy + PartialEq + Eq {
    /// Returns the additive identity (zero)
    fn zero() -> Self;

    /// Returns true if this element is zero
    fn is_zero(&self) -> bool;

    /// Returns the multiplicative identity (one)
    fn one() -> Self;

    /// Computes the multiplicative inverse if it exists
    fn inv(&self) -> Option<Self>;

    /// Computes 2 * self
    fn double(&self) -> Self;

    /// Computes self * rhs
    fn mul(&self, rhs: &Self) -> Self;

    /// Computes self + rhs
    fn add(&self, rhs: &Self) -> Self;

    /// Computes self²
    fn square(&self) -> Self;
}

/// Implement Field for mathlib's `FieldElement<C>`
impl<C: FieldConfig> Field for mathlib::FieldElement<C> {
    fn zero() -> Self {
        mathlib::FieldElement::<C>::zero()
    }

    fn is_zero(&self) -> bool {
        mathlib::FieldElement::is_zero(self)
    }

    fn one() -> Self {
        mathlib::FieldElement::<C>::one()
    }

    fn inv(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(mathlib::FieldElement::inv(self))
        }
    }

    fn double(&self) -> Self {
        mathlib::FieldElement::double(self)
    }

    fn mul(&self, rhs: &Self) -> Self {
        *self * *rhs
    }

    fn add(&self, rhs: &Self) -> Self {
        *self + *rhs
    }

    fn square(&self) -> Self {
        mathlib::FieldElement::square(self)
    }
}
