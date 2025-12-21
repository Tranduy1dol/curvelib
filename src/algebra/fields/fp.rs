//! Field element types for elliptic curve cryptography.
//!
//! This module provides type aliases and implementations for field elements
//! used in curve arithmetic. Field elements are generic over a `FieldConfig`
//! type rather than using lifetime parameters.

use mathlib::{FieldConfig, FieldElement, U1024};

use crate::traits::ToU1024;

/// Type alias for prime field elements.
///
/// `Fp<C>` is an alias for `FieldElement<C>` from mathlib, providing convenient
/// shorthand for working with elements in the prime field defined by `C: FieldConfig`.
///
/// # Example
///
/// ```rust,ignore
/// use curvelib::algebra::fields::Fp;
/// use curvelib::instances::bls6_6::Bls6_6BaseField;
///
/// let a: Fp<Bls6_6BaseField> = Fp::new(U1024::from_u64(5));
/// let b: Fp<Bls6_6BaseField> = Fp::one();
/// let c = a + b;
/// ```
pub type Fp<C> = FieldElement<C>;

/// Implement ToU1024 for `FieldElement<C>`
impl<C: FieldConfig> ToU1024 for FieldElement<C> {
    /// Convert this field element into its canonical `U1024` integer representation.
    fn to_u1024(&self) -> U1024 {
        FieldElement::to_u1024(self)
    }
}
