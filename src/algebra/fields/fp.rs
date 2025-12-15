use std::ops::{Add, Deref, Mul, Neg, Sub};

use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{FieldElement, U1024};

use crate::traits::{Field, ToU1024};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fp<'a>(pub FieldElement<'a>);

impl<'a> Fp<'a> {
    /// Constructs an `Fp` by creating an inner `FieldElement` from the given value and Montgomery parameters.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use curvelib::traits::ToU1024;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = get_tiny_params();
    /// let value: U1024 = U1024::from_u64(42);
    /// let fp = Fp::new(value, params);
    /// assert_eq!(fp.to_u1024(), U1024::from_u64(3));
    /// ```
    pub fn new(value: U1024, params: &'a MontgomeryParams) -> Self {
        Self(FieldElement::new(value, params))
    }
}

impl<'a> From<FieldElement<'a>> for Fp<'a> {
    /// Wraps a `FieldElement` value in the `Fp` new type.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use mathlib::{FieldElement, U1024, BigInt};
    ///
    /// let params = get_tiny_params();
    /// let fe: FieldElement<'_> = FieldElement::new(U1024::from_u64(7), params);
    /// let fp: Fp<'_> = fe.into();
    /// assert_eq!(fp.to_u1024(), U1024::from_u64(7));
    /// ```
    fn from(f: FieldElement<'a>) -> Self {
        Self(f)
    }
}

// Deref to a field element for convenience (accessing .params, .value)
impl<'a> Deref for Fp<'a> {
    type Target = FieldElement<'a>;
    /// Accesses the underlying `FieldElement` by reference.
    ///
    /// Returns a shared reference to the wrapped `FieldElement`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use std::ops::Deref;
    /// use curvelib::algebra::fields::Fp;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = get_tiny_params();
    /// let fp = Fp::new(U1024::from_u64(1), params);
    /// let inner = fp.deref();
    /// assert_eq!(inner.to_u1024(), U1024::from_u64(1));
    /// ```
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<'a> Add for Fp<'a> {
    type Output = Self;
    /// Adds two `Fp` elements.
    ///
    /// Returns an `Fp` containing the sum of `self` and `rhs`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = get_tiny_params();
    /// let a = Fp::new(U1024::from_u64(5), params);
    /// let b = Fp::new(U1024::from_u64(9), params);
    /// let c = a + b;
    /// // tiny field modulus is 13, so 5 + 9 = 14 ≡ 1 (mod 13)
    /// assert_eq!(c.to_u1024(), U1024::from_u64(1));
    /// ```
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl<'a> Sub for Fp<'a> {
    type Output = Self;
    /// Subtracts one `Fp` value from another and returns the difference.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = get_tiny_params();
    /// let a = Fp::new(U1024::from_u64(2), params);
    /// let b = Fp::new(U1024::from_u64(5), params);
    /// let c = a - b;
    /// // 2 - 5 ≡ -3 ≡ 10 (mod 13)
    /// assert_eq!(c.to_u1024(), U1024::from_u64(10));
    /// ```
    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0)
    }
}

impl<'a> Mul for Fp<'a> {
    type Output = Self;
    /// Multiplies two field elements and returns their product wrapped in `Fp`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = get_tiny_params();
    /// let a = Fp::new(U1024::from_u64(2), params);
    /// let b = Fp::new(U1024::from_u64(3), params);
    /// let c = a * b;
    /// assert_eq!(c.to_u1024(), U1024::from_u64(6));
    /// ```
    fn mul(self, rhs: Self) -> Self {
        Self(self.0 * rhs.0)
    }
}

impl<'a> Neg for Fp<'a> {
    type Output = Self;
    /// Computes the additive inverse of the field element.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use curvelib::traits::Field;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = get_tiny_params();
    /// let a = Fp::new(U1024::from_u64(3), params);
    /// let neg_a = -a;
    /// assert_eq!(a + neg_a, Fp::zero(params));
    /// ```
    fn neg(self) -> Self {
        let zero = FieldElement::zero(self.0.params);
        Self(zero - self.0)
    }
}

impl<'a> Field<'a> for Fp<'a> {
    /// Create the additive identity (zero) in the field defined by `params`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use curvelib::traits::Field;
    ///
    /// let params = get_tiny_params();
    /// let z = Fp::zero(params);
    /// assert!(z.is_zero());
    /// ```
    fn zero(params: &'a MontgomeryParams) -> Self {
        Self(FieldElement::zero(params))
    }

    /// Checks whether the element is the additive zero for its Montgomery parameters.
    ///
    /// # Returns
    ///
    /// `true` if the element equals the field zero for its associated `MontgomeryParams`, `false` otherwise.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use curvelib::traits::Field;
    ///
    /// let params = get_tiny_params();
    /// let z = Fp::zero(params);
    /// assert!(z.is_zero());
    ///
    /// let one = Fp::one(params);
    /// assert!(!one.is_zero());
    /// ```
    fn is_zero(&self) -> bool {
        self.0 == FieldElement::zero(self.0.params)
    }

    /// Constructs the multiplicative identity (one) for the given Montgomery parameters.
    ///
    /// # Returns
    ///
    /// The field element representing 1 for `params`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use curvelib::traits::{Field, ToU1024};
    ///
    /// let params = get_tiny_params();
    /// let one = Fp::one(params);
    /// assert_eq!(one.mul(&one).to_u1024(), one.to_u1024());
    /// ```
    fn one(params: &'a MontgomeryParams) -> Self {
        Self(FieldElement::one(params))
    }

    /// Compute the multiplicative inverse of this field element.
    ///
    /// Returns `Some(Self)` containing the inverse when the element is not zero, or `None` when the element is zero.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use curvelib::traits::{Field, ToU1024};
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = get_tiny_params();
    /// let a = Fp::new(U1024::from_u64(5), params);
    /// let inv = a.inv().expect("nonzero element should have an inverse");
    /// assert_eq!(a.mul(&inv).to_u1024(), Fp::one(params).to_u1024());
    ///
    /// let zero = Fp::zero(params);
    /// assert!(zero.inv().is_none());
    /// ```
    fn inv(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(Self(self.0.inv()))
        }
    }

    /// Computes the element added to itself.
    ///
    /// # Returns
    ///
    /// `Self` equal to this element plus itself.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use curvelib::traits::Field;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = get_tiny_params();
    /// let fp = Fp::new(U1024::from_u64(7), params);
    /// let doubled = fp.double();
    /// assert_eq!(doubled, fp + fp);
    /// ```
    fn double(&self) -> Self {
        Self(self.0 + self.0)
    }

    /// Multiplies two field elements and returns their product.
    ///
    /// # Returns
    ///
    /// The product of `self` and `rhs`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use std::ops::Mul;
    /// use curvelib::algebra::fields::Fp;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = get_tiny_params();
    /// let a = Fp::new(U1024::from_u64(3), params);
    /// let b = Fp::new(U1024::from_u64(4), params);
    /// let c = a.mul(b);
    /// assert_eq!(c.to_u1024(), U1024::from_u64(12));
    /// ```
    fn mul(&self, rhs: &Self) -> Self {
        Self(self.0 * rhs.0)
    }

    /// Adds two field elements and returns their sum.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use std::ops::Add;
    /// use curvelib::algebra::fields::Fp;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = get_tiny_params();
    /// let a = Fp::new(U1024::from_u64(8), params);
    /// let b = Fp::new(U1024::from_u64(6), params);
    /// let sum = a.add(b);
    /// assert_eq!(sum, a + b);
    /// ```
    fn add(&self, rhs: &Self) -> Self {
        Self(self.0 + rhs.0)
    }

    /// Computes the square of the field element.
    ///
    /// # Returns
    ///
    /// `Self` equal to this element multiplied by itself.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use curvelib::traits::Field;
    ///
    /// let params = get_tiny_params();
    /// let a = Fp::one(params);
    /// let s = a.square();
    /// assert_eq!(s, a.mul(&a));
    /// ```
    fn square(&self) -> Self {
        Self(self.0 * self.0)
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
    fn to_u1024_val(&self) -> U1024 {
        self.0.to_u1024()
    }
}
