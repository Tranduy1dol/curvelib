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
    /// Multiply two Fp6 field elements.
    ///
    /// The product is computed in the cubic extension Fp6 = Fp2[v]/(v^3 - ξ) using the
    /// relation v^3 = ξ, v^4 = ξ v, v^5 = ξ v^2. This implementation assumes the
    /// cubic non-residue `ξ` equals `1`.
    ///
    /// # Returns
    ///
    /// `Self` representing the product of the two Fp6 elements.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp6;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use curvelib::traits::Field;
    ///
    /// let params = get_tiny_params();
    /// let a = Fp6::one(params);
    /// let b = Fp6::one(params);
    /// let product = a * b;
    /// assert_eq!(product, Fp6::one(params));
    /// ```
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
    /// Constructs the additive identity element of Fp6 using the given Montgomery parameters.
    ///
    /// All three underlying Fp2 components (c0, c1, c2) are set to zero.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp6;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use curvelib::traits::Field;
    ///
    /// let params = get_tiny_params();
    /// let z = Fp6::zero(params);
    /// assert!(z.is_zero());
    /// ```
    fn zero(params: &'a MontgomeryParams) -> Self {
        let z = Fp2::zero(params);
        Self {
            c0: z,
            c1: z,
            c2: z,
        }
    }

    /// Checks whether the element is the additive zero.
    ///
    /// # Returns
    /// `true` if all three Fp2 components (`c0`, `c1`, and `c2`) are zero, `false` otherwise.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp6;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use curvelib::traits::Field;
    ///
    /// let params = get_tiny_params();
    /// let z = Fp6::zero(params);
    /// assert!(z.is_zero());
    ///
    /// let one = Fp6::one(params);
    /// assert!(!one.is_zero());
    /// ```
    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero() && self.c2.is_zero()
    }

    /// Constructs the multiplicative identity element of Fp6.
    ///
    /// The returned element has c0 = Fp2::one(params) and c1 = c2 = Fp2::zero(params).
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::{Fp2, Fp6};
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use curvelib::traits::Field;
    ///
    /// let params = get_tiny_params();
    /// let one = Fp6::one(params);
    /// assert_eq!(one.c0, Fp2::one(params));
    /// assert_eq!(one.c1, Fp2::zero(params));
    /// assert_eq!(one.c2, Fp2::zero(params));
    /// ```
    fn one(params: &'a MontgomeryParams) -> Self {
        let z = Fp2::zero(params);
        let o = Fp2::one(params);
        Self {
            c0: o,
            c1: z,
            c2: z,
        }
    }

    /// Attempts to compute the multiplicative inverse of this Fp6 element.
    ///
    /// The inversion operation is not implemented in this type, and this method always returns `None`.
    /// Currently, this function is too complex to implement efficiently, so it is not implemented.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp6;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use curvelib::traits::Field;
    ///
    /// let params = get_tiny_params();
    /// let a = Fp6::one(params);
    /// assert_eq!(a.inv(), None);
    /// ```
    fn inv(&self) -> Option<Self> {
        None
    }

    /// Compute the sum of the element with itself.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp6;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use curvelib::traits::Field;
    ///
    /// let params = get_tiny_params();
    /// let a = Fp6::zero(params);
    /// let doubled = a.double();
    /// assert!(doubled.is_zero());
    /// ```
    fn double(&self) -> Self {
        *self + *self
    }

    /// Multiplies two Fp6 elements and returns their product.
    ///
    /// # Returns
    ///
    /// The product of `self` and `rhs` as an `Fp6`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp6;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use curvelib::traits::Field;
    ///
    /// let params = get_tiny_params();
    /// let a = Fp6::one(params);
    /// let b = Fp6::one(params);
    /// let c = a.mul(&b);
    /// assert_eq!(c, Fp6::one(params));
    /// ```
    fn mul(&self, rhs: &Self) -> Self {
        *self * *rhs
    }

    /// Adds two Fp6 elements and returns their sum.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp6;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use curvelib::traits::Field;
    ///
    /// let params = get_tiny_params();
    /// let a = Fp6::zero(params);
    /// let b = Fp6::one(params);
    /// let c = a.add(&b);
    /// assert_eq!(c, b);
    /// ```
    fn add(&self, rhs: &Self) -> Self {
        *self + *rhs
    }

    /// Squares the field element.
    ///
    /// # Returns
    ///
    /// The product of the element with itself.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::Fp6;
    /// use curvelib::instances::tiny_jubjub::get_tiny_params;
    /// use curvelib::traits::Field;
    ///
    /// let params = get_tiny_params();
    /// let x = Fp6::one(params);
    /// let y = x.square();
    /// assert_eq!(y, x * x);
    /// ```
    fn square(&self) -> Self {
        *self * *self
    }
}
