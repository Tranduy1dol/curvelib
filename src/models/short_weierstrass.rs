use mathlib::field::{element::FieldElement, montgomery::MontgomeryParams};
use mathlib::{BigInt, U1024};

use crate::traits::{Curve, ProjectivePoint};

#[derive(Clone, Debug)]
pub struct WeierstrassCurve<'a> {
    pub a: FieldElement<'a>,
    pub b: FieldElement<'a>,
    pub params: &'a MontgomeryParams,
    pub scalar_params: &'a MontgomeryParams,
    pub generator_x: U1024,
    pub generator_y: U1024,
}

impl<'a> WeierstrassCurve<'a> {
    /// Creates a WeierstrassCurve with the given curve coefficients and Montgomery parameters.
    pub fn new(
        a: FieldElement<'a>,
        b: FieldElement<'a>,
        params: &'a MontgomeryParams,
        scalar_params: &'a MontgomeryParams,
        generator_x: U1024,
        generator_y: U1024,
    ) -> Self {
        Self {
            a,
            b,
            params,
            scalar_params,
            generator_x,
            generator_y,
        }
    }
}

impl<'a> Curve<'a> for WeierstrassCurve<'a> {
    type Point = SWPoint<'a>;

    /// Creates the identity (point at infinity) for this curve in projective coordinates.
    ///
    /// The returned `SWPoint` has its projective `z` coordinate set to zero, representing the point at infinity.
    ///
    /// # Examples
    ///
    /// ```
    /// use mathlib::{U1024, BigInt};
    /// use mathlib::field::{element::FieldElement, montgomery::MontgomeryParams};
    /// use curvelib::models::short_weierstrass::WeierstrassCurve;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    ///
    /// let p = U1024::from_u64(17);
    /// let params = MontgomeryParams::new(p, U1024::zero());
    /// let a = FieldElement::new(U1024::from_u64(1), &params);
    /// let b = FieldElement::new(U1024::from_u64(1), &params);
    /// let curve = WeierstrassCurve::new(a, b, &params, &params, U1024::from_u64(1), U1024::from_u64(1));
    ///
    /// let id = curve.identity();
    /// assert!(id.is_identity());
    /// ```
    fn identity(&self) -> Self::Point {
        let curve = self.clone();
        let one = FieldElement::one(curve.params);
        let zero = FieldElement::zero(curve.params);
        SWPoint {
            x: one,
            y: one,
            z: zero,
            curve,
        }
    }

    /// Determines whether the given affine point lies on this Weierstrass curve.
    ///
    /// - `x`: Affine x-coordinate of the point.
    /// - `y`: Affine y-coordinate of the point.
    ///
    /// # Returns
    /// `true` if `y^2 = x^3 + a*x + b` holds for this curve's parameters, `false` otherwise.
    fn is_on_curve(&self, x: &FieldElement, y: &FieldElement) -> bool {
        let y2 = *y * *y;
        let x2 = *x * *x;
        let x3 = x2 * *x;
        let ax = self.a * *x;
        let rhs = x3 + ax + self.b;
        y2 == rhs
    }

    fn scalar_params(&self) -> &'a MontgomeryParams {
        self.scalar_params
    }

    fn generator(&self) -> Self::Point {
        let x = FieldElement::new(self.generator_x, self.params);
        let y = FieldElement::new(self.generator_y, self.params);
        let z = FieldElement::one(self.params);
        SWPoint {
            x,
            y,
            z,
            curve: self.clone(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct SWPoint<'a> {
    pub x: FieldElement<'a>,
    pub y: FieldElement<'a>,
    pub z: FieldElement<'a>,
    pub curve: WeierstrassCurve<'a>,
}

impl<'a> PartialEq for SWPoint<'a> {
    fn eq(&self, other: &Self) -> bool {
        if self.is_identity() {
            return other.is_identity();
        }
        if other.is_identity() {
            return self.is_identity();
        }

        let (x1, y1) = self.to_affine();
        let (x2, y2) = other.to_affine();
        x1 == x2 && y1 == y2
    }
}

impl<'a> Eq for SWPoint<'a> {}

impl<'a> SWPoint<'a> {
    /// Constructs a projective curve point from affine coordinates, treating (0, 0) as the point at infinity.
    ///
    /// If `x` and `y` are both zero, the curve identity is returned; otherwise the point is returned with `z = 1`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mathlib::{U1024, BigInt};
    /// use mathlib::field::{element::FieldElement, montgomery::MontgomeryParams};
    /// use curvelib::models::short_weierstrass::{WeierstrassCurve, SWPoint};
    /// use curvelib::traits::{Curve, ProjectivePoint};
    ///
    /// let p = U1024::from_u64(17);
    /// let params = MontgomeryParams::new(p, U1024::zero());
    /// let a = FieldElement::new(U1024::from_u64(1), &params);
    /// let b = FieldElement::new(U1024::from_u64(1), &params);
    /// let curve = WeierstrassCurve::new(a, b, &params, &params, U1024::from_u64(1), U1024::from_u64(1));
    ///
    /// let x = FieldElement::one(&params);
    /// let y = FieldElement::one(&params);
    /// let p = SWPoint::new_affine(x, y, &curve);
    /// // Note: (1,1) might not be on the curve y^2 = x^3 + x + 1 mod 17, but new_affine constructs it anyway.
    /// assert!(!p.is_identity());
    ///
    /// let z0 = FieldElement::zero(&params);
    /// let id = SWPoint::new_affine(z0.clone(), z0, &curve);
    /// assert!(id.is_identity());
    /// ```
    pub fn new_affine(
        x: FieldElement<'a>,
        y: FieldElement<'a>,
        curve: &'a WeierstrassCurve<'a>,
    ) -> Self {
        if x.value == U1024::zero() && y.value == U1024::zero() {
            return curve.identity();
        }
        let z = FieldElement::one(curve.params);
        Self {
            x,
            y,
            z,
            curve: curve.clone(),
        }
    }
}

impl<'a> ProjectivePoint for SWPoint<'a> {
    /// Checks whether this point is the identity (point at infinity).
    ///
    /// The point is considered the identity when its projective `z` coordinate is zero.
    ///
    /// # Returns
    ///
    /// `true` if the point's `z` coordinate equals zero, `false` otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// use mathlib::{U1024, BigInt};
    /// use mathlib::field::{element::FieldElement, montgomery::MontgomeryParams};
    /// use curvelib::models::short_weierstrass::WeierstrassCurve;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    ///
    /// let p = U1024::from_u64(17);
    /// let params = MontgomeryParams::new(p, U1024::zero());
    /// let a = FieldElement::new(U1024::from_u64(1), &params);
    /// let b = FieldElement::new(U1024::from_u64(1), &params);
    /// let curve = WeierstrassCurve::new(a, b, &params, &params, U1024::from_u64(1), U1024::from_u64(1));
    ///
    /// let id = curve.identity();
    /// assert!(id.is_identity());
    /// ```
    fn is_identity(&self) -> bool {
        self.z.value == U1024::zero()
    }

    /// Adds two points on the Weierstrass curve using projective coordinates.
    ///
    /// Performs elliptic-curve point addition in projective form and returns the resulting point.
    /// If either operand is the identity, the other operand is returned. If the points are
    /// inverses of each other, the curve identity is returned.
    ///
    /// # Examples
    ///
    /// ```
    /// use mathlib::{U1024, BigInt};
    /// use mathlib::field::{element::FieldElement, montgomery::MontgomeryParams};
    /// use curvelib::models::short_weierstrass::WeierstrassCurve;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    ///
    /// let p = U1024::from_u64(17);
    /// let params = MontgomeryParams::new(p, U1024::zero());
    /// let a = FieldElement::new(U1024::from_u64(1), &params);
    /// let b = FieldElement::new(U1024::from_u64(1), &params);
    /// let curve = WeierstrassCurve::new(a, b, &params, &params, U1024::from_u64(1), U1024::from_u64(1));
    ///
    /// let p = curve.identity();
    /// let q = curve.identity();
    /// let r = p.add(&q);
    /// assert!(r.is_identity());
    /// ```
    ///
    /// # Returns
    ///
    /// The projective `SWPoint` representing the sum of `self` and `rhs`.
    fn add(&self, rhs: &Self) -> Self {
        if self.is_identity() {
            return rhs.clone();
        }
        if rhs.is_identity() {
            return self.clone();
        }

        let z1z1 = self.z * self.z;
        let z2z2 = rhs.z * rhs.z;

        let u1 = self.x * z2z2;
        let u2 = rhs.x * z1z1;

        let s1 = self.y * (rhs.z * z2z2);
        let s2 = rhs.y * (self.z * z1z1);

        if u1 == u2 {
            return if s1 == s2 {
                self.double()
            } else {
                self.curve.identity()
            };
        }

        let h = u2 - u1;
        let r = s2 - s1;
        let hh = h * h;
        let hhh = hh * h;

        let v = u1 * hh;
        let two = FieldElement::new(U1024::from_u64(2), self.curve.params);

        let x3 = (r * r) - hhh - (two * v);

        let y3 = (r * (v - x3)) - (s1 * hhh);

        let z3 = self.z * rhs.z * h;

        Self {
            x: x3,
            y: y3,
            z: z3,
            curve: self.curve.clone(),
        }
    }

    /// Doubles this point on its associated Weierstrass curve using projective coordinates.
    ///
    /// Returns the point 2P; if this point is the identity (point at infinity), returns a clone of it.
    ///
    /// # Examples
    ///
    /// ```
    /// use mathlib::{U1024, BigInt};
    /// use mathlib::field::{element::FieldElement, montgomery::MontgomeryParams};
    /// use curvelib::models::short_weierstrass::WeierstrassCurve;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    ///
    /// let p = U1024::from_u64(17);
    /// let params = MontgomeryParams::new(p, U1024::zero());
    /// let a = FieldElement::new(U1024::from_u64(1), &params);
    /// let b = FieldElement::new(U1024::from_u64(1), &params);
    /// let curve = WeierstrassCurve::new(a, b, &params, &params, U1024::from_u64(1), U1024::from_u64(1));
    ///
    /// let p = curve.identity();
    /// let r = p.double();
    /// assert!(r.is_identity());
    /// ```
    fn double(&self) -> Self {
        if self.is_identity() {
            return self.clone();
        }

        let xx = self.x * self.x;
        let yy = self.y * self.y;
        let yyyy = yy * yy;
        let zz = self.z * self.z;

        let two = FieldElement::new(U1024::from_u64(2), self.curve.params);
        let s = two * ((self.x * yy) * two);

        let three = FieldElement::new(U1024::from_u64(3), self.curve.params);
        let zzzz = zz * zz;
        let m = (three * xx) + (self.curve.a * zzzz);

        let x_new = (m * m) - (s * two);

        let z_new = (self.y * self.z) * two;

        let eight = FieldElement::new(U1024::from_u64(8), self.curve.params);
        let t = eight * yyyy;
        let y_new = (m * (s - x_new)) - t;

        Self {
            x: x_new,
            y: y_new,
            z: z_new,
            curve: self.curve.clone(),
        }
    }

    /// Convert this projective point to affine coordinates.
    ///
    /// Returns a pair `(x, y)` representing the affine coordinates of the point. If the point is the
    /// identity (point at infinity), it returns `(0, 0)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mathlib::{U1024, BigInt};
    /// use mathlib::field::{element::FieldElement, montgomery::MontgomeryParams};
    /// use curvelib::models::short_weierstrass::WeierstrassCurve;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    ///
    /// let p = U1024::from_u64(17);
    /// let params = MontgomeryParams::new(p, U1024::zero());
    /// let a = FieldElement::new(U1024::from_u64(1), &params);
    /// let b = FieldElement::new(U1024::from_u64(1), &params);
    /// let curve = WeierstrassCurve::new(a, b, &params, &params, U1024::from_u64(1), U1024::from_u64(1));
    ///
    /// let p = curve.identity();
    /// let (x, y) = p.to_affine();
    /// assert_eq!(x, FieldElement::zero(&params));
    /// assert_eq!(y, FieldElement::zero(&params));
    /// ```
    fn to_affine(&self) -> (FieldElement<'a>, FieldElement<'a>) {
        if self.is_identity() {
            let zero = FieldElement::zero(self.curve.params);
            return (zero, zero);
        }

        let z_inv = self.z.inv();
        let z2_inv = z_inv * z_inv;
        let z3_inv = z2_inv * z_inv;

        let x_aff = self.x * z2_inv;
        let y_aff = self.y * z3_inv;

        (x_aff, y_aff)
    }

    /// Multiplies this point by a 1024-bit scalar using the double-and-add algorithm.
    ///
    /// The operation returns k * P where k is `scalar` and P is `self`.
    ///
    /// # Examples
    ///
    /// ```
    /// use mathlib::{U1024, BigInt};
    /// use mathlib::field::{element::FieldElement, montgomery::MontgomeryParams};
    /// use curvelib::models::short_weierstrass::{WeierstrassCurve, SWPoint};
    /// use curvelib::traits::{Curve, ProjectivePoint};
    ///
    /// let p = U1024::from_u64(17);
    /// let params = MontgomeryParams::new(p, U1024::zero());
    /// let a = FieldElement::new(U1024::from_u64(1), &params);
    /// let b = FieldElement::new(U1024::from_u64(1), &params);
    /// let curve = WeierstrassCurve::new(a, b, &params, &params, U1024::from_u64(1), U1024::from_u64(1));
    ///
    /// let p = curve.identity();
    /// let k = U1024::from_u64(3);
    /// let r = p.mul(&k);
    /// assert!(r.is_identity());
    /// ```
    fn mul(&self, scalar: &U1024) -> Self {
        let mut res = self.curve.identity();

        let num_bits = scalar.bits();
        for i in (0..num_bits).rev() {
            res = res.double();
            if scalar.bit(i) {
                res = res.add(self);
            }
        }
        res
    }
}
