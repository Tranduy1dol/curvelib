use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{BigInt, U1024};

use crate::algebra::fields::fp::Fp;
use crate::def_weierstrass_curve;
use crate::traits::{Curve, Field, ProjectivePoint};

def_weierstrass_curve!(WeierstrassCurve, Fp<'a>);

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
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use curvelib::models::short_weierstrass::WeierstrassCurve;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use curvelib::algebra::fields::fp::Fp;
    /// use curvelib::traits::Field;
    ///
    /// let p = U1024::from_u64(17);
    /// let params = MontgomeryParams::new(p, U1024::zero());
    /// let a = Fp::new(U1024::from_u64(1), &params);
    /// let b = Fp::new(U1024::from_u64(1), &params);
    /// let curve = WeierstrassCurve::new(a, b, &params, &params, Fp::new(U1024::from_u64(1), &params), Fp::new(U1024::from_u64(1), &params));
    ///
    /// let id = curve.identity();
    /// assert!(id.is_identity());
    /// ```
    fn identity(&self) -> Self::Point {
        let curve = self.clone();
        let one = Fp::one(self.params);
        let zero = Fp::zero(self.params);
        SWPoint {
            x: one,
            y: one,
            z: zero,
            curve,
        }
    }

    /// Returns whether the affine point (x, y) satisfies this curve's short Weierstrass equation.
    ///
    /// # Returns
    ///
    /// `true` if `y^2 = x^3 + a*x + b` for this curve's parameters, `false` otherwise.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::models::short_weierstrass::WeierstrassCurve;
    /// use curvelib::traits::{Curve, ProjectivePoint, Field};
    /// use curvelib::algebra::fields::fp::Fp;
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let p = U1024::from_u64(43);
    /// let params = MontgomeryParams::new(p, U1024::zero());
    /// let a = Fp::new(U1024::from_u64(0), &params);
    /// let b = Fp::new(U1024::from_u64(6), &params);
    /// let curve = WeierstrassCurve::new(a, b, &params, &params, Fp::new(U1024::from_u64(13), &params), Fp::new(U1024::from_u64(15), &params));
    ///
    /// let g = curve.generator();
    /// let (x, y) = g.to_affine();
    /// assert!(curve.is_on_curve(&x, &y));
    /// ```
    fn is_on_curve(
        &self,
        x: &<Self::Point as ProjectivePoint<'a>>::Field,
        y: &<Self::Point as ProjectivePoint<'a>>::Field,
    ) -> bool {
        let y2 = *y * *y;
        let x2 = *x * *x;
        let x3 = x2 * *x;
        let ax = self.a * *x;
        let rhs = x3 + ax + self.b;
        y2 == rhs
    }

    /// Accesses the curve's Montgomery scalar parameters.
    ///
    /// Returns a reference to the `MontgomeryParams` used for scalar arithmetic on this curve.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::models::short_weierstrass::WeierstrassCurve;
    /// use curvelib::algebra::fields::fp::Fp;
    /// use curvelib::traits::{Curve, Field};
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use mathlib::{BigInt, U1024};
    ///
    /// let p = U1024::from_u64(17);
    /// let params = MontgomeryParams::new(p, U1024::zero());
    /// let curve = WeierstrassCurve::new(
    ///     Fp::new(U1024::from_u64(1), &params),
    ///     Fp::new(U1024::from_u64(1), &params),
    ///     &params,
    ///     &params,
    ///     Fp::one(&params),
    ///     Fp::one(&params),
    /// );
    /// assert_eq!(curve.scalar_params().modulus, U1024::from_u64(17));
    /// ```
    fn scalar_params(&self) -> &'a MontgomeryParams {
        self.scalar_params
    }

    /// Constructs the curve's generator point in projective coordinates (with z = 1).
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::models::short_weierstrass::WeierstrassCurve;
    /// use curvelib::algebra::fields::fp::Fp;
    /// use curvelib::traits::{Curve, ProjectivePoint, Field};
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let p = U1024::from_u64(17);
    /// let params = MontgomeryParams::new(p, U1024::zero());
    /// let a = Fp::new(U1024::from_u64(1), &params);
    /// let b = Fp::new(U1024::from_u64(1), &params);
    /// let curve = WeierstrassCurve::new(a, b, &params, &params, Fp::one(&params), Fp::one(&params));
    ///
    /// let g = curve.generator();
    /// assert!(!g.is_identity());
    /// ```
    fn generator(&self) -> Self::Point {
        let x = self.generator_x;
        let y = self.generator_y;
        let z = Fp::one(self.params);
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
    pub x: Fp<'a>,
    pub y: Fp<'a>,
    pub z: Fp<'a>,
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
    /// Constructs a projective point from affine coordinates, treating (0, 0) as the point at infinity.
    ///
    /// If both `x` and `y` are zero, the curve identity is returned; otherwise returns the point with `z = 1`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::models::short_weierstrass::{WeierstrassCurve, SWPoint};
    /// use curvelib::algebra::fields::fp::Fp;
    /// use curvelib::traits::{Curve, ProjectivePoint, Field};
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = MontgomeryParams::new(U1024::from_u64(17), U1024::zero());
    /// let a = Fp::new(U1024::from_u64(1), &params);
    /// let b = Fp::new(U1024::from_u64(1), &params);
    /// let curve = WeierstrassCurve::new(a, b, &params, &params, Fp::one(&params), Fp::one(&params));
    ///
    /// let p = SWPoint::new_affine(Fp::one(&params), Fp::one(&params), &curve);
    /// assert!(!p.is_identity());
    ///
    /// let id = SWPoint::new_affine(Fp::zero(&params), Fp::zero(&params), &curve);
    /// assert!(id.is_identity());
    /// ```
    pub fn new_affine(x: Fp<'a>, y: Fp<'a>, curve: &'a WeierstrassCurve<'a>) -> Self {
        if x.value == U1024::zero() && y.value == U1024::zero() {
            return curve.identity();
        }
        let z = Fp::one(curve.params);
        Self {
            x,
            y,
            z,
            curve: curve.clone(),
        }
    }
}

impl<'a> ProjectivePoint<'a> for SWPoint<'a> {
    type Field = Fp<'a>;

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
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use curvelib::models::short_weierstrass::WeierstrassCurve;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use curvelib::algebra::fields::fp::Fp;
    /// use curvelib::traits::Field;
    ///
    /// let p = U1024::from_u64(17);
    /// let params = MontgomeryParams::new(p, U1024::zero());
    /// let a = Fp::new(U1024::from_u64(1), &params);
    /// let b = Fp::new(U1024::from_u64(1), &params);
    /// let curve = WeierstrassCurve::new(a, b, &params, &params, Fp::new(U1024::from_u64(1), &params), Fp::new(U1024::from_u64(1), &params));
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
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use curvelib::models::short_weierstrass::WeierstrassCurve;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use curvelib::algebra::fields::fp::Fp;
    /// use curvelib::traits::Field;
    ///
    /// let p = U1024::from_u64(17);
    /// let params = MontgomeryParams::new(p, U1024::zero());
    /// let a = Fp::new(U1024::from_u64(1), &params);
    /// let b = Fp::new(U1024::from_u64(1), &params);
    /// let curve = WeierstrassCurve::new(a, b, &params, &params, Fp::new(U1024::from_u64(1), &params), Fp::new(U1024::from_u64(1), &params));
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
        let two_val = U1024::from_u64(2);
        let two = Fp::new(two_val, self.curve.params);

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

    /// Compute the doubling of this point on its short Weierstrass curve using projective coordinates.
    ///
    /// # Returns
    ///
    /// `Self` representing the point 2P; if this point is the identity (point at infinity), returns a clone of it.
    ///
    /// # Examples
    ///
    /// ```
    /// use mathlib::{U1024, BigInt};
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use curvelib::models::short_weierstrass::WeierstrassCurve;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use curvelib::algebra::fields::fp::Fp;
    ///
    /// let p = U1024::from_u64(17);
    /// let params = MontgomeryParams::new(p, U1024::zero());
    /// let a = Fp::new(U1024::from_u64(1), &params);
    /// let b = Fp::new(U1024::from_u64(1), &params);
    /// let curve = WeierstrassCurve::new(a, b, &params, &params, Fp::new(U1024::from_u64(1), &params), Fp::new(U1024::from_u64(1), &params));
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

        let two_val = U1024::from_u64(2);
        let two = Fp::new(two_val, self.curve.params);
        let s = two * ((self.x * yy) * two);

        let three_val = U1024::from_u64(3);
        let three = Fp::new(three_val, self.curve.params);
        let zzzz = zz * zz;
        let m = (three * xx) + (self.curve.a * zzzz);

        let x_new = (m * m) - (s * two);

        let z_new = (self.y * self.z) * two;

        let eight_val = U1024::from_u64(8);
        let eight = Fp::new(eight_val, self.curve.params);
        let t = eight * yyyy;
        let y_new = (m * (s - x_new)) - t;

        Self {
            x: x_new,
            y: y_new,
            z: z_new,
            curve: self.curve.clone(),
        }
    }

    /// Converts this projective point to affine coordinates.
    ///
    /// Returns `(x, y)` as affine coordinates; returns `(0, 0)` when the point is the identity (point at infinity).
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::models::short_weierstrass::WeierstrassCurve;
    /// use curvelib::algebra::fields::fp::Fp;
    /// use curvelib::traits::{Curve, ProjectivePoint, Field};
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let p = U1024::from_u64(17);
    /// let params = MontgomeryParams::new(p, U1024::zero());
    /// let a = Fp::new(U1024::from_u64(1), &params);
    /// let b = Fp::new(U1024::from_u64(1), &params);
    /// let curve = WeierstrassCurve::new(a, b, &params, &params, Fp::one(&params), Fp::one(&params));
    ///
    /// let id = curve.identity();
    /// let (x, y) = id.to_affine();
    /// assert_eq!(x, Fp::zero(&params));
    /// assert_eq!(y, Fp::zero(&params));
    /// ```
    fn to_affine(&self) -> (Fp<'a>, Fp<'a>) {
        if self.is_identity() {
            let zero = Fp::zero(self.curve.params);
            return (zero, zero);
        }

        let z_inv = self.z.inv().unwrap();
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
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use curvelib::models::short_weierstrass::{WeierstrassCurve, SWPoint};
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use curvelib::algebra::fields::fp::Fp;
    /// use curvelib::traits::Field;
    ///
    /// let p = U1024::from_u64(17);
    /// let params = MontgomeryParams::new(p, U1024::zero());
    /// let a = Fp::new(U1024::from_u64(1), &params);
    /// let b = Fp::new(U1024::from_u64(1), &params);
    /// let curve = WeierstrassCurve::new(a, b, &params, &params, Fp::new(U1024::from_u64(1), &params), Fp::new(U1024::from_u64(1), &params));
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
