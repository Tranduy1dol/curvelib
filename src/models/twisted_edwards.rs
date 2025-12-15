use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{BigInt, U1024};

use crate::algebra::fields::fp::Fp;
use crate::traits::{Curve, Field, ProjectivePoint};

impl<'a> Curve<'a> for EdwardsCurve<'a> {
    type Point = TePoint<'a>;

    /// Returns the identity point on this Edwards curve in extended/projective coordinates.
    ///
    /// The identity has coordinates x = 0, y = 1, z = 1, t = 0 and is bound to this curve.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::fp::Fp;
    /// use curvelib::models::twisted_edwards::EdwardsCurve;
    /// use curvelib::traits::{Curve, ProjectivePoint, Field};
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = MontgomeryParams::new(U1024::from_u64(17), U1024::zero());
    /// let a = Fp::new(U1024::from_u64(1), &params);
    /// let d = Fp::new(U1024::from_u64(2), &params);
    /// let curve = EdwardsCurve::new(a, d, &params, &params, U1024::from_u64(0), U1024::from_u64(1));
    ///
    /// let id = curve.identity();
    /// let (x, y) = id.to_affine();
    /// assert_eq!(x, Fp::zero(&params));
    /// assert_eq!(y, Fp::one(&params));
    /// ```
    fn identity(&self) -> Self::Point {
        let curve = self.clone();
        let zero = Fp::zero(self.params);
        let one = Fp::one(self.params);
        TePoint {
            x: zero,
            y: one,
            z: one,
            t: zero,
            curve,
        }
    }

    /// Determines whether an affine point (x, y) satisfies the Edwards curve equation.
    ///
    /// Returns `true` if the coordinates satisfy a*x^2 + y^2 = 1 + d*x^2*y^2, `false` otherwise.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use mathlib::{U1024, BigInt};
    /// use mathlib::field::{montgomery::MontgomeryParams};
    /// use curvelib::models::twisted_edwards::EdwardsCurve;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use curvelib::traits::Field;
    /// use curvelib::algebra::fields::fp::Fp;
    ///
    /// let p = U1024::from_u64(17);
    /// let params = MontgomeryParams::new(p, U1024::zero());
    /// let a = Fp::new(U1024::from_u64(1), &params);
    /// let d = Fp::new(U1024::from_u64(2), &params);
    /// let curve = EdwardsCurve::new(a, d, &params, &params, U1024::from_u64(0), U1024::from_u64(1));
    ///
    /// let x = Fp::zero(&params);
    /// let y = Fp::one(&params);
    /// assert!(curve.is_on_curve(&x, &y));
    /// ```
    fn is_on_curve(
        &self,
        x: &<Self::Point as ProjectivePoint<'a>>::Field,
        y: &<Self::Point as ProjectivePoint<'a>>::Field,
    ) -> bool {
        let x2 = *x * *x;
        let y2 = *y * *y;
        let lhs = (self.a * x2) + y2;

        let one = Fp::one(self.params);
        let rhs = one + (self.d * x2 * y2);

        lhs == rhs
    }

    /// Returns a reference to the curve's scalar-field Montgomery parameters.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::fp::Fp;
    /// use curvelib::models::twisted_edwards::EdwardsCurve;
    /// use curvelib::traits::Curve;
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = MontgomeryParams::new(U1024::from_u64(17), U1024::zero());
    /// let curve = EdwardsCurve::new(
    ///     Fp::new(U1024::from_u64(1), &params),
    ///     Fp::new(U1024::from_u64(2), &params),
    ///     &params,
    ///     &params,
    ///     U1024::from_u64(0),
    ///     U1024::from_u64(1),
    /// );
    /// let p1 = curve.scalar_params();
    /// let p2 = curve.scalar_params();
    /// assert!(std::ptr::eq(p1, p2));
    /// ```
    fn scalar_params(&self) -> &'a MontgomeryParams {
        self.scalar_params
    }

    /// Constructs the curve's generator point in extended/projective coordinates.
    ///
    /// The returned point uses the curve's stored generator_x and generator_y as affine
    /// coordinates, sets `z` to one, computes `t = x * y`, and clones the curve reference.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::fp::Fp;
    /// use curvelib::models::twisted_edwards::EdwardsCurve;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = MontgomeryParams::new(U1024::from_u64(17), U1024::zero());
    /// let a = Fp::new(U1024::from_u64(1), &params);
    /// let d = Fp::new(U1024::from_u64(2), &params);
    /// let curve = EdwardsCurve::new(a, d, &params, &params, U1024::from_u64(0), U1024::from_u64(1));
    ///
    /// let g = curve.generator();
    /// let (gx, gy) = g.to_affine();
    /// assert!(curve.is_on_curve(&gx, &gy));
    /// ```
    fn generator(&self) -> Self::Point {
        let x = Fp::new(self.generator_x, self.params);
        let y = Fp::new(self.generator_y, self.params);
        let z = Fp::one(self.params);
        let t = x * y;
        TePoint {
            x,
            y,
            z,
            t,
            curve: self.clone(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct EdwardsCurve<'a> {
    pub a: Fp<'a>,
    pub d: Fp<'a>,
    pub params: &'a MontgomeryParams,
    pub scalar_params: &'a MontgomeryParams,
    pub generator_x: U1024,
    pub generator_y: U1024,
}

impl<'a> EdwardsCurve<'a> {
    /// Creates a new EdwardsCurve with the specified curve parameters and generator.
    ///
    /// - `a` and `d` are the Edwards curve parameters in the base field.
    /// - `params` are the Montgomery parameters for the base field.
    /// - `scalar_params` are the Montgomery parameters for the scalar field.
    /// - `generator_x` and `generator_y` are the affine coordinates of the curve generator as `U1024`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::fp::Fp;
    /// use curvelib::models::twisted_edwards::EdwardsCurve;
    /// use curvelib::traits::Field;
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = MontgomeryParams::new(U1024::from_u64(17), U1024::zero());
    /// let scalar_params = &params;
    ///
    /// let a = Fp::one(&params);
    /// let d = Fp::new(U1024::from_u64(5), &params);
    /// let gx = U1024::from_u64(0);
    /// let gy = U1024::from_u64(1);
    ///
    /// let curve = EdwardsCurve::new(a, d, &params, scalar_params, gx, gy);
    /// assert_eq!(curve.generator_x, gx);
    /// assert_eq!(curve.generator_y, gy);
    /// ```
    pub fn new(
        a: Fp<'a>,
        d: Fp<'a>,
        params: &'a MontgomeryParams,
        scalar_params: &'a MontgomeryParams,
        generator_x: U1024,
        generator_y: U1024,
    ) -> Self {
        Self {
            a,
            d,
            params,
            scalar_params,
            generator_x,
            generator_y,
        }
    }
}

#[derive(Clone, Debug)]
pub struct TePoint<'a> {
    pub x: Fp<'a>,
    pub y: Fp<'a>,
    pub z: Fp<'a>,
    pub t: Fp<'a>,
    pub curve: EdwardsCurve<'a>,
}

impl<'a> PartialEq for TePoint<'a> {
    fn eq(&self, other: &Self) -> bool {
        if self.is_identity() {
            return other.is_identity();
        }
        if other.is_identity() {
            return self.is_identity();
        }

        let x1z2 = self.x * other.z;
        let x2z1 = other.x * self.z;

        let y1z2 = self.y * other.z;
        let y2z1 = other.y * self.z;

        x1z2 == x2z1 && y1z2 == y2z1
    }
}
impl<'a> Eq for TePoint<'a> {}

impl<'a> TePoint<'a> {
    /// Creates a projective/extended Edwards curve point from affine coordinates.
    ///
    /// The resulting point has `z` set to the field multiplicative identity and `t` set to `x * y`.
    /// The provided curve reference is cloned and stored with the point.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::fp::Fp;
    /// use curvelib::models::twisted_edwards::{EdwardsCurve, TePoint};
    /// use curvelib::traits::{Curve, ProjectivePoint, Field};
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = MontgomeryParams::new(U1024::from_u64(17), U1024::zero());
    /// let a = Fp::new(U1024::from_u64(1), &params);
    /// let d = Fp::new(U1024::from_u64(2), &params);
    /// let curve = EdwardsCurve::new(a, d, &params, &params, U1024::from_u64(0), U1024::from_u64(1));
    ///
    /// let x = Fp::zero(&params);
    /// let y = Fp::one(&params);
    /// let p = TePoint::new_affine(x, y, &curve);
    /// let (xa, ya) = p.to_affine();
    /// assert_eq!(xa, Fp::zero(&params));
    /// assert_eq!(ya, Fp::one(&params));
    /// ```
    pub fn new_affine(x: Fp<'a>, y: Fp<'a>, curve: &'a EdwardsCurve<'a>) -> Self {
        let z = Fp::one(curve.params);
        let t = x * y;
        Self {
            x,
            y,
            z,
            t,
            curve: curve.clone(),
        }
    }
}

impl<'a> ProjectivePoint<'a> for TePoint<'a> {
    type Field = Fp<'a>;

    /// Determines whether the point is the identity element on its Edwards curve.
    ///
    /// The identity is encoded in extended/projective coordinates as `x == 0` and `y == z`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::fp::Fp;
    /// use curvelib::models::twisted_edwards::{EdwardsCurve, TePoint};
    /// use curvelib::traits::{Curve, ProjectivePoint, Field};
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = MontgomeryParams::new(U1024::from_u64(17), U1024::zero());
    /// let a = Fp::new(U1024::from_u64(1), &params);
    /// let d = Fp::new(U1024::from_u64(2), &params);
    /// let curve = EdwardsCurve::new(a, d, &params, &params, U1024::from_u64(0), U1024::from_u64(1));
    /// let id = TePoint::new_affine(Fp::zero(&params), Fp::one(&params), &curve);
    /// assert!(id.is_identity());
    /// ```
    fn is_identity(&self) -> bool {
        let zero = Fp::zero(self.curve.params);
        self.x == zero && self.y == self.z
    }

    /// Adds two Edwards-curve points in projective/extended coordinates.
    ///
    /// The result is the curve point representing `self + rhs` using the curve's
    /// parameters stored in the points.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use std::ops::Add;
    /// use curvelib::algebra::fields::fp::Fp;
    /// use curvelib::models::twisted_edwards::EdwardsCurve;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = MontgomeryParams::new(U1024::from_u64(17), U1024::zero());
    /// let a = Fp::new(U1024::from_u64(1), &params);
    /// let d = Fp::new(U1024::from_u64(2), &params);
    /// let curve = EdwardsCurve::new(a, d, &params, &params, U1024::from_u64(0), U1024::from_u64(1));
    ///
    /// let p = curve.identity();
    /// let q = curve.identity();
    /// let r = p.add(&q);
    /// assert_eq!(r, q.add(&p));
    /// ```
    fn add(&self, rhs: &Self) -> Self {
        let a = self.x * rhs.x;
        let b = self.y * rhs.y;

        let t1t2 = self.t * rhs.t;
        let c = self.curve.d * t1t2;

        let d = self.z * rhs.z;

        let x1_plus_y1 = self.x + self.y;
        let x2_plus_y2 = rhs.x + rhs.y;
        let e = (x1_plus_y1 * x2_plus_y2) - a - b;

        let f = d - c;
        let g = d + c;

        let h = b - (self.curve.a * a);

        let x3 = e * f;
        let y3 = g * h;
        let t3 = e * h;
        let z3 = f * g;

        Self {
            x: x3,
            y: y3,
            z: z3,
            t: t3,
            curve: self.curve.clone(),
        }
    }

    /// Computes the Edwards-curve point 2P for this point in extended projective coordinates.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::fp::Fp;
    /// use curvelib::models::twisted_edwards::EdwardsCurve;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = MontgomeryParams::new(U1024::from_u64(17), U1024::zero());
    /// let a = Fp::new(U1024::from_u64(1), &params);
    /// let d = Fp::new(U1024::from_u64(2), &params);
    /// let curve = EdwardsCurve::new(a, d, &params, &params, U1024::from_u64(0), U1024::from_u64(1));
    ///
    /// let point = curve.identity();
    /// let doubled = point.double();
    /// assert_eq!(doubled, point.add(&point));
    /// ```
    fn double(&self) -> Self {
        let a = self.x * self.x;
        let b = self.y * self.y;
        let two = Fp::new(U1024::from_u64(2), self.curve.params);
        let c = two * (self.z * self.z);

        let d = self.curve.a * a;

        let x_plus_y = self.x + self.y;
        let e = (x_plus_y * x_plus_y) - a - b;

        let g = d + b;
        let f = g - c;
        let h = d - b;

        let x3 = e * f;
        let y3 = g * h;
        let t3 = e * h;
        let z3 = f * g;

        Self {
            x: x3,
            y: y3,
            z: z3,
            t: t3,
            curve: self.curve.clone(),
        }
    }

    /// Convert this projective/extended point into affine (x, y) coordinates.
    ///
    /// Returns the affine coordinates corresponding to this point. If the point is the projective
    /// identity (z == 0), this returns the affine identity `(0, 1)`.
    ///
    /// # Returns
    ///
    /// A tuple `(x_aff, y_aff)` containing the affine `x` and `y` coordinates.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::fp::Fp;
    /// use curvelib::models::twisted_edwards::EdwardsCurve;
    /// use curvelib::traits::{Curve, ProjectivePoint, Field};
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = MontgomeryParams::new(U1024::from_u64(17), U1024::zero());
    /// let a = Fp::new(U1024::from_u64(1), &params);
    /// let d = Fp::new(U1024::from_u64(2), &params);
    /// let curve = EdwardsCurve::new(a, d, &params, &params, U1024::from_u64(0), U1024::from_u64(1));
    ///
    /// let pt = curve.identity();
    /// let (x, y) = pt.to_affine();
    /// assert_eq!(x, Fp::zero(curve.params));
    /// assert_eq!(y, Fp::one(curve.params));
    /// ```
    fn to_affine(&self) -> (Fp<'a>, Fp<'a>) {
        if self.z.value == U1024::zero() {
            // Access via deref
            let zero = Fp::zero(self.curve.params);
            let one = Fp::one(self.curve.params);
            return (zero, one);
        }

        let z_inv = self.z.inv().unwrap();
        let x_aff = self.x * z_inv;
        let y_aff = self.y * z_inv;
        (x_aff, y_aff)
    }

    /// Multiply the point by a 1024-bit scalar using a binary double-and-add algorithm.
    ///
    /// Returns the scalar multiple of the receiver as a new `TePoint<'a>` on the same curve.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::fp::Fp;
    /// use curvelib::models::twisted_edwards::EdwardsCurve;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use mathlib::field::montgomery::MontgomeryParams;
    /// use mathlib::{U1024, BigInt};
    ///
    /// let params = MontgomeryParams::new(U1024::from_u64(17), U1024::zero());
    /// let a = Fp::new(U1024::from_u64(1), &params);
    /// let d = Fp::new(U1024::from_u64(2), &params);
    /// let curve = EdwardsCurve::new(a, d, &params, &params, U1024::from_u64(0), U1024::from_u64(1));
    ///
    /// let p = curve.identity();
    /// let k = U1024::zero();
    /// let r = p.mul(&k);
    /// assert!(r.is_identity());
    /// ```
    fn mul(&self, scalar: &U1024) -> Self {
        let mut res = self.curve.identity();

        for i in (0..1024).rev() {
            res = res.double();
            if scalar.bit(i) {
                res = res.add(self);
            }
        }
        res
    }
}
