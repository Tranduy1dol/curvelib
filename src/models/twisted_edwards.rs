use mathlib::field::element::FieldElement;
use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{BigInt, U1024};

use crate::traits::{Curve, ProjectivePoint};

impl<'a> Curve<'a> for EdwardsCurve<'a> {
    type Point = TePoint<'a>;

    /// Constructs the identity (neutral) point for this Edwards curve.
    ///
    /// Returns the identity point in projective/extended coordinates: x = 0, y = 1, z = 1, t = 0, tied to this curve.
    ///
    /// # Examples
    ///
    /// ```
    /// // assume `curve` is an `EdwardsCurve` already available in scope
    /// let id = curve.identity();
    /// let (x, y) = id.to_affine();
    /// assert_eq!(x, FieldElement::zero(curve.params));
    /// assert_eq!(y, FieldElement::one(curve.params));
    /// ```
    fn identity(&self) -> Self::Point {
        let curve = self.clone();
        let zero = FieldElement::zero(self.params);
        let one = FieldElement::one(self.params);
        TePoint {
            x: zero,
            y: one,
            z: one,
            t: zero,
            curve,
        }
    }

    /// Checks whether an affine point (x, y) satisfies the Edwards curve equation.
    ///
    /// Returns `true` if the coordinates satisfy a*x^2 + y^2 = 1 + d*x^2*y^2, `false` otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// // Given a curve `curve` and field elements `x`, `y` already constructed for the same field:
    /// // assert!(curve.is_on_curve(&x, &y));
    /// ```
    fn is_on_curve(&self, x: &FieldElement, y: &FieldElement) -> bool {
        let x2 = *x * *x;
        let y2 = *y * *y;
        let lhs = (self.a * x2) + y2;

        let one = FieldElement::one(self.params);
        let rhs = one + (self.d * x2 * y2);

        lhs == rhs
    }
}

#[derive(Clone, Debug)]
pub struct EdwardsCurve<'a> {
    pub a: FieldElement<'a>,
    pub d: FieldElement<'a>,
    pub params: &'a MontgomeryParams,
}

impl<'a> EdwardsCurve<'a> {
    pub fn new(a: FieldElement<'a>, d: FieldElement<'a>, params: &'a MontgomeryParams) -> Self {
        Self { a, d, params }
    }
}

#[derive(Clone, Debug)]
pub struct TePoint<'a> {
    pub x: FieldElement<'a>,
    pub y: FieldElement<'a>,
    pub z: FieldElement<'a>,
    pub t: FieldElement<'a>,
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
    /// ```
    /// // Given existing `x: FieldElement`, `y: FieldElement`, and `curve: EdwardsCurve`
    /// let p = TePoint::new_affine(x, y, &curve);
    /// let (xa, ya) = p.to_affine();
    /// assert_eq!(xa, x);
    /// assert_eq!(ya, y);
    /// ```
    pub fn new_affine(
        x: FieldElement<'a>,
        y: FieldElement<'a>,
        curve: &'a EdwardsCurve<'a>,
    ) -> Self {
        let z = FieldElement::one(curve.params);
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

impl<'a> ProjectivePoint for TePoint<'a> {
    /// Checks whether this point is the identity element on its Edwards curve.
    ///
    /// The identity is encoded in projective/extended coordinates as `x == 0` and `y == z`.
    ///
    /// # Returns
    ///
    /// `true` if the point is the identity element, `false` otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// // Construct an identity point and verify it.
    /// // `curve` must be a valid `EdwardsCurve` and `FieldElement` helpers available in scope.
    /// let id = TePoint::new_affine(FieldElement::zero(curve.params), FieldElement::one(curve.params), curve.clone());
    /// assert!(id.is_identity());
    /// ```
    fn is_identity(&self) -> bool {
        let zero = FieldElement::zero(self.curve.params);
        self.x == zero && self.y == self.z
    }

    /// Adds two Edwards-curve points in projective/extended coordinates.
    ///
    /// The result is the curve point representing `self + rhs` using the curve's
    /// parameters stored in the points.
    ///
    /// # Examples
    ///
    /// ```
    /// // Given two points `p` and `q` on the same Edwards curve:
    /// let r = p.add(&q);
    /// // addition is commutative for curve points
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

    /// Doubles this point on the Edwards curve.
    ///
    /// Computes the curve point 2P using the point's projective/extended representation.
    ///
    /// # Examples
    ///
    /// ```
    /// let doubled = point.double();
    /// assert_eq!(doubled, point.add(&point));
    /// ```
    fn double(&self) -> Self {
        let a = self.x * self.x;
        let b = self.y * self.y;
        let two = FieldElement::new(U1024::from_u64(2), self.curve.params);
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

    /// Converts this projective/extended point to affine coordinates.
    ///
    /// If the point is the projective identity (z == 0) this returns the affine identity `(0, 1)`.
    ///
    /// # Returns
    ///
    /// A tuple `(x_aff, y_aff)` containing the affine x and y coordinates.
    ///
    /// # Examples
    ///
    /// ```
    /// // Given a projective point `p`, obtain its affine coordinates:
    /// let (x, y) = p.to_affine();
    /// // `x` and `y` are FieldElement values representing the affine point.
    /// ```
    fn to_affine(&self) -> (FieldElement<'a>, FieldElement<'a>) {
        if self.z.value == U1024::zero() {
            let zero = FieldElement::zero(self.curve.params);
            let one = FieldElement::one(self.curve.params);
            return (zero, one);
        }

        let z_inv = self.z.inv();
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
    /// ```
    /// // Given a point `p: TePoint<'_>` and a 1024-bit scalar `k: U1024`,
    /// // compute the scalar multiple and verify the zero-scalar case.
    /// let p = /* existing TePoint<'_'> */ unimplemented!();
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