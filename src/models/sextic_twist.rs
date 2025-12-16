use mathlib::{FieldElement, MontgomeryParams, fp, u1024};

use crate::algebra::fields::{Fp, Fp2};
use crate::def_weierstrass_curve;
use crate::traits::{Curve, Field, ProjectivePoint};

def_weierstrass_curve!(SexticTwist, Fp2<'a>);

impl<'a> Curve<'a> for SexticTwist<'a> {
    type Point = STPoint<'a>;

    /// Returns the curve's point at infinity represented in Jacobian coordinates.
    ///
    /// The returned point uses Fp2 coordinates x = (1, 0), y = (1, 0), and z = (0, 0)
    /// and is associated with this curve instance.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::{Fp, Fp2};
    /// use curvelib::models::SexticTwist;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use mathlib::{mont, BigInt, FieldElement, u1024, MontgomeryParams};
    ///
    /// // Build a tiny, self-contained SexticTwist instance over F_13.
    /// let params = mont!(u1024!(13), u1024!(0));
    /// let zero = Fp::from(FieldElement::zero(&params));
    /// let one = Fp::from(FieldElement::one(&params));
    ///
    /// let a = Fp2::new(zero, zero); // a = 0
    /// let b = Fp2::new(one, zero);  // b = 1
    /// let gx = Fp2::new(zero, zero);
    /// let gy = Fp2::new(one, zero);
    ///
    /// let curve = SexticTwist::new(a, b, &params, &params, gx, gy);
    /// let inf = curve.identity();
    /// assert!(inf.is_identity());
    /// ```
    fn identity(&self) -> Self::Point {
        // Create base field elements for constructing Fp2 coordinates
        let zero_fp = FieldElement::zero(self.params);
        let one_fp = FieldElement::one(self.params);

        // Build Fp2 elements: zero = (0, 0), one = (1, 0)
        let zero_fp2 = Fp2::new(Fp::from(zero_fp), Fp::from(zero_fp));
        let one_fp2 = Fp2::new(Fp::from(one_fp), Fp::from(zero_fp)); // 1 + 0u

        // Point at infinity in Jacobian: any point with Z = 0
        // We use (1, 1, 0) as a canonical representation
        STPoint {
            x: one_fp2,
            y: one_fp2,
            z: zero_fp2, // Z = 0 signals point at infinity
            curve: self.clone(),
        }
    }

    /// Checks whether an affine point (x, y) satisfies the curve equation Y^2 = X^3 + aX + b in Fp2.
    ///
    /// # Returns
    ///
    /// `true` if the provided affine coordinates satisfy the curve equation, `false` otherwise.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::{Fp, Fp2};
    /// use curvelib::models::SexticTwist;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use mathlib::{mont, BigInt, FieldElement, u1024, MontgomeryParams};
    ///
    /// let params = mont!(u1024!(13), u1024!(0));
    /// let zero = Fp::from(FieldElement::zero(&params));
    /// let one = Fp::from(FieldElement::one(&params));
    ///
    /// // Curve: y^2 = x^3 + 1, and choose affine point (0, 1).
    /// let a = Fp2::new(zero, zero);
    /// let b = Fp2::new(one, zero);
    /// let gx = Fp2::new(zero, zero);
    /// let gy = Fp2::new(one, zero);
    /// let curve = SexticTwist::new(a, b, &params, &params, gx, gy);
    ///
    /// let (x, y) = curve.generator().to_affine();
    /// assert!(curve.is_on_curve(&x, &y));
    /// ```
    fn is_on_curve(
        &self,
        x: &<Self::Point as ProjectivePoint<'a>>::Field,
        y: &<Self::Point as ProjectivePoint<'a>>::Field,
    ) -> bool {
        // Check affine curve equation: Y^2 = X^3 + aX + b (all operations in Fp2)
        let y2 = *y * *y; // Y^2
        let x2 = *x * *x; // X^2
        let x3 = x2 * *x; // X^3
        let ax = self.a * *x; // aX

        let rhs = x3 + ax + self.b; // X^3 + aX + b
        y2 == rhs // Compare LHS and RHS
    }

    /// Access the curve's Montgomery scalar parameters.
    ///
    /// # Returns
    ///
    /// A reference to the `MontgomeryParams` used for scalar-field arithmetic.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::{Fp, Fp2};
    /// use curvelib::models::SexticTwist;
    /// use curvelib::traits::Curve;
    /// use mathlib::{mont, BigInt, FieldElement, u1024, MontgomeryParams};
    ///
    /// let params = mont!(u1024!(13), u1024!(0));
    /// let zero = Fp::from(FieldElement::zero(&params));
    /// let one = Fp::from(FieldElement::one(&params));
    /// let curve = SexticTwist::new(
    ///     Fp2::new(zero, zero),
    ///     Fp2::new(one, zero),
    ///     &params,
    ///     &params,
    ///     Fp2::new(zero, zero),
    ///     Fp2::new(one, zero),
    /// );
    /// let sp = curve.scalar_params();
    /// assert_eq!(sp.modulus, u1024!(13));
    /// ```
    fn scalar_params(&self) -> &'a MontgomeryParams {
        self.scalar_params
    }

    /// Returns the curve's generator point expressed in Jacobian-like coordinates.
    ///
    /// The returned `STPoint` uses the curve's stored generator X and Y values and a Z coordinate
    /// representing the affine point (Z = 1 in Fp2).
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::{Fp, Fp2};
    /// use curvelib::models::SexticTwist;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use mathlib::{mont, BigInt, FieldElement, u1024, MontgomeryParams};
    ///
    /// let params = mont!(u1024!(13), u1024!(0));
    /// let zero = Fp::from(FieldElement::zero(&params));
    /// let one = Fp::from(FieldElement::one(&params));
    ///
    /// let curve = SexticTwist::new(
    ///     Fp2::new(zero, zero),
    ///     Fp2::new(one, zero),
    ///     &params,
    ///     &params,
    ///     Fp2::new(zero, zero),
    ///     Fp2::new(one, zero),
    /// );
    /// let g = curve.generator();
    /// assert!(!g.is_identity());
    /// ```
    fn generator(&self) -> Self::Point {
        let x = self.generator_x;
        let y = self.generator_y;

        // Z = (1, 0) in Fp2 represents affine coordinates (Z = 1)
        let z = Fp2::new(
            Fp::from(FieldElement::one(self.params)),
            Fp::from(FieldElement::zero(self.params)),
        );

        // Return generator in Jacobian form with Z = 1
        STPoint {
            x,
            y,
            z,
            curve: self.clone(),
        }
    }
}

// --- 2. G2 POINT DEFINITION (Jacobian Coordinates) ---
#[derive(Clone, Debug)]
pub struct STPoint<'a> {
    pub x: Fp2<'a>,
    pub y: Fp2<'a>,
    pub z: Fp2<'a>,
    pub curve: SexticTwist<'a>,
}

impl<'a> PartialEq for STPoint<'a> {
    /// Determines whether two projective points represent the same point on the curve.
    ///
    /// Points with Z == 0 are treated as the identity; otherwise the points are compared by
    /// converting both to affine coordinates and checking that their X and Y components are equal.
    ///
    /// # Returns
    /// `true` if both points are the identity or their affine X and Y coordinates are equal, `false` otherwise.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::{Fp, Fp2};
    /// use curvelib::models::{SexticTwist, STPoint};
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use mathlib::{mont, BigInt, FieldElement, u1024, MontgomeryParams};
    ///
    /// let params = mont!(u1024!(13), u1024!(0));
    /// let zero = Fp::from(FieldElement::zero(&params));
    /// let one = Fp::from(FieldElement::one(&params));
    /// let curve = SexticTwist::new(
    ///     Fp2::new(zero, zero),
    ///     Fp2::new(one, zero),
    ///     &params,
    ///     &params,
    ///     Fp2::new(zero, zero),
    ///     Fp2::new(one, zero),
    /// );
    ///
    /// let p = curve.generator();
    /// let q = STPoint::new(p.x, p.y, p.z, curve.clone());
    /// assert_eq!(p, q);
    /// assert!(p.eq(&q));
    /// ```
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
impl<'a> Eq for STPoint<'a> {}

impl<'a> STPoint<'a> {
    /// Constructs an `STPoint` from the given projective (Jacobian) coordinates and associates it with `curve`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::{Fp, Fp2};
    /// use curvelib::models::{SexticTwist, STPoint};
    /// use mathlib::{mont, BigInt, FieldElement, u1024, MontgomeryParams};
    ///
    /// let params = mont!(u1024!(13), u1024!(0));
    /// let zero = Fp::from(FieldElement::zero(&params));
    /// let one = Fp::from(FieldElement::one(&params));
    /// let curve = SexticTwist::new(
    ///     Fp2::new(zero, zero),
    ///     Fp2::new(one, zero),
    ///     &params,
    ///     &params,
    ///     Fp2::new(zero, zero),
    ///     Fp2::new(one, zero),
    /// );
    ///
    /// let x = Fp2::new(one, zero);
    /// let y = Fp2::new(one, zero);
    /// let z = Fp2::new(zero, zero);
    /// let p = STPoint::new(x, y, z, curve);
    /// let _ = p;
    /// ```
    pub fn new(x: Fp2<'a>, y: Fp2<'a>, z: Fp2<'a>, curve: SexticTwist<'a>) -> Self {
        Self { x, y, z, curve }
    }

    /// Computes the additive inverse of this point on the curve.
    ///
    /// The inverse of a point P = (X:Y:Z) in Jacobian coordinates is -P = (X:-Y:Z).
    /// If `self` is the identity (point at infinity), a clone of `self` is returned.
    ///
    /// # Returns
    ///
    /// A new `STPoint` representing the additive inverse of `self`.
    pub fn neg(&self) -> Self {
        // Identity is its own negation
        if self.is_identity() {
            return self.clone();
        }
        // Negate a point in Jacobian: -P = (X, -Y, Z)
        // Only the Y coordinate changes sign
        let neg_y = -self.y; // Fp2 implements Neg trait
        Self {
            x: self.x,
            y: neg_y,
            z: self.z,
            curve: self.curve.clone(),
        }
    }
}

impl<'a> ProjectivePoint<'a> for STPoint<'a> {
    type Field = Fp2<'a>;

    /// Checks whether this point is the curve identity (point at infinity).
    ///
    /// The point is considered the identity when its projective Z coordinate is zero in both Fp2 components.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::{Fp, Fp2};
    /// use curvelib::models::SexticTwist;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use mathlib::{mont, BigInt, FieldElement, u1024, MontgomeryParams};
    ///
    /// let params = mont!(u1024!(13), u1024!(0));
    /// let zero = Fp::from(FieldElement::zero(&params));
    /// let one = Fp::from(FieldElement::one(&params));
    /// let curve = SexticTwist::new(
    ///     Fp2::new(zero, zero),
    ///     Fp2::new(one, zero),
    ///     &params,
    ///     &params,
    ///     Fp2::new(zero, zero),
    ///     Fp2::new(one, zero),
    /// );
    /// let inf = curve.identity();
    /// assert!(inf.is_identity());
    /// ```
    fn is_identity(&self) -> bool {
        // A point is at infinity iff Z = 0 in Fp2
        // Check both components of Z are zero
        self.z.c0.value == u1024!(0) && self.z.c1.value == u1024!(0)
    }

    /// Adds two points on the curve and returns their sum in Jacobian coordinates.
    ///
    /// Performs elliptic-curve addition with correct handling of the identity element; if the
    /// inputs are inverses the result is the curve identity. The returned point is associated
    /// with the same curve as `self`.
    ///
    /// # Returns
    /// The curve point representing `self + rhs` in Jacobian coordinates.
    fn add(&self, rhs: &Self) -> Self {
        // Handle identity cases
        if self.is_identity() {
            return rhs.clone();
        }
        if rhs.is_identity() {
            return self.clone();
        }

        // Standard Jacobian addition formulas for elliptic curves
        // See "Guide to Elliptic Curve Cryptography" Algorithm 3.22
        let z1z1 = self.z.square(); // Z1^2
        let z2z2 = rhs.z.square(); // Z2^2

        let u1 = self.x * z2z2; // U1 = X1 * Z2^2
        let u2 = rhs.x * z1z1; // U2 = X2 * Z1^2

        let s1 = self.y * (rhs.z * z2z2); // S1 = Y1 * Z2^3
        let s2 = rhs.y * (self.z * z1z1); // S2 = Y2 * Z1^3

        // Check if points have same X coordinate (in affine)
        if u1 == u2 {
            return if s1 == s2 {
                // Same point: use doubling formula
                self.double()
            } else {
                // Inverse points: return identity
                self.curve.identity()
            };
        }

        // Compute new coordinates
        let h = u2 - u1; // H = U2 - U1
        let r = s2 - s1; // R = S2 - S1
        let hh = h.square(); // HH = H^2
        let hhh = hh * h; // HHH = H^3
        let v = u1 * hh; // V = U1 * HH

        // Create constant 2 in Fp2
        let two_val = u1024!(2);
        let zero_fp = FieldElement::zero(self.curve.params);
        let two_fp = fp!(two_val, self.curve.params);
        let two_fp2 = Fp2::new(Fp::from(two_fp), Fp::from(zero_fp));

        // X3 = R^2 - HHH - 2V
        let x3 = r.square() - hhh - (two_fp2 * v);
        // Y3 = R(V - X3) - S1 * HHH
        let y3 = (r * (v - x3)) - (s1 * hhh);
        // Z3 = Z1 * Z2 * H
        let z3 = self.z * rhs.z * h;

        Self {
            x: x3,
            y: y3,
            z: z3,
            curve: self.curve.clone(),
        }
    }

    /// Doubles this point using Jacobian-coordinate doubling on the sextic-twist curve.
    ///
    /// Uses the Jacobian doubling formulas over Fp2 to compute 2P and returns the resulting point
    /// in the same projective representation.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// // Given an existing STPoint `p` (for example obtained from a curve's generator),
    /// // compute its double.
    /// // let curve = SexticTwist::some_constructor(...);
    /// // let p = curve.generator();
    /// // let r = p.double();
    /// // assert_eq!(r, p.add(&p));
    /// ```
    fn double(&self) -> Self {
        // Identity doubled is still identity
        if self.is_identity() {
            return self.clone();
        }

        // Jacobian point doubling formulas (adapted for Fp2 operations)
        // See "Explicit-Formulas Database" - dbl-2007-bl for similar formulas
        let xx = self.x.square(); // XX = X^2
        let yy = self.y.square(); // YY = Y^2
        let yyyy = yy.square(); // YYYY = YY^2
        let zz = self.z.square(); // ZZ = Z^2

        // Construct small integer constants in Fp2
        let zero_fp = FieldElement::zero(self.curve.params);
        let two_val = u1024!(2);
        let two_fp = fp!(two_val, self.curve.params);
        let two_fp2 = Fp2::new(Fp::from(two_fp), Fp::from(zero_fp));

        // S = 4XY^2 = 2 * ((X * YY) * 2)
        let s = two_fp2 * ((self.x * yy) * two_fp2);

        // M = 3X^2 + aZ^4 (slope of tangent line)
        let three_val = u1024!(3);
        let three_fp = fp!(three_val, self.curve.params);
        let three_fp2 = Fp2::new(Fp::from(three_fp), Fp::from(zero_fp));

        let zzzz = zz.square(); // Z^4
        let m = (three_fp2 * xx) + (self.curve.a * zzzz);

        // X' = M^2 - 2S
        let x_new = m.square() - (s * two_fp2);

        // Z' = 2YZ (alternative: (Y+Z)^2 - YY - ZZ)
        let z_new = (self.y * self.z) * two_fp2;

        // Y' = M(S - X') - 8YYYY
        let eight_val = u1024!(8);
        let eight_fp = fp!(eight_val, self.curve.params);
        let eight_fp2 = Fp2::new(Fp::from(eight_fp), Fp::from(zero_fp));

        let t = eight_fp2 * yyyy;
        let y_new = (m * (s - x_new)) - t;

        Self {
            x: x_new,
            y: y_new,
            z: z_new,
            curve: self.curve.clone(),
        }
    }

    /// Converts this Jacobian-coordinate point into affine coordinates over Fp2.
    ///
    /// If the point is the identity (point-at-infinity), it returns (0, 0) in Fp2.
    /// Otherwise, returns the pair (x_aff, y_aff) computed by multiplying X by Z^{-2}
    /// and Y by Z^{-3}, where Z^{-1} is the multiplicative inverse of Z in Fp2.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::{Fp, Fp2};
    /// use curvelib::models::SexticTwist;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use mathlib::{mont, BigInt, FieldElement, u1024, MontgomeryParams};
    ///
    /// let params = mont!(u1024!(13), u1024!(0));
    /// let zero = Fp::from(FieldElement::zero(&params));
    /// let one = Fp::from(FieldElement::one(&params));
    /// let curve = SexticTwist::new(
    ///     Fp2::new(zero, zero),
    ///     Fp2::new(one, zero),
    ///     &params,
    ///     &params,
    ///     Fp2::new(zero, zero),
    ///     Fp2::new(one, zero),
    /// );
    ///
    /// let pt = curve.identity();
    /// let (x_aff, y_aff) = pt.to_affine();
    /// assert_eq!(x_aff, Fp2::new(zero, zero));
    /// assert_eq!(y_aff, Fp2::new(zero, zero));
    /// ```
    fn to_affine(&self) -> (Self::Field, Self::Field) {
        // Point at infinity maps to (0, 0) in affine
        if self.is_identity() {
            let zero = FieldElement::zero(self.curve.params);
            let zero_fp2 = Fp2::new(Fp::from(zero), Fp::from(zero));
            return (zero_fp2, zero_fp2);
        }

        // Convert Jacobian (X:Y:Z) to affine (x, y) where x = X/Z^2, y = Y/Z^3
        let z_inv = self.z.inv().unwrap(); // Compute Z^(-1) in Fp2
        let z2_inv = z_inv.square(); // Z^(-2)
        let z3_inv = z2_inv * z_inv; // Z^(-3)

        let x_aff = self.x * z2_inv; // x = X * Z^(-2)
        let y_aff = self.y * z3_inv; // y = Y * Z^(-3)
        (x_aff, y_aff)
    }

    /// Computes the scalar multiple of this point.
    ///
    /// The result is equivalent to adding this point to itself `scalar` times.
    ///
    /// # Returns
    ///
    /// The point resulting from multiplying this point by `scalar`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::algebra::fields::{Fp, Fp2};
    /// use curvelib::models::SexticTwist;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use mathlib::{mont, BigInt, FieldElement, u1024, MontgomeryParams};
    ///
    /// let params = mont!(u1024!(13), u1024!(0));
    /// let zero = Fp::from(FieldElement::zero(&params));
    /// let one = Fp::from(FieldElement::one(&params));
    /// let curve = SexticTwist::new(
    ///     Fp2::new(zero, zero),
    ///     Fp2::new(one, zero),
    ///     &params,
    ///     &params,
    ///     Fp2::new(zero, zero),
    ///     Fp2::new(one, zero),
    /// );
    ///
    /// let p = curve.generator();
    /// let k = u1024!(3);
    /// let r = p.mul(&k);
    /// // group law sanity: 3P == P + P + P
    /// assert_eq!(r, p.add(&p).add(&p));
    /// ```
    fn mul(&self, scalar: &mathlib::U1024) -> Self {
        // Double-and-add scalar multiplication algorithm
        // Start with identity (0 * P = O)
        let mut res = self.curve.identity();
        let num_bits = scalar.bits();

        // Process bits from most significant to least significant
        for i in (0..num_bits).rev() {
            res = res.double(); // Double accumulated result
            if scalar.bit(i) {
                // If current bit is set
                res = res.add(self); // Add P to accumulator
            }
        }
        res
    }
}
