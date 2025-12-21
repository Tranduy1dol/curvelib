//! Miller loop for pairing computation.
//!
//! This module provides the Miller loop algorithm for computing
//! pairings on BLS-style curves.

use mathlib::{FieldConfig, FieldElement, U1024};

use crate::algebra::fields::{Fp2, Fp6};
use crate::models::{TwistPoint, WeierstrassPoint};
use crate::traits::ProjectivePoint;

/// Type aliases for clarity.
pub type G1Affine<C> = WeierstrassPoint<C>;
pub type G2Projective<C> = TwistPoint<C>;

/// Helper to create small field constants.
fn fp_from_u64<C: FieldConfig>(val: u64) -> FieldElement<C> {
    FieldElement::<C>::new(U1024::from_u64(val))
}

/// Calculate the slope for the line through two G1 points.
///
/// Returns None if the line is vertical.
fn calculate_slope<C: FieldConfig>(t: &G1Affine<C>, p: &G1Affine<C>) -> Option<FieldElement<C>> {
    let (x1, y1) = t.to_affine();
    let (x2, y2) = p.to_affine();

    if x1 == x2 && y1 == y2 {
        // Tangent: λ = (3x₁² + a) / (2y₁)
        let three = fp_from_u64::<C>(3);
        let two = fp_from_u64::<C>(2);

        let num = (x1 * x1 * three) + t.curve.a;
        let den = y1 * two;

        if den.is_zero() {
            return None;
        }

        Some(num * den.inv())
    } else {
        // Chord: λ = (y₂ - y₁) / (x₂ - x₁)
        let num = y2 - y1;
        let den = x2 - x1;

        if den.is_zero() {
            return None;
        }

        Some(num * den.inv())
    }
}

/// Evaluate the line function at a G2 point.
fn evaluate_line<C: FieldConfig>(t: &G1Affine<C>, p: &G1Affine<C>, q: &G2Projective<C>) -> Fp6<C> {
    let (xt, yt) = t.to_affine();

    // Convert Q to affine
    let (xq, yq) = if q.z == Fp2::<C>::one() {
        (q.x, q.y)
    } else {
        let z_inv = q.z.inv().expect("Q z-coordinate is zero");
        let z2 = z_inv.square();
        let z3 = z2 * z_inv;
        (q.x * z2, q.y * z3)
    };

    let zero_fp = FieldElement::<C>::zero();
    let zero_fp2 = Fp2::<C>::zero();

    let lambda_opt = calculate_slope(t, p);
    let xi = Fp2::u();

    match lambda_opt {
        Some(lambda) => {
            // Sextic twist map: ψ(xq, yq) = (xq * v², yq * v³) where v³ = ξ
            // Line evaluation: L(X, Y) = (Y - yt) - λ(X - xt)
            // L(ψ(xq, yq)) = (yq * v³ - yt) - λ(xq * v² - xt)
            //              = (yq * ξ - yt + λ*xt) - λ*xq * v²

            let yt_fp2 = Fp2::new(yt, zero_fp);
            let xt_fp2 = Fp2::new(xt, zero_fp);
            let lambda_fp2 = Fp2::new(lambda, zero_fp);

            let c0 = (yq * xi) - yt_fp2 + (lambda_fp2 * xt_fp2);
            let c1 = zero_fp2;
            let c2 = -(lambda_fp2 * xq);

            Fp6::new(c0, c1, c2)
        }
        None => {
            // Vertical line: X - xt = xq * v² - xt
            let xt_fp2 = Fp2::new(xt, zero_fp);
            let c0 = -xt_fp2;
            let c1 = zero_fp2;
            let c2 = xq;

            Fp6::new(c0, c1, c2)
        }
    }
}

/// Computes the Miller function value f_{r,P}(Q).
///
/// # Arguments
///
/// * `p` - A G1 point
/// * `q` - A G2 point
/// * `r_order` - The group order
///
/// # Returns
///
/// The resulting element of Fp6.
pub fn miller_loop<C: FieldConfig>(p: &G1Affine<C>, q: &G2Projective<C>, r_order: U1024) -> Fp6<C> {
    let mut t = p.clone();
    let mut f = Fp6::<C>::one();

    let num_bits = r_order.bits();

    for i in (0..num_bits - 1).rev() {
        let l_tt = evaluate_line(&t, &t, q);
        f = f * f;
        f = f * l_tt;

        t = t.double();

        if r_order.bit(i) {
            let l_tp = evaluate_line(&t, p, q);
            f = f * l_tp;

            t = t.add(p);
        }
    }
    f
}

/// Creates the Fp2 element ξ = (0, 1).
pub fn generate_xi_fp2<C: FieldConfig>() -> Fp2<C> {
    Fp2::u()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instances::bls6_6::Bls6_6BaseField;

    #[test]
    fn test_generate_xi() {
        let xi = generate_xi_fp2::<Bls6_6BaseField>();
        assert!(xi.c0.is_zero());
        assert!(!xi.c1.is_zero());
    }
}
