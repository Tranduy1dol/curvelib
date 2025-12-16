use mathlib::{MontgomeryParams, U1024, fp, u1024};

use crate::{
    algebra::fields::{Fp, Fp2, Fp6},
    models::{sextic_twist::STPoint as G2Projective, short_weierstrass::SWPoint as G1Affine},
    traits::{Field, ProjectivePoint},
};

/// Compute the slope (λ) of the line through two G1 affine points for use in line evaluations.
///
/// This returns the curve slope for the tangent case when `t == p` (λ = (3*x^2 + a) / (2*y)) or for
/// the chord case when `t != p` (λ = (y2 - y1) / (x2 - x1)). If the denominator is zero (vertical
/// tangent or vertical line), the function returns `None` to signal a vertical line.
fn calculate_slope<'a>(t: &G1Affine<'a>, p: &G1Affine<'a>) -> Option<Fp<'a>> {
    let (x1, y1) = t.to_affine();
    let (x2, y2) = p.to_affine();

    if x1 == x2 && y1 == y2 {
        // Tangent: lambda = (3x1^2 + a) / 2y1
        let three = fp!(u1024!(3), t.curve.params);
        let two = fp!(u1024!(2), t.curve.params);

        let num = (x1 * x1 * three) + t.curve.a;
        let den = y1 * two;

        // If y1 = 0 (vertical tangent) -> return None to signal vertical line
        if den == Fp::zero(t.curve.params) {
            return None;
        }

        Some(num * Field::inv(&den).unwrap())
    } else {
        // Chord: lambda = (y2 - y1) / (x2 - x1)
        let num = y2 - y1;
        let den = x2 - x1;

        // If x1 == x2 (vertical line) -> return None
        if den == Fp::zero(t.curve.params) {
            return None;
        }

        Some(num * Field::inv(&den).unwrap())
    }
}

/// Evaluate the line function l_{T,P} at a G2 point Q and embed the result into Fp6.
///
/// Computes the value of the line defined by T and P (in G1) evaluated at Q (in G2).
/// For a non-vertical line: y_Q - y_T - lambda * (x_Q - x_T).
/// For a vertical line: x_Q - x_T.
/// Coordinates from Fp are lifted into Fp2 to match Q, and the resulting Fp2 value
/// is placed in the c0 coefficient of the returned Fp6 with c1 and c2 set to zero.
/// If Q is projective (z ≠ 1) it is converted to affine first.
fn evaluate_line<'a>(t: &G1Affine<'a>, p: &G1Affine<'a>, q: &G2Projective<'a>) -> Fp6<'a> {
    let (xt_aff, yt_aff) = t.to_affine();
    let xt = xt_aff;
    let yt = yt_aff;

    // Convert Q to affine coordinates
    let (xq, yq) = if q.z == Fp2::one(t.curve.params) {
        (q.x, q.y)
    } else {
        let z_inv =
            q.z.inv()
                .expect("Q z-coordinate is zero in miller loop (point at infinity?)");
        let z2 = z_inv.square();
        let z3 = z2 * z_inv;
        (q.x * z2, q.y * z3)
    };

    let zero = Fp::zero(t.curve.params);
    let z_fp2 = Fp2::new(zero, zero);

    // Pattern-match on the slope to handle vertical vs. non-vertical lines
    let lambda_opt = calculate_slope(t, p);

    let res_fp2 = match lambda_opt {
        Some(lambda) => {
            // Non-vertical line: (yq - yt) - lambda * (xq - xt)
            // Embed Fp -> Fp2 (imaginary part = 0)
            let yt_fp2 = Fp2::new(yt, zero);
            let xt_fp2 = Fp2::new(xt, zero);
            let lambda_fp2 = Fp2::new(lambda, zero);

            let term1 = yq - yt_fp2;
            let term2 = lambda_fp2 * (xq - xt_fp2);
            term1 - term2
        }
        None => {
            // Vertical line: xq - xt
            let xt_fp2 = Fp2::new(xt, zero);
            xq - xt_fp2
        }
    };

    // Embed Fp2 -> Fp6
    // Fp6 = c0 + c1 v + c2 v^2. Embed into c0.
    Fp6::new(res_fp2, z_fp2, z_fp2)
}

/// Computes the Miller function value f_{r,P}(Q) for the given G1 point `p`, G2 point `q`, and scalar `r_order`.
///
/// Returns the resulting element of Fp6 representing f_{r,P}(Q) computed by the Miller loop.
pub fn miller_loop<'a>(p: &G1Affine<'a>, q: &G2Projective<'a>, r_order: U1024) -> Fp6<'a> {
    let params = p.curve.params;

    let mut t = p.clone();
    let mut f = Fp6::one(params);

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

/// Constructs the Fp2 element ξ = (0, 1) using the provided Montgomery parameters.
///
/// The returned value has its first (real) component set to `0` and its second (imaginary) component set to `1`.
///
/// # Examples
///
/// ```
/// use curvelib::{
///     algebra::fields::{Fp, Fp2},
///     protocol::pairing::miller::generate_xi_fp2,
/// };
/// use mathlib::{u1024, mont};
///
/// let p_val = u1024!(43);
/// let params = mont!(p_val, u1024!(0));
///
/// let xi = generate_xi_fp2(&params);
/// let expected = Fp2::new(
///     Fp::zero(&params),
///     Fp::one(&params)
/// );
/// assert_eq!(xi, expected);
/// ```
pub fn generate_xi_fp2(params: &MontgomeryParams) -> Fp2<'_> {
    Fp2::new(Fp::zero(params), Fp::one(params))
}
