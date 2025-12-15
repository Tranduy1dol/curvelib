use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{BigInt, FieldElement, U1024};

use crate::{
    algebra::fields::{Fp, Fp2, Fp6},
    models::{sextic_twist::STPoint as G2Projective, short_weierstrass::SWPoint as G1Affine},
    traits::{Field, ProjectivePoint},
};

fn calculate_slope<'a>(t: &G1Affine<'a>, p: &G1Affine<'a>) -> Fp<'a> {
    let (x1, y1) = t.to_affine();
    let (x2, y2) = p.to_affine();

    if x1 == x2 && y1 == y2 {
        // Tangent: lambda = (3x1^2 + a) / 2y1
        let three = Fp::new(U1024::from_u64(3), t.curve.params);
        let two = Fp::new(U1024::from_u64(2), t.curve.params);

        let num = (x1 * x1 * three) + t.curve.a;
        let den = y1 * two;

        // If y1 = 0 (vertical tangent) -> lambda infinite -> handle separately
        if den == Fp::zero(t.curve.params) {
            return Fp::zero(t.curve.params);
        } // Simplify

        num * den.inv().unwrap()
    } else {
        // Chord: lambda = (y2 - y1) / (x2 - x1)
        let num = y2 - y1;
        let den = x2 - x1;

        if den == Fp::zero(t.curve.params) {
            return Fp::zero(t.curve.params);
        } // Vertical line

        num * den.inv().unwrap()
    }
}

/// Evaluate line function l_{T, P}(Q)
/// Returns an element of Fp6
fn evaluate_line<'a>(t: &G1Affine<'a>, p: &G1Affine<'a>, q: &G2Projective<'a>) -> Fp6<'a> {
    // Line: y - y1 - lambda(x - x1) = 0
    // => Val = y_Q - y_T - lambda * (x_Q - x_T)
    // Note: T, P in G1 (Fp), Q in G2 (Fp2).
    // Lambda in Fp.
    // Result in Fp2 (since Q is Fp2), then embed into Fp6.

    let lambda = calculate_slope(t, p); // Fp

    let (xt_aff, yt_aff) = t.to_affine();
    let xt = xt_aff;
    let yt = yt_aff;

    // Q coordinates (Fp2)
    // Assume G2Point stores Affine, or we convert here
    // For simplicity given G2Point doesn't have affine cache, we use x, y directly assuming z=1 or projective calc
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

    // Calculate: (yq - yt) - lambda * (xq - xt)
    // Need to lift yt, xt, lambda from Fp to Fp2 to compute with xq, yq
    let zero = FieldElement::zero(t.curve.params);

    // Embed Fp -> Fp2 (imaginary part = 0)
    let yt_fp2 = Fp2::new(yt, Fp::from(zero));
    let xt_fp2 = Fp2::new(xt, Fp::from(zero));
    let lambda_fp2 = Fp2::new(lambda, Fp::from(zero));

    let term1 = yq - yt_fp2;
    let term2 = lambda_fp2 * (xq - xt_fp2);
    let res_fp2 = term1 - term2;

    // Embed Fp2 -> Fp6
    // Fp6 = c0 + c1 v + c2 v^2. Embed into c0 (or other pos depending on Twist).
    // For standard Tate, this value is at c0.
    let z_fp2 = Fp2::new(Fp::from(zero), Fp::from(zero));
    Fp6::new(res_fp2, z_fp2, z_fp2)
}

/// Tính f_{r, P}(Q)
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

pub fn generate_xi_fp2(params: &MontgomeryParams) -> Fp2<'_> {
    let z = FieldElement::zero(params);
    let o = FieldElement::one(params);
    Fp2::new(Fp::from(z), Fp::from(o))
}
