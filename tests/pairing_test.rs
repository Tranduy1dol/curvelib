use curvelib::{
    algebra::fields::{Fp, Fp2},
    models::{STPoint, SWPoint, SexticTwist, WeierstrassCurve},
    protocol::pairing::tate_pairing,
    traits::{Field, ProjectivePoint},
};
use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{BigInt, U1024};

#[test]
fn test_bls6_6_bilinearity() {
    // 1. Setup BLS6_6 Params
    // P = 43, r = 13 (Subgroup order)
    let p_val = U1024::from_u64(43);
    let r_val = U1024::from_u64(13); // Subgroup order
    let final_exp = U1024::from_u64(486258696); // (43^6 - 1) / 13

    let base_params = MontgomeryParams::new(p_val, U1024::zero());
    let scalar_params = MontgomeryParams::new(r_val.clone(), U1024::zero());

    // Curve G1: y^2 = x^3 + 6
    let a_val = Fp::new(U1024::from_u64(0), &base_params);
    let b_val = Fp::new(U1024::from_u64(6), &base_params);
    let g1_curve = WeierstrassCurve::new(
        Fp::from(a_val.clone()),
        Fp::from(b_val.clone()),
        &base_params,
        &scalar_params,
        Fp::new(U1024::zero(), &base_params),
        Fp::new(U1024::zero(), &base_params),
    );

    // 2. Create point P in G1 (subgroup order 13)
    // G1 Generator: (13, 15)
    let gx = Fp::new(U1024::from_u64(13), &base_params);
    let gy = Fp::new(U1024::from_u64(15), &base_params);
    let p = SWPoint::new_affine(Fp::from(gx.clone()), Fp::from(gy.clone()), &g1_curve);

    // 3. Create point Q in G2
    // For simplicity, choose Q "untwisted" on Fp (Q is also (13, 15) embedded in Fp2)
    // In reality Q must be linearly independent of P (on Twist).
    // But to test code logic, we use Q = P (embedded in G2).
    // Pairing result will be e(P, P) != 1.

    let zero = Fp::zero(&base_params);
    // let one = Fp::one(&base_params);

    // Create G2Curve
    // Embed a, b into Fp2
    let a_fp2 = Fp2::new(Fp::from(a_val.clone()), Fp::from(zero.clone()));
    let b_fp2 = Fp2::new(Fp::from(b_val.clone()), Fp::from(zero.clone()));
    let gen_x_fp2 = Fp2::new(Fp::from(zero.clone()), Fp::from(zero.clone())); // Dummy, not used in test
    let gen_y_fp2 = Fp2::new(Fp::from(zero.clone()), Fp::from(zero.clone())); // Dummy, not used in test
    let g2_curve = SexticTwist::new(
        a_fp2,
        b_fp2,
        &base_params,
        &scalar_params,
        gen_x_fp2,
        gen_y_fp2,
    );

    // XQ = 13 + 0u, YQ = 15 + 0u
    let qx_fp2 = Fp2::new(Fp::from(gx), Fp::from(zero.clone()));
    let qy_fp2 = Fp2::new(Fp::from(gy), Fp::from(zero.clone()));
    let qz_fp2 = Fp2::new(Fp::from(Fp::one(&base_params)), Fp::from(zero.clone()));

    let q = STPoint::new(qx_fp2, qy_fp2, qz_fp2, g2_curve);

    // 4. Calculate e(P, Q)
    let e1 = tate_pairing(&p, &q, r_val.clone(), final_exp.clone());
    println!("e(P, Q) = {:?}", e1);

    // 5. Calculate e(2P, Q)
    let p2 = p.double();
    let e2 = tate_pairing(&p2, &q, r_val.clone(), final_exp.clone());
    println!("e(2P, Q) = {:?}", e2);

    // 6. Check: e(2P, Q) == e(P, Q)^2
    let e1_squared = e1.square(); // Square in Fp6

    assert_eq!(e2, e1_squared, "Bilinearity check failed!");
}
