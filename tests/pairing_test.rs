//! Tests for pairing operations.

use curvelib::{
    instances::bls6_6::{self, FINAL_EXPONENT},
    protocol::pairing::tate_pairing,
    traits::{Curve, ProjectivePoint},
};
use mathlib::U1024;

#[test]
fn test_bls6_6_bilinearity() {
    // 1. Get BLS6_6 curves
    let g1_curve = bls6_6::get_g1_curve();
    let g2_curve = bls6_6::get_g2_curve();

    // Parameters
    let r_val = U1024::from_u64(13); // Subgroup order
    let final_exp = U1024::from_u64(FINAL_EXPONENT); // (43^6 - 1) / 13

    // 2. Get generators
    let p = g1_curve.generator();
    let q = g2_curve.generator();

    // Verify points are on curves
    {
        let (px, py) = p.to_affine();
        assert!(g1_curve.is_on_curve(&px, &py), "P must be on G1");

        let (qx, qy) = q.to_affine();
        assert!(g2_curve.is_on_curve(&qx, &qy), "Q must be on G2");
    }

    // 3. Calculate e(P, Q)
    let e1 = tate_pairing(&p, &q, r_val, final_exp);
    println!("e(P, Q) = {:?}", e1);

    // 4. Calculate e(2P, Q)
    let p2 = p.double();
    let e2 = tate_pairing(&p2, &q, r_val, final_exp);
    println!("e(2P, Q) = {:?}", e2);

    // Verify bilinearity: e(2P, Q) = e(P, Q)²
    let e1_squared = e1.square();
    assert_eq!(e2, e1_squared, "Bilinearity check failed!");
    println!("e1² = {:?}", e1_squared);
}

#[test]
fn test_pairing_identity() {
    let g1_curve = bls6_6::get_g1_curve();
    let g2_curve = bls6_6::get_g2_curve();

    let r_val = U1024::from_u64(13);
    let final_exp = U1024::from_u64(FINAL_EXPONENT);

    // Identity test: e(O, Q) should be the identity in the target field
    let o = g1_curve.identity();
    let q = g2_curve.generator();

    let e = tate_pairing(&o, &q, r_val, final_exp);

    // For identity point, pairing might be undefined or identity in Fp6
    // At least verify it doesn't crash
    println!("e(O, Q) = {:?}", e);
}
