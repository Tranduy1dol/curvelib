//! KZG integration tests using the PolynomialCommitment trait.
//!
//! Tests the trait-based KZG implementation with mathlib 1.1.0.

use curvelib::{
    instances::bls6_6::{self, Bls6_6BaseField, FINAL_EXPONENT},
    protocol::{PolynomialCommitment, commitment::kzg::Kzg},
    traits::{Curve, ProjectivePoint},
};
use mathlib::{FieldElement, Polynomial, U1024, fp};

#[test]
fn test_kzg_setup_trait() {
    println!("🚀 KZG Setup Test (Trait-based)");

    let g1_curve = bls6_6::get_g1_curve();
    let g2_curve = bls6_6::get_g2_curve();

    let tau = FieldElement::<Bls6_6BaseField>::new(U1024::from_u64(5));
    let max_degree = 3;
    let r_order = U1024::from_u64(13);
    let final_exp = U1024::from_u64(FINAL_EXPONENT);

    let params =
        Kzg::<Bls6_6BaseField>::setup(tau, max_degree, &g1_curve, &g2_curve, r_order, final_exp);

    assert_eq!(params.powers_of_tau.len(), 4);

    // Verify G2 points are on curve
    let (qx, qy) = params.g2.to_affine();
    assert!(g2_curve.is_on_curve(&qx, &qy), "G2 must be on curve");

    println!("✅ KZG Setup passed!");
}

#[test]
fn test_kzg_commit_trait() {
    println!("🔒 KZG Commit Test (Trait-based)");

    let g1_curve = bls6_6::get_g1_curve();
    let g2_curve = bls6_6::get_g2_curve();

    let tau = FieldElement::<Bls6_6BaseField>::new(U1024::from_u64(5));
    let r_order = U1024::from_u64(13);
    let final_exp = U1024::from_u64(FINAL_EXPONENT);

    let params = Kzg::<Bls6_6BaseField>::setup(tau, 3, &g1_curve, &g2_curve, r_order, final_exp);

    // Create polynomial P(x) = 1 + 2x + 3x²
    let poly = Polynomial::new(vec![
        fp!(1u64, Bls6_6BaseField),
        fp!(2u64, Bls6_6BaseField),
        fp!(3u64, Bls6_6BaseField),
    ]);

    // Use trait method
    let commitment = Kzg::commit(&params, &poly);

    let (cx, cy) = commitment.to_affine();
    println!("Commitment: ({:?}, {:?})", cx.to_u1024(), cy.to_u1024());

    println!("✅ KZG Commit passed!");
}

#[test]
fn test_kzg_open_trait() {
    println!("🔓 KZG Open Test (Trait-based)");

    let g1_curve = bls6_6::get_g1_curve();
    let g2_curve = bls6_6::get_g2_curve();

    let tau = FieldElement::<Bls6_6BaseField>::new(U1024::from_u64(5));
    let r_order = U1024::from_u64(13);
    let final_exp = U1024::from_u64(FINAL_EXPONENT);

    let params = Kzg::<Bls6_6BaseField>::setup(tau, 3, &g1_curve, &g2_curve, r_order, final_exp);

    // Create polynomial P(x) = 1 + 2x + 3x²
    let poly = Polynomial::new(vec![
        fp!(1u64, Bls6_6BaseField),
        fp!(2u64, Bls6_6BaseField),
        fp!(3u64, Bls6_6BaseField),
    ]);

    // Open at z = 2: P(2) = 1 + 4 + 12 = 17
    let z = FieldElement::<Bls6_6BaseField>::new(U1024::from_u64(2));

    // Use trait method
    let (_proof, value) = Kzg::open(&params, &poly, &z);

    println!("Evaluation at z=2: {:?}", value.to_u1024());
    assert_eq!(value.to_u1024(), U1024::from_u64(17));

    println!("✅ KZG Open passed!");
}

#[test]
fn test_kzg_full_flow_trait() {
    println!("🔄 KZG Full Flow Test (Trait-based)");

    let g1_curve = bls6_6::get_g1_curve();
    let g2_curve = bls6_6::get_g2_curve();

    let tau = FieldElement::<Bls6_6BaseField>::new(U1024::from_u64(5));
    let r_order = U1024::from_u64(13);
    let final_exp = U1024::from_u64(FINAL_EXPONENT);

    let params = Kzg::<Bls6_6BaseField>::setup(tau, 3, &g1_curve, &g2_curve, r_order, final_exp);

    // Polynomial P(x) = 1 + 2x + 3x²
    let poly = Polynomial::new(vec![
        fp!(1u64, Bls6_6BaseField),
        fp!(2u64, Bls6_6BaseField),
        fp!(3u64, Bls6_6BaseField),
    ]);

    // Commit using trait
    let commitment = Kzg::commit(&params, &poly);

    // Open at z = 2 using trait
    let z = FieldElement::<Bls6_6BaseField>::new(U1024::from_u64(2));
    let (proof, value) = Kzg::open(&params, &poly, &z);

    println!("Committed and opened polynomial");
    println!("Value at z=2: {:?}", value.to_u1024());

    // Note: Verification may not work on toy curve
    let valid = Kzg::verify(&params, &commitment, &z, &value, &proof);
    assert!(valid, "Kzg verify failed");

    println!("✅ KZG Full Flow passed!");
}
