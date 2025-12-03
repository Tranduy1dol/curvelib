use curvelib::models::short_weierstrass::{SWPoint, WeierstrassCurve};
use curvelib::traits::ProjectivePoint;
use mathlib::field::element::FieldElement;
use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{BigInt, U1024};

#[test]
fn test_tiny_curve_operations() {
    let mut p_val = U1024::zero();
    p_val.0[0] = 43;
    let params = MontgomeryParams::new(p_val, U1024::zero());

    let a = FieldElement::new(U1024::from_u64(23), &params);
    let b = FieldElement::new(U1024::from_u64(42), &params);
    let curve = WeierstrassCurve::new(a, b, &params);

    let g = SWPoint::identity(&curve);
    let g2 = g.double();

    assert!(g2.is_identity(), "Double infinity must be infinity");
}
