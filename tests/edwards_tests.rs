//! Tests for Edwards curve operations.

use curvelib::{
    instances::bls6_6::Bls6_6BaseField,
    instances::tiny_jubjub::{self, TinyJubjubBaseField},
    models::{EdwardsPoint, WeierstrassCurve},
    traits::{Curve, ProjectivePoint},
};
use mathlib::{FieldElement, fp};

#[test]
fn test_tiny_jubjub_addition() {
    let curve = tiny_jubjub::get_curve();

    // Identity point: (0, 1)
    let p1 = EdwardsPoint::new_affine(
        FieldElement::<TinyJubjubBaseField>::zero(),
        FieldElement::<TinyJubjubBaseField>::one(),
        curve.clone(),
    );
    assert!(p1.is_identity());

    // Another point on the curve
    let p2 = EdwardsPoint::new_affine(
        fp!(1u64, TinyJubjubBaseField),
        fp!(2u64, TinyJubjubBaseField),
        curve.clone(),
    );

    let sum1 = p1.add(&p2);
    assert_eq!(sum1, p2);

    let double_p2 = p2.double();
    assert!(!double_p2.is_identity());

    let add_p2 = p2.add(&p2);
    assert_eq!(
        double_p2, add_p2,
        "Double and Add(P,P) must be equal on Edwards"
    );

    let (x, y) = double_p2.to_affine();
    println!("2 * (1,2) = ({:?}, {:?})", x.to_u1024(), y.to_u1024());
}

#[test]
fn test_tiny_curve_operations() {
    // Use BLS6_6 for Weierstrass curve test
    let a = fp!(23u64, Bls6_6BaseField);
    let b = fp!(42u64, Bls6_6BaseField);
    let gx = fp!(1u64, Bls6_6BaseField);
    let gy = fp!(1u64, Bls6_6BaseField);

    let curve = WeierstrassCurve::<Bls6_6BaseField>::new(a, b, gx, gy);

    let g = curve.identity();
    let g2 = g.double();

    assert!(g2.is_identity(), "Double infinity must be infinity");
}

#[test]
fn test_edwards_generator_on_curve() {
    let curve = tiny_jubjub::get_curve();
    let g = curve.generator();
    let (x, y) = g.to_affine();
    assert!(curve.is_on_curve(&x, &y), "Generator must be on curve");
}

#[test]
fn test_edwards_identity() {
    let curve = tiny_jubjub::get_curve();
    let id = curve.identity();
    assert!(id.is_identity());

    // Identity has coords (0, 1)
    let (x, y) = id.to_affine();
    assert!(x.is_zero());
    assert_eq!(y, FieldElement::<TinyJubjubBaseField>::one());
}
