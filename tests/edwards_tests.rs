use curvelib::curves::tiny_jubjub;
use curvelib::models::twisted_edwards::TePoint;
use curvelib::traits::ProjectivePoint;
use mathlib::field::element::FieldElement;
use mathlib::{BigInt, U1024};

#[test]
fn test_tiny_jubjub_addition() {
    let curve = tiny_jubjub::get_curve();
    let params = curve.params;

    let p1 = TePoint::new_affine(
        FieldElement::new(U1024::from_u64(0), params),
        FieldElement::new(U1024::from_u64(1), params),
        &curve,
    );
    assert!(p1.is_identity());

    let p2 = TePoint::new_affine(
        FieldElement::new(U1024::from_u64(1), params),
        FieldElement::new(U1024::from_u64(2), params),
        &curve,
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
