use curvelib::{
    algebra::fields::Fp,
    instances::tiny_jubjub,
    models::TePoint,
    traits::{Curve, ProjectivePoint},
};
use mathlib::{BigInt, U1024};

#[test]
fn test_scalar_multiplication() {
    let curve = tiny_jubjub::get_curve();
    let params = curve.params;

    // Valid point on tiny_jubjub: (3, 0)
    // 3*3^2 + 0^2 = 27 = 1 mod 13
    // 1 + 8*3^2*0^2 = 1 mod 13
    let p = TePoint::new_affine(
        Fp::new(U1024::from_u64(3), params),
        Fp::new(U1024::from_u64(0), params),
        &curve,
    );

    assert!(curve.is_on_curve(&p.x, &p.y), "Point must be on the curve");

    let scalar_2 = U1024::from_u64(2);
    let mul_2 = p.mul(&scalar_2);
    let add_2 = p.add(&p);

    assert_eq!(mul_2, add_2, "Scalar mul by 2 failed");

    let scalar_0 = U1024::zero();
    let mul_0 = p.mul(&scalar_0);
    assert!(mul_0.is_identity(), "Scalar mul by 0 must be identity");

    let scalar_1 = U1024::from_u64(1);
    let mul_1 = p.mul(&scalar_1);
    assert_eq!(mul_1, p, "Scalar mul by 1 must be self");

    let scalar_5 = U1024::from_u64(5);
    let mul_5 = p.mul(&scalar_5);

    let p2 = p.double();
    let p4 = p2.double();
    let p5_add = p4.add(&p);

    assert_eq!(mul_5, p5_add, "Scalar mul by 5 failed");

    let (x, y) = mul_5.to_affine();
    println!("5 * P = ({:?}, {:?})", x.to_u1024(), y.to_u1024());
}
