use curvelib::{
    algebra::fields::Fp,
    instances::tiny_jubjub,
    models::{TePoint, WeierstrassCurve},
    traits::{Curve, ProjectivePoint},
};
use mathlib::{BigInt, U1024, field::montgomery::MontgomeryParams};

#[test]
fn test_tiny_jubjub_addition() {
    let curve = tiny_jubjub::get_curve();
    let params = curve.params;

    let p1 = TePoint::new_affine(
        Fp::new(U1024::from_u64(0), params),
        Fp::new(U1024::from_u64(1), params),
        &curve,
    );
    assert!(p1.is_identity());

    let p2 = TePoint::new_affine(
        Fp::new(U1024::from_u64(1), params),
        Fp::new(U1024::from_u64(2), params),
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

/// Checks that doubling the identity point on a Weierstrass curve yields the identity.
///
/// # Examples
///
/// ```
/// use mathlib::{U1024, Fp, field::montgomery::MontgomeryParams};
/// use short_weierstrass::WeierstrassCurve;
///
/// let mut p_val = U1024::zero();
/// p_val.0[0] = 43;
/// let params = MontgomeryParams::new(p_val, U1024::zero());
///
/// let a = Fp::new(U1024::from_u64(23), &params);
/// let b = Fp::new(U1024::from_u64(42), &params);
/// let curve = WeierstrassCurve::new(
///     a,
///     b,
///     &params,
///     &params,
///     Fp::new(U1024::from_u64(1), &params),
///     Fp::new(U1024::from_u64(1), &params),
/// );
///
/// let g = curve.identity();
/// let g2 = g.double();
///
/// assert!(g2.is_identity());
/// ```
#[test]
fn test_tiny_curve_operations() {
    let mut p_val = U1024::zero();
    p_val.0[0] = 43;
    let params = MontgomeryParams::new(p_val, U1024::zero());

    let a = Fp::new(U1024::from_u64(23), &params);
    let b = Fp::new(U1024::from_u64(42), &params);
    let curve = WeierstrassCurve::new(
        a,
        b,
        &params,
        &params,
        Fp::new(U1024::from_u64(1), &params),
        Fp::new(U1024::from_u64(1), &params),
    );

    let g = curve.identity();
    let g2 = g.double();

    assert!(g2.is_identity(), "Double infinity must be infinity");
}
