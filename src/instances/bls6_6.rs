use std::sync::OnceLock;

use mathlib::field::element::FieldElement;
use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{BigInt, U1024};

use crate::algebra::fields::fp::Fp;
use crate::models::short_weierstrass::WeierstrassCurve;

static PARAMS: OnceLock<MontgomeryParams> = OnceLock::new();
static SCALAR_PARAMS: OnceLock<MontgomeryParams> = OnceLock::new();

/// Returns the field parameters for the bls6_6 curve (p = 43).
pub fn get_params() -> &'static MontgomeryParams {
    PARAMS.get_or_init(|| {
        let p = U1024::from_u64(43);
        MontgomeryParams::new(p, U1024::zero())
    })
}

/// Returns the scalar field parameters for the bls6_6 curve (order = 39).
pub fn get_scalar_params() -> &'static MontgomeryParams {
    SCALAR_PARAMS.get_or_init(|| {
        let order = U1024::from_u64(39);
        MontgomeryParams::new(order, U1024::zero())
    })
}

/// Generator point coordinates for bls6_6: (13, 15)
pub fn get_generator_coords() -> (U1024, U1024) {
    (U1024::from_u64(13), U1024::from_u64(15))
}

/// Beta value for the quadratic extension field.
pub fn get_beta() -> U1024 {
    U1024::from_u64(42)
}

/// Returns the bls6_6 Weierstrass curve: y^2 = x^3 + 6 over F_43.
pub fn get_curve() -> WeierstrassCurve<'static> {
    let params = get_params();
    let scalar_params = get_scalar_params();

    // a = 0, b = 6
    let a = Fp::from(FieldElement::zero(params));
    let b = Fp::new(U1024::from_u64(6), params);

    let (gen_x_val, gen_y_val) = get_generator_coords();
    let gen_x = Fp::new(gen_x_val, params);
    let gen_y = Fp::new(gen_y_val, params);

    WeierstrassCurve::new(a, b, params, scalar_params, gen_x, gen_y)
}
