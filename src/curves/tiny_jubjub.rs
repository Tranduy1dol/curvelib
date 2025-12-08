use std::sync::OnceLock;

use mathlib::field::element::FieldElement;
use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{BigInt, U1024};

use crate::models::twisted_edwards::EdwardsCurve;

static PARAMS: OnceLock<MontgomeryParams> = OnceLock::new();
static SCALAR_PARAMS: OnceLock<MontgomeryParams> = OnceLock::new();

/// Returns the field parameters for the tiny jubjub curve (p = 13).
pub fn get_tiny_params() -> &'static MontgomeryParams {
    PARAMS.get_or_init(|| {
        let p = U1024::from_u64(13);
        MontgomeryParams::new(p, U1024::zero())
    })
}

/// Returns the scalar field parameters for the tiny jubjub curve.
pub fn get_scalar_params() -> &'static MontgomeryParams {
    SCALAR_PARAMS.get_or_init(|| {
        let order = U1024::from_u64(5);
        MontgomeryParams::new(order, U1024::zero())
    })
}

/// Generator point coordinates for tiny jubjub: (6, 9)
/// This point has order 5.
pub fn get_generator_coords() -> (U1024, U1024) {
    (U1024::from_u64(6), U1024::from_u64(9))
}

/// Returns the tiny jubjub Edwards curve with a = 3, d = 8 over F_13.
/// Generator point is (1, 2) which is on the curve.
pub fn get_curve() -> EdwardsCurve<'static> {
    let params = get_tiny_params();
    let scalar_params = get_scalar_params();
    let a = FieldElement::new(U1024::from_u64(3), params);
    let d = FieldElement::new(U1024::from_u64(8), params);

    let (gen_x_val, gen_y_val) = get_generator_coords();
    EdwardsCurve::new(a, d, params, scalar_params, gen_x_val, gen_y_val)
}
