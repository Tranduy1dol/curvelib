use std::sync::OnceLock;

use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{BigInt, U1024};

use crate::algebra::fields::Fp;
use crate::models::EdwardsCurve;

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

/// Tiny Jubjub Twisted Edwards curve over F_13 with curve parameters a = 3 and d = 8.
///
/// The curve is constructed using the crate's tiny field and scalar parameters. The
/// generator point used is (6, 9) (as U1024 values) and has order 5.
///
/// # Examples
///
/// ```rust
/// use curvelib::instances::tiny_jubjub::{get_curve, get_generator_coords, get_scalar_params, get_tiny_params};
/// use mathlib::{U1024, BigInt};
///
/// let curve = get_curve();
///
/// // sanity-check the parameters used by the instance
/// assert_eq!(get_tiny_params().modulus, U1024::from_u64(13));
/// assert_eq!(get_scalar_params().modulus, U1024::from_u64(5));
/// assert_eq!(get_generator_coords(), (U1024::from_u64(6), U1024::from_u64(9)));
///
/// // and the curve wires the same params in
/// assert_eq!(curve.params.modulus, U1024::from_u64(13));
/// assert_eq!(curve.scalar_params.modulus, U1024::from_u64(5));
/// ```
pub fn get_curve() -> EdwardsCurve<'static> {
    let params = get_tiny_params();
    let scalar_params = get_scalar_params();
    let a = Fp::new(U1024::from_u64(3), params);
    let d = Fp::new(U1024::from_u64(8), params);

    let (gen_x_val, gen_y_val) = get_generator_coords();
    EdwardsCurve::new(a, d, params, scalar_params, gen_x_val, gen_y_val)
}
