use std::sync::OnceLock;

use mathlib::{FieldElement, MontgomeryParams, U1024, fp, mont, u1024};

use crate::models::WeierstrassCurve;

static PARAMS: OnceLock<MontgomeryParams> = OnceLock::new();
static SCALAR_PARAMS: OnceLock<MontgomeryParams> = OnceLock::new();

/// Field parameters for the bls6_6 base field (prime p = 43).
///
/// The returned value is a static reference to the Montgomery parameters used for
/// arithmetic in the curve's base field.
///
/// # Examples
///
/// ```rust
/// use curvelib::instances::bls6_6::get_params;
/// use mathlib::u1024;
///
/// let params = get_params();
/// assert_eq!(params.modulus, u1024!(43));
/// ```
pub fn get_params() -> &'static MontgomeryParams {
    PARAMS.get_or_init(|| mont!(u1024!(43), u1024!(0)))
}

/// Provides the Montgomery parameters for the scalar field of the bls6_6 curve (modulus = 39).
///
/// # Examples
///
/// ```rust
/// use curvelib::instances::bls6_6::get_scalar_params;
/// use mathlib::u1024;
///
/// let params = get_scalar_params();
/// assert_eq!(params.modulus, u1024!(39));
/// ```
pub fn get_scalar_params() -> &'static MontgomeryParams {
    SCALAR_PARAMS.get_or_init(|| mont!(u1024!(39), u1024!(0)))
}

/// Generator point coordinates for the bls6_6 curve.
///
/// The generator point is returned as an (x, y) pair in the base field, represented as `U1024` values.
///
/// # Examples
///
/// ```rust
/// use curvelib::instances::bls6_6::get_generator_coords;
/// use mathlib::u1024;
///
/// let (x, y) = get_generator_coords();
/// assert_eq!(x, u1024!(13));
/// assert_eq!(y, u1024!(15));
/// ```
pub fn get_generator_coords() -> (U1024, U1024) {
    (u1024!(13), u1024!(15))
}

/// Beta constant for the quadratic extension field.
///
/// Returns the beta value as a `U1024`.
///
/// # Examples
///
/// ```rust
/// use curvelib::instances::bls6_6::get_beta;
/// use mathlib::u1024;
///
/// let beta = get_beta();
/// assert_eq!(beta, u1024!(42));
/// ```
pub fn get_beta() -> U1024 {
    u1024!(42)
}

/// Constructs the BLS6-6 Weierstrass curve over F_43 with equation y^2 = x^3 + 6.
///
/// The curve is created with curve parameters a = 0 and b = 6, uses the base field
/// modulus p = 43 and the scalar field order = 39, and is returned with its generator point
/// (x, y) = (13, 15).
///
/// # Examples
///
/// ```rust
/// use curvelib::instances::bls6_6::{get_curve, get_generator_coords, get_params, get_scalar_params};
/// use mathlib::u1024;
///
/// let curve = get_curve();
///
/// // sanity-check the instance parameters
/// assert_eq!(get_params().modulus, u1024!(43));
/// assert_eq!(get_scalar_params().modulus, u1024!(39));
/// assert_eq!(get_generator_coords(), (u1024!(13), u1024!(15)));
///
/// // and the curve wires the same params in
/// assert_eq!(curve.params.modulus, u1024!(43));
/// assert_eq!(curve.scalar_params.modulus, u1024!(39));
/// ```
pub fn get_curve() -> WeierstrassCurve<'static> {
    let params = get_params();
    let scalar_params = get_scalar_params();

    // a = 0, b = 6
    let a = FieldElement::zero(params);
    let b = fp!(u1024!(6), params);

    let (gen_x_val, gen_y_val) = get_generator_coords();
    let gen_x = fp!(gen_x_val, params);
    let gen_y = fp!(gen_y_val, params);

    WeierstrassCurve::new(a, b, params, scalar_params, gen_x, gen_y)
}
