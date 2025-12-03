use std::sync::OnceLock;

use mathlib::field::element::FieldElement;
use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{BigInt, U1024};

use crate::models::twisted_edwards::EdwardsCurve;

static PARAMS: OnceLock<MontgomeryParams> = OnceLock::new();

pub fn get_tiny_params() -> &'static MontgomeryParams {
    PARAMS.get_or_init(|| {
        let p = U1024::from_u64(13);
        MontgomeryParams::new(p, U1024::zero())
    })
}

pub fn get_curve() -> EdwardsCurve<'static> {
    let params = get_tiny_params();
    let a = FieldElement::new(U1024::from_u64(3), params);
    let d = FieldElement::new(U1024::from_u64(8), params);
    EdwardsCurve::new(a, d, params)
}
