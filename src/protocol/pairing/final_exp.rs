use mathlib::U1024;

use crate::{algebra::fields::Fp6, traits::Field};

/// Compute exponentiation of an Fp6 element by an integer exponent using square-and-multiply.
///
/// # Parameters
/// - `f`: base element in Fp6.
/// - `exponent`: exponent as a `U1024` integer.
///
/// # Returns
/// The value `f` raised to `exponent` in the Fp6 field.
///
/// # Examples
///
/// ```rust
/// use curvelib::algebra::fields::Fp6;
/// use curvelib::instances::bls6_6::get_params;
/// use curvelib::protocol::pairing::final_exp::final_exponentiation;
/// use curvelib::traits::Field;
/// use mathlib::{U1024, BigInt};
///
/// let params = get_params();
///
/// // pick a simple, concrete base element
/// let f = Fp6::one(params);
/// let exp = U1024::from_u64(3);
///
/// let r = final_exponentiation(&f, &exp);
/// assert_eq!(r, f * f * f);
/// ```
pub fn final_exponentiation<'a>(f: &Fp6<'a>, exponent: &U1024) -> Fp6<'a> {
    // Perform f ^ exponent using Square-and-Multiply
    // Note: f is in Fp6.

    let params = f.c0.c0.params;
    // xi generation removed as it's not used in standard mul

    let mut res = Fp6::one(params);
    let mut base = *f;

    let num_bits = exponent.bits();

    for i in 0..num_bits {
        if exponent.bit(i) {
            res = res * base;
        }
        base = base * base;
    }

    res
}
