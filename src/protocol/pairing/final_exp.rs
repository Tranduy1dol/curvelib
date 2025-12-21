//! Final exponentiation for pairing computation.
//!
//! This module provides the final exponentiation step of the Tate pairing.

use mathlib::{FieldConfig, U1024};

use crate::algebra::fields::Fp6;

/// Compute exponentiation of an Fp6 element by an integer exponent.
///
/// Uses square-and-multiply algorithm.
///
/// # Type Parameters
///
/// * `C` - The base field configuration
///
/// # Parameters
///
/// * `f` - Base element in Fp6
/// * `exponent` - Exponent as a `U1024` integer
///
/// # Returns
///
/// The value `f` raised to `exponent` in the Fp6 field.
pub fn final_exponentiation<C: FieldConfig>(f: &Fp6<C>, exponent: &U1024) -> Fp6<C> {
    let mut res = Fp6::<C>::one();
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instances::bls6_6::Bls6_6BaseField;
    use mathlib::u1024;

    #[test]
    fn test_final_exp_one() {
        let one = Fp6::<Bls6_6BaseField>::one();
        let exp = u1024!(5);
        let result = final_exponentiation(&one, &exp);
        assert_eq!(result, one);
    }
}
