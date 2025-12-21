//! Pairing operations for elliptic curves.
//!
//! This module provides the Tate pairing implementation for pairing-friendly curves.

pub mod final_exp;
pub mod miller;

use mathlib::{FieldConfig, U1024};

use crate::algebra::fields::Fp6;
use crate::models::{TwistPoint, WeierstrassPoint};

/// Type aliases for clarity.
pub type G1Affine<C> = WeierstrassPoint<C>;
pub type G2Projective<C> = TwistPoint<C>;

/// Computes the Tate pairing of a G1 point and a G2 point.
///
/// The result is the pairing value in the Fp6 target field computed by:
/// 1. Running the Miller loop
/// 2. Applying the final exponentiation
///
/// # Type Parameters
///
/// * `C` - The base field configuration
///
/// # Parameters
///
/// * `p` - A G1 point
/// * `q` - A G2 point
/// * `r_order` - The curve subgroup order
/// * `final_exp_val` - The final exponentiation exponent
///
/// # Returns
///
/// A Fp6 element representing the pairing value.
pub fn tate_pairing<C: FieldConfig>(
    p: &G1Affine<C>,
    q: &G2Projective<C>,
    r_order: U1024,
    final_exp_val: U1024,
) -> Fp6<C> {
    // 1. Miller Loop
    let f = miller::miller_loop(p, q, r_order);

    // 2. Final Exponentiation
    final_exp::final_exponentiation(&f, &final_exp_val)
}

#[cfg(test)]
mod tests {
    // Tests will be added once all dependencies compile
}
