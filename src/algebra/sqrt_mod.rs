//! Modular square root computation using Tonelli-Shanks algorithm.
//!
//! This module provides the `sqrt_mod` function for computing square roots
//! in prime fields.

use mathlib::{BigInt, FieldConfig, FieldElement, U1024};

/// Compute a modular square root of `n` in its prime field using the Tonelli–Shanks algorithm.
///
/// Attempts to find `r` such that `r * r == n` in the field of `n`. Handles the zero case,
/// uses Euler's criterion to test quadratic residuosity, applies the P ≡ 3 (mod 4) shortcut
/// when applicable, and otherwise runs the Tonelli–Shanks iteration to produce a root.
///
/// # Type Parameters
///
/// * `C` - The field configuration
///
/// # Returns
///
/// `Some(root)` containing a field element `r` such that `r * r == n` when a root exists, `None` otherwise.
pub fn sqrt_mod<C: FieldConfig>(n: &FieldElement<C>) -> Option<FieldElement<C>> {
    let zero = FieldElement::<C>::zero();
    let one = FieldElement::<C>::one();

    // Edge case: √0 = 0
    if *n == zero {
        return Some(zero);
    }

    // Step 1: Test quadratic residuosity using Euler's criterion
    // For prime p, n is a quadratic residue iff n^((p-1)/2) ≡ 1 (mod p)
    let p = C::MODULUS;
    let u_one = U1024::from_u64(1);
    let u_two = U1024::from_u64(2);
    let p_minus_1 = p - u_one;
    let (legendre_exp, _) = p_minus_1.div_rem(&u_two);

    // Check if n^((p-1)/2) == 1; if not, n has no square root
    if n.pow(legendre_exp) != one {
        return None;
    }

    // Step 2: Factor out powers of 2 from (p-1)
    // Write p - 1 = Q * 2^S where Q is odd
    let mut s = 0u32;
    let mut q = p_minus_1;

    // Extract all factors of 2 from q
    loop {
        let (div, rem) = q.div_rem(&u_two);
        if rem != U1024::zero() {
            break;
        }
        q = div;
        s += 1;
    }

    // Step 3: Special case for p ≡ 3 (mod 4)
    if s == 1 {
        let p_plus_1 = p + u_one;
        let u_four = U1024::from_u64(4);
        let (exp, _) = p_plus_1.div_rem(&u_four);
        return Some(n.pow(exp));
    }

    // Step 4: Find a quadratic non-residue z
    let mut z = u_two;
    let neg_one = zero - one;
    let mut z_elem = FieldElement::<C>::new(z);

    loop {
        if z_elem.pow(legendre_exp) == neg_one {
            break;
        }
        z = z + u_one;
        z_elem = FieldElement::<C>::new(z);
    }

    // Step 5: Initialize Tonelli-Shanks variables
    let mut c = z_elem.pow(q);
    let (exp_r, _) = (q + u_one).div_rem(&u_two);

    let mut r = n.pow(exp_r);
    let mut t = n.pow(q);
    let mut m = s;

    // Step 6: Main Tonelli-Shanks iteration
    loop {
        if t == one {
            return Some(r);
        }

        // Find the least i such that t^(2^i) = 1
        let mut i = 0u32;
        let mut temp = t;
        while temp != one && i < m {
            temp = temp * temp;
            i += 1;
        }

        if i == m {
            return None;
        }

        // Compute b = c^(2^(m-i-1))
        let mut b = c;
        for _ in 0..(m - i - 1) {
            b = b * b;
        }

        // Update variables for next iteration
        m = i;
        c = b * b;
        t = t * c;
        r = r * b;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instances::bls6_6::Bls6_6BaseField;
    use mathlib::fp;

    #[test]
    fn test_sqrt_of_zero() {
        let zero = FieldElement::<Bls6_6BaseField>::zero();
        let result = sqrt_mod(&zero);
        assert_eq!(result, Some(zero));
    }

    #[test]
    fn test_sqrt_quadratic_residue() {
        // In F_43, 4 is a quadratic residue (2² = 4)
        let four = fp!(4u64, Bls6_6BaseField);
        let result = sqrt_mod(&four);
        assert!(result.is_some());
        let r = result.unwrap();
        assert_eq!(r * r, four);
    }
}
