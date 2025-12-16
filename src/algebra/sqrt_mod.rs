use mathlib::{BigInt, FieldElement, U1024, fp, u1024};

/// Compute a modular square root of `n` in its prime field using the Tonelli–Shanks algorithm.
///
/// Attempts to find `r` such that `r * r == n` in the field of `n`. Handles the zero case,
/// uses Euler's criterion to test quadratic residuosity, applies the P ≡ 3 (mod 4) shortcut
/// when applicable, and otherwise runs the Tonelli–Shanks iteration to produce a root.
///
/// # Returns
///
/// `Some(root)` containing a field element `r` such that `r * r == n` when a root exists, `None` otherwise.
///
/// # Examples
///
/// ```rust
/// use curvelib::algebra::sqrt_mod::sqrt_mod;
/// use mathlib::{fp, u1024, BigInt, mont};
///
/// // Work in the tiny prime field F_13.
/// let params = mont!(u1024!(13u64), u1024!(0u64));
///
/// // 10 is a quadratic residue mod 13 (since 6^2 = 36 ≡ 10 mod 13).
/// let n = fp!(u1024!(10u64), &params);
/// let r = sqrt_mod(&n).expect("10 should have a square root in F_13");
/// assert_eq!(r * r, n);
///
/// // 2 is a non-residue mod 13, so no square root exists.
/// let n2 = fp!(u1024!(2u64), &params);
/// assert!(sqrt_mod(&n2).is_none());
/// ```
pub fn sqrt_mod<'a>(n: &FieldElement<'a>) -> Option<FieldElement<'a>> {
    let params = n.params;
    // Create field constants we'll need throughout
    let zero = FieldElement::zero(params);
    let one = fp!(u1024!(1u64), params);

    // Edge case: √0 = 0
    if *n == zero {
        return Some(zero);
    }

    // Step 1: Test quadratic residuosity using Euler's criterion
    // For prime p, n is a quadratic residue iff n^((p-1)/2) ≡ 1 (mod p)
    let p = params.modulus;
    let u_one = u1024!(1u64);
    let u_two = u1024!(2u64);
    let p_minus_1 = p - u_one; // p - 1
    let (legendre_exp, _) = p_minus_1.div_rem(&u_two); // (p-1)/2

    // Check if n^((p-1)/2) == 1; if not, n has no square root
    if n.pow(legendre_exp) != one {
        return None; // n is not a quadratic residue
    }

    // Step 2: Factor out powers of 2 from (p-1)
    // Write p - 1 = Q * 2^S where Q is odd
    let mut s = 0u32; // Exponent S (number of times 2 divides p-1)
    let mut q = p_minus_1; // Will become the odd part Q

    // Extract all factors of 2 from q
    loop {
        let (div, rem) = q.div_rem(&u_two);
        if rem != U1024::zero() {
            break; // q is now odd
        }
        q = div; // Continue dividing by 2
        s += 1; // Count the factor of 2
    }

    // Step 3: Special case for p ≡ 3 (mod 4)
    // When S=1, we have p = Q*2 + 1, so p ≡ 3 (mod 4)
    // In this case, r = n^((p+1)/4) is a square root
    if s == 1 {
        let p_plus_1 = p + u_one;
        let u_four = u1024!(4u64);
        let (exp, _) = p_plus_1.div_rem(&u_four); // (p+1)/4
        return Some(n.pow(exp));
    }

    // Step 4: Find a quadratic non-residue z
    // We need some element z where z^((p-1)/2) ≡ -1 (mod p)
    let mut z = u_two; // Start searching from 2
    let neg_one = zero - one; // -1 in the field
    let mut z_elem = fp!(z, params);

    // Search for a non-residue by testing successive integers
    loop {
        if z_elem.pow(legendre_exp) == neg_one {
            break; // Found a non-residue
        }
        z = z + u_one; // Try next integer
        z_elem = fp!(z, params);
    }

    // Step 5: Initialize Tonelli-Shanks variables
    // c is our "generator" raised to an odd power
    let mut c = z_elem.pow(q); // c = z^Q
    let (exp_r, _) = (q + u_one).div_rem(&u_two); // (Q+1)/2

    let mut r = n.pow(exp_r); // r will converge to the square root
    let mut t = n.pow(q); // t tracks n^Q in the iteration
    let mut m = s; // m is the current "order" exponent

    // Step 6: Main Tonelli-Shanks iteration
    // Loop invariant: r^2 * t = n (mod p)
    loop {
        // If t = 1, then r^2 = n, so r is our answer
        if t == one {
            return Some(r);
        }

        // Find the least i such that t^(2^i) = 1
        let mut i = 0u32;
        let mut temp = t;
        while temp != one && i < m {
            temp = temp * temp; // Square repeatedly
            i += 1;
        }

        // Sanity check: i should be less than m
        if i == m {
            return None; // Should not happen if n is a quadratic residue
        }

        // Compute b = c^(2^(m-i-1))
        // This adjusts our approximation r to get closer to the true root
        let mut b = c;
        for _ in 0..(m - i - 1) {
            b = b * b; // Square (m-i-1) times
        }

        // Update variables for next iteration
        m = i; // New "order" is i
        c = b * b; // New c is b^2
        t = t * c; // Update t (gets us closer to 1)
        r = r * b; // Update r (refine our square root approximation)
    }
}
