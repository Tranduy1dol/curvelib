use mathlib::{BigInt, FieldElement, U1024};

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
/// use mathlib::field::montgomery::MontgomeryParams;
/// use mathlib::{FieldElement, U1024, BigInt};
///
/// // Work in the tiny prime field F_13.
/// let params = MontgomeryParams::new(U1024::from_u64(13), U1024::zero());
///
/// // 10 is a quadratic residue mod 13 (since 6^2 = 36 ≡ 10 mod 13).
/// let n = FieldElement::new(U1024::from_u64(10), &params);
/// let r = sqrt_mod(&n).expect("10 should have a square root in F_13");
/// assert_eq!(r * r, n);
///
/// // 2 is a non-residue mod 13, so no square root exists.
/// let n2 = FieldElement::new(U1024::from_u64(2), &params);
/// assert!(sqrt_mod(&n2).is_none());
/// ```
pub fn sqrt_mod<'a>(n: &FieldElement<'a>) -> Option<FieldElement<'a>> {
    let params = n.params;
    let zero = FieldElement::zero(params);
    let one = FieldElement::new(U1024::from_u64(1), params);

    if *n == zero {
        return Some(zero);
    }

    // Euler's criterion: n^((P-1)/2) == 1
    let p = params.modulus;
    let u_one = U1024::from_u64(1);
    let u_two = U1024::from_u64(2);
    let p_minus_1 = p - u_one;
    let (legendre_exp, _) = p_minus_1.div_rem(&u_two);

    if n.pow(legendre_exp) != one {
        return None;
    }

    // P - 1 = Q * 2^S
    let mut s = 0u32;
    let mut q = p_minus_1;

    // Check q % 2 == 0
    loop {
        let (div, rem) = q.div_rem(&u_two);
        if rem != U1024::zero() {
            break;
        }
        q = div;
        s += 1;
    }

    if s == 1 {
        // Case P = 3 mod 4: r = n^((P+1)/4)
        let p_plus_1 = p + u_one;
        let u_four = U1024::from_u64(4);
        let (exp, _) = p_plus_1.div_rem(&u_four);
        return Some(n.pow(exp));
    }

    // Find z such that z is a quadratic non-residue
    let mut z = u_two;
    let neg_one = zero - one;
    let mut z_elem = FieldElement::new(z, params);

    loop {
        if z_elem.pow(legendre_exp) == neg_one {
            break;
        }
        z = z + u_one;
        z_elem = FieldElement::new(z, params);
    }

    let mut c = z_elem.pow(q);
    let (exp_r, _) = (q + u_one).div_rem(&u_two);

    let mut r = n.pow(exp_r);
    let mut t = n.pow(q);
    let mut m = s;

    loop {
        if t == one {
            return Some(r);
        }

        let mut i = 0u32;
        let mut temp = t;
        while temp != one && i < m {
            temp = temp * temp;
            i += 1;
        }

        if i == m {
            return None; // Should not happen
        }

        let mut b = c;
        for _ in 0..(m - i - 1) {
            b = b * b;
        }

        m = i;
        c = b * b;
        t = t * c;
        r = r * b;
    }
}
