//! KZG polynomial commitment scheme.
//!
//! This module provides the KZG (Kate-Zaverucha-Goldberg) polynomial commitment
//! scheme for use in zero-knowledge proofs.
//!
//! Uses `mathlib::Polynomial` for polynomial operations.

use std::marker::PhantomData;
use std::ops::Neg;

use mathlib::{FieldConfig, FieldElement, Polynomial, U1024};

use crate::models::{SexticTwist, TwistPoint, WeierstrassCurve, WeierstrassPoint};
use crate::protocol::pairing::tate_pairing;
use crate::traits::{Curve, ProjectivePoint};

/// Type aliases for clarity.
pub type G1Point<C> = WeierstrassPoint<C>;
pub type G2Point<C> = TwistPoint<C>;

/// KZG Parameters (Structured Reference String).
pub struct KzgParams<C: FieldConfig> {
    /// Powers of tau in G1: [G, τG, τ²G, ...]
    pub powers_of_tau: Vec<G1Point<C>>,
    /// G2 generator
    pub g2: G2Point<C>,
    /// τ·G2
    pub tau_g2: G2Point<C>,
}

/// KZG polynomial commitment scheme.
pub struct Kzg<C: FieldConfig> {
    _marker: PhantomData<C>,
}

impl<C: FieldConfig> Kzg<C> {
    /// Trusted Setup: Generate SRS from secret tau.
    ///
    /// WARNING: In production, tau must come from a secure MPC ceremony!
    pub fn setup(
        tau: FieldElement<C>,
        max_degree: usize,
        g1_curve: &WeierstrassCurve<C>,
        g2_curve: &SexticTwist<C>,
    ) -> KzgParams<C> {
        let g1 = g1_curve.generator();
        let g2 = g2_curve.generator();

        let mut powers_of_tau = Vec::with_capacity(max_degree + 1);
        let mut curr_tau = FieldElement::<C>::one();

        for _ in 0..=max_degree {
            let scalar = curr_tau.to_u1024();
            let point = g1.mul(&scalar);
            powers_of_tau.push(point);
            curr_tau = curr_tau * tau;
        }

        let tau_g2 = g2.mul(&tau.to_u1024());

        KzgParams {
            powers_of_tau,
            g2: g2.clone(),
            tau_g2,
        }
    }

    /// Commit to a polynomial.
    ///
    /// Computes C = Σ cᵢ · \[τⁱ\]G₁
    pub fn commit(params: &KzgParams<C>, poly: &Polynomial<C>) -> G1Point<C> {
        let coeffs = poly.to_vec();
        let mut c = params.powers_of_tau[0].curve.identity();

        for (i, coeff) in coeffs.iter().enumerate() {
            if i < params.powers_of_tau.len() {
                let scalar = coeff.to_u1024();
                let term = params.powers_of_tau[i].mul(&scalar);
                c = c.add(&term);
            }
        }
        c
    }

    /// Open: Create proof for P(z) = y.
    ///
    /// Returns the proof π and the evaluation y.
    pub fn open(
        params: &KzgParams<C>,
        poly: &Polynomial<C>,
        z: &FieldElement<C>,
    ) -> (G1Point<C>, FieldElement<C>) {
        // 1. Evaluate y = P(z)
        let y = poly.evaluate(z);

        // 2. Compute quotient Q(x) = (P(x) - y) / (x - z)
        // P(x) - y
        let p_minus_y = poly.clone() - Polynomial::constant(y);

        // (x - z)
        let divisor = Polynomial::new(vec![-*z, FieldElement::<C>::one()]);

        // Divide using mathlib's polynomial division
        let (quotient, _remainder) = p_minus_y.divide_with_remainder(&divisor);

        // 3. Commit to Q(x) -> Proof π
        let pi = Self::commit(params, &quotient);

        (pi, y)
    }

    /// Verify an opening proof using pairing.
    ///
    /// Checks: e(π, \[τ\]G₂ - \[`z`\]G₂) = e(C - \[`y`\]G₁, G₂)
    pub fn verify(
        params: &KzgParams<C>,
        commitment: &G1Point<C>,
        z: &FieldElement<C>,
        y: &FieldElement<C>,
        proof_pi: &G1Point<C>,
        r_order: U1024,
        final_exp: U1024,
    ) -> bool {
        let z_scalar = z.to_u1024();
        let y_scalar = y.to_u1024();

        // LHS = e(π, τ·G₂ - z·G₂)
        let g2_z = params.g2.mul(&z_scalar);
        let g2_term = params.tau_g2.clone().add(&g2_z.neg());

        let lhs = tate_pairing(proof_pi, &g2_term, r_order, final_exp);

        // RHS = e(C - y·G₁, G₂)
        let g1 = &params.powers_of_tau[0];
        let gy = g1.mul(&y_scalar);
        let g1_term = commitment.add(&gy.neg());

        let rhs = tate_pairing(&g1_term, &params.g2, r_order, final_exp);

        lhs == rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instances::bls6_6::{Bls6_6BaseField, get_g1_curve, get_g2_curve};

    #[test]
    fn test_polynomial_operations() {
        // Create a simple polynomial: P(x) = 2 + 3x + x²
        let coeffs = vec![
            FieldElement::<Bls6_6BaseField>::new(U1024::from_u64(2)),
            FieldElement::<Bls6_6BaseField>::new(U1024::from_u64(3)),
            FieldElement::<Bls6_6BaseField>::new(U1024::from_u64(1)),
        ];
        let poly = Polynomial::new(coeffs);

        // Evaluate at z = 2: P(2) = 2 + 6 + 4 = 12
        let z = FieldElement::<Bls6_6BaseField>::new(U1024::from_u64(2));
        let y = poly.evaluate(&z);
        assert_eq!(y.to_u1024().0[0], 12);
    }

    #[test]
    fn test_kzg_setup() {
        let g1_curve = get_g1_curve();
        let g2_curve = get_g2_curve();
        let tau = FieldElement::<Bls6_6BaseField>::new(U1024::from_u64(5));

        let params = Kzg::<Bls6_6BaseField>::setup(tau, 3, &g1_curve, &g2_curve);

        assert_eq!(params.powers_of_tau.len(), 4); // 0, 1, 2, 3 = 4 powers
    }
}
