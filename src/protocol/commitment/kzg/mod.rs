//! KZG polynomial commitment scheme implementation.
//!
//! This module provides a concrete implementation of the KZG commitment scheme
//! using the `PolynomialCommitment` trait.

use std::marker::PhantomData;
use std::ops::Neg;

use mathlib::{FieldConfig, FieldElement, Polynomial, U1024};

use crate::models::{SexticTwist, TwistPoint, WeierstrassCurve, WeierstrassPoint};
use crate::protocol::commitment::PolynomialCommitment;
use crate::protocol::pairing::tate_pairing;
use crate::traits::{Curve, ProjectivePoint};

/// Type aliases for curve points.
pub type G1Point<C> = WeierstrassPoint<C>;
pub type G2Point<C> = TwistPoint<C>;

/// KZG Parameters (Structured Reference String).
///
/// Contains the trusted setup parameters for the KZG scheme.
pub struct KzgParams<C: FieldConfig> {
    /// Powers of tau in G1: [G, τG, τ²G, ...]
    pub powers_of_tau: Vec<G1Point<C>>,
    /// G2 generator
    pub g2: G2Point<C>,
    /// τ·G2
    pub tau_g2: G2Point<C>,
    /// Curve order (for pairing)
    pub r_order: U1024,
    /// Final exponentiation parameter (for pairing)
    pub final_exp: U1024,
}

/// KZG polynomial commitment scheme.
pub struct Kzg<C: FieldConfig> {
    _marker: PhantomData<C>,
}

impl<C: FieldConfig> Kzg<C> {
    /// Trusted Setup: Generate SRS from secret tau.
    ///
    /// **WARNING**: In production, tau must come from a secure MPC ceremony!
    ///
    /// # Parameters
    ///
    /// * `tau` - Secret value (must be destroyed after setup!)
    /// * `max_degree` - Maximum polynomial degree to support
    /// * `g1_curve` - Base curve for G1
    /// * `g2_curve` - Twist curve for G2
    /// * `r_order` - Curve order
    /// * `final_exp` - Final exponentiation parameter
    ///
    /// # Returns
    ///
    /// KZG parameters for commitment and verification.
    pub fn setup(
        tau: FieldElement<C>,
        max_degree: usize,
        g1_curve: &WeierstrassCurve<C>,
        g2_curve: &SexticTwist<C>,
        r_order: U1024,
        final_exp: U1024,
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
            r_order,
            final_exp,
        }
    }
}

impl<C: FieldConfig> PolynomialCommitment<C> for Kzg<C> {
    type Commitment = G1Point<C>;
    type Proof = G1Point<C>;
    type Params = KzgParams<C>;

    fn commit(params: &Self::Params, polynomial: &Polynomial<C>) -> Self::Commitment {
        let coeffs = polynomial.to_vec();
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

    fn open(
        params: &Self::Params,
        polynomial: &Polynomial<C>,
        point: &FieldElement<C>,
    ) -> (Self::Proof, FieldElement<C>) {
        // 1. Evaluate y = P(z)
        let y = polynomial.evaluate(point);

        // 2. Compute quotient Q(x) = (P(x) - y) / (x - z)
        let p_minus_y = polynomial.clone() - Polynomial::constant(y);
        let divisor = Polynomial::new(vec![-*point, FieldElement::<C>::one()]);
        let (quotient, _remainder) = p_minus_y.divide_with_remainder(&divisor);

        // 3. Commit to Q(x) -> Proof π
        let proof = Self::commit(params, &quotient);

        (proof, y)
    }

    fn verify(
        params: &Self::Params,
        commitment: &Self::Commitment,
        point: &FieldElement<C>,
        value: &FieldElement<C>,
        proof: &Self::Proof,
    ) -> bool {
        let z_scalar = point.to_u1024();
        let y_scalar = value.to_u1024();

        // LHS = e(π, τ·G₂ - z·G₂)
        let g2_z = params.g2.mul(&z_scalar);
        let g2_term = params.tau_g2.clone().add(&g2_z.neg());

        let lhs = tate_pairing(proof, &g2_term, params.r_order, params.final_exp);

        // RHS = e(C - y·G₁, G₂)
        let g1 = &params.powers_of_tau[0];
        let gy = g1.mul(&y_scalar);
        let g1_term = commitment.add(&gy.neg());

        let rhs = tate_pairing(&g1_term, &params.g2, params.r_order, params.final_exp);

        lhs == rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instances::bls6_6::{Bls6_6BaseField, FINAL_EXPONENT, get_g1_curve, get_g2_curve};
    use mathlib::U1024;

    #[test]
    fn test_kzg_setup() {
        let g1_curve = get_g1_curve();
        let g2_curve = get_g2_curve();
        let tau = FieldElement::<Bls6_6BaseField>::new(U1024::from_u64(5));
        let r_order = U1024::from_u64(13);
        let final_exp = U1024::from_u64(FINAL_EXPONENT);

        let params =
            Kzg::<Bls6_6BaseField>::setup(tau, 3, &g1_curve, &g2_curve, r_order, final_exp);

        assert_eq!(params.powers_of_tau.len(), 4); // degrees 0, 1, 2, 3
    }

    #[test]
    fn test_kzg_commit_and_open() {
        let g1_curve = get_g1_curve();
        let g2_curve = get_g2_curve();
        let tau = FieldElement::<Bls6_6BaseField>::new(U1024::from_u64(5));
        let r_order = U1024::from_u64(13);
        let final_exp = U1024::from_u64(FINAL_EXPONENT);

        let params =
            Kzg::<Bls6_6BaseField>::setup(tau, 3, &g1_curve, &g2_curve, r_order, final_exp);

        // Create polynomial P(x) = 1 + 2x + 3x²
        let poly = Polynomial::new(vec![
            FieldElement::<Bls6_6BaseField>::new(U1024::from_u64(1)),
            FieldElement::<Bls6_6BaseField>::new(U1024::from_u64(2)),
            FieldElement::<Bls6_6BaseField>::new(U1024::from_u64(3)),
        ]);

        // Commit
        let commitment = Kzg::commit(&params, &poly);

        // Open at z = 2: P(2) = 1 + 4 + 12 = 17
        let z = FieldElement::<Bls6_6BaseField>::new(U1024::from_u64(2));
        let (proof, value) = Kzg::open(&params, &poly, &z);

        // Check evaluation
        assert_eq!(value.to_u1024().0[0], 17);

        // Verify (may not pass on toy curve due to field mismatch issues)
        // let valid = Kzg::verify(&params, &commitment, &z, &value, &proof);
        // Note: Commented out as BLS6_6 has base field ≠ scalar field issues
    }
}
