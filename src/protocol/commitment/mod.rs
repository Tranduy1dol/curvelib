//! Polynomial commitment scheme traits and implementations.
//!
//! This module provides generic traits for polynomial commitment schemes
//! and concrete implementations like KZG.

pub mod kzg;

use mathlib::{FieldConfig, FieldElement, Polynomial};

/// Core trait for polynomial commitment schemes.
///
/// A polynomial commitment scheme allows a prover to commit to a polynomial
/// and later prove evaluations at specific points without revealing the polynomial.
pub trait PolynomialCommitment<C: FieldConfig> {
    /// Type representing a commitment to a polynomial.
    type Commitment;

    /// Type representing a proof of evaluation.
    type Proof;

    /// Type representing the public parameters (SRS, etc.)
    type Params;

    /// Commit to a polynomial.
    ///
    /// # Parameters
    ///
    /// * `params` - Public parameters for the scheme
    /// * `polynomial` - The polynomial to commit to
    ///
    /// # Returns
    ///
    /// A commitment to the polynomial.
    fn commit(params: &Self::Params, polynomial: &Polynomial<C>) -> Self::Commitment;

    /// Create an opening proof for a polynomial evaluation.
    ///
    /// Proves that `polynomial(point) = value`.
    ///
    /// # Parameters
    ///
    /// * `params` - Public parameters
    /// * `polynomial` - The committed polynomial
    /// * `point` - The evaluation point
    ///
    /// # Returns
    ///
    /// A tuple of `(proof, value)` where `value = polynomial(point)`.
    fn open(
        params: &Self::Params,
        polynomial: &Polynomial<C>,
        point: &FieldElement<C>,
    ) -> (Self::Proof, FieldElement<C>);

    /// Verify an opening proof.
    ///
    /// # Parameters
    ///
    /// * `params` - Public parameters
    /// * `commitment` - The polynomial commitment
    /// * `point` - The evaluation point
    /// * `value` - The claimed evaluation value
    /// * `proof` - The opening proof
    ///
    /// # Returns
    ///
    /// `true` if the proof is valid, `false` otherwise.
    fn verify(
        params: &Self::Params,
        commitment: &Self::Commitment,
        point: &FieldElement<C>,
        value: &FieldElement<C>,
        proof: &Self::Proof,
    ) -> bool;
}

/// Extension trait for batch polynomial commitments.
///
/// Allows committing to multiple polynomials at once or verifying
/// multiple openings simultaneously.
pub trait BatchCommitment<C: FieldConfig>: PolynomialCommitment<C> {
    /// Commit to multiple polynomials at once.
    fn batch_commit(params: &Self::Params, polynomials: &[Polynomial<C>]) -> Vec<Self::Commitment>;

    /// Verify multiple opening proofs at once.
    ///
    /// This can be more efficient than verifying each proof individually.
    fn batch_verify(
        params: &Self::Params,
        commitments: &[Self::Commitment],
        points: &[FieldElement<C>],
        values: &[FieldElement<C>],
        proofs: &[Self::Proof],
    ) -> bool;
}
