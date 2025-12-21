//! Curve instances.
//!
//! This module provides pre-defined curve instances with their configurations.

pub mod bls6_6;
pub mod tiny_jubjub;

// Re-export commonly used items
pub use bls6_6::{
    Bls6_6BaseField, Bls6_6G1Config, Bls6_6ScalarField, FINAL_EXPONENT, Fp43, Fr13, G1Point,
    G2Point,
};
pub use tiny_jubjub::{Fp13, Fr5, TinyJubjubBaseField, TinyJubjubConfig, TinyJubjubScalarField};
