//! Traits for elliptic curve cryptography.
//!
//! This module provides the core trait definitions used throughout curvelib.

pub mod config;
pub mod curve;
pub mod field;
pub mod point;

pub use config::{CurveConfig, ShortWeierstrassConfig, TwistedEdwardsConfig};
pub use curve::{Curve, ToU1024};
pub use field::Field;
pub use point::ProjectivePoint;
