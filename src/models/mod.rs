//! Curve models.
//!
//! This module provides implementations for different elliptic curve models.

pub mod sextic_twist;
pub mod short_weierstrass;
pub mod twisted_edwards;

pub use sextic_twist::{SexticTwist, TwistPoint};
pub use short_weierstrass::{WeierstrassCurve, WeierstrassPoint};
pub use twisted_edwards::{EdwardsCurve, EdwardsPoint};

pub mod sw {
    //! Short Weierstrass curve types.
    pub use super::short_weierstrass::{Affine, Projective};
}

pub mod te {
    //! Twisted Edwards curve types.
    pub use super::twisted_edwards::{Affine, Projective};
}
