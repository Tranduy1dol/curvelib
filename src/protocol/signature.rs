use mathlib::U1024;

/// ECDSA signature with (r, s) components.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Signature {
    pub r: U1024,
    pub s: U1024,
}

impl Signature {
    pub fn new(r: U1024, s: U1024) -> Self {
        Self { r, s }
    }
}
