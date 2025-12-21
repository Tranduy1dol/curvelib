//! ECDSA signature type.

use mathlib::U1024;

/// An ECDSA signature with (r, s) components.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Signature {
    /// The r component.
    r: U1024,
    /// The s component.
    s: U1024,
}

impl Signature {
    /// Create a new signature from r and s components.
    pub fn new(r: U1024, s: U1024) -> Self {
        Self { r, s }
    }

    /// Get the r component.
    pub fn r(&self) -> &U1024 {
        &self.r
    }

    /// Get the s component.
    pub fn s(&self) -> &U1024 {
        &self.s
    }

    /// Convert r to hex string.
    pub fn r_hex(&self) -> String {
        format_u1024_hex(&self.r)
    }

    /// Convert s to hex string.
    pub fn s_hex(&self) -> String {
        format_u1024_hex(&self.s)
    }

    /// Parse from hex strings. Returns None if parsing fails or strings are empty.
    pub fn from_hex(r_hex: &str, s_hex: &str) -> Option<Self> {
        let r_hex = r_hex.strip_prefix("0x").unwrap_or(r_hex);
        let r_hex = r_hex.strip_prefix("0X").unwrap_or(r_hex);
        let s_hex = s_hex.strip_prefix("0x").unwrap_or(s_hex);
        let s_hex = s_hex.strip_prefix("0X").unwrap_or(s_hex);

        // Basic validation: must be non-empty and only contain hex digits
        if r_hex.is_empty() || s_hex.is_empty() {
            return None;
        }

        for c in r_hex.chars().chain(s_hex.chars()) {
            if !c.is_ascii_hexdigit() {
                return None;
            }
        }

        let r = U1024::from_hex(r_hex);
        let s = U1024::from_hex(s_hex);
        Some(Signature::new(r, s))
    }
}

fn format_u1024_hex(val: &U1024) -> String {
    let mut hex = String::new();
    let mut started = false;
    for limb in val.0.iter().rev() {
        if *limb != 0 || started {
            if started {
                hex.push_str(&format!("{:016x}", limb));
            } else {
                hex.push_str(&format!("{:x}", limb));
                started = true;
            }
        }
    }
    if hex.is_empty() {
        hex.push('0');
    }
    format!("0x{}", hex)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_signature_hex() {
        let r = U1024::from_u64(12345);
        let s = U1024::from_u64(67890);
        let sig = Signature::new(r, s);

        let r_hex = sig.r_hex();
        let s_hex = sig.s_hex();

        let recovered = Signature::from_hex(&r_hex, &s_hex).unwrap();
        assert_eq!(sig, recovered);
    }

    #[test]
    fn test_signature_accessors() {
        let r = U1024::from_u64(0x12345);
        let s = U1024::from_u64(0x67890);
        let sig = Signature::new(r, s);

        assert_eq!(*sig.r(), r);
        assert_eq!(*sig.s(), s);
    }

    #[test]
    fn test_signature_invalid_hex() {
        // Empty strings
        assert!(Signature::from_hex("", "0x123").is_none());
        assert!(Signature::from_hex("0x123", "").is_none());

        // Invalid characters
        assert!(Signature::from_hex("0xGHI", "0x123").is_none());
        assert!(Signature::from_hex("0x123", "0xGHI").is_none());
    }
}
