//! Hex encoding/decoding for keys.

use crate::traits::{CurveConfig, ProjectivePoint};
use mathlib::FieldConfig;

use super::keypair::{PrivateKey, PublicKey};

/// Error type for hex parsing.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum HexError {
    /// Invalid hex string length.
    InvalidLength { expected: usize, got: usize },
    /// Invalid hex character.
    InvalidCharacter(char),
    /// Invalid format or prefix.
    InvalidFormat(String),
}

impl std::fmt::Display for HexError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HexError::InvalidLength { expected, got } => {
                write!(f, "invalid length: expected {}, got {}", expected, got)
            }
            HexError::InvalidCharacter(c) => {
                write!(f, "invalid hex character: '{}'", c)
            }
            HexError::InvalidFormat(msg) => {
                write!(f, "invalid format: {}", msg)
            }
        }
    }
}

impl std::error::Error for HexError {}

/// Trait for types that can be encoded to hex.
pub trait ToHex {
    /// Encode to hex string.
    fn to_hex(&self, compressed: bool) -> String;
}

/// Trait for types that can be decoded from hex.
pub trait FromHex: Sized {
    /// Decode from hex string.
    fn from_hex(hex: &str) -> Result<Self, HexError>;
}

impl<C: CurveConfig> ToHex for PrivateKey<C> {
    fn to_hex(&self, _compressed: bool) -> String {
        format!("0x{}", self.to_hex_string())
    }
}

impl<C: CurveConfig> FromHex for PrivateKey<C> {
    fn from_hex(hex_str: &str) -> Result<Self, HexError> {
        // Remove optional 0x prefix
        let hex_str = hex_str.strip_prefix("0x").unwrap_or(hex_str);
        let hex_str = hex_str.strip_prefix("0X").unwrap_or(hex_str);

        // Validate hex characters
        for c in hex_str.chars() {
            if !c.is_ascii_hexdigit() {
                return Err(HexError::InvalidCharacter(c));
            }
        }

        PrivateKey::<C>::from_hex_string(hex_str)
            .ok_or_else(|| HexError::InvalidFormat("failed to parse private key".to_string()))
    }
}

impl<C: CurveConfig> ToHex for PublicKey<C> {
    fn to_hex(&self, compressed: bool) -> String {
        let (ax, ay) = self.point().to_affine();

        if self.point().is_identity() {
            return "0x00".to_string();
        }

        // Compute field byte length from the modulus
        let field_byte_len = compute_field_byte_length::<C>();

        // Get x coordinate as hex with fixed width
        let x = ax.to_u1024();
        let x_hex = format_u1024_fixed_width(&x, field_byte_len);

        if compressed {
            // Compressed: prefix (02/03) + x coordinate (fixed width)
            // 02 if y is even, 03 if y is odd
            let y = ay.to_u1024();
            let prefix = if y.0[0] & 1 == 0 { "02" } else { "03" };
            format!("0x{}{}", prefix, x_hex)
        } else {
            // Uncompressed: 04 + x + y (both fixed width)
            let y = ay.to_u1024();
            let y_hex = format_u1024_fixed_width(&y, field_byte_len);
            format!("0x04{}{}", x_hex, y_hex)
        }
    }
}

/// Compute the field byte length for a given curve configuration.
/// This determines how many bytes are needed to represent a field element.
fn compute_field_byte_length<C: CurveConfig>() -> usize {
    // MODULUS is a constant on FieldConfig
    let modulus = C::BaseField::MODULUS;

    // Find the highest non-zero limb to determine the byte length
    let mut byte_len = 0;
    for (i, &limb) in modulus.0.iter().enumerate().rev() {
        if limb != 0 {
            // This limb contributes up to 8 bytes.
            // We need to find how many bytes are actually used in this limb.
            let bits_used = 64 - u64::leading_zeros(limb);
            let bytes_in_limb = bits_used.div_ceil(8) as usize;
            byte_len = i * 8 + bytes_in_limb;
            break;
        }
    }

    // Ensure at least 1 byte
    byte_len.max(1)
}

/// Format U1024 as hex with fixed width (zero-padded to exactly byte_len * 2 hex chars).
fn format_u1024_fixed_width(val: &mathlib::U1024, byte_len: usize) -> String {
    let mut bytes = Vec::with_capacity(byte_len);

    // Limbs are little-endian: val.0[0] is the least significant 64 bits
    for i in (0..byte_len).rev() {
        let limb_idx = i / 8;
        let byte_idx = i % 8;
        if limb_idx < val.0.len() {
            let byte = ((val.0[limb_idx] >> (byte_idx * 8)) & 0xFF) as u8;
            bytes.push(byte);
        } else {
            bytes.push(0);
        }
    }

    // Convert bytes to hex string (already in big-endian order)
    bytes.iter().map(|b| format!("{:02x}", b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instances::bls6_6::Bls6_6G1Config;
    use mathlib::U1024;

    #[test]
    fn test_private_key_hex_roundtrip() {
        let scalar = U1024::from_u64(5); // Within order 13
        let key = PrivateKey::<Bls6_6G1Config>::new(scalar);
        let hex_str = key.to_hex(false);
        let recovered = PrivateKey::<Bls6_6G1Config>::from_hex(&hex_str).unwrap();
        assert_eq!(key.to_u1024(), recovered.to_u1024());
    }

    #[test]
    fn test_private_key_from_hex_with_prefix() {
        let key = PrivateKey::<Bls6_6G1Config>::from_hex("0x5").unwrap();
        assert_eq!(key.to_u1024(), U1024::from_u64(5));
    }

    #[test]
    fn test_private_key_from_hex_without_prefix() {
        let key = PrivateKey::<Bls6_6G1Config>::from_hex("7").unwrap();
        assert_eq!(key.to_u1024(), U1024::from_u64(7));
    }

    #[test]
    fn test_invalid_hex_character() {
        let result = PrivateKey::<Bls6_6G1Config>::from_hex("0xGHIJ");
        assert!(matches!(result, Err(HexError::InvalidCharacter(_))));
    }

    #[test]
    fn test_public_key_hex_fixed_width() {
        let g = <Bls6_6G1Config as CurveConfig>::generator();
        let pk = PublicKey::<Bls6_6G1Config>::new(g);

        // For BLS6_6, field modulus is 43 (1 byte)
        // Generator x = 13 (0x0d), y = 15 (0x0f)
        // Fixed-width hex should be 2 chars per coordinate

        let uncompressed = pk.to_hex(false);
        // 0x04 + 0d + 0f = 0x040d0f
        assert_eq!(uncompressed, "0x040d0f");

        let compressed = pk.to_hex(true);
        // y=15 is odd -> prefix 03
        // 0x03 + 0d = 0x030d
        assert_eq!(compressed, "0x030d");
    }

    #[test]
    fn test_public_key_identity_hex() {
        // Create an identity point by multiplying generator by 0
        let id_point = <Bls6_6G1Config as CurveConfig>::generator().mul(&U1024::from_u8(0));
        let id_pk = PublicKey::<Bls6_6G1Config>::new(id_point);
        assert_eq!(id_pk.to_hex(false), "0x00");
    }
}
