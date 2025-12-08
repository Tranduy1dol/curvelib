use std::fmt;

/// Errors that can occur during ECDSA signature generation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SignatureError {
    /// The private key provided is invalid (e.g., zero or >= curve order n).
    InvalidPrivateKey,
    /// The nonce provided is invalid (e.g., zero or >= curve order n).
    InvalidNonce,
    /// Signature generation failed due to a transient error (e.g. r or s computed to zero).
    /// This is statistically rare but can happen; retrying usually resolves it.
    SignatureGenerationFailed,
}

impl fmt::Display for SignatureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SignatureError::InvalidPrivateKey => {
                write!(f, "Invalid private key: must be in range [1, n-1]")
            }
            SignatureError::InvalidNonce => write!(f, "Invalid nonce: must be in range [1, n-1]"),
            SignatureError::SignatureGenerationFailed => {
                write!(f, "Signature generation failed (r or s was zero). Retry.")
            }
        }
    }
}

impl std::error::Error for SignatureError {}
