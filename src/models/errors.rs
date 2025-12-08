use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SignatureError {
    InvalidPrivateKey,
    InvalidNonce,
    SignatureGenerationFailed, // For example if r or s is zero
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
