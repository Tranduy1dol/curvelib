//! Signature error types.

/// Error type for signature operations.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SignatureError {
    /// The nonce is invalid (zero or too large).
    InvalidNonce,
    /// The signature is malformed.
    MalformedSignature,
    /// The public key is invalid.
    InvalidPublicKey,
    /// The private key is invalid.
    InvalidPrivateKey,
    /// The scalar is out of range.
    ScalarOutOfRange,
}

impl std::fmt::Display for SignatureError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SignatureError::InvalidNonce => write!(f, "invalid nonce"),
            SignatureError::MalformedSignature => write!(f, "malformed signature"),
            SignatureError::InvalidPublicKey => write!(f, "invalid public key"),
            SignatureError::InvalidPrivateKey => write!(f, "invalid private key"),
            SignatureError::ScalarOutOfRange => write!(f, "scalar out of range"),
        }
    }
}

impl std::error::Error for SignatureError {}
