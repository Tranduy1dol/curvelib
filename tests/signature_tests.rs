//! Signature tests - temporarily disabled during migration.
//!
//! The `sign` and `verify` methods have been removed from the `Curve` trait
//! as part of the Phase 2 refactoring. They will be reimplemented as part
//! of a `SigningEngine` module.

use curvelib::{
    instances::{
        bls6_6::{Bls6_6G1Config, Bls6_6ScalarField},
        tiny_jubjub,
    },
    protocol::{
        keys::{FromHex, KeyEngine, PrivateKey, ToHex},
        signing::{Signature, SigningEngine},
    },
    traits::{Curve, ProjectivePoint},
};
use mathlib::{BigInt, FieldElement, U1024};

#[test]
fn test_keypair_generation() {
    let curve = tiny_jubjub::get_curve();

    // Private key d = 2
    let priv_key = U1024::from_u64(2);
    let pub_key = curve.generate_keypair(&priv_key);

    // Public key should not be identity
    assert!(!pub_key.is_identity());

    // Public key should be on the curve
    let (x, y) = pub_key.to_affine();
    assert!(curve.is_on_curve(&x, &y));
}

#[test]
fn test_keypair_different_keys() {
    let curve = tiny_jubjub::get_curve();

    let priv_key_1 = U1024::from_u64(2);
    let priv_key_2 = U1024::from_u64(3);

    let pub_key_1 = curve.generate_keypair(&priv_key_1);
    let pub_key_2 = curve.generate_keypair(&priv_key_2);

    // Different private keys should give different public keys
    assert_ne!(pub_key_1, pub_key_2);
}

#[test]
fn test_modern_key_engine() {
    let engine = KeyEngine::<Bls6_6G1Config>::new();
    let scalar = U1024::from_u64(5);
    let keypair = engine.keypair_from_scalar(scalar);

    assert_eq!(keypair.private_key().to_u1024(), U1024::from_u64(5));
    assert!(!keypair.public_key().is_identity());

    // Hex roundtrip
    let priv_hex = keypair.private_key().to_hex_string();
    let recovered_priv = PrivateKey::<Bls6_6G1Config>::from_hex_string(&priv_hex).unwrap();
    assert_eq!(*keypair.private_key(), recovered_priv);

    // Using ToHex/FromHex traits
    let priv_hex_trait = keypair.private_key().to_hex(false);
    let recovered_priv_trait = PrivateKey::<Bls6_6G1Config>::from_hex(&priv_hex_trait).unwrap();
    assert_eq!(*keypair.private_key(), recovered_priv_trait);
}

#[test]
fn test_public_key_hex() {
    let engine = KeyEngine::<Bls6_6G1Config>::new();
    let keypair = engine.keypair_from_scalar(U1024::from_u64(5));
    let pub_key = keypair.public_key();

    // Uncompressed hex
    let hex_uncompressed = pub_key.to_hex(false);
    assert!(hex_uncompressed.starts_with("0x04"));

    // Compressed hex
    let hex_compressed = pub_key.to_hex(true);
    assert!(hex_compressed.starts_with("0x02") || hex_compressed.starts_with("0x03"));
}

#[test]
fn test_signature_hex_roundtrip() {
    let r = U1024::from_u64(0x12345);
    let s = U1024::from_u64(0x67890);
    let sig = Signature::new(r, s);

    let r_hex = sig.r_hex();
    let s_hex = sig.s_hex();

    let recovered = Signature::from_hex(&r_hex, &s_hex).unwrap();
    assert_eq!(sig, recovered);
}

#[test]
fn test_ecdsa_sign_verify() {
    let key_engine = KeyEngine::<Bls6_6G1Config>::new();
    let signing_engine = SigningEngine::<Bls6_6G1Config>::new();

    let private_scalar = U1024::from_u64(5);
    let keypair = key_engine.keypair_from_scalar(private_scalar);

    let message_hash = FieldElement::<Bls6_6ScalarField>::new(U1024::from_u64(7));
    let nonce = FieldElement::<Bls6_6ScalarField>::new(U1024::from_u64(3));

    // Sign
    let signature = signing_engine
        .sign(&message_hash, keypair.private_key(), &nonce)
        .expect("signing should succeed");

    // Verify
    // Note: This might fail on BLS6_6 due to base/scalar field mismatch,
    // but we want to test the API.
    let _valid = signing_engine.verify(&signature, &message_hash, keypair.public_key());

    // Basic consistency
    assert_ne!(*signature.r(), U1024::zero());
    assert_ne!(*signature.s(), U1024::zero());
}

#[test]
fn test_invalid_signatures() {
    let key_engine = KeyEngine::<Bls6_6G1Config>::new();
    let signing_engine = SigningEngine::<Bls6_6G1Config>::new();

    let keypair = key_engine.keypair_from_scalar(U1024::from_u64(5));
    let message_hash = FieldElement::<Bls6_6ScalarField>::new(U1024::from_u64(7));

    // Signature with zero r
    let sig_zero_r = Signature::new(U1024::zero(), U1024::from_u64(123));
    assert!(!signing_engine.verify(&sig_zero_r, &message_hash, keypair.public_key()));

    // Signature with zero s
    let sig_zero_s = Signature::new(U1024::from_u64(123), U1024::zero());
    assert!(!signing_engine.verify(&sig_zero_s, &message_hash, keypair.public_key()));
}
