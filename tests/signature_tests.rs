use curvelib::{instances::tiny_jubjub, traits::Curve};
use mathlib::{BigInt, U1024};

#[test]
fn test_ecdsa_signature_verification() {
    let curve = tiny_jubjub::get_curve();

    // 1. Generate Keypair
    // Private key d = 2
    let priv_key = U1024::from_u64(2);
    let pub_key = curve.generate_keypair(&priv_key);

    // 2. Prepare Message (nonce is internal)
    // Message hash z = 11
    let message_hash = U1024::from_u64(11);

    // 3. Sign
    let signature = curve
        .sign(&message_hash, &priv_key)
        .expect("signing failed");

    println!("Signature: r={:?}, s={:?}", signature.r, signature.s);

    // 4. Verify
    let valid = curve.verify(&signature, &message_hash, &pub_key);

    assert!(valid, "Signature verification failed");
}

#[test]
fn test_ecdsa_invalid_signature() {
    let curve = tiny_jubjub::get_curve();
    let priv_key = U1024::from_u64(2);
    let pub_key = curve.generate_keypair(&priv_key);
    let message_hash = U1024::from_u64(11);

    let mut signature = curve
        .sign(&message_hash, &priv_key)
        .expect("signing failed");

    // Tamper with signature
    signature.s = signature.s + U1024::one();

    let valid = curve.verify(&signature, &message_hash, &pub_key);

    assert!(!valid, "Invalid signature should not verify");
}

#[test]
fn test_ecdsa_different_message() {
    let curve = tiny_jubjub::get_curve();
    let priv_key = U1024::from_u64(2);
    let pub_key = curve.generate_keypair(&priv_key);
    let message_hash = U1024::from_u64(11);

    let signature = curve
        .sign(&message_hash, &priv_key)
        .expect("signing failed");

    // Verify against a different message
    let wrong_message = U1024::from_u64(12);
    let valid = curve.verify(&signature, &wrong_message, &pub_key);

    assert!(!valid, "Signature should not verify for wrong message");
}

#[cfg(feature = "test")]
#[test]
fn test_ecdsa_deterministic_failure() {
    // This tests the test-only API for deterministic nonces
    let curve = tiny_jubjub::get_curve();
    let priv_key = U1024::from_u64(2);
    let message_hash = U1024::from_u64(11);

    // Invalid nonce (0)
    let k_zero = U1024::zero();
    let res = curve.sign_with_nonce(&message_hash, &priv_key, &k_zero);
    assert!(res.is_err());

    // Invalid nonce (>= n)
    let scalar_params = curve.scalar_params();
    let n = &scalar_params.modulus;
    let k_ge_n = *n;
    let res = curve.sign_with_nonce(&message_hash, &priv_key, &k_ge_n);
    assert!(res.is_err());
}
