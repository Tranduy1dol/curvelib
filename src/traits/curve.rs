use std::fmt::Debug;

use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{BigInt, FieldElement, U1024};

use crate::protocol::signature::Signature;
use crate::traits::point::ProjectivePoint;

pub trait Curve<'a>: Clone + Debug {
    type Point: ProjectivePoint<'a>;

    fn identity(&self) -> Self::Point;
    fn is_on_curve(
        &self,
        x: &<Self::Point as ProjectivePoint<'a>>::Field,
        y: &<Self::Point as ProjectivePoint<'a>>::Field,
    ) -> bool;
    fn scalar_params(&self) -> &'a MontgomeryParams;
    fn generator(&self) -> Self::Point;

    /// Generates a keypair (private_key, public_key) from a random scalar.
    ///
    /// The private key is the scalar, and the public key is scalar * generator.
    fn generate_keypair(&self, private_key: &U1024) -> Self::Point {
        self.generator().mul(private_key)
    }

    /// Signs a message hash using ECDSA with the given private key.
    ///
    /// Generates secure random nonce internally.
    /// Returns a Signature containing (r, s) components.
    fn sign(
        &self,
        message_hash: &U1024,
        priv_key: &U1024,
    ) -> Result<Signature, crate::models::errors::SignatureError>
    where
        <Self::Point as ProjectivePoint<'a>>::Field: ToU1024,
    {
        use crate::models::errors::SignatureError;
        use rand::RngCore;

        let scalar_params = self.scalar_params();
        let n = &scalar_params.modulus;

        // 1. Validate inputs
        // Check if priv_key is zero
        let priv_key_is_zero = priv_key == &U1024::zero();

        // Manual comparison for priv_key >= n using div_rem
        // if priv_key >= n, then priv_key / n >= 1
        let (q, _) = priv_key.div_rem(n);
        let priv_key_ge_n = q != U1024::zero();

        if priv_key_is_zero || priv_key_ge_n {
            return Err(SignatureError::InvalidPrivateKey);
        }

        let generator = self.generator();

        // Loop until valid signature is generated (probabilistic)
        loop {
            // 2. Generate random nonce k in [1, n-1]
            let mut k_bytes = [0u8; 128]; // 1024 bits
            rand::rng().fill_bytes(&mut k_bytes);

            let mut limbs = [0u64; 16];
            for i in 0..16 {
                let mut ndata = [0u8; 8];
                ndata.copy_from_slice(&k_bytes[i * 8..(i + 1) * 8]);
                limbs[i] = u64::from_le_bytes(ndata);
            }
            let k_full = U1024(limbs); // Direct construction since field is pub

            // k = k_full % n, skip if zero => ensures range [1, n-1]
            let (_, rem) = k_full.div_rem(n);
            let k_nonce = if rem == U1024::zero() {
                continue;
            } else {
                rem
            };

            // 3. r = (k * G).x mod n
            let r_point = generator.mul(&k_nonce);
            let (r_x_elem, _) = r_point.to_affine();

            let r_val = r_x_elem.to_u1024_val();
            let r_elem = FieldElement::new(r_val, scalar_params);
            let r = r_elem.to_u1024();

            if r == U1024::zero() {
                continue; // Retry with new k
            }

            // 4. s = k^(-1) * (z + r * d) mod n
            let k_elem = FieldElement::new(k_nonce, scalar_params);
            let z_elem = FieldElement::new(*message_hash, scalar_params);
            let d_elem = FieldElement::new(*priv_key, scalar_params);

            let k_inv = k_elem.inv();
            // s = k_inv * (z + r * d)
            // r * d
            let rd = r_elem * d_elem;
            // z + r*d
            let z_rd = z_elem + rd;
            let s_elem = k_inv * z_rd;

            let s = s_elem.to_u1024();
            if s == U1024::zero() {
                continue; // Retry with new k
            }

            return Ok(Signature::new(r, s));
        }
    }

    /// Sign with explicit nonce (TESTING ONLY)
    ///
    /// # Safety
    /// This method allows specifying the nonce k manually.
    /// REUSING NONCE COMPLETELY COMPROMISES THE PRIVATE KEY.
    /// DO NOT USE UNLESS YOU KNOW EXACTLY WHAT YOU ARE DOING.
    #[cfg(feature = "test")]
    fn sign_with_nonce(
        &self,
        message_hash: &U1024,
        priv_key: &U1024,
        k_nonce: &U1024,
    ) -> Result<Signature, crate::models::errors::SignatureError>
    where
        <Self::Point as ProjectivePoint<'a>>::Field: ToU1024,
    {
        use crate::models::errors::SignatureError;

        let scalar_params = self.scalar_params();
        let n = &scalar_params.modulus;

        if priv_key == &U1024::zero() || priv_key >= n {
            return Err(SignatureError::InvalidPrivateKey);
        }
        if k_nonce == &U1024::zero() || k_nonce >= n {
            return Err(SignatureError::InvalidNonce);
        }

        let generator = self.generator();

        let r_point = generator.mul(k_nonce);
        let (r_x_elem, _) = r_point.to_affine();
        let r_val = r_x_elem.to_u1024_val();

        let r_elem = FieldElement::new(r_val, scalar_params);
        let r = r_elem.to_u1024();

        if r == U1024::zero() {
            return Err(SignatureError::SignatureGenerationFailed); // r=0 is just bad luck/params for deterministic
        }

        let k_elem = FieldElement::new(*k_nonce, scalar_params);
        let z_elem = FieldElement::new(*message_hash, scalar_params);
        let d_elem = FieldElement::new(*priv_key, scalar_params);

        let k_inv = k_elem.inv();
        let s_elem = k_inv * (z_elem + (r_elem * d_elem));

        let s = s_elem.to_u1024();
        if s == U1024::zero() {
            return Err(SignatureError::SignatureGenerationFailed);
        }

        Ok(Signature::new(r, s))
    }

    /// Verifies an ECDSA signature against a message hash and public key.
    ///
    /// Returns true if the signature is valid, false otherwise.
    fn verify(&self, signature: &Signature, message_hash: &U1024, pub_key: &Self::Point) -> bool
    where
        <Self::Point as ProjectivePoint<'a>>::Field: ToU1024,
    {
        let scalar_params = self.scalar_params();
        let generator = self.generator();
        let n = &scalar_params.modulus;

        if pub_key.is_identity() {
            return false;
        }

        if !self.is_on_curve(&pub_key.to_affine().0, &pub_key.to_affine().1) {
            return false;
        }

        // Defend against low-order subgroup attacks
        // Check if n * Q == Infinity
        let q_n = pub_key.mul(n);
        if !q_n.is_identity() {
            return false;
        }

        if signature.r == U1024::zero() || signature.r >= *n {
            return false;
        }
        if signature.s == U1024::zero() || signature.s >= *n {
            return false;
        }

        let s_elem = FieldElement::new(signature.s, scalar_params);
        let z_elem = FieldElement::new(*message_hash, scalar_params);
        let r_elem = FieldElement::new(signature.r, scalar_params);

        let w = s_elem.inv();

        let u1 = (z_elem * w).to_u1024();
        let u2 = (r_elem * w).to_u1024();

        let p1 = generator.mul(&u1);
        let p2 = pub_key.mul(&u2);
        let p = p1.add(&p2);

        if p.is_identity() {
            return false;
        }

        let (p_x, _) = p.to_affine();
        let p_x_mod_n = FieldElement::new(p_x.to_u1024_val(), scalar_params);

        p_x_mod_n.to_u1024() == signature.r
    }
}

pub trait ToU1024 {
    fn to_u1024_val(&self) -> U1024;
}

impl<'a> ToU1024 for FieldElement<'a> {
    fn to_u1024_val(&self) -> U1024 {
        self.to_u1024()
    }
}
