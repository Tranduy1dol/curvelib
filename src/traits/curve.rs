use std::fmt::Debug;

use mathlib::{MontgomeryParams, U1024, fp, u1024};

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

    /// Computes the public key corresponding to a private scalar by multiplying the curve generator by that scalar.
    ///
    /// The provided `private_key` is treated as the scalar secret; callers must ensure it is a valid private key for the curve (e.g., greater than zero and less than the curve order).
    ///
    /// # Examples
    ///
    /// ```rust
    /// use curvelib::instances::tiny_jubjub;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use mathlib::u1024;
    ///
    /// let curve = tiny_jubjub::get_curve();
    /// let sk = u1024!(3);
    /// let pk = curve.generate_keypair(&sk);
    /// assert!(!pk.is_identity());
    /// ```
    fn generate_keypair(&self, private_key: &U1024) -> Self::Point {
        self.generator().mul(private_key)
    }

    /// Creates an ECDSA-style signature for `message_hash` using `priv_key`.
    ///
    /// The method validates the private key and produces a signature (r, s). It uses
    /// fresh cryptographic nonce for each attempt and will retry internally until
    /// a valid signature is produced or the private key is rejected.
    ///
    /// # Errors
    ///
    /// Returns `SignatureError::InvalidPrivateKey` when `priv_key` is zero or greater
    /// than or equal to the curve order.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # #[cfg(feature = "test")]
    /// # {
    /// use curvelib::instances::tiny_jubjub;
    /// use curvelib::traits::Curve;
    /// use mathlib::u1024;
    ///
    /// let curve = tiny_jubjub::get_curve();
    /// let msg = u1024!(1);
    /// let sk = u1024!(2);      // valid since scalar modulus is 5
    /// let k  = u1024!(1);      // test nonce in [1, n-1]
    ///
    /// let sig = curve.sign_with_nonce(&msg, &sk, &k).expect("failed to sign");
    /// assert_ne!(sig.r, u1024!(0));
    /// assert_ne!(sig.s, u1024!(0));
    /// # }
    /// ```
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
        let priv_key_is_zero = priv_key == &u1024!(0);

        // Manual comparison for priv_key >= n using div_rem
        // if priv_key >= n, then priv_key / n >= 1
        let (q, _) = priv_key.div_rem(n);
        let priv_key_ge_n = q != u1024!(0);

        if priv_key_is_zero || priv_key_ge_n {
            return Err(SignatureError::InvalidPrivateKey);
        }

        let generator = self.generator();

        // Loop until valid signature is generated (probabilistic)
        loop {
            // 2. Generate random nonce k in [1, n-1]
            let mut k_bytes = [0u8; 128]; // 1024 bits
            rand::rng().fill_bytes(&mut k_bytes);
            let k_full = U1024::from_le_bytes(&k_bytes);

            // k = k_full % n, skip if zero => ensures range [1, n-1]
            let (_, rem) = k_full.div_rem(n);
            let k_nonce = if rem == u1024!(0) {
                continue;
            } else {
                rem
            };

            // 3. r = (k * G).x mod n
            let r_point = generator.mul(&k_nonce);
            let (r_x_elem, _) = r_point.to_affine();

            let r_val = r_x_elem.to_u1024();
            let r_elem = fp!(r_val, scalar_params);
            let r = r_elem.to_u1024();

            if r == u1024!(0) {
                continue; // Retry with new k
            }

            // 4. s = k^(-1) * (z + r * d) mod n
            let k_elem = fp!(k_nonce, scalar_params);
            let z_elem = fp!(*message_hash, scalar_params);
            let d_elem = fp!(*priv_key, scalar_params);

            let k_inv = k_elem.inv();
            // s = k_inv * (z + r * d)
            // r * d
            let rd = r_elem * d_elem;
            // z + r*d
            let z_rd = z_elem + rd;
            let s_elem = k_inv * z_rd;

            let s = s_elem.to_u1024();
            if s == u1024!(0) {
                continue; // Retry with new k
            }

            return Ok(Signature::new(r, s));
        }
    }

    /// Generate an ECDSA-like signature for `message_hash` using the provided nonce `k_nonce` (TESTING ONLY).
    ///
    /// # Safety
    /// This method accepts an explicit nonce. Reusing or exposing the same nonce for multiple signatures
    /// completely compromises the private key. Do not use in production or unless you fully understand the risk.
    ///
    /// # Errors
    /// Returns `Err(SignatureError::InvalidPrivateKey)` if `priv_key` is zero or not less than the curve order,
    /// `Err(SignatureError::InvalidNonce)` if `k_nonce` is zero or not less than the curve order,
    /// and `Err(SignatureError::SignatureGenerationFailed)` if the generated `r` or `s` equals zero.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// // Example (conceptual): sign deterministically with a test-only nonce.
    /// // let sig = curve.sign_with_nonce(&message_hash, &private_key, &k_nonce).unwrap();
    /// // assert!(sig.r != 0 && sig.s != 0);
    /// ```ignore
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

        if priv_key == &u1024!(0) || priv_key >= n {
            return Err(SignatureError::InvalidPrivateKey);
        }
        if k_nonce == &u1024!(0) || k_nonce >= n {
            return Err(SignatureError::InvalidNonce);
        }

        let generator = self.generator();

        let r_point = generator.mul(k_nonce);
        let (r_x_elem, _) = r_point.to_affine();
        let r_val = r_x_elem.to_u1024();

        let r_elem = fp!(r_val, scalar_params);
        let r = r_elem.to_u1024();

        if r == u1024!(0) {
            return Err(SignatureError::SignatureGenerationFailed); // r=0 is just bad luck/params for deterministic
        }

        let k_elem = fp!(*k_nonce, scalar_params);
        let z_elem = fp!(*message_hash, scalar_params);
        let d_elem = fp!(*priv_key, scalar_params);

        let k_inv = k_elem.inv();
        let s_elem = k_inv * (z_elem + (r_elem * d_elem));

        let s = s_elem.to_u1024();
        if s == u1024!(0) {
            return Err(SignatureError::SignatureGenerationFailed);
        }

        Ok(Signature::new(r, s))
    }

    /// Verifies an ECDSA-style signature against a message hash and public key.
    ///
    /// Performs all canonical checks: rejects the identity or off-curve public keys, enforces
    /// that the public key lies in the prime-order subgroup (n*Q == identity), and verifies
    /// that r and s are in the valid range before performing the ECDSA verification equations.
    ///
    /// # Returns
    ///
    /// `true` if the signature is valid, `false` otherwise.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # #[cfg(feature = "test")]
    /// # {
    /// use curvelib::instances::tiny_jubjub;
    /// use curvelib::traits::{Curve, ProjectivePoint};
    /// use mathlib::u1024;
    ///
    /// let curve = tiny_jubjub::get_curve();
    /// let msg = u1024!(1);
    /// let sk = u1024!(2);
    /// let k  = u1024!(1);
    ///
    /// let pk = curve.generate_keypair(&sk);
    /// let sig = curve.sign_with_nonce(&msg, &sk, &k).unwrap();
    /// assert!(curve.verify(&sig, &msg, &pk));
    /// # }
    /// ```
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

        if signature.r == u1024!(0) || signature.r >= *n {
            return false;
        }
        if signature.s == u1024!(0) || signature.s >= *n {
            return false;
        }

        let s_elem = fp!(signature.s, scalar_params);
        let z_elem = fp!(*message_hash, scalar_params);
        let r_elem = fp!(signature.r, scalar_params);

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
        let p_x_mod_n = fp!(p_x.to_u1024(), scalar_params);

        p_x_mod_n.to_u1024() == signature.r
    }
}

pub trait ToU1024 {
    fn to_u1024(&self) -> U1024;
}
