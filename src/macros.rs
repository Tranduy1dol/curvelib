//! Macros for defining field extensions and curve types.

/// Define a quadratic extension field Fp2 = Fp\[u\] / (u² + 1)
///
/// Now generic over C: FieldConfig instead of using lifetime 'a.
#[macro_export]
macro_rules! def_fp2 {
    ($name:ident, $base:ty) => {
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        pub struct $name<C: mathlib::FieldConfig> {
            pub c0: mathlib::FieldElement<C>,
            pub c1: mathlib::FieldElement<C>,
        }

        impl<C: mathlib::FieldConfig> $name<C> {
            pub fn new(c0: mathlib::FieldElement<C>, c1: mathlib::FieldElement<C>) -> Self {
                Self { c0, c1 }
            }

            pub fn zero() -> Self {
                Self {
                    c0: mathlib::FieldElement::<C>::zero(),
                    c1: mathlib::FieldElement::<C>::zero(),
                }
            }

            pub fn one() -> Self {
                Self {
                    c0: mathlib::FieldElement::<C>::one(),
                    c1: mathlib::FieldElement::<C>::zero(),
                }
            }
        }

        impl<C: mathlib::FieldConfig> std::ops::Add for $name<C> {
            type Output = Self;
            fn add(self, rhs: Self) -> Self {
                Self {
                    c0: self.c0 + rhs.c0,
                    c1: self.c1 + rhs.c1,
                }
            }
        }

        impl<C: mathlib::FieldConfig> std::ops::Sub for $name<C> {
            type Output = Self;
            fn sub(self, rhs: Self) -> Self {
                Self {
                    c0: self.c0 - rhs.c0,
                    c1: self.c1 - rhs.c1,
                }
            }
        }

        impl<C: mathlib::FieldConfig> std::ops::Neg for $name<C> {
            type Output = Self;
            fn neg(self) -> Self {
                Self {
                    c0: -self.c0,
                    c1: -self.c1,
                }
            }
        }

        impl<C: mathlib::FieldConfig> std::ops::Mul for $name<C> {
            type Output = Self;
            fn mul(self, rhs: Self) -> Self {
                // (a + bu)(c + du) = (ac - bd) + (ad + bc)u
                // where u² = -1
                let ac = self.c0 * rhs.c0;
                let bd = self.c1 * rhs.c1;
                let ad = self.c0 * rhs.c1;
                let bc = self.c1 * rhs.c0;
                Self {
                    c0: ac - bd,
                    c1: ad + bc,
                }
            }
        }
    };
}

/// Define a sextic extension field Fp6 = Fp2\[v\] / (v³ - ξ)
#[macro_export]
macro_rules! def_fp6 {
    ($name:ident, $base:ty) => {
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        pub struct $name<C: mathlib::FieldConfig> {
            pub c0: $base,
            pub c1: $base,
            pub c2: $base,
        }

        impl<C: mathlib::FieldConfig> $name<C> {
            pub fn new(c0: $base, c1: $base, c2: $base) -> Self {
                Self { c0, c1, c2 }
            }
        }

        impl<C: mathlib::FieldConfig> std::ops::Add for $name<C> {
            type Output = Self;
            fn add(self, rhs: Self) -> Self {
                Self {
                    c0: self.c0 + rhs.c0,
                    c1: self.c1 + rhs.c1,
                    c2: self.c2 + rhs.c2,
                }
            }
        }

        impl<C: mathlib::FieldConfig> std::ops::Sub for $name<C> {
            type Output = Self;
            fn sub(self, rhs: Self) -> Self {
                Self {
                    c0: self.c0 - rhs.c0,
                    c1: self.c1 - rhs.c1,
                    c2: self.c2 - rhs.c2,
                }
            }
        }

        impl<C: mathlib::FieldConfig> std::ops::Neg for $name<C> {
            type Output = Self;
            fn neg(self) -> Self {
                Self {
                    c0: -self.c0,
                    c1: -self.c1,
                    c2: -self.c2,
                }
            }
        }
    };
}

/// Define a Short Weierstrass curve: y² = x³ + ax + b
///
/// Now generic over C: FieldConfig instead of using lifetime 'a.
#[macro_export]
macro_rules! def_weierstrass_curve {
    ($name:ident, $field:ty) => {
        #[derive(Clone, Debug)]
        pub struct $name<C: mathlib::FieldConfig> {
            pub a: mathlib::FieldElement<C>,
            pub b: mathlib::FieldElement<C>,
            pub generator_x: mathlib::FieldElement<C>,
            pub generator_y: mathlib::FieldElement<C>,
            _marker: std::marker::PhantomData<C>,
        }

        impl<C: mathlib::FieldConfig> $name<C> {
            pub fn new(
                a: mathlib::FieldElement<C>,
                b: mathlib::FieldElement<C>,
                generator_x: mathlib::FieldElement<C>,
                generator_y: mathlib::FieldElement<C>,
            ) -> Self {
                Self {
                    a,
                    b,
                    generator_x,
                    generator_y,
                    _marker: std::marker::PhantomData,
                }
            }
        }
    };
}

// Keep the old macros for backward compatibility during migration
// These will be removed once migration is complete

#[macro_export]
macro_rules! def_fp2_legacy {
    ($name:ident, $base:ty) => {
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        pub struct $name<'a> {
            pub c0: $base,
            pub c1: $base,
        }

        impl<'a> $name<'a> {
            pub fn new(c0: $base, c1: $base) -> Self {
                Self { c0, c1 }
            }
        }

        impl<'a> std::ops::Add for $name<'a> {
            type Output = Self;
            fn add(self, rhs: Self) -> Self {
                Self {
                    c0: self.c0 + rhs.c0,
                    c1: self.c1 + rhs.c1,
                }
            }
        }

        impl<'a> std::ops::Sub for $name<'a> {
            type Output = Self;
            fn sub(self, rhs: Self) -> Self {
                Self {
                    c0: self.c0 - rhs.c0,
                    c1: self.c1 - rhs.c1,
                }
            }
        }

        impl<'a> std::ops::Neg for $name<'a> {
            type Output = Self;
            fn neg(self) -> Self {
                Self {
                    c0: -self.c0,
                    c1: -self.c1,
                }
            }
        }
    };
}

#[macro_export]
macro_rules! def_fp6_legacy {
    ($name:ident, $base:ty) => {
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        pub struct $name<'a> {
            pub c0: $base,
            pub c1: $base,
            pub c2: $base,
        }

        impl<'a> $name<'a> {
            pub fn new(c0: $base, c1: $base, c2: $base) -> Self {
                Self { c0, c1, c2 }
            }
        }

        impl<'a> std::ops::Add for $name<'a> {
            type Output = Self;
            fn add(self, rhs: Self) -> Self {
                Self {
                    c0: self.c0 + rhs.c0,
                    c1: self.c1 + rhs.c1,
                    c2: self.c2 + rhs.c2,
                }
            }
        }

        impl<'a> std::ops::Sub for $name<'a> {
            type Output = Self;
            fn sub(self, rhs: Self) -> Self {
                Self {
                    c0: self.c0 - rhs.c0,
                    c1: self.c1 - rhs.c1,
                    c2: self.c2 - rhs.c2,
                }
            }
        }

        impl<'a> std::ops::Neg for $name<'a> {
            type Output = Self;
            fn neg(self) -> Self {
                Self {
                    c0: -self.c0,
                    c1: -self.c1,
                    c2: -self.c2,
                }
            }
        }
    };
}

#[macro_export]
macro_rules! def_weierstrass_curve_legacy {
    ($name:ident, $field:ty) => {
        #[derive(Clone, Debug)]
        pub struct $name<'a> {
            pub a: $field,
            pub b: $field,
            pub params: &'a mathlib::MontgomeryContext,
            pub scalar_params: &'a mathlib::MontgomeryContext,
            pub generator_x: $field,
            pub generator_y: $field,
        }

        impl<'a> $name<'a> {
            pub fn new(
                a: $field,
                b: $field,
                params: &'a mathlib::MontgomeryContext,
                scalar_params: &'a mathlib::MontgomeryContext,
                generator_x: $field,
                generator_y: $field,
            ) -> Self {
                Self {
                    a,
                    b,
                    params,
                    scalar_params,
                    generator_x,
                    generator_y,
                }
            }
        }
    };
}
