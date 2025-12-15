#[macro_export]
macro_rules! def_fp2 {
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
macro_rules! def_fp6 {
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
macro_rules! def_weierstrass_curve {
    ($name:ident, $field:ty) => {
        #[derive(Clone, Debug)]
        pub struct $name<'a> {
            pub a: $field,
            pub b: $field,
            pub params: &'a mathlib::field::montgomery::MontgomeryParams,
            pub scalar_params: &'a mathlib::field::montgomery::MontgomeryParams,
            pub generator_x: $field,
            pub generator_y: $field,
        }

        impl<'a> $name<'a> {
            pub fn new(
                a: $field,
                b: $field,
                params: &'a mathlib::field::montgomery::MontgomeryParams,
                scalar_params: &'a mathlib::field::montgomery::MontgomeryParams,
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
