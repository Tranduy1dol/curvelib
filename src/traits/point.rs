use std::fmt::Debug;

use mathlib::U1024;

use crate::traits::field::Field;

pub trait ProjectivePoint<'a>: Sized + Clone + Debug + PartialEq + Eq {
    type Field: Field<'a>;

    fn is_identity(&self) -> bool;

    fn add(&self, rhs: &Self) -> Self;

    fn double(&self) -> Self;

    fn to_affine(&self) -> (Self::Field, Self::Field);

    fn mul(&self, scalar: &U1024) -> Self;
}
