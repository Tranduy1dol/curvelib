use std::fmt::Debug;

use mathlib::field::element::FieldElement;

pub trait CurveParams: Clone + Debug {
    fn zero_point(&self) -> FieldElement<'_>;
    fn is_on_curve(&self, x: &FieldElement, y: &FieldElement) -> bool;
}

pub trait ProjectivePoint: Sized + Clone + Debug + PartialEq + Eq {
    fn identity() -> Self;

    fn is_identity(&self) -> bool;

    fn add(&self, rhs: &Self) -> Self;

    fn double(&self) -> Self;

    fn to_affine(&self) -> (FieldElement<'_>, FieldElement<'_>);
}
