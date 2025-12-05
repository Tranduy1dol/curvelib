use std::fmt::Debug;

use mathlib::U1024;
use mathlib::field::element::FieldElement;

pub trait Curve<'a>: Clone + Debug {
    type Point: ProjectivePoint;

    fn identity(&self) -> Self::Point;
    fn is_on_curve(&self, x: &FieldElement, y: &FieldElement) -> bool;
}

pub trait ProjectivePoint: Sized + Clone + Debug + PartialEq + Eq {
    fn is_identity(&self) -> bool;

    fn add(&self, rhs: &Self) -> Self;

    fn double(&self) -> Self;

    fn to_affine(&self) -> (FieldElement<'_>, FieldElement<'_>);

    fn mul(&self, scalar: &U1024) -> Self;
}
