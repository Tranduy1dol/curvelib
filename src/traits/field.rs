use mathlib::field::montgomery::MontgomeryParams;

pub trait Field<'a>: Sized + Clone + Copy + PartialEq + Eq {
    fn zero(params: &'a MontgomeryParams) -> Self;
    fn is_zero(&self) -> bool;
    fn one(params: &'a MontgomeryParams) -> Self;
    fn inv(&self) -> Option<Self>;
    fn double(&self) -> Self;
    fn mul(&self, rhs: &Self) -> Self;
    fn add(&self, rhs: &Self) -> Self;
    fn square(&self) -> Self;
}
