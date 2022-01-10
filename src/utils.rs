use crate::Fr;
use ark_ec::{
    short_weierstrass_jacobian::{GroupAffine, GroupProjective},
    AffineCurve, ProjectiveCurve, SWModelParameters,
};
use ark_ff::Field;
use std::ops::Add;

pub fn compress_basis<P: SWModelParameters>(
    left: &[GroupAffine<P>],
    right: &[GroupAffine<P>],
    challenge: Fr<P>,
) -> Vec<GroupAffine<P>> {
    assert_eq!(left.len(), right.len());
    let inverse = challenge.inverse().unwrap();
    let left = left.iter().map(|elem| elem.mul(inverse));
    let right = right.iter().map(|elem| elem.mul(challenge));
    left.zip(right)
        .map(|(a, b)| (a + b).into_affine())
        .collect()
}
pub fn compress<P: SWModelParameters>(
    left: &[Fr<P>],
    right: &[Fr<P>],
    challenge: Fr<P>,
) -> Vec<Fr<P>> {
    assert_eq!(left.len(), right.len());
    let inverse = challenge.inverse().unwrap();
    let left = left.iter().map(|elem| *elem * inverse);
    let right = right.iter().map(|elem| *elem * challenge);
    left.zip(right).map(|(a, b)| a + b).collect()
}

pub fn inner_product<P: SWModelParameters>(
    a: &[GroupAffine<P>],
    b: &[Fr<P>],
) -> GroupProjective<P> {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b.iter())
        .map(|(a, b)| a.mul(*b))
        .reduce(Add::add)
        .unwrap()
}
pub fn scalar_inner_product<P: SWModelParameters>(a: &[Fr<P>], b: &[Fr<P>]) -> Fr<P> {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b.iter())
        .map(|(a, b)| *a * *b)
        .reduce(Add::add)
        .unwrap()
}
pub fn split<T>(slice: &[T]) -> (&[T], &[T]) {
    let len = slice.len();
    assert_eq!(len % 2, 0);
    (&slice[0..len / 2], &slice[len / 2..])
}
