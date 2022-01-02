use ark_ec::{
    short_weierstrass_jacobian::GroupAffine, AffineCurve, ModelParameters, ProjectiveCurve,
    SWModelParameters,
};
use ark_ff::{Field, One};
use ark_poly::univariate::DensePolynomial;
use rand::{prelude::StdRng, Rng, SeedableRng};
use std::{
    convert::identity,
    fmt::Debug,
    iter::{repeat, successors},
    ops::Add,
};

use crate::challenges::ChallengeGenerator;

mod challenges;
#[cfg(test)]
mod tests;

type Poly<Fr> = DensePolynomial<Fr>;
type Fr<P> = <GroupAffine<P> as AffineCurve>::ScalarField;
pub struct IpaScheme<P>
where
    P: ModelParameters + SWModelParameters,
{
    basis: Vec<GroupAffine<P>>,
    max_degree: usize,
}

impl<P> IpaScheme<P>
where
    P: ModelParameters + SWModelParameters,
    Fr<P>: One,
    Commitment<P>: Clone + Copy,
{
    pub fn init(init: Init<P>, max_size: u8) -> Self {
        let max_degree = 2_usize.pow(max_size as u32);
        println!("m: {}", max_degree);
        let basis = init.to_elements(max_degree);
        Self { basis, max_degree }
    }
    pub fn commit<'a>(&self, coeffs: Vec<Fr<P>>) -> Commitment<P> {
        assert_eq!(coeffs.len(), self.max_degree);
        let point = self
            .basis
            .iter()
            .zip(coeffs.into_iter())
            .map(|(g, f)| g.mul(f))
            .reduce(Add::add)
            .unwrap();
        Commitment(point.into_affine())
    }
    fn b(&self, z: Fr<P>) -> Vec<Fr<P>> {
        successors(Some(<Fr<P>>::one()), |previous| Some(*previous * z))
            .take(self.basis.len())
            .collect()
    }
    fn round(
        basis: &[GroupAffine<P>],
        a: &[Fr<P>],
        b: &[Fr<P>],
        u: GroupAffine<P>,
    ) -> RoundOutput<P> {
        assert_eq!(basis.len(), a.len());
        assert_eq!(basis.len(), b.len());
        assert!(basis.len() > 1);
        let (a_l, a_r) = split(a);
        let (b_l, b_r) = split(b);
        let (g_l, g_r) = split(basis);

        let lj = inner_product(g_r, a_l) + u.mul(scalar_inner_product::<P>(a_l, b_r)).into_affine();
        let rj = inner_product(g_l, a_r) + u.mul(scalar_inner_product::<P>(a_r, b_l)).into_affine();

        let challenge = <ChallengeGenerator<P>>::round_challenge(&lj, &rj);
        let a = compress::<P>(a_l, a_r, challenge);
        let b = compress::<P>(b_r, b_l, challenge);
        let basis = compress_basis(g_r, g_l, challenge);
        RoundOutput {
            lj,
            rj,
            a,
            b,
            basis,
        }
    }
    fn open_recursive(
        prev: RoundOutput<P>,
        mut rounds: Vec<(GroupAffine<P>, GroupAffine<P>)>,
        point: Fr<P>,
        eval: Fr<P>,
        u: GroupAffine<P>,
    ) -> Opening<P> {
        let RoundOutput {
            lj,
            rj,
            a,
            b,
            basis,
        } = prev;
        println!("lj: {}", lj);
        println!("rj: {}", rj);
        if a.len().is_one() {
            println!("final b: {}", b[0]);
            println!("final basis: {}", basis[0]);
            let comm = basis[0].mul(a[0]) + u.mul(a[0] * b[0]);
            println!("p comm: {}", comm);

            Opening::<P> {
                a: a[0],
                rounds,
                point,
                eval,
            }
        } else {
            let prev = Self::round(&*basis, &*a, &*b, u);
            rounds.push((prev.lj, prev.rj));
            Self::open_recursive(prev, rounds, point, eval, u)
        }
    }
    pub fn open(
        &self,
        commitment: Commitment<P>,
        a: &[Fr<P>],
        point: Fr<P>,
        eval: Fr<P>,
    ) -> Opening<P> {
        let u = ChallengeGenerator::inner_product_basis(&commitment, &point);
        println!("U: {}", u);
        //let p = commitment.0 + u.mul(eval).into_affine();
        let basis = &*self.basis;
        let b = self.b(point);
        let first = Self::round(basis, a, &*b, u);
        let rounds = vec![(first.lj, first.rj)];
        let opening = Self::open_recursive(first, rounds, point, eval, u);
        opening
    }
    pub fn verify(&self, commitment: Commitment<P>, open: Opening<P>) -> Option<Fr<P>> {
        let Opening::<P> {
            point,
            eval,
            a,
            rounds,
        } = open;
        let u = ChallengeGenerator::inner_product_basis(&commitment, &point);
        println!("U: {}", u);
        let p = commitment.0.into_projective() + u.mul(eval);
        let basis = self.basis.clone();
        let b = self.b(point);

        let (final_commit, basis, b) = rounds.iter().fold((p, basis, b), |state, (lj, rj)| {
            let (p, basis, b) = state;
            let challenge = <ChallengeGenerator<P>>::round_challenge(lj, rj);
            let inverse = challenge.inverse().unwrap();
            println!("lj: {}", lj);
            println!("rj: {}", rj);
            let new_commit = p + (rj.mul(challenge.square()) + lj.mul(inverse.square()));

            let (b_l, b_r) = split(&*b);
            let (g_l, g_r) = split(&*basis);
            let b = compress::<P>(b_r, b_l, challenge);
            let basis = compress_basis(g_r, g_l, challenge);

            (new_commit, basis, b)
        });
        println!("final b: {}", b[0]);
        println!("final basis: {}", basis[0]);

        println!("p2 comm: {}", basis[0].mul(a) + u.mul(a * b[0]));
        println!("p3 comm: {}", final_commit);
        let final_check = final_commit == basis[0].mul(a) + u.mul(a * b[0]);
        if final_check {
            Some(eval)
        } else {
            None
        }
    }
}
#[derive(Debug)]
pub struct Opening<P: SWModelParameters> {
    point: Fr<P>,
    eval: Fr<P>,
    rounds: Vec<(GroupAffine<P>, GroupAffine<P>)>,
    a: Fr<P>,
}
struct RoundOutput<P: SWModelParameters> {
    lj: GroupAffine<P>,
    rj: GroupAffine<P>,
    a: Vec<Fr<P>>,
    b: Vec<Fr<P>>,
    basis: Vec<GroupAffine<P>>,
}
fn compress_basis<P: SWModelParameters>(
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
fn compress<P: SWModelParameters>(left: &[Fr<P>], right: &[Fr<P>], challenge: Fr<P>) -> Vec<Fr<P>> {
    assert_eq!(left.len(), right.len());
    let inverse = challenge.inverse().unwrap();
    let left = left.iter().map(|elem| *elem * inverse);
    let right = right.iter().map(|elem| *elem * challenge);
    left.zip(right).map(|(a, b)| a + b).collect()
}
fn inner_product<P: SWModelParameters>(a: &[GroupAffine<P>], b: &[Fr<P>]) -> GroupAffine<P> {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b.iter())
        .map(|(a, b)| a.mul(*b))
        .reduce(Add::add)
        .unwrap()
        .into_affine()
}
fn scalar_inner_product<P: SWModelParameters>(a: &[Fr<P>], b: &[Fr<P>]) -> Fr<P> {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b.iter())
        .map(|(a, b)| *a * *b)
        .reduce(Add::add)
        .unwrap()
}
fn split<T>(slice: &[T]) -> (&[T], &[T]) {
    let len = slice.len();
    assert_eq!(len % 2, 0);
    (&slice[0..len / 2], &slice[len / 2..])
}

#[derive(Clone, Copy, Debug)]
pub struct Commitment<T: SWModelParameters>(GroupAffine<T>)
where
    GroupAffine<T>: Debug;

pub enum Init<T: SWModelParameters> {
    Seed(u64),
    Elements(Vec<GroupAffine<T>>),
}

impl<T: SWModelParameters> Init<T> {
    fn to_elements(self, size: usize) -> Vec<GroupAffine<T>> {
        match self {
            Init::Seed(seed) => {
                let mut rng = StdRng::seed_from_u64(seed);
                let elems = repeat(())
                    .map(|_| {
                        let bytes: [u8; 32] = rng.gen();
                        let x = <T::BaseField as Field>::from_random_bytes(&bytes)?;
                        GroupAffine::<T>::get_point_from_x(x, false)
                    })
                    .filter_map(identity);
                elems.take(size).collect()
            }
            Init::Elements(mut elems) => {
                assert!(elems.len() >= size);
                elems.truncate(size);
                assert_eq!(elems.len(), size);
                elems
            }
        }
    }
}
