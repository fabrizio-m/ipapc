use ark_ec::{
    short_weierstrass_jacobian::{GroupAffine, GroupProjective},
    AffineCurve, ModelParameters, SWModelParameters,
};
use ark_ff::{Field, One, PrimeField};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, Evaluations, Radix2EvaluationDomain,
};
use commit::CommitmentTrait;
use itertools::Itertools;
pub use prove::{Commitment, HidingOpening, Opening, UnsafeHidingCommitment};
use rand::{prelude::StdRng, Rng, SeedableRng};
use std::{
    cell::RefCell,
    convert::identity,
    fmt::Debug,
    iter::{repeat, successors},
};

pub mod amortization;
mod challenges;
mod commit;
mod fft;
mod homomorphism;
mod open;
pub mod prove;
#[cfg(test)]
mod tests;
mod utils;
mod verify;

//type Poly<Fr> = DensePolynomial<Fr>;
type Fr<P> = <GroupAffine<P> as AffineCurve>::ScalarField;
pub struct IpaScheme<P, R>
where
    P: ModelParameters + SWModelParameters,
    R: Rng,
{
    basis: Vec<GroupAffine<P>>,
    ///second basis to commit to evals linearly
    evaluation_basis: Option<Vec<GroupAffine<P>>>,
    blinding_basis: GroupAffine<P>,
    max_degree: usize,
    rng: RefCell<R>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Blinding<P: SWModelParameters>(Fr<P>);

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Init<T: SWModelParameters> {
    Seed(u64),
    Elements(Vec<GroupAffine<T>>, GroupAffine<T>),
}

impl<P, R> IpaScheme<P, R>
where
    P: ModelParameters + SWModelParameters,
    Fr<P>: One,
    R: Rng,
{
    pub fn init(init: Init<P>, max_size: u8, commit_to_evals: bool, rng: R) -> Self {
        let max_degree = 2_usize.pow(max_size as u32);
        let (basis, blinding_basis) = init.to_elements(max_degree);
        let scheme = Self {
            basis,
            evaluation_basis: None,
            max_degree,
            blinding_basis,
            rng: RefCell::new(rng),
        };
        match commit_to_evals {
            true => {
                let eval_basis = scheme
                    .lagrange_commitments()
                    .into_iter()
                    .map(|point| point.0)
                    .collect_vec();
                Self {
                    evaluation_basis: Some(eval_basis),
                    ..scheme
                }
            }
            false => scheme,
        }
    }
    fn commit_simple(&self, poly: impl Into<CoeffsOrEvals<P>>) -> GroupProjective<P> {
        let poly: CoeffsOrEvals<P> = poly.into();
        let (poly, basis) = self.poly_to_msm_vecs(poly);
        debug_assert_eq!(poly.len(), self.max_degree);
        let bases = &*basis;
        let coeffs = poly.into_iter().map(|e| e.into_repr()).collect::<Vec<_>>();
        let scalars = &*coeffs;
        let result = ark_ec::msm::VariableBaseMSM::multi_scalar_mul(bases, scalars);
        result
    }
    fn poly_to_msm_vecs(&self, poly: CoeffsOrEvals<P>) -> (Vec<Fr<P>>, &Vec<GroupAffine<P>>) {
        //let basis = &self.basis;
        match (poly, &self.evaluation_basis) {
            (CoeffsOrEvals::Coeffs(coeffs), _) => (coeffs, &self.basis),
            (CoeffsOrEvals::Evals(evals), None) => {
                //interpolate
                let domain = Radix2EvaluationDomain::<Fr<P>>::new(evals.len()).unwrap();
                let evals = Evaluations::from_vec_and_domain(evals, domain);
                (evals.interpolate().coeffs, &self.basis)
            }
            (CoeffsOrEvals::Evals(evals), Some(basis)) => (evals, basis),
        }
    }
    pub fn commit<C: CommitmentTrait<P, R>>(&self, poly: impl Into<CoeffsOrEvals<P>>) -> C {
        C::commit(self, poly)
    }
    fn b(&self, z: Fr<P>) -> Vec<Fr<P>> {
        successors(Some(<Fr<P>>::one()), |previous| Some(*previous * z))
            .take(self.basis.len())
            .collect()
    }
    ///returns the reference string as commitments
    pub fn string(&self) -> Vec<Commitment<P, false>> {
        self.basis
            .iter()
            .cloned()
            .map(|elem| Commitment(elem))
            .collect()
    }
}

impl<T: SWModelParameters> Init<T> {
    fn to_elements(self, size: usize) -> (Vec<GroupAffine<T>>, GroupAffine<T>) {
        match self {
            Init::Seed(seed) => {
                let mut rng = StdRng::seed_from_u64(seed);
                let mut elems = repeat(())
                    .map(|_| {
                        let bytes: [u8; 32] = rng.gen();
                        let x = <T::BaseField as Field>::from_random_bytes(&bytes)?;
                        GroupAffine::<T>::get_point_from_x(x, false)
                    })
                    .filter_map(identity);
                let blind = elems.next().unwrap();
                (elems.take(size).collect(), blind)
            }
            Init::Elements(mut elems, blinding) => {
                assert!(elems.len() >= size);
                elems.truncate(size);
                assert_eq!(elems.len(), size);
                (elems, blinding)
            }
        }
    }
}

pub struct Assert<const COND: bool> {}

pub trait IsTrue {}
pub trait IsFalse {}

impl IsTrue for Assert<true> {}
impl IsFalse for Assert<false> {}

pub enum CoeffsOrEvals<P: SWModelParameters> {
    Coeffs(Vec<Fr<P>>),
    Evals(Vec<Fr<P>>),
}
impl<P: SWModelParameters> From<Vec<Fr<P>>> for CoeffsOrEvals<P> {
    fn from(coeffs: Vec<Fr<P>>) -> Self {
        Self::Coeffs(coeffs)
    }
}

impl<P: SWModelParameters> From<DensePolynomial<Fr<P>>> for CoeffsOrEvals<P> {
    fn from(poly: DensePolynomial<Fr<P>>) -> Self {
        Self::Coeffs(poly.coeffs)
    }
}

impl<P: SWModelParameters> From<Evaluations<Fr<P>>> for CoeffsOrEvals<P> {
    fn from(evals: Evaluations<Fr<P>>) -> Self {
        Self::Evals(evals.evals)
    }
}
