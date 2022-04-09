use crate::{Commitment, Fr, IpaScheme};
use ark_ec::{
    short_weierstrass_jacobian::{GroupAffine, GroupProjective},
    AffineCurve, ModelParameters, ProjectiveCurve, SWModelParameters,
};
use ark_ff::{Field, One};
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use itertools::{izip, Itertools};

impl<P> IpaScheme<P>
where
    P: ModelParameters + SWModelParameters,
    Fr<P>: One,
{
    pub fn lagrange_commitments(&self) -> Vec<Commitment<P, false>> {
        let basis = self.basis.iter().map(|e| e.into_projective()).collect_vec();
        let commitments = self
            .ifft(basis)
            .into_iter()
            .map(|e| Commitment(e))
            .collect_vec();
        commitments
    }

    fn ifft(&self, coeffs: Vec<GroupProjective<P>>) -> Vec<GroupAffine<P>> {
        let inverted_size = Radix2EvaluationDomain::<Fr<P>>::new(coeffs.len())
            .unwrap()
            .size_as_field_element()
            .inverse()
            .unwrap();
        let mut res = self.fft(coeffs).into_iter();
        let first = res.by_ref().next().unwrap();
        [first]
            .into_iter()
            .chain(res.rev())
            .map(|e| e.into_affine().mul(inverted_size).into_affine())
            .collect_vec()
    }
    fn fft(&self, coeffs: Vec<GroupProjective<P>>) -> Vec<GroupProjective<P>> {
        let len = coeffs.len();
        assert!(len.is_power_of_two());
        if coeffs.len().is_one() {
            return coeffs;
        }
        let (even, odd): (Vec<GroupProjective<P>>, Vec<GroupProjective<P>>) = coeffs
            .iter()
            .cloned()
            .enumerate()
            .partition_map(|(i, item)| match i % 2 == 0 {
                true => itertools::Either::Left(item),
                false => itertools::Either::Right(item),
            });
        let [even, odd] = [even, odd].map(|coeffs| self.fft(coeffs));
        let mut result = coeffs;
        let (left, right) = result.split_at_mut(len / 2);
        let domain = {
            let domain = Radix2EvaluationDomain::<Fr<P>>::new(len).unwrap();
            domain.elements()
        };
        //let domain = successors(Some(domain), |(a, b)| Some((*a * gen, fb * gen)));
        for (even, odd, left, right, domain) in izip!(even, odd, left, right, domain) {
            let rhs = odd.into_affine().mul(domain);
            *left = even + rhs;
            *right = even - rhs;
        }
        result
    }
}

#[test]
fn lagrange_commitment() {
    use ark_pallas::PallasParameters;
    use ark_poly::{EvaluationDomain, Evaluations, GeneralEvaluationDomain};
    let scheme = IpaScheme::<PallasParameters>::init(crate::Init::Seed(1), 3, true);
    let domain = GeneralEvaluationDomain::new(8).unwrap();
    let evals = [1_i32, 0, 0, 0, 0, 0, 0, 0];
    let lcommit = |i| {
        let evals = evals
            .iter()
            .cloned()
            .map(Fr::<PallasParameters>::from)
            .cycle()
            .skip(evals.len() - i)
            .take(evals.len())
            .collect_vec();
        let poly = Evaluations::from_vec_and_domain(evals, domain).interpolate();
        let good_commit = scheme.commit(poly.coeffs);
        good_commit
    };
    let good_commitments = (0..8).map(lcommit).collect_vec();
    let commitments = scheme.lagrange_commitments();
    for (a, b) in good_commitments.iter().zip(commitments.iter()) {
        assert_eq!(a, b);
    }
}
