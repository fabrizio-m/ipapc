use crate::{challenges::ChallengeGenerator, Assert, Commitment, Fr, IpaScheme, IsFalse, Opening};
use ark_ec::{
    short_weierstrass_jacobian::{GroupAffine, GroupProjective},
    AffineCurve, ModelParameters, ProjectiveCurve, SWModelParameters,
};
use ark_ff::{One, Zero};
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial, SparsePolynomial},
    Polynomial, UVPolynomial,
};
use rand::Rng;
use std::{iter::successors, ops::Mul};

#[derive(Clone)]
pub struct MultiOpening<P: SWModelParameters> {
    openings: Vec<(Opening<P>, GroupAffine<P>)>,
    batch_opening: Opening<P>,
}

impl<P, R> IpaScheme<P, R>
where
    P: ModelParameters + SWModelParameters,
    Fr<P>: One,
    R: Rng,
{
    pub fn batch_open<const HIDING: bool>(
        &self,
        opens: Vec<(Commitment<P, HIDING>, &[Fr<P>], Fr<P>, Fr<P>)>,
    ) -> MultiOpening<P>
    where
        Assert<HIDING>: IsFalse,
    {
        let len = opens.len();
        let (openings, amortization, challenges) = opens
            .into_iter()
            .map(|open| {
                let (commitment, a, point, eval) = open;
                let u = ChallengeGenerator::inner_product_basis(&commitment, &point);
                let basis = &*self.basis;
                let b = self.b(point);
                let first = Self::round(basis, a, &*b, u, Some(vec![]));
                let rounds = vec![(first.lj, first.rj)];
                let open = Self::open_recursive(first, rounds, point, eval, u);
                open
            })
            .fold(
                (
                    Vec::with_capacity(len),
                    Vec::with_capacity(len),
                    ChallengeGenerator::new(),
                ),
                |mut state, open| {
                    let (open, s_challenges, basis) = open;
                    state.0.push(open);
                    state.1.push((s_challenges, basis));
                    state.2.digest_for_amortization(basis);
                    state
                },
            );
        let (point, combinations) = challenges.amortization_elements();
        let mut bs = vec![];
        let (_, s_poly, commitment) = amortization.into_iter().fold(
            (
                Fr::<P>::one(),
                DensePolynomial::from_coefficients_slice(&[]),
                GroupProjective::zero(),
            ),
            |acc, val| {
                let (combination, s_poly, commitment) = acc;
                let (s_challenges, basis) = val;
                let s_poly = Self::challenges_to_poly(s_challenges.unwrap(), combination) + s_poly;

                let commitment = basis.mul(combination) + commitment;
                bs.push(basis);

                (combination * combinations, s_poly, commitment)
            },
        );
        let eval = s_poly.evaluate(&point);
        let commitment = Commitment(commitment.into_affine());
        let batch_opening = self.open(commitment, &*s_poly.coeffs, point, eval);
        debug_assert_eq!(openings.len(), bs.len());
        let openings = openings.into_iter().zip(bs.into_iter()).collect();
        MultiOpening {
            openings,
            batch_opening,
        }
    }
    pub fn batch_verify<const HIDING: bool>(
        &self,
        commitments: &[Commitment<P, HIDING>],
        multi_open: MultiOpening<P>,
    ) -> Option<Vec<Fr<P>>>
    where
        Commitment<P, HIDING>: Clone + Copy,
    {
        assert_eq!(commitments.len(), multi_open.openings.len());
        let MultiOpening {
            openings,
            batch_opening,
        } = multi_open;
        let challenge_generator = openings.iter().fold(
            ChallengeGenerator::new(),
            |mut challenge_generator, elem| {
                challenge_generator.digest_for_amortization(elem.1.clone());
                challenge_generator
            },
        );
        let (combination_point, combination_element) = challenge_generator.amortization_elements();

        let a: (Vec<_>, Vec<_>) = openings
            .into_iter()
            .zip(commitments.iter())
            .map(|((open, final_basis), commitment)| {
                let Opening {
                    point,
                    eval,
                    rounds,
                    a,
                } = open;
                let u = ChallengeGenerator::inner_product_basis(&commitment, &point);
                let (final_commit, b_poly) = Self::process_rounds(*commitment, point, eval, rounds);
                let b = Self::eval_b_poly(&b_poly, point);
                let amortization_eval = Self::eval_b_poly(&b_poly, combination_point);
                let basis = final_basis;

                (
                    (final_commit, basis.mul(a) + u.mul(a * b), eval),
                    (amortization_eval, final_basis),
                )
            })
            .unzip();
        let (opens, amorti) = a;
        let combinations = successors(Some(Fr::<P>::one()), |e| Some(*e * combination_element));
        let (amortization_eval, amortization_commitment) = amorti
            .into_iter()
            .zip(combinations)
            .map(|elem| {
                let ((eval, commitment), combination) = elem;
                (eval * combination, commitment.mul(combination))
            })
            .reduce(|(a_e, a_c), (b_e, b_c)| (a_e + b_e, a_c + b_c))
            .unwrap();
        let commitment: Commitment<_, false> = Commitment(amortization_commitment.into_affine());
        let eval = Self::verify(&self, commitment, batch_opening)?;
        if amortization_eval != eval {
            return None;
        }
        let mut all_valid = true;
        let openings = opens
            .into_iter()
            .map(|(left, right, eval)| {
                if left != right {
                    all_valid = false;
                }
                eval
            })
            .collect();
        if all_valid {
            Some(openings)
        } else {
            None
        }
    }
    fn challenges_to_poly(
        challenges: Vec<(Fr<P>, Fr<P>)>,
        combination_element: Fr<P>,
    ) -> DensePolynomial<Fr<P>> {
        let powers = std::iter::successors(Some(1_usize), |a| Some(a * 2));
        let poly = challenges
            .into_iter()
            .rev()
            .zip(powers)
            .map(|((challenge, inverse), exp)| {
                SparsePolynomial::from_coefficients_vec(vec![(0, inverse), (exp, challenge)])
            })
            .reduce(|a, b| a.mul(&b))
            .unwrap();
        let poly: DensePolynomial<_> = DenseOrSparsePolynomial::from(poly).into();
        poly.mul(combination_element)
    }
}

#[test]
fn test_multi() {
    use crate::{tests::commit_and_open, Init};
    use ark_pallas::PallasParameters;
    use rand::thread_rng;
    let scheme = IpaScheme::<PallasParameters, _>::init(Init::Seed(1), 8, true, thread_rng());
    let (commitments, opens) = (0..4)
        .map(|_| {
            let a = commit_and_open(&scheme);
            (a.0.clone(), a)
        })
        .unzip::<_, _, Vec<_>, Vec<_>>();
    let (polys, opens) = opens
        .into_iter()
        .map(|open| {
            let (a, b, c, d) = open;
            ((b), (a, c, d))
        })
        .unzip::<_, _, Vec<_>, Vec<_>>();
    let opens = opens
        .into_iter()
        .zip(polys.iter())
        .map(|(open, poly)| {
            let (commit, point, eval) = open;
            (commit, &**poly, point, eval)
        })
        .collect::<Vec<_>>();

    let multi_open = scheme.batch_open(opens);
    let verif = scheme.batch_verify(&*commitments, multi_open);
    assert!(verif.is_some())
}
