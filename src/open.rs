use crate::{
    challenges::ChallengeGenerator, commit::CommitmentTrait, Commitment, Fr, HidingOpening,
    IpaScheme, Opening, UnsafeHidingCommitment,
};
use ark_ec::{AffineCurve, SWModelParameters};
use rand::Rng;

pub trait OpenTrait<P, R>
where
    P: SWModelParameters,
    Self::Commit: CommitmentTrait<P, R>,
    R: Rng,
{
    type Commit: CommitmentTrait<P, R>;

    fn open(
        scheme: &IpaScheme<P, R>,
        commitment: Self::Commit,
        coeffs: &[Fr<P>],
        point: Fr<P>,
        eval: Fr<P>,
    ) -> Self;
}

impl<P, R> OpenTrait<P, R> for Opening<P>
where
    P: SWModelParameters,
    R: Rng,
{
    type Commit = Commitment<P, false>;

    fn open(
        scheme: &IpaScheme<P, R>,
        commitment: Self::Commit,
        coeffs: &[Fr<P>],
        point: Fr<P>,
        eval: Fr<P>,
    ) -> Self {
        let u = ChallengeGenerator::inner_product_basis(&commitment, &point);
        let basis = &*scheme.basis;
        let b = scheme.b(point);
        //let mut rng = &self.rng;
        let first = IpaScheme::<P, R>::round(basis, coeffs, &*b, u, None);
        let rounds = vec![(first.lj, first.rj)];
        let (opening, _, _) = IpaScheme::<P, R>::open_recursive(first, rounds, point, eval, u);
        opening
    }
}

impl<P, R> OpenTrait<P, R> for HidingOpening<P>
where
    P: SWModelParameters,
    R: Rng,
{
    type Commit = UnsafeHidingCommitment<P>;

    fn open(
        scheme: &IpaScheme<P, R>,
        commitment: Self::Commit,
        coeffs: &[Fr<P>],
        point: Fr<P>,
        eval: Fr<P>,
    ) -> Self {
        //let rng = rng.unwrap();
        let rng = { &mut *scheme.rng.borrow_mut() };
        let UnsafeHidingCommitment(commitment, blinding) = commitment;
        let commitment = Commitment::<_, true>(commitment);
        let u = ChallengeGenerator::inner_product_basis(&commitment, &point);
        let basis = &*scheme.basis;
        let blinding_basis = scheme.blinding_basis;
        let b = scheme.b(point);
        //let mut rng = &self.rng;
        let first = IpaScheme::<P, R>::hiding_round(
            basis,
            coeffs,
            &*b,
            u,
            blinding_basis,
            rng,
            blinding,
            None,
        );
        let rounds = vec![(first.lj, first.rj)];
        let opening = IpaScheme::<P, R>::open_recursive_hiding(
            first,
            rounds,
            point,
            eval,
            u,
            blinding_basis,
            rng,
        );
        opening
    }
}

pub trait VerifTrait<P, R>
where
    P: SWModelParameters,
    R: Rng,
{
    type Commit;

    fn verify(self, scheme: &IpaScheme<P, R>, commitment: Self::Commit) -> Option<Fr<P>>;
}

impl<P, R> VerifTrait<P, R> for Opening<P>
where
    P: SWModelParameters,
    R: Rng,
{
    type Commit = Commitment<P, false>;

    fn verify(self, scheme: &IpaScheme<P, R>, commitment: Self::Commit) -> Option<Fr<P>> {
        let open = self;
        let Opening::<P> {
            point,
            eval,
            a,
            rounds,
        } = open;
        let (final_commit, check) = scheme.general_verify(commitment, point, eval, a, rounds);
        if final_commit == check {
            Some(eval)
        } else {
            None
        }
    }
}

impl<P, R> VerifTrait<P, R> for HidingOpening<P>
where
    P: SWModelParameters,
    R: Rng,
{
    type Commit = Commitment<P, true>;

    fn verify(self, scheme: &IpaScheme<P, R>, commitment: Self::Commit) -> Option<Fr<P>> {
        let open = self;
        let HidingOpening::<P> {
            point,
            eval,
            a,
            rounds,
            blinding_factor,
        } = open;
        let (final_commit, check) = scheme.general_verify(commitment, point, eval, a, rounds);
        if final_commit == check + scheme.blinding_basis.mul(blinding_factor) {
            Some(eval)
        } else {
            None
        }
    }
}
