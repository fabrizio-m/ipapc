use ark_ec::{short_weierstrass_jacobian::GroupAffine, AffineCurve};
use ark_pallas::PallasParameters;
use ark_poly::{Polynomial, UVPolynomial};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use ipapc::{Commitment, Init, IpaScheme};
use rand::thread_rng;

type Fr = <GroupAffine<PallasParameters> as AffineCurve>::ScalarField;
const MAX: usize = 8;
const SIZE: usize = 12;

pub(crate) fn commit_and_open(
    scheme: &IpaScheme<PallasParameters>,
) -> (Commitment<PallasParameters, false>, Vec<Fr>, Fr, Fr) {
    use ark_ff::UniformRand;

    type Poly = ark_poly::univariate::DensePolynomial<Fr>;
    let mut rng = thread_rng();
    let poly = [0; 2_usize.pow(SIZE as u32)]
        .map(|_| Fr::rand(&mut rng))
        .to_vec();
    let commit = scheme.commit(poly.clone());
    let point = Fr::from(5);
    let eval = {
        let poly = Poly::from_coefficients_slice(&*poly);
        poly.evaluate(&point)
    };
    (commit, poly, point, eval)
}
fn sample(
    scheme: &IpaScheme<PallasParameters>,
    size: usize,
) -> (
    Vec<Vec<Fr>>,
    Vec<Commitment<PallasParameters, false>>,
    Vec<(Commitment<PallasParameters, false>, Fr, Fr)>,
) {
    let (commitments, opens) = (0..=size)
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
    (polys, commitments, opens)
}

pub fn batch_verify(c: &mut Criterion) {
    let mut group = c.benchmark_group("batch_verify");
    let scheme = IpaScheme::<PallasParameters>::init(Init::Seed(1), SIZE as u8);

    let (polys, commitments, opens) = sample(&scheme, 2_usize.pow(MAX as u32));
    for size in 0..=MAX {
        let max = 2_usize.pow(size as u32);
        group.throughput(Throughput::Elements(max as u64));

        let opens = opens[0..max].to_vec();
        let opens = opens
            .into_iter()
            .zip(polys.iter())
            .map(|(open, poly)| {
                let (commit, point, eval) = open;
                (commit, &**poly, point, eval)
            })
            .collect::<Vec<_>>();
        let multi_open = scheme.batch_open(opens);
        let commitments = commitments[0..max].to_vec();

        group.bench_with_input(BenchmarkId::from_parameter(max), &size, |b, _size| {
            b.iter_batched(
                || (multi_open.clone(), commitments.clone()),
                |(multi_open, commitments)| {
                    let evals = scheme.batch_verify(&*commitments, multi_open);
                    (evals, commitments)
                },
                criterion::BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = batch_verify
);
criterion_main!(benches);
