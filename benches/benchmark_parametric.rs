use ark_ec::{short_weierstrass_jacobian::GroupAffine, AffineCurve};
use ark_pallas::PallasParameters;
use ark_poly::{Polynomial, UVPolynomial};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use ipapc::{Init, IpaScheme};
use rand::{prelude::ThreadRng, thread_rng, Rng};
use std::{iter::repeat, ops::RangeInclusive};

type Scheme = IpaScheme<PallasParameters>;
type Fr = <GroupAffine<PallasParameters> as AffineCurve>::ScalarField;

const RANGE: RangeInclusive<u8> = 8..=14;

fn sample(size: u8) -> (Scheme, Vec<Fr>, ThreadRng) {
    let scheme = IpaScheme::<PallasParameters>::init(Init::Seed(1), size);
    let mut rng = thread_rng();
    //let poly: [Fr<PallasParameters>; 2_usize.pow(SIZE as u32)] = rng.gen();
    let poly: Vec<Fr> = repeat(())
        .map(|_| rng.gen())
        .take(2_usize.pow(size as u32))
        .collect();
    (scheme, poly, rng)
}

pub fn commit(c: &mut Criterion) {
    let mut group = c.benchmark_group("commit");
    for size in RANGE {
        group.throughput(Throughput::Elements(size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), &size, |b, size| {
            let (scheme, poly, _rng) = sample(*size);

            let coeffs = black_box(poly.to_vec());
            b.iter(|| {
                let _commit = scheme.commit(black_box(coeffs.clone()));
            })
        });
    }
    group.finish();
}
pub fn open(c: &mut Criterion) {
    let mut group = c.benchmark_group("open");
    for size in RANGE {
        group.throughput(Throughput::Elements(size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), &size, |b, size| {
            let (scheme, poly, mut rng) = sample(*size);
            b.iter_batched(
                || {
                    let coeffs = black_box(poly.to_vec());
                    let commit = scheme.commit(coeffs.clone());
                    let point: Fr = rng.gen();
                    let (eval, poly) = {
                        let poly =
                            ark_poly::univariate::DensePolynomial::<Fr>::from_coefficients_slice(
                                &*coeffs,
                            );
                        let eval = poly.evaluate(&point);
                        (eval, poly.coeffs)
                    };
                    (point, commit, poly, eval)
                },
                |(point, commitment, a, eval)| {
                    let open = scheme.open(commitment, &*a, point, eval);
                    (open, a)
                },
                criterion::BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}
pub fn verify(c: &mut Criterion) {
    let mut group = c.benchmark_group("verify");
    for size in RANGE {
        group.throughput(Throughput::Elements(size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), &size, |b, size| {
            let (scheme, poly, mut rng) = sample(*size);
            b.iter_batched(
                || {
                    let coeffs = black_box(poly.to_vec());
                    let commit = scheme.commit(coeffs.clone());
                    let point: Fr = rng.gen();
                    let (eval, poly) = {
                        let poly =
                            ark_poly::univariate::DensePolynomial::<Fr>::from_coefficients_slice(
                                &*coeffs,
                            );
                        let eval = poly.evaluate(&point);
                        (eval, poly.coeffs)
                    };
                    //(point, commit, poly, eval)
                    let open = scheme.open(commit, &*poly, point, eval);
                    (commit, open)
                },
                |(commit, open)| {
                    let open = scheme.verify(commit, open);
                    open
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
    targets = commit, open, verify
);
criterion_main!(benches);