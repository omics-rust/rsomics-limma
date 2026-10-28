use crate::Result;
use crate::input::orthonormal_columns;
use crate::model::{Dataset, FilterThresholds, LibraryState};
use anyhow::ensure;
use nalgebra::DMatrix;

pub fn prepare(dataset: &Dataset) -> Result<LibraryState> {
    let counts = &dataset.counts;
    let n = counts.n_samples();
    let mut raw_library_sizes = vec![0.0; n];
    for row in counts.values.chunks_exact(n) {
        for (library, count) in raw_library_sizes.iter_mut().zip(row) {
            *library += count;
        }
    }
    positive_libraries(&raw_library_sizes, dataset, "raw")?;
    let design = DMatrix::from_fn(n, dataset.design.n_coefs() + 1, |sample, column| {
        if column == 0 {
            1.0
        } else {
            dataset.design.row(sample)[column - 1]
        }
    });
    let basis = orthonormal_columns(&design)?;
    let max_hat = (0..n)
        .map(|sample| {
            basis
                .iter()
                .map(|column| column[sample].powi(2))
                .sum::<f64>()
        })
        .fold(0.0, f64::max);
    let mut min_sample_size = 1.0 / max_hat;
    if min_sample_size > 10.0 {
        min_sample_size = 10.0 + (min_sample_size - 10.0) * 0.7;
    }
    let cpm_cutoff = 10.0 / quantile(&raw_library_sizes, 0.5) * 1e6;
    let keep: Vec<_> = counts
        .values
        .chunks_exact(n)
        .map(|row| {
            let above = row
                .iter()
                .zip(&raw_library_sizes)
                .filter(|(count, library)| **count * 1e6 / **library >= cpm_cutoff)
                .count();
            above as f64 >= min_sample_size - 1e-14 && row.iter().sum::<f64>() >= 15.0 - 1e-14
        })
        .collect();
    let kept_indices: Vec<_> = keep
        .iter()
        .enumerate()
        .filter_map(|(index, &keep)| keep.then_some(index))
        .collect();
    ensure!(
        !kept_indices.is_empty(),
        "expression filter retained no genes"
    );
    let mut postfilter_library_sizes = vec![0.0; n];
    for &gene in &kept_indices {
        for (library, count) in postfilter_library_sizes.iter_mut().zip(counts.row(gene)) {
            *library += count;
        }
    }
    positive_libraries(&postfilter_library_sizes, dataset, "post-filter")?;
    let (tmm_reference_sample, tmm_factors) =
        tmm(dataset, &kept_indices, &postfilter_library_sizes)?;
    let effective_library_sizes: Vec<_> = postfilter_library_sizes
        .iter()
        .zip(&tmm_factors)
        .map(|(library, factor)| library * factor)
        .collect();
    positive_libraries(&effective_library_sizes, dataset, "effective")?;
    let offsets: Vec<_> = effective_library_sizes
        .iter()
        .map(|library| library.ln())
        .collect();
    Ok(LibraryState {
        keep,
        kept_indices,
        raw_library_sizes,
        postfilter_library_sizes,
        tmm_factors,
        effective_library_sizes,
        offsets,
        filter: FilterThresholds {
            cpm_cutoff,
            min_sample_size,
            min_total_count: 15.0,
        },
        tmm_reference_sample,
    })
}

fn positive_libraries(libraries: &[f64], dataset: &Dataset, label: &str) -> Result<()> {
    for (sample, &library) in dataset.counts.sample_ids.iter().zip(libraries) {
        ensure!(
            library.is_finite() && library > 0.0,
            "sample '{sample}' has invalid {label} library size {library}"
        );
    }
    ensure!(
        libraries.iter().sum::<f64>().is_finite(),
        "{label} library total overflowed"
    );
    Ok(())
}

fn tmm(dataset: &Dataset, kept: &[usize], libraries: &[f64]) -> Result<(usize, Vec<f64>)> {
    let n = libraries.len();
    let quartiles: Vec<_> = (0..n)
        .map(|sample| {
            let values: Vec<_> = kept
                .iter()
                .map(|&gene| dataset.counts.row(gene)[sample])
                .collect();
            quantile(&values, 0.75) / libraries[sample]
        })
        .collect();
    let reference = if quantile(&quartiles, 0.5) < 1e-20 {
        let masses: Vec<_> = (0..n)
            .map(|sample| {
                kept.iter()
                    .map(|&gene| dataset.counts.row(gene)[sample].sqrt())
                    .sum::<f64>()
            })
            .collect();
        (1..n).fold(0, |best, sample| {
            if masses[sample] > masses[best] {
                sample
            } else {
                best
            }
        })
    } else {
        let mean = quartiles.iter().sum::<f64>() / n as f64;
        (1..n).fold(0, |best, sample| {
            if (quartiles[sample] - mean).abs() < (quartiles[best] - mean).abs() {
                sample
            } else {
                best
            }
        })
    };
    let reference_values: Vec<_> = kept
        .iter()
        .map(|&gene| dataset.counts.row(gene)[reference])
        .collect();
    let mut observed = vec![0.0; kept.len()];
    let mut factors = Vec::with_capacity(n);
    for sample in 0..n {
        for (value, &gene) in observed.iter_mut().zip(kept) {
            *value = dataset.counts.row(gene)[sample];
        }
        factors.push(tmm_pair(
            &observed,
            &reference_values,
            libraries[sample],
            libraries[reference],
        ));
    }
    let center = (factors.iter().map(|factor| factor.ln()).sum::<f64>() / n as f64).exp();
    for factor in &mut factors {
        *factor /= center;
        ensure!(
            factor.is_finite() && *factor > 0.0,
            "TMM produced an invalid factor"
        );
    }
    Ok((reference, factors))
}

fn tmm_pair(
    observed: &[f64],
    reference: &[f64],
    observed_library: f64,
    reference_library: f64,
) -> f64 {
    let mut ratios = Vec::with_capacity(observed.len());
    let mut abundance = Vec::with_capacity(observed.len());
    let mut variance = Vec::with_capacity(observed.len());
    for (&observed, &reference) in observed.iter().zip(reference) {
        let observed_fraction = observed / observed_library;
        let reference_fraction = reference / reference_library;
        let ratio = (observed_fraction / reference_fraction).log2();
        let absolute = (observed_fraction.log2() + reference_fraction.log2()) / 2.0;
        if ratio.is_finite() && absolute.is_finite() && absolute > -1e10 {
            ratios.push(ratio);
            abundance.push(absolute);
            variance.push(
                (observed_library - observed) / observed_library / observed
                    + (reference_library - reference) / reference_library / reference,
            );
        }
    }
    if ratios.is_empty() || ratios.iter().all(|ratio| ratio.abs() < 1e-6) {
        return 1.0;
    }
    let ratio_ranks = average_rank(&ratios);
    let abundance_ranks = average_rank(&abundance);
    let n = ratios.len() as f64;
    let lower_ratio = (n * 0.3).floor() + 1.0;
    let upper_ratio = n + 1.0 - lower_ratio;
    let lower_abundance = (n * 0.05).floor() + 1.0;
    let upper_abundance = n + 1.0 - lower_abundance;
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for index in 0..ratios.len() {
        if ratio_ranks[index] >= lower_ratio
            && ratio_ranks[index] <= upper_ratio
            && abundance_ranks[index] >= lower_abundance
            && abundance_ranks[index] <= upper_abundance
        {
            let value = ratios[index] / variance[index];
            let weight = 1.0 / variance[index];
            if !value.is_nan() {
                numerator += value;
            }
            if !weight.is_nan() {
                denominator += weight;
            }
        }
    }
    let mean = numerator / denominator;
    2.0_f64.powf(if mean.is_nan() { 0.0 } else { mean })
}

fn average_rank(values: &[f64]) -> Vec<f64> {
    let mut order: Vec<_> = (0..values.len()).collect();
    order.sort_by(|&left, &right| values[left].total_cmp(&values[right]));
    let mut ranks = vec![0.0; values.len()];
    let mut first = 0;
    while first < order.len() {
        let mut end = first + 1;
        while end < order.len() && values[order[first]] == values[order[end]] {
            end += 1;
        }
        let rank = (first + 1 + end) as f64 / 2.0;
        for &index in &order[first..end] {
            ranks[index] = rank;
        }
        first = end;
    }
    ranks
}

fn quantile(values: &[f64], probability: f64) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let position = (sorted.len() - 1) as f64 * probability;
    let lower = position.floor() as usize;
    if lower + 1 == sorted.len() {
        sorted[lower]
    } else {
        sorted[lower] + (position - lower as f64) * (sorted[lower + 1] - sorted[lower])
    }
}
