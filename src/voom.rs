use crate::{
    Result,
    linear::System,
    model::{Dataset, LibraryState, LinearFit, SolverStatus},
    poisson::PreparedDesign,
    smoothing::{self, Curve},
};
use anyhow::{Context, ensure};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Target {
    pub sample: String,
    pub group: String,
    pub lib_size: f64,
    pub norm_factor: f64,
    pub metadata: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BaseFit {
    #[serde(flatten)]
    pub linear: LinearFit,
    pub targets: Vec<Target>,
    pub weights: Vec<Vec<f64>>,
    pub zero_mask: Vec<Vec<bool>>,
    pub zero_probes: Vec<Option<SolverStatus>>,
    pub derived_span: f64,
    pub voom_points: Vec<[f64; 2]>,
    pub voom_line: Curve,
    pub trend_policy: String,
}

pub fn fit(data: &Dataset, lib: &LibraryState) -> Result<BaseFit> {
    let n = data.counts.n_samples();
    let p = data.design.n_coefs();
    let m = lib.kept_indices.len();
    ensure!(m >= 2, "voom needs at least two retained genes");
    let expression: Vec<Vec<f64>> = lib
        .kept_indices
        .par_iter()
        .map(|&gene| {
            data.counts
                .row(gene)
                .iter()
                .zip(&lib.effective_library_sizes)
                .map(|(&y, &size)| ((y + 0.5) / (size + 1.0) * 1e6).log2())
                .collect()
        })
        .collect();
    let unweighted = System::new(&data.design, None, None)?;
    let initial = expression
        .par_iter()
        .map(|y| unweighted.fit(y))
        .collect::<Result<Vec<_>>>()?;
    let amean: Vec<f64> = expression
        .iter()
        .map(|y| y.iter().sum::<f64>() / n as f64)
        .collect();
    let leverage = crate::poisson::leverages(&data.design, &vec![1.0; n]);
    let min_group = 1.0 / leverage.into_iter().fold(0.0, f64::max);
    let poisson = PreparedDesign::new(&data.design)?;
    let probes: Vec<(Vec<bool>, Option<SolverStatus>)> = lib
        .kept_indices
        .par_iter()
        .map(|&gene| {
            let counts = data.counts.row(gene);
            if counts.iter().filter(|&&y| y < 1e-4).count() as f64 <= min_group.max(2.0) - 1e-4 {
                return Ok((vec![false; n], None));
            }
            let fitted = poisson.fit(counts, &lib.offsets, 0.0).with_context(|| {
                format!("zero-loss Poisson fit for {}", data.counts.gene_ids[gene])
            })?;
            Ok((
                counts
                    .iter()
                    .zip(&fitted.means)
                    .map(|(&y, &mu)| y < 1e-4 && mu < 1e-4)
                    .collect(),
                Some(fitted.status),
            ))
        })
        .collect::<Result<_>>()?;
    let any_zero = probes.iter().any(|(mask, _)| mask.iter().any(|x| *x));
    let corrected = expression
        .par_iter()
        .enumerate()
        .map(|(i, y)| {
            if probes[i].0.iter().any(|x| *x) {
                System::new(&data.design, None, Some(&probes[i].0))?.fit(y)
            } else {
                Ok(crate::linear::RowFit {
                    coefficients: Vec::new(),
                    sigma: initial[i].sigma,
                    residual_df: initial[i].residual_df,
                })
            }
        })
        .collect::<Result<Vec<_>>>()?;
    let mean_log_library = lib
        .effective_library_sizes
        .iter()
        .map(|size| (size + 1.0).log2())
        .sum::<f64>()
        / n as f64;
    let mut sx = Vec::new();
    let mut sy = Vec::new();
    let mut df_weights = Vec::new();
    for i in 0..m {
        if corrected[i].residual_df == 0 {
            continue;
        }
        let mean = if probes[i].0.iter().any(|x| *x) {
            let included: Vec<f64> = expression[i]
                .iter()
                .zip(&probes[i].0)
                .filter_map(|(&y, &omit)| (!omit).then_some(y))
                .collect();
            included.iter().sum::<f64>() / included.len() as f64
        } else {
            amean[i]
        };
        sx.push(mean + mean_log_library - 1e6_f64.log2());
        sy.push(corrected[i].sigma.value().sqrt());
        df_weights.push(corrected[i].residual_df as f64);
    }
    ensure!(
        sx.len() >= 2,
        "fewer than two genes have replication for the complete voom workflow"
    );
    let derived_span = (0.3 + 0.7 * (50.0 / m as f64).powf(1.0 / 3.0)).min(1.0);
    let voom_line = smoothing::curve(
        &sx,
        &sy,
        any_zero.then_some(df_weights.as_slice()),
        derived_span,
    )?;
    let weights: Vec<Vec<f64>> = initial
        .par_iter()
        .map(|row| {
            (0..n)
                .map(|sample| {
                    let value: f64 = data
                        .design
                        .row(sample)
                        .iter()
                        .zip(&row.coefficients)
                        .map(|(&x, &beta)| x * beta)
                        .sum();
                    let count = 1e-6 * value.exp2() * (lib.effective_library_sizes[sample] + 1.0);
                    let trend = voom_line.evaluate(count.log2());
                    ensure!(
                        trend.is_finite() && trend > 0.0,
                        "nonpositive fitted voom trend"
                    );
                    Ok(1.0 / trend.powi(4))
                })
                .collect::<Result<_>>()
        })
        .collect::<Result<_>>()?;
    let fitted = expression
        .par_iter()
        .zip(&weights)
        .enumerate()
        .map(|(i, (y, w))| {
            let system = System::new(&data.design, Some(w), None)?;
            ensure!(system.rank == p, "voom weighted design lost rank");
            let mut row = system.fit(y)?;
            if probes[i].0.iter().any(|x| *x) {
                let correction = System::new(&data.design, Some(w), Some(&probes[i].0))?.fit(y)?;
                row.sigma = correction.sigma;
                row.residual_df = correction.residual_df;
            }
            Ok((row, system.unscaled_se))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(BaseFit {
        linear: LinearFit {
            design: data.design.clone(),
            coefficients: fitted
                .iter()
                .map(|(row, _)| row.coefficients.clone())
                .collect(),
            unscaled_se: fitted.iter().map(|(_, se)| se.clone()).collect(),
            sigma: fitted.iter().map(|(row, _)| row.sigma).collect(),
            residual_df: fitted.iter().map(|(row, _)| row.residual_df).collect(),
            covariance: unweighted
                .covariance
                .chunks_exact(p)
                .map(|row| row.to_vec())
                .collect(),
            rank: unweighted.rank,
            pivot: unweighted.pivot.iter().map(|j| j + 1).collect(),
            method: "ls".into(),
            amean,
            expression,
        },
        targets: (0..n)
            .map(|i| Target {
                sample: data.counts.sample_ids[i].clone(),
                group: data
                    .metadata_columns
                    .iter()
                    .position(|x| x == "group")
                    .map_or_else(|| "1".into(), |j| data.metadata[i][j].clone()),
                lib_size: lib.postfilter_library_sizes[i],
                norm_factor: lib.tmm_factors[i],
                metadata: data.metadata[i].clone(),
            })
            .collect(),
        weights,
        zero_mask: probes.iter().map(|(mask, _)| mask.clone()).collect(),
        zero_probes: probes.into_iter().map(|(_, status)| status).collect(),
        derived_span,
        voom_points: sx.into_iter().zip(sy).map(|(x, y)| [x, y]).collect(),
        voom_line,
        trend_policy: if any_zero {
            "weighted-lowess/frequency-df/four-passes"
        } else {
            "ordinary-lowess/three-robustness-updates"
        }
        .into(),
    })
}
