use crate::{
    Result, ebayes,
    number::{Number, RowValues},
    voom::BaseFit,
};
use anyhow::{bail, ensure};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ResultRow {
    pub log_fc: f64,
    pub ave_expr: f64,
    pub unscaled_se: f64,
    pub posterior_variance: f64,
    pub t: f64,
    pub p_value: f64,
    pub fdr: f64,
    pub lods: f64,
    pub df_total: f64,
    pub f_statistic: f64,
    pub f_p_value: f64,
    pub f_denominator_df: Number,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ContrastResult {
    pub contrast: Vec<f64>,
    pub covariance: Vec<Vec<f64>>,
    pub prior_df: Number,
    pub prior_variance: RowValues,
    pub coefficient_prior_variance: f64,
    pub coefficient_prior_fallback: bool,
    pub proportion: f64,
    pub moderation_method: String,
    pub base_rank: usize,
    pub base_pivot: Vec<usize>,
    pub rows: Vec<ResultRow>,
    pub full_base_reused: bool,
    pub linear_model_refit: bool,
}

pub fn derive(base: &BaseFit, contrast: &[f64]) -> Result<ContrastResult> {
    let p = base.design.n_coefs();
    let n = base.coefficients.len();
    ensure!(
        contrast.len() == p
            && contrast.iter().all(|x| x.is_finite())
            && contrast.iter().any(|&x| x != 0.0),
        "invalid complete named contrast"
    );
    let mut covariance = 0.0;
    for j in 0..p {
        for k in 0..p {
            covariance += contrast[j] * base.covariance[j][k] * contrast[k];
        }
    }
    ensure!(
        covariance.is_finite() && covariance > 0.0,
        "contrast covariance is not positive"
    );
    let moderation = ebayes::moderate(&base.sigma, &base.residual_df)?;
    let prior_df = moderation.prior_df.value();
    let pooled = base.residual_df.iter().sum::<usize>() as f64;
    let mut rows = Vec::with_capacity(n);
    for (i, ((beta, se), &post)) in base
        .coefficients
        .iter()
        .zip(&base.unscaled_se)
        .zip(&moderation.posterior_variance)
        .enumerate()
    {
        let coefficient = beta.iter().zip(contrast).map(|(&b, &c)| b * c).sum::<f64>();
        let mut variance = 0.0;
        for j in 0..p {
            for k in 0..p {
                let correlation =
                    base.covariance[j][k] / (base.covariance[j][j] * base.covariance[k][k]).sqrt();
                variance += se[j] * contrast[j] * correlation * se[k] * contrast[k];
            }
        }
        let unscaled = variance.sqrt();
        let statistic = coefficient / unscaled / post.sqrt();
        let uncapped = prior_df + base.residual_df[i] as f64;
        let df_total = uncapped.min(pooled);
        let p_value = two_sided(statistic, df_total)?;
        let f_statistic = statistic * statistic;
        let f_p_value = if uncapped.is_infinite() {
            rsomics_stats::chi2_sf(f_statistic, 1.0)
        } else {
            two_sided(statistic, uncapped)?
        };
        ensure!(
            statistic.is_finite() && f_p_value.is_finite(),
            "nonfinite moderated statistic"
        );
        rows.push(ResultRow {
            log_fc: coefficient,
            ave_expr: base.amean[i],
            unscaled_se: unscaled,
            posterior_variance: post,
            t: statistic,
            p_value,
            fdr: 0.0,
            lods: 0.0,
            df_total,
            f_statistic,
            f_p_value,
            f_denominator_df: Number::from_value(uncapped),
        });
    }
    let adjusted =
        rsomics_stats::bh_adjust(&rows.iter().map(|row| row.p_value).collect::<Vec<_>>())?;
    let prior_scale = match &moderation.prior_variance {
        RowValues::Scalar(x) => *x,
        RowValues::Rows(x) => crate::numerics::quantile(x, 0.5),
    };
    let limits = (0.1_f64.powi(2) / prior_scale, 4.0_f64.powi(2) / prior_scale);
    let target = (0.01 * n as f64 / 2.0).ceil() as usize;
    let proportion = 0.01_f64.max(target as f64 / n as f64);
    let max_df = rows.iter().map(|row| row.df_total).fold(0.0, f64::max);
    let mut order = rows
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mapped = if row.df_total < max_df {
                upper_quantile(log_upper(row.t.abs(), row.df_total)?, max_df)?
            } else {
                row.t.abs()
            };
            Ok((i, mapped))
        })
        .collect::<Result<Vec<_>>>()?;
    order.sort_by(|a, b| b.1.total_cmp(&a.1));
    let mut variance_sum = 0.0;
    for (position, &(i, t)) in order.iter().take(target).enumerate() {
        let p0 = two_sided(t, max_df)?;
        let ptarget = (((position as f64 + 0.5) / n as f64) - (1.0 - proportion) * p0) / proportion;
        let v = if ptarget > p0 {
            let q = upper_quantile((ptarget / 2.0).ln(), max_df)?;
            rows[i].unscaled_se.powi(2) * ((t / q).powi(2) - 1.0)
        } else {
            0.0
        };
        variance_sum += v.clamp(limits.0, limits.1);
    }
    let estimated = variance_sum / target as f64;
    let fallback = !estimated.is_finite();
    let coefficient_prior_variance = if fallback {
        1.0 / prior_scale
    } else {
        estimated
    };
    for (row, adjusted) in rows.iter_mut().zip(adjusted) {
        let r = (row.unscaled_se.powi(2) + coefficient_prior_variance) / row.unscaled_se.powi(2);
        let t2 = row.t * row.t;
        let kernel = if prior_df > 1e6 {
            t2 * (1.0 - 1.0 / r) / 2.0
        } else {
            (1.0 + row.df_total) / 2.0 * ((t2 + row.df_total) / (t2 / r + row.df_total)).ln()
        };
        row.lods = (0.01_f64 / 0.99).ln() - r.ln() / 2.0 + kernel;
        row.fdr = adjusted;
        ensure!(row.lods.is_finite(), "nonfinite posterior log odds");
    }
    Ok(ContrastResult {
        contrast: contrast.to_vec(),
        covariance: vec![vec![covariance]],
        prior_df: moderation.prior_df,
        prior_variance: moderation.prior_variance,
        coefficient_prior_variance,
        coefficient_prior_fallback: fallback,
        proportion: 0.01,
        moderation_method: moderation.method,
        base_rank: base.rank,
        base_pivot: base.pivot.clone(),
        rows,
        full_base_reused: true,
        linear_model_refit: false,
    })
}

fn log_upper(t: f64, df: f64) -> Result<f64> {
    ensure!(
        t >= 0.0 && t.is_finite() && df > 0.0 && df.is_finite(),
        "invalid Student-t tail arguments"
    );
    let probability =
        statrs::function::beta::checked_beta_reg(df / 2.0, 0.5, 1.0 / (1.0 + t * t / df))?;
    Ok(probability.ln() - std::f64::consts::LN_2)
}

fn two_sided(t: f64, df: f64) -> Result<f64> {
    Ok((log_upper(t.abs(), df)? + std::f64::consts::LN_2)
        .exp()
        .min(1.0))
}

fn upper_quantile(log_probability: f64, df: f64) -> Result<f64> {
    ensure!(
        log_probability.is_finite() && log_probability <= -std::f64::consts::LN_2,
        "invalid Student-t upper-tail probability"
    );
    if log_probability == -std::f64::consts::LN_2 {
        return Ok(0.0);
    }
    let mut upper = 1.0;
    let mut bracketed = false;
    for _ in 0..128 {
        if log_upper(upper, df)? <= log_probability {
            bracketed = true;
            break;
        }
        upper *= 2.0;
    }
    if !bracketed {
        bail!("Student-t inverse-tail bracket exhausted");
    }
    let mut lower = 0.0;
    for _ in 0..100 {
        let middle = (lower + upper) / 2.0;
        if log_upper(middle, df)? > log_probability {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    Ok((lower + upper) / 2.0)
}
