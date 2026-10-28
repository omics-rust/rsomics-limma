use crate::{
    Result,
    number::{Number, RowValues},
    numerics::{brent, log_minus_digamma, quantile, trigamma},
};
use anyhow::{bail, ensure};

pub struct Moderated {
    pub prior_df: Number,
    pub prior_variance: RowValues,
    pub posterior_variance: Vec<f64>,
    pub method: String,
}

pub fn moderate(sigma: &[Number], df: &[usize]) -> Result<Moderated> {
    let n = sigma.len();
    let variance: Vec<f64> = sigma.iter().map(|x| x.value().powi(2)).collect();
    ensure!(
        df.iter().any(|&x| x > 0),
        "no residual degrees of freedom for empirical Bayes"
    );
    if n < 3 {
        ensure!(
            variance.iter().all(|x| x.is_finite()),
            "nonfinite unmoderated variance in a tiny fit"
        );
        return Ok(Moderated {
            prior_df: Number::Finite(0.0),
            prior_variance: RowValues::Rows(variance.clone()),
            posterior_variance: variance,
            method: "no-shrink-small-population".into(),
        });
    }
    let positive_df: Vec<usize> = df.iter().copied().filter(|&x| x > 0).collect();
    let legacy = positive_df.iter().min() == positive_df.iter().max();
    let raw: Vec<f64> = variance
        .iter()
        .zip(df)
        .map(|(&x, &d)| if d == 0 { 0.0 } else { x })
        .collect();
    let (scale, prior_df) = if legacy {
        legacy_prior(&raw, df)?
    } else {
        unequal_prior(&raw, df)?
    };
    let posterior_variance: Vec<f64> = raw
        .iter()
        .zip(df)
        .map(|(&x, &d)| {
            if prior_df.is_infinite() {
                scale
            } else {
                (d as f64 * x + prior_df * scale) / (d as f64 + prior_df)
            }
        })
        .collect();
    ensure!(
        posterior_variance.iter().all(|x| x.is_finite() && *x > 0.0),
        "nonfinite or nonpositive posterior variance"
    );
    Ok(Moderated {
        prior_df: Number::from_value(prior_df),
        prior_variance: RowValues::Scalar(scale),
        posterior_variance,
        method: if legacy {
            "legacy-equal-positive-df"
        } else {
            "unequal-df"
        }
        .into(),
    })
}

fn unequal_prior(raw: &[f64], df: &[usize]) -> Result<(f64, f64)> {
    let mut x = raw.to_vec();
    let mut d = Vec::with_capacity(df.len());
    let mut pw = Vec::with_capacity(df.len());
    for i in 0..x.len() {
        let informative = x[i].is_finite() && df[i] > 0;
        if !x[i].is_finite() {
            x[i] = 0.0;
        }
        ensure!(x[i] >= 0.0, "negative variance in unequal-df moderation");
        pw.push(if informative { 1.0 } else { 0.0 });
        d.push(if df[i] > 0 { df[i] as f64 / 2.0 } else { 0.5 });
    }
    let informative: Vec<f64> = x
        .iter()
        .zip(&pw)
        .filter_map(|(&x, &w)| (x > 0.0 && w > 0.0).then_some(x))
        .collect();
    ensure!(
        informative.len() >= 2,
        "could not estimate variance prior: fewer than two informative variances"
    );
    if informative.len() == 2 && pw.contains(&0.0) {
        bail!("weighted two-informative variance source boundary is not qualified");
    }
    let floor = 1e-12 * quantile(&informative, 0.5);
    for value in &mut x {
        *value = value.max(floor);
    }
    let weights: Vec<f64> = d
        .iter()
        .zip(&pw)
        .map(|(&d, &pw)| pw / trigamma(d))
        .collect();
    let emean = x
        .iter()
        .zip(&d)
        .zip(&weights)
        .map(|((&x, &d), &w)| w * (x.ln() + log_minus_digamma(d)))
        .sum::<f64>()
        / weights.iter().sum::<f64>();
    let theta = brent(
        |theta| {
            let delta = theta / (1.0 - theta);
            let denominator = delta * (emean - log_minus_digamma(delta)).exp();
            Ok(-2.0
                * x.iter()
                    .zip(&d)
                    .zip(&pw)
                    .filter(|(_, w)| **w > 0.0)
                    .map(|((&x, &d), &pw)| {
                        pw * (-(d + delta) * (d * x / denominator).ln_1p() - d * denominator.ln()
                            + libm::lgamma(d + delta)
                            - libm::lgamma(delta))
                    })
                    .sum::<f64>())
        },
        0.5,
        0.9998,
        2.0_f64.powi(-13),
    )?;
    let delta = theta / (1.0 - theta);
    Ok(((emean - log_minus_digamma(delta)).exp(), 2.0 * delta))
}

fn legacy_prior(raw: &[f64], df: &[usize]) -> Result<(f64, f64)> {
    let mut values: Vec<f64> = raw
        .iter()
        .zip(df)
        .filter_map(|(&x, &d)| (d > 0 && x.is_finite() && x > -1e-15).then_some(x.max(0.0)))
        .collect();
    ensure!(!values.is_empty(), "no usable variances for equal-df prior");
    if values.len() == 1 {
        return Ok((values[0], 0.0));
    }
    let median = quantile(&values, 0.5);
    let floor = 1e-5 * if median == 0.0 { 1.0 } else { median };
    for x in &mut values {
        *x = x.max(floor);
    }
    let half = *df
        .iter()
        .find(|&&d| d > 0)
        .ok_or_else(|| anyhow::anyhow!("no positive residual df"))? as f64
        / 2.0;
    let correction = log_minus_digamma(half);
    let e: Vec<f64> = values.iter().map(|x| x.ln() + correction).collect();
    let mean = e.iter().sum::<f64>() / e.len() as f64;
    let spread =
        e.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (e.len() - 1) as f64 - trigamma(half);
    if spread <= 0.0 {
        return Ok((
            values.iter().sum::<f64>() / values.len() as f64,
            f64::INFINITY,
        ));
    }
    let mut lower = 0.0;
    let mut upper = (1.0 + 2.0 / spread).max(1.0);
    for _ in 0..128 {
        if trigamma(upper) <= spread {
            break;
        }
        upper *= 2.0;
    }
    ensure!(
        upper.is_finite() && trigamma(upper) <= spread,
        "inverse trigamma bracket exhausted"
    );
    for _ in 0..100 {
        let middle = (lower + upper) / 2.0;
        if trigamma(middle) > spread {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    let delta = (lower + upper) / 2.0;
    Ok(((mean - log_minus_digamma(delta)).exp(), 2.0 * delta))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn automatic_equal_df_prior_and_nonfinite_codec() {
        let moderated = moderate(&[Number::Finite(2.0); 11], &[4; 11]).unwrap();
        assert_eq!(moderated.method, "legacy-equal-positive-df");
        assert_eq!(moderated.prior_df, Number::PosInf);
        assert_eq!(moderated.prior_variance, RowValues::Scalar(4.0));
        assert_eq!(moderated.posterior_variance, vec![4.0; 11]);
        let values = [Number::Na, Number::Nan, Number::PosInf, Number::NegInf];
        let encoded = serde_json::to_vec(&values).unwrap();
        assert_eq!(
            serde_json::from_slice::<Vec<Number>>(&encoded).unwrap(),
            values
        );
    }
}
