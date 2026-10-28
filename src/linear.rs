use crate::Result;
use crate::model::Design;
use crate::number::Number;
use anyhow::ensure;
use nalgebra::{DMatrix, DVector};

pub struct System {
    observations: Vec<usize>,
    square_weights: Vec<f64>,
    design: DMatrix<f64>,
    q: DMatrix<f64>,
    r: DMatrix<f64>,
    pub rank: usize,
    pub pivot: Vec<usize>,
    pub covariance: Vec<f64>,
    pub unscaled_se: Vec<f64>,
}

pub struct RowFit {
    pub coefficients: Vec<f64>,
    pub sigma: Number,
    pub residual_df: usize,
}

impl System {
    pub fn new(design: &Design, weights: Option<&[f64]>, omit: Option<&[bool]>) -> Result<Self> {
        let p = design.n_coefs();
        let observations: Vec<usize> = (0..design.n_samples)
            .filter(|&i| !omit.is_some_and(|mask| mask[i]))
            .collect();
        let n = observations.len();
        let square_weights: Vec<f64> = observations
            .iter()
            .map(|&i| weights.map_or(1.0, |values| values[i]).sqrt())
            .collect();
        ensure!(
            square_weights.iter().all(|x| x.is_finite() && *x > 0.0),
            "invalid linear-fit weight"
        );
        let x = DMatrix::from_fn(n, p, |i, j| {
            design.row(observations[i])[j] * square_weights[i]
        });
        let mut directions: Vec<DVector<f64>> = Vec::new();
        let mut pivot = Vec::new();
        let mut aliased = Vec::new();
        for j in 0..p {
            let column = x.column(j);
            let original = column.norm();
            let mut residual = column.into_owned();
            for _ in 0..2 {
                for direction in &directions {
                    let amount = residual.dot(direction);
                    residual -= direction * amount;
                }
            }
            let norm = residual.norm();
            if original == 0.0 || norm < 1e-7 * original || directions.len() == n {
                aliased.push(j);
            } else {
                directions.push(residual / norm);
                pivot.push(j);
            }
        }
        let rank = pivot.len();
        let selected = DMatrix::from_fn(n, rank, |i, j| x[(i, pivot[j])]);
        let (q, r) = if rank == 0 {
            (DMatrix::zeros(n, 0), DMatrix::zeros(0, 0))
        } else {
            let qr = selected.clone().qr();
            (qr.q(), qr.r())
        };
        let mut covariance = vec![f64::NAN; p * p];
        let mut unscaled_se = vec![f64::NAN; p];
        if rank > 0 {
            let inverse = r
                .clone()
                .try_inverse()
                .ok_or_else(|| anyhow::anyhow!("linear QR triangular factor is singular"))?;
            let cov = &inverse * inverse.transpose();
            for i in 0..rank {
                unscaled_se[pivot[i]] = cov[(i, i)].sqrt();
                for j in 0..rank {
                    covariance[pivot[i] * p + pivot[j]] = cov[(i, j)];
                }
            }
        }
        pivot.extend(aliased);
        Ok(Self {
            observations,
            square_weights,
            design: selected,
            q,
            r,
            rank,
            pivot,
            covariance,
            unscaled_se,
        })
    }

    pub fn fit(&self, response: &[f64]) -> Result<RowFit> {
        let y = DVector::from_iterator(
            self.observations.len(),
            self.observations
                .iter()
                .zip(&self.square_weights)
                .map(|(&i, &w)| response[i] * w),
        );
        let beta = if self.rank == 0 {
            DVector::zeros(0)
        } else {
            self.r
                .solve_upper_triangular(&(self.q.transpose() * &y))
                .ok_or_else(|| anyhow::anyhow!("linear QR solve failed"))?
        };
        let residual = y - &self.design * &beta;
        let residual_df = self.observations.len() - self.rank;
        let sigma = if residual_df == 0 {
            Number::Na
        } else {
            Number::Finite((residual.norm_squared() / residual_df as f64).sqrt())
        };
        let mut coefficients = vec![f64::NAN; self.pivot.len()];
        for j in 0..self.rank {
            coefficients[self.pivot[j]] = beta[j];
        }
        ensure!(
            beta.iter().all(|x| x.is_finite()),
            "linear coefficient overflow"
        );
        Ok(RowFit {
            coefficients,
            sigma,
            residual_df,
        })
    }
}
