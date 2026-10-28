use crate::Result;
use crate::model::{Design, ScalarFit, SolverStatus};
use anyhow::{bail, ensure};
use nalgebra::{DMatrix, DVector};
use serde::{Deserialize, Serialize};

pub fn one_group(counts: &[f64], offsets: &[f64], phi: f64) -> Result<ScalarFit> {
    if counts.iter().all(|&value| value <= 1e-10) {
        return Ok(ScalarFit {
            coefficient: f64::NEG_INFINITY,
            status: SolverStatus::AllZero,
        });
    }
    if phi == 0.0 {
        return Ok(ScalarFit {
            coefficient: (counts.iter().sum::<f64>()
                / offsets.iter().map(|x| x.exp()).sum::<f64>())
            .ln(),
            status: SolverStatus::Converged {
                iterations: 0,
                reason: "poisson".into(),
            },
        });
    }
    let mut coefficient = (counts
        .iter()
        .zip(offsets)
        .map(|(&y, &o)| y / o.exp())
        .sum::<f64>()
        / counts.len() as f64)
        .ln();
    for iteration in 1..=50 {
        let mut score = 0.0;
        let mut information = 0.0;
        for (&y, &offset) in counts.iter().zip(offsets) {
            let mean = (coefficient + offset).exp();
            let denominator = 1.0 + mean * phi;
            score += (y - mean) / denominator;
            information += mean / denominator;
        }
        let step = score / information;
        ensure!(step.is_finite(), "nonfinite one-group scoring step");
        coefficient += step;
        if step.abs() < 1e-10 {
            return Ok(ScalarFit {
                coefficient,
                status: SolverStatus::Converged {
                    iterations: iteration,
                    reason: "coefficient-step".into(),
                },
            });
        }
    }
    bail!("one-group solver exceeded 50 iterations")
}

pub fn unit_deviance(y: f64, mean: f64, phi: f64) -> f64 {
    let y = y + 1e-8;
    let mean = mean + 1e-8;
    let residual = y - mean;
    let value = if phi < 1e-4 {
        2.0 * (y * (y / mean).ln() - residual - 0.5 * residual * residual * phi * (1.0 - phi * y))
    } else if mean * phi > 1e6 {
        2.0 * (residual / mean - (y / mean).ln()) * mean / (1.0 + mean * phi)
    } else {
        let size = 1.0 / phi;
        2.0 * (y * (y / mean).ln() + (y + size) * ((mean + size) / (y + size)).ln())
    };
    value.max(0.0)
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct GeneFit {
    pub coefficients: Vec<f64>,
    pub means: Vec<f64>,
    pub deviance: f64,
    pub status: SolverStatus,
}

#[derive(Clone)]
struct Groups {
    indices: Vec<Vec<usize>>,
    representatives: DMatrix<f64>,
}

pub struct PreparedDesign {
    pub design: Design,
    projection: DVector<f64>,
    groups: Option<Groups>,
}

impl PreparedDesign {
    pub fn new(design: &Design) -> Result<Self> {
        let matrix = design.matrix();
        let qr = matrix.clone().qr();
        let q = qr.q();
        let r = qr.r();
        let response = q.transpose() * DVector::from_element(design.n_samples, 1.0);
        let projection = r
            .solve_upper_triangular(&response)
            .ok_or_else(|| anyhow::anyhow!("design QR factor is singular"))?;
        let mut unique: Vec<usize> = Vec::new();
        let mut indices: Vec<Vec<usize>> = Vec::new();
        for sample in 0..design.n_samples {
            if let Some(index) = unique
                .iter()
                .position(|&other| design.row(sample) == design.row(other))
            {
                indices[index].push(sample);
            } else {
                unique.push(sample);
                indices.push(vec![sample]);
            }
        }
        let groups = if unique.len() == design.n_coefs() {
            let values: Vec<f64> = unique
                .iter()
                .flat_map(|&i| design.row(i).iter().copied())
                .collect();
            Some(Groups {
                indices,
                representatives: DMatrix::from_row_slice(
                    design.n_coefs(),
                    design.n_coefs(),
                    &values,
                ),
            })
        } else {
            None
        };
        Ok(Self {
            design: design.clone(),
            projection,
            groups,
        })
    }

    pub fn fit(&self, counts: &[f64], offsets: &[f64], phi: f64) -> Result<GeneFit> {
        if let Some(groups) = &self.groups {
            let mut group_coefficients = Vec::with_capacity(groups.indices.len());
            let mut means = vec![0.0; counts.len()];
            let mut iterations = 0;
            for group in &groups.indices {
                let y: Vec<f64> = group.iter().map(|&i| counts[i]).collect();
                let o: Vec<f64> = group.iter().map(|&i| offsets[i]).collect();
                let fitted = one_group(&y, &o, phi)?;
                if let SolverStatus::Converged {
                    iterations: count, ..
                } = fitted.status
                {
                    iterations = iterations.max(count);
                }
                let coefficient = fitted.coefficient.max(-1e8);
                for &i in group {
                    means[i] = (coefficient + offsets[i]).exp();
                }
                group_coefficients.push(coefficient);
            }
            let coefficients = groups
                .representatives
                .clone()
                .lu()
                .solve(&DVector::from_vec(group_coefficients))
                .ok_or_else(|| anyhow::anyhow!("one-way representative design is singular"))?
                .as_slice()
                .to_vec();
            let deviance = counts
                .iter()
                .zip(&means)
                .map(|(&y, &mu)| unit_deviance(y, mu, phi))
                .sum();
            return Ok(GeneFit {
                coefficients,
                means,
                deviance,
                status: SolverStatus::Converged {
                    iterations,
                    reason: "one-way".into(),
                },
            });
        }
        self.general_fit(counts, offsets, phi)
    }

    fn general_fit(&self, counts: &[f64], offsets: &[f64], phi: f64) -> Result<GeneFit> {
        let p = self.design.n_coefs();
        let maximum = counts.iter().copied().fold(0.0, f64::max);
        ensure!(
            maximum > 1e-10,
            "all-zero general GLM has unidentified coefficients"
        );
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for (&y, &offset) in counts.iter().zip(offsets) {
            let size = offset.exp();
            let weight = size / (1.0 + phi * size);
            numerator += y * weight / size;
            denominator += weight;
        }
        let log_abundance = (numerator / denominator).ln();
        let mut coefficients: Vec<f64> =
            self.projection.iter().map(|&x| x * log_abundance).collect();
        let mut means = fitted_means(&self.design, &coefficients, offsets);
        let mut deviance: f64 = counts
            .iter()
            .zip(&means)
            .map(|(&y, &m)| unit_deviance(y, m, phi))
            .sum();
        let mut max_information: f64 = 0.0;
        let mut damping = 0.0;
        let mut score = vec![0.0; p];
        let mut information = vec![0.0; p * p];
        for iteration in 1..=250 {
            score.fill(0.0);
            information.fill(0.0);
            for sample in 0..counts.len() {
                let mean = means[sample];
                let denominator = 1.0 + mean * phi;
                let weight = mean / denominator;
                let derivative = (counts[sample] - mean) / denominator;
                let x = self.design.row(sample);
                for j in 0..p {
                    score[j] += x[j] * derivative;
                    for k in 0..=j {
                        information[j * p + k] += x[j] * x[k] * weight;
                    }
                }
            }
            for j in 0..p {
                max_information = max_information.max(information[j * p + j]);
                for k in 0..j {
                    information[k * p + j] = information[j * p + k];
                }
            }
            if iteration == 1 {
                damping = (max_information * 1e-6).max(1e-13);
            }
            let mut accepted = None;
            for retry in 0..128 {
                let mut damped = information.clone();
                for j in 0..p {
                    damped[j * p + j] += damping;
                }
                let Some((step, _)) = cholesky_solve(&damped, &score, p) else {
                    damping *= 10.0;
                    continue;
                };
                let next: Vec<f64> = coefficients
                    .iter()
                    .zip(&step)
                    .map(|(&b, &d)| b + d)
                    .collect();
                let next_means = fitted_means(&self.design, &next, offsets);
                let next_deviance: f64 = counts
                    .iter()
                    .zip(&next_means)
                    .map(|(&y, &m)| unit_deviance(y, m, phi))
                    .sum();
                let low_deviance = next_deviance / maximum < 1e-13;
                if next_deviance.is_finite() && (next_deviance <= deviance || low_deviance) {
                    let divergence: f64 = score.iter().zip(&step).map(|(&g, &d)| g * d).sum();
                    coefficients = next;
                    means = next_means;
                    deviance = next_deviance;
                    accepted = Some((divergence, low_deviance, retry));
                    break;
                }
                damping *= 2.0;
                ensure!(
                    damping.is_finite() && damping / max_information <= 1e13,
                    "GLM damping exhausted"
                );
            }
            let (divergence, low_deviance, retry) =
                accepted.ok_or_else(|| anyhow::anyhow!("GLM damping exceeded 128 retries"))?;
            if low_deviance || divergence < 1e-6 {
                return Ok(GeneFit {
                    coefficients,
                    means,
                    deviance,
                    status: SolverStatus::Converged {
                        iterations: iteration,
                        reason: if low_deviance {
                            "low-deviance"
                        } else {
                            "score-step"
                        }
                        .into(),
                    },
                });
            }
            if retry == 0 {
                damping /= 10.0;
            }
        }
        bail!("GLM exceeded 250 iterations")
    }
}

pub fn fitted_means(design: &Design, coefficients: &[f64], offsets: &[f64]) -> Vec<f64> {
    (0..design.n_samples)
        .map(|i| {
            (offsets[i]
                + design
                    .row(i)
                    .iter()
                    .zip(coefficients)
                    .map(|(&x, &b)| x * b)
                    .sum::<f64>())
            .exp()
        })
        .collect()
}

pub fn cholesky_solve(matrix: &[f64], rhs: &[f64], p: usize) -> Option<(Vec<f64>, f64)> {
    let mut lower = vec![0.0; p * p];
    for i in 0..p {
        for j in 0..=i {
            let mut value = matrix[i * p + j];
            for k in 0..j {
                value -= lower[i * p + k] * lower[j * p + k];
            }
            if i == j {
                if !value.is_finite() || value <= 0.0 {
                    return None;
                }
                lower[i * p + j] = value.sqrt();
            } else {
                lower[i * p + j] = value / lower[j * p + j];
            }
        }
    }
    let mut solution = rhs.to_vec();
    for i in 0..p {
        for j in 0..i {
            solution[i] -= lower[i * p + j] * solution[j];
        }
        solution[i] /= lower[i * p + i];
    }
    for i in (0..p).rev() {
        for j in i + 1..p {
            solution[i] -= lower[j * p + i] * solution[j];
        }
        solution[i] /= lower[i * p + i];
    }
    let half_log_determinant = (0..p).map(|i| lower[i * p + i].ln()).sum();
    Some((solution, half_log_determinant))
}

pub fn leverages(design: &Design, weights: &[f64]) -> Vec<f64> {
    let n = design.n_samples;
    let mut basis: Vec<Vec<f64>> = Vec::new();
    for j in 0..design.n_coefs() {
        let mut column: Vec<f64> = (0..n)
            .map(|i| design.row(i)[j] * weights[i].sqrt())
            .collect();
        let initial_norm = column.iter().map(|x| x * x).sum::<f64>().sqrt();
        for _ in 0..2 {
            for previous in &basis {
                let inner: f64 = column.iter().zip(previous).map(|(&a, &b)| a * b).sum();
                for (value, &direction) in column.iter_mut().zip(previous) {
                    *value -= inner * direction;
                }
            }
        }
        let norm = column.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm > 1e-7 * initial_norm && norm > 0.0 {
            for value in &mut column {
                *value /= norm;
            }
            basis.push(column);
        }
    }
    (0..n)
        .map(|i| basis.iter().map(|column| column[i] * column[i]).sum())
        .collect()
}
