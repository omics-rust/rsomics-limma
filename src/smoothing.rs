use crate::Result;
use anyhow::ensure;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Curve {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub delta: f64,
}

impl Curve {
    fn tied_mean(&self, index: usize) -> f64 {
        let point = self.x[index];
        let mut left = index;
        let mut right = index + 1;
        while left > 0 && self.x[left - 1] == point {
            left -= 1;
        }
        while right < self.x.len() && self.x[right] == point {
            right += 1;
        }
        self.y[left..right].iter().sum::<f64>() / (right - left) as f64
    }
    pub fn evaluate(&self, point: f64) -> f64 {
        let upper = self.x.partition_point(|&x| x <= point);
        if upper == 0 {
            return self.tied_mean(0);
        }
        if upper == self.x.len() {
            return self.tied_mean(self.x.len() - 1);
        }
        let left = upper - 1;
        let fraction = (point - self.x[left]) / (self.x[upper] - self.x[left]);
        let lower = self.tied_mean(left);
        lower + fraction * (self.tied_mean(upper) - lower)
    }
}

pub fn curve(x: &[f64], y: &[f64], weights: Option<&[f64]>, span: f64) -> Result<Curve> {
    ensure!(
        x.len() >= 2,
        "voom trend needs at least two replicated genes"
    );
    let mut order: Vec<usize> = (0..x.len()).collect();
    order.sort_by(|&a, &b| x[a].total_cmp(&x[b]));
    let sorted_x: Vec<f64> = order.iter().map(|&i| x[i]).collect();
    if let Some(weights) = weights {
        let sorted_y: Vec<f64> = order.iter().map(|&i| y[i]).collect();
        let sorted_w: Vec<f64> = order.iter().map(|&i| weights[i]).collect();
        weighted_curve(sorted_x, sorted_y, sorted_w, span)
    } else {
        let fitted = lowess(x, y, span, 3);
        let delta = 0.01 * (sorted_x[sorted_x.len() - 1] - sorted_x[0]);
        Ok(Curve {
            x: sorted_x,
            y: order.iter().map(|&i| fitted[i]).collect(),
            delta,
        })
    }
}

fn weighted_curve(x: Vec<f64>, y: Vec<f64>, w: Vec<f64>, span: f64) -> Result<Curve> {
    let n = x.len();
    let delta = if n <= 200 {
        0.0
    } else {
        let mut gaps: Vec<f64> = x.windows(2).map(|p| p[1] - p[0]).collect();
        gaps.sort_by(f64::total_cmp);
        let mut total = 0.0;
        let cumulative: Vec<f64> = gaps
            .into_iter()
            .map(|gap| {
                total += gap;
                total
            })
            .collect();
        (0..200)
            .map(|clusters| cumulative[n - 2 - clusters] / (200 - clusters) as f64)
            .fold(f64::INFINITY, f64::min)
    };
    let mut seeds = vec![0];
    let mut previous = 0;
    for i in 1..n - 1 {
        if x[i] - x[previous] > delta {
            seeds.push(i);
            previous = i;
        }
    }
    seeds.push(n - 1);
    let total_weight: f64 = w.iter().sum();
    let span_weight = total_weight * span;
    let windows: Vec<(usize, usize, f64)> = seeds
        .iter()
        .map(|&seed| {
            let (mut left, mut right, mut total) = (seed, seed, w[seed]);
            while total < span_weight && (left > 0 || right + 1 < n) {
                if right + 1 == n || (left > 0 && x[seed] - x[left - 1] < x[right + 1] - x[seed]) {
                    left -= 1;
                    total += w[left];
                } else {
                    right += 1;
                    total += w[right];
                }
            }
            let radius = (x[seed] - x[left]).max(x[right] - x[seed]);
            while left > 0 && x[left - 1] == x[left] {
                left -= 1;
            }
            while right + 1 < n && x[right + 1] == x[right] {
                right += 1;
            }
            (left, right, radius)
        })
        .collect();
    let mut robustness = vec![1.0; n];
    let mut fitted = vec![0.0; n];
    let subrange = (x[n - 1] - x[0]) / n as f64;
    for _ in 0..4 {
        for (&seed, &(left, right, radius)) in seeds.iter().zip(&windows) {
            let kernel = |j: usize| {
                if radius < 1e-7 {
                    1.0
                } else {
                    (1.0 - ((x[seed] - x[j]).abs() / radius).powi(3)).powi(3)
                }
            };
            let mut total = 0.0;
            let mut sum_x = 0.0;
            let mut sum_y = 0.0;
            for j in left..=right {
                let weight = w[j] * robustness[j] * kernel(j);
                total += weight;
                sum_x += weight * x[j];
                sum_y += weight * y[j];
            }
            ensure!(
                total > 0.0,
                "weighted voom LOWESS has an empty neighborhood"
            );
            let center = sum_x / total;
            let mean = sum_y / total;
            let mut variance = 0.0;
            let mut covariance = 0.0;
            for j in left..=right {
                let weight = w[j] * robustness[j] * kernel(j);
                variance += weight * (x[j] - center).powi(2);
                covariance += weight * (x[j] - center) * (y[j] - mean);
            }
            fitted[seed] = if radius < 1e-7 || variance < 1e-7 {
                mean
            } else {
                mean + (x[seed] - center) * covariance / variance
            };
        }
        for pair in seeds.windows(2) {
            let (left, right) = (pair[0], pair[1]);
            let distance = x[right] - x[left];
            for j in left + 1..right {
                fitted[j] = if distance > 1e-7 * subrange {
                    let slope = (fitted[right] - fitted[left]) / distance;
                    slope * x[j] + (fitted[right] - slope * x[right])
                } else {
                    0.5 * (fitted[right] + fitted[left])
                };
            }
        }
        let residuals: Vec<f64> = y
            .iter()
            .zip(&fitted)
            .map(|(&a, &b)| (a - b).abs())
            .collect();
        let scale = residuals.iter().sum::<f64>() / n as f64;
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| residuals[a].total_cmp(&residuals[b]));
        let mut cumulative = 0.0;
        let mut cmad = 0.0;
        for (position, &i) in order.iter().enumerate() {
            cumulative += w[i];
            if cumulative == total_weight / 2.0 {
                cmad = 3.0 * (residuals[i] + residuals[order[position + 1]]);
                break;
            }
            if cumulative > total_weight / 2.0 {
                cmad = 6.0 * residuals[i];
                break;
            }
        }
        if cmad <= 1e-7 * scale {
            break;
        }
        for i in 0..n {
            robustness[i] = if residuals[i] < cmad {
                (1.0 - (residuals[i] / cmad).powi(2)).powi(2)
            } else {
                0.0
            };
        }
    }
    ensure!(
        fitted.iter().all(|x| x.is_finite()),
        "nonfinite voom LOWESS trend"
    );
    Ok(Curve {
        x,
        y: fitted,
        delta,
    })
}

pub fn lowess(x: &[f64], y: &[f64], span: f64, robust_iterations: usize) -> Vec<f64> {
    let mut order: Vec<usize> = (0..x.len()).collect();
    order.sort_by(|&a, &b| x[a].total_cmp(&x[b]));
    let sorted_x: Vec<f64> = order.iter().map(|&i| x[i]).collect();
    let sorted_y: Vec<f64> = order.iter().map(|&i| y[i]).collect();
    let fitted = sorted_lowess(&sorted_x, &sorted_y, span, robust_iterations);
    let mut result = vec![0.0; x.len()];
    for (&i, value) in order.iter().zip(fitted) {
        result[i] = value;
    }
    result
}

fn sorted_lowess(x: &[f64], y: &[f64], span: f64, iterations: usize) -> Vec<f64> {
    let n = x.len();
    if n < 2 {
        return y.to_vec();
    }
    let width = ((span * n as f64 + 1e-7) as usize).clamp(2, n);
    let delta = 0.01 * (x[n - 1] - x[0]);
    let mut fitted = vec![0.0; n];
    let mut robust = vec![1.0; n];
    let mut residuals = vec![0.0; n];
    let mut weights = vec![0.0; n];
    for iteration in 0..=iterations {
        let mut left = 0;
        let mut right = width - 1;
        let mut previous: Option<usize> = None;
        let mut point = 0;
        loop {
            while right + 1 < n && x[point] - x[left] > x[right + 1] - x[point] {
                left += 1;
                right += 1;
            }
            let radius = (x[point] - x[left]).max(x[right] - x[point]);
            let mut end = right;
            while end + 1 < n && x[end + 1] == x[right] {
                end += 1;
            }
            weights.fill(0.0);
            let mut total = 0.0;
            for j in left..=end {
                let distance = (x[j] - x[point]).abs();
                let weight = if distance <= 0.001 * radius {
                    1.0
                } else if distance <= 0.999 * radius {
                    (1.0 - (distance / radius).powi(3)).powi(3)
                } else {
                    0.0
                };
                weights[j] = weight * robust[j];
                total += weights[j];
            }
            if total > 0.0 {
                for value in &mut weights[left..=end] {
                    *value /= total;
                }
                let center: f64 = (left..=end).map(|j| weights[j] * x[j]).sum();
                let variance: f64 = (left..=end)
                    .map(|j| weights[j] * (x[j] - center).powi(2))
                    .sum();
                if variance.sqrt() > 0.001 * (x[n - 1] - x[0]) {
                    let slope = (x[point] - center) / variance;
                    for j in left..=end {
                        weights[j] *= 1.0 + slope * (x[j] - center);
                    }
                }
                fitted[point] = (left..=end).map(|j| weights[j] * y[j]).sum();
            } else {
                fitted[point] = y[point];
            }
            if let Some(last) = previous {
                for j in last + 1..point {
                    let fraction = (x[j] - x[last]) / (x[point] - x[last]);
                    fitted[j] = fraction * fitted[point] + (1.0 - fraction) * fitted[last];
                }
            }
            let mut last = point;
            let cutoff = x[point] + delta;
            let mut next = point + 1;
            while next < n && x[next] <= cutoff {
                if x[next] == x[last] {
                    fitted[next] = fitted[last];
                    last = next;
                }
                next += 1;
            }
            previous = Some(last);
            if last == n - 1 {
                break;
            }
            point = (last + 1).max(next.saturating_sub(1));
        }
        if iteration == iterations {
            break;
        }
        for i in 0..n {
            residuals[i] = (y[i] - fitted[i]).abs();
        }
        let scale = residuals.iter().sum::<f64>() / n as f64;
        let mut sorted = residuals.clone();
        sorted.sort_by(f64::total_cmp);
        let center = if n.is_multiple_of(2) {
            3.0 * (sorted[n / 2 - 1] + sorted[n / 2])
        } else {
            6.0 * sorted[n / 2]
        };
        if center < 1e-7 * scale || center == 0.0 {
            break;
        }
        for i in 0..n {
            robust[i] = if residuals[i] <= 0.001 * center {
                1.0
            } else if residuals[i] <= 0.999 * center {
                (1.0 - (residuals[i] / center).powi(2)).powi(2)
            } else {
                0.0
            };
        }
    }
    fitted
}
