use crate::Result;
use anyhow::{bail, ensure};

pub fn quantile(values: &[f64], probability: f64) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let position = (sorted.len() - 1) as f64 * probability;
    let lower = position.floor() as usize;
    let upper = (lower + 1).min(sorted.len() - 1);
    let fraction = position - lower as f64;
    (1.0 - fraction) * sorted[lower] + fraction * sorted[upper]
}

pub fn brent(
    mut function: impl FnMut(f64) -> Result<f64>,
    lower: f64,
    upper: f64,
    tolerance: f64,
) -> Result<f64> {
    let golden = (3.0 - 5.0_f64.sqrt()) / 2.0;
    let mut a = lower;
    let mut b = upper;
    let mut x = a + golden * (b - a);
    let mut w = x;
    let mut v = x;
    let mut fx = function(x)?;
    ensure!(fx.is_finite(), "nonfinite scalar objective");
    let mut fw = fx;
    let mut fv = fx;
    let mut displacement: f64 = 0.0;
    let mut previous: f64 = 0.0;
    for _ in 0..1000 {
        let middle = (a + b) / 2.0;
        let small = f64::EPSILON.sqrt() * x.abs() + tolerance / 3.0;
        let doubled = 2.0 * small;
        if (x - middle).abs() <= doubled - (b - a) / 2.0 {
            return Ok(x);
        }
        let mut parabolic = false;
        if previous.abs() > small {
            let r = (x - w) * (fx - fv);
            let q = (x - v) * (fx - fw);
            let mut numerator = (x - v) * q - (x - w) * r;
            let mut denominator = 2.0 * (q - r);
            if denominator > 0.0 {
                numerator = -numerator;
            }
            denominator = denominator.abs();
            let old_previous = previous;
            previous = displacement;
            if numerator.abs() < (denominator * old_previous / 2.0).abs()
                && numerator > denominator * (a - x)
                && numerator < denominator * (b - x)
            {
                displacement = numerator / denominator;
                let candidate = x + displacement;
                if candidate - a < doubled || b - candidate < doubled {
                    displacement = small.copysign(middle - x);
                }
                parabolic = true;
            }
        }
        if !parabolic {
            previous = if x < middle { b - x } else { a - x };
            displacement = golden * previous;
        }
        let candidate = x + if displacement.abs() >= small {
            displacement
        } else {
            small.copysign(displacement)
        };
        let value = function(candidate)?;
        ensure!(value.is_finite(), "nonfinite scalar objective");
        if value <= fx {
            if candidate < x {
                b = x;
            } else {
                a = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = candidate;
            fx = value;
        } else {
            if candidate < x {
                a = candidate;
            } else {
                b = candidate;
            }
            if value <= fw || w == x {
                v = w;
                fv = fw;
                w = candidate;
                fw = value;
            } else if value <= fv || v == x || v == w {
                v = candidate;
                fv = value;
            }
        }
    }
    bail!("scalar optimizer exceeded 1000 evaluations")
}

pub fn log_minus_digamma(mut x: f64) -> f64 {
    let mut correction = 0.0;
    while x < 16.0 {
        correction += 1.0 / x - (1.0 / x).ln_1p();
        x += 1.0;
    }
    let square = 1.0 / (x * x);
    correction
        + 0.5 / x
        + square
            * (1.0 / 12.0
                + square
                    * (-1.0 / 120.0
                        + square
                            * (1.0 / 252.0
                                + square
                                    * (-1.0 / 240.0
                                        + square * (1.0 / 132.0 - square * 691.0 / 32760.0)))))
}

pub fn trigamma(mut x: f64) -> f64 {
    let mut correction = 0.0;
    while x < 16.0 {
        correction += 1.0 / (x * x);
        x += 1.0;
    }
    let inverse = 1.0 / x;
    let square = inverse * inverse;
    correction
        + inverse
        + square / 2.0
        + inverse
            * square
            * (1.0 / 6.0
                + square
                    * (-1.0 / 30.0
                        + square
                            * (1.0 / 42.0
                                + square
                                    * (-1.0 / 30.0
                                        + square * (5.0 / 66.0 - square * 691.0 / 2730.0)))))
}
