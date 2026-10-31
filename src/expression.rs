use crate::{
    Result,
    linear::System,
    model::{ExpressionDataset, LinearFit},
};
use anyhow::ensure;
use rayon::prelude::*;

pub fn fit(data: &ExpressionDataset) -> Result<LinearFit> {
    let system = System::new(&data.design, None, None)?;
    let p = data.design.n_coefs();
    ensure!(system.rank == p, "ordinary design lost rank");
    let expression: Vec<Vec<f64>> = (0..data.expression.gene_ids.len())
        .map(|gene| data.expression.row(gene).to_vec())
        .collect();
    let rows = expression
        .par_iter()
        .map(|row| system.fit(row))
        .collect::<Result<Vec<_>>>()?;
    Ok(LinearFit {
        design: data.design.clone(),
        coefficients: rows.iter().map(|row| row.coefficients.clone()).collect(),
        unscaled_se: vec![system.unscaled_se.clone(); rows.len()],
        sigma: rows.iter().map(|row| row.sigma).collect(),
        residual_df: rows.iter().map(|row| row.residual_df).collect(),
        covariance: system
            .covariance
            .chunks_exact(p)
            .map(|row| row.to_vec())
            .collect(),
        rank: system.rank,
        pivot: system.pivot.iter().map(|index| index + 1).collect(),
        method: "ls".into(),
        amean: expression
            .iter()
            .map(|row| row.iter().sum::<f64>() / row.len() as f64)
            .collect(),
        expression,
    })
}
