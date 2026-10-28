use crate::{
    Result,
    model::{Dataset, LibraryState},
    number::Number,
    results::ContrastResult,
    voom::BaseFit,
};
use anyhow::{Context, ensure};
use serde::{Deserialize, Serialize};
use std::{fs::File, io::BufReader, path::Path};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Analysis {
    pub schema: u32,
    pub product: String,
    pub product_version: String,
    pub upstream: String,
    pub workers: usize,
    pub parameters: serde_json::Value,
    pub dataset: Dataset,
    pub normalization: LibraryState,
    pub base: BaseFit,
    pub initial_results: ContrastResult,
}

pub fn named(path: &Path) -> Result<()> {
    ensure!(
        path != Path::new("-"),
        "this operation requires a named file, not '-'"
    );
    Ok(())
}

pub fn parameters() -> serde_json::Value {
    serde_json::json!({
        "input_design":"complete named numeric matrix; no formula inference",
        "construction":{"lib.size":null,"norm.factors":null,"group":null,"genes":null,"remove.zeros":false},
        "filtering":{"explicit_design":true,"group":null,"lib.size":null,"min.count":10,"min.total.count":15,"large.n":10,"min.prop":0.7},
        "subset":{"keep.lib.sizes":false},
        "normalization":{"method":"TMM","refColumn":null,"logratioTrim":0.3,"sumTrim":0.05,"doWeighting":true,"Acutoff":-1e10,"p":0.75},
        "voom":{"block":null,"prior.weights":null,"sample.weights":false,"var.design":null,"var.group":null,"prior.n":10,"lib.size":null,"offset":null,"offset.prior":null,"normalize.method":"none","span":0.5,"adaptive.span":true,"plot":false,"save.plot":true,"keep.EList":true,"derived_span_rule":"min(.3+.7*(50/kept_genes)^(1/3),1)"},
        "moderation":{"proportion":0.01,"stdev.coef.lim":[0.1,4.0],"trend":false,"span":null,"robust":false,"winsor.tail.p":[0.05,0.1],"legacy":null},
        "topTable":{"coef":1,"number":"Inf","genelist":"fit$genes","adjust.method":"BH","sort.by":"none","resort.by":null,"p.value":1,"fc":null,"lfc":0,"confint":false},
        "reload":"saved full base reused; contrast/moderation only; no linear-model refit"
    })
}

pub fn validate(analysis: &Analysis) -> Result<()> {
    ensure!(
        analysis.schema == 1 && analysis.product == "rsomics-limma",
        "unsupported analysis product/schema"
    );
    crate::input::validate(&analysis.dataset)?;
    let data = &analysis.dataset;
    let lib = &analysis.normalization;
    let base = &analysis.base;
    let n = data.counts.n_samples();
    let p = data.design.n_coefs();
    let m = lib.kept_indices.len();
    ensure!(analysis.workers > 0, "invalid saved worker count");
    ensure!(
        lib.keep.len() == data.counts.n_genes()
            && m >= 2
            && lib.kept_indices
                == lib
                    .keep
                    .iter()
                    .enumerate()
                    .filter_map(|(i, &keep)| keep.then_some(i))
                    .collect::<Vec<_>>(),
        "invalid retained-feature mapping"
    );
    let mut raw = vec![0.0; n];
    let mut kept = vec![0.0; n];
    for i in 0..data.counts.n_genes() {
        for (j, &value) in data.counts.row(i).iter().enumerate() {
            raw[j] += value;
            if lib.keep[i] {
                kept[j] += value;
            }
        }
    }
    ensure!(
        raw == lib.raw_library_sizes && kept == lib.postfilter_library_sizes,
        "saved counts disagree with library totals"
    );
    for values in [
        &lib.raw_library_sizes,
        &lib.postfilter_library_sizes,
        &lib.tmm_factors,
        &lib.effective_library_sizes,
    ] {
        ensure!(
            values.len() == n && values.iter().all(|x| x.is_finite() && *x > 0.0),
            "invalid saved library state"
        );
    }
    ensure!(
        lib.offsets.len() == n && lib.tmm_reference_sample < n,
        "invalid saved normalization dimensions"
    );
    for i in 0..n {
        let expected = lib.postfilter_library_sizes[i] * lib.tmm_factors[i];
        ensure!(
            (expected - lib.effective_library_sizes[i]).abs() <= 1e-12 * expected
                && (expected.ln() - lib.offsets[i]).abs() <= 1e-12,
            "inconsistent effective library/offset"
        );
    }
    ensure!(
        base.design == data.design
            && base.rank == p
            && base.method == "ls"
            && base.pivot.len() == p,
        "invalid saved design/rank/method"
    );
    let mut pivot = base.pivot.clone();
    ensure!(
        base.targets.len() == n
            && base
                .targets
                .iter()
                .enumerate()
                .all(|(i, t)| t.sample == data.counts.sample_ids[i]
                    && t.lib_size == lib.postfilter_library_sizes[i]
                    && t.norm_factor == lib.tmm_factors[i]
                    && t.metadata == data.metadata[i]),
        "invalid saved sample targets"
    );
    pivot.sort_unstable();
    ensure!(pivot == (1..=p).collect::<Vec<_>>(), "invalid saved pivot");
    for (matrix, columns) in [
        (&base.coefficients, p),
        (&base.unscaled_se, p),
        (&base.expression, n),
        (&base.weights, n),
    ] {
        ensure!(
            matrix.len() == m
                && matrix
                    .iter()
                    .all(|row| row.len() == columns && row.iter().all(|x| x.is_finite())),
            "invalid saved fitted matrix"
        );
    }
    ensure!(
        base.weights.iter().flatten().all(|&x| x > 0.0)
            && base.unscaled_se.iter().flatten().all(|&x| x > 0.0),
        "nonpositive saved weights/SE"
    );
    ensure!(
        base.sigma.len() == m
            && base.residual_df.len() == m
            && base.amean.len() == m
            && base.zero_mask.len() == m
            && base.zero_probes.len() == m,
        "invalid saved gene dimensions"
    );
    for i in 0..m {
        ensure!(
            base.zero_mask[i].len() == n
                && base.residual_df[i] <= n - p
                && base.amean[i].is_finite(),
            "invalid saved zero-loss/abundance state"
        );
        ensure!(
            matches!(base.sigma[i],Number::Finite(x) if x>=0.0 && x.is_finite())
                || (base.sigma[i] == Number::Na && base.residual_df[i] == 0),
            "invalid residual sigma state"
        );
    }
    ensure!(
        base.covariance.len() == p
            && base
                .covariance
                .iter()
                .all(|row| row.len() == p && row.iter().all(|x| x.is_finite())),
        "invalid saved global coefficient covariance"
    );
    ensure!(
        base.voom_line.x.len() >= 2
            && base.voom_line.x.len() == base.voom_line.y.len()
            && base.voom_line.x.windows(2).all(|pair| pair[0] <= pair[1])
            && base
                .voom_line
                .x
                .iter()
                .chain(&base.voom_line.y)
                .all(|x| x.is_finite()),
        "invalid saved voom curve"
    );
    let result = &analysis.initial_results;
    ensure!(
        result.contrast.len() == p
            && result.contrast.iter().all(|x| x.is_finite())
            && result.contrast.iter().any(|&x| x != 0.0)
            && result.rows.len() == m,
        "invalid saved contrast/results"
    );
    ensure!(
        result.full_base_reused && !result.linear_model_refit,
        "invalid saved derivation semantics"
    );
    Ok(())
}

pub fn save(path: &Path, analysis: &Analysis) -> Result<()> {
    named(path)?;
    validate(analysis)?;
    rsomics_common::write_output(Some(path), |writer| {
        serde_json::to_writer(&mut *writer, analysis).map_err(crate::cli::json_error)?;
        writer
            .write_all(b"\n")
            .map_err(rsomics_common::RsomicsError::Io)?;
        Ok(())
    })?;
    Ok(())
}

pub fn load(path: &Path) -> Result<Analysis> {
    named(path)?;
    let file = File::open(path).with_context(|| format!("opening bundle {}", path.display()))?;
    let analysis = serde_json::from_reader(BufReader::new(file))
        .with_context(|| format!("reading bundle {}", path.display()))?;
    validate(&analysis)?;
    Ok(analysis)
}

pub fn write_results(path: &Path, analysis: &Analysis, result: &ContrastResult) -> Result<()> {
    named(path)?;
    rsomics_common::write_output(Some(path), |output| {
        write_results_to(output, analysis, result)
    })?;
    Ok(())
}

pub fn write_results_pair(
    path: &Path,
    state: &Path,
    analysis: &Analysis,
    result: &ContrastResult,
) -> Result<()> {
    use std::io::{BufWriter, Write};
    named(path)?;
    named(state)?;
    let mut encoded = serde_json::to_vec(result)?;
    encoded.push(b'\n');
    rsomics_common::write_atomic_pair(path, state, |table, json| {
        let mut table = BufWriter::new(table);
        write_results_to(&mut table, analysis, result)?;
        table.flush().map_err(rsomics_common::RsomicsError::Io)?;
        json.write_all(&encoded)
            .map_err(rsomics_common::RsomicsError::Io)?;
        Ok(())
    })?;
    Ok(())
}

fn write_results_to(
    output: &mut dyn std::io::Write,
    analysis: &Analysis,
    result: &ContrastResult,
) -> rsomics_common::Result<()> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(output);
    writer
        .write_record(["gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"])
        .map_err(csv_error)?;
    for (&gene, row) in analysis.normalization.kept_indices.iter().zip(&result.rows) {
        writer
            .write_record([
                analysis.dataset.counts.gene_ids[gene].clone(),
                row.log_fc.to_string(),
                row.ave_expr.to_string(),
                row.t.to_string(),
                row.p_value.to_string(),
                row.fdr.to_string(),
                row.lods.to_string(),
            ])
            .map_err(csv_error)?;
    }
    writer.flush().map_err(rsomics_common::RsomicsError::Io)?;
    Ok(())
}

pub fn write_export(path: &Path, analysis: &Analysis, kind: crate::cli::ExportKind) -> Result<()> {
    named(path)?;
    let base = &analysis.base;
    rsomics_common::write_output(Some(path), |output| {
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(output);
        let mut header = vec!["gene".to_string()];
        header.extend(analysis.dataset.counts.sample_ids.iter().cloned());
        writer.write_record(header).map_err(csv_error)?;
        for (i, &gene) in analysis.normalization.kept_indices.iter().enumerate() {
            let mut row = vec![analysis.dataset.counts.gene_ids[gene].clone()];
            for sample in 0..base.design.n_samples {
                let value = match kind {
                    crate::cli::ExportKind::LogCpm => base.expression[i][sample],
                    crate::cli::ExportKind::Weights => base.weights[i][sample],
                    crate::cli::ExportKind::Fitted | crate::cli::ExportKind::Residuals => {
                        let fitted = base
                            .design
                            .row(sample)
                            .iter()
                            .zip(&base.coefficients[i])
                            .map(|(&x, &beta)| x * beta)
                            .sum::<f64>();
                        if matches!(kind, crate::cli::ExportKind::Residuals) {
                            base.expression[i][sample] - fitted
                        } else {
                            fitted
                        }
                    }
                };
                row.push(value.to_string());
            }
            writer.write_record(row).map_err(csv_error)?;
        }
        writer.flush().map_err(rsomics_common::RsomicsError::Io)?;
        Ok(())
    })?;
    Ok(())
}

fn csv_error(error: csv::Error) -> rsomics_common::RsomicsError {
    match error.into_kind() {
        csv::ErrorKind::Io(error) => rsomics_common::RsomicsError::Io(error),
        kind => rsomics_common::RsomicsError::InvalidInput(format!("{kind:?}")),
    }
}
