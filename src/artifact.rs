use crate::{
    Result,
    model::{Dataset, Design, ExpressionDataset, LibraryState, LinearFit},
    number::Number,
    results::ContrastResult,
    voom::BaseFit,
};
use anyhow::{Context, ensure};
use serde::{Deserialize, Serialize};
use std::{fs::File, io::BufReader, path::Path};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CountAnalysis {
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

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ExpressionAnalysis {
    pub schema: u32,
    pub product: String,
    pub product_version: String,
    pub upstream: String,
    pub workers: usize,
    pub input_kind: String,
    pub parameters: serde_json::Value,
    pub source_provenance: serde_json::Value,
    pub dataset: ExpressionDataset,
    pub base: LinearFit,
    pub initial_results: ContrastResult,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
#[serde(untagged)]
pub enum Analysis {
    Counts(Box<CountAnalysis>),
    Expression(Box<ExpressionAnalysis>),
}

impl Analysis {
    pub fn base(&self) -> &LinearFit {
        match self {
            Self::Counts(value) => &value.base.linear,
            Self::Expression(value) => &value.base,
        }
    }

    pub fn design(&self) -> &Design {
        &self.base().design
    }

    pub fn sample_ids(&self) -> &[String] {
        match self {
            Self::Counts(value) => &value.dataset.counts.sample_ids,
            Self::Expression(value) => &value.dataset.expression.sample_ids,
        }
    }

    fn gene_id(&self, row: usize) -> &str {
        match self {
            Self::Counts(value) => {
                &value.dataset.counts.gene_ids[value.normalization.kept_indices[row]]
            }
            Self::Expression(value) => &value.dataset.expression.gene_ids[row],
        }
    }

    pub fn raw_genes(&self) -> usize {
        match self {
            Self::Counts(value) => value.dataset.counts.n_genes(),
            Self::Expression(value) => value.dataset.expression.gene_ids.len(),
        }
    }
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

pub fn expression_parameters() -> serde_json::Value {
    let count = parameters();
    serde_json::json!({
        "input_kind":"log-expression; supplied values already logged and unchanged",
        "input_design":"complete named numeric matrix; no formula inference",
        "lmFit":{"method":"ls","weights":null,"block":null,"ndups":1,"spacing":1,"correlation":"unused"},
        "moderation":count["moderation"], "topTable":count["topTable"],
        "reload":"saved full base reused; contrast/moderation only; no linear-model refit"
    })
}

pub fn expression_provenance() -> serde_json::Value {
    serde_json::json!({
        "reference":"limma3.68.5/825d1c8330d024adb6c54af66d79b0f4ec258530",
        "source_run":41184814043_u64,
        "source_manifest_sha256":"c374d2be76ce0d0ee25aae2eeb14548ddb02d9866a0266e814bd39c5dec5c19f",
        "source_witness_sha256":"ee4a903ce605ac51d059469d738020fb31e02ac45332156ce30ead5abaece4ab",
        "source_reload_witness_sha256":"87ca32ca6d0a199c7b6c729ae99096242fcc223e5ef8eee94bde3e5ae922af5f",
        "feature_identity":"native retains supplied IDs; source one-feature coefficient name loss and table ID 1 are compared through the retained positional mapping",
        "scope":"method source reference, not provenance of user-supplied expression or biological validation"
    })
}

pub fn validate(analysis: &Analysis) -> Result<()> {
    match analysis {
        Analysis::Counts(value) => validate_count(value),
        Analysis::Expression(value) => {
            ensure!(
                value.schema == 2
                    && value.product == "rsomics-limma"
                    && value.input_kind == "log-expression",
                "unsupported ordinary analysis product/schema/input kind"
            );
            crate::input::validate_expression(&value.dataset)?;
            let n = value.dataset.expression.sample_ids.len();
            let m = value.dataset.expression.gene_ids.len();
            ensure!(value.workers > 0, "invalid saved worker count");
            validate_linear(&value.base, &value.dataset.design, m)?;
            ensure!(
                value
                    .base
                    .residual_df
                    .iter()
                    .all(|&df| df == n - value.dataset.design.n_coefs()),
                "ordinary residual df differs from full design"
            );
            ensure!(
                value.base.expression.iter().flatten().copied().eq(value
                    .dataset
                    .expression
                    .values
                    .iter()
                    .copied()),
                "saved supplied expression differs from base observations"
            );
            validate_result(&value.initial_results, value.dataset.design.n_coefs(), m)
        }
    }
}

fn validate_linear(base: &LinearFit, design: &Design, m: usize) -> Result<()> {
    let n = design.n_samples;
    let p = design.n_coefs();
    ensure!(
        base.design == *design && base.rank == p && base.method == "ls" && base.pivot.len() == p,
        "invalid saved design/rank/method"
    );
    let mut pivot = base.pivot.clone();
    pivot.sort_unstable();
    ensure!(pivot == (1..=p).collect::<Vec<_>>(), "invalid saved pivot");
    for (matrix, rows, columns) in [
        (&base.coefficients, m, p),
        (&base.unscaled_se, m, p),
        (&base.expression, m, n),
        (&base.covariance, p, p),
    ] {
        ensure!(
            matrix.len() == rows
                && matrix
                    .iter()
                    .all(|row| row.len() == columns && row.iter().all(|x| x.is_finite())),
            "invalid saved linear matrix"
        );
    }
    ensure!(
        base.unscaled_se.iter().flatten().all(|&value| value > 0.0),
        "nonpositive saved unscaled SE"
    );
    ensure!(
        base.sigma.len() == m && base.residual_df.len() == m && base.amean.len() == m,
        "invalid saved linear gene dimensions"
    );
    for i in 0..m {
        ensure!(
            base.residual_df[i] <= n - p && base.amean[i].is_finite(),
            "invalid saved residual df/abundance"
        );
        ensure!(
            matches!(base.sigma[i], Number::Finite(value) if value.is_finite() && value >= 0.0)
                || (base.sigma[i] == Number::Na && base.residual_df[i] == 0),
            "invalid residual sigma state"
        );
    }
    Ok(())
}

fn validate_result(result: &ContrastResult, p: usize, m: usize) -> Result<()> {
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
    ensure!(
        result.rows.iter().all(|row| row.log_fc.is_finite()
            && row.unscaled_se.is_finite()
            && row.unscaled_se > 0.0),
        "invalid saved finite effect/unscaled SE"
    );
    Ok(())
}

fn validate_count(analysis: &CountAnalysis) -> Result<()> {
    ensure!(
        analysis.schema == 1 && analysis.product == "rsomics-limma",
        "unsupported analysis product/schema"
    );
    crate::input::validate(&analysis.dataset)?;
    let data = &analysis.dataset;
    let lib = &analysis.normalization;
    let base = &analysis.base;
    let fit = &base.linear;
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
        fit.design == data.design && fit.rank == p && fit.method == "ls" && fit.pivot.len() == p,
        "invalid saved design/rank/method"
    );
    let mut pivot = fit.pivot.clone();
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
        (&fit.coefficients, p),
        (&fit.unscaled_se, p),
        (&fit.expression, n),
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
            && fit.unscaled_se.iter().flatten().all(|&x| x > 0.0),
        "nonpositive saved weights/SE"
    );
    ensure!(
        fit.sigma.len() == m
            && fit.residual_df.len() == m
            && fit.amean.len() == m
            && base.zero_mask.len() == m
            && base.zero_probes.len() == m,
        "invalid saved gene dimensions"
    );
    for i in 0..m {
        ensure!(
            base.zero_mask[i].len() == n && fit.residual_df[i] <= n - p && fit.amean[i].is_finite(),
            "invalid saved zero-loss/abundance state"
        );
        ensure!(
            matches!(fit.sigma[i],Number::Finite(x) if x>=0.0 && x.is_finite())
                || (fit.sigma[i] == Number::Na && fit.residual_df[i] == 0),
            "invalid residual sigma state"
        );
    }
    ensure!(
        fit.covariance.len() == p
            && fit
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
    validate_result(&analysis.initial_results, p, m)
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
    let value: serde_json::Value = serde_json::from_reader(BufReader::new(file))
        .with_context(|| format!("reading bundle {}", path.display()))?;
    let analysis = match value.get("schema").and_then(serde_json::Value::as_u64) {
        Some(1) => {
            ensure!(
                value.get("input_kind").is_none(),
                "schema1 count bundle cannot declare a different input kind"
            );
            Analysis::Counts(Box::new(serde_json::from_value(value)?))
        }
        Some(2) => Analysis::Expression(Box::new(serde_json::from_value(value)?)),
        _ => anyhow::bail!("unsupported analysis schema"),
    };
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
    for (gene, row) in result.rows.iter().enumerate() {
        writer
            .write_record([
                analysis.gene_id(gene).to_string(),
                row.log_fc.to_string(),
                row.ave_expr.to_string(),
                row.t.text(),
                row.p_value.text(),
                row.fdr.text(),
                row.lods.text(),
            ])
            .map_err(csv_error)?;
    }
    writer.flush().map_err(rsomics_common::RsomicsError::Io)?;
    Ok(())
}

pub fn write_export(path: &Path, analysis: &Analysis, kind: crate::cli::ExportKind) -> Result<()> {
    named(path)?;
    let base = analysis.base();
    let weights = match (analysis, kind) {
        (Analysis::Counts(value), crate::cli::ExportKind::Weights) => Some(&value.base.weights),
        (
            Analysis::Expression(_),
            crate::cli::ExportKind::Weights | crate::cli::ExportKind::LogCpm,
        ) => anyhow::bail!("log-cpm and weights exports require a count-input bundle"),
        _ => None,
    };
    rsomics_common::write_output(Some(path), |output| {
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(output);
        let mut header = vec!["gene".to_string()];
        header.extend(analysis.sample_ids().iter().cloned());
        writer.write_record(header).map_err(csv_error)?;
        for i in 0..base.coefficients.len() {
            let mut row = vec![analysis.gene_id(i).to_string()];
            for sample in 0..base.design.n_samples {
                let value = match kind {
                    crate::cli::ExportKind::LogCpm | crate::cli::ExportKind::Expression => {
                        base.expression[i][sample]
                    }
                    crate::cli::ExportKind::Weights => {
                        weights.expect("weights variant was validated")[i][sample]
                    }
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
