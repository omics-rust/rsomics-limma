use crate::Result;
use crate::model::{Counts, Dataset, Design};
use anyhow::{Context, anyhow, ensure};
use csv::StringRecord;
use nalgebra::DMatrix;
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

struct Table {
    columns: Vec<String>,
    rows: Vec<StringRecord>,
    sha256: String,
}

fn read_table(path: &Path) -> Result<Table> {
    let bytes = std::fs::read(path).with_context(|| format!("read {}", path.display()))?;
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(bytes.as_slice());
    let header = reader
        .headers()
        .with_context(|| format!("header in {}", path.display()))?
        .clone();
    ensure!(!header.is_empty(), "{} has no ID column", path.display());
    let rows = reader
        .records()
        .collect::<std::result::Result<Vec<_>, _>>()
        .with_context(|| format!("TSV rows in {}", path.display()))?;
    Ok(Table {
        columns: header.iter().skip(1).map(str::to_owned).collect(),
        rows,
        sha256: format!("{:x}", Sha256::digest(&bytes)),
    })
}

fn align(table: Table, samples: &[String], label: &str) -> Result<Vec<Vec<String>>> {
    let mut rows = BTreeMap::new();
    for row in table.rows {
        ensure!(!row[0].is_empty(), "{label} has an empty sample ID");
        let previous = rows.insert(
            row[0].to_owned(),
            row.iter().skip(1).map(str::to_owned).collect(),
        );
        ensure!(previous.is_none(), "duplicate {label} sample '{}'", &row[0]);
    }
    ensure!(
        rows.len() == samples.len(),
        "{label} sample identities differ from counts"
    );
    let aligned = samples
        .iter()
        .map(|sample| {
            rows.remove(sample)
                .ok_or_else(|| anyhow!("{label} is missing count sample '{sample}'"))
        })
        .collect::<Result<Vec<_>>>()?;
    ensure!(rows.is_empty(), "{label} has extra sample identities");
    Ok(aligned)
}

pub fn read_inputs(counts: &Path, samples: &Path, design: &Path) -> Result<Dataset> {
    let counts = read_table(counts)?;
    let samples = read_table(samples)?;
    let design = read_table(design)?;
    let input_sha256 = BTreeMap::from([
        ("counts".into(), counts.sha256.clone()),
        ("samples".into(), samples.sha256.clone()),
        ("design".into(), design.sha256.clone()),
    ]);
    let mut gene_ids = Vec::with_capacity(counts.rows.len());
    let mut values = Vec::with_capacity(counts.rows.len() * counts.columns.len());
    for row in counts.rows {
        gene_ids.push(row[0].to_owned());
        for (sample, value) in counts.columns.iter().zip(row.iter().skip(1)) {
            values.push(
                value
                    .trim()
                    .parse::<f64>()
                    .with_context(|| format!("count for gene '{}' sample '{sample}'", &row[0]))?,
            );
        }
    }
    let metadata_columns = samples.columns.clone();
    let names = design.columns.clone();
    let metadata = align(samples, &counts.columns, "metadata")?;
    let design_rows = align(design, &counts.columns, "design")?;
    let mut design_values = Vec::with_capacity(design_rows.len() * names.len());
    for (sample, row) in counts.columns.iter().zip(design_rows) {
        for (name, value) in names.iter().zip(row) {
            design_values.push(
                value
                    .trim()
                    .parse::<f64>()
                    .with_context(|| format!("design column '{name}' sample '{sample}'"))?,
            );
        }
    }
    let dataset = Dataset {
        design: Design {
            names,
            n_samples: counts.columns.len(),
            values: design_values,
        },
        counts: Counts {
            gene_ids,
            sample_ids: counts.columns,
            values,
        },
        metadata_columns,
        metadata,
        input_sha256,
    };
    validate(&dataset)?;
    Ok(dataset)
}

fn unique(names: &[String], label: &str) -> Result<()> {
    let mut seen = BTreeSet::new();
    for name in names {
        ensure!(!name.is_empty(), "empty {label}");
        ensure!(seen.insert(name), "duplicate {label} '{name}'");
    }
    Ok(())
}

pub fn validate(dataset: &Dataset) -> Result<()> {
    let counts = &dataset.counts;
    let design = &dataset.design;
    let n = counts.n_samples();
    let p = design.n_coefs();
    ensure!(counts.n_genes() > 0, "count matrix has no genes");
    ensure!(
        p >= 1 && n > p,
        "design requires n > p and at least one coefficient"
    );
    ensure!(
        design.n_samples == n && n.checked_mul(p) == Some(design.values.len()),
        "design dimensions differ from count samples"
    );
    ensure!(
        counts.n_genes().checked_mul(n) == Some(counts.values.len()),
        "count dimensions differ from IDs"
    );
    ensure!(
        dataset.metadata.len() == n
            && dataset
                .metadata
                .iter()
                .all(|row| row.len() == dataset.metadata_columns.len()),
        "sample metadata dimensions differ"
    );
    unique(&counts.gene_ids, "gene ID")?;
    unique(&counts.sample_ids, "sample ID")?;
    unique(&design.names, "design column")?;
    unique(&dataset.metadata_columns, "metadata column")?;
    for (index, &value) in counts.values.iter().enumerate() {
        ensure!(
            value.is_finite() && value >= 0.0,
            "gene '{}' sample '{}' requires a finite nonnegative count",
            counts.gene_ids[index / n],
            counts.sample_ids[index % n]
        );
    }
    ensure!(
        design.values.iter().all(|value| value.is_finite()),
        "design values must be finite"
    );
    ensure!(
        orthonormal_columns(&design.matrix())?.len() == p,
        "design is rank deficient at QR relative tolerance 1e-7"
    );
    ensure!(
        dataset.input_sha256.len() == 3
            && ["counts", "samples", "design"].iter().all(|name| {
                dataset.input_sha256.get(*name).is_some_and(|digest| {
                    digest.len() == 64
                        && digest
                            .bytes()
                            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
                })
            }),
        "input SHA-256 identities are incomplete or invalid"
    );
    Ok(())
}

pub(crate) fn orthonormal_columns(matrix: &DMatrix<f64>) -> Result<Vec<Vec<f64>>> {
    let mut basis: Vec<Vec<f64>> = Vec::with_capacity(matrix.ncols());
    for column in matrix.column_iter() {
        let norm = column
            .iter()
            .fold(0.0_f64, |norm, &value| norm.hypot(value));
        ensure!(norm.is_finite(), "design QR norm overflowed");
        if norm == 0.0 {
            continue;
        }
        let mut residual: Vec<_> = column.iter().map(|value| value / norm).collect();
        for _ in 0..2 {
            for direction in &basis {
                let projection = residual
                    .iter()
                    .zip(direction)
                    .map(|(x, q)| x * q)
                    .sum::<f64>();
                for (value, q) in residual.iter_mut().zip(direction) {
                    *value -= projection * q;
                }
            }
        }
        let norm = residual
            .iter()
            .fold(0.0_f64, |norm, &value| norm.hypot(value));
        if norm >= 1e-7 {
            residual.iter_mut().for_each(|value| *value /= norm);
            basis.push(residual);
        }
    }
    Ok(basis)
}

pub fn contrast(names: &[String], terms: &[String]) -> Result<Vec<f64>> {
    unique(names, "contrast coefficient")?;
    ensure!(
        terms.len() == names.len(),
        "contrast must name every design coefficient exactly once"
    );
    let mut weights = BTreeMap::new();
    for term in terms {
        let (name, value) = term
            .rsplit_once('=')
            .ok_or_else(|| anyhow!("contrast requires name=weight: '{term}'"))?;
        let value = value
            .trim()
            .parse::<f64>()
            .with_context(|| format!("contrast weight '{term}'"))?;
        ensure!(
            value.is_finite(),
            "contrast weight must be finite: '{term}'"
        );
        ensure!(
            weights.insert(name, value).is_none(),
            "duplicate contrast coefficient '{name}'"
        );
    }
    let contrast = names
        .iter()
        .map(|name| {
            weights
                .remove(name.as_str())
                .ok_or_else(|| anyhow!("contrast is missing coefficient '{name}'"))
        })
        .collect::<Result<Vec<_>>>()?;
    ensure!(weights.is_empty(), "contrast contains unknown coefficients");
    ensure!(
        contrast.iter().any(|&weight| weight != 0.0),
        "contrast must be nonzero"
    );
    Ok(contrast)
}
