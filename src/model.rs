use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Counts {
    pub gene_ids: Vec<String>,
    pub sample_ids: Vec<String>,
    pub values: Vec<f64>,
}

impl Counts {
    pub fn n_samples(&self) -> usize {
        self.sample_ids.len()
    }

    pub fn n_genes(&self) -> usize {
        self.gene_ids.len()
    }

    pub fn row(&self, gene: usize) -> &[f64] {
        let n = self.n_samples();
        &self.values[gene * n..(gene + 1) * n]
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Design {
    pub names: Vec<String>,
    pub n_samples: usize,
    pub values: Vec<f64>,
}

impl Design {
    pub fn n_coefs(&self) -> usize {
        self.names.len()
    }

    pub fn row(&self, sample: usize) -> &[f64] {
        let p = self.n_coefs();
        &self.values[sample * p..(sample + 1) * p]
    }

    pub fn matrix(&self) -> nalgebra::DMatrix<f64> {
        nalgebra::DMatrix::from_row_slice(self.n_samples, self.n_coefs(), &self.values)
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Dataset {
    pub counts: Counts,
    pub design: Design,
    pub metadata_columns: Vec<String>,
    pub metadata: Vec<Vec<String>>,
    pub input_sha256: BTreeMap<String, String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FilterThresholds {
    pub cpm_cutoff: f64,
    pub min_sample_size: f64,
    pub min_total_count: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LibraryState {
    pub keep: Vec<bool>,
    pub kept_indices: Vec<usize>,
    pub raw_library_sizes: Vec<f64>,
    pub postfilter_library_sizes: Vec<f64>,
    pub tmm_factors: Vec<f64>,
    pub effective_library_sizes: Vec<f64>,
    pub offsets: Vec<f64>,
    pub filter: FilterThresholds,
    pub tmm_reference_sample: usize,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "kebab-case")]
pub enum SolverStatus {
    Converged { iterations: usize, reason: String },
    AllZero,
}

#[derive(Clone, Debug)]
pub struct ScalarFit {
    pub coefficient: f64,
    pub status: SolverStatus,
}
