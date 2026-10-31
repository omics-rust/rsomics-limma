mod artifact;
mod cli;
mod ebayes;
mod expression;
mod input;
mod linear;
mod model;
mod normalization;
mod number;
mod numerics;
mod poisson;
mod results;
mod smoothing;
mod voom;

type Result<T> = anyhow::Result<T>;

fn analyze(
    dataset: model::Dataset,
    workers: usize,
    contrast: Vec<f64>,
) -> Result<artifact::Analysis> {
    input::validate(&dataset)?;
    anyhow::ensure!(workers > 0, "worker count must be positive");
    anyhow::ensure!(
        contrast.len() == dataset.design.n_coefs()
            && contrast.iter().all(|x| x.is_finite())
            && contrast.iter().any(|&x| x != 0.0),
        "invalid complete named contrast"
    );
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(workers)
        .build()?;
    pool.install(|| {
        let normalization = normalization::prepare(&dataset)?;
        let base = voom::fit(&dataset, &normalization)?;
        let initial_results = results::derive(&base.linear, &contrast)?;
        Ok(artifact::Analysis::Counts(Box::new(artifact::CountAnalysis {
            schema: 1,
            product: "rsomics-limma".into(),
            product_version: env!("CARGO_PKG_VERSION").into(),
            upstream: "edgeR4.10.7/27ba7b39d819e68df2ecdea9f32a4b6039b2a0b5;limma3.68.5/825d1c8330d024adb6c54af66d79b0f4ec258530".into(),
            workers,
            parameters: artifact::parameters(),
            dataset,
            normalization,
            base,
            initial_results,
        })))
    })
}

fn analyze_expression(
    dataset: model::ExpressionDataset,
    workers: usize,
    contrast: Vec<f64>,
) -> Result<artifact::Analysis> {
    input::validate_expression(&dataset)?;
    anyhow::ensure!(workers > 0, "worker count must be positive");
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(workers)
        .build()?;
    pool.install(|| {
        let base = expression::fit(&dataset)?;
        let initial_results = results::derive(&base, &contrast)?;
        Ok(artifact::Analysis::Expression(Box::new(
            artifact::ExpressionAnalysis {
                schema: 2,
                product: "rsomics-limma".into(),
                product_version: env!("CARGO_PKG_VERSION").into(),
                upstream: "limma3.68.5/825d1c8330d024adb6c54af66d79b0f4ec258530".into(),
                workers,
                input_kind: "log-expression".into(),
                parameters: artifact::expression_parameters(),
                source_provenance: artifact::expression_provenance(),
                dataset,
                base,
                initial_results,
            },
        )))
    })
}

pub fn run_cli() -> std::process::ExitCode {
    cli::run()
}
