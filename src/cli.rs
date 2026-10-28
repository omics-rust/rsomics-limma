use crate::{artifact, input, results};
use clap::{Args, Parser, Subcommand, ValueEnum};
use rsomics_common::{OutputArgs, RsomicsError, ThreadArgs, ToolMeta};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(
    name = "rsomics-limma",
    version,
    about = "Current voom weighted linear-model count analysis",
    arg_required_else_help = true
)]
struct Cli {
    #[command(flatten)]
    output: OutputArgs,
    #[command(flatten)]
    threads: ThreadArgs,
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Fit counts, filter/TMM and the current voom linear model
    Run(Run),
    /// Test a named contrast by moderating saved full base-fit state
    Results(Derive),
    /// Export voom logCPM, precision weights, fitted values or residuals
    Export(Export),
}

#[derive(Debug, Args)]
struct Contrast {
    /// NAME=WEIGHT for every design column, including zero weights
    #[arg(long = "contrast", required = true)]
    terms: Vec<String>,
}

#[derive(Debug, Args)]
struct Run {
    #[arg(long)]
    counts: PathBuf,
    #[arg(long)]
    samples: PathBuf,
    #[arg(long)]
    design: PathBuf,
    #[arg(long)]
    bundle: PathBuf,
    #[command(flatten)]
    contrast: Contrast,
}

#[derive(Debug, Args)]
struct Derive {
    #[arg(long)]
    bundle: PathBuf,
    #[arg(long)]
    output: PathBuf,
    /// Save complete derived contrast and moderation state as a paired JSON output
    #[arg(long)]
    state: Option<PathBuf>,
    #[command(flatten)]
    contrast: Contrast,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum ExportKind {
    LogCpm,
    Weights,
    Fitted,
    Residuals,
}

#[derive(Debug, Args)]
struct Export {
    #[arg(long)]
    bundle: PathBuf,
    #[arg(long)]
    output: PathBuf,
    #[arg(long, value_enum)]
    kind: ExportKind,
}

pub fn run() -> std::process::ExitCode {
    let cli = rsomics_help::parse::<Cli>();
    rsomics_common::run(
        &cli.output,
        ToolMeta {
            name: "rsomics-limma",
            version: env!("CARGO_PKG_VERSION"),
        },
        || {
            let workers = cli.threads.requested().map(|n| n.get()).unwrap_or_else(|| {
                std::thread::available_parallelism()
                    .map(|n| n.get())
                    .unwrap_or(1)
            });
            execute(cli.command, workers).map_err(process_error)
        },
    )
}

fn execute(command: Command, workers: usize) -> crate::Result<serde_json::Value> {
    match command {
        Command::Run(args) => {
            for path in [&args.counts, &args.samples, &args.design, &args.bundle] {
                artifact::named(path)?;
            }
            rsomics_common::reject_output_alias(
                &args.bundle,
                [
                    args.counts.as_path(),
                    args.samples.as_path(),
                    args.design.as_path(),
                ],
            )?;
            let dataset = input::read_inputs(&args.counts, &args.samples, &args.design)?;
            let contrast = input::contrast(&dataset.design.names, &args.contrast.terms)?;
            let analysis = crate::analyze(dataset, workers, contrast)?;
            artifact::save(&args.bundle, &analysis)?;
            Ok(
                serde_json::json!({"operation":"run","raw_genes":analysis.dataset.counts.n_genes(),"retained_genes":analysis.base.coefficients.len(),"workers":workers}),
            )
        }
        Command::Results(args) => {
            artifact::named(&args.output)?;
            rsomics_common::reject_output_alias(&args.output, [args.bundle.as_path()])?;
            if let Some(state) = &args.state {
                artifact::named(state)?;
                rsomics_common::reject_output_alias(
                    state,
                    [args.bundle.as_path(), args.output.as_path()],
                )?;
            }
            let analysis = artifact::load(&args.bundle)?;
            let contrast = input::contrast(&analysis.dataset.design.names, &args.contrast.terms)?;
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()?;
            let table = pool.install(|| results::derive(&analysis.base, &contrast))?;
            if let Some(state) = &args.state {
                artifact::write_results_pair(&args.output, state, &analysis, &table)?;
            } else {
                artifact::write_results(&args.output, &analysis, &table)?;
            }
            Ok(
                serde_json::json!({"operation":"results","genes":table.rows.len(),"full_base_reused":true,"linear_model_refit":false,"workers":workers}),
            )
        }
        Command::Export(args) => {
            artifact::named(&args.output)?;
            rsomics_common::reject_output_alias(&args.output, [args.bundle.as_path()])?;
            let analysis = artifact::load(&args.bundle)?;
            artifact::write_export(&args.output, &analysis, args.kind)?;
            Ok(
                serde_json::json!({"operation":"export","genes":analysis.base.coefficients.len(),"full_base_reused":true}),
            )
        }
    }
}

pub fn json_error(error: serde_json::Error) -> RsomicsError {
    if error.is_io() {
        RsomicsError::Io(error.into())
    } else {
        RsomicsError::InvalidInput(error.to_string())
    }
}

fn process_error(error: anyhow::Error) -> RsomicsError {
    let error = match error.downcast::<RsomicsError>() {
        Ok(error) => return error,
        Err(error) => error,
    };
    let error = match error.downcast::<serde_json::Error>() {
        Ok(error) => return json_error(error),
        Err(error) => error,
    };
    if let Some(io) = error.downcast_ref::<std::io::Error>() {
        return RsomicsError::Io(std::io::Error::new(io.kind(), format!("{error:#}")));
    }
    if let Some(csv) = error.downcast_ref::<csv::Error>()
        && let csv::ErrorKind::Io(io) = csv.kind()
    {
        return RsomicsError::Io(std::io::Error::new(io.kind(), format!("{error:#}")));
    }
    RsomicsError::InvalidInput(format!("{error:#}"))
}
