use serde_json::Value;
use std::path::Path;
use std::process::{Command, Output};

fn invoke(arguments: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_rsomics-limma"))
        .args(arguments)
        .output()
        .unwrap()
}

fn success(arguments: &[&str]) {
    let output = invoke(arguments);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn read(path: &Path) -> Value {
    serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap()
}

#[test]
fn current_voom_saved_lifecycle() {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let directory = tempfile::tempdir().unwrap();
    let bundle = directory.path().join("analysis.json");
    let parallel = directory.path().join("parallel.json");
    let counts = fixture.join("counts.tsv");
    let samples = fixture.join("samples.tsv");
    let design = fixture.join("design.tsv");
    let contrast = [
        "--contrast=(Intercept)=0",
        "--contrast=typesingle-read=0",
        "--contrast=conditiontreated=1",
    ];
    for (workers, path) in [("1", &bundle), ("4", &parallel)] {
        success(&[
            "--json",
            "--threads",
            workers,
            "run",
            "--counts",
            counts.to_str().unwrap(),
            "--samples",
            samples.to_str().unwrap(),
            "--design",
            design.to_str().unwrap(),
            "--bundle",
            path.to_str().unwrap(),
            contrast[0],
            contrast[1],
            contrast[2],
        ]);
    }
    let analysis = read(&bundle);
    assert_eq!(
        analysis["dataset"]["counts"]["gene_ids"]
            .as_array()
            .unwrap()
            .len(),
        13
    );
    assert_eq!(
        analysis["dataset"]["metadata"][0],
        serde_json::json!(["untreated", "single-read"])
    );
    assert_eq!(
        analysis["normalization"]["kept_indices"],
        serde_json::json!([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    );
    assert_eq!(
        analysis["base"]["residual_df"],
        serde_json::json!([4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2])
    );
    assert_eq!(
        analysis["base"]["zero_mask"][10],
        serde_json::json!([true, true, false, false, true, false, false])
    );
    assert!(analysis["base"]["coefficients"][10][1].as_f64().unwrap() < -3.0);
    assert!(
        analysis["base"]["unscaled_se"][10]
            .as_array()
            .unwrap()
            .iter()
            .all(|v| v.as_f64().unwrap() > 0.0)
    );
    assert_eq!(
        analysis["initial_results"]["moderation_method"],
        "unequal-df"
    );
    let mut parallel = read(&parallel);
    parallel["workers"] = serde_json::json!(1);
    assert_eq!(analysis, parallel);
    let original = std::fs::read(&bundle).unwrap();
    for (name, terms) in [
        ("forward", contrast),
        (
            "reverse",
            [
                "--contrast=(Intercept)=0",
                "--contrast=typesingle-read=0",
                "--contrast=conditiontreated=-1",
            ],
        ),
        (
            "mixed",
            [
                "--contrast=(Intercept)=0",
                "--contrast=typesingle-read=-1",
                "--contrast=conditiontreated=1",
            ],
        ),
    ] {
        let table = directory.path().join(format!("{name}.tsv"));
        let state = directory.path().join(format!("{name}.json"));
        success(&[
            "--json",
            "results",
            "--bundle",
            bundle.to_str().unwrap(),
            "--output",
            table.to_str().unwrap(),
            "--state",
            state.to_str().unwrap(),
            terms[0],
            terms[1],
            terms[2],
        ]);
        let result = read(&state);
        assert_eq!(result["full_base_reused"], true);
        assert_eq!(result["linear_model_refit"], false);
        assert_eq!(result["rows"].as_array().unwrap().len(), 11);
    }
    let forward = read(&directory.path().join("forward.json"));
    assert_eq!(analysis["initial_results"], forward);
    let reverse = read(&directory.path().join("reverse.json"));
    for (a, b) in forward["rows"]
        .as_array()
        .unwrap()
        .iter()
        .zip(reverse["rows"].as_array().unwrap())
    {
        assert_eq!(
            a["log_fc"].as_f64().unwrap(),
            -b["log_fc"].as_f64().unwrap()
        );
        assert_eq!(a["t"].as_f64().unwrap(), -b["t"].as_f64().unwrap());
        for key in ["p_value", "fdr", "lods", "f_statistic"] {
            assert_eq!(a[key], b[key]);
        }
    }
    for kind in ["log-cpm", "weights", "fitted", "residuals"] {
        let path = directory.path().join(format!("{kind}.tsv"));
        success(&[
            "--json",
            "export",
            "--bundle",
            bundle.to_str().unwrap(),
            "--output",
            path.to_str().unwrap(),
            "--kind",
            kind,
        ]);
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(path)
            .unwrap();
        assert_eq!(reader.headers().unwrap().len(), 8);
        assert_eq!(reader.records().count(), 11);
    }
    let table = directory.path().join("forward.tsv");
    let before = std::fs::read(&table).unwrap();
    for (output, state) in [(&bundle, None), (&table, Some(&bundle))] {
        let mut arguments = vec![
            "--json",
            "results",
            "--bundle",
            bundle.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            contrast[0],
            contrast[1],
            contrast[2],
        ];
        if let Some(state) = state {
            arguments.extend(["--state", state.to_str().unwrap()]);
        }
        let result = invoke(&arguments);
        assert!(!result.status.success());
        let error: Value = serde_json::from_slice(&result.stderr).unwrap();
        assert_eq!(error["error"]["kind"], "ConfigError");
    }
    assert_eq!(std::fs::read(table).unwrap(), before);
    let bad = directory.path().join("bad-design.tsv");
    std::fs::write(&bad,"sample\t(Intercept)\ttypesingle-read\tconditiontreated\ns1\t1\t1\t0\ns2\t1\t1\t0\ns3\t1\t1\t0\ns4\t1\t1\t0\ns5\t1\t1\t1\ns6\t1\t1\t1\ns7\t1\t1\t1\n").unwrap();
    let failed_bundle = directory.path().join("failed.json");
    let result = invoke(&[
        "--json",
        "run",
        "--counts",
        counts.to_str().unwrap(),
        "--samples",
        samples.to_str().unwrap(),
        "--design",
        bad.to_str().unwrap(),
        "--bundle",
        failed_bundle.to_str().unwrap(),
        contrast[0],
        contrast[1],
        contrast[2],
    ]);
    assert!(!result.status.success());
    assert!(!failed_bundle.exists());
    assert_eq!(std::fs::read(bundle).unwrap(), original);
}
