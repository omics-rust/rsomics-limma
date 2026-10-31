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

#[test]
fn ordinary_expression_saved_lifecycle_and_tiny_results() {
    let directory = tempfile::tempdir().unwrap();
    let root = directory.path();
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let expression = root.join("expression.tsv");
    let samples = fixture.join("samples.tsv");
    let design = fixture.join("design.tsv");
    std::fs::write(&expression, "gene\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ng1\t-3\t-2\t-1\t-4\t-1\t0\t1\ng2\t4\t5.2\t4.1\t6\t7\t6.9\t8\ng3\t0\t2\t-1\t1\t3\t1.5\t2\ng4\t-8\t-7\t-5\t-6\t-4\t-3\t-2\ng5\t2\t0\t1\t-1\t2.3\t2.9\t3.6\ng6\t9\t8\t10\t11\t8.6\t10\t12\n").unwrap();
    let forward = [
        "--contrast=(Intercept)=0",
        "--contrast=typesingle-read=0",
        "--contrast=conditiontreated=1",
    ];
    let bundle = root.join("ordinary.json");
    let parallel = root.join("parallel.json");
    for (workers, output) in [("1", &bundle), ("4", &parallel)] {
        success(&[
            "--json",
            "--threads",
            workers,
            "run",
            "--expression",
            expression.to_str().unwrap(),
            "--samples",
            samples.to_str().unwrap(),
            "--design",
            design.to_str().unwrap(),
            "--bundle",
            output.to_str().unwrap(),
            forward[0],
            forward[1],
            forward[2],
        ]);
    }
    let analysis = read(&bundle);
    assert_eq!(analysis["schema"], 2);
    assert_eq!(analysis["input_kind"], "log-expression");
    assert_eq!(
        analysis["dataset"]["expression"]["gene_ids"],
        serde_json::json!(["g1", "g2", "g3", "g4", "g5", "g6"])
    );
    assert_eq!(
        analysis["dataset"]["metadata"][0],
        serde_json::json!(["untreated", "single-read"])
    );
    assert_eq!(
        analysis["base"]["residual_df"],
        serde_json::json!([4, 4, 4, 4, 4, 4])
    );
    assert!(analysis.get("normalization").is_none());
    for name in [
        "weights",
        "targets",
        "zero_mask",
        "zero_probes",
        "voom_line",
    ] {
        assert!(analysis["base"].get(name).is_none());
    }
    let mut other = read(&parallel);
    other["workers"] = serde_json::json!(1);
    assert_eq!(analysis, other);
    let original = std::fs::read(&bundle).unwrap();
    for (name, terms) in [
        ("forward", forward),
        (
            "reverse",
            [forward[0], forward[1], "--contrast=conditiontreated=-1"],
        ),
        (
            "mixed",
            [forward[0], "--contrast=typesingle-read=-1", forward[2]],
        ),
    ] {
        let table = root.join(format!("{name}.tsv"));
        let state = root.join(format!("{name}.json"));
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
        assert_eq!(result["rows"].as_array().unwrap().len(), 6);
    }
    let result = read(&root.join("forward.json"));
    assert_eq!(analysis["initial_results"], result);
    let reverse = read(&root.join("reverse.json"));
    for (a, b) in result["rows"]
        .as_array()
        .unwrap()
        .iter()
        .zip(reverse["rows"].as_array().unwrap())
    {
        assert_eq!(
            a["log_fc"].as_f64().unwrap(),
            -b["log_fc"].as_f64().unwrap()
        );
        for key in ["p_value", "fdr", "lods", "f_statistic"] {
            assert_eq!(a[key], b[key]);
        }
    }
    for kind in ["expression", "fitted", "residuals"] {
        let output = root.join(format!("{kind}.tsv"));
        success(&[
            "export",
            "--bundle",
            bundle.to_str().unwrap(),
            "--kind",
            kind,
            "--output",
            output.to_str().unwrap(),
        ]);
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(output)
            .unwrap();
        assert_eq!(
            reader.headers().unwrap().iter().collect::<Vec<_>>(),
            ["gene", "s1", "s2", "s3", "s4", "s5", "s6", "s7"]
        );
        assert_eq!(reader.records().count(), 6);
    }
    for kind in ["log-cpm", "weights"] {
        let output = root.join(format!("invalid-{kind}.tsv"));
        assert!(
            !invoke(&[
                "export",
                "--bundle",
                bundle.to_str().unwrap(),
                "--kind",
                kind,
                "--output",
                output.to_str().unwrap()
            ])
            .status
            .success()
        );
        assert!(!output.exists());
    }
    let counts = fixture.join("counts.tsv");
    let failed = root.join("invalid.json");
    for arguments in [
        vec![
            "--counts",
            counts.to_str().unwrap(),
            "--expression",
            expression.to_str().unwrap(),
        ],
        vec![
            "--expression",
            expression.to_str().unwrap(),
            "--contrast=unknown=1",
        ],
    ] {
        let mut command = vec![
            "run",
            "--samples",
            samples.to_str().unwrap(),
            "--design",
            design.to_str().unwrap(),
            "--bundle",
            failed.to_str().unwrap(),
            forward[0],
            forward[1],
            forward[2],
        ];
        command.extend(arguments);
        assert!(!invoke(&command).status.success());
        assert!(!failed.exists());
    }
    let alias = invoke(&[
        "results",
        "--bundle",
        bundle.to_str().unwrap(),
        "--output",
        bundle.to_str().unwrap(),
        forward[0],
        forward[1],
        forward[2],
    ]);
    assert!(!alias.status.success());
    assert_eq!(std::fs::read(&bundle).unwrap(), original);
    let tiny_design = root.join("tiny-design.tsv");
    let tiny_samples = root.join("tiny-samples.tsv");
    std::fs::write(
        &tiny_design,
        "sample\tIntercept\tgroup\ns4\t1\t1\ns2\t1\t0\ns1\t1\t0\ns3\t1\t1\n",
    )
    .unwrap();
    std::fs::write(
        &tiny_samples,
        "sample\torigin\ns4\tsynthetic\ns3\tsynthetic\ns2\tsynthetic\ns1\tsynthetic\n",
    )
    .unwrap();
    for (name, content, n) in [
        (
            "two",
            "gene\ts1\ts2\ts3\ts4\ng1\t0\t1\t3\t5\ng2\t-3\t-1\t0\t1\n",
            2,
        ),
        ("zero", "gene\ts1\ts2\ts3\ts4\ng1\t0\t0\t0\t0\n", 1),
    ] {
        let input = root.join(format!("{name}-expression.tsv"));
        let saved = root.join(format!("{name}-bundle.json"));
        std::fs::write(&input, content).unwrap();
        success(&[
            "run",
            "--expression",
            input.to_str().unwrap(),
            "--samples",
            tiny_samples.to_str().unwrap(),
            "--design",
            tiny_design.to_str().unwrap(),
            "--bundle",
            saved.to_str().unwrap(),
            "--contrast=Intercept=0",
            "--contrast=group=1",
        ]);
        let fitted = read(&saved);
        assert_eq!(fitted["initial_results"]["prior_df"], 0.0);
        assert_eq!(fitted["initial_results"]["prior_variance"]["shape"], "rows");
        assert_eq!(
            fitted["initial_results"]["rows"].as_array().unwrap().len(),
            n
        );
        let table = root.join(format!("{name}-result.tsv"));
        let state = root.join(format!("{name}-result.json"));
        success(&[
            "results",
            "--bundle",
            saved.to_str().unwrap(),
            "--output",
            table.to_str().unwrap(),
            "--state",
            state.to_str().unwrap(),
            "--contrast=Intercept=0",
            "--contrast=group=1",
        ]);
        assert_eq!(fitted["initial_results"], read(&state));
        if name == "zero" {
            assert_eq!(
                fitted["initial_results"]["coefficient_prior_variance"],
                "Inf"
            );
            assert_eq!(
                fitted["initial_results"]["coefficient_prior_fallback"],
                true
            );
            for key in ["t", "p_value", "fdr", "lods", "f_statistic", "f_p_value"] {
                assert_eq!(fitted["initial_results"]["rows"][0][key], "NaN");
            }
            assert!(
                std::fs::read_to_string(table)
                    .unwrap()
                    .lines()
                    .nth(1)
                    .unwrap()
                    .starts_with("g1\t")
            );
        }
    }
    let boundary_samples = root.join("boundary-samples.tsv");
    let boundary_design = root.join("boundary-design.tsv");
    let boundary_expression = root.join("boundary-expression.tsv");
    std::fs::write(&boundary_samples, "sample\ns1\ns2\ns3\n").unwrap();
    std::fs::write(&boundary_expression, "gene\ts1\ts2\ts3\ng1\t1\t2\t3\n").unwrap();
    let mut boundary_failures = Vec::new();
    for (name, scale) in [("unit", "1"), ("scaled", "1e-100")] {
        std::fs::write(
            &boundary_design,
            format!("sample\tx\ns1\t{scale}\ns2\t0\ns3\t0\n"),
        )
        .unwrap();
        let saved = root.join(format!("{name}-bundle.json"));
        success(&[
            "run",
            "--expression",
            boundary_expression.to_str().unwrap(),
            "--samples",
            boundary_samples.to_str().unwrap(),
            "--design",
            boundary_design.to_str().unwrap(),
            "--bundle",
            saved.to_str().unwrap(),
            "--contrast=x=1",
        ]);
        let fitted = read(&saved);
        let row = &fitted["initial_results"]["rows"][0];
        if !row["t"]
            .as_f64()
            .is_some_and(|value| (value - 0.3922322702763681).abs() < 1e-12)
            || !row["p_value"]
                .as_f64()
                .is_some_and(|value| (value - 0.7327387580875756).abs() < 1e-12)
        {
            boundary_failures.push(format!("{name} scaling changed finite t/P: {row}"));
        }
        let table = root.join(format!("{name}-results.tsv"));
        let state = root.join(format!("{name}-results.json"));
        success(&[
            "results",
            "--bundle",
            saved.to_str().unwrap(),
            "--output",
            table.to_str().unwrap(),
            "--state",
            state.to_str().unwrap(),
            "--contrast=x=1",
        ]);
        assert_eq!(fitted["initial_results"], read(&state));
    }
    std::fs::write(&boundary_expression, "gene\ts1\ts2\ts3\ng1\t1e308\t0\t0\n").unwrap();
    std::fs::write(&boundary_design, "sample\tx\ns1\t1\ns2\t0\ns3\t0\n").unwrap();
    let overflow_bundle = root.join("overflow-bundle.json");
    success(&[
        "run",
        "--expression",
        boundary_expression.to_str().unwrap(),
        "--samples",
        boundary_samples.to_str().unwrap(),
        "--design",
        boundary_design.to_str().unwrap(),
        "--bundle",
        overflow_bundle.to_str().unwrap(),
        "--contrast=x=1",
    ]);
    let before = std::fs::read(&overflow_bundle).unwrap();
    let output = root.join("overflow-results.tsv");
    let state = root.join("overflow-results.json");
    std::fs::write(&output, "table sentinel\n").unwrap();
    std::fs::write(&state, "state sentinel\n").unwrap();
    for arguments in [
        vec![
            "run",
            "--expression",
            boundary_expression.to_str().unwrap(),
            "--samples",
            boundary_samples.to_str().unwrap(),
            "--design",
            boundary_design.to_str().unwrap(),
            "--bundle",
            overflow_bundle.to_str().unwrap(),
            "--contrast=x=2",
        ],
        vec![
            "results",
            "--bundle",
            overflow_bundle.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--state",
            state.to_str().unwrap(),
            "--contrast=x=2",
        ],
    ] {
        std::fs::write(&overflow_bundle, &before).unwrap();
        let failed = invoke(&arguments);
        if failed.status.success() {
            boundary_failures.push(format!("{} accepted effect overflow", arguments[0]));
        }
        if std::fs::read(&overflow_bundle).unwrap() != before {
            boundary_failures.push(format!(
                "{} replaced the bundle after overflow",
                arguments[0]
            ));
        }
    }
    if std::fs::read(output).unwrap() != b"table sentinel\n"
        || std::fs::read(state).unwrap() != b"state sentinel\n"
    {
        boundary_failures.push("overflow replaced existing outputs".to_string());
    }
    assert!(
        boundary_failures.is_empty(),
        "{}",
        boundary_failures.join("\n")
    );
}
