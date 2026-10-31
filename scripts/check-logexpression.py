#!/usr/bin/env python3
"""Compare complete ordinary native lifecycles with the unchanged adopted source witness."""

import argparse
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import time

_spec = importlib.util.spec_from_file_location("limma_count_contract", Path(__file__).with_name("check-workflow.py"))
_shared = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_shared)
sha, table, require, numeric = _shared.sha, _shared.table, _shared.require, _shared.numeric

SOURCE_RUN = 41184814043
MANIFEST = "c374d2be76ce0d0ee25aae2eeb14548ddb02d9866a0266e814bd39c5dec5c19f"
CAPTURE = "ee4a903ce605ac51d059469d738020fb31e02ac45332156ce30ead5abaece4ab"
RELOAD = "87ca32ca6d0a199c7b6c729ae99096242fcc223e5ef8eee94bde3e5ae922af5f"
SOURCE_COMMIT = "825d1c8330d024adb6c54af66d79b0f4ec258530"
SOURCE_ARCHIVE = "40b1d2059a76ad98db4bf80139760acf8e8bae92a8265aaaf16ac199831533f4"
BOUNDS = dict(expression=(1e-10, 1e-10), linear=(1e-7, 1e-6), eb=(1e-6, 1e-4),
              df=(1e-3, 1e-4), odds=(1e-5, 5e-4), probability=(1e-12, 5e-4))
LOG_BOUND = 1e-3
CASES = ("real", "equalvar", "constvar", "two-row", "one-zero")


def save(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def qualify(witness):
    require(sha(witness/"sha256-manifest.tsv") == MANIFEST, "source manifest differs")
    require(sha(witness/"witness.json") == CAPTURE, "source capture witness differs")
    require(sha(witness/"reload-process/witness.json") == RELOAD, "source reload witness differs")
    manifest = table(witness/"sha256-manifest.tsv")
    paths = [row[1] for row in manifest[1:]]
    require(len(paths) == 827 and len(set(paths)) == 827, "source membership count differs")
    require(set(paths) == {str(path.relative_to(witness)) for path in witness.rglob("*")
                           if path.is_file() and path.name != "sha256-manifest.tsv"}, "source membership differs")
    for digest, name in manifest[1:]:
        path = witness/name
        require(not path.is_symlink() and path.resolve().is_relative_to(witness.resolve()) and sha(path) == digest,
                "source member differs: " + name)
    capture = json.loads((witness/"witness.json").read_text())
    reload = json.loads((witness/"reload-process/witness.json").read_text())
    require(capture["status"] == "captured-awaiting-independent-reload" and capture["failure"] is None
            and reload["status"] == "complete", "source capture/reload status differs")
    pin = capture["identity"]["pin"]
    require(pin["version"] == "3.68.5" and pin["source_commit"] == SOURCE_COMMIT
            and pin["source_archive_sha256"] == SOURCE_ARCHIVE, "released source pin differs")
    for name in CASES:
        require(capture["cases"][name]["complete_derivation"] and all(capture["cases"][name]["checks"].values())
                and all(reload["checks"][name].values()), "source case qualification failed: " + name)


def compare_case(binary, witness, output, case, workers):
    output.mkdir(parents=True)
    source = witness/"cases"/case
    inputs = witness/"inputs"/case
    raw = table(inputs/"expression.tsv")
    ids, samples = [row[0] for row in raw[1:]], raw[0][1:]
    contrasts = table(inputs/"contrasts.tsv")
    names = [row[0] for row in contrasts[1:]]
    terms = {name: [f"--contrast={row[0]}={float(row[index])}" for row in contrasts[1:]]
             for index, name in enumerate(contrasts[0][1:], 1)}
    commands, failures, metrics = [], [], {}

    def run(arguments):
        index = len(commands)
        argv = [str(binary), "--json", "--threads", str(workers), *map(str, arguments)]
        record = dict(argv=argv)
        started = time.perf_counter()
        try:
            with (output/f"command-{index}.stdout").open("xb") as stdout, (output/f"command-{index}.stderr").open("xb") as stderr:
                process = subprocess.run(argv, stdout=stdout, stderr=stderr, timeout=1800, check=False)
            record["exit"] = process.returncode
        except subprocess.TimeoutExpired:
            record["timeout"] = True
            raise
        finally:
            record["seconds"] = time.perf_counter()-started
            commands.append(record)
        require(record["exit"] == 0, f"native command failed: {arguments[0]}; command-{index}.stderr retained")

    def point(field, actual, expected, kind, identity=None):
        aclass, a = numeric(actual)
        eclass, e = numeric(expected)
        finite = aclass in ("finite", "zero") and eclass in ("finite", "zero")
        if aclass != eclass and not finite:
            failures.append(dict(field=field, identity=identity, actual_class=aclass, expected_class=eclass))
            return
        if aclass not in ("finite", "zero"):
            return
        if kind == "probability" and e > 0 and a == 0:
            failures.append(dict(field=field, identity=identity, actual=a, expected=e,
                                 error="positive source probability became zero"))
            return
        atol, rtol = BOUNDS[kind]
        scaled = abs(a-e)/(atol+rtol*abs(e))
        log_error = abs(math.log(a)-math.log(e)) if kind == "probability" and a > 0 and e > 0 else 0.0
        metric = metrics.setdefault(field, dict(count=0, failed=0, max_abs=0.0, max_scaled=0.0, max_log_error=0.0))
        metric["count"] += 1
        metric["max_abs"] = max(metric["max_abs"], abs(a-e))
        metric["max_scaled"] = max(metric["max_scaled"], scaled)
        metric["max_log_error"] = max(metric["max_log_error"], log_error)
        if scaled > 1 or log_error > LOG_BOUND:
            metric["failed"] += 1
            failures.append(dict(field=field, identity=identity, actual=a, expected=e, scaled=scaled, log_error=log_error))

    def vectors(relative, values, kinds, mapping=None, identities=None):
        reference = table(source/relative)
        require(len(values) == len(reference)-1 and all(len(row) == len(reference[0])-1 for row in values),
                "matrix dimensions differ: " + relative)
        if mapping:
            mapped = table(source/mapping)
            require(mapped[0] == ["source_row", "input_feature_id", "source_feature_id"]
                    and [row[0] for row in mapped[1:]] == [str(i+1) for i in range(len(ids))]
                    and [row[1] for row in mapped[1:]] == ids, "source positional mapping differs: " + mapping)
            expected_ids = [str(i+1) if row[2] == "NA" else row[2] for i, row in enumerate(mapped[1:])]
        else:
            expected_ids = identities
        require([row[0] for row in reference[1:]] == expected_ids, "source IDs differ: " + relative)
        for i, (actual, expected) in enumerate(zip(values, reference[1:])):
            for column, kind in enumerate(kinds):
                point(relative+":"+reference[0][column+1], actual[column], expected[column+1], kind,
                      expected_ids[i])

    try:
        bundle = output/"fit.json"
        first = next(iter(terms))
        run(["run", "--expression", inputs/"expression.tsv", "--samples", inputs/"samples.tsv",
             "--design", inputs/"design.tsv", "--bundle", bundle, *terms[first]])
        before = sha(bundle)
        for name in terms:
            run(["results", "--bundle", bundle, "--output", output/f"{name}.tsv", "--state", output/f"{name}.json", *terms[name]])
        for kind in ("expression", "fitted", "residuals"):
            run(["export", "--bundle", bundle, "--output", output/f"{kind}.tsv", "--kind", kind])
        require(sha(bundle) == before, "saved bundle bytes changed on reload")
        analysis = json.loads(bundle.read_text())
        require(analysis["schema"] == 2 and analysis["product"] == "rsomics-limma"
                and analysis["input_kind"] == "log-expression", "wrong ordinary bundle schema/input kind")
        require(analysis["source_provenance"]["source_run"] == SOURCE_RUN
                and analysis["source_provenance"]["source_manifest_sha256"] == MANIFEST, "ordinary provenance differs")
        data, base = analysis["dataset"], analysis["base"]
        expression = data["expression"]
        require(expression["gene_ids"] == ids and expression["sample_ids"] == samples
                and expression["values"] == [float(value) for row in raw[1:] for value in row[1:]], "supplied E changed")
        for key in ("expression", "samples", "design"):
            require(data["input_sha256"][key] == sha(inputs/f"{key}.tsv"), "input hash differs: " + key)
        for filename, field in [("samples", "metadata"), ("design", "design")]:
            reference = table(inputs/f"{filename}.tsv")
            joined = {row[0]: row[1:] for row in reference[1:]}
            if field == "metadata":
                require(data["metadata_columns"] == reference[0][1:] and data[field] == [joined[sample] for sample in samples], "metadata join differs")
            else:
                require(data[field]["names"] == names and data[field]["n_samples"] == len(samples)
                        and data[field]["values"] == [float(value) for sample in samples for value in joined[sample]], "named design differs")
        require(base["design"] == data["design"] and not any(key in base for key in
                ("weights", "targets", "zero_mask", "zero_probes", "voom_line", "derived_span"))
                and "normalization" not in analysis, "ordinary fit contains count state")
        descriptor = json.loads((source/"base/state-descriptors.json").read_text())
        require(base["rank"] == descriptor["rank"] and base["pivot"] == descriptor["pivot"]
                and base["method"] == descriptor["method"], "base rank/pivot/method differs")
        parameters = json.loads((witness/"metadata/public-call-parameters.json").read_text())
        for section in ("lmFit", "moderation", "topTable"):
            for key, value in parameters[section].items():
                require(analysis["parameters"][section][key] == value, "source default differs: " + section+":"+key)
        for field, key in [("coefficients", "coefficients"), ("stdev.unscaled", "unscaled_se")]:
            vectors(f"base/{field}.tsv", base[key], ["linear"]*len(names), "base/features-input-mapping.tsv")
        vectors("base/cov.coefficients.tsv", base["covariance"], ["linear"]*len(names), identities=names)
        vectors("base/row-diagnostics.tsv", list(zip(base["amean"], base["sigma"], base["residual_df"])),
                ["expression", "linear", "df"], "base/features-input-mapping.tsv")
        require(base["residual_df"] == [len(samples)-len(names)]*len(ids), "ordinary row residual df differs")
        discoveries = {}
        for name in terms:
            state = json.loads((output/f"{name}.json").read_text())
            require(state["full_base_reused"] and not state["linear_model_refit"]
                    and state["base_rank"] == base["rank"] and state["base_pivot"] == base["pivot"], "result refit/rank semantics differ")
            require(state["contrast"] == [float(row[contrasts[0].index(name)]) for row in contrasts[1:]], "complete named contrast differs")
            if name == first:
                require(state == analysis["initial_results"], "initial/reloaded full contrast state differs")
            rows = state["rows"]
            prefix = f"results/{name}"
            for field, key in [("coefficients", "log_fc"), ("stdev.unscaled", "unscaled_se")]:
                vectors(prefix+f"-contrast/{field}.tsv", [[row[key]] for row in rows], ["linear"], prefix+"-contrast/features-input-mapping.tsv")
            vectors(prefix+"-contrast/cov.coefficients.tsv", state["covariance"], ["linear"], identities=[name])
            vectors(prefix+"-contrast/row-diagnostics.tsv", list(zip(base["amean"], base["sigma"], base["residual_df"])),
                    ["expression", "linear", "df"], prefix+"-contrast/features-input-mapping.tsv")
            moderated = prefix+"-moderated/"
            desc = json.loads((source/(moderated+"state-descriptors.json")).read_text())
            require(desc["shapes"]["df.prior"]["length"] == 1 and desc["shapes"]["s2.prior"]["length"] in (1, len(ids)), "source prior shapes differ")
            prior = state["prior_variance"]
            require(prior["shape"] == ("rows" if len(ids) < 3 else "scalar"), "native prior shape differs")
            scales = prior["values"] if prior["shape"] == "rows" else [prior["values"]]*len(ids)
            require(len(scales) == len(ids), "native prior population differs")
            point(name+":coefficient prior", state["coefficient_prior_variance"], desc["var.prior"], "odds")
            require(state["proportion"] == desc["proportion"], "coefficient mixture proportion differs")
            expected_method = "no-shrink-small-population" if len(ids) < 3 else "legacy-equal-positive-df"
            conditions = table(witness/"metadata/conditions.tsv")
            warned = any(row[1] == f"{case} initial eBayes {name}" and row[2] == "warning"
                         and row[4] == "Estimation of var.prior failed - set to default value" for row in conditions[1:])
            require(state["moderation_method"] == expected_method and state["coefficient_prior_fallback"] == warned, "moderation branch/fallback differs")
            for field, key, kind in [("coefficients", "log_fc", "linear"), ("stdev.unscaled", "unscaled_se", "linear"),
                                     ("t", "t", "eb"), ("p.value", "p_value", "probability"), ("lods", "lods", "odds")]:
                vectors(moderated+field+".tsv", [[row[key]] for row in rows], [kind], moderated+"features-input-mapping.tsv")
            values = [[base["amean"][i], base["sigma"][i], base["residual_df"][i], scales[i], state["prior_df"], row["posterior_variance"],
                       row["df_total"], row["f_statistic"], row["f_p_value"]] for i, row in enumerate(rows)]
            vectors(moderated+"row-diagnostics.tsv", values, ["expression", "linear", "df", "eb", "df", "eb", "df", "eb", "probability"], moderated+"features-input-mapping.tsv")
            for i, row in enumerate(rows):
                point(name+":uncapped F df", row["f_denominator_df"], float(state["prior_df"])+base["residual_df"][i], "df", ids[i])
            fields = ("log_fc", "ave_expr", "t", "p_value", "fdr", "lods")
            vectors(prefix+".tsv", [[row[key] for key in fields] for row in rows],
                    ["linear", "expression", "eb", "probability", "probability", "odds"], prefix+"-input-mapping.tsv")
            native = table(output/f"{name}.tsv")
            require(native[0] == ["gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"]
                    and [row[0] for row in native[1:]] == ids, "native stable table IDs differ")
            for actual, row in zip(native[1:], rows):
                require(all(numeric(value) == numeric(row[key]) for value, key in zip(actual[1:], fields)), "table differs from complete native state")
            reference = table(source/(prefix+".tsv"))
            actual_ids = [ids[i] for i, row in enumerate(native[1:]) if float(row[-2]) <= .05]
            expected_ids = [ids[i] for i, row in enumerate(reference[1:]) if float(row[-2]) <= .05]
            if actual_ids != expected_ids:
                failures.append(dict(field=name+":FDR discovery identities", added=sorted(set(actual_ids)-set(expected_ids)),
                                     missing=sorted(set(expected_ids)-set(actual_ids))))
            discoveries[name] = len(actual_ids)
        for kind, filename, band in [("expression", "supplied_expression", "expression"), ("fitted", "fitted", "linear"), ("residuals", "residuals", "linear")]:
            native = table(output/f"{kind}.tsv")
            require(native[0] == ["gene", *samples] and [row[0] for row in native[1:]] == ids, "native observation export IDs differ")
            values = [row[1:] for row in native[1:]]
            vectors(f"results/{filename}.tsv", values, [band]*len(samples), f"results/{filename}-input-mapping.tsv")
            if kind == "expression":
                require([[float(value) for value in row] for row in values] == base["expression"], "stored/supplied E export differs")
        save(output/"comparison.json", dict(passed=not failures, case=case, workers=workers, metrics=metrics, discoveries=discoveries,
             native_bundle_sha256=before, source_run=SOURCE_RUN))
        return analysis, discoveries
    finally:
        save(output/"commands.json", commands)
        save(output/"failures.json", failures)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--witness", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    binary, witness = args.binary.resolve(strict=True), args.witness.resolve(strict=True)
    require(not args.output.exists(), "output already exists")
    args.output.mkdir(parents=True)
    save(args.output/"contract.json", dict(source_run=SOURCE_RUN, manifest_sha256=MANIFEST, capture_sha256=CAPTURE,
        reload_sha256=RELOAD, source_commit=SOURCE_COMMIT, source_archive_sha256=SOURCE_ARCHIVE,
        binary_sha256=sha(binary), controller_sha256=sha(Path(__file__)), shared_controller_sha256=sha(Path(__file__).with_name("check-workflow.py")),
        bounds=BOUNDS, positive_probability_log_bound=LOG_BOUND, workers=[1, 4], cases=CASES))
    qualify(witness)
    summaries, errors = [], []
    for case in CASES:
        states = []
        for workers in (1, 4):
            output = args.output/case/f"workers-{workers}"
            try:
                analysis, discoveries = compare_case(binary, witness, output, case, workers)
                comparison = json.loads((output/"comparison.json").read_text())
                summaries.append(dict(case=case, workers=workers, passed=comparison["passed"], discoveries=discoveries))
                analysis["workers"] = 1
                states.append(analysis)
            except Exception as error:
                errors.append(dict(case=case, workers=workers, error=type(error).__name__, message=str(error)))
        if len(states) == 2:
            if states[0] != states[1]:
                errors.append(dict(case=case, error="worker full bundle values differ"))
            one, four = args.output/case/"workers-1", args.output/case/"workers-4"
            for path in one.iterdir():
                if path.suffix == ".tsv" or path.suffix == ".json" and path.stem in ("forward", "reverse", "mixed"):
                    if sha(path) != sha(four/path.name):
                        errors.append(dict(case=case, error="worker complete result/export bytes differ", path=path.name))
    passed = not errors and len(summaries) == 10 and all(value["passed"] for value in summaries)
    save(args.output/"summary.json", dict(passed=passed, comparisons=summaries, errors=errors, source_run=SOURCE_RUN))
    require(passed, "ordinary source comparison failed; all original outputs/commands/mismatches retained")
    print(json.dumps(dict(passed=True, source_run=SOURCE_RUN, comparisons=len(summaries))))


if __name__ == "__main__":
    main()
