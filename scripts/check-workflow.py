#!/usr/bin/env python3
"""Check a complete native current-voom lifecycle against the frozen source witness."""

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

SOURCE_RUN = 40606516204
MANIFEST_SHA = "af682fed901644550411e6ae74e22b6bda45a7e567ca15e2dddc2c6b17027e97"
WITNESS_SHA = "091f476ca726487d23e83de7d62b5bffa7cea9021b533a4dda2043c89289a577"
BOUNDS = dict(tmm=(1e-10, 1e-10), library=(1e-8, 1e-10), expression=(1e-10, 1e-10),
              linear=(1e-7, 1e-6), trend=(1e-7, 1e-6), weights=(1e-6, 1e-5),
              eb=(1e-6, 1e-4), df=(1e-3, 1e-4), odds=(1e-5, 5e-4), probability=(1e-12, 5e-4))
LOG_PROBABILITY_BOUND = 1e-3


def sha(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            result.update(block)
    return result.hexdigest()


def save(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def table(path):
    with path.open(newline="") as stream:
        return list(csv.reader(stream, delimiter="\t"))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def numeric(value):
    require(value is not None, "unexplained numeric JSON null")
    if isinstance(value, str):
        value = value.strip()
    if value == "NA":
        return "NA", None
    value = float(value)
    if math.isnan(value):
        return "NaN", None
    if math.isinf(value):
        return "Inf" if value > 0 else "-Inf", None
    if value == 0:
        return "zero", value
    return "finite", value


def compare(binary, witness, output, workers, reuse):
    require(sha(witness/"sha256-manifest.tsv") == MANIFEST_SHA and sha(witness/"witness.json") == WITNESS_SHA,
            "accepted source witness identity differs")
    binary = binary.resolve(strict=True)
    output.mkdir(parents=True, exist_ok=True)
    contract = dict(source_run=SOURCE_RUN, manifest_sha256=MANIFEST_SHA, witness_sha256=WITNESS_SHA,
                    binary_sha256=sha(binary), workers=workers, bounds=BOUNDS,
                    positive_probability_log_bound=LOG_PROBABILITY_BOUND, reuse=reuse)
    save(output/"comparison-contract.json", contract)
    contrasts = table(witness/"inputs/contrasts.tsv")
    names = [row[0] for row in contrasts[1:]]
    terms = {name: [f"--contrast={row[0]}={float(row[j])}" for row in contrasts[1:]]
             for j, name in enumerate(contrasts[0][1:], 1)}
    commands = []

    def run(arguments):
        index = len(commands)
        argv = [str(binary), "--json", "--threads", str(workers), *map(str, arguments)]
        started = time.perf_counter()
        record = dict(argv=argv)
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
        require(record["exit"] == 0, f"native command failed: {arguments[0]}; retained command-{index}.stderr")

    bundle = output/"fit.json"
    try:
        if not reuse:
            require(not bundle.exists(), "fresh output already has a fitted bundle")
            run(["run", "--counts", witness/"inputs/counts.tsv", "--samples", witness/"inputs/samples.tsv",
                 "--design", witness/"inputs/model-matrix.tsv", "--bundle", bundle, *terms["forward"]])
        require(bundle.is_file(), "saved native full base fit is absent")
        bundle_before = sha(bundle)
        for name in terms:
            run(["results", "--bundle", bundle, "--output", output/f"{name}.tsv", "--state", output/f"{name}.json", *terms[name]])
        for kind, name in [("log-cpm", "logCPM"), ("weights", "weights"), ("fitted", "fitted"), ("residuals", "residuals")]:
            run(["export", "--bundle", bundle, "--output", output/f"{name}.tsv", "--kind", kind])
        require(sha(bundle) == bundle_before, "reload mutated the original fitted bundle bytes")
    finally:
        save(output/"commands.json", commands)

    analysis = json.loads(bundle.read_text())
    require(analysis["schema"] == 1 and analysis["product"] == "rsomics-limma", "wrong analysis product/schema")
    data, lib, base = analysis["dataset"], analysis["normalization"], analysis["base"]
    raw = table(witness/"inputs/counts.tsv")
    raw_ids = [row[0] for row in raw[1:]]
    require(data["counts"]["gene_ids"] == raw_ids and data["counts"]["sample_ids"] == raw[0][1:], "raw identities/order differ")
    require(data["counts"]["values"] == [float(value) for row in raw[1:] for value in row[1:]], "raw count values differ")
    for key, name in [("counts", "counts.tsv"), ("samples", "samples.tsv"), ("design", "model-matrix.tsv")]:
        require(data["input_sha256"][key] == sha(witness/"inputs"/name), "input hash differs")
    metadata = table(witness/"inputs/samples.tsv")
    by_sample = {row[0]: row[1:] for row in metadata[1:]}
    require(data["metadata_columns"] == metadata[0][1:] and data["metadata"] == [by_sample[sample] for sample in raw[0][1:]], "joined metadata differs")
    design = table(witness/"inputs/model-matrix.tsv")
    by_sample = {row[0]: row[1:] for row in design[1:]}
    require(data["design"]["names"] == names and data["design"]["n_samples"] == len(raw[0])-1 and
            data["design"]["values"] == [float(value) for sample in raw[0][1:] for value in by_sample[sample]], "raw named design differs")
    require(base["design"] == data["design"], "saved full base design differs")
    mask = table(witness/"fit/filter-mask.tsv")
    require([row[0] for row in mask[1:]] == raw_ids, "source mask identities differ")
    require(lib["keep"] == [row[1] == "TRUE" for row in mask[1:]], "filter decisions differ")
    require(lib["kept_indices"] == [i for i, keep in enumerate(lib["keep"]) if keep], "retained mapping differs")
    ids = [raw_ids[i] for i in lib["kept_indices"]]
    n = len(raw[0])-1
    require(len(base["coefficients"]) == len(ids), "fitted population differs")
    kept = table(witness/"fit/kept-counts.tsv")
    require([row[0] for row in kept[1:]] == ids and [[float(x) for x in row[1:]] for row in kept[1:]] ==
            [data["counts"]["values"][i*n:(i+1)*n] for i in lib["kept_indices"]], "kept count identities/values differ")
    descriptor = json.loads((witness/"fit/base/state-descriptors.json").read_text())
    require(base["rank"] == descriptor["rank"] and base["pivot"] == descriptor["pivot"] and base["method"] == descriptor["method"], "rank/pivot/method differs")
    source_parameters = json.loads((witness/"inputs/public-call-parameters.json").read_text())
    for section in ["construction", "subset", "normalization", "moderation"]:
        require(analysis["parameters"][section] == source_parameters[section], "configured defaults differ: "+section)
    for section in ["filtering", "voom", "topTable"]:
        for field, value in analysis["parameters"][section].items():
            if field == "derived_span_rule":
                continue
            require(source_parameters[section][field] == value, "configured option differs: "+section+":"+field)

    failures, metrics = [], {}

    def point(field, actual, expected, kind, identity=None):
        aclass, a = numeric(actual)
        eclass, e = numeric(expected)
        if aclass != eclass:
            failures.append(dict(field=field, identity=identity, actual_class=aclass, expected_class=eclass))
            return
        if aclass not in ("finite", "zero"):
            return
        atol, rtol = BOUNDS[kind]
        scaled = abs(a-e)/(atol+rtol*abs(e))
        log_error = abs(math.log(a/e)) if kind == "probability" and a > 0 and e > 0 else 0.0
        metric = metrics.setdefault(field, dict(count=0, failed=0, max_abs=0.0, max_scaled=0.0, max_log_probability_error=0.0))
        metric["count"] += 1
        metric["max_abs"] = max(metric["max_abs"], abs(a-e))
        metric["max_scaled"] = max(metric["max_scaled"], scaled)
        metric["max_log_probability_error"] = max(metric["max_log_probability_error"], log_error)
        if scaled > 1 or log_error > LOG_PROBABILITY_BOUND:
            metric["failed"] += 1
            failures.append(dict(field=field, identity=identity, actual=a, expected=e, scaled=scaled, log_probability_error=log_error))

    def vectors(name, actual, kinds, identities=ids):
        reference = table(witness/name)
        require([row[0] for row in reference[1:]] == identities and len(actual) == len(identities), "identity order differs: "+name)
        require(all(len(row) == len(reference[0])-1 for row in actual), "matrix shape differs: "+name)
        for i, (row, expected) in enumerate(zip(actual, reference[1:])):
            for column, (field, kind) in enumerate(zip(reference[0][1:], kinds)):
                point(name+":"+field, row[column], expected[column+1], kind, identities[i])

    vectors("fit/library-state.tsv", list(zip(lib["raw_library_sizes"], lib["postfilter_library_sizes"], lib["tmm_factors"], lib["effective_library_sizes"], lib["offsets"])),
            ["library", "library", "tmm", "library", "tmm"], raw[0][1:])
    targets = table(witness/"fit/targets.tsv")
    require([row[0] for row in targets[1:]] == raw[0][1:] and len(base["targets"]) == n, "sample targets differ")
    for i, (actual, expected) in enumerate(zip(base["targets"], targets[1:])):
        require(actual["sample"] == expected[0] and actual["group"] == expected[1] and actual["metadata"] == expected[4:], "target identity/metadata differs")
        point("targets:lib.size", actual["lib_size"], expected[2], "library", actual["sample"])
        point("targets:norm.factors", actual["norm_factor"], expected[3], "tmm", actual["sample"])
    vectors("fit/base/coefficients.tsv", base["coefficients"], ["linear"]*len(names))
    vectors("fit/base/stdev.unscaled.tsv", base["unscaled_se"], ["linear"]*len(names))
    vectors("fit/base/cov.coefficients.tsv", base["covariance"], ["linear"]*len(names), names)
    vectors("fit/base/row-diagnostics.tsv", list(zip(base["amean"], base["sigma"], base["residual_df"])), ["expression", "trend", "df"])
    diagnostics = table(witness/"fit/base/row-diagnostics.tsv")
    require(base["residual_df"] == [int(float(row[-1])) for row in diagnostics[1:]], "per-gene structural zero residual df differs")
    semantics = json.loads((witness/"fit/state-semantics.json").read_text())
    point("derived_span", base["derived_span"], semantics["derived_span"], "tmm")
    points = table(witness/"fit/voom-points.tsv")
    vectors("fit/voom-points.tsv", base["voom_points"], ["expression", "trend"], [row[0] for row in points[1:]])
    line = table(witness/"fit/voom-line.tsv")
    vectors("fit/voom-line.tsv", list(zip(base["voom_line"]["x"], base["voom_line"]["y"])), ["expression", "trend"], [row[0] for row in line[1:]])
    xs = [float(row[1]) for row in line[1:]]
    cumulative = []
    total = 0.0
    for gap in sorted(b-a for a, b in zip(xs, xs[1:])):
        total += gap
        cumulative.append(total)
    delta = 0.0 if len(xs) <= 200 else min(cumulative[len(xs)-2-k]/(200-k) for k in range(200))
    point("weighted_LOWESS_delta", base["voom_line"]["delta"], delta, "tmm")
    require(base["trend_policy"] == "weighted-lowess/frequency-df/four-passes", "wrong witnessed voom trend route")

    discoveries = {}

    def result_state(label, result, source_label):
        require(result["full_base_reused"] and not result["linear_model_refit"], "reload derivation semantics differ")
        require(result["contrast"] == [float(row[contrasts[0].index(source_label)]) for row in contrasts[1:]], "named contrast weights differ")
        require(result["base_rank"] == base["rank"] and result["base_pivot"] == base["pivot"], "contrast lost base rank/pivot")
        require(result["moderation_method"] == "unequal-df", "wrong default variance-prior route")
        source_base = "results/"+source_label
        vectors(source_base+"-contrast/coefficients.tsv", [[row["log_fc"]] for row in result["rows"]], ["linear"])
        vectors(source_base+"-contrast/stdev.unscaled.tsv", [[row["unscaled_se"]] for row in result["rows"]], ["linear"])
        source_cov = table(witness/(source_base+"-contrast/cov.coefficients.tsv"))
        require(len(result["covariance"]) == 1 and len(result["covariance"][0]) == 1, "contrast covariance shape differs")
        point(label+":covariance", result["covariance"][0][0], source_cov[1][1], "linear")
        vectors(source_base+"-contrast/row-diagnostics.tsv", list(zip(base["amean"], base["sigma"], base["residual_df"])), ["expression", "trend", "df"])
        moderated = source_base+"-moderated/"
        desc = json.loads((witness/(moderated+"state-descriptors.json")).read_text())
        require(desc["shapes"]["df.prior"]["length"] == 1 and desc["shapes"]["s2.prior"]["length"] == 1 and result["prior_variance"]["shape"] == "scalar", "prior scalar/row shapes differ")
        point(label+":var.prior", result["coefficient_prior_variance"], desc["var.prior"], "odds")
        require(result["proportion"] == desc["proportion"], "DE proportion differs")
        require(desc["rank"] == result["base_rank"] and desc["pivot"] == result["base_pivot"], "moderated base rank/pivot differs")
        for field, key, kind in [("coefficients", "log_fc", "linear"), ("stdev.unscaled", "unscaled_se", "linear"), ("t", "t", "eb"), ("p.value", "p_value", "probability"), ("lods", "lods", "odds")]:
            vectors(moderated+field+".tsv", [[row[key]] for row in result["rows"]], [kind])
        scale = result["prior_variance"]["values"]
        values = [[base["amean"][i], base["sigma"][i], base["residual_df"][i], scale, result["prior_df"], row["posterior_variance"], row["df_total"], row["f_statistic"], row["f_p_value"]] for i, row in enumerate(result["rows"])]
        vectors(moderated+"row-diagnostics.tsv", values, ["expression", "trend", "df", "eb", "df", "eb", "df", "eb", "probability"])
        vectors(source_base+".tsv", [[row[key] for key in ["log_fc", "ave_expr", "t", "p_value", "fdr", "lods"]] for row in result["rows"]], ["linear", "expression", "eb", "probability", "probability", "odds"])

    result_state("initial-forward", analysis["initial_results"], "forward")
    for name in terms:
        state = json.loads((output/f"{name}.json").read_text())
        result_state("reloaded-"+name, state, name)
        if name == "forward":
            require(state == analysis["initial_results"], "initial/reloaded complete forward contrast state differs")
        actual, reference = table(output/f"{name}.tsv"), table(witness/f"results/{name}.tsv")
        require(actual[0] == reference[0] and [row[0] for row in actual[1:]] == ids, "reloaded result names/order differ")
        require([[float(x) for x in row[1:]] for row in actual[1:]] == [[row[key] for key in ["log_fc", "ave_expr", "t", "p_value", "fdr", "lods"]] for row in state["rows"]], "TSV differs from complete result state")
        significant = {row[0] for row in actual[1:] if float(row[-2]) <= .05}
        expected = {row[0] for row in reference[1:] if float(row[-2]) <= .05}
        if significant != expected:
            failures.append(dict(field=name+" discovery identities", added=sorted(significant-expected), missing=sorted(expected-significant)))
        discoveries[name] = dict(actual=len(significant), source=len(expected))
    for name, kind in [("logCPM", "expression"), ("weights", "weights"), ("fitted", "linear"), ("residuals", "linear")]:
        actual = table(output/f"{name}.tsv")
        require(actual[0][1:] == raw[0][1:] and [row[0] for row in actual[1:]] == ids, "observation/model export identities differ")
        values = [[float(value) for value in row[1:]] for row in actual[1:]]
        vectors(f"results/{name}.tsv", values, [kind]*n)
        if name in ("logCPM", "weights"):
            require(values == base["expression" if name == "logCPM" else "weights"], "export differs from stored observation state")

    save(output/"comparison.json", dict(passed=not failures, failure_count=len(failures), bounds=BOUNDS,
         metrics=metrics, discoveries=discoveries, raw_genes=len(raw_ids), kept_genes=len(ids),
         residual_df_distribution=dict(Counter(base["residual_df"])), native_bundle_sha256=bundle_before,
         source_run=SOURCE_RUN, source_manifest_sha256=MANIFEST_SHA, source_witness_sha256=WITNESS_SHA))
    save(output/"failures.json", failures)
    save(output/"export-sha256.json", {path.name: sha(path) for path in sorted(output.glob("*.tsv"))})
    require(not failures, f"unchanged whole current-voom comparison failed in {len(failures)} points; all failures retained")
    print(json.dumps(dict(passed=True, kept_genes=len(ids), source_run=SOURCE_RUN, discoveries=discoveries)), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--witness", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, choices=(1, 4), required=True)
    parser.add_argument("--reuse", action="store_true")
    args = parser.parse_args()
    compare(args.binary, args.witness, args.output, args.workers, args.reuse)
