# rsomics-limma

Private, unpublished implementation of the current count-input voom weighted
linear-model lifecycle. Source compatibility, native platforms, performance
and release qualification are separate gates; `publish=false`.

The supported interface is the CLI, versioned analysis/contrast JSON and tables.
The Rust library exports only `run_cli` as the binary bootstrap; statistical,
model, solver and cache internals are not a public Rust API.

```sh
rsomics-limma --threads 4 run --counts counts.tsv --samples samples.tsv --design design.tsv --bundle fit.json \
  '--contrast=(Intercept)=0' --contrast=typesingle-read=0 --contrast=conditiontreated=1
rsomics-limma results --bundle fit.json --output forward.tsv --state forward.json \
  '--contrast=(Intercept)=0' --contrast=typesingle-read=0 --contrast=conditiontreated=1
rsomics-limma export --bundle fit.json --output weights.tsv --kind weights
```

Counts are finite/nonnegative, sample and design rows join by explicit sample
identity, and each contrast names every raw design coefficient including zero.
Numeric design input does not infer factors, formulas, donors or pairing.
Filtering/TMM, current structural-zero voom fitting, default automatic
equal/unequal-df empirical Bayes, per-row precision and mixed-contrast standard
errors belong to this product. It does not depend on another Layer B product.

The bundle preserves raw counts and filter mapping, sample/library/design
identities, full uncontrasted base fit, logCPM/weights/zero-loss diagnostics and
the initial complete contrast/moderation state. Reloaded results transform and
moderate the saved base without filtering, normalization, voom or linear-model
refitting. Optional explicit `--state` saves full contrast/moderation JSON and
TSV as a common atomic pair. Both outputs must differ from the bundle and each
other; no implicit sidecars or large stdout state. Exports are voom logCPM,
observation weights, base-model fitted values or raw residuals. All paths are
named files; `-` is rejected.

NA, NaN, Inf and -Inf numeric states have distinct explicit string tokens;
finite values remain JSON numbers. Scalar versus per-row prior variance is
tagged. Infinite prior df is legitimate in the default equal-df route; a zero
residual df can have an NA sigma. Neither is an unexplained JSON null.

Not exposed: log-expression/missing count input, sample/prior weights, supplied
libraries/offsets, blocks/correlation, other normalizations, forced legacy,
robust/trend/TREAT modes, multi-contrast F, formula engines, plotting UI,
batch/splice/genesets or cross-product fitted objects. No placeholders.

Use `bash scripts/guarded-cargo.sh` for local builds/tests. It verifies external
homes/target/TMP, actual Cargo target, free disks and boot APFS below80% before
the fixed Rust1.91 toolchain. Original source-oracle checks are in
scripts/check-workflow.py; numerical bands were adopted before native output.

Full project license texts are LICENSE-MIT and LICENSE-APACHE-2.0; dependency
license/notice attributions are in THIRD_PARTY_NOTICES.md.
