# rsomics-limma

Private, unpublished implementation of ordinary supplied log-expression and
current count-input voom linear-model lifecycles. Source compatibility, native
platforms, performance and release qualification are separate gates;
`publish=false`.

The supported interface is the CLI, versioned analysis/contrast JSON and tables.
The Rust library exports only `run_cli` as the binary bootstrap; statistical,
model, solver and cache internals are not a public Rust API.

```sh
rsomics-limma --threads 4 run --counts counts.tsv --samples samples.tsv --design design.tsv --bundle fit.json \
  '--contrast=(Intercept)=0' --contrast=typesingle-read=0 --contrast=conditiontreated=1
rsomics-limma results --bundle fit.json --output forward.tsv --state forward.json \
  '--contrast=(Intercept)=0' --contrast=typesingle-read=0 --contrast=conditiontreated=1
rsomics-limma export --bundle fit.json --output weights.tsv --kind weights
rsomics-limma --threads 4 run --expression expression.tsv --samples samples.tsv --design design.tsv --bundle ordinary.json \
  '--contrast=(Intercept)=0' --contrast=typesingle-read=0 --contrast=conditiontreated=1
rsomics-limma export --bundle ordinary.json --output expression-copy.tsv --kind expression
```

Counts are finite/nonnegative, sample and design rows join by explicit sample
identity, and each contrast names every raw design coefficient including zero.
Numeric design input does not infer factors, formulas, donors or pairing.
Filtering/TMM, current structural-zero voom fitting, default automatic
equal/unequal-df empirical Bayes, per-row precision and mixed-contrast standard
errors belong to this product. It does not depend on another Layer B product.

`--counts` and `--expression` are mutually exclusive and one is required.
Expression values are finite, signed and already logged: ordinary fitting does
not filter, transform, normalize or weight them. Both modes require a named
full-column-rank numeric design with more samples than coefficients. Ordinary
fitting reuses one unweighted QR system across all features and retains every
supplied feature ID, including one- and two-feature analyses.

The bundle preserves raw counts and filter mapping, sample/library/design
identities, full uncontrasted base fit, logCPM/weights/zero-loss diagnostics and
the initial complete contrast/moderation state. Reloaded results transform and
moderate the saved base without filtering, normalization, voom or linear-model
refitting. Count bundles retain the schema1 flat base-fit encoding. Ordinary
bundles use explicit schema2/input_kind=log-expression with supplied E,
input digests, metadata/design, the shared full base fit, source-reference
provenance and initial contrast/moderation. They do not contain fabricated
library, voom-weight or structural-zero fields. Source-reference provenance
does not assert the origin or biological validity of user-supplied data.

Optional explicit `--state` saves full contrast/moderation JSON and
TSV as a common atomic pair. Both outputs must differ from the bundle and each
other; no implicit sidecars or large stdout state. Exports are voom logCPM,
observation weights, base-model fitted values or raw residuals. `expression`
exports model expression in either mode; `log-cpm` and `weights` fail explicitly
for an ordinary bundle. All paths are
named files; `-` is rejected.

NA, NaN, Inf and -Inf numeric states have distinct explicit string tokens;
finite values remain JSON numbers. Scalar versus per-row prior variance is
tagged. Infinite prior df is legitimate in the default equal-df route; a zero
residual df can have an NA sigma. Neither is an unexplained JSON null.
Tiny analyses retain the source no-shrink variance model. Zero residual
variance can legitimately produce NaN statistics and an infinite coefficient
prior with recorded fallback; these are not replaced by epsilon variances or
fake significant rows. Native outputs keep supplied IDs; retained source
single-feature name loss is handled only by explicit comparison mappings.
Effects must remain finite and unscaled standard errors finite and positive.
Unsupported arithmetic overflow or degenerate contrast variance fails before
replacing bundle/table/state outputs; it is not encoded as JSON null or treated
as legitimate zero-residual-variance statistics.

Not exposed: missing expression/count input, sample/prior weights, supplied
libraries/offsets, blocks/correlation, other normalizations, forced legacy,
robust/trend/TREAT modes, multi-contrast F, formula engines, plotting UI,
batch/splice/genesets or cross-product fitted objects. No placeholders.

Use `bash scripts/guarded-cargo.sh` for local builds/tests. It verifies external
Cargo home/resolved target, rustup/toolchain, TMPDIR and caches physically,
then invokes the fixed Rust1.91 toolchain. Unrelated boot-disk occupancy does
not block work. Original source-oracle checks are in
scripts/check-workflow.py; numerical bands were adopted before native output.
`scripts/check-logexpression.py` checks all five frozen ordinary source cases
at one/four workers, complete states/tables/exports and unchanged saved bundles
without executing R. Its source identity and numerical bounds are fixed.

Full project license texts are LICENSE-MIT and LICENSE-APACHE-2.0; dependency
license/notice attributions are in THIRD_PARTY_NOTICES.md.
