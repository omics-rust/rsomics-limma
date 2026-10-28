# Implementation provenance

Current behavior references are edgeR4.10.7, commit
27ba7b39d819e68df2ecdea9f32a4b6039b2a0b5, archive
a3cc0f66be29f1a61a74021f9848ec940794e653905c4332327a719dd9ad8393,
and limma3.68.5, commit825d1c8330d024adb6c54af66d79b0f4ec258530,
archive40b1d2059a76ad98db4bf80139760acf8e8bae92a8265aaaf16ac199831533f4.
Actual source was read; new GPL R/C implementation text or fitted tables were
not transplanted. This is engineering provenance, not legal clearance.

Selected team-owned assets were merged/refactored privately from edgeR main
3a97e8f4ea510b6d83d9e3554b71e9eab0b37add: named input/hash/contrast joins,
filter/TMM libraries, Poisson/general scalar solver arithmetic, ordinary LOWESS
and mathematical Brent/g/trigamma. Input accepts upstream fractional counts;
no unused held-AveLogCPM, QL moments, QL coefficients or abundance-EB policy is
carried into limma. Its prior QL comparison failure is not a limma acceptance
or an excuse to change this product's frozen bands.

Historical seed provenance inherited through those assets:

| Asset | Revision | Use |
|---|---|---|
| filter-by-expr | 451580a86a4c36de536605fc08037f2c4b1b83c7 | design-filter/current-source arithmetic |
| tmm-norm | 3f37f84a8935c419a37378bc43dc09930da83067 | rank/reference/trim arithmetic |
| edger-glm-qlf | c2f26cf9d9dc755961a8f7ea529ff28ca65262d8 | NB/Poisson mathematical seed through current reviewed solver, not legacy QL |
| voom-quality-weights | 18ad63de6ed12827da6edc2625fd2f99cd7c36cc | weighted-QR/voom fixture/math reference; old fit state and raw-library shortcuts discarded |
| limma-ebayes | 8fbfa0b8fcb9444dcf5740a4c4feab7a398e9c30 | contrast/mixture fixture/math reference; constant df/SE and positional parsers discarded |
| ebayes-core | historical team-owned local asset | equal-df mathematical/fixture reference; no dependency or revived public API, bounded inverse solve |

Weighted voom LOWESS is an independently written four-pass extension of the
team's local-linear frequency-weight arithmetic. Original df-weight windows,
weighted-MAD robustness, ties, span/delta and interpolation were traced in
limma's actual R/C source; no upstream C implementation text copied. QR/WLS,
scaled-F likelihood/moments, Student-t/beta identities and posterior-odds
formulas are implemented independently with typed product-local state.
statrs0.18 is MIT; existing common/help/stats are concrete consumers. No new
public numerical API or B→B/path dependency was added.

The unchanged compatibility authority is accepted R4.6.1/Bioc3.23 source run
40606516204: original14599x7 integer pasilla counts are raw data carried by
DESeq2, not DESeq fitting. Full source references and formula plan are in
world/docs/plans/2026-10-28-limma-native-count-workflow.md. Equal-df/Inf and
ordinary-no-zero routes are real implemented defaults, but are not claimed
source-covered by this witness. Broad numeric-design/fractional-input and
release qualification require their own proportionate evidence.
