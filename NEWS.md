# PhyloIgSeq 0.4.0

* The MA-plot **change (`M`) axis is now pluggable**, via `change_transform` on
  `get_ma_coordinates()`, `get_ma_plot_data()`, `get_slide_z()`, `plot_ma()` and
  `plot_slide_z()`. `log2(pos/neg)` is a base-2 logit up to an additive constant, and the
  sliding Z-score centers and scales, so any *affine* reparameterization of that axis
  leaves the score bit-identical; only a non-affine axis changes anything.
* New score `"purity_corrected_slide_z"`: the sliding Z-score computed on
  `change_transform = "purity_corrected"`, i.e. `log2(a_t/b_t)` on the populations
  recovered by un-mixing the two sorted fractions with their measured Ig+ frequencies.
  Unlike `purity_corrected_prob_*` it needs **no pre-sort fraction** — the pre-sort
  frequency enters that axis only as an additive constant, which the centering removes —
  so `getPhyloIgSeq()`'s purity gate is now applied per score rather than as one block.
  Measured on a 50-sample cohort, it agrees closely with `"slide_z"` for taxa inside the
  model's admissible range (median |Δz| 0.21) and diverges in the tails (median |Δz| 2.81),
  which is where the un-mixing is curved.
* **Breaking:** the MA geometry underlying a sliding Z-score moved out of `ig_coating`
  into a new `ma_coords` slot, long in `change_transform` — two change axes cannot share
  one set of `obs_change`/`null_change` columns. `obs_change`, `obs_abundance`,
  `null_change`, `null_abundance` and `ellipse_level` are no longer `ig_coating` columns;
  use the new `ma_coords()` accessor. `ellipse_coords` gains a `change_transform` column.
* New `ig_fraction_names()` returns which `ig_coating` columns hold per-fraction
  abundances, derived from the object's own fraction-name slots. Prefer it to subtracting a
  hand-maintained set of identifier/score/diagnostic names from `colnames(ig_coating)`,
  which silently absorbs any new diagnostic column as though it were a fraction.
* **Behaviour change:** the un-mixed populations are now **regularized** — both are
  clamped at 0 and then offset by `pool_prior`, one read's worth of composition by default.
  Bare clamping scored a taxon the un-mixing places entirely in one population at exactly
  0 or 1, whose logit is infinite, so `purity_corrected_prob_ratio` returned `NA` for
  precisely the most strongly coated taxa. Such taxa now saturate: finite, correctly
  ordered, and shrunk toward the null by `O(pool_prior/a_t)`, so a rare taxon that escaped
  on sampling noise is pulled back hard while an abundant, genuinely coated one barely
  moves. The offset is symmetric, so the model's fixed point is untouched. This changes
  `purity_corrected_prob_index`/`_ratio` values; pass `pool_prior = 0` for the old
  behaviour. Falling outside the admissible range is common rather than exceptional — a
  median of 41.5% of taxa per sample on the cohort above.
* **Behaviour change:** `compute_slide_z()` now windows the imputed-taxa population among
  themselves instead of scoring it against one global mean/sd. Imputed taxa span an
  abundance gradient of their own, and pooling discarded it; with fewer of them than
  `window_size` the arithmetic collapses to a single window, i.e. the previous behaviour.
  They remain a population apart — measured on `ps_igseq` under `pseudo_count` their null
  is 2.8x-6.9x wider than that of taxa with no imputed zero, and in some samples the whole
  bottom 60% of the abundance range is imputed, so there is no unimputed neighbour to
  estimate a local null from.
* `compute_slide_z()` gains `estimate_from`, marking taxa that are scored but may not
  inform a window's mean/sd. Out-of-cone taxa on the purity-corrected axis use it: their
  reference value is saturated rather than measured, so it must not enter the estimate,
  but they still get a score on their window's scale. The mask follows whichever pair
  supplies the reference — `null_in_cone` in empirical-null mode, `obs_in_cone` otherwise.
* A window with too few usable taxa to estimate its own null now falls back to the estimate
  pooled over that population, with a `warning()`, instead of emitting `NA` or a scale
  derived from one or two points.
* The two change axes differ in what they read: `"purity_corrected"` is computed on
  relative abundances, because the un-mixing needs compositions outright, whereas
  `"log_ratio"` stays on the abundances exactly as supplied. Normalizing each member of a
  pair by its own total adds a per-pair constant, and the sliding Z-score cancels an
  additive constant only when the observed and the null pair carry the *same* one — which
  they do not, since the two pairs have different totals even under per-sample
  rarefaction. Measured on `ps_igseq` with `zero_treatment = "no_zero"`, normalizing
  `"log_ratio"` shifted whole samples by up to 1.5 Z units and flipped 3.9% of
  `|Z| > 1.96` calls. `"purity_corrected"` cannot avoid this and carries a residual
  per-sample offset `"log_ratio"` does not — a known limitation of that axis.
* `resolve_ig_freqs()` is now exported, so callers that bypass `getPhyloIgSeq()` can obtain
  validated per-fraction Ig+ frequencies without re-implementing the unit conversion and
  range checks.
* New `SLIDE_Z_SCORES`, mapping each windowed score to the change axis it is computed on.
* `getPhyloIgSeq()` now drops `prob_index`/`prob_ratio` with a `warning()` when the pre-sort
  P(Ig+) they weight by is unreadable, instead of returning silent all-NA columns — the same
  treatment the `purity_corrected_*` scores already got. Gated per score, since their inputs
  differ: `prob_index` divides by the pre-sort fraction's own abundance and so needs that
  fraction under either layout, whereas `prob_ratio` divides by the negative fraction and
  needs the pre-sort fraction only under `ig_freq_layout = "long"`, where P(Ig+) is recorded
  per fraction and is therefore readable only off the pre-sort rows.
* `plot_ma()` and `plot_slide_z()` draw out-of-cone taxa **in place but hollow**, rather
  than in the off-axis jitter band: their abundance is real and only their change value is
  saturated. The saturation ceiling grows with abundance, so they trace a rising envelope
  rather than a flat clipped line.
* **Bug fix:** `plot_ma()` decided which taxa to send to the off-axis jitter band from the
  `imputed_taxa` union, so every `zero_treatment` facet banished the same set of taxa —
  including under treatments that impute nothing. `get_ma_plot_data()`'s `plot_data` now
  carries a per-row `imputed` flag, i.e. whether *that* row's treatment derived its
  coordinates from an imputed zero, and `plot_ma()` reads that instead. `imputed_taxa` is
  unchanged (still the union over treatments) for callers that want one set per sample.
* `plot_ma()` and `plot_slide_z()` no longer label the negative part of the log-abundance
  axis. That region holds the off-axis band imputed taxa are parked in, which is not an
  abundance at all, plus taxa carrying less than one read's worth of signal. Only the
  breaks are dropped -- limits are untouched, so nothing is clipped -- and the negatives
  are kept when the axis is legitimately all-negative, as under
  `transform_by_sample = "compositional"`.
* `plot_ma(type = "superposed")` no longer lets the last `zero_treatment` hide all the
  others: the ordinal size scale ran smallest-to-largest in draw order, so the treatment
  drawn last was also the biggest. It now runs largest-to-smallest.

# PhyloIgSeq 0.3.0

* `getPhyloIgSeq()` gains `ig_freq_layout`, controlling what one value of the
  `ig_freq_name` column means. `"wide"` (the default, and the previous behaviour) is one
  Ig+ frequency per biological sample, repeated across its fraction rows. `"long"` reads
  the column as *fraction-specific* — each fraction's row records the Ig+ frequency
  measured in that fraction — yielding the pre-sort P(Ig+) plus the Ig+ fraction's purity
  and the Ig− fraction's impurity.
* The resolved (unit-converted, range-checked) frequencies are now stored per sample in
  the returned object's `sample_data` slot, as `presort_ig_freq`, `pos_ig_freq`,
  `neg1_ig_freq` and `neg2_ig_freq` (the last three are `NA` under `"wide"`). With several
  `positive_fraction_name`s, `pos_ig_freq` follows each row's own positive fraction. New
  `ig_freq_layout` slot records which layout was used.
* Under `"wide"`, a column that varies across a sample's fractions now produces a
  `warning()` pointing at `"long"`, instead of silently collapsing to `NA` and returning
  all-`NA` scores.
* `"purity_corrected_prob_index"`/`"purity_corrected_prob_ratio"` joined `IG_SCORES` and
  are computed by `getPhyloIgSeq()` when their inputs are available (a negative fraction,
  a pre-sort fraction and `ig_freq_layout = "long"`); they are dropped from `scores` with
  a `warning()` otherwise. They remain experimental — the sorting parameter they depend on
  is derived from the three measured frequencies, with the alternatives and caveats
  documented in `compute_ig_score()`'s source.
* **Breaking:** `compute_ig_score()`'s probability arguments were renamed for one
  consistent vocabulary — `ig_freq` → `presort_ig_freq`, `pos_purity` → `pos_ig_freq`,
  `neg_impurity` → `neg_ig_freq` — and `pos_fraction`/`neg_fraction` were removed (the
  quantity they described is now derived internally).

# PhyloIgSeq 0.2.0

* `getPhyloIgSeq()`'s `positive_fraction_name` now accepts a character vector of
  several positive fraction values sharing the same negative/pre-sort fraction(s)
  (e.g. multiple phage-display sort rounds/output pools). Each (sample, positive
  fraction) combination becomes its own row, keyed by a synthetic
  `"<sample_id>_<positive fraction>"` sample_id, with two new `sample_data` columns
  (`original_sample_id`, `positive_fraction_name`) and a coalesced
  `positive_fraction_abundance` column. The negative/pre-sort fraction(s) are
  rarefied once, jointly across all positive fractions of a sample, so comparisons
  stay at a common depth. Behavior with a single positive fraction (the default) is
  unchanged.
* New toy dataset `ps_phage_display`: a synthetic phage-display panning experiment
  (8 campaigns, one shared "input" library plus three sequential positive output
  rounds) demonstrating the multiple-positive-fractions workflow above.

# PhyloIgSeq 0.1.0

* Initial CRAN submission.
