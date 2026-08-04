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
