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
