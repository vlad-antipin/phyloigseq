library(PhyloIgSeq)

# ---- getPhyloIgSeq: positive-fraction-only path ----

# ps_igseq ships a "whole" (pre-sort) fraction and an "ig_pheno" phenotype
# column (see ?ps_igseq) specifically so the positive-only path (positive +
# presort + ig_freq, no negative fraction) can be exercised end-to-end.
make_pos_only_physeq <- function() {
  data("ps_igseq", package = "PhyloIgSeq", envir = environment())
  physeq <- ps_igseq

  sd <- as(phyloseq::sample_data(physeq), "data.frame")
  # sample_6/sample_7 carry intentionally duplicated pos/neg1/neg2 rows (used
  # to exercise group_sorted_samples()'s duplicate-exclusion elsewhere); skip
  # them here so this fixture stays a clean 2-fraction-per-sample case.
  keep <- sd$sorting_fraction %in% c("pos", "whole") &
    !sd$sample_id %in% c("sample_6", "sample_7")
  phyloseq::prune_samples(rownames(sd)[keep], physeq)
}

test_that("getPhyloIgSeq computes prob_index from positive + presort + ig_freq, no negative fraction", {
  physeq <- make_pos_only_physeq()

  pis <- getPhyloIgSeq(
    physeq = physeq,
    sample_id_name = "sample_id",
    fraction_id_name = "sorting_fraction",
    positive_fraction_name = "pos",
    first_negative_fraction_name = NULL,
    presorting_fraction_name = "whole",
    ig_freq_name = "ig_pheno",
    scores = "prob_index",
    empirical_null_distribution = FALSE
  )

  expect_s4_class(pis, "PhyloIgSeq")
  expect_null(pis@first_negative_fraction_name)
  expect_equal(pis@score_names, "prob_index")
  expect_true("prob_index" %in% names(pis@ig_coating))
  expect_false(all(is.na(pis@ig_coating$prob_index)))
})

test_that("getPhyloIgSeq errors when no negative fraction and presort/ig_freq are incomplete", {
  physeq <- make_pos_only_physeq()

  expect_error(
    getPhyloIgSeq(
      physeq = physeq,
      sample_id_name = "sample_id",
      fraction_id_name = "sorting_fraction",
      positive_fraction_name = "pos",
      first_negative_fraction_name = NULL,
      presorting_fraction_name = "whole",
      ig_freq_name = NULL,
      scores = "prob_index"
    ),
    "first_negative_fraction_name"
  )

  expect_error(
    getPhyloIgSeq(
      physeq = physeq,
      sample_id_name = "sample_id",
      fraction_id_name = "sorting_fraction",
      positive_fraction_name = "pos",
      first_negative_fraction_name = NULL,
      presorting_fraction_name = NULL,
      ig_freq_name = "ig_pheno",
      scores = "prob_index"
    ),
    "first_negative_fraction_name"
  )
})

test_that("getPhyloIgSeq drops slide_z with a warning when no negative fraction is furnished", {
  physeq <- make_pos_only_physeq()

  expect_warning(
    pis <- getPhyloIgSeq(
      physeq = physeq,
      sample_id_name = "sample_id",
      fraction_id_name = "sorting_fraction",
      positive_fraction_name = "pos",
      first_negative_fraction_name = NULL,
      presorting_fraction_name = "whole",
      ig_freq_name = "ig_pheno",
      scores = c("slide_z", "prob_index"),
      empirical_null_distribution = FALSE
    ),
    "slide_z"
  )

  expect_false("slide_z" %in% pis@score_names)
  expect_true("prob_index" %in% pis@score_names)
})

# ---- getPhyloIgSeq: ig_freq_units ----

test_that("getPhyloIgSeq's ig_freq_units = 'percent' matches an equivalent 'frequency' value", {
  physeq_freq <- make_pos_only_physeq()

  physeq_pct <- physeq_freq
  sd_pct <- as(phyloseq::sample_data(physeq_pct), "data.frame")
  sd_pct$ig_pheno <- sd_pct$ig_pheno * 100
  phyloseq::sample_data(physeq_pct) <- phyloseq::sample_data(sd_pct)

  # rarefy_by_sample = FALSE -- rarefaction is a random draw, so two
  # independent getPhyloIgSeq() calls with it enabled would differ in pos/pre
  # counts alone, confounding this comparison of ig_freq_units specifically.
  pis_freq <- getPhyloIgSeq(
    physeq = physeq_freq,
    sample_id_name = "sample_id",
    fraction_id_name = "sorting_fraction",
    positive_fraction_name = "pos",
    presorting_fraction_name = "whole",
    ig_freq_name = "ig_pheno",
    ig_freq_units = "frequency",
    rarefy_by_sample = FALSE,
    scores = "prob_index",
    empirical_null_distribution = FALSE
  )
  pis_pct <- getPhyloIgSeq(
    physeq = physeq_pct,
    sample_id_name = "sample_id",
    fraction_id_name = "sorting_fraction",
    positive_fraction_name = "pos",
    presorting_fraction_name = "whole",
    ig_freq_name = "ig_pheno",
    ig_freq_units = "percent",
    rarefy_by_sample = FALSE,
    scores = "prob_index",
    empirical_null_distribution = FALSE
  )

  expect_equal(
    pis_pct@ig_coating[order(pis_pct@ig_coating$sample_id, pis_pct@ig_coating$taxon_id), "prob_index"],
    pis_freq@ig_coating[order(pis_freq@ig_coating$sample_id, pis_freq@ig_coating$taxon_id), "prob_index"]
  )
})

test_that("getPhyloIgSeq treats an out-of-range ig_freq (after unit conversion) as NA with a warning, not an error", {
  physeq <- make_pos_only_physeq()
  sd <- as(phyloseq::sample_data(physeq), "data.frame")
  bad_sample <- sd$sample_id[1]
  sd$ig_pheno[sd$sample_id == bad_sample] <- 1.5 # not a valid [0, 1] frequency
  phyloseq::sample_data(physeq) <- phyloseq::sample_data(sd)

  expect_warning(
    pis <- getPhyloIgSeq(
      physeq = physeq,
      sample_id_name = "sample_id",
      fraction_id_name = "sorting_fraction",
      positive_fraction_name = "pos",
      presorting_fraction_name = "whole",
      ig_freq_name = "ig_pheno",
      ig_freq_units = "frequency",
      scores = "prob_index",
      empirical_null_distribution = FALSE
    ),
    "outside the expected \\[0, 1\\]"
  )

  expect_true(all(is.na(
    pis@ig_coating$prob_index[pis@ig_coating$sample_id == bad_sample]
  )))
  expect_false(all(is.na(
    pis@ig_coating$prob_index[pis@ig_coating$sample_id != bad_sample]
  )))
})

test_that("getPhyloIgSeq treats a non-numeric ig_freq column as NA with a warning, not an error", {
  physeq <- make_pos_only_physeq()
  sd <- as(phyloseq::sample_data(physeq), "data.frame")
  sd$ig_pheno <- sd$sex # a character column, not a numeric probability
  phyloseq::sample_data(physeq) <- phyloseq::sample_data(sd)
  # A single sample is enough to exercise the check; keeps this test to one
  # warning instead of one per sample.
  physeq <- phyloseq::prune_samples(
    rownames(sd)[sd$sample_id == sd$sample_id[1]],
    physeq
  )

  expect_warning(
    pis <- getPhyloIgSeq(
      physeq = physeq,
      sample_id_name = "sample_id",
      fraction_id_name = "sorting_fraction",
      positive_fraction_name = "pos",
      presorting_fraction_name = "whole",
      ig_freq_name = "ig_pheno",
      ig_freq_units = "percent", # would divide by 100 if not caught first
      scores = "prob_index",
      empirical_null_distribution = FALSE
    ),
    "not numeric"
  )

  expect_true(all(is.na(pis@ig_coating$prob_index)))
})
