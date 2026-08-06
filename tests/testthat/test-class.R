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

# ---- getPhyloIgSeq: multiple positive fractions ----

# Full pos/neg1/neg2/whole fixture (unlike make_pos_only_physeq(), keeps
# neg1 too), still skipping sample_6/sample_7's intentionally duplicated
# fraction rows.
make_full_physeq <- function() {
  data("ps_igseq", package = "PhyloIgSeq", envir = environment())
  physeq <- ps_igseq
  sd <- as(phyloseq::sample_data(physeq), "data.frame")
  keep <- !sd$sample_id %in% c("sample_6", "sample_7")
  phyloseq::prune_samples(rownames(sd)[keep], physeq)
}

test_that("getPhyloIgSeq with a single positive fraction is unchanged (no suffix, no new columns)", {
  physeq <- make_full_physeq()

  pis <- getPhyloIgSeq(
    physeq = physeq,
    sample_id_name = "sample_id",
    fraction_id_name = "sorting_fraction",
    positive_fraction_name = "pos",
    first_negative_fraction_name = "neg1",
    scores = "palm"
  )

  expect_false("original_sample_id" %in% names(pis@sample_data))
  expect_false("positive_fraction_name" %in% names(pis@sample_data))
  expect_true("pos" %in% names(pis@ig_coating))
  expect_false("positive_fraction_abundance" %in% names(pis@ig_coating))
  expect_true(all(grepl("^sample_[0-9]+$", unique(pis@ig_coating$sample_id))))
})

test_that("getPhyloIgSeq with multiple positive fractions folds them into the sample dimension", {
  physeq <- make_full_physeq()

  # "neg2" stands in for a second positive fraction here purely to exercise
  # the multi-fraction mechanics on this toy dataset (see the vignette for
  # the same choice, with rationale).
  pis <- getPhyloIgSeq(
    physeq = physeq,
    sample_id_name = "sample_id",
    fraction_id_name = "sorting_fraction",
    positive_fraction_name = c("pos", "neg2"),
    first_negative_fraction_name = "neg1",
    scores = "palm",
    rarefy_by_sample = FALSE
  )

  expect_s4_class(pis, "PhyloIgSeq")

  # One synthetic sample per (original sample, positive fraction) pair
  expect_true(all(grepl("_(pos|neg2)$", unique(pis@ig_coating$sample_id))))
  expect_setequal(unique(pis@sample_data$positive_fraction_name), c("pos", "neg2"))
  expect_setequal(
    unique(pis@sample_data$original_sample_id),
    paste0("sample_", c(1:5, 8))
  )

  # Coalesced abundance column always populated, regardless of which
  # positive fraction the row came from
  expect_true("positive_fraction_abundance" %in% names(pis@ig_coating))
  expect_false("pos" %in% names(pis@ig_coating))
  expect_false("neg2" %in% names(pis@ig_coating))
  expect_false(any(is.na(pis@ig_coating$positive_fraction_abundance)))

  # Shared negative fraction stays fully populated across every row (not
  # sparse like the coalesced positive column would have been without the fix)
  expect_true("neg1" %in% names(pis@ig_coating))
  expect_false(any(is.na(pis@ig_coating$neg1)))
})

test_that("getPhyloIgSeq rarefies the shared negative fraction once across all positive fractions of a sample", {
  physeq <- make_full_physeq()

  pis <- getPhyloIgSeq(
    physeq = physeq,
    sample_id_name = "sample_id",
    fraction_id_name = "sorting_fraction",
    positive_fraction_name = c("pos", "neg2"),
    first_negative_fraction_name = "neg1",
    scores = "palm",
    rarefy_by_sample = TRUE
  )

  ig <- pis@ig_coating
  ig$original_sample_id <- sub("_(pos|neg2)$", "", ig$sample_id)
  ig$fraction <- sub("^.*_", "", ig$sample_id)

  # impute_zeros() runs separately per positive fraction and can legitimately
  # drop a different set of taxa each time (by design -- a zero in one
  # comparison shouldn't drop a taxon from the other), so compare shared
  # neg1 *values* for taxa common to both fractions of a sample, rather than
  # each fraction's neg1 total (which would differ whenever the surviving
  # taxon sets differ, even with a correctly shared rarefaction).
  for (orig_id in unique(ig$original_sample_id)) {
    sub_ig <- ig[ig$original_sample_id == orig_id, ]
    pos_df <- sub_ig[sub_ig$fraction == "pos", c("taxon_id", "neg1")]
    neg2_df <- sub_ig[sub_ig$fraction == "neg2", c("taxon_id", "neg1")]
    merged <- merge(pos_df, neg2_df, by = "taxon_id", suffixes = c("_pos", "_neg2"))
    expect_gt(nrow(merged), 0)
    expect_equal(merged$neg1_pos, merged$neg1_neg2)
  }
})

test_that("show() and plot_slide_z() handle a multi-positive-fraction object without erroring", {
  physeq <- make_full_physeq()

  # suppressWarnings(): this fixture has no second negative fraction or
  # precomputed MA/ellipse coordinates, so getPhyloIgSeq()/plot_slide_z() take
  # their documented fallback paths (observed pos-neg distribution as null, no
  # ellipses) with an informational warning each; plot_slide_z()'s pre-existing
  # geom_jitter()/geom_point() `text` aesthetic (for the app's plotly tooltip)
  # and discrete `size` scale also warn on every call. All expected/harmless
  # here -- muffled so they don't drown out this test's own signal.
  pis <- suppressWarnings(getPhyloIgSeq(
    physeq = physeq,
    sample_id_name = "sample_id",
    fraction_id_name = "sorting_fraction",
    positive_fraction_name = c("pos", "neg2"),
    first_negative_fraction_name = "neg1",
    scores = "slide_z",
    rarefy_by_sample = FALSE
  ))

  expect_output(print(pis), "pos, neg2")
  expect_s3_class(
    suppressWarnings(plot_slide_z(pis, sample_ids = unique(pis@ig_coating$sample_id)[1])),
    "ggplot"
  )
})

# ---- getPhyloIgSeq: ig_freq_layout ----

# ps_igseq's own `ig_pheno` is constant across a sample's fractions (the "wide"
# layout). Real IgSeq metadata often records the Ig+ frequency measured *in each
# fraction* instead, so build that layout here rather than reshaping the shipped
# dataset, which the wide-layout tests above depend on.
make_long_ig_freq_physeq <- function(
  presort = 0.40,
  pos = 0.80,
  neg = 0.10,
  samples = NULL
) {
  data("ps_igseq", package = "PhyloIgSeq", envir = environment())
  physeq <- ps_igseq

  sd <- as(phyloseq::sample_data(physeq), "data.frame")
  # sample_6/sample_7 carry duplicated fraction rows on purpose (see
  # make_pos_only_physeq()); a duplicated fraction has no single Ig+ frequency.
  keep <- !sd$sample_id %in% c("sample_6", "sample_7")
  if (!is.null(samples)) {
    keep <- keep & sd$sample_id %in% samples
  }
  physeq <- phyloseq::prune_samples(rownames(sd)[keep], physeq)

  sd <- as(phyloseq::sample_data(physeq), "data.frame")
  sd$ig_pheno_fraction <- c(
    whole = presort,
    pos = pos,
    neg1 = neg,
    # neg1 and neg2 are two aliquots of the same negative fraction, so they
    # carry the same measurement
    neg2 = neg
  )[sd$sorting_fraction]
  phyloseq::sample_data(physeq) <- phyloseq::sample_data(sd)
  physeq
}

run_long <- function(physeq, ...) {
  getPhyloIgSeq(
    physeq = physeq,
    sample_id_name = "sample_id",
    fraction_id_name = "sorting_fraction",
    positive_fraction_name = "pos",
    first_negative_fraction_name = "neg1",
    second_negative_fraction_name = "neg2",
    presorting_fraction_name = "whole",
    ig_freq_name = "ig_pheno_fraction",
    ig_freq_layout = "long",
    empirical_null_distribution = FALSE,
    ...
  )
}

test_that("ig_freq_layout = 'long' reads each fraction's own Ig+ frequency into sample_data", {
  pis <- run_long(make_long_ig_freq_physeq(), scores = "prob_index")

  expect_equal(pis@ig_freq_layout, "long")
  expect_equal(pis@ig_freq_name, "ig_pheno_fraction")

  sd <- pis@sample_data
  expect_true(all(
    c("presort_ig_freq", "pos_ig_freq", "neg1_ig_freq", "neg2_ig_freq") %in%
      names(sd)
  ))
  expect_true(all(sd$presort_ig_freq == 0.40))
  expect_true(all(sd$pos_ig_freq == 0.80))
  expect_true(all(sd$neg1_ig_freq == 0.10))
  # neg1 and neg2 are aliquots of the same fraction
  expect_equal(sd$neg1_ig_freq, sd$neg2_ig_freq)

  # prob_index must use the *pre-sort* value, not any of the in-fraction ones
  expect_false(all(is.na(pis@ig_coating$prob_index)))
})

test_that("ig_freq_layout = 'wide' still collapses to a single per-sample frequency", {
  data("ps_igseq", package = "PhyloIgSeq", envir = environment())
  physeq <- ps_igseq
  sd <- as(phyloseq::sample_data(physeq), "data.frame")
  physeq <- phyloseq::prune_samples(
    rownames(sd)[!sd$sample_id %in% c("sample_6", "sample_7")],
    physeq
  )

  pis <- getPhyloIgSeq(
    physeq = physeq,
    sample_id_name = "sample_id",
    fraction_id_name = "sorting_fraction",
    positive_fraction_name = "pos",
    first_negative_fraction_name = "neg1",
    presorting_fraction_name = "whole",
    ig_freq_name = "ig_pheno",
    scores = "prob_index",
    empirical_null_distribution = FALSE
  )

  expect_equal(pis@ig_freq_layout, "wide")
  sd_out <- pis@sample_data
  # the sample's single ig_pheno lands in presort_ig_freq; nothing is known
  # about how pure either sorted fraction turned out to be
  expect_false(any(is.na(sd_out$presort_ig_freq)))
  expect_true(all(is.na(sd_out$pos_ig_freq)))
  expect_true(all(is.na(sd_out$neg1_ig_freq)))
})

test_that("ig_freq_layout = 'wide' warns when the column varies across fractions", {
  physeq <- make_long_ig_freq_physeq(samples = "sample_1")

  expect_warning(
    pis <- getPhyloIgSeq(
      physeq = physeq,
      sample_id_name = "sample_id",
      fraction_id_name = "sorting_fraction",
      positive_fraction_name = "pos",
      first_negative_fraction_name = "neg1",
      presorting_fraction_name = "whole",
      ig_freq_name = "ig_pheno_fraction",
      ig_freq_layout = "wide",
      scores = "prob_index",
      empirical_null_distribution = FALSE
    ),
    'ig_freq_layout = "long"'
  )
  expect_true(all(is.na(pis@sample_data$presort_ig_freq)))
})

test_that("ig_freq_units = 'percent' converts every fraction's frequency", {
  physeq <- make_long_ig_freq_physeq(presort = 40, pos = 80, neg = 10)

  pis <- run_long(physeq, ig_freq_units = "percent", scores = "prob_index")

  sd <- pis@sample_data
  expect_true(all(sd$presort_ig_freq == 0.40))
  expect_true(all(sd$pos_ig_freq == 0.80))
  expect_true(all(sd$neg1_ig_freq == 0.10))
})

test_that("getPhyloIgSeq computes the purity-corrected scores under the long layout", {
  pis <- run_long(
    make_long_ig_freq_physeq(),
    scores = c("prob_index", "purity_corrected_prob_index")
  )

  expect_true("purity_corrected_prob_index" %in% names(pis@ig_coating))
  expect_false(all(is.na(pis@ig_coating$purity_corrected_prob_index)))
  # w < 1 here (presort 0.40 < purity 0.80), so the corrected score shrinks
  # towards the Ig- fraction and must differ from plain prob_index
  expect_false(isTRUE(all.equal(
    pis@ig_coating$purity_corrected_prob_index,
    pis@ig_coating$prob_index
  )))
})

test_that("getPhyloIgSeq drops the purity-corrected scores when their inputs are missing", {
  data("ps_igseq", package = "PhyloIgSeq", envir = environment())
  physeq <- ps_igseq
  sd <- as(phyloseq::sample_data(physeq), "data.frame")
  physeq <- phyloseq::prune_samples(
    rownames(sd)[!sd$sample_id %in% c("sample_6", "sample_7")],
    physeq
  )

  # wide layout carries no in-fraction frequencies
  expect_warning(
    pis <- getPhyloIgSeq(
      physeq = physeq,
      sample_id_name = "sample_id",
      fraction_id_name = "sorting_fraction",
      positive_fraction_name = "pos",
      first_negative_fraction_name = "neg1",
      presorting_fraction_name = "whole",
      ig_freq_name = "ig_pheno",
      scores = c("prob_index", "purity_corrected_prob_index"),
      empirical_null_distribution = FALSE
    ),
    "dropping from `scores`"
  )
  expect_equal(pis@score_names, "prob_index")
  expect_false("purity_corrected_prob_index" %in% names(pis@ig_coating))
})

test_that("getPhyloIgSeq drops the Bayesian scores when their prior is unreadable", {
  # `prob_index` divides by the pre-sort fraction's abundance, so it needs that
  # fraction under either layout; `prob_ratio` divides by the negative fraction
  # and only needs the pre-sort fraction under "long", where P(Ig+) is recorded
  # per fraction. Either way, dropped rather than returned as an all-NA column.
  physeq <- make_long_ig_freq_physeq()

  # one warning per score, so both have to be matched
  expect_warning(
    expect_warning(
      pis <- getPhyloIgSeq(
        physeq = physeq,
        sample_id_name = "sample_id",
        fraction_id_name = "sorting_fraction",
        positive_fraction_name = "pos",
        first_negative_fraction_name = "neg1",
        ig_freq_name = "ig_pheno_fraction",
        ig_freq_layout = "long",
        scores = c("prob_index", "prob_ratio", "palm"),
        empirical_null_distribution = FALSE
      ),
      "Cannot compute `prob_index`"
    ),
    "Cannot compute `prob_ratio`"
  )
  expect_equal(pis@score_names, "palm")
  expect_false(any(c("prob_index", "prob_ratio") %in% names(pis@ig_coating)))

  # Under "wide", P(Ig+) is a property of the sample, so `prob_ratio` stays
  # computable without a pre-sort fraction -- only `prob_index` goes.
  expect_warning(
    pis_wide <- getPhyloIgSeq(
      physeq = physeq,
      sample_id_name = "sample_id",
      fraction_id_name = "sorting_fraction",
      positive_fraction_name = "pos",
      first_negative_fraction_name = "neg1",
      ig_freq_name = "ig_pheno",
      ig_freq_layout = "wide",
      scores = c("prob_index", "prob_ratio"),
      empirical_null_distribution = FALSE
    ),
    "Cannot compute `prob_index`"
  )
  expect_equal(pis_wide@score_names, "prob_ratio")
  expect_false(all(is.na(pis_wide@ig_coating$prob_ratio)))
})

test_that("pos_ig_freq follows each row's own positive fraction with several of them", {
  physeq <- make_long_ig_freq_physeq()
  sd <- as(phyloseq::sample_data(physeq), "data.frame")
  # give neg2 a distinguishable frequency, then score it as a second positive
  # fraction so the two synthetic samples must disagree on pos_ig_freq alone
  sd$ig_pheno_fraction[sd$sorting_fraction == "neg2"] <- 0.25
  phyloseq::sample_data(physeq) <- phyloseq::sample_data(sd)

  pis <- getPhyloIgSeq(
    physeq = physeq,
    sample_id_name = "sample_id",
    fraction_id_name = "sorting_fraction",
    positive_fraction_name = c("pos", "neg2"),
    first_negative_fraction_name = "neg1",
    presorting_fraction_name = "whole",
    ig_freq_name = "ig_pheno_fraction",
    ig_freq_layout = "long",
    scores = "prob_index",
    empirical_null_distribution = FALSE
  )

  sd_out <- pis@sample_data
  expect_equal(
    unique(sd_out$pos_ig_freq[sd_out$positive_fraction_name == "pos"]),
    0.80
  )
  expect_equal(
    unique(sd_out$pos_ig_freq[sd_out$positive_fraction_name == "neg2"]),
    0.25
  )
  # the fraction-invariant ones stay the same across both
  expect_equal(unique(sd_out$presort_ig_freq), 0.40)
  expect_equal(unique(sd_out$neg1_ig_freq), 0.10)
})

test_that("getPhyloIgSeq treats a non-existent ig_freq_name column as NA with a warning", {
  physeq <- make_long_ig_freq_physeq(samples = "sample_1")

  expect_warning(
    pis <- getPhyloIgSeq(
      physeq = physeq,
      sample_id_name = "sample_id",
      fraction_id_name = "sorting_fraction",
      positive_fraction_name = "pos",
      first_negative_fraction_name = "neg1",
      presorting_fraction_name = "whole",
      ig_freq_name = "no_such_column",
      ig_freq_layout = "long",
      scores = "prob_index",
      empirical_null_distribution = FALSE
    ),
    "is not a column of sample_data"
  )
  expect_true(all(is.na(pis@ig_coating$prob_index)))
})

# ---- getPhyloIgSeq: two change axes, ma_coords slot, accessors ----

test_that("getPhyloIgSeq stores one ma_coords block per requested change axis", {
  # This fixture's made-up p/q leave a fifth of the taxa outside the admissible
  # cone, which the run reports once; incidental to what is asserted here, and
  # checked on its own below.
  suppressWarnings(
    pis <- run_long(
      make_long_ig_freq_physeq(samples = c("sample_1", "sample_2")),
      scores = c("slide_z", "purity_corrected_slide_z"),
      confidence_levels = c(0.95, 0.99)
    )
  )

  # Both scores are plain columns of ig_coating, so score_names/get_ig_score()
  # and everything keyed off them keep working unchanged.
  expect_true(all(
    c("slide_z", "purity_corrected_slide_z") %in% pis@score_names
  ))
  expect_true(all(
    c("slide_z", "purity_corrected_slide_z") %in% names(pis@ig_coating)
  ))

  coords <- ma_coords(pis)
  expect_setequal(
    unique(coords$change_transform),
    c("log_ratio", "purity_corrected")
  )
  # One row per (sample, taxon) per axis, keyed the same way as ig_coating.
  expect_equal(nrow(coords), 2 * nrow(pis@ig_coating))
  expect_true(all(c("obs_in_cone", "null_in_cone") %in% names(coords)))

  # The cone only exists on the purity-corrected axis.
  log_ratio_block <- ma_coords(pis, "log_ratio")
  purity_block <- ma_coords(pis, "purity_corrected")
  expect_true(all(is.na(log_ratio_block$obs_in_cone)))
  expect_false(all(is.na(purity_block$obs_in_cone)))

  # Ellipses are labelled by axis too, so a consumer can draw the right ones.
  expect_true("change_transform" %in% names(pis@ellipse_coords))
  expect_setequal(
    unique(pis@ellipse_coords$change_transform),
    c("log_ratio", "purity_corrected")
  )
})

test_that("getPhyloIgSeq no longer puts MA geometry in ig_coating", {
  pis <- run_long(
    make_long_ig_freq_physeq(samples = "sample_1"),
    scores = "slide_z"
  )
  expect_false(any(
    c(
      "obs_change",
      "obs_abundance",
      "null_change",
      "null_abundance",
      "ellipse_level"
    ) %in%
      names(pis@ig_coating)
  ))
})

test_that("purity_corrected_slide_z needs no presorting fraction, unlike the prob scores", {
  physeq <- make_long_ig_freq_physeq(samples = c("sample_1", "sample_2"))
  requested <- c(
    "purity_corrected_slide_z",
    "purity_corrected_prob_index",
    "purity_corrected_prob_ratio"
  )

  # The pre-sort Ig+ frequency is only a prior on the purity-corrected change
  # axis, and an additive constant at that, so the Z-score's centering removes
  # it -- no presorting fraction required. The prob scores do need one, and each
  # is gated (and so complains) separately, hence collecting every warning rather
  # than expect_warning()'ing one and letting the other escape the test.
  warnings_seen <- character(0)
  withCallingHandlers(
    pis <- getPhyloIgSeq(
      physeq = physeq,
      sample_id_name = "sample_id",
      fraction_id_name = "sorting_fraction",
      positive_fraction_name = "pos",
      first_negative_fraction_name = "neg1",
      second_negative_fraction_name = "neg2",
      ig_freq_name = "ig_pheno_fraction",
      ig_freq_layout = "long",
      scores = requested,
      empirical_null_distribution = FALSE
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  dropped <- grep("^Cannot compute", warnings_seen, value = TRUE)
  expect_true(any(grepl("purity_corrected_prob_index", dropped)))
  expect_true(any(grepl("purity_corrected_prob_ratio", dropped)))
  # Match the drop message specifically, not the score name anywhere: the
  # out-of-cone report also names this score, and it is not a complaint about
  # missing inputs.
  expect_false(any(grepl("purity_corrected_slide_z", dropped)))

  expect_true("purity_corrected_slide_z" %in% pis@score_names)
  expect_false(any(
    c("purity_corrected_prob_index", "purity_corrected_prob_ratio") %in%
      pis@score_names
  ))
  expect_false(all(is.na(pis@ig_coating$purity_corrected_slide_z)))
})

test_that("purity_corrected_slide_z is dropped without the long ig_freq layout", {
  physeq <- make_long_ig_freq_physeq(samples = "sample_1")

  # This fixture's ig_freq column is deliberately per-fraction, so reading it as
  # "wide" also warns that it varies within the sample -- incidental here, so
  # collect every warning rather than let it escape the test.
  warnings_seen <- character(0)
  withCallingHandlers(
    pis <- getPhyloIgSeq(
      physeq = physeq,
      sample_id_name = "sample_id",
      fraction_id_name = "sorting_fraction",
      positive_fraction_name = "pos",
      first_negative_fraction_name = "neg1",
      presorting_fraction_name = "whole",
      ig_freq_name = "ig_pheno_fraction",
      ig_freq_layout = "wide",
      scores = c("slide_z", "purity_corrected_slide_z"),
      empirical_null_distribution = FALSE
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("purity_corrected_slide_z", warnings_seen)))
  expect_false("purity_corrected_slide_z" %in% pis@score_names)
  # The uncorrected axis survives, so the pipeline still produces geometry.
  expect_equal(unique(ma_coords(pis)$change_transform), "log_ratio")
})

test_that("ig_fraction_names lists the abundance columns and nothing else", {
  pis <- run_long(
    make_long_ig_freq_physeq(samples = "sample_1"),
    scores = c("slide_z", "palm")
  )

  fractions <- ig_fraction_names(pis)
  expect_setequal(fractions, c("pos", "neg1", "neg2", "whole"))
  # The point of the accessor: no identifier, score or diagnostic column leaks in
  # the way a hand-maintained denylist over colnames(ig_coating) would let them.
  expect_false(any(
    c("taxon_id", "sample_id", "slide_z", "palm", "zeros_imputed") %in%
      fractions
  ))
  expect_true(all(fractions %in% names(pis@ig_coating)))
})

test_that("ig_fraction_names follows the multi-positive-fraction rename", {
  data("ps_phage_display", package = "PhyloIgSeq", envir = environment())
  pis <- suppressWarnings(getPhyloIgSeq(
    physeq = ps_phage_display,
    sample_id_name = "sample_id",
    fraction_id_name = "phage_fraction",
    positive_fraction_name = c("round1", "round2"),
    first_negative_fraction_name = "input",
    scores = "palm"
  ))

  fractions <- ig_fraction_names(pis)
  expect_true("positive_fraction_abundance" %in% fractions)
  expect_false(any(c("round1", "round2") %in% fractions))
  expect_true(all(fractions %in% names(pis@ig_coating)))
})

test_that("getPhyloIgSeq reports the out-of-cone rate once, cohort-wide", {
  physeq <- make_long_ig_freq_physeq(
    samples = c("sample_1", "sample_2", "sample_3")
  )

  warnings_seen <- character(0)
  withCallingHandlers(
    pis <- run_long(physeq, scores = "purity_corrected_slide_z"),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  cone <- grep("admissible cone", warnings_seen, value = TRUE)
  # Once for the whole run, not once per sample: the share across the cohort is
  # the number that matters, and per-sample repeats would bury it.
  expect_length(cone, 1)
  expect_match(cone, "^[0-9]+% of the Ig\\+ vs Ig- pairs")
  # It has to say what to do about it, or it is just noise.
  expect_match(cone, "obs_in_cone")

  # The reported share must match the slot it was computed from.
  purity <- ma_coords(pis, "purity_corrected")
  expected <- round(100 * mean(!purity$obs_in_cone, na.rm = TRUE))
  expect_match(cone, paste0("^", expected, "%"))
})

test_that("no out-of-cone report when nothing was clamped", {
  # A near-perfect sort leaves the whole cone open, so there is nothing to report.
  physeq <- make_long_ig_freq_physeq(
    presort = 0.40, pos = 0.999, neg = 0.001, samples = "sample_1"
  )
  warnings_seen <- character(0)
  withCallingHandlers(
    run_long(physeq, scores = "purity_corrected_slide_z"),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(grep("admissible cone", warnings_seen), 0)
})
