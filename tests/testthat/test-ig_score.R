library(PhyloIgSeq)

# ---- compute_ig_score ----

test_that("compute_ig_score computes 'palm' as the ratio of relative abundances", {
  pos <- c(50, 30, 20, 5)
  neg <- c(5, 10, 40, 45)

  result <- compute_ig_score(method = "palm", pos = pos, neg = neg)

  pos_abund <- pos / sum(pos)
  neg_abund <- neg / sum(neg)
  expect_equal(result, pos_abund / neg_abund)
})

test_that("compute_ig_score computes 'kau' per its closed-form formula", {
  pos <- c(50, 30, 20, 5)
  neg <- c(5, 10, 40, 45)

  result <- compute_ig_score(method = "kau", pos = pos, neg = neg)

  pos_abund <- pos / sum(pos)
  neg_abund <- neg / sum(neg)
  expected <- -log2(pos_abund / neg_abund) / log10(pos_abund * neg_abund)
  expect_equal(result, expected)
})

test_that("compute_ig_score computes 'prob_index' from pos and pre abundances", {
  pos <- c(50, 30, 20, 5)
  pre <- c(20, 20, 20, 40)
  presort_ig_freq <- 0.3

  result <- compute_ig_score(
    method = "prob_index",
    pos = pos,
    pre = pre,
    presort_ig_freq = presort_ig_freq
  )

  pos_abund <- pos / sum(pos)
  pre_abund <- pre / sum(pre)
  expect_equal(result, pos_abund * presort_ig_freq / pre_abund)
})

test_that("compute_ig_score computes 'prob_ratio' from pos, neg, and presort_ig_freq", {
  pos <- c(50, 30, 20, 5)
  neg <- c(5, 10, 40, 45)
  presort_ig_freq <- 0.3

  result <- compute_ig_score(
    method = "prob_ratio",
    pos = pos,
    neg = neg,
    presort_ig_freq = presort_ig_freq
  )

  pos_abund <- pos / sum(pos)
  neg_abund <- neg / sum(neg)
  expected <- log2(
    pos_abund * presort_ig_freq / (neg_abund * (1 - presort_ig_freq))
  )
  expect_equal(result, expected)
})

# The purity-corrected scores un-mix the two fractions into the really-Ig+ and
# really-Ig- populations using each fraction's own Ig+ frequency (see
# compute_ig_score()'s Details for the model). These tests re-solve that 2x2
# system here rather than reusing the package's internal helper.
#
# `prior` is the regularization added to both clamped pools. Recomputed here
# rather than taken from the package so the expectation is independent of it:
# one read's worth of composition, mapped through the un-mixing into pool space.
default_pool_prior <- function(pos, neg, p, q) {
  1 / ((p - q) * min(sum(pos), sum(neg)))
}

unmix_of <- function(pos_abund, neg_abund, p, q, prior = 0) {
  list(
    ig_pos = pmax(((1 - q) * pos_abund - (1 - p) * neg_abund) / (p - q), 0) +
      prior,
    ig_neg = pmax((p * neg_abund - q * pos_abund) / (p - q), 0) + prior
  )
}

test_that("compute_ig_score computes 'purity_corrected_prob_index' by un-mixing the fractions", {
  pos <- c(50, 30, 20, 5)
  neg <- c(5, 10, 40, 45)
  pre <- c(20, 20, 20, 40)
  presort_ig_freq <- 0.3
  pos_ig_freq <- 0.9
  neg_ig_freq <- 0.1

  result <- compute_ig_score(
    method = "purity_corrected_prob_index",
    pos = pos,
    neg = neg,
    pre = pre,
    presort_ig_freq = presort_ig_freq,
    pos_ig_freq = pos_ig_freq,
    neg_ig_freq = neg_ig_freq
  )

  pools <- unmix_of(
    pos / sum(pos),
    neg / sum(neg),
    pos_ig_freq,
    neg_ig_freq,
    prior = default_pool_prior(pos, neg, pos_ig_freq, neg_ig_freq)
  )
  expected <- presort_ig_freq *
    pools$ig_pos /
    (presort_ig_freq * pools$ig_pos + (1 - presort_ig_freq) * pools$ig_neg)
  expect_equal(result, expected)
})

test_that("`pool_prior = 0` restores the unregularized clamp", {
  pos <- c(50, 30, 20, 5)
  neg <- c(5, 10, 40, 45)
  presort_ig_freq <- 0.3
  pos_ig_freq <- 0.9
  neg_ig_freq <- 0.1

  result <- compute_ig_score(
    method = "purity_corrected_prob_index",
    pos = pos,
    neg = neg,
    presort_ig_freq = presort_ig_freq,
    pos_ig_freq = pos_ig_freq,
    neg_ig_freq = neg_ig_freq,
    pool_prior = 0
  )

  pools <- unmix_of(pos / sum(pos), neg / sum(neg), pos_ig_freq, neg_ig_freq)
  expected <- presort_ig_freq *
    pools$ig_pos /
    (presort_ig_freq * pools$ig_pos + (1 - presort_ig_freq) * pools$ig_neg)
  expect_equal(result, expected)
  # Taxa clamped entirely into one population sit exactly at the bounds without
  # regularization -- which is what makes their logit infinite, see below.
  expect_true(any(result %in% c(0, 1)))
})

test_that("regularization keeps maximally coated taxa finite instead of NA", {
  # Taxon 1 is far enough above the cone's ceiling (pos/neg > p/q) that the
  # un-mixing places it entirely in the Ig+ population.
  pos <- c(500, 30, 20, 5)
  neg <- c(1, 10, 40, 45)
  args <- list(
    pos = pos,
    neg = neg,
    presort_ig_freq = 0.3,
    pos_ig_freq = 0.9,
    neg_ig_freq = 0.1
  )

  bare <- do.call(
    compute_ig_score,
    c(list(method = "purity_corrected_prob_ratio", pool_prior = 0), args)
  )
  regularized <- do.call(
    compute_ig_score,
    c(list(method = "purity_corrected_prob_ratio"), args)
  )

  # Without regularization the most strongly coated taxon -- the one the score
  # exists to find -- is dropped as NA. With it, it saturates instead.
  expect_true(is.na(bare[1]))
  expect_true(is.finite(regularized[1]))
  # And it is still ranked above every other taxon.
  expect_equal(which.max(regularized), 1L)
})

test_that("regularization shrinks rare taxa harder than abundant ones", {
  # Two taxa with the same fold change but a 100x difference in depth. Both are
  # above the cone ceiling, so both saturate; the rare one is pulled far closer
  # to the null because the prior is a larger share of its pool.
  p <- 0.9
  q <- 0.1
  filler <- rep(1000, 4)
  abundant <- compute_ig_score(
    method = "purity_corrected_prob_ratio",
    pos = c(1000, filler),
    neg = c(1, filler),
    presort_ig_freq = 0.3,
    pos_ig_freq = p,
    neg_ig_freq = q
  )[1]
  rare <- compute_ig_score(
    method = "purity_corrected_prob_ratio",
    pos = c(10, filler),
    neg = c(1, filler),
    presort_ig_freq = 0.3,
    pos_ig_freq = p,
    neg_ig_freq = q
  )[1]

  expect_true(is.finite(abundant) && is.finite(rare))
  expect_gt(abundant, rare)
})

test_that("compute_ig_score warns and gives up regularizing compositional input", {
  pos <- c(50, 30, 20, 5)
  neg <- c(5, 10, 40, 45)

  expect_warning(
    result <- compute_ig_score(
      method = "purity_corrected_prob_ratio",
      pos = pos / sum(pos),
      neg = neg / sum(neg),
      presort_ig_freq = 0.3,
      pos_ig_freq = 0.9,
      neg_ig_freq = 0.1
    ),
    "already compositional"
  )
  # Falling back to pool_prior = 0 means the clamped taxa are NA again.
  expect_true(anyNA(result))
})

test_that("the purity-corrected scores ignore `pre` and stay in range", {
  pos <- c(50, 30, 20, 5)
  neg <- c(5, 10, 40, 45)
  args <- list(
    method = "purity_corrected_prob_index",
    pos = pos,
    neg = neg,
    presort_ig_freq = 0.3,
    pos_ig_freq = 0.9,
    neg_ig_freq = 0.1
  )

  # The pre-sort library is not always the same material as the sorted fractions
  # (it may skip the clean-up they go through), so the corrected scores must not
  # depend on it at all.
  with_pre <- do.call(compute_ig_score, c(args, list(pre = c(20, 20, 20, 40))))
  other_pre <- do.call(compute_ig_score, c(args, list(pre = c(1, 99, 1, 99))))
  expect_equal(with_pre, other_pre)
  expect_equal(with_pre, do.call(compute_ig_score, args))
  expect_true(all(with_pre >= 0 & with_pre <= 1, na.rm = TRUE))
})

test_that("the purity-corrected ratio is exactly the logit of the purity-corrected index", {
  args <- list(
    pos = c(50, 30, 20, 5),
    neg = c(5, 10, 40, 45),
    presort_ig_freq = 0.3,
    pos_ig_freq = 0.9,
    neg_ig_freq = 0.1
  )

  index <- do.call(
    compute_ig_score,
    c(list(method = "purity_corrected_prob_index"), args)
  )
  ratio <- do.call(
    compute_ig_score,
    c(list(method = "purity_corrected_prob_ratio"), args)
  )

  # A taxon clamped out of one population lands at index 0 or 1, whose logit is
  # infinite; compute_ig_score() maps that to NA for every score, so apply the
  # same mapping to the expectation.
  expected <- log2(index / (1 - index))
  expected[is.nan(expected) | is.infinite(expected)] <- NA
  expect_equal(ratio, expected)
})

test_that("the purity-corrected scores reduce to their uncorrected forms at a perfect sort", {
  pos <- c(50, 30, 20, 5)
  neg <- c(5, 10, 40, 45)
  presort_ig_freq <- 0.4
  # At p = 1, q = 0 the fractions ARE the two populations, so the pre-sort
  # composition consistent with the model is the P-weighted mixture of them.
  pre_consistent <- presort_ig_freq * pos / sum(pos) +
    (1 - presort_ig_freq) * neg / sum(neg)

  # `pool_prior = 0` because the reduction is exact only for the bare un-mixing;
  # the default read-scale prior perturbs every taxon by O(prior/pool), which is
  # the point of it. See the separate test that the two agree to within that.
  expect_equal(
    compute_ig_score(
      method = "purity_corrected_prob_index",
      pos = pos,
      neg = neg,
      presort_ig_freq = presort_ig_freq,
      pos_ig_freq = 1,
      neg_ig_freq = 0,
      pool_prior = 0
    ),
    compute_ig_score(
      method = "prob_index",
      pos = pos,
      pre = pre_consistent,
      presort_ig_freq = presort_ig_freq
    )
  )

  expect_equal(
    compute_ig_score(
      method = "purity_corrected_prob_ratio",
      pos = pos,
      neg = neg,
      presort_ig_freq = presort_ig_freq,
      pos_ig_freq = 1,
      neg_ig_freq = 0,
      pool_prior = 0
    ),
    compute_ig_score(
      method = "prob_ratio",
      pos = pos,
      neg = neg,
      presort_ig_freq = presort_ig_freq
    )
  )
})

test_that("the default pool prior perturbs a perfect sort only at the read scale", {
  pos <- c(5000, 3000, 2000, 500)
  neg <- c(500, 1000, 4000, 4500)
  presort_ig_freq <- 0.4

  bare <- compute_ig_score(
    method = "purity_corrected_prob_index",
    pos = pos,
    neg = neg,
    presort_ig_freq = presort_ig_freq,
    pos_ig_freq = 1,
    neg_ig_freq = 0,
    pool_prior = 0
  )
  regularized <- compute_ig_score(
    method = "purity_corrected_prob_index",
    pos = pos,
    neg = neg,
    presort_ig_freq = presort_ig_freq,
    pos_ig_freq = 1,
    neg_ig_freq = 0
  )

  # Every taxon here is well inside the cone and far above one read, so the prior
  # is a negligible share of each pool: the two agree to a few parts in a
  # thousand rather than exactly.
  expect_equal(regularized, bare, tolerance = 0.01)
  expect_false(isTRUE(all.equal(regularized, bare, tolerance = 1e-12)))
})

test_that("a pre-sort frequency outside the fractions' range warns but still scores", {
  pos <- c(50, 30, 20, 5)
  neg <- c(5, 10, 40, 45)
  # P > p is impossible if all three frequencies were purities of the same two
  # populations. Un-mixing uses P only as a prior, so unlike the old sort-recovery
  # clamp the score must stay defined rather than fall back to an uncorrected one.
  expect_warning(
    result <- compute_ig_score(
      method = "purity_corrected_prob_index",
      pos = pos,
      neg = neg,
      presort_ig_freq = 0.6,
      pos_ig_freq = 0.5,
      neg_ig_freq = 0.1
    ),
    "outside the range spanned by the sorted"
  )
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
  expect_false(all(is.na(result)))
})

test_that("compute_ig_score returns NA for purity-corrected scores when a frequency is missing", {
  pos <- c(50, 30, 20, 5)
  neg <- c(5, 10, 40, 45)
  pre <- c(20, 20, 20, 40)

  # No separation between the fractions leaves the recovery undefined, as does
  # an absent in-fraction frequency.
  expect_true(all(is.na(compute_ig_score(
    method = "purity_corrected_prob_index",
    pos = pos, neg = neg, pre = pre,
    presort_ig_freq = 0.3, pos_ig_freq = 0.3, neg_ig_freq = 0.3
  ))))
  expect_true(all(is.na(compute_ig_score(
    method = "purity_corrected_prob_index",
    pos = pos, neg = neg, pre = pre,
    presort_ig_freq = 0.3, pos_ig_freq = NULL, neg_ig_freq = 0.1
  ))))
})

test_that("compute_ig_score converts NaN/Inf results to NA", {
  # A zero-count taxon in neg makes neg_abund 0 for that taxon, so palm's
  # pos_abund / neg_abund is Inf; compute_ig_score() must map that to NA.
  pos <- c(10, 20)
  neg <- c(5, 0)

  result <- compute_ig_score(method = "palm", pos = pos, neg = neg)

  expect_false(is.infinite(result[2]))
  expect_true(is.na(result[2]))
})

test_that("compute_ig_score rejects an unrecognized method via match.arg", {
  expect_error(
    compute_ig_score(method = "bogus", pos = c(1, 2, 3)),
    "should be one of"
  )
})

test_that("compute_ig_score returns all-NA when a required fraction argument is missing", {
  # 'palm' needs neg; omitting it leaves neg_abund as a scalar NA, which
  # recycles to an all-NA score vector of the right length.
  result <- compute_ig_score(method = "palm", pos = c(10, 20, 30))

  expect_length(result, 3)
  expect_true(all(is.na(result)))
})

# ---- to_wider_ig_score ----

make_ig_coating_agglom <- function() {
  data.frame(
    sample_id = rep(paste0("sample_", 1:3), each = 2),
    taxon_id = rep(c("taxon_1", "taxon_2"), times = 3),
    slide_z = c(0.1, 2.2, 0.6, 1.3, -0.1, 1.5),
    palm = c(0.7, 0.4, 0.2, 0.2, 0.7, 0.1),
    not_a_score = 1:6
  )
}

test_that("to_wider_ig_score defaults scores to columns shared with IG_SCORES", {
  result <- to_wider_ig_score(make_ig_coating_agglom())

  expect_named(result, c("slide_z", "palm"))
})

test_that("to_wider_ig_score pivots each requested score to one row per sample_id and one column per taxon_id", {
  result <- to_wider_ig_score(make_ig_coating_agglom(), scores = "slide_z")

  expect_named(result, "slide_z")
  wide <- result$slide_z
  expect_equal(wide$sample_id, c("sample_1", "sample_2", "sample_3"))
  expect_setequal(colnames(wide), c("sample_id", "taxon_1", "taxon_2"))
  expect_equal(wide$taxon_1, c(0.1, 0.6, -0.1))
  expect_equal(wide$taxon_2, c(2.2, 1.3, 1.5))
})

test_that("to_wider_ig_score drops taxa not present in at least shared_by of samples", {
  ig_coating_agglom <- rbind(
    make_ig_coating_agglom(),
    data.frame(
      sample_id = "sample_1",
      taxon_id = "taxon_rare",
      slide_z = 5,
      palm = 0.9,
      not_a_score = 7
    )
  )

  # taxon_rare has a value in only 1/3 samples; a 0.5 threshold should drop it.
  result <- to_wider_ig_score(ig_coating_agglom, scores = "slide_z", shared_by = 0.5)

  expect_false("taxon_rare" %in% colnames(result$slide_z))
})

test_that("to_wider_ig_score keeps every taxon when shared_by is NULL", {
  ig_coating_agglom <- rbind(
    make_ig_coating_agglom(),
    data.frame(
      sample_id = "sample_1",
      taxon_id = "taxon_rare",
      slide_z = 5,
      palm = 0.9,
      not_a_score = 7
    )
  )

  result <- to_wider_ig_score(ig_coating_agglom, scores = "slide_z")

  expect_true("taxon_rare" %in% colnames(result$slide_z))
})

test_that("to_wider_ig_score always excludes an all-NA taxon, warning about it, even when shared_by is NULL/0", {
  ig_coating_agglom <- rbind(
    make_ig_coating_agglom(),
    data.frame(
      sample_id = "sample_1",
      taxon_id = "taxon_all_na",
      slide_z = NA_real_,
      palm = NA_real_,
      not_a_score = 7
    )
  )

  expect_warning(
    result <- to_wider_ig_score(ig_coating_agglom, scores = "slide_z"),
    "taxon_all_na"
  )
  expect_false("taxon_all_na" %in% colnames(result$slide_z))
})

# ---- plot_slide_z helpers ----

test_that(".jitter_offset computes a jitter band below the value range", {
  result <- PhyloIgSeq:::.jitter_offset(c(1, 2, 3, 4, 5))

  expect_equal(result$width, diff(range(1:5)) / 6)
  expect_equal(result$x, min(1:5) - result$width * 3)
})

test_that(".nonnegative_abundance_breaks drops breaks below zero", {
  breaks <- PhyloIgSeq:::.nonnegative_abundance_breaks(c(-4, 9))

  expect_true(all(breaks >= 0))
  expect_gte(length(breaks), 2)
})

test_that(".nonnegative_abundance_breaks keeps them when the axis is all negative", {
  # `transform_by_sample = "compositional"` puts every A below zero legitimately;
  # filtering there would leave the axis with no labels at all.
  breaks <- PhyloIgSeq:::.nonnegative_abundance_breaks(c(-13, -1))

  expect_true(any(breaks < 0))
  expect_gte(length(breaks), 2)
})

test_that("plot_ma and plot_slide_z hide negative x breaks without clipping points", {
  data("ps_igseq", package = "PhyloIgSeq", envir = environment())
  physeq <- ps_igseq
  sd <- as(phyloseq::sample_data(physeq), "data.frame")
  physeq <- phyloseq::prune_samples(
    rownames(sd)[!sd$sample_id %in% c("sample_6", "sample_7")],
    physeq
  )
  # pseudo-count imputation parks taxa below zero on the abundance axis, which is
  # what puts negative breaks on an otherwise count-scaled axis
  grouped <- suppressWarnings(group_sorted_samples(
    physeq = physeq,
    sample_id_name = "sample_id",
    fraction_id_name = "sorting_fraction",
    sample_ids = "sample_1",
    fraction_ids = c("pos", "neg1", "neg2"),
    rarefy_by_sample = TRUE,
    transform_by_sample = "identity"
  ))
  ma_plot_data <- suppressWarnings(get_ma_plot_data(
    sorted_sample_df = grouped[[1]],
    positive_fraction_name = "pos",
    first_negative_fraction_name = "neg1",
    second_negative_fraction_name = "neg2",
    zero_treatments = "pseudo_count"
  ))

  built <- ggplot2::ggplot_build(suppressWarnings(plot_ma(ma_plot_data)))
  panel <- built$layout$panel_params[[1]]
  breaks <- panel$x$breaks
  breaks <- breaks[!is.na(breaks)]
  drawn_x <- unlist(lapply(built$data, function(layer) layer$x))
  drawn_x <- drawn_x[is.finite(drawn_x)]

  expect_true(all(breaks >= 0))
  # the parked band is still drawn, below the lowest break and inside the range
  expect_lt(min(drawn_x), 0)
  expect_true(all(
    drawn_x >= panel$x.range[1] & drawn_x <= panel$x.range[2]
  ))
})

test_that(".truncate_for_tooltip leaves short values and NA untouched", {
  expect_equal(PhyloIgSeq:::.truncate_for_tooltip("short", 10), "short")
  expect_true(is.na(PhyloIgSeq:::.truncate_for_tooltip(NA, 10)))
})

test_that(".truncate_for_tooltip cuts long values short with a trailing ellipsis", {
  long <- strrep("A", 50)

  result <- PhyloIgSeq:::.truncate_for_tooltip(long, 10)

  expect_equal(result, paste0(strrep("A", 10), "..."))
})

test_that(".imputed_taxa_lookup builds one 'sample_id taxon_id' key per imputed pair", {
  lookup <- PhyloIgSeq:::.imputed_taxa_lookup(list(s1 = c(1, 3), s2 = 2))

  expect_setequal(lookup, c("s1 1", "s1 3", "s2 2"))
})

test_that(".imputed_taxa_lookup returns an empty vector for NULL/empty input", {
  expect_equal(PhyloIgSeq:::.imputed_taxa_lookup(NULL), character(0))
  expect_equal(PhyloIgSeq:::.imputed_taxa_lookup(list()), character(0))
  expect_equal(PhyloIgSeq:::.imputed_taxa_lookup(list(s1 = integer(0))), character(0))
})

test_that(".slide_z_tooltip falls back to a slide_z-only line without a tax_table", {
  ig_df <- data.frame(taxon_id = 1:2, slide_z = c(1.2345, -2))

  result <- PhyloIgSeq:::.slide_z_tooltip(ig_df, NULL, max_chars = 40)

  expect_equal(result, c("slide_z: 1.234", "slide_z: -2"))
})

test_that(".slide_z_tooltip includes every tax_table column, matched by taxon_id and truncated", {
  ig_df <- data.frame(taxon_id = c(2, 1), slide_z = c(0.5, -0.5))
  tax_table <- data.frame(
    taxon_id = 1:2,
    taxon_name = c("short", strrep("A", 20)),
    Genus = c("Alpha", "Beta")
  )

  result <- PhyloIgSeq:::.slide_z_tooltip(ig_df, tax_table, max_chars = 5)

  expect_equal(
    result[1],
    paste0("taxon_id: 2<br>taxon_name: AAAAA...<br>Genus: Beta<br>slide_z: 0.5")
  )
  expect_equal(
    result[2],
    "taxon_id: 1<br>taxon_name: short<br>Genus: Alpha<br>slide_z: -0.5"
  )
})

# ---- plot_slide_z ----

make_phyloigseq_fixture <- function(
  with_tax_table = TRUE,
  with_ellipses = TRUE,
  with_null = TRUE,
  change_transform = "log_ratio",
  # NA everywhere is what the "log_ratio" axis produces (it has no cone); set
  # FALSE entries to exercise the hollow saturated-point layer.
  obs_in_cone = NA
) {
  ig_coating <- data.frame(
    taxon_id = c(1, 2, 3, 4, 1, 2, 3, 4),
    sample_id = rep(c("s1", "s2"), each = 4),
    slide_z = c(-3, -0.5, 0.5, 3, -4, 0, 1, 4)
  )
  names(ig_coating)[3] <- names(SLIDE_Z_SCORES)[
    SLIDE_Z_SCORES == change_transform
  ]

  # MA geometry lives in its own slot, long in change_transform, keyed the same
  # way as ig_coating.
  ma_coords <- data.frame(
    sample_id = ig_coating$sample_id,
    taxon_id = ig_coating$taxon_id,
    change_transform = change_transform,
    obs_abundance = c(5, 6, 7, 8, 5, 6, 7, 8),
    obs_change = c(-1.2, -0.2, 0.2, 1.2, -1.5, 0, 0.5, 1.5),
    null_abundance = NA_real_,
    null_change = NA_real_,
    obs_in_cone = obs_in_cone,
    null_in_cone = NA,
    ellipse_level = NA_character_
  )
  if (with_null) {
    ma_coords$null_change <- c(-0.1, 0.05, -0.05, 0.1, 0, 0.1, -0.1, 0)
    ma_coords$null_abundance <- c(5, 6, 7, 8, 5, 6, 7, 8)
  }

  tax_table <- if (with_tax_table) {
    data.frame(
      taxon_id = 1:4,
      taxon_name = c("taxon_1", strrep("A", 60), "taxon_3", "taxon_4"),
      Kingdom = "Bacteria",
      Genus = c("Alpha", "Beta", "Gamma", "Delta")
    )
  } else {
    NULL
  }

  ellipse_coords <- if (with_ellipses) {
    data.frame(
      sample_id = "s1",
      x = c(1, 2, 3),
      y = c(1, 2, 1),
      ellipse_level = "0.95",
      change_transform = change_transform
    )
  } else {
    data.frame()
  }

  new(
    "PhyloIgSeq",
    ig_coating = ig_coating,
    score_names = names(ig_coating)[3],
    ma_coords = ma_coords,
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1",
    second_negative_fraction_name = "Neg2",
    ellipse_coords = ellipse_coords,
    tax_table = tax_table,
    # taxon 4 is only flagged as imputed for s1
    imputed_taxa = list(s1 = 4, s2 = integer(0))
  )
}

test_that("plot_slide_z returns a ggplot", {
  pis <- make_phyloigseq_fixture()

  expect_s3_class(suppressWarnings(plot_slide_z(pis)), "ggplot")
})

# plot_slide_z's `text` aes (for the plotly tooltip, see plotly::ggplotly(tooltip = "text") in
# the app) and its discrete size scale are intentional and trigger their own harmless ggplot2
# cosmetic warnings on every call, alongside whichever warning a given test is actually checking
# for; muffle just those two known ones so expect_warning()'s regexp isn't drowned out below.
quiet_plot_cosmetics <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (
        grepl(
          "Ignoring unknown aesthetics|for a discrete variable is not advised",
          conditionMessage(w)
        )
      ) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

test_that("plot_slide_z falls back to the observed distribution with a warning when there are no null columns", {
  pis <- make_phyloigseq_fixture(with_null = FALSE)

  expect_warning(
    quiet_plot_cosmetics(plot_slide_z(pis, ellipses = FALSE)),
    "empirical null distribution"
  )
})

test_that("plot_slide_z disables ellipses with a warning when none are furnished", {
  pis <- make_phyloigseq_fixture(with_ellipses = FALSE)

  expect_warning(quiet_plot_cosmetics(plot_slide_z(pis)), "No ellipse coordinates")
})

test_that("plot_slide_z restricts and relevels to the requested sample_ids", {
  pis <- make_phyloigseq_fixture()

  plt <- suppressWarnings(plot_slide_z(pis, sample_ids = "s2", ellipses = FALSE))

  expect_true(all(as.character(plt$data$sample_id) == "s2"))
})

test_that("plot_slide_z routes imputed taxa to the jittered layer instead of the main scatter", {
  pis <- make_phyloigseq_fixture()

  plt <- suppressWarnings(plot_slide_z(pis, ellipses = FALSE))

  expect_false(any(plt$data$sample_id == "s1" & plt$data$taxon_id == 4))
  imputed_layer_data <- plt$layers[[1]]$data
  expect_true(
    any(imputed_layer_data$sample_id == "s1" & imputed_layer_data$taxon_id == 4)
  )
})

# The cone legend has to be judged on the *rendered* plot: `scale_shape_manual(guide =)`
# alone does not settle it, because the later `guides()` call overrides the scale's own
# guide. Collect every piece of legend text from the built gtable and look for the labels.
legend_text <- function(plt) {
  gtable <- ggplot2::ggplotGrob(plt)
  boxes <- gtable$grobs[grepl("guide-box", gtable$layout$name)]
  labels <- character(0)
  collect <- function(grob) {
    if (!is.null(grob$label)) {
      labels <<- c(labels, as.character(grob$label))
    }
    for (child in c(grob$children, grob$grobs)) collect(child)
  }
  for (box in boxes) collect(box)
  labels
}

test_that("plot_slide_z hides the admissible-cone legend when nothing was clamped", {
  pis <- make_phyloigseq_fixture()

  plt <- suppressWarnings(plot_slide_z(pis, ellipses = FALSE))

  labels <- legend_text(plt)
  expect_true(any(grepl("Ig+ vs Ig-", labels, fixed = TRUE)))
  expect_false(any(grepl("cone", labels, fixed = TRUE)))
})

test_that("plot_slide_z shows the admissible-cone legend once the un-mixing clamped a taxon", {
  # Taxon 3 of s1, not taxon 4: taxon 4 is the fixture's imputed one, and imputed
  # taxa are drawn in the jitter band rather than hollow in place.
  pis <- make_phyloigseq_fixture(
    change_transform = "purity_corrected",
    obs_in_cone = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)
  )

  plt <- suppressWarnings(
    plot_slide_z(pis, ellipses = FALSE, change_transform = "purity_corrected")
  )

  labels <- legend_text(plt)
  expect_true(any(grepl("inside admissible cone", labels, fixed = TRUE)))
  expect_true(any(grepl("outside cone", labels, fixed = TRUE)))
})

test_that("plot_slide_z truncates tax_table values in the tooltip via tooltip_max_chars", {
  pis <- make_phyloigseq_fixture()

  plt <- suppressWarnings(
    plot_slide_z(pis, tooltip_max_chars = 5, ellipses = FALSE)
  )

  all_tooltips <- c(plt$data$tooltip, plt$layers[[1]]$data$tooltip)
  expect_true(any(grepl("taxon_name: AAAAA\\.\\.\\.", all_tooltips)))
})

# ---- plot_ig_score ----

test_that(".ig_score_boundary returns the documented boundary for each known score", {
  expect_equal(
    PhyloIgSeq:::.ig_score_boundary("slide_z", z_alpha2 = 1.96),
    list(left_lim = -1.96, right_lim = 1.96, midpoint = 0, left_boundary = -Inf, right_boundary = Inf)
  )
  expect_equal(
    PhyloIgSeq:::.ig_score_boundary("kau", z_alpha2 = 1.96),
    PhyloIgSeq:::.ig_score_boundary("prob_ratio", z_alpha2 = 1.96)
  )
  expect_equal(
    PhyloIgSeq:::.ig_score_boundary("palm", z_alpha2 = 1.96),
    list(left_lim = 1, right_lim = 1, midpoint = 1, left_boundary = 0, right_boundary = Inf)
  )
  expect_equal(
    PhyloIgSeq:::.ig_score_boundary("prob_index", z_alpha2 = 1.96),
    list(left_lim = 0.5, right_lim = 0.5, midpoint = 0.5, left_boundary = 0, right_boundary = 1)
  )
})

test_that(".ig_score_boundary uses the supplied z_alpha2 for slide_z", {
  result <- PhyloIgSeq:::.ig_score_boundary("slide_z", z_alpha2 = 3)
  expect_equal(c(result$left_lim, result$right_lim), c(-3, 3))
})

test_that(".ig_score_boundary errors on an unsupported score_name instead of leaving limits unbound", {
  expect_error(
    PhyloIgSeq:::.ig_score_boundary("not_a_real_score", z_alpha2 = 1.96),
    "no known plotting boundary"
  )
})

test_that(".ig_score_agglomerate: 'both' and 'taxon' agree, 'sample' diverges when sample sizes are unequal", {
  # 3 taxa nested in s1, a single (larger) taxon in s2 - a classic case where "median of
  # per-sample medians" (sample-first) differs from both the raw median ("both") and
  # "median of per-taxon medians" (taxon-first, which degenerates to the raw values here
  # since each taxon appears in only one sample).
  plot_data <- data.frame(
    sample_id = c("s1", "s1", "s1", "s2"),
    taxon_id = c(1, 2, 3, 4),
    Genus = "G",
    group = "A",
    score = c(1, 2, 3, 1000)
  )

  agg <- function(mode) {
    PhyloIgSeq:::.ig_score_agglomerate(
      plot_data,
      score_name = "score",
      score_agglom_fn = "median",
      taxrank_score = "Genus",
      taxrank_facet = NULL,
      group_score = "group",
      group_facet = NULL,
      first_score_agglom_for_each = mode
    )
  }

  expect_equal(unique(agg("both")$agglom_score), 2.5)
  expect_equal(unique(agg("taxon")$agglom_score), 2.5)
  expect_equal(unique(agg("sample")$agglom_score), 501)
})

test_that(".ig_score_agglomerate 'both'/'sample'/'taxon' agree when score_agglom_fn is 'mean'", {
  plot_data <- data.frame(
    sample_id = c("s1", "s1", "s1", "s2"),
    taxon_id = c(1, 2, 3, 4),
    Genus = "G",
    group = "A",
    score = c(1, 2, 3, 1000)
  )
  agg <- function(mode) {
    PhyloIgSeq:::.ig_score_agglomerate(
      plot_data,
      score_name = "score",
      score_agglom_fn = "mean",
      taxrank_score = "Genus",
      taxrank_facet = NULL,
      group_score = "group",
      group_facet = NULL,
      first_score_agglom_for_each = mode
    )
  }
  expect_equal(unique(agg("both")$agglom_score), mean(c(1, 2, 3, 1000)))
  expect_equal(unique(agg("taxon")$agglom_score), mean(c(1, 2, 3, 1000)))
  expect_equal(unique(agg("sample")$agglom_score), mean(c(mean(c(1, 2, 3)), 1000)))
})

test_that(".ig_score_valid_comparisons pairs up every level with >= 2 points when there are no facets", {
  plot_data <- data.frame(
    taxrank = c("A", "A", "B", "B"),
    agglom_score = 1:4
  )
  result <- PhyloIgSeq:::.ig_score_valid_comparisons(
    plot_data,
    taxrank_score = "taxrank",
    taxrank_facet = NULL,
    group_facet = NULL
  )
  expect_equal(result, list(c("A", "B")))
})

test_that(".ig_score_valid_comparisons returns NULL when fewer than 2 levels qualify", {
  plot_data <- data.frame(
    taxrank = c("A", "B", "B"),
    agglom_score = 1:3
  )
  result <- PhyloIgSeq:::.ig_score_valid_comparisons(
    plot_data,
    taxrank_score = "taxrank",
    taxrank_facet = NULL,
    group_facet = NULL
  )
  expect_null(result)
})

test_that(".ig_score_valid_comparisons excludes a level with < 2 points in one facet panel even though its overall total is >= 2", {
  # taxrank "A" has 1 point in facet f1 and 3 in f2 (4 overall); taxrank "B" has 2 in both.
  plot_data <- data.frame(
    taxrank = c("A", "A", "A", "A", "B", "B", "B", "B"),
    facet = c("f1", "f2", "f2", "f2", "f1", "f1", "f2", "f2"),
    agglom_score = 1:8
  )
  result <- PhyloIgSeq:::.ig_score_valid_comparisons(
    plot_data,
    taxrank_score = "taxrank",
    taxrank_facet = "facet",
    group_facet = NULL
  )
  # only "B" qualifies (min per-facet count 2 >= 2); "A"'s worst facet (f1) has only 1 point.
  expect_null(result)
})

test_that(".ig_score_valid_comparisons doesn't penalize a level that is simply absent from a facet panel", {
  plot_data <- data.frame(
    taxrank = c("A", "A", "B", "B", "B", "B"),
    facet = c("f2", "f2", "f1", "f1", "f2", "f2"),
    agglom_score = 1:6
  )
  result <- PhyloIgSeq:::.ig_score_valid_comparisons(
    plot_data,
    taxrank_score = "taxrank",
    taxrank_facet = "facet",
    group_facet = NULL
  )
  # "A" is entirely absent from f1 (not "too few points", just not there) and has 2 in f2.
  expect_equal(result, list(c("A", "B")))
})

make_ig_score_fixture <- function() {
  ig_coating <- data.frame(
    taxon_id = rep(1:4, times = 3),
    sample_id = rep(c("s1", "s2", "s3"), each = 4),
    slide_z = c(
      3, -0.5, 2.0, -2.5, # s1: taxa 1..4
      2.5, 0.3, 2.2, -2.0, # s2
      -3, 0.1, 1.8, -1.5 # s3
    ),
    palm = c(
      3, 1, 2, 0.2,
      2.5, 1, 2, 0.1,
      0.2, 1, 2, 0.3
    )
  )

  sample_data <- data.frame(
    sample_id = c("s1", "s2", "s3"),
    group = c("group_a", "group_a", "group_b"),
    batch = c("b1", "b2", "b1")
  )

  tax_table <- data.frame(
    taxon_id = 1:4,
    Genus = c("G1", "G1", "G2", "G2"),
    Phylum = "P1"
  )

  new(
    "PhyloIgSeq",
    ig_coating = ig_coating,
    score_names = c("slide_z", "palm"),
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1",
    sample_data = sample_data,
    tax_table = tax_table
  )
}

test_that("plot_ig_score returns a ggplot for each plot_type", {
  pis <- make_ig_score_fixture()

  for (pt in c("boxplot", "violin", "bubbleplot")) {
    plt <- suppressWarnings(
      plot_ig_score(pis, plot_type = pt, taxrank_score = "Genus", group_score = "group")
    )
    expect_s3_class(plt, "ggplot")
  }
})

test_that("plot_ig_score rejects an unrecognized plot_type/score_agglom_fn/first_score_agglom_for_each via match.arg", {
  pis <- make_ig_score_fixture()

  expect_error(plot_ig_score(pis, plot_type = "not_a_type"), "should be one of")
  expect_error(plot_ig_score(pis, score_agglom_fn = "not_a_fn"), "should be one of")
  expect_error(
    plot_ig_score(pis, first_score_agglom_for_each = "not_a_mode"),
    "should be one of"
  )
})

test_that("plot_ig_score errors for a score_name with no known plotting boundary", {
  pis <- make_ig_score_fixture()
  pis@ig_coating$custom_score <- 1

  expect_error(
    plot_ig_score(pis, score_name = "custom_score"),
    "no known plotting boundary"
  )
})

test_that("plot_ig_score applies coord_flip when transpose = TRUE", {
  pis <- make_ig_score_fixture()

  plt <- suppressWarnings(plot_ig_score(pis, taxrank_score = "Genus", group_score = "group"))
  plt_transposed <- suppressWarnings(
    plot_ig_score(pis, taxrank_score = "Genus", group_score = "group", transpose = TRUE)
  )

  expect_false(inherits(plt$coordinates, "CoordFlip"))
  expect_true(inherits(plt_transposed$coordinates, "CoordFlip"))
})

test_that("plot_ig_score adds significance brackets when at least 2 taxrank_score levels qualify", {
  pis <- make_ig_score_fixture()

  plt <- suppressWarnings(
    plot_ig_score(
      pis,
      plot_type = "boxplot",
      taxrank_score = "Genus",
      group_score = "group",
      add_stats = TRUE
    )
  )

  stat_layers <- vapply(
    plt$layers,
    function(l) inherits(l$stat, "StatSignif") || inherits(l$stat, "StatCompareMeans"),
    logical(1)
  )
  expect_true(any(stat_layers))
})

test_that("plot_ig_score silently disables add_stats (no error) when fewer than 2 taxrank_score levels have data", {
  pis <- make_ig_score_fixture()
  pis@ig_coating <- pis@ig_coating[pis@ig_coating$taxon_id %in% c(1, 2), ]

  plt <- suppressWarnings(
    plot_ig_score(
      pis,
      plot_type = "boxplot",
      taxrank_score = "Genus",
      group_score = "group",
      add_stats = TRUE
    )
  )
  expect_s3_class(plt, "ggplot")

  stat_layers <- vapply(
    plt$layers,
    function(l) inherits(l$stat, "StatSignif") || inherits(l$stat, "StatCompareMeans"),
    logical(1)
  )
  expect_false(any(stat_layers))
})

test_that("plot_ig_score accepts taxrank_facet/group_facet without error", {
  pis <- make_ig_score_fixture()

  plt <- suppressWarnings(
    plot_ig_score(
      pis,
      plot_type = "bubbleplot",
      taxrank_score = "Genus",
      taxrank_facet = "Phylum",
      group_score = "group",
      group_facet = "batch"
    )
  )
  expect_s3_class(plt, "ggplot")
})

# ---- .central_tendency ----

test_that(".central_tendency computes mean/median with NAs dropped", {
  x <- c(1, 2, 3, NA)
  expect_equal(.central_tendency(x, "mean"), mean(x, na.rm = TRUE))
  expect_equal(.central_tendency(x, "median"), median(x, na.rm = TRUE))
})

test_that(".central_tendency weight_by_abund computes a weighted mean and requires weights", {
  x <- c(1, 3)
  w <- c(10, 20)
  expect_equal(.central_tendency(x, "weight_by_abund", weights = w), weighted.mean(x, w))
  expect_error(
    .central_tendency(x, "weight_by_abund"),
    "Need abundance fraction"
  )
})

test_that(".central_tendency errors on an unrecognized method", {
  expect_error(.central_tendency(1:3, "bogus"), "wrong agglomeration method")
})

# ---- agglomPhyloIgSeq ----

make_agglom_fixture <- function(total_reads = NULL) {
  ig_coating <- data.frame(
    taxon_id = rep(1:4, times = 2),
    sample_id = rep(c("s1", "s2"), each = 4),
    # taxa 1-2 -> Genus G1, taxa 3-4 -> Genus G2 (see tax_table below)
    slide_z = c(1, 3, 5, 7, 2, 4, 6, 8),
    palm = NA_real_, # all-NA -> should be dropped from score_names/output
    Pos = c(10, 20, 30, 5, 15, 25, 35, 10)
  )
  tax_table <- data.frame(
    Kingdom = "Bacteria",
    Genus = c("G1", "G1", "G2", "G2"),
    taxon_id = 1:4
  )
  new(
    "PhyloIgSeq",
    ig_coating = ig_coating,
    score_names = c("slide_z", "palm"),
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1",
    tax_table = tax_table,
    total_reads = total_reads
  )
}

test_that("agglomPhyloIgSeq drops score columns that are all-NA from score_names and ig_coating", {
  agg <- agglomPhyloIgSeq(make_agglom_fixture(), taxrank = "Genus", agglom_method = "median")

  expect_equal(agg@score_names, "slide_z")
  expect_false("palm" %in% colnames(agg@ig_coating))
})

test_that("agglomPhyloIgSeq agglomerates scores per (taxrank, sample_id) via median/mean", {
  agg_median <- agglomPhyloIgSeq(
    make_agglom_fixture(),
    taxrank = "Genus",
    agglom_method = "median"
  )
  d <- as.data.frame(agg_median@ig_coating)
  d <- d[order(d$sample_id, d$taxon_id), ]

  expect_equal(d$taxon_id, c("G1", "G2", "G1", "G2"))
  expect_equal(d$sample_id, c("s1", "s1", "s2", "s2"))
  # medians of {1, 3}, {5, 7}, {2, 4}, {6, 8}
  expect_equal(d$slide_z, c(2, 6, 3, 7))

  agg_mean <- agglomPhyloIgSeq(make_agglom_fixture(), taxrank = "Genus", agglom_method = "mean")
  d_mean <- as.data.frame(agg_mean@ig_coating)
  d_mean <- d_mean[order(d_mean$sample_id, d_mean$taxon_id), ]
  expect_equal(d_mean$slide_z, c(2, 6, 3, 7))
})

test_that("agglomPhyloIgSeq weight_by_abund weights by the pre-agglomeration abundance_fraction", {
  agg <- agglomPhyloIgSeq(
    make_agglom_fixture(),
    taxrank = "Genus",
    abundance_fraction = "Pos",
    agglom_method = "weight_by_abund"
  )
  d <- as.data.frame(agg@ig_coating)
  d <- d[order(d$sample_id, d$taxon_id), ]

  # G1/s1: values c(1, 3), weights c(10, 20); G2/s1: values c(5, 7), weights c(30, 5)
  expect_equal(
    d$slide_z,
    c(
      weighted.mean(c(1, 3), c(10, 20)),
      weighted.mean(c(5, 7), c(30, 5)),
      weighted.mean(c(2, 4), c(15, 25)),
      weighted.mean(c(6, 8), c(35, 10))
    )
  )
  # abundance_fraction column itself is agglomerated by sum, not by agglom_method
  expect_equal(d$Pos, c(30, 35, 40, 45))
})

test_that("agglomPhyloIgSeq rejects an unrecognized agglom_method", {
  expect_error(
    agglomPhyloIgSeq(make_agglom_fixture(), taxrank = "Genus", agglom_method = "bogus"),
    "wrong agglomeration method"
  )
})

test_that("agglomPhyloIgSeq defaults agglom_method to median with a warning when NULL", {
  expect_warning(
    agg <- agglomPhyloIgSeq(make_agglom_fixture(), taxrank = "Genus"),
    "agglomeration method is set to median"
  )
  d <- as.data.frame(agg@ig_coating)
  d <- d[order(d$sample_id, d$taxon_id), ]
  expect_equal(d$slide_z, c(2, 6, 3, 7))
})

test_that("agglomPhyloIgSeq drops taxa below abundance_quantile", {
  agg <- agglomPhyloIgSeq(
    make_agglom_fixture(),
    taxrank = "Genus",
    abundance_fraction = "Pos",
    agglom_method = "mean",
    abundance_quantile = 0.5
  )
  d <- as.data.frame(agg@ig_coating)

  # per sample, only the higher-abundance Genus (G2: 35 vs G1: 30 for s1; 45 vs 40 for s2)
  # clears the per-sample median quantile threshold
  expect_true(all(d$taxon_id == "G2"))
  expect_equal(sort(d$sample_id), c("s1", "s2"))
})

test_that("agglomPhyloIgSeq's min_rel_abundance uses the matching total_reads slot over the fallback sum", {
  # total_reads covers "Pos" (the chosen abundance_fraction) with values far above the
  # per-taxon abundances, so every taxon fails a 0.5 threshold -> 0 rows
  pis_with_total_reads <- make_agglom_fixture(
    total_reads = data.frame(sample_id = c("s1", "s2"), Pos = c(1000, 2000))
  )
  agg_matching <- agglomPhyloIgSeq(
    pis_with_total_reads,
    taxrank = "Genus",
    abundance_fraction = "Pos",
    agglom_method = "mean",
    min_rel_abundance = 0.5
  )
  expect_equal(nrow(agg_matching@ig_coating), 0)

  # without a matching total_reads, the fallback sums the (already taxrank-agglomerated)
  # abundance_fraction column itself, so G2 (35/45) clears 0.5 * (30+35)/(40+45) while G1 doesn't
  agg_fallback <- agglomPhyloIgSeq(
    make_agglom_fixture(),
    taxrank = "Genus",
    abundance_fraction = "Pos",
    agglom_method = "mean",
    min_rel_abundance = 0.5
  )
  d <- as.data.frame(agg_fallback@ig_coating)
  expect_true(all(d$taxon_id == "G2"))
})

test_that("agglomPhyloIgSeq restricts total_reads to samples remaining in the agglomerated ig_coating", {
  pis_with_total_reads <- make_agglom_fixture(
    total_reads = data.frame(sample_id = c("s1", "s2"), Pos = c(1000, 2000))
  )
  agg <- agglomPhyloIgSeq(
    pis_with_total_reads,
    taxrank = "Genus",
    abundance_fraction = "Pos",
    agglom_method = "mean",
    min_rel_abundance = 0.5 # drops every taxon for both samples, see test above
  )

  expect_equal(nrow(agg@total_reads), 0)
})

test_that("agglomPhyloIgSeq defaults taxrank to taxon_id (no cross-taxon agglomeration)", {
  agg <- agglomPhyloIgSeq(make_agglom_fixture(), agglom_method = "mean")
  d <- as.data.frame(agg@ig_coating)

  expect_equal(nrow(d), 8)
  expect_setequal(d$taxon_id, as.character(1:4))
})

test_that("agglomPhyloIgSeq's make_unique_taxonomy disambiguates colliding taxon names before agglomerating", {
  tax_table <- data.frame(
    Kingdom = c("Bacteria", "Archaea"), # different lineages, same Genus name below
    Genus = c("Foo", "Foo"),
    taxon_id = 1:2
  )
  pis <- new(
    "PhyloIgSeq",
    ig_coating = data.frame(
      taxon_id = c(1, 2),
      sample_id = c("s1", "s1"),
      slide_z = c(1, 5)
    ),
    score_names = "slide_z",
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1",
    tax_table = tax_table
  )

  agg_unique <- agglomPhyloIgSeq(
    pis,
    taxrank = "Genus",
    agglom_method = "mean",
    make_unique_taxonomy = TRUE
  )
  expect_setequal(agg_unique@ig_coating$taxon_id, c("Foo", "Foo.1"))

  agg_not_unique <- agglomPhyloIgSeq(
    pis,
    taxrank = "Genus",
    agglom_method = "mean",
    make_unique_taxonomy = FALSE
  )
  expect_equal(nrow(agg_not_unique@ig_coating), 1)
  expect_equal(agg_not_unique@ig_coating$slide_z, mean(c(1, 5)))
})
