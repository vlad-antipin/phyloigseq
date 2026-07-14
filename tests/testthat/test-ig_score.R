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
  ig_freq <- 0.3

  result <- compute_ig_score(method = "prob_index", pos = pos, pre = pre, ig_freq = ig_freq)

  pos_abund <- pos / sum(pos)
  pre_abund <- pre / sum(pre)
  expect_equal(result, pos_abund * ig_freq / pre_abund)
})

test_that("compute_ig_score computes 'prob_ratio' from pos, neg, and ig_freq", {
  pos <- c(50, 30, 20, 5)
  neg <- c(5, 10, 40, 45)
  ig_freq <- 0.3

  result <- compute_ig_score(method = "prob_ratio", pos = pos, neg = neg, ig_freq = ig_freq)

  pos_abund <- pos / sum(pos)
  neg_abund <- neg / sum(neg)
  expected <- log2(pos_abund * ig_freq / (neg_abund * (1 - ig_freq)))
  expect_equal(result, expected)
})

test_that("compute_ig_score computes 'purity_corrected_prob_index' from all purity/fraction inputs", {
  pos <- c(50, 30, 20, 5)
  neg <- c(5, 10, 40, 45)
  pre <- c(20, 20, 20, 40)

  result <- compute_ig_score(
    method = "purity_corrected_prob_index",
    pos = pos,
    neg = neg,
    pre = pre,
    pos_purity = 0.9,
    neg_impurity = 0.1,
    pos_fraction = 0.4,
    neg_fraction = 0.6
  )

  pos_abund <- pos / sum(pos)
  neg_abund <- neg / sum(neg)
  pre_abund <- pre / sum(pre)
  expected <- (pos_abund * 0.9 * 0.4 + neg_abund * 0.1 * 0.6) / pre_abund
  expect_equal(result, expected)
})

test_that("compute_ig_score computes 'purity_corrected_prob_ratio' from all purity/fraction inputs", {
  pos <- c(50, 30, 20, 5)
  neg <- c(5, 10, 40, 45)

  result <- compute_ig_score(
    method = "purity_corrected_prob_ratio",
    pos = pos,
    neg = neg,
    pos_purity = 0.9,
    neg_impurity = 0.1,
    pos_fraction = 0.4,
    neg_fraction = 0.6
  )

  pos_abund <- pos / sum(pos)
  neg_abund <- neg / sum(neg)
  prob <- pos_abund * 0.9 * 0.4 + neg_abund * 0.1 * 0.6
  expect_equal(result, log2(prob / (1 - prob)))
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

# ---- plot_slide_z helpers ----

test_that(".jitter_offset computes a jitter band below the value range", {
  result <- PhyloIgSeq:::.jitter_offset(c(1, 2, 3, 4, 5))

  expect_equal(result$width, diff(range(1:5)) / 6)
  expect_equal(result$x, min(1:5) - result$width * 3)
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
  with_null = TRUE
) {
  ig_coating <- data.frame(
    taxon_id = c(1, 2, 3, 4, 1, 2, 3, 4),
    sample_id = rep(c("s1", "s2"), each = 4),
    slide_z = c(-3, -0.5, 0.5, 3, -4, 0, 1, 4),
    obs_change = c(-1.2, -0.2, 0.2, 1.2, -1.5, 0, 0.5, 1.5),
    obs_abundance = c(5, 6, 7, 8, 5, 6, 7, 8)
  )
  if (with_null) {
    ig_coating$null_change <- c(-0.1, 0.05, -0.05, 0.1, 0, 0.1, -0.1, 0)
    ig_coating$null_abundance <- c(5, 6, 7, 8, 5, 6, 7, 8)
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
      ellipse_level = "0.95"
    )
  } else {
    data.frame()
  }

  new(
    "PhyloIgSeq",
    ig_coating = ig_coating,
    score_names = "slide_z",
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

test_that("plot_slide_z falls back to the observed distribution with a warning when there are no null columns", {
  pis <- make_phyloigseq_fixture(with_null = FALSE)

  expect_warning(
    plot_slide_z(pis, ellipses = FALSE),
    "empirical null distribution"
  )
})

test_that("plot_slide_z disables ellipses with a warning when none are furnished", {
  pis <- make_phyloigseq_fixture(with_ellipses = FALSE)

  expect_warning(plot_slide_z(pis), "No ellipse coordinates")
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

test_that("plot_slide_z truncates tax_table values in the tooltip via tooltip_max_chars", {
  pis <- make_phyloigseq_fixture()

  plt <- suppressWarnings(
    plot_slide_z(pis, tooltip_max_chars = 5, ellipses = FALSE)
  )

  all_tooltips <- c(plt$data$tooltip, plt$layers[[1]]$data$tooltip)
  expect_true(any(grepl("taxon_name: AAAAA\\.\\.\\.", all_tooltips)))
})
