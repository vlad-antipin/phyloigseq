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
