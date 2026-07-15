library(PhyloIgSeq)

# ---- smart_facet_ncol ----

test_that("smart_facet_ncol returns n unchanged for n <= 3", {
  expect_equal(PhyloIgSeq:::smart_facet_ncol(1), 1)
  expect_equal(PhyloIgSeq:::smart_facet_ncol(2), 2)
  expect_equal(PhyloIgSeq:::smart_facet_ncol(3), 3)
})

test_that("smart_facet_ncol picks the smallest perfect-square-ish divisor for composite n", {
  # 12 = 3 x 4 (also 2 x 6, 1 x 12); the tied-widest perfect factorization
  # with ncol >= sqrt(12) is 4, which is also the smallest such divisor.
  expect_equal(PhyloIgSeq:::smart_facet_ncol(12), 4)
  # 9 = 3 x 3, a perfect square.
  expect_equal(PhyloIgSeq:::smart_facet_ncol(9), 3)
})

test_that("smart_facet_ncol falls back to minimizing empty cells for prime/awkward n", {
  ncol <- PhyloIgSeq:::smart_facet_ncol(7)
  nrow <- ceiling(7 / ncol)
  empties <- nrow * ncol - 7
  # No other ncol in [ceiling(sqrt(7)), ceiling(7/2)] should produce fewer
  # empty cells, and ties should prefer the widest (largest) ncol.
  candidates <- seq(ceiling(sqrt(7)), ceiling(7 / 2))
  all_empties <- ceiling(7 / candidates) * candidates - 7
  expect_equal(empties, min(all_empties))
  expect_equal(ncol, max(candidates[all_empties == min(all_empties)]))
})

test_that("smart_facet_ncol never produces empty rows/columns beyond what's needed", {
  for (n in 4:40) {
    ncol <- PhyloIgSeq:::smart_facet_ncol(n)
    expect_true(ncol >= 1 && ncol <= n)
    expect_true(ceiling(n / ncol) * ncol - n < ncol)
  }
})

# ---- IG_SCORES ----

test_that("IG_SCORES is the expected character vector", {
  expect_type(IG_SCORES, "character")
  expect_setequal(
    IG_SCORES,
    c("slide_z", "palm", "kau", "prob_index", "prob_ratio")
  )
})

# ---- reverseASV ----

test_that("reverseASV transposes a taxa-are-rows phyloseq object", {
  mat <- matrix(1:6, nrow = 2, dimnames = list(c("t1", "t2"), c("s1", "s2", "s3")))
  ps <- phyloseq::phyloseq(phyloseq::otu_table(mat, taxa_are_rows = TRUE))
  result <- PhyloIgSeq:::reverseASV(ps)
  expect_false(phyloseq::taxa_are_rows(result))
  expect_equal(phyloseq::ntaxa(result), 2)
  expect_equal(phyloseq::nsamples(result), 3)
})

test_that("reverseASV leaves a taxa-are-columns phyloseq object unchanged", {
  mat <- matrix(1:6, nrow = 2, dimnames = list(c("t1", "t2"), c("s1", "s2", "s3")))
  ps <- phyloseq::phyloseq(phyloseq::otu_table(mat, taxa_are_rows = FALSE))
  result <- PhyloIgSeq:::reverseASV(ps)
  expect_false(phyloseq::taxa_are_rows(result))
  expect_identical(as(phyloseq::otu_table(result), "matrix"), as(phyloseq::otu_table(ps), "matrix"))
})

# ---- geom_jitter ----

test_that("geom_jitter defaults height to 0 when position is not supplied", {
  layer <- PhyloIgSeq:::geom_jitter()
  expect_s3_class(layer, "Layer")
  expect_equal(layer$position$height, 0)
})

test_that("geom_jitter respects an explicitly supplied position", {
  layer <- PhyloIgSeq:::geom_jitter(position = ggplot2::position_jitter(height = 0.5))
  expect_s3_class(layer, "Layer")
  expect_equal(layer$position$height, 0.5)
})

# ---- %in_interval% ----

test_that("%in_interval% detects values inside/outside a sorted interval", {
  expect_true(3.5 %in_interval% c(3, 4))
  expect_true(3 %in_interval% c(3, 4))
  expect_true(4 %in_interval% c(3, 4))
  expect_false(2.9 %in_interval% c(3, 4))
  expect_false(4.1 %in_interval% c(3, 4))
})

test_that("%in_interval% normalizes an unsorted interval", {
  expect_true(3.5 %in_interval% c(4, 3))
})

test_that("%in_interval% is vectorized over x", {
  expect_equal(c(1, 3.5, 10) %in_interval% c(3, 4), c(FALSE, TRUE, FALSE))
})

# ---- transform_abundances ----

test_that("transform_abundances applies the compositional transform column-wise", {
  mat <- matrix(c(1, 3, 6, 2, 2, 4), nrow = 2, byrow = TRUE)
  result <- PhyloIgSeq:::transform_abundances(mat, transform = "compositional")
  expect_equal(colSums(result), c(1, 1, 1))
  expect_equal(result[1, ], c(1 / 3, 3 / 5, 6 / 10))
})

test_that("transform_abundances applies the hellinger transform column-wise", {
  mat <- matrix(c(1, 3, 6, 2, 2, 4), nrow = 2, byrow = TRUE)
  result <- PhyloIgSeq:::transform_abundances(mat, transform = "hellinger")
  expect_equal(colSums(result^2), c(1, 1, 1))
})

test_that("transform_abundances transposes for taxa_are_rows = FALSE and transposes back", {
  mat <- matrix(c(1, 3, 6, 2, 2, 4), nrow = 2, byrow = TRUE)
  by_row <- PhyloIgSeq:::transform_abundances(mat, transform = "compositional")
  by_col <- PhyloIgSeq:::transform_abundances(t(mat), transform = "compositional", taxa_are_rows = FALSE)
  expect_equal(by_col, t(by_row))
})

test_that("transform_abundances defaults to compositional when transform is omitted", {
  mat <- matrix(c(1, 3, 6, 2, 2, 4), nrow = 2, byrow = TRUE)
  result <- PhyloIgSeq:::transform_abundances(mat)
  expect_equal(result, PhyloIgSeq:::transform_abundances(mat, transform = "compositional"))
})

test_that("transform_abundances errors clearly on an invalid transform", {
  mat <- matrix(c(1, 3, 6, 2, 2, 4), nrow = 2, byrow = TRUE)
  expect_error(
    PhyloIgSeq:::transform_abundances(mat, transform = "clr"),
    "should be one of"
  )
})

# ---- plot_rarefaction / plot_seq_depth ----

make_seq_depth_ps <- function() {
  mat <- matrix(
    c(
      1, 20, 30, 40, 0,
      5, 15, 25, 35, 0,
      8, 12, 22, 32, 0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("t1", "t2", "t3"), c("s1", "s2", "s3", "s4", "s5"))
  )
  sdata <- data.frame(
    group = c("A", "A", "B", "B", "B"),
    row.names = colnames(mat)
  )
  phyloseq::phyloseq(
    phyloseq::otu_table(mat, taxa_are_rows = TRUE),
    phyloseq::sample_data(sdata)
  )
}

test_that("plot_rarefaction returns a ggplot excluding zero-count samples", {
  ps <- make_seq_depth_ps()
  p <- plot_rarefaction(ps, step = 5)
  expect_s3_class(p, "ggplot")
  expect_setequal(unique(p$data$Sample), c("s1", "s2", "s3", "s4"))
})

test_that("plot_rarefaction toggles the legend via show_legend", {
  ps <- make_seq_depth_ps()
  p_legend <- plot_rarefaction(ps, step = 5, show_legend = TRUE)
  p_no_legend <- plot_rarefaction(ps, step = 5, show_legend = FALSE)
  expect_equal(p_legend$theme$legend.position, "right")
  expect_equal(p_no_legend$theme$legend.position, "none")
})

test_that("plot_rarefaction does not open a graphics device as a side effect", {
  ps <- make_seq_depth_ps()
  devs_before <- grDevices::dev.list()
  plot_rarefaction(ps, step = 5)
  expect_identical(grDevices::dev.list(), devs_before)
})

test_that("plot_seq_depth bar chart reports per-sample sequencing depth", {
  ps <- make_seq_depth_ps()
  p <- plot_seq_depth(ps)
  expect_s3_class(p, "ggplot")
  expect_equal(sort(p$data$Depth), sort(unname(phyloseq::sample_sums(ps))))
})

test_that("plot_seq_depth box type requires x_var", {
  ps <- make_seq_depth_ps()
  expect_error(plot_seq_depth(ps, type = "box"), "x_var")
})

test_that("plot_seq_depth box type facets when facet_var is supplied", {
  ps <- make_seq_depth_ps()
  p <- plot_seq_depth(ps, type = "box", x_var = "group")
  expect_true(is.null(p$facet) || inherits(p$facet, "FacetNull"))
  p_facet <- plot_seq_depth(ps, type = "box", x_var = "group", facet_var = "group")
  expect_s3_class(p_facet$facet, "FacetWrap")
})

test_that("plot_seq_depth errors on an invalid type", {
  ps <- make_seq_depth_ps()
  expect_error(plot_seq_depth(ps, type = "pie"))
})

# ---- rarefy_abundances ----

make_rarefy_matrix <- function() {
  matrix(
    c(
      10, 20, 30, 0,
      5, 15, 10, 0,
      15, 25, 20, 0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("t1", "t2", "t3"), c("s1", "s2", "s3", "s4"))
  )
}

test_that("rarefy_abundances rarefies every retained sample to common_count_sum", {
  set.seed(1)
  mat <- make_rarefy_matrix()
  result <- rarefy_abundances(mat, silent_warnings = TRUE)
  expect_true(all(colSums(result) == 30))
})

test_that("rarefy_abundances defaults common_count_sum to the smallest nonzero sample total", {
  set.seed(1)
  mat <- make_rarefy_matrix()
  result <- rarefy_abundances(mat, silent_warnings = TRUE)
  expect_setequal(colnames(result), c("s1", "s2", "s3"))
})

test_that("rarefy_abundances drops under-depth samples by default", {
  set.seed(1)
  mat <- make_rarefy_matrix()
  result <- rarefy_abundances(mat, silent_warnings = TRUE)
  expect_false("s4" %in% colnames(result))
})

test_that("rarefy_abundances keeps under-depth samples when trim_samples = FALSE", {
  set.seed(1)
  mat <- make_rarefy_matrix()
  result <- rarefy_abundances(mat, trim_samples = FALSE, silent_warnings = TRUE)
  expect_true("s4" %in% colnames(result))
  expect_equal(unname(colSums(result)["s4"]), 0)
})

test_that("rarefy_abundances honors an explicit common_count_sum", {
  set.seed(1)
  mat <- make_rarefy_matrix()
  result <- rarefy_abundances(mat, common_count_sum = 20, silent_warnings = TRUE)
  expect_true(all(colSums(result) == 20))
  expect_setequal(colnames(result), c("s1", "s2", "s3"))
})

test_that("rarefy_abundances trims taxa left all-zero after rarefaction by default", {
  set.seed(1)
  mat <- make_rarefy_matrix()
  result <- rarefy_abundances(mat, silent_warnings = TRUE)
  expect_true(all(rowSums(result) > 0))
})

test_that("rarefy_abundances trim_taxa = FALSE preserves every taxon row", {
  set.seed(1)
  mat <- make_rarefy_matrix()
  result <- rarefy_abundances(mat, trim_taxa = FALSE, silent_warnings = TRUE)
  expect_equal(nrow(result), nrow(mat))
})

test_that("rarefy_abundances handles taxa_are_rows = FALSE", {
  set.seed(1)
  mat <- t(make_rarefy_matrix())
  result <- rarefy_abundances(mat, taxa_are_rows = FALSE, silent_warnings = TRUE)
  expect_true(all(rowSums(result) == 30))
  expect_false("s4" %in% rownames(result))
})

test_that("rarefy_abundances warns about trimmed samples/taxa unless silenced", {
  set.seed(1)
  mat <- make_rarefy_matrix()
  expect_warning(rarefy_abundances(mat), "Trimmed")
  expect_no_warning(rarefy_abundances(mat, silent_warnings = TRUE))
})

test_that("rarefy_abundances leaves an all-zero sample all-zero rather than up-sampling it", {
  set.seed(1)
  mat <- make_rarefy_matrix()
  result <- rarefy_abundances(mat, trim_samples = FALSE, silent_warnings = TRUE)
  expect_equal(unname(result[, "s4"]), c(0, 0, 0))
})
