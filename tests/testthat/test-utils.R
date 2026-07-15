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
