library(PhyloIgSeq)

# ---- filter_reads ----

# 2 taxa x 3 samples, engineered so that filtering order changes the result:
# taxa_sums = t1: 205, t2: 110; sample_sums (both taxa) = 105/105/105 (tied).
# min_taxa_sum = 150 always drops t2. With taxa_first = TRUE, dropping t2
# first changes sample sums to 100/100/5 (from only t1), so min_sample_sum =
# 50 then also drops s3. With taxa_first = FALSE, sample sums are checked
# against the *original* (both-taxa) totals, all tied at 105, so no sample is
# dropped before t2 is removed.
make_filter_ps <- function() {
  mat <- matrix(
    c(
      100, 100, 5,
      5, 5, 100
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("t1", "t2"), c("s1", "s2", "s3"))
  )
  sdata <- data.frame(dummy = c("a", "b", "c"), row.names = colnames(mat))
  # phyloseq() with a single component returns it unwrapped (not a
  # "phyloseq"-class object), so a sample_data component is required here.
  phyloseq::phyloseq(
    phyloseq::otu_table(mat, taxa_are_rows = TRUE),
    phyloseq::sample_data(sdata)
  )
}

test_that("filter_reads with taxa_first = TRUE filters taxa then recomputed sample sums", {
  ps <- make_filter_ps()
  out <- filter_reads(ps, min_sample_sum = 50, min_taxa_sum = 150, taxa_first = TRUE)
  expect_equal(phyloseq::taxa_names(out), "t1")
  expect_equal(phyloseq::sample_names(out), c("s1", "s2"))
})

test_that("filter_reads with taxa_first = FALSE checks samples against original sums", {
  ps <- make_filter_ps()
  out <- filter_reads(ps, min_sample_sum = 50, min_taxa_sum = 150, taxa_first = FALSE)
  expect_equal(phyloseq::taxa_names(out), "t1")
  expect_equal(phyloseq::sample_names(out), c("s1", "s2", "s3"))
})

test_that("filter_reads warns and returns NULL when every taxon is filtered out", {
  ps <- make_filter_ps()
  expect_warning(
    result <- filter_reads(ps, min_sample_sum = 0, min_taxa_sum = 1000),
    "All taxa filtered out"
  )
  expect_null(result)
})

test_that("filter_reads warns and returns NULL when every sample is filtered out", {
  ps <- make_filter_ps()
  expect_warning(
    result <- filter_reads(ps, min_sample_sum = 1000, min_taxa_sum = 0),
    "All samples filtered out"
  )
  expect_null(result)
})

test_that("filter_reads skips filtering entirely when either threshold is NULL", {
  ps <- make_filter_ps()
  out1 <- filter_reads(ps, min_sample_sum = NULL, min_taxa_sum = 150)
  out2 <- filter_reads(ps, min_sample_sum = 50, min_taxa_sum = NULL)
  expect_equal(phyloseq::ntaxa(out1), 2)
  expect_equal(phyloseq::nsamples(out1), 3)
  expect_equal(phyloseq::ntaxa(out2), 2)
  expect_equal(phyloseq::nsamples(out2), 3)
})

test_that("filter_reads reorients a taxa-are-rows input to taxa-are-columns", {
  ps <- make_filter_ps()
  expect_true(phyloseq::taxa_are_rows(ps))
  out <- filter_reads(ps, min_sample_sum = 0, min_taxa_sum = 0)
  expect_false(phyloseq::taxa_are_rows(out))
})

test_that("filter_reads preserves a sparse_otu_table through actual filtering", {
  ps <- as_sparse_phyloseq(make_filter_ps())
  out <- filter_reads(ps, min_sample_sum = 50, min_taxa_sum = 150, taxa_first = TRUE)
  expect_s4_class(phyloseq::otu_table(out), "sparse_otu_table")
})

test_that("filter_reads returns NULL for non-phyloseq input", {
  expect_null(filter_reads(list(a = 1)))
  expect_null(filter_reads(matrix(1:4, 2)))
})

# ---- plot_reads ----

test_that("plot_reads returns a faceted ggplot with no threshold lines by default", {
  ps <- make_filter_ps()
  p <- plot_reads(ps)
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  # geom_histogram() is implemented on top of GeomBar, there is no GeomHistogram
  expect_true("GeomBar" %in% layer_classes)
  expect_false("GeomVline" %in% layer_classes)
})

test_that("plot_reads draws threshold reference lines when both thresholds are given", {
  ps <- make_filter_ps()
  p <- plot_reads(ps, min_sample_sum = 50, min_taxa_sum = 150)
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomVline" %in% layer_classes)
})

test_that("plot_reads omits threshold lines when only one threshold is given", {
  ps <- make_filter_ps()
  p <- plot_reads(ps, min_sample_sum = 50, min_taxa_sum = NA)
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_false("GeomVline" %in% layer_classes)
})

test_that("plot_reads treats NULL thresholds the same as NA", {
  ps <- make_filter_ps()
  p <- plot_reads(ps, min_sample_sum = NULL, min_taxa_sum = NULL)
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_false("GeomVline" %in% layer_classes)
})

test_that("plot_reads returns NULL for non-phyloseq input", {
  expect_null(plot_reads(list(a = 1)))
})
