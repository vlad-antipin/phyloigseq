library(PhyloIgSeq)

# 6 taxa x 6 samples; Rank1 groups taxa into 3 pairs (A/A, B/B, C/C) so
# tax_glom(taxrank = "Rank1") has something to agglomerate. t5 is all-zero so
# its coefficient of variation is undefined (mean = 0) -- used to exercise
# the CV guard. No phy_tree, so unifrac/wunifrac always hit the bray
# fallback. "fraction" splits samples into two groups of 3 for
# fraction_id_name/fraction_ids. "na_var" has one NA to exercise
# vars_to_remove_na.
make_heatmap_ps <- function() {
  mat <- matrix(
    c(
      50, 40, 30, 20, 10, 5,
      10, 12, 14, 16, 18, 20,
      5, 5, 5, 5, 5, 5,
      100, 90, 80, 5, 5, 5,
      0, 0, 0, 0, 0, 0,
      30, 25, 20, 15, 10, 5
    ),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(paste0("t", 1:6), paste0("s", 1:6))
  )
  tax <- matrix(
    c("A", "A", "B", "B", "C", "C"),
    nrow = 6,
    ncol = 1,
    dimnames = list(paste0("t", 1:6), "Rank1")
  )
  sdata <- data.frame(
    group = c("x", "x", "x", "y", "y", "y"),
    numeric_var = c(1, 2, 3, 4, 5, 6),
    fraction = c("f1", "f1", "f1", "f2", "f2", "f2"),
    na_var = c(NA, "a", "a", "a", "a", "a"),
    row.names = paste0("s", 1:6)
  )
  phyloseq::phyloseq(
    phyloseq::otu_table(mat, taxa_are_rows = TRUE),
    phyloseq::tax_table(tax),
    phyloseq::sample_data(sdata)
  )
}

# ---- get_phylo_heatmap ----

test_that("get_phylo_heatmap returns the expected structure", {
  ps <- make_heatmap_ps()
  hd <- get_phylo_heatmap(ps)
  expect_named(hd, c("heat_matrix", "taxa_sorted_by_abundance", "sample_data", "dendrogram"))
  expect_true(is.matrix(hd$heat_matrix))
  expect_equal(dim(hd$heat_matrix), c(6, 6))
  expect_setequal(hd$taxa_sorted_by_abundance, paste0("t", 1:6))
  expect_s3_class(hd$sample_data, "data.frame")
  expect_s3_class(hd$dendrogram, "dendrogram")
})

test_that("get_phylo_heatmap orients heat_matrix with taxa as rows regardless of input orientation", {
  ps <- make_heatmap_ps()
  hd_cols <- get_phylo_heatmap(ps) # otu_table is taxa_are_rows = TRUE on construction
  hd_rows <- get_phylo_heatmap(phyloseq::t(ps))
  expect_equal(hd_cols$heat_matrix, hd_rows$heat_matrix)
})

test_that("get_phylo_heatmap sorts taxa by decreasing total abundance", {
  ps <- make_heatmap_ps()
  hd <- get_phylo_heatmap(ps)
  totals <- phyloseq::taxa_sums(ps)[hd$taxa_sorted_by_abundance]
  expect_equal(unname(totals), sort(unname(totals), decreasing = TRUE))
})

test_that("get_phylo_heatmap filters samples by fraction_id_name/fraction_ids", {
  ps <- make_heatmap_ps()
  hd <- get_phylo_heatmap(ps, fraction_id_name = "fraction", fraction_ids = "f1")
  expect_equal(nrow(hd$sample_data), 3)
  expect_equal(ncol(hd$heat_matrix), 3)
})

test_that("get_phylo_heatmap agglomerates the heatmap matrix to taxrank_for_heatmap", {
  ps <- make_heatmap_ps()
  hd <- get_phylo_heatmap(ps, taxrank_for_heatmap = "Rank1")
  expect_equal(nrow(hd$heat_matrix), 3)
})

test_that("taxrank_for_hclust agglomerates clustering independently of taxrank_for_heatmap", {
  ps <- make_heatmap_ps()
  hd <- get_phylo_heatmap(ps, taxrank_for_hclust = "Rank1")
  # heatmap matrix unaffected (taxrank_for_heatmap is NULL)...
  expect_equal(nrow(hd$heat_matrix), 6)
  # ...but the dendrogram is still produced without error, from the
  # Rank1-agglomerated distance matrix.
  expect_s3_class(hd$dendrogram, "dendrogram")
})

test_that("get_phylo_heatmap drops samples with NA in vars_to_remove_na", {
  ps <- make_heatmap_ps()
  hd <- get_phylo_heatmap(ps, vars_to_remove_na = "na_var")
  expect_equal(nrow(hd$sample_data), 5)
  expect_equal(ncol(hd$heat_matrix), 5)
  expect_false(anyNA(hd$sample_data$na_var))
})

test_that("get_phylo_heatmap falls back to bray when distance is NULL", {
  ps <- make_heatmap_ps()
  expect_message(get_phylo_heatmap(ps, distance = NULL), "bray-curtis selected by default")
})

test_that("get_phylo_heatmap falls back to bray when unifrac is requested without a phy_tree", {
  ps <- make_heatmap_ps()
  expect_null(phyloseq::access(ps, "phy_tree"))
  expect_message(
    get_phylo_heatmap(ps, distance = "unifrac"),
    "requires a phylogenetic tree"
  )
})

test_that("get_phylo_heatmap requires a phyloseq object", {
  expect_error(get_phylo_heatmap(list(a = 1)), "phyloseq")
  expect_error(get_phylo_heatmap(matrix(1:4, 2)), "phyloseq")
})

test_that("get_phylo_heatmap requires physeq to be supplied", {
  expect_error(get_phylo_heatmap())
})

# ---- plot_phylo_heatmap ----

test_that("plot_phylo_heatmap returns a Heatmap sorted by abundance by default", {
  ps <- make_heatmap_ps()
  hd <- get_phylo_heatmap(ps)
  ht <- plot_phylo_heatmap(hd)
  expect_s4_class(ht, "Heatmap")
  expect_equal(rownames(ht@matrix), hd$taxa_sorted_by_abundance)
})

test_that("plot_phylo_heatmap caps rows to nb_top_taxa, including the nb_top_taxa = 1 edge case", {
  ps <- make_heatmap_ps()
  hd <- get_phylo_heatmap(ps)
  ht <- plot_phylo_heatmap(hd, nb_top_taxa = 3)
  expect_equal(nrow(ht@matrix), 3)

  ht1 <- plot_phylo_heatmap(hd, nb_top_taxa = 1)
  expect_true(is.matrix(ht1@matrix))
  expect_equal(nrow(ht1@matrix), 1)
})

test_that("plot_phylo_heatmap sorts by coefficient of variation and warns for zero-mean taxa", {
  ps <- make_heatmap_ps()
  hd <- get_phylo_heatmap(ps)
  expect_warning(
    ht <- plot_phylo_heatmap(hd, sort_taxa_by_diff_abundance = TRUE, nb_top_taxa = 6),
    "Coefficient of variation is undefined.*t5"
  )
  # t5 (mean = 0, undefined CV) is sorted last via na.last = TRUE
  expect_equal(rownames(ht@matrix)[6], "t5")
})

test_that("plot_phylo_heatmap sorts by Kruskal-Wallis effect size for a character variable", {
  ps <- make_heatmap_ps()
  hd <- get_phylo_heatmap(ps)
  expect_message(
    ht <- plot_phylo_heatmap(
      hd,
      sort_taxa_by_diff_abundance = TRUE,
      var_for_diff_abundance = "group",
      nb_top_taxa = 6
    ),
    "Kruskal-Wallis"
  )
  expect_s4_class(ht, "Heatmap")
})

test_that("plot_phylo_heatmap sorts by Spearman effect size for a numeric variable", {
  ps <- make_heatmap_ps()
  hd <- get_phylo_heatmap(ps)
  suppressWarnings(
    expect_message(
      ht <- plot_phylo_heatmap(
        hd,
        sort_taxa_by_diff_abundance = TRUE,
        var_for_diff_abundance = "numeric_var",
        nb_top_taxa = 6
      ),
      "Spearman"
    )
  )
  expect_s4_class(ht, "Heatmap")
})

test_that("plot_phylo_heatmap errors clearly for an unsupported var_for_diff_abundance type", {
  ps <- make_heatmap_ps()
  hd <- get_phylo_heatmap(ps)
  hd$sample_data$logical_var <- c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
  expect_error(
    plot_phylo_heatmap(hd, sort_taxa_by_diff_abundance = TRUE, var_for_diff_abundance = "logical_var"),
    "must be a character, factor, or numeric"
  )
})

test_that("plot_phylo_heatmap builds top/bottom annotations only when requested", {
  ps <- make_heatmap_ps()
  hd <- get_phylo_heatmap(ps)
  ht_none <- plot_phylo_heatmap(hd)
  expect_null(ht_none@top_annotation)
  expect_null(ht_none@bottom_annotation)

  ht_ann <- plot_phylo_heatmap(
    hd,
    top_annotation_vars = c("group", "fraction"),
    bottom_annotation_var = "numeric_var"
  )
  expect_false(is.null(ht_ann@top_annotation))
  expect_false(is.null(ht_ann@bottom_annotation))
})

test_that("plot_phylo_heatmap ignores an absent bottom_annotation_var", {
  ps <- make_heatmap_ps()
  hd <- get_phylo_heatmap(ps)
  ht <- plot_phylo_heatmap(hd, bottom_annotation_var = "not_a_real_var")
  expect_null(ht@bottom_annotation)
})

test_that("plot_phylo_heatmap scales columns when scale_cols = TRUE", {
  ps <- make_heatmap_ps()
  hd <- get_phylo_heatmap(ps)
  ht_raw <- plot_phylo_heatmap(hd)
  ht_scaled <- plot_phylo_heatmap(hd, scale_cols = TRUE)
  expect_false(isTRUE(all.equal(ht_raw@matrix, ht_scaled@matrix)))
  # each column (sample) is z-scored: mean ~0, sd ~1
  expect_equal(unname(colMeans(ht_scaled@matrix)), rep(0, ncol(ht_scaled@matrix)), tolerance = 1e-10)
})
