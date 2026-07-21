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

# ---- is_count_like ----

test_that("is_count_like returns TRUE for a non-negative integer matrix", {
  expect_true(is_count_like(matrix(1:6, nrow = 2)))
})

test_that("is_count_like returns FALSE for non-integer values", {
  expect_false(is_count_like(matrix(c(0.5, 1, 2, 3), nrow = 2)))
})

test_that("is_count_like returns FALSE for negative values", {
  expect_false(is_count_like(matrix(c(-1, 2, 3, 4), nrow = 2)))
})

test_that("is_count_like returns FALSE for non-numeric input", {
  expect_false(is_count_like(matrix(c("a", "b", "c", "d"), nrow = 2)))
})

test_that("is_count_like ignores NAs by default and rejects them when allow_na = FALSE", {
  mat <- matrix(c(1, NA, 3, 4), nrow = 2)
  expect_true(is_count_like(mat))
  expect_false(is_count_like(mat, allow_na = FALSE))
})

test_that("is_count_like works on a phyloseq otu_table", {
  ot <- phyloseq::otu_table(
    matrix(1:6, nrow = 2, dimnames = list(c("t1", "t2"), c("s1", "s2", "s3"))),
    taxa_are_rows = TRUE
  )
  expect_true(is_count_like(ot))
})

test_that("is_count_like checks only stored (non-zero) values of a sparse_otu_table", {
  mat <- matrix(
    c(1L, 0L, 0L, 3L),
    nrow = 2,
    dimnames = list(c("t1", "t2"), c("s1", "s2"))
  )
  sp <- sparse_otu_table(phyloseq::otu_table(mat, taxa_are_rows = TRUE))
  expect_true(is_count_like(sp))
})

test_that("is_count_like's consider_small_part subsamples large matrices; consider_small_part = FALSE checks everything", {
  mat <- matrix(1L, nrow = 150, ncol = 150)
  mat[150, 150] <- 0.5 # violation outside the checked 100x100 corner
  expect_true(is_count_like(mat, consider_small_part = TRUE))
  expect_false(is_count_like(mat, consider_small_part = FALSE))
})

# ---- make_unique_taxa_table ----

test_that("make_unique_taxa_table leaves 'good' duplication (identical full lineage) untouched", {
  taxa_table <- data.frame(
    Genus = c("Bacteroides", "Bacteroides"),
    Species = c("caccae", "caccae")
  )
  result <- make_unique_taxa_table(taxa_table)
  expect_equal(result$Species, c("caccae", "caccae"))
})

test_that("make_unique_taxa_table disambiguates 'bad' duplication (same name, different lineage)", {
  taxa_table <- data.frame(
    Genus = c("Bacteroides", "Anaerostipes"),
    Species = c("caccae", "caccae")
  )
  result <- make_unique_taxa_table(taxa_table)
  expect_equal(result$Species, c("caccae", "caccae.1"))
})

test_that("make_unique_taxa_table replaces NA entries with the literal string 'NA'", {
  taxa_table <- data.frame(
    Genus = c("Bacteroides", NA),
    Species = c("caccae", "unknown")
  )
  result <- make_unique_taxa_table(taxa_table)
  expect_equal(result$Genus, c("Bacteroides", "NA"))
})

test_that("make_unique_taxa_table cascades correctly across three rank columns", {
  taxa_table <- data.frame(
    Phylum = c("P1", "P1", "P2"),
    Genus = c("Bacteroides", "Bacteroides", "Anaerostipes"),
    Species = c("caccae", "caccae", "caccae")
  )
  result <- make_unique_taxa_table(taxa_table)
  expect_equal(result$Genus, c("Bacteroides", "Bacteroides", "Anaerostipes"))
  expect_equal(result$Species, c("caccae", "caccae", "caccae.1"))
})

# ---- reorder_taxonomy_columns ----

test_that("reorder_taxonomy_columns leaves a well-formed table unchanged", {
  taxa_table <- data.frame(
    Kingdom = c("Bacteria", "Bacteria", "Bacteria", "Bacteria"),
    Genus = c("Bacteroides", "Bacteroides", "Anaerostipes", "Anaerostipes"),
    Species = c("caccae", "vulgatus", "hadrus", "caccae")
  )
  result <- reorder_taxonomy_columns(taxa_table)
  expect_equal(result$taxa_table, taxa_table)
  expect_equal(result$moved, character(0))
})

test_that("reorder_taxonomy_columns relocates a leading fully-unique ASV id column", {
  taxa_table <- data.frame(
    ASV = paste0("asv", 1:4),
    Genus = c("Bacteroides", "Bacteroides", "Anaerostipes", "Anaerostipes"),
    Species = c("caccae", "vulgatus", "hadrus", "caccae"),
    stringsAsFactors = FALSE
  )
  result <- reorder_taxonomy_columns(taxa_table)
  expect_equal(colnames(result$taxa_table), c("Genus", "Species", "ASV"))
  expect_equal(result$moved, "ASV")
  expect_equal(result$taxa_table$Genus, taxa_table$Genus)
})

test_that("reorder_taxonomy_columns relocates two leading id columns of different cardinality", {
  taxa_table <- data.frame(
    ASV = paste0("asv", 1:6),
    OTU_id = c("otu1", "otu2", "otu3", "otu4", "otu5", "otu5"),
    Genus = c(
      "Bacteroides",
      "Bacteroides",
      "Bacteroides",
      "Anaerostipes",
      "Anaerostipes",
      "Anaerostipes"
    ),
    Species = c("caccae", "caccae", "vulgatus", "hadrus", "hadrus", "caccae"),
    stringsAsFactors = FALSE
  )
  result <- reorder_taxonomy_columns(taxa_table)
  expect_equal(colnames(result$taxa_table), c("Genus", "Species", "ASV", "OTU_id"))
  expect_equal(result$moved, c("ASV", "OTU_id"))
})

test_that("reorder_taxonomy_columns still relocates a leading id column tied with a later column", {
  # ASV (4 uniques) happens to tie with Species (also fully unique) further right, but
  # the anomaly is local (ASV vs its immediate neighbour Genus), so it's still flagged.
  taxa_table <- data.frame(
    ASV = paste0("asv", 1:4),
    Genus = c("Bacteroides", "Bacteroides", "Anaerostipes", "Anaerostipes"),
    Species = paste0("sp", 1:4),
    stringsAsFactors = FALSE
  )
  result <- reorder_taxonomy_columns(taxa_table)
  expect_equal(colnames(result$taxa_table), c("Genus", "Species", "ASV"))
  expect_equal(result$moved, "ASV")
})

test_that("reorder_taxonomy_columns does not flag exactly tied adjacent leading columns (known limitation)", {
  # ASV and OTU_id both have 4 (fully unique) distinct values -- an exact tie at the
  # very first comparison, which the local/adjacent check can't disambiguate from a
  # legitimate rank column, so it conservatively leaves the table untouched.
  taxa_table <- data.frame(
    ASV = paste0("asv", 1:4),
    OTU_id = paste0("otu", 1:4),
    Genus = c("Bacteroides", "Bacteroides", "Anaerostipes", "Anaerostipes"),
    Species = c("caccae", "caccae", "hadrus", "hadrus"),
    stringsAsFactors = FALSE
  )
  result <- reorder_taxonomy_columns(taxa_table)
  expect_equal(result$taxa_table, taxa_table)
  expect_equal(result$moved, character(0))
})

test_that("reorder_taxonomy_columns ignores a dip occurring later in an otherwise normal hierarchy", {
  # Class (8) dips below Phylum (15) -- a plausible real-world sparse-classification
  # artifact, not a leading id column -- so nothing before it should be touched.
  taxa_table <- data.frame(
    Kingdom = rep("Bacteria", 8),
    Phylum = paste0("phy", c(1, 2, 3, 4, 5, 6, 7, 8)),
    Class = paste0("cls", c(1, 1, 2, 2, 3, 3, 4, 4)),
    Species = paste0("sp", 1:8),
    stringsAsFactors = FALSE
  )
  result <- reorder_taxonomy_columns(taxa_table)
  expect_equal(result$taxa_table, taxa_table)
  expect_equal(result$moved, character(0))
})

test_that("reorder_taxonomy_columns detects a leading id column even with a highly unique final rank", {
  # Regression test: a legitimately high-cardinality trailing rank (e.g. an implicit
  # ASV/species-name rank) must not mask an anomaly anchored at the front.
  taxa_table <- data.frame(
    ASV_id = paste0("id", 1:20),
    Phylum = rep(c("Bacteroidota", "Firmicutes"), 10),
    Class = paste0("cls", rep(1:4, 5)),
    Species = paste0("sp", 1:20),
    stringsAsFactors = FALSE
  )
  result <- reorder_taxonomy_columns(taxa_table)
  expect_equal(colnames(result$taxa_table), c("Phylum", "Class", "Species", "ASV_id"))
  expect_equal(result$moved, "ASV_id")
})

test_that("reorder_taxonomy_columns leaves single-column tables unchanged", {
  taxa_table <- data.frame(Genus = c("Bacteroides", "Anaerostipes"))
  result <- reorder_taxonomy_columns(taxa_table)
  expect_equal(result$taxa_table, taxa_table)
  expect_equal(result$moved, character(0))
})

# ---- impute_with_central_tendency ----

test_that("impute_with_central_tendency fills numeric NAs with the column mean/median", {
  df <- data.frame(a = c(1, NA, 3, 4))
  expect_equal(impute_with_central_tendency(df, "mean")$a, c(1, 8 / 3, 3, 4))
  expect_equal(impute_with_central_tendency(df, "median")$a, c(1, 3, 3, 4))
})

test_that("impute_with_central_tendency mode imputation uses the most frequent value", {
  df <- data.frame(a = c(1, 1, 2, NA))
  result <- impute_with_central_tendency(df, "mode")
  expect_equal(result$a, c(1, 1, 2, 1))
})

test_that("impute_with_central_tendency imputes factor columns with their mode, ignoring character columns", {
  df <- data.frame(
    fct = factor(c("x", "x", NA, "y")),
    chr = c("a", NA, "c", "d"),
    stringsAsFactors = FALSE
  )
  result <- impute_with_central_tendency(df)
  expect_equal(as.character(result$fct), c("x", "x", "x", "y"))
  expect_true(is.na(result$chr[2]))
})

test_that("impute_with_central_tendency errors clearly on an invalid central_tendency", {
  df <- data.frame(a = c(1, NA, 3))
  expect_error(impute_with_central_tendency(df, "bogus"), "should be one of")
})

# ---- dataImpute ----

test_that("dataImpute 'Replace NA with 0' fills NAs with zero and passes exceptions through unchanged", {
  df <- data.frame(id = 1:3, a = c(1, NA, 3), b = c(NA, 5, 6))
  result <- dataImpute(df, exceptions = "id", method = "Replace NA with 0")
  expect_equal(colnames(result), c("id", "a", "b"))
  expect_equal(result$id, 1:3)
  expect_equal(result$a, c(1, 0, 3))
  expect_equal(result$b, c(0, 5, 6))
})

test_that("dataImpute 'Central Tendency' delegates to impute_with_central_tendency", {
  df <- data.frame(a = c(1, NA, 3, 4))
  result <- dataImpute(df, method = "Central Tendency", central_tendency = "median")
  expect_equal(result$a, c(1, 3, 3, 4))
})

test_that("dataImpute 'KNN' fills NAs using neighboring rows", {
  set.seed(1)
  df <- data.frame(
    a = c(1, 2, 1, 2, NA, 1, 2, 1, 2, 1),
    b = c(10, 20, 10, 20, 10, 10, 20, 10, 20, 10)
  )
  result <- dataImpute(df, method = "KNN", nb_neighbors = 3)
  expect_false(anyNA(result$a))
})

test_that("dataImpute errors clearly (not VIM's cryptic error) when method = 'KNN' hits an all-NA column", {
  df <- data.frame(a = c(1, 2, 3), all_na_col = c(NA_real_, NA_real_, NA_real_))
  expect_error(dataImpute(df, method = "KNN"), "all_na_col")
})

test_that("dataImpute's all-NA KNN guard ignores exceptions columns", {
  df <- data.frame(
    id = rep(NA_real_, 5),
    a = c(1, 2, NA, 4, 5),
    b = c(10, 9, 8, 7, 6)
  )
  expect_error(
    dataImpute(df, exceptions = "id", method = "KNN", nb_neighbors = 2),
    NA
  )
})

test_that("dataImpute errors clearly on an invalid method", {
  df <- data.frame(a = c(1, NA, 3))
  expect_error(dataImpute(df, method = "bogus"), "should be one of")
})

# ---- plot_phylo_tree ----

make_phylo_tree_ps <- function() {
  n_taxa <- 6
  taxa <- paste0("ASV", seq_len(n_taxa))
  set.seed(1)
  otu <- matrix(
    sample(1:20, n_taxa * 4, replace = TRUE),
    nrow = n_taxa,
    dimnames = list(taxa, paste0("sample", 1:4))
  )
  tax <- matrix(
    c(
      rep("Bacteria", n_taxa),
      rep(c("Firmicutes", "Bacteroidetes"), each = n_taxa / 2)
    ),
    nrow = n_taxa,
    dimnames = list(taxa, c("Kingdom", "Phylum"))
  )
  sam <- data.frame(
    fraction = rep(c("pos", "neg"), 2),
    row.names = paste0("sample", 1:4)
  )
  tree <- ape::rtree(n_taxa, tip.label = taxa)
  phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = TRUE),
    phyloseq::tax_table(tax),
    phyloseq::sample_data(sam),
    tree
  )
}

test_that("plot_phylo_tree errors on a non-phyloseq input", {
  expect_error(plot_phylo_tree(NULL), "Need a phyloseq object")
  expect_error(plot_phylo_tree(list(a = 1)), "Need a phyloseq object")
})

test_that("plot_phylo_tree errors when physeq has no phy_tree", {
  ps <- make_phylo_tree_ps()
  ps_no_tree <- phyloseq::phyloseq(
    phyloseq::otu_table(ps),
    phyloseq::tax_table(ps),
    phyloseq::sample_data(ps)
  )
  expect_error(plot_phylo_tree(ps_no_tree), "has to contain a tree")
})

test_that("plot_phylo_tree errors on an invalid ladderize value", {
  ps <- make_phylo_tree_ps()
  expect_error(
    suppressWarnings(plot_phylo_tree(ps, ladderize = "up")),
    '"left", or "right"'
  )
})

test_that("plot_phylo_tree returns a ggtree object for a plain default call", {
  ps <- make_phylo_tree_ps()
  p <- suppressWarnings(plot_phylo_tree(ps))
  expect_s3_class(p, "ggtree")
})

test_that("plot_phylo_tree accepts ladderize = NULL (no reordering)", {
  ps <- make_phylo_tree_ps()
  p <- suppressWarnings(plot_phylo_tree(ps, ladderize = NULL))
  expect_s3_class(p, "ggtree")
})

test_that("plot_phylo_tree's ladderize = \"left\"/\"right\" produce different tip orders", {
  ps <- make_phylo_tree_ps()
  p_left <- suppressWarnings(plot_phylo_tree(ps, ladderize = "left"))
  p_right <- suppressWarnings(plot_phylo_tree(ps, ladderize = "right"))
  order_left <- p_left$data$label[order(p_left$data$y)]
  order_right <- p_right$data$label[order(p_right$data$y)]
  expect_false(identical(order_left, order_right))
})

test_that("plot_phylo_tree agglomerates to the requested taxrank", {
  ps <- make_phylo_tree_ps()
  p <- suppressWarnings(plot_phylo_tree(ps, taxrank = "Phylum"))
  expect_equal(sum(p$data$isTip), 2)
})

test_that("plot_phylo_tree(label_tips = TRUE) truncates long tip names and adds labels", {
  ps <- make_phylo_tree_ps()
  # Distinguishing digit up front so truncation to 25 chars doesn't collapse
  # every tip to the same prefix (which would force make.unique() to pad the
  # truncated names back out past 25 chars, defeating the point of this test).
  phyloseq::taxa_names(ps) <- paste0(
    seq_len(phyloseq::ntaxa(ps)),
    strrep("x", 30)
  )
  phyloseq::phy_tree(ps)$tip.label <- phyloseq::taxa_names(ps)
  p <- suppressWarnings(plot_phylo_tree(ps, label_tips = TRUE))
  tip_labels <- p$data$label[p$data$isTip]
  expect_true(all(nchar(tip_labels) <= 25))
  has_tiplab_layer <- any(vapply(
    p$layers,
    function(l) inherits(l$geom, "GeomText") || inherits(l$geom, "GeomLabel"),
    logical(1)
  ))
  expect_true(has_tiplab_layer)
})

test_that("plot_phylo_tree(label_tips = FALSE) doesn't truncate tip names", {
  ps <- make_phylo_tree_ps()
  long_names <- paste0(strrep("x", 30), seq_len(phyloseq::ntaxa(ps)))
  phyloseq::taxa_names(ps) <- long_names
  phyloseq::phy_tree(ps)$tip.label <- long_names
  p <- suppressWarnings(plot_phylo_tree(ps, label_tips = FALSE))
  tip_labels <- p$data$label[p$data$isTip]
  expect_true(any(nchar(tip_labels) > 25))
})

test_that("plot_phylo_tree colors tips by a valid tax_table column", {
  ps <- make_phylo_tree_ps()
  p <- suppressWarnings(plot_phylo_tree(ps, tip_color = "Phylum"))
  expect_true("Phylum" %in% colnames(p$data))
})

test_that("plot_phylo_tree silently ignores an unknown tip_color column", {
  ps <- make_phylo_tree_ps()
  p <- suppressWarnings(plot_phylo_tree(ps, tip_color = "not_a_column"))
  expect_s3_class(p, "ggtree")
  expect_false("not_a_column" %in% colnames(p$data))
})

test_that("plot_phylo_tree restricts samples via fraction_id_name/fraction_ids without error", {
  ps <- make_phylo_tree_ps()
  p <- suppressWarnings(
    plot_phylo_tree(
      ps,
      fraction_id_name = "fraction",
      fraction_ids = "pos"
    )
  )
  expect_s3_class(p, "ggtree")
})

test_that("plot_phylo_tree forwards ... to ggtree::ggtree()", {
  ps <- make_phylo_tree_ps()
  p <- suppressWarnings(plot_phylo_tree(ps, branch.length = "none"))
  expect_s3_class(p, "ggtree")
})

.quiet_ggtree_cosmetics <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("must be used|deprecated", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

test_that("plot_phylo_tree midpoint-roots and warns on an unrooted tree", {
  ps <- make_phylo_tree_ps()
  phyloseq::phy_tree(ps) <- ape::unroot(phyloseq::phy_tree(ps))
  expect_false(ape::is.rooted(phyloseq::phy_tree(ps)))

  expect_warning(
    p <- .quiet_ggtree_cosmetics(plot_phylo_tree(ps)),
    "unrooted, midpoint was set as root"
  )
  expect_s3_class(p, "ggtree")
  expect_setequal(
    p$data$label[p$data$isTip],
    phyloseq::taxa_names(ps)
  )
})

test_that("plot_phylo_tree doesn't warn when the tree is already rooted", {
  ps <- make_phylo_tree_ps()
  expect_true(ape::is.rooted(phyloseq::phy_tree(ps)))
  expect_no_warning(.quiet_ggtree_cosmetics(plot_phylo_tree(ps)))
})
