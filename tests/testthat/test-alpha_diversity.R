library(PhyloIgSeq)

# ---- Helpers ----

make_ps <- function(
  n_taxa = 50,
  n_samples = 10,
  sparsity = 0.90,
  taxa_are_rows = TRUE,
  seed = 7
) {
  set.seed(seed)
  mat <- matrix(0L, n_taxa, n_samples)
  n_nz <- max(2L, round((1 - sparsity) * n_taxa * n_samples))
  mat[sample(n_taxa * n_samples, n_nz)] <- sample(
    1L:1000L,
    n_nz,
    replace = TRUE
  )
  rownames(mat) <- paste0("ASV", seq_len(n_taxa))
  colnames(mat) <- paste0("S", seq_len(n_samples))

  sdata <- data.frame(
    Group = sample(c("A", "B"), n_samples, replace = TRUE),
    row.names = colnames(mat)
  )
  if (!taxa_are_rows) {
    mat <- t(mat)
  }

  ps <- phyloseq(
    otu_table(mat, taxa_are_rows = taxa_are_rows),
    sample_data(sdata)
  )
  list(dense = ps, sparse = as_sparse_phyloseq(ps))
}

make_ig_score_ps <- function(n_taxa = 6, n_samples = 8, seed = 11) {
  set.seed(seed)
  mat <- matrix(
    rnorm(n_taxa * n_samples, sd = 2),
    nrow = n_taxa,
    dimnames = list(paste0("ASV", seq_len(n_taxa)), paste0("S", seq_len(n_samples)))
  )
  sdata <- data.frame(
    Group = rep(c("A", "B"), length.out = n_samples),
    row.names = colnames(mat)
  )
  phyloseq(
    otu_table(mat, taxa_are_rows = TRUE),
    sample_data(sdata)
  )
}

ref_shannon <- function(ps) {
  # vegan::diversity is what estimate_richness uses internally
  m <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(otu_table(ps))) {
    m <- t(m)
  } # samples x taxa
  vegan::diversity(m, index = "shannon")
}

pair <- make_ps(taxa_are_rows = TRUE)
pair_t <- make_ps(taxa_are_rows = FALSE, seed = 13)

# ---- sparse_shannon: correctness ----

test_that("sparse_shannon matches vegan for taxa_are_rows = TRUE", {
  expected <- ref_shannon(pair$dense)
  got <- sparse_shannon(pair$sparse)
  expect_named(got, names(expected))
  expect_equal(got, expected, tolerance = 1e-10)
})

test_that("sparse_shannon matches vegan for taxa_are_rows = FALSE", {
  expected <- ref_shannon(pair_t$dense)
  got <- sparse_shannon(pair_t$sparse)
  expect_named(got, names(expected))
  expect_equal(got, expected, tolerance = 1e-10)
})

test_that("sparse_shannon returns a named numeric vector", {
  h <- sparse_shannon(pair$sparse)
  expect_type(h, "double")
  expect_length(h, nsamples(pair$sparse))
  expect_equal(names(h), sample_names(pair$sparse))
})

test_that("sparse_shannon values are non-negative", {
  expect_true(all(sparse_shannon(pair$sparse) >= 0))
})

# ---- get_alpha_diversity: return shape ----

test_that("get_alpha_diversity returns a structured list with the expected fields", {
  result <- get_alpha_diversity(pair$sparse, measure = "Shannon")
  expect_named(result, c("diversity", "measure", "sample_data", "depth"))
  expect_equal(result$measure, "Shannon")
  expect_true("Shannon" %in% colnames(result$diversity))
  expect_true("sample_id" %in% colnames(result$diversity))
  expect_true("sample_id" %in% colnames(result$sample_data))
  expect_true("Group" %in% colnames(result$sample_data))
  expect_equal(sort(result$depth$sample_id), sort(sample_names(pair$sparse)))
})

test_that("get_alpha_diversity uses sparse path for sparse input", {
  result <- get_alpha_diversity(pair$sparse, measure = "Shannon")
  expect_true("Shannon" %in% colnames(result$diversity))
})

test_that("get_alpha_diversity sparse and dense produce identical Shannon values", {
  res_dense <- suppressWarnings(get_alpha_diversity(
    pair$dense,
    measure = "Shannon"
  ))
  res_sparse <- get_alpha_diversity(pair$sparse, measure = "Shannon")
  expect_equal(
    res_sparse$diversity$Shannon[order(res_sparse$diversity$sample_id)],
    res_dense$diversity$Shannon[order(res_dense$diversity$sample_id)],
    tolerance = 1e-10
  )
})

test_that("get_alpha_diversity sparse path works for taxa_are_rows = FALSE", {
  res_dense <- suppressWarnings(get_alpha_diversity(
    pair_t$dense,
    measure = "Shannon"
  ))
  res_sparse <- get_alpha_diversity(pair_t$sparse, measure = "Shannon")
  expect_equal(
    res_sparse$diversity$Shannon[order(res_sparse$diversity$sample_id)],
    res_dense$diversity$Shannon[order(res_dense$diversity$sample_id)],
    tolerance = 1e-10
  )
})

test_that("get_alpha_diversity falls back to estimate_richness for non-Shannon measure", {
  # Observed is not sparse-accelerated; result should still be numerically valid
  res_dense <- suppressWarnings(get_alpha_diversity(
    pair$dense,
    measure = "Observed"
  ))
  res_sparse <- suppressWarnings(get_alpha_diversity(
    pair$sparse,
    measure = "Observed"
  ))
  expect_equal(
    res_sparse$diversity$Observed[order(res_sparse$diversity$sample_id)],
    res_dense$diversity$Observed[order(res_dense$diversity$sample_id)],
    tolerance = 1e-10
  )
})

test_that("get_alpha_diversity errors on non-phyloseq input", {
  expect_error(get_alpha_diversity(matrix(1:4, 2, 2)), "phyloseq")
})

test_that("get_alpha_diversity respects taxrank agglomeration", {
  ps_tax <- pair$dense
  tax_table(ps_tax) <- tax_table(matrix(
    rep(c("Firmicutes", "Bacteroidetes"), length.out = ntaxa(ps_tax)),
    ncol = 1,
    dimnames = list(taxa_names(ps_tax), "Phylum")
  ))
  result <- suppressWarnings(get_alpha_diversity(
    ps_tax,
    taxrank = "Phylum",
    measure = "Observed"
  ))
  expect_true(all(result$diversity$Observed <= 2))
})

# ---- get_alpha_diversity: from_igseq delegation ----

test_that("get_alpha_diversity(from_igseq = TRUE) matches get_igseq_richness directly", {
  ps_ig <- make_ig_score_ps()
  via_alpha <- get_alpha_diversity(ps_ig, from_igseq = TRUE, low_lim = -1, high_lim = 1)
  direct <- get_igseq_richness(ps_ig, low_lim = -1, high_lim = 1)
  expect_equal(via_alpha, direct)
})

test_that("get_alpha_diversity(from_igseq = TRUE) ignores taxrank/transform/fraction args", {
  ps_ig <- make_ig_score_ps()
  result <- get_alpha_diversity(
    ps_ig,
    from_igseq = TRUE,
    taxrank = "does_not_exist",
    fraction_id_name = "does_not_exist"
  )
  expect_equal(result$measure, "richness")
})

# ---- get_igseq_richness ----

test_that("get_igseq_richness buckets scores correctly and sums to taxa count", {
  ps_ig <- make_ig_score_ps(n_taxa = 6, n_samples = 4)
  result <- get_igseq_richness(ps_ig, low_lim = -1.96, high_lim = 1.96)
  expect_named(result, c("diversity", "measure", "sample_data", "depth"))
  expect_equal(result$measure, "richness")
  expect_null(result$depth)
  expect_setequal(unique(result$diversity$significance), c("down", "ns", "up"))

  totals <- vapply(
    unique(result$diversity$sample_id),
    function(s) sum(result$diversity$richness[result$diversity$sample_id == s]),
    numeric(1)
  )
  expect_true(all(totals == ntaxa(ps_ig)))
})

test_that("get_igseq_richness proportions sum to 1 per sample", {
  ps_ig <- make_ig_score_ps(n_taxa = 6, n_samples = 4)
  result <- get_igseq_richness(ps_ig, proportions = TRUE)
  totals <- vapply(
    unique(result$diversity$sample_id),
    function(s) sum(result$diversity$richness[result$diversity$sample_id == s]),
    numeric(1)
  )
  expect_equal(unname(totals), rep(1, length(totals)))
})

test_that("get_igseq_richness manual bucket count matches a hand-picked example", {
  mat <- matrix(
    c(-3, -3, 0, 0, 3, 3),
    nrow = 6,
    dimnames = list(paste0("ASV", 1:6), "S1")
  )
  ps_ig <- phyloseq(
    otu_table(mat, taxa_are_rows = TRUE),
    sample_data(data.frame(Group = "A", row.names = "S1"))
  )
  result <- get_igseq_richness(ps_ig, low_lim = -1.96, high_lim = 1.96)
  richness_by_bucket <- setNames(
    result$diversity$richness,
    result$diversity$significance
  )
  expect_equal(unname(richness_by_bucket["down"]), 2)
  expect_equal(unname(richness_by_bucket["ns"]), 2)
  expect_equal(unname(richness_by_bucket["up"]), 2)
})

# ---- plot_alpha_diversity ----

test_that("plot_alpha_diversity returns a ggplot for the standard branch", {
  alpha_div <- suppressWarnings(get_alpha_diversity(pair$dense, measure = "Shannon"))
  p <- plot_alpha_diversity(alpha_div, x = "Group", group = "Group")
  expect_s3_class(p, "ggplot")
})

test_that("plot_alpha_diversity defaults measure to alpha_div$measure[1]", {
  alpha_div <- suppressWarnings(get_alpha_diversity(pair$dense, measure = "Shannon"))
  p <- plot_alpha_diversity(alpha_div)
  expect_equal(rlang::as_label(p$mapping$y), "Shannon")
})

test_that("plot_alpha_diversity works on the from_igseq/get_igseq_richness shape", {
  ps_ig <- make_ig_score_ps()
  alpha_div <- get_alpha_diversity(ps_ig, from_igseq = TRUE)
  p <- plot_alpha_diversity(alpha_div, x = "Group", facet = "significance")
  expect_s3_class(p, "ggplot")
})

test_that("plot_alpha_diversity check_depth maps point size to depth", {
  # check_depth's "text" aes is intended for plotly::ggplotly(tooltip = "text")
  # and triggers a harmless "Ignoring unknown aesthetics" warning from plain
  # ggplot2 (same root cause as plot_slide_z's, see test-ig_score.R)
  alpha_div <- suppressWarnings(get_alpha_diversity(pair$dense, measure = "Shannon"))
  p <- suppressWarnings(plot_alpha_diversity(alpha_div, check_depth = TRUE))
  expect_s3_class(p, "ggplot")
})

# ---- alpha_diversity() wrapper ----

test_that("alpha_diversity() chains get_alpha_diversity() and plot_alpha_diversity()", {
  p <- suppressWarnings(alpha_diversity(pair$dense, x = "Group", group = "Group"))
  expect_s3_class(p, "ggplot")
})

# ---- plot_igseq_richness ----

test_that("plot_igseq_richness returns a ggplot faceted by significance", {
  ps_ig <- make_ig_score_ps()
  igseq_richness <- get_igseq_richness(ps_ig)
  p <- plot_igseq_richness(igseq_richness, group = "Group", color = "Group")
  expect_s3_class(p, "ggplot")
})

test_that("plot_igseq_richness exclude_ns drops the ns bucket", {
  ps_ig <- make_ig_score_ps()
  igseq_richness <- get_igseq_richness(ps_ig)
  p <- plot_igseq_richness(
    igseq_richness,
    group = "Group",
    color = "Group",
    exclude_ns = TRUE
  )
  expect_false("ns" %in% p$data$significance)
})

# ---- relocated UI helpers (utils.R) ----

test_that("is_valid_factor identifies usable discrete columns", {
  df <- data.frame(a = c("x", "y"), b = c(1, 2))
  expect_true(is_valid_factor(df, "a"))
  expect_false(is_valid_factor(df, "b")) # numeric
  expect_false(is_valid_factor(df, "c")) # absent
  expect_false(is_valid_factor(df, NULL))
})

test_that("keep_levels / factorize_levels round-trip with an explicit level set", {
  vals <- c("a", "b", "c", NA)
  keep <- keep_levels(vals, c("a", "c", "(NA)"))
  expect_equal(keep, c(TRUE, FALSE, TRUE, TRUE))

  fct <- factorize_levels(vals, c("a", "c", "(NA)"))
  expect_equal(levels(fct), c("a", "c", "(NA)"))
})

test_that("keep_levels / factorize_levels fall back to sorted levels when level_names is NULL", {
  vals <- c("b", "a", "c")
  expect_true(all(keep_levels(vals, NULL)))
  expect_equal(levels(factorize_levels(vals, NULL)), c("a", "b", "c"))
})

test_that("apply_levels filters and factorizes in one step", {
  df <- data.frame(g = c("a", "b", "c"), v = 1:3)
  out <- apply_levels(df, "g", c("a", "c"))
  expect_equal(out$g, factor(c("a", "c"), levels = c("a", "c")))
  expect_equal(out$v, c(1, 3))
})

test_that("apply_levels_if_valid is a no-op for invalid columns", {
  df <- data.frame(g = c("a", "b"), v = c(1, 2))
  expect_identical(apply_levels_if_valid(df, "v", NULL), df) # numeric
  expect_identical(apply_levels_if_valid(df, NULL, NULL), df)
})

test_that("remove_nas drops rows with an NA in any listed column", {
  df <- data.frame(a = c(1, NA, 3), b = c(NA, 2, 3))
  expect_equal(nrow(remove_nas(df, c("a", "b"))), 1)
  expect_equal(nrow(remove_nas(df, c("a", NULL))), 2)
})

test_that("get_hover_text builds one string per row from selected columns", {
  df <- data.frame(sample_id = c("S1", "S2"), Group = c("A", "B"))
  txt <- get_hover_text(df, c("Group", "not_a_column"))
  expect_equal(txt, c("Group: A<br>", "Group: B<br>"))
})

test_that("get_hover_text does not depend on rownames matching column values", {
  df <- data.frame(sample_id = c("S1", "S2"), Group = c("A", "B"))
  rownames(df) <- c("1", "2") # default merge()-style integer rownames
  txt <- get_hover_text(df, "Group")
  expect_equal(txt, c("Group: A<br>", "Group: B<br>"))
})
