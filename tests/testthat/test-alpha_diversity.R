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

# ---- get_alpha_diversity: sparse path ----

test_that("get_alpha_diversity uses sparse path for sparse input", {
  result <- get_alpha_diversity(pair$sparse, measure = "Shannon")
  expect_true("Shannon" %in% colnames(result))
})

test_that("get_alpha_diversity sparse and dense produce identical Shannon values", {
  res_dense <- suppressWarnings(get_alpha_diversity(
    pair$dense,
    measure = "Shannon"
  ))
  res_sparse <- get_alpha_diversity(pair$sparse, measure = "Shannon")
  expect_equal(res_sparse$Shannon, res_dense$Shannon, tolerance = 1e-10)
})

test_that("get_alpha_diversity sparse path works for taxa_are_rows = FALSE", {
  res_dense <- suppressWarnings(get_alpha_diversity(
    pair_t$dense,
    measure = "Shannon"
  ))
  res_sparse <- get_alpha_diversity(pair_t$sparse, measure = "Shannon")
  expect_equal(res_sparse$Shannon, res_dense$Shannon, tolerance = 1e-10)
})

test_that("get_alpha_diversity result has correct row names and sample_data columns", {
  res <- get_alpha_diversity(pair$sparse, measure = "Shannon")
  expect_equal(rownames(res), sample_names(pair$sparse))
  expect_true("Group" %in% colnames(res))
  expect_true("depth" %in% colnames(res))
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
  expect_equal(res_sparse$Observed, res_dense$Observed, tolerance = 1e-10)
})
