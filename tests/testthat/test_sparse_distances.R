library(testthat)
library(phyloseq)
library(Matrix)
library(vegan)

source("../R/sparse_otu_table.R")
source("../R/sparse_distances.R")

# ---- Helpers ----

make_pair <- function(
  n_taxa = 100,
  n_samples = 20,
  sparsity = 0.95,
  taxa_are_rows = TRUE,
  seed = 42
) {
  set.seed(seed)
  mat <- matrix(0L, n_taxa, n_samples)
  n_nz <- max(2L, round((1 - sparsity) * n_taxa * n_samples))
  mat[sample(n_taxa * n_samples, n_nz)] <- sample(
    1L:10000L,
    n_nz,
    replace = TRUE
  )
  rownames(mat) <- paste0("ASV", seq_len(n_taxa))
  colnames(mat) <- paste0("S", seq_len(n_samples))

  if (!taxa_are_rows) {
    mat <- t(mat)
  }

  snames <- if (taxa_are_rows) colnames(mat) else rownames(mat)
  sdata <- data.frame(id = seq_len(n_samples), row.names = snames)
  ps_dense <- phyloseq(
    otu_table(mat, taxa_are_rows = taxa_are_rows),
    sample_data(sdata)
  )
  ps_sparse <- as_sparse_phyloseq(ps_dense)
  list(dense = ps_dense, sparse = ps_sparse)
}

# Reference: vegan Bray-Curtis (samples × taxa)
vegan_bray <- function(ps) {
  m <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(otu_table(ps))) {
    m <- t(m)
  }
  vegan::vegdist(m, method = "bray")
}

pair <- make_pair()
pair_t <- make_pair(taxa_are_rows = FALSE)

# ---- 1. Return type ----

test_that("bray_curtis_sparse returns a dist object", {
  d <- bray_curtis_sparse(pair$sparse)
  expect_s3_class(d, "dist")
})

test_that("dist has correct size (n*(n-1)/2 elements)", {
  n <- nsamples(pair$sparse)
  d <- bray_curtis_sparse(pair$sparse)
  expect_equal(length(d), n * (n - 1L) / 2L)
})

test_that("dist labels match sample names", {
  d <- bray_curtis_sparse(pair$sparse)
  expect_identical(attr(d, "Labels"), sample_names(pair$sparse))
})

# ---- 2. Values are in [0, 1] ----

test_that("all distances are in [0, 1]", {
  d <- bray_curtis_sparse(pair$sparse)
  expect_true(all(as.numeric(d) >= 0))
  expect_true(all(as.numeric(d) <= 1))
})

# ---- 3. Correctness vs vegan ----

test_that("matches vegan::vegdist (tar = TRUE)", {
  dm_sparse <- bray_curtis_sparse(pair$sparse)
  dm_vegan <- vegan_bray(pair$dense)
  expect_equal(as.numeric(dm_sparse), as.numeric(dm_vegan), tolerance = 1e-8)
})

test_that("matches vegan::vegdist (tar = FALSE)", {
  dm_sparse <- bray_curtis_sparse(pair_t$sparse)
  dm_vegan <- vegan_bray(pair_t$dense)
  expect_equal(as.numeric(dm_sparse), as.numeric(dm_vegan), tolerance = 1e-8)
})

test_that("matches phyloseq::distance (dense reference)", {
  dm_sparse <- bray_curtis_sparse(pair$sparse)
  dm_phyloseq <- phyloseq::distance(pair$dense, method = "bray")
  expect_equal(as.numeric(dm_sparse), as.numeric(dm_phyloseq), tolerance = 1e-8)
})

# ---- 4. Dense phyloseq fallback ----

test_that("accepts plain dense phyloseq (no sparse_otu_table)", {
  dm_dense <- bray_curtis_sparse(pair$dense)
  dm_vegan <- vegan_bray(pair$dense)
  expect_equal(as.numeric(dm_dense), as.numeric(dm_vegan), tolerance = 1e-8)
})

# ---- 5. Hand-computed two-sample case ----

test_that("two-sample Bray-Curtis is correct by formula", {
  # s1 = [10, 0, 5], s2 = [0, 3, 5]
  # C = min(10,0) + min(0,3) + min(5,5) = 5
  # BC = 1 - 2*5 / (15 + 8) = 13/23
  mat <- matrix(c(10L, 0L, 5L, 0L, 3L, 5L), nrow = 3, ncol = 2)
  rownames(mat) <- paste0("t", 1:3)
  colnames(mat) <- c("s1", "s2")
  ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
  otu_table(ps) <- sparse_otu_table(otu_table(ps))

  d <- as.numeric(bray_curtis_sparse(ps))
  expect_equal(d, 13 / 23, tolerance = 1e-12)
})

# ---- 6. Identical samples → distance 0 ----

test_that("identical samples have distance 0", {
  mat <- matrix(c(1L, 2L, 3L, 1L, 2L, 3L), nrow = 3)
  rownames(mat) <- paste0("t", 1:3)
  colnames(mat) <- c("s1", "s2")
  ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
  otu_table(ps) <- sparse_otu_table(otu_table(ps))

  expect_equal(as.numeric(bray_curtis_sparse(ps)), 0, tolerance = 1e-12)
})

# ---- 7. Completely disjoint samples → distance 1 ----

test_that("completely disjoint samples have distance 1", {
  # 2×2 diagonal matrices become ddiMatrix in Matrix; use dense ps to avoid it
  mat <- matrix(c(5L, 0L, 0L, 3L), nrow = 2)
  rownames(mat) <- c("t1", "t2")
  colnames(mat) <- c("s1", "s2")
  ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))

  expect_equal(as.numeric(bray_curtis_sparse(ps)), 1, tolerance = 1e-12)
})

# ---- 8. Single-sample phyloseq ----

test_that("single-sample phyloseq returns empty dist of length 0", {
  mat <- matrix(c(1L, 2L, 3L), nrow = 3)
  rownames(mat) <- paste0("t", 1:3)
  colnames(mat) <- "s1"
  # use dense ps: sparse_otu_table may produce ddiMatrix for a 1-column matrix
  ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))

  d <- bray_curtis_sparse(ps)
  expect_s3_class(d, "dist")
  expect_equal(length(d), 0L)
})

# ---- 9. Symmetry ----

test_that("distance matrix is symmetric", {
  d <- bray_curtis_sparse(pair$sparse)
  m <- as.matrix(d)
  expect_equal(m, t(m))
})

# ---- 10. Diagonal is zero ----

test_that("diagonal of full distance matrix is zero", {
  d <- bray_curtis_sparse(pair$sparse)
  expect_true(all(diag(as.matrix(d)) == 0))
})

# ---- 11. Low-sparsity (near-dense) matrix ----

test_that("works correctly at low sparsity (50%)", {
  p <- make_pair(n_taxa = 50, n_samples = 10, sparsity = 0.50, seed = 7)
  dm_sparse <- bray_curtis_sparse(p$sparse)
  dm_vegan <- vegan_bray(p$dense)
  expect_equal(as.numeric(dm_sparse), as.numeric(dm_vegan), tolerance = 1e-8)
})
