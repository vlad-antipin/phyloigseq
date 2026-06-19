library(phyloseq)
library(PhyloIgSeq)

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

pair   <- make_pair()
pair_t <- make_pair(taxa_are_rows = FALSE)

# ---- SPARSE_DISTANCE_METHODS metadata ----

test_that("SPARSE_DISTANCE_METHODS contains 'bray'", {
  expect_true("bray" %in% SPARSE_DISTANCE_METHODS)
})

test_that("unsupported method falls back to phyloseq::distance with a warning", {
  expect_warning(
    sparse_distance(pair$dense, method = "jaccard"),
    regexp = "No sparse version"
  )
})

# ---- Per-method generic tests ----
# Runs for every method in SPARSE_DISTANCE_METHODS.
# local() captures the loop variable so each block closes over its own `m`.

for (method in SPARSE_DISTANCE_METHODS) {
  local({
    m <- method

    test_that(paste0("[", m, "] returns a dist object"), {
      expect_s3_class(sparse_distance(pair$sparse, method = m), "dist")
    })

    test_that(paste0("[", m, "] dist has correct size (n*(n-1)/2)"), {
      n <- nsamples(pair$sparse)
      expect_equal(length(sparse_distance(pair$sparse, method = m)), n * (n - 1L) / 2L)
    })

    test_that(paste0("[", m, "] dist labels match sample names"), {
      expect_identical(
        attr(sparse_distance(pair$sparse, method = m), "Labels"),
        sample_names(pair$sparse)
      )
    })

    test_that(paste0("[", m, "] matches phyloseq::distance (tar = TRUE)"), {
      expect_equal(
        as.numeric(sparse_distance(pair$sparse, method = m)),
        as.numeric(phyloseq::distance(pair$dense, method = m)),
        tolerance = 1e-8
      )
    })

    test_that(paste0("[", m, "] matches phyloseq::distance (tar = FALSE)"), {
      expect_equal(
        as.numeric(sparse_distance(pair_t$sparse, method = m)),
        as.numeric(phyloseq::distance(pair_t$dense, method = m)),
        tolerance = 1e-8
      )
    })

    test_that(paste0("[", m, "] accepts plain dense phyloseq"), {
      expect_equal(
        as.numeric(sparse_distance(pair$dense, method = m)),
        as.numeric(phyloseq::distance(pair$dense, method = m)),
        tolerance = 1e-8
      )
    })

    test_that(paste0("[", m, "] identical samples have distance 0"), {
      mat <- matrix(c(1L, 2L, 3L, 1L, 2L, 3L), nrow = 3)
      rownames(mat) <- paste0("t", 1:3)
      colnames(mat) <- c("s1", "s2")
      ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
      otu_table(ps) <- sparse_otu_table(otu_table(ps))
      expect_equal(as.numeric(sparse_distance(ps, method = m)), 0, tolerance = 1e-12)
    })

    test_that(paste0("[", m, "] single-sample returns empty dist"), {
      mat <- matrix(c(1L, 2L, 3L), nrow = 3)
      rownames(mat) <- paste0("t", 1:3)
      colnames(mat) <- "s1"
      ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
      d <- sparse_distance(ps, method = m)
      expect_s3_class(d, "dist")
      expect_equal(length(d), 0L)
    })

    test_that(paste0("[", m, "] distance matrix is symmetric"), {
      dm <- as.matrix(sparse_distance(pair$sparse, method = m))
      expect_equal(dm, t(dm))
    })

    test_that(paste0("[", m, "] diagonal of distance matrix is zero"), {
      dm <- as.matrix(sparse_distance(pair$sparse, method = m))
      expect_true(all(diag(dm) == 0))
    })

    test_that(paste0("[", m, "] works correctly at low sparsity (50%)"), {
      p <- make_pair(n_taxa = 50, n_samples = 10, sparsity = 0.50, seed = 7)
      expect_equal(
        as.numeric(sparse_distance(p$sparse, method = m)),
        as.numeric(phyloseq::distance(p$dense, method = m)),
        tolerance = 1e-8
      )
    })
  })
}

# ---- Bray-Curtis specific ----

test_that("[bray] all distances are in [0, 1]", {
  d <- sparse_distance(pair$sparse, method = "bray")
  expect_true(all(as.numeric(d) >= 0))
  expect_true(all(as.numeric(d) <= 1))
})

test_that("[bray] two-sample hand-computed formula (13/23)", {
  mat <- matrix(c(10L, 0L, 5L, 0L, 3L, 5L), nrow = 3, ncol = 2)
  rownames(mat) <- paste0("t", 1:3)
  colnames(mat) <- c("s1", "s2")
  ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
  otu_table(ps) <- sparse_otu_table(otu_table(ps))
  expect_equal(as.numeric(sparse_distance(ps, method = "bray")), 13 / 23, tolerance = 1e-12)
})

test_that("[bray] completely disjoint samples have distance 1", {
  mat <- matrix(c(5L, 0L, 0L, 3L), nrow = 2)
  rownames(mat) <- c("t1", "t2")
  colnames(mat) <- c("s1", "s2")
  ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
  expect_equal(as.numeric(sparse_distance(ps, method = "bray")), 1, tolerance = 1e-12)
})
