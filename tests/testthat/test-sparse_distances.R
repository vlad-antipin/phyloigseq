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

# Reference distance: phyloseq::distance for methods it supports;
# vegan decostand + vegdist for chord and hellinger.
ref_distance <- function(ps, method) {
  if (method %in% c("chord", "hellinger")) {
    # vegan expects samples x taxa; respect the otu_table orientation
    m  <- as(otu_table(ps), "matrix")
    ot <- if (taxa_are_rows(otu_table(ps))) t(m) else m
    tr <- switch(method,
      chord     = vegan::decostand(ot, method = "normalize"),
      hellinger = vegan::decostand(ot, method = "hellinger")
    )
    vegan::vegdist(tr, method = "euclidean")
  } else {
    phyloseq::distance(ps, method = method)
  }
}

pair   <- make_pair()
pair_t <- make_pair(taxa_are_rows = FALSE)

# UniFrac methods need a phy_tree, which make_pair() doesn't attach. Reuse a
# pair and bolt a random tree onto its taxa (named "ASV1".."ASVn" by
# make_pair()) so sparse_distance()'s dispatch to sparse_unifrac() can be
# exercised here; full sparse_unifrac() correctness/edge-case coverage lives
# in test-sparse_unifrac.R.
make_tree_pair <- function(n_taxa = 100, ...) {
  p <- make_pair(n_taxa = n_taxa, ...)
  tree <- ape::rtree(n_taxa, tip.label = paste0("ASV", seq_len(n_taxa)))
  phy_tree(p$dense) <- tree
  phy_tree(p$sparse) <- tree
  p
}

pair_tree <- make_tree_pair()

# Distance methods that require a phy_tree and so are excluded from the
# generic per-method loop below, which runs against tree-less fixtures.
TREE_METHODS <- c("unifrac", "wunifrac")

# ---- SPARSE_DISTANCE_METHODS metadata ----

test_that("SPARSE_DISTANCE_METHODS contains 'bray'", {
  expect_true("bray" %in% SPARSE_DISTANCE_METHODS)
})

test_that("SPARSE_DISTANCE_METHODS contains all expected methods", {
  expected <- c("bray", "jaccard", "kulczynski", "manhattan",
                "euclidean", "canberra", "horn", "chord", "hellinger",
                "unifrac", "wunifrac")
  expect_true(all(expected %in% SPARSE_DISTANCE_METHODS))
})

test_that("unsupported method falls back to phyloseq::distance with a warning", {
  warns <- character(0)
  withCallingHandlers(
    sparse_distance(pair$dense, method = "gower"),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("No sparse version", warns)))
})

# ---- Per-method generic tests ----
# Runs for every method in SPARSE_DISTANCE_METHODS except TREE_METHODS (see
# the dedicated UniFrac section below).
# local() captures the loop variable so each block closes over its own `m`.

for (method in setdiff(SPARSE_DISTANCE_METHODS, TREE_METHODS)) {
  local({
    m <- method

    test_that(paste0("[", m, "] returns a dist object"), {
      expect_s3_class(sparse_distance(pair$sparse, method = m), "dist")
    })

    test_that(paste0("[", m, "] dist has correct size (n*(n-1)/2)"), {
      n <- nsamples(pair$sparse)
      expect_equal(
        length(sparse_distance(pair$sparse, method = m)),
        n * (n - 1L) / 2L
      )
    })

    test_that(paste0("[", m, "] dist labels match sample names"), {
      expect_identical(
        attr(sparse_distance(pair$sparse, method = m), "Labels"),
        sample_names(pair$sparse)
      )
    })

    test_that(paste0("[", m, "] matches reference distance (tar = TRUE)"), {
      expect_equal(
        as.numeric(sparse_distance(pair$sparse, method = m)),
        as.numeric(ref_distance(pair$dense, method = m)),
        tolerance = 1e-8
      )
    })

    test_that(paste0("[", m, "] matches reference distance (tar = FALSE)"), {
      expect_equal(
        as.numeric(sparse_distance(pair_t$sparse, method = m)),
        as.numeric(ref_distance(pair_t$dense, method = m)),
        tolerance = 1e-8
      )
    })

    test_that(paste0("[", m, "] accepts plain dense phyloseq"), {
      expect_equal(
        as.numeric(sparse_distance(pair$dense, method = m)),
        as.numeric(ref_distance(pair$dense, method = m)),
        tolerance = 1e-8
      )
    })

    test_that(paste0("[", m, "] identical samples have distance 0"), {
      mat <- matrix(c(1L, 2L, 3L, 1L, 2L, 3L), nrow = 3)
      rownames(mat) <- paste0("t", 1:3)
      colnames(mat) <- c("s1", "s2")
      ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
      otu_table(ps) <- sparse_otu_table(otu_table(ps))
      expect_equal(
        as.numeric(sparse_distance(ps, method = m)),
        0,
        tolerance = 1e-12
      )
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
        as.numeric(ref_distance(p$dense, method = m)),
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
  expect_equal(
    as.numeric(sparse_distance(ps, method = "bray")),
    13 / 23,
    tolerance = 1e-12
  )
})

test_that("[bray] completely disjoint samples have distance 1", {
  mat <- matrix(c(5L, 0L, 0L, 3L), nrow = 2)
  rownames(mat) <- c("t1", "t2")
  colnames(mat) <- c("s1", "s2")
  ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
  expect_equal(
    as.numeric(sparse_distance(ps, method = "bray")),
    1,
    tolerance = 1e-12
  )
})

# ---- Jaccard specific ----

test_that("[jaccard] all distances are in [0, 1]", {
  d <- sparse_distance(pair$sparse, method = "jaccard")
  expect_true(all(as.numeric(d) >= 0))
  expect_true(all(as.numeric(d) <= 1 + 1e-10))
})

test_that("[jaccard] two-sample hand-computed formula (13/18)", {
  # s1=(10,0,5) A=15, s2=(0,3,5) B=8, C=min(10,0)+min(0,3)+min(5,5)=5
  # jaccard = 1 - 5/(15+8-5) = 13/18
  mat <- matrix(c(10L, 0L, 5L, 0L, 3L, 5L), nrow = 3, ncol = 2)
  rownames(mat) <- paste0("t", 1:3)
  colnames(mat) <- c("s1", "s2")
  ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
  otu_table(ps) <- sparse_otu_table(otu_table(ps))
  expect_equal(
    as.numeric(sparse_distance(ps, method = "jaccard")),
    13 / 18,
    tolerance = 1e-12
  )
})

test_that("[jaccard] disjoint samples have distance 1", {
  mat <- matrix(c(5L, 0L, 0L, 3L), nrow = 2)
  rownames(mat) <- c("t1", "t2")
  colnames(mat) <- c("s1", "s2")
  ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
  expect_equal(
    as.numeric(sparse_distance(ps, method = "jaccard")),
    1,
    tolerance = 1e-12
  )
})

# ---- Kulczynski specific ----

test_that("[kulczynski] all distances are in [0, 1]", {
  d <- sparse_distance(pair$sparse, method = "kulczynski")
  expect_true(all(as.numeric(d) >= -1e-10))
  expect_true(all(as.numeric(d) <= 1 + 1e-10))
})

test_that("[kulczynski] two-sample hand-computed formula (25/48)", {
  # s1=(10,0,5) A=15, s2=(0,3,5) B=8, C=5
  # kulczynski = 1 - 0.5*(5/15 + 5/8) = 1 - 0.5*(1/3 + 5/8) = 25/48
  mat <- matrix(c(10L, 0L, 5L, 0L, 3L, 5L), nrow = 3, ncol = 2)
  rownames(mat) <- paste0("t", 1:3)
  colnames(mat) <- c("s1", "s2")
  ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
  otu_table(ps) <- sparse_otu_table(otu_table(ps))
  expect_equal(
    as.numeric(sparse_distance(ps, method = "kulczynski")),
    25 / 48,
    tolerance = 1e-12
  )
})

# ---- Manhattan specific ----

test_that("[manhattan] all distances are non-negative", {
  d <- sparse_distance(pair$sparse, method = "manhattan")
  expect_true(all(as.numeric(d) >= 0))
})

test_that("[manhattan] two-sample hand-computed formula (13)", {
  # s1=(10,0,5), s2=(0,3,5): |10-0|+|0-3|+|5-5| = 13
  mat <- matrix(c(10L, 0L, 5L, 0L, 3L, 5L), nrow = 3, ncol = 2)
  rownames(mat) <- paste0("t", 1:3)
  colnames(mat) <- c("s1", "s2")
  ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
  otu_table(ps) <- sparse_otu_table(otu_table(ps))
  expect_equal(
    as.numeric(sparse_distance(ps, method = "manhattan")),
    13,
    tolerance = 1e-12
  )
})

# ---- Euclidean specific ----

test_that("[euclidean] all distances are non-negative", {
  d <- sparse_distance(pair$sparse, method = "euclidean")
  expect_true(all(as.numeric(d) >= 0))
})

test_that("[euclidean] two-sample hand-computed formula (1)", {
  # s1=(1,2,3), s2=(1,2,4): sqrt((0)^2+(0)^2+(1)^2) = 1
  mat <- matrix(c(1L, 2L, 3L, 1L, 2L, 4L), nrow = 3, ncol = 2)
  rownames(mat) <- paste0("t", 1:3)
  colnames(mat) <- c("s1", "s2")
  ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
  otu_table(ps) <- sparse_otu_table(otu_table(ps))
  expect_equal(
    as.numeric(sparse_distance(ps, method = "euclidean")),
    1,
    tolerance = 1e-12
  )
})

# ---- Canberra specific ----

test_that("[canberra] all distances are in [0, 1]", {
  d <- sparse_distance(pair$sparse, method = "canberra")
  expect_true(all(as.numeric(d) >= 0))
  expect_true(all(as.numeric(d) <= 1 + 1e-10))
})

test_that("[canberra] two-sample hand-computed formula (5/6)", {
  # s1=(5,0,3), s2=(0,2,1)
  # |5-0|/(5+0)=1, |0-2|/(0+2)=1, |3-1|/(3+1)=0.5; n_active=3
  # canberra = (1+1+0.5)/3 = 5/6
  mat <- matrix(c(5L, 0L, 3L, 0L, 2L, 1L), nrow = 3, ncol = 2)
  rownames(mat) <- paste0("t", 1:3)
  colnames(mat) <- c("s1", "s2")
  ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
  otu_table(ps) <- sparse_otu_table(otu_table(ps))
  expect_equal(
    as.numeric(sparse_distance(ps, method = "canberra")),
    5 / 6,
    tolerance = 1e-12
  )
})

test_that("[canberra] disjoint samples have distance 1", {
  mat <- matrix(c(5L, 0L, 0L, 3L), nrow = 2)
  rownames(mat) <- c("t1", "t2")
  colnames(mat) <- c("s1", "s2")
  ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
  expect_equal(
    as.numeric(sparse_distance(ps, method = "canberra")),
    1,
    tolerance = 1e-12
  )
})

# ---- Horn specific ----

test_that("[horn] all distances are in [0, 1]", {
  d <- sparse_distance(pair$sparse, method = "horn")
  expect_true(all(as.numeric(d) >= -1e-10))
  expect_true(all(as.numeric(d) <= 1 + 1e-10))
})

test_that("[horn] completely disjoint samples have distance 1", {
  mat <- matrix(c(5L, 0L, 0L, 3L), nrow = 2)
  rownames(mat) <- c("t1", "t2")
  colnames(mat) <- c("s1", "s2")
  ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
  expect_equal(
    as.numeric(sparse_distance(ps, method = "horn")),
    1,
    tolerance = 1e-12
  )
})

# ---- Chord specific ----

test_that("[chord] all distances are in [0, sqrt(2)]", {
  d <- sparse_distance(pair$sparse, method = "chord")
  expect_true(all(as.numeric(d) >= -1e-10))
  expect_true(all(as.numeric(d) <= sqrt(2) + 1e-10))
})

test_that("[chord] two-sample hand-computed formula (sqrt(2)/5)", {
  # s1=(3,4), s2=(4,3): ||s1||=||s2||=5, dot=3*4+4*3=24
  # chord = sqrt(2 - 2*24/25) = sqrt(2/25) = sqrt(2)/5
  mat <- matrix(c(3L, 4L, 4L, 3L), nrow = 2, ncol = 2)
  rownames(mat) <- c("t1", "t2")
  colnames(mat) <- c("s1", "s2")
  ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
  otu_table(ps) <- sparse_otu_table(otu_table(ps))
  expect_equal(
    as.numeric(sparse_distance(ps, method = "chord")),
    sqrt(2) / 5,
    tolerance = 1e-12
  )
})

test_that("[chord] orthogonal samples have max distance sqrt(2)", {
  mat <- matrix(c(5L, 0L, 0L, 3L), nrow = 2)
  rownames(mat) <- c("t1", "t2")
  colnames(mat) <- c("s1", "s2")
  ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
  expect_equal(
    as.numeric(sparse_distance(ps, method = "chord")),
    sqrt(2),
    tolerance = 1e-12
  )
})

# ---- Hellinger specific ----

test_that("[hellinger] all distances are in [0, sqrt(2)]", {
  d <- sparse_distance(pair$sparse, method = "hellinger")
  expect_true(all(as.numeric(d) >= -1e-10))
  expect_true(all(as.numeric(d) <= sqrt(2) + 1e-10))
})

test_that("[hellinger] two-sample hand-computed formula (sqrt(2))", {
  # s1=(4,0), s2=(0,1): disjoint support
  # h1=(1,0), h2=(0,1); euclidean = sqrt(2)
  mat <- matrix(c(4L, 0L, 0L, 1L), nrow = 2, ncol = 2)
  rownames(mat) <- c("t1", "t2")
  colnames(mat) <- c("s1", "s2")
  ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE))
  otu_table(ps) <- sparse_otu_table(otu_table(ps))
  expect_equal(
    as.numeric(sparse_distance(ps, method = "hellinger")),
    sqrt(2),
    tolerance = 1e-12
  )
})

# ---- UniFrac / weighted UniFrac (dispatch to sparse_unifrac) ----
# These check that sparse_distance() dispatches correctly and matches
# phyloseq::distance() when a tree is present, and warns + returns NULL when
# it isn't. Full sparse_unifrac() correctness (rooted/unrooted trees, with/
# without edge lengths, etc.) is covered in test-sparse_unifrac.R.

for (method in TREE_METHODS) {
  local({
    m <- method

    test_that(paste0("[", m, "] sparse_distance matches phyloseq::distance"), {
      expect_equal(
        as.numeric(sparse_distance(pair_tree$sparse, method = m)),
        as.numeric(phyloseq::distance(pair_tree$dense, method = m)),
        tolerance = 1e-8
      )
    })

    test_that(paste0("[", m, "] sparse_distance warns and returns NULL without a tree"), {
      expect_warning(
        result <- sparse_distance(pair$dense, method = m),
        "phy_tree slot is empty"
      )
      expect_null(result)
    })
  })
}
