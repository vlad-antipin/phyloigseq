library(PhyloIgSeq)

# ---- Helpers ----

# Random sparse count matrix + a phylo tree, wrapped in a phyloseq object.
# Mirrors make_pair()/make_ps() in the sibling distance test files, but with
# a phy_tree attached so sparse_unifrac() has something to walk.
make_tree_pair <- function(
  n_taxa = 15,
  n_samples = 8,
  sparsity = 0.6,
  taxa_are_rows = TRUE,
  seed = 42,
  tree = NULL
) {
  set.seed(seed)
  mat <- matrix(0L, n_taxa, n_samples)
  n_nz <- max(2L, round((1 - sparsity) * n_taxa * n_samples))
  mat[sample(n_taxa * n_samples, n_nz)] <- sample(
    1L:1000L,
    n_nz,
    replace = TRUE
  )
  taxa <- paste0("ASV", seq_len(n_taxa))
  rownames(mat) <- taxa
  colnames(mat) <- paste0("S", seq_len(n_samples))

  if (is.null(tree)) {
    tree <- ape::rtree(n_taxa, tip.label = taxa)
  }

  if (!taxa_are_rows) {
    mat <- t(mat)
  }

  ps_dense <- phyloseq(
    otu_table(mat, taxa_are_rows = taxa_are_rows),
    phy_tree(tree)
  )
  ps_sparse <- as_sparse_phyloseq(ps_dense)
  list(dense = ps_dense, sparse = ps_sparse, tree = tree)
}

ref_unifrac <- function(ps, method) phyloseq::distance(ps, method = method)

METHODS <- c("unifrac", "wunifrac")

pair <- make_tree_pair()
pair_t <- make_tree_pair(taxa_are_rows = FALSE)

# ---- Generic structural checks (rooted tree, real branch lengths) ----
# local() captures the loop variable so each block closes over its own `m`.

for (method in METHODS) {
  local({
    m <- method

    test_that(paste0("[", m, "] returns a dist object"), {
      expect_s3_class(sparse_unifrac(pair$sparse, method = m), "dist")
    })

    test_that(paste0("[", m, "] dist has correct size (n*(n-1)/2)"), {
      n <- nsamples(pair$sparse)
      expect_equal(
        length(sparse_unifrac(pair$sparse, method = m)),
        n * (n - 1L) / 2L
      )
    })

    test_that(paste0("[", m, "] dist labels match sample names"), {
      expect_identical(
        attr(sparse_unifrac(pair$sparse, method = m), "Labels"),
        sample_names(pair$sparse)
      )
    })

    test_that(paste0("[", m, "] distance matrix is symmetric"), {
      dm <- as.matrix(sparse_unifrac(pair$sparse, method = m))
      expect_equal(dm, t(dm))
    })

    test_that(paste0("[", m, "] diagonal of distance matrix is zero"), {
      dm <- as.matrix(sparse_unifrac(pair$sparse, method = m))
      expect_true(all(diag(dm) == 0))
    })

    test_that(paste0("[", m, "] matches phyloseq::distance (tar = TRUE)"), {
      expect_equal(
        as.numeric(sparse_unifrac(pair$sparse, method = m)),
        as.numeric(ref_unifrac(pair$dense, method = m)),
        tolerance = 1e-8
      )
    })

    test_that(paste0("[", m, "] matches phyloseq::distance (tar = FALSE)"), {
      expect_equal(
        as.numeric(sparse_unifrac(pair_t$sparse, method = m)),
        as.numeric(ref_unifrac(pair_t$dense, method = m)),
        tolerance = 1e-8
      )
    })

    test_that(paste0("[", m, "] accepts plain dense phyloseq"), {
      expect_equal(
        as.numeric(sparse_unifrac(pair$dense, method = m)),
        as.numeric(ref_unifrac(pair$dense, method = m)),
        tolerance = 1e-8
      )
    })

    test_that(paste0("[", m, "] identical samples have distance 0"), {
      mat <- matrix(c(1L, 2L, 3L, 1L, 2L, 3L), nrow = 3)
      rownames(mat) <- paste0("ASV", 1:3)
      colnames(mat) <- c("s1", "s2")
      tr <- ape::rtree(3, tip.label = rownames(mat))
      ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE), phy_tree(tr))
      expect_equal(
        as.numeric(sparse_unifrac(ps, method = m)),
        0,
        tolerance = 1e-10
      )
    })

    test_that(
      paste0(
        "[",
        m,
        "] single-sample errors the same way as phyloseq::distance"
      ),
      {
        # Neither implementation can form a pair from one sample; both
        # bottom out in combn(samples, 2) with the same "n < m" error, so
        # this documents parity rather than a bug.
        mat <- matrix(c(1L, 2L, 3L), nrow = 3)
        rownames(mat) <- paste0("ASV", 1:3)
        colnames(mat) <- "s1"
        tr <- ape::rtree(3, tip.label = rownames(mat))
        ps <- phyloseq(otu_table(mat, taxa_are_rows = TRUE), phy_tree(tr))
        expect_error(sparse_unifrac(ps, method = m), "n < m")
        expect_error(ref_unifrac(ps, method = m), "n < m")
      }
    )

    test_that(
      paste0(
        "[",
        m,
        "] result is unaffected by OTU-table/tree-tip ordering mismatch"
      ),
      {
        shuffled <- pair$sparse
        otu_table(shuffled) <- otu_table(pair$sparse)[
          sample(taxa_names(pair$sparse)),
        ]
        expect_equal(
          as.numeric(sparse_unifrac(shuffled, method = m)),
          as.numeric(sparse_unifrac(pair$sparse, method = m)),
          tolerance = 1e-10
        )
      }
    )
  })
}

# ---- Range sanity checks ----

test_that("[unifrac] all distances are in [0, 1]", {
  d <- as.numeric(sparse_unifrac(pair$sparse, method = "unifrac"))
  expect_true(all(d >= -1e-10))
  expect_true(all(d <= 1 + 1e-10))
})

test_that("[wunifrac] all distances are non-negative", {
  d <- as.numeric(sparse_unifrac(pair$sparse, method = "wunifrac"))
  expect_true(all(d >= -1e-10))
})

# ---- No tree present ----

for (method in METHODS) {
  local({
    m <- method

    test_that(
      paste0("[", m, "] warns and returns NULL when phy_tree is absent"),
      {
        ps <- phyloseq(otu_table(pair$sparse))
        expect_warning(
          result <- sparse_unifrac(ps, method = m),
          "phy_tree slot is empty"
        )
        expect_null(result)
      }
    )
  })
}

# ---- Rooted tree, no informative branch lengths (unit edge lengths) ----

test_that("[unifrac] matches phyloseq::distance with unit edge lengths", {
  tr <- pair$tree
  tr$edge.length <- rep(1, length(tr$edge.length))
  p <- make_tree_pair(tree = tr)
  expect_equal(
    as.numeric(sparse_unifrac(p$sparse, method = "unifrac")),
    as.numeric(ref_unifrac(p$dense, method = "unifrac")),
    tolerance = 1e-8
  )
})

test_that("[wunifrac] matches phyloseq::distance with unit edge lengths", {
  tr <- pair$tree
  tr$edge.length <- rep(1, length(tr$edge.length))
  p <- make_tree_pair(tree = tr)
  expect_equal(
    as.numeric(sparse_unifrac(p$sparse, method = "wunifrac")),
    as.numeric(ref_unifrac(p$dense, method = "wunifrac")),
    tolerance = 1e-8
  )
})

# ---- Unrooted tree ----
# phyloseq::UniFrac roots an unrooted tree at a *random* tip, so its result
# on the original unrooted tree isn't a valid reference. Instead we check:
#  (1) sparse_unifrac() warns and midpoint-roots the tree, and
#  (2) that midpoint-rooted result is self-consistent with calling
#      sparse_unifrac() directly on the manually midpoint-rooted tree, and
#  (3) matches phyloseq::distance() computed on that *same* rooted tree,
#      which validates the core algorithm independent of rooting policy.

test_that("[unifrac] unrooted tree is midpoint-rooted with a warning", {
  tr_unrooted <- ape::unroot(pair$tree)
  expect_false(ape::is.rooted(tr_unrooted))
  p_unrooted <- make_tree_pair(tree = tr_unrooted)

  expect_warning(
    sparse_unifrac(p_unrooted$sparse, method = "unifrac"),
    "unrooted"
  )
})

for (method in METHODS) {
  local({
    m <- method

    test_that(
      paste0(
        "[",
        m,
        "] unrooted-tree result matches explicit midpoint-rooting"
      ),
      {
        tr_unrooted <- ape::unroot(pair$tree)
        p_unrooted <- make_tree_pair(tree = tr_unrooted)
        res_auto <- suppressWarnings(
          sparse_unifrac(p_unrooted$sparse, method = m)
        )

        tr_mid <- phytools::midpoint_root(tr_unrooted)
        p_mid <- make_tree_pair(tree = tr_mid)
        res_manual <- sparse_unifrac(p_mid$sparse, method = m)

        expect_equal(as.numeric(res_auto), as.numeric(res_manual))
      }
    )

    test_that(
      paste0(
        "[",
        m,
        "] midpoint-rooted result matches phyloseq::distance on same tree"
      ),
      {
        tr_unrooted <- ape::unroot(pair$tree)
        tr_mid <- phytools::midpoint_root(tr_unrooted)
        p_mid <- make_tree_pair(tree = tr_mid)

        expect_equal(
          as.numeric(sparse_unifrac(p_mid$sparse, method = m)),
          as.numeric(ref_unifrac(p_mid$dense, method = m)),
          tolerance = 1e-8
        )
      }
    )
  })
}
