library(PhyloIgSeq)

# ---- Helpers ----

# Small taxonomy table with two deliberate edge cases:
# - ASV1/ASV2 share identical taxonomy at every real rank (different ASVs)
# - ASV4/ASV5 are NA at Species but already differ at Genus, so a naive
#   "NA matches NA" comparison would understate their distance
make_tax_ps <- function() {
  tox <- rbind(
    c(
      Kingdom = "Bacteria",
      Phylum = "Firmicutes",
      Genus = "Bacteroides",
      Species = "caccae"
    ),
    c(
      Kingdom = "Bacteria",
      Phylum = "Firmicutes",
      Genus = "Bacteroides",
      Species = "caccae"
    ),
    c(
      Kingdom = "Bacteria",
      Phylum = "Firmicutes",
      Genus = "Prevotella",
      Species = "copri"
    ),
    c(
      Kingdom = "Bacteria",
      Phylum = "Bacteroidetes",
      Genus = "Bacteroides",
      Species = NA
    ),
    c(
      Kingdom = "Bacteria",
      Phylum = "Bacteroidetes",
      Genus = "Prevotella",
      Species = NA
    ),
    c(
      Kingdom = "Archaea",
      Phylum = "Euryarchaeota",
      Genus = "Methanobrevibacter",
      Species = "smithii"
    )
  )
  rownames(tox) <- paste0("ASV", 1:6)
  otu <- matrix(
    sample(0:10, nrow(tox) * 3, replace = TRUE),
    nrow = nrow(tox),
    dimnames = list(rownames(tox), paste0("S", 1:3))
  )
  phyloseq(otu_table(otu, taxa_are_rows = TRUE), tax_table(tox))
}

# Taxonomy engineered so one rank (Phylum) covers almost all taxa, which is
# exactly the shape that overflows build_taxonomy_distance_rankwise()'s
# sparse crossprod (see the guard it raises).
make_overflow_ps <- function(n = 60000) {
  tox <- cbind(
    Kingdom = rep("Bacteria", n),
    Phylum = rep("Firmicutes", n),
    Genus = paste0("G", seq_len(n))
  )
  rownames(tox) <- paste0("ASV", seq_len(n))
  phyloseq(
    otu_table(
      matrix(1, nrow = n, ncol = 1, dimnames = list(rownames(tox), "S1")),
      taxa_are_rows = TRUE
    ),
    tax_table(tox)
  )
}

DIST_FUNS <- list(
  rankwise = build_taxonomy_distance_rankwise,
  longest_prefix = build_taxonomy_distance_longest_prefix
)

# ---- Generic structural checks shared by both distance-matrix methods ----

for (method_name in names(DIST_FUNS)) {
  local({
    m <- method_name
    fn <- DIST_FUNS[[m]]
    ps <- make_tax_ps()

    test_that(paste0("[", m, "] returns a dist object with correct size"), {
      d <- fn(ps)
      expect_s3_class(d, "dist")
      expect_equal(length(d), ntaxa(ps) * (ntaxa(ps) - 1L) / 2L)
    })

    test_that(paste0("[", m, "] labels match taxa names"), {
      expect_identical(attr(fn(ps), "Labels"), taxa_names(ps))
    })

    test_that(paste0("[", m, "] distances lie in [0, 1] with zero diagonal"), {
      dm <- as.matrix(fn(ps))
      expect_true(all(dm >= -1e-10 & dm <= 1 + 1e-10))
      expect_true(all(diag(dm) == 0))
    })

    test_that(
      paste0(
        "[",
        m,
        "] identical taxonomy gives a small nonzero distance, not 0"
      ),
      {
        dm <- as.matrix(fn(ps))
        d12 <- dm["ASV1", "ASV2"]
        expect_gt(d12, 0)
        expect_lt(d12, 0.5)
      }
    )

    test_that(
      paste0(
        "[",
        m,
        "] shared NA across different lineages is not a false match"
      ),
      {
        # ASV4 and ASV5 are both NA at Species but already differ at Genus;
        # their distance should be at least as large as two taxa that only
        # differ at the (real) Genus rank, e.g. ASV1 vs ASV3.
        dm <- as.matrix(fn(ps))
        expect_gte(dm["ASV4", "ASV5"], dm["ASV1", "ASV3"] - 1e-10)
      }
    )

    test_that(
      paste0(
        "[",
        m,
        "] maximally different taxa (differ at Kingdom) hit distance 1"
      ),
      {
        dm <- as.matrix(fn(ps))
        expect_equal(dm["ASV1", "ASV6"], 1, tolerance = 1e-10)
      }
    )

    test_that(paste0("[", m, "] rank_weight_base changes the result"), {
      expect_false(isTRUE(all.equal(
        as.numeric(fn(ps, rank_weight_base = 1.1)),
        as.numeric(fn(ps, rank_weight_base = 3))
      )))
    })
  })
}

test_that("rankwise and longest_prefix agree exactly once NA/duplicate labels are disambiguated", {
  ps <- make_tax_ps()
  d_rw <- build_taxonomy_distance_rankwise(ps)
  d_lp <- build_taxonomy_distance_longest_prefix(ps)
  expect_equal(as.numeric(d_rw), as.numeric(d_lp), tolerance = 1e-10)
})

test_that("build_taxonomy_distance_rankwise errors instead of crashing on overflow-prone taxonomy", {
  ps <- make_overflow_ps()
  expect_error(build_taxonomy_distance_rankwise(ps), "hierarchy")
})

# ---- build_taxonomy_tree_hierarchy() ----

test_that("build_taxonomy_tree_hierarchy returns a valid rooted binary phylo object", {
  ps <- make_tax_ps()
  tree <- build_taxonomy_tree_hierarchy(ps)
  expect_s3_class(tree, "phylo")
  expect_true(ape::is.rooted(tree))
  expect_true(ape::is.binary(tree))
  expect_setequal(tree$tip.label, taxa_names(ps))
})

test_that("build_taxonomy_tree_hierarchy patristic distances match the rankwise heuristic (up to the 2x tree-path factor)", {
  ps <- make_tax_ps()
  tree <- build_taxonomy_tree_hierarchy(ps)
  cophenetic_dist <- ape::cophenetic.phylo(tree)[taxa_names(ps), taxa_names(ps)]
  rankwise_dist <- as.matrix(build_taxonomy_distance_rankwise(ps))

  expect_equal(cophenetic_dist / 2, rankwise_dist, tolerance = 1e-8)
})

test_that("build_taxonomy_tree_hierarchy identical taxonomy gives nonzero (not zero) patristic distance", {
  ps <- make_tax_ps()
  tree <- build_taxonomy_tree_hierarchy(ps)
  cophenetic_dist <- ape::cophenetic.phylo(tree)
  expect_gt(cophenetic_dist["ASV1", "ASV2"], 0)
})

test_that("build_taxonomy_tree_hierarchy rank_weight_base changes branch lengths", {
  ps <- make_tax_ps()
  tree_default <- build_taxonomy_tree_hierarchy(ps)
  tree_base3 <- build_taxonomy_tree_hierarchy(ps, rank_weight_base = 3)
  expect_false(isTRUE(all.equal(
    tree_default$edge.length,
    tree_base3$edge.length
  )))
})

test_that("build_taxonomy_tree_hierarchy works with phyloseq::UniFrac", {
  ps <- make_tax_ps()
  phy_tree(ps) <- build_taxonomy_tree_hierarchy(ps)
  uf <- UniFrac(ps, weighted = TRUE)
  expect_s3_class(uf, "dist")
  expect_equal(length(uf), nsamples(ps) * (nsamples(ps) - 1L) / 2L)
})

test_that("build_taxonomy_tree_hierarchy scales to a large number of taxa without error", {
  skip_on_cran()
  ps <- make_overflow_ps(n = 20000)
  tree <- build_taxonomy_tree_hierarchy(ps)
  expect_s3_class(tree, "phylo")
  expect_true(ape::is.binary(tree))
  expect_equal(length(tree$tip.label), ntaxa(ps))
})

# ---- get_taxonomy_tree() ----

METHODS <- c("hierarchy", "rankwise", "longest_prefix")

for (method_name in METHODS) {
  local({
    m <- method_name
    ps <- make_tax_ps()

    test_that(
      paste0(
        "get_taxonomy_tree(method = '",
        m,
        "') returns a tree matching taxa names"
      ),
      {
        tree <- get_taxonomy_tree(ps, method = m)
        expect_s3_class(tree, "phylo")
        expect_setequal(tree$tip.label, taxa_names(ps))
      }
    )
  })
}

test_that("get_taxonomy_tree defaults to the hierarchy method", {
  ps <- make_tax_ps()
  expect_equal(
    get_taxonomy_tree(ps),
    get_taxonomy_tree(ps, method = "hierarchy")
  )
})

test_that("get_taxonomy_tree forwards rank_weight_base to the chosen method", {
  ps <- make_tax_ps()
  tree_default <- get_taxonomy_tree(ps, method = "hierarchy")
  tree_base3 <- get_taxonomy_tree(
    ps,
    method = "hierarchy",
    rank_weight_base = 3
  )
  expect_false(isTRUE(all.equal(
    tree_default$edge.length,
    tree_base3$edge.length
  )))
})
