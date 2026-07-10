library(PhyloIgSeq)

# ---- Synthetic PhyloIgSeq fixture ----

make_pis_bd <- function(n_samples = 10, n_taxa = 20, seed = 789) {
  set.seed(seed)
  samp_nms <- paste0("S", seq_len(n_samples))
  taxa_nms <- paste0("T", seq_len(n_taxa))

  anchor_i <- sample(n_samples, n_taxa, replace = TRUE)
  anchor_j <- seq_len(n_taxa)
  extra_i <- sample(n_samples, 15L, replace = TRUE)
  extra_j <- sample(n_taxa, 15L, replace = TRUE)
  i <- c(anchor_i, extra_i)
  j <- c(anchor_j, extra_j)
  ij_key <- paste(i, j, sep = "_")
  keep <- !duplicated(ij_key)
  i <- i[keep]
  j <- j[keep]
  n_obs <- length(i)
  vals <- rnorm(n_obs, sd = 2)

  ig_coating <- data.frame(
    sample_id = samp_nms[i],
    taxon_id = taxa_nms[j],
    slide_z = vals,
    stringsAsFactors = FALSE
  )

  sample_data <- data.frame(
    sample_id = samp_nms,
    group = rep(c("A", "B"), length.out = n_samples),
    stringsAsFactors = FALSE
  )

  tax_table <- data.frame(
    taxon_id = c(taxa_nms, NA_character_),
    Phylum = c(paste0("P", seq_len(n_taxa)), "P_NA"),
    Family = c(paste0("F", seq_len(n_taxa)), "F_NA"),
    stringsAsFactors = FALSE
  )

  pis <- new(
    "PhyloIgSeq",
    ig_coating = ig_coating,
    positive_fraction_name = "pos",
    first_negative_fraction_name = "neg",
    sample_data = sample_data,
    tax_table = tax_table
  )

  list(
    pis = pis,
    samp_nms = samp_nms,
    taxa_nms = taxa_nms,
    n_samples = n_samples,
    n_taxa = n_taxa
  )
}

f <- make_pis_bd()
ps <- PhyloIgSeq_to_phyloseq(
  f$pis,
  score_name = "slide_z",
  imputation_method = "SVD",
  svd_rank = 5L
)
ot <- phyloseq::otu_table(ps)

# ---- Test 6: Distance matrix properties ----

test_that("SVD embedding distance matrix is symmetric", {
  svd_fit <- ot@svd_fit
  emb <- svd_fit$u %*% diag(svd_fit$d)
  dm <- as.matrix(dist(emb))
  expect_true(isSymmetric(dm))
})

test_that("SVD embedding distance matrix has zero diagonal", {
  svd_fit <- ot@svd_fit
  emb <- svd_fit$u %*% diag(svd_fit$d)
  dm <- as.matrix(dist(emb))
  expect_true(all(diag(dm) == 0))
})

test_that("SVD embedding distance matrix is non-negative", {
  svd_fit <- ot@svd_fit
  emb <- svd_fit$u %*% diag(svd_fit$d)
  dm <- as.matrix(dist(emb))
  expect_true(all(dm >= 0))
})

test_that("SVD embedding distance matrix has correct dimensions", {
  svd_fit <- ot@svd_fit
  emb <- svd_fit$u %*% diag(svd_fit$d)
  dm <- as.matrix(dist(emb))
  expect_equal(dim(dm), c(f$n_samples, f$n_samples))
})

# ---- Test 7: get_beta_diversity end-to-end ----

result <- get_beta_diversity(ps, method = "PCoA")

test_that("get_beta_diversity returns a list without error", {
  expect_true(is.list(result))
})

test_that("get_beta_diversity coords[[1]] has n_samples rows", {
  expect_equal(nrow(result$coords[[1]]), f$n_samples)
})

test_that("get_beta_diversity coords[[1]] rownames match sample names", {
  expect_setequal(rownames(result$coords[[1]]), f$samp_nms)
})

test_that("get_beta_diversity coords[[1]] has at least 1 column", {
  expect_gte(ncol(result$coords[[1]]), 1L)
})

test_that("get_beta_diversity loadings[[1]] has n_taxa rows", {
  expect_equal(nrow(result$loadings[[1]]), f$n_taxa)
})

test_that("get_beta_diversity loadings[[1]] rownames match taxa names", {
  expect_setequal(rownames(result$loadings[[1]]), f$taxa_nms)
})

test_that("get_beta_diversity loadings[[1]] col count matches coords col count", {
  expect_equal(ncol(result$loadings[[1]]), ncol(result$coords[[1]]))
})

test_that("get_beta_diversity loadings[[1]] colnames match coords colnames", {
  expect_identical(colnames(result$loadings[[1]]), colnames(result$coords[[1]]))
})

test_that("get_beta_diversity returns non-null eigen.values for PCoA", {
  expect_false(is.null(result$eigen.values))
})

# ---- Test 7b: fit.filter subset path ----

result_ff <- get_beta_diversity(
  ps,
  method = "PCoA",
  fit.filter.name = "group",
  fit.filter.values = "A"
)

test_that("get_beta_diversity with fit.filter returns coords for ALL samples", {
  expect_equal(nrow(result_ff$coords[[1]]), f$n_samples)
  expect_setequal(rownames(result_ff$coords[[1]]), f$samp_nms)
})

test_that("get_beta_diversity with fit.filter marks fit samples correctly", {
  expect_true(".is.fit.sample" %in% colnames(result_ff$sample.data))
})
