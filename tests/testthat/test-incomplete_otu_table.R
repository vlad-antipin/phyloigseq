library(PhyloIgSeq)

# ---- Synthetic test-data helpers ----

# Build an Incomplete object + softImpute fit from random triplets.
# Returns a list with everything needed to construct/verify the class.
make_inc_data <- function(n_samples = 12, n_taxa = 25, seed = 42) {
  set.seed(seed)
  n_obs <- 35L
  i <- sample(n_samples, n_obs, replace = TRUE)
  j <- sample(n_taxa,   n_obs, replace = TRUE)
  # deduplicate: keep first occurrence of each (i,j) pair
  ij_key <- paste(i, j, sep = "_")
  keep   <- !duplicated(ij_key)
  i      <- i[keep]; j <- j[keep]; n_obs <- sum(keep)
  vals   <- rnorm(n_obs, sd = 2)

  # Anchor (n_samples, n_taxa) so Incomplete() infers the full dimensions.
  # Remove any existing entry at that position first to avoid a duplicate.
  anchor_ij <- paste(n_samples, n_taxa, sep = "_")
  drop       <- paste(i, j, sep = "_") == anchor_ij
  i <- c(i[!drop], n_samples)
  j <- c(j[!drop], n_taxa)
  vals <- c(vals[!drop], rnorm(1L, sd = 2))
  n_obs <- length(i)

  samp_names <- paste0("S", seq_len(n_samples))
  tax_names  <- paste0("T", seq_len(n_taxa))

  X_inc <- softImpute::Incomplete(i, j, vals)
  dimnames(X_inc) <- list(samp_names, tax_names)

  # Dense matrix with NA for unobserved positions
  X_dense <- matrix(NA_real_, n_samples, n_taxa,
                    dimnames = list(samp_names, tax_names))
  for (k in seq_len(n_obs)) X_dense[i[k], j[k]] <- vals[k]

  col_means            <- colMeans(X_dense, na.rm = TRUE)
  col_means[is.nan(col_means)] <- 0
  X_centered           <- sweep(X_dense, 2, col_means, "-")

  r   <- 3L
  fit <- softImpute::softImpute(X_centered, rank.max = r, lambda = 1,
                                type = "svd", maxit = 200L)

  list(
    X_inc      = X_inc,
    svd_fit    = list(u = fit$u, d = fit$d, v = fit$v),
    col_means  = col_means,
    i          = i, j = j, vals = vals,
    n_samples  = n_samples, n_taxa = n_taxa, r = r,
    samp_names = samp_names, tax_names = tax_names
  )
}

make_ot <- function(d = make_inc_data()) {
  incomplete_otu_table(
    X_inc     = d$X_inc,
    svd_fit   = d$svd_fit,
    col_means = d$col_means
  )
}

# Canonical objects for most tests
d  <- make_inc_data()
ot <- make_ot(d)

# ---- Test 1: Class hierarchy ----

test_that("incomplete_otu_table has correct S4 class hierarchy", {
  expect_s4_class(ot, "incomplete_otu_table")
  expect_s4_class(ot, "sparse_otu_table")
  expect_s4_class(ot, "otu_table")
  expect_true(is(ot, "matrix"))
})

# ---- Test 2: Dimensions & names ----

test_that("dim returns (n_samples, n_taxa) for taxa_are_rows = FALSE", {
  expect_equal(dim(ot), c(d$n_samples, d$n_taxa))
  expect_equal(nrow(ot), d$n_samples)
  expect_equal(ncol(ot), d$n_taxa)
})

test_that("taxa_names and sample_names are correct for taxa_are_rows = FALSE", {
  expect_identical(taxa_names(ot),   d$tax_names)
  expect_identical(sample_names(ot), d$samp_names)
  expect_equal(ntaxa(ot),    d$n_taxa)
  expect_equal(nsamples(ot), d$n_samples)
})

test_that("rownames / colnames are sample / taxon names respectively", {
  expect_identical(rownames(ot), d$samp_names)
  expect_identical(colnames(ot), d$tax_names)
})

# ---- Test 3: is.na ----

test_that("is.na has correct dimensions", {
  na_mat <- is.na(ot)
  expect_identical(dim(na_mat), dim(ot))
  expect_true(is.logical(na_mat))
})

test_that("is.na returns FALSE for every observed entry", {
  na_mat <- is.na(ot)
  for (k in seq_along(d$i)) {
    expect_false(na_mat[d$i[k], d$j[k]])
  }
})

test_that("is.na: total observed count matches number of unique triplets", {
  na_mat <- is.na(ot)
  expect_equal(sum(!na_mat), length(d$i))
})

test_that("is.na returns TRUE for a known unobserved position", {
  na_mat <- is.na(ot)
  # Find a position that was NOT observed
  obs_set <- paste(d$i, d$j, sep = "_")
  for (ri in seq_len(d$n_samples)) {
    for (ci in seq_len(d$n_taxa)) {
      if (!paste(ri, ci, sep = "_") %in% obs_set) {
        expect_true(na_mat[ri, ci])
        return()
      }
    }
  }
})

# ---- Test 4: as(x, "matrix") ----

test_that("as(x, 'matrix') produces correct dimensions and dimnames", {
  mat <- as(ot, "matrix")
  expect_identical(dim(mat), dim(ot))
  expect_identical(dimnames(mat), dimnames(ot))
})

test_that("as(x, 'matrix') has NA at every unobserved position", {
  mat    <- as(ot, "matrix")
  na_mat <- is.na(ot)
  expect_true(all(is.na(mat[na_mat])))
})

test_that("as(x, 'matrix') has the original value at every observed position", {
  mat <- as(ot, "matrix")
  for (k in seq_along(d$i)) {
    expect_equal(mat[d$i[k], d$j[k]], d$vals[k])
  }
})

test_that("explicit observed zero is not NA in as(x, 'matrix')", {
  # Build a tiny object with an explicit zero as an observed score
  X_inc0 <- softImpute::Incomplete(c(1L, 2L), c(1L, 2L), c(0.0, 1.0))
  dimnames(X_inc0) <- list(c("S1","S2"), c("T1","T2"))
  X_d0 <- matrix(NA_real_, 2L, 2L); X_d0[1,1] <- 0; X_d0[2,2] <- 1
  fit0 <- softImpute::softImpute(X_d0, rank.max = 1L, lambda = 0.1, type = "svd")
  ot0  <- incomplete_otu_table(X_inc0, list(u=fit0$u, d=fit0$d, v=fit0$v))
  mat0 <- as(ot0, "matrix")
  expect_false(is.na(mat0[1L, 1L]))   # observed 0.0 → not NA
  expect_equal(mat0[1L, 1L], 0.0)
  expect_true(is.na(mat0[1L, 2L]))    # unobserved → NA
  expect_true(is.na(mat0[2L, 1L]))    # unobserved → NA
})

test_that("as.matrix() (S3) is identical to as(x, 'matrix')", {
  expect_identical(as.matrix(ot), as(ot, "matrix"))
})

# ---- Test 5: SVD fit ----

test_that("svd_fit slots have correct shapes", {
  svd <- ot@svd_fit
  expect_equal(dim(svd$u), c(d$n_samples, d$r))
  expect_equal(length(svd$d), d$r)
  expect_equal(dim(svd$v), c(d$n_taxa, d$r))
  expect_true(all(svd$d >= 0))
})

test_that("svd_fit reconstruction is close at observed positions", {
  svd     <- ot@svd_fit
  X_hat   <- svd$u %*% diag(svd$d) %*% t(svd$v)            # centered approx
  X_hat_orig <- sweep(X_hat, 2L, ot@col_means, "+")          # de-center
  mat     <- as(ot, "matrix")
  obs     <- !is.na(mat)
  mse     <- mean((mat[obs] - X_hat_orig[obs])^2)
  expect_lt(mse, 20)   # generous bound; real data RMSE ~ 1-2
})

test_that("col_means slot has correct length and names", {
  expect_length(ot@col_means, d$n_taxa)
  expect_named(ot@col_means)
  expect_identical(names(ot@col_means), d$tax_names)
})

test_that("sample coordinates from svd_fit have correct shape", {
  coords <- ot@svd_fit$u %*% diag(ot@svd_fit$d)
  expect_equal(dim(coords), c(d$n_samples, d$r))
})

# ---- Test 9: rowSums / colSums warning ----

test_that("rowSums emits a warning and returns named vector of length n_samples", {
  expect_warning(rs <- rowSums(ot))
  expect_named(rs)
  expect_length(rs, d$n_samples)
  expect_identical(names(rs), d$samp_names)
})

test_that("colSums emits a warning and returns named vector of length n_taxa", {
  expect_warning(cs <- colSums(ot))
  expect_named(cs)
  expect_length(cs, d$n_taxa)
  expect_identical(names(cs), d$tax_names)
})

test_that("sample_sums emits a warning and returns named vector", {
  expect_warning(ss <- sample_sums(ot))
  expect_named(ss)
  expect_length(ss, d$n_samples)
  expect_identical(names(ss), d$samp_names)
})

test_that("taxa_sums emits a warning and returns named vector", {
  expect_warning(ts <- taxa_sums(ot))
  expect_named(ts)
  expect_length(ts, d$n_taxa)
  expect_identical(names(ts), d$tax_names)
})

# rowSums / colSums on sparse_otu_table (non-incomplete) must still NOT warn
test_that("rowSums on sparse_otu_table does not warn", {
  m <- matrix(c(1,0,0,2), 2, 2)
  sp <- sparse_otu_table(phyloseq::otu_table(m, taxa_are_rows = TRUE))
  expect_no_warning(rowSums(sp))
})

# ---- Test 10: Subsetting drops SVD fit ----

test_that("[ row-subset returns standard otu_table (not incomplete_otu_table)", {
  sub <- ot[1:3, ]
  expect_s4_class(sub, "otu_table")
  expect_false(is(sub, "incomplete_otu_table"))
  expect_equal(nrow(sub), 3L)
  expect_identical(rownames(sub), d$samp_names[1:3])
})

test_that("[ column-subset returns standard otu_table with NAs preserved", {
  sub <- ot[, 1:5]
  expect_s4_class(sub, "otu_table")
  expect_equal(ncol(sub), 5L)
  expect_identical(colnames(sub), d$tax_names[1:5])
  # Values at observed positions in the subset match original
  mat <- as(ot, "matrix")
  expect_equal(as(sub, "matrix"), mat[, 1:5])
})

test_that("[ flat extraction via logical index returns plain vector", {
  na_mat <- is.na(ot)
  vals   <- ot[!na_mat]
  expect_true(is.numeric(vals))
  expect_false(is.matrix(vals))
  expect_equal(sort(vals), sort(d$vals))
})
