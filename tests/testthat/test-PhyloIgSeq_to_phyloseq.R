library(PhyloIgSeq)

# ---- Synthetic PhyloIgSeq fixture ----

make_pis <- function(n_samples = 8, n_taxa = 15, seed = 123) {
  set.seed(seed)
  samp_nms <- paste0("S", seq_len(n_samples))
  taxa_nms <- paste0("T", seq_len(n_taxa))

  # One observation per taxon (anchors) + extra random ones; deduplicate.
  anchor_i <- sample(n_samples, n_taxa, replace = TRUE)
  anchor_j <- seq_len(n_taxa)
  extra_i <- sample(n_samples, 12L, replace = TRUE)
  extra_j <- sample(n_taxa, 12L, replace = TRUE)
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
  # Include some NA-score rows to verify na.omit is applied
  ig_coating <- rbind(
    ig_coating,
    data.frame(
      sample_id = samp_nms[1],
      taxon_id = "T_NA",
      slide_z = NA_real_,
      stringsAsFactors = FALSE
    )
  )

  sample_data <- data.frame(
    sample_id = samp_nms,
    group = rep(c("A", "B"), length.out = n_samples),
    stringsAsFactors = FALSE
  )

  tax_table <- data.frame(
    taxon_id = c(taxa_nms, NA_character_), # one NA taxon_id row → filtered out
    Phylum = c(paste0("P", seq_len(n_taxa)), "P_NA"),
    Family = c(paste0("F", seq_len(n_taxa)), "F_NA"),
    stringsAsFactors = FALSE
  )

  list(
    pis = new(
      "PhyloIgSeq",
      ig_coating = ig_coating,
      positive_fraction_name = "pos",
      first_negative_fraction_name = "neg",
      sample_data = sample_data,
      tax_table = tax_table
    ),
    samp_nms = samp_nms,
    taxa_nms = taxa_nms,
    i = i,
    j = j,
    vals = vals,
    n_samples = n_samples,
    n_taxa = n_taxa
  )
}

f <- make_pis()
ps <- PhyloIgSeq_to_phyloseq(
  f$pis,
  score_name = "slide_z",
  imputation_method = "SVD",
  svd_rank = 3L
)
ot <- phyloseq::otu_table(ps)

# ---- Test 1: return type ----

test_that("SVD path returns a phyloseq object", {
  expect_s4_class(ps, "phyloseq")
})

test_that("OTU table is an incomplete_otu_table", {
  expect_s4_class(ot, "incomplete_otu_table")
  expect_s4_class(ot, "sparse_otu_table")
  expect_s4_class(ot, "otu_table")
})

# ---- Test 2: dimensions & names ----

test_that("sample_names match PhyloIgSeq sample_data$sample_id", {
  expect_identical(sort(phyloseq::sample_names(ps)), sort(f$samp_nms))
})

test_that("taxa_names are the unique taxon_ids from ig_coating (NA rows excluded)", {
  expect_setequal(phyloseq::taxa_names(ps), f$taxa_nms)
})

test_that("nsamples and ntaxa are correct", {
  expect_equal(phyloseq::nsamples(ps), f$n_samples)
  expect_equal(phyloseq::ntaxa(ps), f$n_taxa)
})

test_that("OTU table dim matches (n_samples × n_taxa) for taxa_are_rows = FALSE", {
  expect_equal(nrow(ot), f$n_samples)
  expect_equal(ncol(ot), f$n_taxa)
})

# ---- Test 3: is.na ----

test_that("is.na has correct dimensions", {
  na_mat <- is.na(ot)
  expect_identical(dim(na_mat), dim(ot))
})

test_that("is.na returns FALSE for known observed entries from ig_coating", {
  na_mat <- is.na(ot)
  for (k in seq_along(f$i)) {
    si <- f$samp_nms[f$i[k]]
    ti <- f$taxa_nms[f$j[k]]
    expect_false(
      na_mat[si, ti],
      label = paste0("is.na[", si, ",", ti, "] should be FALSE")
    )
  }
})

test_that("total observed count equals unique (sample, taxon) pairs", {
  na_mat <- is.na(ot)
  expect_equal(sum(!na_mat), length(f$i))
})

test_that("is.na returns TRUE for a position not in ig_coating", {
  na_mat <- is.na(ot)
  obs_set <- paste(f$samp_nms[f$i], f$taxa_nms[f$j], sep = "_")
  found <- FALSE
  for (ri in seq_len(f$n_samples)) {
    for (ci in seq_len(f$n_taxa)) {
      key <- paste(f$samp_nms[ri], f$taxa_nms[ci], sep = "_")
      if (!key %in% obs_set) {
        expect_true(na_mat[f$samp_nms[ri], f$taxa_nms[ci]])
        found <- TRUE
        break
      }
    }
    if (found) break
  }
  expect_true(found, label = "Expected at least one unobserved position")
})

# ---- Test 4: as(x, "matrix") round-trips original values ----

test_that("as(x, 'matrix') has NA at unobserved positions", {
  mat <- as(ot, "matrix")
  na_mat <- is.na(ot)
  expect_true(all(is.na(mat[na_mat])))
})

test_that("as(x, 'matrix') has original ig_coating values at observed positions", {
  mat <- as(ot, "matrix")
  for (k in seq_along(f$i)) {
    si <- f$samp_nms[f$i[k]]
    ti <- f$taxa_nms[f$j[k]]
    expect_equal(
      mat[si, ti],
      f$vals[k],
      label = paste0("mat[", si, ",", ti, "]")
    )
  }
})

test_that("as.matrix() (S3) matches as(x, 'matrix')", {
  expect_identical(as.matrix(ot), as(ot, "matrix"))
})

# ---- Test 5: SVD fit slots ----

test_that("svd_fit slot is present with u / d / v components", {
  svd <- ot@svd_fit
  expect_true(is.list(svd))
  expect_true(all(c("u", "d", "v") %in% names(svd)))
})

test_that("svd_fit shapes are (n_samples × r), r, (n_taxa × r)", {
  svd <- ot@svd_fit
  r <- length(svd$d)
  expect_equal(nrow(svd$u), f$n_samples)
  expect_equal(ncol(svd$u), r)
  expect_equal(nrow(svd$v), f$n_taxa)
  expect_equal(ncol(svd$v), r)
  expect_true(all(svd$d >= 0))
})

test_that("col_means has correct length and names", {
  cm <- ot@col_means
  expect_length(cm, f$n_taxa)
  expect_named(cm)
  expect_setequal(names(cm), f$taxa_nms)
})

# ---- Test 8: legacy imputation methods still work ----

# `f$pis` includes one all-NA-taxon row ("T_NA") specifically to exercise `na.omit()`
# on the SVD path (see above); on the dense path that used to produce an all-NA
# otu_table column, which VIM::kNN() errored on ungracefully ("subscript out of
# bounds") rather than skipping (known issue #8). `to_wider_ig_score()` now excludes
# any all-NA-score taxon up front (with a `warning()` naming it, see its own tests in
# test-ig_score.R), and `.PhyloIgSeq_to_phyloseq_dense()` keeps `tax_table` in sync
# with that exclusion — so `f$pis` (T_NA included) can be used directly everywhere
# below instead of a manually-pruned fixture, and the fixture-level workaround that
# used to live here is gone.

test_that("imputation_method = NULL returns standard phyloseq with dense otu_table", {
  expect_warning(
    ps_null <- PhyloIgSeq_to_phyloseq(
      f$pis,
      score_name = "slide_z",
      imputation_method = NULL
    ),
    "T_NA"
  )
  ot_null <- phyloseq::otu_table(ps_null)
  expect_s4_class(ot_null, "otu_table")
  expect_false(is(ot_null, "incomplete_otu_table"))
})

test_that("imputation_method = 'Replace NA with 0' returns standard dense otu_table", {
  expect_warning(
    ps_r0 <- PhyloIgSeq_to_phyloseq(
      f$pis,
      score_name = "slide_z",
      imputation_method = "Replace NA with 0"
    ),
    "T_NA"
  )
  ot_r0 <- phyloseq::otu_table(ps_r0)
  expect_s4_class(ot_r0, "otu_table")
  expect_false(is(ot_r0, "incomplete_otu_table"))
  # All NAs should be replaced with 0
  mat_r0 <- as(ot_r0, "matrix")
  expect_false(anyNA(mat_r0))
})

test_that("imputation_method = 'KNN' returns standard dense otu_table with no NAs", {
  expect_warning(
    ps_knn <- PhyloIgSeq_to_phyloseq(
      f$pis,
      score_name = "slide_z",
      imputation_method = "KNN",
      nb_neighbors = 3
    ),
    "T_NA"
  )
  ot_knn <- phyloseq::otu_table(ps_knn)
  expect_s4_class(ot_knn, "otu_table")
  expect_false(is(ot_knn, "incomplete_otu_table"))
  mat_knn <- as(ot_knn, "matrix")
  expect_false(anyNA(mat_knn))
  # Observed entries are untouched by kNN imputation
  for (k in seq_along(f$i)) {
    si <- f$samp_nms[f$i[k]]
    ti <- f$taxa_nms[f$j[k]]
    expect_equal(mat_knn[si, ti], f$vals[k])
  }
})

test_that("the all-NA taxon dropped by to_wider_ig_score() is also absent from tax_table", {
  ps_knn <- suppressWarnings(PhyloIgSeq_to_phyloseq(
    f$pis,
    score_name = "slide_z",
    imputation_method = "KNN",
    nb_neighbors = 3
  ))
  expect_false("T_NA" %in% phyloseq::taxa_names(ps_knn))
  expect_false("T_NA" %in% rownames(phyloseq::tax_table(ps_knn)))
})

test_that("imputation_method = 'Central Tendency' fills NAs with the per-taxon median", {
  expect_warning(
    ps_ct <- PhyloIgSeq_to_phyloseq(
      f$pis,
      score_name = "slide_z",
      imputation_method = "Central Tendency",
      central_tendency = "median"
    ),
    "T_NA"
  )
  ot_ct <- phyloseq::otu_table(ps_ct)
  expect_s4_class(ot_ct, "otu_table")
  expect_false(is(ot_ct, "incomplete_otu_table"))
  mat_ct <- as(ot_ct, "matrix")
  expect_false(anyNA(mat_ct))

  # Every unobserved (sample, taxon) cell gets filled with that taxon's observed median
  observed_by_taxon <- split(f$vals, f$taxa_nms[f$j])
  for (ti in f$taxa_nms) {
    expected_median <- median(observed_by_taxon[[ti]])
    unobserved_samples <- setdiff(
      f$samp_nms,
      f$samp_nms[f$i[f$taxa_nms[f$j] == ti]]
    )
    for (si in unobserved_samples) {
      expect_equal(
        mat_ct[si, ti],
        expected_median,
        label = paste0("mat_ct[", si, ",", ti, "]")
      )
    }
  }
})
