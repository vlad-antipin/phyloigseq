library(PhyloIgSeq)

# ---- scree_plot() ----

test_that("scree_plot warns and returns NULL for NULL input", {
  expect_warning(result <- scree_plot(NULL), "No eigen values")
  expect_null(result)
})

test_that("scree_plot warns and returns NULL for empty input", {
  expect_warning(result <- scree_plot(numeric(0)), "No eigen values")
  expect_null(result)
})

test_that("scree_plot returns a ggplot for named eigenvalues", {
  eigen_values <- c(Axis.1 = 5, Axis.2 = 3, Axis.3 = 1, Axis.4 = 0.5)
  plt <- scree_plot(eigen_values)
  expect_s3_class(plt, "ggplot")
})

test_that("scree_plot computes percent variability against positive eigenvalues only", {
  eigen_values <- c(Axis.1 = 6, Axis.2 = 2, Axis.3 = 2)
  plt <- scree_plot(eigen_values)
  built <- ggplot2::layer_data(plt)
  # total_var = 10, so proportions are 60/20/20
  expect_equal(sort(built$y), c(20, 20, 60))
})

test_that("scree_plot excludes negative eigenvalues from the variance total but still plots them", {
  eigen_values <- c(Axis.1 = 8, Axis.2 = 2, Axis.3 = -1)
  plt <- scree_plot(eigen_values)
  built <- ggplot2::layer_data(plt)
  # total_var = 10 (only positive values count), Axis.3 still plotted at -10%
  expect_equal(built$y, c(80, 20, -10))
})

test_that("scree_plot truncates to max_nb_comp leading axes", {
  eigen_values <- c(Axis.1 = 5, Axis.2 = 4, Axis.3 = 3, Axis.4 = 2, Axis.5 = 1)
  plt <- scree_plot(eigen_values, max_nb_comp = 2)
  built <- ggplot2::layer_data(plt)
  expect_equal(nrow(built), 2)
})

test_that("scree_plot falls back to integer dim labels when eigen_values is unnamed", {
  eigen_values <- c(5, 3, 1)
  plt <- scree_plot(eigen_values)
  expect_equal(levels(plt$data$dim), c("1", "2", "3"))
})

test_that("scree_plot dim factor levels preserve input order", {
  eigen_values <- c(Axis.2 = 3, Axis.1 = 5, Axis.3 = 1)
  plt <- scree_plot(eigen_values)
  expect_equal(levels(plt$data$dim), c("Axis.2", "Axis.1", "Axis.3"))
})

# ---- get_beta_diversity() ----

make_bd_ps <- function(n_samples = 10, n_taxa = 12, seed = 42) {
  set.seed(seed)
  samp_nms <- paste0("S", seq_len(n_samples))
  taxa_nms <- paste0("T", seq_len(n_taxa))

  mat <- matrix(
    rpois(n_samples * n_taxa, lambda = 20),
    nrow = n_taxa,
    ncol = n_samples,
    dimnames = list(taxa_nms, samp_nms)
  )
  mat[1, ] <- mat[1, ] + 5 # avoid an all-zero taxon
  mat[, 1] <- mat[, 1] + 5 # avoid an all-zero sample

  sdata <- data.frame(
    group = rep(c("A", "B"), length.out = n_samples),
    batch = rep(c("x", "y", "z"), length.out = n_samples),
    row.names = samp_nms
  )

  taxtab <- matrix(
    c(
      rep(c("Firmicutes", "Bacteroidetes", "Proteobacteria"), length.out = n_taxa),
      taxa_nms
    ),
    nrow = n_taxa,
    dimnames = list(taxa_nms, c("Phylum", "taxon_id"))
  )

  phyloseq(
    otu_table(mat, taxa_are_rows = TRUE),
    sample_data(sdata),
    tax_table(taxtab)
  )
}

test_that("get_beta_diversity rejects non-phyloseq input", {
  expect_error(get_beta_diversity(list()), "Need a phyloseq")
})

test_that("get_beta_diversity errors on an invalid method", {
  ps <- make_bd_ps()
  expect_error(get_beta_diversity(ps, method = "bogus"), "Invalid method")
})

test_that("get_beta_diversity requires a model for constrained methods", {
  ps <- make_bd_ps()
  expect_error(
    get_beta_diversity(ps, method = "CCA"),
    "Model is required"
  )
})

test_that("get_beta_diversity has dropped the unused species argument", {
  ps <- make_bd_ps()
  expect_error(
    get_beta_diversity(ps, method = "PCoA", species = TRUE),
    "unused argument"
  )
})

test_that("get_beta_diversity messages (not prints) the unifrac-without-tree fallback", {
  ps <- make_bd_ps()
  expect_message(
    get_beta_diversity(ps, method = "PCoA", dist = "unifrac"),
    "bray-curtis distance selected"
  )
})

test_that("get_beta_diversity dispatches distance-based methods (PCoA/tSNE/UMAP)", {
  ps <- make_bd_ps()
  for (m in c("PCoA", "tSNE", "UMAP")) {
    bd <- get_beta_diversity(
      ps,
      method = m,
      dist = "bray",
      perplexity = 3,
      nb_neighbors = 4
    )
    expect_equal(nrow(bd$coords[[1]]), 10L, info = m)
    expect_equal(bd$method, m, info = m)
    expect_equal(bd$covariates, list(), info = m)
  }
  expect_false(is.null(get_beta_diversity(ps, method = "PCoA")$eigen_values))
})

test_that("get_beta_diversity dispatches abundance-based methods (NMDS/PCA/CA/DCA)", {
  ps <- make_bd_ps()
  for (m in c("NMDS", "PCA", "CA", "DCA")) {
    bd <- suppressMessages(get_beta_diversity(ps, method = m))
    expect_equal(nrow(bd$coords[[1]]), 10L, info = m)
    expect_equal(bd$covariates, list(), info = m)
  }
  expect_null(suppressMessages(get_beta_diversity(ps, method = "NMDS"))$eigen_values)
  expect_false(is.null(get_beta_diversity(ps, method = "PCA")$eigen_values))
})

test_that("get_beta_diversity dispatches constrained methods (CCA/RDA/dbRDA)", {
  ps <- make_bd_ps()
  for (m in c("CCA", "RDA", "dbRDA")) {
    bd <- get_beta_diversity(ps, method = m, model = "group", dist = "bray")
    expect_equal(nrow(bd$coords[[1]]), 10L, info = m)
    expect_gt(length(bd$covariates), 0L, label = m)
    expect_equal(bd$model, "group", info = m)
  }
})

test_that("get_beta_diversity's confounders augment the returned model string", {
  ps <- make_bd_ps()
  bd <- get_beta_diversity(
    ps,
    method = "CCA",
    model = "group",
    confounders = "batch"
  )
  expect_equal(bd$model, "group + Condition(`batch`)")
})

test_that("get_beta_diversity with fit_filter warns and NA-fills for methods without projection", {
  ps <- make_bd_ps(n_samples = 14)
  expect_warning(
    bd <- get_beta_diversity(
      ps,
      method = "dbRDA",
      model = "group",
      dist = "bray",
      fit_filter_name = "batch",
      fit_filter_values = c("x", "y")
    ),
    "Projection of non-fit samples is not supported for dbRDA"
  )
  expect_equal(nrow(bd$coords[[1]]), 14L)
  expect_true(all(is.na(bd$coords[[1]][!bd$sample_data$.is_fit_sample, ])))
})

test_that("get_beta_diversity with fit_filter projects samples for methods that support it (CCA)", {
  ps <- make_bd_ps(n_samples = 14)
  bd <- get_beta_diversity(
    ps,
    method = "CCA",
    model = "group",
    fit_filter_name = "batch",
    fit_filter_values = c("x", "y")
  )
  expect_equal(nrow(bd$coords[[1]]), 14L)
  expect_false(any(is.na(bd$coords[[1]])))
  expect_true(".is_fit_sample" %in% colnames(bd$sample_data))
})

test_that("get_beta_diversity agglomerates by taxrank before ordination", {
  ps <- make_bd_ps()
  bd <- get_beta_diversity(ps, method = "PCoA", taxrank = "Phylum")
  expect_equal(nrow(bd$tax_table), 3L)
  expect_equal(nrow(bd$loadings[[1]]), 3L)
})
