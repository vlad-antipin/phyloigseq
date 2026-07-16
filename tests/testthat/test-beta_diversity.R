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

# ---- stat_beta_diversity() ----

make_bd_fit <- function(n_samples = 24, seed = 42) {
  ps <- make_bd_ps(n_samples = n_samples, seed = seed)
  sdata <- as(phyloseq::sample_data(ps), "data.frame")
  # facet1 is blocked in pairs so it doesn't align with group's period-2
  # alternation (which would make group constant within a facet1 level).
  sdata$facet1 <- rep(c("F1", "F1", "F2", "F2"), length.out = n_samples)
  sdata$facet2 <- rep(c("G1", "G2"), each = n_samples / 2)
  set.seed(seed)
  sdata$cont_var <- stats::rnorm(n_samples)
  phyloseq::sample_data(ps) <- phyloseq::sample_data(sdata)
  get_beta_diversity(ps, method = "PCoA", dist = "bray")
}

test_that("stat_beta_diversity warns and returns NULL for an unknown label name", {
  fit <- make_bd_fit()
  expect_warning(
    result <- stat_beta_diversity(fit, label_name = "nope"),
    "Wrong label name"
  )
  expect_null(result)
})

test_that("stat_beta_diversity runs PERMANOVA for a categorical label, no facet", {
  fit <- make_bd_fit()
  result <- stat_beta_diversity(fit, label_name = "group")
  expect_equal(result$stat, "permanova")
  expect_match(result$p_value[[1]], "PERMANOVA p")
  expect_null(result$p_value_df)
  expect_length(result$p_value_raw, 1L)
  expect_null(names(result$p_value_raw)) # unfaceted: unnamed, like the pre-refactor behavior
})

test_that("stat_beta_diversity runs envfit for a continuous label", {
  fit <- make_bd_fit()
  result <- stat_beta_diversity(fit, label_name = "cont_var")
  expect_equal(result$stat, "envfit")
  expect_match(result$p_value[[1]], "Correlation p")
})

test_that("stat_beta_diversity wrap-mode facets, one PERMANOVA per facet level", {
  fit <- make_bd_fit()
  result <- stat_beta_diversity(
    fit,
    facet_mode = "wrap",
    facet = "facet1",
    label_name = "group"
  )
  expect_setequal(names(result$p_value_raw), c("F1", "F2"))
  expect_length(result$test_result, 2L)
})

test_that("stat_beta_diversity grid mode returns a p_value_df with one row per cell", {
  fit <- make_bd_fit()
  # suppressMessages(): grid cells are small enough that vegan::adonis2() falls
  # back to complete permutation enumeration and notes so via message().
  result <- suppressMessages(stat_beta_diversity(
    fit,
    facet_mode = "grid",
    facet_row = "facet1",
    facet_col = "facet2",
    label_name = "group"
  ))
  expect_equal(nrow(result$p_value_df), 4L)
  expect_setequal(result$p_value_df$facet_row, c("F1", "F2"))
  expect_setequal(result$p_value_df$facet_col, c("G1", "G2"))
  expect_null(result$p_value)
})

test_that("stat_beta_diversity strata restricts the PERMANOVA permutation scheme", {
  fit <- make_bd_fit()
  result <- stat_beta_diversity(
    fit,
    label_name = "group",
    strata_name = "batch"
  )
  expect_match(result$p_value[[1]], "restricted by batch")
})

test_that("stat_beta_diversity pairwise adds a BH-adjusted pairwise table", {
  fit <- make_bd_fit()
  result <- stat_beta_diversity(fit, label_name = "batch", pairwise = TRUE)
  expect_setequal(
    colnames(result$pairwise_df),
    c("group1", "group2", "R2", "p_raw", "p_adj")
  )
  expect_equal(nrow(result$pairwise_df), 3L) # choose(3, 2) levels of batch
})

test_that("stat_beta_diversity comp = NULL uses every retained ordination axis", {
  fit <- make_bd_fit()
  result <- stat_beta_diversity(fit, comp = NULL, label_name = "group")
  expect_equal(result$dim_used, paste("all", ncol(fit$coords[[1]])))
})

test_that("stat_beta_diversity has dropped the dead facet.name back-compat alias", {
  fit <- make_bd_fit()
  expect_error(
    stat_beta_diversity(fit, facet.name = "facet1", label_name = "group"),
    "unused argument"
  )
})

test_that("stat_beta_diversity wrap-mode skips only the degenerate facet, keeping the others (bug fix)", {
  fit <- make_bd_fit()
  sdata <- fit$sample_data
  # Force facet2's "G1" level to a single, constant group value.
  sdata$group[sdata$facet2 == "G1"] <- "A"
  fit$sample_data <- sdata
  expect_warning(
    result <- stat_beta_diversity(
      fit,
      facet_mode = "wrap",
      facet = "facet2",
      label_name = "group"
    ),
    "doesn't vary"
  )
  expect_false(is.null(result)) # no longer aborts entirely
  expect_equal(names(result$p_value_raw), "G2")
  expect_length(result$test_result, 1L)
})

test_that("stat_beta_diversity grid mode skips only the degenerate cell, keeping the others", {
  fit <- make_bd_fit()
  sdata <- fit$sample_data
  sdata$group[sdata$facet1 == "F1" & sdata$facet2 == "G1"] <- "A"
  fit$sample_data <- sdata
  # suppressMessages(): grid cells are small enough that vegan::adonis2() falls
  # back to complete permutation enumeration and notes so via message().
  expect_warning(
    result <- suppressMessages(stat_beta_diversity(
      fit,
      facet_mode = "grid",
      facet_row = "facet1",
      facet_col = "facet2",
      label_name = "group"
    )),
    "doesn't vary"
  )
  expect_equal(nrow(result$p_value_df), 3L)
})

# ---- beta_diversity() ----
#
# suppressWarnings(): plot_beta_diversity()'s pre-existing geom_point(point_aes, ...)
# call passes an unused `text` aesthetic (for the app's plotly tooltip), triggering a
# harmless "Ignoring unknown aesthetics" cosmetic warning on every call. Out of scope
# for this batch (plot_beta_diversity() itself is Batch 14e/14f); muffled here so it
# doesn't drown out these tests' own signal.

test_that("beta_diversity returns a ggplot", {
  plt <- suppressWarnings(beta_diversity(make_bd_ps(), group = "group"))
  expect_s3_class(plt, "ggplot")
})

test_that("beta_diversity colors the plot by group even when stat = FALSE (bug fix)", {
  plt <- suppressWarnings(beta_diversity(make_bd_ps(), group = "group"))
  expect_setequal(as.character(unique(plt$layers[[1]]$data$label)), c("A", "B"))
})

test_that("beta_diversity forwards comp to the plot, not just to stat_beta_diversity (bug fix)", {
  ps <- make_bd_ps()
  fit <- get_beta_diversity(ps, method = "PCoA", dist = "bray")
  plt_default <- suppressWarnings(beta_diversity(ps, group = "group"))
  plt_23 <- suppressWarnings(beta_diversity(ps, group = "group", comp = c(2, 3)))
  expect_equal(plt_default$layers[[1]]$data$Comp1, unname(fit$coords[[1]][, 1]))
  expect_equal(plt_23$layers[[1]]$data$Comp1, unname(fit$coords[[1]][, 2]))
})

test_that("beta_diversity with stat = TRUE runs without erroring (regression: stat.beta.dispersion[[1]]$test.result indexing bug)", {
  plt <- suppressWarnings(beta_diversity(make_bd_ps(), group = "group", stat = TRUE))
  expect_s3_class(plt, "ggplot")
})

test_that("beta_diversity has dropped the dead species argument and the full_ prefix", {
  expect_false("species" %in% names(formals(beta_diversity)))
  expect_false(exists("full_beta_diversity", where = asNamespace("PhyloIgSeq"), inherits = FALSE))
})

test_that("beta_diversity forwards ... to plot_beta_diversity", {
  plt <- suppressWarnings(beta_diversity(make_bd_ps(), group = "group", ellipses = TRUE))
  expect_true(any(vapply(plt$layers, function(l) inherits(l$stat, "StatEllipse"), logical(1))))
})

# ---- animate_by_variable() ----

test_that("animate_by_variable returns a gganim/gif_image object by default", {
  ps <- make_bd_ps()
  plt <- suppressWarnings(beta_diversity(
    ps,
    group = "group",
    animation.variable.name = "group"
  ))
  anim <- animate_by_variable(
    plt,
    "group",
    nframes = 6,
    fps = 3,
    width = 200,
    height = 200,
    res = 72
  )
  expect_true(inherits(anim, "gif_image") || inherits(anim, "gganim"))
})

test_that("animate_by_variable returns NULL when return_anim = FALSE", {
  ps <- make_bd_ps()
  plt <- suppressWarnings(beta_diversity(
    ps,
    group = "group",
    animation.variable.name = "group"
  ))
  anim <- animate_by_variable(
    plt,
    "group",
    return_anim = FALSE,
    nframes = 6,
    fps = 3,
    width = 200,
    height = 200,
    res = 72
  )
  expect_null(anim)
})

test_that("animate_by_variable writes a file to save_path when given", {
  ps <- make_bd_ps()
  plt <- suppressWarnings(beta_diversity(
    ps,
    group = "group",
    animation.variable.name = "group"
  ))
  out <- tempfile(fileext = ".gif")
  on.exit(unlink(out), add = TRUE)
  animate_by_variable(
    plt,
    "group",
    save_path = out,
    nframes = 6,
    fps = 3,
    width = 200,
    height = 200,
    res = 72
  )
  expect_true(file.exists(out))
})
