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

test_that("get_beta_diversity's dbRDA loadings are real per-taxon values, not vegan's all-NA species scores", {
  # capscale() is never given comm =, so scores(fit, display = "species")
  # would just be an all-NA "Dim1".."DimN" placeholder -- .get_beta_diversity_
  # constrained() falls back to correlating raw abundances against the site
  # scores instead (same approach the standalone PCoA branch already uses).
  ps <- make_bd_ps()
  bd <- get_beta_diversity(ps, method = "dbRDA", model = "group", dist = "bray")
  expect_false(any(is.na(bd$loadings[[1]])))
  expect_setequal(rownames(bd$loadings[[1]]), taxa_names(ps))
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

test_that("get_beta_diversity with fit_filter projects dbRDA's constrained and residual axes correctly", {
  ps <- make_bd_ps(n_samples = 14)
  for (d in c("bray", "euclidean")) {
    bd <- get_beta_diversity(
      ps,
      method = "dbRDA",
      model = "group",
      dist = d,
      fit_filter_name = "batch",
      fit_filter_values = c("x", "y")
    )
    expect_equal(nrow(bd$coords[[1]]), 14L, info = d)
    expect_false(any(is.na(bd$coords[[1]])), info = d)
    # Projected (fit-subset) rows must reproduce capscale's own site scores
    # exactly, for both scaling 1 and 2 -- regression test for
    # .dbrda_projected_scores(), vegan's predict(type = "wa") being
    # unconditionally blocked for capscale objects.
    own_1 <- vegan::scores(bd$fit, display = "sites", choices = seq_len(ncol(bd$coords[[1]])), scaling = 1)
    own_2 <- vegan::scores(bd$fit, display = "sites", choices = seq_len(ncol(bd$coords[[2]])), scaling = 2)
    fit_names <- bd$fit_sample_names
    expect_equal(unname(bd$coords[[1]][fit_names, ]), unname(own_1), tolerance = 1e-6, ignore_attr = TRUE, info = d)
    expect_equal(unname(bd$coords[[2]][fit_names, ]), unname(own_2), tolerance = 1e-6, ignore_attr = TRUE, info = d)
    expect_false(isTRUE(all.equal(bd$coords[[1]], bd$coords[[2]])), info = d)
  }
})

test_that("get_beta_diversity with fit_filter projects confounders-only RDA/dbRDA (no free predictor)", {
  ps <- make_bd_ps(n_samples = 14)
  for (m in c("RDA", "dbRDA")) {
    bd <- get_beta_diversity(
      ps,
      method = m,
      confounders = "group",
      dist = "bray",
      fit_filter_name = "batch",
      fit_filter_values = c("x", "y")
    )
    expect_equal(nrow(bd$coords[[1]]), 14L, info = m)
    expect_false(any(is.na(bd$coords[[1]])), info = m)
    # Every axis here is "residual" (post-confounder) -- .pcca_fitted_all()'s
    # hand-derived fit, in place of predict(type = "lc")/"working", which
    # only support model = "CCA"/"CA" (i.e. a free predictor, not
    # confounders-only). Regression check: fit-subset rows must still
    # reproduce the fit's own site scores exactly.
    own_1 <- vegan::scores(bd$fit, display = "sites", choices = seq_len(ncol(bd$coords[[1]])), scaling = 1)
    fit_names <- bd$fit_sample_names
    expect_equal(unname(bd$coords[[1]][fit_names, ]), unname(own_1), tolerance = 1e-6, ignore_attr = TRUE, info = m)
  }
})

test_that("get_beta_diversity with fit_filter warns and NA-fills CCA when there's no free predictor to project with", {
  ps <- make_bd_ps(n_samples = 14)
  expect_warning(
    bd <- get_beta_diversity(
      ps,
      method = "CCA",
      confounders = "group",
      dist = "bray",
      fit_filter_name = "batch",
      fit_filter_values = c("x", "y")
    ),
    "Projection of non-fit samples is not supported for CCA"
  )
  expect_equal(nrow(bd$coords[[1]]), 14L)
  expect_true(all(is.na(bd$coords[[1]][!bd$sample_data$.is_fit_sample, ])))
})

test_that("get_beta_diversity with fit_filter warns and NA-fills a fit-subset sample with a missing model variable, instead of erroring", {
  ps <- make_bd_ps(n_samples = 14)
  sd <- as(sample_data(ps), "data.frame")
  sd$group[2] <- NA # S2, in the fit subset (batch == "y")
  sample_data(ps) <- sample_data(sd)

  warnings_seen <- character(0)
  bd <- withCallingHandlers(
    get_beta_diversity(
      ps,
      method = "RDA",
      model = "group",
      fit_filter_name = "batch",
      fit_filter_values = c("x", "y")
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_match(
    warnings_seen,
    "fit-subset sample\\(s\\) excluded from the model fit",
    all = FALSE
  )
  expect_match(
    warnings_seen,
    "sample\\(s\\) excluded from projection",
    all = FALSE
  )
  expect_equal(nrow(bd$coords[[1]]), 14L)
  # S2 is excluded from the fit itself, but its "wa" (constrained-axis) score
  # doesn't need the predictor value at all -- only its residual axis, which
  # does, ends up NA (same partial-NA pattern as a held-out sample missing
  # the same variable).
  expect_false(is.na(bd$coords[[1]]["S2", "RDA1"]))
  expect_true(all(is.na(bd$coords[[1]]["S2", -1])))
  expect_false(any(is.na(bd$coords[[1]][setdiff(rownames(bd$coords[[1]]), "S2"), ])))
})

test_that("get_beta_diversity with fit_filter partially NA-fills a held-out sample with a missing model variable", {
  ps <- make_bd_ps(n_samples = 14)
  sd <- as(sample_data(ps), "data.frame")
  sd$group[3] <- NA # S3, held out (batch == "z")
  sample_data(ps) <- sample_data(sd)

  expect_warning(
    bd <- get_beta_diversity(
      ps,
      method = "RDA",
      model = "group",
      fit_filter_name = "batch",
      fit_filter_values = c("x", "y")
    ),
    "sample\\(s\\) excluded from projection"
  )
  # "wa" (constrained-axis) scores don't need the predictor value at all, so
  # S3's constrained axis is still real; only its residual axis (which does
  # need it, to residualize against the fitted constrained part) is NA.
  expect_false(is.na(bd$coords[[1]]["S3", "RDA1"]))
  expect_true(all(is.na(bd$coords[[1]]["S3", -1])))
  expect_false(any(is.na(bd$coords[[1]][setdiff(rownames(bd$coords[[1]]), "S3"), ])))
})

test_that("get_beta_diversity with fit_filter NA-fills a held-out sample missing a confounder, for confounders-only RDA/dbRDA", {
  ps <- make_bd_ps(n_samples = 14)
  sd <- as(sample_data(ps), "data.frame")
  sd$group[3] <- NA # S3, held out; "group" is the sole confounder here
  sample_data(ps) <- sample_data(sd)

  for (m in c("RDA", "dbRDA")) {
    expect_warning(
      bd <- get_beta_diversity(
        ps,
        method = m,
        confounders = "group",
        dist = "bray",
        fit_filter_name = "batch",
        fit_filter_values = c("x", "y")
      ),
      "sample\\(s\\) excluded from projection",
      info = m
    )
    expect_true(all(is.na(bd$coords[[1]]["S3", ])), info = m)
    expect_false(any(is.na(bd$coords[[1]][setdiff(rownames(bd$coords[[1]]), "S3"), ])), info = m)
  }
})

test_that("get_beta_diversity with fit_filter projects samples for methods that support it (CCA)", {
  ps <- make_bd_ps(n_samples = 14)
  # CCA's constrained axis is exactly projectable for every sample; its
  # residual axes are chi-square weighted and have no equivalent closed-form
  # projection (unlike RDA's, see .rda_residual_scores()), so those stay NA
  # (with a warning) for out-of-fit samples.
  expect_warning(
    bd <- get_beta_diversity(
      ps,
      method = "CCA",
      model = "group",
      fit_filter_name = "batch",
      fit_filter_values = c("x", "y")
    ),
    "Projection of non-fit samples is not supported for CCA"
  )
  expect_equal(nrow(bd$coords[[1]]), 14L)
  expect_false(any(is.na(bd$coords[[1]][, "CCA1"])))
  expect_true(all(is.na(bd$coords[[1]][!bd$sample_data$.is_fit_sample, -1])))
  expect_true(".is_fit_sample" %in% colnames(bd$sample_data))
})

test_that("get_beta_diversity with fit_filter projects RDA's residual axes correctly too", {
  ps <- make_bd_ps(n_samples = 14)
  bd <- get_beta_diversity(
    ps,
    method = "RDA",
    model = "group",
    fit_filter_name = "batch",
    fit_filter_values = c("x", "y")
  )
  expect_equal(nrow(bd$coords[[1]]), 14L)
  expect_false(any(is.na(bd$coords[[1]])))
  # Projected (fit-subset) rows must reproduce the fit's own site scores
  # exactly, for both scaling 1 and 2 -- regression test for a bug where the
  # projection branch ignored the `scaling` argument entirely.
  own_1 <- vegan::scores(bd$fit, display = "sites", choices = seq_len(ncol(bd$coords[[1]])), scaling = 1)
  own_2 <- vegan::scores(bd$fit, display = "sites", choices = seq_len(ncol(bd$coords[[2]])), scaling = 2)
  fit_names <- bd$fit_sample_names
  expect_equal(unname(bd$coords[[1]][fit_names, ]), unname(own_1), tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(unname(bd$coords[[2]][fit_names, ]), unname(own_2), tolerance = 1e-8, ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(bd$coords[[1]], bd$coords[[2]])))
})

test_that("get_beta_diversity agglomerates by taxrank before ordination", {
  ps <- make_bd_ps()
  bd <- get_beta_diversity(ps, method = "PCoA", taxrank = "Phylum")
  expect_equal(nrow(bd$tax_table), 3L)
  expect_equal(nrow(bd$loadings[[1]]), 3L)
})

# ---- get_beta_diversity(dist = "euclidean"): sparse/dense density gating ----
# make_bd_ps() itself is unusable here: rpois(lambda = 20) essentially never
# returns 0, so it's already ~100% dense regardless of the sparsity this
# section wants to control.

make_bd_ps_density <- function(
  n_samples = 10,
  n_taxa = 60,
  sparsity = 0.9,
  seed = 42
) {
  set.seed(seed)
  samp_nms <- paste0("S", seq_len(n_samples))
  taxa_nms <- paste0("T", seq_len(n_taxa))
  mat <- matrix(0L, n_taxa, n_samples, dimnames = list(taxa_nms, samp_nms))
  n_nz <- max(2L, round((1 - sparsity) * n_taxa * n_samples))
  mat[sample(n_taxa * n_samples, n_nz)] <- sample(1L:5000L, n_nz, replace = TRUE)
  mat[1, ] <- mat[1, ] + 5 # avoid an all-zero taxon
  mat[, 1] <- mat[, 1] + 5 # avoid an all-zero sample

  sdata <- data.frame(
    group = rep(c("A", "B"), length.out = n_samples),
    row.names = samp_nms
  )
  taxtab <- matrix(taxa_nms, nrow = n_taxa, dimnames = list(taxa_nms, "taxon_id"))

  phyloseq(
    otu_table(mat, taxa_are_rows = TRUE),
    sample_data(sdata),
    tax_table(taxtab)
  )
}

# Reference PCoA-on-standardized-euclidean, independent of which internal
# branch get_beta_diversity() takes (dense decostand()+vegdist() vs the
# sparse euclidean_sparse(standardize = TRUE) kernel).
ref_pcoa_euclidean_coords <- function(ps, n_axes) {
  otu_mat <- as(otu_table(reverseASV(ps)), "matrix")
  ref <- vegan::vegdist(
    vegan::decostand(
      otu_mat[, apply(otu_mat, 2, var, na.rm = TRUE) > 0],
      "standardize"
    ),
    method = "euclidean"
  )
  fit <- vegan::wcmdscale(ref, eig = TRUE)
  vegan::scores(fit, display = "sites", choices = seq_len(n_axes))
}

test_that(".otu_table_density counts non-zero cells for a dense otu_table", {
  ps <- make_bd_ps_density(sparsity = 0.9)
  mat <- as(otu_table(ps), "matrix")
  expect_equal(.otu_table_density(otu_table(ps)), sum(mat != 0) / length(mat))
})

test_that(".otu_table_density reads nnz directly from a sparse_otu_table", {
  ps <- make_bd_ps_density(sparsity = 0.9)
  ps_sparse <- as_sparse_phyloseq(ps)
  mat <- as(otu_table(ps), "matrix")
  expect_equal(
    .otu_table_density(otu_table(ps_sparse)),
    sum(mat != 0) / length(mat)
  )
})

test_that("get_beta_diversity(dist='euclidean') matches the dense reference below the sparse-density threshold", {
  ps <- make_bd_ps_density(sparsity = 0.9)
  expect_lt(.otu_table_density(otu_table(ps)), .EUCLIDEAN_SPARSE_DENSITY_THRESHOLD)
  bd <- get_beta_diversity(ps, method = "PCoA", dist = "euclidean")
  n_axes <- ncol(bd$coords[[1]])
  ref_coords <- ref_pcoa_euclidean_coords(ps, n_axes)
  cors <- sapply(
    seq_len(n_axes),
    function(k) abs(cor(bd$coords[[1]][, k], ref_coords[, k]))
  )
  expect_equal(cors, rep(1, n_axes), tolerance = 1e-6)
})

test_that("get_beta_diversity(dist='euclidean') matches the dense reference above the sparse-density threshold (e.g. post-CLR)", {
  ps <- make_bd_ps_density(sparsity = 0) # fully dense, mimics a CLR-transformed table
  expect_gt(.otu_table_density(otu_table(ps)), .EUCLIDEAN_SPARSE_DENSITY_THRESHOLD)
  bd <- get_beta_diversity(ps, method = "PCoA", dist = "euclidean")
  n_axes <- ncol(bd$coords[[1]])
  ref_coords <- ref_pcoa_euclidean_coords(ps, n_axes)
  cors <- sapply(
    seq_len(n_axes),
    function(k) abs(cor(bd$coords[[1]][, k], ref_coords[, k]))
  )
  expect_equal(cors, rep(1, n_axes), tolerance = 1e-6)
})

# `sign_invariant_cor()` is used below to compare two ordinations' site
# scores axis-by-axis: eigen-decomposition-based methods (PCA/PCoA, RDA/dbRDA)
# are only defined up to an arbitrary sign flip per axis, so a plain
# `expect_equal()` on the raw coordinates would be a false negative.
sign_invariant_cor <- function(a, b, n_axes) {
  sapply(seq_len(n_axes), function(k) abs(cor(a[, k], b[, k])))
}

test_that("PCA and PCoA + Euclidean agree exactly, up to axis sign (same math, different route)", {
  ps <- make_bd_ps_density(sparsity = 0.9)
  expect_lt(.otu_table_density(otu_table(ps)), .EUCLIDEAN_SPARSE_DENSITY_THRESHOLD)
  bd_pca <- get_beta_diversity(ps, method = "PCA")
  bd_pcoa <- get_beta_diversity(ps, method = "PCoA", dist = "euclidean")
  n_axes <- min(ncol(bd_pca$coords[[1]]), ncol(bd_pcoa$coords[[1]]))
  cors <- sign_invariant_cor(bd_pca$coords[[1]], bd_pcoa$coords[[1]], n_axes)
  expect_equal(cors, rep(1, n_axes), tolerance = 1e-6)
})

test_that("RDA and dbRDA + Euclidean agree exactly, up to axis sign (RDA's scale=TRUE matches dbRDA's standardized euclidean)", {
  ps <- make_bd_ps_density(sparsity = 0.9)
  expect_lt(.otu_table_density(otu_table(ps)), .EUCLIDEAN_SPARSE_DENSITY_THRESHOLD)
  bd_rda <- get_beta_diversity(ps, method = "RDA", model = "group")
  bd_dbrda <- get_beta_diversity(ps, method = "dbRDA", model = "group", dist = "euclidean")
  n_axes <- min(ncol(bd_rda$coords[[1]]), ncol(bd_dbrda$coords[[1]]))
  cors <- sign_invariant_cor(bd_rda$coords[[1]], bd_dbrda$coords[[1]], n_axes)
  expect_equal(cors, rep(1, n_axes), tolerance = 1e-6)
  expect_equal(
    as.numeric(bd_rda$eigen_values[seq_len(n_axes)]),
    as.numeric(bd_dbrda$eigen_values[seq_len(n_axes)]),
    tolerance = 1e-6
  )
})

test_that("RDA and dbRDA + Euclidean still agree above the sparse-density threshold (e.g. post-CLR)", {
  ps <- make_bd_ps_density(sparsity = 0) # fully dense, mimics a CLR-transformed table
  expect_gt(.otu_table_density(otu_table(ps)), .EUCLIDEAN_SPARSE_DENSITY_THRESHOLD)
  bd_rda <- get_beta_diversity(ps, method = "RDA", model = "group")
  bd_dbrda <- get_beta_diversity(ps, method = "dbRDA", model = "group", dist = "euclidean")
  n_axes <- min(ncol(bd_rda$coords[[1]]), ncol(bd_dbrda$coords[[1]]))
  cors <- sign_invariant_cor(bd_rda$coords[[1]], bd_dbrda$coords[[1]], n_axes)
  expect_equal(cors, rep(1, n_axes), tolerance = 1e-6)
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

# ---- plot_beta_diversity() ----
# Kept suppressWarnings() around the same pre-existing "unknown aesthetics"
# geom_point() warning noted above.

make_bd_constrained_fit <- function(n_samples = 14, seed = 42) {
  ps <- make_bd_ps(n_samples = n_samples, seed = seed)
  get_beta_diversity(ps, method = "RDA", model = "group", dist = "bray")
}

test_that("plot_beta_diversity plots Comp1/Comp2 from the first scaling by default", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(fit))
  expect_s3_class(plt, "ggplot")
  expect_equal(plt$layers[[1]]$data$Comp1, unname(fit$coords[[1]][, 1]))
  expect_equal(plt$layers[[1]]$data$Comp2, unname(fit$coords[[1]][, 2]))
})

test_that("plot_beta_diversity selects the requested scaling", {
  fit <- make_bd_constrained_fit()
  plt <- suppressWarnings(plot_beta_diversity(fit, scaling = 2))
  expect_equal(plt$layers[[1]]$data$Comp1, unname(fit$coords[[2]][, 1]))
})

muffle_unknown_aes <- function(w) {
  if (grepl("unknown aesthetics", conditionMessage(w))) {
    invokeRestart("muffleWarning")
  }
}

test_that("plot_beta_diversity warns and falls back to scaling 1 for an out-of-range scaling", {
  fit <- make_bd_fit()
  withCallingHandlers(
    expect_warning(
      plt <- plot_beta_diversity(fit, scaling = 99),
      "Wrong scaling, scaling 1 is used"
    ),
    warning = muffle_unknown_aes
  )
  expect_equal(plt$layers[[1]]$data$Comp1, unname(fit$coords[[1]][, 1]))
})

test_that("plot_beta_diversity warns and resets to c(1, 2) for invalid comp", {
  fit <- make_bd_fit()
  withCallingHandlers(
    expect_warning(
      plt <- plot_beta_diversity(fit, comp = c(1, 2, 3)),
      "Wrong components, forced to first two components"
    ),
    warning = muffle_unknown_aes
  )
  expect_equal(plt$layers[[1]]$data$Comp1, unname(fit$coords[[1]][, 1]))
})

test_that("plot_beta_diversity plots the requested comp axes", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(fit, comp = c(2, 3)))
  expect_equal(plt$layers[[1]]$data$Comp1, unname(fit$coords[[1]][, 2]))
  expect_equal(plt$layers[[1]]$data$Comp2, unname(fit$coords[[1]][, 3]))
})

test_that("plot_beta_diversity filters samples by label_levels and reorders factor levels", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(
    fit,
    label_name = "facet1",
    label_levels = c("F2", "F1")
  ))
  expect_equal(levels(plt$layers[[1]]$data$label), c("F2", "F1"))
  expect_equal(nrow(plt$layers[[1]]$data), nrow(fit$sample_data))
})

test_that("plot_beta_diversity filters samples by shape_levels", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(
    fit,
    shape_name = "facet1",
    shape_levels = "F1"
  ))
  expect_equal(nrow(plt$layers[[1]]$data), sum(fit$sample_data$facet1 == "F1"))
})

test_that("plot_beta_diversity filters samples by facet_levels in wrap mode", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(
    fit,
    facet = "facet2",
    facet_levels = "G1"
  ))
  expect_equal(nrow(plt$layers[[1]]$data), sum(fit$sample_data$facet2 == "G1"))
})

test_that("plot_beta_diversity filters samples by facet_row_levels/facet_col_levels in grid mode", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(
    fit,
    facet_mode = "grid",
    facet_row = "facet1",
    facet_row_levels = "F1",
    facet_col = "facet2"
  ))
  expect_equal(nrow(plt$layers[[1]]$data), sum(fit$sample_data$facet1 == "F1"))
})

test_that("plot_beta_diversity removes samples with NA in a graphical variable when remove_na_from_plot = TRUE", {
  fit <- make_bd_fit()
  fit$sample_data$facet1[1] <- NA
  plt <- suppressWarnings(plot_beta_diversity(
    fit,
    label_name = "facet1",
    remove_na_from_plot = TRUE
  ))
  expect_equal(nrow(plt$layers[[1]]$data), nrow(fit$sample_data) - 1L)
})

test_that("plot_beta_diversity negates the first component's coordinates when reverse_dim1 = TRUE", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(fit, reverse_dim1 = TRUE))
  expect_equal(plt$layers[[1]]$data$Comp1, -unname(fit$coords[[1]][, 1]))
  expect_equal(plt$layers[[1]]$data$Comp2, unname(fit$coords[[1]][, 2]))
})

test_that("plot_beta_diversity negates the second component's coordinates when reverse_dim2 = TRUE", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(fit, reverse_dim2 = TRUE))
  expect_equal(plt$layers[[1]]$data$Comp2, -unname(fit$coords[[1]][, 2]))
})

test_that("plot_beta_diversity adds a StatEllipse layer when ellipses = TRUE and label is a factor", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(fit, label_name = "facet1", ellipses = TRUE))
  expect_true(any(vapply(plt$layers, function(l) inherits(l$stat, "StatEllipse"), logical(1))))
})

test_that("plot_beta_diversity does not add ellipses when label is numeric", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(fit, label_name = "cont_var", ellipses = TRUE))
  expect_false(any(vapply(plt$layers, function(l) inherits(l$stat, "StatEllipse"), logical(1))))
})

test_that("plot_beta_diversity fills ellipses with color when fill_ellipses = TRUE", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(
    fit,
    label_name = "facet1",
    ellipses = TRUE,
    fill_ellipses = TRUE
  ))
  ellipse_layer <- Filter(function(l) inherits(l$stat, "StatEllipse"), plt$layers)[[1]]
  expect_equal(ellipse_layer$geom_params$outline.type %||% NULL, NULL)
  expect_s3_class(ellipse_layer$geom, "GeomPolygon")
})

test_that("plot_beta_diversity facet_wraps by facet in wrap mode", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(fit, facet = "facet2"))
  expect_s3_class(plt$facet, "FacetWrap")
})

test_that("plot_beta_diversity facet_grids by facet_row x facet_col in grid mode", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(
    fit,
    facet_mode = "grid",
    facet_row = "facet1",
    facet_col = "facet2"
  ))
  expect_s3_class(plt$facet, "FacetGrid")
})

test_that("plot_beta_diversity populates hover_text from hover_variables", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(fit, hover_variables = "facet1"))
  expect_true(all(grepl("^facet1: (F1|F2)<br>$", plt$layers[[1]]$data$hover_text)))
})

test_that("plot_beta_diversity maps size_name to the size aesthetic and rescales point size", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(fit, size_name = "cont_var"))
  # size aesthetic is rescaled to a point-size range by scale_size(), so check
  # rank order is preserved against the raw variable rather than exact values.
  expect_equal(rank(ggplot2::layer_data(plt)$size), rank(fit$sample_data$cont_var))
  expect_true(any(vapply(plt$scales$scales, function(s) "size" %in% s$aesthetics, logical(1))))
})

test_that("plot_beta_diversity maps shape_name to the shape aesthetic for a factor column", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(fit, shape_name = "facet1"))
  expect_true("shape" %in% names(plt$layers[[1]]$data))
  expect_setequal(as.character(plt$layers[[1]]$data$shape), c("F1", "F2"))
})

test_that("plot_beta_diversity embeds a single p-value in the subtitle when there's no facet", {
  fit <- make_bd_fit()
  stat <- stat_beta_diversity(fit, label_name = "facet1")
  plt <- suppressWarnings(plot_beta_diversity(fit, stat_beta_dispersion = stat))
  expect_true(grepl(strsplit(stat$p_value[[1]], "\n")[[1]][1], plt$labels$subtitle, fixed = TRUE))
})

test_that("plot_beta_diversity overrides label_name with stat_beta_dispersion$label_name", {
  fit <- make_bd_fit()
  stat <- stat_beta_diversity(fit, label_name = "facet1")
  plt <- suppressWarnings(plot_beta_diversity(
    fit,
    label_name = "facet2",
    stat_beta_dispersion = stat
  ))
  expect_setequal(as.character(unique(plt$layers[[1]]$data$label)), c("F1", "F2"))
})

test_that("plot_beta_diversity annotates grid-mode facets with p_value_df text", {
  fit <- make_bd_fit()
  stat <- suppressMessages(stat_beta_diversity(
    fit,
    label_name = "group",
    facet_mode = "grid",
    facet_row = "facet1",
    facet_col = "facet2"
  ))
  plt <- suppressWarnings(plot_beta_diversity(
    fit,
    facet_mode = "grid",
    facet_row = "facet1",
    facet_col = "facet2",
    stat_beta_dispersion = stat
  ))
  text_layer <- Filter(
    function(l) inherits(l$geom, "GeomText") && "p_label" %in% names(l$data),
    plt$layers
  )
  expect_length(text_layer, 1L)
  expect_equal(nrow(text_layer[[1]]$data), nrow(stat$p_value_df))
})

test_that("plot_beta_diversity draws two point layers with different alpha when a fit_filter is active", {
  ps <- make_bd_ps(n_samples = 14)
  fit <- suppressWarnings(get_beta_diversity(
    ps,
    method = "CCA",
    model = "group",
    fit_filter_name = "batch",
    fit_filter_values = c("x", "y")
  ))
  plt <- suppressWarnings(plot_beta_diversity(fit, projected_alpha = 0.1, point_alpha = 0.9))
  point_layers <- Filter(function(l) inherits(l$geom, "GeomPoint"), plt$layers)
  expect_length(point_layers, 2L)
  expect_setequal(
    vapply(point_layers, function(l) l$aes_params$alpha, numeric(1)),
    c(0.1, 0.9)
  )
})

test_that("plot_beta_diversity draws loading arrows and filters short ones by arrow_cutoff_load", {
  fit <- make_bd_fit()
  plt_all <- suppressWarnings(plot_beta_diversity(fit, biplot_loadings = TRUE))
  segment_all <- Filter(function(l) inherits(l$geom, "GeomSegment"), plt_all$layers)
  expect_length(segment_all, 1L)
  expect_equal(nrow(segment_all[[1]]$data), sum(!is.na(fit$loadings[[1]][, 1])))

  plt_cut <- suppressWarnings(plot_beta_diversity(
    fit,
    biplot_loadings = TRUE,
    arrow_cutoff_load = 0.99
  ))
  segment_cut <- Filter(function(l) inherits(l$geom, "GeomSegment"), plt_cut$layers)
  expect_lt(nrow(segment_cut[[1]]$data), nrow(segment_all[[1]]$data))
})

test_that("plot_beta_diversity draws covariate arrows for a constrained fit", {
  fit <- make_bd_constrained_fit()
  plt <- suppressWarnings(plot_beta_diversity(fit, biplot_covariates = TRUE))
  segment_layers <- Filter(function(l) inherits(l$geom, "GeomSegment"), plt$layers)
  expect_length(segment_layers, 1L)
  expect_equal(nrow(segment_layers[[1]]$data), nrow(fit$covariates[[1]]))
})

test_that("plot_beta_diversity labels arrows with geom_text when arrow_labels = TRUE and repel = FALSE", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(
    fit,
    biplot_loadings = TRUE,
    arrow_labels = TRUE,
    repel = FALSE
  ))
  text_layers <- Filter(function(l) inherits(l$geom, "GeomText"), plt$layers)
  expect_length(text_layers, 1L)
  expect_setequal(text_layers[[1]]$data$Names, rownames(na.omit(fit$loadings[[1]])))
})

test_that("plot_beta_diversity labels arrows with ggrepel when repel = TRUE", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(
    fit,
    biplot_loadings = TRUE,
    arrow_labels = TRUE,
    repel = TRUE
  ))
  expect_true(any(vapply(
    plt$layers,
    function(l) inherits(l$geom, "GeomTextRepel"),
    logical(1)
  )))
})

test_that("plot_beta_diversity labels arrows using taxonomy ranks when arrow_taxonomy_labels is given", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(
    fit,
    biplot_loadings = TRUE,
    arrow_labels = TRUE,
    arrow_taxonomy_labels = "Phylum"
  ))
  text_layers <- Filter(function(l) inherits(l$geom, "GeomText"), plt$layers)
  expect_true(all(text_layers[[1]]$data$Names %in% fit$tax_table$Phylum))
})

test_that("plot_beta_diversity labels the secondary axis 'Loadings' (raw.loadings distinction dropped)", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(fit, biplot_loadings = TRUE))
  x_scale <- Find(function(s) "x" %in% s$aesthetics, plt$scales$scales)
  expect_match(x_scale$secondary.axis$name, "^Loadings ")
})

test_that("plot_beta_diversity has dropped the dead raw.loadings/raw_loadings argument", {
  expect_false("raw_loadings" %in% names(formals(plot_beta_diversity)))
  expect_false("raw.loadings" %in% names(formals(plot_beta_diversity)))
})

test_that("plot_beta_diversity filters short covariate arrows by arrow_cutoff_covar", {
  ps <- make_bd_ps(n_samples = 14)
  fit <- get_beta_diversity(ps, method = "RDA", model = "group + batch", dist = "bray")
  plt_all <- suppressWarnings(plot_beta_diversity(fit, biplot_covariates = TRUE))
  segment_all <- Filter(function(l) inherits(l$geom, "GeomSegment"), plt_all$layers)
  expect_equal(nrow(segment_all[[1]]$data), 3L)

  plt_cut <- suppressWarnings(plot_beta_diversity(
    fit,
    biplot_covariates = TRUE,
    arrow_cutoff_covar = 0.5
  ))
  segment_cut <- Filter(function(l) inherits(l$geom, "GeomSegment"), plt_cut$layers)
  expect_lt(nrow(segment_cut[[1]]$data), nrow(segment_all[[1]]$data))
})

test_that("plot_beta_diversity maps loading arrows to individual colors when color_arrows_by_taxa = TRUE", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(
    fit,
    biplot_loadings = TRUE,
    color_arrows_by_taxa = TRUE
  ))
  segment_layer <- Filter(function(l) inherits(l$geom, "GeomSegment"), plt$layers)[[1]]
  expect_true("colour" %in% names(segment_layer$mapping))
  expect_null(segment_layer$aes_params$colour)
})

test_that("plot_beta_diversity fixes loading arrows to a single darkgrey color when color_arrows_by_taxa = FALSE", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(fit, biplot_loadings = TRUE))
  segment_layer <- Filter(function(l) inherits(l$geom, "GeomSegment"), plt$layers)[[1]]
  expect_false("colour" %in% names(segment_layer$mapping))
  expect_equal(segment_layer$aes_params$colour, "darkgrey")
})

test_that("plot_beta_diversity omits the secondary axis when both loadings and covariates arrows are drawn", {
  fit <- make_bd_constrained_fit()
  plt <- suppressWarnings(plot_beta_diversity(
    fit,
    biplot_loadings = TRUE,
    biplot_covariates = TRUE
  ))
  x_scale <- Find(function(s) "x" %in% s$aesthetics, plt$scales$scales)
  expect_null(x_scale)
})

test_that("plot_beta_diversity omits the secondary axis when marginal_plot is set", {
  fit <- make_bd_fit()
  # No label_name: the marginal-plot condition (is.factor(label)) isn't met,
  # so ggMarginal() never wraps plt -- lets us inspect $scales directly while
  # still exercising the arrow-plot's own is.null(marginal_plot) check.
  plt <- suppressWarnings(plot_beta_diversity(
    fit,
    biplot_loadings = TRUE,
    marginal_plot = "boxplot"
  ))
  expect_s3_class(plt, "ggplot")
  x_scale <- Find(function(s) "x" %in% s$aesthetics, plt$scales$scales)
  expect_null(x_scale)
})

test_that("plot_beta_diversity wraps the plot in ggMarginal when marginal_plot is set, no facet, no grid", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(
    fit,
    label_name = "facet1",
    marginal_plot = "boxplot"
  ))
  expect_s3_class(plt, "ggExtraPlot")
})

test_that("plot_beta_diversity adds the animation variable as a plot_df column", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(fit, animation_variable_name = "facet1"))
  expect_true("facet1" %in% names(plt$layers[[1]]$data))
})

test_that("plot_beta_diversity's title includes the method and, for constrained fits, the model formula", {
  fit <- make_bd_constrained_fit()
  plt <- suppressWarnings(plot_beta_diversity(fit))
  expect_match(plt$labels$title, "RDA", fixed = TRUE)
  expect_match(plt$labels$title, "group", fixed = TRUE)
})

test_that("plot_beta_diversity's subtitle reports taxrank as 'none' when NULL", {
  fit <- make_bd_fit()
  plt <- suppressWarnings(plot_beta_diversity(fit))
  expect_match(plt$labels$subtitle, "taxa agglom: none", fixed = TRUE)
})

# ---- animate_by_variable() ----

test_that("animate_by_variable returns a gganim/gif_image object by default", {
  ps <- make_bd_ps()
  plt <- suppressWarnings(beta_diversity(
    ps,
    group = "group",
    animation_variable_name = "group"
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
    animation_variable_name = "group"
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
    animation_variable_name = "group"
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
