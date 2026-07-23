#' Fill in NA coordinates for samples a model can't project, with a warning
#'
#' Shared tail of the "project fit-subset scores onto all samples" pattern
#' used by every ordination method in [get_beta_diversity()] that cannot
#' project new samples onto an existing fit (tSNE, UMAP, NMDS, DCA, dbRDA).
#' Builds an all-sample matrix of `NA`, fills in the fit-subset rows, and
#' warns that the rest are unprojected.
#'
#' @param method_label Method name to name in the warning message.
#' @param all_sample_names Character vector of every sample name, used as the
#'   row names (and row order) of the returned matrix.
#' @param fit_sample_names Character vector of the fit-subset sample names.
#' @param fit_scores Numeric matrix of scores for the fit-subset samples only
#'   (rows named by `fit_sample_names`); its `colnames` are carried over.
#' @return A numeric matrix with one row per `all_sample_names`, `NA` outside
#'   `fit_sample_names`.
#' @noRd
.warn_and_na_fill <- function(
  method_label,
  all_sample_names,
  fit_sample_names,
  fit_scores
) {
  warning(paste0(
    "Projection of non-fit samples is not supported for ",
    method_label,
    "; coordinates will be NA."
  ))
  all_scores <- matrix(
    NA,
    nrow = length(all_sample_names),
    ncol = ncol(fit_scores),
    dimnames = list(all_sample_names, colnames(fit_scores))
  )
  all_scores[fit_sample_names, ] <- fit_scores
  all_scores
}

#' Fraction of non-zero cells in an `otu_table`
#'
#' Cheap for a [sparse_otu_table-class]: reads the backing `dgCMatrix`'s `@x`
#' length directly, no densification. For a plain (dense) `otu_table`, counts
#' non-zero cells in the matrix the caller already has in hand.
#'
#' Used to decide whether the sparse-matrix euclidean kernel
#' (`euclidean_sparse()`) is actually worth using: a table with a
#' zero-destroying transform applied upstream (e.g. CLR, which replaces true
#' zeros with a small pseudocount) has effectively no sparsity left to
#' exploit, and representing it as a `dgCMatrix` just adds sparse-matrix
#' bookkeeping overhead on top of what would otherwise be a plain dense BLAS
#' `crossprod()` -- benchmarked to be ~40-70% slower than the dense path once
#' density exceeds ~0.7 at realistic sample/taxa counts (100-300 samples,
#' 2000-5000 taxa).
#'
#' @param ot An `otu_table` (sparse or dense).
#' @return Numeric in `[0, 1]`.
#' @noRd
.otu_table_density <- function(ot) {
  if (is(ot, "sparse_otu_table")) {
    length(ot@sparse_data@x) / prod(dim(ot@sparse_data))
  } else {
    mat <- as(ot, "matrix")
    sum(mat != 0, na.rm = TRUE) / length(mat)
  }
}

# Below this density, the sparse euclidean kernel (dgCMatrix crossprod) beats
# the dense decostand()+vegdist() reference with a comfortable margin (>=2x
# at both benchmarked sizes); above it, dense wins and can be meaningfully
# faster once density approaches 1 (e.g. after a CLR transform).
.EUCLIDEAN_SPARSE_DENSITY_THRESHOLD <- 0.5

#' Standardized (z-scored per taxon) euclidean distance, sparse when it's
#' actually worth it
#'
#' Shared by every place in [get_beta_diversity()] where `"euclidean"` is
#' offered as a distance metric (the unconstrained PCoA/tSNE/UMAP branch and
#' dbRDA) -- unlike raw [sparse_distance()], euclidean is always standardized
#' here so that high-abundance taxa don't dominate the distance in
#' ordination, matching `PCA`'s and `RDA`'s hardcoded `scale = TRUE`
#' (`.get_beta_diversity_from_abundance()`'s `pca` branch and
#' `.get_beta_diversity_constrained()`'s `rda` branch): `PCA`/`RDA` and their
#' distance-based counterparts `PCoA`/`dbRDA` + Euclidean are meant to agree
#' exactly (up to axis sign).
#'
#' Routes through the sparse kernel (`euclidean_sparse()`) when
#' `.otu_table_density()` says it's worth it, else falls back to the dense
#' `decostand("standardize")` + `vegdist("euclidean")` reference the sparse
#' kernel was verified against.
#'
#' @param physeq A `phyloseq` object (already reversed to `taxa_are_rows =
#'   FALSE` by the caller).
#' @return A [stats::dist] object.
#' @noRd
.standardized_euclidean <- function(physeq) {
  if (.otu_table_density(otu_table(physeq)) <= .EUCLIDEAN_SPARSE_DENSITY_THRESHOLD) {
    euclidean_sparse(physeq, standardize = TRUE)
  } else {
    otu_mat <- as(otu_table(physeq), "matrix")
    vegdist(
      decostand(
        otu_mat[, apply(otu_mat, 2, var, na.rm = TRUE) > 0],
        method = "standardize"
      ),
      method = "euclidean"
    )
  }
}

#' Distance-based ordination methods (PCoA, tSNE, UMAP)
#'
#' Handles the `get_beta_diversity()` methods that first reduce `physeq` to a
#' sample x sample distance matrix, then ordinate that matrix. Distance
#' matrices are the natural input for [sparse_distance()], so this is the one
#' branch that can stay sparse until the ordination step itself.
#'
#' @inheritParams get_beta_diversity
#' @param physeq_fit `physeq` pruned to the fit-subset samples (or identical
#'   to `physeq` if no `fit_filter_name`/`fit_filter_values` was given).
#' @param fit_sample_names Character vector of fit-subset sample names, or
#'   `NULL` if all samples are used to fit.
#' @return A list with elements `fit`, `coords`, `loadings`, `covariates`
#'   (always `list()` here — these are unconstrained methods), `eigen_values`.
#' @noRd
.get_beta_diversity_from_distance <- function(
  physeq,
  physeq_fit,
  dist,
  method,
  ndims,
  fit_sample_names,
  pca,
  perplexity,
  nb_neighbors,
  min_dist
) {
  # Compute distance matrix for ALL samples first (needed for projection)
  # NOTE: tSNE and UMAP can take as well a precomputed distance metric
  if (is(otu_table(physeq), "incomplete_otu_table")) {
    svd_emb <- otu_table(physeq)@svd_fit
    emb <- svd_emb$u %*% diag(svd_emb$d)
    rownames(emb) <- sample_names(physeq)
    dist_matrix <- as.matrix(dist(emb))
  } else if (is(dist, "dist")) {
    # if distance matrix furnished directly (not recommended)
    dist_matrix <- dist
  } else {
    if (
      is.null(access(physeq, "phy_tree")) && dist %in% c("unifrac", "wunifrac")
    ) {
      message(
        "No phylogenetic tree in this phyloseq object, bray-curtis distance selected."
      )
      dist_matrix <- phyloseq::distance(physeq, method = "bray")
    } else {
      # NOTE: euclidean and manhattan use z-score standardization (decostand "standardize")
      # so that high-abundance taxa don't dominate the distance in ordination.
      # This produces STANDARDIZED euclidean/manhattan, intentionally different from
      # sparse_distance(ps, "euclidean"/"manhattan") which returns raw distances.
      # The zero-variance filter (apply var > 0) prevents division-by-zero during
      # z-score scaling; if sparse data leaves very few taxa, results may be degenerate.
      # NOTE: if transform_abundances was applied above (e.g. CLR, Hellinger), this
      # z-score step stacks on top of it. For CLR that is still useful (removes
      # scale differences between taxa); for Hellinger it is redundant but harmless.
      # Euclidean goes through .standardized_euclidean() (sparse when worth it, see
      # .otu_table_density()); manhattan has no equivalent sparse identity and always
      # takes the dense path below.
      if (!is.null(dist) && dist == "euclidean") {
        dist_matrix <- .standardized_euclidean(physeq)
      } else if (!is.null(dist) && dist == "manhattan") {
        otu_mat <- as(otu_table(physeq), "matrix")
        dist_matrix <- vegdist(
          decostand(
            otu_mat[, apply(otu_mat, 2, var, na.rm = TRUE) > 0],
            method = "standardize"
          ),
          method = "manhattan"
        )
      } else {
        dist_matrix <- sparse_distance(physeq, method = dist)
      }
    }
  }
  dist_matrix <- as.matrix(dist_matrix)

  # Subset distance matrix to fitting samples only
  dist_matrix_fit <- if (!is.null(fit_sample_names)) {
    dist_matrix[fit_sample_names, fit_sample_names]
  } else {
    dist_matrix
  }

  coords <- list()
  loadings <- list()
  eigen_values <- NULL

  if (method == "pcoa") {
    fit <- vegan::wcmdscale(dist_matrix_fit, eig = TRUE)

    eigen_values <- vegan::eigenvals(fit)
    n_axes <- min(ndims, length(eigen_values), ncol(fit$points))

    # Project all samples using Gower's formula; falls back to scores() when no filter.
    # Implemented directly (instead of predict.wcmdscale) to avoid vegan version issues
    # and because we already have the distance matrices in scope.
    # Derivation: score_i = B[i,] %*% V / sqrt(lambda)
    #           = -0.5 * (d^2_i - delta_plus) %*% fit$points / lambda
    # where delta_plus = colMeans(D^2_fit). Row-centering terms vanish when projected
    # onto eigenvectors of the centered B matrix (they are orthogonal to the
    # all-ones vector).
    site_scores <- if (!is.null(fit_sample_names)) {
      d_fit_sq <- dist_matrix_fit^2
      delta <- colMeans(d_fit_sq)
      d_all <- dist_matrix[, fit_sample_names, drop = FALSE]
      f <- -0.5 * sweep(d_all^2, 2, delta, "-")
      eig_k <- fit$eig[1:n_axes]
      sc <- sweep(f %*% fit$points[, 1:n_axes, drop = FALSE], 2, eig_k, "/")
      rownames(sc) <- rownames(dist_matrix)
      colnames(sc) <- colnames(fit$points)[1:n_axes]
      sc
    } else {
      scores(fit, display = "sites", choices = 1:n_axes, scaling = 1)
    }
    coords[[1]] <- site_scores
    coords[[2]] <- site_scores # PCoA site scores are identical across scalings

    # NOTE: we don't get loadings from the model but we can nevertheless
    #       correlate abundances against axes - be careful when interpreting.
    # We used to get these via envfit(), which fits a multiple regression of
    # each taxon onto ALL requested axes. PCoA axes are orthogonal though, so
    # that multiple regression reduces to a plain per-axis correlation -
    # envfit's QR-based implementation doesn't exploit this and becomes very
    # slow as the number of axes/taxa grows, so just correlate directly.
    if (is(otu_table(physeq), "incomplete_otu_table")) {
      svd_ld <- otu_table(physeq)@svd_fit
      U <- svd_ld$u
      V <- svd_ld$v
      S <- coords[[1]]
      B <- solve(crossprod(U), crossprod(U, S))
      ld <- V %*% B
      rownames(ld) <- taxa_names(physeq)
      colnames(ld) <- colnames(S)
      loadings[[1]] <- ld
    } else {
      loadings[[1]] <- suppressWarnings(cor(
        as(otu_table(physeq), "matrix"),
        coords[[1]]
      ))
    }
    # TODO:for now same thing for both "scalings"
    # loadings[[2]] = wascores(fit$points, otu_table(physeq)) ≈ scaling 2 ?
  } else if (method == "tsne") {
    max_perplexity <- floor((nrow(dist_matrix_fit) - 1) / 3)

    if (is.null(perplexity)) {
      perplexity <- 30
    }

    if (perplexity > max_perplexity) {
      warning(
        "Perplexity must not exceed 3 * perplexity < nrow(X) - 1, set to this value"
      )
      perplexity <- max_perplexity
    }

    fit <- Rtsne(
      dist_matrix_fit,
      pca = pca,
      perplexity = perplexity,
      is_distance = TRUE
    )
    site_scores <- if (!is.null(fit_sample_names)) {
      .warn_and_na_fill(
        "tSNE",
        rownames(dist_matrix),
        fit_sample_names,
        fit$Y
      )
    } else {
      fit$Y
    }
    coords[[1]] <- site_scores
    coords[[2]] <- site_scores
  } else if (method == "umap") {
    umap_config <- umap::umap.defaults
    if (!is.null(nb_neighbors)) {
      umap_config$n_neighbors <- nb_neighbors
    }

    if (!is.null(min_dist)) {
      umap_config$min_dist <- min_dist
    }

    umap_config$input <- "dist" # distance matrix as input

    fit <- umap::umap(dist_matrix_fit, config = umap_config)

    site_scores <- if (!is.null(fit_sample_names)) {
      .warn_and_na_fill(
        "UMAP",
        rownames(dist_matrix),
        fit_sample_names,
        fit$layout
      )
    } else {
      fit$layout
    }
    coords[[1]] <- site_scores
    coords[[2]] <- site_scores
  }

  list(
    fit = fit,
    coords = coords,
    loadings = loadings,
    covariates = list(),
    eigen_values = eigen_values
  )
}

#' Unconstrained ordination directly on the abundance matrix (NMDS/PCA/CA/DCA)
#'
#' Handles the `get_beta_diversity()` methods that hand the (dense) abundance
#' matrix straight to a `vegan` ordination, without a separate
#' distance-matrix step.
#'
#' @inheritParams .get_beta_diversity_from_distance
#' @noRd
.get_beta_diversity_from_abundance <- function(
  physeq,
  physeq_fit,
  dist,
  method,
  ndims,
  fit_sample_names
) {
  otu_fit <- as(otu_table(reverseASV(physeq_fit)), "matrix")

  if (method == "nmds") {
    fit <- vegan::metaMDS(otu_fit, distance = dist, trace = 0)
  } else if (method == "pca") {
    # RDA without constraint = PCA
    fit <- vegan::rda(otu_fit, scale = TRUE)
  } else if (method == "ca") {
    fit <- vegan::cca(otu_fit)
  } else if (method == "dca") {
    fit <- vegan::decorana(otu_fit)
  }
  # NMDS has no eigenvalues; its dimension count comes from fit$points directly
  if (method == "nmds") {
    eigen_values <- NULL
    n_axes <- ncol(fit$points)
  } else {
    eigen_values <- vegan::eigenvals(fit)
    n_axes <- min(ndims, length(eigen_values))
  }

  # Project all samples onto the fitted model if fit.filter is active.
  # RDA (PCA) and CCA support predict(); NMDS and DCA do not.
  get_site_coords <- function(scaling) {
    fit_scores <- scores(
      fit,
      display = "sites",
      choices = 1:n_axes,
      scaling = scaling
    )
    if (is.null(fit_sample_names)) {
      return(fit_scores)
    }
    if (inherits(fit, c("rda", "cca"))) {
      otu_all <- as(otu_table(reverseASV(physeq)), "matrix")
      wa <- predict(fit, newdata = otu_all, type = "wa")
      return(wa[, 1:min(n_axes, ncol(wa)), drop = FALSE])
    }
    .warn_and_na_fill(method, sample_names(physeq), fit_sample_names, fit_scores)
  }

  coords <- list(get_site_coords(1), get_site_coords(2))

  # NMDS/DCA may have no species scores; vegan 2.7.5 returns list() instead of
  # NULL/error when scores are unavailable — validate result is a matrix.
  valid_scores <- function(sc) {
    if (is.matrix(sc) || is.data.frame(sc)) sc else NULL
  }
  loadings <- list(
    valid_scores(tryCatch(
      scores(fit, display = "species", choices = 1:n_axes, scaling = 1),
      error = function(e) NULL
    )),
    valid_scores(tryCatch(
      scores(fit, display = "species", choices = 1:n_axes, scaling = 2),
      error = function(e) NULL
    ))
  )

  list(
    fit = fit,
    coords = coords,
    loadings = loadings,
    covariates = list(),
    eigen_values = eigen_values
  )
}

#' Compute a `newdata`-based prediction only for complete-case rows,
#' `NA`-filling (with a warning) any row missing a needed variable
#'
#' `predict(fit, newdata, type = "lc"/"working", ...)` hard-errors for the
#' *whole* call (`na.fail.default(...)`, from `vegan`'s internal
#' `ordiParseFormula()`) if even one row of `newdata` has an `NA` in a model
#' variable, and `model.matrix()` (used by `.pcca_fitted_all()`) silently
#' *drops* incomplete rows instead -- causing a row-count mismatch a few
#' lines later when the result is combined with something that still has
#' every row. Both defeat "project every sample that can be projected, leave
#' the rest `NA`" (the same treatment already given to samples a *method*
#' can't project, e.g. `.warn_and_na_fill()`). This restricts `f_complete`
#' to the complete-case rows only, then reinserts `NA` rows for the rest, so
#' ordinary matrix arithmetic downstream (subtracting/projecting) propagates
#' `NA` to just those samples' final coordinates instead of the whole
#' pipeline erroring or misaligning.
#'
#' @param newdata Data frame to check for completeness and pass to
#'   `f_complete` (one row per sample).
#' @param relevant_cols Character vector of `newdata` column names that
#'   actually matter for completeness (the model/confounder variables --
#'   `NA` in an unrelated column doesn't exclude a sample).
#' @param f_complete Function taking the complete-case-only `newdata` subset
#'   and returning a matrix with as many rows, in the same row order.
#' @return Matrix with one row per `nrow(newdata)` (row names from
#'   `rownames(newdata)`), `NA` for incomplete rows.
#' @noRd
.predictor_complete_cases <- function(newdata, relevant_cols, f_complete) {
  usable <- stats::complete.cases(newdata[, relevant_cols, drop = FALSE])
  if (all(usable)) {
    return(f_complete(newdata))
  }
  warning(paste0(
    sum(!usable),
    " sample(s) excluded from projection: missing value(s) in a model or ",
    "confounder variable."
  ))
  out_complete <- f_complete(newdata[usable, , drop = FALSE])
  out <- matrix(
    NA_real_,
    nrow = nrow(newdata),
    ncol = ncol(out_complete),
    dimnames = list(rownames(newdata), colnames(out_complete))
  )
  out[usable, ] <- out_complete
  out
}

#' Fitted values for held-out samples from a partial (`Condition()`) model's
#' `pCCA` component, computed by hand
#'
#' Neither `predict()`'s `"lc"` nor `"working"` type accepts `model = "pCCA"`
#' (only `"CCA"`/`"CA"` are valid) -- there is no `vegan`-native way to get a
#' confounders-only model's fitted values for new predictor data (and
#' `"lc"`/model `"CCA"` doesn't apply at all when there's no free predictor:
#' `predict()` silently falls back to `model = "CA"` when `object$CCA` is
#' `NULL`, and `"lc"` scores are explicitly unavailable for that model).
#' This derives the fitted values directly instead: `fit$pCCA$QR` is the
#' very QR decomposition `vegan` fit the confounders against internally
#' (`vegan:::ordPartial()`), so its own regression coefficients
#' (`qr.coef()`), applied to a freshly built confounder design matrix
#' centered the same way `vegan` centered its own (`fit$pCCA$envcentre`),
#' reproduce `qr.fitted(fit$pCCA$QR, fit$Ybar)` exactly for the training rows
#' (verified numerically, including with a confounder level that `vegan`
#' itself dropped as redundant -- the `[, names(fit$pCCA$envcentre)]`
#' reindex below is what keeps that case aligned) and extend cleanly to any
#' row count.
#'
#' @param fit An `rda`- or `capscale`-class object with a non-`NULL` `$pCCA`
#'   (i.e. `confounders` was supplied to [get_beta_diversity()]).
#' @param confounders Character vector of confounder column names, as
#'   originally passed to [get_beta_diversity()].
#' @param df_all `sample_data(physeq)` as a data frame, all samples.
#' @return Numeric matrix, all samples x `ncol(fit$Ybar)`, in the same
#'   standardized units as `fit$Ybar`.
#' @noRd
.pcca_fitted_all <- function(fit, confounders, df_all) {
  z_formula <- as.formula(paste(
    "~",
    paste(paste0("`", confounders, "`"), collapse = " + ")
  ))
  z_all <- model.matrix(z_formula, data = df_all)[, -1, drop = FALSE]
  z_all <- z_all[, names(fit$pCCA$envcentre), drop = FALSE]
  z_all <- sweep(z_all, 2, fit$pCCA$envcentre, "-")
  z_all %*% qr.coef(fit$pCCA$QR, fit$Ybar)
}

#' Fitted values for held-out samples from an RDA/dbRDA fit's free-predictor
#' (`CCA`) component, computed by hand
#'
#' `predict(fit, newdata, type = "working"/"lc", model = "CCA")` hard-errors
#' (`"factor ... has new levels ..."`, from `model.frame()`'s own
#' `xlev`-checking against `fit$terms`) whenever `newdata` has a factor level
#' -- in *either* the free predictor or a `Condition()` confounder -- that
#' wasn't present in the fit-subset data. This happens routinely with
#' `fit_filter`: fitting on a data subset (e.g. two groups out of many) and
#' projecting the rest is exactly "newdata has levels the fit never saw".
#' `.pcca_fitted_all()` above never hits this because it builds
#' `model.matrix()` fresh from `newdata` instead of routing through
#' `predict()`'s stored `xlevels`; this does the same for the free-predictor
#' component. `fit$CCA$QR` turns out to be the QR decomposition of the
#' *combined* `[Condition(...) | model]` design (confirmed via its own
#' column names), but applying only the coefficient rows belonging to the
#' free predictor's own columns to the free predictor's own (freshly built,
#' `newdata`-derived) design matrix -- never touching the confounder columns
#' at all -- reproduces `predict(..., type = "working"/"lc", model = "CCA")`
#' exactly (verified to float tolerance, including with `Condition()`
#' confounders present in the fit). A `newdata` level never seen while
#' fitting simply has no dummy column in `fit$CCA$envcentre` to align to, so
#' it's silently dropped by the reindex below -- equivalent to treating that
#' sample as the reference level for this predictor, the same convention
#' `.pcca_fitted_all()` already established for unseen confounder levels.
#'
#' @param fit An `rda`- or `capscale`-class object with a non-`NULL` `$CCA`
#'   (i.e. a free predictor, not just confounders, was fit).
#' @param model String with the free predictor's own model formula RHS, as
#'   originally passed to [get_beta_diversity()] -- *not* the
#'   `Condition(...)`-augmented version used to build `fit`'s own formula.
#' @param df_all `sample_data(physeq)` as a data frame, all samples.
#' @return Numeric matrix, all samples x `ncol(fit$Ybar)`, in the same
#'   standardized units as `fit$Ybar`.
#' @noRd
.cca_working_fitted <- function(fit, model, df_all) {
  x_formula <- as.formula(paste("~", model))
  x_all <- model.matrix(x_formula, data = df_all)[, -1, drop = FALSE]
  x_all <- x_all[, names(fit$CCA$envcentre), drop = FALSE]
  x_all <- sweep(x_all, 2, fit$CCA$envcentre, "-")
  coefs <- qr.coef(fit$CCA$QR, fit$Ybar)
  x_all %*% coefs[names(fit$CCA$envcentre), , drop = FALSE]
}

#' Project new samples onto the residual (unconstrained `PC`) axes of an RDA
#' fit
#'
#' `vegan`'s own `predict(fit, newdata, type = "wa", model = "CA")` does not
#' residualize `newdata` against the fitted constrained (`CCA`) component
#' before projecting it onto the residual species scores -- verified to
#' diverge from `scores(fit, ...)`'s residual-axis output by 20-70% (and
#' sometimes in sign), even when `newdata` is exactly the fit-subset data the
#' model was trained on. This reimplements the projection correctly: since
#' RDA's residual axes are a plain (unweighted) PCA of the residuals, the
#' constrained fit's prediction (`.cca_working_fitted()`, reproducing
#' `predict(..., type = "working", model = "CCA")` but tolerant of `newdata`
#' levels the fit never saw, in the same standardized units as `fit$Ybar`)
#' can be subtracted from `otu_all` directly, and the remainder projected
#' onto `fit$CA$v` --
#' reproduces `scores(fit, display = "sites", scaling = scaling)`'s residual
#' columns exactly (to float tolerance) for the fit subset, and extends that
#' correctly to out-of-fit samples. CCA's residual axes are chi-square
#' weighted rather than a plain PCA, so this does not carry over to `method =
#' "cca"` -- see the caller.
#'
#' Also handles a confounders-only fit (no free predictor: `fit$CCA` is
#' `NULL`, only `fit$pCCA`) -- there, every requested axis is "residual"
#' (there's no constrained part to compute at all), and the quantity
#' subtracted before projecting onto `fit$CA$v` is the confounders' own
#' fitted values (`.pcca_fitted_all()`) instead of the free predictor's.
#'
#' @param fit An `rda`-class object from `vegan::rda()`.
#' @param otu_all Numeric matrix, all samples (fit + out-of-fit) x taxa,
#'   samples-are-rows.
#' @param df_all `sample_data(physeq)` as a data frame, all samples.
#' @param n_axes Number of residual axes needed.
#' @param scaling `1` or `2` (vegan scaling convention).
#' @param model String with the free predictor's own model formula RHS
#'   (*not* the `Condition(...)`-augmented version used to build `fit`);
#'   only used (and required) when `fit` has a free predictor (`fit$CCA` is
#'   not `NULL`) -- passed through to `.cca_working_fitted()`.
#' @param confounders Character vector of confounder column names; only used
#'   (and required) when `fit` has no free predictor (`fit$CCA` is `NULL`).
#' @param predictor_vars Character vector of every model/confounder column
#'   name `df_all` needs to be complete in -- passed through to
#'   `.predictor_complete_cases()`, which `NA`-fills (rather than errors on)
#'   any sample missing one of them.
#' @return Numeric matrix, all samples x `n_axes`, columns named `PC1`,
#'   `PC2`, ...
#' @noRd
.rda_residual_scores <- function(
  fit,
  otu_all,
  df_all,
  n_axes,
  scaling,
  model = NULL,
  confounders = NULL,
  predictor_vars
) {
  cent <- attr(fit$Ybar, "scaled:center")
  scal <- attr(fit$Ybar, "scaled:scale")
  otu_all <- otu_all[, names(cent), drop = FALSE]
  nr <- nobs(fit) - 1

  ybar_all <- sweep(otu_all, 2, cent, "-") / sqrt(nr)
  if (!is.null(scal)) {
    nz <- scal > 0
    ybar_all[, nz] <- sweep(ybar_all[, nz, drop = FALSE], 2, scal[nz], "/")
  }
  yhat_all <- if (!is.null(fit$CCA) && fit$CCA$rank > 0) {
    .predictor_complete_cases(df_all, predictor_vars, function(d) {
      .cca_working_fitted(fit, model, d)
    })
  } else {
    .predictor_complete_cases(df_all, predictor_vars, function(d) {
      .pcca_fitted_all(fit, confounders, d)
    })
  }
  resid_all <- ybar_all - yhat_all

  v_ca <- fit$CA$v[, 1:n_axes, drop = FALSE]
  eig_ca <- fit$CA$eig[1:n_axes]
  scores_out <- resid_all %*% v_ca %*% diag(1 / sqrt(eig_ca), nrow = n_axes)
  colnames(scores_out) <- paste0("PC", seq_len(n_axes))

  if (scaling) {
    slam <- sqrt(eig_ca / fit$tot.chi)
    const <- sqrt(sqrt((nobs(fit) - 1) * fit$tot.chi))
    lam <- list(slam, 1, sqrt(slam))[[abs(scaling)]]
    scores_out <- const * sweep(scores_out, 2, lam, "*")
  }
  scores_out
}

#' Project new samples onto a dbRDA (`capscale`) fit's constrained and
#' residual axes
#'
#' `vegan`'s own `predict(fit, newdata, type = "wa")` is unconditionally
#' blocked for `capscale` objects (`stop("'wa' scores not available in
#' capscale with 'newdata'")`), and the alternative it does allow, `type =
#' "lc"`, is not a usable substitute: it's purely the model's fitted value
#' from the predictors, so every sample sharing the same predictor
#' combination collapses to one identical point, discarding all of that
#' sample's actual community signal (verified: with a single 2-level factor
#' predictor, all samples projected to just 2 distinct points). This
#' reimplements the projection from scratch by treating `capscale` as what
#' it actually is internally (see `vegan:::ordConstrained()`'s `"capscale"`
#' branch and `vegan:::initCAP()`): classical scaling (PCoA, via
#' `vegan::wcmdscale()`) of the training distance matrix, mean-centered but
#' *not* z-scored (`initCAP()`, unlike `rda()`'s `initPCA()`), followed by a
#' plain unweighted RDA of that embedding against the predictors. Held-out
#' samples are embedded into the exact same classical-scaling space via
#' Gower's interpolation formula (the same one
#' `.get_beta_diversity_from_distance()` already uses for standalone PCoA
#' projection), replicating `capscale()`'s own pre-embedding "adjust"
#' rescaling (`vegan::capscale`'s `if (max(X) >= 4) X <- X / sqrt(k)` step)
#' so the interpolation lands in the same units. From there, projection onto
#' the constrained and residual axes follows the same "wa" logic as
#' `.rda_residual_scores()` -- except `capscale`'s own species/embedding-axis
#' scores (`fit$CCA$v`/`fit$CA$v`) are only populated when a `comm =`
#' argument was supplied to `capscale()` (this package never supplies one,
#' since there's no single meaningful "community matrix" once the ordination
#' is already distance-based), so `v` is reconstructed by hand from `u`,
#' `eig`, and the model's own QR decomposition (`fit$CCA$QR`) instead --
#' verified to exactly reproduce `vegan`'s own `v` when it IS available
#' (i.e. cross-checked against a plain `rda()` fit). The whole pipeline
#' reproduces `scores(fit, display = "sites", scaling = scaling)` exactly
#' (to float tolerance) for the fit subset.
#'
#' @param fit A `capscale`-class object from `vegan::capscale()`.
#' @param dist_matrix_fit Sample x sample distance matrix (plain matrix, not
#'   `dist`) for the fit-subset samples only, using the exact distance `fit`
#'   itself was built from.
#' @param dist_matrix_all Same distance metric and (for `euclidean`)
#'   standardization reference as `dist_matrix_fit`, but covering every
#'   sample (fit + out-of-fit); only the `[, fit_sample_names]` cross-block
#'   is actually used, i.e. distances from every sample to each fit-subset
#'   sample.
#' @param df_all `sample_data(physeq)` as a data frame, all samples.
#' @param n_axes Number of axes requested; split automatically between
#'   constrained and residual based on `fit`'s own rank.
#' @param scaling `1` or `2` (vegan scaling convention).
#' @param model String with the free predictor's own model formula RHS
#'   (*not* the `Condition(...)`-augmented version used to build `fit`);
#'   only used (and required) when `fit` has a free predictor (`fit$CCA` is
#'   not `NULL`) -- passed through to `.cca_working_fitted()`, which this
#'   function's `has_free_predictor` branch uses in place of `predict(type =
#'   "lc")` (unusable as-is: it hard-errors on any `newdata` factor level
#'   the fit-subset data never saw, which `fit_filter` runs into routinely).
#' @param confounders Character vector of confounder column names; only used
#'   (and required) when `fit` has no free predictor (`fit$CCA` is `NULL`) --
#'   see `.pcca_fitted_all()`, the `capscale` analogue of which this
#'   function's confounders-only branch uses in place of `predict(type =
#'   "lc")`.
#' @param predictor_vars Character vector of every model/confounder column
#'   name `df_all` needs to be complete in -- passed through to
#'   `.predictor_complete_cases()`, which `NA`-fills (rather than errors on)
#'   any sample missing one of them.
#' @return Numeric matrix, all samples (`rownames(dist_matrix_all)`) x up to
#'   `n_axes` columns, named like `scores(fit, display = "sites")` (`CAP1`,
#'   `CAP2`, ..., then `MDS1`, `MDS2`, ...) -- or just `MDS1`, `MDS2`, ... for
#'   a confounders-only fit, which has no constrained axis at all.
#' @noRd
.dbrda_projected_scores <- function(
  fit,
  dist_matrix_fit,
  dist_matrix_all,
  df_all,
  n_axes,
  scaling,
  model = NULL,
  confounders = NULL,
  predictor_vars
) {
  fit_names <- rownames(dist_matrix_fit)
  k <- length(fit_names) - 1
  adjust <- if (max(dist_matrix_fit) >= 4 + sqrt(.Machine$double.eps)) {
    sqrt(k)
  } else {
    1
  }
  d_fit_adj <- as.dist(dist_matrix_fit / adjust)
  d_all_adj <- dist_matrix_all / adjust

  emb <- wcmdscale(d_fit_adj, eig = TRUE, x.ret = TRUE, add = FALSE)
  n_dims <- ncol(emb$points)

  # Gower's classical-scaling interpolation formula (see the PCoA branch of
  # .get_beta_diversity_from_distance() for the same derivation, spelled out
  # in full there).
  delta <- colMeans(as.matrix(d_fit_adj)^2)
  cross <- d_all_adj[, fit_names, drop = FALSE]
  f <- -0.5 * sweep(cross^2, 2, delta, "-")
  z_all <- sweep(f %*% emb$points, 2, emb$eig[1:n_dims], "/")
  rownames(z_all) <- rownames(dist_matrix_all)

  cent <- attr(fit$Ybar, "scaled:center")
  ybar_all <- sweep(z_all, 2, cent, "-")

  has_free_predictor <- !is.null(fit$CCA) && fit$CCA$rank > 0

  if (has_free_predictor) {
    take_cca <- fit$CCA$rank
    n_cca_take <- min(n_axes, take_cca)
    n_ca_take <- min(max(0, n_axes - take_cca), fit$CA$rank)

    u_cca <- fit$CCA$u[, 1:take_cca, drop = FALSE]
    eig_cca <- fit$CCA$eig[1:take_cca]
    yhat_fit <- qr.fitted(fit$CCA$QR, fit$Ybar)
    v_cca <- t(yhat_fit) %*% u_cca %*% diag(1 / sqrt(eig_cca), nrow = take_cca)

    u_lc_all <- .predictor_complete_cases(df_all, predictor_vars, function(d) {
      yhat_d <- .cca_working_fitted(fit, model, d)
      yhat_d %*% v_cca %*% diag(1 / sqrt(eig_cca), nrow = take_cca)
    })
    yhat_all <- u_lc_all %*% diag(sqrt(eig_cca), nrow = take_cca) %*% t(v_cca)

    wa_cca <- sweep(ybar_all %*% v_cca, 2, sqrt(eig_cca), "/")
    wa_cca <- wa_cca[, 1:n_cca_take, drop = FALSE]
    colnames(wa_cca) <- paste0("CAP", seq_len(n_cca_take))
    eig_used <- eig_cca[1:n_cca_take]
  } else {
    # Confounders-only fit: every requested axis is "residual" (there's no
    # constrained part), and the quantity subtracted before projecting onto
    # fit$CA$v is the confounders' own fitted values instead of the free
    # predictor's -- see .pcca_fitted_all().
    n_ca_take <- min(n_axes, fit$CA$rank)
    yhat_fit <- qr.fitted(fit$pCCA$QR, fit$Ybar)
    yhat_all <- .predictor_complete_cases(df_all, predictor_vars, function(d) {
      .pcca_fitted_all(fit, confounders, d)
    })
    wa_cca <- NULL
    eig_used <- numeric(0)
  }

  if (n_ca_take > 0) {
    u_ca <- fit$CA$u[, 1:n_ca_take, drop = FALSE]
    eig_ca <- fit$CA$eig[1:n_ca_take]
    resid_fit <- fit$Ybar - yhat_fit
    resid_all <- ybar_all - yhat_all
    v_ca <- t(resid_fit) %*% u_ca %*% diag(1 / sqrt(eig_ca), nrow = n_ca_take)
    wa_ca <- sweep(resid_all %*% v_ca, 2, sqrt(eig_ca), "/")
    colnames(wa_ca) <- paste0("MDS", seq_len(n_ca_take))
    result <- if (is.null(wa_cca)) wa_ca else cbind(wa_cca, wa_ca)
    eig_used <- c(eig_used, eig_ca)
  } else {
    result <- wa_cca
  }

  if (scaling) {
    tot.chi <- fit$tot.chi
    const <- sqrt(sqrt((nobs(fit) - 1) * tot.chi))
    slam <- sqrt(eig_used / tot.chi)
    lam <- list(slam, 1, sqrt(slam))[[abs(scaling)]]
    result <- const * sweep(result, 2, lam, "*")
  }
  result
}

#' Constrained/supervised ordination (CCA/RDA/dbRDA)
#'
#' Handles the `get_beta_diversity()` methods that fit the abundance matrix
#' against a `sample_data` model formula.
#'
#' @inheritParams .get_beta_diversity_from_distance
#' @param model String with the right-hand side of the model formula (e.g.
#'   `"Var1+Var2"`).
#' @param confounders Character vector of `sample_data` column names to
#'   partial out via `vegan::Condition()`.
#' @noRd
.get_beta_diversity_constrained <- function(
  physeq,
  physeq_fit,
  dist,
  method,
  model,
  confounders,
  ndims,
  fit_sample_names
) {
  # vegan::rda()/cca()/capscale()'s na.action = na.exclude is unreliable for
  # this package's purposes: an NA in a fit-subset sample's model/confounder
  # value can crash the fit itself (verified: even a bare vegan::rda(... ,
  # na.action = na.exclude) call errors internally, in ordiNAexclude()'s own
  # predict() call, trying to pad the excluded row back in). Pruning those
  # samples out of physeq_fit before fitting sidesteps that entirely -- the
  # model just never sees an NA. If fit_sample_names was NULL (no fit_filter
  # requested), pruning here means physeq_fit is now smaller than physeq, so
  # fit_sample_names is set to the (now-pruned) actual fit set purely to
  # make get_constrained_coords() take its "project every physeq sample"
  # branches below instead of assuming physeq_fit == physeq -- that's what
  # NA-fills (with a warning) the pruned-out samples' coordinates instead of
  # silently returning a shorter, misaligned coords matrix.
  predictor_vars <- unique(c(
    if (!is.null(model)) all.vars(stats::as.formula(paste("~", model))),
    confounders
  ))
  if (length(predictor_vars) > 0) {
    fit_sdata <- as(sample_data(physeq_fit), "data.frame")
    fit_usable <- stats::complete.cases(fit_sdata[, predictor_vars, drop = FALSE])
    if (!all(fit_usable)) {
      warning(paste0(
        sum(!fit_usable),
        " fit-subset sample(s) excluded from the model fit: missing ",
        "value(s) in a model or confounder variable."
      ))
      physeq_fit <- prune_samples(fit_usable, physeq_fit)
      if (is.null(fit_sample_names)) {
        fit_sample_names <- sample_names(physeq_fit)
      }
    }
  }

  # .cca_working_fitted() needs the free predictor's own formula RHS, without
  # the Condition(...) term about to be appended below -- Condition() is
  # vegan formula-parsing syntax, not something model.matrix() understands.
  free_predictor_model <- model

  if (!is.null(confounders)) {
    model <- paste0(
      model,
      " + ",
      "Condition(",
      paste(paste0("`", confounders, "`"), collapse = " + "),
      ")"
    )
  }

  df_fit <- as(sample_data(physeq_fit), "data.frame")

  if (method == "cca") {
    formula <- as.formula(paste(
      "as(otu_table(reverseASV(physeq_fit)), 'matrix')",
      "~",
      model
    ))
    fit <- vegan::cca(formula, data = df_fit, na.action = na.exclude)
  } else if (method == "rda") {
    formula <- as.formula(paste(
      "as(otu_table(reverseASV(physeq_fit)), 'matrix')",
      "~",
      model
    ))
    # scale = TRUE (z-score each taxon) matches PCA's hardcoded scale = TRUE
    # above and dbRDA's standardized euclidean below, so RDA and dbRDA +
    # Euclidean agree exactly (up to axis sign) -- see .standardized_euclidean().
    fit <- vegan::rda(formula, data = df_fit, scale = TRUE, na.action = na.exclude)
  } else if (method == "dbrda") {
    # Unlike RDA/CCA, dbRDA is fundamentally distance-based: capscale()
    # accepts a pre-computed dist object as the formula LHS and skips its own
    # (dense) vegdist() call entirely when given one. Euclidean goes through
    # .standardized_euclidean() (matching RDA's scale = TRUE and staying
    # sparse whenever the table actually is, see .otu_table_density()); every
    # other metric goes through sparse_distance() directly, keeping dbRDA off
    # the dense path -- verified numerically identical (site scores,
    # eigenvalues, permutation F-stat) to passing the raw matrix + dist=
    # method name.
    dist_obj <- if (dist == "euclidean") {
      .standardized_euclidean(physeq_fit)
    } else {
      sparse_distance(physeq_fit, method = dist)
    }
    # For projecting out-of-fit samples later (.dbrda_projected_scores()):
    # the same metric extended to every sample. Ordinary dissimilarities
    # (bray, jaccard, unifrac, ...) are pairwise-only, so the fit-subset
    # block of the full matrix is identical to dist_obj regardless; only
    # euclidean needs care, since its z-scoring reference (taxon sd) must
    # come from physeq_fit specifically, not physeq, to stay consistent with
    # dist_obj -- .standardized_euclidean(physeq) would z-score against the
    # wrong (full-dataset) statistics.
    dbrda_dist_all <- if (!is.null(fit_sample_names)) {
      if (dist == "euclidean") {
        otu_fit_mat <- as(otu_table(reverseASV(physeq_fit)), "matrix")
        otu_all_mat <- as(otu_table(reverseASV(physeq)), "matrix")
        keep <- apply(otu_fit_mat, 2, var, na.rm = TRUE) > 0
        taxon_sd <- apply(otu_fit_mat[, keep, drop = FALSE], 2, sd, na.rm = TRUE)
        # Euclidean distance is translation-invariant, so centering (unlike
        # scaling) can't change the result and is skipped -- see
        # euclidean_sparse()'s same observation.
        z_all <- sweep(otu_all_mat[, keep, drop = FALSE], 2, taxon_sd, "/")
        as.matrix(vegdist(z_all, method = "euclidean"))
      } else {
        as.matrix(sparse_distance(physeq, method = dist))
      }
    }
    formula <- as.formula(paste("dist_obj", "~", model))
    fit <- vegan::capscale(formula, data = df_fit, na.action = na.exclude)
  } else {
    stop("Invalid method")
  }

  eigen_values <- vegan::eigenvals(fit)
  n_axes <- min(ndims, length(eigen_values))

  # For RDA/CCA/dbRDA project all samples using WA scores.
  # capscale (dbRDA) inherits from "rda"/"cca" so check method name, not class.
  get_constrained_coords <- function(scaling) {
    # No free (unconditioned) predictor — e.g. confounders-only partial RDA/CCA
    # — leaves fit$CCA NULL, so there are no "wa" scores to predict(); fall
    # through to the scores()-based branch below instead.
    has_free_predictor <- !is.null(fit$CCA) && fit$CCA$rank > 0
    if (!is.null(fit_sample_names) && method %in% c("rda", "cca") && has_free_predictor) {
      otu_all <- as(otu_table(reverseASV(physeq)), "matrix")
      # predict(type="wa") defaults to model="CCA" — constrained axes only.
      # Exact for every sample, fit or not: "wa" scores are, by definition,
      # weighted averages of the raw community data against the model's
      # species scores, so passing `scaling` through here (it was previously
      # dropped, silently forcing scaling="none" regardless of what the
      # caller asked for -- invisible for PCoA, whose site scores don't
      # depend on scaling, but wrong for RDA/CCA, where they do) makes this
      # match scores(fit, ...) exactly for the fit subset.
      wa_cca <- predict(fit, newdata = otu_all, type = "wa", model = "CCA", scaling = scaling)
      n_cca <- ncol(wa_cca)
      wa_cca <- wa_cca[, 1:min(n_axes, n_cca), drop = FALSE]
      n_ca_need <- max(0, n_axes - n_cca)
      if (n_ca_need > 0 && !is.null(fit$CA) && fit$CA$rank > 0) {
        n_ca_take <- min(n_ca_need, fit$CA$rank)
        if (method == "rda") {
          df_all <- as(sample_data(physeq), "data.frame")
          wa_ca <- .rda_residual_scores(
            fit,
            otu_all,
            df_all,
            n_ca_take,
            scaling,
            model = free_predictor_model,
            predictor_vars = predictor_vars
          )
          return(cbind(wa_cca, wa_ca))
        }
        # CCA's residual axes are chi-square weighted (unlike RDA's plain
        # PCA-of-residuals, see .rda_residual_scores()), so there's no
        # equivalent closed-form projection here, and vegan's own
        # predict(type="wa", model="CA") doesn't residualize newdata against
        # the constrained fit first -- verified wrong even for fit-subset
        # samples (diverges from the fit's own scores() by 20-70%). NA-fill +
        # warn instead, matching the dbRDA no-projection path below, but only
        # for these residual columns / the out-of-fit rows: the constrained
        # columns above and the fit subset's own residual scores stay exact.
        ca_scores <- scores(
          fit,
          display = "sites",
          choices = n_cca + seq_len(n_ca_take),
          scaling = scaling
        )
        wa_ca <- .warn_and_na_fill("CCA", sample_names(physeq), fit_sample_names, ca_scores)
        return(cbind(wa_cca, wa_ca))
      }
      return(wa_cca)
    }
    if (
      !is.null(fit_sample_names) &&
        method == "rda" &&
        !has_free_predictor &&
        !is.null(fit$pCCA) &&
        fit$CA$rank > 0
    ) {
      # Confounders-only RDA (no free predictor: "Model" empty, only
      # "Confounders"): every requested axis is a residual axis after
      # partialling out the confounders -- same math as the residual-axis
      # case above, just subtracting the confounders' own fit
      # (.pcca_fitted_all()) instead of the free predictor's.
      otu_all <- as(otu_table(reverseASV(physeq)), "matrix")
      df_all <- as(sample_data(physeq), "data.frame")
      n_ca_take <- min(n_axes, fit$CA$rank)
      return(.rda_residual_scores(
        fit,
        otu_all,
        df_all,
        n_ca_take,
        scaling,
        confounders = confounders,
        predictor_vars = predictor_vars
      ))
    }
    if (
      !is.null(fit_sample_names) &&
        method == "dbrda" &&
        (has_free_predictor || (!is.null(fit$pCCA) && fit$CA$rank > 0))
    ) {
      df_all <- as(sample_data(physeq), "data.frame")
      dist_matrix_fit <- as.matrix(dist_obj)
      return(.dbrda_projected_scores(
        fit,
        dist_matrix_fit,
        dbrda_dist_all,
        df_all,
        n_axes,
        scaling,
        model = free_predictor_model,
        confounders = confounders,
        predictor_vars = predictor_vars
      ))
    }
    if (!is.null(fit_sample_names)) {
      # Reached for: CCA (with or without a free predictor -- its residual
      # axes are chi-square weighted, so there's no equivalent closed-form
      # projection, see .rda_residual_scores()'s docs), or any constrained
      # method/model combination with neither a free predictor nor
      # confounders to project from.
      sc <- scores(fit, display = "sites", choices = 1:n_axes, scaling = scaling)
      method_label <- c(rda = "RDA", cca = "CCA", dbrda = "dbRDA")[[method]]
      return(.warn_and_na_fill(method_label, sample_names(physeq), fit_sample_names, sc))
    }
    scores(fit, display = "sites", choices = 1:n_axes, scaling = scaling)
  }

  # get_constrained_coords() warns (once) via .warn_and_na_fill() on the dbRDA
  # no-projection path; scaling 1/2 genuinely differ (unlike PCoA, so they
  # can't just share one computed result), but the warning reason doesn't --
  # suppress the identical repeat on the second call.
  coords <- list(
    get_constrained_coords(1),
    suppressWarnings(get_constrained_coords(2))
  )
  loadings <- if (method == "dbrda") {
    # capscale() is never given comm= above (there's no single meaningful
    # community matrix once the ordination is already distance-based), so
    # vegan sets species scores to NA unconditionally (its own else-branch:
    # `sol$CA$v[] <- NA; if (!is.null(sol$CCA)) sol$CCA$v[] <- NA`) --
    # scores(fit, display = "species", ...) is just an all-NA placeholder,
    # which .plot_beta_diversity_arrow_group()'s na.omit() then reduces to
    # zero taxa arrows, always (verified unchanged from before this
    # package's dbRDA fit_filter support was added). Falls back to the same
    # per-taxon correlation-against-axes approach the standalone PCoA branch
    # already uses (.get_beta_diversity_from_distance()) for the identical
    # problem, using coords[[1]] so it covers every projected sample, not
    # just the fit subset. use = "complete.obs" (unlike PCoA's plain cor(),
    # which never needs it: PCoA always projects every sample) guards
    # against a handful of NA rows in coords[[1]] -- e.g. from
    # .predictor_complete_cases() -- poisoning every taxon's correlation.
    otu_mat <- as(otu_table(reverseASV(physeq)), "matrix")
    ld <- suppressWarnings(cor(otu_mat, coords[[1]], use = "complete.obs"))
    list(ld, ld) # same caveat as PCoA: no separate formula for scaling 2
  } else {
    list(
      scores(fit, display = "species", choices = 1:n_axes, scaling = 1),
      scores(fit, display = "species", choices = 1:n_axes, scaling = 2)
    )
  }
  covariates <- list(
    scores(fit, display = "bp", choices = 1:n_axes, scaling = 1),
    scores(fit, display = "bp", choices = 1:n_axes, scaling = 2)
  )

  list(
    fit = fit,
    coords = coords,
    loadings = loadings,
    covariates = covariates,
    eigen_values = eigen_values,
    model = model # possibly confounders-augmented; caller reports this, not its own copy
  )
}

#' Compute a Beta-Diversity Ordination
#'
#' Ordinates a `phyloseq` object's samples by inter-sample dissimilarity, via
#' one of three families of `vegan`/`Rtsne`/`umap` methods: distance-based
#' (`"PCoA"`, `"tSNE"`, `"UMAP"`), direct-on-abundance (`"NMDS"`, `"PCA"`,
#' `"CA"`, `"DCA"`), or constrained/supervised on a `sample_data` model
#' formula (`"CCA"`, `"RDA"`, `"dbRDA"`). Optionally fits the model on a
#' subset of samples (`fit_filter_name`/`fit_filter_values`) and projects the
#' rest onto that fit; some methods have no such projection and fall back to
#' `NA` coordinates with a `warning()` for the non-fit samples.
#'
#' @param physeq A `phyloseq` object.
#' @param taxrank Taxonomic rank to agglomerate to before ordination (e.g.
#'   `"Phylum"`), via [tax_glom()]. Default `NULL` (no agglomeration).
#' @param fraction_id_name,fraction_ids Optionally restrict to the samples
#'   where `sample_data(physeq)[[fraction_id_name]] %in% fraction_ids` (e.g.
#'   selecting a single IgSeq sort fraction) before ordination. Both must be
#'   supplied together; default `NULL` (no restriction).
#' @param fit_filter_name,fit_filter_values Optionally fit the ordination
#'   model only on the samples where
#'   `sample_data(physeq)[[fit_filter_name]] %in% fit_filter_values`, then
#'   project the remaining samples onto that fit. Both must be supplied
#'   together; default `NULL` (fit on every sample).
#' @param transform_abundances Abundance transform applied via
#'   [microbiome::transform()] before computing distances/ordination (e.g.
#'   `"compositional"`, `"clr"`, `"hellinger"`). Default `"identity"` (no
#'   transform).
#' @param dist Distance metric. Used directly for `method = "PCoA"`/`"tSNE"`/
#'   `"UMAP"` (also accepts a precomputed `dist` object) and for
#'   `method = "dbRDA"`; ignored by the other methods. Default `"bray"`.
#' @param method Ordination method: `"PCoA"`, `"tSNE"`, `"UMAP"` (distance-
#'   based); `"NMDS"`, `"PCA"`, `"CA"`, `"DCA"` (direct on abundances); or
#'   `"CCA"`, `"RDA"`, `"dbRDA"` (constrained, require `model`).
#'   Case-insensitive. Default `"PCoA"`.
#' @param ndims Maximum number of ordination axes to retain in `coords`/
#'   `loadings`/`covariates`. Models can have up to `nsamples - 1` axes, but
#'   only a couple are ever plotted — keeping them all makes correlating
#'   every taxon against every axis for `loadings` (and storing the
#'   resulting matrices) far slower than necessary. Raise this if you need
#'   to plot axes beyond the default range. Default `10`.
#' @param model For constrained methods, a string with the right-hand side
#'   of the model formula (e.g. `"Var1+Var2"`). Required when `method` is
#'   `"CCA"`/`"RDA"`/`"dbRDA"`; ignored otherwise. Default `NULL`.
#' @param confounders Character vector of `sample_data` column names to
#'   partial out via `vegan::Condition()` in constrained models. Default
#'   none.
#' @param pca For `method = "tSNE"`, whether `Rtsne()` should run an initial
#'   PCA step. Default `TRUE`.
#' @param perplexity For `method = "tSNE"`, the perplexity parameter.
#'   Default `NULL` (30, capped to `floor((n - 1) / 3)` with a `warning()`
#'   if exceeded).
#' @param nb_neighbors For `method = "UMAP"`, the number of neighbors.
#'   Default `15`.
#' @param min_dist For `method = "UMAP"`, the minimum distance. Default
#'   `0.1`.
#'
#' @return A list:
#'   \describe{
#'     \item{fit}{The underlying model object returned by the `vegan`/
#'       `Rtsne`/`umap` fitting function.}
#'     \item{taxrank}{The `taxrank` argument, passed through.}
#'     \item{coords}{Sample (site) scores: a length-2 list, `[[1]]` for
#'       `vegan` scaling 1 (inter-sample distances) and `[[2]]` for scaling
#'       2 (inter-taxon angles). Identical for methods with only one natural
#'       scaling (PCoA, tSNE, UMAP).}
#'     \item{loadings}{Taxon (species) scores, same `[[1]]`/`[[2]]` scaling
#'       convention as `coords`. May be `NULL` per-scaling for methods that
#'       don't provide them (e.g. NMDS, DCA).}
#'     \item{covariates}{Constraint (biplot) scores, same scaling
#'       convention. `list()` for unconstrained methods.}
#'     \item{eigen_values}{Named numeric vector of eigenvalues (for scree
#'       plots / percent-explained-variance), or `NULL` for methods without
#'       them (NMDS).}
#'     \item{dist}{The `dist` argument, passed through.}
#'     \item{transform_abundances}{The `transform_abundances` argument,
#'       passed through.}
#'     \item{method}{`method`, in its original (non-lowercased) casing.}
#'     \item{model}{The (possibly `confounders`-augmented) model formula
#'       string, for constrained methods; `NULL` otherwise.}
#'     \item{fit_filter}{`list(name = fit_filter_name, values =
#'       fit_filter_values)`, or `NULL` if no fit filter was used.}
#'     \item{fit_sample_names}{Character vector of the fit-subset sample
#'       names, or `NULL` if all samples were used to fit.}
#'     \item{sample_data}{`sample_data(physeq)` as a `data.frame`. Gains a
#'       `.is_fit_sample` logical column when a fit filter was used.}
#'     \item{tax_table}{`tax_table(physeq)` as a `data.frame`.}
#'   }
#'
#' @seealso [scree_plot()], [stat_beta_diversity()], [plot_beta_diversity()]
#'
#' @examples
#' data(ps_16s_refinement)
#' bd <- get_beta_diversity(ps_16s_refinement, method = "PCoA", dist = "bray")
#' head(bd$coords[[1]])
#'
#' # Constrained ordination against a sample_data variable, agglomerated to
#' # genus level for speed
#' bd_cca <- get_beta_diversity(
#'   ps_16s_refinement,
#'   taxrank = "Genus",
#'   method = "CCA",
#'   model = "Mechanical_Lyse"
#' )
#' bd_cca$eigen_values
#'
#' @export
get_beta_diversity <- function(
  physeq,
  taxrank = NULL,
  fraction_id_name = NULL,
  fraction_ids = NULL,
  fit_filter_name = NULL,
  fit_filter_values = NULL,
  transform_abundances = "identity",
  dist = "bray",
  method = "PCoA",
  ndims = 10,
  model = NULL,
  confounders = c(),
  pca = TRUE,
  perplexity = NULL,
  nb_neighbors = 15,
  min_dist = 0.1
) {
  method_orig <- method # keep method name with the original case (used later for plot title)
  method <- tolower(method)

  .check_phyloseq(physeq)
  physeq <- reverseASV(physeq)

  physeq <- .filter_glom_transform_physeq(
    physeq,
    fraction_id_name = fraction_id_name,
    fraction_ids = fraction_ids,
    taxrank = taxrank,
    transform_abundances = transform_abundances
  )

  # Determine samples used to fit the ordination model (after all preprocessing)
  fit_sample_names <- NULL
  if (!is.null(fit_filter_name) && !is.null(fit_filter_values)) {
    fit_mask <- sample_data(physeq)[[fit_filter_name]] %in% fit_filter_values
    physeq_fit <- prune_samples(fit_mask, physeq)
    fit_sample_names <- sample_names(physeq_fit)
  } else {
    physeq_fit <- physeq
  }

  if (method %in% c("pcoa", "tsne", "umap")) {
    ord <- .get_beta_diversity_from_distance(
      physeq = physeq,
      physeq_fit = physeq_fit,
      dist = dist,
      method = method,
      ndims = ndims,
      fit_sample_names = fit_sample_names,
      pca = pca,
      perplexity = perplexity,
      nb_neighbors = nb_neighbors,
      min_dist = min_dist
    )
  } else if (method %in% c("nmds", "pca", "ca", "dca")) {
    ord <- .get_beta_diversity_from_abundance(
      physeq = physeq,
      physeq_fit = physeq_fit,
      dist = dist,
      method = method,
      ndims = ndims,
      fit_sample_names = fit_sample_names
    )
  } else if (!is.null(model) || length(confounders) > 0) {
    ord <- .get_beta_diversity_constrained(
      physeq = physeq,
      physeq_fit = physeq_fit,
      dist = dist,
      method = method,
      model = model,
      confounders = confounders,
      ndims = ndims,
      fit_sample_names = fit_sample_names
    )
  } else if (method %in% c("cca", "rda", "dbrda")) {
    stop("Model is required for constrained models")
  } else {
    stop("Invalid method")
  }

  eigen_values <- ord$eigen_values
  if (length(eigen_values) == 1 && is.na(eigen_values)) {
    eigen_values <- NULL
  }

  sample_data_df <- sample_data(physeq) %>% as("data.frame")
  if (!is.null(fit_sample_names)) {
    sample_data_df$.is_fit_sample <- rownames(sample_data_df) %in%
      fit_sample_names
  }

  return(list(
    fit = ord$fit,
    taxrank = taxrank,

    # lists: [[1]] - scaling 1, [[2]] - scaling 2 (see vegan's scaling)
    coords = ord$coords, # = sample (site) scores
    loadings = ord$loadings, # = variable (species) scores
    covariates = ord$covariates, # = covariate scores (if constrained model)

    eigen_values = eigen_values, # NULL if no eigen values provided for this method

    dist = dist, # distance metric
    transform_abundances = transform_abundances, # abundance transform applied before ordination
    method = method_orig, # name of ordination method
    model = if (method %in% c("rda", "cca", "dbrda")) {
      ord$model # confounders-augmented, from .get_beta_diversity_constrained()
    } else {
      NULL
    }, # formula with covariates
    fit_filter = if (!is.null(fit_filter_name)) {
      list(name = fit_filter_name, values = fit_filter_values)
    } else {
      NULL
    },
    fit_sample_names = fit_sample_names,
    sample_data = sample_data_df,
    tax_table = tax_table(physeq) %>% as.data.frame()
  ))
}

#' Run PERMANOVA or envfit on a single beta-diversity data slice
#'
#' Shared per-facet worker for [stat_beta_diversity()]: a factor/character
#' `Label` column runs [vegan::adonis2()] (PERMANOVA); a numeric `Label`
#' column runs [vegan::envfit()] (correlation-style test).
#'
#' @param stat_data A data frame with `Comp1`/`Comp2` (or one column per
#'   retained ordination axis), a `Label` column, and optionally a
#'   `strata_name` column. Must not contain the internal `.facet_row`/
#'   `.facet_col` marker columns.
#' @param strata_name Column of `stat_data` to use as `strata` in
#'   [vegan::adonis2()], or `NULL`.
#' @return A list with elements `test_result`, `p_raw`, `p_text`, `stat`
#'   (`"permanova"` or `"envfit"`).
#' @noRd
.stat_beta_diversity_test <- function(stat_data, strata_name) {
  if (is.factor(stat_data$Label) || is.character(stat_data$Label)) {
    stat_data <- na.omit(stat_data)
    if (!is.null(strata_name)) {
      strata <- stat_data[[strata_name]]
      stat_data[[strata_name]] <- NULL
    } else {
      strata <- NULL
    }
    test_result <- vegan::adonis2(
      stat_data[, colnames(stat_data) != "Label"] ~ Label,
      data = stat_data,
      strata = strata,
      permutations = 999,
      method = "euclidean"
    )
    p_raw <- test_result$`Pr(>F)`[1]
    p_text <- if (!is.null(strata)) {
      paste0(
        "restricted by ",
        strata_name,
        " PERMANOVA p = ",
        format(signif(p_raw, digits = 2), scientific = TRUE)
      )
    } else {
      paste0(
        "PERMANOVA p = ",
        format(signif(p_raw, digits = 2), scientific = TRUE)
      )
    }
    list(test_result = test_result, p_raw = p_raw, p_text = p_text, stat = "permanova")
  } else {
    test_result <- vegan::envfit(
      stat_data[, colnames(stat_data) != "Label"] ~ Label,
      data = stat_data,
      permutations = 999
    )
    p_raw <- test_result$vectors$pvals[1]
    p_text <- paste0(
      "Correlation p = ",
      format(signif(p_raw, digits = 2), scientific = TRUE)
    )
    list(test_result = test_result, p_raw = p_raw, p_text = p_text, stat = "envfit")
  }
}

#' Pairwise PERMANOVA between every pair of label levels
#'
#' Per-facet worker backing [stat_beta_diversity()]'s `pairwise = TRUE` path:
#' runs [vegan::adonis2()] on every pair of `Label` levels in `stat_data` and
#' BH-adjusts the resulting p-values across pairs.
#'
#' @inheritParams .stat_beta_diversity_test
#' @return A data frame with columns `group1`, `group2`, `R2`, `p_raw`,
#'   `p_adj`; `NULL` if `Label` is not categorical, has fewer than 2 levels,
#'   or no pair yields a fittable model.
#' @noRd
.stat_beta_diversity_pairwise <- function(stat_data, strata_name) {
  if (!is.factor(stat_data$Label) && !is.character(stat_data$Label)) {
    return(NULL)
  }
  labels <- unique(na.omit(as.character(stat_data$Label)))
  if (length(labels) < 2) {
    return(NULL)
  }
  pairs <- combn(labels, 2, simplify = FALSE)
  pair_results <- lapply(pairs, function(pair) {
    sub <- stat_data[as.character(stat_data$Label) %in% pair, ]
    sub <- na.omit(sub)
    if (length(unique(as.character(sub$Label))) < 2) {
      return(NULL)
    }
    sub_strata <- NULL
    if (!is.null(strata_name) && strata_name %in% colnames(sub)) {
      sub_strata <- sub[[strata_name]]
      sub[[strata_name]] <- NULL
    }
    res <- tryCatch(
      vegan::adonis2(
        sub[, colnames(sub) != "Label"] ~ Label,
        data = sub,
        strata = sub_strata,
        permutations = 999,
        method = "euclidean"
      ),
      error = function(e) NULL
    )
    if (is.null(res)) {
      return(NULL)
    }
    data.frame(
      group1 = pair[1],
      group2 = pair[2],
      R2 = res$R2[1],
      p_raw = res$`Pr(>F)`[1],
      stringsAsFactors = FALSE
    )
  })
  pair_results <- do.call(rbind, Filter(Negate(is.null), pair_results))
  if (is.null(pair_results) || nrow(pair_results) == 0) {
    return(NULL)
  }
  pair_results$p_adj <- p.adjust(pair_results$p_raw, method = "BH")
  pair_results
}

#' Permutation Test for Beta-Diversity Ordination
#'
#' Tests whether `label_name` explains variation in the ordination
#' coordinates returned by [get_beta_diversity()]: a factor/character
#' `label_name` runs a PERMANOVA ([vegan::adonis2()]); a numeric one runs a
#' correlation-style test ([vegan::envfit()]). Optionally repeats the test
#' once per facet level (wrap- or grid-mode faceting, mirroring
#' [plot_beta_diversity()]'s own faceting) and/or adds pairwise,
#' BH-adjusted comparisons between every pair of `label_name` levels.
#'
#' @param beta_dispersion_fit A list as returned by [get_beta_diversity()].
#' @param facet_mode `"wrap"` (facet by the single `facet` variable) or
#'   `"grid"` (facet by `facet_row` x `facet_col`; only takes effect when
#'   both are supplied). Defaults to `"wrap"` if `NULL`.
#' @param facet Wrap-mode faceting variable: a column name of
#'   `beta_dispersion_fit$sample_data`, or `NULL` for no faceting.
#' @param facet_row,facet_col Grid-mode row/column faceting variable names
#'   (column names of `beta_dispersion_fit$sample_data`).
#' @param comp Length-2 integer vector selecting which two ordination axes
#'   (columns of `beta_dispersion_fit$coords[[1]]`) to test, or `NULL` to
#'   test using every retained axis at once.
#' @param strata_name Optional column of `beta_dispersion_fit$sample_data`
#'   passed as `strata` to [vegan::adonis2()], restricting permutations to
#'   within each stratum (paired/blocked designs).
#' @param label_name Column of `beta_dispersion_fit$sample_data` to test
#'   against the ordination coordinates.
#' @param pairwise If `TRUE` and `label_name` is categorical, also compute
#'   pairwise PERMANOVA between every pair of `label_name` levels (per
#'   facet, if faceted), BH-adjusted across pairs.
#'
#' @return `NULL` (with a `warning()`) if `label_name` is not a column of
#'   `beta_dispersion_fit$sample_data`. Otherwise a list (facets whose
#'   `label_name` is constant are individually skipped with a `warning()`,
#'   not fatal -- `test_result`/`p_value`-family elements simply omit them,
#'   and can end up empty if every facet was skipped):
#'   \describe{
#'     \item{stat}{`"permanova"` or `"envfit"`.}
#'     \item{label_name, facet_mode, facet, facet_row, facet_col}{Echo the
#'       corresponding arguments (`facet`/`facet_row`/`facet_col` may differ
#'       from the input when only one grid dimension was supplied).}
#'     \item{facet_name}{The active single faceting variable.}
#'     \item{dim_used}{The `comp` argument, or `"all <n>"` when
#'       `comp = NULL`.}
#'     \item{test_result}{A list of raw [vegan::adonis2()]/[vegan::envfit()]
#'       objects, one per facet level tested (unnamed if unfaceted); `NULL`
#'       if every facet was skipped for having a constant `label_name`.}
#'     \item{p_value}{Wrap-mode/no-facet only: a character vector of
#'       human-readable p-value strings, one per facet level (unnamed if
#'       unfaceted).}
#'     \item{p_value_df}{Grid-mode only: a data frame with columns
#'       `facet_row`, `facet_col`, `p_label`, `p_raw`, one row per grid
#'       cell tested.}
#'     \item{p_value_raw}{Numeric p-values, named by facet level (wrap
#'       mode) or unnamed (grid mode uses row order matching
#'       `p_value_df`).}
#'     \item{pairwise_df}{If `pairwise = TRUE`: a data frame with columns
#'       `group1`, `group2`, `R2`, `p_raw`, `p_adj`, and (if faceted) the
#'       faceting variable(s); `NULL` if `pairwise = FALSE` or no pair was
#'       fittable.}
#'   }
#'
#' @examples
#' taxa_nms <- paste0("T", 1:12)
#' samp_nms <- paste0("S", 1:10)
#' ps <- phyloseq::phyloseq(
#'   phyloseq::otu_table(
#'     matrix(rpois(120, 20), nrow = 12, dimnames = list(taxa_nms, samp_nms)),
#'     taxa_are_rows = TRUE
#'   ),
#'   phyloseq::sample_data(data.frame(
#'     group = rep(c("A", "B"), 5),
#'     row.names = samp_nms
#'   )),
#'   phyloseq::tax_table(matrix(taxa_nms, dimnames = list(taxa_nms, "taxon_id")))
#' )
#' bd <- get_beta_diversity(ps, method = "PCoA", dist = "bray")
#' stat_beta_diversity(bd, label_name = "group")
#'
#' @export
stat_beta_diversity <- function(
  beta_dispersion_fit,
  facet_mode = "wrap", # "grid", "wrap"
  facet = NULL, # wrap-mode facet variable name
  facet_row = NULL, # grid-mode row facet variable name
  facet_col = NULL, # grid-mode col facet variable name
  comp = c(1, 2), # if NULL, consider all components
  strata_name = NULL, # in case of paired design
  label_name,
  pairwise = FALSE
) {
  if (is.null(facet_mode)) {
    facet_mode <- "wrap"
  }

  grid_mode <- facet_mode == "grid" && !is.null(facet_row) && !is.null(facet_col)
  active_facet <- if (facet_mode == "wrap") {
    facet
  } else if (!is.null(facet_row) && is.null(facet_col)) {
    facet_row
  } else if (is.null(facet_row) && !is.null(facet_col)) {
    facet_col
  } else {
    NULL
  }

  if (!label_name %in% colnames(beta_dispersion_fit$sample_data)) {
    warning("Wrong label name")
    return(NULL)
  }

  if (!is.null(comp)) {
    stat_data_full <- data.frame(
      Comp1 = beta_dispersion_fit$coords[[1]][, comp[1]],
      Comp2 = beta_dispersion_fit$coords[[1]][, comp[2]],
      Label = beta_dispersion_fit$sample_data[[label_name]]
    )
  } else {
    stat_data_full <- data.frame(
      beta_dispersion_fit$coords[[1]],
      Label = beta_dispersion_fit$sample_data[[label_name]]
    )
  }

  if (!is.null(strata_name)) {
    stat_data_full[[strata_name]] <- beta_dispersion_fit$sample_data[[strata_name]]
  }

  # Embed facet variables into stat_data_full before na.omit so subsetting is consistent
  if (
    grid_mode &&
      facet_row %in% colnames(beta_dispersion_fit$sample_data)
  ) {
    stat_data_full[[".facet_row"]] <- beta_dispersion_fit$sample_data[[facet_row]]
  }
  if (
    grid_mode &&
      facet_col %in% colnames(beta_dispersion_fit$sample_data)
  ) {
    stat_data_full[[".facet_col"]] <- beta_dispersion_fit$sample_data[[facet_col]]
  }
  if (
    !grid_mode &&
      !is.null(active_facet) &&
      active_facet %in% colnames(beta_dispersion_fit$sample_data)
  ) {
    stat_data_full[[".facet_row"]] <- beta_dispersion_fit$sample_data[[active_facet]]
  }

  stat_data_full <- na.omit(stat_data_full)

  dim_used <- if (!is.null(comp)) {
    comp
  } else {
    paste("all", ncol(beta_dispersion_fit$coords[[1]]))
  }

  # Build the facet combinations to iterate: grid mode crosses facet_row x
  # facet_col; wrap/no-facet mode is the same shape with a single, unused
  # "col" column (NA sentinel -- no legitimate facet value can be NA here,
  # since stat_data_full was just na.omit()'d).
  combinations <- if (grid_mode) {
    expand.grid(
      row = unique(stat_data_full[[".facet_row"]]),
      col = unique(stat_data_full[[".facet_col"]]),
      stringsAsFactors = FALSE
    )
  } else if (".facet_row" %in% colnames(stat_data_full)) {
    data.frame(
      row = unique(stat_data_full[[".facet_row"]]),
      col = NA,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(row = NA, col = NA, stringsAsFactors = FALSE)
  }

  cell_results <- list()
  for (i in seq_len(nrow(combinations))) {
    row_val <- combinations$row[i]
    col_val <- combinations$col[i]

    cell_data <- stat_data_full
    if (!is.na(row_val)) {
      cell_data <- cell_data[cell_data[[".facet_row"]] == row_val, ]
    }
    if (grid_mode) {
      cell_data <- cell_data[cell_data[[".facet_col"]] == col_val, ]
    }
    cell_data[[".facet_row"]] <- NULL
    cell_data[[".facet_col"]] <- NULL

    if (length(unique(cell_data$Label)) < 2) {
      warning(if (grid_mode) {
        paste0("Facet [", row_val, ", ", col_val, "]: label doesn't vary, skipping.")
      } else if (!is.na(row_val)) {
        paste0("Facet [", row_val, "]: label doesn't vary, skipping.")
      } else {
        "Variable doesn't vary, no stats performed..."
      })
      next
    }

    res <- .stat_beta_diversity_test(cell_data, strata_name)
    pw <- if (pairwise && res$stat == "permanova") {
      .stat_beta_diversity_pairwise(cell_data, strata_name)
    } else {
      NULL
    }
    if (!is.null(pw)) {
      if (grid_mode) {
        pw[[facet_row]] <- as.character(row_val)
        pw[[facet_col]] <- as.character(col_val)
      } else if (!is.na(row_val) && !is.null(active_facet)) {
        pw[[active_facet]] <- as.character(row_val)
      }
    }

    cell_results[[length(cell_results) + 1]] <- list(
      row = row_val,
      col = col_val,
      res = res,
      pairwise = pw
    )
  }

  stat_type <- if (length(cell_results) > 0) {
    cell_results[[length(cell_results)]]$res$stat
  } else {
    NULL
  }
  cell_keys <- if (grid_mode) {
    vapply(cell_results, function(c) paste(c$row, c$col, sep = "//"), character(1))
  } else if (length(cell_results) > 0 && !is.na(cell_results[[1]]$row)) {
    vapply(cell_results, function(c) as.character(c$row), character(1))
  } else {
    NULL
  }
  test_result_all <- if (length(cell_results) > 0) {
    stats::setNames(lapply(cell_results, function(c) c$res$test_result), cell_keys)
  } else {
    NULL
  }
  pairwise_list <- Filter(Negate(is.null), lapply(cell_results, `[[`, "pairwise"))
  pairwise_df <- if (pairwise && length(pairwise_list) > 0) {
    do.call(rbind, pairwise_list)
  } else {
    NULL
  }

  if (grid_mode) {
    p_value_df <- data.frame(
      facet_row = vapply(cell_results, function(c) as.character(c$row), character(1)),
      facet_col = vapply(cell_results, function(c) as.character(c$col), character(1)),
      p_label = vapply(cell_results, function(c) c$res$p_text, character(1)),
      p_raw = vapply(cell_results, function(c) c$res$p_raw, numeric(1)),
      stringsAsFactors = FALSE
    )
    return(list(
      stat = stat_type,
      label_name = label_name,
      facet_name = facet_row,
      facet_mode = facet_mode,
      facet = NULL,
      facet_row = facet_row,
      facet_col = facet_col,
      dim_used = dim_used,
      test_result = test_result_all,
      p_value = NULL,
      p_value_df = p_value_df,
      p_value_raw = p_value_df$p_raw,
      pairwise_df = pairwise_df
    ))
  }

  p_value_raw <- if (length(cell_results) > 0) {
    stats::setNames(vapply(cell_results, function(c) c$res$p_raw, numeric(1)), cell_keys)
  } else {
    numeric(0)
  }
  p_value_all <- if (length(cell_results) > 0) {
    stats::setNames(vapply(cell_results, function(c) c$res$p_text, character(1)), cell_keys)
  } else {
    character(0)
  }

  return(list(
    stat = stat_type,
    label_name = label_name,
    facet_name = active_facet,
    facet_mode = facet_mode,
    facet = facet,
    facet_row = facet_row,
    facet_col = facet_col,
    dim_used = dim_used,
    test_result = test_result_all,
    p_value = p_value_all,
    p_value_df = NULL,
    p_value_raw = p_value_raw,
    pairwise_df = pairwise_df
  ))
}

#' Compute and Plot Beta-Diversity
#'
#' Thin convenience wrapper chaining [get_beta_diversity()], optionally
#' [stat_beta_diversity()], and [plot_beta_diversity()] in one call.
#'
#' @inheritParams get_beta_diversity
#' @param group `NULL` (default) or the name of a `sample_data` column used
#'   to color/group the plot and, when `stat = TRUE`, as the
#'   [stat_beta_diversity()] test label.
#' @param stat Logical; run [stat_beta_diversity()] on `group` and annotate
#'   the plot with its result. Default `FALSE`.
#' @param comp Length-2 integer vector of ordination axes to plot and (when
#'   `stat = TRUE`) to test. Default `c(1, 2)`.
#' @param ... Passed on to [plot_beta_diversity()] (e.g. `facet`,
#'   `hover_variables`, `ellipses`; see its documentation for the full list).
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @seealso [get_beta_diversity()], [stat_beta_diversity()],
#'   [plot_beta_diversity()]
#'
#' @examples
#' data(ps_16s_refinement)
#' beta_diversity(ps_16s_refinement, group = "Protocol", stat = TRUE)
#'
#' @export
beta_diversity <- function(
  physeq,
  taxrank = NULL,
  fraction_id_name = NULL,
  fraction_ids = NULL,
  method = "PCoA",
  model = NULL,
  dist = "bray",
  group = NULL,
  stat = FALSE,
  comp = c(1, 2),
  ...
) {
  beta_dispersion_fit <- get_beta_diversity(
    physeq = physeq,
    taxrank = taxrank,
    fraction_id_name = fraction_id_name,
    fraction_ids = fraction_ids,
    dist = dist,
    method = method,
    model = model
  )

  stat_beta_dispersion <- if (isTRUE(stat)) {
    stat_beta_diversity(
      beta_dispersion_fit,
      comp = comp,
      label_name = group
    )
  } else {
    NULL
  }

  plot_beta_diversity(
    beta_dispersion_fit = beta_dispersion_fit,
    stat_beta_dispersion = stat_beta_dispersion,
    comp = comp,
    label_name = group,
    ...
  )
}

#' Scree Plot from Eigenvalues
#'
#' Bar chart of the percentage of variability explained by each ordination
#' axis, for judging how many components/axes to retain from an ordination
#' (e.g. the `eigen.values` element returned by [get_beta_diversity()]).
#'
#' @param eigen_values Numeric vector of eigenvalues, one per ordination
#'   axis, typically named (e.g. `"Axis.1"`, `"Axis.2"`, ...). Percentages
#'   are computed against the sum of *positive* eigenvalues only; negative
#'   eigenvalues (as can occur with PCoA on non-Euclidean distances) are
#'   still plotted but excluded from that total.
#' @param max_nb_comp Integer; maximum number of leading axes to plot.
#'   Default `10`.
#'
#' @return A `ggplot` object (bar chart of percent variability explained per
#'   axis), or `NULL` (with a `warning()`) if `eigen_values` is `NULL` or
#'   empty.
#'
#' @export
#'
#' @examples
#' eigen_values <- c(Axis.1 = 5, Axis.2 = 3, Axis.3 = 1, Axis.4 = 0.5)
#' scree_plot(eigen_values)
scree_plot <- function(eigen_values, max_nb_comp = 10) {
  if (is.null(eigen_values) || length(eigen_values) == 0) {
    warning("No eigen values are furnished")
    return(NULL)
  }
  total_var <- sum(eigen_values[eigen_values > 0])
  eigen_values <- eigen_values[1:min(max_nb_comp, length(eigen_values))]
  plot_data <- data.frame(
    prop_var = eigen_values / total_var * 100,
    dim = if (!is.null(names(eigen_values))) {
      names(eigen_values)
    } else {
      1:length(eigen_values)
    }
  )
  plot_data$dim <- factor(plot_data$dim, levels = plot_data$dim)

  plt <-
    ggplot(plot_data, aes(x = dim, y = prop_var)) +
    geom_bar(stat = "identity", fill = "skyblue", color = "black") +
    xlab("Dimension") +
    ylab("% of variability explained") +
    ggtitle("Scree plot") +
    .plot_title_theme()
  return(plt)
}

#' Pick one scaling's element out of a `coords`/`loadings`/`covariates` list
#'
#' Shared tail of the "resolve a `get_beta_diversity()` scaling list to one
#' matrix" pattern, used identically for `coords`, `loadings`, and
#' `covariates` in [plot_beta_diversity()]. Falls back to the first scaling
#' with a `warning()` if `scaling` is out of range; passes through unchanged
#' (without validating `scaling`) when the list has 0 or 1 elements, since
#' there's nothing to choose between.
#'
#' @param scaling_list A `coords`/`loadings`/`covariates`-style length-0/1/2
#'   list from [get_beta_diversity()].
#' @param scaling Requested scaling (`1` or `2`).
#' @return `list(value =, scaling =)`: the resolved matrix (or `NULL` if
#'   `scaling_list` was empty) and the (possibly corrected) `scaling`.
#' @noRd
.plot_beta_diversity_pick_scaling <- function(scaling_list, scaling) {
  if (length(scaling_list) == 0) {
    return(list(value = NULL, scaling = scaling))
  }
  if (length(scaling_list) == 1) {
    return(list(value = scaling_list[[1]], scaling = scaling))
  }
  if (0 < scaling && scaling <= length(scaling_list)) {
    return(list(value = scaling_list[[scaling]], scaling = scaling))
  }
  warning("Wrong scaling, scaling 1 is used")
  list(value = scaling_list[[1]], scaling = 1)
}

#' Resolve `coords`/`loadings`/`covariates` and the plotted axes
#'
#' For [plot_beta_diversity()]: selects the requested `scaling` out of
#' `beta_dispersion_fit`'s `coords`/`loadings`/`covariates` lists (falling
#' back to scaling 1 with a `warning()` if invalid — the scaling correction
#' from `coords` cascades into the `loadings`/`covariates` lookups, matching
#' `get_beta_diversity()`'s own scaling-list ordering), validates `comp`
#' against the resolved `coords` (falling back to `c(1, 2)` with a
#' `warning()`), drops `loadings`/`covariates` that don't cover `comp` (with
#' a `warning()`), and applies `reverse_dim1`/`reverse_dim2` axis negation.
#' None of this depends on sample-level filtering (it only touches matrix
#' columns/axes, never rows/samples), so it can run before or after
#' .plot_beta_diversity_filter_levels() with an identical result.
#'
#' @param beta_dispersion_fit A list as returned by [get_beta_diversity()].
#' @param scaling,comp,reverse_dim1,reverse_dim2 See [plot_beta_diversity()].
#' @return `list(coords =, loadings =, covariates =, scaling =, comp =,
#'   dim_names =)`.
#' @noRd
.plot_beta_diversity_select_scaling <- function(
  beta_dispersion_fit,
  scaling,
  comp,
  reverse_dim1,
  reverse_dim2
) {
  coords_sel <- .plot_beta_diversity_pick_scaling(
    beta_dispersion_fit$coords,
    scaling
  )
  coords <- coords_sel$value
  scaling <- coords_sel$scaling

  loadings <- .plot_beta_diversity_pick_scaling(
    beta_dispersion_fit$loadings,
    scaling
  )$value
  covariates <- .plot_beta_diversity_pick_scaling(
    beta_dispersion_fit$covariates,
    scaling
  )$value

  dim_names <- if (!is.null(colnames(coords))) {
    colnames(coords)
  } else {
    paste0("Dim ", seq_len(ncol(coords)))
  }

  if (
    !is.numeric(comp) ||
      length(comp) != 2 ||
      !all(comp %in% seq_len(ncol(coords)))
  ) {
    warning("Wrong components, forced to first two components")
    comp <- c(1, 2)
  }

  if (!is.null(loadings) && !all(comp %in% seq_len(ncol(loadings)))) {
    warning("No species score provided for these components")
    loadings <- NULL
  }
  if (!is.null(covariates) && !all(comp %in% seq_len(ncol(covariates)))) {
    warning("No species score provided for these components")
    covariates <- NULL
  }

  if (reverse_dim1) {
    if (!is.null(coords)) {
      coords[, comp[1]] <- -coords[, comp[1]]
    }
    if (!is.null(loadings)) {
      loadings[, comp[1]] <- -loadings[, comp[1]]
    }
    if (!is.null(covariates)) {
      covariates[, comp[1]] <- -covariates[, comp[1]]
    }
  }
  if (reverse_dim2) {
    if (!is.null(coords)) {
      coords[, comp[2]] <- -coords[, comp[2]]
    }
    if (!is.null(loadings)) {
      loadings[, comp[2]] <- -loadings[, comp[2]]
    }
    if (!is.null(covariates)) {
      covariates[, comp[2]] <- -covariates[, comp[2]]
    }
  }

  list(
    coords = coords,
    loadings = loadings,
    covariates = covariates,
    scaling = scaling,
    comp = comp,
    dim_names = dim_names
  )
}

#' Filter `sample_data`/`coords` to a variable's kept levels
#'
#' For [plot_beta_diversity()]: shared tail of the "optionally subset
#' samples to `keep_levels()`, then `factorize_levels()`" pattern, applied
#' identically to `label_name`, `shape_name`, `facet`, `facet_row`,
#' `facet_col`, and `animation_variable_name`. A no-op when `var_name` is
#' `NULL` or already numeric. When `var_levels` is `NULL`: a no-op for
#' `label_name`/`shape_name`/`animation_variable_name` (`auto_factor =
#' FALSE`), but still coerces to a plain (naturally-ordered) factor for
#' `facet`/`facet_row`/`facet_col` (`auto_factor = TRUE`), matching each
#' variable's pre-refactor behavior.
#'
#' @param sample_data,coords Current (possibly already filtered by an
#'   earlier call) `sample_data`/`coords`.
#' @param var_name `NULL` or a `sample_data` column name.
#' @param var_levels `NULL` or the subset/order of levels to keep, per
#'   `keep_levels()`/`factorize_levels()`.
#' @param auto_factor Logical; coerce `sample_data[[var_name]]` to a factor
#'   even when `var_levels` is `NULL`. Default `FALSE`.
#' @return `list(sample_data =, coords =)`.
#' @noRd
.plot_beta_diversity_filter_levels <- function(
  sample_data,
  coords,
  var_name,
  var_levels,
  auto_factor = FALSE
) {
  if (is.null(var_name) || is.numeric(sample_data[[var_name]])) {
    return(list(sample_data = sample_data, coords = coords))
  }

  if (!is.null(var_levels)) {
    keep <- keep_levels(sample_data[[var_name]], var_levels)
    sample_data <- sample_data[keep, , drop = FALSE]
    if (!is.null(coords)) {
      coords <- coords[keep, , drop = FALSE]
    }
    sample_data[[var_name]] <- factorize_levels(
      sample_data[[var_name]],
      var_levels
    )
  } else if (auto_factor) {
    sample_data[[var_name]] <- as.factor(sample_data[[var_name]])
  }

  list(sample_data = sample_data, coords = coords)
}

#' Extract p-value(s) for the subtitle/facet strips/facet annotations
#'
#' For [plot_beta_diversity()]: reads a [stat_beta_diversity()] result and
#' picks the right shape of p-value to embed: a full grid's `p_value_df`
#' (rendered later as `geom_text()` annotations by
#' .plot_beta_diversity_facets()), a single subtitle-embedded p-value (no
#' facet, or a facet whose levels don't need individual p-values), or a
#' named-by-level vector (one p-value per facet-strip label).
#'
#' @param stat_beta_dispersion `NULL` or a [stat_beta_diversity()] result.
#' @param sample_data Current (post-level-filtering) `sample_data`.
#' @param grid_mode,facet,facet_row,facet_col See [plot_beta_diversity()].
#' @return `list(p_value =, p_value_df =)`.
#' @noRd
.plot_beta_diversity_pvalues <- function(
  stat_beta_dispersion,
  sample_data,
  grid_mode,
  facet,
  facet_row,
  facet_col
) {
  p_value <- NULL
  p_value_df <- NULL

  if (!is.null(stat_beta_dispersion$label_name)) {
    if (
      grid_mode &&
        !is.null(facet_row) &&
        !is.null(facet_col) &&
        !is.null(stat_beta_dispersion$p_value_df)
    ) {
      # Full grid mode: p-values as data frame, rendered as geom_text annotations
      p_value_df <- stat_beta_dispersion$p_value_df
    } else if (!is.null(stat_beta_dispersion$p_value)) {
      active_facet <- if (grid_mode) facet_row %||% facet_col else facet
      if (
        is.null(active_facet) && is.null(names(stat_beta_dispersion$p_value))
      ) {
        # No facets: embed single p-value in subtitle
        p_value <- paste(
          strsplit(stat_beta_dispersion$p_value[[1]], "\n")[[1]],
          collapse = " "
        )
      } else if (
        !is.null(active_facet) &&
          !is.null(names(stat_beta_dispersion$p_value)) &&
          all(
            levels(sample_data[[active_facet]]) %in%
              names(stat_beta_dispersion$p_value)
          )
      ) {
        # Single-facet: embed per-facet p-value into facet strip labels
        p_value <- stat_beta_dispersion$p_value[levels(sample_data[[
          active_facet
        ]])]
      }
    }
  }

  list(p_value = p_value, p_value_df = p_value_df)
}

#' Assemble `plot_beta_diversity()`'s point-layer data frame
#'
#' Builds the `Comp1`/`Comp2` base data frame and adds every optional
#' column: `label`, the animation variable (under its own dynamic column
#' name), grid/wrap facet column(s) (with per-facet p-values embedded in
#' strip labels when `p_value` is named), `shape` (only for a
#' character/factor `shape_name`), and hover text (via `get_hover_text()`,
#' `utils.R`). `size` is resolved but not added as a `plot_df` column, since
#' it's mapped via the *local* `size` value directly in
#' .plot_beta_diversity_base_plot()'s `aes()` — same treatment as `label`/
#' `shape`.
#'
#' @param sample_data,coords,comp Post-filtering `sample_data`/`coords`, and
#'   the resolved `comp` (from .plot_beta_diversity_select_scaling()).
#' @param label `NULL` or the resolved (factorized, if character) label
#'   vector.
#' @param animation_variable_name,animation_variable `NULL`/`NULL`, or the
#'   animation column name and its resolved value vector.
#' @param grid_mode,facet,facet_row,facet_col See [plot_beta_diversity()].
#' @param p_value From .plot_beta_diversity_pvalues(); embedded into facet
#'   strip labels when it's named by facet level.
#' @param shape_name,size_name,hover_variables See [plot_beta_diversity()].
#' @return `list(plot_df =, shape =, size =)`.
#' @noRd
.plot_beta_diversity_build_df <- function(
  sample_data,
  coords,
  comp,
  label,
  animation_variable_name,
  animation_variable,
  grid_mode,
  facet,
  facet_row,
  facet_col,
  p_value,
  shape_name,
  size_name,
  hover_variables
) {
  plot_df <- data.frame(Comp1 = coords[, comp[1]], Comp2 = coords[, comp[2]])

  if (!is.null(label)) {
    plot_df$label <- label
  }

  if (!is.null(animation_variable)) {
    plot_df[[animation_variable_name]] <- animation_variable
  }

  if (grid_mode) {
    if (!is.null(facet_row)) {
      facet_row_data <- sample_data[[facet_row]]
      if (!is.factor(facet_row_data)) {
        facet_row_data <- factor(facet_row_data)
      }
      if (is.null(facet_col) && length(p_value) > 1) {
        plot_df$facet_row <- factor(
          paste0(
            facet_row,
            " = ",
            as.character(facet_row_data),
            "\n",
            p_value[as.character(facet_row_data)]
          ),
          levels = paste0(
            facet_row,
            " = ",
            as.character(levels(facet_row_data)),
            "\n",
            p_value[as.character(levels(facet_row_data))]
          )
        )
      } else {
        plot_df$facet_row <- factor(
          paste0(facet_row, " = ", as.character(facet_row_data)),
          levels = paste0(
            facet_row,
            " = ",
            as.character(levels(facet_row_data))
          )
        )
      }
    }
    if (!is.null(facet_col)) {
      facet_col_data <- sample_data[[facet_col]]
      if (!is.factor(facet_col_data)) {
        facet_col_data <- factor(facet_col_data)
      }
      if (is.null(facet_row) && length(p_value) > 1) {
        plot_df$facet_col <- factor(
          paste0(
            facet_col,
            " = ",
            as.character(facet_col_data),
            "\n",
            p_value[as.character(facet_col_data)]
          ),
          levels = paste0(
            facet_col,
            " = ",
            as.character(levels(facet_col_data)),
            "\n",
            p_value[as.character(levels(facet_col_data))]
          )
        )
      } else {
        plot_df$facet_col <- factor(
          paste0(facet_col, " = ", as.character(facet_col_data)),
          levels = paste0(
            facet_col,
            " = ",
            as.character(levels(facet_col_data))
          )
        )
      }
    }
  } else if (!is.null(facet)) {
    facet_data <- sample_data[[facet]]
    if (is.character(facet_data) || is.factor(facet_data)) {
      if (!is.factor(facet_data)) {
        facet_data <- as.factor(facet_data)
      }
      plot_df$facet <- factor(
        paste0(
          facet,
          " = ",
          as.character(facet_data),
          if (length(p_value) > 1) {
            paste0("\n", p_value[as.character(facet_data)])
          }
        ),
        levels = paste0(
          facet,
          " = ",
          as.character(levels(facet_data)),
          if (length(p_value) > 1) {
            paste0("\n", p_value[as.character(levels(facet_data))])
          }
        )
      )
    }
  }

  shape <- NULL
  if (!is.null(shape_name)) {
    shape_vals <- sample_data[[shape_name]]
    if (is.character(shape_vals) || is.factor(shape_vals)) {
      shape <- if (is.factor(shape_vals)) shape_vals else as.factor(shape_vals)
      plot_df$shape <- shape
    }
  }

  size <- if (!is.null(size_name)) sample_data[[size_name]] else NULL

  plot_df$hover_text <- get_hover_text(sample_data, hover_variables)

  if (".is_fit_sample" %in% colnames(sample_data)) {
    plot_df$.is_fit_sample <- sample_data$.is_fit_sample
  }

  list(plot_df = plot_df, shape = shape, size = size)
}

#' Build `plot_beta_diversity()`'s base point layer(s) and ellipses
#'
#' Draws sample points, split into two `geom_point()` layers (at
#' `projected_alpha`/`point_alpha`) when `plot_df` has a `.is_fit_sample`
#' column (i.e. `beta_dispersion_fit` used a `fit_filter`), otherwise one
#' layer at `point_alpha`; `size_name` additionally picks a
#' `ggplot2::scale_size()` range vs. a fixed `point_size`. Adds a
#' [ggplot2::stat_ellipse()] per level of `label` when `ellipses = TRUE` and
#' `label` is discrete.
#'
#' `label`/`shape`/`size` are mapped via `aes()`'s lazy (quosure) evaluation
#' against these *local* values, not `plot_df` columns — required so the
#' mapping still resolves correctly when a value is `NULL` (no `plot_df`
#' column to fall back on).
#'
#' @param plot_df,shape,size From .plot_beta_diversity_build_df().
#' @param label `NULL` or the resolved label vector (see
#'   .plot_beta_diversity_build_df()'s `label` argument).
#' @param label_name,shape_name,size_name,point_alpha,point_size,
#'   projected_alpha,ellipses,fill_ellipses See [plot_beta_diversity()].
#' @return A `ggplot` object.
#' @noRd
.plot_beta_diversity_base_plot <- function(
  plot_df,
  label,
  shape,
  size,
  label_name,
  shape_name,
  size_name,
  point_alpha,
  point_size,
  projected_alpha,
  ellipses,
  fill_ellipses
) {
  point_aes <- aes(
    x = Comp1,
    y = Comp2,
    color = label, # OK if NULL
    text = hover_text,
    shape = shape, # OK if NULL
    size = size # OK if NULL
  )

  if (".is_fit_sample" %in% colnames(plot_df)) {
    if (!is.null(size_name)) {
      plt <- ggplot() +
        geom_point(
          point_aes,
          plot_df[!plot_df$.is_fit_sample, ],
          alpha = projected_alpha
        ) +
        geom_point(
          point_aes,
          plot_df[plot_df$.is_fit_sample, ],
          alpha = point_alpha
        ) +
        ggplot2::scale_size(range = c(point_size * 0.5, point_size * 3)) +
        labs(color = label_name, shape = shape_name, size = size_name)
    } else {
      plt <- ggplot() +
        geom_point(
          point_aes,
          plot_df[!plot_df$.is_fit_sample, ],
          alpha = projected_alpha,
          size = point_size
        ) +
        geom_point(
          point_aes,
          plot_df[plot_df$.is_fit_sample, ],
          alpha = point_alpha,
          size = point_size
        ) +
        labs(color = label_name, shape = shape_name, size = size_name)
    }
  } else {
    if (!is.null(size_name)) {
      plt <- ggplot() +
        geom_point(point_aes, plot_df, alpha = point_alpha) +
        ggplot2::scale_size(range = c(point_size * 0.5, point_size * 3)) +
        labs(color = label_name, shape = shape_name, size = size_name)
    } else {
      plt <- ggplot() +
        geom_point(point_aes, plot_df, alpha = point_alpha, size = point_size) +
        labs(color = label_name, shape = shape_name, size = size_name)
    }
  }

  if (ellipses && is.factor(label)) {
    if (fill_ellipses) {
      plt <- plt +
        stat_ellipse(
          aes(x = Comp1, y = Comp2, color = label, fill = label),
          plot_df,
          geom = "polygon",
          alpha = 0.2
        )
    } else {
      plt <- plt +
        stat_ellipse(aes(x = Comp1, y = Comp2, color = label), plot_df)
    }
  }

  plt
}

#' Add axis labels, title, subtitle, and theme to a `plot_beta_diversity()` plot
#'
#' Builds the x/y axis labels (with percent-variance-explained, when
#' available), the title (method, and model formula for constrained
#' methods), and the subtitle (taxrank, fit-filter subset, NA-removal,
#' distance, scaling, biplot arrow cutoffs, and a single embedded p-value
#' when there's no facet or only one).
#'
#' @param plt A `ggplot` object (from .plot_beta_diversity_base_plot()).
#' @param beta_dispersion_fit A list as returned by [get_beta_diversity()].
#' @param dim_names,comp,scaling From .plot_beta_diversity_select_scaling().
#' @param prop_var_explained `NULL` or a named numeric vector (percent of
#'   variance explained per axis, from `beta_dispersion_fit$eigen_values`).
#' @param remove_na_from_plot,biplot,biplot_loadings,biplot_covariates,
#'   arrow_cutoff_load,arrow_cutoff_covar See [plot_beta_diversity()].
#' @param loadings,covariates From .plot_beta_diversity_select_scaling();
#'   only their presence/absence is used here (arrow cutoff subtitle text).
#' @param p_value From .plot_beta_diversity_pvalues().
#' @return `plt`, with labs/theme layers added.
#' @noRd
.plot_beta_diversity_labs <- function(
  plt,
  beta_dispersion_fit,
  dim_names,
  comp,
  prop_var_explained,
  remove_na_from_plot,
  scaling,
  biplot,
  biplot_loadings,
  loadings,
  biplot_covariates,
  covariates,
  arrow_cutoff_load,
  arrow_cutoff_covar,
  p_value
) {
  plt +
    labs(
      x = paste0(
        dim_names[comp[1]],
        if (!is.null(prop_var_explained)) {
          paste0(" (", round(prop_var_explained[comp[1]]), "%)")
        } else {
          NULL
        }
      ),
      y = paste0(
        dim_names[comp[2]],
        if (!is.null(prop_var_explained)) {
          paste0(" (", round(prop_var_explained[comp[2]]), "%)")
        } else {
          NULL
        }
      ),
      # get rid of backticks in the model (formula) for the title
      title = paste0(
        "Beta-Diversity",
        if (is.null(beta_dispersion_fit$model)) {
          ""
        } else {
          paste0(" ~ ", gsub("`", " ", beta_dispersion_fit$model))
        },
        " (",
        beta_dispersion_fit$method,
        ")"
      ),

      subtitle = paste0(
        "taxa agglom: ",
        if (!is.null(beta_dispersion_fit$taxrank)) {
          beta_dispersion_fit$taxrank
        } else {
          "none"
        },
        paste0(
          "  transform: ",
          beta_dispersion_fit$transform_abundances %||% "identity"
        ),
        if (!is.null(beta_dispersion_fit$fit_filter)) {
          paste0(
            "  fit subset: ",
            beta_dispersion_fit$fit_filter$name,
            " in {",
            paste(beta_dispersion_fit$fit_filter$values, collapse = ", "),
            "}"
          )
        },
        if (remove_na_from_plot) {
          "  NA's removed  "
        },
        if (!is.null(beta_dispersion_fit$dist)) {
          paste0(" dist=", beta_dispersion_fit$dist)
        },
        if (
          is.null(beta_dispersion_fit$method) ||
            !tolower(beta_dispersion_fit$method) %in% c("pcoa", "tsne", "umap")
        ) {
          paste0("  scaling: ", scaling)
        } else {
          NULL
        },
        if (biplot && (!is.null(loadings) || !is.null(covariates))) {
          "  arrow cut-offs:"
        } else {
          NULL
        },
        if (biplot_loadings && !is.null(loadings)) {
          paste0("  taxa=", arrow_cutoff_load)
        } else {
          NULL
        },
        if (biplot_covariates && !is.null(covariates)) {
          paste0("  predictors=", arrow_cutoff_covar)
        } else {
          NULL
        },
        # add stats to the title only if there's no facet or only one facet
        if (length(p_value) == 1) {
          paste0("\n", p_value)
        } else {
          NULL
        }
      )
    ) +
    .plot_title_theme() +
    ggsci::scale_fill_npg()
}

#' Apply grid/wrap faceting (and grid p-value annotations) to a
#' `plot_beta_diversity()` plot
#'
#' Grid mode facets by `facet_row`/`facet_col` (whichever are non-`NULL`),
#' additionally annotating each panel with `p_value_df`'s per-cell p-value
#' text when a full grid was tested; wrap mode facets by `facet` (only if
#' .plot_beta_diversity_build_df() actually added a `facet` column, i.e.
#' the faceting variable was character/factor).
#'
#' @param plt A `ggplot` object.
#' @param plot_df From .plot_beta_diversity_build_df(); supplies the
#'   `facet`/`facet_row`/`facet_col` columns faceted on.
#' @param grid_mode,facet,facet_row,facet_col See [plot_beta_diversity()].
#' @param p_value_df From .plot_beta_diversity_pvalues().
#' @return `plt`, with a facet layer (and grid annotations) added.
#' @noRd
.plot_beta_diversity_facets <- function(
  plt,
  plot_df,
  grid_mode,
  facet,
  facet_row,
  facet_col,
  p_value_df
) {
  if (grid_mode) {
    if (!is.null(facet_row) && !is.null(facet_col)) {
      plt <- plt + facet_grid(facet_row ~ facet_col, scales = "fixed")
      if (!is.null(p_value_df) && nrow(p_value_df) > 0) {
        annot_df <- data.frame(
          facet_row = factor(
            paste0(facet_row, " = ", p_value_df$facet_row),
            levels = levels(plot_df$facet_row)
          ),
          facet_col = factor(
            paste0(facet_col, " = ", p_value_df$facet_col),
            levels = levels(plot_df$facet_col)
          ),
          p_label = p_value_df$p_label,
          x = Inf,
          y = Inf
        )
        plt <- plt +
          geom_text(
            data = annot_df,
            aes(x = x, y = y, label = p_label),
            hjust = 1.05,
            vjust = 1.5,
            size = 3,
            inherit.aes = FALSE
          )
      }
    } else if (!is.null(facet_row)) {
      plt <- plt + facet_grid(rows = vars(facet_row), scales = "fixed")
    } else if (!is.null(facet_col)) {
      plt <- plt + facet_grid(cols = vars(facet_col), scales = "fixed")
    }
  } else if (!is.null(facet) && "facet" %in% colnames(plot_df)) {
    plt <- plt +
      facet_wrap(
        . ~ facet,
        scales = "fixed",
        ncol = smart_facet_ncol(nlevels(plot_df$facet))
      )
  }
  plt
}

#' Filter and scale one biplot arrow group (loadings or covariates)
#'
#' For [plot_beta_diversity()]'s biplot: shared tail of the "drop NA/absent
#' rows, filter out arrows shorter than a length cutoff, build a
#' `Comp1`/`Comp2`/`type` data frame, and scale it onto `coords`' range"
#' pattern applied identically to loadings and covariates by
#' .plot_beta_diversity_arrow_data(). Hover text/label columns are added by
#' the caller instead, since those differ between the two arrow types.
#'
#' @param scores A `loadings` or `covariates` matrix (rows = taxa/
#'   predictors, columns = ordination axes).
#' @param type `"load"` or `"covar"`, recorded in the returned data frame's
#'   `type` column.
#' @param comp,coords See [plot_beta_diversity()]; `coords` sets the range
#'   arrow length is scaled to.
#' @param cutoff Arrows shorter than `cutoff` times the longest arrow in this
#'   group are dropped (`arrow_cutoff_load`/`arrow_cutoff_covar`, see
#'   [plot_beta_diversity()]).
#' @return `list(arrow =, scale_factor =)`; both `NULL` if nothing survives
#'   filtering.
#' @noRd
.plot_beta_diversity_arrow_group <- function(scores, type, comp, coords, cutoff) {
  # Account for absent scores: taxa/predictors with zero variance across the
  # fitted samples have an undefined correlation/score (e.g. PCoA's loadings
  # are cor(otu_matrix, coords), see .get_beta_diversity_from_distance()) and
  # are dropped here rather than plotted as NA arrows.
  filtered <- na.omit(scores[, comp, drop = FALSE])

  # Filter arrows by length (get rid of smaller arrows)
  arrow_lengths <- apply(filtered, 1, function(x) sqrt(sum(x^2)))
  filtered <- filtered[arrow_lengths / max(arrow_lengths) > cutoff, , drop = FALSE]

  if (nrow(filtered) == 0) {
    return(list(arrow = NULL, scale_factor = NULL))
  }

  arrow <- data.frame(Comp1 = filtered[, 1], Comp2 = filtered[, 2], type = type)
  rownames(arrow) <- rownames(filtered)

  # Scale the secondary axes so that arrows and sample scores are on a
  # comparable scale. NOTE: both components are scaled by the same factor,
  # so arrows point in the same direction as before, only their size changes.
  # na.rm = TRUE: coords can have NA rows for samples excluded from
  # projection (e.g. .predictor_complete_cases()/.warn_and_na_fill()) --
  # without it, a single NA row poisons max() and NAs out every arrow.
  lim <- max(abs(arrow$Comp1), abs(arrow$Comp2))
  scale_factor <- max(abs(coords), na.rm = TRUE) / lim
  arrow$Comp1 <- arrow$Comp1 * scale_factor
  arrow$Comp2 <- arrow$Comp2 * scale_factor

  list(arrow = arrow, scale_factor = scale_factor)
}

#' Build the biplot arrow data frame and length-scale factors
#'
#' For [plot_beta_diversity()]: turns `loadings`/`covariates` into the arrow
#' data frame drawn by .plot_beta_diversity_arrow_plot(), gated by
#' `biplot_loadings`/`biplot_covariates`. Loading arrows and covariate
#' arrows share the filter/scale steps
#' (.plot_beta_diversity_arrow_group()), but differ in hover text/label
#' source: loadings use `tax_table`/`arrow_taxonomy_labels` (species scores
#' are taxa), covariates always use their own (backtick-stripped) name,
#' since they're `sample_data` column names, not taxa.
#'
#' @param loadings,covariates `NULL` or a scores matrix, as resolved by
#'   .plot_beta_diversity_select_scaling().
#' @param comp,coords,biplot_loadings,biplot_covariates,arrow_cutoff_load,
#'   arrow_cutoff_covar,arrow_taxonomy_labels See [plot_beta_diversity()].
#' @param tax_table `beta_dispersion_fit$tax_table`, or `NULL`.
#' @return `list(arrow_df =, scale_factor_load =, scale_factor_covar =)`;
#'   `arrow_df` has 0 rows when neither arrow group survives filtering.
#' @noRd
.plot_beta_diversity_arrow_data <- function(
  loadings,
  covariates,
  comp,
  coords,
  biplot_loadings,
  biplot_covariates,
  arrow_cutoff_load,
  arrow_cutoff_covar,
  tax_table,
  arrow_taxonomy_labels
) {
  arrow_df <- data.frame()
  scale_factor_load <- NULL
  scale_factor_covar <- NULL

  # Plot species arrows = loadings here
  if (!is.null(loadings) && biplot_loadings) {
    grp <- .plot_beta_diversity_arrow_group(
      loadings,
      "load",
      comp,
      coords,
      arrow_cutoff_load
    )
    if (!is.null(grp$arrow)) {
      arrow_load <- grp$arrow
      scale_factor_load <- grp$scale_factor

      # If phyloseq object provided, get hover info from it - taxonomy for each species (arrow)
      if (!is.null(tax_table)) {
        arrow_load$hover_text <- get_hover_text(
          tax_table[rownames(arrow_load), , drop = FALSE],
          colnames(tax_table)
        )
      } else {
        arrow_load$hover_text <- rownames(arrow_load)
      }

      if (!is.null(arrow_taxonomy_labels)) {
        # use taxonomy as label instead of taxa names
        arrow_load$Names <- apply(
          tax_table[rownames(arrow_load), ],
          1,
          function(row) {
            gsub("NA", " ", paste(row[arrow_taxonomy_labels], collapse = ":"))
          }
        )
      } else {
        arrow_load$Names <- rownames(arrow_load)
      }
      arrow_df <- rbind(arrow_df, arrow_load)
    }
  }

  # Similar for covariates (predictors in case of constrained model), scaled
  # independently so both groups are of comparable size on the plot
  if (!is.null(covariates) && biplot_covariates) {
    grp <- .plot_beta_diversity_arrow_group(
      covariates,
      "covar",
      comp,
      coords,
      arrow_cutoff_covar
    )
    if (!is.null(grp$arrow)) {
      arrow_covar <- grp$arrow
      scale_factor_covar <- grp$scale_factor
      # Hover text is just the name of covariate
      arrow_covar$hover_text <- rownames(arrow_covar)
      arrow_covar$Names <- gsub("`", "", rownames(arrow_covar))
      arrow_df <- rbind(arrow_df, arrow_covar)
    }
  }

  list(
    arrow_df = arrow_df,
    scale_factor_load = scale_factor_load,
    scale_factor_covar = scale_factor_covar
  )
}

#' Draw biplot arrow layers (loadings/covariates) and their secondary axis
#'
#' For [plot_beta_diversity()]: adds `geom_segment()` arrow layers for
#' whichever of `"load"`/`"covar"` rows are present in `arrow_df`, optional
#' arrow labels (`geom_text()` or `ggrepel::geom_text_repel()`), and a
#' secondary axis scaled back to the arrows' pre-scaling units. The
#' secondary axis is only added when exactly one of the two arrow groups is
#' present and no `marginal_plot` is set, to avoid an ambiguous/cluttered
#' axis (loadings and covariates use independent scale factors, and
#' `marginal_plot` has no room for it).
#'
#' @param plt The base `ggplot` object (already carrying the sample points).
#' @param arrow_df As returned by .plot_beta_diversity_arrow_data(); must
#'   have at least 1 row.
#' @param scale_factor_load,scale_factor_covar As returned by
#'   .plot_beta_diversity_arrow_data().
#' @param color_arrows_by_taxa,arrow_labels,repel,max_overlaps,marginal_plot
#'   See [plot_beta_diversity()].
#' @param dim_names,comp As resolved by
#'   .plot_beta_diversity_select_scaling().
#' @return `plt` with the arrow layers/secondary axis added.
#' @noRd
.plot_beta_diversity_arrow_plot <- function(
  plt,
  arrow_df,
  scale_factor_load,
  scale_factor_covar,
  color_arrows_by_taxa,
  arrow_labels,
  repel,
  max_overlaps,
  marginal_plot,
  dim_names,
  comp
) {
  arrow_params <- arrow(
    type = "closed",
    angle = 20,
    length = unit(0.1, "inches")
  )

  if ("load" %in% arrow_df$type) {
    if (color_arrows_by_taxa) {
      plt <- plt +
        geom_segment(
          data = arrow_df[arrow_df$type == "load", ],
          aes(
            x = 0,
            y = 0,
            xend = Comp1,
            yend = Comp2,
            text = hover_text,
            color = Names
          ), # Arrow from (0,0) to each point
          arrow = arrow_params,
          size = 0.7, # Size of the arrows
          alpha = 0.7 # Transparency of the arrows
        )
    } else {
      plt <- plt +
        geom_segment(
          data = arrow_df[arrow_df$type == "load", ],
          aes(x = 0, y = 0, xend = Comp1, yend = Comp2, text = hover_text), # Arrow from (0,0) to each point
          arrow = arrow_params,
          color = "darkgrey",
          size = 0.7, # Size of the arrows
          alpha = 0.7 # Transparency of the arrows
        )
    }
  }

  if ("covar" %in% arrow_df$type) {
    plt <- plt +
      geom_segment(
        data = arrow_df[arrow_df$type == "covar", ],
        aes(x = 0, y = 0, xend = Comp1, yend = Comp2, text = hover_text), # Arrow from (0,0) to each point
        arrow = arrow_params,
        color = "darkred",
        size = 0.7, # Size of the arrows
        alpha = 0.7 # Transparency of the arrows
      )
  }
  if (arrow_labels) {
    if (repel) {
      plt <- plt +
        ggrepel::geom_text_repel(
          data = arrow_df,
          aes(x = Comp1, y = Comp2, label = Names),
          max.overlaps = max_overlaps,
          show.legend = FALSE
        )
    } else {
      plt <- plt +
        geom_text(
          data = arrow_df,
          aes(x = Comp1, y = Comp2, label = Names, text = hover_text)
        )
    }
  }

  # if there are marginal plots or covariates and loadings at the same time,
  # don't show secondary axis names to avoid visual mess
  if (
    is.null(marginal_plot) &&
      xor(is.null(scale_factor_load), is.null(scale_factor_covar)) # <=> only loadings or covariates, not both
  ) {
    scale_factor <- c(scale_factor_load, scale_factor_covar) # = the one which is not NULL
    plt <- plt +
      scale_x_continuous(
        sec.axis = sec_axis(
          trans = ~ . / scale_factor,
          name = paste0("Loadings ", dim_names[comp[1]])
        )
      ) +
      scale_y_continuous(
        sec.axis = sec_axis(
          trans = ~ . / scale_factor,
          name = paste0("Loadings ", dim_names[comp[2]])
        )
      )
  }

  plt
}

#' Plot a Beta-Diversity Ordination
#'
#' Scatter plot of ordination sample scores (`beta_dispersion_fit$coords`,
#' as returned by [get_beta_diversity()]), optionally colored/shaped/sized/
#' faceted/animated by `sample_data` columns, annotated with a
#' [stat_beta_diversity()] result, and/or overlaid with a biplot of taxon
#' (`loadings`) and/or predictor (`covariates`) arrows.
#'
#' @param beta_dispersion_fit A list as returned by [get_beta_diversity()].
#' @param hover_variables Character vector of `sample_data` column names to
#'   include in each point's hover text (via the `text` aesthetic, consumed
#'   by `plotly::ggplotly(tooltip = "text")`). Columns not present in
#'   `sample_data` are silently skipped. Default `NULL` (no hover text).
#' @param scaling Which of `beta_dispersion_fit`'s two `vegan` scalings to
#'   plot (`1`: distances between samples are interpretable; `2`: angles
#'   between taxa/predictor arrows are interpretable). Falls back to `1`
#'   with a `warning()` if out of range. Default `1`.
#' @param comp Length-2 integer vector of ordination axes to plot. Falls
#'   back to `c(1, 2)` with a `warning()` if invalid (e.g. an axis beyond
#'   what `beta_dispersion_fit` retained). Default `c(1, 2)`.
#' @param label_name,label_levels `sample_data` column used to color points
#'   (only if non-numeric), and optionally the subset/order of its levels
#'   to keep (samples outside `label_levels` are dropped from the plot; see
#'   `keep_levels()`/`factorize_levels()`, `utils.R`). Overridden by
#'   `stat_beta_dispersion$label_name` when `stat_beta_dispersion` is
#'   supplied. Default `NULL` (no coloring).
#' @param facet_mode `"wrap"` (facet by `facet`, via
#'   [ggplot2::facet_wrap()]) or `"grid"` (facet by `facet_row`/`facet_col`,
#'   via [ggplot2::facet_grid()]). Default `"wrap"`.
#' @param facet,facet_levels Wrap-mode faceting variable (only used if
#'   character/factor) and, optionally, the subset/order of its levels to
#'   keep. Ignored when `facet_mode = "grid"`. Default `NULL` (no facet).
#' @param facet_row,facet_row_levels,facet_col,facet_col_levels Grid-mode
#'   row/column faceting variables and their optional level subset/order.
#'   Ignored when `facet_mode = "wrap"`. Default `NULL` (no facet).
#' @param shape_name,shape_levels `sample_data` column mapped to point shape
#'   (only if character/factor), and its optional level subset/order.
#'   Default `NULL` (no shape mapping).
#' @param size_name `sample_data` column mapped to point size (rescaled to
#'   `[point_size * 0.5, point_size * 3]`). Default `NULL` (fixed size,
#'   `point_size`).
#' @param animation_variable_name,animation_variable_levels `sample_data`
#'   column to embed in the returned plot's data (for a later
#'   [animate_by_variable()] call), and its optional level subset/order.
#'   Default `NULL` (no animation variable embedded).
#' @param remove_na_from_plot Logical; drop samples with `NA` in any of
#'   `label_name`/`shape_name`/`size_name`/`facet`/`facet_row`/`facet_col`/
#'   `animation_variable_name` before plotting. Default `FALSE`.
#' @param ellipses Logical; draw a [ggplot2::stat_ellipse()] per level of
#'   `label_name`, only when `label_name` is discrete. Default `FALSE`.
#' @param fill_ellipses Logical; fill the ellipses with color instead of
#'   just outlining them. Only used when `ellipses = TRUE`. Default
#'   `FALSE`.
#' @param stat_beta_dispersion `NULL` (default) or the result of
#'   [stat_beta_diversity()]: its p-value(s) are embedded in the plot's
#'   subtitle, facet strip labels, or (full-grid mode) per-panel
#'   annotations, and its `label_name` overrides this call's own
#'   `label_name`.
#' @param biplot_loadings,biplot_covariates Logical; draw taxon
#'   (`biplot_loadings`, from `beta_dispersion_fit$loadings`) and/or
#'   predictor (`biplot_covariates`, from `beta_dispersion_fit$covariates`)
#'   arrows. Default `FALSE` (no biplot).
#' @param arrow_labels Logical; label each drawn arrow with its taxon/
#'   predictor name (taxa use `arrow_taxonomy_labels` instead, when given).
#'   Default `FALSE`.
#' @param arrow_taxonomy_labels `NULL` (default; label taxa by their taxon
#'   name) or a character vector of `tax_table` ranks (e.g.
#'   `c("Genus", "Species")`) to build taxon arrow labels from instead
#'   (e.g. `"Escherichia:Coli"`).
#' @param color_arrows_by_taxa Logical; color taxon arrows individually by
#'   `arrow_taxonomy_labels`/taxon name, instead of a single fixed color.
#'   Default `FALSE`.
#' @param arrow_cutoff_load,arrow_cutoff_covar Minimum taxon/predictor arrow
#'   length, normalized to the longest arrow of that type (range `[0, 1]`),
#'   below which the arrow is hidden. Default `0` (show every arrow).
#' @param repel Logical; label arrows with [ggrepel::geom_text_repel()]
#'   instead of [ggplot2::geom_text()]. Incompatible with `plotly`. Default
#'   `FALSE`.
#' @param max_overlaps Passed to [ggrepel::geom_text_repel()]'s
#'   `max.overlaps`, when `repel = TRUE`. Default `10`.
#' @param marginal_plot `NULL` (default; no marginal plot) or a
#'   [ggExtra::ggMarginal()] `type` (e.g. `"boxplot"`, `"density"`),
#'   rendered along both axes. Only applied when `label_name` is discrete,
#'   there's no facet, and `facet_mode = "wrap"`. Incompatible with
#'   `plotly`.
#' @param point_alpha,point_size Alpha/size of sample points. When
#'   `beta_dispersion_fit` used a `fit_filter`, only fit-subset points use
#'   `point_alpha` (the rest use `projected_alpha`); `size_name` overrides
#'   `point_size` with a rescaled range instead of a fixed value. Default
#'   `1`/`3`.
#' @param projected_alpha Alpha of samples that were projected onto (but not
#'   used to fit) the ordination, when `beta_dispersion_fit` used a
#'   `fit_filter`. Unused otherwise. Default `0.4`.
#' @param reverse_dim1,reverse_dim2 Logical; negate the first/second plotted
#'   component (`comp[1]`/`comp[2]`, respectively) of `coords`, `loadings`,
#'   and `covariates` alike — e.g. to match another ordination's arbitrary
#'   axis sign/orientation. Default `FALSE`.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object, or — when `marginal_plot`
#'   is used — the `ggExtraPlot`/`gtable` object returned by
#'   [ggExtra::ggMarginal()] instead.
#'
#' @seealso [get_beta_diversity()], [stat_beta_diversity()],
#'   [beta_diversity()] (a thin wrapper chaining all three),
#'   [animate_by_variable()]
#'
#' @examples
#' data(ps_16s_refinement)
#' bd <- get_beta_diversity(ps_16s_refinement, method = "PCoA", dist = "bray")
#' plot_beta_diversity(bd, label_name = "Protocol", ellipses = TRUE)
#'
#' @export
plot_beta_diversity <- function(
  beta_dispersion_fit,
  hover_variables = NULL,
  scaling = 1,
  comp = c(1, 2),
  label_name = NULL,
  label_levels = NULL,
  facet_mode = "wrap",
  facet = NULL,
  facet_levels = NULL,
  facet_row = NULL,
  facet_row_levels = NULL,
  facet_col = NULL,
  facet_col_levels = NULL,
  shape_name = NULL,
  shape_levels = NULL,
  size_name = NULL,
  animation_variable_name = NULL,
  animation_variable_levels = NULL,
  remove_na_from_plot = FALSE,
  ellipses = FALSE,
  fill_ellipses = FALSE,
  stat_beta_dispersion = NULL,
  biplot_loadings = FALSE,
  biplot_covariates = FALSE,
  arrow_labels = FALSE,
  arrow_taxonomy_labels = NULL,
  color_arrows_by_taxa = FALSE,
  arrow_cutoff_load = 0,
  arrow_cutoff_covar = 0,
  repel = FALSE,
  max_overlaps = 10,
  marginal_plot = NULL,
  point_alpha = 1,
  point_size = 3,
  projected_alpha = 0.4,
  reverse_dim1 = FALSE,
  reverse_dim2 = FALSE
) {
  sample_data <- beta_dispersion_fit$sample_data
  grid_mode <- facet_mode == "grid"

  # for arrow hover info/taxonomy labels only:
  tax_table <- beta_dispersion_fit$tax_table

  if (is.null(ellipses)) {
    ellipses <- FALSE
  }

  # NOTE: species scores indicated by `vegan` = loadings here
  #       covariates (predictors) scores = covariates (constrained models only)
  # -> coordinates of arrows for biplot

  # NOTE: Comp (component) and Dim (dimension) are used interchangeably here
  if (is.null(biplot_loadings)) {
    biplot_loadings <- FALSE
  }
  if (is.null(biplot_covariates)) {
    biplot_covariates <- FALSE
  }
  biplot <- biplot_loadings || biplot_covariates

  axes <- .plot_beta_diversity_select_scaling(
    beta_dispersion_fit,
    scaling,
    comp,
    reverse_dim1,
    reverse_dim2
  )
  coords <- axes$coords
  loadings <- axes$loadings
  covariates <- axes$covariates
  scaling <- axes$scaling
  comp <- axes$comp
  dim_names <- axes$dim_names

  if (!is.null(beta_dispersion_fit$eigen_values)) {
    prop_var_explained <- beta_dispersion_fit$eigen_values /
      sum(beta_dispersion_fit$eigen_values[
        beta_dispersion_fit$eigen_values > 0
      ]) *
      100
  } else {
    prop_var_explained <- NULL
  }

  filtered <- .plot_beta_diversity_filter_levels(
    sample_data,
    coords,
    label_name,
    label_levels
  )
  sample_data <- filtered$sample_data
  coords <- filtered$coords

  filtered <- .plot_beta_diversity_filter_levels(
    sample_data,
    coords,
    shape_name,
    shape_levels
  )
  sample_data <- filtered$sample_data
  coords <- filtered$coords

  filtered <- .plot_beta_diversity_filter_levels(
    sample_data,
    coords,
    facet,
    facet_levels,
    auto_factor = TRUE
  )
  sample_data <- filtered$sample_data
  coords <- filtered$coords

  if (grid_mode) {
    filtered <- .plot_beta_diversity_filter_levels(
      sample_data,
      coords,
      facet_row,
      facet_row_levels,
      auto_factor = TRUE
    )
    sample_data <- filtered$sample_data
    coords <- filtered$coords

    filtered <- .plot_beta_diversity_filter_levels(
      sample_data,
      coords,
      facet_col,
      facet_col_levels,
      auto_factor = TRUE
    )
    sample_data <- filtered$sample_data
    coords <- filtered$coords
  }

  filtered <- .plot_beta_diversity_filter_levels(
    sample_data,
    coords,
    animation_variable_name,
    animation_variable_levels
  )
  sample_data <- filtered$sample_data
  coords <- filtered$coords

  pvalues <- .plot_beta_diversity_pvalues(
    stat_beta_dispersion,
    sample_data,
    grid_mode,
    facet,
    facet_row,
    facet_col
  )
  p_value <- pvalues$p_value
  p_value_df <- pvalues$p_value_df

  # Remove all NA's from plot data (labels, facets or shape) by removing samples
  # having NA for at least one of graphical parameters (label, shape, facet, animation)
  # remove these samples from sample data AND from corresponding coordinates!
  if (remove_na_from_plot) {
    samples_wo_na <- rep(TRUE, nrow(sample_data))
    for (var_name in c(
      label_name,
      shape_name,
      size_name,
      facet,
      facet_row,
      facet_col,
      animation_variable_name
    )) {
      if (!is.null(sample_data[[var_name]])) {
        samples_wo_na <- samples_wo_na & !is.na(sample_data[[var_name]])
      }
    }
    sample_data <- sample_data[samples_wo_na, ]
    coords <- coords[samples_wo_na, ]
  }

  # If stats on group are provided, overwrite the label_name by the one of this group
  if (!is.null(stat_beta_dispersion$label_name)) {
    label_name <- stat_beta_dispersion$label_name
  }

  label <- if (!is.null(label_name)) sample_data[[label_name]] else NULL
  animation_variable <- if (!is.null(animation_variable_name)) {
    sample_data[[animation_variable_name]]
  } else {
    NULL
  }

  if (is.character(label)) {
    label <- factor(label, levels = unique(label))
  }

  built <- .plot_beta_diversity_build_df(
    sample_data,
    coords,
    comp,
    label,
    animation_variable_name,
    animation_variable,
    grid_mode,
    facet,
    facet_row,
    facet_col,
    p_value,
    shape_name,
    size_name,
    hover_variables
  )
  plot_df <- built$plot_df
  shape <- built$shape
  size <- built$size

  plt <- .plot_beta_diversity_base_plot(
    plot_df,
    label,
    shape,
    size,
    label_name,
    shape_name,
    size_name,
    point_alpha,
    point_size,
    projected_alpha,
    ellipses,
    fill_ellipses
  )

  plt <- .plot_beta_diversity_labs(
    plt,
    beta_dispersion_fit,
    dim_names,
    comp,
    prop_var_explained,
    remove_na_from_plot,
    scaling,
    biplot,
    biplot_loadings,
    loadings,
    biplot_covariates,
    covariates,
    arrow_cutoff_load,
    arrow_cutoff_covar,
    p_value
  )

  if (biplot && (!is.null(loadings) || !is.null(covariates))) {
    arrow_data <- .plot_beta_diversity_arrow_data(
      loadings,
      covariates,
      comp,
      coords,
      biplot_loadings,
      biplot_covariates,
      arrow_cutoff_load,
      arrow_cutoff_covar,
      tax_table,
      arrow_taxonomy_labels
    )

    # If there are loadings and/or covariates, plot arrows
    if (nrow(arrow_data$arrow_df) > 0) {
      plt <- .plot_beta_diversity_arrow_plot(
        plt,
        arrow_data$arrow_df,
        arrow_data$scale_factor_load,
        arrow_data$scale_factor_covar,
        color_arrows_by_taxa,
        arrow_labels,
        repel,
        max_overlaps,
        marginal_plot,
        dim_names,
        comp
      )
    }
  }

  plt <- .plot_beta_diversity_facets(
    plt,
    plot_df,
    grid_mode,
    facet,
    facet_row,
    facet_col,
    p_value_df
  )

  # Marginal plot in case of factor variable - boxplot, density...
  if (
    !is.null(marginal_plot) && is.factor(label) && is.null(facet) && !grid_mode
  ) {
    plt <- plt +
      theme(
        legend.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "left"
      ) # place legend to the left so that it doesn't interfere with the marginal plot
    plt <- ggExtra::ggMarginal(
      plt,
      type = marginal_plot,
      groupColour = TRUE,
      groupFill = TRUE
    )
  } else if (!is.null(label)) {
    plt <- plt +
      theme(legend.title = element_text(face = "bold", hjust = 0.5))
  }

  return(plt)
}


# Inspired from animateDimReduction() in Feature Selector
#' Animate a Ggplot Object by a Variable
#'
#' Adds a [gganimate::transition_states()] animation over the levels of
#' `animation_variable_name` to an existing `ggplot` object (typically the
#' output of [plot_beta_diversity()]) and renders it to a GIF via
#' [gifski::gifski()].
#'
#' @param ggplot_obj A `ggplot` object to animate.
#' @param animation_variable_name Name of the (typically discrete) variable
#'   in `ggplot_obj`'s underlying data to animate over; one animation frame
#'   group per level.
#' @param return_anim Logical; return the rendered animation object. Default
#'   `TRUE`.
#' @param save_path `NULL` (default; don't save) or a file path to save the
#'   rendered animation to (e.g. a `.gif` path), via
#'   [gganimate::anim_save()].
#' @param nframes Total number of animation frames. Default `150`.
#' @param fps Frames per second. Default `25`.
#' @param width,height Rendered animation size in pixels. Default `900`/
#'   `700`.
#' @param res Rendered animation resolution (dpi). Default `200`.
#'
#' @return The rendered `gganim` animation object if `return_anim = TRUE`,
#'   otherwise `NULL`. Either way, the animation is written to `save_path`
#'   first (as a side effect) if one was given.
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(ps_16s_refinement)
#' # animation_variable_name must also be passed to plot_beta_diversity()
#' # (via `...`) so the variable is embedded in the plot's data first.
#' bd <- beta_diversity(
#'   ps_16s_refinement,
#'   group = "Protocol",
#'   animation_variable_name = "Protocol"
#' )
#' animate_by_variable(bd, "Protocol", nframes = 10, fps = 5)
#' }
animate_by_variable <- function(
  ggplot_obj,
  animation_variable_name,
  return_anim = TRUE,
  save_path = NULL,
  nframes = 150,
  fps = 25,
  width = 900,
  height = 700,
  res = 200
) {
  plt <-
    ggplot_obj +
    gganimate::transition_states(
      get(animation_variable_name),
      transition_length = 2,
      state_length = 3
    ) +
    labs(subtitle = paste0(animation_variable_name, " = {closest_state}")) +
    ease_aes('linear')
  anim <-
    animate(
      plt,
      nframes = nframes,
      fps = fps,
      width = width,
      height = height,
      res = res,
      renderer = gifski_renderer(loop = TRUE)
    )

  if (!is.null(save_path)) {
    anim_save(save_path, anim)
  }

  if (return_anim) {
    return(anim)
  } else {
    return(NULL)
  }
}
