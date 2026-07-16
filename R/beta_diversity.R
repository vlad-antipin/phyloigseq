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
      # TODO: consider routing through sparse_distance() once standardization can be
      # done sparsely, to avoid materializing the full dense OTU matrix here.
      if (!is.null(dist) && dist %in% c("euclidean", "manhattan")) {
        otu_mat <- as(otu_table(physeq), "matrix")
        dist_matrix <- vegdist(
          decostand(
            otu_mat[, apply(otu_mat, 2, var, na.rm = TRUE) > 0],
            method = "standardize"
          ),
          method = dist
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
    fit <- vegan::metaMDS(otu_fit, distance = dist)
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
  if (!is.null(confounders)) {
    model <- paste0(
      model,
      " + ",
      "Condition(",
      paste(paste0("`", confounders, "`"), collapse = " + "),
      ")"
    )
  }

  formula <- as.formula(paste(
    "as(otu_table(reverseASV(physeq_fit)), 'matrix')",
    "~",
    model
  ))
  df_fit <- as(sample_data(physeq_fit), "data.frame")

  if (method == "cca") {
    fit <- vegan::cca(formula, data = df_fit, na.action = na.exclude)
  } else if (method == "rda") {
    fit <- vegan::rda(formula, data = df_fit, na.action = na.exclude)
  } else if (method == "dbrda") {
    fit <- vegan::capscale(
      formula,
      data = df_fit,
      na.action = na.exclude,
      dist = dist
    )
  } else {
    stop("Invalid method")
  }

  eigen_values <- vegan::eigenvals(fit)
  n_axes <- min(ndims, length(eigen_values))

  # For RDA/CCA project all samples using WA scores; dbRDA has no standard projection.
  # capscale (dbRDA) inherits from "rda"/"cca" so check method name, not class.
  get_constrained_coords <- function(scaling) {
    if (!is.null(fit_sample_names) && method %in% c("rda", "cca")) {
      otu_all <- as(otu_table(reverseASV(physeq)), "matrix")
      # predict(type="wa") defaults to model="CCA" — constrained axes only.
      # Explicitly fetch residual (CA) axes too and combine up to n_axes.
      wa_cca <- predict(fit, newdata = otu_all, type = "wa", model = "CCA")
      n_cca <- ncol(wa_cca)
      n_ca_need <- max(0, n_axes - n_cca)
      if (n_ca_need > 0 && !is.null(fit$CA) && fit$CA$rank > 0) {
        wa_ca <- predict(fit, newdata = otu_all, type = "wa", model = "CA")
        wa_ca <- wa_ca[, 1:min(n_ca_need, ncol(wa_ca)), drop = FALSE]
        return(cbind(wa_cca, wa_ca))
      }
      return(wa_cca[, 1:min(n_axes, n_cca), drop = FALSE])
    }
    if (!is.null(fit_sample_names)) {
      sc <- scores(fit, display = "sites", choices = 1:n_axes, scaling = scaling)
      return(.warn_and_na_fill("dbRDA", sample_names(physeq), fit_sample_names, sc))
    }
    scores(fit, display = "sites", choices = 1:n_axes, scaling = scaling)
  }

  coords <- list(get_constrained_coords(1), get_constrained_coords(2))
  loadings <- list(
    scores(fit, display = "species", choices = 1:n_axes, scaling = 1),
    scores(fit, display = "species", choices = 1:n_axes, scaling = 2)
  )
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

  if (!is(physeq, "phyloseq")) {
    stop("Need a phyloseq or PhyloIgSeq object")
  }

  if (!is.null(fraction_id_name) && !is.null(fraction_ids)) {
    physeq <- prune_samples(
      sample_data(physeq)[[fraction_id_name]] %in% fraction_ids,
      physeq
    )
  }

  # Force taxa to be columns of otu table - like in internal phyloseq function phyloseq:::veganifyOTU()
  if (taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }

  # NOTE: agglomerate taxa BEFORE transforming the data
  if (!is.null(taxrank)) {
    physeq <- tax_glom(physeq = physeq, taxrank = taxrank)
    taxa_names(physeq) <- make.unique(tax_table(physeq)[, taxrank]) # NOTE: names of taxons are not necesserily unique!
  }

  if (!is.null(transform_abundances) && transform_abundances != "identity") {
    physeq <- microbiome::transform(
      physeq,
      transform = transform_abundances,
      target = "OTU", # TODO: and still, clr will scale over samples
      shift = 0, # pseudocount added (shifts baseline)
      scale = 1, # if transform is "scale"
      log10 = TRUE,
      reference = 1
    )
  }

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
  } else if (!is.null(model)) {
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

#' Permutation Test for Beta-Diversity Ordination
#' @export
stat_beta_diversity <- function(
  beta.dispersion.fit,
  facet.name = NULL, # backward compat: treated as facet in wrap mode
  facet.mode = "wrap", # "grid", "wrap"
  facet = NULL, # wrap-mode facet variable name
  facet.row = NULL, # grid-mode row facet variable name
  facet.col = NULL, # grid-mode col facet variable name
  comp = c(1, 2), # if NULL, consider all components
  strata.name = NULL, # in case of paired design
  label.name,
  pairwise = FALSE
) {
  if (!is.null(facet.name) && is.null(facet)) {
    facet <- facet.name
  }
  if (is.null(facet.mode)) {
    facet.mode <- "wrap"
  }

  full.grid.mode <- facet.mode == "grid" &&
    !is.null(facet.row) &&
    !is.null(facet.col)
  active.facet <- if (facet.mode == "wrap") {
    facet
  } else if (!is.null(facet.row) && is.null(facet.col)) {
    facet.row
  } else if (is.null(facet.row) && !is.null(facet.col)) {
    facet.col
  } else {
    NULL
  }

  if (!label.name %in% colnames(beta.dispersion.fit$sample_data)) {
    warning("Wrong label name")
    return(NULL)
  }

  if (!is.null(comp)) {
    stat.data.full <- data.frame(
      Comp1 = beta.dispersion.fit$coords[[1]][, comp[1]],
      Comp2 = beta.dispersion.fit$coords[[1]][, comp[2]],
      Label = beta.dispersion.fit$sample_data[[label.name]]
    )
  } else {
    stat.data.full <- data.frame(
      beta.dispersion.fit$coords[[1]],
      Label = beta.dispersion.fit$sample_data[[label.name]]
    )
  }

  if (!is.null(strata.name)) {
    stat.data.full[[strata.name]] <- beta.dispersion.fit$sample_data[[
      strata.name
    ]]
  }

  # Embed facet variables into stat.data.full before na.omit so subsetting is consistent
  if (
    full.grid.mode &&
      !is.null(facet.row) &&
      facet.row %in% colnames(beta.dispersion.fit$sample_data)
  ) {
    stat.data.full[[".facet.row"]] <- beta.dispersion.fit$sample_data[[
      facet.row
    ]]
  }
  if (
    full.grid.mode &&
      !is.null(facet.col) &&
      facet.col %in% colnames(beta.dispersion.fit$sample_data)
  ) {
    stat.data.full[[".facet.col"]] <- beta.dispersion.fit$sample_data[[
      facet.col
    ]]
  }
  if (
    !full.grid.mode &&
      !is.null(active.facet) &&
      active.facet %in% colnames(beta.dispersion.fit$sample_data)
  ) {
    stat.data.full[[".facet.row"]] <- beta.dispersion.fit$sample_data[[
      active.facet
    ]]
  }

  stat.data.full <- na.omit(stat.data.full)

  # Shared helper: run PERMANOVA or envfit on a single data slice
  run_test <- function(stat.data) {
    stat.data[[".facet.row"]] <- NULL
    stat.data[[".facet.col"]] <- NULL
    if (is.factor(stat.data$Label) || is.character(stat.data$Label)) {
      stat.data <- na.omit(stat.data)
      if (!is.null(strata.name)) {
        strata <- stat.data[[strata.name]]
        stat.data[[strata.name]] <- NULL
      } else {
        strata <- NULL
      }
      test.result <- vegan::adonis2(
        stat.data[, colnames(stat.data) != "Label"] ~ Label,
        data = stat.data,
        strata = strata,
        permutations = 999,
        method = "euclidean"
      )
      p.raw <- test.result$`Pr(>F)`[1]
      p.text <- if (!is.null(strata)) {
        paste0(
          "restricted by ",
          strata.name,
          " PERMANOVA p = ",
          format(signif(p.raw, digits = 2), scientific = TRUE)
        )
      } else {
        paste0(
          "PERMANOVA p = ",
          format(signif(p.raw, digits = 2), scientific = TRUE)
        )
      }
      list(
        test.result = test.result,
        p.raw = p.raw,
        p.text = p.text,
        stat = "permanova"
      )
    } else {
      test.result <- vegan::envfit(
        stat.data[, colnames(stat.data) != "Label"] ~ Label,
        data = stat.data,
        permutations = 999
      )
      p.raw <- test.result$vectors$pvals[1]
      p.text <- paste0(
        "Correlation p = ",
        format(signif(p.raw, digits = 2), scientific = TRUE)
      )
      list(
        test.result = test.result,
        p.raw = p.raw,
        p.text = p.text,
        stat = "envfit"
      )
    }
  }

  run_pairwise <- function(stat.data) {
    stat.data[[".facet.row"]] <- NULL
    stat.data[[".facet.col"]] <- NULL
    if (!is.factor(stat.data$Label) && !is.character(stat.data$Label)) {
      return(NULL)
    }
    labels <- unique(na.omit(as.character(stat.data$Label)))
    if (length(labels) < 2) {
      return(NULL)
    }
    pairs <- combn(labels, 2, simplify = FALSE)
    pair_results <- lapply(pairs, function(pair) {
      sub <- stat.data[as.character(stat.data$Label) %in% pair, ]
      sub <- na.omit(sub)
      if (length(unique(as.character(sub$Label))) < 2) {
        return(NULL)
      }
      sub_strata <- NULL
      if (!is.null(strata.name) && strata.name %in% colnames(sub)) {
        sub_strata <- sub[[strata.name]]
        sub[[strata.name]] <- NULL
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
        p.raw = res$`Pr(>F)`[1],
        stringsAsFactors = FALSE
      )
    })
    pair_results <- do.call(rbind, Filter(Negate(is.null), pair_results))
    if (is.null(pair_results) || nrow(pair_results) == 0) {
      return(NULL)
    }
    pair_results$p.adj <- p.adjust(pair_results$p.raw, method = "BH")
    pair_results
  }

  dim_used <- if (!is.null(comp)) {
    comp
  } else {
    paste("all", ncol(beta.dispersion.fit$coords))
  }

  if (full.grid.mode) {
    row.vals <- unique(stat.data.full[[".facet.row"]])
    col.vals <- unique(stat.data.full[[".facet.col"]])
    combinations <- expand.grid(
      row = row.vals,
      col = col.vals,
      stringsAsFactors = FALSE
    )

    test.result.all <- list()
    p.value.df <- data.frame(
      facet.row = character(0),
      facet.col = character(0),
      p.label = character(0),
      p.raw = numeric(0),
      stringsAsFactors = FALSE
    )
    stat.type <- NULL
    pairwise.list <- list()

    for (i in seq_len(nrow(combinations))) {
      row.val <- combinations$row[i]
      col.val <- combinations$col[i]
      stat.data <- stat.data.full[
        stat.data.full[[".facet.row"]] == row.val &
          stat.data.full[[".facet.col"]] == col.val,
      ]
      if (length(unique(stat.data$Label)) < 2) {
        warning(paste0(
          "Facet [",
          row.val,
          ", ",
          col.val,
          "]: label doesn't vary, skipping."
        ))
        next
      }
      res <- run_test(stat.data)
      key <- paste(row.val, col.val, sep = "//")
      test.result.all[[key]] <- res$test.result
      stat.type <- res$stat
      p.value.df <- rbind(
        p.value.df,
        data.frame(
          facet.row = as.character(row.val),
          facet.col = as.character(col.val),
          p.label = res$p.text,
          p.raw = res$p.raw,
          stringsAsFactors = FALSE
        )
      )
      if (pairwise && res$stat == "permanova") {
        pw <- run_pairwise(stat.data)
        if (!is.null(pw)) {
          if (!is.null(facet.row)) {
            pw[[facet.row]] <- as.character(row.val)
          }
          if (!is.null(facet.col)) {
            pw[[facet.col]] <- as.character(col.val)
          }
          pairwise.list[[length(pairwise.list) + 1]] <- pw
        }
      }
    }

    return(list(
      stat = stat.type,
      label.name = label.name,
      facet.name = facet.row,
      facet.mode = facet.mode,
      facet = NULL,
      facet.row = facet.row,
      facet.col = facet.col,
      dim_used = dim_used,
      test.result = if (length(test.result.all) > 0) test.result.all else NULL,
      p.value = NULL,
      p.value.df = p.value.df,
      p.value.raw = p.value.df$p.raw,
      pairwise.df = if (pairwise && length(pairwise.list) > 0) {
        do.call(rbind, pairwise.list)
      } else {
        NULL
      }
    ))
  }

  # Single-facet or no-facet mode
  facet.levels.iter <- if (".facet.row" %in% colnames(stat.data.full)) {
    unique(stat.data.full[[".facet.row"]])
  } else {
    NA
  }

  test.result.all <- list()
  p.value.all <- c()
  p.value.raw <- c()
  stat.type <- NULL
  pairwise.list <- list()

  for (facet.val in facet.levels.iter) {
    if (is.na(facet.val)) {
      stat.data <- stat.data.full
      idx <- 1
    } else {
      stat.data <- stat.data.full[stat.data.full[[".facet.row"]] == facet.val, ]
      idx <- facet.val
    }

    if (length(unique(stat.data$Label)) == 1) {
      warning("Variable doesn't vary, no stats performed...")
      return(NULL)
    }

    res <- run_test(stat.data)
    test.result.all[[idx]] <- res$test.result
    stat.type <- res$stat
    p.value.raw[idx] <- res$p.raw
    p.value.all[idx] <- res$p.text
    if (pairwise && res$stat == "permanova") {
      pw <- run_pairwise(stat.data)
      if (!is.null(pw)) {
        if (!is.na(facet.val) && !is.null(active.facet)) {
          pw[[active.facet]] <- as.character(facet.val)
        }
        pairwise.list[[length(pairwise.list) + 1]] <- pw
      }
    }
  }

  return(list(
    stat = stat.type,
    label.name = label.name,
    facet.name = active.facet,
    facet.mode = facet.mode,
    facet = facet,
    facet.row = facet.row,
    facet.col = facet.col,
    dim_used = dim_used,
    test.result = if (length(test.result.all) > 0) test.result.all else NULL,
    p.value = p.value.all,
    p.value.df = NULL,
    p.value.raw = p.value.raw,
    pairwise.df = if (pairwise && length(pairwise.list) > 0) {
      do.call(rbind, pairwise.list)
    } else {
      NULL
    }
  ))
}

#' Get and Plot Beta-Diversity from Phyloseq Object
#' @export
full_beta_diversity <- function(
  physeq,
  taxrank = NULL,
  fraction_id_name = NULL,
  fraction_ids = NULL,
  method = "PCoA",
  model = NULL,
  dist = "bray",
  group = NULL,
  stat = "none",
  species = FALSE,
  axis_x = 1,
  axis_y = 2,
  # nf= 5,
  # type= "boxplot",
  # color_vector= c("cyan4","brown","deepskyblue", "black","red"),
  # legend_title= NULL,
  # lwd=1,
  # conf=0.9,
  # cex=2,
  # font=2,
  # pch=20,
  # draw= "lines",
  # ylimits="auto",
  # xlimits= "auto",
  # text=FALSE,
  # ncol=1,
  # x.intersp = 1,
  # y.intersp=0.5,
  # where="topleft",
  # inset=0.2,
  # pca=TRUE,
  # stat.cex= 2,
  # legend.cex=2,
  # widths= c(1,1),
  # heights= c(1,1),
  # margins= c(1,1,1,1),
  ...
) {
  # NOTE: `species` is no longer forwarded — get_beta_diversity() dropped its
  # own unused `species` param in the 14b cleanup pass. This function's own
  # `species` arg is now fully dead; left for the 14d argument-hygiene pass.
  beta.dispersion.fit <- get_beta_diversity(
    physeq = physeq,
    taxrank = taxrank,
    fraction_id_name = fraction_id_name,
    fraction_ids = fraction_ids,
    dist = dist,
    method = method,
    model = model
  )
  if (!is.null(stat) && stat != FALSE && stat != "none") {
    stat.beta.dispersion <- stat_beta_diversity(
      beta.dispersion.fit,
      comp = c(axis_x, axis_y),
      label.name = group
    )
  } else {
    stat.beta.dispersion <- NULL
  }

  plt <- plot_beta_diversity(
    beta.dispersion.fit = beta.dispersion.fit,
    stat.beta.dispersion = stat.beta.dispersion,
    ...
  )
  print(plt)

  return(list(
    reduction = beta.dispersion.fit$fit,
    stat = stat.beta.dispersion[[1]]$test.result
  ))
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
    theme_minimal() +
    ggplot2::theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.title = element_text(face = "bold", hjust = 0.5)
    )
  return(plt)
}


#' Plot Beta-Diversity Ordination
#' @export
plot_beta_diversity <- function(
  beta.dispersion.fit, # output of get_beta_diversity()
  hover.variables = NULL, # all columns of sample data if NULL
  scaling = 1, # see vegan's scalings 1 and 2
  comp = c(1, 2), # components to plot (loadings might be lacking components beyond 2)
  label.name = NULL,
  label.levels = NULL,
  facet.mode = "wrap", # "grid", "wrap"
  facet = NULL, # wrap-mode facet variable name
  facet.levels = NULL,
  facet.row = NULL, # grid-mode row facet variable name
  facet.row.levels = NULL,
  facet.col = NULL, # grid-mode col facet variable name
  facet.col.levels = NULL,
  shape.name = NULL,
  shape.levels = NULL,
  size.name = NULL,
  animation.variable.name = NULL,
  animation.variable.levels = NULL,
  remove.na.from.plot = FALSE,
  ellipses = FALSE, # only for factor lable
  fill.ellipses = FALSE, # fill ellipses with color (requires ellipses = TRUE)
  stat.beta.dispersion = NULL, # in case of stats performed between groups
  # can contain  label.name, overwriting the parameter `label.name`
  biplot.loadings = FALSE, # aka species
  biplot.covariates = FALSE, # aka environemental variables

  arrow.labels = FALSE, # display species names near the arrows
  arrow.taxonomy.labels = NULL, # taxranks to display as text label for
  # each species arrow e.g. c("Genus", "Species") will give
  # Escherichia : Coli
  color.arrows.by.taxa = FALSE, # color arrows according to taxonomy labels
  arrow.cutoff.load = 0, # [0,1] arrow length normalized to the longest arrow
  arrow.cutoff.covar = 0, # to hide shorter arrows
  repel = FALSE, # repel arrow text labels
  # NOTE: repel is incompatible with plotly!
  max.overlaps = 10, # in case of repel
  marginal.plot = NULL, # only for factor label
  # NOTE: marginal plot is incompatible with plotly!
  raw.loadings = TRUE, # TODO: think if you should implement correlations here as an alternative
  point.alpha = 1,
  point.size = 3,
  projected.alpha = 0.4, # alpha for samples not in the fit subset (when fit.filter is used)
  reverse.dim1 = FALSE,
  reverse.dim2 = FALSE
) {
  sample.data <- beta.dispersion.fit$sample_data

  grid.mode <- facet.mode == "grid"

  # for the hover info only:
  tax.table <- beta.dispersion.fit$tax_table
  taxrank <- beta.dispersion.fit$taxrank

  if (is.null(ellipses)) {
    ellipses <- FALSE
  }

  # NOTE: species scores indicated by `vegan` = loadings here
  #       covariates (predictors) scores = covariates (constrained models only)
  # -> coordinates of arrows for biplot

  # NOTE: Comp (component) and Dim (dimension) are used interchangeably here
  if (is.null(biplot.loadings)) {
    biplot.loadings <- FALSE
  }
  if (is.null(biplot.covariates)) {
    biplot.covariates <- FALSE
  }
  biplot <- biplot.loadings | biplot.covariates

  # If no coordinates provided
  if (length(beta.dispersion.fit$coords) == 0) {
    coords <- NULL
    # If only one scaling provided
  } else if (length(beta.dispersion.fit$coords) == 1) {
    coords <- beta.dispersion.fit$coords[[1]]
    # Select the scaling
  } else if (0 < scaling & scaling <= length(beta.dispersion.fit$coords)) {
    coords <- beta.dispersion.fit$coords[[scaling]]
  } else {
    warning("Wrong scaling, scaling 1 is used")
    coords <- beta.dispersion.fit$coords[[1]]
    scaling <- 1
  }

  # Same as for coords
  if (length(beta.dispersion.fit$loadings) == 0) {
    loadings <- NULL
  } else if (length(beta.dispersion.fit$loadings) == 1) {
    loadings <- beta.dispersion.fit$loadings[[1]]
  } else if (0 < scaling & scaling <= length(beta.dispersion.fit$loadings)) {
    loadings <- beta.dispersion.fit$loadings[[scaling]]
  } else {
    warning("Wrong scaling, scaling 1 is used")
    loadings <- beta.dispersion.fit$loadings[[1]]
    scaling <- 1
  }

  # Same as for coords
  if (length(beta.dispersion.fit$covariates) == 0) {
    covariates <- NULL
  } else if (length(beta.dispersion.fit$covariates) == 1) {
    covariates <- beta.dispersion.fit$covariates[[1]]
  } else if (0 < scaling & scaling <= length(beta.dispersion.fit$covariates)) {
    covariates <- beta.dispersion.fit$covariates[[scaling]]
  } else {
    warning("Wrong scaling, scaling 1 is used")
    covariates <- beta.dispersion.fit$covariates[[1]]
    scaling <- 1
  }

  if (!is.null(beta.dispersion.fit$eigen_values)) {
    prop.var.explained <- beta.dispersion.fit$eigen_values /
      sum(beta.dispersion.fit$eigen_values[
        beta.dispersion.fit$eigen_values > 0
      ]) *
      100
  } else {
    prop.var.explained <- NULL
  }

  dist <- beta.dispersion.fit$dist
  method <- beta.dispersion.fit$method
  model <- beta.dispersion.fit$model

  if (
    !is.null(label.name) &&
      !is.null(label.levels) &&
      !is.numeric(sample.data[[label.name]])
  ) {
    keep <- keep_levels(sample.data[[label.name]], label.levels)
    sample.data <- sample.data[keep, , drop = FALSE]
    if (!is.null(coords)) {
      coords <- coords[keep, , drop = FALSE]
    }
    sample.data[[label.name]] <- factorize_levels(
      sample.data[[label.name]],
      label.levels
    )
  }

  if (
    !is.null(shape.name) &&
      !is.null(shape.levels) &&
      !is.numeric(sample.data[[shape.name]])
  ) {
    keep <- keep_levels(sample.data[[shape.name]], shape.levels)
    sample.data <- sample.data[keep, , drop = FALSE]
    if (!is.null(coords)) {
      coords <- coords[keep, , drop = FALSE]
    }
    sample.data[[shape.name]] <- factorize_levels(
      sample.data[[shape.name]],
      shape.levels
    )
  }

  if (!is.null(facet) && !is.numeric(sample.data[[facet]])) {
    if (!is.null(facet.levels)) {
      keep <- keep_levels(sample.data[[facet]], facet.levels)
      sample.data <- sample.data[keep, , drop = FALSE]
      if (!is.null(coords)) {
        coords <- coords[keep, , drop = FALSE]
      }
      sample.data[[facet]] <- factorize_levels(
        sample.data[[facet]],
        facet.levels
      )
    } else {
      sample.data[[facet]] <- as.factor(sample.data[[facet]])
    }
  }

  if (
    grid.mode && !is.null(facet.row) && !is.numeric(sample.data[[facet.row]])
  ) {
    if (!is.null(facet.row.levels)) {
      keep <- keep_levels(sample.data[[facet.row]], facet.row.levels)
      sample.data <- sample.data[keep, , drop = FALSE]
      if (!is.null(coords)) {
        coords <- coords[keep, , drop = FALSE]
      }
      sample.data[[facet.row]] <- factorize_levels(
        sample.data[[facet.row]],
        facet.row.levels
      )
    } else {
      sample.data[[facet.row]] <- as.factor(sample.data[[facet.row]])
    }
  }

  if (
    grid.mode && !is.null(facet.col) && !is.numeric(sample.data[[facet.col]])
  ) {
    if (!is.null(facet.col.levels)) {
      keep <- keep_levels(sample.data[[facet.col]], facet.col.levels)
      sample.data <- sample.data[keep, , drop = FALSE]
      if (!is.null(coords)) {
        coords <- coords[keep, , drop = FALSE]
      }
      sample.data[[facet.col]] <- factorize_levels(
        sample.data[[facet.col]],
        facet.col.levels
      )
    } else {
      sample.data[[facet.col]] <- as.factor(sample.data[[facet.col]])
    }
  }

  if (
    !is.null(animation.variable.name) &&
      !is.null(animation.variable.levels) &&
      !is.numeric(sample.data[[animation.variable.name]])
  ) {
    keep <- keep_levels(
      sample.data[[animation.variable.name]],
      animation.variable.levels
    )
    sample.data <- sample.data[keep, , drop = FALSE]
    if (!is.null(coords)) {
      coords <- coords[keep, , drop = FALSE]
    }
    sample.data[[animation.variable.name]] <- factorize_levels(
      sample.data[[animation.variable.name]],
      animation.variable.levels
    )
  }
  # (split by newline character to get one line of text for the subtitle)
  # if no label name is provided, stats are ignored
  p.value <- NULL
  p.value.df <- NULL

  if (!is.null(stat.beta.dispersion$label.name)) {
    if (
      grid.mode &&
        !is.null(facet.row) &&
        !is.null(facet.col) &&
        !is.null(stat.beta.dispersion$p.value.df)
    ) {
      # Full grid mode: p-values as data frame, rendered as geom_text annotations
      p.value.df <- stat.beta.dispersion$p.value.df
    } else if (!is.null(stat.beta.dispersion$p.value)) {
      active.facet <- if (grid.mode) facet.row %||% facet.col else facet
      if (
        is.null(active.facet) && is.null(names(stat.beta.dispersion$p.value))
      ) {
        # No facets: embed single p-value in subtitle
        p.value <- paste(
          strsplit(stat.beta.dispersion$p.value[[1]], "\n")[[1]],
          collapse = " "
        )
      } else if (
        !is.null(active.facet) &&
          !is.null(names(stat.beta.dispersion$p.value)) &&
          all(
            levels(sample.data[[active.facet]]) %in%
              names(stat.beta.dispersion$p.value)
          )
      ) {
        # Single-facet: embed per-facet p-value into facet strip labels
        p.value <- stat.beta.dispersion$p.value[levels(sample.data[[
          active.facet
        ]])]
      }
    }
  }

  # Remove all NA's from plot data (labels, facets or shape) by removing samples
  # having NA for at least one of graphical parameters (label, shape, facet, animation)
  # remove these samples from sample data AND from corresponding coordinates!
  if (remove.na.from.plot) {
    samples.wo.na <- rep(TRUE, nrow(sample.data))
    for (var.name in c(
      label.name,
      shape.name,
      size.name,
      facet,
      facet.row,
      facet.col,
      animation.variable.name
    )) {
      if (!is.null(sample.data[[var.name]])) {
        samples.wo.na <- samples.wo.na & !is.na(sample.data[[var.name]])
      }
    }
    sample.data <- sample.data[samples.wo.na, ]
    coords <- coords[samples.wo.na, ]
  }

  # If stats on group are provided, overwrite the label.name by the one of this group
  if (!is.null(stat.beta.dispersion$label.name)) {
    label.name <- stat.beta.dispersion$label.name
  }

  if (!is.null(label.name)) {
    label <- sample.data[[label.name]]
  } else {
    label <- NULL
  }

  if (!is.null(animation.variable.name)) {
    animation.variable <- sample.data[[animation.variable.name]]
  } else {
    animation.variable <- NULL
  }

  if (is.character(label)) {
    label <- factor(label, levels = unique(label))
  }

  if (!is.null(colnames(coords))) {
    dim.names <- colnames(coords)
  } else {
    dim.names <- paste0("Dim ", 1:ncol(coords))
  }

  # Check if components are valid
  if (!is.numeric(comp) | length(comp) != 2 | !all(comp %in% 1:ncol(coords))) {
    warning("Wrong components, forced to first two components")
    comp <- c(1, 2)
  }

  # Check if there are loadings and covariates for given components
  if (!is.null(loadings)) {
    if (!all(comp %in% 1:ncol(loadings))) {
      warning("No species score provided for these components")
      loadings <- NULL
    }
  }

  if (!is.null(covariates)) {
    if (!all(comp %in% 1:ncol(covariates))) {
      warning("No species score provided for these components")
      covariates <- NULL
    }
  }

  if (reverse.dim1 && !is.null(coords)) {
    coords[, comp[1]] <- -coords[, comp[1]]
  }
  if (reverse.dim2 && !is.null(coords)) {
    coords[, comp[2]] <- -coords[, comp[2]]
  }
  if (reverse.dim1 && !is.null(loadings)) {
    loadings[, comp[1]] <- -loadings[, comp[1]]
  }
  if (reverse.dim2 && !is.null(loadings)) {
    loadings[, comp[2]] <- -loadings[, comp[2]]
  }
  if (reverse.dim1 && !is.null(covariates)) {
    covariates[, comp[1]] <- -covariates[, comp[1]]
  }
  if (reverse.dim2 && !is.null(covariates)) {
    covariates[, comp[2]] <- -covariates[, comp[2]]
  }

  plot.df <- data.frame(Comp1 = coords[, comp[1]], Comp2 = coords[, comp[2]])

  if (!is.null(label)) {
    plot.df$label <- label
  }

  if (!is.null(animation.variable)) {
    plot.df[[animation.variable.name]] <- animation.variable
  }

  if (grid.mode) {
    if (!is.null(facet.row)) {
      facet.row.data <- sample.data[[facet.row]]
      if (!is.factor(facet.row.data)) {
        facet.row.data <- factor(facet.row.data)
      }
      if (is.null(facet.col) && length(p.value) > 1) {
        plot.df$facet.row <- factor(
          paste0(
            facet.row,
            " = ",
            as.character(facet.row.data),
            "\n",
            p.value[as.character(facet.row.data)]
          ),
          levels = paste0(
            facet.row,
            " = ",
            as.character(levels(facet.row.data)),
            "\n",
            p.value[as.character(levels(facet.row.data))]
          )
        )
      } else {
        plot.df$facet.row <- factor(
          paste0(facet.row, " = ", as.character(facet.row.data)),
          levels = paste0(
            facet.row,
            " = ",
            as.character(levels(facet.row.data))
          )
        )
      }
    }
    if (!is.null(facet.col)) {
      facet.col.data <- sample.data[[facet.col]]
      if (!is.factor(facet.col.data)) {
        facet.col.data <- factor(facet.col.data)
      }
      if (is.null(facet.row) && length(p.value) > 1) {
        plot.df$facet.col <- factor(
          paste0(
            facet.col,
            " = ",
            as.character(facet.col.data),
            "\n",
            p.value[as.character(facet.col.data)]
          ),
          levels = paste0(
            facet.col,
            " = ",
            as.character(levels(facet.col.data)),
            "\n",
            p.value[as.character(levels(facet.col.data))]
          )
        )
      } else {
        plot.df$facet.col <- factor(
          paste0(facet.col, " = ", as.character(facet.col.data)),
          levels = paste0(
            facet.col,
            " = ",
            as.character(levels(facet.col.data))
          )
        )
      }
    }
  } else if (!is.null(facet)) {
    facet.data <- sample.data[[facet]]
    if (is.character(facet.data) | is.factor(facet.data)) {
      if (!is.factor(facet.data)) {
        facet.data <- as.factor(facet.data)
      }
      plot.df$facet <- factor(
        paste0(
          facet,
          " = ",
          as.character(facet.data),
          if (length(p.value) > 1) {
            paste0("\n", p.value[as.character(facet.data)])
          }
        ),
        levels = paste0(
          facet,
          " = ",
          as.character(levels(facet.data)),
          if (length(p.value) > 1) {
            paste0("\n", p.value[as.character(levels(facet.data))])
          }
        )
      )
    }
  }

  if (!is.null(shape.name)) {
    shape <- sample.data[[shape.name]]
    if (is.character(shape) | is.factor(shape)) {
      if (!is.factor(shape)) {
        shape <- as.factor(shape)
      }
      plot.df$shape <- shape
    } else {
      shape <- NULL
    }
  } else {
    shape <- NULL
  }

  if (!is.null(size.name)) {
    size <- sample.data[[size.name]]
  } else {
    size <- NULL
  }

  hover.variables <- colnames(sample.data)[
    colnames(sample.data) %in% hover.variables
  ]

  hover.text.all <- rep(NA, nrow(sample.data))
  for (i in seq_len(nrow(sample.data))) {
    sample <- rownames(sample.data)[i]
    hover.text <- ""
    values <- sample.data[sample, , drop = FALSE] %>% as("data.frame")
    for (variable in hover.variables) {
      value <- values[[variable]]
      hover.text <- paste0(hover.text, variable, ": ", value, "<br>")
    }
    hover.text.all[i] <- hover.text
  }

  plot.df$hover.text <- hover.text.all

  if (".is_fit_sample" %in% colnames(sample.data)) {
    plot.df$.is_fit_sample <- sample.data$.is_fit_sample
  }

  point_aes <- aes(
    x = Comp1,
    y = Comp2,
    color = label, # OK if NULL
    text = hover.text,
    shape = shape, # OK if NULL
    size = size # OK if NULL
  )

  if (".is_fit_sample" %in% colnames(plot.df)) {
    if (!is.null(size.name)) {
      plt <- ggplot() +
        geom_point(
          point_aes,
          plot.df[!plot.df$.is_fit_sample, ],
          alpha = projected.alpha
        ) +
        geom_point(
          point_aes,
          plot.df[plot.df$.is_fit_sample, ],
          alpha = point.alpha
        ) +
        ggplot2::scale_size(range = c(point.size * 0.5, point.size * 3)) +
        labs(color = label.name, shape = shape.name, size = size.name)
    } else {
      plt <- ggplot() +
        geom_point(
          point_aes,
          plot.df[!plot.df$.is_fit_sample, ],
          alpha = projected.alpha,
          size = point.size
        ) +
        geom_point(
          point_aes,
          plot.df[plot.df$.is_fit_sample, ],
          alpha = point.alpha,
          size = point.size
        ) +
        labs(color = label.name, shape = shape.name, size = size.name)
    }
  } else {
    if (!is.null(size.name)) {
      plt <- ggplot() +
        geom_point(point_aes, plot.df, alpha = point.alpha) +
        ggplot2::scale_size(range = c(point.size * 0.5, point.size * 3)) +
        labs(color = label.name, shape = shape.name, size = size.name)
    } else {
      plt <- ggplot() +
        geom_point(point_aes, plot.df, alpha = point.alpha, size = point.size) +
        labs(color = label.name, shape = shape.name, size = size.name)
    }
  }

  if (ellipses & is.factor(label)) {
    if (fill.ellipses) {
      plt <- plt +
        stat_ellipse(
          aes(x = Comp1, y = Comp2, color = label, fill = label),
          plot.df,
          geom = "polygon",
          alpha = 0.2
        )
    } else {
      plt <- plt +
        stat_ellipse(aes(x = Comp1, y = Comp2, color = label), plot.df)
    }
  }

  plt <- plt +
    labs(
      x = paste0(
        dim.names[comp[1]],
        if (!is.null(prop.var.explained)) {
          paste0(" (", round(prop.var.explained[comp[1]]), "%)")
        } else {
          NULL
        }
      ),
      y = paste0(
        dim.names[comp[2]],
        if (!is.null(prop.var.explained)) {
          paste0(" (", round(prop.var.explained[comp[2]]), "%)")
        } else {
          NULL
        }
      ),
      # get rid of backticks in the model (formula) for the title
      title = paste0(
        "Beta-Diversity",
        if (is.null(model)) {
          ""
        } else {
          paste0(" ~ ", gsub("`", " ", model))
        },
        " (",
        beta.dispersion.fit$method,
        ")"
      ),

      subtitle = paste0(
        "taxa agglom: ",
        if (!is.null(taxrank)) {
          taxrank
        } else {
          "none"
        },
        if (!is.null(beta.dispersion.fit$fit_filter)) {
          paste0(
            "  fit subset: ",
            beta.dispersion.fit$fit_filter$name,
            " in {",
            paste(beta.dispersion.fit$fit_filter$values, collapse = ", "),
            "}"
          )
        },
        if (remove.na.from.plot) {
          "  NA's removed  "
        },
        if (!is.null(beta.dispersion.fit$dist)) {
          paste0(" dist=", beta.dispersion.fit$dist)
        },
        if (
          is.null(beta.dispersion.fit$method) ||
            !tolower(beta.dispersion.fit$method) %in% c("pcoa", "tsne", "umap")
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
        if (biplot.loadings && !is.null(loadings)) {
          paste0("  taxa=", arrow.cutoff.load)
        } else {
          NULL
        },
        if (biplot.covariates && !is.null(covariates)) {
          paste0("  predictors=", arrow.cutoff.covar)
        } else {
          NULL
        },
        # add stats to the title only if there no or only one facet
        if (length(p.value) == 1) {
          paste0("\n", p.value)
        } else {
          NULL
        }
      )
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5)
    ) +
    ggsci::scale_fill_npg()

  if (biplot & (!is.null(loadings) | !is.null(covariates))) {
    scale.factor.load <- NULL
    scale.factor.covar <- NULL

    arrow.df <- data.frame()

    # Plot species arrows = loadings here
    if (!is.null(loadings) & biplot.loadings) {
      # Account for absent loadings
      filtered.loadings <- na.omit(loadings[, comp, drop = FALSE]) # TODO: how come that some loadings are NaN?
      #       is it for all models?
      # Filter arrows by length (get rid of smaller arrows)
      arrow.lengths <- apply(filtered.loadings, 1, function(x) sqrt(sum(x^2)))
      filtered.loadings <- filtered.loadings[
        arrow.lengths / max(arrow.lengths) > arrow.cutoff.load,
        ,
        drop = FALSE
      ]

      # If there are loadings left...
      if (nrow(filtered.loadings) > 0) {
        arrow.load <- data.frame(
          Comp1 = filtered.loadings[, 1],
          Comp2 = filtered.loadings[, 2],
          type = "load"
        )

        rownames(arrow.load) <- row.names(filtered.loadings)

        # Scale the secondary axes for loadings so that loadings and scores are on
        # comparable scale
        # NOTE: both components are scaled by the same factor, so arrows point
        # in the same direction as before, only their size changes!
        lim.load <- max(abs(arrow.load$Comp1), abs(arrow.load$Comp2))
        scale.factor.load <- max(abs(coords)) / lim.load
        arrow.load$Comp1 <- arrow.load$Comp1 * scale.factor.load
        arrow.load$Comp2 <- arrow.load$Comp2 * scale.factor.load

        # If phyloseq object provided, get hover info from it - taxonomy for each species (arrow)
        if (!is.null(tax.table)) {
          hover.text.all <- rep(NA, nrow(arrow.load))
          for (i in 1:nrow(arrow.load)) {
            # it makes sure that taxa names match
            taxon <- rownames(arrow.load)[i]
            hover.text <- ""
            values <- tax.table[taxon, , drop = FALSE]
            for (variable in colnames(values)) {
              value <- values[[variable]]
              hover.text <- paste0(hover.text, variable, ": ", value, "<br>")
            }
            hover.text.all[i] <- hover.text
          }

          arrow.load$hover.text <- hover.text.all
        } else {
          arrow.load$hover.text <- rownames(arrow.load)
        }

        if (!is.null(arrow.taxonomy.labels)) {
          # use taxonomy as label instead of taxa names
          arrow.load$Names <- apply(
            tax.table[rownames(arrow.load), ],
            1,
            function(row) {
              gsub("NA", " ", paste(row[arrow.taxonomy.labels], collapse = ":"))
            }
          )
        } else {
          arrow.load$Names <- rownames(arrow.load)
        }
        arrow.df <- rbind(arrow.df, arrow.load)
      }
    }
    # Similar for covariates (predictors in case of constrained model)
    if (!is.null(covariates) & biplot.covariates) {
      # Filter out absent data or short arrows
      filtered.covariates <- na.omit(covariates[, comp, drop = FALSE])
      arrow.lengths <- apply(filtered.covariates, 1, function(x) sqrt(sum(x^2)))
      filtered.covariates <- filtered.covariates[
        arrow.lengths / max(arrow.lengths) > arrow.cutoff.covar,
        ,
        drop = FALSE
      ]
      # If some are left...
      if (nrow(filtered.covariates) > 0) {
        arrow.covar <- data.frame(
          Comp1 = filtered.covariates[, 1],
          Comp2 = filtered.covariates[, 2],
          type = "covar"
        )
        rownames(arrow.covar) <- row.names(filtered.covariates)

        # Scale in the same manner as for loading but notice that they are scaled
        # by a different factor, so that all are of comparable size on plot
        lim.covar <- max(abs(arrow.covar$Comp1), abs(arrow.covar$Comp2))
        scale.factor.covar <- max(abs(coords)) / lim.covar
        arrow.covar$Comp1 <- arrow.covar$Comp1 * scale.factor.covar
        arrow.covar$Comp2 <- arrow.covar$Comp2 * scale.factor.covar
        # Hover text is just the name of covariate
        arrow.covar$hover.text <- rownames(arrow.covar)
        arrow.covar$Names <- gsub("`", "", rownames(arrow.covar))
        arrow.df <- rbind(arrow.df, arrow.covar)
      }
    }

    # If there are loadings and/or covariates, plot arrows
    if (nrow(arrow.df) > 0) {
      arrow.params <- arrow(
        type = "closed",
        angle = 20,
        length = unit(0.1, "inches")
      )

      if ("load" %in% arrow.df$type) {
        if (color.arrows.by.taxa) {
          plt <- plt +
            geom_segment(
              data = arrow.df[arrow.df$type == "load", ],
              aes(
                x = 0,
                y = 0,
                xend = Comp1,
                yend = Comp2,
                text = hover.text,
                color = Names
              ), # Arrow from (0,0) to each point
              arrow = arrow.params,
              size = 0.7, # Size of the arrows
              alpha = 0.7 # Transparency of the arrows
            )
        } else {
          plt <- plt +
            geom_segment(
              data = arrow.df[arrow.df$type == "load", ],
              aes(x = 0, y = 0, xend = Comp1, yend = Comp2, text = hover.text), # Arrow from (0,0) to each point
              arrow = arrow.params,
              color = "darkgrey",
              size = 0.7, # Size of the arrows
              alpha = 0.7 # Transparency of the arrows
            )
        }
      }

      if ("covar" %in% arrow.df$type) {
        plt <- plt +
          geom_segment(
            data = arrow.df[arrow.df$type == "covar", ],
            aes(x = 0, y = 0, xend = Comp1, yend = Comp2, text = hover.text), # Arrow from (0,0) to each point
            arrow = arrow.params,
            color = "darkred",
            size = 0.7, # Size of the arrows
            alpha = 0.7 # Transparency of the arrows
          )
      }
      if (arrow.labels) {
        if (repel) {
          plt <- plt +
            ggrepel::geom_text_repel(
              data = arrow.df,
              aes(x = Comp1, y = Comp2, label = Names),
              max.overlaps = max.overlaps,
              show.legend = FALSE
            )
        } else {
          plt <- plt +
            geom_text(
              data = arrow.df,
              aes(x = Comp1, y = Comp2, label = Names, text = hover.text)
            )
        }
      }

      # if there are marginal plots or covariates and loadings at the same time,
      # don't show secondary axis names to avoid visual mess
      # TODO: or show them finally?
      if (
        is.null(marginal.plot) &
          xor(is.null(scale.factor.load), is.null(scale.factor.covar)) # <=> only loadings or covariates, not both
      ) {
        scale.factor <- c(scale.factor.load, scale.factor.covar) # = the one which is not NULL
        plt <- plt +
          scale_x_continuous(
            sec.axis = sec_axis(
              trans = ~ . / scale.factor,
              name = paste0(
                ifelse(raw.loadings, 'Loadings ', 'Correlations '),
                dim.names[comp[1]]
              )
            )
          ) +
          scale_y_continuous(
            sec.axis = sec_axis(
              trans = ~ . / scale.factor,
              name = paste0(
                ifelse(raw.loadings, 'Loadings ', 'Correlations '),
                dim.names[comp[2]]
              )
            )
          )
      }
    }
  }
  if (grid.mode) {
    if (!is.null(facet.row) && !is.null(facet.col)) {
      plt <- plt + facet_grid(facet.row ~ facet.col, scales = "fixed")
      if (!is.null(p.value.df) && nrow(p.value.df) > 0) {
        annot.df <- data.frame(
          facet.row = factor(
            paste0(facet.row, " = ", p.value.df$facet.row),
            levels = levels(plot.df$facet.row)
          ),
          facet.col = factor(
            paste0(facet.col, " = ", p.value.df$facet.col),
            levels = levels(plot.df$facet.col)
          ),
          p.label = p.value.df$p.label,
          x = Inf,
          y = Inf
        )
        plt <- plt +
          geom_text(
            data = annot.df,
            aes(x = x, y = y, label = p.label),
            hjust = 1.05,
            vjust = 1.5,
            size = 3,
            inherit.aes = FALSE
          )
      }
    } else if (!is.null(facet.row)) {
      plt <- plt + facet_grid(rows = vars(facet.row), scales = "fixed")
    } else if (!is.null(facet.col)) {
      plt <- plt + facet_grid(cols = vars(facet.col), scales = "fixed")
    }
  } else if (!is.null(facet) && "facet" %in% colnames(plot.df)) {
    plt <- plt +
      facet_wrap(
        . ~ facet,
        scales = "fixed",
        ncol = smart_facet_ncol(nlevels(plot.df$facet))
      )
  }

  # Marginal plot in case of factor variable - boxplot, density...
  if (
    !is.null(marginal.plot) & is.factor(label) & is.null(facet) & !grid.mode
  ) {
    plt <- plt +
      theme(
        legend.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "left"
      ) # place legend to the left so that it doesn't interfere with the marginal plot
    plt <- ggExtra::ggMarginal(
      plt,
      type = marginal.plot,
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
# TODO: this function is better since it doesn't assume any facets for animation,
# implement it in feature selector
#' Animate a Ggplot Object by a Variable
#' @export
animate_by_variable <- function(
  ggplot.obj, # ggplot2 object
  animation.variable.name,
  return.anim = TRUE,
  save.path = NULL # gif file name
) {
  plt <- ggplot.obj

  plt <-
    plt +
    gganimate::transition_states(
      get(animation.variable.name),
      transition_length = 2,
      state_length = 3
    ) +
    labs(subtitle = paste0(animation.variable.name, " = {closest_state}")) +
    ease_aes('linear')
  anim <-
    animate(
      plt,
      nframes = 150, # More frames for better smoothness
      fps = 25, # Higher FPS for fluid animation
      width = 900,
      height = 700,
      res = 200, # Higher resolution
      renderer = gifski_renderer(loop = TRUE) # High-quality GIF output
    )

  if (!is.null(save.path)) {
    anim_save(save.path, anim)
  }

  if (return.anim) {
    return(anim)
  } else {
    return(NULL)
  }
}
