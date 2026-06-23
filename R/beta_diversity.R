#' Get Beta-Diversity from Phyloseq Object
#' @export
get_beta_dispersion = function(
  physeq,
  taxrank = NULL, # taxrank to agglomerate by (e.g."Phylum")
  fraction_id_name = NULL,
  fraction_ids = NULL,
  transform.abundances = "identity",
  dist = "bray", # for PCoA, NMDS or dbRDA
  method = "PCoA",
  ndims = 10, # number of ordination axes to retain in coords/loadings/
  # covariates. Models can have up to nsamples - 1 axes, but only a couple
  # are ever plotted - keeping them all makes correlating every taxon
  # against every axis for `loadings` (and storing the resulting matrices)
  # far slower than necessary. Raise this if you need to plot axes beyond
  # the default range.
  model = NULL, # in case of constrained model,
  # string with a part of the formula
  # for covariates e.g. "Var1+Var2"
  confounders = c(),
  species = FALSE, # NOTE: this argument is not used, kept from the original function
  # In case of tSNE method
  pca = TRUE, # whether initial pca step should be performed
  perplexity = NULL,
  # In case of UMAP method
  nb.neighbors = 15,
  min.dist = 0.1
) {
  method.orig = method # keep method name with the original case (used later for plot title)
  method = tolower(method)
  # Lists by scaling: [[1]] - vegan's scaling 1 -> interpretable distances between samples,
  #                   [[2]] - interpretable angles between arrows (as correlations)
  # Initiate species coordinates (arrows) to NULL still sometimes they're not availible
  loadings = list() # = species scores
  # NOTE: term loadings is not really appropriate for all models,
  # but it's used here for convenience

  covariates = list() # = covariate scores for biplot in case of constrained models

  coords = list()

  eigen.values = NULL # eventually for scree plots and % of explained variance

  if (class(physeq) != "phyloseq") {
    stop("Need a phyloseq or PhyloIgSeq object")
  }

  if (!is.null(fraction_id_name) & !is.null(fraction_ids)) {
    physeq = prune_samples(
      sample_data(physeq)[[fraction_id_name]] %in% fraction_ids,
      physeq
    )
  }

  # Force taxa to be columns of otu table - like in internal phyloseq function phyloseq:::veganifyOTU()
  if (taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }

  # NOTE: agglomerate taxa BEFORE transforming the data

  # Agglomerate taxa by taxrank
  if (!is.null(taxrank)) {
    physeq = tax_glom(physeq = physeq, taxrank = taxrank)
    taxa_names(physeq) = make.unique(tax_table(physeq)[, taxrank]) # NOTE: names of taxons are not necesserily unique!
  }

  if (!is.null(transform.abundances) && transform.abundances != "identity") {
    physeq = microbiome::transform(
      physeq,
      transform = transform.abundances,
      target = "OTU", # TODO: and still, clr will scale over samples
      shift = 0, # pseudocount added (shifts baseline)
      scale = 1, # if transform is "scale"
      log10 = TRUE,
      reference = 1
    )
  }

  # Unconstrained ordination models

  if (method %in% c("pcoa", "tsne", "umap")) {
    # Compute distance matrix myself
    # NOTE: tSNE and UMAP can take as well a precomputed distance metric
    if (class(dist) == "dist") {
      # if distance matrix furnished directly (not recommended)
      dist.matrix = dist
    } else {
      if (is.null(access(physeq, "phy_tree")) & dist == "unifrac") {
        # FIXME: it was Unifrac before - check the spelling
        print(
          "No phylogenetic tree in this phyloseq object, bray-curtis distance selected."
        )
        dist.matrix = phyloseq::distance(physeq, method = "bray")
      } else {
        # NOTE: euclidean and manhattan use z-score standardization (decostand "standardize")
        # so that high-abundance taxa don't dominate the distance in ordination.
        # This produces STANDARDIZED euclidean/manhattan, intentionally different from
        # sparse_distance(ps, "euclidean"/"manhattan") which returns raw distances.
        # The zero-variance filter (apply var > 0) prevents division-by-zero during
        # z-score scaling; if sparse data leaves very few taxa, results may be degenerate.
        # NOTE: if transform.abundances was applied above (e.g. CLR, Hellinger), this
        # z-score step stacks on top of it. For CLR that is still useful (removes
        # scale differences between taxa); for Hellinger it is redundant but harmless.
        # TODO: consider routing through sparse_distance once standardization can be
        # done sparsely, to avoid materializing the full dense OTU matrix here.
        if (!is.null(dist) && dist %in% c("euclidean", "manhattan")) {
          otu.mat = as(otu_table(physeq), "matrix")
          dist.matrix = vegdist(
            decostand(
              otu.mat[, apply(otu.mat, 2, var, na.rm = TRUE) > 0],
              method = "standardize"
            ),
            method = dist
          )
        } else {
          dist.matrix = sparse_distance(physeq, method = dist)
        }
      }
    }
    dist.matrix = as.matrix(dist.matrix)
    if (method == "pcoa") {
      fit = vegan::wcmdscale(dist.matrix, eig = TRUE)
      eigen.values = vegan::eigenvals(fit)
      n.axes = min(ndims, length(eigen.values))

      coords[[1]] = scores(
        fit,
        display = "sites",
        choices = 1:n.axes,
        scaling = 1
      )
      coords[[2]] = scores(
        fit,
        display = "sites",
        choices = 1:n.axes,
        scaling = 2
      )

      # NOTE: we don't get loadings from the model but we can nevertheless
      #       correlate abundances against axes - be careful when interpreting.
      # We used to get these via envfit(), which fits a multiple regression of
      # each taxon onto ALL requested axes. PCoA axes are orthogonal though, so
      # that multiple regression reduces to a plain per-axis correlation -
      # envfit's QR-based implementation doesn't exploit this and becomes very
      # slow as the number of axes/taxa grows, so just correlate directly.
      loadings[[1]] = suppressWarnings(cor(
        as(otu_table(physeq), "matrix"),
        coords[[1]]
      ))
      # TODO:for now same thing for both "scalings"
      # loadings[[2]] = wascores(fit$points, otu_table(physeq)) ≈ scaling 2 ?
    } else if (method == "tsne") {
      max.perplexity = floor((nsamples(physeq) - 1) / 3)

      if (is.null(perplexity)) {
        perplexity = 30
      }

      if (perplexity > max.perplexity) {
        warning(
          "Perplexity must not exceed 3 * perplexity < nrow(X) - 1, set to this value"
        )
        perplexity = max.perplexity
      }

      fit = Rtsne(
        dist.matrix,
        pca = pca,
        perplexity = perplexity,
        is_distance = TRUE
      )
      coords[[1]] = fit$Y
      coords[[2]] = fit$Y
    } else if (method == "umap") {
      umap.config = umap::umap.defaults
      if (!is.null(nb.neighbors)) {
        umap.config$n_neighbors = nb.neighbors
      }

      if (!is.null(min.dist)) {
        umap.config$min_dist = min.dist
      }

      umap.config$input = "dist" # distance matrix as input

      fit = umap::umap(dist.matrix, config = umap.config)

      coords[[1]] = fit$layout
      coords[[2]] = fit$layout
    }
  } else if (method %in% c("nmds", "pca", "ca", "dca")) {
    otu = as(otu_table(reverseASV(physeq)), 'matrix')

    if (method == "nmds") {
      fit = vegan::metaMDS(otu, distance = dist)
    } else if (method == "pca") {
      # RDA without constraint = PCA
      fit = vegan::rda(otu, scale = TRUE)
    } else if (method == "ca") {
      fit = vegan::cca(otu)
    } else if (method == "dca") {
      fit = vegan::decorana(otu)
    }
    eigen.values = vegan::eigenvals(fit)
    n.axes = min(ndims, length(eigen.values))

    coords[[1]] = scores(
      fit,
      display = "sites",
      choices = 1:n.axes,
      scaling = 1
    )
    coords[[2]] = scores(
      fit,
      display = "sites",
      choices = 1:n.axes,
      scaling = 2
    )
    loadings[[1]] = scores(
      fit,
      display = "species",
      choices = 1:n.axes,
      scaling = 1
    )
    loadings[[2]] = scores(
      fit,
      display = "species",
      choices = 1:n.axes,
      scaling = 2
    )
  } else if (!is.null(model)) {
    # Constrained ordination models

    if (!is.null(confounders)) {
      model = paste0(
        model,
        " + ",
        "Condition(",
        paste(paste0("`", confounders, "`"), collapse = " + "),
        ")"
      )
    }

    formula = as.formula(paste(
      "as(otu_table(reverseASV(physeq)), 'matrix')",
      "~",
      model
    ))
    df = as(sample_data(physeq), "data.frame")

    if (method == "cca") {
      fit = vegan::cca(formula, data = df, na.action = na.exclude)
    } else if (method == "rda") {
      fit = vegan::rda(formula, data = df, na.action = na.exclude)
    } else if (method == "dbrda") {
      fit = vegan::capscale(
        formula,
        data = df,
        na.action = na.exclude,
        dist = dist
      )
    } else {
      stop("Invalid method")
    }

    eigen.values = vegan::eigenvals(fit)
    n.axes = min(ndims, length(eigen.values))

    coords[[1]] = scores(
      fit,
      display = "sites",
      choices = 1:n.axes,
      scaling = 1
    )
    coords[[2]] = scores(
      fit,
      display = "sites",
      choices = 1:n.axes,
      scaling = 2
    )
    loadings[[1]] = scores(
      fit,
      display = "species",
      choices = 1:n.axes,
      scaling = 1
    )
    loadings[[2]] = scores(
      fit,
      display = "species",
      choices = 1:n.axes,
      scaling = 2
    )
    covariates[[1]] = scores(
      fit,
      display = "bp",
      choices = 1:n.axes,
      scaling = 1
    )
    covariates[[2]] = scores(
      fit,
      display = "bp",
      choices = 1:n.axes,
      scaling = 2
    )
  } else {
    if (method %in% c("cca", "rda", "dbrda")) {
      stop("Model is required for constrained models")
    } else {
      stop("Invalid method")
    }
  }

  if (length(eigen.values) == 1 && is.na(eigen.values)) {
    eigen.values = NULL
  }

  return(list(
    fit = fit,
    taxrank = taxrank,

    # lists: [[1]] - scaling 1, [[2]] - sclaing 2 (see vegan's scaling)
    coords = coords, # = sample (site) scores
    loadings = loadings, # = variable (species) scores
    covariates = covariates, # = covariate scores (if constrained model)

    eigen.values = eigen.values, # NULL if no eigen values provided for this method

    dist = dist, # distance metric
    method = method.orig, # name of ordination method
    model = if (method %in% c("rda", "cca", "dbrda")) {
      model
    } else {
      NULL
    }, # formula with covariates
    sample.data = sample_data(physeq) %>% as("data.frame"),
    tax.table = tax_table(physeq) %>% as.data.frame()
  ))
}

#' Permutation Test for Beta-Diversity Ordination
#' @export
stat_beta_dispersion = function(
  beta.dispersion.fit,
  facet.name = NULL, # single facet (backward compat alias for facet.row.name)
  facet.row.name = NULL, # row facet variable name
  facet.col.name = NULL, # col facet variable name; when set alongside facet.row.name, enables grid mode
  comp = c(1, 2), # if NULL, consider all components
  strata.name = NULL, # in case of paired design
  label.name
) {
  if (!is.null(facet.name) && is.null(facet.row.name)) {
    facet.row.name = facet.name
  }

  grid.mode = !is.null(facet.row.name) && !is.null(facet.col.name)

  if (!label.name %in% colnames(beta.dispersion.fit$sample.data)) {
    warning("Wrong label name")
    return(NULL)
  }

  if (!is.null(comp)) {
    stat.data.full = data.frame(
      Comp1 = beta.dispersion.fit$coords[[1]][, comp[1]],
      Comp2 = beta.dispersion.fit$coords[[1]][, comp[2]],
      Label = beta.dispersion.fit$sample.data[[label.name]]
    )
  } else {
    stat.data.full = data.frame(
      beta.dispersion.fit$coords[[1]],
      Label = beta.dispersion.fit$sample.data[[label.name]]
    )
  }

  if (!is.null(strata.name)) {
    stat.data.full[[strata.name]] = beta.dispersion.fit$sample.data[[
      strata.name
    ]]
  }

  # Embed facet variables into stat.data.full before na.omit so subsetting is consistent
  if (
    !is.null(facet.row.name) &&
      facet.row.name %in% colnames(beta.dispersion.fit$sample.data)
  ) {
    stat.data.full[[".facet.row"]] = beta.dispersion.fit$sample.data[[
      facet.row.name
    ]]
  }
  if (
    grid.mode && facet.col.name %in% colnames(beta.dispersion.fit$sample.data)
  ) {
    stat.data.full[[".facet.col"]] = beta.dispersion.fit$sample.data[[
      facet.col.name
    ]]
  }

  stat.data.full = na.omit(stat.data.full)

  # Shared helper: run PERMANOVA or envfit on a single data slice
  run_test = function(stat.data) {
    stat.data[[".facet.row"]] = NULL
    stat.data[[".facet.col"]] = NULL
    if (is.factor(stat.data$Label) || is.character(stat.data$Label)) {
      stat.data = na.omit(stat.data)
      if (!is.null(strata.name)) {
        strata = stat.data[[strata.name]]
        stat.data[[strata.name]] = NULL
      } else {
        strata = NULL
      }
      test.result = vegan::adonis2(
        stat.data[, colnames(stat.data) != "Label"] ~ Label,
        data = stat.data,
        strata = strata,
        permutations = 999,
        method = "euclidean"
      )
      p.raw = test.result$`Pr(>F)`[1]
      p.text = if (!is.null(strata)) {
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
      test.result = vegan::envfit(
        stat.data[, colnames(stat.data) != "Label"] ~ Label,
        data = stat.data,
        permutations = 999
      )
      p.raw = test.result$vectors$pvals[1]
      p.text = paste0(
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

  dim_used = if (!is.null(comp)) {
    comp
  } else {
    paste("all", ncol(beta.dispersion.fit$coords))
  }

  if (grid.mode) {
    row.vals = unique(stat.data.full[[".facet.row"]])
    col.vals = unique(stat.data.full[[".facet.col"]])
    combinations = expand.grid(
      row = row.vals,
      col = col.vals,
      stringsAsFactors = FALSE
    )

    test.result.all = list()
    p.value.df = data.frame(
      facet.row = character(0),
      facet.col = character(0),
      p.label = character(0),
      p.raw = numeric(0),
      stringsAsFactors = FALSE
    )
    stat.type = NULL

    for (i in seq_len(nrow(combinations))) {
      row.val = combinations$row[i]
      col.val = combinations$col[i]
      stat.data = stat.data.full[
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
      res = run_test(stat.data)
      key = paste(row.val, col.val, sep = "//")
      test.result.all[[key]] = res$test.result
      stat.type = res$stat
      p.value.df = rbind(
        p.value.df,
        data.frame(
          facet.row = as.character(row.val),
          facet.col = as.character(col.val),
          p.label = res$p.text,
          p.raw = res$p.raw,
          stringsAsFactors = FALSE
        )
      )
    }

    return(list(
      stat = stat.type,
      label.name = label.name,
      facet.name = facet.row.name,
      facet.row.name = facet.row.name,
      facet.col.name = facet.col.name,
      dim_used = dim_used,
      test.result = if (length(test.result.all) > 0) test.result.all else NULL,
      p.value = NULL,
      p.value.df = p.value.df,
      p.value.raw = p.value.df$p.raw
    ))
  }

  # Single-facet or no-facet mode
  facets = if (".facet.row" %in% colnames(stat.data.full)) {
    unique(stat.data.full[[".facet.row"]])
  } else {
    NA
  }

  test.result.all = list()
  p.value.all = c()
  p.value.raw = c()
  stat.type = NULL

  for (facet in facets) {
    if (is.na(facet)) {
      stat.data = stat.data.full
      idx = 1
    } else {
      stat.data = stat.data.full[stat.data.full[[".facet.row"]] == facet, ]
      idx = facet
    }

    if (length(unique(stat.data$Label)) == 1) {
      warning("Variable doesn't vary, no stats performed...")
      return(NULL)
    }

    res = run_test(stat.data)
    test.result.all[[idx]] = res$test.result
    stat.type = res$stat
    p.value.raw[idx] = res$p.raw
    p.value.all[idx] = res$p.text
  }

  return(list(
    stat = stat.type,
    label.name = label.name,
    facet.name = facet.row.name,
    facet.row.name = facet.row.name,
    facet.col.name = NULL,
    dim_used = dim_used,
    test.result = if (length(test.result.all) > 0) test.result.all else NULL,
    p.value = p.value.all,
    p.value.df = NULL,
    p.value.raw = p.value.raw
  ))
}

#' Plot Beta-Diversity
#' @export
plot_beta_dispersion = function(
  fraction_id_name = NULL,
  fraction_ids = NULL,
  scaling = 1, # see vegan's scalings 1 and 2
  beta.dispersion.fit, # a list of p, coords, species.coords.x, species.coords.y
  stat.beta.dispersion, # a list of stat, group, res, p.val
  axis_x = 1,
  axis_y = 2,
  nf = 5,
  type = "boxplot",
  color_vector = c("cyan4", "brown", "deepskyblue", "black", "red"),
  legend_title = NULL,
  lwd = 1,
  conf = 0.9,
  cex = 2,
  font = 2,
  pch = 20,
  draw = "lines",
  ylimits = "auto",
  xlimits = "auto",
  text = FALSE,
  ncol = 1,
  species = FALSE,
  x.intersp = 1,
  y.intersp = 0.5,
  where = "topleft",
  inset = 0.2,
  pca = TRUE,
  stat.cex = 2,
  legend.cex = 2,
  widths = c(1, 1),
  heights = c(1, 1),
  margins = c(1, 1, 1, 1),
  ...
) {
  # /!\ This function is used only to imitate the behavior of the original function
  # by assembling it with other function in full_beta_dispersion()
  # It will not be used in the PhyloIgSeq app

  fit = beta.dispersion.fit$fit

  if (length(beta.dispersion.fit$coords) != 0) {
    coords = beta.dispersion.fit$coords[[scaling]]
  } else {
    coords = NULL
  }

  if (length(beta.dispersion.fit$loadings) != 0) {
    loadings = beta.dispersion.fit$loadings[[scaling]]
  } else {
    loadings = NULL
  }

  species.coords.x = loadings[, 1]
  species.coords.y = loadings[, 2]
  dist = beta.dispersion.fit$dist
  method = beta.dispersion.fit$method

  label = as.factor(beta.dispersion.fit$sample.data[[group]])

  if (!is.null(stat.beta.dispersion$label.name)) {
    stat = stat.beta.dispersion$stat
    group = stat.beta.dispersion$label.name
  } else {
    stat = NULL
    group = NULL
  }

  if (!is.null(stat.beta.dispersion$test.result)) {
    test.result = stat.beta.dispersion$test.result[[1]]
  } else {
    test.result = NULL
  }

  p.value = stat.beta.dispersion$p.value

  if (is.null(group)) {
    # TODO: account for when the group=NULL
    stop("Need factor to segregate result.")
  }

  # Set up the plot itself

  if (xlimits == "auto") {
    xlimits = c(coords[, axis_x] %>% max / 0.9, coords[, 1] %>% min / 0.9)
  }

  if (ylimits == "auto") {
    ylimits = c(coords[, axis_y] %>% max / 0.9, coords[, 2] %>% min / 0.9)
  } else {
    ylimits = ylim
    xlimits = ylim
  }

  col1 = color_vector[unique(label)]
  col2 = color_vector[label]

  if (type == "boxplot") {
    # prepare the layout matrix to display boxplots on the margins
    layout(
      matrix(
        c(2, 2, 2, 4, 1, 1, 1, 3, 1, 1, 1, 3, 1, 1, 1, 3),
        nrow = 4,
        ncol = 4,
        byrow = TRUE
      ),
      widths = widths,
      heights = heights
    )

    par(mar = margins)
    on.exit(layout(matrix(c(1, 1)))) # reset the layout to one plotting region
    # Plot coordinates (scores)
    plot(
      coords[, axis_x],
      coords[, axis_y],
      bg = col2,
      axes = FALSE,
      xlab = "",
      ylab = "",
      las = 2,
      pch = 21,
      cex = cex,
      ylim = ylimits,
      xlim = xlimits,
      ...
    )
    # display groups
    disp.groups = ordispider(
      coords,
      groups = label,
      col = adjustcolor(color_vector, alpha = 0.3),
      lwd = lwd,
      ylim = ylimits,
      xlim = xlimits
    )
    # add ellipse
    # TODO: make limits of the plot englobe all the ellipse
    ordiellipse(
      coords,
      groups = label,
      conf = conf,
      col = color_vector,
      lwd = lwd,
      draw = draw,
      ylim = ylimits,
      xlim = xlimits
    )
    # TODO: this line is redundant?
    points(
      coords[, axis_x],
      coords[, axis_y],
      bg = col2,
      xlab = "",
      ylab = "",
      las = 2,
      pch = 21,
      cex = cex,
      ylim = ylimits,
      xlim = xlimits,
      ...
    )

    if (text) {
      # TODO: might be a problem with labels and levels
      text(
        x = unique(disp.groups[, 1:2]),
        labels = unique(group),
        col = "black",
        cex = cex,
        font = font
      )
    }

    if (!is.null(stat) && (stat == "permanova" | stat == "envfit")) {
      # plot statistics
      legend(where, inset = 0.2, legend = p.value, bty = "n", cex = stat.cex)
    }

    # Plot marginal boxplots
    boxplot(
      coords[, axis_x] ~ label,
      data = coords,
      horizontal = TRUE,
      axes = FALSE,
      xlab = NULL,
      ylab = NULL,
      col = adjustcolor(color_vector, alpha = 0.7),
      xaxt = "n",
      lwd = lwd / 2,
      ylim = xlimits
    )

    stripchart(
      coords[, axis_x] ~ label,
      data = coords,
      method = "jitter",
      vertical = FALSE,
      add = TRUE,
      pch = 21,
      col = "black",
      bg = "gray",
      lwd = lwd / 2,
      ylim = xlimits
    )

    boxplot(
      coords[, axis_y] ~ label,
      data = coords,
      axes = FALSE,
      ylab = "",
      xlab = "",
      col = adjustcolor(color_vector, alpha = 0.7),
      lwd = lwd / 2,
      ylim = ylimits
    )
    stripchart(
      coords[, axis_y] ~ label,
      data = coords,
      method = "jitter",
      vertical = TRUE,
      add = TRUE,
      pch = 21,
      col = "black",
      bg = "gray",
      lwd = lwd / 2,
      ylim = ylimits
    )

    # Plot legend
    # TODO: what is it for?
    plot(
      NULL,
      xlim = c(0, 1),
      ylim = c(0, 1),
      ylab = "",
      xlab = "",
      axes = FALSE
    )
    legend(
      "center",
      legend = unique(label),
      col = col1,
      title = legend_title,
      pch = 20,
      cex = legend.cex,
      bty = "n",
      ncol = ncol,
      y.intersp = y.intersp
    )

    # return(head(p$eig/sum(p$eig)*100,5))
  }

  if (type == "pure") {
    par(mar = c(0.5, 0.5, 1, 0.5))
    plot(
      coords[, axis_x],
      coords[, axis_y],
      bg = col2,
      axes = FALSE,
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "",
      las = 2,
      pch = 21,
      cex = cex,
      ylim = ylimits,
      xlim = xlimits,
      ...
    )

    disp.groups = ordispider(
      coords,
      groups = label,
      col = adjustcolor(color_vector, alpha = 0.3),
      lwd = lwd,
      ylim = ylimits,
      xlim = xlimits,
    )
    ordiellipse(
      coords,
      groups = label,
      conf = conf,
      col = adjustcolor(color_vector, alpha = 0.3),
      lwd = lwd,
      draw = draw,
      ylim = ylimits,
      xlim = xlimits
    )
    points(
      coords[, axis_x],
      coords[, axis_y],
      bg = col2,
      xlab = "",
      ylab = "",
      las = 2,
      pch = 21,
      cex = cex,
      ylim = ylimits,
      xlim = xlimits,
      ...
    )

    legend(
      where,
      legend = unique(label),
      col = col1,
      title = legend_title,
      pch = 20,
      cex = legend.cex,
      bty = "n",
      ncol = ncol,
      y.intersp = y.intersp
    )

    if (stat == "permanova" | stat == "envfit") {
      legend(where, inset = 0.2, legend = p.value, bty = "n", cex = stat.cex)
    }
  }

  if (type == "arrows") {
    # TODO: and for others it's okay?
    if (method == "PCoA") {
      stop("Can't plot loadings for a PCoA, use NMDS for that purpose")
    }

    plot(fit, type = "n", axes = FALSE, xlab = "", ylab = "", bty = "n")
    abline(h = 0, v = 0, col = "white", lwd = 3)

    disp.groups = ordiellipse(
      fit,
      groups = label,
      conf = conf,
      col = adjustcolor(color_vector, alpha = 0.3),
      lwd = lwd,
      draw = draw,
      ylim = ylimits,
      xlim = xlimits
    )
    disp.groups = ordispider(
      fit,
      groups = label,
      col = adjustcolor(color_vector, alpha = 0.3),
      lwd = lwd,
      xlim = c(min(species.coords.x), max(species.coords.x)),
      ylim = c(min(species.coords.y), max(species.coords.y))
    )

    points(
      fit,
      display = "sites",
      bg = col2,
      xlab = "",
      ylab = "",
      las = 2,
      pch = 21,
      cex = cex
    )

    arrows(
      x0 = 0,
      x1 = species.coords.x,
      y0 = 0,
      y1 = species.coords.y,
      lwd = lwd / 1.5
    )

    legend(
      where,
      legend = unique(label),
      col = col1,
      title = legend_title,
      pch = 20,
      cex = legend.cex,
      bty = "n",
      ncol = ncol,
      y.intersp = y.intersp
    )

    if (text) {
      text(
        x = unique(disp.groups[, 1:2]),
        labels = unique(label),
        col = "black",
        cex = cex,
        font = font
      )
    }
    if (stat == "permanova" | stat == "envfit") {
      legend(where, inset = 0.2, legend = p.value, bty = "n", cex = stat.cex)
    }
  }

  return(NULL) # builds plot as a side-effect
}


#' Get and Plot Beta-Diversity from Phyloseq Object
#' @export
full_beta_dispersion = function(
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
  beta.dispersion.fit = get_beta_dispersion(
    physeq = physeq,
    taxrank = taxrank,
    fraction_id_name = fraction_id_name,
    fraction_ids = fraction_ids,
    dist = dist,
    method = method,
    model = model,
    species = species
  )
  if (!is.null(stat) && stat != FALSE && stat != "none") {
    stat.beta.dispersion = stat_beta_dispersion(
      beta.dispersion.fit,
      comp = c(axis_x, axis_y),
      label.name = group
    )
  } else {
    stat.beta.dispersion = NULL
  }

  plot_beta_dispersion(
    fraction_id_name = fraction_id_name,
    fraction_ids = fraction_ids,
    beta.dispersion.fit = beta.dispersion.fit,
    stat.beta.dispersion = stat.beta.dispersion,
    ...
  )

  return(list(
    reduction = beta.dispersion.fit$fit,
    stat = stat.beta.dispersion[[1]]$test.result
  ))
}

#' Scree Plot from Eigenvalues
#' @export
scree_plot = function(eigen.values, max.nb.comp = 10) {
  if (is.null(eigen.values)) {
    warning("No eigen values are furnished")
    return(NULL)
  }
  total_var = sum(eigen.values[eigen.values > 0])
  eigen.values = eigen.values[1:min(max.nb.comp, length(eigen.values))]
  plot.data = data.frame(
    prop_var = eigen.values / total_var * 100,
    dim = if (!is.null(names(eigen.values))) {
      names(eigen.values)
    } else {
      1:length(eigen.values)
    }
  )
  plot.data$dim = factor(plot.data$dim, levels = plot.data$dim)

  plot =
    ggplot(plot.data, aes(x = dim, y = prop_var)) +
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
  return(plot)
}


#' Ggplot Beta-Diversity Ordination
#' @export
ggplot_beta_dispersion = function(
  beta.dispersion.fit, # output of get_beta_dispersion()
  hover.variables = NULL, # all columns of sample data if NULL
  scaling = 1, # see vegan's scalings 1 and 2
  comp = c(1, 2), # components to plot (loadings might be lacking components beyond 2)
  label.name = NULL,
  label.levels = NULL,
  facet.name = NULL, # single facet (backward compat alias for facet.row.name)
  facet.levels = NULL, # levels for facet.name
  facet.row.name = NULL, # row facet variable
  facet.row.levels = NULL,
  facet.col.name = NULL, # col facet variable; when set alongside facet.row.name, enables grid mode
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
  reverse.dim1 = FALSE,
  reverse.dim2 = FALSE
) {
  fit = beta.dispersion.fit$fit
  sample.data = beta.dispersion.fit$sample.data

  if (!is.null(facet.name) && is.null(facet.row.name)) {
    facet.row.name = facet.name
    facet.row.levels = facet.levels
  }
  grid.mode = !is.null(facet.row.name) && !is.null(facet.col.name)

  # for the hover info only:
  tax.table = beta.dispersion.fit$tax.table
  taxrank = beta.dispersion.fit$taxrank

  if (is.null(ellipses)) {
    ellipses = FALSE
  }

  # NOTE: species scores indicated by `vegan` = loadings here
  #       covariates (predictors) scores = covariates (constrained models only)
  # -> coordinates of arrows for biplot

  # NOTE: Comp (component) and Dim (dimension) are used interchangeably here
  if (is.null(biplot.loadings)) {
    biplot.loadings = FALSE
  }
  if (is.null(biplot.covariates)) {
    biplot.covariates = FALSE
  }
  biplot = biplot.loadings | biplot.covariates

  # If no coordinates provided
  if (length(beta.dispersion.fit$coords) == 0) {
    coords = NULL
    # If only one scaling provided
  } else if (length(beta.dispersion.fit$coords) == 1) {
    coords = beta.dispersion.fit$coords[[1]]
    # Select the scaling
  } else if (0 < scaling & scaling <= length(beta.dispersion.fit$coords)) {
    coords = beta.dispersion.fit$coords[[scaling]]
  } else {
    warning("Wrong scaling, scaling 1 is used")
    coords = beta.dispersion.fit$coords[[1]]
    scaling = 1
  }

  # Same as for coords
  if (length(beta.dispersion.fit$loadings) == 0) {
    loadings = NULL
  } else if (length(beta.dispersion.fit$loadings) == 1) {
    loadings = beta.dispersion.fit$loadings[[1]]
  } else if (0 < scaling & scaling <= length(beta.dispersion.fit$loadings)) {
    loadings = beta.dispersion.fit$loadings[[scaling]]
  } else {
    warning("Wrong scaling, scaling 1 is used")
    loadings = beta.dispersion.fit$loadings[[1]]
    scaling = 1
  }

  # Same as for coords
  if (length(beta.dispersion.fit$covariates) == 0) {
    covariates = NULL
  } else if (length(beta.dispersion.fit$covariates) == 1) {
    covariates = beta.dispersion.fit$covariates[[1]]
  } else if (0 < scaling & scaling <= length(beta.dispersion.fit$covariates)) {
    covariates = beta.dispersion.fit$covariates[[scaling]]
  } else {
    warning("Wrong scaling, scaling 1 is used")
    covariates = beta.dispersion.fit$covariates[[1]]
    scaling = 1
  }

  if (!is.null(beta.dispersion.fit$eigen.values)) {
    prop.var.explained = beta.dispersion.fit$eigen.values /
      sum(beta.dispersion.fit$eigen.values[
        beta.dispersion.fit$eigen.values > 0
      ]) *
      100
  } else {
    prop.var.explained = NULL
  }

  dist = beta.dispersion.fit$dist
  method = beta.dispersion.fit$method
  model = beta.dispersion.fit$model

  if (
    !is.null(label.name) &&
      !is.null(label.levels) &&
      !is.numeric(sample.data[[label.name]])
  ) {
    sample.data[[label.name]] = factor(
      sample.data[[label.name]],
      levels = label.levels
    )
  }

  if (
    !is.null(shape.name) &&
      !is.null(shape.levels) &&
      !is.numeric(sample.data[[shape.name]])
  ) {
    sample.data[[shape.name]] = factor(
      sample.data[[shape.name]],
      levels = shape.levels
    )
  }

  if (!is.null(facet.row.name) && !is.numeric(sample.data[[facet.row.name]])) {
    if (!is.null(facet.row.levels)) {
      sample.data[[facet.row.name]] = factor(
        sample.data[[facet.row.name]],
        levels = facet.row.levels
      )
    } else {
      sample.data[[facet.row.name]] = as.factor(sample.data[[facet.row.name]])
    }
  }

  if (
    grid.mode &&
      !is.null(facet.col.name) &&
      !is.numeric(sample.data[[facet.col.name]])
  ) {
    if (!is.null(facet.col.levels)) {
      sample.data[[facet.col.name]] = factor(
        sample.data[[facet.col.name]],
        levels = facet.col.levels
      )
    } else {
      sample.data[[facet.col.name]] = as.factor(sample.data[[facet.col.name]])
    }
  }

  if (
    !is.null(animation.variable.name) &&
      !is.null(animation.variable.levels) &&
      !is.numeric(sample.data[[animation.variable.name]])
  ) {
    sample.data[[animation.variable.name]] = factor(
      sample.data[[animation.variable.name]],
      levels = animation.variable.levels
    )
  }
  # (split by newline character to get one line of text for the subtitle)
  # if no label name is provided, stats are ignored
  p.value = NULL
  p.value.df = NULL

  if (!is.null(stat.beta.dispersion$label.name)) {
    if (grid.mode && !is.null(stat.beta.dispersion$p.value.df)) {
      # Grid mode: p-values come as a data frame, rendered as geom_text annotations
      p.value.df = stat.beta.dispersion$p.value.df
    } else if (!is.null(stat.beta.dispersion$p.value)) {
      if (
        is.null(facet.row.name) && is.null(names(stat.beta.dispersion$p.value))
      ) {
        # No facets: embed single p-value in subtitle
        p.value = paste(
          strsplit(stat.beta.dispersion$p.value[[1]], "\n")[[1]],
          collapse = " "
        )
      } else if (
        !is.null(facet.row.name) &&
          !is.null(names(stat.beta.dispersion$p.value)) &&
          all(
            levels(sample.data[[facet.row.name]]) %in%
              names(stat.beta.dispersion$p.value)
          )
      ) {
        # Single-facet: embed per-facet p-value into facet strip labels
        p.value = stat.beta.dispersion$p.value[levels(sample.data[[
          facet.row.name
        ]])]
      }
    }
  }

  # Remove all NA's from plot data (labels, facets or shape) by removing samples
  # having NA for at least one of graphical parameters (label, shape, facet, animation)
  # remove these samples from sample data AND from corresponding coordinates!
  if (remove.na.from.plot) {
    samples.wo.na = rep(TRUE, nrow(sample.data))
    for (var.name in c(
      label.name,
      shape.name,
      size.name,
      facet.row.name,
      facet.col.name,
      animation.variable.name
    )) {
      if (!is.null(sample.data[[var.name]])) {
        samples.wo.na = samples.wo.na & !is.na(sample.data[[var.name]])
      }
    }
    sample.data = sample.data[samples.wo.na, ]
    coords = coords[samples.wo.na, ]
  }

  # If stats on group are provided, overwrite the label.name by the one of this group
  if (!is.null(stat.beta.dispersion$label.name)) {
    label.name = stat.beta.dispersion$label.name
  }

  if (!is.null(label.name)) {
    label = sample.data[[label.name]]
  } else {
    label = NULL
  }

  if (!is.null(animation.variable.name)) {
    animation.variable = sample.data[[animation.variable.name]]
  } else {
    animation.variable = NULL
  }

  if (is.character(label)) {
    label = factor(label, levels = unique(label))
  }

  if (!is.null(colnames(coords))) {
    dim.names = colnames(coords)
  } else {
    dim.names = paste0("Dim ", 1:ncol(coords))
  }

  # Check if components are valid
  if (!is.numeric(comp) | length(comp) != 2 | !all(comp %in% 1:ncol(coords))) {
    warning("Wrong components, forced to first two components")
    comp = c(1, 2)
  }

  # Check if there are loadings and covariates for given components
  if (!is.null(loadings)) {
    if (!all(comp %in% 1:ncol(loadings))) {
      warning("No species score provided for these components")
      loadings = NULL
    }
  }

  if (!is.null(covariates)) {
    if (!all(comp %in% 1:ncol(covariates))) {
      warning("No species score provided for these components")
      covariates = NULL
    }
  }

  plot.df = data.frame(Comp1 = coords[, comp[1]], Comp2 = coords[, comp[2]])

  if (!is.null(label)) {
    plot.df$label = label
  }

  if (!is.null(animation.variable)) {
    plot.df[[animation.variable.name]] = animation.variable
  }

  if (grid.mode) {
    facet.row = sample.data[[facet.row.name]]
    facet.col = sample.data[[facet.col.name]]
    plot.df$facet.row = factor(
      paste0(facet.row.name, " = ", as.character(facet.row)),
      levels = paste0(facet.row.name, " = ", as.character(levels(facet.row)))
    )
    plot.df$facet.col = factor(
      paste0(facet.col.name, " = ", as.character(facet.col)),
      levels = paste0(facet.col.name, " = ", as.character(levels(facet.col)))
    )
    facet = NULL
  } else if (!is.null(facet.row.name)) {
    facet = sample.data[[facet.row.name]]
    if (is.character(facet) | is.factor(facet)) {
      if (!is.factor(facet)) {
        facet = as.factor(facet)
      }
      plot.df$facet = factor(
        paste0(
          paste0(facet.row.name, " = "),
          as.character(facet),
          if (length(p.value) > 1) {
            paste0("\n", p.value[facet])
          }
        ),
        levels = paste0(
          paste0(facet.row.name, " = "),
          as.character(levels(facet)),
          if (length(p.value) > 1) {
            paste0("\n", p.value[levels(facet)])
          }
        )
      )
    } else {
      facet = NULL
    }
  } else {
    facet = NULL
  }

  if (!is.null(shape.name)) {
    shape = sample.data[[shape.name]]
    if (is.character(shape) | is.factor(shape)) {
      if (!is.factor(shape)) {
        shape = as.factor(shape)
      }
      plot.df$shape = shape
    } else {
      shape = NULL
    }
  } else {
    shape = NULL
  }

  if (!is.null(size.name)) {
    size = sample.data[[size.name]]
  } else {
    size = NULL
  }

  hover.variables = colnames(sample.data)[
    colnames(sample.data) %in% hover.variables
  ]

  hover.text.all = rep(NA, nrow(sample.data))
  for (i in 1:nrow(sample.data)) {
    sample = rownames(sample.data)[i]
    hover.text = ""
    values = sample.data[sample, , drop = FALSE] %>% as("data.frame")
    for (variable in hover.variables) {
      value = values[[variable]]
      hover.text = paste0(hover.text, variable, ": ", value, "<br>")
    }
    hover.text.all[i] = hover.text
  }

  plot.df$hover.text = hover.text.all

  plt = ggplot() +
    geom_point(
      aes(
        x = Comp1,
        y = Comp2,
        color = label, # OK if NULL
        text = hover.text,
        shape = shape, # OK if NULL
        size = size # OK if NULL
      ),
      plot.df,
      alpha = point.alpha
    ) +
    labs(color = label.name, shape = shape.name, size = size.name)

  if (ellipses & is.factor(label)) {
    if (fill.ellipses) {
      plt = plt +
        stat_ellipse(
          aes(x = Comp1, y = Comp2, color = label, fill = label),
          plot.df,
          geom = "polygon",
          alpha = 0.2
        )
    } else {
      plt = plt +
        stat_ellipse(aes(x = Comp1, y = Comp2, color = label), plot.df)
    }
  }

  plt = plt +
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
    scale.factor.load = NULL
    scale.factor.covar = NULL

    arrow.df = data.frame()

    # Plot species arrows = loadings here
    if (!is.null(loadings) & biplot.loadings) {
      # Account for absent loadings
      filtered.loadings = na.omit(loadings[, comp, drop = FALSE]) # TODO: how come that some loadings are NaN?
      #       is it for all models?
      # Filter arrows by length (get rid of smaller arrows)
      arrow.lengths = apply(filtered.loadings, 1, function(x) sqrt(sum(x^2)))
      filtered.loadings = filtered.loadings[
        arrow.lengths / max(arrow.lengths) > arrow.cutoff.load,
        ,
        drop = FALSE
      ]

      # If there are loadings left...
      if (nrow(filtered.loadings) > 0) {
        arrow.load = data.frame(
          Comp1 = filtered.loadings[, 1],
          Comp2 = filtered.loadings[, 2],
          type = "load"
        )

        rownames(arrow.load) = row.names(filtered.loadings)

        # Scale the secondary axes for loadings so that loadings and scores are on
        # comparable scale
        # NOTE: both components are scaled by the same factor, so arrows point
        # in the same direction as before, only their size changes!
        lim.load = max(abs(arrow.load$Comp1), abs(arrow.load$Comp2))
        scale.factor.load = max(abs(coords)) / lim.load
        arrow.load$Comp1 = arrow.load$Comp1 * scale.factor.load
        arrow.load$Comp2 = arrow.load$Comp2 * scale.factor.load

        # If phyloseq object provided, get hover info from it - taxonomy for each species (arrow)
        if (!is.null(tax.table)) {
          hover.text.all = rep(NA, nrow(arrow.load))
          for (i in 1:nrow(arrow.load)) {
            # it makes sure that taxa names match
            taxon = rownames(arrow.load)[i]
            hover.text = ""
            values = tax.table[taxon, , drop = FALSE]
            for (variable in colnames(values)) {
              value = values[[variable]]
              hover.text = paste0(hover.text, variable, ": ", value, "<br>")
            }
            hover.text.all[i] = hover.text
          }

          arrow.load$hover.text = hover.text.all
        } else {
          arrow.load$hover.text = rownames(arrow.load)
        }

        if (!is.null(arrow.taxonomy.labels)) {
          # use taxonomy as label instead of taxa names
          arrow.load$Names = apply(
            tax.table[rownames(arrow.load), ],
            1,
            function(row) {
              gsub("NA", " ", paste(row[arrow.taxonomy.labels], collapse = ":"))
            }
          )
        } else {
          arrow.load$Names = rownames(arrow.load)
        }
        arrow.df = rbind(arrow.df, arrow.load)
      }
    }
    # Similar for covariates (predictors in case of constrained model)
    if (!is.null(covariates) & biplot.covariates) {
      # Filter out absent data or short arrows
      filtered.covariates = na.omit(covariates[, comp, drop = FALSE])
      arrow.lengths = apply(filtered.covariates, 1, function(x) sqrt(sum(x^2)))
      filtered.covariates = filtered.covariates[
        arrow.lengths / max(arrow.lengths) > arrow.cutoff.covar,
        ,
        drop = FALSE
      ]
      # If some are left...
      if (nrow(filtered.covariates) > 0) {
        arrow.covar = data.frame(
          Comp1 = filtered.covariates[, 1],
          Comp2 = filtered.covariates[, 2],
          type = "covar"
        )
        rownames(arrow.covar) = row.names(filtered.covariates)

        # Scale in the same manner as for loading but notice that they are scaled
        # by a different factor, so that all are of comparable size on plot
        lim.covar = max(abs(arrow.covar$Comp1), abs(arrow.covar$Comp2))
        scale.factor.covar = max(abs(coords)) / lim.covar
        arrow.covar$Comp1 = arrow.covar$Comp1 * scale.factor.covar
        arrow.covar$Comp2 = arrow.covar$Comp2 * scale.factor.covar
        # Hover text is just the name of covariate
        arrow.covar$hover.text = rownames(arrow.covar)
        arrow.covar$Names = gsub("`", "", rownames(arrow.covar))
        arrow.df = rbind(arrow.df, arrow.covar)
      }
    }

    # If there are loadings and/or covariates, plot arrows
    if (nrow(arrow.df) > 0) {
      arrow.params = arrow(
        type = "closed",
        angle = 20,
        length = unit(0.1, "inches")
      )

      if ("load" %in% arrow.df$type) {
        if (color.arrows.by.taxa) {
          plt = plt +
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
          plt = plt +
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
        plt = plt +
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
          plt = plt +
            ggrepel::geom_text_repel(
              data = arrow.df,
              aes(x = Comp1, y = Comp2, label = Names),
              max.overlaps = max.overlaps,
              show.legend = FALSE
            )
        } else {
          plt = plt +
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
        scale.factor = c(scale.factor.load, scale.factor.covar) # = the one which is not NULL
        plt = plt +
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
    plt = plt +
      facet_grid(facet.row ~ facet.col, scales = "fixed")

    if (!is.null(p.value.df) && nrow(p.value.df) > 0) {
      annot.df = data.frame(
        facet.row = factor(
          paste0(facet.row.name, " = ", p.value.df$facet.row),
          levels = levels(plot.df$facet.row)
        ),
        facet.col = factor(
          paste0(facet.col.name, " = ", p.value.df$facet.col),
          levels = levels(plot.df$facet.col)
        ),
        p.label = p.value.df$p.label,
        x = Inf,
        y = Inf
      )
      plt = plt +
        geom_text(
          data = annot.df,
          aes(x = x, y = y, label = p.label),
          hjust = 1.05,
          vjust = 1.5,
          size = 3,
          inherit.aes = FALSE
        )
    }
  } else if (!is.null(facet)) {
    plt = plt +
      facet_wrap(
        . ~ facet,
        scales = "fixed",
        ncol = smart_facet_ncol(nlevels(plot.df$facet))
      )
  }

  # Marginal plot in case of factor variable - boxplot, density...
  if (!is.null(marginal.plot) & is.factor(label) & is.null(facet)) {
    plt = plt +
      theme(
        legend.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "left"
      ) # place legend to the left so that it doesn't interfere with the marginal plot
    plt = ggExtra::ggMarginal(
      plt,
      type = marginal.plot,
      groupColour = TRUE,
      groupFill = TRUE
    )
  } else if (!is.null(label)) {
    plt = plt +
      theme(legend.title = element_text(face = "bold", hjust = 0.5))
  }

  if (reverse.dim1) {
    plt = plt + scale_x_reverse()
  }

  if (reverse.dim2) {
    plt = plt + scale_y_reverse()
  }

  return(plt)
}


# Inspired from animateDimReduction() in Feature Selector
# TODO: this function is better since it doesn't assume any facets for animation,
# implement it in feature selector
#' Animate a Ggplot Object by a Variable
#' @export
animate_by_variable = function(
  ggplot.obj, # ggplot2 object
  animation.variable.name,
  return.anim = TRUE,
  save.path = NULL # gif file name
) {
  plt = ggplot.obj

  plt =
    plt +
    gganimate::transition_states(
      get(animation.variable.name),
      transition_length = 2,
      state_length = 3
    ) +
    labs(subtitle = paste0(animation.variable.name, " = {closest_state}")) +
    ease_aes('linear')
  anim =
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
