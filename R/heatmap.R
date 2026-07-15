#' Prepare Data for a Phylogenetic Heatmap
#'
#' Builds the abundance matrix, sample metadata, taxa ordering, and sample
#' dendrogram consumed by [plot_phylo_heatmap()].
#'
#' @param physeq A `phyloseq` object.
#' @param fraction_id_name Name of a `sample_data` variable used, together
#'   with `fraction_ids`, to keep only a subset of samples (e.g. a
#'   sorting-fraction indicator). `NULL` (default) keeps all samples.
#' @param fraction_ids Values of `fraction_id_name` to keep. Ignored if
#'   `fraction_id_name` is `NULL`.
#' @param taxrank_for_heatmap Taxonomic rank to agglomerate to (via
#'   [tax_glom()]) before building the heatmap matrix. `NULL` (default) keeps
#'   the original taxa.
#' @param taxrank_for_hclust Taxonomic rank to agglomerate to before
#'   computing the sample distance/dendrogram; independent of
#'   `taxrank_for_heatmap` so clustering and display can use different
#'   resolutions. `NULL` (default) keeps the original taxa.
#' @param transform_abundances Abundance transform applied before both the
#'   heatmap matrix and the clustering distance are computed, forwarded to
#'   [microbiome::transform()] (e.g. `"clr"`, `"hellinger"`). `"identity"`
#'   (default) applies no transform.
#' @param distance Distance method forwarded to [sparse_distance()] for the
#'   sample dendrogram. Falls back to `"bray"` (with a message) if `NULL`, or
#'   if a UniFrac-family method is requested without a `phy_tree` in
#'   `physeq`.
#' @param vars_to_remove_na Character vector of `sample_data` variable names;
#'   samples with `NA` in any of them are dropped before clustering. This is
#'   done here rather than at plot time because the dendrogram would
#'   otherwise need to be rebuilt whenever samples are excluded.
#'
#' @return A list with:
#'   \item{heat_matrix}{Abundance matrix, taxa as rows, samples as columns.}
#'   \item{taxa_sorted_by_abundance}{Taxa names ordered by decreasing total
#'     abundance.}
#'   \item{sample_data}{Sample metadata as a `data.frame`.}
#'   \item{dendrogram}{A [stats::as.dendrogram()] object from clustering
#'     samples.}
#'
#' @examples
#' data(ps_16s_refinement)
#' heatmap_data <- get_phylo_heatmap(ps_16s_refinement)
#' dim(heatmap_data$heat_matrix)
#' length(heatmap_data$taxa_sorted_by_abundance)
#'
#' @export
get_phylo_heatmap <- function(
  physeq,
  fraction_id_name = NULL,
  fraction_ids = NULL,
  taxrank_for_heatmap = NULL,
  taxrank_for_hclust = NULL,
  transform_abundances = "identity",
  distance = "bray",
  vars_to_remove_na = c()
) {
  if (!is(physeq, "phyloseq")) {
    stop("Need a phyloseq object")
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

  distance <- .resolve_heatmap_distance(physeq, distance)

  # Remove all NA's from plot data by removing samples having NA for at least one
  # of graphical parameters ( heatmap annotation variables eventually applied for plotting)
  # NOTE: here, it will affect the analysis - hierarchical clustering
  samples_wo_na <- rep(TRUE, nrow(sample_data(physeq)))
  for (var_name in vars_to_remove_na) {
    if (!is.null(sample_data(physeq)[[var_name]])) {
      samples_wo_na <- samples_wo_na & !is.na(sample_data(physeq)[[var_name]])
    }
  }
  physeq <- prune_samples(samples_wo_na, physeq)

  # Clustering of samples is performed on the original (agglomerated and transformed) data,
  # taxrank for agglomeration can be different for clustering (to balance better between noise and information)
  # and not to affect the heatmap visualisation
  dendrogram <- .compute_sample_dendrogram(
    physeq,
    taxrank_for_hclust = taxrank_for_hclust,
    transform_abundances = transform_abundances,
    distance = distance
  )

  heat <- .build_heat_matrix(
    physeq,
    taxrank_for_heatmap = taxrank_for_heatmap,
    transform_abundances = transform_abundances
  )

  list(
    heat_matrix = heat$heat_matrix,
    taxa_sorted_by_abundance = heat$taxa_sorted_by_abundance,
    sample_data = as(sample_data(physeq), "data.frame"),
    dendrogram = dendrogram
  )
}

#' Resolve the Heatmap Distance Method
#'
#' Falls back to `"bray"` (with a message) when `distance` is `NULL`, or when
#' a UniFrac-family method is requested but `physeq` has no `phy_tree`.
#'
#' @return A single distance-method string.
#' @noRd
.resolve_heatmap_distance <- function(physeq, distance) {
  if (is.null(distance)) {
    message("No distance provided, bray-curtis selected by default.")
    return("bray")
  }
  if (
    is.null(access(physeq, "phy_tree")) && distance %in% c("unifrac", "wunifrac")
  ) {
    message(
      "Using \"", distance, "\" distance requires a phylogenetic tree in the ",
      "phyloseq object; bray-curtis distance selected instead."
    )
    return("bray")
  }
  distance
}

#' Compute the Sample Dendrogram for a Heatmap
#'
#' Agglomerates/transforms `physeq` as requested, then hierarchically
#' clusters samples (`ward.D2`) on the resulting distance matrix.
#'
#' @return A [stats::as.dendrogram()] object.
#' @noRd
.compute_sample_dendrogram <- function(
  physeq,
  taxrank_for_hclust,
  transform_abundances,
  distance
) {
  if (!is.null(taxrank_for_hclust)) {
    physeq <- tax_glom(physeq, taxrank = taxrank_for_hclust)
  }
  if (!is.null(transform_abundances) && transform_abundances != "identity") {
    physeq <- microbiome::transform(physeq, transform = transform_abundances)
  }

  dist_mat <- sparse_distance(physeq, method = distance)
  cluster_fit <- hclust(dist_mat, method = "ward.D2")
  as.dendrogram(cluster_fit)
}

#' Build the Heatmap Abundance Matrix
#'
#' Agglomerates to `taxrank_for_heatmap` (if given), records the
#' pre-transform abundance ordering, then applies `transform_abundances` and
#' extracts a taxa-as-rows matrix.
#'
#' @return A list with `heat_matrix` and `taxa_sorted_by_abundance`.
#' @noRd
.build_heat_matrix <- function(physeq, taxrank_for_heatmap, transform_abundances) {
  if (!is.null(taxrank_for_heatmap)) {
    physeq <- tax_glom(physeq, taxrank = taxrank_for_heatmap)
    # NOTE: names of taxa are not necessarily unique!
    taxa_names(physeq) <- make.unique(tax_table(physeq)[, taxrank_for_heatmap])
  }

  # Sort by abundance before applying the transformation (most abundant will be on top)
  taxa_sorted_by_abundance <- names(sort(taxa_sums(physeq), decreasing = TRUE))

  # NOTE: agglomerate taxa BEFORE transforming the data
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

  heat_matrix <- as(otu_table(physeq), "matrix")

  # Taxa must be rows for heatmap
  if (!taxa_are_rows(physeq)) {
    heat_matrix <- t(heat_matrix)
    message("Abundance matrix transposed so that taxa are rows.")
  }
  rownames(heat_matrix) <- taxa_names(physeq)

  list(heat_matrix = heat_matrix, taxa_sorted_by_abundance = taxa_sorted_by_abundance)
}

#' Plot a Phylogenetic Heatmap
#'
#' Renders the abundance matrix prepared by [get_phylo_heatmap()] as a
#' [ComplexHeatmap::Heatmap()], with samples clustered by the precomputed
#' dendrogram and taxa ordered/limited to the top `nb_top_taxa`.
#'
#' @param heatmap_data A list as returned by [get_phylo_heatmap()].
#' @param sort_taxa_by_diff_abundance If `TRUE`, order taxa by variation
#'   (coefficient of variation) or by effect size against
#'   `var_for_diff_abundance`, instead of by total abundance.
#' @param var_for_diff_abundance Name of a `sample_data` variable to sort
#'   taxa by effect size against: Kruskal-Wallis epsilon-squared for a
#'   character/factor variable, absolute Spearman correlation for a numeric
#'   one. `NULL` (default) sorts by coefficient of variation instead. Ignored
#'   if `sort_taxa_by_diff_abundance` is `FALSE`.
#' @param scale_cols If `TRUE`, scale (z-score) each sample (column) after
#'   selecting the top taxa.
#' @param nb_top_taxa Number of top-ranked taxa (per the sort above) to keep.
#' @param top_annotation_vars Character vector of `sample_data` variable
#'   names to show as top column annotations.
#' @param bottom_annotation_var Name of a single `sample_data` variable to
#'   show as a bottom text annotation.
#' @param split If a number `>= 2`, the number of clusters to split columns
#'   into (`column_split`); `NULL`/`< 2` disables splitting.
#' @param color_vector Vector of colors forwarded to `Heatmap(col = ...)`.
#' @param taxa_names_par A [grid::gpar()] object for row (taxon) name text.
#' @param label_names_par A [grid::gpar()] object for annotation label text.
#'
#' @return A [ComplexHeatmap::Heatmap-class] object.
#'
#' @examples
#' data(ps_16s_refinement)
#' heatmap_data <- get_phylo_heatmap(ps_16s_refinement)
#' heatmap <- plot_phylo_heatmap(heatmap_data, nb_top_taxa = 10)
#' methods::is(heatmap, "Heatmap")
#'
#' @export
plot_phylo_heatmap <- function(
  heatmap_data,
  sort_taxa_by_diff_abundance = FALSE,
  var_for_diff_abundance = NULL,
  scale_cols = FALSE,
  nb_top_taxa = 30,
  top_annotation_vars = NULL,
  bottom_annotation_var = NULL,
  split = NULL,
  color_vector = c("white", "#88CCAA", "#771122"),
  taxa_names_par = gpar(fontsize = 15, fontface = "bold.italic"),
  label_names_par = gpar(fontsize = 15, fontface = "bold")
) {
  # TODO: add option to display number of reads
  heat_matrix <- heatmap_data$heat_matrix # rows are taxa (sorted by abundance), columns are samples
  sample_data <- heatmap_data$sample_data
  dendrogram <- heatmap_data$dendrogram

  taxa_order <- if (sort_taxa_by_diff_abundance) {
    .sort_taxa_by_diff_abundance(heat_matrix, sample_data, var_for_diff_abundance)
  } else {
    message("Taxa sorted by abundance.")
    heatmap_data$taxa_sorted_by_abundance
  }
  heat_matrix <- heat_matrix[taxa_order, , drop = FALSE]

  nb_top_taxa <- min(nb_top_taxa, nrow(heat_matrix))
  heat_matrix <- heat_matrix[seq_len(nb_top_taxa), , drop = FALSE]

  # Scale by column (sample) - optional
  # TODO: does it even make sense if abundances were transformed?
  if (scale_cols) {
    heat_matrix <- scale(heat_matrix)
  }

  top_annotation <- .build_top_annotation(sample_data, top_annotation_vars)
  bottom_annotation <- .build_bottom_annotation(sample_data, bottom_annotation_var)

  Heatmap(
    heat_matrix,
    cluster_rows = FALSE,
    cluster_columns = dendrogram,
    top_annotation = top_annotation,
    bottom_annotation = bottom_annotation,
    row_names_gp = gpar(taxa_names_par),
    show_column_names = FALSE,
    col = color_vector,
    column_dend_height = unit(3, "cm"),
    name = "Abundance",
    column_split = if (is.null(split) || split < 2) {
      NULL
    } else {
      split
    },
    column_gap = unit(5, "mm"),
    column_title = NULL
  )
}

#' Order Heatmap Taxa by Variation or Effect Size
#'
#' Dispatches to coefficient-of-variation ordering (`var_for_diff_abundance`
#' is `NULL`), Kruskal-Wallis effect-size ordering (character/factor
#' variable), or Spearman effect-size ordering (numeric variable).
#'
#' @return A character vector of taxa names, ordered by decreasing effect.
#' @noRd
.sort_taxa_by_diff_abundance <- function(heat_matrix, sample_data, var_for_diff_abundance) {
  if (is.null(var_for_diff_abundance)) {
    return(.cv_taxa_order(heat_matrix))
  }
  var_for_diff <- sample_data[[var_for_diff_abundance]]
  if (is.character(var_for_diff) || is.factor(var_for_diff)) {
    .kruskal_effect_size_order(heat_matrix, var_for_diff)
  } else if (is.numeric(var_for_diff)) {
    .spearman_effect_size_order(heat_matrix, var_for_diff)
  } else {
    stop(
      "`var_for_diff_abundance` (\"", var_for_diff_abundance, "\") must be a ",
      "character, factor, or numeric sample_data variable, not ",
      class(var_for_diff)[1], "."
    )
  }
}

#' Order Taxa by Coefficient of Variation
#'
#' Taxa whose mean abundance is 0 have an undefined coefficient of variation;
#' these are placed last and reported via `warning()` rather than silently
#' producing `Inf`/`NaN`.
#'
#' @return A character vector of taxa names, ordered by decreasing CV.
#' @noRd
.cv_taxa_order <- function(heat_matrix) {
  message("Taxa sorted by coefficient of variation.")
  cv <- apply(heat_matrix, 1, function(row) {
    row_mean <- mean(row)
    if (row_mean == 0) {
      return(NA_real_)
    }
    sd(row) / row_mean
  })
  if (any(is.na(cv))) {
    warning(
      "Coefficient of variation is undefined (mean abundance is 0) for taxa: ",
      paste(names(cv)[is.na(cv)], collapse = ", "),
      call. = FALSE
    )
  }
  names(sort(cv, decreasing = TRUE, na.last = TRUE))
}

#' Order Taxa by Kruskal-Wallis Effect Size
#'
#' @return A character vector of taxa names, ordered by decreasing
#'   epsilon-squared.
#' @noRd
.kruskal_effect_size_order <- function(heat_matrix, group_var) {
  message("Taxa sorted by effect size (Kruskal-Wallis).")
  eps2 <- apply(heat_matrix, 1, function(row) {
    h_stat <- kruskal.test(row ~ group_var)$statistic
    k <- length(unique(group_var))
    n <- length(group_var)
    (h_stat - k + 1) / (n - k)
  })
  names(sort(eps2, decreasing = TRUE))
}

#' Order Taxa by Spearman Effect Size
#'
#' @return A character vector of taxa names, ordered by decreasing absolute
#'   Spearman correlation.
#' @noRd
.spearman_effect_size_order <- function(heat_matrix, numeric_var) {
  message("Taxa sorted by effect size (Spearman correlation).")
  rho <- apply(heat_matrix, 1, function(row) {
    abs(cor.test(row, numeric_var, method = "spearman")$estimate)
  })
  names(sort(rho, decreasing = TRUE))
}

#' Build the Heatmap Top Annotation
#'
#' @return A [ComplexHeatmap::HeatmapAnnotation()], or `NULL` if none of
#'   `top_annotation_vars` are present in `sample_data`.
#' @noRd
.build_top_annotation <- function(sample_data, top_annotation_vars) {
  args_top <- list()
  for (var in top_annotation_vars) {
    if (var %in% names(sample_data)) {
      args_top[[var]] <- sample_data[[var]]
    }
  }
  if (length(args_top) == 0) {
    return(NULL)
  }
  args_top <- c(args_top, list(show_annotation_name = TRUE))
  do.call(HeatmapAnnotation, args_top)
}

#' Build the Heatmap Bottom Annotation
#'
#' @return A [ComplexHeatmap::HeatmapAnnotation()] text annotation, or `NULL`
#'   if `bottom_annotation_var` is `NULL` or absent from `sample_data`.
#' @noRd
.build_bottom_annotation <- function(sample_data, bottom_annotation_var) {
  if (is.null(bottom_annotation_var) || !bottom_annotation_var %in% names(sample_data)) {
    return(NULL)
  }
  HeatmapAnnotation(
    Name = anno_text(sample_data[[bottom_annotation_var]], gp = gpar(fontsize = 9))
  )
}
