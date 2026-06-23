#' Prepare Data for Heatmap from Phyloseq Object
#' @export
get_phylo_heatmap = function(
  physeq = physeq,
  fraction_id_name = NULL,
  fraction_ids = NULL,
  # user might want to use different agglomeration
  # for visualization on heatmap and for clustering samples
  taxrank.for.heatmap = NULL,
  taxrank.for.hclust = NULL,
  transform.abundances = "identity",
  distance = "bray", # for hclust
  vars.to.remove.na = c() # a vector of variable names for which no NA's allowed
  # done here and not at plotting level because dendrogram has
  # to be be redone if some samples are excluded
) {
  if (class(physeq) != "phyloseq") {
    stop("Need a phyloseq object")
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

  if (is.null(distance)) {
    print("no distance provided, bray-curtis selected by default")
    distance = "bray"
  }
  if (
    is.null(access(physeq, "phy_tree")) & distance %in% c("unifrac", "wunifrac")
  ) {
    # TODO: check the spelling for unifrac
    print(paste0(
      "Uning ",
      distance,
      "distance requires phylogenetic tree in phyloseq object, bray-curtis distance selected instead."
    ))
    distance = "bray"
  }

  # Remove all NA's from plot data by removing samples having NA for at least one
  # of graphical parameters ( heatmap annotation variables eventually applied for plotting)
  # NOTE: here, it will affect the analysis - hierarchical clustering
  samples.wo.na = rep(TRUE, nrow(sample_data(physeq)))
  for (var.name in vars.to.remove.na) {
    if (!is.null(sample_data(physeq)[[var.name]])) {
      samples.wo.na = samples.wo.na & !is.na(sample_data(physeq)[[var.name]])
    }
  }
  samples.wo.na.global <<- samples.wo.na
  physeq = subset_samples(physeq, samples.wo.na.global)

  # Clustering of samples is performed on the original (agglomerated and transformed) data,
  # taxrank for agglomeration can be different for clustering (to balance better between noise and information)
  # and not to affect the heatmap visualisation
  dist_mat = sparse_distance(
    physeq %>%
      {
        if (!is.null(taxrank.for.hclust)) {
          tax_glom(., taxrank = taxrank.for.hclust)
        } else {
          .
        }
      } %>%
      {
        if (
          !is.null(transform.abundances) &
            transform.abundances != "identity"
        ) {
          microbiome::transform(., transform = transform.abundances)
        } else {
          .
        }
      },
    method = distance
  )
  cluster.fit = hclust(
    dist_mat,
    "ward.D2"
  )

  dendrogram = as.dendrogram(cluster.fit)

  if (!is.null(taxrank.for.heatmap)) {
    physeq = tax_glom(physeq, taxrank = taxrank.for.heatmap)
    taxa_names(physeq) = make.unique(tax_table(physeq)[, taxrank.for.heatmap]) # NOTE: names of taxa are not necessarily unique!
  }

  # Sort by abundance before applying the transformation
  # (most abundant will be on top)
  taxa.sorted.by.abundance = names(sort(taxa_sums(physeq), decreasing = TRUE))

  # NOTE: agglomerate taxa BEFORE transforming the data
  if (!is.null(transform.abundances) & transform.abundances != "identity") {
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

  heat.matrix = as(otu_table(physeq), "matrix")

  # Taxa must be rows for heatmap
  if (!taxa_are_rows(physeq)) {
    heat.matrix = t(heat.matrix)
    print("abundance matrix is reversed so that taxa are rows")
  }

  row.names(heat.matrix) = taxa_names(physeq)

  return(list(
    heat.matrix = heat.matrix,
    taxa.sorted.by.abundance = taxa.sorted.by.abundance,
    sample.data = sample_data(physeq) %>% as("data.frame"),
    dendrogram = dendrogram
  ))
}


#' Plot ComplexHeatmap
#' @export
plot_phylo_heatmap = function(
  heatmap.data,
  sort_taxa_by_diff_abundance = FALSE, # by abundance otherwise by variation
  # or by effect size (Kruskal Wallis) - group separation
  # TODO: if variation, decide between sd and CV (coefficient of variation)
  var_for_diff_abundance = NULL, # if NULL
  scale_cols = FALSE, # scale by column (sample) after selecting top taxa
  nb_top_taxa = 30,
  top_annotation_vars = NULL, # can be multiple variables in a vector
  bottom_annotation_var = NULL, # one variable
  split = NULL,
  color_vector = c("white", "#88CCAA", "#771122"),
  taxa_names_par = gpar(fontsize = 15, fontface = "bold.italic"),
  label_names_par = gpar(fontsize = 15, fontface = "bold")
) {
  # TODO: add option to display number of reads
  heat.matrix = heatmap.data$heat.matrix # rows are taxa (sorted by abundance), columns are samples
  sample.data = heatmap.data$sample.data
  dendrogram = heatmap.data$dendrogram

  if (sort_taxa_by_diff_abundance) {
    if (is.null(var_for_diff_abundance)) {
      print("taxa sorted by coefficient of variation")
      heat.matrix = heat.matrix[
        names(sort(
          apply(heat.matrix, 1, function(row) {
            sd(row) / mean(row)
          }),
          decreasing = TRUE
        )),
      ]
    } else {
      print("taxa sorted by effect size")
      var.for.diff = sample.data[[var_for_diff_abundance]]
      if (is.character(var.for.diff) | is.factor(var.for.diff)) {
        # sort by effect size of Kruskal Wallis test
        heat.matrix = heat.matrix[
          names(sort(
            apply(heat.matrix, 1, function(row) {
              H = kruskal.test(row ~ var.for.diff)$statistic
              k = length(unique(var.for.diff))
              n = length(var.for.diff)
              eps2 = (H - k + 1) / (n - k)
              return(eps2)
            }),
            decreasing = TRUE
          )),
        ]
      } else if (is.numeric(var.for.diff)) {
        # sort by effect size of Spearman correlation test
        heat.matrix = heat.matrix[
          names(sort(
            apply(heat.matrix, 1, function(row) {
              rho = cor.test(row, var.for.diff, method = "spearman")$estimate

              return(abs(rho))
            }),
            decreasing = TRUE
          )),
        ]
      }
    }
  } else {
    print("taxa sorted by abundance")
    heat.matrix = heat.matrix[heatmap.data$taxa.sorted.by.abundance, ]
  }

  nb_top_taxa = min(nb_top_taxa, nrow(heat.matrix))
  heat.matrix = heat.matrix[1:nb_top_taxa, ]

  # Scale by column (sample) - optional
  # TODO: does it even make sense if abundances were transformed?
  if (scale_cols) {
    heat.matrix = scale(heat.matrix)
  }

  # Variable(s) to put at the top of the heatmap
  args.top = list()
  for (var in top_annotation_vars) {
    if (var %in% names(sample.data)) {
      args.top[[var]] = sample.data[[var]]
    }
  }
  if (length(args.top) > 0) {
    args.top = c(args.top, list(show_annotation_name = TRUE))
    top.annotation = do.call(HeatmapAnnotation, args.top)
  } else {
    top.annotation = NULL
  }

  # Variable to use for bottom annotation
  if (
    !is.null(bottom_annotation_var) &&
      bottom_annotation_var %in% names(sample.data)
  ) {
    bottom.annotation = HeatmapAnnotation(
      Name = anno_text(
        sample.data[[bottom_annotation_var]],
        gp = gpar(fontsize = 9)
      )
    )
  } else {
    bottom.annotation = NULL
  }

  heat.map = Heatmap(
    heat.matrix,
    cluster_rows = FALSE,
    cluster_columns = dendrogram,
    top_annotation = top.annotation,
    bottom_annotation = bottom.annotation,
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

  return(heat.map)
}
