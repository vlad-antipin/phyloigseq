# for a function foo() from ImmuMicrobiome package
# get_foo(): performs calculation and outputs data
# plot_foo(): builds plot based on the result of get_foo()
# full_foo(): combines get_foo() and plot_foo() to restore the original foo()

#' Sparse-aware Shannon alpha diversity
#'
#' Computes Shannon entropy per sample directly from the sparse matrix slot,
#' without materialising the full dense OTU table.
#'
#' @param physeq A \code{phyloseq} object whose OTU table is a
#'   \code{\link{sparse_otu_table-class}}.
#' @return A named numeric vector of Shannon entropy values (nats), one per sample.
#' @export
sparse_shannon <- function(physeq) {
  ot <- phyloseq::otu_table(physeq)
  sp <- ot@sparse_data
  if (phyloseq::taxa_are_rows(ot)) {
    sp <- Matrix::t(sp)
  }
  # sp: samples x taxa (rows = samples)
  N <- Matrix::rowSums(sp)
  sp_norm <- Matrix::Diagonal(x = 1 / N) %*% sp
  log_nz <- sp_norm
  log_nz@x <- log(sp_norm@x)
  h <- -Matrix::rowSums(sp_norm * log_nz)
  setNames(as.numeric(h), rownames(sp))
}

#' Get Alpha-Diversity from Phyloseq Object
#'
#' Identity (raw counts) is suitable for all metrics,
# compositional (proportions) are not suitable for Chao1, ACE and Fisher
# all other transformations are not suitable for alpha div. metrics
# CCL: use identity all the time
#' @export
get_alpha_diversity <- function(
  physeq,
  from_igseq = FALSE,
  transform.abundances = "identity",
  # if from Ig Seq
  proportions = FALSE,
  low_lim = -1.96,
  high_lim = 1.96,
  # if not from Ig Seq
  taxrank = NULL,
  fraction_id_name = NULL,
  fraction_ids = NULL,
  measure = "Shannon"
) {
  if (class(physeq) != "phyloseq") {
    stop("Need a phyloseq object")
  }

  # Force taxa to be columns of otu table - like in internal phyloseq function phyloseq:::veganifyOTU()
  if (taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }

  if (!from_igseq) {
    if (!is.null(fraction_id_name) & !is.null(fraction_ids)) {
      physeq <- prune_samples(
        sample_data(physeq)[[fraction_id_name]] %in% fraction_ids,
        physeq
      )
    }

    # NOTE: agglomerate taxa BEFORE transforming the data
    # Agglomerate taxa up to a certain taxrank
    if (!is.null(taxrank)) {
      physeq <- tax_glom(physeq = physeq, taxrank = taxrank)
      taxa_names(physeq) <- make.unique(tax_table(physeq)[, taxrank])
    }

    if (!is.null(transform.abundances) & transform.abundances != "identity") {
      physeq <- microbiome::transform(
        physeq,
        transform = transform.abundances,
        target = "OTU", # TODO: and still, clr will scale over samples
        shift = 0, # pseudocount added (shifts baseline)
        scale = 1, # if transform is "scale"
        log10 = TRUE,
        reference = 1
      )
    }

    # alpha diversity is in the first column
    if (
      is(otu_table(physeq), "sparse_otu_table") && identical(measure, "Shannon")
    ) {
      shannon_vals <- sparse_shannon(physeq)
      alpha.diversity <- data.frame(
        Shannon = shannon_vals,
        row.names = names(shannon_vals)
      )
    } else {
      alpha.diversity <- estimate_richness(physeq, measures = measure)
    }
    full.sample.data <- cbind(
      alpha.diversity,
      as(sample_data(physeq), "data.frame")
    )

    full.sample.data$depth <- sample_sums(physeq) # nreads for each sample
    rownames(full.sample.data) <- phyloseq::sample_names(physeq)
  } else {
    full.sample.data <-
      get_igseq_richness(
        ps_ig_score = physeq,
        proportions = proportions,
        low_lim = low_lim,
        high_lim = high_lim
      )
  }

  return(full.sample.data)
}

is_valid_factor <- function(df, name) {
  return(!is.null(name) && name %in% names(df) && !is.numeric(df[[name]]))
}

keep_levels <- function(vals, level_names) {
  if (!is.null(level_names) && length(level_names) > 0) {
    if ("(NA)" %in% level_names) {
      vals <- as.character(vals)
      vals[is.na(vals)] <- "(NA)"
    }
    as.character(vals) %in% level_names
  } else {
    rep(TRUE, length(vals))
  }
}

factorize_levels <- function(vals, level_names) {
  if (!is.null(level_names) && length(level_names) > 0) {
    if ("(NA)" %in% level_names) {
      vals <- as.character(vals)
      vals[is.na(vals)] <- "(NA)"
    }
    factor(vals, levels = level_names)
  } else {
    factor(vals, levels = gtools::mixedsort(unique(vals[!is.na(vals)])))
  }
}

apply_levels <- function(df, name, level_names) {
  df <- df[keep_levels(df[[name]], level_names), , drop = FALSE]
  df[[name]] <- factorize_levels(df[[name]], level_names)
  df
}

remove_nas <- function(df, var.names) {
  samples.wo.na <- rep(TRUE, nrow(df))
  for (var.name in var.names) {
    if (!is.null(df[[var.name]])) {
      samples.wo.na <- samples.wo.na & !is.na(df[[var.name]])
    }
  }
  return(df[samples.wo.na, ])
}

get_hover_text <- function(df, hover.variables) {
  hover.variables <- colnames(df)[colnames(df) %in% hover.variables]
  hover.text.all <- rep(NA, nrow(df))

  for (i in seq_len(nrow(df))) {
    sample <- rownames(df)[i]
    hover.text <- ""
    values <- df[sample, , drop = FALSE] %>% as("data.frame")
    for (variable in hover.variables) {
      value <- values[[variable]]
      hover.text <- paste0(hover.text, variable, ": ", value, "<br>")
    }
    hover.text.all[i] <- hover.text
  }

  return(hover.text.all)
}

#' Plot Alpha-Diversity
#' @export
plot_alpha_diversity <- function(
  full.sample.data, # containing alpha diversity, measure and depth
  hover.variables = NULL,
  x = NULL,
  x.levels = NULL,
  group = NULL,
  group.levels = NULL,
  facet.mode = "wrap", # "grid", "wrap"
  facet = NULL,
  facet.levels = NULL,
  facet.row = NULL,
  facet.row.levels = NULL,
  facet.col = NULL,
  facet.col.levels = NULL,
  facet_labeller = "label_value", #"label_both" or "label_value"
  shape = NULL,
  shape.levels = NULL,
  size = NULL,
  point.size = 1.5,
  remove.na.from.plot = FALSE,
  plot_type = NULL, # automatic if NULL
  color_vector = c("brown", "darkgreen", "orange", "violet"),
  stat = FALSE,
  check_depth = FALSE, # NOTE: not compatible with plotly
  alpha = 1
) {
  # Assuming alpha diversity is in the first column !
  measure <- colnames(full.sample.data)[1] # contains measure name

  if (is.null(facet.mode)) {
    facet.mode <- "wrap"
  }
  # Set proper data types

  if (is_valid_factor(full.sample.data, group)) {
    full.sample.data <- apply_levels(full.sample.data, group, group.levels)
  }

  if (is_valid_factor(full.sample.data, shape)) {
    full.sample.data <- apply_levels(full.sample.data, shape, shape.levels)
  } else {
    shape <- NULL
  }

  if (is_valid_factor(full.sample.data, x)) {
    full.sample.data <- apply_levels(full.sample.data, x, x.levels)
  }

  is_valid_facet <- is_valid_factor(full.sample.data, facet) &&
    facet.mode == "wrap"

  if (is_valid_facet) {
    full.sample.data <- apply_levels(full.sample.data, facet, facet.levels)
  }

  is_valid_facet_row <- is_valid_factor(full.sample.data, facet.row) &&
    facet.mode == "grid"

  if (is_valid_facet_row) {
    full.sample.data <- apply_levels(
      full.sample.data,
      facet.row,
      facet.row.levels
    )
  }

  is_valid_facet_col <- is_valid_factor(full.sample.data, facet.col) &&
    facet.mode == "grid"

  if (is_valid_facet_col) {
    full.sample.data <- apply_levels(
      full.sample.data,
      facet.col,
      facet.col.levels
    )
  }

  # Check whether depth is present
  check_depth <- check_depth && !is.null(full.sample.data$depth)

  # Remove all NA's from plot data (labels, facets or shape) by removing samples
  # having NA for at least one of graphical parameters (x-axis variable and group variable)
  # remove these samples from sample data
  if (remove.na.from.plot) {
    full.sample.data <- remove_nas(
      full.sample.data,
      c(x, group, facet, facet.row, facet.col, shape, size)
    )
  }

  # Prepare hover information about samples based on sample_data
  hover.variables <- c(
    hover.variables,
    if (check_depth) {
      "depth"
    }
  )
  full.sample.data$hover.text <-
    get_hover_text(full.sample.data, hover.variables)

  # Handle various plot types

  is_continuous_x <- !is.null(x) && is.numeric(full.sample.data[[x]])

  if (is.null(plot_type)) {
    if (is_continuous_x) {
      plot_type <- "scatter"
    } else {
      plot_type <- "boxplot"
    }
  }

  if (is.null(x)) {
    mapping <- aes(x = 0, y = .data[[measure]])
  } else {
    mapping <- aes(x = .data[[x]], y = .data[[measure]])
  }

  point_mapping <- aes()

  if (!is.null(group)) {
    point_mapping <- modifyList(
      point_mapping,
      aes(color = .data[[group]])
    )
  }

  if (length(hover.variables) > 0) {
    point_mapping <- modifyList(
      point_mapping,
      aes(text = .data[["hover.text"]])
    )
  }

  if (!is.null(shape)) {
    point_mapping <- modifyList(
      point_mapping,
      aes(shape = .data[[shape]])
    )
  }

  has_size_aes <- FALSE
  if (check_depth) {
    point_mapping <- modifyList(
      point_mapping,
      aes(size = .data[["depth"]])
    )
    has_size_aes <- TRUE
  } else if (!is.null(size)) {
    point_mapping <- modifyList(
      point_mapping,
      aes(size = .data[[size]])
    )
    has_size_aes <- TRUE
  }

  if (plot_type == "boxplot") {
    plot_layer_fn <- ggplot2::geom_boxplot
  }

  if (plot_type == "violin") {
    plot_layer_fn <- ggplot2::geom_violin
  }

  if (plot_type == "scatter") {
    plot_layer_fn <- if (has_size_aes) {
      function() {
        list(
          ggplot2::geom_point(point_mapping),
          ggplot2::scale_size(range = c(point.size * 0.5, point.size * 3))
        )
      }
    } else {
      function() {
        ggplot2::geom_point(point_mapping, size = point.size)
      }
    }
  }

  if (plot_type %in% c("boxplot", "violin")) {
    plt <- ggplot(full.sample.data, mapping) +
      plot_layer_fn(outlier.shape = NA)
  } else {
    plt <- ggplot(full.sample.data, mapping) +
      plot_layer_fn()
  }

  if (plot_type != "scatter") {
    if (has_size_aes) {
      plt <- plt +
        ggplot2::geom_jitter(point_mapping, alpha = alpha) +
        ggplot2::scale_size(range = c(point.size * 0.5, point.size * 3))
    } else {
      plt <- plt +
        ggplot2::geom_jitter(point_mapping, alpha = alpha, size = point.size)
    }
  } else {
    plt <- plt +
      stat_smooth(aes_string(color = group), method = "lm", alpha = 0.1)
  }

  if (stat) {
    if (plot_type != "scatter") {
      plt <- plt +
        stat_compare_means()
    } else {
      plt <- plt +
        ggpubr::stat_cor()
    }
  }

  plt <- plt + theme_minimal()

  if (is_valid_facet_row && is_valid_facet_col) {
    plt <- plt +
      facet_grid(
        rows = vars(!!sym(facet.row)),
        cols = vars(!!sym(facet.col)),
        labeller = facet_labeller
      )
    #facet_grid(.data[[facet.row]] ~ .data[[facet.col]])
  } else if (is_valid_facet_row) {
    plt <- plt +
      facet_grid(
        rows = vars(!!sym(facet.row)),
        labeller = facet_labeller
      )
  } else if (is_valid_facet_col) {
    plt <- plt +
      facet_grid(
        cols = vars(!!sym(facet.col)),
        labeller = facet_labeller
      )
  } else if (is_valid_facet) {
    plt <- plt +
      facet_wrap(
        ~ .data[[facet]],
        ncol = smart_facet_ncol(
          nlevels(factor(full.sample.data[[facet]]))
        ),
        labeller = facet_labeller
      )
  }

  # TODO: ignored it since ggarrange is incompatible with plotly after
  # else{
  #   if(check_depth){
  #
  #     plt.depth = ggplot(full.sample.data, aes_string(y=measure, x = "depth"))+
  #       geom_point(aes_string(color=group, shape=shape), size = size, alpha = alpha)+
  #       stat_smooth(aes_string(color= group), method = "lm", alpha=0.1 )+
  #       theme_minimal()+
  #       labs(x= "depth", y = measure)
  #     if(stat){
  #       plt.depth = plt.depth+
  #         ggpubr::stat_cor()
  #     }
  #     plt = ggarrange(plt, plt.depth, common.legend = TRUE)
  #   }
  # }

  plt <- plt +
    labs(title = "Alpha Diversity") +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold", hjust = 0.5)
    )
  if (!is.null(group) && !is.numeric(full.sample.data[[group]])) {
    plt <- plt +
      ggsci::scale_fill_npg()
  }

  if (is.null(x)) {
    plt <- plt +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }

  return(plt)
}


#' Get and Plot Alpha-Diversity from Phyloseq Object
#' @export
full_alpha_diversity <- function(
  physeq,
  taxrank = NULL,
  fraction_id_name = NULL,
  fraction_ids = NULL,
  measure = "Shannon",
  x,
  group = NULL,
  plot_type = "boxplot",
  hover.variables = hover.variables,
  color_vector = c("brown", "darkgreen", "orange", "violet"),
  stat = FALSE,
  check_depth = FALSE,
  size = 1.5,
  alpha = 1
) {
  full.sample.data <- get_alpha_diversity(
    physeq = physeq,
    taxrank = taxrank,
    fraction_id_name = fraction_id_name,
    fraction_ids = fraction_ids,
    measure = measure
  )
  p <- plot_alpha_diversity(
    full.sample.data = full.sample.data,
    hover.variables = hover.variables,
    x = x,
    group = group,
    plot_type = plot_type,
    color_vector = color_vector,
    stat = stat,
    check_depth = check_depth,
    size = size,
    alpha = alpha
  )

  return(p)
}

#' Get Richness Based on Ig Score Significance
#' @export
get_igseq_richness <-
  function(ps_ig_score, proportions = FALSE, low_lim = -5, high_lim = 5) {
    counts <- matrix(
      NA,
      nrow = nsamples(ps_ig_score),
      ncol = 3,
      dimnames = list(sample_names(ps_ig_score), c("down", "ns", "up"))
    )
    for (sample_id in sample_names(ps_ig_score)) {
      counts[sample_id, "down"] <- sum(
        otu_table(ps_ig_score)[sample_id, ] < low_lim,
        na.rm = TRUE
      )
      counts[sample_id, "up"] <- sum(
        otu_table(ps_ig_score)[sample_id, ] > high_lim,
        na.rm = TRUE
      )
      counts[sample_id, "ns"] <- sum(
        otu_table(ps_ig_score)[sample_id, ] >= low_lim &
          otu_table(ps_ig_score)[sample_id, ] <= high_lim,
        na.rm = TRUE
      )
    }

    if (proportions) {
      counts[, c("down", "ns", "up")] <- t(apply(
        counts[, c("down", "ns", "up")],
        1,
        function(row) {
          row / sum(row, na.rm = TRUE)
        }
      ))
    }

    igseq_richness_df <- cbind(counts, sample_data(ps_ig_score)) %>%
      as.data.frame() %>%
      mutate(sample_id = rownames(counts), .before = 1) %>%
      pivot_longer(
        cols = c("down", "ns", "up"),
        values_to = "richness",
        names_to = "significance"
      ) %>%
      as.data.frame()

    #rownames(igseq_richness_df) = igseq_richness_df$sample_id

    igseq_richness_df <- igseq_richness_df[, unique(c(
      "richness",
      colnames(igseq_richness_df)
    ))]

    return(igseq_richness_df)
  }

#' Plot Richness Based on Ig Score Significance
#' @export
plot_igseq_richness <-
  function(igseq_richness_df, group, color, exclude_ns = FALSE) {
    if (exclude_ns) {
      igseq_richness_df <- igseq_richness_df[
        igseq_richness_df$significance != "ns",
      ]
    }

    igseq_richness_df[[group]] <- factor(
      igseq_richness_df[[group]],
      levels = gtools::mixedsort(unique(igseq_richness_df[[group]]))
    )

    plt <-
      ggplot(
        igseq_richness_df,
        aes(x = .data[[group]], y = richness, color = .data[[color]])
      ) +
      geom_jitter() +
      geom_violin(alpha = 0.1) +
      facet_grid(~significance) +
      theme_minimal()
    return(plt)
  }

# get_... - gets the fit object and coordinates, parameters
# stat_... - compares groups
# plot_... - plots the result with plot() function as a sideeffect
# ggplot_... - outputs a ggplot instead
# full_... - imitates the original function

# constrained_beta_dispersion() is very similar, so it's in this function as well
# BUT: original functions for constrained and unconstrained beta-dispersion
# don't return the same object, so this one returns the format of original
