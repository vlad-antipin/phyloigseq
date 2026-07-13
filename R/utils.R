#' PhyloIgSeq
#'
#' Package containing wrappers mainly around \pkg{phyloseq} and \pkg{vegan} packages
#' for downstream analysis of 16S amplicon sequencing results and allowing to further
#' analyze the Ig-Seq data, computing and plotting various Ig scores. The functions
#' from this package are used in PhyloIgSeq Shiny Web application on
#'  \url{https://www.immulab.fr/cms/index.php/team/tools/lab-tools}.
"_PACKAGE"

# Compute the optimal ncol for facet_wrap: minimize empty cells, prefer wider
# (more columns) layouts. Avoids degenerate single-row/column results for n > 3.
smart_facet_ncol <- function(n) {
  if (n <= 3) {
    return(n)
  }
  divisors <- which(n %% seq_len(n) == 0)
  perfect_wide <- divisors[divisors >= sqrt(n) & divisors <= ceiling(n / 2)]
  if (length(perfect_wide) > 0) {
    return(min(perfect_wide))
  }
  nc_range <- seq(ceiling(sqrt(n)), ceiling(n / 2))
  empties <- ceiling(n / nc_range) * nc_range - n
  max(nc_range[empties == min(empties)])
}


#' Ig-Coating Score Names
#'
#' Character vector of the Ig-coating score names computable via \code{\link{compute_ig_score}}
#' and requestable through \code{\link{getPhyloIgSeq}}'s \code{scores} argument.
#'
#' @import openxlsx
#' @import DT
#' @import ade4
#' @import vegan
#' @import Rtsne
#' @import umap
#' @import VIM
#' @import dplyr
#' @import tidyr
#' @import ggpubr
#' @import tools
#' @import htmlwidgets
#' @import rstatix
#' @import scales
#' @import phyloseq
#' @import ggplot2
#' @import methods
#' @importFrom microbiome transform
#' @import zCompositions
#' @importFrom magrittr %>%
#' @importFrom gtools mixedsort
#' @import rlang
#' @import ggExtra
#' @importFrom car dataEllipse
#' @import sp
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation anno_text
#' @import ragg
#' @import gganimate
#' @import gifski
#' @importFrom utils globalVariables
#'
#' @examples
#' IG_SCORES
#'
#' @export
IG_SCORES <- c("slide_z", "palm", "kau", "prob_index", "prob_ratio")
# scores to come:
# "purity_corrected_prob_index", "purity_corrected_prob_ratio")

# Prevent R CMD check NOTE about undefined global variable
utils::globalVariables("IG_SCORES")

#' Similar to veganifyOTU from phyloseq.
#' @keywords internal
#' @export
reverseASV <- function(physeq) {
  if (taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  return(physeq)
}

# Fix vertical jitter issue:
geom_jitter <- function(..., height = 0) {
  # avoid conflicts with `position` argument if it was explicitly specified
  if (hasArg(position)) {
    ggplot2::geom_jitter(...)
  } else {
    ggplot2::geom_jitter(..., height = height)
  }
}


# Custom function, returns TRUE or FALSE depending on whether the number is in the
# the interval given is a vector, e.g. 3.5 %in_interval% c(3,4) gives TRUE
`%in_interval%` <- function(x, interval) {
  interval <- range(interval, na.rm = TRUE)
  x >= interval[1] & x <= interval[2]
}

transform_abundances <-
  function(
    abundance_table, # matrix
    transform = c("compositional", "hellinger"), # ,"clr", "log10", "log10p"),
    taxa_are_rows = TRUE
  ) {
    # assumes taxa are rows
    if (!taxa_are_rows) {
      abundance_table <- t(abundance_table)
    }

    if (transform == "compositional") {
      abundance_table <- sweep(
        abundance_table,
        2,
        colSums(abundance_table),
        "/"
      )
    } else if (transform == "hellinger") {
      abundance_table <- sqrt(sweep(
        abundance_table,
        2,
        colSums(abundance_table),
        "/"
      ))
    } else {
      # TODO:
      # else if(transform == "clr"){
      #
      # }else if(transform == "log10"){
      #
      # }else if(transform == "log10p"){
      #
      # }
      stop("Wrong transformation method")
    }

    if (!taxa_are_rows) {
      return(t(abundance_table))
    } else {
      return(abundance_table)
    }
  }

#' Plot Rarefaction Curves
#'
#' Computes and plots rarefaction curves for all samples in a \code{phyloseq}
#' object using \code{\link[vegan]{rarecurve}}.
#'
#' @param ps A \code{\link[phyloseq]{phyloseq}} object with an OTU table and
#'   sample data. Samples with zero total counts are automatically excluded.
#' @param step Integer. Sampling-depth increment passed to
#'   \code{\link[vegan]{rarecurve}}. Smaller values yield smoother curves at
#'   the cost of computation time. Default \code{100}.
#' @param show_legend Logical. Whether to display the sample legend.
#'   Default \code{TRUE}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object with one rarefaction curve
#'   per sample, coloured by sample name.
#'
#' @seealso \code{\link[vegan]{rarecurve}}, \code{\link{plot_seq_depth}}
#'
#' @export
plot_rarefaction <- function(ps, step = 100, show_legend = TRUE) {
  # Extract OTU table
  otu <- as(otu_table(ps), "matrix")

  # Ensure samples are rows
  if (taxa_are_rows(ps)) {
    otu <- t(otu)
  }

  # Remove samples with zero counts
  otu <- otu[rowSums(otu) > 0, ]

  # Compute rarefaction curves
  rare_df <- rarecurve(
    otu,
    step = step,
    xlab = "Reads",
    ylab = "Richness",
    tidy = TRUE,
    label = FALSE
  )
  rare_df <- rare_df %>%
    rename(Reads = Sample, Richness = Species, Sample = Site)

  # names(rare_list) <- row.names(ps@sam_data)
  # # Convert rarecurve output to dataframe
  # rare_df <- lapply(names(rare_list), function(sample) {
  #   data.frame(
  #     Sample = sample,
  #     Reads = attr(rare_list[[sample]], "Subsample"),
  #     Richness = rare_list[[sample]]
  #   )
  # }) %>%
  #   bind_rows()

  # print(head(rare_df))

  # Plot
  p <- ggplot(
    rare_df,
    aes(x = Reads, y = Richness, group = Sample, color = Sample)
  ) +
    geom_line(alpha = 0.8) +
    theme_bw() +
    labs(
      x = "Sequencing depth",
      y = "Observed richness",
      title = "Rarefaction Curves"
    )

  p <- p +
    theme_minimal() +
    ggplot2::theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.title = if (show_legend) {
        element_text(face = "bold", hjust = 0.5)
      },
      legend.position = if (show_legend) {
        "right"
      } else {
        "none"
      }
    )

  return(p)
}

#' Plot Sequencing Depth
#'
#' Visualises per-sample sequencing depth from a \code{\link[phyloseq]{phyloseq}}
#' object as either a bar chart (one bar per sample) or a grouped boxplot with
#' jittered points.
#'
#' @param ps A \code{\link[phyloseq]{phyloseq}} object with an OTU table and
#'   sample data.
#' @param type Character. Plot type: \code{"bar"} (default) for a bar chart
#'   ordered by sample name, or \code{"box"} for a boxplot grouped by a
#'   metadata variable.
#' @param x_var Character. Name of a sample-data column to use as the x-axis
#'   grouping variable. Required when \code{type = "box"}, ignored otherwise.
#' @param facet_var Character or \code{NULL}. Name of a sample-data column used
#'   to facet the plot with \code{\link[ggplot2]{facet_wrap}}. Only applied
#'   when \code{type = "box"}. Default \code{NULL} (no faceting).
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @seealso \code{\link{plot_rarefaction}}
#'
#' @export
plot_seq_depth <- function(
  ps,
  type = c("bar", "box"),
  x_var = NULL,
  facet_var = NULL
) {
  type <- match.arg(type)

  # Compute sequencing depth
  depth_df <- data.frame(
    Sample = sample_names(ps),
    Depth = sample_sums(ps)
  )

  # Add metadata
  meta <- as.data.frame(sample_data(ps))
  meta$Sample <- rownames(meta)

  depth_df <- left_join(depth_df, meta, by = "Sample")

  # -----------------------
  # BARPLOT
  # -----------------------
  if (type == "bar") {
    p <- ggplot(depth_df, aes(x = Sample, y = Depth)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      labs(x = "Sample", y = "Sequencing depth") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }

  # -----------------------
  # BOXPLOT + JITTER
  # -----------------------
  if (type == "box") {
    if (is.null(x_var)) {
      stop("For boxplot, please provide x_var (independent variable).")
    }

    p <- ggplot(depth_df, aes_string(x = x_var, y = "Depth")) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.7) +
      theme_bw() +
      labs(x = x_var, y = "Sequencing depth")

    if (!is.null(facet_var)) {
      p <- p + facet_wrap(as.formula(paste("~", facet_var)))
    }
  }

  return(p)
}

#' Rarefy Abundances to Same Depth by multinomial resampling
#' @export
rarefy_abundances <-
  function(
    abundance_table,
    taxa_are_rows = TRUE,
    common_count_sum = NULL,
    trim_taxa = TRUE,
    trim_samples = TRUE,
    silent_warnings = FALSE
  ) {
    if (taxa_are_rows) {
      sample_margin <- 2
      taxa_margin <- 1
      sampleSums <- colSums
      taxaSums <- rowSums
      nsamples <- ncol
      # assign_to_sample = function(sample, x) abundance_table[,sample ] <<- x
      # sample_slice = function(sample) abundance_table[,sample]
    } else {
      sample_margin <- 1
      taxa_margin <- 2
      sampleSums <- rowSums
      taxaSums <- colSums
      nsamples <- nrow
      # assign_to_sample = function(sample, x) abundance_table[sample, ] <<- x
      # sample_slice = function(sample) abundance_table[sample,]
    }

    if (is.null(common_count_sum)) {
      common_count_sum <- min(setdiff(sampleSums(abundance_table), 0))
    }

    dims_orig <- dim(abundance_table)

    if (trim_samples) {
      if (taxa_are_rows) {
        abundance_table <- abundance_table[,
          sampleSums(abundance_table) >= common_count_sum,
          drop = FALSE
        ]
      } else {
        abundance_table <- abundance_table[
          sampleSums(abundance_table) >= common_count_sum, ,
          drop = FALSE
        ]
      }
    }

    # For each sample, draw counts from multinomial distribution over taxa
    # for (sample in seq_len(nsamples(abundance_table))) {
    #     sample_sum = sum(sample_slice(sample))
    #     assign_to_sample( sample, if( sample_sum == 0){ 0 }else{
    #       rmultinom(1, size = common_count_sum, prob = sample_slice(sample) / sample_sum)[, 1]
    #       })
    #
    # }

    abundance_table <-
      apply(abundance_table, sample_margin, function(sample_counts) {
        if (sum(sample_counts) == 0) {
          return(rep(0, length(sample_counts)))
        }
        rmultinom(
          1,
          size = common_count_sum,
          prob = sample_counts / sum(sample_counts)
        )[, 1]
      })

    if (sample_margin == 1) {
      # apply() to rows puts them in columns, so have to transoose the matrix back

      abundance_table <- t(abundance_table)
    }

    if (trim_taxa) {
      if (taxa_are_rows) {
        abundance_table <- abundance_table[
          taxaSums(abundance_table) > 0, ,
          drop = FALSE
        ]
      } else {
        abundance_table <- abundance_table[,
          taxaSums(abundance_table) > 0,
          drop = FALSE
        ]
      }
    }

    dims_rare <- dim(abundance_table)

    if (!all(dims_rare == dims_orig) & !silent_warnings) {
      warning(paste0(
        "Trimmed ",
        dims_orig[sample_margin] - dims_rare[sample_margin],
        " samples (with count sum < common_count_sum) and ",
        dims_orig[taxa_margin] - dims_rare[taxa_margin],
        " taxa (zero counts after rarefaction)",
        "\n"
      ))
    }

    return(abundance_table)
  }

#' Check whether abundances look like counts
#' @export
is.count.like <- function(x, allow.na = TRUE, consider.small.part = TRUE) {
  if (!is.numeric(x)) {
    return(FALSE)
  }

  if (is(x, "sparse_otu_table")) {
    # Only check stored (non-zero) values; zero is trivially count-like.
    vals <- x@sparse_data@x
    if (consider.small.part) vals <- head(vals, 10000L)
  } else {
    # Covers both otu_table and any subclass that is not sparse_otu_table.
    if (is(x, "otu_table")) {
      x <- as(x, "matrix")
    }
    if (consider.small.part) {
      x <- x[seq_len(min(100L, nrow(x))), seq_len(min(100L, ncol(x)))]
    }
    vals <- c(x)
  }

  if (allow.na) {
    vals <- vals[!is.na(vals)]
  } else if (anyNA(vals)) {
    return(FALSE)
  }

  all(vals >= 0 & is.finite(vals) & vals == floor(vals))
}

#' Make taxonomy names unique at all taxonomic rank levels
#' @export
make_unique_taxa_table <- function(taxa_table) {
  taxa_table[is.na(taxa_table)] <- "NA"

  # "good duplication" example:
  # low level duplications correspond to whole taxonomy duplications
  #  ASV   | ... |    Genus     | Species
  # __________________________________
  #  ASV1  | ... | Bacteroides  | caccae
  #  ASV2  | ... | Bacteroides  | caccae

  # "bad duplication" example:
  #  ASV   | ... |    Genus     | Species
  # ___________________________________
  #  ASV1  | ... | Bacteroides  | caccae
  #  ASV2  | ... | Anaerostipes | caccae

  # Will be changed to :
  #  ASV   | ... |    Genus     | Species
  # ___________________________________
  #  ASV1  | ... | Bacteroides  | caccae
  #  ASV2  | ... | Anaerostipes | caccae.1

  # Go down taxonomy hierarchy (assuming it goes from left to right in the table)
  # the topmost level of taxonomy has only "good" duplications
  for (i in 2:ncol(taxa_table)) {
    # taxonomy if defined by the current level + all higher levels
    # (here, by just one higher level since it's ensured to be unique)
    taxonomies <- apply(taxa_table[, (i - 1):i], 1, function(row) {
      paste(row, collapse = "#")
    })
    # get rid of "good duplications"
    unique_taxonomies <- unique(taxonomies)

    # identify "bad duplications" and make them unique
    duplicated_names <- sapply(
      strsplit(unique_taxonomies, "#"),
      function(taxons) {
        tail(taxons, 1)
      }
    )
    unique_names <- make.unique(duplicated_names)

    # skip if no "bad duplications" - most often the case
    if (all(unique_names == duplicated_names)) {
      next
    }

    # replace "bad duplications" (vectorized: a loop doing
    # `taxonomies == unique_taxonomies[j]` per name is O(n * length(unique_names))
    # and can take a minute-plus on large tables)
    taxonomies <- unique_names[match(taxonomies, unique_taxonomies)]

    taxa_table[, i] <- taxonomies
  }

  return(taxa_table)
}

#' Impute with central tendency
#' @export
impute_with_central_tendency <- function(df, central.tendency = "median") {
  df <- as.data.frame(df)
  # Function to calculate the mode
  get.mode <- function(v) {
    v <- v[!is.na(v)]
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }

  # Loop through each column in the data frame
  for (col in colnames(df)) {
    # Impute numeric columns
    if (is.numeric(df[[col]])) {
      if (central.tendency == "mean") {
        df[[col]][is.na(df[[col]])] <- mean(df[[col]], na.rm = TRUE)
      } else if (central.tendency == "median") {
        df[[col]][is.na(df[[col]])] <- median(df[[col]], na.rm = TRUE)
      } else if (central.tendency == "mode") {
        df[[col]][is.na(df[[col]])] <- get.mode(df[[col]])
      } else {
        stop("Invalid method. Choose 'mean', 'median', or 'mode'.")
      }
      # Impute categorical columns, character columns are not affected /!\
    } else if (is.factor(df[[col]])) {
      df[[col]][is.na(df[[col]])] <- get.mode(df[[col]])
    }
  }
  return(df)
}


#' Impute Data with KNN, Central Tendency or Zero
#' @export
dataImpute <- function(
  data.tmp,
  exceptions = NULL,
  method = "KNN", # "KNN" or "Central Tendency"
  central.tendency = "median", # "mean", "median" or "mode"
  nb.neighbors = 5, # /!\ Find optimal usually it's between 3 and 10
  add.imputation.indicators = FALSE
) {
  if (method == "KNN") {
    data.imputed <-
      VIM::kNN(
        data.tmp[, !colnames(data.tmp) %in% exceptions],
        k = nb.neighbors,
        imp_var = add.imputation.indicators
      )
  } else if (method == "Central Tendency") {
    data.imputed <-
      impute_with_central_tendency(
        data.tmp[, !colnames(data.tmp) %in% exceptions],
        central.tendency = central.tendency
      )
  } else if (method == "Replace NA with 0") {
    data.imputed <- data.tmp
    data.imputed[is.na(data.imputed)] <- 0
  }

  data.imputed <- cbind(
    data.tmp[, colnames(data.tmp) %in% exceptions, drop = FALSE],
    data.imputed
  )

  return(data.imputed)
}

#' Plot Phylogenetic Tree from Phyloseq Object
#' @export
plot_phylo_tree <- function(
  physeq,
  taxrank = NULL,
  fraction_id_name = NULL,
  fraction_ids = NULL,
  layout = "rectangular",
  tip_color = NULL,
  label.tips = NULL,
  label_size = 2.5,
  ladderize = "left",
  ...
) {
  if (is.null(physeq) || class(physeq) != "phyloseq") {
    stop("Need a phyloseq object")
  }
  if (is.null(access(physeq, "phy_tree"))) {
    stop("Phyloseq object has to contain a tree")
  }

  if (!is.null(fraction_id_name) && !is.null(fraction_ids)) {
    physeq <- prune_samples(
      sample_data(physeq)[[fraction_id_name]] %in% fraction_ids,
      physeq
    )
  }

  if (!is.null(taxrank)) {
    physeq <- tax_glom(physeq = physeq, taxrank = taxrank)
    taxa_names(physeq) <- make.unique(tax_table(physeq)[, taxrank])
  }

  # Truncate tip labels so long ASV hashes don't clutter the plot;
  # make.unique() prevents duplicate names after truncation
  if (!is.null(label.tips)) {
    taxa_names(physeq) <- make.unique(substr(taxa_names(physeq), 1, 25))
  }

  tree_obj <- phy_tree(physeq)
  if (!is.null(ladderize)) {
    tree_obj <- ape::ladderize(tree_obj, right = (ladderize != "left"))
  }

  p <- ggtree::ggtree(tree_obj, layout = layout)

  if (!is.null(label.tips)) {
    p <- p + ggtree::geom_tiplab(size = label_size, align = FALSE)
  }

  if (!is.null(tip_color) && tip_color %in% colnames(tax_table(physeq))) {
    tax_df <- as.data.frame(tax_table(physeq))
    tax_df$.join_label <- taxa_names(physeq)
    # %<+% can fail when p$data has a duplicate "label" column in some ggtree
    # versions; merge directly into p$data to avoid the conflict
    p$data <- dplyr::left_join(p$data, tax_df, by = c("label" = ".join_label"))
    p <- p +
      ggtree::geom_tippoint(
        ggplot2::aes(color = .data[[tip_color]]),
        size = 2
      ) +
      ggsci::scale_color_igv(na.value = "grey50", name = tip_color)
  }

  p
}
