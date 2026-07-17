# Implements the idea of MA plot (from microarray analysis) to abundance analysis

#' Compute MA-Plot Coordinates Between Ig-Seq Fractions
#'
#' Computes MA-plot-style coordinates (borrowed from microarray analysis: `M`
#' = log-ratio, `A` = log-abundance) between a positive Ig-Seq fraction and a
#' negative fraction, per taxon. If a second negative fraction is supplied,
#' the same coordinates are also computed between the two negative
#' fractions, to serve as an empirical null distribution (see
#' [get_slide_z()]).
#'
#' @param sorted_sample_df A data frame for one biological sample, as
#'   produced by [group_sorted_samples()], with one row per taxon and one
#'   column per fraction.
#' @param positive_fraction_name Name of the column in `sorted_sample_df`
#'   holding the positive (Ig+) fraction abundance.
#' @param first_negative_fraction_name Name of the column holding the
#'   primary negative fraction abundance (e.g. 9/10 of the whole negative
#'   fraction for IgSeq).
#' @param second_negative_fraction_name Name of the column holding a second
#'   negative fraction abundance (e.g. the remaining 1/10), used to build an
#'   empirical null distribution, or `NULL` (default) to skip it.
#'
#' @return A data frame with one row per taxon: `taxon_id`, `sample_id`, the
#'   original fraction abundances (`pos`, `neg1`, `neg2`; `neg2` is `NA` when
#'   `second_negative_fraction_name` is `NULL`), the observed coordinates
#'   (`obs_abundance`, `obs_change`, positive vs. first negative fraction),
#'   and, when a second negative fraction is supplied, the empirical-null
#'   coordinates (`null_abundance`, `null_change`, first vs. second negative
#'   fraction); otherwise `NA`.
#'
#' @examples
#' data(ps_igseq)
#' grouped <- group_sorted_samples(
#'   physeq = ps_igseq,
#'   sample_id_name = "sample_id",
#'   sample_ids = c("sample_1", "sample_2"),
#'   fraction_id_name = "sorting_fraction",
#'   fraction_ids = c("Pos", "Neg1", "Neg2")
#' )
#' ma_coords <- get_ma_coordinates(
#'   sorted_sample_df = grouped[["sample_1"]],
#'   positive_fraction_name = "Pos",
#'   first_negative_fraction_name = "Neg1",
#'   second_negative_fraction_name = "Neg2"
#' )
#' head(ma_coords)
#'
#' @export
get_ma_coordinates <- function(
  sorted_sample_df, # dataframe from the result of group_sorted_samples()
  positive_fraction_name,
  first_negative_fraction_name, # 9/10 of the whole negative fraction for IgSeq
  second_negative_fraction_name = NULL # 1/10
) {
  # Retrieve taxa abundances for each fraction
  pos <- sorted_sample_df[, positive_fraction_name]
  neg1 <- sorted_sample_df[, first_negative_fraction_name]
  if (!is.null(second_negative_fraction_name)) {
    neg2 <- sorted_sample_df[, second_negative_fraction_name]
    empirical_null <- TRUE
  } else {
    neg2 <- rep(NA, nrow(sorted_sample_df))
    empirical_null <- FALSE
  }

  # Compute MA coordinates (M = log-ratio, A = log-abundance), see MA plot.
  # Alternatives to log2/log10 (e.g. asin(sqrt(x/sum(x))), sqrt(x/sum(x)))
  # are worth revisiting if 0-handling becomes an issue, since they don't
  # blow up at x = 0 the way a log transform does.
  transform_M <- function(x) log2(x)
  transform_A <- function(x) log10(x)

  # Obs - observed, null - empirical null distribution (control vs control)
  compute_ma_pair <- function(a, b) {
    abundance <- transform_A(a) + transform_A(b)
    change <- transform_M(a) - transform_M(b)
    abundance[is.nan(abundance) | is.infinite(abundance)] <- NA
    change[is.nan(change) | is.infinite(change)] <- NA
    list(abundance = abundance, change = change)
  }

  obs <- compute_ma_pair(pos, neg1)
  obs_abundance <- obs$abundance
  obs_change <- obs$change

  if (empirical_null) {
    null <- compute_ma_pair(neg1, neg2)
    null_abundance <- null$abundance
    null_change <- null$change
  } else {
    null_abundance <- rep(NA, nrow(sorted_sample_df))
    null_change <- rep(NA, nrow(sorted_sample_df))
  }

  ma_coords <-
    data.frame(
      taxon_id = sorted_sample_df$taxon_id,
      sample_id = sorted_sample_df$sample_id,
      pos,
      neg1,
      neg2,
      obs_abundance,
      obs_change,
      null_abundance,
      null_change
    )

  return(ma_coords)
}

#' Build MA-Plot Data Across Zero-Handling Strategies
#'
#' Wraps [get_ma_coordinates()] to build a plotting-ready long-format MA-plot
#' data frame for one biological sample, across one or more zero-count
#' handling strategies (see [impute_zeros()]). Comparing several
#' `zero_treatments` side by side (e.g. via [plot_ma()]) shows how sensitive
#' the resulting MA-plot is to the choice of zero-handling method.
#'
#' @param sorted_sample_df A data frame for one biological sample, as
#'   produced by [group_sorted_samples()], with one row per taxon and one
#'   column per fraction.
#' @param positive_fraction_name Name of the column in `sorted_sample_df`
#'   holding the positive (Ig+) fraction abundance.
#' @param first_negative_fraction_name Name of the column holding the
#'   primary negative fraction abundance (e.g. 9/10 of the whole negative
#'   fraction for IgSeq).
#' @param second_negative_fraction_name Name of the column holding a second
#'   negative fraction abundance (e.g. the remaining 1/10), used to build an
#'   empirical null distribution, or `NULL` (default) to skip it.
#' @param zero_treatments Character vector of [impute_zeros()] `method`
#'   values to compute MA coordinates for; each contributes one
#'   `zero_treatment` level to the returned `plot_data`. Default
#'   `"keep_zeros"` (zeros left as-is, no imputation).
#'
#' @return A list with:
#'   \describe{
#'     \item{`sample_id`}{The biological sample id
#'       (`sorted_sample_df$sample_id[1]`).}
#'     \item{`nb_zero_taxa`}{Number of taxa with a zero abundance in at
#'       least one compared fraction, before imputation.}
#'     \item{`plot_data`}{A long-format data frame (`M`, `A`, `comparison`,
#'       `taxon_id`, `zero_treatment`) ready for [plot_ma()].}
#'     \item{`imputed_taxa`}{The union, across all `zero_treatments`, of
#'       `taxon_id`s imputed under that treatment (see [impute_zeros()]).}
#'   }
#'
#' @examples
#' data(ps_igseq)
#' grouped <- group_sorted_samples(
#'   physeq = ps_igseq,
#'   sample_id_name = "sample_id",
#'   sample_ids = c("sample_1", "sample_2"),
#'   fraction_id_name = "sorting_fraction",
#'   fraction_ids = c("Pos", "Neg1", "Neg2")
#' )
#' ma_plot_data <- get_ma_plot_data(
#'   sorted_sample_df = grouped[["sample_1"]],
#'   positive_fraction_name = "Pos",
#'   first_negative_fraction_name = "Neg1",
#'   second_negative_fraction_name = "Neg2",
#'   zero_treatments = c("keep_zeros", "pseudo_count")
#' )
#' head(ma_plot_data$plot_data)
#'
#' @export
get_ma_plot_data <- function(
  sorted_sample_df, # dataframe from the result of group_sorted_samples()
  positive_fraction_name,
  first_negative_fraction_name,
  second_negative_fraction_name = NULL,
  zero_treatments = c("keep_zeros")
) {
  pos <- sorted_sample_df[, positive_fraction_name]
  neg1 <- sorted_sample_df[, first_negative_fraction_name]

  if (!is.null(second_negative_fraction_name)) {
    neg2 <- sorted_sample_df[, second_negative_fraction_name]
    empirical_null <- TRUE
  } else {
    empirical_null <- FALSE
  }

  nb_zero_taxa <-
    if (empirical_null) {
      sum(pos == 0 | neg1 == 0 | neg2 == 0)
    } else {
      sum(pos == 0 | neg1 == 0)
    }

  sample_id <- sorted_sample_df$sample_id[1]

  plot_data <- data.frame()
  imputed_taxa <- NULL

  for (zero_treatment in zero_treatments) {
    zero_imputation_result <-
      impute_zeros(
        data = sorted_sample_df,
        # Don't impute zeros in other fractions!
        fraction_names = c(
          positive_fraction_name,
          first_negative_fraction_name,
          second_negative_fraction_name
        ),
        method = zero_treatment
      )

    sorted_sample_df_imputed <- zero_imputation_result$data
    # Union across treatments so a taxon flagged as imputed under any one
    # of them is still highlighted in its facet, regardless of the order
    # `zero_treatments` was given in.
    imputed_taxa <- union(imputed_taxa, zero_imputation_result$imputed_taxa)

    ma_coords <- get_ma_coordinates(
      sorted_sample_df = sorted_sample_df_imputed,
      positive_fraction_name = positive_fraction_name,
      first_negative_fraction_name = first_negative_fraction_name,
      second_negative_fraction_name = second_negative_fraction_name
    )

    # Convert to longer format
    # coordinates for pos vs neg1
    ma_coords_long <-
      data.frame(
        M = ma_coords$obs_change,
        A = ma_coords$obs_abundance,
        comparison = paste0(
          positive_fraction_name,
          " vs ",
          first_negative_fraction_name
        ),
        taxon_id = ma_coords$taxon_id
      )

    # coordinates for neg1 vs neg2
    if (empirical_null) {
      ma_coords_long <- rbind(
        ma_coords_long,

        data.frame(
          M = ma_coords$null_change,
          A = ma_coords$null_abundance,
          comparison = paste0(
            first_negative_fraction_name,
            " vs ",
            second_negative_fraction_name
          ),
          taxon_id = ma_coords$taxon_id
        )
      )
    }

    plot_data <- rbind(
      plot_data,
      cbind(
        ma_coords_long,
        data.frame(
          zero_treatment = rep(
            gsub("_", " ", zero_treatment),
            nrow(ma_coords_long)
          )
        )
      )
    )
  }

  plot_data$zero_treatment <- factor(
    plot_data$zero_treatment,
    levels = gsub("_", " ", zero_treatments)
  )

  return(list(
    sample_id = sample_id,
    nb_zero_taxa = nb_zero_taxa,
    plot_data = plot_data,
    imputed_taxa = imputed_taxa
  ))
}


#' Plot an Ig-Seq MA Plot
#'
#' Draws the MA-plot coordinates produced by [get_ma_plot_data()] for one
#' biological sample: log-abundance (`A`) on the x-axis, log-ratio (`M`) on
#' the y-axis. Taxa with a zero abundance in at least one compared fraction
#' (`ma_plot_data$imputed_taxa`) are drawn separately, jittered just past the
#' plot's x-range rather than placed at their imputed/undefined `A` value.
#'
#' @param ma_plot_data A list as returned by [get_ma_plot_data()].
#' @param ellipses If `TRUE`, overlay a 95% normal-theory confidence ellipse
#'   per `comparison` (only applies when `type = "facet"`). Default `FALSE`.
#' @param type One of `"facet"` (default; facet by `zero_treatment`, color
#'   by `comparison`) or `"superposed"` (facet by `comparison`, color by
#'   `zero_treatment`, comparing zero-handling strategies within the same
#'   panel).
#'
#' @return A `ggplot` object.
#'
#' @examples
#' data(ps_igseq)
#' grouped <- group_sorted_samples(
#'   physeq = ps_igseq,
#'   sample_id_name = "sample_id",
#'   sample_ids = c("sample_1", "sample_2"),
#'   fraction_id_name = "sorting_fraction",
#'   fraction_ids = c("Pos", "Neg1", "Neg2")
#' )
#' ma_plot_data <- get_ma_plot_data(
#'   sorted_sample_df = grouped[["sample_1"]],
#'   positive_fraction_name = "Pos",
#'   first_negative_fraction_name = "Neg1",
#'   second_negative_fraction_name = "Neg2",
#'   zero_treatments = c("keep_zeros", "pseudo_count")
#' )
#' plot_ma(ma_plot_data)
#' plot_ma(ma_plot_data, type = "superposed")
#'
#' @export
plot_ma <-
  function(ma_plot_data, ellipses = FALSE, type = c("facet", "superposed")) {
    type <- match.arg(type)

    ma_non_imputed <- ma_plot_data$plot_data[
      !ma_plot_data$plot_data$taxon_id %in% ma_plot_data$imputed_taxa,
    ]

    ma_imputed <- ma_plot_data$plot_data[
      ma_plot_data$plot_data$taxon_id %in% ma_plot_data$imputed_taxa,
    ]

    jitter <- .jitter_offset(ma_non_imputed$A)
    jitter_width <- jitter$width
    jitter_x <- jitter$x

    plt <- ggplot(ma_non_imputed, aes(x = A, y = M))

    if (type == "facet") {
      plt <- plt +
        geom_point(aes(color = comparison), alpha = 0.8) +
        geom_jitter(
          data = ma_imputed,
          aes(x = jitter_x, y = M, color = comparison),
          alpha = 0.8,
          position = position_jitter(width = jitter_width, height = 0)
        ) +
        facet_wrap(. ~ zero_treatment, scales = "fixed")

      if (ellipses) {
        plt <- plt +
          stat_ellipse(
            aes(color = comparison),
            type = "norm",
            level = 0.95,
            alpha = 0.7
          )
      }
    } else if (type == "superposed") {
      plt <- plt +
        geom_point(
          aes(color = zero_treatment, size = zero_treatment),
          alpha = 0.5
        ) +
        geom_jitter(
          data = ma_imputed,
          aes(
            x = jitter_x,
            y = M,
            color = zero_treatment,
            size = zero_treatment
          ),
          alpha = 0.5,
          position = position_jitter(width = jitter_width, height = 0)
        ) +
        facet_wrap(. ~ comparison, scales = "fixed")
    }

    plt <- plt +
      labs(
        title = paste0("Ig-Seq MA plot for sample ", ma_plot_data$sample_id),
        subtitle = paste0(
          ma_plot_data$nb_zero_taxa,
          " taxa with zero abundances in at least one fraction"
        ),
        x = latex2exp::TeX(
          "Log-Abundance: $\\log_{10}\\left(\\fraction_{1}\\cdot\\fraction_{2}\\right)$"
        ),
        y = latex2exp::TeX(
          "Log-Ratio: $\\log_{2}\\left(\\frac{\\fraction_{1}}{\\fraction_{2}}\\right)$"
        )
      ) +
      theme_minimal() +
      .plot_title_theme() +
      scale_color_manual(values = c("darkgray", ggsci::pal_npg()(4)))
    #ggsci::scale_color_npg()

    return(plt)
  }
