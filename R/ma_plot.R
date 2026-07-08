# analogous to ImmuMicrobiome::log_ratio()

#' Get MA Coordinates
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

  # Compute MA coordinates
  # see MA plot

  transform_M <- function(x) {
    log2(x) # DEFAULT
    #asin(sqrt(x/sum(x)))
    #sqrt(x/sum(x))
    #x
  }
  transform_A <- function(x) {
    log10(x) # DEFAULT
    #asin(sqrt(x/sum(x)))
    #sqrt(x/sum(x))
    #x
  }
  # Obs - observed, null - empirical null distribution (control vs control)
  obs_abundance <- transform_A(pos) + transform_A(neg1) #log10(pos * neg1)
  obs_change <- transform_M(pos) - transform_M(neg1) #log2(pos/neg1)
  obs_abundance[is.nan(obs_abundance) | is.infinite(obs_abundance)] <- NA
  obs_change[is.nan(obs_change) | is.infinite(obs_change)] <- NA

  if (empirical_null) {
    null_abundance <- transform_A(neg1) + transform_A(neg2) #log10(neg1 * neg2)
    null_change <- transform_M(neg1) - transform_M(neg2) #log2(neg1/neg2)
    null_abundance[is.nan(null_abundance) | is.infinite(null_abundance)] <- NA
    null_change[is.nan(null_change) | is.infinite(null_change)] <- NA
  } else {
    null_abundance <- rep(NA, nrow(sorted_sample_df))
    null_change <- rep(NA, nrow(sorted_sample_df))
  }

  ma_coords <-
    data.frame(
      taxon_id = sorted_sample_df$taxon_id,
      sample_id = sorted_sample_df$sample_id,
      # TODO: maybe no need for original pos, neg1, neg2??
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

#' Get Data for Ig-Seq MA plot
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
    imputed_taxa = zero_imputation_result$imputed_taxa
  ))
}


#' Ig-Seq MA plot
#' @export
plot_ma <-
  function(ma_plot_data, ellipses = FALSE, type = "facet") {
    ma_non_imputed <- ma_plot_data$plot_data[
      !ma_plot_data$plot_data$taxon_id %in% ma_plot_data$imputed_taxa,
    ]

    ma_imputed <- ma_plot_data$plot_data[
      ma_plot_data$plot_data$taxon_id %in% ma_plot_data$imputed_taxa,
    ]

    jitter_width <- diff(range(ma_non_imputed$A)) / 6
    jitter_x <- min(ma_non_imputed$A) - jitter_width * 3

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
      theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        legend.title = element_text(face = "bold", hjust = 0.5)
      ) +
      scale_color_manual(values = c("darkgray", ggsci::pal_npg()(4)))
    #ggsci::scale_color_npg()

    return(plt)
  }
