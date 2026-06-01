#' Get Slide Z score from Ig Fractions Data
#' @export
get_slide_z =
  function(
    sorted_sample_df, # dataframe from the result of group_sorted_samples()
    positive_fraction_name = "pos",
    first_negative_fraction_name = "neg1",
    second_negative_fraction_name = "neg2",
    window_size = 50,
    empirical_null_distribution = TRUE,
    confidence_levels = NULL, # c(0.95, 0.99, 0.999)
    imputed_taxa = NULL # calculate slide z differently on those taxa
  ) {
    if (empirical_null_distribution & is.null(second_negative_fraction_name)) {
      warning(
        "No second negative fraction furnished, cannot model empirical null ( Ig-.1 vs Ig-.2) distribution...\n"
      )
      empirical_null_distribution = FALSE
    }

    if (
      !all(
        c(
          positive_fraction_name,
          first_negative_fraction_name,
          second_negative_fraction_name
        ) %in%
          colnames(sorted_sample_df)
      )
    ) {
      warning(paste0(
        "Sample ",
        sorted_sample_df$sample_id[1],
        " lacks fraction(s)\n"
      ))
      return(list(
        slide_z = NA,
        ma_coords = data.frame(),
        ellipse_level = NA,
        ellipse_coords = data.frame()
      ))
    }

    ma_coords =
      get_ma_coordinates(
        sorted_sample_df = sorted_sample_df,
        positive_fraction_name = positive_fraction_name,
        first_negative_fraction_name = first_negative_fraction_name,
        second_negative_fraction_name = second_negative_fraction_name
      )

    slide_z = compute_slide_z(
      ma_coords = ma_coords,
      was_imputed = ma_coords$taxon_id %in% imputed_taxa,
      window_size = window_size,
      empirical_null_distribution = empirical_null_distribution
    )
    if (!is.null(confidence_levels)) {
      ellipse_data = get_ellipse_data(
        sorted_sample_df = ma_coords,
        imputed_taxa = imputed_taxa,
        empirical_null_distribution = empirical_null_distribution,
        confidence_levels = confidence_levels
      )
    } else {
      ellipse_data = list(coords = data.frame())
    }

    return(list(
      slide_z = slide_z,
      ma_coords = ma_coords,
      ellipse_level = ellipse_data$levels,
      ellipse_coords = ellipse_data$coords
    ))
  }

#' Compute Slide Z from MA coordinates
#' @export
compute_slide_z = function(
  ma_coords,
  was_imputed = NULL, # boolean vector - whether this row (taxon)
  # has imputed zero(s) and should
  # be considered separately
  window_size = 50,
  empirical_null_distribution = TRUE
) {
  if (is.logical(was_imputed) & length(was_imputed) > 0) {
    imputed_coords = ma_coords[was_imputed, ]
    ma_coords = ma_coords[!was_imputed, ]
  } else {
    imputed_coords = NULL
  }

  # Treat imputed coordinates as a slice apart:
  if (!is.null(imputed_coords)) {
    if (empirical_null_distribution) {
      # aka "slide_z_modern"
      # center and scale each pos vs neg ratio with empirical null (neg vs neg) mean and sd
      slide_z_imputed = (imputed_coords$obs_change -
        mean(imputed_coords$null_change, na.rm = TRUE)) /
        sd(imputed_coords$null_change, na.rm = TRUE)
    } else {
      # aka "slide_z_standard"
      # or center and scale based on the same pos vs neg distribution
      slide_z_imputed = (imputed_coords$obs_change -
        mean(imputed_coords$obs_change, na.rm = TRUE)) /
        sd(imputed_coords$obs_change, na.rm = TRUE)
    }
  } else {
    slide_z_imputed = NULL
  }

  taxa_order = order(ma_coords$obs_abundance, decreasing = TRUE)
  # Make sure that coordinates are sorted by observed abundance
  # TODO: slice will not concern obs_abundance and null_abundance's of the same
  # rank then - probably algorithm should be changed!
  ma_coords = ma_coords[taxa_order, ]

  if (nrow(ma_coords) == 0) {
    warning(paste0("No taxa in MA coordinates, cannot compute slide z score\n"))
    return(NULL)
  }

  # Overlap by the half of the size of the window
  overlap = window_size %/% 2

  # Obtain number of windows by integer division
  n_last_window = nrow(ma_coords) %/% window_size

  # If what's left exceeds the overlap or the window size is bigger then number of rows,
  # add one more window
  # FIXME: what happens if what's left is less then overlap
  if ((nrow(ma_coords) %% window_size > overlap) | (n_last_window == 0)) {
    n_last_window = n_last_window + 1
  }

  slide_z_all = c()

  # Loop through windows
  for (n_window in 1:n_last_window) {
    # Get indices of window with overlap : window_start/_end, to obtain a slice
    # and get indices of window w/o overlap INSIDE THE SLICE: slice_window_start/_end
    # In case of the first window (and if there's more than one window overall)
    if (n_window == 1 && n_window != n_last_window) {
      # overlap only on the right, start at first row
      window_start = 1
      window_end = window_size + overlap
      slice_window_start = 1
      slice_window_end = window_size
    } else if (n_window > 1 && n_window < n_last_window) {
      # If it's the window in between
      # start and end in window_start n_window-th window +- overlap
      # (n_window - 1) * window_size is the end of the previous window
      # so (n_window - 1) * window_size + 1 = start of the current window w/o overlap
      window_start = (n_window - 1) * window_size + 1 - overlap
      window_end = n_window * window_size + overlap
      slice_window_start = overlap + 1
      slice_window_end = overlap + window_size
    } else if (n_window > 1 && n_window == n_last_window) {
      # If it's the last window (and if there's more than one window overall)
      # overlap only on the left, end at the last row
      window_start = (n_window - 1) * window_size + 1 - overlap
      window_end = nrow(ma_coords)
      slice_window_start = overlap + 1
      slice_window_end = window_end - window_start + 1
    } else if (n_window == 1 && n_window == n_last_window) {
      # If there's only one window
      # start with the first and end with the last row
      window_start = 1
      window_end = nrow(ma_coords)
      slice_window_start = 1
      slice_window_end = window_end - window_start + 1
    }

    # taxa_slice = slice window +- overlap, use it to estimate mean and sd
    taxa_slice = ma_coords[window_start:window_end, ]

    if (empirical_null_distribution) {
      # aka "slide_z_modern"
      # center and scale each pos vs neg ratio with empirical null (neg vs neg) mean and sd
      slide_z_slice = (taxa_slice$obs_change -
        mean(taxa_slice$null_change, na.rm = TRUE)) /
        sd(taxa_slice$null_change, na.rm = TRUE)
    } else {
      # aka "slide_z_standard"
      # or center and scale based on the same pos vs neg distribution
      slide_z_slice = (taxa_slice$obs_change -
        mean(taxa_slice$obs_change, na.rm = TRUE)) /
        sd(taxa_slice$obs_change, na.rm = TRUE)
    }

    # We use window +- overlap only to estimate distribution parameters, but we don't keep z scores
    # for taxa from overlap i.e  we apply the distribution only to the window w/o overlap

    # take only Z scores from the window w/o overlap
    slide_z_slice = slide_z_slice[slice_window_start:slice_window_end]

    # APPEND the window result to final Z score vector
    slide_z_all = c(slide_z_all, slide_z_slice)
  }
  # set taxa back to initial order
  # FIXME: verify if it's ok
  slide_z_all = slide_z_all[order(taxa_order)]

  # merge with scores on imputed taxa

  if (!is.null(slide_z_imputed)) {
    slide_z_merged = rep(NA, length(was_imputed))
    slide_z_merged[!was_imputed] = slide_z_all
    slide_z_merged[was_imputed] = slide_z_imputed
    return(slide_z_merged)
  } else {
    return(slide_z_all)
  }
}

#' Get Confidence Ellipse Coordinates and Taxa Confidence Levels
#' @export
get_ellipse_data =
  function(
    sorted_sample_df,
    imputed_taxa = NULL,
    empirical_null_distribution = TRUE,
    confidence_levels # a vector of confidence levels
  ) {
    if (
      all(is.null(sorted_sample_df$null_abundance)) ||
        all(is.na(sorted_sample_df$null_abundance)) ||
        all(is.null(sorted_sample_df$null_change)) ||
        all(is.na(sorted_sample_df$null_change))
    ) {
      if (empirical_null_distribution) {
        warning(paste0(
          "No MA coordinates to model empirical null distribution for ",
          unique(sorted_sample_df$sample_id),
          "...\n"
        ))
      }
      empirical_null_distribution = FALSE
    }

    valid_taxa_obs = !is.na(sorted_sample_df$obs_abundance) &
      !is.na(sorted_sample_df$obs_change) &
      !sorted_sample_df$taxon_id %in% imputed_taxa

    if (empirical_null_distribution) {
      valid_taxa_null = !is.na(sorted_sample_df$null_abundance) &
        !is.na(sorted_sample_df$null_change) &
        !sorted_sample_df$taxon_id %in% imputed_taxa
      # based on neg1 vs neg2, construct Ellipses
      abund_coords = sorted_sample_df$null_abundance[valid_taxa_null]
      change_coords = sorted_sample_df$null_change[valid_taxa_null]
    } else {
      # based on pos vs neg1, construct Ellipses
      abund_coords = sorted_sample_df$obs_abundance[valid_taxa_obs]
      change_coords = sorted_sample_df$obs_change[valid_taxa_obs]
    }

    min_nb_points = 2
    if (length(abund_coords) > min_nb_points) {
      ellipse_data = car::dataEllipse(
        abund_coords,
        change_coords,
        levels = confidence_levels,
        draw = FALSE
      )
    } else {
      warning(paste0(
        "Cannot build ellipse with only ",
        length(abund_coords),
        " points\n"
      ))
      return(list(levels = NULL, coords = data.frame()))
    }

    if (length(confidence_levels) == 1) {
      ellipse_list = list()
      ellipse_list[[as.character(confidence_levels)]] = ellipse_data
    } else {
      ellipse_list = ellipse_data
    }

    ellipse_coords = data.frame()
    for (ellipse_level in names(ellipse_list)) {
      ellipse_coords = rbind(
        ellipse_coords,
        cbind(
          ellipse_list[[ellipse_level]],
          data.frame(
            ellipse_level = rep(
              ellipse_level,
              nrow(ellipse_list[[ellipse_level]])
            )
          ),
          data.frame(
            sample_id = rep(
              sorted_sample_df$sample_id[1],
              nrow(ellipse_list[[ellipse_level]])
            )
          )
        )
      )
    }

    # Get a boolean indicator whether the point from pos vs neg1 is in the ellipse

    is_outside_ellipse = data.frame("ns" = rep(TRUE, nrow(sorted_sample_df)))

    for (confidence_level in as.character(sort(confidence_levels))) {
      is_outside_ellipse[[confidence_level]] =
        !sp::point.in.polygon(
          sorted_sample_df$obs_abundance,
          sorted_sample_df$obs_change,
          ellipse_list[[confidence_level]][, 1],
          ellipse_list[[confidence_level]][, 2]
        )
    }
    # Get a maximum confidence level for each taxon
    max_confidence_level = factor(
      names(is_outside_ellipse)[rowSums(is_outside_ellipse)],
      # put "ns" as the first level
      levels = unique(c("ns", names(is_outside_ellipse)))
    )

    names(max_confidence_level) = rownames(sorted_sample_df)
    max_confidence_level[!valid_taxa_obs] = NA

    # # Treat imputed taxa separately - instead of ellipses, build confidence intervals
    #
    # valid_taxa_obs = !is.na(sorted_sample_df$obs_abundance) & !is.na(sorted_sample_df$obs_change) &
    #   sorted_sample_df$taxon_id %in% imputed_taxa
    #
    # if(empirical_null_distribution){
    #   valid_taxa_null = !is.na(sorted_sample_df$null_abundance) & !is.na(sorted_sample_df$null_change) &
    #     sorted_sample_df$taxon_id %in% imputed_taxa
    #
    #   change_coords = sorted_sample_df$null_change[valid_taxa_null]
    #
    # }else{
    #   change_coords = sorted_sample_df$obs_change[valid_taxa_obs]
    # }
    #
    # min_nb_points = 2
    # if(length(abund_coords) > min_nb_points){
    #   z_score = ( sorted_sample_df$obs_change[valid_taxa_obs] - mean(change_coords, na.rm = TRUE))/sd(change_coords, na.rm = TRUE)
    #   max_confidence_level[valid_taxa_obs] =
    #
    # }else{
    #   warning(paste0("Cannot build confidence intervals with only ", length(abund_coords), " points\n"))
    # }

    return(list(levels = max_confidence_level, coords = ellipse_coords))
  }
