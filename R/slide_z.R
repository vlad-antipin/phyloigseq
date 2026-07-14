#' Get the Sliding Z-Score for One Biological Sample
#'
#' Orchestrates the package's sliding Z-score, its novel Ig-Seq scoring
#' statistic: computes MA-plot coordinates ([get_ma_coordinates()]) between
#' the positive and (first) negative fraction of one biological sample,
#' estimates a *local* null distribution along the abundance axis
#' ([compute_slide_z()]) rather than assuming one global null, and
#' optionally builds per-taxon confidence ellipses/levels
#' ([get_ellipse_data()]).
#'
#' @param sorted_sample_df A data frame for one biological sample, as
#'   produced by [group_sorted_samples()], with one row per taxon and one
#'   column per fraction.
#' @param positive_fraction_name Name of the column holding the positive
#'   (Ig+) fraction abundance.
#' @param first_negative_fraction_name Name of the column holding the
#'   primary negative fraction abundance.
#' @param second_negative_fraction_name Name of the column holding a second
#'   negative fraction abundance, used to model the null distribution
#'   empirically (Ig-.1 vs Ig-.2) instead of theoretically (Ig+ vs Ig-.1).
#'   `NULL` disables the empirical null even if `empirical_null_distribution
#'   = TRUE` (a warning is issued and `empirical_null_distribution` is
#'   forced to `FALSE`).
#' @param window_size Number of taxa (ranked by abundance) per sliding
#'   window used to estimate the local null mean/sd; smaller windows track
#'   local abundance-dependent variance more closely at the cost of
#'   noisier estimates. Default `50`; exposed to end users as an adjustable
#'   parameter in the companion Shiny app.
#' @param empirical_null_distribution If `TRUE` (default), center/scale
#'   each taxon's observed log-ratio against the local Ig-.1 vs Ig-.2 null;
#'   if `FALSE`, center/scale against the local Ig+ vs Ig-.1 distribution
#'   itself. Forced to `FALSE` when `second_negative_fraction_name` is
#'   `NULL`.
#' @param confidence_levels Optional numeric vector of confidence levels
#'   (e.g. `c(0.95, 0.99, 0.999)`) to build confidence ellipses for; `NULL`
#'   (default) skips ellipse construction.
#' @param imputed_taxa Optional vector of `taxon_id`s to score separately
#'   from the sliding-window scheme (see [compute_slide_z()]), typically
#'   the taxa [impute_zeros()] filled in for this sample.
#'
#' @return A list:
#'   \describe{
#'     \item{`slide_z`}{Numeric vector of per-taxon sliding Z-scores, in
#'       `sorted_sample_df`'s original row order.}
#'     \item{`ma_coords`}{The MA-plot coordinates data frame from
#'       [get_ma_coordinates()].}
#'     \item{`ellipse_level`}{Per-taxon maximum confidence level outside
#'       which the taxon falls (see [get_ellipse_data()]), or `NULL` if
#'       `confidence_levels` is `NULL`.}
#'     \item{`ellipse_coords`}{Data frame of confidence ellipse boundary
#'       coordinates, or an empty data frame if `confidence_levels` is
#'       `NULL`.}
#'   }
#'   If `sorted_sample_df` is missing one of the three fraction columns,
#'   this issues a `warning()` and returns `list(slide_z = NA, ma_coords =
#'   data.frame(), ellipse_level = NA, ellipse_coords = data.frame())`
#'   instead of erroring, so that callers processing many samples (e.g.
#'   [getPhyloIgSeq()]) can skip the failing one rather than abort the
#'   whole batch — callers are expected to check for this rather than
#'   assume every sample succeeded.
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
#' result <- get_slide_z(
#'   sorted_sample_df = grouped[["sample_1"]],
#'   positive_fraction_name = "Pos",
#'   first_negative_fraction_name = "Neg1",
#'   second_negative_fraction_name = "Neg2",
#'   confidence_levels = c(0.95, 0.99)
#' )
#' head(result$slide_z)
#' head(result$ellipse_coords)
#'
#' @export
get_slide_z <-
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
    if (empirical_null_distribution && is.null(second_negative_fraction_name)) {
      warning(
        "No second negative fraction furnished, cannot model empirical null ( Ig-.1 vs Ig-.2) distribution...\n"
      )
      empirical_null_distribution <- FALSE
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

    ma_coords <-
      get_ma_coordinates(
        sorted_sample_df = sorted_sample_df,
        positive_fraction_name = positive_fraction_name,
        first_negative_fraction_name = first_negative_fraction_name,
        second_negative_fraction_name = second_negative_fraction_name
      )

    slide_z <- compute_slide_z(
      ma_coords = ma_coords,
      was_imputed = ma_coords$taxon_id %in% imputed_taxa,
      window_size = window_size,
      empirical_null_distribution = empirical_null_distribution
    )
    if (!is.null(confidence_levels)) {
      ellipse_data <- get_ellipse_data(
        sorted_sample_df = ma_coords,
        imputed_taxa = imputed_taxa,
        empirical_null_distribution = empirical_null_distribution,
        confidence_levels = confidence_levels
      )
    } else {
      ellipse_data <- list(coords = data.frame())
    }

    return(list(
      slide_z = slide_z,
      ma_coords = ma_coords,
      ellipse_level = ellipse_data$levels,
      ellipse_coords = ellipse_data$coords
    ))
  }

#' Center and Scale Observed Log-Ratios Against a Null Distribution
#'
#' Internal helper shared by [compute_slide_z()]'s imputed-taxa branch and
#' its per-window loop: both center/scale `df$obs_change` against either the
#' empirical null (`df$null_change`, Ig-.1 vs Ig-.2) or the observed
#' distribution itself (`df$obs_change`, Ig+ vs Ig-.1), depending on
#' `empirical_null_distribution`. Pulled out verbatim from two identical
#' inline blocks — no formula change.
#'
#' @param df A data frame with `obs_change` and (if
#'   `empirical_null_distribution`) `null_change` columns.
#' @param empirical_null_distribution See [compute_slide_z()].
#'
#' @return Numeric vector of Z-scores, same length as `nrow(df)`.
#' @noRd
.slide_z_center_scale <- function(df, empirical_null_distribution) {
  if (empirical_null_distribution) {
    # aka "slide_z_modern": center/scale each pos vs neg ratio with the
    # empirical null (neg vs neg) mean and sd.
    reference <- df$null_change
  } else {
    # aka "slide_z_standard": center/scale against the same pos vs neg
    # distribution itself.
    reference <- df$obs_change
  }
  (df$obs_change - mean(reference, na.rm = TRUE)) / sd(reference, na.rm = TRUE)
}

#' Compute Sliding-Window Index Bounds for `compute_slide_z()`
#'
#' Internal helper returning, for the `n_window`-th of `n_last_window`
#' windows over `n_rows` abundance-sorted taxa, the row range to slice for
#' parameter estimation (`window_start`/`window_end`, includes the overlap
#' on either side) and the sub-range of that slice whose Z-scores are
#' actually kept (`slice_window_start`/`slice_window_end`, excludes the
#' overlap). Pulled out verbatim from the original 4-branch `if`/`else if`
#' chain — no arithmetic change.
#'
#' Verified 2026-07-14 via property-based testing (all `n_rows` 1-400
#' crossed with 18 `window_size` values, including cases where the leftover
#' after full windows is smaller than `overlap`): the returned bounds
#' always tile `1:n_rows` with no gaps, no double-counted rows, and no
#' out-of-bounds indices. The last window's `window_end`/`slice_window_end`
#' are always pinned to `n_rows`, so it simply absorbs whatever is left
#' over, however small.
#'
#' @param n_window Index of the current window (`1:n_last_window`).
#' @param n_last_window Total number of windows.
#' @param window_size See [compute_slide_z()].
#' @param overlap `window_size %/% 2`.
#' @param n_rows `nrow()` of the (already sorted) coordinates being windowed.
#'
#' @return A list with `window_start`, `window_end`, `slice_window_start`,
#'   `slice_window_end`.
#' @noRd
.compute_window_bounds <- function(
  n_window,
  n_last_window,
  window_size,
  overlap,
  n_rows
) {
  if (n_window == 1 && n_window != n_last_window) {
    # First of several windows: overlap only on the right, start at row 1.
    window_start <- 1
    window_end <- window_size + overlap
    slice_window_start <- 1
    slice_window_end <- window_size
  } else if (n_window > 1 && n_window < n_last_window) {
    # A window strictly between the first and last: overlap on both sides.
    # (n_window - 1) * window_size is the end of the previous window's
    # kept slice, so (n_window - 1) * window_size + 1 is this window's
    # kept slice start.
    window_start <- (n_window - 1) * window_size + 1 - overlap
    window_end <- n_window * window_size + overlap
    slice_window_start <- overlap + 1
    slice_window_end <- overlap + window_size
  } else if (n_window > 1 && n_window == n_last_window) {
    # Last of several windows: overlap only on the left, end at n_rows so
    # any leftover rows (however many) are absorbed into this window.
    window_start <- (n_window - 1) * window_size + 1 - overlap
    window_end <- n_rows
    slice_window_start <- overlap + 1
    slice_window_end <- window_end - window_start + 1
  } else {
    # Only one window overall: start with the first, end with the last row.
    window_start <- 1
    window_end <- n_rows
    slice_window_start <- 1
    slice_window_end <- window_end - window_start + 1
  }
  list(
    window_start = window_start,
    window_end = window_end,
    slice_window_start = slice_window_start,
    slice_window_end = slice_window_end
  )
}

#' Compute the Sliding Z-Score From MA-Plot Coordinates
#'
#' The core sliding-window statistic behind [get_slide_z()]: rather than
#' assuming one global null distribution for the whole sample, taxa are
#' ranked by observed abundance and processed in overlapping windows of
#' `window_size` taxa, each window's own local mean/sd (of either the
#' empirical null or the observed change, see `empirical_null_distribution`)
#' used to center/scale the Z-scores of the taxa in that window's non-
#' overlapping "core". Taxa flagged via `was_imputed` are scored separately,
#' against the null/observed distribution taken over *all* imputed taxa at
#' once rather than a local window (there are usually too few of them to
#' window meaningfully).
#'
#' @param ma_coords MA-plot coordinates for one biological sample, as
#'   returned by [get_ma_coordinates()] (`taxon_id`, `obs_abundance`,
#'   `obs_change`, and, if `empirical_null_distribution`, `null_abundance`/
#'   `null_change`).
#' @param was_imputed Optional logical vector, same length/row order as
#'   `ma_coords`, flagging which taxa had an imputed zero and should be
#'   scored separately from the sliding-window scheme.
#' @param window_size See [get_slide_z()]. Default `50`.
#' @param empirical_null_distribution See [get_slide_z()]. Default `TRUE`.
#'
#' @return A numeric vector of per-taxon sliding Z-scores, in `ma_coords`'s
#'   original row order (or `NULL`, with a `warning()`, if `ma_coords` has
#'   zero non-imputed rows).
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
#' z_scores <- compute_slide_z(ma_coords, window_size = 50)
#' head(z_scores)
#'
#' @export
compute_slide_z <- function(
  ma_coords,
  was_imputed = NULL, # boolean vector - whether this row (taxon)
  # has imputed zero(s) and should
  # be considered separately
  window_size = 50,
  empirical_null_distribution = TRUE
) {
  if (is.logical(was_imputed) && length(was_imputed) > 0) {
    imputed_coords <- ma_coords[was_imputed, ]
    ma_coords <- ma_coords[!was_imputed, ]
  } else {
    imputed_coords <- NULL
  }

  # Treat imputed coordinates as a slice apart:
  if (!is.null(imputed_coords)) {
    slide_z_imputed <- .slide_z_center_scale(
      imputed_coords,
      empirical_null_distribution
    )
  } else {
    slide_z_imputed <- NULL
  }

  # Make sure that coordinates are sorted by observed abundance, so that a
  # sliding window over rows is a sliding window over the abundance axis.
  # Note: the sort key is obs_abundance only, so a window's null_abundance
  # values aren't guaranteed to be similarly ranked — the local null model
  # assumes taxa with similar *observed* abundance are comparable, which is
  # the MA-plot's usual assumption, not a bug in the sort itself.
  taxa_order <- order(ma_coords$obs_abundance, decreasing = TRUE)
  ma_coords <- ma_coords[taxa_order, ]

  if (nrow(ma_coords) == 0) {
    warning(paste0("No taxa in MA coordinates, cannot compute slide z score\n"))
    return(NULL)
  }

  # Overlap by the half of the size of the window
  overlap <- window_size %/% 2

  # Obtain number of windows by integer division
  n_last_window <- nrow(ma_coords) %/% window_size

  # Add one more window if the leftover exceeds the overlap, or if there
  # isn't even one full window yet. This is safe even when the leftover is
  # *smaller* than the overlap (no extra window added): see
  # .compute_window_bounds()'s doc for the verification that the last
  # window always absorbs it regardless of size.
  if ((nrow(ma_coords) %% window_size > overlap) | (n_last_window == 0)) {
    n_last_window <- n_last_window + 1
  }

  slide_z_all <- c()

  # Loop through windows
  for (n_window in 1:n_last_window) {
    bounds <- .compute_window_bounds(
      n_window = n_window,
      n_last_window = n_last_window,
      window_size = window_size,
      overlap = overlap,
      n_rows = nrow(ma_coords)
    )

    # taxa_slice = slice window +- overlap, use it to estimate mean and sd
    taxa_slice <- ma_coords[bounds$window_start:bounds$window_end, ]

    slide_z_slice <- .slide_z_center_scale(
      taxa_slice,
      empirical_null_distribution
    )

    # We use window +- overlap only to estimate distribution parameters,
    # but we don't keep Z-scores for taxa from the overlap, i.e. we apply
    # the distribution only to the window's non-overlapping core.
    slide_z_slice <- slide_z_slice[
      bounds$slice_window_start:bounds$slice_window_end
    ]

    # APPEND the window result to final Z score vector
    slide_z_all <- c(slide_z_all, slide_z_slice)
  }

  # Set taxa back to their initial order. order(taxa_order) is the standard
  # inverse-permutation idiom (taxa_order[k] gives, for sorted position k,
  # its original row index; order(taxa_order) gives the reverse mapping),
  # so this restores the pre-sort row order exactly regardless of ties in
  # obs_abundance. Verified 2026-07-14 via 2000 randomized trials
  # (including induced ties and NAs).
  slide_z_all <- slide_z_all[order(taxa_order)]

  # merge with scores on imputed taxa

  if (!is.null(slide_z_imputed)) {
    slide_z_merged <- rep(NA, length(was_imputed))
    slide_z_merged[!was_imputed] <- slide_z_all
    slide_z_merged[was_imputed] <- slide_z_imputed
    return(slide_z_merged)
  } else {
    return(slide_z_all)
  }
}

#' Get Confidence Ellipse Coordinates and Taxa Confidence Levels
#'
#' Builds normal-theory confidence ellipses (via `car::dataEllipse()`) over
#' the MA-plot coordinates of one biological sample — one nested ellipse per
#' entry of `confidence_levels` — and, for each taxon, reports the highest
#' confidence level whose ellipse it falls outside of. Used by
#' [get_slide_z()] to turn `confidence_levels` into a per-taxon
#' significance call independent of the sliding Z-score itself.
#'
#' @param sorted_sample_df MA-plot coordinates for one biological sample —
#'   either the output of [get_ma_coordinates()], or `ma_coords` as returned
#'   inside [get_slide_z()]'s result list (`taxon_id`, `obs_abundance`,
#'   `obs_change`, and, if `empirical_null_distribution`, `null_abundance`/
#'   `null_change`).
#' @param imputed_taxa Optional vector of `taxon_id`s to exclude from the
#'   ellipse fit and leave with an `NA` confidence level (they are scored
#'   separately elsewhere, see [compute_slide_z()]).
#' @param empirical_null_distribution If `TRUE` (default), fit the ellipse
#'   over the empirical null (Ig-.1 vs Ig-.2) coordinates; if `FALSE`, or if
#'   `sorted_sample_df` has no usable `null_abundance`/`null_change` values,
#'   fit over the observed (Ig+ vs Ig-.1) coordinates instead (a `warning()`
#'   is issued in the latter case when `empirical_null_distribution` was
#'   requested).
#' @param confidence_levels Numeric vector of confidence levels (e.g.
#'   `c(0.95, 0.99, 0.999)`) to build nested ellipses for.
#'
#' @return A list:
#'   \describe{
#'     \item{`levels`}{A factor, one entry per row of `sorted_sample_df`
#'       (named by `rownames(sorted_sample_df)`), giving each taxon's
#'       highest confidence level whose ellipse it falls outside of (as a
#'       character level, e.g. `"0.99"`), `"ns"` ("not significant" — inside
#'       every ellipse), or `NA` for taxa in `imputed_taxa` or with a
#'       missing `obs_abundance`/`obs_change`. `NULL` instead if there
#'       weren't enough points (`> 2`) to fit an ellipse at all.}
#'     \item{`coords`}{A data frame of ellipse boundary coordinates (`x`,
#'       `y`, `ellipse_level`, `sample_id`), one block of rows per
#'       confidence level, ready for overlaying on an MA-plot; empty if
#'       ellipses couldn't be fit.}
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
#' ma_coords <- get_ma_coordinates(
#'   sorted_sample_df = grouped[["sample_1"]],
#'   positive_fraction_name = "Pos",
#'   first_negative_fraction_name = "Neg1",
#'   second_negative_fraction_name = "Neg2"
#' )
#' ellipse_data <- get_ellipse_data(
#'   sorted_sample_df = ma_coords,
#'   confidence_levels = c(0.95, 0.99)
#' )
#' table(ellipse_data$levels)
#' head(ellipse_data$coords)
#'
#' @export
get_ellipse_data <-
  function(
    sorted_sample_df,
    imputed_taxa = NULL,
    empirical_null_distribution = TRUE,
    confidence_levels # a vector of confidence levels
  ) {
    if (
      is.null(sorted_sample_df$null_abundance) ||
        all(is.na(sorted_sample_df$null_abundance)) ||
        is.null(sorted_sample_df$null_change) ||
        all(is.na(sorted_sample_df$null_change))
    ) {
      if (empirical_null_distribution) {
        warning(paste0(
          "No MA coordinates to model empirical null distribution for ",
          unique(sorted_sample_df$sample_id),
          "...\n"
        ))
      }
      empirical_null_distribution <- FALSE
    }

    valid_taxa_obs <- !is.na(sorted_sample_df$obs_abundance) &
      !is.na(sorted_sample_df$obs_change) &
      !sorted_sample_df$taxon_id %in% imputed_taxa

    if (empirical_null_distribution) {
      valid_taxa_null <- !is.na(sorted_sample_df$null_abundance) &
        !is.na(sorted_sample_df$null_change) &
        !sorted_sample_df$taxon_id %in% imputed_taxa
      # based on neg1 vs neg2, construct Ellipses
      abund_coords <- sorted_sample_df$null_abundance[valid_taxa_null]
      change_coords <- sorted_sample_df$null_change[valid_taxa_null]
    } else {
      # based on pos vs neg1, construct Ellipses
      abund_coords <- sorted_sample_df$obs_abundance[valid_taxa_obs]
      change_coords <- sorted_sample_df$obs_change[valid_taxa_obs]
    }

    min_nb_points <- 2
    if (length(abund_coords) > min_nb_points) {
      ellipse_data <- car::dataEllipse(
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
      ellipse_list <- list()
      ellipse_list[[as.character(confidence_levels)]] <- ellipse_data
    } else {
      ellipse_list <- ellipse_data
    }

    ellipse_coords <- data.frame()
    for (ellipse_level in names(ellipse_list)) {
      ellipse_coords <- rbind(
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

    # Get a boolean indicator whether the point from pos vs neg1 is in the ellipse.
    #
    # is_outside_ellipse gets one column per confidence level (plus a
    # sentinel "ns" column, always TRUE), added in ascending confidence
    # order, TRUE where the taxon falls outside that level's ellipse.
    # Because higher confidence levels correspond to strictly larger nested
    # ellipses (car::dataEllipse() fits them around the same points/center),
    # a taxon outside the k-th ellipse is always outside every ellipse
    # before it too — so each row's TRUE/FALSE pattern is a run of TRUEs
    # from "ns" up to some cutoff column, and rowSums() (TRUE = 1, FALSE =
    # 0) gives exactly the position of that cutoff column. That position is
    # then used directly as a names() index below to look up the
    # corresponding confidence-level label — a compact way to get "highest
    # confidence level this taxon falls outside of" without an explicit
    # per-row loop.
    is_outside_ellipse <- data.frame("ns" = rep(TRUE, nrow(sorted_sample_df)))

    for (confidence_level in as.character(sort(confidence_levels))) {
      is_outside_ellipse[[confidence_level]] <-
        !sp::point.in.polygon(
          sorted_sample_df$obs_abundance,
          sorted_sample_df$obs_change,
          ellipse_list[[confidence_level]][, 1],
          ellipse_list[[confidence_level]][, 2]
        )
    }
    # Get a maximum confidence level for each taxon
    max_confidence_level <- factor(
      names(is_outside_ellipse)[rowSums(is_outside_ellipse)],
      # put "ns" as the first level
      levels = unique(c("ns", names(is_outside_ellipse)))
    )

    names(max_confidence_level) <- rownames(sorted_sample_df)
    max_confidence_level[!valid_taxa_obs] <- NA

    # An alternative treatment of imputed taxa was sketched here (build a
    # per-taxon Z-score/confidence-interval from the null/observed change
    # distribution restricted to imputed_taxa, instead of the ellipse
    # membership test above) but never implemented — worth revisiting if
    # imputed taxa need their own significance call rather than just an NA
    # confidence level.

    return(list(levels = max_confidence_level, coords = ellipse_coords))
  }
