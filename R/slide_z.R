#' Get the Sliding Z-Score for One Biological Sample
#'
#' Orchestrates the package's sliding Z-score, its novel IgSeq scoring
#' statistic: computes MA-plot coordinates ([get_ma_coordinates()]) between
#' the positive and (first) negative fraction of one biological sample,
#' estimates a *local* null distribution along the abundance axis
#' ([compute_slide_z()]) rather than assuming one global null, and
#' optionally builds per-taxon confidence ellipses/levels
#' ([get_ellipse_data()]).
#'
#' @section Saturated taxa under the purity-corrected axis:
#' A taxon the un-mixing has to clamp (`obs_in_cone = FALSE`) gets a change value of
#' `±(log2(pool) - log2(pool_prior))` -- one pool is 0, so what is left is the
#' regularization scale, not a measurement. Two consequences, both measured on a real
#' cohort rather than argued:
#'
#' * On an MA plot these taxa fall on two straight lines forming a **V opening to the
#'   right**, with slopes `±log2(10)/2 ≈ ±1.66` (observed `+1.77`/`-1.76`, R² 0.89),
#'   because the change is a log of one pool while the abundance axis is a log of the
#'   product of two. Taxa inside the cone show no such structure (slope 0.00, R² 0.000).
#'   The V is the saturation boundary, drawn by the model's own geometry -- it is not a
#'   feature of the data and should not be read as one.
#' * Their `|Z|` is large **by construction**: a change value fixed by `pool_prior`,
#'   divided by a null that measures only technical-replicate noise. On that cohort
#'   *every* out-of-cone taxon cleared `|Z| > 1.96` (median 7.3, versus 1.5 and a 37%
#'   significant rate for taxa inside the cone), and scaling `pool_prior` over three
#'   orders of magnitude moved the largest `|Z|` around without shrinking it.
#'
#' So the *ranking* among saturated taxa is meaningful and they are deliberately kept
#' rather than dropped, but their magnitude is a floor and a `|Z|` threshold does not
#' test anything for them. Split them out with `obs_in_cone` before thresholding.
#' [plot_slide_z()] draws them hollow for the same reason. Note this is a property of
#' point-estimating an unbounded quantity from noisy counts, not of the clamp: shrinking
#' the underlying enrichment toward the null instead would address it, and is not
#' implemented.
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
#' @param empirical_null_distribution If `TRUE` (default), scale each taxon's
#'   observed log-ratio by the local Ig-.1 vs Ig-.2 null's sd; if `FALSE`, by the
#'   local Ig+ vs Ig-.1 distribution's own sd. Forced to `FALSE` when
#'   `second_negative_fraction_name` is `NULL`. This chooses the *scale* only --
#'   `center_on` chooses the location, and the two are separate decisions.
#' @param center_on Which distribution supplies the local *location*: `"reference"`
#'   (default, the mean of whichever distribution supplies the scale — at matched depth
#'   this keeps `0` meaning "coated at the community's own Ig+ frequency") or
#'   `"observed"` (the median Ig+ vs Ig-.1 change, which re-anchors on the typical
#'   taxon). Also shifts the confidence ellipses to match, so both significance calls
#'   agree. See [compute_slide_z()]'s "Where the centre comes from" — the two answer
#'   different questions and `"observed"` gives up the tie to the pre-sort frequency.
#' @param confidence_levels Optional numeric vector of confidence levels
#'   (e.g. `c(0.95, 0.99, 0.999)`) to build confidence ellipses for; `NULL`
#'   (default) skips ellipse construction.
#' @param imputed_taxa Optional vector of `taxon_id`s to score as a population apart
#'   from the sliding-window scheme (see [compute_slide_z()]), typically
#'   the taxa [impute_zeros()] filled in for this sample.
#' @param change_transform Which change (`M`) axis to score on: `"log_ratio"`
#'   (default) or `"purity_corrected"`. See [get_ma_coordinates()] for the two axes.
#'   Because the Z-score centers and scales, only a *non-affine* axis changes the
#'   result -- so `"purity_corrected"` differs from the default in the tails, where
#'   the un-mixing is curved, and barely at all near the null, where its local slope
#'   is the constant `1/(pos_ig_freq - neg_ig_freq)`. **Read the note on saturated
#'   taxa below before thresholding its output.**
#' @param pos_ig_freq,neg_ig_freq,cone_policy,pool_prior Passed to
#'   [get_ma_coordinates()]; only read under `change_transform = "purity_corrected"`,
#'   which needs the two Ig+ frequencies but *not* a pre-sort fraction.
#' @param was_zero Passed to [get_ma_coordinates()]: the `was_zero` element of
#'   [impute_zeros()]'s result, marking which counts were really zero before
#'   imputation. **Pass it whenever `sorted_sample_df` has been imputed.** Taxa
#'   whose Ig+ and Ig-.1 were both zero then score `NA` instead of receiving one
#'   shared, manufactured Z-score, and taxa whose two negative fractions were both
#'   zero stop contributing an exact 0 to the null. Both effects are large: on
#'   `ps_igseq` under `pseudo_count`, 13-43% of the imputed population's reference
#'   values were such manufactured zeros, and removing them widens the pooled null
#'   sd by 7-31%.
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
#'     \item{`change_transform`}{Which change axis was used, so callers storing
#'       several can tell them apart.}
#'   }
#'   If `sorted_sample_df` is missing one of the three fraction columns,
#'   this issues a `warning()` and returns `list(slide_z = NA, ma_coords =
#'   data.frame(), ellipse_level = NA, ellipse_coords = data.frame(),
#'   change_transform = <requested>)`
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
#'   fraction_ids = c("pos", "neg1", "neg2")
#' )
#' result <- get_slide_z(
#'   sorted_sample_df = grouped[["sample_1"]],
#'   positive_fraction_name = "pos",
#'   first_negative_fraction_name = "neg1",
#'   second_negative_fraction_name = "neg2",
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
    empirical_null_distribution = TRUE, # picks the SCALE
    center_on = c("reference", "observed"), # picks the LOCATION
    confidence_levels = NULL, # c(0.95, 0.99, 0.999)
    imputed_taxa = NULL, # calculate slide z differently on those taxa
    change_transform = c("log_ratio", "purity_corrected"),
    pos_ig_freq = NULL, # "p", required by "purity_corrected"
    neg_ig_freq = NULL, # "q", -//-
    cone_policy = c("regularize", "na"),
    pool_prior = NULL,
    was_zero = NULL # pre-imputation zero pattern, see impute_zeros()
  ) {
    change_transform <- match.arg(change_transform)
    cone_policy <- match.arg(cone_policy)
    center_on <- match.arg(center_on)
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
        ellipse_coords = data.frame(),
        change_transform = change_transform
      ))
    }

    ma_coords <-
      get_ma_coordinates(
        sorted_sample_df = sorted_sample_df,
        positive_fraction_name = positive_fraction_name,
        first_negative_fraction_name = first_negative_fraction_name,
        second_negative_fraction_name = second_negative_fraction_name,
        change_transform = change_transform,
        pos_ig_freq = pos_ig_freq,
        neg_ig_freq = neg_ig_freq,
        cone_policy = cone_policy,
        pool_prior = pool_prior,
        was_zero = was_zero
      )

    slide_z <- compute_slide_z(
      ma_coords = ma_coords,
      was_imputed = ma_coords$taxon_id %in% imputed_taxa,
      window_size = window_size,
      empirical_null_distribution = empirical_null_distribution,
      center_on = center_on
    )
    if (!is.null(confidence_levels)) {
      ellipse_data <- get_ellipse_data(
        sorted_sample_df = ma_coords,
        imputed_taxa = imputed_taxa,
        empirical_null_distribution = empirical_null_distribution,
        center_on = center_on,
        confidence_levels = confidence_levels
      )
    } else {
      ellipse_data <- list(coords = data.frame())
    }

    return(list(
      slide_z = slide_z,
      ma_coords = ma_coords,
      ellipse_level = ellipse_data$levels,
      ellipse_coords = ellipse_data$coords,
      change_transform = change_transform
    ))
  }

#' Center and Scale Observed Log-Ratios Against a Null Distribution
#'
#' Internal helper shared by [compute_slide_z()]'s imputed-taxa branch and
#' its per-window loop: both center/scale `df$obs_change` against a location and a
#' scale estimated by `.slide_z_params()`.
#'
#' @param df A data frame with `obs_change` and (if
#'   `empirical_null_distribution`) `null_change` columns.
#' @param empirical_null_distribution,center_on See [compute_slide_z()].
#' @param estimate_from Optional logical vector, `nrow(df)` long, restricting which
#'   rows may contribute to the scale. The numerator is never restricted, so every
#'   row still gets a score on the same scale. `NULL` (default) uses every row.
#' @param center_from As `estimate_from`, but for the location; `NULL` reuses
#'   `estimate_from`.
#'
#' @return Numeric vector of Z-scores, same length as `nrow(df)`.
#' @noRd
.slide_z_center_scale <- function(
  df,
  empirical_null_distribution,
  estimate_from = NULL,
  center_on = "reference",
  center_from = NULL
) {
  params <- .slide_z_params(
    df,
    empirical_null_distribution,
    estimate_from,
    center_on,
    center_from
  )
  (df$obs_change - params$center) / params$scale
}

#' Estimate a Slide-Z Null's Location and Scale
#'
#' The single point where slide-Z's reference distribution is estimated. Split out from
#' `.slide_z_center_scale()` so callers can inspect how many rows actually contributed
#' before deciding whether the estimate is usable (see `.slide_z_windowed()`).
#'
#' Location and scale are estimated **separately**, because they are not properties of
#' the same distribution -- see [compute_slide_z()]'s "Where the centre comes from"
#' section for why, and what each choice assumes.
#'
#' @param df,empirical_null_distribution,estimate_from,center_on,center_from See
#'   `.slide_z_center_scale()`.
#'
#' @return A list with `center`, `scale`, `n` (non-`NA` values behind the scale) and
#'   `n_center` (non-`NA` values behind the location). Under
#'   `center_on = "reference"` the two counts are equal by construction: both come
#'   from the same vector.
#' @noRd
.slide_z_params <- function(
  df,
  empirical_null_distribution,
  estimate_from = NULL,
  center_on = c("reference", "observed"),
  center_from = NULL
) {
  center_on <- match.arg(center_on)

  usable <- function(values, mask) {
    if (is.logical(mask) && length(mask) == length(values)) {
      values <- values[mask]
    }
    values[!is.na(values)]
  }

  if (empirical_null_distribution) {
    # aka "slide_z_modern": scale each pos vs neg ratio by the empirical null
    # (neg vs neg) sd.
    reference <- df$null_change
  } else {
    # aka "slide_z_standard": scale against the same pos vs neg distribution
    # itself.
    reference <- df$obs_change
  }
  reference <- usable(reference, estimate_from)

  if (center_on == "reference") {
    # Legacy behaviour: location and scale both read off the reference, so under
    # the empirical null the centre is the mean of the Ig-.1 vs Ig-.2 pair.
    location <- reference
    center <- if (length(location) > 0) mean(location) else NA_real_
  } else {
    # The location belongs to the comparison being scored, not to the pair that
    # supplies the scale. Median rather than mean so the coated taxa the score is
    # meant to find cannot drag the point they are measured from.
    location <- usable(
      df$obs_change,
      if (is.null(center_from)) estimate_from else center_from
    )
    center <- if (length(location) > 0) stats::median(location) else NA_real_
  }

  list(
    center = center,
    scale = if (length(reference) > 1) sd(reference) else NA_real_,
    n = length(reference),
    n_center = length(location)
  )
}

#' Score One Population With the Sliding-Window Scheme
#'
#' The sliding-window engine, applied to the taxa with no imputed zero -- the only
#' population [compute_slide_z()] windows, since imputed taxa get a pooled estimate
#' instead (see there). Ranks rows by observed abundance, tiles them into overlapping
#' windows via `.compute_window_bounds()`, estimates each window's reference from its
#' `estimate_from` rows only, and scores every row in the window's non-overlapping core
#' against it.
#'
#' A window whose location or scale rests on fewer than `min_estimable` usable values
#' falls back, for that quantity alone, to the estimate pooled over this population's
#' whole `estimate_from`/`center_from` set, rather than emitting `NA` or a scale
#' derived from one or two points. The two fall back independently because they can
#' rest on different vectors (see `.slide_z_params()`); under
#' `center_on = "reference"` they rest on the same one, so they always fall back
#' together and the behaviour is the historical one.
#'
#' @param coords MA-plot coordinates for one population, in any row order.
#' @param estimate_from Logical vector, `nrow(coords)` long: which rows may contribute
#'   to a window's scale.
#' @param center_from As `estimate_from`, but for the location; `NULL` reuses
#'   `estimate_from`.
#' @param window_size,empirical_null_distribution,center_on See [compute_slide_z()].
#' @param min_estimable Minimum usable values for a window's own estimate to be
#'   trusted.
#'
#' @return A list with `slide_z` (numeric, `coords`'s original row order) and
#'   `n_fallback`, how many windows had to use the pooled estimate.
#' @noRd
.slide_z_windowed <- function(
  coords,
  estimate_from,
  window_size,
  empirical_null_distribution,
  center_on = "reference",
  center_from = NULL,
  min_estimable = 3L
) {
  n_rows <- nrow(coords)
  if (n_rows == 0) {
    return(list(slide_z = numeric(0), n_fallback = 0L))
  }
  if (is.null(center_from)) {
    center_from <- estimate_from
  }

  # Make sure that coordinates are sorted by observed abundance, so that a
  # sliding window over rows is a sliding window over the abundance axis.
  # Note: the sort key is obs_abundance only, so a window's null_abundance
  # values aren't guaranteed to be similarly ranked — the local null model
  # assumes taxa with similar *observed* abundance are comparable, which is
  # the MA-plot's usual assumption, not a bug in the sort itself.
  taxa_order <- order(coords$obs_abundance, decreasing = TRUE)
  coords <- coords[taxa_order, ]
  estimate_from <- estimate_from[taxa_order]
  center_from <- center_from[taxa_order]

  # Overlap by the half of the size of the window
  overlap <- window_size %/% 2

  # Obtain number of windows by integer division
  n_last_window <- n_rows %/% window_size

  # Add one more window if the leftover exceeds the overlap, or if there
  # isn't even one full window yet. This is safe even when the leftover is
  # *smaller* than the overlap (no extra window added): see
  # .compute_window_bounds()'s doc for the verification that the last
  # window always absorbs it regardless of size.
  if ((n_rows %% window_size > overlap) | (n_last_window == 0)) {
    n_last_window <- n_last_window + 1
  }

  # Fallback for windows too thin to estimate from on their own. Note this is
  # also exactly what a population smaller than `window_size` gets anyway: the
  # arithmetic above yields a single window spanning every row, so a small
  # population is pooled without needing a special case.
  pooled <- .slide_z_params(
    coords,
    empirical_null_distribution,
    estimate_from,
    center_on,
    center_from
  )

  slide_z_all <- c()
  n_fallback <- 0L

  # Loop through windows
  for (n_window in 1:n_last_window) {
    bounds <- .compute_window_bounds(
      n_window = n_window,
      n_last_window = n_last_window,
      window_size = window_size,
      overlap = overlap,
      n_rows = n_rows
    )

    # taxa_slice = slice window +- overlap, use it to estimate mean and sd
    slice_idx <- bounds$window_start:bounds$window_end
    taxa_slice <- coords[slice_idx, ]

    params <- .slide_z_params(
      taxa_slice,
      empirical_null_distribution,
      estimate_from[slice_idx],
      center_on,
      center_from[slice_idx]
    )
    fell_back <- FALSE
    if (is.na(params$scale) || params$n < min_estimable) {
      params$scale <- pooled$scale
      fell_back <- TRUE
    }
    if (is.na(params$center) || params$n_center < min_estimable) {
      params$center <- pooled$center
      fell_back <- TRUE
    }
    # With a single window the pooled set and the window's own set are the same
    # rows, so this is not a fallback and warning about one would be noise --
    # it would fire for every population smaller than `window_size`.
    if (fell_back && n_last_window > 1) {
      n_fallback <- n_fallback + 1L
    }

    slide_z_slice <- (taxa_slice$obs_change - params$center) / params$scale

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
  list(
    slide_z = slide_z_all[order(taxa_order)],
    n_fallback = n_fallback
  )
}

#' Which Taxa May Contribute to a Slide-Z Null
#'
#' A taxon is usable as a reference observation when the *pair* that supplies the
#' reference is admissible: `null_in_cone` in empirical-null mode (where the reference
#' is `null_change`, from Ig-.1 vs Ig-.2), `obs_in_cone` otherwise. The two can differ
#' — a taxon whose Ig+ vs Ig- pair left the admissible cone may still have a perfectly
#' good negative-vs-negative pair, and vice versa.
#'
#' `NA` counts as usable: the `"log_ratio"` change axis has no cone, so
#' [get_ma_coordinates()] leaves both flags `NA` and every taxon is admissible.
#'
#' Taxa flagged by `was_imputed` are usable regardless of the cone, because inside that
#' population the flag no longer reports measurement quality -- see the note in
#' [compute_slide_z()]'s Details.
#'
#' This mask covers admissibility only. The *other* exclusion -- pairs with both
#' fractions really zero, which carry no information at all -- is not handled here:
#' [get_ma_coordinates()] already marks those `NA`, and `.slide_z_params()` drops `NA`
#' references on its own. Do not re-encode it as a mask; `NA` also removes such a taxon
#' from the numerator, which a mask would not, and that is the correct behaviour for an
#' untestable observed pair.
#'
#' Which pair's flag applies is therefore a property of the *quantity* being
#' estimated, not of the call: under `center_on = "observed"` the scale still comes
#' from the null pair while the location comes from the observed one, so the two need
#' different masks. `.slide_z_cone_mask()` builds either; `.slide_z_estimate_mask()`
#' is the thin wrapper picking whichever pair supplies the scale.
#'
#' @param ma_coords MA-plot coordinates, optionally carrying `obs_in_cone`/`null_in_cone`.
#' @param pair Which pair's admissibility to read: `"observed"` or `"null"`.
#' @param empirical_null_distribution See [compute_slide_z()].
#' @param was_imputed Optional logical vector, `nrow(ma_coords)` long, flagging the
#'   imputed population; `NULL` (default) applies the cone to every row.
#'
#' @return Logical vector, `nrow(ma_coords)` long.
#' @noRd
.slide_z_cone_mask <- function(
  ma_coords,
  pair = c("observed", "null"),
  was_imputed = NULL
) {
  pair <- match.arg(pair)
  flag <- if (pair == "null") ma_coords$null_in_cone else ma_coords$obs_in_cone
  if (is.null(flag)) {
    return(rep(TRUE, nrow(ma_coords)))
  }
  mask <- ifelse(is.na(flag), TRUE, as.logical(flag))
  if (is.logical(was_imputed) && length(was_imputed) == length(mask)) {
    mask[was_imputed] <- TRUE
  }
  mask
}

#' Mask for Whichever Pair Supplies the Slide-Z Scale
#'
#' See `.slide_z_cone_mask()`, of which this is the wrapper that reads
#' `empirical_null_distribution` instead of a `pair`.
#'
#' @param ma_coords,empirical_null_distribution,was_imputed See `.slide_z_cone_mask()`.
#' @return Logical vector, `nrow(ma_coords)` long.
#' @noRd
.slide_z_estimate_mask <- function(
  ma_coords,
  empirical_null_distribution,
  was_imputed = NULL
) {
  .slide_z_cone_mask(
    ma_coords,
    pair = if (empirical_null_distribution) "null" else "observed",
    was_imputed = was_imputed
  )
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
#' `window_size` taxa, each window's own local location and scale used to
#' center/scale the Z-scores of the taxa in that window's non-overlapping
#' "core".
#'
#' @section Where the centre comes from:
#' Location and scale are estimated from **different** distributions, because they are
#' not properties of the same thing. `empirical_null_distribution` chooses the scale:
#' the sd of the Ig-.1 vs Ig-.2 pair measures replicate noise and carries no coating
#' signal, which is what a yardstick should be. `center_on` chooses the location.
#'
#' **Zero on the change axis is not an arbitrary point.** Write `p_i` for the fraction
#' of taxon `i`'s cells that are Ig-coated and `n_i` for its true abundance. The Ig+
#' pool is `n_i p_i`, the Ig- pool `n_i (1 - p_i)`, and with both fractions sequenced to
#' a common depth,
#'
#' \deqn{E[\mathrm{obs\_change}_i] = \mathrm{logit}_2(p_i) + \log_2\frac{\sum n(1-p)}{\sum n p}}
#'
#' The second term is one constant fixed by closure, and since `sum(n p) / sum(n)` is
#' `P`, the community's overall Ig+ frequency -- the very quantity the assay measures on
#' the pre-sort fraction -- that constant is exactly `-logit2(P)`. So
#' `E[obs_change_i] = logit2(p_i) - logit2(P)`: the raw log-ratio already sits on an
#' anchored scale where **0 means "coated at the community's own rate `P`"**. (Verified
#' numerically to 6 decimals; a corollary is that a community where *every* taxon is
#' coated at the same rate has `obs_change = 0` for every taxon, whatever that rate --
#' uniform coating is invisible to this design, and no choice of centre changes that.)
#'
#' * `"reference"` (default) -- the mean of whichever distribution supplies the scale,
#'   i.e. the mean of Ig-.1 vs Ig-.2 under the empirical null. At matched depth that is
#'   ~0 by construction, so this **preserves the anchor**: `Z_i` estimates
#'   `(logit2(p_i) - logit2(P))` over the local replicate sd. It answers "is this taxon
#'   coated at a higher rate than the community as a whole?"
#' * `"observed"` -- the **median** of `obs_change` over the window, which re-anchors on
#'   the typical taxon instead. It answers "is this taxon coated at a higher rate than
#'   the median taxon in its abundance band?" -- a legitimate but different question,
#'   and one that gives up the tie to `P`.
#'
#' On simulations with known `p_i`, `"reference"` recovers the sign of the truth
#' (whether taxon `i` is coated above `P`) in 100% / 97.2% / 83.2% of taxa across a
#' coated minority, a coated majority and an all-coated community, against 68.2% /
#' 65.2% / 83.5% for `"observed"`. Worse, when coating is confounded with abundance --
#' the top decile coated at 0.95, the rest at 0.05 -- windowing the median puts the
#' coated taxa in a window whose median *is* the coated level, collapsing
#' `Spearman(Z, p_i)` from 0.62 to 0.22. `"observed"` centring removes signal in exactly
#' the case where the signal is strongest.
#'
#' Prefer `"observed"` only when the tie to `P` is not wanted or not available, e.g. no
#' pre-sort fraction and no reason to trust that the two sorted fractions were sequenced
#' at comparable depth. Which brings up the one real hazard of `"reference"`:
#'
#' @section The anchor depends on matched depth:
#' `"reference"` subtracts the mean of the Ig-.1 vs Ig-.2 pair, which is ~0 only while
#' the fractions are at a common depth. It is not a general-purpose estimate of the
#' anchor: with `rarefy_by_sample = FALSE` in [group_sorted_samples()] it subtracts
#' `log2(S_neg1/S_neg2)` from a quantity carrying `log2(S_pos/S_neg1)`, two different
#' constants. On `ps_igseq` that translated whole samples by -1.77 to +1.20 Z units and
#' **inverted the conclusion**: sample_4 read 5 enriched / 29 depleted rarefied against
#' 44 / 5 unrarefied, sample_5 1 / 10 against 13 / 0. Keep `rarefy_by_sample = TRUE`
#' unless you have corrected the depth ratio yourself; the anchor argument above assumes
#' it throughout.
#'
#' Taxa flagged via `was_imputed` form a **separate population**, scored against a
#' single pooled mean/sd taken over all of them at once rather than mixed into the
#' sliding windows. Two separate reasons, both measured on `ps_igseq` under
#' `pseudo_count`: their null is 2.8x-6.9x wider than that of taxa with no imputed
#' zero and they barely overlap the rest in abundance, so windowing them *with* the
#' rest would inflate their Z-scores several-fold; and their abundance axis is largely
#' an artefact of the imputation -- 75-83% of them share an `obs_abundance` with
#' another taxon -- so windowing them *among themselves* estimates the scale from tied,
#' arbitrary slices and deflates it by up to 5x, which inflates the Z-scores just as
#' badly in the other direction. Pooling is what neither failure mode reaches.
#'
#' Independently of that split, `estimate_from` marks rows that are scored but may not
#' *inform* a window's mean/sd. This is how out-of-cone taxa on the
#' `"purity_corrected"` change axis are handled: their reference value is saturated
#' rather than measured, so it must not enter the estimate, but they still receive a
#' score on their window's scale. Unlike imputed taxa they are not segregated, because
#' measurement shows they sit at the same abundances as the rest and their null spread
#' is within noise of it.
#'
#' That cone exclusion is deliberately **not** applied inside the imputed population,
#' where the flag stops reporting measurement quality and starts reporting the
#' imputation pattern. A taxon zeroed in only one of the two negative fractions gets an
#' extreme ratio and falls outside the cone; one zeroed in neither stays inside. Masking
#' on the cone there therefore keeps only the taxa whose negative fractions were both
#' measured and discards every censored pair -- and censored pairs are the dominant mode
#' of genuine replicate noise in this assay, not an edge case (see
#' [get_ma_coordinates()]'s Uninformative pairs section). Measured on `ps_igseq` under
#' `pseudo_count`, masking halves the reference sd (pooled 2.2-3.6 versus masked
#' 1.0-1.8, ratio 0.47-0.64), which would roughly double every imputed taxon's `|Z|`.
#'
#' The pairs that genuinely carry no information -- *both* fractions zero, so the change
#' is one substituted value over itself -- are cut upstream instead, by
#' [get_ma_coordinates()], which marks them `NA` on the affected axis; `.slide_z_params()`
#' then drops them from the estimate, and from the numerator too when they are the
#' observed pair. That is a cleaner cut than the cone, because it targets exactly the
#' fabricated values without touching the censored ones. Before it existed those pairs
#' were reported as a change of exactly 0: on `ps_igseq` under `pseudo_count` they made
#' up 13-43% of this population's reference values and deflated its pooled sd by 7-31%,
#' and every taxon whose Ig+ and Ig-.1 were both zero received the same single
#' manufactured Z-score as all the others.
#'
#' Even so, neither remaining subset of this population is clean -- the censored ones
#' depend on the pseudocount for their magnitude, the out-of-cone ones are saturated --
#' so pooling over all of them is the least-wrong estimate, not a good one. Their `|Z|`
#' remains a floor rather than a test, exactly as for the saturated taxa described in
#' [get_slide_z()].
#'
#' @param ma_coords MA-plot coordinates for one biological sample, as
#'   returned by [get_ma_coordinates()] (`taxon_id`, `obs_abundance`,
#'   `obs_change`, and, if `empirical_null_distribution`, `null_abundance`/
#'   `null_change`).
#' @param was_imputed Optional logical vector, same length/row order as
#'   `ma_coords`, flagging which taxa had an imputed zero and should be
#'   scored as a population apart.
#' @param window_size See [get_slide_z()]. Default `50`.
#' @param empirical_null_distribution See [get_slide_z()]. Default `TRUE`. Chooses the
#'   *scale* only; see the "Where the centre comes from" section and `center_on`.
#' @param center_on Which distribution supplies each window's location:
#'   `"reference"` (default) the mean of whichever distribution
#'   `empirical_null_distribution` picked for the scale, `"observed"` the median of
#'   `obs_change`. See the "Where the centre comes from" section -- these are different
#'   statistics answering different questions, not a tuning knob.
#' @param estimate_from Optional logical vector, same length/row order as `ma_coords`,
#'   marking which taxa may contribute to a window's estimate. `NULL` (default) derives
#'   it per quantity from `ma_coords`'s `obs_in_cone`/`null_in_cone` columns -- the
#'   scale from whichever pair supplies it, the location always from `obs_in_cone`
#'   under `center_on = "observed"` -- treating `NA` (the `"log_ratio"` axis, which has
#'   no cone) as usable. Supplying it explicitly overrides **both**.
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
#'   fraction_ids = c("pos", "neg1", "neg2")
#' )
#' ma_coords <- get_ma_coordinates(
#'   sorted_sample_df = grouped[["sample_1"]],
#'   positive_fraction_name = "pos",
#'   first_negative_fraction_name = "neg1",
#'   second_negative_fraction_name = "neg2"
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
  empirical_null_distribution = TRUE, # picks the SCALE
  center_on = c("reference", "observed"), # picks the LOCATION
  estimate_from = NULL # which rows may inform a window's estimate
) {
  center_on <- match.arg(center_on)
  has_imputed <- is.logical(was_imputed) &&
    length(was_imputed) == nrow(ma_coords) &&
    any(was_imputed)

  if (is.null(estimate_from)) {
    estimate_from <- .slide_z_estimate_mask(
      ma_coords,
      empirical_null_distribution,
      was_imputed = if (has_imputed) was_imputed else NULL
    )
    # The location is read off `obs_change`, so it is the OBSERVED pair's
    # admissibility that governs it -- not the null pair's, which may differ (see
    # `.slide_z_cone_mask()`). A saturated observed value must not shift the median
    # even when its negative-vs-negative pair is perfectly good.
    center_from <- if (center_on == "observed") {
      .slide_z_cone_mask(
        ma_coords,
        pair = "observed",
        was_imputed = if (has_imputed) was_imputed else NULL
      )
    } else {
      estimate_from
    }
  } else {
    center_from <- estimate_from
  }

  if (!has_imputed) {
    result <- .slide_z_windowed(
      coords = ma_coords,
      estimate_from = estimate_from,
      window_size = window_size,
      empirical_null_distribution = empirical_null_distribution,
      center_on = center_on,
      center_from = center_from
    )
    if (nrow(ma_coords) == 0) {
      warning(paste0(
        "No taxa in MA coordinates, cannot compute slide z score\n"
      ))
      return(NULL)
    }
    .warn_slide_z_fallback(result$n_fallback)
    return(result$slide_z)
  }

  # Imputed taxa are scored as their own population, not merged into the windows
  # of the rest. Measured on ps_igseq under pseudo_count imputation, their null
  # distribution is 2.8x-6.9x wider than that of taxa with no imputed zero, and
  # the two barely overlap in abundance -- in some samples the entire bottom 60%
  # of the abundance range is imputed, so there is no non-imputed taxon nearby to
  # estimate a null from in the first place. Scoring them against the
  # non-imputed null would inflate their Z-scores several-fold and call almost
  # every imputed taxon significant.
  #
  # Within that population the estimate is POOLED, not windowed, because the
  # abundance axis there is largely an artefact of the imputation and so cannot
  # order taxa into meaningful neighbourhoods: a zero treatment substitutes the
  # same value for every zero, which on ps_igseq left 75-83% of imputed taxa
  # sharing an `obs_abundance` with another and collapsed 943 taxa onto as few as
  # 174 distinct `null_change` values. Windows cut from that tied ordering are
  # arbitrary slices rather than local neighbourhoods, and their sd underestimates
  # the population spread by up to 5x (pooled 3.45 vs 0.71 in the thinnest window),
  # which divides the Z-scores up rather than localizing them -- it doubled the
  # per-sample sd and turned 320 significant imputed taxa into 1527.
  # Verified 2026-08-06 against f86e051. Don't window this population without
  # first showing its abundance axis carries real information.
  clean_coords <- ma_coords[!was_imputed, ]
  imputed_coords <- ma_coords[was_imputed, ]

  if (nrow(clean_coords) == 0) {
    warning(paste0("No taxa in MA coordinates, cannot compute slide z score\n"))
    return(NULL)
  }

  clean_result <- .slide_z_windowed(
    coords = clean_coords,
    estimate_from = estimate_from[!was_imputed],
    window_size = window_size,
    empirical_null_distribution = empirical_null_distribution,
    center_on = center_on,
    center_from = center_from[!was_imputed]
  )
  .warn_slide_z_fallback(clean_result$n_fallback)

  slide_z_merged <- rep(NA_real_, length(was_imputed))
  slide_z_merged[!was_imputed] <- clean_result$slide_z
  # This population gets its own location too, not the clean population's: its
  # `obs_change` values are shifted by the zero treatment (a substituted count on
  # one side of a pair moves the ratio), so borrowing the clean centre would score
  # that shift as coating.
  slide_z_merged[was_imputed] <- .slide_z_center_scale(
    imputed_coords,
    empirical_null_distribution,
    estimate_from = estimate_from[was_imputed],
    center_on = center_on,
    center_from = center_from[was_imputed]
  )
  slide_z_merged
}

#' Warn About Windows That Could Not Estimate Their Own Null
#'
#' @param n_fallback How many windows fell back to the pooled estimate.
#' @return `invisible(NULL)`, called for the `warning()`.
#' @noRd
.warn_slide_z_fallback <- function(n_fallback) {
  if (n_fallback > 0) {
    warning(
      n_fallback,
      " sliding window(s)",
      " had too few usable taxa to estimate a local null and fell back to the ",
      "pooled estimate. Consider a larger `window_size`, or check how many taxa ",
      "the zero treatment and the admissible cone are excluding.\n"
    )
  }
  invisible(NULL)
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
#' @param center_on Where the fitted envelope sits on the change axis:
#'   `"observed"` (default) translates it so its centre is the median observed
#'   change, `"reference"` leaves it centred on the fitted cloud's own mean (the
#'   behaviour before 2026-08-07). The same choice [compute_slide_z()] makes for the
#'   Z-score, applied here so the two significance calls do not disagree — see its
#'   "Where the centre comes from" section. Only the change axis is moved; the
#'   abundance axis is left on the fitted cloud, which is a separate mismatch (an
#'   observed point can lie outside the envelope on `A` alone, measured at 0-4 taxa
#'   per sample on `ps_igseq`).
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
#'       missing `obs_abundance`/`obs_change`. That last case now also covers taxa
#'       whose observed pair was zero in both fractions, which [get_ma_coordinates()]
#'       marks `NA`: they are reported as untested rather than as `"ns"`, which is
#'       what they used to be — a fabricated negative. `NULL` instead if there
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
#'   fraction_ids = c("pos", "neg1", "neg2")
#' )
#' ma_coords <- get_ma_coordinates(
#'   sorted_sample_df = grouped[["sample_1"]],
#'   positive_fraction_name = "pos",
#'   first_negative_fraction_name = "neg1",
#'   second_negative_fraction_name = "neg2"
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
    center_on = c("reference", "observed"),
    confidence_levels # a vector of confidence levels
  ) {
    center_on <- match.arg(center_on)
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

    # Taxa the un-mixing had to clamp must not shape the envelope: their change
    # value is saturated rather than measured. They are still *tested* against it
    # below -- their abundance is real and their saturated change is genuinely
    # extreme -- so this restricts the fit only, matching how compute_slide_z()
    # keeps them out of a window's mean/sd while still scoring them. All-`NA`
    # under the "log_ratio" axis, which has no cone.
    estimable_taxa <- .slide_z_estimate_mask(
      sorted_sample_df,
      empirical_null_distribution
    )

    if (empirical_null_distribution) {
      valid_taxa_null <- !is.na(sorted_sample_df$null_abundance) &
        !is.na(sorted_sample_df$null_change) &
        !sorted_sample_df$taxon_id %in% imputed_taxa
      # based on neg1 vs neg2, construct Ellipses
      fit_taxa <- valid_taxa_null & estimable_taxa
      abund_coords <- sorted_sample_df$null_abundance[fit_taxa]
      change_coords <- sorted_sample_df$null_change[fit_taxa]
    } else {
      # based on pos vs neg1, construct Ellipses
      fit_taxa <- valid_taxa_obs & estimable_taxa
      abund_coords <- sorted_sample_df$obs_abundance[fit_taxa]
      change_coords <- sorted_sample_df$obs_change[fit_taxa]
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

    # `car::dataEllipse()` centres every level on the fitted cloud's own mean. Under
    # `center_on = "observed"` the envelope keeps that cloud's shape, spread and tilt
    # -- which is what the Ig-.1 vs Ig-.2 pair is there to measure -- but slides along
    # the change axis to sit on the observed cloud's median, the same location
    # compute_slide_z() subtracts. Without this the ellipse and the Z-score would
    # disagree about which taxa are significant whenever the two clouds do not
    # coincide, which on ps_igseq is always.
    if (center_on == "observed") {
      center_taxa <- valid_taxa_obs &
        .slide_z_cone_mask(sorted_sample_df, pair = "observed")
      observed_center <- stats::median(
        sorted_sample_df$obs_change[center_taxa],
        na.rm = TRUE
      )
      change_shift <- observed_center - mean(change_coords)
      if (is.finite(change_shift)) {
        ellipse_list <- lapply(ellipse_list, function(ellipse) {
          ellipse[, 2] <- ellipse[, 2] + change_shift
          ellipse
        })
      }
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
