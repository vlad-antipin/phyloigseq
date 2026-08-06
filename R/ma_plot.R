# Implements the idea of MA plot (from microarray analysis) to abundance analysis

#' Compute MA-Plot Coordinates Between IgSeq Fractions
#'
#' Computes MA-plot-style coordinates (borrowed from microarray analysis: `M`
#' = change, `A` = log-abundance) between a positive IgSeq fraction and a
#' negative fraction, per taxon. If a second negative fraction is supplied,
#' the same coordinates are also computed between the two negative
#' fractions, to serve as an empirical null distribution (see
#' [get_slide_z()]).
#'
#' The change (`M`) axis is selectable via `change_transform`, because the default
#' `log2(pos/neg)` is a base-2 logit up to an additive constant -- it equals
#' `"prob_ratio"` minus `log2(P/(1-P))`. Any *affine* reparameterization of it is
#' annihilated by [compute_slide_z()]'s centering and scaling, so only a non-affine
#' change axis alters a sliding Z-score:
#'
#' * `"log_ratio"` (default) -- `log2(a) - log2(b)`, the classical MA-plot `M`, on the
#'   fraction abundances exactly as supplied.
#' * `"purity_corrected"` -- `log2(a_t/b_t)` on the populations recovered by un-mixing
#'   the two fractions with their measured Ig+ frequencies (see [compute_ig_score()]'s
#'   Details for the model). Requires `pos_ig_freq` and `neg_ig_freq`, but *not* a
#'   pre-sort fraction: the prior term `log2(P/(1-P))` is an additive per-sample
#'   constant, so it is omitted, leaving an axis centred at 0 under the null and
#'   directly overlayable on `"log_ratio"`. Nonlinear in `pos/neg`, with local slope
#'   `1/(pos_ig_freq - neg_ig_freq)` at the null, so it stretches the tails while
#'   leaving taxa near the null essentially where they were.
#'
#' The `A` axis is `log10(a) + log10(b)` for both: it is a precision proxy (how much
#' evidence a taxon carries) and should not follow the change transform.
#'
#' Both the observed and the empirical-null pair go through the same map. Under the
#' null hypothesis both are two samples of the same underlying composition at matched
#' depth, so their ratios are identically distributed and identically transformed
#' values are too -- which is what makes the null a valid reference for the observed
#' axis, and why the second negative fraction's own Ig+ frequency is never needed.
#'
#' That equivalence is what the two axes trade against each other. `"log_ratio"` reads
#' the abundances as given, so both pairs carry the same (zero, under matched depth)
#' offset and [compute_slide_z()] cancels it exactly. `"purity_corrected"` must convert
#' each fraction to a composition before un-mixing, which normalizes the two pairs
#' independently; because the un-mixing is non-linear the resulting per-pair offsets do
#' not cancel, so that axis carries a residual shift `"log_ratio"` does not. Do not
#' "harmonize" the two by normalizing `"log_ratio"` as well -- that reintroduces the
#' same non-cancelling shift on the default axis, and it is not visible in the scores,
#' only in where the significance threshold lands.
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
#' @param change_transform Which change (`M`) axis to compute: `"log_ratio"`
#'   (default) or `"purity_corrected"`. See Details.
#' @param pos_ig_freq,neg_ig_freq Single probabilities in `[0, 1]`, the Ig+ frequencies
#'   measured inside the positive and first negative fractions (`p` and `q`). Required
#'   by `change_transform = "purity_corrected"`, ignored otherwise.
#' @param cone_policy How to treat taxa outside the admissible cone
#'   `pos/neg \in [(1-q')/(1-p'), p'/q']` under `"purity_corrected"`, where the un-mixing
#'   returns a negative population. `"regularize"` (default) offsets both pools by
#'   `pool_prior`, so such taxa saturate at a finite, correctly-ordered value.
#'   `"na"` leaves them `NA`; this is for diagnostics only -- it discards a median of
#'   41.5% of taxa per sample on the SMILE cohort, including the most strongly coated
#'   ones.
#' @param pool_prior Passed to the un-mixing under `"purity_corrected"`; `NULL`
#'   (default) derives one read's worth of composition from the fraction depth. See
#'   [compute_ig_score()].
#'
#' @return A data frame with one row per taxon: `taxon_id`, `sample_id`, the
#'   original fraction abundances (`pos`, `neg1`, `neg2`; `neg2` is `NA` when
#'   `second_negative_fraction_name` is `NULL`), the observed coordinates
#'   (`obs_abundance`, `obs_change`, positive vs. first negative fraction),
#'   and, when a second negative fraction is supplied, the empirical-null
#'   coordinates (`null_abundance`, `null_change`, first vs. second negative
#'   fraction); otherwise `NA`. Plus `change_transform`, recording which axis was
#'   computed, and `obs_in_cone`/`null_in_cone`, logical flags marking taxa the
#'   un-mixing had to clamp (`NA` throughout under `"log_ratio"`, which has no cone).
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
#' head(ma_coords)
#'
#' @export
get_ma_coordinates <- function(
  sorted_sample_df, # dataframe from the result of group_sorted_samples()
  positive_fraction_name,
  first_negative_fraction_name, # 9/10 of the whole negative fraction for IgSeq
  second_negative_fraction_name = NULL, # 1/10
  change_transform = c("log_ratio", "purity_corrected"),
  pos_ig_freq = NULL, # "p", required by "purity_corrected"
  neg_ig_freq = NULL, # "q", -//-
  cone_policy = c("regularize", "na"),
  pool_prior = NULL
) {
  change_transform <- match.arg(change_transform)
  cone_policy <- match.arg(cone_policy)

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

  if (change_transform == "purity_corrected") {
    if (
      !isTRUE(is.finite(pos_ig_freq)) || !isTRUE(is.finite(neg_ig_freq))
    ) {
      warning(
        "`change_transform = \"purity_corrected\"` needs both `pos_ig_freq` and ",
        "`neg_ig_freq` (the Ig+ frequencies measured inside the sorted ",
        "fractions); the change axis will be all-NA.\n"
      )
    } else if (
      cone_policy == "regularize" &&
        is.null(pool_prior) &&
        is.na(.resolve_pool_prior(NULL, pos, neg1, pos_ig_freq, neg_ig_freq))
    ) {
      # Compositional input carries no read scale, so the default prior cannot be
      # derived. Falling back to "na" is honest: it reports the out-of-cone taxa
      # as missing rather than saturating them at an arbitrary offset.
      warning(
        "Fraction abundances are already compositional, so the read scale needed ",
        "for the default `pool_prior` is unavailable; falling back to ",
        "`cone_policy = \"na\"`, which drops out-of-cone taxa from the change ",
        "axis. Pass counts (`transform_by_sample = \"identity\"`), or set ",
        "`pool_prior` explicitly.\n"
      )
      cone_policy <- "na"
    }
  }

  # Compute MA coordinates (M = change, A = log-abundance), see MA plot.
  # Alternatives to log2/log10 (e.g. asin(sqrt(x/sum(x))), sqrt(x/sum(x)))
  # are worth revisiting if 0-handling becomes an issue, since they don't
  # blow up at x = 0 the way a log transform does.
  transform_A <- function(x) log10(x)

  # A stays log10(a) + log10(b) whatever the change axis: it is a precision proxy
  # (how much evidence this taxon carries) and should not follow the Y transform.
  #
  # The change axis is pluggable because log2(a/b) is a base-2 logit up to an
  # additive constant, and any *affine* reparameterization of it is annihilated by
  # compute_slide_z()'s centering and scaling. Only a non-affine change axis --
  # i.e. "purity_corrected" -- actually alters a sliding Z-score.
  compute_change <- switch(
    change_transform,
    log_ratio = function(a, b, a_rel, b_rel) {
      # Raw values, NOT a_rel/b_rel: see compute_ma_pair() for why per-pair
      # normalization is not neutral here.
      list(change = log2(a) - log2(b), in_cone = rep(NA, length(a)))
    },
    purity_corrected = function(a, b, a_rel, b_rel) {
      prior <- if (cone_policy == "regularize") {
        .resolve_pool_prior(pool_prior, a, b, pos_ig_freq, neg_ig_freq)
      } else {
        0
      }
      if (is.na(prior)) prior <- 0
      pools <- .unmix_ig_pools(
        a_rel,
        b_rel,
        pos_ig_freq,
        neg_ig_freq,
        pool_prior = prior
      )
      # log2(a_t/b_t): the purity-corrected log2 fold change. This is
      # `purity_corrected_prob_ratio` minus log2(P/(1-P)); the prior term is
      # dropped because it is an additive per-sample constant that centering
      # removes anyway, which keeps this axis independent of the pre-sort
      # fraction and centred at 0 under the null.
      list(
        change = log2(pools$ig_pos) - log2(pools$ig_neg),
        in_cone = pools$in_cone
      )
    }
  )

  # Obs - observed, null - empirical null distribution (control vs control).
  # Both pairs go through the *same* map: under the null hypothesis they are two
  # samples of the same underlying composition at matched depth, so their ratios
  # are identically distributed and identically transformed values are too. That
  # is what makes the null usable as a reference for the observed axis, and it is
  # why the second negative fraction's own Ig+ frequency is not needed here.
  compute_ma_pair <- function(a, b) {
    # A stays on the raw count scale: more reads means more precision, which is
    # what this axis stands in for, and it is only ever used for ranking taxa into
    # windows and for display.
    abundance <- transform_A(a) + transform_A(b)
    # Relative abundances are offered to the change transforms, but ONLY
    # "purity_corrected" takes them, because the un-mixing needs compositions
    # outright (the two fractions enter the 2x2 system with different
    # coefficients, so their totals do not cancel).
    #
    # "log_ratio" must NOT use them, however tempting the symmetry. Dividing each
    # member of a pair by its own total adds a per-pair constant: log2(sum(b)/sum(a))
    # for this pair. compute_slide_z() removes an additive constant only when the
    # observed and the null pair carry the SAME one -- and they do not. The observed
    # pair (pos, neg1) picks up log2(S_neg1/S_pos) while the null pair (neg1, neg2)
    # picks up log2(S_neg2/S_neg1), so under the default
    # `empirical_null_distribution = TRUE`, where the null supplies the center, the
    # residual log2(S_neg1^2 / (S_pos * S_neg2)) survives and translates every
    # Z-score in the sample by that amount over the local sd.
    #
    # These sums are not equal even under `rarefy_by_sample = TRUE`: rarefaction
    # equalizes depth across ALL taxa, but impute_zeros() then drops or fills taxa,
    # and get_ma_coordinates() only ever sees the surviving subset. Measured on
    # ps_igseq with the default `zero_treatment = "no_zero"`, this shifted whole
    # samples by up to 1.5 Z units (worst in the sample where only 16 taxa
    # survived) and flipped 3.9% of |Z| > 1.96 calls. Verified 2026-08-06 against
    # f86e051; raw values reproduce it bit-for-bit.
    #
    # The pipeline's matched-depth assumption (see rarefy_abundances()) is what
    # makes the raw log-ratio the right scale: both pairs then carry the constant 0
    # and cancel exactly. Note "purity_corrected" cannot have this property, since
    # it needs compositions; its two pairs are normalized independently and the
    # un-mixing is non-linear, so its residual offset does not cancel either. That
    # is a known limitation of that axis, not something raw values could fix.
    res <- compute_change(
      a,
      b,
      a_rel = a / sum(a, na.rm = TRUE),
      b_rel = b / sum(b, na.rm = TRUE)
    )
    change <- res$change
    abundance[is.nan(abundance) | is.infinite(abundance)] <- NA
    change[is.nan(change) | is.infinite(change)] <- NA
    list(abundance = abundance, change = change, in_cone = res$in_cone)
  }

  obs <- compute_ma_pair(pos, neg1)
  obs_abundance <- obs$abundance
  obs_change <- obs$change
  obs_in_cone <- obs$in_cone

  if (empirical_null) {
    null <- compute_ma_pair(neg1, neg2)
    null_abundance <- null$abundance
    null_change <- null$change
    null_in_cone <- null$in_cone
  } else {
    null_abundance <- rep(NA, nrow(sorted_sample_df))
    null_change <- rep(NA, nrow(sorted_sample_df))
    null_in_cone <- rep(NA, nrow(sorted_sample_df))
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
      null_change,
      change_transform = change_transform,
      obs_in_cone,
      null_in_cone
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
#'   Comparing these on `change_transform = "purity_corrected"` is much less
#'   informative than it looks: the un-mixing regularizes both pools by
#'   `pool_prior`, one read's worth of composition by default, and imputed
#'   abundances are typically *below* that (3.6x below under `"pseudo_count"`,
#'   347x under `"bayesian_inference"`, on `ps_igseq`), so the prior rather than
#'   the imputed value sets the change of every imputed taxon. The treatments
#'   then agree to within ~2% of the spread and render as one band. That is the
#'   regularization behaving as intended -- an imputed abundance below one read
#'   carries less than one read of evidence -- but it means an apparent
#'   insensitivity to zero handling on that axis is not evidence of one. Compare
#'   on `"log_ratio"`, where the same taxa separate by 29x more.
#' @param change_transform Character vector of change (`M`) axes to compute --
#'   `"log_ratio"` and/or `"purity_corrected"`, see [get_ma_coordinates()]. Like
#'   `zero_treatments`, several may be given: each contributes one
#'   `change_transform` level to `plot_data`, so a caller can switch axes by
#'   filtering rather than recomputing. Default `"log_ratio"`.
#' @param pos_ig_freq,neg_ig_freq,cone_policy,pool_prior Passed to
#'   [get_ma_coordinates()]; only read when `"purity_corrected"` is among
#'   `change_transform`.
#'
#' @return A list with:
#'   \describe{
#'     \item{`sample_id`}{The biological sample id
#'       (`sorted_sample_df$sample_id[1]`).}
#'     \item{`nb_zero_taxa`}{Number of taxa with a zero abundance in at
#'       least one compared fraction, before imputation.}
#'     \item{`plot_data`}{A long-format data frame (`M`, `A`, `comparison`,
#'       `taxon_id`, `in_cone`, `imputed`, `zero_treatment`, `change_transform`)
#'       ready for [plot_ma()]. `imputed` is per-row, i.e. per `zero_treatment`:
#'       whether *that* treatment derived this taxon's coordinates from an
#'       imputed zero.}
#'     \item{`imputed_taxa`}{The union, across all `zero_treatments`, of
#'       `taxon_id`s imputed under that treatment (see [impute_zeros()]). Kept
#'       for callers that want one set for the whole sample; use `plot_data$imputed`
#'       when the treatment matters.}
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
#' ma_plot_data <- get_ma_plot_data(
#'   sorted_sample_df = grouped[["sample_1"]],
#'   positive_fraction_name = "pos",
#'   first_negative_fraction_name = "neg1",
#'   second_negative_fraction_name = "neg2",
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
  zero_treatments = c("keep_zeros"),
  # Accepts several, like `zero_treatments`, so one call yields every axis a
  # caller might want to switch between without recomputing.
  change_transform = c("log_ratio"),
  pos_ig_freq = NULL,
  neg_ig_freq = NULL,
  cone_policy = c("regularize", "na"),
  pool_prior = NULL
) {
  cone_policy <- match.arg(cone_policy)
  change_transform <- unique(match.arg(
    change_transform,
    choices = c("log_ratio", "purity_corrected"),
    several.ok = TRUE
  ))

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

    # Zero handling depends only on the treatment, so impute once and reuse it
    # across change axes -- the axes differ purely in how M is derived.
    for (this_transform in change_transform) {
      ma_coords <- get_ma_coordinates(
        sorted_sample_df = sorted_sample_df_imputed,
        positive_fraction_name = positive_fraction_name,
        first_negative_fraction_name = first_negative_fraction_name,
        second_negative_fraction_name = second_negative_fraction_name,
        change_transform = this_transform,
        pos_ig_freq = pos_ig_freq,
        neg_ig_freq = neg_ig_freq,
        cone_policy = cone_policy,
        pool_prior = pool_prior
      )

      # Per-treatment, NOT the union below: whether a taxon's coordinates are
      # model-derived is a property of the treatment that produced them, and it
      # is the only thing that distinguishes the facets. Marking a taxon imputed
      # everywhere just because one treatment imputed it is what made every
      # facet render the same set of points the same way.
      imputed <- ma_coords$taxon_id %in% zero_imputation_result$imputed_taxa

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
          taxon_id = ma_coords$taxon_id,
          in_cone = ma_coords$obs_in_cone,
          imputed = imputed
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
            taxon_id = ma_coords$taxon_id,
            in_cone = ma_coords$null_in_cone,
            imputed = imputed
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
            ),
            change_transform = rep(this_transform, nrow(ma_coords_long))
          )
        )
      )
    }
  }

  plot_data$zero_treatment <- factor(
    plot_data$zero_treatment,
    levels = gsub("_", " ", zero_treatments)
  )
  plot_data$change_transform <- factor(
    plot_data$change_transform,
    levels = change_transform
  )

  return(list(
    sample_id = sample_id,
    nb_zero_taxa = nb_zero_taxa,
    plot_data = plot_data,
    imputed_taxa = imputed_taxa
  ))
}


#' Plot an IgSeq MA Plot
#'
#' Draws the MA-plot coordinates produced by [get_ma_plot_data()] for one
#' biological sample: log-abundance (`A`) on the x-axis, change (`M`) on
#' the y-axis. Taxa whose coordinates come from an imputed zero are drawn
#' separately, jittered just past the plot's x-range rather than placed at their
#' imputed/undefined `A` value.
#'
#' Which taxa those are is read per `zero_treatment`, from `plot_data$imputed`,
#' not from the `imputed_taxa` union: a treatment that imputes nothing (`"no_zero"`,
#' which drops the taxa instead) must not inherit another treatment's imputed set,
#' or every facet renders the same points the same way and the comparison the
#' facets exist for shows nothing.
#'
#' Under `change_transform = "purity_corrected"`, taxa the un-mixing had to clamp
#' (`in_cone = FALSE`) are drawn **in place but hollow**. Unlike imputed taxa, whose
#' `A` is an artefact of the pseudocount, these have a perfectly real `A` -- only
#' their `M` is saturated -- so moving them off-axis would discard true information.
#' They trace a rising envelope rather than a flat clipped line, because the
#' saturation ceiling grows with abundance; a thick band riding that envelope means
#' `pos_ig_freq` and `neg_ig_freq` are too close for the correction to say much.
#'
#' @param ma_plot_data A list as returned by [get_ma_plot_data()].
#' @param ellipses If `TRUE`, overlay a 95% normal-theory confidence ellipse
#'   per `comparison` (only applies when `type = "facet"`). Default `FALSE`.
#' @param type One of `"facet"` (default; facet by `zero_treatment`, color
#'   by `comparison`) or `"superposed"` (facet by `comparison`, color by
#'   `zero_treatment`, comparing zero-handling strategies within the same
#'   panel).
#' @param change_transform Which change axis to draw when `ma_plot_data` carries
#'   more than one. `NULL` (default) draws the only one present, or `"log_ratio"`
#'   when several are. Errors if the requested axis was never computed.
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
#'   fraction_ids = c("pos", "neg1", "neg2")
#' )
#' ma_plot_data <- get_ma_plot_data(
#'   sorted_sample_df = grouped[["sample_1"]],
#'   positive_fraction_name = "pos",
#'   first_negative_fraction_name = "neg1",
#'   second_negative_fraction_name = "neg2",
#'   zero_treatments = c("keep_zeros", "pseudo_count")
#' )
#' plot_ma(ma_plot_data)
#' plot_ma(ma_plot_data, type = "superposed")
#'
#' @export
plot_ma <-
  function(
    ma_plot_data,
    ellipses = FALSE,
    type = c("facet", "superposed"),
    change_transform = NULL
  ) {
    type <- match.arg(type)

    plot_data <- ma_plot_data$plot_data
    # `get_ma_plot_data()` may carry several change axes at once so callers can
    # switch between them without recomputing; draw exactly one of them.
    change_transform <- .pick_change_transform(
      change_transform,
      plot_data$change_transform
    )
    if (!is.null(change_transform) && !is.null(plot_data$change_transform)) {
      plot_data <- plot_data[
        as.character(plot_data$change_transform) == change_transform,
      ]
    }

    # `imputed` is per zero-treatment; the union in `imputed_taxa` is the
    # fallback for data built before that column existed.
    is_imputed <- if (is.null(plot_data$imputed)) {
      plot_data$taxon_id %in% ma_plot_data$imputed_taxa
    } else {
      !is.na(plot_data$imputed) & plot_data$imputed
    }
    # Out-of-cone taxa under the purity-corrected axis. Their M is saturated but
    # their A is real, unlike imputed taxa whose A comes from the pseudocount --
    # so they are drawn in place, hollow, rather than banished to the jitter band.
    is_saturated <- !is_imputed &
      !is.null(plot_data$in_cone) &
      !is.na(plot_data$in_cone) &
      !plot_data$in_cone

    ma_non_imputed <- plot_data[!is_imputed, ]
    ma_clean <- plot_data[!is_imputed & !is_saturated, ]
    ma_saturated <- plot_data[is_saturated, ]
    ma_imputed <- plot_data[is_imputed, ]

    jitter <- .jitter_offset(ma_non_imputed$A)
    jitter_width <- jitter$width
    jitter_x <- jitter$x

    plt <- ggplot(ma_clean, aes(x = A, y = M))

    if (type == "facet") {
      plt <- plt +
        geom_point(aes(color = comparison), alpha = 0.8) +
        geom_point(
          data = ma_saturated,
          aes(color = comparison),
          alpha = 0.8,
          shape = 1
        ) +
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
        geom_point(
          data = ma_saturated,
          aes(color = zero_treatment, size = zero_treatment),
          alpha = 0.5,
          shape = 1
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
        facet_wrap(. ~ comparison, scales = "fixed") +
        # Earliest treatment largest, so each later one draws smaller on top of
        # it instead of hiding it: the default ordinal scale runs the other way,
        # which let the last treatment cover every earlier one outright.
        scale_size_ordinal(range = c(3, 1))
    }

    plt <- plt +
      labs(
        title = paste0("IgSeq MA plot for sample ", ma_plot_data$sample_id),
        subtitle = paste0(
          ma_plot_data$nb_zero_taxa,
          " taxa with zero abundances in at least one fraction",
          if (nrow(ma_saturated) > 0) {
            paste0(
              "; ",
              nrow(ma_saturated),
              " point(s) outside the admissible cone, drawn hollow"
            )
          }
        ),
        x = latex2exp::TeX(
          "Log-Abundance: $\\log_{10}\\left(\\fraction_{1}\\cdot\\fraction_{2}\\right)$"
        ),
        y = .change_axis_label(change_transform)
      ) +
      scale_x_continuous(breaks = .nonnegative_abundance_breaks) +
      .plot_title_theme() +
      scale_color_manual(values = c("darkgray", ggsci::pal_npg()(4)))
    #ggsci::scale_color_npg()

    return(plt)
  }


#' Resolve Which Change Axis to Draw
#'
#' Shared by [plot_ma()] and [plot_slide_z()]: both accept data that may carry more than
#' one change axis and must draw exactly one. `NULL` means "pick for me": the only axis
#' present, or `"log_ratio"` when several are.
#'
#' @param requested The caller's `change_transform`, or `NULL`.
#' @param available The `change_transform` column of the data (may be `NULL` for data
#'   produced before the column existed).
#'
#' @return A single change-transform name, or `NULL` when none is recorded in the data.
#' @noRd
.pick_change_transform <- function(requested, available) {
  present <- unique(as.character(available[!is.na(available)]))
  if (length(present) == 0) {
    # Data predates the column; nothing to filter on.
    return(if (is.null(requested)) NULL else requested)
  }
  if (is.null(requested)) {
    return(if (length(present) == 1) present else "log_ratio")
  }
  requested <- match.arg(
    requested,
    choices = c("log_ratio", "purity_corrected")
  )
  if (!requested %in% present) {
    stop(
      "`change_transform = \"",
      requested,
      "\"` is not present in the data; available: ",
      paste0("\"", present, "\"", collapse = ", "),
      ". Recompute with that axis requested.",
      call. = FALSE
    )
  }
  requested
}


#' Y-Axis Label for a Change Axis
#'
#' @param change_transform A change-transform name, or `NULL` for the historical default.
#' @return A `latex2exp` expression.
#' @noRd
.change_axis_label <- function(change_transform) {
  if (identical(change_transform, "purity_corrected")) {
    latex2exp::TeX(
      "Purity-corrected: $\\log_{2}\\left(\\frac{a_{t}}{b_{t}}\\right)$"
    )
  } else {
    latex2exp::TeX(
      "Log-Ratio: $\\log_{2}\\left(\\frac{\\fraction_{1}}{\\fraction_{2}}\\right)$"
    )
  }
}
