library(PhyloIgSeq)

# ---- compute_slide_z ----
#
# ma_coords-shaped fixtures, built directly (not via get_ma_coordinates) so
# obs_change/null_change can be pinned to values whose expected Z-scores are
# easy to hand-verify.

make_ma_coords <- function(obs_abundance, obs_change, null_change = NA_real_) {
  n <- length(obs_abundance)
  data.frame(
    taxon_id = paste0("taxon_", seq_len(n)),
    sample_id = "sample_1",
    obs_abundance = obs_abundance,
    obs_change = obs_change,
    null_abundance = rep(NA_real_, n),
    null_change = rep(null_change, length.out = n)
  )
}

test_that("compute_slide_z matches a hand-computed Z-score for a single window", {
  # window_size >= nrow(ma_coords), so this is the "only one window" branch:
  # every taxon is centered/scaled against the same overall mean/sd.
  obs_change <- c(1, 2, 3, 4, 5, 6, 7)
  ma_coords <- make_ma_coords(obs_abundance = 7:1, obs_change = obs_change)

  result <- compute_slide_z(
    ma_coords,
    window_size = 50,
    empirical_null_distribution = FALSE
  )

  expected <- (obs_change - mean(obs_change)) / sd(obs_change)
  expect_equal(result, expected)
})

test_that("compute_slide_z restores the caller's original row order after internally sorting by abundance", {
  obs_abundance <- 10:1
  obs_change <- 1:10
  ma_coords_sorted <- make_ma_coords(obs_abundance, obs_change)

  set.seed(42)
  shuffle <- sample(nrow(ma_coords_sorted))
  ma_coords_shuffled <- ma_coords_sorted[shuffle, ]

  result_sorted <- compute_slide_z(
    ma_coords_sorted,
    window_size = 5,
    empirical_null_distribution = FALSE
  )
  result_shuffled <- compute_slide_z(
    ma_coords_shuffled,
    window_size = 5,
    empirical_null_distribution = FALSE
  )

  # Match by taxon_id rather than assuming any particular row order.
  names(result_sorted) <- ma_coords_sorted$taxon_id
  names(result_shuffled) <- ma_coords_shuffled$taxon_id
  expect_equal(
    result_shuffled[ma_coords_sorted$taxon_id],
    result_sorted[ma_coords_sorted$taxon_id],
    ignore_attr = TRUE
  )
})

test_that("compute_slide_z uses the empirical null instead of obs_change when requested", {
  obs_change <- c(1, 2, 3, 4, 5)
  null_change <- c(100, 102, 104, 106, 108)
  ma_coords <- make_ma_coords(
    obs_abundance = 5:1,
    obs_change = obs_change,
    null_change = null_change
  )

  result <- compute_slide_z(
    ma_coords,
    window_size = 50,
    empirical_null_distribution = TRUE
  )

  expected <- (obs_change - mean(null_change)) / sd(null_change)
  expect_equal(result, expected)
})

test_that("compute_slide_z scores imputed taxa from their own change distribution, not the sliding window", {
  obs_change <- c(1, 2, 100, 4, 5, 6, 7, 200, 9, 10)
  ma_coords <- make_ma_coords(obs_abundance = 10:1, obs_change = obs_change)
  was_imputed <- c(
    FALSE, FALSE, TRUE, FALSE, FALSE,
    FALSE, FALSE, TRUE, FALSE, FALSE
  )

  result <- compute_slide_z(
    ma_coords,
    was_imputed = was_imputed,
    window_size = 5,
    empirical_null_distribution = FALSE
  )

  imputed_vals <- obs_change[was_imputed]
  expected_imputed <- (imputed_vals - mean(imputed_vals)) / sd(imputed_vals)

  expect_length(result, nrow(ma_coords))
  expect_equal(result[was_imputed], expected_imputed)
  expect_false(anyNA(result[!was_imputed]))
})

test_that("compute_slide_z warns and returns NULL when every taxon is imputed", {
  ma_coords <- make_ma_coords(obs_abundance = 3:1, obs_change = c(1, 2, 3))

  expect_warning(
    result <- compute_slide_z(
      ma_coords,
      was_imputed = c(TRUE, TRUE, TRUE),
      window_size = 5
    ),
    "No taxa in MA coordinates"
  )
  expect_null(result)
})

# A denser fixture (window_size << nrow) to exercise the multi-window path
# (first/middle/last branches of the sliding window) end to end.
test_that("compute_slide_z tiles multiple windows without dropping or duplicating any taxon", {
  set.seed(7)
  n <- 137
  ma_coords <- make_ma_coords(
    obs_abundance = rnorm(n),
    obs_change = rnorm(n)
  )
  result <- compute_slide_z(
    ma_coords,
    window_size = 20,
    empirical_null_distribution = FALSE
  )
  expect_length(result, n)
  expect_false(anyNA(result))
})

# ---- compute_slide_z: the estimation mask and the two populations ----

test_that("compute_slide_z is invariant to an affine change of the Y axis", {
  # The property the pluggable change axis rests on: because the score centers
  # and scales, rescaling and shifting obs_change (and the null it is measured
  # against) cannot move any Z-score. Only a NON-affine axis, i.e.
  # change_transform = "purity_corrected", can.
  set.seed(11)
  n <- 80
  ma_coords <- make_ma_coords(
    obs_abundance = rnorm(n),
    obs_change = rnorm(n),
    null_change = rnorm(n)
  )
  shifted <- ma_coords
  shifted$obs_change <- 3.7 * ma_coords$obs_change - 12.5
  shifted$null_change <- 3.7 * ma_coords$null_change - 12.5

  for (empirical in c(TRUE, FALSE)) {
    expect_equal(
      compute_slide_z(
        shifted,
        window_size = 20,
        empirical_null_distribution = empirical
      ),
      compute_slide_z(
        ma_coords,
        window_size = 20,
        empirical_null_distribution = empirical
      )
    )
  }
})

test_that("compute_slide_z ignores an all-in-cone mask (regression guard)", {
  # Adding the in_cone columns must not perturb anything when nothing is out of
  # cone -- which is always the case on the "log_ratio" axis, where they are NA.
  set.seed(13)
  n <- 60
  bare <- make_ma_coords(
    obs_abundance = rnorm(n),
    obs_change = rnorm(n),
    null_change = rnorm(n)
  )
  flagged <- bare
  flagged$obs_in_cone <- NA
  flagged$null_in_cone <- NA
  all_true <- bare
  all_true$obs_in_cone <- TRUE
  all_true$null_in_cone <- TRUE

  reference <- compute_slide_z(bare, window_size = 15)
  expect_equal(compute_slide_z(flagged, window_size = 15), reference)
  expect_equal(compute_slide_z(all_true, window_size = 15), reference)
})

test_that("compute_slide_z keeps out-of-cone taxa out of the null but still scores them", {
  obs_change <- c(1, 2, 3, 4, 5, 6, 7)
  ma_coords <- make_ma_coords(
    obs_abundance = 7:1,
    obs_change = obs_change,
    null_change = c(0.1, 0.2, 999, 0.3, 0.1, 0.2, 0.3)
  )
  # Taxon 3's null pair left the cone, so its (saturated) null value must not
  # inform the mean/sd -- but the taxon still gets a score.
  ma_coords$obs_in_cone <- TRUE
  ma_coords$null_in_cone <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)

  result <- compute_slide_z(
    ma_coords,
    window_size = 50,
    empirical_null_distribution = TRUE
  )

  clean_null <- ma_coords$null_change[ma_coords$null_in_cone]
  expected <- (obs_change - mean(clean_null)) / sd(clean_null)
  expect_equal(result, expected)
  expect_false(anyNA(result))
})

test_that("compute_slide_z does not apply the cone mask inside the imputed population", {
  # Regression guard. Inside the imputed population the cone flag reports the
  # imputation pattern, not measurement quality: a taxon zeroed in BOTH negative
  # fractions gets the same substituted value twice, so its ratio is exactly 1,
  # it lands in the cone, and its null_change is a manufactured 0. A taxon zeroed
  # in only one gets an extreme ratio and falls outside. Masking on the cone here
  # therefore keeps the fabricated zeros and drops the real spread -- on ps_igseq
  # that shrank the reference sd 3x and inflated every imputed |Z| ~2.5x.
  n_clean <- 6
  imputed_null <- c(rep(0, 6), c(-3.1, 2.8, -2.4, 3.3))
  in_cone <- c(rep(TRUE, 6), rep(FALSE, 4)) # exactly the fabricated-zero split
  ma_coords <- make_ma_coords(
    obs_abundance = c(20:15, 10:1),
    obs_change = c(seq_len(n_clean), seq(2, 20, length.out = 10))
  )
  ma_coords$null_change <- c(c(0.4, -0.3, 0.5, -0.2, 0.1, -0.5), imputed_null)
  ma_coords$obs_in_cone <- TRUE
  ma_coords$null_in_cone <- c(rep(TRUE, n_clean), in_cone)
  was_imputed <- c(rep(FALSE, n_clean), rep(TRUE, 10))

  result <- compute_slide_z(
    ma_coords,
    was_imputed = was_imputed,
    window_size = 50,
    empirical_null_distribution = TRUE
  )

  # Every imputed taxon informs the estimate, in-cone or not.
  expected <- (ma_coords$obs_change[was_imputed] - mean(imputed_null)) /
    sd(imputed_null)
  expect_equal(result[was_imputed], expected)

  # ...and masking on the cone would have used only the fabricated zeros, whose
  # sd is far too small -- which is what this guards against.
  fabricated <- imputed_null[in_cone]
  expect_lt(sd(fabricated), sd(imputed_null) / 2)

  # The cone mask still applies to the non-imputed population.
  expect_false(anyNA(result))
})

test_that("compute_slide_z's mask follows whichever pair supplies the reference", {
  ma_coords <- make_ma_coords(
    obs_abundance = 6:1,
    obs_change = c(1, 2, 3, 4, 5, 6),
    null_change = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
  )
  # Flags disagree: taxon 1 is out of cone on the observed pair only, taxon 6 on
  # the null pair only.
  ma_coords$obs_in_cone <- c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)
  ma_coords$null_in_cone <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)

  empirical <- compute_slide_z(
    ma_coords,
    window_size = 50,
    empirical_null_distribution = TRUE
  )
  # Empirical mode reads null_change, so it must drop taxon 6, not taxon 1.
  ref_null <- ma_coords$null_change[c(1:5)]
  expect_equal(
    empirical,
    (ma_coords$obs_change - mean(ref_null)) / sd(ref_null)
  )

  theoretical <- compute_slide_z(
    ma_coords,
    window_size = 50,
    empirical_null_distribution = FALSE
  )
  # Theoretical mode reads obs_change, so it must drop taxon 1 instead.
  ref_obs <- ma_coords$obs_change[2:6]
  expect_equal(
    theoretical,
    (ma_coords$obs_change - mean(ref_obs)) / sd(ref_obs)
  )
})

test_that("compute_slide_z pools the imputed population instead of windowing it", {
  # Regression guard. Windowing this population is tempting -- it spans an
  # abundance range like any other -- but that range is largely an artefact of
  # the zero treatment, which substitutes one value for every zero: on ps_igseq
  # 75-83% of imputed taxa then share an obs_abundance with another. Windows cut
  # from that tied ordering are arbitrary slices whose sd underestimates the
  # population spread by up to 5x, inflating |Z| rather than localizing it (it
  # turned 320 significant imputed taxa into 1527). Enough imputed taxa here to
  # fill several windows, so the two schemes are distinguishable.
  n_per <- 12
  obs_abundance <- c(seq(100, 89), seq(10, -1))
  obs_change <- c(seq(1, 12), seq(100, 1200, length.out = n_per))
  ma_coords <- make_ma_coords(obs_abundance, obs_change)
  was_imputed <- rep(TRUE, 2 * n_per)
  # Keep a handful of non-imputed taxa so the function doesn't bail out.
  ma_coords <- rbind(
    make_ma_coords(c(50, 49, 48, 47), c(1, 2, 3, 4)),
    ma_coords
  )
  ma_coords$taxon_id <- paste0("taxon_", seq_len(nrow(ma_coords)))
  was_imputed <- c(rep(FALSE, 4), was_imputed)

  slide_z <- compute_slide_z(
    ma_coords,
    was_imputed = was_imputed,
    window_size = n_per,
    empirical_null_distribution = FALSE
  )
  # One mean/sd over every imputed taxon, whatever its abundance.
  imputed_vals <- ma_coords$obs_change[was_imputed]
  pooled <- (imputed_vals - mean(imputed_vals)) / sd(imputed_vals)

  expect_length(slide_z, nrow(ma_coords))
  expect_false(anyNA(slide_z[was_imputed]))
  expect_equal(slide_z[was_imputed], pooled)
  # The non-imputed taxa are unaffected: they are still windowed among themselves.
  expect_false(anyNA(slide_z[!was_imputed]))
})

test_that("compute_slide_z pools a two-taxon imputed population", {
  # The pooled estimate does not depend on how many imputed taxa there are;
  # two of them get the same treatment as the several-windows-worth above.
  obs_change <- c(1, 2, 100, 4, 5, 6, 7, 200, 9, 10)
  ma_coords <- make_ma_coords(obs_abundance = 10:1, obs_change = obs_change)
  was_imputed <- c(
    FALSE, FALSE, TRUE, FALSE, FALSE,
    FALSE, FALSE, TRUE, FALSE, FALSE
  )

  result <- compute_slide_z(
    ma_coords,
    was_imputed = was_imputed,
    window_size = 50,
    empirical_null_distribution = FALSE
  )
  imputed_vals <- obs_change[was_imputed]
  expect_equal(
    result[was_imputed],
    (imputed_vals - mean(imputed_vals)) / sd(imputed_vals)
  )
})

# ---- get_slide_z ----

make_sorted_sample_df <- function() {
  data.frame(
    sample_id = rep("sample_1", 5),
    taxon_id = paste0("taxon_", 1:5),
    Pos = c(10, 20, 5, 8, 12),
    Neg1 = c(8, 15, 6, 4, 9),
    Neg2 = c(9, 18, 4, 5, 10)
  )
}

test_that("get_slide_z warns and returns an NA placeholder list when a fraction column is missing", {
  df <- make_sorted_sample_df()

  expect_warning(
    result <- get_slide_z(
      sorted_sample_df = df,
      positive_fraction_name = "Pos",
      first_negative_fraction_name = "Neg1",
      second_negative_fraction_name = "DoesNotExist"
    ),
    "lacks fraction"
  )
  expect_true(is.na(result$slide_z))
  expect_equal(result$ma_coords, data.frame())
  expect_true(is.na(result$ellipse_level))
  expect_equal(result$ellipse_coords, data.frame())
})

test_that("get_slide_z returns per-taxon Z-scores and MA coordinates for a normal call", {
  df <- make_sorted_sample_df()

  result <- get_slide_z(
    sorted_sample_df = df,
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1",
    second_negative_fraction_name = "Neg2",
    window_size = 5
  )

  expect_length(result$slide_z, nrow(df))
  expect_equal(nrow(result$ma_coords), nrow(df))
  expect_null(result$ellipse_level)
  expect_equal(result$ellipse_coords, data.frame())
})

test_that("get_slide_z falls back to non-empirical scoring with a warning when there is no second negative fraction", {
  df <- make_sorted_sample_df()

  expect_warning(
    result <- get_slide_z(
      sorted_sample_df = df,
      positive_fraction_name = "Pos",
      first_negative_fraction_name = "Neg1",
      second_negative_fraction_name = NULL,
      empirical_null_distribution = TRUE,
      window_size = 5
    ),
    "No second negative fraction"
  )
  expect_length(result$slide_z, nrow(df))
  expect_false(anyNA(result$slide_z))
})

test_that("get_slide_z includes confidence ellipse data when confidence_levels is supplied", {
  df <- make_sorted_sample_df()

  result <- get_slide_z(
    sorted_sample_df = df,
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1",
    second_negative_fraction_name = "Neg2",
    window_size = 5,
    confidence_levels = 0.95
  )

  expect_length(result$ellipse_level, nrow(df))
  expect_gt(nrow(result$ellipse_coords), 0)
})

# ---- get_ellipse_data ----

n_inliers <- 30

make_ellipse_fixture <- function() {
  set.seed(123)
  data.frame(
    sample_id = "sample_1",
    taxon_id = paste0("taxon_", seq_len(n_inliers + 1)),
    obs_abundance = c(rnorm(n_inliers, mean = 5, sd = 1), 5),
    # One clear outlier: typical abundance, wildly extreme change.
    obs_change = c(rnorm(n_inliers, mean = 0, sd = 1), 50),
    null_abundance = NA_real_,
    null_change = NA_real_
  )
}

test_that("get_ellipse_data flags a clear outlier as outside every confidence ellipse", {
  df <- make_ellipse_fixture()

  result <- get_ellipse_data(
    sorted_sample_df = df,
    empirical_null_distribution = FALSE,
    confidence_levels = c(0.95, 0.99)
  )

  expect_equal(as.character(result$levels[n_inliers + 1]), "0.99")
  expect_true(any(result$levels[seq_len(n_inliers)] == "ns"))
})

test_that("get_ellipse_data warns and returns NULL levels with too few points", {
  df <- make_ellipse_fixture()[1:2, ]

  expect_warning(
    result <- get_ellipse_data(
      sorted_sample_df = df,
      empirical_null_distribution = FALSE,
      confidence_levels = 0.95
    ),
    "Cannot build ellipse"
  )
  expect_null(result$levels)
  expect_equal(result$coords, data.frame())
})

test_that("get_ellipse_data falls back to observed coordinates with a warning when no empirical null is available", {
  df <- make_ellipse_fixture()

  expect_warning(
    result_requested_empirical <- get_ellipse_data(
      sorted_sample_df = df,
      empirical_null_distribution = TRUE,
      confidence_levels = c(0.95, 0.99)
    ),
    "No MA coordinates to model empirical null"
  )
  result_observed <- get_ellipse_data(
    sorted_sample_df = df,
    empirical_null_distribution = FALSE,
    confidence_levels = c(0.95, 0.99)
  )

  expect_equal(result_requested_empirical$levels, result_observed$levels)
  expect_equal(result_requested_empirical$coords, result_observed$coords)
})

test_that("get_ellipse_data assigns NA confidence level to imputed taxa", {
  df <- make_ellipse_fixture()

  result <- get_ellipse_data(
    sorted_sample_df = df,
    imputed_taxa = df$taxon_id[1],
    empirical_null_distribution = FALSE,
    confidence_levels = c(0.95, 0.99)
  )

  expect_true(is.na(result$levels[1]))
})
