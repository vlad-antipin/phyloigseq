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
