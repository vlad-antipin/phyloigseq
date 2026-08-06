library(PhyloIgSeq)

# ---- Fixture ----
#
# A small "sorted_sample_df"-shaped data frame, as produced by
# group_sorted_samples(): one row per taxon, columns for sample/taxon id and
# one column per fraction. taxon_5 has a zero in Neg1 (but not Pos/Neg2) so
# zero-handling/imputation paths get exercised.

make_sorted_sample_df <- function() {
  data.frame(
    sample_id = rep("sample_1", 5),
    taxon_id = paste0("taxon_", 1:5),
    Pos = c(10, 20, 5, 8, 12),
    Neg1 = c(8, 15, 6, 4, 0),
    Neg2 = c(9, 18, 4, 5, 3)
  )
}

df <- make_sorted_sample_df()

# ---- get_ma_coordinates ----

test_that("get_ma_coordinates computes obs coordinates and keeps raw fraction columns", {
  ma_coords <- get_ma_coordinates(
    sorted_sample_df = df,
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1"
  )
  expect_equal(nrow(ma_coords), nrow(df))
  expect_setequal(
    colnames(ma_coords),
    c(
      "taxon_id",
      "sample_id",
      "pos",
      "neg1",
      "neg2",
      "obs_abundance",
      "obs_change",
      "null_abundance",
      "null_change",
      "change_transform",
      "obs_in_cone",
      "null_in_cone"
    )
  )
  expect_equal(ma_coords$pos, df$Pos)
  expect_equal(ma_coords$neg1, df$Neg1)
  expect_equal(unique(ma_coords$change_transform), "log_ratio")
  # The "log_ratio" axis has no admissible cone, so both flags stay NA.
  expect_true(all(is.na(ma_coords$obs_in_cone)))
  # taxon_5 has Neg1 == 0 (log2(0) = -Inf, cleaned to NA) -- checked
  # separately below; compare only the well-defined rows here.
  non_zero <- df$Neg1 != 0
  expect_equal(
    ma_coords$obs_change[non_zero],
    log2(df$Pos[non_zero]) - log2(df$Neg1[non_zero])
  )
})

test_that("get_ma_coordinates' log-ratio axis carries no per-pair normalization", {
  # Regression guard. Dividing each fraction by its own total before the log
  # gives the observed pair (Pos, Neg1) the constant log2(S_Neg1/S_Pos) and the
  # null pair (Neg1, Neg2) the constant log2(S_Neg2/S_Neg1). compute_slide_z()
  # cancels an additive constant only when both pairs carry the SAME one, so
  # under the default empirical null the residual translates every Z-score in
  # the sample -- measured at up to 1.5 Z units on ps_igseq, with no change in
  # spread to make it visible. Both pairs must stay on the values as supplied.
  ma_coords <- get_ma_coordinates(
    sorted_sample_df = df,
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1",
    second_negative_fraction_name = "Neg2"
  )
  defined <- df$Pos != 0 & df$Neg1 != 0 & df$Neg2 != 0
  obs_offset <- ma_coords$obs_change[defined] -
    (log2(df$Pos[defined]) - log2(df$Neg1[defined]))
  null_offset <- ma_coords$null_change[defined] -
    (log2(df$Neg1[defined]) - log2(df$Neg2[defined]))
  expect_equal(obs_offset, rep(0, sum(defined)))
  expect_equal(null_offset, rep(0, sum(defined)))
})

test_that("get_ma_coordinates leaves null_* columns NA without a second negative fraction", {
  ma_coords <- get_ma_coordinates(
    sorted_sample_df = df,
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1"
  )
  expect_true(all(is.na(ma_coords$null_abundance)))
  expect_true(all(is.na(ma_coords$null_change)))
  expect_true(all(is.na(ma_coords$neg2)))
})

test_that("get_ma_coordinates fills in null_* coordinates with a second negative fraction", {
  ma_coords <- get_ma_coordinates(
    sorted_sample_df = df,
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1",
    second_negative_fraction_name = "Neg2"
  )
  expect_equal(ma_coords$neg2, df$Neg2)
  # taxon_5 has Neg1 == 0 (log2(0) = -Inf, cleaned to NA); compare only the
  # well-defined rows here.
  non_zero <- df$Neg1 != 0
  expect_equal(
    ma_coords$null_change[non_zero],
    log2(df$Neg1[non_zero]) - log2(df$Neg2[non_zero])
  )
})

test_that("get_ma_coordinates sets NA for a zero-count fraction (log of zero)", {
  ma_coords <- get_ma_coordinates(
    sorted_sample_df = df,
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1"
  )
  zero_row <- which(df$Neg1 == 0)
  expect_true(is.na(ma_coords$obs_abundance[zero_row]))
  expect_true(is.na(ma_coords$obs_change[zero_row]))
})

# ---- get_ma_coordinates: the purity-corrected change axis ----
#
# Non-canonical fraction names and a deliberately imperfect sort (p = 0.7,
# q = 0.15 => admissible cone pos/neg in [(1-p)/(1-q), p/q] = [0.353, 4.667]).

purity_coords <- function(data = df, p = 0.7, q = 0.15, ...) {
  get_ma_coordinates(
    sorted_sample_df = data,
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1",
    second_negative_fraction_name = "Neg2",
    change_transform = "purity_corrected",
    pos_ig_freq = p,
    neg_ig_freq = q,
    ...
  )
}

test_that("the purity-corrected axis keeps A identical and records the transform", {
  log_ratio <- get_ma_coordinates(
    sorted_sample_df = df,
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1",
    second_negative_fraction_name = "Neg2"
  )
  purity <- purity_coords()

  # A is a precision proxy (how much evidence a taxon carries) and must not
  # follow the change transform.
  expect_equal(purity$obs_abundance, log_ratio$obs_abundance)
  expect_equal(purity$null_abundance, log_ratio$null_abundance)
  expect_equal(unique(purity$change_transform), "purity_corrected")
  expect_false(isTRUE(all.equal(purity$obs_change, log_ratio$obs_change)))
})

test_that("the purity-corrected axis is 0 exactly where the compositions match", {
  # Equal relative abundances in both fractions => a_t == b_t => change 0, for
  # any valid p > q and any pool prior, because the prior is symmetric.
  equal_df <- df
  equal_df$Neg1 <- equal_df$Pos
  equal_df$Neg2 <- equal_df$Pos

  for (params in list(c(0.7, 0.15), c(0.55, 0.45), c(0.99, 0.01))) {
    coords <- purity_coords(equal_df, p = params[1], q = params[2])
    expect_equal(coords$obs_change, rep(0, nrow(equal_df)))
    expect_equal(coords$null_change, rep(0, nrow(equal_df)))
  }
})

test_that("the purity-corrected axis reduces to the log-ratio at a perfect sort", {
  perfect <- purity_coords(p = 1, q = 0, pool_prior = 0)
  log_ratio <- get_ma_coordinates(
    sorted_sample_df = df,
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1",
    second_negative_fraction_name = "Neg2"
  )
  # Equal up to an additive constant, not identical: the un-mixing needs
  # compositions, so the purity-corrected axis normalizes each fraction by its
  # own total while the log-ratio axis reads the values as supplied. At p = 1,
  # q = 0 the un-mixing is the identity on those compositions, so all that is
  # left between the two axes is log2(sum(Neg1)/sum(Pos)).
  defined <- !is.na(perfect$obs_change) & !is.na(log_ratio$obs_change)
  offset <- perfect$obs_change[defined] - log_ratio$obs_change[defined]
  expect_equal(offset, rep(log2(sum(df$Neg1) / sum(df$Pos)), sum(defined)))
})

test_that("the purity-corrected axis is strictly increasing in pos, past the cone ceiling", {
  # p/q = 4.667, so the last few of these are above the ceiling and saturate --
  # the ordering must survive there too, since that is what slide_z ranks on.
  changes <- vapply(
    c(1, 5, 20, 60, 100, 200, 400, 1000, 5000),
    function(k) {
      one <- data.frame(
        sample_id = "s",
        taxon_id = c("a", "b"),
        Pos = c(k, 1000),
        Neg1 = c(100, 1000),
        Neg2 = c(100, 1000)
      )
      purity_coords(one)$obs_change[1]
    },
    numeric(1)
  )
  expect_true(all(diff(changes) > 0))
  expect_true(all(is.finite(changes)))
})

test_that("in_cone marks exactly the taxa outside [(1-p)/(1-q), p/q]", {
  p <- 0.7
  q <- 0.15
  spread <- data.frame(
    sample_id = "s",
    taxon_id = paste0("t", 1:4),
    # Chosen so the relative-abundance ratio brackets both cone bounds.
    Pos = c(1, 250, 250, 999),
    Neg1 = c(999, 250, 250, 1),
    Neg2 = c(999, 250, 250, 1)
  )
  coords <- purity_coords(spread, p = p, q = q)
  ratio <- (spread$Pos / sum(spread$Pos)) / (spread$Neg1 / sum(spread$Neg1))
  expected_in_cone <- ratio >= (1 - p) / (1 - q) & ratio <= p / q
  expect_equal(coords$obs_in_cone, expected_in_cone)

  # cone_policy = "na" blanks exactly those, and nothing else.
  na_coords <- purity_coords(spread, p = p, q = q, cone_policy = "na")
  expect_equal(is.na(na_coords$obs_change), !expected_in_cone)
})

test_that("the purity-corrected axis is all-NA for a singular sort (p <= q)", {
  expect_warning(
    coords <- get_ma_coordinates(
      sorted_sample_df = df,
      positive_fraction_name = "Pos",
      first_negative_fraction_name = "Neg1",
      change_transform = "purity_corrected",
      pos_ig_freq = 0.2,
      neg_ig_freq = 0.2
    ),
    NA
  )
  expect_true(all(is.na(coords$obs_change)))
})

test_that("the purity-corrected axis warns when the Ig+ frequencies are missing", {
  expect_warning(
    coords <- get_ma_coordinates(
      sorted_sample_df = df,
      positive_fraction_name = "Pos",
      first_negative_fraction_name = "Neg1",
      change_transform = "purity_corrected"
    ),
    "needs both `pos_ig_freq` and `neg_ig_freq`"
  )
  expect_true(all(is.na(coords$obs_change)))
})

# ---- get_ma_plot_data ----

test_that("get_ma_plot_data reports sample_id and the pre-imputation zero-taxa count", {
  ma_plot_data <- get_ma_plot_data(
    sorted_sample_df = df,
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1"
  )
  expect_equal(ma_plot_data$sample_id, "sample_1")
  expect_equal(ma_plot_data$nb_zero_taxa, sum(df$Pos == 0 | df$Neg1 == 0))
})

test_that("get_ma_plot_data's plot_data has one comparison row per taxon per treatment (no empirical null)", {
  ma_plot_data <- get_ma_plot_data(
    sorted_sample_df = df,
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1",
    zero_treatments = c("keep_zeros", "pseudo_count")
  )
  expect_equal(nrow(ma_plot_data$plot_data), nrow(df) * 2)
  expect_setequal(
    levels(ma_plot_data$plot_data$zero_treatment),
    c("keep zeros", "pseudo count")
  )
})

test_that("get_ma_plot_data doubles rows per treatment when an empirical null is requested", {
  ma_plot_data <- get_ma_plot_data(
    sorted_sample_df = df,
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1",
    second_negative_fraction_name = "Neg2",
    zero_treatments = "keep_zeros"
  )
  expect_equal(nrow(ma_plot_data$plot_data), nrow(df) * 2)
  expect_setequal(
    unique(ma_plot_data$plot_data$comparison),
    c("Pos vs Neg1", "Neg1 vs Neg2")
  )
})

test_that("get_ma_plot_data unions imputed_taxa across zero_treatments regardless of order", {
  # "no_zero" resets its own imputed-taxa tracking to NULL (those taxa are
  # dropped entirely under that method); placing it last used to make it
  # silently clobber the imputed_taxa found by the earlier "pseudo_count"
  # treatment. taxon_5 (zero in Neg1) should still show up as imputed.
  ma_plot_data <- get_ma_plot_data(
    sorted_sample_df = df,
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1",
    zero_treatments = c("pseudo_count", "no_zero")
  )
  expect_true("taxon_5" %in% ma_plot_data$imputed_taxa)
})

test_that("get_ma_plot_data flags `imputed` per treatment, not as the union", {
  # The union in `imputed_taxa` is deliberately order-independent, but it must
  # not be what a treatment's own rows are labelled with: "no_zero" drops
  # taxon_5 outright and imputes nothing, so none of its rows are imputed.
  ma_plot_data <- get_ma_plot_data(
    sorted_sample_df = df,
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1",
    zero_treatments = c("pseudo_count", "no_zero")
  )
  plot_data <- ma_plot_data$plot_data
  expect_true(all(
    plot_data$imputed[plot_data$zero_treatment == "pseudo count"] ==
      (plot_data$taxon_id[plot_data$zero_treatment == "pseudo count"] ==
        "taxon_5")
  ))
  expect_false(any(plot_data$imputed[plot_data$zero_treatment == "no zero"]))
})

# ---- plot_ma ----

ma_plot_data <- get_ma_plot_data(
  sorted_sample_df = df,
  positive_fraction_name = "Pos",
  first_negative_fraction_name = "Neg1",
  second_negative_fraction_name = "Neg2",
  zero_treatments = c("keep_zeros", "pseudo_count")
)

test_that("plot_ma returns a ggplot for both supported types", {
  expect_s3_class(plot_ma(ma_plot_data, type = "facet"), "ggplot")
  expect_s3_class(plot_ma(ma_plot_data, type = "superposed"), "ggplot")
})

test_that("plot_ma errors on an unrecognized type instead of silently dropping geoms", {
  expect_error(plot_ma(ma_plot_data, type = "not_a_real_type"))
})

# Every x/y actually rendered, pooled across layers and panels.
drawn_points <- function(plt) {
  # A layer with no rows is built without x/y columns at all, so skip it.
  layers <- lapply(ggplot2::ggplot_build(plt)$data, function(d) {
    if (nrow(d) == 0) NULL else d[, c("x", "y")]
  })
  drawn <- do.call(rbind, layers)
  drawn[!is.na(drawn$x) & !is.na(drawn$y), , drop = FALSE]
}

test_that("plot_ma jitters imputed taxa into the off-axis band", {
  ma_plot_data <- get_ma_plot_data(
    sorted_sample_df = df,
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1",
    zero_treatments = "pseudo_count"
  )
  plot_data <- ma_plot_data$plot_data
  expect_true(any(plot_data$imputed)) # fixture must exercise the path
  drawn <- drawn_points(plot_ma(ma_plot_data))
  # The band sits below every real A, so an imputed taxon is never drawn among
  # the taxa that kept their own log-abundance.
  expect_lt(min(drawn$x), min(plot_data$A[!plot_data$imputed], na.rm = TRUE))
})

test_that("plot_ma reads `imputed` per treatment, so a treatment that imputes nothing has an empty band", {
  # The `imputed_taxa` union names taxon_5 because "pseudo_count" imputed it;
  # under "no_zero" that taxon is dropped outright and nothing is imputed, so
  # nothing may be banished off-axis in that facet. Marking it imputed
  # everywhere is what made every facet render the same points the same way.
  ma_plot_data <- get_ma_plot_data(
    sorted_sample_df = df,
    positive_fraction_name = "Pos",
    first_negative_fraction_name = "Neg1",
    zero_treatments = c("pseudo_count", "no_zero")
  )
  expect_true("taxon_5" %in% ma_plot_data$imputed_taxa)
  no_zero_rows <- ma_plot_data$plot_data[
    ma_plot_data$plot_data$zero_treatment == "no zero",
  ]
  expect_gt(nrow(no_zero_rows), 0)
  expect_false(any(no_zero_rows$imputed))
})
