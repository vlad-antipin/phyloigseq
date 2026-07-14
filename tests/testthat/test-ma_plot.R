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
      "null_change"
    )
  )
  expect_equal(ma_coords$pos, df$Pos)
  expect_equal(ma_coords$neg1, df$Neg1)
  # taxon_5 has Neg1 == 0 (log2(0) = -Inf, cleaned to NA) -- checked
  # separately below; compare only the well-defined rows here.
  non_zero <- df$Neg1 != 0
  expect_equal(
    ma_coords$obs_change[non_zero],
    log2(df$Pos[non_zero]) - log2(df$Neg1[non_zero])
  )
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
