library(PhyloIgSeq)

# ---- Fixture ----
#
# A small phyloseq object mimicking IgSeq sort-fraction data: several
# "biological" samples (bio_sample), each split across sort fractions
# (fraction), with one instance of each of group_sorted_samples()'s
# exclusion cases baked in (duplicated fraction, single fraction, an
# all-zero fraction).

make_igseq_ps <- function() {
  taxa <- paste0("ASV", 1:5)

  counts <- cbind(
    # sample_1: clean, 3 fractions, has some per-taxon zero cells
    s1_pos = c(10, 0, 5, 20, 3),
    s1_neg1 = c(8, 2, 6, 15, 4),
    s1_neg2 = c(9, 1, 4, 18, 2),
    # sample_2: clean, 3 fractions
    s2_pos = c(12, 3, 0, 22, 5),
    s2_neg1 = c(11, 2, 1, 19, 6),
    s2_neg2 = c(10, 4, 2, 20, 7),
    # sample_dup: "Pos" fraction duplicated -> excluded
    dup_a = c(5, 5, 5, 5, 5),
    dup_b = c(6, 6, 6, 6, 6),
    dup_c = c(7, 7, 7, 7, 7),
    # sample_single: only one fraction present -> excluded
    single_a = c(4, 4, 4, 4, 4),
    # sample_zero: "Neg2" fraction has zero reads for every taxon -> excluded
    zero_pos = c(9, 1, 2, 3, 4),
    zero_neg1 = c(8, 1, 2, 3, 4),
    zero_neg2 = c(0, 0, 0, 0, 0)
  )
  rownames(counts) <- taxa

  sdata <- data.frame(
    bio_sample = c(
      "sample_1",
      "sample_1",
      "sample_1",
      "sample_2",
      "sample_2",
      "sample_2",
      "sample_dup",
      "sample_dup",
      "sample_dup",
      "sample_single",
      "sample_zero",
      "sample_zero",
      "sample_zero"
    ),
    fraction = c(
      "Pos",
      "Neg1",
      "Neg2",
      "Pos",
      "Neg1",
      "Neg2",
      "Pos",
      "Pos",
      "Neg1",
      "Pos",
      "Pos",
      "Neg1",
      "Neg2"
    ),
    row.names = colnames(counts)
  )

  phyloseq(
    otu_table(counts, taxa_are_rows = TRUE),
    sample_data(sdata)
  )
}

ps <- make_igseq_ps()

# ---- group_sorted_samples ----

test_that("group_sorted_samples keeps only well-formed samples", {
  grouped <- suppressWarnings(group_sorted_samples(
    physeq = ps,
    sample_id_name = "bio_sample",
    fraction_id_name = "fraction",
    fraction_ids = c("Pos", "Neg1", "Neg2"),
    rarefy_by_sample = FALSE
  ))
  expect_named(grouped, c("sample_1", "sample_2"))
})

test_that("group_sorted_samples warns and excludes duplicated-fraction samples", {
  expect_warning(
    group_sorted_samples(
      physeq = ps,
      sample_id_name = "bio_sample",
      sample_ids = "sample_dup",
      fraction_id_name = "fraction",
      fraction_ids = c("Pos", "Neg1", "Neg2"),
      rarefy_by_sample = FALSE
    ),
    "duplicated fraction"
  )
})

test_that("group_sorted_samples warns and excludes single-fraction samples", {
  expect_warning(
    group_sorted_samples(
      physeq = ps,
      sample_id_name = "bio_sample",
      sample_ids = "sample_single",
      fraction_id_name = "fraction",
      fraction_ids = c("Pos", "Neg1", "Neg2"),
      rarefy_by_sample = FALSE
    ),
    "only one or no fraction"
  )
})

test_that("group_sorted_samples warns and excludes samples with a zero-read fraction", {
  expect_warning(
    group_sorted_samples(
      physeq = ps,
      sample_id_name = "bio_sample",
      sample_ids = "sample_zero",
      fraction_id_name = "fraction",
      fraction_ids = c("Pos", "Neg1", "Neg2"),
      rarefy_by_sample = FALSE
    ),
    "no reads for at least one fraction"
  )
})

test_that("group_sorted_samples output has one row per taxon and a column per fraction", {
  grouped <- suppressWarnings(group_sorted_samples(
    physeq = ps,
    sample_id_name = "bio_sample",
    sample_ids = "sample_1",
    fraction_id_name = "fraction",
    fraction_ids = c("Pos", "Neg1", "Neg2"),
    rarefy_by_sample = FALSE
  ))
  df <- grouped[["sample_1"]]
  expect_setequal(
    colnames(df),
    c("sample_id", "taxon_id", "Pos", "Neg1", "Neg2")
  )
  expect_equal(nrow(df), 5)
  expect_true(all(df$sample_id == "sample_1"))
})

test_that("group_sorted_samples rarefaction equalizes per-fraction totals", {
  grouped <- suppressWarnings(group_sorted_samples(
    physeq = ps,
    sample_id_name = "bio_sample",
    sample_ids = "sample_1",
    fraction_id_name = "fraction",
    fraction_ids = c("Pos", "Neg1", "Neg2"),
    rarefy_by_sample = TRUE
  ))
  df <- grouped[["sample_1"]]
  totals <- colSums(df[, c("Pos", "Neg1", "Neg2")])
  expect_equal(unname(totals["Pos"]), unname(totals["Neg1"]))
  expect_equal(unname(totals["Neg1"]), unname(totals["Neg2"]))
})

test_that("group_sorted_samples compositional transform makes fractions sum to 1", {
  grouped <- suppressWarnings(group_sorted_samples(
    physeq = ps,
    sample_id_name = "bio_sample",
    sample_ids = "sample_1",
    fraction_id_name = "fraction",
    fraction_ids = c("Pos", "Neg1", "Neg2"),
    rarefy_by_sample = FALSE,
    transform_by_sample = "compositional"
  ))
  totals <- colSums(grouped[["sample_1"]][, c("Pos", "Neg1", "Neg2")])
  expect_equal(unname(totals), c(1, 1, 1), tolerance = 1e-10)
})

test_that("group_sorted_samples errors on a non-phyloseq input", {
  expect_error(
    group_sorted_samples(
      physeq = data.frame(),
      sample_id_name = "bio_sample",
      fraction_id_name = "fraction"
    ),
    "phyloseq"
  )
})

# ---- impute_zeros ----

grouped_s1 <- suppressWarnings(group_sorted_samples(
  physeq = ps,
  sample_id_name = "bio_sample",
  sample_ids = "sample_1",
  fraction_id_name = "fraction",
  fraction_ids = c("Pos", "Neg1", "Neg2"),
  rarefy_by_sample = FALSE
))[["sample_1"]]

test_that("impute_zeros 'no_zero' drops every taxon with a zero anywhere", {
  result <- impute_zeros(
    data = grouped_s1,
    fraction_names = c("Pos", "Neg1", "Neg2"),
    method = "no_zero"
  )
  expect_true(all(result$data[, c("Pos", "Neg1", "Neg2")] != 0))
  expect_null(result$imputed_taxa)
})

test_that("impute_zeros 'pseudo_count' removes remaining zeros without dropping taxa", {
  result <- impute_zeros(
    data = grouped_s1,
    fraction_names = c("Pos", "Neg1", "Neg2"),
    method = "pseudo_count"
  )
  expect_equal(nrow(result$data), nrow(grouped_s1))
  expect_true(all(result$data[, c("Pos", "Neg1", "Neg2")] > 0))
  expect_true(length(result$imputed_taxa) >= 1)
})

test_that("impute_zeros 'random_pseudo_count' removes remaining zeros without dropping taxa", {
  set.seed(1)
  result <- impute_zeros(
    data = grouped_s1,
    fraction_names = c("Pos", "Neg1", "Neg2"),
    method = "random_pseudo_count"
  )
  expect_equal(nrow(result$data), nrow(grouped_s1))
  expect_true(all(result$data[, c("Pos", "Neg1", "Neg2")] > 0))
})

test_that("impute_zeros 'keep_zeros' leaves zero counts untouched", {
  result <- impute_zeros(
    data = grouped_s1,
    fraction_names = c("Pos", "Neg1", "Neg2"),
    method = "keep_zeros"
  )
  expect_equal(nrow(result$data), nrow(grouped_s1))
  expect_true(any(result$data[, c("Pos", "Neg1", "Neg2")] == 0))
})

test_that("impute_zeros always drops taxa that are zero in every fraction", {
  data_all_zero <- grouped_s1
  data_all_zero[1, c("Pos", "Neg1", "Neg2")] <- 0
  result <- impute_zeros(
    data = data_all_zero,
    fraction_names = c("Pos", "Neg1", "Neg2"),
    method = "keep_zeros"
  )
  expect_equal(nrow(result$data), nrow(data_all_zero) - 1)
})

test_that("impute_zeros 'bayesian_inference' removes remaining zeros without dropping taxa", {
  result <- impute_zeros(
    data = grouped_s1,
    fraction_names = c("Pos", "Neg1", "Neg2"),
    method = "bayesian_inference"
  )
  expect_equal(nrow(result$data), nrow(grouped_s1))
  expect_true(all(result$data[, c("Pos", "Neg1", "Neg2")] > 0))
})

test_that("impute_zeros errors on an unrecognized method", {
  expect_error(
    impute_zeros(
      data = grouped_s1,
      fraction_names = c("Pos", "Neg1", "Neg2"),
      method = "not_a_real_method"
    ),
    "Wrong 'method' argument"
  )
})

test_that("impute_zeros returns NULL imputed_taxa when data has no taxon_id column", {
  no_taxon_id <- grouped_s1[, c("sample_id", "Pos", "Neg1", "Neg2")]
  result <- impute_zeros(
    data = no_taxon_id,
    fraction_names = c("Pos", "Neg1", "Neg2"),
    method = "pseudo_count"
  )
  expect_null(result$imputed_taxa)
})
