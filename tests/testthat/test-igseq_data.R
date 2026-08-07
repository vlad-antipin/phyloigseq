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

# ---- .resolve_ig_freq_value ----

test_that(".resolve_ig_freq_value passes through a valid frequency and converts percent", {
  expect_equal(
    .resolve_ig_freq_value(0.42, "frequency", sample_id = "s1"),
    0.42
  )
  expect_equal(
    .resolve_ig_freq_value(42, "percent", sample_id = "s1"),
    0.42
  )
  # the boundaries themselves are legal probabilities
  expect_equal(.resolve_ig_freq_value(0, "frequency", sample_id = "s1"), 0)
  expect_equal(.resolve_ig_freq_value(1, "frequency", sample_id = "s1"), 1)
})

test_that(".resolve_ig_freq_value returns NA for a missing value without warning", {
  expect_silent(result <- .resolve_ig_freq_value(NA, "frequency", sample_id = "s1"))
  expect_true(is.na(result))
  # length 0: the column doesn't exist, or the sample has no row for the fraction
  expect_silent(result <- .resolve_ig_freq_value(numeric(0), "frequency", sample_id = "s1"))
  expect_true(is.na(result))
})

test_that(".resolve_ig_freq_value rejects non-numeric and out-of-range values with a warning", {
  expect_warning(
    result <- .resolve_ig_freq_value("M", "percent", sample_id = "s1"),
    "not numeric"
  )
  expect_true(is.na(result))

  expect_warning(
    result <- .resolve_ig_freq_value(1.5, "frequency", sample_id = "s1"),
    "outside the expected \\[0, 1\\]"
  )
  expect_true(is.na(result))

  # 150% is out of range only after the conversion
  expect_warning(
    result <- .resolve_ig_freq_value(150, "percent", sample_id = "s1"),
    "outside the expected \\[0, 1\\]"
  )
  expect_true(is.na(result))
})

test_that(".resolve_ig_freq_value warns when a value resolves to more than one row", {
  expect_warning(
    result <- .resolve_ig_freq_value(c(0.2, 0.3), "frequency", sample_id = "s1"),
    "resolves to 2 values"
  )
  expect_true(is.na(result))
})

# ---- resolve_ig_freqs ----

make_ig_freq_metadata <- function(values = c(whole = 0.4, pos = 0.8, neg1 = 0.1, neg2 = 0.1)) {
  data.frame(
    sample_id = "s1",
    sorting_fraction = names(values),
    ig_pheno = unname(values),
    stringsAsFactors = FALSE
  )
}

resolve <- function(metadata, layout, positive = "pos", ...) {
  resolve_ig_freqs(
    sam_metadata_df = metadata,
    fraction_id_name = "sorting_fraction",
    ig_freq_name = "ig_pheno",
    ig_freq_units = "frequency",
    ig_freq_layout = layout,
    positive_fraction_name = positive,
    first_negative_fraction_name = "neg1",
    second_negative_fraction_name = "neg2",
    presorting_fraction_name = "whole",
    sample_id = "s1",
    ...
  )
}

test_that("resolve_ig_freqs reads one value per fraction under the long layout", {
  result <- resolve(make_ig_freq_metadata(), "long")

  expect_equal(result$presort_ig_freq, 0.4)
  expect_equal(result$pos_ig_freq, 0.8)
  expect_equal(result$neg1_ig_freq, 0.1)
  expect_equal(result$neg2_ig_freq, 0.1)
})

test_that("resolve_ig_freqs follows the positive fraction it is asked for", {
  metadata <- make_ig_freq_metadata(
    c(whole = 0.4, pos = 0.8, neg1 = 0.1, neg2 = 0.25)
  )
  expect_equal(resolve(metadata, "long", positive = "pos")$pos_ig_freq, 0.8)
  expect_equal(resolve(metadata, "long", positive = "neg2")$pos_ig_freq, 0.25)
})

test_that("resolve_ig_freqs leaves a fraction NA when it has no row", {
  metadata <- make_ig_freq_metadata(c(whole = 0.4, pos = 0.8, neg1 = 0.1))
  result <- resolve(metadata, "long")

  expect_equal(result$presort_ig_freq, 0.4)
  expect_true(is.na(result$neg2_ig_freq))
})

test_that("resolve_ig_freqs collapses to presort only under the wide layout", {
  metadata <- make_ig_freq_metadata(
    c(whole = 0.3, pos = 0.3, neg1 = 0.3, neg2 = 0.3)
  )
  result <- resolve(metadata, "wide")

  expect_equal(result$presort_ig_freq, 0.3)
  # a single overall P(Ig+) says nothing about either sorted fraction's purity
  expect_true(is.na(result$pos_ig_freq))
  expect_true(is.na(result$neg1_ig_freq))
  expect_true(is.na(result$neg2_ig_freq))
})

test_that("resolve_ig_freqs warns and gives up when a wide column varies by fraction", {
  expect_warning(
    result <- resolve(make_ig_freq_metadata(), "wide"),
    'ig_freq_layout = "long"'
  )
  expect_true(all(is.na(unlist(result))))
})

test_that("resolve_ig_freqs returns all NA for a NULL or unknown ig_freq_name", {
  metadata <- make_ig_freq_metadata()

  result <- resolve_ig_freqs(
    sam_metadata_df = metadata,
    fraction_id_name = "sorting_fraction",
    ig_freq_name = NULL,
    ig_freq_units = "frequency",
    ig_freq_layout = "long",
    positive_fraction_name = "pos"
  )
  expect_true(all(is.na(unlist(result))))

  expect_warning(
    result <- resolve_ig_freqs(
      sam_metadata_df = metadata,
      fraction_id_name = "sorting_fraction",
      ig_freq_name = "no_such_column",
      ig_freq_units = "frequency",
      ig_freq_layout = "long",
      positive_fraction_name = "pos"
    ),
    "is not a column of sample_data"
  )
  expect_true(all(is.na(unlist(result))))
})

# ---- impute_zeros(): the pre-imputation zero pattern ----

test_that("impute_zeros reports which counts were really zero, per fraction", {
  df <- data.frame(
    sample_id = rep("s1", 4),
    taxon_id = paste0("taxon_", 1:4),
    Pos = c(10, 0, 5, 0),
    Neg1 = c(8, 0, 0, 0),
    Neg2 = c(9, 4, 0, 0)
  )
  res <- impute_zeros(df, c("Pos", "Neg1", "Neg2"), "pseudo_count")

  # taxon_4 is zero everywhere and is dropped before anything else happens.
  expect_equal(nrow(res$data), 3)
  expect_equal(nrow(res$was_zero), 3)
  expect_equal(colnames(res$was_zero), c("Pos", "Neg1", "Neg2"))
  expect_equal(res$was_zero$Pos, c(FALSE, TRUE, FALSE))
  expect_equal(res$was_zero$Neg1, c(FALSE, TRUE, TRUE))
  expect_equal(res$was_zero$Neg2, c(FALSE, FALSE, TRUE))
  expect_identical(rownames(res$was_zero), rownames(res$data))
  # The imputed frame itself no longer carries the information, which is the
  # whole reason this is returned separately.
  expect_false(any(res$data[, c("Pos", "Neg1", "Neg2")] == 0))
})

test_that("impute_zeros keeps was_zero aligned when the method drops rows too", {
  df <- data.frame(
    sample_id = rep("s1", 4),
    taxon_id = paste0("taxon_", 1:4),
    Pos = c(10, 0, 5, 3),
    Neg1 = c(8, 7, 0, 4),
    Neg2 = c(9, 4, 2, 6)
  )
  res <- impute_zeros(df, c("Pos", "Neg1", "Neg2"), "no_zero")
  # "no_zero" drops every taxon with a zero, so the pattern must be all-FALSE
  # over exactly the survivors rather than over the input rows.
  expect_equal(res$data$taxon_id, c("taxon_1", "taxon_4"))
  expect_equal(nrow(res$was_zero), 2)
  expect_false(any(as.matrix(res$was_zero)))
})
