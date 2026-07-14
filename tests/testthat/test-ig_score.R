library(PhyloIgSeq)

# ---- compute_ig_score ----

test_that("compute_ig_score computes 'palm' as the ratio of relative abundances", {
  pos <- c(50, 30, 20, 5)
  neg <- c(5, 10, 40, 45)

  result <- compute_ig_score(method = "palm", pos = pos, neg = neg)

  pos_abund <- pos / sum(pos)
  neg_abund <- neg / sum(neg)
  expect_equal(result, pos_abund / neg_abund)
})

test_that("compute_ig_score computes 'kau' per its closed-form formula", {
  pos <- c(50, 30, 20, 5)
  neg <- c(5, 10, 40, 45)

  result <- compute_ig_score(method = "kau", pos = pos, neg = neg)

  pos_abund <- pos / sum(pos)
  neg_abund <- neg / sum(neg)
  expected <- -log2(pos_abund / neg_abund) / log10(pos_abund * neg_abund)
  expect_equal(result, expected)
})

test_that("compute_ig_score computes 'prob_index' from pos and pre abundances", {
  pos <- c(50, 30, 20, 5)
  pre <- c(20, 20, 20, 40)
  ig_freq <- 0.3

  result <- compute_ig_score(method = "prob_index", pos = pos, pre = pre, ig_freq = ig_freq)

  pos_abund <- pos / sum(pos)
  pre_abund <- pre / sum(pre)
  expect_equal(result, pos_abund * ig_freq / pre_abund)
})

test_that("compute_ig_score computes 'prob_ratio' from pos, neg, and ig_freq", {
  pos <- c(50, 30, 20, 5)
  neg <- c(5, 10, 40, 45)
  ig_freq <- 0.3

  result <- compute_ig_score(method = "prob_ratio", pos = pos, neg = neg, ig_freq = ig_freq)

  pos_abund <- pos / sum(pos)
  neg_abund <- neg / sum(neg)
  expected <- log2(pos_abund * ig_freq / (neg_abund * (1 - ig_freq)))
  expect_equal(result, expected)
})

test_that("compute_ig_score computes 'purity_corrected_prob_index' from all purity/fraction inputs", {
  pos <- c(50, 30, 20, 5)
  neg <- c(5, 10, 40, 45)
  pre <- c(20, 20, 20, 40)

  result <- compute_ig_score(
    method = "purity_corrected_prob_index",
    pos = pos,
    neg = neg,
    pre = pre,
    pos_purity = 0.9,
    neg_impurity = 0.1,
    pos_fraction = 0.4,
    neg_fraction = 0.6
  )

  pos_abund <- pos / sum(pos)
  neg_abund <- neg / sum(neg)
  pre_abund <- pre / sum(pre)
  expected <- (pos_abund * 0.9 * 0.4 + neg_abund * 0.1 * 0.6) / pre_abund
  expect_equal(result, expected)
})

test_that("compute_ig_score computes 'purity_corrected_prob_ratio' from all purity/fraction inputs", {
  pos <- c(50, 30, 20, 5)
  neg <- c(5, 10, 40, 45)

  result <- compute_ig_score(
    method = "purity_corrected_prob_ratio",
    pos = pos,
    neg = neg,
    pos_purity = 0.9,
    neg_impurity = 0.1,
    pos_fraction = 0.4,
    neg_fraction = 0.6
  )

  pos_abund <- pos / sum(pos)
  neg_abund <- neg / sum(neg)
  prob <- pos_abund * 0.9 * 0.4 + neg_abund * 0.1 * 0.6
  expect_equal(result, log2(prob / (1 - prob)))
})

test_that("compute_ig_score converts NaN/Inf results to NA", {
  # A zero-count taxon in neg makes neg_abund 0 for that taxon, so palm's
  # pos_abund / neg_abund is Inf; compute_ig_score() must map that to NA.
  pos <- c(10, 20)
  neg <- c(5, 0)

  result <- compute_ig_score(method = "palm", pos = pos, neg = neg)

  expect_false(is.infinite(result[2]))
  expect_true(is.na(result[2]))
})

test_that("compute_ig_score rejects an unrecognized method via match.arg", {
  expect_error(
    compute_ig_score(method = "bogus", pos = c(1, 2, 3)),
    "should be one of"
  )
})

test_that("compute_ig_score returns all-NA when a required fraction argument is missing", {
  # 'palm' needs neg; omitting it leaves neg_abund as a scalar NA, which
  # recycles to an all-NA score vector of the right length.
  result <- compute_ig_score(method = "palm", pos = c(10, 20, 30))

  expect_length(result, 3)
  expect_true(all(is.na(result)))
})

# ---- to_wider_ig_score ----

make_ig_coating_agglom <- function() {
  data.frame(
    sample_id = rep(paste0("sample_", 1:3), each = 2),
    taxon_id = rep(c("taxon_1", "taxon_2"), times = 3),
    slide_z = c(0.1, 2.2, 0.6, 1.3, -0.1, 1.5),
    palm = c(0.7, 0.4, 0.2, 0.2, 0.7, 0.1),
    not_a_score = 1:6
  )
}

test_that("to_wider_ig_score defaults scores to columns shared with IG_SCORES", {
  result <- to_wider_ig_score(make_ig_coating_agglom())

  expect_named(result, c("slide_z", "palm"))
})

test_that("to_wider_ig_score pivots each requested score to one row per sample_id and one column per taxon_id", {
  result <- to_wider_ig_score(make_ig_coating_agglom(), scores = "slide_z")

  expect_named(result, "slide_z")
  wide <- result$slide_z
  expect_equal(wide$sample_id, c("sample_1", "sample_2", "sample_3"))
  expect_setequal(colnames(wide), c("sample_id", "taxon_1", "taxon_2"))
  expect_equal(wide$taxon_1, c(0.1, 0.6, -0.1))
  expect_equal(wide$taxon_2, c(2.2, 1.3, 1.5))
})

test_that("to_wider_ig_score drops taxa not present in at least shared_by of samples", {
  ig_coating_agglom <- rbind(
    make_ig_coating_agglom(),
    data.frame(
      sample_id = "sample_1",
      taxon_id = "taxon_rare",
      slide_z = 5,
      palm = 0.9,
      not_a_score = 7
    )
  )

  # taxon_rare has a value in only 1/3 samples; a 0.5 threshold should drop it.
  result <- to_wider_ig_score(ig_coating_agglom, scores = "slide_z", shared_by = 0.5)

  expect_false("taxon_rare" %in% colnames(result$slide_z))
})

test_that("to_wider_ig_score keeps every taxon when shared_by is NULL", {
  ig_coating_agglom <- rbind(
    make_ig_coating_agglom(),
    data.frame(
      sample_id = "sample_1",
      taxon_id = "taxon_rare",
      slide_z = 5,
      palm = 0.9,
      not_a_score = 7
    )
  )

  result <- to_wider_ig_score(ig_coating_agglom, scores = "slide_z")

  expect_true("taxon_rare" %in% colnames(result$slide_z))
})
