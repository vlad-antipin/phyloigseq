library(phyloseq)
library(PhyloIgSeq)

# ---- Synthetic test data helpers ----

# n_taxa × n_samples phyloseq with a Genus-level tax_table and sample_data.
# taxa_are_rows controls OTU table orientation.
make_ps <- function(
  n_taxa = 50,
  n_samples = 20,
  sparsity = 0.95,
  taxa_are_rows = TRUE,
  seed = 42
) {
  set.seed(seed)
  mat <- matrix(0L, n_taxa, n_samples)
  n_nz <- max(2L, round((1 - sparsity) * n_taxa * n_samples))
  mat[sample(n_taxa * n_samples, n_nz)] <- sample(
    1L:1000L,
    n_nz,
    replace = TRUE
  )
  rownames(mat) <- paste0("ASV", seq_len(n_taxa))
  colnames(mat) <- paste0("S", seq_len(n_samples))

  n_genus <- max(5L, n_taxa %/% 5L)
  genus_ids <- sample(seq_len(n_genus), n_taxa, replace = TRUE)
  # speedyseq::tax_glom groups by the FULL taxonomic path up to taxrank.
  # Derive all higher ranks from genus_id so every ASV sharing a genus also
  # shares the same Kingdom/Phylum/.../Family — otherwise ASVs with the same
  # Genus but different higher ranks end up in separate groups and glom does
  # not reduce the taxon count as expected.
  tax <- cbind(
    Kingdom = rep("Bacteria", n_taxa),
    Phylum = paste0(
      "Phy",
      (genus_ids - 1L) %/% max(1L, ceiling(n_genus / 3L)) + 1L
    ),
    Class = paste0(
      "Class",
      (genus_ids - 1L) %/% max(1L, ceiling(n_genus / 5L)) + 1L
    ),
    Order = paste0(
      "Order",
      (genus_ids - 1L) %/% max(1L, ceiling(n_genus / 8L)) + 1L
    ),
    Family = paste0("Family", (genus_ids - 1L) %/% 2L + 1L),
    Genus = paste0("Genus", genus_ids)
  )
  rownames(tax) <- rownames(mat)

  sdata <- data.frame(
    Group = sample(c("A", "B"), n_samples, replace = TRUE),
    row.names = colnames(mat)
  )

  if (!taxa_are_rows) {
    mat <- t(mat)
  }

  phyloseq(
    otu_table(mat, taxa_are_rows = taxa_are_rows),
    tax_table(tax),
    sample_data(sdata)
  )
}

sparsify <- function(ps) {
  otu_table(ps) <- sparse_otu_table(otu_table(ps))
  ps
}

# Canonical test objects (taxa_are_rows = TRUE and FALSE)
ps_d <- make_ps() # dense,  tar = TRUE
ps_s <- sparsify(ps_d) # sparse, tar = TRUE
ps_dt <- make_ps(taxa_are_rows = FALSE) # dense,  tar = FALSE
ps_st <- sparsify(ps_dt) # sparse, tar = FALSE

# ---- 1. Class ----

test_that("sparse_otu_table has correct S4 class hierarchy", {
  ot <- otu_table(ps_s)
  expect_s4_class(ot, "sparse_otu_table")
  expect_s4_class(ot, "otu_table")
  expect_true(is(ot, "matrix"))
  expect_s4_class(ps_s, "phyloseq")
})

# ---- 2. Dimensions / names ----

test_that("dim matches original", {
  expect_identical(dim(otu_table(ps_s)), dim(otu_table(ps_d)))
  expect_identical(dim(otu_table(ps_st)), dim(otu_table(ps_dt)))
})

test_that("taxa_names / sample_names match original", {
  expect_identical(taxa_names(ps_s), taxa_names(ps_d))
  expect_identical(sample_names(ps_s), sample_names(ps_d))
  expect_identical(ntaxa(ps_s), ntaxa(ps_d))
  expect_identical(nsamples(ps_s), nsamples(ps_d))
})

test_that("taxa_are_rows preserved for both orientations", {
  expect_true(taxa_are_rows(otu_table(ps_s)))
  expect_false(taxa_are_rows(otu_table(ps_st)))
})

test_that("length() returns total element count", {
  expect_identical(length(otu_table(ps_s)), length(otu_table(ps_d)))
})

# ---- 3. Materialization ----

test_that("as(x, 'matrix') equals original matrix (values; sparse returns double, original integer)", {
  expect_equal(as(otu_table(ps_s), "matrix"), as(otu_table(ps_d), "matrix"))
  expect_equal(
    as(otu_table(ps_st), "matrix"),
    as(otu_table(ps_dt), "matrix")
  )
})

test_that("as.matrix() (S3) is identical to as(x, 'matrix')", {
  expect_identical(as.matrix(otu_table(ps_s)), as(otu_table(ps_s), "matrix"))
  expect_identical(as.matrix(otu_table(ps_st)), as(otu_table(ps_st), "matrix"))
})

test_that("as.data.frame() equals dense version (values; column type may differ double vs integer)", {
  expect_equal(
    as.data.frame(otu_table(ps_s)),
    as.data.frame(otu_table(ps_d))
  )
  expect_equal(
    as.data.frame(otu_table(ps_st)),
    as.data.frame(otu_table(ps_dt))
  )
})

test_that("as.data.frame() round-trips through as.matrix correctly", {
  df <- as.data.frame(otu_table(ps_s))
  m <- as(otu_table(ps_s), "matrix")
  expect_identical(df, as.data.frame(m))
})

# ---- 4. Sums (sparse-native path) ----

test_that("taxa_sums matches dense (tar = TRUE)", {
  expect_equal(taxa_sums(ps_s), taxa_sums(ps_d))
})

test_that("taxa_sums matches dense (tar = FALSE)", {
  expect_equal(taxa_sums(ps_st), taxa_sums(ps_dt))
})

test_that("sample_sums matches dense (tar = TRUE)", {
  expect_equal(sample_sums(ps_s), sample_sums(ps_d))
})

test_that("sample_sums matches dense (tar = FALSE)", {
  expect_equal(sample_sums(ps_st), sample_sums(ps_dt))
})

# ---- 5. Subsetting ----

test_that("[ row subset returns valid otu_table with correct values", {
  keep <- taxa_names(ps_s)[1:10]
  ot_sub <- otu_table(ps_s)[keep, ]
  expect_s4_class(ot_sub, "otu_table")
  expect_equal(nrow(ot_sub), 10L)
  expect_identical(rownames(ot_sub), keep)
  expect_equal(
    as(ot_sub, "matrix"),
    as(otu_table(ps_d)[keep, ], "matrix")
  )
})

test_that("[ column subset returns valid otu_table", {
  keep <- sample_names(ps_s)[1:5]
  ot_sub <- otu_table(ps_s)[, keep]
  expect_s4_class(ot_sub, "otu_table")
  expect_equal(ncol(ot_sub), 5L)
  expect_identical(colnames(ot_sub), keep)
})

test_that("[ two-index subset has correct dimensions and values", {
  kr <- taxa_names(ps_s)[1:10]
  kc <- sample_names(ps_s)[1:5]
  ot_sub <- otu_table(ps_s)[kr, kc]
  expect_s4_class(ot_sub, "otu_table")
  expect_equal(dim(ot_sub), c(10L, 5L))
  expect_equal(
    as(ot_sub, "matrix"),
    as(otu_table(ps_d)[kr, kc], "matrix")
  )
})

test_that("[ preserves taxa_are_rows flag in result", {
  keep <- taxa_names(ps_s)[1:5]
  ot_sub <- otu_table(ps_s)[keep, ]
  expect_true(taxa_are_rows(ot_sub))

  ot_sub_t <- otu_table(ps_st)[, keep] # tar=FALSE: samples are rows, taxa are cols
  expect_false(taxa_are_rows(ot_sub_t))
})

# ---- 5b. Integer index extraction ----

test_that("[ row subset by integer returns correct otu_table", {
  ot <- otu_table(ps_s)
  ot_sub <- ot[1:10, ]
  expect_s4_class(ot_sub, "otu_table")
  expect_equal(nrow(ot_sub), 10L)
  expect_identical(rownames(ot_sub), taxa_names(ps_s)[1:10])
  expect_equal(as(ot_sub, "matrix"), as(otu_table(ps_d)[1:10, ], "matrix"))
})

test_that("[ column subset by integer returns correct otu_table", {
  ot <- otu_table(ps_s)
  ot_sub <- ot[, 1:5]
  expect_s4_class(ot_sub, "otu_table")
  expect_equal(ncol(ot_sub), 5L)
  expect_identical(colnames(ot_sub), sample_names(ps_s)[1:5])
  expect_equal(as(ot_sub, "matrix"), as(otu_table(ps_d)[, 1:5], "matrix"))
})

test_that("[ two-index integer subset has correct dimensions and values", {
  ot <- otu_table(ps_s)
  ot_sub <- ot[1:10, 1:5]
  expect_s4_class(ot_sub, "otu_table")
  expect_equal(dim(ot_sub), c(10L, 5L))
  expect_equal(as(ot_sub, "matrix"), as(otu_table(ps_d)[1:10, 1:5], "matrix"))
})

test_that("[ single-row integer subset returns one-row otu_table", {
  ot <- otu_table(ps_s)
  ot_sub <- ot[3L, ]
  expect_s4_class(ot_sub, "otu_table")
  expect_equal(nrow(ot_sub), 1L)
  expect_identical(rownames(ot_sub), taxa_names(ps_s)[3L])
})

test_that("[ mixed character row and integer column subset", {
  ot <- otu_table(ps_s)
  row_names <- taxa_names(ps_s)[1:5]
  ot_sub <- ot[row_names, 1:3]
  expect_s4_class(ot_sub, "otu_table")
  expect_equal(dim(ot_sub), c(5L, 3L))
  expect_identical(rownames(ot_sub), row_names)
  expect_identical(colnames(ot_sub), sample_names(ps_s)[1:3])
})

# ---- 5c. Flat element extraction ----

test_that("is.na() returns correctly-sized FALSE matrix (not 0x0 from .Data stub)", {
  ot <- otu_table(ps_s)
  na_mat <- is.na(ot)
  expect_identical(dim(na_mat), dim(ot))
  expect_true(is.logical(na_mat))
  expect_false(any(na_mat)) # count data has no NAs
})

test_that("x[!is.na(x)] returns a plain vector of all values (flat extraction)", {
  ot <- otu_table(ps_s)
  vals <- ot[!is.na(ot)]
  expect_true(is.numeric(vals))
  expect_false(is.matrix(vals))
  expect_equal(sort(vals), sort(as(ot, "matrix")[!is.na(as(ot, "matrix"))]))
})

test_that("x[matrix_index] returns plain vector (flat extraction path)", {
  ot <- otu_table(ps_s)
  m <- as(ot, "matrix")
  idx <- m > 0 # correctly-sized logical matrix via dense path
  expect_equal(sort(ot[idx]), sort(m[idx]))
})

# ---- 6. Transpose ----

test_that("t() flips taxa_are_rows and swaps dimensions", {
  ot <- otu_table(ps_s)
  ott <- t(ot)
  expect_s4_class(ott, "sparse_otu_table")
  expect_false(taxa_are_rows(ott))
  expect_equal(dim(ott), rev(dim(ot)))
  expect_equal(as(ott, "matrix"), t(as(ot, "matrix")))
})

test_that("t(t(x)) round-trips to original", {
  ot <- otu_table(ps_s)
  expect_equal(as(t(t(ot)), "matrix"), as(ot, "matrix"))
  expect_true(taxa_are_rows(t(t(ot))))
})

# ---- 7. phyloseq integration: prune_taxa / prune_samples ----

test_that("prune_taxa produces valid phyloseq with correct taxa", {
  # Use taxa_sums > median to prune roughly half the taxa
  thresh <- median(taxa_sums(ps_s))
  ps_pruned <- prune_taxa(taxa_sums(ps_s) > thresh, ps_s)
  expect_s4_class(ps_pruned, "phyloseq")
  expect_lte(ntaxa(ps_pruned), ntaxa(ps_s))
  # Counts for kept taxa must match the original
  kept <- taxa_names(ps_pruned)
  expect_equal(
    as(otu_table(ps_pruned), "matrix"),
    as(otu_table(ps_d)[kept, ], "matrix")
  )
})

test_that("prune_taxa(taxa_sums > 0) keeps all taxa that have reads", {
  ps_pruned <- prune_taxa(taxa_sums(ps_s) > 0, ps_s)
  expect_s4_class(ps_pruned, "phyloseq")
  # All surviving taxa should have reads in dense version too
  expect_true(all(taxa_sums(ps_pruned) > 0))
})

test_that("prune_samples produces valid phyloseq", {
  ps_pruned <- prune_samples(sample_sums(ps_s) > 0, ps_s)
  expect_s4_class(ps_pruned, "phyloseq")
  expect_lte(nsamples(ps_pruned), nsamples(ps_s))
})

test_that("prune_taxa works with tar = FALSE orientation", {
  ps_pruned <- prune_taxa(taxa_sums(ps_st) > 0, ps_st)
  expect_s4_class(ps_pruned, "phyloseq")
})

# ---- 8. transform_sample_counts ----

test_that("transform_sample_counts (relative abundance) sums to 1 per sample", {
  ps_rel <- transform_sample_counts(ps_s, function(x) x / sum(x))
  expect_s4_class(ps_rel, "phyloseq")
  # Samples with zero reads yield NaN after 0/0; test only samples that had reads
  has_reads <- sample_names(ps_s)[sample_sums(ps_s) > 0]
  ss <- sample_sums(ps_rel)
  expect_equal(
    unname(ss[has_reads]),
    rep(1.0, length(has_reads)),
    tolerance = 1e-10
  )
})

test_that("transform_sample_counts result matches dense", {
  ps_rel_s <- transform_sample_counts(ps_s, function(x) x / sum(x))
  ps_rel_d <- transform_sample_counts(ps_d, function(x) x / sum(x))
  expect_equal(
    as(otu_table(ps_rel_s), "matrix"),
    as(otu_table(ps_rel_d), "matrix")
  )
})

# ---- 9. tax_glom (speedyseq) ----

test_that("tax_glom reduces number of taxa", {
  ps_glom <- tax_glom(ps_s, taxrank = "Genus")
  expect_s4_class(ps_glom, "phyloseq")
  expect_lt(ntaxa(ps_glom), ntaxa(ps_s))
})

test_that("tax_glom on sparse gives same taxon count as dense", {
  ps_glom_s <- tax_glom(ps_s, taxrank = "Genus")
  ps_glom_d <- speedyseq::tax_glom(ps_d, taxrank = "Genus")
  expect_equal(ntaxa(ps_glom_s), ntaxa(ps_glom_d))
})

test_that("tax_glom genus-level sums match (sparse vs dense)", {
  ps_glom_s <- tax_glom(ps_s, taxrank = "Genus")
  ps_glom_d <- speedyseq::tax_glom(ps_d, taxrank = "Genus")
  # Total read counts must be identical regardless of glom path
  expect_equal(sum(taxa_sums(ps_glom_s)), sum(taxa_sums(ps_glom_d)))
})

# ---- 10. rowSums / colSums ----

test_that("rowSums on sparse matches dense (tar = TRUE)", {
  expect_equal(rowSums(otu_table(ps_s)), rowSums(otu_table(ps_d)))
})

test_that("colSums on sparse matches dense (tar = TRUE)", {
  expect_equal(colSums(otu_table(ps_s)), colSums(otu_table(ps_d)))
})

test_that("rowSums on sparse matches dense (tar = FALSE)", {
  expect_equal(rowSums(otu_table(ps_st)), rowSums(otu_table(ps_dt)))
})

test_that("colSums on sparse matches dense (tar = FALSE)", {
  expect_equal(colSums(otu_table(ps_st)), colSums(otu_table(ps_dt)))
})

test_that("rowSums result is named", {
  rs <- rowSums(otu_table(ps_s))
  expect_named(rs)
  expect_identical(names(rs), taxa_names(ps_s))
})

test_that("colSums result is named", {
  cs <- colSums(otu_table(ps_s))
  expect_named(cs)
  expect_identical(names(cs), sample_names(ps_s))
})

test_that("rowSums falls back to base for plain matrix", {
  m <- matrix(1:6, 2, 3)
  expect_identical(rowSums(m), base::rowSums(m))
})

test_that("colSums falls back to base for plain matrix", {
  m <- matrix(1:6, 2, 3)
  expect_identical(colSums(m), base::colSums(m))
})

# ---- 11. Memory ----

test_that("sparse phyloseq uses less memory than dense for sparse data", {
  # Use a larger object for this check so the difference is reliable
  ps_d_big <- make_ps(n_taxa = 200, n_samples = 50, sparsity = 0.97)
  ps_s_big <- sparsify(ps_d_big)
  expect_lt(
    as.numeric(object.size(ps_s_big)),
    as.numeric(object.size(ps_d_big))
  )
})
