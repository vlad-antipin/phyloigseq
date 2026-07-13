# NOTE: the right move is to keep phyloseq's default orientation (taxa x samples) with dgCMatrix,
# where each sample is a column. This is both cache-optimal for your distance kernel and maximally
# compatible with the ecosystem.
#
# WARNING: ps must not contain samples with zero total count (all-zero OTU vectors). Such samples
# produce denom == 0 in several formulas, yielding NaN or degenerate distances that can cause
# downstream ordination (wcmdscale/eigen) to fail. Filter empty samples before calling these
# functions, or they are silently treated as distance 0 to each other via the ifelse guard below.

# ============================================================
# Internal helpers
# ============================================================

get_sp_taxa_by_samples <- function(ps) {
  ot <- phyloseq::otu_table(ps)
  sp <- if (is(ot, "sparse_otu_table")) {
    ot@sparse_data
  } else {
    as(
      as(
        Matrix::Matrix(methods::as(ot, "matrix"), sparse = TRUE),
        "generalMatrix"
      ),
      "CsparseMatrix"
    )
  }
  if (!phyloseq::taxa_are_rows(ot)) {
    sp <- Matrix::t(sp)
  }
  sp
}

# C[i,j] = sum_k min(x_ki, x_kj), where k indexes taxa and i,j index samples.
# Iterates over the CSC structure of t(sp): each column is a taxon, non-zero entries
# are the samples that carry it.  Only pairs of samples sharing the same taxon contribute.
compute_min_matrix <- function(sp) {
  n <- ncol(sp)
  spt <- Matrix::t(sp)
  p <- spt@p
  ridx <- spt@i + 1L
  rval <- spt@x
  C <- matrix(0.0, n, n)
  for (k in seq_len(ncol(spt))) {
    s <- p[k] + 1L
    e <- p[k + 1L]
    if (e < s + 1L) {
      next
    } # 0 or 1 non-zero sample: no pairwise contribution
    si <- ridx[s:e]
    sv <- rval[s:e]
    C[si, si] <- C[si, si] + outer(sv, sv, pmin)
  }
  C
}

# Shared tail for every *_sparse kernel below: zero the diagonal, attach
# sample-name dimnames, and coerce to a dist object.
finalize_sparse_dist <- function(mat, nms) {
  diag(mat) <- 0
  dimnames(mat) <- list(nms, nms)
  stats::as.dist(mat)
}

# ============================================================
# Distance implementations (not exported; called via sparse_distance)
# ============================================================

bray_curtis_sparse <- function(ps) {
  sp <- get_sp_taxa_by_samples(ps)
  ss <- Matrix::colSums(sp)
  nms <- colnames(sp)
  C <- compute_min_matrix(sp)
  denom <- outer(ss, ss, `+`)
  bc <- ifelse(denom == 0, 0.0, 1.0 - 2.0 * C / denom)
  finalize_sparse_dist(bc, nms)
}

# Abundance-weighted Jaccard: 1 - C / (A + B - C)
# Equivalent to vegan's "jaccard" for quantitative data (Ruzicka index).
jaccard_sparse <- function(ps) {
  sp <- get_sp_taxa_by_samples(ps)
  ss <- Matrix::colSums(sp)
  nms <- colnames(sp)
  C <- compute_min_matrix(sp)
  n <- length(ss)
  A <- matrix(ss, n, n, byrow = FALSE)
  B <- matrix(ss, n, n, byrow = TRUE)
  denom <- A + B - C
  jac <- ifelse(denom == 0, 0.0, 1.0 - C / denom)
  finalize_sparse_dist(jac, nms)
}

# Kulczynski: 1 - 0.5 * (C/A + C/B)
kulczynski_sparse <- function(ps) {
  sp <- get_sp_taxa_by_samples(ps)
  ss <- Matrix::colSums(sp)
  nms <- colnames(sp)
  C <- compute_min_matrix(sp)
  n <- length(ss)
  A <- matrix(ss, n, n, byrow = FALSE)
  B <- matrix(ss, n, n, byrow = TRUE)
  sim <- 0.5 * (ifelse(A == 0, 0, C / A) + ifelse(B == 0, 0, C / B))
  kul <- 1 - sim
  finalize_sparse_dist(kul, nms)
}

# L1 (Manhattan) distance for non-negative data: A + B - 2*C
# Valid because |xi - yi| = (xi + yi) - 2*min(xi, yi) for xi, yi >= 0.
manhattan_sparse <- function(ps) {
  sp <- get_sp_taxa_by_samples(ps)
  ss <- Matrix::colSums(sp)
  nms <- colnames(sp)
  C <- compute_min_matrix(sp)
  man <- outer(ss, ss, `+`) - 2 * C
  man[man < 0] <- 0 # floating-point guard
  finalize_sparse_dist(man, nms)
}

# Euclidean via the dot-product identity: ||xi - xj||^2 = ||xi||^2 + ||xj||^2 - 2*(xi.xj)
euclidean_sparse <- function(ps) {
  sp <- get_sp_taxa_by_samples(ps)
  nms <- colnames(sp)
  D <- as.matrix(Matrix::crossprod(sp)) # n x n dot-product matrix
  norm2 <- diag(D)
  n <- ncol(sp)
  eu_sq <- outer(norm2, rep(1, n)) + outer(rep(1, n), norm2) - 2 * D
  eu_sq[eu_sq < 0] <- 0 # floating-point guard
  eu <- sqrt(eu_sq)
  finalize_sparse_dist(eu, nms)
}

# Canberra: (1/n_active) * sum_{union(xi,xj)>0} |xi-yi|/(xi+yi)
# For non-negative data and the k-th taxon:
#   - xi>0, xj=0: term = 1
#   - xi=0, xj>0: term = 1
#   - xi>0, xj>0: term = 1 - 2*min(xi,xj)/(xi+xj)
# Summed: n_union - 2 * sum_{both>0} min(xi,xj)/(xi+xj)
# Normalized by n_union = nnz_i + nnz_j - n_both (union of non-zero taxa).
canberra_sparse <- function(ps) {
  sp <- get_sp_taxa_by_samples(ps)
  n <- ncol(sp)
  nms <- colnames(sp)

  spt <- Matrix::t(sp)
  p <- spt@p
  ridx <- spt@i + 1L
  rval <- spt@x

  n_both <- matrix(0L, n, n)
  harm_min <- matrix(0.0, n, n)

  for (k in seq_len(ncol(spt))) {
    s <- p[k] + 1L
    e <- p[k + 1L]
    if (e < s) {
      next
    } # taxon absent from all samples
    si <- ridx[s:e]
    sv <- rval[s:e]
    n_both[si, si] <- n_both[si, si] + 1L
    harm_min[si, si] <- harm_min[si, si] +
      outer(sv, sv, function(a, b) pmin(a, b) / (a + b))
  }

  nnz <- Matrix::colSums(sp != 0) # non-zero taxa count per sample
  n_union <- outer(nnz, rep(1, n)) + outer(rep(1, n), nnz) - n_both

  can <- ifelse(n_union == 0, 0.0, 1.0 - 2 * harm_min / n_union)
  finalize_sparse_dist(can, nms)
}

# Morisita-Horn: 1 - 2*(xi.xj) / (da*N1*N2 + db*N1*N2)
# da = sum(xi^2)/N1^2, db = sum(xj^2)/N2^2
# Denominator simplifies to outer(sq/ss, ss) + outer(ss, sq/ss).
horn_sparse <- function(ps) {
  sp <- get_sp_taxa_by_samples(ps)
  nms <- colnames(sp)
  ss <- Matrix::colSums(sp)
  sq <- Matrix::colSums(sp^2)
  D <- as.matrix(Matrix::crossprod(sp))
  sqss <- ifelse(ss == 0, 0, sq / ss)
  denom <- outer(sqss, ss) + outer(ss, sqss)
  horn <- ifelse(denom == 0, 0.0, 1.0 - 2 * D / denom)
  finalize_sparse_dist(horn, nms)
}

# Chord: Euclidean distance after L2-normalizing each sample.
# chord[i,j] = sqrt(2 - 2*(xi.xj) / (||xi|| * ||xj||))
chord_sparse <- function(ps) {
  sp <- get_sp_taxa_by_samples(ps)
  nms <- colnames(sp)
  D <- as.matrix(Matrix::crossprod(sp))
  norm2 <- diag(D)
  denom <- outer(sqrt(norm2), sqrt(norm2))
  inner <- ifelse(denom == 0, 0, 2 - 2 * D / denom)
  inner[inner < 0] <- 0 # floating-point guard
  ch <- sqrt(inner)
  finalize_sparse_dist(ch, nms)
}

# Hellinger: Euclidean distance after Hellinger transformation (sqrt of relative abundances).
# hellinger[i,j] = sqrt(2 - 2 * DH[i,j] / sqrt(ss_i * ss_j))
# where DH[i,j] = sum_k sqrt(x_ki) * sqrt(x_kj)  (sparse: zero entries stay zero)
hellinger_sparse <- function(ps) {
  sp <- get_sp_taxa_by_samples(ps)
  nms <- colnames(sp)
  ss <- Matrix::colSums(sp)
  sp_sqrt <- sp
  sp_sqrt@x <- sqrt(sp_sqrt@x)
  DH <- as.matrix(Matrix::crossprod(sp_sqrt))
  denom <- outer(sqrt(ss), sqrt(ss))
  inner <- ifelse(denom == 0, 0, 2 - 2 * DH / denom)
  inner[inner < 0] <- 0 # floating-point guard
  hel <- sqrt(inner)
  finalize_sparse_dist(hel, nms)
}

# ============================================================
# Exported interface
# ============================================================

#' Distance methods with native sparse implementations
#'
#' Character vector of distance method names that have a native sparse-matrix
#' implementation in \pkg{PhyloIgSeq}. Pass one of these to
#' \code{\link{sparse_distance}} to avoid converting the OTU table to a dense
#' matrix.
#'
#' @examples
#' SPARSE_DISTANCE_METHODS
#' "bray" %in% SPARSE_DISTANCE_METHODS
#'
#' @export
SPARSE_DISTANCE_METHODS <- c(
  "bray",
  "jaccard",
  "kulczynski",
  "manhattan",
  "euclidean",
  "canberra",
  "horn",
  "chord",
  "hellinger",
  "unifrac",
  "wunifrac"
)

#' Compute pairwise sample distances with sparse OTU table support
#'
#' Dispatches to a native sparse-matrix implementation when \code{method}
#' appears in \code{\link{SPARSE_DISTANCE_METHODS}}, and falls back to
#' \code{\link[phyloseq]{distance}} otherwise.
#'
#' @param ps A \code{\link[phyloseq]{phyloseq}} object. The OTU table may be
#'   a \code{\link{sparse_otu_table-class}} or a standard dense
#'   \code{\link[phyloseq]{otu_table}}.
#' @param method A single character string naming the distance metric
#'   (e.g. \code{"bray"}, \code{"jaccard"}, \code{"kulczynski"},
#'   \code{"manhattan"}, \code{"euclidean"}, \code{"canberra"}, \code{"horn"},
#'   \code{"chord"}, \code{"hellinger"}).
#'
#' @return A \code{\link[stats]{dist}} object of pairwise sample distances.
#'
#' @examples
#' data(ps_16s_refinement)
#' ps_sparse <- as_sparse_phyloseq(ps_16s_refinement)
#' bc <- sparse_distance(ps_sparse, method = "bray")
#' class(bc)
#'
#' # Falls back to phyloseq::distance() (with a warning) for methods without
#' # a native sparse implementation.
#' gower <- sparse_distance(ps_sparse, method = "gower")
#'
#' @export
sparse_distance <- function(ps, method) {
  if (method %in% SPARSE_DISTANCE_METHODS) {
    switch(method,
      bray = bray_curtis_sparse(ps),
      jaccard = jaccard_sparse(ps),
      kulczynski = kulczynski_sparse(ps),
      manhattan = manhattan_sparse(ps),
      euclidean = euclidean_sparse(ps),
      canberra = canberra_sparse(ps),
      horn = horn_sparse(ps),
      chord = chord_sparse(ps),
      hellinger = hellinger_sparse(ps),
      unifrac = sparse_unifrac(ps, "unifrac"),
      wunifrac = sparse_unifrac(ps, "wunifrac")
    )
  } else {
    warning(
      "No sparse version for '",
      method,
      "', falling back to phyloseq::distance"
    )
    phyloseq::distance(ps, method = method)
  }
}
