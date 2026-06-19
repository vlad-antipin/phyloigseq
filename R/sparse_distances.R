# NOTE: the right move is to keep phyloseq's default orientation (taxa x samples) with dgCMatrix,
# where each sample is a column. This is both cache-optimal for your distance kernel and maximally
# compatible with the ecosystem.
#
# WARNING: ps must not contain samples with zero total count (all-zero OTU vectors). Such samples
# produce denom == 0 in the Bray-Curtis formula (0+0), yielding NaN distances that cause
# downstream ordination (wcmdscale/eigen) to fail. Filter empty samples before calling these
# functions, or they are silently treated as distance 0 to each other via the ifelse guard below.

bray_curtis_sparse <- function(ps) {
  ot <- phyloseq::otu_table(ps)
  sp <- if (is(ot, "sparse_otu_table")) {
    ot@sparse_data
  } else {
    as(as(Matrix::Matrix(methods::as(ot, "matrix"), sparse = TRUE), "generalMatrix"), "CsparseMatrix")
  }
  if (!phyloseq::taxa_are_rows(ot)) {
    sp <- Matrix::t(sp)
  } # ensure taxa × samples

  n <- ncol(sp)
  ss <- Matrix::colSums(sp)
  nms <- colnames(sp)

  # t(sp) = samples × taxa as dgCMatrix; taxa are now columns → easy to iterate
  spt <- Matrix::t(sp)
  p <- spt@p
  ridx <- spt@i + 1L # 1-based sample indices
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

  denom <- outer(ss, ss, `+`)
  bc <- ifelse(denom == 0, 0.0, 1.0 - 2.0 * C / denom)
  diag(bc) <- 0
  dimnames(bc) <- list(nms, nms)
  stats::as.dist(bc)
}

#' Distance methods with native sparse implementations
#'
#' Character vector of distance method names that have a native sparse-matrix
#' implementation in \pkg{PhyloIgSeq}. Pass one of these to
#' \code{\link{sparse_distance}} to avoid converting the OTU table to a dense
#' matrix. Currently supported: \code{"bray"} (Bray-Curtis dissimilarity).
#'
#' @export
SPARSE_DISTANCE_METHODS <- c("bray")

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
#'   (e.g. \code{"bray"}).
#'
#' @return A \code{\link[stats]{dist}} object of pairwise sample distances.
#' @export
sparse_distance <- function(ps, method) {
  if (method %in% SPARSE_DISTANCE_METHODS) {
    switch(method,
      bray = bray_curtis_sparse(ps)
    )
  } else {
    warning("No sparse version for '", method, "', falling back to phyloseq::distance")
    phyloseq::distance(ps, method = method)
  }
}
