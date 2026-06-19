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
    methods::as(Matrix::Matrix(methods::as(ot, "matrix"), sparse = TRUE), "dgCMatrix")
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

#' @export
sparse_distance = function(ps, method) {
  if (method == "bray") {
    dist.matrix = bray_curtis_sparse(ps)
  } else {
    warning("No sparse version for this distance metric")
    dist.matrix = phyloseq::distance(ps, method = method)
  }
  dist.matrix
}
