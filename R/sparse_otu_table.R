# Clean up stale sparse_otu_table method definitions
suppressMessages({
  for (.m in c("initialize", "dim", "dimnames", "length", "[", "t")) {
    if (existsMethod(.m, "sparse_otu_table")) {
      removeMethod(.m, "sparse_otu_table")
    }
  }
  if (isClass("sparse_otu_table")) {
    removeClass("sparse_otu_table")
  }
})

#' Sparse OTU table backed by a dgCMatrix
#'
#' An S4 class extending \code{\link[phyloseq]{otu_table}} that stores
#' abundance counts in a sparse matrix (\code{dgCMatrix}) instead of a dense R
#' matrix, reducing memory use for typical microbiome datasets where most
#' taxon-sample counts are zero. All standard phyloseq operations
#' (\code{prune_taxa}, \code{sample_sums}, subsetting with \code{[}, etc.)
#' work transparently.
#'
#' @details
#' Low-level coercion via \code{as(x, "matrix")} / \code{as(x, "data.frame")}
#' is also defined and is what \code{\link{as.matrix.sparse_otu_table}} /
#' \code{\link{as.data.frame.sparse_otu_table}} call into; most users will
#' prefer the \code{as.matrix()}/\code{as.data.frame()} spelling below.
#'
#' @slot sparse_data A \code{dgCMatrix} holding the actual abundance counts.
#'   The \code{.Data} slot inherited from \code{otu_table} is kept as a 0x0
#'   stub to satisfy S4 validity without allocating a full dense copy.
#'
#' @name sparse_otu_table-class
#' @aliases sparse_otu_table-class
#' @exportClass sparse_otu_table
setClass(
  "sparse_otu_table",
  contains = "otu_table",
  representation(sparse_data = "dgCMatrix")
)

# otu_table validity reads .Data dimensions directly, bypassing our dim override,
# and fails on an empty stub. Skip validObject entirely — phyloseq's own validity
# check (which runs when otu_table(ps) <- ...) uses taxa_names() which DOES
# dispatch through our dimnames override and sees the correct names.
setMethod(
  "initialize",
  "sparse_otu_table",
  function(.Object, .Data, taxa_are_rows, sparse_data) {
    .Object@.Data <- .Data
    .Object@taxa_are_rows <- taxa_are_rows
    .Object@sparse_data <- sparse_data
    .Object # no callNextMethod() → no validObject
  }
)

#' Create a sparse OTU table
#'
#' Converts a \code{\link[phyloseq]{otu_table}} into a
#' \code{\link{sparse_otu_table-class}} backed by a \code{dgCMatrix}, reducing
#' memory use when most counts are zero.
#'
#' @param otu An \code{\link[phyloseq]{otu_table}} object.
#'
#' @return A \code{\link{sparse_otu_table-class}} object.
#'
#' @seealso \code{\link{as_sparse_phyloseq}} to convert a whole
#'   \code{phyloseq} object at once.
#' @examples
#' data(ps_16s_refinement)
#' ot_sparse <- sparse_otu_table(phyloseq::otu_table(ps_16s_refinement))
#' class(ot_sparse)
#' dim(ot_sparse)
#'
#' @export
sparse_otu_table <- function(otu) {
  tar <- phyloseq::taxa_are_rows(otu)
  sp <- as(
    as(Matrix::Matrix(as(otu, "matrix"), sparse = TRUE), "generalMatrix"),
    "CsparseMatrix"
  )
  new(
    "sparse_otu_table",
    .Data = matrix(integer(0), 0L, 0L),
    taxa_are_rows = tar,
    sparse_data = sp
  )
}

# dim + dimnames are the root: rownames/colnames/nrow/ncol/ntaxa/nsamples
# all derive from these via base R and phyloseq generics
#' @rdname sparse_otu_table-class
#' @param x A \code{sparse_otu_table} object.
setMethod("dim", "sparse_otu_table", function(x) dim(x@sparse_data))

#' @rdname sparse_otu_table-class
setMethod("dimnames", "sparse_otu_table", function(x) dimnames(x@sparse_data))

#' @rdname sparse_otu_table-class
setMethod("length", "sparse_otu_table", function(x) length(x@sparse_data))

#' @rdname sparse_otu_table-class
#' @param value Replacement dimnames (a list of two character vectors).
setReplaceMethod("dimnames", "sparse_otu_table", function(x, value) {
  dimnames(x@sparse_data) <- value
  x
})

# is.na on the 0x0 .Data stub returns a 0x0 matrix, which later causes
# "subscript too long" in code that does x[!is.na(x)].  Override to read
# sparse_data instead — for count OTU data all entries are non-NA.
#' @rdname sparse_otu_table-class
setMethod("is.na", "sparse_otu_table", function(x) {
  is.na(as(x@sparse_data, "matrix"))
})

#' @rdname sparse_otu_table-class
#' @param i Row index (integer, character, or logical).
#' @param j Column index (integer, character, or logical).
#' @param ... Unused; required by the generic signature.
#' @param drop Ignored; included for S4 generic compatibility.
setMethod("[", "sparse_otu_table", function(x, i, j, ..., drop = FALSE) {
  if (!missing(i) && is.character(i)) {
    i <- match(i, rownames(x@sparse_data))
  }
  if (!missing(j) && is.character(j)) {
    j <- match(j, colnames(x@sparse_data))
  }
  if (!missing(i) && missing(j)) {
    nr <- nrow(x@sparse_data)
    # Flat element extraction: x[!is.na(x)], x[x > 0], x[matrix_index].
    # A matrix index or a logical vector longer than nrow means the caller
    # wants a plain vector of selected values, not a row-subset otu_table.
    if (is.matrix(i) || (is.logical(i) && length(i) != nr)) {
      return(as(x, "matrix")[i])
    }
  }
  sp <- if (missing(i) && missing(j)) {
    x@sparse_data
  } else if (missing(i)) {
    x@sparse_data[, j, drop = FALSE]
  } else if (missing(j)) {
    x@sparse_data[i, , drop = FALSE]
  } else {
    x@sparse_data[i, j, drop = FALSE]
  }
  # Must return otu_table (not bare matrix): prune_taxa does
  # otu_table(ps) <- x[taxa, ], which requires otu_tableOrNULL
  phyloseq::otu_table(as(sp, "matrix"), taxa_are_rows = x@taxa_are_rows)
})

setAs("sparse_otu_table", "matrix", function(from) {
  as(from@sparse_data, "matrix")
})

setAs("sparse_otu_table", "data.frame", function(from) {
  as.data.frame(as(from, "matrix"))
})

# S3 methods so as.matrix() / as.data.frame() route here instead of
# *.default, which checks is.matrix() at C-level, sees the 0x0 .Data stub,
# and returns the sparse_otu_table unchanged (silently wrong).
#' @rdname sparse_otu_table-class
#' @exportS3Method base::as.matrix
as.matrix.sparse_otu_table <- function(x, ...) as(x, "matrix")
#' @rdname sparse_otu_table-class
#' @exportS3Method base::as.data.frame
as.data.frame.sparse_otu_table <- function(x, ...) {
  as.data.frame(as(x, "matrix"), ...)
}

#' @rdname sparse_otu_table-class
setMethod("t", "sparse_otu_table", function(x) {
  new(
    "sparse_otu_table",
    matrix(integer(0), 0L, 0L),
    taxa_are_rows = !x@taxa_are_rows,
    sparse_data = t(x@sparse_data)
  )
})

# Plain functions in .GlobalEnv shadow phyloseq's taxa_sums/sample_sums for
# unqualified user-level calls (bench::mark, direct console calls, other
# packages that search the global path). phyloseq's own namespace is NOT
# patched, so an explicitly-qualified phyloseq::taxa_sums()/sample_sums()
# call still runs phyloseq's dense-only implementation and will silently
# return wrong results (or error) on a sparse_otu_table — always call
# PhyloIgSeq::taxa_sums()/PhyloIgSeq::sample_sums() explicitly from code
# that isn't relying on search-path shadowing (see CLAUDE.md's app
# integration note).
#' Sparse-aware taxa and sample sum functions
#'
#' Drop-in replacements for \code{\link[phyloseq]{taxa_sums}} and
#' \code{\link[phyloseq]{sample_sums}} that read directly from the sparse
#' \code{dgCMatrix} slot when the OTU table is a
#' \code{\link{sparse_otu_table-class}}, avoiding materialisation of the full
#' dense matrix. Falls back to the phyloseq implementations for standard OTU
#' tables. For an \code{\link{incomplete_otu_table-class}} (structurally
#' missing, not just zero, entries), the sum is taken only over the
#' \emph{observed} taxon-sample pairs and a warning is issued, since
#' unobserved entries cannot be assumed to be zero.
#'
#' @param physeq A \code{\link[phyloseq]{phyloseq}} object or an
#'   \code{\link[phyloseq]{otu_table}}.
#' @param ... Additional arguments passed to \code{phyloseq::taxa_sums} or
#'   \code{phyloseq::sample_sums} for non-sparse tables.
#'
#' @return A named numeric vector of per-taxon (\code{taxa_sums}) or
#'   per-sample (\code{sample_sums}) abundance sums.
#' @examples
#' data(ps_16s_refinement)
#' ps_sparse <- as_sparse_phyloseq(ps_16s_refinement)
#' head(taxa_sums(ps_sparse))
#' head(sample_sums(ps_sparse))
#' @name sparse-sums
NULL

# Shared 3-way dispatch for taxa_sums/sample_sums: incomplete_otu_table (warn,
# then sum only the observed entries), sparse_otu_table (direct sparse sum),
# or fall back to the phyloseq implementation for standard OTU tables.
.sparse_margin_sums <- function(physeq, want_taxa, warn_label, phyloseq_fn, ...) {
  ot <- if (is(physeq, "phyloseq")) phyloseq::otu_table(physeq) else physeq
  if (is(ot, "incomplete_otu_table")) {
    warning(
      warn_label, " on an incomplete_otu_table sums only observed Ig score ",
      "entries; unobserved taxon-sample pairs are excluded."
    )
  } else if (!is(ot, "sparse_otu_table")) {
    return(phyloseq_fn(physeq, ...))
  }
  use_row <- phyloseq::taxa_are_rows(ot) == want_taxa
  if (use_row) Matrix::rowSums(ot@sparse_data) else Matrix::colSums(ot@sparse_data)
}

#' @rdname sparse-sums
#' @export
taxa_sums <- function(physeq, ...) {
  .sparse_margin_sums(
    physeq,
    want_taxa = TRUE,
    warn_label = "taxa_sums",
    phyloseq_fn = phyloseq::taxa_sums,
    ...
  )
}
#' @rdname sparse-sums
#' @export
sample_sums <- function(physeq, ...) {
  .sparse_margin_sums(
    physeq,
    want_taxa = FALSE,
    warn_label = "sample_sums",
    phyloseq_fn = phyloseq::sample_sums,
    ...
  )
}

# rowSums / colSums: plain .GlobalEnv functions shadow base:: for all user-level
# calls (base::rowSums reads .Data stub directly, returning wrong 0x0 results).
# Shared 3-way dispatch, mirroring .sparse_margin_sums above but operating
# directly on the otu_table/matrix rather than a phyloseq wrapper.
.sparse_matrix_sums <- function(x, na.rm, dims, sparse_fn, base_fn, warn_label, ...) {
  if (is(x, "incomplete_otu_table")) {
    warning(
      warn_label, " on an incomplete_otu_table sums only observed Ig score ",
      "entries; unobserved taxon-sample pairs are excluded."
    )
    return(sparse_fn(x@sparse_data, na.rm = na.rm))
  }
  if (is(x, "sparse_otu_table")) {
    return(sparse_fn(x@sparse_data, na.rm = na.rm))
  }
  base_fn(x, na.rm = na.rm, dims = dims, ...)
}

rowSums <- function(x, na.rm = FALSE, dims = 1L, ...) {
  .sparse_matrix_sums(x, na.rm, dims, Matrix::rowSums, base::rowSums, "rowSums", ...)
}

colSums <- function(x, na.rm = FALSE, dims = 1L, ...) {
  .sparse_matrix_sums(x, na.rm, dims, Matrix::colSums, base::colSums, "colSums", ...)
}

# Convert a phyloseq object to use a sparse-backed OTU table.
# Safe to call repeatedly: returns the object unchanged if the OTU table is
# already a sparse_otu_table.  All other phyloseq components (tax_table,
# sample_data, phy_tree, refseq) are preserved.
#' Convert a phyloseq object to use a sparse OTU table
#'
#' Replaces the OTU table of a \code{\link[phyloseq]{phyloseq}} object with a
#' \code{\link{sparse_otu_table-class}} backed by a \code{dgCMatrix}. Safe to
#' call repeatedly: returns the object unchanged if the OTU table is already
#' sparse. All other phyloseq components (tax table, sample data, phylogenetic
#' tree, reference sequences) are preserved.
#'
#' @param ps A \code{\link[phyloseq]{phyloseq}} object.
#'
#' @return The same \code{phyloseq} object with its OTU table replaced by a
#'   \code{\link{sparse_otu_table-class}}.
#'
#' @seealso \code{\link{sparse_otu_table}} to convert an OTU table directly.
#' @examples
#' data(ps_16s_refinement)
#' ps_sparse <- as_sparse_phyloseq(ps_16s_refinement)
#' class(phyloseq::otu_table(ps_sparse))
#' @export
as_sparse_phyloseq <- function(ps) {
  stopifnot(is(ps, "phyloseq"))
  ot <- phyloseq::otu_table(ps)
  if (is(ot, "sparse_otu_table")) {
    return(ps)
  }
  ps@otu_table <- sparse_otu_table(ot)
  ps
}

# ── incomplete_otu_table ───────────────────────────────────────────────────────

# Clean up stale incomplete_otu_table method definitions before (re-)loading
suppressMessages({
  for (.m in c("initialize", "is.na", "[")) {
    if (existsMethod(.m, "incomplete_otu_table")) {
      removeMethod(.m, "incomplete_otu_table")
    }
  }
  for (.to in c("matrix", "data.frame")) {
    if (existsMethod("coerce", c("incomplete_otu_table", .to))) {
      removeMethod("coerce", c("incomplete_otu_table", .to))
    }
  }
  if (isClass("incomplete_otu_table")) {
    removeClass("incomplete_otu_table")
  }
})

#' Incomplete OTU table backed by a softImpute matrix factorization
#'
#' An S4 class extending \code{\link{sparse_otu_table-class}} for IgA-Seq Ig
#' score matrices where most entries are \code{NA} (taxon not observed in that
#' sample — genuinely missing, not zero). Observed entries are stored as an
#' \code{Incomplete} object (from the \code{softImpute} package, which extends
#' \code{dgCMatrix}), and a low-rank SVD factorization is fitted eagerly so
#' that sample-to-sample distances and taxon loadings can be computed without
#' reconstructing the full dense matrix.
#'
#' @slot svd_fit A named list with elements \code{$u} (n_samples × r),
#'   \code{$d} (length-r singular values), and \code{$v} (n_taxa × r) from
#'   \code{\link[softImpute]{softImpute}}.
#' @slot col_means Named numeric vector (length n_taxa) of per-taxon means
#'   computed from observed entries before centering. \code{numeric(0)} when
#'   centering was skipped.
#'
#' @name incomplete_otu_table-class
#' @aliases incomplete_otu_table-class
#' @exportClass incomplete_otu_table
setClass(
  "incomplete_otu_table",
  contains = "sparse_otu_table",
  representation(svd_fit = "list", col_means = "numeric")
)

setMethod(
  "initialize",
  "incomplete_otu_table",
  function(.Object, .Data, taxa_are_rows, sparse_data, svd_fit, col_means) {
    .Object@.Data <- .Data
    .Object@taxa_are_rows <- taxa_are_rows
    .Object@sparse_data <- sparse_data
    .Object@svd_fit <- svd_fit
    .Object@col_means <- col_means
    .Object
  }
)

#' Create an incomplete OTU table with an embedded SVD factorization
#'
#' Wraps a \code{softImpute::Incomplete} object and a fitted SVD factorization
#' into an \code{\link{incomplete_otu_table-class}}. In normal use this is
#' called from \code{PhyloIgSeq_to_phyloseq(imputation_method = "SVD")} rather
#' than directly.
#'
#' @param X_inc An \code{Incomplete} object (from
#'   \code{\link[softImpute]{Incomplete}}) with \code{dimnames} set to
#'   \code{list(sample_names, taxa_names)}.
#' @param svd_fit A list with elements \code{$u}, \code{$d}, \code{$v} as
#'   returned by \code{\link[softImpute]{softImpute}}.
#' @param col_means Named numeric vector of per-taxon column means used for
#'   centering before fitting; \code{numeric(0)} if centering was not applied.
#' @param taxa_are_rows Logical; \code{FALSE} (default) means rows = samples,
#'   columns = taxa, which matches the layout of the \code{Incomplete} object
#'   built from \code{ig_coating} triplets.
#'
#' @return An \code{\link{incomplete_otu_table-class}} object.
#' @export
incomplete_otu_table <- function(
  X_inc,
  svd_fit,
  col_means = numeric(0),
  taxa_are_rows = FALSE
) {
  new(
    "incomplete_otu_table",
    .Data = matrix(integer(0), 0L, 0L),
    taxa_are_rows = taxa_are_rows,
    sparse_data = X_inc,
    svd_fit = svd_fit,
    col_means = col_means
  )
}

# as.matrix: NA-filled dense matrix where only the positions in the
# dgCMatrix sparsity pattern (the observed entries) are filled in.
# This differs from the parent's setAs which converts structural zeros
# (= unobserved positions in Incomplete) to 0.
setAs("incomplete_otu_table", "matrix", function(from) {
  sp <- from@sparse_data
  dims <- dim(sp)
  mat <- matrix(
    NA_real_,
    nrow = dims[1L],
    ncol = dims[2L],
    dimnames = dimnames(sp)
  )
  if (length(sp@x) > 0L) {
    col_idx <- rep(seq_len(dims[2L]), diff(sp@p))
    row_idx <- sp@i + 1L # 0-based → 1-based
    mat[cbind(row_idx, col_idx)] <- sp@x
  }
  mat
})

setAs("incomplete_otu_table", "data.frame", function(from) {
  as.data.frame(as(from, "matrix"))
})

#' @rdname incomplete_otu_table-class
#' @exportS3Method base::as.matrix
as.matrix.incomplete_otu_table <- function(x, ...) as(x, "matrix")

#' @rdname incomplete_otu_table-class
#' @exportS3Method base::as.data.frame
as.data.frame.incomplete_otu_table <- function(x, ...) {
  as.data.frame(as(x, "matrix"), ...)
}

# is.na: TRUE for every unobserved position (structural zero in the dgCMatrix),
# FALSE for every observed position — the inverse of standard dgCMatrix
# semantics where stored zeros are "non-missing".
setMethod("is.na", "incomplete_otu_table", function(x) {
  sp <- x@sparse_data
  dims <- dim(sp)
  mat <- matrix(TRUE, nrow = dims[1L], ncol = dims[2L], dimnames = dimnames(sp))
  if (length(sp@x) > 0L) {
    col_idx <- rep(seq_len(dims[2L]), diff(sp@p))
    row_idx <- sp@i + 1L
    mat[cbind(row_idx, col_idx)] <- FALSE
  }
  mat
})

# [ subset: materialise to NA-filled dense matrix, then return a standard
# otu_table.  The SVD fit is not meaningful after arbitrary subsetting, so
# it is not carried over.
setMethod("[", "incomplete_otu_table", function(x, i, j, ..., drop = FALSE) {
  if (!missing(i) && is.character(i)) {
    i <- match(i, rownames(x@sparse_data))
  }
  if (!missing(j) && is.character(j)) {
    j <- match(j, colnames(x@sparse_data))
  }
  if (!missing(i) && missing(j)) {
    nr <- nrow(x@sparse_data)
    if (is.matrix(i) || (is.logical(i) && length(i) != nr)) {
      return(as(x, "matrix")[i])
    }
  }
  mat <- as(x, "matrix")
  sub <- if (missing(i) && missing(j)) {
    mat
  } else if (missing(i)) {
    mat[, j, drop = FALSE]
  } else if (missing(j)) {
    mat[i, , drop = FALSE]
  } else {
    mat[i, j, drop = FALSE]
  }
  phyloseq::otu_table(sub, taxa_are_rows = x@taxa_are_rows)
})

# speedyseq::merge_taxa_vec calls phyloseq::taxa_sums(otu_table(x)) internally,
# which falls through to base::rowSums — a C routine that reads .Data directly
# and errors on the 0x0 stub. Materialise to dense before delegating.
#' Taxonomic agglomeration with sparse OTU table support
#'
#' A wrapper around \code{\link[speedyseq]{tax_glom}} that temporarily
#' materialises the OTU table to a dense matrix when the input uses a
#' \code{\link{sparse_otu_table-class}}, then delegates to
#' \code{speedyseq::tax_glom}. This is necessary because
#' \code{speedyseq::merge_taxa_vec} calls \code{phyloseq::taxa_sums}
#' internally, which falls through to a C-level routine that reads the stub
#' \code{.Data} slot directly and errors on the 0x0 placeholder.
#'
#' @param physeq A \code{\link[phyloseq]{phyloseq}} object.
#' @param taxrank A single character string naming the taxonomic rank to
#'   agglomerate at (must be a column of the tax table).
#' @param ... Additional arguments passed to \code{speedyseq::tax_glom}.
#'
#' @return A \code{\link[phyloseq]{phyloseq}} object agglomerated at
#'   \code{taxrank}, with a standard dense OTU table.
#' @export
tax_glom <- function(physeq, taxrank, ...) {
  if (
    is(physeq, "phyloseq") &&
      is(phyloseq::otu_table(physeq), "sparse_otu_table")
  ) {
    ot <- phyloseq::otu_table(physeq)
    phyloseq::otu_table(physeq) <- phyloseq::otu_table(
      as(ot, "matrix"),
      taxa_are_rows = phyloseq::taxa_are_rows(ot)
    )
  }
  speedyseq::tax_glom(physeq, taxrank = taxrank, ...)
}
