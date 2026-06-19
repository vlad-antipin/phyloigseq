# Clean up stale sparse_otu_table method definitions
suppressMessages({
  for (.m in c("initialize", "dim", "dimnames", "length", "[", "t")) {
    if (existsMethod(.m, "sparse_otu_table")) removeMethod(.m, "sparse_otu_table")
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
#' @export
sparse_otu_table <- function(otu) {
  tar <- phyloseq::taxa_are_rows(otu)
  sp <- as(as(Matrix::Matrix(as(otu, "matrix"), sparse = TRUE), "generalMatrix"), "CsparseMatrix")
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
#' @exportS3Method base::as.matrix
as.matrix.sparse_otu_table <- function(x, ...) as(x, "matrix")
#' @exportS3Method base::as.data.frame
as.data.frame.sparse_otu_table <- function(x, ...) {
  as.data.frame(as(x, "matrix"), ...)
}

setMethod("t", "sparse_otu_table", function(x) {
  new(
    "sparse_otu_table",
    matrix(integer(0), 0L, 0L),
    taxa_are_rows = !x@taxa_are_rows,
    sparse_data = t(x@sparse_data)
  )
})

# Plain functions in .GlobalEnv shadow phyloseq's taxa_sums/sample_sums for
# user-level calls (bench::mark, direct console calls, other packages that
# search the global path).  phyloseq's own copies are also patched below via
# assignInNamespace so that explicit phyloseq::sample_sums() calls from third-
# party packages (e.g. PhyloIgSeq) go through the sparse-aware path too.
#' Sparse-aware taxa and sample sum functions
#'
#' Drop-in replacements for \code{\link[phyloseq]{taxa_sums}} and
#' \code{\link[phyloseq]{sample_sums}} that read directly from the sparse
#' \code{dgCMatrix} slot when the OTU table is a
#' \code{\link{sparse_otu_table-class}}, avoiding materialisation of the full
#' dense matrix. Falls back to the phyloseq implementations for standard OTU
#' tables.
#'
#' @param physeq A \code{\link[phyloseq]{phyloseq}} object or an
#'   \code{\link[phyloseq]{otu_table}}.
#' @param ... Additional arguments passed to \code{phyloseq::taxa_sums} or
#'   \code{phyloseq::sample_sums} for non-sparse tables.
#'
#' @return A named numeric vector of per-taxon (\code{taxa_sums}) or
#'   per-sample (\code{sample_sums}) abundance sums.
#' @name sparse-sums
NULL

#' @rdname sparse-sums
#' @export
taxa_sums <- function(physeq, ...) {
  ot <- if (is(physeq, "phyloseq")) phyloseq::otu_table(physeq) else physeq
  if (is(ot, "sparse_otu_table")) {
    sp <- ot@sparse_data
    if (phyloseq::taxa_are_rows(ot)) {
      Matrix::rowSums(sp)
    } else {
      Matrix::colSums(sp)
    }
  } else {
    phyloseq::taxa_sums(physeq, ...)
  }
}
#' @rdname sparse-sums
#' @export
sample_sums <- function(physeq, ...) {
  ot <- if (is(physeq, "phyloseq")) phyloseq::otu_table(physeq) else physeq
  if (is(ot, "sparse_otu_table")) {
    sp <- ot@sparse_data
    if (phyloseq::taxa_are_rows(ot)) {
      Matrix::colSums(sp)
    } else {
      Matrix::rowSums(sp)
    }
  } else {
    phyloseq::sample_sums(physeq, ...)
  }
}

# rowSums / colSums: plain .GlobalEnv functions shadow base:: for all user-level
# calls (base::rowSums reads .Data stub directly, returning wrong 0x0 results).
rowSums <- function(x, na.rm = FALSE, dims = 1L, ...) {
  if (is(x, "sparse_otu_table")) {
    Matrix::rowSums(x@sparse_data, na.rm = na.rm)
  } else {
    base::rowSums(x, na.rm = na.rm, dims = dims, ...)
  }
}

colSums <- function(x, na.rm = FALSE, dims = 1L, ...) {
  if (is(x, "sparse_otu_table")) {
    Matrix::colSums(x@sparse_data, na.rm = na.rm)
  } else {
    base::colSums(x, na.rm = na.rm, dims = dims, ...)
  }
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
  if (is(physeq, "phyloseq") && is(phyloseq::otu_table(physeq), "sparse_otu_table")) {
    ot <- phyloseq::otu_table(physeq)
    phyloseq::otu_table(physeq) <- phyloseq::otu_table(
      as(ot, "matrix"),
      taxa_are_rows = phyloseq::taxa_are_rows(ot)
    )
  }
  speedyseq::tax_glom(physeq, taxrank = taxrank, ...)
}
