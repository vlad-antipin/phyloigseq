# Fast UniFrac (weighted/unweighted) for phyloseq objects, exploiting OTU-table
# sparsity. API mirrors `phyloseq::distance(physeq, method = ...)` for the
# "unifrac" / "wunifrac" methods and returns a `dist` object.
#
# Unlike phyloseq's own implementation, the per-edge descendant-abundance
# totals are computed with a single bottom-up (postorder) pass over the tree
# instead of re-summing the full descendant-tip set for every edge from
# scratch.
#
# Both methods compute every sample pair at once via matrix algebra (matrix
# multiplication for "unifrac", stats::dist() for "wunifrac") instead of an
# R-level loop over pairs -- see the comments in each branch below for the
# derivation. That's fast enough on its own that there's no parallel path.

# Descendant-tip / edge-abundance computation, shared by both the weighted
# and unweighted branches below.
#
# Computes edge_abundance[e, sample] = descendant-tip abundance total for
# the edge leading into tree$edge[e, 2] (i.e. the "child end" of edge e),
# via a single bottom-up (postorder) pass over the tree instead of
# re-summing the full descendant-tip set for every edge from scratch.
.unifrac_edge_abundance <- function(tree, count_mat) {
  n_tips <- length(tree$tip.label)
  n_nodes <- n_tips + tree$Nnode

  # Descendant tip set per node, built bottom-up (postorder: every edge's
  # child subtree is fully resolved before it's added to its parent). This
  # touches only tree topology, not sample data, so its cost is O(sum of
  # descendant-set sizes) -- O(n_tips log n_tips) for a reasonably balanced
  # tree -- independent of nsamples, unlike summing actual abundances edge
  # by edge would be. `index.only = TRUE` returns just the row-permutation
  # of `tree$edge` needed to visit edges child-before-parent, without
  # touching `tree` itself.
  postorder_edge_idx <- ape::reorder.phylo(
    tree,
    order = "postorder",
    index.only = TRUE
  )
  descendant_tips <- vector("list", n_nodes)
  for (tip in seq_len(n_tips)) {
    descendant_tips[[tip]] <- tip
  }
  for (i in postorder_edge_idx) {
    parent_node <- tree$edge[i, 1]
    child_node <- tree$edge[i, 2]
    descendant_tips[[parent_node]] <- c(
      descendant_tips[[parent_node]],
      descendant_tips[[child_node]]
    )
  }

  # ancestor_indicator[node, tip] = 1 if `tip` descends from `node` (a tip
  # descends from itself). Sparse: nnz is the same sum-of-descendant-set-
  # sizes total as above, a tiny fraction of n_nodes * n_tips for any
  # non-degenerate tree.
  descendant_counts <- lengths(descendant_tips)
  ancestor_indicator <- Matrix::sparseMatrix(
    i = rep(seq_len(n_nodes), descendant_counts),
    j = unlist(descendant_tips, use.names = FALSE),
    x = 1,
    dims = c(n_nodes, n_tips)
  )

  # subtree_abundance[node, sample] = total abundance of all tips descending
  # from `node`, for that sample. One sparse matrix multiplication replaces
  # the per-edge accumulation loop above (which, done directly on
  # abundances, would cost O(Nedge * nsamples) -- the dominant cost at
  # scale) and stays sparse -- count_mat is never densified, and neither is
  # this result. Most nodes are only "on" for a fraction of samples even
  # when a handful near the root are nonzero for nearly all of them, so
  # each branch below decides for itself, edge-group by edge-group, when
  # (if ever) densifying is worth it, instead of eagerly allocating a dense
  # n_nodes x nsamples matrix up front regardless of how it'll be used.
  subtree_abundance <- ancestor_indicator %*% count_mat
  colnames(subtree_abundance) <- colnames(count_mat)

  subtree_abundance[tree$edge[, 2], , drop = FALSE]
}

# Weighted UniFrac branch: abundance-weighted numerator/denominator computed
# for every sample pair at once.
.weighted_unifrac <- function(tree, count_mat, edge_abundance, sample_depths) {
  n_tips <- length(tree$tip.label)

  # Cumulative branch length from the root to each tip, used as the
  # per-tip weight in the weighted-UniFrac denominator (sum of branch
  # lengths times average tip abundance across the two samples).
  tip_root_distance <- ape::node.depth.edgelength(tree)
  tip_root_distance <- tip_root_distance[1:n_tips]
  names(tip_root_distance) <- tree$tip.label
  tip_root_distance <- tip_root_distance[rownames(count_mat)]

  inv_sample_depth <- Matrix::Diagonal(x = 1 / sample_depths)
  tip_proportions <- count_mat %*% inv_sample_depth # stays sparse
  colnames(tip_proportions) <- colnames(count_mat)
  edge_proportions <- edge_abundance %*% inv_sample_depth # stays sparse

  # Weighted UniFrac numerator for every sample pair at once: sum_e
  # edge_length_e * |p_ea - p_eb|. Unlike the presence/absence
  # intersection unweighted UniFrac uses below, this L1 (Manhattan) sum
  # has no bilinear form, so it can't be turned into a single BLAS matrix
  # multiplication -- stats::dist(method = "manhattan") is the direct way
  # to get it, but it's a plain single-threaded C loop with none of the
  # SIMD/multi-threaded speedup BLAS gets for the unweighted branch.
  #
  # Edges split into two groups, handled differently:
  #  - A handful near the tree's root are nonzero for close to every
  #    sample; for those a per-edge loop is no better than dist(), so
  #    they're densified (just this small subset, not the full edge set)
  #    and batched into one dist() call.
  #  - Most edges (leaf-ward ones) have small nnz in real sparse OTU data.
  #    For these, |p_ea - p_eb| = p_ea + p_eb - 2*min(p_ea, p_eb), and
  #    min(p_ea, p_eb) is 0 whenever *either* is 0 -- so the "one zero,
  #    one nonzero" pairs need no work at all, only pairs where both
  #    samples are nonzero on that edge contribute to the min term. That
  #    avoids ever materializing which samples are zero on a sparse edge
  #    (an O(nsamples) scan per edge that would otherwise dominate: with
  #    many sparse edges, O(Nedge_sparse * nsamples) reintroduces the
  #    exact cost this split was meant to avoid). The p_ea + p_eb part is
  #    just a per-sample column-sum over the sparse edges.
  weighted_edge <- edge_proportions * tree$edge.length
  n_samples <- ncol(weighted_edge)
  nnz_per_edge <- Matrix::rowSums(weighted_edge != 0)
  # 10% of nsamples is comfortably inside the plateau where the split
  # point barely matters (checked empirically from 5-15%); it just needs
  # to separate the long tail of near-full edges from the sparse bulk.
  dense_cutoff <- 0.1 * n_samples
  dense_edges <- which(nnz_per_edge > dense_cutoff)
  sparse_edges <- which(nnz_per_edge > 0 & nnz_per_edge <= dense_cutoff)

  numerator_mat <- matrix(0, n_samples, n_samples)
  if (length(dense_edges) > 0) {
    dense_block <- as.matrix(weighted_edge[dense_edges, , drop = FALSE])
    numerator_mat <- numerator_mat +
      as.matrix(stats::dist(t(dense_block), method = "manhattan"))
  }
  if (length(sparse_edges) > 0) {
    sparse_block <- weighted_edge[sparse_edges, , drop = FALSE]
    sparse_depth <- Matrix::colSums(sparse_block)
    numerator_mat <- numerator_mat + outer(sparse_depth, sparse_depth, "+")

    # One bulk triplet extraction (not one Matrix `[` call per edge, which
    # would pay S4 dispatch overhead O(Nedge_sparse) times), then a plain
    # R loop over edges using only base vectors from here on.
    triplet <- Matrix::summary(sparse_block)
    nz_by_row <- split(triplet$j, triplet$i)
    val_by_row <- split(triplet$x, triplet$i)
    min_sum <- matrix(0, n_samples, n_samples)
    for (local_row in seq_along(nz_by_row)) {
      nz <- nz_by_row[[local_row]]
      if (length(nz) > 1) {
        vals <- val_by_row[[local_row]]
        min_sum[nz, nz] <- min_sum[nz, nz] + outer(vals, vals, pmin)
      }
    }
    numerator_mat <- numerator_mat - 2 * min_sum
  }
  dimnames(numerator_mat) <- list(
    colnames(weighted_edge),
    colnames(weighted_edge)
  )
  numerator <- stats::as.dist(numerator_mat)

  # Denominator for every pair at once: sum_t tip_root_distance_t *
  # (tip_prop_ta + tip_prop_tb) splits into a per-sample scalar
  # depth[s] = sum_t tip_root_distance_t * tip_prop_ts, so
  # denominator[a, b] = depth[a] + depth[b] -- no per-pair tip-level work
  # needed either.
  tip_weighted_depth <- as.vector(tip_root_distance %*% tip_proportions)
  names(tip_weighted_depth) <- colnames(tip_proportions)
  denominator <- stats::as.dist(
    outer(tip_weighted_depth, tip_weighted_depth, "+")
  )

  numerator / denominator
}

# Unweighted UniFrac branch: presence/absence numerator/denominator computed
# for every sample pair at once, via matrix multiplication instead of a
# per-pair intersect/union loop. With present_e,s the 0/1 indicator of
# whether edge e has nonzero descendant abundance in sample s:
#   shared[a, b] = sum_e edge_length_e * present_e,a * present_e,b
#     (branch length shared by both samples)
#   total[s]     = shared[s, s] = sum_e edge_length_e * present_e,s
#     (branch length present in sample s alone)
#   union[a, b]  = total[a] + total[b] - shared[a, b]
# shared is exactly the cross-product of the (branch-length-scaled)
# presence matrix with itself, so one matrix multiplication produces every
# pairwise numerator at once; union then only needs the diagonal of shared
# plus an outer sum, both O(nsamples^2) but allocation-only (no tree walk).
.unweighted_unifrac <- function(tree, edge_abundance) {
  # That cross-product can be done as a sparse Matrix `%*%` (skips zero
  # entries, one S4 call rather than one per pair so no per-pair dispatch
  # tax) or as a plain dense `%*%` (routed to BLAS). edge_present starts
  # sparse (edge_abundance was never densified above); a bulk sparse
  # multiply wins once the matrix is both large and sparse enough that
  # skipping zeros outweighs its CSC bookkeeping overhead, but BLAS's raw
  # throughput wins for small/moderately-sparse inputs, so it's worth
  # densifying just for those. Checked empirically: dense wins under ~5M
  # edge*sample cells or ~15% edge-presence density; past that, sparse pulls
  # ahead (e.g. 2.7x faster on a 16000-edge x 2000-sample, 2.6%-dense case).
  edge_present <- edge_abundance > 0
  edge_cells <- nrow(edge_present) * ncol(edge_present)
  edge_density <- Matrix::nnzero(edge_present) / edge_cells
  if (!(edge_cells > 5e6 && edge_density < 0.15)) {
    edge_present <- as.matrix(edge_present)
  }
  scaled_edge_present <- edge_present * tree$edge.length
  shared_length <- as.matrix(t(scaled_edge_present) %*% edge_present)
  total_length <- diag(shared_length)
  union_length <- outer(total_length, total_length, "+") - shared_length
  stats::as.dist(1 - shared_length / union_length)
}

#' Fast UniFrac distances via a sparse OTU table
#'
#' Computes unweighted or weighted UniFrac distances between samples,
#' exploiting OTU-table sparsity for speed. Results are numerically
#' equivalent to \code{\link[phyloseq]{distance}(physeq, method = method)},
#' but per-edge descendant-abundance totals are computed with a single
#' bottom-up (postorder) pass over the tree instead of re-summing the full
#' descendant-tip set for every edge from scratch. Idea inspired from
#' \url{https://github.com/joey711/phyloseq/issues/524}.
#'
#' If the tree in \code{physeq} is unrooted, it is midpoint-rooted here,
#' whereas \code{\link[phyloseq]{UniFrac}}/\code{\link[phyloseq]{distance}}
#' root it at a random tip instead; the two are therefore only guaranteed to
#' agree when both are run on the same already-rooted tree, not when handed
#' the same unrooted tree.
#'
#' @param physeq A \code{\link[phyloseq]{phyloseq}} object with a
#'   \code{\link[phyloseq]{phy_tree}} slot. The OTU table may be a
#'   \code{\link{sparse_otu_table-class}} or a standard dense
#'   \code{\link[phyloseq]{otu_table}}. If the tree is unrooted, it is
#'   midpoint-rooted (with a warning) before distances are computed.
#' @param method Either \code{"unifrac"} (unweighted, presence/absence) or
#'   \code{"wunifrac"} (abundance-weighted).
#'
#' @return A \code{\link[stats]{dist}} object of pairwise UniFrac distances
#'   between samples, or \code{NULL} with a warning if \code{physeq} has no
#'   \code{phy_tree}.
#'
#' @examples
#' data(ps_16s_refinement)
#' ps_sparse <- as_sparse_phyloseq(ps_16s_refinement)
#'
#' # Tree is unrooted, so this midpoint-roots it first (with a warning)
#' d_unweighted <- sparse_unifrac(ps_sparse, method = "unifrac")
#' class(d_unweighted)
#'
#' d_weighted <- sparse_unifrac(ps_sparse, method = "wunifrac")
#'
#' # Also works directly on a dense (non-sparse) phyloseq object
#' d_dense <- sparse_unifrac(ps_16s_refinement, method = "unifrac")
#'
#' @export
sparse_unifrac <- function(physeq, method = c("unifrac", "wunifrac")) {
  method <- match.arg(method)
  weighted <- method == "wunifrac"

  tree <- phy_tree(physeq, errorIfNULL = FALSE)
  if (is.null(tree)) {
    warning("phy_tree slot is empty, UniFrac requires a tree")
    return(NULL)
  }

  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }

  if (!ape::is.rooted(tree)) {
    tree <- phytools::midpoint_root(tree)
    warning("Tree is unrooted, midpoint was set as root")
  }

  # Pull the count matrix out as a sparse (taxa x samples) CsparseMatrix,
  # whether the OTU table is already our sparse S4 class or a dense phyloseq
  # otu_table.
  otu_tab <- otu_table(physeq)
  count_mat <- if (methods::is(otu_tab, "sparse_otu_table")) {
    otu_tab@sparse_data
  } else {
    methods::as(
      Matrix::Matrix(as(otu_tab, "matrix"), sparse = TRUE),
      "CsparseMatrix"
    )
  }

  sample_names_vec <- sample_names(physeq)
  sample_depths <- sample_sums(physeq)
  # Neither branch below builds a sample-pairs list (both are fully
  # vectorized), so this reproduces the "n < m" error phyloseq::distance()
  # gets for free from combn() when there are fewer than 2 samples.
  if (length(sample_names_vec) < 2) {
    stop("n < m")
  }

  # Match OTU names in OTU table to tree edges: the bottom-up pass below
  # seeds tip rows by position, so the count-matrix row order must line up
  # with `tree$tip.label` order.
  if (!all(rownames(count_mat) == taxa_names(tree))) {
    count_mat <- count_mat[taxa_names(tree), ]
  }

  edge_abundance <- .unifrac_edge_abundance(tree, count_mat)

  if (weighted) {
    return(.weighted_unifrac(tree, count_mat, edge_abundance, sample_depths))
  }

  .unweighted_unifrac(tree, edge_abundance)
}
