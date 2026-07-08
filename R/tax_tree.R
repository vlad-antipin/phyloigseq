#' Taxonomic distance from longest common prefix of shared ranks
#'
#' Computes a weighted taxonomy distance matrix using the
#' longest-common-prefix (shared ancestry) principle: two taxa are similar
#' if they share many consecutive taxonomic ranks from the top (Kingdom ->
#' Phylum -> ...), and higher ranks receive larger weights (sharing Kingdom
#' and Phylum but differing at Class is more informative than sharing only
#' Kingdom).
#'
#' `distance = 1 - (shared_weight / total_weight)`, where `shared_weight` is
#' the sum of weights for the longest run of ranks shared from Kingdom
#' downward -- once a rank differs, all deeper ranks are treated as
#' unshared regardless of whether they coincidentally match. A distance of
#' 1 means the two taxa already differ at the highest rank. Taxon names are
#' appended as an implicit finest rank, so distance between two different
#' taxa is never exactly 0 even when their annotated taxonomy is identical
#' -- two different taxa are known to be different sequences (otherwise
#' they would already have been merged upstream), so collapsing them to
#' distance 0 would discard that information.
#'
#' Requires an n x n comparison (`NOTE: very slow`), so it does not scale to
#' large datasets -- use [build_taxonomy_tree_hierarchy()] instead for those.
#'
#' @param phyloseq_object A phyloseq object with a taxonomy table.
#' @param rank_weights Optional numeric vector of weights, one per rank
#'   (including the appended taxon-name rank), highest rank first. Defaults
#'   to `rank_weight_base^((number_of_ranks - 1):0)`, normalized to sum to 1.
#' @param rank_weight_base Base of the default exponential rank-weight
#'   decay (see [get_taxonomy_tree()] for its effect). Ignored if
#'   `rank_weights` is supplied.
#' @return A `dist` object with one entry per taxon pair.
#' @export
build_taxonomy_distance_longest_prefix <- function(
  phyloseq_object,
  rank_weights = NULL,
  rank_weight_base = 1.1
) {
  # make_unique_taxa_table() prevents a real correctness bug, not just an
  # NA-comparison nuisance: without it, two taxa that are both NA (or share
  # any other repeated placeholder label) at some rank, but already differ
  # at a coarser rank, would spuriously "match" at that rank and understate
  # their distance. Called here (not just in get_taxonomy_tree()) so this
  # function is correct even when called directly.
  taxonomy_matrix <- as.matrix(PhyloIgSeq::make_unique_taxa_table(
    phyloseq_object@tax_table
  ))
  taxonomy_matrix <- cbind(
    taxonomy_matrix,
    ASV = taxa_names(phyloseq_object)
  )

  number_of_taxa <- nrow(taxonomy_matrix)
  number_of_ranks <- ncol(taxonomy_matrix)

  # Default weights:
  # Kingdom gets the highest weight,
  # the appended ASV rank gets the lowest weight.
  if (is.null(rank_weights)) {
    rank_weights <- rank_weight_base^((number_of_ranks - 1):0)
  }

  # Normalize so distances lie between 0 and 1
  rank_weights <- rank_weights / sum(rank_weights)

  # Preallocate distance matrix
  distance_matrix <- matrix(
    0,
    nrow = number_of_taxa,
    ncol = number_of_taxa,
    dimnames = list(
      rownames(taxonomy_matrix),
      rownames(taxonomy_matrix)
    )
  )

  # Compute only upper triangle
  for (taxon_i in seq_len(number_of_taxa)) {
    for (taxon_j in taxon_i:number_of_taxa) {
      # Determine which ranks are identical
      ranks_match <- taxonomy_matrix[taxon_i, ] == taxonomy_matrix[taxon_j, ]

      # Longest common prefix:
      # once a rank differs, deeper ranks are ignored
      shared_prefix <- cumprod(ranks_match)

      # Total weight of shared ancestry
      shared_weight <- sum(shared_prefix * rank_weights)

      # Convert similarity into distance
      taxonomy_distance <- 1 - shared_weight

      distance_matrix[taxon_i, taxon_j] <- taxonomy_distance
      distance_matrix[taxon_j, taxon_i] <- taxonomy_distance
    }
  }

  as.dist(distance_matrix)
}

#' Taxonomic distance from rankwise shared weight
#'
#' Computes a heuristic taxonomic distance between taxa by comparing them
#' rank-by-rank (Kingdom, Phylum, ..., plus an implicit finest rank made of
#' the taxon names) and summing the weight of every rank where they match,
#' regardless of whether shallower ranks also matched. Higher ranks get
#' larger weights, since differences at broad taxonomic levels represent
#' larger taxonomic divergence.
#'
#' `distance = 1 - similarity`, where
#' `similarity = sum(weight of ranks where taxonomy is identical)`. As with
#' [build_taxonomy_distance_longest_prefix()], appending taxon names as
#' an implicit finest rank guarantees distance is never exactly 0 between
#' two different taxa. After [make_unique_taxa_table()] has disambiguated
#' repeated labels across lineages, this produces exactly the same result
#' as [build_taxonomy_distance_longest_prefix()] (a rankwise match can
#' never occur past a real point of divergence), but is vectorized via a
#' sparse one-hot encoding + crossprod rather than an explicit double loop.
#'
#' Requires an n x n similarity matrix, so it does not scale to large
#' datasets -- use [build_taxonomy_tree_hierarchy()] instead for those. This
#' function raises an informative error (rather than crashing) if the
#' sparse crossprod would overflow.
#'
#' @inheritParams build_taxonomy_distance_longest_prefix
#' @return A `dist` object with one entry per taxon pair.
#' @export
build_taxonomy_distance_rankwise <- function(
  phyloseq_object,
  rank_weights = NULL,
  rank_weight_base = 1.1
) {
  # See build_taxonomy_distance_longest_prefix() for why
  # make_unique_taxa_table() is required here (prevents NA/placeholder
  # labels from creating false rank matches across different lineages)
  # and why taxon names are appended as the finest rank (guarantees a
  # nonzero distance between taxa that are known to be different
  # sequences even when their annotated taxonomy is identical).
  taxonomy_matrix <- as.matrix(PhyloIgSeq::make_unique_taxa_table(
    phyloseq_object@tax_table
  ))
  taxonomy_matrix <- cbind(
    taxonomy_matrix,
    ASV = taxa_names(phyloseq_object)
  )

  # Taxon names (rows of taxonomy table)
  taxon_names <- rownames(taxonomy_matrix)

  # Ensure taxa names exist and are unique
  if (is.null(taxon_names)) {
    stop("Taxonomy table has no row names (taxon names).")
  }

  if (anyDuplicated(taxon_names)) {
    stop("Taxon names are not unique.")
  }

  number_of_taxa <- nrow(taxonomy_matrix)
  number_of_ranks <- ncol(taxonomy_matrix)

  # Default decreasing weights:
  # first ranks (Kingdom/Phylum) contribute more
  if (is.null(rank_weights)) {
    rank_weights <- rank_weight_base^((number_of_ranks - 1):0)
  }

  # Normalize weights to sum to 1
  rank_weights <- rank_weights / sum(rank_weights)

  # Guard against Matrix::tcrossprod() silently overflowing (this is what
  # produces the "bus error: invalid alignment" crash on large datasets).
  # A rank where one level covers a large fraction of all taxa turns that
  # rank's one-hot block into a near-dense n x n contribution once crossed
  # with itself, and dgCMatrix uses 32-bit nonzero indices internally, so
  # once the *sum* of (group size)^2 across ranks passes 2^31-1 the sparse
  # multiply can corrupt memory instead of raising a normal R error.
  estimated_nnz <- 0
  for (rank_index in seq_len(number_of_ranks)) {
    current_rank <- taxonomy_matrix[, rank_index, drop = TRUE]
    group_sizes <- tabulate(match(current_rank, unique(current_rank)))
    estimated_nnz <- estimated_nnz + sum(as.double(group_sizes)^2)
  }
  if (estimated_nnz > .Machine$integer.max) {
    stop(
      "build_taxonomy_distance_rankwise() would need an estimated ",
      format(estimated_nnz, scientific = TRUE),
      " nonzero entries in an ",
      "intermediate sparse matrix, past the 32-bit limit Matrix::tcrossprod() ",
      "can safely index -- this is the cause of the bus-error crash on large ",
      "datasets. Use get_taxonomy_tree(ps, method = \"hierarchy\") instead: it ",
      "builds the tree directly from the taxonomy's rank structure without ",
      "ever forming an n x n matrix, so it has no such limit."
    )
  }

  # Build a sparse one-hot encoding of all ranks at once, with each
  # rank's block of columns scaled by sqrt(weight). Then
  # encoded_matrix %*% t(encoded_matrix) sums, in a single sparse
  # matrix multiplication, exactly the weighted rank-match totals
  # that the previous version computed via one dense n x n outer()
  # comparison per rank (O(ranks * n^2) R-level work). This avoids
  # materializing `number_of_ranks` dense n x n logical matrices.
  block_list <- vector("list", number_of_ranks)
  for (rank_index in seq_len(number_of_ranks)) {
    # Extract rank as a vector
    current_rank <- as.vector(taxonomy_matrix[, rank_index, drop = TRUE])

    rank_levels <- unique(current_rank)
    level_codes <- match(current_rank, rank_levels)

    block_list[[rank_index]] <- Matrix::sparseMatrix(
      i = seq_len(number_of_taxa),
      j = level_codes,
      x = sqrt(rank_weights[rank_index]),
      dims = c(number_of_taxa, length(rank_levels))
    )
  }

  encoded_matrix <- do.call(cbind, block_list)
  similarity_matrix <- as.matrix(Matrix::tcrossprod(encoded_matrix))
  dimnames(similarity_matrix) <- list(taxon_names, taxon_names)

  # Convert similarity into distance
  distance_matrix <- 1 - similarity_matrix

  # Restore taxon names
  dimnames(distance_matrix) <- list(
    taxon_names,
    taxon_names
  )

  # Return a dist object compatible with clustering/ordination
  return(as.dist(distance_matrix))
}

#' Build a phylogenetic tree directly from the taxonomy hierarchy
#'
#' Builds a tree directly from the nested rank structure of a taxonomy
#' table (Kingdom -> Phylum -> ... -> finest rank -> taxon name), without
#' ever forming a pairwise distance matrix. This is the scalable
#' alternative to [build_taxonomy_distance_rankwise()] /
#' [build_taxonomy_distance_longest_prefix()] + `hclust()`, which need an
#' n x n matrix and so are O(n^2) in both time and memory (and can
#' overflow `Matrix`'s 32-bit sparse indices for large n). This walks each
#' rank once to assign every taxon an integer "clade code" identifying its
#' lineage at that depth (vectorized paste+match, no per-row recursion or
#' subsetting), then builds the tree's Newick representation bottom-up
#' from those codes: `O(n * number_of_ranks)`, with no quadratic case.
#'
#' Each rank (including the appended taxon-name rank) contributes a branch
#' length equal to its (normalized) weight, so the patristic (tip-to-tip
#' tree path) distance between two taxa equals `2 *` the
#' [build_taxonomy_distance_rankwise()] heuristic distance between them
#' -- the factor of 2 is just the usual up-then-down tree-path duplication
#' and is irrelevant for (w)UniFrac, which normalizes by total tree length.
#'
#' The result is resolved into a strictly bifurcating rooted tree (via
#' [ape::collapse.singles()] then [ape::multi2di()]), since
#' `phyloseq::UniFrac()` requires one and fails with an opaque error
#' otherwise; neither step changes any tip-to-tip patristic distance.
#'
#' @inheritParams build_taxonomy_distance_longest_prefix
#' @return A rooted, bifurcating object of class `phylo`.
#' @export
build_taxonomy_tree_hierarchy <- function(
  phyloseq_object,
  rank_weights = NULL,
  rank_weight_base = 1.1
) {
  # See build_taxonomy_distance_longest_prefix() for why
  # make_unique_taxa_table() is required (prevents NA/placeholder labels
  # from creating false shared-ancestry across different lineages) and why
  # taxon names are appended as the finest rank (taxa are known to be
  # different sequences even when their annotated taxonomy is identical,
  # so they should not collapse to distance 0).
  taxonomy_matrix <- as.matrix(PhyloIgSeq::make_unique_taxa_table(
    phyloseq_object@tax_table
  ))
  taxon_names <- taxa_names(phyloseq_object)
  taxonomy_matrix <- cbind(taxonomy_matrix, ASV = taxon_names)

  number_of_taxa <- nrow(taxonomy_matrix)
  number_of_ranks <- ncol(taxonomy_matrix)

  if (is.null(rank_weights)) {
    rank_weights <- rank_weight_base^((number_of_ranks - 1):0)
  }
  rank_weights <- rank_weights / sum(rank_weights)

  # For each rank depth, assign an integer code identifying the clade
  # (full lineage from the root down to that rank) each taxon belongs to,
  # plus a lookup from that code to its parent's (one rank coarser) code.
  clade_codes <- matrix(0L, nrow = number_of_taxa, ncol = number_of_ranks)
  parent_lookup <- vector("list", number_of_ranks)
  parent_code <- rep(1L, number_of_taxa) # implicit single root above rank 1
  for (rank_index in seq_len(number_of_ranks)) {
    lineage_key <- paste(parent_code, taxonomy_matrix[, rank_index], sep = "\r")
    codes <- match(lineage_key, unique(lineage_key))

    first_occurrence <- !duplicated(codes)
    lookup <- integer(max(codes))
    lookup[codes[first_occurrence]] <- parent_code[first_occurrence]
    parent_lookup[[rank_index]] <- lookup

    clade_codes[, rank_index] <- codes
    parent_code <- codes
  }

  # Bottom-up: start from tip labels carrying their own pendant length (the
  # ASV pseudo-rank's weight -- ASV values are unique, so this level is
  # always a singleton "group" and needs no wrapping node of its own),
  # then repeatedly wrap each clade at the current depth into a Newick
  # subtree with that rank's branch length, moving to the parent depth
  # using the lookup above. Note: wrapping a bare label in parens and
  # appending a length attaches that length to the *enclosing* node, not
  # to the label itself, which is why the ASV level is handled separately
  # here instead of going through the same loop.
  current_labels <- paste0(taxon_names, ":", rank_weights[number_of_ranks])
  current_codes <- clade_codes[, number_of_ranks - 1]

  for (rank_index in (number_of_ranks - 1):1) {
    grouped <- split(current_labels, current_codes)
    subtree_strings <- vapply(
      grouped,
      function(members) paste0("(", paste(members, collapse = ","), ")"),
      character(1)
    )
    group_codes <- as.integer(names(grouped))
    subtree_strings <- paste0(subtree_strings, ":", rank_weights[rank_index])

    current_labels <- unname(subtree_strings)
    current_codes <- parent_lookup[[rank_index]][group_codes]
  }

  newick <- paste0("(", paste(current_labels, collapse = ","), ");")
  taxonomy_tree <- ape::read.tree(text = newick)

  # Every rank contributes its own node even where it doesn't actually
  # branch, so the tree is full of single-child ("knuckle") nodes; on top
  # of that, ranks with >2 children (or several taxa sharing the exact
  # same non-ASV lineage) are genuine polytomies. phyloseq::UniFrac()/
  # fastUniFrac() require a strictly bifurcating tree and fail with an
  # opaque "subscript out of bounds" otherwise, so collapse the knuckle
  # nodes (summing their branch lengths into the parent edge) and then
  # resolve remaining polytomies with zero-length binary splits. Neither
  # step changes any tip-to-tip patristic distance.
  taxonomy_tree <- ape::collapse.singles(taxonomy_tree)
  taxonomy_tree <- ape::multi2di(taxonomy_tree)

  stopifnot(
    setequal(taxonomy_tree$tip.label, taxon_names)
  )

  taxonomy_tree
}

#' Get a phylogenetic tree from a taxonomy table
#'
#' Dispatches to one of the taxonomy-based tree-building methods and
#' returns a tree whose tip labels match `taxa_names(ps)`.
#'
#' @param ps A phyloseq object with a taxonomy table.
#' @param method `"hierarchy"` (default) builds the tree directly from the
#'   taxonomy's rank structure via [build_taxonomy_tree_hierarchy()] -- no
#'   pairwise distance matrix, no `hclust()`, scales to arbitrarily large
#'   datasets. `"rankwise"` and `"longest_prefix"` go through an n x n
#'   distance matrix ([build_taxonomy_distance_rankwise()] /
#'   [build_taxonomy_distance_longest_prefix()]) + UPGMA instead; they are
#'   kept for reference/small datasets but do not scale, and `"rankwise"`
#'   raises an informative error rather than crashing if the dataset is too
#'   large for its sparse-matrix approach.
#' @param rank_weight_base Controls how quickly rank weights decay from
#'   Kingdom (highest) to the finest rank (lowest): weight at depth `d`
#'   (from the deepest) is `rank_weight_base^d`, normalized to sum to 1. A
#'   base close to 1 weights all ranks near-evenly; a larger base makes
#'   coarse-rank differences dominate the distance. Passed through to
#'   whichever `method` is selected.
#' @return An object of class `phylo`.
#' @export
get_taxonomy_tree <- function(
  ps,
  method = c("hierarchy", "rankwise", "longest_prefix"),
  rank_weight_base = 1.1
) {
  method <- match.arg(method)

  taxonomy_tree <- switch(
    method,
    hierarchy = build_taxonomy_tree_hierarchy(
      ps,
      rank_weight_base = rank_weight_base
    ),
    rankwise = {
      taxonomy_distance <- build_taxonomy_distance_rankwise(
        ps,
        rank_weight_base = rank_weight_base
      )
      taxonomy_hclust <- hclust(
        taxonomy_distance,
        method = "average" # UPGMA; appropriate for distance-based taxonomy
      )
      tree <- ape::as.phylo(taxonomy_hclust)
      tree$tip.label <- labels(taxonomy_distance)
      tree
    },
    longest_prefix = {
      taxonomy_distance <- build_taxonomy_distance_longest_prefix(
        ps,
        rank_weight_base = rank_weight_base
      )
      taxonomy_hclust <- hclust(taxonomy_distance, method = "average")
      tree <- ape::as.phylo(taxonomy_hclust)
      tree$tip.label <- labels(taxonomy_distance)
      tree
    }
  )

  # Check compatibility
  stopifnot(
    setequal(taxonomy_tree$tip.label, taxa_names(ps))
  )

  return(taxonomy_tree)
}
