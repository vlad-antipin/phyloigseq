setClassUnion("data.frameOrNULL", c("data.frame", "NULL"))
setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))
#' PhyloIgSeq class
#'
#' An S4 class to represent the results of an IgSeq experiment, including Ig coating scores,
#' fractions, and taxonomic or sample-level metadata.
#'
#' @details
#' \code{ig_coating} always carries the identifier columns \code{taxon_id}/\code{sample_id},
#' followed by the actual Ig score columns (named in \code{score_names}), followed by the
#' per-fraction abundance columns and \code{zeros_imputed}. Only the columns listed in
#' \code{score_names} are Ig-coating scores; use \code{\link{get_ig_score}} to pull one out, and
#' \code{\link{ig_fraction_names}} to enumerate the abundance columns, rather than inferring
#' either from what is left over.
#'
#' The MA-plot geometry a sliding Z-score is derived from lives in its own \code{ma_coords} slot
#' rather than alongside the scores, because there can be more than one change axis per run (see
#' \code{\link{get_ma_coordinates}}) and they would otherwise collide on one set of
#' \code{obs_change}/\code{null_change} column names. That slot is long in
#' \code{change_transform}: one block of rows per axis, so adding an axis costs rows rather than
#' columns. \code{\link{ma_coords}} is the accessor.
#'
#' @slot ig_coating A data.frame containing per-taxon/sample Ig scores plus supporting metadata
#'   (see Details)
#' @slot score_names Character. Names of the \code{ig_coating} columns that are actual Ig scores
#'   (as opposed to fraction/diagnostic metadata columns)
#' @slot ma_coords A data.frame or NULL. Per-taxon/sample MA-plot geometry underlying the sliding
#'   Z-scores, long in \code{change_transform}: \code{sample_id}, \code{taxon_id},
#'   \code{change_transform}, \code{obs_abundance}, \code{obs_change}, \code{null_abundance},
#'   \code{null_change}, \code{ellipse_level}, \code{obs_in_cone}, \code{null_in_cone}
#' @slot positive_fraction_name Character. Name(s) of the positive Ig-coated fraction(s) used
#'   to build this object. A single value in the common case; when [getPhyloIgSeq()] was given
#'   more than one, this holds the full set folded into the sample dimension (see
#'   [getPhyloIgSeq()]'s `positive_fraction_name` documentation) -- informational only, not an
#'   indexable per-row value (use the `sample_data` column of the same name for that).
#' @slot first_negative_fraction_name Character or NULL. Name of the main negative fraction (e.g., 90%)
#' @slot second_negative_fraction_name Character or NULL. Name of the secondary negative fraction (e.g., 10%)
#' @slot presorting_fraction_name Character or NULL. Name of the pre-sorting (whole community) fraction
#' @slot ig_freq_name Character or NULL. Name of the `sample_data` column the Ig+ frequencies were
#'   read from. The resolved per-fraction values themselves live in the `sample_data` slot, as the
#'   `presort_ig_freq`/`pos_ig_freq`/`neg1_ig_freq`/`neg2_ig_freq` columns.
#' @slot ig_freq_layout Character or NULL. Which layout `ig_freq_name` was read with, `"wide"`
#'   (one Ig+ frequency per sample) or `"long"` (one per sort fraction) -- see [getPhyloIgSeq()].
#' @slot ellipse_coords A data.frame or NULL. Stores coordinates for sliding Z-score ellipses,
#'   with a \code{change_transform} column identifying which change axis each block belongs to
#' @slot sample_data A data.frame or NULL. Optional metadata for each sample
#' @slot tax_table A data.frame or NULL. Taxonomic information
#' @slot total_reads A data.frame or NULL. Total read counts per sample and fraction before rarefaction
#' @slot imputed_taxa List or NULL. Taxa that had zeros imputed, stored per sample
#'
#' @exportClass PhyloIgSeq
setClass(
  Class = "PhyloIgSeq",
  slots = list(
    ig_coating = "data.frame",
    score_names = "character",
    ma_coords = "data.frameOrNULL",
    positive_fraction_name = "character",
    first_negative_fraction_name = "characterOrNULL", # 9/10 of the whole negative fraction for IgSeq
    second_negative_fraction_name = "characterOrNULL", # 1/10 -//-
    presorting_fraction_name = "characterOrNULL", # before sorting
    ig_freq_name = "characterOrNULL",
    ig_freq_layout = "characterOrNULL",
    ellipse_coords = "data.frameOrNULL",
    sample_data = "data.frameOrNULL",
    tax_table = "data.frameOrNULL",
    total_reads = "data.frameOrNULL",
    imputed_taxa = "listOrNULL"
  ),
  prototype = list(score_names = character(0))
)

#' @rdname PhyloIgSeq-class
#' @param object A PhyloIgSeq object.
setMethod("show", "PhyloIgSeq", function(object) {
  n_samples <- length(unique(object@ig_coating$sample_id))
  n_taxa <- length(unique(object@ig_coating$taxon_id))
  scores_label <- if (length(object@score_names) > 0) {
    paste(object@score_names, collapse = ", ")
  } else {
    "(none computed)"
  }

  fraction_names <- c(
    positive = paste(object@positive_fraction_name, collapse = ", "),
    neg1 = object@first_negative_fraction_name,
    neg2 = object@second_negative_fraction_name,
    presort = object@presorting_fraction_name
  )

  cat("PhyloIgSeq-class Ig-coating scoring result\n")
  cat(sprintf(
    "ig_coating       Ig scores:      [ %s ] across %d sample(s), %d taxa\n",
    scores_label,
    n_samples,
    n_taxa
  ))
  if (length(fraction_names) > 0) {
    cat(
      "Fractions        ",
      paste(sprintf("%s=\"%s\"", names(fraction_names), fraction_names), collapse = ", "),
      "\n"
    )
  }
  if (!is.null(object@sample_data)) {
    cat(sprintf(
      "sample_data()    Sample Data:    [ %d samples by %d sample variables ]\n",
      nrow(object@sample_data),
      ncol(object@sample_data)
    ))
  }
  if (!is.null(object@tax_table)) {
    cat(sprintf(
      "tax_table()      Taxonomy Table: [ %d taxa by %d taxonomic ranks ]\n",
      nrow(object@tax_table),
      ncol(object@tax_table)
    ))
  }
  invisible(NULL)
})

#' Collapse a List of PhyloIgSeq objects
#'
#' Row-binds the \code{ig_coating}, \code{ellipse_coords}, \code{sample_data} and
#' \code{imputed_taxa} of a list of per-sample \code{\link{PhyloIgSeq-class}} objects (as produced
#' internally by \code{\link{getPhyloIgSeq}}) into a single object.
#'
#' @param phyloigseq_list A list of \code{PhyloIgSeq} objects. Elements that are not a
#'   \code{PhyloIgSeq} object are skipped.
#'
#' @return A single \code{PhyloIgSeq} object. Its \code{score_names} is the union of
#'   \code{score_names} across all input objects; \code{tax_table} and the fraction-name slots are
#'   not combined here (left at their defaults) and are set by the caller.
#'
#' @examples
#' pis_1 <- new(
#'   "PhyloIgSeq",
#'   ig_coating = data.frame(taxon_id = 1:2, sample_id = "s1", slide_z = c(0.5, -0.2)),
#'   score_names = "slide_z",
#'   positive_fraction_name = "pos",
#'   first_negative_fraction_name = "neg"
#' )
#' pis_2 <- new(
#'   "PhyloIgSeq",
#'   ig_coating = data.frame(taxon_id = 1:2, sample_id = "s2", slide_z = c(1.1, 0.3)),
#'   score_names = "slide_z",
#'   positive_fraction_name = "pos",
#'   first_negative_fraction_name = "neg"
#' )
#' collapsePhyloIgSeq(list(pis_1, pis_2))
#'
#' @export
collapsePhyloIgSeq <- function(phyloigseq_list) {
  ig_coating <- data.frame()
  sample_data <- data.frame()
  ellipse_coords <- data.frame()
  ma_coords <- data.frame()
  imputed_taxa <- list()
  score_names <- character(0)
  for (phyloigseq_obj in phyloigseq_list) {
    if (is(phyloigseq_obj, "PhyloIgSeq")) {
      # bind_rows() matches columns by name, fills in NA for missing columns
      ig_coating <- bind_rows(ig_coating, phyloigseq_obj@ig_coating)
      ellipse_coords <- bind_rows(ellipse_coords, phyloigseq_obj@ellipse_coords)
      ma_coords <- bind_rows(ma_coords, phyloigseq_obj@ma_coords)
      sample_data <- bind_rows(sample_data, phyloigseq_obj@sample_data)
      imputed_taxa <- c(imputed_taxa, phyloigseq_obj@imputed_taxa)
      score_names <- union(score_names, phyloigseq_obj@score_names)
    }
  }

  return(new(
    Class = "PhyloIgSeq",
    ig_coating = ig_coating,
    score_names = score_names,
    ma_coords = ma_coords,
    ellipse_coords = ellipse_coords,
    sample_data = sample_data,
    tax_table = NULL,
    imputed_taxa = imputed_taxa
  ))
}


#' Report How Much of a Run the Un-Mixing Had to Clamp
#'
#' Warns once per [getPhyloIgSeq()] call when a material share of taxa fall outside the
#' admissible cone on the purity-corrected change axis. Their change value is then set by
#' `pool_prior` rather than measured, so their sliding Z-score's magnitude -- and its
#' significance call -- are properties of the regularization, not of the data.
#'
#' @param ma_coords The `ma_coords` slot of the assembled object.
#' @param threshold Share of out-of-cone taxa above which to warn.
#'
#' @return `invisible(NULL)`, called for the `warning()`.
#' @noRd
.warn_out_of_cone_rate <- function(ma_coords, threshold = 0.1) {
  if (is.null(ma_coords) || prod(dim(ma_coords)) == 0) {
    return(invisible(NULL))
  }
  purity <- ma_coords[
    as.character(ma_coords$change_transform) == "purity_corrected",
  ]
  if (nrow(purity) == 0) {
    return(invisible(NULL))
  }
  for (which_pair in c("obs", "null")) {
    flag <- purity[[paste0(which_pair, "_in_cone")]]
    flag <- flag[!is.na(flag)]
    if (length(flag) == 0) next
    share <- mean(!flag)
    if (share <= threshold) next
    warning(
      sprintf("%.0f%%", 100 * share),
      " of the ",
      if (which_pair == "obs") "Ig+ vs Ig-" else "Ig-.1 vs Ig-.2",
      " pairs fall outside the admissible cone of the purity-corrected change ",
      "axis, so the un-mixing had to clamp them. Their change value is then set ",
      "by `pool_prior` rather than measured: their `purity_corrected_slide_z` ",
      "ranking is still meaningful, but its magnitude and any |Z| threshold ",
      "applied to it are not. Use the `obs_in_cone`/`null_in_cone` columns of ",
      "the `ma_coords` slot to separate or exclude them",
      if (which_pair == "null") {
        ", and note that a high rate here also distorts the null every Z-score is divided by"
      },
      ". A high rate means `pos_ig_freq` and `neg_ig_freq` are close together, ",
      "or that capture efficiency varies by taxon in a way the model does not ",
      "represent.\n"
    )
  }
  invisible(NULL)
}

#' MA-Plot Geometry Behind the Sliding Z-Scores
#'
#' Accessor for the \code{ma_coords} slot: the per-taxon/sample MA-plot coordinates each
#' sliding Z-score was derived from. Long in \code{change_transform}, one block of rows per
#' change axis (see \code{\link{get_ma_coordinates}}), so requesting both
#' \code{"slide_z"} and \code{"purity_corrected_slide_z"} from \code{\link{getPhyloIgSeq}}
#' yields two blocks and a consumer can switch between them by filtering instead of
#' recomputing.
#'
#' @param phyloigseq_obj A \code{\link{PhyloIgSeq-class}} object.
#' @param change_transform Optional. Restrict to one change axis (\code{"log_ratio"} or
#'   \code{"purity_corrected"}); \code{NULL} (the default) returns every axis present.
#'
#' @return A data frame with \code{sample_id}, \code{taxon_id}, \code{change_transform},
#'   \code{obs_abundance}, \code{obs_change}, \code{null_abundance}, \code{null_change},
#'   \code{obs_in_cone}, \code{null_in_cone} and \code{ellipse_level}; or \code{NULL} when no
#'   sliding Z-score was computed.
#'
#' @examples
#' data(ps_igseq)
#' pis <- getPhyloIgSeq(
#'   physeq = ps_igseq,
#'   sample_ids = c("sample_1", "sample_2"),
#'   sample_id_name = "sample_id",
#'   fraction_id_name = "sorting_fraction",
#'   positive_fraction_name = "pos",
#'   first_negative_fraction_name = "neg1",
#'   second_negative_fraction_name = "neg2",
#'   scores = "slide_z"
#' )
#' head(ma_coords(pis))
#'
#' @export
ma_coords <- function(phyloigseq_obj, change_transform = NULL) {
  coords <- phyloigseq_obj@ma_coords
  if (is.null(coords) || prod(dim(coords)) == 0) {
    return(coords)
  }
  if (!is.null(change_transform)) {
    change_transform <- match.arg(
      change_transform,
      choices = c("log_ratio", "purity_corrected")
    )
    coords <- coords[
      as.character(coords$change_transform) == change_transform,
    ]
  }
  coords
}


#' Names of the Per-Fraction Abundance Columns in \code{ig_coating}
#'
#' Returns which \code{ig_coating} columns hold per-fraction abundances, derived from the
#' object's own fraction-name slots. Use this rather than subtracting a hand-maintained set of
#' identifier/score/diagnostic names from \code{colnames(ig_coating)}: that denylist approach
#' silently absorbs any new diagnostic column as though it were a fraction, which is how a
#' column like \code{zeros_imputed} ends up offered as an abundance to weight by.
#'
#' @param phyloigseq_obj A \code{\link{PhyloIgSeq-class}} object.
#'
#' @return A character vector of column names present in \code{ig_coating}. In
#'   multiple-positive-fraction mode the positive fraction's column is
#'   \code{"positive_fraction_abundance"} rather than the fraction's own name (see
#'   \code{\link{getPhyloIgSeq}}), and that is what is returned.
#'
#' @examples
#' data(ps_igseq)
#' pis <- getPhyloIgSeq(
#'   physeq = ps_igseq,
#'   sample_ids = c("sample_1", "sample_2"),
#'   sample_id_name = "sample_id",
#'   fraction_id_name = "sorting_fraction",
#'   positive_fraction_name = "pos",
#'   first_negative_fraction_name = "neg1",
#'   second_negative_fraction_name = "neg2",
#'   scores = "palm"
#' )
#' ig_fraction_names(pis)
#'
#' @export
ig_fraction_names <- function(phyloigseq_obj) {
  candidates <- c(
    phyloigseq_obj@positive_fraction_name,
    "positive_fraction_abundance",
    phyloigseq_obj@first_negative_fraction_name,
    phyloigseq_obj@second_negative_fraction_name,
    phyloigseq_obj@presorting_fraction_name
  )
  intersect(unique(candidates), colnames(phyloigseq_obj@ig_coating))
}

#' Compute Ig Scores from a Phyloseq Object
#'
#' This function computes various Ig scores based on sample fraction data.
#'
#' @param physeq A `phyloseq` object containing raw count data.
#' @param taxrank Character or NULL. Taxonomic rank to agglomerate to before computing scores.
#' @param sample_id_name Name of the column identifying unique samples.
#' @param sample_ids Optional. A character vector of sample IDs to subset.
#' @param fraction_id_name Name of the column indicating fraction (e.g., pos, neg).
#' @param rarefy_by_sample Logical. Rarefy read counts across fractions within each sample.
#' @param transform_by_sample Transformation method (e.g., "identity", "log").
#' @param positive_fraction_name Name of the positive fraction. Usually a single value, but can
#'   be a character vector of several positive fraction values sharing the same negative/pre-sort
#'   fraction(s) below (e.g. multiple phage-display sort rounds/output pools sequenced against one
#'   shared negative/whole-community pool for the same samples). When more than one is given:
#'   \itemize{
#'     \item each (sample, positive fraction) combination becomes its own row, keyed by a synthetic
#'       `sample_id` of the form `"<original sample_id>_<positive fraction>"`;
#'     \item the negative/pre-sort fraction(s) are rarefied \emph{once}, jointly across every
#'       positive fraction of a sample, so they're compared against each positive fraction at the
#'       same common depth (rather than independently re-rarefied per positive fraction, which
#'       would add spurious noise between otherwise-comparable positive fractions);
#'     \item `sample_data` gains two columns, `original_sample_id` and `positive_fraction_name`,
#'       recording which original sample and positive fraction each row came from;
#'     \item the raw positive-fraction abundance column (normally named after the fraction's own
#'       value, e.g. `"pos"`) is instead named `positive_fraction_abundance` and is always
#'       populated, holding each row's own positive fraction's abundance.
#'   }
#'   With a single positive fraction (the default/common case), output is unchanged from previous
#'   versions: no suffix, no extra columns, original fraction-named abundance column.
#'
#'   Rows from the same original sample across different positive fractions are \strong{not}
#'   independent replicates -- they share the same underlying biology and the same negative/pre-sort
#'   fraction.
#' @param first_negative_fraction_name Name of the first negative fraction, or `NULL`
#'   to run without a negative fraction. In that case, `presorting_fraction_name` and
#'   `ig_freq_name` must both be supplied instead, and only `"prob_index"` (of `scores`)
#'   can be computed — every other score, and `slide_z`'s MA-plot geometry, is a
#'   positive-vs-negative comparison and requires a negative fraction.
#' @param second_negative_fraction_name Optional. Name of a second negative fraction.
#' @param presorting_fraction_name Optional. Name of the presorting fraction.
#' @param ig_freq_name Optional. Column name (in `physeq`'s `sample_data`) for the
#'   Ig+ frequency phenotype, if precomputed. How it is read depends on `ig_freq_layout`.
#' @param ig_freq_units Units `ig_freq_name`'s values are recorded in: `"frequency"`
#'   (default, already a probability in `[0, 1]`) or `"percent"` (`[0, 100]`, divided by
#'   100 before use). After conversion, a value outside `[0, 1]` for a given sample is
#'   treated as `NA` for that sample (with a `warning()`) rather than fed into
#'   `compute_ig_score()` — matching this function's general "skip what can't be
#'   computed, don't abort the batch" behavior (see [get_slide_z()]).
#' @param ig_freq_layout What one `ig_freq_name` value means:
#'   \describe{
#'     \item{`"wide"` (default)}{One Ig+ frequency per biological *sample*, i.e. the same value
#'       repeated on each of its fraction rows — the historical layout. Only the overall
#'       P(Ig+) is recoverable, which is all `"prob_index"`/`"prob_ratio"` need. A column that
#'       is *not* constant across a sample's fractions can't be used and produces a `warning()`
#'       pointing at `"long"`.}
#'     \item{`"long"`}{The frequency is *fraction-specific*: each fraction's row records the Ig+
#'       frequency measured in that fraction. This yields the pre-sort P(Ig+) plus the
#'       "positive fraction purity" and "negative fraction impurity" that the
#'       `"purity_corrected_*"` scores require.}
#'   }
#'   Either way the resolved values are stored per sample in the returned object's `sample_data`
#'   slot, as `presort_ig_freq`, `pos_ig_freq`, `neg1_ig_freq` and `neg2_ig_freq` (the last three
#'   are `NA` under `"wide"`). With several `positive_fraction_name`s, `pos_ig_freq` follows each
#'   row's own positive fraction.
#' @param zero_treatment How to handle zeros ("no_zero", "pseudocount", etc.).
#' @param window_size Integer. Window size for smoothing.
#' @param empirical_null_distribution Logical. Whether to scale the sliding Z-score by
#'   the empirical (Ig-.1 vs Ig-.2) null's sd rather than the observed distribution's.
#'   Chooses the scale only; `center_on` chooses the location. See [get_slide_z()].
#' @param center_on Which distribution locates the sliding Z-score: `"reference"`
#'   (default, the mean of whatever supplies the scale, which at matched depth keeps `0`
#'   meaning "coated at the community's own Ig+ frequency") or `"observed"` (the local
#'   median observed change, which re-anchors on the typical taxon). Passed to
#'   [get_slide_z()]; read [compute_slide_z()]'s "Where the centre comes from" before
#'   changing it.
#' @param confidence_levels Optional. Confidence levels for scoring.
#' @param scores Vector of score names to compute, from [IG_SCORES] (the default is all of
#'   them). A requested score whose inputs are not available -- a fraction it divides by, or
#'   the Ig+ frequency it weights with -- is dropped with a `warning()` rather than returned
#'   as an all-NA column, so `score_names` on the result lists what was actually computed.
#' @param taxon_id_source How to derive the `taxon_id` used throughout `ig_coating`/`tax_table`:
#'   `"sequential"` (default) renumbers taxa as fresh sequential integers, recoverable via
#'   `tax_table$taxon_name`; `"original"` uses `physeq`'s own (possibly `taxrank`-agglomerated,
#'   and not necessarily ASV-level if `physeq` was already agglomerated upstream) taxa name
#'   directly as `taxon_id`, matching the identifiers shown by [group_sorted_samples()] when
#'   called directly (e.g. for a single-sample preview).
#'
#' @return A \code{\link{PhyloIgSeq-class}} object -- always a single object, whether one or
#'   several \code{positive_fraction_name} values were given (see that parameter's documentation
#'   for the schema differences in the latter case). Its \code{ig_coating} slot holds one row per
#'   taxon/sample with the requested \code{scores} as columns (also recorded in \code{score_names})
#'   plus supporting fraction/diagnostic columns; see \code{\link{get_ig_score}} to retrieve a
#'   single score without dealing with the rest of \code{ig_coating}.
#'
#' @examples
#' data(ps_igseq)
#' pis <- getPhyloIgSeq(
#'   physeq = ps_igseq,
#'   sample_ids = c("sample_1", "sample_2", "sample_3"),
#'   sample_id_name = "sample_id",
#'   fraction_id_name = "sorting_fraction",
#'   positive_fraction_name = "pos",
#'   first_negative_fraction_name = "neg1",
#'   second_negative_fraction_name = "neg2",
#'   scores = c("slide_z", "palm", "kau")
#' )
#' pis
#' get_ig_score(pis, score_name = "slide_z", sample_ids = "sample_1")
#'
#' @export
getPhyloIgSeq <- function(
  physeq, # containing raw counts
  taxrank = NULL,
  sample_id_name, # name identifying each sample_id (individual for
  # which sorting was performed)
  sample_ids = NULL, # if NULL, all samples taken
  fraction_id_name, # name of the fraction column
  rarefy_by_sample = TRUE, # inside each sample, rarefy abundances for each
  # fraction (so that fractions of the same sample have the same total sum of reads)
  transform_by_sample = "identity",
  positive_fraction_name = "pos",
  first_negative_fraction_name = NULL,
  second_negative_fraction_name = NULL,
  presorting_fraction_name = NULL, # before sorting
  ig_freq_name = NULL,
  ig_freq_units = c("frequency", "percent"),
  ig_freq_layout = c("wide", "long"),
  zero_treatment = "no_zero",
  window_size = 50,
  empirical_null_distribution = TRUE,
  center_on = c("reference", "observed"),
  confidence_levels = NULL,
  scores = IG_SCORES,
  taxon_id_source = c("sequential", "original")
) {
  taxon_id_source <- match.arg(taxon_id_source)
  ig_freq_units <- match.arg(ig_freq_units)
  ig_freq_layout <- match.arg(ig_freq_layout)
  center_on <- match.arg(center_on)
  if (
    is.null(first_negative_fraction_name) &&
      (is.null(presorting_fraction_name) || is.null(ig_freq_name))
  ) {
    stop(
      "Either `first_negative_fraction_name` must be supplied, or both ",
      "`presorting_fraction_name` and `ig_freq_name` must be supplied ",
      "(required to compute `prob_index` without a negative fraction)."
    )
  }
  if (
    is.null(first_negative_fraction_name) &&
      any(names(SLIDE_Z_SCORES) %in% scores)
  ) {
    warning(
      "No negative fraction furnished, cannot compute ",
      paste0(
        "`",
        intersect(names(SLIDE_Z_SCORES), scores),
        "`",
        collapse = " / "
      ),
      " (or the MA-plot geometry they depend on); dropping from `scores`.\n"
    )
    scores <- setdiff(scores, names(SLIDE_Z_SCORES))
  }
  # The purity-corrected scores need the Ig+ frequency measured *inside* each
  # sorted fraction, i.e. a negative fraction and a "long" ig_freq column. Drop
  # them rather than emit all-NA columns -- they are part of the default `scores`
  # since they joined IG_SCORES, so most calls would otherwise pick them up
  # silently.
  #
  # Their requirements are NOT the same, so they are gated per score rather than
  # as one block: `purity_corrected_prob_index`/`_ratio` use the pre-sort Ig+
  # frequency as a Bayesian prior and so need a pre-sort fraction to read it off
  # under the long layout, whereas `purity_corrected_slide_z` does not -- that
  # prior is an additive constant on its change axis, which the Z-score's
  # centering removes. Don't "simplify" this back into a single check.
  purity_requirements <- list(
    purity_corrected_prob_index = TRUE,
    purity_corrected_prob_ratio = TRUE,
    purity_corrected_slide_z = FALSE # needs no pre-sort fraction
  )
  for (score in intersect(names(purity_requirements), scores)) {
    missing_inputs <- c(
      if (is.null(first_negative_fraction_name)) "a negative fraction",
      if (purity_requirements[[score]] && is.null(presorting_fraction_name)) {
        "a presorting fraction"
      },
      if (is.null(ig_freq_name)) {
        "`ig_freq_name`"
      } else if (!identical(ig_freq_layout, "long")) {
        "`ig_freq_layout = \"long\"` (the per-fraction Ig+ frequencies)"
      }
    )
    if (length(missing_inputs) > 0) {
      warning(
        "Cannot compute `",
        score,
        "` without ",
        paste(missing_inputs, collapse = " and "),
        "; dropping from `scores`.\n"
      )
      scores <- setdiff(scores, score)
    }
  }
  # `prob_index`/`prob_ratio` are Bayesian: both weight the sorted fractions by
  # the pre-sort P(Ig+), which comes from `ig_freq_name`. Gated for the same
  # reason as the purity-corrected block above -- they are in the default
  # `scores`, so a call that never asked for them would otherwise carry a silent
  # all-NA column -- and, like that block, per score rather than as one, since
  # their inputs differ:
  #
  #   * `prob_index` is P(taxon | Ig+) P(Ig+) / P(taxon), so the denominator is
  #     the pre-sort fraction's own abundance: it needs that fraction under
  #     either layout.
  #   * `prob_ratio` divides by the negative fraction instead, so it needs a
  #     negative fraction, and a pre-sort fraction only under
  #     `ig_freq_layout = "long"`, where P(Ig+) is recorded per fraction and is
  #     therefore readable only off the pre-sort rows. Under `"wide"` P(Ig+) is a
  #     property of the sample and any fraction's row carries it.
  prior_requirements <- list(
    prob_index = "presort",
    prob_ratio = if (identical(ig_freq_layout, "long")) "presort" else "negative"
  )
  for (score in intersect(names(prior_requirements), scores)) {
    missing_inputs <- c(
      if (
        identical(prior_requirements[[score]], "presort") &&
          is.null(presorting_fraction_name)
      ) {
        "a presorting fraction"
      },
      if (
        identical(prior_requirements[[score]], "negative") &&
          is.null(first_negative_fraction_name)
      ) {
        "a negative fraction"
      },
      if (is.null(ig_freq_name)) "`ig_freq_name`"
    )
    if (length(missing_inputs) > 0) {
      warning(
        "Cannot compute `",
        score,
        "` without ",
        paste(missing_inputs, collapse = " and "),
        if (
          identical(score, "prob_ratio") && identical(ig_freq_layout, "long")
        ) {
          paste0(
            " (under `ig_freq_layout = \"long\"` the pre-sort Ig+ frequency it ",
            "uses as a prior is recorded on the presorting fraction's own rows)"
          )
        },
        "; dropping from `scores`.\n"
      )
      scores <- setdiff(scores, score)
    }
  }
  # Only relevant when a sliding Z-score is actually being computed -- otherwise
  # empirical_null_distribution is never read.
  if (
    any(names(SLIDE_Z_SCORES) %in% scores) &&
      is.null(second_negative_fraction_name) &&
      empirical_null_distribution
  ) {
    warning(
      "No second negative fraction furnished, cannot model empirical null ( Ig-.1 vs Ig-.2) distribution...\n"
    )
    empirical_null_distribution <- FALSE
  }
  if (!is.null(taxrank)) {
    physeq <- tax_glom(physeq = physeq, taxrank = taxrank)
    taxa_names(physeq) <- make.unique(tax_table(physeq)[, taxrank])
  }

  original_taxa_names <- taxa_names(physeq)
  if (taxon_id_source == "sequential") {
    # To make matching and computation easier, don't carry long strings as taxa names
    taxon_ids <- seq_along(original_taxa_names)
    names(original_taxa_names) <- taxon_ids # keep the mapping
    taxa_names(physeq) <- taxon_ids
  } else {
    # "original": keep physeq's own (possibly taxrank-agglomerated) taxa_names
    # as taxon_id directly instead of renumbering — original_taxa_names becomes
    # an identity mapping so downstream lookups (e.g. `taxon_name = original_taxa_names[taxon_ids_to_keep]`
    # below) still work unchanged.
    names(original_taxa_names) <- original_taxa_names
  }

  all_fraction_names <- c(
    positive_fraction_name,
    first_negative_fraction_name,
    second_negative_fraction_name,
    presorting_fraction_name
  )

  grouped_data <-
    group_sorted_samples(
      physeq = physeq,
      taxrank = NULL, # already done above
      sample_id_name = sample_id_name,
      sample_ids = sample_ids,
      fraction_id_name = fraction_id_name,
      fraction_ids = all_fraction_names,
      rarefy_by_sample = rarefy_by_sample,
      # NOTE: only pos, neg1 and neg2 are rarefied
      fractions_to_rarefy = c(
        positive_fraction_name,
        first_negative_fraction_name,
        second_negative_fraction_name
      ),
      transform_by_sample = transform_by_sample
    )

  if (is.null(sample_ids)) {
    sample_ids <- names(grouped_data)
  } else {
    sample_ids <- sample_ids[sample_ids %in% names(grouped_data)]
  }

  if (length(sample_ids) == 0) {
    warning("No suitable samples left, returning NULL\n")
    return(NULL)
  }

  metadata <- phyloseq::sample_data(physeq) %>% as("data.frame")
  # to avoid potential problems if there is already a column called "sample_id"
  # which is not a name of smaple_id column indicated by a user
  names(metadata)[
    names(metadata) != sample_id_name & names(metadata) == "sample_id"
  ] <-
    paste0(
      "original___",
      names(metadata)[
        names(metadata) != sample_id_name & names(metadata) == "sample_id"
      ]
    )

  names(metadata)[names(metadata) == sample_id_name] <- "sample_id"
  # Same collision guard for the per-fraction Ig+ frequency columns this
  # function adds to `sample_data` below -- a user column of the same name would
  # otherwise be silently overwritten.
  IG_FREQ_COLUMNS <- c(
    "presort_ig_freq",
    "pos_ig_freq",
    "neg1_ig_freq",
    "neg2_ig_freq"
  )
  clashing <- names(metadata) %in% IG_FREQ_COLUMNS
  # `ig_freq_name` itself may be one of the clashing columns, in which case the
  # lookup below has to follow the rename -- while the slot keeps recording the
  # name the user actually passed.
  ig_freq_column <- ig_freq_name
  if (!is.null(ig_freq_column) && ig_freq_column %in% names(metadata)[clashing]) {
    ig_freq_column <- paste0("original___", ig_freq_column)
  }
  names(metadata)[clashing] <- paste0("original___", names(metadata)[clashing])
  # put sample_id on the first place
  metadata <- metadata[, c("sample_id", setdiff(names(metadata), "sample_id"))]

  phyloigseq_list <- list()

  total_reads <- data.frame()

  # When more than one positive fraction is given, each (sample x positive
  # fraction) combination becomes its own synthetic "sample" below, so that
  # every existing sample_id-keyed consumer (agglomeration, wide pivots,
  # phyloseq conversion, plots) gets correct per-fraction resolution for
  # free. A single positive fraction keeps today's sample_id untouched.
  multi_fraction <- length(positive_fraction_name) > 1

  for (sample_id in sample_ids) {
    # Fraction-agnostic, computed once per sample and reused for every
    # positive fraction below.

    # Keep track of total presorting read counts for each sample (used to compute relative abundances later)
    sample_total_reads_value <- if (!is.null(presorting_fraction_name)) {
      if (all(is.na(grouped_data[[sample_id]][[presorting_fraction_name]]))) {
        NA
      } else {
        sum(
          grouped_data[[sample_id]][[presorting_fraction_name]],
          na.rm = TRUE
        )
      }
    }

    # Retrieve sample  metadata
    sam_metadata_df <-
      metadata[
        metadata[["sample_id"]] == sample_id &
          metadata[[fraction_id_name]] %in% all_fraction_names, ,
        drop = FALSE
      ]
    sam_metadata_row <- data.frame(matrix(NA, nrow = 1, ncol = ncol(metadata)))
    names(sam_metadata_row) <- names(metadata)

    for (var_name in names(metadata)) {
      unique_values <- unique(sam_metadata_df[[var_name]])

      if (length(unique_values) == 1) {
        sam_metadata_row[[var_name]] <- unique_values
      }
    }
    sam_metadata_row <- sam_metadata_row[
      ,
      names(sam_metadata_row) != fraction_id_name
    ]

    # Fraction-specific: looped once per positive fraction, reusing the
    # sample-level grouping/rarefaction already computed above in
    # `grouped_data[[sample_id]]` (shared across all positive fractions of
    # this sample).
    for (pos in positive_fraction_name) {
      synthetic_id <- if (multi_fraction) {
        paste0(sample_id, "_", pos)
      } else {
        sample_id
      }

      if (!is.null(presorting_fraction_name)) {
        total_reads <- rbind(
          total_reads,
          data.frame(sample_id = synthetic_id, total_reads = sample_total_reads_value)
        )
      }

      # Resolve the Ig+ frequency of each fraction to a [0, 1] probability.
      # Inside the positive-fraction loop because `pos_ig_freq` is a property of
      # the positive fraction being scored, so it has to follow `pos` -- the
      # other three are the same for every iteration.
      ig_freqs <- resolve_ig_freqs(
        sam_metadata_df = sam_metadata_df,
        fraction_id_name = fraction_id_name,
        ig_freq_name = ig_freq_column,
        ig_freq_units = ig_freq_units,
        ig_freq_layout = ig_freq_layout,
        positive_fraction_name = pos,
        first_negative_fraction_name = first_negative_fraction_name,
        second_negative_fraction_name = second_negative_fraction_name,
        presorting_fraction_name = presorting_fraction_name,
        sample_id = sample_id
      )

      # Handle zeros
      this_fraction_names <- c(
        pos,
        first_negative_fraction_name,
        second_negative_fraction_name,
        presorting_fraction_name
      )
      present_fraction_names <- this_fraction_names[
        this_fraction_names %in% colnames(grouped_data[[sample_id]])
      ]

      zero_imputation_result <-
        impute_zeros(
          data = grouped_data[[sample_id]],
          # Don't impute zeros in other fractions!
          fraction_names = intersect(
            present_fraction_names,
            c(pos, first_negative_fraction_name, second_negative_fraction_name)
          ),
          method = zero_treatment
        )

      ig_coating <- zero_imputation_result$data %>%
        select(all_of(unique(c("taxon_id", "sample_id", present_fraction_names))))
      ig_coating$sample_id <- synthetic_id

      ig_coating$zeros_imputed <- ig_coating$taxon_id %in%
        zero_imputation_result$imputed_taxa

      if (nrow(ig_coating) == 0) {
        warning(paste0(
          synthetic_id,
          " excluded: no taxa left after zero treatment\n"
        ))
        next
      }

      # Compute scores.
      #
      # Every requested sliding Z-score is computed here, one per change axis (see
      # SLIDE_Z_SCORES). Only the score itself goes into `ig_coating`; the MA
      # geometry it came from goes to the `ma_coords` slot, long in
      # `change_transform`, since two axes cannot share one set of
      # obs_change/null_change columns.
      ellipse_coords <- data.frame()
      ma_coords_this <- data.frame()
      for (slide_score in intersect(names(SLIDE_Z_SCORES), scores)) {
        this_transform <- SLIDE_Z_SCORES[[slide_score]]
        slide_z_result <-
          get_slide_z(
            sorted_sample_df = ig_coating,
            positive_fraction_name = pos,
            first_negative_fraction_name = first_negative_fraction_name,
            second_negative_fraction_name = second_negative_fraction_name,
            window_size = window_size,
            empirical_null_distribution = empirical_null_distribution,
            center_on = center_on,
            confidence_levels = confidence_levels,
            imputed_taxa = zero_imputation_result$imputed_taxa,
            change_transform = this_transform,
            # NA outside `ig_freq_layout = "long"`; only the purity-corrected axis
            # reads them, and its score is dropped from `scores` above when the
            # long layout isn't in use.
            pos_ig_freq = ig_freqs$pos_ig_freq,
            neg_ig_freq = ig_freqs$neg1_ig_freq,
            # `ig_coating` is `zero_imputation_result$data` with columns selected,
            # so the two are still row-aligned. Without this the imputed counts
            # look like measured ones and every pair that was zero on both sides
            # would be scored off a manufactured change of exactly 0.
            was_zero = zero_imputation_result$was_zero
          )

        ig_coating[[slide_score]] <- slide_z_result$slide_z

        if (prod(dim(slide_z_result$ma_coords)) != 0) {
          ma_coords_this <- rbind(
            ma_coords_this,
            data.frame(
              sample_id = synthetic_id,
              taxon_id = slide_z_result$ma_coords$taxon_id,
              change_transform = this_transform,
              slide_z_result$ma_coords[, c(
                "obs_abundance",
                "obs_change",
                "null_abundance",
                "null_change",
                "obs_in_cone",
                "null_in_cone",
                "obs_estimable",
                "null_estimable"
              )],
              ellipse_level = if (is.null(slide_z_result$ellipse_level)) {
                NA
              } else {
                as.character(slide_z_result$ellipse_level)
              },
              row.names = NULL
            )
          )
        }
        if (prod(dim(slide_z_result$ellipse_coords)) != 0) {
          ellipse_coords <- rbind(
            ellipse_coords,
            data.frame(
              slide_z_result$ellipse_coords,
              change_transform = this_transform,
              row.names = NULL
            )
          )
        }
      }

      # Other Ig scores:
      # Both purity-corrected scores derive the same sort recovery from the same
      # three Ig+ frequencies, so an inconsistent set warns identically once per
      # score. That reads like the sample was processed twice; report each
      # distinct complaint once per (sample x positive fraction) instead.
      warned_here <- character(0)
      for (score in setdiff(scores, names(SLIDE_Z_SCORES))) {
        ig_coating[[score]] <-
          withCallingHandlers(
            compute_ig_score(
              method = score,
              pos = ig_coating[[pos]],
              neg = if (!is.null(first_negative_fraction_name)) {
                ig_coating[[first_negative_fraction_name]]
              },
              pre = if (!is.null(presorting_fraction_name)) {
                ig_coating[[presorting_fraction_name]]
              },
              presort_ig_freq = ig_freqs$presort_ig_freq,
              # NA outside `ig_freq_layout = "long"`; only the purity-corrected
              # scores read them, and those are dropped from `scores` above when
              # the long layout isn't in use.
              pos_ig_freq = ig_freqs$pos_ig_freq,
              neg_ig_freq = ig_freqs$neg1_ig_freq
            ),
            warning = function(w) {
              message_text <- conditionMessage(w)
              if (message_text %in% warned_here) {
                invokeRestart("muffleWarning")
              }
              warned_here <<- c(warned_here, message_text)
            }
          )
      }

      # Ig scores are the columns users look for first; keep them right after the
      # taxon_id/sample_id identifiers, ahead of fraction/diagnostic columns.
      score_names_present <- intersect(scores, names(ig_coating))
      ig_coating <- ig_coating %>%
        relocate(all_of(score_names_present), .after = "sample_id")

      sam_metadata_row_this <- sam_metadata_row
      sam_metadata_row_this$sample_id <- synthetic_id

      # Record the resolved (unit-converted, range-checked) Ig+ frequencies per
      # sample, so they can be inspected, exported, and used to colour/facet
      # downstream plots. These are the *measured* phenotypes only -- nothing
      # derived from them is stored.
      sam_metadata_row_this[names(ig_freqs)] <- ig_freqs

      if (multi_fraction) {
        sam_metadata_row_this$original_sample_id <- sample_id
        sam_metadata_row_this$positive_fraction_name <- pos

        # The raw positive-fraction abundance column is named after `pos`
        # itself, so after combining all positive fractions' rows it would
        # only be populated for the rows that came from that one fraction
        # (NA elsewhere) -- silently breaking anything that filters/weights
        # by a single named abundance column (e.g. agglomPhyloIgSeq()'s
        # `abundance_fraction`). Coalesce it into one universally-populated
        # column instead: every row's own positive fraction's abundance,
        # regardless of which fraction that was.
        if (pos %in% names(ig_coating)) {
          names(ig_coating)[names(ig_coating) == pos] <- "positive_fraction_abundance"
        }
      }

      imputed_taxa <- list()
      imputed_taxa[[synthetic_id]] <- zero_imputation_result$imputed_taxa

      phyloigseq_list[[synthetic_id]] <-
        new(
          Class = "PhyloIgSeq",
          ig_coating = ig_coating,
          score_names = score_names_present,
          ma_coords = ma_coords_this,
          positive_fraction_name = pos,
          first_negative_fraction_name = first_negative_fraction_name,
          second_negative_fraction_name = second_negative_fraction_name,
          ellipse_coords = ellipse_coords,
          sample_data = sam_metadata_row_this,
          tax_table = NULL,
          imputed_taxa = imputed_taxa
        )
    }
  }

  phyloigseq_obj <- collapsePhyloIgSeq(phyloigseq_list)

  if (prod(dim(total_reads)) != 0) {
    names(total_reads)[2] <- presorting_fraction_name
    phyloigseq_obj@total_reads <- total_reads[
      total_reads$sample_id %in% phyloigseq_obj@ig_coating$sample_id,
    ]
  } else {
    phyloigseq_obj@total_reads <- NULL
  }

  phyloigseq_obj@positive_fraction_name <- positive_fraction_name
  phyloigseq_obj@first_negative_fraction_name <- first_negative_fraction_name
  phyloigseq_obj@second_negative_fraction_name <- second_negative_fraction_name
  phyloigseq_obj@presorting_fraction_name <- presorting_fraction_name
  phyloigseq_obj@ig_freq_name <- ig_freq_name
  phyloigseq_obj@ig_freq_layout <- if (!is.null(ig_freq_name)) ig_freq_layout

  # Cohort-level saturation report for the purity-corrected axis. Warned once
  # here rather than per sample, because the number that matters is the share
  # across the run, and 50 identical per-sample warnings would bury it.
  #
  # This is not a cosmetic complaint. A saturated taxon's change value is set by
  # `pool_prior` rather than measured, while the null it is divided by is
  # technical-replicate noise, so its |Z| is large by construction: on the SMILE
  # cohort every single out-of-cone taxon clears |Z| > 1.96, and multiplying
  # `pool_prior` by 1000 moves the largest |Z| around without shrinking it. The
  # ranking among such taxa is meaningful; the magnitude and the significance
  # call are not.
  .warn_out_of_cone_rate(phyloigseq_obj@ma_coords)

  taxon_ids_to_keep <- unique(phyloigseq_obj@ig_coating$taxon_id)
  tax_table <- as.matrix(phyloseq::tax_table(physeq)@.Data)[
    taxon_ids_to_keep,
  ] %>%
    as.data.frame()

  # In the hierarchical order
  tax_table <- cbind(
    tax_table,
    data.frame(
      taxon_name = original_taxa_names[taxon_ids_to_keep],
      taxon_id = taxon_ids_to_keep
    )
  )
  rownames(tax_table) <- NULL
  # TODO: keep only taxa that are left
  phyloigseq_obj@tax_table <- tax_table

  return(phyloigseq_obj)
}

#' Get an Ig Score from a PhyloIgSeq Object
#'
#' Retrieves a single Ig score from \code{ig_coating}, optionally restricted to a subset of taxa
#' and/or samples, without having to know which of \code{ig_coating}'s other columns are metadata.
#'
#' @param phyloigseq_obj A \code{\link{PhyloIgSeq-class}} object.
#' @param score_name Character. Name of the score to retrieve; must be one of
#'   \code{phyloigseq_obj@score_names}.
#' @param taxa_ids Optional. A vector of taxon IDs to restrict to; \code{NULL} (the default) keeps
#'   all taxa.
#' @param sample_ids Optional. A vector of sample IDs to restrict to; \code{NULL} (the default)
#'   keeps all samples.
#'
#' @return A data frame with columns \code{taxon_id}, \code{sample_id} and \code{score_name},
#'   filtered to \code{taxa_ids}/\code{sample_ids} if given.
#'
#' @examples
#' data(ps_igseq)
#' pis <- getPhyloIgSeq(
#'   physeq = ps_igseq,
#'   sample_ids = c("sample_1", "sample_2", "sample_3"),
#'   sample_id_name = "sample_id",
#'   fraction_id_name = "sorting_fraction",
#'   positive_fraction_name = "pos",
#'   first_negative_fraction_name = "neg1",
#'   second_negative_fraction_name = "neg2",
#'   scores = c("slide_z", "palm", "kau")
#' )
#' get_ig_score(pis, score_name = "palm", sample_ids = c("sample_1", "sample_2"))
#'
#' @export
get_ig_score <- function(
  phyloigseq_obj,
  score_name,
  taxa_ids = NULL,
  sample_ids = NULL
) {
  if (!is(phyloigseq_obj, "PhyloIgSeq")) {
    stop("`phyloigseq_obj` must be a PhyloIgSeq object")
  }
  if (!score_name %in% phyloigseq_obj@score_names) {
    stop(
      "`score_name` must be one of: ",
      paste(phyloigseq_obj@score_names, collapse = ", ")
    )
  }

  result <- phyloigseq_obj@ig_coating[, c("taxon_id", "sample_id", score_name)]
  if (!is.null(taxa_ids)) {
    result <- result[result$taxon_id %in% taxa_ids, ]
  }
  if (!is.null(sample_ids)) {
    result <- result[result$sample_id %in% sample_ids, ]
  }
  rownames(result) <- NULL
  result
}

# Mimics seq_table() function with these differences:
# - naming of variables and columns
# - preprocessing (rarefaction + transformation) integrated in the function
# - only those taxa that have zero counts for all fractions are excluded
# - no restriction of names or number of fractions (but they have to be unique for each sample)
