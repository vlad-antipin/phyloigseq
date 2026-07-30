# Generates the ps_phage_display toy dataset: an entirely synthetic
# phyloseq object simulating a phage-display panning experiment, intended
# as a toy dataset for trying out getPhyloIgSeq()'s support for multiple
# positive fractions (see ?ps_phage_display and the "Multiple positive
# fractions" section of vignette("igseq-scoring", package = "PhyloIgSeq")).
#
# Story: 8 independent panning campaigns, each sequenced at 4 stages -- the
# naive/unselected "input" library (shared negative/pre-sort fraction) and
# three sequential output rounds ("round1"/"round2"/"round3"), each round
# progressively enriched for a handful of campaign-specific "responder"
# clones (monotonic, compounding enrichment across rounds, as in real
# panning). Taxa represent individual phage clone sequences rather than
# taxonomic units, so tax_table carries no real taxonomy -- just a single
# sequence_id column identical to the taxon name.
#
# Re-run this script (after adjusting as needed) to regenerate
# data/ps_phage_display.rda; devtools::document() picks up the roxygen
# block in R/data.R.

set.seed(20260730)

n_clones <- 400
n_samples <- 8
fractions <- c("input", "round1", "round2", "round3")

clone_ids <- sprintf("clone_%03d", seq_len(n_clones))
sample_ids <- sprintf("campaign_%d", seq_len(n_samples))

target_antigens <- c("antigen_A", "antigen_B", "antigen_C", "antigen_D")
library_types <- c("naive", "synthetic", "immune")
host_species <- c("human", "camelid", "mouse")
operators <- c("operator_1", "operator_2", "operator_3")
seq_batches <- c("batch_1", "batch_2")

# Campaign-level metadata, constant across all fraction rows of a sample_id
# (like sex/age in ps_igseq). Balanced round-robin + shuffle (rather than iid
# sample(..., replace = TRUE)) so every level gets a reasonable number of
# campaigns even with only 8 of them -- otherwise a small level could end up
# with 0-1 campaigns by chance alone, of little use for group comparisons.
balanced_levels <- function(levels, n) sample(rep_len(levels, n))

sample_meta <- data.frame(
  sample_id = sample_ids,
  target_antigen = balanced_levels(target_antigens, n_samples),
  library_type = balanced_levels(library_types, n_samples),
  host_species = balanced_levels(host_species, n_samples),
  # Synthetic independent-assay phenotype (e.g. ELISA-estimated fraction of
  # target-specific phage) -- usable as ig_freq_name together with
  # presorting_fraction_name = "input" to exercise the positive-fraction-only
  # workflow, the same role ps_igseq's ig_pheno column plays.
  pct_target_specific = round(runif(n_samples, 0.01, 0.35), 3),
  stringsAsFactors = FALSE
)

# Baseline (input-library) log-abundance shared across campaigns, like a
# typical naive-library composition: most clones rare, a handful more
# common. Per-campaign jitter is added below so no two campaigns' input
# fractions are identical.
baseline_log_mean <- rnorm(n_clones, mean = 0, sd = 1.2)

row_meta <- list()
counts <- list()

for (i in seq_len(n_samples)) {
  sid <- sample_ids[i]

  # A handful of clones specific to this campaign's target, enriched over
  # rounds; drawn per campaign so different campaigns enrich different
  # clones (target_antigen-dependent signal a user can explore).
  n_responders <- sample(3:8, 1)
  responders <- sample(clone_ids, n_responders)

  jitter <- rnorm(n_clones, 0, 0.3)
  rel <- exp(baseline_log_mean + jitter)
  rel <- rel / sum(rel)

  # Per-round multiplicative enrichment factor (and wash-count step) for
  # this campaign, drawn once and reused each round so enrichment -- and
  # wash stringency -- compound/increase monotonically (input -> round1 ->
  # round2 -> round3), as in a real panning protocol, rather than being
  # redrawn (and so potentially non-monotonic) each round.
  enrich_per_round <- 3 + rgamma(1, shape = 4, rate = 1) # roughly 3-10x/round
  wash_step <- sample(3:5, 1)

  for (fr in fractions) {
    round_index <- match(fr, fractions) - 1L # 0 for input, 1/2/3 for rounds

    if (round_index > 0) {
      rel[clone_ids %in% responders] <- rel[clone_ids %in% responders] * enrich_per_round
      rel <- rel / sum(rel)
    }

    depth <- sample(8000:30000, 1)
    cts <- as.integer(stats::rmultinom(1, size = depth, prob = rel))

    row_id <- paste0(sid, "_", fr)
    counts[[row_id]] <- cts

    # Output titer (log10 cfu/mL) typically climbs with successive rounds of
    # enrichment; wash stringency (number of washes) is usually increased
    # round over round too -- both vary per fraction row within a campaign,
    # like batch/operator in ps_igseq.
    row_meta[[row_id]] <- data.frame(
      sample_id = sid,
      phage_fraction = fr,
      operator = sample(operators, 1),
      sequencing_batch = sample(seq_batches, 1),
      output_titer_log10 = round(3.5 + round_index * runif(1, 0.6, 1.1) + rnorm(1, 0, 0.15), 2),
      wash_count = 3L + round_index * wash_step,
      stringsAsFactors = FALSE
    )
  }
}

sample_data_df <- do.call(rbind, row_meta)
sample_data_df <- merge(sample_data_df, sample_meta, by = "sample_id", sort = FALSE)
rownames(sample_data_df) <- paste0(sample_data_df$sample_id, "_", sample_data_df$phage_fraction)
sample_data_df <- sample_data_df[names(counts), ] # merge() can reorder rows

# samples x taxa (taxa_are_rows = FALSE), matching ps_igseq/ps_16s_refinement's
# orientation
otu_mat <- t(do.call(cbind, counts))
rownames(otu_mat) <- names(counts)
colnames(otu_mat) <- clone_ids

tax_mat <- matrix(clone_ids, ncol = 1, dimnames = list(clone_ids, "sequence_id"))

ps_phage_display <- phyloseq::phyloseq(
  phyloseq::otu_table(otu_mat, taxa_are_rows = FALSE),
  phyloseq::sample_data(sample_data_df),
  phyloseq::tax_table(tax_mat)
)

save(
  ps_phage_display,
  file = file.path("data", "ps_phage_display.rda"),
  compress = "xz"
)
