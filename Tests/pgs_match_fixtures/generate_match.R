#!/usr/bin/env Rscript
# Generate reference fixtures for PgsMatchOne from the installed bigsnpr::snp_match.
#
# Run once: Rscript Tests/pgs_match_fixtures/generate_match.R
# Output:   Tests/pgs_match_fixtures/match_reference.txt
#
# bigsnpr is only needed to (re-)generate fixtures; the harness reads the
# committed text file and needs no R, matching the arrangement used by
# score_regression_golden/ and the statsmodels logistic fixtures.

suppressMessages(library(bigsnpr))

cat("bigsnpr version:", as.character(packageVersion("bigsnpr")), "\n")

# ---------------------------------------------------------------------------
# Helper: classify each source variant from the snp_match output.
#
# snp_match returns only matched rows.  We must derive the status of the
# unmatched rows ourselves by comparing the source table to the output.
# ---------------------------------------------------------------------------
classify_variants <- function(sumstats_src, snp_match_out,
                              target_df, strand_flip_enabled = TRUE) {
  n_src <- nrow(sumstats_src)
  status     <- character(n_src)
  out_a0     <- character(n_src)
  out_a1     <- character(n_src)
  out_beta   <- numeric(n_src)

  # Helper: is a pair ambiguous?
  is_ambiguous <- function(a0, a1) {
    paste(a0, a1) %in% c("A T", "T A", "C G", "G C")
  }

  # Does (chr, pos) exist anywhere in the target?
  target_key <- paste(target_df$chr, target_df$pos)
  src_key    <- paste(sumstats_src$chr, sumstats_src$pos)

  # Build matched lookup keyed by source _NUM_ID_.
  # snp_match suffixes the shared _NUM_ID_ column: _NUM_ID_.ss for source,
  # _NUM_ID_ for target.  We need the SOURCE column (_NUM_ID_.ss).
  src_id_col <- if ("_NUM_ID_.ss" %in% names(snp_match_out)) "_NUM_ID_.ss" else "_NUM_ID_"
  matched_ids <- if (!is.null(snp_match_out) && nrow(snp_match_out) > 0) {
    snp_match_out[[src_id_col]]
  } else {
    integer(0)
  }

  for (i in seq_len(n_src)) {
    row <- sumstats_src[i, ]
    sid <- row[["_NUM_ID_"]]

    # 1. Position absent from target?
    if (!src_key[i] %in% target_key) {
      status[i]   <- "NO_MATCH_POS"
      out_a0[i]   <- "NA"
      out_a1[i]   <- "NA"
      out_beta[i] <- NA_real_
      next
    }

    # 2. Ambiguous (and strand-flip is enabled)?
    if (strand_flip_enabled && is_ambiguous(row$a0, row$a1)) {
      status[i]   <- "AMBIGUOUS"
      out_a0[i]   <- "NA"
      out_a1[i]   <- "NA"
      out_beta[i] <- NA_real_
      next
    }

    # 3. Was it matched?
    if (sid %in% matched_ids) {
      m <- snp_match_out[snp_match_out[[src_id_col]] == sid, ]
      # If duplicated on (chr, pos), snp_match removes all of them; treat
      # first as DUPLICATE (bigsnpr removes both; we label them DUPLICATE).
      # In practice our generator avoids this case in the hand-picked set.
      flip <- m[["_FLIP_"]]
      rev  <- m[["_REV_"]]
      if (flip && rev) {
        status[i] <- "FLIPPED_REVERSED"
      } else if (flip) {
        status[i] <- "FLIPPED"
      } else if (rev) {
        status[i] <- "REVERSED"
      } else {
        status[i] <- "EXACT"
      }
      out_a0[i]   <- m$a0
      out_a1[i]   <- m$a1
      out_beta[i] <- m$beta
    } else {
      # Position present but no candidate matched: NO_MATCH_ALLELE
      status[i]   <- "NO_MATCH_ALLELE"
      out_a0[i]   <- "NA"
      out_a1[i]   <- "NA"
      out_beta[i] <- NA_real_
    }
  }

  data.frame(status = status, out_a0 = out_a0, out_a1 = out_a1,
             out_beta = out_beta, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# SECTION 1: Hand-designed case
#
# Target: chr1, positions 100..700 in steps of 100
#   pos=100: A/G   pos=200: A/G   pos=300: A/G   pos=400: A/G
#   pos=500: A/T   pos=600: AT/A  pos=700: C/G
#
# Source (same positions):
#   pos=100: A/G  beta=0.1  -> EXACT
#   pos=200: G/A  beta=0.2  -> REVERSED  (beta -> -0.2)
#   pos=300: T/C  beta=0.3  -> FLIPPED
#   pos=400: C/T  beta=0.4  -> FLIPPED_REVERSED  (beta -> -0.4)
#   pos=500: A/T  beta=0.5  -> AMBIGUOUS (dropped)
#   pos=600: AT/A beta=0.6  -> EXACT (indel, not flippable)
#   pos=700: C/G  beta=0.7  -> AMBIGUOUS (dropped)
#
# Extra rows (not at these positions):
#   pos=800: A/G  beta=0.8  -> NO_MATCH_POS
#   pos=100: A/G  beta=0.9  -> (same pos as 100 in target, different alleles forced):
#     We use pos=900, alleles C/A where target has no C/A -> NO_MATCH_ALLELE
#     Actually add pos=100 to target but force alleles C/T there.
# Simpler: pos=900 is absent -> NO_MATCH_POS; for NO_MATCH_ALLELE use
#   source alleles that cannot match the target variant at an existing position.
# ---------------------------------------------------------------------------

# Target table for hand-designed case
tgt_hand <- data.frame(
  chr = 1L,
  pos = c(100L, 200L, 300L, 400L, 500L, 600L, 700L, 999L),
  a0  = c("A",  "A",  "A",  "A",  "A",  "AT", "C",  "C"),
  a1  = c("G",  "G",  "G",  "G",  "T",  "A",  "G",  "T"),
  stringsAsFactors = FALSE
)
tgt_hand[["_NUM_ID_"]] <- seq_len(nrow(tgt_hand))

# Source table for hand-designed case
#   pos=800: absent from target -> NO_MATCH_POS
#   pos=999: present in target (C/T), source says A/C -> NO_MATCH_ALLELE
#     (A/G would match via flip+rev: flip(A)=T, flip(G)=C, swap -> C/T = target)
#     A/C cannot match C/T by any of the four candidate forms.
src_hand <- data.frame(
  chr  = 1L,
  pos  = c(100L, 200L, 300L, 400L, 500L, 600L, 700L, 800L, 999L),
  a0   = c("A",  "G",  "T",  "C",  "A",  "AT", "C",  "A",  "A"),
  a1   = c("G",  "A",  "C",  "T",  "T",  "A",  "G",  "G",  "C"),
  beta = c(0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9),
  stringsAsFactors = FALSE
)
src_hand[["_NUM_ID_"]] <- seq_len(nrow(src_hand))

cat("\n=== Hand-designed case ===\n")
out_hand <- tryCatch(
  snp_match(src_hand, tgt_hand, return_flip_and_rev = TRUE),
  error = function(e) { cat("snp_match error:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(out_hand)) {
  cat("Matched output:\n")
  print(out_hand[, c("chr", "pos", "a0", "a1", "beta", "_FLIP_", "_REV_")])
}

class_hand <- classify_variants(src_hand, out_hand, tgt_hand)

# Sanity: verify ground truth
expected_status <- c("EXACT","REVERSED","FLIPPED","FLIPPED_REVERSED",
                      "AMBIGUOUS","EXACT","AMBIGUOUS","NO_MATCH_POS","NO_MATCH_ALLELE")
# Note: pos=600 (AT/A indel) is EXACT; it cannot be strand-flipped (indel alleles
# are multi-char) so the flip candidates produce NA and only the exact match applies.
expected_beta   <- c(0.1, -0.2, 0.3, -0.4, NA_real_, 0.6, NA_real_, NA_real_, NA_real_)

mismatches <- 0
for (i in seq_along(expected_status)) {
  got_s <- class_hand$status[i]
  exp_s <- expected_status[i]
  if (got_s != exp_s) {
    cat(sprintf("MISMATCH row %d: status got=%s expected=%s\n", i, got_s, exp_s))
    mismatches <- mismatches + 1
  }
  got_b <- class_hand$out_beta[i]
  exp_b <- expected_beta[i]
  if (is.na(exp_b) && is.na(got_b)) next
  if (is.na(exp_b) != is.na(got_b) || abs(got_b - exp_b) > 1e-12) {
    cat(sprintf("MISMATCH row %d: beta got=%.6g expected=%.6g\n", i, got_b, exp_b))
    mismatches <- mismatches + 1
  }
}
if (mismatches == 0) {
  cat("Hand-designed case: all statuses and betas match expected values.\n")
} else {
  stop(sprintf("%d mismatch(es) in hand-designed case", mismatches))
}

# ---------------------------------------------------------------------------
# SECTION 2: Randomised case
# ---------------------------------------------------------------------------
set.seed(20260802L)

# Allele alphabet for random SNPs
snp_allele_pairs <- list(
  c("A","G"), c("G","A"), c("T","C"), c("C","T"),
  c("A","C"), c("C","A"), c("T","G"), c("G","T"),
  c("A","T"), c("T","A"), c("C","G"), c("G","C")
)
indel_pairs <- list(c("AT","A"), c("A","AT"), c("GT","G"), c("G","GT"),
                    c("ACT","A"), c("A","ACT"))

N_TGT <- 300L
N_SRC <- 350L

# Build target
rand_pos <- sort(sample(seq(1000L, 99000L, by = 100L), N_TGT))
rand_tgt_alleles <- lapply(seq_len(N_TGT), function(i) {
  if (runif(1) < 0.05) indel_pairs[[sample(length(indel_pairs), 1)]]
  else snp_allele_pairs[[sample(12, 1)]]
})
tgt_rand <- data.frame(
  chr = sample(c(1L, 2L, 3L), N_TGT, replace = TRUE),
  pos = rand_pos,
  a0  = sapply(rand_tgt_alleles, `[`, 1),
  a1  = sapply(rand_tgt_alleles, `[`, 2),
  stringsAsFactors = FALSE
)
# Ensure (chr, pos) uniqueness in target
tgt_rand <- tgt_rand[!duplicated(paste(tgt_rand$chr, tgt_rand$pos)), ]
tgt_rand[["_NUM_ID_"]] <- seq_len(nrow(tgt_rand))
N_TGT_ACTUAL <- nrow(tgt_rand)

# Build source variants, mixing all outcome types
# ~60% overlap by position with the target
n_overlap <- as.integer(N_TGT_ACTUAL * 0.6)
overlap_idx <- sample(seq_len(N_TGT_ACTUAL), n_overlap)
n_pos_only <- N_SRC - n_overlap  # no-pos variants

src_rows <- vector("list", N_SRC)
for (k in seq_len(n_overlap)) {
  ti <- overlap_idx[k]
  trow <- tgt_rand[ti, ]
  # Decide outcome type randomly
  otype <- sample(c("exact","reversed","flipped","flip_rev","ambiguous","no_allele"),
                  1, prob = c(0.25, 0.20, 0.15, 0.15, 0.10, 0.15))

  flip_strand_single <- function(a) {
    switch(a, A = "T", C = "G", T = "A", G = "C", NA_character_)
  }

  if (otype == "exact") {
    a0 <- trow$a0; a1 <- trow$a1
  } else if (otype == "reversed") {
    a0 <- trow$a1; a1 <- trow$a0
  } else if (otype == "flipped") {
    fa0 <- flip_strand_single(trow$a0)
    fa1 <- flip_strand_single(trow$a1)
    if (is.na(fa0) || is.na(fa1)) {
      otype <- "exact"; a0 <- trow$a0; a1 <- trow$a1
    } else {
      # Check not ambiguous after flip
      if (paste(fa0, fa1) %in% c("A T", "T A", "C G", "G C")) {
        otype <- "exact"; a0 <- trow$a0; a1 <- trow$a1
      } else {
        a0 <- fa0; a1 <- fa1
      }
    }
  } else if (otype == "flip_rev") {
    fa0 <- flip_strand_single(trow$a0)
    fa1 <- flip_strand_single(trow$a1)
    if (is.na(fa0) || is.na(fa1)) {
      otype <- "reversed"; a0 <- trow$a1; a1 <- trow$a0
    } else {
      if (paste(fa1, fa0) %in% c("A T", "T A", "C G", "G C")) {
        otype <- "reversed"; a0 <- trow$a1; a1 <- trow$a0
      } else {
        a0 <- fa1; a1 <- fa0
      }
    }
  } else if (otype == "ambiguous") {
    # Pick a genuine ambiguous pair
    amb <- list(c("A","T"), c("T","A"), c("C","G"), c("G","C"))[[sample(4,1)]]
    a0 <- amb[1]; a1 <- amb[2]
  } else {  # no_allele
    # Pick alleles that cannot match by any of the four candidates
    repeat {
      cand <- snp_allele_pairs[[sample(12, 1)]]
      # Check that none of the 4 candidates match target
      fa0t <- flip_strand_single(trow$a0)
      fa1t <- flip_strand_single(trow$a1)
      pairs_to_avoid <- list(
        c(trow$a0, trow$a1), c(trow$a1, trow$a0)
      )
      if (!is.na(fa0t) && !is.na(fa1t)) {
        pairs_to_avoid <- c(pairs_to_avoid,
                            list(c(fa0t, fa1t), c(fa1t, fa0t)))
      }
      avoid <- any(sapply(pairs_to_avoid, function(p)
        cand[1] == p[1] && cand[2] == p[2]))
      if (!avoid) { a0 <- cand[1]; a1 <- cand[2]; break }
    }
  }
  src_rows[[k]] <- list(chr = trow$chr, pos = trow$pos,
                         a0 = a0, a1 = a1,
                         beta = round(rnorm(1), 4))
}

# Add position-absent variants
all_target_pos <- paste(tgt_rand$chr, tgt_rand$pos)
absent_rows <- 0
for (k in (n_overlap + 1):N_SRC) {
  repeat {
    ch  <- sample(1:3, 1)
    pos <- sample(c(1L:999L, 100000L:200000L), 1) * 100L
    if (!paste(ch, pos) %in% all_target_pos) break
  }
  src_rows[[k]] <- list(chr = ch, pos = pos,
                         a0 = snp_allele_pairs[[sample(12,1)]][1],
                         a1 = snp_allele_pairs[[sample(12,1)]][2],
                         beta = round(rnorm(1), 4))
  absent_rows <- absent_rows + 1
}

src_rand <- do.call(rbind, lapply(src_rows, as.data.frame, stringsAsFactors = FALSE))
src_rand$chr  <- as.integer(src_rand$chr)
src_rand$pos  <- as.integer(src_rand$pos)
src_rand$beta <- as.numeric(src_rand$beta)
src_rand[["_NUM_ID_"]] <- seq_len(nrow(src_rand))

cat("\n=== Randomised case ===\n")
cat("Target:", nrow(tgt_rand), "variants\n")
cat("Source:", nrow(src_rand), "variants\n")

out_rand <- tryCatch(
  snp_match(src_rand, tgt_rand, return_flip_and_rev = TRUE),
  error = function(e) { cat("snp_match error:", conditionMessage(e), "\n"); NULL }
)

if (is.null(out_rand)) out_rand <- src_rand[integer(0), ]
cat("Matched:", nrow(out_rand), "\n")

class_rand <- classify_variants(src_rand, out_rand, tgt_rand)
cat("Status breakdown:\n"); print(table(class_rand$status))

# ---------------------------------------------------------------------------
# Write fixture file
# ---------------------------------------------------------------------------
out_path <- "Tests/pgs_match_fixtures/match_reference.txt"

write_section <- function(f, section_name, src_df, class_df) {
  cat(sprintf("# section %s\n", section_name), file = f)
  cat(sprintf("# n_variants %d\n", nrow(src_df)), file = f)
  cat("# columns: chr pos src_a0 src_a1 src_beta status out_a0 out_a1 out_beta\n",
      file = f)
  for (i in seq_len(nrow(src_df))) {
    ob <- class_df$out_beta[i]
    ob_str <- if (is.na(ob)) "NA" else sprintf("%.17g", ob)
    cat(sprintf("%d %d %s %s %.17g %s %s %s %s\n",
                src_df$chr[i], src_df$pos[i],
                src_df$a0[i],  src_df$a1[i],
                src_df$beta[i],
                class_df$status[i],
                class_df$out_a0[i], class_df$out_a1[i],
                ob_str),
        file = f)
  }
}

f <- file(out_path, "w")
cat("# adjutant pgs_match reference -- generated by generate_match.R, do not edit\n",
    file = f)
cat(sprintf("# bigsnpr %s\n", as.character(packageVersion("bigsnpr"))), file = f)
cat("# Two sections: hand (hand-designed) and rand (randomised, set.seed(20260802)).\n",
    file = f)
cat("# Statuses: EXACT REVERSED FLIPPED FLIPPED_REVERSED AMBIGUOUS ",
    "NO_MATCH_POS NO_MATCH_ALLELE\n", sep = "", file = f)
write_section(f, "hand", src_hand, class_hand)
write_section(f, "rand", src_rand, class_rand)
close(f)

cat(sprintf("\nWrote %s\n", out_path))
cat("Fixture generation complete.\n")
