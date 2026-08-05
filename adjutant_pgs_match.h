#ifndef __ADJUTANT_PGS_MATCH_H__
#define __ADJUTANT_PGS_MATCH_H__

// Per-variant harmonisation of a PGS weight file against a target .bim.
//
// Reimplements bigsnpr::snp_match() as a pure function so it can be
// exercised directly from pgs_eval_test against committed fixtures without
// any file I/O.
//
// Reference: Privé et al., bigsnpr R package, snp_match().

#include <stdint.h>

#include "adjutant_mhc.h"
#include "include/plink2_base.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// -------------------------------------------------------------------
// Status codes for a single source variant after harmonisation.
//
// kPgsMatchNoPos and kPgsMatchDuplicate are decided by the caller
// (position lookup and dedup); they live in the enum but are NOT
// returned by PgsMatchOne.
//
// kPgsMatchResolvedExact/Reversed are deliberately distinct from plain
// Exact/Reversed: a palindrome resolved by comparing allele frequencies
// is a weaker inference than an allele match, and the report should not
// hide that.
// -------------------------------------------------------------------
typedef enum {
  kPgsMatchExact,
  kPgsMatchReversed,
  kPgsMatchFlipped,
  kPgsMatchFlippedReversed,
  kPgsMatchAmbiguous,
  kPgsMatchResolvedExact,
  kPgsMatchResolvedReversed,
  kPgsMatchAmbiguousUnresolved,
  kPgsMatchNoPos,
  kPgsMatchNoAllele,
  kPgsMatchDuplicate
} PgsMatchStatus;

// -------------------------------------------------------------------
// Frequency inputs for resolving one strand-ambiguous (palindromic)
// variant.  Passed as nullptr when --resolve-ambiguous is absent.
//
//   maf_threshold  refuse to resolve when min(f, 1 - f) reaches this for
//                  EITHER frequency, because the hypotheses are then
//                  indistinguishable.  The signal is |2f - 1| on each
//                  side, so a target frequency near 0.5 is exactly as
//                  uninformative as a source frequency near 0.5.
//   src_freq       effect-allele (a1) frequency in the SOURCE study.
//   tgt_a1_freq    frequency of the TARGET's a1 allele.  NaN when the
//                  target frequency is unavailable for this variant.
// -------------------------------------------------------------------
typedef struct PgsAmbigResolveStruct {
  double maf_threshold;
  double src_freq;
  double tgt_a1_freq;
} PgsAmbigResolve;

// -------------------------------------------------------------------
// Harmonise one source variant against one target variant.
//
// Pure: no I/O, no allocation.
//
// Convention (matching bigsnpr):
//   a0 = other allele, a1 = effect allele, beta is w.r.t. a1.
//
// Four candidate representations are tried in order:
//   1. (src_a0, src_a1)                  -> kPgsMatchExact,           reversed=0
//   2. (src_a1, src_a0)                  -> kPgsMatchReversed,        reversed=1
//   3. (comp(src_a0), comp(src_a1))      -> kPgsMatchFlipped,         reversed=0
//   4. (comp(src_a1), comp(src_a0))      -> kPgsMatchFlippedReversed, reversed=1
// The first candidate equal to (tgt_a0, tgt_a1) wins.
// comp() returns a failure marker for anything that is not a single
// A/C/T/G (case-sensitive), so indels never match through flipped forms,
// replicating flip_strand which returns NA for anything else.
//
// Ambiguous check: src allele pair in {AT, TA, CG, GC} when
// !allow_ambiguous -> kPgsMatchAmbiguous (no candidate tried).
//
// When `resolve` is non-null the ambiguous pair is instead decided by
// frequency.  For a palindrome the exact and flipped+reversed forms
// carry identical allele labels, as do reversed and flipped, so there
// are exactly TWO distinguishable outcomes: the source effect allele is
// the target's a1, or it is the target's a0.
//   - allele pair sets differ  -> kPgsMatchNoAllele  (an A/T source
//     against a C/G target is palindromic on both sides but shares no
//     allele; frequencies must not be allowed to "resolve" it)
//   - either frequency absent           -> kPgsMatchAmbiguousUnresolved
//   - min(f, 1-f) >= threshold for the SOURCE or the TARGET frequency
//                                       -> kPgsMatchAmbiguousUnresolved
//   - src_freq nearer tgt_a1_freq       -> kPgsMatchResolvedExact
//   - src_freq nearer 1 - tgt_a1_freq   -> kPgsMatchResolvedReversed
//
// The threshold comparison carries a few ulps of slack, so that the two
// decimal spellings of one MAF (0.4 and 0.6) decide the same way through
// a parser that is only accurate to ~2 bits; see FreqTooCloseToHalf.
//
// The two distances are equal only when a frequency is exactly 0.5, and
// a MAF of 0.5 is at or above every permitted threshold, so the two
// guards above make the tie unreachable -- there is no input for which
// the sign is decided on zero evidence.
//
// *out_reversed is set to 1 for the reversed outcomes (including
// kPgsMatchResolvedReversed); the CALLER negates beta.  The sign
// decision is kept at one visible call site so it cannot be applied
// twice.
//
// Returns kPgsMatchNoAllele when no candidate matches.
// -------------------------------------------------------------------
PgsMatchStatus PgsMatchOne(const char* src_a0, const char* src_a1,
                           const char* tgt_a0, const char* tgt_a1,
                           uint32_t allow_ambiguous,
                           const PgsAmbigResolve* resolve,
                           uint32_t* out_reversed);

// -------------------------------------------------------------------
// Full I/O: read a .bim and a PGS weight file, harmonise, write
// <out_prefix>.matched and <out_prefix>.match_report.
//
// Parameters:
//   bim_fname         path to the target .bim (CHR ID CM POS A1 A2, no header)
//   pgs_model_fname   path to the weight file (possibly .zst)
//   col_chr           1-based column index for CHR in the weight file
//   col_pos           1-based column index for POS in the weight file
//   col_effect        1-based column index for the effect (a1) allele
//   col_other         1-based column index for the other (a0) allele
//   col_weight        1-based column index for the beta/weight
//   col_freq          1-based column index for the source effect-allele
//                     frequency, or 0 when absent
//   allow_ambiguous   if nonzero, keep A/T and C/G variants
//   resolve_ambig     if nonzero, resolve palindromes by frequency;
//                     requires col_freq and read_freq_fname
//   maf_threshold     refuse to resolve at or above this MAF, on the
//                     source and on the target frequency alike
//   read_freq_fname   target .afreq, or nullptr
//   outname           output filename prefix (NUL-terminated)
//   outname_end       pointer past the last character of the prefix
// -------------------------------------------------------------------
PglErr RunPgsModel(const char* bim_fname,
                   const char* pgs_model_fname,
                   uint32_t col_chr,
                   uint32_t col_pos,
                   uint32_t col_effect,
                   uint32_t col_other,
                   uint32_t col_weight,
                   uint32_t col_freq,
                   uint32_t allow_ambiguous,
                   uint32_t resolve_ambig,
                   double maf_threshold,
                   const char* read_freq_fname,
                   MhcRegion mhc_region,
                   char* outname,
                   char* outname_end);

// -------------------------------------------------------------------
// Merge several PGS models into one matrix.
//
// Reads a file of model paths, harmonises every one against the SAME
// .bim through the same matching code --pgs-model uses, and writes
// <out_prefix>.matched as
//
//   ID  A1  <model1>  <model2> ...
//
// one row per variant matched by at least one model, in .bim order,
// with 0 where a model does not carry that variant.
//
// Harmonising before merging is the point of the feature, not an
// implementation detail: two models may name opposite effect alleles
// for the same variant, so merging raw weight files on variant ID
// would put one column's beta on the wrong allele and nothing
// downstream would notice.  Every column here is w.r.t. the .bim's A1.
//
// Reorienting a weight shifts the absolute score by a constant, and it
// is worth knowing why that is acceptable rather than assuming it.  For
// a biallelic variant dosage(other) == 2 - dosage(effect), so
//
//     beta * dosage(effect) == 2*beta - beta * dosage(other)
//
// and writing -beta against the other allele reproduces the original
// contribution MINUS 2*beta.  Summed over the reversed variants that is
// one constant, identical for every sample: rankings, correlations, R^2
// and any regression with an intercept are unchanged, and only raw
// SCORE_SUM values move.  Measured on 5 variants at |beta| = 0.4, every
// sample's sum differed by exactly 4.0 with r = 1.0 over 200 samples.
//
// The identity holds for ADDITIVE models only, which is all adjutant
// scores.  It fails for dominant/recessive effect types -- which is why
// tools supporting those must split a model set across several files
// instead of reorienting, and why one merged file always suffices here.
//
// <out_prefix>.match_report is the --pgs-model report with a leading
// MODEL column, concatenated across models.
//
// Parameters are as RunPgsModel, except:
//   list_fname   file of model paths, one per line.  Lines whose first
//                non-space character is '#', and whitespace-only lines,
//                are skipped.  A COMPLETELY empty interior line is a
//                hard error, refused by the text reader itself rather
//                than here (plink2_text tolerates trailing empty lines
//                but not interior ones, by deliberate API choice) -- the
//                message names the file, so it is not a silent failure.
//                The column indices apply to EVERY file.
//
// Model names are the basename up to the first '.', so
// "PGS000123.txt.gz" -> "PGS000123".  Two paths yielding the same name
// is an error.  A model matching nothing is a warning and its column
// is dropped, since an all-zero column is a constant predictor.
// -------------------------------------------------------------------
PglErr RunPgsModelsList(const char* bim_fname,
                        const char* list_fname,
                        uint32_t col_chr,
                        uint32_t col_pos,
                        uint32_t col_effect,
                        uint32_t col_other,
                        uint32_t col_weight,
                        uint32_t col_freq,
                        uint32_t allow_ambiguous,
                        uint32_t resolve_ambig,
                        double maf_threshold,
                        const char* read_freq_fname,
                        MhcRegion mhc_region,
                        char* outname,
                        char* outname_end);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __ADJUTANT_PGS_MATCH_H__
