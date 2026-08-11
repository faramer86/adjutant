// This file is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "include/SFMT.h"
#include "include/pgenlib_misc.h"
#include "include/pgenlib_read.h"
#include "include/plink2_base.h"
#include "include/plink2_bits.h"
#include "include/plink2_float.h"
#include "include/plink2_string.h"
#include "include/plink2_thread.h"
#include "include/plink2_zstfile.h"
#include "adjutant_cmdline.h"
#include "adjutant_common.h"
#include "adjutant_compress_stream.h"
#include "adjutant_data.h"
#include "adjutant_filter.h"
#include "adjutant_help.h"
#include "adjutant_matrix.h"
#include "adjutant_matrix_calc.h"
#include "adjutant_mhc.h"
#include "adjutant_pgs_eval.h"
#include "adjutant_pgs_match.h"
#include "adjutant_psam.h"
#include "adjutant_pvar.h"

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>  // time()
#include <unistd.h>  // unlink()

#ifdef __cplusplus
namespace plink2 {
#endif

static PREFER_CONSTEXPR char ver_str[] = "adjutant v1.0.0"
#ifdef NOLAPACK
  "NL"
#elif defined(LAPACK_ILP64)
  "LM"
#endif

#ifdef __LP64__
#  ifdef USE_CUDA
    " CUDA"
#  elif defined(USE_AVX2)
    " AVX2"
#  elif defined(__APPLE__) && !defined(__x86_64__)
    " M1"
#  elif defined(USE_SSE42)
    " SSE4.2"
#  else
    " 64-bit"
#  endif
#else
  " 32-bit"
#endif

#ifdef USE_MKL
  " Intel"
#elif defined(USE_AOCL)
  " AMD"
#endif
  " (3 Jun 2026)";
static PREFER_CONSTEXPR char ver_str2[] =
  // include leading space if day < 10, so character length stays the same
  " "

#ifdef NOLAPACK
#elif defined(LAPACK_ILP64)
#else
  "  "
#endif

#ifdef USE_MKL
#elif defined(USE_AOCL)
  "  "
#else
  "      "
#endif

#ifdef __LP64__
#  ifdef USE_CUDA
    "  "
#  elif defined(USE_AVX2)
    "  "
#  elif defined(__APPLE__) && !defined(__x86_64__)
    "    "
#  endif
#endif

  "    PLINK v2.0.0-a.7.1.a derivative\n"
  "(C) 2005-2026 Shaun Purcell, Christopher Chang    GNU General Public License v3\n";
#ifdef HAS_CONSTEXPR
static_assert(CompileTimeSlen(ver_str) + CompileTimeSlen(ver_str2) == 160, "ver_str/ver_str2 must be updated");
#endif
static const char errstr_append[] = "For more info, try \"" PROG_NAME_STR " --help <flag name>\" or \"" PROG_NAME_STR " --help | more\".\n";

static const char notestr_null_calc2[] = "The only command in this build is --score.\n\n\"" PROG_NAME_STR " --help | more\" describes all functions.\n";

// multiallelics-already-joined + terminating null
CONSTI32(kMaxFlagBlen, 29);

FLAGSET_DEF_START()
  kfLoadParams0,
  kfLoadParamsPgen = (1 << 0),
  kfLoadParamsPsam = (1 << 1),
  kfLoadParamsPvar = (1 << 2),
  kfLoadParamsPfileAll = (kfLoadParamsPgen | kfLoadParamsPsam | kfLoadParamsPvar)
FLAGSET_DEF_END(LoadParams);

FLAGSET_DEF_START()
  kfXload0,
  kfXloadVcf = (1 << 0),
  kfXloadBcf = (1 << 1),
  kfXloadOxSample = (1 << 2),
  kfXloadOxGen = (1 << 3),
  kfXloadOxBgen = (1 << 4),
  kfXloadOxHaps = (1 << 5),
  kfXloadOxLegend = (1 << 6),
  kfXloadPlink1Dosage = (1 << 7),
  kfXloadMap = (1 << 8),
  kfXloadGenDummy = (1 << 9),
  kfXloadPed = (1 << 10),
  kfXloadTped = (1 << 11),
  kfXloadEigGeno = (1 << 12),
  kfXloadEigInd = (1 << 13),
  kfXloadEigSnp = (1 << 14)
FLAGSET_DEF_END(Xload);


// assume for now that .pgen must always be accompanied by both .pvar and .psam
// currently does double-duty in tracking dependencies
FLAGSET64_DEF_START()
  kfFilter0,
  kfFilterAllReq = (1 << 0),
  kfFilterPvarReq = (1 << 1),
  kfFilterPsamReq = (1 << 2),
  kfFilterOpportunisticPgen = (1 << 3),
  kfFilterNonrefFlagsNeeded = (1 << 4),
  kfFilterNonrefFlagsNeededSet = kfFilterNonrefFlagsNeeded | kfFilterOpportunisticPgen,
  kfFilterExclFemales = (1 << 6),
  kfFilterExclMales = (1 << 7),
  kfFilterExclNosex = (1 << 8),
  kfFilterExclFounders = (1 << 9),
  kfFilterExclNonfounders = (1 << 10),
  kfFilterSnpsOnlyJustAcgt = (1 << 11),
  kfFilterSelectSidRepresentatives = (1 << 18),
  kfFilterMendel = (1 << 19)
FLAGSET64_DEF_END(FilterFlags);

FLAGSET64_DEF_START()
  kfCommand10,
  kfCommand1MakePlink2 = (1 << 0),
  kfCommand1Exportf = (1 << 1),
  kfCommand1AlleleFreq = (1 << 6),
  kfCommand1GenoCounts = (1 << 7),
  kfCommand1Hardy = (1 << 8),
  kfCommand1Pca = (1 << 10),
  kfCommand1Score = (1 << 15),
  kfCommand1WriteCovar = (1 << 16),
  kfCommand1PgenInfo = (1 << 19),
  kfCommand1Pmerge = (1 << 26),
  kfCommand1Vcor = (1 << 29),
  kfCommand1CheckOrImputeSex = (1U << 31),
  kfCommand1MendelReport = (1LLU << 32)
FLAGSET64_DEF_END(Command1Flags);

void PgenInfoPrint(const char* pgenname, const PgenFileInfo* pgfip, PgenExtensionLl* header_exts, PgenHeaderCtrl header_ctrl, uint32_t max_allele_ct) {
  logprintfww("--pgen-info on %s:\n", pgenname);
  logprintf("  Variants: %u\n", pgfip->raw_variant_ct);
  logprintf("  Samples: %u\n", pgfip->raw_sample_ct);
  const uint32_t nonref_flags_status = header_ctrl >> 6;
  if (!nonref_flags_status) {
    logputs("  REF allele known/provisional status not stored in .pgen\n");
  } else if (nonref_flags_status == 1) {
    logputs("  REF alleles are all known\n");
  } else if (nonref_flags_status == 2) {
    logputs("  REF alleles are all provisional\n");
  } else {
    // could report exact counts of each
    logputs("  REF alleles are a mix of known and provisional\n");
  }
  if (max_allele_ct >= UINT32_MAXM1) {
    if (max_allele_ct == UINT32_MAX) {
      logputs("  Maximum allele count for a single variant: >2, not explicitly stored\n");
    } else {
      logputs("  Maximum allele count for a single variant: not explicitly stored\n");
    }
  } else {
    logprintf("  Maximum allele count for a single variant: %u\n", max_allele_ct);
  }
  if (pgfip->gflags & kfPgenGlobalHardcallPhasePresent) {
    logputs("  Explicitly phased hardcalls present\n");
  } else {
    logputs("  No hardcalls are explicitly phased\n");
  }
  if (pgfip->gflags & kfPgenGlobalDosagePresent) {
    if (pgfip->gflags & kfPgenGlobalDosagePhasePresent) {
      logputs("  Explicitly phased dosages present\n");
    } else {
      logputs("  Dosage present, none explicitly phased\n");
    }
  } else {
    logputs("  No dosages present\n");
  }
  assert(header_exts && (header_exts->type_idx == 1));
  if (header_exts->size != ~0LLU) {
    char* write_iter = strcpya_k(g_logbuf, "  Writer: ");
    write_iter = memcpya(write_iter, header_exts->contents, header_exts->size);
    strcpy_k(write_iter, "\n");
    WordWrapB(0);
    logputsb();
  }
}

PglErr PgenInfoStandalone(const char* pgenname, const char* pginame) {
  PgenFileInfo pgfi;
  PglErr reterr = kPglRetSuccess;
  PreinitPgfi(&pgfi);
  PgenExtensionLl ext_slot; // shouldn't have shorter lifetime than pgfi
  ext_slot.contents = nullptr;
  {
    PgenHeaderCtrl header_ctrl;
    uintptr_t cur_alloc_cacheline_ct;
    reterr = PgfiInitPhase1(pgenname, pginame, UINT32_MAX, UINT32_MAX, &header_ctrl, &pgfi, &cur_alloc_cacheline_ct, g_logbuf);
    if (unlikely(reterr)) {
      if ((reterr == kPglRetSampleMajorBed) || (reterr == kPglRetImproperFunctionCall)) {
        logerrputs("Warning: Skipping --pgen-info since a .bed file was provided.\n");
        reterr = kPglRetSuccess;
      } else {
        WordWrapB(0);
        logerrputsb();
      }
      goto PgenInfoStandalone_ret_1;
    }
    const uint32_t raw_variant_ct = pgfi.raw_variant_ct;
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    unsigned char* pgfi_alloc;
    if (unlikely(bigstack_alloc_uc(cur_alloc_cacheline_ct * kCacheline, &pgfi_alloc) ||
                 bigstack_alloc_w(raw_variant_ct + 1, &pgfi.allele_idx_offsets) ||
                 bigstack_alloc_w(raw_variant_ctl, &pgfi.nonref_flags))) {
      goto PgenInfoStandalone_ret_NOMEM;
    }
    ext_slot.next = nullptr;
    ext_slot.size = ~0LLU;
    ext_slot.type_idx = 1;
    PgenExtensionLl* header_exts = &ext_slot;
    uintptr_t pgr_alloc_cacheline_ct = 0;
    uint32_t max_vrec_width;
    reterr = PgfiInitPhase2Ex(header_ctrl, 0, 0, 1, 0, raw_variant_ct, &max_vrec_width, &pgfi, pgfi_alloc, header_exts, nullptr, &pgr_alloc_cacheline_ct, g_logbuf);
    if (unlikely(reterr)) {
      WordWrapB(0);
      logerrputsb();
      goto PgenInfoStandalone_ret_1;
    }
    if (header_exts->size != ~0LLU) {
#ifndef __LP64__
      if (unlikely(header_exts->size >= (1LLU << 31))) {
        goto PgenInfoStandalone_ret_NOMEM;
      }
#endif
      if (unlikely(pgl_malloc(header_exts->size, &(header_exts->contents)))) {
        goto PgenInfoStandalone_ret_NOMEM;
      }
      reterr = PgfiInitLoadExts(header_ctrl, &pgfi, header_exts, nullptr, g_logbuf);
      if (unlikely(reterr)) {
        WordWrapB(0);
        logerrputsb();
        goto PgenInfoStandalone_ret_1;
      }
    }
    uint32_t max_allele_ct = 2;
    if (pgfi.gflags & kfPgenGlobalMultiallelicHardcallFound) {
      max_allele_ct = UINT32_MAX;
    } else if (pgfi.gflags & kfPgenGlobalDosagePresent) {
      max_allele_ct = UINT32_MAXM1;
    }
    PgenInfoPrint(pgenname, &pgfi, header_exts, header_ctrl, max_allele_ct);
  }
  while (0) {
  PgenInfoStandalone_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
 PgenInfoStandalone_ret_1:
  free_cond(ext_slot.contents);
  CleanupPgfi2(pgenname, &pgfi, &reterr);
  // no BigstackReset() needed?
  return reterr;
}

typedef struct Plink2CmdlineStruct {
  NONCOPYABLE(Plink2CmdlineStruct);
  MiscFlags misc_flags;

  // filter_flags tracks some info about flags which may cause the .pgen part
  // of --make-pgen's output to change (affects sample_include/variant_include,
  // etc.)  other .psam/.pvar/.pgen dependency info is now tracked by
  // dependency_flags (which is set to a superset of filter_flags at the end of
  // command-line parsing).
  // (filter_flags & kfFilterPsamReq) and (filter_flags & kfFilterPvarReq) can
  // now be used to detect whether sample_include or variant_include may be
  // modified.
  FilterFlags filter_flags;
  FilterFlags dependency_flags;

  // .log should show where --chr, --exclude-if-info, etc. are applied
  LoadFilterLogFlags load_filter_log_flags;

  Command1Flags command_flags1;
  PvarPsamFlags pvar_psam_flags;

  // --exclude-mhc-hg19 / --exclude-xmhc-hg19 / --exclude-mhc-hg38 /
  // --exclude-xmhc-hg38.  At most one; kMhcNone when none was given.
  MhcRegion mhc_region;

  RangeList snps_range_list;
  RangeList exclude_snps_range_list;
  RangeList pheno_range_list;
  RangeList covar_range_list;
  FamCol fam_cols;
  STD_ARRAY_DECL(FreqFilterMode, 4, filter_modes);
  SelectSidMissingnessMode select_sid_missingness_mode;
  SelectSidTiebreakMode select_sid_tiebreak_mode;
  ScoreInfo score_info;
  CmpExpr keep_if_expr;
  CmpExpr remove_if_expr;
  CmpExpr extract_if_info_expr;
  CmpExpr exclude_if_info_expr;
  ExtractColCondInfo extract_col_cond_info;
  float var_min_qual;
  uint32_t splitpar_bound1;
  uint32_t splitpar_bound2;
  uint32_t new_variant_id_max_allele_slen;
  // maybe support BGEN v1.2-style variable-precision dosages later, at which
  // point these should be floating-point numbers; but let's first see what we
  // gain from v1.1 fixed-point arithmetic
  uint32_t hard_call_thresh;
  uint32_t dosage_erase_thresh;
  double geno_thresh;
  double hwe_ln_thresh;
  double hwe_sample_size_term; // -1 if unspecified
  double mach_r2_max;
  double minimac3_r2_max;
  double af_pseudocount;
  double min_maf;
  double max_maf;
  double pgs_prevalence; // < 0 means not supplied
  uint32_t pgs_folds;      // --folds; 1 means no cross-validation
  uint32_t pgs_fold_seed;  // --seed, kept for reproducible fold assignment
  uint64_t min_allele_ddosage;
  uint64_t max_allele_ddosage;
  int32_t missing_pheno;
  int32_t window_bp;
  uint32_t xchr_model;
  uint32_t max_thread_ct;
  uint32_t keep_col_match_num;
  uint32_t filter_min_allele_ct;
  uint32_t filter_max_allele_ct;
  char input_missing_geno_char;
  // might want input_missing_catname and/or output_missing_catname later, but
  // let's start with the simplest implementation
  char missing_catname[kMaxMissingPhenostrBlen];
  char output_missing_pheno[kMaxMissingPhenostrBlen];
  char legacy_output_missing_pheno[kMaxMissingPhenostrBlen];

  char* pginame;
  char* var_filter_exceptions_flattened;
  char* varid_template_str;
  char* varid_multi_template_str;
  char* varid_multi_nonsnp_template_str;
  char* missing_varid_match;
  char* varid_from;
  char* varid_to;
  char* varid_snp;
  char* varid_exclude_snp;
  char* pheno_fname;
  char* covar_fname;
  char* extract_fnames;
  char* extract_intersect_fnames;
  char* exclude_fnames;
  char* keep_fnames;
  char* keepfam_fnames;
  char* remove_fnames;
  char* removefam_fnames;
  char* sample_sort_fname;
  char* freq_ref_binstr;
  char* freq_alt1_binstr;
  char* glm_local_covar_fname;
  char* glm_local_pvar_fname;
  char* glm_local_psam_fname;
  char* read_freq_fname;
  char* within_fname;
  char* catpheno_name;
  char* family_missing_catname;
  char* keep_cats_fname;
  char* keep_cat_names_flattened;
  char* keep_cat_phenoname;
  char* remove_cats_fname;
  char* remove_cat_names_flattened;
  char* remove_cat_phenoname;
  char* split_cat_phenonames_flattened;
  char* require_pheno_flattened;
  char* require_covar_flattened;
  char* vstd_flattened;
  char* quantnorm_flattened;
  char* covar_quantnorm_flattened;
  char* loop_cats_phenoname;
  char* fa_fname;
  char* require_info_flattened;
  char* require_no_info_flattened;
  char* keep_col_match_fname;
  char* keep_col_match_flattened;
  char* keep_col_match_name;
  char* update_sample_ids_fname;
  char* update_parental_ids_fname;
  char* recover_var_ids_fname;
  char* vscore_fname;
  char* sscore_fname;  // --sscore: user-supplied pre-computed score file
  char* bim_fname;         // --bim: target .bim for --pgs-model
  char* pgs_model_fname;   // --pgs-model: weight file
  // --pgs-models-list: file of weight-file paths, merged into one matrix.
  // Alternative to pgs_model_fname; the column indices below serve both,
  // and so does every --bim/--allow-ambiguous/--resolve-ambiguous guard.
  char* pgs_models_list_fname;
  char* pgs_pheno_pairs_fname;  // --pgs-pheno-pairs
  uint32_t pgs_model_by_id;   // --pgs-model[s-list] id ...
  uint32_t pgs_model_col_id;
  uint32_t pgs_model_col_chr;     // 1-based column indices
  uint32_t pgs_model_col_pos;
  uint32_t pgs_model_col_effect;
  uint32_t pgs_model_col_other;
  uint32_t pgs_model_col_weight;
  uint32_t pgs_model_col_freq;    // optional 6th index; 0 when absent
  uint32_t pgs_model_allow_ambiguous;  // --allow-ambiguous flag
  uint32_t pgs_model_resolve_ambig;    // --resolve-ambiguous flag
  double pgs_model_maf_threshold;      // --resolve-ambiguous argument
  char* indep_preferred_fname;
  char* not_pheno_flattened;
  char* not_covar_flattened;
  char* indv_str;
  char* rename_chrs_fname;
  char* zero_cluster_fname;
  char* zero_cluster_phenoname;
  TwoColParams* ref_allele_flag;
  TwoColParams* alt_allele_flag;
  TwoColParams* update_chr_flag;
  TwoColParams* update_map_flag;
  TwoColParams* update_name_flag;
} Plink2Cmdline;

// er, probably time to just always initialize this...
uint32_t SingleVariantLoaderIsNeeded(Command1Flags command_flags1, double hwe_ln_thresh) {
  return (command_flags1 & kfCommand1Score) || (hwe_ln_thresh != -DBL_MAX);
}


uint32_t DecentAlleleFreqsAreNeeded(Command1Flags command_flags1, ScoreFlags score_flags) {
  // --score mean-imputes missing genotypes unless told otherwise, and
  // center/variance-standardize both need the frequencies too.
  return (command_flags1 & kfCommand1Score) && ((!(score_flags & kfScoreNoMeanimpute)) || (score_flags & (kfScoreCenter | kfScoreVarianceStandardize)));
}

// only needs to cover cases not captured by DecentAlleleFreqsAreNeeded()
uint32_t IndecentAlleleFreqsAreNeeded(double min_maf, double max_maf) {
  return (min_maf != 0.0) || (max_maf != 1.0);
}


uint32_t GetFirstHaploidUidx(const ChrInfo* cip, UnsortedVar vpos_sortstatus) {
  // returns 0x7fffffff if no X/haploid chromosomes present
  if (!(vpos_sortstatus & kfUnsortedVarSplitChr)) {
    const uint32_t chr_ct = cip->chr_ct;
    for (uint32_t chr_fo_idx = 0; chr_fo_idx != chr_ct; ++chr_fo_idx) {
      const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      if (IsSet(cip->haploid_mask, chr_idx)) {
        return cip->chr_fo_vidx_start[chr_fo_idx];
      }
    }
  }
  return 0x7fffffff;
}

uint32_t AlleleDosagesAreNeeded(Command1Flags command_flags1, MiscFlags misc_flags, uint32_t afreq_needed, uint64_t min_allele_ddosage, uint64_t max_allele_ddosage, uint32_t* regular_freqcounts_neededp) {
  if (!(misc_flags & kfMiscNonfounders)) {
    return 0;
  }
  if ((command_flags1 & kfCommand1AlleleFreq) ||
      (misc_flags & kfMiscMajRef) ||
      min_allele_ddosage ||
      (max_allele_ddosage != (~0LLU))) {
    // Keep this in sync with --error-on-freq-calc.
    *regular_freqcounts_neededp = 1;
    return 1;
  }
  return afreq_needed;
}

uint32_t FounderAlleleDosagesAreNeeded(Command1Flags command_flags1, MiscFlags misc_flags, uint32_t afreq_needed, uint64_t min_allele_ddosage, uint64_t max_allele_ddosage, uint32_t* regular_freqcounts_neededp) {
  if (misc_flags & kfMiscNonfounders) {
    return 0;
  }
  if ((command_flags1 & kfCommand1AlleleFreq) ||
      (misc_flags & kfMiscMajRef) ||
      min_allele_ddosage ||
      (max_allele_ddosage != (~0LLU))) {
    // Keep this in sync with --error-on-freq-calc.
    *regular_freqcounts_neededp = 1;
    return 1;
  }
  return afreq_needed;
}
uint32_t VariantMissingHcCtsAreNeeded(MiscFlags misc_flags, double geno_thresh) {
  return (geno_thresh != 1.0) && (!(misc_flags & kfMiscGenoDosage));
}

uint32_t VariantHethapCtsAreNeeded(MiscFlags misc_flags, double geno_thresh, uint32_t first_hap_uidx) {
  // bugfix (20 Oct 2024): --geno hh-missing condition was backwards
  return (first_hap_uidx != 0x7fffffff) &&
    (geno_thresh != 1.0) && (misc_flags & kfMiscGenoHhMissing);
}

uint32_t VariantMissingDosageCtsAreNeeded(MiscFlags misc_flags, double geno_thresh) {
  return (geno_thresh != 1.0) && (misc_flags & kfMiscGenoDosage);
}

// can simplify --geno-counts all-biallelic case, but let's first make sure the
// general case works for multiallelic variants
uint32_t RawGenoCtsAreNeeded(Command1Flags command_flags1, MiscFlags misc_flags, double hwe_ln_thresh) {
  // Keep this in sync with --error-on-freq-calc.
  return (command_flags1 & kfCommand1GenoCounts) ||
    ((misc_flags & kfMiscNonfounders) && ((command_flags1 & kfCommand1Hardy) || (hwe_ln_thresh != -DBL_MAX)));
}

uint32_t FounderRawGenoCtsAreNeeded(Command1Flags command_flags1, MiscFlags misc_flags, double hwe_ln_thresh) {
  // Keep this in sync with --error-on-freq-calc.
  return (!(misc_flags & kfMiscNonfounders)) && ((command_flags1 & kfCommand1Hardy) || (hwe_ln_thresh != -DBL_MAX));
}
uint32_t AlleleUniquenessCheckIsNeeded(const char* read_freq_fname, Command1Flags command_flags1) {
  return (command_flags1 & kfCommand1Score) || read_freq_fname;
}
// The interval-file modifiers --extract/--exclude used to accept.
//
// adjutant has no positional variant filtering: --chr, --from-bp/--to-bp,
// --from/--to and --snp/--snps are all gone, and the BED-interval code went
// with them.  What survived here was the parsing, so 'bed1' was accepted,
// echoed into "Options in effect", and then filtered nothing -- an --extract
// restricting a score to one region silently returned the genome-wide score.
// A filter that passes everything is worse than one that refuses, because the
// run looks like it worked.  Refuse.
uint32_t IsIntervalBedModifier(const char* modif, uint32_t slen) {
  return strequal_k(modif, "bed0", slen) || strequal_k(modif, "ibed0", slen) ||
         strequal_k(modif, "bed1", slen) || strequal_k(modif, "ibed1", slen) ||
         strequal_k(modif, "range", slen);
}

// Chromosome 6.  --chr-set/--autosome-num are not parsed, so the chromosome set
// is always the human default and autosomes are coded by their own number.
static const uint32_t kMhcChrIdx = 6;

// Drop chromosome 6 variants inside the selected MHC interval.
//
// A linear scan, not a binary search over sorted positions.  --from-bp searched
// and so had to reject an unsorted .pvar; walking one chromosome against a fixed
// interval costs nothing worth reclaiming, and in exchange this imposes no
// ordering requirement on the input at all.
//
// Composes with --extract/--exclude by construction: it only ever clears bits,
// so whatever ran before it still holds.
void ApplyMhcExclusion(MhcRegion region, const ChrInfo* cip, const uint32_t* variant_bps, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  if ((region == kMhcNone) || (!(*variant_ct_ptr))) {
    return;
  }
  uint32_t removed_ct = 0;
  for (uint32_t chr_fo_idx = 0; chr_fo_idx != cip->chr_ct; ++chr_fo_idx) {
    if (cip->chr_file_order[chr_fo_idx] != kMhcChrIdx) {
      continue;
    }
    const uint32_t vidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
    for (uint32_t vidx = cip->chr_fo_vidx_start[chr_fo_idx]; vidx != vidx_end; ++vidx) {
      if (IsSet(variant_include, vidx) && MhcExcludesBp(region, variant_bps[vidx])) {
        ClearBit(vidx, variant_include);
        ++removed_ct;
      }
    }
  }
  *variant_ct_ptr -= removed_ct;
  // Print the interval, not just the count, so a log says which coordinates the
  // run actually used rather than leaving the reader to trust the build flag.
  logprintf("%s: chr%u:%u-%u, %u variant%s excluded.\n", MhcFlagName(region), kMhcChrIdx, MhcStartBp(region), MhcEndBp(region), removed_ct, (removed_ct == 1)? "" : "s");
}

void UpdateSampleSubsets(const uintptr_t* sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t* founder_info, uint32_t* founder_ct_ptr, uintptr_t* sex_nm, uintptr_t* sex_male, uint32_t* male_ct_ptr, uint32_t* nosex_ct_ptr) {
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  BitvecAnd(sample_include, raw_sample_ctl, founder_info);
  *founder_ct_ptr = PopcountWords(founder_info, raw_sample_ctl);
  BitvecAnd(sample_include, raw_sample_ctl, sex_male);
  *male_ct_ptr = PopcountWords(sex_male, raw_sample_ctl);
  BitvecAnd(sample_include, raw_sample_ctl, sex_nm);
  *nosex_ct_ptr = sample_ct - PopcountWords(sex_nm, raw_sample_ctl);
}

// command_flags2 will probably be needed before we're done
static_assert(kPglMaxAlleleCt == 255, "Plink2Core() --maj-ref needs to be updated.");
PglErr Plink2Core(const Plink2Cmdline* pcp, MakePlink2Flags make_plink2_flags, char* pgenname, char* psamname, char* pvarname, char* outname, char* outname_end, ChrInfo* cip, sfmt_t* sfmtp) {
  PhenoCol* pheno_cols = nullptr;
  PhenoCol* covar_cols = nullptr;
  PhenoCol* loop_cats_pheno_col = nullptr;
  char* pheno_names = nullptr;
  char* covar_names = nullptr;
  uint32_t pheno_ct = 0;
  uint32_t covar_ct = 0;
  PglErr reterr = kPglRetSuccess;
  PgenFileInfo pgfi;
  PgenReader simple_pgr;
  PreinitPgfi(&pgfi);
  PreinitPgr(&simple_pgr);
  PgenExtensionLl ext_slot; // shouldn't have shorter lifetime than pgfi
  ext_slot.contents = nullptr;
  {
    uint32_t pvar_renamed = 0;
    if (make_plink2_flags & (kfMakeBed | kfMakePgen)) {
      uint32_t fname_slen;
#ifdef _WIN32
      fname_slen = GetFullPathName(pgenname, kPglFnamesize, g_textbuf, nullptr);
      if (unlikely((!fname_slen) || (fname_slen > kPglFnamesize)))
#else
      if (unlikely(!realpath(pgenname, g_textbuf)))
#endif
      {
        logerrprintfww(kErrprintfFopen, pgenname, strerror(errno));
        goto Plink2Core_ret_OPEN_FAIL;
      }
      uint32_t pgen_rename = 0;
      if (make_plink2_flags & kfMakePgen) {
        snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
        pgen_rename = RealpathIdentical(outname, g_textbuf, &(g_textbuf[kPglFnamesize + 64]));
      }
      if (pgen_rename) {
        logerrprintf("Warning: --make-%s input and output filenames match.  Appending '~' to input\nfilenames.\n", (make_plink2_flags & kfMakeBed)? "bed" : ((make_plink2_flags & kfMakePvar)? "pgen" : "bpgen"));
        fname_slen = strlen(pgenname);
        memcpy(g_textbuf, pgenname, fname_slen + 1);
        snprintf(&(pgenname[fname_slen]), 2, "~");
        if (unlikely(rename(g_textbuf, pgenname))) {
          logerrprintf("Error: Failed to append '~' to input .bed/.pgen filename: %s.\n", strerror(errno));
          goto Plink2Core_ret_OPEN_FAIL;
        }
        // if/when --make-pgen can generate .pgen.pgi files, change pginame to
        // a Plink2Core() parameter, and conditionally rename the file here
        fname_slen = strlen(pvarname);
        memcpy(g_textbuf, pvarname, fname_slen + 1);
        snprintf(&(pvarname[fname_slen]), 2, "~");
        if (unlikely(rename(g_textbuf, pvarname))) {
          logerrprintf("Error: Failed to append '~' to input .bim/.pvar filename: %s.\n", strerror(errno));
          goto Plink2Core_ret_OPEN_FAIL;
        }
        pvar_renamed = 1;
        fname_slen = strlen(psamname);
        memcpy(g_textbuf, psamname, fname_slen + 1);
        snprintf(&(psamname[fname_slen]), 2, "~");
        if (unlikely(rename(g_textbuf, psamname))) {
          logerrprintf("Error: Failed to append '~' to input .fam/.psam filename: %s.\n", strerror(errno));
          goto Plink2Core_ret_OPEN_FAIL;
        }
      }
    }
    uint32_t raw_sample_ct = 0;
    uintptr_t* sample_include = nullptr;
    uintptr_t* sex_nm = nullptr;
    uintptr_t* sex_male = nullptr;
    uintptr_t* founder_info = nullptr;
    uintptr_t max_pheno_name_blen = 0;
    uint32_t raw_sample_ctl = 0;
    uint32_t sample_ct = 0;

    uint32_t update_sample_ids_empty = 0;
    uint32_t update_parental_ids_empty = 0;

    // There's a tradeoff between using structs like this and passing each of
    // the components separately: large structs make for shorter parameter
    // lists and hence slightly prettier code, but going too far in that
    // direction (e.g. having a single huge struct that gets passed to all
    // commands) makes it far harder to see what the actual data dependencies
    // are.  Right now the codebase errs, to an almost comical degree, in the
    // direction of making data dependencies clear and minimizing passing of
    // unneeded information.  I'll scale this back, but I'll try to avoid doing
    // it in ways that obscure data flows.
    // For instance, sample IDs and parental IDs could be stored in the same
    // struct, but if we pass parental IDs separately it's easier to see which
    // functions actually need them.
    // In lower-level functions, a reasonable rule of thumb is that aggregation
    // of parameters that normally appear together is worthwhile when there are
    // more than 6 parameters, since the x86-64 C calling convention uses
    // registers for the first 6 parameters.
    PedigreeIdInfo pii;
    InitPedigreeIdInfo(pcp->misc_flags, &pii);
    if (psamname[0]) {
      if (pcp->update_sample_ids_fname) {
        if (reterr) {
          if (unlikely(reterr != kPglRetEof)) {
            goto Plink2Core_ret_1;
          }
          reterr = kPglRetSuccess;
          update_sample_ids_empty = 1;
        }
      }
      if (pcp->update_parental_ids_fname) {
        if (reterr) {
          if (unlikely(reterr != kPglRetEof)) {
            goto Plink2Core_ret_1;
          }
          reterr = kPglRetSuccess;
          update_parental_ids_empty = 1;
        }
      }
      // update (26 Nov 2017): change --no-pheno to also apply to .psam file
      const uint32_t ignore_psam_phenos = (!(pcp->fam_cols & kfFamCol6)) || (pcp->pheno_fname && pcp->pheno_range_list.name_ct);
      reterr = LoadPsam(psamname, pcp->pheno_fname? nullptr : &(pcp->pheno_range_list), pcp->missing_catname, pcp->fam_cols, ignore_psam_phenos? 0 : 0x7fffffff, pcp->missing_pheno, (pcp->misc_flags / kfMiscAffection01) & 1, (pcp->misc_flags / kfMiscNoCategorical) & 1, (pcp->misc_flags / kfMiscNeg9PhenoReallyMissing) & 1, pcp->max_thread_ct, &pii, &sample_include, &founder_info, &sex_nm, &sex_male, &pheno_cols, &pheno_names, &raw_sample_ct, &pheno_ct, &max_pheno_name_blen);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }

      // todo: add option to discard loaded SIDs
      raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
      sample_ct = PopcountWords(sample_include, raw_sample_ctl);
      const uint32_t known_sex_ct = PopcountWords(sex_nm, raw_sample_ctl);
      const uint32_t male_ct = PopcountWords(sex_male, raw_sample_ctl);
      const uint32_t female_ct = known_sex_ct - male_ct;
      const uint32_t founder_ct = PopcountWords(founder_info, raw_sample_ctl);
      if (known_sex_ct == sample_ct) {
        logprintfww("%u sample%s (%u female%s, %u male%s; %u founder%s) loaded from %s.\n", sample_ct, (sample_ct == 1)? "" : "s", female_ct, (female_ct == 1)? "" : "s", male_ct, (male_ct == 1)? "" : "s", founder_ct, (founder_ct == 1)? "" : "s", psamname);
      } else {
        const uint32_t unknown_sex_ct = sample_ct - known_sex_ct;
        logprintfww("%u sample%s (%u female%s, %u male%s, %u ambiguous; %u founder%s) loaded from %s.\n", sample_ct, (sample_ct == 1)? "" : "s", female_ct, (female_ct == 1)? "" : "s", male_ct, (male_ct == 1)? "" : "s", unknown_sex_ct, founder_ct, (founder_ct == 1)? "" : "s", psamname);
      }
    }

    uint32_t max_variant_id_slen = 1;
    uint32_t info_reload_slen = 0;  // --make-pgen/--export/--rm-dup are gone
    uintptr_t* allele_idx_offsets = nullptr;
    uint32_t raw_variant_ct = 0;
    uint32_t variant_ct = 0;
    char* xheader = nullptr;
    uintptr_t xheader_blen = 0;
    uintptr_t* variant_include = nullptr;
    uint32_t* variant_bps = nullptr;
    char** variant_ids_mutable = nullptr;

    // This is actually a hybrid of const char** and char**: length >1
    // pointed-to strings can be modified in-place by --update-alleles, but we
    // must be able to assign const pointers to length-1 strings to it.  So we
    // set the official type to char**, but there's more const_cast-ing than
    // usual.
    char** allele_storage_mutable = nullptr;

    uintptr_t* pvar_qual_present = nullptr;
    float* pvar_quals = nullptr;
    uintptr_t* pvar_filter_present = nullptr;
    uintptr_t* pvar_filter_npass = nullptr;
    const char* const* pvar_filter_storage = nullptr;
    uintptr_t* nonref_flags = nullptr;
    InfoFlags info_flags = kfInfo0;
    uint32_t max_allele_ct = 2;
    uint32_t max_allele_slen = 0;
    uint32_t max_filter_slen = 0;
    UnsortedVar vpos_sortstatus = kfUnsortedVar0;
    double* variant_cms = nullptr;
    ChrIdx* chr_idxs = nullptr;  // split-chromosome case only
    if (pvarname[0]) {
      char** pvar_filter_storage_mutable = nullptr;

      // LoadPvar() uses pvar_psam_flags to determine what's needed for .pvar
      // export.  These booleans are just for tracking requirements beyond
      // that.
      const uint32_t xheader_needed = 0;  // only --export vcf/bcf needed the .pvar header

      const uint32_t qualfilter_needed = xheader_needed;
      uint32_t neg_bp_seen = 0;
      reterr = LoadPvar(pvarname, pcp->var_filter_exceptions_flattened, pcp->varid_template_str, pcp->varid_multi_template_str, pcp->varid_multi_nonsnp_template_str, pcp->missing_varid_match, pcp->require_info_flattened, pcp->require_no_info_flattened, &(pcp->extract_if_info_expr), &(pcp->exclude_if_info_expr), pcp->misc_flags, pcp->pvar_psam_flags, pcp->load_filter_log_flags, xheader_needed, qualfilter_needed, pcp->var_min_qual, pcp->splitpar_bound1, pcp->splitpar_bound2, pcp->new_variant_id_max_allele_slen, (pcp->filter_flags / kfFilterSnpsOnlyJustAcgt) & 1, /* split_chr_ok = */ 1, pcp->filter_min_allele_ct, pcp->filter_max_allele_ct, pcp->input_missing_geno_char, pcp->max_thread_ct, cip, &max_variant_id_slen, &info_reload_slen, &vpos_sortstatus, &xheader, &variant_include, &variant_bps, &variant_ids_mutable, &allele_idx_offsets, K_CAST(const char***, &allele_storage_mutable), &pvar_qual_present, &pvar_quals, &pvar_filter_present, &pvar_filter_npass, &pvar_filter_storage_mutable, &nonref_flags, &variant_cms, &chr_idxs, &raw_variant_ct, &variant_ct, &neg_bp_seen, &max_allele_ct, &max_allele_slen, &xheader_blen, &info_flags, &max_filter_slen);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
      LoadFilterLogFlags load_filter_log_flags = pcp->load_filter_log_flags;
      if (load_filter_log_flags & kfLoadFilterLogImportMergeAlreadyApplied) {
        load_filter_log_flags &= ~kfLoadFilterLogImportMergeMask;
      }
      if (unlikely(!variant_ct)) {
        // conditionally permit this?
        if (!raw_variant_ct) {
          logerrprintfww("Error: No variants in %s.\n", pvarname);
          goto Plink2Core_ret_MALFORMED_INPUT;
        }
        char* write_iter = strcpya_k(g_logbuf, "Error: All ");
        write_iter = u32toa(raw_variant_ct, write_iter);
        write_iter = strcpya_k(write_iter, " variant");
        if (raw_variant_ct != 1) {
          *write_iter++ = 's';
        }
        write_iter = strcpya_k(write_iter, " in ");
        write_iter = strcpya(write_iter, pvarname);
        write_iter = strcpya_k(write_iter, " excluded by ");
        if (load_filter_log_flags) {
          AppendLoadFilterFlagnames(load_filter_log_flags, &write_iter);
        }
        if (neg_bp_seen) {
          if (load_filter_log_flags) {
            write_iter = strcpya_k(write_iter, " and/or ");
          }
          write_iter = strcpya_k(write_iter, "negative bp coordinates");
        }
        strcpy_k(write_iter, ".\n");
        goto Plink2Core_ret_INCONSISTENT_INPUT_WW;
      }
      pvar_filter_storage = TO_CONSTCPCONSTP(pvar_filter_storage_mutable);
      {
        char* write_iter = u32toa(raw_variant_ct, g_logbuf);
        write_iter = strcpya_k(write_iter, " variant");
        if (raw_variant_ct != 1) {
          *write_iter++ = 's';
        }
        if (raw_variant_ct != variant_ct) {
          write_iter = strcpya_k(write_iter, " in ");
          write_iter = strcpya(write_iter, pvarname);
          write_iter = strcpya_k(write_iter, "; ");
          write_iter = u32toa(raw_variant_ct - variant_ct, write_iter);
          write_iter = strcpya_k(write_iter, " excluded by ");
          if (load_filter_log_flags) {
            AppendLoadFilterFlagnames(load_filter_log_flags, &write_iter);
          }
          if (neg_bp_seen) {
            if (load_filter_log_flags) {
              write_iter = strcpya_k(write_iter, " and/or ");
            }
            write_iter = strcpya_k(write_iter, "negative bp coordinates");
          }
          write_iter = strcpya_k(write_iter, ", ");
          write_iter = u32toa(variant_ct, write_iter);
          write_iter = strcpya_k(write_iter, " remaining");
        } else {
          write_iter = strcpya_k(write_iter, " loaded from ");
          write_iter = strcpya(write_iter, pvarname);
          if (load_filter_log_flags) {
            write_iter = strcpya_k(write_iter, " (");
            AppendLoadFilterFlagnames(load_filter_log_flags, &write_iter);
            write_iter = strcpya_k(write_iter, " had no effect)");
          }
        }
        strcpy_k(write_iter, ".\n");
        WordWrapB(0);
        logputsb();
      }
      if (info_reload_slen && (make_plink2_flags & (kfMakeBim | kfMakePvar)) && (!pvar_renamed)) {
        // need to be careful with .pvar in this case
        uint32_t fname_slen;
#ifdef _WIN32
        fname_slen = GetFullPathName(pvarname, kPglFnamesize, g_textbuf, nullptr);
        if (unlikely((!fname_slen) || (fname_slen > kPglFnamesize)))
#else
        if (unlikely(!realpath(pvarname, g_textbuf)))
#endif
        {
          logerrprintfww(kErrprintfFopen, pvarname, strerror(errno));
          goto Plink2Core_ret_OPEN_FAIL;
        }
        if (make_plink2_flags & kfMakeBim) {
          OutnameZstSet(".bim", make_plink2_flags & kfMakeBimZs, outname_end);
          pvar_renamed = RealpathIdentical(outname, g_textbuf, &(g_textbuf[kPglFnamesize + 64]));
          if (pvar_renamed) {
            logerrputs("Warning: .bim input and output filenames match.  Appending '~' to input\nfilename.\n");
            fname_slen = strlen(pvarname);
            memcpy(g_textbuf, pvarname, fname_slen + 1);
            snprintf(&(pvarname[fname_slen]), 2, "~");
            if (unlikely(rename(g_textbuf, pvarname))) {
              logerrprintf("Error: Failed to append '~' to input .bim filename: %s.\n", strerror(errno));
              goto Plink2Core_ret_OPEN_FAIL;
            }
          }
        }
        if ((!pvar_renamed) && (make_plink2_flags & kfMakePvar)) {
          OutnameZstSet(".pvar", pcp->pvar_psam_flags & kfPvarZs, outname_end);
          // pvar_renamed = RealpathIdentical();
          if (RealpathIdentical(outname, g_textbuf, &(g_textbuf[kPglFnamesize + 64]))) {
            logerrputs("Warning: .pvar input and output filenames match.  Appending '~' to input\nfilename.\n");
            fname_slen = strlen(pvarname);
            memcpy(g_textbuf, pvarname, fname_slen + 1);
            snprintf(&(pvarname[fname_slen]), 2, "~");
            if (unlikely(rename(g_textbuf, pvarname))) {
              logerrprintf("Error: Failed to append '~' to input .pvar filename: %s.\n", strerror(errno));
              goto Plink2Core_ret_OPEN_FAIL;
            }
          }
        }
      }
      if ((pcp->dependency_flags & (kfFilterAllReq | kfFilterNonrefFlagsNeeded)) == kfFilterNonrefFlagsNeeded) {
        // We need INFO/PR, but we don't really care about anything else
        // potentially in the .pgen.
        // So if INFO was present in the variant file, we don't need to look at
        // the .pgen at all (even if there was no PR key; we infer that all
        // reference alleles are accurate in that case).
        // Conversely, if no .pgen filename was provided, and the variant file
        // doesn't contain an INFO field, print a warning.
        if (info_flags & kfInfoAll) {
          pgenname[0] = '\0';
        } else if (!pgenname[0]) {
          logerrputs("Warning: Variant file does not distinguish between provisional and trusted REF\nalleles, and no .pgen input file was provided.\n");
          if (info_flags & kfInfoPrNonrefDefault) {
            logerrputs("Assuming all REF alleles are provisional.\n");
          } else {
            logerrputs("Assuming all REF alleles are trusted.\n");
          }
        }
      }
    }

    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uintptr_t pgr_alloc_cacheline_ct = 0;
    if (pgenname[0]) {
      PgenHeaderCtrl header_ctrl;
      uintptr_t cur_alloc_cacheline_ct;
      while (1) {
        reterr = PgfiInitPhase1(pgenname, pcp->pginame, raw_variant_ct, raw_sample_ct, &header_ctrl, &pgfi, &cur_alloc_cacheline_ct, g_logbuf);
        if (!reterr) {
          break;
        }
        // detect and autoconvert plink 1 sample-major files, instead of
        // failing (don't bother supporting plink 0.99 files any more)
        if (unlikely(reterr != kPglRetSampleMajorBed)) {
          WordWrapB(0);
          logerrputsb();
          goto Plink2Core_ret_1;
        }
        char* pgenname_end = memcpya(pgenname, outname, outname_end - outname);
        pgenname_end = strcpya_k(pgenname_end, ".pgen");
        const uint32_t no_vmaj_ext = (pcp->command_flags1 & kfCommand1MakePlink2) && (!pcp->filter_flags) && ((make_plink2_flags & (kfMakePgen | kfMakePlink2MJoin | (kfMakePgenFormatBase * 3) | kfMakePgenWriterVer)) == kfMakePgen);
        if (no_vmaj_ext) {
          *pgenname_end = '\0';
          make_plink2_flags &= ~kfMakePgen;
          // no --make-just-pgen command, so we'll never entirely skip the
          // MakePlink2 operation
        } else {
          snprintf(pgenname_end, kMaxOutfnameExtBlen - 5, ".vmaj");
        }
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
        fclose(pgfi.shared_ff);
        pgfi.shared_ff = nullptr;
      }
      pgfi.allele_idx_offsets = allele_idx_offsets;
      pgfi.max_allele_ct = max_allele_ct;
      unsigned char* pgfi_alloc;
      if (unlikely(bigstack_alloc_uc(cur_alloc_cacheline_ct * kCacheline, &pgfi_alloc))) {
        goto Plink2Core_ret_NOMEM;
      }
      const uint32_t nonref_flags_already_loaded = (nonref_flags != nullptr);
      if (!nonref_flags_already_loaded) {
        const uint32_t nonref_flags_status_shifted = header_ctrl & 192;
        if (nonref_flags_status_shifted == 192) {
          if (unlikely(bigstack_alloc_w(raw_variant_ctl, &nonref_flags))) {
            goto Plink2Core_ret_NOMEM;
          }
        } else if ((pgfi.const_vrtype != kPglVrtypePlink1) && (!nonref_flags_status_shifted)) {
          // This should never happen with a chain of plink2-generated .pgen
          // files, but it's permitted by the specification and makes future
          // merge unsafe.
          logerrputs("Warning: .pgen file indicates that provisional-REF information is stored in the\ncompanion .pvar, but that .pvar lacks either an INFO/PR header or an INFO\ncolumn.  Assuming all REF alleles are correct.\n");
        }
      }
      pgfi.nonref_flags = nonref_flags;
      PgenExtensionLl* header_exts = nullptr;
      uint32_t max_vrec_width;
      // only practical effect of setting use_blockload to zero here is that
      // pgr_alloc_cacheline_ct is overestimated by
      // DivUp(max_vrec_width, kCacheline).
      reterr = PgfiInitPhase2Ex(header_ctrl, 1, nonref_flags_already_loaded, 1, 0, raw_variant_ct, &max_vrec_width, &pgfi, pgfi_alloc, header_exts, nullptr, &pgr_alloc_cacheline_ct, g_logbuf);
      if (unlikely(reterr)) {
        WordWrapB(0);
        logerrputsb();
        goto Plink2Core_ret_1;
      }
      if (unlikely((!allele_idx_offsets) && (pgfi.gflags & kfPgenGlobalMultiallelicHardcallFound))) {
        logerrputs("Error: .pgen file contains multiallelic variants, while .pvar does not.\n");
        goto Plink2Core_ret_INCONSISTENT_INPUT;
      }
      if (pcp->misc_flags & kfMiscRealRefAlleles) {
        if (unlikely(nonref_flags && (!AllBitsAreOne(nonref_flags, raw_variant_ct)))) {
          // To reduce the ease of foot-shooting, we don't allow this to
          // clobber arbitrary nonref_flags arrays.
          logerrputs("Error: --real-ref-alleles must be used on a plink1 or similar fileset.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }

        // wasteful if nonref_flags was allocated, but nonref_flags isn't that
        // large, and --real-ref-alleles + --make-pgen can be run separately
        // from anything truly memory-limited
        nonref_flags = nullptr;
        pgfi.nonref_flags = nullptr;

        pgfi.gflags &= ~kfPgenGlobalAllNonref;
      }
      if (header_exts && (header_exts->size != ~0LLU)) {
#ifndef __LP64__
        if (unlikely(header_exts->size >= (1LLU << 31))) {
          goto Plink2Core_ret_NOMEM;
        }
#endif
        if (unlikely(pgl_malloc(header_exts->size, &(header_exts->contents)))) {
          goto Plink2Core_ret_NOMEM;
        }
        reterr = PgfiInitLoadExts(header_ctrl, &pgfi, header_exts, nullptr, g_logbuf);
        if (unlikely(reterr)) {
          WordWrapB(0);
          logerrputsb();
          goto Plink2Core_ret_1;
        }
      }
      if (SingleVariantLoaderIsNeeded(pcp->command_flags1, pcp->hwe_ln_thresh)) {
        // ugly kludge, probably want to add pgenlib_internal support for this
        // hybrid use pattern
        FILE* shared_ff_copy = pgfi.shared_ff;
        pgfi.shared_ff = nullptr;
        unsigned char* simple_pgr_alloc;
        if (unlikely(bigstack_alloc_uc((pgr_alloc_cacheline_ct + DivUp(max_vrec_width, kCacheline)) * kCacheline, &simple_pgr_alloc))) {
          goto Plink2Core_ret_NOMEM;
        }
        reterr = PgrInit(pgenname, max_vrec_width, &pgfi, &simple_pgr, simple_pgr_alloc);
        if (unlikely(reterr)) {
          if (reterr == kPglRetOpenFail) {
            logerrprintfww(kErrprintfFopen, pgenname, strerror(errno));
          } else {
            assert(reterr == kPglRetReadFail);
            logerrprintfww(kErrprintfFread, pgenname, rstrerror(errno));
          }
          goto Plink2Core_ret_1;
        }
        pgfi.shared_ff = shared_ff_copy;
        }
      // any functions using blockload must perform its own PgrInit(), etc.
      } else {
      // bugfix (10-11 Feb 2018): these variables may still be accessed
      pgfi.gflags = S_CAST(PgenGlobalFlags, ((info_flags / kfInfoPrNonrefDefault) & 1) * kfPgenGlobalAllNonref);
      pgfi.nonref_flags = nonref_flags;
      // dirty hack (3 Oct 2023): this is needed for --sort-vars
      // --make-just-pvar to work
      GET_PRIVATE(simple_pgr, m).fi.gflags = pgfi.gflags;
      GET_PRIVATE(simple_pgr, m).fi.nonref_flags = nonref_flags;
    }
    if (pcp->misc_flags & kfMiscMakeFoundersFirst) {
      reterr = MakeFounders(sample_include, raw_sample_ct, sample_ct, (pcp->misc_flags / kfMiscMakeFoundersRequire2Missing) & 1, &pii, founder_info);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
    }
    if (pcp->pheno_fname) {
      reterr = LoadPhenos(pcp->pheno_fname, &(pcp->pheno_range_list), sample_include, &pii.sii, pcp->missing_catname, raw_sample_ct, sample_ct, pcp->missing_pheno, (pcp->misc_flags / kfMiscAffection01) & 1, (pcp->misc_flags / kfMiscNoCategorical) & 1, (pcp->misc_flags / kfMiscPhenoIidOnly) & 1, (pcp->misc_flags / kfMiscPhenoColNums) & 1, (pcp->misc_flags / kfMiscNeg9PhenoReallyMissing) & 1, pcp->max_thread_ct, &pheno_cols, &pheno_names, &pheno_ct, &max_pheno_name_blen);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
    }
    if (pcp->not_pheno_flattened) {
      reterr = IgnorePhenosOrCovars(pcp->not_pheno_flattened, 0, &pheno_cols, &pheno_names, &pheno_ct, &max_pheno_name_blen);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
    }

    // move processing of PLINK 1.x cluster-loading/filtering flags here, since
    // they're now under the categorical-phenotype umbrella
    if ((pcp->misc_flags & kfMiscCatPhenoFamily) || pcp->within_fname) {
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
    }

    if (!pheno_ct) {
      logputs("Note: No phenotype data present.\n");
    } else {
      if (pheno_ct == 1) {
        const uint32_t obs_ct = PopcountWords(pheno_cols[0].nonmiss, raw_sample_ctl);
        if (pheno_cols[0].type_code == kPhenoDtypeCc) {
          const uint32_t case_ct = PopcountWords(pheno_cols[0].data.cc, raw_sample_ctl);
          const uint32_t ctrl_ct = obs_ct - case_ct;
          logprintf("1 binary phenotype loaded (%u case%s, %u control%s).\n", case_ct, (case_ct == 1)? "" : "s", ctrl_ct, (ctrl_ct == 1)? "" : "s");
        } else if (pheno_cols[0].type_code == kPhenoDtypeQt) {
          logprintf("1 quantitative phenotype loaded (%u value%s).\n", obs_ct, (obs_ct == 1)? "" : "s");
        } else {
          logprintf("1 categorical phenotype loaded (%u value%s).\n", obs_ct, (obs_ct == 1)? "" : "s");
        }
      } else {
        uint32_t cc_ct = 0;
        uint32_t qt_ct = 0;
        for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
          const PhenoDtype cur_type_code = pheno_cols[pheno_idx].type_code;
          if (pheno_cols[pheno_idx].type_code == kPhenoDtypeCc) {
            ++cc_ct;
          } else if (cur_type_code == kPhenoDtypeQt) {
            ++qt_ct;
          }
        }
        uint32_t cat_ct = pheno_ct - cc_ct - qt_ct;
        // just brute force this for now
        if (!cc_ct) {
          if (!qt_ct) {
            logprintf("%u categorical phenotypes loaded.\n", pheno_ct);
          } else if (!cat_ct) {
            logprintf("%u quantitative phenotypes loaded.\n", pheno_ct);
          } else {
            logprintf("%u phenotypes loaded (%u quantitative, %u categorical).\n", pheno_ct, qt_ct, cat_ct);
          }
        } else if (!qt_ct) {
          if (!cat_ct) {
            logprintf("%u binary phenotypes loaded.\n", pheno_ct);
          } else {
            logprintf("%u phenotypes loaded (%u binary, %u categorical).\n", pheno_ct, cc_ct, cat_ct);
          }
        } else if (!cat_ct) {
          logprintf("%u phenotypes loaded (%u binary, %u quantitative).\n", pheno_ct, cc_ct, qt_ct);
        } else {
          logprintfww("%u phenotypes loaded (%u binary, %u quantitative, %u categorical).\n", pheno_ct, cc_ct, qt_ct, cat_ct);
        }
      }
    }
    if (variant_ct && AlleleUniquenessCheckIsNeeded(pcp->read_freq_fname, pcp->command_flags1)) {
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
    }
    // If something like --snps is combined with a position-based filter which
    // may remove some of the named variants, we need to apply --snps first.
    // Otherwise, it may be very advantageous to apply the position-based
    // filters before constructing the variant ID hash table.  So we split this
    // into two cases.
    // (Actually, hash table construction is now usually fast enough, even with
    // >80m variants, that this optimization isn't a big deal.  But it's
    // already implemented, so let's keep it around until/unless it creates
    // some sort of conflict.)

    // Additional minor optimization: --recover-var-ids doesn't actually need
    // the ID hash table; it's just positioned in this code block anyway
    // because it would be too weird for it to be in a different position than
    // --update-name in the order of operations.
    const uint32_t htable_needed_early = variant_ct && (pcp->varid_from || pcp->varid_to || pcp->varid_snp || pcp->varid_exclude_snp || pcp->snps_range_list.name_ct || pcp->exclude_snps_range_list.name_ct);
    const uint32_t full_variant_id_htable_needed = variant_ct && (htable_needed_early || pcp->update_map_flag || pcp->update_name_flag || pcp->extract_col_cond_info.params);

    // Every --extract/--exclude is by ID now; the interval forms are rejected
    // during parsing.
    const uint32_t extract_exclude_by_id = pcp->extract_fnames || pcp->extract_intersect_fnames || pcp->exclude_fnames;
    if (variant_ct && (full_variant_id_htable_needed || extract_exclude_by_id || pcp->recover_var_ids_fname)) {
      if (vpos_sortstatus & kfUnsortedVarBp) {
        if (unlikely(pcp->varid_from || pcp->varid_to)) {
          logerrputs("Error: --from/--to require a sorted .pvar/.bim.  Retry this command after using\n--make-pgen/--make-bed + --sort-vars to sort your data.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
        if (unlikely(pcp->window_bp != -1)) {
          logerrputs("Error: --window requires a sorted .pvar/.bim.  Retry this command after using\n--make-pgen/--make-bed + --sort-vars to sort your data.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
        if (unlikely(pcp->recover_var_ids_fname)) {
          logerrputs("Error: --recover-var-ids requires a sorted .pvar/.bim.  Retry this command\nafter using --make-pgen/--make-bed + --sort-vars to sort your data.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
      }
      // don't bother with having different allow_dups vs. no allow_dups hash
      // table structures, just check specific IDs for duplication in the
      // no-duplicates-allowed cases
      unsigned char* bigstack_mark = g_bigstack_base;
      uint32_t* variant_id_htable = nullptr;
      uint32_t* htable_dup_base = nullptr;
      uint32_t dup_ct = 0;
      uint32_t variant_id_htable_size = 0;
      if ((!pcp->recover_var_ids_fname) || htable_needed_early) {
        reterr = AllocAndPopulateIdHtableMt(variant_include, TO_CONSTCPCONSTP(variant_ids_mutable), variant_ct, bigstack_left() / 8, pcp->max_thread_ct, &variant_id_htable, &htable_dup_base, &variant_id_htable_size, &dup_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->varid_from || pcp->varid_to) {
        reterr = FromToFlag(TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, htable_dup_base, pcp->varid_from, pcp->varid_to, raw_variant_ct, max_variant_id_slen, variant_id_htable_size, variant_include, cip, &variant_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }

      if (variant_ct) {
        if (pcp->update_map_flag) {
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        } else if (pcp->update_name_flag || pcp->recover_var_ids_fname) {
          if (pcp->update_name_flag) {
          } else {
          }
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
          }

        if (pcp->extract_fnames) {
          reterr = TokenExtractExclude(TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, htable_dup_base, pcp->extract_fnames, "extract", raw_variant_ct, max_variant_id_slen, variant_id_htable_size, kVfilterExtract, pcp->max_thread_ct, variant_include, &variant_ct);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        }
        if (pcp->extract_intersect_fnames) {
          reterr = TokenExtractExclude(TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, htable_dup_base, pcp->extract_intersect_fnames, "extract-intersect", raw_variant_ct, max_variant_id_slen, variant_id_htable_size, kVfilterExtractIntersect, pcp->max_thread_ct, variant_include, &variant_ct);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        }
        if (pcp->exclude_fnames) {
          reterr = TokenExtractExclude(TO_CONSTCPCONSTP(variant_ids_mutable), variant_id_htable, htable_dup_base, pcp->exclude_fnames, "exclude", raw_variant_ct, max_variant_id_slen, variant_id_htable_size, kVfilterExclude, pcp->max_thread_ct, variant_include, &variant_ct);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        }
        }

      BigstackReset(bigstack_mark);

      // todo: --attrib; although it isn't a "standard" file format, it can be
      // more convenient than forcing users to generate full-blown sites-only
      // VCF files, etc.
    }

    // After the ID-based filters, so the count in the log is what the MHC
    // exclusion actually removed rather than what it would have removed from
    // the unfiltered file.
    ApplyMhcExclusion(pcp->mhc_region, cip, variant_bps, variant_include, &variant_ct);

    // variant_ids[] is fixed from this point on.
    const char* const* variant_ids = TO_CONSTCPCONSTP(variant_ids_mutable);

    // SetRefalt1FromFile() can alter pointers-to-missing in
    // allele_storage[] so we can't make this a const char* const*, but nothing
    // alters the pointed-to strings after this point.
    const char** allele_storage = K_CAST(const char**, allele_storage_mutable);

    uint32_t* sample_missing_dosage_cts = nullptr;
    uint32_t* sample_missing_hc_cts = nullptr;
    uint32_t* sample_hethap_cts = nullptr;
    uintptr_t max_covar_name_blen = 0;
    if (psamname[0]) {
      // xid_mode may vary between these operations in a single run, and
      // sample-sort is relatively cheap, so we abandon plink 1.9's "construct
      // sample ID map only once" optimization.
      if (pcp->update_sample_ids_fname) {
        if (update_sample_ids_empty) {
          logputs("--update-ids: 0 samples updated.\n");
        } else {
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        }
      } else {
        if (pcp->update_parental_ids_fname) {
          if (update_parental_ids_empty) {
            logputs("--update-parents: 0 samples updated.\n");
          } else {
            if (unlikely(reterr)) {
              goto Plink2Core_ret_1;
            }
            // bugfix (16 Jan 2021)
            pii.sii.flags |= kfSampleIdParentsPresent;
          }
        }
        // --update-parents goes here
        }
      if ((pii.sii.flags & kfSampleIdFidPresent) && (pii.sii.flags & kfSampleIdNoIdHeaderIidOnly)) {
        for (uint32_t sample_idx = 0; sample_idx != raw_sample_ct; ++sample_idx) {
          const char* cur_sample_id = &(pii.sii.sample_ids[sample_idx * pii.sii.max_sample_id_blen]);
          if (unlikely(!memequal_sk(cur_sample_id, "0\t"))) {
            logerrputs("Error: 'iid-only' modifier can only be used when FIDs are missing or all-0.\n");
            goto Plink2Core_ret_INCONSISTENT_INPUT;
          }
        }
      }
      if (pcp->keepfam_fnames) {
        reterr = KeepOrRemove(pcp->keepfam_fnames, &pii.sii, raw_sample_ct, kfKeepFam, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->keep_fnames) {
        reterr = KeepOrRemove(pcp->keep_fnames, &pii.sii, raw_sample_ct, kfKeep0, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->removefam_fnames) {
        reterr = KeepOrRemove(pcp->removefam_fnames, &pii.sii, raw_sample_ct, kfKeepRemove | kfKeepFam, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->remove_fnames) {
        reterr = KeepOrRemove(pcp->remove_fnames, &pii.sii, raw_sample_ct, kfKeepRemove, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->indv_str) {
        KeepOneId(pcp->indv_str, &pii.sii, raw_sample_ct, (pcp->misc_flags / kfMiscIidSid) & 1, sample_include, &sample_ct);
      }

      // todo: --attrib-indiv

      if (pcp->misc_flags & kfMiscRequirePheno) {
        reterr = RequirePheno(pheno_cols, pheno_names, pcp->require_pheno_flattened, raw_sample_ct, pheno_ct, max_pheno_name_blen, 0, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->filter_flags & (kfFilterExclFemales | kfFilterExclMales | kfFilterExclNosex)) {
        if (pcp->filter_flags & kfFilterExclFemales) {
          for (uint32_t widx = 0; widx != raw_sample_ctl; ++widx) {
            sample_include[widx] &= (~sex_nm[widx]) | sex_male[widx];
          }
        }
        if (pcp->filter_flags & kfFilterExclMales) {
          BitvecInvmask(sex_male, raw_sample_ctl, sample_include);
        }
        if (pcp->filter_flags & kfFilterExclNosex) {
          BitvecAnd(sex_nm, raw_sample_ctl, sample_include);
        }
        const uint32_t old_sample_ct = sample_ct;
        sample_ct = PopcountWords(sample_include, raw_sample_ctl);
        const uint32_t removed_ct = old_sample_ct - sample_ct;
        logprintf("%u sample%s removed due to sex filter(s).\n", removed_ct, (removed_ct == 1)? "" : "s");
      }
      if (pcp->filter_flags & (kfFilterExclFounders | kfFilterExclNonfounders)) {
        const uint32_t keep_founders = (pcp->filter_flags / kfFilterExclNonfounders) & 1;
        if (keep_founders) {
          BitvecAnd(founder_info, raw_sample_ctl, sample_include);
        } else {
          BitvecInvmask(founder_info, raw_sample_ctl, sample_include);
        }
        const uint32_t old_sample_ct = sample_ct;
        sample_ct = PopcountWords(sample_include, raw_sample_ctl);
        const uint32_t removed_ct = old_sample_ct - sample_ct;
        logprintf("--keep-%sfounders: %u sample%s removed.\n", keep_founders? "" : "non", removed_ct, (removed_ct == 1)? "" : "s");
      }

      if (pcp->filter_flags & kfFilterSelectSidRepresentatives) {
        reterr = SelectSidRepresentatives(sex_nm, sex_male, (pcp->select_sid_missingness_mode == kSelectSidMissingnessDosage)? sample_missing_dosage_cts : sample_missing_hc_cts, (pcp->select_sid_missingness_mode == kSelectSidMissingnessHhMissing)? sample_hethap_cts : nullptr, &pii, pcp->select_sid_missingness_mode, pcp->select_sid_tiebreak_mode, (pcp->misc_flags / kfMiscSelectSidParentsOnly) & 1, raw_sample_ct, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      // this results in a small "memory leak" when a regular missingness
      // report is requested, not a big deal
      if (pcp->misc_flags & kfMiscMakeFoundersNotfirst) {
        reterr = MakeFounders(sample_include, raw_sample_ct, sample_ct, (pcp->misc_flags / kfMiscMakeFoundersRequire2Missing) & 1, &pii, founder_info);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      if (pcp->covar_fname || pcp->covar_range_list.name_ct) {
        const char* cur_covar_fname = pcp->covar_fname? pcp->covar_fname : (pcp->pheno_fname? pcp->pheno_fname : psamname);
        reterr = LoadPhenos(cur_covar_fname, &(pcp->covar_range_list), sample_include, &pii.sii, pcp->missing_catname, raw_sample_ct, sample_ct, pcp->missing_pheno, 2, (pcp->misc_flags / kfMiscNoCategorical) & 1, (pcp->misc_flags / kfMiscCovarIidOnly) & 1, (pcp->misc_flags / kfMiscCovarColNums) & 1, (pcp->misc_flags / kfMiscNeg9PhenoReallyMissing) & 1, pcp->max_thread_ct, &covar_cols, &covar_names, &covar_ct, &max_covar_name_blen);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
        logprintf("%u covariate%s loaded from %s.\n", covar_ct, (covar_ct == 1)? "" : "s", cur_covar_fname);

        // do we still want to clear some main phenotype values here if some
        // covariate values are missing?  (don't think there's a point to
        // preserving that behavior, let the regression functions do it to
        // their local phenotype copies on their own.)
      }
      if (pcp->not_covar_flattened) {
        reterr = IgnorePhenosOrCovars(pcp->not_covar_flattened, 1, &covar_cols, &covar_names, &covar_ct, &max_covar_name_blen);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }


      if (pcp->misc_flags & kfMiscRequireCovar) {
        reterr = RequirePheno(covar_cols, covar_names, pcp->require_covar_flattened, raw_sample_ct, covar_ct, max_covar_name_blen, 1, sample_include, &sample_ct);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
      }
      // meow
    }

    // Validate --give-best-pgs requirements before allele-frequency calculation
    // and scoring, so a misconfigured run fails immediately rather than after
    // wasting scoring work.
    if (pcp->score_info.flags & kfScoreGiveBestPgs) {
      // Temporary storage for the usable-phenotype index list; we only care
      // whether validation passes here, not which phenotypes are usable.
      uint32_t* tmp_usable_idxs = S_CAST(uint32_t*, malloc(pheno_ct? pheno_ct * sizeof(uint32_t) : sizeof(uint32_t)));
      if (!tmp_usable_idxs) {
        goto Plink2Core_ret_NOMEM;
      }
      uint32_t tmp_usable_ct = 0;
      reterr = BestPgsValidate(sample_include, sample_ct,
                               pheno_cols, pheno_names, max_pheno_name_blen,
                               pheno_ct, covar_cols, covar_ct, pcp->pgs_folds,
                               tmp_usable_idxs, &tmp_usable_ct);
      free(tmp_usable_idxs);
      if (unlikely(reterr)) {
        goto Plink2Core_ret_1;
      }
    }

    const uint32_t nonfounders = (pcp->misc_flags / kfMiscNonfounders) & 1;
    uint32_t founder_ct = 0;
    uint32_t male_ct = 0;
    uint32_t nosex_ct = 0;
    if (psamname[0]) {
      if (unlikely(!sample_ct)) {
        logerrputs("Error: No samples remaining after main filters.\n");
        goto Plink2Core_ret_DEGENERATE_DATA;
      }
      UpdateSampleSubsets(sample_include, raw_sample_ct, sample_ct, founder_info, &founder_ct, sex_nm, sex_male, &male_ct, &nosex_ct);
      if (pcp->filter_flags & kfFilterPsamReq) {
        const uint32_t female_ct = sample_ct - male_ct - nosex_ct;
        if (!nosex_ct) {
          logprintfww("%u sample%s (%u female%s, %u male%s; %u founder%s) remaining after main filters.\n", sample_ct, (sample_ct == 1)? "" : "s", female_ct, (female_ct == 1)? "" : "s", male_ct, (male_ct == 1)? "" : "s", founder_ct, (founder_ct == 1)? "" : "s");
        } else {
          logprintfww("%u sample%s (%u female%s, %u male%s, %u ambiguous; %u founder%s) remaining after main filters.\n", sample_ct, (sample_ct == 1)? "" : "s", female_ct, (female_ct == 1)? "" : "s", male_ct, (male_ct == 1)? "" : "s", nosex_ct, founder_ct, (founder_ct == 1)? "" : "s");
        }
        if (pheno_ct == 1) {
          const PhenoDtype pheno_type_code = pheno_cols[0].type_code;
          const uint32_t obs_ct = PopcountWordsIntersect(pheno_cols[0].nonmiss, sample_include, raw_sample_ctl);
          if (pheno_type_code == kPhenoDtypeCc) {
            const uint32_t case_ct = PopcountWordsIntersect(pheno_cols[0].data.cc, sample_include, raw_sample_ctl);
            const uint32_t ctrl_ct = obs_ct - case_ct;
            logprintf("%u case%s and %u control%s remaining after main filters.\n", case_ct, (case_ct == 1)? "" : "s", ctrl_ct, (ctrl_ct == 1)? "" : "s");
          } else {
            logprintf("%u %s phenotype value%s remaining after main filters.\n", obs_ct, (pheno_type_code == kPhenoDtypeQt)? "quantitative" : "categorical", (obs_ct == 1)? "" : "s");
          }
        }
      }
    }
    // quantnorm before variance-standardize, since at least that has a minor
    // effect, whereas the other order is pointless
    #ifndef NOLAPACK
    #endif

    char* loop_cats_outname_endp1_backup = &(outname_end[1]);
    uintptr_t* loop_cats_sample_include_backup = nullptr;
    uintptr_t* loop_cats_founder_info_backup = nullptr;
    uintptr_t* loop_cats_sex_nm_backup = nullptr;
    uintptr_t* loop_cats_sex_male_backup = nullptr;
    uintptr_t* loop_cats_variant_include_backup = nullptr;
    uintptr_t* loop_cats_cat_include = nullptr;
    uint32_t loop_cats_sample_ct = 0;
    uint32_t loop_cats_variant_ct = variant_ct;
    uint32_t loop_cats_uidx = 0;
    uint32_t loop_cats_ct = 1;
    if (pcp->loop_cats_phenoname) {
      // 1. check phenotype names, verify it's a categorical pheno; fall back
      //    on checking covariates if no phenotype has same name
      // 2. remove phenotype/covariate, back up sample_include
      // 3. iterate through positive integer category numbers, skip empty
      //    categories
      const char* loop_cats_phenoname = pcp->loop_cats_phenoname;
      const uintptr_t name_blen = 1 + strlen(loop_cats_phenoname);
      if (name_blen <= max_pheno_name_blen) {
        // this boilerplate may belong in its own function
        for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
          if (memequal(loop_cats_phenoname, &(pheno_names[pheno_idx * max_pheno_name_blen]), name_blen)) {
            PhenoCol* cur_pheno_col = &(pheno_cols[pheno_idx]);
            if (unlikely(cur_pheno_col->type_code != kPhenoDtypeCat)) {
              logerrprintfww("Error: '%s' is not a categorical phenotype.\n", loop_cats_phenoname);
              goto Plink2Core_ret_INCONSISTENT_INPUT;
            }
            // bugfix (13 Dec 2019): we need to copy the contents of
            // pheno_cols[pheno_idx] when pheno_idx != pheno_ct - 1.
            // Also note that this makes pheno_ct == 0 and pheno_cols !=
            // nullptr possible simultaneously.
            --pheno_ct;
            if (pheno_idx < pheno_ct) {
              PhenoCol pheno_col_copy = *cur_pheno_col;
              for (uint32_t pheno_idx2 = pheno_idx; pheno_idx2 != pheno_ct; ++pheno_idx2) {
                pheno_cols[pheno_idx2] = pheno_cols[pheno_idx2 + 1];
                memcpy(&(pheno_names[pheno_idx2 * max_pheno_name_blen]), &(pheno_names[(pheno_idx2 + 1) * max_pheno_name_blen]), max_pheno_name_blen);
              }
              pheno_cols[pheno_ct] = pheno_col_copy;
            }
            loop_cats_pheno_col = &(pheno_cols[pheno_ct]);
            break;
          }
        }
      }
      if ((!loop_cats_pheno_col) && (name_blen <= max_covar_name_blen)) {
        for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
          if (memequal(loop_cats_phenoname, &(covar_names[covar_idx * max_covar_name_blen]), name_blen)) {
            PhenoCol* cur_covar_col = &(covar_cols[covar_idx]);
            if (unlikely(cur_covar_col->type_code != kPhenoDtypeCat)) {
              logerrprintfww("Error: '%s' is not a categorical covariate.\n", loop_cats_phenoname);
              goto Plink2Core_ret_INCONSISTENT_INPUT;
            }
            --covar_ct;
            if (covar_idx < covar_ct) {
              PhenoCol covar_col_copy = *cur_covar_col;
              for (uint32_t covar_idx2 = covar_idx; covar_idx2 != covar_ct; ++covar_idx2) {
                covar_cols[covar_idx2] = covar_cols[covar_idx2 + 1];
                memcpy(&(covar_names[covar_idx2 * max_covar_name_blen]), &(covar_names[(covar_idx2 + 1) * max_covar_name_blen]), max_covar_name_blen);
              }
              covar_cols[covar_ct] = covar_col_copy;
            }
            loop_cats_pheno_col = &(covar_cols[covar_ct]);
            break;
          }
        }
      }
      if (unlikely(!loop_cats_pheno_col)) {
        logerrprintfww("Error: --loop-cats phenotype '%s' not loaded.\n", loop_cats_phenoname);
        goto Plink2Core_ret_INCONSISTENT_INPUT;
      }
      if (unlikely(bigstack_alloc_w(raw_sample_ctl, &loop_cats_sample_include_backup) ||
                   bigstack_alloc_w(raw_sample_ctl, &loop_cats_founder_info_backup) ||
                   bigstack_alloc_w(raw_sample_ctl, &loop_cats_sex_nm_backup) ||
                   bigstack_alloc_w(raw_sample_ctl, &loop_cats_sex_male_backup) ||
                   bigstack_alloc_w(1 + (loop_cats_pheno_col->nonnull_category_ct / kBitsPerWord), &loop_cats_cat_include))) {
        goto Plink2Core_ret_NOMEM;
      }
      if (variant_ct != raw_variant_ct) {
        if (unlikely(bigstack_alloc_w(raw_variant_ctl, &loop_cats_variant_include_backup))) {
          goto Plink2Core_ret_NOMEM;
        }
        memcpy(loop_cats_variant_include_backup, variant_include, raw_variant_ctl * sizeof(intptr_t));
      }
      BitvecAndCopy(sample_include, loop_cats_pheno_col->nonmiss, raw_sample_ctl, loop_cats_sample_include_backup);
      loop_cats_sample_ct = PopcountWords(loop_cats_sample_include_backup, raw_sample_ctl);
      loop_cats_ct = IdentifyRemainingCats(sample_include, loop_cats_pheno_col, sample_ct, loop_cats_cat_include);
      // bugfix (13 Dec 2019): --loop-cats does not iterate over the null
      // category.  Don't include it in the count.
      loop_cats_ct -= loop_cats_cat_include[0] & 1;
      if (unlikely(!loop_cats_ct)) {
        logerrputs("Error: All --loop-cats categories are empty.\n");
        goto Plink2Core_ret_INCONSISTENT_INPUT;
      }
      logprintf("--loop-cats: %u categor%s present.\n", loop_cats_ct, (loop_cats_ct == 1)? "y" : "ies");
      BitvecAndCopy(sample_include, founder_info, raw_sample_ctl, loop_cats_founder_info_backup);
      BitvecAndCopy(sample_include, sex_nm, raw_sample_ctl, loop_cats_sex_nm_backup);
      BitvecAndCopy(sample_include, sex_male, raw_sample_ctl, loop_cats_sex_male_backup);
      *outname_end = '.';
    }
    const uintptr_t raw_allele_ct = allele_idx_offsets? allele_idx_offsets[raw_variant_ct] : (2 * raw_variant_ct);
    unsigned char* bigstack_mark_varfilter = g_bigstack_base;
    unsigned char* bigstack_end_mark_varfilter = g_bigstack_end;
    for (uint32_t loop_cats_idx = 0; loop_cats_idx != loop_cats_ct; ++loop_cats_idx) {
      BigstackDoubleReset(bigstack_mark_varfilter, bigstack_end_mark_varfilter);
      if (loop_cats_pheno_col) {
        loop_cats_uidx = AdvTo1Bit(loop_cats_cat_include, loop_cats_uidx + 1);
        const char* catname = loop_cats_pheno_col->category_names[loop_cats_uidx];
        const uint32_t catname_slen = strlen(catname);
        if (unlikely(catname_slen + S_CAST(uintptr_t, loop_cats_outname_endp1_backup - outname) > (kPglFnamesize - kMaxOutfnameExtBlen))) {
          logerrputs("Error: --loop-cats category name too long.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
        outname_end = memcpya(loop_cats_outname_endp1_backup, catname, catname_slen);
        sample_ct = GetCatSamples(loop_cats_sample_include_backup, loop_cats_pheno_col, raw_sample_ctl, loop_cats_sample_ct, loop_cats_uidx, sample_include);
        memcpy(founder_info, loop_cats_founder_info_backup, raw_sample_ctl * sizeof(intptr_t));
        memcpy(sex_nm, loop_cats_sex_nm_backup, raw_sample_ctl * sizeof(intptr_t));
        memcpy(sex_male, loop_cats_sex_male_backup, raw_sample_ctl * sizeof(intptr_t));
        UpdateSampleSubsets(sample_include, raw_sample_ct, sample_ct, founder_info, &founder_ct, sex_nm, sex_male, &male_ct, &nosex_ct);
        variant_ct = loop_cats_variant_ct;
        if (loop_cats_variant_include_backup) {
          memcpy(variant_include, loop_cats_variant_include_backup, raw_variant_ctl * sizeof(intptr_t));
        } else {
          SetAllBits(variant_ct, variant_include);
        }
        logprintf("--loop-cats: Processing category '%s' (%u sample%s).\n", catname, sample_ct, (sample_ct == 1)? "" : "s");
      }

      uintptr_t* allele_presents = nullptr;

      // dosages are currently in 32768ths, hence the extra 'd'
      // same indexes as allele_storage.  We can't omit the last allele (in the
      // way allele_freqs does) because the sum isn't constant (missing
      // values).
      uint64_t* allele_ddosages = nullptr;
      uint64_t* founder_allele_ddosages = nullptr;

      double* allele_freqs = nullptr;
      STD_ARRAY_PTR_DECL(uint32_t, 3, raw_geno_cts) = nullptr;
      STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts) = nullptr;
      unsigned char* bigstack_mark_allele_ddosages = g_bigstack_base;

      if (pgenname[0]) {
        const uint32_t decent_afreqs_needed = DecentAlleleFreqsAreNeeded(pcp->command_flags1, pcp->score_info.flags);
        if (decent_afreqs_needed || IndecentAlleleFreqsAreNeeded(pcp->min_maf, pcp->max_maf)) {
          if (unlikely((!pcp->read_freq_fname) && ((sample_ct < 50) || (founder_ct < 50)) && decent_afreqs_needed)) {
            // --freq, --nonfounders and --bad-freqs no longer exist in this
            // build, so --read-freq is the only remedy left to suggest.
            logerrputs("Error: --score needs allele frequencies to mean-impute missing genotypes, but\nfewer than 50 samples/founders are available to estimate them from.\nSupply frequencies from a larger reference dataset with --read-freq, or add the\n'no-mean-imputation' modifier to --score to drop missing calls instead.\n");
            goto Plink2Core_ret_DEGENERATE_DATA;
          }
          //   allele_freqs[allele_idx_offsets[variant_uidx] - variant_uidx]
          // stores the frequency estimate for the reference allele; if there's
          // more than 1 alt allele, next element stores alt1 freq, etc.  To
          // save memory, we omit the last alt.
          if (unlikely(bigstack_alloc_d(raw_allele_ct - raw_variant_ct, &allele_freqs))) {
            goto Plink2Core_ret_NOMEM;
          }
        }
        if ((sample_ct != founder_ct) && (pcp->min_allele_ddosage || (pcp->max_allele_ddosage != (~0LLU))) && (!nonfounders) && (!(pcp->misc_flags & kfMiscAcFounders))) {
          logerrputs("Error: --mac/--max-mac/\"--freq counts\" specified, but with neither\n--ac-founders nor --nonfounders; and nonfounders are present.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
        uint32_t x_start = 0;
        uint32_t x_len = 0;
        uint32_t hwe_x_probs_needed = 0;
        uint32_t x_code;
        if (XymtExists(cip, kChrOffsetX, &x_code) && (!(vpos_sortstatus & kfUnsortedVarSplitChr))) {
          const uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[x_code];
          x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
          const uint32_t x_end = cip->chr_fo_vidx_start[x_chr_fo_idx + 1];
          x_len = x_end - x_start;
          if (x_len && ((pcp->command_flags1 & kfCommand1Hardy) || (pcp->hwe_ln_thresh != -DBL_MAX)) && (!AllBitsAreZero(variant_include, x_start, x_end))) {
            if (nonfounders) {
              hwe_x_probs_needed = (sample_ct > nosex_ct);
            } else {
              // at least one founder with known gender?
              hwe_x_probs_needed = !IntersectionIsEmpty(founder_info, sex_nm, raw_sample_ctl);
            }
          }
        }
        bigstack_mark_allele_ddosages = g_bigstack_base;
        const uint32_t first_hap_uidx = GetFirstHaploidUidx(cip, vpos_sortstatus);
        if (make_plink2_flags & kfMakePlink2TrimAlts) {
          if (unlikely(bigstack_alloc_w(BitCtToWordCt(raw_allele_ct), &allele_presents))) {
            goto Plink2Core_ret_NOMEM;
          }
        }
        uint32_t regular_freqcounts_needed = (allele_presents != nullptr);
        if (AlleleDosagesAreNeeded(pcp->command_flags1, pcp->misc_flags, (allele_freqs != nullptr), pcp->min_allele_ddosage, pcp->max_allele_ddosage, &regular_freqcounts_needed)) {
          if (unlikely(bigstack_alloc_u64(raw_allele_ct, &allele_ddosages))) {
            goto Plink2Core_ret_NOMEM;
          }
        }
        uint32_t overlapping_allele_ddosages = 0;
        if (FounderAlleleDosagesAreNeeded(pcp->command_flags1, pcp->misc_flags, (allele_freqs != nullptr), pcp->min_allele_ddosage, pcp->max_allele_ddosage, &regular_freqcounts_needed)) {
          if (allele_ddosages && (founder_ct == sample_ct)) {
            overlapping_allele_ddosages = 1;
          } else {
            if (unlikely(bigstack_alloc_u64(raw_allele_ct, &founder_allele_ddosages))) {
              goto Plink2Core_ret_NOMEM;
            }
          }
        }
        double* imp_r2_vals = nullptr;
        const uint32_t is_minimac3_r2 = (pcp->minimac3_r2_max != 0.0);
        if (is_minimac3_r2 || (pcp->mach_r2_max != 0.0)) {
          if (unlikely(bigstack_alloc_d(raw_variant_ct, &imp_r2_vals))) {
            goto Plink2Core_ret_NOMEM;
          }
        }

        unsigned char* bigstack_mark_geno_cts = g_bigstack_base;

        // no longer includes hethaps by default
        uint32_t* variant_missing_hc_cts = nullptr;
        uint32_t* variant_hethap_cts = nullptr;
        if (VariantMissingHcCtsAreNeeded(pcp->misc_flags, pcp->geno_thresh)) {
          if (unlikely(bigstack_alloc_u32(raw_variant_ct, &variant_missing_hc_cts))) {
            goto Plink2Core_ret_NOMEM;
          }
          if (VariantHethapCtsAreNeeded(pcp->misc_flags, pcp->geno_thresh, first_hap_uidx)) {
            // first_hap_uidx offset can save an entire GB...
            if (unlikely(bigstack_alloc_u32(raw_variant_ct - first_hap_uidx, &variant_hethap_cts))) {
              goto Plink2Core_ret_NOMEM;
            }
          }
        }
        uint32_t* variant_missing_dosage_cts = nullptr;
        if (VariantMissingDosageCtsAreNeeded(pcp->misc_flags, pcp->geno_thresh)) {
          if ((!variant_missing_hc_cts) || (pgfi.gflags & kfPgenGlobalDosagePresent)) {
            if (unlikely(bigstack_alloc_u32(raw_variant_ct, &variant_missing_dosage_cts))) {
              goto Plink2Core_ret_NOMEM;
            }
          } else {
            variant_missing_dosage_cts = variant_missing_hc_cts;
          }
        }
        STD_ARRAY_PTR_DECL(uint32_t, 3, x_male_geno_cts) = nullptr;
        STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts) = nullptr;
        STD_ARRAY_PTR_DECL(uint32_t, 3, x_nosex_geno_cts) = nullptr;
        STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts) = nullptr;
        // [0] = homref ct, [1] = het ref-altx total, [2] = nonref diploid
        //   total
        // use unfiltered indexes, since we remove more variants later
        if (RawGenoCtsAreNeeded(pcp->command_flags1, pcp->misc_flags, pcp->hwe_ln_thresh)) {
          if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(uint32_t, 3, raw_variant_ct, &raw_geno_cts))) {
            goto Plink2Core_ret_NOMEM;
          }
          if (x_len) {
            if (male_ct) {
              if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(uint32_t, 3, x_len, &x_male_geno_cts))) {
                goto Plink2Core_ret_NOMEM;
              }
            }
            if (nosex_ct && hwe_x_probs_needed && nonfounders) {
              if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(uint32_t, 3, x_len, &x_nosex_geno_cts))) {
                goto Plink2Core_ret_NOMEM;
              }
            }
          }
        }
        if (FounderRawGenoCtsAreNeeded(pcp->command_flags1, pcp->misc_flags, pcp->hwe_ln_thresh)) {
          if ((founder_ct == sample_ct) && raw_geno_cts) {
            founder_raw_geno_cts = raw_geno_cts;
            founder_x_male_geno_cts = x_male_geno_cts;
          } else {
            if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(uint32_t, 3, raw_variant_ct, &founder_raw_geno_cts))) {
              goto Plink2Core_ret_NOMEM;
            }
            if (x_len && male_ct) {
              if (!IntersectionIsEmpty(founder_info, sex_male, raw_sample_ctl)) {
                if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(uint32_t, 3, x_len, &founder_x_male_geno_cts))) {
                  goto Plink2Core_ret_NOMEM;
                }
              }
            }
          }
          if (nosex_ct && hwe_x_probs_needed && (!nonfounders)) {
            const uint32_t founder_knownsex_ct = PopcountWordsIntersect(founder_info, sex_nm, raw_sample_ctl);
            if (founder_knownsex_ct < founder_ct) {
              if ((founder_ct == sample_ct) && x_nosex_geno_cts) {
                // shouldn't be possible for now, since x_nosex_geno_cts can
                // only be non-null when nonfounders is true
                assert(0);
                // founder_x_nosex_geno_cts = x_nosex_geno_cts;
              } else {
                if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(uint32_t, 3, x_len, &founder_x_nosex_geno_cts))) {
                  goto Plink2Core_ret_NOMEM;
                }
              }
            }
          }
        }
        unsigned char* bigstack_mark_read_freqs = g_bigstack_base;
        regular_freqcounts_needed = regular_freqcounts_needed || variant_missing_hc_cts || variant_missing_dosage_cts || variant_hethap_cts || raw_geno_cts || founder_raw_geno_cts || imp_r2_vals;
        uintptr_t* variant_afreqcalc = nullptr;
        uint32_t afreqcalc_variant_ct = 0;
        if (allele_freqs) {
          if (pcp->read_freq_fname) {
            reterr = ReadAlleleFreqs(variant_include, variant_ids, allele_idx_offsets, allele_storage, pcp->read_freq_fname, raw_variant_ct, variant_ct, max_allele_ct, max_variant_id_slen, max_allele_slen, pcp->af_pseudocount, pcp->max_thread_ct, allele_freqs, &variant_afreqcalc);
            if (unlikely(reterr)) {
              goto Plink2Core_ret_1;
            }
            if (variant_ct) {
              afreqcalc_variant_ct = PopcountWords(variant_afreqcalc, raw_variant_ctl);
            }
          } else {
            variant_afreqcalc = variant_include;
            afreqcalc_variant_ct = variant_ct;
          }
        } else if (pcp->read_freq_fname) {
          logerrputs("Warning: Ignoring --read-freq since no command would use the frequencies.\n");
        }
        if (regular_freqcounts_needed || afreqcalc_variant_ct) {
          if (unlikely(pcp->misc_flags & kfMiscErrorOnFreqCalc)) {
            // This needs to be kept in sync with DecentAlleleFreqsAreNeeded(),
            // IndecentAlleleFreqsAreNeeded(),
            // [Founder]AlleleDosagesAreNeeded(),
            // VariantMissingHcCtsAreNeeded(),
            // VariantMissingDosageCtsAreNeeded(),
            // [Founder]RawGenoCtsAreNeeded(), TrimAlts, and is_minimac3_r2.
            logerrputs("Error: --error-on-freq-calc specified, but allele frequency calculation is\nneeded.\nFlags which may invoke the allele frequency calculation include --check-sex,\n--freq, --geno, --geno-counts, --genotyping-rate, --glm (unless 'omit-ref' is\nspecified), --hardy, --het (unless 'small-sample' is specified), --hwe,\n--impute-sex, --indep-pairwise, the 'trim-alts' modifier of\n--make-[b]pgen/--make-bed, --make-grm-{bin,list}, --make-rel, --[max-]mac,\n--[max-]maf, --minimac3-r2-filter, --missing, --pca, --score[-list] (unless\n'no-mean-imputation' is specified, and neither 'center' nor\n'variance-standardize' are), and --variant-score.\n");
            goto Plink2Core_ret_INVALID_CMDLINE;
          }
          // note that --geno depends on different handling of X/Y than --maf.

          // possible todo: "free" these arrays early in some cases
          // todo: oblig-missing

          // bugfix (22 Nov 2017): given "--genotyping-rate dosage" on a file
          // with no dosage data, variant_missing_hc_cts may be nullptr while
          // variant_missing_dosage_cts is zero-initialized.  In this case,
          // it's fine to pass variant_missing_dosage_cts in the
          // hardcall-missing-count slot... and it's NOT fine to pass in
          // nullptrs for both missing-count arrays...
          const uint32_t dosageless_file = !(pgfi.gflags & kfPgenGlobalDosagePresent);
          reterr = LoadAlleleAndGenoCounts(sample_include, founder_info, sex_nm, sex_male, regular_freqcounts_needed? variant_include : variant_afreqcalc, cip, allele_idx_offsets, raw_sample_ct, sample_ct, founder_ct, male_ct, nosex_ct, raw_variant_ct, regular_freqcounts_needed? variant_ct : afreqcalc_variant_ct, first_hap_uidx, is_minimac3_r2, (pcp->misc_flags / kfMiscYNosexMissingStats) & 1, pcp->max_thread_ct, pgr_alloc_cacheline_ct, &pgfi, allele_presents, allele_ddosages, founder_allele_ddosages, ((!variant_missing_hc_cts) && dosageless_file)? variant_missing_dosage_cts : variant_missing_hc_cts, dosageless_file? nullptr : variant_missing_dosage_cts, variant_hethap_cts, raw_geno_cts, founder_raw_geno_cts, x_male_geno_cts, founder_x_male_geno_cts, x_nosex_geno_cts, founder_x_nosex_geno_cts, imp_r2_vals);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
          if (overlapping_allele_ddosages) {
            founder_allele_ddosages = allele_ddosages;
          }
          }
        if (allele_freqs) {
          if (afreqcalc_variant_ct) {
            ComputeAlleleFreqs(variant_afreqcalc, allele_idx_offsets, nonfounders? allele_ddosages : founder_allele_ddosages, afreqcalc_variant_ct, pcp->af_pseudocount, allele_freqs);
          }
        }
        BigstackReset(bigstack_mark_read_freqs);


        if ((pcp->command_flags1 & kfCommand1Hardy) || (pcp->hwe_ln_thresh != -DBL_MAX)) {
          if (cip->haploid_mask[0] & 1) {
            if (unlikely(pcp->command_flags1 & kfCommand1Hardy)) {
              logerrputs("Error: --hardy is pointless on an all-haploid genome.\n");
              goto Plink2Core_ret_INCONSISTENT_INPUT;
            }
            // could make hwe_ln_thresh non-const and set it to -DBL_MAX
            // earlier on in this case
            logerrputs("Warning: --hwe has no effect since entire genome is haploid.\n");
          } else {
            STD_ARRAY_PTR_DECL(uint32_t, 3, hwe_x_male_geno_cts) = nonfounders? x_male_geno_cts : founder_x_male_geno_cts;

            // For multiallelic variants, perform one 'biallelic' test per
            // allele.  This is much, much faster than e.g. triallelic exact
            // tests, and the resulting set of p-values are arguably more
            // useful anyway.

            // One entry per autosomal multiallelic alt allele.
            // Only two cells per entry, since third can be inferred by
            // subtracting from the nonmissing genotype count (this is
            // potentially a memory hog).
            STD_ARRAY_PTR_DECL(uint32_t, 2, autosomal_xgeno_cts) = nullptr;
            // one entry per chrX multiallelic alt allele.
            STD_ARRAY_PTR_DECL(uint32_t, 2, x_knownsex_xgeno_cts) = nullptr;
            STD_ARRAY_PTR_DECL(uint32_t, 2, x_male_xgeno_cts) = nullptr;
            uint32_t hwe_x_ct = 0;
            if (hwe_x_probs_needed) {
              hwe_x_ct = CountChrVariantsUnsafe(variant_include, cip, cip->xymt_codes[kChrOffsetX]);
              // hwe_x_ct == 0 possible when hwe_x_probs_needed set, if --geno
              // filters out all chrX variants
            }
            uintptr_t autosomal_xallele_ct = 0;
            uintptr_t x_xallele_ct = 0;
            if (allele_idx_offsets) {
              // bugfix (3 Mar 2023): if !hwe_x_probs_needed, need to subtract
              // chrX variant count here
              const uint32_t chr_ct = cip->chr_ct;
              for (uint32_t chr_fo_idx = 0; chr_fo_idx != chr_ct; ++chr_fo_idx) {
                const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
                if ((chr_idx != x_code) && (!IsSet(cip->haploid_mask, chr_idx))) {
                  autosomal_xallele_ct += CountExtraAlleles(variant_include, allele_idx_offsets, cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1], 1);
                }
              }
              if (hwe_x_ct) {
                x_xallele_ct = CountExtraAlleles(variant_include, allele_idx_offsets, x_start, x_start + x_len, 1);
              }
              if (autosomal_xallele_ct || x_xallele_ct) {
                if (autosomal_xallele_ct) {
                  if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(uint32_t, 2, autosomal_xallele_ct, &autosomal_xgeno_cts))) {
                    goto Plink2Core_ret_NOMEM;
                  }
                }
                if (x_xallele_ct) {
                  if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(uint32_t, 2, x_xallele_ct, &x_knownsex_xgeno_cts))) {
                    goto Plink2Core_ret_NOMEM;
                  }
                  if (hwe_x_male_geno_cts) {
                    if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(uint32_t, 2, x_xallele_ct, &x_male_xgeno_cts))) {
                      goto Plink2Core_ret_NOMEM;
                    }
                  }
                }
                if (unlikely(reterr)) {
                  goto Plink2Core_ret_1;
                }
              }
            }

            }
        }
        // raw_geno_cts/founder_raw_geno_cts/hwe_x_ln_pvals no longer needed
        BigstackReset(bigstack_mark_geno_cts);


        }


      if (pcp->filter_flags & kfFilterPvarReq) {
        if (unlikely(!variant_ct)) {
          // do we want this to be conditionally acceptable?
          logerrputs("Error: No variants remaining after main filters.\n");
          goto Plink2Core_ret_DEGENERATE_DATA;
        }
        logprintf("%u variant%s remaining after main filters.\n", variant_ct, (variant_ct == 1)? "" : "s");
      }

      if ((pcp->command_flags1 & kfCommand1MendelReport) || (pcp->filter_flags & kfFilterMendel)) {
        if (vpos_sortstatus & kfUnsortedVarBp) {
          logerrputs("Error: --me and --mendel require a sorted .pvar/.bim.  Retry this command after\nusing --make-pgen/--make-bed + --sort-vars to sort your data.\n");
          goto Plink2Core_ret_INCONSISTENT_INPUT;
        }
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
        if (pcp->filter_flags & kfFilterMendel) {
          sample_ct = PopcountWords(sample_include, raw_sample_ctl);
          variant_ct = PopcountWords(variant_include, raw_variant_ctl);
          logprintf("--me: %u sample%s and %u variant%s remaining.\n", sample_ct, (sample_ct == 1)? "" : "s", variant_ct, (variant_ct == 1)? "" : "s");
          UpdateSampleSubsets(sample_include, raw_sample_ct, sample_ct, founder_info, &founder_ct, sex_nm, sex_male, &male_ct, &nosex_ct);
        }
      }

      #ifndef NOLAPACK
      #endif

      if (pcp->command_flags1 & (kfCommand1MakePlink2 | kfCommand1Exportf | kfCommand1WriteCovar)) {
        // If non-null, allele_permute[allele_idx_offsets[j]] stores the old
        // allele-index of the new ref allele,
        // allele_permute[allele_idx_offsets[j] + 1] stores the old
        // allele-index of the new alt1 allele, etc.
        // Operations which instantiate this (--maj-ref, --ref-allele,
        // --alt[1]-allele) are only usable with fileset creation commands.  No
        // more pass-marker_reverse-to-everything nonsense.
        // (Technically, I could also drop support for --export, but that would
        // force too many real-world jobs to require two plink2 runs instead of
        // one.)

        uint32_t* variant_bps_backup = nullptr;
        const char** allele_storage_backup = nullptr;
        uint32_t max_allele_slen_backup = max_allele_slen;
        uintptr_t* nonref_flags_backup = nullptr;
        uint32_t nonref_flags_was_null = (nonref_flags == nullptr);
        AlleleCode* allele_permute = nullptr;
        if (allele_presents) {
          if (!allele_permute) {
            if (unlikely(bigstack_end_alloc_ac(raw_allele_ct, &allele_permute))) {
              goto Plink2Core_ret_NOMEM;
            }
            reterr = InitAllelePermuteUnsafe(allele_idx_offsets, raw_variant_ct, pcp->max_thread_ct, allele_permute);
            if (unlikely(reterr)) {
              goto Plink2Core_ret_1;
            }
          }
          reterr = TrimAllelePermute(variant_include, allele_idx_offsets, allele_presents, variant_ct, pcp->max_thread_ct, allele_permute);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        }
        BigstackReset(bigstack_mark_allele_ddosages);

        // update (18 Mar 2018): permit no-covariate --write-covar when used to
        // write phenotypes.
        if (variant_bps_backup) {
          memcpy(variant_bps, variant_bps_backup, raw_variant_ct * sizeof(int32_t));
        }
        if (allele_storage_backup) {
          max_allele_slen = max_allele_slen_backup;
          memcpy(allele_storage, allele_storage_backup, raw_allele_ct * sizeof(intptr_t));
        }
        if (nonref_flags_backup) {
          memcpy(nonref_flags, nonref_flags_backup, raw_variant_ctl * sizeof(intptr_t));
        } else if (nonref_flags_was_null) {
          nonref_flags = nullptr;
          pgfi.nonref_flags = nullptr;
          GET_PRIVATE(simple_pgr, m).fi.nonref_flags = nullptr;
        }
      }
      BigstackReset(bigstack_mark_allele_ddosages);

      if (pcp->command_flags1 & kfCommand1Score) {
        reterr = ScoreReport(sample_include, &pii.sii, sex_nm, sex_male, pheno_cols, pheno_names, variant_include, cip, variant_ids, allele_idx_offsets, allele_storage, allele_freqs, &(pcp->score_info), pcp->output_missing_pheno, raw_sample_ct, sample_ct, nosex_ct, pheno_ct, max_pheno_name_blen, raw_variant_ct, variant_ct, max_variant_id_slen, pcp->xchr_model, pcp->max_thread_ct, &simple_pgr, outname, outname_end);
        if (unlikely(reterr)) {
          goto Plink2Core_ret_1;
        }
        const uint32_t give_best_pgs = (pcp->score_info.flags / kfScoreGiveBestPgs) & 1;
        const uint32_t corr_matrix = (pcp->score_info.flags / kfScoreGiveCorrMatrix) & 1;
        if (give_best_pgs || corr_matrix) {
          // ScoreReport has written <outname>.sscore[.zst]; read it back.  Its
          // internal score matrix is bigstack scratch, freed before the
          // function returns, so the file is the only durable form.
          // --q-score-range is not supported by either report (multiple
          // output files).
          if (pcp->score_info.qsr_range_fname) {
            logerrprintf("Error: %s cannot be combined with --q-score-range.\n",
                         give_best_pgs? "--give-best-pgs" : "--give-corr-matrix");
            reterr = kPglRetInvalidCmdline;
            goto Plink2Core_ret_1;
          }
          const uint32_t output_zst = (pcp->score_info.flags / kfScoreZs) & 1;
          OutnameZstSet(".sscore", output_zst, outname_end);
          reterr = BestPgsReport(outname, sample_include, sample_ct,
                                 pheno_cols, pheno_names, max_pheno_name_blen,
                                 pheno_ct, covar_cols, covar_names,
                                 max_covar_name_blen, covar_ct,
                                 pcp->pgs_prevalence, pcp->pgs_folds,
                                 pcp->pgs_fold_seed,
                                 give_best_pgs, corr_matrix,
                                 (pcp->score_info.flags / kfScoreBestPgsNoNa) & 1,
                                 pcp->pgs_pheno_pairs_fname,
                                 outname, outname_end);
          if (unlikely(reterr)) {
            goto Plink2Core_ret_1;
          }
        }
      }
      // eventually check for nonzero pheno_ct here?

      }
    }
  while (0) {
  Plink2Core_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  Plink2Core_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  Plink2Core_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  Plink2Core_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  Plink2Core_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  Plink2Core_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  Plink2Core_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 Plink2Core_ret_1:
  if (loop_cats_pheno_col) {
    // Current implementation requires this to happen before CleanupPhenoCols()
    // on pheno_cols/covar_cols, since loop_cats_pheno_col actually points to
    // the last element of one of those arrays.
    vecaligned_free_cond(loop_cats_pheno_col->nonmiss);
  }
  CleanupPhenoCols(covar_ct, covar_cols);
  CleanupPhenoCols(pheno_ct, pheno_cols);
  free_cond(covar_names);
  free_cond(pheno_names);
  CleanupPgr2(".pgen file", &simple_pgr, &reterr);
  free_cond(ext_slot.contents);
  CleanupPgfi2(".pgen file", &pgfi, &reterr);
  assert(pgfi.block_base == nullptr);
  // no BigstackReset() needed?
  return reterr;
}
// flagname_p only needed when check_file_existence true
void PrintVer() {
  fputs(ver_str, stdout);
  fputs(ver_str2, stdout);
}
static_assert(sizeof(int) == sizeof(int32_t), "main() assumes int and int32_t are synonymous.");
static_assert(!kChrOffsetX, "--autosome-num/--chr-set/--cow/etc. assume kChrOffsetX == 0.");
static_assert(kChrOffsetY == 1, "--chr-set/--cow/... assume kChrOffsetY == 1.");
static_assert(kChrOffsetXY == 2, "--chr-set/--cow/... assume kChrOffsetXY == 2.");
static_assert(kChrOffsetMT == 3, "--chr-set/--cow/... assume kChrOffsetMT == 3.");
#ifdef __cplusplus
}  // namespace plink2
#endif

#if defined(CPU_CHECK_SSE42) || defined(CPU_CHECK_AVX2)
int RealMain(int argc, char** argv) {
#else
int main(int argc, char** argv) {
#endif
#ifdef __cplusplus
  using namespace plink2;
#endif
  flush_denormals();
  // special case, since it may dump to stdout
  if (argc > 1) {
    const char* argv1 = argv[1];
    const uint32_t argv1_slen = strlen(argv1);
      }

  unsigned char* bigstack_ua = nullptr;
  Plink2CmdlineMeta pcm;
  PreinitPlink2CmdlineMeta(&pcm);
  const char* flagname_p = nullptr;
  char* const_fid = nullptr;
  char* import_single_chr_str = nullptr;
  char* ox_missing_code = nullptr;
  char* vcf_dosage_import_field = nullptr;
  uint32_t* rseeds = nullptr;
  LlStr* file_delete_list = nullptr;
  uint32_t arg_idx = 0;
  uint32_t print_end_time = 0;
  uint32_t warning_errcode = 0;
  PglErr reterr = kPglRetSuccess;
  Plink2Cmdline pc;
  pc.filter_flags = kfFilter0;
  pc.dependency_flags = kfFilter0;
  pc.load_filter_log_flags = kfLoadFilterLog0;
  pc.pginame = nullptr;
  pc.var_filter_exceptions_flattened = nullptr;
  pc.varid_template_str = nullptr;
  pc.varid_multi_template_str = nullptr;
  pc.varid_multi_nonsnp_template_str = nullptr;
  pc.missing_varid_match = nullptr;
  pc.varid_from = nullptr;
  pc.varid_to = nullptr;
  pc.varid_snp = nullptr;
  pc.varid_exclude_snp = nullptr;
  pc.pheno_fname = nullptr;
  pc.covar_fname = nullptr;
  pc.sample_sort_fname = nullptr;
  pc.keep_fnames = nullptr;
  pc.keepfam_fnames = nullptr;
  pc.remove_fnames = nullptr;
  pc.removefam_fnames = nullptr;
  pc.extract_fnames = nullptr;
  pc.extract_intersect_fnames = nullptr;
  pc.exclude_fnames = nullptr;
  pc.freq_ref_binstr = nullptr;
  pc.freq_alt1_binstr = nullptr;
  pc.glm_local_covar_fname = nullptr;
  pc.glm_local_pvar_fname = nullptr;
  pc.glm_local_psam_fname = nullptr;
  pc.read_freq_fname = nullptr;
  pc.within_fname = nullptr;
  pc.catpheno_name = nullptr;
  pc.family_missing_catname = nullptr;
  pc.keep_cats_fname = nullptr;
  pc.keep_cat_names_flattened = nullptr;
  pc.keep_cat_phenoname = nullptr;
  pc.remove_cats_fname = nullptr;
  pc.remove_cat_names_flattened = nullptr;
  pc.remove_cat_phenoname = nullptr;
  pc.split_cat_phenonames_flattened = nullptr;
  pc.require_pheno_flattened = nullptr;
  pc.require_covar_flattened = nullptr;
  pc.vstd_flattened = nullptr;
  pc.quantnorm_flattened = nullptr;
  pc.covar_quantnorm_flattened = nullptr;
  pc.loop_cats_phenoname = nullptr;
  pc.fa_fname = nullptr;
  pc.require_info_flattened = nullptr;
  pc.require_no_info_flattened = nullptr;
  pc.keep_col_match_fname = nullptr;
  pc.keep_col_match_flattened = nullptr;
  pc.keep_col_match_name = nullptr;
  pc.ref_allele_flag = nullptr;
  pc.alt_allele_flag = nullptr;
  pc.update_chr_flag = nullptr;
  pc.update_map_flag = nullptr;
  pc.update_name_flag = nullptr;
  pc.update_sample_ids_fname = nullptr;
  pc.update_parental_ids_fname = nullptr;
  pc.recover_var_ids_fname = nullptr;
  pc.vscore_fname = nullptr;
  pc.sscore_fname = nullptr;
  pc.bim_fname = nullptr;
  pc.pgs_model_fname = nullptr;
  pc.pgs_models_list_fname = nullptr;
  pc.pgs_pheno_pairs_fname = nullptr;
  pc.pgs_model_col_chr     = 0;
  pc.pgs_model_col_pos     = 0;
  pc.pgs_model_col_effect  = 0;
  pc.pgs_model_col_other   = 0;
  pc.pgs_model_col_weight  = 0;
  pc.pgs_model_col_freq    = 0;
  pc.pgs_model_allow_ambiguous = 0;
  pc.pgs_model_resolve_ambig   = 0;
  pc.pgs_model_maf_threshold   = 0.4;
  pc.indep_preferred_fname = nullptr;
  pc.not_pheno_flattened = nullptr;
  pc.not_covar_flattened = nullptr;
  pc.indv_str = nullptr;
  pc.rename_chrs_fname = nullptr;
  pc.zero_cluster_fname = nullptr;
  pc.zero_cluster_phenoname = nullptr;
  InitRangeList(&pc.snps_range_list);
  InitRangeList(&pc.exclude_snps_range_list);
  InitRangeList(&pc.pheno_range_list);
  InitRangeList(&pc.covar_range_list);
  InitScore(&pc.score_info);
  InitCmpExpr(&pc.keep_if_expr);
  InitCmpExpr(&pc.remove_if_expr);
  InitCmpExpr(&pc.extract_if_info_expr);
  InitCmpExpr(&pc.exclude_if_info_expr);
  InitExtractColCond(&pc.extract_col_cond_info);
  ChrInfo chr_info;
  if (unlikely(InitChrInfo(&chr_info))) {
    goto main_ret_NOMEM_NOLOG;
  }

  {
    // standardize strtod() behavior
    // (not necessary any more since we use dtoa_g() instead)
    // setlocale(LC_NUMERIC, "C");

    uint32_t first_arg_idx;
    uint32_t flag_ct;
    reterr = CmdlineParsePhase1(ver_str, ver_str2, PROG_NAME_STR, notestr_null_calc2, kCmdlineFormatStr, errstr_append, kMaxFlagBlen, DispHelp, &argc, &argv, &pcm, &first_arg_idx, &flag_ct);
    if (reterr) {
      if (reterr == kPglRetHelp) {
        reterr = kPglRetSuccess;
      }
      goto main_ret_NOLOG;
    }
    if (!flag_ct) {
      goto main_ret_NULL_CALC_0;
    }
    // + (kBytesPerWord - 1) to support overread during sort
    if (unlikely(pgl_malloc(flag_ct * kMaxFlagBlen + kBytesPerWord - 1, &pcm.flag_buf) ||
                 pgl_malloc(flag_ct * sizeof(int32_t), &pcm.flag_map))) {
      goto main_ret_NOMEM_NOLOG2;
    }

    // No modifications to argv past this point.
    const char* const* argvk = TO_CONSTCPCONSTP(argv);

    char* flagname_write_iter = pcm.flag_buf;
    uint32_t cur_flag_idx = 0;
    for (arg_idx = first_arg_idx; arg_idx != S_CAST(uint32_t, argc); ++arg_idx) {
      flagname_p = IsCmdlineFlagStart(argvk[arg_idx]);
      if (flagname_p) {
        const uint32_t flag_slen = strlen(flagname_p);
        // Flag aliases.  Only the two --threads spellings survive; every other
        // alias in stock plink2 pointed at a command this build does not have.
        if (strequal_k(flagname_p, "num_threads", flag_slen) ||
            strequal_k(flagname_p, "thread-num", flag_slen)) {
          snprintf(flagname_write_iter, kMaxFlagBlen, "threads");
        } else {
          memcpy(flagname_write_iter, flagname_p, flag_slen + 1);
        }
        flagname_write_iter = &(flagname_write_iter[kMaxFlagBlen]);
        pcm.flag_map[cur_flag_idx++] = arg_idx;
      }
    }
    char outname[kPglFnamesize];
    memcpy_k(outname, PROG_NAME_STR, sizeof(PROG_NAME_STR) - 1);
    char* outname_end = nullptr;
    char range_delim;
    uint32_t strict_extra_chr;
    int32_t known_procs;
    reterr = CmdlineParsePhase2(ver_str, errstr_append, argvk, sizeof(PROG_NAME_STR) - 1, kMaxFlagBlen, argc, flag_ct, &pcm, outname, &outname_end, &range_delim, &strict_extra_chr, &known_procs, &pc.max_thread_ct);
    if (unlikely(reterr)) {
      goto main_ret_NOLOG;
    }

    char pgenname[kPglFnamesize];
    char psamname[kPglFnamesize];
    char pvarname[kPglFnamesize];
    pgenname[0] = '\0';
    psamname[0] = '\0';
    pvarname[0] = '\0';
    snprintf(pc.missing_catname, kMaxMissingPhenostrBlen, "NONE");
    snprintf(pc.output_missing_pheno, kMaxMissingPhenostrBlen, "NA");
    snprintf(pc.legacy_output_missing_pheno, kMaxMissingPhenostrBlen, "-9");
    cur_flag_idx = 0;
    pc.command_flags1 = kfCommand10;
    // uint64_t command_flags2 = 0;
    pc.misc_flags = strict_extra_chr? kfMiscProhibitExtraChr : kfMisc0;
    pc.pvar_psam_flags = kfPvarPsam0;
    pc.fam_cols = kfFamCol13456;
    pc.select_sid_missingness_mode = kSelectSidMissingness0;
    pc.select_sid_tiebreak_mode = kSelectSidTiebreak0;
    for (uint32_t uii = 0; uii != 4; ++uii) {
      pc.filter_modes[uii] = kFreqFilterNonmajor;
    }

    // Default value is 1638 = 32768 / 20, and that's applied to imported
    // dosages when --hard-call-threshold is not specified.
    // However, when --make-[b]pgen is run on a dosage-containing dataset,
    // explicit --hard-call-threshold will cause the hardcall set to be
    // regenerated, and that won't happen without --hard-call-threshold.  So we
    // need to distinguish between --hard-call-threshold 0.1 and no flag.
    pc.hard_call_thresh = UINT32_MAX;

    pc.dosage_erase_thresh = 0;
    pc.geno_thresh = 1.0;
    pc.hwe_ln_thresh = -DBL_MAX;
    pc.mach_r2_max = 0.0;
    pc.minimac3_r2_max = 0.0;
    pc.af_pseudocount = 0.0;
    pc.min_maf = 0.0;
    pc.max_maf = 1.0;
    pc.pgs_prevalence = -1.0;  // < 0 means not supplied
    pc.pgs_folds = 1;
    // Fold assignment defaults to a fixed seed rather than the clock: every
    // other number this tool prints is a deterministic function of its inputs,
    // and a cross-validated R^2 that moved between identical runs would not be.
    // --seed overrides it.
    pc.pgs_fold_seed = 0;
    pc.min_allele_ddosage = 0;
    pc.max_allele_ddosage = (~0LLU);
    pc.var_min_qual = -1;
    pc.new_variant_id_max_allele_slen = 23;
    pc.splitpar_bound1 = 0;
    pc.splitpar_bound2 = 0;
    pc.missing_pheno = -9;
    pc.mhc_region = kMhcNone;
    pc.window_bp = -1;
    pc.xchr_model = 2;
    pc.keep_col_match_num = 0;
    pc.filter_min_allele_ct = 0;
    pc.filter_max_allele_ct = UINT32_MAX;
    pc.input_missing_geno_char = '0';
    uintptr_t malloc_size_mib = 0;
    LoadParams load_params = kfLoadParams0;
    Xload xload = kfXload0;
    uint32_t rseed_ct = 0;
    MakePlink2Flags make_plink2_flags = kfMake0;
    char id_delim = '\0';
    char idspace_to = '\0';
    ImportFlags import_flags = kfImport0;
    uint32_t delete_pmerge_result = 0;
    uint32_t score_col_nums_present = 0;
    uint32_t r2_required = 0;
    uint32_t memory_require = 0;
#ifdef USE_MKL
    uint32_t mkl_native = 0;
#endif
    uint32_t pmerge_required = 0;
    uint32_t allow_misleading_out_arg = 0;
    uint32_t import_max_allele_ct = 0x7ffffffe;
    // In retrospect, plink 1.9's per-flag perm-count modifiers were a bit
    // silly; reasonable to drop that bit of flexibility so we can organize
    // permutation-test formatting under a single --perm-format flag.
    //
    // We error out on the very-unlikely case where that flexibility was
    // actually used (i.e. multiple association tests requested, at least one
    // with perm-count and at least one without, and --perm-count not
    // specified), and otherwise just treat perm-count modifiers (and the
    // --perm-count flag) like "--perm-format counts".
    uint32_t explicit_perm_count_modifier = 0;
    uint32_t implicit_perm_count_omission = 0;
    do {
      flagname_p = &(pcm.flag_buf[cur_flag_idx * kMaxFlagBlen]);
      const char* flagname_p2 = &(flagname_p[1]);
      arg_idx = pcm.flag_map[cur_flag_idx];
      uint32_t param_ct = GetParamCt(argvk, argc, arg_idx);
      switch (*flagname_p) {
      case '1':
        if (likely(*flagname_p2 == '\0')) {
          pc.misc_flags |= kfMiscAffection01;
          goto main_param_zero;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'a':
        if (strequal_k_unsafe(flagname_p2, "llow-ambiguous")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 0))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          pc.pgs_model_allow_ambiguous = 1;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'b':
        if (strequal_k_unsafe(flagname_p2, "im")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.bim_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "file")) {
          if (unlikely(xload)) {
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_modif_idx = 1;
          if (param_ct == 2) {
            if (unlikely(CheckExtraParam(&(argvk[arg_idx]), "vzs", &fname_modif_idx))) {
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          const char* fname_prefix = argvk[arg_idx + fname_modif_idx];
          const uint32_t slen = strlen(fname_prefix);
          // pathological case bugfix (1 Feb 2018): need to subtract 1 more, to
          // avoid buffer overflow in the case we rename and append '~'.
          if (unlikely(slen > (kPglFnamesize - 10))) {
            // could use kPglFnamesize - 2 - 3 * param_ct, but that's pointless
            logerrputs("Error: --bfile argument too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          snprintf(memcpya(pgenname, fname_prefix, slen), 9, ".bed");
          snprintf(memcpya(psamname, fname_prefix, slen), 9, ".fam");
          char* bimname_end = memcpya(pvarname, fname_prefix, slen);
          bimname_end = Stpcpy(bimname_end, ".bim");
          if (param_ct == 2) {
            snprintf(bimname_end, 5, ".zst");
          }
          load_params |= kfLoadParamsPfileAll;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'c':
        if (strequal_k_unsafe(flagname_p2, "ovar")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_idx = 1;
          if (param_ct == 2) {
            if (unlikely(CheckExtraParam(&(argvk[arg_idx]), "iid-only", &fname_idx))) {
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.misc_flags |= kfMiscCovarIidOnly;
          }
          reterr = AllocFname(argvk[arg_idx + fname_idx], flagname_p, &pc.covar_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "ovar-name")) {
          if (unlikely(pc.covar_range_list.name_ct)) {
            logerrputs("Error: --covar-name can't be used with --covar-col-nums.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          // can now be used without --covar
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 0, range_delim, &pc.covar_range_list);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPsamReq;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'd':
        if (unlikely(*flagname_p2)) {
          // --d is preprocessed
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'e':
        if (strequal_k_unsafe(flagname_p2, "xtract")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(IsIntervalBedModifier(cur_modif, strlen(cur_modif)))) {
            logerrprintfww("Error: '--extract %s' is not supported; --extract takes variant ID lists only. Convert your intervals to a variant ID list by filtering the .bim/.pvar on position first.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_A;
          }
          reterr = AllocAndFlattenFnames(&(argvk[arg_idx + 1]), flagname_p, param_ct, &pc.extract_fnames);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else if (strequal_k_unsafe(flagname_p2, "xclude")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          const char* cur_modif = argvk[arg_idx + 1];
          if (unlikely(IsIntervalBedModifier(cur_modif, strlen(cur_modif)))) {
            logerrprintfww("Error: '--exclude %s' is not supported; --exclude takes variant ID lists only. Convert your intervals to a variant ID list by filtering the .bim/.pvar on position first.\n", cur_modif);
            goto main_ret_INVALID_CMDLINE_A;
          }
          reterr = AllocAndFlattenFnames(&(argvk[arg_idx + 1]), flagname_p, param_ct, &pc.exclude_fnames);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPvarReq;
        } else {
          // The four MHC exclusions share a shape, so match them off a table
          // rather than as four near-identical branches.
          const char* const mhc_suffixes[] = {
            "xclude-mhc-hg19", "xclude-xmhc-hg19",
            "xclude-mhc-hg38", "xclude-xmhc-hg38"
          };
          const MhcRegion mhc_regions[] = {
            kMhcHg19, kMhcXHg19, kMhcHg38, kMhcXHg38
          };
          uint32_t mhc_idx = 4;
          for (uint32_t uii = 0; uii != 4; ++uii) {
            if (!strcmp(flagname_p2, mhc_suffixes[uii])) {
              mhc_idx = uii;
              break;
            }
          }
          if (unlikely(mhc_idx == 4)) {
            goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 0))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          // Two of these together is always a mistake: same build is redundant
          // (xMHC contains MHC), different builds is a contradiction.
          if (unlikely(pc.mhc_region != kMhcNone)) {
            snprintf(g_logbuf, kLogbufSize, "Error: %s cannot be used with %s.\n", MhcFlagName(mhc_regions[mhc_idx]), MhcFlagName(pc.mhc_region));
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          pc.mhc_region = mhc_regions[mhc_idx];
          pc.filter_flags |= kfFilterPvarReq;
        }
        break;

      case 'f':
        if (strequal_k_unsafe(flagname_p2, "olds")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          // ScanPosintCappedx, not ScanUintCappedx: 0 folds is wrong on its
          // face, so it is refused here rather than after the pgen, phenotypes
          // and covariates have loaded -- which is where the range test
          // against the complete-case count genuinely has to happen, since
          // BestPgsValidate is the first place that count exists.
          if (unlikely(ScanPosintCappedx(argvk[arg_idx + 1], 0x7fffffff, &pc.pgs_folds))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --folds argument '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'g':
        if (strequal_k_unsafe(flagname_p2, "ive-best-pgs")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            const char* modif = argvk[arg_idx + 1];
            if (unlikely(!strequal_k(modif, "no-na", strlen(modif)))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --give-best-pgs argument '%s'; 'no-na' is the only one accepted.\n", modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            pc.score_info.flags |= kfScoreBestPgsNoNa;
          }
          pc.score_info.flags |= kfScoreGiveBestPgs;
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "ive-corr-matrix")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 0))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          pc.score_info.flags |= kfScoreGiveCorrMatrix;
          pc.dependency_flags |= kfFilterPsamReq;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'h':
        goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        break;

      case 'i':
        goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        break;

      case 'k':
        if (strequal_k_unsafe(flagname_p2, "eep")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlattenFnames(&(argvk[arg_idx + 1]), flagname_p, param_ct, &pc.keep_fnames);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'l':
        goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        break;

      case 'm':
        if (strequal_k_unsafe(flagname_p2, "emory")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t mb_modif_idx = 1;
          if (param_ct == 2) {
            if (unlikely(CheckExtraParam(&(argvk[arg_idx]), "require", &mb_modif_idx))) {
              goto main_ret_INVALID_CMDLINE_A;
            }
            memory_require = 1;
          }
          const char* mb_modif = argvk[arg_idx + mb_modif_idx];
          if (unlikely(ScanPosintptrx(mb_modif, &malloc_size_mib))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --memory argument '%s'.\n", mb_modif);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (unlikely(malloc_size_mib < kBigstackMinMib)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --memory argument '%s' (minimum %u).\n", mb_modif, kBigstackMinMib);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
#ifndef __LP64__
          if (unlikely(malloc_size_mib > kMalloc32bitMibMax)) {
            logerrprintf("Error: --memory argument too large for 32-bit version (max %u).\n", kMalloc32bitMibMax);
            goto main_ret_INVALID_CMDLINE;
          }
#endif
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'n':
        goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        break;

      case 'o':
        if (unlikely(!strequal_k_unsafe(flagname_p2, "ut"))) {
          // --out is a special case due to logging
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'p':
        if (strequal_k_unsafe(flagname_p2, "file")) {
          if (unlikely(load_params || xload)) {
            // currently only possible with --bcf, --bfile, --pfile
            goto main_ret_INVALID_CMDLINE_INPUT_CONFLICT;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_modif_idx = 1;
          if (param_ct == 2) {
            if (unlikely(CheckExtraParam(&(argvk[arg_idx]), "vzs", &fname_modif_idx))) {
              goto main_ret_INVALID_CMDLINE_A;
            }
          }
          const char* fname_prefix = argvk[arg_idx + fname_modif_idx];
          const uint32_t slen = strlen(fname_prefix);
          if (unlikely(slen > (kPglFnamesize - 11))) {
            logerrputs("Error: --pfile argument too long.\n");
            goto main_ret_OPEN_FAIL;
          }
          snprintf(memcpya(pgenname, fname_prefix, slen), 10, ".pgen");
          snprintf(memcpya(psamname, fname_prefix, slen), 10, ".psam");
          char* pvarname_end = memcpya(pvarname, fname_prefix, slen);
          pvarname_end = Stpcpy(pvarname_end, ".pvar");
          if (param_ct == 2) {
            snprintf(pvarname_end, 5, ".zst");
          }
          load_params |= kfLoadParamsPfileAll;
        } else if (strequal_k_unsafe(flagname_p2, "heno")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 2))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          uint32_t fname_idx = 1;
          if (param_ct == 2) {
            if (unlikely(CheckExtraParam(&(argvk[arg_idx]), "iid-only", &fname_idx))) {
              goto main_ret_INVALID_CMDLINE_A;
            }
            pc.misc_flags |= kfMiscPhenoIidOnly;
          }
          reterr = AllocFname(argvk[arg_idx + fname_idx], flagname_p, &pc.pheno_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "heno-col-nums")) {
          if (unlikely(!pc.pheno_fname)) {
            logerrputs("Error: --pheno-col-nums must be used with --pheno.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.pheno_range_list.name_ct)) {
            logerrputs("Error: --pheno-col-nums cannot be used with --pheno-name.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 1, '-', &pc.pheno_range_list);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.misc_flags |= kfMiscPhenoColNums;
        } else if (strequal_k_unsafe(flagname_p2, "heno-name")) {
          if (unlikely(pc.pheno_range_list.name_ct)) {
            logerrputs("Error: --pheno-name cannot be used with --pheno-col-nums.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          // can now be used without --pheno
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 0, range_delim, &pc.pheno_range_list);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "revalence")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          double dxx;
          if (unlikely((!ScantokDouble(argvk[arg_idx + 1], &dxx)) ||
                       (dxx <= 0.0) || (dxx >= 1.0))) {
            logerrputs("Error: Invalid --prevalence argument (must be in (0, 1)).\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          pc.pgs_prevalence = dxx;
        } else if (strequal_k_unsafe(flagname_p2, "gs-pheno-pairs")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.pgs_pheno_pairs_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "gs-model") ||
                   strequal_k_unsafe(flagname_p2, "gs-models-list")) {
          // --pgs-model       <file> <chr> <pos> <effect> <other> <weight> [<freq>]
          // --pgs-models-list <list> <chr> <pos> <effect> <other> <weight> [<freq>]
          //
          // The indices mean the same thing for both, so they are parsed
          // once here; --pgs-models-list applies them to every file it
          // lists.  Only the destination filename differs.
          const uint32_t is_list =
            strequal_k_unsafe(flagname_p2, "gs-models-list");
          // A leading literal 'id' selects ID matching, which takes one key
          // column instead of chr+pos.  It cannot be mistaken for a column
          // index, so the two forms never overlap.
          const uint32_t by_id = (param_ct >= 2) &&
            strequal_k(argvk[arg_idx + 2], "id", strlen(argvk[arg_idx + 2]));
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct,
                                           by_id? 6 : 6, by_id? 7 : 7))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          pc.pgs_model_by_id = by_id;
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p,
                              is_list ? &pc.pgs_models_list_fname
                                      : &pc.pgs_model_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          // Parse the 1-based column indices: five required, the sixth
          // (source effect-allele frequency) optional.
          const uint32_t pgs_col_ct = param_ct - 1 - by_id;
          uint32_t pgs_cols[6];
          for (uint32_t pidx = 0; pidx < pgs_col_ct; ++pidx) {
            const char* cur_modif = argvk[arg_idx + 2 + by_id + pidx];
            if (unlikely(ScanPosintCappedx(cur_modif, 65535, &pgs_cols[pidx]))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s column index '%s' (must be a positive integer).\n", flagname_p, cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          // Check they are all distinct.
          for (uint32_t pi = 0; pi < pgs_col_ct; ++pi) {
            for (uint32_t pj = pi + 1; pj < pgs_col_ct; ++pj) {
              if (unlikely(pgs_cols[pi] == pgs_cols[pj])) {
                snprintf(g_logbuf, kLogbufSize, "Error: --%s column indices must all be distinct.\n", flagname_p);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
            }
          }
          if (by_id) {
            pc.pgs_model_col_id     = pgs_cols[0];
            pc.pgs_model_col_effect = pgs_cols[1];
            pc.pgs_model_col_other  = pgs_cols[2];
            pc.pgs_model_col_weight = pgs_cols[3];
            if (pgs_col_ct == 5) {
              pc.pgs_model_col_freq = pgs_cols[4];
            }
          } else {
            pc.pgs_model_col_chr    = pgs_cols[0];
            pc.pgs_model_col_pos    = pgs_cols[1];
            pc.pgs_model_col_effect = pgs_cols[2];
            pc.pgs_model_col_other  = pgs_cols[3];
            pc.pgs_model_col_weight = pgs_cols[4];
            if (pgs_col_ct == 6) {
              pc.pgs_model_col_freq = pgs_cols[5];
            }
          }
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'q':
        goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        break;

      case 'r':
        if (strequal_k_unsafe(flagname_p2, "emove")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocAndFlattenFnames(&(argvk[arg_idx + 1]), flagname_p, param_ct, &pc.remove_fnames);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.filter_flags |= kfFilterPsamReq;
        } else if (strequal_k_unsafe(flagname_p2, "ead-freq")) {
          if (unlikely(pc.command_flags1 & kfCommand1AlleleFreq)) {
            // --read-freq can't promise OBS_CTs, so simplest to just continue
            // disallowing this
            logerrputs("Error: --freq and --read-freq cannot be used together.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(pc.misc_flags & kfMiscMajRef)) {
            logerrputs("Error: --read-freq cannot be used with --maj-ref.  (We recommend preprocessing\nthe --read-freq file to generate a file that can be used with\n--ref-allele/--alt-allele.)\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.read_freq_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "esolve-ambiguous")) {
          // --resolve-ambiguous [<maf-threshold>]
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (param_ct) {
            double dxx;
            // A MAF is at most 0.5, so a larger threshold would "resolve"
            // coin flips; 0 would resolve nothing.
            if (unlikely((!ScantokDouble(argvk[arg_idx + 1], &dxx)) ||
                         (dxx <= 0.0) || (dxx > 0.5))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --resolve-ambiguous argument '%s' (must be in (0, 0.5]).\n", argvk[arg_idx + 1]);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
            pc.pgs_model_maf_threshold = dxx;
          }
          pc.pgs_model_resolve_ambig = 1;
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 's':
        if (strequal_k_unsafe(flagname_p2, "eed")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 0x7fffffff))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          rseed_ct = param_ct;
          if (unlikely(pgl_malloc(param_ct * sizeof(int32_t), &rseeds))) {
            goto main_ret_NOMEM;
          }
          for (uint32_t param_idx = 1; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            if (unlikely(ScanUintCappedx(cur_modif, UINT32_MAX, &(rseeds[param_idx - 1])))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --seed argument '%s'.\n", cur_modif);
              goto main_ret_INVALID_CMDLINE_WWA;
            }
          }
          // rseeds is consumed and freed by the sfmt initialisation below,
          // so --folds keeps its own copy.  The first value is enough: fold
          // assignment needs one reproducible 64-bit seed, not sfmt's state.
          pc.pgs_fold_seed = rseeds[0];
        } else if (strequal_k_unsafe(flagname_p2, "core")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 12))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.score_info.input_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          uint32_t numeric_param_ct = 0;
          uint32_t score_cols[4];
          for (uint32_t param_idx = 2; param_idx <= param_ct; ++param_idx) {
            const char* cur_modif = argvk[arg_idx + param_idx];
            const uint32_t cur_modif_slen = strlen(cur_modif);
            if (strequal_k(cur_modif, "header", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreHeaderIgnore;
            } else if (strequal_k(cur_modif, "header-read", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreHeaderRead;
            } else if (strequal_k(cur_modif, "no-mean-imputation", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreNoMeanimpute;
            } else if (strequal_k(cur_modif, "dominant", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreDominant;
            } else if (strequal_k(cur_modif, "recessive", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreRecessive;
            } else if (strequal_k(cur_modif, "center", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreCenter;
            } else if (strequal_k(cur_modif, "variance-standardize", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreVarianceStandardize;
            } else if (strequal_k(cur_modif, "variance-normalize", cur_modif_slen)) {
              logerrputs("Note: --score[-list]'s 'variance-normalize' modifier has been renamed to the\nmore precise 'variance-standardize'.\n");
              pc.score_info.flags |= kfScoreVarianceStandardize;
            } else if (strequal_k(cur_modif, "se", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreSe;
            } else if (strequal_k(cur_modif, "zs", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreZs;
            } else if (strequal_k(cur_modif, "ignore-dup-ids", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreIgnoreDupIds;
            } else if (strequal_k(cur_modif, "list-variants", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreListVariants;
            } else if (strequal_k(cur_modif, "list-variants-zs", cur_modif_slen)) {
              pc.score_info.flags |= kfScoreListVariants | kfScoreListVariantsZs;
            } else if (StrStartsWith(cur_modif, "cols=", cur_modif_slen)) {
              if (unlikely(pc.score_info.flags & kfScoreColAll)) {
                logerrputs("Error: Multiple --score[-list] cols= modifiers.\n");
                goto main_ret_INVALID_CMDLINE;
              }
              reterr = ParseColDescriptor(&(cur_modif[5]), "maybefid\0fid\0maybesid\0sid\0pheno1\0phenos\0nallele\0denom\0dosagesum\0scoreavgs\0scoresums\0", "score", kfScoreColMaybefid, kfScoreColDefault, 1, &pc.score_info.flags);
              if (unlikely(reterr)) {
                goto main_ret_1;
              }
            } else {
              if (unlikely(ScanPosintCappedx(cur_modif, kMaxLongLine / 2, &(score_cols[numeric_param_ct])))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid --score[-list] argument '%s'.\n", cur_modif);
                goto main_ret_INVALID_CMDLINE_WWA;
              }
              if (unlikely(numeric_param_ct == 3)) {
                logerrputs("Error: --score[-list] takes at most three numeric arguments.\n");
                goto main_ret_INVALID_CMDLINE_A;
              }
              for (uint32_t uii = 0; uii != numeric_param_ct; ++uii) {
                if (unlikely(score_cols[uii] == score_cols[numeric_param_ct])) {
                  logerrputs("Error: Identical --score[-list] column indexes.\n");
                  goto main_ret_INVALID_CMDLINE_A;
                }
              }
              ++numeric_param_ct;
            }
          }
          if (unlikely((pc.score_info.flags & (kfScoreHeaderIgnore | kfScoreHeaderRead)) == (kfScoreHeaderIgnore | kfScoreHeaderRead))) {
            logerrputs("Error: --score[-list] 'header' and 'header-read' modifiers cannot be used\ntogether.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          uint32_t model_flags_u = S_CAST(uint32_t, pc.score_info.flags & (kfScoreDominant | kfScoreRecessive | kfScoreCenter | kfScoreVarianceStandardize));
          if (unlikely(model_flags_u & (model_flags_u - 1))) {
            logerrputs("Error: --score[-list] 'dominant', 'recessive', 'center', and\n'variance-standardize' modifiers are mutually exclusive.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          if (!(pc.score_info.flags & kfScoreColAll)) {
            pc.score_info.flags |= kfScoreColDefault;
          }
          if (numeric_param_ct) {
            pc.score_info.varid_col_p1 = score_cols[0];
          }
          if (numeric_param_ct > 1) {
            pc.score_info.allele_col_p1 = score_cols[1];
          } else {
            pc.score_info.allele_col_p1 = pc.score_info.varid_col_p1 + 1;
          }
          if (numeric_param_ct == 3) {
            // a bit artificial, but it works
            const uint32_t col_idx = score_cols[2];
            const uint32_t col_idx_blen = 1 + UintSlen(col_idx);
            char* new_buf;
            if (unlikely(pgl_malloc(col_idx_blen + 1, &new_buf))) {
              goto main_ret_NOMEM;
            }
            pc.score_info.input_col_idx_range_list.names = new_buf;
            pc.score_info.input_col_idx_range_list.name_max_blen = col_idx_blen;
            pc.score_info.input_col_idx_range_list.name_ct = 1;
            u32toa_x(col_idx, '\0', new_buf);
            new_buf[col_idx_blen] = '\0';
            pc.score_info.input_col_idx_range_list.starts_range = R_CAST(unsigned char*, &(new_buf[col_idx_blen]));
          }
          pc.command_flags1 |= kfCommand1Score;
          pc.dependency_flags |= kfFilterAllReq;
        } else if (strequal_k_unsafe(flagname_p2, "core-col-nums")) {
          if (unlikely(pc.score_info.input_col_idx_range_list.name_ct)) {
            logerrputs("Error: --score-col-nums cannot be used when three numeric arguments are\nprovided to --score.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
          reterr = ParseNameRanges(&(argvk[arg_idx]), errstr_append, param_ct, 1, '-', &pc.score_info.input_col_idx_range_list);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
          score_col_nums_present = 1;
        } else if (strequal_k_unsafe(flagname_p2, "core-partitions")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.score_info.partition_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (strequal_k_unsafe(flagname_p2, "score")) {
          // --sscore <filename>: user-supplied pre-computed score file.
          // Mutually exclusive with --pfile, --bfile, --score.
          // Requires --give-best-pgs and --pheno (checked post-parse).
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          reterr = AllocFname(argvk[arg_idx + 1], flagname_p, &pc.sscore_fname);
          if (unlikely(reterr)) {
            goto main_ret_1;
          }
        } else if (unlikely(!strequal_k_unsafe(flagname_p2, "ilent"))) {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 't':
        if (strequal_k_unsafe(flagname_p2, "hreads")) {
          if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
            goto main_ret_INVALID_CMDLINE_2A;
          }
          if (unlikely(ScanPosintDefcapx(argvk[arg_idx + 1], &pc.max_thread_ct))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --threads argument '%s'.\n", argvk[arg_idx + 1]);
            goto main_ret_INVALID_CMDLINE_WWA;
          }
          if (pc.max_thread_ct > kMaxThreads) {
            logprintfww("Note: Reducing --threads argument to %u.  (If this is not large enough,\nrecompile with a larger kMaxThreads setting.)\n", kMaxThreads);
            pc.max_thread_ct = kMaxThreads;
          } else if (known_procs == -1) {
            // trigger BLAS/LAPACK warning?
            known_procs = 0;
          }
        } else {
          goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        }
        break;

      case 'u':
        goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        break;

      case 'v':
        goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        break;

      case 'w':
        goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        break;

      case 'x':
        goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        break;

      case 'y':
        goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        break;

      case 'z':
        goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
        break;

      default:
        goto main_ret_INVALID_CMDLINE_UNRECOGNIZED;
      main_param_zero:
        if (unlikely(param_ct)) {
          snprintf(g_logbuf, kLogbufSize, "Error: You provided %u argument%s to --%s, which has no parameters.\n", param_ct, (param_ct == 1)? "" : "s", flagname_p);
          goto main_ret_INVALID_CMDLINE_WWA;
        }
      }
    } while ((++cur_flag_idx) < flag_ct);
    if (!outname_end) {
      outname_end = &(outname[sizeof(PROG_NAME_STR) - 1]);
    } else if (!allow_misleading_out_arg) {
      char* outname_last_dot = Memrchr(outname, '.', outname_end - outname);
      if (outname_last_dot) {
        char* outname_ext = &(outname_last_dot[1]);
        const uint32_t outname_ext_slen = outname_end - outname_ext;
        if (strequal_k(outname_ext, "bcf", outname_ext_slen) ||
            strequal_k(outname_ext, "bed", outname_ext_slen) ||
            strequal_k(outname_ext, "bgen", outname_ext_slen) ||
            strequal_k(outname_ext, "gz", outname_ext_slen) ||
            strequal_k(outname_ext, "ped", outname_ext_slen) ||
            strequal_k(outname_ext, "pgen", outname_ext_slen) ||
            strequal_k(outname_ext, "vcf", outname_ext_slen)) {
          logerrputs("Error: --out argument is supposed to be a filename *prefix*; the value you've\nprovided looks like it is intended to be an entire filename.\nIf the current --out argument is really what you want, use\n--allow-misleading-out-arg to suppress this error.\n");
          goto main_ret_INVALID_CMDLINE_A;
        }
      }
    }

    pc.dependency_flags |= pc.filter_flags;

    // --pgs-model / --pgs-models-list / --bim validation.
    //
    // Every guard below is about the harmonisation MODE, not about which
    // of the two flags spelled it, so they test this one name.  Giving
    // the list route its own copy is precisely how the --pfile and
    // --sscore phenotype checks drifted apart, and shipped a wrong answer.
    if (unlikely(pc.pgs_model_fname && pc.pgs_models_list_fname)) {
      logerrputs("Error: --pgs-model cannot be used with --pgs-models-list.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    {
      const char* harmonise_fname =
        pc.pgs_model_fname ? pc.pgs_model_fname : pc.pgs_models_list_fname;
      const char* harmonise_flag =
        pc.pgs_model_fname ? "--pgs-model" : "--pgs-models-list";

    if (harmonise_fname || pc.bim_fname) {
      // The two flags require each other.
      if (unlikely(!harmonise_fname)) {
        logerrputs("Error: --bim requires --pgs-model or --pgs-models-list.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
      if (unlikely(!pc.bim_fname)) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s requires --bim.\n", harmonise_flag);
        goto main_ret_INVALID_CMDLINE_WWA;
      }
      // Mutually exclusive with genotype input and scoring commands.
      if (unlikely(load_params)) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s cannot be used with --pfile or --bfile.\n", harmonise_flag);
        goto main_ret_INVALID_CMDLINE_WWA;
      }
      if (unlikely(pc.score_info.input_fname)) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s cannot be used with --score.\n", harmonise_flag);
        goto main_ret_INVALID_CMDLINE_WWA;
      }
      if (unlikely(pc.sscore_fname)) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s cannot be used with --sscore.\n", harmonise_flag);
        goto main_ret_INVALID_CMDLINE_WWA;
      }
      if (unlikely(pc.score_info.flags & kfScoreGiveBestPgs)) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s cannot be used with --give-best-pgs.\n", harmonise_flag);
        goto main_ret_INVALID_CMDLINE_WWA;
      }
      if (unlikely(pc.score_info.flags & kfScoreGiveCorrMatrix)) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s cannot be used with --give-corr-matrix.\n", harmonise_flag);
        goto main_ret_INVALID_CMDLINE_WWA;
      }
    }
    if (unlikely(pc.pgs_model_allow_ambiguous && !harmonise_fname)) {
      logerrputs("Error: --allow-ambiguous requires --pgs-model or --pgs-models-list.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (pc.pgs_model_resolve_ambig) {
      if (unlikely(!harmonise_fname)) {
        logerrputs("Error: --resolve-ambiguous requires --pgs-model or --pgs-models-list.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
      // --allow-ambiguous keeps palindromes blindly; --resolve-ambiguous
      // resolves them by frequency and drops what it cannot.  Mixing
      // blindly-kept and frequency-resolved palindromes in one output
      // would hide a real difference in evidence.
      if (unlikely(pc.pgs_model_allow_ambiguous)) {
        logerrputs("Error: --resolve-ambiguous cannot be used with --allow-ambiguous.  --allow-ambiguous\nkeeps ambiguous variants without resolving their strand; --resolve-ambiguous\nresolves them from allele frequencies and drops the ones it cannot.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
      if (unlikely(!pc.read_freq_fname)) {
        logerrputs("Error: --resolve-ambiguous requires --read-freq (target allele frequencies).\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
      if (unlikely(!pc.pgs_model_col_freq)) {
        // Name the flag the user actually typed: both harmonisation routes
        // take the same sixth index, and "--pgs-model" would be wrong advice
        // for someone running --pgs-models-list.
        //
        // Built here and emitted with logerrputs rather than through the
        // word-wrapping path: wrapping re-flows the text at the terminal
        // width, which splits "column index" across a line break and makes
        // the message unmatchable by anything reading it back.
        char sixth_buf[256];
        snprintf(sixth_buf, sizeof(sixth_buf),
                 "Error: --resolve-ambiguous requires the optional sixth %s column index,\nthe effect-allele frequency in the source study.\n",
                 harmonise_flag);
        logerrputs(sixth_buf);
        goto main_ret_INVALID_CMDLINE_A;
      }
    } else if (harmonise_fname && pc.read_freq_fname) {
      // --resolve-ambiguous is the only thing in the harmonisation path that
      // reads the file.  Without it the path checks that the name exists
      // and then never opens it, so a wrong -- or entirely non-.afreq --
      // file produced output identical to a run without the flag.
      //
      // Scoped to the harmonisation mode in the wording as well as in the
      // guard: this branch cannot be reached alongside --score, which does
      // read the frequencies (to mean-impute missing genotypes), and an
      // unqualified "--read-freq is ignored" invites exactly the wrong
      // generalisation.
      logerrputs("Warning: in harmonisation mode --read-freq is ignored without --resolve-ambiguous;\nnothing else in the harmonisation path reads allele frequencies.\n");
    }
    }  // harmonise_fname scope

    // --sscore validation: mutually exclusive with --pfile/--bfile/--score.
    if (pc.sscore_fname) {
      if (unlikely(load_params)) {
        logerrputs("Error: --sscore cannot be used with --pfile or --bfile.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
      if (unlikely(pc.score_info.input_fname)) {
        logerrputs("Error: --sscore cannot be used with --score.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
      if (unlikely(!(pc.score_info.flags & (kfScoreGiveBestPgs | kfScoreGiveCorrMatrix)))) {
        logerrputs("Error: --sscore requires --give-best-pgs or --give-corr-matrix.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
      // A correlation matrix is a property of the scores alone, so --pheno is
      // required only by the part of the run that fits models.
      if (unlikely((pc.score_info.flags & kfScoreGiveBestPgs) && (!pc.pheno_fname))) {
        logerrputs("Error: --sscore requires --pheno (there is no .psam to fall back on).\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
      // --sscore consumes a score matrix that is already computed; there are no
      // variants left to exclude, so accepting the flag would silently do
      // nothing.
      if (unlikely(pc.mhc_region != kMhcNone)) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s cannot be used with --sscore, which reads scores that have already been computed.\n", MhcFlagName(pc.mhc_region));
        goto main_ret_INVALID_CMDLINE_WWA;
      }
    }
    if ((pc.score_info.flags & (kfScoreGiveBestPgs | kfScoreGiveCorrMatrix)) && (!pc.score_info.input_fname) && (!pc.sscore_fname)) {
      logerrprintf("Error: %s must be used with --score or --sscore.\n",
                   (pc.score_info.flags & kfScoreGiveBestPgs)? "--give-best-pgs" : "--give-corr-matrix");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (pc.pgs_pheno_pairs_fname) {
      if (unlikely(!(pc.score_info.flags & kfScoreGiveBestPgs))) {
        logerrputs("Error: --pgs-pheno-pairs must be used with --give-best-pgs.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
      // Two sources of truth for one question is exactly how the --pfile and
      // --sscore phenotype checks drifted apart; the pairs file already says
      // which phenotypes matter.
      if (unlikely(pc.pheno_range_list.name_ct)) {
        logerrprintf("Error: --pgs-pheno-pairs cannot be used with %s; the pairs file already selects\nthe phenotypes.\n",
                     (pc.misc_flags & kfMiscPhenoColNums)? "--pheno-col-nums" : "--pheno-name");
        goto main_ret_INVALID_CMDLINE_A;
      }
    }
    if ((pc.pgs_prevalence >= 0.0) && (!(pc.score_info.flags & kfScoreGiveBestPgs))) {
      logerrputs("Error: --prevalence must be used with --give-best-pgs.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if ((pc.pgs_folds != 1) && (!(pc.score_info.flags & kfScoreGiveBestPgs))) {
      logerrputs("Error: --folds must be used with --give-best-pgs.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    // --score is the only command; without it there is nothing to do.  (This
    // guard was lost when the other commands were removed, which left the tool
    // loading the whole fileset and then exiting 0 in silence.)
    if (pc.score_info.partition_fname) {
      if (unlikely(!(pc.command_flags1 & kfCommand1Score))) {
        logerrputs("Error: --score-partitions must be used with --score.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
      // No --q-score-range guard here, though the two fill the same membership
      // bitarray: adjutant has no --q-score-range parser case at all, so
      // qsr_range_fname is permanently null and the check would be dead code.
      // (Three pre-existing reads of it are dead for the same reason.)
    }
    if (!pc.command_flags1 && !pc.sscore_fname && !pc.pgs_model_fname &&
        !pc.pgs_models_list_fname) {
      goto main_ret_NULL_CALC;
    }
    if (unlikely((xload & kfXloadMap) && (!(xload & (kfXloadPed | kfXloadPlink1Dosage))))) {
      logerrputs("Error: --map must be used with --import-dosage or --ped.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely((xload & kfXloadOxGen) && (!(xload & kfXloadOxSample)))) {
      // could permit .fam/.psam, but unless Oxford software supports that mode
      // it's pointless
      logerrputs("Error: --gen must be used with --sample or --data.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely((xload & kfXloadEigGeno) && ((xload & (kfXloadEigInd | kfXloadEigSnp)) != (kfXloadEigInd | kfXloadEigSnp)))) {
      logerrputs("Error: --eiggeno must be used with either --eigfile, or --eigind + --eigsnp.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    // may as well permit merge here
    if (unlikely((make_plink2_flags & (kfMakePlink2MMask | kfMakePlink2TrimAlts | kfMakePgenErasePhase | kfMakePgenEraseDosage)) && (pc.command_flags1 & (~(kfCommand1MakePlink2 | kfCommand1Pmerge))))) {
      logerrputs("Error: When the 'multiallelics=', 'trim-alts', and/or 'erase-...' modifier is\npresent, --make-bed/--make-[b]pgen cannot be combined with other commands.\n(Other filters are fine.)\n");
      goto main_ret_INVALID_CMDLINE;
    }
    if (make_plink2_flags & kfMakePlink2MMask) {
      if (unlikely(make_plink2_flags & kfMakePlink2TrimAlts)) {
        logerrputs("Error: --make-bed/--make-[b]pgen 'multiallelics=' cannot be used with\n'trim-alts'.\n");
        goto main_ret_INVALID_CMDLINE;
      }
    }
    if (pc.command_flags1 & (~(kfCommand1MakePlink2 | kfCommand1Pmerge))) {
      if (unlikely(make_plink2_flags & (kfMakePlink2SetInvalidHaploidMissing | kfMakePlink2SetMixedMtMissing | kfMakePlink2SetMeMissing | kfMakePlink2FillMissingWithRef))) {
        logerrputs("Errpr: --set-invalid-haploid-missing/--set-mixed-mt-missing/--set-me-missing/\n--fill-missing-with-ref must be used with --make-[b]pgen/--make-bed and no\nother non-merge commands.\n");
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely(pc.zero_cluster_fname)) {
        logerrputs("Error: --zero-cluster must be used with --make-[b]pgen/--make-bed and no other\nnon-merge commands.\n");
        goto main_ret_INVALID_CMDLINE;
      }
      if (pc.command_flags1 & (~(kfCommand1MakePlink2 | kfCommand1Exportf | kfCommand1Pmerge))) {
      }
    }
    if (pc.keep_cat_phenoname && (!pc.keep_cat_names_flattened) && (!pc.keep_cats_fname)) {
      logerrputs("Error: --keep-cat-pheno must be used with --keep-cats and/or --keep-cat-names.\n");
    }
    if (pc.remove_cat_phenoname && (!pc.remove_cat_names_flattened) && (!pc.remove_cats_fname)) {
      logerrputs("Error: --remove-cat-pheno must be used with --remove-cats and/or\n--remove-cat-names.\n");
    }
    if ((pc.hard_call_thresh == UINT32_MAX) && (xload & (kfXloadVcf | kfXloadBcf | kfXloadOxGen | kfXloadOxBgen))) {
      if (unlikely(pc.dosage_erase_thresh > (kDosageMid / 10))) {
        logerrputs("Error: --dosage-erase-threshold value cannot be larger than (default)\n--hard-call-threshold value.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
    } else {
      if (unlikely(pc.dosage_erase_thresh > pc.hard_call_thresh)) {
        logerrputs("Error: --dosage-erase-threshold value cannot be larger than\n--hard-call-threshold value.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
    }
    {
    }
    if (unlikely(!strcmp(pc.missing_catname, pc.output_missing_pheno))) {
      logerrputs("Error: --missing-catname and --output-missing-phenotype strings can't match.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely((pc.misc_flags & kfMiscChrOverrideCmdline) && (!chr_info.chrset_source))) {
      logerrputs("Error: --chr-override requires an explicit chromosome set.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely((xload & kfXloadPlink1Dosage) && (!(load_params & kfLoadParamsPsam)))) {
      logerrputs("Error: --import-dosage requires a .fam file.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely((pc.misc_flags & kfMiscCovarColNums) && (!pc.covar_fname) && (!pc.pheno_fname))) {
      logerrputs("Error: --covar-col-nums requires --covar or --pheno.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely(pc.score_info.qsr_range_fname && (!pc.score_info.input_fname))) {
      logerrputs("Error: --q-score-range cannot be used without --score[-list].\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely(score_col_nums_present && (!(pc.command_flags1 & kfCommand1Score)) && (!pc.sscore_fname))) {
      logerrputs("Error: --score-col-nums must be used with --score[-list] or --sscore.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (unlikely(r2_required && (!(pc.command_flags1 & kfCommand1Vcor)))) {
      logerrputs("Error: --ld-window.../--ld-snp... must be used with --r[2]-[un]phased.\n");
      goto main_ret_INVALID_CMDLINE_A;
    }
    if (pc.extract_col_cond_info.params) {
      if (unlikely((!pc.extract_col_cond_info.match_substr) && pc.extract_col_cond_info.match_flattened && pc.extract_col_cond_info.mismatch_flattened)) {
        logerrputs("Error: --extract-col-cond-match and --extract-col-cond-mismatch can only be\nused together when --extract-col-cond-substr is specified.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
      if (unlikely(pc.extract_col_cond_info.max < pc.extract_col_cond_info.min)) {
        logerrputs("Error: --extract-col-cond-max value can't be smaller than\n--extract-col-cond-min value.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
    }
    if (explicit_perm_count_modifier || implicit_perm_count_omission) {
      // at least one permutation test requested
      if (explicit_perm_count_modifier) {
      }
    }
    // Upstream arbitrated between seven variant inclusion filters here, since
    // it is not obvious whether the union or the intersection is meant.  Only
    // --extract survives in adjutant, so there is nothing left to arbitrate.
    if (unlikely(pmerge_required && (!(pc.command_flags1 & kfCommand1Pmerge)))) {
      if (delete_pmerge_result) {
        logerrputs("Error: --delete-pmerge-result must be used with --pmerge[-list].\n");
        logerrputs("Error: --merge-cm-mode must be used with --pmerge[-list].\n");
        logerrputs("Error: --merge-filter-mode must be used with --pmerge[-list].\n");
        logerrputs("Error: --merge-info-mode must be used with --pmerge[-list].\n");
        logerrputs("Error: --merge-info-sort must be used with --pmerge[-list].\n");
        logerrputs("Error: --merge-max-allele-ct must be used with --pmerge[-list].\n");
        logerrputs("Error: --merge-mode must be used with --pmerge[-list].\n");
        logerrputs("Error: --merge-parents-mode must be used with --pmerge[-list].\n");
        logerrputs("Error: --merge-pheno-mode must be used with --pmerge[-list].\n");
        logerrputs("Error: --merge-pheno-sort must be used with --pmerge[-list].\n");
        logerrputs("Error: --merge-qual-mode must be used with --pmerge[-list].\n");
        logerrputs("Error: --merge-sids must be used with --pmerge[-list].\n");
        logerrputs("Error: --merge-xheader-mode must be used with --pmerge[-list].\n");
        logerrputs("Error: --multiallelics-already-joined must be used with --pmerge[-list].\n");
        logerrputs("Error: --pheno-inner-join must be used with --pmerge[-list].\n");
      } else {
        assert(0);
      }
      goto main_ret_INVALID_CMDLINE;
    }
    if (delete_pmerge_result) {
      if (unlikely(pc.command_flags1 == kfCommand1Pmerge)) {
        logerrputs("Error: --delete-pmerge-result doesn't make sense when --pmerge[-list] is run\nwithout any other commands.\n");
        goto main_ret_INVALID_CMDLINE_A;
      }
    }

#ifdef USE_MKL
    if (!mkl_native) {
#  ifdef USE_SSE42
#    ifdef USE_AVX2
      int status = mkl_cbwr_set(MKL_CBWR_AVX2);
#    else
      int status = mkl_cbwr_set(MKL_CBWR_SSE4_2);
#    endif
      if (status != MKL_CBWR_SUCCESS) {
        // could happen for AMD processors.  This is the least-bad option?
        mkl_cbwr_set(MKL_CBWR_COMPATIBLE);
      }
#  else
      mkl_cbwr_set(MKL_CBWR_COMPATIBLE);
#  endif
    }
#endif
    sfmt_t main_sfmt;
    if (!rseeds) {
      uint32_t seed = S_CAST(uint32_t, time(nullptr));
      snprintf(g_logbuf, kLogbufSize, "Random number seed: %u\n", seed);
      logputs_silent(g_logbuf);
      sfmt_init_gen_rand(&main_sfmt, seed);
    } else {
      if (rseed_ct == 1) {
        sfmt_init_gen_rand(&main_sfmt, rseeds[0]);
      } else {
        sfmt_init_by_array(&main_sfmt, rseeds, rseed_ct);
      }
      free(rseeds);
      rseeds = nullptr;
    }

    if (unlikely(CmdlineParsePhase3(0, malloc_size_mib, memory_require, &pcm, &bigstack_ua))) {
      goto main_ret_NOMEM;
    }

    if (pc.max_thread_ct > 8) {
      logprintf("Using up to %u threads (change this with --threads).\n", pc.max_thread_ct);
    } else {
      // "x compute threads" instead of "x threads" since e.g. when
      // max_thread_ct == 2, some code will use one I/O thread and two
      // compute threads.  Not worth the trouble of writing special-case code
      // to avoid that: with 2 cores, the I/O thread usually isn't busy enough
      // to justify only 1 compute thread.
      logprintf("Using %s%u compute thread%s.\n", (pc.max_thread_ct > 1)? "up to " : "", pc.max_thread_ct, (pc.max_thread_ct == 1)? "" : "s");
    }
    print_end_time = 1;


    // --pgs-model standalone mode: harmonise a PGS weight file against a
    // target .bim.  No genotypes loaded.  Bypasses Plink2Core entirely.
    if (pc.pgs_model_fname) {
      reterr = RunPgsModel(pc.bim_fname, pc.pgs_model_fname,
                           pc.pgs_model_by_id,
                           pc.pgs_model_col_id,
                           pc.pgs_model_col_chr,
                           pc.pgs_model_col_pos,
                           pc.pgs_model_col_effect,
                           pc.pgs_model_col_other,
                           pc.pgs_model_col_weight,
                           pc.pgs_model_col_freq,
                           pc.pgs_model_allow_ambiguous,
                           pc.pgs_model_resolve_ambig,
                           pc.pgs_model_maf_threshold,
                           pc.read_freq_fname,
                           pc.mhc_region,
                           outname, outname_end);
      goto main_ret_1;
    }

    // --pgs-models-list: same harmonisation, several models, merged into
    // one matrix.
    if (pc.pgs_models_list_fname) {
      reterr = RunPgsModelsList(pc.bim_fname, pc.pgs_models_list_fname,
                                pc.pgs_model_by_id,
                                pc.pgs_model_col_id,
                                pc.pgs_model_col_chr,
                                pc.pgs_model_col_pos,
                                pc.pgs_model_col_effect,
                                pc.pgs_model_col_other,
                                pc.pgs_model_col_weight,
                                pc.pgs_model_col_freq,
                                pc.pgs_model_allow_ambiguous,
                                pc.pgs_model_resolve_ambig,
                                pc.pgs_model_maf_threshold,
                                pc.read_freq_fname,
                                pc.mhc_region,
                                outname, outname_end);
      goto main_ret_1;
    }

    // --sscore standalone mode: evaluate a pre-computed score matrix without
    // loading any genotype data.  Bypasses Plink2Core entirely.
    if (pc.sscore_fname) {
      BLAS_SET_NUM_THREADS(1);
      // Without --give-best-pgs there is nothing to fit, so no phenotype file
      // is passed even if one was supplied: --give-corr-matrix must give the
      // same matrix either way, and the surest way to guarantee that is for
      // the phenotype not to reach the computation at all.
      reterr = BestPgsFromSscore(pc.sscore_fname,
                                 pc.score_info.input_col_idx_range_list.name_ct ?
                                     &pc.score_info.input_col_idx_range_list : nullptr,
                                 (pc.score_info.flags & kfScoreGiveBestPgs)?
                                     pc.pheno_fname : nullptr,
                                 pc.pheno_range_list.name_ct ?
                                     &pc.pheno_range_list : nullptr,
                                 pc.covar_fname,
                                 pc.pgs_prevalence, pc.pgs_folds,
                                 pc.pgs_fold_seed,
                                 (pc.score_info.flags / kfScoreGiveCorrMatrix) & 1,
                                 (pc.score_info.flags / kfScoreBestPgsNoNa) & 1,
                                 pc.pgs_pheno_pairs_fname,
                                 outname, outname_end);
      goto main_ret_1;
    }

    // nonstandard cases (CNV, etc.) here
    if ((pc.command_flags1 == kfCommand1PgenInfo) && (load_params != kfLoadParamsPfileAll) && (!xload)) {
      // special case: don't require .psam/.pvar file
      if (unlikely(!(load_params & kfLoadParamsPgen))) {
        logerrputs("Error: --pgen-info requires a .pgen file.\n");
        goto main_ret_INVALID_CMDLINE;
      }
      reterr = PgenInfoStandalone(pgenname, pc.pginame);
    } else {
      // --real-ref-alleles is an exception since it applies to merge.
      if (unlikely(pc.dependency_flags && (!(pc.command_flags1 & (~kfCommand1Pmerge))))) {
        logerrputs("Error: Basic file conversions do not support regular filter or transform\noperations.  Rerun your command with --make-bed/--make-[b]pgen.\n");
        goto main_ret_INVALID_CMDLINE;
      }
      if ((xload && pc.command_flags1 && (import_flags & kfImportKeepAutoconv)) || ((pc.command_flags1 & kfCommand1Pmerge) & (pc.command_flags1 & (~kfCommand1Pmerge)))) {
        // pfile-affecting input/output settings must be consistent since we
        // need to reload what we write.
        if (unlikely(pc.misc_flags & kfMiscAffection01)) {
          logerrputs("Error: --1 cannot be used with --keep-autoconv or non-standalone\n--pmerge[-list].\n");
          goto main_ret_INVALID_CMDLINE_A;
        }
      }
      if (xload) {
        char* convname_end = outname_end;
        if (pc.command_flags1) {
          if (!(import_flags & kfImportKeepAutoconv)) {
            convname_end = Stpcpy(convname_end, "-temporary");
          }
        } else {
          import_flags |= kfImportKeepAutoconv;
        }
        const uint32_t convname_slen = convname_end - outname;
        const LoadFilterLogFlags load_filter_log_import_flags = pc.load_filter_log_flags & kfLoadFilterLogImportMask;
        uint32_t pgen_generated = 1;
        uint32_t pgi_generated = 0;
        uint32_t psam_generated = 1;

        // Compress by default for VCF/BCF due to potentially large INFO
        // section; otherwise only do it in "--keep-autoconv vzs" case.
        uint32_t pvar_is_compressed;
        if (xload & (kfXloadVcf | kfXloadBcf)) {
          const uint32_t no_samples_ok = !((pc.dependency_flags & (kfFilterAllReq | kfFilterPsamReq)) || (pc.command_flags1 & kfCommand1Pmerge));
          const uint32_t is_vcf = (xload / kfXloadVcf) & 1;
          pvar_is_compressed = ((import_flags & (kfImportKeepAutoconv | kfImportKeepAutoconvVzs)) != kfImportKeepAutoconv);
          if (no_samples_ok && is_vcf && (!(import_flags & (kfImportKeepAutoconv | kfImportVcfRefNMissing))) && pc.command_flags1) {
            // special case: just treat the VCF as a .pvar file
            strcpy(pvarname, pgenname);
            pgenname[0] = '\0';
            goto main_reinterpret_vcf_instead_of_converting;
          }
          // Default to compression level 1 for temporary .pvar files for now.
          //
          // Level 1 may actually be best in a much wider variety of scenarios
          // as of this writing, but I won't try to tune any other compression
          // defaults for now.  Interestingly, the current setup is actually a
          // bit backwards given observed zstd behavior: the long INFO fields
          // motivating automatic compression here are best handled with level
          // 3, while the other simpler text files generated by plink2 are
          // likely to compress *better*, not just faster, with level 1 for
          // some reason.
          const uint32_t zst_level = g_zst_level;
          if (!(import_flags & kfImportKeepAutoconv)) {
            g_zst_level = 1;
          }
          if (is_vcf) {
          } else {
          }
          g_zst_level = zst_level;
        } else {
          pvar_is_compressed = (import_flags / kfImportKeepAutoconvVzs) & 1;
          // We have previously verified that --1 and either --keep-autoconv
          // or non-standalone --pmerge[-list] haven't been specified
          // simultaneously.  Therefore:
          // - If --1 is specified, and we're executing either --pmerge[-list]
          //   or at least one regular plink2 command, the temporary .psam file
          //   should use 0/1 coding for binary phenotypes.
          // - However, if --1 is specified and we aren't executing any regular
          //   plink2 commands, the output .psam is not temporary and should
          //   use 1/2 coding.
          if (xload & kfXloadPed) {
          } else if (xload & kfXloadTped) {
          } else if (xload & kfXloadOxGen) {
          } else if (xload & kfXloadOxBgen) {
          } else if (xload & kfXloadOxHaps) {
          } else if (xload & kfXloadEigGeno) {
          } else if (xload & kfXloadPlink1Dosage) {
          } else if (likely(xload & kfXloadGenDummy)) {
          } else {
            // We should have errored out with a better message before this
            // point.
            assert(0);
            logerrputs("Error: Invalid input flag combination.\n");
            goto main_ret_INVALID_CMDLINE;
          }
        }
        if (reterr || (!pc.command_flags1)) {
          goto main_ret_1;
        }

        if (load_filter_log_import_flags) {
          pc.load_filter_log_flags |= kfLoadFilterLogImportMergeAlreadyApplied;
        }

        pc.hard_call_thresh = UINT32_MAX;
        pc.dosage_erase_thresh = 0;

        if (pgen_generated) {
          snprintf(memcpya(pgenname, outname, convname_slen), kMaxOutfnameExtBlen - 10, ".pgen");
        }
        snprintf(memcpya(pvarname, outname, convname_slen), kMaxOutfnameExtBlen - 10, pvar_is_compressed? ".pvar.zst" : ".pvar");
        if (psam_generated) {
          snprintf(memcpya(psamname, outname, convname_slen), kMaxOutfnameExtBlen - 10, ".psam");
        }
        if (!(import_flags & kfImportKeepAutoconv)) {
          if (pgen_generated) {
            if (unlikely(PushLlStr(pgenname, &file_delete_list))) {
              goto main_ret_NOMEM;
            }
            if (pgi_generated) {
              char* pgenname_end = &(pgenname[convname_slen + strlen(".pgen")]);
              strcpy_k(pgenname_end, ".pgi");
              if (unlikely(PushLlStr(pgenname, &file_delete_list))) {
                goto main_ret_NOMEM;
              }
              pgenname_end[0] = '\0';
            }
          }
          if (unlikely(PushLlStr(pvarname, &file_delete_list))) {
            goto main_ret_NOMEM;
          }
          if (psam_generated) {
            if (unlikely(PushLlStr(psamname, &file_delete_list))) {
              goto main_ret_NOMEM;
            }
          }
        }
        *outname_end = '\0';
      }
    main_reinterpret_vcf_instead_of_converting:
      if ((pc.dependency_flags & kfFilterOpportunisticPgen) && (pgenname[0] != '\0')) {
        pc.dependency_flags |= kfFilterAllReq;
      }
      if ((pc.dependency_flags & kfFilterAllReq) || (pc.command_flags1 & kfCommand1Pmerge)) {
        if (unlikely((!xload) && (load_params != kfLoadParamsPfileAll))) {
          logerrputs("Error: A full fileset (.pgen/.bed + .pvar/.bim + .psam/.fam) is required for\nthis.\n");
          goto main_ret_INVALID_CMDLINE_A;
        }
      } else {
        // no genotype file required
        if (!(pc.dependency_flags & kfFilterOpportunisticPgen)) {
          pgenname[0] = '\0';
        }

        if (pc.dependency_flags & kfFilterPvarReq) {
          if (unlikely((!xload) && (!(load_params & kfLoadParamsPvar)))) {
            logerrputs("Error: A .pvar/.bim file is required for this.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
        } else {
          pvarname[0] = '\0';
        }
        if (pc.dependency_flags & kfFilterPsamReq) {
          if (unlikely((!xload) && (!(load_params & kfLoadParamsPsam)))) {
            logerrputs("Error: A .psam/.fam file is required for this.\n");
            goto main_ret_INVALID_CMDLINE_A;
          }
        } else {
          psamname[0] = '\0';
        }
      }

      if (pc.command_flags1 & kfCommand1Pmerge) {
        const LoadFilterLogFlags load_filter_log_merge_flags = pc.load_filter_log_flags & kfLoadFilterLogMergeMask;
        char* merge_outname_end = outname_end;
        if (make_plink2_flags & (kfMakePgen | kfMakePvar | kfMakePsam)) {
          merge_outname_end = strcpya_k(merge_outname_end, "-merge");
        }
        if (unlikely(reterr)) {
          goto main_ret_1;
        }
        pc.command_flags1 ^= kfCommand1Pmerge;
        if (!pc.command_flags1) {
          goto main_ret_1;
        }
        if (delete_pmerge_result) {
          if (unlikely(PushLlStr(pgenname, &file_delete_list) ||
                       PushLlStr(pvarname, &file_delete_list) ||
                       PushLlStr(psamname, &file_delete_list))) {
            goto main_ret_NOMEM;
          }
        }
        // Don't need to repeat --set-{all,missing}-var-ids.
        if (pc.varid_template_str) {
          free(pc.varid_template_str);
          pc.varid_template_str = nullptr;
        }
        if (load_filter_log_merge_flags) {
          pc.load_filter_log_flags |= kfLoadFilterLogImportMergeAlreadyApplied;
        }
        // --real-ref-alleles communicates to --pmerge[-list] that all input
        // .bed filesets have real REF alleles.
        // However, we then need to clear this flag since Plink2Core() would
        // otherwise error out (since --real-ref-alleles does not make sense on
        // the merged .pgen).
        pc.misc_flags &= ~kfMiscRealRefAlleles;
      }


      BLAS_SET_NUM_THREADS(1);
      reterr = Plink2Core(&pc, make_plink2_flags, pgenname, psamname, pvarname, outname, outname_end, &chr_info, &main_sfmt);
    }
  }
  while (0) {
  main_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  main_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  main_ret_INVALID_CMDLINE_UNRECOGNIZED:
    InvalidArg(argv[arg_idx]);
    logerrputsb();
    logerrputs(errstr_append);
    reterr = kPglRetInvalidCmdline;
    break;
  main_ret_INVALID_CMDLINE_INPUT_CONFLICT:
    logerrprintf("Error: --%s conflicts with another input flag.\n%s", flagname_p, errstr_append);
    reterr = kPglRetInvalidCmdline;
    break;
  main_ret_INVALID_CMDLINE_WWA:
    WordWrapB(0);
  main_ret_INVALID_CMDLINE_2A:
    logerrputsb();
  main_ret_INVALID_CMDLINE_A:
    logerrputs(errstr_append);
  main_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  main_ret_NULL_CALC:
    if (pc.dependency_flags) {
      logerrputs("Warning: No output requested.  (Did you forget --make-bed/--make-[b]pgen?)\nExiting.\n");
    } else {
      logerrputs("Warning: No output requested.  Exiting.\n");
    }
  main_ret_NULL_CALC_0:
    fputs(kCmdlineFormatStr, stdout);
    fputs(notestr_null_calc2, stdout);
    reterr = kPglRetSkipped;
    break;
  }
 main_ret_1:
  if (reterr == kPglRetNomemCustomMsg) {
    if (g_failed_alloc_attempt_size) {
      logerrprintf("Failed allocation size: %" PRIu64 "\n", g_failed_alloc_attempt_size);
    }
    reterr = kPglRetNomem;
  } else {
    DispExitMsg(reterr);
  }
  while (0) {
  main_ret_NOMEM_NOLOG:
    PrintVer();
  main_ret_NOMEM_NOLOG2:
    fputs(kErrstrNomem, stderr);
    if (g_failed_alloc_attempt_size) {
      fprintf(stderr, "Failed allocation size: %" PRIu64 "\n", g_failed_alloc_attempt_size);
    }
    reterr = kPglRetNomem;
    break;
  }
 main_ret_NOLOG:
  free_cond(vcf_dosage_import_field);
  free_cond(ox_missing_code);
  free_cond(import_single_chr_str);
  free_cond(const_fid);
  free_cond(rseeds);
  CleanupPlink2CmdlineMeta(&pcm);
  free_cond(pc.zero_cluster_phenoname);
  free_cond(pc.zero_cluster_fname);
  free_cond(pc.rename_chrs_fname);
  free_cond(pc.indv_str);
  free_cond(pc.not_covar_flattened);
  free_cond(pc.not_pheno_flattened);
  free_cond(pc.indep_preferred_fname);
  free_cond(pc.vscore_fname);
  free_cond(pc.sscore_fname);
  free_cond(pc.bim_fname);
  free_cond(pc.pgs_model_fname);
  free_cond(pc.pgs_models_list_fname);
  free_cond(pc.pgs_pheno_pairs_fname);
  free_cond(pc.recover_var_ids_fname);
  free_cond(pc.update_parental_ids_fname);
  free_cond(pc.update_sample_ids_fname);
  free_cond(pc.update_name_flag);
  free_cond(pc.update_map_flag);
  free_cond(pc.update_chr_flag);
  free_cond(pc.alt_allele_flag);
  free_cond(pc.ref_allele_flag);
  free_cond(pc.keep_col_match_name);
  free_cond(pc.keep_col_match_flattened);
  free_cond(pc.keep_col_match_fname);
  free_cond(pc.require_no_info_flattened);
  free_cond(pc.require_info_flattened);
  free_cond(pc.fa_fname);
  free_cond(pc.loop_cats_phenoname);
  free_cond(pc.covar_quantnorm_flattened);
  free_cond(pc.quantnorm_flattened);
  free_cond(pc.vstd_flattened);
  free_cond(pc.require_covar_flattened);
  free_cond(pc.require_pheno_flattened);
  free_cond(pc.split_cat_phenonames_flattened);
  free_cond(pc.remove_cat_phenoname);
  free_cond(pc.remove_cat_names_flattened);
  free_cond(pc.remove_cats_fname);
  free_cond(pc.keep_cat_phenoname);
  free_cond(pc.keep_cat_names_flattened);
  free_cond(pc.keep_cats_fname);
  free_cond(pc.family_missing_catname);
  free_cond(pc.catpheno_name);
  free_cond(pc.within_fname);
  free_cond(pc.read_freq_fname);
  free_cond(pc.glm_local_covar_fname);
  free_cond(pc.glm_local_pvar_fname);
  free_cond(pc.glm_local_psam_fname);
  free_cond(pc.freq_alt1_binstr);
  free_cond(pc.freq_ref_binstr);
  free_cond(pc.removefam_fnames);
  free_cond(pc.remove_fnames);
  free_cond(pc.keepfam_fnames);
  free_cond(pc.keep_fnames);
  free_cond(pc.exclude_fnames);
  free_cond(pc.extract_intersect_fnames);
  free_cond(pc.extract_fnames);
  free_cond(pc.sample_sort_fname);
  free_cond(pc.covar_fname);
  free_cond(pc.pheno_fname);
  free_cond(pc.varid_exclude_snp);
  free_cond(pc.varid_snp);
  free_cond(pc.varid_to);
  free_cond(pc.varid_from);
  free_cond(pc.missing_varid_match);
  free_cond(pc.varid_multi_nonsnp_template_str);
  free_cond(pc.varid_multi_template_str);
  free_cond(pc.varid_template_str);
  free_cond(pc.var_filter_exceptions_flattened);
  free_cond(pc.pginame);
  if (file_delete_list) {
    do {
      LlStr* llstr_ptr = file_delete_list->next;
      if (unlikely(unlink(file_delete_list->str))) {
        logerrprintfww("Error: Failed to delete %s : %s.\n", file_delete_list->str, strerror(errno));
        if (reterr == kPglRetSuccess) {
          reterr = kPglRetWriteFail;
        }
      }
      free(file_delete_list);
      file_delete_list = llstr_ptr;
    } while (file_delete_list);
  }
  CleanupChrInfo(&chr_info);
  CleanupExtractColCond(&pc.extract_col_cond_info);
  CleanupCmpExpr(&pc.exclude_if_info_expr);
  CleanupCmpExpr(&pc.extract_if_info_expr);
  CleanupCmpExpr(&pc.remove_if_expr);
  CleanupCmpExpr(&pc.keep_if_expr);
  CleanupScore(&pc.score_info);
  CleanupRangeList(&pc.covar_range_list);
  CleanupRangeList(&pc.pheno_range_list);
  CleanupRangeList(&pc.exclude_snps_range_list);
  CleanupRangeList(&pc.snps_range_list);
  if (warning_errcode && g_stderr_written_to && (!reterr)) {
    logerrputs("--warning-errcode: One or more warnings in this run; exiting with code 61.\n");
    reterr = kPglRetWarningErrcode;
  }
  if (CleanupLogfile(print_end_time) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  free_cond(bigstack_ua);
  return S_CAST(int32_t, reterr);
}
