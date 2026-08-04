#ifndef __ADJUTANT_FILTER_H__
#define __ADJUTANT_FILTER_H__

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
#include "adjutant_cmdline.h"
#include "adjutant_common.h"

#ifdef __cplusplus
namespace plink2 {
#endif

typedef struct ExtractColCondStruct {
  NONCOPYABLE(ExtractColCondStruct);
  TwoColParams* params;
  char* match_flattened;
  char* mismatch_flattened;
  uint32_t match_substr;
  double min;
  double max;
} ExtractColCondInfo;

void InitExtractColCond(ExtractColCondInfo* efip);

void CleanupExtractColCond(ExtractColCondInfo* efip);

PglErr FromToFlag(const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const char* varid_from, const char* varid_to, uint32_t raw_variant_ct, uint32_t max_variant_id_slen, uintptr_t variant_id_htable_size, uintptr_t* variant_include, ChrInfo* cip, uint32_t* variant_ct_ptr);




PglErr TokenExtractExclude(const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const char* fnames, const char* flagname, uint32_t raw_variant_ct, uint32_t max_variant_id_slen, uintptr_t variant_id_htable_size, VfilterType vft, uint32_t max_thread_ct, uintptr_t* variant_include, uint32_t* variant_ct_ptr);
FLAGSET_DEF_START()
  kfKeep0,
  kfKeepRemove = (1 << 0),
  kfKeepFam = (1 << 1)
FLAGSET_DEF_END(KeepFlags);

PglErr KeepOrRemove(const char* fnames, const SampleIdInfo* siip, uint32_t raw_sample_ct, KeepFlags flags, uintptr_t* sample_include, uint32_t* sample_ct_ptr);

void KeepOneId(const char* sample_id_flattened, const SampleIdInfo* siip, uint32_t raw_sample_ct, uint32_t iid_sid, uintptr_t* sample_include, uint32_t* sample_ct_ptr);


PglErr RequirePheno(const PhenoCol* pheno_cols, const char* pheno_names, const char* require_pheno_flattened, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t is_covar, uintptr_t* sample_include, uint32_t* sample_ct_ptr);



// this mode is 0 in SID-only case
ENUM_U31_DEF_START()
  kSelectSidMissingness0,
  kSelectSidMissingnessHhMissing,
  kSelectSidMissingnessHc,
  kSelectSidMissingnessDosage
ENUM_U31_DEF_END(SelectSidMissingnessMode);

ENUM_U31_DEF_START()
  kSelectSidTiebreak0,
  kSelectSidTiebreakFirst,
  kSelectSidTiebreakFirstAscii,
  kSelectSidTiebreakLast,
  kSelectSidTiebreakLastAscii
ENUM_U31_DEF_END(SelectSidTiebreakMode);

void ComputeAlleleFreqs(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const uint64_t* founder_allele_ddosages, uint32_t variant_ct, double af_pseudocount, double* allele_freqs);

PglErr ReadAlleleFreqs(const uintptr_t* variant_include, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const char* read_freq_fname, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, double af_pseudocount, uint32_t max_thread_ct, double* allele_freqs, uintptr_t** variant_afreqcalcp);



PglErr SelectSidRepresentatives(const uintptr_t* sex_nm, const uintptr_t* sex_male, const uint32_t* sample_missing_cts, const uint32_t* sample_hethap_cts, const PedigreeIdInfo* piip, SelectSidMissingnessMode missingness_mode, SelectSidTiebreakMode tiebreak_mode, uint32_t parents_only, uint32_t raw_sample_ct, uintptr_t* sample_include, uint32_t* sample_ct_ptr);



// main loop has a bitwise Nref-or-Alt1 check
ENUM_U31_DEF_START()
  kFreqFilterNonmajor,
  kFreqFilterNref,
  kFreqFilterAlt1,
  kFreqFilterMinor = 4
ENUM_U31_DEF_END(FreqFilterMode);
PglErr MakeFounders(const uintptr_t* sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t require_two, PedigreeIdInfo* piip, uintptr_t* founder_info);



#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __ADJUTANT_FILTER_H__
