#ifndef __ADJUTANT_MATRIX_CALC_H__
#define __ADJUTANT_MATRIX_CALC_H__

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
// --make-rel, --make-grm-list, --make-grm-bin, --make-grm-sparse
//   --make-rel iff kfGrmMatrixShapemask nonzero
//   --make-grm-list iff kfGrmList
//   --make-grm-bin iff kfGrmBin
//   --make-grm-sparse iff kfGrmSparse
//   (no GRM output possible in --pca case)
FLAGSET_DEF_START()
  kfScore0,
  kfScoreHeaderIgnore = (1 << 0),
  kfScoreHeaderRead = (1 << 1),
  kfScoreNoMeanimpute = (1 << 2),
  kfScoreDominant = (1 << 3),
  kfScoreRecessive = (1 << 4),
  kfScoreCenter = (1 << 5),
  kfScoreVarianceStandardize = (1 << 6),
  kfScoreSe = (1 << 7),
  kfScoreZs = (1 << 8),
  kfScoreIgnoreDupIds = (1 << 9),
  kfScoreListVariants = (1 << 10),
  kfScoreListVariantsZs = (1 << 11),
  kfScoreMultiInput = (1 << 12),

  kfScoreQsrHeader = (1 << 13),
  kfScoreQsrMin = (1 << 14),

  kfScoreColMaybefid = (1 << 15),
  kfScoreColFid = (1 << 16),
  kfScoreColMaybesid = (1 << 17),
  kfScoreColSid = (1 << 18),
  kfScoreColPheno1 = (1 << 19),
  kfScoreColPhenos = (1 << 20),
  kfScoreColNallele = (1 << 21),
  kfScoreColDenom = (1 << 22),
  kfScoreColDosageSum = (1 << 23),
  kfScoreColScoreAvgs = (1 << 24),
  kfScoreColScoreSums = (1 << 25),
  kfScoreGiveBestPgs = (1 << 26),
  kfScoreGiveCorrMatrix = (1 << 27),
  // --give-best-pgs no-na: omit the .best_pgs columns the run's configuration
  // cannot fill.  Must stay clear of kfScoreColAll above, which is derived from
  // kfScoreColScoreSums and so covers only bits 16..25.
  kfScoreBestPgsNoNa = (1 << 28),
  kfScoreListColDefault = (kfScoreColMaybefid | kfScoreColMaybesid | kfScoreColPhenos | kfScoreColScoreAvgs),
  kfScoreColDefault = (kfScoreListColDefault | kfScoreColNallele | kfScoreColDosageSum),
  kfScoreColAll = ((kfScoreColScoreSums * 2) - kfScoreColMaybefid)
FLAGSET_DEF_END(ScoreFlags);
typedef struct ScoreInfoStruct {
  NONCOPYABLE(ScoreInfoStruct);
  ScoreFlags flags;
  uint32_t varid_col_p1;
  uint32_t allele_col_p1;
  char* input_fname;
  RangeList input_col_idx_range_list;

  // --q-score-range
  char* qsr_range_fname;
  char* qsr_data_fname;
  uint32_t qsr_varid_col_p1;
  uint32_t qsr_val_col_p1;

  // --score-partitions.  Fills the same qsr_include membership bitarray as
  // --q-score-range (hence mutually exclusive with it) from an explicit 0/1
  // mask, and writes one merged file instead of one file per subset.
  char* partition_fname;
} ScoreInfo;
void InitScore(ScoreInfo* score_info_ptr);

void CleanupScore(ScoreInfo* score_info_ptr);


CONSTI32(kMaxPc, 8000);
#ifndef NOLAPACK
#endif

PglErr ScoreReport(const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uintptr_t* variant_include, const ChrInfo* cip, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const double* allele_freqs, const ScoreInfo* score_info_ptr, const char* output_missing_pheno, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t nosex_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_variant_id_slen, uint32_t xchr_model, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end);


#ifndef NOLAPACK
#endif

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __ADJUTANT_MATRIX_CALC_H__
