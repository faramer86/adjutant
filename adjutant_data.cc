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

#include "adjutant_data.h"

#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>  // time()
#include <unistd.h>  // unlink()

#include "include/plink2_bgzf.h"
#include "include/plink2_bits.h"
#include "include/plink2_htable.h"
#include "include/plink2_simd.h"
#include "include/plink2_string.h"
#include "include/plink2_thread.h"
#include "adjutant_compress_stream.h"
#include "adjutant_decompress.h"
#include "adjutant_pvar.h"

#ifdef __cplusplus
namespace plink2 {
#endif
// May use all remaining workspace memory.
// fileformat, fileDate, source
/*
// Note that the order-of-operations page lists this as happening right after
// the filtering performed by LoadPvar().  Which is effectively true, since we
// ignore variant_include (this is safe since LoadPvar() always initializes
// all variant_bps[] and allele_storage[] entries appropriately).
// possible todo: ChrInfo can have a length field, which is initialized by the
// ##contig header line when possible, but when that doesn't exist LoadPvar()
// can conditionally detect INFO/END and take that into account.  (Or a reason
// to keep the entire info_end array in memory may emerge.)
//
// Update (10 Apr 2022): Discontinuing lower-bound length= entries in exported
// VCF/BCF headers, since --fa can now be used whenever length values are
// needed.
*/
// Read: assumes hkvline_iter, idval, and id_slen are return values from
//       HkvlineId(), and closing_gt points to the closing '>'.
// Write: assumes ID key=value has been written, but NOT trailing comma, and
//        finishes in an analogous state.
/*
// unaligned accesses here.
*/

// Now assumes trailing bits of genovec have been zeroed out.
// writebuf need not be word-aligned.
// Revised phaseraw:
//   4 byte het_ct, 4 byte explicit_phasepresent_ct
//   first half, up to (1 + (het_ct / kBitsPerWord)) words
//   second half, rounded up to vector boundary
// these also work on dphaseraw
// Requires trailing bits of genovec to be zeroed out.
// "Flat" = don't separate one_cts and two_cts.
void GetMFlatCounts64(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const PgenVariant* pgvp, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t allele_ct, STD_ARRAY_REF(uint32_t, 4) genocounts, uint64_t* all_dosages) {
  if (sample_ct == raw_sample_ct) {
    GenoarrCountFreqsUnsafe(pgvp->genovec, sample_ct, genocounts);
  } else {
    GenoarrCountSubsetFreqs(pgvp->genovec, sample_include_interleaved_vec, raw_sample_ct, sample_ct, genocounts);
  }
  all_dosages[0] = 2 * genocounts[0] + genocounts[1];
  all_dosages[1] = 2 * genocounts[2] + genocounts[1];
  ZeroU64Arr(allele_ct - 2, &(all_dosages[2]));
  const AlleleCode* patch_01_vals = pgvp->patch_01_vals;
  const AlleleCode* patch_10_vals = pgvp->patch_10_vals;
  const uint32_t patch_01_ct = pgvp->patch_01_ct;
  const uint32_t patch_10_ct = pgvp->patch_10_ct;
  if (sample_ct == raw_sample_ct) {
    all_dosages[1] -= patch_01_ct + 2 * patch_10_ct;
    for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
      all_dosages[patch_01_vals[uii]] += 1;
    }
    const uint32_t patch_10_ct_x2 = patch_10_ct * 2;
    for (uint32_t uii = 0; uii != patch_10_ct_x2; ++uii) {
      all_dosages[patch_10_vals[uii]] += 1;
    }
  } else {
    if (patch_01_ct) {
      const uintptr_t* patch_01_set = pgvp->patch_01_set;
      uintptr_t sample_widx = 0;
      uintptr_t patch_01_bits = patch_01_set[0];
      uint32_t subsetted_patch_01_ct = 0;
      for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
        const uintptr_t lowbit = BitIter1y(patch_01_set, &sample_widx, &patch_01_bits);
        if (sample_include[sample_widx] & lowbit) {
          all_dosages[patch_01_vals[uii]] += 1;
          ++subsetted_patch_01_ct;
        }
      }
      all_dosages[1] -= subsetted_patch_01_ct;
    }
    if (patch_10_ct) {
      const uintptr_t* patch_10_set = pgvp->patch_10_set;
      uintptr_t sample_widx = 0;
      uintptr_t patch_10_bits = patch_10_set[0];
      uint32_t subsetted_patch_10_ct = 0;
      for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
        const uintptr_t lowbit = BitIter1y(patch_10_set, &sample_widx, &patch_10_bits);
        if (sample_include[sample_widx] & lowbit) {
          all_dosages[patch_10_vals[2 * uii]] += 1;
          all_dosages[patch_10_vals[2 * uii + 1]] += 1;
          ++subsetted_patch_10_ct;
        }
      }
      all_dosages[1] -= 2 * subsetted_patch_10_ct;
    }
  }
}

void GetMCounts64(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const PgenVariant* pgvp, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t allele_ct, STD_ARRAY_REF(uint32_t, 4) genocounts, uint64_t* __restrict one_cts, uint64_t* __restrict two_cts) {
  // This mirrors GetMultiallelicCountsAndDosage16s().
  if (sample_ct == raw_sample_ct) {
    GenoarrCountFreqsUnsafe(pgvp->genovec, sample_ct, genocounts);
  } else {
    GenoarrCountSubsetFreqs(pgvp->genovec, sample_include_interleaved_vec, raw_sample_ct, sample_ct, genocounts);
  }
  one_cts[0] = genocounts[1];
  one_cts[1] = genocounts[1];
  ZeroU64Arr(allele_ct - 2, &(one_cts[2]));
  two_cts[0] = genocounts[0];
  two_cts[1] = genocounts[2];
  ZeroU64Arr(allele_ct - 2, &(two_cts[2]));
  const AlleleCode* patch_01_vals = pgvp->patch_01_vals;
  const AlleleCode* patch_10_vals = pgvp->patch_10_vals;
  const uint32_t patch_01_ct = pgvp->patch_01_ct;
  const uint32_t patch_10_ct = pgvp->patch_10_ct;
  if (sample_ct == raw_sample_ct) {
    one_cts[1] -= patch_01_ct;
    for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
      one_cts[patch_01_vals[uii]] += 1;
    }
    two_cts[1] -= patch_10_ct;
    const AlleleCode* patch_10_vals_iter = patch_10_vals;
    for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
      const AlleleCode code_lo = *patch_10_vals_iter++;
      const AlleleCode code_hi = *patch_10_vals_iter++;
      if (code_lo == code_hi) {
        two_cts[code_lo] += 1;
      } else {
        one_cts[code_lo] += 1;
        one_cts[code_hi] += 1;
      }
    }
  } else {
    if (patch_01_ct) {
      const uintptr_t* patch_01_set = pgvp->patch_01_set;
      uintptr_t sample_widx = 0;
      uintptr_t patch_01_bits = patch_01_set[0];
      uint32_t subsetted_patch_01_ct = 0;
      for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
        const uintptr_t lowbit = BitIter1y(patch_01_set, &sample_widx, &patch_01_bits);
        if (sample_include[sample_widx] & lowbit) {
          one_cts[patch_01_vals[uii]] += 1;
          ++subsetted_patch_01_ct;
        }
      }
      one_cts[1] -= subsetted_patch_01_ct;
    }
    if (patch_10_ct) {
      const uintptr_t* patch_10_set = pgvp->patch_10_set;
      uintptr_t sample_widx = 0;
      uintptr_t patch_10_bits = patch_10_set[0];
      uint32_t subsetted_patch_10_ct = 0;
      for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
        const uintptr_t lowbit = BitIter1y(patch_10_set, &sample_widx, &patch_10_bits);
        if (sample_include[sample_widx] & lowbit) {
          ++subsetted_patch_10_ct;
          const AlleleCode code_lo = patch_10_vals[2 * uii];
          const AlleleCode code_hi = patch_10_vals[2 * uii + 1];
          if (code_lo == code_hi) {
            two_cts[code_lo] += 1;
          } else {
            one_cts[code_lo] += 1;
            one_cts[code_hi] += 1;
          }
        }
      }
      two_cts[1] -= subsetted_patch_10_ct;
    }
  }
}

typedef struct LoadAlleleAndGenoCountsCtxStruct {
  const uintptr_t* variant_include;
  const ChrInfo* cip;
  const uintptr_t* allele_idx_offsets;
  const uintptr_t* sample_include;
  uintptr_t* sample_include_interleaved_vec;
  uint32_t* sample_include_cumulative_popcounts;
  const uintptr_t* sex_male;
  const uintptr_t* sex_nonfemale;
  uintptr_t* sex_male_interleaved_vec;
  uintptr_t* sex_nonfemale_interleaved_vec;
  uint32_t* sex_nonfemale_cumulative_popcounts;
  uintptr_t* nosex_interleaved_vec;
  const uintptr_t* founder_info;
  uintptr_t* founder_info_interleaved_vec;
  uint32_t* founder_info_cumulative_popcounts;
  uintptr_t* founder_male;
  uintptr_t* founder_nonfemale;
  uintptr_t* founder_male_interleaved_vec;
  uintptr_t* founder_nonfemale_interleaved_vec;
  uint32_t* founder_nonfemale_cumulative_popcounts;
  uintptr_t* founder_nosex_interleaved_vec;
  uint32_t raw_sample_ct;
  uint32_t sample_ct;
  uint32_t founder_ct;
  uint32_t male_ct;
  uint32_t nosex_ct;
  uint32_t chry_missingstat_sample_ct;
  uint32_t founder_male_ct;
  uint32_t founder_nosex_ct;
  uint32_t founder_chry_missingstat_sample_ct;
  uint32_t first_hap_uidx;
  uint32_t is_minimac3_r2;

  PgenReader** pgr_ptrs;

  uintptr_t** genovecs;
  uintptr_t** thread_read_mhc;
  uintptr_t** dosage_presents;
  Dosage** dosage_mains;
  uint64_t** all_dosages;
  uint32_t* read_variant_uidx_starts;

  // shouldn't need array, or errno storage, since kPglRetMalformedInput is the
  // only possible error for now
  uint64_t err_info;

  uint32_t cur_block_size;

  unsigned char* allele_presents_bytearr;
  uint64_t* allele_ddosages;
  STD_ARRAY_PTR_DECL(uint32_t, 3, raw_geno_cts);
  uint32_t* variant_missing_hc_cts;
  uint32_t* variant_missing_dosage_cts;
  uint32_t* variant_hethap_cts;
  uint64_t* founder_allele_ddosages;
  STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts);
  STD_ARRAY_PTR_DECL(uint32_t, 3, x_male_geno_cts);
  STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts);
  STD_ARRAY_PTR_DECL(uint32_t, 3, x_nosex_geno_cts);
  STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts);
  double* imp_r2_vals;
} LoadAlleleAndGenoCountsCtx;

THREAD_FUNC_DECL LoadAlleleAndGenoCountsThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  LoadAlleleAndGenoCountsCtx* ctx = S_CAST(LoadAlleleAndGenoCountsCtx*, arg->sharedp->context);

  const uintptr_t* variant_include = ctx->variant_include;
  const ChrInfo* cip = ctx->cip;
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const uint32_t thread_ct = GetThreadCt(arg->sharedp);
  const uint32_t subset_ct = (ctx->founder_info != nullptr) + 1;
  const uint32_t raw_sample_ct = ctx->raw_sample_ct;
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t first_hap_uidx = ctx->first_hap_uidx;
  const uint32_t is_minimac3_r2 = ctx->is_minimac3_r2;
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  PgenVariant pgv;
  pgv.genovec = ctx->genovecs[tidx];
  SetPgvThreadMhcNull(raw_sample_ct, tidx, ctx->thread_read_mhc, &pgv);
  pgv.dosage_present = nullptr;
  pgv.dosage_main = nullptr;
  if (ctx->dosage_presents) {
    pgv.dosage_present = ctx->dosage_presents[tidx];
    pgv.dosage_main = ctx->dosage_mains[tidx];
  }
  uint64_t* all_dosages = nullptr;
  if (ctx->all_dosages) {
    all_dosages = ctx->all_dosages[tidx];
  }
  uint32_t is_y = 0;
  uint32_t is_nonxy_haploid = 0;
  uint32_t x_start = 0;
  uint32_t x_code;
  if (XymtExists(cip, kChrOffsetX, &x_code)) {
    const uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[x_code];
    x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
  }
  uint32_t allele_ct = 2;
  uint64_t new_err_info;
  do {
    const uintptr_t cur_block_size = ctx->cur_block_size;
    // no overflow danger since cur_block_size <= 2^16, tidx < (2^16 - 1)
    const uint32_t cur_idx_end = ((tidx + 1) * cur_block_size) / thread_ct;
    const uintptr_t* sample_include = ctx->sample_include;
    const uintptr_t* sample_include_interleaved_vec = ctx->sample_include_interleaved_vec;
    const uint32_t* sample_include_cumulative_popcounts = ctx->sample_include_cumulative_popcounts;
    const uintptr_t* sex_male = ctx->sex_male;
    const uintptr_t* sex_male_interleaved_vec = ctx->sex_male_interleaved_vec;
    const uintptr_t* sex_nonfemale = ctx->sex_nonfemale;
    const uintptr_t* sex_nonfemale_interleaved_vec = ctx->sex_nonfemale_interleaved_vec;
    const uint32_t* sex_nonfemale_cumulative_popcounts = ctx->sex_nonfemale_cumulative_popcounts;
    const uintptr_t* nosex_interleaved_vec = ctx->nosex_interleaved_vec;
    uint32_t sample_ct = ctx->sample_ct;
    uint32_t male_ct = ctx->male_ct;
    uint32_t nosex_ct = ctx->nosex_ct;
    uint32_t nonfemale_ct = male_ct + nosex_ct;
    uint32_t chry_missingstat_sample_ct = ctx->chry_missingstat_sample_ct;
    unsigned char* allele_presents_bytearr = ctx->allele_presents_bytearr;
    uint64_t* allele_ddosages = ctx->allele_ddosages;
    STD_ARRAY_PTR_DECL(uint32_t, 3, raw_geno_cts) = ctx->raw_geno_cts;
    uint32_t* variant_missing_hc_cts = ctx->variant_missing_hc_cts;
    uint32_t* variant_missing_dosage_cts = ctx->variant_missing_dosage_cts;
    uint32_t* variant_hethap_cts = ctx->variant_hethap_cts;
    STD_ARRAY_PTR_DECL(uint32_t, 3, x_male_geno_cts) = ctx->x_male_geno_cts;
    STD_ARRAY_PTR_DECL(uint32_t, 3, x_nosex_geno_cts) = ctx->x_nosex_geno_cts;
    double* imp_r2_vals = ctx->imp_r2_vals;
    pgv.dosage_ct = 0;
    for (uint32_t subset_idx = 0; ; ) {
      // bugfix (29 Dec 2019): this boolean can change with subset_idx
      const uint32_t no_multiallelic_branch = (!variant_hethap_cts) && (!allele_presents_bytearr) && (!allele_ddosages) && (!imp_r2_vals);
      PgrSampleSubsetIndex pssi;
      PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, pgrp, &pssi);
      uint32_t cur_idx = (tidx * cur_block_size) / thread_ct;
      uintptr_t variant_uidx_base;
      uintptr_t variant_include_bits;
      BitIter1Start(variant_include, ctx->read_variant_uidx_starts[tidx], &variant_uidx_base, &variant_include_bits);
      uint32_t chr_end = 0;
      uint32_t is_x_or_y = 0;
      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      STD_ARRAY_DECL(uint32_t, 4, sex_specific_genocounts);
      for (; cur_idx != cur_idx_end; ++cur_idx) {
        const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        if (variant_uidx >= chr_end) {
          const uint32_t chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx);
          const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
          is_y = 0;
          is_nonxy_haploid = 0;
          if (chr_idx == x_code) {
            is_x_or_y = 1;
            PgrClearSampleSubsetIndex(pgrp, &pssi);
          } else if (chr_idx == y_code) {
            is_x_or_y = 1;
            is_y = 1;
            // ugh
            if ((nonfemale_ct == chry_missingstat_sample_ct) && ((!allele_presents_bytearr) || (sample_ct == nonfemale_ct))) {
              PgrSetSampleSubsetIndex(sex_nonfemale_cumulative_popcounts, pgrp, &pssi);
            } else {
              PgrClearSampleSubsetIndex(pgrp, &pssi);
            }
          } else {
            if (is_x_or_y) {
              PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, pgrp, &pssi);
            }
            is_x_or_y = 0;
            // true for MT
            is_nonxy_haploid = IsSet(cip->haploid_mask, chr_idx);
          }
        }
        uintptr_t cur_allele_idx_offset;
        if (!allele_idx_offsets) {
          cur_allele_idx_offset = 2 * variant_uidx;
        } else {
          cur_allele_idx_offset = allele_idx_offsets[variant_uidx];
          allele_ct = allele_idx_offsets[variant_uidx + 1] - cur_allele_idx_offset;
        }
        uint32_t hethap_ct;
        if ((allele_ct == 2) || no_multiallelic_branch) {
          uint64_t cur_dosages[2];
          if (!is_x_or_y) {
            const PglErr reterr = PgrGetDCounts(sample_include, sample_include_interleaved_vec, pssi, sample_ct, variant_uidx, is_minimac3_r2, pgrp, imp_r2_vals? (&(imp_r2_vals[variant_uidx])) : nullptr, genocounts, cur_dosages);
            if (unlikely(reterr)) {
              new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
              goto LoadAlleleAndGenoCountsThread_err;
            }
            if (allele_presents_bytearr) {
              if (cur_dosages[0]) {
                allele_presents_bytearr[cur_allele_idx_offset] = 128;
              }
              if (cur_dosages[1]) {
                allele_presents_bytearr[cur_allele_idx_offset + 1] = 128;
              }
            }
            if (!is_nonxy_haploid) {
              hethap_ct = 0;
              if (allele_ddosages) {
                // ...but save all allele counts here.
                allele_ddosages[cur_allele_idx_offset] = cur_dosages[0] * 2;
                allele_ddosages[cur_allele_idx_offset + 1] = cur_dosages[1] * 2;
              }
            } else {
              // this hethap_ct can be inaccurate in multiallelic case
              hethap_ct = genocounts[1];
              if (imp_r2_vals && (!is_minimac3_r2)) {
                // Assuming the input data isn't malformed "phased haploid",
                // minimac3-r2 is independent of haploid/diploid state; only
                // mach-r2 requires a haploid correction.
                imp_r2_vals[variant_uidx] *= 0.5;
              }
              if (allele_ddosages) {
                allele_ddosages[cur_allele_idx_offset] = cur_dosages[0];
                allele_ddosages[cur_allele_idx_offset + 1] = cur_dosages[1];
              }
            }
          } else if (is_y) {
            if ((nonfemale_ct == chry_missingstat_sample_ct) && ((!allele_presents_bytearr) || (sample_ct == nonfemale_ct))) {
              const PglErr reterr = PgrGetDCounts(sex_nonfemale, sex_nonfemale_interleaved_vec, pssi, nonfemale_ct, variant_uidx, 0, pgrp, imp_r2_vals? (&(imp_r2_vals[variant_uidx])) : nullptr, genocounts, cur_dosages);
              if (unlikely(reterr)) {
                new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
                goto LoadAlleleAndGenoCountsThread_err;
              }
              hethap_ct = genocounts[1];
              if (imp_r2_vals && (!is_minimac3_r2)) {
                // note that female is not counted here
                imp_r2_vals[variant_uidx] *= 0.5;
              }
              if (allele_presents_bytearr) {
                if (cur_dosages[0]) {
                  allele_presents_bytearr[cur_allele_idx_offset] = 128;
                }
                if (cur_dosages[1]) {
                  allele_presents_bytearr[cur_allele_idx_offset + 1] = 128;
                }
              }
              if (allele_ddosages) {
                allele_ddosages[cur_allele_idx_offset] = cur_dosages[0];
                allele_ddosages[cur_allele_idx_offset + 1] = cur_dosages[1];
              }
              if (raw_geno_cts) {
                STD_ARRAY_REF(uint32_t, 3) cur_raw_geno_cts = raw_geno_cts[variant_uidx];
                cur_raw_geno_cts[0] = genocounts[0];
                cur_raw_geno_cts[1] = genocounts[1];
                cur_raw_geno_cts[2] = genocounts[2];
              }
            } else {
              // females need to be counted for allele_presents, and ignored
              // elsewhere.
              // unknown-sex might need to be excluded from missing-genotype
              // and hethap counts.
              const PglErr reterr = PgrGetD(nullptr, pssi, raw_sample_ct, variant_uidx, pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &pgv.dosage_ct);
              if (unlikely(reterr)) {
                new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
                goto LoadAlleleAndGenoCountsThread_err;
              }
              const uint32_t dosage_is_relevant = pgv.dosage_ct && ((sample_ct == raw_sample_ct) || (!IntersectionIsEmpty(sample_include, pgv.dosage_present, raw_sample_ctl)));
              if (allele_presents_bytearr) {
                if (dosage_is_relevant) {
                  // at least one dosage value is present, that's all we need
                  // to know
                  allele_presents_bytearr[cur_allele_idx_offset] = 128;
                  allele_presents_bytearr[cur_allele_idx_offset + 1] = 128;
                } else {
                  // only hardcalls matter
                  // bugfix (31 Jul 2018): forgot to initialize genocounts here
                  // possible todo: use a specialized function which just
                  // checks which alleles exist
                  if (sample_ct == raw_sample_ct) {
                    ZeroTrailingNyps(raw_sample_ct, pgv.genovec);
                    GenoarrCountFreqsUnsafe(pgv.genovec, sample_ct, genocounts);
                  } else {
                    GenoarrCountSubsetFreqs(pgv.genovec, sample_include_interleaved_vec, raw_sample_ct, sample_ct, genocounts);
                  }
                  if (genocounts[0] || genocounts[1]) {
                    allele_presents_bytearr[cur_allele_idx_offset] = 128;
                  }
                  if (genocounts[1] || genocounts[2]) {
                    allele_presents_bytearr[cur_allele_idx_offset + 1] = 128;
                  }
                }
              }
              GenoarrCountSubsetFreqs(pgv.genovec, sex_nonfemale_interleaved_vec, raw_sample_ct, nonfemale_ct, genocounts);
              hethap_ct = genocounts[1];
              // x2, x4 since this is haploid
              uintptr_t alt1_ct_x2 = genocounts[2] * 2 + hethap_ct;
              uintptr_t alt1_sq_sum_x4 = genocounts[2] * (4 * k1LU) + hethap_ct;
              uint64_t alt1_ddosage = 0;  // in 32768ths
              uint64_t alt1_ddosage_sq_sum = 0;
              uint32_t additional_dosage_ct = 0;
              if (dosage_is_relevant) {
                uintptr_t sample_widx = 0;
                uintptr_t dosage_present_bits = pgv.dosage_present[0];
                uint32_t sample_uidx = 0;
                for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ++dosage_idx) {
                  const uintptr_t lowbit = BitIter1y(pgv.dosage_present, &sample_widx, &dosage_present_bits);
                  if (sample_include[sample_widx] & lowbit) {
                    const uintptr_t cur_dosage_val = pgv.dosage_main[dosage_idx];
                    alt1_ddosage += cur_dosage_val;
                    alt1_ddosage_sq_sum += cur_dosage_val * cur_dosage_val;
                    const uintptr_t hardcall_code = GetNyparrEntry(pgv.genovec, sample_uidx);
                    if (hardcall_code != 3) {
                      alt1_ct_x2 -= hardcall_code;
                      alt1_sq_sum_x4 -= hardcall_code * hardcall_code;
                    } else {
                      ++additional_dosage_ct;
                    }
                  }
                }
              }
              const uintptr_t obs_ct = nonfemale_ct + additional_dosage_ct - genocounts[3];
              alt1_ddosage += alt1_ct_x2 * S_CAST(uint64_t, kDosageMid);
              alt1_ddosage_sq_sum += alt1_sq_sum_x4 * 0x10000000LLU;
              cur_dosages[0] = obs_ct * S_CAST(uint64_t, kDosageMax) - alt1_ddosage;
              cur_dosages[1] = alt1_ddosage;
              if (imp_r2_vals) {
                // minimac3-r2 and mach-r2 are identical in haploid case
                const double dosage_sumd = u63tod(alt1_ddosage);
                const double dosage_avg = dosage_sumd / u31tod(obs_ct);
                const double dosage_variance = u63tod(alt1_ddosage_sq_sum) - dosage_sumd * dosage_avg;
                imp_r2_vals[variant_uidx] = dosage_variance / (dosage_sumd * (32768 - dosage_avg));
              }
              if (allele_ddosages) {
                allele_ddosages[cur_allele_idx_offset] = cur_dosages[0];
                allele_ddosages[cur_allele_idx_offset + 1] = alt1_ddosage;
              }
              // argh, these values must be saved before potential male-only
              // clobber.
              if (raw_geno_cts) {
                STD_ARRAY_REF(uint32_t, 3) cur_raw_geno_cts = raw_geno_cts[variant_uidx];
                cur_raw_geno_cts[0] = genocounts[0];
                cur_raw_geno_cts[1] = genocounts[1];
                cur_raw_geno_cts[2] = genocounts[2];
              }
              if (chry_missingstat_sample_ct != nonfemale_ct) {
                // By default, we include unknown-sex samples in
                // allele-frequency and imputation-r2 computations, but exclude
                // them from missing-rate stats (and hethap stats because they
                // can't be cleanly separated out)
                GenoarrCountSubsetFreqs(pgv.genovec, sex_male_interleaved_vec, raw_sample_ct, male_ct, genocounts);
                hethap_ct = genocounts[1];
                // cur_dosages[0] + cur_dosages[1] used later to compute
                // missing-dosage rate.
                uint32_t male_missing_ct = genocounts[3];
                if (dosage_is_relevant && male_missing_ct) {
                  const uint32_t raw_sample_ctl2 = NypCtToWordCt(raw_sample_ct);
                  const uintptr_t* genovec = pgv.genovec;
                  const Halfword* sex_male_hwalias = DowncastKWToHW(sex_male);
                  const Halfword* dosage_present_hwalias = DowncastWToHW(pgv.dosage_present);
                  for (uint32_t widx = 0; widx != raw_sample_ctl2; ++widx) {
                    const Halfword hw = Pack11ToHalfword(genovec[widx]) & sex_male_hwalias[widx] & dosage_present_hwalias[widx];
                    if (hw) {
                      male_missing_ct -= PopcountHW(hw);
                    }
                  }
                }
                const uintptr_t male_dosage_obs_ct = male_ct - male_missing_ct;
                cur_dosages[0] = male_dosage_obs_ct * S_CAST(uint64_t, kDosageMax);
                cur_dosages[1] = 0;
              }
            }
          } else {
            // chrX
            const PglErr reterr = PgrGetD(nullptr, pssi, raw_sample_ct, variant_uidx, pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &pgv.dosage_ct);
            if (unlikely(reterr)) {
              new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
              goto LoadAlleleAndGenoCountsThread_err;
            }
            if (sample_ct == raw_sample_ct) {
              ZeroTrailingNyps(raw_sample_ct, pgv.genovec);
              GenoarrCountFreqsUnsafe(pgv.genovec, sample_ct, genocounts);
            } else {
              GenoarrCountSubsetFreqs(pgv.genovec, sample_include_interleaved_vec, raw_sample_ct, sample_ct, genocounts);
            }
            GenoarrCountSubsetFreqs(pgv.genovec, sex_male_interleaved_vec, raw_sample_ct, male_ct, sex_specific_genocounts);
            hethap_ct = sex_specific_genocounts[1];
            // Could compute imputation r2 iff there are no unknown-sex
            // samples, but probably not worth it since larger datasets could
            // have a small number of Klinefelter syndrome cases, etc. coded as
            // unknown-sex, and we don't want to discourage their inclusion;
            // let's delegate that chrX filter to other software for now.

            if (allele_presents_bytearr) {
              if (pgv.dosage_ct && ((sample_ct == raw_sample_ct) || (!IntersectionIsEmpty(sample_include, pgv.dosage_present, raw_sample_ctl)))) {
                // at least one dosage value is present, that's all we need to
                // know
                allele_presents_bytearr[cur_allele_idx_offset] = 128;
                allele_presents_bytearr[cur_allele_idx_offset + 1] = 128;
              } else {
                // only hardcalls matter
                if (genocounts[0] || genocounts[1]) {
                  allele_presents_bytearr[cur_allele_idx_offset] = 128;
                }
                if (genocounts[1] || genocounts[2]) {
                  allele_presents_bytearr[cur_allele_idx_offset + 1] = 128;
                }
              }
            }
            if (allele_ddosages) {
              uintptr_t alt1_ct = 4 * genocounts[2] + 2 * genocounts[1] - 2 * sex_specific_genocounts[2] - hethap_ct;  // nonmales count twice
              uint64_t alt1_ddosage = 0;  // in 32768ths, nonmales count twice
              uint32_t additional_dosage_ct = 0;  // missing hardcalls only; nonmales count twice
              // bugfix (12 Jul 2018): dosage_present may be null if dosage_ct
              // == 0
              if (pgv.dosage_ct) {
                uintptr_t sample_uidx_base = 0;
                uintptr_t dosage_present_bits = pgv.dosage_present[0];
                if (sample_ct == raw_sample_ct) {
                  for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ++dosage_idx) {
                    const uintptr_t sample_uidx = BitIter1(pgv.dosage_present, &sample_uidx_base, &dosage_present_bits);
                    const uintptr_t cur_dosage_val = pgv.dosage_main[dosage_idx];
                    const uintptr_t sex_multiplier = 2 - IsSet(sex_male, sample_uidx);
                    alt1_ddosage += cur_dosage_val * sex_multiplier;

                    // could call GenoarrCountSubsetIntersectFreqs() twice
                    // instead, but since we've already manually extracted the
                    // sex bit it probably doesn't help?
                    const uintptr_t hardcall_code = GetNyparrEntry(pgv.genovec, sample_uidx);
                    if (hardcall_code != 3) {
                      alt1_ct -= hardcall_code * sex_multiplier;
                    } else {
                      additional_dosage_ct += sex_multiplier;
                    }
                  }
                } else {
                  for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ++dosage_idx) {
                    const uintptr_t sample_uidx = BitIter1(pgv.dosage_present, &sample_uidx_base, &dosage_present_bits);
                    if (IsSet(sample_include, sample_uidx)) {
                      const uintptr_t cur_dosage_val = pgv.dosage_main[dosage_idx];
                      const uintptr_t sex_multiplier = 2 - IsSet(sex_male, sample_uidx);
                      alt1_ddosage += cur_dosage_val * sex_multiplier;
                      const uintptr_t hardcall_code = GetNyparrEntry(pgv.genovec, sample_uidx);
                      if (hardcall_code != 3) {
                        alt1_ct -= hardcall_code * sex_multiplier;
                      } else {
                        additional_dosage_ct += sex_multiplier;
                      }
                    }
                  }
                }
              }
              alt1_ddosage += alt1_ct * S_CAST(uint64_t, kDosageMid);

              // bugfix (14 May 2018): this didn't correctly distinguish
              // between missing vs. 'replaced' hardcalls
              const uintptr_t weighted_obs_ct = (2 * (sample_ct - genocounts[3]) - male_ct + sex_specific_genocounts[3] + additional_dosage_ct) * (2 * k1LU);

              allele_ddosages[cur_allele_idx_offset] = weighted_obs_ct * S_CAST(uint64_t, kDosageMid) - alt1_ddosage;
              allele_ddosages[cur_allele_idx_offset + 1] = alt1_ddosage;
            }
            if (x_male_geno_cts) {
              STD_ARRAY_REF(uint32_t, 3) cur_x_male_geno_cts = x_male_geno_cts[variant_uidx - x_start];
              cur_x_male_geno_cts[0] = sex_specific_genocounts[0];
              cur_x_male_geno_cts[1] = sex_specific_genocounts[1];
              cur_x_male_geno_cts[2] = sex_specific_genocounts[2];
              if (x_nosex_geno_cts) {
                GenoarrCountSubsetFreqs(pgv.genovec, nosex_interleaved_vec, raw_sample_ct, nosex_ct, sex_specific_genocounts);
                STD_ARRAY_REF(uint32_t, 3) cur_nosex_geno_cts = x_nosex_geno_cts[variant_uidx - x_start];
                cur_nosex_geno_cts[0] = sex_specific_genocounts[0];
                cur_nosex_geno_cts[1] = sex_specific_genocounts[1];
                cur_nosex_geno_cts[2] = sex_specific_genocounts[2];
              }
            }
          }
          if (variant_missing_dosage_cts) {
            uint32_t missing_dosage_ct;
            if (!is_x_or_y) {
              missing_dosage_ct = sample_ct - ((cur_dosages[0] + cur_dosages[1]) / kDosageMax);
            } else if (is_y) {
              missing_dosage_ct = chry_missingstat_sample_ct - ((cur_dosages[0] + cur_dosages[1]) / kDosageMax);
            } else {
              if (pgv.dosage_ct) {
                ZeroTrailingNyps(raw_sample_ct, pgv.genovec);
                missing_dosage_ct = GenoarrCountMissingInvsubsetUnsafe(pgv.genovec, pgv.dosage_present, raw_sample_ct);
              } else {
                missing_dosage_ct = genocounts[3];
              }
            }
            variant_missing_dosage_cts[variant_uidx] = missing_dosage_ct;
          }
        } else {
          // multiallelic cases
          if (!is_x_or_y) {
            const PglErr reterr = PgrGetMDCounts(sample_include, sample_include_interleaved_vec, pssi, sample_ct, variant_uidx, is_minimac3_r2, pgrp, imp_r2_vals? (&(imp_r2_vals[variant_uidx])) : nullptr, &hethap_ct, genocounts, all_dosages);
            if (unlikely(reterr)) {
              new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
              goto LoadAlleleAndGenoCountsThread_err;
            }
            if (allele_presents_bytearr) {
              for (uintptr_t aidx = 0; aidx != allele_ct; ++aidx) {
                if (all_dosages[aidx]) {
                  allele_presents_bytearr[cur_allele_idx_offset + aidx] = 128;
                }
              }
            }
            if (!is_nonxy_haploid) {
              hethap_ct = 0;
              if (allele_ddosages) {
                for (uintptr_t aidx = 0; aidx != allele_ct; ++aidx) {
                  allele_ddosages[cur_allele_idx_offset + aidx] = all_dosages[aidx] * 2;
                }
              }
            } else {
              if (imp_r2_vals && (!is_minimac3_r2)) {
                imp_r2_vals[variant_uidx] *= 0.5;
              }
              if (allele_ddosages) {
                memcpy(&(allele_ddosages[cur_allele_idx_offset]), all_dosages, allele_ct * sizeof(int64_t));
              }
            }
          } else if (is_y) {
            if ((nonfemale_ct == chry_missingstat_sample_ct) && ((!allele_presents_bytearr) || (sample_ct == nonfemale_ct))) {
              const PglErr reterr = PgrGetMDCounts(sex_nonfemale, sex_nonfemale_interleaved_vec, pssi, nonfemale_ct, variant_uidx, 0, pgrp, imp_r2_vals? (&(imp_r2_vals[variant_uidx])) : nullptr, &hethap_ct, genocounts, all_dosages);
              if (unlikely(reterr)) {
                new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
                goto LoadAlleleAndGenoCountsThread_err;
              }
              if (imp_r2_vals && (!is_minimac3_r2)) {
                imp_r2_vals[variant_uidx] *= 0.5;
              }
              if (allele_presents_bytearr) {
                for (uintptr_t aidx = 0; aidx != allele_ct; ++aidx) {
                  if (all_dosages[aidx]) {
                    allele_presents_bytearr[cur_allele_idx_offset + aidx] = 128;
                  }
                }
              }
              if (allele_ddosages) {
                memcpy(&(allele_ddosages[cur_allele_idx_offset]), all_dosages, allele_ct * sizeof(int64_t));
              }
              if (raw_geno_cts) {
                STD_ARRAY_REF(uint32_t, 3) cur_raw_geno_cts = raw_geno_cts[variant_uidx];
                cur_raw_geno_cts[0] = genocounts[0];
                cur_raw_geno_cts[1] = genocounts[1];
                cur_raw_geno_cts[2] = genocounts[2];
              }
            } else {
              // females need to be counted for allele_presents, and ignored
              // elsewhere.
              // unknown-sex might need to be excluded from missing-genotype
              // and hethap counts.
              const PglErr reterr = PgrGetM(nullptr, pssi, raw_sample_ct, variant_uidx, pgrp, &pgv);
              if (unlikely(reterr)) {
                new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
                goto LoadAlleleAndGenoCountsThread_err;
              }
              // possible todo: use a specialized function which just checks
              // which alleles exist
              ZeroTrailingNyps(raw_sample_ct, pgv.genovec);
              GetMFlatCounts64(sample_include, sample_include_interleaved_vec, &pgv, raw_sample_ct, sample_ct, allele_ct, genocounts, all_dosages);
              if (allele_presents_bytearr) {
                for (uintptr_t aidx = 0; aidx != allele_ct; ++aidx) {
                  if (all_dosages[aidx]) {
                    allele_presents_bytearr[cur_allele_idx_offset + aidx] = 128;
                  }
                }
              }

              uint64_t* two_cts = &(all_dosages[allele_ct]);
              GetMCounts64(sex_nonfemale, sex_nonfemale_interleaved_vec, &pgv, raw_sample_ct, nonfemale_ct, allele_ct, genocounts, all_dosages, two_cts);
              if (allele_ddosages) {
                for (uintptr_t aidx = 0; aidx != allele_ct; ++aidx) {
                  allele_ddosages[cur_allele_idx_offset + aidx] = all_dosages[aidx] * kDosageMid + two_cts[aidx] * kDosageMax;
                }
              }
              if (imp_r2_vals) {
                for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
                  const uint64_t one_ct = allele_ddosages[aidx];
                  const uint64_t two_ct = two_cts[aidx];
                  // now sums
                  allele_ddosages[aidx] = one_ct * kDosageMid + two_ct * kDosageMax;
                  // now ssqs
                  two_cts[aidx] = one_ct * kDosageMid * kDosageMid + two_ct * kDosageMax * kDosageMax;
                }
                imp_r2_vals[variant_uidx] = 0.5 * MultiallelicDiploidMachR2(all_dosages, two_cts, nonfemale_ct - genocounts[3], allele_ct);
              }
              // argh
              if (raw_geno_cts) {
                STD_ARRAY_REF(uint32_t, 3) cur_raw_geno_cts = raw_geno_cts[variant_uidx];
                cur_raw_geno_cts[0] = genocounts[0];
                cur_raw_geno_cts[1] = genocounts[1];
                cur_raw_geno_cts[2] = genocounts[2];
              }
              if (nonfemale_ct != chry_missingstat_sample_ct) {
                GetMCounts64(sex_male, sex_male_interleaved_vec, &pgv, raw_sample_ct, male_ct, allele_ct, genocounts, all_dosages, two_cts);
              }
              uintptr_t hethap_x2 = 0;
              for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
                hethap_x2 += all_dosages[aidx];
              }
              hethap_ct = hethap_x2 / 2;
            }
          } else {
            // chrX
            // multiallelic dosages not supported yet
            const PglErr reterr = PgrGetM(nullptr, pssi, raw_sample_ct, variant_uidx, pgrp, &pgv);
            if (unlikely(reterr)) {
              new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
              goto LoadAlleleAndGenoCountsThread_err;
            }
            ZeroTrailingNyps(raw_sample_ct, pgv.genovec);
            // We don't attempt to compute imp_r2 on chrX, so flat counts are
            // fine.
            GetMFlatCounts64(sample_include, sample_include_interleaved_vec, &pgv, raw_sample_ct, sample_ct, allele_ct, genocounts, all_dosages);

            // Double all counts, then subtract male counts.
            for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
              all_dosages[aidx] *= 2;
            }
            GenoarrCountSubsetFreqs(pgv.genovec, sex_male_interleaved_vec, raw_sample_ct, male_ct, sex_specific_genocounts);
            hethap_ct = sex_specific_genocounts[1];
            if (male_ct) {
              all_dosages[0] -= 2 * sex_specific_genocounts[0] + hethap_ct;

              // may underflow
              all_dosages[1] -= 2 * sex_specific_genocounts[2] + hethap_ct;

              if (pgv.patch_01_ct) {
                uintptr_t sample_widx = 0;
                uintptr_t patch_01_bits = pgv.patch_01_set[0];
                uint32_t male_patch_01_ct = 0;
                for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                  const uintptr_t lowbit = BitIter1y(pgv.patch_01_set, &sample_widx, &patch_01_bits);
                  if (sex_male[sample_widx] & lowbit) {
                    ++male_patch_01_ct;
                    all_dosages[pgv.patch_01_vals[uii]] -= 1;
                  }
                }
                all_dosages[1] += male_patch_01_ct;
              }
              if (pgv.patch_10_ct) {
                uintptr_t sample_widx = 0;
                uintptr_t patch_10_bits = pgv.patch_10_set[0];
                uint32_t male_patch_10_ct = 0;
                for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                  const uintptr_t lowbit = BitIter1y(pgv.patch_10_set, &sample_widx, &patch_10_bits);
                  if (sex_male[sample_widx] & lowbit) {
                    ++male_patch_10_ct;
                    const AlleleCode code_lo = pgv.patch_10_vals[2 * uii];
                    const AlleleCode code_hi = pgv.patch_10_vals[2 * uii + 1];
                    all_dosages[code_lo] -= 1;
                    all_dosages[code_hi] -= 1;
                    hethap_ct += (code_lo != code_hi);
                  }
                }
                all_dosages[1] += male_patch_10_ct * 2;
              }
            }
            if (allele_presents_bytearr) {
              for (uintptr_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                if (all_dosages[allele_idx]) {
                  allele_presents_bytearr[cur_allele_idx_offset + allele_idx] = 128;
                }
              }
            }
            if (allele_ddosages) {
              for (uintptr_t aidx = 0; aidx != allele_ct; ++aidx) {
                allele_ddosages[cur_allele_idx_offset + aidx] = all_dosages[aidx] * kDosageMid;
              }
            }
            if (x_male_geno_cts) {
              STD_ARRAY_REF(uint32_t, 3) cur_x_male_geno_cts = x_male_geno_cts[variant_uidx - x_start];
              cur_x_male_geno_cts[0] = sex_specific_genocounts[0];
              cur_x_male_geno_cts[1] = sex_specific_genocounts[1];
              cur_x_male_geno_cts[2] = sex_specific_genocounts[2];
              if (x_nosex_geno_cts) {
                GenoarrCountSubsetFreqs(pgv.genovec, nosex_interleaved_vec, raw_sample_ct, nosex_ct, sex_specific_genocounts);
                STD_ARRAY_REF(uint32_t, 3) cur_nosex_geno_cts = x_nosex_geno_cts[variant_uidx - x_start];
                cur_nosex_geno_cts[0] = sex_specific_genocounts[0];
                cur_nosex_geno_cts[1] = sex_specific_genocounts[1];
                cur_nosex_geno_cts[2] = sex_specific_genocounts[2];
              }
            }
          }
          if (variant_missing_dosage_cts) {
            // multiallelic dosage not supported yet
            variant_missing_dosage_cts[variant_uidx] = genocounts[3];
          }
        }
        if (raw_geno_cts && (!is_y)) {
          STD_ARRAY_REF(uint32_t, 3) cur_raw_geno_cts = raw_geno_cts[variant_uidx];
          cur_raw_geno_cts[0] = genocounts[0];
          cur_raw_geno_cts[1] = genocounts[1];
          cur_raw_geno_cts[2] = genocounts[2];
        }
        if (variant_missing_hc_cts) {
          variant_missing_hc_cts[variant_uidx] = genocounts[3];
          if (variant_hethap_cts && (variant_uidx >= first_hap_uidx)) {
            variant_hethap_cts[variant_uidx - first_hap_uidx] = hethap_ct;
          }
        }
      }
      if (++subset_idx == subset_ct) {
        break;
      }
      sample_include = ctx->founder_info;
      sample_include_interleaved_vec = ctx->founder_info_interleaved_vec;
      sample_include_cumulative_popcounts = ctx->founder_info_cumulative_popcounts;
      sex_male = ctx->founder_male;
      sex_male_interleaved_vec = ctx->founder_male_interleaved_vec;
      sex_nonfemale = ctx->founder_nonfemale;
      sex_nonfemale_interleaved_vec = ctx->founder_nonfemale_interleaved_vec;
      sex_nonfemale_cumulative_popcounts = ctx->founder_nonfemale_cumulative_popcounts;

      nosex_interleaved_vec = ctx->founder_nosex_interleaved_vec;

      sample_ct = ctx->founder_ct;
      male_ct = ctx->founder_male_ct;
      nosex_ct = ctx->founder_nosex_ct;
      nonfemale_ct = male_ct + nosex_ct;
      chry_missingstat_sample_ct = ctx->founder_chry_missingstat_sample_ct;
      allele_presents_bytearr = nullptr;
      allele_ddosages = ctx->founder_allele_ddosages;
      variant_missing_hc_cts = nullptr;
      variant_missing_dosage_cts = nullptr;
      raw_geno_cts = ctx->founder_raw_geno_cts;
      x_male_geno_cts = ctx->founder_x_male_geno_cts;
      x_nosex_geno_cts = ctx->founder_x_nosex_geno_cts;
      imp_r2_vals = nullptr;
    }
    while (0) {
    LoadAlleleAndGenoCountsThread_err:
      UpdateU64IfSmaller(new_err_info, &ctx->err_info);
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr LoadAlleleAndGenoCounts(const uintptr_t* sample_include, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uintptr_t* allele_idx_offsets, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t founder_ct, uint32_t male_ct, uint32_t nosex_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t first_hap_uidx, uint32_t is_minimac3_r2, uint32_t y_nosex_missing_stats, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, uintptr_t* allele_presents, uint64_t* allele_ddosages, uint64_t* founder_allele_ddosages, uint32_t* variant_missing_hc_cts, uint32_t* variant_missing_dosage_cts, uint32_t* variant_hethap_cts, STD_ARRAY_PTR_DECL(uint32_t, 3, raw_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, x_male_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, x_nosex_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts), double* imp_r2_vals) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  LoadAlleleAndGenoCountsCtx ctx;
  {
    if (!variant_ct) {
      goto LoadAlleleAndGenoCounts_ret_1;
    }

    // four cases:
    // 1. allele_ddosages, raw_geno_cts, and/or variant_missing_{hc,dosage}_cts
    //    required, and that's it
    // 2. founder_allele_ddosages and/or founder_raw_geno_cts required, and
    //    that's it
    // 3. both required, and founder_ct != sample_ct.
    // 4. both required, and founder_ct == sample_ct.  caller is expected to
    //    make founder_allele_ddosages and allele_ddosages point to the same
    //    memory, ditto for founder_raw_geno_cts/raw_geno_cts.
    const uint32_t only_founder_cts_required = (!allele_presents) && (!allele_ddosages) && (!raw_geno_cts) && (!variant_missing_hc_cts) && (!variant_missing_dosage_cts);
    const uint32_t two_subsets_required = (founder_ct != sample_ct) && (!only_founder_cts_required) && (founder_allele_ddosages || founder_raw_geno_cts);
    ctx.cip = cip;
    ctx.sample_include = only_founder_cts_required? founder_info : sample_include;
    ctx.raw_sample_ct = raw_sample_ct;
    ctx.sample_ct = only_founder_cts_required? founder_ct : sample_ct;
    ctx.allele_ddosages = only_founder_cts_required? founder_allele_ddosages : allele_ddosages;
    ctx.raw_geno_cts = only_founder_cts_required? founder_raw_geno_cts : raw_geno_cts;
    ctx.x_male_geno_cts = only_founder_cts_required? founder_x_male_geno_cts : x_male_geno_cts;
    ctx.x_nosex_geno_cts = only_founder_cts_required? founder_x_nosex_geno_cts : x_nosex_geno_cts;
    ctx.imp_r2_vals = imp_r2_vals;
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t raw_sample_ctv = BitCtToVecCt(raw_sample_ct);
    if (unlikely(bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &ctx.sample_include_interleaved_vec) ||
                 bigstack_alloc_u32(raw_sample_ctl, &ctx.sample_include_cumulative_popcounts) ||
                 bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &ctx.sex_male_interleaved_vec) ||
                 bigstack_alloc_u32(raw_sample_ctl, &ctx.sex_nonfemale_cumulative_popcounts))) {
      goto LoadAlleleAndGenoCounts_ret_NOMEM;
    }
    FillInterleavedMaskVec(ctx.sample_include, raw_sample_ctv, ctx.sample_include_interleaved_vec);
    FillCumulativePopcounts(ctx.sample_include, raw_sample_ctl, ctx.sample_include_cumulative_popcounts);
    if ((founder_ct == sample_ct) || (!only_founder_cts_required)) {
      ctx.sex_male = sex_male;
    } else {
      // no nonfounder counts required
      uintptr_t* new_sex_male;
      if (unlikely(bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &new_sex_male))) {
        goto LoadAlleleAndGenoCounts_ret_NOMEM;
      }
      BitvecAndCopy(sex_male, founder_info, raw_sample_ctl, new_sex_male);
      ZeroTrailingWords(raw_sample_ctl, new_sex_male);
      ctx.sex_male = new_sex_male;
    }
    FillInterleavedMaskVec(ctx.sex_male, raw_sample_ctv, ctx.sex_male_interleaved_vec);
    // slightly wasteful, but founder_male_ct isn't precomputed
    ctx.male_ct = PopcountWords(ctx.sex_male, raw_sample_ctl);
    ctx.nosex_ct = 0;
    if (nosex_ct) {
      uintptr_t* sex_nonfemale;
      if (unlikely(bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &sex_nonfemale) ||
                   bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &ctx.sex_nonfemale_interleaved_vec))) {
        goto LoadAlleleAndGenoCounts_ret_NOMEM;
      }
      if ((founder_ct == sample_ct) || (!only_founder_cts_required)) {
        // bugfix (30 May 2025): sex_nonfemale initialization was incorrect
        // when sample_ct != raw_sample_ct
        BitvecXor3Copy(ctx.sample_include, ctx.sex_male, sex_nm, raw_sample_ctl, sex_nonfemale);
      } else {
        BitvecInvmaskCopy(founder_info, sex_nm, raw_sample_ctl, sex_nonfemale);
        BitvecOr(ctx.sex_male, raw_sample_ctl, sex_nonfemale);
      }
      ZeroTrailingWords(raw_sample_ctl, sex_nonfemale);
      ctx.sex_nonfemale = sex_nonfemale;
      const uint32_t nonfemale_ct = PopcountWords(sex_nonfemale, raw_sample_ctl);
      ctx.nosex_ct = nonfemale_ct - ctx.male_ct;
      FillInterleavedMaskVec(sex_nonfemale, raw_sample_ctv, ctx.sex_nonfemale_interleaved_vec);
    } else {
      ctx.sex_nonfemale = ctx.sex_male;
      ctx.sex_nonfemale_interleaved_vec = ctx.sex_male_interleaved_vec;
    }
    FillCumulativePopcounts(ctx.sex_nonfemale, raw_sample_ctl, ctx.sex_nonfemale_cumulative_popcounts);
    ctx.chry_missingstat_sample_ct = y_nosex_missing_stats? (ctx.male_ct + ctx.nosex_ct) : ctx.male_ct;
    ctx.nosex_interleaved_vec = nullptr;
    uintptr_t* nosex_buf = nullptr;
    if (x_nosex_geno_cts || founder_x_nosex_geno_cts) {
      if (unlikely(bigstack_end_alloc_w(raw_sample_ctl, &nosex_buf) ||
                   bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &ctx.nosex_interleaved_vec))) {
        goto LoadAlleleAndGenoCounts_ret_NOMEM;
      }
      BitvecInvmaskCopy(ctx.sample_include, sex_nm, raw_sample_ctl, nosex_buf);
      ZeroTrailingWords(raw_sample_ctl, nosex_buf);
      FillInterleavedMaskVec(nosex_buf, raw_sample_ctv, ctx.nosex_interleaved_vec);
    }

    ctx.variant_missing_hc_cts = variant_missing_hc_cts;
    ctx.variant_missing_dosage_cts = variant_missing_dosage_cts;
    ctx.variant_hethap_cts = variant_hethap_cts;
    ctx.first_hap_uidx = first_hap_uidx;
    ctx.is_minimac3_r2 = is_minimac3_r2;

    ctx.founder_info = nullptr;
    ctx.founder_info_interleaved_vec = nullptr;
    ctx.founder_info_cumulative_popcounts = nullptr;
    ctx.founder_male = nullptr;
    ctx.founder_male_interleaved_vec = nullptr;
    ctx.founder_nonfemale = nullptr;
    ctx.founder_nonfemale_interleaved_vec = nullptr;
    ctx.founder_nonfemale_cumulative_popcounts = nullptr;
    ctx.founder_nosex_interleaved_vec = nullptr;
    ctx.founder_ct = 0;
    ctx.founder_male_ct = 0;
    ctx.founder_nosex_ct = 0;
    ctx.founder_chry_missingstat_sample_ct = 0;
    ctx.founder_allele_ddosages = nullptr;
    ctx.founder_raw_geno_cts = nullptr;
    ctx.founder_x_male_geno_cts = nullptr;
    ctx.founder_x_nosex_geno_cts = nullptr;
    if (two_subsets_required) {
      if (founder_ct) {
        ctx.founder_info = founder_info;
        if (unlikely(bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &ctx.founder_info_interleaved_vec) ||
                     bigstack_alloc_u32(raw_sample_ctl, &ctx.founder_info_cumulative_popcounts) ||
                     bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &ctx.founder_male) ||
                     bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &ctx.founder_male_interleaved_vec) ||
                     bigstack_alloc_u32(raw_sample_ctl, &ctx.founder_nonfemale_cumulative_popcounts))) {
          goto LoadAlleleAndGenoCounts_ret_NOMEM;
        }
        FillInterleavedMaskVec(founder_info, raw_sample_ctv, ctx.founder_info_interleaved_vec);
        FillCumulativePopcounts(founder_info, raw_sample_ctl, ctx.founder_info_cumulative_popcounts);
        BitvecAndCopy(sex_male, founder_info, raw_sample_ctl, ctx.founder_male);
        ZeroTrailingWords(raw_sample_ctl, ctx.founder_male);
        FillInterleavedMaskVec(ctx.founder_male, raw_sample_ctv, ctx.founder_male_interleaved_vec);
        ctx.founder_male_ct = PopcountWords(ctx.founder_male, raw_sample_ctl);
        if (ctx.nosex_ct) {
          if (unlikely(bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &ctx.founder_nonfemale) ||
                       bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &ctx.founder_nonfemale_interleaved_vec))) {
            goto LoadAlleleAndGenoCounts_ret_NOMEM;
          }
          BitvecInvmaskCopy(founder_info, sex_nm, raw_sample_ctl, ctx.founder_nonfemale);
          BitvecOr(ctx.founder_male, raw_sample_ctl, ctx.founder_nonfemale);
          ZeroTrailingWords(raw_sample_ctl, ctx.founder_nonfemale);
          const uint32_t founder_nonfemale_ct = PopcountWords(ctx.founder_nonfemale, raw_sample_ctl);
          ctx.founder_nosex_ct = founder_nonfemale_ct - ctx.founder_male_ct;
          FillInterleavedMaskVec(ctx.founder_nonfemale, raw_sample_ctv, ctx.founder_nonfemale_interleaved_vec);
        } else {
          ctx.founder_nonfemale = ctx.founder_male;
          ctx.founder_nonfemale_interleaved_vec = ctx.founder_male_interleaved_vec;
        }
        FillCumulativePopcounts(ctx.founder_nonfemale, raw_sample_ctl, ctx.founder_nonfemale_cumulative_popcounts);
        ctx.founder_chry_missingstat_sample_ct = y_nosex_missing_stats? (ctx.founder_male_ct + ctx.founder_nosex_ct) : ctx.founder_male_ct;
        if (ctx.founder_nosex_ct) {
          // caller currently responsible for ensuring that when
          // founder_nosex_ct is zero, founder_x_nosex_geno_cts == nullptr
          if (founder_x_nosex_geno_cts) {
            if (unlikely(bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &ctx.founder_nosex_interleaved_vec))) {
              goto LoadAlleleAndGenoCounts_ret_NOMEM;
            }
            BitvecAnd(founder_info, raw_sample_ctl, nosex_buf);
            ZeroTrailingWords(raw_sample_ctl, nosex_buf);
            FillInterleavedMaskVec(nosex_buf, raw_sample_ctv, ctx.founder_nosex_interleaved_vec);
            ctx.founder_x_nosex_geno_cts = founder_x_nosex_geno_cts;
          }
        }
        ctx.founder_ct = founder_ct;
        ctx.founder_allele_ddosages = founder_allele_ddosages;
        ctx.founder_raw_geno_cts = founder_raw_geno_cts;
        ctx.founder_x_male_geno_cts = founder_x_male_geno_cts;
      } else {
        if (founder_allele_ddosages) {
          ZeroU64Arr(allele_idx_offsets? allele_idx_offsets[raw_variant_ct] : (2 * raw_variant_ct), founder_allele_ddosages);
        }
        if (founder_raw_geno_cts) {
          memset(founder_raw_geno_cts, 0, raw_variant_ct * (3 * sizeof(int32_t)));
        }
      }
    } else if (founder_ct == sample_ct) {
      // bugfix: some founder and some nonfounder counts required
      if ((!ctx.allele_ddosages) && founder_allele_ddosages) {
        ctx.allele_ddosages = founder_allele_ddosages;
      }
      if ((!ctx.raw_geno_cts) && founder_raw_geno_cts) {
        ctx.raw_geno_cts = founder_raw_geno_cts;
      }
      if ((!ctx.x_male_geno_cts) && founder_x_male_geno_cts) {
        ctx.x_male_geno_cts = founder_x_male_geno_cts;
      }
      if ((!ctx.x_nosex_geno_cts) && founder_x_nosex_geno_cts) {
        ctx.x_nosex_geno_cts = founder_x_nosex_geno_cts;
      }
    }
    const uintptr_t raw_allele_ct = allele_idx_offsets? allele_idx_offsets[raw_variant_ct] : (2 * raw_variant_ct);
    if (!ctx.sample_ct) {
      if (allele_presents) {
        ZeroWArr(BitCtToWordCt(raw_allele_ct), allele_presents);
      }
      if (ctx.allele_ddosages) {
        ZeroU64Arr(raw_allele_ct, ctx.allele_ddosages);
      }
      if (ctx.raw_geno_cts) {
        memset(ctx.raw_geno_cts, 0, raw_variant_ct * (3 * sizeof(int32_t)));
      }
      // early exit
      goto LoadAlleleAndGenoCounts_ret_1;
    }
    BigstackEndReset(bigstack_end_mark);  // free nosex_buf
    if (allele_presents) {
      const uintptr_t raw_allele_ct_a64 = RoundUpPow2(raw_allele_ct, kCacheline);
      if (unlikely(bigstack_left() < raw_allele_ct_a64)) {
        goto LoadAlleleAndGenoCounts_ret_NOMEM;
      }
      // fill byte-array instead of bitarray so multithreading works
      ctx.allele_presents_bytearr = S_CAST(unsigned char*, bigstack_alloc_raw(raw_allele_ct_a64));
      memset(ctx.allele_presents_bytearr, 0, raw_allele_ct_a64);
    } else {
      ctx.allele_presents_bytearr = nullptr;
    }

    uint32_t unused_chr_code;
    uint32_t unused_chr_code2;
    uint32_t xy_complications_present = ((allele_presents || allele_ddosages || founder_allele_ddosages || variant_missing_dosage_cts) && XymtExists(cip, kChrOffsetX, &unused_chr_code)) || (((ctx.chry_missingstat_sample_ct != ctx.male_ct + ctx.nosex_ct) || (allele_presents && (sample_ct != male_ct))) && XymtExists(cip, kChrOffsetY, &unused_chr_code2));
    const uint32_t xy_dosages_needed = (pgfip->gflags & kfPgenGlobalDosagePresent) && xy_complications_present;

    // todo: check when this saturates
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    const uint32_t max_allele_ct = pgfip->max_allele_ct;
    uint32_t mhc_needed = 0;
    ctx.thread_read_mhc = nullptr;
    if ((max_allele_ct > 2) && (variant_hethap_cts || allele_presents || allele_ddosages || founder_allele_ddosages || imp_r2_vals)) {
      if (unlikely(bigstack_alloc_u64p(calc_thread_ct, &ctx.all_dosages))) {
        goto LoadAlleleAndGenoCounts_ret_NOMEM;
      }
      mhc_needed = (xy_complications_present || ((variant_hethap_cts || imp_r2_vals) && XymtExists(cip, kChrOffsetX, &unused_chr_code)));
      for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
        // double allocation size, to leave room for chrY ssqs
        if (unlikely(bigstack_alloc_u64(max_allele_ct * 2, &(ctx.all_dosages[tidx])))) {
          goto LoadAlleleAndGenoCounts_ret_NOMEM;
        }
      }
    } else {
      ctx.all_dosages = nullptr;
    }
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    // defensive
    ctx.dosage_presents = nullptr;
    ctx.dosage_mains = nullptr;
    uint32_t read_block_size;
    // todo: check if raw_sample_ct should be replaced with sample_ct here
    if (unlikely(PgenMtLoadInit(variant_include, raw_sample_ct, variant_ct, bigstack_left(), pgr_alloc_cacheline_ct, 0, 0, 0, pgfip, &calc_thread_ct, &ctx.genovecs, mhc_needed? (&ctx.thread_read_mhc) : nullptr, nullptr, nullptr, xy_dosages_needed? (&ctx.dosage_presents) : nullptr, xy_dosages_needed? (&ctx.dosage_mains) : nullptr, nullptr, nullptr, &read_block_size, nullptr, main_loadbufs, &ctx.pgr_ptrs, &ctx.read_variant_uidx_starts))) {
      goto LoadAlleleAndGenoCounts_ret_NOMEM;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto LoadAlleleAndGenoCounts_ret_NOMEM;
    }
    ctx.variant_include = variant_include;
    ctx.allele_idx_offsets = allele_idx_offsets;
    ctx.err_info = (~0LLU) << 32;
    SetThreadFuncAndData(LoadAlleleAndGenoCountsThread, &ctx, &tg);

    logputs("Calculating allele frequencies... ");
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t pct = 0;

    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
    for (uint32_t variant_idx = 0; ; ) {
      const uint32_t cur_block_size = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
      if (unlikely(reterr)) {
        goto LoadAlleleAndGenoCounts_ret_PGR_FAIL;
      }
      if (variant_idx) {
        JoinThreads(&tg);
        reterr = S_CAST(PglErr, ctx.err_info);
        if (unlikely(reterr)) {
          PgenErrPrintNV(reterr, ctx.err_info >> 32);
          goto LoadAlleleAndGenoCounts_ret_1;
        }
      }
      if (!IsLastBlock(&tg)) {
        ctx.cur_block_size = cur_block_size;
        ComputeUidxStartPartition(variant_include, cur_block_size, calc_thread_ct, read_block_idx * read_block_size, ctx.read_variant_uidx_starts);
        PgrCopyBaseAndOffset(pgfip, calc_thread_ct, ctx.pgr_ptrs);
        if (variant_idx + cur_block_size == variant_ct) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto LoadAlleleAndGenoCounts_ret_THREAD_CREATE_FAIL;
        }
      }

      parity = 1 - parity;
      if (variant_idx == variant_ct) {
        break;
      }
      if (variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
      }

      ++read_block_idx;
      variant_idx += cur_block_size;
      // crucially, this is independent of the PgenReader block_base
      // pointers
      pgfip->block_base = main_loadbufs[parity];
    }
    if (allele_presents) {
      const uintptr_t raw_allele_ctl = BitCtToWordCt(raw_allele_ct);
      allele_presents[raw_allele_ctl - 1] = 0;
#ifdef __LP64__
      // todo: better ARM implementation
      const uintptr_t vec_ct = DivUp(raw_allele_ct, kBytesPerVec);
      VecUc* bytearr_alias = R_CAST(VecUc*, ctx.allele_presents_bytearr);
      Vec8thUint* allele_presents_alias = R_CAST(Vec8thUint*, allele_presents);
      for (uintptr_t vec_idx = 0; vec_idx != vec_ct; ++vec_idx) {
        allele_presents_alias[vec_idx] = vecuc_movemask(bytearr_alias[vec_idx]);
      }
#else
      const uintptr_t twovec_ct = DivUp(raw_allele_ct, 8);
      uintptr_t* bytearr_iter = R_CAST(uintptr_t*, ctx.allele_presents_bytearr);
      unsigned char* allele_presents_iter = R_CAST(unsigned char*, allele_presents);
      unsigned char* allele_presents_stop = &(allele_presents_iter[twovec_ct]);
      for (; allele_presents_iter != allele_presents_stop; ++allele_presents_iter) {
        // 31,23,15,7 -> 3,2,1,0: multiply by number with bits 0,7,14,21 set,
        // then right-shift
        uintptr_t cur_word = ((*bytearr_iter++) * 0x204081) >> 28;
        cur_word |= ((*bytearr_iter++) * 0x204081) >> 24;
        *allele_presents_iter = cur_word;
      }
#endif
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
  }
  while (0) {
  LoadAlleleAndGenoCounts_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LoadAlleleAndGenoCounts_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  LoadAlleleAndGenoCounts_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 LoadAlleleAndGenoCounts_ret_1:
  CleanupThreads(&tg);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  pgfip->block_base = nullptr;
  return reterr;
}
// Join behavior:
// - Require sorted input .pvar for now, though it won't be difficult to lift
//   this restriction later.  Ok for input to contain multiallelic variants.
// - Don't need to do anything different when chr:pos only appears once in a
//   biallelic variant, or in a multiallelic variant in "+any" mode.  For
//   multiallelic variants in +both/+snps mode, error out if the variant is
//   mixed SNP/non-SNP.  In +both case, also error out if the variant is mixed
//   symbolic/non-symbolic, or it's symbolic and does not satisfy the REF or
//   INFO/END constraints.
// - When multiple variants have the same chr:pos:
//   - In +snps mode, also don't need to do anything different for non-SNPs.
//   - Otherwise, create up to three linked lists of input variant records
//     (need to distinguish SNPs from non-SNPs in +both mode, and within the
//     non-SNP category, symbolic alleles are separate from non-symbolic).
//     The variant with symbolic alleles has additional constraints: a warning
//     is printed if INFO/END isn't defined, and an error occurs if either REF
//     is multi-character, or there's an INFO/END mismatch.
//   - Error out if REF alleles aren't all consistent, or any ALT allele is
//     duplicated (note that --pmerge must support the latter).
//   - For joined not-entirely-SNP non-symbolic variants, the final REF is the
//     longest of the original REFs; ALT alleles have bases added to the end if
//     necessary.  (Yes, this causes SNPs to stop being visible to a strlen ==
//     1 check in +any mode, which is why + is interpreted as +both instead.)
//   - Final ALT allele order is based on allele frequency (highest first),
//     with ties broken by natural-sort.
//   - ID: 1. If --set-all-var-ids specified, apply template.
//         2. Otherwise, keep original variant ID if all sources identical and
//            nonmissing.
//         3. Otherwise, if vid-join specified, check if all original IDs are
//            nonmissing and contain exactly one ';' per extra ALT allele.  If
//            so, join on ';' in final ALT allele order.
//         4. Otherwise, set to --set-missing-var-ids template or missing code.
//   - QUAL is minimum of inputs ('.' treated as positive infinity).  FILTER is
//     natural-sorted union of non-PASS values; if all inputs were '.'/PASS,
//     output is PASS unless all inputs were missing.
//   - INFO join is based on the Number field in the key's header line.  (If
//     there's no header line corresponding to a key, we error out.)  For
//     Number=A, we join in the obvious manner.  For Number=R (or Number=G in
//     the haploid case), we replace the reference allele entry with . iff
//     there's any string mismatch (e.g. '13' and '13.0' will be treated as
//     unequal), and print a warning (with more than 3 warnings, later warnings
//     are only written to log file).  For diploid Number=G, we do the same for
//     the hom-ref entry.  For Number=0, we error out if the header line Type
//     isn't Flag, and the final flag is set iff any of the original variants
//     have the flag set.  For the other cases (Number=<fixed constant> or
//     '.'), we replace the value with '.' if there's any string mismatch.
//     The key won't appear at all iff it doesn't appear in any of the original
//     variants (not even a '.').
//   - We error out if the joined variant has total ALT dosage > ~2.02, or
//     either side of the total ALT phased dosage > ~1.01.  (We scale the
//     components down if there's a <1% overflow.)
// - Missing ALTs... ugh.  Have to permit this in biallelic case, but don't
//   allow it elsewhere.  This forces a SNP (if REF is single-char) and/or
//   non-SNP (if REF is multichar) entry to be written when no corresponding
//   regular ALT is present, genotype writing/merging will error out if a
//   missing allele has any dosage, and QUAL/FILTER/INFO is merged as usual
//   when there are multiple missing-ALT same-type variants at the same
//   position.  However, when a regular ALT is present, the missing-ALT
//   variants are completely ignored.

// The simplest design involves precomputing the entire (new variant idx, new
// allele idx) -> (old variant uidx, old allele idx) mapping here, and
// referring to that in both the .pvar and .pgen writers.  Unfortunately, that
// has a rather high memory requirement of 5 bytes per allele (assuming
// sizeof(AlleleCode) == 1).  While that's smaller than the 8 bytes/allele we
// pay for allele_storage[], and also practically always smaller than the
// 21+[variant ID len] per variant we pay for variant_bps + variant_ids +
// allele_idx_offsets, it's still worth some effort to avoid; in particular, we
// want an 8 GiB workspace to be sufficient for most operations on the full
// 1000 Genomes phase 3 variant set (~84.8 million), and that's barely true
// right now, so a bit of additional complexity to avoid losing ~850 MB is
// justified.
//
// Thus, we only save the number of alleles in each new variant here, and force
// the which-allele-comes-from-where computation (as well as SNP vs. non-SNP
// vs. symbolic handling) to be repeated in the .pvar and .pgen writers.  This
// sucks, but being forced to perform an ordinary analysis on a remote machine
// rather than locally sucks a bit more.
//
// Incidentally, another place to look, if it's important to further reduce
// memory requirements, is internal representation of the FILTER field.  The
// current design gains a bit of extra speed by simply storing the non-./PASS
// strings without parsing them further; but we could put them into a temporary
// storage location and then convert to bitarray + string table at the end of
// LoadPvar().  (Note that we already use only bitarrays when all FILTER values
// are ./PASS, though.)
// *write_allele_idx_offsetsp is assumed to be initialized to nullptr.
// *max_missalt_ctp is assumed to be initialized to 0.
// Returns 1 iff there are exactly (allele_ct - 2) semicolons in
// orig_variant_id, and no two are adjacent (or leading/trailing).
// Similar to WriteMapOrBim(), but there are enough small differences to
// justify making this a separate function instead of clogging the original
// with more conditionals.
// info_keys[] entries point to the (variable-size) key[] member of InfoVtype
// structs.  We use [const_]container_of(x)->num to look up the associated
// Number= value.
// Assumes info_has_g_key initialized to 0.
// Sets INFO/PR Flag num to kInfoVtypeSkip.
// Final filter_keys is natural-sorted.
// Return values are allocated on bottom of bigstack.
// Caller must initialize all return values to correspond to the null table.
/*
*/
// Returned arrays are allocated at end of bigstack.
// possible todo: use more compact representation.  PLINK 1.9's setdefs are one
// possibility, there are others.
// Returned arrays are allocated at bottom of bigstack.
// May set extraneous phaseinfo bits.
// todo: add .pgen-specific stuff
// initialized mcp fields: cip, fip, sex_male_collapsed,
// sex_male_collapsed_interleaved, sex_female_collapsed_interleaved,
// write_allele_permute, flip_subset_*, raw_sample_ct, sample_ct, male_ct,
// female_ct, plink2_write_flags
// One-thread-per-vblock is sensible for possibly-phased biallelic data, where
// subsetting and LD-compression are a substantial fraction of processing time,
// and memory requirements tend to be low enough that it's actually reasonable
// for each thread job to comprise 64k variants.
// Beyond that... the VCF/.pgen division of labor looks nice, but far too much
// of the work is usually being done in the initial PgrGetRaw() call, so just
// fall back on single-threaded invocation of the same function; only
// difference is that the worker thread owns the writer object.
//
// Possible todo: it doesn't come up frequently, but sorting a .pgen that is
// thoroughly out of order can be very slow due to all the nonsequential I/O.
// That could be accelerated with the pvar-INFO temporary-file strategy (split
// into remaining-workspace-sized original-order shard files, then load one
// entire shard at a time and write it out in the correct order).
// May want to benchmark on Linux first to see if this is a real problem; the
// horrible benchmark result motivating this comment was on macOS, which is
// known to have a sketchy filesystem.
// Single-output-thread implementation.  Allows variants to be unsorted.
// (Note that MakePlink2NoVsort() currently requires enough memory for 64k * 2
// variants per output thread, due to LD compression.  This is faster in the
// common case, but once you have 150k+ samples with dosage data...)
//
// initialized mcp fields: cip, fip, sex_male_collapsed,
// sex_male_collapsed_interleaved, sex_female_collapsed_interleaved,
// write_allele_permute, flip_subset_*, raw_sample_ct, sample_ct,
// plink2_write_flags
// hybrid of WriteMapOrBim() and write_pvar_resorted()
// allele_presents must be nullptr unless we're trimming alt alleles.
//
// The annoying part of this is handling a sequence of INFO strings that don't
// fit in memory; we use a multipass approach for that.  File creation,
// allocation of buffers, and generating the header line occurs directly in
// this function, while loading the next pvar_info_strs batch and writing the
// next .pvar line batch are one level down.
#ifdef __cplusplus
}  // namespace plink2
#endif
