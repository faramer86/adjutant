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

#include "adjutant_matrix_calc.h"

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "include/plink2_bits.h"
#include "include/plink2_float.h"
#include "include/plink2_htable.h"
#include "include/plink2_simd.h"
#include "include/plink2_string.h"
#include "include/plink2_text.h"
#include "include/plink2_thread.h"
#include "adjutant_compress_stream.h"
#include "adjutant_decompress.h"
#include "adjutant_matrix.h"

#ifdef USE_CUDA
#  include "cuda/plink2_matrix_cuda.h"
#endif

#include <unistd.h>  // unlink()

#ifdef __cplusplus
namespace plink2 {
#endif

void InitScore(ScoreInfo* score_info_ptr) {
  score_info_ptr->flags = kfScore0;
  score_info_ptr->varid_col_p1 = 1;
  score_info_ptr->allele_col_p1 = 0;  // defensive
  score_info_ptr->input_fname = nullptr;
  InitRangeList(&(score_info_ptr->input_col_idx_range_list));

  score_info_ptr->qsr_range_fname = nullptr;
  score_info_ptr->qsr_data_fname = nullptr;
  score_info_ptr->qsr_varid_col_p1 = 1;
  score_info_ptr->qsr_val_col_p1 = 0;  // defensive

  score_info_ptr->partition_fname = nullptr;
}

void CleanupScore(ScoreInfo* score_info_ptr) {
  free_cond(score_info_ptr->input_fname);
  CleanupRangeList(&(score_info_ptr->input_col_idx_range_list));

  free_cond(score_info_ptr->qsr_range_fname);
  free_cond(score_info_ptr->qsr_data_fname);
  free_cond(score_info_ptr->partition_fname);
}
// Cost function for thread load-balancing: (max - min) * ((max + min)/2)
// i.e. we don't waste any time processing entries in the irrelevant
// upper-right triangle
// set align to 1 for no alignment
// Cost function for thread load-balancing: (max - min) * max
// i.e. when we can't avoid processing entries in the irrelevant upper-right
// triangle
// Returns 0 if cells_avail is insufficient.
#ifdef USE_SSE42
CONSTI32(kKingMultiplex, 1024);
CONSTI32(kKingMultiplexWords, kKingMultiplex / kBitsPerWord);
#else  // !USE_SSE42
#  ifdef __LP64__
CONSTI32(kKingMultiplex, 1536);
#  else
CONSTI32(kKingMultiplex, 960);
#  endif
static_assert(kKingMultiplex % (3 * kBitsPerVec) == 0, "Invalid kKingMultiplex value.");
CONSTI32(kKingMultiplexWords, kKingMultiplex / kBitsPerWord);
// expensive PopcountWord().  Use Lauradoux/Walisch accumulators, since
// Harley-Seal requires too many variables.
#endif
static_assert(!(kKingMultiplexWords % 2), "kKingMultiplexWords must be even for safe bit-transpose.");
/*
*/

// '2' refers to the larger index here
// could also return pointer to end?
#ifdef USE_SSE42
#else
#endif
// assumes trailing bits of genovec are zeroed out
// This breaks the "don't pass pssi between functions" rule since it's a thin
// wrapper around PgrGetInv1D().
double ComputeDiploidMultiallelicVariance(const double* cur_allele_freqs, uint32_t cur_allele_ct) {
  const uint32_t cur_allele_ct_m1 = cur_allele_ct - 1;
  double variance = 0.0;
  double freq_sum = 0.0;
  for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct_m1; ++allele_idx) {
    const double cur_allele_freq = cur_allele_freqs[allele_idx];
    variance += cur_allele_freq * (1.0 - cur_allele_freq);
    freq_sum += cur_allele_freq;
  }
  if (freq_sum < 1.0 - kSmallEpsilon) {
    const double last_allele_freq = 1.0 - freq_sum;
    variance += freq_sum * last_allele_freq;
  }
  return variance;
}

// Assumes trailing bits of pgvp->genovec have been zeroed out.
// This function handles error-message-logging.
// turns out dsyrk_ does exactly what we want here
// can't use dsyrk_, so we manually partition the GRM piece we need to compute
// into an appropriate number of sub-pieces
// should be able to remove NOLAPACK later since we already have a non-LAPACK
// SVD implementation
#ifndef NOLAPACK
// This seems to be better than 256 (due to avoidance of cache critical
// stride?)
// (still want this to be a multiple of 8, for cleaner multithreading)
// TODO: try to allow run to proceed with a smaller value when 240 induces an
// out-of-memory error.  It's cutting it close for UK Biobank.  Though we don't
// want to reduce the default value, I've tested that and it degrades
// performance.
#endif

PglErr ScoreListLoad(const char* list_fname, TextStream* score_txsp, LlStr** fname_iterp, uint32_t* infile_ctp) {
  PglErr reterr = kPglRetSuccess;
  {
    LlStr** llstr_endp = fname_iterp;
    uintptr_t infile_ct = 0;
    while (1) {
      const char* fname = TextGet(score_txsp);
      if (!fname) {
        break;
      }
      ++infile_ct;
      const char* fname_end = CurTokenEnd(fname);
      const uint32_t fname_slen = fname_end - fname;
      LlStr* new_entry;
      if (unlikely(bigstack_end_alloc_llstr(fname_slen + 1, &new_entry))) {
        goto ScoreListLoad_ret_NOMEM;
      }
      new_entry->next = nullptr;
      memcpyx(new_entry->str, fname, fname_slen, '\0');
      *llstr_endp = new_entry;
      llstr_endp = &(new_entry->next);
    }
    if (unlikely(!infile_ct)) {
      logerrprintfww("Error: --score-list: %s is empty.\n", list_fname);
      goto ScoreListLoad_ret_MALFORMED_INPUT;
    }
#ifdef __LP64__
    if (unlikely(infile_ct > 0xffffffffU)) {
      logerrprintfww("Error: --score-list: Too many entries in %s.\n", list_fname);
      reterr = kPglRetNotYetSupported;
      goto ScoreListLoad_ret_1;
    }
#endif
    reterr = TextRetarget((*fname_iterp)->str, score_txsp);
    if (unlikely(reterr)) {
      goto ScoreListLoad_ret_TSTREAM_FAIL;
    }
    *infile_ctp = infile_ct;
  }
  while (0) {
  ScoreListLoad_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ScoreListLoad_ret_TSTREAM_FAIL:
    TextStreamErrPrint(list_fname, score_txsp);
    break;
  ScoreListLoad_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
#ifdef __LP64__
 ScoreListLoad_ret_1:
#endif
  return reterr;
}

// to test: do we actually want cur_dosage_ints to be uint64_t* instead of
// uint32_t*?
void FillDdosageInts(const uintptr_t* genovec_buf, const uintptr_t* dosage_present, const Dosage* dosage_main_buf, uint32_t sample_ct, uint32_t dosage_ct, uint32_t is_diploid, uint64_t* cur_ddosage_ints) {
  uint64_t lookup_table[32] ALIGNV16;
  lookup_table[0] = 0;
  lookup_table[2] = (kDosageMid * k1LU) << is_diploid;
  lookup_table[4] = (kDosageMax * k1LU) << is_diploid;
  lookup_table[6] = 0;
  InitLookup16x8bx2(lookup_table);
  GenoarrLookup16x8bx2(genovec_buf, lookup_table, sample_ct, cur_ddosage_ints);
  if (!dosage_ct) {
    return;
  }
  uintptr_t sample_idx_base = 0;
  uintptr_t cur_bits = dosage_present[0];
  for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
    const uintptr_t sample_idx = BitIter1(dosage_present, &sample_idx_base, &cur_bits);
    cur_ddosage_ints[sample_idx] = (dosage_main_buf[dosage_idx] * k1LU) << is_diploid;
  }
}

CONSTI32(kScoreVariantBlockSize, 240);

// must be a power of 2, >= kBitsPerWord
#ifdef USE_AVX2
CONSTI32(kScoreSampleShardSizeAlign, 256);
#else
CONSTI32(kScoreSampleShardSizeAlign, 128);
#endif

// The previous implementation had poor parallelism; the main thread was
// responsible for too many per-sample uint64_t and floating-point operations.
// We now shard those operations across all worker threads, and leave only
// much faster bitvector and similar operations in the main thread.
typedef struct CalcScoreCtxStruct {
  const uintptr_t* variant_include;
  const uint32_t* variant_include_cumulative_popcounts;
  const uintptr_t* qsr_include;
  const uintptr_t* sex_nonmale_collapsed;
  const uintptr_t* sex_female_collapsed;
  uint32_t score_final_col_ct;
  uint32_t sample_shard_size;
  uint32_t sample_ct;
  uint32_t max_difflist_len;
  ScoreFlags flags;
  uint32_t qsr_ct;

  uint32_t* variant_uidxs[2];
  uintptr_t* genovecs[2];
  uintptr_t* raregenos[2];
  uint32_t* difflist_lens[2];
  uint32_t* difflist_sample_ids[2];
  int8_t* difflist_common_genos[2];
  uintptr_t* dosage_presents[2];
  Dosage* dosage_mains[2];
  uint32_t* dosage_cts[2];
  uintptr_t* missing_bitvecs[2];
  uintptr_t* missing_male_bitvecs[2];
  uintptr_t* missing_nonfemale_bitvecs[2];
  double* allele_freqs[2];
  double* geno_slopes[2];
  double* geno_intercepts[2];
  double* score_dense_coefs_cmaj[2];
  double* score_sparse_coefs_vmaj[2];
  unsigned char is_nonx_haploids[2][kScoreVariantBlockSize];
  unsigned char is_relevant_xs[2][kScoreVariantBlockSize];
  unsigned char is_ys[2][kScoreVariantBlockSize];
  unsigned char ploidy_m1s[2][kScoreVariantBlockSize];
  uint32_t cur_variant_batch_size;

  // per-thread buffers
  uint64_t** ddosage_incrs;
  double** dosages_vmaj;

  // final results
  uint64_t* ddosage_sums;
  double** sharded_final_scores_cmaj;
} CalcScoreCtx;

THREAD_FUNC_DECL CalcScoreThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  CalcScoreCtx* ctx = S_CAST(CalcScoreCtx*, arg->sharedp->context);

  const uintptr_t* variant_include = ctx->variant_include;
  const uint32_t* variant_include_cumulative_popcounts = ctx->variant_include_cumulative_popcounts;
  const uintptr_t* qsr_include = ctx->qsr_include;

  const uint32_t score_final_col_ct = ctx->score_final_col_ct;
  uint32_t sample_shard_size = ctx->sample_shard_size;
  const uint32_t sample_idx_start = tidx * sample_shard_size;
  double* shard_final_scores_cmaj = ctx->sharded_final_scores_cmaj[tidx];
  const uint32_t sample_ct = ctx->sample_ct;
  if (sample_shard_size + sample_idx_start > sample_ct) {
    sample_shard_size = sample_ct - sample_idx_start;
  }
  const uint32_t sample_idx_startl2 = NypCtToWordCt(sample_idx_start);
  const uint32_t sample_idx_startl = BitCtToWordCt(sample_idx_start);
  const uint32_t shard_sizel = BitCtToWordCt(sample_shard_size);
  const uint32_t sample_ctaw2 = NypCtToAlignedWordCt(sample_ct);
  const uint32_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
  const uint32_t max_difflist_len = ctx->max_difflist_len;
  const uintptr_t raregeno_stride = NypCtToAlignedWordCt(max_difflist_len);
  const uintptr_t difflist_sample_ids_stride = RoundUpPow2(max_difflist_len, kInt32PerVec);
  const uint32_t dosage_main_stride = RoundUpPow2(sample_ct, kDosagePerVec);
  const uintptr_t* shard_sex_nonmale_collapsed = &(ctx->sex_nonmale_collapsed[sample_idx_startl]);
  const uintptr_t* shard_sex_female_collapsed = &(ctx->sex_female_collapsed[sample_idx_startl]);
  const uint32_t shard_nonmale_ct = PopcountWords(shard_sex_nonmale_collapsed, shard_sizel);
  const uint32_t shard_male_ct = sample_shard_size - shard_nonmale_ct;
  const uint32_t shard_female_ct = PopcountWords(shard_sex_female_collapsed, shard_sizel);
  const uint32_t shard_nonfemale_ct = sample_shard_size - shard_female_ct;
  const ScoreFlags flags = ctx->flags;
  const uint32_t model_dominant = (flags / kfScoreDominant) & 1;
  const uint32_t domrec = model_dominant || (flags & kfScoreRecessive);
  const uint32_t no_meanimpute = (flags / kfScoreNoMeanimpute) & 1;
  const uint32_t se_mode = (flags / kfScoreSe) & 1;
  const uint32_t qsr_ct = ctx->qsr_ct;

  uint64_t* ddosage_incrs = ctx->ddosage_incrs[tidx];
  double* dosages_vmaj = ctx->dosages_vmaj[tidx];

  uint64_t* shard_ddosage_sums = &(ctx->ddosage_sums[sample_idx_start]);

  double cur_allele_freq = 0.0;
  double geno_slope = kRecipDosageMax;
  double geno_intercept = 0.0;
  uint32_t parity = 0;
  do {
    const uint32_t cur_variant_batch_size = ctx->cur_variant_batch_size;
    if (cur_variant_batch_size) {
      const uint32_t* variant_uidxs = ctx->variant_uidxs[parity];
      const uintptr_t* genovec_iter = &(ctx->genovecs[parity][sample_idx_startl2]);
      const uintptr_t* raregeno_iter = ctx->raregenos[parity];
      const uint32_t* difflist_lens = ctx->difflist_lens[parity];
      const uint32_t* difflist_sample_ids_iter = ctx->difflist_sample_ids[parity];
      const int8_t* difflist_common_genos = ctx->difflist_common_genos[parity];
      const uintptr_t* dosage_present_iter = ctx->dosage_presents[parity];
      const Dosage* dosage_main_iter = ctx->dosage_mains[parity];
      const uint32_t* dosage_cts = ctx->dosage_cts[parity];
      const uintptr_t* missing_bitvec_iter = &(ctx->missing_bitvecs[parity][sample_idx_startl]);
      const uintptr_t* missing_male_bitvec_iter = &(ctx->missing_male_bitvecs[parity][sample_idx_startl]);
      const uintptr_t* missing_nonfemale_bitvec_iter = &(ctx->missing_nonfemale_bitvecs[parity][sample_idx_startl]);
      const double* allele_freqs = ctx->allele_freqs[parity];
      const double* geno_slopes = ctx->geno_slopes[parity];
      const double* geno_intercepts = ctx->geno_intercepts[parity];
      const double* score_sparse_coefs_vmaj_iter = ctx->score_sparse_coefs_vmaj[parity];
      const unsigned char* is_nonx_haploids = &(ctx->is_nonx_haploids[parity][0]);
      const unsigned char* is_relevant_xs = &(ctx->is_relevant_xs[parity][0]);
      const unsigned char* is_ys = &(ctx->is_ys[parity][0]);
      const unsigned char* ploidy_m1s = &(ctx->ploidy_m1s[parity][0]);
      const uintptr_t* shard_dosage_present = nullptr;
      const Dosage* shard_dosage_main = nullptr;
      double* dosages_vmaj_iter = dosages_vmaj;

      uint32_t sparse_vidx = 0;
      uint32_t shard_dosage_ct = 0;
      for (uint32_t vidx = 0; vidx != cur_variant_batch_size; ++vidx) {
        const uint32_t is_nonx_haploid = is_nonx_haploids[vidx];
        const uint32_t is_relevant_x = is_relevant_xs[vidx];
        const uint32_t is_y = is_ys[vidx];
        if (dosage_cts) {
          // defensive
          shard_dosage_present = nullptr;
          shard_dosage_main = nullptr;
          // bugfix (3 Jun 2022)
          shard_dosage_ct = 0;
          const uint32_t dosage_ct = dosage_cts[vidx];
          if (dosage_ct) {
            const uint32_t dosage_main_offset = PopcountWords(dosage_present_iter, sample_idx_startl);
            if (dosage_main_offset != dosage_ct) {
              shard_dosage_ct = PopcountWords(&(dosage_present_iter[sample_idx_startl]), shard_sizel);
              if (shard_dosage_ct) {
                shard_dosage_present = &(dosage_present_iter[sample_idx_startl]);
                shard_dosage_main = &(dosage_main_iter[dosage_main_offset]);
              }
            }
          }
        } else if (!domrec) {
          // For now, the sparse-genotype optimization is only considered if
          // there are no dosages in the entire .pgen, and
          // 'dominant'/'recessive' wasn't specified.
          // It's actually straightforward to support domrec
          // (modifying difflist_common_geno and raregeno upfront is almost
          // enough), but let's get this working first.
          const uint32_t difflist_common_geno = S_CAST(int32_t, difflist_common_genos[vidx]);
          if (difflist_common_geno != UINT32_MAX) {
            const uint32_t difflist_len = difflist_lens[sparse_vidx];
            if (difflist_len) {
              uint32_t difflist_start_idx = 0;
              if (sample_idx_start) {
                difflist_start_idx = LowerBoundNonemptyU32(difflist_sample_ids_iter, difflist_len, sample_idx_start);
              }
              if (difflist_len > difflist_start_idx) {
                const uint32_t shard_difflist_len = LowerBoundNonemptyU32(&(difflist_sample_ids_iter[difflist_start_idx]), difflist_len - difflist_start_idx, sample_idx_start + sample_shard_size);
                if (shard_difflist_len) {
                  // Shard contains at least one uncommon genotype, so we have
                  // some work to do.
                  if (allele_freqs) {
                    cur_allele_freq = allele_freqs[vidx];
                    geno_slope = geno_slopes[vidx];
                    geno_intercept = geno_intercepts[vidx];
                  }
                  uint64_t ddosage_incr_lookup_table[4];
                  ddosage_incr_lookup_table[0] = 0;
                  ddosage_incr_lookup_table[1] = kDosageMax;
                  ddosage_incr_lookup_table[2] = kDosageMax * 2LL;
                  ddosage_incr_lookup_table[3] = 0;
                  const uint64_t common_ddosage = ddosage_incr_lookup_table[difflist_common_geno];
                  if (common_ddosage) {
                    for (uint32_t uii = 0; uii != 4; ++uii) {
                      // deliberate underflow
                      ddosage_incr_lookup_table[uii] -= common_ddosage;
                    }
                  }
                  double dosage_incr_lookup_table[4];
                  dosage_incr_lookup_table[0] = geno_intercept;
                  const double geno_incr = kDosageMax * geno_slope;
                  dosage_incr_lookup_table[1] = geno_incr + geno_intercept;
                  dosage_incr_lookup_table[2] = dosage_incr_lookup_table[1] + geno_incr;
                  dosage_incr_lookup_table[3] = 0.0;
                  if (!no_meanimpute) {
                    dosage_incr_lookup_table[3] = kDosageMax * 2LL * cur_allele_freq * geno_slope;
                  }
                  if (se_mode) {
                    for (uint32_t uii = 0; uii != 4; ++uii) {
                      dosage_incr_lookup_table[uii] *= dosage_incr_lookup_table[uii];
                    }
                  }
                  const double common_dosage = dosage_incr_lookup_table[difflist_common_geno];
                  for (uint32_t uii = 0; uii != 4; ++uii) {
                    dosage_incr_lookup_table[uii] -= common_dosage;
                  }

                  const uint32_t difflist_end_idx = difflist_start_idx + shard_difflist_len;
                  const uint32_t widx_last = (difflist_end_idx - 1) / kBitsPerWordD2;
                  uint32_t widx = difflist_start_idx / kBitsPerWordD2;
                  const uint32_t* cur_difflist_sample_ids = &(difflist_sample_ids_iter[widx * kBitsPerWordD2]);
                  uint32_t difflist_idx_lowbits = difflist_start_idx % kBitsPerWordD2;
                  uint32_t difflist_idx_lowbits_end = kBitsPerWordD2;
                  if (widx == widx_last) {
                    difflist_idx_lowbits_end = ModNz(difflist_end_idx, kBitsPerWordD2);
                  }
                  uintptr_t raregeno_word = raregeno_iter[widx] >> (2 * difflist_idx_lowbits);
                  while (1) {
                    for (; difflist_idx_lowbits != difflist_idx_lowbits_end; ++difflist_idx_lowbits) {
                      const uintptr_t cur_geno = raregeno_word & 3;
                      const uint32_t shard_sample_idx = cur_difflist_sample_ids[difflist_idx_lowbits] - sample_idx_start;
                      uint64_t* ddosage_sum_col_base = &(shard_ddosage_sums[shard_sample_idx]);
                      const uint64_t cur_ddosage_incr = ddosage_incr_lookup_table[cur_geno];
                      if (!qsr_ct) {
                        ddosage_sum_col_base[0] += cur_ddosage_incr;
                      } else {
                        const uint32_t variant_uidx = variant_uidxs[vidx];
                        const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
                        for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                          if (IsSet(qsr_include, qsr_idx + bit_idx_base)) {
                            ddosage_sum_col_base[qsr_idx * sample_ct] += cur_ddosage_incr;
                          }
                        }
                      }
                      double* final_score_col_iter = &(shard_final_scores_cmaj[shard_sample_idx]);
                      const double cur_dosage_incr = dosage_incr_lookup_table[cur_geno];
                      for (uint32_t uii = 0; uii != score_final_col_ct; ++uii) {
                        *final_score_col_iter = prefer_fma(cur_dosage_incr, score_sparse_coefs_vmaj_iter[uii], *final_score_col_iter);
                        final_score_col_iter = &(final_score_col_iter[sample_shard_size]);
                      }
                      raregeno_word = raregeno_word >> 2;
                    }
                    ++widx;
                    if (widx >= widx_last) {
                      if (widx > widx_last) {
                        break;
                      }
                      difflist_idx_lowbits_end = ModNz(difflist_end_idx, kBitsPerWordD2);
                    }
                    cur_difflist_sample_ids = &(cur_difflist_sample_ids[kBitsPerWordD2]);
                    difflist_idx_lowbits = 0;
                    raregeno_word = raregeno_iter[widx];
                  }
                }
              }
            }

            genovec_iter = &(genovec_iter[sample_ctaw2]);
            missing_bitvec_iter = &(missing_bitvec_iter[sample_ctaw]);
            missing_male_bitvec_iter = &(missing_male_bitvec_iter[sample_ctaw]);
            missing_nonfemale_bitvec_iter = &(missing_nonfemale_bitvec_iter[sample_ctaw]);
            score_sparse_coefs_vmaj_iter = &(score_sparse_coefs_vmaj_iter[score_final_col_ct]);
            ++sparse_vidx;
            raregeno_iter = &(raregeno_iter[raregeno_stride]);
            difflist_sample_ids_iter = &(difflist_sample_ids_iter[difflist_sample_ids_stride]);
            continue;
          }
        }
        FillDdosageInts(genovec_iter, shard_dosage_present, shard_dosage_main, sample_shard_size, shard_dosage_ct, 1 - is_nonx_haploid, ddosage_incrs);
        if (is_y) {
          if (shard_female_ct) {
            uintptr_t shard_sample_idx_base = 0;
            uintptr_t shard_sex_female_collapsed_bits = shard_sex_female_collapsed[0];
            for (uint32_t shard_female_idx = 0; shard_female_idx != shard_female_ct; ++shard_female_idx) {
              const uintptr_t shard_sample_idx = BitIter1(shard_sex_female_collapsed, &shard_sample_idx_base, &shard_sex_female_collapsed_bits);
              ddosage_incrs[shard_sample_idx] = 0;
            }
          }
        } else if (is_relevant_x) {
          uintptr_t shard_sample_idx_base = 0;
          uintptr_t shard_sex_nonmale_collapsed_inv_bits = ~shard_sex_nonmale_collapsed[0];
          for (uint32_t shard_male_idx = 0; shard_male_idx != shard_male_ct; ++shard_male_idx) {
            const uintptr_t shard_sample_idx = BitIter0(shard_sex_nonmale_collapsed, &shard_sample_idx_base, &shard_sex_nonmale_collapsed_inv_bits);
            ddosage_incrs[shard_sample_idx] /= 2;
          }
        }
        if (domrec) {
          if (model_dominant) {
            for (uint32_t shard_sample_idx = 0; shard_sample_idx != sample_shard_size; ++shard_sample_idx) {
              if (ddosage_incrs[shard_sample_idx] > kDosageMax) {
                ddosage_incrs[shard_sample_idx] = kDosageMax;
              }
            }
          } else {
            // recessive
            for (uint32_t shard_sample_idx = 0; shard_sample_idx != sample_shard_size; ++shard_sample_idx) {
              uint64_t cur_ddosage_incr = ddosage_incrs[shard_sample_idx];
              if (cur_ddosage_incr <= kDosageMax) {
                cur_ddosage_incr = 0;
              } else {
                cur_ddosage_incr -= kDosageMax;
              }
              ddosage_incrs[shard_sample_idx] = cur_ddosage_incr;
            }
          }
        }
        if (!qsr_ct) {
          for (uint32_t shard_sample_idx = 0; shard_sample_idx != sample_shard_size; ++shard_sample_idx) {
            shard_ddosage_sums[shard_sample_idx] += ddosage_incrs[shard_sample_idx];
          }
        } else {
          const uint32_t variant_uidx = variant_uidxs[vidx];
          const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
          for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
            if (IsSet(qsr_include, qsr_idx + bit_idx_base)) {
              uint64_t* cur_ddosage_sums = &(shard_ddosage_sums[qsr_idx * sample_ct]);
              for (uint32_t shard_sample_idx = 0; shard_sample_idx != sample_shard_size; ++shard_sample_idx) {
                cur_ddosage_sums[shard_sample_idx] += ddosage_incrs[shard_sample_idx];
              }
            }
          }
        }
        if (allele_freqs) {
          cur_allele_freq = allele_freqs[vidx];
          geno_slope = geno_slopes[vidx];
          geno_intercept = geno_intercepts[vidx];
        }
        if ((!shard_dosage_ct) && (!is_nonx_haploid) && (!is_relevant_x)) {
          // Fast path for common no-dosage case.
          // (We wait till this point to branch since the shard_ddosage_sums
          // update isn't that expensive.)
          double lookup_table[32] ALIGNV16;
          lookup_table[0] = geno_intercept;
          const double geno_incr = kDosageMax * geno_slope;
          const double geno_one = geno_incr + geno_intercept;
          if (!domrec) {
            lookup_table[2] = geno_one;
            lookup_table[4] = geno_incr + geno_one;
          } else {
            lookup_table[2] = model_dominant? geno_one : geno_intercept;
            lookup_table[4] = geno_one;
          }
          lookup_table[6] = 0.0;
          if (!no_meanimpute) {
            double missing_effect = kDosageMax * cur_allele_freq * geno_slope;
            if (!domrec) {
              missing_effect *= 2;
            }
            lookup_table[6] = missing_effect;
          }
          if (se_mode) {
            lookup_table[0] *= lookup_table[0];
            lookup_table[2] *= lookup_table[2];
            lookup_table[4] *= lookup_table[4];
            lookup_table[6] *= lookup_table[6];
          }
          InitLookup16x8bx2(lookup_table);
          GenoarrLookup16x8bx2(genovec_iter, lookup_table, sample_shard_size, dosages_vmaj_iter);
        } else {
          const uint32_t shard_missing_ct = PopcountWords(missing_bitvec_iter, shard_sizel);
          if (shard_missing_ct) {
            double missing_effect = 0.0;
            if (!no_meanimpute) {
              missing_effect = kDosageMax * cur_allele_freq * geno_slope;
            }
            if (is_y) {
              ZeroDArr(sample_shard_size, dosages_vmaj_iter);
              if ((!no_meanimpute) && shard_nonfemale_ct) {
                for (uint32_t shard_widx = 0; shard_widx != shard_sizel; ++shard_widx) {
                  uintptr_t cur_missing_nonfemale_bits = missing_nonfemale_bitvec_iter[shard_widx];
                  if (!cur_missing_nonfemale_bits) {
                    continue;
                  }
                  double* cur_dosages_vmaj_iter = &(dosages_vmaj_iter[shard_widx * kBitsPerWord]);
                  do {
                    const uint32_t sample_idx_lowbits = ctzw(cur_missing_nonfemale_bits);
                    cur_dosages_vmaj_iter[sample_idx_lowbits] = missing_effect;
                    cur_missing_nonfemale_bits &= cur_missing_nonfemale_bits - 1;
                  } while (cur_missing_nonfemale_bits);
                }
              }
            } else if (is_relevant_x) {
              ZeroDArr(sample_shard_size, dosages_vmaj_iter);
              if (!no_meanimpute) {
                if (shard_male_ct) {
                  for (uint32_t shard_widx = 0; shard_widx != shard_sizel; ++shard_widx) {
                    uintptr_t cur_missing_male_bits = missing_male_bitvec_iter[shard_widx];
                    if (!cur_missing_male_bits) {
                      continue;
                    }
                    double* cur_dosages_vmaj_iter = &(dosages_vmaj_iter[shard_widx * kBitsPerWord]);
                    do {
                      const uint32_t sample_idx_lowbits = ctzw(cur_missing_male_bits);
                      cur_dosages_vmaj_iter[sample_idx_lowbits] = missing_effect;
                      cur_missing_male_bits &= cur_missing_male_bits - 1;
                    } while (cur_missing_male_bits);
                  }
                }
                if (shard_nonmale_ct) {
                  missing_effect *= 2;
                  for (uint32_t shard_widx = 0; shard_widx != shard_sizel; ++shard_widx) {
                    uintptr_t cur_missing_nonmale_bits = missing_bitvec_iter[shard_widx] & shard_sex_nonmale_collapsed[shard_widx];
                    if (!cur_missing_nonmale_bits) {
                      continue;
                    }
                    double* cur_dosages_vmaj_iter = &(dosages_vmaj_iter[shard_widx * kBitsPerWord]);
                    do {
                      const uint32_t sample_idx_lowbits = ctzw(cur_missing_nonmale_bits);
                      cur_dosages_vmaj_iter[sample_idx_lowbits] = missing_effect;
                      cur_missing_nonmale_bits &= cur_missing_nonmale_bits - 1;
                    } while (cur_missing_nonmale_bits);
                  }
                }
              }
            } else {
              const double ploidy_d = ploidy_m1s[vidx]? 2.0 : 1.0;
              missing_effect *= ploidy_d;
              uintptr_t shard_sample_idx_base = 0;
              uintptr_t missing_bits = missing_bitvec_iter[0];
              for (uint32_t missing_idx = 0; missing_idx != shard_missing_ct; ++missing_idx) {
                const uintptr_t shard_sample_idx = BitIter1(missing_bitvec_iter, &shard_sample_idx_base, &missing_bits);
                dosages_vmaj_iter[shard_sample_idx] = missing_effect;
              }
            }
          }
          const uint32_t shard_nm_sample_ct = sample_shard_size - shard_missing_ct;
          uintptr_t shard_sample_idx_base = 0;
          uintptr_t missing_inv_bits = ~missing_bitvec_iter[0];
          for (uint32_t shard_nm_sample_idx = 0; shard_nm_sample_idx != shard_nm_sample_ct; ++shard_nm_sample_idx) {
            const uintptr_t shard_sample_idx = BitIter0(missing_bitvec_iter, &shard_sample_idx_base, &missing_inv_bits);
            dosages_vmaj_iter[shard_sample_idx] = u63tod(ddosage_incrs[shard_sample_idx]) * geno_slope + geno_intercept;
          }
          if (se_mode) {
            // Suppose our score coefficients are drawn from independent
            // Gaussians.  Then the variance of the final score average is the
            // sum of the variances of the individual terms, divided by (T^2)
            // where T is the number of terms.  These individual variances are
            // of the form (<genotype value> * <stdev>)^2.
            //
            // Thus, we can use the same inner loop to compute standard errors,
            // as long as
            //   1. we square the genotypes and the standard errors before
            //      matrix multiplication, and
            //   2. we take the square root of the sums at the end.
            for (uint32_t shard_sample_idx = 0; shard_sample_idx != sample_shard_size; ++shard_sample_idx) {
              dosages_vmaj_iter[shard_sample_idx] *= dosages_vmaj_iter[shard_sample_idx];
            }
          }
        }
        dosages_vmaj_iter = &(dosages_vmaj_iter[sample_shard_size]);

        genovec_iter = &(genovec_iter[sample_ctaw2]);
        missing_bitvec_iter = &(missing_bitvec_iter[sample_ctaw]);
        missing_male_bitvec_iter = &(missing_male_bitvec_iter[sample_ctaw]);
        missing_nonfemale_bitvec_iter = &(missing_nonfemale_bitvec_iter[sample_ctaw]);
        if (dosage_present_iter) {
          dosage_present_iter = &(dosage_present_iter[sample_ctaw]);
          dosage_main_iter = &(dosage_main_iter[dosage_main_stride]);
        }
      }

      const uint32_t dense_variant_ct = cur_variant_batch_size - sparse_vidx;
      if (dense_variant_ct) {
        const double* score_dense_coefs_cmaj = ctx->score_dense_coefs_cmaj[parity];
        RowMajorMatrixMultiplyStridedIncr(score_dense_coefs_cmaj, dosages_vmaj, score_final_col_ct, kScoreVariantBlockSize, sample_shard_size, sample_shard_size, dense_variant_ct, sample_shard_size, shard_final_scores_cmaj);
      }
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

typedef struct ParsedQscoreRangeStruct {
  char* range_name;
  double lbound;
  double ubound;
} ParsedQscoreRange;

PglErr ScoreReport(const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uintptr_t* variant_include, const ChrInfo* cip, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const double* allele_freqs, const ScoreInfo* score_info_ptr, const char* output_missing_pheno, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t nosex_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_variant_id_slen, uint32_t xchr_model, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uintptr_t line_idx = 0;
  const char* cur_input_fname = nullptr;
  char* cswritep = nullptr;
  FILE* score_tmpfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  TextStream score_txs;
  ThreadGroup tg;
  CompressStreamState css;
  PreinitTextStream(&score_txs);
  PreinitThreads(&tg);
  PreinitCstream(&css);
  {
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    if (!xchr_model) {
      uint32_t x_code;
      if (XymtExists(cip, kChrOffsetX, &x_code)) {
        uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[x_code];
        uint32_t x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
        uint32_t x_end = cip->chr_fo_vidx_start[x_chr_fo_idx + 1];
        if (!AllBitsAreZero(variant_include, x_start, x_end)) {
          uintptr_t* variant_include_no_x;
          if (unlikely(bigstack_alloc_w(raw_variant_ctl, &variant_include_no_x))) {
            goto ScoreReport_ret_NOMEM;
          }
          memcpy(variant_include_no_x, variant_include, raw_variant_ctl * sizeof(intptr_t));
          ClearBitsNz(x_start, x_end, variant_include_no_x);
          variant_include = variant_include_no_x;
        }
      }
    } else if (xchr_model == 2) {
      xchr_model = 0;
    }
    // now xchr_model is set iff it's 1

    const ScoreFlags flags = score_info_ptr->flags;
    const uint32_t multi_input = (flags / kfScoreMultiInput) & 1;
    if (multi_input) {
      assert(!(flags & (kfScoreListVariants | kfScoreColNallele | kfScoreColDenom | kfScoreColDosageSum)));
    }
    const uint32_t output_zst = (flags / kfScoreZs) & 1;
    uint32_t* variant_id_htable = nullptr;
    uint32_t variant_id_htable_size = 0;
    uint32_t* variant_include_cumulative_popcounts = nullptr;
    uintptr_t* qsr_include = nullptr;
    char** range_names = nullptr;
    uintptr_t qsr_ct = 0;
    if (multi_input || score_info_ptr->qsr_range_fname || score_info_ptr->partition_fname) {
      // Limit this to ~1/8 of available memory, since memory may be tight.
      variant_id_htable_size = GetHtableFastSize(variant_ct);
      const uintptr_t htable_size_limit = bigstack_left() / (8 * sizeof(int32_t));
      if (variant_id_htable_size > htable_size_limit) {
        variant_id_htable_size = htable_size_limit;
        const uint32_t htable_size_min = GetHtableMinSize(variant_ct);
        if (htable_size_min > variant_id_htable_size) {
          variant_id_htable_size = htable_size_min;
        }
      }
      if (unlikely(bigstack_alloc_u32(variant_id_htable_size, &variant_id_htable))) {
        goto ScoreReport_ret_NOMEM;
      }
      reterr = PopulateIdHtableMt(nullptr, variant_include, variant_ids, variant_ct, 0, variant_id_htable_size, max_thread_ct, nullptr, variant_id_htable, nullptr);
      if (unlikely(reterr)) {
        goto ScoreReport_ret_1;
      }
    }
    if (score_info_ptr->qsr_range_fname) {
      // Strictly speaking, textFILE would be more appropriate for the range
      // file since it should be tiny, but it doesn't really matter.
      // We still reserve bigstack_left() / 8 for the line-buffer since the
      // data file usually contains allele codes, and we use TextRetarget()
      // below to use the buffer allocated here for the data file too (and for
      // the score file later).
      reterr = SizeAndInitTextStream(score_info_ptr->qsr_range_fname, bigstack_left() / 8, 1, &score_txs);
      if (unlikely(reterr)) {
        goto ScoreReport_ret_QSR_RANGE_TSTREAM_FAIL;
      }
      unsigned char* bigstack_mark2 = g_bigstack_base;
      // strlen("<prefix>.<range name>.sscore[.zst]") < kPglFnamesize
      const uint32_t max_name_slen = kPglFnamesize - S_CAST(uintptr_t, outname_end - outname) - 9 - output_zst * 4;
      ParsedQscoreRange* parsed_qscore_ranges = R_CAST(ParsedQscoreRange*, g_bigstack_base);
      unsigned char* tmp_alloc_end = g_bigstack_end;
      uintptr_t miss_ct = 0;
      while (1) {
        ++line_idx;
        const char* line_start = TextGet(&score_txs);
        if (!line_start) {
          if (likely(!TextStreamErrcode2(&score_txs, &reterr))) {
            break;
          }
          goto ScoreReport_ret_QSR_RANGE_TSTREAM_FAIL;
        }
        // range name, p-value lower bound, p-value upper bound
        const char* range_name_end = CurTokenEnd(line_start);
        const char* lbound_start = FirstNonTspace(range_name_end);
        double lbound;
        // PLINK 1.9 documentation promises that lines with too few entries or
        // nonnumeric values in the second and third column are ignored.
        const char* lbound_end = ScantokDouble(lbound_start, &lbound);
        if (!lbound_end) {
          continue;
        }
        const char* ubound_start = FirstNonTspace(lbound_end);
        double ubound;
        const char* ubound_end = ScantokDouble(ubound_start, &ubound);
        if (!ubound_end) {
          continue;
        }
        if (unlikely(lbound > ubound)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Upper bound < lower bound on line %" PRIuPTR " of --q-score-range range file.\n", line_idx);
          goto ScoreReport_ret_MALFORMED_INPUT_WW;
        }
        const uint32_t name_slen = range_name_end - line_start;
        if (name_slen > max_name_slen) {
          snprintf(g_logbuf, kLogbufSize, "Error: Name too long on line %" PRIuPTR " of --q-score-range range file.\n", line_idx);
        }
        unsigned char* tmp_alloc_base = R_CAST(unsigned char*, &(parsed_qscore_ranges[qsr_ct]));
        if (S_CAST(uintptr_t, tmp_alloc_end - tmp_alloc_base) <= name_slen + sizeof(ParsedQscoreRange)) {
          goto ScoreReport_ret_NOMEM;
        }
        tmp_alloc_end -= name_slen + 1;
        char* stored_name = R_CAST(char*, tmp_alloc_end);
        memcpyx(stored_name, line_start, name_slen, '\0');
        parsed_qscore_ranges[qsr_ct].range_name = stored_name;
        parsed_qscore_ranges[qsr_ct].lbound = lbound;
        parsed_qscore_ranges[qsr_ct].ubound = ubound;
        ++qsr_ct;
      }
      if (unlikely(!qsr_ct)) {
        logerrputs("Error: Empty --q-score-range range file.\n");
        goto ScoreReport_ret_INCONSISTENT_INPUT;
      }
      BigstackEndSet(tmp_alloc_end);
      if (unlikely(BigstackBaseSetChecked(&(parsed_qscore_ranges[qsr_ct])))) {
        goto ScoreReport_ret_NOMEM;
      }
#ifndef LAPACK_ILP64
      if (unlikely(qsr_ct > (0x7fffffff / kScoreVariantBlockSize))) {
        logerrputs("Error: --q-score-range range count too large for this " PROG_NAME_STR " build.  If this is\nreally the computation you want, use a " PROG_NAME_STR " build with large-matrix support.\n");
        goto ScoreReport_ret_INCONSISTENT_INPUT;
      }
#  ifndef __LP64__
      const uint64_t bit_ct = S_CAST(uint64_t, qsr_ct) * variant_ct;
      if (unlikely(bit_ct > 0xffffffffU)) {
        goto ScoreReport_ret_NOMEM;
      }
#  endif
#endif
      if (unlikely(bigstack_end_alloc_u32(raw_variant_ctl, &variant_include_cumulative_popcounts) ||
                   bigstack_end_calloc_w(BitCtToWordCt(S_CAST(uint64_t, qsr_ct) * variant_ct), &qsr_include) ||
                   bigstack_end_alloc_cp(qsr_ct, &range_names))) {
        goto ScoreReport_ret_NOMEM;
      }
      for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
        range_names[qsr_idx] = parsed_qscore_ranges[qsr_idx].range_name;
      }
      const uint32_t variant_ctl = BitCtToWordCt(variant_ct);
      uintptr_t* already_seen;
      if (unlikely(bigstack_calloc_w(variant_ctl, &already_seen))) {
        goto ScoreReport_ret_NOMEM;
      }
      FillCumulativePopcounts(variant_include, raw_variant_ctl, variant_include_cumulative_popcounts);
      reterr = TextRetarget(score_info_ptr->qsr_data_fname, &score_txs);
      if (unlikely(reterr)) {
        goto ScoreReport_ret_QSR_RANGE_TSTREAM_FAIL;
      }
      const uint32_t colid_first = (score_info_ptr->qsr_varid_col_p1 < score_info_ptr->qsr_val_col_p1);
      uint32_t colmin;
      uint32_t coldiff;
      if (colid_first) {
        colmin = score_info_ptr->qsr_varid_col_p1 - 1;
        coldiff = score_info_ptr->qsr_val_col_p1 - score_info_ptr->qsr_varid_col_p1;
      } else {
        colmin = score_info_ptr->qsr_val_col_p1 - 1;
        coldiff = score_info_ptr->qsr_varid_col_p1 - score_info_ptr->qsr_val_col_p1;
      }
      line_idx = 0;
      miss_ct = 0;
      if (flags & kfScoreQsrHeader) {
        ++line_idx;
        if (unlikely(!TextGet(&score_txs))) {
          if (!TextStreamErrcode2(&score_txs, &reterr)) {
            logerrputs("Error: Empty --q-score-range data file.\n");
            goto ScoreReport_ret_MALFORMED_INPUT;
          }
          goto ScoreReport_ret_QSR_DATA_TSTREAM_FAIL;
        }
      }
      double* min_vals = nullptr;
      if (flags & kfScoreQsrMin) {
        // something like this is needed to handle --glm output for
        // multiallelic variants.
        // (possible todo: --glm modifier which requests all-allele joint tests
        // for multiallelic variants)
        if (unlikely(bigstack_alloc_d(variant_ct, &min_vals))) {
          goto ScoreReport_ret_NOMEM;
        }
      }
      while (1) {
        ++line_idx;
        const char* line_start = TextGet(&score_txs);
        if (!line_start) {
          if (likely(!TextStreamErrcode2(&score_txs, &reterr))) {
            break;
          }
          goto ScoreReport_ret_QSR_DATA_TSTREAM_FAIL;
        }
        const char* colid_ptr;
        const char* colval_ptr;
        if (colid_first) {
          colid_ptr = NextTokenMult0(line_start, colmin);
          colval_ptr = NextTokenMult(colid_ptr, coldiff);
          if (unlikely(!colval_ptr)) {
            goto ScoreReport_ret_QSR_DATA_MISSING_TOKENS;
          }
        } else {
          colval_ptr = NextTokenMult0(line_start, colmin);
          colid_ptr = NextTokenMult(colval_ptr, coldiff);
          if (unlikely(!colid_ptr)) {
            goto ScoreReport_ret_QSR_DATA_MISSING_TOKENS;
          }
        }
        const uint32_t varid_slen = strlen_se(colid_ptr);
        const uint32_t variant_uidx = VariantIdDupflagHtableFind(colid_ptr, variant_ids, variant_id_htable, varid_slen, variant_id_htable_size, max_variant_id_slen);
        if ((variant_uidx >> 31) || (!IsSet(variant_include, variant_uidx))) {
          ++miss_ct;
          continue;
        }
        double cur_val;
        if (!ScantokDouble(colval_ptr, &cur_val)) {
          // Tolerate NA without erroring out.  (Could count this as seen, but
          // that would force the min_vals logic to be more complicated.)
          const char* colval_end = CurTokenEnd(colval_ptr);
          if (likely(IsNanStr(colval_ptr, colval_end - colval_ptr))) {
            continue;
          }
          *K_CAST(char*, colval_end) = '\0';
          logerrprintfww("Error: Invalid value '%s' on line %" PRIuPTR " of --q-score-range data file.\n", colval_ptr, line_idx);
          goto ScoreReport_ret_MALFORMED_INPUT;
        }
        const uint32_t variant_idx = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx);
        const uintptr_t bit_idx_base = variant_idx * qsr_ct;
        if (min_vals) {
          if (IsSet(already_seen, variant_idx)) {
            if (min_vals[variant_idx] <= cur_val) {
              continue;
            }
            ClearBitsNz(bit_idx_base, bit_idx_base + qsr_ct, qsr_include);
          }
          min_vals[variant_idx] = cur_val;
        } else {
          if (IsSet(already_seen, variant_idx)) {
            logerrprintfww("Error: Duplicate variant ID '%s' in --q-score-range data file. (Add the 'min' modifier if this is a multiallelic variant that you want to use the minimum p-value for.)\n", variant_ids[variant_uidx]);
            goto ScoreReport_ret_MALFORMED_INPUT;
          }
        }
        SetBit(variant_idx, already_seen);
        for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
          if ((cur_val < parsed_qscore_ranges[qsr_idx].lbound) || (cur_val > parsed_qscore_ranges[qsr_idx].ubound)) {
            continue;
          }
          SetBit(bit_idx_base + qsr_idx, qsr_include);
        }
      }
      const uint32_t qsr_variant_ct = PopcountWords(already_seen, variant_ctl);
      if (unlikely(!qsr_variant_ct)) {
        logerrputs("Error: No valid entries in --q-score-range data file.\n");
        goto ScoreReport_ret_INCONSISTENT_INPUT;
      }
      logprintf("--q-score-range: %" PRIuPTR " range%s and %u variant%s loaded.\n", qsr_ct, (qsr_ct == 1)? "" : "s", qsr_variant_ct, (qsr_variant_ct == 1)? "" : "s");
      if (miss_ct) {
        logerrprintf("Warning: %" PRIuPTR " line%s skipped in --q-score-range data file.\n", miss_ct, (miss_ct == 1)? "" : "s");
      }
      // possible todo: replace variant_include with already_seen, and compact
      // qsr_include.
      // but for now, we just free already_seen, and in the common use cases
      // this should be fine.
      reterr = TextRetarget(score_info_ptr->input_fname, &score_txs);
      if (unlikely(reterr)) {
        cur_input_fname = score_info_ptr->input_fname;
        goto ScoreReport_ret_TSTREAM_FAIL;
      }
      BigstackReset(bigstack_mark2);
      line_idx = 0;
    } else if (score_info_ptr->partition_fname) {
      // --score-partitions: fill the SAME qsr_include membership bitarray from
      // an explicit 0/1 mask.  Everything downstream -- the per-partition
      // denominators in allele_ct_bases[], the [sample x partition x model]
      // accumulation, the missingness accounting -- is the --q-score-range
      // machinery unchanged.  Only membership and the output shape differ.
      //
      // Partition 0 is the synthetic TOTAL: every variant set, so it reproduces
      // the unpartitioned score through the same code path rather than a
      // special case.  range_names[0] is null, which is how the header writer
      // knows to emit "<model>_AVG" with no partition suffix.
      reterr = SizeAndInitTextStream(score_info_ptr->partition_fname, bigstack_left() / 8, 1, &score_txs);
      if (unlikely(reterr)) {
        logerrprintfww("Error: --score-partitions: failed to open %s .\n", score_info_ptr->partition_fname);
        goto ScoreReport_ret_1;
      }
      unsigned char* bigstack_mark2 = g_bigstack_base;
      {
        char* header_line = TextGet(&score_txs);
        if (unlikely(!header_line)) {
          if (!TextStreamErrcode2(&score_txs, &reterr)) {
            logerrputs("Error: --score-partitions: file is empty.\n");
            reterr = kPglRetMalformedInput;
          } else {
            TextStreamErrPrint(score_info_ptr->partition_fname, &score_txs);
          }
          goto ScoreReport_ret_1;
        }
        // Copy the header before parsing: TextGet's buffer is reused by the
        // next line, and the partition names have to outlive the whole read.
        const uint32_t header_slen = S_CAST(uint32_t, strlen(header_line));
        char* name_store;
        if (unlikely(bigstack_end_alloc_c(header_slen + 1, &name_store))) {
          goto ScoreReport_ret_NOMEM;
        }
        memcpy(name_store, header_line, header_slen + 1);

        // Token 0 names the variant-ID column and is ignored; the rest name
        // partitions.
        uint32_t header_tok_ct = 0;
        {
          char* iter = FirstNonTspace(name_store);
          while (!IsEolnKns(*iter)) {
            ++header_tok_ct;
            iter = FirstNonTspace(CurTokenEnd(iter));
          }
        }
        if (unlikely(header_tok_ct < 2)) {
          logerrputs("Error: --score-partitions: no partition columns; expected a variant ID column\nfollowed by one 0/1 column per partition.\n");
          reterr = kPglRetMalformedInput;
          goto ScoreReport_ret_1;
        }
        const uint32_t part_ct = header_tok_ct - 1;
        qsr_ct = part_ct + 1;  // +1 for the total

        if (unlikely(bigstack_end_alloc_u32(raw_variant_ctl, &variant_include_cumulative_popcounts) ||
                     bigstack_end_calloc_w(BitCtToWordCt(S_CAST(uint64_t, qsr_ct) * variant_ct), &qsr_include) ||
                     bigstack_end_alloc_cp(qsr_ct, &range_names))) {
          goto ScoreReport_ret_NOMEM;
        }
        // Null-terminate each partition name in place within name_store.
        range_names[0] = nullptr;  // the total
        {
          char* iter = FirstNonTspace(name_store);
          iter = FirstNonTspace(CurTokenEnd(iter));  // skip the ID column name
          for (uint32_t part_idx = 0; part_idx != part_ct; ++part_idx) {
            char* tok_end = CurTokenEnd(iter);
            const uint32_t tok_slen = S_CAST(uint32_t, tok_end - iter);
            char* next = FirstNonTspace(tok_end);
            *tok_end = '\0';
            range_names[part_idx + 1] = iter;
            for (uint32_t prev = 0; prev != part_idx; ++prev) {
              if (unlikely(!strcmp(range_names[prev + 1], iter))) {
                logerrprintfww("Error: --score-partitions: duplicate partition name '%s'.\n", iter);
                reterr = kPglRetMalformedInput;
                goto ScoreReport_ret_1;
              }
            }
            if (unlikely(!tok_slen)) {
              logerrputs("Error: --score-partitions: empty partition name in the header.\n");
              reterr = kPglRetMalformedInput;
              goto ScoreReport_ret_1;
            }
            iter = next;
          }
        }
        FillCumulativePopcounts(variant_include, raw_variant_ctl, variant_include_cumulative_popcounts);

        // The total covers every variant the scoring pass will see.
        for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
          SetBit(S_CAST(uintptr_t, variant_idx) * qsr_ct, qsr_include);
        }

        const uint32_t variant_ctl = BitCtToWordCt(variant_ct);
        uintptr_t* already_seen;
        uint32_t* part_variant_cts;
        if (unlikely(bigstack_calloc_w(variant_ctl, &already_seen) ||
                     bigstack_alloc_u32(part_ct, &part_variant_cts))) {
          goto ScoreReport_ret_NOMEM;
        }
        ZeroU32Arr(part_ct, part_variant_cts);
        uintptr_t miss_ct = 0;
        line_idx = 1;
        while (1) {
          ++line_idx;
          char* line_start = TextGet(&score_txs);
          if (!line_start) {
            break;
          }
          line_start = FirstNonTspace(line_start);
          if (IsEolnKns(*line_start) || (*line_start == '#')) {
            continue;
          }
          char* varid_end = CurTokenEnd(line_start);
          const uint32_t varid_slen = S_CAST(uint32_t, varid_end - line_start);
          const uint32_t variant_uidx = VariantIdDupflagHtableFind(line_start, variant_ids, variant_id_htable, varid_slen, variant_id_htable_size, max_variant_id_slen);
          const uint32_t known = !((variant_uidx >> 31) || (!IsSet(variant_include, variant_uidx)));
          if (!known) {
            ++miss_ct;
          }
          // Parse and validate every row, known variant or not: a malformed
          // mask value is a malformed file regardless of whether that variant
          // happens to be in this dataset.
          uintptr_t bit_idx_base = 0;
          if (known) {
            const uint32_t variant_idx = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx);
            if (unlikely(IsSet(already_seen, variant_idx))) {
              *varid_end = '\0';
              logerrprintfww("Error: --score-partitions: duplicate variant ID '%s' on line %" PRIuPTR ".\n", line_start, line_idx);
              reterr = kPglRetMalformedInput;
              goto ScoreReport_ret_1;
            }
            SetBit(variant_idx, already_seen);
            bit_idx_base = S_CAST(uintptr_t, variant_idx) * qsr_ct;
          }
          char* iter = FirstNonTspace(varid_end);
          for (uint32_t part_idx = 0; part_idx != part_ct; ++part_idx) {
            if (unlikely(IsEolnKns(*iter))) {
              // Count partition columns, not total columns: comparing "3
              // columns" against "5 partitions" reads as an off-by-one when
              // the difference is really the leading variant-ID column.
              logerrprintfww("Error: --score-partitions: line %" PRIuPTR " of %s has %u partition column%s; the header names %u.\n", line_idx, score_info_ptr->partition_fname, part_idx, (part_idx == 1)? "" : "s", part_ct);
              reterr = kPglRetMalformedInput;
              goto ScoreReport_ret_1;
            }
            char* tok_end = CurTokenEnd(iter);
            // Strictly 0 or 1.  "nonzero means include" would let a stray value
            // change a score without saying so.
            const uint32_t tok_slen = S_CAST(uint32_t, tok_end - iter);
            const uint32_t is_one = (tok_slen == 1) && (*iter == '1');
            if (unlikely(!is_one && !((tok_slen == 1) && (*iter == '0')))) {
              *tok_end = '\0';
              logerrprintfww("Error: --score-partitions: value '%s' on line %" PRIuPTR " is not 0 or 1.\n", iter, line_idx);
              reterr = kPglRetMalformedInput;
              goto ScoreReport_ret_1;
            }
            if (is_one && known) {
              SetBit(bit_idx_base + part_idx + 1, qsr_include);
              part_variant_cts[part_idx] += 1;
            }
            iter = FirstNonTspace(tok_end);
          }
          if (unlikely(!IsEolnKns(*iter))) {
            logerrprintfww("Error: --score-partitions: line %" PRIuPTR " of %s has more columns than the header's %u partitions.\n", line_idx, score_info_ptr->partition_fname, part_ct);
            reterr = kPglRetMalformedInput;
            goto ScoreReport_ret_1;
          }
        }
        if (TextStreamErrcode2(&score_txs, &reterr)) {
          if (!strcmp(TextStreamError(&score_txs), kShortErrInteriorEmptyLine)) {
            logerrprintfww("Error: --score-partitions: blank line inside %s. Blank lines are allowed only at the end of the file.\n", score_info_ptr->partition_fname);
            reterr = kPglRetMalformedInput;
          } else {
            TextStreamErrPrint(score_info_ptr->partition_fname, &score_txs);
          }
          goto ScoreReport_ret_1;
        }
        // An empty partition would occupy a column and score nothing.
        for (uint32_t part_idx = 0; part_idx != part_ct; ++part_idx) {
          if (unlikely(!part_variant_cts[part_idx])) {
            logerrprintfww("Error: --score-partitions: partition '%s' contains no variant present in this dataset.\n", range_names[part_idx + 1]);
            reterr = kPglRetInconsistentInput;
            goto ScoreReport_ret_1;
          }
        }
        logprintf("--score-partitions: %u partition%s and %u variant%s loaded.\n", part_ct, (part_ct == 1)? "" : "s", PopcountWords(already_seen, variant_ctl), (PopcountWords(already_seen, variant_ctl) == 1)? "" : "s");
        if (miss_ct) {
          logerrprintf("Warning: %" PRIuPTR " variant ID%s in --score-partitions %s not in the dataset.\n", miss_ct, (miss_ct == 1)? "" : "s", (miss_ct == 1)? "is" : "are");
        }
      }
      reterr = TextRetarget(score_info_ptr->input_fname, &score_txs);
      if (unlikely(reterr)) {
        cur_input_fname = score_info_ptr->input_fname;
        goto ScoreReport_ret_TSTREAM_FAIL;
      }
      BigstackReset(bigstack_mark2);
      line_idx = 0;
    } else {
      reterr = SizeAndInitTextStream(score_info_ptr->input_fname, bigstack_left() / 8, 1, &score_txs);
      if (unlikely(reterr)) {
        cur_input_fname = score_info_ptr->input_fname;
        goto ScoreReport_ret_TSTREAM_FAIL;
      }
    }

    uint32_t varid_col_idx = score_info_ptr->varid_col_p1 - 1;
    uint32_t allele_col_idx = score_info_ptr->allele_col_p1 - 1;
    uint32_t relevant_col_idx_end = MAXV(varid_col_idx, allele_col_idx);
    uint32_t* score_col_idx_deltas = nullptr;
    uintptr_t score_col_ct = 1;
    if (!score_info_ptr->input_col_idx_range_list.name_ct) {
      // catch edge case
      const uint32_t score_col_idx = allele_col_idx + 1;
      if (unlikely(score_col_idx == varid_col_idx)) {
        logerrprintf("Error: --score%s variant ID column index matches a coefficient column index.\n", multi_input? "-list" : "");
        goto ScoreReport_ret_INVALID_CMDLINE;
      }
      relevant_col_idx_end += (score_col_idx > varid_col_idx);
      if (unlikely(bigstack_alloc_u32(1, &score_col_idx_deltas))) {
        goto ScoreReport_ret_NOMEM;
      }
      score_col_idx_deltas[0] = score_col_idx;
    } else {
      const uint32_t range_list_max = NumericRangeListMax(&(score_info_ptr->input_col_idx_range_list));
      if (range_list_max > relevant_col_idx_end) {
        relevant_col_idx_end = range_list_max;
      }
      unsigned char* bigstack_end_mark2 = g_bigstack_end;
      const uint32_t last_col_idxl = BitCtToWordCt(relevant_col_idx_end);
      uintptr_t* score_col_bitarr;
      if (unlikely(bigstack_end_calloc_w(last_col_idxl, &score_col_bitarr))) {
        goto ScoreReport_ret_NOMEM;
      }
      if (unlikely(NumericRangeListToBitarr(&(score_info_ptr->input_col_idx_range_list), relevant_col_idx_end, 1, 0, score_col_bitarr))) {
        goto ScoreReport_ret_MISSING_TOKENS;
      }
      if (unlikely(IsSet(score_col_bitarr, varid_col_idx))) {
        logerrprintf("Error: --score%s variant ID column index matches a coefficient column index.\n", multi_input? "-list" : "");
        goto ScoreReport_ret_INVALID_CMDLINE;
      }
      if (unlikely(IsSet(score_col_bitarr, allele_col_idx))) {
        logerrprintf("Error: --score%s allele column index matches a coefficient column index.\n", multi_input? "-list" : "");
        goto ScoreReport_ret_INVALID_CMDLINE;
      }
      score_col_ct = PopcountWords(score_col_bitarr, last_col_idxl);
      if (unlikely(bigstack_alloc_u32(score_col_ct, &score_col_idx_deltas))) {
        goto ScoreReport_ret_NOMEM;
      }
      uintptr_t col_uidx_base = 0;
      uintptr_t score_col_bitarr_bits = score_col_bitarr[0];
      for (uintptr_t score_col_idx = 0; score_col_idx != score_col_ct; ++score_col_idx) {
        const uint32_t col_uidx = BitIter1(score_col_bitarr, &col_uidx_base, &score_col_bitarr_bits);
        score_col_idx_deltas[score_col_idx] = col_uidx;
      }
      // now convert to deltas
      for (uintptr_t score_col_idx = score_col_ct - 1; score_col_idx; --score_col_idx) {
        score_col_idx_deltas[score_col_idx] -= score_col_idx_deltas[score_col_idx - 1];
      }
      BigstackEndReset(bigstack_end_mark2);
    }

    LlStr* fname_iter = nullptr;
    uint32_t infile_ct = 1;
    if (multi_input) {
      reterr = ScoreListLoad(score_info_ptr->input_fname, &score_txs, &fname_iter, &infile_ct);
      if (unlikely(reterr)) {
        goto ScoreReport_ret_1;
      }
    } else {
      const uint32_t fname_blen = strlen(score_info_ptr->input_fname) + 1;
      if (unlikely(bigstack_alloc_llstr(fname_blen, &fname_iter))) {
        goto ScoreReport_ret_NOMEM;
      }
      fname_iter->next = nullptr;
      memcpy(fname_iter->str, score_info_ptr->input_fname, fname_blen);
    }

    uintptr_t score_final_col_ct = score_col_ct;
    if (qsr_ct) {
      const uint64_t prod = S_CAST(uint64_t, qsr_ct) * score_col_ct;
      if (prod > 0x7fffffff) {
        // little point in supporting this even in large-matrix build
        logerrputs("Error: <--score[-list] column count> * <--q-score-range range count> too large.\n");
        goto ScoreReport_ret_INCONSISTENT_INPUT;
      }
#ifndef LAPACK_ILP64
      if (unlikely(prod > (0x7fffffff / kScoreVariantBlockSize))) {
        logerrputs("Error: <--score[-list] column count> * <--q-score-range range count> too large\nfor this " PROG_NAME_STR " build.  If this is really the computation you want, use a\n" PROG_NAME_STR "\nbuild with large-matrix support.\n");
        goto ScoreReport_ret_INCONSISTENT_INPUT;
      }
#endif
      score_final_col_ct = qsr_ct * score_col_ct;
#ifndef LAPACK_ILP64
    } else {
      if (unlikely(score_final_col_ct > (0x7fffffff / kScoreVariantBlockSize))) {
        logerrputs("Error: --score[-list] column count too large for this " PROG_NAME_STR " build.  If this is\nreally the computation you want, use a " PROG_NAME_STR " build with large-matrix support.\n");
        goto ScoreReport_ret_INCONSISTENT_INPUT;
      }
#endif
    }
    const uint32_t calc_thread_ct = MAXV(1, max_thread_ct - 1);
    uint32_t sample_shard_size;
    uint32_t sample_shard_ct;
    if (sample_ct <= calc_thread_ct * kScoreSampleShardSizeAlign) {
      sample_shard_size = kScoreSampleShardSizeAlign;
      sample_shard_ct = DivUp(sample_ct, kScoreSampleShardSizeAlign);
    } else {
      sample_shard_size = kScoreSampleShardSizeAlign * DivUp(sample_ct, calc_thread_ct * kScoreSampleShardSizeAlign);
      sample_shard_ct = DivUp(sample_ct, sample_shard_size);
    }
#ifndef LAPACK_ILP64
    if (unlikely(score_final_col_ct * sample_shard_size > 0x7fffffff)) {
      logerrputs("Error: <--score[-list] column count> * <sample count> too large for this " PROG_NAME_STR "\nbuild.  If this is really the computation you want, use a " PROG_NAME_STR " build with\nlarge-matrix support.\n");
      goto ScoreReport_ret_INCONSISTENT_INPUT;
    }
#endif
    const uint32_t sample_shard_ct_m1 = sample_shard_ct - 1;
    CalcScoreCtx ctx;
    ctx.variant_include = qsr_ct? variant_include : nullptr;
    ctx.variant_include_cumulative_popcounts = variant_include_cumulative_popcounts;
    ctx.qsr_include = qsr_include;
    ctx.score_final_col_ct = score_final_col_ct;
    ctx.sample_shard_size = sample_shard_size;
    ctx.sample_ct = sample_ct;
    ctx.flags = flags;
    ctx.qsr_ct = qsr_ct;
    ctx.cur_variant_batch_size = kScoreVariantBlockSize;
    if (unlikely(SetThreadCt(sample_shard_ct, &tg))) {
      goto ScoreReport_ret_NOMEM;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    const uint32_t sample_ctaw2 = NypCtToAlignedWordCt(sample_ct);
    const uint32_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
    const uintptr_t dosage_main_stride = RoundUpPow2(sample_ct, kDosagePerVec);
    const uint32_t acc1_vec_ct = BitCtToVecCt(sample_ct);
    const uint32_t acc4_vec_ct = acc1_vec_ct * 4;
    const uint32_t acc8_vec_ct = acc1_vec_ct * 8;
    const uint32_t write_score_avgs = (flags / kfScoreColScoreAvgs) & 1;
    const uint32_t write_score_sums = (flags / kfScoreColScoreSums) & 1;
    if (unlikely(S_CAST(uint64_t, score_col_ct) * infile_ct > 0x7fffff)) {
      goto ScoreReport_ret_TOO_MANY_COLUMNS;
    }
    const uintptr_t global_score_col_ct = score_col_ct * infile_ct;
    // --score-partitions writes every partition on one row, so the row holds
    // qsr_ct times as many values as a per-partition file would.  Sizing this
    // for one partition overruns the buffer.
    const uintptr_t merge_partitions = (score_info_ptr->partition_fname != nullptr);
    const uintptr_t out_col_mult = merge_partitions? qsr_ct : 1;
    const uintptr_t overflow_buf_size = RoundUpPow2((global_score_col_ct * out_col_mult * (write_score_avgs + write_score_sums) + pheno_ct) * 16 + 3 * kMaxIdSlen + kCompressStreamBlock + 64, kCacheline);
    if (unlikely(overflow_buf_size > S_CAST(uintptr_t, kMaxLongLine) + S_CAST(uintptr_t, kCompressStreamBlock))) {
      goto ScoreReport_ret_TOO_MANY_COLUMNS;
    }
    uintptr_t overflow_buf_alloc = overflow_buf_size;
    if (flags & (kfScoreZs | kfScoreListVariantsZs)) {
      overflow_buf_alloc += CstreamWkspaceReq(overflow_buf_size);
    }
    uintptr_t raw_allele_ct = 2 * raw_variant_ct;
    if (allele_idx_offsets) {
      raw_allele_ct = allele_idx_offsets[raw_variant_ct];
    }
    const uintptr_t raw_allele_ctl = BitCtToWordCt(raw_allele_ct);
    const uintptr_t qsr_ct_nz = qsr_ct + (qsr_ct == 0);
    uint32_t* sample_include_cumulative_popcounts = nullptr;
    uintptr_t* sex_nonmale_collapsed = nullptr;
    VecW* missing_diploid_accx = nullptr;
    VecW* missing_haploid_accx = nullptr;
    uint32_t* allele_ct_bases;
    uint32_t* male_allele_ct_decrs;
    uint32_t* nonfemale_allele_ct_incrs;
    uint32_t* variant_ct_rems;
    uint32_t* variant_hap_ct_rems;
    uintptr_t* already_seen_variants;
    uintptr_t* already_seen_alleles;
    double* tmpfile_buf;
    char* overflow_buf = nullptr;
    if (unlikely(bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw2, &ctx.genovecs[0]) ||
                 bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw2, &ctx.genovecs[1]) ||
                 bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw, &ctx.missing_bitvecs[0]) ||
                 bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw, &ctx.missing_bitvecs[1]) ||
                 bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw, &ctx.missing_male_bitvecs[0]) ||
                 bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw, &ctx.missing_male_bitvecs[1]) ||
                 bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw, &ctx.missing_nonfemale_bitvecs[0]) ||
                 bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw, &ctx.missing_nonfemale_bitvecs[1]) ||
                 bigstack_alloc_d(kScoreVariantBlockSize * score_final_col_ct, &(ctx.score_dense_coefs_cmaj[0])) ||
                 bigstack_alloc_d(kScoreVariantBlockSize * score_final_col_ct, &(ctx.score_dense_coefs_cmaj[1])) ||
                 bigstack_alloc_dp(sample_shard_ct, &ctx.sharded_final_scores_cmaj) ||
                 bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
                 bigstack_alloc_w(sample_ctl, &sex_nonmale_collapsed) ||
                 bigstack_end_clalloc_v(11 * acc4_vec_ct * qsr_ct_nz, &missing_diploid_accx) ||
                 bigstack_end_clalloc_v(11 * acc4_vec_ct * qsr_ct_nz, &missing_haploid_accx) ||
                 bigstack_alloc_u32(qsr_ct_nz, &allele_ct_bases) ||
                 bigstack_alloc_u32(qsr_ct_nz, &male_allele_ct_decrs) ||
                 bigstack_alloc_u32(qsr_ct_nz, &nonfemale_allele_ct_incrs) ||
                 bigstack_alloc_u32(qsr_ct_nz * 2, &variant_ct_rems) ||
                 bigstack_alloc_u32(qsr_ct_nz * 2, &variant_hap_ct_rems) ||
                 bigstack_alloc_u64p(sample_shard_ct, &ctx.ddosage_incrs) ||
                 bigstack_alloc_dp(sample_shard_ct, &ctx.dosages_vmaj) ||
                 bigstack_alloc_u64(qsr_ct_nz * sample_ct, &ctx.ddosage_sums) ||
                 bigstack_alloc_w(raw_variant_ctl, &already_seen_variants) ||
                 bigstack_alloc_w(raw_allele_ctl, &already_seen_alleles) ||
                 bigstack_end_clalloc_d(score_col_ct, &tmpfile_buf) ||
                 bigstack_end_clalloc_c(overflow_buf_alloc, &overflow_buf))) {
      goto ScoreReport_ret_NOMEM;
    }
    {
      uintptr_t cur_sample_shard_size = sample_shard_size;
      for (uint32_t shard_idx = 0; ; ++shard_idx) {
        if (shard_idx >= sample_shard_ct_m1) {
          if (shard_idx > sample_shard_ct_m1) {
            break;
          }
          cur_sample_shard_size = sample_ct - sample_shard_ct_m1 * sample_shard_size;
        }
        if (unlikely(bigstack_alloc_d(score_final_col_ct * cur_sample_shard_size, &ctx.sharded_final_scores_cmaj[shard_idx]) ||
                     bigstack_alloc_u64(cur_sample_shard_size, &ctx.ddosage_incrs[shard_idx]) ||
                     bigstack_alloc_d(kScoreVariantBlockSize * cur_sample_shard_size, &ctx.dosages_vmaj[shard_idx]))) {
          goto ScoreReport_ret_NOMEM;
        }
      }
    }
    if (qsr_ct) {
      if (unlikely(bigstack_alloc_u32(kScoreVariantBlockSize, &ctx.variant_uidxs[0]) ||
                   bigstack_alloc_u32(kScoreVariantBlockSize, &ctx.variant_uidxs[1]))) {
        goto ScoreReport_ret_NOMEM;
      }
    } else {
      ctx.variant_uidxs[0] = nullptr;
      ctx.variant_uidxs[1] = nullptr;
    }
    if (PgrGetGflags(simple_pgrp) & kfPgenGlobalDosagePresent) {
      // We may want to buffer less than 2x240 variants at a time in this case
      // if there are >500k samples...
      if (unlikely(bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw, &ctx.dosage_presents[0]) ||
                   bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw, &ctx.dosage_presents[1]) ||
                   bigstack_alloc_dosage(kScoreVariantBlockSize * dosage_main_stride, &ctx.dosage_mains[0]) ||
                   bigstack_alloc_dosage(kScoreVariantBlockSize * dosage_main_stride, &ctx.dosage_mains[1]) ||
                   bigstack_alloc_u32(kScoreVariantBlockSize, &ctx.dosage_cts[0]) ||
                   bigstack_alloc_u32(kScoreVariantBlockSize, &ctx.dosage_cts[1]))) {
        goto ScoreReport_ret_NOMEM;
      }
    } else {
      ctx.dosage_presents[0] = nullptr;
      ctx.dosage_presents[1] = nullptr;
      ctx.dosage_mains[0] = nullptr;
      ctx.dosage_mains[1] = nullptr;
      ctx.dosage_cts[0] = nullptr;
      ctx.dosage_cts[1] = nullptr;
    }
    const uint32_t domrec = !!(flags & (kfScoreDominant | kfScoreRecessive));
    uint32_t* common_geno_sum_incrs = nullptr;
    double* common_score_incrs = nullptr;
    uint32_t max_difflist_len = 0;
    uintptr_t raregeno_stride = 0;
    uintptr_t difflist_sample_ids_stride = 0;
    if ((!ctx.dosage_presents[0]) && (!domrec)) {
      max_difflist_len = (score_final_col_ct == 1)? (sample_ct / 16) : (sample_ct / 32);
      raregeno_stride = NypCtToAlignedWordCt(max_difflist_len);
      difflist_sample_ids_stride = RoundUpPow2(max_difflist_len, kInt32PerVec);
      if (unlikely(bigstack_alloc_u32(qsr_ct_nz, &common_geno_sum_incrs) ||
                   bigstack_alloc_d(score_final_col_ct, &common_score_incrs) ||
                   bigstack_alloc_w(kScoreVariantBlockSize * raregeno_stride, &ctx.raregenos[0]) ||
                   bigstack_alloc_w(kScoreVariantBlockSize * raregeno_stride, &ctx.raregenos[1]) ||
                   bigstack_alloc_u32(kScoreVariantBlockSize, &ctx.difflist_lens[0]) ||
                   bigstack_alloc_u32(kScoreVariantBlockSize, &ctx.difflist_lens[1]) ||
                   bigstack_alloc_u32(kScoreVariantBlockSize * difflist_sample_ids_stride, &ctx.difflist_sample_ids[0]) ||
                   bigstack_alloc_u32(kScoreVariantBlockSize * difflist_sample_ids_stride, &ctx.difflist_sample_ids[1]) ||
                   bigstack_alloc_i8(kScoreVariantBlockSize, &ctx.difflist_common_genos[0]) ||
                   bigstack_alloc_i8(kScoreVariantBlockSize, &ctx.difflist_common_genos[1]) ||
                   bigstack_alloc_d(kScoreVariantBlockSize * score_final_col_ct, &ctx.score_sparse_coefs_vmaj[0]) ||
                   bigstack_alloc_d(kScoreVariantBlockSize * score_final_col_ct, &ctx.score_sparse_coefs_vmaj[1]))) {
        goto ScoreReport_ret_NOMEM;
      }
    } else {
      ctx.raregenos[0] = nullptr;
      ctx.raregenos[1] = nullptr;
      ctx.difflist_lens[0] = nullptr;
      ctx.difflist_lens[1] = nullptr;
      ctx.difflist_sample_ids[0] = nullptr;
      ctx.difflist_sample_ids[1] = nullptr;
      ctx.difflist_common_genos[0] = nullptr;
      ctx.difflist_common_genos[1] = nullptr;
      ctx.score_sparse_coefs_vmaj[0] = nullptr;
      ctx.score_sparse_coefs_vmaj[1] = nullptr;
    }
    ctx.max_difflist_len = max_difflist_len;
    if (allele_freqs) {
      if (unlikely(bigstack_alloc_d(kScoreVariantBlockSize, &ctx.allele_freqs[0]) ||
                   bigstack_alloc_d(kScoreVariantBlockSize, &ctx.allele_freqs[1]) ||
                   bigstack_alloc_d(kScoreVariantBlockSize, &ctx.geno_slopes[0]) ||
                   bigstack_alloc_d(kScoreVariantBlockSize, &ctx.geno_slopes[1]) ||
                   bigstack_alloc_d(kScoreVariantBlockSize, &ctx.geno_intercepts[0]) ||
                   bigstack_alloc_d(kScoreVariantBlockSize, &ctx.geno_intercepts[1]))) {
        goto ScoreReport_ret_NOMEM;
      }
    } else {
      ctx.allele_freqs[0] = nullptr;
      ctx.allele_freqs[1] = nullptr;
      ctx.geno_slopes[0] = nullptr;
      ctx.geno_slopes[1] = nullptr;
      ctx.geno_intercepts[0] = nullptr;
      ctx.geno_intercepts[1] = nullptr;
    }
    SetThreadFuncAndData(CalcScoreThread, &ctx, &tg);

    FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
    CopyBitarrSubset(sex_male, sample_include, sample_ct, sex_nonmale_collapsed);
    AlignedBitarrInvert(sample_ct, sex_nonmale_collapsed);
    ctx.sex_nonmale_collapsed = sex_nonmale_collapsed;
    const uintptr_t* sex_nonfemale = sex_male;
    uintptr_t* sex_female_collapsed = sex_nonmale_collapsed;
    if (nosex_ct) {
      uintptr_t* sex_female_tmp;
      if (unlikely(bigstack_alloc_w(raw_sample_ctl, &sex_female_tmp) ||
                   bigstack_alloc_w(sample_ctl, &sex_female_collapsed))) {
        goto ScoreReport_ret_NOMEM;
      }
      BitvecInvmaskCopy(sex_nm, sex_male, raw_sample_ctl, sex_female_tmp);
      CopyBitarrSubset(sex_female_tmp, sample_include, sample_ct, sex_female_collapsed);
      // now invert to nonfemale
      BitvecInvert(raw_sample_ctl, sex_female_tmp);
      // (not actually necessary here, but let's keep sex_nonfemale definition
      // consistent since failure to do so has resulted in a bug)
      BitvecAnd(sample_include, raw_sample_ctl, sex_female_tmp);
      sex_nonfemale = sex_female_tmp;
    }
    ctx.sex_female_collapsed = sex_nonmale_collapsed;

    if (!variant_id_htable) {
      reterr = AllocAndPopulateIdHtableMt(variant_include, variant_ids, variant_ct, 0, max_thread_ct, &variant_id_htable, nullptr, &variant_id_htable_size, nullptr);
      if (unlikely(reterr)) {
        goto ScoreReport_ret_1;
      }
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".sscore.tmp");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &score_tmpfile))) {
      goto ScoreReport_ret_OPEN_FAIL;
    }

    const uint32_t ignore_dup_ids = (flags / kfScoreIgnoreDupIds) & 1;
    const uint32_t list_variants = (flags / kfScoreListVariants) & 1;
    if (list_variants) {
      const uint32_t list_variants_zst = (flags / kfScoreListVariantsZs) & 1;
      OutnameZstSet(".sscore.vars", list_variants_zst, outname_end);
      reterr = InitCstream(outname, 0, list_variants_zst, 1, overflow_buf_size, overflow_buf, R_CAST(unsigned char*, &(overflow_buf[overflow_buf_size])), &css);
      if (unlikely(reterr)) {
        goto ScoreReport_ret_1;
      }
      cswritep = overflow_buf;
    }

    char** score_col_names;
    if (unlikely(bigstack_end_clalloc_cp(global_score_col_ct, &score_col_names))) {
      goto ScoreReport_ret_NOMEM;
    }
    if (!(flags & kfScoreHeaderRead)) {
      // just perform this check once, rather than once per file.
      // errored out with "Too many score columns" earlier if LHS > 2^31.
      if (global_score_col_ct * 16 > bigstack_left()) {
        goto ScoreReport_ret_NOMEM;
      }
    }

#ifdef USE_MTBLAS
    if (sample_shard_ct < calc_thread_ct) {
      const uint32_t matrix_multiply_thread_ct = 1 + ((calc_thread_ct - 1) / sample_shard_ct);
      BLAS_SET_NUM_THREADS(matrix_multiply_thread_ct);
    }
#endif

    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi);
    const uint32_t lines_to_skip_p1 = 1 + ((flags / kfScoreHeaderIgnore) & 1);
    // Stricter threshold for multi-score since standard matrix-multiply
    // workflow is more likely to pay off
    const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
    const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
    const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    const uint32_t variance_standardize = (flags / kfScoreVarianceStandardize) & 1;
    const uint32_t center = variance_standardize || (flags & kfScoreCenter);
    const uint32_t no_meanimpute = (flags / kfScoreNoMeanimpute) & 1;
    const uint32_t y_prohibited = nosex_ct && (!no_meanimpute);
    const uint32_t se_mode = (flags / kfScoreSe) & 1;
    const uintptr_t tmpfile_buf_byte_ct = score_col_ct * sizeof(double);
    // If we're not attempting the sparse-genotype optimization, this is always
    // UINT32_MAX.
    uint32_t difflist_common_geno = UINT32_MAX;
    uint32_t cur_allele_ct = 2;
    uint32_t dosage_ct = 0;
    double geno_slope = kRecipDosageMax;
    double geno_intercept = 0.0;
    char* score_name_write_iter = R_CAST(char*, g_bigstack_end);
    for (uint32_t file_idx1 = 1; file_idx1 <= infile_ct; ++file_idx1, fname_iter = fname_iter->next) {
      cur_input_fname = fname_iter->str;
      if (file_idx1 > 1) {
        reterr = TextRetarget(cur_input_fname, &score_txs);
        if (unlikely(reterr)) {
          goto ScoreReport_ret_TSTREAM_FAIL;
        }
      }
      ZeroVecArr(11 * acc4_vec_ct * qsr_ct_nz, missing_diploid_accx);
      ZeroVecArr(11 * acc4_vec_ct * qsr_ct_nz, missing_haploid_accx);
      ZeroU32Arr(qsr_ct_nz, allele_ct_bases);
      ZeroU32Arr(qsr_ct_nz, male_allele_ct_decrs);
      ZeroU32Arr(qsr_ct_nz, nonfemale_allele_ct_incrs);
      ZeroU64Arr(qsr_ct_nz * sample_ct, ctx.ddosage_sums);
      ZeroWArr(raw_variant_ctl, already_seen_variants);
      ZeroWArr(raw_allele_ctl, already_seen_alleles);
      {
        uintptr_t cur_sample_shard_size = sample_shard_size;
        for (uint32_t shard_idx = 0; ; ++shard_idx) {
          if (shard_idx >= sample_shard_ct_m1) {
            if (shard_idx > sample_shard_ct_m1) {
              break;
            }
            cur_sample_shard_size = sample_ct - sample_shard_ct_m1 * sample_shard_size;
          }
          ZeroDArr(score_final_col_ct * cur_sample_shard_size, ctx.sharded_final_scores_cmaj[shard_idx]);
        }
      }
      if (common_geno_sum_incrs) {
        ZeroU32Arr(qsr_ct_nz, common_geno_sum_incrs);
        ZeroDArr(score_final_col_ct, common_score_incrs);
      }
      for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct_nz; ++qsr_idx) {
        variant_ct_rems[2 * qsr_idx] = 15;
        variant_ct_rems[2 * qsr_idx + 1] = 17;
        variant_hap_ct_rems[2 * qsr_idx] = 15;
        variant_hap_ct_rems[2 * qsr_idx + 1] = 17;
      }
      ctx.cur_variant_batch_size = kScoreVariantBlockSize;
      line_idx = 0;
      char* line_start;
      for (uint32_t uii = 0; uii != lines_to_skip_p1; ++uii) {
        ++line_idx;
        line_start = TextGet(&score_txs);
        if (unlikely(!line_start)) {
          if (!TextStreamErrcode2(&score_txs, &reterr)) {
            logerrprintf("Error: --score%s: %s is empty.\n", multi_input? "-list" : "", cur_input_fname);
            goto ScoreReport_ret_MALFORMED_INPUT;
          }
          goto ScoreReport_ret_TSTREAM_FAIL;
        }
      }
      {
        const uintptr_t global_score_col_idx_start = (file_idx1 - 1) * score_col_ct;
        if (flags & kfScoreHeaderRead) {
          char* line_end = TextLineEnd(&score_txs);
          if (S_CAST(uintptr_t, line_end - line_start) > bigstack_left()) {
            goto ScoreReport_ret_NOMEM;
          }
          char** cur_score_col_names = &(score_col_names[global_score_col_idx_start]);
          char* read_iter = line_start;
          for (uintptr_t score_col_idx = 0; score_col_idx != score_col_ct; ++score_col_idx) {
            read_iter = NextTokenMult0(read_iter, score_col_idx_deltas[score_col_idx]);
            if (unlikely(!read_iter)) {
              goto ScoreReport_ret_MISSING_TOKENS;
            }
            char* token_end = CurTokenEnd(read_iter);
            const uint32_t slen = token_end - read_iter;
            score_name_write_iter -= 1 + slen;
            cur_score_col_names[score_col_idx] = score_name_write_iter;
            memcpyx(score_name_write_iter, read_iter, slen, '\0');
          }
          if (file_idx1 == infile_ct) {
            BigstackEndSet(score_name_write_iter);
            bigstack_end_clalign();
            unsigned char* bigstack_mark2 = g_bigstack_base;
            uintptr_t* col_subset_mask;
            if (unlikely(bigstack_alloc_w(BitCtToWordCt(score_col_ct), &col_subset_mask))) {
              goto ScoreReport_ret_NOMEM;
            }
            SetAllBits(score_col_ct, col_subset_mask);
            uint32_t dup_found;
            reterr = CheckIdUniqueness(g_bigstack_base, g_bigstack_end, col_subset_mask, TO_CONSTCPCONSTP(score_col_names), score_col_ct, max_thread_ct, &dup_found);
            if (unlikely(reterr)) {
              goto ScoreReport_ret_1;
            }
            if (unlikely(dup_found)) {
              logerrprintf("Error: --score%s: Score IDs are not unique.\n", multi_input? "-list" : "");
              goto ScoreReport_ret_MALFORMED_INPUT;
            }
            BigstackReset(bigstack_mark2);
          }
        } else {
          const uintptr_t col_idx_stop = global_score_col_idx_start + score_col_ct;
          for (uintptr_t global_score_col_idx = global_score_col_idx_start; global_score_col_idx != col_idx_stop; ++global_score_col_idx) {
            const uint32_t str_blen = 6 + UintSlen(global_score_col_idx + 1);
            score_name_write_iter -= str_blen;
            score_col_names[global_score_col_idx] = score_name_write_iter;
            char* tmp_write_iter = strcpya_k(score_name_write_iter, "SCORE");
            u32toa_x(global_score_col_idx + 1, '\0', tmp_write_iter);
          }
          if (file_idx1 == infile_ct) {
            BigstackEndSet(score_name_write_iter);
            bigstack_end_clalign();
          }
        }
      }

      ReinitThreads(&tg);
      uint32_t block_vidx = 0;
      uint32_t sparse_vidx = 0;
      uint32_t parity = 0;
      uint32_t prev_variant_uidx = UINT32_MAX;
      uintptr_t* genovec_iter = ctx.genovecs[0];
      uintptr_t* raregeno_iter = ctx.raregenos[0];
      uint32_t* difflist_lens = ctx.difflist_lens[0];
      uint32_t* difflist_sample_ids_iter = ctx.difflist_sample_ids[0];
      int8_t* difflist_common_genos = ctx.difflist_common_genos[0];
      uintptr_t* dosage_present_iter = ctx.dosage_presents[0];
      Dosage* dosage_main_iter = ctx.dosage_mains[0];
      uint32_t* dosage_cts = ctx.dosage_cts[0];
      uintptr_t* missing_bitvec_iter = ctx.missing_bitvecs[0];
      uintptr_t* missing_male_bitvec_iter = ctx.missing_male_bitvecs[0];
      uintptr_t* missing_nonfemale_bitvec_iter = ctx.missing_nonfemale_bitvecs[0];
      double* cur_allele_freqs = ctx.allele_freqs[0];
      double* geno_slopes = ctx.geno_slopes[0];
      double* geno_intercepts = ctx.geno_intercepts[0];
      double* score_dense_coefs_cmaj = ctx.score_dense_coefs_cmaj[0];
      double* score_sparse_coefs_vmaj = ctx.score_sparse_coefs_vmaj[0];
      uint32_t* variant_uidxs = ctx.variant_uidxs[0];
      unsigned char* is_nonx_haploids = &(ctx.is_nonx_haploids[0][0]);
      unsigned char* is_relevant_xs = &(ctx.is_relevant_xs[0][0]);
      unsigned char* is_ys = &(ctx.is_ys[0][0]);
      unsigned char* ploidy_m1s = &(ctx.ploidy_m1s[0][0]);
      uint32_t valid_variant_ct = 0;
      uint32_t difflist_len = 0;
      uintptr_t missing_var_id_ct = 0;
      uintptr_t missing_allele_code_ct = 0;
      uintptr_t duplicated_var_id_ct = 0;
      if (flags & kfScoreHeaderRead) {
        ++line_idx;
        line_start = TextGet(&score_txs);
      }
      for (; line_start; ++line_idx, line_start = TextGet(&score_txs)) {
        // varid_col_idx and allele_col_idx will almost always be very small
        char* variant_id_start = NextTokenMult0(line_start, varid_col_idx);
        if (unlikely(!variant_id_start)) {
          goto ScoreReport_ret_MISSING_TOKENS;
        }
        char* variant_id_token_end = CurTokenEnd(variant_id_start);
        const uint32_t variant_id_slen = variant_id_token_end - variant_id_start;
        uint32_t variant_uidx = VariantIdDupflagHtableFind(variant_id_start, variant_ids, variant_id_htable, variant_id_slen, variant_id_htable_size, max_variant_id_slen);
        if (variant_uidx >> 31) {
          ++missing_var_id_ct;
          if (variant_uidx != UINT32_MAX) {
            if (unlikely(!ignore_dup_ids)) {
              snprintf(g_logbuf, kLogbufSize, "Error: --score%s variant ID '%s' appears multiple times in main dataset.\n", multi_input? "-list" : "", variant_ids[variant_uidx & 0x7fffffff]);
              goto ScoreReport_ret_INCONSISTENT_INPUT_WW;
            }
            ++duplicated_var_id_ct;
            // subtract this from missing_var_id_ct later
          }
          continue;
        }
        char* allele_start = NextTokenMult0(line_start, allele_col_idx);
        if (unlikely(!allele_start)) {
          goto ScoreReport_ret_MISSING_TOKENS;
        }
        uintptr_t allele_idx_offset_base;
        if (!allele_idx_offsets) {
          allele_idx_offset_base = variant_uidx * 2;
        } else {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
        }
        char* allele_end = CurTokenEnd(allele_start);
        const char allele_end_char = *allele_end;
        *allele_end = '\0';
        const uint32_t allele_blen = 1 + S_CAST(uintptr_t, allele_end - allele_start);
        const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);

        uint32_t cur_aidx = 0;
        for (; cur_aidx != cur_allele_ct; ++cur_aidx) {
          if (memequal(allele_start, cur_alleles[cur_aidx], allele_blen)) {
            break;
          }
        }
        // compiler is smart enough to avoid repeating this test
        if (cur_aidx == cur_allele_ct) {
          ++missing_allele_code_ct;
          *allele_end = allele_end_char;
          continue;
        }
        const uintptr_t allele_idx = allele_idx_offset_base + cur_aidx;
        if (unlikely(IsSet(already_seen_alleles, allele_idx))) {
          char* errwrite_iter = strcpya_k(g_logbuf, "Error: --score");
          if (multi_input) {
            errwrite_iter = strcpya_k(errwrite_iter, "-list");
          }
          errwrite_iter = strcpya_k(errwrite_iter, ": ");
          // Don't write allele code, since it might be too long for the
          // buffer.
          if (!cur_aidx) {
            errwrite_iter = strcpya_k(errwrite_iter, "REF");
          } else {
            errwrite_iter = strcpya_k(errwrite_iter, "ALT");
            errwrite_iter = u32toa(cur_aidx, errwrite_iter);
          }
          errwrite_iter = strcpya_k(errwrite_iter, " allele for variant '");
          errwrite_iter = strcpya(errwrite_iter, variant_ids[variant_uidx]);
          errwrite_iter = strcpya_k(errwrite_iter, "' appears multiple times in ");
          errwrite_iter = strcpya(errwrite_iter, cur_input_fname);
          if (multi_input) {
            errwrite_iter = strcpya_k(errwrite_iter, "-list");
          }
          strcpy_k(errwrite_iter, " file.\n");
          goto ScoreReport_ret_MALFORMED_INPUT_WW;
        }
        SetBit(allele_idx, already_seen_alleles);
        const uint32_t is_new_variant = 1 - IsSet(already_seen_variants, variant_uidx);
        SetBit(variant_uidx, already_seen_variants);

        // okay, the variant and allele are in our dataset.
        const uint32_t chr_idx = GetVariantChr(cip, variant_uidx);
        uint32_t is_nonx_haploid = IsSet(cip->haploid_mask, chr_idx);
        if (unlikely(domrec && is_nonx_haploid)) {
          fputs("\n", stdout);
          logerrputs("Error: --score[-list] 'dominant' and 'recessive' modifiers cannot be used with\nhaploid chromosomes.\n");
          goto ScoreReport_ret_INCONSISTENT_INPUT;
        }
        uint32_t is_relevant_x = (chr_idx == x_code);
        if (unlikely(variance_standardize && (is_relevant_x || (chr_idx == mt_code)))) {
          fputs("\n", stdout);
          logerrputs("Error: --score[-list] 'variance-standardize' cannot be used with chrX or MT.\n");
          goto ScoreReport_ret_INCONSISTENT_INPUT;
        }
        is_nonx_haploid = (!is_relevant_x) && is_nonx_haploid;

        // only if --xchr-model 1 (which is no longer the default)
        is_relevant_x = is_relevant_x && xchr_model;

        const uint32_t is_y = (chr_idx == y_code);
        if (unlikely(is_y && y_prohibited)) {
          fputs("\n", stdout);
          logerrputs("Error: When both chrY variants and unknown-sex samples are present,\n--score[-list] can only be run with the 'no-mean-imputation' modifier.\n");
          goto ScoreReport_ret_INCONSISTENT_INPUT;
        }
        uint32_t ploidy_m1;
        if ((variant_uidx == prev_variant_uidx) && (cur_allele_ct == 2)) {
          // Previous record is the exact inverse of what we now want, so just
          // copy and invert it.
          // This is worth the trouble since PgrGet1D() is usually the main
          // bottleneck, and the repeated-variant case definitely comes up from
          // at least the recommended PCA-projection workflow.
          if (difflist_common_geno != UINT32_MAX) {
            memcpy(raregeno_iter, &(raregeno_iter[-S_CAST(intptr_t, raregeno_stride)]), NypCtToWordCt(difflist_len) * sizeof(intptr_t));
            memcpy(difflist_sample_ids_iter, &(difflist_sample_ids_iter[-S_CAST(intptr_t, difflist_sample_ids_stride)]), difflist_len * sizeof(int32_t));
            GenovecInvertUnsafe(difflist_len, raregeno_iter);
            difflist_common_geno = (6 - difflist_common_geno) & 3;
          } else {
            memcpy(genovec_iter, &(genovec_iter[-S_CAST(intptr_t, sample_ctaw2)]), sample_ctaw2 * sizeof(intptr_t));
            GenovecInvertUnsafe(sample_ct, genovec_iter);
            if (dosage_present_iter) {
              memcpy(dosage_present_iter, &(dosage_present_iter[-S_CAST(intptr_t, sample_ctaw)]), sample_ctaw * sizeof(intptr_t));
              dosage_ct = dosage_cts[block_vidx - 1];
              memcpy(dosage_main_iter, &(dosage_main_iter[-S_CAST(intptr_t, dosage_main_stride)]), dosage_ct * sizeof(Dosage));
              BiallelicDosage16Invert(dosage_ct, dosage_main_iter);
            }
          }
        } else if (dosage_present_iter || (cur_allele_ct > 2) || domrec || is_nonx_haploid || is_relevant_x) {
          difflist_common_geno = UINT32_MAX;
          reterr = PgrGet1D(sample_include, pssi, sample_ct, variant_uidx, cur_aidx, simple_pgrp, genovec_iter, dosage_present_iter, dosage_main_iter, &dosage_ct);
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, variant_uidx);
            goto ScoreReport_ret_1;
          }
        } else {
          // Sparse-genotype optimization.  We limit this to the no-dosage,
          // all-diploid, biallelic, no 'dominant'/'recessive' case to keep
          // code size under control.
          reterr = PgrGetDifflistOrGenovec(sample_include, pssi, sample_ct, max_difflist_len, variant_uidx, simple_pgrp, genovec_iter, &difflist_common_geno, raregeno_iter, difflist_sample_ids_iter, &difflist_len);
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, variant_uidx);
            goto ScoreReport_ret_1;
          }
          if (!cur_aidx) {
            if (difflist_common_geno == UINT32_MAX) {
              GenovecInvertUnsafe(sample_ct, genovec_iter);
            } else {
              GenovecInvertUnsafe(difflist_len, raregeno_iter);
              difflist_common_geno = (6 - difflist_common_geno) & 3;
            }
          }
        }
        if (difflist_common_geno != UINT32_MAX) {
          ZeroTrailingNyps(difflist_len, raregeno_iter);
          if (is_new_variant) {
            uint32_t missing_exists = (difflist_common_geno == 3);
            if ((!missing_exists) && difflist_len) {
              STD_ARRAY_DECL(uint32_t, 4, genocounts);
              GenoarrCountFreqsUnsafe(raregeno_iter, difflist_len, genocounts);
              missing_exists = (genocounts[3] != 0);
            }
            // This must be kept in sync with how the main branch processes
            // missingness.
            if (missing_exists) {
              SparseToMissingness(raregeno_iter, difflist_sample_ids_iter, sample_ct, difflist_common_geno, difflist_len, missing_bitvec_iter);
            }
            if (!qsr_ct) {
              allele_ct_bases[0] += 2;
              if (missing_exists) {
                VerticalCounterUpdate(missing_bitvec_iter, acc1_vec_ct, variant_ct_rems, missing_diploid_accx);
              }
            } else {
              const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
              for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                if (IsSet(qsr_include, qsr_idx + bit_idx_base)) {
                  allele_ct_bases[qsr_idx] += 2;
                  if (missing_exists) {
                    VerticalCounterUpdate(missing_bitvec_iter, acc1_vec_ct, &(variant_ct_rems[2 * qsr_idx]), &(missing_diploid_accx[acc4_vec_ct * 11 * qsr_idx]));
                  }
                }
              }
            }
          }
          // Enable sparse ddosage_sums update.
          if ((difflist_common_geno == 1) || (difflist_common_geno == 2)) {
            if (!qsr_ct) {
              common_geno_sum_incrs[0] += difflist_common_geno;
            } else {
              const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
              for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                if (IsSet(qsr_include, qsr_idx + bit_idx_base)) {
                  common_geno_sum_incrs[qsr_idx] += difflist_common_geno;
                }
              }
            }
          }
          ploidy_m1 = 1;
        } else {
          ZeroTrailingNyps(sample_ct, genovec_iter);
          GenoarrToMissingnessUnsafe(genovec_iter, sample_ct, missing_bitvec_iter);
          if (dosage_ct) {
            BitvecInvmask(dosage_present_iter, sample_ctl, missing_bitvec_iter);
          }
          if (is_nonx_haploid) {
            if (is_y) {
              if (is_new_variant) {
                if (!qsr_ct) {
                  nonfemale_allele_ct_incrs[0] += 1;
                } else {
                  const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
                  for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                    if (IsSet(qsr_include, qsr_idx + bit_idx_base)) {
                      nonfemale_allele_ct_incrs[qsr_idx] += 1;
                    }
                  }
                }
              }
              BitvecInvmask(sex_female_collapsed, sample_ctl, missing_bitvec_iter);
            } else if (is_new_variant) {
              if (!qsr_ct) {
                allele_ct_bases[0] += 1;
              } else {
                const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
                for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                  if (IsSet(qsr_include, qsr_idx + bit_idx_base)) {
                    allele_ct_bases[qsr_idx] += 1;
                  }
                }
              }
            }
            if (is_new_variant) {
              if (!qsr_ct) {
                VerticalCounterUpdate(missing_bitvec_iter, acc1_vec_ct, variant_hap_ct_rems, missing_haploid_accx);
              } else {
                const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
                for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                  if (IsSet(qsr_include, qsr_idx + bit_idx_base)) {
                    VerticalCounterUpdate(missing_bitvec_iter, acc1_vec_ct, &(variant_ct_rems[2 * qsr_idx]), &(missing_haploid_accx[acc4_vec_ct * 11 * qsr_idx]));
                  }
                }
              }
            }
            if (is_y) {
              memcpy(missing_nonfemale_bitvec_iter, missing_bitvec_iter, sample_ctl * sizeof(intptr_t));
              BitvecOr(sex_female_collapsed, sample_ctl, missing_bitvec_iter);
            }
            ploidy_m1 = 0;
          } else {
            // !is_nonx_haploid
            if (is_relevant_x) {
              // Handle nonmales first.  Remove males from missing_bitvec_iter,
              // to be restored at the end of this block.
              BitvecInvmaskCopy(missing_bitvec_iter, sex_nonmale_collapsed, sample_ctl, missing_male_bitvec_iter);
              BitvecAnd(sex_nonmale_collapsed, sample_ctl, missing_bitvec_iter);
            }
            if (is_new_variant) {
              if (!qsr_ct) {
                allele_ct_bases[0] += 2;
                VerticalCounterUpdate(missing_bitvec_iter, acc1_vec_ct, variant_ct_rems, missing_diploid_accx);
              } else {
                const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
                for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                  if (IsSet(qsr_include, qsr_idx + bit_idx_base)) {
                    allele_ct_bases[qsr_idx] += 2;
                    VerticalCounterUpdate(missing_bitvec_iter, acc1_vec_ct, &(variant_ct_rems[2 * qsr_idx]), &(missing_diploid_accx[acc4_vec_ct * 11 * qsr_idx]));
                  }
                }
              }
            }
            if (is_relevant_x) {
              // Now count males.
              if (is_new_variant) {
                if (!qsr_ct) {
                  male_allele_ct_decrs[0] += 1;
                  VerticalCounterUpdate(missing_male_bitvec_iter, acc1_vec_ct, variant_hap_ct_rems, missing_haploid_accx);
                } else {
                  const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
                  for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                    if (IsSet(qsr_include, qsr_idx + bit_idx_base)) {
                      male_allele_ct_decrs[qsr_idx] += 1;
                      VerticalCounterUpdate(missing_male_bitvec_iter, acc1_vec_ct, &(variant_hap_ct_rems[2 * qsr_idx]), &(missing_haploid_accx[acc4_vec_ct * 11 * qsr_idx]));
                    }
                  }
                }
              }
              // Restore missing males to missing_bitvec_iter.
              BitvecOr(missing_male_bitvec_iter, sample_ctl, missing_bitvec_iter);
            }
            ploidy_m1 = domrec? 0 : 1;
          }
        }
        const double ploidy_d = ploidy_m1? 2.0 : 1.0;
        if (allele_freqs) {
          const double cur_allele_freq = GetAlleleFreq(&(allele_freqs[allele_idx_offset_base - variant_uidx]), cur_aidx, cur_allele_ct);
          if (center) {
            // Note that we prohibit 'dominant'/'recessive' from being used
            // with 'center'/'variance-standardize'.
            if (variance_standardize) {
              const double variance = ploidy_d * 0.5 * ComputeDiploidMultiallelicVariance(&(allele_freqs[allele_idx_offset_base - variant_uidx]), cur_allele_ct);
              if (!(variance > kSmallEpsilon)) {
                // ZeroTrailingNyps(sample_ct, genovec_buf);
                STD_ARRAY_DECL(uint32_t, 4, genocounts);
                if (difflist_common_geno != UINT32_MAX) {
                  GenoarrCountFreqsUnsafe(raregeno_iter, difflist_len, genocounts);
                  genocounts[difflist_common_geno] = sample_ct - difflist_len;
                } else {
                  GenoarrCountFreqsUnsafe(genovec_iter, sample_ct, genocounts);
                }
                if (unlikely(dosage_ct || genocounts[1] || genocounts[2])) {
                  snprintf(g_logbuf, kLogbufSize, "Error: --score[-list] variance-standardize failure for variant '%s': estimated allele frequency is zero or NaN, but not all dosages are zero. (This is possible when e.g. allele frequencies are estimated from founders, but the allele is only observed in nonfounders.)\n", variant_ids[variant_uidx]);
                  goto ScoreReport_ret_DEGENERATE_DATA_WW;
                }
                geno_slope = 0.0;
              } else {
                geno_slope = kRecipDosageMax / sqrt(variance);
              }
            }
            // (ploidy * cur_allele_freq * kDosageMax) * geno_slope +
            //   geno_intercept == 0
            // bugfix: must use "-1.0 *" instead of - to avoid unsigned int
            //   wraparound
            geno_intercept = (-1.0 * kDosageMax) * ploidy_d * cur_allele_freq * geno_slope;
          }
          cur_allele_freqs[block_vidx] = cur_allele_freq;
          geno_slopes[block_vidx] = geno_slope;
          geno_intercepts[block_vidx] = geno_intercept;
        }

        *allele_end = allele_end_char;
        const uint32_t dense_vidx = block_vidx - sparse_vidx;
        double* score_coefs_iter = (difflist_common_geno == UINT32_MAX)? (&(score_dense_coefs_cmaj[dense_vidx])) : (&(score_sparse_coefs_vmaj[sparse_vidx * score_final_col_ct]));
        double* common_score_incrs_iter = common_score_incrs;
        const char* read_iter = line_start;
        for (uint32_t score_col_idx = 0; score_col_idx != score_col_ct; ++score_col_idx) {
          read_iter = NextTokenMult0(read_iter, score_col_idx_deltas[score_col_idx]);
          if (unlikely(!read_iter)) {
            goto ScoreReport_ret_MISSING_TOKENS;
          }
          double raw_coef;
          const char* token_end = ScantokDouble(read_iter, &raw_coef);
          if (unlikely(!token_end)) {
            char* errwrite_iter = strcpya_k(g_logbuf, "Error: --score");
            if (multi_input) {
              errwrite_iter = strcpya_k(errwrite_iter, "-list");
            }
            errwrite_iter = strcpya_k(errwrite_iter, ": Invalid coefficient ");
            // don't want to write this unconditionally, since it can
            // theoretically overflow the buffer.
            token_end = CurTokenEnd(read_iter);
            const uint32_t token_slen = token_end - read_iter;
            if (token_slen <= 77) {
              *errwrite_iter++ = '\'';
              errwrite_iter = memcpya(errwrite_iter, read_iter, token_slen);
              errwrite_iter = strcpya_k(errwrite_iter, "' ");
            }
            errwrite_iter = strcpya_k(errwrite_iter, "on line ");
            errwrite_iter = wtoa(line_idx, errwrite_iter);
            errwrite_iter = strcpya_k(errwrite_iter, " of ");
            errwrite_iter = strcpya(errwrite_iter, cur_input_fname);
            strcpy_k(errwrite_iter, " .\n");
            goto ScoreReport_ret_MALFORMED_INPUT_WW;
          }
          if (difflist_common_geno == UINT32_MAX) {
            // dense case: fill column for later matrix-multiply
            if (!qsr_ct) {
              *score_coefs_iter = raw_coef;
              score_coefs_iter = &(score_coefs_iter[kScoreVariantBlockSize]);
            } else {
              const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
              for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                double cur_coef = raw_coef * u31tod(IsSet(qsr_include, qsr_idx + bit_idx_base));
                *score_coefs_iter = cur_coef;
                score_coefs_iter = &(score_coefs_iter[kScoreVariantBlockSize]);
              }
            }
          } else {
            // sparse case: use variant-major representation instead
            if (!qsr_ct) {
              score_coefs_iter[0] = raw_coef;
            } else {
              const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
              for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                const double cur_coef = raw_coef * u31tod(IsSet(qsr_include, qsr_idx + bit_idx_base));
                score_coefs_iter[qsr_idx] = cur_coef;
              }
            }
            double common_dosage = 0.0;
            if (difflist_common_geno != 3) {
              common_dosage = prefer_fma(kDosageMax * geno_slope, difflist_common_geno, geno_intercept);
            } else if (!no_meanimpute) {
              common_dosage = kDosageMax * 2LL * cur_allele_freqs[block_vidx] * geno_slope;
            }
            if (!se_mode) {
              for (uint32_t uii = 0; uii != qsr_ct_nz; ++uii) {
                common_score_incrs_iter[uii] = prefer_fma(common_dosage, score_coefs_iter[uii], common_score_incrs_iter[uii]);
              }
            } else {
              for (uint32_t uii = 0; uii != qsr_ct_nz; ++uii) {
                score_coefs_iter[uii] *= score_coefs_iter[uii];
              }
              const double common_dosage_sq = common_dosage * common_dosage;
              for (uint32_t uii = 0; uii != qsr_ct_nz; ++uii) {
                common_score_incrs_iter[uii] = prefer_fma(common_dosage_sq, score_coefs_iter[uii], common_score_incrs_iter[uii]);
              }
            }
            score_coefs_iter = &(score_coefs_iter[qsr_ct_nz]);
            common_score_incrs_iter = &(common_score_incrs_iter[qsr_ct_nz]);
          }
          read_iter = token_end;
        }
        if (difflist_common_genos) {
          difflist_common_genos[block_vidx] = S_CAST(int8_t, difflist_common_geno);
          if (difflist_common_geno != UINT32_MAX) {
            difflist_lens[sparse_vidx] = difflist_len;
            raregeno_iter = &(raregeno_iter[raregeno_stride]);
            difflist_sample_ids_iter = &(difflist_sample_ids_iter[difflist_sample_ids_stride]);
            ++sparse_vidx;
          }
        }
        if (is_new_variant) {
          if (list_variants) {
            cswritep = strcpya(cswritep, variant_ids[variant_uidx]);
            AppendBinaryEoln(&cswritep);
            if (unlikely(Cswrite(&css, &cswritep))) {
              goto ScoreReport_ret_WRITE_FAIL;
            }
          }
          ++valid_variant_ct;
          if (!(valid_variant_ct % 10000)) {
            // Technically still processing a few of these variants, but this
            // is less misleading than the previous 'loaded' message.
            const uint32_t k_ct = valid_variant_ct / 1000;
            if (multi_input) {
              printf("\r--score-list file %u/%u: %uk variants processed.", file_idx1, infile_ct, k_ct);
            } else {
              printf("\r--score: %uk variants processed.", k_ct);
            }
            fflush(stdout);
          }
        }
        if (variant_uidxs) {
          variant_uidxs[block_vidx] = variant_uidx;
        }
        is_nonx_haploids[block_vidx] = is_nonx_haploid;
        is_relevant_xs[block_vidx] = is_relevant_x;
        is_ys[block_vidx] = is_y;
        ploidy_m1s[block_vidx] = ploidy_m1;
        prev_variant_uidx = variant_uidx;
        genovec_iter = &(genovec_iter[sample_ctaw2]);
        missing_bitvec_iter = &(missing_bitvec_iter[sample_ctaw]);
        missing_male_bitvec_iter = &(missing_male_bitvec_iter[sample_ctaw]);
        missing_nonfemale_bitvec_iter = &(missing_nonfemale_bitvec_iter[sample_ctaw]);
        if (dosage_present_iter) {
          dosage_present_iter = &(dosage_present_iter[sample_ctaw]);
          dosage_main_iter = &(dosage_main_iter[dosage_main_stride]);
          dosage_cts[block_vidx] = dosage_ct;
        }
        ++block_vidx;
        if (block_vidx == kScoreVariantBlockSize) {
          if (se_mode) {
            const uint32_t dense_ct = block_vidx - sparse_vidx;
            if (dense_ct) {
              for (uintptr_t score_final_col_idx = 0; score_final_col_idx != score_final_col_ct; ++score_final_col_idx) {
                double* score_dense_coefs_row = &(score_dense_coefs_cmaj[score_final_col_idx * kScoreVariantBlockSize]);
                for (uint32_t uii = 0; uii != dense_ct; ++uii) {
                  score_dense_coefs_row[uii] *= score_dense_coefs_row[uii];
                }
              }
            }
            for (uintptr_t ulii = 0; ulii != sparse_vidx * score_final_col_ct; ++ulii) {
              score_sparse_coefs_vmaj[ulii] *= score_sparse_coefs_vmaj[ulii];
            }
          }
          parity = 1 - parity;
          const uint32_t is_not_first_block = ThreadsAreActive(&tg);
          if (is_not_first_block) {
            JoinThreads(&tg);
            // CalcScoreThread() never errors out
          }
          if (unlikely(SpawnThreads(&tg))) {
            goto ScoreReport_ret_THREAD_CREATE_FAIL;
          }
          // we could instead have prev_genovec, prev_dosage_present, etc.
          // pointers, but this should be good enough
          prev_variant_uidx = UINT32_MAX;

          variant_uidxs = ctx.variant_uidxs[parity];
          genovec_iter = ctx.genovecs[parity];
          raregeno_iter = ctx.raregenos[parity];
          difflist_lens = ctx.difflist_lens[parity];
          difflist_sample_ids_iter = ctx.difflist_sample_ids[parity];
          difflist_common_genos = ctx.difflist_common_genos[parity];
          dosage_present_iter = ctx.dosage_presents[parity];
          dosage_main_iter = ctx.dosage_mains[parity];
          dosage_cts = ctx.dosage_cts[parity];
          missing_bitvec_iter = ctx.missing_bitvecs[parity];
          missing_male_bitvec_iter = ctx.missing_male_bitvecs[parity];
          missing_nonfemale_bitvec_iter = ctx.missing_nonfemale_bitvecs[parity];
          cur_allele_freqs = ctx.allele_freqs[parity];
          geno_slopes = ctx.geno_slopes[parity];
          geno_intercepts = ctx.geno_intercepts[parity];
          score_dense_coefs_cmaj = ctx.score_dense_coefs_cmaj[parity];
          score_sparse_coefs_vmaj = ctx.score_sparse_coefs_vmaj[parity];
          is_nonx_haploids = &(ctx.is_nonx_haploids[parity][0]);
          is_relevant_xs = &(ctx.is_relevant_xs[parity][0]);
          is_ys = &(ctx.is_ys[parity][0]);
          ploidy_m1s = &(ctx.ploidy_m1s[parity][0]);
          block_vidx = 0;
          sparse_vidx = 0;
        }
      }
      if (unlikely(TextStreamErrcode2(&score_txs, &reterr))) {
        goto ScoreReport_ret_TSTREAM_FAIL;
      }
      for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct_nz; ++qsr_idx) {
        VecW* missing_diploid_acc4 = &(missing_diploid_accx[acc4_vec_ct * 11 * qsr_idx]);
        VecW* missing_diploid_acc8 = &(missing_diploid_acc4[acc4_vec_ct]);
        VecW* missing_diploid_acc32 = &(missing_diploid_acc8[acc8_vec_ct]);
        VecW* missing_haploid_acc4 = &(missing_haploid_accx[acc4_vec_ct * 11 * qsr_idx]);
        VecW* missing_haploid_acc8 = &(missing_haploid_acc4[acc4_vec_ct]);
        VecW* missing_haploid_acc32 = &(missing_haploid_acc8[acc8_vec_ct]);
        VcountIncr4To8(missing_diploid_acc4, acc4_vec_ct, missing_diploid_acc8);
        VcountIncr8To32(missing_diploid_acc8, acc8_vec_ct, missing_diploid_acc32);
        VcountIncr4To8(missing_haploid_acc4, acc4_vec_ct, missing_haploid_acc8);
        VcountIncr8To32(missing_haploid_acc8, acc8_vec_ct, missing_haploid_acc32);
      }
      const uint32_t is_not_first_block = ThreadsAreActive(&tg);
      putc_unlocked('\r', stdout);
      if (missing_var_id_ct || missing_allele_code_ct || duplicated_var_id_ct) {
        missing_var_id_ct -= duplicated_var_id_ct;
        if (!missing_var_id_ct) {
          if (missing_allele_code_ct) {
            if (missing_allele_code_ct == 1) {
              snprintf(g_logbuf, kLogbufSize, "Warning: --score%s: 1 entry in %s was skipped due to a mismatching allele code.\n", multi_input? "-list" : "", cur_input_fname);
            } else {
              snprintf(g_logbuf, kLogbufSize, "Warning: --score%s: %" PRIuPTR " entries in %s were skipped due to mismatching allele codes.\n", multi_input? "-list" : "", missing_allele_code_ct, cur_input_fname);
            }
          }
        } else if (!missing_allele_code_ct) {
          if (missing_var_id_ct == 1) {
            snprintf(g_logbuf, kLogbufSize, "Warning: --score%s: 1 entry in %s was skipped due to a missing variant ID.\n", multi_input? "-list" : "", cur_input_fname);
          } else {
            snprintf(g_logbuf, kLogbufSize, "Warning: --score%s: %" PRIuPTR " entries in %s were skipped due to missing variant IDs.\n", multi_input? "-list" : "", missing_var_id_ct, cur_input_fname);
          }
        } else {
          char* warning_write_iter = strcpya_k(g_logbuf, "Warning: --score");
          if (multi_input) {
            warning_write_iter = strcpya_k(warning_write_iter, "-list");
          }
          warning_write_iter = strcpya_k(warning_write_iter, ": ");
          if (missing_var_id_ct == 1) {
            warning_write_iter = strcpya_k(warning_write_iter, "1 entry in ");
            warning_write_iter = strcpya(warning_write_iter, cur_input_fname);
            warning_write_iter = strcpya_k(warning_write_iter," was skipped due to a missing variant ID");
          } else {
            warning_write_iter = i64toa(missing_var_id_ct, warning_write_iter);
            warning_write_iter = strcpya_k(warning_write_iter, " entries in ");
            warning_write_iter = strcpya(warning_write_iter, cur_input_fname);
            warning_write_iter = strcpya_k(warning_write_iter, " were skipped due to missing variant IDs");
          }
          warning_write_iter = strcpya_k(warning_write_iter, ", and ");
          if (missing_allele_code_ct == 1) {
            warning_write_iter = strcpya_k(warning_write_iter, "1 was skipped due to a mismatching allele code");
          } else {
            warning_write_iter = i64toa(missing_allele_code_ct, warning_write_iter);
            warning_write_iter = strcpya_k(warning_write_iter, " were skipped due to mismatching allele codes");
          }
          strcpy_k(warning_write_iter, ".\n");
        }
        WordWrapB(0);
        logerrputsb();
        if (duplicated_var_id_ct) {
          if (duplicated_var_id_ct == 1) {
            snprintf(g_logbuf, kLogbufSize, "Warning: --score%s: 1 entry in %s was skipped since its variant ID appears multiple times in the main dataset.\n", multi_input? "-list" : "", cur_input_fname);
          } else {
            snprintf(g_logbuf, kLogbufSize, "Warning: --score%s: %" PRIuPTR " entries in %s were skipped since their variant IDs appear multiple times in the main dataset.\n", multi_input? "-list" : "", duplicated_var_id_ct, cur_input_fname);
          }
          WordWrapB(0);
          logerrputsb();
        }
        if (!list_variants) {
          logerrputs("(Add the 'list-variants' modifier to see which variants were actually used for\nscoring.)\n");
        }
      }
      if (block_vidx) {
        if (is_not_first_block) {
          JoinThreads(&tg);
        }
      } else if (unlikely(!valid_variant_ct)) {
        logerrprintf("Error: --score%s: No valid variants in %s.\n", multi_input? "-list" : "", cur_input_fname);
        goto ScoreReport_ret_DEGENERATE_DATA;
      } else {
        JoinThreads(&tg);
      }
      DeclareLastThreadBlock(&tg);
      ctx.cur_variant_batch_size = block_vidx;
      if (se_mode) {
        const uint32_t dense_ct = block_vidx - sparse_vidx;
        if (dense_ct) {
          for (uintptr_t score_final_col_idx = 0; score_final_col_idx != score_final_col_ct; ++score_final_col_idx) {
            double* score_dense_coefs_row = &(score_dense_coefs_cmaj[score_final_col_idx * kScoreVariantBlockSize]);
            for (uint32_t uii = 0; uii != dense_ct; ++uii) {
              score_dense_coefs_row[uii] *= score_dense_coefs_row[uii];
            }
          }
        }
        for (uintptr_t ulii = 0; ulii != sparse_vidx * score_final_col_ct; ++ulii) {
          score_sparse_coefs_vmaj[ulii] *= score_sparse_coefs_vmaj[ulii];
        }
      }
      if (unlikely(SpawnThreads(&tg))) {
        goto ScoreReport_ret_THREAD_CREATE_FAIL;
      }
      JoinThreads(&tg);
      // Sparse-optimization postprocessing.
      if (common_geno_sum_incrs) {
        uint64_t* ddosage_sums_iter = ctx.ddosage_sums;
        for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
          for (uint32_t qsr_idx = 0; qsr_idx != qsr_ct_nz; ++qsr_idx) {
            *ddosage_sums_iter += common_geno_sum_incrs[qsr_idx] * S_CAST(uint64_t, kDosageMax);
            ++ddosage_sums_iter;
          }
        }
        uintptr_t cur_sample_shard_size = sample_shard_size;
        for (uint32_t shard_idx = 0; ; ++shard_idx) {
          if (shard_idx >= sample_shard_ct_m1) {
            if (shard_idx > sample_shard_ct_m1) {
              break;
            }
            cur_sample_shard_size = sample_ct - sample_shard_ct_m1 * sample_shard_size;
          }
          double* shard_final_scores_cmaj_iter = ctx.sharded_final_scores_cmaj[shard_idx];
          for (uintptr_t ulii = 0; ulii != score_final_col_ct; ++ulii) {
            const double cur_score_incr = common_score_incrs[ulii];
            for (uint32_t shard_sample_idx = 0; shard_sample_idx != cur_sample_shard_size; ++shard_sample_idx) {
              *shard_final_scores_cmaj_iter += cur_score_incr;
              ++shard_final_scores_cmaj_iter;
            }
          }
        }
      }
      if (se_mode) {
        // sample_ct * score_final_col_ct
        uintptr_t cur_sample_shard_size = sample_shard_size;
        for (uint32_t shard_idx = 0; ; ++shard_idx) {
          if (shard_idx >= sample_shard_ct_m1) {
            if (shard_idx > sample_shard_ct_m1) {
              break;
            }
            cur_sample_shard_size = sample_ct - sample_shard_ct_m1 * sample_shard_size;
          }
          double* shard_final_scores_cmaj = ctx.sharded_final_scores_cmaj[shard_idx];
          for (uintptr_t ulii = 0; ulii != cur_sample_shard_size * score_final_col_ct; ++ulii) {
            shard_final_scores_cmaj[ulii] = sqrt(shard_final_scores_cmaj[ulii]);
          }
        }
      }
      // Now compute and spill scores to temporary file.
      for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct_nz; ++qsr_idx) {
        const uint32_t* scrambled_missing_diploid_cts = &(R_CAST(uint32_t*, missing_diploid_accx)[acc4_vec_ct * (3 + 11 * qsr_idx) * kInt32PerVec]);
        const uint32_t* scrambled_missing_haploid_cts = &(R_CAST(uint32_t*, missing_haploid_accx)[acc4_vec_ct * (3 + 11 * qsr_idx) * kInt32PerVec]);

        uint32_t cur_shard_start = 0;
        uintptr_t cur_sample_shard_size = sample_shard_size;
        uintptr_t sample_uidx_base = 0;
        uintptr_t sample_include_bits = sample_include[0];
        for (uint32_t shard_idx = 0; ; ++shard_idx) {
          if (shard_idx >= sample_shard_ct_m1) {
            if (shard_idx > sample_shard_ct_m1) {
              break;
            }
            cur_sample_shard_size = sample_ct - cur_shard_start;
          }
          const double* shard_final_scores_cmaj = ctx.sharded_final_scores_cmaj[shard_idx];
          for (uint32_t shard_sample_idx = 0; shard_sample_idx != cur_sample_shard_size; ++shard_sample_idx) {
            const uint32_t sample_idx = shard_sample_idx + cur_shard_start;
            const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
            const uint32_t scrambled_idx = VcountScramble1(sample_idx);
            uint32_t denom = allele_ct_bases[qsr_idx] + IsSet(sex_nonfemale, sample_uidx) * nonfemale_allele_ct_incrs[qsr_idx] - IsSet(sex_male, sample_uidx) * male_allele_ct_decrs[qsr_idx];
            if (no_meanimpute) {
              const uint32_t nallele = denom - 2 * scrambled_missing_diploid_cts[scrambled_idx] - scrambled_missing_haploid_cts[scrambled_idx];
              denom = nallele;
            }
            const double* final_score_col = &(shard_final_scores_cmaj[shard_sample_idx]);
            if (write_score_avgs) {
              const double denom_recip = 1.0 / S_CAST(double, denom);
              double* tmpfile_write_iter = tmpfile_buf;
              for (uintptr_t score_final_col_idx = qsr_idx; score_final_col_idx < score_final_col_ct; score_final_col_idx += qsr_ct_nz) {
                *tmpfile_write_iter++ = final_score_col[score_final_col_idx * cur_sample_shard_size] * denom_recip;
              }
              if (unlikely(fwrite_checked(tmpfile_buf, tmpfile_buf_byte_ct, score_tmpfile))) {
                goto ScoreReport_ret_WRITE_FAIL;
              }
            }
            if (write_score_sums) {
              double* tmpfile_write_iter = tmpfile_buf;
              for (uint32_t score_final_col_idx = qsr_idx; score_final_col_idx < score_final_col_ct; score_final_col_idx += qsr_ct_nz) {
                *tmpfile_write_iter++ = final_score_col[score_final_col_idx * cur_sample_shard_size];
              }
              if (unlikely(fwrite_checked(tmpfile_buf, tmpfile_buf_byte_ct, score_tmpfile))) {
                goto ScoreReport_ret_WRITE_FAIL;
              }
            }
          }
          cur_shard_start += cur_sample_shard_size;
        }
      }
      if (multi_input) {
        logprintf("--score-list file %u/%u: %u variant%s processed.\n", file_idx1, infile_ct, valid_variant_ct, (valid_variant_ct == 1)? "" : "s");
      } else {
        logprintf("--score: %u variant%s processed.\n", valid_variant_ct, (valid_variant_ct == 1)? "" : "s");
      }
    }
    if (unlikely(fclose_null(&score_tmpfile))) {
      goto ScoreReport_ret_WRITE_FAIL;
    }
    if (list_variants) {
      if (unlikely(CswriteCloseNull(&css, cswritep))) {
        goto ScoreReport_ret_WRITE_FAIL;
      }
      cswritep = nullptr;
      logprintf("--score: Variant list written to %s .\n", outname);
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".sscore.tmp");
    if (unlikely(fopen_checked(outname, FOPEN_RB, &score_tmpfile))) {
      goto ScoreReport_ret_OPEN_FAIL;
    }

    // major: file_idx
    // then qsr_idx (out of qsr_ct_nz)
    // then sample_idx (out of sample_ct)
    // then score_avgs vs. score_sums (out of write_type_ct)
    // then score_col_idx within a single input (out of score_col_ct doubles)
    const uint32_t write_type_ct = write_score_avgs + write_score_sums;
    const uint64_t sample_idx_stride = write_type_ct * score_col_ct;
    const uint64_t file_idx_stride = sample_idx_stride * sample_ct * qsr_ct_nz;
    const uint32_t* scrambled_missing_diploid_cts = nullptr;
    const uint32_t* scrambled_missing_haploid_cts = nullptr;
    double* saved_scores = nullptr;
    double* cur_saved_scores = tmpfile_buf;
    // Under --score-partitions this loop runs once and every partition lands in
    // that single file; range_names is non-null but must not reach the filename.
    const uintptr_t out_file_ct = merge_partitions? 1 : qsr_ct_nz;
    for (uintptr_t qsr_idx = 0; qsr_idx != out_file_ct; ++qsr_idx) {
      BigstackReset(bigstack_mark);
      char* outname_end2 = outname_end;
      if (range_names && (!merge_partitions)) {
        *outname_end2++ = '.';
        outname_end2 = strcpya(outname_end2, range_names[qsr_idx]);
      }
      OutnameZstSet(".sscore", output_zst, outname_end2);
      reterr = InitCstream(outname, 0, output_zst, max_thread_ct, overflow_buf_size, overflow_buf, R_CAST(unsigned char*, &(overflow_buf[overflow_buf_size])), &css);
      if (unlikely(reterr)) {
        goto ScoreReport_ret_1;
      }
      cswritep = overflow_buf;

      if (!qsr_idx) {
        const uint64_t saved_scores_byte_ct = S_CAST(uint64_t, score_final_col_ct) * infile_ct * sample_ct * write_type_ct * sizeof(double);
        if (saved_scores_byte_ct <= bigstack_left()) {
          // Guaranteed to not exhaust memory, since there are no allocations
          // after this point, and if qsr_ct > 1 the subsequent InitCstream()
          // calls will allocate the same memory as the first one.
          saved_scores = S_CAST(double*, bigstack_end_alloc_raw_rd(saved_scores_byte_ct));
          if (unlikely(fread_checked(saved_scores, saved_scores_byte_ct, score_tmpfile))) {
            goto ScoreReport_ret_READ_TMP_FAIL;
          }
        }
      }
      const uint32_t write_fid = FidColIsRequired(siip, flags / kfScoreColMaybefid);
      const char* sample_ids = siip->sample_ids;
      const char* sids = siip->sids;
      const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
      const uintptr_t max_sid_blen = siip->max_sid_blen;
      const uint32_t write_sid = SidColIsRequired(sids, flags / kfScoreColMaybesid);
      const uint32_t write_empty_pheno = (flags & kfScoreColPheno1) && (!pheno_ct);
      const uint32_t write_phenos = (flags & (kfScoreColPheno1 | kfScoreColPhenos)) && pheno_ct;
      if (write_phenos && (!(flags & kfScoreColPhenos))) {
        pheno_ct = 1;
      }
      *cswritep++ = '#';
      if (write_fid) {
        cswritep = strcpya_k(cswritep, "FID\t");
      }
      cswritep = strcpya_k(cswritep, "IID");
      if (write_sid) {
        cswritep = strcpya_k(cswritep, "\tSID");
      }
      if (write_phenos) {
        for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
          *cswritep++ = '\t';
          cswritep = strcpya(cswritep, &(pheno_names[pheno_idx * max_pheno_name_blen]));
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto ScoreReport_ret_WRITE_FAIL;
          }
        }
      } else if (write_empty_pheno) {
        cswritep = strcpya_k(cswritep, "\tPHENO1");
      }
      const uint32_t write_nallele = (flags / kfScoreColNallele) & 1;
      if (write_nallele) {
        cswritep = strcpya_k(cswritep, "\tALLELE_CT");
      }
      const uint32_t write_denom = (flags / kfScoreColDenom) & 1;
      if (write_denom) {
        cswritep = strcpya_k(cswritep, "\tDENOM");
      }
      const uint32_t write_dosage_sum = (flags / kfScoreColDosageSum) & 1;
      if (write_dosage_sum) {
        cswritep = strcpya_k(cswritep, "\tNAMED_ALLELE_DOSAGE_SUM");
      }
      // With --score-partitions each model contributes out_col_mult columns:
      // the total (range_names[0] == nullptr, so no suffix) then one per
      // partition, grouped by model so a model's columns stay adjacent.
      for (uint32_t write_type_idx = 0; write_type_idx != write_type_ct; ++write_type_idx) {
        const char* type_suffix = (write_score_avgs && (write_type_idx == 0))? "_AVG" : "_SUM";
        for (uint32_t global_score_col_idx = 0; global_score_col_idx != global_score_col_ct; ++global_score_col_idx) {
          for (uintptr_t out_qsr_idx = 0; out_qsr_idx != out_col_mult; ++out_qsr_idx) {
            *cswritep++ = '\t';
            cswritep = strcpya(cswritep, score_col_names[global_score_col_idx]);
            if (merge_partitions && range_names[out_qsr_idx]) {
              *cswritep++ = '_';
              cswritep = strcpya(cswritep, range_names[out_qsr_idx]);
            }
            cswritep = strcpya(cswritep, type_suffix);
            if (unlikely(Cswrite(&css, &cswritep))) {
              goto ScoreReport_ret_WRITE_FAIL;
            }
          }
        }
      }
      AppendBinaryEoln(&cswritep);
      if (!multi_input) {
        scrambled_missing_diploid_cts = &(R_CAST(uint32_t*, missing_diploid_accx)[acc4_vec_ct * (3 + 11 * qsr_idx) * kInt32PerVec]);
        scrambled_missing_haploid_cts = &(R_CAST(uint32_t*, missing_haploid_accx)[acc4_vec_ct * (3 + 11 * qsr_idx) * kInt32PerVec]);
      }
      const uint32_t omp_slen = strlen(output_missing_pheno);

      uintptr_t sample_uidx_base = 0;
      uintptr_t sample_include_bits = sample_include[0];
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
        cswritep = AppendXid(sample_ids, sids, write_fid, write_sid, max_sample_id_blen, max_sid_blen, sample_uidx, cswritep);
        if (write_phenos) {
          // er, this probably belongs in its own function
          for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
            const PhenoCol* cur_pheno_col = &(pheno_cols[pheno_idx]);
            const PhenoDtype type_code = cur_pheno_col->type_code;
            *cswritep++ = '\t';
            if (type_code <= kPhenoDtypeQt) {
              if (!IsSet(cur_pheno_col->nonmiss, sample_uidx)) {
                cswritep = memcpya(cswritep, output_missing_pheno, omp_slen);
              } else if (type_code == kPhenoDtypeCc) {
                *cswritep++ = '1' + IsSet(cur_pheno_col->data.cc, sample_uidx);
              } else {
                cswritep = dtoa_g(cur_pheno_col->data.qt[sample_uidx], cswritep);
              }
            } else {
              // category index guaranteed to be zero for missing values
              cswritep = strcpya(cswritep, cur_pheno_col->category_names[cur_pheno_col->data.cat[sample_uidx]]);
              if (unlikely(Cswrite(&css, &cswritep))) {
                goto ScoreReport_ret_WRITE_FAIL;
              }
            }
          }
        } else if (write_empty_pheno) {
          *cswritep++ = '\t';
          cswritep = memcpya(cswritep, output_missing_pheno, omp_slen);
        }

        if (!multi_input) {
          const uint32_t scrambled_idx = VcountScramble1(sample_idx);
          uint32_t denom = allele_ct_bases[qsr_idx] + IsSet(sex_nonfemale, sample_uidx) * nonfemale_allele_ct_incrs[qsr_idx] - IsSet(sex_male, sample_uidx) * male_allele_ct_decrs[qsr_idx];
          const uint32_t nallele = denom - 2 * scrambled_missing_diploid_cts[scrambled_idx] - scrambled_missing_haploid_cts[scrambled_idx];
          if (write_nallele) {
            *cswritep++ = '\t';
            cswritep = u32toa(nallele, cswritep);
          }
          if (no_meanimpute) {
            denom = nallele;
          }
          if (write_denom) {
            *cswritep++ = '\t';
            cswritep = u32toa(denom, cswritep);
          }
          if (write_dosage_sum) {
            *cswritep++ = '\t';
            cswritep = ddosagetoa(ctx.ddosage_sums[sample_idx + qsr_idx * sample_ct], cswritep);
          }
        }
        // Value order must match the header exactly: write_type, then model,
        // then partition.  The saved-score layout is [file][qsr][sample][type]
        // [col], so the partition loop is the strided one.
        for (uint32_t write_type_idx = 0; write_type_idx != write_type_ct; ++write_type_idx) {
          const uint64_t type_offset = (write_type_idx == 1)? score_col_ct : 0;
          for (uintptr_t score_col_idx = 0; score_col_idx != score_col_ct; ++score_col_idx) {
            for (uintptr_t out_qsr_idx = 0; out_qsr_idx != out_col_mult; ++out_qsr_idx) {
              const uintptr_t cur_qsr_idx = merge_partitions? out_qsr_idx : qsr_idx;
              const uint64_t qsr_sample_idx_offset =
                  (cur_qsr_idx * sample_ct + sample_idx) * sample_idx_stride + type_offset;
              for (uintptr_t file_idx = 0; file_idx != infile_ct; ++file_idx) {
                const uint64_t tmp_idx = file_idx * file_idx_stride + qsr_sample_idx_offset;
                if (saved_scores) {
                  cur_saved_scores = &(saved_scores[tmp_idx]);
                } else {
                  if (unlikely(fseeko(score_tmpfile, tmp_idx * sizeof(double), SEEK_SET) ||
                               fread_checked(cur_saved_scores, tmpfile_buf_byte_ct, score_tmpfile))) {
                    goto ScoreReport_ret_READ_TMP_FAIL;
                  }
                }
                *cswritep++ = '\t';
                cswritep = dtoa_g(cur_saved_scores[score_col_idx], cswritep);
              }
            }
          }
        }
        AppendBinaryEoln(&cswritep);
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto ScoreReport_ret_WRITE_FAIL;
        }
      }
      if (unlikely(CswriteCloseNull(&css, cswritep))) {
        goto ScoreReport_ret_WRITE_FAIL;
      }
    }
    if ((!qsr_ct) || merge_partitions) {
      // --score-partitions writes the one ordinary .sscore, so it reports the
      // path it actually wrote rather than a "<range name>" pattern.
      logprintfww("--score%s%s: Results written to %s .\n", multi_input? "-list" : "", merge_partitions? " + --score-partitions" : "", outname);
    } else {
      *outname_end = '\0';
      logprintfww("--score%s + --q-score-range: Results written to %s.<range name>.sscore%s .\n", multi_input? "-list" : "", outname, output_zst? ".zst" : "");
    }
    if (unlikely(fclose_null(&score_tmpfile))) {
      goto ScoreReport_ret_READ_TMP_FAIL;
    }
    snprintf(outname_end, kMaxOutfnameExtBlen, ".sscore.tmp");
    if (unlink(outname)) {
      logerrprintfww("Error: Failed to delete %s : %s.\n", outname, strerror(errno));
      goto ScoreReport_ret_WRITE_FAIL;
    }
  }
  while (0) {
  ScoreReport_ret_TSTREAM_FAIL:
    TextStreamErrPrint(cur_input_fname, &score_txs);
    break;
  ScoreReport_ret_QSR_RANGE_TSTREAM_FAIL:
    TextStreamErrPrint("--q-score-range range file", &score_txs);
    break;
  ScoreReport_ret_QSR_DATA_TSTREAM_FAIL:
    TextStreamErrPrint("--q-score-range data file", &score_txs);
    break;
  ScoreReport_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ScoreReport_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ScoreReport_ret_READ_TMP_FAIL:
    logerrprintfww(kErrprintfFread, "--score[-list] temporary file", rstrerror(errno));
    reterr = kPglRetReadFail;
    break;
  ScoreReport_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ScoreReport_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  ScoreReport_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    fputs("\n", stdout);
    logerrputsb();
  ScoreReport_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ScoreReport_ret_QSR_DATA_MISSING_TOKENS:
    logerrprintfww("Error: Line %" PRIuPTR " of --q-score-range data file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetInconsistentInput;
    break;
  ScoreReport_ret_MISSING_TOKENS:
    fputs("\n", stdout);
    logerrprintfww("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, cur_input_fname);
    reterr = kPglRetInconsistentInput;
    break;
  ScoreReport_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    fputs("\n", stdout);
    logerrputsb();
  ScoreReport_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  ScoreReport_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  ScoreReport_ret_DEGENERATE_DATA_WW:
    WordWrapB(0);
    fputs("\n", stdout);
    logerrputsb();
  ScoreReport_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  ScoreReport_ret_TOO_MANY_COLUMNS:
    logerrprintf("Error: --score%s: Too many score columns.\n", (score_info_ptr->flags & kfScoreMultiInput)? "-list" : "");
    reterr = kPglRetNotYetSupported;
    break;
  }
 ScoreReport_ret_1:
  fclose_cond(score_tmpfile);
  CswriteCloseCond(&css, cswritep);
  CleanupThreads(&tg);
  BLAS_SET_NUM_THREADS(1);
  CleanupTextStream2("--score[-list] file", &score_txs, &reterr);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}
// This setting seems optimal on my Mac (smaller doesn't take full advantage of
// AVX, larger creates cache problems?).
#ifndef NOLAPACK
// may need to increase textbuf size if this stops being true
static_assert(kMaxPc < 100000, "PhenoSvd() needs to be updated.");
#endif

#ifdef __cplusplus
}  // namespace plink2
#endif
