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

#include "adjutant_filter.h"

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <string.h>

#include "include/plink2_bits.h"
#include "include/plink2_float.h"
#include "include/plink2_htable.h"
#include "include/plink2_simd.h"
#include "include/plink2_string.h"
#include "include/plink2_text.h"
#include "include/plink2_thread.h"
#include "adjutant_decompress.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void InitExtractColCond(ExtractColCondInfo* eccip) {
  eccip->params = nullptr;
  eccip->match_flattened = nullptr;
  eccip->mismatch_flattened = nullptr;
  eccip->match_substr = 0;
  eccip->min = 0.0;
  eccip->max = DBL_MAX;
}

void CleanupExtractColCond(ExtractColCondInfo* eccip) {
  free_cond(eccip->params);
  free_cond(eccip->match_flattened);
  free_cond(eccip->mismatch_flattened);
}

PglErr FromToFlag(const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const char* varid_from, const char* varid_to, uint32_t raw_variant_ct, uint32_t max_variant_id_slen, uintptr_t variant_id_htable_size, uintptr_t* variant_include, ChrInfo* cip, uint32_t* variant_ct_ptr) {
  PglErr reterr = kPglRetSuccess;
  {
    if (!(*variant_ct_ptr)) {
      goto FromToFlag_ret_1;
    }
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t variant_uidx_start = UINT32_MAX;
    if (varid_from) {
      uint32_t cur_llidx;
      variant_uidx_start = VariantIdDupHtableFind(varid_from, variant_ids, variant_id_htable, htable_dup_base, strlen(varid_from), variant_id_htable_size, max_variant_id_slen, &cur_llidx);
      if (unlikely(variant_uidx_start == UINT32_MAX)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --from variant '%s' not found.\n", varid_from);
        goto FromToFlag_ret_INCONSISTENT_INPUT_WW;
      }
      // do *not* check variant_include here.  variant ID uniqueness should not
      // be dependent on the order in which filters are applied.
      if (unlikely(cur_llidx != UINT32_MAX)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --from variant ID '%s' appears multiple times.\n", varid_from);
        goto FromToFlag_ret_INCONSISTENT_INPUT_WW;
      }
      chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx_start);
    }
    uint32_t variant_uidx_end = 0;
    if (varid_to) {
      uint32_t cur_llidx;
      variant_uidx_end = VariantIdDupHtableFind(varid_to, variant_ids, variant_id_htable, htable_dup_base, strlen(varid_to), variant_id_htable_size, max_variant_id_slen, &cur_llidx);
      if (unlikely(variant_uidx_end == UINT32_MAX)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --to variant '%s' not found.\n", varid_to);
        goto FromToFlag_ret_INCONSISTENT_INPUT_WW;
      }
      if (unlikely(cur_llidx != UINT32_MAX)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --to variant ID '%s' appears multiple times.\n", varid_to);
        goto FromToFlag_ret_INCONSISTENT_INPUT_WW;
      }
      uint32_t chr_fo_idx2 = GetVariantChrFoIdx(cip, variant_uidx_end);
      if (variant_uidx_start == UINT32_MAX) {
        chr_fo_idx = chr_fo_idx2;
        variant_uidx_start = cip->chr_fo_vidx_start[chr_fo_idx];
      } else {
        if (unlikely(chr_fo_idx != chr_fo_idx2)) {
          logerrputs("Error: --from and --to variants are not on the same chromosome.\n");
          goto FromToFlag_ret_INCONSISTENT_INPUT;
        }
        if (variant_uidx_start > variant_uidx_end) {
          // permit order to be reversed
          uint32_t uii = variant_uidx_start;
          variant_uidx_start = variant_uidx_end;
          variant_uidx_end = uii;
        }
      }
      ++variant_uidx_end;  // convert to half-open interval
    } else {
      variant_uidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
    }
    if (variant_uidx_start) {
      ClearBitsNz(0, variant_uidx_start, variant_include);
    }
    if (variant_uidx_end < raw_variant_ct) {
      ClearBitsNz(variant_uidx_end, raw_variant_ct, variant_include);
    }
    ZeroWArr(kChrMaskWords, cip->chr_mask);
    SetBit(cip->chr_file_order[chr_fo_idx], cip->chr_mask);
    const uint32_t new_variant_ct = PopcountBitRange(variant_include, variant_uidx_start, variant_uidx_end);
    logprintf("--from/--to: %u variant%s remaining.\n", new_variant_ct, (new_variant_ct == 1)? "" : "s");
    *variant_ct_ptr = new_variant_ct;
  }
  while (0) {
  FromToFlag_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  FromToFlag_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 FromToFlag_ret_1:
  return reterr;
}
void ExtractExcludeProcessTokens(const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const char* shard_start, const char* shard_end, uint32_t variant_id_htable_size, uint32_t max_variant_id_slen, uintptr_t* already_seen) {
  const char* shard_iter = shard_start;
  while (1) {
    shard_iter = FirstPostspaceBounded(shard_iter, shard_end);
    if (shard_iter == shard_end) {
      return;
    }
    const char* token_end = CurTokenEnd(shard_iter);
    uint32_t cur_llidx;
    uint32_t variant_uidx = VariantIdDupHtableFind(shard_iter, variant_ids, variant_id_htable, htable_dup_base, token_end - shard_iter, variant_id_htable_size, max_variant_id_slen, &cur_llidx);
    shard_iter = token_end;
    if (variant_uidx == UINT32_MAX) {
      continue;
    }
    if (IsSet(already_seen, variant_uidx)) {
      continue;
    }
    for (; ; cur_llidx = htable_dup_base[cur_llidx + 1]) {
      SetBit(variant_uidx, already_seen);
      if (cur_llidx == UINT32_MAX) {
        break;
      }
      variant_uidx = htable_dup_base[cur_llidx];
    }
  }
}

CONSTI32(kMaxExtractExcludeThreads, 8);

typedef struct ExtractExcludeCtxStruct {
  const char* const* variant_ids;
  const uint32_t* variant_id_htable;
  const uint32_t* htable_dup_base;
  uintptr_t variant_id_htable_size;
  uint32_t max_variant_id_slen;

  char* shard_boundaries[kMaxExtractExcludeThreads + 1];
  uintptr_t* already_seens[kMaxExtractExcludeThreads];
} ExtractExcludeCtx;

// probable todo: rename and move into plink2_cmdline (or plink2_htable if
// VariantIdDupHtableFind is also moved).
THREAD_FUNC_DECL ExtractExcludeThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx_p1 = arg->tidx + 1;
  ExtractExcludeCtx* ctx = S_CAST(ExtractExcludeCtx*, arg->sharedp->context);

  const char* const* variant_ids = ctx->variant_ids;
  const uint32_t* variant_id_htable = ctx->variant_id_htable;
  const uint32_t* htable_dup_base = ctx->htable_dup_base;
  const uintptr_t variant_id_htable_size = ctx->variant_id_htable_size;
  const uint32_t max_variant_id_slen = ctx->max_variant_id_slen;
  uintptr_t* already_seen = ctx->already_seens[tidx_p1];
  do {
    ExtractExcludeProcessTokens(variant_ids, variant_id_htable, htable_dup_base, ctx->shard_boundaries[tidx_p1], ctx->shard_boundaries[tidx_p1 + 1], variant_id_htable_size, max_variant_id_slen, already_seen);
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr TokenExtractExclude(const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const char* fnames, const char* flagname, uint32_t raw_variant_ct, uint32_t max_variant_id_slen, uintptr_t variant_id_htable_size, VfilterType vft, uint32_t max_thread_ct, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  const char* fname_tks = nullptr;
  PglErr reterr = kPglRetSuccess;
  TokenStream tks;
  ThreadGroup tg;
  PreinitTokenStream(&tks);
  PreinitThreads(&tg);
  ExtractExcludeCtx ctx;
  {
    const uint32_t calc_thread_ct_m1 = MINV(max_thread_ct, kMaxExtractExcludeThreads) - 1;
    if (unlikely(SetThreadCt0(calc_thread_ct_m1, &tg))) {
      goto TokenExtractExclude_ret_NOMEM;
    }
    if (!(*variant_ct_ptr)) {
      goto TokenExtractExclude_ret_1;
    }
    uint32_t decompress_thread_ct = 1;
    if (max_thread_ct > calc_thread_ct_m1 + 2) {
      decompress_thread_ct = max_thread_ct - calc_thread_ct_m1 - 1;
    }
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    for (uint32_t tidx = 0; tidx <= calc_thread_ct_m1; ++tidx) {
      if (unlikely(bigstack_calloc_w(raw_variant_ctl, &(ctx.already_seens[tidx])))) {
        goto TokenExtractExclude_ret_NOMEM;
      }
    }
    if (calc_thread_ct_m1) {
      ctx.variant_ids = variant_ids;
      ctx.variant_id_htable = variant_id_htable;
      ctx.htable_dup_base = htable_dup_base;
      ctx.variant_id_htable_size = variant_id_htable_size;
      ctx.max_variant_id_slen = max_variant_id_slen;
      SetThreadFuncAndData(ExtractExcludeThread, &ctx, &tg);
    }
    const char* fnames_iter = fnames;
    do {
      if (fnames_iter == fnames) {
        fname_tks = fnames_iter;
        reterr = InitTokenStream(fnames_iter, decompress_thread_ct, &tks);
        if (unlikely(reterr)) {
          goto TokenExtractExclude_ret_TKSTREAM_FAIL;
        }
      } else {
        reterr = TokenStreamRetarget(fnames_iter, &tks);
        if (unlikely(reterr)) {
          goto TokenExtractExclude_ret_TKSTREAM_FAIL;
        }
        fname_tks = fnames_iter;
      }
      while (1) {
        reterr = TksNext(&tks, calc_thread_ct_m1 + 1, ctx.shard_boundaries);
        if (reterr) {
          break;
        }
        if (calc_thread_ct_m1) {
          if (unlikely(SpawnThreads(&tg))) {
            goto TokenExtractExclude_ret_THREAD_CREATE_FAIL;
          }
        }
        ExtractExcludeProcessTokens(variant_ids, variant_id_htable, htable_dup_base, ctx.shard_boundaries[0], ctx.shard_boundaries[1], variant_id_htable_size, max_variant_id_slen, ctx.already_seens[0]);
        JoinThreads0(&tg);
      }
      if (unlikely(reterr != kPglRetEof)) {
        goto TokenExtractExclude_ret_TKSTREAM_FAIL;
      }
      if (vft == kVfilterExtractIntersect) {
        for (uint32_t tidx = 1; tidx <= calc_thread_ct_m1; ++tidx) {
          BitvecOr(ctx.already_seens[tidx], raw_variant_ctl, ctx.already_seens[0]);
          ZeroWArr(raw_variant_ctl, ctx.already_seens[tidx]);
        }
        BitvecAnd(ctx.already_seens[0], raw_variant_ctl, variant_include);
        ZeroWArr(raw_variant_ctl, ctx.already_seens[0]);
      }
      fnames_iter = strnul(fnames_iter);
      ++fnames_iter;
    } while (*fnames_iter);
    reterr = kPglRetSuccess;
    if (vft == kVfilterExclude) {
      for (uint32_t tidx = 1; tidx <= calc_thread_ct_m1; ++tidx) {
        BitvecOr(ctx.already_seens[tidx], raw_variant_ctl, ctx.already_seens[0]);
      }
      BitvecInvmask(ctx.already_seens[0], raw_variant_ctl, variant_include);
    } else if (vft == kVfilterExtract) {
      for (uint32_t tidx = 1; tidx <= calc_thread_ct_m1; ++tidx) {
        BitvecOr(ctx.already_seens[tidx], raw_variant_ctl, ctx.already_seens[0]);
      }
      BitvecAnd(ctx.already_seens[0], raw_variant_ctl, variant_include);
    }
    const uint32_t new_variant_ct = PopcountWords(variant_include, raw_variant_ctl);
    logprintf("--%s: %u variant%s remaining.\n", flagname, new_variant_ct, (new_variant_ct == 1)? "" : "s");
    *variant_ct_ptr = new_variant_ct;
  }
  while (0) {
  TokenExtractExclude_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  TokenExtractExclude_ret_TKSTREAM_FAIL:
    TokenStreamErrPrint(fname_tks, &tks);
    break;
  TokenExtractExclude_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 TokenExtractExclude_ret_1:
  CleanupThreads(&tg);
  if (fname_tks) {
    CleanupTokenStream2(fname_tks, &tks, &reterr);
  }
  BigstackReset(bigstack_mark);
  return reterr;
}
// could permit split-chromosome here
static const char kKeepRemoveFlagStrs[4][11] = {"keep", "remove", "keep-fam", "remove-fam"};

PglErr KeepOrRemove(const char* fnames, const SampleIdInfo* siip, uint32_t raw_sample_ct, KeepFlags flags, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t orig_sample_ct = *sample_ct_ptr;
    if (!orig_sample_ct) {
      goto KeepOrRemove_ret_1;
    }
    const char* flag_name = kKeepRemoveFlagStrs[flags % 4];
    const LoadSampleIdsFlags load_sample_ids_flags = (flags & kfKeepFam)? (kfLoadSampleIdsMultifile | kfLoadSampleIdsFamOnly) : kfLoadSampleIdsMultifile;
    uintptr_t* seen_uidxs;
    uint32_t duplicate_ct;
    reterr = LoadSampleIds(fnames, sample_include, siip, flag_name, raw_sample_ct, orig_sample_ct, load_sample_ids_flags, &seen_uidxs, &duplicate_ct);
    if (unlikely(reterr)) {
      goto KeepOrRemove_ret_1;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    if (flags & kfKeepRemove) {
      BitvecInvmask(seen_uidxs, raw_sample_ctl, sample_include);
    } else {
      memcpy(sample_include, seen_uidxs, raw_sample_ctl * sizeof(intptr_t));
    }
    const uint32_t sample_ct = PopcountWords(sample_include, raw_sample_ctl);
    *sample_ct_ptr = sample_ct;
    logprintf("--%s: %u sample%s remaining.\n", flag_name, sample_ct, (sample_ct == 1)? "" : "s");
    if (duplicate_ct) {
      // "At least" since this does not count duplicate IDs absent from the
      // .fam.
      logerrprintf("Warning: At least %" PRIuPTR " duplicate ID%s in --%s file(s).\n", duplicate_ct, (duplicate_ct == 1)? "" : "s", flag_name);
    }
  }
 KeepOrRemove_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}

static_assert(kTextbufSize >= 2 * kMaxIdBlen, "KeepOneId() needs to be updated.");
void KeepOneId(const char* sample_id_flattened, const SampleIdInfo* siip, uint32_t raw_sample_ct, uint32_t iid_sid, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  // Assumes sample_id_flattened is a multistr (sequence of nonempty
  // null-terminated strings, end-of-sequence is denoted by an empty string)
  // with 1-3 parts.
  // If it has 1 part, interpret as an IID, setting FID to 0.
  // If it has 2 parts, interpret as FID-IID unless iid_sid is true.
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  {
    const char* cur_sid = nullptr;
    uintptr_t cur_sid_blen = 0;
    const char* cur_fid_iid;
    uintptr_t cur_fid_iid_blen;
    {
      const char* part1_start = sample_id_flattened;
      const uint32_t part1_blen = strlen(part1_start) + 1;
      const char* part2_start = &(part1_start[part1_blen]);
      char* fid_iid_write_iter = g_textbuf;
      if (!(*part2_start)) {
        fid_iid_write_iter = strcpya_k(fid_iid_write_iter, "0\t");
        fid_iid_write_iter = memcpya(fid_iid_write_iter, part1_start, part1_blen);
      } else {
        const uint32_t part2_blen = strlen(part2_start) + 1;
        const char* sid_start = &(part2_start[part2_blen]);
        if (!(*sid_start)) {
          if (!iid_sid) {
            fid_iid_write_iter = memcpyax(fid_iid_write_iter, part1_start, part1_blen - 1, '\t');
            fid_iid_write_iter = memcpya(fid_iid_write_iter, part2_start, part2_blen);

          } else {
            fid_iid_write_iter = strcpya_k(fid_iid_write_iter, "0\t");
            fid_iid_write_iter = memcpya(fid_iid_write_iter, part1_start, part1_blen);
            cur_sid = part2_start;
            cur_sid_blen = part2_blen;
          }
        } else {
          fid_iid_write_iter = memcpyax(fid_iid_write_iter, part1_start, part1_blen - 1, '\t');
          fid_iid_write_iter = memcpya(fid_iid_write_iter, part2_start, part2_blen);
          cur_sid = sid_start;
          cur_sid_blen = strlen(sid_start) + 1;
        }
      }
      cur_fid_iid_blen = fid_iid_write_iter - g_textbuf;
      cur_fid_iid = g_textbuf;
    }
    const char* sids = siip->sids;
    if (cur_sid) {
      if (!sids) {
        // No per-sample SID comparisons.  Either SID is always ok, or it never
        // is and we can immediately return sample_ct = 0.
        if ((siip->flags & kfSampleIdStrictSid0) && (!strequal_k_unsafe(cur_sid, "0"))) {
          goto KeepOneId_match_none;
        }
        cur_sid = nullptr;
        cur_sid_blen = 0;
      }
    } else {
      if (sids) {
        cur_sid = &(g_one_char_strs[96]);  // "0"
        cur_sid_blen = 2;
      }
    }
    const char* sample_ids = siip->sample_ids;
    const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
    const uintptr_t max_sid_blen = siip->max_sid_blen;
    if ((cur_fid_iid_blen > max_sample_id_blen) || (cur_sid_blen > max_sid_blen)) {
      goto KeepOneId_match_none;
    }
    uint32_t new_sample_ct = 0;
    for (uint32_t widx = 0; widx != raw_sample_ctl; ++widx) {
      uintptr_t cur_word = sample_include[widx];
      if (cur_word) {
        const uintptr_t sample_uidx_base = widx * kBitsPerWord;
        uintptr_t new_word = 0;
        do {
          const uint32_t sample_uidx_lowbits = ctzw(cur_word);
          const uintptr_t sample_uidx = sample_uidx_base + sample_uidx_lowbits;
          if (memequal(&(sample_ids[sample_uidx * max_sample_id_blen]), cur_fid_iid, cur_fid_iid_blen) &&
              ((!cur_sid) || memequal(&(sids[sample_uidx * max_sid_blen]), cur_sid, cur_sid_blen))) {
            new_word |= k1LU << sample_uidx_lowbits;
            ++new_sample_ct;
          }
          cur_word &= cur_word - 1;
        } while (cur_word);
        sample_include[widx] = new_word;
      }
    }
    logprintf("--indv: %u sample%s remaining.\n", new_sample_ct, (new_sample_ct == 1)? "" : "s");
    *sample_ct_ptr = new_sample_ct;
    return;
  }
 KeepOneId_match_none:
  ZeroWArr(raw_sample_ctl, sample_include);
  *sample_ct_ptr = 0;
  logputs("--indv: 0 samples remaining.\n");
}

// Minor extension of PLINK 1.x --filter.  (Renamed since --filter is not
// sufficiently self-describing; PLINK has lots of other filters on both
// samples and variants.  --filter is automatically converted to
// --keep-col-match for backward compatibility, though.)
PglErr RequirePheno(const PhenoCol* pheno_cols, const char* pheno_names, const char* require_pheno_flattened, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t is_covar, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t orig_sample_ct = *sample_ct_ptr;
    if (!orig_sample_ct) {
      goto RequirePheno_ret_1;
    }
    uint32_t required_pheno_ct = 0;
    uintptr_t max_required_pheno_blen = 2;
    uintptr_t* matched_phenos = nullptr;
    char* sorted_required_pheno_names = nullptr;
    if (require_pheno_flattened) {
      if (unlikely(MultistrToStrboxDedupAlloc(require_pheno_flattened, &sorted_required_pheno_names, &required_pheno_ct, &max_required_pheno_blen))) {
        goto RequirePheno_ret_NOMEM;
      }
      if (unlikely(bigstack_calloc_w(1 + (required_pheno_ct / kBitsPerWord), &matched_phenos))) {
        goto RequirePheno_ret_NOMEM;
      }
    } else {
      if (!pheno_ct) {
        logerrputs(is_covar? "Warning: No covariates loaded; ignoring --require-covar.\n" : "Warning: No phenotypes loaded; ignoring --require-pheno.\n");
        goto RequirePheno_ret_1;
      }
      required_pheno_ct = pheno_ct;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
      if (sorted_required_pheno_names) {
        const char* cur_pheno_name = &(pheno_names[pheno_idx * max_pheno_name_blen]);
        const int32_t ii = bsearch_strbox(cur_pheno_name, sorted_required_pheno_names, strlen(cur_pheno_name), max_required_pheno_blen, required_pheno_ct);
        if (ii == -1) {
          continue;
        }
        SetBitI(ii, matched_phenos);
      }
      BitvecAnd(pheno_cols[pheno_idx].nonmiss, raw_sample_ctl, sample_include);
    }
    if (matched_phenos) {
      const uint32_t first_unmatched_idx = AdvTo0Bit(matched_phenos, 0);
      if (unlikely(first_unmatched_idx < required_pheno_ct)) {
        logerrprintfww("Error: --require-%s '%s' not loaded.\n", is_covar? "covar covariate" : "pheno phenotype", &(sorted_required_pheno_names[first_unmatched_idx * max_required_pheno_blen]));
        goto RequirePheno_ret_INCONSISTENT_INPUT;
      }
    }
    const uint32_t new_sample_ct = PopcountWords(sample_include, raw_sample_ctl);
    const uint32_t removed_sample_ct = orig_sample_ct - new_sample_ct;
    logprintf("--require-%s: %u sample%s removed.\n", is_covar? "covar" : "pheno", removed_sample_ct, (removed_sample_ct == 1)? "" : "s");
    *sample_ct_ptr = new_sample_ct;
  }
  while (0) {
  RequirePheno_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  RequirePheno_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 RequirePheno_ret_1:
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}
void ComputeAlleleFreqs(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const uint64_t* founder_allele_ddosages, uint32_t variant_ct, double af_pseudocount, double* allele_freqs) {
  // ok for maj_alleles or allele_freqs to be nullptr
  // note that founder_allele_ddosages is in 32768ths
  // could multithread this, but not a high priority since it would only tend
  // to reduce wall-clock time by a fraction of a second
  const double af_pseudocount_ddosage = af_pseudocount * u31tod(kDosageMax);
  uint32_t cur_allele_ct = 2;
  uintptr_t variant_uidx_base = 0;
  uintptr_t cur_bits = variant_include[0];
  for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
    const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
    uintptr_t allele_idx_offset_base;
    if (!allele_idx_offsets) {
      allele_idx_offset_base = 2 * variant_uidx;
    } else {
      allele_idx_offset_base = allele_idx_offsets[variant_uidx];
      cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
    }
    const uint64_t* cur_founder_allele_ddosages = &(founder_allele_ddosages[allele_idx_offset_base]);
    uint64_t tot_ddosage = 0;
    for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct; ++allele_idx) {
      tot_ddosage += cur_founder_allele_ddosages[allele_idx];
    }
    double* cur_allele_freqs_base = &(allele_freqs[allele_idx_offset_base - variant_uidx]);
    const uint32_t cur_allele_ct_m1 = cur_allele_ct - 1;
    if (!tot_ddosage) {
      const double cur_allele_ct_m1_recip = 1.0 / u31tod(cur_allele_ct);
      for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct_m1; ++allele_idx) {
        cur_allele_freqs_base[allele_idx] = cur_allele_ct_m1_recip;
      }
    } else {
      const double adj_tot_ddosage_recip = 1.0 / (u63tod(tot_ddosage) + af_pseudocount_ddosage * u31tod(cur_allele_ct));
      for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct_m1; ++allele_idx) {
        const double cur_ddosage = u63tod(cur_founder_allele_ddosages[allele_idx]) + af_pseudocount_ddosage;
        cur_allele_freqs_base[allele_idx] = cur_ddosage * adj_tot_ddosage_recip;
      }
    }
  }
}

CONSTI32(kMaxReadFreqAlleles, 255);

// relevant column types:
// 0: variant ID
// 1: ref allele code
// 2: all alt allele codes (potentially just alt1)
//
// (3-4 are --freq only)
// 3: ref freq/count
// 4: either all freqs/counts, or all-but-ref
//
// 5: obs ct (only relevant for --freq, but can be in --geno-counts)
//
// (6-11 are --geno-counts/--freqx only)
// 6: hom-ref count
// 7: het ref-alt counts (worst case, just ref-alt1)
// 8: altx-alty counts (worst case, just hom-alt1), or all pairs
// 9: hap-ref count
// 10: hap-alt counts (worst case, just hap-alt1), or all hap counts
// 11: --geno-counts numeq (if present, ignore 6..10)
//
// overrideable:
// 12->2: ALT1
// 13->4: ALT1_FREQ/ALT1_CT
// 14->7: HET_REF_ALT1_CT
// 15->8: HOM_ALT1_CT
// 16->10: HAP_ALT1_CT
ENUM_U31_DEF_START()
  kReadFreqColVarId = 0,
  kReadFreqColRefAllele,
  kReadFreqColAltAlleles,

  kReadFreqColRefFreq,
  kReadFreqColAltFreqs,

  kReadFreqColObsCt,

  kReadFreqColHomRefCt,
  kReadFreqColHetRefAltCts,
  kReadFreqColNonrefDiploidCts,
  kReadFreqColHapRefCt,
  kReadFreqColHapAltCts,
  kReadFreqColGenoCtNumeq,

  kReadFreqColAlt1Allele,
  kReadFreqColAlt1Freq,
  kReadFreqColHetRefAlt1Ct,
  kReadFreqColHomAlt1Ct,
  kReadFreqColHapAlt1Ct,

  kReadFreqColNull
ENUM_U31_DEF_END(ReadFreqColidx);

FLAGSET_DEF_START()
  kfReadFreqColset0,
  kfReadFreqColsetVarId = (1 << kReadFreqColVarId),
  kfReadFreqColsetRefAllele = (1 << kReadFreqColRefAllele),
  kfReadFreqColsetAltAlleles = (1 << kReadFreqColAltAlleles),
  kfReadFreqColsetBase = (kfReadFreqColsetVarId | kfReadFreqColsetRefAllele | kfReadFreqColsetAltAlleles),

  kfReadFreqColsetRefFreq = (1 << kReadFreqColRefFreq),
  kfReadFreqColsetAltFreqs = (1 << kReadFreqColAltFreqs),
  kfReadFreqColsetAfreqOnly = (kfReadFreqColsetRefFreq | kfReadFreqColsetAltFreqs),

  kfReadFreqColsetObsCt = (1 << kReadFreqColObsCt),

  kfReadFreqColsetHomRefCt = (1 << kReadFreqColHomRefCt),
  kfReadFreqColsetHetRefAltCts = (1 << kReadFreqColHetRefAltCts),
  kfReadFreqColsetNonrefDiploidCts = (1 << kReadFreqColNonrefDiploidCts),
  kfReadFreqColsetHapRefCt = (1 << kReadFreqColHapRefCt),
  kfReadFreqColsetHapAltCts = (1 << kReadFreqColHapAltCts),
  kfReadFreqColsetGcountDefault = ((kfReadFreqColsetHapAltCts * 2) - kfReadFreqColsetHomRefCt),

  kfReadFreqColsetGenoCtNumeq = (1 << kReadFreqColGenoCtNumeq),
  kfReadFreqColsetGcountOnly = (kfReadFreqColsetGcountDefault | kfReadFreqColsetGenoCtNumeq),

  kfReadFreqColsetAlt1Allele = (1 << kReadFreqColAlt1Allele),
  kfReadFreqColsetAlt1Freq = (1 << kReadFreqColAlt1Freq),
  kfReadFreqColsetHetRefAlt1Ct = (1 << kReadFreqColHetRefAlt1Ct),
  kfReadFreqColsetHomAlt1Ct = (1 << kReadFreqColHomAlt1Ct),
  kfReadFreqColsetHapAlt1Ct = (1 << kReadFreqColHapAlt1Ct)
FLAGSET_DEF_END(ReadFreqColFlags);

// Support exact reconstruction of original allele frequencies in --freq counts
// case.
static inline double ForceCountToDosage(double raw_count) {
  return u63tod(S_CAST(int64_t, raw_count * kDosageMax + 0.5)) * kRecipDosageMax;
}

PglErr ReadAlleleFreqs(const uintptr_t* variant_include, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const char* read_freq_fname, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, double af_pseudocount, uint32_t max_thread_ct, double* allele_freqs, uintptr_t** variant_afreqcalcp) {
  // support PLINK 1.9 --freq/--freqx, and 2.0 --freq/--geno-counts.
  // GCTA-format no longer supported since it inhibits the allele consistency
  // check.
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream read_freq_txs;
  PreinitTextStream(&read_freq_txs);
  {
    if (!variant_ct) {
      logerrputs("Warning: Skipping --read-freq since no variants remain.\n");
      goto ReadAlleleFreqs_ret_1;
    }
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    double* cur_allele_freqs;
    uintptr_t* matched_loaded_alleles = nullptr; // spurious g++ 4.8 warning
    uintptr_t* matched_internal_alleles = nullptr;
    uint32_t* loaded_to_internal_allele_idx;
    uintptr_t* already_seen;
    if (unlikely(bigstack_calloc_w(raw_variant_ctl, variant_afreqcalcp) ||
                 bigstack_calloc_d(kMaxReadFreqAlleles, &cur_allele_freqs) ||
                 bigstack_alloc_w(BitCtToWordCt(kMaxReadFreqAlleles), &matched_loaded_alleles) ||
                 bigstack_alloc_w(BitCtToWordCt(max_allele_ct), &matched_internal_alleles) ||
                 bigstack_alloc_u32(kMaxReadFreqAlleles, &loaded_to_internal_allele_idx) ||
                 bigstack_calloc_w(raw_variant_ctl, &already_seen))) {
      goto ReadAlleleFreqs_ret_NOMEM;
    }
    bigstack_mark = R_CAST(unsigned char*, cur_allele_freqs);

    reterr = SizeAndInitTextStream(read_freq_fname, bigstack_left() / 8, MAXV(max_thread_ct - 1, 1), &read_freq_txs);
    if (unlikely(reterr)) {
      goto ReadAlleleFreqs_ret_TSTREAM_FAIL;
    }
    uint32_t* variant_id_htable = nullptr;
    uint32_t variant_id_htable_size;
    reterr = AllocAndPopulateIdHtableMt(variant_include, variant_ids, variant_ct, 0, max_thread_ct, &variant_id_htable, nullptr, &variant_id_htable_size, nullptr);
    if (unlikely(reterr)) {
      goto ReadAlleleFreqs_ret_1;
    }
    char* line_start;
    do {
      ++line_idx;
      line_start = TextGet(&read_freq_txs);
      if (unlikely(!line_start)) {
        if (!TextStreamErrcode2(&read_freq_txs, &reterr)) {
          logerrputs("Error: Empty --read-freq file.\n");
          goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
        }
        goto ReadAlleleFreqs_ret_TSTREAM_FAIL;
      }
      // automatically skip header lines that start with '##' or '# '
    } while ((*line_start == '#') && (ctou32(line_start[1]) <= '#'));
    uint32_t col_skips[kReadFreqColNull];
    ReadFreqColidx col_types[kReadFreqColNull];
    uint32_t overrideable_pos[kReadFreqColNull - kReadFreqColAlt1Allele];
    uint32_t geno_counts = 0;
    uint32_t main_eq = 0;
    uint32_t is_numeq = 0;
    uint32_t use_obs_ct = 0;
    uint32_t infer_one_freq = 0;
    uint32_t infer_freq_loaded_idx = 0;
    uint32_t relevant_col_ct = 0;

    // interpretation of ColAltAlleles
    uint32_t allele_list_just_alt1 = 1;

    uint32_t is_frac = 0;  // if true, one frequency can be missing
    // could add consistency check (can't mix FREQ and CT)

    // interpretation of ColAltFreqs and ColNonrefDiploidCts
    uint32_t main_allele_idx_start = 1;
    uint32_t main_list_just_alt1 = 1;

    // interpretation of ColHapAltCts
    uint32_t hap_allele_idx_start = 1;
    uint32_t hap_list_just_alt1 = 1;

    uint32_t het_list_just_alt1 = 1;  // ColHetRefAltCts

    uint32_t biallelic_only = 0;

    ReadFreqColFlags header_cols = kfReadFreqColset0;
    uint32_t skip_header = 1;
    if (*line_start == '#') {
      // PLINK 2.0
      // guaranteed nonspace
      const char* linebuf_iter = &(line_start[1]);

      uint32_t col_idx = 0;
      while (1) {
        const char* token_end = CurTokenEnd(linebuf_iter);
        const uint32_t token_slen = token_end - linebuf_iter;
        ReadFreqColidx cur_colidx = kReadFreqColNull;
        if (token_slen <= 4) {
          if (strequal_k(linebuf_iter, "ID", token_slen)) {
            cur_colidx = kReadFreqColVarId;
          } else if (token_slen == 3) {
            if (memequal_sk(linebuf_iter, "REF")) {
              cur_colidx = kReadFreqColRefAllele;
            } else if (memequal_sk(linebuf_iter, "ALT")) {
              cur_colidx = kReadFreqColAltAlleles;
              if (allele_list_just_alt1) {
                header_cols &= ~kfReadFreqColsetAlt1Allele;
                allele_list_just_alt1 = 0;
              }
            } else if (memequal_sk(linebuf_iter, "CTS")) {
              goto ReadAlleleFreqs_freqmain_found1;
            }
          } else if (strequal_k(linebuf_iter, "ALT1", token_slen) && allele_list_just_alt1) {
            cur_colidx = kReadFreqColAlt1Allele;
          }
        } else if (strequal_k(linebuf_iter, "REF_FREQ", token_slen) ||
                   strequal_k(linebuf_iter, "REF_CT", token_slen)) {
          cur_colidx = kReadFreqColRefFreq;
          if (linebuf_iter[4] == 'F') {
            is_frac = 1;
          }
        } else if ((strequal_k(linebuf_iter, "ALT1_FREQ", token_slen) ||
                    strequal_k(linebuf_iter, "ALT1_CT", token_slen)) &&
                   main_list_just_alt1) {
          cur_colidx = kReadFreqColAlt1Freq;
          if (linebuf_iter[5] == 'F') {
            is_frac = 1;
          }
        } else if (strequal_k(linebuf_iter, "ALT_FREQS", token_slen) ||
                   strequal_k(linebuf_iter, "ALT_CTS", token_slen)) {
          if (linebuf_iter[4] == 'F') {
            is_frac = 1;
          }
          goto ReadAlleleFreqs_freqmain_found2;
        } else if (strequal_k(linebuf_iter, "FREQS", token_slen)) {
          is_frac = 1;
          goto ReadAlleleFreqs_freqmain_found1;
        } else if (strequal_k(linebuf_iter, "ALT_NUM_FREQS", token_slen) ||
                   strequal_k(linebuf_iter, "ALT_NUM_CTS", token_slen)) {
          is_numeq = 1;
          goto ReadAlleleFreqs_freqmain_found2;
        } else if (strequal_k(linebuf_iter, "NUM_FREQS", token_slen) ||
                   strequal_k(linebuf_iter, "NUM_CTS", token_slen)) {
          is_numeq = 1;
        ReadAlleleFreqs_freqmain_found1:
          main_allele_idx_start = 0;
        ReadAlleleFreqs_freqmain_found2:
          cur_colidx = kReadFreqColAltFreqs;
          if (main_list_just_alt1) {
            header_cols &= ~kfReadFreqColsetAlt1Freq;
            main_list_just_alt1 = 0;
          }
        } else if (strequal_k(linebuf_iter, "OBS_CT", token_slen)) {
          cur_colidx = kReadFreqColObsCt;
        } else if (strequal_k(linebuf_iter, "HOM_REF_CT", token_slen)) {
          cur_colidx = kReadFreqColHomRefCt;
        } else if (strequal_k(linebuf_iter, "HET_REF_ALT1_CT", token_slen) && het_list_just_alt1) {
          cur_colidx = kReadFreqColHetRefAlt1Ct;
        } else if (strequal_k(linebuf_iter, "HET_REF_ALT_CTS", token_slen)) {
          cur_colidx = kReadFreqColHetRefAltCts;
          if (het_list_just_alt1) {
            header_cols &= ~kfReadFreqColsetHetRefAlt1Ct;
            het_list_just_alt1 = 0;
          }
        } else if (strequal_k(linebuf_iter, "HOM_ALT1_CT", token_slen) && main_list_just_alt1) {
          cur_colidx = kReadFreqColHomAlt1Ct;
        } else if (strequal_k(linebuf_iter, "TWO_ALT_GENO_CTS", token_slen) ||
                   strequal_k(linebuf_iter, "NONREF_DIPLOID_GENO_CTS", token_slen)) {
          goto ReadAlleleFreqs_countmain_found;
        } else if (strequal_k(linebuf_iter, "DIPLOID_GENO_CTS", token_slen)) {
          main_allele_idx_start = 0;
        ReadAlleleFreqs_countmain_found:
          cur_colidx = kReadFreqColNonrefDiploidCts;
          if (main_list_just_alt1) {
            header_cols &= ~kfReadFreqColsetHomAlt1Ct;
            // could make this use a different variable than FREQS does
            main_list_just_alt1 = 0;
          }
        } else if (strequal_k(linebuf_iter, "HAP_REF_CT", token_slen)) {
          cur_colidx = kReadFreqColHapRefCt;
        } else if (strequal_k(linebuf_iter, "HAP_ALT1_CT", token_slen) && hap_list_just_alt1) {
          cur_colidx = kReadFreqColHapAlt1Ct;
        } else if (strequal_k(linebuf_iter, "HAP_ALT_CTS", token_slen)) {
          goto ReadAlleleFreqs_hapmain_found;
        } else if (strequal_k(linebuf_iter, "HAP_CTS", token_slen)) {
          hap_allele_idx_start = 0;
        ReadAlleleFreqs_hapmain_found:
          cur_colidx = kReadFreqColHapAltCts;
          if (hap_list_just_alt1) {
            header_cols &= ~kfReadFreqColsetHapAlt1Ct;
            hap_list_just_alt1 = 0;
          }
        } else if (strequal_k(linebuf_iter, "GENO_NUM_CTS", token_slen)) {
          cur_colidx = kReadFreqColGenoCtNumeq;
          is_numeq = 1;
        }
        if (cur_colidx != kReadFreqColNull) {
          const ReadFreqColFlags cur_colset = S_CAST(ReadFreqColFlags, 1U << cur_colidx);
          if (unlikely(header_cols & cur_colset)) {
            logerrputs("Error: Conflicting columns in header line of --read-freq file.\n");
            goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
          }
          if (cur_colidx >= kReadFreqColAlt1Allele) {
            overrideable_pos[cur_colidx - kReadFreqColAlt1Allele] = relevant_col_ct;
          }
          header_cols |= cur_colset;
          col_skips[relevant_col_ct] = col_idx;
          col_types[relevant_col_ct++] = cur_colidx;
        }
        linebuf_iter = FirstNonTspace(token_end);
        if (IsEolnKns(*linebuf_iter)) {
          break;
        }
        ++col_idx;
      }
      ReadFreqColFlags semifinal_header_cols = header_cols;
      if (header_cols & kfReadFreqColsetAlt1Allele) {
        header_cols ^= kfReadFreqColsetAltAlleles | kfReadFreqColsetAlt1Allele;
        col_types[overrideable_pos[0]] = kReadFreqColAltAlleles;
      }
      if (header_cols & kfReadFreqColsetAlt1Freq) {
        header_cols ^= kfReadFreqColsetAltFreqs | kfReadFreqColsetAlt1Freq;
        col_types[overrideable_pos[kReadFreqColAlt1Freq - kReadFreqColAlt1Allele]] = kReadFreqColAltFreqs;
      }
      if (header_cols & kfReadFreqColsetHetRefAlt1Ct) {
        header_cols ^= kfReadFreqColsetHetRefAltCts | kfReadFreqColsetHetRefAlt1Ct;
        col_types[overrideable_pos[kReadFreqColHetRefAlt1Ct - kReadFreqColAlt1Allele]] = kReadFreqColHetRefAltCts;
      }
      if (header_cols & kfReadFreqColsetHomAlt1Ct) {
        header_cols ^= kfReadFreqColsetNonrefDiploidCts | kfReadFreqColsetHomAlt1Ct;
        col_types[overrideable_pos[kReadFreqColHomAlt1Ct - kReadFreqColAlt1Allele]] = kReadFreqColNonrefDiploidCts;
      }
      if (header_cols & kfReadFreqColsetHapAlt1Ct) {
        header_cols ^= kfReadFreqColsetHapAltCts | kfReadFreqColsetHapAlt1Ct;
        col_types[overrideable_pos[kReadFreqColHapAlt1Ct - kReadFreqColAlt1Allele]] = kReadFreqColHapAltCts;
      }
      if ((semifinal_header_cols != header_cols) && (!(header_cols & kfReadFreqColsetGenoCtNumeq))) {
        // we're treating at least one ALT1 column as if it spoke for all ALT
        // alleles
        biallelic_only = 1;
      }

      main_eq = is_numeq;
      semifinal_header_cols = header_cols;
      if (header_cols & kfReadFreqColsetAfreqOnly) {
        if (unlikely(header_cols & kfReadFreqColsetGcountOnly)) {
          logerrputs("Error: Conflicting columns in header line of --read-freq file (--freq and\n--geno-counts values mixed together).\n");
          goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
        }
        ReadFreqColFlags header_cols_exempt = kfReadFreqColset0;
        if ((header_cols & kfReadFreqColsetAltFreqs) && (!is_numeq)) {
          // [ALT_]FREQS can be formatted as either
          //   0.5,0,0.2
          // or
          //   A=0.5,G=0.2
          // Look at the first nonheader line to distinguish between these two.
          ++line_idx;
          line_start = TextGet(&read_freq_txs);
          if (unlikely(!line_start)) {
            if (!TextStreamErrcode2(&read_freq_txs, &reterr)) {
              logerrputs("Error: Empty --read-freq file.\n");
              goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
            }
            goto ReadAlleleFreqs_ret_TSTREAM_FAIL;
          }
          linebuf_iter = line_start;
          const char* alt_freq_str = nullptr;
          for (uint32_t relevant_col_idx = 0; relevant_col_idx != relevant_col_ct; ++relevant_col_idx) {
            if (col_types[relevant_col_idx] == kReadFreqColAltFreqs) {
              alt_freq_str = NextTokenMult0(linebuf_iter, col_skips[relevant_col_idx]);
              break;
            }
          }
          if (unlikely(!alt_freq_str)) {
            goto ReadAlleleFreqs_ret_MISSING_TOKENS;
          }
          const uint32_t alt_freq_slen = CurTokenEnd(alt_freq_str) - alt_freq_str;
          // bare '.' can only appear in eq formats
          main_eq = ((alt_freq_slen == 1) && (*alt_freq_str == '.')) || (memchr(alt_freq_str, '=', alt_freq_slen) != nullptr);
          if (main_eq) {
            header_cols_exempt = kfReadFreqColsetAltAlleles;
            if (!main_allele_idx_start) {
              header_cols_exempt |= kfReadFreqColsetRefAllele;
            }
            header_cols &= ~header_cols_exempt;
          }
          skip_header = 0;
        }
        if (unlikely(((header_cols & kfReadFreqColsetBase) | header_cols_exempt) != kfReadFreqColsetBase)) {
          logerrputs("Error: Missing column(s) in --read-freq file (ID, REF, ALT[1] usually\nrequired).\n");
          goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
        }
        if (!main_allele_idx_start) {
          header_cols &= ~kfReadFreqColsetRefFreq;
        } else {
          if ((header_cols & (kfReadFreqColsetRefFreq | kfReadFreqColsetAltFreqs)) != (kfReadFreqColsetRefFreq | kfReadFreqColsetAltFreqs)) {
            if (main_list_just_alt1) {
              biallelic_only = 1;
            }
            infer_one_freq = 1;
            infer_freq_loaded_idx = (header_cols / kfReadFreqColsetRefFreq) & 1;
            if (!is_frac) {
              if (unlikely(!(header_cols & kfReadFreqColsetObsCt))) {
                logerrputs("Error: Missing column(s) in --read-freq file (at least two of {REF_CT, ALT1_CT,\nALT_CTS, OBS_CT} must be present).\n");
                goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
              }
              use_obs_ct = 1;
            }
          }
        }
        logputs("--read-freq: PLINK 2 --freq file detected.\n");
      } else if (likely(header_cols & kfReadFreqColsetGcountOnly)) {
        if (unlikely((header_cols & kfReadFreqColsetBase) != kfReadFreqColsetBase)) {
          logerrputs("Error: Missing column(s) in --read-freq file (ID, REF, ALT[1] required).\n");
          goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
        }
        // possible todo: allow one frequency/count to be missing.  (not really
        // necessary since PLINK 1.9 --freqx does not leave anything out,
        // unlike PLINK 1.x --freq)
        if (header_cols & kfReadFreqColsetGenoCtNumeq) {
          // don't need anything but GENO_NUM_CTS
          header_cols &= ~kfReadFreqColsetGcountDefault;
        } else {
          // require both diploid and haploid columns for now.  (could
          // conditionally drop one of these requirements later.)
          if (unlikely(!(header_cols & kfReadFreqColsetNonrefDiploidCts))) {
            logerrputs("Error: Missing column(s) in --read-freq file (HOM_ALT1_CT,\nTWO_ALT_GENO_CTS, or DIPLOID_GENO_CTS required).\n");
            goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
          }
          if (!main_allele_idx_start) {
            header_cols &= ~(kfReadFreqColsetHomRefCt | kfReadFreqColsetHetRefAltCts);
          } else if (unlikely((header_cols & (kfReadFreqColsetHomRefCt | kfReadFreqColsetHetRefAltCts)) != (kfReadFreqColsetHomRefCt | kfReadFreqColsetHetRefAltCts))) {
            logerrputs("Error: Missing column(s) in --read-freq file (HOM_REF_CT, HET_REF_ALT1_CT, or\nHET_REF_ALT_CTS required unless {DIPLOID_}GENO_CTS present).\n");
            goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
          }
          if (unlikely(!(header_cols & kfReadFreqColsetHapAltCts))) {
            logerrputs("Error: Missing column(s) in --read-freq file (HAP_ALT1_CT, HAP_ALT_CTS, or\nHAP_CTS required).\n");
            goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
          }
          if (!hap_allele_idx_start) {
            header_cols &= ~kfReadFreqColsetHapRefCt;
          } else if (unlikely(!(header_cols & kfReadFreqColsetHapRefCt))) {
            logerrputs("Error: Missing column(s) in --read-freq file (HAP_REF_CT required unless\nHAP_CTS or GENO_CTS present).\n");
            goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
          }
        }
        geno_counts = 1;
        logputs("--read-freq: PLINK 2 --geno-counts file detected.\n");
      } else {
        logerrputs("Error: Missing column(s) in --read-freq file (no frequencies/counts).\n");
        goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
      }
      if (!use_obs_ct) {
        header_cols &= ~kfReadFreqColsetObsCt;
      }
      if (semifinal_header_cols != header_cols) {
        // remove redundant columns
        uint32_t relevant_col_idx_read = 0;
        while ((S_CAST(uint32_t, header_cols) >> col_types[relevant_col_idx_read]) & 1) {
          ++relevant_col_idx_read;
        }
        uint32_t relevant_col_idx_write = relevant_col_idx_read++;
        for (; relevant_col_idx_read != relevant_col_ct; ++relevant_col_idx_read) {
          const ReadFreqColidx cur_colidx = col_types[relevant_col_idx_read];
          if ((S_CAST(uint32_t, header_cols) >> cur_colidx) & 1) {
            col_types[relevant_col_idx_write] = cur_colidx;
            col_skips[relevant_col_idx_write] = col_skips[relevant_col_idx_read];
            ++relevant_col_idx_write;
          }
        }
        relevant_col_ct = relevant_col_idx_write;
      }
      for (uint32_t uii = relevant_col_ct - 1; uii; --uii) {
        col_skips[uii] -= col_skips[uii - 1];
      }
    } else {
      // PLINK 1.x
      // .frq:       CHR  SNP  A1  A2  MAF        NCHROBS
      // .frq.count: CHR  SNP  A1  A2  C1         C2       G0
      // .frqx:      CHR  SNP  A1  A2  C(HOM A1)  C(HET)   C(HOM A2)  C(HAP A1)
      //   C(HAP A2)  C(MISSING)
      // (yeah, the spaces in the .frqx header were a mistake, should have used
      // underscores.  oh well, live and learn.)
      col_skips[0] = 1;
      col_skips[1] = 1;
      col_skips[2] = 1;
      col_skips[3] = 1;

      col_types[0] = kReadFreqColVarId;
      // doesn't matter if we treat A1 or A2 as ref
      col_types[1] = kReadFreqColRefAllele;
      col_types[2] = kReadFreqColAltAlleles;
      biallelic_only = 1;
      if (StrStartsWithUnsafe(line_start, "CHR\tSNP\tA1\tA2\tC(HOM A1)\tC(HET)\tC(HOM A2)\tC(HAP A1)\tC(HAP A2)\tC(MISSING)")) {
        col_skips[4] = 1;
        col_skips[5] = 1;
        col_skips[6] = 1;
        col_skips[7] = 1;

        col_types[3] = kReadFreqColHomRefCt;
        col_types[4] = kReadFreqColHetRefAltCts;
        col_types[5] = kReadFreqColNonrefDiploidCts;
        col_types[6] = kReadFreqColHapRefCt;
        col_types[7] = kReadFreqColHapAltCts;
        header_cols = kfReadFreqColsetBase | kfReadFreqColsetGcountOnly;
        geno_counts = 1;
        relevant_col_ct = 8;
        logputs("--read-freq: PLINK 1.9 --freqx file detected.\n");
      } else {
        if (unlikely(!tokequal_k(line_start, "CHR"))) {
          goto ReadAlleleFreqs_ret_UNRECOGNIZED_HEADER;
        }
        const char* linebuf_iter = FirstNonTspace(&(line_start[3]));
        if (unlikely(!tokequal_k(linebuf_iter, "SNP"))) {
          goto ReadAlleleFreqs_ret_UNRECOGNIZED_HEADER;
        }
        linebuf_iter = FirstNonTspace(&(linebuf_iter[3]));
        if (unlikely(!tokequal_k(linebuf_iter, "A1"))) {
          goto ReadAlleleFreqs_ret_UNRECOGNIZED_HEADER;
        }
        linebuf_iter = FirstNonTspace(&(linebuf_iter[2]));
        if (unlikely(!tokequal_k(linebuf_iter, "A2"))) {
          goto ReadAlleleFreqs_ret_UNRECOGNIZED_HEADER;
        }
        linebuf_iter = FirstNonTspace(&(linebuf_iter[2]));
        col_types[3] = kReadFreqColRefFreq;
        if (tokequal_k(linebuf_iter, "MAF")) {
          is_frac = 1;
          infer_one_freq = 1;
          infer_freq_loaded_idx = 1;
          header_cols = kfReadFreqColsetBase | kfReadFreqColsetRefFreq;
          relevant_col_ct = 4;
          logputs("--read-freq: PLINK 1.x --freq file detected.\n");
        } else {
          if (unlikely(!tokequal_k(linebuf_iter, "C1"))) {
            goto ReadAlleleFreqs_ret_UNRECOGNIZED_HEADER;
          }
          linebuf_iter = FirstNonTspace(&(linebuf_iter[2]));
          if (unlikely(!tokequal_k(linebuf_iter, "C2"))) {
            goto ReadAlleleFreqs_ret_UNRECOGNIZED_HEADER;
          }
          col_skips[4] = 1;
          col_types[4] = kReadFreqColAltFreqs;
          header_cols = kfReadFreqColsetBase | kfReadFreqColsetAfreqOnly;
          relevant_col_ct = 5;
          logputs("--read-freq: PLINK 1.x '--freq counts' file detected.\n");
        }
      }
    }
    assert(relevant_col_ct <= 8);

    double freq_max = 4294967295.0;
    if (is_frac) {
      af_pseudocount = 0.0;
      freq_max = 1.0;
    }
    uintptr_t skipped_variant_ct = 0;
    uint32_t loaded_variant_ct = 0;
    uint32_t cur_allele_ct = 2;
    uint32_t variant_uidx = 0;
    char* line_iter = line_start;
    if (skip_header) {
      line_iter = TextLineEnd(&read_freq_txs);
      ++line_idx;
    }
    for (; TextGetUnsafe2(&read_freq_txs, &line_iter); ++line_idx) {
      {
        // not const since tokens may be null-terminated or comma-terminated
        // later
        char* token_ptrs[12];
        uint32_t token_slens[12];
        line_iter = TokenLex0(line_iter, R_CAST(uint32_t*, col_types), col_skips, relevant_col_ct, token_ptrs, token_slens);
        if (unlikely(!line_iter)) {
          goto ReadAlleleFreqs_ret_MISSING_TOKENS;
        }
        // characters may be modified, so better find the \n now just in case
        // it gets clobbered
        line_iter = AdvPastDelim(line_iter, '\n');
        const char* variant_id_start = token_ptrs[kReadFreqColVarId];
        const uint32_t variant_id_slen = token_slens[kReadFreqColVarId];
        variant_uidx = VariantIdDupflagHtableFind(variant_id_start, variant_ids, variant_id_htable, variant_id_slen, variant_id_htable_size, max_variant_id_slen);
        if (variant_uidx >> 31) {
          if (likely(variant_uidx == UINT32_MAX)) {
            ++skipped_variant_ct;
            continue;
          }
          snprintf(g_logbuf, kLogbufSize, "Error: --read-freq variant ID '%s' appears multiple times in main dataset.\n", variant_ids[variant_uidx & 0x7fffffff]);
          goto ReadAlleleFreqs_ret_MALFORMED_INPUT_WW;
        }
        if (unlikely(IsSet(already_seen, variant_uidx))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Variant ID '%s' appears multiple times in --read-freq file.\n", variant_ids[variant_uidx]);
          goto ReadAlleleFreqs_ret_MALFORMED_INPUT_WW;
        }
        SetBit(variant_uidx, already_seen);

        uintptr_t allele_idx_offset_base;
        if (!allele_idx_offsets) {
          allele_idx_offset_base = variant_uidx * 2;
        } else {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
          if (biallelic_only) {
            goto ReadAlleleFreqs_skip_variant;
          }
        }
        ZeroWArr(BitCtToWordCt(cur_allele_ct), matched_internal_alleles);
        const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
        uint32_t loaded_allele_idx_end = 0;
        {
          uint32_t unmatched_allele_ct = cur_allele_ct;
          uint32_t cur_loaded_allele_code_slen;
          char* cur_loaded_allele_code;
          char* loaded_allele_code_iter;
          char* loaded_allele_code_end;
          if (header_cols & kfReadFreqColsetRefAllele) {
            cur_loaded_allele_code_slen = token_slens[kReadFreqColRefAllele];
            cur_loaded_allele_code = token_ptrs[kReadFreqColRefAllele];
            // No special handling of missing allele code needed (if the
            // fileset is valid).
            cur_loaded_allele_code[cur_loaded_allele_code_slen] = '\0';
            if (header_cols & kfReadFreqColsetAltAlleles) {
              loaded_allele_code_iter = token_ptrs[kReadFreqColAltAlleles];
              loaded_allele_code_end = &(loaded_allele_code_iter[token_slens[kReadFreqColAltAlleles]]);
              *loaded_allele_code_end++ = ',';
            } else {
              // special case: with --freq alteq or alteqz column, we only need
              // to scrape REF here
              loaded_allele_code_iter = &(cur_loaded_allele_code[cur_loaded_allele_code_slen + 1]);
              loaded_allele_code_end = loaded_allele_code_iter;
            }
          } else {
            // no REF, only ALT.  skip the first iteration of the loop.
            loaded_allele_code_iter = token_ptrs[kReadFreqColAltAlleles];
            loaded_allele_code_end = &(loaded_allele_code_iter[token_slens[kReadFreqColAltAlleles]]);
            *loaded_allele_code_end++ = ',';

            cur_loaded_allele_code = loaded_allele_code_iter;
            loaded_allele_code_iter = S_CAST(char*, memchr(loaded_allele_code_iter, ',', loaded_allele_code_end - cur_loaded_allele_code));
            cur_loaded_allele_code_slen = loaded_allele_code_iter - cur_loaded_allele_code;
            *loaded_allele_code_iter++ = '\0';
            loaded_allele_idx_end = 1;
            matched_loaded_alleles[0] = 0;
          }
          uint32_t widx = 0;
          while (1) {
            if (!(loaded_allele_idx_end % kBitsPerWord)) {
              widx = loaded_allele_idx_end / kBitsPerWord;
              matched_loaded_alleles[widx] = 0;
            }
            if (cur_loaded_allele_code_slen <= max_allele_slen) {
              const uint32_t cur_blen = cur_loaded_allele_code_slen + 1;
              uintptr_t internal_allele_idx_base = 0;
              uintptr_t cur_inv_bits = ~matched_internal_alleles[0];
              for (uint32_t unmatched_allele_idx = 0; unmatched_allele_idx != unmatched_allele_ct; ++unmatched_allele_idx) {
                const uintptr_t internal_allele_idx = BitIter0(matched_internal_alleles, &internal_allele_idx_base, &cur_inv_bits);
                if (memequal(cur_loaded_allele_code, cur_alleles[internal_allele_idx], cur_blen)) {
                  if (unlikely(IsSet(matched_internal_alleles, internal_allele_idx))) {
                    snprintf(g_logbuf, kLogbufSize, "Error: Duplicate allele code on line %" PRIuPTR " of --read-freq file.\n", line_idx);
                    goto ReadAlleleFreqs_ret_MALFORMED_INPUT_2;
                  }
                  SetBit(internal_allele_idx, matched_internal_alleles);
                  SetBit(loaded_allele_idx_end, matched_loaded_alleles);
                  loaded_to_internal_allele_idx[loaded_allele_idx_end] = internal_allele_idx;
                  break;
                }
              }
            }
            ++loaded_allele_idx_end;
            if (loaded_allele_code_iter == loaded_allele_code_end) {
              break;
            }
            if (unlikely(loaded_allele_idx_end == kMaxReadFreqAlleles)) {
              snprintf(g_logbuf, kLogbufSize, "Error: --read-freq file entry for variant ID '%s' has more than %u ALT alleles.\n", variant_ids[variant_uidx], kMaxReadFreqAlleles - 1);
              goto ReadAlleleFreqs_ret_MALFORMED_INPUT_WW;
            }
            cur_loaded_allele_code = loaded_allele_code_iter;
            loaded_allele_code_iter = S_CAST(char*, memchr(loaded_allele_code_iter, ',', loaded_allele_code_end - cur_loaded_allele_code));
            cur_loaded_allele_code_slen = loaded_allele_code_iter - cur_loaded_allele_code;
            *loaded_allele_code_iter++ = '\0';
          }
        }
        // bugfix (31 May 2023): forgot to skip variant when not enough alleles
        // were loaded.
        {
          const uint32_t matched_allele_ct = PopcountWords(matched_loaded_alleles, BitCtToWordCt(loaded_allele_idx_end));
          if (matched_allele_ct != cur_allele_ct) {
            assert(matched_allele_ct < cur_allele_ct);
            goto ReadAlleleFreqs_skip_variant;
          }
        }

        double* allele_freqs_write = &(allele_freqs[allele_idx_offset_base - variant_uidx]);
        if (geno_counts) {
          ZeroDArr(cur_allele_ct, cur_allele_freqs);
          if (is_numeq) {
            const uint32_t full_slen = token_slens[kReadFreqColGenoCtNumeq];
            char* geno_num_cts = token_ptrs[kReadFreqColGenoCtNumeq];
            if (full_slen > 1) {
              geno_num_cts[full_slen] = ',';
              const char* geno_num_cts_iter = geno_num_cts;
              const char* geno_num_cts_end = &(geno_num_cts[full_slen]);
#ifndef __LP64__
              const uint32_t cap_div_10 = (loaded_allele_idx_end - 1) / 10;
              const uint32_t cap_mod_10 = (loaded_allele_idx_end - 1) % 10;
#endif
              while (1) {
                uint32_t second_loaded_allele_idx = UINT32_MAX;
                uint32_t first_loaded_allele_idx;
#ifdef __LP64__
                if (unlikely(ScanmovUintCapped(loaded_allele_idx_end - 1, &geno_num_cts_iter, &first_loaded_allele_idx))) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                if (*geno_num_cts_iter == '/') {
                  ++geno_num_cts_iter;
                  if (unlikely(ScanmovUintCapped(loaded_allele_idx_end - 1, &geno_num_cts_iter, &second_loaded_allele_idx))) {
                    goto ReadAlleleFreqs_ret_INVALID_FREQS;
                  }
                }
#else
                if (unlikely(ScanmovUintCapped32(cap_div_10, cap_mod_10, &geno_num_cts_iter, &first_loaded_allele_idx))) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                if (*geno_num_cts_iter == '/') {
                  ++geno_num_cts_iter;
                  if (unlikely(ScanmovUintCapped32(cap_div_10, cap_mod_10, &geno_num_cts_iter, &second_loaded_allele_idx))) {
                    goto ReadAlleleFreqs_ret_INVALID_FREQS;
                  }
                }
#endif
                if (unlikely(*geno_num_cts_iter != '=')) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                ++geno_num_cts_iter;
                double dxx;
                const char* cur_ct_end = ScanadvDouble(geno_num_cts_iter, &dxx);
                if (unlikely((!cur_ct_end) || (*cur_ct_end != ',') || (dxx < 0.0) || (dxx > 4294967295.0))) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                dxx = ForceCountToDosage(dxx);
                if (IsSet(matched_loaded_alleles, first_loaded_allele_idx)) {
                  cur_allele_freqs[loaded_to_internal_allele_idx[first_loaded_allele_idx]] += dxx;
                }
                if ((second_loaded_allele_idx != UINT32_MAX) && IsSet(matched_loaded_alleles, second_loaded_allele_idx)) {
                  cur_allele_freqs[loaded_to_internal_allele_idx[second_loaded_allele_idx]] += dxx;
                }
                geno_num_cts_iter = cur_ct_end;
                if (geno_num_cts_iter == geno_num_cts_end) {
                  break;
                }
                ++geno_num_cts_iter;
              }
            } else if (unlikely(*geno_num_cts != '.')) {
              goto ReadAlleleFreqs_ret_INVALID_FREQS;
            }
          } else {
            const uint32_t internal0 = IsSet(matched_loaded_alleles, 0)? loaded_to_internal_allele_idx[0] : UINT32_MAX;
            if (header_cols & kfReadFreqColsetHomRefCt) {
              if (internal0 != UINT32_MAX) {
                const char* hom_ref_str = token_ptrs[kReadFreqColHomRefCt];
                double dxx;
                const char* hom_ref_end = ScantokDouble(hom_ref_str, &dxx);
                if (unlikely((!hom_ref_end) || (dxx < 0.0) || (dxx > 4294967295.0))) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                dxx = ForceCountToDosage(dxx);
                cur_allele_freqs[internal0] += 2 * dxx;
              }

              char* het_refalt = token_ptrs[kReadFreqColHetRefAltCts];
              const uint32_t het_refalt_slen = token_slens[kReadFreqColHetRefAltCts];
              het_refalt[het_refalt_slen] = ',';
              const char* het_refalt_iter = het_refalt;
              const char* het_refalt_end = &(het_refalt[het_refalt_slen]);
              for (uint32_t alt_allele_idx = 1; alt_allele_idx != cur_allele_ct; ++alt_allele_idx) {
                if (unlikely(het_refalt_iter >= het_refalt_end)) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                double dxx;
                const char* cur_entry_end = ScanadvDouble(het_refalt_iter, &dxx);
                if (unlikely((!cur_entry_end) || (*cur_entry_end != ',') || (dxx < 0.0) || (dxx > 4294967295.0))) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                dxx = ForceCountToDosage(dxx);
                if (internal0 != UINT32_MAX) {
                  cur_allele_freqs[internal0] += dxx;
                }
                if (IsSet(matched_loaded_alleles, alt_allele_idx)) {
                  cur_allele_freqs[loaded_to_internal_allele_idx[alt_allele_idx]] += dxx;
                }
                het_refalt_iter = &(cur_entry_end[1]);
              }
            }
            // ColNonrefDiploidCts required
            char* diploid_cts = token_ptrs[kReadFreqColNonrefDiploidCts];
            const uint32_t diploid_cts_slen = token_slens[kReadFreqColNonrefDiploidCts];
            diploid_cts[diploid_cts_slen] = ',';
            const char* diploid_cts_iter = diploid_cts;
            const char* diploid_cts_end = &(diploid_cts[diploid_cts_slen]);
            for (uint32_t second_allele_idx = main_allele_idx_start; second_allele_idx != cur_allele_ct; ++second_allele_idx) {
              uint32_t internalx = UINT32_MAX;
              if (IsSet(matched_loaded_alleles, second_allele_idx)) {
                internalx = loaded_to_internal_allele_idx[second_allele_idx];
              }
              // 1/1, 1/2, 2/2, 1/3, ...
              for (uint32_t first_allele_idx = main_allele_idx_start; first_allele_idx <= second_allele_idx; ++first_allele_idx) {
                if (unlikely(diploid_cts_iter >= diploid_cts_end)) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                double dxx;
                const char* cur_entry_end = ScanadvDouble(diploid_cts_iter, &dxx);
                if (unlikely((!cur_entry_end) || (*cur_entry_end != ',') || (dxx < 0.0) || (dxx > 4294967295.0))) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                dxx = ForceCountToDosage(dxx);
                if (IsSet(matched_loaded_alleles, first_allele_idx)) {
                  cur_allele_freqs[loaded_to_internal_allele_idx[first_allele_idx]] += dxx;
                }
                if (internalx != UINT32_MAX) {
                  cur_allele_freqs[internalx] += dxx;
                }
                diploid_cts_iter = &(cur_entry_end[1]);
              }
            }

            if ((header_cols & kfReadFreqColsetHapRefCt) && (internal0 != UINT32_MAX)) {
              const char* hap_ref_str = token_ptrs[kReadFreqColHapRefCt];
              double dxx;
              const char* hap_ref_end = ScantokDouble(hap_ref_str, &dxx);
              if (unlikely((!hap_ref_end) || (dxx < 0.0) || (dxx > 4294967295.0))) {
                goto ReadAlleleFreqs_ret_INVALID_FREQS;
              }
              dxx = ForceCountToDosage(dxx);
              cur_allele_freqs[internal0] += dxx;
            }
            // ColHapAltCts required
            char* hap_alt = token_ptrs[kReadFreqColHapAltCts];
            const uint32_t hap_alt_slen = token_slens[kReadFreqColHapAltCts];
            hap_alt[hap_alt_slen] = ',';
            const char* hap_alt_iter = hap_alt;
            const char* hap_alt_end = &(hap_alt[hap_alt_slen]);
            for (uint32_t alt_allele_idx = 1; alt_allele_idx != cur_allele_ct; ++alt_allele_idx) {
              if (unlikely(hap_alt_iter >= hap_alt_end)) {
                goto ReadAlleleFreqs_ret_INVALID_FREQS;
              }
              double dxx;
              const char* cur_entry_end = ScanadvDouble(hap_alt_iter, &dxx);
              if (unlikely((!cur_entry_end) || (*cur_entry_end != ',') || (dxx < 0.0) || (dxx > 4294967295.0))) {
                goto ReadAlleleFreqs_ret_INVALID_FREQS;
              }
              if (IsSet(matched_loaded_alleles, alt_allele_idx)) {
                dxx = ForceCountToDosage(dxx);
                cur_allele_freqs[loaded_to_internal_allele_idx[alt_allele_idx]] += dxx;
              }
              hap_alt_iter = &(cur_entry_end[1]);
            }
          }
        } else {
          if ((header_cols & kfReadFreqColsetRefFreq) && IsSet(matched_loaded_alleles, 0)) {
            const char* ref_freq_str = token_ptrs[kReadFreqColRefFreq];
            double dxx;
            if (!ScantokDouble(ref_freq_str, &dxx)) {
              if (likely(IsNanStr(ref_freq_str, token_slens[kReadFreqColRefFreq]))) {
                goto ReadAlleleFreqs_skip_variant;
              }
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid REF frequency/count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
              goto ReadAlleleFreqs_ret_MALFORMED_INPUT_WW;
            }
            if (unlikely((dxx < 0.0) || (dxx > freq_max))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid REF frequency/count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
              goto ReadAlleleFreqs_ret_MALFORMED_INPUT_WW;
            }
            if (!is_frac) {
              dxx = ForceCountToDosage(dxx);
            }
            cur_allele_freqs[loaded_to_internal_allele_idx[0]] = dxx;
          }
          if (header_cols & kfReadFreqColsetAltFreqs) {
            const uint32_t full_slen = token_slens[kReadFreqColAltFreqs];
            char* alt_freq_start = token_ptrs[kReadFreqColAltFreqs];
            alt_freq_start[full_slen] = ',';
            if (!main_eq) {
              const char* alt_freq_iter = alt_freq_start;
              const char* alt_freq_end = &(alt_freq_start[full_slen]);
              for (uint32_t allele_idx = main_allele_idx_start; allele_idx != loaded_allele_idx_end; ++allele_idx, ++alt_freq_iter) {
                if (unlikely(alt_freq_iter >= alt_freq_end)) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                if (!IsSet(matched_loaded_alleles, allele_idx)) {
                  alt_freq_iter = AdvToDelim(alt_freq_iter, ',');
                  continue;
                }
                double dxx;
                const char* cur_freq_end = ScanadvDouble(alt_freq_iter, &dxx);
                if (!cur_freq_end) {
                  cur_freq_end = AdvToDelim(alt_freq_iter, ',');
                  if (likely(IsNanStr(alt_freq_iter, cur_freq_end - alt_freq_iter))) {
                    goto ReadAlleleFreqs_skip_variant;
                  }
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                if (unlikely((*cur_freq_end != ',') || (dxx < 0.0) || (dxx > freq_max))) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                if (!is_frac) {
                  dxx = ForceCountToDosage(dxx);
                }
                alt_freq_iter = cur_freq_end;
                cur_allele_freqs[loaded_to_internal_allele_idx[allele_idx]] = dxx;
              }
            } else {
              ZeroDArr(cur_allele_ct, cur_allele_freqs);
              if ((full_slen > 1) || (*alt_freq_start != '.')) {
                if (is_numeq) {
                  const char* alt_freq_iter = alt_freq_start;
                  const char* alt_freq_end = &(alt_freq_start[full_slen]);
#ifndef __LP64__
                  const uint32_t cap_div_10 = (loaded_allele_idx_end - 1) / 10;
                  const uint32_t cap_mod_10 = (loaded_allele_idx_end - 1) % 10;
#endif
                  while (1) {
                    const char* cur_entry_end = AdvToDelim(alt_freq_iter, ',');
                    uint32_t loaded_allele_idx;
#ifdef __LP64__
                    if (unlikely(ScanmovUintCapped(loaded_allele_idx_end - 1, &alt_freq_iter, &loaded_allele_idx))) {
                      goto ReadAlleleFreqs_ret_INVALID_FREQS;
                    }
#else
                    if (unlikely(ScanmovUintCapped32(cap_div_10, cap_mod_10, &alt_freq_iter, &loaded_allele_idx))) {
                      goto ReadAlleleFreqs_ret_INVALID_FREQS;
                    }
#endif
                    if (unlikely(*alt_freq_iter != '=')) {
                      goto ReadAlleleFreqs_ret_INVALID_FREQS;
                    }
                    if (IsSet(matched_loaded_alleles, loaded_allele_idx)) {
                      const uint32_t internal_allele_idx = loaded_to_internal_allele_idx[loaded_allele_idx];
                      if (unlikely(cur_allele_freqs[internal_allele_idx] != 0.0)) {
                        snprintf(g_logbuf, kLogbufSize, "Error: Duplicate entry on line %" PRIuPTR " of --read-freq file.\n", line_idx);
                        goto ReadAlleleFreqs_ret_MALFORMED_INPUT_2;
                      }
                      ++alt_freq_iter;
                      double dxx;
                      const char* cur_freq_end = ScantokDouble(alt_freq_iter, &dxx);
                      if (!cur_freq_end) {
                        if (likely(IsNanStr(alt_freq_iter, cur_entry_end - alt_freq_iter))) {
                          goto ReadAlleleFreqs_skip_variant;
                        }
                        goto ReadAlleleFreqs_ret_INVALID_FREQS;
                      }
                      if (unlikely((dxx < 0.0) || (dxx > freq_max))) {
                        goto ReadAlleleFreqs_ret_INVALID_FREQS;
                      }
                      if (!is_frac) {
                        dxx = ForceCountToDosage(dxx);
                      }
                      cur_allele_freqs[internal_allele_idx] = dxx;
                    }
                    alt_freq_iter = cur_entry_end;
                    if (alt_freq_iter == alt_freq_end) {
                      break;
                    }
                    ++alt_freq_iter;
                  }
                } else {
                  char* alt_freq_iter = alt_freq_start;
                  char* alt_freq_end = &(alt_freq_start[full_slen]);
                  while (1) {
                    char* cur_entry_end = AdvToDelim(alt_freq_iter, ',');
                    const uint32_t cur_entry_slen = cur_entry_end - alt_freq_iter;
                    char* eq_ptr = S_CAST(char*, memchr(alt_freq_iter, '=', cur_entry_slen));
                    if (unlikely(!eq_ptr)) {
                      goto ReadAlleleFreqs_ret_INVALID_FREQS;
                    }

                    // necessary for string comparison
                    *eq_ptr++ = '\0';

                    const uint32_t cur_blen = eq_ptr - alt_freq_iter;
                    // O(n^2), may want to replace with O(n log n)
                    for (uint32_t internal_allele_idx = 0; internal_allele_idx != cur_allele_ct; ++internal_allele_idx) {
                      if (memequal(alt_freq_iter, cur_alleles[internal_allele_idx], cur_blen)) {
                        if (unlikely(cur_allele_freqs[internal_allele_idx] != 0.0)) {
                          snprintf(g_logbuf, kLogbufSize, "Error: Duplicate entry on line %" PRIuPTR " of --read-freq file.\n", line_idx);
                          goto ReadAlleleFreqs_ret_MALFORMED_INPUT_2;
                        }
                        alt_freq_iter = eq_ptr;
                        double dxx;
                        const char* cur_freq_end = ScantokDouble(alt_freq_iter, &dxx);
                        if (!cur_freq_end) {
                          if (likely(IsNanStr(alt_freq_iter, cur_entry_end - alt_freq_iter))) {
                            goto ReadAlleleFreqs_skip_variant;
                          }
                          goto ReadAlleleFreqs_ret_INVALID_FREQS;
                        }
                        if (unlikely((dxx < 0.0) || (dxx > freq_max))) {
                          goto ReadAlleleFreqs_ret_INVALID_FREQS;
                        }
                        if (!is_frac) {
                          dxx = ForceCountToDosage(dxx);
                        }
                        cur_allele_freqs[internal_allele_idx] = dxx;
                        break;
                      }
                    }
                    alt_freq_iter = cur_entry_end;
                    if (alt_freq_iter == alt_freq_end) {
                      break;
                    }
                    ++alt_freq_iter;
                  }
                }
              }
            }
          }
        }
        // The IsSet condition is false when the following happens:
        // - --read-freq file has multiple ALT alleles.
        // - These ALT alleles collectively match *all* alleles in the main
        //   dataset's version of the variant, including the main dataset's
        //   REF.  The --read-freq file's REF allele is not in the main
        //   dataset at all.
        // - So, while we normally infer the REF frequency by subtracting the
        //   sum of the other frequencies from 1, in this subcase there is
        //   nothing to infer.
        if (infer_one_freq && IsSet(matched_loaded_alleles, infer_freq_loaded_idx)) {
          double adj_obs_ct_recip = 1.0;
          if (header_cols & kfReadFreqColsetObsCt) {
            uint32_t obs_ct_raw;
            if (unlikely(ScanUintCapped(token_ptrs[kReadFreqColObsCt], UINT32_MAX, &obs_ct_raw))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid allele count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
              goto ReadAlleleFreqs_ret_MALFORMED_INPUT_2;
            }
            if ((!obs_ct_raw) && (af_pseudocount == 0.0)) {
              goto ReadAlleleFreqs_skip_variant;
            }
            adj_obs_ct_recip = 1.0 / (u63tod(obs_ct_raw) + af_pseudocount * cur_allele_ct);
          }
          const uint32_t infer_freq_internal_idx = loaded_to_internal_allele_idx[infer_freq_loaded_idx];
          if (cur_allele_ct == 2) {
            // optimize common case
            double known_freq_d = cur_allele_freqs[1 - infer_freq_internal_idx] + af_pseudocount;
            double known_scaled_freq = known_freq_d * adj_obs_ct_recip;
            if (known_scaled_freq <= 1.0) {
              if (infer_freq_internal_idx) {
                allele_freqs_write[0] = known_scaled_freq;
              } else {
                allele_freqs_write[0] = 1.0 - known_scaled_freq;
              }
            } else if (likely(known_scaled_freq < ((1.0 + kSmallEpsilon) / 0.99))) {
              if (infer_freq_internal_idx) {
                allele_freqs_write[0] = 1.0;
              } else {
                allele_freqs_write[0] = 0.0;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Frequency/count too large on line %" PRIuPTR " of --read-freq file.\n", line_idx);
              goto ReadAlleleFreqs_ret_MALFORMED_INPUT_2;
            }
          } else {
            if (af_pseudocount != 0.0) {
              for (uint32_t internal_allele_idx = 0; internal_allele_idx != cur_allele_ct; ++internal_allele_idx) {
                cur_allele_freqs[internal_allele_idx] += af_pseudocount;
              }
            }
            cur_allele_freqs[infer_freq_internal_idx] = 0.0;
            double known_freq_sum_d = 0.0;
            for (uint32_t internal_allele_idx = 0; internal_allele_idx != cur_allele_ct; ++internal_allele_idx) {
              known_freq_sum_d += cur_allele_freqs[internal_allele_idx];
            }
            double known_scaled_freq_sum = known_freq_sum_d * adj_obs_ct_recip;
            if (likely(known_scaled_freq_sum < ((1.0 + kSmallEpsilon) / 0.99))) {
              if (known_scaled_freq_sum > 1.0) {
                // possible rounding error, rescale
                adj_obs_ct_recip = 1.0 / known_scaled_freq_sum;
                known_scaled_freq_sum = 1.0;
              }
              const uint32_t cur_allele_ct_m1 = cur_allele_ct - 1;
              for (uint32_t internal_allele_idx = 0; internal_allele_idx != cur_allele_ct_m1; ++internal_allele_idx) {
                double dxx;
                if (internal_allele_idx == infer_freq_internal_idx) {
                  dxx = 1.0 - known_scaled_freq_sum;
                } else {
                  dxx = adj_obs_ct_recip * cur_allele_freqs[internal_allele_idx];
                }
                allele_freqs_write[internal_allele_idx] = dxx;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Frequency/count too large on line %" PRIuPTR " of --read-freq file.\n", line_idx);
              goto ReadAlleleFreqs_ret_MALFORMED_INPUT_2;
            }
          }
        } else {
          // complete frequency or count data
          if (af_pseudocount != 0.0) {
            for (uint32_t internal_allele_idx = 0; internal_allele_idx != cur_allele_ct; ++internal_allele_idx) {
              cur_allele_freqs[internal_allele_idx] += af_pseudocount;
            }
          }
          double tot_freq = 0.0;
          for (uint32_t internal_allele_idx = 0; internal_allele_idx != cur_allele_ct; ++internal_allele_idx) {
            tot_freq += cur_allele_freqs[internal_allele_idx];
          }
          if (tot_freq == 0.0) {
            goto ReadAlleleFreqs_skip_variant;
          }
          const double tot_freq_recip = 1.0 / tot_freq;
          const uint32_t cur_allele_ct_m1 = cur_allele_ct - 1;
          for (uint32_t internal_allele_idx = 0; internal_allele_idx != cur_allele_ct_m1; ++internal_allele_idx) {
            allele_freqs_write[internal_allele_idx] = tot_freq_recip * cur_allele_freqs[internal_allele_idx];
          }
        }

        ++loaded_variant_ct;
        if (!(loaded_variant_ct % 10000)) {
          printf("\r--read-freq: Frequencies for %uk variants loaded.", loaded_variant_ct / 1000);
          fflush(stdout);
        }
      }
      while (0) {
      ReadAlleleFreqs_skip_variant:
        SetBit(variant_uidx, *variant_afreqcalcp);
        ++skipped_variant_ct;
      }
    }
    if (unlikely(TextStreamErrcode2(&read_freq_txs, &reterr))) {
      goto ReadAlleleFreqs_ret_TSTREAM_FAIL;
    }
    // Set variant_afreqcalc to the set of variants in variant_include, but
    // without loaded frequencies.
    BitvecInvertAndMask(variant_include, raw_variant_ctl, already_seen);
    // already_seen is now the set of remaining variants that *weren't* seen,
    // and variant_afreqcalc is the set of variants which were seen but
    // skipped due to incomplete information.
    BitvecOr(already_seen, raw_variant_ctl, *variant_afreqcalcp);
    putc_unlocked('\r', stdout);
    logprintf("--read-freq: Frequencies for %u variant%s loaded.\n", loaded_variant_ct, (loaded_variant_ct == 1)? "" : "s");
    if (skipped_variant_ct) {
      logerrprintfww("Warning: %" PRIuPTR " entr%s skipped due to missing variant IDs, mismatching allele codes, and/or zero observations.\n", skipped_variant_ct, (skipped_variant_ct == 1)? "y" : "ies");
    }
  }
  while (0) {
  ReadAlleleFreqs_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--read-freq file", &read_freq_txs);
    break;
  ReadAlleleFreqs_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ReadAlleleFreqs_ret_UNRECOGNIZED_HEADER:
    logerrputs("Error: Unrecognized header line in --read-freq file.\n");
    reterr = kPglRetMalformedInput;
    break;
  ReadAlleleFreqs_ret_MISSING_TOKENS:
    logerrprintfww("Error: Line %" PRIuPTR " of --read-freq file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  ReadAlleleFreqs_ret_INVALID_FREQS:
    snprintf(g_logbuf, kLogbufSize, "Error: Invalid frequencies/counts on line %" PRIuPTR " of --read-freq file.\n", line_idx);
  ReadAlleleFreqs_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
  ReadAlleleFreqs_ret_MALFORMED_INPUT_2:
    logputs("\n");
    logerrputsb();
  ReadAlleleFreqs_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
 ReadAlleleFreqs_ret_1:
  CleanupTextStream2("--read-freq file", &read_freq_txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}
PglErr SelectSidRepresentatives(const uintptr_t* sex_nm, const uintptr_t* sex_male, const uint32_t* sample_missing_cts, const uint32_t* sample_hethap_cts, const PedigreeIdInfo* piip, SelectSidMissingnessMode missingness_mode, SelectSidTiebreakMode tiebreak_mode, uint32_t parents_only, uint32_t raw_sample_ct, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t orig_sample_ct = *sample_ct_ptr;
    if (!orig_sample_ct) {
      goto SelectSidRepresentatives_ret_1;
    }
    const SampleIdInfo* siip = &(piip->sii);
    if (!siip->sids) {
      logputs("Skipping --select-sid-representatives: all FID+IID pairs are unique.\n");
      goto SelectSidRepresentatives_ret_1;
    }
    const uint32_t use_nsort = (tiebreak_mode == kSelectSidTiebreakFirst) || (tiebreak_mode == kSelectSidTiebreakLast);
    char* sorted_3idbox;
    uint32_t* xid_map;
    uintptr_t max_3id_blen;
    reterr = SortedXidboxInitAlloc(sample_include, siip, orig_sample_ct, kfXidModeFidIidSid, use_nsort, &sorted_3idbox, &xid_map, &max_3id_blen);
    if (unlikely(reterr)) {
      goto SelectSidRepresentatives_ret_1;
    }
    char* idbuf;
    if (unlikely(bigstack_alloc_c(siip->max_sample_id_blen + 1, &idbuf))) {
      goto SelectSidRepresentatives_ret_NOMEM;
    }
    uintptr_t* all_parents = nullptr;
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    if (parents_only) {
      if (unlikely(bigstack_calloc_w(raw_sample_ctl, &all_parents))) {
        goto SelectSidRepresentatives_ret_NOMEM;
      }
      MarkParents(&(piip->parental_id_info), sorted_3idbox, xid_map, orig_sample_ct, max_3id_blen, use_nsort, all_parents, idbuf);
      if (AllWordsAreZero(all_parents, raw_sample_ctl)) {
        logputs("Skipping --select-sid-representatives parents-only: no parents remaining.\n");
        goto SelectSidRepresentatives_ret_1;
      }
    }
    if (missingness_mode == kSelectSidMissingness0) {
      sample_missing_cts = nullptr;
    }
    const int32_t tiebreak_dir = (tiebreak_mode >= kSelectSidTiebreakLast)? -1 : 1;
    uint32_t xid_idx_end = 0;
    do {
      const uint32_t xid_idx_start = xid_idx_end;
      const char* fid_start = &(sorted_3idbox[xid_idx_start * max_3id_blen]);
      const char* iid_end = AdvToDelim(AdvPastDelim(fid_start, '\t'), '\t');
      const uint32_t sample_2id_slen = iid_end - fid_start;
      memcpy(idbuf, fid_start, sample_2id_slen);
      idbuf[sample_2id_slen] = ' ';
      const uint32_t sample_2id_blen = sample_2id_slen + 1;
      idbuf[sample_2id_blen] = '\0';
      if (use_nsort) {
        xid_idx_end = ExpsearchNsortStrLb(idbuf, sorted_3idbox, max_3id_blen, orig_sample_ct, xid_idx_start + 1);
      } else {
        xid_idx_end = ExpsearchStrLb(idbuf, sorted_3idbox, sample_2id_blen, max_3id_blen, orig_sample_ct, xid_idx_start + 1);
      }
      if (xid_idx_end == xid_idx_start + 1) {
        continue;
      }
      uint32_t sample_uidx = xid_map[xid_idx_start];
      if (all_parents && (!IsSet(all_parents, sample_uidx))) {
        continue;
      }
      // intentional underflow on xid_idx_start == 0
      uint32_t xid_idx = (tiebreak_dir == 1)? (xid_idx_start - 1) : xid_idx_end;
      const uint32_t xid_idx_stop = (tiebreak_dir == 1)? (xid_idx_end - 1) : xid_idx_start;
      const uint32_t sex_nm_plus_male = IsSet(sex_nm, sample_uidx) + IsSet(sex_male, sample_uidx);
      uint32_t best_missing_ct = UINT32_MAX;
      // keep correct sample in sid-only case
      uint32_t cur_missing_ct = UINT32_MAXM1;
      uint32_t best_sample_uidx = 0;
      do {
        xid_idx += tiebreak_dir;
        sample_uidx = xid_map[xid_idx];
        if (unlikely(sex_nm_plus_male != IsSet(sex_nm, sample_uidx) + IsSet(sex_male, sample_uidx))) {
          idbuf[sample_2id_slen] = '\0';
          TabsToSpaces(idbuf);
          logerrprintf("Error: Inconsistent sex information for sample \"%s\".\n", idbuf);
          goto SelectSidRepresentatives_ret_INCONSISTENT_INPUT;
        }
        ClearBit(sample_uidx, sample_include);
        if (sample_missing_cts) {
          cur_missing_ct = sample_missing_cts[sample_uidx];
          if (sample_hethap_cts) {
            cur_missing_ct += sample_hethap_cts[sample_uidx];
          }
        }
        if (cur_missing_ct >= best_missing_ct) {
          continue;
        }
        best_missing_ct = cur_missing_ct;
        best_sample_uidx = sample_uidx;
      } while (xid_idx != xid_idx_stop);
      SetBit(best_sample_uidx, sample_include);
    } while (xid_idx_end != orig_sample_ct);
    const uint32_t sample_ct = PopcountWords(sample_include, raw_sample_ctl);
    const uint32_t removed_ct = orig_sample_ct - sample_ct;
    logprintf("--select-sid-representatives: %u sample%s removed.\n", removed_ct, (removed_ct == 1)? "" : "s");
    *sample_ct_ptr = sample_ct;
  }
  while (0) {
  SelectSidRepresentatives_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  SelectSidRepresentatives_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 SelectSidRepresentatives_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}
// Generalization of plink 1.9 load_ax_alleles().
//
// Note that, when a variant has 2 (or more) nonmissing allele codes, we print
// a warning if the loaded allele doesn't match one of them (and the variant is
// unchanged).  However, if the variant is biallelic with one or two missing
// allele codes, a missing allele code is filled in.  This maintains
// compatibility with plink 1.9 --a1-allele/--a2-allele's behavior, and is good
// enough for most real-world use cases; however, it's a bit annoying to be
// unable to e.g. add the ref allele when dealing with a single-sample file
// with a het alt1/alt2 call.
// (Workaround for that case when merge is implemented: generate a
// single-sample file with all the right reference alleles, and merge with
// that.)
PglErr MakeFounders(const uintptr_t* sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t require_two, PedigreeIdInfo* piip, uintptr_t* founder_info) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uintptr_t max_sample_id_blen = piip->sii.max_sample_id_blen;
    char* idbuf;
    uintptr_t* nf_bitvec;
    if (bigstack_alloc_c(max_sample_id_blen, &idbuf) ||
        bigstack_alloc_w(raw_sample_ctl, &nf_bitvec)) {
      goto MakeFounders_ret_NOMEM;
    }
    BitvecInvmaskCopy(sample_include, founder_info, raw_sample_ctl, nf_bitvec);
    uint32_t sample_uidx = AdvBoundedTo1Bit(nf_bitvec, 0, raw_sample_ct);
    if (sample_uidx == raw_sample_ct) {
      logputs("Note: Skipping --make-founders since there are no nonfounders.\n");
      goto MakeFounders_ret_1;
    }
    const char* sample_ids = piip->sii.sample_ids;
    uint32_t* id_htable;
    uint32_t id_htable_size;
    if (unlikely(HtableGoodSizeAlloc(sample_ct, bigstack_left(), &id_htable, &id_htable_size))) {
      goto MakeFounders_ret_NOMEM;
    }
    PopulateStrboxSubsetHtable(sample_ids, sample_include, sample_ct, max_sample_id_blen, 1, id_htable_size, id_htable);
    const uintptr_t max_paternal_id_blen = piip->parental_id_info.max_paternal_id_blen;
    const uintptr_t max_maternal_id_blen = piip->parental_id_info.max_maternal_id_blen;
    char* paternal_ids = piip->parental_id_info.paternal_ids;
    char* maternal_ids = piip->parental_id_info.maternal_ids;
    uint32_t new_founder_ct = 0;
    do {
      const char* fid_start = &(sample_ids[sample_uidx * max_sample_id_blen]);
      const char* fid_postend = AdvPastDelim(fid_start, '\t');
      const uint32_t fid_blen = fid_postend - fid_start;
      char* iid_write_start = memcpya(idbuf, fid_start, fid_blen);

      uint32_t missing_parent_ct = 0;
      char* pat_ptr = &(paternal_ids[sample_uidx * max_paternal_id_blen]);
      const uint32_t pat_slen = strlen(pat_ptr);
      if (fid_blen + pat_slen >= max_sample_id_blen) {
        ++missing_parent_ct;
      } else {
        memcpy(iid_write_start, pat_ptr, pat_slen + 1);
        const uint32_t uii = StrboxHtableFind(idbuf, sample_ids, id_htable, max_sample_id_blen, fid_blen + pat_slen, id_htable_size);
        if (uii == UINT32_MAX) {
          ++missing_parent_ct;
        }
      }
      char* mat_ptr = &(maternal_ids[sample_uidx * max_maternal_id_blen]);
      const uint32_t mat_slen = strlen(mat_ptr);
      if (fid_blen + mat_slen >= max_sample_id_blen) {
        ++missing_parent_ct;
      } else {
        memcpy(iid_write_start, mat_ptr, mat_slen + 1);
        const uint32_t uii = StrboxHtableFind(idbuf, sample_ids, id_htable, max_sample_id_blen, fid_blen + mat_slen, id_htable_size);
        if (uii == UINT32_MAX) {
          ++missing_parent_ct;
        }
      }
      if (missing_parent_ct > require_two) {
        SetBit(sample_uidx, founder_info);
        strcpy_k(pat_ptr, "0");
        strcpy_k(mat_ptr, "0");
        ++new_founder_ct;
      }
      sample_uidx = AdvBoundedTo1Bit(nf_bitvec, sample_uidx + 1, raw_sample_ct);
    } while (sample_uidx != raw_sample_ct);
    logprintf("--make-founders: %u sample%s affected.\n", new_founder_ct, (new_founder_ct == 1)? "" : "s");
  }
  while (0) {
  MakeFounders_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
 MakeFounders_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}

static_assert(sizeof(Dosage) + sizeof(AlleleCode) <= 4, "MajRef() needs to be updated.");
#ifdef __cplusplus
}  // namespace plink2
#endif
