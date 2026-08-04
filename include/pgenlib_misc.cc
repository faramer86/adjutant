// This library is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.

#include "pgenlib_misc.h"

#include <limits.h>

#ifdef __cplusplus
namespace plink2 {
#endif

#ifndef NDEBUG
char* g_pgl_errbuf = nullptr;
char* g_pgl_errbuf_write_iter = nullptr;
char* g_pgl_errbuf_end = nullptr;
#endif

#ifdef USE_AVX2
void CopyNyparrNonemptySubset(const uintptr_t* __restrict raw_nyparr, const uintptr_t* __restrict subset_mask, uint32_t raw_nyparr_entry_ct, uint32_t subset_entry_ct, uintptr_t* __restrict output_nyparr) {
  if (subset_entry_ct == raw_nyparr_entry_ct) {
    memcpy(output_nyparr, raw_nyparr, DivUp(subset_entry_ct, kBitsPerWordD2) * sizeof(intptr_t));
    ZeroTrailingNyps(subset_entry_ct, output_nyparr);
    return;
  }
  assert(subset_entry_ct);
  uintptr_t cur_output_word = 0;

  uintptr_t* output_nyparr_iter = output_nyparr;

  uintptr_t* output_nyparr_last = &(output_nyparr[subset_entry_ct / kBitsPerWordD2]);
  const uint32_t word_write_shift_end = 2 * (subset_entry_ct % kBitsPerWordD2);
  uint32_t word_write_shift = 0;
  for (uint32_t subset_mask_widx = 0; ; ++subset_mask_widx) {
    const uintptr_t cur_include_word = subset_mask[subset_mask_widx];
    if (cur_include_word) {
      uint32_t cur_include_halfword = S_CAST(Halfword, cur_include_word);
      for (uint32_t wordhalf_idx = 0; ; ++wordhalf_idx) {
        if (cur_include_halfword) {
          uintptr_t extracted_bits = raw_nyparr[subset_mask_widx * 2 + wordhalf_idx];
          uint32_t set_bit_ct = kBitsPerWord;
          if (cur_include_halfword != UINT32_MAX) {
            const uintptr_t pext_mask = 3 * UnpackHalfwordToWord(cur_include_halfword);
            extracted_bits = _pext_u64(extracted_bits, pext_mask);
            set_bit_ct = PopcountWord(pext_mask);
          }
          cur_output_word |= extracted_bits << word_write_shift;
          word_write_shift += set_bit_ct;
          if (word_write_shift >= kBitsPerWord) {
            *output_nyparr_iter++ = cur_output_word;
            word_write_shift -= kBitsPerWord;
            cur_output_word = 0;
            if (word_write_shift) {
              cur_output_word = extracted_bits >> (set_bit_ct - word_write_shift);
            }
          }
        }
        if (wordhalf_idx) {
          break;
        }
        cur_include_halfword = cur_include_word >> kBitsPerWordD2;
      }
      if (output_nyparr_iter == output_nyparr_last) {
        if (word_write_shift == word_write_shift_end) {
          if (word_write_shift_end) {
            *output_nyparr_last = cur_output_word;
          }
          return;
        }
      }
    }
  }
}

// bit_idx_start assumed to be < kBitsPerWord
// Variant of ExpandBytearr() which is based off a target 2-bit value instead
// of single expand_mask bits.  expand_size must be the number of instances of
// the target value in genovec.
void ExpandBytearrFromGenoarr(const void* __restrict compact_bitarr, const uintptr_t* __restrict genoarr, uintptr_t match_word, uint32_t genoword_ct, uint32_t expand_size, uint32_t read_start_bit, uintptr_t* __restrict target) {
  const uint32_t expand_sizex_m1 = expand_size + read_start_bit - 1;
  const uint32_t leading_byte_ct = 1 + (expand_sizex_m1 % kBitsPerWord) / CHAR_BIT;
  const uint32_t genoword_ct_m1 = genoword_ct - 1;
  uintptr_t compact_word = SubwordLoad(compact_bitarr, leading_byte_ct) >> read_start_bit;
  const unsigned char* compact_bitarr_biter = &(S_CAST(const unsigned char*, compact_bitarr)[leading_byte_ct]);
  uint32_t compact_idx_lowbits = read_start_bit + CHAR_BIT * (sizeof(intptr_t) - leading_byte_ct);
  for (uint32_t widx = 0; ; widx += 2) {
    uintptr_t mask_word;
    if (widx >= genoword_ct_m1) {
      if (widx > genoword_ct_m1) {
        return;
      }
      mask_word = 0;
    } else {
      const uintptr_t geno_word1 = genoarr[widx + 1] ^ match_word;
      mask_word = PackWordToHalfwordMask5555(~(geno_word1 | (geno_word1 >> 1)));
      mask_word = mask_word << 32;
    }
    const uintptr_t geno_word0 = genoarr[widx] ^ match_word;
    mask_word |= PackWordToHalfwordMask5555(~(geno_word0 | (geno_word0 >> 1)));
    uintptr_t write_word = 0;
    if (mask_word) {
      const uint32_t mask_set_ct = PopcountWord(mask_word);
      uint32_t next_compact_idx_lowbits = compact_idx_lowbits + mask_set_ct;
      if (next_compact_idx_lowbits <= kBitsPerWord) {
        write_word = _pdep_u64(compact_word, mask_word);
        if (mask_set_ct != kBitsPerWord) {
          compact_word = compact_word >> mask_set_ct;
        } else {
          // avoid nasal demons
          compact_word = 0;
        }
      } else {
        uintptr_t next_compact_word;
        CopyFromUnalignedIncrW(&next_compact_word, &compact_bitarr_biter);
        next_compact_idx_lowbits -= kBitsPerWord;
        compact_word |= next_compact_word << (kBitsPerWord - compact_idx_lowbits);
        write_word = _pdep_u64(compact_word, mask_word);
        if (next_compact_idx_lowbits != kBitsPerWord) {
          compact_word = next_compact_word >> next_compact_idx_lowbits;
        } else {
          compact_word = 0;
        }
      }
      compact_idx_lowbits = next_compact_idx_lowbits;
    }
    target[widx / 2] = write_word;
  }
}
#else  // !USE_AVX2
void CopyNyparrNonemptySubset(const uintptr_t* __restrict raw_nyparr, const uintptr_t* __restrict subset_mask, uint32_t raw_nyparr_entry_ct, uint32_t subset_entry_ct, uintptr_t* __restrict output_nyparr) {
  if (subset_entry_ct == raw_nyparr_entry_ct) {
    // subset_mask may be nullptr in this case
    memcpy(output_nyparr, raw_nyparr, DivUp(subset_entry_ct, kBitsPerWordD2) * sizeof(intptr_t));
    ZeroTrailingNyps(subset_entry_ct, output_nyparr);
    return;
  }
  assert(subset_entry_ct);
  assert(raw_nyparr_entry_ct >= subset_entry_ct);
  uintptr_t cur_output_word = 0;

  uintptr_t* output_nyparr_iter = output_nyparr;

  uintptr_t* output_nyparr_last = &(output_nyparr[subset_entry_ct / kBitsPerWordD2]);
  const uint32_t word_write_halfshift_end = subset_entry_ct % kBitsPerWordD2;
  uint32_t word_write_halfshift = 0;
  // if <= 2/3-filled, use sparse copy algorithm
  // (tried CopyBitarrSubset() approach, that actually worsened things)
  if (subset_entry_ct * (3 * k1LU) <= raw_nyparr_entry_ct * (2 * k1LU)) {
    for (uint32_t subset_mask_widx = 0; ; ++subset_mask_widx) {
      const uintptr_t cur_include_word = subset_mask[subset_mask_widx];
      if (cur_include_word) {
        uint32_t cur_include_halfword = S_CAST(Halfword, cur_include_word);
        for (uint32_t wordhalf_idx = 0; ; ++wordhalf_idx) {
          if (cur_include_halfword) {
            uintptr_t raw_nyparr_word = raw_nyparr[subset_mask_widx * 2 + wordhalf_idx];
            do {
              uint32_t rqa_idx_lowbits = ctzu32(cur_include_halfword);
              cur_output_word |= ((raw_nyparr_word >> (rqa_idx_lowbits * 2)) & 3) << (word_write_halfshift * 2);
              if (++word_write_halfshift == kBitsPerWordD2) {
                *output_nyparr_iter++ = cur_output_word;
                word_write_halfshift = 0;
                cur_output_word = 0;
              }
              cur_include_halfword &= cur_include_halfword - 1;
            } while (cur_include_halfword);
          }
          if (wordhalf_idx) {
            break;
          }
          cur_include_halfword = cur_include_word >> kBitsPerWordD2;
        }
        if (output_nyparr_iter == output_nyparr_last) {
          if (word_write_halfshift == word_write_halfshift_end) {
            if (word_write_halfshift_end) {
              *output_nyparr_last = cur_output_word;
            }
            return;
          }
        }
      }
    }
  }
  // blocked copy
  const uintptr_t* raw_nyparr_iter = raw_nyparr;
  for (; ; ++subset_mask) {
    const uintptr_t cur_include_word = *subset_mask;
    uintptr_t cur_include_halfword = S_CAST(Halfword, cur_include_word);
    for (uint32_t wordhalf_idx = 0; ; ++wordhalf_idx) {
      uintptr_t raw_nyparr_word = *raw_nyparr_iter++;
      while (cur_include_halfword) {
        const uint32_t rqa_idx_lowbits = ctzw(cur_include_halfword);

        // TAOCP, 7.1.3, (43).
        const uintptr_t bottom_block_remover = (cur_include_halfword | (cur_include_halfword - 1)) + 1;

        const uintptr_t raw_nyparr_curblock_unmasked = raw_nyparr_word >> (rqa_idx_lowbits * 2);
        const uint32_t rqa_block_len = ctzw(bottom_block_remover) - rqa_idx_lowbits;
        const uint32_t block_len_limit = kBitsPerWordD2 - word_write_halfshift;
        cur_output_word |= raw_nyparr_curblock_unmasked << (2 * word_write_halfshift);
        if (rqa_block_len < block_len_limit) {
          word_write_halfshift += rqa_block_len;
          cur_output_word = bzhi(cur_output_word, word_write_halfshift * 2);
        } else {
          // no need to mask, extra bits vanish off the high end
          *output_nyparr_iter++ = cur_output_word;
          word_write_halfshift = rqa_block_len - block_len_limit;
          if (word_write_halfshift) {
            cur_output_word = bzhi(raw_nyparr_curblock_unmasked >> (2 * block_len_limit), 2 * word_write_halfshift);
          } else {
            // avoid potential right-shift-[word length]
            cur_output_word = 0;
          }
        }
        cur_include_halfword &= bottom_block_remover;
      }
      if (wordhalf_idx) {
        break;
      }
      cur_include_halfword = cur_include_word >> kBitsPerWordD2;
    }
    if (output_nyparr_iter == output_nyparr_last) {
      if (word_write_halfshift == word_write_halfshift_end) {
        if (word_write_halfshift_end) {
          *output_nyparr_last = cur_output_word;
        }
        return;
      }
    }
  }
}
void ExpandBytearrFromGenoarr(const void* __restrict compact_bitarr, const uintptr_t* __restrict genoarr, uintptr_t match_word, uint32_t genoword_ct, uint32_t expand_size, uint32_t read_start_bit, uintptr_t* __restrict target) {
  Halfword* target_alias = DowncastWToHW(target);
  ZeroHwArr(RoundUpPow2(genoword_ct, 2), target_alias);
  const unsigned char* compact_bitarr_alias = S_CAST(const unsigned char*, compact_bitarr);
  const uint32_t expand_sizex_m1 = expand_size + read_start_bit - 1;
  const uint32_t compact_widx_last = expand_sizex_m1 / kBitsPerWord;
  uint32_t compact_idx_lowbits = read_start_bit;
  uint32_t loop_len = kBitsPerWord;
  uintptr_t write_hwidx = 0;
  uintptr_t genomatch_bits = genoarr[0] ^ match_word;
  genomatch_bits = (~(genomatch_bits | (genomatch_bits >> 1))) & kMask5555;
  for (uint32_t compact_widx = 0; ; ++compact_widx) {
    uintptr_t compact_word;
    if (compact_widx >= compact_widx_last) {
      if (compact_widx > compact_widx_last) {
        return;
      }
      loop_len = 1 + (expand_sizex_m1 % kBitsPerWord);
      // avoid possible segfault
      compact_word = SubwordLoad(&(compact_bitarr_alias[compact_widx * kBytesPerWord]), DivUp(loop_len, CHAR_BIT));
    } else {
      CopyFromUnalignedOffsetW(&compact_word, compact_bitarr_alias, compact_widx);
    }
    for (; compact_idx_lowbits != loop_len; ++compact_idx_lowbits) {
      while (!genomatch_bits) {
        genomatch_bits = genoarr[++write_hwidx] ^ match_word;
        genomatch_bits = (~(genomatch_bits | (genomatch_bits >> 1))) & kMask5555;
      }
      if (compact_word & (k1LU << compact_idx_lowbits)) {
        const uint32_t lowbit_idx = ctzw(genomatch_bits);
        target_alias[write_hwidx] |= 1U << (lowbit_idx / 2);
      }
      genomatch_bits &= genomatch_bits - 1;
    }
    compact_idx_lowbits = 0;
  }
}
#endif

// Harley-Seal algorithm only works for bitarrays, not nyparrays, so don't
// add an AVX2 specialization here.
// ...unless something like the interleaved_vec strategy is used?  hmm.  should
// test this on basic frequency counter.
/*
void Count2FreqVec3(const VecW* geno_vvec, uint32_t vec_ct, uint32_t* __restrict alt1_plus_bothset_ctp, uint32_t* __restrict bothset_ctp) {
  assert(!(vec_ct % 3));
  // Increments bothset_ct by the number of 0b11 in the current block, and
  // alt1_ct by twice the number of 0b10 plus the number of 0b01.
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW* geno_vvec_iter = geno_vvec;
  uint32_t alt1_plus_bothset_ct = 0;
  uint32_t bothset_ct = 0;

  while (1) {
    UniVec acc_alt1_plus_bothset;
    UniVec acc_bothset;
    acc_alt1_plus_bothset.vw = vecw_setzero();
    acc_bothset.vw = vecw_setzero();
    const VecW* geno_vvec_stop;
    if (vec_ct < 30) {
      if (!vec_ct) {
        *alt1_plus_bothset_ctp = alt1_plus_bothset_ct;
        *bothset_ctp = bothset_ct;
        return;
      }
      geno_vvec_stop = &(geno_vvec_iter[vec_ct]);
      vec_ct = 0;
    } else {
      geno_vvec_stop = &(geno_vvec_iter[30]);
      vec_ct -= 30;
    }
    do {
      VecW cur_geno_vword1 = *geno_vvec_iter++;
      // process first two vwords simultaneously to minimize linear dependence
      VecW cur_geno_vword2 = *geno_vvec_iter++;
      VecW cur_geno_vword_low_lshifted1 = vecw_slli(cur_geno_vword1 & m1, 1);
      VecW cur_geno_vword_low_lshifted2 = vecw_slli(cur_geno_vword2 & m1, 1);

      // 00 -> 00; 01 -> 01; 10 -> 10; 11 -> 01
      // note that _mm_andnot_si128 flips the *first* argument before the AND
      // operation.
      VecW alt1_plus_bothset1 = vecw_and_notfirst(cur_geno_vword_low_lshifted1, cur_geno_vword1);
      VecW alt1_plus_bothset2 = vecw_and_notfirst(cur_geno_vword_low_lshifted2, cur_geno_vword2);

      VecW bothset1 = vecw_srli(cur_geno_vword_low_lshifted1 & cur_geno_vword1, 1);
      VecW bothset2 = vecw_srli(cur_geno_vword_low_lshifted2 & cur_geno_vword2, 1);

      cur_geno_vword1 = *geno_vvec_iter++;
      alt1_plus_bothset1 = (alt1_plus_bothset1 & m2) + (vecw_srli(alt1_plus_bothset1, 2) & m2);
      bothset2 = bothset1 + bothset2;
      alt1_plus_bothset2 = (alt1_plus_bothset2 & m2) + (vecw_srli(alt1_plus_bothset2, 2) & m2);
      cur_geno_vword_low_lshifted1 = vecw_slli(cur_geno_vword1 & m1, 1);

      alt1_plus_bothset2 = alt1_plus_bothset1 + alt1_plus_bothset2;
      // alt1_plus_bothset2 now contains 4-bit values from 0-8, while bothset2
      // contains 2-bit values from 0-2
      // (todo: check whether this is faster if we use double_bothsetx
      // variables instead of bothset1/bothset2)
      bothset1 = vecw_srli(cur_geno_vword_low_lshifted1 & cur_geno_vword1, 1);
      alt1_plus_bothset1 = vecw_and_notfirst(cur_geno_vword_low_lshifted1, cur_geno_vword1);
      bothset2 = bothset1 + bothset2;
      alt1_plus_bothset1 = (alt1_plus_bothset1 & m2) + (vecw_srli(alt1_plus_bothset1, 2) & m2);

      bothset2 = (bothset2 & m2) + (vecw_srli(bothset2, 2) & m2);
      alt1_plus_bothset2 = alt1_plus_bothset1 + alt1_plus_bothset2;
      // alt1_plus_bothset2 now contains 4-bit values from 0-12, while bothset2
      // contains 4-bit values from 0-6.  aggregate both into 8-bit values.
      bothset2 = (bothset2 & m4) + (vecw_srli(bothset2, 4) & m4);
      alt1_plus_bothset2 = (alt1_plus_bothset2 & m4) + (vecw_srli(alt1_plus_bothset2, 4) & m4);

      acc_bothset.vw = acc_bothset.vw + bothset2;
      acc_alt1_plus_bothset.vw = acc_alt1_plus_bothset.vw + alt1_plus_bothset2;
    } while (geno_vvec_iter < geno_vvec_stop);
    const VecW m8 = VCONST_W(kMask00FF);
    acc_bothset.vw = (acc_bothset.vw + vecw_srli(acc_bothset.vw, 8)) & m8;
    acc_alt1_plus_bothset.vw = (acc_alt1_plus_bothset.vw & m8) + (vecw_srli(acc_alt1_plus_bothset.vw, 8) & m8);
    bothset_ct += UniVecHsum16(acc_bothset);
    alt1_plus_bothset_ct += UniVecHsum16(acc_alt1_plus_bothset);
  }
}
*/

// todo: benchmark against general-purpose counter
void Count12Vec6(const VecW* geno_vvec, uint32_t vec_ct, uint32_t* __restrict raw_01_ctp, uint32_t* __restrict raw_both_ctp) {
  assert(!(vec_ct % 6));
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW* geno_vvec_iter = geno_vvec;
  VecW acc_01 = vecw_setzero();
  VecW acc_both = vecw_setzero();
  uintptr_t cur_incr = 60;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 60) {
      if (!vec_ct) {
        *raw_01_ctp = HsumW(acc_01);
        *raw_both_ctp = HsumW(acc_both);
        return;
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc_01 = vecw_setzero();
    VecW inner_acc_both = vecw_setzero();
    const VecW* geno_vvec_stop = &(geno_vvec_iter[cur_incr]);
    do {
      VecW cur_vvec = *geno_vvec_iter++;
      VecW vvec_rshift = vecw_srli(cur_vvec, 1);
      VecW nyp_both = m1 & (cur_vvec ^ vvec_rshift);
      VecW nyp_01 = nyp_both & cur_vvec;

      cur_vvec = *geno_vvec_iter++;
      vvec_rshift = vecw_srli(cur_vvec, 1);
      VecW vvec_both = m1 & (cur_vvec ^ vvec_rshift);
      nyp_01 = nyp_01 + (vvec_both & cur_vvec);
      nyp_both = nyp_both + vvec_both;

      cur_vvec = *geno_vvec_iter++;
      vvec_rshift = vecw_srli(cur_vvec, 1);
      vvec_both = m1 & (cur_vvec ^ vvec_rshift);
      nyp_01 = nyp_01 + (vvec_both & cur_vvec);
      nyp_both = nyp_both + vvec_both;

      VecW nybble_01 = (nyp_01 & m2) + (vecw_srli(nyp_01, 2) & m2);
      VecW nybble_both = (nyp_both & m2) + (vecw_srli(nyp_both, 2) & m2);

      cur_vvec = *geno_vvec_iter++;
      vvec_rshift = vecw_srli(cur_vvec, 1);
      nyp_both = m1 & (cur_vvec ^ vvec_rshift);
      nyp_01 = nyp_both & cur_vvec;

      cur_vvec = *geno_vvec_iter++;
      vvec_rshift = vecw_srli(cur_vvec, 1);
      vvec_both = m1 & (cur_vvec ^ vvec_rshift);
      nyp_01 = nyp_01 + (vvec_both & cur_vvec);
      nyp_both = nyp_both + vvec_both;

      cur_vvec = *geno_vvec_iter++;
      vvec_rshift = vecw_srli(cur_vvec, 1);
      vvec_both = m1 & (cur_vvec ^ vvec_rshift);
      nyp_01 = nyp_01 + (vvec_both & cur_vvec);
      nyp_both = nyp_both + vvec_both;

      nybble_01 = nybble_01 + (nyp_01 & m2) + (vecw_srli(nyp_01, 2) & m2);
      nybble_both = nybble_both + (nyp_both & m2) + (vecw_srli(nyp_both, 2) & m2);

      inner_acc_01 = inner_acc_01 + (nybble_01 & m4) + (vecw_srli(nybble_01, 4) & m4);
      inner_acc_both = inner_acc_both + (nybble_both & m4) + (vecw_srli(nybble_both, 4) & m4);
    } while (geno_vvec_iter < geno_vvec_stop);
    const VecW m0 = vecw_setzero();
    acc_01 = acc_01 + vecw_bytesum(inner_acc_01, m0);
    acc_both = acc_both + vecw_bytesum(inner_acc_both, m0);
  }
}

void GenovecCount12Unsafe(const uintptr_t* genovec, uint32_t sample_ct, uint32_t* __restrict raw_01_ctp, uint32_t* __restrict raw_10_ctp) {
  // assumes trailing bits of last genovec word are zeroed out.
  // sample_ct == 0 ok.
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  uint32_t raw_01_ct;
  uint32_t raw_both_ct;
  uint32_t word_idx = sample_ctl2 - (sample_ctl2 % (6 * kWordsPerVec));
  assert(IsVecAligned(genovec));
  Count12Vec6(R_CAST(const VecW*, genovec), word_idx / kWordsPerVec, &raw_01_ct, &raw_both_ct);
  for (; word_idx != sample_ctl2; ++word_idx) {
    const uintptr_t cur_geno_word = genovec[word_idx];
    const uintptr_t cur_rshift = cur_geno_word >> 1;
    const uintptr_t cur_both = kMask5555 & (cur_geno_word ^ cur_rshift);
    raw_01_ct += Popcount01Word(cur_both & cur_geno_word);
    raw_both_ct += Popcount01Word(cur_both);
  }
  *raw_01_ctp = raw_01_ct;
  *raw_10_ctp = raw_both_ct - raw_01_ct;
}

void Count3FreqVec6(const void* genoarr, uint32_t vec_ct, uint32_t* __restrict even_ctp, uint32_t* __restrict odd_ctp, uint32_t* __restrict bothset_ctp) {
  assert(!(vec_ct % 6));
  // Sets even_ct to the number of set low bits in the current block, odd_ct to
  // the number of set high bits, and bothset_ct by the number of 0b11s.
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const unsigned char* geno_vvec_biter = S_CAST(const unsigned char*, genoarr);
  VecW acc_even = vecw_setzero();
  VecW acc_odd = vecw_setzero();
  VecW acc_bothset = vecw_setzero();
  uintptr_t cur_incr = 60;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 60) {
      if (!vec_ct) {
        *even_ctp = HsumW(acc_even);
        *odd_ctp = HsumW(acc_odd);
        *bothset_ctp = HsumW(acc_bothset);
        return;
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc_even = vecw_setzero();
    VecW inner_acc_odd = vecw_setzero();
    VecW inner_acc_bothset = vecw_setzero();
    const unsigned char* geno_vvec_stop = &(geno_vvec_biter[cur_incr * kBytesPerVec]);
    do {
      // hmm, this seems to have more linear dependence than I'd want, but the
      // reorderings I tried just made the code harder to read without helping,
      // so I'll leave this alone
      VecW cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      VecW odd1 = m1 & vecw_srli(cur_geno_vword, 1);
      VecW even1 = m1 & cur_geno_vword;
      VecW bothset1 = odd1 & cur_geno_vword;

      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      VecW cur_geno_vword_high = m1 & vecw_srli(cur_geno_vword, 1);
      even1 = even1 + (m1 & cur_geno_vword);
      odd1 = odd1 + cur_geno_vword_high;
      bothset1 = bothset1 + (cur_geno_vword_high & cur_geno_vword);

      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      cur_geno_vword_high = m1 & vecw_srli(cur_geno_vword, 1);
      even1 = even1 + (m1 & cur_geno_vword);
      odd1 = odd1 + cur_geno_vword_high;
      bothset1 = bothset1 + (cur_geno_vword_high & cur_geno_vword);

      even1 = (even1 & m2) + (vecw_srli(even1, 2) & m2);
      odd1 = (odd1 & m2) + (vecw_srli(odd1, 2) & m2);
      bothset1 = (bothset1 & m2) + (vecw_srli(bothset1, 2) & m2);

      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      VecW odd2 = m1 & vecw_srli(cur_geno_vword, 1);
      VecW even2 = m1 & cur_geno_vword;
      VecW bothset2 = odd2 & cur_geno_vword;

      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      cur_geno_vword_high = m1 & vecw_srli(cur_geno_vword, 1);
      even2 = even2 + (m1 & cur_geno_vword);
      odd2 = odd2 + cur_geno_vword_high;
      bothset2 = bothset2 + (cur_geno_vword_high & cur_geno_vword);

      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      cur_geno_vword_high = m1 & vecw_srli(cur_geno_vword, 1);
      even2 = even2 + (m1 & cur_geno_vword);
      odd2 = odd2 + cur_geno_vword_high;
      bothset2 = bothset2 + (cur_geno_vword_high & cur_geno_vword);

      even1 = even1 + (even2 & m2) + (vecw_srli(even2, 2) & m2);
      odd1 = odd1 + (odd2 & m2) + (vecw_srli(odd2, 2) & m2);
      bothset1 = bothset1 + (bothset2 & m2) + (vecw_srli(bothset2, 2) & m2);
      // these now contain 4-bit values from 0-12

      inner_acc_even = inner_acc_even + (even1 & m4) + (vecw_srli(even1, 4) & m4);
      inner_acc_odd = inner_acc_odd + (odd1 & m4) + (vecw_srli(odd1, 4) & m4);
      inner_acc_bothset = inner_acc_bothset + (bothset1 & m4) + (vecw_srli(bothset1, 4) & m4);
    } while (geno_vvec_biter < geno_vvec_stop);
    const VecW m0 = vecw_setzero();
    acc_even = acc_even + vecw_bytesum(inner_acc_even, m0);
    acc_odd = acc_odd + vecw_bytesum(inner_acc_odd, m0);
    acc_bothset = acc_bothset + vecw_bytesum(inner_acc_bothset, m0);
  }
}

void FillInterleavedMaskVec(const uintptr_t* __restrict subset_mask, uint32_t base_vec_ct, uintptr_t* interleaved_mask_vec) {
#ifdef USE_SSE2
  // This is a cousin of github.com/KWillets/simd_interleave , which was
  // written in response to
  //   https://lemire.me/blog/2018/01/09/how-fast-can-you-bit-interleave-32-bit-integers-simd-edition/
  // AVX2 implementation takes ~40% less time than before, and SSE4.2 takes
  // ~65% less.  This also avoids the Ryzen-screwing _pdep_u64()/_pext_u64()
  // operations.
  const VecW m4 = VCONST_W(kMask0F0F);
#  ifdef USE_SHUFFLE8
  const VecW lookup0 = vecw_setr8(
      0, 1, 4, 5, 16, 17, 20, 21,
      64, 65, 68, 69, 80, 81, 84, 85);
  const VecW lookup1 = vecw_slli(lookup0, 1);
#  else
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
#  endif
  const VecW* subset_mask_valias = R_CAST(const VecW*, subset_mask);
  VecW* interleaved_mask_valias = R_CAST(VecW*, interleaved_mask_vec);

  for (uint32_t vidx = 0; vidx != base_vec_ct; ++vidx) {
    // I'll assume the compiler can handle this register allocation job.
    VecW vec_lo;
    VecW vec_hi;
    vecw_lo_and_hi_nybbles(subset_mask_valias[vidx], m4, &vec_lo, &vec_hi);
#  ifdef USE_SHUFFLE8
    vec_lo = vecw_shuffle8(lookup0, vec_lo);
    vec_hi = vecw_shuffle8(lookup1, vec_hi);
#  else
    vec_lo = (vec_lo | vecw_slli(vec_lo, 2)) & m2;
    vec_hi = (vec_hi | vecw_slli(vec_hi, 2)) & m2;
    vec_lo = (vec_lo | vecw_slli(vec_lo, 1)) & m1;
    vec_hi = (vec_hi | vecw_slli(vec_hi, 1)) & m1;
    vec_hi = vecw_slli(vec_hi, 1);
#  endif
    interleaved_mask_valias[vidx] = vec_lo | vec_hi;
  }
#else  // !USE_SSE2
  for (uint32_t widx = 0; widx != base_vec_ct; ++widx) {
    const uintptr_t orig_word = subset_mask[widx];
    uintptr_t ww_even = S_CAST(Halfword, orig_word);
    uintptr_t ww_odd = orig_word >> kBitsPerWordD2;
    ww_even = UnpackHalfwordToWord(ww_even);
    ww_odd = UnpackHalfwordToWord(ww_odd);
    interleaved_mask_vec[widx] = ww_even | (ww_odd << 1);
  }
#endif
}

/*
void GenovecAlleleCtsUnsafe(const uintptr_t* genovec, uint32_t sample_ct, uint32_t* __restrict allele_cts, uint32_t* __restrict bothset_ctp) {
  // assumes trailing bits of last genovec word are zeroed out.
  // sets allele_cts[0] to the number of observed ref alleles, and
  // allele_cts[1] to the number of observed alt1s.
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  uint32_t word_idx = sample_ctl2 - (sample_ctl2 % (3 * kWordsPerVec));
  uint32_t alt1_plus_bothset_ct;
  uint32_t bothset_ct;
  assert(IsVecAligned(genovec));
  Count2FreqVec3(R_CAST(const VecW*, genovec), word_idx / kWordsPerVec, &alt1_plus_bothset_ct, &bothset_ct);
  for (; word_idx != sample_ctl2; ++word_idx) {
    const uintptr_t cur_geno_word = genovec[word_idx];
    const uintptr_t cur_geno_word_low_lshifted = (cur_geno_word & kMask5555) << 1;
    alt1_plus_bothset_ct += NypsumWord((~cur_geno_word_low_lshifted) & cur_geno_word);
    bothset_ct += NypsumWord(cur_geno_word_low_lshifted & cur_geno_word);
  }
  const uint32_t alt1_ct = alt1_plus_bothset_ct - bothset_ct;
  allele_cts[0] = (sample_ct - bothset_ct) * 2 - alt1_ct;
  allele_cts[1] = alt1_ct;
  *bothset_ctp = bothset_ct;
}
*/

void GenoarrCountFreqsUnsafe(const uintptr_t* genoarr, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  // fills genocounts[0] with the number of 00s, genocounts[1] with the number
  // of 01s, etc.
  // assumes trailing bits of last genoarr word are zeroed out.
  // sample_ct == 0 ok.
  // no longer requires vector-alignment.
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  uint32_t even_ct;
  uint32_t odd_ct;
  uint32_t bothset_ct;
  uint32_t word_idx = sample_ctl2 - (sample_ctl2 % (6 * kWordsPerVec));
  Count3FreqVec6(genoarr, word_idx / kWordsPerVec, &even_ct, &odd_ct, &bothset_ct);
  for (; word_idx != sample_ctl2; ++word_idx) {
    const uintptr_t cur_geno_word = genoarr[word_idx];
    const uintptr_t cur_geno_word_high = kMask5555 & (cur_geno_word >> 1);
    even_ct += Popcount01Word(cur_geno_word & kMask5555);
    odd_ct += Popcount01Word(cur_geno_word_high);
    bothset_ct += Popcount01Word(cur_geno_word & cur_geno_word_high);
  }
  genocounts[0] = sample_ct + bothset_ct - even_ct - odd_ct;
  genocounts[1] = even_ct - bothset_ct;
  genocounts[2] = odd_ct - bothset_ct;
  genocounts[3] = bothset_ct;
}
// genoarr now allowed to be unaligned.
void CountSubset3FreqVec6(const void* __restrict genoarr, const VecW* __restrict interleaved_mask_vvec, uint32_t vec_ct, uint32_t* __restrict even_ctp, uint32_t* __restrict odd_ctp, uint32_t* __restrict bothset_ctp) {
  assert(!(vec_ct % 6));
  // Sets even_ct to the number of set low bits in the current block after
  // subsetting, odd_ct to the number of set high bits, and bothset_ct by the
  // number of 0b11s.
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const unsigned char* geno_vvec_biter = S_CAST(const unsigned char*, genoarr);
  const VecW* interleaved_mask_vvec_iter = interleaved_mask_vvec;
  VecW acc_even = vecw_setzero();
  VecW acc_odd = vecw_setzero();
  VecW acc_bothset = vecw_setzero();
  uintptr_t cur_incr = 60;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 60) {
      if (!vec_ct) {
        *even_ctp = HsumW(acc_even);
        *odd_ctp = HsumW(acc_odd);
        *bothset_ctp = HsumW(acc_bothset);
        return;
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc_even = vecw_setzero();
    VecW inner_acc_odd = vecw_setzero();
    VecW inner_acc_bothset = vecw_setzero();
    const unsigned char* geno_vvec_stop = &(geno_vvec_biter[cur_incr * kBytesPerVec]);
    do {
      VecW interleaved_mask_vword = *interleaved_mask_vvec_iter++;
      VecW cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      VecW cur_mask = interleaved_mask_vword & m1;
      VecW odd1 = cur_mask & vecw_srli(cur_geno_vword, 1);
      VecW even1 = cur_mask & cur_geno_vword;
      VecW bothset1 = odd1 & cur_geno_vword;

      cur_mask = vecw_srli(interleaved_mask_vword, 1) & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      VecW cur_geno_vword_high_masked = cur_mask & vecw_srli(cur_geno_vword, 1);
      even1 = even1 + (cur_mask & cur_geno_vword);
      odd1 = odd1 + cur_geno_vword_high_masked;
      bothset1 = bothset1 + (cur_geno_vword_high_masked & cur_geno_vword);

      interleaved_mask_vword = *interleaved_mask_vvec_iter++;
      cur_mask = interleaved_mask_vword & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      cur_geno_vword_high_masked = cur_mask & vecw_srli(cur_geno_vword, 1);
      even1 = even1 + (cur_mask & cur_geno_vword);
      odd1 = odd1 + cur_geno_vword_high_masked;
      bothset1 = bothset1 + (cur_geno_vword_high_masked & cur_geno_vword);

      even1 = (even1 & m2) + (vecw_srli(even1, 2) & m2);
      odd1 = (odd1 & m2) + (vecw_srli(odd1, 2) & m2);
      bothset1 = (bothset1 & m2) + (vecw_srli(bothset1, 2) & m2);

      cur_mask = vecw_srli(interleaved_mask_vword, 1) & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      VecW odd2 = cur_mask & vecw_srli(cur_geno_vword, 1);
      VecW even2 = cur_mask & cur_geno_vword;
      VecW bothset2 = odd2 & cur_geno_vword;

      interleaved_mask_vword = *interleaved_mask_vvec_iter++;
      cur_mask = interleaved_mask_vword & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      cur_geno_vword_high_masked = cur_mask & vecw_srli(cur_geno_vword, 1);
      even2 = even2 + (cur_mask & cur_geno_vword);
      odd2 = odd2 + cur_geno_vword_high_masked;
      bothset2 = bothset2 + (cur_geno_vword_high_masked & cur_geno_vword);

      cur_mask = vecw_srli(interleaved_mask_vword, 1) & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      cur_geno_vword_high_masked = cur_mask & vecw_srli(cur_geno_vword, 1);
      even2 = even2 + (cur_mask & cur_geno_vword);
      odd2 = odd2 + cur_geno_vword_high_masked;
      bothset2 = bothset2 + (cur_geno_vword_high_masked & cur_geno_vword);

      even1 = even1 + (even2 & m2) + (vecw_srli(even2, 2) & m2);
      odd1 = odd1 + (odd2 & m2) + (vecw_srli(odd2, 2) & m2);
      bothset1 = bothset1 + (bothset2 & m2) + (vecw_srli(bothset2, 2) & m2);
      // these now contain 4-bit values from 0-12

      inner_acc_even = inner_acc_even + (even1 & m4) + (vecw_srli(even1, 4) & m4);
      inner_acc_odd = inner_acc_odd + (odd1 & m4) + (vecw_srli(odd1, 4) & m4);
      inner_acc_bothset = inner_acc_bothset + (bothset1 & m4) + (vecw_srli(bothset1, 4) & m4);
    } while (geno_vvec_biter < geno_vvec_stop);
    const VecW m0 = vecw_setzero();
    acc_even = acc_even + vecw_bytesum(inner_acc_even, m0);
    acc_odd = acc_odd + vecw_bytesum(inner_acc_odd, m0);
    acc_bothset = acc_bothset + vecw_bytesum(inner_acc_bothset, m0);
  }
}

void GenoarrCountSubsetFreqs(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict sample_include_interleaved_vec, uint32_t raw_sample_ct, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  // fills genocounts[0] with the number of 00s, genocounts[1] with the number
  // of 01s, etc.
  // {raw_}sample_ct == 0 ok.
  const uint32_t raw_sample_ctv2 = NypCtToVecCt(raw_sample_ct);
  uint32_t even_ct;
  uint32_t odd_ct;
  uint32_t bothset_ct;
#ifdef USE_SSE2
  uint32_t vec_idx = raw_sample_ctv2 - (raw_sample_ctv2 % 6);
  CountSubset3FreqVec6(genoarr, R_CAST(const VecW*, sample_include_interleaved_vec), vec_idx, &even_ct, &odd_ct, &bothset_ct);
  const uintptr_t* genoarr_iter = &(genoarr[kWordsPerVec * vec_idx]);
  const uintptr_t* interleaved_mask_iter = &(sample_include_interleaved_vec[(kWordsPerVec / 2) * vec_idx]);
#  ifdef USE_AVX2
  uintptr_t mask_base1 = 0;
  uintptr_t mask_base2 = 0;
  uintptr_t mask_base3 = 0;
  uintptr_t mask_base4 = 0;
  for (; vec_idx != raw_sample_ctv2; ++vec_idx) {
    uintptr_t mask_word1;
    uintptr_t mask_word2;
    uintptr_t mask_word3;
    uintptr_t mask_word4;
    if (!(vec_idx % 2)) {
      mask_base1 = *interleaved_mask_iter++;
      mask_base2 = *interleaved_mask_iter++;
      mask_base3 = *interleaved_mask_iter++;
      mask_base4 = *interleaved_mask_iter++;
      mask_word1 = mask_base1 & kMask5555;
      mask_word2 = mask_base2 & kMask5555;
      mask_word3 = mask_base3 & kMask5555;
      mask_word4 = mask_base4 & kMask5555;
    } else {
      mask_word1 = (mask_base1 >> 1) & kMask5555;
      mask_word2 = (mask_base2 >> 1) & kMask5555;
      mask_word3 = (mask_base3 >> 1) & kMask5555;
      mask_word4 = (mask_base4 >> 1) & kMask5555;
    }
    for (uint32_t vechalf_idx = 0; ; ++vechalf_idx) {
      const uintptr_t cur_geno_word1 = *genoarr_iter++;
      const uintptr_t cur_geno_word2 = *genoarr_iter++;
      const uintptr_t cur_geno_word1_high_masked = mask_word1 & (cur_geno_word1 >> 1);
      const uintptr_t cur_geno_word2_high_masked = mask_word2 & (cur_geno_word2 >> 1);
      even_ct += PopcountWord(((cur_geno_word1 & mask_word1) << 1) | (cur_geno_word2 & mask_word2));
      odd_ct += PopcountWord((cur_geno_word1_high_masked << 1) | cur_geno_word2_high_masked);
      bothset_ct += PopcountWord(((cur_geno_word1 & cur_geno_word1_high_masked) << 1) | (cur_geno_word2 & cur_geno_word2_high_masked));
      if (vechalf_idx) {
        break;
      }
      mask_word1 = mask_word3;
      mask_word2 = mask_word4;
    }
  }
#  else  // not USE_AVX2
  uintptr_t mask_base1 = 0;
  uintptr_t mask_base2 = 0;
  for (; vec_idx != raw_sample_ctv2; ++vec_idx) {
    uintptr_t mask_word1;
    uintptr_t mask_word2;
    if (!(vec_idx % 2)) {
      mask_base1 = *interleaved_mask_iter++;
      mask_base2 = *interleaved_mask_iter++;
      mask_word1 = mask_base1 & kMask5555;
      mask_word2 = mask_base2 & kMask5555;
    } else {
      mask_word1 = (mask_base1 >> 1) & kMask5555;
      mask_word2 = (mask_base2 >> 1) & kMask5555;
    }
    const uintptr_t cur_geno_word1 = *genoarr_iter++;
    const uintptr_t cur_geno_word2 = *genoarr_iter++;
    const uintptr_t cur_geno_word1_high_masked = mask_word1 & (cur_geno_word1 >> 1);
    const uintptr_t cur_geno_word2_high_masked = mask_word2 & (cur_geno_word2 >> 1);
#    ifdef USE_SSE42
    even_ct += PopcountWord(((cur_geno_word1 & mask_word1) << 1) | (cur_geno_word2 & mask_word2));
    odd_ct += PopcountWord((cur_geno_word1_high_masked << 1) | cur_geno_word2_high_masked);
    bothset_ct += PopcountWord(((cur_geno_word1 & cur_geno_word1_high_masked) << 1) | (cur_geno_word2 & cur_geno_word2_high_masked));
#    else
    even_ct += NypsumWord((cur_geno_word1 & mask_word1) + (cur_geno_word2 & mask_word2));
    odd_ct += NypsumWord(cur_geno_word1_high_masked + cur_geno_word2_high_masked);
    bothset_ct += NypsumWord((cur_geno_word1 & cur_geno_word1_high_masked) + (cur_geno_word2 & cur_geno_word2_high_masked));
#    endif
  }
#  endif  // not USE_AVX2
#else  // not USE_SSE2
  uint32_t word_idx = raw_sample_ctv2 - (raw_sample_ctv2 % 6);
  CountSubset3FreqVec6(genoarr, R_CAST(const VecW*, sample_include_interleaved_vec), word_idx, &even_ct, &odd_ct, &bothset_ct);
  const uintptr_t* interleaved_mask_iter = &(sample_include_interleaved_vec[word_idx / 2]);
  uintptr_t mask_base = 0;
  for (; word_idx != raw_sample_ctv2; ++word_idx) {
    uintptr_t mask_word;
    if (!(word_idx % 2)) {
      mask_base = *interleaved_mask_iter++;
      mask_word = mask_base & kMask5555;
    } else {
      mask_word = (mask_base >> 1) & kMask5555;
    }
    const uintptr_t cur_geno_word = genoarr[word_idx];
    const uintptr_t cur_geno_word_high_masked = mask_word & (cur_geno_word >> 1);
    even_ct += Popcount01Word(cur_geno_word & mask_word);
    odd_ct += Popcount01Word(cur_geno_word_high_masked);
    bothset_ct += Popcount01Word(cur_geno_word & cur_geno_word_high_masked);
  }
#endif
  genocounts[0] = sample_ct + bothset_ct - even_ct - odd_ct;
  genocounts[1] = even_ct - bothset_ct;
  genocounts[2] = odd_ct - bothset_ct;
  genocounts[3] = bothset_ct;
}

void GenoarrCountSubsetFreqs2(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  // slower GenoarrCountSubsetFreqs() which does not require
  // sample_include_interleaved_vec to be precomputed.
  // {raw_}sample_ct == 0 ok.
  const uint32_t raw_sample_ctl2 = NypCtToWordCt(raw_sample_ct);
  const uint32_t fullword_ct = raw_sample_ctl2 / 2;
  uint32_t even_ct = 0;
  uint32_t odd_ct = 0;
  uint32_t bothset_ct = 0;
  for (uint32_t widx = 0; widx != fullword_ct; ++widx) {
    // possible todo: try vectorizing this loop in SSE4.2+ high-sample-ct case
    // with shuffle-based dynamic unpacking of sample_include?
    const uintptr_t mask_word = sample_include[widx];
    if (mask_word) {
      uintptr_t geno_word = genoarr[2 * widx];
      uintptr_t geno_even = PackWordToHalfwordMask5555(geno_word);
      uintptr_t geno_odd = PackWordToHalfwordMaskAAAA(geno_word);
      geno_word = genoarr[2 * widx + 1];
      geno_even |= S_CAST(uintptr_t, PackWordToHalfwordMask5555(geno_word)) << kBitsPerWordD2;
      geno_odd |= S_CAST(uintptr_t, PackWordToHalfwordMaskAAAA(geno_word)) << kBitsPerWordD2;
      const uintptr_t geno_even_masked = geno_even & mask_word;
      even_ct += PopcountWord(geno_even_masked);
      odd_ct += PopcountWord(geno_odd & mask_word);
      bothset_ct += PopcountWord(geno_odd & geno_even_masked);
    }
  }
  if (raw_sample_ctl2 % 2) {
    const uintptr_t mask_hw = sample_include[fullword_ct];
    if (mask_hw) {
      const uintptr_t geno_word = genoarr[2 * fullword_ct];
      // todo: benchmark main loop unpack vs. pack
      const uintptr_t mask_word = UnpackHalfwordToWord(mask_hw);
      const uintptr_t geno_word_shifted = geno_word >> 1;
      const uintptr_t geno_word_masked = geno_word & mask_word;
      even_ct += Popcount01Word(geno_word_masked);
      odd_ct += Popcount01Word(geno_word_shifted & mask_word);
      bothset_ct += Popcount01Word(geno_word_masked & geno_word_shifted);
    }
  }
  genocounts[0] = sample_ct + bothset_ct - even_ct - odd_ct;
  genocounts[1] = even_ct - bothset_ct;
  genocounts[2] = odd_ct - bothset_ct;
  genocounts[3] = bothset_ct;
}
void GenoarrCountSubsetIntersectFreqs(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict subset1, const uintptr_t* __restrict subset2, uint32_t raw_sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  // {raw_}sample_ct == 0 ok.
  const uint32_t raw_sample_ctl2 = NypCtToWordCt(raw_sample_ct);
  const uint32_t fullword_ct = raw_sample_ctl2 / 2;
  uint32_t subset_intersect_ct = 0;
  uint32_t even_ct = 0;
  uint32_t odd_ct = 0;
  uint32_t bothset_ct = 0;
  for (uint32_t widx = 0; widx != fullword_ct; ++widx) {
    // hmm, there may be little point to vectorizing this
    const uintptr_t mask_word = subset1[widx] & subset2[widx];
    if (mask_word) {
      uintptr_t geno_word = genoarr[2 * widx];
      uintptr_t geno_even = PackWordToHalfwordMask5555(geno_word);
      uintptr_t geno_odd = PackWordToHalfwordMaskAAAA(geno_word);
      geno_word = genoarr[2 * widx + 1];
      geno_even |= S_CAST(uintptr_t, PackWordToHalfwordMask5555(geno_word)) << kBitsPerWordD2;
      geno_odd |= S_CAST(uintptr_t, PackWordToHalfwordMaskAAAA(geno_word)) << kBitsPerWordD2;
      const uintptr_t geno_even_masked = geno_even & mask_word;
      subset_intersect_ct += PopcountWord(mask_word);
      even_ct += PopcountWord(geno_even_masked);
      odd_ct += PopcountWord(geno_odd & mask_word);
      bothset_ct += PopcountWord(geno_odd & geno_even_masked);
    }
  }
  if (raw_sample_ctl2 % 2) {
    const uintptr_t mask_hw = subset1[fullword_ct] & subset2[fullword_ct];
    if (mask_hw) {
      const uintptr_t geno_word = genoarr[fullword_ct * 2];
      const uintptr_t mask_word = UnpackHalfwordToWord(mask_hw);
      const uintptr_t geno_word_shifted = geno_word >> 1;
      const uintptr_t geno_word_masked = geno_word & mask_word;
      subset_intersect_ct += Popcount01Word(mask_word);
      even_ct += Popcount01Word(geno_word_masked);
      odd_ct += Popcount01Word(geno_word_shifted & mask_word);
      bothset_ct += Popcount01Word(geno_word_masked & geno_word_shifted);
    }
  }
  genocounts[0] = subset_intersect_ct + bothset_ct - even_ct - odd_ct;
  genocounts[1] = even_ct - bothset_ct;
  genocounts[2] = odd_ct - bothset_ct;
  genocounts[3] = bothset_ct;
}

void GenovecInvertUnsafe(uint32_t sample_ct, uintptr_t* genovec) {
  // flips 0 to 2 and vice versa.
  // "unsafe" because trailing bits are not zeroed out.
  const uint32_t vec_ct = NypCtToVecCt(sample_ct);
  assert(IsVecAligned(genovec));
  const VecW not_m1 = VCONST_W(kMaskAAAA);
  VecW* vptr = R_CAST(VecW*, genovec);
  for (uint32_t vidx = 0; vidx != vec_ct; ++vidx) {
    VecW cur_vec = vptr[vidx];
    // flip high bit iff low bit is unset
    vptr[vidx] = cur_vec ^ vecw_and_notfirst(vecw_slli(cur_vec, 1), not_m1);
  }
}
#ifdef USE_SSE2
static_assert(kPglNypTransposeBatch == S_CAST(uint32_t, kNypsPerCacheline), "TransposeNypblock64() needs to be updated.");
#  ifdef CACHELINE64
#  else
#    ifndef CACHELINE128
#      error "CACHELINE64 or CACHELINE128 expected."
#    endif
// assumes USE_SHUFFLE8, !USE_AVX2
#  endif
#else  // !USE_SSE2
#  ifdef __LP64__
static_assert(kWordsPerVec == 1, "TransposeNypblock64() needs to be updated.");
#  else
static_assert(kWordsPerVec == 1, "TransposeNypblock32() needs to be updated.");
#  endif  // !__LP64__
#endif  // !USE_SSE2

void BiallelicDosage16Invert(uint32_t dosage_ct, uint16_t* dosage_main) {
  // replace each x with (32768 - x).
  // compiler is smart enough to vectorize this.
  for (uint32_t uii = 0; uii != dosage_ct; ++uii) {
    dosage_main[uii] = 32768 - dosage_main[uii];
  }
}
#if defined(USE_SSE2) && !defined(USE_AVX2)
void PackWordsToHalfwordsInvmatch(const uintptr_t* __restrict genoarr, uintptr_t inv_match_word, uint32_t inword_ct, uintptr_t* __restrict dst) {
  // In shuffle8 case, this takes ~30% less time than a
  // PackWordToHalfwordMask5555 loop.
  const uint32_t out_fullvec_ct = inword_ct / (kWordsPerVec * 2);
  const VecW xor_vec = vecw_set1(inv_match_word);
  const VecW m1 = VCONST_W(kMask5555);
#  ifdef USE_SHUFFLE8
  const VecW swap12 = vecw_setr8(
      0, 1, 4, 5, 2, 3, 6, 7,
      8, 9, 12, 13, 10, 11, 14, 15);
#  else
  const VecW m2 = VCONST_W(kMask3333);
#  endif
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW m8 = VCONST_W(kMask00FF);
  for (uintptr_t vidx = 0; vidx != out_fullvec_ct; ++vidx) {
    VecW vec_lo = vecw_loadu(&(genoarr[2 * kWordsPerVec * vidx])) ^ xor_vec;
    VecW vec_hi = vecw_loadu(&(genoarr[2 * kWordsPerVec * vidx + kWordsPerVec])) ^ xor_vec;
    vec_lo = vec_lo & vecw_srli(vec_lo, 1) & m1;
    vec_hi = vec_hi & vecw_srli(vec_hi, 1) & m1;
#  ifdef USE_SHUFFLE8
    vec_lo = (vec_lo | vecw_srli(vec_lo, 3)) & m4;
    vec_hi = (vec_hi | vecw_srli(vec_hi, 3)) & m4;
    vec_lo = vecw_shuffle8(swap12, vec_lo);
    vec_hi = vecw_shuffle8(swap12, vec_hi);
#  else
    vec_lo = (vec_lo | vecw_srli(vec_lo, 1)) & m2;
    vec_hi = (vec_hi | vecw_srli(vec_hi, 1)) & m2;
    vec_lo = (vec_lo | vecw_srli(vec_lo, 2)) & m4;
    vec_hi = (vec_hi | vecw_srli(vec_hi, 2)) & m4;
#  endif
    vec_lo = vec_lo | vecw_srli(vec_lo, 4);
    vec_hi = vec_hi | vecw_srli(vec_hi, 4);
    const VecW vec_packed = vecw_gather_even(vec_lo, vec_hi, m8);
    vecw_storeu(&(dst[kWordsPerVec * vidx]), vec_packed);
  }
  Halfword* dst_alias = DowncastWToHW(dst);
  uint32_t widx = RoundDownPow2(inword_ct, kWordsPerVec * 2);
  for (; widx != inword_ct; ++widx) {
    const uintptr_t cur_word = genoarr[widx] ^ inv_match_word;
    const Halfword hw = PackWordToHalfwordMask5555(cur_word & (cur_word >> 1));
    dst_alias[widx] = hw;
  }
}
void MaskWordsToHalfwordsInvmatch(const uintptr_t* __restrict genoarr, uintptr_t inv_match_word, uint32_t inword_ct, uintptr_t* src, uintptr_t* dst) {
  const uint32_t out_fullvec_ct = inword_ct / (kWordsPerVec * 2);
  const VecW xor_vec = vecw_set1(inv_match_word);
  const VecW m1 = VCONST_W(kMask5555);
#  ifdef USE_SHUFFLE8
  const VecW swap12 = vecw_setr8(
      0, 1, 4, 5, 2, 3, 6, 7,
      8, 9, 12, 13, 10, 11, 14, 15);
#  else
  const VecW m2 = VCONST_W(kMask3333);
#  endif
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW m8 = VCONST_W(kMask00FF);
  for (uintptr_t vidx = 0; vidx != out_fullvec_ct; ++vidx) {
    VecW vec_lo = vecw_loadu(&(genoarr[2 * kWordsPerVec * vidx])) ^ xor_vec;
    VecW vec_hi = vecw_loadu(&(genoarr[2 * kWordsPerVec * vidx + kWordsPerVec])) ^ xor_vec;
    VecW src_vec = vecw_loadu(&(src[kWordsPerVec * vidx]));
    vec_lo = vec_lo & vecw_srli(vec_lo, 1) & m1;
    vec_hi = vec_hi & vecw_srli(vec_hi, 1) & m1;
#  ifdef USE_SHUFFLE8
    vec_lo = (vec_lo | vecw_srli(vec_lo, 3)) & m4;
    vec_hi = (vec_hi | vecw_srli(vec_hi, 3)) & m4;
    vec_lo = vecw_shuffle8(swap12, vec_lo);
    vec_hi = vecw_shuffle8(swap12, vec_hi);
#  else
    vec_lo = (vec_lo | vecw_srli(vec_lo, 1)) & m2;
    vec_hi = (vec_hi | vecw_srli(vec_hi, 1)) & m2;
    vec_lo = (vec_lo | vecw_srli(vec_lo, 2)) & m4;
    vec_hi = (vec_hi | vecw_srli(vec_hi, 2)) & m4;
#  endif
    vec_lo = vec_lo | vecw_srli(vec_lo, 4);
    vec_hi = vec_hi | vecw_srli(vec_hi, 4);
    const VecW vec_packed = vecw_gather_even(vec_lo, vec_hi, m8);
    vecw_storeu(&(dst[kWordsPerVec * vidx]), src_vec & vec_packed);
  }
  Halfword* src_alias = DowncastWToHW(src);
  Halfword* dst_alias = DowncastWToHW(dst);
  uint32_t widx = RoundDownPow2(inword_ct, kWordsPerVec * 2);
  for (; widx != inword_ct; ++widx) {
    const uintptr_t cur_word = genoarr[widx] ^ inv_match_word;
    const Halfword hw = PackWordToHalfwordMask5555(cur_word & (cur_word >> 1));
    dst_alias[widx] = src_alias[widx] & hw;
  }
}
#else
void PackWordsToHalfwordsInvmatch(const uintptr_t* __restrict genoarr, uintptr_t inv_match_word, uint32_t inword_ct, uintptr_t* __restrict dst) {
  Halfword* dst_alias = DowncastWToHW(dst);
  for (uint32_t widx = 0; widx != inword_ct; ++widx) {
    const uintptr_t cur_word = genoarr[widx] ^ inv_match_word;
    const Halfword hw = PackWordToHalfwordMask5555(cur_word & (cur_word >> 1));
    dst_alias[widx] = hw;
  }
}
void MaskWordsToHalfwordsInvmatch(const uintptr_t* __restrict genoarr, uintptr_t inv_match_word, uint32_t inword_ct, uintptr_t* src, uintptr_t* dst) {
  Halfword* src_alias = DowncastWToHW(src);
  Halfword* dst_alias = DowncastWToHW(dst);
  for (uint32_t widx = 0; widx != inword_ct; ++widx) {
    const uintptr_t cur_word = genoarr[widx] ^ inv_match_word;
    const Halfword hw = PackWordToHalfwordMask5555(cur_word & (cur_word >> 1));
    dst_alias[widx] = src_alias[widx] & hw;
  }
}
#endif

void SparseToMissingness(const uintptr_t* __restrict raregeno, const uint32_t* difflist_sample_ids, uint32_t sample_ct, uint32_t difflist_common_geno, uint32_t difflist_len, uintptr_t* __restrict missingness) {
  if (difflist_common_geno != 3) {
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    ZeroWArr(sample_ctl, missingness);
    if (!difflist_len) {
      return;
    }
    const uint32_t raregeno_word_ct = NypCtToWordCt(difflist_len);
    for (uint32_t widx = 0; widx != raregeno_word_ct; ++widx) {
      const uintptr_t raregeno_word = raregeno[widx];
      uintptr_t raregeno_11 = raregeno_word & (raregeno_word >> 1) & kMask5555;
      if (raregeno_11) {
        const uint32_t* cur_difflist_sample_ids = &(difflist_sample_ids[widx * kBitsPerWordD2]);
        do {
          const uint32_t sample_idx_lowbits = ctzw(raregeno_11) / 2;
          const uint32_t cur_sample_id = cur_difflist_sample_ids[sample_idx_lowbits];
          SetBit(cur_sample_id, missingness);
          raregeno_11 &= raregeno_11 - 1;
        } while (raregeno_11);
      }
    }
  } else {
    SetAllBits(sample_ct, missingness);
    // Don't need to look at raregeno, all cases are nonmissing.
    for (uint32_t uii = 0; uii != difflist_len; ++uii) {
      const uint32_t cur_sample_id = difflist_sample_ids[uii];
      ClearBit(cur_sample_id, missingness);
    }
  }
}

#if defined(USE_SSE2) && !defined(USE_AVX2)
#else
#endif
#if defined(USE_SSE2) && !defined(USE_AVX2)
#else
#endif

#ifdef USE_SHUFFLE8
// ~70% less time than per-byte lookup.
#else
#endif

#ifndef NO_UNALIGNED
#else
// plink2_glm_logistic does not guarantee 8 byte table alignment.
#endif

// this might be important for genovec -> AlleleCode expansion
#ifdef USE_SSE2
void GenoarrLookup16x8bx2(const uintptr_t* genoarr, const void* table16x8bx2, uint32_t sample_ct, void* __restrict result) {
  const __m128i* table_alias = S_CAST(const __m128i*, table16x8bx2);
  unsigned char* result_biter = S_CAST(unsigned char*, result);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
        if (sample_ct % 2) {
          memcpy(result_biter, &(table_alias[geno_word & 3]), 8);
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != loop_len; ++uii) {
      const uintptr_t cur_2geno = geno_word & 15;
      _mm_storeu_si128(R_CAST(__m128i*, result_biter), table_alias[cur_2geno]);
      result_biter += 16;
      geno_word >>= 4;
    }
  }
}
#else // !USE_SSE2
void GenoarrLookup16x8bx2(const uintptr_t* genoarr, const void* table16x8bx2, uint32_t sample_ct, void* __restrict result) {
  const uint64_t* table_alias = S_CAST(const uint64_t*, table16x8bx2);
  unsigned char* result_biter = S_CAST(unsigned char*, result);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
        if (sample_ct % 2) {
          CopyToUnalignedU64(result_biter, &(table_alias[(geno_word & 3) * 2]));
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != loop_len; ++uii) {
      const uintptr_t cur_2geno = geno_word & 15;
      memcpy(result_biter, &(table_alias[cur_2geno * 2]), 16);
      result_biter += 16;
      geno_word >>= 4;
    }
  }
}
#endif
void InitLookup16x8bx2(void* table16x8bx2) {
  uint64_t* table_iter = S_CAST(uint64_t*, table16x8bx2);
  uint64_t vals[4];
  vals[0] = table_iter[0];
  table_iter[1] = vals[0];
  vals[1] = table_iter[2];
  table_iter[3] = vals[0];
  vals[2] = table_iter[4];
  table_iter[5] = vals[0];
  vals[3] = table_iter[6];
  table_iter[7] = vals[0];
  table_iter = &(table_iter[8]);
  for (uint32_t high_idx = 1; high_idx != 4; ++high_idx) {
    // bugfix (20 Jun 2018): cur_high needs to be a uint64_t, not a uint32_t
    const uint64_t cur_high = vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
}
#ifdef USE_SSE2
#else // !USE_SSE2
#endif
// bits 0..3: two genotypes
// bits 4..5: two (phasepresent | sex_male) bits
// bits 1,3: unpacked phaseinfo xor
#ifdef USE_SSE2
#else // !USE_SSE2
#endif
double MultiallelicDiploidMinimac3R2(const uint64_t* __restrict sums, const uint64_t* __restrict hap_ssqs_x2, uint32_t nm_sample_ct, uint32_t allele_ct, uint32_t extra_phased_het_ct) {
  // sums[k] == sum_i [left_dosage_{ik} + right_dosage_{ik}]
  // hap_ssqs_x2[k] ==
  //   2 * sum_i [(left_dosage_{ik})^2 + (right_dosage_{ik})^2]
  //   This may be odd, since it's computed as
  //     (left + right)^2 + (left - right)^2
  //   and the parities of the two integers can be different.
  // For phased hardcalls, it is fine for the hap_ssqs_x2[k] values to
  // correspond to unphased hardcalls iff extra_phased_het_ct is the number of
  // phased-hets that weren't accounted for in hap_ssqs_x2[]; this makes it
  // straightforward for GetMultiallelicCountsAndDosage16s to stick to the
  // custom internal multiallelic-count functions.
  if (!nm_sample_ct) {
    return (0.0 / 0.0);
  }
  // Allelic phased-dosages are on a (k-1)-dimensional simplex; embed this in
  // R^k as the (1, 0, ..., 0), (0, 1, ..., 0), ..., (0, 0, ..., 1) polytope.
  // Define
  //   m_k := (1/2n) * sum_i [left_dosage_{ik} + right_dosage_{ik}]
  // Minimac3-r2 is defined as empirical phased-dosage variance divided by
  // expected-under-allele-frequencies variance.
  // Expected sum-of-squared-Euclidean-distance with perfect imputation is
  //   2n("1"^2 - sum_k ((m_k)^2))
  // and observed sum-of-squared-distance is
  //   sum_k (sum_i [(left_dosage_{ik})^2 + (right_dosage_{ik})^2] -
  //          2n((m_k)^2))

  // ssq_sum_x2 can be as large as 2^31 * nm_sample_ct; meansq_sum can cancel
  // as little as (1 / allele_ct) of it
  if (nm_sample_ct < 92682) {
    uint64_t ssq_sum_x2 = extra_phased_het_ct * 0x20000000LLU;
    uint64_t meansq_sum = 0;
    for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
      const uint64_t cur_allele_dosage = sums[allele_idx];
      ssq_sum_x2 += hap_ssqs_x2[allele_idx];
      // cur_allele_dosage == 2n * m_k
      // -> meansq_sum becomes 2n * sum_k [2n((m_k)^2)]
      meansq_sum += cur_allele_dosage * cur_allele_dosage;
    }
    const uint64_t observed_variance_times_2n = ssq_sum_x2 * nm_sample_ct - meansq_sum;
    // "1"^2 -> 16384^2 in our units.  So 2n * 2n * "1"^2 is equal to
    //   n * n * 16384^2 * 4.
    const uint64_t expected_variance_times_2n = nm_sample_ct * 0x40000000LLU * nm_sample_ct - meansq_sum;
    // mach_r2 == 1 test cases:
    // - AA, AB, BB: 1, 4, 4
    //   sums[0] = 6 * 2^14
    //   sums[1] = 12 * 2^14
    //   ssqs[0] = 8 * 2^28
    //   ssqs[1] = 20 * 2^28
    //   ssq_sum = (8 + 20) * 2^28
    //   meansq_sum = (6 * 6 + 12 * 12) * 2^28
    //   observed_variance = 28 * 9 * 2^28 - 180 * 2^28
    //   expected_variance = (9 * 9 * 4 * 2^28 - 180 * 2^28) / 2
    // - AA, AB, BB, AC, BC, CC: 1, 4, 4, 6, 12, 9
    //   sums[0] = 12 * 2^14
    //   sums[1] = 24 * 2^14
    //   sums[2] = 36 * 2^14
    //   ssqs[0] = 14 * 2^28
    //   ssqs[1] = 32 * 2^28
    //   ssqs[2] = 54 * 2^28
    //   ssq_sum = (14 + 32 + 54) * 2^28
    //   meansq_sum = (12 * 12 + 24 * 24 + 36 * 36) * 2^28
    //   observed_variance = 100 * 36 * 2^28 - 56 * 36 * 2^28
    //   expected_variance = (36 * 36 * 4 * 2^28 - 56 * 36 * 2^28) / 2
    return S_CAST(double, observed_variance_times_2n) / S_CAST(double, expected_variance_times_2n);
  }
  uint64_t ssq_sum_x2 = extra_phased_het_ct * 0x20000000LLU;
  uint64_t meansq_sum_lo = 0;
  uint64_t meansq_sum_hi = 0;
  for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
    const uint64_t cur_allele_dosage = sums[allele_idx];
    const uint64_t cur_ssq_x2 = hap_ssqs_x2[allele_idx];
    ssq_sum_x2 += cur_ssq_x2;
    uint64_t incr_hi;
    uint64_t incr_lo = multiply64to128(cur_allele_dosage, cur_allele_dosage, &incr_hi);
    meansq_sum_lo += incr_lo;
    meansq_sum_hi += incr_hi + (meansq_sum_lo < incr_lo);
  }
  uint64_t prod_hi;
  uint64_t prod_lo = multiply64to128(ssq_sum_x2, nm_sample_ct, &prod_hi);
  const uint64_t observed_variance_times_2n_lo = prod_lo - meansq_sum_lo;
  const uint64_t observed_variance_times_2n_hi = prod_hi - meansq_sum_hi - (prod_lo < meansq_sum_lo);
  const double observed_variance_times_2n = u127tod(observed_variance_times_2n_hi, observed_variance_times_2n_lo);

  prod_lo = multiply64to128(nm_sample_ct * 0x40000000LLU, nm_sample_ct, &prod_hi);
  const uint64_t expected_variance_times_2n_lo = prod_lo - meansq_sum_lo;
  const uint64_t expected_variance_times_2n_hi = prod_hi - meansq_sum_hi - (prod_lo < meansq_sum_lo);
  const double expected_variance_times_2n = u127tod(expected_variance_times_2n_hi, expected_variance_times_2n_lo);
  return observed_variance_times_2n / expected_variance_times_2n;
}

void PgrDifflistToGenovecUnsafe(const uintptr_t* __restrict raregeno, const uint32_t* difflist_sample_ids, uintptr_t difflist_common_geno, uint32_t sample_ct, uint32_t difflist_len, uintptr_t* __restrict genovec) {
  // Ok for trailing bits of raregeno to be nonzero.  Does not zero out
  // trailing bits of genovec.
  const uint32_t vec_ct = NypCtToVecCt(sample_ct);
  // could just memset up to word boundary; this should be a bit more
  // vector-instruction-friendly, though
  vecset(genovec, difflist_common_geno * kMask5555, vec_ct);
  const uintptr_t* raregeno_incr = raregeno;
  uint32_t difflist_idx = 0;
  uint32_t difflist_idx_stop = kBitsPerWordD2;
  if (!difflist_common_geno) {
    // faster inner loop since there's no existing value to mask out
    // todo: check if this should just be deleted since the code bloat causes
    // too many more cache misses
    for (; ; difflist_idx_stop += kBitsPerWordD2) {
      if (difflist_idx_stop > difflist_len) {
        if (difflist_idx == difflist_len) {
          return;
        }
        difflist_idx_stop = difflist_len;
      }
      uintptr_t raregeno_word = *raregeno_incr++;
      for (; difflist_idx != difflist_idx_stop; ++difflist_idx) {
        const uint32_t cur_sample_idx = difflist_sample_ids[difflist_idx];
        genovec[cur_sample_idx / kBitsPerWordD2] |= (raregeno_word & 3) << (2 * (cur_sample_idx % kBitsPerWordD2));
        raregeno_word >>= 2;
      }
    }
  }
  for (; ; difflist_idx_stop += kBitsPerWordD2) {
    if (difflist_idx_stop > difflist_len) {
      if (difflist_idx == difflist_len) {
        return;
      }
      difflist_idx_stop = difflist_len;
    }
    uintptr_t raregeno_word = *raregeno_incr++;
    for (; difflist_idx != difflist_idx_stop; ++difflist_idx) {
      const uint32_t cur_sample_idx = difflist_sample_ids[difflist_idx];
      AssignNyparrEntry(cur_sample_idx, raregeno_word & 3, genovec);
      raregeno_word >>= 2;
    }
  }
}

const uint16_t kHcToAlleleCodes[1024] = QUAD_TABLE256(0, 0x100, 0x101, 0xffff);

static_assert(sizeof(AlleleCode) == 1, "PglMultiallelicSparseToDenseMiss() needs to be updated.");
// Ok for nybble_vvec to be unaligned.
uint32_t CountNybbleVec(const unsigned char* nybble_vvec_biter, uintptr_t nybble_word, uint32_t vec_ct) {
  const VecW m0 = vecw_setzero();
  const VecW alld15 = VCONST_W(kMask1111);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW xor_vvec = vecw_set1(nybble_word);
  VecW prev_sad_result = vecw_setzero();
  VecW acc = vecw_setzero();
  uintptr_t cur_incr = 15;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 15) {
      if (!vec_ct) {
        acc = acc + prev_sad_result;
        return HsumW(acc);
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc = vecw_setzero();
    const unsigned char* nybble_vvec_stop = &(nybble_vvec_biter[cur_incr * kBytesPerVec]);
    do {
      VecW loader = vecw_loadu(nybble_vvec_biter) ^ xor_vvec;
      nybble_vvec_biter += kBytesPerVec;
      // DetectAllZeroNybbles() followed by right-shift-3 is the same number of
      // operations, can see if that's any faster in practice
      loader = vecw_srli(loader, 1) | loader;
      loader = vecw_srli(loader, 2) | loader;
      inner_acc = inner_acc + vecw_and_notfirst(loader, alld15);
    } while (nybble_vvec_biter < nybble_vvec_stop);
    inner_acc = (inner_acc & m4) + (vecw_srli(inner_acc, 4) & m4);
    acc = acc + prev_sad_result;
    prev_sad_result = vecw_bytesum(inner_acc, m0);
  }
}

uint32_t CountNybble(const void* nybblearr, uintptr_t nybble_word, uintptr_t nybble_ct) {
  const unsigned char* nybblearr_uc = S_CAST(const unsigned char*, nybblearr);
  const uint32_t fullword_ct = nybble_ct / kBitsPerWordD4;
  uint32_t tot = CountNybbleVec(nybblearr_uc, nybble_word, fullword_ct / kWordsPerVec);
#ifdef __LP64__
  for (uint32_t word_idx = RoundDownPow2(fullword_ct, kWordsPerVec); word_idx != fullword_ct; ++word_idx) {
    uintptr_t cur_word;
    CopyFromUnalignedOffsetW(&cur_word, nybblearr_uc, word_idx);
    cur_word ^= nybble_word;
    cur_word = cur_word | (cur_word >> 1);
    cur_word = cur_word | (cur_word >> 2);
    tot += Popcount0001Word((~cur_word) & kMask1111);
  }
#endif
  const uint32_t trailing_nybble_ct = nybble_ct % kBitsPerWordD4;
  if (trailing_nybble_ct) {
    const uint32_t trailing_byte_ct = DivUp(trailing_nybble_ct, (CHAR_BIT / 4));
    uintptr_t cur_word = SubwordLoad(&(nybblearr_uc[fullword_ct * kBytesPerWord]), trailing_byte_ct) ^ nybble_word;
    cur_word = cur_word | (cur_word >> 1);
    cur_word = cur_word | (cur_word >> 2);
    cur_word = bzhi((~cur_word) & kMask1111, trailing_nybble_ct * 4);
#if defined(USE_SSE42) || !defined(__LP64__)
    tot += Popcount0001Word(cur_word);
#else
    // minor optimization, can't overflow
    tot += (cur_word * kMask1111) >> 60;
#endif
  }
  return tot;
}

uint64_t PglHeaderBaseEndOffset(uint32_t variant_ct, uintptr_t vrec_len_byte_ct, uint32_t phase_or_dosage_present, uint32_t explicit_nonref_flags) {
  const uint32_t vblock_ct = DivUp(variant_ct, kPglVblockSize);
  uint64_t offset = 12 + vblock_ct * sizeof(int64_t) + variant_ct * vrec_len_byte_ct;
  if (phase_or_dosage_present) {
    // 8-bit vrtypes
    offset += variant_ct;
  } else {
    // 4-bit vrtypes
    offset += DivUp(variant_ct, 2);
  }
  if (explicit_nonref_flags) {
    offset += DivUp(variant_ct, CHAR_BIT);
  }
  return offset;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
