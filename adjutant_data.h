#ifndef __ADJUTANT_DATA_H__
#define __ADJUTANT_DATA_H__

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

#include "include/pgenlib_misc.h"
#include "include/pgenlib_read.h"
#include "include/plink2_base.h"
#include "include/plink2_text.h"
#include "adjutant_cmdline.h"
#include "adjutant_common.h"

#ifdef __cplusplus
namespace plink2 {
#endif

FLAGSET_DEF_START()
  kfMake0,
  kfMakeBed = (1 << 0),
  kfMakeBim = (1 << 1),
  kfMakeFam = (1 << 2),
  kfMakePgen = (1 << 3),
  kfMakePvar = (1 << 4),
  kfMakePsam = (1 << 5),
  kfMakeBimZs = (1 << 6),
  kfMakePlink2MSplitBase = (1 << 7), // three bits for multiallelic mode
  kfMakePlink2MSplitAll = kfMakePlink2MSplitBase,
  kfMakePlink2MSplitSnps = 2 * kfMakePlink2MSplitBase,
  kfMakePlink2MJoin = (1 << 9),
  kfMakePlink2MJoinBoth = kfMakePlink2MJoin,
  kfMakePlink2MJoinSnps = kfMakePlink2MJoin + kfMakePlink2MSplitBase,
  kfMakePlink2MJoinAny = kfMakePlink2MJoin + 2 * kfMakePlink2MSplitBase,
  // don't support e.g. '+indels' for now due to lack of standardization re:
  // handling of MNP/'other' classes
  kfMakePlink2TrimAlts = (1 << 10),
  kfMakePlink2MMask = kfMakePlink2TrimAlts - kfMakePlink2MSplitBase,
  kfMakePlink2VaridSemicolon = (1 << 11),
  kfMakePlink2VaridDup = (1 << 12),

  kfMakePlink2SetInvalidHaploidMissing = (1 << 13),
  kfMakePlink2SetInvalidHaploidMissingKeepDosage = (1 << 14),
  kfMakePlink2SetMixedMtMissing = (1 << 15),
  kfMakePlink2SetMixedMtMissingKeepDosage = (1 << 16),
  kfMakePlink2SetMeMissing = (1 << 17),
  kfMakePlink2FillMissingWithRef = (1 << 18),
  kfMakePgenFormatBase = (1 << 19), // two bits
  kfMakePgenErasePhase = (1 << 21),
  kfMakePgenEraseDosage = (1 << 22),
  kfMakePgenFillMissingFromDosage = (1 << 23),
  kfMakePgenWriterVer = (1 << 24)
FLAGSET_DEF_END(MakePlink2Flags);

FLAGSET_DEF_START()
  kfImport0,
  kfImportKeepAutoconv = (1 << 0),
  kfImportKeepAutoconvVzs = (1 << 1),
  kfImportDoubleId = (1 << 2),
  kfImportVcfRefNMissing = (1 << 3),
  kfImportLaxChrX = (1 << 4),
  kfImportPolyploidMissing = (1 << 5),
  kfImportPolyploidExplicitError = (1 << 6),
  kfImportVcfAllowNoNonvar = (1 << 7),
  kfImportLaxBgen = (1 << 8),
  kfImportEigNohash = (1 << 9)
FLAGSET_DEF_END(ImportFlags);
CONSTI32(kMaxInfoKeySlen, kMaxIdSlen);
#define MAX_INFO_KEY_SLEN_STR MAX_ID_SLEN_STR

// We only need to distinguish between the following INFO-value-type cases:
// Number=0 (flag), Number=<positive integer>, Number=., Number=A, and
// Number=R.  (Number=G is being treated as Number=., since otherwise chrX is a
// nightmare; see https://github.com/samtools/hts-specs/issues/272 for some
// discussion.)  We use negative numbers to represent the last 3 cases in
// InfoVtype.
CONSTI32(kInfoVtypeUnknown, -1);
CONSTI32(kInfoVtypeA, -2);
CONSTI32(kInfoVtypeR, -3);
// Used in plink2_data (though not plink2_merge) for INFO/PR flag.
CONSTI32(kInfoVtypeSkip, -4);
// Main fixed data structure when splitting/joining INFO is a hashmap of keys.
// Behavior when splitting:
// - Field order in the original variant is retained.
// - Number >= 0 and Number=. don't require any special handling, just copy the
//   entire key=value pair (or lone key, in the Flag case).
// - Number=A and Number=R require splitting the value on ',' and verifying the
//   comma count is correct, but is otherwise straightforward since alleles
//   can't be permuted.
// - Number=G is being treated as Number=.
// When joining:
// - Field order is determined by header line order.
// - Number=./G and Number>0 just require a buffer of size ~info_reload_slen,
//   and a boolean indicating whether no mismatch has been found.
// - Number=0 (Flag) requires a single boolean, we perform an or operation.
// - Number=A/R are the messy ones: we need to have enough space for
//   max_write_allele_ct (or that minus 1) comma-separated values in the =A and
//   =R cases.
// uint32_t ChrLenLbound(const ChrInfo* cip, const uint32_t* variant_bps, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uint32_t* new_variant_idx_to_old, uint32_t chr_fo_idx, uint32_t max_allele_slen, UnsortedVar vpos_sortstatus);




// fileformat, fileDate, source
// assumes adequate buffer space











PglErr LoadAlleleAndGenoCounts(const uintptr_t* sample_include, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uintptr_t* allele_idx_offsets, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t founder_ct, uint32_t male_ct, uint32_t nosex_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t first_hap_uidx, uint32_t is_minimac3_r2, uint32_t y_nosex_missing_stats, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, uintptr_t* allele_presents, uint64_t* allele_ddosages, uint64_t* founder_allele_ddosages, uint32_t* variant_missing_hc_cts, uint32_t* variant_missing_dosage_cts, uint32_t* variant_hethap_cts, STD_ARRAY_PTR_DECL(uint32_t, 3, raw_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, x_male_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, x_nosex_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts), double* imp_r2_vals);








#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __ADJUTANT_DATA_H__
