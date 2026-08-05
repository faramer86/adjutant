// Variant harmonisation for --pgs-model.
//
// Reimplements bigsnpr::snp_match() allele logic as a pure function.
// See adjutant_pgs_match.h for the full contract.

#include "adjutant_pgs_match.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <cmath>

#include "adjutant_cmdline.h"
#include "adjutant_decompress.h"
#include "include/plink2_htable.h"
#include "include/plink2_string.h"
#include "include/plink2_text.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// -------------------------------------------------------------------
// Strand complement of a single character.
//
// Returns the complement for A/C/T/G (case-sensitive).
// Returns '\0' for anything else, so flipped candidates never match.
// -------------------------------------------------------------------
static char CompAlleleChar(char c) {
  switch (c) {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'C': return 'G';
    case 'G': return 'C';
    default:  return '\0';
  }
}

static void FlipSingleAllele(const char* allele, char out[2]) {
  if (allele[0] != '\0' && allele[1] == '\0') {
    out[0] = CompAlleleChar(allele[0]);
    out[1] = '\0';
  } else {
    out[0] = '\0';
    out[1] = '\0';
  }
}

static uint32_t IsAmbiguous(const char* a0, const char* a1) {
  if (a0[0] == '\0' || a0[1] != '\0') return 0;
  if (a1[0] == '\0' || a1[1] != '\0') return 0;
  const char c0 = a0[0];
  const char c1 = a1[0];
  return (c0 == 'A' && c1 == 'T') ||
         (c0 == 'T' && c1 == 'A') ||
         (c0 == 'C' && c1 == 'G') ||
         (c0 == 'G' && c1 == 'C');
}

// Do {a0, a1} and {b0, b1} name the same unordered allele pair?
static uint32_t AllelePairSetEqual(const char* a0, const char* a1,
                                   const char* b0, const char* b1) {
  return (strcmp(a0, b0) == 0 && strcmp(a1, b1) == 0) ||
         (strcmp(a0, b1) == 0 && strcmp(a1, b0) == 0);
}

// Is `f` too close to 0.5 for the two strand hypotheses to be told apart?
//
// The slack is not decoration.  ScanadvDouble, which parses every
// frequency here, says outright that it is "ok to lose last ~2 bits of
// precision", and it uses them: "0.6" comes back one ulp high, so 1 - it
// lands just under a threshold of 0.4 and resolves, while the identical
// MAF spelled "0.4" is refused.  1e-15 is 18 ulp at the default
// threshold of 0.4 -- several times the parser's documented ~2 bits, and
// nine orders of magnitude below the 1e-6 resolution of a frequency
// estimated from even a million samples.  So the two spellings of one
// MAF agree, and what error remains goes toward refusing rather than
// toward a guessed sign.
static uint32_t FreqTooCloseToHalf(double f, double maf_threshold) {
  const double maf = (f <= 0.5) ? f : (1.0 - f);
  return maf >= maf_threshold - 1e-15;
}

// -------------------------------------------------------------------
// Decide a palindromic variant from allele frequencies.  See the
// contract in adjutant_pgs_match.h.
// -------------------------------------------------------------------
static PgsMatchStatus ResolveAmbiguous(const char* src_a0, const char* src_a1,
                                       const char* tgt_a0, const char* tgt_a1,
                                       const PgsAmbigResolve* r,
                                       uint32_t* out_reversed) {
  // The pair must actually be the same pair before frequencies get a
  // vote.  A source A/T against a target C/G is palindromic on both
  // sides yet shares no allele.  For a palindrome the strand complement
  // of the pair is that same pair, so set equality is the whole test.
  if (!AllelePairSetEqual(src_a0, src_a1, tgt_a0, tgt_a1)) {
    return kPgsMatchNoAllele;
  }
  // NaN fails every comparison, so this also rejects a missing or
  // unparseable frequency on either side.
  if (!(r->src_freq >= 0.0 && r->src_freq <= 1.0)) {
    return kPgsMatchAmbiguousUnresolved;
  }
  if (!(r->tgt_a1_freq >= 0.0 && r->tgt_a1_freq <= 1.0)) {
    return kPgsMatchAmbiguousUnresolved;
  }
  // The threshold applies to BOTH frequencies.  What separates the two
  // hypotheses is |2f - 1| on each side, so a target frequency near 0.5
  // is exactly as uninformative as a source frequency near 0.5 --
  // guarding only the source would let a target of 0.5 decide the sign
  // on no evidence at all.
  if (FreqTooCloseToHalf(r->src_freq, r->maf_threshold) ||
      FreqTooCloseToHalf(r->tgt_a1_freq, r->maf_threshold)) {
    return kPgsMatchAmbiguousUnresolved;
  }
  const double d_same = fabs(r->src_freq - r->tgt_a1_freq);
  const double d_opp  = fabs(r->src_freq - (1.0 - r->tgt_a1_freq));
  // d_same == d_opp only when src_freq or tgt_a1_freq is exactly 0.5,
  // and a MAF of 0.5 is at or above every permitted threshold -- the CLI
  // caps --resolve-ambiguous at 0.5 -- so both were refused above and the
  // tie cannot reach here.
  if (d_same < d_opp) {
    return kPgsMatchResolvedExact;
  }
  // The source effect allele is the target's a0, so the beta has to be
  // re-expressed w.r.t. a1.
  *out_reversed = 1;
  return kPgsMatchResolvedReversed;
}

PgsMatchStatus PgsMatchOne(const char* src_a0, const char* src_a1,
                           const char* tgt_a0, const char* tgt_a1,
                           uint32_t allow_ambiguous,
                           const PgsAmbigResolve* resolve,
                           uint32_t* out_reversed) {
  *out_reversed = 0;

  if (IsAmbiguous(src_a0, src_a1)) {
    if (resolve) {
      return ResolveAmbiguous(src_a0, src_a1, tgt_a0, tgt_a1,
                              resolve, out_reversed);
    }
    if (!allow_ambiguous) {
      return kPgsMatchAmbiguous;
    }
  }

  char f_a0[2], f_a1[2];
  FlipSingleAllele(src_a0, f_a0);
  FlipSingleAllele(src_a1, f_a1);

  // Candidate 1: exact
  if (strcmp(src_a0, tgt_a0) == 0 && strcmp(src_a1, tgt_a1) == 0) {
    return kPgsMatchExact;
  }
  // Candidate 2: reversed -- caller negates beta
  if (strcmp(src_a1, tgt_a0) == 0 && strcmp(src_a0, tgt_a1) == 0) {
    *out_reversed = 1;
    return kPgsMatchReversed;
  }
  // Candidate 3: strand-flipped
  if (f_a0[0] != '\0' && f_a1[0] != '\0' &&
      strcmp(f_a0, tgt_a0) == 0 && strcmp(f_a1, tgt_a1) == 0) {
    return kPgsMatchFlipped;
  }
  // Candidate 4: strand-flipped + reversed
  if (f_a0[0] != '\0' && f_a1[0] != '\0' &&
      strcmp(f_a1, tgt_a0) == 0 && strcmp(f_a0, tgt_a1) == 0) {
    *out_reversed = 1;
    return kPgsMatchFlippedReversed;
  }
  return kPgsMatchNoAllele;
}

// -------------------------------------------------------------------
// Chromosome normaliser (not thread-safe; single-threaded use only).
// Strips "chr" prefix, maps M -> MT.
// Returns pointer to static buffer, or nullptr if empty.
// -------------------------------------------------------------------
static const char* NormChr(const char* s) {
  if (!s || !s[0]) return nullptr;
  if (s[0] == 'c' && s[1] == 'h' && s[2] == 'r') s += 3;
  if (!s[0]) return nullptr;
  static char buf[8];
  if (s[0] == 'M' && s[1] == '\0') {
    buf[0] = 'M'; buf[1] = 'T'; buf[2] = '\0';
    return buf;
  }
  size_t n = strlen(s);
  if (n > 7) n = 7;
  memcpy(buf, s, n);
  buf[n] = '\0';
  return buf;
}

// True when this chr/pos text pair falls inside the excluded MHC interval.
//
// This route has text tokens rather than a ChrInfo, so it establishes "is
// chromosome 6" by normalised name; the interval itself comes from the one
// shared table in adjutant_mhc.cc, same as the scoring route.
//
// An unparseable position is not excluded.  Deciding it is malformed is the
// matcher's job, not this predicate's, and it already reports such a row as
// kPgsMatchNoPos -- silently swallowing it here would hide a broken input
// behind a flag the user set for an unrelated reason.
static uint32_t MhcExcludesChrPos(MhcRegion region, const char* chr,
                                  const char* pos_str) {
  if (region == kMhcNone) return 0;
  const char* nchr = NormChr(chr);
  if (!nchr || strcmp(nchr, "6")) return 0;
  uint32_t bp;
  if (ScanUintDefcap(pos_str, &bp)) return 0;
  return MhcExcludesBp(region, bp);
}

// Build hash key "CHR\tPOS" into caller buffer.  Returns key length or 0.
static uint32_t BuildChrPosKey(const char* chr, const char* pos_str,
                                char* key_buf, uint32_t key_cap) {
  const char* nchr = NormChr(chr);
  if (!nchr) return 0;
  uint32_t chr_len = (uint32_t)strlen(nchr);
  uint32_t pos_len = (uint32_t)strlen(pos_str);
  uint32_t total = chr_len + 1 + pos_len;
  if (total >= key_cap) return 0;
  memcpy(key_buf, nchr, chr_len);
  key_buf[chr_len] = '\t';
  memcpy(key_buf + chr_len + 1, pos_str, pos_len);
  key_buf[total] = '\0';
  return total;
}

// -------------------------------------------------------------------
// BIM record.
// -------------------------------------------------------------------
typedef struct BimEntryStruct {
  char* varid;
  char* chr;  // column 1
  char* pos;  // column 4, kept as text -- it is only ever printed
  char* a1;   // column 5 (effect allele in bim convention)
  char* a0;   // column 6 (other allele)
  // Inside the --exclude-*mhc-* interval.  Flagged rather than dropped: the
  // ID-matching route has no position on the source side, so filtering there
  // is impossible, and dropping from the index alone would report an excluded
  // variant as NO_MATCH -- "absent from your target" rather than "you asked me
  // to drop this".  Matching against a flagged row drops it honestly, in both
  // modes, and lets the source-side filter go away entirely.
  unsigned char mhc_excluded;
} BimEntry;

// -------------------------------------------------------------------
// Source variant record (one per row of the weight file).
// -------------------------------------------------------------------
typedef struct SrcVariantStruct {
  char*  src_id;      // ID mode only; nullptr under chr:pos matching
  char*  chr;         // chr:pos mode only; nullptr under ID matching
  char*  pos_str;
  char*  src_a1;      // effect allele (source col_effect)
  char*  src_a0;      // other allele  (source col_other)
  double src_beta;
  double src_freq;    // effect-allele freq in the source; NaN if absent
  PgsMatchStatus status;
  char*  bim_id;      // nullptr if unmatched
  char*  bim_chr;     // the matched .bim row's coordinates, for the report
  char*  bim_pos;
  char*  out_a1;      // nullptr if unmatched
  char*  out_a2;      // the target's other allele, for .matched
  double out_beta;
  // Row of the matched .bim entry, UINT32_MAX when unmatched.  --pgs-model
  // never needs it, but merging several models does: it keys the union and
  // orders it by .bim position without comparing ID strings.
  uint32_t bim_idx;
} SrcVariant;

// The key a source variant joins on, in whichever mode is active.  One helper
// so the lookup and the dedup pass cannot drift apart -- they must agree, or a
// variant could match a .bim row and then fail to be recognised as its own
// duplicate.
static uint32_t SrcKey(const SrcVariant* sv, uint32_t match_by_id,
                       char* key_buf, uint32_t key_cap) {
  if (!match_by_id) {
    return BuildChrPosKey(sv->chr, sv->pos_str, key_buf, key_cap);
  }
  const uint32_t slen = (uint32_t)strlen(sv->src_id);
  if (slen >= key_cap) return 0;
  memcpy(key_buf, sv->src_id, slen + 1);
  return slen;
}

// A variant is written to .matched exactly when its status is one of these.
static uint32_t PgsMatchIsMatched(PgsMatchStatus s) {
  return (s == kPgsMatchExact) || (s == kPgsMatchReversed) ||
         (s == kPgsMatchFlipped) || (s == kPgsMatchFlippedReversed) ||
         (s == kPgsMatchResolvedExact) || (s == kPgsMatchResolvedReversed);
}

static char* DupStr(const char* s, uint32_t len) {
  char* p = (char*)malloc(len + 1);
  if (p) { memcpy(p, s, len); p[len] = '\0'; }
  return p;
}

// Longest variant ID either side will key on.  Shared so that an over-long ID
// truncates identically in the .bim and in the weight file, and therefore still
// joins rather than silently failing to.
enum { kPgsMaxIdBlen = 256 };

// Bump allocator for a model's per-variant strings.
//
// Each SrcVariant carries six of them -- chr, pos_str, src_a1, src_a0, bim_id,
// out_a1 -- and a strdup apiece meant a 600k-variant model made 3.6M
// allocations of one to twelve characters, every one paying 24-32 bytes of
// allocator bookkeeping to hold them.  That was 227 bytes per source variant
// measured, about 11x the size of the input file, and freeing them left the
// heap too fragmented for the next model to reuse: peak memory climbed with
// every model in a --pgs-models-list run instead of staying flat.
//
// Chunked and never realloc'd, because SrcVariant points INTO these blocks and
// the backing store must not move.  Releasing a model is then a few free()s of
// megabyte blocks, which the allocator can actually hand back.
typedef struct StrArenaBlockStruct {
  struct StrArenaBlockStruct* next;
  size_t used;
  size_t cap;
} StrArenaBlock;

typedef struct StrArenaStruct {
  StrArenaBlock* head;   // every block, in allocation order
  StrArenaBlock* cur;    // the one being filled
} StrArena;

static const size_t kStrArenaBlockBytes = 1U << 20;

// Data lives immediately after the header, so a block is one allocation
// without relying on flexible array members.
static char* StrArenaData(StrArenaBlock* b) {
  return (char*)(b + 1);
}

static char* StrArenaDup(StrArena* a, const char* s, uint32_t len) {
  const size_t need = (size_t)len + 1;
  StrArenaBlock* b = a->cur;
  while (b && (b->cap - b->used < need)) b = b->next;   // reuse a spent block
  if (!b) {
    size_t cap = kStrArenaBlockBytes;
    if (cap < need) cap = need;
    b = (StrArenaBlock*)malloc(sizeof(StrArenaBlock) + cap);
    if (!b) return nullptr;
    b->next = nullptr;
    b->used = 0;
    b->cap  = cap;
    if (a->cur) {
      StrArenaBlock* tail = a->cur;
      while (tail->next) tail = tail->next;
      tail->next = b;
    } else {
      a->head = b;
    }
  }
  a->cur = b;
  char* dst = StrArenaData(b) + b->used;
  memcpy(dst, s, len);
  dst[len] = '\0';
  b->used += need;
  return dst;
}

// Hand every block back to the arena without returning it to the allocator.
// Models are processed one at a time, so the next one refills the same memory
// instead of making the heap churn 100k-plus allocations per model.
static void StrArenaReset(StrArena* a) {
  for (StrArenaBlock* b = a->head; b; b = b->next) b->used = 0;
  a->cur = a->head;
}

static void StrArenaFree(StrArena* a) {
  StrArenaBlock* b = a->head;
  while (b) {
    StrArenaBlock* nx = b->next;
    free(b);
    b = nx;
  }
  a->head = nullptr;
  a->cur  = nullptr;
}

// -------------------------------------------------------------------
// Flat-string store: appends NUL-terminated strings to a growing heap
// buffer and returns the pointer to the stored copy.
// -------------------------------------------------------------------
typedef struct FlatStrStoreStruct {
  char*    buf;
  uint32_t used;
  uint32_t cap;
  const char** ptrs;   // parallel array of pointers into buf
  uint32_t ptr_ct;
  uint32_t ptr_cap;
} FlatStrStore;

static FlatStrStore FlatStrStoreInit(uint32_t init_buf, uint32_t init_ptrs) {
  FlatStrStore st;
  st.buf     = (char*)malloc(init_buf);
  st.used    = 0;
  st.cap     = st.buf ? init_buf : 0;
  st.ptrs    = (const char**)malloc(init_ptrs * sizeof(const char*));
  st.ptr_ct  = 0;
  st.ptr_cap = st.ptrs ? init_ptrs : 0;
  return st;
}

// Returns the stored pointer, or nullptr on OOM.
static const char* FlatStrStoreAdd(FlatStrStore* st, const char* s, uint32_t slen) {
  uint32_t needed = st->used + slen + 1;
  if (needed > st->cap) {
    uint32_t new_cap = st->cap * 2;
    if (new_cap < needed) new_cap = needed * 2;
    char* nb = (char*)realloc(st->buf, new_cap);
    if (!nb) return nullptr;
    if (nb != st->buf) {
      ptrdiff_t delta = nb - st->buf;
      for (uint32_t i = 0; i < st->ptr_ct; ++i) st->ptrs[i] += delta;
    }
    st->buf = nb;
    st->cap = new_cap;
  }
  char* dst = st->buf + st->used;
  memcpy(dst, s, slen);
  dst[slen] = '\0';
  st->used += slen + 1;

  // Grow ptrs array.
  if (st->ptr_ct == st->ptr_cap) {
    uint32_t new_cap = st->ptr_cap * 2;
    if (!new_cap) new_cap = 64;
    const char** np = (const char**)realloc(st->ptrs, new_cap * sizeof(const char*));
    if (!np) return nullptr;
    st->ptrs = np;
    st->ptr_cap = new_cap;
  }
  st->ptrs[st->ptr_ct++] = dst;
  return dst;
}

static void FlatStrStoreFree(FlatStrStore* st) {
  free(st->buf);  st->buf  = nullptr;
  free(st->ptrs); st->ptrs = nullptr;
  st->used = st->cap = st->ptr_ct = st->ptr_cap = 0;
}

// -------------------------------------------------------------------
// Target allele frequencies, read straight out of a --read-freq .afreq.
//
// This path is standalone -- no .pgen and no .psam -- so plink2's
// ReadAlleleFreqs machinery, which indexes by loaded variant, is not
// available.  The file is read here and joined on the variant ID.
//
// Rows whose ALT or ALT_FREQS lists more than one allele (multi-allelic)
// are NOT stored: which ALT the .bim's allele corresponds to would be a
// guess, and a wrong guess inverts an effect.  Such a variant simply has
// no target frequency and comes out AMBIGUOUS_UNRESOLVED.  Rows that are
// too short or carry an unparseable frequency are dropped the same way,
// and counted so the caller can say how many.
//
// A repeated variant ID keeps the first row -- see the insertion loop in
// ReadAfreqTable, which is where that is actually decided.  Orientation
// is resolved per-row against the .bim's own allele (see
// AfreqLookupAlleleFreq), so whichever row wins still yields the
// frequency of the right allele.
// -------------------------------------------------------------------
typedef struct AfreqTableStruct {
  FlatStrStore ids;
  FlatStrStore refs;
  FlatStrStore alts;
  double*      alt_freqs;
  uint32_t     ct;
  uint32_t     cap;
  uint32_t*    htable;
  uint32_t     htsize;
  uint32_t     skipped_ct;   // rows dropped: short, multi-allelic, or unparseable
} AfreqTable;

static void AfreqTableFree(AfreqTable* t) {
  FlatStrStoreFree(&t->ids);
  FlatStrStoreFree(&t->refs);
  FlatStrStoreFree(&t->alts);
  free(t->alt_freqs); t->alt_freqs = nullptr;
  free(t->htable);    t->htable    = nullptr;
  t->ct = t->cap = t->htsize = t->skipped_ct = 0;
}

// 1-based index of the header field named `want`, or 0 if absent.
// A leading '#' on the first field is ignored (#CHROM, #ID, ...).
static uint32_t AfreqHeaderCol(const char* header, const char* want) {
  const uint32_t want_slen = (uint32_t)strlen(want);
  uint32_t idx = 1;
  for (const char* tok = header; tok; tok = NextToken(tok), ++idx) {
    const char* tok_end = CurTokenEnd(tok);
    const char* name = (tok[0] == '#') ? (tok + 1) : tok;
    if (((uint32_t)(tok_end - name) == want_slen) &&
        (memcmp(name, want, want_slen) == 0)) {
      return idx;
    }
  }
  return 0;
}

// 1-based token lookup; nullptr when the line is too short.
static const char* TokenAt(const char* line, uint32_t idx) {
  const char* tok = line;
  for (uint32_t ii = 1; (ii < idx) && tok; ++ii) {
    tok = NextToken(tok);
  }
  return tok;
}

static uint32_t TokenHasComma(const char* tok) {
  for (const char* p = tok; p < CurTokenEnd(tok); ++p) {
    if (*p == ',') return 1;
  }
  return 0;
}

static PglErr ReadAfreqTable(const char* fname, AfreqTable* t) {
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  reterr = SizeAndInitTextStream(fname, bigstack_left() / 8, 1, &txs);
  if (unlikely(reterr)) {
    logerrprintfww("Error: --read-freq: failed to open '%s'.\n", fname);
    goto ReadAfreqTable_ret_1;
  }

  {
    char* header = TextGet(&txs);
    if (unlikely(!header)) {
      logerrprintfww("Error: --read-freq: '%s' is empty.\n", fname);
      reterr = kPglRetMalformedInput;
      goto ReadAfreqTable_ret_1;
    }
    const uint32_t col_id   = AfreqHeaderCol(header, "ID");
    const uint32_t col_ref  = AfreqHeaderCol(header, "REF");
    const uint32_t col_alt  = AfreqHeaderCol(header, "ALT");
    const uint32_t col_freq = AfreqHeaderCol(header, "ALT_FREQS");
    if (unlikely(!col_id || !col_ref || !col_alt || !col_freq)) {
      logerrprintfww(
        "Error: --read-freq: '%s' needs ID, REF, ALT and ALT_FREQS header"
        " columns (a plink2 .afreq).\n", fname);
      reterr = kPglRetMalformedInput;
      goto ReadAfreqTable_ret_1;
    }

    t->cap        = 65536;
    t->skipped_ct = 0;
    t->ids        = FlatStrStoreInit(1 << 20, t->cap);
    t->refs       = FlatStrStoreInit(1 << 18, t->cap);
    t->alts       = FlatStrStoreInit(1 << 18, t->cap);
    t->alt_freqs  = (double*)malloc(t->cap * sizeof(double));
    if (unlikely(!t->ids.buf || !t->refs.buf || !t->alts.buf ||
                 !t->ids.ptrs || !t->refs.ptrs || !t->alts.ptrs ||
                 !t->alt_freqs)) {
      reterr = kPglRetNomem;
      goto ReadAfreqTable_ret_1;
    }

    char* line;
    while ((line = TextGet(&txs)) != nullptr) {
      if (line[0] == '#') continue;
      const char* id_tok   = TokenAt(line, col_id);
      const char* ref_tok  = TokenAt(line, col_ref);
      const char* alt_tok  = TokenAt(line, col_alt);
      const char* freq_tok = TokenAt(line, col_freq);
      // Every reason to drop a row is counted.  A file where a SUBSET of
      // rows is malformed used to be accepted in silence -- only an
      // all-bad file tripped "no usable rows" -- so a truncated or
      // mis-delimited .afreq looked exactly like one that simply had no
      // entry for those variants.
      if (!id_tok || !ref_tok || !alt_tok || !freq_tok) {
        ++t->skipped_ct; continue;
      }
      // Multi-allelic: unavailable rather than guessed at.
      if (TokenHasComma(alt_tok) || TokenHasComma(freq_tok)) {
        ++t->skipped_ct; continue;
      }
      double alt_freq;
      if (!ScantokDouble(freq_tok, &alt_freq)) { ++t->skipped_ct; continue; }

      if (t->ct == t->cap) {
        t->cap *= 2;
        double* nf = (double*)realloc(t->alt_freqs, t->cap * sizeof(double));
        if (unlikely(!nf)) { reterr = kPglRetNomem; goto ReadAfreqTable_ret_1; }
        t->alt_freqs = nf;
      }
      const char* s_id =
          FlatStrStoreAdd(&t->ids, id_tok,
                          (uint32_t)(CurTokenEnd(id_tok) - id_tok));
      const char* s_ref =
          FlatStrStoreAdd(&t->refs, ref_tok,
                          (uint32_t)(CurTokenEnd(ref_tok) - ref_tok));
      const char* s_alt =
          FlatStrStoreAdd(&t->alts, alt_tok,
                          (uint32_t)(CurTokenEnd(alt_tok) - alt_tok));
      if (unlikely(!s_id || !s_ref || !s_alt)) {
        reterr = kPglRetNomem;
        goto ReadAfreqTable_ret_1;
      }
      t->alt_freqs[t->ct] = alt_freq;
      ++t->ct;
    }
    if (TextStreamErrcode2(&txs, &reterr)) {
      TextStreamErrPrint(fname, &txs);
      goto ReadAfreqTable_ret_1;
    }
    if (unlikely(!t->ct)) {
      logerrprintfww("Error: --read-freq: no usable rows in '%s'.\n", fname);
      reterr = kPglRetMalformedInput;
      goto ReadAfreqTable_ret_1;
    }
    if (t->skipped_ct) {
      logerrprintfww(
        "Warning: --read-freq: %u of %u rows in '%s' skipped (too few columns,"
        " multi-allelic, or unparseable ALT_FREQS); those variants have no"
        " target frequency.\n",
        t->skipped_ct, t->skipped_ct + t->ct, fname);
    }

    t->htsize = GetHtableFastSize(t->ct);
    t->htable = (uint32_t*)malloc(t->htsize * sizeof(uint32_t));
    if (unlikely(!t->htable)) { reterr = kPglRetNomem; goto ReadAfreqTable_ret_1; }
    SetAllU32Arr(t->htsize, t->htable);
    // File order.  HtableAddNondup parks the first value inserted for a
    // key at the head of that key's probe sequence and IdHtableFind
    // returns the head, so ascending -- and only ascending -- is what
    // makes the FIRST row of a repeated ID the one that wins.  Inserting
    // descending hands the win to the last row, which for two rows of
    // the same variant at different frequencies is the difference
    // between +beta and -beta.  Same idiom as the .bim table below.
    for (uint32_t ii = 0; ii < t->ct; ++ii) {
      const char* k = t->ids.ptrs[ii];
      HtableAddNondup(k, (uint32_t)strlen(k), t->htsize, ii, t->htable);
    }
  }

 ReadAfreqTable_ret_1:
  CleanupTextStream2(fname, &txs, &reterr);
  return reterr;
}

// Frequency of `allele` in the target, or NaN when unavailable.
//
// ALT_FREQS is the frequency of ALT, and nothing guarantees that the
// .afreq's REF/ALT orientation matches the .bim's A1/A2 -- so the
// correspondence is established explicitly against the allele asked for
// rather than assumed.
//
// *out_id_found reports whether the ID join hit at all, separately from
// whether a frequency came back: a .bim and an .afreq on different ID
// schemes match nothing, and without this the only symptom is a resolved
// count of zero.
static double AfreqLookupAlleleFreq(const AfreqTable* t, const char* varid,
                                    const char* allele,
                                    uint32_t* out_id_found) {
  const double unavailable = nan("");
  *out_id_found = 0;
  if (!t->ct) return unavailable;
  const uint32_t idx = IdHtableFind(varid, t->ids.ptrs, t->htable,
                                    (uint32_t)strlen(varid), t->htsize);
  if (idx == UINT32_MAX) return unavailable;
  *out_id_found = 1;
  if (strcmp(t->alts.ptrs[idx], allele) == 0) return t->alt_freqs[idx];
  if (strcmp(t->refs.ptrs[idx], allele) == 0) return 1.0 - t->alt_freqs[idx];
  return unavailable;   // .afreq names neither of the .bim's alleles
}

// -------------------------------------------------------------------
// RunPgsModel implementation.
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// The target .bim, indexed by (chr, pos).
//
// Loaded once per run and shared by every model.  --pgs-models-list
// harmonises N weight files against ONE target, so re-reading a
// million-line .bim per model would dominate the run for no benefit.
// -------------------------------------------------------------------
typedef struct PgsBimIndexStruct {
  BimEntry*    entries;
  uint32_t     ct;
  uint32_t*    htable;
  uint32_t     htsize;
  FlatStrStore keys;
  uint32_t     mhc_excluded_ct;  // .bim rows inside --exclude-*mhc-*
  uint32_t     no_id_ct;         // rows skipped for a '.' variant ID
  uint32_t     dup_id_ct;        // .bim rows dropped for a repeated ID
} PgsBimIndex;

static void PgsBimIndexPreinit(PgsBimIndex* bx) {
  bx->entries = nullptr;
  bx->ct      = 0;
  bx->htable  = nullptr;
  bx->htsize  = 0;
  bx->mhc_excluded_ct = 0;
  bx->no_id_ct = 0;
  bx->dup_id_ct = 0;
  bx->keys    = {nullptr, 0, 0, nullptr, 0, 0};
}

static void PgsBimIndexFree(PgsBimIndex* bx) {
  if (bx->entries) {
    for (uint32_t ii = 0; ii < bx->ct; ++ii) {
      free(bx->entries[ii].varid);
      free(bx->entries[ii].chr);
      free(bx->entries[ii].pos);
      free(bx->entries[ii].a1);
      free(bx->entries[ii].a0);
    }
    free(bx->entries);
    bx->entries = nullptr;
  }
  free(bx->htable);
  bx->htable = nullptr;
  FlatStrStoreFree(&bx->keys);
}

// Read the .bim (CHR ID CM POS A1 A2, no header) and build the
// (chr \t pos) -> row hash table.
static PglErr LoadBimIndex(const char* bim_fname, MhcRegion mhc_region, uint32_t match_by_id, const char* dup_fname, PgsBimIndex* bx) {
  PglErr reterr = kPglRetSuccess;
  uint32_t bim_cap = 65536;

  TextStream bim_txs;
  PreinitTextStream(&bim_txs);
  reterr = SizeAndInitTextStream(bim_fname, bigstack_left() / 8, 1, &bim_txs);
  if (unlikely(reterr)) {
    logerrprintfww("Error: --bim: failed to open '%s'.\n", bim_fname);
    CleanupTextStream2(bim_fname, &bim_txs, &reterr);
    return reterr;
  }

  bx->entries = (BimEntry*)malloc(bim_cap * sizeof(BimEntry));
  bx->keys    = FlatStrStoreInit(1 << 20, 65536);
  if (!bx->entries || !bx->keys.buf || !bx->keys.ptrs) {
    reterr = kPglRetNomem;
    CleanupTextStream2(bim_fname, &bim_txs, &reterr);
    return kPglRetNomem;
  }

  char* line;
  char key_buf[64];
  while ((line = TextGet(&bim_txs)) != nullptr) {
    if (line[0] == '#') continue;
    const char* t1 = line;
    const char* t2 = NextToken(t1);   // ID
    const char* t3 = NextToken(t2);   // CM
    const char* t4 = NextToken(t3);   // POS
    const char* t5 = NextToken(t4);   // A1
    const char* t6 = NextToken(t5);   // A2
    if (!t6) {
      logerrprintfww("Error: --bim: malformed line in '%s'.\n", bim_fname);
      reterr = kPglRetMalformedInput;
      CleanupTextStream2(bim_fname, &bim_txs, &reterr);
      return kPglRetMalformedInput;
    }

    char chr_tok[32]; const char* chr_end = CurTokenEnd(t1);
    uint32_t chr_slen = (uint32_t)(chr_end - t1);
    if (chr_slen > 31) chr_slen = 31;
    memcpy(chr_tok, t1, chr_slen); chr_tok[chr_slen] = '\0';

    char pos_tok[32]; const char* pos_end = CurTokenEnd(t4);
    uint32_t pos_slen = (uint32_t)(pos_end - t4);
    if (pos_slen > 31) pos_slen = 31;
    memcpy(pos_tok, t4, pos_slen); pos_tok[pos_slen] = '\0';

    uint32_t id_slen = (uint32_t)(CurTokenEnd(t2) - t2);
    char id_tok[kPgsMaxIdBlen];
    if (id_slen >= kPgsMaxIdBlen) id_slen = kPgsMaxIdBlen - 1;
    memcpy(id_tok, t2, id_slen); id_tok[id_slen] = '\0';

    uint32_t key_len;
    if (match_by_id) {
      // '.' is the conventional missing variant ID.  Left in, every such row
      // would collide with every other and first-wins would attach one
      // arbitrary weight to one arbitrary variant.
      if ((id_slen == 1) && (id_tok[0] == '.')) { ++bx->no_id_ct; continue; }
      key_len = id_slen;
      memcpy(key_buf, id_tok, id_slen + 1);
    } else {
      key_len = BuildChrPosKey(chr_tok, pos_tok, key_buf, sizeof(key_buf));
    }
    if (!key_len) continue;

    const uint32_t is_mhc = MhcExcludesChrPos(mhc_region, chr_tok, pos_tok);

    if (bx->ct == bim_cap) {
      bim_cap *= 2;
      BimEntry* nb = (BimEntry*)realloc(bx->entries, bim_cap * sizeof(BimEntry));
      if (!nb) {
        CleanupTextStream2(bim_fname, &bim_txs, &reterr);
        return kPglRetNomem;
      }
      bx->entries = nb;
    }

    const char* stored_key = FlatStrStoreAdd(&bx->keys, key_buf, key_len);
    if (!stored_key) {
      CleanupTextStream2(bim_fname, &bim_txs, &reterr);
      return kPglRetNomem;
    }

    uint32_t a1_slen = (uint32_t)(CurTokenEnd(t5) - t5);
    uint32_t a0_slen = (uint32_t)(CurTokenEnd(t6) - t6);
    BimEntry* be = &bx->entries[bx->ct];
    be->varid = DupStr(t2, id_slen);
    be->chr   = DupStr(chr_tok, (uint32_t)strlen(chr_tok));
    be->pos   = DupStr(pos_tok, (uint32_t)strlen(pos_tok));
    be->a1    = DupStr(t5, a1_slen);
    be->a0    = DupStr(t6, a0_slen);
    be->mhc_excluded = (unsigned char)(is_mhc != 0);
    bx->mhc_excluded_ct += (is_mhc != 0);
    if (!be->varid || !be->a1 || !be->a0) {
      CleanupTextStream2(bim_fname, &bim_txs, &reterr);
      return kPglRetNomem;
    }
    ++bx->ct;
  }
  if (TextStreamErrcode2(&bim_txs, &reterr)) {
    TextStreamErrPrint(bim_fname, &bim_txs);
    CleanupTextStream2(bim_fname, &bim_txs, &reterr);
    return reterr;
  }
  CleanupTextStream2(bim_fname, &bim_txs, &reterr);
  if (unlikely(reterr)) return reterr;

  if (!bx->ct) {
    logerrprintfww("Error: --bim: no variants in '%s'.\n", bim_fname);
    return kPglRetMalformedInput;
  }

  bx->htsize = GetHtableFastSize(bx->ct);
  bx->htable = (uint32_t*)malloc(bx->htsize * sizeof(uint32_t));
  if (!bx->htable) return kPglRetNomem;
  SetAllU32Arr(bx->htsize, bx->htable);

  // HtableAddNondup keeps the first insertion, so a duplicate key is simply
  // not reachable.  That is the long-standing policy for chr:pos and it stays,
  // but duplicate rsIDs are far more common than duplicate positions, so in ID
  // mode they are counted and listed rather than passed over in silence.
  FILE* dup_f = nullptr;
  for (uint32_t ii = 0; ii < bx->ct; ++ii) {
    const char* k = bx->keys.ptrs[ii];
    const uint32_t klen = (uint32_t)strlen(k);
    if (match_by_id) {
      const uint32_t prev = IdHtableFind(k, bx->keys.ptrs, bx->htable, klen, bx->htsize);
      if (prev != UINT32_MAX) {
        ++bx->dup_id_ct;
        if (dup_fname) {
          if (!dup_f) {
            dup_f = fopen(dup_fname, "w");
            if (dup_f) fputs("#ID\tKEPT_ROW\tDROPPED_ROW\n", dup_f);
          }
          // 1-based rows, so they line up with the .bim as a text file.
          if (dup_f) fprintf(dup_f, "%s\t%u\t%u\n", k, prev + 1, ii + 1);
        }
        continue;
      }
    }
    HtableAddNondup(k, klen, bx->htsize, ii, bx->htable);
  }
  if (dup_f) fclose(dup_f);
  return kPglRetSuccess;
}

// What loading the target turned up.  Shared by --pgs-model and
// --pgs-models-list so the two cannot describe the same file differently.
static void LogBimIndexSummary(const char* flag, const PgsBimIndex* bx,
                               MhcRegion mhc_region, uint32_t match_by_id,
                               const char* dup_fname) {
  if (mhc_region != kMhcNone) {
    logprintf("%s: chr6:%u-%u, %u .bim variant%s inside the interval.\n",
              MhcFlagName(mhc_region), MhcStartBp(mhc_region),
              MhcEndBp(mhc_region), bx->mhc_excluded_ct,
              bx->mhc_excluded_ct == 1 ? "" : "s");
  }
  if (!match_by_id) return;
  if (bx->no_id_ct) {
    logprintf("%s: %u .bim variant%s skipped for a missing ('.') ID.\n",
              flag, bx->no_id_ct, bx->no_id_ct == 1 ? "" : "s");
  }
  if (bx->dup_id_ct) {
    // Count here, list in the file: a warning long enough to bury the rest of
    // the log gets skipped wholesale, and the full set is more use greppable.
    logerrprintfww("Warning: %s: %u .bim row%s unreachable because an earlier row carries the same ID; the first of each wins. Full list in %s. Consider re-IDing both files as chr_pos_ref_alt, which is unique by construction.\n",
                   flag, bx->dup_id_ct, bx->dup_id_ct == 1 ? "" : "s",
                   dup_fname? dup_fname : "<none>");
  }
}

// -------------------------------------------------------------------
// One model's harmonisation outcome: the source rows with their
// statuses resolved, plus the counters the summary reports.
// -------------------------------------------------------------------
typedef struct PgsModelStatsStruct {
  uint32_t n_found_pos;
  uint32_t n_ambiguous;
  uint32_t n_matched;
  uint32_t n_flipped;
  uint32_t n_reversed;
  uint32_t n_no_allele;
  uint32_t n_resolved;
  uint32_t n_unresolved;
  uint32_t n_afreq_id_hit;
  uint32_t n_duplicate;
  uint32_t n_mhc_excluded;  // rows dropped for matching an MHC .bim row
  uint32_t n_no_id;         // source rows skipped for a '.' variant ID
} PgsModelStats;

typedef struct PgsModelResultStruct {
  SrcVariant*   vars;
  uint32_t      ct;
  PgsModelStats stats;
  // Backs every string in vars[]; freed in one go.
  StrArena      strs;
  // Retained across models so the vars[] buffer is grown once, to the largest
  // model, rather than allocated and freed for every one.
  uint32_t      cap;
} PgsModelResult;

static void PgsModelResultPreinit(PgsModelResult* r) {
  r->vars = nullptr;
  r->ct   = 0;
  memset(&r->stats, 0, sizeof(r->stats));
  r->strs.head = nullptr;
  r->strs.cur  = nullptr;
  r->cap = 0;
}

static void PgsModelResultFree(PgsModelResult* r) {
  free(r->vars);
  r->vars = nullptr;
  r->ct = 0;
  r->cap = 0;
  StrArenaFree(&r->strs);
}

// Between models: keep the buffers, drop the contents.
static void PgsModelResultRecycle(PgsModelResult* r) {
  r->ct = 0;
  memset(&r->stats, 0, sizeof(r->stats));
  StrArenaReset(&r->strs);
}

// -------------------------------------------------------------------
// Read one weight file, match every row against the .bim index, and
// dedup on (chr, pos).
//
// This is the whole of the matching contract, and it has exactly one
// implementation on purpose: --pgs-model is the one-model case and
// --pgs-models-list calls it per model.  A second copy is what let the
// --pfile and --sscore phenotype checks drift apart.
//
// `flag` names the calling flag for diagnostics, so messages read
// "--pgs-models-list: ..." when that is what the user typed.  `label`
// names what is being matched, for the progress line.
// -------------------------------------------------------------------
static PglErr HarmoniseModel(const PgsBimIndex* bx,
                             const char* model_fname,
                             uint32_t match_by_id,
                             uint32_t col_id,
                             uint32_t col_chr, uint32_t col_pos,
                             uint32_t col_effect, uint32_t col_other,
                             uint32_t col_weight, uint32_t col_freq,
                             uint32_t allow_ambiguous,
                             uint32_t resolve_ambig,
                             double maf_threshold,
                             const AfreqTable* afreq_tab,
                             MhcRegion mhc_region,
                             const char* flag,
                             const char* label,
                             PgsModelResult* out) {
  PglErr reterr = kPglRetSuccess;
  uint32_t src_cap = 4096;
  uint32_t*    dedup_htable = nullptr;
  uint32_t     dedup_htsize = 0;
  FlatStrStore dedup_store  = {nullptr, 0, 0, nullptr, 0, 0};

  TextStream txs;
  PreinitTextStream(&txs);

  reterr = SizeAndInitTextStream(model_fname, bigstack_left() / 8, 1, &txs);
  if (unlikely(reterr)) {
    logerrprintfww("Error: %s: failed to open '%s'.\n", flag, model_fname);
    goto HarmoniseModel_ret_1;
  }

  // Reuse whatever the previous model left behind; only grow.
  if (out->cap < src_cap) {
    SrcVariant* nv = (SrcVariant*)realloc(out->vars, src_cap * sizeof(SrcVariant));
    if (!nv) { reterr = kPglRetNomem; goto HarmoniseModel_ret_1; }
    out->vars = nv;
    out->cap  = src_cap;
  } else {
    src_cap = out->cap;
  }

  {
    double beta_sum = 0.0;
    uint32_t beta_ct = 0;
    uint32_t any_nonpositive = 0;

    uint32_t max_col = col_chr;
    if (col_pos    > max_col) max_col = col_pos;
    if (col_effect > max_col) max_col = col_effect;
    if (col_other  > max_col) max_col = col_other;
    if (col_weight > max_col) max_col = col_weight;
    if (col_freq   > max_col) max_col = col_freq;

    char* line;
    while ((line = TextGet(&txs)) != nullptr) {
      if (line[0] == '#') continue;

      const char* col_ptrs[6] = {nullptr, nullptr, nullptr,
                                 nullptr, nullptr, nullptr};
      // Slots 0 and 1 hold the key columns.  Under ID matching both name the
      // same column, so the "all required fields present" test below stays one
      // expression rather than branching on the mode.
      uint32_t wanted[6] = {match_by_id? col_id : col_chr,
                            match_by_id? col_id : col_pos,
                            col_effect, col_other, col_weight, col_freq};
      const char* tok = line;
      uint32_t cur = 1;
      while (tok && cur <= max_col) {
        for (uint32_t fi = 0; fi < 6; ++fi) {
          if (wanted[fi] == cur) col_ptrs[fi] = tok;
        }
        if (cur < max_col) tok = NextToken(tok);
        ++cur;
      }
      if (!col_ptrs[0] || !col_ptrs[1] || !col_ptrs[2] ||
          !col_ptrs[3] || !col_ptrs[4]) continue;

      char chr_tok[32] = {0};
      char pos_tok[32] = {0};
      char id_tok[kPgsMaxIdBlen] = {0};
      uint32_t id_slen = 0;
      if (match_by_id) {
        id_slen = (uint32_t)(CurTokenEnd(col_ptrs[0]) - col_ptrs[0]);
        if (id_slen >= kPgsMaxIdBlen) id_slen = kPgsMaxIdBlen - 1;
        memcpy(id_tok, col_ptrs[0], id_slen); id_tok[id_slen] = '\0';
        // Same reasoning as on the .bim side: '.' is missing, not a key.
        if ((id_slen == 1) && (id_tok[0] == '.')) {
          ++out->stats.n_no_id;
          continue;
        }
      } else {
        uint32_t chr_slen = (uint32_t)(CurTokenEnd(col_ptrs[0]) - col_ptrs[0]);
        if (chr_slen > 31) chr_slen = 31;
        memcpy(chr_tok, col_ptrs[0], chr_slen); chr_tok[chr_slen] = '\0';

        uint32_t pos_slen = (uint32_t)(CurTokenEnd(col_ptrs[1]) - col_ptrs[1]);
        if (pos_slen > 31) pos_slen = 31;
        memcpy(pos_tok, col_ptrs[1], pos_slen); pos_tok[pos_slen] = '\0';
      }

      uint32_t eff_slen = (uint32_t)(CurTokenEnd(col_ptrs[2]) - col_ptrs[2]);
      uint32_t oth_slen = (uint32_t)(CurTokenEnd(col_ptrs[3]) - col_ptrs[3]);

      double beta;
      if (!ScantokDouble(col_ptrs[4], &beta)) continue;

      double src_freq = nan("");
      if (col_ptrs[5]) {
        double parsed;
        if (ScantokDouble(col_ptrs[5], &parsed)) src_freq = parsed;
      }

      if (out->ct == src_cap) {
        src_cap *= 2;
        SrcVariant* ns = (SrcVariant*)realloc(out->vars,
                                              src_cap * sizeof(SrcVariant));
        if (!ns) { reterr = kPglRetNomem; goto HarmoniseModel_ret_1; }
        out->vars = ns;
        out->cap  = src_cap;
      }

      SrcVariant* sv = &out->vars[out->ct];
      sv->src_id  = match_by_id? StrArenaDup(&out->strs, id_tok, id_slen) : nullptr;
      sv->chr     = match_by_id? nullptr : StrArenaDup(&out->strs, chr_tok, (uint32_t)strlen(chr_tok));
      sv->pos_str = match_by_id? nullptr : StrArenaDup(&out->strs, pos_tok, (uint32_t)strlen(pos_tok));
      sv->src_a1  = StrArenaDup(&out->strs, col_ptrs[2], eff_slen);
      sv->src_a0  = StrArenaDup(&out->strs, col_ptrs[3], oth_slen);
      sv->src_beta = beta;
      sv->src_freq = src_freq;
      sv->status  = kPgsMatchNoPos;
      sv->bim_id  = nullptr;
      sv->bim_chr = nullptr;
      sv->bim_pos = nullptr;
      sv->out_a1  = nullptr;
      sv->out_a2  = nullptr;
      sv->out_beta = 0.0;
      sv->bim_idx = UINT32_MAX;
      const uint32_t key_alloc_ok = match_by_id? (sv->src_id != nullptr)
                                              : (sv->chr && sv->pos_str);
      if (!key_alloc_ok || !sv->src_a1 || !sv->src_a0) {
        reterr = kPglRetNomem;
        goto HarmoniseModel_ret_1;
      }

      beta_sum += beta;
      ++beta_ct;
      if (beta <= 0.0) any_nonpositive = 1;
      ++out->ct;
    }
    if (TextStreamErrcode2(&txs, &reterr)) {
      TextStreamErrPrint(model_fname, &txs);
      goto HarmoniseModel_ret_1;
    }
    CleanupTextStream2(model_fname, &txs, &reterr);
    if (unlikely(reterr)) goto HarmoniseModel_ret_1;

    // Odds-ratio guard.
    if (beta_ct && !any_nonpositive) {
      double mean = beta_sum / (double)beta_ct;
      if (mean >= 0.5 && mean <= 2.0) {
        logputs("Warning: ");
        logputs(flag);
        logputs(": all weights are positive with mean ");
        char tmp[64];
        snprintf(tmp, sizeof(tmp), "%.4g", mean);
        logputs(tmp);
        logputs(
          ", which looks like\nodds ratios rather than betas.  Harmonisation"
          " negates weights for reversed\nalleles, which is only valid for"
          " betas -- supply log(OR) if so.\n");
      }
    }
  }

  if (!out->ct) {
    logerrprintfww("Error: %s: no variants in '%s'.\n", flag, model_fname);
    reterr = kPglRetMalformedInput;
    goto HarmoniseModel_ret_1;
  }

  {
    char key_buf[64];
    // Matching is the part with no output of its own, and on a large weight
    // file against a large target it is where the wall-clock goes.  Same
    // idiom as "Calculating allele frequencies...": the text goes through
    // logputs so the log keeps it, the percentage is written straight to
    // stdout and backspaced away, so the log does not fill with 0%..99%.
    logprintf("%s: matching %s... ", flag, label);
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t pct = 0;
    uint32_t next_print_idx = (out->ct + 99) / 100;

    for (uint32_t ii = 0; ii < out->ct; ++ii) {
      if (ii >= next_print_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (ii * 100LLU) / out->ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_idx = (pct * S_CAST(uint64_t, out->ct) + 99) / 100;
      }
      SrcVariant* sv = &out->vars[ii];
      uint32_t key_len = SrcKey(sv, match_by_id, key_buf,
                                (uint32_t)sizeof(key_buf));
      if (!key_len) { sv->status = kPgsMatchNoPos; continue; }

      uint32_t bim_idx = IdHtableFind(key_buf, bx->keys.ptrs,
                                      bx->htable, key_len, bx->htsize);
      if (bim_idx == UINT32_MAX) { sv->status = kPgsMatchNoPos; continue; }

      const BimEntry* be = &bx->entries[bim_idx];
      // Excluded on the target's authority, not the source's.  Tested after the
      // lookup because that is the only point where the position is known in
      // both modes -- the ID route never sees one on the source side.
      if (be->mhc_excluded) {
        ++out->stats.n_mhc_excluded;
        sv->status = kPgsMatchExcluded;
        continue;
      }
      ++out->stats.n_found_pos;
      sv->bim_chr = StrArenaDup(&out->strs, be->chr, (uint32_t)strlen(be->chr));
      sv->bim_pos = StrArenaDup(&out->strs, be->pos, (uint32_t)strlen(be->pos));
      PgsAmbigResolve resolve;
      const PgsAmbigResolve* resolve_ptr = nullptr;
      if (resolve_ambig) {
        uint32_t id_found = 0;
        resolve.maf_threshold = maf_threshold;
        resolve.src_freq      = sv->src_freq;
        resolve.tgt_a1_freq   = AfreqLookupAlleleFreq(afreq_tab, be->varid,
                                                      be->a1, &id_found);
        out->stats.n_afreq_id_hit += id_found;
        resolve_ptr = &resolve;
      }
      uint32_t out_reversed = 0;
      PgsMatchStatus ms = PgsMatchOne(sv->src_a0, sv->src_a1,
                                      be->a0, be->a1,
                                      allow_ambiguous, resolve_ptr,
                                      &out_reversed);
      sv->status = ms;
      if (ms == kPgsMatchAmbiguous) { ++out->stats.n_ambiguous; }
      else if (ms == kPgsMatchAmbiguousUnresolved) { ++out->stats.n_unresolved; }
      else if (ms == kPgsMatchNoAllele) { ++out->stats.n_no_allele; }
      else {
        if (ms == kPgsMatchResolvedExact || ms == kPgsMatchResolvedReversed) {
          ++out->stats.n_resolved;
        }
        // Matched: negate beta HERE (exactly once) for reversed outcomes.
        double out_beta = sv->src_beta;
        if (out_reversed) { out_beta = -out_beta; ++out->stats.n_reversed; }
        if (ms == kPgsMatchFlipped || ms == kPgsMatchFlippedReversed) {
          ++out->stats.n_flipped;
        }
        sv->bim_id  = StrArenaDup(&out->strs, be->varid, (uint32_t)strlen(be->varid));
      sv->bim_chr = StrArenaDup(&out->strs, be->chr, (uint32_t)strlen(be->chr));
      sv->bim_pos = StrArenaDup(&out->strs, be->pos, (uint32_t)strlen(be->pos));
        sv->out_a1  = StrArenaDup(&out->strs, be->a1, (uint32_t)strlen(be->a1));
        sv->out_a2  = StrArenaDup(&out->strs, be->a0, (uint32_t)strlen(be->a0));
        sv->out_beta = out_beta;
        sv->bim_idx  = bim_idx;
        if (!sv->bim_id || !sv->out_a1) {
          reterr = kPglRetNomem;
          goto HarmoniseModel_ret_1;
        }
        ++out->stats.n_matched;
      }
    }

    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");

    // Dedup: keep the first occurrence of each (chr, pos); later matched
    // occurrences become DUPLICATE.
    const uint32_t nm = out->stats.n_matched ? out->stats.n_matched : 1;
    dedup_store  = FlatStrStoreInit(1 << 16, nm);
    dedup_htsize = GetHtableFastSize(nm);
    dedup_htable = (uint32_t*)malloc(dedup_htsize * sizeof(uint32_t));
    if (!dedup_store.buf || !dedup_store.ptrs || !dedup_htable) {
      reterr = kPglRetNomem;
      goto HarmoniseModel_ret_1;
    }
    SetAllU32Arr(dedup_htsize, dedup_htable);

    for (uint32_t ii = 0; ii < out->ct; ++ii) {
      SrcVariant* sv = &out->vars[ii];
      if (!PgsMatchIsMatched(sv->status)) continue;
      uint32_t key_len = SrcKey(sv, match_by_id, key_buf,
                                (uint32_t)sizeof(key_buf));
      uint32_t dedup_idx = IdHtableFind(key_buf, dedup_store.ptrs,
                                        dedup_htable, key_len, dedup_htsize);
      if (dedup_idx != UINT32_MAX) {
        sv->status = kPgsMatchDuplicate;
        // Arena-owned: drop the reference, never free it here.
        sv->bim_id  = nullptr;
        sv->bim_chr = nullptr;
        sv->bim_pos = nullptr;
        sv->out_a1  = nullptr;
        sv->out_a2  = nullptr;
        sv->out_beta = 0.0;
        sv->bim_idx  = UINT32_MAX;
        ++out->stats.n_duplicate;
        --out->stats.n_matched;
      } else {
        const char* stored = FlatStrStoreAdd(&dedup_store, key_buf, key_len);
        if (!stored) { reterr = kPglRetNomem; goto HarmoniseModel_ret_1; }
        uint32_t idx = dedup_store.ptr_ct - 1;
        HtableAddNondup(key_buf, key_len, dedup_htsize, idx, dedup_htable);
      }
    }
  }

 HarmoniseModel_ret_1:
  CleanupTextStream2(model_fname, &txs, &reterr);
  free(dedup_htable);
  FlatStrStoreFree(&dedup_store);
  return reterr;
}

// The status spelling used by .match_report, in one place so the two
// report writers cannot disagree about it.
static const char* PgsStatusStr(PgsMatchStatus s) {
  switch (s) {
    case kPgsMatchExact:            return "EXACT";
    case kPgsMatchReversed:         return "REVERSED";
    case kPgsMatchFlipped:          return "FLIPPED";
    case kPgsMatchFlippedReversed:  return "FLIPPED_REVERSED";
    case kPgsMatchAmbiguous:        return "AMBIGUOUS";
    case kPgsMatchResolvedExact:    return "RESOLVED_EXACT";
    case kPgsMatchResolvedReversed: return "RESOLVED_REVERSED";
    case kPgsMatchAmbiguousUnresolved: return "AMBIGUOUS_UNRESOLVED";
    case kPgsMatchNoPos:            return "NO_MATCH_POS";
    case kPgsMatchNoAllele:         return "NO_MATCH_ALLELE";
    case kPgsMatchDuplicate:        return "DUPLICATE";
    case kPgsMatchExcluded:         return "EXCLUDED";
    default:                        return "UNKNOWN";
  }
}

// One model's rows of .match_report.  `model_name` is nullptr for
// --pgs-model, which has no MODEL column, and the model's name for
// --pgs-models-list, which prepends one.
// Everything the merge still needs from a model once it has been harmonised
// and its report rows written: which .bim row, and the weight w.r.t. that
// row's A1.
//
// This exists because holding whole PgsModelResults for every model did not
// scale.  A SrcVariant is 88 bytes plus six separate strdups (chr, pos_str,
// src_a1, src_a0, bim_id, out_a1), each paying allocator overhead to store one
// or two characters -- 227 bytes per source variant, measured, or about 11x the
// size of the input file.  249 models totalling 5.9 GB needed ~65 GB and were
// OOM-killed.  At 16 bytes a hit that becomes a few GB.
typedef struct PgsModelHitStruct {
  uint32_t bim_idx;
  double   out_beta;
} PgsModelHit;

// Ascending by .bim row, so the writer can walk every model with a cursor
// instead of materialising a dense union_ct x kept_ct matrix.  Dedup already
// guarantees at most one hit per .bim row per model, so no tie-break is needed
// and the sort does not have to be stable.
static int PgsModelHitCmp(const void* aa, const void* bb) {
  const uint32_t x = ((const PgsModelHit*)aa)->bim_idx;
  const uint32_t y = ((const PgsModelHit*)bb)->bim_idx;
  return (x > y) - (x < y);
}

// Both sides of the join, because the question this file answers is "my source
// said chr6:31000000 -- what did that align to?".  In position mode both pairs
// are populated and a build mismatch shows up directly as a disagreement
// between them; in ID mode the source columns are '.' and SRC_ID identifies the
// row.  An unmatched row is always identifiable by one pair or the other.
static void WriteMatchReportRows(FILE* f, const PgsModelResult* r,
                                 const char* model_name) {
  for (uint32_t ii = 0; ii < r->ct; ++ii) {
    const SrcVariant* sv = &r->vars[ii];
    const char* ss = PgsStatusStr(sv->status);
    if (model_name) fprintf(f, "%s\t", model_name);
    fprintf(f, "%s\t%s\t%s\t",
            sv->src_id? sv->src_id : ".",
            sv->chr? sv->chr : ".",
            sv->pos_str? sv->pos_str : ".");
    fprintf(f, "%s\t%s\t", sv->bim_chr? sv->bim_chr : ".",
            sv->bim_pos? sv->bim_pos : ".");
    if (sv->bim_id) {
      fprintf(f, "%s\t%s\t%.15g\t%s\t%s\t%s\t%.15g\n",
              sv->src_a1, sv->src_a0, sv->src_beta,
              ss, sv->bim_id, sv->out_a1, sv->out_beta);
    } else {
      fprintf(f, "%s\t%s\t%.15g\t%s\tNA\tNA\tNA\n",
              sv->src_a1, sv->src_a0, sv->src_beta, ss);
    }
  }
}

// The per-model counts.  `label` is what precedes each line: the flag
// name for --pgs-model, and "<flag>: <model>" for a list entry.
static void LogModelSummary(const char* label, const PgsModelResult* r,
                            uint32_t allow_ambiguous, uint32_t resolve_ambig,
                            double maf_threshold, uint32_t match_by_id) {
  const PgsModelStats* st = &r->stats;
  logprintf("%s: %u source variant%s, %u found in the .bim by %s.\n",
            label, r->ct, r->ct == 1 ? "" : "s", st->n_found_pos,
            match_by_id? "ID" : "position");
  // Stated separately from r->ct: these rows never became source variants, so
  // they are absent from the .match_report rather than reported as NO_MATCH.
  if (match_by_id) {
    // Prominent by design: ID matching trades a build-mismatch problem for an
    // rsID-versioning one, and unlike a build mismatch that failure is
    // invisible in the coordinates.  A low rate is the only signal that the two
    // files disagree about what an rsID means.
    logprintf("%s: %.1f%% of source variants found in the .bim by ID.\n",
              label, r->ct? (100.0 * S_CAST(double, st->n_found_pos) / S_CAST(double, r->ct)) : 0.0);
  }
  if (st->n_no_id) {
    logprintf("%s: %u source variant%s skipped for a missing ('.') ID.\n",
              label, st->n_no_id, st->n_no_id == 1 ? "" : "s");
  }
  if (st->n_mhc_excluded) {
    logprintf("%s: %u source variant%s dropped for matching an excluded MHC row.\n",
              label, st->n_mhc_excluded, st->n_mhc_excluded == 1 ? "" : "s");
  }
  if (!allow_ambiguous && st->n_ambiguous) {
    logprintf("%s: %u ambiguous (A/T, C/G) dropped; "
              "use --allow-ambiguous to keep them.\n", label, st->n_ambiguous);
  }
  if (resolve_ambig) {
    // The .bim and the .afreq are joined on variant ID, and the two
    // routinely use different ID schemes (rsID vs chr_pos_ref_alt).
    // If not one lookup hit, nothing could ever resolve, and
    // "0 resolved" alone does not say why.
    if (st->n_found_pos && !st->n_afreq_id_hit) {
      logerrputs(
        "Warning: --resolve-ambiguous: none of the .bim variant IDs appear in the"
        " --read-freq\nfile, so no target frequency is available for any"
        " variant.  The two files\nprobably use different ID schemes (rsID vs"
        " chr_pos_ref_alt).\n");
    }
    logprintf("%s: %u ambiguous; %u resolved by frequency "
              "(--resolve-ambiguous, MAF < %g on both sides),\n"
              "             %u unresolved and dropped.\n",
              label, st->n_resolved + st->n_unresolved, st->n_resolved,
              maf_threshold, st->n_unresolved);
  }
  logprintf("%s: %u matched (%u strand-flipped, %u reversed), "
            "%u unmatched on alleles.\n",
            label, st->n_matched, st->n_flipped, st->n_reversed,
            st->n_no_allele);
  if (st->n_duplicate) {
    logprintf("%s: %u duplicate position%s dropped.\n",
              label, st->n_duplicate, st->n_duplicate == 1 ? "" : "s");
  }
}

PglErr RunPgsModel(const char* bim_fname,
                   const char* pgs_model_fname,
                   uint32_t match_by_id,
                   uint32_t col_id,
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
                   char* outname_end) {
  PglErr reterr = kPglRetSuccess;

  PgsBimIndex bx;
  PgsBimIndexPreinit(&bx);
  PgsModelResult res;
  PgsModelResultPreinit(&res);
  AfreqTable afreq_tab;
  memset(&afreq_tab, 0, sizeof(afreq_tab));

  FILE* matched_f = nullptr;
  FILE* report_f  = nullptr;

  memcpy(outname_end, ".duplicate_ids", 15);
  reterr = LoadBimIndex(bim_fname, mhc_region, match_by_id,
                        match_by_id? outname : nullptr, &bx);
  *outname_end = '\0';
  if (unlikely(reterr)) goto RunPgsModel_ret_1;
  logprintf("--pgs-model: %u variant%s loaded from %s.\n",
            bx.ct, bx.ct == 1 ? "" : "s", bim_fname);
  memcpy(outname_end, ".duplicate_ids", 15);
  LogBimIndexSummary("--pgs-model", &bx, mhc_region, match_by_id,
                     match_by_id? outname : nullptr);
  *outname_end = '\0';

  // Target allele frequencies, only when --resolve-ambiguous asked for them.
  // The CLI has already refused --resolve-ambiguous without --read-freq, so
  // read_freq_fname is non-null here.
  if (resolve_ambig) {
    reterr = ReadAfreqTable(read_freq_fname, &afreq_tab);
    if (unlikely(reterr)) goto RunPgsModel_ret_1;
  }

  reterr = HarmoniseModel(&bx, pgs_model_fname, match_by_id, col_id,
                          col_chr, col_pos, col_effect, col_other,
                          col_weight, col_freq, allow_ambiguous,
                          resolve_ambig, maf_threshold, &afreq_tab, mhc_region,
                          "--pgs-model", pgs_model_fname, &res);
  if (unlikely(reterr)) goto RunPgsModel_ret_1;

  memcpy(outname_end, ".matched", 9);
  matched_f = fopen(outname, "w");
  if (!matched_f) {
    logerrprintfww("Error: --pgs-model: failed to create '%s': %s\n",
                   outname, strerror(errno));
    reterr = kPglRetOpenFail;
    goto RunPgsModel_ret_1;
  }
  memcpy(outname_end, ".match_report", 14);
  report_f = fopen(outname, "w");
  if (!report_f) {
    logerrprintfww("Error: --pgs-model: failed to create '%s': %s\n",
                   outname, strerror(errno));
    reterr = kPglRetOpenFail;
    goto RunPgsModel_ret_1;
  }

  // Every identifying field comes from the matched .bim row: the result is
  // scored against that target, so it carries the target's coordinates and
  // alleles rather than the source's, which may be a different build.
  fputs("CHR\tPOS\tID\tA1\tA2\tBETA\n", matched_f);
  for (uint32_t ii = 0; ii != res.ct; ++ii) {
    const SrcVariant* sv = &res.vars[ii];
    if (sv->bim_id) {
      fprintf(matched_f, "%s\t%s\t%s\t%s\t%s\t%.15g\n",
              sv->bim_chr, sv->bim_pos, sv->bim_id, sv->out_a1,
              sv->out_a2, sv->out_beta);
    }
  }

  fputs("#SRC_ID\tSRC_CHR\tSRC_POS\tBIM_CHR\tBIM_POS\t"
        "SRC_EFFECT\tSRC_OTHER\tSRC_BETA\tSTATUS\tBIM_ID\tOUT_A1\tOUT_BETA\n",
        report_f);
  WriteMatchReportRows(report_f, &res, nullptr);

  fclose(matched_f); matched_f = nullptr;
  fclose(report_f);  report_f  = nullptr;

  *outname_end = '\0';
  LogModelSummary("--pgs-model", &res, allow_ambiguous, resolve_ambig,
                  maf_threshold, match_by_id);
  logprintfww("--pgs-model: results written to %s.matched and %s.match_report .\n",
              outname, outname);

 RunPgsModel_ret_1:
  if (matched_f) { fclose(matched_f); matched_f = nullptr; }
  if (report_f)  { fclose(report_f);  report_f  = nullptr; }
  PgsModelResultFree(&res);
  PgsBimIndexFree(&bx);
  AfreqTableFree(&afreq_tab);
  return reterr;
}

// -------------------------------------------------------------------
// Column name for one model: basename with everything from the first
// '.' removed, so "PGS000123.txt.gz" -> "PGS000123".
//
// Returns 0 when nothing is left (a path like "/a/.txt"), which the
// caller reports rather than emitting a nameless column.
// -------------------------------------------------------------------
static uint32_t PgsModelNameFromPath(const char* path, char* out,
                                     uint32_t out_cap) {
  const char* base = path;
  for (const char* p = path; *p; ++p) {
    if (*p == '/') base = p + 1;
  }
  uint32_t n = 0;
  while (base[n] && base[n] != '.' && n + 1 < out_cap) ++n;
  memcpy(out, base, n);
  out[n] = '\0';
  return n;
}

PglErr RunPgsModelsList(const char* bim_fname,
                        const char* list_fname,
                        uint32_t match_by_id,
                        uint32_t col_id,
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
                        char* outname_end) {
  PglErr reterr = kPglRetSuccess;

  PgsBimIndex bx;
  PgsBimIndexPreinit(&bx);
  AfreqTable afreq_tab;
  memset(&afreq_tab, 0, sizeof(afreq_tab));

  char**  paths     = nullptr;
  char**  names     = nullptr;
  uint32_t model_ct = 0;
  uint32_t model_cap = 0;

  // One model is live at a time; `hits` keeps only the compact residue.
  PgsModelResult res;
  PgsModelResultPreinit(&res);
  PgsModelHit** hits      = nullptr;   // per model, sorted by bim_idx
  uint32_t*     hit_ct    = nullptr;
  uint32_t      hits_ct   = 0;         // models whose hits[] slot is live

  uint32_t* keep_col  = nullptr;   // model -> output column, or UINT32_MAX
  uint32_t* union_row = nullptr;   // bim row -> output row, or UINT32_MAX
  uint32_t* cursor    = nullptr;   // per kept column, position in its hits[]
  double*   rowbuf    = nullptr;   // one output row
  uint32_t  kept_ct   = 0;
  uint32_t  union_ct  = 0;

  FILE* matched_f = nullptr;
  FILE* report_f  = nullptr;

  // ---- Step 1: read the list of model paths. ----
  {
    TextStream list_txs;
    PreinitTextStream(&list_txs);
    reterr = SizeAndInitTextStream(list_fname, bigstack_left() / 8, 1, &list_txs);
    if (unlikely(reterr)) {
      logerrprintfww("Error: --pgs-models-list: failed to open '%s'.\n",
                     list_fname);
      CleanupTextStream2(list_fname, &list_txs, &reterr);
      goto RunPgsModelsList_ret_1;
    }
    model_cap = 16;
    paths = (char**)malloc(model_cap * sizeof(char*));
    names = (char**)malloc(model_cap * sizeof(char*));
    if (!paths || !names) {
      CleanupTextStream2(list_fname, &list_txs, &reterr);
      reterr = kPglRetNomem;
      goto RunPgsModelsList_ret_1;
    }
    char* line;
    while ((line = TextGet(&list_txs)) != nullptr) {
      // Skip comments and whitespace-only lines, so a list can be
      // annotated.  A completely empty interior line never reaches here:
      // plink2_text refuses those outright (trailing ones are fine), so
      // this cannot silently drop a path.
      const char* s = line;
      while (*s == ' ' || *s == '\t') ++s;
      if (*s == '\0' || *s == '#') continue;
      const char* end = CurTokenEnd(s);
      const uint32_t slen = (uint32_t)(end - s);

      if (model_ct == model_cap) {
        model_cap *= 2;
        char** np = (char**)realloc(paths, model_cap * sizeof(char*));
        char** nn = (char**)realloc(names, model_cap * sizeof(char*));
        if (np) paths = np;
        if (nn) names = nn;
        if (!np || !nn) {
          CleanupTextStream2(list_fname, &list_txs, &reterr);
          reterr = kPglRetNomem;
          goto RunPgsModelsList_ret_1;
        }
      }
      paths[model_ct] = DupStr(s, slen);
      if (!paths[model_ct]) {
        CleanupTextStream2(list_fname, &list_txs, &reterr);
        reterr = kPglRetNomem;
        goto RunPgsModelsList_ret_1;
      }
      char namebuf[256];
      if (!PgsModelNameFromPath(paths[model_ct], namebuf, sizeof(namebuf))) {
        logerrprintfww("Error: --pgs-models-list: cannot derive a model name "
                       "from '%s'.\n", paths[model_ct]);
        free(paths[model_ct]);
        CleanupTextStream2(list_fname, &list_txs, &reterr);
        reterr = kPglRetMalformedInput;
        goto RunPgsModelsList_ret_1;
      }
      names[model_ct] = DupStr(namebuf, (uint32_t)strlen(namebuf));
      if (!names[model_ct]) {
        free(paths[model_ct]);
        CleanupTextStream2(list_fname, &list_txs, &reterr);
        reterr = kPglRetNomem;
        goto RunPgsModelsList_ret_1;
      }
      ++model_ct;
    }
    if (TextStreamErrcode2(&list_txs, &reterr)) {
      TextStreamErrPrint(list_fname, &list_txs);
      CleanupTextStream2(list_fname, &list_txs, &reterr);
      goto RunPgsModelsList_ret_1;
    }
    CleanupTextStream2(list_fname, &list_txs, &reterr);
    if (unlikely(reterr)) goto RunPgsModelsList_ret_1;
  }

  if (!model_ct) {
    logerrprintfww("Error: --pgs-models-list: no model paths in '%s'.\n",
                   list_fname);
    reterr = kPglRetMalformedInput;
    goto RunPgsModelsList_ret_1;
  }

  // Duplicate names are fatal, not silently disambiguated: the names
  // become column headers that the user addresses with --score-col-nums,
  // so handing back a header they did not write would be worse than
  // stopping.  Quadratic, but model counts are small and this runs once.
  for (uint32_t ii = 0; ii != model_ct; ++ii) {
    for (uint32_t jj = ii + 1; jj != model_ct; ++jj) {
      if (strcmp(names[ii], names[jj]) == 0) {
        logerrprintfww("Error: --pgs-models-list: '%s' and '%s' both give the "
                       "model name '%s'.\n", paths[ii], paths[jj], names[ii]);
        reterr = kPglRetMalformedInput;
        goto RunPgsModelsList_ret_1;
      }
    }
  }

  // ---- Step 2: the target, read once for every model. ----
  memcpy(outname_end, ".duplicate_ids", 15);
  reterr = LoadBimIndex(bim_fname, mhc_region, match_by_id,
                        match_by_id? outname : nullptr, &bx);
  *outname_end = '\0';
  if (unlikely(reterr)) goto RunPgsModelsList_ret_1;
  logprintf("--pgs-models-list: %u variant%s loaded from %s.\n",
            bx.ct, bx.ct == 1 ? "" : "s", bim_fname);
  memcpy(outname_end, ".duplicate_ids", 15);
  LogBimIndexSummary("--pgs-models-list", &bx, mhc_region, match_by_id,
                     match_by_id? outname : nullptr);
  *outname_end = '\0';

  if (resolve_ambig) {
    reterr = ReadAfreqTable(read_freq_fname, &afreq_tab);
    if (unlikely(reterr)) goto RunPgsModelsList_ret_1;
  }

  // ---- Step 3: harmonise every model against it. ----
  //
  // One model is held at a time.  Its report rows are written as soon as it is
  // harmonised and it is then reduced to (bim_idx, out_beta) pairs, so peak
  // memory tracks the LARGEST model plus the compact residue of the rest,
  // rather than the sum of every model's parsed records.
  hits   = (PgsModelHit**)calloc(model_ct, sizeof(PgsModelHit*));
  hit_ct = (uint32_t*)calloc(model_ct, sizeof(uint32_t));
  if (!hits || !hit_ct) { reterr = kPglRetNomem; goto RunPgsModelsList_ret_1; }
  hits_ct = model_ct;

  // Opened before the loop because the report is now streamed into it.  A run
  // that fails midway therefore leaves a partial .match_report; it also exits
  // nonzero and says why, which is the signal to trust.
  memcpy(outname_end, ".match_report", 14);
  report_f = fopen(outname, "w");
  if (!report_f) {
    logerrprintfww("Error: --pgs-models-list: failed to create '%s': %s\n",
                   outname, strerror(errno));
    reterr = kPglRetOpenFail;
    goto RunPgsModelsList_ret_1;
  }
  fputs("#MODEL\tSRC_ID\tSRC_CHR\tSRC_POS\tBIM_CHR\tBIM_POS\t"
        "SRC_EFFECT\tSRC_OTHER\tSRC_BETA\tSTATUS\tBIM_ID\tOUT_A1\tOUT_BETA\n", report_f);

  keep_col = (uint32_t*)malloc(model_ct * sizeof(uint32_t));
  if (!keep_col) { reterr = kPglRetNomem; goto RunPgsModelsList_ret_1; }

  for (uint32_t mi = 0; mi != model_ct; ++mi) {
    // HarmoniseModel sizes its text stream from bigstack_left(), and bigstack
    // is a stack allocator: CleanupTextStream2 closes the stream but only a
    // BigstackReset returns the buffer.  Without this mark every model in the
    // list took a fresh slice and kept it, which is what actually made peak
    // memory climb with model count -- not the parsed records, which are freed.
    unsigned char* model_bigstack_mark = g_bigstack_base;

    char prog_label[288];
    snprintf(prog_label, sizeof(prog_label), "%s (%u of %u)",
             names[mi], mi + 1, model_ct);
    reterr = HarmoniseModel(&bx, paths[mi], match_by_id, col_id,
                            col_chr, col_pos, col_effect, col_other,
                            col_weight, col_freq, allow_ambiguous,
                            resolve_ambig, maf_threshold, &afreq_tab, mhc_region,
                            "--pgs-models-list", prog_label, &res);
    if (unlikely(reterr)) goto RunPgsModelsList_ret_1;
    const PgsModelStats* st = &res.stats;
    logprintf("--pgs-models-list: %s: %u source variant%s, %u matched "
              "(%u strand-flipped, %u reversed).\n",
              names[mi], res.ct, res.ct == 1 ? "" : "s",
              st->n_matched, st->n_flipped, st->n_reversed);

    WriteMatchReportRows(report_f, &res, names[mi]);

    // bim_idx is set only on a matched row and cleared again on a duplicate,
    // so this count equals n_matched and no .bim row repeats within a model.
    if (st->n_matched) {
      hits[mi] = (PgsModelHit*)malloc(st->n_matched * sizeof(PgsModelHit));
      if (!hits[mi]) { reterr = kPglRetNomem; goto RunPgsModelsList_ret_1; }
      uint32_t nh = 0;
      for (uint32_t ii = 0; ii != res.ct; ++ii) {
        const uint32_t bi = res.vars[ii].bim_idx;
        if (bi == UINT32_MAX) continue;
        hits[mi][nh].bim_idx  = bi;
        hits[mi][nh].out_beta = res.vars[ii].out_beta;
        ++nh;
      }
      assert(nh == st->n_matched);
      hit_ct[mi] = nh;
      qsort(hits[mi], nh, sizeof(PgsModelHit), PgsModelHitCmp);
      keep_col[mi] = kept_ct++;
    } else {
      // An all-zero column is a constant predictor: degenerate as a PGS,
      // and it would occupy a --score-col-nums slot that predicts nothing.
      logerrprintfww("Warning: --pgs-models-list: '%s' matched no variants; "
                     "its column is omitted.\n", paths[mi]);
      keep_col[mi] = UINT32_MAX;
    }

    PgsModelResultRecycle(&res);
    // Nothing kept above lives on bigstack -- vars, the string arena and hits
    // are all malloc'd -- so the model's stream buffer can go back now.
    BigstackReset(model_bigstack_mark);
  }

  // ---- Step 4: build the union. ----
  if (!kept_ct) {
    logerrputs("Error: --pgs-models-list: no model matched any variant.\n");
    reterr = kPglRetMalformedInput;
    goto RunPgsModelsList_ret_1;
  }

  union_row = (uint32_t*)malloc(bx.ct * sizeof(uint32_t));
  if (!union_row) { reterr = kPglRetNomem; goto RunPgsModelsList_ret_1; }
  SetAllU32Arr(bx.ct, union_row);
  for (uint32_t mi = 0; mi != model_ct; ++mi) {
    for (uint32_t ii = 0; ii != hit_ct[mi]; ++ii) {
      union_row[hits[mi][ii].bim_idx] = 0;   // mark present
    }
  }
  // Number the marked rows in .bim order, so the output comes out sorted by
  // target position without a sort.
  for (uint32_t bi = 0; bi != bx.ct; ++bi) {
    if (union_row[bi] == 0) union_row[bi] = union_ct++;
    else union_row[bi] = UINT32_MAX;
  }

  // ---- Step 5: write.  Absent (model, variant) is 0. ----
  //
  // Written a row at a time by advancing one cursor per model through its
  // bim_idx-sorted hits, rather than materialising a union_ct x kept_ct matrix.
  // That matrix was the second scaling term: 249 models against an 8M-variant
  // .bim would be 16 GB of doubles on its own.
  cursor = (uint32_t*)calloc(model_ct, sizeof(uint32_t));
  rowbuf = (double*)malloc(kept_ct * sizeof(double));
  if (!cursor || !rowbuf) { reterr = kPglRetNomem; goto RunPgsModelsList_ret_1; }

  memcpy(outname_end, ".matched", 9);
  matched_f = fopen(outname, "w");
  if (!matched_f) {
    logerrprintfww("Error: --pgs-models-list: failed to create '%s': %s\n",
                   outname, strerror(errno));
    reterr = kPglRetOpenFail;
    goto RunPgsModelsList_ret_1;
  }

  fputs("CHR\tPOS\tID\tA1\tA2", matched_f);
  for (uint32_t mi = 0; mi != model_ct; ++mi) {
    if (keep_col[mi] != UINT32_MAX) fprintf(matched_f, "\t%s", names[mi]);
  }
  fputc('\n', matched_f);
  for (uint32_t bi = 0; bi != bx.ct; ++bi) {
    if (union_row[bi] == UINT32_MAX) continue;
    for (uint32_t cc = 0; cc != kept_ct; ++cc) rowbuf[cc] = 0.0;
    for (uint32_t mi = 0; mi != model_ct; ++mi) {
      const uint32_t col = keep_col[mi];
      if (col == UINT32_MAX) continue;
      // Monotone: bi only increases, and hits[mi] is sorted by bim_idx, so
      // every hit is examined at most once across the whole write.
      uint32_t* cp = &(cursor[mi]);
      while ((*cp < hit_ct[mi]) && (hits[mi][*cp].bim_idx < bi)) ++(*cp);
      if ((*cp < hit_ct[mi]) && (hits[mi][*cp].bim_idx == bi)) {
        rowbuf[col] = hits[mi][*cp].out_beta;
        ++(*cp);
      }
    }
    const BimEntry* obe = &bx.entries[bi];
    fprintf(matched_f, "%s\t%s\t%s\t%s\t%s",
            obe->chr, obe->pos, obe->varid, obe->a1, obe->a0);
    for (uint32_t cc = 0; cc != kept_ct; ++cc) {
      fprintf(matched_f, "\t%.15g", rowbuf[cc]);
    }
    fputc('\n', matched_f);
  }

  fclose(report_f);  report_f  = nullptr;
  fclose(matched_f); matched_f = nullptr;

  *outname_end = '\0';
  logprintf("--pgs-models-list: %u model%s, %u variant%s in the union.\n",
            kept_ct, kept_ct == 1 ? "" : "s",
            union_ct, union_ct == 1 ? "" : "s");
  logprintfww("--pgs-models-list: results written to %s.matched and "
              "%s.match_report .\n", outname, outname);

 RunPgsModelsList_ret_1:
  if (matched_f) { fclose(matched_f); matched_f = nullptr; }
  if (report_f)  { fclose(report_f);  report_f  = nullptr; }
  free(rowbuf);
  free(cursor);
  free(union_row);
  free(keep_col);
  PgsModelResultFree(&res);
  if (hits) {
    for (uint32_t mi = 0; mi != hits_ct; ++mi) free(hits[mi]);
    free(hits);
  }
  free(hit_ct);
  if (paths) {
    for (uint32_t mi = 0; mi != model_ct; ++mi) free(paths[mi]);
    free(paths);
  }
  if (names) {
    for (uint32_t mi = 0; mi != model_ct; ++mi) free(names[mi]);
    free(names);
  }
  PgsBimIndexFree(&bx);
  AfreqTableFree(&afreq_tab);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
