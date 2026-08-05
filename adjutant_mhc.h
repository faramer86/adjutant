#ifndef __ADJUTANT_MHC_H__
#define __ADJUTANT_MHC_H__

// Exclusion of the major histocompatibility complex on chromosome 6.
//
// The MHC's long-range LD and extreme polymorphism make it behave unlike the
// rest of the genome, so PGS work routinely drops it.  Doing that by hand means
// writing an interval file and hoping its coordinates match the build the data
// is actually in.
//
// This is the single definition of those coordinates.  Both the scoring route
// and the --bim/--pgs-model harmonisation route call in here; neither carries
// its own copy.  That is deliberate -- the --score/--sscore phenotype
// validators drifted apart precisely because each route owned a private copy of
// the rule, and a coordinate table is worse in that respect than a validator,
// because a stale copy produces plausible output instead of an error.

#include <stdint.h>

#include "include/plink2_base.h"

#ifdef __cplusplus
namespace plink2 {
#endif

typedef enum {
  kMhcNone = 0,
  kMhcHg19,
  kMhcXHg19,
  kMhcHg38,
  kMhcXHg38
} MhcRegion;

// Interval bounds, inclusive at both ends, 1-based -- the same convention as
// the position column of a .bim/.pvar.  Both are 0 for kMhcNone.
uint32_t MhcStartBp(MhcRegion region);
uint32_t MhcEndBp(MhcRegion region);

// True when a chromosome 6 variant at this position falls in the region.
// The caller establishes that the variant IS on chromosome 6: the scoring route
// has a ChrInfo and the harmonisation route has a text token, so each does it in
// its own idiom and this stays free of both.  Always false for kMhcNone.
uint32_t MhcExcludesBp(MhcRegion region, uint32_t bp);

// The flag that selected this region, e.g. "--exclude-mhc-hg19", for log and
// error messages.  Returns nullptr for kMhcNone.
const char* MhcFlagName(MhcRegion region);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __ADJUTANT_MHC_H__
