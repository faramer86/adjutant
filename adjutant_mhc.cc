#include "adjutant_mhc.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Coordinates, inclusive, 1-based, chromosome 6.
//
// MHC is the Genome Reference Consortium's MHC region:
//   https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC
// xMHC is the extended region, HIST1H2AA to RPL12P1:
//   https://cloufield.github.io/gwaslab/HLA/
//
// Every one of these eight numbers is published for its own build.  None is
// lifted over, and none may be "corrected" to look consistent with another:
// the xMHC hg38 start (25726063) is LOWER than the hg19 start (25726291), so
// the region does not shift uniformly between builds and any offset applied to
// convert one into the other would be wrong at one end.
//
// The check that these are the same region in both builds is the span, not the
// endpoints: MHC is 4970558 bp in hg19 and 4970458 in hg38, xMHC 7642131 and
// 7674582.
typedef struct MhcIntervalStruct {
  uint32_t start_bp;
  uint32_t end_bp;
  const char* flag_name;
} MhcInterval;

static const MhcInterval kMhcIntervals[] = {
  {0, 0, nullptr},                              // kMhcNone
  {28477797, 33448354, "--exclude-mhc-hg19"},   // kMhcHg19
  {25726291, 33368421, "--exclude-xmhc-hg19"},  // kMhcXHg19
  {28510120, 33480577, "--exclude-mhc-hg38"},   // kMhcHg38
  {25726063, 33400644, "--exclude-xmhc-hg38"}   // kMhcXHg38
};

static const MhcInterval* MhcLookup(MhcRegion region) {
  const uint32_t idx = S_CAST(uint32_t, region);
  if (idx >= (sizeof(kMhcIntervals) / sizeof(kMhcIntervals[0]))) {
    return &(kMhcIntervals[0]);
  }
  return &(kMhcIntervals[idx]);
}

uint32_t MhcStartBp(MhcRegion region) {
  return MhcLookup(region)->start_bp;
}

uint32_t MhcEndBp(MhcRegion region) {
  return MhcLookup(region)->end_bp;
}

uint32_t MhcExcludesBp(MhcRegion region, uint32_t bp) {
  const MhcInterval* ivp = MhcLookup(region);
  // Explicit, not folded into the range test: kMhcNone is {0, 0}, and a bare
  // (bp >= 0 && bp <= 0) would exclude a variant at position 0 when no region
  // was selected at all.  Position 0 occurs in real filesets.
  if (!ivp->end_bp) {
    return 0;
  }
  return (bp >= ivp->start_bp) && (bp <= ivp->end_bp);
}

const char* MhcFlagName(MhcRegion region) {
  return MhcLookup(region)->flag_name;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
