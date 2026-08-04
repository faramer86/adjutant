#ifndef __PLINK2_MEMORY_H__
#define __PLINK2_MEMORY_H__

// This library is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.


// Double-ended stack allocator, and a few support functions (e.g.
// GetMemAvailableKib to estimate how much memory we can afford to use).

#include <assert.h>

#include "plink2_base.h"
#include "plink2_simd.h"

#ifdef __cplusplus
namespace plink2 {
#endif

uint64_t DetectMib();

// On more recent Linux builds, the MemAvailable value reported by
// "cat /proc/meminfo" is a pretty accurate estimate of how much we can afford
// to allocate before introducing substantial risk of OOM-killer action.
// Returns (~0LLU) if this estimate is unavailable.
// Requires a text buffer large enough to fit longest line in /proc/meminfo.
uint64_t GetMemAvailableKib(uint32_t textbuf_blen, char* textbuf);

// ensure vector-alignment
CONSTI32(kEndAllocAlign, MAXV(kBytesPerVec, 16));

HEADER_INLINE void* arena_alloc(unsigned char* arena_top, uintptr_t size, unsigned char** arena_bottom_ptr) {
  size = RoundUpPow2(size, kCacheline);
  if (S_CAST(uintptr_t, arena_top - (*arena_bottom_ptr)) < size) {
    g_failed_alloc_attempt_size = size;
    return nullptr;
  }
  return arena_alloc_raw(size, arena_bottom_ptr);
}
HEADER_INLINE BoolErr arena_alloc_d(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, double** d_arr_ptr) {
  *d_arr_ptr = S_CAST(double*, arena_alloc(arena_top, ct * sizeof(double), arena_bottom_ptr));
  return !(*d_arr_ptr);
}
HEADER_INLINE BoolErr arena_alloc_u32(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, uint32_t** u32_arr_ptr) {
  *u32_arr_ptr = S_CAST(uint32_t*, arena_alloc(arena_top, ct * sizeof(int32_t), arena_bottom_ptr));
  return !(*u32_arr_ptr);
}
HEADER_INLINE BoolErr arena_alloc_kcp(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, const char*** kcp_arr_ptr) {
  *kcp_arr_ptr = S_CAST(const char**, arena_alloc(arena_top, ct * sizeof(intptr_t), arena_bottom_ptr));
  return !(*kcp_arr_ptr);
}


HEADER_INLINE void ArenaBaseSet(const void* unaligned_base, unsigned char** arena_bottom_ptr) {
  *arena_bottom_ptr = R_CAST(unsigned char*, RoundUpPow2(R_CAST(uintptr_t, unaligned_base), kCacheline));
}
HEADER_INLINE void* arena_end_alloc_raw(uintptr_t size, unsigned char** arena_top_ptr) {
  assert(!(size % kEndAllocAlign));
  unsigned char* alloc_ptr = *arena_top_ptr;
  alloc_ptr -= size;
  *arena_top_ptr = alloc_ptr;
  return alloc_ptr;
}

HEADER_INLINE void* arena_end_alloc(unsigned char* arena_bottom, uintptr_t size, unsigned char** arena_top_ptr) {
  size = RoundUpPow2(size, kEndAllocAlign);
  if (unlikely(S_CAST(uintptr_t, (*arena_top_ptr) - arena_bottom) < size)) {
    g_failed_alloc_attempt_size = size;
    return nullptr;
  }
  return arena_end_alloc_raw(size, arena_top_ptr);
}
HEADER_INLINE BoolErr arena_end_alloc_u32(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, uint32_t** u32_arr_ptr) {
  *u32_arr_ptr = S_CAST(uint32_t*, arena_end_alloc(arena_bottom, ct * sizeof(int32_t), arena_top_ptr));
  return !(*u32_arr_ptr);
}
#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_MEMORY_H__
