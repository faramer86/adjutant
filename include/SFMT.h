#pragma once
/**
 * @file SFMT.h
 *
 * @brief SIMD oriented Fast Mersenne Twister(SFMT) pseudorandom
 * number generator using C structure.
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (The University of Tokyo)
 *
 * Copyright (C) 2006, 2007 Mutsuo Saito, Makoto Matsumoto and Hiroshima
 * University.
 * Copyright (C) 2012 Mutsuo Saito, Makoto Matsumoto, Hiroshima
 * University and The University of Tokyo.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software:
 *
 * * * * * *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * - Neither the names of Hiroshima University, The University of Tokyo nor the
 *   names of its contributors may be used to endorse or promote products
 *   derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * * * * * *
 *
 * @note We assume that your system has inttypes.h.  If your system
 * doesn't have inttypes.h, you have to typedef uint32_t and uint64_t,
 * and you have to define PRIu64 and PRIx64 in this file as follows:
 * @verbatim
 typedef unsigned int uint32_t
 typedef unsigned long long uint64_t
 #define PRIu64 "llu"
 #define PRIx64 "llx"
@endverbatim
 * uint32_t must be exactly 32-bit unsigned integer type (no more, no
 * less), and uint64_t must be exactly 64-bit unsigned integer type.
 * PRIu64 and PRIx64 are used for printf function to print 64-bit
 * unsigned int and 64-bit unsigned int in hexadecimal format.
 */

#define __USE_MINGW_ANSI_STDIO 1

#ifdef __LP64__
#  ifdef __x86_64__
#    include <emmintrin.h>  // IWYU pragma: export
#  else
#    define SIMDE_ENABLE_NATIVE_ALIASES
#    ifdef IGNORE_BUNDLED_SIMDE
#      include <simde/x86/sse2.h>  // IWYU pragma: export
#    else
#      include "../simde/x86/sse2.h"  // IWYU pragma: export
#    endif
#  endif
#endif

#ifndef SFMTST_H
#define SFMTST_H
#if defined(__cplusplus)
extern "C" {
#endif

#include <assert.h>

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
  #include <inttypes.h>
#elif defined(_MSC_VER) || defined(__BORLANDC__)
  typedef unsigned int uint32_t;
  typedef unsigned __int64 uint64_t;
  #define inline __inline
#else
  #include <inttypes.h>
  #if defined(__GNUC__)
    #define inline __inline__
  #endif
#endif

#ifndef PRIu64
  #define PRIu64 "llu"
  #define PRIx64 "llx"
#endif


  // fold in SFMT-params.h, etc. to reduce compilation time a bit
#if !defined(SFMT_MEXP)
  #define SFMT_MEXP 19937
#endif
/*-----------------
  BASIC DEFINITIONS
  -----------------*/
/** Mersenne Exponent. The period of the sequence
 *  is a multiple of 2^MEXP-1.
 * #define SFMT_MEXP 19937 */
/** SFMT generator has an internal state array of 128-bit integers,
 * and N is its size. */
#define SFMT_N (SFMT_MEXP / 128 + 1)
/** N32 is the size of internal state array when regarded as an array
 * of 32-bit integers.*/
#define SFMT_N32 (SFMT_N * 4)
/** N64 is the size of internal state array when regarded as an array
 * of 64-bit integers.*/
#define SFMT_N64 (SFMT_N * 2)
#define SFMT_POS1	122
#define SFMT_SL1	18
#define SFMT_SL2	1
#define SFMT_SR1	11
#define SFMT_SR2	1
#define SFMT_MSK1	0xdfffffefU
#define SFMT_MSK2	0xddfecb7fU
#define SFMT_MSK3	0xbffaffffU
#define SFMT_MSK4	0xbffffff6U
#define SFMT_PARITY1	0x00000001U
#define SFMT_PARITY2	0x00000000U
#define SFMT_PARITY3	0x00000000U
#define SFMT_PARITY4	0x13c9e684U
#define SFMT_IDSTR	"SFMT-19937:122-18-1-11-1:dfffffef-ddfecb7f-bffaffff-bffffff6"


/*------------------------------------------
  128-bit SIMD like data type for standard C
  ------------------------------------------*/
#ifdef __LP64__
/** 128-bit data structure */
union W128_T {
    uint32_t u[4];
    uint64_t u64[2];
    __m128i si;
};
#else
/** 128-bit data structure */
union W128_T {
    uint32_t u[4];
    uint64_t u64[2];
};
#endif

/** 128-bit data type */
typedef union W128_T w128_t;

/**
 * SFMT internal state
 */
struct SFMT_T {
    /** the 128-bit internal state array */
    w128_t state[SFMT_N];
    /** index counter to the 32-bit internal state array */
    int idx;
};

typedef struct SFMT_T sfmt_t;

void sfmt_init_gen_rand(sfmt_t * sfmt, uint32_t seed);
void sfmt_init_by_array(sfmt_t * sfmt, uint32_t * init_key, int key_length);

/**
 * This function generates and returns 32-bit pseudorandom number.
 * init_gen_rand or init_by_array must be called before this function.
 * @param sfmt SFMT internal state
 * @return 32-bit pseudorandom number
 */
/**
 * This function generates and returns 64-bit pseudorandom number.
 * init_gen_rand or init_by_array must be called before this function.
 * The function gen_rand64 should not be called after gen_rand32,
 * unless an initialization is again executed.
 * @param sfmt SFMT internal state
 * @return 64-bit pseudorandom number
 */
/* =================================================
   The following real versions are due to Isaku Wada
   ================================================= */
/**
 * converts an unsigned 32-bit number to a double on [0,1]-real-interval.
 * @param v 32-bit unsigned integer
 * @return double on [0,1]-real-interval
 */
/**
 * generates a random number on [0,1]-real-interval
 * @param sfmt SFMT internal state
 * @return double on [0,1]-real-interval
 */
/**
 * converts an unsigned 32-bit integer to a double on [0,1)-real-interval.
 * @param v 32-bit unsigned integer
 * @return double on [0,1)-real-interval
 */
/**
 * generates a random number on [0,1)-real-interval
 * @param sfmt SFMT internal state
 * @return double on [0,1)-real-interval
 */
/**
 * converts an unsigned 32-bit integer to a double on (0,1)-real-interval.
 * @param v 32-bit unsigned integer
 * @return double on (0,1)-real-interval
 */
/**
 * generates a random number on (0,1)-real-interval
 * @param sfmt SFMT internal state
 * @return double on (0,1)-real-interval
 */
/**
 * converts an unsigned 64-bit integer to double on [0,1)
 * with 53-bit resolution.
 * @param v 64-bit unsigned integer
 * @return double on [0,1)-real-interval with 53-bit resolution.
 */
/**
 * generates a random number on [0,1) with 53-bit resolution
 * @param sfmt SFMT internal state
 * @return double on [0,1) with 53-bit resolution
 */
/* =================================================
   The following function are added by Saito.
   ================================================= */
/**
 * generates a random number on [0,1) with 53-bit resolution from two
 * 32 bit integers
 */
/**
 * generates a random number on [0,1) with 53-bit resolution
 * using two 32bit integers.
 * @param sfmt SFMT internal state
 * @return double on [0,1) with 53-bit resolution
 */
#if defined(__cplusplus)
}
#endif

#endif
