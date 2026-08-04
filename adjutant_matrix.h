#ifndef __ADJUTANT_MATRIX_H__
#define __ADJUTANT_MATRIX_H__

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


// Wrappers for frequent BLAS/LAPACK calls (sometimes with no-BLAS/LAPACK
// fallbacks).
//
// Currently supports Accelerate, AOCL, ATLAS, MKL, and OpenBLAS backends.
//
// BLAS functions are referred to as cblas_<fname>(), LAPACK functions are
// referred to as LAPACK_<fname>(), and the integer index type is lapack_int.

#include "include/plink2_base.h"
#include "include/plink2_simd.h"

#ifdef NOLAPACK
typedef double MatrixInvertBuf1;
CONSTI32(kMatrixInvertBuf1ElemAlloc, 2 * sizeof(double));
CONSTI32(kMatrixInvertBuf1CheckedAlloc, 2 * sizeof(double));
#  define lapack_int int
#  define BLAS_SET_NUM_THREADS(num)
#else  // not NOLAPACK

// 1. Define LAPACK_ILP64 and USE_MKL consistently.
#  ifdef __APPLE__
// Make -DLAPACK_ILP64 and -DACCELERATE_LAPACK_ILP64 have the same effect.
#    if defined(LAPACK_ILP64) && !defined(ACCELERATE_LAPACK_ILP64)
#      define ACCELERATE_LAPACK_ILP64
#    endif
#    if defined(ACCELERATE_LAPACK_ILP64) && !defined(LAPACK_ILP64)
#      define LAPACK_ILP64
#    endif
#  endif

#  if defined(LAPACK_ILP64) && !defined(__LP64__)
#    error "Invalid compile flags."
#  endif

#  ifdef DYNAMIC_MKL
#    define USE_MKL
#  endif

#  ifdef __APPLE__
// 2a. Apple Accelerate.
#    ifdef USE_AOCL
#      error "plink2 cannot use AOCL on macOS."
#    endif
#    ifdef USE_MKL
#      error "plink2 cannot currently use MKL on macOS."
#    endif
#    ifdef USE_OPENBLAS
#      error "plink2 cannot currently use OpenBLAS on macOS."
#    endif

#    if defined(LAPACK_ILP64) && !defined(ACCELERATE_NEW_LAPACK)
#      error "LAPACK_ILP64 requires ACCELERATE_NEW_LAPACK on macOS"
#    endif
#    include <Accelerate/Accelerate.h>  // IWYU pragma: export
#    if defined(ACCELERATE_NEW_LAPACK)
#      define USE_MTBLAS
HEADER_INLINE void BLAS_SET_NUM_THREADS(__attribute__((unused)) int num_threads) {
  if (__builtin_available(macOS 15.0, *)) {
    if (num_threads > 1) {
      BLASSetThreading(BLAS_THREADING_MULTI_THREADED);
    } else {
      BLASSetThreading(BLAS_THREADING_SINGLE_THREADED);
    }
  }
}
#    else
#      define BLAS_SET_NUM_THREADS(num)

#    endif

#  elif defined(USE_MKL)
// 2b. Intel MKL.
#    if defined(USE_AOCL)
#      error "USE_AOCL and USE_MKL cannot both be defined."
#    endif
#    if defined(USE_OPENBLAS)
#      error "USE_MKL and USE_OPENBLAS cannot both be defined."
#    endif
#    ifdef _WIN32
#      error "plink2 cannot currently use MKL when compiled for Windows."
#    endif

#    ifdef LAPACK_ILP64
#      define MKL_ILP64
#    endif
#    ifdef DYNAMIC_MKL
#      include <mkl_service.h>  // IWYU pragma: export
#      include <mkl_cblas.h>  // IWYU pragma: export
#      include <mkl_lapack.h>  // IWYU pragma: export
#    else
// If this isn't initially found, use the compiler's -I option to specify the
// appropriate include-file directory.  As of this writing, it is expected to
// be something like /opt/intel/oneapi/mkl/<version>/include .
#      include "mkl_service.h"  // IWYU pragma: export
#      include "mkl_cblas.h"  // IWYU pragma: export
#      include "mkl_lapack.h"  // IWYU pragma: export
#    endif
#    define USE_MTBLAS
// This technically doesn't have to be a macro, but it's surrounded by other
// things which do have to be macros, so changing this to a namespaced function
// arguably *decreases* overall readability...
#    define BLAS_SET_NUM_THREADS mkl_set_num_threads

#  else
// 2c. AOCL, ATLAS, OpenBLAS.
#    if defined(USE_AOCL) && defined(USE_OPENBLAS)
#      error "USE_AOCL and USE_OPENBLAS cannot both be defined."
#    endif
#    if defined(__cplusplus) && !defined(__LP64__) && !defined(_WIN32)
// only subcase that seems to need this; win32 requires it *not* be here
extern "C" {
#    endif
#    include <cblas.h>  // IWYU pragma: export
#    include <lapacke.h>  // IWYU pragma: export
#    if defined(__cplusplus) && !defined(__LP64__) && !defined(_WIN32)
}  // extern "C"
#    endif
#    ifdef USE_OPENBLAS
#      define USE_MTBLAS
#      define BLAS_SET_NUM_THREADS openblas_set_num_threads
#    else
#      define BLAS_SET_NUM_THREADS(num)
#    endif
#  endif

// This will probably be moved the next time CUDA work is done.
#  ifdef USE_CUDA
#    include "cuda/plink2_matrix_cuda.h"
#  endif

// 3. Make lapack_int and LAPACK_<fname>() work everywhere.
#  ifndef lapack_int
#    ifdef LAPACK_ILP64
#      ifdef __APPLE__
// argh
#        define lapack_int long
#      elif defined(USE_MKL)
#        define lapack_int long long
#      else
#        define lapack_int int64_t
#      endif
#    else
#      define lapack_int int32_t
#    endif
#  endif

#  if defined(__APPLE__) || defined(USE_MKL)
#    define LAPACK_dpotri dpotri_
#  endif

typedef lapack_int MatrixInvertBuf1;
// need to be careful about >= 2^32?
CONSTI32(kMatrixInvertBuf1ElemAlloc, sizeof(lapack_int));
// invert_matrix_checked() usually requires a larger buffer
CONSTI32(kMatrixInvertBuf1CheckedAlloc, 2 * sizeof(lapack_int));

#endif  // !NOLAPACK

#ifdef __cplusplus
namespace plink2 {
#endif

static const double kMatrixSingularRcond = 1e-14;

// Returns -1 if no inf/nan found.
// May move this to a more central location if there are other users.


// Copies (C-order) lower-left to upper right.
void ReflectMatrix(uint32_t dim, double* matrix);



// If dim < stride, this zeroes out the trailing elements of each row.


// AddDVec and AddFVec execute main += arg, under the assumption that both
// vectors have been zero-padded.
// "ctav" is the float-length of both vectors, rounded up to a vector boundary.
// (With the current implementation, it's harmless to just pass in ct instead,
// but ctav is better practice since it indicates awareness of what this
// function is actually doing.)
#ifdef __LP64__
// FillDVec sets all elements of dst to dxx, and zero-fills any trailing
// elements (w.r.t. vector boundaries).
// bugfix (21 Nov 2023): forgot about trailing zero-fill.
#else
#endif

#ifdef __LP64__
// FillFVec sets all elements of dst to fxx, and zero-fills any trailing
// elements (w.r.t. vector boundaries).
#else
#endif
// todo: benchmark again, last round was before Spectre/Meltdown mitigations
CONSTI32(kDotprodDThresh, 17);
CONSTI32(kDotprodFThresh, 15);

#ifdef NOLAPACK
BoolErr InvertMatrix(int32_t dim, double* matrix, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf);
HEADER_INLINE BoolErr InvertSymmdefMatrixChecked(int32_t dim, double* matrix, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf) {
  ReflectMatrix(dim, matrix);
  return InvertMatrix(dim, matrix, dbl_1d_buf, dbl_2d_buf);
}
// first half computes either LU or singular value decomposition, and
//   determinant
// second half actually inverts matrix, assuming 1d_buf and 2d_buf have results
//   from first half


#else // !NOLAPACK
// not worthwhile for ct < 16.
// extra if-statement in DotprodxF() seems disproportionally expensive in
// test?... guess I won't have auto-branch for now.

BoolErr InvertMatrix(lapack_int dim, double* matrix, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf);



// InvertSymmdef... functions only assume (C-order) lower left of matrix is
// filled, and only the lower left of the return matrix is valid.

BoolErr InvertSymmdefStridedMatrix(lapack_int dim, lapack_int stride, double* matrix, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf);

// dbl_2d_buf must have room for at least max(dim, 3) * dim elements.
BoolErr InvertSymmdefMatrixChecked(lapack_int dim, double* matrix, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf);

BoolErr InvertSymmdefMatrixFirstHalf(lapack_int dim, lapack_int stride, double* matrix, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf);




/*
HEADER_INLINE double HalfInvertedDet(__maybe_unused const double* half_inverted_iter, __maybe_unused const MatrixInvertBuf1* int_1d_buf, uint32_t dim) {
  uint32_t dim_p1 = dim + 1;
  double det_u = *half_inverted_iter;
  for (uint32_t uii = 1; uii != dim; ++uii) {
    half_inverted_iter = &(half_inverted_iter[dim_p1]);
    det_u *= (*half_inverted_iter);
  }
  return fabs(det_u);
}
*/

#endif // !NOLAPACK
// this is essentially a full-blown dgemm wrapper, only missing the alpha
// parameter now
void ColMajorMatrixMultiplyStridedAddassign(const double* inmatrix1, const double* inmatrix2, lapack_int row1_ct, lapack_int stride1, lapack_int col2_ct, lapack_int stride2, lapack_int common_ct, lapack_int stride3, double beta, double* outmatrix);

HEADER_INLINE void RowMajorMatrixMultiplyStridedIncr(const double* inmatrix1, const double* inmatrix2, lapack_int row1_ct, lapack_int stride1, lapack_int col2_ct, lapack_int stride2, lapack_int common_ct, lapack_int stride3, double* outmatrix) {
  return ColMajorMatrixMultiplyStridedAddassign(inmatrix2, inmatrix1, col2_ct, stride2, row1_ct, stride1, common_ct, stride3, 1.0, outmatrix);
}

// out^T := V^T * M
// out := M^T * V

// out := M * V

// out := M * V

// out^T := V^T * M
// out := M^T * V




// A(A^T), where A is row-major; result is dim x dim
// ONLY UPDATES LOWER TRIANGLE OF result[].




// (A^T)A

#ifndef NOLAPACK

// currently a wrapper for LAPACK_dgesvd().

// currently a wrapper for LAPACK_dsyevr().  Matrix is expected to be
// lower-triangular from C's perspective.
// reverse_eigvecs is eigenvector-major, but the vectors are in order of
// *increasing* eigenvalue.
#endif

// Computes inverse of
//   [ A   b^T ]
//   [ b   c   ]
// given precomputed A^{-1} (must be fully filled out, not just lower left).
// See e.g.
//   https://en.wikipedia.org/wiki/Invertible_matrix#Blockwise_inversion .
// Only fills lower left of outmatrix.
// insert_idx specifies the zero-based row/column number of b/c in outmatrix.

// When you only need the diagonal from InvertRank1Symm().  insert_idx
// assumed to be orig_dim.

//   [ A   B^T ]
//   [ B   D   ]
// B is row-major.


#ifdef NOLAPACK
#else
#endif

// now assumes xtx_inv is predictors_pmaj * transpose on input
// just for debugging
#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __ADJUTANT_MATRIX_H__
