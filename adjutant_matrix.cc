// This file is part of PLINK 2.0, copyright (C) 2005-2022 Shaun Purcell,
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

#include "adjutant_matrix.h"

#include <assert.h>
#include <math.h>  // fabs(), sqrt()
#include <string.h>

#include "include/plink2_float.h"

#ifndef NOLAPACK
#  if defined(__APPLE__) || defined(USE_MKL)
#    define LAPACK_dgecon dgecon_
#    define LAPACK_dgesvd dgesvd_
#    define LAPACK_dgetrf dgetrf_
#    define LAPACK_dgetri dgetri_
#    define LAPACK_dlange dlange_
#    define LAPACK_dlansy dlansy_
#    define LAPACK_dpocon dpocon_
#    define LAPACK_dpotrf dpotrf_
#    define LAPACK_dpotrs dpotrs_
#    define LAPACK_dsyevr dsyevr_
#  endif
#endif  // !NOLAPACK


#ifdef __cplusplus
namespace plink2 {
#endif
void ReflectMatrix(uint32_t dim, double* matrix) {
  const uintptr_t dim_p1l = dim + 1;
  double* write_row = matrix;
  for (uint32_t row_idx = 0; row_idx != dim; ++row_idx) {
    double* read_col_iter = &(matrix[dim_p1l * row_idx + dim]);
    for (uint32_t col_idx = row_idx + 1; col_idx != dim; ++col_idx) {
      write_row[col_idx] = *read_col_iter;
      read_col_iter = &(read_col_iter[dim]);
    }
    write_row = &(write_row[dim]);
  }
}
static inline double Sqr(const double aa) {
  return aa * aa;
}

double Pythag(const double a, const double b) {
  // PLINK stats.cpp Pythag().
  double absa,absb;

  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) {
    return absa*sqrt(1.0+Sqr(absb/absa));
  }
  return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+Sqr(absa/absb)));
}

#ifdef NOLAPACK
static inline double MultiplyBySgn(double aa, double bb) {
  // PLINK helper.h MultiplyBySgn() template specialized to doubles.
  return (bb >= 0)? (aa >= 0 ? aa : -aa) : (aa >= 0 ? -aa : aa);
}

uint32_t SvdcmpC(int32_t m, double* a, double* w, double* v) {
  // C port of PLINK stats.cpp svdcmp().
  //   m: (in) number of rows/columns
  //   a: (in/out) input matrix, becomes first part of decomposition
  //   w: (out) singular values
  //   v: (out) final part of decomposition
  //
  // Now thread-safe.
  double* rv1 = &(w[S_CAST(uint32_t, m)]);
  int32_t n = m;
  int32_t flag;
  int32_t l = 0;  // suppress compile warning
  int32_t i,its,j,jj,k,nm;
  double anorm,c,f,g,h,s,scale,x,y,z;
  double temp;

  g=scale=anorm=0.0;
  for (i=0;i!=n;++i) {
    l=i+2;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for (k=i;k!=m;++k) scale += fabs(a[k * m + i]);
      if (scale != 0.0) {
        for (k=i;k!=m;++k) {
          a[k * m + i] /= scale;
          s += a[k * m + i]*a[k * m + i];
        }
        f=a[i * m + i];
        g = -MultiplyBySgn(sqrt(s),f);
        h=f*g-s;
        a[i * m + i]=f-g;
        for (j=l-1;j!=n;++j) {
          for (s=0.0,k=i;k!=m;++k) s += a[k * m + i]*a[k * m + j];
          f=s/h;
          for (k=i;k!=m;++k) a[k * m + j] += f*a[k * m + i];
        }
        for (k=i;k!=m;++k) a[k * m + i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i+1 <= m && i+1 != n) {
      for (k=l-1;k!=n;++k) scale += fabs(a[i * m + k]);
      if (scale != 0.0) {
        for (k=l-1;k!=n;++k) {
          a[i * m + k] /= scale;
          s += a[i * m + k]*a[i * m + k];
        }
        f=a[i * m + l-1];
        g = -MultiplyBySgn(sqrt(s),f);
        h=f*g-s;
        a[i * m + l-1]=f-g;
        for (k=l-1;k!=n;++k) rv1[k]=a[i * m + k]/h;
        for (j=l-1;j!=m;++j) {
          for (s=0.0,k=l-1;k!=n;++k) s += a[j * m + k]*a[i * m + k];
          for (k=l-1;k!=n;++k) a[j * m + k] += s*rv1[k];
        }
        for (k=l-1;k!=n;++k) a[i * m + k] *= scale;
      }
    }
    anorm=MAXV(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n-1;i>=0;i--) {
    if (i < n-1) {
      if (g != 0.0) {
        for (j=l;j!=n;++j)
          v[j * m + i]=(a[i * m + j]/a[i * m + l])/g;
        for (j=l;j!=n;++j) {
          for (s=0.0,k=l;k!=n;++k) s += a[i * m + k]*v[k * m + j];
          for (k=l;k!=n;++k) v[k * m + j] += s*v[k * m + i];
        }
      }
      for (j=l;j!=n;++j) v[i * m + j]=v[j * m + i]=0.0;
    }
    v[i * m + i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=MINV(m,n)-1;i>=0;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j!=n;++j) a[i * m + j]=0.0;
    if (g != 0.0) {
      g=1.0/g;
      for (j=l;j!=n;++j) {
        for (s=0.0,k=l;k<m;++k) s += a[k * m + i]*a[k * m + j];
        f=(s/a[i * m + i])*g;
        for (k=i;k!=m;++k) a[k * m + j] += f*a[k * m + i];
      }
      for (j=i;j!=m;++j) a[j * m + i] *= g;
    } else for (j=i;j!=m;++j) a[j * m + i]=0.0;
    ++a[i * m + i];
  }
  for (k=n-1;k>=0;k--) {
    for (its=0;its!=30;++its) {
      flag=1;
      nm=0;  // maybe-uninitialized warning
      for (l=k;l>=0;l--) {
        nm=l-1;
        temp=fabs(rv1[l])+anorm;
        if (temp == anorm) {
          flag=0;
          break;
        }
        temp=fabs(w[nm])+anorm;
        if (temp == anorm) break;
      }
      if (flag) {
        c=0.0;
        s=1.0;
        for (i=l;i!=k+1;++i) {
          f=s*rv1[i];
          rv1[i]=c*rv1[i];
          temp = fabs(f)+anorm;
          if (temp == anorm) break;
          g=w[i];
          h=Pythag(f,g);
          w[i]=h;
          h=1.0/h;
          c=g*h;
          s = -f*h;
          for (j=0;j!=m;++j) {
            y=a[j * m + nm];
            z=a[j * m + i];
            a[j * m + nm]=y*c+z*s;
            a[j * m + i]=z*c-y*s;
          }
        }
      }
      z=w[k];
      if (l == k) {
        if (z < 0.0) {
          w[k] = -z;
          for (j=0;j!=n;++j) v[j * m + k] = -v[j * m + k];
        }
        break;
      }
      if (its == 29)
        return 0;  // cannot converge: multi-collinearity?
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=Pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+MultiplyBySgn(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
        i=j+1;
        g=rv1[i];
        y=w[i];
        h=s*g;
        g=c*g;
        z=Pythag(f,h);
        rv1[j]=z;
        c=f/z;
        s=h/z;
        f=x*c+g*s;
        g=g*c-x*s;
        h=y*s;
        y *= c;
        for (jj=0;jj!=n;++jj) {
          x=v[jj * m + j];
          z=v[jj * m + i];
          v[jj * m + j]=x*c+z*s;
          v[jj * m + i]=z*c-x*s;
        }
        z=Pythag(f,h);
        w[j]=z;
        if (z != 0.0) {
          z=1.0/z;
          c=f*z;
          s=h*z;
        }
        f=c*g+s*y;
        x=c*y-s*g;
        for (jj=0;jj!=m;++jj) {
          y=a[jj * m + j];
          z=a[jj * m + i];
          a[jj * m + j]=y*c+z*s;
          a[jj * m + i]=z*c-y*s;
        }
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  return 1;
}

BoolErr InvertMatrix(int32_t dim, double* matrix, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf) {
  // C port of PLINK stats.cpp's svd_inverse() function.
  assert(dim > 0);

  // w -> dbl_1d_buf
  // v -> dbl_2d_buf
  const double eps = 1e-24;
  // not unlikely(), matrix inversion failure can be common
  if (!SvdcmpC(dim, matrix, dbl_1d_buf, dbl_2d_buf)) {
    return 1;
  }

  // Look for singular values
  double wmax = 0;
  int32_t i;
  for (i=0; i!=dim; ++i) {
    wmax = dbl_1d_buf[i] > wmax ? dbl_1d_buf[i] : wmax;
  }
  double wmin = wmax * eps;
  for (i=0; i!=dim; ++i) {
    dbl_1d_buf[i] = dbl_1d_buf[i] < wmin ? 0 : (1 / dbl_1d_buf[i]);
  }

  int32_t j;
  for (i=0; i!=dim; ++i) {
    for (j=0; j!=dim; ++j) {
      matrix[i * dim + j] = matrix[i * dim + j] * dbl_1d_buf[j];
    }
  }

  int32_t k;
  // [nxn].[t(v)]
  for (i=0; i!=dim; ++i) {
    ZeroDArr(dim, dbl_1d_buf);
    for (j=0; j!=dim; ++j) {
      for (k=0; k!=dim; ++k) {
        dbl_1d_buf[j] += matrix[i * dim + k] * dbl_2d_buf[j * dim + k];
      }
    }
    for (j = 0; j != dim; ++j) {
      matrix[i * dim + j] = dbl_1d_buf[j];
    }
  }
  for (i=1; i!=dim; ++i) {
    for(j=0; j!=i; ++j) {
      swap_f64(&(matrix[i * dim + j]), &(matrix[j * dim + i]));
    }
  }
  return 0;
}
#else  // !NOLAPACK
BoolErr InvertMatrix(lapack_int dim, double* matrix, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf) {
  // InvertSymmdefMatrix() is noticeably faster in the symmetric
  // positive-semidefinite case.
  lapack_int info;
  LAPACK_dgetrf(&dim, &dim, matrix, &dim, int_1d_buf, &info);
  if (info) {
    return 1;
  }
  lapack_int lwork = dim * dim;
  LAPACK_dgetri(&dim, matrix, &dim, int_1d_buf, dbl_2d_buf, &lwork, &info);
  assert(info == 0);
  return 0;
}
BoolErr InvertSymmdefStridedMatrix(lapack_int dim, lapack_int stride, double* matrix, __maybe_unused MatrixInvertBuf1* int_1d_buf, __maybe_unused double* dbl_2d_buf) {
  char uplo = 'U';
  lapack_int info;
  LAPACK_dpotrf(&uplo, &dim, matrix, &stride, &info);
  if (info) {
    return 1;
  }
  LAPACK_dpotri(&uplo, &dim, matrix, &stride, &info);
  assert(info == 0);
  return 0;
}

BoolErr InvertSymmdefMatrixChecked(lapack_int dim, double* matrix, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf) {
  char cc = '1';
  char uplo = 'U';
  double norm = LAPACK_dlansy(&cc, &uplo, &dim, matrix, &dim, dbl_2d_buf);
  lapack_int info;
  LAPACK_dpotrf(&uplo, &dim, matrix, &dim, &info);
  if (info > 0) {
    return 1;
  }
  double rcond;
  LAPACK_dpocon(&uplo, &dim, matrix, &dim, &norm, &rcond, dbl_2d_buf, int_1d_buf, &info);
  if (rcond < kMatrixSingularRcond) {
    return 1;
  }
  LAPACK_dpotri(&uplo, &dim, matrix, &dim, &info);
  return 0;
}

BoolErr InvertSymmdefMatrixFirstHalf(lapack_int dim, lapack_int stride, double* matrix, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf) {
  char cc = '1';
  char uplo = 'U';
  double norm = LAPACK_dlansy(&cc, &uplo, &dim, matrix, &stride, dbl_2d_buf);
  lapack_int info;
  LAPACK_dpotrf(&uplo, &dim, matrix, &stride, &info);
  if (info > 0) {
    return 1;
  }
  double rcond;
  LAPACK_dpocon(&uplo, &dim, matrix, &stride, &norm, &rcond, dbl_2d_buf, int_1d_buf, &info);
  return (rcond < kMatrixSingularRcond);
}

// quasi-bugfix (20 Sep 2017): give up on doing this with single-precision
// numbers.  Instead, convert to double-precision, then perform the usual
// inversion, then downcode back to single-precision.
#endif  // !NOLAPACK
void ColMajorMatrixMultiplyStridedAddassign(const double* inmatrix1, const double* inmatrix2, lapack_int row1_ct, lapack_int stride1, lapack_int col2_ct, lapack_int stride2, lapack_int common_ct, lapack_int stride3, double beta, double* outmatrix) {
  // stride1 should be close to row1_ct
  // stride2 should be close to common_ct
  // output matrix uses stride3, which should be close to row1_ct
#ifdef NOLAPACK
  const uintptr_t row1_ct_l = row1_ct;
  const uintptr_t col2_ct_l = col2_ct;
  const uintptr_t common_ct_l = common_ct;
  // not optimized, no beta == 0 special case
  for (uintptr_t col_idx = 0; col_idx != col2_ct_l; ++col_idx) {
    double* outmatrix_row_iter = &(outmatrix[col_idx * stride3]);
    for (uintptr_t row_idx = 0; row_idx != row1_ct_l; ++row_idx) {
      double cur_entry = 0.0;
      const double* col2_iter = &(inmatrix2[col_idx * stride2]);
      for (uintptr_t com_idx = 0; com_idx != common_ct_l; com_idx++) {
        cur_entry += (*col2_iter++) * inmatrix1[com_idx * stride1 + row_idx];
      }
      *outmatrix_row_iter = (*outmatrix_row_iter) * beta + cur_entry;
      ++outmatrix_row_iter;
    }
  }
#else
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, row1_ct, col2_ct, common_ct, 1.0, inmatrix1, stride1, inmatrix2, stride2, beta, outmatrix, stride3);
#endif  // !NOLAPACK
}
// er, should make this _addassign for consistency...
// Briefly experimented with trying to speed this up, didn't make any progress.
// Tried vDSP_mtransD() on macOS and mkl_domatcopy() on Linux, neither of those
// helped, so it's likely there's no progress to be made.
/*
*/

// A(A^T), where A is row-major; result is dim x dim
// ONLY UPDATES LOWER TRIANGLE OF result[].
#ifndef NOLAPACK
// dsyevr_ takes ~30% less time than dsyevd_ on OS X dev machine.  todo: retest
// for Linux 64-bit MKL.
#endif  // !NOLAPACK
// now assumes xtx_inv is predictors_pmaj * transpose on input
#ifdef NOLAPACK
#else
#endif  // !NOLAPACK

#ifdef __cplusplus
}  // namespace plink2
#endif
