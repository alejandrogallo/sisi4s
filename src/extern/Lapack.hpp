/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LAPACK_DEFINED
#define LAPACK_DEFINED

#include <math/Float.hpp>
#include <math/Complex.hpp>

// TODO: use define for name mangling: underscore or not

extern "C" {
  void dgeev_(
    const char *jobvLeft,
    const char *jobvRight,
    const int *n,
    cc4s::Float64 *a,
    const int *lda,
    cc4s::Float64 *wReal,
    cc4s::Float64 *wImag,
    cc4s::Float64 *vLeft,
    const int *ldvLeft,
    cc4s::Float64 *vRight,
    const int *ldvRight,
    cc4s::Float64 *work,
    const int *workSize,
    int *info
  );
  void zgeev_(
    const char *jobvLeft,
    const char *jobvRight,
    const int *n,
    const cc4s::Complex64 *a,
    const int *lda,
    cc4s::Complex64 *w,
    cc4s::Complex64 *vLeft,
    const int *ldvLeft,
    cc4s::Complex64 *vRight,
    const int *ldvRight,
    cc4s::Complex64 *work,
    const int *workSize,
    cc4s::Float64 *rwork,
    int *info
  );
  void dsyev_(
    const char *jobz,
    const char *uplo,
    const int *n,
    cc4s::Float64 *a,
    const int *lda,
    cc4s::Float64 *w,
    cc4s::Float64 *work,
    const int *lwork,
    int *info
  );
  void dsysv_(
    const char *uplo,
    const int *n,
    const int *m,
    cc4s::Float64 *a,
    const int *lda,
    const int *ipiv,
    cc4s::Float64 *b,
    const int *ldb,
    cc4s::Float64 *work,
    const int *lwork,
    const int *info
  );
  void dgetrf_(
    const int *m,
    const int *n,
    cc4s::Float64 *a,
    const int *lda,
    const int *rowPermutation,
    int *info
  );
  void zgetrf_(
    const int *m,
    const int *n,
    cc4s::Complex64 *a,
    const int *lda,
    const int *rowPermutation,
    int *info
  );
  void dgetri_(
    const int *n,
    cc4s::Float64 *a,
    const int *lda,
    const int *rowPermutation,
    cc4s::Float64 *work,
    const int *workSize,
    int *info
  );
  void zgetri_(
    const int *n,
    cc4s::Complex64 *a,
    const int *lda,
    const int *rowPermutation,
    cc4s::Complex64 *work,
    const int *workSize,
    int *info
  );
  void dgemm_(
    const char *transa,
    const char *transb,
    const int *m,
    const int *n,
    const int *k,
    double *alpha,
    const double *A,
    const int *lda,
    const double *B,
    const int *ldb,
    double *beta,
    double *C,
    const int *ldc
  );
  void dger_(
    const int *M,
    const int *N,
    const double *alpha,
    const double *X,
    const int *incx,
    const double *Y,
    const int *incy,
    double *A,
    const int *lda
  );

};

#endif

