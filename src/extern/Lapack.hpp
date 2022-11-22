/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LAPACK_DEFINED
#define LAPACK_DEFINED

#include <math/Float.hpp>
#include <math/Complex.hpp>

// TODO: use define for name mangling: underscore or not

extern "C" {
void dgeev_(const char *jobvLeft,
            const char *jobvRight,
            const int *n,
            sisi4s::Float64 *a,
            const int *lda,
            sisi4s::Float64 *wReal,
            sisi4s::Float64 *wImag,
            sisi4s::Float64 *vLeft,
            const int *ldvLeft,
            sisi4s::Float64 *vRight,
            const int *ldvRight,
            sisi4s::Float64 *work,
            const int *workSize,
            int *info);
void zgeev_(const char *jobvLeft,
            const char *jobvRight,
            const int *n,
            const sisi4s::Complex64 *a,
            const int *lda,
            sisi4s::Complex64 *w,
            sisi4s::Complex64 *vLeft,
            const int *ldvLeft,
            sisi4s::Complex64 *vRight,
            const int *ldvRight,
            sisi4s::Complex64 *work,
            const int *workSize,
            sisi4s::Float64 *rwork,
            int *info);
void dsyev_(const char *jobz,
            const char *uplo,
            const int *n,
            sisi4s::Float64 *a,
            const int *lda,
            sisi4s::Float64 *w,
            sisi4s::Float64 *work,
            const int *lwork,
            int *info);
void dsysv_(const char *uplo,
            const int *n,
            const int *m,
            sisi4s::Float64 *a,
            const int *lda,
            const int *ipiv,
            sisi4s::Float64 *b,
            const int *ldb,
            sisi4s::Float64 *work,
            const int *lwork,
            const int *info);
void dgetrf_(const int *m,
             const int *n,
             sisi4s::Float64 *a,
             const int *lda,
             const int *rowPermutation,
             int *info);
void zgetrf_(const int *m,
             const int *n,
             sisi4s::Complex64 *a,
             const int *lda,
             const int *rowPermutation,
             int *info);
void dgetri_(const int *n,
             sisi4s::Float64 *a,
             const int *lda,
             const int *rowPermutation,
             sisi4s::Float64 *work,
             const int *workSize,
             int *info);
void zgetri_(const int *n,
             sisi4s::Complex64 *a,
             const int *lda,
             const int *rowPermutation,
             sisi4s::Complex64 *work,
             const int *workSize,
             int *info);
// void dgemm_(
//   const char *transa,
//   const char *transb,
//   const int *m,
//   const int *n,
//   const int *k,
//   double *alpha,
//   const double *A,
//   const int *lda,
//   const double *B,
//   const int *ldb,
//   double *beta,
//   double *C,
//   const int *ldc
// );
void dger_(const int *M,
           const int *N,
           const double *alpha,
           const double *X,
           const int *incx,
           const double *Y,
           const int *incy,
           double *A,
           const int *lda);
}

#endif
