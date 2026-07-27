/* S0.3 proof of concept: supply zheevx_ / zhegvx_ on top of routines R does
 * export (zheevd_, zpotrf_, ztrsm_, zhemm_), so a vendored PRIMME links against
 * R's Rblas/Rlapack without an external LAPACK.
 *
 * PRIMME uses these only for the SMALL dense projected eigenproblem, whose order
 * is the Krylov basis size, so computing the whole spectrum and then selecting is
 * cheap. Correctness, not asymptotics, is the constraint here.
 *
 * Spike code. Not shipped. Range selection implements RANGE = 'A', 'V', 'I'.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct { double r, i; } zcplx;

extern void zheevd_(const char *jobz, const char *uplo, const int *n, zcplx *a,
                    const int *lda, double *w, zcplx *work, const int *lwork,
                    double *rwork, const int *lrwork, int *iwork,
                    const int *liwork, int *info);
extern void zpotrf_(const char *uplo, const int *n, zcplx *a, const int *lda, int *info);
extern void ztrsm_(const char *side, const char *uplo, const char *transa,
                   const char *diag, const int *m, const int *n, const zcplx *alpha,
                   const zcplx *a, const int *lda, zcplx *b, const int *ldb);

static int lsame(const char *a, char b) {
   char c = *a;
   if (c >= 'a' && c <= 'z') c = (char)(c - 'a' + 'A');
   return c == b;
}

void zheevx_(const char *jobz, const char *range, const char *uplo, const int *n,
             zcplx *a, const int *lda, const double *vl, const double *vu,
             const int *il, const int *iu, const double *abstol, int *m,
             double *w, zcplx *z, const int *ldz, zcplx *work, const int *lwork,
             double *rwork, int *iwork, int *ifail, int *info)
{
   int N = *n, i, j, lo, hi, cnt, wantz = lsame(jobz, 'V');
   int lw = -1, lrw = -1, liw = -1, i1;
   zcplx wq; double rwq; int iwq;
   zcplx *A = NULL, *W2 = NULL; double *ww = NULL, *rw = NULL; int *iw = NULL;

   (void)abstol; (void)lwork; (void)work; (void)rwork; (void)iwork;
   *info = 0;

   /* work on a copy: zheevd overwrites A with the eigenvectors */
   A = (zcplx *)malloc((size_t)N * N * sizeof(zcplx));
   ww = (double *)malloc((size_t)N * sizeof(double));
   if (!A || !ww) { *info = -1; free(A); free(ww); return; }
   for (j = 0; j < N; j++)
      for (i = 0; i < N; i++) A[i + (size_t)j * N] = a[i + (size_t)j * (*lda)];

   /* workspace query, then allocate */
   zheevd_(wantz ? "V" : "N", uplo, &N, A, &N, ww, &wq, &lw, &rwq, &lrw, &iwq, &liw, info);
   if (*info != 0) { free(A); free(ww); return; }
   lw = (int)wq.r; lrw = (int)rwq; liw = iwq;
   W2 = (zcplx *)malloc((size_t)lw * sizeof(zcplx));
   rw = (double *)malloc((size_t)lrw * sizeof(double));
   iw = (int *)malloc((size_t)liw * sizeof(int));
   if (!W2 || !rw || !iw) { *info = -1; goto done; }

   zheevd_(wantz ? "V" : "N", uplo, &N, A, &N, ww, W2, &lw, rw, &lrw, iw, &liw, info);
   if (*info != 0) goto done;

   /* select the requested range; zheevd returns w in ascending order */
   lo = 0; hi = N - 1;
   if (lsame(range, 'I')) { lo = *il - 1; hi = *iu - 1; }
   else if (lsame(range, 'V')) {
      lo = N; hi = -1;
      for (i = 0; i < N; i++)
         if (ww[i] > *vl && ww[i] <= *vu) { if (i < lo) lo = i; if (i > hi) hi = i; }
      if (hi < lo) { lo = 0; hi = -1; }
   }
   cnt = hi - lo + 1;
   if (cnt < 0) cnt = 0;
   *m = cnt;
   for (i = 0; i < cnt; i++) w[i] = ww[lo + i];
   if (wantz)
      for (j = 0; j < cnt; j++)
         for (i = 0; i < N; i++) z[i + (size_t)j * (*ldz)] = A[i + (size_t)(lo + j) * N];
   if (ifail) for (i1 = 0; i1 < cnt; i1++) ifail[i1] = 0;

done:
   free(A); free(ww); free(W2); free(rw); free(iw);
}

/* Generalized: reduce A x = lambda B x to standard form via Cholesky of B. */
void zhegvx_(const int *itype, const char *jobz, const char *range, const char *uplo,
             const int *n, zcplx *a, const int *lda, zcplx *b, const int *ldb,
             const double *vl, const double *vu, const int *il, const int *iu,
             const double *abstol, int *m, double *w, zcplx *z, const int *ldz,
             zcplx *work, const int *lwork, double *rwork, int *iwork,
             int *ifail, int *info)
{
   int N = *n, i, j;
   const zcplx one = {1.0, 0.0};
   zcplx *B = (zcplx *)malloc((size_t)N * N * sizeof(zcplx));
   zcplx *A = (zcplx *)malloc((size_t)N * N * sizeof(zcplx));
   if (!A || !B) { *info = -1; free(A); free(B); return; }

   for (j = 0; j < N; j++) for (i = 0; i < N; i++) {
      B[i + (size_t)j * N] = b[i + (size_t)j * (*ldb)];
      A[i + (size_t)j * N] = a[i + (size_t)j * (*lda)];
   }
   zpotrf_(uplo, &N, B, &N, info);
   if (*info != 0) { free(A); free(B); return; }

   /* itype 1: C = L^-1 A L^-H  (upper: U^-H A U^-1). Only itype 1 is used here. */
   (void)itype;
   if (lsame(uplo, 'U')) {
      ztrsm_("L", "U", "C", "N", &N, &N, &one, B, &N, A, &N);
      ztrsm_("R", "U", "N", "N", &N, &N, &one, B, &N, A, &N);
   } else {
      ztrsm_("R", "L", "C", "N", &N, &N, &one, B, &N, A, &N);
      ztrsm_("L", "L", "N", "N", &N, &N, &one, B, &N, A, &N);
   }

   zheevx_(jobz, range, uplo, &N, A, &N, vl, vu, il, iu, abstol, m, w, z, ldz,
           work, lwork, rwork, iwork, ifail, info);
   if (*info != 0) { free(A); free(B); return; }

   /* back-transform eigenvectors: x = L^-H y  (upper: U^-1 y) */
   if (lsame(jobz, 'V') && *m > 0) {
      if (lsame(uplo, 'U')) ztrsm_("L", "U", "N", "N", &N, m, &one, B, &N, z, ldz);
      else                  ztrsm_("L", "L", "C", "N", &N, m, &one, B, &N, z, ldz);
   }
   free(A); free(B);
}
