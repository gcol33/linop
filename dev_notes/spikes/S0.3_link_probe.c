/* S0.3 link probe: real-double PRIMME path against R's BLAS/LAPACK. */
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "primme.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* 1-D Dirichlet Laplacian, tridiag(-1,2,-1): eigenvalues 4 sin^2(k pi/(2(n+1))) */
static void matvec(void *vx, PRIMME_INT *ldx, void *vy, PRIMME_INT *ldy,
                   int *blockSize, primme_params *primme, int *ierr) {
   double *x = (double *)vx, *y = (double *)vy;
   int n = (int)primme->n, b, i;
   for (b = 0; b < *blockSize; b++) {
      double *xv = x + b * (*ldx), *yv = y + b * (*ldy);
      for (i = 0; i < n; i++) {
         yv[i] = 2.0 * xv[i];
         if (i > 0)     yv[i] -= xv[i - 1];
         if (i < n - 1) yv[i] -= xv[i + 1];
      }
   }
   *ierr = 0;
}

int main(void) {
   primme_params primme;
   int n = 200, k = 4, i, ret;
   double *evals, *evecs, *rnorms;

   primme_initialize(&primme);
   primme.matrixMatvec = matvec;
   primme.n            = n;
   primme.numEvals     = k;
   primme.target       = primme_smallest;
   primme.eps          = 1e-10;
   primme_set_method(PRIMME_DEFAULT_MIN_TIME, &primme);

   evals  = (double *)malloc(k * sizeof(double));
   rnorms = (double *)malloc(k * sizeof(double));
   evecs  = (double *)malloc((size_t)n * k * sizeof(double));

   ret = dprimme(evals, evecs, rnorms, &primme);
   if (ret != 0) { printf("dprimme failed: %d\n", ret); return 1; }

   printf("dprimme OK, %d smallest eigenvalues of tridiag(-1,2,-1), n=%d\n", k, n);
   for (i = 0; i < k; i++) {
      double exact = 4.0 * pow(sin((i + 1) * M_PI / (2.0 * (n + 1))), 2.0);
      printf("  lambda[%d] = %.12f   exact = %.12f   |diff| = %.3e   rnorm = %.3e\n",
             i, evals[i], exact, fabs(evals[i] - exact), rnorms[i]);
   }
   printf("matvecs: %lld\n", (long long)primme.stats.numMatvecs);
   primme_free(&primme);
   free(evals); free(evecs); free(rnorms);
   return 0;
}
