#ifndef GMRES_H
#define GMRES_H

typedef int (*GmresMatVecFn)(double *alpha, double *x, double *beta, double *y);
typedef int (*GmresPrecondFn)(double *x, double *b);

int gmres(int n, double *b, double *x, int restrt, double *work, int ldw,
          double *h, int ldh, int *iter, double *resid,
          GmresMatVecFn matvec, GmresPrecondFn psolve, int *info);

#endif
