/*
 *  fmm.c
 *    routines related to fmm-style calculations
 *
 *    for piecewise constant elements
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gkGlobal.h"
#include "gk.h"

#define STOREM2L 0
#define SETUPONLY 0

/* blas: matrix times vector */
void dgemv_(char *tr, int *m, int *n, double *alpha, double *A, int *lda,
          double *x, int *incx, double *beta, double *y, int *incy);
void dgetrf_(int *n, int *m, double *A, int *lda, int *ipiv, int *incy);
void dgetrs_(char *tr, int *n, int *nrhs, double *A, int *lda, int *ipiv,
          double *b, int *lbd, int *incy);

void kernelKER4( double *x, double *y);
double *panelIA0(panel *pnlX, panel *pnlY );

void setupDerivs(int order, double *x );

extern int normErr;
extern void (*kernel)();
extern double kappa;
extern double epsilon;

double **Gp;            /* translation matrices */
double **Q2PK1, **Q2PK2, **Q2PK3, **Q2PK4;
double **Q2M0, **Q2M1;
double **L2P0, **L2P1;
double *matrixA;

/*
 * direct solver, using LU method
 */
void directSolver( ssystem *sys, double *sgm ) {
  panel *pnlX, *pnlY;
  double *KER;
  int i, j, ii, nPnls=sys->nPnls, n, one, *ipiv, inc=1;

  kernel = kernelKER4;

  CALLOC(matrixA, 4*nPnls*nPnls, double, ON, AQ2P);
  CALLOC(ipiv, 2*nPnls, int, ON, AQ2P);

  for ( i=0,pnlX=sys->pnlLst; i<nPnls; i++,pnlX=pnlX->nextC ) {
    for ( j=0,pnlY=sys->pnlLst; j<nPnls; j++,pnlY=pnlY->nextC ) {
      KER = panelIA0(pnlX, pnlY);
//      printf("%f %f\n",KER[0],KER[3]);
      matrixA[i*2*nPnls+j] = -KER[1];
      matrixA[i*2*nPnls+j+nPnls] = -KER[0];
      matrixA[(i+nPnls)*2*nPnls+j] = -KER[3];
      matrixA[(i+nPnls)*2*nPnls+j+nPnls] = -KER[2];
//      printf("%f %f %f %f\n",KER[1],KER[0],KER[3],KER[2]);
//      printf("%f %f\n",matrixA[(j+nPnls)*2*nPnls+i],matrixA[(j+nPnls)*2*nPnls+i+nPnls]);
    }
    //exit(0);
  }

  for ( i=0,pnlX=sys->pnlLst; i<nPnls; i++,pnlX=pnlX->nextC ) {
    matrixA[i*2*nPnls+i] += (1.0+epsilon)/2.0*pnlX->area;
    matrixA[(i+nPnls)*2*nPnls+i+nPnls] += (1.0+1.0/epsilon)/2.0*pnlX->area;
    //matrixA[i*2*nPnls+i] += (1.0+epsilon)/2.0;
    //matrixA[(i+nPnls)*2*nPnls+i+nPnls] += (1.0+1.0/epsilon)/2.0;

    //ipiv[i] = i;
    //ipiv[i+nPnls] = i+nPnls;
  }

  n = 2*nPnls;
  one = 1;

  dgetrf_( &n, &n, matrixA, &n, ipiv, &inc );
  /* hChr is the transpose, since the dgetrs stores in column */
  dgetrs_( &hChr, &n, &one, matrixA, &n, ipiv, sgm, &n, &inc );

}/* directSolver */

/*
 *  direct sum for test
 */
void directSum( ssystem *sys ) {
  panel *pnlX, *pnlY;
  double *KER;
  kernel = kernelKER4;

  for ( pnlX=sys->pnlLst; pnlX!=NULL; pnlX=pnlX->nextC ) {
    for ( pnlY=sys->pnlLst; pnlY!=NULL; pnlY=pnlY->nextC ) {
      //pnlY=pnlX->next;
      KER = panelIA0(pnlX, pnlY);
      printf("%f %f\n",KER[0],KER[1]);
    }
    exit(0);
  }
}/* directSum */

/*
 *  compute the matrices (piecewise constant elements) for direct sum
 */
void setupDirect( ssystem *sys ) {
  panel *pnlX, *pnlY;
  cube *cb1, *cb;
  int nPnls=sys->nPnls;
  int i, j, k, ii, inbr, idx;
  double *KER;

  /* set up kernel */
  kernel = kernelKER4;

  /* allocate four terms for Coulomb and screened Coulomb */
/*
  CALLOC(Q2PK1, nPnls, double*, ON, AQ2P);
  CALLOC(Q2PK2, nPnls, double*, ON, AQ2P);
  CALLOC(Q2PK3, nPnls, double*, ON, AQ2P);
  CALLOC(Q2PK4, nPnls, double*, ON, AQ2P);
  for ( i=0; i<nPnls; i++ ) {
    CALLOC(Q2PK1[i], nPnls, double, ON, AQ2P);
    CALLOC(Q2PK2[i], nPnls, double, ON, AQ2P);
    CALLOC(Q2PK3[i], nPnls, double, ON, AQ2P);
    CALLOC(Q2PK4[i], nPnls, double, ON, AQ2P);
  }

  for ( i=0,pnlX=sys->pnlOLst; pnlX!=NULL; pnlX=pnlX->next,i++ ) {
    for ( j=0,pnlY=sys->pnlOLst; pnlY!=NULL; pnlY=pnlY->next,j++ ) {
      KER = panelIA0(pnlX, pnlY);
      Q2PK1[i][j] = KER[0];
      Q2PK2[i][j] = KER[1];
      Q2PK3[i][j] = KER[2];
      Q2PK4[i][j] = KER[3];
    }
  }
*/

  CALLOC(matrixA, 4*nPnls*nPnls, double, ON, AQ2P);

  for ( i=0,pnlX=sys->pnlLst; i<nPnls; i++,pnlX=pnlX->nextC ) {
    for ( j=0,pnlY=sys->pnlLst; j<nPnls; j++,pnlY=pnlY->nextC ) {
      KER = panelIA0(pnlX, pnlY);
      matrixA[i*2*nPnls+j] = -KER[1];
      matrixA[i*2*nPnls+j+nPnls] = -KER[0];
      matrixA[(i+nPnls)*2*nPnls+j] = -KER[3];
      matrixA[(i+nPnls)*2*nPnls+j+nPnls] = -KER[2];
    }
//    exit(0);
  }

  for ( i=0,pnlX=sys->pnlLst; i<nPnls; i++,pnlX=pnlX->nextC ) {
    matrixA[i*2*nPnls+i] += (1.0+epsilon)/2.0*pnlX->area;
    matrixA[(i+nPnls)*2*nPnls+i+nPnls] += (1.0+1.0/epsilon)/2.0*pnlX->area;
    //matrixA[i*2*nPnls+i] += (1.0+epsilon)/2.0;
    //matrixA[(i+nPnls)*2*nPnls+i+nPnls] += (1.0+1.0/epsilon)/2.0;
  }
} /* setupDirect */

/*
 * Add the nearfield.
 * Pcw constant case. The nearfield coefficients are computed
 * by setupNearfield0().
 * Same parameters as applyFMM().
 */
void applyDirect(ssystem *sys, double *sgm, double *pot) {
  panel *pnlX, *pnlY;
  double *KER;
  double *x_pot, *y_pot, *x_dpdn, *y_dpdn;
  double tp, tpd;
  int nPnls=sys->nPnls;
  int i, j, n, inc=1;

  n = 2*nPnls;

//  dgemv_(&hChr, &n, &n, &one, matrixA, &n, sgm, &inc, &zero, pot, &inc);

  x_pot = &(sgm[0]);
  x_dpdn = &(sgm[nPnls]);
  y_pot = &(pot[0]);
  y_dpdn = &(pot[nPnls]);

  kernel = kernelKER4;

  for ( i=0,pnlX=sys->pnlLst; pnlX!=NULL; pnlX=pnlX->nextC,i++ ) {
    tp = 0.0;
    tpd = 0.0;
    for ( j=0,pnlY=sys->pnlLst; pnlY!=NULL; pnlY=pnlY->nextC,j++ ) {
      KER = panelIA0(pnlX, pnlY);
      tp += KER[0]*x_dpdn[j]+KER[1]*x_pot[j];
      tpd += KER[2]*x_dpdn[j]+KER[3]*x_pot[j];
//      printf("%.12e %.12e %.12e %.12e\n",KER[1],KER[0],KER[3],KER[2]);
//      printf("%.12e %.12e %.12e %.12e\n",KER[0]*x_dpdn[j],KER[1]*x_pot[j],KER[2]*x_dpdn[j],KER[3]*x_pot[j]);
//      printf("%.12e %.12e\n",KER[2]*x_dpdn[j],KER[3]*x_pot[j]);
//      printf("%.12e %.12e\n",tp,x_dpdn[j]);

    }
//    exit(0);
    //y_pot[i] = (1.0+epsilon)/2.0*x_pot[i] - tp;
    //y_dpdn[i] = (1.0+1.0/epsilon)/2.0*x_dpdn[i] - tpd;
    y_pot[i] = (1.0+epsilon)/2.0*x_pot[i]*pnlX->area - tp;
    y_dpdn[i] = (1.0+1.0/epsilon)/2.0*x_dpdn[i]*pnlX->area - tpd;
//    printf("%f %f\n",tp,tpd);
  }
//  exit(0);


} /* applyNearfield0 */
