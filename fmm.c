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

void kernelKER4( double *x, double *y);
double *panelIA0(panel *pnlX, panel *pnlY );

double *calcMoments0(ssystem *sys, int order, cube *cb, int job);
void transM2M(ssystem *sys, cube *cbIn, cube *cbOut);
void transM2L(ssystem *sys, double *G0, double *Gk, cube *cbIn, cube *cbOut );
void transL2L(ssystem *sys, cube *cbIn, cube *cbOut );
void kernelDC0( double r, int p, double *G );
void kernelDS0( double r, int p, double *G );
void setupDerivs(int order, double *x );

double curvature(double *x, double *h20, double *h11, double *h02);
double paramEllip( panel *pnl, double x, double y, double *r, double *nrm );

extern double **dG0;     /* workspace for setupDerivs */
extern double **dGk;     /* workspace for setupDerivs */
extern int normErr;
extern void (*kernel)(), (*kernelDC)(), (*kernelDS)();
extern double kappa;
extern double epsilon;

double **Gp;            /* translation matrices */
double **Q2PK1, **Q2PK2, **Q2PK3, **Q2PK4;
double **Q2M0, **Q2M1;
double **L2P0, **L2P1;

/*
 *  compute the nearfield matrices (piecewise constant elements)
 */
void setupNearfield0( ssystem *sys ) {
  panel *pnlX, *pnlY;
  cube *cb1, *cb;
  int nNbrs, nPnl, nPnl1;
  int i, j, k, ii, inbr, idx;
  double *KER, *KER1, *KER2, *KER3, *KER4;

  /* set up kernel */
  kernel = kernelKER4;

  /* allocate eight terms for Coulomb and screened Coulomb */
  for ( idx=0, cb=sys->cubeList[sys->depth]; cb != NULL; cb=cb->next ) {
    idx += cb->nNbrs;
  }
  CALLOC(Q2PK1, idx, double*, ON, AQ2P);
  CALLOC(Q2PK2, idx, double*, ON, AQ2P);
  CALLOC(Q2PK3, idx, double*, ON, AQ2P);
  CALLOC(Q2PK4, idx, double*, ON, AQ2P);

  for ( idx=0, cb=sys->cubeList[sys->depth]; cb!=NULL; cb=cb->next ) {
    nNbrs = cb->nNbrs;
    nPnl = cb->nPnls;
//    printf("nNbrs:%d nPnl:%d idx:%d\n",nNbrs, nPnl, idx);
    for ( inbr=0; inbr<nNbrs; inbr++ ) {
      cb1 = cb->nbrs[inbr];
      nPnl1 = cb1->nPnls;

      CALLOC(KER1, nPnl*nPnl1, double, ON, AQ2P);
      Q2PK1[idx] = KER1;
      CALLOC(KER2, nPnl*nPnl1, double, ON, AQ2P);
      Q2PK2[idx] = KER2;
      CALLOC(KER3, nPnl*nPnl1, double, ON, AQ2P);
      Q2PK3[idx] = KER3;
      CALLOC(KER4, nPnl*nPnl1, double, ON, AQ2P);
      Q2PK4[idx++] = KER4;
      for ( ii=j=0, pnlY=cb1->pnls; j<nPnl1; j++, pnlY=pnlY->nextC ) {
        for ( i=0, pnlX=cb->pnls; i<nPnl; i++, pnlX=pnlX->nextC, ii++ ) {
          KER = panelIA0(pnlX, pnlY);
//          printf("%d %d %d %d\n",ii,pnlX->nSurf,pnlY->nSurf,pnlY->idx);

          KER1[ii] = KER[0];
          KER2[ii] = KER[1];
          KER3[ii] = KER[2];
          KER4[ii] = KER[3];
        }
      }
    }
  }
} /* setupNearfield0 */



/*
 * setup everything related to FMM-style matrix-vector multiply
 * note that the orders used in the FMM are given by sys->ordM2L[]
 */
void setupFMM(ssystem *sys) {
  cube *cb, *cb1, *kid;
  int depth=sys->depth, height=sys->height;
  int lev, idx, i, j, k;
  int inbr, nNbrs;
  int order, nMoments, nCubesL;
  double r[3];
  time_t time1, time2;

  /*
   * allocate and compute Q2P's or selfterms
   */
  time1 = time(&time1);
  setupNearfield0(sys);
  time2 = time(&time2);
  setupQ2PTime = difftime(time2,time1);

  /*
   * allocate and calculate Q2M's and L2P's
   * moments, local expansion coefficients
   */
  time1 = time(&time1);

  for ( nCubesL=0, cb=sys->cubeList[depth]; cb != NULL; cb=cb->next ) {
    nCubesL++;
  }
  CALLOC(Q2M0, nCubesL, double*, ON, AQ2M);
  CALLOC(Q2M1, nCubesL, double*, ON, AQ2M);
  CALLOC(L2P0, nCubesL, double*, ON, AQ2M);
  CALLOC(L2P1, nCubesL, double*, ON, AQ2M);

  order=sys->ordMom[depth];

  /* moments depend on which layer operator we have */
  for ( idx=0, cb=sys->cubeList[depth]; cb!=NULL; cb=cb->next, idx++ ) {
    L2P0[idx] = Q2M0[idx] = calcMoments0(sys, order, cb, 0);
    L2P1[idx] = Q2M1[idx] = calcMoments0(sys, order, cb, 1);
  }

  time2 = time(&time2);
  setupQ2MTime = difftime(time2,time1);

  for ( lev=sys->depth; lev>=sys->height; lev-- ) {
    order = sys->ordMom[lev];
    nMoments  = sys->nMom[order];
    for ( cb=sys->cubeList[lev]; cb != NULL; cb=cb->next ) {
      CALLOC(cb->mom_pot, nMoments, double, ON, ACUBES);
      CALLOC(cb->mom_dpdn, nMoments, double, ON, ACUBES);
      CALLOC(cb->lec_k1, nMoments, double, ON, ACUBES);
      CALLOC(cb->lec_k2, nMoments, double, ON, ACUBES);
      CALLOC(cb->lec_k3, nMoments, double, ON, ACUBES);
      CALLOC(cb->lec_k4, nMoments, double, ON, ACUBES);
    }
  }

    kernelDC = kernelDC0;
    kernelDS = kernelDS0;

#if STOREM2L
  /*
   * allocate and compute M2L's
   */
  time1 = time(&time1);

  for ( idx=0, lev=depth; lev>=height; lev-- ) {
    for ( cb=sys->cubeList[lev]; cb != NULL; cb=cb->next ) {
      idx += cb->n2Nbrs - cb->nNbrs;
    }
  }
  if ( idx==0 ) {
    printf("\n*****Warning: No multipole acceleration for panels\n");
  }
  else {
    CALLOC(Gp, idx, double*, ON, AM2L);
    for ( idx=0, lev=depth; lev>=height; lev-- ) {
      order=sys->ordM2L[lev];
      nMoments=sys->nMom[order];
      for ( cb=sys->cubeList[lev]; cb != NULL; cb=cb->next ) {
        nNbrs = cb->nNbrs;
        for ( inbr=cb->n2Nbrs-1; inbr>=nNbrs; inbr--, idx++ ) {
          cb1 = cb->nbrs[inbr];
          for ( k=0; k<3; k++ ) r[k] = cb->x[k] - cb1->x[k];
          setupDerivs(order, r);
          CALLOC(Gp[idx], nMoments, double, ON, AM2L);
          memcpy(Gp[idx], dG[0], nMoments*sizeof(double));
        }
      }
    }
  }

  time2 = time(&time2);
  setupM2LTime = difftime(time2,time1);
#endif

#if SETUPONLY
  dumpStats(sys);
  exit(1);
#endif

} /* setupFMM */


/*
 * re-do the necessary operations of setupFMM in case we have alreay computed the
 * adjoint and now want to do the single-layer
 */
void setupAdjoint2Single(ssystem *sys){
  cube *cb;
  int depth=sys->depth;
  int lev, idx, i, j, k;
  int inbr, nNbrs;
  time_t time1, time2;

//  sys->layer = SINGLE;

  time1 = time(&time1);
  setupNearfield0(sys);
  time2 = time(&time2);
  setupQ2PTime += difftime(time2,time1);

  for ( idx=0, cb=sys->cubeList[depth]; cb!=NULL; cb=cb->next, idx++ ) {
//    L2P[idx] = Q2M[idx];
  }

}  /* setupAdjoint2Single */




/*
 * Add the nearfield.
 * Pcw constant case. The nearfield coefficients are computed
 * by setupNearfield0().
 * Same parameters as applyFMM().
 */
void applyNearfield0(ssystem *sys, double *sgm, double *pot) {
  cube *cb, *cb1;
  double *x_pot, *y_pot, *x_dpdn, *y_dpdn;
  int inbr, nNbrs, idx, nPnls=sys->nPnls;
  int i, k, n, n1, inc=1;

  for ( idx=0, cb=sys->cubeList[sys->depth]; cb!=NULL; cb=cb->next ) {
    y_pot = &(pot[cb->pnls->idx]);
    y_dpdn = &(pot[cb->pnls->idx+nPnls]);
    n = cb->nPnls;
    nNbrs = cb->nNbrs;
    for ( inbr=0; inbr<nNbrs; inbr++, idx++ ) {
      cb1 = cb->nbrs[inbr];
      n1 = cb1->nPnls;
      x_pot = &(sgm[cb1->pnls->idx]);
      x_dpdn = &(sgm[cb1->pnls->idx+nPnls]);
      dgemv_(&nChr, &n, &n1, &one, Q2PK1[idx], &n, x_dpdn, &inc, &one, y_pot, &inc);
      dgemv_(&nChr, &n, &n1, &one, Q2PK2[idx], &n, x_pot, &inc, &one, y_pot, &inc);
      dgemv_(&nChr, &n, &n1, &one, Q2PK3[idx], &n, x_dpdn, &inc, &one, y_dpdn, &inc);
      dgemv_(&nChr, &n, &n1, &one, Q2PK4[idx], &n, x_pot, &inc, &one, y_dpdn, &inc);
    }
  }
} /* applyNearfield0 */




/*
 * FMM-style matrix-vector multiply for panels
 * Parameters
 *    sgm      input density (function on panels)
 *    pot      output potential (function on panels)
 */
void applyFMM(ssystem *sys, double *sgm, double *pot) {
  cube *cb, *cb1;
  double *x, *y, *lec;
  double r[3], *self;
  int depth=sys->depth, height=sys->height, nPnls=sys->nPnls;
  int nKid, nKid1, iNbr, iPnl, nNbrs, idx, nMom, order;
  int i, k, lev, n, n1, inc = 1;

  /* zero out mom's and lec's */
  for ( lev=depth; lev>=height; lev-- ) {
    nMom  = sys->nMom[sys->ordMom[lev]];
    for ( cb=sys->cubeList[lev]; cb != NULL; cb=cb->next ) {
      for ( k=0; k<nMom;  k++ ) {
        cb->mom_pot[k] = 0.0;
        cb->mom_dpdn[k] = 0.0;
        cb->lec_k1[k] = 0.0;
        cb->lec_k2[k] = 0.0;
        cb->lec_k3[k] = 0.0;
        cb->lec_k4[k] = 0.0;
      }
    }
  }

  /* Q2M transformations */
  nMom = sys->nMom[sys->ordMom[depth]];
  for ( idx=0, cb=sys->cubeList[depth]; cb != NULL; cb=cb->next, idx++ ) {
    x = &(sgm[cb->pnls->idx]);
    n = cb->nPnls;
    y = cb->mom_pot;
    dgemv_(&nChr, &nMom, &n, &one, Q2M1[idx], &nMom, x, &inc, &one, y, &inc);
    x = &(sgm[cb->pnls->idx+nPnls]);
    y = cb->mom_dpdn;
    dgemv_(&nChr, &nMom, &n, &one, Q2M0[idx], &nMom, x, &inc, &one, y, &inc);

  }

  /* upward pass */
  for ( lev=depth-1; lev>=height; lev-- ) {
    for ( cb=sys->cubeList[lev]; cb != NULL; cb=cb->next ) {
      for ( nKid=0; nKid<cb->nKids; nKid++ ) {
        transM2M(sys, cb->kids[nKid], cb);
      }
    }
  }

  /* Interaction phase */
  for ( idx=0, lev=depth; lev>=height; lev-- ) {
    for ( cb=sys->cubeList[lev]; cb != NULL; cb=cb->next ) {
      nNbrs = cb->nNbrs;
#if !STOREM2L
      order=sys->ordM2L[lev];
#endif
      for ( iNbr=cb->n2Nbrs-1; iNbr>=nNbrs; iNbr--, idx++ ) {
#if STOREM2L
        transM2L(sys, Gp[idx], Gp[idx], cb->nbrs[iNbr], cb);
#else
        cb1 = cb->nbrs[iNbr];
        for ( k=0; k<3; k++ ) r[k] = cb->x[k] - cb1->x[k];
        setupDerivs(order, r);
        transM2L(sys, dG0[0], dGk[0], cb1, cb);
#endif
      }
    }
  }

  /* downward pass */
  for ( lev=height; lev<depth; lev++ ) {
    for ( cb=sys->cubeList[lev]; cb != NULL; cb=cb->next ) {
      for ( nKid=0; nKid<cb->nKids; nKid++ ) {
        transL2L(sys, cb, cb->kids[nKid]);
      }
    }
  }

  /* L2P transformations */
  nMom = sys->nMom[sys->ordMom[depth]];
  for ( idx=0, cb=sys->cubeList[depth]; cb != NULL; cb=cb->next, idx++ ) {
    y = &(pot[cb->pnls->idx]);
    n = cb->nPnls;
    x = cb->lec_k1;
    dgemv_(&hChr, &nMom, &n, &one, L2P0[idx], &nMom, x, &inc, &zero, y, &inc);
    x = cb->lec_k2;
    dgemv_(&hChr, &nMom, &n, &one, L2P0[idx], &nMom, x, &inc, &one, y, &inc);

    y = &(pot[cb->pnls->idx+nPnls]);
    x = cb->lec_k3;
    dgemv_(&hChr, &nMom, &n, &one, L2P1[idx], &nMom, x, &inc, &zero, y, &inc);
    x = cb->lec_k4;
    dgemv_(&hChr, &nMom, &n, &one, L2P1[idx], &nMom, x, &inc, &one, y, &inc);
  }

  applyNearfield0(sys, sgm, pot);
  //for ( i=0; i<nPnls; i++ ) {
    //printf("%f %f\n",sgm[i],pot[i]);
    //pot[i] = 0.0;
    //pot[i+nPnls] = 0.0;
  //}
  //exit(0);

} /* applyFMM */
