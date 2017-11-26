/*
 *  routines to solve the linear system with cg()
 *
 *  Copyright: J. Tuasch
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gkGlobal.h"
#include "gk.h"

#define epsmac 1.e-15
#define DUMPSOL 0
extern int verbose;
extern double twoPi;

double *ss=NULL;
double **hh, *c, *s, *rs, *spc;
time_t time1, time2;

/* blas: return x^H*y */
double ddot_(int *n, double *x, int *incx, double *y, int *incy);
void daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);


void MtV(ssystem *sys, double *sgm, double *pot);
void PtV(ssystem *sys, double *sgm, double *pot);

void dumpSoln(ssystem *sys, double *sgm);

/*
 *  solve a linear system with conjugate gradients
 *  this is straight from Golub/VanLoan pg. 523
 *  stop the iteration when the initial residal is reduced by the factor eps
 *  returns in maxits the number of iterations and in eps the ratio
 *        ||final residual|| / ||initital residual||
 *
 *  jtdeb: hasn't been updated to work with hermitian matrices
 */
void cg(int nn, ssystem *sys, double *rhs, double *sol, int *maxits, double *eps) {
  double *res, *p, *y;
  double beta, alpha;
  double rnorm, rnormPrev, eps2, nrmy;
  int i, its, err_cde=0;

  time1 = time(&time1);

  CALLOC(res, nn, double, ON, ASOLVER);
  CALLOC(p, nn, double, ON, ASOLVER);
  CALLOC(y, nn, double, ON, ASOLVER);

  for ( rnorm=0.0, i=0; i<nn; i++ ) {
    res[i] = rhs[i];
    sol[i] = 0.0;
    rnorm += FABS2(res[i]);
  }

  eps2 = SQR(*eps)*rnorm;

  its=0;

  if ( verbose ) {
    printf(" **cg**  its=%d, ro=%lf\n", its, sqrt(rnorm) );
    fflush(stdout);
  }

  while ( rnorm > eps2 && its <= *maxits ) {
    its++;
    if ( its == 1 ) {
      for ( i=0; i<nn; i++ )
        p[i] = res[i];
    }
    else {
      beta = rnorm/rnormPrev;
      for ( i=0; i<nn; i++ )
        p[i] = res[i] + beta*p[i];
    }

    MtV(sys, p, y);

    for ( alpha=0.0, i=0; i<nn; i++ )
      alpha += p[i]*y[i];
    if ( (double) alpha < 0 ) {   /* jtdb: fix this for hermitian A */
      for ( nrmy=0.0, i=0; i<nn; i++ ) nrmy += FABS2(y[i]);
      nrmy = sqrt(nrmy);
      printf("\n**cg-error** Matrix not +def its=%d alpha=%lg\n",
             its, (double)alpha/nrmy);
      fflush(stdout);
    }
    alpha = rnorm/alpha;
    rnormPrev = rnorm;

    for ( rnorm=0.0,i=0; i<nn; i++ ) {
      sol[i] = sol[i] + alpha*p[i];
      res[i] = res[i] - alpha*y[i];
      rnorm += FABS2(res[i]);
    }
    if ( verbose ) {
#if DUMPSOL
      dumpSoln(sys, sol);
#endif
      printf(" **cg**  its=%d, ro=%lf, factor=%lf\n",
             its, sqrt(rnorm), sqrt(rnorm/rnormPrev) );
      fflush(stdout);
    }
  } /* while */
  *eps = sqrt(rnorm/eps2)*(*eps);
  *maxits = its;

  time2 = time(&time2);
  solveTimeNoPC = difftime(time2,time1);

} /* cg */

/*
 *  solve a linear system with preconditioned conjugate gradients
 *  this is straight from Golub/VanLoan pg. 529
 *  stop the iteration when the initial residal is reduced by the factor eps
 *  returns in maxits the number of iterations and in eps the ratio
 *        ||final residual|| / ||initital residual||
 *  jtdeb: hasn't been updated to work with hermitian matrices
 */
void pcg(int nn, ssystem *sys, double *rhs, double *sol, int *maxits, double *eps) {
  double *res, *p, *y, *z;
  double beta, alpha;
  double rnorm, rnormPrev, rTz, rTzPrev, eps2, nrmy;
  int i, its, err_cde=0;

  time1 = time(&time1);

  CALLOC(res, nn, double, ON, ASOLVER);
  CALLOC(p, nn, double, ON, ASOLVER);
  CALLOC(y, nn, double, ON, ASOLVER);
  CALLOC(z, nn, double, ON, ASOLVER);

  for ( rnorm=0.0, i=0; i<nn; i++ ) {
    res[i] = rhs[i];
    sol[i] = 0.0;
    rnorm += FABS2(res[i]);
  }

  eps2 = SQR(*eps)*rnorm;

  its=0;

  if ( verbose ) {
    printf(" **pcg**  its=%d, ro=%lf\n", its, sqrt(rnorm) );
    fflush(stdout);
  }

  while ( rnorm > eps2 && its <= *maxits ) {
    its++;
    PtV(sys, res, z);
    for ( rTz=0.0, i=0; i<nn; i++ )
      rTz += res[i]*z[i];

    if ( its == 1 ) {
      for ( i=0; i<nn; i++ )
        p[i] = z[i];
    }
    else {
      beta = rTz/rTzPrev;
      for ( i=0; i<nn; i++ )
        p[i] = z[i] + beta*p[i];
    }

    MtV(sys, p, y);

    for ( alpha=0.0, i=0; i<nn; i++ )
      alpha += p[i]*y[i];
    if ( (double)alpha < 0 ) {  /* jtdeb: fix for hermitian matrices */
      for ( nrmy=0.0, i=0; i<nn; i++ ) nrmy += SQR(y[i]);
      nrmy = sqrt(nrmy);
      printf("\n**pcg-error** Matrix not +def its=%d alpha=%lg\n",its, alpha/nrmy);
      fflush(stdout);
    }
    alpha = rTz/alpha;
    rnormPrev = rnorm;
    rTzPrev = rTz;

    for ( rnorm=0.0,i=0; i<nn; i++ ) {
      sol[i] = sol[i] + alpha*p[i];
      res[i] = res[i] - alpha*y[i];
      rnorm += FABS2(res[i]);
    }
    if ( verbose ) {
#if DUMPSOL
      dumpSoln(sys, sol);
#endif
      printf(" **pcg**  its=%d, ro=%lf, factor=%lf\n",
             its, sqrt(rnorm), sqrt(rnorm/rnormPrev) );
      fflush(stdout);
    }

  } /* while */
  *eps = sqrt(rnorm/eps2)*(*eps);
  *maxits = its;

  time2 = time(&time2);
  solveTimePC = difftime(time2,time1);

} /* pcg */



/*
 * Vector times Vector
*/
double VtV( int n_rows, double *x, double *y) {
  int i;
  double t;

  t = 0.0;
  for (i=0; i<n_rows; i++)
    t += x[i]*y[i];
  return t;
} /* VtV */


/*
 * This is a translation of gmrd.f, the original Saad version
 * Parameter list:
 * n	size of the problem
 * sys  system
 * im	size of Krylow subspace: should not exceed 50 in this version
 *	(can be reset in code. Looking at comment below)
 * rhs	right hand side
 * sol	initial guess on input, approximate solution on output
 * eps	tolerance for stopping criterion; process is stopped as soon as
 *	|| current residual || / || initial residual ||  <=  eps
 *	on output,  eps  =  || final residual || / || initial residual ||
 * precond flag for preconditioner
 * maxits maximum number of iterations allowed on input
 *	on output,  maxits = total number of iterations
 *
*/
void gmres(int n, ssystem *sys, int im, double *rhs, double *sol, int precond,
           int *maxits, double *eps ) {
  int its, i, j, k, ii, i1, k1, pos_i, pos_j, inc=1;
  double ro, ro1, gam, factor, avgfactor, eps1=1.0, tR;
  double t;

  time1 = time(&time1);

  if ( ss==NULL ) { /* first call: allocate space */
    CALLOC(ss, n*(im+1), double , ON, ASOLVER);
    CALLOC(hh, (im+2), double* , ON, ASOLVER);
    for ( i=0; i<=im+1; i++ ) {
      CALLOC( hh[i], (im+1), double, ON, ASOLVER);
    }
    CALLOC(c, (im+1), double, ON, ASOLVER);
    CALLOC(s, im, double, ON, ASOLVER);
    CALLOC(rs, (im+1), double, ON, ASOLVER);
  }
  if ( precond && spc==NULL ) CALLOC(spc, n, double , ON, ASOLVER);
  its = 0;
  avgfactor = 0.0;

  if ( precond ) {
    PtV(sys, rhs, spc);
    for ( i=0; i<n; i++ ) {
      rhs[i] = spc[i];
    }
  }

  ro = sqrt( ddot_(&n, sol, &inc, sol, &inc) );
  if ( ro <= epsmac ) goto ret;
  eps1 = *eps*ro;

  /*
   * outer loop starts here
   */
           l_10:
  /*    compute initial residual vectors */
  if ( precond ) {
    MtV(sys, sol, spc);
    PtV(sys, spc, ss);
  }
  else {
    MtV(sys, sol, ss);
  }

  for ( j=0; j<n; j++ ) {
    ss[j] =  rhs[j] - ss[j];
  } /* 21 */
  ro1 = ro;

  ro = sqrt( ddot_(&n, ss, &inc, ss, &inc) );
  if ( ro <= epsmac ) goto ret;

  for ( j=0; j<n; j++ ) {
    ss[j] =  ss[j]/ro;
  } /* 210 */
  if ( verbose ) {
    printf("\n ** gmres, restart:\t its=%d ro=%lf",its,ro);
    fflush(stdout);
  }
  /* initialize 1-st term of Hessenberg system */
  rs[1] = ro;
  i = 0;
           l_4:
  pos_i = n*i;
  i++;
  its++;
  i1 = i + 1;
  if ( precond ) {
    MtV(sys, ss+pos_i, spc);
    PtV(sys, spc, ss+pos_i+n);
  }
  else {
    MtV(sys, ss+pos_i, ss+pos_i+n);
  }
  /*
   * modified Gram - Schmidt
   */
  pos_i = i*n;
  pos_j = 0;
  for ( j=1; j<=i; j++ ) {
    t = ddot_(&n, ss+pos_j, &inc, ss+pos_i, &inc);
    hh[j][i] = t;
    t = -t;
    daxpy_(&n, &t, ss+pos_j, &inc, ss+pos_i, &inc);
    pos_j += n;
  } /* 55 */

  tR = sqrt( ddot_(&n, ss+pos_i, &inc, ss+pos_i, &inc) );
  hh[i1][i] = tR;
  for ( k=0; k<n; k++ ) {
    ss[pos_i+k] /= tR; /* 57 */
  }

  /* done with modified Schmidt and Arnoldi step
   * now update factorization of hh
  */
  if ( i != 1 )
    /* perform previous transformations on i-th column of h */
    for ( k=2; k<=i; k++ ) {
      k1 = k - 1;
      t = hh[k1][i];
      hh[k1][i] = c[k1]*t + s[k1]*hh[k][i];
      hh[k][i] = -s[k1]*t + c[k1]*hh[k][i];
    } /* 66 */
  gam = sqrt( FABS2(hh[i][i]) + FABS2(hh[i1][i]) );
  if ( gam == 0.0 ) gam = epsmac;

  /* determine next plane rotation */
  c[i] = hh[i][i]/gam;
  s[i] = hh[i1][i]/gam;
  rs[i1] = -s[i]*rs[i];
  rs[i] = c[i]*rs[i];

  /* determine residual norm and test for convergence */
  hh[i][i] = c[i]*hh[i][i] + s[i]*hh[i1][i];
  ro1 = ro;
  ro = FABS(rs[i1]);
  factor =  ro/ro1;
  if ( verbose ) {
    printf("\n ** gmres: \tits=%d ro=%lf factor=%lf",its,ro,factor);
    fflush(stdout);
  }
  avgfactor += factor;

  if ( (i<im) && (ro>eps1) && (its<*maxits) ) goto l_4;

  /*
   *  now compute solution. First solve upper triangular system
  */
  rs[i] = rs[i]/hh[i][i];
  for ( ii=2; ii<=i; ii++ ) {
    k = i - ii + 1;
    k1 = k + 1;
    t = rs[k];
    for (j=k1;j<=i;j++) {
      t = t - hh[k][j]*rs[j];
    }  /* 40 */
    rs[k] = t/hh[k][k];
  } /* 30 */

  /*
   * Form a linear combination to get solution
  */
  pos_j = 0;
  for ( j=1; j<=i; j++) {
    t = rs[j];
    daxpy_(&n, &t, ss+pos_j, &inc, sol, &inc);
    pos_j += n;
  } /* 16 */

  /* restart outer loop when necessary */
  if ( (ro > eps1) && (its<*maxits) ) goto l_10;

           ret:
  *eps = (ro/eps1)*(*eps);
  *maxits = its;
  if  ( verbose && (its>0) ) {
    avgfactor = avgfactor/(double)its;
    printf("\n **gmres** \tits = %d\t eps = %lf\t af = %lf",its,*eps,avgfactor);

  }

  time2 = time(&time2);
  if ( precond ) {
    solveTimePC = difftime(time2,time1);
  }
  else {
    solveTimeNoPC = difftime(time2,time1);
  }

} /* gmres */
