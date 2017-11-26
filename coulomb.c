/*
 * coulomb.c: main driver
 * This program computes the screened Coulomb potential with fmm method
 *         exp(-K|x-y|)
 * \phi = --------------
 *            |x-y|
 * usage:
 *   coulomb [options] panelfile [options]
 *
 * Copyright Jiahui Chen
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gkGlobal.h"
#include "gk.h"

/* global variables */
int verbose=0, orderMom=0;
double kappa=0.12572533685666004, epsilon, epsilon1=1.0, epsilon2=40.0;
int doProtein=0, doEllipsoid=0, multiSpheres=0;
double Ra=50.0, l2Inv[3]; /* axes of sphere */

/* function pointers to kernel routines */
void (*kernel)(), (*kernelDC)(), (*kernelDS)();

/* routines used by the main routine */
panel *mkIco(int lev, int *nPanels, ssystem *sys);
void liftOntoEllipsoid(panel *snglist, double a, double b, double c, int nSurf);
panel *loadPanel(char *panelfile, char *density, int *numSing, ssystem *sys);
void makeImage(panel *pnlList, int *nPnls);
void gkInit(ssystem *sys, panel *pnlList, int order, int orderMom);
void dispCube(cube *cb);
void setupFMM(ssystem *sys);
void gmres(int n, ssystem *sys, int im, double *rhs, double *sol, int precond,
           int *maxits, double *eps );
void applyFMM(ssystem *sys, double *sgm, double *pot);
void directSum( ssystem *sys );
void setupDirect( ssystem *sys );
void applyDirect( ssystem *sys, double *sgm, double *pot );
void directSolver( ssystem *sys, double *sgm ) ;
void potential_molecule( ssystem *sys, double *r, double *sgm, double *potential );
double *panelRHS(int qOrder, panel *pnlX, double *chrY );

/* ub is solution on the boundary */
double ub(double *x, double *grad){
  double xBar=0.0, yBar=0.0, zBar=0.0;
  double rI, r2, r3I;
  double y;


} /* ub */

/* u1 is solution inside the boundary */
double u1(double *x, double *grad){
  /* double xBar=0.1, yBar=0.5, zBar=-0.1; */
  double xBar=0.0, yBar=0.0, zBar=0.0;
  double rI, r2, r3I;
  double y;

  r2 = SQR(x[0]-xBar) + SQR(x[1]-yBar) + SQR(x[2]-zBar);
  rI = 1.0/sqrt(r2);
  r3I = rI/r2;
  r3I *= fourPiI;
  grad[0] = -(x[0]-xBar)*r3I;
  grad[1] = -(x[1]-yBar)*r3I;
  grad[2] = -(x[2]-zBar)*r3I;

  y = 1.0/epsilon/(1.0+kappa*Ra)/Ra-1.0/Ra+rI;
  return y*fourPiI;
} /* u1 */

/* u2 is solution outside the boundary */
double u2(double *x, double *grad){
  double xBar=0.0, yBar=0.0, zBar=0.0;
  double r, rI, r2, expr, coef, fac;
  double y;

  r2 = SQR(x[0]-xBar) + SQR(x[1]-yBar) + SQR(x[2]-zBar);
  r = sqrt(r2);
  rI = 1.0/r;
  coef = fourPiI*exp( kappa*Ra )/epsilon/(1+kappa*Ra);
  expr = exp( -kappa*r );
  fac = (kappa + rI)*expr*coef/r2;
  grad[0] = -(x[0]-xBar)*fac;
  grad[1] = -(x[1]-yBar)*fac;
  grad[2] = -(x[2]-zBar)*fac;

  y = coef*expr*rI;
  return y;
} /* u2 */

/*
 *  setup right hand side (exterior Neumann problem)
 */
void setupRHS(ssystem *sys, double *sgm) {
  int i, j;
  int nPnls = sys->nPnls, nChar = sys->nChar, qOrder=sys->maxQuadOrder;
  double fac, x[3], ri, r2, r3i, ip, dudn;
  double *intgr;
  panel *pnl;

  fac = fourPiI/epsilon1;

  /* triangles order for Direct */
//  for ( i=0, pnl=sys->pnlOLst; pnl!=NULL; pnl=pnl->next, i++ ) {

  for ( i=0, pnl=sys->pnlLst; pnl!=NULL; pnl=pnl->nextC, i++ ) {
    sgm[i] = 0.0; sgm[nPnls+i] = 0.0;
    for ( j=0; j<nChar; j++ ) {
      /*
      x[0] = pnl->x[0] - sys->pos[3*j];
      x[1] = pnl->x[1] - sys->pos[3*j+1];
      x[2] = pnl->x[2] - sys->pos[3*j+2];
      r2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
      ri = 1.0/sqrt(r2);
      r3i = ri/r2;
      ip = pnl->normal[0]*x[0] + pnl->normal[1]*x[1] + pnl->normal[2]*x[2];
      dudn = -ip*r3i;
      sgm[i] += sys->chr[j]*ri;
      sgm[i+nPnls] += sys->chr[j]*dudn;
      */
      intgr=panelRHS(qOrder, pnl, &sys->pos[3*j]);
      sgm[i] += sys->chr[j]*intgr[0];
      sgm[i+nPnls] += sys->chr[j]*intgr[1];
    }
    //sgm[i] *= fac*pnl->area;
    //sgm[nPnls+i] *= fac*pnl->area;
    sgm[i] *= fac;
    sgm[nPnls+i] *= fac;
  }
  //exit(0);
} /* setupRHS */

/* sphereTest */
void sphereTest(ssystem *sys) {
  double dist;
  int i;
  panel *pnl;

  for ( pnl=sys->pnlLst; pnl!=NULL; pnl=pnl->nextC ) {
    dist = sqrt(SQR(pnl->x[0])+SQR(pnl->x[1])+SQR(pnl->x[2]));
    for ( i=0; i<3; i++ ) {
      pnl->x[i] = pnl->x[i]/dist*Ra;
      pnl->normal[i] = pnl->x[i]/Ra;
    }
  }

}/* sphereTest */

int main(int nargs, char *argv[]){
  char panelfile[80], density[80];
  int order=-1, image=0, refineLev=0, numSurfOne=1;
  int i, j, k, n, nPnls, nChar;
  int numItr, arnoldiSz;
  extern long memcount;
  panel *inputLst, *pnl;
  cube *cb;
  double tolpar=1.0e-6, tol;
  double *sgm, *pot, *pos, *chr;
  double ptl;
  ssystem *sys;

  double abserr, relerr, relinf_err, absinf_err;
  double abserr1, relerr1, relinf_err1, absinf_err1, temp1;
  double abserr2, relerr2, relinf_err2, absinf_err2, temp2;
  double realsol1, realsol2;
  double dist;

  CALLOC(sys, 1, ssystem, ON, AMISC);
  sys->height = 2;
  sys->depth = -1;
  sys->maxSepRatio = 0.8;
  sys->maxQuadOrder = 1;
  sys->nKerl = 4;
  sys->depth = 4;
  sprintf(density,"1");
  epsilon = epsilon2/epsilon1;

  /* parse the command line */
  panelfile[0] = 0;
  for ( i=1; i<nargs; i++ )
    if ( argv[i][0] == '-' )
      switch ( argv[i][1] ) {
        case 'E': doEllipsoid = 1;
          refineLev = atoi( argv[i]+2 );
          sprintf(panelfile, "E%d", refineLev);
          break;
        case 'S': sys->maxSepRatio = atof( argv[i]+3 );
          break;
        case 'e': tolpar = atof( argv[i]+3 );
          break;
        case 'h': sys->height = atoi( argv[i]+3 );
          break;
        case 'p':
          if ( argv[i][2] == '=' ) order = atoi( argv[i]+3 );
          if ( argv[i][2] == 'm' ) orderMom = atoi( argv[i]+4 );
          break;
        case 'q':
          sys->maxQuadOrder = atoi( argv[i]+3 );
          break;
        case 't': sys->depth = atoi( argv[i]+3 );
          break;
        case 'v': verbose = atoi( argv[i]+3 );
          break;
        case 'd': strcpy(density,argv[i]+3);
      }
    else
      strcpy(panelfile,argv[i]);

  if ( panelfile[0] == 0 ) {
    printf("\n Name of the panel file > ");
    scanf("%s",panelfile);
  }
  if ( sys->depth < 0 ) {
    printf("Select tree depth > ");
    scanf("%d", &sys->depth);
    if( sys->depth < 1 ) {
      printf("bad tree depth: %d\n", sys->depth );
      exit(0);
    }
  }

  /*
   * get panels by msms from pqr
   */
  if ( doEllipsoid ) {
    l2Inv[0] = 1.0/SQR(Ra);
    l2Inv[1] = 1.0/SQR(Ra);
    l2Inv[2] = 1.0/SQR(Ra);
    inputLst = mkIco(refineLev, &nPnls, sys);
    if ( verbose > 0 ) {
      printf("Icosahedral grid %d panels ", nPnls);
      printf("A=%3.1lf B=%3.1lf C=%3.1lf\n", Ra, Ra, Ra );
    }
    liftOntoEllipsoid(inputLst, Ra, Ra, Ra, 1);
  }
  else {
    inputLst = loadPanel(panelfile, density, &nPnls, sys);
  }
  if ( image ) makeImage(inputLst, &nPnls);

  printf("--- %s den=%s nPnls=%d nLev=%d ord=%d SepRat=%lg qOrd=%d\n",
         panelfile, density, nPnls, sys->depth, order, sys->maxSepRatio, sys->maxQuadOrder );

  sys->pnlOLst = inputLst;
  gkInit(sys, inputLst, order, orderMom);

//  dispCube(sys->cubeList[0]);
//  for (cb=sys->cubeList[2]; cb!=NULL; cb=cb->next){
//    dispCube(cb);
//  }
//  exit(1);

  CALLOC(sgm, 2*nPnls, double, ON, AMISC);
  CALLOC(pot, 2*nPnls, double, ON, AMISC);

  //setupFMM(sys);
  //sphereTest(sys);
  setupRHS(sys, sgm);
  setupFMM(sys);

  //setupDirect(sys);
  //exit(0);
  //directSum(sys);

  for ( i=0; i<2*nPnls; i++ ) pot[i] = sgm[i];

  /*
  //for ( i=0; i<nPnls; i++ ) {
    //sgm[i] = 1.0;
    //sgm[i+nPnls] = 1.0;
  //}

  applyFMM(sys, sgm, pot);
  for (  i=0, pnl=sys->pnlLst; pnl!=NULL; pnl=pnl->nextC, i++ ) {
    pot[i] = (1.0+epsilon)/2.0*pnl->area*sgm[i]-pot[i];
    pot[i+nPnls] = (1.0+1.0/epsilon)/2.0*pnl->area*sgm[i+nPnls]-pot[i+nPnls];
  }

  for ( i=0; i<nPnls; i++ ) {
    sgm[i] = 1.0;
    sgm[i+nPnls] = 1.0;
  }
  setupRHS(sys, sgm);

  double *pot1;
  CALLOC(pot1, 2*nPnls, double, ON, AMISC);
  applyDirect(sys, sgm, pot1);

  abserr=0.0; relerr=0.0; relinf_err=0.0; absinf_err=0.0;
  for ( i=0; i<nPnls; i++ ) {
    //printf("%.15f %.15f\n",sgm[i],sgm[i+nPnls]);
    temp1 = fabs(pot1[i]-pot[i]);
    relerr += temp1*temp1;
    abserr += pot[i]*pot[i];
    if ( temp1 > relinf_err ) relinf_err = temp1;
    temp1 = fabs(pot1[i]);
    if ( temp1 > absinf_err ) absinf_err = temp1;

    temp1 = fabs(pot1[i+nPnls]-pot[i+nPnls]);
    relerr += temp1*temp1;
    abserr += pot[i+nPnls]*pot[i+nPnls];
    if ( temp1 > relinf_err ) relinf_err = temp1;
    temp1 = fabs(pot1[i+nPnls]);
    if ( temp1 > absinf_err ) absinf_err = temp1;
  }
  relerr = sqrt(relerr/abserr);
  relinf_err =  relinf_err/absinf_err;
  printf("the Relative L2 and Inf norm: %.16f %.16f\n",relerr,relinf_err);
  exit(0);
  */

  //directSolver( sys, sgm );
  //for ( i=0; i<nPnls; i++ ) printf("%f\n",sgm[i]);


  n = 2*nPnls;
  tol = tolpar;
  arnoldiSz = 50;
  numItr = 2*arnoldiSz;
  gmres(n, sys, arnoldiSz, pot, sgm, 0, &numItr, &tol);
  printf("\ngmres-its=%d\n", numItr);


  double potential = 0.0, para=332.0716;
  for ( i=0; i<sys->nChar; i++ ) {
    potential_molecule( sys, &sys->pos[3*i], sgm, &ptl );
    potential += sys->chr[i]*twoPi*para*ptl;
  }
  //for ( i=0; i<80; i++ ) printf("%f\n",sgm[i]);
  printf("pot: %f\n",potential);

  abserr1=0.0; relerr1=0.0; relinf_err1=0.0; absinf_err1=0.0;
  abserr2=0.0; relerr2=0.0; relinf_err2=0.0; absinf_err2=0.0;
  sphereTest(sys);
  for ( i=0, pnl=sys->pnlLst; pnl!=NULL; pnl=pnl->nextC, i++ ) {
    //printf("%.15f %.15f\n",sgm[i],sgm[i+nPnls]);
    dist = sqrt(SQR(pnl->x[0])+SQR(pnl->x[1])+SQR(pnl->x[2]));
    //dist = 50.0;
    //printf("%.15f\n",dist);
    if ( dist - Ra > 1.0e-10 ) {
      realsol1 = fourPiI*exp(kappa*(Ra-dist))/epsilon/(1.0+kappa*Ra)
                /dist*sys->chr[0];
      realsol2 = -fourPiI/dist/dist/epsilon*sys->chr[0];
      //printf("%f %f\n",realsol1,realsol2);
    } else {
      realsol1 = fourPiI*(1.0/epsilon/(1.0+kappa*Ra)/Ra
               - 1.0/Ra + 1.0/dist)*sys->chr[0];
      realsol2 = -fourPiI/dist/dist*sys->chr[0];
    }
    //printf("%.16f %.16f %.16f\n",dist,realsol1,realsol2);
    /*//relative error
    temp1 = fabs(realsol1-sgm[i]);
    relerr1 += temp1*temp1;
    abserr1 += realsol1*realsol1;
    if ( temp1 > relinf_err1 ) relinf_err1 = temp1;
    temp1 = fabs(realsol1);
    if ( temp1 > absinf_err1 ) absinf_err1 = temp1;

    temp1 = fabs(realsol2-sgm[i+nPnls]);
    relerr2 += temp1*temp1;
    abserr2 += realsol2*realsol2;
    if ( temp1 > relinf_err2 ) relinf_err2 = temp1;
    temp1 = fabs(realsol2);
    if ( temp1 > absinf_err2 ) absinf_err2 = temp1;

    //printf("%.15f %.15f %.15f\n",dist,realsol1,realsol2);
    */
    temp1 = fabs(realsol1-sgm[i]);
    temp2 = temp1/fabs(realsol1);
    relerr1 += temp1*temp1; //e3
    //abserr1 += realsol1*realsol1;
    if ( temp1 > relinf_err1 ) relinf_err1 = temp1;
    if ( temp2 > absinf_err1 ) absinf_err1 = temp2;
    //printf("%.16f %.16f\n",relinf_err1,absinf_err1);

    temp1 = fabs(realsol2-sgm[i+nPnls]);
    temp2 = temp1/fabs(realsol2);
    relerr2 += temp1*temp1;
    if ( temp1 > relinf_err2 ) relinf_err2 = temp1;
    if ( temp2 > absinf_err2 ) absinf_err2 = temp1;

    sgm[i] = realsol1;
    sgm[i+nPnls] = realsol2;
  }
  printf("Error linf:%e %e l2:%e\n",relinf_err1,absinf_err1,sqrt(relerr1/nPnls));
  printf("Error: %e %e %e\n",relinf_err2,absinf_err2,sqrt(relerr2/nPnls));
  /*
  relerr1 = sqrt(relerr1/abserr1);
  relinf_err1 =  relinf_err1/absinf_err1;
  printf("the Relative L2 and Inf norm: %f %f\n",relerr1,relinf_err1);
  relerr2 = sqrt(relerr2/abserr2);
  relinf_err2 =  relinf_err2/absinf_err2;
  printf("the Relative L2 and Inf norm: %f %f\n",relerr2,relinf_err2);
  */

  potential = 0.0, para=332.0716;
  for ( i=0; i<sys->nChar; i++ ) {
    potential_molecule( sys, &sys->pos[3*i], sgm, &ptl );
    potential += sys->chr[i]*twoPi*para*ptl;
  }
  //for ( i=0; i<80; i++ ) printf("%f\n",sgm[i]);
  printf("pot: %f\n",potential);


}

void potential_molecule( ssystem *sys, double *r, double *sgm, double *potential ) {
  panel *pnlX;
  int i, nPnls=sys->nPnls;
  double dist, tp1, tp2, v[3], *nrmX, *ctr;
  double expKa, cos_theta, G0, Gk, G1, G2, H1, H2;

  *potential = 0.0;

  for ( i=0,pnlX=sys->pnlLst; i<nPnls; i++,pnlX=pnlX->nextC ) {
    ctr = pnlX->x;
    nrmX = pnlX->normal;
    v[0] = ctr[0]-r[0];
    v[1] = ctr[1]-r[1];
    v[2] = ctr[2]-r[2];
    dist = sqrt(SQR(v[0])+SQR(v[1])+SQR(v[2]));
//    printf("%f\n",dist);

    G0 = fourPiI/dist;
    expKa = exp(-kappa*dist);
    Gk = expKa*G0;

    cos_theta = (nrmX[0]*v[0]+nrmX[1]*v[1]+nrmX[2]*v[2])/dist;

    tp1 = G0/dist;
    tp2 = (1.0+kappa*dist)*expKa;

    G1 = cos_theta*tp1;
    G2 = tp2*G1;

    H1=G1-epsilon*G2;
    H2=G0-Gk;

    *potential+=pnlX->area*H1*sgm[i];
    *potential+=pnlX->area*H2*sgm[nPnls+i];

  }
  //exit(0);
}

/*
 * Matrix times Vector, subroutine of the iterative solver
 * the vector sgm and the result pot are ordered contiguously within cubes
*/
void MtV(ssystem *sys, double *sgm, double *pot) {
  int i, lev, inc=1, nPnls = sys->nPnls;
  cube *cb;
  panel *pnl;
  double output;

  applyFMM(sys, sgm, pot);
  //applyDirect(sys, sgm, pot);
  //output = 0.0;
//  for ( i=0;i<20;i++ ) {
//    output += pot[i];
//    printf("%f %f\n",sgm[i],pot[i]);
//  }
  for (  i=0, pnl=sys->pnlLst; pnl!=NULL; pnl=pnl->nextC, i++ ) {
    pot[i] = (1.0+epsilon)/2.0*pnl->area*sgm[i]-pot[i];
    pot[i+nPnls] = (1.0+1.0/epsilon)/2.0*pnl->area*sgm[i+nPnls]-pot[i+nPnls];
    //pot[i] = (1.0+epsilon)/2.0*sgm[i]-pot[i];
    //pot[i+nPnls] = (1.0+1.0/epsilon)/2.0*sgm[i+nPnls]-pot[i+nPnls];
  }
//  printf("\n%f",output);
  //exit(0);
} /* MtV */


/*
 * dummy routine, not implemented for FMM
*/
void PtV(ssystem *sys, double *sgm, double *pot) {
  int i, lev;
  cube *cb;

} /* PtV */
