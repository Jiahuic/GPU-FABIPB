/*
 * gkDisplay.c
 *   routines that print stuff for debugging
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "gkGlobal.h"
#include "gk.h"

#define MEG 1048576

extern double minx, miny, minz;         /* corner of level-0 cube */
extern double *sideLen;                 /* sidelengths of uniform cubes */

/*
 * print synopsis of memory usage, non-zero matrix elmts, etc
 */
void dumpStats(ssystem *sys) {
  double ratio;
  long total;
  int lev;
  cube *cb;


  printf("\nLevel:             ");
  for ( lev=sys->depth; lev>=sys->height; lev-- ) printf(" %3d",lev);
  printf("\nOrder of Moments:  ");
  for ( lev=sys->depth; lev>=sys->height; lev-- ) printf(" %3d",sys->ordMom[lev]);
  printf("\nOrder of M2L's:    ");
  for ( lev=sys->depth; lev>=sys->height; lev-- ) printf(" %3d",sys->ordM2L[lev]);

  printf("\nNumber panels=%d, vertices=%d",
         sys->nPnls, sys->nVtxs);
  printf("\nMemory: Total                 %lg MB\n", ((double)memcount)/MEG);
  printf("        Panels/Vertices/Edges %lg MB\n", ((double)memPVE)/MEG);
  printf("        Cubes                 %lg MB\n", ((double)memCUBES)/MEG);
  printf("        Q2P transforms        %lg MB\n", ((double)memQ2P)/MEG);
  printf("        Q2M/L2P transforms    %lg MB\n", ((double)memQ2M)/MEG);
  printf("        M2L transforms        %lg MB\n", ((double)memM2L)/MEG);
  printf("        Solver                %lg MB\n", ((double)memSOLVER)/MEG);
  printf("        Misc                  %lg MB\n", ((double)memMISC)/MEG);
  printf("Time for solving w/o PC            %6.0lf sec\n", solveTimeNoPC);
  printf("Time for setup (Q2P+Q2M+M2L=total) %6.0lf %6.0lf %6.0lf %6.0lf sec\n",
         setupQ2PTime, setupQ2MTime, setupM2LTime,
         setupQ2PTime+setupQ2MTime+setupM2LTime );
#if NUMKERNEL
  printf("Number of kernel evaluations:  Quad=%lg / M2L=%lg\n",
         (double)numKernRealEval, (double)numKernCplxEval);
#endif

}

/*
 * dump panels to stdout (in input order)
 */
void dumpPanels(panel *pnlLst) {
  panel *pnl;
  int i;


  for ( pnl=pnlLst; pnl!=NULL; pnl=pnl->next  ) {

    if ( pnl->shape==4 ) { /* quadrilateral */
      printf("sorry, no quadrilaterals\n");
    }
    else {
      if ( pnl->shape==3 ) { /* triangle */
#if 0
        printf("T %d %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", pnl->nSurf,
               pnl->vtx[0]->x[0], pnl->vtx[0]->x[1], pnl->vtx[0]->x[2],
               pnl->vtx[1]->x[0], pnl->vtx[1]->x[1], pnl->vtx[1]->x[2],
               pnl->vtx[2]->x[0], pnl->vtx[2]->x[1], pnl->vtx[2]->x[2] );
#else
        printf("T %d %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", pnl->nSurf,
               pnl->vtx[0][0], pnl->vtx[0][1], pnl->vtx[0][2],
               pnl->vtx[1][0], pnl->vtx[1][1], pnl->vtx[1][2],
               pnl->vtx[2][0], pnl->vtx[2][1], pnl->vtx[2][2] );
#endif
      }
      else {
        fprintf(stderr,"\nonly triangles or quadrilaterals, you dodo!\n");
        exit(1);
      }
    }
  }
} /* dumpPanels */



/*
 * print a matrix
 */
void dispMax( int nRows, int nCols, double *A ) {
  int i,j;

  for (i=0; i<nRows; i++ ) {
    for (j=0; j<nCols; j++ )
      printf("%lf ", A[j*nRows + i]);
    printf("\n");
  }

} /* dispMax */




/*
 * print info about one cube
 */
void dispCube(cube *cb) {
  int i, lev;
  cube *nbr;
  double x, y, z;
  double dist, maxDist;

  lev = cb->level;
  printf("\n-------\t lev=%d (i,j,k)=(%d,%d,%d)-------\n", lev, cb->i, cb->j, cb->k );
  x =  sideLen[lev]*((double) cb->i + 0.5) + minx;
  y =  sideLen[lev]*((double) cb->j + 0.5) + miny;
  z =  sideLen[lev]*((double) cb->k + 0.5) + minz;
  printf("  ChebCtr=(%lf,%lf,%lf) ChebRad=%lf\n",
         cb->x[0], cb->x[1], cb->x[2], cb->eRad);

  printf("  Ctr=(%lf,%lf,%lf) SideLen=%lf nPnls=%d \n",
         x, y, z, sideLen[lev], cb->nPnls);
  /* Kids */
#if 1
  printf("  %d kids:", cb->nKids);
  for (i=0; i<cb->nKids; i++) {
    printf("(%d,%d,%d),", cb->kids[i]->i,cb->kids[i]->j,cb->kids[i]->k );
  }
#endif

#if 1
  /* First neighbors */
  printf("\n  %d 1st nbrs:", cb->nNbrs);
  maxDist = 0;
  for (i=0; i<cb->nNbrs; i++) {
    nbr = cb->nbrs[i];
    dist = sqrt(SQR(cb->x[0]-nbr->x[0]) + SQR(cb->x[1]-nbr->x[1]) + SQR(cb->x[2]-nbr->x[2]));
    if ( dist > maxDist ) maxDist = dist;
    printf("(%d,%d,%d),", nbr->i, nbr->j,nbr->k );
  }
  printf(" max-distance=%lf", maxDist );
  /* Interaction list */
  printf("\n  %d 2nd nbrs:", cb->n2Nbrs);
  maxDist = 0;
  for ( i=cb->n2Nbrs-1; i>=cb->nNbrs; i-- ) {
    nbr = cb->nbrs[i];
    dist = sqrt(SQR(cb->x[0]-nbr->x[0]) + SQR(cb->x[1]-nbr->x[1]) + SQR(cb->x[2]-nbr->x[2]));
    if ( dist > maxDist ) maxDist = dist;
    printf("(%d,%d,%d),", nbr->i, nbr->j,nbr->k );
  }
  printf(" max-distance=%lf\n", maxDist );
#endif
} /* dispCube */


#if 0
/*
 * print a list of all moments
 */
void dispMoments(ssystem *sys) {
  int lev, i, numMoments;
  cube *cb;

  for ( lev=sys->depth; lev>=sys->height; lev-- ) {
    printf("\n  ======================= \n");
    printf("\t lev=%d", lev);
    printf("\n  ======================= \n");
    numMoments = sys->nMom[sys->ordMom[lev]];
    for ( cb=sys->cubeList[lev]; cb!=NULL; cb=cb->next ) {
      if ( cb->MomentsL != NULL ) {
        printf("\n------- lev=%d ctr=(%lf,%lf,%lf) nr=(%d,%d,%d) left -------\n",
               cb->level, cb->x[0], cb->x[1], cb->x[2], cb->i, cb->j, cb->k);
        dispMax(numMoments,cb->nPhiL,cb->MomentsL);
      }
      if ( cb->MomentsR != NULL ) {
        if ( sys->layer==DOUBLE || sys->layer==ADJOINT ) {
          printf("\n------- lev=%d ctr=(%lf,%lf,%lf) nr=(%d,%d,%d) right -------\n",
                 cb->level, cb->x[0], cb->x[1], cb->x[2], cb->i, cb->j, cb->k);
          dispMax(numMoments,cb->nPhiR,cb->MomentsR);
        }
      }
    }
  }
} /* dispMoments */
#endif

/* dump Panel information */
void dispPanel(panel *pnl) {
  int i;
  double c[4][3], *vtx;

  printf("shape=%d area=%g", pnl->shape, pnl->area);

  printf(" centr=[%g %g %g]\n", pnl->x[0], pnl->x[1], pnl->x[2]);
  printf("normal=[%g %g %g]\n", pnl->normal[0], pnl->normal[1], pnl->normal[2]);
  printf("vertices (or edges)\n");
  for (i=0; i < pnl->shape; i++) {
    printf("\t (%lf, %lf, %lf)\n", pnl->vtx[0][0], pnl->vtx[0][1], pnl->vtx[0][2] );
    printf("\t (%lf, %lf, %lf)\n", pnl->vtx[1][0], pnl->vtx[1][1], pnl->vtx[1][2] );
    printf("\t (%lf, %lf, %lf)\n", pnl->vtx[2][0], pnl->vtx[2][1], pnl->vtx[2][2] );
  }

  printf("\n");
} /* dispPanel */

/*
 * print all cubes below the cube cb
 */
void dispCubeTree(cube *cb){
  int i,j;
  for ( j=0; j<cb->level; j++ ) printf(" ");
  printf("lev=%d [%d %d %d] nPnls=%d nKids=%d\n",
         cb->level, cb->i,cb->j,cb->k,cb->nPnls,cb->nKids);
  for ( i=0; i<cb->nKids; i++ ){
    dispCubeTree(cb->kids[i]);
  }
}
