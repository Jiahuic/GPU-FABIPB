/*
 * input.c
 *   various input/output routines for gk package
 *
 *   copyright by Johannes Tausch
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gkGlobal.h"
#include "gk.h"
#define MAXCOND 50
#define DENSITY 10    /* for mkMultiSpheres() */

int nSurf = 0;

/*
 * calculate area and normal
 */
double triangle_area(double v[3][3]){
  int i;
  double a[3], b[3], c[3], aa, bb, cc, ss, t_area;
  for (i=0;i<3;i++){
    a[i] = v[0][i]-v[1][i];
    b[i] = v[0][i]-v[2][i];
    c[i] = v[1][i]-v[2][i];
  }
  aa = sqrt(SQR(a[0])+SQR(a[1])+SQR(a[2]));
  bb = sqrt(SQR(b[0])+SQR(b[1])+SQR(b[2]));
  cc = sqrt(SQR(c[0])+SQR(c[1])+SQR(c[2]));
  ss = 0.5*(aa+bb+cc);
  t_area = sqrt(ss*(ss-aa)*(ss-bb)*(ss-cc));
  return(t_area);
}

/*
 * loadpanel returns a list of panel structs derived from passed data:
 * shape, vertices, and type.
 */
panel *loadPanel(char *panelfile, char *density, int *numSing, ssystem *sys) {
  int i, j, k, ii, shape, type, nSurf;
  panel *pnlList, *pnl;
  char fpath[256], fname[256];
  FILE *fp, *wfp;

  char c,c1[10],c2[10],c3[10],c4[10],c5[10];
  double a1,a2,a3,b1,b2,b3;//a_norm,r0_norm,v0_norm;
  int i1,i2,i3,j1,j2,j3,ierr,iface[3],jface[3],ialert;
  double den,prob_rds,xx[3],yy[3],face[3][3],tface[3][3],s_area;
  double **sptpos, **sptnrm;
  int **extr_v, **extr_f, *nvert;
  double dist_local, area_local, cpuf;
  int nspt, natm, nface;

  /* read in vertices */
  sys->nChar = 0;
  sprintf(fpath,"test_proteins/");
  sprintf(fname,"%s%s.pqr",fpath,panelfile);
  fp=fopen(fname,"r");
  sprintf(fname,"%s%s.xyzr",fpath,panelfile);
  wfp=fopen(fname,"w");
  while(fscanf(fp,"%s %s %s %s %s %lf %lf %lf %lf %lf",c1,c2,c3,
               c4,c5,&a1,&a2,&a3,&b1,&b2) != EOF){
    fprintf(wfp,"%f %f %f %f\n",a1,a2,a3,b2);
    sys->nChar++;
  }
  fclose(fp);
  fclose(wfp);

  CALLOC(sys->pos, 3*sys->nChar, double, ON, APVE);
  CALLOC(sys->chr, sys->nChar, double, ON, APVE);
  sprintf(fname,"%s%s.pqr",fpath,panelfile);
  fp=fopen(fname,"r");
  for (i=0;i<sys->nChar;i++){
    ierr=fscanf(fp,"%s %s %s %s %s %lf %lf %lf %lf %lf",c1,c2,c3,
                 c4,c5,&a1,&a2,&a3,&b1,&b2);
    sys->pos[3*i]=a1;
    sys->pos[3*i+1]=a2;
    sys->pos[3*i+2]=a3;
    sys->chr[i]=b1;
  }
  fclose(fp);

  /* run msms */
  sprintf(fname,"msms -if %s%s.xyzr -prob 1.4 -de %s -of %s%s > msms.output",
                fpath, panelfile, density, fpath, panelfile);
  printf("%s\n",fname);
  ierr=system(fname);
  sprintf(fname,"rm msms.output");
  ierr=system(fname);
  /*======================================================================*/

  /* read in vert */
  sprintf(fname, "%s%s.vert",fpath,panelfile);
  printf("%s\n",fname);

  /* open the file and read through the first two rows */
  fp=fopen(fname,"r");
  for (i=1;i<=2;i++){
    while (c=getc(fp)!='\n'){
   }
  }

  ierr=fscanf(fp,"%d %d %lf %lf ",&nspt,&natm,&den,&prob_rds);
  //printf("nspt=%d, natm=%d, den=%lf, prob=%lf\n", nspt,natm,den,prob_rds);

  /*allocate variables for vertices file*/
  CALLOC(sptpos, 3, double*, ON, APVE);
  CALLOC(sptnrm, 3, double*, ON, APVE);
  for (i=0;i<3;i++) {
    CALLOC(sptpos[i], nspt, double, ON, APVE);
    CALLOC(sptnrm[i], nspt, double, ON, APVE);
  }

  for (i=0;i<=nspt-1;i++){
    ierr=fscanf(fp,"%lf %lf %lf %lf %lf %lf %d %d %d",&a1,&a2,&a3,&b1,&b2,&b3,&i1,&i2,&i3);

    sptpos[0][i]=a1;
    sptpos[1][i]=a2;
    sptpos[2][i]=a3;
    sptnrm[0][i]=b1;
    sptnrm[1][i]=b2;
    sptnrm[2][i]=b3;

  }
  fclose(fp);
  printf("finish reading vertices file\n");

  /* read in faces */
  ierr=sprintf(fname, "%s%s.face",fpath,panelfile);
  fp=fopen(fname,"r");
  for (i=1;i<=2;i++){ while (c=getc(fp)!='\n'){} }

  ierr=fscanf(fp,"%d %d %lf %lf ",&nface,&natm,&den,&prob_rds);
  //printf("nface=%d, natm=%d, den=%lf, prob=%lf\n", nface,natm,den,prob_rds);

  /* allocate variables for vertices file */
  CALLOC(nvert, 3*nface, int, ON, APVE);

  for (i=0;i<=nface-1;i++){
    ierr=fscanf(fp,"%d %d %d %d %d",&j1,&j2,&j3,&i1,&i2);
    nvert[3*i]=j1;
    nvert[3*i+1]=j2;
    nvert[3*i+2]=j3;
  }
  fclose(fp);
  printf("finish reading face file\n");

  /* we delete ill performence triangles */
  s_area = 0.0;
  printf("Number of surfaces = %d\n",nface);
  printf("Number of surface points = %d\n",nspt);
  time_t cpu1 = time(NULL);

  pnlList = NULL;
  *numSing = 0;
  for (i=0;i<nface;i++){

    for (j=0;j<3;j++){
  iface[j]=nvert[3*i+j];
  xx[j]=0.0;
    }
    for (j=0; j<3; j++){
  for (k=0; k<3; k++){
    face[k][j]=sptpos[j][iface[k]-1];
	  xx[j] += 1.0/3.0*(face[k][j]);
  }
    }
    /* compute the area of each triangule */
    area_local=triangle_area(face);
    /* if the point is too close with 10 points infront, delete it */
    ialert=0;
    for (j=i-10; (j>=0&&j<i);j++){ /* like k=max(1,i-10), i-1 in fortran */
  for (k=0; k<3; k++) {
    jface[k] = nvert[3*j+k];
    yy[k]=0.0;
  }
  for (k=0;k<3;k++){
    for (ii=0;ii<3;ii++){
      tface[ii][k] = sptpos[k][jface[ii]-1];
      yy[k] += 1.0/3.0*(tface[ii][k]);
    }
  }
  dist_local = 0.0;
  for(k=0;k<3;k++) dist_local +=(xx[k]-yy[k])*(xx[k]-yy[k]);/* dot_product */
  dist_local=sqrt(dist_local);
  if (dist_local<1.0e-6) {
    ialert=1;
    printf("particles %d and %d are too close: %e\n", i,j,dist_local);
  }
    }

    /* allocate the panel */
    if(area_local>=1.0e-5 && ialert == 0){
  (*numSing)++;
  if (pnlList == NULL) {
    CALLOC(pnlList, 1, panel, ON, APVE);
    pnl = pnlList;
  }
  else{
    CALLOC(pnl->next, 1, panel, ON, APVE);
    pnl = pnl->next;
  }
  /* Fill in corners. */
  for (j=0;j<3;j++){
    for (k=0;k<3;k++){
      pnl->vtx[j][k] = sptpos[k][iface[j]-1];
    }
  }
  pnl->shape = 3;
  pnl->nSurf = *numSing;
  s_area += area_local;
    }
  }

  printf("total area = %f\n",s_area);

  time_t cpu2 = time(NULL);
  cpuf = ((double)(cpu2-cpu1));
  printf("total MSMS post-processing time = %f\n",cpuf);

  sprintf(fname,"rm %s%s.xyzr",fpath,panelfile);
  ierr=system(fname);
  sprintf(fname,"rm %s%s.vert",fpath,panelfile);
  ierr=system(fname);
  sprintf(fname,"rm %s%s.face",fpath,panelfile);
  ierr=system(fname);

  for (i=0;i<3;i++){
    free(sptpos[i]);
    free(sptnrm[i]);
  }
  free(sptpos);
  free(sptnrm);
  free(nvert);

  return pnlList;
} /* loadPanel */

/*
 * reads a panellist, and appends a list of their mirror image
 * across the plane x_3 = 0. This is useful when we have
 * zero Dirichlet boundary condition on this plane.
 * An alternative is to use the method of images, but this would be difficult
 * in this implementation.
 * Caution: we do not check whether the original domain is above the plane.
 */
void makeImage(panel *pnlList, int *nSing) {
  panel *pnlList2 = NULL;
  panel *pnl, *pnlNew;
  int j, nSurf=0;

  /* nSurf is the number of surfaces in the original domain */
  pnl = pnlList;
  while ( pnl != NULL ) {
    if ( pnl->nSurf > nSurf ) nSurf = pnl->nSurf;
    pnl = pnl->next;
  }

  pnl=pnlList;
  while ( 1 ) {
    if(pnlList2 == NULL) {
      CALLOC(pnlList2, 1, panel, ON, APVE);
      pnlNew = pnlList2;
    }
    else {
      CALLOC(pnlNew->next, 1, panel, ON, APVE);
      pnlNew = pnlNew->next;
    }
    /* Fill in the corners: just flip the z-component
     * note that the orientation must be changed to ensure the normal
     * of the flipped panel points out of the flipped domain
     * jtdeb: only implemented for triangles.
     */
    pnlNew->vtx[0][0] = pnl->vtx[0][0];
    pnlNew->vtx[0][1] = pnl->vtx[0][1];
    pnlNew->vtx[0][2] = -pnl->vtx[0][2];

    pnlNew->vtx[1][0] = pnl->vtx[2][0];
    pnlNew->vtx[1][1] = pnl->vtx[2][1];
    pnlNew->vtx[1][2] = -pnl->vtx[2][2];

    pnlNew->vtx[2][0] = pnl->vtx[1][0];
    pnlNew->vtx[2][1] = pnl->vtx[1][1];
    pnlNew->vtx[2][2] = -pnl->vtx[1][2];

    /* Fill in descriptors. */
    pnlNew->shape = pnl->shape;
    pnlNew->nSurf = pnl->nSurf+nSurf;
    if ( pnl->next == NULL ) {
      pnl->next = pnlList2; /* append flipped panels to the old list */
      pnlNew->next = NULL;
      *nSing *= 2;
      return;
    }
    else {
      pnl = pnl->next;
    }
  }

} /* makeImage */

/*
 * Uniform refinement: Each panel is refined into four panels of equal size.
 * The refined panels are inserted into the panel list.
 * Only triangles are supported
 */
void unifRefine(panel *pnlList, int *nPanels, double raduis) {
   panel *pnl, *pnlNext;
   int j, k, nSurf;
   double v0[3], v1[3], v2[3], w0[3], w1[3], w2[3];
   double leng1, leng2, leng3;

   pnl=pnlList;

   while ( pnl != NULL ) {

     pnlNext = pnl->next;
     if ( pnl->shape != 3 ) {
       printf("\n**error** unifRefine(): only triangles\n");
       exit(1);
     }

     leng1 = leng2 = leng3 = 0.0;
     for ( k=0; k<3; k++ ) {
       v0[k] = pnl->vtx[0][k];
       v1[k] = pnl->vtx[1][k];
       v2[k] = pnl->vtx[2][k];
       w0[k] = 0.5*(v1[k] + v2[k]);
       w1[k] = 0.5*(v2[k] + v0[k]);
       w2[k] = 0.5*(v0[k] + v1[k]);
       leng1 += SQR(w0[k]);
       leng2 += SQR(w0[k]);
       leng3 += SQR(w0[k]);
     }
     leng1 = sqrt(leng1);
     leng2 = sqrt(leng2);
     leng3 = sqrt(leng3);

     for ( k=0; k<3; k++ ) {
       w0[k] *= raduis/leng1;
       w1[k] *= raduis/leng2;
       w2[k] *= raduis/leng3;
       pnl->vtx[0][k] = v0[k];
       pnl->vtx[1][k] = w2[k];
       pnl->vtx[2][k] = w1[k];
     }

     CALLOC(pnl->next, 1, panel, ON, APVE);
     pnl = pnl->next;
     for ( k=0; k<3; k++ ) {
       pnl->vtx[0][k] = w2[k];
       pnl->vtx[1][k] = v1[k];
       pnl->vtx[2][k] = w0[k];
     }
     pnl->nSurf = nSurf++;
     //printf("%d\n",pnl->nSurf);
     pnl->shape = 3;

     CALLOC(pnl->next, 1, panel, ON, APVE);
     pnl = pnl->next;
     for ( k=0; k<3; k++ ) {
       pnl->vtx[0][k] = w0[k];
       pnl->vtx[1][k] = v2[k];
       pnl->vtx[2][k] = w1[k];
     }
     pnl->nSurf = nSurf++;
     //printf("%d\n",pnl->nSurf);
     pnl->shape = 3;

     CALLOC(pnl->next, 1, panel, ON, APVE);
     pnl = pnl->next;
     for ( k=0; k<3; k++ ) {
       pnl->vtx[0][k] = w0[k];
       pnl->vtx[1][k] = w1[k];
       pnl->vtx[2][k] = w2[k];
     }
     pnl->nSurf = nSurf++;
     //printf("%d\n",pnl->nSurf);
     pnl->shape = 3;

     pnl->next = pnlNext;
     pnl = pnlNext;
   }

   *nPanels *= 4;

}/* unifRefine */

/*
 * Make the icosahedron.
 * The panels are oriented such that the normal determined by the right hand
 * rule points into the exterior of the cube.
 *
 * Parameters:
 *  lev       refinement level
 *  nPanels   (output) number of panels:
 *                     lev=0=>nPanels=48, lev=1=>nPanels=192,...
 */
panel *mkIco(int lev, int *nPanels, ssystem *sys) {
  panel *pnlList=NULL, *pnl;

  int cnt, j, k, level;
  /* 20 panels of initial triangulation */
  double vtx[20][9] = {
    {0.85065080835203993218, 0, 0.525731112119133606025,
      0.85065080835203993218, 0, -0.525731112119133606025,
      0.525731112119133606025, 0.85065080835203993218},
    {0.85065080835203993218, 0, 0.525731112119133606025,
      0.525731112119133606025, -0.85065080835203993218, 0,
      0.85065080835203993218, 0, -0.525731112119133606025},
    {-0.85065080835203993218, 0, -0.525731112119133606025,
      -0.85065080835203993218, 0, 0.525731112119133606025,
      -0.525731112119133606025, 0.85065080835203993218},
    {-0.85065080835203993218, 0, 0.525731112119133606025,
      -0.85065080835203993218, 0, -0.525731112119133606025,
      -0.525731112119133606025, -0.85065080835203993218},
    {0.525731112119133606025, 0.85065080835203993218, 0,
      -0.525731112119133606025, 0.85065080835203993218, 0,
      0, 0.525731112119133606025, 0.85065080835203993218},
    {-0.525731112119133606025, 0.85065080835203993218, 0,
      0.525731112119133606025, 0.85065080835203993218, 0,
      0, 0.525731112119133606025, -0.85065080835203993218 },
    {0, -0.525731112119133606025, -0.85065080835203993218,
      0, 0.525731112119133606025, -0.85065080835203993218,
      0.85065080835203993218, 0, -0.525731112119133606025 },
    {0, 0.525731112119133606025, -0.85065080835203993218,
      0, -0.525731112119133606025, -0.85065080835203993218,
      -0.85065080835203993218, 0, -0.525731112119133606025},
    {0.525731112119133606025, -0.85065080835203993218, 0,
      -0.525731112119133606025, -0.85065080835203993218, 0,
      0, -0.525731112119133606025, -0.85065080835203993218},
    {-0.525731112119133606025, -0.85065080835203993218, 0,
      0.525731112119133606025, -0.85065080835203993218, 0,
      0, -0.525731112119133606025, 0.85065080835203993218},
    {0, 0.525731112119133606025, 0.85065080835203993218, 0,
      -0.525731112119133606025, 0.85065080835203993218,
      0.85065080835203993218, 0, 0.525731112119133606025 },
    {0, -0.525731112119133606025, 0.85065080835203993218,
      0, 0.525731112119133606025, 0.85065080835203993218,
      -0.85065080835203993218, 0, 0.525731112119133606025},
    {0.525731112119133606025, 0.85065080835203993218, 0,
      0.85065080835203993218, 0, -0.525731112119133606025,
      0, 0.525731112119133606025, -0.85065080835203993218},
    {0.85065080835203993218, 0, 0.525731112119133606025,
      0.525731112119133606025, 0.85065080835203993218, 0,
      0, 0.525731112119133606025, 0.85065080835203993218 },
    {-0.85065080835203993218, 0, -0.525731112119133606025,
      -0.525731112119133606025, 0.85065080835203993218, 0,
      0, 0.525731112119133606025, -0.85065080835203993218 },
    {-0.525731112119133606025, 0.85065080835203993218, 0,
      -0.85065080835203993218, 0, 0.525731112119133606025,
      0, 0.525731112119133606025, 0.85065080835203993218 },
    {0.85065080835203993218, 0, -0.525731112119133606025,
      0.525731112119133606025, -0.85065080835203993218, 0,
      0, -0.525731112119133606025, -0.85065080835203993218},
    {0.525731112119133606025, -0.85065080835203993218, 0,
      0.85065080835203993218, 0, 0.525731112119133606025,
      0, -0.525731112119133606025, 0.85065080835203993218 },
    {-0.85065080835203993218, 0, -0.525731112119133606025,
      0, -0.525731112119133606025, -0.85065080835203993218,
      -0.525731112119133606025, -0.85065080835203993218, 0},
    {-0.85065080835203993218, 0, 0.525731112119133606025,
      -0.525731112119133606025, -0.85065080835203993218, 0,
      0, -0.525731112119133606025, 0.85065080835203993218}
  };

  /* The unit charge is located at the center of the sphere with some raduis */
  sys->nChar = 1;
  CALLOC(sys->pos, 3*sys->nChar, double, ON, APVE);
  CALLOC(sys->chr, sys->nChar, double, ON, APVE);
  sys->pos[0] = 0.0;
  sys->pos[1] = 0.0;
  sys->pos[2] = 0.0;
  sys->chr[0] = 50.0;
  *nPanels = 20;
  nSurf = 0;

  for ( cnt=0; cnt<20; cnt++ ) {
    /* Allocate panel struct to fill in. */
    if(pnlList == NULL) {
      CALLOC(pnlList, 1, panel, ON, APVE);
      pnl = pnlList;
    }
    else {
      CALLOC(pnl->next, 1, panel, ON, APVE);
      pnl = pnl->next;
    }

    /* Fill in corners. */
    for(j=0;j<3;j++) {
      for(k=0;k<3;k++) {
        pnl->vtx[j][k] = vtx[cnt][j*3+k];
      }
    }
    pnl->shape = 3;
    pnl->nSurf = nSurf++;
    //printf("mkIco%d\n",pnl->nSurf);
  }

  for ( level=1; level<=lev; level++ ){
    unifRefine(pnlList, nPanels, 1.0);
  }

  return pnlList;
} /* mkIco */
/*
 * maps input panles on the ellipsoid
 *      (x/a)^2 + (y/b)^2 + (z/c)^2 = 1
 *
 */
void liftOntoEllipsoid(panel *pnlList, double a, double b, double c, int nSurf) {
  panel *pnl;
  double len;
  int i;

  for ( pnl=pnlList; pnl != NULL; pnl=pnl->next ) {
    for ( i=0; i<pnl->shape; i++ ) {
      /* this gives a better projection */
      len = SQR(pnl->vtx[i][0]) + SQR(pnl->vtx[i][1]) + SQR(pnl->vtx[i][2]);
      len = sqrt(len);
      pnl->vtx[i][0] *= a/len;
      pnl->vtx[i][1] *= b/len;
      pnl->vtx[i][2] *= c/len;
    }
  }

/*
  FILE *fp;
  char fname[256];
  sprintf(fname,"triangles.txt");
  fp=fopen(fname,"w");
  for (pnl=inputLst; pnl!=NULL; pnl=pnl->next) {
    fprintf(fp,"%d %f %f %f\n",i++,pnl->vtx[0][0],pnl->vtx[0][1],pnl->vtx[0][2]);
    fprintf(fp,"%d %f %f %f\n",i++,pnl->vtx[1][0],pnl->vtx[1][1],pnl->vtx[1][2]);
    fprintf(fp,"%d %f %f %f\n",i++,pnl->vtx[2][0],pnl->vtx[2][1],pnl->vtx[2][2]);
  }
  fclose(fp);
*/
} /* liftOntoEllipsoid */
