#include "mshmet.h"
#include <float.h>


int scaleMesh(pMesh mesh,pSol sol) {
  pPoint    ppt;
  Info     *info;
  double    dd;
  int       i,k,iadr;

  /* compute bounding box */
  info = &mesh->info;
  for (i=0; i<mesh->dim; i++) {
    info->min[i] =  FLT_MAX;
    info->max[i] = -FLT_MAX;
  }

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    for (i=0; i<mesh->dim; i++) {
      if ( ppt->c[i] > info->max[i] )  info->max[i] = ppt->c[i];
      if ( ppt->c[i] < info->min[i] )  info->min[i] = ppt->c[i];
    }
  }
  info->delta = 0.0;
  for (i=0; i<mesh->dim; i++) {
    dd = fabs(info->max[i]-info->min[i]);
    if ( dd > info->delta )  info->delta = dd;
  }
  if ( info->delta < EPS1 ) {
    fprintf(stdout,"  ## Unable to scale mesh\n");
    return(0);
  }

  /* compute solmax */
  if ( mesh->info.ls ) {
    sol->umin =  1.0e20;
    sol->umax = -1.0e20;
    for (k=1; k<=sol->np; k++) {
      iadr = (k-1) * sol->size + 1;
      for (i=0; i<sol->type; i++) {
        if ( sol->typtab[i] == GmfSca ) {
		      sol->umin = MS_MIN(sol->sol[iadr],sol->umin);
		      sol->umax = MS_MAX(sol->sol[iadr],sol->umax);
          iadr++;
        }
	    }
	  }
	}
	
  return(1);
}


int unscaleMesh(pMesh mesh,pSol sol) {
  pPoint     ppt;
  Info      *info;
  double     dd;
  int        i,k;

  /* de-normalize coordinates */
  info = &mesh->info;
  dd   = info->delta / (double)PRECI;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    for (i=0; i<mesh->dim; i++)
      ppt->c[i] = ppt->c[i] * dd + info->min[i];
  }

  return(1);
}



