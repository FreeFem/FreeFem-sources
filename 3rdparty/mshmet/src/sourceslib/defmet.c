#include "mshmet.h"


/* create 3d metric */
int defmet_3d(pMesh mesh,pSol sol,int is) {
  double    lmax,lambda[3],vp[3][3],m[6],mr[6],hlmin,hlmax;
  double   *tmet,*hes,dd;
  int       i,j,k,l,iadr,iasol;

  hlmin = 1.0 / (mesh->info.hmax*mesh->info.hmax);
  hlmax = 1.0 / (mesh->info.hmin*mesh->info.hmin);

  for (k=1; k<=mesh->np; k++) {
		iadr = (k-1)*6 + 1;
		hes  = &sol->hes[iadr];
    if ( !eigenv(1,hes,lambda,vp) ) {
      if ( mesh->info.imprim < 0 ) 
        fprintf(stdout,"  %%%% Unable to define metric %d.\n",k);
      return(0);
    }

    /* truncation */
		lmax = 0.0;
		for (i=0; i<3; i++) {
      dd        = MS_MAX(fabs(lambda[i]),hlmin);
      lambda[i] = MS_MIN(dd,hlmax);
			if ( lambda[i] > lmax )  lmax = lambda[i];
    }

    if ( mesh->info.iso ) {
      tmet = &sol->met[k];
      if ( is == 0 && !mesh->info.metric )
        tmet[0] = 1.0 / sqrt(lmax);
      else
        tmet[0] = MS_MIN(1.0/sqrt(lmax),tmet[0]);
    }

    else {
      iasol = (k-1)*6 + 1;
      tmet  = &sol->met[iasol];

      /* compute Mat = B M Bt */
      for (l=0,i=0; i<3; i++)
        for (j=i; j<3; j++)
          m[l++] = lambda[0]*vp[0][i]*vp[0][j] + lambda[1]*vp[1][i]*vp[1][j] \
                 + lambda[2]*vp[2][i]*vp[2][j];

      if ( is == 0 && !mesh->info.metric ) {
				for (i=0; i<6; i++)  tmet[i] = m[i];
      }
      else {
        if ( !redsim(tmet,m,mr) ) {
          fprintf(stdout,"  %%%% Metric problem (point: %d)\n",k);
          if ( mesh->info.ddebug ) { 
            fprintf(stdout,"  lambda: %e %e %e\n",lambda[0],lambda[1],lambda[2]);
					  fprintf(stdout,"  mat %e %e %e\n%e %e\n %e\n",m[0],m[1],m[2],m[3],m[4],m[5]);
          }
          return(0);
	      } 
        else
				  for (i=0; i<6; i++)  tmet[i] = mr[i];
      }
    }
  }

	return(1);
}


int defmet_2d(pMesh mesh,pSol sol,int is) {
  double    lambda[2],vp[2][2],mat[3],mr[3],dd,hlmin,hlmax;
  double   *tmet,*hes;
  int       i,k,iasol,iadr;

  hlmin = 1.0 / (mesh->info.hmax*mesh->info.hmax);
  hlmax = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  for (k=1; k<=mesh->np; k++) {
		iadr = (k-1)*3 + 1;
		hes  = &sol->hes[iadr];
    if ( !eigen2(hes,lambda,vp) ) {
      if ( mesh->info.imprim < 0 ) 
        fprintf(stdout,"  %%%% Unable to define metric %d.\n",k);
      return(0);
    }

    /* truncation */
    dd        = MS_MAX(fabs(lambda[0]),hlmin);
    lambda[0] = MS_MIN(dd,hlmax);
    dd        = MS_MAX(fabs(lambda[1]),hlmin);
    lambda[1] = MS_MIN(dd,hlmax);

    if ( mesh->info.iso ) {
      tmet = &sol->met[k];
      if ( lambda[0] < lambda[1] ) 
        lambda[0] = lambda[1];
      else
        lambda[1] = lambda[0];
      if ( is == 0 && !mesh->info.metric )
        tmet[0] = 1.0 / sqrt(lambda[0]);
      else
        tmet[0] = MS_MIN(1.0/sqrt(lambda[0]),tmet[0]);
    }

    else {
      iasol = (k-1)*3 + 1;
      tmet  = &sol->met[iasol];

      /* compute Mat = B M Bt */
      mat[0] = lambda[0]*vp[0][0]*vp[0][0] + lambda[1]*vp[1][0]*vp[1][0];
      mat[1] = lambda[0]*vp[0][0]*vp[0][1] + lambda[1]*vp[1][0]*vp[1][1];
      mat[2] = lambda[0]*vp[0][1]*vp[0][1] + lambda[1]*vp[1][1]*vp[1][1];

      if ( is == 0 && !mesh->info.metric ) {
        tmet[0] = mat[0];
	      tmet[1] = mat[1];
	      tmet[2] = mat[2];
      }
      else {  /* intersection */
        redsim(tmet,mat,mr);
				for (i=0; i<3; i++)  tmet[i] = mr[i];
      }
    }
  }

  return(1);
}


int defmet_s(pMesh mesh,pSol sol,int is) {
  double    lmax,lambda[3],vp[3][3],m[6],mr[6],hlmin,hlmax;
  double   *tmet,*hes,*nn,dd,d1,d2,d3;
  int       i,j,k,l,iadr,iasol;

  hlmin = 1.0 / (mesh->info.hmax*mesh->info.hmax);
  hlmax = 1.0 / (mesh->info.hmin*mesh->info.hmin);

  for (k=1; k<=mesh->np; k++) {
		iadr = (k-1)*6 + 1;
		hes  = &sol->hes[iadr];
    if ( !eigenv(1,hes,lambda,vp) ) {
      if ( mesh->info.imprim < 0 ) 
        fprintf(stdout,"  %%%% Unable to define metric %d.\n",k);
      return(0);
    }

		iadr = (k-1)*sol->dim + 1;
		nn = &sol->nn[iadr];

		d1 = fabs(nn[0]*vp[0][0] + nn[1]*vp[0][1] + nn[2]*vp[0][2]);
		d2 = fabs(nn[0]*vp[1][0] + nn[1]*vp[1][1] + nn[2]*vp[1][2]);
		d3 = fabs(nn[0]*vp[2][0] + nn[1]*vp[2][1] + nn[2]*vp[2][2]);
		
    if ( d1 > d2 ) {
      if ( d1 > d3 ) {
			  //memcpy(vp[0],nn,3*sizeof(double));
				lambda[0] = hlmin;
			}
      else {
        //memcpy(vp[2],nn,3*sizeof(double));
				lambda[2] = hlmin;
			}
		}
		else if ( d2 > d3 ) {
		  //memcpy(vp[1],nn,3*sizeof(double));
			lambda[1] = hlmin;
		}
		else {
		  //memcpy(vp[2],nn,3*sizeof(double));
		  lambda[2] = hlmin;
		}

    /* truncation */
		lmax = 0.0;
		for (i=0; i<3; i++) {
      dd        = MS_MAX(fabs(lambda[i]),hlmin);
      lambda[i] = MS_MIN(dd,hlmax);
			if ( lambda[i] > lmax )  lmax = lambda[i];
    }
  
    if ( mesh->info.iso ) {
      tmet = &sol->met[k];
      if ( is == 0 && !mesh->info.metric )
        tmet[0] = 1.0 / sqrt(lmax);
      else
        tmet[0] = MS_MIN(1.0/sqrt(lmax),tmet[0]);
    }

    else {
      iasol = (k-1)*6 + 1;
      tmet  = &sol->met[iasol];

      /* compute Mat = B M Bt */
      for (l=0,i=0; i<3; i++)
        for (j=i; j<3; j++)
          m[l++] = lambda[0]*vp[0][i]*vp[0][j] + lambda[1]*vp[1][i]*vp[1][j] \
                 + lambda[2]*vp[2][i]*vp[2][j];

      if ( is == 0 && !mesh->info.metric ) {
				for (i=0; i<6; i++)  tmet[i] = m[i];
      }
      else {
        redsim(tmet,m,mr);
				for (i=0; i<6; i++)  tmet[i] = mr[i];
      }
    }
  }

  return(1);
}
