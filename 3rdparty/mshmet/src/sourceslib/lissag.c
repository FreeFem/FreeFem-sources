#include "mshmet.h"

#define lambda 0.04
#define mu     0.03


int lissag_3d(pMesh mesh,pSol sol,int is,int it) {
  pPoint    p0,p1;
	pTetra    pt;
  double   *nu,u,v;
  int       i,j,k,iadr,nsdep,ip1,nb,list[LONMAX+2],ilist;

	nu = calloc(mesh->np+1,sizeof(double));
	assert(nu);

	for (j=1; j<=it; j++) {

    /* 1st stage: new value */
    for (k=1; k<=mesh->np; k++) {
  		p0    = &mesh->point[k];
      nsdep = p0->s;
      assert(nsdep);
  
  	  pt = &mesh->tetra[nsdep];
      i  = 0;
      if ( pt->v[1] == k )      i = 1;
      else if ( pt->v[2] == k ) i = 2;
      else if ( pt->v[3] == k ) i = 3;
  
  		ilist = boulep(mesh,nsdep,i,list);
      u  = getSol(sol,k,is);
  		v  = 0.0;
  		nb = 0;
  	
      for (i=1; i<=ilist; i++) {
        ip1 = list[i];
        p1  = &mesh->point[ip1];
        v  += getSol(sol,ip1,is);
  			nb++;
      }
  		v     = v / nb;
  		nu[k] = u + lambda*(v-u);
    }
  
    /* 2nd stage: update value */
    for (k=1; k<=mesh->np; k++) {
  		p0    = &mesh->point[k];
      nsdep = p0->s;
  
  	  pt = &mesh->tetra[nsdep];
      i  = 0;
      if ( pt->v[1] == k )      i = 1;
      else if ( pt->v[2] == k ) i = 2;
      else if ( pt->v[3] == k ) i = 3;
  
  		ilist = boulep(mesh,nsdep,i,list);
      u  = nu[k];
  		v  = 0.0;
  		nb = 0;
      for (i=1; i<=ilist; i++) {
        ip1 = list[i];
        p1  = &mesh->point[ip1];
        v  += nu[ip1];
  			nb++;
      }
  		v = v / nb;
  
      iadr = (k-1)*sol->size + is+1;
  		sol->sol[iadr] = u - mu*(u-v);
    }

  	if ( j < it )  memset(nu,0,(mesh->np+1)*sizeof(double));
  }

	free(nu);
	return(1);
}


/* regularisation of sol is */
int lissag_2d(pMesh mesh,pSol sol,int is,int it) {
  pPoint    p0,p1;
  pTria     pt;
  double   *nu,u,v;
  int       i,j,k,iadr,nsdep,ip1,nb,list[LONMAX+2],ilist;

	nu = calloc(mesh->np+1,sizeof(double));
	assert(nu);

	for (j=1; j<=it; j++) {

    /* 1st stage: new value */
    for (k=1; k<=mesh->np; k++) {
		  p0    = &mesh->point[k];
      nsdep = p0->s;
      assert(nsdep);

      pt = &mesh->tria[nsdep];
      i  = 0;
      if ( pt->v[1] == k )      i = 1;
      else if ( pt->v[2] == k ) i = 2;
  
  		ilist = boulep(mesh,nsdep,i,list);
  		if ( ilist < 1 )  return(0);
      u  = getSol(sol,k,is);
  		v  = 0.0;
  		nb = 0;
  
      for (i=1; i<=ilist; i++) {
        ip1 = list[i];
        p1  = &mesh->point[ip1];
        v  += getSol(sol,ip1,is);
  			nb++;
      }
  
  		v = v / nb;
  		nu[k] = u + lambda*(v-u);
    }
  
    /* 2nd stage: update value */
    for (k=1; k<=mesh->np; k++) {
  		p0    = &mesh->point[k];
      nsdep = p0->s;
  
      pt = &mesh->tria[nsdep];
      i  = 0;
      if ( pt->v[1] == k )      i = 1;
      else if ( pt->v[2] == k ) i = 2;
  
  		ilist = boulep(mesh,nsdep,i,list);
  		if ( ilist < 1 )  return(0);
      u  = nu[k];
  		v  = 0.0;
  		nb = 0;
  		for (i=1; i<=ilist; i++) {
        ip1 = list[i];
        p1  = &mesh->point[ip1];
        v  += nu[ip1];
  			nb++;
      }
  		v = v / nb;
  
      iadr = (k-1)*sol->size + is+1;
  		sol->sol[iadr] = u - mu*(u-v);
    }

		if ( j < it )  memset(nu,0,(mesh->np+1)*sizeof(double));
  }
 
	free(nu);
  return(1);
}

