#include "mshmet.h"


/* return unique sol is at vertex ip, convert vector-tensor into scalar */
double getSol_3d(pSol sol,int ip,int is) {
  double  ux,uy,uz,u,m[6],lambda[3],vp[3][3];
  int     i,iadr;

  u = 0.0;
  iadr = (ip-1)*sol->size + is+1;
  switch(sol->typtab[is]) {
  case GmfSca:
    u = sol->sol[iadr];
    break;
  case GmfVec:
    ux = sol->sol[iadr];
    uy = sol->sol[iadr+1];
    uz = sol->sol[iadr+2];
    u  = sqrt(ux*ux + uy*uy + uz*uz);
    break;
  case GmfSymMat:
    for (i=0; i<6; i++)  m[i] = sol->sol[iadr+i];
    eigenv(1,m,lambda,vp);
    u = MS_MAX(lambda[0],MS_MAX(lambda[1],lambda[2]));
  }

  return(u);
}


double getSol_2d(pSol sol,int ip,int is) {
  double  ux,uy,u,m[3],lambda[2],vp[2][2];
  int     i,iadr;

  u = 0.0;
  iadr = (ip-1)*sol->size + is+1;
  switch(sol->typtab[is]) {
  case GmfSca:
    u = sol->sol[iadr];
    break;
  case GmfVec:
    ux = sol->sol[iadr];
    uy = sol->sol[iadr+1];
    u  = sqrt(ux*ux + uy*uy);
    break;
  case GmfSymMat:
    for (i=0; i<3; i++)  m[i] = sol->sol[iadr+i];
    eigen2(m,lambda,vp);
    u = MS_MAX(lambda[0],lambda[1]);
  }

  return(u);
}


/* compute gradient via least square approx 
   in:  ip = vertex num   is = sol num 
   out: grd[3] gradient */
int gradLS_3d(pMesh mesh,pSol sol,int ip,int is) {
  pPoint    p0,p1;
	pTetra    pt;
  double   *grd,a[6],b[3],ax,ay,az,du,u,u1,aa,bb,cc,dd,ee,ff;
  int       i,iadr,nsdep,ip1,list[LONMAX+2],ilist;

  p0     = &mesh->point[ip];
  p0->nv = 0;
  nsdep  = p0->s;
	assert(nsdep);

	i  = 0;
	pt = &mesh->tetra[nsdep];
  if ( pt->v[1] == ip )      i = 1;
  else if ( pt->v[2] == ip ) i = 2;
  else if ( pt->v[3] == ip ) i = 3;

  ilist = boulep(mesh,nsdep,i,list);
	if ( ilist < 1 )  return(0);

  memset(a,0,6*sizeof(double));
  memset(b,0,3*sizeof(double));
  u = getSol(sol,ip,is);

	iadr = (ip-1)*sol->dim + 1;
	grd  = &sol->grd[iadr];

  /* Gradient: ui = u + <grad u,PPi> */
	for (i=1; i<=ilist; i++) {
		ip1 = list[i];
		p1  = &mesh->point[ip1];
    u1  = getSol(sol,ip1,is);
    du  = u1 - u;

    /* M = At*A symmetric definite positive */
    ax = p1->c[0] - p0->c[0];
    ay = p1->c[1] - p0->c[1];
    az = p1->c[2] - p0->c[2];
    a[0] += ax * ax;
    a[1] += ax * ay;
    a[2] += ax * az;
    a[3] += ay * ay;
    a[4] += ay * az;
    a[5] += az * az;

    /* b = At*du */
    b[0] += ax * du;
    b[1] += ay * du;
    b[2] += az * du;
    p0->nv++;
	}

  /* solution of A(3,3)*grad(1,3) = b(1,3) */ 
	aa = a[3]*a[5] - a[4]*a[4];
  bb = a[4]*a[2] - a[1]*a[5];
  cc = a[1]*a[4] - a[2]*a[3];
  du = aa*a[0] + bb*a[1] + cc*a[2];
  if ( fabs(du) == 0.0 ) {
    fprintf(stdout," Invalid matrix (%d). Exit\n",ip);
    return(-1);
  }
  du = 1.0 / du;
  dd = a[0]*a[5] - a[2]*a[2];
  ee = a[1]*a[2] - a[0]*a[4];
  ff = a[0]*a[3] - a[1]*a[1];

  grd[0] = (aa*b[0] + bb*b[1] + cc*b[2]) * du;
  grd[1] = (bb*b[0] + dd*b[1] + ee*b[2]) * du;
  grd[2] = (cc*b[0] + ee*b[1] + ff*b[2]) * du;

  /* normalize */
  if ( mesh->info.ls ) {
    du = grd[0]*grd[0] + grd[1]*grd[1] + grd[2]*grd[2];
    if ( du > EPS1 ) {
      du      = 1.0 / sqrt(du);
      grd[0] *= du;
      grd[1] *= du;
      grd[2] *= du;
    }
  }

  return(1);
}


int gradLS_2d(pMesh mesh,pSol sol,int ip,int is) {
  pPoint    p0,p1;
  pTria     pt;
  double   *grd,a[3],b[2],ax,ay,du,u,u1;
  int       i,iadr,nsdep,ip1,list[LONMAX+2],ilist;

  p0     = &mesh->point[ip];
  p0->nv = 0;
  nsdep  = p0->s;
	assert(nsdep);

  pt   = &mesh->tria[nsdep];
  if ( pt->v[0] == ip )      i = 0;
  else if ( pt->v[1] == ip ) i = 1;
  else                       i = 2;

  ilist = boulep(mesh,nsdep,i,list);
	if ( ilist < 1 )  return(0);

  memset(a,0,3*sizeof(double));
  memset(b,0,2*sizeof(double));
  u = getSol(sol,ip,is);

	iadr = (ip-1)*sol->dim + 1;
	grd  = &sol->grd[iadr];

  /* Gradient: ui = u + <grad u,PPi> */
  for (i=1; i<=ilist; i++) {
    ip1 = list[i];
    p1  = &mesh->point[ip1];
    u1  = getSol(sol,ip1,is);

    /* M = At*A symmetric definite positive */
    ax    = p1->c[0] - p0->c[0];
    ay    = p1->c[1] - p0->c[1];
    a[0] += ax * ax;
    a[1] += ax * ay;
    a[2] += ay * ay;

    /* b = At*du */
    du    = u1 - u;
    b[0] += ax * du;
    b[1] += ay * du;
    p0->nv++;
  }

  /* solution of A(2,2)*grad(1,2) = b(1,2) */ 
  du = a[0]*a[2] - a[1]*a[1];
  if ( fabs(du) == 0.0 ) {
    fprintf(stdout," Invalid matrix (%E). Exit\n",du);
    return(0);
  }
  du  = 1.0 / du;
  grd[0] = (a[2]*b[0] - a[1]*b[1]) * du;
  grd[1] = (a[0]*b[1] - a[1]*b[0]) * du;

  /* normalize */
  if ( mesh->info.ls ) {
    du = grd[0]*grd[0] + grd[1]*grd[1];
    if ( du > EPS1 ) {
      du      = 1.0 / sqrt(du);
      grd[0] *= du;
      grd[1] *= du;
    }
  }

  return(1);
}


/* for surface mesh */
int gradLS_s(pMesh mesh,pSol sol,int ip,int is) {
  pPoint    p0,p1;
	pTria     pt;
  double   *grd,*nn,a[6],b[3],ax,ay,az,du,u,u1,aa,bb,cc,dd,ee,ff;
  int       i,iadr,nsdep,ip1,list[LONMAX+2],ilist;

  p0     = &mesh->point[ip];
  p0->nv = 0;
  nsdep  = p0->s;
	assert(nsdep);

	i  = 0;
	pt = &mesh->tria[nsdep];
  if ( pt->v[1] == ip )      i = 1;
  else if ( pt->v[2] == ip ) i = 2;

  ilist = boulep(mesh,nsdep,i,list);
	if ( ilist < 1 )  return(0);

  memset(a,0,6*sizeof(double));
  memset(b,0,3*sizeof(double));
  u = getSol(sol,ip,is);

	iadr = (ip-1)*sol->dim + 1;
	grd  = &sol->grd[iadr];
	nn   = &sol->nn[iadr];

  /* Gradient: ui = u + <grad u,PPi> */
	for (i=1; i<=ilist; i++) {
		ip1 = list[i];
		p1  = &mesh->point[ip1];
    u1  = getSol(sol,ip1,is);
    du  = u1 - u;

    ax = p1->c[0] - p0->c[0];
    ay = p1->c[1] - p0->c[1];
    az = p1->c[2] - p0->c[2];

    /* M = At*A symmetric definite positive */
    a[0] += ax * ax;
    a[1] += ax * ay;
    a[2] += ax * az;
    a[3] += ay * ay;
    a[4] += ay * az;
    a[5] += az * az;

    /* b = At*du */
    b[0] += ax * du;
    b[1] += ay * du;
    b[2] += az * du;
    p0->nv++;
	}

  /* solution of A(3,3)*grad(1,3) = b(1,3) */ 
	aa = a[3]*a[5] - a[4]*a[4];
  bb = a[4]*a[2] - a[1]*a[5];
  cc = a[1]*a[4] - a[2]*a[3];
  du = aa*a[0] + bb*a[1] + cc*a[2];
  if ( fabs(du) == 0.0 ) {
    fprintf(stdout," Invalid matrix (%E). Exit\n",du);
    return(0);
  }
  du = 1.0 / du;
  dd = a[0]*a[5] - a[2]*a[2];
  ee = a[1]*a[2] - a[0]*a[4];
  ff = a[0]*a[3] - a[1]*a[1];

  grd[0] = (aa*b[0] + bb*b[1] + cc*b[2]) * du;
  grd[1] = (bb*b[0] + dd*b[1] + ee*b[2]) * du;
  grd[2] = (cc*b[0] + ee*b[1] + ff*b[2]) * du;

  /* align w normal */
  du = grd[0]*grd[0] + grd[1]*grd[1] + grd[2]*grd[2];
  grd[0] = du * nn[0];
  grd[1] = du * nn[1];
  grd[2] = du * nn[2];

  /* normalize */
  if ( mesh->info.ls ) {
    du = grd[0]*grd[0] + grd[1]*grd[1] + grd[2]*grd[2];
    if ( du > EPS1 ) {
      du      = 1.0 / sqrt(du);
      grd[0] *= du;
      grd[1] *= du;
      grd[2] *= du;
    }
  }

	return(1);
}

