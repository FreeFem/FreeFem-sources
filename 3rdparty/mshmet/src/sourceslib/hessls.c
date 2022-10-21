#include "mshmet.h"


/* solve Mx=b using Gauss */
static int gauss(int n,double m[n][n],double *x,double *b,char dbg) {
  double    nn,dd;
  int       i,j,k,ip;

  nn = m[0][0];
  for (i=0; i<n; i++) 
    for (j=0; j<n; j++)
      nn = MS_MAX(nn,m[i][j]);

  if ( fabs(nn) < EPS1 ) {
    if ( dbg )  fprintf(stdout,"  %%%% Null matrix\n");
    return(0);
  }
  nn = 1.0 / nn;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++)
      m[i][j] *= nn;
    b[i] *= nn; 
  }

  /* partial pivoting, column j */
  for (j=0; j<n-1; j++) {
    /* find line pivot */
    ip = j;
    for (i=j+1; i<n; i++)
      if ( fabs(m[i][j]) > fabs(m[ip][j]) )  ip = i;

    /* swap i <-> ip */
    if ( i != ip ) {
      for (k=j; k<n; k++) {
        dd       = m[j][k];
        m[j][k]  = m[ip][k];
        m[ip][k] = dd;
      }
      dd    = b[j];
      b[j]  = b[ip];
      b[ip] = dd;
    }
    if ( fabs(m[j][j]) < EPS1 ) {
      if ( dbg )  fprintf(stdout,"  %%%% Null pivot: %e.\n",m[j][i]);
      return(0);
    }

    /* elimination */
    for (i=j+1; i<n; i++) {
      dd = m[i][j] / m[j][j];
      m[i][j] = 0.0;
      for (k=j+1; k<n; k++)
        m[i][k] -= dd * m[j][k]; 
      b[i] -= dd * b[j];
    }
  }

  if ( fabs(m[n-1][n-1]) < EPS1 ) {
    if ( dbg )  fprintf(stdout,"  %%%% Null pivot.\n");
    return(0);
  }

  x[n-1] = b[n-1] / m[n-1][n-1];
  for (i=n-2; i>=0; i--) {
    dd = 0.0;
    for (j=i+1; j<n; j++)
      dd += m[i][j] * x[j];
    x[i] = (b[i] - dd) / m[i][i];
  }

  if ( dbg ) {
    for (i=0; i<n; i++) {
      dd = 0.;
      for (j=0; j<n; j++)  dd += m[i][j] * x[j];
      if ( fabs(dd-b[i]) > EPS ) {
        fprintf(stdout,"  Ax[%d] = %f   b[%d] = %E\n",i,dd,i,b[i]);
        exit(1);
      }
    }
  }

  return(1);
}


/* least square approximation of hessian */
int hessLS_3d(pMesh mesh,pSol sol,int ip,int is) {
  pPoint    p0,p1;
  pTetra    pt;
  double    a[6],b[6],m[6][6],ax,ay,az;
  double   *grd,*hes,du,u,u1;
  int       i,j,l,iadr,ier,nsdep,ip1,list[LONMAX+2],ilist;

  p0 = &mesh->point[ip];
  if ( p0->nv < 6 )  return(-1);
  nsdep = p0->s;
  assert(nsdep);

  i  = 0;
  pt = &mesh->tetra[nsdep];
  if ( pt->v[1] == ip )      i = 1;
  else if ( pt->v[2] == ip ) i = 2;
  else if ( pt->v[3] == ip ) i = 3;

  ilist = boulep(mesh,nsdep,i,list);
  if ( ilist < 1 )  return(0);

  memset(b,0,6*sizeof(double));
  memset(m,0,36*sizeof(double));
  u = getSol(sol,ip,is);

  iadr = (ip-1)*sol->dim + 1;
  grd  = &sol->grd[iadr];

  /* Hessian: Ui = U + <gradU,PPi> + 0.5*(tPPi.Hess.PPi) */
  for (i=1; i<=ilist; i++) {
    ip1 = list[i];
    p1  = &mesh->point[ip1];
    u1  = getSol(sol,ip1,is);

    ax = p1->c[0] - p0->c[0];
    ay = p1->c[1] - p0->c[1];
    az = p1->c[2] - p0->c[2];

    a[0] = 0.5 * ax*ax;
    a[1] =       ax*ay;
    a[2] =       ax*az;
    a[3] = 0.5 * ay*ay;
    a[4] =       ay*az;
    a[5] = 0.5 * az*az;
    du   = (u1-u) - (ax*grd[0] + ay*grd[1] + az*grd[2]);

    /* M = At*A symmetric positive definite */
    for (j=0; j<6; j++)
      for (l=j; l<6; l++) {
        m[j][l] += a[j]*a[l];
        m[l][j] += a[l]*a[j];
      }

    /* c = At*b */
    for (j=0; j<6; j++)
      b[j] += a[j]*du;
  }

  /* solve m(6,6)*x(6) = b(6) */
  iadr = 6*(ip-1) + 1;
  hes  = &sol->hes[iadr];

  ier = gauss(6,m,hes,b,mesh->info.ddebug);
  if ( !ier ) {
	  if ( mesh->info.ddebug )  fprintf(stdout," Ill cond'ed matrix (%d, %d).\n",ip,ier);
    return(-1);
  }

  p0->h = 1;
  return(1);
}


/* least square approximation of hessian */
int hessLS_2d(pMesh mesh,pSol sol,int ip,int is) {
  pPoint    p0,p1;
  pTria     pt;
  double    a[3],ma[6],mb[3],ax,ay;
  double    det,aa,bb,cc,dd,ee,ff;
  double   *grd,*hes,u,u1;
  int       i,iadr,nsdep,ip1,list[LONMAX+2],ilist;

  p0 = &mesh->point[ip];
  if ( p0->nv < 3 )  return(-1);
  nsdep = p0->s;
  assert(nsdep);

  pt   = &mesh->tria[nsdep];
  if ( pt->v[0] == ip )      i = 0;
  else if ( pt->v[1] == ip ) i = 1;
  else                       i = 2;

  ilist = boulep(mesh,nsdep,i,list);
  if ( ilist < 1 )  return(0);

  memset(ma,0,6*sizeof(double));
  memset(mb,0,3*sizeof(double));
  u = getSol(sol,ip,is);
  
  iadr = (ip-1)*sol->dim + 1;
  grd  = &sol->grd[iadr];

  /* Hessian: Ui = U + <gradU,PPi> + 0.5*(tPPi.Hess.PPi) */
  for (i=1; i<=ilist; i++) {
    ip1 = list[i];
    p1  = &mesh->point[ip1];
    u1  = getSol(sol,ip1,is);

    ax    = p1->c[0] - p0->c[0];
    ay    = p1->c[1] - p0->c[1];
    a[0]  = 0.5 * ax*ax;
    a[1]  =       ax*ay;
    a[2]  = 0.5 * ay*ay;
    dd    = (u1-u) - (ax*grd[0] + ay*grd[1]);

    /* M = At*A symmetric definite positive */
    ma[0] += a[0]*a[0];
    ma[1] += a[0]*a[1];
    ma[2] += a[1]*a[1];
    ma[3] += a[0]*a[2];
    ma[4] += a[1]*a[2];
    ma[5] += a[2]*a[2];

    /* c = At*b */
    mb[0] += a[0]*dd;
    mb[1] += a[1]*dd;
    mb[2] += a[2]*dd;
  }

  if ( fabs(ma[0]) < EPS1 ) {
    if ( mesh->info.ddebug )  fprintf(stdout," Ill cond'ed matrix (%d, %E).\n",ip,ma[0]);
    return(-1);
  }

  /* direct solving */
  aa = ma[2]*ma[5] - ma[4]*ma[4];
  bb = ma[4]*ma[3] - ma[1]*ma[5];
  cc = ma[1]*ma[4] - ma[2]*ma[3];
  dd = ma[0]*aa + ma[1]*bb + ma[3]*cc;

  /* singular matrix */
  if ( fabs(dd) == 0.0 ) {
    if ( mesh->info.ddebug )  fprintf(stderr," Singular matrix (%d, %E).\n",ip,dd);
    return(-1);
  }
  det = 1.0 / dd;
  dd = ma[0]*ma[5] - ma[3]*ma[3];
  ee = ma[1]*ma[3] - ma[0]*ma[4];
  ff = ma[0]*ma[2] - ma[1]*ma[1];

  iadr = (ip-1)*3 + 1;
  hes  = &sol->hes[iadr];

  hes[0] = (mb[0]*aa + mb[1]*bb + mb[2]*cc) * det;
  hes[1] = (mb[0]*bb + mb[1]*dd + mb[2]*ee) * det;
  hes[2] = (mb[0]*cc + mb[1]*ee + mb[2]*ff) * det;

  if ( !p0->b && p0->nv == 4 )  return(-1);

  p0->h = 1;
  return(1);
}


int hessLS_s(pMesh mesh,pSol sol,int ip,int is) {
  pPoint    p0,p1;
  pTria     pt;
  double    a[6],b[6],m[6][6],ax,ay,az;
  double   *grd,*hes,du,u,u1;
  int       i,j,l,iadr,ier,nsdep,ip1,list[LONMAX+2],ilist;

  p0 = &mesh->point[ip];
  if ( p0->nv < 3 )  return(-1);
  nsdep = p0->s;
  assert(nsdep);

  i  = 0;
  pt = &mesh->tria[nsdep];
  if ( pt->v[1] == ip )      i = 1;
  else if ( pt->v[2] == ip ) i = 2;

  ilist = boulep(mesh,nsdep,i,list);
  if ( ilist < 1 )  return(0);

  memset(b,0,6*sizeof(double));
  memset(m,0,36*sizeof(double));
  u = getSol(sol,ip,is);

  iadr = (ip-1)*sol->dim + 1;
  grd  = &sol->grd[iadr];

  /* Hessian: Ui = U + <gradU,PPi> + 0.5*(tPPi.Hess.PPi) */
  for (i=1; i<=ilist; i++) {
    ip1 = list[i];
    p1  = &mesh->point[ip1];
    u1  = getSol(sol,ip1,is);

    ax = p1->c[0] - p0->c[0];
    ay = p1->c[1] - p0->c[1];
    az = p1->c[2] - p0->c[2];
 
    a[0] = 0.5 * ax*ax;
    a[1] =       ax*ay;
    a[2] =       ax*az;
    a[3] = 0.5 * ay*ay;
    a[4] =       ay*az;
    a[5] = 0.5 * az*az;
    du   = (u1-u) - (ax*grd[0] + ay*grd[1] + az*grd[2]);

    /* M = At*A symmetric positive definite */
    for (j=0; j<6; j++)
      for (l=j; l<6; l++) {
        m[j][l] += a[j]*a[l];
        m[l][j] += a[l]*a[j];
      }

    /* c = At*b */
    for (j=0; j<6; j++)
      b[j] += a[j]*du;
  }

  /* solve m(6,6)*x(6) = b(6) */
  iadr = (ip-1)*6 + 1;
  hes  = &sol->hes[iadr];

  ier = gauss(6,m,hes,b,mesh->info.ddebug);
  if ( !ier ) {
    if ( mesh->info.ddebug )  fprintf(stdout," Ill cond'ed matrix (%d, %d).\n",ip,ier);
    return(-1);
  }

  p0->h = 1;
  return(1);
}


/* average values of neighbors */
int avgval_3d(pMesh mesh,pSol sol,int ip) {
  pPoint    p0,p1;
  pTetra    pt;
  double   *hes,dd;
  int       i,j,iadr,nsdep,ip1,nb,list[LONMAX+2],ilist;

  p0    = &mesh->point[ip];
  nsdep = p0->s;
  assert(nsdep);

  i = 0;
  pt = &mesh->tetra[nsdep];
  if ( pt->v[1] == ip )      i = 1;
  else if ( pt->v[2] == ip ) i = 2;
  else if ( pt->v[3] == ip ) i = 3;

  ilist = boulep(mesh,nsdep,i,list);
  if ( ilist < 1 )  return(0);

  iadr = (ip-1)*6 + 1;
  hes  = &sol->hes[iadr];
  memset(hes,0,6*sizeof(double));
  nb = 0;
  for (i=1; i<=ilist; i++) {
    ip1 = list[i];
    p1  = &mesh->point[ip1];
    if ( p1->h > 0 ) {
      iadr = (ip1-1)*6 + 1;
      for (j=0; j<6; j++)
        hes[j] += sol->hes[iadr+j];
      nb++;
    }
  }

  if ( nb > 0 ) {
    dd = 1.0 / (double)nb;
    for (j=0; j<6; j++)
      hes[j] = hes[j] * dd;
    p0->h = 1;
    return(1);
  }

  return(0);
}


int avgval_2d(pMesh mesh,pSol sol,int ip) {
  pPoint    p0,p1;
  pTria     pt;
  double   *hes,dd;
  int       i,j,iadr,nsdep,ip1,nb,list[LONMAX+2],ilist;

  p0    = &mesh->point[ip];
  nsdep = p0->s;
  assert(nsdep);

  pt   = &mesh->tria[nsdep];
  if ( pt->v[0] == ip )      i = 0;
  else if ( pt->v[1] == ip ) i = 1;
  else                       i = 2;

  ilist = boulep(mesh,nsdep,i,list);
  if ( ilist < 1 )  return(0);

  iadr = (ip-1)*3 + 1;
  hes  = &sol->hes[iadr];
  memset(hes,0,3*sizeof(double));
  nb = 0;
  for (i=1; i<=ilist; i++) {
    ip1 = list[i];
    p1    = &mesh->point[ip1];
    if ( p1->h > 0 ) {
      iadr = (ip1-1)*3 +1;
      for (j=0; j<3; j++)
        hes[j] += sol->hes[iadr+j];
      nb++;
    }
  }

  if ( nb > 0 ) {
    dd = 1.0 / (double)nb;
    for (j=0; j<3; j++)
      hes[j] = hes[j] * dd;
    p0->h = 1;
    return(1);
  }

  return(0);
}


/* assign value of closest vertex */
int clsval_3d(pMesh mesh,pSol sol,int ip) {
  pPoint    p0,p1;
  pTetra    pt;
  double   *hes;
  int       i,j,iadr,nsdep,ip1;

  p0    = &mesh->point[ip];
  nsdep = p0->s;
  assert(nsdep);

  pt = &mesh->tetra[nsdep];
  iadr = (ip-1)*6 + 1;
  hes  = &sol->hes[iadr];
  for (i=0; i<4; i++) {
    ip1 = pt->v[i];
    p1  = &mesh->point[ip1];
    if ( p1->h ) {
      iadr = (ip1-1)*6 + 1;
      for (j=0; j<6; j++)
        hes[j] = sol->hes[iadr+j];
      p0->h = 1;
      return(1);
    }
  }

  return(0);
}


int clsval_2d(pMesh mesh,pSol sol,int ip) {
  pPoint    p0,p1;
  pTria     pt;
  double   *hes;
  int       i,j,iadr,nsdep,ip1;

  p0    = &mesh->point[ip];
  nsdep = p0->s;
  assert(nsdep);

  pt = &mesh->tria[nsdep];
  iadr = (ip-1)*3 + 1;
  hes  = &sol->hes[iadr];
  for (i=0; i<3; i++) {
    ip1 = pt->v[i];
    p1  = &mesh->point[ip1];
    if ( p1->h ) {
      iadr = (ip1-1)*3 + 1;
      for (j=0; j<3; j++)
        hes[j] = sol->hes[iadr+j];
      p0->h = 1;
      return(1);
    }
  }
  
  return(0);
}


int nrmhes_3d(pMesh mesh,pSol sol,int is) {
  pPoint   p0;
  Info     info;
  double   err,err1,errs,u,norm;
  int      i,k,iadr;

  info = mesh->info;

  if ( info.nnu > 0 || mesh->info.ls ) {
    for (k=1; k<=mesh->np; k++) {
      u = fabs(getSol(sol,k,is));
      sol->umax = MS_MAX(u,sol->umax);
    }
  }
  
  switch(info.nnu) {
  /* no norm */
  case 0:
    err = CTE3D / info.eps;
    for (k=1; k<=mesh->np; k++) {
      iadr = (k-1)*6 + 1;
      for (i=0; i<6; i++)  sol->hes[iadr+i] *= err;
    }
    break;

  /* relative value: M(u)= |H(u)| / (err*||u||_inf) */
  case 1:
  default:
    if ( sol->umax < 1e-30 )  return(1);
    err =  CTE3D / (info.eps * sol->umax);
    for (k=1; k<=mesh->np; k++) {
      iadr = (k-1)*6 + 1;
      for (i=0; i<6; i++)  sol->hes[iadr+i] *= err;
    }
    break;

  /* local norm: M(u)= |H(u)| / err*|u| */
  case 2:
    errs = sol->umax > 0.0 ? sol->umax*0.01 : 0.01;
    for (k=1; k<=mesh->np; k++) {
      u    = fabs(getSol(sol,k,is));
      err  = CTE3D / MS_MAX(errs,u);
      iadr = (k-1)*6 + 1;
      for (i=0; i<6; i++)  sol->hes[iadr+i] *= err;
    }
    break;

  /* local norm: M(u)= |H(u)| / err*(e1*|u|+e2*h*|du|) */
  case 3:
      puts("A CODER");
      exit(1);
/*
    maxu = 0.0;
    maxg = 0.0;
    for (k=1; k<=mesh->np; k++) {
      u    = fabs(getSol(sol,k,is));
      maxu = MS_MAX(maxu,u);
      norm = der[k].grd[0]*der[k].grd[0] + der[k].grd[1]*der[k].grd[1] + der[k].grd[2]*der[k].grd[2];
      norm = sqrt(norm);
      if ( norm > maxg )  maxg = norm;
    }
    if ( maxu == 0.0 )
      errs = 0.01;
    else
      errs = maxu * 0.01;
    err1 = maxg * 0.01;

    for (k=1; k<=mesh->np; k++) {
      p0    = &mesh->point[k];
      u     = getSol(sol,k,is);
      norm  = der[k].grd[0]*der[k].grd[0] + der[k].grd[1]*der[k].grd[1] + der[k].grd[2]*der[k].grd[2];
      norm  = sqrt(norm);
      norm  = MS_MAX(err1,norm);
      //norm *= p0->rins / p0->nv;
      err1  = MS_MAX(errs,fabs(u));
      
      /* variable weight /
      err1  = 0.01*err1 + 0.01*norm;
      maxu  = CTE3D / (info.eps*err1);
      for (i=0; i<6; i++)  der[k].hes[i] *= maxu;
    } 
*/
    break;
  }

  return(1);
}


int nrmhes_2d(pMesh mesh,pSol sol,int is) {
  pPoint   p0;
  Info     info;
  double   err,err1,errs,u,norm;
  int      i,k,iadr;

  info = mesh->info;

  if ( info.nnu > 0 || mesh->info.ls ) {
    sol->umax = 0.0;
    for (k=1; k<=mesh->np; k++) {
      u = fabs(getSol(sol,k,is));
      sol->umax = MS_MAX(u,sol->umax);
    }
  }
  
  switch(info.nnu) {  
  /* no norm */
  case 0:
    err = CTE2D / info.eps;
    for (k=1; k<=mesh->np; k++) {
      iadr = (k-1)*3 + 1;
      for (i=0; i<3; i++)  sol->hes[iadr+i] *= err;
    }
    break;

  /* relative value: M(u)= |H(u)| / err*||u||_inf */
  case 1:
  default:
    if ( fabs(sol->umax) < 1e-30 )  return(1);
    err =  CTE2D / (info.eps * sol->umax);
    for (k=1; k<=mesh->np; k++) {
      iadr = (k-1)*3 + 1;
      for (i=0; i<3; i++)  sol->hes[iadr+i] *= err;
    }
    break;

  /* local norm: M(u)= |H(u)| / err*|u| */
  case 2:
    errs = sol->umax > 0.0 ? sol->umax*0.01 : 0.01;
    for (k=1; k<=mesh->np; k++) {
      u    = fabs(getSol(sol,k,is));
      err  = CTE2D / MS_MAX(errs,u);
      iadr = (k-1)*3 + 1;
      for (i=0; i<3; i++)  sol->hes[iadr+i] *= err;
    }
    break;

  /* local norm: M(u)= |H(u)| / err*(e1*|u|+e2*h*|du|) */
  case 3:
      puts("A CODER");
      exit(1);
/*
    maxg = 0.0;
    for (k=1; k<=mesh->np; k++) {
      u    = fabs(getSol(sol,k,is));
      norm = der[k].grd[0]*der[k].grd[0] + der[k].grd[1]*der[k].grd[1];
      norm = sqrt(norm);
      if ( norm > maxg )  maxg = norm;
    }
    if ( sol->umax < 1.0e-30 )
      errs = 0.01;
    else
      errs = sol->umax * 0.01;
    err1 = maxg * 0.01;

    for (k=1; k<=mesh->np; k++) {
      p0    = &mesh->point[k];
      u     = getSol(sol,k,is);
      norm  = der[k].grd[0]*der[k].grd[0] + der[k].grd[1]*der[k].grd[1];
      norm  = sqrt(norm);
      norm  = MS_MAX(err1,norm);
      puts("A CODER");
      exit(1);
      //norm *= p0->rins / p0->nv;
      err1  = MS_MAX(errs,fabs(u));
      
      /* variable weight *
      err1  = 0.01*err1 + 0.01*norm;
      maxu  = CTE2D / (info.eps*err1);
      der[k].hes[0] *= maxu;
      der[k].hes[1] *= maxu;
      der[k].hes[2] *= maxu;
    } */
    break;
  }

  return(1);
}


/* build metric for levelsets */
static double fsize(pMesh mesh,double u) {
  double   dd;

  if ( u < mesh->info.width )
    dd = mesh->info.hmin;
  else
    dd = mesh->info.hmin + u * (mesh->info.hmax - mesh->info.hmin);
  
  return(dd);
}


/* compute ansotropic metric for levelset */
int metrLS_3d(pMesh mesh,pSol sol) {
  double   *tmet,*grd,*hes,mat[6],mr[6],e1[3],e2[3],e3[3],dd,u,urel,rap;
  double    tail,kappa,hh,hhmin,hhmax,lambda1,lambda2,lambda3,dx2,dy2,dz2,dxy,dxz,dyz;
  int       i,j,k,l,ias,iadr;

  hhmin = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  hhmax = 1.0 / (mesh->info.hmax*mesh->info.hmax);

  for (k=1; k<=mesh->np; k++) {
    u       = fabs(getSol(sol,k,0));
    urel    = u / sol->umax;
    hh      = fsize(mesh,urel);
    lambda3 = 1.0 / (hh*hh);

    iadr = (k-1)*3 + 1;
    grd  = &sol->grd[iadr];
    iadr = (k-1)*6 + 1;
    hes  = &sol->hes[iadr];

    /* gradient norm */
    dx2 = grd[0]*grd[0];
    dxy = grd[0]*grd[1];
    dxz = grd[0]*grd[2];
    dy2 = grd[1]*grd[1];
    dyz = grd[1]*grd[2];
    dz2 = grd[2]*grd[2];
    dd  = dx2 + dy2 + dz2;

    /* curvature */
    kappa   = 0.0;
    if ( dd > EPS1 ) {
      kappa = dx2*hes[3] + dz2*hes[3] + dy2*hes[0] + dz2*hes[0] + dx2*hes[5] + dy2*hes[5] \
            - 2.0 * (dxy*hes[1] + dxz*hes[2] + dyz*hes[4]);
      kappa = fabs(kappa) / 2.0 * sqrt(dd*dd*dd);
    }

    /* vicinity of isovalue */
    if ( urel < mesh->info.width ) {
      lambda1 = kappa / mesh->info.eps; //kappa*kappa*eps; //(eps*hh);
      lambda1 = MS_MIN(lambda1,hhmin);
      lambda1 = MS_MAX(lambda1,hhmax);
      lambda2 = lambda1;
    }
    else 
      lambda1 = lambda2 = hhmax;

    /* isotropic */
    if ( mesh->info.iso ) {
      if ( lambda3 > lambda1 ) 
        tail = 1.0 / sqrt(lambda3);
      else
        tail = 1.0 / sqrt(lambda1);
      if ( !mesh->info.metric ) {
        sol->met[k] = tail;
      }
      else
        sol->met[k] = MS_MIN(sol->met[k],tail);
    }
    else {
      /* lambda_i = sizes, truncation */
      if ( mesh->info.ani > 0.0 ) {
        if ( lambda1 > lambda3 ) {
          rap = lambda1 / lambda3;
          if ( rap*rap > mesh->info.ani ) {
            lambda3 = lambda1 /(mesh->info.ani * mesh->info.ani);
            lambda3 = MS_MIN(lambda3,hhmin);
            lambda3 = MS_MAX(lambda3,hhmax);
          }
        }
        else {
          rap = lambda3 / lambda1;
          if ( rap*rap > mesh->info.ani ) {
            lambda1 = lambda3 /(mesh->info.ani * mesh->info.ani);
            lambda1 = MS_MIN(lambda1,hhmin);
            lambda1 = MS_MAX(lambda1,hhmax);
						lambda2 = lambda1;
          }
        }
      }

      /* compute local basis */
      dd = sqrt(dd);
      if ( fabs(dd) > EPS1 ) {
        memcpy(e3,grd,3*sizeof(double));
        e3[0] /= dd;
        e3[1] /= dd;
        e3[2] /= dd;
        if ( fabs(grd[0]) > EPS1 ) {
          e2[0] = -(grd[1]+grd[2]) / grd[0];
          e2[1] = e2[2] = 1.0;
        }
        else if ( fabs(grd[1]) > EPS1 ) {
          e2[1] = -(grd[0]+grd[2]) / grd[1];
          e2[0] = e2[2] = 1.0;
        }
        else if ( fabs(grd[2]) > EPS1 ) {
          e2[2] = -(grd[0]+grd[1]) / grd[2];
          e2[0] = e2[1] = 1.0;
        }
        else {
          e2[0] = 0.0; e2[1] = 1.0; e2[2] = 0.0;
        }
        dd = sqrt(e2[0]*e2[0]+e2[1]*e2[1]+e2[2]*e2[2]);
        if ( fabs(dd) > EPS1 ) { 
          e2[0] /= dd;
          e2[1] /= dd;
          e2[2] /= dd;
        }
        e1[0] = e2[1]*e3[2] - e2[2]*e3[1];
        e1[1] = e2[2]*e3[0] - e2[0]*e3[2];
        e1[2] = e2[0]*e3[1] - e2[1]*e3[0];
      }
      else {
        e1[0] = 1.0; e1[1] = 0.0; e1[2] = 0.0;
        e2[0] = 0.0; e2[1] = 1.0; e2[2] = 0.0;
        e3[0] = 0.0; e3[1] = 0.;  e3[2] = 1.0;
      }

      ias  = (k-1)*6 + 1;
      tmet = &sol->met[ias];

      /* set coeffs directly */
      if ( !mesh->info.metric ) {
        for (l=0,i=0; i<3; i++)
          for (j=i; j<3; j++)
           tmet[l++] = lambda1*e1[i]*e1[j] + lambda2*e2[i]*e2[j] + lambda3*e3[i]*e3[j];
      }
      else {
        for (l=0,i=0; i<3; i++) {
          for (j=i; j<3; j++)
            mat[l++] = lambda1*e1[i]*e1[j] + lambda2*e2[i]*e2[j] + lambda3*e3[i]*e3[j];
        }
        if ( !redsim(tmet,mat,mr) ) {
          fprintf(stdout,"  %%%% Metric problem (point: %d)\n",k);
          if ( mesh->info.ddebug ) { 
            fprintf(stdout,"  lambda: %e %e %e\n",lambda1,lambda2,lambda3);
            fprintf(stdout,"  mat %e %e %e %e %e %e\n",mat[0],mat[1],mat[2],mat[3],mat[4],mat[5]);
            fprintf(stdout,"  mat %e %e %e %e %e %e\n",tmet[0],tmet[1],tmet[2],tmet[3],tmet[4],tmet[5]);
          }
          /*memset(tmet,0,6*sizeof(double));
          tmet[0] = tmet[3] = tmet[5] = hhmax;*/
        }
        else
          for (i=0; i<6 ; i++)  tmet[i] = mr[i];
      }
    }
  }

  return(1);
}


int metrLS_2d(pMesh mesh,pSol sol) {
  double   *tmet,*hes,*grd,mat[3],mr[3],u,urel;
  double    rap,tail,kappa,hh,hhmin,hhmax,lambda1,lambda2,dd,dx2,dy2,dxy;
  int       k,ias,iadr;

  hhmin = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  hhmax = 1.0 / (mesh->info.hmax*mesh->info.hmax);

  for (k=1; k<=mesh->np; k++) {
    u       = fabs(getSol(sol,k,0));
    urel    = u / sol->umax;
    hh      = fsize(mesh,urel);
    lambda1 = 1.0 / (hh*hh);
    
    iadr = (k-1)*2 + 1;
    grd  = &sol->grd[iadr];
    iadr = (k-1)*3 + 1;
    hes  = &sol->hes[iadr];

    dx2 = grd[0]*grd[0];
    dxy = grd[0]*grd[1];
    dy2 = grd[1]*grd[1];
    dd  = dx2 + dy2;
    kappa = dx2*hes[2] - 2.0*dxy*hes[1] + dy2*hes[0];
    kappa = fabs(kappa) / 2.0 * sqrt(dd*dd*dd);

    if ( urel < mesh->info.width ) {
      lambda2 = kappa / mesh->info.eps;
      lambda2 = MS_MAX(lambda2,hhmin);
      lambda2 = MS_MIN(lambda2,hhmax);
    }
    else {
      lambda2 = hhmax;
    }
    lambda2 = MS_MAX(lambda2,hhmin);   /* <-- chgt le 16/03/09 */
    lambda2 = MS_MIN(lambda2,hhmax);
  
    /* isotropic */
    if ( mesh->info.iso ) {
      if ( lambda1 > lambda2 ) 
        tail = 1.0 / sqrt(lambda1);
      else
        tail = 1.0 / sqrt(lambda2);
      if ( !mesh->info.metric )
        sol->met[k] = tail;
      else
        sol->met[k] = MS_MIN(sol->met[k],tail);
    }
    else {
      /* lambda_i = sizes, truncation */
      if ( mesh->info.ani > 0.0 ) {
        if ( lambda1 > lambda2 ) {
          rap = lambda1 / lambda2;
          if ( rap*rap > mesh->info.ani ) {
            lambda2 = lambda1 /(mesh->info.ani * mesh->info.ani);
            lambda2 = MS_MIN(lambda2,hhmin);
            lambda2 = MS_MAX(lambda2,hhmax);
          }
        }
        else {
          rap = lambda2 / lambda1;
          if ( rap*rap > mesh->info.ani ) {
            lambda1 = lambda2 /(mesh->info.ani * mesh->info.ani);
            lambda1 = MS_MIN(lambda1,hhmin);
            lambda1 = MS_MAX(lambda1,hhmax);
          }
        }
      }
      ias  = (k-1)*3 + 1;
      tmet = &sol->met[ias];
      if ( !mesh->info.metric ) {
        if ( fabs(hes[0]) < EPS1 && fabs(hes[1]) < EPS1 && fabs(hes[2]) < EPS1 ) {
          tmet[0] = lambda1;
          tmet[1] = 0.0;
          tmet[2] = lambda2;
        }
        else {
          tmet[0] = lambda1*dx2 + lambda2*dy2;
          tmet[1] = (lambda1-lambda2)*dxy;
          tmet[2] = lambda2*dx2 + lambda1*dy2;
        }
      }
      else {
        if ( fabs(hes[0]) < EPS1 && fabs(hes[1]) < EPS1 && fabs(hes[2]) < EPS1 ) {
          mat[0] = lambda1;
          mat[1] = 0.0;
          mat[2] = lambda2;
        }
        else {
          mat[0] = lambda1*dx2 + lambda2*dy2;
          mat[1] = (lambda1-lambda2)*dxy;
          mat[2] = lambda2*dx2 + lambda1*dy2;
        }
        if ( !redsim(tmet,mat,mr) )  return(0);
        tmet[0] = mr[0];
        tmet[1] = mr[1];
        tmet[2] = mr[2];
      }
    }
  }

  return(1);
}
