#include "mshmet.h"

#define tuMv(u,m,v) (\
  (u[0]*m[0]+u[1]*m[1])*v[0] + \
  (u[0]*m[1]+u[1]*m[2])*v[1])


/* invert 3x3 symmetric matrix */
static int invmat(double m[6],double mi[6]) {
  double  aa,bb,cc,det,vmin,vmax,maxx;
  int     k;

  /* check diagonal matrices */
  vmax = fabs(m[1]);
  maxx = fabs(m[2]);
  if( maxx > vmax ) vmax = maxx;
  maxx = fabs(m[4]);
  if( maxx > vmax ) vmax = maxx;
  if ( vmax < EPS ) {
    mi[0]  = 1.0/m[0];
    mi[3]  = 1.0/m[3];
    mi[5]  = 1.0/m[5];
    mi[1] = mi[2] = mi[4] = 0.0;
    return(1);
  }

  /* check ill-conditionned matrix */
  vmin = vmax = fabs(m[0]);
  for (k=1; k<6; k++) {
    maxx = fabs(m[k]);
    if ( maxx < vmin )  vmin = maxx;
    else if ( maxx > vmax )  vmax = maxx;
  }
	if ( vmax == 0.0 )  return(0);

  /* compute sub-dets */
  aa  = m[3]*m[5] - m[4]*m[4];
  bb  = m[4]*m[2] - m[1]*m[5];
  cc  = m[1]*m[4] - m[2]*m[3];
  det = m[0]*aa + m[1]*bb + m[2]*cc;
	if ( fabs(det) < EPS1 )  return(0);
  det = 1.0 / det;

  mi[0] = aa*det;
  mi[1] = bb*det;
  mi[2] = cc*det;
  mi[3] = (m[0]*m[5] - m[2]*m[2])*det;
  mi[4] = (m[1]*m[2] - m[0]*m[4])*det;
  mi[5] = (m[0]*m[3] - m[1]*m[1])*det;

  return(1);
}


/* invert 3x3 non-symmetric matrix */
int invmatg(double m[9],double mi[9]) {
  double  aa,bb,cc,det,vmin,vmax,maxx;
  int     k;

  /* check ill-conditionned matrix */
  vmin = vmax = fabs(m[0]);
  for (k=1; k<9; k++) {
    maxx = fabs(m[k]);
    if ( maxx < vmin )  vmin = maxx;
    else if ( maxx > vmax )  vmax = maxx;
  }
  if ( vmax == 0.0 )  return(0);

  /* compute sub-dets */
  aa = m[4]*m[8] - m[5]*m[7];
  bb = m[5]*m[6] - m[3]*m[8];
  cc = m[3]*m[7] - m[4]*m[6];
  det = m[0]*aa + m[1]*bb + m[2]*cc;
  if ( fabs(det) < EPS1 )  return(0);
  det = 1.0f / det;

  mi[0] = aa*det;
  mi[3] = bb*det;
  mi[6] = cc*det;
  mi[1] = (m[2]*m[7] - m[1]*m[8])*det;
  mi[4] = (m[0]*m[8] - m[2]*m[6])*det;
  mi[7] = (m[1]*m[6] - m[0]*m[7])*det;
  mi[2] = (m[1]*m[5] - m[2]*m[4])*det;
  mi[5] = (m[2]*m[3] - m[0]*m[5])*det;
  mi[8] = (m[0]*m[4] - m[1]*m[3])*det;

  return(1);
}


/* simultaneous reduction  mr = m1 cap m2 */
int redsim_3d(double *m1,double *m2,double *mr) {
  double  dd,lambda[3],v[3][3],hh[3];
  double  maxd1,maxd2,ex,ey,ez,m1i[6],n[9],p[9],pi[9];
  int     i,j,k,order;

  maxd1 = fabs(m1[1]);
  maxd1 = MS_MAX(maxd1,fabs(m1[2]));
  maxd1 = MS_MAX(maxd1,fabs(m1[4]));
  
  maxd2 = fabs(m2[1]);
  maxd2 = MS_MAX(maxd2,fabs(m2[2]));
  maxd2 = MS_MAX(maxd2,fabs(m2[4]));

  /* check diag matrices */
  if ( maxd1 < EPS1 && maxd2 < EPS1 ) {
    mr[0] = MS_MAX(m1[0],m2[0]);
    mr[3] = MS_MAX(m1[3],m2[3]);
    mr[5] = MS_MAX(m1[5],m2[5]);
    mr[1] = mr[2] = mr[4] = 0.0;
    return(1);
  }
	if ( !invmat(m1,m1i) )  return(0);

  /* n = (m1)^-1*m2 */
  n[0] = m1i[0]*m2[0] + m1i[1]*m2[1] + m1i[2]*m2[2];
  n[1] = m1i[0]*m2[1] + m1i[1]*m2[3] + m1i[2]*m2[4];
  n[2] = m1i[0]*m2[2] + m1i[1]*m2[4] + m1i[2]*m2[5];
  n[3] = m1i[1]*m2[0] + m1i[3]*m2[1] + m1i[4]*m2[2];
  n[4] = m1i[1]*m2[1] + m1i[3]*m2[3] + m1i[4]*m2[4];
  n[5] = m1i[1]*m2[2] + m1i[3]*m2[4] + m1i[4]*m2[5];
  n[6] = m1i[2]*m2[0] + m1i[4]*m2[1] + m1i[5]*m2[2];
  n[7] = m1i[2]*m2[1] + m1i[4]*m2[3] + m1i[5]*m2[4];
  n[8] = m1i[2]*m2[2] + m1i[4]*m2[4] + m1i[5]*m2[5];

  /* eigenvectors */
  order = eigenv(0,n,lambda,v);

  if ( order == 0 ) 
    return(0);
  else if ( order == 3 ) {
    mr[0] = mr[3] = mr[5] = lambda[0];
    mr[1] = mr[2] = mr[4] = 0.0;
    return(1);
  }
  else {
    /* matrix of passage */
    for (i=0,k=0; i<3; i++)
      for (j=0; j<3; j++)
        p[k++] = v[j][i];
    if ( !invmatg(p,pi) )  return(0);

    for (i=0; i<3; i++) {
      ex = v[i][0];
      ey = v[i][1];
      ez = v[i][2];
      maxd1 = ex*(m1[0]*ex+m1[1]*ey+m1[2]*ez) 
            + ey*(m1[1]*ex+m1[3]*ey+m1[4]*ez)
            + ez*(m1[2]*ex+m1[4]*ey+m1[5]*ez);
      maxd2 = ex*(m2[0]*ex+m2[1]*ey+m2[2]*ez) 
            + ey*(m2[1]*ex+m2[3]*ey+m2[4]*ez)
            + ez*(m2[2]*ex+m2[4]*ey+m2[5]*ez);
      hh[i] = MS_MAX(maxd1,maxd2);
    }

    /* compose matrix tP^-1*lambda*P^-1 */
    mr[0] = pi[0]*hh[0]*pi[0] + pi[3]*hh[1]*pi[3] + pi[6]*hh[2]*pi[6];
    mr[1] = pi[0]*hh[0]*pi[1] + pi[3]*hh[1]*pi[4] + pi[6]*hh[2]*pi[7];
    mr[2] = pi[0]*hh[0]*pi[2] + pi[3]*hh[1]*pi[5] + pi[6]*hh[2]*pi[8];
    mr[3] = pi[1]*hh[0]*pi[1] + pi[4]*hh[1]*pi[4] + pi[7]*hh[2]*pi[7];
    mr[4] = pi[1]*hh[0]*pi[2] + pi[4]*hh[1]*pi[5] + pi[7]*hh[2]*pi[8];
    mr[5] = pi[2]*hh[0]*pi[2] + pi[5]*hh[1]*pi[5] + pi[8]*hh[2]*pi[8];

    /* for safety... */
    maxd1 = fabs(mr[0]);
    for (i=1; i<6; i++)
      maxd1 = MS_MAX(maxd1,fabs(mr[i]));
    if ( maxd1 < EPS1 ) {
      fprintf(stderr,"  %%%% Null matrix.\n");
      exit(1);
    }
    dd = 1.0 / maxd1;
    for (i=0; i<6; i++) {
      mr[i] *= dd;
      /*if ( fabs(m[i]) < EPS )  m[i] = 0.0;*/
    }

    maxd2 = fabs(mr[1]);
    maxd2 = MS_MAX(maxd2,fabs(mr[2]));
    maxd2 = MS_MAX(maxd2,fabs(mr[4]));

    if ( maxd2 < EPS ) {
      mr[1] = mr[2] = mr[4] = 0.0;
      mr[0] = fabs(mr[0]) * maxd1;
      mr[3] = fabs(mr[3]) * maxd1;
      mr[5] = fabs(mr[5]) * maxd1;
    }
    else {
      order = eigenv(1,mr,lambda,v);
      for (i=0; i<6; i++)   mr[i] *= maxd1;

      if ( lambda[0] < -EPS1 || lambda[1] < -EPS1 || lambda[2] < -EPS1 ) {
        fprintf(stderr,"  %%%% Metric inconsistency\n");
        fprintf(stderr,"  %.6f %.6f %.6f\n%.6f %.6f\n%.6f\n",
                       mr[0],mr[1],mr[2],mr[3],mr[4],mr[5]);
        fprintf(stderr,"  Lambda %f %f %f\n",lambda[0],lambda[1],lambda[2]);
	      return(0);
      }
    }
  }
  return(1);
}


/* find non orthogonal basis (v1,v2) where m1,m2 are diag */
int redsim_2d(double *m1,double *m2,double *mr) {
  double     n[4],v1[2],v2[2],lambda[2],vnorm,dd,det1,inv1,tracn;
  double     lambda1[2],lambda2[2],ip[4],tip[4],disc,det,vp1,vp2;

  /* eigenvect of n = m1^(-1) m2 */
  det1 = m1[0]*m1[2]-m1[1]*m1[1];
  if ( fabs(det1) < EPS1 ) {
	  return(0);
  }
  inv1 = 1.0 / det1;

  /* verif idem vp tq : M1v = vp*M2v <=> (M1-vp*M2)v = 0 */ 
  n[0] = (m1[2]*m2[0] - m1[1]*m2[1]) * inv1;
  n[1] = (m1[2]*m2[1] - m1[1]*m2[2]) * inv1;
  n[2] = (m1[0]*m2[1] - m1[1]*m2[0]) * inv1;
  n[3] = (m1[0]*m2[2] - m1[1]*m2[1]) * inv1;

  /* delta >=0 as N diag in IR */
  dd    = n[0] - n[3];
  disc  = sqrt(fabs(dd*dd + 4.0*n[1]*n[2]));
  tracn = n[0] + n[3];

  lambda[0] = 0.5 * (tracn + disc);
  if ( fabs(lambda[0]) < EPS1 ) {
    mr[0] = m1[0];
  	mr[1] = m1[1];
  	mr[2] = m1[2];
		return(1);
  }

  /* check for ellipse included in other */
  if ( disc / lambda[0] < EPS ) {
    /* delta<<Tr => delta~0 : eigenvalue double, P not invertible */
    if ( lambda[0] <= 1.0 ) {
  	   mr[0] = m1[0];
  	   mr[1] = m1[1];
  	   mr[2] = m1[2];
    }
    else {
  	  mr[0] = m2[0];
  	  mr[1] = m2[1];
  	  mr[2] = m2[2];
    }
    return(1);
  }
  else {
    lambda[1] = 0.5 * (tracn - disc);

    v1[0] = -n[1];
    v1[1] =  n[0] - lambda[0];
    vnorm = sqrt(v1[0]*v1[0] + v1[1]*v1[1]);
    if ( vnorm < EPS1 ) {
       v1[0] = lambda[0] - n[3];
       v1[1] = n[2];   
       vnorm = sqrt(v1[0]*v1[0] + v1[1]*v1[1]);         
    }
    vnorm = 1.0 / vnorm;
    v1[0] *= vnorm;
    v1[1] *= vnorm;

    v2[0] = -n[1];
    v2[1] =  n[0] - lambda[1];
    vnorm = sqrt(v2[0]*v2[0] + v2[1]*v2[1]);
    if ( vnorm < EPS1 ) {
       v2[0] = lambda[1] - n[3];
       v2[1] = n[2];   
       vnorm = sqrt(v2[0]*v2[0] + v2[1]*v2[1]);         
    }
    vnorm = 1.0 / vnorm;
    v2[0] *= vnorm;
    v2[1] *= vnorm;

    lambda1[0] = tuMv(v1,m1,v1);
    lambda2[0] = tuMv(v1,m2,v1);
    lambda1[1] = tuMv(v2,m1,v2);
    lambda2[1] = tuMv(v2,m2,v2);

    /* intersect (m1,m2)= tP^(-1)&*L&*P^(-1) */
    det = v1[0]*v2[1] - v2[0]*v1[1];
		if ( fabs(det) < EPS1 )  return(0);
    det = 1.0 / det;

    ip[0] =  v2[1]*det;
    ip[1] = -v2[0]*det;
    ip[2] = -v1[1]*det;
    ip[3] =  v1[0]*det;

    /* tP^(-1)&*L */
    vp1 = MS_MAX(lambda1[0],lambda2[0]);
    vp2 = MS_MAX(lambda1[1],lambda2[1]);

    tip[0] = ip[0] * vp1;
    tip[1] = ip[2] * vp2;
    tip[2] = ip[1] * vp1;
    tip[3] = ip[3] * vp2;

    mr[0] = tip[0]*ip[0] + tip[1]*ip[2];
    mr[1] = tip[0]*ip[1] + tip[1]*ip[3];
    mr[2] = tip[2]*ip[1] + tip[3]*ip[3]; 
  }

  return(1);
}
