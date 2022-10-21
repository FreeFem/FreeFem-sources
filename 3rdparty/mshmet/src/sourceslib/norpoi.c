#include "mshmet.h"


int norpoi(pMesh mesh,pSol sol) {
  pTria     pt;
	pPoint    p0,p1,p2;
	double   *np,ux,uy,uz,vx,vy,vz,d1,d2,d,nx,ny,nz,poi,sum;
  int       i,k,iadr,nsdep,ip1,ip2,list[LONMAX+2],ilist;

  assert(mesh->nt);
	sol->nn = (double*)M_calloc(sol->np,sol->dim*sizeof(double),"nn");
	assert(sol->nn);

  for (k=1; k<=mesh->np; k++) {
    p0    = &mesh->point[k];
		nsdep = p0->s;

		pt = &mesh->tria[nsdep];
    if ( pt->v[0] == k )      i = 0;
    else if ( pt->v[1] == k ) i = 1;
    else                      i = 2;

		ilist = boulep(mesh,nsdep,i,list);
		if ( ilist < 1 )  return(0);

 		iadr = (k-1)*sol->dim + 1;
 		np   = &sol->nn[iadr];
		sum  = 0.0;

		ip1 = list[1];
		p1  = &mesh->point[ip1];
		for (i=2; i<=ilist; i++) {
      ip2 = list[i];
      p2  = &mesh->point[ip2];
 
			ux = p1->c[0] - p0->c[0];
			uy = p1->c[1] - p0->c[1];
			uz = p1->c[2] - p0->c[2];
			d1 = ux*ux + uy*uy + uz*uz;

			vx = p2->c[0] - p0->c[0];
			vy = p2->c[1] - p0->c[1];
			vz = p2->c[2] - p0->c[2];
			d2 = vx*vx + vy*vy + vz*vz;
			
			nx = uy*vz - uz*vy;
      ny = uz*vx - ux*vz;
      nz = ux*vy - uy*vx;
      d  = nx*nx + ny*ny + nz*nz;
      if ( d > 0.0 ) {
				d = 1.0 / sqrt(d);
	      poi = (ux*vx + uy*vy + uz*vz) / sqrt(d1*d2);
        poi = acos(poi) * d;
        np[0] += poi * nx;
        np[1] += poi * ny;
        np[2] += poi * nz;
        sum   += poi;
      }
			p1 = p2;
    }
    
		d = np[0]*np[0] + np[1]*np[1] + np[2]*np[2];
		if ( d > 0.0 ) {
			d = 1.0 / sqrt(d);
			np[0] *= d;
			np[1] *= d;
			np[2] *= d;
		}
  }

  return(1);
}
