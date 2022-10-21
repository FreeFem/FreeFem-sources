#include "mshmet.h"

#define SQRT3DIV2   0.8660254037844386


/* control mesh gradation (aniso) */
static int hgrada_3d(pMesh mesh,pSol sol) {
  pPoint      p1,p2;
  pHash       hash;
  hedge      *pht;
  double      ux,uy,uz,dd,d1,d2,dh,dd1,dd2,ma1[6],mb1[6];
  double      logh,tail,rap,logs,coef,coef1,coef2;
  double     *ma,*mb;
  int         i,k,a,b,nc,nbc,it,maxt;

  /* reset color */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].mark = mesh->mark;

  /* build edge list */
  hash = hashEdge_3d(mesh);
  assert(hash);

  /* default values */
  logh = log(mesh->info.hgrad);
  logs = 0.001 + logh;
  it   = 0;
  nbc  = 0;
  maxt = 100;

  /* scan edge table */
  do {
		mesh->mark++;
    nc = 0;
    for (k=1; k<=hash->hnxt; k++) {
      pht = &hash->item[k];
      if ( !pht->min )  continue;
      a  = pht->min;
      b  = pht->max;
      p1 = &mesh->point[a];
      p2 = &mesh->point[b];
			if ( p1->mark < mesh->mark-1 && p1->mark < mesh->mark-1 )  continue;

      ma = &sol->met[6*(a-1)+1];
      mb = &sol->met[6*(b-1)+1];
      
      /* compute edge lengths */
      ux = p2->c[0] - p1->c[0];
      uy = p2->c[1] - p1->c[1];
      uz = p2->c[2] - p1->c[2];
      d1 =      ma[0]*ux*ux + ma[3]*uy*uy + ma[5]*uz*uz \
         + 2.0*(ma[1]*ux*uy + ma[2]*ux*uz + ma[4]*uy*uz);
      if ( d1 <= EPS1 )  d1 = EPS1;
      dd1 = MS_MAX(EPS1,sqrt(d1));
      d2 =      mb[0]*ux*ux + mb[3]*uy*uy + mb[5]*uz*uz \
         + 2.0*(mb[1]*ux*uy + mb[2]*ux*uz + mb[4]*uy*uz);
      if ( d2 <= EPS1 )  d2 = EPS1;
      dd2 = MS_MAX(EPS1,sqrt(d2));
      
      /* swap vertices */
      if ( dd1 > dd2 ) {
        p1  = &mesh->point[b];
        p2  = &mesh->point[a];
        ma = &sol->met[6*(b-1)+1];
        mb = &sol->met[6*(a-1)+1];
        dd  = dd1;
        dd1 = dd2;
        dd2 = dd;
      }
      
      /* check edge ratio */
      rap = dd2 / dd1;
      dh = rap - 1.0;
      if ( fabs(dh) > EPS1 ) {
        tail = (dd1+dd2+4*sqrt(0.5*(d1+d2))) / 6.0;
        coef = log(rap) / tail;
      
        /* update sizes */
        if ( coef > logs ) {
          coef      = exp(tail*logh);
          p1->mark = mesh->mark;
          p2->mark = mesh->mark;
      
          coef1 = 1.0 + logh*dd2;
          coef1 = 1.0 / (coef1*coef1);
          coef2 = 1.0 + logh*dd1;
          coef2 = 1.0 / (coef2*coef2);
      
          /* metric intersection */
          for (i=0; i<6; i++) {
            ma1[i] = coef2 * ma[i];
            mb1[i] = coef1 * mb[i];
          }
          if ( !redsim_3d(ma,mb1,ma) )
            for (i=0; i<6; i++)  ma[i] = SQRT3DIV2 * (ma[i]+mb1[i]);
          if ( !redsim_3d(ma1,mb,mb) )
            for (i=0; i<6; i++)  mb[i] = SQRT3DIV2 * (mb[i]+ma1[i]);
          nc++;
        }
      }
    }
    nbc += nc;
    if ( nc > 0 )  fprintf(stdout,"  -- %d sizes adapted\n",nc);
  } 
  while ( nc > 0 && ++it <= maxt );

  return(1);
}


/* control mesh gradation (iso) */
static int hgradi_3d(pMesh mesh,pSol sol) {
  pPoint      p1,p2;
  pHash       hash;
  hedge      *pht;
  double      ux,uy,uz,dd,dh;
  double      logh,tail,rap,logs,lograp,coef;
  double     *ma,*mb;
  int         k,a,b,nc,nbc,it,maxt;

  /* reset color */
	mesh->mark = 0;
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].mark = mesh->mark;

  /* build edge list */
  hash = hashEdge_3d(mesh);
  assert(hash);

  /* default values */
  logh = log(mesh->info.hgrad);
  logs = 0.001 + logh;
  it   = 0;
  nbc  = 0;
  maxt = 100;

  /* scan edge table */
  do {
    mesh->mark++;
    nc = 0;
    for (k=1; k<=hash->hnxt; k++) {
      pht = &hash->item[k];
      if ( !pht->min )  continue;
      a  = pht->min;
      b  = pht->max;
      p1 = &mesh->point[a];
      p2 = &mesh->point[b];
			if ( p1->mark < mesh->mark-1 && p1->mark < mesh->mark-1 )  continue;

      /* compute edge lengths */
      ux = p2->c[0] - p1->c[0];
      uy = p2->c[1] - p1->c[1];
      uz = p2->c[2] - p1->c[2];
      dd = sqrt(ux*ux + uy*uy + uz*uz);
      ma = &sol->met[a];
      mb = &sol->met[b];
      if ( ma[0] > mb[0] ) {
        p1 = &mesh->point[b];
        p2 = &mesh->point[a];
        ma = &sol->met[b];
        mb = &sol->met[a];
	    }
      rap = mb[0] / ma[0];
      dh  = rap - 1.0;
      if ( fabs(dh) > EPS ) {
	      lograp = log(rap);
	      tail   = dd * dh / (mb[0]*lograp);
        coef   = lograp / tail;
	      /* update sizes */
	      if ( coef > logs ) {
          p2->mark = mesh->mark;
	        mb[0] = ma[0] * exp(tail*logh);
	        nc++;
	      }
      }
    }
    nbc += nc;
    if ( nc > 0 )  fprintf(stdout,"  -- %d sizes adapted\n",nc);
  } 
  while ( nc > 0 && ++it <= maxt );

  return(1);
}


/* control mesh gradation (aniso) */
static int hgrada_2d(pMesh mesh,pSol sol) {
  pPoint      p1,p2;
  pHash       hash;
  hedge      *pht;
  double      ux,uy,dd,d1,d2,dh,dd1,dd2,ma1[6],mb1[6];
  double      logh,tail,rap,logs,coef,coef1,coef2;
  double     *ma,*mb;
  int         i,k,a,b,nc,nbc,it,maxt;

  /* reset color */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].mark = mesh->mark;

  /* build edge list */
  hash = hashEdge_2d(mesh);
  assert(hash);

  /* default values */
  logh = log(mesh->info.hgrad);
  logs = 0.001 + logh;
  it   = 0;
  nbc  = 0;
  maxt = 100;

  /* scan edge table */
  do {
    mesh->mark++;
    nc = 0;
    for (k=1; k<=hash->hnxt; k++) {
      pht = &hash->item[k];
      if ( !pht->min )  continue;
      a  = pht->min;
      b  = pht->max;
      p1 = &mesh->point[a];
      p2 = &mesh->point[b];
			if ( p1->mark < mesh->mark-1 && p1->mark < mesh->mark-1 )  continue;

      ma = &sol->met[3*(a-1)+1];
      mb = &sol->met[3*(b-1)+1];
      
      /* compute edge lengths */
      ux = p2->c[0] - p1->c[0];
      uy = p2->c[1] - p1->c[1];
      d1 = ma[0]*ux*ux + ma[2]*uy*uy + 2.0*ma[1]*ux*uy;
      if ( d1 <= EPS1 )  d1 = EPS1;
      dd1 = MS_MAX(EPS1,sqrt(d1));
      d2 = mb[0]*ux*ux + mb[2]*uy*uy + 2.0*mb[1]*ux*uy;
      if ( d2 <= EPS1 )  d2 = EPS1;
      dd2 = MS_MAX(EPS1,sqrt(d2));
      
      /* swap vertices */
      if ( dd1 > dd2 ) {
        p1  = &mesh->point[b];
        p2  = &mesh->point[a];
        ma = &sol->met[3*(b-1)+1];
        mb = &sol->met[3*(a-1)+1];
        dd  = dd1;
        dd1 = dd2;
        dd2 = dd;
      }
      
      /* check edge ratio */
      rap = dd2 / dd1;
      dh = rap - 1.0;
      if ( fabs(dh) > EPS1 ) {
        tail = (dd1+dd2+4*sqrt(0.5*(d1+d2))) / 6.0;
        coef = log(rap) / tail;
      
        /* update sizes */
        if ( coef > logs ) {
          coef      = exp(tail*logh);
          p1->mark = mesh->mark;
          p2->mark = mesh->mark;
      
          coef1 = 1.0 + logh*dd2;
          coef1 = 1.0 / (coef1*coef1);
          coef2 = 1.0 + logh*dd1;
          coef2 = 1.0 / (coef2*coef2);
      
          /* metric intersection */
          for (i=0; i<6; i++) {
            ma1[i] = coef2 * ma[i];
            mb1[i] = coef1 * mb[i];
          }
          if ( !redsim_2d(ma,mb1,ma) )
            for (i=0; i<3; i++)  ma[i] = SQRT3DIV2 * (ma[i]+mb1[i]);
          if ( !redsim_2d(ma1,mb,mb) )
            for (i=0; i<3; i++)  mb[i] = SQRT3DIV2 * (mb[i]+ma1[i]);
          nc++;
        }
      }
      
    }
    nbc += nc;
    if ( nc > 0 )  fprintf(stdout,"  -- %d sizes adapted\n",nc);
  } 
  while ( nc > 0 && ++it <= maxt );

  return(1);
}


/* control mesh gradation (aniso) */
static int hgradi_2d(pMesh mesh,pSol sol) {
  pPoint      p1,p2;
  pHash       hash;
  hedge      *pht;
  double      ux,uy,dd,dh;
  double      logh,tail,rap,logs,lograp,coef;
  double     *ma,*mb;
  int         k,a,b,nc,nbc,it,maxt;

  /* reset color */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].mark = mesh->mark;

  /* build edge list */
  hash = hashEdge_2d(mesh);
  assert(hash);

  /* default values */
  logh = log(mesh->info.hgrad);
  logs = 0.001 + logh;
  it   = 0;
  nbc  = 0;
  maxt = 100;

  /* scan edge table */
  do {
    mesh->mark++;
    nc = 0;
    for (k=1; k<=hash->hnxt; k++) {
      pht = &hash->item[k];
      if ( !pht->min )  continue;
      a  = pht->min;
      b  = pht->max;
      p1 = &mesh->point[a];
      p2 = &mesh->point[b];
			if ( p1->mark < mesh->mark-1 && p1->mark < mesh->mark-1 )  continue;

      /* compute edge lengths */
      ux = p2->c[0] - p1->c[0];
      uy = p2->c[1] - p1->c[1];
      dd = sqrt(ux*ux + uy*uy);
      ma = &sol->met[a];
      mb = &sol->met[b];
      if ( ma[0] > mb[0] ) {
        ma = &sol->met[b];
        mb = &sol->met[a];
        p1 = &mesh->point[b];
        p2 = &mesh->point[a];
	    }
      rap = mb[0] / ma[0];
      dh  = rap - 1.0;
      if ( fabs(dh) > EPS ) {
	      lograp = log(rap);
	      tail   = dd * dh / (mb[0]*lograp);
        coef   = lograp / tail;
	      /* update sizes */
	      if ( coef > logs ) {
          p1->mark = mesh->mark;
          p2->mark = mesh->mark;
	        mb[0]  = ma[0] * exp(tail*logh);
	        nc++;
	      }
      }
      
    }
    nbc += nc;
    if ( nc > 0 )  fprintf(stdout,"  -- %d sizes adapted\n",nc);
  } 
  while ( nc > 0 && ++it <= maxt );

  return(1);
}


int mshme1(pMesh mesh,pSol sol) {
  pPoint   ppt;
  int      j,k,siz,ier,it,nex,ney;

  /* memory alloc */
  sol->grd = (double*)M_calloc(sol->np+1,sol->dim*sizeof(double),"grd");
  assert(sol->grd);
  siz = sol->dim == 2 ? 3 : 6;
  sol->hes = (double*)M_calloc(sol->np+1,siz*sizeof(double),"hes");
  assert(sol->hes);

  if ( mesh->info.metric == 0 ) {
    if ( mesh->info.iso )
      sol->met = (double*)M_calloc(sol->np+1,sizeof(double),"mshme1");
    else {
      siz = sol->dim == 2 ? 3 : 6;
      sol->met = (double*)M_calloc(sol->np+1,siz*sizeof(double),"mshme1");
    }
    assert(sol->met);
  }
  if ( mesh->info.ls )  sol->type = 1;

  for (j=0; j<sol->type; j++) {
    if ( mesh->info.nsol > -1 && mesh->info.nsol != j )  continue;
    if ( mesh->info.imprim < 0 )
      fprintf(stdout,"     Solution %d: %d vertices\n",j+1,mesh->np);

    /* solution smoothing */
    if ( mesh->info.nlis )  lissag(mesh,sol,j,mesh->info.nlis);

    /* compute gradient */
    nex = ney = 0;
    for (k=1; k<=mesh->np; k++) {
			ier = gradLS(mesh,sol,k,j);
			if ( !ier ) {
        fprintf(stdout,"  %%%% Unable to evaluate gradient %d\n.",k);
        return(0);
      }
    }

    /* compute hessian */
    for (k=1; k<=mesh->np; k++) {
      ier = hessLS(mesh,sol,k,j);
      if ( !ier )  return(0);
      else if ( ier < 0 )  nex++;
    }

    /* correction */
    if ( nex ) {
      ney = nex;
      it  = 0;
      do {
        nex = 0;
        for (k=1; k<=mesh->np; k++) {
          ppt = &mesh->point[k];
          if ( !ppt->h ) {
            ier = avgval(mesh,sol,k);
            if ( !ier )  ier = clsval(mesh,sol,k);
            if ( !ier )  nex++;
          }  
        }
      }
      while ( nex > 0 && ++it < 100 );
    }

    if ( mesh->info.imprim < 0 && ney )
      fprintf(stdout,"  %%%% %d corrected  %d unknowns\n",ney,nex);

    /* norm + avg hessian */
    if ( !nrmhes(mesh,sol,j) )  return(0);
    //if ( !laplac(mesh,der) )    return(0);

    /* define metric */
    if ( mesh->info.ddebug )  outder(mesh,sol);
    if ( mesh->info.ls )
      break;
    else if ( !defmet(mesh,sol,j) )  return(0);
  }

  /* levelset metric */
  if ( mesh->info.ls )
    if ( !metrLS(mesh,sol) )  return(0);

  M_free(sol->grd);
  M_free(sol->hes);
  if ( sol->dim == 3 && !mesh->ne )  M_free(sol->nn);

  if ( mesh->info.hgrad > 0.0 ) {
    if ( mesh->dim == 2 )
      ier = mesh->info.iso ? hgradi_2d(mesh,sol) : hgrada_2d(mesh,sol);
    else
      ier = mesh->info.iso ? hgradi_3d(mesh,sol) : hgrada_3d(mesh,sol);
  }
  
  return(1);
}
