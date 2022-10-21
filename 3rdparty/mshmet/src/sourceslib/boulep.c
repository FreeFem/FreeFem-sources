#include "mshmet.h"

unsigned char inxt[3] = {1,2,0};
unsigned char iprv[3] = {2,0,1};


/* find all vertices connected to P
   in:  start : tetrahedron containing p 
        ip    : index of p in start
        list  : dynamic list structure (allocated outside)
   out: list  : list of tets */
int boulep_3d(pMesh mesh,int start,int ip,int *list) {
  pTetra  pt,pt1;
  pPoint  ppt;
  int    *adja,adj,i,j,kk,indp,iel,iadr,base,ilist,iilist,nump;
  int     llist[LONMAX+2],vois[4];

  pt   = &mesh->tetra[start];
  if ( !pt->v[0] )  return(0);
  nump = pt->v[ip];
  ppt  = &mesh->point[nump];

  /* store initial tetra */
  base     = ++mesh->mark;
  pt->mark = base;
  ilist    = 0;
  iilist   = 1;
  llist[1] = 4*start + ip;
  /* store first 3 vertices */
  for (j=0; j<4; j++) {
    if ( j != ip ) {
      ilist++;
      list[ilist] = pt->v[j];
      mesh->point[ pt->v[j] ].mark = base;
    }
  }

  /* store 3 neighbors sharing P */
  iadr = (start-1)*4 + 1;
  adja = &mesh->adja[iadr];
  vois[0]  = adja[0] >> 2;
  vois[1]  = adja[1] >> 2;
  vois[2]  = adja[2] >> 2;
  vois[3]  = adja[3] >> 2;
  for (i=0; i<4; i++) {
    if ( i == ip )  continue;
    adj = vois[i];
    if ( adj ) {
      pt1 = &mesh->tetra[adj];
      if ( pt1->mark != base ) {
        /* not stored yet */
        pt1->mark = base;
        for (j=0; j<4; j++)
          if ( pt1->v[j] == nump )  break;
        iilist++;
        llist[iilist] = 4*adj + j;
      }
    }
  }
  if ( iilist < 2 )  return(ilist);

  /* explore list of neighbors */
  indp = 2;
  do {
    iel  = llist[indp] >> 2;
    pt   = &mesh->tetra[iel];
    iadr = (iel-1)*4 + 1;
    adja = &mesh->adja[iadr];
    vois[0]  = adja[0] >> 2;
    vois[1]  = adja[1] >> 2;
    vois[2]  = adja[2] >> 2;
    vois[3]  = adja[3] >> 2;

    for (i=0; i<4; i++) {
      if ( pt->v[i] == nump )  continue;
      adj = vois[i];
      if ( adj ) {
        pt1 = &mesh->tetra[adj];
        if ( pt1->mark != base ) {
          pt1->mark = base;
          for (j=0; j<4; j++)
            if ( pt1->v[j] == nump )  break;
          iilist++;
          llist[iilist] = 4*adj + j;
        }
      }
    }
    /* overflow */
    if ( iilist > LONMAX-3 )  return(-ilist);
  }
  while ( ++indp <= iilist );

  /* store vertices from tetra list */
  for (i=2; i<=iilist; i++) {
    kk = llist[i] / 4;
    pt = &mesh->tetra[kk];
    for (j=0; j<4; j++) {
      if ( pt->v[j] != nump ) {
        ppt = &mesh->point[ pt->v[j] ];
        if ( ppt->mark < base ) {
          ilist++;
          list[ilist] = pt->v[j];
          ppt->mark = base;
        }
      }
    }
  }

  return(ilist);
}


/* store neighboring vertices, return < 0 if overflow */
int boulep_2d(pMesh mesh,int start,int ip,int *list) {
  pTria    pt;
  int     *adja,i1,iadr,nump,voy,ilist,iel;

  if ( start < 1 )  return(0);
  pt = &mesh->tria[start];
  if ( !pt->v[0] )  return(0);
  nump = pt->v[ip];

  /* init list */
  i1      = inxt[ip];
  ilist   = 1;
  list[1] = pt->v[i1];

  iadr = (start-1)*3 + 1;
  adja = &mesh->adja[iadr];
  iel  =  adja[i1] / 3;
  while ( iel && (iel != start) ) {
    pt  = &mesh->tria[iel]; 
    voy = adja[i1] % 3;
    i1  = iprv[voy];
    ++ilist;
    if ( ilist > LONMAX-2 )  return(-ilist);
    list[ilist] = pt->v[i1];
    iadr = (iel-1)*3 + 1;
    adja = &mesh->adja[iadr];
    iel  = adja[i1] / 3;
  }

  /* reverse loop */
  if ( iel != start ) {
    voy  = inxt[i1];
    ++ilist;
    list[ilist] = pt->v[voy];

    pt   = & mesh->tria[start];
    i1   = iprv[ip];
    iadr = (start-1)*3 + 1;
    adja = &mesh->adja[iadr];
    iel =  adja[i1] / 3;
    while ( iel && (iel != start) ) {
      pt  = &mesh->tria[iel];
      voy = adja[i1] % 3;
      ++ilist;
      if ( ilist > LONMAX-2 )  return(-ilist);
      list[ilist] = pt->v[voy];
      iadr = (iel-1)*3 + 1;
      adja = &mesh->adja[iadr];
      i1   = inxt[voy];
      iel  = adja[i1] / 3;
    }
  }

  return(ilist);
}



