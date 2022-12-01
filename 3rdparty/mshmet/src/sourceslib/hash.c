#include "mshmet.h"

#define KA     31
#define KB     57
#define KC     79

#define KTA     7
#define KTB    11

static int idirt[7] = {0,1,2,3,0,1,2};


/* hash mesh edges */
pHash hashEdge_3d(pMesh mesh) {
  pTetra   pt;
  hedge   *ph;
	pHash    hash;
  int      i,j,k,ia,ib,mins,maxs,key;

  /* adjust hash table params */
	hash = malloc(1*sizeof(Hash));
	assert(hash);
  hash->item = (hedge*)calloc(9*mesh->np,sizeof(hedge));
  assert(hash->item);
  hash->size  = mesh->np;
  hash->nhmax = 9*mesh->np - 1;
  hash->hnxt  = mesh->np;
  memset(hash->item,0,hash->nhmax*sizeof(hedge));
  for (k=mesh->np; k<hash->nhmax; k++)
    hash->item[k].nxt = k+1;

  /* scan tetra edges */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    for (i=0; i<4; i++) {
      ia = pt->v[i];
      for (j=i+1; j<4; j++) {
        ib = pt->v[j];
        mins = ia;
        maxs = ib;
        if ( ia > ib ) {
          mins = ib;
          maxs = ia;
        }
        key = (KA*mins + KB*maxs) % hash->size;
        ph  = &hash->item[key];

        if ( !ph->min ) {
          ph->min = mins;
          ph->max = maxs;
          ph->nxt = 0;
	      }
	      else {
          while ( ph->nxt && ph->nxt < hash->nhmax ) {
						if ( ph->min == mins && ph->max == maxs )  break;
            ph = &hash->item[ph->nxt];
          }
					if ( ph->min == mins && ph->max == maxs )  continue;
					ph->nxt = hash->hnxt;
          ph      = &hash->item[hash->hnxt];
          ph->min = mins;
          ph->max = maxs;
          ph->nxt = 0;
          ++hash->hnxt;
					assert(hash->hnxt < hash->nhmax );
        }
      }
    }
  }      

  return(hash);
}


/* hash mesh edges */
pHash hashEdge_2d(pMesh mesh) {
  pTria    pt;
  hedge   *ph;
	pHash    hash;
  int      i,j,k,ia,ib,mins,maxs,key;

  /* adjust hash table params */
	hash = malloc(1*sizeof(Hash));
	assert(hash);
  hash->item = (hedge*)calloc(4*mesh->np,sizeof(hedge));
  assert(hash->item);
  hash->size  = mesh->np;
  hash->nhmax = 4*mesh->np - 1;
  hash->hnxt  = mesh->np;
  memset(hash->item,0,hash->nhmax*sizeof(hedge));
  for (k=mesh->np; k<hash->nhmax; k++)
    hash->item[k].nxt = k+1;

  /* scan tetra edges */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    for (i=0; i<3; i++) {
      ia = pt->v[i];
      for (j=i+1; j<3; j++) {
        ib = pt->v[j];
        mins = ia;
        maxs = ib;
        if ( ia > ib ) {
          mins = ib;
          maxs = ia;
        }
        key = (KA*mins + KB*maxs) % hash->size;
        ph  = &hash->item[key];

        if ( !ph->min ) {
          ph->min = mins;
          ph->max = maxs;
          ph->nxt = 0;
	      }
	      else {
          while ( ph->nxt && ph->nxt < hash->nhmax ) {
						if ( ph->min == mins && ph->max == maxs )  break;
            ph = &hash->item[ph->nxt];
          }
					if ( ph->min == mins && ph->max == maxs )  continue;
					ph->nxt = hash->hnxt;
          ph      = &hash->item[hash->hnxt];
          ph->min = mins;
          ph->max = maxs;
          ph->nxt = 0;
          ++hash->hnxt;
					assert(hash->hnxt < hash->nhmax );
        }
      }
    }
  }      

  return(hash);
}


int hashel_3d(pMesh mesh) {
  pTetra    pt,pt1;
  pPoint    ppt;
  int       k,kk,pp,l,ll,mins,mins1,maxs,maxs1,sum,sum1,iadr;
  int      *hcode,*link,*adja,inival,hsize;
  unsigned char   i,ii,i1,i2,i3;
  unsigned int    key;

  /* memory alloc */
  hcode = (int*)M_calloc(mesh->ne+1,sizeof(int),"hash");
  assert(hcode);
  link  = mesh->adja;
  hsize = mesh->ne;

  /* init */
	inival = 2147483647;
  for (k=0; k<=mesh->ne; k++)
    hcode[k] = -inival;

  /* build hash table */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !pt->v[0] )  continue;
    for (i=0; i<4; i++) {
      i1 = idirt[i+1];
      i2 = idirt[i+2];
      i3 = idirt[i+3];
      mins = MS_MIN(pt->v[i1],pt->v[i2]);
      mins = MS_MIN(mins,pt->v[i3]);
      maxs = MS_MAX(pt->v[i1],pt->v[i2]);
      maxs = MS_MAX(maxs,pt->v[i3]);

      /* compute key */
      sum = pt->v[i1] + pt->v[i2] + pt->v[i3];
      key = KA*mins + KB*maxs + KC*sum;
      key = key % hsize + 1;

      /* insert */
      iadr = 4*(k-1) + i+1;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }

  /* set adjacency */
  for (l=4*mesh->ne; l>0; l--) {
    if ( link[l] >= 0 )  continue;
    k = (l-1) / 4 + 1;
    i = (l-1) % 4;
    i1 = idirt[i+1];
    i2 = idirt[i+2];
    i3 = idirt[i+3];
    pt = &mesh->tetra[k];

    sum  = pt->v[i1] + pt->v[i2] + pt->v[i3];
    mins = MS_MIN(pt->v[i1],pt->v[i2]);
    mins = MS_MIN(mins,pt->v[i3]);
    maxs = MS_MAX(pt->v[i1],pt->v[i2]);
    maxs = MS_MAX(maxs,pt->v[i3]);

    /* accross link */
    ll = -link[l];
    pp = 0;
    link[l] = 0;
    while ( ll != inival ) {
      kk = (ll-1) / 4 + 1;
      ii = (ll-1) % 4;
      i1 = idirt[ii+1];
      i2 = idirt[ii+2];
      i3 = idirt[ii+3];
      pt1  = &mesh->tetra[kk];
      sum1 = pt1->v[i1] + pt1->v[i2] + pt1->v[i3];
      if ( sum1 == sum ) {
        mins1 = MS_MIN(pt1->v[i1],pt1->v[i2]);
        mins1 = MS_MIN(mins1,pt1->v[i3]);
        if ( mins1 == mins ) {
          maxs1 = MS_MAX(pt1->v[i1],pt1->v[i2]);
          maxs1 = MS_MAX(maxs1,pt1->v[i3]);
          if ( maxs1 == maxs ) {
            /* adjacent found */
            if ( pp != 0 )  link[pp] = link[ll];
            link[l] = 4*kk + ii;
            link[ll]= 4*k + i;
            break;
          }
        }
      }
      pp = ll;
      ll = -link[ll];
    }
  }
  M_free(hcode);

  /* set seed */
  for (k=1; k<=mesh->ne; k++) {
    pt   = &mesh->tetra[k];
    iadr = 4*(k-1)+1;
    adja = &mesh->adja[iadr];
    if ( !adja[0] )  mesh->point[pt->v[1]].s = k;
    if ( !adja[1] )  mesh->point[pt->v[2]].s = k;
    if ( !adja[2] )  mesh->point[pt->v[3]].s = k;
    if ( !adja[3] )  mesh->point[pt->v[0]].s = k;
  }
	
  for (k=1; k<=mesh->ne; k++) {
	  pt = &mesh->tetra[k];
    for (i=0; i<4; i++) {
	    ppt = &mesh->point[pt->v[i]];
	    if ( !ppt->s )  ppt->s = k;	
    }
  }

  return(1);
}


int hashel_2d(pMesh mesh) {
  pTria     pt,pt1;
  int       k,kk,pp,l,ll,mins,mins1,maxs,maxs1,iadr;
  int      *hcode,*link,*adj,inival,hsize,iad;
  unsigned char  *hvoy,i,ii,i1,i2;
  unsigned int    key;

  if ( !mesh->nt )  return(0);

  /* memory alloc */
  hcode = (int*)M_calloc(mesh->nt+1,sizeof(int),"hash");
  assert(hcode);
  link  = mesh->adja;
  hsize = mesh->nt;
  hvoy  = (unsigned char*)hcode;

  /* init */
	inival = 2147483647;
  for (k=0; k<=mesh->nt; k++)
    hcode[k] = -inival;

  /* build hash table */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] )  continue;
    for (i=0; i<3; i++) {
      i1 = idir[i+1];
      i2 = idir[i+2];
      if ( pt->v[i1] < pt->v[i2] ) {
        mins = pt->v[i1];
        maxs = pt->v[i2];
      }
      else {
        mins = pt->v[i2];
        maxs = pt->v[i1];
      }

      /* compute key */
      key = KTA*mins + KTB*maxs;
      key = key % hsize + 1;

      /* insert */
      iadr = 3*(k-1) + i+1;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }

  /* set adjacency */
  for (l=3*mesh->nt; l>0; l--) {
    if ( link[l] >= 0 )  continue;
    k = (l-1) / 3 + 1;
    i = (l-1) % 3;
    i1 = idir[i+1];
    i2 = idir[i+2];
    pt = &mesh->tria[k];

    mins = MS_MIN(pt->v[i1],pt->v[i2]);
    maxs = MS_MAX(pt->v[i1],pt->v[i2]);

    /* accross link */
    ll = -link[l];
    pp = 0;
    link[l] = 0;
    hvoy[l] = 0;
    while ( ll != inival ) {
      kk = (ll-1) / 3 + 1;
      ii = (ll-1) % 3;
      i1 = idir[ii+1];
      i2 = idir[ii+2];
      pt1  = &mesh->tria[kk];
      if ( pt1->v[i1] < pt1->v[i2] ) {
        mins1 = pt1->v[i1];
        maxs1 = pt1->v[i2];
      }
      else {
        mins1 = pt1->v[i2];
        maxs1 = pt1->v[i1];
      }
      
      if ( mins1 == mins  && maxs1 == maxs ) {
        /* adjacent found */
        if ( pp != 0 )  link[pp] = link[ll];
        link[l] = 3*kk + ii;
        link[ll]= 3*k + i;
        break;
      }
      pp = ll;
      ll = -link[ll];
    }
  }
  M_free(hcode);

  for (k=1; k<=mesh->nt; k++) {
    pt  = &mesh->tria[k];
    iad = 3*(k-1)+1;
    adj = &mesh->adja[iad];
    if ( !adj[0] )  mesh->point[pt->v[1]].s = k;
    if ( !adj[1] )  mesh->point[pt->v[2]].s = k;
    if ( !adj[2] )  mesh->point[pt->v[0]].s = k;
  }

  return(1);
}

