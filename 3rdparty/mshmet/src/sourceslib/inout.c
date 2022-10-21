#include "mshmet.h"

extern int eigenv(int symmat,double *mat,double lambda[3],double v[3][3]);
extern int eigen2(double *mm,double *lambda,double vp[2][2]);


/* read mesh data */
int loadMesh(pMesh mesh,char *filename) {
  pPoint       ppt;
  pTetra       pt;
  pTria        pt1;
  float        fp1,fp2,fp3;
  int          i,k,ref;
  int64_t inm;
  char        *ptr,data[128];

  strcpy(data,filename);
  ptr = strstr(data,".mesh");
  mesh->info.bin = 0;
  if ( !ptr ) {
    strcat(data,".meshb");
    mesh->info.bin = 1;
    if ( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if ( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
        return(0);
      }
      mesh->info.bin = 0;
    }
  }
  else if ( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
    fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"  -- READING DATA FILE %s\n",data);

  mesh->np = GmfStatKwd(inm,GmfVertices);
  mesh->nt = GmfStatKwd(inm,GmfTriangles);
  mesh->ne = GmfStatKwd(inm,GmfTetrahedra);
  if ( !mesh->np || mesh->ne+mesh->nt == 0 ) {
    fprintf(stdout,"  ** MISSING DATA\n");
    return(0);
  }

  /* mem alloc */
  mesh->point = (pPoint)M_calloc(mesh->np+1,sizeof(Point),"point");
  assert(mesh->point);
  if ( mesh->ne ) {
    mesh->tetra = (pTetra)M_calloc(mesh->ne+1,sizeof(Tetra),"tetra");
    assert(mesh->tetra);
    mesh->adja = (int*)M_calloc(4*mesh->ne+5,sizeof(int),"adja");
    assert(mesh->adja);
    mesh->nt = 0;
  }
  else if ( mesh->nt ) {
    mesh->tria  = (pTria)M_calloc(mesh->nt+1,sizeof(Tria),"tria");
    assert(mesh->tria);
    mesh->adja = (int*)M_calloc(3*mesh->nt+5,sizeof(int),"adja");
    assert(mesh->adja);
  }

  /* read mesh vertices */
  GmfGotoKwd(inm,GmfVertices);
  if ( mesh->dim == 2 ) {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( mesh->ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ref);
    }
  }
  else {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( mesh->ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&fp3,&ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
        ppt->c[2] = fp3;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->c[2],&ref);
    }
  }

  /* read mesh triangles */
  GmfGotoKwd(inm,GmfTriangles);
  for (k=1; k<=mesh->nt; k++) {
    pt1 = &mesh->tria[k];
    GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&ref);
    for (i=0; i<3; i++) {    
      ppt = &mesh->point[pt1->v[i]];
      if ( !ppt->s )  ppt->s = k;
    }
  }

  /* read mesh tetrahedra */
  GmfGotoKwd(inm,GmfTetrahedra);
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    GmfGetLin(inm,GmfTetrahedra,&pt->v[0],&pt->v[1],&pt->v[2],&pt->v[3],&ref);
    for (i=0; i<4; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( !ppt->s )  ppt->s = k;
    }
  }

  GmfCloseMesh(inm);
  return(1);
}


/* load solution (metric) */
int loadSol(pSol sol,Info *info,char *filename) {
  double       dbuf[ GmfMaxTyp ];
  float        fbuf[ GmfMaxTyp ];
  int          k,i,ia;
  int64_t inm;
  char        *ptr,data[128];

  strcpy(data,filename);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  strcat(data,".solb");
  if (!(inm = GmfOpenMesh(data,GmfRead,&sol->ver,&sol->dim)) ) {
    ptr  = strstr(data,".solb");
    *ptr = '\0';
    strcat(data,".sol");
    if (!(inm = GmfOpenMesh(data,GmfRead,&sol->ver,&sol->dim)) ) {
      fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
      return(0);
    }
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  if ( abs(info->imprim) > 3 )
    fprintf(stdout,"  -- READING DATA FILE %s\n",data);

  sol->np = GmfStatKwd(inm,GmfSolAtVertices,&sol->type,&sol->size,sol->typtab);
  if ( !sol->np ) {
    fprintf(stdout,"  ** MISSING DATA.\n");
    return(0);
  }

  /* mem alloc */
  sol->sol = (double*)M_calloc(sol->np+1,sol->size*sizeof(double),"inout");
  assert(sol->sol);

  /* read mesh solutions */
  GmfGotoKwd(inm,GmfSolAtVertices);
  for (k=1; k<=sol->np; k++) {
    if ( sol->ver == GmfFloat ) {
      GmfGetLin(inm,GmfSolAtVertices,fbuf);
      ia = (k-1)*sol->size + 1;
      for (i=0; i<sol->size; i++)
        sol->sol[ia+i] = fbuf[i];   // < -0.0001 ? -1.0 : fbuf[i];
    }
    else {
      GmfGetLin(inm,GmfSolAtVertices,dbuf);
      ia = (k-1)*sol->size + 1;
      for (i=0; i<sol->size; i++)
        sol->sol[ia+i] = dbuf[i];
    }
  }

  GmfCloseMesh(inm);
  return(1);  
}


/* load solution (metric) */
int loadMetric(pSol sol,Info *info,char *filename) {
  double   dbuf[ GmfMaxTyp ],tmpd;
  float    fbuf[ GmfMaxTyp ],tmpf;
  int      i,k,ia,np,ver,dim,type,size,typtab[ GmfMaxTyp ];
  int64_t inm;
  if ( !(inm = GmfOpenMesh(filename,GmfRead,&ver,&dim)) ) {
    fprintf(stderr,"  ** %s  NOT FOUND.\n",filename);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",filename);

  if ( abs(info->imprim) > 3 )
    fprintf(stdout,"  -- READING METRIC FILE %s\n",filename);

  np = GmfStatKwd(inm,GmfSolAtVertices,&type,&size,typtab);
  if ( !np || np != sol->np ) {
    fprintf(stdout,"  ** INCOMPLETE DATA.\n");
    return(0);
  }

  /* mem alloc */
  sol->met = (double*)M_calloc(np+1,size*sizeof(double),"inout");
  assert(sol->met);

  /* read mesh solutions */
  GmfGotoKwd(inm,GmfSolAtVertices);
  for (k=1; k<=np; k++) {
    if ( ver == GmfFloat ) {
      GmfGetLin(inm,GmfSolAtVertices,fbuf);
      ia = (k-1)*size + 1;
      if ( sol->dim == 3 ) {
        tmpf    = fbuf[3];
        fbuf[3] = fbuf[2];
        fbuf[2] = tmpf;
      }
      for (i=0; i<size; i++)
        sol->met[ia+i] = fbuf[i];
    }
    else {
      GmfGetLin(inm,GmfSolAtVertices,dbuf);
      ia = (k-1)*size + 1;
      if ( sol->dim == 3 ) {
        tmpd    = dbuf[3];
        dbuf[3] = dbuf[2];
        dbuf[2] = tmpd;
      }
      for (i=0; i<size; i++)
        sol->met[ia+i] = dbuf[i];
    }
  }

  GmfCloseMesh(inm);
  return(1);  
}


int saveMet(pSol sol,Info *info,char *filename) {
  double       dbuf[ GmfMaxTyp ],tmpd,hmin,lambda[3],vp2[2][2],vp3[3][3];
  float        fbuf[ GmfMaxTyp ],tmpf;
  int          k,i,ia,size;
  int64_t outm;
  char        *ptr,data[128];

  strcpy(data,filename);
  ptr = strstr(data,".mesh");
  if ( ptr )  {
    *ptr = '\0';
    strcat(data,".solb");
  }
  else {
    ptr = strstr(data,".sol");
    if ( !ptr )  strcat(data,".sol");
  }
  if (!(outm = GmfOpenMesh(data,GmfWrite,sol->ver,sol->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /* write sol */
  sol->type = 1;
  if ( info->iso ) {
    sol->typtab[0] = GmfSca;
    sol->size = 1;
  }
  else {
    sol->typtab[0] = GmfSymMat;
    sol->size = sol->dim == 2 ? 3 : 6;
  }
  GmfSetKwd(outm,GmfSolAtVertices,sol->np,sol->type,sol->typtab);
  for (k=1; k<=sol->np; k++) {
    ia = (k-1)*sol->size + 1;
    if ( sol->ver == GmfFloat ) {
      for (i=0; i<sol->size; i++)
        fbuf[i] = sol->met[ia+i];
      if ( sol->dim == 3 ) {
        tmpf    = fbuf[3];
        fbuf[3] = fbuf[2];
        fbuf[2] = tmpf;
      }
      GmfSetLin(outm,GmfSolAtVertices,fbuf);
    }
    else {
      for (i=0; i<sol->size; i++)
        dbuf[i] = sol->met[ia+i];
      if ( sol->dim == 3 ) {
        tmpd    = dbuf[3];
        dbuf[3] = dbuf[2];
        dbuf[2] = tmpd;
      }
      GmfSetLin(outm,GmfSolAtVertices,dbuf);
    }
  }
  GmfCloseMesh(outm);

	if ( !sol->mapname )  return(1);
	
	/* save sizemap iso */ 
	strcpy(data,sol->mapname);
  ptr = strstr(data,".mesh");
  if ( ptr )  {
    *ptr = '\0';
    strcat(data,".solb");
  }
  else {
    ptr = strstr(data,".sol");
    if ( !ptr )  strcat(data,".sol");
  }
  if ( !(outm = GmfOpenMesh(data,GmfWrite,sol->ver,sol->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  sol->type = 1;
  sol->typtab[0] = GmfSca;
  size = sol->dim*(sol->dim+1)/2;
  GmfSetKwd(outm,GmfSolAtVertices,sol->np,sol->type,sol->typtab);
  if ( sol->dim == 2 ) {
    for (k=1; k<=sol->np; k++) {
      ia = (k-1)*size + 1;
			eigen2(&sol->met[ia],lambda,vp2);
			lambda[0] = 1.0 / sqrt(lambda[0]);
			lambda[1] = 1.0 / sqrt(lambda[1]);
			hmin = MS_MIN(lambda[0],lambda[1]);
      if ( sol->ver == GmfFloat ) {
				fbuf[0] = hmin;
        GmfSetLin(outm,GmfSolAtVertices,fbuf);
      }
      else
        GmfSetLin(outm,GmfSolAtVertices,&hmin);
    }
	}
	else {
    for (k=1; k<=sol->np; k++) {
      ia = (k-1)*size + 1;
      eigenv(1,&sol->met[ia],lambda,vp3);
			lambda[0] = 1.0 / sqrt(lambda[0]);
			lambda[1] = 1.0 / sqrt(lambda[1]);
			lambda[2] = 1.0 / sqrt(lambda[2]);
			hmin = MS_MIN(lambda[0],lambda[1]);
			hmin = MS_MIN(hmin,lambda[2]);
      if ( sol->ver == GmfFloat ) {
        fbuf[0] = hmin;
        GmfSetLin(outm,GmfSolAtVertices,fbuf);
      }
      else
        GmfSetLin(outm,GmfSolAtVertices,&hmin);
    }
  }

  GmfCloseMesh(outm);
  return(1);
}


int saveSol(pSol sol,Info *info,char *filename) {
  double       dbuf[ GmfMaxTyp ];
  float        fbuf[ GmfMaxTyp ];
  int          k,i,ia;
  int64_t outm;
  char        *ptr,data[128];

  strcpy(data,filename);
  ptr = strstr(data,".mesh");
  if ( ptr )  {
    *ptr = '\0';
    strcat(data,".solb");
  }
  else {
    ptr = strstr(data,".sol");
    if ( !ptr )  strcat(data,".sol");
  }

  if (!(outm = GmfOpenMesh(data,GmfWrite,sol->ver,sol->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
    return(0);
  }

  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /* write sol */
  sol->type = 1;
  sol->typtab[0] = GmfSca;
  sol->size = 1;
  GmfSetKwd(outm,GmfSolAtVertices,sol->np,sol->type,sol->typtab);
  for (k=1; k<=sol->np; k++) {
    ia = (k-1)*sol->size + 1;
    if ( sol->ver == GmfFloat ) {
      for (i=0; i<sol->size; i++)
        fbuf[i] = sol->sol[ia+i];
      GmfSetLin(outm,GmfSolAtVertices,fbuf);
    }
    else {
      for (i=0; i<sol->size; i++)
        dbuf[i] = sol->sol[ia+i];
      GmfSetLin(outm,GmfSolAtVertices,dbuf);
    }
  }

  GmfCloseMesh(outm);
  return(1);
}


int outder(pMesh mesh,pSol sol) { 
  double    *grd,*hes,*np;
  int        k,iadr,ver,type,siz,typtab[GmfMaxTyp];
  int64_t outm;
  ver = GmfDouble;

  /* gradients */
  if ( !(outm = GmfOpenMesh("gradient.sol",GmfWrite,ver,sol->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN gradient\n");
    return(0);
  }
  type = 1;
  typtab[0] = GmfVec;
  GmfSetKwd(outm,GmfSolAtVertices,mesh->np,type,typtab);
  for (k=1; k<=mesh->np; k++) {
	  iadr = (k-1)*sol->dim + 1;
	  grd  = &sol->grd[iadr];
    GmfSetLin(outm,GmfSolAtVertices,grd);
  }
  GmfCloseMesh(outm);

  /* hessian */
  if ( !(outm = GmfOpenMesh("hessien.sol",GmfWrite,ver,sol->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN hessien\n");
    return(0);
  }
  type = 1;
  typtab[0] = GmfSymMat;
	siz = sol->dim == 2 ? 3 : (mesh->ne > 0 ? 6 : 3);

  GmfSetKwd(outm,GmfSolAtVertices,mesh->np,type,typtab);
  for (k=1; k<=mesh->np; k++) {
    iadr = (k-1)*siz + 1;
    hes = &sol->hes[iadr];
    GmfSetLin(outm,GmfSolAtVertices,hes);
  }
  GmfCloseMesh(outm);

  if ( mesh->dim == 3 && mesh->ne == 0 ) {
    if ( !(outm = GmfOpenMesh("normal.sol",GmfWrite,ver,sol->dim)) ) {
      fprintf(stderr,"  ** UNABLE TO OPEN gradient\n");
      return(0);
    }
    type = 1;
    typtab[0] = GmfVec;
    GmfSetKwd(outm,GmfSolAtVertices,mesh->np,type,typtab);
    for (k=1; k<=mesh->np; k++) {
	    iadr = (k-1)*sol->dim + 1;
	    np   = &sol->nn[iadr];
      GmfSetLin(outm,GmfSolAtVertices,np);
    }
    GmfCloseMesh(outm);
  }

	return(1);
}


