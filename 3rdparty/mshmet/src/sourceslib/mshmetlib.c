/* mshmetlib.c

 mshmetlib(int option, ...) to use mshmet via a library
 * compute metric based on hessian

 * j.morice LJLL 2010
 * Copyright (c) LJLL, 2010.
*/

#include "mshmet.h"
#include "compil.date"
extern long verbosity;
char     idir[5]     = {0,1,2,0,1};
mytime   mshmet_ctim[TIMEMAX];
int   (*boulep)(pMesh ,int ,int ,int *);
int   (*hashel)(pMesh );
int   (*gradLS)(pMesh ,pSol ,int ,int );
int   (*hessLS)(pMesh ,pSol ,int ,int );
int   (*avgval)(pMesh ,pSol ,int );
int   (*clsval)(pMesh ,pSol ,int );
int   (*nrmhes)(pMesh ,pSol ,int );
int   (*redsim)(double *,double *,double *);
int   (*defmet)(pMesh ,pSol ,int );
double (*getSol)(pSol ,int ,int );
int   (*metrLS)(pMesh mesh,pSol );
int   (*lissag)(pMesh ,pSol , int ,int );


static void mshmet_excfun(int sigid) {
  fprintf(stdout,"\n Unexpected error:");  fflush(stdout);
  switch(sigid) {
    case SIGABRT:
      fprintf(stdout,"  Abnormal stop\n");  exit(1);
    case SIGFPE:
      fprintf(stdout,"  Floating-point exception\n"); exit(1);
    case SIGILL:
      fprintf(stdout,"  Illegal instruction\n"); exit(1);
    case SIGSEGV:
      fprintf(stdout,"  Segmentation fault\n");  exit(1);
    case SIGTERM:
    case SIGINT:
      //fprintf(stdout,"  Program killed\n");  exit(1);
      fprintf(stdout," Abnormal end\n");  exit(1);
  }
  exit(1);
}

/*
static void usage(char *prog) {
  fprintf(stdout,"\n usage: %s filein[.mesh] [solin[.sol]] [fileout.sol] -eps x -hmin y -hmax z -v -iso -norm\n",prog);
  
  fprintf(stdout,"\n** Generic options :\n");
  fprintf(stdout,"-d      Turn on debug mode\n");
  fprintf(stdout,"-h      Print this message\n");
  fprintf(stdout,"-ls     Build levelset metric\n");
  fprintf(stdout,"-v [n]  Tune level of verbosity\n");
  fprintf(stdout,"-m file Use metric file\n");

  fprintf(stdout,"\n** Specific options : \n");
  fprintf(stdout,"  -eps :  tolerance\n");
  fprintf(stdout,"  -hmin:  min size\n");
  fprintf(stdout,"  -hmax:  max size\n");
  fprintf(stdout,"  -iso :  isotropic\n");
  fprintf(stdout,"  -w   :  relative width for LS (0<w<1)\n");
  fprintf(stdout,"  -n[i]:  normalization (level i), default=0\n");
  fprintf(stdout,"  -s n :  consider solution n (only)\n");

  fprintf(stdout,"\n** DEFAULT.hmet file allows to store parameter values\n");
  fprintf(stdout,"eps   x\n");
  fprintf(stdout,"hmin  y\n");
  fprintf(stdout,"hmax  z\n");
  fprintf(stdout,"iso\n");
  exit(1);
}
*/

/*
int parsop(pMesh mesh,pSol sol) {
  char    *ptr,data[256],key[256];
  float    dummy;
  int      i,ret;
  FILE    *in;

  strcpy(data,sol->name);
  ptr = strstr(data,".sol");
  if ( ptr )  *ptr = '\0';
  strcat(data,".mhes");

  in = fopen(data,"r");
  if ( !in ) {
    strcpy(data,"DEFAULT.hmet");
    in = fopen(data,"r");
    if ( !in )  {
      if ( mesh->info.imprim < 0 )
        fprintf(stdout,"  %%%% DEFAULT VALUES (%g %g %g)\n",
                mesh->info.eps,mesh->info.hmin,mesh->info.hmax);
      return(1);
    }
  }
  fprintf(stdout,"  %%%% %s FOUND\n",data);

  while ( !feof(in) ) {
    ret = fscanf(in,"%s",key);
    if ( !ret || feof(in) )  break;
    for (i=0; i<strlen(key); i++) key[i] = tolower(key[i]);

    if ( !strcmp(key,"hmin") ) {
      fscanf(in,"%f",&dummy);
      mesh->info.hmin = dummy;
    }
    else if ( !strcmp(key,"hmax") ) {
      fscanf(in,"%f",&dummy);
      mesh->info.hmax = dummy;
    }
    else if ( !strcmp(key,"eps") ) {
      fscanf(in,"%f",&dummy);
      mesh->info.eps = dummy;
    }
    else if ( !strcmp(key,"iso") ) {
      mesh->info.iso = 1;
    }
    else if ( !strcmp(key,"norm") ) {
      fscanf(in,"%d",&mesh->info.nnu);
    }
    else if ( key[0] == '#' ) {
      fgets(key,255,in);
    }
    else fprintf(stderr,"  unrecognized keyword : %s\n",key);
  }

  fclose(in);
  return(1);
}
*/

static void mshmet_stats(pMesh mesh,pSol sol) {
  fprintf(stdout,"     NUMBER OF GIVEN VERTICES   %8d\n",mesh->np);
  if ( mesh->nt )
    fprintf(stdout,"     NUMBER OF GIVEN TRIANGLES  %8d\n",mesh->nt);
  if ( mesh->ne )
    fprintf(stdout,"     NUMBER OF GIVEN TETRAHEDRA %8d\n",mesh->ne);
  fprintf(stdout,"     NUMBER OF GIVEN DATA       %8d\n",sol->np);
}


static void mshmet_endcod() {
  double   ttot,ttim[TIMEMAX];
  int      k,call[TIMEMAX];

  chrono(OFF,&mshmet_ctim[0]);
  for (k=0; k<TIMEMAX; k++) {
    call[k] = mshmet_ctim[k].call;
    ttim[k] = mshmet_ctim[k].call ? gttime(mshmet_ctim[k]) : 0.0;
  }
  ttot    = ttim[1]+ttim[2]+ttim[3]+ttim[4];
  ttim[0] = MS_MAX(ttim[0],ttot);

  if ( ttot > 0.01 ) {
    fprintf(stdout,"\n  -- CPU REQUIREMENTS\n");
    fprintf(stdout,"  in/out    %8.2f %%    %3d. calls,   %7.2f sec/call\n",
        100.*ttim[1]/ttim[0],call[1],ttim[1]/(float)call[1]);
    fprintf(stdout,"  analysis  %8.2f %%    %3d. calls,   %7.2f sec/call\n",
        100.*ttim[2]/ttim[0],call[2],ttim[2]/(float)call[2]);
    fprintf(stdout,"  metric    %8.2f %%    %3d. calls,   %7.2f sec/call\n",
        100.*ttim[3]/ttim[0],call[3],ttim[3]/(float)call[3]);
    fprintf(stdout,"  total     %8.2f %%    %3d. calls,   %7.2f sec/call\n",
        100.*ttot/ttim[0],call[0],ttot/(float)call[0]);
  }
  fprintf(stdout,"\n   ELAPSED TIME  %.2f SEC.  (%.2f)\n",ttim[0],ttot);
}


/* set function pointers */
/* set function pointers */
void MSHMET_setfunc(pMesh mesh) {
  if ( mesh->dim == 2 ) {
    boulep = boulep_2d;
    hashel = hashel_2d;
    gradLS = gradLS_2d;
    hessLS = hessLS_2d;
    getSol = getSol_2d;
    avgval = avgval_2d;
    clsval = clsval_2d;
    nrmhes = nrmhes_2d;
    defmet = defmet_2d;
    redsim = redsim_2d;
    metrLS = metrLS_2d;
    lissag = lissag_2d;
  }
  else {
    if ( mesh->ne > 0 ) { /* 3d */
      boulep = boulep_3d;
      hashel = hashel_3d;
      gradLS = gradLS_3d;
      hessLS = hessLS_3d;
      getSol = getSol_3d;
      avgval = avgval_3d;
      clsval = clsval_3d;
      nrmhes = nrmhes_3d;
      defmet = defmet_3d;
      redsim = redsim_3d;
      metrLS = metrLS_3d;
			lissag = lissag_3d;
    }
    else { /* surface mesh */
      boulep = boulep_2d;
      hashel = hashel_2d;
      lissag = lissag_2d;
      avgval = avgval_3d;
      clsval = clsval_3d;
      nrmhes = nrmhes_3d;
      getSol = getSol_3d;
      redsim = redsim_3d;
      gradLS = gradLS_s;
      hessLS = hessLS_s;
      defmet = defmet_s;

      metrLS = metrLS_3d;
    }
  }
}



int MSHMET_mshmet(int intopt[7], double fopt[4], pMesh mesh, pSol sol){
  Info *info;
  if ( intopt[4] ) {
    fprintf(stdout,"  -- MSHMET, Release %s (%s) \n",MS_VER,MS_REL);
    fprintf(stdout,"     %s\n",MS_CPY);
    fprintf(stdout,"    %s\n",COMPIL);
  }

  /* trap exceptions */
  signal(SIGABRT,mshmet_excfun);
  signal(SIGFPE,mshmet_excfun);
  signal(SIGILL,mshmet_excfun);
  signal(SIGSEGV,mshmet_excfun);
  signal(SIGTERM,mshmet_excfun);
  signal(SIGINT,mshmet_excfun);
  //atexit(mshmet_endcod);

  tminit(mshmet_ctim,TIMEMAX);
  chrono(ON,&mshmet_ctim[0]);
  chrono(ON,&mshmet_ctim[1]);
  /* default */
  info = &mesh->info;
  info->hmin   = (float) fopt[0]; // 0.01;
  info->hmax   = (float) fopt[1]; // 1.0;
  info->eps    = (float) fopt[2]; // 0.01;
  info->width  = (float) fopt[3]; // 0.05;

  info->nnu    = intopt[0]; //  0;
  info->iso    = intopt[1]; //  0;
  info->ls     = intopt[2]; //  0;
  info->ddebug = intopt[3]; //  0;
  info->imprim = intopt[4]; // 10;
  info->nlis   = intopt[5]; //  0;
  
  info->bin    =   1;          // pas besoin c'est pour le fichier
  info->nsol   =  -1; //-1;    // pas besoin ==> on peut prendre plusieurs solutions en meme temps ???
  info->metric =  intopt[6];        // 0; // metric given besoin ???
 
  MSHMET_setfunc(mesh);
  chrono(OFF,&mshmet_ctim[1]);
  if ( mesh->info.imprim ) {
    mshmet_stats(mesh,sol);
    fprintf(stdout,"  -- DATA READING COMPLETED.     %.2f sec.\n",gttime(mshmet_ctim[1]));

    fprintf(stdout,"\n  %s\n   MODULE MSHMET-LJLL : %s (%s)\n  %s\n",
            MS_STR,MS_VER,MS_REL,MS_STR);
  }

  /* analysis */
  chrono(ON,&mshmet_ctim[2]);
  if ( mesh->info.imprim )   fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  if ( abs(mesh->info.imprim) > 4 ) {
    fprintf(stdout,"  ** SETTING ADJACENCIES\n");
    fflush(stdout);
  }
  if ( !scaleMesh(mesh,sol) )  return(1);
  if ( !hashel(mesh) )         exit(1);
  chrono(OFF,&mshmet_ctim[2]);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %.2f sec.\n",gttime(mshmet_ctim[2]));

  /* metric */
  chrono(ON,&mshmet_ctim[3]);
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- PHASE 2 : METRIC\n");
  if ( !mshme1(mesh,sol) )  exit(1);
  
  chrono(OFF,&mshmet_ctim[3]);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %.2f sec.\n",gttime(mshmet_ctim[3]));
  
  if ( mesh->info.imprim )
    fprintf(stdout,"\n  %s\n   END OF MODULE MSHMET \n  %s\n",MS_STR,MS_STR);
  /*
  sol->outn="zzzz";
  if ( !saveMet(sol,&mesh->info,sol->outn) )  exit(1);
  */
  if ( mesh->info.imprim ) {
    mshmet_endcod();

    fprintf(stdout,"\n  %s\n   END OF MODULE MSHMET \n  %s\n",MS_STR,MS_STR);
  }
  if ( mesh->info.imprim < -4 || mesh->info.ddebug )  M_memDump();

  return(0);
}
