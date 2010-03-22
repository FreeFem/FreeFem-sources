#define  __YAMSLIB

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>

#include "yams.h"
#include "defines.h"
#include "sproto.h"


/* globals (see globals.h) */
Error       yerr;
Info        info;
Options     opts;
pHashtable  hash;
mytime      ctim[TIMEMAX];

long      nhmax,hnext,hsize;
int       out; //,yams_idir[5] = {0,1,2,0,1},yams_idirq[7] = {0,1,2,3,0,1,2};
short     imprim;
ubyte     ddebug;
ubyte     ecp;


static void yams_excfun(int sigid) {
  switch(sigid){
  case SIGFPE:
    fprintf(stderr,"  ## FP EXCEPTION. STOP\n");
    break;
  case SIGILL:
    fprintf(stderr,"  ## ILLEGAL INSTRUCTION. STOP\n");
    break;
  case SIGSEGV:
    fprintf(stderr,"  ## SEGMENTATION FAULT. STOP\n");
    break;
  case SIGABRT:
  case SIGTERM:
  case SIGINT:
    fprintf(stderr,"  ## ABNORMAL END. STOP\n");
    break;
  }
  out = 0;
  exit(1);
}

static void yams_endcod() {
  chrono(OFF,&ctim[0]);
  chrono(OFF,&ctim[1]);
  E_dump();
  if ( out <= 0 ) {
    prierr(WAR,8002);
    fprintf(stdout,"\n   ELAPSED TIME  %.2f SEC.\n",gttime(ctim[0]));
  }
}

static void yams_inival(){
  /* initialize data */
  E_put("inival");
  info.dmin  = (double)FLT_MAX; 
  info.dmax  = (double)FLT_MIN;
  info.xmin  = info.ymin = info.zmin = (double)FLT_MAX;
  info.xmax  = info.ymax = info.zmax = (double)-FLT_MAX/2.;
  info.nedg  = info.nrid = info.ndang = 0;
  info.ncoi  = info.nreq = info.nvus  = 0;
  info.cc    = info.flip = 0;
  info.nulp  = info.nulf = info.nuln = 0;
  info.qpire = 0;
  info.manifold = TRUE;

  /* set default values for options */
  opts.hmin   = -2.0;    
  opts.hmax   = -2.0;
  opts.shock  =   1.3;       /* default mesh gradation     */
  opts.eps    =   0.01;      /* geometric approximation    */
  opts.iso    =   0.0;
  opts.declic =   1.0 / BETAC;
  opts.lambda =   -1.0;
  opts.mu     =   -1.0;
  opts.ridge  =   cos(RIDG*M_PI/180.);
  opts.geom   =   cos(GEOM*M_PI/180.);
  opts.walton =   COS45DEG;  /* Walton limitation          */
  opts.bande  =   -2;       /* default = 1 unit           */
  opts.degrad =   QUALCOE;   /* quality degradation        */
  opts.ctrl   =   REL | ISO;
  opts.iter   =   -1;
  opts.check  = 1;
  opts.alpha  = sqrt(opts.eps*(2.0-opts.eps));
  opts.gap    = 1 - opts.eps;
  
  opts.minnp  =   -1;
  opts.alpha  =   sqrt(opts.eps*(2.0-opts.eps));
  opts.gap    =   1.0 - opts.eps;

  E_pop();
}

void yams_printval() {
  /* set default values for options */
  printf("opts.hmin %f\n",opts.hmin);
  printf("opts.hmax %f\n",opts.hmax);
  printf("opts.kmin %f\n",opts.kmin);
  printf("opts.kmax %f\n",opts.kmax);
  printf("opts.eps %f\n",opts.eps);
  printf("opts.iso %f\n",opts.iso);
  printf(" opts.alpha %f\n", opts.alpha );
  printf(" opts.gap %f\n", opts.gap );
  printf(" opts.degrad %f\n", opts.degrad);
  printf("opts.ridge %f\n", opts.ridge);
  printf(" opts.geom %f\n", opts.geom);
  printf("opts.shock %f\n",opts.shock);
  printf(" opts.bande %f\n", opts.bande );
  printf(" opts.walton %f\n", opts.walton);
  printf("opts.declic %f\n", opts.declic);
  printf("opts.lambda %f\n",opts.lambda);
  printf("opts.mu %f\n",opts.mu);
   
  printf(" opts.ctrl %d\n", opts.ctrl );
  printf(" opts.iter %d\n", opts.iter );
  printf(" opts.choix %d\n", opts.choix );
  printf(" opts.minnp %d\n", opts.minnp );
  
  printf(" opts.check %X\n", (unsigned char) opts.check);
  printf(" opts.ptmult %X\n",  (unsigned char) opts.ptmult);
  printf(" opts.noreff %X\n",  (unsigned char) opts.noreff);
  printf(" opts.ffem %X\n",  (unsigned char) opts.ffem );
 
}


int yams_main(pSurfMesh sm, int intopt[23], double fopt[14], int infondang, int infocc ) {
  float       declic;
  float       ridge=RIDG;
  int         option,absopt,ret,memory,corr;
  int         choix;
  short       phase;
  int k;
  /* trap exceptions */
  signal(SIGABRT,yams_excfun);
  signal(SIGFPE,yams_excfun);
  signal(SIGILL,yams_excfun);
  signal(SIGSEGV,yams_excfun);
  signal(SIGTERM,yams_excfun);
  signal(SIGINT,yams_excfun);
  atexit(yams_endcod);

  /* init time and calls */
  tminit(ctim,TIMEMAX);
  chrono(ON,&ctim[0]);

  

  /* assign default values */
  yerr.lerror = FALSE;
  yerr.coderr = 0;
  phase  = 0;
  ret    = TRUE;
  out    = -1;
  memory = -1;
  imprim = -99;
  option = -99;
  choix  = option;
  ddebug = FALSE;
  declic = 0.009;
  ecp    = 0;
  
  // assigne option and surfacemesh
  
  /* setting defaults */
  sm->infile  = NULL;
  sm->outfile = NULL;
  sm->type    = M_SMOOTH | M_QUERY | M_DETECT | M_BINARY | M_OUTPUT;
  yams_inival();

  for (k=1; k<=sm->npfixe; k++) {
    pPoint ppt = &sm->point[k];
    /* find extrema coordinates */
    if ( ppt->c[0] < info.xmin ) info.xmin = ppt->c[0];
    if ( ppt->c[0] > info.xmax ) info.xmax = ppt->c[0];
    if ( ppt->c[1] < info.ymin ) info.ymin = ppt->c[1];
    if ( ppt->c[1] > info.ymax ) info.ymax = ppt->c[1];
    if ( ppt->c[2] < info.zmin ) info.zmin = ppt->c[2];
    if ( ppt->c[2] > info.zmax ) info.zmax = ppt->c[2];
  }

  // info nuln et nulp
  info.nuln = 0; 
  for (k=1; k<=sm->nvfixe; k++) {
    pGeomSupp g0 = &sm->geom[ k ];    
    double dd = g0->vn[0]*g0->vn[0] + g0->vn[1]*g0->vn[1] + g0->vn[2]*g0->vn[2];
    if ( dd < 0.0 ) 
      info.nuln++;
  }
  info.nulp = 0;

  /* mark used vertices */
  for (k=1; k<=sm->nefixe; k++) {
    pTriangle pt1 = &sm->tria[k];
    int i;
    if ( pt1->v[0] )
      for (i=0; i<3; i++) {
        pPoint ppt = &sm->point[pt1->v[i]];
        ppt->tag &= ~M_UNUSED;
      }
  }

  /* count unused vertices */
  for (k=1; k<=sm->npfixe; k++) {
    pPoint ppt;
    ppt = &sm->point[k];
    if ( ppt->tag & M_UNUSED )  info.nulp++;
  }


  /* get decimation parameters */
  opts.noreff = 0;
  opts.ffem   = 1;
  opts.ptmult = 0;

  /*
    intopt : 0  !! anisotropie
             1  !! ecp 
             2  !! extended out put file
	     3  !! FE correction 
	     4  !! Formatted (ascii) output file
	     5  !! save metric file
	     6  !! msh2
	     7  !! Split multiple connected points
	     8  !! memory
	     9  !! connected component
	    10  !! vrml 
	    11  !! imprim
	    12  !! nm : Create point on straight edge (no mapping)
	    13  !! nc : No validity check during smoothing (opt. 9)
	    14  !! np : Specify number of points desired
	    15  !! nit : Nb Iter
	    16  !! nq  : Output quads
	    17  !! nr  : No ridge detection
	    18  !! ns  : No point smoothing
	    19  !! no  : No output file
	    20  !! ref : Ignore face references
	    // rajouter lors de l'ouverture du fichiers yams
	    21  !! absolute : opts.ctrl &= ~REL;
	    22  !! set optim option

    fopt   : 0  !! iso 
             1  !! eps 
	     2  // pas de valeur
	     3  !! opts.lambda
	     4  !! opts.mu
	     5  // pas de valeur
	     6  !! hgrad  :: opts.shock
	     7  !! hmin   :: opts.hmin
	     8  !! hmax   :: opts.hmax
	     // rajouter lors de l'ouverture du fichiers yams
	     9  !! tolerance :: opts.bande
	     10 !! degrad :: opts.degrad
	     11 !! declic :: opts.declic 
	     12 !! walton :: opts.walton = cos(dummy/180.0*M_PI);
	     13 !! ridge  :: opts.ridge
   */
  if( intopt[0] == 1)
    opts.ctrl ^= ISO;
  opts.iso = fopt[0];  
  if( intopt[1] == 1 ) { ecp = 1;   sm->type &= ~M_BINARY; }
  opts.eps = fopt[1];
  if( intopt[2] == 1 )  sm->type |= M_EXTEND;
  if( intopt[3] == 1 )  opts.ffem = 0;
  if( intopt[4] == 1 )  sm->type &= ~M_BINARY;
  if( intopt[5] == 1 )  sm->type |= M_METRIC;
  if( intopt[6] == 1 ){
    sm->type |=  M_MSH2;
    sm->type &= ~M_BINARY;
    sm->type &= ~M_EXTEND; 
  }
  if( intopt[7] == 1 ){
    opts.ptmult = 1;
  }
  memory = intopt[8]; 
  sm->connex = intopt[9]; // a initialiser à -1 par défault
  if( intopt[10] == 1 ){
    sm->type |=  M_VRML;
    sm->type &= ~M_BINARY;
    sm->type &= ~M_EXTEND;
  }
  imprim = intopt[11]; 
  // parsar -n
  if( intopt[12] == 1 ) sm->type &= ~M_QUERY;
  if( intopt[13] == 1 ) opts.check = 0;
  opts.minnp = intopt[14];
  opts.iter = intopt[15];
  if( intopt[16] == 1 ) sm->type |= M_QUADS;
  if( intopt[17] == 1 ) sm->type &= ~M_DETECT;
  if( intopt[18] == 1 ) sm->type &= ~M_SMOOTH;
  //if( intopt[19] == 1 ) sm->type &= ~M_OUTPUT;
  sm->type &= ~M_OUTPUT;
  // parsar -r 
  if( intopt[20] == 1 ) opts.noreff = 1;
  // parsar -l
  opts.lambda  =  fopt[3];
  opts.mu      =  fopt[4];
  // parsar -O
  option = intopt[22];
  choix  = intopt[22];
  // parsar -h
  opts.shock = fopt[6];
  opts.hmin = fopt[7];
  opts.hmax = fopt[8];

  // fin parsar
  opts.choix = option;

  // yams0
  
  /* check option */
  if ( (option) > 0 )
    sm->type |= M_ENRICH;
  else
    memory = -1;
  /*
  if ( abs(*choix) > 4 && !(sm->type & M_QUADS) )
    sm->type &= ~M_SMOOTH;
  */
  if ( !(opts.ctrl & ISO) && abs(option) != 1 && abs(option) != 6 )
    opts.ctrl ^= ISO;

  if ( imprim )   fprintf(stdout,"  -- INPUT DATA\n");
  chrono(ON,&ctim[5]);
  
  opts.bande =  fopt[9];
  opts.degrad = fopt[10];
  if( intopt[21] == 1) opts.ctrl &= ~REL;

  // parsop check
  /* check parameters consistency */
  ridge = fopt[13];
  if ( ridge < 0.0 || !(sm->type & M_DETECT) )
    opts.ridge = -1.0;
  else
    opts.ridge  = cos(ridge*M_PI / 180.0);
  opts.degrad = min(opts.degrad,1.0);
  opts.degrad = max(opts.degrad,0.001);

  /* bound values */
  opts.alpha = sqrt(opts.eps * (2.-opts.eps));
  opts.gap   = 1.0 - opts.eps;
  if ( opts.walton < COS45DEG )  opts.walton = COS45DEG;

  // end assignement mesh and options
  //int bb = loadSol(sm,sm->infile);
  //sm->nmfixe = bb ? sm->npfixe : 0;
  absopt = abs(option);

  chrono(OFF,&ctim[5]);
  if ( imprim ) {
    fprintf(stdout,"     NUMBER OF GIVEN VERTICES    %8d\n",sm->npfixe);
    fprintf(stdout,"     NUMBER OF GIVEN TRIANGLES   %8d\n",sm->nefixe);
    fprintf(stdout,"  -- DATA READING COMPLETED.     %.2f sec.\n",
            gttime(ctim[5]));
    if ( imprim < -4 )  priopt(choix);
  }


  yams_printval();


  /* set adjacencies  */
  chrono(ON,&ctim[1]);
  chrono(ON,&ctim[2]);
  ret = tabvo2(sm,declic);
  chrono(OFF,&ctim[2]);
  if ( !ret ) {
    prierr(ERR,yerr.coderr);
    exit(1);
  }

  /* print surface quality */
  if ( imprim ) {
    if ( opts.ctrl & ISO ) 
      priqua(sm);
    else if ( sm->metric )
      priqua_a(sm);
    primsg(0000);
    if ( abs(imprim) > 1 ) {
      yerr.inderr[0] = sm->npmax;
      yerr.inderr[1] = sm->nemax;
      primsg(0002);
    }
  }

  /* pre-processing stage */
  yerr.inderr[0] = ++phase;
  out = 0;
  if ( abs(imprim) > 1 )  primsg(1000);
  chrono(ON,&ctim[2]);
  corr = sm->type & M_DETECT ? 1 : 0;
  if ( !setvoi(sm,corr) )   exit(1);
  if ( !ptmult(sm) )   exit(1);
  if ( absopt < 6 ) { 
    declic = 0.038;
		declic = opts.ctrl & ISO ? 1e-6 : 1.e-8;
    if ( !sident(sm,corr) )    exit(1);
    if ( !delnul(sm,declic) )  exit(1);
    if ( !optedg(sm) )         exit(1);
  }
  if ( sm->type & M_DETECT && !sident(sm,1) )  exit(1);

  /* smoothing */
  if ( absopt == 9 ) {
    if ( !noshrk(sm,opts.check) )  exit(1);
	  //if ( !hilbert(sm) )  exit(1);
    //if ( !denois(sm) )  exit(1);
  }
  else {
    if ( opts.iter < 0 )  opts.iter = 5;
    if ( absopt < 5 ) {
      if ( !norpoi(sm,0,corr) )  exit(1);
      if ( !tgepoi(sm,0,corr) )  exit(1);
    }
  }
  chrono(OFF,&ctim[2]);

  yerr.inderr[0] = phase;
  yerr.cooerr[0] = gttime(ctim[2]);
  if ( abs(imprim) > 1 ) {
    primsg(1001);
    if ( imprim < -4 ) {
      bilan(sm);
      prigap(sm);
    }
  }

  printf("absopt= %d\n", absopt);
  printf("imprim= %d\n", imprim);
  printf("sm->np %d\n", sm->np);
  printf("sm->dim %d\n", sm->dim);
  /* surface remeshing */
  yerr.inderr[0] = ++phase;
  if ( absopt && absopt <= 6 ) {
    if ( abs(imprim) > 1 )  primsg(1000);
    chrono(ON,&ctim[4]);

    /* geometry enrichment */
    if ( option > 0 ) {
      chrono(ON,&ctim[6]);
      if ( option == 4 )
        ret = yams4(sm);
      else if ( option == 6 )
        ret = yams6(sm);
       else
        ret = yams3(sm);
      chrono(OFF,&ctim[6]);
      if ( !ret )  exit(1);
    }

    /* surface simplification */
    if ( absopt == 1 )
      ret = yams1(sm);
    else if ( absopt == 2 ) {
      if ( opts.minnp < 0 )
        ret = yams2(sm);
      else
        ret = yams22(sm);
    }
    else if ( absopt == 5 && sm->type & M_METRIC ) 
      ret = calmet(sm);
    chrono(OFF,&ctim[4]);
    if ( !ret )  exit(1);

    yerr.inderr[0] = phase;
    yerr.cooerr[0] = gttime(ctim[4]);
    if ( abs(imprim) > 1 ) {
      primsg(1001);
      if ( imprim < -4 ) {
        if ( opts.ctrl & ISO )
          priqua(sm);
        else
          priqua_a(sm);
        prilen(sm);
      }
    }
  }

  /* mesh optimization */
  yerr.inderr[0] = ++phase;
  if ( absopt < 4 && absopt != 2 && yerr.coderr != 4000 ) {
    if ( abs(imprim) > 1 )  primsg(1000);

    chrono(ON,&ctim[3]);
    if ( sm->type & M_SMOOTH && yerr.coderr != 4000 ) {
      ret = optra4(sm,option);
      if ( !ret )  exit(1);
    }
    if ( absopt < 2 && opts.ffem && !optfem(sm) )  exit(1);
    chrono(OFF,&ctim[3]);
    yerr.inderr[0] = phase;
    yerr.cooerr[0] = gttime(ctim[3]);
    if ( abs(imprim) > 1 ) primsg(1001);
  }

  /* convert to quads (09-2003) */
  if ( sm->type & M_QUADS ) {
    yerr.inderr[0] = ++phase;
    if ( abs(imprim) > 1 )  primsg(1000);
    chrono(ON,&ctim[4]);

    if ( !yamsq(sm) )  exit(1);
 
    yerr.inderr[0] = phase;
    yerr.cooerr[0] = gttime(ctim[4]);    
    if ( abs(imprim) > 1 )  primsg(1001);
  }
  chrono(OFF,&ctim[1]);

  /* evaluation histograms */
  if ( abs(imprim) > 1 && absopt < 10 ) {
    if ( sm->type & M_QUADS )
      outqua_q(sm);
    else {
      if ( absopt == 1 )  prilen(sm);
      if ( opts.ctrl & ISO )
        outqua(sm);
      else {
	outqua_a(sm);
        outqua1_a(sm);
      }
      if ( sm->connex && info.cc > 1 )  rchsub(sm);
    }
  }
  if ( abs(imprim) > 1 )  primsg(0001);

  /* write resulting mesh */
    // a voir 
  if ( sm->type & M_OUTPUT ) {
    printf("freefem++:: outputfile yams\n");
    chrono(ON,&ctim[5]);
    out = yams8(sm,sm->outfile,absopt);
    chrono(OFF,&ctim[5]);
  }
  else {
    if ( imprim )  priout(sm);
    out=1;
  }

  yams_printval();

  /* print CPU requirements */
  chrono(OFF,&ctim[0]);
  if ( imprim ) {
    if ( imprim < 0 ) primem(sm->npmax);
    pritim(sm,option);
  }

  M_free(hash);

  ///* check for mem leaks */
  //if ( imprim < 0 && M_memLeak() )  M_memDump();

#ifdef DISTRIB
  /* free token */
  if ( !IsKeyCodeProtected(keycode) )
    free_token(&token);
#endif
  
  infondang = info.ndang;
  infocc = info.cc;
  
  return(0);
}

