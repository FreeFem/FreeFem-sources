// ORIG-DATE:     Fev 2010
// -*- Mode : c++ -*-
//
// SUMMARY  : liaison medit freefem++ : adaptmesh in 3d 
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
// AUTHOR   : Jacques Morice
// E-MAIL   : jacques.morice@ann.jussieu.fr
//
//   for automatic  compilation with ff-c++
//ff-c++-LIBRARY-dep:   yams.cpp
//ff-c++-cpp-dep: 
//  

/* 
 This file is part of Freefem++
 
 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.
 
 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

 Thank to the ARN ()  FF2A3 grant
 ref:ANR-07-CIS7-002-01 
 */

/*
./ff-c++ yams.cpp -I/Users/morice/work/postdoc/freefem++prod/src/libMesh/ -I/Users/morice/Desktop/adaptmesh3d/yams2.2010.02.22/sourcesnew -L/Users/morice/Desktop/adaptmesh3d/yams2.2010.02.22/objects/i386/ -lyams2 -L/Users/morice/work/postdoc/freefem++prod/src/libMesh/ -lMesh
*/

// ./ff-c++ yams.cpp -I../src/libMesh/ -I../download/include/yams/ -L../download/lib/yams/ -lyams2 -L/Users/morice/work/postdoc/freefem++prod/src/libMesh/ -lMesh

#include "ff++.hpp" 
#include "msh3.hpp"
//#define ADAPTLIBRARY
#include "memory.h"
#include "yamslib.h"
#include "eigenv.h"  // include dans libMesh

using namespace  Fem2D;
using namespace  yams;

// 3d mesh function

void mesh3_to_yams_pSurfMesh( const Mesh3 &Th3 , int memory, int choix, 
			      yams_pSurfMesh meshyams){

  /*
    Mesh3  :: maillage initiale
    memory :: memoire pour yams
    choix  :: option du remaillage
    ref    :: 
   */ 
  int k;
  int npinit,neinit;

  meshyams->dim = 3;
  meshyams->npfixe = Th3.nv;
  meshyams->nefixe = Th3.nbe;
  meshyams->ntet = Th3.nt;
  meshyams->nafixe  = 0; // Edges
  meshyams->nvfixe  = 0; // Normals 
  meshyams->ntfixe  = 0; // Tangents 
  npinit = meshyams->npfixe;
  neinit = meshyams->nefixe;
  // cette fonction change la taille des tableaux en fonctions des options : choix, memory, sm->type
  zaldy1( meshyams->nefixe, meshyams->npfixe, meshyams->nvfixe, memory, meshyams, choix);

  
  yams_pPoint ppt;
  for (k=1; k<=npinit; k++) {
    ppt = &meshyams->point[k];
    ppt->c[0] = Th3.vertices[k-1].x;
    ppt->c[1] = Th3.vertices[k-1].y;
    ppt->c[2] = Th3.vertices[k-1].z;
    ppt->ref  = Th3.vertices[k-1].lab & 0x7fff;
   
    ppt->tag  = M_UNUSED;
    ppt->color= 0;
    ppt->size = -1.;
    ppt->tge  = 0;
    ppt->geom = M_CURVE;
  }
  meshyams->npfixe = npinit;
  
  /* read mesh triangles */
  yams_pTriangle ptriangle;
  for (k=1; k<=neinit; k++) {
    const Triangle3 & K(Th3.be(k-1));
    ptriangle = &meshyams->tria[k];
    ptriangle->v[0] = Th3.operator()(K[0])+1;
    ptriangle->v[1] = Th3.operator()(K[1])+1;
    ptriangle->v[2] = Th3.operator()(K[2])+1;
    ptriangle->ref = K.lab& 0x7fff; 
  }

  /* tetrahedra */
  if( meshyams->ntet ){
    yams_pTetra ptetra;
    meshyams->tetra = (yams_Tetra*)calloc((meshyams->ntet+1)*sizeof(yams_Tetra));
    assert(meshyams->tetra);
    
    for (k=1; k<=meshyams->ntet; k++) {
      const Tet & K(Th3.elements[k-1]);
      ptetra = &meshyams->tetra[k];
      ptetra->v[0] = Th3.operator()(K[0])+1;
      ptetra->v[1] = Th3.operator()(K[1])+1;
      ptetra->v[2] = Th3.operator()(K[2])+1;
      ptetra->v[3] = Th3.operator()(K[3])+1;
      ptetra->ref = K.lab & 0x7fff;
    }
  }

  
  meshyams->ne = meshyams->nefixe;
  meshyams->np = meshyams->npfixe;
}

Mesh3 * yams_pSurfMesh_to_mesh3( yams_pSurfMesh sm, int infondang, int infocc ){

  /*
    Mesh3  :: maillage initiale
    memory :: memoire pour yams
    choix  :: option du remaillage
    ref    :: 
   */ 
  // variable a enlever par la suite
  yams_pGeomSupp    gs;
  yams_pGeomtge     gt;
  yams_pPoint       ppt;
  yams_pTriangle    pt1;
  yams_pTetra       ptt;
  yams_pEdge        pte;
  int i,k,np,ne,nn,nt,nav,natv,tatv,nbl;
  int          nedge,nridge,ndang,nrequis;
  int          is1,is2,ncorner,prequis;
   
  // freefempp variable
  int ff_nv, ff_nt, ff_nbe;
  
  
  /* mark connected component */
  ne = 0;
  for (k=1; k<=sm->npmax; k++) {
    ppt = & sm->point[k];
    ppt->tag |= M_UNUSED;
    ppt->flag = ppt->color = 0;
  }
  printf("sm->connex %d\n",sm->connex);
  if ( sm->connex > 0 ) {
    for (k=1; k<=sm->ne; k++) {
      pt1 = &sm->tria[k];
      if ( pt1->v[0] > 0 && pt1->cc == sm->connex ) {
        ne++;
        for (i=0; i<3; i++) {
          ppt = &sm->point[pt1->v[i]];
          ppt->tag &= ~M_UNUSED;
        }
      }
    }
  }
  else {
    /* mark used faces */
    for (k=1; k<=sm->ne; k++) {
      pt1 = &sm->tria[k];
      if ( !pt1->v[0] )  continue;
      ++ne;
      for (i=0; i<3; i++) {
        ppt = &sm->point[pt1->v[i]];
        ppt->tag &= ~M_UNUSED;
      }
    }
  }

  if ( sm->ntet ) {
    for (k=1; k<=sm->ntet; k++) {
      ptt = &sm->tetra[k];
      if ( !ptt->v[0] )  continue;
      for (i=0; i<4; i++) {
        ppt = &sm->point[ptt->v[i]];
        ppt->tag &= ~M_UNUSED;
      }
    }
  }

  /* mark used vertices */
  np = nav = 0;
  ncorner = prequis = 0;
  for (k=1; k<=sm->npmax; k++) {
    ppt = &sm->point[k];
    if ( ppt->tag & M_UNUSED )  continue;
    ppt->tmp = ++np;
    if ( ppt->tag == M_NOTAG )  nav++;
  }

  ff_nv = np;   // number of vertex
  //
  Vertex3 *ff_v = new Vertex3[ff_nv];
  int kk=0;
  for(k=1; k<=sm->npmax; k++) {
    ppt = &sm->point[k];
    if ( ppt->tag & M_UNUSED )  continue;
    ff_v[kk].x = ppt->c[0];
    ff_v[kk].y = ppt->c[1];
    ff_v[kk].z = ppt->c[2];
    ff_v[kk].lab = ppt->ref;
    kk++;
    if (ppt->tag & M_CORNER)    ncorner++;
    if (ppt->tag & M_REQUIRED ) prequis++;    
  }
  assert(kk==ff_nv);
  // write triangle
  nedge  = sm->dim == 3 ? infondang : 0;
  nridge = nrequis = nn = nt = natv = tatv = 0;
  
  for (k=1; k<=sm->ne; k++) {
    pt1  = &sm->tria[k];
    if ( !pt1->v[0] )  continue;
    else if ( sm->connex > 0 && pt1->cc != sm->connex ) continue;
    nt++;
  }
  
  ff_nbe = nt;
  Triangle3 *ff_b = new Triangle3[ff_nbe];
  Triangle3 *ff_bb = ff_b;


  for (k=1; k<=sm->ne; k++) {
    int iv[3],lab;
    pt1  = &sm->tria[k];
    if ( !pt1->v[0] )  continue;
    else if ( sm->connex > 0 && pt1->cc != sm->connex ) continue;
    iv[0] = sm->point[pt1->v[0]].tmp-1;
    iv[1] = sm->point[pt1->v[1]].tmp-1;
    iv[2] = sm->point[pt1->v[2]].tmp-1;
    lab = (int)(sm->connex < 0 ? pt1->cc : pt1->ref);
 
    (*ff_bb++).set( ff_v, iv, lab);
    
    for (i=0; i<3; i++) {
      ppt = &sm->point[pt1->v[i]];
      gs  = &sm->geom[pt1->vn[i]];
      gt  = &sm->tgte[ppt->tge];
      if ( ppt->tag > M_NOTAG ) {
        natv++;
        if ( ppt->tag & M_CORNER )  tatv++;
      }
      if ( !gs->newnum )  gs->newnum = ++nn;
      if ( !gt->newnum )  gt->newnum = ++nt;
      if ( !pt1->edg[i] && pt1->tag[i] == M_NOTAG )  continue;
      else if ( pt1->adj[i] && (k > pt1->adj[i]) )   continue;
      nedge++;
      if ( pt1->tag[i] & M_RIDGE_GEO )  nridge++;
      if ( pt1->tag[i] & M_REQUIRED )   nrequis++;
    } 
    
  }

  // les autres avoir par la suite
  cout << " nv " << ff_nv << " nbe" << ff_nbe << endl;
  Mesh3 *TH3_T = new Mesh3(ff_nv,ff_nbe,ff_v,ff_b);
  return TH3_T;
}

void solyams_pSurfMesh( yams_pSurfMesh sm, const int &type, const KN<double> & tabsol, float hmin, float hmax){
  yams_pPoint ppt;
  yams_pMetric  pm;
  int i,k;
  double   sizeh,m[6],lambda[3],vp[2][2],vp3[3][3];
  hmin =  FLT_MAX;
  hmax = -FLT_MAX;
  

  if(type == 1){
    for (k=1; k<=sm->npfixe; k++) {
      ppt = &sm->point[k];
      ppt->size = (float) tabsol[k];
      hmin = min(ppt->size,hmin);
      hmax = max(ppt->size,hmax);
    }
  }
  else if( type == 3 ){
    if( !sm->metric ){
      cerr << " we give metric solution bug" << endl;
    }
    
    for (k=1; k<=sm->npfixe; k++) {
      ppt = &sm->point[k];    
      pm  = &sm->metric[k];
      memset(pm->m,6*sizeof(float),0.);
      for (i=0; i<6; i++)
	m[i] = (float) tabsol[(k-1)*6+i];
      
      pm->m[0] = m[0];
      pm->m[1] = m[1];
      pm->m[2] = m[3];
      pm->m[3] = m[2];
      pm->m[4] = m[4];
      pm->m[5] = m[5];
      pm->k1   = pm->k2 = (float)FLT_MAX;
      for (i=0; i<6; i++)  m[i] = pm->m[i];
      if ( !eigenv(1,m,lambda,vp3) ) {
	fprintf(stderr,"  ## ERR 9201, inbbf, Not a metric tensor. Discarded\n");
	free(sm->metric);
	sm->metric = 0;
	exit(1);
      }
      sizeh     = max(max(lambda[0],lambda[1]),lambda[2]);
      ppt->size = max(1.0 / sqrt(sizeh),EPS);
      hmin = min(ppt->size,hmin);
      hmax = max(ppt->size,hmax);
    }
  }
     
}

void yams_inival(int intopt[23],double fopt[14]){
/*
    intopt : 0  !! anisotropie
             1  !! ecp  // enl
             2  !! extended out put file 
	     3  !! FE correction 
	     4  !! Formatted (ascii) output file
	     5  !! save metric file
	     6  !! msh2  // enl
	     7  !! Split multiple connected points
	     8  !! memory
	     9  !! connected component
	    10  !! vrml  //enl
	    11  !! imprim  
	    12  !! nm : Create point on straight edge (no mapping)
	    13  !! nc : No validity check during smoothing (opt. 9)
	    14  !! np : Specify number of points desired
	    15  !! nit : Nb Iter
	    16  !! nq  : Output quads // enl
	    17  !! nr  : No ridge detection
	    18  !! ns  : No point smoothing
	    19  !! no  : No output file  // enl
	    20  !! ref : Ignore face references // enl
	    // rajouter lors de l'ouverture du fichiers yams
	    21  !! absolute : opts.ctrl &= ~REL;
	    22  !! set optim option

    fopt   : 0  !! iso 
             1  !! eps 
	     3  !! opts.lambda
	     4  !! opts.mu
	     5  !! set optim option
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

/* set default values for options */
  // fopt 5,
  fopt[7]   = -2.0;    
  fopt[8]   = -2.0;
  fopt[6]   =   1.3;       /* default mesh gradation     */
  fopt[1]   =   0.01;      /* geometric approximation    */
  fopt[0]   =   0.0;
  fopt[11]  =   1.0 / BETAC;
  fopt[3]   =   -1.0;
  fopt[4]   =   -1.0;
  fopt[13]  =   45.;  // default RIDG = 45.
  //opts.ridge  =   cos(RIDG*M_PI/180.);
  //opts.geom   =   cos(GEOM*M_PI/180.);
  fopt[12] =   COS45DEG;  /* Walton limitation          */
  fopt[9]  =   -2;       /* default = 1 unit           */
  fopt[10] =   QUALCOE;   /* quality degradation        */
  //opts.ctrl   =   REL | ISO;  initialisation by default
  
  

  // intopt :: 3,7,13,14,15,20
  intopt[15]  = -1;
  intopt[13]  = 0;
  intopt[14]  = -1;
  
  /* get decimation parameters */
  intopt[20] = 0;
  intopt[3]  = 0;
  intopt[7]  = 1;
  intopt[22] = 0;

  // demander P. Frey
  intopt[0] = 0;
  intopt[1] = 0; 
  intopt[2] = 0;
    
  intopt[4] = 0;
  intopt[5] = 0;
  intopt[6] = 0;

  intopt[8] = -1;
  intopt[9] =  0; // par default   // a verifier
  intopt[10] = 0;
  intopt[11] = -99;
  intopt[12] = 0;   // par default

  intopt[16] = 0;
  intopt[17] = 0;
  intopt[18] = 0;
  intopt[19] = 0;

  intopt[21] = 0;
  
}

class yams_Op: public E_F0mps 
{
public:
  typedef pmesh3  Result;
  Expression eTh;
  int nbsol;
  int nbsolsize;
  int type;
  int dim;
  vector<Expression> sol;

  static const int n_name_param = 3; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  KN_<int>  arg(int i,Stack stack,KN_<int> a ) const
  { return nargs[i] ? GetAny<KN_<int> >( (*nargs[i])(stack) ): a;}
  KN_<double>  arg(int i,Stack stack,KN_<double> a ) const
  { return nargs[i] ? GetAny<KN_<double> >( (*nargs[i])(stack) ): a;}
  double  arg(int i,Stack stack,double a ) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  int  arg(int i,Stack stack, int a ) const{ return nargs[i] ? GetAny< int >( (*nargs[i])(stack) ): a;}
  
  
public:
  yams_Op(const basicAC_F0 &  args) : sol( args.size()-1 )
  {
    
    cout << "yams"<< endl;
    args.SetNameParam(n_name_param,name_param,nargs);
    eTh=to<pmesh3>(args[0]); 
    dim=3;
    nbsol = args.size()-1;
    if(nbsol >1) 
      CompileError(" yams accept only one solution ");
    int ksol=0; 
      
    if(nbsol == 1){
      int i=1;
      if (args[i].left()==atype<E_Array>())
	{
	  const E_Array * a = dynamic_cast<const E_Array *>(args[i].LeftValue());
	  ffassert(a);
	  ksol+=a->size(); 
	}
      else
	ksol++;
      sol.resize(ksol); 
      
      // type :: 1 sca, 2 vector, 3 symtensor
      
      ksol=0; 
      nbsolsize=0;
      type = 0;
          
      if (args[i].left()==atype<E_Array>())
	{
	  const E_Array * a = dynamic_cast<const E_Array *>(args[i].LeftValue());
	  ffassert(a);
	  int N=a->size();
	  nbsolsize=nbsolsize+N;
	  switch (N){
	    /*
	      case 3 :
		type[i-1]=2; 
		for (int j=0;j<N;j++)             
		sol[ksol++]=to<double>((*a)[j]);
		break;
	    */
	  case 6 :
	    type=3;   
	    for (int j=0;j<N;j++)             
	      sol[ksol++]=to<double>((*a)[j]); 
	    break;
	  default :
	    CompileError(" 3D solution for yams is a scalar (1 comp) or a symetric tensor (6 comp)");
	    break;
	  }
	}
      else 
	{
	  type=1;
	  nbsolsize=nbsolsize+1;
	  sol[ksol++]=to<double>(args[i]);
	} 

      if( nargs[2]  ) 
	CompileError(" we give two metric for yams ");
    }

  }
    
  static ArrayOfaType  typeargs() { return  ArrayOfaType( atype< pmesh3 >(), true); }// all type
  static  E_F0 * f(const basicAC_F0 & args) { return new yams_Op(args);} 
  AnyType operator()(Stack stack)  const ;
  operator aType () const { return atype< pmesh3 >();} 
};


basicAC_F0::name_and_type  yams_Op::name_param[]= {
  {  "loptions", &typeid(KN_<long>)},
  {  "doptions", &typeid(KN_<double>)},
  {  "metric", &typeid(KN_<double>)}
};

AnyType yams_Op::operator()(Stack stack)  const 
{
  // initialisation
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh3 * pTh= GetAny<Mesh3 *>((*eTh)(stack));
  ffassert( pTh );
  Mesh3 &Th3=*pTh;
  int nv=Th3.nv;
  int nt=Th3.nt;
  int nbe=Th3.nbe;

  KN<int> defaultintopt(23);
  KN<double> defaultfopt(14);
  defaultintopt = 0;
  defaultfopt   = 0.;
  yams_inival( defaultintopt, defaultfopt);

  KN<int> intopt(arg(0,stack,defaultintopt));
  assert( intopt.N() == 23 );
  KN<double> fopt(arg(1,stack,defaultfopt));
  assert( fopt.N() == 14 );
  KN<double> metric;

  int mtype=type;
  if( nargs[2]  ){ 
    metric = GetAny<KN_<double> >( (*nargs[2])(stack) );
    if(metric.N()==Th3.nv){
      mtype=1;
      intopt[1]=0;
    }
    else if(metric.N()==6*Th3.nv){ 
      intopt[1]=1; 
      mtype=3;
    }
    else 
      cerr << "sizeof vector metric is incorrect, size will be Th.nv or 6*Th.nv" << endl;
  }
  else if(nbsol>0){
    if( type == 1 ){
      intopt[1]=0;
      metric.resize(Th3.nv);
      metric=0.;
    }
    else if( type ==3 ){
      intopt[1]=1;
      metric.resize(6*Th3.nv);
      metric=0.;
    }
  }
  else{
    if( intopt[1]==0 ){ metric.resize(Th3.nv); metric=0.;}
    else if ( intopt[1]==1 ){ metric.resize(6*Th3.nv); metric=0.;}
  }
  // mesh for yams
  yams_pSurfMesh yamsmesh;
  yamsmesh = (yams_pSurfMesh)calloc(1,sizeof(yams_SurfMesh));
  if ( !yamsmesh ){
    cerr << "allocation error for SurfMesh for yams" << endl;
  }
  yamsmesh->infile  = NULL;
  yamsmesh->outfile = NULL;
  yamsmesh->type    = M_SMOOTH | M_QUERY | M_DETECT | M_BINARY | M_OUTPUT;
 

  mesh3_to_yams_pSurfMesh( Th3 , intopt[8], intopt[22], yamsmesh);
    
  
  // solution for yams2
  if(nbsol){
    MeshPoint *mp3(MeshPointStack(stack)); 
    
    KN<bool> takemesh(nv);
    takemesh=false;
    for(int it=0;it<nt;it++){
      for(int iv=0;iv<4;iv++){
	int i=Th3(it,iv);
	
	if(takemesh[i]==false){
	  mp3->setP(&Th3,it,iv);

	  for(int ii=0;ii<nbsolsize;ii++){
	    metric[i*nbsolsize+ii] = GetAny< double >( (*sol[ii])(stack) );
	  }
	  takemesh[i] = true; 
	}
      }
    }
  }
  if( nargs[2] || (nbsol > 0) ){ 
    float hmin,hmax;
    solyams_pSurfMesh( yamsmesh, mtype, metric, hmin, hmax);
    yamsmesh->nmfixe = yamsmesh->npfixe;
    fopt[7]=hmin;
    fopt[8]=hmax;
  }
  else{
    yamsmesh->nmfixe = 0;
  }
  int infondang=0, infocc=0;
  int res = yams_main(yamsmesh, intopt, fopt, infondang, infocc );

  cout << " yamsmesh->dim " << yamsmesh->dim << endl;
  if( res > 0){
    cout << " problem with yams :: error " <<  res << endl; 
    exit(1);
  }
  
  Mesh3 *Th3_T = yams_pSurfMesh_to_mesh3( yamsmesh, infondang, infocc );
  
  // recuperer la solution ????
  cout << &yamsmesh->point << " " << &yamsmesh->tria << " "  <<&yamsmesh->geom << " "  << &yamsmesh->tgte << endl;
  cout << &yamsmesh << endl;
  M_free(yamsmesh->point);
  M_free(yamsmesh->tria);
  M_free(yamsmesh->geom);
  M_free(yamsmesh->tgte);
  if ( yamsmesh->metric ) M_free(yamsmesh->metric);
  if ( yamsmesh->edge )   M_free(yamsmesh->edge);
  if ( yamsmesh->tetra )   M_free(yamsmesh->tetra);
  M_free(yamsmesh);

  /* check for mem leaks */
  if ( M_memLeak() )  M_memDump();
  *mp=mps;
  return SetAny<pmesh3>(Th3_T);
}



class Init1 { public:
  Init1();
};

static Init1 init1;  //  une variable globale qui serat construite  au chargement dynamique 

Init1::Init1(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 
  //typedef Mesh3 *pmesh3;
  if(verbosity) cout << " load: yams  " << endl;
  
  Global.Add("yams","(",new OneOperatorCode<yams_Op>);
 
}


#define  WITH_NO_INIT
#include "msh3.hpp" 

