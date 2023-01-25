//ff-c++-LIBRARY-dep: mmg scotch
//ff-c++-cpp-dep:

#include "ff++.hpp"
#include "memory.h"
#include "mmg/libmmg.h"
#include "GenericMesh.hpp"

using namespace Fem2D;

template<class ffmesh> int ffmesh_to_MMG5_pMesh(const ffmesh &, MMG5_pMesh&);

template<>
int ffmesh_to_MMG5_pMesh<Mesh3>(const Mesh3 &Th, MMG5_pMesh& mesh) {
  int nVertices       = Th.nv;
  int nTetrahedra     = Th.nt;
  int nPrisms         = 0;
  int nTriangles      = Th.nbe;
  int nQuadrilaterals = 0;
  int nEdges          = 0;

    if ( MMG3D_Set_meshSize(mesh,nVertices,nTetrahedra,nPrisms,nTriangles,
                           nQuadrilaterals,nEdges) != 1 ) {
      exit(EXIT_FAILURE);
    }

    for (int k = 0; k < Th.nv; k++) {
      if ( MMG3D_Set_vertex(mesh,Th.vertices[k].x,Th.vertices[k].y,
                           Th.vertices[k].z, Th.vertices[k].lab, k+1) != 1 ) { 
        exit(EXIT_FAILURE);
      }
    }

    for (int k = 0; k < Th.nt; k++) {
      const Tet &K(Th.elements[k]);
      if ( MMG3D_Set_tetrahedron(mesh,Th.operator()(K[0])+1,Th.operator()(K[1])+1,
                                Th.operator()(K[2])+1,Th.operator()(K[3])+1,K.lab,k+1) != 1 ) {
        exit(EXIT_FAILURE);
      }
    }

    for (int k = 0; k < Th.nbe; k++) {
      const Triangle3 &K(Th.be(k));
      if ( MMG3D_Set_triangle(mesh,Th.operator()(K[0])+1,Th.operator()(K[1])+1,Th.operator()(K[2])+1,
                             K.lab,k+1) != 1 ) {
        exit(EXIT_FAILURE);
      }
    }

  return 0;
}

template<>
int ffmesh_to_MMG5_pMesh<MeshS>(const MeshS &Th, MMG5_pMesh& mesh) {
  int nVertices       = Th.nv;
  int nTetrahedra     = 0;
  int nPrisms         = 0;
  int nTriangles      = Th.nt;
  int nQuadrilaterals = 0;
  int nEdges          = Th.nbe;

    if ( MMGS_Set_meshSize(mesh,nVertices,nTriangles,nEdges) != 1 ) {
      exit(EXIT_FAILURE);
    }

    for (int k = 0; k < Th.nv; k++) {
      if ( MMGS_Set_vertex(mesh,Th.vertices[k].x,Th.vertices[k].y,
                           Th.vertices[k].z, Th.vertices[k].lab, k+1) != 1 ) { 
        exit(EXIT_FAILURE);
      }
    }

    for (int k = 0; k < Th.nt; k++) {
      const TriangleS &K(Th.elements[k]);
      if ( MMGS_Set_triangle(mesh,Th.operator()(K[0])+1,Th.operator()(K[1])+1,Th.operator()(K[2])+1,
                             K.lab,k+1) != 1 ) {
        exit(EXIT_FAILURE);
      }
    }

    for (int k = 0; k < Th.nbe; k++) {
      const BoundaryEdgeS &K(Th.be(k));
      if ( MMG3D_Set_edge(mesh,Th.operator()(K[0])+1,Th.operator()(K[1])+1,
                             K.lab,k+1) != 1 ) {
        exit(EXIT_FAILURE);
      }
    }

  return 0;
}

template<class ffmesh> int MMG5_pMesh_to_ffmesh(const MMG5_pMesh&, ffmesh *&);

template<>
int MMG5_pMesh_to_ffmesh<Mesh3>(const MMG5_pMesh& mesh, Mesh3 *&T_TH3) {
    int ier;

    int nVertices   = 0;
    int nTetrahedra = 0;
    int nTriangles  = 0;
    int nEdges      = 0;

    if ( MMG3D_Get_meshSize(mesh,&nVertices,&nTetrahedra,NULL,&nTriangles,NULL,
                       &nEdges) !=1 ) { 
      ier = MMG5_STRONGFAILURE;
    }

    Vertex3 *v = new Vertex3[nVertices];
    Tet *t = new Tet[nTetrahedra];
    Tet *tt = t;
    Triangle3 *b = new Triangle3[nTriangles];
    Triangle3 *bb = b;
    int k;

    int corner, required;

    for (k = 0; k < nVertices; k++) {
      if ( MMG3D_Get_vertex(mesh,&(v[k].x),&(v[k].y),&(v[k].z),
                                   &(v[k].lab),&(corner),&(required)) != 1 ) {
        cout << "Unable to get mesh vertex " << k << endl;
        ier = MMG5_STRONGFAILURE;
      }
    }

    for ( k=0; k<nTetrahedra; k++ ) {
      int iv[4], lab;
      if ( MMG3D_Get_tetrahedron(mesh,
                               &(iv[0]),&(iv[1]),
                               &(iv[2]),&(iv[3]),
                               &(lab),&(required)) != 1 ) {
        cout << "Unable to get mesh tetra " << k << endl;
        ier = MMG5_STRONGFAILURE;
      }
      for (int i=0; i<4; i++)
        iv[i]--;
      tt++->set(v, iv, lab);
    }

    for ( k=0; k<nTriangles; k++ ) {
      int iv[3], lab;
      if ( MMG3D_Get_triangle(mesh,
                               &(iv[0]),&(iv[1]),&(iv[2]),
                               &(lab),&(required)) != 1 ) {
        cout << "Unable to get mesh triangle " << k << endl;
        ier = MMG5_STRONGFAILURE;
      }
      for (int i=0; i<3; i++)
        iv[i]--;
      bb++->set(v, iv, lab);
    }

    T_TH3 = new Mesh3(nVertices, nTetrahedra, nTriangles, v, t, b, 1);

    if (verbosity > 1) {
      cout << "transformation maillage --> mesh3 " << endl;
      cout << "vertices =" << nVertices << endl;
      cout << "tetrahedrons =" << nTetrahedra << endl;
      cout << "triangles =" << nTriangles << endl;
      cout << "T_TH3" << T_TH3->nv << " " << T_TH3->nt << " " << T_TH3->nbe << endl;
    }

  return 0;
}

template<>
int MMG5_pMesh_to_ffmesh<MeshS>(const MMG5_pMesh& mesh, MeshS *&T_TH3) {
    int ier;

    int nVertices   = 0;
    int nTriangles  = 0;
    int nEdges      = 0;

    if ( MMGS_Get_meshSize(mesh,&nVertices,&nTriangles,&nEdges) !=1 ) { 
      ier = MMG5_STRONGFAILURE;
    }

    Vertex3 *v = new Vertex3[nVertices];
    TriangleS *t = new TriangleS[nTriangles];
    TriangleS *tt = t;
    BoundaryEdgeS *b = new BoundaryEdgeS[nEdges];
    BoundaryEdgeS *bb = b;
    int k;

    int corner, required, ridge;

    for (k = 0; k < nVertices; k++) {
      if ( MMGS_Get_vertex(mesh,&(v[k].x),&(v[k].y),&(v[k].z),
                                   &(v[k].lab),&(corner),&(required)) != 1 ) {
        cout << "Unable to get mesh vertex " << k << endl;
        ier = MMG5_STRONGFAILURE;
      }
    }

    for ( k=0; k<nTriangles; k++ ) {
      int iv[3], lab;
      if ( MMGS_Get_triangle(mesh,
                               &(iv[0]),&(iv[1]),&(iv[2]),
                               &(lab),&(required)) != 1 ) {
        cout << "Unable to get mesh triangle " << k << endl;
        ier = MMG5_STRONGFAILURE;
      }
      for (int i=0; i<3; i++)
        iv[i]--;
      tt++->set(v, iv, lab);
    }

    for ( k=0; k<nEdges; k++ ) {
      int iv[2], lab;
      if ( MMGS_Get_edge(mesh,
                               &(iv[0]),&(iv[1]),
                               &(lab), &(ridge), &(required)) != 1 ) {
        cout << "Unable to get mesh edge " << k << endl;
        ier = MMG5_STRONGFAILURE;
      }
      for (int i=0; i<2; i++)
        iv[i]--;
      bb++->set(v, iv, lab);
    }

    T_TH3 = new MeshS(nVertices, nTriangles, nEdges, v, t, b, 1);

    if (verbosity > 1) {
      cout << "transformation maillage --> meshS " << endl;
      cout << "vertices =" << nVertices << endl;
      cout << "triangles =" << nTriangles << endl;
      cout << "edges =" << nEdges << endl;
      cout << "T_TH3" << T_TH3->nv << " " << T_TH3->nt << " " << T_TH3->nbe << endl;
    }

  return 0;
}

template<class ffmesh>
class mmg_Op : public E_F0mps {
 public:
  Expression eTh;
  static const int n_name_param = std::is_same<ffmesh,Mesh3>::value ? 27 : 20;
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];

  KN_< long > arg(int i, Stack stack, KN_< long > a) const {
    return nargs[i] ? GetAny< KN_< long > >((*nargs[i])(stack)) : a;
  }

  KN_< double > arg(int i, Stack stack, KN_< double > a) const {
    return nargs[i] ? GetAny< KN_< double > >((*nargs[i])(stack)) : a;
  }

  double arg(int i, Stack stack, double a) const {
    return nargs[i] ? GetAny< double >((*nargs[i])(stack)) : a;
  }

  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

  long arg(int i, Stack stack, bool a) const {
    return nargs[i] ? GetAny< bool >((*nargs[i])(stack)) : a;
  }

 public:
  mmg_Op(const basicAC_F0 &args, Expression tth) : eTh(tth) {
    if (verbosity > 1) {
      cout << "mmg" << endl;
    }

    args.SetNameParam(n_name_param, name_param, nargs);

  }

  AnyType operator( )(Stack stack) const;
};

template<>
basicAC_F0::name_and_type mmg_Op<Mesh3>::name_param[] = {
{"metric"            , &typeid(KN< double > *)},
{"verbose"           , &typeid(long)},/*!< [-1..10], Tune level of verbosity */
{"mem"               , &typeid(long)},/*!< [n/-1], Set memory size to n Mbytes or keep the default value */
{"debug"             , &typeid(bool)},/*!< [1/0], Turn on/off debug mode */
{"angle"             , &typeid(bool)},/*!< [1/0], Turn on/off angle detection */
{"iso"               , &typeid(bool)},/*!< [1/0], Level-set meshing */
{"nofem"             , &typeid(bool)},/*!< [1/0], Generate a non finite element mesh */
{"opnbdy"            , &typeid(bool)},/*!< [1/0], Preserve triangles at interface of 2 domains with same reference */
{"lag"               , &typeid(long)},/*!< [-1/0/1/2], Lagrangian option */
{"optim"             , &typeid(bool)},/*!< [1/0], Optimize mesh keeping its initial edge sizes */
{"optimLES"          , &typeid(bool)},/*!< [1/0], Strong mesh optimization for Les computations */
{"noinsert"          , &typeid(bool)},/*!< [1/0], Avoid/allow point insertion */
{"noswap"            , &typeid(bool)},/*!< [1/0], Avoid/allow edge or face flipping */
{"nomove"            , &typeid(bool)},/*!< [1/0], Avoid/allow point relocation */
{"nosurf"            , &typeid(bool)},/*!< [1/0], Avoid/allow surface modifications */
//{"nreg"              , &typeid(bool)},/*!< [0/1], Enable normal regularization */
{"renum"             , &typeid(bool)},/*!< [1/0], Turn on/off point relocation with Scotch */
{"anisosize"         , &typeid(bool)},/*!< [1/0], Turn on/off anisotropic metric creation when no metric is provided */
{"octree"            , &typeid(long)},/*!< [n], Specify the max number of points per PROctree cell (DELAUNAY) */
{"angleDetection"    , &typeid(double)},/*!< [val], Value for angle detection */
{"hmin"              , &typeid(double)},/*!< [val], Minimal mesh size */
{"hmax"              , &typeid(double)},/*!< [val], Maximal mesh size */
{"hsiz"              , &typeid(double)},/*!< [val], Constant mesh size */
{"hausd"             , &typeid(double)},/*!< [val], Control global Hausdorff distance (on all the boundary surfaces of the mesh) */
{"hgrad"             , &typeid(double)},/*!< [val], Control gradation */
{"ls"                , &typeid(double)},/*!< [val], Value of level-set */
//{"rmc"               , &typeid(double)},/*!< [-1/val], Remove small connex componants in level-set mode */
{"requiredTriangle"  , &typeid(KN<long>*)},/*!< [val], References of surfaces with required triangles */
{"localParameter"    , &typeid(KNM<double>*)}/*!< [val], Local parameters on given surfaces */
};

template<>
basicAC_F0::name_and_type mmg_Op<MeshS>::name_param[] = {
{"metric"            , &typeid(KN< double > *)},
{"verbose"           , &typeid(long)},/*!< [-1..10], Tune level of verbosity */
{"mem"               , &typeid(long)},/*!< [n/-1], Set memory size to n Mbytes or keep the default value */
{"debug"             , &typeid(bool)},/*!< [1/0], Turn on/off debug mode */
{"angle"             , &typeid(bool)},/*!< [1/0], Turn on/off angle detection */
{"iso"               , &typeid(bool)},/*!< [1/0], Level-set meshing */
{"keepRef"           , &typeid(bool)},/*!< [1/0], Preserve the initial domain references in level-set mode */
//{"optim"             , &typeid(bool)},/*!< [1/0], Optimize mesh keeping its initial edge sizes */
{"noinsert"          , &typeid(bool)},/*!< [1/0], Avoid/allow point insertion */
{"noswap"            , &typeid(bool)},/*!< [1/0], Avoid/allow edge or face flipping */
{"nomove"            , &typeid(bool)},/*!< [1/0], Avoid/allow point relocation */
{"nreg"              , &typeid(bool)},/*!< [0/1], Disabled/enabled normal regularization */
{"renum"             , &typeid(bool)},/*!< [1/0], Turn on/off point relocation with Scotch */
{"angleDetection"    , &typeid(double)},/*!< [val], Value for angle detection */
{"hmin"              , &typeid(double)},/*!< [val], Minimal mesh size */
{"hmax"              , &typeid(double)},/*!< [val], Maximal mesh size */
{"hsiz"              , &typeid(double)},/*!< [val], Constant mesh size */
{"hausd"             , &typeid(double)},/*!< [val], Control global Hausdorff distance (on all the boundary surfaces of the mesh) */
{"hgrad"             , &typeid(double)},/*!< [val], Control gradation */
{"ls"                , &typeid(double)},/*!< [val], Value of level-set */
{"requiredEdge"      , &typeid(KN<long>*)}/*!< [val], References of boundaries with required edges */
};

template<class ffmesh>
class mmg_ff : public OneOperator {
 public:
    mmg_ff( ) : OneOperator(atype< const ffmesh* >( ), atype< const ffmesh* >( )) {pref=10;}
    // to remove ambiguity with mmg3-v4

  E_F0 *code(const basicAC_F0 &args) const { return new mmg_Op<ffmesh>(args, t[0]->CastTo(args[0])); }
};

template<>
AnyType mmg_Op<Mesh3>::operator( )(Stack stack) const {
  // initialisation
  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
  Mesh3 *pTh = GetAny< Mesh3 * >((*eTh)(stack));

  ffassert(pTh);
  Mesh3 &Th = *pTh;
  int nv = Th.nv;
  int nt = Th.nt;
  int nbe = Th.nbe;

  KN< double > *pmetric = 0;
  KN< long > *prequiredTriangle = 0;
  KNM< double > *plocalParameter = 0;

  if (nargs[0]) {
    pmetric = GetAny< KN< double > * >((*nargs[0])(stack));
  }
  if (nargs[25]) {
    prequiredTriangle = GetAny< KN< long > * >((*nargs[25])(stack));
  }
  if (nargs[26]) {
    plocalParameter = GetAny< KNM< double > * >((*nargs[26])(stack));
  }

  MMG5_pMesh mesh;
  MMG5_pSol sol;
  MMG5_pSol met;

  mesh = nullptr;
  sol = nullptr;
  met = nullptr;

  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&sol,
                  MMG5_ARG_end);

  ffmesh_to_MMG5_pMesh(Th, mesh);
  
    if (pmetric && pmetric->N( ) > 0) {
      const KN< double > &metric = *pmetric;
      if (metric.N( ) == Th.nv) {
        if ( MMG3D_Set_solSize(mesh,sol,MMG5_Vertex,Th.nv,MMG5_Scalar) != 1 ) { 
          printf("Unable to allocate the metric array.\n");
          exit(EXIT_FAILURE);
        }
        if ( MMG3D_Set_scalarSols(sol,metric) != 1 ) { 
          printf("Unable to set metric.\n");
          exit(EXIT_FAILURE);
        }
      }
      else {
        if ( MMG3D_Set_solSize(mesh,sol,MMG5_Vertex,Th.nv,MMG5_Tensor) != 1 ) { 
          printf("Unable to allocate the metric array.\n");
          exit(EXIT_FAILURE);
        }
        static const int perm[6] = {0, 1, 3, 2, 4, 5};
        for (int k=0; k<Th.nv; k++) {
          if ( MMG3D_Set_tensorSol(sol, metric[6*k+perm[0]], metric[6*k+perm[1]], metric[6*k+perm[2]], 
                                  metric[6*k+perm[3]], metric[6*k+perm[4]], metric[6*k+perm[5]], k+1) != 1 ) { 
            printf("Unable to set metric.\n");
            exit(EXIT_FAILURE);
          }
        }
      }
    }
    if (prequiredTriangle && prequiredTriangle->N( ) > 0) {
      const KN< long > &requiredTriangle = *prequiredTriangle;
      std::sort(requiredTriangle + 0, requiredTriangle + requiredTriangle.N());
      int nt;
      if ( MMG3D_Get_meshSize(mesh,NULL,NULL,NULL,&nt,NULL,NULL) !=1 ) {
        exit(EXIT_FAILURE);
      }
      for (int k=1; k<=nt; k++) {
        int ref, dummy;
        if ( MMG3D_Get_triangle(mesh,&dummy,&dummy,&dummy,
                    &ref,NULL) != 1 ) {
          exit(EXIT_FAILURE);
        }
        if (std::binary_search(requiredTriangle + 0, requiredTriangle + requiredTriangle.N(), ref)) {
          if ( MMG3D_Set_requiredTriangle(mesh,k) != 1 ) {
            exit(EXIT_FAILURE);
          }
        }
      }
    }
    if (plocalParameter && plocalParameter->M( ) > 0) {
      const KNM< double > &localParameter = *plocalParameter;
      ffassert(localParameter.N() == 4);
      if ( MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_numberOfLocalParam,localParameter.M()) != 1 ) {
        exit(EXIT_FAILURE);
      }
      for(int j = 0; j < localParameter.M(); ++j) {
        if ( MMG3D_Set_localParameter(mesh,sol,MMG5_Triangle,localParameter(0,j),localParameter(1,j),localParameter(2,j),localParameter(3,j)) != 1 ) {
          exit(EXIT_FAILURE);
        }
      }
    }

  int i=1;
  if (nargs[i]) MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_verbose,       arg(i,stack,0L));    i++;   /*!< [-1..10], Tune level of verbosity */
  if (nargs[i]) MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_mem,           arg(i,stack,0L));    i++;   /*!< [n/-1], Set memory size to n Mbytes or keep the default value */
  if (nargs[i]) MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_debug,         arg(i,stack,false)); i++;   /*!< [1/0], Turn on/off debug mode */
  if (nargs[i]) MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_angle,         arg(i,stack,false)); i++;   /*!< [1/0], Turn on/off angle detection */
  if (nargs[i]) MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_iso,           arg(i,stack,false)); i++;   /*!< [1/0], Level-set meshing */
  if (nargs[i]) MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_nofem,         arg(i,stack,false)); i++;   /*!< [1/0], Generate a non finite element mesh */
  if (nargs[i]) MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_opnbdy,        arg(i,stack,false)); i++;   /*!< [1/0], Preserve triangles at interface of 2 domains with same reference */
  if (nargs[i]) MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_lag,           arg(i,stack,0L));    i++;   /*!< [-1/0/1/2], Lagrangian option */
  if (nargs[i]) MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_optim,         arg(i,stack,false)); i++;   /*!< [1/0], Optimize mesh keeping its initial edge sizes */
  if (nargs[i]) MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_optimLES,      arg(i,stack,false)); i++;   /*!< [1/0], Strong mesh optimization for Les computations */
  if (nargs[i]) MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_noinsert,      arg(i,stack,false)); i++;   /*!< [1/0], Avoid/allow point insertion */
  if (nargs[i]) MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_noswap,        arg(i,stack,false)); i++;   /*!< [1/0], Avoid/allow edge or face flipping */
  if (nargs[i]) MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_nomove,        arg(i,stack,false)); i++;   /*!< [1/0], Avoid/allow point relocation */
  if (nargs[i]) MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_nosurf,        arg(i,stack,false)); i++;   /*!< [1/0], Avoid/allow surface modifications */
  //if (nargs[i]) MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_nreg,          arg(i,stack,false)); i++;   /*!< [0/1], Enable normal regularization */
  if (nargs[i]) MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_renum,         arg(i,stack,false)); i++;   /*!< [1/0], Turn on/off point relocation with Scotch */
  if (nargs[i]) MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_anisosize,     arg(i,stack,false)); i++;   /*!< [1/0], Turn on/off anisotropic metric creation when no metric is provided */
  if (nargs[i]) MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_octree,        arg(i,stack,0L));    i++;   /*!< [n], Specify the max number of points per PROctree cell (DELAUNAY) */
  if (nargs[i]) MMG3D_Set_dparameter(mesh,sol,MMG3D_DPARAM_angleDetection,arg(i,stack,0.));    i++;   /*!< [val], Value for angle detection */
  if (nargs[i]) MMG3D_Set_dparameter(mesh,sol,MMG3D_DPARAM_hmin,          arg(i,stack,0.));    i++;   /*!< [val], Minimal mesh size */
  if (nargs[i]) MMG3D_Set_dparameter(mesh,sol,MMG3D_DPARAM_hmax,          arg(i,stack,0.));    i++;   /*!< [val], Maximal mesh size */
  if (nargs[i]) MMG3D_Set_dparameter(mesh,sol,MMG3D_DPARAM_hsiz,          arg(i,stack,0.));    i++;   /*!< [val], Constant mesh size */
  if (nargs[i]) MMG3D_Set_dparameter(mesh,sol,MMG3D_DPARAM_hausd,         arg(i,stack,0.));    i++;   /*!< [val], Control global Hausdorff distance (on all the boundary surfaces of the mesh) */
  if (nargs[i]) MMG3D_Set_dparameter(mesh,sol,MMG3D_DPARAM_hgrad,         arg(i,stack,0.));    i++;   /*!< [val], Control gradation */
  if (nargs[i]) MMG3D_Set_dparameter(mesh,sol,MMG3D_DPARAM_ls,            arg(i,stack,0.));    i++;   /*!< [val], Value of level-set */
  //if (nargs[i]) MMG3D_Set_dparameter(mesh,sol,MMG3D_DPARAM_rmc,           arg(i,stack,0.));    i++;   /*!< [-1/val], Remove small connex componants in level-set mode */

  bool bls = MMG3D_Get_iparameter(mesh,MMG3D_IPARAM_iso);

  int ier;
  if (!bls)
    ier = MMG3D_mmg3dlib(mesh,sol);
  else
    ier = MMG3D_mmg3dls(mesh,sol,met);
  
  Mesh3 *Th_T = nullptr;
  
  MMG5_pMesh_to_ffmesh(mesh,Th_T);

  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&sol,
                 MMG5_ARG_end);

  Th_T->BuildGTree();
  
  Add2StackOfPtr2FreeRC(stack, Th_T);
  return Th_T;
}

template<>
AnyType mmg_Op<MeshS>::operator( )(Stack stack) const {
  // initialisation
  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
  MeshS *pTh = GetAny< MeshS * >((*eTh)(stack));

  ffassert(pTh);
  MeshS &Th = *pTh;
  int nv = Th.nv;
  int nt = Th.nt;
  int nbe = Th.nbe;

  KN< double > *pmetric = 0;
  KN< long > *prequiredEdge = 0;

  if (nargs[0]) {
    pmetric = GetAny< KN< double > * >((*nargs[0])(stack));
  }
  if (nargs[19]) {
    prequiredEdge = GetAny< KN< long > * >((*nargs[19])(stack));
  }

  MMG5_pMesh mesh;
  MMG5_pSol sol;

  mesh = nullptr;
  sol = nullptr;

  MMGS_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&sol,
                  MMG5_ARG_end);

  ffmesh_to_MMG5_pMesh(Th, mesh);
  
    if (pmetric && pmetric->N( ) > 0) {
      const KN< double > &metric = *pmetric;
      if (metric.N( ) == Th.nv) {
        if ( MMGS_Set_solSize(mesh,sol,MMG5_Vertex,Th.nv,MMG5_Scalar) != 1 ) { 
          printf("Unable to allocate the metric array.\n");
          exit(EXIT_FAILURE);
        }
        if ( MMGS_Set_scalarSols(sol,metric) != 1 ) { 
          printf("Unable to set metric.\n");
          exit(EXIT_FAILURE);
        }
      }
      else {
        if ( MMGS_Set_solSize(mesh,sol,MMG5_Vertex,Th.nv,MMG5_Tensor) != 1 ) { 
          printf("Unable to allocate the metric array.\n");
          exit(EXIT_FAILURE);
        }
        static const int perm[6] = {0, 1, 3, 2, 4, 5};
        for (int k=0; k<Th.nv; k++) {
          if ( MMGS_Set_tensorSol(sol, metric[6*k+perm[0]], metric[6*k+perm[1]], metric[6*k+perm[2]], 
                                  metric[6*k+perm[3]], metric[6*k+perm[4]], metric[6*k+perm[5]], k+1) != 1 ) { 
            printf("Unable to set metric.\n");
            exit(EXIT_FAILURE);
          }
        }
      }
    }
    if (prequiredEdge && prequiredEdge->N( ) > 0) {
      const KN< long > &requiredEdge = *prequiredEdge;
      std::sort(requiredEdge + 0, requiredEdge + requiredEdge.N());
      int na;
      if ( MMGS_Get_meshSize(mesh,NULL,NULL,&na) !=1 ) {
        exit(EXIT_FAILURE);
      }
      for (int k=1; k<=na; k++) {
        int ref, dummy;
        if ( MMGS_Get_edge(mesh, &dummy, &dummy, &ref,
                                  &dummy, &dummy) != 1 ) {
          exit(EXIT_FAILURE);
        }
        if (std::binary_search(requiredEdge + 0, requiredEdge + requiredEdge.N(), ref)) {
          if ( MMG3D_Set_requiredEdge(mesh,k) != 1 ) {
            exit(EXIT_FAILURE);
          }
        }
      }
    }

  long iso=0L;

  int i=1;
  if (nargs[i]) MMGS_Set_iparameter(mesh,sol,MMGS_IPARAM_verbose,       arg(i,stack,0L));    i++;   /*!< [-1..10], Tune level of verbosity */
  if (nargs[i]) MMGS_Set_iparameter(mesh,sol,MMGS_IPARAM_mem,           arg(i,stack,0L));    i++;   /*!< [n/-1], Set memory size to n Mbytes or keep the default value */
  if (nargs[i]) MMGS_Set_iparameter(mesh,sol,MMGS_IPARAM_debug,         arg(i,stack,false)); i++;   /*!< [1/0], Turn on/off debug mode */
  if (nargs[i]) MMGS_Set_iparameter(mesh,sol,MMGS_IPARAM_angle,         arg(i,stack,false)); i++;   /*!< [1/0], Turn on/off angle detection */
  if (nargs[i]) {iso = arg(i,stack,false); MMGS_Set_iparameter(mesh,sol,MMGS_IPARAM_iso,iso);} i++; /*!< [1/0], Level-set meshing */
  if (nargs[i]) MMGS_Set_iparameter(mesh,sol,MMGS_IPARAM_keepRef,       arg(i,stack,false)); i++;   /*!< [1/0], Preserve the initial domain references in level-set mode */
  //if (nargs[i]) MMGS_Set_iparameter(mesh,sol,MMGS_IPARAM_optim,       arg(i,stack,false));   i++;   /*!< [1/0], Optimize mesh keeping its initial edge sizes */
  if (nargs[i]) MMGS_Set_iparameter(mesh,sol,MMGS_IPARAM_noinsert,      arg(i,stack,false)); i++;   /*!< [1/0], Avoid/allow point insertion */
  if (nargs[i]) MMGS_Set_iparameter(mesh,sol,MMGS_IPARAM_noswap,        arg(i,stack,false)); i++;   /*!< [1/0], Avoid/allow edge or face flipping */
  if (nargs[i]) MMGS_Set_iparameter(mesh,sol,MMGS_IPARAM_nomove,        arg(i,stack,false)); i++;   /*!< [1/0], Avoid/allow point relocation */
  if (nargs[i]) MMGS_Set_iparameter(mesh,sol,MMGS_IPARAM_nreg,          arg(i,stack,false)); i++;   /*!< [0/1], Disabled/enabled normal regularization */
  if (nargs[i]) MMGS_Set_iparameter(mesh,sol,MMGS_IPARAM_renum,         arg(i,stack,false)); i++;   /*!< [1/0], Turn on/off point relocation with Scotch */
  if (nargs[i]) MMGS_Set_dparameter(mesh,sol,MMGS_DPARAM_angleDetection,arg(i,stack,0.));    i++;   /*!< [val], Value for angle detection */
  if (nargs[i]) MMGS_Set_dparameter(mesh,sol,MMGS_DPARAM_hmin,          arg(i,stack,0.));    i++;   /*!< [val], Minimal mesh size */
  if (nargs[i]) MMGS_Set_dparameter(mesh,sol,MMGS_DPARAM_hmax,          arg(i,stack,0.));    i++;   /*!< [val], Maximal mesh size */
  if (nargs[i]) MMGS_Set_dparameter(mesh,sol,MMGS_DPARAM_hsiz,          arg(i,stack,0.));    i++;   /*!< [val], Constant mesh size */
  if (nargs[i]) MMGS_Set_dparameter(mesh,sol,MMGS_DPARAM_hausd,         arg(i,stack,0.));    i++;   /*!< [val], Control global Hausdorff distance (on all the boundary surfaces of the mesh) */
  if (nargs[i]) MMGS_Set_dparameter(mesh,sol,MMGS_DPARAM_hgrad,         arg(i,stack,0.));    i++;   /*!< [val], Control gradation */
  if (nargs[i]) MMGS_Set_dparameter(mesh,sol,MMGS_DPARAM_ls,            arg(i,stack,0.));    i++;   /*!< [val], Value of level-set */

  int ier;
  if (!iso)
    ier = MMGS_mmgslib(mesh,sol);
  else
    ier = MMGS_mmgsls(mesh,sol,NULL);

  MeshS *Th_T = nullptr;
  
  MMG5_pMesh_to_ffmesh(mesh,Th_T);

  MMGS_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&sol,
                 MMG5_ARG_end);

  Th_T->BuildGTree();
  
  Add2StackOfPtr2FreeRC(stack, Th_T);
  return Th_T;
}

static void Load_Init( ) {
  if (verbosity && mpirank == 0) {
    cout << " load: mmg " << endl;
  }

  Global.Add("mmg3d", "(", new mmg_ff<Mesh3>);
  Global.Add("mmgs", "(", new mmg_ff<MeshS>);
}

LOADFUNC(Load_Init)
