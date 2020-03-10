//ff-c++-LIBRARY-dep: parmmg mmg mpi metis scotch
//ff-c++-cpp-dep:

#include "ff++.hpp"
#include "memory.h"
#include "parmmg/libparmmg.h"
#include "GenericMesh.hpp"

#if 0
extern "C" int PMMG_grp_to_saveMesh(PMMG_pParMesh, int, char*);
#endif

using namespace Fem2D;

int ffmesh_to_PMMG_pParMesh(const Mesh3 &Th, PMMG_pParMesh& mesh, bool distributed) {
  int nVertices       = Th.nv;
  int nTetrahedra     = Th.nt;
  int nPrisms         = 0;
  int nTriangles      = Th.nbe;
  int nQuadrilaterals = 0;
  int nEdges          = 0;

  if(mesh->myrank == mesh->info.root || distributed) {
    if ( PMMG_Set_meshSize(mesh,nVertices,nTetrahedra,nPrisms,nTriangles,
                           nQuadrilaterals,nEdges) != 1 ) {
      exit(EXIT_FAILURE);
    }

    for (int k = 0; k < Th.nv; k++) {
      if ( PMMG_Set_vertex(mesh,Th.vertices[k].x,Th.vertices[k].y,
                           Th.vertices[k].z, Th.vertices[k].lab, k+1) != 1 ) {
        exit(EXIT_FAILURE);
      }
    }

    for (int k = 0; k < Th.nt; k++) {
      const Tet &K(Th.elements[k]);
      if ( PMMG_Set_tetrahedron(mesh,Th.operator()(K[0])+1,Th.operator()(K[1])+1,
                                Th.operator()(K[2])+1,Th.operator()(K[3])+1,K.lab,k+1) != 1 ) {
        exit(EXIT_FAILURE);
      }
    }

    for (int k = 0; k < Th.nbe; k++) {
      const Triangle3 &K(Th.be(k));
      if ( PMMG_Set_triangle(mesh,Th.operator()(K[0])+1,Th.operator()(K[1])+1,Th.operator()(K[2])+1,
                             K.lab,k+1) != 1 ) {
        exit(EXIT_FAILURE);
      }
    }
  }
  return 0;
}

int PMMG_pParMesh_to_ffmesh(const PMMG_pParMesh& mesh, Mesh3 *&T_TH3, bool distributed) {
    int ier;

    int nVertices   = 0;
    int nTetrahedra = 0;
    int nTriangles  = 0;
    int nEdges      = 0;

    if(mesh->myrank == mesh->info.root || distributed) {
      if ( PMMG_Get_meshSize(mesh,&nVertices,&nTetrahedra,NULL,&nTriangles,NULL,
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
        if ( PMMG_Get_vertex(mesh,&(v[k].x),&(v[k].y),&(v[k].z),
                                     &(v[k].lab),&(corner),&(required)) != 1 ) {
          cout << "Unable to get mesh vertex " << k << endl;
          ier = MMG5_STRONGFAILURE;
        }
      }

      for ( k=0; k<nTetrahedra; k++ ) {
        int iv[4], lab;
        if ( PMMG_Get_tetrahedron(mesh,
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
        if ( PMMG_Get_triangle(mesh,
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
    }
    return 0;
}

class parmmg_Op : public E_F0mps {
 public:
  Expression eTh, xx, yy, zz;
  static const int n_name_param = 31;
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

 public:
  parmmg_Op(const basicAC_F0 &args, Expression tth) : eTh(tth) {
    if (verbosity > 1) {
      cout << "parmmg" << endl;
    }

    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type parmmg_Op::name_param[] = {
{"metric"            , &typeid(KN< double > *)},
{"comm"              , &typeid(pcommworld)},
{"verbose"           , &typeid(long)},/*!< [-10..10], Tune level of verbosity */
{"mmgVerbose"        , &typeid(long)},/*!< [-10..10], Tune level of verbosity of Mmg */
{"mem"               , &typeid(long)},/*!< [n/-1], Set memory size to n Mbytes or keep the default value */
{"debug"             , &typeid(bool)},/*!< [1/0], Turn on/off debug mode */
{"mmgDebug"          , &typeid(bool)},/*!< [1/0], Turn on/off debug mode */
{"angle"             , &typeid(bool)},/*!< [1/0], Turn on/off angle detection */
{"iso"               , &typeid(bool)},/*!< [1/0], Level-set meshing */
{"lag"               , &typeid(long)},/*!< [-1/0/1/2], Lagrangian option */
{"optim"             , &typeid(bool)},/*!< [1/0], Optimize mesh keeping its initial edge sizes */
{"optimLES"          , &typeid(bool)},/*!< [1/0], Strong mesh optimization for Les computations */
{"noinsert"          , &typeid(bool)},/*!< [1/0], Avoid/allow point insertion */
{"noswap"            , &typeid(bool)},/*!< [1/0], Avoid/allow edge or face flipping */
{"nomove"            , &typeid(bool)},/*!< [1/0], Avoid/allow point relocation */
{"nosurf"            , &typeid(bool)},/*!< [1/0], Avoid/allow surface modifications */
{"anisosize"         , &typeid(bool)},/*!< [1/0], Turn on/off anisotropic metric creation when no metric is provided */
{"octree"            , &typeid(long)},/*!< [n], Specify the max number of points per octree cell (DELAUNAY) */
{"meshSize"          , &typeid(long)},/*!< [n], Target mesh size of Mmg (advanced use) */
{"metisRatio"        , &typeid(long)},/*!< [n], wanted ratio # mesh / # metis super nodes (advanced use) */
{"ifcLayers"         , &typeid(long)},/*!< [n], Number of layers of interface displacement */
{"groupsRatio"       , &typeid(double)},/*!< [val], Allowed imbalance between current and desired groups size */
{"niter"             , &typeid(long)},/*!< [n], Set the number of remeshing iterations */
{"angleDetection"    , &typeid(double)},/*!< [val], Value for angle detection */
{"hmin"              , &typeid(double)},/*!< [val], Minimal mesh size */
{"hmax"              , &typeid(double)},/*!< [val], Maximal mesh size */
{"hsiz"              , &typeid(double)},/*!< [val], Constant mesh size */
{"hausd"             , &typeid(double)},/*!< [val], Control global Hausdorff distance (on all the boundary surfaces of the mesh) */
{"hgrad"             , &typeid(double)},/*!< [val], Control gradation */
{"ls"                , &typeid(double)},/*!< [val], Value of level-set */
{"nodeCommunicators" , &typeid(KN<KN<long>>*)}
};

class parmmg_ff : public OneOperator {
 public:
  parmmg_ff( ) : OneOperator(atype< const Mesh3* >( ), atype< const Mesh3* >( )) {}

  E_F0 *code(const basicAC_F0 &args) const { return new parmmg_Op(args, t[0]->CastTo(args[0])); }
};

AnyType parmmg_Op::operator( )(Stack stack) const {
  // initialisation
  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
  Mesh3 *pTh = GetAny< Mesh3 * >((*eTh)(stack));

  ffassert(pTh);
  Mesh3 &Th = *pTh;
  int nv = Th.nv;
  int nt = Th.nt;
  int nbe = Th.nbe;

  KN< double > *pmetric = 0;

  if (nargs[0]) {
    pmetric = GetAny< KN< double > * >((*nargs[0])(stack));
  }

  pcommworld pcomm = 0;

  if (nargs[1]) {
    pcomm = GetAny<pcommworld>((*nargs[1])(stack));
  }

  MPI_Comm comm = pcomm ? *(MPI_Comm*)pcomm : MPI_COMM_WORLD;

  PMMG_pParMesh mesh;
  MMG5_pSol sol;
  MMG5_pSol met;

  mesh = nullptr;
  sol = nullptr;
  met = nullptr;

  PMMG_Init_parMesh(PMMG_ARG_start,
                    PMMG_ARG_ppParMesh,&mesh,
                    PMMG_ARG_pMesh,PMMG_ARG_pMet,
                    PMMG_ARG_dim,3,PMMG_ARG_MPIComm,comm,
                    PMMG_ARG_end);

  KN< KN< long > >* communicators = nargs[30] ? GetAny< KN< KN< long > >* >((*nargs[30])(stack)) : 0;
  ffmesh_to_PMMG_pParMesh(Th, mesh, communicators != NULL);

  int root = mesh->info.root;
  int myrank = mesh->myrank;

  if(myrank == root || communicators != NULL) {
    if (pmetric && pmetric->N( ) > 0) {
      const KN< double > &metric = *pmetric;
      if (metric.N( ) == Th.nv) {
        if ( PMMG_Set_metSize(mesh,MMG5_Vertex,Th.nv,MMG5_Scalar) != 1 ) {
          printf("Unable to allocate the metric array.\n");
          exit(EXIT_FAILURE);
        }
        if ( PMMG_Set_scalarMets(mesh,metric) != 1 ) {
          printf("Unable to set metric.\n");
          exit(EXIT_FAILURE);
        }
      }
      else {
        if ( PMMG_Set_metSize(mesh,MMG5_Vertex,Th.nv,MMG5_Tensor) != 1 ) {
          printf("Unable to allocate the metric array.\n");
          exit(EXIT_FAILURE);
        }
        static const int perm[6] = {0, 1, 3, 2, 4, 5};
        for (int k=0; k<Th.nv; k++) {
          if ( PMMG_Set_tensorMet(mesh, metric[6*k+perm[0]], metric[6*k+perm[1]], metric[6*k+perm[2]],
                                  metric[6*k+perm[3]], metric[6*k+perm[4]], metric[6*k+perm[5]], k+1) != 1 ) {
            printf("Unable to set metric.\n");
            exit(EXIT_FAILURE);
          }
        }
      }
    }
  }

  int i=2;

  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_verbose,       arg(i,stack,0L)); i++;   /*!< [-10..10], Tune level of verbosity */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_mmgVerbose,    arg(i,stack,0L)); i++;   /*!< [-10..10], Tune level of verbosity of Mmg */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_mem,           arg(i,stack,0L)); i++;   /*!< [n/-1], Set memory size to n Mbytes or keep the default value */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_debug,         arg(i,stack,0L)); i++;   /*!< [1/0], Turn on/off debug mode */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_mmgDebug,      arg(i,stack,0L)); i++;   /*!< [1/0], Turn on/off debug mode */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_angle,         arg(i,stack,0L)); i++;   /*!< [1/0], Turn on/off angle detection */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_iso,           arg(i,stack,0L)); i++;   /*!< [1/0], Level-set meshing */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_lag,           arg(i,stack,0L)); i++;   /*!< [-1/0/1/2], Lagrangian option */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_optim,         arg(i,stack,0L)); i++;   /*!< [1/0], Optimize mesh keeping its initial edge sizes */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_optimLES,      arg(i,stack,0L)); i++;   /*!< [1/0], Strong mesh optimization for Les computations */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_noinsert,      arg(i,stack,0L)); i++;   /*!< [1/0], Avoid/allow point insertion */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_noswap,        arg(i,stack,0L)); i++;   /*!< [1/0], Avoid/allow edge or face flipping */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_nomove,        arg(i,stack,0L)); i++;   /*!< [1/0], Avoid/allow point relocation */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_nosurf,        arg(i,stack,0L)); i++;   /*!< [1/0], Avoid/allow surface modifications */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_anisosize,     arg(i,stack,0L)); i++;   /*!< [1/0], Turn on/off anisotropic metric creation when no metric is provided */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_octree,        arg(i,stack,0L)); i++;   /*!< [n], Specify the max number of points per octree cell (DELAUNAY) */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_meshSize,      arg(i,stack,0L)); i++;   /*!< [n], Target mesh size of Mmg (advanced use) */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_metisRatio,    arg(i,stack,0L)); i++;   /*!< [n], wanted ratio # mesh / # metis super nodes (advanced use) */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_ifcLayers,     arg(i,stack,0L)); i++;   /*!< [n], Number of layers of interface displacement */
  if (nargs[i]) PMMG_Set_dparameter(mesh,PMMG_DPARAM_groupsRatio,   arg(i,stack,0.)); i++;   /*!< [val], Allowed imbalance between current and desired groups size */
  if (nargs[i]) PMMG_Set_iparameter(mesh,PMMG_IPARAM_niter,         arg(i,stack,0L)); i++;   /*!< [n], Set the number of remeshing iterations */
  if (nargs[i]) PMMG_Set_dparameter(mesh,PMMG_DPARAM_angleDetection,arg(i,stack,0.)); i++;   /*!< [val], Value for angle detection */
  if (nargs[i]) PMMG_Set_dparameter(mesh,PMMG_DPARAM_hmin,          arg(i,stack,0.)); i++;   /*!< [val], Minimal mesh size */
  if (nargs[i]) PMMG_Set_dparameter(mesh,PMMG_DPARAM_hmax,          arg(i,stack,0.)); i++;   /*!< [val], Maximal mesh size */
  if (nargs[i]) PMMG_Set_dparameter(mesh,PMMG_DPARAM_hsiz,          arg(i,stack,0.)); i++;   /*!< [val], Constant mesh size */
  if (nargs[i]) PMMG_Set_dparameter(mesh,PMMG_DPARAM_hausd,         arg(i,stack,0.)); i++;   /*!< [val], Control global Hausdorff distance (on all the boundary surfaces of the mesh) */
  if (nargs[i]) PMMG_Set_dparameter(mesh,PMMG_DPARAM_hgrad,         arg(i,stack,0.)); i++;   /*!< [val], Control gradation */
  if (nargs[i]) PMMG_Set_dparameter(mesh,PMMG_DPARAM_ls,            arg(i,stack,0.)); i++;   /*!< [val], Value of level-set */

  if(communicators != NULL) {
    /* Set API mode */
    if( !PMMG_Set_iparameter( mesh, PMMG_IPARAM_APImode, PMMG_APIDISTRIB_nodes ) ) {
      exit(EXIT_FAILURE);
    }

    /* Set the number of interfaces */
    PMMG_Set_numberOfNodeCommunicators(mesh, communicators->operator[](0).N());

    /* Loop on each interface (proc pair) seen by the current rank) */
    for(int icomm=0; icomm<communicators->operator[](0).N(); icomm++ ) {

      /* Set nb. of entities on interface and rank of the outward proc */
      PMMG_Set_ithNodeCommunicatorSize(mesh, icomm,
                                       communicators->operator[](0)[icomm],
                                       communicators->operator[](1 + 2 * icomm).N());

      /* Set local and global index for each entity on the interface */
      KN<int> local = communicators->operator[](1 + 2 * icomm);
      KN<int> global = communicators->operator[](2 + 2 * icomm);
      PMMG_Set_ithNodeCommunicator_nodes(mesh, icomm,
                                         local.operator int*(),
                                         global.operator int*(), 1);
    }
  }
#if 0
  char filemesh[48];
  sprintf(filemesh,"cube_in.%d.mesh",mpirank);
  // PMMG_grp_to_saveMesh(mesh, mpirank, filemesh);
#endif
  int ier = communicators == NULL ? PMMG_parmmglib_centralized(mesh) : PMMG_parmmglib_distributed(mesh);

  Mesh3 *Th_T = nullptr;

  PMMG_pParMesh_to_ffmesh(mesh, Th_T, communicators != NULL);

  PMMG_Free_all(PMMG_ARG_start,
              PMMG_ARG_ppParMesh, &mesh,
              PMMG_ARG_end);
  if(communicators == NULL) {
    Serialize *buf = 0;
    long nbsize = 0;
    if (myrank == root) {
      buf = new Serialize((*Th_T).serialize());
      nbsize = buf->size();
    }

    MPI_Bcast(reinterpret_cast<void*>(&nbsize), 1, MPI_LONG, root, comm);
    if (myrank != root)
      buf = new Serialize(nbsize, Fem2D::GenericMesh_magicmesh);
    MPI_Bcast(reinterpret_cast<void*>((char *)(*buf)), nbsize, MPI_BYTE, root, comm);
    if (myrank != root)
      Th_T = new Fem2D::Mesh3(*buf);
    delete buf;
  }

  Th_T->BuildGTree();

  Add2StackOfPtr2FreeRC(stack, Th_T);
  return Th_T;
}

static void Load_Init( ) {
  if (verbosity) {
    cout << " load: parmmg  " << endl;
  }

  Global.Add("parmmg3d", "(", new parmmg_ff);
}

LOADFUNC(Load_Init)
