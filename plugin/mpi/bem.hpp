#ifndef BEM_HPP_
#define BEM_HPP_

#include <htool/htool.hpp>

#if defined(WITH_metis)
extern "C" {
#include <metis.h>
}

#ifdef METIS_VER_MAJOR
extern "C" {
real_t libmetis__ComputeElementBalance(idx_t ne, idx_t nparts, idx_t *where);
}
#else
typedef idxtype idx_t;
#endif

template< class FESPACE, int NO, typename R >
KN< R > *partmetis( KN< R > *const &part, const FESPACE * pVh, long const &lparts) {
  ffassert(pVh);
  const FESPACE &Vh(*pVh);
   int nve = Vh[0].NbDoF( );
  const typename FESPACE::Mesh & Th = Vh.Th;
    idx_t nt = Th.nt, nv = Vh.NbOfDF;

  KN< idx_t > eptr(nt + 1), elmnts(nve * nt), epart(nt), npart(nv);
  if(lparts > 1) {
      for (idx_t k = 0, i = 0; k < nt; ++k) {
          eptr[k] = i;
          for (idx_t j = 0; j < nve; j++) {
              elmnts[i++] = Vh(k, j);
          }
          eptr[k + 1] = i;
      }

      idx_t numflag = 0;
      idx_t nparts = lparts;
      idx_t edgecut;
      idx_t etype = nve - 2;    // triangle or tet .  change FH fevr 2010
      idx_t ncommon = 1;
#ifdef METIS_VER_MAJOR
      if (NO == 0) {
          METIS_PartMeshNodal(&nt, &nv, eptr, (idx_t *)elmnts, 0, 0, &nparts, 0, 0, &edgecut,
                  (idx_t *)epart, (idx_t *)npart);
      } else {
          METIS_PartMeshDual(&nt, &nv, eptr, (idx_t *)elmnts, 0, 0, &ncommon, &nparts, 0, 0, &edgecut,
                  (idx_t *)epart, (idx_t *)npart);
      }

      if (verbosity) {
          printf("  --metis: %d-way Edge-Cut: %7d, Balance: %5.2f Nodal=0/Dual %d\n", nparts, nve,
                  libmetis__ComputeElementBalance(nt, nparts, epart), NO);
      }

#else
      if (NO == 0) {
          METIS_PartMeshNodal(&nt, &nv, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart);
      } else {
          METIS_PartMeshDual(&nt, &nv, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart);
      }

      if (verbosity) {
          printf("  --metis: %d-way Edge-Cut: %7d, Balance: %5.2f Nodal=0/Dual %d\n", nparts, nve,
                  ComputeElementBalance(nt, nparts, epart), NO);
      }

#endif
  } else npart = 0;
  part->resize(nv);
  *part = npart;
  return part;
}
#endif

double ff_htoolEta=10., ff_htoolEpsilon=1e-3;
long ff_htoolMinclustersize=10, ff_htoolMintargetdepth=0, ff_htoolMinsourcedepth=0;

template<class K>
class HMatrixVirt {
public:
    virtual const std::map<std::string, std::string>& get_infos() const = 0;
    virtual void mvprod_global(const K* const in, K* const out,const int& mu=1) const = 0;
    virtual int nb_rows() const = 0;
    virtual int nb_cols() const = 0;
    virtual void cluster_to_target_permutation(const K* const in, K* const out) const = 0;
    virtual void source_to_cluster_permutation(const K* const in, K* const out) const = 0;
    virtual MPI_Comm get_comm() const = 0;
    virtual int get_rankworld() const = 0;
    virtual int get_sizeworld() const = 0;
    virtual int get_local_size() const = 0;
    virtual int get_local_offset() const = 0;
    virtual std::vector<const htool::HMatrix<K>*> get_dense_blocks() const = 0;
    virtual std::vector<const htool::HMatrix<K>*> get_low_rank_blocks() const = 0;
    virtual std::vector<const htool::HMatrix<K>*> get_diagonal_blocks() const = 0;
    virtual htool::Matrix<K> get_local_dense() const = 0;
    virtual int get_permt(int) const = 0;
    virtual int get_perms(int) const = 0;

    virtual ~HMatrixVirt() {};
};


template<class K>
class HMatrixImpl : public HMatrixVirt<K> {
private:
    std::shared_ptr<htool::Cluster<double>>  target_cluster;
    std::shared_ptr<htool::Cluster<double>>  source_cluster;
    htool::DistributedOperatorFromHMatrix<K> distributed_operator_holder;


public:
    const htool::HMatrix<K>& H;
    const htool::DistributedOperator<K>& distributed_operator;
    std::map<std::string, std::string> infos;
    HMatrixImpl(htool::VirtualGenerator<K> &mat, std::shared_ptr<htool::Cluster<double>> t, std::shared_ptr<htool::Cluster<double>> s,const htool::HMatrixTreeBuilder<K,double>& hmatrix_tree_builder,MPI_Comm comm) : target_cluster(t), source_cluster(s), distributed_operator_holder(mat,*target_cluster,*source_cluster,hmatrix_tree_builder,comm),distributed_operator(distributed_operator_holder.distributed_operator), H(distributed_operator_holder.hmatrix), infos(htool::get_distributed_hmatrix_information(H,distributed_operator.get_comm())){}
    const std::map<std::string, std::string>& get_infos() const { return infos; }
    void mvprod_global(const K* const in, K* const out,const int& mu=1) const {return distributed_operator.matrix_product_global_to_global(in,out,mu);}
    int nb_rows() const { return distributed_operator.get_target_partition().get_global_size();}
    int nb_cols() const { return distributed_operator.get_source_partition().get_global_size();}
    void cluster_to_target_permutation(const K* const in, K* const out) const {return htool::cluster_to_user(*target_cluster,in,out);}
    void source_to_cluster_permutation(const K* const in, K* const out) const {return htool::user_to_cluster(*source_cluster,in,out);}
    MPI_Comm get_comm() const {return distributed_operator.get_comm();}
    int get_rankworld() const {int rankworld; MPI_Comm_rank(distributed_operator.get_comm(),&rankworld); return rankworld;}
    int get_sizeworld() const {int sizeworld; MPI_Comm_size(distributed_operator.get_comm(),&sizeworld); return sizeworld;}
    int get_local_size() const {int rankworld; MPI_Comm_rank(distributed_operator.get_comm(),&rankworld); return distributed_operator.get_target_partition().get_size_of_partition(rankworld);}
    int get_local_offset() const {int rankworld; MPI_Comm_rank(distributed_operator.get_comm(),&rankworld); return distributed_operator.get_target_partition().get_offset_of_partition(rankworld);}
    std::vector<const htool::HMatrix<K>*> get_dense_blocks()const{    
        std::vector<const htool::HMatrix<K>*> dense_blocks;
        preorder_tree_traversal(
        H,
        [&dense_blocks](const htool::HMatrix<K> &current_hmatrix) {
            if (current_hmatrix.is_leaf() && current_hmatrix.is_dense()) {
                dense_blocks.push_back(&current_hmatrix);
            }
        });
        return dense_blocks;}
    std::vector<const htool::HMatrix<K>*> get_low_rank_blocks()const{    
        std::vector<const htool::HMatrix<K>*> low_rank_blocks;
        preorder_tree_traversal(
        H,
        [&low_rank_blocks](const htool::HMatrix<K> &current_hmatrix) {
            if (current_hmatrix.is_leaf() && current_hmatrix.is_low_rank()) {
                low_rank_blocks.push_back(&current_hmatrix);
            }
        });
        return low_rank_blocks;}
    std::vector<const htool::HMatrix<K>*> get_diagonal_blocks() const{    
        std::vector<const htool::HMatrix<K>*> diagonal_blocks;
        preorder_tree_traversal(
        H,
        [&diagonal_blocks](const htool::HMatrix<K> &current_hmatrix) {
            if (current_hmatrix.is_dense() && current_hmatrix.get_target_cluster().get_offset()==current_hmatrix.get_source_cluster().get_offset()) {
                diagonal_blocks.push_back(&current_hmatrix);
            }
        });
        return diagonal_blocks;}
    htool::Matrix<K> get_local_dense() const {htool::Matrix<K> dense(H.nb_rows(),H.nb_cols());htool::copy_to_dense(H,dense.data());return dense;}
    int get_permt(int i) const {return H.get_target_cluster().get_permutation()[i];}
    int get_perms(int i) const {return H.get_source_cluster().get_permutation()[i];}
};


struct Data_Bem_Solver
: public Data_Sparse_Solver {
    double eta;
    int minclustersize,mintargetdepth,minsourcedepth;
    string compressor, initialclustering;
       
    Data_Bem_Solver()
       : Data_Sparse_Solver(),
       eta(ff_htoolEta),
       minclustersize(ff_htoolMinclustersize),
       mintargetdepth(ff_htoolMintargetdepth),
       minsourcedepth(ff_htoolMinsourcedepth),
       compressor("partialACA"),
       initialclustering("pca")
    
        {epsilon=ff_htoolEpsilon;}
     
    template<class R>
       void Init_sym_positive_var();
    
    private:
        Data_Bem_Solver(const Data_Bem_Solver& ); // pas de copie
        
};

template<class R>
void Data_Bem_Solver::Init_sym_positive_var()
{
    
    Data_Sparse_Solver::Init_sym_positive_var<R>(-1);
    
}

template<class R>
inline void SetEnd_Data_Bem_Solver(Stack stack,Data_Bem_Solver & ds,Expression const *nargs ,int n_name_param)
{
    
    {
        bool unset_eps=true;
        ds.initmat=true;
        ds.factorize=0;
        int kk = n_name_param-(NB_NAME_PARM_MAT+NB_NAME_PARM_HMAT)-1;
        if (nargs[++kk]) ds.initmat= ! GetAny<bool>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.solver= * GetAny<string*>((*nargs[kk])(stack));
        ds.Init_sym_positive_var<R>();//  set def value of sym and posi
        if (nargs[++kk]) ds.epsilon= GetAny<double>((*nargs[kk])(stack)),unset_eps=false;
        if (nargs[++kk])
        {// modif FH fev 2010 ...
            const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(nargs[kk]);
            if(op)
            {
                ds.precon = op->Find("(",ArrayOfaType(atype<KN<R>* >(),false)); // strange bug in g++ is R become a double
                ffassert(ds.precon);
            } // add miss
        }
        
        if (nargs[++kk]) ds.NbSpace= GetAny<long>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.tgv= GetAny<double>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.factorize= GetAny<long>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.strategy = GetAny<long>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.tol_pivot = GetAny<double>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.tol_pivot_sym = GetAny<double>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.itmax = GetAny<long>((*nargs[kk])(stack)); //  frev 2007 OK
        if (nargs[++kk]) ds.data_filename = *GetAny<string*>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.lparams = GetAny<KN_<long> >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.dparams = GetAny<KN_<double> >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.smap = GetAny<MyMap<String,String> *>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.perm_r = GetAny<KN_<long> >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.perm_c = GetAny<KN_<long> >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.scale_r = GetAny<KN_<double> >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.scale_c = GetAny<KN_<double> >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.sparams = *GetAny<string*>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.commworld = GetAny<pcommworld>((*nargs[kk])(stack));
#ifdef VDATASPARSESOLVER
        if (nargs[++kk]) ds.master = GetAny<long>((*nargs[kk])(stack));
#else
        ++kk;
#endif
        if (nargs[++kk]) ds.rinfo = GetAny<KN<double>* >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.info = GetAny<KN<long>* >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.kerneln = GetAny< KNM<double>* >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.kernelt = GetAny< KNM<double>* >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.kerneldim = GetAny<long * >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.verb = GetAny<long  >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.x0 = GetAny<bool>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.veps= GetAny<double*>((*nargs[kk])(stack));
        if( unset_eps && ds.veps) ds.epsilon = *ds.veps;//  if veps  and no def value  => veps def value of epsilon.
        if (nargs[++kk]) ds.rightprecon= GetAny<bool>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.sym= GetAny<bool>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.positive= GetAny<bool>((*nargs[kk])(stack));
        if (nargs[++kk])  { ds.getnbiter= GetAny<long*>((*nargs[kk])(stack));
                   if( ds.getnbiter) *ds.getnbiter=-1; //undef
               }
        if(ds.solver == "") { // SET DEFAULT SOLVER TO HRE ...
            if( ds.sym && ds.positive ) ds.solver=*def_solver_sym_dp;
            else if( ds.sym ) ds.solver=*def_solver_sym;
            else  ds.solver=*def_solver;
            if(mpirank==0 && verbosity>4) cout << "  **Warning: set default solver to " << ds.solver << endl;
        }
        if (nargs[++kk]) ds.eta = GetAny<double>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.minclustersize = GetAny<int>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.mintargetdepth = GetAny<int>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.minsourcedepth = GetAny<int>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.compressor = *GetAny<string*>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.initialclustering = *GetAny<string*>((*nargs[kk])(stack));
        ffassert(++kk == n_name_param);
    }   
}

template<class K>
class HMatrixVirt;

template<class Type, class K>
KNM<K>* HMatrixVirtToDense(HMatrixVirt<K> **Hmat)
{
    ffassert(Hmat && *Hmat);
    HMatrixVirt<K>& H = **Hmat;
    if(std::is_same<Type, KNM<K>>::value) {

        htool::Matrix<K> mdense = H.get_local_dense();
        KNM<K>* M = new KNM<K>( (long) H.nb_rows(), (long) H.nb_cols() );
        for (int i=0; i< H.nb_rows(); i++)
            for (int j=0; j< H.nb_cols(); j++)
                (*M)(i,j) = 0;
        for (int i=0; i< mdense.nb_rows(); i++)
            for (int j=0; j< mdense.nb_cols(); j++)
                (*M)(H.get_permt(i+H.get_local_offset()),H.get_perms(j)) = mdense(i,j);
        return M;
    }
    else {
        std::cerr << " Error of value type HMatrixVirt<K> and KNM<K> " << std::endl;
        ffassert(0);
    }
}

template< class FESPACE>
std::shared_ptr<htool::Cluster<double>> build_clustering(int n, const FESPACE * Uh, const std::vector<double> & p, const Data_Bem_Solver& ds, MPI_Comm comm) {
    htool::ClusterTreeBuilder<double> cluster_builder;
    std::shared_ptr<htool::Cluster<double>> cluster;

    cluster_builder.set_direction_computation_strategy(std::make_shared<htool::ComputeLargestExtent<double>>());
    cluster_builder.set_splitting_strategy(std::make_shared<htool::RegularSplitting<double>>());
    cluster_builder.set_minclustersize(ds.minclustersize);
    // std::shared_ptr<htool::Cluster<htool::PCARegularClustering>> c = std::make_shared<htool::Cluster<htool::PCARegularClustering>>();

    int sizeWorld;
    MPI_Comm_size(comm, &sizeWorld);
    // c->set_minclustersize(ds.minclustersize);
    if (ds.initialclustering == "" || ds.initialclustering == "pca") {
        cluster = std::make_shared<htool::Cluster<double>>(cluster_builder.create_cluster_tree(n,3,p.data(),2,sizeWorld));
    }
    else if (ds.initialclustering == "metis") {
    #if defined(HTOOL_VERSION_GE)
    #if HTOOL_VERSION_GE(0,8,1)
    #if defined(WITH_metis)
    KN<double> parts(n);
    int npart;
    int rank;
    MPI_Comm_size(comm, &npart);
    MPI_Comm_rank(comm, &rank);
    std::vector<int> part(n);

    if (npart > 1) {
        if (rank == 0) {
            partmetis<FESPACE,0,double>(&parts, Uh, npart);
            for (int i=0; i<n; i++)
                part[i] = parts[i];
        }
        MPI_Bcast(part.data(), n, MPI_INT, 0, comm);
    }
    else
        std::fill(part.begin(),part.end(),0);
    cluster = std::make_shared<htool::Cluster<double>>(cluster_builder.create_cluster_tree(n,3,p.data(),2,sizeWorld,part.data()));
    #else
    if (mpirank == 0) std::cerr << "Error: cannot use metis for initial htool clustering ; no metis library" << std::endl;
    ffassert(0);
    #endif
    #else
    if (mpirank == 0) std::cerr << "Error: cannot use metis for initial htool clustering ; need htool 0.8.1 or later" << std::endl;
    ffassert(0);
    #endif
    #else
    if (mpirank == 0) std::cerr << "Error: cannot use metis for initial htool clustering ; need htool 0.8.1 or later" << std::endl;
    ffassert(0);
    #endif
    }
    else {
        if (mpirank == 0) std::cerr << "Error: unknown initial clustering \"" << ds.initialclustering << "\", please use \"pca\" or \"metis\"" << std::endl;
        ffassert(0);
    }
    return cluster;
}

template <class R>
void buildHmat(HMatrixVirt<R>** Hmat, htool::VirtualGenerator<R>* generatorP,const Data_Bem_Solver& data,
                std::shared_ptr<htool::Cluster<double>> t, std::shared_ptr<htool::Cluster<double>> s,
                vector<double> &pt,vector<double> &ps,MPI_Comm comm) {

    int rankWorld;
    MPI_Comm_rank(comm, &rankWorld);
    auto hmatrix_builder = htool::HMatrixTreeBuilder<R, double>(*t, *s, data.epsilon, data.eta, data.sym?'S':'N',data.sym?'U':'N', -1,rankWorld,rankWorld);
    std::shared_ptr<htool::VirtualLowRankGenerator<R,double>> LowRankGenerator = nullptr;
    if (data.compressor=="" || data.compressor == "partialACA")
        LowRankGenerator = std::make_shared<htool::partialACA<R>>();
    else if (data.compressor == "fullACA")
        LowRankGenerator = std::make_shared<htool::fullACA<R>>();
    else if (data.compressor == "SVD")
        LowRankGenerator = std::make_shared<htool::SVD<R>>();
    else {
        cerr << "Error: unknown htool compressor \""+data.compressor+"\"" << endl;
        ffassert(0);
    }

    hmatrix_builder.set_low_rank_generator(LowRankGenerator);
    hmatrix_builder.set_minimal_target_depth(data.mintargetdepth);
    hmatrix_builder.set_minimal_source_depth(data.minsourcedepth);
    hmatrix_builder.set_block_tree_consistency(false);
    HMatrixImpl<R>* test =new HMatrixImpl<R>(*generatorP, t,s,hmatrix_builder,comm);
    *Hmat = test;
}

template<class R, class MMesh, class FESpace1, class FESpace2> 
void creationHMatrixtoBEMForm(const FESpace1 * Uh, const FESpace2 * Vh, const int & VFBEM, 
                             const list<C_F0> & largs, Stack stack, const Data_Bem_Solver &ds, HMatrixVirt<R>** Hmat){


    typedef typename FESpace1::Mesh SMesh;
    typedef typename FESpace2::Mesh TMesh;
    typedef typename SMesh::RdHat SRdHat;
    typedef typename TMesh::RdHat TRdHat;

    typedef typename std::conditional<SMesh::RdHat::d==1,bemtool::Mesh1D,bemtool::Mesh2D>::type MeshBemtool;
    typedef typename std::conditional<SMesh::RdHat::d==1,bemtool::P0_1D,bemtool::P0_2D>::type P0;
    typedef typename std::conditional<SMesh::RdHat::d==1,bemtool::P1_1D,bemtool::P1_2D>::type P1;
    typedef typename std::conditional<SMesh::RdHat::d==1,bemtool::P2_1D,bemtool::P2_2D>::type P2;

    // size of the matrix
    int m=Uh->NbOfDF;
    int n=Vh->NbOfDF;

    /*
    int VFBEM = typeVFBEM(largs,stack);
    if (mpirank == 0 && verbosity>5)
        cout << "test VFBEM type (1 kernel / 2 potential) "  << VFBEM << endl;
    */

    // compression infogfg
    
    MPI_Comm comm = ds.commworld ? *(MPI_Comm*)ds.commworld : MPI_COMM_WORLD;
    
    // source/target meshes
    const SMesh & ThU =Uh->Th; // line
    const TMesh & ThV =Vh->Th; // colunm
    bool samemesh = (void*)&Uh->Th == (void*)&Vh->Th;  // same Fem2D::Mesh     +++ pot or kernel
 
    bemtool::Geometry node; MeshBemtool mesh;
    Mesh2Bemtool(ThU, node, mesh);


    vector<double> pt(3*n);
    vector<double> ps(3*m);
    Fem2D::R3 pp;
    bemtool::R3 p;
    SRdHat pbs;
    TRdHat pbt;
    pbs[0] = 1./(SRdHat::d+1);
    pbs[1] = 1./(SRdHat::d+1);
    if (SRdHat::d == 3) pbs[2] = 1./(SRdHat::d+1);
    pbt[0] = 1./(TRdHat::d+1);
    pbt[1] = 1./(TRdHat::d+1);
    if (TRdHat::d == 3) pbt[2] = 1./(TRdHat::d+1);

    int Snbv = Uh->TFE[0]->ndfonVertex;
    int Snbe = Uh->TFE[0]->ndfonEdge;
    int Snbt = Uh->TFE[0]->ndfonFace;
    bool SP0 = SRdHat::d == 1 ? (Snbv == 0) && (Snbe == 1) && (Snbt == 0) : (Snbv == 0) && (Snbe == 0) && (Snbt == 1);
    bool SP1 = (Snbv == 1) && (Snbe == 0) && (Snbt == 0);
    bool SP2 = (Snbv == 1) && (Snbe == 1) && (Snbt == 0);
    bool SRT0 = (SRdHat::d == 2) && (Snbv == 0) && (Snbe == 1) && (Snbt == 0);

    if (SP2) {
        bemtool::Dof<P2> dof(mesh,true);
        for (int i=0; i<m; i++) {
            const std::vector<bemtool::N2>& jj = dof.ToElt(i);
            p = dof(jj[0][0])[jj[0][1]];
            ps[3*i+0] = p[0];
            ps[3*i+1] = p[1];
            ps[3*i+2] = p[2];
        }
    }
    else if (SRT0) {
        bemtool::Dof<bemtool::RT0_2D> dof(mesh);
        for (int i=0; i<m; i++) {
            const std::vector<bemtool::N2>& jj = dof.ToElt(i);
            p = dof(jj[0][0])[jj[0][1]];
            ps[3*i+0] = p[0];
            ps[3*i+1] = p[1];
            ps[3*i+2] = p[2];
        }
    }
    else
    for (int i=0; i<m; i++) {
        if (SP1)
            pp = ThU.vertices[i];
        else if (SP0)
            pp = ThU[i](pbs);
        else {
            if (mpirank == 0) std::cerr << "ff-BemTool error: only P0, P1 and P2 discretizations are available for now." << std::endl;
            ffassert(0);
        }
        ps[3*i+0] = pp.x;
        ps[3*i+1] = pp.y;
        ps[3*i+2] = pp.z;
    }

    std::shared_ptr<htool::Cluster<double>> t, s;
    s = build_clustering(m, Uh, ps, ds, comm);

    if(!samemesh) {
        if( Vh->TFE[0]->N == 1){
            // case the target FE is scalar
            for (int i=0; i<n; i++) {
                if (Vh->MaxNbNodePerElement == TRdHat::d + 1)
                    pp = ThV.vertices[i];
                else if (Vh->MaxNbNodePerElement == 1)
                    pp = ThV[i](pbt);
                else {
                    if (mpirank == 0) std::cerr << "ff-BemTool error: only P0 and P1 FEspaces are available for reconstructions." << std::endl;
                    ffassert(0);
                }
                pt[3*i+0] = pp.x;
                pt[3*i+1] = pp.y;
                pt[3*i+2] = pp.z;
            }
        }
        else{
            // hack for Maxwell case to have one Hmatrix to avoid one Hmatrix by direction
            ffassert(SRT0 && SRdHat::d == 2 && VFBEM==2);

            // Dans un espace vectoriel, [P1,P1,P1] pour les targets, on a:
            //        n correspond au nombre de dof du FEM space
            // Or dans ce cas, on veut que n = mesh_Target.nv
            // 
            // ==> on n'a pas besoin de resize les points p2
            int nnn= Vh->TFE[0]->N; // the size of the vector FESpace. For [P1,P1,P1], nnn=3;

            int nDofScalar = n/nnn; // computation of the dof of one component

            for (int i=0; i<nDofScalar; i++) {
                if (Vh->MaxNbNodePerElement == TRdHat::d + 1)
                    pp = ThV.vertices[i];
                else if (Vh->MaxNbNodePerElement == 1)
                    pp = ThV[i](pbt);
                else {
                    if (mpirank == 0) std::cerr << "ff-BemTool error: only P0 and P1 FEspaces are available for reconstructions." << std::endl;
                    ffassert(0);
                }

                for(int iii=0; iii<nnn; iii++){
                    ffassert( nnn*3*i+3*iii+2 < nnn*3*n );
                    pt[nnn*3*i+3*iii+0] = pp.x;
                    pt[nnn*3*i+3*iii+1] = pp.y;
                    pt[nnn*3*i+3*iii+2] = pp.z;
                }
            }
        }
        t = build_clustering(n, Vh, pt, ds, comm);
    }
    else{
        pt=ps;
        t=s;
    }

    // creation of the generator for htool and creation the matrix
    if (VFBEM==1) {
        // info kernel
        pair<BemKernel*,Complex> kernel = getBemKernel(stack,largs);
        BemKernel *Ker = kernel.first;
        Complex alpha = kernel.second;

        // check for Maxwell case
        if (SRT0 && SRdHat::d == 2) {
            if( Ker->typeKernel[1] >0 ){
                cerr << "vector BEM Maxwell is not valid for combined field formulation."<< endl;
                ffassert(0);
            }
            // BemKernel->typeKernel[0] == 6 :: MA_SL
            ffassert( Ker->typeKernel[0] == 6 ); // check MA_SL
        }

        if(mpirank == 0 && verbosity >5) {
            int nk=-1;
            iscombinedKernel(Ker) ? nk=2 : nk=1;
            //for(int i=0;i<nk;i++)
            //    cout << " kernel info... i: " << i << " typeKernel: " << Ker->typeKernel[i] << " wave number: " << Ker->wavenum[i]  << " coeffcombi: " << Ker->coeffcombi[i] <<endl;
        }
        htool::VirtualGenerator<R>* generator;

        if(mpirank == 0 && verbosity>5)
            std::cout << "Creating dof" << std::endl;

        if (SP0) {
            bemtool::Dof<P0> dof(mesh);
            ff_BIO_Generator<R,P0,SMesh>(generator,Ker,dof,alpha);
        }
        else if (SP1) {
            bemtool::Dof<P1> dof(mesh,true);
            ff_BIO_Generator<R,P1,SMesh>(generator,Ker,dof,alpha);
        }
        else if (SP2) {
            bemtool::Dof<P2> dof(mesh,true);
            ff_BIO_Generator<R,P2,SMesh>(generator,Ker,dof,alpha);
        }
        else if (SRT0 && SRdHat::d == 2) {
            // BemKernel->typeKernel[0] == 6 :: MA_SL
            bemtool::Dof<bemtool::RT0_2D> dof(mesh);
            ff_BIO_Generator_Maxwell<R>(generator,Ker,dof,alpha);
        }
        else
            ffassert(0);

        // build the Hmat
        buildHmat(Hmat, generator, ds, t, s, pt, ps, comm);

        delete generator;
    }
    else if (VFBEM==2) {
        BemPotential *Pot = getBemPotential(stack,largs);
        bemtool::Geometry node_output;

        if (SRT0 && SRdHat::d == 2) {
            // check if we have the good FE space for vector BEM
            ffassert( Pot->typePotential == 4 ); // check MA_SL
        }

        if (Vh->MaxNbNodePerElement == TRdHat::d + 1)
            Mesh2Bemtool(ThV,node_output);
        else if (Vh->MaxNbNodePerElement == 1) {
            for (int i=0; i<n; i++) {
                pp = ThV[i](pbt);
                p[0]=pp.x; p[1]=pp.y; p[2]=pp.z;
                node_output.setnodes(p);
            }
        }
        else
            ffassert(0);

        htool::VirtualGenerator<R>* generator;
        if (SP0) {
            bemtool::Dof<P0> dof(mesh);
            ff_POT_Generator<R,P0,MeshBemtool,SMesh>(generator,Pot,dof,mesh,node_output);
        }
        else if (SP1) {
            bemtool::Dof<P1> dof(mesh,true);
            ff_POT_Generator<R,P1,MeshBemtool,SMesh>(generator,Pot,dof,mesh,node_output);
        }
        else if (SP2) {
            bemtool::Dof<P2> dof(mesh,true);
            ff_POT_Generator<R,P2,MeshBemtool,SMesh>(generator,Pot,dof,mesh,node_output);
        }
        else if (SRT0 && SRdHat::d == 2) {
            bemtool::Dof<bemtool::RT0_2D> dof(mesh);
            ff_POT_Generator_Maxwell<R,bemtool::RT0_2D>(generator,Pot,dof,mesh,node_output);
        }
        else
            ffassert(0);

        buildHmat(Hmat, generator, ds, t, s, pt, ps, comm);

        delete generator;
    }
}

template<> void creationHMatrixtoBEMForm<double, MeshS, FESpaceS, FESpaceS>(const FESpaceS * Uh, const FESpaceS * Vh, const int & VFBEM, 
                             const std::list<C_F0> & largs, Stack stack, const Data_Bem_Solver &ds, HMatrixVirt<double> **Hmat){
                                 cerr << "we can't use bemtool with Real type." << endl;
                                 ffassert(0);
                             }

template<> void creationHMatrixtoBEMForm<double, MeshL, FESpaceL, FESpaceL>(const FESpaceL * Uh, const FESpaceL * Vh, const int & VFBEM, 
                             const std::list<C_F0> & largs, Stack stack, const Data_Bem_Solver &ds, HMatrixVirt<double> **Hmat){
                                 cerr << "we can't use bemtool with Real type." << endl;
                                 ffassert(0);
                             }

#endif
;
