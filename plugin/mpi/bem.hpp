#ifndef BEM_HPP_
#define BEM_HPP_

double ff_htoolEta=10., ff_htoolEpsilon=1e-3;
long ff_htoolMinclustersize=10, ff_htoolMaxblocksize=1000000, ff_htoolMintargetdepth=0, ff_htoolMinsourcedepth=0;

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
    virtual const std::vector<htool::SubMatrix<K>*>& get_MyNearFieldMats() const = 0;
    virtual const htool::LowRankMatrix<K>& get_MyFarFieldMats(int i) const = 0;
    virtual int get_MyFarFieldMats_size() const = 0;
    virtual const std::vector<htool::SubMatrix<K>*>& get_MyStrictlyDiagNearFieldMats() const = 0;
    virtual htool::Matrix<K> get_local_dense() const = 0;
    virtual int get_permt(int) const = 0;
    virtual int get_perms(int) const = 0;
    virtual char get_symmetry_type() const = 0;
    virtual void set_compression(std::shared_ptr<htool::VirtualLowRankGenerator<K>> compressor) = 0;

    // Getters/setters for parameters
    virtual double get_epsilon() const                            = 0;
    virtual double get_eta() const                                = 0;
    virtual int get_minsourcedepth() const                        = 0;
    virtual int get_mintargetdepth() const                        = 0;
    virtual int get_maxblocksize() const                          = 0;
    virtual void set_epsilon(double epsilon0)                     = 0;
    virtual void set_eta(double eta0)                             = 0;
    virtual void set_minsourcedepth(unsigned int minsourcedepth0) = 0;
    virtual void set_mintargetdepth(unsigned int mintargetdepth0) = 0;
    virtual void set_maxblocksize(unsigned int maxblocksize)      = 0;

    // Build
    virtual void build(htool::VirtualGenerator<K> &mat, double* xt, double* xs) = 0;
    virtual void build(htool::VirtualGenerator<K> &mat, double* xt)             = 0;

    virtual ~HMatrixVirt() {};
};


template<class K>
class HMatrixImpl : public HMatrixVirt<K> {
public:
    htool::HMatrix<K> H;
public:
    HMatrixImpl(std::shared_ptr<htool::VirtualCluster> t, std::shared_ptr<htool::VirtualCluster> s, double epsilon=1e-6 ,double eta=10, char symmetry='N', char UPLO='N',const int& reqrank=-1, MPI_Comm comm=MPI_COMM_WORLD) : H(t,s,epsilon,eta,symmetry,UPLO,reqrank,comm){}
    const std::map<std::string, std::string>& get_infos() const {return H.get_infos();}
    void mvprod_global(const K* const in, K* const out,const int& mu=1) const {return H.mvprod_global_to_global(in,out,mu);}
    int nb_rows() const { return H.nb_rows();}
    int nb_cols() const { return H.nb_cols();}
    void cluster_to_target_permutation(const K* const in, K* const out) const {return H.cluster_to_target_permutation(in,out);}
    void source_to_cluster_permutation(const K* const in, K* const out) const {return H.source_to_cluster_permutation(in,out);}
    MPI_Comm get_comm() const {return H.get_comm();}
    int get_rankworld() const {return H.get_rankworld();}
    int get_sizeworld() const {return H.get_sizeworld();}
    int get_local_size() const {return H.get_local_size();}
    int get_local_offset() const {return H.get_local_offset();}
    const std::vector<htool::SubMatrix<K>*>& get_MyNearFieldMats() const {return H.get_MyNearFieldMats();}
    const htool::LowRankMatrix<K>& get_MyFarFieldMats(int i) const {return *(H.get_MyFarFieldMats()[i]);}
    int get_MyFarFieldMats_size() const {return H.get_MyFarFieldMats().size();}
    const std::vector<htool::SubMatrix<K>*>& get_MyStrictlyDiagNearFieldMats() const {return H.get_MyStrictlyDiagNearFieldMats();}
    htool::Matrix<K> get_local_dense() const {return H.get_local_dense();}
    int get_permt(int i) const {return H.get_permt(i);}
    int get_perms(int i) const {return H.get_perms(i);}
    char get_symmetry_type() const {return H.get_symmetry_type();}
    void set_compression(std::shared_ptr<htool::VirtualLowRankGenerator<K>> compressor) {H.set_compression(compressor);};

    // Getters/setters for parameters
    double get_epsilon() const  {return H.get_epsilon();}                          
    double get_eta() const    {return H.get_epsilon();}                            
    int get_minsourcedepth() const {return H.get_minsourcedepth();}              
    int get_mintargetdepth() const {return H.get_mintargetdepth();}            
    int get_maxblocksize() const   {return H.get_maxblocksize();}              
    void set_epsilon(double epsilon0) {H.set_epsilon(epsilon0);}             
    void set_eta(double eta0)        {H.set_eta(eta0);}                        
    void set_minsourcedepth(unsigned int minsourcedepth0) {H.set_minsourcedepth(minsourcedepth0);}   
    void set_mintargetdepth(unsigned int mintargetdepth0) {H.set_mintargetdepth(mintargetdepth0);}   
    void set_maxblocksize(unsigned int maxblocksize0)   {H.set_maxblocksize(maxblocksize0);}      

    // Build
    void build(htool::VirtualGenerator<K> &mat, const std::vector<htool::R3> &xt, const std::vector<double> &rt, const std::vector<int> &tabt, const std::vector<double> &gt, const std::vector<htool::R3> &xs, const std::vector<double> &rs, const std::vector<int> &tabs, const std::vector<double> &gs) {H.build(mat,xt,tabt,xs,tabs);}


    void build(htool::VirtualGenerator<K> &mat, double* xt, double* xs) {H.build(mat,xt,xs);}
    void build(htool::VirtualGenerator<K> &mat, double* xt) {H.build(mat,xt);}                       
};


struct Data_Bem_Solver
: public Data_Sparse_Solver {
    double eta;
    int minclustersize,maxblocksize,mintargetdepth,minsourcedepth;
    string compressor;
       
    Data_Bem_Solver()
       : Data_Sparse_Solver(),
       eta(ff_htoolEta),
       minclustersize(ff_htoolMinclustersize),
       maxblocksize(ff_htoolMaxblocksize),
       mintargetdepth(ff_htoolMintargetdepth),
       minsourcedepth(ff_htoolMinsourcedepth),
       compressor("partialACA")
    
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
        if (nargs[++kk]) ds.maxblocksize = GetAny<int>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.mintargetdepth = GetAny<int>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.minsourcedepth = GetAny<int>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.compressor = *GetAny<string*>((*nargs[kk])(stack));
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
        KNM<K>* M = new KNM<K>( (long) H.nb_cols(), (long) H.nb_rows() );
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


template <class R>
void builHmat(HMatrixVirt<R>** Hmat, htool::VirtualGenerator<R>* generatorP,const Data_Bem_Solver& data,vector<double> &p1,vector<double> &p2,MPI_Comm comm)  {
    std::shared_ptr<htool::Cluster<htool::PCARegularClustering>> t =std::make_shared<htool::Cluster<htool::PCARegularClustering>>();
    std::shared_ptr<htool::Cluster<htool::PCARegularClustering>> s =std::make_shared<htool::Cluster<htool::PCARegularClustering>>();
    t->set_minclustersize(data.minclustersize);
    s->set_minclustersize(data.minclustersize);
    t->build(generatorP->nb_rows(),p2.data(),2,comm);
    s->build(generatorP->nb_cols(),p1.data(),2,comm);

    *Hmat = new HMatrixImpl<R>(t,s,data.epsilon,data.eta,data.sym?'S':'N',data.sym?'U':'N',-1,comm);
    std::shared_ptr<htool::VirtualLowRankGenerator<R>> LowRankGenerator = nullptr;
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

    (*Hmat)->set_compression(LowRankGenerator);
    (*Hmat)->set_maxblocksize(data.maxblocksize);
    (*Hmat)->set_mintargetdepth(data.mintargetdepth);
    (*Hmat)->set_minsourcedepth(data.minsourcedepth);
    // TODO set options
    (*Hmat)->build(*generatorP,p2.data(),p1.data());

}

template<class Matrix, class R, typename std::enable_if< std::is_same< Matrix, HMatrixVirt<R>* >::value >::type* = nullptr>
void Assembly(Matrix* A, htool::VirtualGenerator<R>* generator, const Data_Bem_Solver& data,vector<double> &p1,vector<double> &p2,MPI_Comm comm,int dim,bool = false) {
    builHmat<R>(A,generator,data,p1,p2,comm);
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
    int n=Uh->NbOfDF;
    int m=Vh->NbOfDF;

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


    vector<double> p1(3*n);
    vector<double> p2(3*m);
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
        for (int i=0; i<n; i++) {
            const std::vector<bemtool::N2>& jj = dof.ToElt(i);
            p = dof(jj[0][0])[jj[0][1]];
            p1[3*i+0] = p[0];
            p1[3*i+1] = p[1];
            p1[3*i+2] = p[2];
        }
    }
    else if (SRT0) {
        bemtool::Dof<bemtool::RT0_2D> dof(mesh);
        for (int i=0; i<n; i++) {
            const std::vector<bemtool::N2>& jj = dof.ToElt(i);
            p = dof(jj[0][0])[jj[0][1]];
            p1[3*i+0] = p[0];
            p1[3*i+1] = p[1];
            p1[3*i+2] = p[2];
        }
    }
    else
    for (int i=0; i<n; i++) {
        if (SP1)
            pp = ThU.vertices[i];
        else if (SP0)
            pp = ThU[i](pbs);
        else {
            if (mpirank == 0) std::cerr << "ff-BemTool error: only P0, P1 and P2 discretizations are available for now." << std::endl;
            ffassert(0);
        }
        p1[3*i+0] = pp.x;
        p1[3*i+1] = pp.y;
        p1[3*i+2] = pp.z;
    }

    std::shared_ptr<htool::Cluster<htool::PCARegularClustering>> t =std::make_shared<htool::Cluster<htool::PCARegularClustering>>();
    std::shared_ptr<htool::Cluster<htool::PCARegularClustering>> s =std::make_shared<htool::Cluster<htool::PCARegularClustering>>();
    t->set_minclustersize(ds.minclustersize);
    s->set_minclustersize(ds.minclustersize);
    t->build(n,p1.data(),2,comm);

    if(!samemesh) {
        if( Vh->TFE[0]->N == 1){
            // case the target FE is scalar
            for (int i=0; i<m; i++) {
                if (Vh->MaxNbNodePerElement == TRdHat::d + 1)
                    pp = ThV.vertices[i];
                else if (Vh->MaxNbNodePerElement == 1)
                    pp = ThV[i](pbt);
                else {
                    if (mpirank == 0) std::cerr << "ff-BemTool error: only P0 and P1 FEspaces are available for reconstructions." << std::endl;
                    ffassert(0);
                }
                p2[3*i+0] = pp.x;
                p2[3*i+1] = pp.y;
                p2[3*i+2] = pp.z;
            }
        }
        else{
            // hack for Maxwell case to have one Hmatrix to avoid one Hmatrix by direction
            ffassert(SRT0 && SRdHat::d == 2 && VFBEM==2);

            // Dans un espace verctoriel, [P1,P1,P1] pour les targets, on a:
            //        m correspond au nombre de dof du FEM space
            // Or dans ce cas, on veut que m = mesh_Target.nv 
            // 
            // ==> on n'a pas besoin de resize les points p2
            int nnn= Vh->TFE[0]->N; // the size of the vector FESpace. For [P1,P1,P1], nnn=3;

            int mDofScalar = m/nnn; // computation of the dof of one component 

            for (int i=0; i<mDofScalar; i++) {
                if (Vh->MaxNbNodePerElement == TRdHat::d + 1)
                    pp = ThV.vertices[i];
                else if (Vh->MaxNbNodePerElement == 1)
                    pp = ThV[i](pbt);
                else {
                    if (mpirank == 0) std::cerr << "ff-BemTool error: only P0 and P1 FEspaces are available for reconstructions." << std::endl;
                    ffassert(0);
                }

                for(int iii=0; iii<nnn; iii++){
                    ffassert( nnn*3*i+3*iii+2 < nnn*3*m );
                    p2[nnn*3*i+3*iii+0] = pp.x;
                    p2[nnn*3*i+3*iii+1] = pp.y;
                    p2[nnn*3*i+3*iii+2] = pp.z;
                }
            }
        }
        s->build( m,p2.data(),2,comm);  
    }
    else{
        p2=p1;
        s=t;
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
        *Hmat = new HMatrixImpl<R>(t,s,ds.epsilon,ds.eta,ds.sym?'S':'N',ds.sym?'U':'N',-1,comm);
        std::shared_ptr<htool::VirtualLowRankGenerator<R>> compressor = nullptr;
        if ( ds.compressor == "" || ds.compressor == "partialACA")
            compressor = std::make_shared<htool::partialACA<R>>();
        else if (ds.compressor == "fullACA")
            compressor = std::make_shared<htool::fullACA<R>>();
        else if (ds.compressor == "SVD")
            compressor = std::make_shared<htool::SVD<R>>();
        else {
            cerr << "Error: unknown htool compressor \""+ds.compressor+"\"" << endl;
            ffassert(0);
        }

        (*Hmat)->set_compression(compressor);
        (*Hmat)->set_maxblocksize(ds.maxblocksize);
        (*Hmat)->set_mintargetdepth(ds.mintargetdepth);
        (*Hmat)->set_minsourcedepth(ds.minsourcedepth);
        (*Hmat)->build(*generator,p1.data());


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
            for (int i=0; i<m; i++) {
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
        Assembly(Hmat,generator,ds,p1,p2,comm,MMesh::RdHat::d+1);
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
