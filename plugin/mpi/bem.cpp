//ff-c++-LIBRARY-dep: cxx11 [mkl|blas] [petsccomplex] mpi pthread htool bemtool boost
//ff-c++-cpp-dep:
// for def  M_PI under windows in <cmath>
#define _USE_MATH_DEFINES
#include <ff++.hpp>
#include <AFunction_ext.hpp>
#include <lgfem.hpp>
#include <R3.hpp>

#include <htool/clustering/DDM_cluster.hpp>
#include <htool/lrmat/partialACA.hpp>
#include <htool/lrmat/fullACA.hpp>
#include <htool/lrmat/SVD.hpp>
#include <htool/types/matrix.hpp>
#include <htool/types/hmatrix.hpp>

// include the bemtool library .... path define in where library
//#include <bemtool/operator/block_op.hpp>
#include <bemtool/tools.hpp>
#include <bemtool/fem/dof.hpp>
#include <bemtool/operator/operator.hpp>
#include <bemtool/miscellaneous/htool_wrap.hpp>
#include "PlotStream.hpp"

#ifdef WITH_petsccomplex
#include <petsc.h>
#include "PETSc.hpp"
typedef PETSc::DistributedCSR< HpSchwarz< PetscScalar > > Dmat;
#endif

#include "common.hpp"

extern FILE *ThePlotStream;

using namespace std;
using namespace htool;
using namespace bemtool;

class BemKernel;
class BemPotential;
//class FoperatorKBEM;
// new type for bem
typedef const BemKernel *pBemKernel;
typedef const BemPotential *pBemPotential;

typedef  const BemKernel fkernel;
typedef  const BemPotential fpotential;

template< class K >
AnyType AddIncrement(Stack stack, const AnyType &a) {
    K m = GetAny< K >(a);
    m->increment( );
    Add2StackOfPtr2FreeRC(stack, m);
    if (mpirank==0 && verbosity > 1) cout << "AddIncrement:: increment + Add2StackOfPtr2FreeRC " << endl;
    return a;
}


template<class K>
class HMatrixVirt {
public:
    virtual const std::map<std::string, std::string>& get_infos() const = 0;
    virtual void mvprod_global(const K* const in, K* const out,const int& mu=1) const = 0;
    virtual int nb_rows() const = 0;
    virtual int nb_cols() const = 0;
    virtual void cluster_to_target_permutation(const K* const in, K* const out) const = 0;
    virtual void source_to_cluster_permutation(const K* const in, K* const out) const = 0;
    virtual const MPI_Comm& get_comm() const = 0;
    virtual int get_rankworld() const = 0;
    virtual int get_sizeworld() const = 0;
    virtual int get_local_size() const = 0;
    virtual const std::vector<SubMatrix<K>*>& get_MyNearFieldMats() const = 0;
    virtual const LowRankMatrix<K,GeometricClusteringDDM>& get_MyFarFieldMats(int i) const = 0;
    virtual int get_MyFarFieldMats_size() const = 0;
    virtual const std::vector<SubMatrix<K>*>& get_MyStrictlyDiagNearFieldMats() const = 0;
    virtual Matrix<K> to_dense_perm() const = 0;
    virtual Matrix<K> to_local_dense() const = 0;

    virtual ~HMatrixVirt() {};
};


template<class K, template<class,class> class LR>
class HMatrixImpl : public HMatrixVirt<K> {
private:
    HMatrix<K,LR,GeometricClusteringDDM> H;
public:
    HMatrixImpl(IMatrix<K>& I, const std::vector<htool::R3>& xt, bool symmetry=false,const int& reqrank=-1, MPI_Comm comm=MPI_COMM_WORLD) : H(I,xt,symmetry,reqrank,comm){}
    HMatrixImpl(IMatrix<K>& I, const std::vector<htool::R3>& xt, const std::vector<htool::R3>& xs, const int& reqrank=-1, MPI_Comm comm=MPI_COMM_WORLD) : H(I,xt,xs,reqrank,comm){}
    const std::map<std::string, std::string>& get_infos() const {return H.get_infos();}
    void mvprod_global(const K* const in, K* const out,const int& mu=1) const {return H.mvprod_global(in,out,mu);}
    int nb_rows() const { return H.nb_rows();}
    int nb_cols() const { return H.nb_cols();}
    void cluster_to_target_permutation(const K* const in, K* const out) const {return H.cluster_to_target_permutation(in,out);}
    void source_to_cluster_permutation(const K* const in, K* const out) const {return H.source_to_cluster_permutation(in,out);}
    const MPI_Comm& get_comm() const {return H.get_comm();}
    int get_rankworld() const {return H.get_rankworld();}
    int get_sizeworld() const {return H.get_sizeworld();}
    int get_local_size() const {return H.get_local_size();}
    const std::vector<SubMatrix<K>*>& get_MyNearFieldMats() const {return H.get_MyNearFieldMats();}
    const LowRankMatrix<K,GeometricClusteringDDM>& get_MyFarFieldMats(int i) const {return *(H.get_MyFarFieldMats()[i]);}
    int get_MyFarFieldMats_size() const {return H.get_MyFarFieldMats().size();}
    const std::vector<SubMatrix<K>*>& get_MyStrictlyDiagNearFieldMats() const {return H.get_MyStrictlyDiagNearFieldMats();}
    Matrix<K> to_dense_perm() const {return H.to_dense_perm();}
    Matrix<K> to_local_dense() const {return H.to_local_dense();}
};



template<class ffmesh, class bemtoolmesh>
void Mesh2Bemtool(const ffmesh &Th, Geometry &node, bemtoolmesh &mesh ) {
    
    typedef typename ffmesh::RdHat RdHat;
    typedef typename ffmesh::Element E;
    const int dHat =  RdHat::d;
    
    // create the geometry;
    
    bemtool::R3 p;
    for(int iv=0 ; iv<Th.nv ; iv++){
        p[0]=Th.vertices[iv].x;p[1]=Th.vertices[iv].y;p[2]=Th.vertices[iv].z;
        node.setnodes(p);
    }
    
    node.initEltData();
    
    if(mpirank==0 && verbosity>10) std::cout << "Creating mesh domain (nodes)" << std::endl;
    
    mesh.set_elt(node);
    bemtool::array<dHat+1,int> I;
    if(mpirank==0 && verbosity>10) std::cout << "End creating mesh domain mesh" << std::endl;
    
    if(mpirank==0 && verbosity>10) std::cout << "Creating geometry domain (elements)" << std::endl;
    for(int it=0; it<Th.nt; it++){
        const E &K(Th[it]);
        for(int j=0;j<dHat+1;j++)
            I[j]=Th.operator () (K[j]);
        mesh.setOneElt(node,I);
    }
    
    //mesh = unbounded;
    //Orienting(mesh);
    Normal<dHat> N(mesh);
    for(int it=0; it<Th.nt; it++){
        const E &K(Th[it]);
        Fem2D::R3 nn = K.NFrenetUnitaire();
        bemtool::R3 mm; mm[0]=nn.x; mm[1]=nn.y; mm[2]=nn.z;
        N.set(it, mm);
    }
    mesh.Orienting(N);
    
    if(mpirank==0 && verbosity>10) std::cout << "end creating geometry domain" << std::endl;
}

template<class ffmesh>
void Mesh2Bemtool(const ffmesh &Th, Geometry &node) {
    if(mpirank==0 && verbosity>10) std::cout << "Creating mesh output" << std::endl;
    bemtool::R3 p;
    Fem2D::R3 pp;
    for(int iv=0 ; iv<Th.nv ; iv++){
        pp = Th.vertices[iv];
        p[0]=pp.x; p[1]=pp.y; p[2]=pp.z;
        node.setnodes(p);
    }
}

static PetscErrorCode s2c(PetscContainer ctx, PetscScalar* in, PetscScalar* out) {
    HMatrixVirt<PetscScalar>** Hmat;

    PetscFunctionBegin;
    PetscContainerGetPointer(ctx, (void**)&Hmat);
    (*Hmat)->source_to_cluster_permutation(in, out);
    PetscFunctionReturn(0);
}

static PetscErrorCode c2s(PetscContainer ctx, PetscScalar* in, PetscScalar* out) {
    HMatrixVirt<PetscScalar>** Hmat;

    PetscFunctionBegin;
    PetscContainerGetPointer(ctx, (void**)&Hmat);
    (*Hmat)->cluster_to_target_permutation(in, out);
    PetscFunctionReturn(0);
}

template<class Type, class K>
AnyType To(Stack stack,Expression emat,Expression einter,int init)
{
    ffassert(einter);
    HMatrixVirt<K>** Hmat = GetAny<HMatrixVirt<K>** >((*einter)(stack));
    ffassert(Hmat && *Hmat);
    HMatrixVirt<K>& H = **Hmat;
    if(std::is_same<Type, KNM<K>>::value) {
        Matrix<K> mdense = H.to_dense_perm();
        const std::vector<K>& vdense = mdense.get_mat();
        KNM<K>* M = GetAny<KNM<K>*>((*emat)(stack));
        for (int i=0; i< mdense.nb_rows(); i++)
            for (int j=0; j< mdense.nb_cols(); j++)
                (*M)(i,j) = mdense(i,j);
        return M;
    }
    else {
#ifndef WITH_petsccomplex
        ffassert(0);
#else
        Matrix<K> mdense = H.to_local_dense();
        Dmat* dense = GetAny<Dmat*>((*emat)(stack));
        dense->dtor();
        MatCreate(H.get_comm(), &dense->_petsc);
        MatSetType(dense->_petsc, MATMPIDENSE);
        MatSetSizes(dense->_petsc, H.get_local_size(), PETSC_DECIDE, H.nb_rows(), H.nb_cols());
        MatMPIDenseSetPreallocation(dense->_petsc, PETSC_NULL);
        PetscScalar* array;
        MatDenseGetArray(dense->_petsc, &array);
        if (array) {
           for (int i = 0; i < mdense.nb_rows(); ++i)
               for (int j = 0; j < mdense.nb_cols(); ++j)
                   array[i + j * mdense.nb_rows()] = mdense(i,j);
        }
        MatDenseRestoreArray(dense->_petsc, &array);
        MatAssemblyBegin(dense->_petsc, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(dense->_petsc, MAT_FINAL_ASSEMBLY);
        Mat M;
        MatCreate(H.get_comm(), &M);
        MatSetSizes(M, H.get_local_size(), H.get_local_size(), H.nb_rows(), H.nb_cols());
        MatSetType(M, MATMPIAIJ);
        MatMPIAIJSetPreallocation(M, H.get_local_size(), NULL, H.nb_cols() - H.get_local_size(), NULL);
        MatSetUp(M);

        MatSetOption(M, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
        MatSetOption(M, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);
        MatConvert(dense->_petsc, MATMPIAIJ, MAT_REUSE_MATRIX, &M);
        MatHeaderReplace(dense->_petsc, &M);
        PetscContainer ptr;
        PetscContainerCreate(H.get_comm(), &ptr);
        PetscContainerSetPointer(ptr, Hmat);
        PetscObjectComposeFunction((PetscObject)ptr, "s2c_C", s2c);
        PetscObjectComposeFunction((PetscObject)ptr, "c2s_C", c2s);
        PetscObjectCompose((PetscObject)dense->_petsc, "Hmat", (PetscObject)ptr);
        PetscContainerDestroy(&ptr);
        return dense;
#endif
    }
}

template<class Type, class K, int init>
AnyType To(Stack stack,Expression emat,Expression einter)
{ return To<Type, K>(stack,emat,einter,init);}

template<class V, class K>
class Prod {
public:
    const HMatrixVirt<K>* h;
    const V u;
    Prod(HMatrixVirt<K>** v, V w) : h(*v), u(w) {}
    
    void prod(V x) const {h->mvprod_global(*(this->u), *x);};
    
    static V mv(V Ax, Prod<V, K> A) {
        *Ax = K();
        A.prod(Ax);
        return Ax;
    }
    static V init(V Ax, Prod<V, K> A) {
        Ax->init(A.u->n);
        return mv(Ax, A);
    }
    
};


// post treatment for HMatrix

template<class K>
std::map<std::string, std::string>* get_infos(HMatrixVirt<K>** const& H) {
    return new std::map<std::string, std::string>((*H)->get_infos());
}

string* get_info(std::map<std::string, std::string>* const& infos, string* const& s){
    return new string((*infos)[*s]);
}

ostream & operator << (ostream &out, const std::map<std::string, std::string> & infos)
{
    for (std::map<std::string,std::string>::const_iterator it = infos.begin() ; it != infos.end() ; ++it){
        out<<it->first<<"\t"<<it->second<<std::endl;
    }
    out << std::endl;
    return out;
}

template<class A>
struct PrintPinfos: public binary_function<ostream*,A,ostream*> {
    static ostream* f(ostream* const  & a,const A & b)  {  *a << *b;
        return a;
    }
};

template<class K>
class plotHMatrix : public OneOperator {
public:
    
    class Op : public E_F0info {
    public:
        Expression a;
        
        static const int n_name_param = 2;
        static basicAC_F0::name_and_type name_param[] ;
        Expression nargs[n_name_param];
        bool arg(int i,Stack stack,bool a) const{ return nargs[i] ? GetAny<bool>( (*nargs[i])(stack) ): a;}
        long argl(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
        
    public:
        Op(const basicAC_F0 &  args,Expression aa) : a(aa) {
            args.SetNameParam(n_name_param,name_param,nargs);
        }
        
        AnyType operator()(Stack stack) const{
            
            bool wait = arg(0,stack,false);
            long dim = argl(1,stack,2);
            
            HMatrixVirt<K>** H =GetAny<HMatrixVirt<K>** >((*a)(stack));
            
            PlotStream theplot(ThePlotStream);
            
            if (mpirank == 0 && ThePlotStream) {
                theplot.SendNewPlot();
                theplot << 3L;
                theplot <= wait;
                theplot << 17L;
                theplot <= dim;
                theplot << 4L;
                theplot <= true;
                theplot.SendEndArgPlot();
                theplot.SendMeshes();
                theplot << 0L;
                theplot.SendPlots();
                theplot << 1L;
                theplot << 31L;
            }
            
            if (!H || !(*H)) {
                if (mpirank == 0&& ThePlotStream) {
                    theplot << 0;
                    theplot << 0;
                    theplot << 0L;
                    theplot << 0L;
                }
            }
            else {
                const std::vector<SubMatrix<K>*>& dmats = (*H)->get_MyNearFieldMats();
                
                int nbdense = dmats.size();
                int nblr = (*H)->get_MyFarFieldMats_size();
                
                int sizeworld = (*H)->get_sizeworld();
                int rankworld = (*H)->get_rankworld();
                
                int nbdenseworld[sizeworld];
                int nblrworld[sizeworld];
                MPI_Allgather(&nbdense, 1, MPI_INT, nbdenseworld, 1, MPI_INT, (*H)->get_comm());
                MPI_Allgather(&nblr, 1, MPI_INT, nblrworld, 1, MPI_INT, (*H)->get_comm());
                int nbdenseg = 0;
                int nblrg = 0;
                for (int i=0; i<sizeworld; i++) {
                    nbdenseg += nbdenseworld[i];
                    nblrg += nblrworld[i];
                }
                
                int* buf = new int[4*(mpirank==0?nbdenseg:nbdense) + 5*(mpirank==0?nblrg:nblr)];
                
                for (int i=0;i<nbdense;i++) {
                    const SubMatrix<K>& l = *(dmats[i]);
                    buf[4*i] = l.get_offset_i();
                    buf[4*i+1] = l.get_offset_j();
                    buf[4*i+2] = l.nb_rows();
                    buf[4*i+3] = l.nb_cols();
                }
                
                int displs[sizeworld];
                int recvcounts[sizeworld];
                displs[0] = 0;
                
                for (int i=0; i<sizeworld; i++) {
                    recvcounts[i] = 4*nbdenseworld[i];
                    if (i > 0)    displs[i] = displs[i-1] + recvcounts[i-1];
                }
                MPI_Gatherv(rankworld==0?MPI_IN_PLACE:buf, recvcounts[rankworld], MPI_INT, buf, recvcounts, displs, MPI_INT, 0, (*H)->get_comm());
                
                int* buflr = buf + 4*(mpirank==0?nbdenseg:nbdense);
                double* bufcomp = new double[mpirank==0?nblrg:nblr];
                
                for (int i=0;i<nblr;i++) {
                    const LowRankMatrix<K,GeometricClusteringDDM>& l = (*H)->get_MyFarFieldMats(i);
                    buflr[5*i] = l.get_offset_i();
                    buflr[5*i+1] = l.get_offset_j();
                    buflr[5*i+2] = l.nb_rows();
                    buflr[5*i+3] = l.nb_cols();
                    buflr[5*i+4] = l.rank_of();
                    bufcomp[i] = l.compression();
                }
                
                for (int i=0; i<sizeworld; i++) {
                    recvcounts[i] = 5*nblrworld[i];
                    if (i > 0)    displs[i] = displs[i-1] + recvcounts[i-1];
                }
                
                MPI_Gatherv(rankworld==0?MPI_IN_PLACE:buflr, recvcounts[rankworld], MPI_INT, buflr, recvcounts, displs, MPI_INT, 0, (*H)->get_comm());
                
                for (int i=0; i<sizeworld; i++) {
                    recvcounts[i] = nblrworld[i];
                    if (i > 0)    displs[i] = displs[i-1] + recvcounts[i-1];
                }
                
                MPI_Gatherv(rankworld==0?MPI_IN_PLACE:bufcomp, recvcounts[rankworld], MPI_DOUBLE, bufcomp, recvcounts, displs, MPI_DOUBLE, 0, (*H)->get_comm());
                
                if (mpirank == 0 && ThePlotStream ) {
                    
                    int si = (*H)->nb_rows();
                    int sj = (*H)->nb_cols();
                    
                    theplot << si;
                    theplot << sj;
                    theplot << (long)nbdenseg;
                    theplot << (long)nblrg;
                    
                    for (int i=0;i<nbdenseg;i++) {
                        theplot << buf[4*i];
                        theplot << buf[4*i+1];
                        theplot << buf[4*i+2];
                        theplot << buf[4*i+3];
                    }
                    
                    for (int i=0;i<nblrg;i++) {
                        theplot << buflr[5*i];
                        theplot << buflr[5*i+1];
                        theplot << buflr[5*i+2];
                        theplot << buflr[5*i+3];
                        theplot << buflr[5*i+4];
                        theplot << bufcomp[i];
                    }
                    
                    theplot.SendEndPlot();
                    
                }
                delete [] buf;
                delete [] bufcomp;
                
            }
            
            return 0L;
        }
    };
    
    plotHMatrix() : OneOperator(atype<long>(),atype<HMatrixVirt<K> **>()) {}
    
    E_F0 * code(const basicAC_F0 & args) const
    {
        return  new Op(args,t[0]->CastTo(args[0]));
    }
};

template<class K>
basicAC_F0::name_and_type  plotHMatrix<K>::Op::name_param[]= {
    {  "wait", &typeid(bool)},
    {  "dim", &typeid(long)}
};

template<class T, class U, class K, char trans>
class HMatrixInv {
public:
    const T t;
    const U u;
    
    struct HMatVirt: CGMatVirt<int,K> {
        const T tt;
        
        HMatVirt(T ttt) : tt(ttt), CGMatVirt<int,K>((*ttt)->nb_rows()) {}
        K*  addmatmul(K* x,K* Ax) const { (*tt)->mvprod_global(x, Ax); return Ax;}
    };
    
    struct HMatVirtPrec: CGMatVirt<int,K> {
        const T tt;
        std::vector<K> invdiag;
        
        HMatVirtPrec(T ttt) : tt(ttt), CGMatVirt<int,K>((*ttt)->nb_rows()), invdiag((*ttt)->nb_rows(),0) {
            std::vector<SubMatrix<K>*> diagblocks = (*tt)->get_MyStrictlyDiagNearFieldMats();
            std::vector<K> tmp((*ttt)->nb_rows(),0);
            for (int j=0;j<diagblocks.size();j++){
                SubMatrix<K>& submat = *(diagblocks[j]);
                int local_nr = submat.nb_rows();
                int local_nc = submat.nb_cols();
                int offset_i = submat.get_offset_i();
                int offset_j = submat.get_offset_j();
                for (int i=offset_i;i<offset_i+std::min(local_nr,local_nc);i++){
                    tmp[i] = 1./submat(i-offset_i,i-offset_i);
                }
            }
            (*tt)->cluster_to_target_permutation(tmp.data(),invdiag.data());
            MPI_Allreduce(MPI_IN_PLACE, &(invdiag[0]), (*ttt)->nb_rows(), wrapper_mpi<K>::mpi_type(), MPI_SUM, (*tt)->get_comm());
        }
        
        K*  addmatmul(K* x,K* Ax) const {
            for (int i=0; i<(*tt)->nb_rows(); i++)
                Ax[i] = invdiag[i] * x[i];
            return Ax;
        }
    };
    
    HMatrixInv(T v, U w) : t(v), u(w) {}
    
    void solve(U out) const {
        HMatVirt A(t);
        HMatVirtPrec P(t);
        double eps =1e-6;
        int niterx=2000;
        bool res=fgmres(A,P,1,(K*)*u,(K*)*out,eps,niterx,200,(mpirank==0)*verbosity);
    }
    
    static U inv(U Ax, HMatrixInv<T, U, K, trans> A) {
        A.solve(Ax);
        return Ax;
    }
    static U init(U Ax, HMatrixInv<T, U, K, trans> A) {
        Ax->init(A.u->n);
        return inv(Ax, A);
    }
};

template<class K>
class CompressMat : public OneOperator {
        public:
        class Op : public E_F0info {
                public:
                Expression a,b,c,d;

                static const int n_name_param = 8;
                static basicAC_F0::name_and_type name_param[] ;
                Expression nargs[n_name_param];
                long argl(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
                string* args(int i,Stack stack,string* a) const{ return nargs[i] ? GetAny<string*>( (*nargs[i])(stack) ): a;}
                double arg(int i,Stack stack,double a) const{ return nargs[i] ? GetAny<double>( (*nargs[i])(stack) ): a;}
                pcommworld argc(int i,Stack stack,pcommworld a ) const{ return nargs[i] ? GetAny<pcommworld>( (*nargs[i])(stack) ): a;}

                Op(const basicAC_F0 &  args,Expression aa,Expression bb,Expression cc, Expression dd) : a(aa),b(bb),c(cc),d(dd) {
                        args.SetNameParam(n_name_param,name_param,nargs);
                }
        };

        CompressMat() : OneOperator(atype<const typename CompressMat<K>::Op *>(),
                                    atype<KNM<K>*>(),
                                    atype<KN<double>*>(),
                                    atype<KN<double>*>(),
                                    atype<KN<double>*>()) {}

        E_F0 * code(const basicAC_F0 & args) const {
          return  new Op(args,t[0]->CastTo(args[0]),
                              t[1]->CastTo(args[1]),
                              t[2]->CastTo(args[2]),
                              t[3]->CastTo(args[3]));
        }
};

template<class K>
basicAC_F0::name_and_type  CompressMat<K>::Op::name_param[]= {
  {  "eps", &typeid(double)},
  {  "commworld", &typeid(pcommworld)},
  {  "eta", &typeid(double)},
  {  "minclustersize", &typeid(long)},
  {  "maxblocksize", &typeid(long)},
  {  "mintargetdepth", &typeid(long)},
  {  "minsourcedepth", &typeid(long)},
  {  "compressor", &typeid(string*)}
};

template<class K>
class MyMatrix: public IMatrix<K>{
        const KNM<K> &M;

public:
        MyMatrix(const KNM<K> &mat):IMatrix<K>(mat.N(),mat.M()),M(mat) {}

        K get_coef(const int& i, const int& j)const {return M(i,j);}
};

template<class K, int init>
AnyType SetCompressMat(Stack stack,Expression emat,Expression einter)
{ return SetCompressMat<K>(stack,emat,einter,init);}

template<class K>
AnyType SetCompressMat(Stack stack,Expression emat,Expression einter,int init)
{
  HMatrixVirt<K>** Hmat =GetAny<HMatrixVirt<K>** >((*emat)(stack));
  const typename CompressMat<K>::Op * mi(dynamic_cast<const typename CompressMat<K>::Op *>(einter));

  double epsilon=mi->arg(0,stack,1e-3);
  pcommworld pcomm=mi->argc(1,stack,nullptr);
  double eta=mi->arg(2,stack,10.);
  int minclustersize=mi->argl(3,stack,10);
  int maxblocksize=mi->argl(4,stack,1000000);
  int mintargetdepth=mi->argl(5,stack,htool::Parametres::mintargetdepth);
  int minsourcedepth=mi->argl(6,stack,htool::Parametres::minsourcedepth);
  string* pcompressor=mi->args(7,stack,0);

  SetMaxBlockSize(maxblocksize);
  SetMinClusterSize(minclustersize);
  SetEpsilon(epsilon);
  SetEta(eta);
  SetMinTargetDepth(mintargetdepth);
  SetMinSourceDepth(minsourcedepth);

  string compressor = pcompressor ? *pcompressor : "partialACA";

  MPI_Comm comm = pcomm ? *(MPI_Comm*)pcomm : MPI_COMM_WORLD;

  ffassert(einter);
  KNM<K> * pM = GetAny< KNM<K> * >((* mi->a)(stack));
  KNM<K> & M = *pM;

  KN<double> * px = GetAny< KN<double> * >((* mi->b)(stack));
  KN<double> & xx = *px;

  KN<double> * py = GetAny< KN<double> * >((* mi->c)(stack));
  KN<double> & yy = *py;

  KN<double> * pz = GetAny< KN<double> * >((* mi->d)(stack));
  KN<double> & zz = *pz;

  MyMatrix<K> A(M);

  ffassert(xx.n == M.N());
  ffassert(yy.n == M.N());
  ffassert(zz.n == M.N());

  std::vector<htool::R3> p(xx.n);
  for (int i=0; i<xx.n; i++)
    p[i] = {xx[i], yy[i], zz[i]};

  //cout << M.N() << " " << xx.M() << " " << yy.n << " " << zz.n<< endl;
  if (init) delete *Hmat;

  if ( compressor == "" || compressor == "partialACA")
       *Hmat = new HMatrixImpl<K,partialACA>(A,p);
   else if (compressor == "fullACA")
       *Hmat = new HMatrixImpl<K,fullACA>(A,p);
   else if (compressor == "SVD")
       *Hmat = new HMatrixImpl<K,SVD>(A,p);
   else {
       cerr << "Error: unknown htool compressor \""+compressor+"\"" << endl;
       ffassert(0);
   }

  return Hmat;
}

template<class K>
void addHmat() {
    Dcl_Type<HMatrixVirt<K>**>(Initialize<HMatrixVirt<K>*>, Delete<HMatrixVirt<K>*>);
    Dcl_TypeandPtr<HMatrixVirt<K>*>(0,0,::InitializePtr<HMatrixVirt<K>*>,::DeletePtr<HMatrixVirt<K>*>);
    //atype<HMatrix<LR ,K>**>()->Add("(","",new OneOperator2_<string*, HMatrix<LR ,K>**, string*>(get_infos<LR,K>));
    
    Add<HMatrixVirt<K>**>("infos",".",new OneOperator1_<std::map<std::string, std::string>*, HMatrixVirt<K>**>(get_infos));
    
    Dcl_Type<Prod<KN<K>*, K>>();
    TheOperators->Add("*", new OneOperator2<Prod<KN<K>*, K>, HMatrixVirt<K>**, KN<K>*>(Build));
    TheOperators->Add("=", new OneOperator2<KN<K>*, KN<K>*, Prod<KN<K>*, K>>(Prod<KN<K>*, K>::mv));
    TheOperators->Add("<-", new OneOperator2<KN<K>*, KN<K>*, Prod<KN<K>*, K>>(Prod<KN<K>*, K>::init));
    
    addInv<HMatrixVirt<K>*, HMatrixInv, KN<K>, K>();
    
    Global.Add("display","(",new plotHMatrix<K>);
    
    // to dense:
    TheOperators->Add("=",
                      new OneOperator2_<KNM<K>*, KNM<K>*, HMatrixVirt<K>**,E_F_StackF0F0>(To<KNM<K>, K, 1>));
    TheOperators->Add("<-",
                      new OneOperator2_<KNM<K>*, KNM<K>*, HMatrixVirt<K>**,E_F_StackF0F0>(To<KNM<K>, K, 0>));
#ifdef WITH_petsccomplex
    if(std::is_same<K, PetscScalar>::value) {
        TheOperators->Add("=",
                new OneOperator2_<Dmat*, Dmat*, HMatrixVirt<K>**,E_F_StackF0F0>(To<Dmat, K, 1>));
        TheOperators->Add("<-",
                new OneOperator2_<Dmat*, Dmat*, HMatrixVirt<K>**,E_F_StackF0F0>(To<Dmat, K, 0>));
    }
#endif
    Dcl_Type<const typename CompressMat<K>::Op *>();
    //Add<const typename assembleHMatrix<LR, K>::Op *>("<-","(", new assembleHMatrix<LR, K>);

    TheOperators->Add("=",
    new OneOperator2_<HMatrixVirt<K>**,HMatrixVirt<K>**,const typename CompressMat<K>::Op*,E_F_StackF0F0>(SetCompressMat<K, 1>));
    TheOperators->Add("<-",
    new OneOperator2_<HMatrixVirt<K>**,HMatrixVirt<K>**,const typename CompressMat<K>::Op*,E_F_StackF0F0>(SetCompressMat<K, 0>));

    Global.Add("compress","(",new CompressMat<K>);
}




/// Domain integration - operator



class CPartBemDI: public E_F0mps {
public:
    
    static const int n_name_param =3; //12;
    static basicAC_F0::name_and_type name_param[] ;
    Expression nargs [n_name_param];
    enum typeofkind  { int1dx1d=0, int2dx2d=1, int2dx1d=2, int1dx2d=3  } ;
    typeofkind  kind; //  0
    int d, dHat; // 3d
    typedef const CPartBemDI* Result;
    Expression Th;
    vector<Expression> what;
    vector<int> whatis; // 0 -> long , 1 -> array ???
    
    
    CPartBemDI(const basicAC_F0 & args,typeofkind b=int1dx1d,int ddim=3, int ddimHat=1) // always ddim=3d
    :kind(b),d(ddim),dHat(ddimHat),
    Th(0), what(args.size()-1), whatis(args.size()-1)
    
    {
        args.SetNameParam(n_name_param,name_param,nargs);
        
        if(d==3 && dHat==1)
            Th=CastTo<pmeshL>(args[0]);
        else if(d==3 && dHat==2)
            Th=CastTo<pmeshS>(args[0]);
        else ffassert(0); // a faire
        
        int n=args.size();
        for (int i=1;i<n;i++)
            if(!BCastTo<KN_<long> >(args[i]) ) {
                whatis[i-1]=0;
                what[i-1]=CastTo<long>(args[i]);
            }
            else {
                whatis[i-1]=1;
                what[i-1]=CastTo<KN_<long> >(args[i]);
            }
        
    }
    static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmesh>(), true);} // all type
    AnyType operator()(Stack ) const  { return SetAny<const CPartBemDI *>(this);}
    
    operator aType () const { return atype<const CPartBemDI *>();}
    
    static  E_F0 * f(const basicAC_F0 & args) { return new CPartBemDI(args);}
    
    const Fem2D::QuadratureFormular & FIT(Stack) const ;
    const Fem2D::QuadratureFormular1d & FIE(Stack) const ;
    const Fem2D::GQuadratureFormular<Fem2D::R3> & FIV(Stack) const ;  // 3d
    
};



class CBemDomainOfIntegration: public E_F0mps {
public:
    
    static const int n_name_param =3; //12;
    static basicAC_F0::name_and_type name_param[] ;
    Expression nargs_t [n_name_param];
    // typeofkind  kind; //  0
    int d_s, dHat_s, d_t, dHat_t; // 3d
    typedef const CBemDomainOfIntegration* Result;
    Expression Th_s, Th_t, BemPartDI;
    vector<Expression> what_s, what_t;
    vector<int> whatis_s, whatis_t;
    //CBemDomainOfIntegration();
    CBemDomainOfIntegration(const basicAC_F0 & args_t, int ddim=3, int ddimHat=1) // always ddim=3d
    :d_s(0),dHat_s(0),d_t(ddim),dHat_t(ddimHat),
    Th_s(0), what_s(0), whatis_s(0),
    Th_t(0), what_t(args_t.size()-1), whatis_t(args_t.size()-1)
    
    {
        args_t.SetNameParam(n_name_param,name_param,nargs_t);
        // acces to the value of the fist di Th_s
        const CPartBemDI * sourceDI(dynamic_cast<const CPartBemDI*>((Expression) args_t[0]));
        
        CPartBemDI::typeofkind kind_s= sourceDI->kind; // int1dx1d=0, int2dx2d=1, int2dx1d=2, int1dx2d=3
        d_s=sourceDI->d;
        dHat_s=sourceDI->dHat;
        
        // check the integral operator integral
        if(kind_s==0 || kind_s==2)
            Th_t=CastTo<pmeshL>(args_t[1]);
        else if(kind_s==1 || kind_s==3)
            Th_t=CastTo<pmeshS>(args_t[1]);
        else if(kind_s==0 || kind_s==3)
            Th_s=sourceDI->Th; //Th_s=CastTo<pmeshL>(sourceDI->Th);
        else if(kind_s==1 || kind_s==2)
            Th_s=sourceDI->Th;
        else ffassert(0); // a faire
        
        if (mpirank==0 && verbosity >5)
            cout << " CBemDomainOfIntegration " << kind_s << " Th_s: " << &Th_s << " d_s= " << d_s << " dHat_s= " << dHat_s << " " <<
            " Th_t: " << &Th_t << " d_t= " << d_t << " dHat_t= " << dHat_t << " " << endl;
        
        // read the argument (Th_t,.....)
        int n_t=args_t.size();
        for (int i=2;i<n_t;i++)
            if(!BCastTo<KN_<long> >(args_t[i]) ) {
                whatis_t[i-1]=0;
                what_t[i-1]=CastTo<long>(args_t[i]);
            }
            else {
                whatis_t[i-1]=1;
                what_t[i-1]=CastTo<KN_<long> >(args_t[i]);
            }
        int n_s=(sourceDI->what).size();
        for (int i=1;i<n_s;i++) {
            
            whatis_s[i-1]=sourceDI->whatis[i-1];
            what_s[i-1]=sourceDI->what[i-1];
        }
        
    }
    static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<const CPartBemDI *>(), true);} // all type
    AnyType operator()(Stack ) const  { return SetAny<const CBemDomainOfIntegration *>(this);}
    
    operator aType () const { return atype<const CBemDomainOfIntegration *>();}
    
    //static  E_F0 * f(const basicAC_F0 & args_s,const basicAC_F0 & args_t) { return new CBemDomainOfIntegration(args_s,args_t);}
    static  E_F0 * f(const basicAC_F0 & args_s) { return new CBemDomainOfIntegration(args_s);}
    
    const Fem2D::QuadratureFormular & FIT(Stack) const ;
    const Fem2D::QuadratureFormular1d & FIE(Stack) const ;
    const Fem2D::GQuadratureFormular<Fem2D::R3> & FIV(Stack) const ;  // 3d
    
};



basicAC_F0::name_and_type  CPartBemDI::name_param[]= {
    { "qft", &typeid(const Fem2D::QuadratureFormular *)},    //0
    { "qfe", &typeid(const Fem2D::QuadratureFormular1d *)},
    { "qforder",&typeid(long)},     // 2
    //{ "qfnbpT",&typeid(long)},
    //{ "qfnbpE",&typeid(long)},
    //{ "optimize",&typeid(long)},
    //{ "binside",&typeid(double)},
    //{ "mortar",&typeid(bool)},
    { "qfV", &typeid(const Fem2D::GQuadratureFormular<Fem2D::R3> *)},    // 8 ->  3
    //{ "levelset",&typeid(double)},
    //{ "mapt",&typeid(E_Array)},
    //{ "mapu",&typeid(E_Array)}
    
    
};
basicAC_F0::name_and_type  CBemDomainOfIntegration::name_param[]= {
    { "qft", &typeid(const Fem2D::QuadratureFormular *)},    //0
    { "qfe", &typeid(const Fem2D::QuadratureFormular1d *)},
    { "qforder",&typeid(long)},     // 2
    //{ "qfnbpT",&typeid(long)},
    //{ "qfnbpE",&typeid(long)},
    //{ "optimize",&typeid(long)},
    //{ "binside",&typeid(double)},
    //{ "mortar",&typeid(bool)},
    { "qfV", &typeid(const Fem2D::GQuadratureFormular<Fem2D::R3> *)},    // 8 ->  3
    //{ "levelset",&typeid(double)},
    //{ "mapt",&typeid(E_Array)},
    //{ "mapu",&typeid(E_Array)}
    
    
};



class CPartBemDI1d1d: public CPartBemDI {
public:
    CPartBemDI1d1d( const basicAC_F0 & args_s) :CPartBemDI(args_s,int1dx1d,3,1) {}
    static  E_F0 * f(const basicAC_F0 & args_s) { return new CPartBemDI(args_s,int1dx1d,3,1);}
    static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmeshL>(), true);} // all type
};

class CPartBemDI2d2d: public CPartBemDI {
public:
    CPartBemDI2d2d(const basicAC_F0 & args_s) :CPartBemDI(args_s,int2dx2d,3,2) {}
    static  E_F0 * f(const basicAC_F0 & args_s) { return new CPartBemDI(args_s,int2dx2d,3,2);}
    static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmeshS>(), true);} // all type
};

class CPartBemDI1d2d: public CPartBemDI {
public:
    CPartBemDI1d2d(const basicAC_F0 & args_s) :CPartBemDI(args_s,int1dx2d,3,1) {}
    static  E_F0 * f(const basicAC_F0 & args_s) { return new CPartBemDI(args_s,int2dx2d,3,1);}
    static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmeshS>(), true);} // all type
};
class CPartBemDI2d1d: public CPartBemDI {
public:
    CPartBemDI2d1d(const basicAC_F0 & args_s) :CPartBemDI(args_s,int2dx1d,3,2) {}
    static  E_F0 * f(const basicAC_F0 & args_s) { return new CPartBemDI(args_s,int2dx2d,3,2);}
    static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmeshS>(), true);} // all type
};

//// begin type BEM kernel / potential

class BemKernel : public RefCounter {
public:
    int typeKernel[2]={0,0};  // Laplace, Helmholtz
    // typeKernel ={SL, DL, HS, TDL} and determine equation Laplace, Helmholtz if k==0 or not
    std::complex<double> wavenum[2]={0,0}; // parameter to Helmholtz
    std::complex<double> coeffcombi[2]={0,0};
    BemKernel(){}
    BemKernel(const string tkernel, Complex alpha=1 , Complex k=0) {
        
        coeffcombi[0]=alpha;
        wavenum[0]=k;
        
        if(!tkernel.compare("SL"))
            typeKernel[0] = 1;
        else if(!tkernel.compare("DL"))
            typeKernel[0] = 2;
        else if(!tkernel.compare("HS"))
            typeKernel[0] = 3;
        else if(!tkernel.compare("TDL"))
            typeKernel[0] = 4;
        else if(!tkernel.compare("CST"))
            typeKernel[0] = 5;
        else
            ExecError("unknow BEM kernel type ");
        
        if(mpirank==0 && verbosity>5)
            cout << "type BEM kernel " << tkernel <<": " << typeKernel[0] << " coeff combi " << coeffcombi[0] << " wave number "<< wavenum[0] << endl;
    }
    ~BemKernel() {}
    
    BemKernel(const BemKernel &Bk) {
        for(int i=0;i<2;i++) {typeKernel[i]=Bk.typeKernel[i]; wavenum[i]=Bk.wavenum[i]; coeffcombi[i]=Bk.coeffcombi[i]; } } ;
    
    
private:
    //BemKernel(const BemKernel &);
    void operator=(const BemKernel &);
};


class listBemKernel {
public:
    list<BemKernel const  *> *lbk;
    void init()  { lbk=new list<BemKernel const  *>;}
    void destroy() { delete lbk;}
    listBemKernel(Stack s,BemKernel const*bk) : lbk(Add2StackOfPtr2Free(s,new list<BemKernel const*>)) { lbk->push_back(bk);}
    listBemKernel(Stack s,BemKernel const*const bka,BemKernel const* const bkb) : lbk(Add2StackOfPtr2Free(s,new list<BemKernel const*>)) { lbk->push_back(bka);lbk->push_back(bkb);}
    listBemKernel(Stack s,const listBemKernel &l,BemKernel const*const bk) : lbk(Add2StackOfPtr2Free(s,new list<BemKernel const*>(*l.lbk))) { lbk->push_back(bk);}
    listBemKernel(){};
};


template<class RR,class AA=RR,class BB=AA>
struct Op_addBemKernel: public binary_function<AA,BB,RR> {
    static RR f(Stack s,const AA & a,const BB & b) {
        if (mpirank==0 && verbosity>10) cout << "test " <<typeid(RR).name() << " " << typeid(AA).name() << " " << typeid(BB).name() <<endl;
        return RR(s,a,b);}
};


template<bool INIT,class RR,class AA=RR,class BB=AA>
struct Op_setBemKernel: public binary_function<AA,BB,RR> {
    static RR f(Stack stack, const AA & a,const BB & b)
    {
        ffassert(a);
        const pBemKernel p=combKernel(b);
        
        if (!INIT && *a)
            (**a).destroy( );
        *a = p;
        return a;
    }
};


class BemPotential : public RefCounter {
public:
    
    int typePotential;  // Laplace, Helmholtz
    // typePotential ={SL=0, DL=1, HS=2, TDL=3} and determine equation Laplace, Helmholtz if k==0 or not
    std::complex<double> wavenum; // parameter to Helmholtz
    BemPotential(){}
    BemPotential(const string tpotential, Complex k=0) : typePotential(-1), wavenum(k) {
        
        if(!tpotential.compare("SL"))
            typePotential = 1;
        else if(!tpotential.compare("DL"))
            typePotential = 2;
        else if(!tpotential.compare("CST"))
            typePotential = 3;
        else
            ExecError("unknow BEM Potential type ");
        
        if(mpirank==0 && verbosity>5)
            cout << "type BEM Potential " << tpotential <<": " << typePotential << " wave number "<< wavenum << endl;
    }
    
    ~BemPotential() {}
    
private:
    BemPotential(const BemPotential &);
    void operator=(const BemPotential &);
};

// fusion of Bem Kernel in case combined kernels
BemKernel *combKernel (listBemKernel const &lbemker){
    int kk=0;
    bool LaplaceK=false, HelmholtzK=false;
    const list< const BemKernel * > lbk(*lbemker.lbk);
    BemKernel *combBemKernel=new BemKernel();
    for (list< const BemKernel * >::const_iterator i = lbk.begin( ); i != lbk.end( ); i++) {
        if (!*i) continue;
        const BemKernel &bkb(**i);
        combBemKernel->typeKernel[kk] = bkb.typeKernel[0];
        // test same equation kernel
        (combBemKernel->wavenum[kk]!=0.) ? HelmholtzK=true : LaplaceK=true;
        
        //  ExecError(" combined kernel have to be the same type equation Laplace or Helmholtz");
        combBemKernel->coeffcombi[kk] = bkb.coeffcombi[0];
        combBemKernel->wavenum[kk] = bkb.wavenum[0];
        if(kk>4) ExecError(" combined kernel: 4 max kernels  ");
        kk++;
    }
    // check the same wave number ?
    if( HelmholtzK== LaplaceK) ExecError(" combined kernel: must be same equation kernels Laplace or Helmholtz");
    if (mpirank==0 && verbosity>5)
        for (int i=0;i<kk;i++) cout << "combined type BEM kernel " << combBemKernel->typeKernel[i] << " coeff combi " <<
            combBemKernel->coeffcombi[i] << " wave number "<< combBemKernel->wavenum[i] << endl;
    
    return combBemKernel;
}

inline pBemKernel *  initKernel_Helmholtz(pBemKernel * const & p, string * const & s,std::complex<double> const & alpha,std::complex<double> const & k) {
    BemKernel * m;
    if(mpirank==0 && verbosity > 5)
        cout << " initBemKernel " << *s << endl;
    *p= m =new BemKernel(*s,alpha,k);
    return p;
}

inline pBemKernel *  initKernel_Laplace(pBemKernel * const & p, string * const & s,std::complex<double> const & alpha) {
    return initKernel_Helmholtz(p,s,alpha,0);
    
}

inline pBemKernel *  initKernel_default(pBemKernel * const & p, string * const & s) {
    return initKernel_Helmholtz(p,s,1,0);
    
}

inline pBemPotential *  initPotential_Helmholtz(pBemPotential * const & p, string * const & s, std::complex<double> const & k) {
    BemPotential * m;
    if(mpirank==0 && verbosity > 5)
        cout << " initPotential " << *s << endl;
    *p= m =new BemPotential(*s,k);
    return p;
}

inline pBemPotential *  initPotential_default(pBemPotential * const & p, string * const & s) {
    return initPotential_Helmholtz(p,s,0);
    
}


// BEM variational form

// define a bilinear form for BEM
class FoperatorKBEM : public E_F0mps { public:
    typedef const FoperatorKBEM* Result;
    static const int n_name_param = 1;
    typedef finconnue * Fi;
    typedef ftest * Ft;
    typedef fkernel * KBem;
    
    Fi fi;
    Ft ft;
    Expression kbem;
    
    FoperatorKBEM(const basicAC_F0 & args) :fi(0),ft(0),kbem(0) {
        ffassert(args.size()==3);
 
        kbem= CastTo<KBem>(args[0]);
        fi= dynamic_cast<Fi>(CastTo<Fi>(args[1]));
        ft= dynamic_cast<Ft>(CastTo<Ft>(args[2]));
        ffassert(kbem && fi && ft);
    }
    
    AnyType operator()(Stack ) const { return SetAny<Result>(this);}
    operator aType () const { return atype<Result>();}
    FoperatorKBEM(const FoperatorKBEM & fk) : kbem(fk.kbem),fi(fk.fi),ft(fk.ft) {}
    
};

// define a bilinear form for BEM



class BemFormBilinear : virtual public E_F0mps { public:
    int type=-1;
    static  E_F0 * f(const basicAC_F0 & args);
};


class BemKFormBilinear : public BemFormBilinear { public:
    typedef const BemFormBilinear* Result;
    typedef const CBemDomainOfIntegration * A;
    typedef const FoperatorKBEM * B;
    A  di;
    FoperatorKBEM * b;
    
    BemKFormBilinear(const basicAC_F0 & args) {
        di= dynamic_cast<A>(CastTo<A>(args[0]));
        B Kb= dynamic_cast<B>(CastTo<B>(args[1]));
        b= new FoperatorKBEM(*Kb);
        ffassert(di && Kb);
        type=1;
        
        
    };
    
    static  E_F0 * f(const basicAC_F0 & args) { return new BemKFormBilinear(args);}
    static ArrayOfaType  typeargs() { return ArrayOfaType(atype<A>(),atype<B>());}
    AnyType operator()(Stack ) const { return SetAny<Result>(this);}
    operator aType () const { return atype<Result>();}
    
    
    BemKFormBilinear(A a,Expression bb) : di(a),b(new FoperatorKBEM(*dynamic_cast<FoperatorKBEM *>(bb)))
    {ffassert(b);}
    BemKFormBilinear operator-() const { return  BemKFormBilinear(di,C_F0(TheOperators,"-",C_F0(b,atype<FoperatorKBEM>())));}
    
    BemKFormBilinear(const BemKFormBilinear & fb) : di(fb.di),b(new FoperatorKBEM(*fb.b) ) {}
    
};


class TypeFormBEM: public ForEachType<const BemFormBilinear*> {
public:
    TypeFormBEM() : ForEachType<const BemFormBilinear*>(0,0) {}
    void SetArgs(const ListOfId *lid) const {
        SetArgsFormLinear(lid,2);    }
    
    Type_Expr SetParam(const C_F0 & c,const ListOfId *l,size_t & top) const
    { return Type_Expr(this,CastTo(c));}
    
    
    C_F0 Initialization(const Type_Expr & e) const
    {
        return C_F0(); }
    Type_Expr construct(const Type_Expr & e) const
    {
        return e; }
    
};

// define the function BEM(k,u,v)
class FormalKBEMcode : public OneOperator{
public:
    
    FormalKBEMcode( ): OneOperator(atype<C_F0>(),atype<pBemKernel>(), atype<finconnue*>(), atype<ftest*>()) {}
    FormalKBEMcode(int  ): OneOperator(atype<C_F0>(),atype<pBemKernel>()) {}
    E_F0 *  code(const basicAC_F0 & ) const {ffassert(0);}
    C_F0  code2(const basicAC_F0 &args) const {
        Expression e=new FoperatorKBEM(args);
        aType r=atype<const FoperatorKBEM *>();
        return C_F0(e,r) ;}
    
    AnyType operator()(Stack s)  const {ffassert(0);return 0L;}
    
};



// define the function POT(k,u,v)
class FoperatorPBEM : public E_F0mps { public:
    typedef const FoperatorPBEM* Result;
    typedef finconnue * Fi;
    typedef ftest * Ft;
    typedef fpotential * Pot;
    
    Fi fi;
    Ft ft;
    Expression pot;
    
    FoperatorPBEM(const basicAC_F0 & args) :fi(0),ft(0),pot(0) {
        ffassert(args.size()==3);
        
        pot= CastTo<Pot>(args[0]);
        fi= dynamic_cast<Fi>(CastTo<Fi>(args[1]));
        ft= dynamic_cast<Ft>(CastTo<Ft>(args[2]));
        
        ffassert(pot && fi && ft);
    }
    
    AnyType operator()(Stack ) const { return SetAny<Result>(this);}
    operator aType () const { return atype<Result>();}
    FoperatorPBEM(const FoperatorPBEM & fk) : pot(fk.pot),fi(fk.fi),ft(fk.ft){}
    
};


// define a bilinear form for BEM
class BemPFormBilinear : public BemFormBilinear { public:
    typedef const BemFormBilinear* Result;
    typedef const CDomainOfIntegration * A;
    typedef const FoperatorPBEM * B;
    A  di;
    FoperatorPBEM * b;
    BemPFormBilinear(const basicAC_F0 & args) {
        di= dynamic_cast<A>(CastTo<A>(args[0]));
        B Kb= dynamic_cast<B>(CastTo<B>(args[1]));
        b= new FoperatorPBEM(*Kb);
        ffassert(di && Kb);
        type=2;
    }
    
    static  E_F0 * f(const basicAC_F0 & args) { return new BemPFormBilinear(args);}
    static ArrayOfaType  typeargs() { return ArrayOfaType(atype<A>(),atype<B>());}// all type
    AnyType operator()(Stack ) const { return SetAny<Result>(this);}
    operator aType () const { return atype<Result>();}
    
    
    BemPFormBilinear(A a,Expression bb) : di(a),b(new FoperatorPBEM(*dynamic_cast<FoperatorPBEM *>(bb))/*->Optimize(currentblock) FH1004 */)
    {ffassert(b);}
    BemPFormBilinear operator-() const { return  BemPFormBilinear(di,C_F0(TheOperators,"-",C_F0(b,atype<FoperatorPBEM>())));}
    
    BemPFormBilinear(const BemPFormBilinear & fb) : di(fb.di),b(new FoperatorPBEM(*fb.b) ) {}
    
};

class FormalPBEMcode : public OneOperator{
public:
    
    FormalPBEMcode( ): OneOperator(atype<C_F0>(),atype<pBemPotential>(), atype<finconnue*>(), atype<ftest*>()) {}
    FormalPBEMcode(int  ): OneOperator(atype<C_F0>(),atype<pBemPotential>()) {}
    E_F0 *  code(const basicAC_F0 & ) const {ffassert(0);}
    C_F0  code2(const basicAC_F0 &args) const {
        Expression e=new FoperatorPBEM(args);
        aType r=atype<const FoperatorPBEM *>();
        return C_F0(e,r) ;}
    
    AnyType operator()(Stack s)  const {ffassert(0);return 0L;}
    
};


//// end type BEM kernel / potential



bool C_args::IsBemBilinearOperator() const {
for (const_iterator i=largs.begin(); i != largs.end();i++) {
    C_F0  c= *i;
    aType r=c.left();
    if  ( r!= atype<const class BemFormBilinear *>() )
         return false;
}
return true;}

// bem + fem
bool C_args::IsMixedBilinearOperator() const {
if ( this->IsBilinearOperator() && this->IsBemBilinearOperator()  ) return true;
    
return false;}


template<class R, class v_fes1, class v_fes2>
struct OpHMatrixtoBEMForm
: public OneOperator
{
    typedef typename Call_FormBilinear<v_fes1,v_fes2>::const_iterator const_iterator;
    int init;
    class Op : public E_F0mps {
    public:
        Call_FormBilinear<v_fes1,v_fes2> *b;
        Expression a;
        int init;
        AnyType operator()(Stack s)  const ;
        
        Op(Expression aa,Expression  bb,int initt)
        : b(new Call_FormBilinear<v_fes1,v_fes2>(* dynamic_cast<const Call_FormBilinear<v_fes1,v_fes2> *>(bb))),a(aa),init(initt)
        { assert(b && b->nargs);
        }
        operator aType () const { return atype<HMatrixVirt<R> **>();}
        
    };
    
    E_F0 * code(const basicAC_F0 & args) const
    {
        Expression p=args[1];
        Call_FormBilinear<v_fes1,v_fes2> *t( new Call_FormBilinear<v_fes1,v_fes2>(* dynamic_cast<const Call_FormBilinear<v_fes1,v_fes2> *>(p))) ;
        return  new Op(to<HMatrixVirt<R> **>(args[0]),args[1],init);}
  
     OpHMatrixtoBEMForm(int initt=0) :
     OneOperator(atype<HMatrixVirt<R> **>(),atype<HMatrixVirt<R> **>(),atype<const Call_FormBilinear<v_fes1,v_fes2>*>()),
        init(initt)
        {}
    
};



pair<BemKernel*,double> getBemKernel(Stack stack, const list<C_F0> & largs)  {
    list<C_F0>::const_iterator ii,ib=largs.begin(),ie=largs.end();
    
    BemKernel* K;
    double alpha=0.;
    
    bool haveBemBilinearOperator=false, haveBilinearOperator=false;
    
    for (ii=ib;ii != ie;ii++) {
        Expression e=ii->LeftValue();
        aType r = ii->left();
        
        if (r==atype<const  BemFormBilinear *>() && !haveBemBilinearOperator) {
            BemKFormBilinear * bb=new BemKFormBilinear(*dynamic_cast<const BemKFormBilinear *>(e));
            FoperatorKBEM * b=const_cast<  FoperatorKBEM *>(bb->b);
            if (b == NULL) {
                if(mpirank == 0) cout << "dynamic_cast error" << endl; }
            else
                K=GetAny<BemKernel*>((*b->kbem)(stack));
            haveBemBilinearOperator=true;
        }
        
        else if (r==atype<const  FormBilinear *>() && !haveBilinearOperator) {
            
            const  FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
            const CDomainOfIntegration & di= *bb->di;
            // check the integration (keyword)
            ffassert( (di.kind == CDomainOfIntegration::int1d && di.isMeshL) || (di.kind == CDomainOfIntegration::int2d && di.isMeshS) );  //check only necessary in surface case
            
            BilinearOperator * Op=const_cast<  BilinearOperator *>(bb->b);
            if (Op == NULL) {
                if(mpirank == 0) cout << "dynamic_cast error" << endl; }
            
            
            
            BilinearOperator::const_iterator l=Op->v.begin();
            
            BilinearOperator::K ll(*l);   //  LinearComb<pair<MGauche,MDroit>,C_F0> BilinearOperator;
            pair<int,int> finc(ll.first.first),ftest(ll.first.second);
                                   
            alpha = GetAny<double>(ll.second.eval(stack));
            if(mpirank == 0 && verbosity>5) cout << " test coeff mass matrix " << alpha << endl;
          
            if(mpirank == 0 && verbosity>5) {
                cout << "FormBilinear: number of unknow finc=" << finc.first << " ,ftest= " << ftest.first << endl;
                cout << "FormBilinear: operator order finc=" << finc.second << " ,ftest= " << ftest.second << endl;      // ordre   only op_id=0
            }
            ffassert(finc.first==0 && ftest.first==0);
            ffassert(finc.second==0 && ftest.second==0);
            haveBilinearOperator=true; // or ffassert(Op->v.size()==1); // check size 1
        }
        else
            ffassert(0);
    }
    return std::make_pair(K, alpha); //K;
}


BemPotential* getBemPotential(Stack stack, const list<C_F0> & largs)  {
    list<C_F0>::const_iterator ii,ib=largs.begin(),ie=largs.end();
    
    BemPotential* P;
    
    for (ii=ib;ii != ie;ii++) {
        Expression e=ii->LeftValue();
        aType r = ii->left();
        ffassert (r==atype<const  BemFormBilinear *>());
        
        BemPFormBilinear * bb=new BemPFormBilinear(*dynamic_cast<const BemPFormBilinear *>(e));
        FoperatorPBEM * b=const_cast<  FoperatorPBEM *>(bb->b);
        if (b == NULL) {
            if(mpirank == 0) cout << "dynamic_cast error" << endl; }
        else
            P=GetAny<BemPotential*>((*b->pot)(stack));
    }
    return P;
}

struct Data_Bem_Solver
: public Data_Sparse_Solver {
    double eta;
       int minclustersize,maxblocksize,mintargetdepth,minsourcedepth;
       string compressor;
       
       Data_Bem_Solver()
       : Data_Sparse_Solver(),
       eta(10.),
       minclustersize(10),
       maxblocksize(1000000),
       mintargetdepth(htool::Parametres::mintargetdepth),
       minsourcedepth(htool::Parametres::minsourcedepth),
       compressor("partialACA")
    
      {epsilon=1e-3;}
     
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


bool iscombinedKernel(BemKernel *K ) {
    bool iscombined=false;
    for (int i=0;i<2;i++)
        iscombined= K->coeffcombi[i]!=0. ? 1 : 0;
    return iscombined;
}


const EquationEnum whatEquationEnum(BemKernel *K,int i) {
    bool wavenum=false;
    int nk=2;
    for (int i=0;i<nk;i++)
        wavenum= K->wavenum[i]!=0. ? 1 : 0;
    
    const EquationEnum equEnum = wavenum ? HE : LA ;
    return equEnum;
    
}

BIOpKernelEnum whatTypeEnum(BemKernel *K,int i) {
    BIOpKernelEnum pKernel;
    switch(K->typeKernel[i]) {
        case 1: pKernel=SL_OP ; break;
        case 2: pKernel=DL_OP ; break;
        case 3: pKernel=HS_OP ; break;
        case 4: pKernel=TDL_OP ; break;
    }
    const BIOpKernelEnum cpKernel=pKernel;
    return cpKernel;
}


PotKernelEnum whatTypeEnum(BemPotential *P) {
    PotKernelEnum pPotential;
    switch(P->typePotential) {
        case 1: pPotential=SL_POT ; break;
        case 2: pPotential=DL_POT ; break;
        case 3: pPotential=CST_POT ; break;
            //case 4: pPotential=TDL_POT ; break;
    }
    const PotKernelEnum cpPotential=pPotential;
    return cpPotential;
}


int typeVFBEM(const list<C_F0> & largs, Stack stack)
{
    list<C_F0>::const_iterator ii,ib=largs.begin(),ie=largs.end();
    
    int VVFBEM =-1, ik=-1;;
    for (ii=ib;ii != ie;ii++) {
        Expression e=ii->LeftValue();
        aType r = ii->left();
        
        if (r==atype<const  BemFormBilinear *>()) {
            
            BemFormBilinear * bb= GetAny<BemFormBilinear *>((*e)(0));
            VVFBEM = bb->type;
        }
        ffassert(ik);
    }
    return VVFBEM;
}


template <class R, typename P, class MMesh>
void ff_BIO_Generator(HMatrixVirt<R>** Hmat, BemKernel *typeKernel, Dof<P>& dof, double alpha,string &compressor,vector<htool::R3> &p1,MPI_Comm &comm) {
    IMatrix<R>* generator;
    
    BIOpKernelEnum ker1 = whatTypeEnum(typeKernel,0), ker2 = whatTypeEnum(typeKernel,1);;
    double kappa1 = typeKernel->wavenum[0].real(), kappa2 = typeKernel->wavenum[1].real();
    bool iscombined = iscombinedKernel(typeKernel);
    if(iscombined) ffassert(kappa1==kappa2);
    std::complex<double> coeff1=typeKernel->coeffcombi[0], coeff2=typeKernel->coeffcombi[1];
    
    
    // Eq=HE LA
    // Ker= SL_OP DL_OP HS_OP TDL_OP
    if (kappa1 && !iscombined && !kappa2 && !alpha) {
        switch (ker1) {
            case SL_OP : generator=new BIO_Generator<BIOpKernel<HE,SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1);
                if(mpirank == 0 && verbosity>5) cout << " call bemtool func BIOpKernel<HE,SL_OP ..." << endl; break;
            case DL_OP : generator=new BIO_Generator<BIOpKernel<HE,DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<HE,DL_OP ..." << endl; break;
            case HS_OP : generator=new BIO_Generator<BIOpKernel<HE,HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<HE,HS_OP ..." << endl; break;
            case TDL_OP : generator=new BIO_Generator<BIOpKernel<HE,TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<HE,TDL_OP ..." << endl; break;
            case CST_OP :  if(verbosity>5) cout << " no tested BIOpKernel<HE,CST_OP" << endl; ffassert(0); break;
        }
    }
    else if (!kappa1 && !iscombined && !kappa2 && !alpha ){
        switch (ker1) {
            case SL_OP : generator=new BIO_Generator<BIOpKernel<LA,SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<LA,SL_OP ..." << endl; break;
            case DL_OP : generator=new BIO_Generator<BIOpKernel<LA,DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<LA,TDL_OP ..." << endl; break;
            case HS_OP : generator=new BIO_Generator<BIOpKernel<LA,HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<LA,HS_OP ..." << endl; break;
            case TDL_OP : generator=new BIO_Generator<BIOpKernel<LA,TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<LA,TDL_OP ..." << endl; break;
            case CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested BIOpKernel<LA,CST_OP" << endl; ffassert(0); break;
        }
    }
    
    
    //BIO_Generator_w_mass
    
    else if (kappa1 && !iscombined && !kappa2 && alpha) {
        switch (ker1) {
            case SL_OP : generator=new BIO_Generator_w_mass<BIOpKernel<HE,SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,alpha);
                if(mpirank == 0 && verbosity>5) cout << " call bemtool func BIO_Generator_w_mass<HE,SL_OP ...alpha coeff mass matrix=" << alpha << endl; break;
            case DL_OP : generator=new BIO_Generator_w_mass<BIOpKernel<HE,DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,alpha);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<HE,DL_OP ...alpha coeff mass matrix=" << alpha << endl; break;
            case HS_OP : generator=new BIO_Generator_w_mass<BIOpKernel<HE,HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,alpha);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<HE,HS_OP ...alpha coeff mass matrix=" << alpha << endl; break;
            case TDL_OP : generator=new BIO_Generator_w_mass<BIOpKernel<HE,TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,alpha);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<HE,TDL_OP ...alpha coeff mass matrix=" << alpha << endl; break;
            case CST_OP :  if(mpirank == 0 && verbosity>5) cout << " no tested BIO_Generator_w_mass<HE,CST_OP" << endl; ffassert(0); break;
        }
    }
    else if (!kappa1 && !iscombined && !kappa2 && alpha){
        switch (ker1) {
            case SL_OP : generator=new BIO_Generator_w_mass<BIOpKernel<LA,SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,alpha);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<LA,SL_OP ...alpha coeff mass matrix=" << alpha << endl; break;
            case DL_OP : generator=new BIO_Generator_w_mass<BIOpKernel<LA,DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,alpha);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<LA,TDL_OP ...alpha coeff mass matrix=" << alpha << endl; break;
            case HS_OP : generator=new BIO_Generator_w_mass<BIOpKernel<LA,HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,alpha);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<LA,HS_OP ...alpha coeff mass matrix=" << alpha << endl; break;
            case TDL_OP : generator=new BIO_Generator_w_mass<BIOpKernel<LA,TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,alpha);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<LA,TDL_OP ...alpha coeff mass matrix=" << alpha << endl; break;
            case CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested BIOpKernel<LA,CST_OP" << endl; ffassert(0); break;
        }
    }
    
    
    // combined kernel
    // alpha coeff mass matrix  --> no define
    // eq HE -> ker1 + coeff2 ker2
    else if (kappa1 && kappa2 && iscombined && alpha) {
        
        
        switch (ker1) {
            case SL_OP :
                switch (ker2) {
                    case SL_OP : generator= new Combined_BIO_Generator<BIOpKernel<HE,SL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<HE,SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,SL_OP ... HE,SL_OP" << endl; break;
                    case DL_OP : generator= new Combined_BIO_Generator<BIOpKernel<HE,SL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<HE,DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,SL_OP ... HE,DL_OP" << endl; break;
                    case HS_OP : generator= new Combined_BIO_Generator<BIOpKernel<HE,SL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<HE,HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,SL_OP ... HE,HS_OP" << endl; break;
                    case TDL_OP : generator= new Combined_BIO_Generator<BIOpKernel<HE,SL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<HE,TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,SL_OP ... HE,TDL_OP" << endl; break;
                    case CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested Combined_BIO_Generator<HE,SL_OP ... HE,CST_OP" << endl; ffassert(0); break;
                }
                break;
            case DL_OP :
                switch (ker2) {
                    case SL_OP : generator= new Combined_BIO_Generator<BIOpKernel<HE,DL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<HE,SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,DL_OP ... HE,SL_OP" << endl; break;
                    case DL_OP : generator= new Combined_BIO_Generator<BIOpKernel<HE,DL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<HE,DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,DL_OP ... HE,DL_OP" << endl; break;
                    case HS_OP : generator= new Combined_BIO_Generator<BIOpKernel<HE,DL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<HE,HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,DL_OP ... HE,HS_OP" << endl; break;
                    case TDL_OP : generator= new Combined_BIO_Generator<BIOpKernel<HE,DL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<HE,TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,DL_OP ... HE,TDL_OP" << endl; break;
                    case CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested Combined_BIO_Generator<HE,DL_OP ... HE,CST_OP" << endl; ffassert(0); break;
                }
                    break;
            case HS_OP :
                switch (ker2) {
                    case SL_OP : generator= new Combined_BIO_Generator<BIOpKernel<HE,HS_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<HE,SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,HS_OP ... HE,SL_OP" << endl; break;
                    case DL_OP : generator= new Combined_BIO_Generator<BIOpKernel<HE,HS_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<HE,DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,HS_OP ... HE,DL_OP" << endl; break;
                    case HS_OP : generator= new Combined_BIO_Generator<BIOpKernel<HE,HS_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<HE,HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,HS_OP ... HE,HS_OP" << endl; break;
                    case TDL_OP : generator= new Combined_BIO_Generator<BIOpKernel<HE,HS_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<HE,TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,HS_OP ... HE,TDL_OP" << endl; break;
                    case CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested Combined_BIO_Generator<HE,HS_OP ... HE,CST_OP" << endl; ffassert(0); break;
                }
                    break;
            case TDL_OP :
                switch (ker2) {
                    case SL_OP : generator= new Combined_BIO_Generator<BIOpKernel<HE,TDL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<HE,SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,TDL_OP ... HE,SL_OP" << endl; break;
                    case DL_OP : generator= new Combined_BIO_Generator<BIOpKernel<HE,TDL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<HE,DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,TDL_OP ... HE,DL_OP" << endl; break;
                    case HS_OP : generator= new Combined_BIO_Generator<BIOpKernel<HE,TDL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<HE,HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,TDL_OP ... HE,HS_OP" << endl; break;
                    case TDL_OP : generator= new Combined_BIO_Generator<BIOpKernel<HE,TDL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<HE,TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,TDL_OP ... HE,TDL_OP" << endl; break;
                    case CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested Combined_BIO_Generator<HE,TDL_OP ... CST_OP" << endl; ffassert(0); break;
                }
                    break;
            case CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested BIOpKernel<HE,CST_OP" << endl; ffassert(0); break;
                
        }
    }
    
    // Eq LAPLACE
    // alpha coeff mass matrix  --> no define
    // eq LA -> ker1 + coeff2 ker2
    else if (!kappa1 && !kappa2 && iscombined && alpha) {
        
        
        switch (ker1) {
            case SL_OP :
                switch (ker2) {
                    case SL_OP : generator= new Combined_BIO_Generator<BIOpKernel<LA,SL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<LA,SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,SL_OP ... LA,SL_OP" << endl; break;
                    case DL_OP : generator= new Combined_BIO_Generator<BIOpKernel<LA,SL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<LA,DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,SL_OP ... LA,DL_OP" << endl; break;
                    case HS_OP : generator= new Combined_BIO_Generator<BIOpKernel<LA,SL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<LA,HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,SL_OP ... LA,HS_OP" << endl; break;
                    case TDL_OP : generator= new Combined_BIO_Generator<BIOpKernel<LA,SL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<LA,TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,SL_OP ... LA,TDL_OP" << endl; break;
                    case CST_OP : if(mpirank == 0 && verbosity>5) cout << " no testedCombined_BIO_Generator<LA,SL_OP ... LA,CST_OP" << endl; ffassert(0); break;
                }
            case DL_OP :
                switch (ker2) {
                    case SL_OP : generator= new Combined_BIO_Generator<BIOpKernel<LA,DL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<LA,SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,DL_OP ... LA,SL_OP" << endl; break;
                    case DL_OP : generator= new Combined_BIO_Generator<BIOpKernel<LA,DL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<LA,DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,DL_OP ... LA,DL_OP" << endl; break;
                    case HS_OP : generator= new Combined_BIO_Generator<BIOpKernel<LA,DL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<LA,HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,DL_OP ... LA,HS_OP" << endl; break;
                    case TDL_OP : generator= new Combined_BIO_Generator<BIOpKernel<LA,DL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<LA,TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,DL_OP ... LA,TDL_OP" << endl; break;
                    case CST_OP :  if(mpirank == 0 && verbosity>5) cout << " no testedCombined_BIO_Generator<LA,DL_OP ... LA,CST_OP" << endl; ffassert(0); break;
                }
            case HS_OP :
                switch (ker2) {
                    case SL_OP : generator= new Combined_BIO_Generator<BIOpKernel<LA,HS_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<LA,SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,HS_OP ... LA,SL_OP" << endl; break;
                    case DL_OP : generator= new Combined_BIO_Generator<BIOpKernel<LA,HS_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<LA,DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,HS_OP ... LA,DL_OP" << endl; break;
                    case HS_OP : generator= new Combined_BIO_Generator<BIOpKernel<LA,HS_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<LA,HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,HS_OP ... LA,HS_OP" << endl; break;
                    case TDL_OP : generator= new Combined_BIO_Generator<BIOpKernel<LA,HS_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<LA,TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,HS_OP ... LA,TDL_OP" << endl; break;
                    case CST_OP :  if(mpirank == 0 && verbosity>5) cout << " no testedCombined_BIO_Generator<LA,HS_OP ... LA,CST_OP" << endl; ffassert(0); break;
                }
            case TDL_OP :
                switch (ker2) {
                    case SL_OP : generator= new Combined_BIO_Generator<BIOpKernel<LA,TDL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<LA,SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,TDL_OP ... LA,SL_OP" << endl; break;
                    case DL_OP : generator= new Combined_BIO_Generator<BIOpKernel<LA,TDL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<LA,DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,TDL_OP ... LA,DL_OP" << endl; break;
                    case HS_OP : generator= new Combined_BIO_Generator<BIOpKernel<LA,TDL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<LA,HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,TDL_OP ... LA,HS_OP" << endl; break;
                    case TDL_OP : generator= new Combined_BIO_Generator<BIOpKernel<LA,TDL_OP,MMesh::RdHat::d+1,P,P>,BIOpKernel<LA,TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappa1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,TDL_OP ... LA,TDL_OP" << endl; break;
                    case CST_OP :  if(mpirank == 0 && verbosity>5) cout << " no testedCombined_BIO_Generator<LA,TDL_OP ... LA,CST_OP" << endl; ffassert(0); break;
                }
            case CST_OP : if(mpirank == 0 && verbosity>5) cout << " no testedCombined_BIO_Generator<LA,CST_OP ... LA,CST_OP" << endl; ffassert(0); break;
        }
    }
    else {
        if(mpirank == 0) cout << "kernel definition error" << endl; ffassert(0);}
   // build the Hmat
   if ( compressor == "" || compressor == "partialACA")
        *Hmat = new HMatrixImpl<R,partialACA>(*generator,p1,false,-1,comm);
    
    else if (compressor == "fullACA")
        *Hmat = new HMatrixImpl<R,fullACA>(*generator,p1,false,-1,comm);
    else if (compressor == "SVD")
        *Hmat = new HMatrixImpl<R,SVD>(*generator,p1,false,-1,comm);
    else {
        cerr << "Error: unknown htool compressor \""+compressor+"\"" << endl;
        ffassert(0);
    }
    delete generator;
    
    
}


template <class R>
void builHmat(HMatrixVirt<R>** Hmat, IMatrix<R>* generatorP,string compressor,vector<htool::R3> &p1,vector<htool::R3> &p2,MPI_Comm comm)  {
    
    if (compressor=="" || compressor == "partialACA")
        *Hmat = new HMatrixImpl<R,partialACA>(*generatorP,p2,p1,-1,comm);
    else if (compressor == "fullACA")
        *Hmat = new HMatrixImpl<R,fullACA>(*generatorP,p2,p1,-1,comm);
    else if (compressor == "SVD")
        *Hmat = new HMatrixImpl<R,SVD>(*generatorP,p2,p1,-1,comm);
    else {
        cerr << "Error: unknown htool compressor \""+compressor+"\"" << endl;
        ffassert(0);
    }
    
}



template <class R, typename P, typename MeshBemtool, class MMesh>
void ff_POT_Generator(HMatrixVirt<R>** Hmat,BemPotential *typePot, Dof<P> &dof, MeshBemtool &mesh, Geometry &node_output, string compressor,vector<htool::R3> &p1,vector<htool::R3> &p2,MPI_Comm comm) {
    
    PotKernelEnum pot = whatTypeEnum(typePot);
    double kappa = typePot->wavenum.real();
    
    switch (pot) {
        case SL_POT :
            if (kappa) {
                Potential<PotKernel<HE,SL_POT,MMesh::RdHat::d+1,P>> t(mesh,kappa);
                POT_Generator<PotKernel<HE,SL_POT,MMesh::RdHat::d+1,P>,P> generator(t,dof,node_output);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func POT_Generator<HE,SL_POT ..." << endl;
                builHmat<R>(Hmat,&generator,compressor,p1,p2,comm);
            }
            else {
                Potential<PotKernel<LA,SL_POT,MMesh::RdHat::d+1,P>> t(mesh,kappa);
                POT_Generator<PotKernel<LA,SL_POT,MMesh::RdHat::d+1,P>,P> generator(t,dof,node_output);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func POT_Generator<LA,SL_POT ..." << endl;
                builHmat<R>(Hmat,&generator,compressor,p1,p2,comm);
            }
            break;
            
        case DL_POT :
            if (kappa) {
                Potential<PotKernel<HE,DL_POT,MMesh::RdHat::d+1,P>> t(mesh,kappa);
                POT_Generator<PotKernel<HE,DL_POT,MMesh::RdHat::d+1,P>,P> generator(t,dof,node_output);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func POT_Generator<HE,DL_POT ..." << endl;
                builHmat<R>(Hmat,&generator,compressor,p1,p2,comm);
            }
            
            else {
                Potential<PotKernel<LA,DL_POT,MMesh::RdHat::d+1,P>> t(mesh,kappa);
                POT_Generator<PotKernel<LA,DL_POT,MMesh::RdHat::d+1,P>,P> generator(t,dof,node_output);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func POT_Generator<LA,DL_POT ..." << endl;
                builHmat<R>(Hmat,&generator,compressor,p1,p2,comm);
            }
            break;
            
        case CST_POT :
            if (kappa) {
                if(mpirank == 0 && verbosity>5) cout << " no POT_Generator<HE,CST_POT ... LA" << endl; ffassert(0); break; }
            else {
                if(mpirank == 0 && verbosity>5) cout << " no POT_Generator<LA,CST_POT ... LA" << endl; ffassert(0); break; }
            
    }
}

// the operator
template<class R,class v_fes1,class v_fes2>
AnyType OpHMatrixtoBEMForm<R,v_fes1,v_fes2>::Op::operator()(Stack stack)  const
{
    typedef typename v_fes1::pfes pfes1;
    typedef typename v_fes2::pfes pfes2;
    typedef typename v_fes1::FESpace FESpace1;
    typedef typename v_fes2::FESpace FESpace2;
    typedef typename FESpace1::Mesh SMesh;
    typedef typename FESpace2::Mesh TMesh;
    
    
    typedef typename std::conditional<SMesh::RdHat::d==1,Mesh1D,Mesh2D>::type MeshBemtool;
    typedef typename std::conditional<SMesh::RdHat::d==1,P1_1D,P1_2D>::type P1;
    
    assert(b && b->nargs);
    const list<C_F0> & largs=b->largs;
    
    // FE space
    pfes1  * pUh= GetAny<pfes1 *>((*b->euh)(stack));
    FESpace1 * Uh = **pUh;
    int NUh =Uh->N;
    pfes2  * pVh= GetAny<pfes2 *>((*b->evh)(stack));
    FESpace2 * Vh = **pVh;
    int NVh =Vh->N;
    ffassert(Vh);
    ffassert(Uh);
    
    int n=Uh->NbOfDF;
    int m=Vh->NbOfDF;
    
    // VFBEM =1 kernel VF   =2 Potential VF
    int VFBEM = typeVFBEM(largs,stack);
    if (mpirank == 0 && verbosity>5)
        cout << "test VFBEM type (1 kernel / 2 potential) "  << VFBEM << endl;
    
 
    HMatrixVirt<R>** Hmat =GetAny<HMatrixVirt<R>** >((*a)(stack));
    
    // info about HMatrix and type solver
    Data_Bem_Solver ds;
    ds.factorize=0;
    ds.initmat=true;
    SetEnd_Data_Bem_Solver<R>(stack,ds, b->nargs,OpCall_FormBilinear_np::n_name_param);  // LIST_NAME_PARM_HMAT
    WhereStackOfPtr2Free(stack)=new StackOfPtr2Free(stack);
    
    // compression infogfg
    SetMaxBlockSize(ds.maxblocksize);
    SetMinClusterSize(ds.minclustersize);
    SetEpsilon(ds.epsilon);
    SetEta(ds.eta);
    SetMinTargetDepth(ds.mintargetdepth);
    SetMinSourceDepth(ds.minsourcedepth);
    
    MPI_Comm comm = ds.commworld ? *(MPI_Comm*)ds.commworld : MPI_COMM_WORLD;
    
    // source/target meshes
    const SMesh & ThU =Uh->Th; // line
    const TMesh & ThV =Vh->Th; // colunm
    bool samemesh = (void*)&Uh->Th == (void*)&Vh->Th;  // same Fem2D::Mesh     +++ pot or kernel
 
    if (VFBEM==1)
        ffassert (samemesh);
     if(init)
        *Hmat =0;
      *Hmat =0;
    if( *Hmat)
            delete *Hmat;
    
    *Hmat =0;
    
    Geometry node; MeshBemtool mesh;
    Mesh2Bemtool(ThU, node, mesh);
    if(mpirank == 0 && verbosity>5)
        std::cout << "Creating dof" << std::endl;
    Dof<P1> dof(mesh,true);
    // now the list of dof is known -> can acces to global num triangle and the local num vertice assiciated
    
    vector<htool::R3> p1(n);
    vector<htool::R3> p2(m);
    Fem2D::R3 pp;
    for (int i=0; i<n; i++) {
        pp = ThU.vertices[i];
        p1[i] = {pp.x, pp.y, pp.z};
    }
    
    
    if(!samemesh) {
        for (int i=0; i<m; i++) {
            pp = ThV.vertices[i];
            p2[i] = {pp.x, pp.y, pp.z};
        }
    }
    else
        p2=p1;
    
    if (VFBEM==1) {
        // info kernel
        pair<BemKernel*,double> kernel = getBemKernel(stack,largs);
        BemKernel *Ker = kernel.first;
        double alpha = kernel.second;
        
        if(mpirank == 0 && verbosity >5) {
            int nk=2;
            for(int i=0;i<nk;i++)
                cout << " kernel info... i: " << i << " typeKernel: " << Ker->typeKernel[i] << " wave number: " << Ker->wavenum[i]  << " coeffcombi: " << Ker->coeffcombi[i] <<endl;
        }
        ff_BIO_Generator<R,P1,SMesh>(Hmat,Ker,dof,alpha,ds.compressor,p1,comm);
    }
    else if (VFBEM==2) {
        BemPotential *Pot = getBemPotential(stack,largs);
        Geometry node_output;
        Mesh2Bemtool(ThV,node_output);
        ff_POT_Generator<R,P1,MeshBemtool,SMesh>(Hmat,Pot,dof,mesh,node_output, ds.compressor, p1,p2,comm);
    }
    return Hmat;
    
}


static void Init_Bem() {
    
    map_type[typeid(const BemFormBilinear *).name( )] = new TypeFormBEM;
    
    map_type[typeid(const BemKFormBilinear *).name( )] = new ForEachType< BemKFormBilinear >;
    map_type[typeid(const BemPFormBilinear *).name( )] = new ForEachType< BemPFormBilinear >;
    
    basicForEachType *t_BEM = atype< const C_args * >( ); //atype< const BemFormBilinear * >( );
    basicForEachType *t_fbem = atype< const BemFormBilinear * >( );
    
    aType t_C_args = map_type[typeid(const C_args *).name( )];
    atype< const C_args * >( )->AddCast(new OneOperatorCode< C_args >(t_C_args, t_fbem) );    // bad
    
    
    typedef  const BemKernel fkernel;
    typedef  const BemPotential fpotential;
    // new type for bem
    typedef const BemKernel *pBemKernel;
    typedef const BemPotential *pBemPotential;
    
    
    Dcl_Type< fkernel * >( );  // a bem kernel
    Dcl_Type< fpotential * >( ); // a bem potential
    Dcl_TypeandPtr< pBemKernel >(0, 0, ::InitializePtr< pBemKernel >, ::DestroyPtr< pBemKernel >,
                                 AddIncrement< pBemKernel >, NotReturnOfthisType);
    // pBemPotential initialize
    Dcl_TypeandPtr< pBemPotential >(0, 0, ::InitializePtr< pBemPotential >, ::DestroyPtr< pBemPotential >,
                                    AddIncrement< pBemPotential >, NotReturnOfthisType);
 
    
    zzzfff->Add("BemKernel", atype< pBemKernel * >( ));
    zzzfff->Add("BemPotential", atype< pBemPotential * >( ));
    
     // pBemKernel initialize
    atype<pBemKernel>()->AddCast( new E_F1_funcT<pBemKernel,pBemKernel*>(UnRef<pBemKernel>));
    // BemPotential
    atype<pBemPotential>()->AddCast( new E_F1_funcT<pBemPotential,pBemPotential*>(UnRef<pBemPotential>));
   
    
    // simplified type/function to define varf bem
    Dcl_Type< const FoperatorKBEM * >( );
    Dcl_Type< const FoperatorPBEM * >( );
    Dcl_Type<std::map<std::string, std::string>*>( );
        
    TheOperators->Add("<<",new OneBinaryOperator<PrintPinfos<std::map<std::string, std::string>*>>);
    Add<std::map<std::string, std::string>*>("[","",new OneOperator2_<string*, std::map<std::string, std::string>*, string*>(get_info));

    addHmat<double>();
    addHmat<std::complex<double>>();

    //BemKernel
    TheOperators->Add("<-", new OneOperator4_<pBemKernel*,pBemKernel*,string*,std::complex<double>,std::complex<double> >(&initKernel_Helmholtz));
    TheOperators->Add("<-", new OneOperator3_<pBemKernel*,pBemKernel*,string*,std::complex<double> >(&initKernel_Laplace));
    TheOperators->Add("<-", new OneOperator2_<pBemKernel*,pBemKernel*,string* >(&initKernel_default));
    //BemPotential
    TheOperators->Add("<-", new OneOperator3_<pBemPotential*,pBemPotential*,string*,std::complex<double> >(&initPotential_Helmholtz));
    TheOperators->Add("<-", new OneOperator2_<pBemPotential*,pBemPotential*,string* >(&initPotential_default));
    
    zzzfff->Add("HMatrix", atype<HMatrixVirt<double> **>());
    map_type_of_map[make_pair(atype<HMatrixVirt<double>**>(), atype<double*>())] = atype<HMatrixVirt<double>**>();
    map_type_of_map[make_pair(atype<HMatrixVirt<double>**>(), atype<Complex*>())] = atype<HMatrixVirt<std::complex<double> >**>();
        
    TheOperators->Add("<-", new OpHMatrixtoBEMForm< std::complex<double>, v_fesS, v_fesS > (1) );
    TheOperators->Add("=", new OpHMatrixtoBEMForm< std::complex<double>, v_fesS, v_fesS > );
    TheOperators->Add("<-", new OpHMatrixtoBEMForm< std::complex<double>, v_fesL, v_fesL > (1) );
    TheOperators->Add("=", new OpHMatrixtoBEMForm< std::complex<double>, v_fesL, v_fesL > );
       

    TheOperators->Add("<-", new OpHMatrixtoBEMForm< std::complex<double>, v_fesL, v_fes > (1) );
    TheOperators->Add("=", new OpHMatrixtoBEMForm< std::complex<double>, v_fesL, v_fes > );
    TheOperators->Add("<-", new OpHMatrixtoBEMForm< std::complex<double>, v_fesL, v_fesS > (1) );
    TheOperators->Add("=", new OpHMatrixtoBEMForm< std::complex<double>, v_fesL, v_fesS > );

    TheOperators->Add("<-", new OpHMatrixtoBEMForm< std::complex<double>, v_fesS, v_fes > (1));
    TheOperators->Add("=", new OpHMatrixtoBEMForm< std::complex<double>, v_fesS, v_fes > );

    // BemKernel
    
    Dcl_Type<listBemKernel> ();
    TheOperators->Add("+",new OneBinaryOperator_st< Op_addBemKernel<listBemKernel,pBemKernel,pBemKernel> >);
    TheOperators->Add("+",new OneBinaryOperator_st< Op_addBemKernel<listBemKernel,listBemKernel,pBemKernel> >);
    TheOperators->Add("=",new OneBinaryOperator_st< Op_setBemKernel<false,pBemKernel*,pBemKernel*,listBemKernel> >);
    TheOperators->Add("<-", new OneBinaryOperator_st< Op_setBemKernel<true,pBemKernel*,pBemKernel*,listBemKernel> >);
      
    Dcl_Type< const CBemDomainOfIntegration * >( );
    Dcl_Type< const CPartBemDI * >( );
    
    Add< const CPartBemDI * >("(", "", new OneOperatorCode< CBemDomainOfIntegration >);
    
    Add< const CBemDomainOfIntegration * >("(", "", new OneOperatorCode< BemKFormBilinear >);
    Add< const CDomainOfIntegration * >("(", "", new OneOperatorCode< BemPFormBilinear >);
    
   // zzzfff->AddF("varfbem", t_BEM);
    
    Global.Add("BEM","(",new FormalKBEMcode);
    Global.Add("POT","(",new FormalPBEMcode);
    // 2D / 1D Bem
    Global.Add("int2dx2d","(",new OneOperatorCode<CPartBemDI2d2d>);
    Global.Add("int1dx1d","(",new OneOperatorCode<CPartBemDI1d1d>);
    Global.Add("int1dx2d","(",new OneOperatorCode<CPartBemDI1d2d>);
    Global.Add("int2dx1d","(",new OneOperatorCode<CPartBemDI2d1d>);
    
}

LOADFUNC(Init_Bem)

