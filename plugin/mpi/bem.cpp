//ff-c++-LIBRARY-dep: cxx11 [mkl|blas] mpi pthread htool bemtool boost
//ff-c++-cpp-dep:
// for def  M_PI under windows in <cmath>
#define _USE_MATH_DEFINES

#define BOOST_NO_CXX17_IF_CONSTEXPR
#include <ff++.hpp>
#include <AFunction_ext.hpp>
#include <lgfem.hpp>
#include <R3.hpp>

#include <htool/htool.hpp>

// include the bemtool library .... path define in where library
//#include <bemtool/operator/block_op.hpp>
#include <bemtool/tools.hpp>
#include <bemtool/fem/dof.hpp>
#include <bemtool/operator/operator.hpp>
#include <bemtool/miscellaneous/htool_wrap.hpp>
#include "PlotStream.hpp"

#include "common.hpp"

extern FILE *ThePlotStream;

using namespace std;
using namespace htool;
using namespace bemtool;

#include <type_traits>
#include "bem.hpp"
  
double ff_htoolEta=10., ff_htoolEpsilon=1e-3;
long ff_htoolMinclustersize=10, ff_htoolMaxblocksize=1000000, ff_htoolMintargetdepth=0, ff_htoolMinsourcedepth=0;

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
    virtual MPI_Comm get_comm() const = 0;
    virtual int get_rankworld() const = 0;
    virtual int get_sizeworld() const = 0;
    virtual int get_local_size() const = 0;
    virtual int get_local_offset() const = 0;
    virtual const std::vector<SubMatrix<K>*>& get_MyNearFieldMats() const = 0;
    virtual const LowRankMatrix<K>& get_MyFarFieldMats(int i) const = 0;
    virtual int get_MyFarFieldMats_size() const = 0;
    virtual const std::vector<SubMatrix<K>*>& get_MyStrictlyDiagNearFieldMats() const = 0;
    virtual Matrix<K> get_local_dense() const = 0;
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
    virtual void build(VirtualGenerator<K> &mat, double* xt, double* xs) = 0;
    virtual void build(VirtualGenerator<K> &mat, double* xt)             = 0;

    virtual ~HMatrixVirt() {};
};


template<class K>
class HMatrixImpl : public HMatrixVirt<K> {
private:
    HMatrix<K> H;
public:
    HMatrixImpl(std::shared_ptr<VirtualCluster> t, std::shared_ptr<VirtualCluster> s, double epsilon=1e-6 ,double eta=10, char symmetry='N', char UPLO='N',const int& reqrank=-1, MPI_Comm comm=MPI_COMM_WORLD) : H(t,s,epsilon,eta,symmetry,UPLO,reqrank,comm){}
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
    const std::vector<SubMatrix<K>*>& get_MyNearFieldMats() const {return H.get_MyNearFieldMats();}
    const LowRankMatrix<K>& get_MyFarFieldMats(int i) const {return *(H.get_MyFarFieldMats()[i]);}
    int get_MyFarFieldMats_size() const {return H.get_MyFarFieldMats().size();}
    const std::vector<SubMatrix<K>*>& get_MyStrictlyDiagNearFieldMats() const {return H.get_MyStrictlyDiagNearFieldMats();}
    Matrix<K> get_local_dense() const {return H.get_local_dense();}
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
    void build(VirtualGenerator<K> &mat, const std::vector<htool::R3> &xt, const std::vector<double> &rt, const std::vector<int> &tabt, const std::vector<double> &gt, const std::vector<htool::R3> &xs, const std::vector<double> &rs, const std::vector<int> &tabs, const std::vector<double> &gs) {H.build(mat,xt,tabt,xs,tabs);}


    void build(VirtualGenerator<K> &mat, double* xt, double* xs) {H.build(mat,xt,xs);}
    void build(VirtualGenerator<K> &mat, double* xt) {H.build(mat,xt);}                       
};


template<class Type, class K>
AnyType To(Stack stack,Expression emat,Expression einter,int init)
{
    ffassert(einter);
    HMatrixVirt<K>** Hmat = GetAny<HMatrixVirt<K>** >((*einter)(stack));
    ffassert(Hmat && *Hmat);
    HMatrixVirt<K>& H = **Hmat;
    if(std::is_same<Type, KNM<K>>::value) {
        Matrix<K> mdense = H.get_local_dense();
        KNM<K>* M = GetAny<KNM<K>*>((*emat)(stack));
        for (int i=0; i< H.nb_rows(); i++)
            for (int j=0; j< H.nb_cols(); j++)
                (*M)(i,j) = 0;
        for (int i=0; i< mdense.nb_rows(); i++)
            for (int j=0; j< mdense.nb_cols(); j++)
                (*M)(H.get_permt(i+H.get_local_offset()),H.get_perms(j)) = mdense(i,j);

        return M;
    }
    else {
        ffassert(0);
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
                    const LowRankMatrix<K>& l = (*H)->get_MyFarFieldMats(i);
                    buflr[5*i] = l.get_offset_i();
                    buflr[5*i+1] = l.get_offset_j();
                    buflr[5*i+2] = l.nb_rows();
                    buflr[5*i+3] = l.nb_cols();
                    buflr[5*i+4] = l.rank_of();
                    bufcomp[i] = l.space_saving();
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
        int niterx=3000;
        bool res=fgmres(A,P,1,(K*)*u,(K*)*out,eps,niterx,niterx,(mpirank==0)*verbosity);
        //bool res=fgmres(A,P,1,(K*)*u,(K*)*out,eps,niterx,200,(mpirank==0)*verbosity);
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
class MyMatrix: public VirtualGenerator<K>{
        const KNM<K> &M;

public:
        MyMatrix(const KNM<K> &mat):VirtualGenerator<K>(mat.N(),mat.M()),M(mat) {}

        K get_coef(const int& i, const int& j)const {return M(i,j);}
        void copy_submatrix(int m, int n, const int *const rows, const int *const cols, K *ptr) const {
            std::fill_n(ptr,m*n,0);
            for (int i=0;i<m;i++){
                for (int j=0;j<n;j++){
                    ptr[i+m*j]=M(rows[i],cols[j]);
                }

            }
        }
};

template<class K, int init>
AnyType SetCompressMat(Stack stack,Expression emat,Expression einter)
{ return SetCompressMat<K>(stack,emat,einter,init);}

template<class K>
AnyType SetCompressMat(Stack stack,Expression emat,Expression einter,int init)
{
  HMatrixVirt<K>** Hmat =GetAny<HMatrixVirt<K>** >((*emat)(stack));
  const typename CompressMat<K>::Op * mi(dynamic_cast<const typename CompressMat<K>::Op *>(einter));

  double epsilon=mi->arg(0,stack,ff_htoolEpsilon);
  pcommworld pcomm=mi->argc(1,stack,nullptr);
  double eta=mi->arg(2,stack,ff_htoolEta);
  int minclustersize=mi->argl(3,stack,ff_htoolMinclustersize);
  int maxblocksize=mi->argl(4,stack,ff_htoolMaxblocksize);
  int mintargetdepth=mi->argl(5,stack,ff_htoolMintargetdepth);
  int minsourcedepth=mi->argl(6,stack,ff_htoolMinsourcedepth);
  string* pcompressor=mi->args(7,stack,0);

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

  std::vector<double> p(3*xx.n);
  for (int i=0; i<xx.n; i++) {
    p[3*i+0] = xx[i];
    p[3*i+1] = yy[i];
    p[3*i+2] = zz[i];
  }

  std::shared_ptr<Cluster<PCARegularClustering>> t = std::make_shared<Cluster<PCARegularClustering>>();
  t->set_minclustersize(minclustersize);
  t->build(xx.n,p.data(),2,comm);

  //cout << M.N() << " " << xx.M() << " " << yy.n << " " << zz.n<< endl;
  if (init) delete *Hmat;
    
  *Hmat = new HMatrixImpl<K>(t,t,epsilon,eta);
  std::shared_ptr<htool::VirtualLowRankGenerator<K>> LowRankGenerator = nullptr;

  if ( compressor == "" || compressor == "partialACA")
    LowRankGenerator = std::make_shared<htool::partialACA<K>>();
   else if (compressor == "fullACA")
    LowRankGenerator = std::make_shared<htool::fullACA<K>>();
   else if (compressor == "SVD")
    LowRankGenerator = std::make_shared<htool::SVD<K>>();
   else {
       cerr << "Error: unknown htool compressor \""+compressor+"\"" << endl;
       ffassert(0);
   }
  (*Hmat)->set_compression(LowRankGenerator);
  (*Hmat)->set_maxblocksize(maxblocksize);
  (*Hmat)->set_mintargetdepth(mintargetdepth);
  (*Hmat)->set_minsourcedepth(minsourcedepth);
  (*Hmat)->build(A,p.data());
  return Hmat;
}

template<class K>
void addHmat() {
   // Dcl_Type<HMatrixVirt<K>**>(Initialize<HMatrixVirt<K>*>, Delete<HMatrixVirt<K>*>);
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
    Dcl_Type<const typename CompressMat<K>::Op *>();
    //Add<const typename assembleHMatrix<LR, K>::Op *>("<-","(", new assembleHMatrix<LR, K>);

    TheOperators->Add("=",
    new OneOperator2_<HMatrixVirt<K>**,HMatrixVirt<K>**,const typename CompressMat<K>::Op*,E_F_StackF0F0>(SetCompressMat<K, 1>));
    TheOperators->Add("<-",
    new OneOperator2_<HMatrixVirt<K>**,HMatrixVirt<K>**,const typename CompressMat<K>::Op*,E_F_StackF0F0>(SetCompressMat<K, 0>));

    Global.Add("compress","(",new CompressMat<K>);
}

template<class R, class MMesh, class v_fes1, class v_fes2>
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



template <class R>
void builHmat(HMatrixVirt<R>** Hmat, VirtualGenerator<R>* generatorP,const Data_Bem_Solver& data,vector<double> &p1,vector<double> &p2,MPI_Comm comm)  {
    std::shared_ptr<Cluster<PCARegularClustering>> t =std::make_shared<Cluster<PCARegularClustering>>();
    std::shared_ptr<Cluster<PCARegularClustering>> s =std::make_shared<Cluster<PCARegularClustering>>();
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

// the operator
template<class R,class MMesh, class v_fes1,class v_fes2>
AnyType OpHMatrixtoBEMForm<R,MMesh,v_fes1,v_fes2>::Op::operator()(Stack stack)  const
{
    typedef typename v_fes1::pfes pfes1;
    typedef typename v_fes2::pfes pfes2;
    typedef typename v_fes1::FESpace FESpace1;
    typedef typename v_fes2::FESpace FESpace2;
    typedef typename FESpace1::Mesh SMesh;
    typedef typename FESpace2::Mesh TMesh;
    typedef typename SMesh::RdHat SRdHat;
    typedef typename TMesh::RdHat TRdHat;
    
    
    typedef typename std::conditional<SMesh::RdHat::d==1,Mesh1D,Mesh2D>::type MeshBemtool;
    typedef typename std::conditional<SMesh::RdHat::d==1,P0_1D,P0_2D>::type P0;
    typedef typename std::conditional<SMesh::RdHat::d==1,P1_1D,P1_2D>::type P1;
    typedef typename std::conditional<SMesh::RdHat::d==1,P2_1D,P2_2D>::type P2;

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

/*
    for (int i=0; i<ThU.nt;i++) {
        auto &K = ThU[i];
        for (int e=0;e<3;e++) {
            R signe = K.EdgeOrientation(e) ;
            cout << i << " " << e << " " << signe << endl;
        }
    }
*/

    vector<double> p1(3*n);
    vector<double> p2(3*m);
    Fem2D::R3 pp;
    bemtool::R3 p;
    SRdHat pbs;
    TRdHat pbt;
    pbs[0] = 1./(SRdHat::d+1);
    pbs[1] = 1./(SRdHat::d+1);
    if (SRdHat::d == 2) pbs[2] = 1./(SRdHat::d+1);
    pbt[0] = 1./(TRdHat::d+1);
    pbt[1] = 1./(TRdHat::d+1);
    if (TRdHat::d == 2) pbt[2] = 1./(TRdHat::d+1);

    int Snbv = Uh->TFE[0]->ndfonVertex;
    int Snbe = Uh->TFE[0]->ndfonEdge;
    int Snbt = Uh->TFE[0]->ndfonFace;
    bool SP0 = SRdHat::d == 1 ? (Snbv == 0) && (Snbe == 1) && (Snbt == 0) : (Snbv == 0) && (Snbe == 0) && (Snbt == 1);
    bool SP1 = (Snbv == 1) && (Snbe == 0) && (Snbt == 0);
    bool SP2 = (Snbv == 1) && (Snbe == 1) && (Snbt == 0);
    bool SRT0 = (SRdHat::d == 2) && (Snbv == 0) && (Snbe == 1) && (Snbt == 0);

    if (SP2) {
        Dof<P2> dof(mesh,true);
        for (int i=0; i<n; i++) {
            const std::vector<N2>& jj = dof.ToElt(i);
            p = dof(jj[0][0])[jj[0][1]];
            p1[3*i+0] = p[0];
            p1[3*i+1] = p[1];
            p1[3*i+2] = p[2];
        }
    }
    else if (SRT0) {
        Dof<RT0_2D> dof(mesh);
        for (int i=0; i<n; i++) {
            const std::vector<N2>& jj = dof.ToElt(i);
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

    std::shared_ptr<Cluster<PCARegularClustering>> t =std::make_shared<Cluster<PCARegularClustering>>();
    std::shared_ptr<Cluster<PCARegularClustering>> s = make_shared<Cluster<PCARegularClustering>>();
    t->set_minclustersize(ds.minclustersize);
    s->set_minclustersize(ds.minclustersize);
    t->build(n,p1.data(),2,comm);
    // 
    // In this method :OpHMatrixtoBEMForm<R,MMesh,v_fes1,v_fes2>::Op::operator() 
    // be careful version du 30 mai 2022: - s : correponding to target
    //                                    - t : corresponding to source
    //     

    if(!samemesh) {
        if( Vh->TFE[0]->N == 1){
            // case the targer FE is scalar
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
            int nnn= Vh->TFE[0]->N;

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
                //if( i%2 ==0) cout << i << ",  pp.x =" << pp.x << endl;
                for(int iii=0; iii<nnn; iii++){
                    ffassert( nnn*3*i+3*iii+2 < nnn*3*m );
                    //if(i== m-1) cout << "pp.x =" << pp.x << endl;
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
    
    if (VFBEM==1) {
        // info kernel
        pair<BemKernel*,double> kernel = getBemKernel(stack,largs);
        BemKernel *Ker = kernel.first;
        double alpha = kernel.second;
        
        if(mpirank == 0 && verbosity >5) {
            int nk=-1;
            iscombinedKernel(Ker) ? nk=2 : nk=1;
            for(int i=0;i<nk;i++)
                cout << " kernel info... i: " << i << " typeKernel: " << Ker->typeKernel[i] << " wave number: " << Ker->wavenum[i]  << " coeffcombi: " << Ker->coeffcombi[i] <<endl;
        }
        VirtualGenerator<R>* generator;

        if(mpirank == 0 && verbosity>5)
            std::cout << "Creating dof" << std::endl;

        if (SP0) {
            Dof<P0> dof(mesh);
            ff_BIO_Generator<R,P0,SMesh>(generator,Ker,dof,alpha);
        }
        else if (SP1) {
            Dof<P1> dof(mesh,true);
            ff_BIO_Generator<R,P1,SMesh>(generator,Ker,dof,alpha);
        }
        else if (SP2) {
            Dof<P2> dof(mesh,true);
            ff_BIO_Generator<R,P2,SMesh>(generator,Ker,dof,alpha);
        }
        else if (SRT0 && SRdHat::d == 2) {
            ff_BIO_Generator_Maxwell<R>(generator,Ker,mesh,alpha);
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
        Geometry node_output;
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

        VirtualGenerator<R>* generator;
        if (SP0) {
            Dof<P0> dof(mesh);
            ff_POT_Generator<R,P0,MeshBemtool,SMesh>(generator,Pot,dof,mesh,node_output);
        }
        else if (SP1) {
            Dof<P1> dof(mesh,true);
            ff_POT_Generator<R,P1,MeshBemtool,SMesh>(generator,Pot,dof,mesh,node_output);
        }
        else if (SP2) {
            Dof<P2> dof(mesh,true);
            ff_POT_Generator<R,P2,MeshBemtool,SMesh>(generator,Pot,dof,mesh,node_output);
        }
        else if (SRT0 && SRdHat::d == 2) {
            
            ff_POT_Generator_Maxwell<R,bemtool::RT0_2D>(generator,Pot,mesh,node_output);
        }
        else
            ffassert(0);
        Assembly(Hmat,generator,ds,p1,p2,comm,MMesh::RdHat::d+1);
        delete generator;
    }
    return Hmat;
    
}

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






EquationEnum whatEquationEnum(BemKernel *K,int i) {
    int nk=2; // restriction, max 2 kernels for combined
    int type[2]={-1,-1};
    EquationEnum equation[3] = {LA,HE,YU};
    double wavenumRe, wavenumImag;
        
    for (int i=0;i<nk;i++) {
      wavenumRe=K->wavenum[i].real(); wavenumImag=K->wavenum[i].imag();
      if(wavenumRe==0 && wavenumImag==0) type[i] = 0;
      else if(wavenumRe>0 && wavenumImag==0) type[i] = 1;
      else if(wavenumRe==0 && wavenumImag>0) type[i] = 2;
    }
    if(type[0]-type[1]) {
      cerr << "Error, the combined kernels must be the same equation " << equation[type[0]] << " " << equation[type[1]] << endl;
      throw(ErrorExec("exit",1));}
    
    const EquationEnum equEnum = equation[0];
    return equEnum;
    
}
//  Begin Array of HMatrix[int]
template< class A >
inline AnyType DestroyKNmat(Stack, const AnyType &x) {
  KN< A > *a = GetAny< KN< A > * >(x);
  for (int i = 0; i < a->N( ); i++)
    if ((*a)[i]) delete (*a)[i];
  a->destroy( );
  return Nothing;
}

template< class RR, class A, class B >
RR *get_elementp_(const A &a, const B &b) {
  if (b < 0 || a->N( ) <= b) {
    cerr << " Out of bound  0 <=" << b << " < " << a->N( ) << " array type = " << typeid(A).name( )
         << endl;
    ExecError("Out of bound in operator []");
  }
  return &((*a)[b]);
}

template< class R >
R *set_initinit(R *const &a, const long &n) {
  SHOWVERB(cout << " set_init " << typeid(R).name( ) << " " << n << endl);
  a->init(n);
  for (int i = 0; i < n; i++) (*a)[i] = 0;
  return a;
}

template<class R>
void ArrayofHmat()
{
    typedef HMatrixVirt< R > *Mat;
    typedef Mat *PMat;
    typedef KN< Mat > AMat;

    Dcl_Type< AMat * >(0, ::DestroyKNmat< Mat >);
    // to declare HMatrix[int]
    map_type_of_map[make_pair(atype< long >( ), atype< PMat >( )->right())] = atype< AMat * >( );
    atype<  AMat * >( )->Add(
      "[", "",
      new OneOperator2_<  PMat, AMat *, long >(get_elementp_< Mat, AMat *, long >));

    TheOperators->Add("<-", new OneOperator2_< AMat*, AMat* , long >(&set_initinit));
  

    // resize mars 2006 v2.4-1
    Dcl_Type< Resize< AMat > >( );
    Add< AMat *>("resize", ".",
                         new OneOperator1< Resize< AMat >, AMat * >(to_Resize));
    Add< Resize< AMat > >(
      "(", "", new OneOperator2_< AMat*, Resize< AMat >, long >(resizeandclean1));

}
//  End Array of HMatrix[int]

template<class R, class v_fes1, class v_fes2>
class OpHMatrixUser : public OneOperator
{
    public:
        typedef typename v_fes1::pfes pfes1;
        typedef typename v_fes2::pfes pfes2;
        class Op : public E_F0info {
            public:
                Expression g, uh1, uh2;
                static const int n_name_param = 8;
                static basicAC_F0::name_and_type name_param[] ;
                Expression nargs[n_name_param];
                long argl(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
                string* args(int i,Stack stack,string* a) const{ return nargs[i] ? GetAny<string*>( (*nargs[i])(stack) ): a;}
                double arg(int i,Stack stack,double a) const{ return nargs[i] ? GetAny<double>( (*nargs[i])(stack) ): a;}
                pcommworld argc(int i,Stack stack,pcommworld a ) const{ return nargs[i] ? GetAny<pcommworld>( (*nargs[i])(stack) ): a;}
                Op(const basicAC_F0 &  args, Expression  bb, Expression  cc, Expression  dd) : g(bb),uh1(cc), uh2(dd) {
                    args.SetNameParam(n_name_param,name_param,nargs);
                }
        };

        E_F0 * code(const basicAC_F0 & args) const {
            return  new Op(args,t[0]->CastTo(args[0]),t[1]->CastTo(args[1]),t[2]->CastTo(args[2]));
        }

        OpHMatrixUser() :
        OneOperator(atype<const typename OpHMatrixUser<R,v_fes1,v_fes2>::Op*>(),atype<VirtualGenerator<R>**>(),atype<pfes1*>(),atype<pfes2*>()) {}
};

template<class K, class v_fes1, class v_fes2>
basicAC_F0::name_and_type  OpHMatrixUser<K,v_fes1,v_fes2>::Op::name_param[]= {
  {  "eps", &typeid(double)},
  {  "commworld", &typeid(pcommworld)},
  {  "eta", &typeid(double)},
  {  "minclustersize", &typeid(long)},
  {  "maxblocksize", &typeid(long)},
  {  "mintargetdepth", &typeid(long)},
  {  "minsourcedepth", &typeid(long)},
  {  "compressor", &typeid(string*)}
};

template<class R, class v_fes1,class v_fes2, int init>
AnyType SetOpHMatrixUser(Stack stack,Expression emat, Expression eop)
{
    typedef typename v_fes1::pfes pfes1;
    typedef typename v_fes2::pfes pfes2;
    typedef typename v_fes1::FESpace FESpace1;
    typedef typename v_fes2::FESpace FESpace2;
    typedef typename FESpace1::Mesh SMesh;
    typedef typename FESpace2::Mesh TMesh;
    typedef typename SMesh::RdHat SRdHat;
    typedef typename TMesh::RdHat TRdHat;

    typedef typename OpHMatrixUser<R,v_fes1,v_fes2>::Op UOp;
    const UOp * op(dynamic_cast<const UOp *>(eop));
    pfes1  * pUh= GetAny<pfes1 *>((*op->uh1)(stack));
    FESpace1 * Uh = **pUh;
    int NUh =Uh->N;
    pfes2  * pVh= GetAny<pfes2 *>((*op->uh2)(stack));
    FESpace2 * Vh = **pVh;
    int NVh =Vh->N;
    ffassert(Vh);
    ffassert(Uh);

    int n=Uh->NbOfDF;
    int m=Vh->NbOfDF;

    HMatrixVirt<R>** Hmat =GetAny<HMatrixVirt<R>** >((*emat)(stack));

    Data_Bem_Solver ds;
    ds.factorize=0;
    ds.initmat=true;

    ds.epsilon = op->arg(0,stack,ds.epsilon);
    ds.commworld = op->argc(1,stack,ds.commworld);
    ds.eta = op->arg(2,stack,ds.eta);
    ds.minclustersize = op->argl(3,stack,ds.minclustersize);
    ds.maxblocksize = op->argl(4,stack,ds.maxblocksize);
    ds.mintargetdepth = op->argl(5,stack,ds.mintargetdepth);
    ds.minsourcedepth = op->argl(6,stack,ds.minsourcedepth);
    ds.compressor = *(op->args(7,stack,&ds.compressor));

    const SMesh & ThU =Uh->Th;
    const TMesh & ThV =Vh->Th;
    bool samemesh = (void*)&Uh->Th == (void*)&Vh->Th;

     if(init)
        *Hmat =0;
      *Hmat =0;
    if( *Hmat)
            delete *Hmat;

    *Hmat =0;

    vector<double> p1(3*n);
    vector<double> p2(3*m);
    Fem2D::R3 pp;
    bemtool::R3 p;
    SRdHat pbs;
    TRdHat pbt;
    pbs[0] = 1./(SRdHat::d+1);
    pbs[1] = 1./(SRdHat::d+1);
    if (SRdHat::d == 2) pbs[2] = 1./(SRdHat::d+1);
    pbt[0] = 1./(TRdHat::d+1);
    pbt[1] = 1./(TRdHat::d+1);
    if (TRdHat::d == 2) pbt[2] = 1./(TRdHat::d+1);

    int Snbv = Uh->TFE[0]->ndfonVertex;
    int Snbe = Uh->TFE[0]->ndfonEdge;
    int Snbt = Uh->TFE[0]->ndfonFace;
    bool SP0 = SRdHat::d == 1 ? (Snbv == 0) && (Snbe == 1) && (Snbt == 0) : (Snbv == 0) && (Snbe == 0) && (Snbt == 1);
    bool SP1 = (Snbv == 1) && (Snbe == 0) && (Snbt == 0);

    int Tnbv = Vh->TFE[0]->ndfonVertex;
    int Tnbe = Vh->TFE[0]->ndfonEdge;
    int Tnbt = Vh->TFE[0]->ndfonFace;
    bool TP0 = TRdHat::d == 1 ? (Tnbv == 0) && (Tnbe == 1) && (Tnbt == 0) : (Tnbv == 0) && (Tnbe == 0) && (Tnbt == 1);
    bool TP1 = (Tnbv == 1) && (Tnbe == 0) && (Tnbt == 0);

    for (int i=0; i<n; i++) {
        if (SP1)
            pp = ThU.vertices[i];
        else if (SP0)
            pp = ThU[i](pbs);
        else {
            if (mpirank == 0) std::cerr << "ff-Htool error: only P0 and P1 discretizations are available for now." << std::endl;
            ffassert(0);
        }
        p1[3*i+0] = pp.x;
        p1[3*i+1] = pp.y;
        p1[3*i+2] = pp.z;
    }

    if(!samemesh) {
        for (int i=0; i<m; i++) {
            if (TP1)
                pp = ThV.vertices[i];
            else if (TP0)
                pp = ThV[i](pbt);
            else {
                if (mpirank == 0) std::cerr << "ff-Htool error: only P0 and P1 discretizations are available for now." << std::endl;
                ffassert(0);
            }
            p2[3*i+0] = pp.x;
            p2[3*i+1] = pp.y;
            p2[3*i+2] = pp.z;
        }
    }
    else{
        p2=p1;
    }

    VirtualGenerator<R>** generator = GetAny<VirtualGenerator<R>**>((*op->g)(stack));

    MPI_Comm comm = ds.commworld ? *(MPI_Comm*)ds.commworld : MPI_COMM_WORLD;
    builHmat<R>(Hmat,*generator,ds,p1,p2,comm);

    return Hmat;
}

template<class R, class v_fes1, class v_fes2>
void AddHMatrixUser() {
    typedef typename OpHMatrixUser< R, v_fes1, v_fes2 >::Op OpT;
    Dcl_Type<const OpT*>();
    Global.Add("Build","(", new OpHMatrixUser< R, v_fes1, v_fes2 >);
    Add<const OpT*>("<-","(", new OpHMatrixUser< R, v_fes1, v_fes2 >);
    TheOperators->Add("<-", new OneOperator2_<HMatrixVirt<R>**,HMatrixVirt<R>**,const OpT*,E_F_StackF0F0>(SetOpHMatrixUser<R,v_fes1,v_fes2,0>));
    TheOperators->Add("=", new OneOperator2_<HMatrixVirt<R>**,HMatrixVirt<R>**,const OpT*,E_F_StackF0F0>(SetOpHMatrixUser<R,v_fes1,v_fes2,1>));
}

static void Init_Bem() {
    if(mpirank==0) cout << "\nInit_Bem\n";

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
    
    
   // Dcl_Type< fkernel * >( );  // a bem kernel
   // Dcl_Type< fpotential * >( ); // a bem potential
    Dcl_Type< const OP_MakeBemKernelFunc::Op * >( );
    Dcl_Type< const OP_MakeBemPotentialFunc::Op * >( );
    
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
    TheOperators->Add("<-", new OneOperatorCode< OP_MakeBemKernel >);
    //BemPotential
    TheOperators->Add("<-", new OneOperatorCode< OP_MakeBemPotential >);
   
    zzzfff->Add("HMatrix", atype<HMatrixVirt<double> **>());
    map_type_of_map[make_pair(atype<HMatrixVirt<double>**>(), atype<double*>())] = atype<HMatrixVirt<double>**>();
    map_type_of_map[make_pair(atype<HMatrixVirt<double>**>(), atype<Complex*>())] = atype<HMatrixVirt<std::complex<double> >**>();
    // bem integration space/target space must be review Axel 08/2020    
    TheOperators->Add("<-", new OpHMatrixtoBEMForm< std::complex<double>, MeshS, v_fesS, v_fesS > (1) );
    TheOperators->Add("=", new OpHMatrixtoBEMForm< std::complex<double>, MeshS, v_fesS, v_fesS > );
    TheOperators->Add("<-", new OpHMatrixtoBEMForm< std::complex<double>, MeshS, v_fesS, v_fes3 > (1) );
    TheOperators->Add("=", new OpHMatrixtoBEMForm< std::complex<double>, MeshS, v_fesS, v_fes3 > );
    TheOperators->Add("<-", new OpHMatrixtoBEMForm< std::complex<double>, MeshL, v_fesL, v_fesL > (1) );
    TheOperators->Add("=", new OpHMatrixtoBEMForm< std::complex<double>, MeshL, v_fesL, v_fesL > );
       
    TheOperators->Add("<-", new OpHMatrixtoBEMForm< std::complex<double>, MeshL, v_fesL, v_fes > (1) );
    TheOperators->Add("=", new OpHMatrixtoBEMForm< std::complex<double>, MeshL, v_fesL, v_fes > );
    TheOperators->Add("<-", new OpHMatrixtoBEMForm< std::complex<double>, MeshL, v_fesL, v_fesS > (1) );
    TheOperators->Add("=", new OpHMatrixtoBEMForm< std::complex<double>, MeshL, v_fesL, v_fesS > );

    TheOperators->Add("<-", new OpHMatrixtoBEMForm< std::complex<double>, MeshS, v_fesS, v_fes > (1));
    TheOperators->Add("=", new OpHMatrixtoBEMForm< std::complex<double>, MeshS, v_fesS, v_fes > );

    // operation on BemKernel
    Dcl_Type<listBemKernel> ();
    TheOperators->Add("+",new OneBinaryOperator_st< Op_addBemKernel<listBemKernel,pBemKernel,pBemKernel> >);
    //TheOperators->Add("+",new OneBinaryOperator_st< Op_addBemKernel<listBemKernel,listBemKernel,pBemKernel> >); // no need is the combinaison is only with 2 kernels
    TheOperators->Add("=",new OneBinaryOperator_st< Op_setBemKernel<false,pBemKernel*,pBemKernel*,listBemKernel> >);
    TheOperators->Add("<-", new OneBinaryOperator_st< Op_setBemKernel<true,pBemKernel*,pBemKernel*,listBemKernel> >);

    TheOperators->Add("<-", new OneOperator2_<pBemKernel*,pBemKernel*,pBemKernel >(&set_copy_incr));
    TheOperators->Add("=", new OneOperator2<pBemKernel*,pBemKernel*,pBemKernel >(&set_eqdestroy_incr));
    TheOperators->Add("*",new OneBinaryOperator_st< Op_coeffBemKernel1<pBemKernel,Complex,pBemKernel> >);
    
    Dcl_Type< const CBemDomainOfIntegration * >( );
    Dcl_Type< const CPartBemDI * >( );
    
    Add< const CPartBemDI * >("(", "", new OneOperatorCode< CBemDomainOfIntegration >);
    
    Add< const CBemDomainOfIntegration * >("(", "", new OneOperatorCode< BemKFormBilinear >);
    Add< const CDomainOfIntegration * >("(", "", new OneOperatorCode< BemPFormBilinear >);
        
    Global.Add("BEM","(",new FormalKBEMcode);
    Global.Add("POT","(",new FormalPBEMcode);
 //   Global.Add("Kernel","(",new FormalBemKernel);
    Add< pBemKernel >("<--","(",new FormalBemKernel);
//    Global.Add("Potential","(",new FormalBemPotential);
    Add< pBemPotential >("<--","(",new FormalBemPotential);

    Global.Add("int2dx2d","(",new OneOperatorCode<CPartBemDI2d2d>);
    Global.Add("int1dx1d","(",new OneOperatorCode<CPartBemDI1d1d>);
    Global.Add("int1dx2d","(",new OneOperatorCode<CPartBemDI1d2d>);
    Global.Add("int2dx1d","(",new OneOperatorCode<CPartBemDI2d1d>);

    Global.New("htoolEta",CPValue<double>(ff_htoolEta));
    Global.New("htoolEpsilon",CPValue<double>(ff_htoolEpsilon));
    Global.New("htoolMinclustersize",CPValue<long>(ff_htoolMinclustersize));
    Global.New("htoolMaxblocksize",CPValue<long>(ff_htoolMaxblocksize));
    Global.New("htoolMintargetdepth",CPValue<long>(ff_htoolMintargetdepth));
    Global.New("htoolMinsourcedepth",CPValue<long>(ff_htoolMinsourcedepth));
    ArrayofHmat<double>();
    ArrayofHmat<complex<double>>();

    // Build HMatrix from custom user generator
    Dcl_TypeandPtr< VirtualGenerator<double>* >(0, 0,::InitializePtr< VirtualGenerator<double>* >, ::DeletePtr< VirtualGenerator<double>*> );
    Dcl_TypeandPtr< VirtualGenerator<std::complex<double>>* >(0, 0,::InitializePtr< VirtualGenerator<std::complex<double>>* >, ::DeletePtr< VirtualGenerator<std::complex<double>>*> );
    zzzfff->Add("Generator", atype<VirtualGenerator<double>** >( ));
    map_type_of_map[make_pair(atype<VirtualGenerator<double>**>(), atype<double*>())] = atype<VirtualGenerator<double>**>();
    map_type_of_map[make_pair(atype<VirtualGenerator<double>**>(), atype<Complex*>())] = atype<VirtualGenerator<std::complex<double> >**>();

    AddHMatrixUser<double,v_fesL, v_fesL>();
    AddHMatrixUser<double,v_fesS, v_fesS>();
    AddHMatrixUser<double,v_fes3, v_fes3>();
    AddHMatrixUser<std::complex<double>,v_fesL, v_fesL>();
    AddHMatrixUser<std::complex<double>,v_fesS, v_fesS>();
    AddHMatrixUser<std::complex<double>,v_fes3, v_fes3>();
    AddHMatrixUser<double,v_fesL, v_fesS>();
    AddHMatrixUser<double,v_fesS, v_fes3>();
    AddHMatrixUser<std::complex<double>,v_fesL, v_fesS>();
    AddHMatrixUser<std::complex<double>,v_fesS, v_fes3>();
}

LOADFUNC(Init_Bem)

