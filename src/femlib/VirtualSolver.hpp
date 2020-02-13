#ifndef __VirtualSolver_HPP__
#define __VirtualSolver_HPP__
#include <map>
#include <vector>
#include <string>
extern double ff_tgv;
#define VDATASPARSESOLVER  1

/*struct TypeSolveMat {
    enum TSolveMat { NONESQUARE=0, LU=1, CROUT=2, CHOLESKY=3, GC = 4 , GMRES = 5, SparseSolver=6, SparseSolverSym=7 };
    TSolveMat t;
    bool sym;
    bool profile;
    TypeSolveMat(TSolveMat tt=LU) :t(tt),
    sym(t == CROUT || t ==CHOLESKY  ||  t==GC || t==SparseSolverSym ),
    profile(t == CROUT || t ==CHOLESKY  || t ==LU ) {}
    bool operator==(const TypeSolveMat & a) const { return t == a.t;}
    bool operator!=(const TypeSolveMat & a) const { return t != a.t;}
    static TSolveMat defaultvalue;
};
*/
// add FH , JM  avril 2009
template<class K,class V> class MyMap;
class String;
typedef void *    pcommworld; // to get the pointeur to the comm word ... in mpi

int Data_Sparse_Solver_version() ; //{ return VDATASPARSESOLVER;}


struct Data_Sparse_Solver {
    
    static std::map<std::string,int> *  mds;
    typedef std::map<std::string,int>::const_iterator IMDS;

    Expression einitmat;
    bool initmat;
    string  solver;
    double epsilon;
    const void * precon;
    int NbSpace;
    int strategy;
    double tgv;
    int factorize;
    double tol_pivot;
    double tol_pivot_sym;
    int itmax ;
    string data_filename;
    KN<long> lparams;  //  copy arry more secure ...
    KN<double> dparams;
    
    MyMap<String,String> * smap;
    
    KN<long> perm_r;
    KN<long> perm_c;
    KN<double> scale_r;
    KN<double> scale_c;
    string sparams;
    pcommworld commworld;  // pointeur sur le commworld
    int master; //  master rank in comm add FH 02/2013 for MUMPS ... => VDATASPARSESOLVER exist
    // array for return information for mumps ...
    KN<double> * rinfo;
    KN<long> * info;
    
    KNM<double>* kerneln;
    KNM<double> * kernelt;
    long *kerneldim;
    long  verb;
    bool x0; //  init by 0 the inital data the solution
    double * veps; //    to get and set value of eps
    bool rightprecon;
    bool sym;// know is set or not
    bool positive;
    
    Data_Sparse_Solver()
    :
    einitmat(0),
    solver(""),
    epsilon(1e-6),
    precon(0),
    NbSpace(1000),
    strategy(0),
    tgv(ff_tgv),
    factorize(0),
    tol_pivot(-1),
    tol_pivot_sym(-1),
    itmax(0),
     smap(0) ,
    commworld(0),
    master(0),
    rinfo(0),
    info(0),
    kerneln(0), kernelt(0), kerneldim(0),verb(verbosity) ,x0(true),veps(0),
    rightprecon(true),
    sym(false),
    positive(false)
    {}
    void Update(Stack s)
    {
        if(einitmat)
            initmat= GetAny<bool>((*einitmat)(s));
    }
    template<class R> 
    void Init_sym_positive_var(int syma);

    
    static   std::map<std::string,int> * Set_mds()
    {
        using namespace std;
        static int init=0;
        static map<string,int> ms;
        
        if( init==0)
        {
            ms["eps"]=1;
            ms["precon"]=3;
            ms["nbkrilov"]=4;
            ms["nbspace"]=4;
            ms["strategy"]=5;
            ms["tgv"]=6;
            ms["tol_pivot"]=8;
            ms["tol_pivot_sym"]=9;
            ms["itmax"]=10;
            ms["data_filename"]=11;
            ms["verbosity"]=12;
        }
        return & ms;
    }
    void  set(va_list ap)
    {
        int k=0;
        const char *what;
        
        while (( what=va_arg(ap,const  char *) ))
        {
            k++;
            assert(k<20);
            if( mds==0) mds = Set_mds();
            if( what==0 ) break;
            IMDS iw = mds->find(what);
            if(iw == mds->end()) {break;}
            int cas = iw->second;
            switch (cas)
            {
                case 1: epsilon=va_arg(ap,double);
                    cout << " ds : eps = " << epsilon << endl; break;
                case 3: precon=va_arg(ap,void *);
                    cout << " ds : precon = " << precon << endl; break;
                case 4: NbSpace=va_arg(ap,int);
                    cout << " ds : nbspace (krilov)  = " << NbSpace << endl; break;
                case 5: strategy=va_arg(ap,int);
                    cout << " ds : strategy  = " << strategy << endl; break;
                case 6: tgv=va_arg(ap,double);
                     cout << " ds : tgv  = " << tgv  << endl; break;
                case 8: tol_pivot=va_arg(ap,double);  cout << " ds : tol_pivot  = " << tol_pivot << endl; break;
                case 9: tol_pivot_sym=va_arg(ap,double); cout << " ds : tol_pivot_sym  = " << tol_pivot_sym << endl;  break;
                case 10: itmax=va_arg(ap,int);
                    cout << " ds : itmax  = " << itmax << endl; break;
                case 12: verbosity=va_arg(ap,int);
                    cout << " ds : itmax  = " << itmax << endl; break;
                    // case 11: data_filename=va_arg(ap,string); break;
            }
        }
    }
private:
    Data_Sparse_Solver(const Data_Sparse_Solver& ); // pas de copie
    
};


//template<class K>  using VirtualSolverN=typename VirtualMatrix<int,K>::VSolver;
//typedef VirtualSolverN<double>  * (*SparseRMatSolve)(DCL_ARG_SPARSE_SOLVER(double,A) );
//typedef VirtualSolverN<Complex> * (*SparseCMatSolve)(DCL_ARG_SPARSE_SOLVER(Complex,A) );

/*
template<class R,int sympos> struct DefSparseSolverNew {
    typedef VirtualSolverN<R>  *
    (*SparseMatSolver)(DCL_ARG_SPARSE_SOLVER(R,A) );
    static SparseMatSolver solver;
    static SparseMatSolver solverdef;
    static   VirtualSolverN<R> * Build( DCL_ARG_SPARSE_SOLVER(R,A) )
    {
        VirtualSolverN<R> *ret=0;
        if(solver)
            ret =(solver)(ARG_SPARSE_SOLVER(A));
        return ret;
    }
};
*/

template<class I, class R>
class VirtualSolver : public VirtualMatrix<I,R>::VSolver  {
public:
    typedef I INDEX;
    typedef R SCALAR;
    
    int state,defMatType,MatType;
    long codeini,codesym,codenum;
    long codeininew,codesymnew,codenumnew;// previous code
    VirtualSolver(int dmt=0) : state(0),codeini(0),codesym(0),codenum(0),defMatType(dmt),MatType(dmt) {}
 
    virtual void dosolver(R *x,R*b,int N=0,int trans=0) =0;
    
    virtual void fac_init(){}  // n, nzz fixe
    virtual void fac_symbolic(){} //  i,j fixe
    virtual void fac_numeric(){}   // a fixe
    virtual void UpdateState(){}// to see the code change to redo fac ????

    //  MatType :  & 2  symetric or not
    //             & 4  positive of not
    //             & 8  mpi or not
    //             & 16  operator or not ????
    
    
    virtual int GetDefautMatType() { return defMatType;}// 0 LU, 1 sym , 4 positive  8 MPI
    virtual void SetDefautMatType(int MMatType ) { ffassert(MMatType==defMatType);  MatType=MMatType;}

    void ChangeCodeState(long ci=0,long cs=0, long cn=0)
    {
        if(ci) codeininew=ci;
        if(cs) codesymnew=ci;
        if(cn) codenumnew=ci;
        if(codeininew != codeini) state=0; // redo init
        else if(codesymnew != codesym) state=1; // redo init
        else if(codenumnew != codenumnew) state=2; // redo init
    };
    void factorize(int st)
    {
        UpdateState();
        if(verbosity>9) cout << " VirtualSolver :: factorize state:" << state << " st= "<< st << endl;

        if( (state==0) && (state <st)) {codeini=codeininew;fac_init(); state=1;}
        if( (state==1) && (state <st)) {codesym=codesymnew;fac_symbolic(); state=2;}
        if( (state==2) && (state <st)) {codenum=codenumnew;fac_numeric();state=3;}
    }
    R* solve(R *x,R *b,int N=1,int trans=0)
    {
        factorize(3);
        dosolver(x,b,N,trans);
        return x;
    }
    
    virtual ~VirtualSolver(){}
};

#endif
