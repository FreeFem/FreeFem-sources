#ifndef __SparseLinearSolver_HPP__
#define __SparseLinearSolver_HPP__
#include <cstdarg>
#include "VirtualSolverCG.hpp"
#include "VirtualSolverSparseSuite.hpp"
#include "VirtualSolverSkyLine.hpp"
#include <map>
#include <vector>
#include <cstring>

#include "ffstack.hpp"

extern string *def_solver, *def_solver_sym, *def_solver_sym_dp;
inline void setptrstring( string * & ds, const string & s)
{
    if(ds) delete ds;
    ds = new string(s);
}
template<class K,class V> class MyMap;
class String;
typedef void *    pcommworld; // to get the pointeur to the comm word ... in mpi
//  to build
/*
template<class Z, class K>  struct BuidVSolver<Z,K> {
 virtual   VirtualMatrix<Z,K>::VSolver *build(HashMatrix<Z,K> &A,const char *solver,const  Data_Sparse_Solver & ds,Stack stack )=0;
};
 

template<class Z, class K> map<string, BuidVSolver<Z,K> > VSolverff ;// Choise  of New Solver
 */

template<class Z, class K>
struct TheFFSolver {
    
    struct OneFFSlver {
        int p;// preference ????`
        int orTypeSol; // type
        OneFFSlver(int pp,int tt):p(pp),orTypeSol(tt) {}
        virtual VirtualSolver<Z,K> * create(HashMatrix<Z,K> &A, const Data_Sparse_Solver & ds,Stack stack )=0;
        virtual ~OneFFSlver() {}
        
    };
    
    template<class VS>
    struct OneFFSlverVS:public OneFFSlver  {
        OneFFSlverVS(int pp) :OneFFSlver(pp,VS::orTypeSol) {}
        virtual VirtualSolver<Z,K> * create(HashMatrix<Z,K> &A, const Data_Sparse_Solver & ds,Stack stack )
        { return new VS(A,ds,stack);}
    };
    
    typedef pair<string,OneFFSlver *>  MValue;
    typedef  map<string,OneFFSlver *> MAPSF;
    static MAPSF ffsolver;

    template<class VS>
    static void addsolver (const char* name,int pp,int ts,const  VS* pvs)
    {
        string sn(name);
        std::transform(sn.begin(), sn.end(), sn.begin(), static_cast<int(*)(int)>(std::toupper));
        ffassert( ffsolver.find(sn) == ffsolver.end());
        MValue vm(sn,new OneFFSlverVS<VS>(pp));
        auto ii=ffsolver.insert(vm);
        ffassert( ii.second == true);
    }

    static typename VirtualMatrix<Z,K>::VSolver * Find(HashMatrix<Z,K> &A, const Data_Sparse_Solver & ds,Stack stack );
};



template<class TS> void addsolver(const char *nm,int p)
{
    typedef typename  TS::INDEX ZZ;
    typedef typename  TS::SCALAR KK;
    int ots=TS::orTypeSol;
    TheFFSolver<ZZ,KK>::addsolver(nm,p,ots, (TS*) 0); // trick (TS*) 0 to call the corre ct case ..
}

template<class Z, class K> void InitSolver()
{
   addsolver<VirtualSolverUMFPACK<Z,K>>("UMFPACK",100);
   addsolver<VirtualSolverCHOLMOD<Z,K>>("CHOLMOD",99);
   addsolver<SolverCG<Z,K>>("CG",10);
   addsolver<SolverGMRES<Z,K>>("GMRES",10);
   addsolver<VirtualSolverSkyLine<Z,K>>("LU",10);
   addsolver<VirtualSolverSkyLine<Z,K>>("CROUT",9);
   addsolver<VirtualSolverSkyLine<Z,K>>("CHOLESKY",9);

}
template<class Z, class K>
typename VirtualMatrix<Z,K>::VSolver * NewVSolver(HashMatrix<Z,K> &A, const Data_Sparse_Solver & ds,Stack stack )
{
    typename VirtualMatrix<Z,K>::VSolver *thesolver= TheFFSolver<Z,K>::Find(A,ds,stack);
    ffassert(thesolver);
    return thesolver;
    /*
     VirtualSolver<Z,K> *thesolver=0;
     if(strncmp("UMFPACK",solver,6)==0)
     thesolver = new VirtualSolverUMFPACK<Z,K> (A,solver,ds,stack); //.strategy,ds.tol_pivot,ds.tol_pivot_sym,ds.verb);
     else if(strncmp("CHOLMOD",solver,7)==0)
     thesolver = new VirtualSolverCHOLMOD<Z,K> (A,solver,ds,stack);
     else if(strncmp("CG",solver,2)==0)
        thesolver = new SolverCG<Z,K> (A,solver,ds,stack);
    else if(strncmp("GMRES",solver,5)==0)
        thesolver = new SolverGMRES<Z,K> (A,solver,ds,stack);
    else if(strncmp("LU",solver,2)==0)
        thesolver=new VirtualSolverSkyLine<Z,K> (&A,solver,ds,stack);
    else if(strncmp("CROUT",solver,5)==0)
        thesolver=new VirtualSolverSkyLine<Z,K> (&A,solver,ds,stack);
    else if(strncmp("CHOLESKY",solver,8)==0)
        thesolver=new VirtualSolverSkyLine<Z,K> (&A,solver,ds,stack);
    else if(strncmp("MUMPS",solver,5)==0)
        thesolver=0;
    else if(strncmp("SUPERLU",solver,6)==0)
        thesolver=0;
    if(thesolver ==0) { std::cerr << "\n ******* Solver linear inconnue " << solver << " => UMFPACK " << std::endl;
        {
            thesolver =new VirtualSolverUMFPACK<Z,K> (A,solver,ds,stack);
        }
    }
    */
    return thesolver;
}

template<class Z, class K>
class SparseLinearSolver  {
public:
    typedef HashMatrix<Z,K>  HMat;
    typedef VirtualSolver<Z,K> VR;
    VR *thesolver;
    SparseLinearSolver(HMat &A,const char *solver, ...)
    : thesolver(0)
    {
        /*
         SolverCG<int,R> AGC(A,eps,1,10);
         SolverGMRES<int,R> AGMRES(A,eps,1,10,std::max(50,int(sqrt(n)+10)));
         VirtualSolverUMFPACK<int,R> AUM(A);
         VirtualSolverCHOLMOD<int,R> ACH(A);

         */
        Data_Sparse_Solver ds;
        va_list val;
        va_start(val,solver);
        ds.set(val);
        va_end(val);
        ds.solver=solver;
        // int strategy=-1,double tol_pivot=-1.,double tol_pivot_sym=-1.
        if(strncmp("UMFPACK",solver,6)==0)
         thesolver = new VirtualSolverUMFPACK<Z,K> (A,ds,0);
        else if(strncmp("CHOLMOD",solver,7)==0)
         thesolver = new VirtualSolverCHOLMOD<Z,K> (A,ds,0);

        else if(strncmp("CG",solver,2)==0)
          thesolver = new SolverCG<Z,K> (A,ds,0);
        else if(strncmp("GMRES",solver,5)==0)
            thesolver = new SolverGMRES<Z,K> (A,ds,0);
        else if(strncmp("LU",solver,2)==0)
            thesolver=new VirtualSolverSkyLine<Z,K> (A,ds,0);
        else if(strncmp("CROUT",solver,5)==0)
            thesolver=new VirtualSolverSkyLine<Z,K> (A,ds,0);
        else if(strncmp("CHOLESKY",solver,8)==0)
            thesolver=new VirtualSolverSkyLine<Z,K> (A,ds,0);
        else if(strncmp("MUMPS",solver,5)==0)
            thesolver=0;
        else if(strncmp("SUPERLU",solver,6)==0)
            thesolver=0;
        if(thesolver ==0) { std::cerr << " Solver linear inconnue " << solver << " => UMFPACK " << std::endl;
            {
                thesolver = new  VirtualSolverUMFPACK<Z,K> (A,ds,0);
            }
            
        }

    }
     void solve(K *u,K *b,int nrhs=1){thesolver->solve(u,b,nrhs);}
    ~SparseLinearSolver(){ delete thesolver;}
};

void init_SparseLinearSolver();
/*
template<class TypeIndex=int,class TypeScalar=double>
inline double * ProduitMatVec(const SparseLinearSolver<TypeIndex,TypeScalar> *A,TypeScalar *x, TypeScalar *Ax) { return VR->A->matmul(x,Ax);}
template<class TypeIndex=int,class TypeScalar=double>
inline double * ProduitMatVec(const SparseLinearSolver<TypeIndex,TypeScalar> &A,TypeScalar *x, TypeScalar *Ax) { return VR->A.matmul(x,Ax);}
*/
#endif
