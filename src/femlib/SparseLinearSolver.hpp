#ifndef __SparseLinearSolver_HPP__
#define __SparseLinearSolver_HPP__
#include <cstdarg>
#include "VirtualSolverCG.hpp"
#include "VirtualSolverSparseSuite.hpp"
#include "VirtualSolverSkyLine.hpp"
#include <map>
#include <vector>
#include <cstring> 

template<class K,class V> class MyMap;
class String;
typedef void *    pcommworld; // to get the pointeur to the comm word ... in mpi
//  to build

template<class Z, class K>
typename VirtualMatrix<Z,K>::VSolver * NewVSolver(HashMatrix<Z,K> &A,const char *solver, Data_Sparse_Solver & ds )
{
    VirtualSolver<Z,K> *thesolver;
    if(strncmp("UMFPACK",solver,6)==0)
        thesolver = new VirtualSolverUMFPACK<Z,K> (A,ds.strategy,ds.tol_pivot,ds.tol_pivot_sym);
    else if(strncmp("CHOLMOD",solver,7)==0)
        thesolver = new VirtualSolverCHOLMOD<Z,K> (A);
    else if(strncmp("CG",solver,2)==0)
        thesolver = new SolverCG<Z,K> (A,ds);
    else if(strncmp("GMRES",solver,5)==0)
        thesolver = new SolverGMRES<Z,K> (A,ds);
    else if(strncmp("LU",solver,2)==0)
        thesolver=new VirtualSolverSkyLine<Z,K> (&A,3,ds.tol_pivot);
    else if(strncmp("CROUT",solver,5)==0)
        thesolver=new VirtualSolverSkyLine<Z,K> (&A,2,ds.tol_pivot);
    else if(strncmp("CHOLESKY",solver,8)==0)
        thesolver=new VirtualSolverSkyLine<Z,K> (&A,1,ds.tol_pivot);
    else if(strncmp("MUMPS",solver,5)==0)
        thesolver=0;
    else if(strncmp("SUPERLU",solver,6)==0)
        thesolver=0;
    if(thesolver ==0) { std::cerr << " Solver linear inconnue " << solver << " => UMFPACK " << std::endl;
        {
            thesolver = new VirtualSolverUMFPACK<Z,K> (A);
        }
    }
    
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
        // int strategy=-1,double tol_pivot=-1.,double tol_pivot_sym=-1.
        if(strncmp("UMFPACK",solver,6)==0)
         thesolver = new VirtualSolverUMFPACK<Z,K> (A,ds.strategy,ds.tol_pivot,ds.tol_pivot_sym);
        else if(strncmp("CHOLMOD",solver,7)==0)
         thesolver = new VirtualSolverCHOLMOD<Z,K> (A);

        else if(strncmp("CG",solver,2)==0)
          thesolver = new SolverCG<Z,K> (A,ds);
        else if(strncmp("GMRES",solver,5)==0)
            thesolver = new SolverGMRES<Z,K> (A,ds);
        else if(strncmp("LU",solver,2)==0)
            thesolver=new VirtualSolverSkyLine<Z,K> (&A,3,ds.tol_pivot);
        else if(strncmp("CROUT",solver,5)==0)
            thesolver=new VirtualSolverSkyLine<Z,K> (&A,2,ds.tol_pivot);
        else if(strncmp("CHOLESKY",solver,8)==0)
            thesolver=new VirtualSolverSkyLine<Z,K> (&A,1,ds.tol_pivot);
        else if(strncmp("MUMPS",solver,5)==0)
            thesolver=0;
        else if(strncmp("SUPERLU",solver,6)==0)
            thesolver=0;
        if(thesolver ==0) { std::cerr << " Solver linear inconnue " << solver << " => UMFPACK " << std::endl;
            {
                thesolver = new VirtualSolverUMFPACK<Z,K> (A);
            }
            
        }

    }
     void solve(K *u,K *b,int nrhs=1){thesolver->solve(u,b,nrhs);}
    ~SparseLinearSolver(){ delete thesolver;}
};

/*
template<class TypeIndex=int,class TypeScalar=double>
inline double * ProduitMatVec(const SparseLinearSolver<TypeIndex,TypeScalar> *A,TypeScalar *x, TypeScalar *Ax) { return VR->A->matmul(x,Ax);}
template<class TypeIndex=int,class TypeScalar=double>
inline double * ProduitMatVec(const SparseLinearSolver<TypeIndex,TypeScalar> &A,TypeScalar *x, TypeScalar *Ax) { return VR->A.matmul(x,Ax);}
*/
#endif
