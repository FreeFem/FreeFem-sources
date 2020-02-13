/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-unviersite.fr

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:  mumps parmetis [ptscotch scotch]  scalapack blas  mpifc  fc mpi  pthread
//ff-c++-cpp-dep:
// *INDENT-ON* //

// F. Hecht  december 2011
// ----------------------------
// file to add MUMPS sequentiel interface for sparce linear solver with dynamic load.
#include <mpi.h>
#ifdef _WIN32
__declspec(dllexport) int toto;
MPI_Fint* _imp__MPI_F_STATUS_IGNORE;
MPI_Fint* _imp__MPI_F_STATUSES_IGNORE;
//__declspec(dllexport) void __guard_check_icall_fptr(unsigned long ptr) { }
#endif
#include  <iostream>
using namespace std;

#include "ff++.hpp"


#include <dmumps_c.h>
#include <zmumps_c.h>

const int JOB_INIT = -1;
const int JOB_END = -2;
const int JOB_ANA = 1;
const int JOB_FAC = 2;
const int JOB_ANA_FAC = 4;
const int JOB_SOLVE = 3;
const int USE_COMM_WORLD = -987654;

template<typename RR> struct MUMPS_STRUC_TRAIT {typedef void MUMPS;  typedef void R; };
template<> struct MUMPS_STRUC_TRAIT<double>  {typedef DMUMPS_STRUC_C MUMPS; typedef double R;};
template<> struct MUMPS_STRUC_TRAIT<Complex>  {typedef ZMUMPS_STRUC_C MUMPS; typedef ZMUMPS_COMPLEX R;};
void mumps_c(DMUMPS_STRUC_C *id) { dmumps_c(id);}
void mumps_c(ZMUMPS_STRUC_C *id) { zmumps_c(id);}

template<class T> struct MPI_TYPE {static MPI_Datatype  TYPE(){return MPI_BYTE;}};;
template<> struct MPI_TYPE<long>      {static MPI_Datatype  TYPE(){return MPI_LONG;}};
template<> struct MPI_TYPE<int>      {static MPI_Datatype TYPE(){return MPI_INT;}};
template<> struct MPI_TYPE<double>    {static MPI_Datatype TYPE(){return MPI_DOUBLE;}};
template<> struct MPI_TYPE<char>    {static MPI_Datatype TYPE(){return MPI_BYTE;}};
template<> struct MPI_TYPE<Complex>    {static MPI_Datatype TYPE(){return MPI_DOUBLE_COMPLEX;}};


static std::string analysis[] = {"AMD", "", "AMF", "SCOTCH", "PORD", "METIS", "QAMD", "automatic sequential", "automatic parallel", "PT-SCOTCH", "ParMetis"};

//template<typename R>
template<class R=double>
class SolveMUMPS_mpi: public  VirtualSolver<int,R>
{
public:
    //  1 unsym , 2 sym, 4 pos , 8 nopos, 16  seq, 32  ompi, 64 mpi ,
    static const int orTypeSol = 1|2|4|8|16;
    typedef HashMatrix<int,R>  HMat;
    typedef R K; //
    HMat &A;
    
    
    // typedef double R;
    long verb;
    double eps;
    double tgv;
    int cn,cs;
    typedef typename MUMPS_STRUC_TRAIT<R>::R MR;
    mutable typename MUMPS_STRUC_TRAIT<R>::MUMPS id;
    KN<double> *rinfog;
    KN<long> *infog;
    mutable unsigned char   strategy;
    bool distributed;
    MPI_Comm comm;
    int mpirank;
    int matrank;
    // int distributed;
    
    int&    ICNTL (int i) const {return id.icntl[i - 1];}
    double& CNTL  (int i) const {return id.cntl[i - 1];}
    int&    INFO  (int i) const {return id.info[i - 1];}
    double& RINFO (int i) const {return id.rinfo[i - 1];}
    int&    INFOG (int i) const {return id.infog[i - 1];}
    double& RINFOG (int i) const {return id.rinfog[i - 1];}
    
    void SetVerb () const {
        ICNTL(1) = 6;//   output stream for error messages.
        ICNTL(2) = 6;//  stream for diagnostic printing, statistics, and warning messages.
        ICNTL(3) = 6;//  output stream global information, collected on the host.
        ICNTL(4) = min(max(verb-2,1L),4L); // the level of printing for error, warning, and diag
        if(verb ==0 )ICNTL(4) =0;
        ICNTL(11)=0; // noerroranalysisisperformed(nostatistics).
        if( id.job ==JOB_SOLVE && verb >99)
        { //computes statistics related to an error analysis of the linear system
            if( verb > 999) ICNTL(11)=1; // All Stat (veryexpensive)
            else ICNTL(11)=2;// compute main statistics
        }
        
        
    }
    void Clean ()
    {
        delete [] id.irn;
        delete [] id.jcn;
        delete [] id.a;
        id.irn=0;
        id.jcn=0;
        id.a =0;
    }
    void to_mumps_mat()
    {
        Clean ();
        
        id.nrhs = 0;//
        int n = A.n;
        int nz = A.nnz;
        ffassert(A.n == A.m);
        if( distributed || (mpirank == matrank) )
        {
            int *irn = new int[nz];
            int *jcn = new int[nz];
            R *a = new R[nz];
            A.COO();
            
            for (int k = 0; k < nz; ++k) {
                {
                    irn[k] = A.i[k]+1;
                    jcn[k] = A.j[k] + 1;
                    a[k] = A.aij[k];
                }
            }
            
            id.n = n;
            
            if(!distributed)
            {
                if(mpirank == matrank)
                {
                  id.nz = nz;
                  id.irn = irn;
                  id.jcn = jcn;
                  id.a = (MR *)(void *)a;
                }
                else
                { //  no matrix
                    id.nz=0;
                    id.a =0;
                    id.irn = 0;;
                    id.jcn = 0;
                    
                }
                
            }
            else
            {
                id.nz_loc = nz;
                id.irn_loc = irn;
                id.jcn_loc = jcn;
                id.a_loc = (MR *)(void *)a;
                
            }
            id.rhs = 0;
            ffassert( A.half == (id.sym != 0) );//
            ICNTL(5) = 0;    // input matrix type
            ICNTL(7) = 7;    // NUMBERING ...
            
            ICNTL(9) = 1;    // 1: A x = b, !1 : tA x = b  during slove phase
            ICNTL(18) = 0;
            if(strategy > 0 && strategy < 9 && strategy != 2)
            {
                ICNTL(28) = 1;             // 1: sequential analysis
                ICNTL(7)  = strategy - 1; //     0: AMD
            }
            //     1:
            //     2: AMF
            //     3: SCOTCH
            //     4: PORD
            //     5: METIS
            //     6: QAMD
            //     7: automatic
            else
            {
                ICNTL(28) = 1;
                ICNTL(7)  = 7;
            }
            if(strategy > 8 && strategy < 12)
            {
                ICNTL(28) = 2;              // 2: parallel analysis
                ICNTL(29) = strategy - 9;  //     0: automatic
            }                                   //     1: PT-STOCH
            //     2: ParMetis
            ICNTL(9)  = 1;
            ICNTL(11) = 0;                 // verbose level
            ICNTL(18) = distributed ? 3: 0;        // centralized matrix input if !distributed
            ICNTL(20) = 0;                 // dense RHS
            ICNTL(14) = 30;
        }
    }
    void Check (const char *msg = "mumps_mpi")
    {
        if (INFO(1) != 0) {
            cout << " Erreur Mumps mpi: number " << INFO(1) << endl;
            cout << " Fatal Erreur  " << msg << endl;
            Clean ();
            id.job = JOB_END;
            mumps_c(&id);	/* Terminate instance */
            ErrorExec(msg, INFO(1));
        }
    }
    void CopyInfo()
    {
        if (rinfog) {
            // copy rinfog
            if (rinfog->N() < 40) {rinfog->resize(40);}
            
            for (int i = 0; i < 40; ++i) {
                (*rinfog)[i] = RINFOG(i + 1);
            }
        }
        
        if (infog) {
            // copy ginfo
            if (infog->N() < 40) {infog->resize(40);}
            
            for (int i = 0; i < 40; ++i) {
                (*infog)[i] = INFOG(i + 1);
            }
        }
    }
    SolveMUMPS_mpi (HMat  &AA, const Data_Sparse_Solver & ds,Stack stack )
    : A(AA), verb(ds.verb),
    eps(ds.epsilon),
    tgv(ds.tgv),cn(0),cs(0),
    rinfog(ds.rinfo), infog(ds.info),
    matrank(ds.master),distributed(ds.master<0),
    strategy(ds.strategy)
    {
        
        if(ds.commworld)
            MPI_Comm_dup(*((MPI_Comm*)ds.commworld), &comm);
	else
	    MPI_Comm_dup(MPI_COMM_WORLD, &comm);

        MPI_Comm_rank(comm, &mpirank);
        int master = mpirank==matrank;
        int myid = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        
        id.irn=0;
        id.jcn=0;
        id.a =0;
        
        id.job = JOB_INIT;
        id.par = 1;
        id.sym = A.half;
        id.comm_fortran = MPI_Comm_c2f(comm);
        SetVerb();
        mumps_c(&id);
        
        
        Check("MUMPS_mpi build/init");
        if (verbosity > 3 && master) {
            cout << "  -- MUMPS   n=  " << id.n << ", peak Mem: " << INFOG(22) << " Mb" << " sym: " << id.sym << endl;
        }
        
        
    }
    
    
    
    ~SolveMUMPS_mpi () {
        Clean ();
        id.job = JOB_END;
        SetVerb () ;
        mumps_c(&id);	/* Terminate instance */
        /*int ierr = */
        MPI_Comm_free(&comm);
    }
    
    
    void dosolver(K *x,K*b,int N,int trans)
    {
        size_t  nN=id.n*N;
        if (verbosity > 1 && mpirank==0) {
            cout << " -- MUMPS solve,  peak Mem : " << INFOG(22) << " Mb,   n = "
            << id.n << " sym =" << id.sym <<" trans = " << trans  << endl;
        }
        ICNTL(9) = trans == 0;    // 1: A x = b, !1 : tA x = b  during slove phase
        id.nrhs = N;
        // x = b;
        if(distributed)
        {
            MPI_Reduce( (void *) b,(void *)  x  , nN , MPI_TYPE<R>::TYPE(),MPI_SUM,0,comm);
        }
        else if(mpirank==0)  std::copy(b,b+nN,x);
        id.rhs = (MR *)(void *)(R *)x;
        id.job = JOB_SOLVE;    // performs the analysis. and performs the factorization.
        SetVerb();
        mumps_c(&id);
        Check("MUMPS_mpi dosolver");
        if(distributed) // send the solution ...
            MPI_Bcast(reinterpret_cast<void*> (x),nN, MPI_TYPE<R>::TYPE(), 0,comm);
        
        
        if (verb  > 9 && mpirank==0) {
            
            for(int j=0; j<N; ++j)
            {
                KN_<R> B(b+j*id.n,id.n);
                cout << j <<"   b linfty " << B.linfty()  << endl;
            }
        }
        
        if (verb > 2) {
            
            for(int j=0; j<N; ++j)
            {   KN_<R> B(x+j*id.n,id.n);
                cout << "   x  " << j <<"  linfty " << B.linfty() << endl;
            }
        }
        CopyInfo();
        
    }
    
    void fac_init(){
        to_mumps_mat();
    }  // n, nzz fixe
    void fac_symbolic(){
        id.job = JOB_ANA;
        SetVerb ();
        mumps_c(&id);
        Check("MUMPS_mpi Analyse");
        CopyInfo();
    }
    void fac_numeric(){
        id.job = JOB_FAC;
        SetVerb () ;
        mumps_c(&id);
        Check("MUMPS_mpi Factorize");
        CopyInfo();
    }
    void UpdateState(){
        if( A.GetReDoNumerics() ) cn++;
        if( A.GetReDoSymbolic() ) cs++;
        this->ChangeCodeState(A.n,cs,cn);
    }
    
};


static void Load_Init()
{
    
    addsolver<SolveMUMPS_mpi<double>>("MUMPS",50,1);
    addsolver<SolveMUMPS_mpi<Complex>>("MUMPS",50,1);
    addsolver<SolveMUMPS_mpi<double>>("MUMPSMPI",50,1);
    addsolver<SolveMUMPS_mpi<Complex>>("MUMPSMPI",50,1);
    setptrstring(def_solver,"MUMPSMPI");
}
LOADFUNC(Load_Init)
