#include "SparseLinearSolver.hpp"
#include <complex>

template<class I,class R> typename TheFFSolver<I,R>::MAPSF TheFFSolver<I,R>::ffsolver;

 void setptrstring( string * & ds, const string & s)
{
    if(ds) delete ds;
    ds = new string(s);
}
template<class R>
void Data_Sparse_Solver::Init_sym_positive_var(int syma)
{
    //  put the solver name in UPPER CASE
   std::transform(solver.begin(), solver.end(), solver.begin(), static_cast<int(*)(int)>(std::toupper));
    auto i=  TheFFSolver<int,R>::ffsolver.find(solver);
    if ( i != TheFFSolver<int,R>::ffsolver.end())
    {
        //  1 unsym , 2 sym, 4 pos , 8 nopos, 16  seq, 32  ompi, 64 mpi ,
        int ts = i->second->orTypeSol ;
        if ( (ts & 3) == 3 ) sym= syma <0 ? 0 : syma ; //  both sym/unsym -> prev. value
        else sym =  (ts & 2) == 2 ;
        positive = (ts & (4+8)) == 4;
        if(verbosity>4)
            cout <<  "  The solver "<< solver << " need sym "<< sym << " and  positif def "
                  << positive << " matrix ( prev sym"<< syma <<" ts " << ts << " )  \n";
    }
    else
    {
        if( solver == "CHOLESKY") {sym = true; positive = true;}
        if( solver == "CROUT") {sym = true;}
        if( solver == "CG") {sym = true;positive=true;}
        if( solver == "SPARSESOLVERSYM") {sym=true;}
        if( solver == "CHOLMOD") {sym=true;}
    }
}

template<class Z,class K>
  int TypeOfMat( Data_Sparse_Solver & ds)
{
    string sn = ds.solver;
    typedef  VirtualMatrix<Z,K> VM;
    auto i=  TheFFSolver<Z,K>::ffsolver.find(sn);
    int sym=ds.sym, pos = ds.positive;
    if( sn == "CHOLESKY") {sym = true; pos = true;}
    if( sn == "CROUT") {sym = true;}
    if( sn == "CG") {sym = true;pos=true;}
    if( sn == "SPARSESOLVERSYM") {sym=true;}
    return sym + pos*2;
}

template<class Z,class K>
 typename VirtualMatrix<Z,K>::VSolver * TheFFSolver<Z,K>::Find(HashMatrix<Z,K> &A, const Data_Sparse_Solver & ds,Stack stack )
{
    //  1 unsym , 2 sym, 4 pos , 8 nopos, 16  seq, 32  ompi, 64 mpi ,
    //    static const int  TS_unsym=1, TS_sym=2, TS_def_positif=4,  TS_not_def_positif=8, TS_sequental = 16, TS_mpi = 32;
    typedef  VirtualMatrix<Z,K> VM;
    int sym=ds.sym, pos = ds.positive;
    if(verbosity>3) cout << " ** Search solver "<< ds.solver << " sym = " << sym << " pos.  " << pos << " half "<< A.half <<endl;

    string sn = ds.solver;
    std::transform(sn.begin(), sn.end(), sn.begin(), static_cast<int(*)(int)>(std::toupper));

    int typesolve = (sym ? VM::TS_sym: VM::TS_unsym ) + (pos ? VM::TS_def_positif : VM::TS_not_def_positif );
    auto i=  ffsolver.find(sn);
    auto ii=i;
    int pp=-1;
    if( i == ffsolver.end()) // choose the best ???
        for ( auto j  = ffsolver.begin() ; j != ffsolver.end() ; ++j)
        {
            int ts=j->second->orTypeSol;
            int p=j->second->p;
            if( (ts  & typesolve) == typesolve ) // compiatile
            {

                if( pp < p)
                {
                    i=j; // the last   ???
                    pp=p;
                    if(verbosity>9)
                        cout << " find solver " << i->first << " "<< p << " / " << typesolve << " in "<< ts <<  endl;
                }

            }
            else
                if(verbosity>9)
                    cout << " not find solver " << i->first << " "<< p << " / " << typesolve << " in "<< ts <<  endl;


        }

    if( i != ffsolver.end())
    {
        if(verbosity>2) cout << " ** Find solver "<< i->first << " ts: "<< i->second->orTypeSol
                             << " sym = " << sym << " pos.  " << pos << " half "<< A.half <<endl;

        return i->second->create(A,ds,stack);
    }
    else
    {
        cerr << " TheFFSolver<Z,K>::FATAL ERROR  impossible to find solver  \n"
        << " want "<< ds.solver << " type " << typesolve << endl;
        cerr << " dispo \n";
          for ( auto j  = ffsolver.begin() ; j != ffsolver.end() ; ++j)
          {
              cerr << j->first << "  orTypeSol: " << j->second->orTypeSol <<endl;
          }
        ExecError(" No Solver ????");
        return 0;
    }

}


template<class R>
void SetSolver(Stack stack,bool VF,VirtualMatrix<int,R> & A,const  Data_Sparse_Solver & ds)
{
    using namespace Fem2D;
    typename  VirtualMatrix<int,R>::VSolver * solver=0;
    HashMatrix<int,R> * AH(dynamic_cast<HashMatrix<int,R> *>(&A));
    ffassert(AH);
    solver = NewVSolver<int,R>(*AH,ds,stack);
    if(solver)
    {
        A.SetSolver(solver,true);
        if(ds.factorize)
        {
            solver->factorize(ds.factorize);// full factorization  mars 2019 FH/PHT
        }
    }
    else
        CompileError("SetSolver: type resolution unknown");

}


template<class R>
void DefSolver(Stack stack, VirtualMatrix<int,R>  & A,const Data_Sparse_Solver & ds)
{
    typename  VirtualMatrix<int,R>::VSolver * solver=0;
    HashMatrix<int,R>* AH(dynamic_cast<HashMatrix<int,R> *>(&A));
    ffassert(AH);

    solver = NewVSolver<int,R>(*AH,ds,stack);

    if(solver)
        A.SetSolver(solver,true);
    else
        CompileError("SetSolver: type resolution unknown");


}

typedef double R;
typedef complex<double> C;

// explicit instentition of solver ...
std::map<std::string,int> * Data_Sparse_Solver::mds = Data_Sparse_Solver::Set_mds();

void init_SparseLinearSolver()
{
    InitSolver<int,R>();
    InitSolver<int,C>();
}

template class SparseLinearSolver<int,R>;
template class SparseLinearSolver<int,C>;

template class TheFFSolver<int,R>;
template class TheFFSolver<int,C>;

template int TypeOfMat<int,R>( Data_Sparse_Solver & ds);
template  int TypeOfMat<int,C>( Data_Sparse_Solver & ds);

template void Data_Sparse_Solver::Init_sym_positive_var<R>(int );
template void Data_Sparse_Solver::Init_sym_positive_var<C>(int );

template void SetSolver(Stack stack,bool VF,VirtualMatrix<int,R> & A, const Data_Sparse_Solver & ds);
template void SetSolver(Stack stack,bool VF,VirtualMatrix<int,C> & A, const     Data_Sparse_Solver & ds);

template void DefSolver(Stack stack, VirtualMatrix<int,R>  & A,const Data_Sparse_Solver & ds);
template void DefSolver(Stack stack, VirtualMatrix<int,C> & A,const  Data_Sparse_Solver & ds);
