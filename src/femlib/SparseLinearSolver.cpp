#include "SparseLinearSolver.hpp"
#include <complex>

template<class I,class R> typename TheFFSolver<I,R>::MAPSF TheFFSolver<I,R>::ffsolver;

template<class R>
void Data_Sparse_Solver::Init_sym_positive_var()
{
    //  put the solver name in UPPER CASE
   std::transform(solver.begin(), solver.end(), solver.begin(), static_cast<int(*)(int)>(std::toupper));
    if( solver == "CHOLESKY") {sym = true; positive = true;}
    if( solver == "CROUT") {sym = true;}
    if( solver == "CG") {sym = true;positive=true;}
    if( solver == "SPARSESOLVERSYM") {sym=true;}
    if( solver == "CHOLMOD") {sym=true;}
    
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
    string sn = ds.solver;
    std::transform(sn.begin(), sn.end(), sn.begin(), static_cast<int(*)(int)>(std::toupper));

    int typesolve = (sym ? VM::TS_sym: VM::TS_unsym ) + (pos ? VM::TS_def_positif : VM::TS_not_def_positif );
    auto i=  ffsolver.find(sn);
    int pp=-1;
    if( i == ffsolver.end()) // choose the best ???
        for ( auto j  = ffsolver.begin() ; j != ffsolver.end() ; ++j)
        {
            int ts=i->second->orTypeSol;
            int p=i->second->p;
            if( (ts  & typesolve) == typesolve ) // compiatile
            {
                
                if( pp < p)
                {
                    i=j; // the last   ???
                    pp=p;
                    if(verbosity>999)
                        cout << " find solver " << i->first << " "<< p << " / " << typesolve << " in "<< ts <<  endl;
                }
                
            }
            else
                if(verbosity>999)
                    cout << " not find solver " << i->first << " "<< p << " / " << typesolve << " in "<< ts <<  endl;
            
            
        }
    
    if( i != ffsolver.end())
        return i->second->create(A,ds,stack);
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







// explicit instentition of solver ...
std::map<std::string,int> * Data_Sparse_Solver::mds = Data_Sparse_Solver::Set_mds();

void init_SparseLinearSolver()
{
    InitSolver<int,double>();
    InitSolver<int,std::complex<double> >();
}

template class SparseLinearSolver<int,double>;
template class SparseLinearSolver<int,std::complex<double> >;

template class TheFFSolver<int,double>;
template class TheFFSolver<int,std::complex<double> >;

template int TypeOfMat<int,double>( Data_Sparse_Solver & ds);
template  int TypeOfMat<int,std::complex<double> >( Data_Sparse_Solver & ds);

template void Data_Sparse_Solver::Init_sym_positive_var<double>();
template void Data_Sparse_Solver::Init_sym_positive_var<complex<double> >();
