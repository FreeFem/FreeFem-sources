#include "SparseLinearSolver.hpp"
#include <complex>
template class SparseLinearSolver<int,double>;
template class SparseLinearSolver<int,std::complex<double> >;

std::map<std::string,int> * Data_Sparse_Solver::mds = Data_Sparse_Solver::Set_mds();
