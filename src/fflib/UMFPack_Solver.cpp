//  file to add UMFPACK solver with dynamic load.

#include  <iostream>
using namespace std;

#include "rgraph.hpp"
#include "error.hpp"
#include "AFunction.hpp"
#include "MatriceCreuse_tpl.hpp"
#include "lgsolver.hpp"

#ifdef HAVE_LIBUMFPACK

void init_UMFPack_solver() {
    setptrstring(def_solver,"UMFPACK");
    setptrstring(def_solver_sym,"CHOLMOD");
    setptrstring(def_solver_sym_dp,"CHOLMOD");
   Global.New("HaveUMFPACK",CConstant<bool>(true));
 }
        
#else
void init_UMFPack_solver() {
    setptrstring(def_solver,"LU");
    setptrstring(def_solver_sym,"CROUT");
    setptrstring(def_solver_sym_dp,"CHOLESKY");

   if(verbosity&& (mpirank==0))
        cout << " no UMFPACK -> replace by LU or GMRES  ";
 
    Global.New("HaveUMFPACK",CConstant<bool>(false));
}

#endif
    
