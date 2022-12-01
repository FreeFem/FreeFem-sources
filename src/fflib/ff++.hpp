#ifndef FF___HPP_
#define FF___HPP_
#if defined(__clang__) && defined(__has_warning)
#if __has_warning("-Wundefined-var-template")
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wundefined-var-template"
#endif
#endif
// a not to bad list of freefem++ include 
// to simplify like of programmer.
// FH. 
#include <fstream>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <complex>
#include<cstdlib>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
#include "ufunction.hpp"
using namespace std;  
#include "rgraph.hpp"
#include "RNM.hpp"
// after RNM   otherwise 
// trouble with index in RNM (I do no understander FH)
#include <set>
#include <vector>
#include <map>

#include "fem.hpp"
#include "FESpacen.hpp" 
#include "FESpace.hpp" 
#include "HashMatrix.hpp"
#include "SparseLinearSolver.hpp"       
#include "MatriceCreuse_tpl.hpp"
#include "MatriceCreuse.hpp"
#include "MatriceCreuse_tpl.hpp"
#include "VirtualSolverCG.hpp"
//#include "VirtualSolverSparseSuite.hpp"
//#include "VirtualSolverSkyLine.hpp"
//#include "SparseLinearSolver.hpp"
#include "lgsolver.hpp"

#include "MeshPoint.hpp"
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "Operator.hpp" 
#include "lex.hpp"
#include "libmeshb7.h"
#include "lgfem.hpp"
#include "lgmesh3.hpp"

#include "problem.hpp"

#include "Mesh2.h"
#include "BamgFreeFem.hpp"
#include "ffapi.hpp" 
#if defined(__clang__) && defined(__has_warning)
#if __has_warning("-Wundefined-var-template")
#pragma clang diagnostic pop
#endif
#endif
#endif
