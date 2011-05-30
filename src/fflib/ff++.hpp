#ifndef FF___HPP_
#define FF___HPP_
// a not to bad list of freefem++ include 
// to simplify like of programmer.
// FH. 
#include <fstream>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <complex>
#include<stdlib.h>
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

#include "MatriceCreuse_tpl.hpp"
#include "MeshPoint.hpp"
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "Operator.hpp" 
#include "lex.hpp"
#include "libmesh5.h"
#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"

#include "Mesh2.h"
#include "BamgFreeFem.hpp"

#endif
