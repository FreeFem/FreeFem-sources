#include "ff++.hpp"


static void Load_Init( ) {    // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
  // if (verbosity)
 if(mpirank==0)
    cout << " load: msh3 is obsolete (in kernel of freefem 4/nov/2024 FH"  << endl;
}

LOADFUNC(Load_Init)
