// empty parallele interface if no MPI to build dll
// with or without mpi..
//#include "parallelepmi.hpp" 
extern  void (*initparallele)(int &argc, char **& argv) ;
extern void (*init_lgparallele)();
extern void (*end_parallele)();

void init_ptr_parallelepmi();
void init_ptr_parallelepmi(){};


