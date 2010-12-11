// empty parallele interface if no MPI to build dll
// with or without mpi..

extern void (*initparallele)(int &argc, char **& argv) ;
extern void (*init_lgparallele)();
extern void (*end_parallele)();
extern void Set_pparallele();

