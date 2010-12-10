// empty parallele interface if no MPI to build dll
// with or without mpi..

void initparallele(int &argc, char **& argv);
void init_lgparallele();
void end_parallele();


void initparallele(int &argc, char **& argv){}

void init_lgparallele() { }
void end_parallele(){}

