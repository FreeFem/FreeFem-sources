int mainff (int  argc, char **argv);

void (*initparallele)(int &argc, char **& argv)=0 ;
void (*init_lgparallele)()=0;
void (*end_parallele)()=0;

extern void init_ptr_parallelepmi();

int mymain (int  argc, char **argv)
{

 init_ptr_parallelepmi();
 return mainff(argc,argv);
}
