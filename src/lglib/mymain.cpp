int mainff (int  argc, char **argv);


extern void init_ptr_parallelepmi();

int mymain (int  argc, char **argv)
{

 init_ptr_parallelepmi();
 return mainff(argc,argv);
}
