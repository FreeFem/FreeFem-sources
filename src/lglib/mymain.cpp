/// \file


int mainff (int  argc, char **argv);
namespace ffapi { void init();}
extern void init_ptr_parallelepmi();

/// called by platform-dependent main() in src/Graphics/sansrgraph.cpp and others.
int mymain (int  argc, char **argv)
{
 ffapi::init(); 
 init_ptr_parallelepmi();
 return mainff(argc,argv); // [[file:lg.ypp::mainff]]
}
