/// \file

int mainff (int  argc, char **argv);
namespace ffapi { void init();}
extern void init_ptr_parallelepmi();

/// called by platform-dependent main() in [[file:../Graphics/sansrgraph.cpp::calling_mymain]]
int mymain (int  argc, char **argv)
{
  ffapi::init(); // [[file:~/ff/src/fflib/ffapi.cpp::init]]

  // Calls either [[file:~/ff/src/mpi/parallelempi.cpp::init_ptr_parallelepmi]] or
  // [[file:~/ff/src/mpi/parallelempi-empty.cpp::init_ptr_parallelepmi]]
  init_ptr_parallelepmi();

  return mainff(argc,argv); // [[file:lg.ypp::mainff]]
}
