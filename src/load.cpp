#include  <iostream>
#include  <map>
#include "AFunction.hpp"
using namespace std;
#include "lex.hpp"
#define LOAD 1
#if defined(__INTEL__) || defined(__MWERKS__)
#undef LOAD
#endif

#ifdef LOAD
#include <dlfcn.h>
#endif

bool load(string s)
{

  void * handle = 0;
#ifdef LOAD  
  handle = dlopen (s.c_str(), RTLD_LAZY ); 
  if ( handle == 0 ) 
    cerr  <<   "load: " <<s << "\n \t fail : " << dlerror() << endl;
  else 
    cout << "load: dlopen(" << s << ") = " << handle << endl;
#else
  cout << "load: sorry no dlopen on this system " << s << " \n" ;
#endif  
 
  return handle ;
}

