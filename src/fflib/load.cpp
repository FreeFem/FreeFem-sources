#include "config-wrapper.h" // needed for HAVE_DLFCN_H

#include  <iostream>
#include  <map>
#include "AFunction.hpp"
using namespace std;
#include "lex.hpp"
#define LOAD 1
#if defined(__INTEL__) || defined(__MWERKS__) || !defined(HAVE_DLFCN_H)
#undef LOAD
#endif

#ifdef LOAD
#include <dlfcn.h>
#elif WIN32
#include <windows.h>
#endif
bool load(string s)
{
  bool ret=false;
  void * handle = 0;
#ifdef LOAD  
  handle = dlopen (s.c_str(), RTLD_LAZY ); 
  ret= handle !=0;
  if ( handle == 0 ) 
    cerr  <<   "load: " <<s << "\n \t fail : " << dlerror() << endl;
  else 
    cout << "load: dlopen(" << s << ") = " << handle << endl;
#elif WIN32
  {
   HINSTANCE mod=  LoadLibrary(s.c_str());
   if(mod==0) 
     {
       DWORD merr = GetLastError();
    cerr  <<   "loadLibary : " <<s << "\n \t fail : " << merr << endl;
     }
  else 
    cout << "load: LoadLibrary (" << s << ") = " << mod << endl;

  }
#else
  cout << "load: sorry no dlopen on this system " << s << " \n" ;
#endif  
 
  return handle ;
}

