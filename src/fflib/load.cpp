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
bool load(string ss)
{
  bool ret=false;
  void * handle = 0;
  const int nbprefix=2,nbsuffix=2;
  string  prefix[nbprefix];
  string suffix[nbsuffix] ;
  
  prefix[0]="";
  prefix[1]="./";
  
  suffix[0]="";
  suffix[1]=".so";
#ifdef  _APPLE_
  suffix[1]=".dylib";
#endif  
#ifdef WIN32  
  suffix[1]=".dll";
#endif 
 int i,j; 
 for ( i= 0;i< nbprefix;i++)
 for ( j= 0;j< nbsuffix;j++)
  {
    string s= prefix[i]+ss+suffix[i];
#ifdef LOAD  
  handle = dlopen (s.c_str(), RTLD_LAZY ); 
  ret= handle !=0;
#elif WIN32
  {
   HINSTANCE mod=  LoadLibrary(s.c_str());
   if(mod==0) 
     {
       DWORD merr = GetLastError();
       cerr  <<   "\n try loadLibary : " <<s << "\n \t fail : " << merr << endl;
     }
  else 
    cout << "\nload: LoadLibrary (" << s << ") = " << mod << endl;
   
  }
#else
  cout << "------------------------------------   \n" ;
  cout << "  load: sorry no dlopen on this system " << s << " \n" ;
  cout << "------------------------------------   \n" ;
  
#endif  
 }
  if ( ! ret ) 
    cerr  <<   "load error : [" << prefix[1]<<"]" <<ss <<"[" << suffix[1]<<"]"<< "\n \t fail : "  << endl;
  else 
    cout << "\nload: dlopen(" <<prefix[i] << ss << suffix[j] << ") = " << handle << endl;
  
  return handle ;
}

