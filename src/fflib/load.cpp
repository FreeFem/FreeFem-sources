#include "config-wrapper.h" // needed for HAVE_DLFCN_H

#include  <iostream>
#include  <map>
#include "AFunction.hpp"
#include "environment.hpp"
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
  list<string> prefix(ffenvironment["load"]);
  if(prefix.empty())
    {
      prefix.push_back("");
      prefix.push_back("./");
    }

  string suffix[nbsuffix] ;
  
  suffix[0]="";
  suffix[1]=".so";
#ifdef  __APPLE__
  suffix[1]=".dylib";
#endif  
#ifdef WIN32  
  suffix[1]=".dll";
#endif 
  int j; 
  for (list<string>::const_iterator i= prefix.begin();i !=prefix.end();++i)
    for ( j= 0;j< nbsuffix;++j)
      {
	string s= *i+ss+suffix[j];
	
#ifdef LOAD  
	handle = dlopen (s.c_str(), RTLD_LAZY ); 
	if (verbosity>9) cout << " test dlopen(" << s << ")= " << handle <<  endl;
	ret= handle !=0;
	if (  ret ) 
	  {
	    if(verbosity)
	      cout << "\nload: dlopen(" << s << ") = " << handle << endl;
	    return handle;
	  }
	
#elif WIN32
	{
	  HINSTANCE mod=  LoadLibrary(s.c_str());
	  if (verbosity>9) cout << " test LoadLibrary(" << s << ")= " << mod <<  endl;
	  if(mod==0) 
	    {
	      DWORD merr = GetLastError();
	      if(verbosity>19)
		cerr  <<   "\n try loadLibary : " <<s << "\n \t fail : " << merr << endl;
	    }
	  else 
	    {
	      if(verbosity)
		cout << "\nload: loadLibary(" <<  s <<  ") = " << handle << endl;
	      return mod;
	    }
	}
#else
	cout << "------------------------------------   \n" ;
	cout << "  load: sorry no dlopen on this system " << s << " \n" ;
	cout << "------------------------------------   \n" ;
	CompileError("Error load");
	return 0;
#endif  
      }
  cerr  <<   "\nload error : " << ss << "\n \t fail : "  << endl;
  cerr << "list  prefix: " ;
  for (list<string>::const_iterator i= prefix.begin();i !=prefix.end();++i)
    cerr <<"'"<<*i<<"' ";
  cerr << "list  suffix : '"<< suffix[0] << "' , '"  << suffix[1] << "' "; 

  cerr << endl;
  CompileError("Error load");
  
  return 0 ;
}

