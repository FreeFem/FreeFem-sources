// -*- Mode : c++ -*-
//
// SUMMARY  :      
// USAGE    :        
// ORG      : 
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
//

/*
 
 This file is part of Freefem++
 
 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.
 
 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
#include "config-wrapper.h" // needed for HAVE_DLFCN_H

#include  <iostream>
#include  <map>
#include  <set>
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

set<string> SetLoadFile;

bool load(string ss)
{
  if(SetLoadFile.find(ss) != SetLoadFile.end())
    { 
      if( (mpirank==0)&& verbosity)
	cout << " (already loaded : " <<  ss << " ) " ;
    }
    else
      {
	SetLoadFile.insert(ss);
	bool ret=false;
	void * handle = 0;
	const int /*nbprefix=2,*/nbsuffix=2;
	list<string> prefix(ffenvironment["loadpath"]);
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
		  if(verbosity && (mpirank ==0))
		    cout << " (load: dlopen " << s << " " << handle << ") ";
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
		    if(verbosity&& (mpirank ==0))
		      cout << "(load: loadLibary " <<  s <<  " = " << handle << ")";
		    return mod;
	    }
	      }
#else
	      if((mpirank ==0))
		{
		  cout << "------------------------------------   \n" ;
		  cout << "  load: sorry no dlopen on this system " << s << " \n" ;
		  cout << "------------------------------------   \n" ;
		}
	CompileError("Error load");
	return 0;
#endif  
	    }
	if((mpirank ==0))
	  {
	    cerr  <<   "\nload error : " << ss << "\n \t fail : "  << endl;
	    cerr << "list  prefix: " ;
	    for (list<string>::const_iterator i= prefix.begin();i !=prefix.end();++i)
	      cerr <<"'"<<*i<<"' ";
	    cerr << "list  suffix : '"<< suffix[0] << "' , '"  << suffix[1] << "' "; 
	    
	    cerr << endl;
	  }
	CompileError("Error load");
      }
  return 0 ;
}

