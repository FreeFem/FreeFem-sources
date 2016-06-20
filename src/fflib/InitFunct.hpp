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

#ifndef INITSFUNCT_HPP_
#define INITSFUNCT_HPP_
#include "ffapi.hpp"
 
// [[file:~/ff/src/fflib/InitFunct.cpp::addInitFunct]]
void  addInitFunct(int i,void  (* f)(),const char *name) ;

// [[file:InitFunct.cpp::callInitsFunct]]
void  callInitsFunct();

class  addingInitFunct {  public:
  addingInitFunct(int i,void  (* f)(),const char *name="") {
    // [[file:~/ff/src/fflib/InitFunct.cpp::addInitFunct]]
    addInitFunct(i,f,name);
  }
} ;

// <<LOADINITIO>> In Emscripten manipulating input and output streams is not allowed (cf
// [[file:~/fflib/Makefile::NO_STREAM_REDIRECT]])

#ifdef NO_STREAM_REDIRECT
#define LOADINITIO
#else
#if _WIN32
#define LOADINITIO {					\
    streambuf * so =ffapi::cout()->rdbuf() ;		\
    streambuf * si =ffapi::cin()->rdbuf() ;		\
    streambuf * se =ffapi::cerr()->rdbuf() ;		\
    if( so &&  cout.rdbuf() != so ) cout.rdbuf(so);	\
    if( si &&  cin.rdbuf() != si ) cin.rdbuf(si);	\
    if( se &&  cerr.rdbuf() != se ) cerr.rdbuf(se);	\
} 
#else // [[_WIN32]]
#define LOADINITIO {					\
    streambuf * so =ffapi::cout()->rdbuf() ;		\
    streambuf * si =ffapi::cin()->rdbuf() ;		\
    streambuf * se =ffapi::cerr()->rdbuf() ;		\
    if( so &&  cout.rdbuf() != so ) cout.rdbuf(so);	\
    if( si &&  cin.rdbuf() != si ) cin.rdbuf(si);	\
    if( se &&  cerr.rdbuf() != se ) cerr.rdbuf(se);	\
    stdout = ffapi::ffstdout();				\
    stderr = ffapi::ffstderr();				\
    stdin = ffapi::ffstdin();				\
} 
#endif // [[_WIN32]]
#endif // [[NO_STREAM_REDIRECT]]

#define LOADINITNM(EXEC,NM)					\
								\
  static  void  AutoLoadInit() {				\
								\
    /* [[LOADINITIO]] */					\
    LOADINITIO;							\
    if(verbosity>9) cout << "\n loadfile " NM  "\n" ;		\
    EXEC;							\
  }								\
								\
  static int DoLoadInit() {					\
    if(verbosity>9) cout << " ****  " << NM  <<  " ****\n" ;	\
								\
    /* <<calling_addInitFunct>> */				\
    /* [[file:~/ff/src/fflib/InitFunct.cpp::addInitFunct]] */	\
    addInitFunct(10000,&AutoLoadInit,NM);			\
    return 2;							\
  }								\
								\
  static int callDoLoadInit=DoLoadInit();				

#define LOADINIT(TI) LOADINITNM(TI init obsolete,__FILE__)			     
#define LOADFUNC(FC) LOADINITNM(FC() ,__FILE__)			     

#endif // [[INITSFUNCT_HPP_]]
