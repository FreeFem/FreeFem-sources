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
 
void  addInitFunct(int i,void  (* f)(),const char *name) ;
void  callInitsFunct() ;

class  addingInitFunct {  public:
  addingInitFunct(int i,void  (* f)(),const char *name="") { addInitFunct(i,f,name);}
} ;

//

#define LOADINITIO {					\
    streambuf * so =ffapi::cout()->rdbuf() ;		\
    streambuf * si =ffapi::cin()->rdbuf() ;		\
    streambuf * se =ffapi::cerr()->rdbuf() ;		\
    if( so &&  cout.rdbuf() != so ) cout.rdbuf(so);	\
    if( si &&  cin.rdbuf() != si ) cin.rdbuf(si);	\
    if( se &&  cerr.rdbuf() != se ) cerr.rdbuf(se);	\
} 
  
#define LOADINITNM(EXEC,NM)						\
  static  void  AutoLoadInit() { LOADINITIO ;				\
    if(verbosity>9) cout << "\n loadfile " NM  "\n" ;			\
    EXEC; }								\
  int DoLoadInit() {							\
    if(verbosity>9)							\
      cout << " ****  " << NM  <<  " ****\n" ;				\
    addInitFunct(10000,&AutoLoadInit,NM);				\
    return 2;}								\
									\
  static int callDoLoadInit=DoLoadInit();				

#define LOADINIT(TI) LOADINITNM(TI init,__FILE__)			     
#define LOADFUNC(FC) LOADINITNM(FC() ,__FILE__)			     


#endif
