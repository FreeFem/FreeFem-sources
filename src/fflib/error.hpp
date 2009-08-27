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
#ifndef ERROR_H
#define ERROR_H
#include <cassert>
#include <string>
#include "throwassert.hpp"
#include <exception>

extern int TheCurrentLine; 

#if defined(__GNUC__) && __GNUC__+0 < 3
#include <strstream.h>
 typedef  istrstream istringstream   ;
 typedef  ostrstream ostringstream   ;
#define ENDS << '\0'
#define OLDCPP 1 
#else
// car ostringstream n'est pas encore defin sous g++
//  
#include <sstream>
#define ENDS 

#endif

  using std::exception;

extern long mpirank;

class Error : public exception
{ public:
  enum CODE_ERROR { NONE, COMPILE_ERROR,EXEC_ERROR, MEM_ERROR,MESH_ERROR,ASSERT_ERROR,INTERNAL_ERROR, UNKNOWN };
  
  
private: 
  std::string  message;
  const CODE_ERROR code;
protected:
  Error(CODE_ERROR c,const char * t1,const char * t2,const char * t3=0,
	int n=0,const char * t4=0,const char * t5=0,const char * t6=0,
	const char * t7=0,const char * t8=0,const char * t9=0)     
    : message(),code(c)
  {
    using namespace std;
    ostringstream mess;
    if(t1)  mess << t1;
    if(t2)  mess << t2;
    if(t3)  mess << t3 << n ;
    if(t4)  mess << t4;
    if(t5)  mess << t5;
    if(t6)  mess << t6;
    if(t7)  mess << t7;
    if(t8)  mess << t8;
    if(t9)  mess << t9;
    message = mess.str();
    extern void ShowDebugStack();
    ShowDebugStack();
    if (c!=NONE && mpirank==0) cerr  << message << endl; // cerr << " at exec line  " << TheCurrentLine << endl; 
    }
public:
  virtual int errcode() const {return code;} 
  virtual const char *  what() const   throw () { return message.c_str(); } 
  virtual  ~Error() throw () {}      
};

class ErrorCompile : public Error
{
 public:
  ErrorCompile(const char * Text,int l,const char * t2="") : 
    Error(COMPILE_ERROR,"Compile error : ",Text,"\n\tline number :",l,", ", t2) {}
};

class ErrorExec : public Error
{  
 public:
  ErrorExec(const char * Text,int l) :
    Error(UNKNOWN,"Exec error : ",Text, "\n   -- number :", l)  {}
};

class ErrorInternal : public Error
{  
 public:
  ErrorInternal(const char * Text,int l,const char * t2="") :
    Error(INTERNAL_ERROR,"Internal error : ",Text, "\n\tline  :",l,", in file ", t2)  {}
};
class ErrorAssert : public Error
{  
 public:
  ErrorAssert(const char * Text,const char *file,const int line) :
    Error(ASSERT_ERROR,"Assertion fail : (",Text, ")\n\tline :", line,", in file ",file)  {}
};


class ErrorMemory : public Error
{ public:
  ErrorMemory(const char * Text,int l=0) : 
    Error(MEM_ERROR,"Memory Error : ",Text," number: ",l)  {}
};

class ErrorExit : public Error
{
  int codeexit;
public:
  ErrorExit(const char * ,int l) :
    Error(NONE,"exit","(","",l,")"), codeexit(l)   {}
    // the exit code fo freefem++ is given by l 
  int errcode() const{return codeexit;}
};

#endif
