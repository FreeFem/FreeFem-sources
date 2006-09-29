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

class CodeAlloc { public:

  static size_t nb,nbt,lg, nbdl,nbpx, chunk ;   
  static CodeAlloc ** mem;
  static bool cleanning;
  static void * lgmax;
  static bool sort;
  static bool isdel(int i)
  {
    return  ((char *) (void *)  mem[i] - (char *) 0) % 2 == 1;
  }
      
  static void setdel(int i)
  {
    mem[i] = (CodeAlloc *) (void *) (((char *) mem[i] - (char *) 0)+ 1);
  }
  static void resize(); 
  
 static  void * Add2CleanAtEnd(void * p)
  {
    if(p) {
      if(nbt>=nbpx) resize();
      if(nbt>0) sort = sort && mem[nbt-1] < p;
      nb++; 
      mem[nbt++]=(CodeAlloc*)p;  }
    return p;
  }
  
  void *operator new(size_t ll ) {
    lg+=ll;
    return Add2CleanAtEnd(::operator new(ll));} 

  
  static void Sort_mem();
  static  void clear();
  static void ErrorDel(void *pp);
  void operator delete(void * pp);
  virtual ~CodeAlloc() {}
  
};

template<class T> class CodeAllocT: public CodeAlloc{
  T * p;
  public: 
  CodeAllocT(int n): p(new T[n]) { assert(p);}
  static T * New(int n) { return (new CodeAllocT(n))->p;}
  ~CodeAllocT(){ delete [] p;}
};
