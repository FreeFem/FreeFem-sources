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
#ifndef SERIALEZE_HPP_
#define SERIALEZE_HPP_
#include <cstring>
#include "endian.hpp"
struct MPIrank;
 class Serialize {  
   // we store a refcounter in the pointer p a adresse p-sizeof(long)
   // so we can use the copy constructor
   private: 
   
  size_t lg;
  const char *what;
  char * p; 
  public: 
  Serialize(size_t lgg,const char * wht): 
    lg(lgg), what(wht) , p((new char[lg+sizeof(long)])+sizeof(long)) 
   { //cout << " begin count()=0 " << endl;
    count()=0; }
  
  ~Serialize(){ if(count()--==0) delete [] (p-sizeof(long));}
  size_t size() const { return lg;}
  //  mpi routine
  void mpisend(const MPIrank &,long tag);
  Serialize(const MPIrank &,const char * wht,long tag); 
  //  end mpi routine 
  operator void *() { return p;} 
  operator char *() { return p;} 
  bool samewhat(const char * w) const { return strncmp(what,w,8)==0; }

  Serialize(const Serialize & s) :
    lg(s.lg),
    what(s.what) ,  
    p(s.p)
   { count()++; } 
  
 template<typename T>  inline void get(size_t & k,T & x) const
   { 
     T xx;//= r_endian(x);
     assert(k<=lg+sizeof(T));
     memcpy(&xx,p+ k,sizeof(T));
     k +=  sizeof(T);
     x=r_endian(xx);
   }
   template<typename T>  inline void put(size_t & k,const T & x) 
   { 
     if ( !(k<=lg+sizeof(T)) )
       {
	 cout << " assert put " << k << " <=" << lg + sizeof(T) << endl;
	 assert((k<=lg+sizeof(T)));
       }
     T xx= w_endian(x);
     memcpy( p + k,&xx,sizeof(T));
     k += sizeof(T);
   }
   
 private:
   long & count() const  { return * (long*) (void*) (p-sizeof(long));}
   void operator=(Serialize & s) ; // no affectation

 };
#endif
