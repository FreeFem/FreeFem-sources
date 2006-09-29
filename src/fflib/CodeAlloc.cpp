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
#include <iostream>
#include <typeinfo>
#include <cstddef>
#include <cstdlib>
#include <cassert>
#include "error.hpp"
using namespace std;

#include "CodeAlloc.hpp"

/*
size_t CodeAlloc::nb=0, CodeAlloc::lg=0,CodeAlloc::nbpx=0,CodeAlloc::chunk=2048; 
size_t CodeAlloc::nbt,CodeAlloc::nbdl=0;
CodeAlloc ** CodeAlloc::mem=0;
bool CodeAlloc::sort=true;
bool  CodeAlloc::cleanning=false;
*/
static long kerr=0;
static long nbsort =0;
template<class T>
static void  HeapSort(T *c,long n)
  {
    long l,j,r,i;
    T crit;
    c--; // on decale de 1 pour que le tableau commence a 1
    if( n <= 1) return;
    l = n/2 + 1;
    r = n;
    while (1) { // label 2
      if(l <= 1 ) { // label 20
	crit = c[r];
	c[r--] = c[1];
	if ( r == 1 ) { c[1]=crit; return;}
      } else  crit = c[--l]; 
      j=l;
      while (1) {// label 4
	i=j;
	j=2*j;
	if  (j>r) {c[i]=crit;break;} // L8 -> G2
	if ((j<r) && (c[j] < c[j+1])) j++; // L5
	if (crit < c[j]) c[i]=c[j]; // L6+1 G4
	else {c[i]=crit;break;} //L8 -> G2
      }
    }
  }

 void CodeAlloc::resize() 
{
  Sort_mem();
  if( nbt*1.5 +10 >= nbpx )
    {
      nbpx = chunk;
      mem=(CodeAlloc **)realloc(mem,chunk*sizeof(void*));
      if(mem)  nbpx=chunk;
      else ErrorExec("Alloc problem",0);
      chunk *= 3; 
      chunk /= 2;
      assert(chunk > nbpx);
    }
}

     
     
   void CodeAlloc::Sort_mem()
  { 
    int i,j;
    if(nbt ==0) return;
    nbsort++; 
    HeapSort(mem,nbt);
    
    for( i=0,j=0;i<nbt;i++)
      if ( ! isdel(i) )
        mem[j++]=mem[i];
    nbt=j;
    if(nbt != nb) 
      assert(nbt==nb);
  }
    void CodeAlloc::clear()
  {
     return ;

    cleanning=true;
    Sort_mem();
    cout << " Clear Alloc Mem: "<< nb << endl; ;
    for(int i=0;i<nbt;i++)
      {
        if( ! isdel(i)) 
	  delete mem[i];	
      }
    cout << " nb " << nb <<  " ==  " << 0 << " , nbsort: " << nbsort <<endl;
    nb=0;
    lg=0;
    nbpx=0;
    chunk=2048; 
    sort=true;
    free(mem);
    mem=0;
    nbdl=0;
    if(kerr)
    cout << " CodeAlloc: nb err delete  " << kerr << endl;
    kerr =0; 
     cleanning=false;
  }
   void CodeAlloc::ErrorDel(void *pp)
  {
//    cerr << " Alloc:: Sorry the pointeur " << pp << " is not allocated \n";    
     kerr++;  
  }
  void CodeAlloc::operator delete(void * pp) {
    //  return ; // FH 
    if( ! sort  &&  (nb*2 > nbt)  )  Sort_mem();
    int ib=0,ie=nbt-1,im;
    if(pp< mem[ib]) ErrorDel(pp);
    if(pp> mem[ie]) ErrorDel(pp);
    int p=-1;
    while( ib < ie )
      {
	im = (ib + ie)/2;
	  {
	    if(  pp < mem[im]) ie=im-1; // in [ib , im-1] 
	    else if (mem[im] == pp) { p=im;break;}
	    else ib =im+1;// in [im+1 , ie]
	  }
      }
    if( p <0 && mem[ib] == pp) p=ib;
    
    if(p<0) 
      ErrorDel(pp);
    else {
      setdel(p); //le tableau est detruit   
      nbdl++;nb--; 
      ::operator delete(pp);
    }
  }
