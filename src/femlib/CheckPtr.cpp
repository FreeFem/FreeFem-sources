//#undef NCHECKPTR
#ifndef NCHECKPTR
#define DEBUGUNALLOC 1 
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

void debugunalloc()
{ static long count=0;
 count++;}

void debugalloc()
{ static long count=0;
 count++;}

// AUTHOR:   F. Hecht,    
// ORG    :  UPMC  
// E-MAIL :   Frederic.Hecht@Inria.fr   
// ------------
// ORIG-DATE:     Dec 97
// Modif:         Juin 2001  for debuging  missing  delete point
//  --  THE MACRO 
// TO SHOW ALL ALLOCATION
// #define SHOWALLOC
//  TO DEBUG ALL UN DELETE POINETUR

#include <stdlib.h>
#include <cerrno>
#include <stdio.h>
#include "error.hpp"

int ShowAlloc(char *s,size_t & lg); 
        
//inline void *operator new(size_t, void *place) { return place; }

static int kerr=0;
void * mymalloc(size_t l)
{
  char *p = (char*)malloc(l+16);
  if (!p) return p;   
  for( int i = 0; i < 8 ; i++)
    p[i] = 'a'+i,p[i+l+8] = 'z'-i; // put a marque before 
  return (void*) (p+8);
}
void myfree(char *p,size_t l=0,int nordre=0) 
{
  if(p) {
    p -= 8;
    int k =0;
    for( int i = 0; i < 8 ; i++)
      {
	if (p[i] != 'a' +i)     k++;
	if(l && (p[i+l+8] != 'z' -i)) k++;       
      }
    for (size_t i=0;i<l;i++)
      p[i+8]=127;
    if(!k) free(p);
   else {
     debugalloc();
     if (kerr++<20) 
       printf("@@@@@@@@@@@@@@@@@ Erreur myfree p= %p   l=%d n=%d\n",p,(int) l,nordre);
     
     //throw(ErrorExec("exit",1));
   }
  }
}
void *operator new(size_t) throw (std::bad_alloc);
void operator delete(void * pp )throw () ;



const int N100 = 100;
class  AllocData;

class AllocExtern {
  public:
class OneAlloc {public:
  void * p;
  size_t l;
  long n;
  bool is_array ;
  bool operator<(const OneAlloc & a) const { return n<a.n;}
};
  
  void  HeapSort(OneAlloc **c,long n)
  {
    long l,j,r,i;
    OneAlloc* crit;
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
      if ((j<r) && (*c[j] < *c[j+1])) j++; // L5
      if ( *crit < *c[j]) c[i]=c[j]; // L6+1 G4
      else {c[i]=crit;break;} //L8 -> G2
      }
    }
  }
  
  class AllocData {public:
    OneAlloc *a;
    AllocData * next;
    AllocData();
    ~AllocData();
private:
    AllocData(const AllocData&);
    void operator=(const AllocData&);
  };

private:
  static const long Maxundelptr = 2048;
  static size_t AllocSize ;
  static size_t MaxUsedSize;
  static AllocData * AllocHead ;  
  static long NbAlloc;
  static long NbAllocShow;
  static long NbPtr;
  static void * NextFree;
  static long NbuDelPtr;
  static long uDelPtr[Maxundelptr];
  static bool after_end; 

  AllocData * NewAllocData();
  OneAlloc *Alloc();
public:
  
  void * MyNewOperator(size_t ll,bool is_array );
  void MyDeleteOperator(void * pp,bool is_array);
  AllocExtern();
  ~AllocExtern();
  void init();  
  int ShowAlloc(char *s,size_t & lg); 
  bool IsUnDelPtr(long nn) { // dichotomic find 
    long i=0;
    long j=NbuDelPtr-1;
    while (i<=j) {
      long k = (i+j)/2, kn=uDelPtr[k];
      if ( nn<kn) j=k-1;
      else if ( nn > kn) i = k+1;
      else
        return  true;}
    return false;
  }
};

static AllocExtern AllocExternData;

size_t AllocExtern::AllocSize =0;
size_t AllocExtern::MaxUsedSize =0;
AllocExtern::AllocData * AllocExtern::AllocHead =0;  
long AllocExtern::NbAlloc =0;
long AllocExtern::NbAllocShow=0;
long AllocExtern::NbPtr =0;
void * AllocExtern::NextFree =0;
long AllocExtern::NbuDelPtr =0;
long AllocExtern::uDelPtr[Maxundelptr];
bool AllocExtern::after_end =false;


AllocExtern::AllocData * AllocExtern::NewAllocData()
{
  
  AllocExtern::AllocData * ad = (AllocData *) mymalloc(sizeof(AllocData));
  ad->a = (OneAlloc*) mymalloc(sizeof(OneAlloc)*N100);
  for (int i=0;i<N100;i++)
    ad->a[i].l=0,ad->a[i].p=NextFree,NextFree = & ad->a[i];
  ad->next = AllocHead;
  AllocHead = ad;
#ifdef SHOWALLOC    
  printf("\t\tCheckPtr: OneAlloc[100] %lx\n",this);
#endif    
  return ad;
}



AllocExtern::OneAlloc * AllocExtern::Alloc()
{
  OneAlloc * f =  (OneAlloc *) NextFree;
  if (!f) 
    AllocHead = NewAllocData();
  f =(OneAlloc *) NextFree;
  if (!f) throw(ErrorExec("exit",1));
  NextFree =   f->p;
  return f;
}


void * AllocExtern::MyNewOperator(size_t ll,bool is_array)
{ 
  if(after_end) return malloc(ll);
  init();
  AllocExtern::OneAlloc * a = Alloc();
  a->p = mymalloc(ll);
  a->l = ll+1; // pour les allocation null
  a->n = ++NbAlloc;
  a->is_array = is_array;
  NbPtr++;
  AllocSize += ll;
#ifdef DEBUGUNALLOC
  if ( (IsUnDelPtr(a->n) && (a->n >= DEBUGUNALLOC) )) 
    debugunalloc();
#endif
#ifdef SHOWALLOC    
  printf( "\t%d\tCheckPtr: New Alloc %ld %lx when %ld\n ",a->n, ll, a->p, a->n);
#endif
  MaxUsedSize = AllocSize < MaxUsedSize ? MaxUsedSize :  AllocSize;
  if( !ll &&  !a->p)
    {
      printf("\t\tCheckPtrMem Full Exit(10) New Alloc %ld %p when %ld\n ",(long) ll, a->p, a->n);
      printf ("\t\tCheckPtr:Max Memory used %10.3f kbytes " ,  MaxUsedSize/1024. );
      printf (" Memory undelete %ld \n" , AllocSize);
      throw(ErrorExec("exit",10));
    }
  return (void*) ((char*)a->p);
}

void AllocExtern::MyDeleteOperator(void * pp,bool is_array)
{
  if(after_end) { /*free(pp)*/; return;}
  init();
  if (AllocHead)
    {
      AllocExtern::AllocData *p = AllocHead;
      while (p)
	{
	  for (int i=0;i<N100;i++)
	    if((p->a[i].l > 0) && (p->a[i].p == pp))
	      {
#ifdef SHOWALLOC    	  
		printf("\t%d\tCheckPtr: delete  Alloc %ld %lx when %ld \n",p->a[i].n,p->a[i].l-1,  p->a[i].p, p->a[i].n);
#endif
		size_t ll = p->a[i].l-1;
		for (size_t kkk=0;kkk<ll;kkk++) 
		  ((char *) pp)[kkk]=18;
		
		myfree((char*)pp,ll,p->a[i].n);
		
		AllocSize -= ll;
		NbPtr--;
		p->a[i].l=0;
		p->a[i].p = NextFree;
		p->a[i].n =0;
		if (p->a[i].is_array != is_array)
		  printf("\t\tCheckPtr:  erreur delete [] ");
		if( p->a[i].n < NbAllocShow )
		  debugalloc();      
		NextFree = & p->a[i].p;
		return;}
	  p = p->next;
	}
      if(pp) 
	printf( "\t\tCheckPtr: delete of bad pointer %p -----------\n",pp);
      
    } else 
      myfree((char*)pp); 
}
void AllocExtern::init()
{
   static int count=0;
   if(0== (count++)) 
    {
      AllocSize =0;
      MaxUsedSize =0;
      AllocHead =0;  
      NbAlloc =0;
      NbPtr =0;
      NextFree =0;
      NbuDelPtr =0;
      NbuDelPtr = 0;
      after_end = false;

      FILE *file=fopen("ListOfUnAllocPtr.bin","rb");
      
      if (file) 
	{
	  fread(&NbuDelPtr,sizeof(long),1,file);
	  fread(uDelPtr,sizeof(long),NbuDelPtr,file);
	  fclose(file);
	}  
      else
	{ // printf("fopen ListOfUnAllocPtr errno = %d\n",errno);
	}   
    }
}
AllocExtern::AllocExtern()
{
  init();
  
}

AllocExtern::~AllocExtern()
{
     OneAlloc *  list[Maxundelptr];
     
     AllocData * p=AllocHead;
     int k=0,kk=0;
     int lln=0;
     
     while (p) {int i=N100;
     while(i--)
       if (p->a[i].l >0  )
	 {
	   if ( p->a[i].n >=  p->a[i].n) lln = p->a[i].n;
	   if ( p->a[i].n <= NbAllocShow )
	     k++;
	   else
	     if (kk<Maxundelptr)
	       list[kk++]=p->a+i;
	 }
     // myfree((char*)p->a);
     //     AllocData * pold = p;
      p = p->next;
      //   myfree((char*)pold);
     }
     k+=kk;
    kk=kk < Maxundelptr ? kk : Maxundelptr;
    HeapSort(list,kk);
    for (int i= kk-10<0 ? 0 : kk-10 ;i<kk;i++)
      {
        printf ("\t\tCheckPtr:Undelete pointer  %p size %ld  when %ld\n", list[i]->p,list[i]->l,list[i]->n);        
      }
    if (kk)
      {
	FILE *file=fopen("ListOfUnAllocPtr.bin","wb");
	if (file) 
	  {
	    NbuDelPtr=kk;
	    for (int i=0;i<kk;i++)
	      uDelPtr[i]=list[i]->n;
	    fwrite(&NbuDelPtr,sizeof(long),1,file);
	    fwrite(uDelPtr,sizeof(long),NbuDelPtr,file);
	    fclose(file);
	  }
      }
      else {}
    
    if(k)  printf ("\t\tCheckPtr:Nb of undelete pointer is %d last %d\n",k,lln);
    printf ("\t\tCheckPtr:Max Memory used %10.3f kbytes " ,  MaxUsedSize/1024. );
    printf (" Memory undelete %ld \n" , AllocSize);
    
    //   clean store pointer      
    p=AllocHead;    
    while (p) 
      {//int i=N100;
      myfree((char*)p->a);
      AllocData * pold = p;
      p = p->next;
      myfree((char*)pold);
      }     
    AllocHead=0;
    after_end=true; 
}
// ------------------


void *operator new(size_t ll ) throw (std::bad_alloc)
{ void * p =  AllocExternData.MyNewOperator(ll,false);
 if (ll && !p) { printf("EXIT BECAUSE MEMORY FULL \n");
 throw(ErrorExec("exit",-1));};
 return p;}
void *operator new[](size_t ll ) throw (std::bad_alloc)
{ void * p =  AllocExternData.MyNewOperator(ll,true);
 if (ll && !p) { printf("EXIT BECAUSE MEMORY FULL \n");
 throw(ErrorExec("exit",-1));};
 return p;}
  
void operator delete(void * pp)throw ()
{  AllocExternData.MyDeleteOperator(pp,false);}
void operator delete[](void * pp)throw ()
{  AllocExternData.MyDeleteOperator(pp,true);}

int AllocExtern::ShowAlloc(char *s,size_t & lg) {
  // AllocExtern::AllocData * p=AllocExtern::AllocHead;
  if (!NbAllocShow) NbAllocShow=NbAlloc;
  //int i=N100-1;
  printf ("----------CheckPtr:-----%s------ NbUndelPtr  %ld  Alloc: %ld  NbPtr %ld \n",s,NbPtr,AllocSize,NbAlloc);
  lg = AllocSize;
  return NbPtr;
}
int ShowAlloc(char *s,size_t & lg) 
{  return  AllocExternData.ShowAlloc(s,lg);}
#else
#include <stdlib.h>
int ShowAlloc(char *s,size_t & lg)
{lg=0; return 0;}
#endif
