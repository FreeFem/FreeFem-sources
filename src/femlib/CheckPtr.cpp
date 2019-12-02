//#define SHOWALLOC
#define NCHECKPTR // BUg with MPI ??? FH
/*
#if __cplusplus >= 201103L
#define PARAM_THROW(i) (i)
#define PARAM_NOTHROW(i) (i,nothrow) noexcept
#define DCLPARAM_NOTHROW(i) (i,const std::nothrow_t& nothrow_constant) noexcept
#define DLDCLPARAM_NOTHROW(i) (i,const std::nothrow_t& nothrow_constant) noexcept
#define DCLPARAM_THROW(i) (i) noexcept
#else
#define PARAM_THROW(i) (i) throw(std::bad_alloc);
#define PARAM_NOTHROW(i) (i,nothrow) throw()
#define DCLPARAM_NOTHROW(i) (i,const std::nothrow_t& nothrow_constant) throw()
#define DLDCLPARAM_NOTHROW(i) (i,const std::nothrow_t& nothrow_constant) throw()
#define DCLPARAM_THROW(i) (i) throw(std::bad_alloc)
#endif
 98: // nothrow
 void* operator new (std::size_t size) throw (std::bad_alloc);
 void* operator new (std::size_t size, const std::nothrow_t& nothrow_value) throw();
 void operator delete[] (void* ptr) throw();
 void operator delete[] (void* ptr, const std::nothrow_t& nothrow_constant) throw();

   11:
 void* operator new (std::size_t size);
 void* operator new (std::size_t size, const std::nothrow_t& nothrow_value) noexcept;
 p
 */
#if __APPLE__
#include <malloc/malloc.h>
#elif __FreeBSD__
#include <stdlib.h>
#else
#include <malloc.h>
#endif

static long verbosity;


static long StorageUsed()
{
#if MALLOC_ZONE_SPECIFIC_FLAGS
    struct mstats mem1;
    mem1 = mstats();
    return mem1.bytes_used;
#elif M_MMAP_THRESHOLD
    struct mallinfo mem1;
    mem1=mallinfo();
    return mem1.uordblks;
#else
    return 0;
#endif

}
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
#include <cstdlib>
#include <cerrno>
#include <cstdio>
#include <new>

void debugalloc()
{ }

void debugunalloc()
{ static long count=0;
  //  debugalloc();
 count++;}


void exitalloc(int i)
{ static long count=0;
 count++;
 exit(1);
}

// Modif:         Juin 2001  for debuging  missing  delete point
//  --  THE MACRO
// TO SHOW ALL ALLOCATION
// #define SHOWALLOC
//  TO DEBUG ALL UN DELETE POINETUR


int UnShowAlloc =1;

int ShowAlloc(const char *s,size_t & lg);

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
       printf("@@@@@@@@@@@@@@@@@ Erreur myfree p= %p   l=%ul n=%d\n",p,(unsigned int) l,nordre);

   }
  }
}
/*
void* operator new(std::size_t ) throw(std::bad_alloc) ;
void* operator new[](std::size_t ) throw(std::bad_alloc) ;
void* operator new(std::size_t ,const std::nothrow_t& nothrow_constant) throw() ;
void* operator new[](std::size_t ,const std::nothrow_t& nothrow_constant) throw() ;

void operator delete[](void * ) throw();
void operator delete(void * ) throw();
void operator delete(void * ,const std::nothrow_t& nothrow_constant) throw();
void operator delete[](void * ,const std::nothrow_t& nothrow_constant) throw();
*/
const int N100 = 10000;
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
    AllocData(const AllocExtern::AllocData&);
    void operator=(const AllocExtern::AllocData&);
  };

private:

  static const long Maxundelptr = 2048;
    static size_t StorageUsage;
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
  static char filename[128];
  AllocData * NewAllocData();
  OneAlloc *Alloc();
  static   const int nblastalloc=20;
    static OneAlloc *  LastAlloc[nblastalloc];
public:
  static  bool CheckAlloc ;

  void * MyNewOperator(size_t ll,bool is_array );
  void MyDeleteOperator(void * pp,bool is_array);
  AllocExtern();
  ~AllocExtern();
  void init();
  int ShowAlloc( const char *s,size_t & lg);
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
size_t AllocExtern::StorageUsage=0;
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
bool AllocExtern::CheckAlloc = true;
char AllocExtern::filename[128] ="ListOfUnAllocPtr.bin";
AllocExtern::OneAlloc * AllocExtern::LastAlloc[nblastalloc];
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
  if (!f) exitalloc(1);
  NextFree =   f->p;
  return f;
}


void * AllocExtern::MyNewOperator(size_t ll,bool is_array)
{
  if(after_end || !CheckAlloc)  return malloc(ll);
  init();
  AllocExtern::OneAlloc * a = Alloc();
  a->p = mymalloc(ll);
  a->l = ll+1; // pour les allocation null
  a->n = ++NbAlloc;
  a->is_array = is_array;
  LastAlloc[NbPtr%nblastalloc] = a;

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
        if(verbosity>2) {
      printf("\t\tCheckPtrMem Full Exit(10) New Alloc %ld %p when %ld\n ", ll, a->p, a->n);
      printf ("\t\tCheckPtr:Max Memory used %10.3f kbytes " ,  MaxUsedSize/1024. );
      printf (" Memory undelete %ld \n" , AllocSize);
        }
      exitalloc(1);
    }
  return (void*) ((char*)a->p);
}

void AllocExtern::MyDeleteOperator(void * pp,bool is_array)
{
  if(after_end) { /*free(pp)*/; return;}
  if( !AllocExtern::CheckAlloc) {free(pp); return;}

  init();
    for( int i=0; i< nblastalloc; ++i)
    {
        if( (NbPtr-i) <0) break;
        int k =(NbPtr-i)%nblastalloc;
         if (LastAlloc[k]  && LastAlloc[k]->l>0 && LastAlloc[k]->p ==pp)
         {
             AllocExtern::OneAlloc * pa=LastAlloc[k];
             LastAlloc[k]=0;
             size_t ll = pa->l-1;
             for (size_t kkk=0;kkk<ll;kkk++)
                 ((char *) pp)[kkk]=18;

             myfree((char*)pp,ll,pa->n);
#ifdef SHOWALLOC
             printf("\t%d\tCheckPtr: fast delete  Alloc %ld %lx when %ld \n",pa->n,pa->l-1,  pa->p, pa->n);
#endif

             AllocSize -= ll;
             NbPtr--;
             pa->l=0;
             pa->p = NextFree;
             pa->n =0;
             if (pa->is_array != is_array)
                 printf("\t\tCheckPtr:  erreur delete [] \n");
             //if( pa->n < NbAllocShow )		  debugalloc();
             NextFree = & pa->p;
             return;
         }
    }
  if (AllocHead)
    {
      AllocExtern::AllocData *p = AllocHead;
        int kloop=0;
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
		  printf("\t\tCheckPtr:  erreur delete [] \n");
		//if( p->a[i].n < NbAllocShow )		  debugalloc();
		NextFree = & p->a[i].p;
		return;}
                if(p == p->next ||  kloop++ > 10000)
                  {
                   printf("\n\n ****CheckPTR Big BUG ??????? loop %d\n",kloop);
                   break;
                  }

	  p = p->next;
	}
      if(pp)
	{
	  printf( "\t\tCheckPtr: delete of bad pointer %p  -----------\n",pp);
	  debugalloc();
	}

    } else
      myfree((char*)pp);
}
void AllocExtern::init()
{
   static int count=0;
   if(0== (count++))
    {
      sprintf(filename,"ListOfAllocPtr-%d.bin",(int) sizeof(void*));
      StorageUsage=0;
      AllocSize =0;
      MaxUsedSize =0;
      AllocHead =0;
      NbAlloc =0;
      NbPtr =0;
      NextFree =0;
      NbuDelPtr =0;
      NbuDelPtr = 0;

      after_end = false;
      for(int i=0; i<nblastalloc; i++)
          LastAlloc[i]=0;
      FILE *file=fopen(filename,"rb");

      if (file)
	{
	  fread(&NbuDelPtr,sizeof(long),1,file);
	  fread(uDelPtr,sizeof(long),NbuDelPtr,file);
	  if(NbuDelPtr> 100000000 && NbuDelPtr <0)
	    {
	      printf("Fatal error in the file %s is wrong (please remove)",filename);
	      exit(1);
	    }
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
    if(UnShowAlloc==0) return;
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
      p = p->next;
     }
     k+=kk;
    kk=kk < Maxundelptr ? kk : Maxundelptr;
    HeapSort(list,kk);
    if(verbosity > 2)
    for (int i= kk-10<0 ? 0 : kk-10 ;i<kk;i++)
      {
        printf ("\t\tCheckPtr:Undelete pointer  %p size %ld  when %ld\n", list[i]->p,list[i]->l,list[i]->n);
      }
    if (kk)
      {
	FILE *file=fopen(filename,"wb");
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
    
      if(verbosity>2) {
    if(k)  printf ("\t\tCheckPtr:Nb of undelete pointer is %d last %d\n",k,lln);
    printf ("\t\tCheckPtr:Max Memory used %10.3f kbytes " ,  MaxUsedSize/1024. );
    printf (" Memory undelete %ld \n" , AllocSize);
      }

    //   clean store pointer
    p=AllocHead;
    while (p)
      {
	myfree((char*)p->a);
	AllocData * pold = p;
	p = p->next;
	myfree((char*)pold);
      }
    AllocHead=0;
    after_end=true;
}
// ------------------


void * operator new(std::size_t ll,const std::nothrow_t& nothrow_constant) throw()
{ void * p =  AllocExternData.MyNewOperator(ll,false);
 if (ll && !p) { printf("EXIT BECAUSE MEMORY FULL \n");
 exitalloc(1); };
 return p;}

void * operator new[](std::size_t ll,const std::nothrow_t& nothrow_constant) throw()
{ void * p =  AllocExternData.MyNewOperator(ll,true);
 if (ll && !p) { printf("EXIT BECAUSE MEMORY FULL \n");
 exitalloc(1); };
 return p;}

void *operator new(std::size_t ll) throw(std::bad_alloc)
{ void * p =  AllocExternData.MyNewOperator(ll,false);
    if (ll && !p) { printf("THROW BECAUSE MEMORY FULL \n");
        std::bad_alloc exception;
        throw exception; };
    return p;}

void *operator new[] (std::size_t ll) throw(std::bad_alloc)
{ void * p =  AllocExternData.MyNewOperator(ll,true);
    if (ll && !p) {
        printf("THROW BECAUSE MEMORY FULL \n");
        std::bad_alloc exception;
        throw exception; };
    return p;}


void operator delete(void * pp,const std::nothrow_t& nothrow_constant) throw()
{  AllocExternData.MyDeleteOperator(pp,false);}
void operator delete[] (void * pp,const std::nothrow_t& nothrow_constant) throw()
{  AllocExternData.MyDeleteOperator(pp,true);}
void operator delete(void * pp) throw()
{  AllocExternData.MyDeleteOperator(pp,false);}
void operator delete[](void * pp) throw()
{  AllocExternData.MyDeleteOperator(pp,true);}

int AllocExtern::ShowAlloc(const char *s,size_t & lg) {
    size_t m =StorageUsage;
    StorageUsage =StorageUsed();
    if (!NbAllocShow) {NbAllocShow=NbAlloc;}
    if(verbosity > 2)
  printf ("----------CheckPtr:-----%s------ NbUndelPtr  %ld  Alloc: %ld  NbPtr %ld  Mem Usage: %zu diff: %ld\n",s,NbPtr,AllocSize,NbAlloc,StorageUsage,(long)(StorageUsage-m));
  lg = AllocSize;
  return NbPtr;
}
int ShowAlloc(const char *s,size_t & lg)
{  return  AllocExternData.ShowAlloc(s,lg);}
void WithoutCheckPtr(){AllocExtern::CheckAlloc=false;}
void WithCheckPtr(){AllocExtern::CheckAlloc=true;}

#else
#define XXXX
#ifdef XXXX
#include <cstdlib>
#include <cerrno>
#include <cstdio>
#include <new>
#include <iostream>

long  CheckPtr___nbptr=0;
size_t CheckPtr___memoryusage =0;

void * operator new(std::size_t size,const std::nothrow_t& nothrow_constant) throw()
{
    CheckPtr___nbptr++;
    void *p = malloc( size );
    if(verbosity > 1000000 )
        std::cout << " CheckPtr: new nothrow" << CheckPtr___nbptr << " " << size
                  << " p =" << p <<std::endl;

    return p;
}

void * operator new[](std::size_t size,const std::nothrow_t& nothrow_constant) throw()
{
    void *p = malloc(size);
     CheckPtr___nbptr++;
    if(verbosity > 1000000 )
        std::cout << " CheckPtr: new[] nothrow" << CheckPtr___nbptr << " " << size << " p =" << p <<std::endl;
    return p;
}

void *operator new(std::size_t size) throw(std::bad_alloc)
{
    CheckPtr___nbptr++;
    void *p = malloc( size );
    if(verbosity > 1000000 )
        std::cout << " CheckPtr: new throw " << CheckPtr___nbptr << " " << size
        << " p =" << p <<std::endl;

    return p;
}

void *operator new[](std::size_t size) throw(std::bad_alloc)
{
    void *p = malloc(size);
    CheckPtr___nbptr++;
    if(verbosity > 1000000 )
        std::cout << " CheckPtr: new[]  throw" << CheckPtr___nbptr << " " << size << " p =" << p <<std::endl;
    return p;
}

void operator delete(void * p,const std::nothrow_t& nothrow_constant) throw()
{
    if(verbosity > 1000000 )
        std::cout << " CheckPtr: free nothrow " << CheckPtr___nbptr-1 <<  " p =" << p <<std::endl;

    free(p);

    CheckPtr___nbptr--;

}

void operator delete[](void * p,const std::nothrow_t& nothrow_constant) throw()
{
    if(verbosity > 1000000 )
        std::cout << " CheckPtr: free[] thow" << CheckPtr___nbptr-1 <<  " p =" << p <<std::endl;
    free(p);
    CheckPtr___nbptr--;

}
void operator delete(void * p) throw()
{
    if(verbosity > 1000000 )
        std::cout << " CheckPtr: free throw " << CheckPtr___nbptr-1 <<  " p =" << p <<std::endl;

    free(p);

    CheckPtr___nbptr--;

}

void operator delete[](void * p) throw()
{
    if(verbosity > 1000000 )
        std::cout << " CheckPtr: free[] throw " << CheckPtr___nbptr-1 <<  " p =" << p <<std::endl;
    free(p);
    CheckPtr___nbptr--;

}

int ShowAlloc(const char *s,size_t & lg)
{
    size_t m =StorageUsed();
    long diff = m-CheckPtr___memoryusage;
    if(verbosity > 0 && CheckPtr___memoryusage!=0 && m != CheckPtr___memoryusage)
        printf("CheckPtr:  Warning memory leak with malloc = %ld \n ",diff);
    CheckPtr___memoryusage=m;
    lg = 0; return CheckPtr___nbptr;}
void WithoutCheckPtr(){}
void WithCheckPtr(){}

int UnShowAlloc =0;
#else
#include <cstdlib>

int ShowAlloc(const char *s,size_t & lg)
{lg=0; return 0;}
int UnShowAlloc =0;
        void WithoutCheckPtr(){}
        void WithCheckPtr(){}

#endif
#endif
