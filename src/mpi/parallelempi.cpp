#include <fstream>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <complex>
#include<stdlib.h>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
#include "ufunction.hpp"
using namespace std;  
#include "rgraph.hpp"
#include "RNM.hpp"
// after RNM   otherwise 
// trouble with index in RNM (I do no understander FH)
#include <set>
#include <vector>
#include <map>

#include "fem.hpp"


#include "FESpacen.hpp" 
#include "FESpace.hpp" 

#include "MatriceCreuse_tpl.hpp"
#include "MeshPoint.hpp"
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "Operator.hpp" 
#include "lex.hpp"
#include "libmesh5.h"
#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"




//  pas terrible .... 

#undef MPICH_SKIP_MPICXX
#define  MPICH_SKIP_MPICXX
#undef MPICH_IGNORE_CXX_SEEK
#define MPICH_IGNORE_CXX_SEEK
#include "mpi.h"
//  July change of cluster.cica.es parallele compute
//  where mpicxx.h is wrong and penmpi/ompi/mpi/cxx/mpicxx.h is good !!!
// F. Hecht 
#ifdef XXXXXXXXXXXXXXXXXXXXXXXXXXXXXZZZZ
#ifdef HAVE_OPENMPI_OMPI_MPI_CXX_MPICXX_H
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#elif HAVE_MPI_CXX_MPICXX_H
#include <mpi/cxx/mpicxx.h>
#elif HAVE_OMPI_MPI_CXX_MPICXX_H
#include <ompi/mpi/cxx/mpicxx.h>
#elif HAVE_MPICXX_H
#include <mpicxx.h>
#elif  HAVE_MPI___H
#include <mpi++.h>
#else
#error "no mpixx.h or mpi++.h file" 
#endif 
#endif
/*
 MPI::INT, MPI::LONG, MPI::SHORT, 
 MPI::UNSIGNED SHORT, MPI::UNSIGNED, 
 MPI::UNSIGNED LONG, MPI::SIGNED CHAR, 
 MPI::UNSIGNED CHAR 
 Fortran integer: MPI::INTEGER 
 Floating point: MPI::FLOAT, MPI::DOUBLE, MPI::REAL, 
 MPI::DOUBLE PRECISION, 
 MPI::LONG DOUBLE 
 Logical: MPI::LOGICAL, MPI::BOOL 
 Complex: MPI::F COMPLEX, MPI::COMPLEX, 
 MPI::F DOUBLE COMPLEX, 
 MPI::DOUBLE COMPLEX, 
 MPI::LONG DOUBLE COMPLEX 
 Byte: MPI::BYTE 
 
 */

// to send a sparse matrix we send header, line array ,colmun array, value array.
//  end afer the fist resquest we need to do allocation.
// so in this cas the communacatio is done in
// 2 step  
//    1 the header, 
//       alloc time 
//    2 the  message 
//  a couple request, pointer. 
//    Not use of IPROBE because probelem of wait. 

template<class T> struct MPI_TYPE {};
template<> struct MPI_TYPE<long>      {static const MPI_Datatype  TYPE(){return MPI_LONG;}};
template<> struct MPI_TYPE<int>      {static const MPI_Datatype TYPE(){return MPI_INT;}};
template<> struct MPI_TYPE<double>    {static const MPI_Datatype TYPE(){return MPI_DOUBLE;}};
template<> struct MPI_TYPE<char>    {static const MPI_Datatype TYPE(){return MPI_BYTE;}};

#ifdef MPI_DOUBLE_COMPLEX_
template<> struct MPI_TYPE<Complex>   {static const MPI_Datatype TYPE(){return MPI_DOUBLE_COMPLEX;}};
#endif
template<class T> struct MPI_WHAT {};
template<> struct MPI_WHAT<long>      {static const int WHAT=101;};
template<> struct MPI_WHAT<double>    {static const int WHAT=102;};
template<> struct MPI_WHAT<Complex>   {static const int WHAT=103;};
template<> struct MPI_WHAT<KN<long> >   {static const int WHAT=104;};
template<> struct MPI_WHAT<KN<double> >   {static const int WHAT=105;};
template<> struct MPI_WHAT<KN<Complex> >   {static const int WHAT=106;};

template<class T> struct MPI_TAG {};
template<> struct MPI_TAG<long>   {static const int TAG=5;};
template<> struct MPI_TAG<double>    {static const int TAG=4;};
template<> struct MPI_TAG<Complex >   {static const int TAG=6;};
template<> struct MPI_TAG<KN<long> >   {static const int TAG=11;};
template<> struct MPI_TAG<KN<double> >   {static const int TAG=12;};
template<> struct MPI_TAG<KN<Complex> >   {static const int TAG=13;};
template<> struct MPI_TAG<Mesh *>   {static const int TAG=1000;};
template<> struct MPI_TAG<Mesh3 *>   {static const int TAG=1010;};
template<> struct MPI_TAG<Matrice_Creuse<double> >   {static const int TAG=1020;};
template<> struct MPI_TAG<Matrice_Creuse<Complex> >   {static const int TAG=1030;};

  void initparallele(int &, char **&);
  void init_lgparallele();  
  
 extern long mpirank ;
 extern long mpisize ;


const size_t sizempibuf = 1024*32;

template<class R> 
void  WSend( R * v,int l,int who,int tag,MPI_Comm comm,MPI_Request *rq)
{
    MPI_Request rq0,*request=&rq0;
    if(verbosity>100)
    cout << mpirank<< " send to " << who << " tag " << tag << " " << rq << " " <<  comm <<endl;
    MPI_Isend((void *) v,l, MPI_TYPE<R>::TYPE() , who, tag,comm,request);
    if(rq) *rq=*request;
    else MPI_Request_free(request); 
}

template<> 
 void  WSend<Complex> ( Complex * v,int n,int who,int tag,MPI_Comm comm,MPI_Request *rq)
{
    MPI_Request rq0,*request=&rq0;
    if(verbosity>100)
    cout << mpirank<< " send to " << who << " tag " << tag << " " << rq << " " <<  comm <<  endl;

#ifdef MPI_DOUBLE_COMPLEX
    MPI_Isend(reinterpret_cast<void*> (v) , n, MPI_DOUBLE_COMPLEX, who, tag,comm,request);
#else
    n *=2;
    MPI_Isend(reinterpret_cast<void*> (v), n, MPI_DOUBLE, who, tag,comm,request);
    n /= 2;
#endif
    if(rq) *rq=*request;
    else MPI_Request_free(request);
}

template<class R> 
void  WRecv(R * v,int n,int who,int tag,MPI_Comm comm,MPI_Request *rq)
{
    MPI_Status status;
    if(rq) 
	MPI_Irecv(reinterpret_cast<void*> (v),n, MPI_TYPE<R>::TYPE() , who, tag,comm,rq);
    else 
        MPI_Recv(reinterpret_cast<void*> (v),n, MPI_TYPE<R>::TYPE() , who, tag,comm,&status);
}

template<> 
void  WRecv<Complex> (Complex * v,int n,int who,int tag,MPI_Comm comm,MPI_Request *rq)
{
    MPI_Status status;
#ifdef MPI_DOUBLE_COMPLEX_
    if(rq) 
      MPI_Irecv(reinterpret_cast<void*> (v), n, MPI_DOUBLE_COMPLEX, who, tag,comm,rq);
    else 
      MPI_Recv(reinterpret_cast<void*> (v), n, MPI_DOUBLE_COMPLEX, who, tag,comm,&status);
#else
    n *=2;
    if(rq) 
        MPI_Irecv(reinterpret_cast<void*> (v), n, MPI_DOUBLE, who, tag,comm,rq);
    else 
       MPI_Recv(reinterpret_cast<void*> (v), n, MPI_DOUBLE, who, tag,comm,&status);
    n /= 2;
#endif
}
			 
template<class R> 
void  WBcast(R * v,int  n,int who,MPI_Comm comm)  
{
    assert(v && n>0);
  MPI_Bcast(reinterpret_cast<void*> (v), n, MPI_TYPE<R>::TYPE(), who,comm);
}

template<> 
void  WBcast<Complex>(Complex * v,int  n,int who,MPI_Comm comm)  
{
   assert(v && n>0);
#ifdef MPI_DOUBLE_COMPLEX_
    MPI_Bcast(reinterpret_cast<void*> (v), n, MPI_TYPE<R>::TYPE(), who,comm);
#else
    n *=2;
    MPI_Bcast(reinterpret_cast<void*> (v), n, MPI_DOUBLE, who,comm);
    n /= 2;
#endif
}

			 

struct MPIrank {
    
   int who; 
   MPI_Comm comm; 
   MPI_Request *rq; 

   MPIrank(int i=0,MPI_Comm com=MPI_COMM_WORLD, MPI_Request *rqq=0) : who(i) , comm(com),rq(rqq) {ffassert(i>=-1 && i < mpisize);} 
    
    
   
   const MPIrank & operator<<(double a)const  {
     WSend(&a, 1, who, 4,comm,rq);
      return *this;
   }
   const MPIrank & Bcast(double & a) const {
      WBcast(&a, 1, who,comm);
      return *this;
   }
   const MPIrank & operator>>(double & a) const {
     WRecv(&a, 1, who, 4,comm,rq);
      return *this;
   }
   const MPIrank & operator<<(long a) const {
     WSend(&a, 1, who, 5,comm,rq);
      return *this;
   }
   const MPIrank & Bcast(long & a) const {
     WBcast(&a, 1,  who,comm);
      return *this;
   }
   const MPIrank & operator>>(long & a) const {
      WRecv(&a, 1,  who, 5,comm,rq);
      return *this;
   }
   

    template<class R>
    const MPIrank & operator>>(KN<R> & a) const {
	assert(&a);
	if(verbosity>99)
	  cout << " ---- " << who  << "  >> " << & a << " " << a.N() << " " << MPI_TAG<KN<R> >::TAG 
	       <<" from " << mpirank << "  "<<  (R *) a << endl;
	int n= a.N();
	WRecv((R *) a, n, who, MPI_TAG<KN<R> >::TAG ,comm,rq);
	if(verbosity>99)
	  cout << " ++++ " << who  << "  >> " << & a << " " << a.N() << " " << MPI_TAG<KN<R> >::TAG 
	       <<" from  " << mpirank << "  "<<  (R *) a << endl;
	ffassert(a.N()==n);
	return *this;
    }
     
    template<class R>
    const MPIrank & operator<<(const KN<R> *aa)const  {
	const KN<R> & a=*aa;
	ffassert(&a); 
	int n= a.N();
	if(verbosity>99)
	  cout << " .... " << who  << "  >> " << & a << " " << a.N() << " " << MPI_TAG<KN<R> >::TAG 
	       <<" from  " << mpirank << "  "<<  (R *) a << endl;
	WSend((R *) a,n,who,MPI_TAG<KN<R> >::TAG,comm,rq);
	return *this;
    }
    
   template<class R> 
   const MPIrank & Bcast(const KN<R> &a) const  {
    //const KN<R> & a=*aa;
      assert(&a); 
      int n= a.N();
      WBcast((R *) a, n, who,comm);
      ffassert(a.N()==n);
      return *this;
   }

    
   const MPIrank & Bcast(Fem2D::Mesh *&  a) const {
     if(verbosity>1) 
       cout << " MPI Bcast  (mesh *) " << a << endl;
      Serialize  *buf=0;
      int nbsize=0;
      if(  who == mpirank)  
        {
          buf =new Serialize((*a).serialize());
          nbsize =  buf->size();
        }
      WBcast( &nbsize, 1,  who,comm);
       if (who != mpirank)
          buf= new Serialize(nbsize,Fem2D::Mesh::magicmesh);
       assert(nbsize);
       if(verbosity>2) 
         cout << " size to bcast : " << nbsize << " mpirank : " << mpirank << endl;
       
       WBcast( (char *)(*buf),nbsize,  who,comm);     
        
       if(who != mpirank)
        {
         if (a) (*a).decrement();
          a= new Fem2D::Mesh(*buf);
         Fem2D::R2 Pn,Px;
         a->BoundingBox(Pn,Px);
         a->quadtree=new Fem2D::FQuadTree(a,Pn,Px,a->nv);
        }   
      delete buf;      
      return *this;
   }
 
    const MPIrank & Bcast(Fem2D::Mesh3 *&  a) const {
	if(verbosity>1) 
	    cout << " MPI Bcast  (mesh3 *) " << a << endl;
	Serialize  *buf=0;
	int nbsize=0;
	if(  who == mpirank)  
	  {
	      buf =new Serialize((*a).serialize());
	      nbsize =  buf->size();
	  }
	WBcast( &nbsize, 1,  who,comm);
	if (who != mpirank)
	    buf= new Serialize(nbsize,Fem2D::GenericMesh_magicmesh);
	assert(nbsize);
	if(verbosity>2) 
	    cout << " size to bcast : " << nbsize << " mpirank : " << mpirank << endl;
	
	WBcast( (char *)(*buf),nbsize,  who,comm);     
        
	if(who != mpirank)
	  {
	      if (a) (*a).decrement();
	      a= new Fem2D::Mesh3(*buf);
	      a->BuildGTree();
	  }   
	delete buf;      
	return *this;
    }
    
    template<class R>
    const MPIrank & Bcast(Matrice_Creuse<R> &  a) const 
    {
	if(verbosity>1) 
	    cout << mpirank <<  ":  MPI Bcast " << who << "  (Matrice_Creuse &) " << &a << " " << a.A << endl;
	MatriceMorse<R> *mA=0;
	int ldata[4]={0,0,0,0};
	if(  who == mpirank)  
	  {
	      if(a.A)
		{
	      mA= a.A->toMatriceMorse();
	      ldata[0]=mA->n;
	      ldata[1]=mA->m;
	      ldata[2]=mA->nbcoef;
	      ldata[3]=mA->symetrique;
	     // cout << mpirank << " ldata " << ldata[0] << " " << ldata[1] <<" " << ldata[2] << " " <<ldata[3] << endl;
		}
	  }
        int n4=4;
	WBcast( ldata,n4, who,comm); 
	//cout << mpirank << " after 4 " " ldata " << ldata[0] << " " << ldata[1] <<" " << ldata[2] << " " <<ldata[3] << endl;
	int n1= ldata[0]+1;
	if(  who != mpirank && ldata[0] )
	      mA= new MatriceMorse<R>(ldata[0],ldata[1],ldata[2],ldata[3]); 
	if(ldata[0]) 
	  {
	  // cout << mpirank << " " << who << " lg  " << mA->lg << " " << n1 << endl;
	   WBcast(  mA->lg,n1, who,comm);     
	   //cout << mpirank << " " << who << " cl  " << mA->cl << " " <<  mA->nbcoef << endl;
	   WBcast(  mA->cl,mA->nbcoef, who,comm);     
	   //cout << mpirank << " " << who << " a  " << mA->a << " " <<  mA->nbcoef << endl;
	   WBcast( mA->a,mA->nbcoef , who,comm);  
	  }
	if(  who != mpirank) 
          a.A.master(mA);
	else 
	  delete mA;      
	return *this;
    }
    
    
    template<class R>
    const MPIrank & operator<<(Matrice_Creuse<R> * const &  a) const 
    {
	//if(verbosity>1) 
	    cout << " MPI << (Matrice_Creuse *) " << a << endl;
	ffassert(rq==0) ; // 
	int tag = MPI_TAG<Matrice_Creuse<R> >::TAG;		       
	MatriceMorse<R> *mA=a->A->toMatriceMorse();
	int ldata[4];
	ldata[0]=mA->n;
	ldata[1]=mA->m;
	ldata[2]=mA->nbcoef;
	ldata[3]=mA->symetrique;
    
	cout << " ldata " << ldata[0] << " " << ldata[1] <<" " << ldata[2] << " " <<ldata[3] << endl;
	WSend( ldata,4, who, tag,comm,rq);
	WSend(  mA->lg,mA->n+1,  who, tag+1,comm,rq);     
	WSend( mA->cl,mA->nbcoef,  who, tag+2,comm,rq);     
	WSend(  mA->a,mA->nbcoef,   who, tag+3,comm,rq);  
	delete mA;
	return *this;
    }
    template<class R>
    const MPIrank & operator>>(Matrice_Creuse<R> &  a) const 
    {
	if(verbosity>1) 
	    cout << " MPI << (Matrice_Creuse *) " << a << endl;
	ffassert(rq==0) ; // 
	int tag =  MPI_TAG<Matrice_Creuse<R> >::TAG;		       
	int ldata[4];	
	WRecv( ldata,4, who, tag,comm,rq);
	MatriceMorse<R> *mA= new MatriceMorse<R>(ldata[0],ldata[1],ldata[2],ldata[3]); 
	WRecv(  mA->lg,mA->n+1,   who, tag+1,comm,rq);     
	WRecv(  mA->cl,mA->nbcoef,  who, tag+2,comm,rq);     
	WRecv(  mA->a,mA->nbcoef,  who, tag+3,comm,rq);  
	a.A.master(mA);
	return *this;
    }
    
   const MPIrank & operator<<(Fem2D::Mesh *  a) const {
     if(verbosity>1) 
     cout << " MPI << (mesh *) " << a << endl;
      ffassert(a);
       Serialize  buf=(*a).serialize();       
       buf.mpisend(*this,MPI_TAG<Mesh *>::TAG,static_cast<const void *>(this));
      return *this;
   }
    const MPIrank & operator<<(Fem2D::Mesh3 *  a) const {
	if(verbosity>1) 
	    cout << " MPI << (mesh3 *) " << a << endl;
	ffassert(a);
	Serialize  buf=(*a).serialize();       
	buf.mpisend(*this,MPI_TAG<Mesh3 *>::TAG,static_cast<const void *>(this));
	return *this;
    }
    
   const MPIrank & operator>>(Fem2D::Mesh *& a) const {
     if(verbosity>1) 
     cout << " MPI >> (mesh *) &" << a << endl;
       Serialize buf(*this,Fem2D::Mesh::magicmesh,MPI_TAG<Mesh *>::TAG,static_cast<const void *>(this));
      if (a) (*a).decrement();
      a= new Fem2D::Mesh(buf);
      //  add 3 line FH 08/12/2003  forget build quadtree sorry      
      Fem2D::R2 Pn,Px;
      a->BoundingBox(Pn,Px);
      a->quadtree=new Fem2D::FQuadTree(a,Pn,Px,a->nv);
      return *this;
   }
    
    const MPIrank & operator>>(Fem2D::Mesh3 *& a) const {
	if(verbosity>1) 
	    cout << " MPI >> (mesh3 *) &" << a << endl;
	Serialize buf(*this,Fem2D::GenericMesh_magicmesh,MPI_TAG<Mesh3 *>::TAG,static_cast<const void *>(this));
	if (a) (*a).decrement();
	a= new Fem2D::Mesh3(buf);
	//  add 3 line FH 08/12/2003  forget build quadtree sorry      
	a->BuildGTree();
	return *this;
    }
    
  operator int () const { return who;}     
};

void Serialize::mpisend(const MPIrank & rank,long tag,const void * vmpirank)
{
     const MPIrank * mpirank=static_cast<const MPIrank *> (vmpirank);
     MPI_Comm comm=mpirank->comm;
     MPI_Request *rq=mpirank->rq;
    ffassert(rq==0);
       char * pp = p-sizeof(long);
      long countsave=count(); // save count 
      count()=lg; // store length in count 
      int l=lg+sizeof(long);
      if(verbosity>1) 
         cout << " -- send from  " << mpirank << " to " << rank << " serialized " << what 
              <<   ", l=" << l << ", tag=" << tag << endl;
      if (l <=sizempibuf)
        WSend(pp,l, rank, tag,comm,rq);
      else {
         WSend(pp,sizempibuf,  rank, tag,comm,rq);
         WSend(pp+sizempibuf,l-sizempibuf, rank, tag+1,comm,rq);
      }
      if(verbosity>1) 
         cout << "    ok send is arrived " << endl;      
      count()=countsave; // restore count 
}

Serialize::Serialize(const MPIrank & rank,const char * wht,long tag,const void * vmpirank)
 :what(wht) 
{
    const MPIrank * mpirank=static_cast<const MPIrank *> (vmpirank);
    MPI_Comm comm=mpirank->comm;
    MPI_Request *rq=mpirank->rq;

      if(verbosity>1) 
         cout << " -- waiting " << mpirank << " from  " << rank << " serialized " << what 
              << " tag = " << tag <<  endl;
    ffassert(rq==0);
   char * buf= new char [sizempibuf];
   WRecv(buf, sizempibuf,  rank, tag,comm,rq);
   lg = * (long *) (void *) buf;
   int l=lg+sizeof(long);
   char * pp= new char[l]  ;
   if ( l <= sizempibuf) 
      memcpy(pp,buf,l);
   else 
      {
        memcpy(pp,buf,sizempibuf);
       WRecv(pp+sizempibuf,l-sizempibuf,  rank, tag+1,comm,rq)  ;       
      }
    
   if(verbosity>1) 
       cout << "    " << mpirank << " recived from " << rank << " serialized " << what <<   ", l=" 
            << l << ", tag=" << tag << endl;
   delete [] buf;
   p=pp+sizeof(long);
   count()=0;
   
}

template<class A>
struct Op_Readmpi : public binary_function<MPIrank,A*,MPIrank> {
  static MPIrank  f(MPIrank const  & f,A *  const  & a)  
   { 
     f >> *a;
     return f;
   }
};


template<class A>
struct Op_Writempi : public binary_function<MPIrank,A,MPIrank> {
  static MPIrank  f(MPIrank const  & f,A   const  &  a)  
   { 
     f << a;
     return f;
   }
};

template<class A>
struct Op_Bcastmpi : public binary_function<MPIrank,A*,MPIrank> {
  static MPIrank  f(MPIrank const  & f,A *  const  & a)  
   { 
     f.Bcast(*a);
     return f;
   }
};

template<class R>
struct Op_All2All : public binary_function<KN_<R>,KN_<R>,long> {
    static long  f( KN_<R>  const  & s, KN_<R>  const  &r)  
    { 
	MPI_Comm comm=MPI_COMM_WORLD;
	int mpisizew;
	MPI_Comm_size(comm, &mpisizew); /* local */ 
	int chunk = s.N()/mpisizew;
	ffassert(s.N()==mpisizew*chunk && r.N()==s.N());
	
	return MPI_Alltoall( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			     (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(), comm);	
    }
};

template<class R>
struct Op_Allgather : public binary_function<KN_<R>,KN_<R>,long> {
    static long  f( KN_<R>  const  & s, KN_<R>  const  &r)  
    { 
	MPI_Comm comm=MPI_COMM_WORLD;
	int mpisizew;
	MPI_Comm_size(comm, &mpisizew); /* local */ 
	int chunk = s.N()/mpisizew;
	ffassert(s.N()==mpisizew*chunk && r.N()==s.N());
	
	return MPI_Allgather( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			    (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(), comm);	
    }
};

template<class R>
struct Op_All2All3 : public ternary_function<KN_<R>,KN_<R>,MPI_Comm,long> {
    static long  f(Stack, KN_<R>  const  & s, KN_<R>  const  &r,MPI_Comm const & cmm )  
    { 
	MPI_Comm comm=cmm;
	int mpisizew;
	MPI_Comm_size(comm, &mpisizew); /* local */ 
	int chunk = s.N()/mpisizew;
	ffassert(s.N()==mpisizew*chunk && r.N()==s.N());
	
	return MPI_Alltoall( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			    (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(), comm);	
    }
};

template<class R>
struct Op_Allgather3 : public ternary_function<KN_<R>,KN_<R>,MPI_Comm,long> {
    static long  f(Stack, KN_<R>  const  & s, KN_<R>  const  &r,MPI_Comm const & cmm)  
    { 
	MPI_Comm comm=cmm;
	int mpisizew;
	MPI_Comm_size(comm, &mpisizew); /* local */ 
	int chunk = s.N()/mpisizew;
	ffassert(s.N()==mpisizew*chunk && r.N()==s.N());
	
	return MPI_Allgather( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			     (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(), comm);	
    }
};


template<class R>
struct Op_Scatter3 : public   ternary_function<KN_<R>,KN_<R>,MPIrank,long> {
    static long  f(Stack, KN_<R>  const  & s, KN_<R>  const  &r,  MPIrank const & root)  
    { 
	
	int mpisizew;
	MPI_Comm_size(root.comm, &mpisizew); 
	int chunk = s.N()/mpisizew;
	ffassert(s.N()==mpisizew*chunk && r.N()==s.N());
	
	return MPI_Scatter( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			     (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(),root.who,root.comm);	
    }
};
template<class R>
struct Op_Gather3 : public   ternary_function<KN_<R>,KN_<R>,MPIrank,long> {
    static long  f(Stack, KN_<R>  const  & s, KN_<R>  const  &r,  MPIrank const & root)  
    { 
	
	int mpisizew;
	MPI_Comm_size(root.comm, &mpisizew); 
	int chunk = s.N()/mpisizew;
	ffassert(s.N()==mpisizew*chunk && r.N()==s.N());
	
	return MPI_Gather( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			   (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(),root.who,root.comm);	
    }
};

			     

MPIrank mpiwho(long i) { return MPIrank(i);}
MPIrank mpiwho(long i,MPI_Comm comm) { return MPIrank(i,comm);}
MPIrank mpiwho_(const long &i,const MPI_Comm &comm,MPI_Request * const &rq) { return MPIrank(i,comm,rq);}
MPIrank mpiwho_(const long &i,MPI_Request * const &rq) { return MPIrank(i, MPI_COMM_WORLD ,rq);}

long mpiWait(MPI_Request * rq) { 
    MPI_Status status;
    return MPI_Wait(rq,&status);
}


long mpiWaitAny(KN<MPI_Request>* rq) { 
    MPI_Status status;
    int index;
    cout << "mpiWaitAny " <<rq->N() << " in "  <<  endl ;
    MPI_Waitany(rq->N(),*rq,&index,&status);
    cout << "mpiWaitAny " <<rq->N() << " out " << index <<  endl ;
     return index;
}

 MPIrank * set_copympi( MPIrank* const & a,const MPIrank & b)
 { *a=b;return a;}
 
  void initparallele(int &argc, char **& argv)
  {
     MPI_Init(&argc, &argv);
    
      int mpirank1,mpisize1;
     MPI_Comm_rank(MPI_COMM_WORLD, &mpirank1); /* local */ 
     MPI_Comm_size(MPI_COMM_WORLD, &mpisize1); /* local */ 
       
      mpirank = mpirank1;//MPI::COMM_WORLD.Get_rank();
      mpisize =mpisize1;// MPI::COMM_WORLD.Get_size();
     cout << "initparallele rank " <<  mpirank << " on " << mpisize << endl;
  }

long mpiSize(MPI_Comm  cmm) { 
    int s;
 //   MPI_Comm_rank(MPI_COMM_WORLD, &s); /* local */ 
    MPI_Comm_size(MPI_COMM_WORLD, &s); /* local */ 
    return s;
}
long mpiRank(MPI_Comm  cmm) { 
    int s;
    MPI_Comm_rank(MPI_COMM_WORLD, &s); /* local */ 
 //   MPI_Comm_size(MPI_COMM_WORLD, &s); /* local */ 
    return s;
}

AnyType InitializeGroup(Stack stack,const AnyType &x){
    MPI_Group *g=PGetAny<MPI_Group>(x);
    *g=0;
    MPI_Comm_group(MPI_COMM_WORLD, g);
    return  g;
}
AnyType DeleteGroup(Stack stack,const AnyType &x){
    MPI_Group *g=PGetAny<MPI_Group>(x);
    MPI_Group_free(g);
    return  Nothing;
}
AnyType InitializeComm(Stack stack,const AnyType &x){
    MPI_Comm *comm=PGetAny<MPI_Comm>(x);
    *comm=0;
    MPI_Comm_dup(MPI_COMM_WORLD, comm);
    return  comm;
}
AnyType DeleteComm(Stack stack,const AnyType &x){
    MPI_Comm *comm=PGetAny<MPI_Comm>(x);
    MPI_Comm_free(comm);
    return  Nothing;
}
AnyType InitializeRequest(Stack stack,const AnyType &x){
    MPI_Request *comm=PGetAny<MPI_Request>(x);
    *comm=0;
    
    return  comm;
}
AnyType DeleteRequest(Stack stack,const AnyType &x){
    MPI_Request *comm=PGetAny<MPI_Request>(x);
    if(comm) MPI_Request_free(comm);
    return  Nothing;
}
//  Hack to Bypass a bug in freefem FH  ... 
template<> 
class ForEachType<MPI_Group>:  public basicForEachType{public:// correction july 2009..... FH  Hoooo....  (Il y a un bug DUR DUR FH  ...) 
    ForEachType(Function1 iv=0,Function1 id=0,Function1 OOnReturn=0):basicForEachType(typeid(MPI_Group),sizeof(MPI_Group),0,0,iv,id,OOnReturn) { }
};

template<> 
class ForEachType<MPI_Comm>:  public basicForEachType{public:// coorection july 2009..... FH  Hoooo....  (Il y a un bug DUR DUR FH  ...) 
    ForEachType(Function1 iv=0,Function1 id=0,Function1 OOnReturn=0):basicForEachType(typeid(MPI_Comm),sizeof(MPI_Comm),0,0,iv,id,OOnReturn) {}
};

template<> 
class ForEachType<MPI_Request>:  public basicForEachType{public:// correction july 2009..... FH  Hoooo....  (Il y a un bug DUR DUR FH  ...) 
    ForEachType(Function1 iv=0,Function1 id=0,Function1 OOnReturn=0):basicForEachType(typeid(MPI_Request),sizeof(MPI_Request),0,0,iv,id,OOnReturn) {}
};
// end Hack  ... 
MPI_Group* def_group( MPI_Group* const & a, KN_<long>  const & b)
{
    cout << b.N() <<endl;
    for(int i=0;i<b.N();++i)
	cout << b[i] << endl;
    // ici def a .. 
  //  ffassert(0); //   A AFAIRE  //  pour arete le programm 
    return a;
}

template<typename K,typename KK>
AnyType ClearReturnpKK(Stack stack, const AnyType & a)
{
    // a ne faire que pour les variables local au return...
    //  pour l'instant on copie pour fqire mqrche 
    // a repense  FH  mqi 2009....
    KK * m = GetAny<KK * >(a);
    //   KN<K> *cm=new KN<K>(true, *m); bug quant KN est une variable global
    // KN<K> *cm=new KN<K>( *m); // on duplique le tableau comme en C++  (dur dur ?????? FH)
    m->increment();
    Add2StackOfPtr2FreeRC(stack,m);
    if(verbosity>400)
	cout << "ClearReturnpKK:: increment + Add2StackOfPtr2FreeRC nb ref  " <<  -m->next  << endl;
    return m;
}

template<typename K,typename KK,typename KK_>
AnyType ClearReturnpKK_(Stack stack, const AnyType & a)
{
    // il faut faire un copie du tableau 
    KK_ * m = GetAny<KK_ * >(a);
    KK *cm=new KK(*m); 
    
    Add2StackOfPtr2Free(stack,cm);// detruire la copie 
    if(verbosity>400)
	cout << "ClearReturnpKK_:: copie  Add2StackOfPtr2Free "  << endl;
    return (KK_ *) cm;
}
template<typename K,typename KK,typename KK_>
AnyType ClearReturnKK_(Stack stack, const AnyType & a)
{
    // il faut faire un copie du tableau 
    KK_  m = GetAny<KK_>(a);
    KK *cm=new KK(m); 
    
    Add2StackOfPtr2Free(stack,cm);// detruire la copie 
    if(verbosity>400)
	cout << "ClearReturnKK_:: copie  Add2StackOfPtr2Free   "  << endl;
    return SetAny<KK_>(*cm);
}
template<class RR,class A,class B>  
RR * get_elementp_(const A & a,const B & b){ 
    if( b<0 || a->N() <= b) 
      { cerr << " Out of bound  0 <=" << b << " < "  << a->N() << " array type = " << typeid(A).name() << endl;
      ExecError("Out of bound in operator []");}
return  &((*a)[b]);}

KN<MPI_Request> * set_init0( KN<MPI_Request> * const & a,const long & b)
   { 
   a->init(b);
      for(int i=0;i<b;++i)
	  (*a)[i]=MPI_REQUEST_NULL;
   return a;
   }

void init_lgparallele()
  {
    if(verbosity) cout << "parallelempi ";
     using namespace Fem2D;
     Dcl_TypeandPtr<MPIrank>(0);
      
     Dcl_TypeandPtr<MPI_Group>(0,0,InitializeGroup,DeleteGroup); 
     Dcl_TypeandPtr<MPI_Comm>(0,0,InitializeComm,DeleteComm);  
     zzzfff->Add("mpiGroup",atype<MPI_Group*>());
     zzzfff->Add("mpiComm",atype<MPI_Comm*>());

     map_type[typeid(MPIrank).name()]->AddCast(new E_F1_funcT<MPIrank,MPIrank*>(UnRef<MPIrank>));
     map_type[typeid(MPI_Group).name()]->AddCast(new E_F1_funcT<MPI_Group,MPI_Group*>(UnRef<MPI_Group>));
     map_type[typeid(MPI_Comm).name()]->AddCast(new E_F1_funcT<MPI_Comm,MPI_Comm*>(UnRef<MPI_Comm>));
       
     TheOperators->Add("<-", 
       new OneOperator2_<MPIrank*,MPIrank*,MPIrank>(&set_copympi));
      
      // constructor example ... 
     TheOperators->Add("<-", 
			new OneOperator2_<MPI_Group*,MPI_Group*,KN_<long> >(&def_group)); 
/*  code edp
     int[int] procs=[1,2,3];
     mpiGroup toto(procs);
 
 */
      
     Global.Add("processor","(",new OneOperator1<MPIrank,long>(mpiwho));
    Global.Add("processor","(",new OneOperator2<MPIrank,long,MPI_Comm>(mpiwho));
     TheOperators->Add(">>",
		       new OneBinaryOperator<Op_Readmpi<double> >,
		       new OneBinaryOperator<Op_Readmpi<long> > ,
		       new OneBinaryOperator<Op_Readmpi<KN<double> > > ,
		       new OneBinaryOperator<Op_Readmpi<KN<long> > > ,
		       new OneBinaryOperator<Op_Readmpi<KN<Complex> > > ,
		       new OneBinaryOperator<Op_Readmpi<Mesh *> > ,
		       new OneBinaryOperator<Op_Readmpi<Mesh3 *> > ,
		       new OneBinaryOperator<Op_Readmpi<Matrice_Creuse<R> > > ,
		       new OneBinaryOperator<Op_Readmpi<Matrice_Creuse<Complex> > > 
       );
     TheOperators->Add("<<",
       new OneBinaryOperator<Op_Writempi<double> >,
       new OneBinaryOperator<Op_Writempi<long> > ,
       new OneBinaryOperator<Op_Writempi<KN<double> * > > ,
       new OneBinaryOperator<Op_Writempi<KN<long> * > > ,
       new OneBinaryOperator<Op_Writempi<KN<Complex> * > > ,
       new OneBinaryOperator<Op_Writempi<Mesh *> > ,
       new OneBinaryOperator<Op_Writempi<Mesh3 *> > ,
       new OneBinaryOperator<Op_Writempi<Matrice_Creuse<R> * > > ,
       new OneBinaryOperator<Op_Writempi<Matrice_Creuse<Complex>* > > 
		       
       );
       
    Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<double> >);
    Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<long> >);
    Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<KN<double> > >);
    Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<KN<long> > >);
    Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<KN<Complex> > >);
    Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<Mesh *> >);
    Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<Mesh3 *> >);
      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<Matrice_Creuse<R> > >);
      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<Matrice_Creuse<Complex> > >);
      Global.Add("mpiAlltoall","(",new OneBinaryOperator<Op_All2All< long > >);
      Global.Add("mpiAlltoall","(",new OneBinaryOperator<Op_All2All< double > >);
      Global.Add("mpiAllgather","(",new OneBinaryOperator<Op_Allgather< long > >);
      Global.Add("mpiAllgather","(",new OneBinaryOperator<Op_Allgather< double > >);
      Global.Add("mpiAlltoall","(",new OneTernaryOperator3<Op_All2All3< long > >);
      Global.Add("mpiAlltoall","(",new OneTernaryOperator3<Op_All2All3< double > >);
      Global.Add("mpiAllgather","(",new OneTernaryOperator3<Op_Allgather3< long > >);
      Global.Add("mpiAllgather","(",new OneTernaryOperator3<Op_Allgather3< double > >);
      
      Global.Add("mpiScatter","(",new OneTernaryOperator3<Op_Scatter3< long > >);
      Global.Add("mpiScatter","(",new OneTernaryOperator3<Op_Scatter3< double > >);
      Global.Add("mpiGather","(",new OneTernaryOperator3<Op_Scatter3< long > >);
      Global.Add("mpiGather","(",new OneTernaryOperator3<Op_Gather3< double > >);
      
//      Global.Add("mpiAlltoall","(",new OneBinaryOperator<Op_All2Allm< long > >);
  //    Global.Add("mpiAlltoall","(",new OneBinaryOperator<Op_All2Allm< double > >);
      
 
     Global.New("mpirank",CPValue<long>(mpirank));
     Global.New("mpisize",CPValue<long>(mpisize));
     static MPI_Comm mpiWorld=MPI_COMM_WORLD;
     Global.New("mpiCommWorld",CPValue<MPI_Comm>(mpiWorld));   
      // add FH 
    
    Dcl_TypeandPtr<MPI_Request>(0,0,InitializeRequest,DeleteRequest); // bof bof ... 
    Dcl_TypeandPtr_<KN_<MPI_Request> ,KN<MPI_Request>*  > 
       (0,0,0,::Destroy<KN<MPI_Request> >,
	::ClearReturnKK_<MPI_Request,KN<MPI_Request>,KN_<MPI_Request> >,
	::ClearReturnpKK<MPI_Request,KN<MPI_Request> >);
      
    map_type_of_map[make_pair(atype<long>(),atype<MPI_Request>())]=atype<KN<MPI_Request>*>(); // vector
      
     zzzfff->Add("mpiRequest",atype<MPI_Request*>());
     Global.Add("mpiWtime","(",new OneOperator0<double>(MPI_Wtime));    
     Global.Add("mpiWtick","(",new OneOperator0<double>(MPI_Wtick));    
      Global.Add("processor","(",new OneOperator3_<MPIrank,long,MPI_Comm,MPI_Request*>(mpiwho_));
     Global.Add("processor","(",new OneOperator2_<MPIrank,long,MPI_Request*>(mpiwho_));
     Global.Add("mpiWait","(",new OneOperator1<long,MPI_Request*>(mpiWait));
     Global.Add("mpiWaitAny","(",new OneOperator1<long,KN<MPI_Request>*>(mpiWaitAny));
      Global.Add("mpiSize","(",new OneOperator1<long,MPI_Comm>(mpiSize)); 
      Global.Add("mpiRank","(",new OneOperator1<long,MPI_Comm>(mpiRank)); 
    
      TheOperators->Add("<-", 
			new OneOperator2_<KN<MPI_Request> *,KN<MPI_Request> *,long>(&set_init0)
						);
   atype<KN<MPI_Request>* >()->Add("[","",new OneOperator2_<MPI_Request*,KN<MPI_Request>*,long >(get_elementp_<MPI_Request,KN<MPI_Request>*,long>));    
      
  }
  void end_parallele()
   {
    MPI_Finalize();
    if(mpirank==0) cout << "FreeFem++-mpi finalize correctly .\n" << flush ; 
    else if(verbosity>5)  cout << '.' << endl ;
   }
//   MPI::COMM_WORLD.Recv(&msg, 1, MPI::INT, from, MPI::ANY_TAG);
//    MPI::COMM_WORLD.Isend(&msg, 1, MPI::INT, to, 4);
