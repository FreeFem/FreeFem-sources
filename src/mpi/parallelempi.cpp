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
MPI_Request  WSend( R * v,int l,int who,int tag,MPI_Comm comm)
{
    MPI_Request request;
    MPI_Isend((void *) v,l, MPI_TYPE<R>::TYPE() , who, tag,comm,&request);
    return request;
}

template<> 
 MPI_Request  WSend<Complex> ( Complex * v,int n,int who,int tag,MPI_Comm comm)
{
   MPI_Request request;
#ifdef MPI_DOUBLE_COMPLEX
    MPI_Isend(reinterpret_cast<void*> (v) , n, MPI_DOUBLE_COMPLEX, who, tag,comm,&request);
#else
    n *=2;
    MPI_Isend(reinterpret_cast<void*> (v), n, MPI_DOUBLE, who, tag,comm,&request);
    n /= 2;
#endif
     return request;
}

template<class R> 
void  WRecv(R * v,int n,int who,int tag,MPI_Comm comm)
{
    MPI_Status status;
    MPI_Recv(reinterpret_cast<void*> (v),n, MPI_TYPE<R>::TYPE() , who, tag,comm,&status);
}

template<> 
void  WRecv<Complex> (Complex * v,int n,int who,int tag,MPI_Comm comm)
{
    MPI_Status status;
#ifdef MPI_DOUBLE_COMPLEX_
    MPI_Recv(reinterpret_cast<void*> (v), n, MPI_DOUBLE_COMPLEX, who, tag,comm,&status);
#else
    n *=2;
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
   MPIrank(int i=0,MPI_Comm com=MPI_COMM_WORLD) : who(i) , comm(com){ffassert(i>=0 && i < mpisize);} 
    
    
   
   const MPIrank & operator<<(double a)const  {
     WSend(&a, 1, who, 4,comm);
      return *this;
   }
   const MPIrank & Bcast(double & a) const {
      WBcast(&a, 1, who,comm);
      return *this;
   }
   const MPIrank & operator>>(double & a) const {
     WRecv(&a, 1, who, 4,comm);
      return *this;
   }
   const MPIrank & operator<<(long a) const {
     WSend(&a, 1, who, 5,comm);
      return *this;
   }
   const MPIrank & Bcast(long & a) const {
     WBcast(&a, 1,  who,comm);
      return *this;
   }
   const MPIrank & operator>>(long & a) const {
      WRecv(&a, 1,  who, 5,comm);
      return *this;
   }
   

    template<class R>
    const MPIrank & operator>>(KN<R> & a) const {
	assert(&a);
	if(verbosity>9)
	  cout << " ---- " << who  << "  >> " << & a << " " << a.N() << " " << MPI_TAG<KN<R> >::TAG 
	       <<" from " << mpirank << "  "<<  (R *) a << endl;
	int n= a.N();
	WRecv((R *) a, n, who, MPI_TAG<KN<R> >::TAG ,comm);
	if(verbosity>9)
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
	if(verbosity>9)
	  cout << " .... " << who  << "  >> " << & a << " " << a.N() << " " << MPI_TAG<KN<R> >::TAG 
	       <<" from  " << mpirank << "  "<<  (R *) a << endl;
	WSend((R *) a,n,who,MPI_TAG<KN<R> >::TAG,comm);
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
	
	int tag = MPI_TAG<Matrice_Creuse<R> >::TAG;		       
	MatriceMorse<R> *mA=a->A->toMatriceMorse();
	int ldata[4];
	ldata[0]=mA->n;
	ldata[1]=mA->m;
	ldata[2]=mA->nbcoef;
	ldata[3]=mA->symetrique;
    
	cout << " ldata " << ldata[0] << " " << ldata[1] <<" " << ldata[2] << " " <<ldata[3] << endl;
	WSend( ldata,4, who, tag,comm);
	WSend(  mA->lg,mA->n+1,  who, tag+1,comm);     
	WSend( mA->cl,mA->nbcoef,  who, tag+2,comm);     
	WSend(  mA->a,mA->nbcoef,   who, tag+3,comm);  
	delete mA;
	return *this;
    }
    template<class R>
    const MPIrank & operator>>(Matrice_Creuse<R> &  a) const 
    {
	if(verbosity>1) 
	    cout << " MPI << (Matrice_Creuse *) " << a << endl;
	
	int tag =  MPI_TAG<Matrice_Creuse<R> >::TAG;		       
	int ldata[4];	
	WRecv( ldata,4, who, tag,comm);
	MatriceMorse<R> *mA= new MatriceMorse<R>(ldata[0],ldata[1],ldata[2],ldata[3]); 
	WRecv(  mA->lg,mA->n+1,   who, tag+1,comm);     
	WRecv(  mA->cl,mA->nbcoef,  who, tag+2,comm);     
	WRecv(  mA->a,mA->nbcoef,  who, tag+3,comm);  
	a.A.master(mA);
	return *this;
    }
    
   const MPIrank & operator<<(Fem2D::Mesh *  a) const {
     if(verbosity>1) 
     cout << " MPI << (mesh *) " << a << endl;
      ffassert(a);
      Serialize  buf=(*a).serialize();       
       buf.mpisend(*this,MPI_TAG<Mesh *>::TAG,comm);
      return *this;
   }
    const MPIrank & operator<<(Fem2D::Mesh3 *  a) const {
	if(verbosity>1) 
	    cout << " MPI << (mesh3 *) " << a << endl;
	ffassert(a);
	Serialize  buf=(*a).serialize();       
	buf.mpisend(*this,MPI_TAG<Mesh3 *>::TAG,comm);
	return *this;
    }
    
   const MPIrank & operator>>(Fem2D::Mesh *& a) const {
     if(verbosity>1) 
     cout << " MPI >> (mesh *) &" << a << endl;
       Serialize buf(*this,Fem2D::Mesh::magicmesh,MPI_TAG<Mesh *>::TAG,comm);
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
	Serialize buf(*this,Fem2D::GenericMesh_magicmesh,MPI_TAG<Mesh3 *>::TAG,comm);
	if (a) (*a).decrement();
	a= new Fem2D::Mesh3(buf);
	//  add 3 line FH 08/12/2003  forget build quadtree sorry      
	a->BuildGTree();
	return *this;
    }
    
  operator int () const { return who;}     
};

void Serialize::mpisend(const MPIrank & rank,long tag,void * vcomm)
{
    MPI_Comm comm=static_cast<MPI_Comm> (vcomm);
      char * pp = p-sizeof(long);
      long countsave=count(); // save count 
      count()=lg; // store length in count 
      int l=lg+sizeof(long);
      if(verbosity>1) 
         cout << " -- send from  " << mpirank << " to " << rank << " serialized " << what 
              <<   ", l=" << l << ", tag=" << tag << endl;
      if (l <=sizempibuf)
        WSend(pp,l, rank, tag,comm);
      else {
         WSend(pp,sizempibuf,  rank, tag,comm);
         WSend(pp+sizempibuf,l-sizempibuf, rank, tag+1,comm);
      }
      if(verbosity>1) 
         cout << "    ok send is arrived " << endl;      
      count()=countsave; // restore count 
}

Serialize::Serialize(const MPIrank & rank,const char * wht,long tag,void * vcomm)
 :what(wht) 
{
    MPI_Comm comm=static_cast<MPI_Comm> (vcomm);

      if(verbosity>1) 
         cout << " -- waiting " << mpirank << " from  " << rank << " serialized " << what 
              << " tag = " << tag <<  endl;
   char * buf= new char [sizempibuf];
   WRecv(buf, sizempibuf,  rank, tag,comm);
   lg = * (long *) (void *) buf;
   int l=lg+sizeof(long);
   char * pp= new char[l]  ;
   if ( l <= sizempibuf) 
      memcpy(pp,buf,l);
   else 
      {
        memcpy(pp,buf,sizempibuf);
       WRecv(pp+sizempibuf,l-sizempibuf,  rank, tag+1,comm)  ;       
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

MPIrank mpiwho(long i) { return MPIrank(i);}
MPIrank mpiwho(long i,MPI_Comm comm) { return MPIrank(i,comm);}


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
//  Hack to Bypass a bug in freefem FH  ... 
template<> 
class ForEachType<MPI_Group>:  public basicForEachType{public:// correction july 2009..... FH  Hoooo....  (Il y a un bug DUR DUR FH  ...) 
    ForEachType(Function1 iv=0,Function1 id=0,Function1 OOnReturn=0):basicForEachType(typeid(MPI_Group),sizeof(MPI_Group),0,0,iv,id,OOnReturn) {
	//T i= 0.0;
    }
};

template<> 
class ForEachType<MPI_Comm>:  public basicForEachType{public:// coorection july 2009..... FH  Hoooo....  (Il y a un bug DUR DUR FH  ...) 
    ForEachType(Function1 iv=0,Function1 id=0,Function1 OOnReturn=0):basicForEachType(typeid(MPI_Comm),sizeof(MPI_Comm),0,0,iv,id,OOnReturn) {
	//T i= 0.0;
    }
};
// end Hack  ... 
MPI_Group* def_group( MPI_Group* const & a, KN_<long>  const & b)
{
    cout << b.N() <<endl;
    for(int i=0;i<b.N();++i)
	cout << b[i] << endl;
    // ici def a .. 
    ffassert(0); //   A AFAIRE  //  pour arete le programm 
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
    
 
     Global.New("mpirank",CPValue<long>(mpirank));
     Global.New("mpisize",CPValue<long>(mpisize));
    
  }
  void end_parallele()
   {
    MPI_Finalize();
    if(mpirank==0) cout << "FreeFem++-mpi finalize correctly .\n" << flush ; 
    else if(verbosity>5)  cout << '.' << endl ;
   }
//   MPI::COMM_WORLD.Recv(&msg, 1, MPI::INT, from, MPI::ANY_TAG);
//    MPI::COMM_WORLD.Isend(&msg, 1, MPI::INT, to, 4);
