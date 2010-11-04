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
// FH: I have move AFunction_ext.hpp  to fflib dir.
#include "AFunction_ext.hpp" 
// Add J. Morice for AFunction_ext.hpp 
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


#undef MPICH_SKIP_MPICXX
#define  MPICH_SKIP_MPICXX
#undef MPICH_IGNORE_CXX_SEEK
#define MPICH_IGNORE_CXX_SEEK
#include "mpi.h"

// Remark on mipich  MPI_Comm, MPI_Resquest, MPI_Group, MPI_Op are int 
//  => encapsulation


template<class MPI_type,int DIFF>
struct fMPI { 
  MPI_type v; 
  operator  MPI_type &() {return v;}
  operator  MPI_type *() {return &v;}
  operator  MPI_type () const  {return v;}
  
  // MPI_type * operator  &() {return &v;}
  void operator=(const MPI_type vv) { v=vv;}
  fMPI(const MPI_type vv=0) : v(vv){}
  bool operator!=(MPI_type vv) const {return vv !=v;}
  bool operator==(MPI_type vv) const {return vv ==v;}
  
};

// the encapsulation for the for MPI type  (int on mpich )

typedef fMPI<MPI_Comm,1> fMPI_Comm;
typedef fMPI<MPI_Group,2> fMPI_Group;
typedef fMPI<MPI_Request,3> fMPI_Request;
typedef fMPI<MPI_Op,4> fMPI_Op;



// end of encapsulation ..

// to send a sparse matrix we send header, line array ,colmun array, value array.
//  end afer the fist resquest we need to do allocation.
// so in this cas the communacatio is done in
// 2 step  
//    1 the header, 
//       alloc time 
//    2 the  message 
//  a couple request, pointer. 
//    Not use of IPROBE because probelem of wait. 
class MPIrank; 
class DoOnWaitMPI_Request ; 

map<MPI_Request*,DoOnWaitMPI_Request *> ToDoOnWaitMPI_Request;

void GetPendingWait() ;


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
template<> struct MPI_WHAT<KN<long> *>   {static const int WHAT=104;};
template<> struct MPI_WHAT<KN<double>* >   {static const int WHAT=105;};
template<> struct MPI_WHAT<KN<Complex>* >   {static const int WHAT=106;};

template<class T> struct MPI_TAG {};
template<> struct MPI_TAG<long>   {static const int TAG=5;};
template<> struct MPI_TAG<double>    {static const int TAG=4;};
template<> struct MPI_TAG<Complex >   {static const int TAG=6;};
template<> struct MPI_TAG<KN<long>* >   {static const int TAG=11;};
template<> struct MPI_TAG<KN<double>* >   {static const int TAG=12;};
template<> struct MPI_TAG<KN<Complex>* >   {static const int TAG=13;};
template<> struct MPI_TAG<Mesh *>   {static const int TAG=1000;};
template<> struct MPI_TAG<Mesh3 *>   {static const int TAG=1010;};
template<> struct MPI_TAG<Matrice_Creuse<double> *>   {static const int TAG=1020;};
template<> struct MPI_TAG<Matrice_Creuse<Complex> *>   {static const int TAG=1030;};

void initparallele(int &, char **&);
void init_lgparallele();  

extern long mpirank ;
extern long mpisize ;

//  for syncro communication
MPI_Request *  Syncro_block = reinterpret_cast<MPI_Request * > (1); 

const size_t sizempibuf = 1024*32;

template<class R> 
long  WSend( R * v,int l,int who,int tag,MPI_Comm comm,MPI_Request *rq)
{
  long ret=0;
  MPI_Request rq0,*request=&rq0;
  if(verbosity>100)
    cout << mpirank<< " send to " << who << " tag " << tag << " " << rq << " " <<  comm << " syncro "<<  (rq == Syncro_block) <<endl;
  if(rq == Syncro_block) 
    ret=MPI_Send((void *) v,l, MPI_TYPE<R>::TYPE() , who, tag,comm);
  else
    {
      ret=MPI_Isend((void *) v,l, MPI_TYPE<R>::TYPE() , who, tag,comm,request);
      if(rq) *rq=*request;
      else MPI_Request_free(request); 
    }
}

template<> 
long  WSend<Complex> ( Complex * v,int n,int who,int tag,MPI_Comm comm,MPI_Request *rq)
{
  long ret=0;
  MPI_Request rq0,*request=&rq0;
  if(verbosity>100)
    cout << mpirank<< " send to " << who << " tag " << tag << " " << rq << " " <<  comm << " syncro "<<  (rq == Syncro_block) << endl;
  if(rq == Syncro_block) 
    {
#ifdef MPI_DOUBLE_COMPLEX
      ret=MPI_Send(reinterpret_cast<void*> (v) , n, MPI_DOUBLE_COMPLEX, who, tag,comm);
#else
      n *=2;
      ret=  MPI_Send(reinterpret_cast<void*> (v), n, MPI_DOUBLE, who, tag,comm);
#endif	  
    }
  else
    {
#ifdef MPI_DOUBLE_COMPLEX
      ret=MPI_Isend(reinterpret_cast<void*> (v) , n, MPI_DOUBLE_COMPLEX, who, tag,comm,request);
#else
      n *=2;
      ret=MPI_Isend(reinterpret_cast<void*> (v), n, MPI_DOUBLE, who, tag,comm,request);
      n /= 2;
#endif
      if(rq) *rq=*request;
      else MPI_Request_free(request);
      } 
  return ret;
}

template<class R> 
long  WRecv(R * v,int n,int who,int tag,MPI_Comm comm,MPI_Request *rq)
{
  MPI_Status status;
  if(rq && (rq != Syncro_block)) 
    return MPI_Irecv(reinterpret_cast<void*> (v),n, MPI_TYPE<R>::TYPE() , who, tag,comm,rq);
  else 
    return MPI_Recv(reinterpret_cast<void*> (v),n, MPI_TYPE<R>::TYPE() , who, tag,comm,&status);
}

template<> 
long  WRecv<Complex> (Complex * v,int n,int who,int tag,MPI_Comm comm,MPI_Request *rq)
{
  MPI_Status status;
#ifdef MPI_DOUBLE_COMPLEX_
  if(rq && (rq != Syncro_block)) 
    return MPI_Irecv(reinterpret_cast<void*> (v), n, MPI_DOUBLE_COMPLEX, who, tag,comm,rq);
  else 
    return MPI_Recv(reinterpret_cast<void*> (v), n, MPI_DOUBLE_COMPLEX, who, tag,comm,&status);
#else
  n *=2;
  if(rq && (rq != Syncro_block)) 
    return  MPI_Irecv(reinterpret_cast<void*> (v), n, MPI_DOUBLE, who, tag,comm,rq);
  else 
       return MPI_Recv(reinterpret_cast<void*> (v), n, MPI_DOUBLE, who, tag,comm,&status);
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
#endif
}

			 



struct MPIrank {
  
  int who; 
  MPI_Comm comm; 
  MPI_Request *rq; 
  // mutable bool block;  
  
  MPIrank(int i=0,MPI_Comm com=MPI_COMM_WORLD, MPI_Request *rqq=0) 
    : who(i) , comm(com),rq(rqq) 
  {
	int n;
	MPI_Comm_size(comm, &n);
	// cout <<" who = " << who << " ******** " << n << " "<< ((-2 < who) && (who < n))  << endl;
	// ffassert( (-2 < who) && (who < n) ); // plant sans raison ...
  } 
  
  
  
  long Send(double a)const  {return WSend(&a, 1, who, MPI_TAG< double >::TAG,comm,rq); }
  long Send(long a)const  {return WSend(&a, 1, who, MPI_TAG< long >::TAG,comm,rq); }
  long Send(Complex a)const  {return WSend(&a, 1, who, MPI_TAG< Complex >::TAG,comm,rq); }
  
  long Recv(double & a) const { return  WRecv(&a, 1,  who, MPI_TAG< double >::TAG ,comm,rq);}
  long Recv(long & a) const { return  WRecv(&a, 1,  who, MPI_TAG< long >::TAG ,comm,rq);}
  long Recv(Complex & a) const { return  WRecv(&a, 1,  who, MPI_TAG< Complex >::TAG ,comm,rq);}
  
  const MPIrank & Bcast(double & a) const {WBcast(&a, 1, who,comm); return *this; }
  const MPIrank & Bcast(long & a) const {WBcast(&a, 1, who,comm); return *this; }
  const MPIrank & Bcast(Complex & a) const {WBcast(&a, 1, who,comm); return *this; }
  
  
  
  template<class R>
  long Recv(KN<R> & a) const {
    assert(&a);
    if(verbosity>99)
      cout << " ---- " << who  << "  >> " << & a << " " << a.N() << " " << MPI_TAG<KN<R>* >::TAG 
	       <<" from " << mpirank << "  "<<  (R *) a << endl;
    int n= a.N();
    long ll=WRecv((R *) a, n, who, MPI_TAG<KN<R>* >::TAG ,comm,rq);
    if(verbosity>99)
      cout << " ++++ " << who  << "  >> " << & a << " " << a.N() << " " << MPI_TAG<KN<R>* >::TAG 
	   <<" from  " << mpirank << "  "<<  (R *) a << endl;
    ffassert(a.N()==n);
    return ll;
  }
  
  template<class R>
  long Send(const KN<R> *aa)const  {
    const KN<R> & a=*aa;
    ffassert(&a); 
    int n= a.N();
    if(verbosity>99)
	  cout << " .... " << who  << "  >> " << & a << " " << a.N() << " " << MPI_TAG<KN<R>* >::TAG 
	       <<" from  " << mpirank << "  "<<  (R *) a << endl;
    return WSend((R *) a,n,who,MPI_TAG<KN<R>* >::TAG,comm,rq);
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
    long nbsize=0;
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
    long  nbsize=0;
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
  
  // version asyncrone or syncrone  Now 2010 ...
  template<class R>   long Send(Matrice_Creuse<R> * const &  a) const ;
  template<class R>   long Recv(Matrice_Creuse<R>  &  a) const ;
  long Send(Fem2D::Mesh *  a) const ;
  long Send (Fem2D::Mesh3 *  a) const ;
  long Recv(Fem2D::Mesh *& a)  const;  
  long Recv(Fem2D::Mesh3 *& a) const; 

  operator int () const { return who;}     
};




// for MPI_WAIT_resquets (complex MPI asyncrone MPI recv request ) ..
class DoOnWaitMPI_Request :public   MPIrank
{  
    
public:
    bool sync;
    DoOnWaitMPI_Request( MPIrank  mpr) : MPIrank(mpr),sync((rq==0 || rq == Syncro_block)) {}
    virtual  bool Do(MPI_Request *rrq) =0; // false -> end 
    bool  DoSR() { // do the  Send/Recv Op. 
	bool ret=false;
	if(verbosity>100)
	  cout << mpirank << "   --- Do Send/Recv :  "  << " " << rq << " " << sync <<  endl;
	if(sync) //  wait ...
	  { bool c=1; 
	    if(verbosity>100)
	      cout << mpirank << "   --- Do way :  " << c << " " << rq << endl;
	    while (c)
	      {
		c=Do(rq); 
		if(verbosity>100)
		  cout << mpirank << "   --- Do return :  " << c << " " << rq << endl;
	      }
	    
	    ret=true;// clean 
	  } 
	else 
	  ToDoOnWaitMPI_Request[rq]=this; // add request for WAIT ..
	return ret;
    }
    virtual ~DoOnWaitMPI_Request(){}
private:
    DoOnWaitMPI_Request(const DoOnWaitMPI_Request & );
     DoOnWaitMPI_Request & operator=( DoOnWaitMPI_Request & );
};



void DoOnWaitMPIRequest(MPI_Request *rq)
{
  if( rq  )
    { 
	map<MPI_Request*,DoOnWaitMPI_Request *>:: iterator drd = ToDoOnWaitMPI_Request.find(rq) ;
	if(drd != ToDoOnWaitMPI_Request.end())
	  {
	    if(verbosity>100)
	      cout << " Do on DoOnWaitMPIRequest " << rq  << " "  << endl; 
	    if( !drd->second->Do(rq) )
	      {
		delete drd->second;
		ToDoOnWaitMPI_Request.erase(drd); // finish ... 
	      }

	  }
	
    }

}

void DeSerialize(Serialize * sTh,Fem2D::Mesh ** ppTh)
{
      if (*ppTh) (**ppTh).decrement();
      Fem2D::Mesh * pTh= new Fem2D::Mesh(*sTh);
      *ppTh=pTh;
      Fem2D::R2 Pn,Px;
      pTh->BoundingBox(Pn,Px);
      pTh->quadtree=new Fem2D::FQuadTree(pTh,Pn,Px,pTh->nv);
}

void DeSerialize(Serialize * sTh,Fem2D::Mesh3 ** ppTh)
{
      if (*ppTh) (**ppTh).decrement();
      Fem2D::Mesh3 * pTh= new Fem2D::Mesh3(*sTh);
      pTh->BuildGTree();
      *ppTh=pTh;
}


template<class R>
class RevcWMatd : public DoOnWaitMPI_Request
{
public:  
  typedef Matrice_Creuse<R> Mat;
  Matrice_Creuse<R> * pmat;
  MatriceMorse<R> *mA;
  int state;
  int ldata[4];
  RevcWMatd(const MPIrank *mpirank,Mat * pm)
    : DoOnWaitMPI_Request(*mpirank),
      pmat(pm),mA(0),state(0)
  {
    int tag = MPI_TAG<Matrice_Creuse<R>* >::TAG;
    int ll=WRecv( ldata,4, who, tag,comm,rq);
    ffassert(ll == MPI_SUCCESS);

  }
  
  bool  Do(MPI_Request *rrq)
  {
    state++;
    int tag=MPI_TAG<Mat *>::TAG;
    if(verbosity>100)
      cout << mpirank << "  ---R: ldata " << ldata[0] << " " << ldata[1] <<" " << ldata[2] << " " <<ldata[3] << " " << state << endl;
    
    int ll=0;
    switch (state)
      {
      case 1:
	mA =  new MatriceMorse<R>(ldata[0],ldata[1],ldata[2],ldata[3]); 
	ll=WRecv(  mA->lg,mA->n+1,   who, tag+1,comm,rq);
	break;
      case 2:
	ll=WRecv(  mA->cl,mA->nbcoef,  who, tag+2,comm,rq);
	break;
      case 3:
	ll=WRecv(  mA->a,mA->nbcoef,  who, tag+3,comm,rq);
	break;
      default:
	pmat->A.master(mA);
	mA=0;
	return false;
	break;
      }
    ffassert(ll == MPI_SUCCESS);
    return true; // OK 
  }
  ~RevcWMatd() {
    if(mA) delete mA; 
  }	
  
};

template<class R>
class SendWMatd : public DoOnWaitMPI_Request
{
public:  
  typedef Matrice_Creuse<R> Mat;
  Matrice_Creuse<R> * pmat;
  MatriceMorse<R> *mA;
  int state;
  int ldata[4];
  SendWMatd(const MPIrank *mpirank,Mat * pm)
    : DoOnWaitMPI_Request(*mpirank),
      pmat(pm),mA(0),state(0)
  {
    mA=pmat->A->toMatriceMorse();
    ldata[0]=mA->n;
    ldata[1]=mA->m;
    ldata[2]=mA->nbcoef;
    ldata[3]=mA->symetrique;
    int tag = MPI_TAG<Matrice_Creuse<R>* >::TAG;
    int ll=WSend( ldata,4, who, tag,comm,rq);
    ffassert(ll == MPI_SUCCESS) ;
  }
  bool  Do(MPI_Request *rrq)
  {
    state++;
    int tag=MPI_TAG<Mat *>::TAG;
    if(verbosity>100)
      cout << mpirank << "  ---S  ldata " << ldata[0] << " " << ldata[1] <<" " << ldata[2] << " " <<ldata[3] << endl;
    
    int ll=0;
    switch (state)
      {
      case 1:
	ll=WSend(  mA->lg,mA->n+1,  who, tag+1,comm,rq);     
	break;
      case 2:
	ll=WSend( mA->cl,mA->nbcoef,  who, tag+2,comm,rq);  
	break;
      case 3:
	ll=WSend(  mA->a,mA->nbcoef,   who, tag+3,comm,rq);  
	break;
      default:
	delete mA;
	mA=0;
	return false;
	break;
      }
    ffassert(ll == MPI_SUCCESS);
    return true; // OK 
  }
  ~SendWMatd() {
    if(mA) delete mA; 
  }	
  
};


template<class Mesh>
class RevcWMeshd : public DoOnWaitMPI_Request,Serialize  
{
public:  
  Mesh ** ppTh;
  int state;
  RevcWMeshd(const MPIrank *mpirank,Mesh ** ppThh)
    : DoOnWaitMPI_Request(*mpirank),Serialize(sizempibuf,Fem2D::Mesh::magicmesh),
      ppTh(ppThh),state(0)
  {
    int tag=MPI_TAG<Mesh *>::TAG;
    if(verbosity>100)
      cout << " -- RevcWMeshd   " << rq << " " << comm << " " << p << endl; 
    char * pp = p-sizeof(long);
    int ll=WRecv(pp, sizempibuf,  who, tag,comm,rq); // wait first part ..
    // cout << mpirank << " ++ ll= " << ll << " pp= " << pp << endl;
  }
  
  bool  Do(MPI_Request *rrq)
  {
    int tag=MPI_TAG<Mesh *>::TAG;
    ffassert(rq == rrq);
    long l = * (long *) (void *) p ;
    long l1 = l -( sizempibuf-sizeof(long));
    if(verbosity>100)
      cout << mpirank << " Do RevcWMeshd " <<  l  <<" " << state << "  cont  : " <<  (l1 >0)  << " " << rq << " " << comm << endl; 
    
    if(0==state++ &&  l1>0 ) // recv first part ..
      {
	if(verbosity>100)
	  cout << mpirank << " + Do RevcWMeshd " <<  l  <<" " << state << "  cont  : " <<  ( l > sizempibuf) <<  " " << rq << " " << l-sizempibuf << " p = " << (void *) p <<  endl; 
	resize(l);
	int ll=WRecv(p-sizeof(long)+sizempibuf,l1,  who, tag+state,comm,rq);
	return true;// continue .. 	
      }
    else resize(l);
    // we have the all buffer => DeSerialize
    DeSerialize(this,ppTh);      
    count()=0; 
    if(verbosity>100) 
      cout << "    " << mpirank << " recived from " << who << " serialized " << what <<   ", l=" 
	   << l << ", tag=" << tag << " rq = " << rq << " "  << *ppTh << endl;
    
    return false; // OK 
  }
  ~RevcWMeshd() {count()=0;}	
  
};

template<class Mesh>
class SendWMeshd : public DoOnWaitMPI_Request,Serialize  
{
public:  
  Mesh ** ppTh;
  int state;
  SendWMeshd(const MPIrank *mpirank,Mesh ** ppThh)
    : DoOnWaitMPI_Request(*mpirank),Serialize((**ppThh).serialize()),
      ppTh(ppThh),state(0)
  {
    int tag=MPI_TAG<Mesh *>::TAG;
    if(verbosity>100)
      cout << " -- SendWMeshd   " << rq << " " << comm << " " << p << endl; 
    char * pp = p-sizeof(long);
    count()=lg; // store length in count 
    size_t ls=lg+sizeof(long);
    if (ls<=sizempibuf)
      WSend(pp,ls, who, tag,comm,rq);
    else 
      WSend(pp,sizempibuf,who, tag,comm,rq);
  }
  
  bool  Do(MPI_Request *rrq)
  {
    int tag=MPI_TAG<Mesh *>::TAG;
    char * pp = p-sizeof(long);
    long l1 = lg -(sizempibuf- sizeof(long));
    if(verbosity>100)
      cout << mpirank << " Do SendWMeshd " <<  lg  <<" " << state << "  cont  : " <<  (l1 >0)  << " " << rq << " " << comm << endl; 
    
    if(0==state++ &&  l1>0 ) // send the second part 
      {
	int ll=WSend(pp+sizempibuf,l1,  who, tag+state,comm,rq);
	return true;// Fini
      }
    return false; // OK 
  }

  ~SendWMeshd() {count()=0;}
  
};


template<class R>
  long MPIrank::Send(Matrice_Creuse<R> * const &  a) const 
  {
    if(0)
      {
	if(verbosity>100) 
	  cout << " MPI << (Matrice_Creuse *) " << a << endl;
	ffassert(rq==0 || rq == Syncro_block) ; // 
	int tag = MPI_TAG<Matrice_Creuse<R>* >::TAG;		       
	MatriceMorse<R> *mA=a->A->toMatriceMorse();
	int ldata[4];
	ldata[0]=mA->n;
	ldata[1]=mA->m;
	ldata[2]=mA->nbcoef;
	ldata[3]=mA->symetrique;
	
	if(verbosity>100)
	  cout << " ldata " << ldata[0] << " " << ldata[1] <<" " << ldata[2] << " " <<ldata[3] << endl;
	int ll=0;
	ll=WSend( ldata,4, who, tag,comm,rq);
	if(ll == MPI_SUCCESS) 
	  ll=WSend(  mA->lg,mA->n+1,  who, tag+1,comm,rq);     
	if(ll == MPI_SUCCESS) 
	  ll=WSend( mA->cl,mA->nbcoef,  who, tag+2,comm,rq);  
	if(ll == MPI_SUCCESS) 
	  ll=WSend(  mA->a,mA->nbcoef,   who, tag+3,comm,rq);  
	delete mA;
	return ll;
      }
    else
      {
	SendWMatd<R> *rwm= new SendWMatd<R>(this,a);
	if( rwm->DoSR() ) delete rwm;
	return MPI_SUCCESS;
      }
  }

  template<class R>
  long MPIrank::Recv(Matrice_Creuse<R>  &  a) const 
  {
    if(0)
      {
	if(verbosity>100) 
	  cout << " MPI << (Matrice_Creuse ) " << a << endl;
	ffassert(rq==0 || rq == Syncro_block) ; // 
	int tag =  MPI_TAG<Matrice_Creuse<R>* >::TAG;		       
	int ldata[4];	
	int ll=0;
	ll=WRecv( ldata,4, who, tag,comm,rq);
	MatriceMorse<R> *mA= new MatriceMorse<R>(ldata[0],ldata[1],ldata[2],ldata[3]); 
	if(ll == MPI_SUCCESS) 
	  ll=WRecv(  mA->lg,mA->n+1,   who, tag+1,comm,rq);     
	if(ll == MPI_SUCCESS) 
	  ll=WRecv(  mA->cl,mA->nbcoef,  who, tag+2,comm,rq);     
	if(ll == MPI_SUCCESS) 
	  ll=WRecv(  mA->a,mA->nbcoef,  who, tag+3,comm,rq);  
	a.A.master(mA);
	return ll;
      }
    else
      {
	RevcWMatd<R> *rwm= new RevcWMatd<R>(this,&a);
        if( rwm->DoSR() ) delete rwm;
        return MPI_SUCCESS;	
      }
  }


long MPIrank::Send(Fem2D::Mesh *  a) const {
    if(verbosity>100) 
      cout << " MPI << (mesh *) " << a << endl;
    ffassert(a);
    SendWMeshd<Mesh> *rwm= new SendWMeshd<Mesh>(this,&a);
    cout << " ... "<< endl;
    if( rwm->DoSR() ) delete rwm;
    return MPI_SUCCESS;
  }
long MPIrank::Send (Fem2D::Mesh3 *  a) const {
    if(verbosity>100) 
      cout << " MPI << (mesh3 *) " << a << endl;
    ffassert(a);
    SendWMeshd<Mesh3> *rwm= new SendWMeshd<Mesh3>(this,&a);
    if( rwm->DoSR() ) delete rwm;
    return MPI_SUCCESS;
  }

/*
long MPIrank::Send(Fem2D::Mesh *  a) const {
    if(verbosity>100) 
      cout << " MPI << (mesh *) " << a << endl;
    ffassert(a);
    Serialize  buf=(*a).serialize();       
    buf.mpisend(*this,MPI_TAG<Mesh *>::TAG,static_cast<const void *>(this));
    return MPI_SUCCESS;
  }
long MPIrank::Send (Fem2D::Mesh3 *  a) const {
    if(verbosity>100) 
      cout << " MPI << (mesh3 *) " << a << endl;
    ffassert(a);
    Serialize  buf=(*a).serialize();       
    buf.mpisend(*this,MPI_TAG<Mesh3 *>::TAG,static_cast<const void *>(this));
    return MPI_SUCCESS;
  }
*/

// new version asyncrone ...  Now 2010 ... 
long MPIrank::Recv(Fem2D::Mesh *& a) const  {
    if(verbosity>100) 
	cout << " MPI >> (mesh *) &" << a << " " << &a << endl;
    RevcWMeshd<Mesh> *rwm= new RevcWMeshd<Mesh>(this,&a);
    if( rwm->DoSR() ) delete rwm;
    if((rq==0 || rq == Syncro_block))  
      ffassert( a );
    return MPI_SUCCESS;
}

long MPIrank::Recv(Fem2D::Mesh3 *& a) const  {
    if(verbosity>100) 
      cout << " MPI >> (mesh3 *) &" << a << " " << &a << endl;
    RevcWMeshd<Mesh3> *rwm= new RevcWMeshd<Mesh3>(this,&a);
    if( rwm->DoSR() ) delete rwm;
    if((rq==0 || rq == Syncro_block))  
      ffassert( a );
    return MPI_SUCCESS;
}


void Serialize::mpisend(const MPIrank & rank,long tag,const void * vmpirank)
{
  const MPIrank * mpirank=static_cast<const MPIrank *> (vmpirank);
  MPI_Comm comm=mpirank->comm;
  MPI_Request *rq=mpirank->rq;
  ffassert(rq==0 || rq == Syncro_block);
  char * pp = p-sizeof(long);
  long countsave=count(); // save count 
  count()=lg; // store length in count 
  int l=lg+sizeof(long);
  if(verbosity>100) 
    cout << " -- send from  " << mpirank << " to " << rank << " serialized " << what 
	 <<   ", l=" << l << ", tag=" << tag << " " << (l < sizempibuf) << endl;
  if (l <=sizempibuf)
    WSend(pp,l, rank, tag,comm,rq);
  else {
    WSend(pp,sizempibuf,  rank, tag,comm,rq);
    WSend(pp+sizempibuf,l-sizempibuf, rank, tag+1,comm,rq);
  }
  if(verbosity>100) 
    cout << "    ok send is arrived " << endl;      
  count()=countsave; // restore count 
}


Serialize::Serialize(const MPIrank & rank,const char * wht,long tag,const void * vmpirank)
  :what(wht) 
{
  const MPIrank * mpirank=static_cast<const MPIrank *> (vmpirank);
  MPI_Comm comm=mpirank->comm;
  MPI_Request *rq=mpirank->rq;
  
  if(verbosity>100) 
    cout << " -- waiting " << mpirank << " from  " << rank << " serialized " << what 
	 << " tag = " << tag <<  endl;
  if(!(rq==0 || rq == Syncro_block))
    {
      ExecError("Not async recv of complex  objet!  Sorry to hard to code (FH!).");
      ffassert(rq==0 || rq == Syncro_block); 
    }
      
 
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
  
  if(verbosity>100) 
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
    f.Recv(*a);
    return f;
  }
};

template<class A>
struct Op_Recvmpi : public binary_function<MPIrank,A*,long> {
  static MPIrank  f(MPIrank const  & f,A *  const  & a)  
  { 
    ffassert(f.rq ==0 || f.rq == Syncro_block); // Block 
    return f.Recv(*a);
    
  }
};
template<class A>
struct Op_IRecvmpi : public binary_function<MPIrank,A*,long> {
  static MPIrank  f(MPIrank const  & f,A *  const  & a)  
  { 
    ffassert(f.rq !=0 || f.rq != Syncro_block); // no Block 
    return f.Recv(*a);
    
  }
};


template<class A>
struct Op_Writempi : public binary_function<MPIrank,A,MPIrank> {
  static MPIrank  f(MPIrank const  & f,A   const  &  a)  
  { 
    f.Send(a);
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

template<class A>
struct Op_ISendmpi : public binary_function<MPIrank,A,long> {
  static MPIrank  f(MPIrank const  & f,A   const  & a)  
  { 
    ffassert(f.rq != Syncro_block); 
    return f.Send(a);
  }
};
template<class A>
struct Op_Sendmpi : public binary_function<MPIrank,A,long> {
  static MPIrank  f(MPIrank const  & f,A  const  & a)  
  { 
	MPIrank ff(f.who,f.comm,Syncro_block); 
	return ff.Send(a);
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
struct Op_Allgather1 : public binary_function<R*,KN_<R>,long> {
  static long  f( R*  const  & s, KN_<R>  const  &r)  
    { 
      MPI_Comm comm=MPI_COMM_WORLD;
      int mpisizew;
      MPI_Comm_size(comm, &mpisizew); /* local */ 
      int chunk = 1;
      ffassert(r.N()==mpisizew);
      
      return MPI_Allgather( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
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
      int chunk = r.N()/mpisizew;
      if( ! 	     (r.N()==mpisizew*chunk && chunk==s.N()) )
	{
	  cout << " ???? Error size buf  r.N " << r.N() << " s.N = " << s.N() << " mpisizew " << mpisizew << endl; 
	  ffassert(r.N()==mpisizew*chunk && chunk==s.N());
	}
      return MPI_Allgather( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			    (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(), comm);	
    }
};

template<class R>
struct Op_All2All3 : public ternary_function<KN_<R>,KN_<R>,fMPI_Comm,long> {
  static long  f(Stack, KN_<R>  const  & s, KN_<R>  const  &r,fMPI_Comm const & cmm )  
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
struct Op_Allgather3 : public ternary_function<KN_<R>,KN_<R>,fMPI_Comm,long> {
  static long  f(Stack, KN_<R>  const  & s, KN_<R>  const  &r,fMPI_Comm const & cmm)  
  { 
    MPI_Comm comm=cmm;
    int mpisizew;
    MPI_Comm_size(comm, &mpisizew); /* local */ 
    int chunk = r.N()/mpisizew;    // bug corrected by J. Morice
    //ffassert(s.N()==mpisizew*chunk && r.N()==s.N());
    ffassert(s.N()==chunk && r.N()==s.N()*mpisizew);
	
    return MPI_Allgather( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			  (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(), comm);	
  }
};


template<class R>
struct Op_Allgather13 : public ternary_function<R*,KN_<R>,fMPI_Comm,long> {
  static long  f(Stack, R*  const  & s, KN_<R>  const  &r,fMPI_Comm const & cmm)  
  { 
    MPI_Comm comm=cmm;
    int mpisizew;
    MPI_Comm_size(comm, &mpisizew); /* local */ 
    int chunk = 1;    // bug corrected by J. Morice
    //ffassert(s.N()==mpisizew*chunk && r.N()==s.N());
    ffassert( r.N()==mpisizew);
	
    return MPI_Allgather( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			  (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(), comm);	
  }
};
// Add J. Morice

template<class R>
long  Op_All2Allv( KN_<R>  const  & s, KN_<R>  const  &r, KN_<long> const &sendcnts, KN_<long> const &sdispls, KN_<long> const &recvcnts, KN_<long> const &rdispls)  
{ 
  MPI_Comm comm=MPI_COMM_WORLD;
  int mpirankv=MPI_UNDEFINED;
  MPI_Comm_rank(comm, &mpirankv); 
  
  int mpisizew;
  MPI_Comm_size(comm, &mpisizew); /* local */ 
  ffassert( sendcnts.N() == sdispls.N() && sendcnts.N() == recvcnts.N() && sendcnts.N() == rdispls.N() && sendcnts.N() == mpisizew );
  
  KN<int> INTsendcnts(sendcnts.N());
  KN<int> INTsdispls(sdispls.N());
  KN<int> INTrecvcnts(recvcnts.N());
  KN<int> INTrdispls(rdispls.N());
  
  for(int ii=0; ii< sendcnts.N(); ii++){
    INTsendcnts[ii] = sendcnts[ii];
    INTsdispls[ii]  = sdispls[ii];
    INTrecvcnts[ii]  = recvcnts[ii];
    INTrdispls[ii]  = rdispls[ii];
  }
  
  return MPI_Alltoallv( (void *) (R*) s, INTsendcnts, INTsdispls, MPI_TYPE<R>::TYPE(),
			(void *) (R*) r, INTrecvcnts, INTrdispls, MPI_TYPE<R>::TYPE(), comm);	
}



template<class R>
struct Op_Allgatherv : public quad_function<KN_<R>,KN_<R>,KN_<long>,KN_<long>,long> {
  static long f( Stack ,KN_<R>  const  & s, KN_<R>  const  &r, KN_<long> const & recvcount, KN_<long> const & displs)  
  { 
    MPI_Comm comm=MPI_COMM_WORLD;
    int mpisizew;
    MPI_Comm_size(comm, &mpisizew); 
    ffassert( recvcount.N() == displs.N() && recvcount.N() == mpisizew);
    long sum=0;
    for(int ii=0; ii< recvcount.N(); ii++)
      sum+=recvcount[ii];
    ffassert( sum == r.N() );
    
    KN<int> INTrecvcount(recvcount.N());
    KN<int> INTdispls(displs.N());
    for(int ii=0; ii< recvcount.N(); ii++){
      INTrecvcount[ii]= recvcount[ii];
      INTdispls[ii]= displs[ii];
    }
    return MPI_Allgatherv( (void *) (R*) s, s.N(), MPI_TYPE<R>::TYPE(),
			   (void *) (R*) r, INTrecvcount, INTdispls,MPI_TYPE<R>::TYPE(), comm);	
  }
};

template<class R>
long  Op_All2All3v(KN_<R>  const  & s, KN_<R>  const  &r,fMPI_Comm const & cmm, KN_<long> const &sendcnts, KN_<long> const &sdispls, KN_<long> const &recvcnts, KN_<long> const &rdispls )  
{ 
  MPI_Comm comm=cmm;
  int mpirankv=MPI_UNDEFINED;
  MPI_Comm_rank(comm, &mpirankv); 
      
  int mpisizew;
  MPI_Comm_size(comm, &mpisizew); /* local */ 
  ffassert( sendcnts.N() == sdispls.N() && sendcnts.N() == recvcnts.N() && sendcnts.N() == rdispls.N() && sendcnts.N() == mpisizew );
      
  //ffassert(s.N()==sendcnts[mpirankv] && r.N()==recvbuf[mpirankv]);
      
  KN<int> INTsendcnts(sendcnts.N());
  KN<int> INTsdispls(sdispls.N());
  KN<int> INTrecvcnts(recvcnts.N());
  KN<int> INTrdispls(rdispls.N());
      
  for(int ii=0; ii< sendcnts.N(); ii++){
    INTsendcnts[ii] = sendcnts[ii];
    INTsdispls[ii]  = sdispls[ii];
    INTrecvcnts[ii]  = recvcnts[ii];
    INTrdispls[ii]  = rdispls[ii];
  }

  return MPI_Alltoallv( (void *) (R*) s, INTsendcnts, INTsdispls, MPI_TYPE<R>::TYPE(),
			(void *) (R*) r, INTrecvcnts, INTrdispls, MPI_TYPE<R>::TYPE(), comm);	
}


template<class R>
long Op_Allgatherv3(KN_<R>  const  & s, KN_<R>  const  &r,fMPI_Comm const & cmm, KN_<long> const & recvcount, KN_<long> const & displs)
{ 
  MPI_Comm comm=cmm;
  int mpisizew;
  MPI_Comm_size(comm, &mpisizew); 
  ffassert( recvcount.N() == displs.N() && recvcount.N() == mpisizew);
  long sum=0;
  for(int ii=0; ii< recvcount.N(); ii++)
    sum+=recvcount[ii];
  ffassert( sum == r.N() );
  KN<int> INTrecvcount(recvcount.N());
  KN<int> INTdispls(displs.N());
  for(int ii=0; ii< recvcount.N(); ii++){
    INTrecvcount[ii]= recvcount[ii];
    INTdispls[ii]= displs[ii];
  }

  return MPI_Allgatherv( (void *) (R*) s, s.N(), MPI_TYPE<R>::TYPE(),
			 (void *) (R*) r, INTrecvcount, INTdispls,MPI_TYPE<R>::TYPE(), comm);	
}


template<class R>
struct Op_Scatter1 : public   ternary_function<KN_<R>, R* ,MPIrank,long> {
  static long  f(Stack, KN_<R>  const  & s, R*  const  &r,  MPIrank const & root)  
  { 
    
    int mpisizew;
    MPI_Comm_size(root.comm, &mpisizew); 
    int chunk = 1;
    ffassert(s.N()==mpisizew*chunk);
    
    return MPI_Scatter( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			(void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(),root.who,root.comm);	
  }
};


// Fin add J. Morice


template<class R>
struct Op_Scatter3 : public   ternary_function<KN_<R>,KN_<R>,MPIrank,long> {
  static long  f(Stack, KN_<R>  const  & s, KN_<R>  const  &r,  MPIrank const & root)  
  { 
    
    int mpisizew;
    MPI_Comm_size(root.comm, &mpisizew); 
    int chunk = s.N()/mpisizew;
    ffassert(s.N()==mpisizew*chunk && r.N()==chunk);
    
    return MPI_Scatter( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			(void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(),root.who,root.comm);	
  }
};

// Add J. Morice
template<class R>
//struct Op_Scatterv3 : public   penta_function< KN_<R>, KN_<R>, MPIrank, KN_<long>, KN_<long>, long> {
long Op_Scatterv3( KN_<R>  const  & s, KN_<R>  const  &r,  MPIrank const & root, KN_<long> const &sendcnts, KN_<long> const &displs)  
{ 
  
  int mpirankv=MPI_UNDEFINED;
  if(root.comm != MPI_COMM_NULL)
    MPI_Comm_rank(root.comm, &mpirankv); 
  
  int mpisizew;
  MPI_Comm_size(root.comm, &mpisizew); /* local */ 
  ffassert( sendcnts.N() == displs.N() && sendcnts.N() == mpisizew );
  // size control 
  ffassert( r.N() == sendcnts[mpirankv] );
  long sumsize=0;
  for(int ii=0; ii<sendcnts.N(); ii++){
    sumsize += sendcnts[ii];
  }
  ffassert( s.N() == sumsize );
  
  KN<int> INTsendcnts(sendcnts.N());
  KN<int> INTdispls(displs.N());
  for(int ii=0; ii< sendcnts.N(); ii++){
    INTsendcnts[ii]= sendcnts[ii];
    INTdispls[ii]= displs[ii];
  }
  
  return MPI_Scatterv( (void *) (R*) s, INTsendcnts, INTdispls, MPI_TYPE<R>::TYPE(),
		       (void *) (R*) r, r.N(), MPI_TYPE<R>::TYPE(),root.who,root.comm);
}

// fin J. Morice

template<class R>
struct Op_Reduce  : public   quad_function<KN_<R>,KN_<R>,MPIrank,fMPI_Op,long> {
  static long  f(Stack, KN_<R>  const  & s, KN_<R>  const  &r,  MPIrank const & root, fMPI_Op const &op)  
  { 
    int chunk = s.N();
    ffassert(chunk==r.N());
    
    return MPI_Reduce( (void *) (R*) s,(void *) (R*) r, chunk , MPI_TYPE<R>::TYPE(),op,root.who,root.comm);	
  }
};

template<class R>
struct Op_AllReduce  : public   quad_function<KN_<R>,KN_<R>,fMPI_Comm,fMPI_Op,long> {
  static long  f(Stack, KN_<R>  const  & s, KN_<R>  const  &r,  fMPI_Comm const & comm,fMPI_Op const &op)  
  { 
    int chunk = s.N();
    ffassert(chunk==r.N());
    return MPI_Allreduce( (void *) (R*) s,(void *) (R*) r, chunk , MPI_TYPE<R>::TYPE(),op,comm);	
  }
};
/*
template<class R>
struct Op_Reducescatter  : public   quad_function<KN_<R>,KN_<R>,fMPI_Comm,fMPI_Op,long> {
    static long  f(Stack, KN_<R>  const  & s, KN_<R>  const  &r,  fMPI_Comm const & comm,fMPI_Op const &op)  
    { 
	int chunk = s.N();
	ffassert(chunk==r.N());
        //    chunk est un tableau ????
	MPI_Op oop = reinterpret_cast<MPI_Op> (op);
	return MPI_Reduce_scatter( (void *) (R*) s,(void *) (R*) r, chunk , MPI_TYPE<R>::TYPE(),op,comm);	
    }
};*/

template<class R>
struct Op_Reduce1  : public   quad_function<R*,R*,MPIrank,fMPI_Op,long> {
  static long  f(Stack, R*  const  & s, R*  const  &r,  MPIrank const & root, fMPI_Op const &op)  
  { 
    int chunk = 1;
    return MPI_Reduce( (void *) (R*) s,(void *) (R*) r, chunk , MPI_TYPE<R>::TYPE(),op,root.who,root.comm);	
  }
};

// Add J. Morice
template<class R>
struct Op_Gather1 : public   ternary_function<R*,KN_<R>,MPIrank,long> {
    static long  f(Stack, R* const  & s, KN_<R>  const  &r,  MPIrank const & root)  
  { 
    
    int mpisizew;
    MPI_Comm_size(root.comm, &mpisizew); 
    int chunk = 1;
    ffassert(r.N()==mpisizew*chunk );
    
    return MPI_Gather( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			   (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(),root.who,root.comm);	
  }
};

// Fin Add J. Morice


template<class R>
struct Op_Gather3 : public   ternary_function<KN_<R>,KN_<R>,MPIrank,long> {
    static long  f(Stack, KN_<R>  const  & s, KN_<R>  const  &r,  MPIrank const & root)  
  { 
    
    int mpisizew;
    MPI_Comm_size(root.comm, &mpisizew); 
    int chunk = r.N()/mpisizew;
    ffassert(r.N()==mpisizew*chunk && chunk==s.N());
    
    return MPI_Gather( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			   (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(),root.who,root.comm);	
  }
};

// Add by J. Morice

template<class R>
//struct Op_Gatherv3 : public penta_function<KN_<R>,KN_<R>, MPIrank, KN_<long>, KN_<long>, long> {
long  Op_Gatherv3(KN_<R>  const  & s, KN_<R>  const  &r,  MPIrank const & root, KN_<long> const & recvcount, KN_<long> const & displs)  
{ 
    
  int mpirankw;
  MPI_Comm_rank(root.comm, &mpirankw); 
  int mpisizew;
  MPI_Comm_size(root.comm, &mpisizew); 
  ffassert( recvcount.N() == displs.N() && recvcount.N() == mpisizew);
  
  if( mpirankw == root.who){
    long sum=0;
    for(int ii=0; ii< recvcount.N(); ii++)
      sum+=recvcount[ii];
    ffassert( sum == r.N() );
  }
  KN<int> INTrecvcount(recvcount.N());
  KN<int> INTdispls(displs.N());
  for(int ii=0; ii< recvcount.N(); ii++){
    INTrecvcount[ii]= recvcount[ii];
    INTdispls[ii]= displs[ii];
  }
  
  return MPI_Gatherv( (void *) (R*) s, s.N(), MPI_TYPE<R>::TYPE(),
		      (void *) (R*) r, INTrecvcount, INTdispls,MPI_TYPE<R>::TYPE(),root.who,root.comm);	
}

// Fin Add J. Morice


// Add J. Morice communications entre processeurs complex

template<>
struct Op_All2All<Complex> : public binary_function<KN_<Complex>,KN_<Complex>,long> {
  static long  f( KN_<Complex>  const  & s, KN_<Complex>  const  &r)  
  { 
    MPI_Comm comm=MPI_COMM_WORLD;
    int mpisizew;
    MPI_Comm_size(comm, &mpisizew); /* local */ 
    int chunk = s.N()/mpisizew;
    ffassert(s.N()==mpisizew*chunk && r.N()==s.N());
   
#ifdef MPI_DOUBLE_COMPLEX 
    return MPI_Alltoall( (void *) (Complex*) s, chunk, MPI_DOUBLE_COMPLEX,
			 (void *) (Complex*) r, chunk, MPI_DOUBLE_COMPLEX, comm);	
#else
    chunk*=2;
    return MPI_Alltoall( reinterpret_cast<void*> (s), chunk, MPI_DOUBLE,
			 reinterpret_cast<void*> (r), chunk, MPI_DOUBLE, comm);	
#endif
  }
};

template<>
struct Op_Allgather1<Complex> : public binary_function<Complex *,KN_<Complex>,long> {
  static long  f( Complex *  const  & s, KN_<Complex>  const  &r)  
  { 
    MPI_Comm comm=MPI_COMM_WORLD;
    int mpisizew;
    MPI_Comm_size(comm, &mpisizew); /* local */ 
    int chunk = 1;    
    ffassert( r.N()== mpisizew);
#ifdef MPI_DOUBLE_COMPLEX
    return MPI_Allgather( (void *) (Complex*) s, chunk, MPI_DOUBLE_COMPLEX,
			  (void *) (Complex*) r, chunk, MPI_DOUBLE_COMPLEX, comm);
#else
    chunk*2=;
    return MPI_Allgather( reinterpret_cast<void*> (s), chunk, MPI_DOUBLE,
			  reinterpret_cast<void*> (r), chunk, MPI_DOUBLE, comm);
#endif
  }
};


template<>
struct Op_Allgather<Complex> : public binary_function<KN_<Complex>,KN_<Complex>,long> {
  static long  f( KN_<Complex>  const  & s, KN_<Complex>  const  &r)  
    { 
      MPI_Comm comm=MPI_COMM_WORLD;
      int mpisizew;
      MPI_Comm_size(comm, &mpisizew); /* local */ 
      int chunk = r.N()/mpisizew;
      ffassert( r.N()==chunk*mpisizew && chunk==s.N() );
#ifdef MPI_DOUBLE_COMPLEX 
      return MPI_Allgather( (void *) (Complex*) s, chunk, MPI_DOUBLE_COMPLEX,
			    (void *) (Complex*) r, chunk, MPI_DOUBLE_COMPLEX, comm);
#else
      chunk*=2;
      return MPI_Allgather( reinterpret_cast<void*> (s), chunk, MPI_DOUBLE,
			     reinterpret_cast<void*> (r), chunk, MPI_DOUBLE, comm);
#endif
    }
};

template<>
struct Op_All2All3<Complex> : public ternary_function<KN_<Complex>,KN_<Complex>,fMPI_Comm,long> {
  static long  f(Stack, KN_<Complex>  const  & s, KN_<Complex>  const  &r,fMPI_Comm const & cmm )  
    { 
      MPI_Comm comm=cmm;
      int mpisizew;
      MPI_Comm_size(comm, &mpisizew); /* local */ 
      int chunk = s.N()/mpisizew;
      ffassert(s.N()==mpisizew*chunk && r.N()==s.N());
#ifdef MPI_DOUBLE_COMPLEX       
      return MPI_Alltoall( (void *) (Complex*) s, chunk, MPI_DOUBLE_COMPLEX,
			   (void *) (Complex*) r, chunk, MPI_DOUBLE_COMPLEX, comm);
#else
      chunk*=2;
      return MPI_Alltoall( reinterpret_cast<void*> (s), chunk, MPI_DOUBLE,
			   reinterpret_cast<void*> (r), chunk, MPI_DOUBLE, comm);
#endif	
    }
};

template<>
struct Op_Allgather3<Complex> : public ternary_function<KN_<Complex>,KN_<Complex>,fMPI_Comm,long> {
  static long  f(Stack, KN_<Complex>  const  & s, KN_<Complex>  const  &r,fMPI_Comm const & cmm)  
  { 
    MPI_Comm comm=cmm;
    int mpisizew;
    MPI_Comm_size(comm, &mpisizew); /* local */ 
    int chunk = r.N()/mpisizew;    // bug corrected by J. Morice
    //ffassert(s.N()==mpisizew*chunk && r.N()==s.N());
    ffassert(s.N()==chunk && r.N()==s.N()*mpisizew);
#ifdef MPI_DOUBLE_COMPLEX
    return MPI_Allgather( (void *) (Complex*) s, chunk, MPI_DOUBLE_COMPLEX,
			  (void *) (Complex*) r, chunk, MPI_DOUBLE_COMPLEX, comm);
#else
    chunk*2=;
    return MPI_Allgather( reinterpret_cast<void*> (s), chunk, MPI_DOUBLE,
			  reinterpret_cast<void*> (r), chunk, MPI_DOUBLE, comm);
#endif
  }
};

template<>
struct Op_Allgather13<Complex> : public ternary_function<Complex *,KN_<Complex>,fMPI_Comm,long> {
  static long  f(Stack, Complex *  const  & s, KN_<Complex>  const  &r,fMPI_Comm const & cmm)  
  { 
    MPI_Comm comm=cmm;
    int mpisizew;
    MPI_Comm_size(comm, &mpisizew); /* local */ 
    int chunk = 1;    
    ffassert( r.N()==mpisizew);
#ifdef MPI_DOUBLE_COMPLEX
    return MPI_Allgather( (void *) (Complex*) s, chunk, MPI_DOUBLE_COMPLEX,
			  (void *) (Complex*) r, chunk, MPI_DOUBLE_COMPLEX, comm);
#else
    chunk*2=;
    return MPI_Allgather( reinterpret_cast<void*> (s), chunk, MPI_DOUBLE,
			  reinterpret_cast<void*> (r), chunk, MPI_DOUBLE, comm);
#endif
  }
};


template<>
long  Op_All2Allv<Complex>( KN_<Complex>  const  & s, KN_<Complex>  const  &r, KN_<long> const &sendcnts, KN_<long> const &sdispls, KN_<long> const &recvcnts, KN_<long> const &rdispls)  
{ 
  MPI_Comm comm=MPI_COMM_WORLD;
  int mpirankv=MPI_UNDEFINED;
  MPI_Comm_rank(comm, &mpirankv); 
  
  int mpisizew;
  MPI_Comm_size(comm, &mpisizew); /* local */ 
  ffassert( sendcnts.N() == sdispls.N() && sendcnts.N() == recvcnts.N() && sendcnts.N() == rdispls.N() && sendcnts.N() == mpisizew );
  
  KN<int> INTsendcnts(sendcnts.N());
  KN<int> INTsdispls(sdispls.N());
  KN<int> INTrecvcnts(recvcnts.N());
  KN<int> INTrdispls(rdispls.N());
  
  
#ifdef MPI_DOUBLE_COMPLEX
  for(int ii=0; ii< sendcnts.N(); ii++){
    INTsendcnts[ii] = sendcnts[ii];
    INTsdispls[ii]  = sdispls[ii];
    INTrecvcnts[ii]  = recvcnts[ii];
    INTrdispls[ii]  = rdispls[ii];
  }
  return MPI_Alltoallv( (void *) (Complex*) s, INTsendcnts, INTsdispls, MPI_DOUBLE_COMPLEX,
			(void *) (Complex*) r, INTrecvcnts, INTrdispls, MPI_DOUBLE_COMPLEX, comm);	
#else
  for(int ii=0; ii< sendcnts.N(); ii++){
    INTsendcnts[ii] = 2*sendcnts[ii];
    INTsdispls[ii]  = 2*sdispls[ii];
    INTrecvcnts[ii] = 2*recvcnts[ii];
    INTrdispls[ii]  = 2*rdispls[ii];
  }
  return MPI_Alltoallv( reinterpret_cast<void*> (s), INTsendcnts, INTsdispls, MPI_DOUBLE,
			reinterpret_cast<void*> (r), INTrecvcnts, INTrdispls, MPI_DOUBLE, comm);
#endif
}



template<>
struct Op_Allgatherv<Complex> : public quad_function<KN_<Complex>,KN_<Complex>,KN_<long>,KN_<long>,long> {
  static long f( Stack ,KN_<Complex>  const  & s, KN_<Complex>  const  &r, KN_<long> const & recvcount, KN_<long> const & displs)  
  { 
    MPI_Comm comm=MPI_COMM_WORLD;
    int mpisizew;
    MPI_Comm_size(comm, &mpisizew); 
    ffassert( recvcount.N() == displs.N() && recvcount.N() == mpisizew);
    long sum=0;
    for(int ii=0; ii< recvcount.N(); ii++)
      sum+=recvcount[ii];
    ffassert( sum == r.N() );
  
  
    KN<int> INTrecvcount(recvcount.N());
    KN<int> INTdispls(displs.N());
#ifdef MPI_DOUBLE_COMPLEX
    for(int ii=0; ii< recvcount.N(); ii++){
      INTrecvcount[ii]= recvcount[ii];
      INTdispls[ii]= displs[ii];
    }
    return MPI_Allgatherv( (void *) (Complex*)s, s.N(), MPI_DOUBLE_COMPLEX,
			   (void *) (Complex*)r, INTrecvcount, INTdispls,MPI_DOUBLE_COMPLEX, comm);
#else
    for(int ii=0; ii< recvcount.N(); ii++){
      INTrecvcount[ii]= 2*recvcount[ii];
      INTdispls[ii]= 2*displs[ii];
    }
    return MPI_Allgatherv( reinterpret_cast<void*> (s), 2*s.N(), MPI_DOUBLE,
			   reinterpret_cast<void*> (r), INTrecvcount, INTdispls,MPI_DOUBLE, comm);
#endif	
  }
};

template<>
long  Op_All2All3v<Complex>(KN_<Complex>  const  & s, KN_<Complex>  const  &r,fMPI_Comm const & cmm, KN_<long> const &sendcnts, KN_<long> const &sdispls, KN_<long> const &recvcnts, KN_<long> const &rdispls )  
{ 
  MPI_Comm comm=cmm;
  int mpirankv=MPI_UNDEFINED;
  MPI_Comm_rank(comm, &mpirankv); 
      
  int mpisizew;
  MPI_Comm_size(comm, &mpisizew); /* local */ 
  ffassert( sendcnts.N() == sdispls.N() && sendcnts.N() == recvcnts.N() && sendcnts.N() == rdispls.N() && sendcnts.N() == mpisizew );
      
  //ffassert(s.N()==sendcnts[mpirankv] && r.N()==recvbuf[mpirankv]);
      
  KN<int> INTsendcnts(sendcnts.N());
  KN<int> INTsdispls(sdispls.N());
  KN<int> INTrecvcnts(recvcnts.N());
  KN<int> INTrdispls(rdispls.N());
    
#ifdef MPI_DOUBLE_COMPLEX  
  for(int ii=0; ii< sendcnts.N(); ii++){
    INTsendcnts[ii] = sendcnts[ii];
    INTsdispls[ii]  = sdispls[ii];
    INTrecvcnts[ii]  = recvcnts[ii];
    INTrdispls[ii]  = rdispls[ii];
  }

  return MPI_Alltoallv( (void *) (Complex*)s, INTsendcnts, INTsdispls, MPI_DOUBLE_COMPLEX,
			(void *) (Complex*)r, INTrecvcnts, INTrdispls, MPI_DOUBLE_COMPLEX, comm);	
#else
  for(int ii=0; ii< sendcnts.N(); ii++){
    INTsendcnts[ii] = 2*sendcnts[ii];
    INTsdispls[ii]  = 2*sdispls[ii];
    INTrecvcnts[ii]  = 2*recvcnts[ii];
    INTrdispls[ii]  = 2*rdispls[ii];
  }

  return MPI_Alltoallv( reinterpret_cast<void*> (s), INTsendcnts, INTsdispls, MPI_DOUBLE,
			reinterpret_cast<void*> (r), INTrecvcnts, INTrdispls, MPI_DOUBLE, comm);
#endif
}


template<>
long Op_Allgatherv3<Complex>(KN_<Complex>  const  & s, KN_<Complex>  const  &r,fMPI_Comm const & cmm, KN_<long> const & recvcount, KN_<long> const & displs)
{ 
  MPI_Comm comm=cmm;
  int mpisizew;
  MPI_Comm_size(comm, &mpisizew); 
  ffassert( recvcount.N() == displs.N() && recvcount.N() == mpisizew);
  long sum=0;
  for(int ii=0; ii< recvcount.N(); ii++)
    sum+=recvcount[ii];
  ffassert( sum == r.N() );
  KN<int> INTrecvcount(recvcount.N());
  KN<int> INTdispls(displs.N());
  
#ifdef MPI_DOUBLE_COMPLEX 
  for(int ii=0; ii< recvcount.N(); ii++){
    INTrecvcount[ii]= recvcount[ii];
    INTdispls[ii]= displs[ii];
  }

  return MPI_Allgatherv( (void *) (Complex*)s, s.N(), MPI_DOUBLE_COMPLEX,
			 (void *) (Complex*)r, INTrecvcount, INTdispls,MPI_DOUBLE_COMPLEX, comm);	

#else
  for(int ii=0; ii< recvcount.N(); ii++){
    INTrecvcount[ii]= 2*recvcount[ii];
    INTdispls[ii]= 2*displs[ii];
  }

  return MPI_Allgatherv( reinterpret_cast<void*> (s), 2*s.N(), MPI_DOUBLE,
			 reinterpret_cast<void*> (r), INTrecvcount, INTdispls,MPI_DOUBLE, comm);
#endif

}

template<>
struct Op_Scatter1<Complex> : public   ternary_function<KN_<Complex>,Complex *,MPIrank,long> {
  static long  f(Stack, KN_<Complex> const  & s, Complex *  const  &r,  MPIrank const & root)  
  { 
    
    int mpisizew;
    MPI_Comm_size(root.comm, &mpisizew); 
    int chunk = 1;
    ffassert(s.N()==mpisizew*chunk );
#ifdef MPI_DOUBLE_COMPLEX      
    return MPI_Scatter( (void *) (Complex*)s, chunk, MPI_DOUBLE_COMPLEX,
			(void *) (Complex*)r, chunk, MPI_DOUBLE_COMPLEX,root.who,root.comm);	
#else
    chunk*=2;
    return MPI_Scatter( reinterpret_cast<void*> (s), chunk, MPI_DOUBLE,
			reinterpret_cast<void*> (r), chunk, MPI_DOUBLE,root.who,root.comm);
#endif
  }
};


template<>
struct Op_Scatter3<Complex> : public   ternary_function<KN_<Complex>,KN_<Complex>,MPIrank,long> {
  static long  f(Stack, KN_<Complex>  const  & s, KN_<Complex>  const  &r,  MPIrank const & root)  
  { 
    
    int mpisizew;
    MPI_Comm_size(root.comm, &mpisizew); 
    int chunk = s.N()/mpisizew;
    ffassert(s.N()==mpisizew*chunk && r.N()==chunk);
#ifdef MPI_DOUBLE_COMPLEX      
    return MPI_Scatter( (void *) (Complex*)s, chunk, MPI_DOUBLE_COMPLEX,
			(void *) (Complex*)r, chunk, MPI_DOUBLE_COMPLEX,root.who,root.comm);	
#else
    chunk*=2;
    return MPI_Scatter( reinterpret_cast<void*> (s), chunk, MPI_DOUBLE,
			reinterpret_cast<void*> (r), chunk, MPI_DOUBLE,root.who,root.comm);
#endif
  }
};

template<>
long Op_Scatterv3<Complex>( KN_<Complex>  const  & s, KN_<Complex>  const  &r,  MPIrank const & root, KN_<long> const &sendcnts, KN_<long> const &displs)  
{ 
  
  int mpirankv=MPI_UNDEFINED;
  if(root.comm != MPI_COMM_NULL)
    MPI_Comm_rank(root.comm, &mpirankv); 
  
  int mpisizew;
  MPI_Comm_size(root.comm, &mpisizew); /* local */ 
  ffassert( sendcnts.N() == displs.N() && sendcnts.N() == mpisizew );
  // size control 
  ffassert( r.N() == sendcnts[mpirankv] );
  long sumsize=0;
  for(int ii=0; ii<sendcnts.N(); ii++){
    sumsize += sendcnts[ii];
  }
  ffassert( s.N() == sumsize );
  
  KN<int> INTsendcnts(sendcnts.N());
  KN<int> INTdispls(displs.N());
#ifdef MPI_DOUBLE_COMPLEX  
  for(int ii=0; ii< sendcnts.N(); ii++){
    INTsendcnts[ii]= sendcnts[ii];
    INTdispls[ii]= displs[ii];
  }

  return MPI_Scatterv( (void *) (Complex*)s, INTsendcnts, INTdispls, MPI_DOUBLE_COMPLEX,
		       (void *) (Complex*)r, r.N(), MPI_DOUBLE_COMPLEX,root.who,root.comm);
#else
  for(int ii=0; ii< sendcnts.N(); ii++){
    INTsendcnts[ii]= 2*sendcnts[ii];
    INTdispls[ii]= 2*displs[ii];
  }

  return MPI_Scatterv( reinterpret_cast<void*> (s), INTsendcnts, INTdispls, MPI_DOUBLE,
		       reinterpret_cast<void*> (r), 2*r.N(), MPI_DOUBLE,root.who,root.comm);
#endif

}


template<>
struct Op_Reduce<Complex>  : public   quad_function<KN_<Complex>,KN_<Complex>,MPIrank,fMPI_Op,long> {
  static long  f(Stack, KN_<Complex>  const  & s, KN_<Complex>  const  &r,  MPIrank const & root, fMPI_Op const &op)  
  { 
    int chunk = s.N();
    ffassert(chunk==r.N());
#ifdef MPI_DOUBLE_COMPLEX     
    return MPI_Reduce( (void *) (Complex*)s,(void *) (Complex*)r, chunk , MPI_DOUBLE_COMPLEX,op,root.who,root.comm);
#else
    chunk*=2;
    return MPI_Reduce( reinterpret_cast<void*> (s), reinterpret_cast<void*> (r), chunk , MPI_DOUBLE,op,root.who,root.comm);
#endif	
  }
};

template<>
struct Op_AllReduce<Complex>  : public   quad_function<KN_<Complex>,KN_<Complex>,fMPI_Comm,fMPI_Op,long> {
  static long  f(Stack, KN_<Complex>  const  & s, KN_<Complex>  const  &r,  fMPI_Comm const & comm,fMPI_Op const &op)  
  { 
    int chunk = s.N();
    ffassert(chunk==r.N());
#ifdef MPI_DOUBLE_COMPLEX
    return MPI_Allreduce( (void *) (Complex*)s,(void *) (Complex*)r, chunk , MPI_DOUBLE_COMPLEX,op,comm);
#else
    chunk *=2;
    return MPI_Allreduce( reinterpret_cast<void*> (s), reinterpret_cast<void*> (r), chunk , MPI_DOUBLE,op,comm);
#endif	
  }
};
// /*
// template<>
// struct Op_Reducescatter  : public   quad_function<KN_<Complex>,KN_<Complex>,fMPI_Comm,fMPI_Op,long> {
//     static long  f(Stack, KN_<Complex>  const  & s, KN_<Complex>  const  &r,  fMPI_Comm const & comm,fMPI_Op const &op)  
//     { 
// 	int chunk = s.N();
// 	ffassert(chunk==r.N());
//         //    chunk est un tableau ????
// 	MPI_Op oop = reinterpret_cast<MPI_Op> (op);
// 	return MPI_Reduce_scatter( (void *) (Complex*)s,(void *) (Complex*)r, chunk , MPI_DOUBLE_COMPLEX,op,comm);	
//     }
// };*/

template<>
struct Op_Reduce1<Complex>  : public   quad_function<R*,R*,MPIrank,fMPI_Op,long> {
  static long  f(Stack, R*  const  & s, R*  const  &r,  MPIrank const & root, fMPI_Op const &op)  
  { 
#ifdef MPI_DOUBLE_COMPLEX
    int chunk = 1;
    return MPI_Reduce( (void *) (Complex*)s, (void *) (Complex*)r, chunk , MPI_DOUBLE_COMPLEX,op,root.who,root.comm);
#else
    int chunk = 2;
    return MPI_Reduce( reinterpret_cast<void*> (s), reinterpret_cast<void*> (r), chunk , MPI_DOUBLE,op,root.who,root.comm);
#endif
  }
};

// Add J. Morice
template<>
struct Op_Gather1<Complex> : public   ternary_function<Complex* ,KN_<Complex>,MPIrank,long> {
    static long  f(Stack, Complex * const  & s, KN_<Complex>  const  &r,  MPIrank const & root)  
  { 
    
    int mpisizew;
    MPI_Comm_size(root.comm, &mpisizew); 
    int chunk = 1;
    ffassert(r.N()==mpisizew*chunk );
#ifdef MPI_DOUBLE_COMPLEX  
   
    return MPI_Gather( (void *) (Complex*) s, chunk, MPI_DOUBLE_COMPLEX,
			   (void *) (Complex*) r, chunk, MPI_DOUBLE_COMPLEX, root.who, root.comm);
#else
    chunk = 2;
    return MPI_Gather( reinterpret_cast<void*> (s), chunk, MPI_DOUBLE,
			   reinterpret_cast<void*> (r), chunk, MPI_DOUBLE, root.who, root.comm);
#endif	
  }
};

// Fin Add J. Morice


template<>
struct Op_Gather3<Complex> : public   ternary_function<KN_<Complex>,KN_<Complex>,MPIrank,long> {
    static long  f(Stack, KN_<Complex>  const  & s, KN_<Complex>  const  &r,  MPIrank const & root)  
  { 
    
    int mpisizew;
    MPI_Comm_size(root.comm, &mpisizew); 
    int chunk = r.N()/mpisizew;
    ffassert(r.N()==mpisizew*chunk && chunk==s.N());
#ifdef MPI_DOUBLE_COMPLEX     
    return MPI_Gather( (void *) (Complex*)s, chunk, MPI_DOUBLE_COMPLEX,
			   (void *) (Complex*)r, chunk, MPI_DOUBLE_COMPLEX,root.who,root.comm);	
#else
    chunk *= 2; 
    return MPI_Gather( reinterpret_cast<void*> (s), chunk, MPI_DOUBLE,
		       reinterpret_cast<void*> (r), chunk, MPI_DOUBLE,root.who,root.comm);
#endif
  }
};


template<>
//struct Op_Gatherv3 : public penta_function<KN_<Complex>,KN_<Complex>, MPIrank, KN_<long>, KN_<long>, long> {
long  Op_Gatherv3<Complex>(KN_<Complex>  const  & s, KN_<Complex>  const  &r,  MPIrank const & root, KN_<long> const & recvcount, KN_<long> const & displs)  
{ 
    
  int mpirankw;
  MPI_Comm_rank(root.comm, &mpirankw); 
  int mpisizew;
  MPI_Comm_size(root.comm, &mpisizew); 
  ffassert( recvcount.N() == displs.N() && recvcount.N() == mpisizew);
  
  if( mpirankw == root.who){
    long sum=0;
    for(int ii=0; ii< recvcount.N(); ii++)
      sum+=recvcount[ii];
    ffassert( sum == r.N() );
  }
  KN<int> INTrecvcount(recvcount.N());
  KN<int> INTdispls(displs.N());
#ifdef MPI_DOUBLE_COMPLEX   
  for(int ii=0; ii< recvcount.N(); ii++){
    INTrecvcount[ii]= recvcount[ii];
    INTdispls[ii]= displs[ii];
  }
  return MPI_Gatherv( (void *) (Complex*)s, s.N(), MPI_DOUBLE_COMPLEX,
		      (void *) (Complex*)r, INTrecvcount, INTdispls,MPI_DOUBLE_COMPLEX,root.who,root.comm);
#else
  for(int ii=0; ii< recvcount.N(); ii++){
    INTrecvcount[ii]= 2*recvcount[ii];
    INTdispls[ii]= 2*displs[ii];
  }
  return MPI_Gatherv( reinterpret_cast<void*> (s), 2*s.N(), MPI_DOUBLE,
		      reinterpret_cast<void*> (r), INTrecvcount, INTdispls,MPI_DOUBLE,root.who,root.comm);
#endif
}

// Fin Add J. Morice communication entre complexe

MPIrank mpiwho(long i) { return MPIrank(i);}
MPIrank mpiwho(long i,fMPI_Comm comm) { return MPIrank(i,comm,Syncro_block);}
MPIrank mpiwho(fMPI_Comm comm,long i) { return MPIrank(i,comm,Syncro_block);}
MPIrank mpiwhob(long i) { return MPIrank(i);}
MPIrank mpiwhob(long i,fMPI_Comm comm) { return MPIrank(i,comm,Syncro_block);}

MPIrank mpiwho_(const long &i,const fMPI_Comm &comm,fMPI_Request * const &rq) { return MPIrank(i,comm,*rq);}
MPIrank mpiwho_(const long &i,fMPI_Request * const &rq) { return MPIrank(i, MPI_COMM_WORLD ,*rq);}

long mpiWait(fMPI_Request * frq) { 
 MPI_Request * rq= *frq; 
 MPI_Status status;
 long res=MPI_SUCCESS;
  while(rq && *rq!=MPI_REQUEST_NULL)
    {
      res == MPI_Wait(rq,&status); 
      DoOnWaitMPIRequest(rq);
    }
   

}

long mpiBarrier(fMPI_Comm * comm) 
{    
  return MPI_Barrier(*comm);
}



long mpiWaitAny(KN<MPI_Request>* rq)
{ 
  MPI_Status status;
  int index;
  //cout << "mpiWaitAny " <<rq->N() << " in "  <<  endl ;
  do {
      MPI_Waitany(rq->N(),*rq,&index,&status);
      DoOnWaitMPIRequest(&(*rq)[index]);
  }  while( (MPI_UNDEFINED!= index) &&  ((*rq)[index]!=MPI_REQUEST_NULL));
    
  //cout << "mpiWaitAny " <<rq->N() << " out " << index <<  endl ;
  return index;
}

MPIrank * set_copympi( MPIrank* const & a,const MPIrank & b){ *a=b;return a;}
 

long mpiSize(fMPI_Comm  cmm) { 
    int s=0;
 //   fMPI_Comm_rank(MPI_COMM_WORLD, &s); /* local */ 
    if(cmm != MPI_COMM_NULL)
    MPI_Comm_size(cmm, &s); /* local */ 
    return s;
}
long mpiRank(fMPI_Comm  cmm) { 
    int s=MPI_UNDEFINED;
    if(cmm != MPI_COMM_NULL)
    MPI_Comm_rank(cmm, &s); /* local */ 
 //   MPI_Comm_size(MPI_COMM_WORLD, &s); /* local */ 
    return s;
}

AnyType InitializeGroup(Stack stack,const AnyType &x){
    MPI_Group *g=*PGetAny<fMPI_Group>(x);
    *g=0;
    MPI_Comm_group(MPI_COMM_WORLD, g);
    return  g;
}
AnyType DeleteGroup(Stack stack,const AnyType &x){
    MPI_Group *g=*PGetAny<fMPI_Group>(x);
    MPI_Group_free(g);
    return  Nothing;
}
AnyType InitializeComm(Stack stack,const AnyType &x){
    MPI_Comm *comm= *PGetAny<fMPI_Comm>(x);
    *comm=0;
    MPI_Comm_dup(MPI_COMM_WORLD, comm);
    return  comm;
}
AnyType DeleteComm(Stack stack,const AnyType &x){
    MPI_Comm *comm= *PGetAny<fMPI_Comm>(x);
    if(comm && (*comm != MPI_COMM_NULL))
    MPI_Comm_free(comm);
    return  Nothing;
}
AnyType InitializeRequest(Stack stack,const AnyType &x){
    MPI_Request *comm=*PGetAny<fMPI_Request>(x);
    *comm=0;
    
    return  comm;
}
AnyType DeleteRequest(Stack stack,const AnyType &x){
    MPI_Request *comm=*PGetAny<fMPI_Request>(x);
    if(comm) MPI_Request_free(comm);
    return  Nothing;
}
//  Hack to Bypass a bug in freefem FH  ... 
template<> 
class ForEachType<MPI_Group>:  public basicForEachType{public:// correction july 2009..... FH  Hoooo....  (Il y a un bug DUR DUR FH  ...) 
    ForEachType(Function1 iv=0,Function1 id=0,Function1 OOnReturn=0):basicForEachType(typeid(MPI_Group),sizeof(MPI_Group),0,0,iv,id,OOnReturn) { }
};

template<> 
class ForEachType<fMPI_Comm>:  public basicForEachType{public:// coorection july 2009..... FH  Hoooo....  (Il y a un bug DUR DUR FH  ...) 
    ForEachType(Function1 iv=0,Function1 id=0,Function1 OOnReturn=0):basicForEachType(typeid(fMPI_Comm),sizeof(fMPI_Comm),0,0,iv,id,OOnReturn) {}
};

template<> 
class ForEachType<fMPI_Request>:  public basicForEachType{public:// correction july 2009..... FH  Hoooo....  (Il y a un bug DUR DUR FH  ...) 
    ForEachType(Function1 iv=0,Function1 id=0,Function1 OOnReturn=0):basicForEachType(typeid(fMPI_Request),sizeof(fMPI_Request),0,0,iv,id,OOnReturn) {}
};
// end Hack  ... 


fMPI_Group* def_group( fMPI_Group* const & a,fMPI_Comm * const &comm, KN_<long>  const & b)
{
  MPI_Group group;
  MPI_Comm_group(*comm,& group); 
  KN<int> ranks(b);
  MPI_Group_incl(group, ranks.N(),(int *) ranks, *a);
  MPI_Group_free(&group) ;
  // ici def a .. 
  //  ffassert(0); //   A AFAIRE  //  pour arete le programm 
  return a;}

fMPI_Group* def_group( fMPI_Group* const & a, KN_<long>  const & b)
{
    MPI_Group group;
    MPI_Comm comm=MPI_COMM_WORLD;
    MPI_Comm_group(comm,& group); 
    KN<int> ranks(b);
    MPI_Group_incl(group, ranks.N(),(int *) ranks, *a);
    MPI_Group_free(&group) ;
    // ici def a .. 
    //  ffassert(0); //   A AFAIRE  //  pour arete le programm 
return a;}

fMPI_Group* def_group( fMPI_Group* const & a,fMPI_Comm * const &comm)
{
  MPI_Comm_group(*comm,*a); 
  return a;}

fMPI_Comm* def_comm( fMPI_Comm* const & a,fMPI_Group* const & g)
{
    int ok=MPI_Comm_create(MPI_COMM_WORLD,*g,*a); 
    return a;
}

fMPI_Comm* def_comm( fMPI_Comm* const & a,fMPI_Comm* const & b,fMPI_Group* const & g)
{
    MPI_Comm_create(*b,*g,*a); 
    return a;
}

fMPI_Group* def_group( fMPI_Group* const & a,fMPI_Group  * const & group,KN_<long>  const & b)
{
  
  KN<int> ranks(b);
  MPI_Group_incl(*group, ranks.N(), (int *) ranks, *a);
  // ici def a .. 
  //  ffassert(0); //   A AFAIRE  //  pour arete le programm 
  return a;
}

struct Def_def_Commsplit  : public  quad_function<fMPI_Comm*,fMPI_Comm*,long,long, fMPI_Comm*> 
{
    static fMPI_Comm * f(Stack,fMPI_Comm* const & a,fMPI_Comm* const & comm,const long &color ,const long &key  )
    {    
	MPI_Comm_split(*comm,color, key, *a);
	return a;
    }
};

fMPI_Comm * mpiCommsplit(fMPI_Comm* const & a,const MPIrank &p1,const long &key  )
{    
    MPI_Comm_split(p1.comm, p1.who, key, *a);
    return a;
}




struct Def_def_Intercommcreate  : public  quad_function<fMPI_Comm*,MPIrank,MPIrank,long, fMPI_Comm*> 
{
static fMPI_Comm * f(Stack,fMPI_Comm* const & a, MPIrank const & p1, MPIrank const & p2, long const & tag )
{
    int err;
    err=MPI_Intercomm_create(p1.comm, p1.who, p2.comm,p2.who,tag, *a);     
    return a;	
}
};

fMPI_Comm * def_intercommmerge(fMPI_Comm* const & a,fMPI_Comm* const & b, const long & high)
{
    MPI_Intercomm_merge(*b, high, *a);
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
//template<class RR,class A,class B>  fMPI_Request*,KN<MPI_Request>*,long
fMPI_Request * get_elementp_( KN<MPI_Request> * const & a,const long & b){ 
  if( b<0 || a->N() <= b) 
    { cerr << " Out of bound  0 <=" << b << " < "  << a->N() << " KN<MPI_Request> * " << endl;
      ExecError("Out of bound in operator []");}
  return  reinterpret_cast<fMPI_Request *> (&((*a)[b]));}// bofBof ... 

KN<MPI_Request> * set_init0( KN<MPI_Request> * const & a,const long & b)
   { 
     a->init(b);
     for(int i=0;i<b;++i)
       (*a)[i]=MPI_REQUEST_NULL;
     return a;
   }
bool toBool(fMPI_Comm *comm)
{
  return (comm && (*comm !=MPI_COMM_NULL)); 
}
void * topVoid(fMPI_Comm *comm) {    return comm; }

template<typename T>
class Quad_Op : public E_F0 {
  typedef typename T::result_type R;
  typedef typename T::first_argument_type A;
  typedef typename T::second_argument_type B;
  typedef typename T::third_argument_type C;
  typedef typename T::fourth_argument_type D;
  
  typedef  typename T::result_type Result;
  Expression a,b,c,d;
public:
  AnyType operator()(Stack s)  const 
  {return  SetAny<R>(static_cast<R>(T::f(s, GetAny<A>((*a)(s)) ,
					 GetAny<B>((*b)(s)) ,
					 GetAny<C>((*c)(s)),
					   GetAny<D>((*d)(s))
					 )));}
  Quad_Op(Expression aa,Expression bb,Expression cc,Expression dd) : a(aa),b(bb),c(cc),d(dd) {} 
  bool MeshIndependent() const {
    return a->MeshIndependent() && b->MeshIndependent() && c->MeshIndependent()  && d->MeshIndependent();}
};


// Fin add J. Morice
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

void init_lgparallele()
  {
    if(verbosity && mpirank == 0) cout << "parallelempi ";
    using namespace Fem2D;
    Dcl_TypeandPtr<MPIrank>(0);
    
    Dcl_TypeandPtr<fMPI_Group>(0,0,InitializeGroup,DeleteGroup); 
    Dcl_TypeandPtr<fMPI_Comm>(0,0,InitializeComm,DeleteComm);  
    Dcl_Type<fMPI_Op>(); 
    Dcl_TypeandPtr<fMPI_Request>(0,0,InitializeRequest,DeleteRequest); // bof bof ... 
    Dcl_TypeandPtr_<KN_<MPI_Request> ,KN<MPI_Request>*  > 
      (0,0,0,::Destroy<KN<MPI_Request> >,
       ::ClearReturnKK_<MPI_Request,KN<MPI_Request>,KN_<MPI_Request> >,
       ::ClearReturnpKK<MPI_Request,KN<MPI_Request> >);
      
    zzzfff->Add("mpiGroup",atype<fMPI_Group*>());
    zzzfff->Add("mpiComm",atype<fMPI_Comm*>());
    zzzfff->Add("mpiRequest",atype<fMPI_Request*>());
    
    map_type_of_map[make_pair(atype<long>(),atype<fMPI_Request>())]=atype<KN<MPI_Request>*>(); // vector
    
    
    map_type[typeid(MPIrank).name()]->AddCast(new E_F1_funcT<MPIrank,MPIrank*>(UnRef<MPIrank>));
    map_type[typeid(fMPI_Group).name()]->AddCast(new E_F1_funcT<fMPI_Group,fMPI_Group*>(UnRef<fMPI_Group>));
    map_type[typeid(fMPI_Comm).name()]->AddCast(new E_F1_funcT<fMPI_Comm,fMPI_Comm*>(UnRef<fMPI_Comm>));
    map_type[typeid(bool).name()]->AddCast(new OneOperator1<bool,fMPI_Comm*>(toBool));
    map_type[typeid(void*).name()]->AddCast(new OneOperator1<void *,fMPI_Comm*>(topVoid));
      
      
    TheOperators->Add("<-", 
		      new OneOperator2_<MPIrank*,MPIrank*,MPIrank>(&set_copympi));
    
    // constructor example ... 
    TheOperators->Add("<-", 
		      new OneOperator2_<fMPI_Group*,fMPI_Group*,KN_<long> >(&def_group),
		      new OneOperator3_<fMPI_Group*,fMPI_Group*,fMPI_Group*,KN_<long> >(&def_group),
		      new OneOperator3_<fMPI_Group*,fMPI_Group*,fMPI_Comm*,KN_<long> >(&def_group),
		      new OneOperator2_<fMPI_Group*,fMPI_Group*,fMPI_Comm*>(&def_group));
      
      
       /*
       fMPI_Comm * mpiCommsplit(fMPI_Comm* const & a,const MPIrank &p1,const long &rk  )
       fMPI_Comm * def_intercommmerge(fMPI_Comm* const & a,fMPI_Comm* const & b, const long & high)
	fMPI_Comm * def_Intercommcreate(fMPI_Comm* const & a, MPIrank const & p1, MPIrank const & p2, long const & tag )
	 quad_function<fMPI_Comm*,MPIrank,MPIrank,long, fMPI_Comm* >
       
       */
      TheOperators->Add("<-", 
			new OneOperator2_<fMPI_Comm*,fMPI_Comm*,fMPI_Group* >(&def_comm),
			new OneOperator3_<fMPI_Comm*,fMPI_Comm*,fMPI_Comm*,fMPI_Group* >(&def_comm),
			new OneOperator3_<fMPI_Comm*,fMPI_Comm*,MPIrank,long >(&mpiCommsplit),
			new OneOperator3_<fMPI_Comm*,fMPI_Comm*,fMPI_Comm*,long >(&def_intercommmerge),
			new OneQuadOperator< Def_def_Intercommcreate, Quad_Op<Def_def_Intercommcreate>  >,
			new OneQuadOperator< Def_def_Commsplit, Quad_Op<Def_def_Commsplit>  >
			
			); 
      
    /*  code edp
	int[int] procs=[1,2,3];
	mpiGroup toto(procs);
	mpiComm comm(toto);
	
    */
      
    
    Global.Add("processor","(",new OneOperator1<MPIrank,long>(mpiwho));
    Global.Add("processor","(",new OneOperator2<MPIrank,long,fMPI_Comm>(mpiwho));
      Global.Add("processor","(",new OneOperator2<MPIrank,fMPI_Comm,long>(mpiwho));
   
    Global.Add("processorblock","(",new OneOperator1<MPIrank,long>(mpiwhob));
    Global.Add("processorblock","(",new OneOperator2<MPIrank,long,fMPI_Comm>(mpiwhob));
    
    TheOperators->Add(">>",
		      new OneBinaryOperator<Op_Readmpi<double> >,
		      new OneBinaryOperator<Op_Readmpi<Complex> >,
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
		       new OneBinaryOperator<Op_Writempi<Complex> >,
		       new OneBinaryOperator<Op_Writempi<long> > ,
		       new OneBinaryOperator<Op_Writempi<KN<double> * > > ,
		       new OneBinaryOperator<Op_Writempi<KN<long> * > > ,
		       new OneBinaryOperator<Op_Writempi<KN<Complex> * > > ,
		       new OneBinaryOperator<Op_Writempi<Mesh *> > ,
		       new OneBinaryOperator<Op_Writempi<Mesh3 *> > ,
		       new OneBinaryOperator<Op_Writempi<Matrice_Creuse<R> * > > ,
		       new OneBinaryOperator<Op_Writempi<Matrice_Creuse<Complex>* > > 
		       
		       );
     
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<double> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<long> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<Complex> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<KN<double> *> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<KN<long> *> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<KN<Complex> *> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<Mesh *> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<Mesh3 *> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<Matrice_Creuse<R> *> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<Matrice_Creuse<Complex> *> >);
     
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<double> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<long> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<Complex> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<KN<double> *> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<KN<long> *> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<KN<Complex> *> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<Mesh *> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<Mesh3 *> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<Matrice_Creuse<R> *> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<Matrice_Creuse<Complex> *> >);
     
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<double> >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<long> >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<Complex> >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<KN<double> > >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<KN<long> > >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<KN<Complex> > >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<Mesh *> >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<Mesh3 *> >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<Matrice_Creuse<R> > >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<Matrice_Creuse<Complex> > >);
      
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<double> >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<long> >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<Complex> >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<KN<double> > >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<KN<long> > >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<KN<Complex> > >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<Mesh *> >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<Mesh3 *> >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<Matrice_Creuse<R> > >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<Matrice_Creuse<Complex> > >);
      
    
      
      
      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<double> >);
      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<Complex> >);
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
      
      Global.Add("mpiAllgather","(",new OneBinaryOperator<Op_Allgather1< long > >);   // Add J. Morice
      Global.Add("mpiAllgather","(",new OneBinaryOperator<Op_Allgather1< double > >); // Add J. Morice

      Global.Add("mpiAllgather","(",new OneTernaryOperator3<Op_Allgather13< long > >);  // Add J. Morice
      Global.Add("mpiAllgather","(",new OneTernaryOperator3<Op_Allgather13< double > >);// Add J. Morice
      
      Global.Add("mpiScatter","(",new OneTernaryOperator3<Op_Scatter1< long > >);    // Add J. Morice
      Global.Add("mpiScatter","(",new OneTernaryOperator3<Op_Scatter1< double > >);  // Add J. Morice
      Global.Add("mpiScatter","(",new OneTernaryOperator3<Op_Scatter3< long > >);
      Global.Add("mpiScatter","(",new OneTernaryOperator3<Op_Scatter3< double > >);

      Global.Add("mpiGather","(",new OneTernaryOperator3<Op_Gather1< long > >);   // Add J. Morice
      Global.Add("mpiGather","(",new OneTernaryOperator3<Op_Gather1< double > >); // Add J. Morice
      Global.Add("mpiGather","(",new OneTernaryOperator3<Op_Gather3< long > >); // correction J. Morice Scatter --> Gather
      Global.Add("mpiGather","(",new OneTernaryOperator3<Op_Gather3< double > >);
      
      // Add J. Morice communication with vector of different size    
      Global.Add("mpiAlltoallv","(",new OneOperator6_<long, KN_<double>, KN_<double>, KN_<long>, KN_<long>, KN_<long>, KN_<long> >( Op_All2Allv<double> ) );
      Global.Add("mpiAlltoallv","(",new OneOperator6_<long, KN_<long>, KN_<long>, KN_<long>, KN_<long>, KN_<long>, KN_<long> >( Op_All2Allv<long> ) );
      
      Global.Add("mpiAlltoallv","(",new OneOperator7_<long, KN_<long>, KN_<long>, fMPI_Comm, KN_<long>, KN_<long>, KN_<long>, KN_<long> >( Op_All2All3v<long> ) ); 
      Global.Add("mpiAlltoallv","(",new OneOperator7_<long, KN_<double>, KN_<double>, fMPI_Comm, KN_<long>, KN_<long>, KN_<long>, KN_<long> >( Op_All2All3v<double> ) ); 

      Global.Add("mpiAllgatherv","(",new OneQuadOperator<Op_Allgatherv< long >, Quad_Op<Op_Allgatherv< long > > > );
      Global.Add("mpiAllgatherv","(",new OneQuadOperator<Op_Allgatherv< double >, Quad_Op<Op_Allgatherv< double > > >);
      Global.Add("mpiAllgatherv","(",new OneOperator5_<long, KN_<long>, KN_<long>, fMPI_Comm, KN_<long>, KN_<long> >(Op_Allgatherv3< long >) );
      Global.Add("mpiAllgatherv","(",new OneOperator5_<long, KN_<double>, KN_<double>, fMPI_Comm, KN_<long>, KN_<long> >(Op_Allgatherv3< double >) );

      Global.Add("mpiScatterv","(",new OneOperator5_<long, KN_<long>, KN_<long>, MPIrank, KN_<long>, KN_<long> >(Op_Scatterv3< long >) );
      Global.Add("mpiScatterv","(",new OneOperator5_<long, KN_<double>, KN_<double>, MPIrank, KN_<long>, KN_<long> >(Op_Scatterv3< double >) );

      Global.Add("mpiGatherv","(",new OneOperator5_<long, KN_<long>, KN_<long>, MPIrank, KN_<long>, KN_<long> >( Op_Gatherv3< long > ) );
      Global.Add("mpiGatherv","(",new OneOperator5_<long, KN_<double>, KN_<double>, MPIrank, KN_<long>, KN_<long> >( Op_Gatherv3< double > ) );
      // Fin Add J. Morice

      Global.Add("mpiReduce","(",new OneQuadOperator<Op_Reduce< double >, Quad_Op<Op_Reduce< double > > >);
      Global.Add("mpiReduce","(",new OneQuadOperator<Op_Reduce1< double >, Quad_Op<Op_Reduce1< double > > >);
      Global.Add("mpiAllReduce","(",new OneQuadOperator<Op_AllReduce< double >, Quad_Op<Op_AllReduce< double > > >);
      //    Global.Add("mpiReduceScatter","(",new OneQuadOperator<Op_Reducescatter< double >, Quad_Op<Op_Reducescatter< double > > >);

      // Add J. Morice
      Global.Add("mpiReduce","(",new OneQuadOperator<Op_Reduce< long >, Quad_Op<Op_Reduce< long > > >);
      Global.Add("mpiReduce","(",new OneQuadOperator<Op_Reduce1< long >, Quad_Op<Op_Reduce1< long > > >);
      Global.Add("mpiAllReduce","(",new OneQuadOperator<Op_AllReduce< long >, Quad_Op<Op_AllReduce< long > > >);
      //    Global.Add("mpiReduceScatter","(",new OneQuadOperator<Op_Reducescatter< long >, Quad_Op<Op_Reducescatter< long > > >);
      // fin Add J. Morice

      // Add J. Morice :: complex communication between processor 
      Global.Add("mpiAlltoall","(",new OneBinaryOperator<Op_All2All< Complex > >);
      Global.Add("mpiAllgather","(",new OneBinaryOperator<Op_Allgather< Complex > >);    
      Global.Add("mpiAlltoall","(",new OneTernaryOperator3<Op_All2All3< Complex > >);     
      Global.Add("mpiAllgather","(",new OneTernaryOperator3<Op_Allgather3< Complex > >);
        
      Global.Add("mpiAllgather","(",new OneBinaryOperator<Op_Allgather1< Complex > >);
      Global.Add("mpiAllgather","(",new OneTernaryOperator3<Op_Allgather13< Complex > >);

      Global.Add("mpiScatter","(",new OneTernaryOperator3<Op_Scatter3< Complex > >);
      Global.Add("mpiGather","(",new OneTernaryOperator3<Op_Gather3< Complex > >);

      Global.Add("mpiScatter","(",new OneTernaryOperator3<Op_Scatter1< Complex > >);
      Global.Add("mpiGather","(",new OneTernaryOperator3<Op_Gather1< Complex > >);
      
      // Add J. Morice communication with vector of different size    
      Global.Add("mpiAlltoallv","(",new OneOperator6_<long, KN_<Complex>, KN_<Complex>, KN_<long>, KN_<long>, KN_<long>, KN_<long> >( Op_All2Allv<Complex> ) );   
      Global.Add("mpiAlltoallv","(",new OneOperator7_<long, KN_<Complex>, KN_<Complex>, fMPI_Comm, KN_<long>, KN_<long>, KN_<long>, KN_<long> >( Op_All2All3v<Complex> ) ); 
      Global.Add("mpiAllgatherv","(",new OneQuadOperator<Op_Allgatherv< Complex>, Quad_Op<Op_Allgatherv< Complex > > >);
   
      Global.Add("mpiAllgatherv","(",new OneOperator5_<long, KN_<Complex>, KN_<Complex>, fMPI_Comm, KN_<long>, KN_<long> >(Op_Allgatherv3< Complex >) );
      Global.Add("mpiScatterv","(",new OneOperator5_<long, KN_<Complex>, KN_<Complex>, MPIrank, KN_<long>, KN_<long> >(Op_Scatterv3< Complex >) );
      Global.Add("mpiGatherv","(",new OneOperator5_<long, KN_<Complex>, KN_<Complex>, MPIrank, KN_<long>, KN_<long> >( Op_Gatherv3< Complex > ) );
     
      Global.Add("mpiReduce","(",new OneQuadOperator<Op_Reduce< Complex >, Quad_Op<Op_Reduce< Complex > > >);
      Global.Add("mpiReduce","(",new OneQuadOperator<Op_Reduce1< Complex >, Quad_Op<Op_Reduce1< Complex > > >);
      Global.Add("mpiAllReduce","(",new OneQuadOperator<Op_AllReduce< Complex >, Quad_Op<Op_AllReduce< Complex > > >);
      // Fin Add J. Morice :: complex communication between processor 
      
      Global.New("mpirank",CConstant<long>(mpirank));
      Global.New("mpisize",CConstant<long>(mpisize));
     static long mpiUndefined=MPI_UNDEFINED, mpiAnySource =  MPI_ANY_SOURCE,mpiAnyTag=MPI_ANY_TAG ;
     static fMPI_Comm fmpiWorld=MPI_COMM_WORLD;
     
     Global.New("mpiUndefined",CConstant<long>(mpiUndefined));
     Global.New("mpiAnySource",CConstant<long>(mpiAnySource));
      Global.New("mpiAnyTag",CConstant<long>(mpiAnyTag));
      
      
      
      Global.New("mpiCommWorld",CConstant<fMPI_Comm*>(&fmpiWorld));   
      // add FH 
      
      
     Global.Add("mpiWtime","(",new OneOperator0<double>(MPI_Wtime));    
     Global.Add("mpiWtick","(",new OneOperator0<double>(MPI_Wtick));    
     Global.Add("processor","(",new OneOperator3_<MPIrank,long,fMPI_Comm,fMPI_Request*>(mpiwho_));
     Global.Add("processor","(",new OneOperator2_<MPIrank,long,fMPI_Request*>(mpiwho_));
     Global.Add("mpiWait","(",new OneOperator1<long,fMPI_Request*>(mpiWait));
     Global.Add("mpiWaitAny","(",new OneOperator1<long,KN<MPI_Request>*>(mpiWaitAny));
     Global.Add("mpiSize","(",new OneOperator1<long,fMPI_Comm>(mpiSize)); 
     Global.Add("mpiRank","(",new OneOperator1<long,fMPI_Comm>(mpiRank)); 
      Global.Add("mpiBarrier","(",new OneOperator1<long,fMPI_Comm*>(mpiBarrier)); 
      
      static   fMPI_Op op_max(MPI_MAX);
      static   fMPI_Op op_min(MPI_MIN);
      static   fMPI_Op op_sum(MPI_SUM);
      static   fMPI_Op op_prod(MPI_PROD);
      static   fMPI_Op op_land(MPI_LAND);
      
      static   fMPI_Op op_lor(MPI_LOR);
      static   fMPI_Op op_lxor(MPI_LXOR);
      static   fMPI_Op op_band(MPI_BAND);
      static   fMPI_Op op_bor(MPI_BOR);
      static   fMPI_Op op_bxor(MPI_BXOR);
      static   fMPI_Op op_maxloc(MPI_MAXLOC);
      static   fMPI_Op op_minloc(MPI_MINLOC);
      
      Global.New("mpiMAX",CConstant<fMPI_Op>(op_max));
      Global.New("mpiMIN",CConstant<fMPI_Op>(op_min));
      Global.New("mpiSUM",CConstant<fMPI_Op>(op_sum));
      Global.New("mpiPROD",CConstant<fMPI_Op>(op_prod));
      Global.New("mpiLAND",CConstant<fMPI_Op>(op_land));
      
      Global.New("mpiLOR",CConstant<fMPI_Op>(op_lor));
      Global.New("mpiLXOR",CConstant<fMPI_Op>(op_lxor));
      Global.New("mpiBAND",CConstant<fMPI_Op>(op_band));
      Global.New("mpiBXOR",CConstant<fMPI_Op>(op_bxor));
      //  sur des pair bof bof ...
      Global.New("mpiMAXLOC",CConstant<fMPI_Op>(op_maxloc));
      Global.New("mpiMINLOC",CConstant<fMPI_Op>(op_minloc));
      
      
      TheOperators->Add("<-", 
			new OneOperator2_<KN<MPI_Request> *,KN<MPI_Request> *,long>(&set_init0)
			);
      atype<KN<MPI_Request>* >()->Add("[","",new OneOperator2_<fMPI_Request*,KN<MPI_Request>*,long >(get_elementp_));    
      
  }
void end_parallele()
{
    MPI_Finalize();
    if(mpirank==0) cout << "FreeFem++-mpi finalize correctly .\n" << flush ; 
    else if(verbosity>5)  cout << '.' << endl ;
}
//   MPI::COMM_WORLD.Recv(&msg, 1, MPI::INT, from, MPI::ANY_TAG);
//    MPI::COMM_WORLD.Isend(&msg, 1, MPI::INT, to, 4);
