#include <complex>
#include "rgraph.hpp"
#include "AFunction.hpp"
#include "RNM.hpp"
#include "Operator.hpp"
#include "fem.hpp"
//  pas terrible .... 

#undef MPICH_SKIP_MPICXX
#define  MPICH_SKIP_MPICXX
#undef MPICH_IGNORE_CXX_SEEK
#define MPICH_IGNORE_CXX_SEEK
#include "mpi.h"

#ifdef HAVE_MPICXX_H
#include <mpicxx.h>
#elif HAVE_MPI_CXX_MPICXX_H
#include <mpi/cxx/mpicxx.h>
#elif HAVE_OMPI_MPI_CXX_MPICXX_H
#include <ompi/mpi/cxx/mpicxx.h>
#elif HAVE_OPENMPI_OMPI_MPI_CXX_MPICXX_H
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#elif  HAVE_MPI___H
#include <mpi++.h>
#else
#error "no mpixx.h or mpi++.h file" 
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
template<> struct MPI_TYPE<long>      {static const MPI::Datatype TYPE(){return MPI::LONG;}};
template<> struct MPI_TYPE<double>    {static const MPI::Datatype TYPE(){return MPI::DOUBLE;}};
#ifdef MPI_DOUBLE_COMPLEX_
template<> struct MPI_TYPE<Complex>   {static const MPI::Datatype TYPE(){return MPI::DOUBLE_COMPLEX;}};
#endif
template<class T> struct MPI_WHAT {};
template<> struct MPI_WHAT<long>      {static const int WHAT=101;};
template<> struct MPI_WHAT<double>    {static const int WHAT=102;};
template<> struct MPI_WHAT<Complex>   {static const int WHAT=103;};
template<> struct MPI_WHAT<KN<long> >   {static const int WHAT=104;};
template<> struct MPI_WHAT<KN<double> >   {static const int WHAT=105;};
template<> struct MPI_WHAT<KN<Complex> >   {static const int WHAT=106;};

  void initparallele(int &, char **&);
  void init_lgparallele();  
  
 extern long mpirank ;
 extern long mpisize ;


const size_t sizempibuf = 1024*32;
const int Meshtag =1000;
struct MPIrank {
    
   int who; 
   MPIrank(int i=0) : who(i) {ffassert(i>=0 && i < mpisize);} 
   
   const MPIrank & operator<<(double a)const  {
     MPI::COMM_WORLD.Isend(&a, 1, MPI::DOUBLE, who, 4);
      return *this;
   }
   const MPIrank & Bcast(double & a) const {
     (void)  MPI::COMM_WORLD.Bcast(&a, 1, MPI::DOUBLE, who);
      return *this;
   }
   const MPIrank & operator>>(double & a) const {
      MPI::COMM_WORLD.Recv(&a, 1, MPI::DOUBLE, who, 4);
      return *this;
   }
   const MPIrank & operator<<(long a) const {
     MPI::COMM_WORLD.Isend(&a, 1, MPI::LONG, who, 5);
      return *this;
   }
   const MPIrank & Bcast(long & a) const {
     (void)  MPI::COMM_WORLD.Bcast(&a, 1, MPI::LONG, who);
      return *this;
   }
   const MPIrank & operator>>(long & a) const {
      MPI::COMM_WORLD.Recv(&a, 1, MPI::LONG, who, 5);
      return *this;
   }
   
   const MPIrank & operator>>(KN<double> & a) const {
       assert(&a);
      int n= a.N();
      MPI::COMM_WORLD.Recv((double *) a, n, MPI::DOUBLE, who, 10);
      ffassert(a.N()==n);
      return *this;
   }
    const MPIrank & operator>>(KN<long> & a) const {
	assert(&a);
	int n= a.N();
	MPI::COMM_WORLD.Recv((long *) a, n, MPI::LONG, who, 11);
	ffassert(a.N()==n);
	return *this;
    }

    const MPIrank & operator>>(KN<Complex> & a) const {
	assert(&a);
	int n= a.N();
#ifdef MPI_DOUBLE_COMPLEX_
	MPI::COMM_WORLD.Recv((Complex *) a, n, MPI::DOUBLE_COMPLEX, who, 12);
#else
	n*=2;
	MPI::COMM_WORLD.Recv((Complex *) a, n, MPI::DOUBLE , who, 12);
	n /= 2;
#endif
	ffassert(a.N()==n);
	return *this;
    }
    
   const MPIrank & Bcast(KN<double> & a) const {
       assert(&a);
      int n= a.N();
      (void)  MPI::COMM_WORLD.Bcast((double *) a, n, MPI::DOUBLE, who);
     ffassert(a.N()==n);
      return *this;
   }

    const MPIrank & Bcast(KN<Complex> & a) const {
	assert(&a);
	int n= a.N();
#ifdef MPI_DOUBLE_COMPLEX_
	(void)  MPI::COMM_WORLD.Bcast((Complex *) a, n, MPI::DOUBLE_COMPLEX, who);
#else
	n*=2;
	(void)  MPI::COMM_WORLD.Bcast((Complex *) a, n, MPI::DOUBLE, who);
	n /= 2;
#endif
	ffassert(a.N()==n);
	return *this;
    }

    const MPIrank & Bcast(KN<long> & a) const {
	assert(&a);
	int n= a.N();
	(void)  MPI::COMM_WORLD.Bcast((long *) a, n, MPI::LONG, who);
	ffassert(a.N()==n);
	return *this;
    }
    
   const MPIrank & operator<<(const KN<double> *aa)const  {
     const KN<double> & a=*aa;
      ffassert(a); 
      int n= a.N();
      MPI::COMM_WORLD.Isend((double *) a, n, MPI::DOUBLE, who, 10);
      return *this;
   }
 
    const MPIrank & operator<<(const KN<long> *aa)const  {
	const KN<long> & a=*aa;
	ffassert(a); 
	int n= a.N();
	MPI::COMM_WORLD.Isend((long *) a, n, MPI::LONG, who, 11);
	return *this;
    }
 
    const MPIrank & operator<<(const KN<Complex> *aa)const  {
	const KN<Complex> & a=*aa;
	ffassert(a); 
	int n= a.N();
#ifdef MPI_DOUBLE_COMPLEX_
	MPI::COMM_WORLD.Isend((Complex *) a, n, MPI::DOUBLE_COMPLEX, who, 12);
#else
	n *=2;
	MPI::COMM_WORLD.Isend((Complex *) a, n, MPI::DOUBLE, who, 12);
	n /= 2;
#endif
	return *this;
    }
    
    
   const MPIrank & Bcast(const KN<double> *aa)const  {
      const KN<double> & a=*aa;
      assert(a); 
      int n= a.N();
      (void) MPI::COMM_WORLD.Bcast((double *) a, n, MPI::DOUBLE, who);
      ffassert(a.N()==n);
      return *this;
   }

    const MPIrank & Bcast(const KN<long> *aa)const  {
	const KN<long> & a=*aa;
	assert(a); 
	int n= a.N();
	(void) MPI::COMM_WORLD.Bcast((long *) a, n, MPI::LONG, who);
	ffassert(a.N()==n);
	return *this;
    }

    const MPIrank & Bcast(const KN<Complex> *aa)const  {
	const KN<Complex> & a=*aa;
	assert(a); 
	int n= a.N();
#ifdef MPI_DOUBLE_COMPLEX_
	(void) MPI::COMM_WORLD.Bcast((Complex *) a, n, MPI::DOUBLE_COMPLEX, who);
#else
	n *=2;
	(void) MPI::COMM_WORLD.Bcast((Complex *) a, n, MPI::DOUBLE, who);
	n /= 2;
#endif
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
      (void) MPI::COMM_WORLD.Bcast( &nbsize, 1, MPI::LONG, who);
       if (who != mpirank)
          buf= new Serialize(nbsize,Fem2D::Mesh::magicmesh);
       assert(nbsize);
       if(verbosity>2) 
         cout << " size to bcast : " << nbsize << " mpirank : " << mpirank << endl;
       
       MPI::COMM_WORLD.Bcast( (char *)(*buf),nbsize, MPI::BYTE, who);     
        
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
   
   const MPIrank & operator<<(Fem2D::Mesh *  a) const {
     if(verbosity>1) 
     cout << " MPI << (mesh *) " << a << endl;
      ffassert(a);
      Serialize  buf=(*a).serialize();       
      buf.mpisend(*this,Meshtag);
      return *this;
   }
   
   const MPIrank & operator>>(Fem2D::Mesh *& a) const {
     if(verbosity>1) 
     cout << " MPI >> (mesh *) &" << a << endl;
      Serialize buf(*this,Fem2D::Mesh::magicmesh,Meshtag);
      if (a) (*a).decrement();
      a= new Fem2D::Mesh(buf);
//  add 3 line FH 08/12/2003  forget build quadtree sorry      
      Fem2D::R2 Pn,Px;
      a->BoundingBox(Pn,Px);
      a->quadtree=new Fem2D::FQuadTree(a,Pn,Px,a->nv);
      return *this;
   }
  operator int () const { return who;}     
};

void Serialize::mpisend(const MPIrank & rank,long tag)
{
      char * pp = p-sizeof(long);
      long countsave=count(); // save count 
      count()=lg; // store length in count 
      int l=lg+sizeof(long);
      if(verbosity>1) 
         cout << " -- send from  " << mpirank << " to " << rank << " serialized " << what 
              <<   ", l=" << l << ", tag=" << tag << endl;
      if (l <=sizempibuf)
        MPI::COMM_WORLD.Isend(pp,l, MPI::BYTE, rank, tag);
      else {
         MPI::COMM_WORLD.Isend(pp,sizempibuf, MPI::BYTE, rank, tag);
         MPI::COMM_WORLD.Isend(pp+sizempibuf,l-sizempibuf, MPI::BYTE, rank, tag+1);
      }
      if(verbosity>1) 
         cout << "    ok send is arrived " << endl;      
      count()=countsave; // restore count 
}

Serialize::Serialize(const MPIrank & rank,const char * wht,long tag)
 :what(wht) 
{
      if(verbosity>1) 
         cout << " -- waiting " << mpirank << " from  " << rank << " serialized " << what 
              << " tag = " << tag <<  endl;
   char * buf= new char [sizempibuf];
   MPI::COMM_WORLD.Recv(buf, sizempibuf, MPI::BYTE, rank, tag);
   lg = * (long *) (void *) buf;
   int l=lg+sizeof(long);
   char * pp= new char[l]  ;
   if ( l <= sizempibuf) 
      memcpy(pp,buf,l);
   else 
      {
        memcpy(pp,buf,sizempibuf);
        MPI::COMM_WORLD.Recv(pp+sizempibuf,l-sizempibuf, MPI::BYTE, rank, tag+1)  ;       
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


 MPIrank * set_copympi( MPIrank* const & a,const MPIrank & b)
 { *a=b;return a;}
 
  void initparallele(int &argc, char **& argv)
  {
     MPI::Init(argc, argv);
     mpirank = MPI::COMM_WORLD.Get_rank();
     mpisize = MPI::COMM_WORLD.Get_size();
     cout << "initparallele rank " <<  mpirank << " on " << mpisize << endl;
  }
  
     
void init_lgparallele()
  {
    if(verbosity) cout << "parallelempi ";
     using namespace Fem2D;
     Dcl_TypeandPtr<MPIrank>(0);
     map_type[typeid(MPIrank).name()]->AddCast(
       new E_F1_funcT<MPIrank,MPIrank*>(UnRef<MPIrank>));
       
     TheOperators->Add("<-", 
       new OneOperator2_<MPIrank*,MPIrank*,MPIrank>(&set_copympi));
     Global.Add("processor","(",new OneOperator1<MPIrank,long>(mpiwho));

     TheOperators->Add(">>",
		       new OneBinaryOperator<Op_Readmpi<double> >,
		       new OneBinaryOperator<Op_Readmpi<long> > ,
		       new OneBinaryOperator<Op_Readmpi<KN<double> > > ,
		       new OneBinaryOperator<Op_Readmpi<KN<long> > > ,
		       new OneBinaryOperator<Op_Readmpi<KN<Complex> > > ,
		       new OneBinaryOperator<Op_Readmpi<Mesh *> > 
       );
     TheOperators->Add("<<",
       new OneBinaryOperator<Op_Writempi<double> >,
       new OneBinaryOperator<Op_Writempi<long> > ,
       new OneBinaryOperator<Op_Writempi<KN<double> * > > ,
       new OneBinaryOperator<Op_Writempi<KN<long> * > > ,
       new OneBinaryOperator<Op_Writempi<KN<Complex> * > > ,
       new OneBinaryOperator<Op_Writempi<Mesh *> > 
       );
       
    Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<double> >);
    Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<long> >);
    Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<KN<double> > >);
    Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<KN<long> > >);
    Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<KN<Complex> > >);
    Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<Mesh *> >);
    
 
     Global.New("mpirank",CPValue<long>(mpirank));
     Global.New("mpisize",CPValue<long>(mpisize));
    
  }
  void end_parallele()
   {
    MPI::Finalize();
    if(mpirank==0) cout << "FreeFem++-mpi finalize correctly .\n" << flush ; 
    else   cout << '.' << flush ;
   }
//   MPI::COMM_WORLD.Recv(&msg, 1, MPI::INT, from, MPI::ANY_TAG);
//    MPI::COMM_WORLD.Isend(&msg, 1, MPI::INT, to, 4);
