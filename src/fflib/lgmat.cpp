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
#ifdef __MWERKS__
#pragma optimization_level 0
//#pragma inline_depth(0) 
#endif
/*
#include "error.hpp"
#include "RNM.hpp"
#include <set>
#include <vector>
#include <cstdio>
#include <cstring>
#include <complex>

#include  <cmath>
#include  <iostream>
using namespace std;
#include "FESpacen.hpp"

#include "AFunction.hpp"
#include "rgraph.hpp"
#include "MatriceCreuse_tpl.hpp"

//#include "fem3.hpp"
#include "MeshPoint.hpp"
#include "Operator.hpp" 

#include "lex.hpp"
#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"
*/
#include "ff++.hpp" 
#include "array_resize.hpp" 
#include "AFunction_ext.hpp"
#include "CGNL.hpp"

namespace bamg { class Triangles; }
namespace Fem2D { void DrawIsoT(const R2 Pt[3],const R ff[3],const RN_ & Viso); }

extern const basicForEachType *aatypeknlongp; //// for  compilation error with g++ 3.2.2
#include "BamgFreeFem.hpp"



// Debut FH Houston -------- avril 2004 ---------
//  class for operator on sparse matrix 
//   ------------------------------------
//  pour le addition de matrice ----matrice creuse in french)
//  MatriceCreuse    real class for matrix sparse
//  Matrice_Creuse   class for matrix sparse +  poiteur on FE space def 
//         to recompute matrice in case of mesh change
//  list<triplet<R,MatriceCreuse<R> *,bool> > * is liste of 
//  \sum_i a_i*A_i  where a_i is a scalare and A_i is a sparse matrix
//
 


template<class R> 
list<triplet<R,MatriceCreuse<R> *, bool > > * to(Matrice_Creuse<R> * M)
{
  list<triplet<R,MatriceCreuse<R> *,bool> >  * l=new list<triplet<R,MatriceCreuse<R> *,bool> >;
   l ->push_back(make_triplet<R,MatriceCreuse<R> *,bool>(1,M->A,false));
  return  l;
}
template<class R> 
list<triplet<R,MatriceCreuse<R> *, bool > > * to(Matrice_Creuse_Transpose<R>  M)
{
  list<triplet<R,MatriceCreuse<R> *,bool> >  * l=new list<triplet<R,MatriceCreuse<R> *,bool> >;
   l ->push_back(make_triplet<R,MatriceCreuse<R> *,bool>(1,M.A->A,true));
  return  l;
}

template<class R>
AnyType M2L3 (Stack , const AnyType & pp)
{
    return to(GetAny<Matrice_Creuse<R> *>(pp));
}


template<class R>
AnyType tM2L3 (Stack , const AnyType & pp)
{
    return to(GetAny<Matrice_Creuse_Transpose<R> >(pp));
}


template<class R> 
struct Op2_ListCM: public binary_function<R,Matrice_Creuse<R> *,list<triplet<R,MatriceCreuse<R> *,bool> > *> 
 { 
   typedef triplet<R,MatriceCreuse<R> *,bool>  P;
   typedef list<P> L;
   typedef L * RR;
   typedef R  AA;
   typedef Matrice_Creuse<R> * BB;
   
 static  RR f(const   AA & a,const   BB & b)  
  {
    RR r=  new list<P> ;
    r ->push_back(make_triplet<R,MatriceCreuse<R> *>(a,b->A,false));
    return r;}
};

template<class R> 
struct Op2_ListMC: public binary_function<Matrice_Creuse<R> *,R,list<triplet<R,MatriceCreuse<R> *,bool> > *> 
 { 
   typedef triplet<R,MatriceCreuse<R> *,bool>  P;
   typedef list<P> L;
   typedef L * RR;
   typedef R  AA;
   typedef Matrice_Creuse<R> * BB;
   
 static  RR f(const   BB & b,const   AA & a)  
  {
    RR r=  new list<P> ;
    r ->push_back(make_triplet<R,MatriceCreuse<R> *>(a,b->A,false));
    return r;}
};
//  ADD FH 16/02/2007

template<class R> 
struct Op2_ListCMt: public binary_function<R,Matrice_Creuse_Transpose<R> ,list<triplet<R,MatriceCreuse<R> *,bool> > *> 
{ 
    typedef triplet<R,MatriceCreuse<R> *,bool>  P;
    typedef list<P> L;
    typedef L * RR;
    typedef R  AA;
    typedef Matrice_Creuse_Transpose<R>  BB;
    
    static  RR f(const   AA & a,const   BB & b)  
    {
	RR r=  new list<P> ;
	r ->push_back(make_triplet<R,MatriceCreuse<R> *>(a,b.A->A,true));
	return r;}
};

template<class R> 
struct Op2_ListMtC: public binary_function<Matrice_Creuse_Transpose<R> ,R,list<triplet<R,MatriceCreuse<R> *,bool> > *> 
{ 
    typedef triplet<R,MatriceCreuse<R> *,bool>  P;
    typedef list<P> L;
    typedef L * RR;
    typedef R  AA;
    typedef Matrice_Creuse_Transpose<R> BB;
    
    static  RR f(const   BB & b,const   AA & a)  
    {
	RR r=  new list<P> ;
	r ->push_back(make_triplet<R,MatriceCreuse<R> *>(a,b.A->A,true));
	return r;}
};
// FIN ADD 16/02/2007



template<class R> 
struct Op1_LCMd: public unary_function<list<triplet<R,MatriceCreuse<R> *,bool> > *,
list<triplet<R,MatriceCreuse<R> *,bool> > *  >
{  //  - ...
    typedef triplet<R,MatriceCreuse<R> *,bool>  P;
    typedef list<P> L;
    typedef L * RR;
    
    static   RR f(const RR & l)  
    { 
	typedef typename list<triplet<R,MatriceCreuse<R> *,bool> >::iterator lci;
	for (lci i= l->begin();i !=l->end();++i)
	    i->first *= R(-1);
	return l;
    }    
    
};

template<class R> 
struct Op2_ListCMCMadd: public binary_function<list<triplet<R,MatriceCreuse<R> *,bool> > *,
                                               list<triplet<R,MatriceCreuse<R> *,bool> > *,
                                               list<triplet<R,MatriceCreuse<R> *,bool> > *  >
{  //  ... + ...
   typedef triplet<R,MatriceCreuse<R> *,bool>  P;
   typedef list<P> L;
   typedef L * RR;

  static   RR f(const RR & a,const RR & b)  
  { 
    a->insert(a->end(),b->begin(),b->end());
      
    delete b;
    return a;
  }    
    
};

template<class R> 
struct Op2_ListMCMadd: public binary_function<Matrice_Creuse<R> *,
                                              list<triplet<R,MatriceCreuse<R> *,bool> > *,                                               
                                               list<triplet<R,MatriceCreuse<R> *,bool> > *  >
{  //  M + ....
   typedef triplet<R,MatriceCreuse<R> *,bool> P;
   typedef list<P> L;
   typedef L * RR;
   typedef Matrice_Creuse<R> * MM;

  static   RR f(const MM & a,const RR & b)  
  { 
    
    b->push_front(make_triplet<R,MatriceCreuse<R> *>(R(1.),a->A,false));
    return b;
  }
    
    
};

template<class R> 
struct Op2_ListCMMadd: public binary_function< list<triplet<R,MatriceCreuse<R> *,bool> > *,                                                                                              
                                               Matrice_Creuse<R> * ,
                                               list<triplet<R,MatriceCreuse<R> *,bool> > *>
{  //   .... + M
   typedef triplet<R,MatriceCreuse<R> *,bool> P;
   typedef list<P> L;
   typedef L * RR;
   typedef Matrice_Creuse<R> * MM;

  static   RR f(const RR & a,const MM & b)  
  { 
    
    a->push_back(make_triplet<R,MatriceCreuse<R> *,bool>(R(1.),b->A,false));
    return a;
  }
    
    
};

template<class R> 
struct Op2_ListMMadd: public binary_function< Matrice_Creuse<R> *,
                                              Matrice_Creuse<R> * ,
                                              list<triplet<R,MatriceCreuse<R> *,bool> > *>
{  //  M + M
   typedef triplet<R,MatriceCreuse<R> *,bool> P;
   typedef list<P> L;
   typedef L * RR;
   typedef Matrice_Creuse<R> * MM;

  static   RR f(const MM & a,const MM & b)  
  { 
    L * l=to(a);
    l->push_back(make_triplet<R,MatriceCreuse<R> *>(R(1.),b->A,false));
    return l;
  }
    
    
};
// Fin Add FH Houston --------

// for Jolivet  to build restriction  jan 2014
// t[int] I= restrict(VCh,VGh,IPG); //  ou
template<class pfes>
class RestrictArray : public OneOperator { public:
    template< typename T > struct Base { typedef  T B; };
    template< typename T > struct Base< T* >{ typedef T B;};

    typedef  typename Base<pfes>::B::FESpace FESpace;
    typedef typename FESpace::FElement FElement;
  //  typedef  FESpace FESpace;
    
    class Op : public E_F0info { public:  // passe de l'info ..
        typedef pfes * A;
        Expression a,b,c,d;
        // if c = 0 => a,b FESpace
        // if c != a FESpace et b,c KN_<double>
        static const int n_name_param =0;
     //   static basicAC_F0::name_and_type name_param[] ;
        Expression nargs[n_name_param];
        bool arg(int i,Stack stack,bool a) const{ return nargs[i] ? GetAny<bool>( (*nargs[i])(stack) ): a;}
        long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
        KN_<long>  arg(int i,Stack stack,KN_<long> a ) const{ return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
        
    public:
        Op(const basicAC_F0 &  args,Expression aa,Expression bb,Expression cc) : a(aa),b(bb),c(cc),d(0) {
          args.SetNameParam(n_name_param,0,nargs);
          
        }
        /*
        AnyType operator()(Stack stack)  const
          {
              if( verbosity>9) cout << " -- RestrictArray  "<< endl;
              pfes * pCUh = GetAny< pfes * >((* a)(stack));
              pfes * pFVh = GetAny<  pfes * >((* b)(stack));
              // verif same  FE.
              KN_<long> ncf=  GetAny<  KN_<long>  >((* c)(stack));
              FESpace * pVCh = **pCUh;
              FESpace * pVFh = **pFVh;
              FESpace  VCh = *pVCh;
              FESpace  VFh = *pVFh;
              long neC = VCh.NbOfElements   ;
              long neF = VFh.NbOfElements   ;
              long ndfC = VCh.NbOfDF   ;
              long ndfF = VFh.NbOfDF   ;
              
              KN_<long> nc2f= ncf;
              
              KN<long> *pnum=new KN<long>(ndfC);
              KN<long> inj=*pnum;
              inj = -1; // un set ..
              if( verbosity>9) cout<< " ne =" << neC << " " << neF << endl;
              
              for(int kc=0; kc <VCh.NbOfElements; kc++)
              {
                  
                  int kf = nc2f(kc);
                  FElement KC(pVCh,kc);
                  FElement KF(pVFh,kf);
                  
                  int ndofKC = KC.NbDoF() ;
                  int ndofKF =  KF.NbDoF() ;
                  if( verbosity>99) cout << kc << " " << kf << " : " <<ndofKC << " " << ndofKF << endl;
                  ffassert(ndofKC== ndofKF );
                  ffassert( kf >= 0 && kf < neF);
                  ffassert( kc >= 0 && kc< neC);
 
                  for(int df=0; df<ndofKC; df++)
                  {
                      int dfC =KC(df), dfF=KF(df);
                      if( verbosity>99) cout << dfC <<" ->  "<< dfF << endl;
                      assert(dfC >= 0 && dfC < ndfC);
                      inj[dfC] = dfF;
                  }
                 
              }
              
              if( verbosity>9) cout << " restrict:  Inject= " << inj << endl;
               ffassert(inj.min() != -1);
            //Add2StackOfPtr2FreeRC(stack,pnum);
              
            return pnum;
          }
	*/
    };
    RestrictArray() : OneOperator(  atype<const typename RestrictArray<pfes>::Op *>(),//atype<KN<long>* >(),
                                        atype<pfes *>(),
                                        atype<pfes *>(),
                                        atype<KN_<long> >()) {}
    
    E_F0 * code(const basicAC_F0 & args) const
    {
        if(args.size()!=3)CompileError("Bug in RestrictArray code nb of args !=  3 ????? bizarre" );
      return  new Op(args,t[0]->CastTo(args[0]),
                           t[1]->CastTo(args[1]),
                           t[2]->CastTo(args[2]));
    }
};
// end restrict ...
template< typename T > struct Base { typedef  T B; };
template< typename T > struct Base< T* >{ typedef T B;};

template<class pfes, int INIT>
AnyType SetRestrict(Stack stack,Expression einj,Expression erest)
{
    
    typedef  typename Base<pfes>::B::FESpace FESpace;
    typedef typename FESpace::FElement FElement;

   KN<long>  * pinj =GetAny<KN<long>*>((*einj)(stack));
   const typename RestrictArray<pfes>::Op * ar(dynamic_cast<const typename RestrictArray<pfes>::Op *>(erest));
    ffassert(ar);
    if( verbosity>9) cout << " -- RestrictArray  "<< endl;
    pfes * pCUh = GetAny< pfes * >((* ar->a)(stack));
    pfes * pFVh = GetAny<  pfes * >((* ar->b)(stack));
    // verif same  FE.
    KN_<long> ncf=  GetAny<  KN_<long>  >((* ar->c)(stack));
    FESpace * pVCh = **pCUh;
    FESpace * pVFh = **pFVh;
    FESpace & VCh = *pVCh;
    FESpace & VFh = *pVFh;
    long neC = VCh.NbOfElements   ;
    long neF = VFh.NbOfElements   ;
    long ndfC = VCh.NbOfDF   ;
    long ndfF = VFh.NbOfDF   ;
    
    KN_<long> nc2f= ncf;
    if(INIT==0)
        pinj->init(ndfC);
    else pinj->resize(ndfC);
  //  KN<long> *pnum= INIT ==0 ? new KN<long>(ndfC): pnum;
    KN<long> & inj=*pinj;
    inj = -1; // un set ..
    if( verbosity>9) cout<< " ne =" << neC << " " << neF << endl;
    
    for(int kc=0; kc <VCh.NbOfElements; kc++)
    {
        
        int kf = nc2f(kc);
        FElement KC(pVCh,kc);
        FElement KF(pVFh,kf);
        
        int ndofKC = KC.NbDoF() ;
        int ndofKF =  KF.NbDoF() ;
        if( verbosity>99) cout << kc << " " << kf << " : " <<ndofKC << " " << ndofKF << endl;
        ffassert(ndofKC== ndofKF );
        ffassert( kf >= 0 && kf < neF);
        ffassert( kc >= 0 && kc< neC);
        
        for(int df=0; df<ndofKC; df++)
        {
            int dfC =KC(df), dfF=KF(df);
            if( verbosity>99) cout << dfC <<" ->  "<< dfF << endl;
            assert(dfC >= 0 && dfC < ndfC);
            inj[dfC] = dfF;
        }
        
    }
    if( verbosity>9) cout << " restrict:  Inject= " << inj << endl;
    ffassert(inj.min() != -1);
   
  return pinj;
}

// Fin Add FH Houston --------
template<class pfes>
class MatrixInterpolation : public OneOperator { public:

    class Op : public E_F0info { public:
       typedef pfes * A;
       Expression a,b,c,d;
       // if c = 0 => a,b FESpace
       // if c != a FESpace et b,c KN_<double>
       static const int n_name_param =5;
       static basicAC_F0::name_and_type name_param[] ;
        Expression nargs[n_name_param];
     bool arg(int i,Stack stack,bool a) const{ return nargs[i] ? GetAny<bool>( (*nargs[i])(stack) ): a;}
     long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
     KN_<long>  arg(int i,Stack stack,KN_<long> a ) const{ return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}

       public:
       Op(const basicAC_F0 &  args,Expression aa,Expression bb) : a(aa),b(bb),c(0),d(0) {
         args.SetNameParam(n_name_param,name_param,nargs);  }
       Op(const basicAC_F0 &  args,Expression aa,Expression bb,Expression cc) : a(aa),b(bb),c(cc),d(0) {
         args.SetNameParam(n_name_param,name_param,nargs); } 
	Op(const basicAC_F0 &  args,Expression aa,Expression bb,Expression cc,Expression dd) : a(aa),b(bb),c(cc),d(dd) {
	    args.SetNameParam(n_name_param,name_param,nargs); } 
	
	
    };
    // interpolation(Vh,Vh)
   MatrixInterpolation() : OneOperator(atype<const typename MatrixInterpolation<pfes>::Op *>(),
                                       atype<pfes *>(),
                                       atype<pfes *>()) {}
   // interpolation(Vh,xx,yy) // 2d 
   MatrixInterpolation(int bidon) : OneOperator(atype<const typename MatrixInterpolation<pfes>::Op *>(),
                                                atype<pfes *>(),atype<KN_<double> >(),atype<KN_<double> >()) {}
 
    // interpolation(Vh,xx,yy,zz) // 3d 
    MatrixInterpolation(int bidon,int bidon2) : OneOperator(atype<const typename MatrixInterpolation<pfes>::Op *>(),
						 atype<pfes *>(),atype<KN_<double> >(),atype<KN_<double> >(),atype<KN_<double> >()) {}
    
  
    E_F0 * code(const basicAC_F0 & args) const 
     { 
       if(args.size()==2)
       return  new Op(args,t[0]->CastTo(args[0]),
                           t[1]->CastTo(args[1]));
       else if(args.size()==3)
       return  new Op(args,t[0]->CastTo(args[0]),
                           t[1]->CastTo(args[1]),
                           t[2]->CastTo(args[2]));
       else if(args.size()==4)
	   return  new Op(args,t[0]->CastTo(args[0]),
			  t[1]->CastTo(args[1]),
			  t[2]->CastTo(args[2]),
			  t[2]->CastTo(args[3])
			  );
       else CompileError("Bug in MatrixInterpolation code nb != 2 or 3 ????? bizarre" );
       return 0;
     }
};
template<>
basicAC_F0::name_and_type  MatrixInterpolation<pfes>::Op::name_param[]= {
   {  "t", &typeid(bool)}, 
   {  "op", &typeid(long)},
   {  "inside",&typeid(bool)},
   {  "composante",&typeid(long)},
   {  "U2Vc",&typeid(KN_<long>)}

};

template<>
basicAC_F0::name_and_type  MatrixInterpolation<pfes3>::Op::name_param[]= {
    {  "t", &typeid(bool)}, 
    {  "op", &typeid(long)},
    {  "inside",&typeid(bool)},
    {  "composante",&typeid(long)},
    {  "U2Vc",&typeid(KN_<long>)}
    
};




template<class R>
   class SetMatrix_Op : public E_F0mps { public:
       Expression a; 
       
       static  aType btype;
       static const int n_name_param =NB_NAME_PARM_MAT; //  add nbiter FH 30/01/2007 11 -> 12  //add var MUMPS+autre
       static basicAC_F0::name_and_type name_param[] ;
       Expression nargs[n_name_param];
       const OneOperator * precon;
       bool arg(int i,Stack stack,bool a) const{ return nargs[i] ? GetAny<bool>( (*nargs[i])(stack) ): a;}
       long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}

       public:
       SetMatrix_Op(const basicAC_F0 &  args,Expression aa) : a(aa) {
         args.SetNameParam(n_name_param,name_param,nargs);
         precon = 0; //  a changer 
         if ( nargs[3])
          {
           const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(nargs[3]);
           assert(op);
           precon = op->Find("(",ArrayOfaType(atype<KN<R>* >(),false)); // strange bug in g++ is R become a double
          }
         
       } 
       AnyType operator()(Stack stack)  const ;
    };


template<class R>
class SetMatrix : public OneOperator { public:  

 
  // SetMatrix() : OneOperator(atype<const typename SetMatrix<R>::Op *>(),atype<Matrice_Creuse<R> *>() ) {}
   SetMatrix() : OneOperator(SetMatrix_Op<R>::btype,atype<Matrice_Creuse<R> *>() ) {}
  
    E_F0 * code(const basicAC_F0 & args) const 
     { 
       return  new SetMatrix_Op<R>(args,t[0]->CastTo(args[0])); 
     }
};

template<class R>
       aType  SetMatrix_Op<R>::btype=0;

template <class R> 
basicAC_F0::name_and_type  SetMatrix_Op<R>::name_param[]= {
 LIST_NAME_PARM_MAT



/*
   {   "paramint",&typeid(KN<int>)}, // Add J. Morice 02/09 
   {   "paramdouble",&typeid(KN<double>)},

   {   "paramstring",&typeid(string*)},
   {   "permrow",&typeid(KN_<long>)},
   {   "permcol",&typeid(KN_<long>)},
   {   "fileparamint",&typeid(string*)}, // Add J. Morice 02/09 
   {   "fileparamdouble",&typeid(string*)},
   {   "fileparamstring",&typeid(string* )},
   {   "filepermrow",&typeid(string*)},
   {   "filepermcol",&typeid(string*)} //22
 */
};

template<class R>
AnyType SetMatrix_Op<R>::operator()(Stack stack)  const 
{
   Matrice_Creuse<R> *  A= GetAny<Matrice_Creuse<R> *>((*a)(stack));
   
    ffassert(A);
    if( !A->A) A->A.master(new MatriceMorse<R>());//  set to empty matrix .. mars 2014 FH ..
    Data_Sparse_Solver ds;
    bool VF=false;
   // bool factorize=false;
    ds.factorize=false;
    /*
  long NbSpace = 50; 
  long itmax=0; 
  double epsilon=1e-6;
//  bool VF=false;
//  VF=isVF(op->largs);
 // assert(!VF); 
  double tgv = 1e30;
  double tol_pivot=-1;
  double tol_pivot_sym=-1;
  
  KN<int>  param_int;
  KN<double> param_double; 
  string *param_char=NULL;
  KN<int> perm_r; 
  KN<int> perm_c;
  string *file_param_int;  // Add J. Morice 02/09 
  string *file_param_double; 
  string* file_param_char;
  string* file_param_perm_r;
  string* file_param_perm_c;  
*/
// type de matrice par default 
  TypeSolveMat tmat= TypeSolveMat::defaultvalue; 
  if(   tmat !=  TypeSolveMat::SparseSolver  )    
    tmat=TypeSolveMat::GMRES;
 ds.typemat=&tmat; 
 SetEnd_Data_Sparse_Solver<R>(stack,ds,nargs,n_name_param);  
/*     
  TypeSolveMat    *typemat=&tmat;
  bool initmat=true;
  int strategy=0; 
   
  if (nargs[0]) ds.initmat= ! GetAny<bool>((*nargs[0])(stack));
  if (nargs[1]) ds.typemat= GetAny<TypeSolveMat *>((*nargs[1])(stack));
  if (nargs[2]) ds.epsilon= GetAny<double>((*nargs[2])(stack));
  // 3 precon 
  if (nargs[4]) ds.NbSpace= GetAny<long>((*nargs[4])(stack));
  if (nargs[6]) ds.tgv= GetAny<double>((*nargs[6])(stack));
  if (nargs[7]) ds.factorize= GetAny<bool>((*nargs[7])(stack));
  
  if (nargs[8]) ds.strategy = GetAny<long>((*nargs[8])(stack)); 
  if (nargs[9]) ds.tol_pivot = GetAny<double>((*nargs[9])(stack)); 
  if (nargs[10]) ds.tol_pivot_sym = GetAny<double>((*nargs[10])(stack)); 
  if (nargs[11]) ds.itmax = GetAny<long>((*nargs[11])(stack)); //  frev 2007 OK
   

  if (nargs[12]) ds.param_int= GetAny< KN<int> >((*nargs[12])(stack));  // Add J. Morice 02/09 
  if (nargs[13]) ds.param_double= GetAny< KN<double> >((*nargs[13])(stack));
  if (nargs[14]) ds.param_char= GetAny< string * >((*nargs[14])(stack));  //
  if (nargs[15]) ds.perm_r = GetAny< KN<int > >((*nargs[15])(stack));
  if (nargs[16]) ds.perm_c = GetAny< KN<int> >((*nargs[16])(stack));  //
  if (nargs[17]) ds.file_param_int= GetAny< string* >((*nargs[17])(stack));  // Add J. Morice 02/09 
  if (nargs[18]) ds.file_param_double= GetAny< string* >((*nargs[18])(stack));
  if (nargs[19]) ds.file_param_char= GetAny< string* >((*nargs[19])(stack));  //
  if (nargs[20]) ds.file_param_perm_r = GetAny< string* >((*nargs[20])(stack));
  if (nargs[21]) ds.file_param_perm_c = GetAny< string* >((*nargs[21])(stack));  //
*/


   if(A->typemat.profile != ds.typemat->profile) 
   {
     cerr << " type of matrix " << A->typemat<<endl;
     cerr << " type of matrix for solver " <<*ds.typemat<<endl;
     
     ExecError(" Set incompatibility between solver and type of matrix");
   }
  if( ds.factorize ) {
    MatriceProfile<R> * pf = dynamic_cast<MatriceProfile<R> *>((MatriceCreuse<R> *) A->A);
    assert(pf);
    switch (ds.typemat->t) {
    case TypeSolveMat::LU: pf->LU(Abs(ds.epsilon));break;
    case TypeSolveMat::CROUT: pf->crout(Abs(ds.epsilon));break;
    case TypeSolveMat::CHOLESKY: pf->cholesky(Abs(ds.epsilon));break;
    default: ExecError("Sorry no factorization for this type for matrix"); 
    }
    
  }    
  SetSolver<R>(stack,VF,*A->A,ds) ;/*stack,*A->A,typemat,VF,epsilon,NbSpace,itmax,precon,strategy,tgv,tol_pivot,tol_pivot_sym, 
	       param_int, param_double, param_char, perm_r, perm_c, file_param_int, file_param_double, 
	       file_param_char, file_param_perm_r, file_param_perm_c);*/

  return Nothing; 
}



bool SetDefaultSolver()
{

#ifdef HAVE_LIBUMFPACKXXXXXX
    if(verbosity>1)
	cout << " SetDefault sparse solver to UMFPACK" << endl;
    DefSparseSolver<double>::solver  =BuildSolverUMFPack;
    DefSparseSolver<Complex>::solver =BuildSolverUMFPack;
    DefSparseSolverSym<double>::solver  =BuildSolverGMRES<double>;
    DefSparseSolverSym<Complex>::solver =BuildSolverGMRES<Complex>;
    
#else
    if(verbosity>1)
	cout << " SetDefault sparse solver to GMRES (no UMFPACK)" << endl;
    DefSparseSolver<double>::solver  =BuildSolverGMRES<double>;
    DefSparseSolver<Complex>::solver =BuildSolverGMRES<Complex>;
    DefSparseSolverSym<double>::solver  =BuildSolverGMRES<double>;
    DefSparseSolverSym<Complex>::solver =BuildSolverGMRES<Complex>;
    
#endif
    return true;

}

template<int init>
AnyType SetMatrixInterpolation(Stack,Expression ,Expression);
template<int init>
AnyType SetMatrixInterpolation3(Stack,Expression ,Expression);

//------

template<class R>
void BuildCombMat(map< pair<int,int>, R> & mij,const KNM_<R> & A, int ii00=0,int jj00=0,R coef=R(1.),bool cnj=false)
{
  double eps0=numeric_limits<double>::min();
  int i,j;
  int n = A.N(),m=A.M();
  for ( i=0;i<n;i++)
   for ( j=0;j<m;j++)
          {
           R cij=coef*A(i,j);
           if (cnj)  cij = RNM::conj(cij); 
           if(Fem2D::norm(cij) >eps0)
             mij[ij_mat(false,ii00,jj00,i,j)] += cij;
         
   }

}

void buildInterpolationMatrix(MatriceMorse<R> * m,const FESpace & Uh,const FESpace & Vh,void *data)
{  //  Uh = Vh 

  int op=op_id; //  value of the function
  bool transpose=false;
  bool inside=false;
   int * iU2V=0;  
  if (data)
   {
     int * idata=static_cast<int*>(data);
     transpose=idata[0];
     op=idata[1];
     inside=idata[2];
       iU2V= idata + 5;
     ffassert(op>=0 && op < 4);
   }
  if(verbosity>2) 
    {
      cout << "  -- buildInterpolationMatrix   transpose =" << transpose << endl
           << "              value, dx , dy          op = " << op << endl
           << "                            just  inside = " << inside << endl;
    }
  using namespace Fem2D;
  int n=Uh.NbOfDF;
  int mm=Vh.NbOfDF;
  if(transpose) Exchange(n,mm);
  m->symetrique = false;
  m->dummy=false;
  m->a=0;
  m->lg=0;
  m->cl=0;
  m->nbcoef=0;
  m->n=n;
  m->m=mm;
  int n1=n+1;
  const  Mesh & ThU =Uh.Th; // line 
  const  Mesh & ThV =Vh.Th; // colunm
  bool samemesh =  &Uh.Th == &Vh.Th;  // same Mesh
  int thecolor =0;
  //  int nbn_u = Uh.NbOfNodes;
  //int nbn_v = Vh.NbOfNodes;
  
  int nbcoef =0;
  
  KN<int> color(ThV.nt);
  KN<int> mark(n);
  mark=0;
  
  int *lg = new int [n1];
  int * cl = 0;
  double *a=0;
  
  color=thecolor++;
  FElement Uh0 = Uh[0];
  FElement Vh0 = Vh[0];
  
  FElement::aIPJ ipjU(Uh0.Pi_h_ipj()); 
  FElement::aR2  PtHatU(Uh0.Pi_h_R2()); 
  
  //  FElement::aIPJ ipjV(Vh0.Pi_h_ipj()); 
  // FElement::aR2  PtHatV(Vh0.Pi_h_R2()); 
  
  int nbdfVK= Vh0.NbDoF();
  //  int nbdfUK= Uh0.NbDoF();
  int NVh= Vh0.N;
  //int NUh= Uh0.N;
  
  //ffassert(NVh==NUh); 
  
  
  int nbp= PtHatU.N(); // 
  //  int nbc= ipjU.N(); // 
  KN<R2> PV(nbp);   //  the PtHat in ThV mesh
  KN<int> itV(nbp); // the Triangle number
  KN<bool> intV(nbp); // ouside or not 
  KN<R>   AipjU(ipjU.N());
  
  KNM<R> aaa(nbp,nbdfVK); 


   const R eps = 1.0e-10;
   const int sfb1=Vh0.N*last_operatortype*Vh0.NbDoF();
   KN<R> kv(sfb1*nbp);
   R * v = kv;
    KN<int> ik(nbp); // the Triangle number
//   KNMK_<R> fb(v+ik[i],VhV0.NbDoF,VhV0.N,1); //  the value for basic fonction

   bool whatd[last_operatortype];
   for (int i=0;i<last_operatortype;i++) 
     whatd[i]=false;
   whatd[op]=true; // the value of function
   KN<bool> fait(Uh.NbOfDF);
   fait=false;
  map< pair<int,int> , double > sij; 
  
  for (int step=0;step<2;step++) 
   {
      
	  for (int it=0;it<ThU.nt;it++)
	    {
	      thecolor++; //  change the current color
	      const Triangle & TU(ThU[it]);
	      FElement KU(Uh[it]);
	      KU.Pi_h(AipjU);
	      int nbkU = 0;
	      if (samemesh)
	        {
	          nbkU = 1; 
	          PV = PtHatU;
	          itV = it;
		  intV= false;// add July 2009 (unset varaible FH)
	        }
	      else 
	       {  
	          const Triangle *ts=0;
	          for (int i=0;i<nbp;i++)
	            {
	              bool outside;
	              ts=ThV.Find(TU(PtHatU[i]),PV[i],outside,ts); 
		      if(outside && verbosity>9 ) 
		        cout << it << " " << i << " :: " << TU(PtHatU[i]) << "  -- "<< outside << PV[i] << " " << ThV(ts) << " ->  " <<  (*ts)(PV[i]) <<endl;
	              itV[i]= ThV(ts);
	              intV[i]=outside && inside; //  ouside and inside flag 
	            }
	       }
	       
	      
	      for (int p=0;p<nbp;p++)
	         { 

	              KNMK_<R> fb(v+p*sfb1,nbdfVK,NVh,last_operatortype); // valeur de fonction de base de Vh 
	              // ou:   fb(idf,j,0) valeur de la j composante de la fonction idf 
	              Vh0.tfe->FB(whatd,ThV,ThV[itV[p]],PV[p],fb);  
	          }

	      for (int i=0;i<ipjU.N();i++) 
	          { // pour tous le terme 
	           const FElement::IPJ &ipj_i(ipjU[i]);
	          // assert(ipj_i.j==0); // car  Vh.N=0
	           int dfu = KU(ipj_i.i); // le numero de df global 
	           if(fait[dfu]) continue;
	           int jU = ipj_i.j; // la composante dans U
	           int p=ipj_i.p;  //  le points
	           if (intV[p]) continue; //  ouside and inside flag => next 
	           R aipj = AipjU[i];
	           FElement KV(Vh[itV[p]]);
		      int jV=jU;
		      if(iU2V) jV=iU2V[jU];
		    
		    if(jV>=0 && jV<NVh)
			{
			    KNMK_<R> fb(v+p*sfb1,nbdfVK,NVh,last_operatortype); 
			    KN_<R> fbj(fb('.',jV,op)); 
			    
			    for (int idfv=0;idfv<nbdfVK;idfv++) 
				if (Abs(fbj[idfv])>eps) 
				  {
				      int dfv=KV(idfv);
				      int ii=dfu, jj=dfv;
				      if(transpose) Exchange(ii,jj);
				      // le term dfu,dfv existe dans la matrice
				      R c= fbj[idfv]*aipj;
				      //  cout << " Mat inter " << i << " , "<< j << " = " << c << " " <<step << " " << it << " " <<  endl; 
				      if(Abs(c)>eps)
					  //
					  sij[make_pair(ii,jj)] += c;// correct du to periodic BC. June 2011
				      /*   
				       if(step==0)
				       sij.insert(make_pair(i,j));
				       else	                  
				       (*m)(i,j)=c;
				       */
				  }
			}
	                      
	          }
	          
	      for (int df=0;df<KU.NbDoF();df++) 
	          {  
	           int dfu = KU(df); // le numero de df global 
	           fait[dfu]=true;
	           }
	       
	       
	    }
	    if (step==0)
	     {
	         nbcoef = sij.size();
	         cl = new int[nbcoef];
	         a = new double[nbcoef];
	         int k=0;
	         for(int i=0;i<n1;i++)
	           lg[i]=0;
	          
	         for (map<pair<int,int>, double >::iterator kk=sij.begin();kk!=sij.end();++kk)
	          { 
	            int i= kk->first.first;
	            int j= kk->first.second;
	           // cout << " Mat inter " << i << " , "<< j  << endl;
	            cl[k]=j;
	            a[k]= kk->second;	            
	            lg[i+1]=++k;
	          }
	         assert(k==nbcoef);
	         //  on bouche les ligne vide   lg[i]=0; 
	         //  lg est un tableau croissant  =>
	         for(int i=1;i<n1;i++)
	              lg[i]=max(lg[i-1],lg[i]) ;
	         m->lg=lg;
             m->cl=cl;
             m->a=a;
             m->nbcoef=nbcoef;
             fait=false;
	     }
    }
    sij.clear();
  //assert(0); // a faire to do
}

void buildInterpolationMatrix(MatriceMorse<R> * m,const FESpace3 & Uh,const FESpace3 & Vh,void *data)
{  //  Uh = Vh 
    typedef FESpace3::Mesh Mesh;
    typedef FESpace3::FElement  FElement;
    typedef Mesh::Element  Element;
    typedef FESpace3::Rd  Rd;
    typedef Element::RdHat  RdHat;
    
    int op=op_id; //  value of the function
    bool transpose=false;
    bool inside=false;
    int * iU2V=0;  
    if (data)
      {
	int * idata=static_cast<int*>(data);
	transpose=idata[0];
	op=idata[1];
	inside=idata[2];
	iU2V= idata + 5;
	ffassert(op>=0 && op < 4);
      }
    if(verbosity>2) 
      {
	cout << "  -- buildInterpolationMatrix   transpose =" << transpose << endl
	<< "              value, dx , dy          op = " << op << endl
	<< "                            just  inside = " << inside << endl;
      }
    using namespace Fem2D;
    int n=Uh.NbOfDF;
    int mm=Vh.NbOfDF;
    if(transpose) Exchange(n,mm);
    m->symetrique = false;
    m->dummy=false;
    m->a=0;
    m->lg=0;
    m->cl=0;
    m->nbcoef=0;
    m->n=n;
    m->m=mm;
    int n1=n+1;
    const  Mesh & ThU =Uh.Th; // line 
    const  Mesh & ThV =Vh.Th; // colunm
    bool samemesh =  &Uh.Th == &Vh.Th;  // same Mesh
    int thecolor =0;
    //  int nbn_u = Uh.NbOfNodes;
    //int nbn_v = Vh.NbOfNodes;
    
    int nbcoef =0;
    
    KN<int> color(ThV.nt);
    KN<int> mark(n);
    mark=0;
    
    int *lg = new int [n1];
    int * cl = 0;
    double *a=0;
    
    color=thecolor++;
    FElement Uh0 = Uh[0];
    FElement Vh0 = Vh[0];
    
   // FElement::aIPJ ipjU(Uh0.Pi_h_ipj()); 
   // FElement::aR2  PtHatU(Uh0.Pi_h_R2()); 
    
    //  FElement::aIPJ ipjV(Vh0.Pi_h_ipj()); 
    // FElement::aR2  PtHatV(Vh0.Pi_h_R2()); 
    
    int nbdfVK= Vh0.NbDoF();
    //  int nbdfUK= Uh0.NbDoF();
    int NVh= Vh0.N;
   // int NUh= Uh0.N;
    
    //ffassert(NVh==NUh); 
    
    InterpolationMatrix<RdHat> ipmat(Uh);    
  
    int nbp=ipmat.np; // 
    //  int nbc= ipjU.N(); // 
    KN<RdHat> PV(nbp);   //  the PtHat in ThV mesh
    KN<int> itV(nbp); // the Triangle number
    KN<bool> intV(nbp); // ouside or not 
   // KN<R>   AipjU(ipjU.N());
    
    KNM<R> aaa(nbp,nbdfVK); 
    
    
    const R eps = 1.0e-10;
    const int sfb1=Vh0.N*last_operatortype*Vh0.NbDoF();
    KN<R> kv(sfb1*nbp);
    R * v = kv;
    KN<int> ik(nbp); // the Triangle number
    //   KNMK_<R> fb(v+ik[i],VhV0.NbDoF,VhV0.N,1); //  the value for basic fonction
    op= op==3 ? op_dz : op; //   renumber op ????  dec 2010 FH. 
    What_d whatd= 1<< op;
    KN<bool> fait(Uh.NbOfDF);
    fait=false;
    map< pair<int,int> , double > sij; 
    for (int step=0;step<2;step++) 
      {
	
	for (int it=0;it<ThU.nt;it++)
	  {
	    thecolor++; //  change the current color
	    const Element & TU(ThU[it]);
	    FElement KU(Uh[it]);
	    ipmat.set(KU);
	    //KU.Pi_h(AipjU);
	    int nbkU = 0;
	    if (samemesh)
	      {
		nbkU = 1; 
		PV = ipmat.P;
		itV = it;
		intV= false;// add July 2009 (unset varaible FH)
	      }
	    else 
	      {  
	          const Element *ts=0;
	          for (int i=0;i<nbp;i++)
	            {
	              bool outside;
	              ts=ThV.Find(TU(ipmat.P[i]),PV[i],outside,ts); 
		      if(outside && verbosity>9 ) 
			  cout << it << " " << i << " :: " << TU(ipmat.P[i]) << "  -- "<< outside << PV[i] << " " << ThV(ts) << " ->  " <<  (*ts)(PV[i]) <<endl;
	              itV[i]= ThV(ts);
	              intV[i]=outside && inside; //  ouside and inside flag 
	            }
	      }
	    
	    
	    for (int p=0;p<nbp;p++)
	      { 
		  
		  KNMK_<R> fb(v+p*sfb1,nbdfVK,NVh,last_operatortype); // valeur de fonction de base de Vh 
		  // ou:   fb(idf,j,0) valeur de la j composante de la fonction idf 
		  Vh0.tfe->FB(whatd,ThV,ThV[itV[p]],PV[p],fb);  
	      }
	    
	    for (int i=0;i<ipmat.ncoef;i++) 
	      { // pour tous le terme 
		  
	          // assert(ipj_i.j==0); // car  Vh.N=0
		  int dfu = KU(ipmat.dofe[i]); // le numero de df global 
		  if(fait[dfu]) continue;
		  int jU = ipmat.comp[i]; // la composante dans U
		  int p=ipmat.p[i];  //  le point
		  if (intV[p]) continue; //  ouside and inside flag => next 
		  R aipj = ipmat.coef[i];
		  FElement KV(Vh[itV[p]]);
		  int jV=jU;
		  if(iU2V) jV=iU2V[jU];
		  
		  if(jV>=0 && jV<NVh)
		    {
		      KNMK_<R> fb(v+p*sfb1,nbdfVK,NVh,last_operatortype); 
		      KN_<R> fbj(fb('.',jV,op)); 
		      
		      for (int idfv=0;idfv<nbdfVK;idfv++) 
			  if (Abs(fbj[idfv])>eps) 
			    {
			      int dfv=KV(idfv);
			      int ii=dfu, jj=dfv;
			      if(transpose) Exchange(ii,jj);
			      // le term dfu,dfv existe dans la matrice
			      R c= fbj[idfv]*aipj;
			      //  cout << " Mat inter " << i << " , "<< j << " = " << c << " " <<step << " " << it << " " <<  endl; 
			      if(Abs(c)>eps)
				  //
				  sij[make_pair(ii,jj)] += c;// correct du to periodic BC. June 2011
			      /*   
			       if(step==0)
			       sij.insert(make_pair(i,j));
			       else	                  
			       (*m)(i,j)=c;
			       */
			    }
		    }
		  
	      }
	    
	    for (int df=0;df<KU.NbDoF();df++) 
	      {  
		  int dfu = KU(df); // le numero de df global 
		  fait[dfu]=true;
	      }
	    
	    
	  }
	if (step==0)
	  {
	    nbcoef = sij.size();
	    cl = new int[nbcoef];
	    a = new double[nbcoef];
	    int k=0;
	    for(int i=0;i<n1;i++)
		lg[i]=0;
	    
	    for (map<pair<int,int>, double >::iterator kk=sij.begin();kk!=sij.end();++kk)
	      { 
		  int i= kk->first.first;
		  int j= kk->first.second;
		  // cout << " Mat inter " << i << " , "<< j  << endl;
		  cl[k]=j;
		  a[k]= kk->second;	            
		  lg[i+1]=++k;
	      }
	    assert(k==nbcoef);
	    //  on bouche les ligne vide   lg[i]=0; 
	    //  lg est un tableau croissant  =>
	    for(int i=1;i<n1;i++)
		lg[i]=max(lg[i-1],lg[i]) ;
	    m->lg=lg;
	    m->cl=cl;
	    m->a=a;
	    m->nbcoef=nbcoef;
	    fait=false;
	  }
      }
    sij.clear();
    //assert(0); // a faire to do
}


MatriceMorse<R> *  buildInterpolationMatrix1(const FESpace & Uh,const KN_<double> & xx,const KN_<double> & yy ,int *data)
{  //  Uh = Vh 

  int op=op_id; //  value of the function
  int icomp=0;
  bool transpose=false;
  bool inside=false;
  if (data)
   {
     transpose=data[0];
     op=data[1];
     inside=data[2];
     icomp = data[3];
     ffassert(op>=0 && op < 4);
     
   }
  if(verbosity>2) 
    {
      cout << "  -- buildInterpolationMatrix   transpose =" << transpose << endl
           << "              value, dx , dy          op = " << op << endl
           << "              composante                 = " << icomp << endl
           << "                            just  inside = " << inside << endl;
    }
  using namespace Fem2D;
  int n=Uh.NbOfDF;
  int mm=xx.N();
  int nbxx= mm;
  if(transpose) Exchange(n,mm);
  const  Mesh & ThU =Uh.Th; // line 
  
   
  FElement Uh0 = Uh[0];
  
   
   
 int nbdfUK= Uh0.NbDoF();
 int NUh= Uh0.N;
 
 ffassert(icomp < NUh && icomp >=0); 
    

   const int sfb1=Uh0.N*last_operatortype*Uh0.NbDoF();
   KN<R> kv(sfb1);
   R * v = kv;
   const R eps = 1.0e-10;

   bool whatd[last_operatortype];
   for (int i=0;i<last_operatortype;i++) 
     whatd[i]=false;
   whatd[op]=true; // the value of function
   KN<bool> fait(Uh.NbOfDF);
   fait=false;
   map< pair<int,int> , double > sij; 
   R2 Phat;
   bool outside;
   
   for(int ii=0;ii<nbxx;ii++)
   {
     const Triangle *ts=ThU.Find(R2(xx[ii],yy[ii]),Phat,outside); 
     if(outside && !inside) continue;
     int it = ThU(ts); 
     FElement KU(Uh[it]);
     KNMK_<R> fb(v,nbdfUK,NUh,last_operatortype);
     Uh0.tfe->FB(whatd,ThU,ThU[it],Phat,fb);  
     KN_<R> Fwi(fb('.',icomp,op));
     for (int idfu=0;idfu<nbdfUK;idfu++) 
       {
        int  j = ii;
        int  i = KU(idfu);
        if(transpose) Exchange(i,j);
        R c = Fwi(idfu);
	if(Abs(c)>eps)
	  sij[make_pair(i,j)] += c;// correct du to periodic BC. June 2011
      }
      }
  
   MatriceMorse<R> * m = new MatriceMorse<R>(n,mm,sij,false);
    sij.clear();
   return m;
  //assert(0); // a faire to do
}

MatriceMorse<R> *  buildInterpolationMatrix1(const FESpace3 & Uh,const KN_<double> & xx,const KN_<double> & yy ,const KN_<double> & zz,int *data)
{  //  Uh = Vh 
    typedef FESpace3::Mesh Mesh;
    typedef FESpace3::FElement  FElement;
    typedef Mesh::Element  Element;
    typedef FESpace3::Rd  Rd;
    
    
    int op=op_id; //  value of the function
    int icomp=0;
    bool transpose=false;
    bool inside=false;
    if (data)
      {
	transpose=data[0];
	op=data[1];
	inside=data[2];
	icomp = data[3];
	ffassert(op>=0 && op < 4);
        if(op==3) op=op_dz;//  correct missing
      }
    if(verbosity>2) 
      {
	cout << "  -- buildInterpolationMatrix   transpose =" << transpose << endl
	<< "              value, dx , dy          op = " << op << endl
	<< "              composante                 = " << icomp << endl
	<< "                            just  inside = " << inside << endl;
      }
    using namespace Fem2D;
    int n=Uh.NbOfDF;
    int mm=xx.N();
    int nbxx= mm;
    if(transpose) Exchange(n,mm);
    const  Mesh & ThU =Uh.Th; // line 
    
    
    FElement Uh0 = Uh[0];
    
    
    
    int nbdfUK= Uh0.NbDoF();
    int NUh= Uh0.N;
    
    ffassert(icomp < NUh && icomp >=0); 
    
    
    const int sfb1=Uh0.N*last_operatortype*Uh0.NbDoF();
    KN<R> kv(sfb1);
    R * v = kv;
    const R eps = 1.0e-10;
    
    What_d whatd= 1 <<op;
    KN<bool> fait(Uh.NbOfDF);
    fait=false;
    map< pair<int,int> , double > sij; 
    Rd Phat;
    bool outside;
    
    for(int ii=0;ii<nbxx;ii++)
      {
	const Element *ts=ThU.Find(Rd(xx[ii],yy[ii],zz[ii]),Phat,outside); 
	if(outside && !inside) continue;
	int it = ThU(ts); 
	FElement KU(Uh[it]);
	KNMK_<R> fb(v,nbdfUK,NUh,last_operatortype);
	Uh0.tfe->FB(whatd,ThU,ThU[it],Phat,fb);  
	KN_<R> Fwi(fb('.',icomp,op));
	for (int idfu=0;idfu<nbdfUK;idfu++) 
	  {
	    int  j = ii;
	    int  i = KU(idfu);
	    if(transpose) Exchange(i,j);
	    R c = Fwi(idfu);
	    if(Abs(c)>eps)
		sij[make_pair(i,j)] += c;// correct du to periodic BC. June 2011
	  }
      }
    
    MatriceMorse<R> * m = new MatriceMorse<R>(n,mm,sij,false);
    sij.clear();
    return m;
    //assert(0); // a faire to do
}


AnyType SetMatrixInterpolation1(Stack stack,Expression emat,Expression einter,int init)
{
  using namespace Fem2D;
  
  Matrice_Creuse<R> * sparse_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
  const MatrixInterpolation<pfes>::Op * mi(dynamic_cast<const MatrixInterpolation<pfes>::Op *>(einter));
  ffassert(einter);
  int data[ MatrixInterpolation<pfes>::Op::n_name_param+100];
  data[0]=mi->arg(0,stack,false); // transpose not
  data[1]=mi->arg(1,stack,(long) op_id); ; // get just value
  data[2]=mi->arg(2,stack,false); ; // get just value
  data[3]=mi->arg(3,stack,0L); ; // get just value
  KN<long> U2Vc;
  U2Vc= mi->arg(4,stack,U2Vc); ;
  if( mi->c==0)
  { // old cas 
  pfes * pUh = GetAny< pfes * >((* mi->a)(stack));
  pfes * pVh = GetAny<  pfes * >((* mi->b)(stack));
  FESpace * Uh = **pUh;
  FESpace * Vh = **pVh;
  int NVh =Vh->N;
  int NUh =Uh->N;
  ffassert(NUh< 100-MatrixInterpolation<pfes>::Op::n_name_param);

      for(int i=0;i<NUh;++i)    
        data[5+i]=i;// 
      for(int i=0;i<min(NUh,(int) U2Vc.size());++i)    
	  data[5+i]= U2Vc[i];//
  if(verbosity>3)
	for(int i=0;i<NUh;++i)
	  {
	    cout << "The Uh componante " << i << " -> " << data[5+i] << "  Componante of Vh  " <<endl;
	  }
	  for(int i=0;i<NUh;++i)    
	if(data[5+i]>=NVh)
	  {
	      cout << "The Uh componante " << i << " -> " << data[5+i] << " >= " << NVh << " number of Vh Componante " <<endl;
	      ExecError("Interpolation incompability beetween componante ");
	  }
      
  ffassert(Vh);
  ffassert(Uh);
  
  //  sparse_mat->pUh=pUh;
  //sparse_mat->pVh=pVh;
  if(!init) sparse_mat->init();
  sparse_mat->typemat=TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
  sparse_mat->A.master(new MatriceMorse<R>(*Uh,*Vh,buildInterpolationMatrix,data));
  //  sparse_mat->A.master(new MatriceMorse<R>(*Uh,*Vh,buildInterpolationMatrix,data));
  }
  else 
  {  // new cas mars 2006
  pfes * pUh = GetAny< pfes * >((* mi->a)(stack));
  KN_<double>  xx = GetAny<  KN_<double>  >((* mi->b)(stack));
  KN_<double>  yy = GetAny<  KN_<double>  >((* mi->c)(stack));
  ffassert( xx.N() == yy.N()); 
  FESpace * Uh = **pUh;
  ffassert(Uh);
  
  //  sparse_mat->pUh=0;
  //  sparse_mat->pVh=0;
  if(!init) sparse_mat->init();
  sparse_mat->typemat=TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
  sparse_mat->A.master(buildInterpolationMatrix1(*Uh,xx,yy,data));
  }
   return sparse_mat; // Warning .. no correct gestion of temp ptr ..
 // return Add2StackOfPtr2Free(stack,sparse_mat);
}

AnyType SetMatrixInterpolation31(Stack stack,Expression emat,Expression einter,int init)
{
    using namespace Fem2D;
    
    Matrice_Creuse<R> * sparse_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
    const MatrixInterpolation<pfes3>::Op * mi(dynamic_cast<const MatrixInterpolation<pfes3>::Op *>(einter));
    ffassert(einter);
    int data[ MatrixInterpolation<pfes3>::Op::n_name_param+100];
    data[0]=mi->arg(0,stack,false); // transpose not
    data[1]=mi->arg(1,stack,(long) op_id); ; // get just value
    data[2]=mi->arg(2,stack,false); ; // get just value
    data[3]=mi->arg(3,stack,0L); ; // get just value
    KN<long> U2Vc;
    U2Vc= mi->arg(4,stack,U2Vc); ;
    if( mi->c==0)
      { // old cas 
	  pfes3 * pUh = GetAny< pfes3 * >((* mi->a)(stack));
	  pfes3 * pVh = GetAny<  pfes3 * >((* mi->b)(stack));
	  FESpace3 * Uh = **pUh;
	  FESpace3 * Vh = **pVh;
	  int NVh =Vh->N;
	  int NUh =Uh->N;
	  ffassert(NUh< 100-MatrixInterpolation<pfes3>::Op::n_name_param);
	  
	  for(int i=0;i<NUh;++i)    
	      data[5+i]=i;// 
	  for(int i=0;i<min(NUh,(int) U2Vc.size());++i)    
	      data[5+i]= U2Vc[i];//
	  if(verbosity>3)
	      for(int i=0;i<NUh;++i)
		{
		  cout << "The Uh componante " << i << " -> " << data[5+i] << "  Componante of Vh  " <<endl;
		}
	  for(int i=0;i<NUh;++i)    
	      if(data[5+i]>=NVh)
		{
		  cout << "The Uh componante " << i << " -> " << data[5+i] << " >= " << NVh << " number of Vh Componante " <<endl;
		  ExecError("Interpolation incompability beetween componante ");
		}
	  
	  ffassert(Vh);
	  ffassert(Uh);
	  if(!init) sparse_mat->init();
	  //  sparse_mat->pUh=pUh;
	  //sparse_mat->pVh=pVh;
	  sparse_mat->typemat=TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
	  sparse_mat->A.master(new MatriceMorse<R>(*Uh,*Vh,buildInterpolationMatrix,data));
	  //  sparse_mat->A.master(new MatriceMorse<R>(*Uh,*Vh,buildInterpolationMatrix,data));
      }
    else 
      {  // new cas mars 2006
	  pfes3 * pUh = GetAny< pfes3 * >((* mi->a)(stack));
	  KN_<double>  xx = GetAny<  KN_<double>  >((* mi->b)(stack));
	  KN_<double>  yy = GetAny<  KN_<double>  >((* mi->c)(stack));
	  KN_<double>  zz = GetAny<  KN_<double>  >((* mi->d)(stack));
	  ffassert( xx.N() == yy.N()); 
	  ffassert( xx.N() == zz.N()); 
	  FESpace3 * Uh = **pUh;
	  ffassert(Uh);
	  if(!init) sparse_mat->init();
	  //  sparse_mat->pUh=0;
	  //  sparse_mat->pVh=0;
	  sparse_mat->typemat=TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
	  sparse_mat->A.master(buildInterpolationMatrix1(*Uh,xx,yy,zz,data));
      }
    return sparse_mat;
}

template<int init>
AnyType SetMatrixInterpolation(Stack stack,Expression emat,Expression einter)
{ return SetMatrixInterpolation1(stack,emat,einter,init);}
template<int init>
AnyType SetMatrixInterpolation3(Stack stack,Expression emat,Expression einter)
{ return SetMatrixInterpolation31(stack,emat,einter,init);}


template<class RA,class RB,class RAB,int init>
AnyType ProdMat(Stack stack,Expression emat,Expression prodmat)
{
  using namespace Fem2D;
  
  Matrice_Creuse<RAB> * sparse_mat =GetAny<Matrice_Creuse<RA>* >((*emat)(stack));
  const Matrix_Prod<RA,RB>  AB = GetAny<Matrix_Prod<RA,RB> >((*prodmat)(stack));
  //  sparse_mat->pUh=AB.A->pUh;
  //sparse_mat->pVh=AB.B->pVh;
  MatriceMorse<RA> *mA= AB.A->A->toMatriceMorse(AB.ta);
  MatriceMorse<RB> *mB= AB.B->A->toMatriceMorse(AB.tb);
  if( !mA && ! mB) ExecError(" Sorry error: in MatProd,  pb trans in MorseMat");
  if( mA->m != mB->n) {
    cerr << "  -- Error dim ProdMat A*B : tA =" << AB.ta << " = tB " << AB.tb << endl;
    cerr << "  --MatProd " << mA->n<< " "<< mA->m << " x " << mB->n<< " "<< mB->m <<  endl;
    ExecError(" Wrong mat dim in MatProd");
  }
  MatriceMorse<RAB> *mAB=new MatriceMorse<RA>();
  mA->prod(*mB,*mAB);
  if(!init) sparse_mat->init();
  sparse_mat->typemat=(mA->n == mB->m) ? TypeSolveMat(TypeSolveMat::GMRES) : TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
  sparse_mat->A.master(mAB);
  delete mA;
  delete mB;
  return sparse_mat;
}

template<class R,int init>
AnyType CombMat(Stack stack,Expression emat,Expression combMat)
{
  using namespace Fem2D;
  
  Matrice_Creuse<R> * sparse_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
  list<triplet<R,MatriceCreuse<R> *,bool> > *  lcB = GetAny<list<triplet<R,MatriceCreuse<R> *,bool> >*>((*combMat)(stack));
  //  sparse_mat->pUh=0;
  // sparse_mat->pVh=0;
   MatriceCreuse<R> * AA=BuildCombMat<R>(*lcB,false,0,0);
   if(!init) sparse_mat->init();
  sparse_mat->A.master(AA);
  sparse_mat->typemat=(AA->n == AA->m) ? TypeSolveMat(TypeSolveMat::GMRES) : TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
  delete lcB;
  return sparse_mat;
}


template<class R>
AnyType MatriceCreuse2map(Stack , const AnyType & mat)
{

  using namespace Fem2D;

  Matrice_Creuse<R> * sparse_mat =GetAny<Matrice_Creuse<R>* >(mat);
  ffassert(sparse_mat);
  int n=sparse_mat->N(),m=sparse_mat->M();
  map<pair<int,int>,R> *M=new map<pair<int,int>,R>;
  if (n >0 && m>0 && sparse_mat->A) 
    {
      sparse_mat->A->addMatTo(R(1.),*M);
      // hack 
      (*M)[make_pair(n-1,m-1)]+=R();
    }
  return M;
}


template<class R,int init>
AnyType DiagMat(Stack stack,Expression emat,Expression edia)
{
  using namespace Fem2D;
  KN<R> * diag=GetAny<KN<R>* >((*edia)(stack));
  Matrice_Creuse<R> * sparse_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
  //  sparse_mat->pUh=0;
  // sparse_mat->pVh=0;
  if(!init) sparse_mat->init();
  sparse_mat->typemat=TypeSolveMat(TypeSolveMat::GC); //  none square matrice (morse)
  sparse_mat->A.master(new MatriceMorse<R>((int) diag->N(),(const R*) *diag));
  return sparse_mat;
}



template<class Rin,class Rout>
 struct  ChangeMatriceMorse {
 static  MatriceMorse<Rout> *f(MatriceMorse<Rin> *mr)
 {
    MatriceMorse<Rout>*  mrr=new MatriceMorse<Rout>(*mr);
    delete mr;
    return mrr;
 }
 };

template<class R> 
 struct  ChangeMatriceMorse<R,R> {
 static MatriceMorse<R>* f(MatriceMorse<R>* mr)
 {
   return mr;
 }
 };
template<class R,class RR,int init>
AnyType CopyMat_tt(Stack stack,Expression emat,Expression eA,bool transp)
{
  using namespace Fem2D;
  Matrice_Creuse<R> * Mat;
 
  if(transp)
   {
    Matrice_Creuse_Transpose<R>  tMat=GetAny<Matrice_Creuse_Transpose<R> >((*eA)(stack));
    Mat=tMat; 
   }
  else   Mat =GetAny<Matrice_Creuse<R>*>((*eA)(stack));
  MatriceMorse<R> * mr=Mat->A->toMatriceMorse(transp,false);
  MatriceMorse<RR> * mrr = ChangeMatriceMorse<R,RR>::f(mr);
  
  Matrice_Creuse<RR> * sparse_mat =GetAny<Matrice_Creuse<RR>* >((*emat)(stack));
  //  sparse_mat->pUh=Mat->pUh;
  // sparse_mat->pVh=Mat->pUh;;
  //  cout << " CopyMat_tt " << init << " "<<transp << " " << sparse_mat << endl;
  if(!init) sparse_mat->init() ;
  sparse_mat->typemat=TypeSolveMat(TypeSolveMat::GC); //  none square matrice (morse)
  sparse_mat->A.master(mrr);
  //delete mr;
  return sparse_mat;
}

template<class R,class RR,int init>
AnyType CopyTrans(Stack stack,Expression emat,Expression eA)
{
 return CopyMat_tt<R,RR,init>(stack,emat,eA,true);
}
template<class R,class RR,int init>
AnyType CopyMat(Stack stack,Expression emat,Expression eA)
{
 return CopyMat_tt<R,RR,init>(stack,emat,eA,false);
}


template<class R,int init>
AnyType MatFull2Sparse(Stack stack,Expression emat,Expression eA)
{
  KNM<R> * A=GetAny<KNM<R>* >((*eA)(stack));
  Matrice_Creuse<R> * sparse_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
  //  sparse_mat->pUh=0;
  // sparse_mat->pVh=0;
  if(!init) sparse_mat->init() ;
  sparse_mat->typemat=TypeSolveMat(TypeSolveMat::GMRES); //  none square matrice (morse)
  sparse_mat->A.master(new MatriceMorse<R>((KNM_<R> &)*A,0.0));
  
 return sparse_mat;
}

template<class R,int init>
AnyType MatMap2Sparse(Stack stack,Expression emat,Expression eA)
{
   map< pair<int,int>, R> * A=GetAny< map< pair<int,int>, R> * >((*eA)(stack));
   int n=0,m=0;
   // hack:  the last element must exist in the map  to set matrix size 

   typename map< pair<int,int>, R>::const_iterator last= --A->end(); // le last element
   
   if( last != A->end() )
   
        { 
                n = last->first.first+1; 
                m=last->first.second+1;
        } 
    
  Matrice_Creuse<R> * sparse_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
     if(!init) sparse_mat->init() ;
  //  sparse_mat->pUh=0;
  // sparse_mat->pVh=0;
  sparse_mat->typemat=TypeSolveMat(TypeSolveMat::GMRES); //  none square matrice (morse)  
  sparse_mat->A.master(new MatriceMorse<R>(n,m,*A,false));
  delete A; 
 return sparse_mat;
}

template<class RR,class AA=RR,class BB=AA> 
struct Op2_pair: public binary_function<AA,BB,RR> { 
  static RR f(const AA & a,const BB & b)  
  { return RR( a, b);} 
};


template<class R>
long get_mat_n(Matrice_Creuse<R> * p)
 { ffassert(p ) ;  return p->A ?p->A->n: 0  ;}

template<class R>
long get_mat_m(Matrice_Creuse<R> * p)
 { ffassert(p ) ;  return p->A ?p->A->m: 0  ;}

template<class R>
long get_mat_nbcoef(Matrice_Creuse<R> * p)
 { ffassert(p ) ;  return p->A ?p->A->NbCoef(): 0  ;}

template<class R>
pair<long,long> get_NM(const list<triplet<R,MatriceCreuse<R> *,bool> > & lM)
{
      typedef typename list<triplet<R,MatriceCreuse<R> *,bool> >::const_iterator lconst_iterator;
    
    lconst_iterator begin=lM.begin();
    lconst_iterator end=lM.end();
    lconst_iterator i;
    
    long n=0,m=0;
    for(i=begin;i!=end;i++++)
     {
       ffassert(i->second);
       MatriceCreuse<R> & M=*i->second;
       bool t=i->third;
       int nn= M.n,mm=M.m;
       if (t) swap(nn,mm);
       if ( n==0) n =  nn;
       if ( m==0) m = mm;
       if (n != 0) ffassert(nn == n);
       if (m != 0) ffassert(mm == m);
      }
   return make_pair(n,m);    
}

template<class R>
long get_diag(Matrice_Creuse<R> * p, KN<R> * x)
 { ffassert(p && x ) ;  return p->A ?p->A->getdiag(*x): 0  ;}
template<class R>
long set_diag(Matrice_Creuse<R> * p, KN<R> * x)
 { ffassert(p && x ) ;  return p->A ?p->A->setdiag(*x): 0  ;}

 
template<class R>  
R * get_elementp2mc(Matrice_Creuse<R> * const  & ac,const long & b,const long & c){ 
   MatriceCreuse<R> * a= ac ? ac->A:0 ;
  if(  !a || a->n <= b || c<0 || a->m <= c  ) 
   { cerr << " Out of bound  0 <=" << b << " < "  << a->n << ",  0 <= " << c << " < "  << a->m
           << " Matrix type = " << typeid(ac).name() << endl;
     cerr << ac << " " << a << endl;
     ExecError("Out of bound in operator Matrice_Creuse<R> (,)");}
   R *  p =a->pij(b,c);
   if( !p) { if(verbosity) cerr << "Error: the coef a(" << b << ","   << c << ")  do'nt exist in sparse matrix "
           << " Matrix  type = " << typeid(ac).name() << endl;
       ExecError("Use of unexisting coef in sparse matrix operator a(i,j) ");}
    return  p;}

template<class R>
R  get_element2mc(Matrice_Creuse<R> * const  & ac,const long & b,const long & c){
    MatriceCreuse<R> * a= ac ? ac->A:0 ;
    R r=R();
    if(  !a || a->n <= b || c<0 || a->m <= c  )
    { cerr << " Out of bound  0 <=" << b << " < "  << a->n << ",  0 <= " << c << " < "  << a->m
        << " Matrix type = " << typeid(ac).name() << endl;
        cerr << ac << " " << a << endl;
        ExecError("Out of bound in operator Matrice_Creuse<R> (,)");}
    R *  p =a->pij(b,c);
    if(p) r=*p;
    return  r;}

template<class RR,class AA=RR,class BB=AA> 
struct Op2_mulAv: public binary_function<AA,BB,RR> { 
  static RR f(const AA & a,const BB & b)  
  { return (*a->A * *b );} 
};

template<class RR,class AA=RR,class BB=AA> 
struct Op2_mulvirtAv: public binary_function<AA,BB,RR> { 
  static RR f(const AA & a,const BB & b)  
  { return RR( (*a).A, b );} 
};

class Matrice_Creuse_C2R  { public:
    typedef Complex K;  
    Matrice_Creuse<K> * A;
    int cas; //  0 re , 1 im 
    Matrice_Creuse_C2R(Matrice_Creuse<K> * AA,int cass) : A(AA),cas(cass) {assert(A);}
    operator MatriceCreuse<K> & () const {return *A->A;}
    operator Matrice_Creuse<K> * () const {return A;}
};

template<int cas> 
Matrice_Creuse_C2R Build_Matrice_Creuse_C2R(Matrice_Creuse<Complex> * pAA)
{
    
  return Matrice_Creuse_C2R(pAA,cas);
}
 
template<class K>
class OneBinaryOperatorA_inv : public OneOperator { public:  
  OneBinaryOperatorA_inv() : OneOperator(atype<Matrice_Creuse_inv<K> >(),atype<Matrice_Creuse<K> *>(),atype<long>()) {}
    E_F0 * code(const basicAC_F0 & args) const 
     { Expression p=args[1];
       if ( ! p->EvaluableWithOutStack() ) 
        { 
          bool bb=p->EvaluableWithOutStack();
          cout << bb << " " <<  * p <<  endl;
          CompileError(" A^p, The p must be a constant == -1, sorry");}
       long pv = GetAny<long>((*p)(NullStack));
        if (pv !=-1)   
         { char buf[100];
           sprintf(buf," A^%ld, The pow must be  == -1, sorry",pv);
           CompileError(buf);}     
       return  new E_F_F0<Matrice_Creuse_inv<K>,Matrice_Creuse<K> *>(Build<Matrice_Creuse_inv<K>,Matrice_Creuse<K> *>,t[0]->CastTo(args[0])); 
    }
};




template<class K>
class Psor :  public E_F0 { public:  
 
   typedef double  Result;
   Expression mat;
   Expression xx,gmn,gmx,oomega;  
   Psor(const basicAC_F0 & args) 
    {   
      args.SetNameParam();
      mat=to<Matrice_Creuse<K> *>(args[0]); 
      gmn=to<KN<K>*>(args[1]); 
      gmx=to<KN<K>*>(args[2]); 
      xx=to<KN<K>*>(args[3]); 
      oomega=to<double>(args[4]); 
      
   }   
    static ArrayOfaType  typeargs() { 
      return  ArrayOfaType( atype<double>(),
                            atype<Matrice_Creuse<K> *>(),
                            atype<KN<K>*>(),
                            atype<KN<K>*>(),
                            atype<KN<K>*>(),
                            atype<double>(),false);}
                            
    static  E_F0 * f(const basicAC_F0 & args){ return new Psor(args);} 
    
    AnyType operator()(Stack s) const {
      Matrice_Creuse<K>* A= GetAny<Matrice_Creuse<K>* >( (*mat)(s) );
      KN<K>* gmin = GetAny<KN<K>* >( (*gmn)(s) );
      KN<K>* gmax = GetAny<KN<K>* >( (*gmx)(s) );
      KN<K>* x = GetAny<KN<K>* >( (*xx)(s) );
      double omega = GetAny<double>((*oomega)(s));
      return A->A->psor(*gmin,*gmax,*x,omega);
    }
  
};
template <class R>
 struct TheDiagMat {
  Matrice_Creuse<R> * A; 
  TheDiagMat(Matrice_Creuse<R> * AA) :A(AA) {ffassert(A);}
  void   get_mat_daig( KN_<R> & x) { ffassert(A && A->A && x.N() == A->A->n  && A->A->n == A->A->m );
     A->A->getdiag(x);}
  void   set_mat_daig(const  KN_<R> & x) { ffassert(A && A->A && x.N() == A->A->n  && A->A->n == A->A->m );
     A->A->setdiag(x);}
 };
 
 template <class R>
 struct TheCoefMat {
  Matrice_Creuse<R> * A; 
  TheCoefMat(Matrice_Creuse<R> * AA) :A(AA) {ffassert(A);}
  void   get_mat_coef( KN_<R> & x) { ffassert(A && A->A && x.N() == A->A->NbCoef()  );
     A->A->getcoef(x);}
  void   set_mat_coef(const  KN_<R> & x) { ffassert(A && A->A && x.N() == A->A->NbCoef() );
     A->A->setcoef(x);}
 };
 
template<class R>
TheDiagMat<R> thediag(Matrice_Creuse<R> * p)
 {  return  TheDiagMat<R>(p);}

template<class R>
TheCoefMat<R> thecoef(Matrice_Creuse<R> * p)
 {  return  TheCoefMat<R>(p);}
 
template<class R>
TheDiagMat<R> set_mat_daig(TheDiagMat<R> dm,KN<R> * x)
{
  dm.set_mat_daig(*x);
  return dm;
}
template<class R>
KN<R> * get_mat_daig(KN<R> * x,TheDiagMat<R> dm)
{
  dm.get_mat_daig(*x);
  return x;
}

template<class R>
TheCoefMat<R> set_mat_coef(TheCoefMat<R> dm,KN<R> * x)
{
  dm.set_mat_coef(*x);
  return dm;
}
template<class R>
KN<R> * get_mat_coef(KN<R> * x,TheCoefMat<R> dm)
{
  dm.get_mat_coef(*x);
  return x;
}

template<class R> 
bool IsRawMat(const basicAC_F0 & args) 
{   
    
    const E_Array * pee= dynamic_cast<const E_Array*>((Expression) args[1]);
    if (!pee) return 0;
    const E_Array &ee=*pee;
    int N=ee.size();
    if (N==1)
    {
	C_F0 c0(ee[0]);
	return 
	    atype<KN_<R> >()->CastingFrom(ee[0].left());
	    
    }
    else if (N==3)
    {
	C_F0 c0(ee[0]),c1(ee[1]),c2(ee[2]);
	return 
	    atype<KN_<long> >()->CastingFrom(ee[0].left())
	    && 	    atype<KN_<long> >()->CastingFrom(ee[1].left())
	    &&      atype<KN_<R> >()->CastingFrom(ee[2].left());
	    
    }
    return 0;
}


template<typename R>
class RawMatrix :  public E_F0 { public:
    int init;
    typedef Matrice_Creuse<R> * Result;
    Expression emat; 
    Expression coef,col,lig;
    RawMatrix(const basicAC_F0 & args,int initt) ;
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<Matrice_Creuse<R>*>(),atype<E_Array>());}
    AnyType operator()(Stack s) const ;    
};

 template<typename R>
 class BlockMatrix :  public E_F0 { public: 
   typedef Matrice_Creuse<R> * Result;
   int N,M;
   int init;
   Expression emat; 
   Expression ** e_Mij;
   int ** t_Mij;
   BlockMatrix(const basicAC_F0 & args,int iinit=0) ;
   ~BlockMatrix() ;
      
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<Matrice_Creuse<R>*>(),atype<E_Array>());}
    static  E_F0 * f(const basicAC_F0 & args){
	if(IsRawMat<R>(args)) return new RawMatrix<R>(args,0);
	else return new BlockMatrix(args,0);
    }  
    AnyType operator()(Stack s) const ;
    
};
template<typename R>
class BlockMatrix1 :  public BlockMatrix<R> { public:
    BlockMatrix1(const basicAC_F0 & args): BlockMatrix<R>(args,1) {}
    static  E_F0 * f(const basicAC_F0 & args){
        if(IsRawMat<R>(args)) return new RawMatrix<R>(args,1);
        else return new BlockMatrix<R>(args,1);
    }  

};

template<typename R>  
map< pair<int,int>, R> *Matrixfull2mapIJ_inv (KNM<R>   * const & pa,const Inv_KN_long & iii,const Inv_KN_long & jjj)
{
   const  KN_<long> &ii(iii), &jj(jjj);
   const KNM<R> & a(*pa);
   int N=a.N(),M=a.M();
   long n = ii(SubArray(N)).max()+1;
   long m= jj(SubArray(M)).max()+1;
   
/*
  long minn = ii(SubArray(N)).min()+1;
  long minm= jj(SubArray(M)).min()+1;
  if ( !(0 <= minn && 0 <=  minm) ) 
  {
  cerr << " Out of Bound  in A(I^-1,J^1) :  "<< minn << " " << minm <<" =>  negative value!! " << endl;
  ExecError("Out of Bound Error");
  }
*/
   
  // cout << "  ### n m " << n << " " << m << endl; 
   map< pair<int,int>, R> *pA= new map< pair<int,int>, R>;
   map< pair<int,int>, R> & A(*pA);
   A[make_pair(n-1,m-1)] = R(); // Hack to be sure that the last term existe 
  
   for (long i=0;i<N;++i)
    for (long j=0;j<M;++j)
     { R aij=a(i,j);
       //cout << i << " " << j << " :: " << ii[i] << " " << jj[j] << " = " << aij << endl;
       if(ii[i]>=0 && jj[j]>=0 && Fem2D::norm(aij)>1e-40) 
         A[make_pair(ii[i],jj[j])] += aij;
     }
      
  return pA;
}

template<typename R>  
map< pair<int,int>, R> *Matrixfull2mapIJ (KNM<R>   * const & pa,const KN_<long> & ii,const  KN_<long> & jj)
{
   const KNM<R> & a(*pa);
   int N=a.N(),M=a.M();
   long n = ii.N();
   long m= jj.N();
  // cout << "  ### n m " << n << " " << m << endl; 
   map< pair<int,int>, R> *pA= new map< pair<int,int>, R>;
   map< pair<int,int>, R> & A(*pA);
   A[make_pair(n-1,m-1)] += R(); // Hack to be sure that the last term existe 
  
   for (long il=0;il<N;++il)
    for (long jl=0;jl<M;++jl)
     { 
       long i = ii[il];
       long j = jj[jl];
       if( i>=0 && j >=0) {
          if ( !(0 <= i && i < N && 0 <= j && j < M ) )
            {
              cerr << " Out of Bound  in A(I,J) : " << i << " " << j << " not in " << "[0,"<<N<<"[x[0," << M << "[ \n";
              ExecError("Out of Bound Error");
             }
       
          R aij=a(i,j);
       //cout << i << " " << j << " :: " << ii[i] << " " << jj[j] << " = " << aij << endl;
	  if (Fem2D::norm(aij)>1e-40) 
           A[make_pair(il,jl)] += aij;
       }
     }
      
  return pA;
}

template<class R>
AnyType Matrixfull2map (Stack , const AnyType & pp)
{
   const KNM<R> & a(*GetAny<KNM<R> *>(pp));
   int N=a.N(),M=a.M();
   int n = N;
   int m= M;
   map< pair<int,int>, R> *pA= new map< pair<int,int>, R>;
   map< pair<int,int>, R> & A(*pA);
   A[make_pair(n-1,m-1)] = R(); // Hack to be sure that the last term existe 
  
   for (int i=0;i<N;++i)
    for (int j=0;j<M;++j)
     { R aij=a(i,j);
       if (Fem2D::norm(aij)>1e-40) 
      A[make_pair(i,j)] += aij;
     }
      
  return pA;
}


template<class R>
map< pair<int,int>, R> *Matrixoutp2mapIJ_inv (outProduct_KN_<R>   * const & pop,const Inv_KN_long & iii,const Inv_KN_long & jjj)
{
   const KN_<long> &ii(iii), &jj(jjj);
   const outProduct_KN_<R> & op(*pop);
   long  N=op.a.N(),M=op.b.N();
   long  n = ii(SubArray(N)).max()+1;
   long m= jj(SubArray(M)).max()+1;
/*
   long minn = ii(SubArray(N)).min()+1;
   long minm= jj(SubArray(M)).min()+1;
     if ( !(0 <= minn && 0 <=  minm) ) 
        {
            cerr << " Out of Bound  in A(I^-1,J^1) :  "<< minn << " " << minm <<" =>  negative value!! " << endl;
            ExecError("Out of Bound Error");
        }
 */
   map< pair<int,int>, R> *pA= new map< pair<int,int>, R>;
   map< pair<int,int>, R> & A(*pA);
   A[make_pair(n-1,m-1)] = R(); // Hack to be sure that the last term existe 
  
   for (int i=0;i<N;++i)
    for (int j=0;j<M;++j)
     { 
       R aij=op.a[i]*RNM::conj(op.b[j]);
       if(ii[i]>=0 && jj[j]>=0 && Fem2D::norm(aij)>1e-40) 
//       if (norm(aij)>1e-40 &) 
          A[make_pair(ii[i],jj[j])] += aij;
     }   
  delete pop;
    
  return pA;
}


template<class R>
map< pair<int,int>, R> *Matrixmapp2mapIJ1 (map< pair<int,int>, R> *const &  B,const Inv_KN_long & iii,const Inv_KN_long  & jjj)
{
    const KN_<long> &ii(iii), &jj(jjj);  
    typedef typename map< pair<int,int>, R>::const_iterator It;
    
    map< pair<int,int>, R> *pA= new map< pair<int,int>, R>;
    map< pair<int,int>, R> & A(*pA);
    int n=0,m=0;
    // hack:  the last element must exist in the map  to set matrix size 
/*    
    It lastb= --B->end(); // le last element
    
    if( lastb != B->end() )	
    { 
	nb = last->first.first+1; 
	mb=last->first.second+1;
    } */
    int N=ii.N(),M=jj.N();
   // A[make_pair(n-1,m-1)] = R(); // Hack to be sure that the last term existe 
    
    for (It k=B->begin();k!=B->end();++k)
    {
	int il =  k->first.first;
	int jl =  k->first.second;
	if ( !( 0 <= il && il < N && 0 <= jl && jl < M )  )
	{
	    cerr << " Out of Bound  in (Map)(I,J) : " << il << " " << jl << " not in " << "[0,"<<N<<"[x[0," << M << "[ \n";
	    ExecError("Out of Bound Error");
	}
	int i=ii(il);
	int j=jj(jl);
	n=max(i,n);
	m=max(j,m);
	R aij =k->second;
	if(i >=0 && j>=0) 
	  A[make_pair(i,j)] += aij;
    }
     A[make_pair(n,m)] += R(); // Hack to be sure that the last term existe 
     delete B;
    
    return pA;
}

template<class R>
map< pair<int,int>, R> *Matrixmapp2mapIJ (map< pair<int,int>, R> *const &  B,const KN_<long> & ii,const KN_<long>  & jj)
{
    
    typedef typename map< pair<int,int>, R>::const_iterator It;
    typedef typename multimap< int,int>::iterator  MI;
    
    map< pair<int,int>, R> *pA= new map< pair<int,int>, R>;
    map< pair<int,int>, R> & A(*pA);
    multimap< int,int > I,J;
    int N=ii.N(),M=jj.N();
    for (int i=0;i<N;++i)
	if(ii[i]>=0)
	  I.insert(make_pair(ii[i],i));
    for (int j=0;j<M;++j)
	if(jj[j]>=0)
	    J.insert(make_pair(jj[j],j));
    int n=N-1,m=M-1;// change FH  sep 2009 to have the correct size..

    for (It k=B->begin();k!=B->end();++k)
    {
	int il =  k->first.first;
	int jl =  k->first.second;
	R aij =k->second;
	pair<MI,MI> PPi=I.equal_range(il);
	pair<MI,MI> PPj=J.equal_range(jl);
	for(MI pi=PPi.first ; pi !=  PPi.second; ++pi)
	{
	    int i=pi->second;
	    for(MI pj=PPj.first ; pj !=  PPj.second; ++pj)
	    { 
		int j=pj->second;
		n=max(i,n);
	        m=max(j,m);	
	       if(i >=0 && j>=0) 
	         A[make_pair(i,j)] += aij;
	    }
	}   
    }
    A[make_pair(n,m)] += R(); // Hack to be sure that the last term existe 
    delete B;
    
    return pA;
}

template<class R>
map< pair<int,int>, R> *Matrixoutp2mapIJ (outProduct_KN_<R>   * const & pop,const KN_<long> & ii,const KN_<long>  & jj)
{
   const outProduct_KN_<R> & op(*pop);
   long N=op.a.N(),M=op.b.N();
   long n=ii.N(),m=jj.N();
   
   map< pair<int,int>, R> *pA= new map< pair<int,int>, R>;
   map< pair<int,int>, R> & A(*pA);
   A[make_pair(n-1,m-1)] = R(); // Hack to be sure that the last term existe 
   
   for (long il=0;il<n;++il)
    for (long jl=0;jl<m;++jl)
     { 
       long i = ii[il];
       long j = jj[jl];
       if(i>=0 && j >=0)
        {
               if ( !( 0 <= i && i < N && 0 <= j && j < M )  )
                {
                    cerr << " Out of Bound  in (a*b')(I,J) : " << i << " " << j << " not in " << "[0,"<<N<<"[x[0," << M << "[ \n";
                    ExecError("Out of Bound Error");
                }
               R aij=op.a[i]*RNM::conj(op.b[j]);
               if (Fem2D::norm(aij)>1e-40) 
                  A[make_pair(il,jl)] += aij;
               }
     }   
  delete pop;
    
  return pA;
}


template<class R>
AnyType Matrixoutp2map (Stack , const AnyType & pp)
{
   const outProduct_KN_<R> & op(*GetAny<outProduct_KN_<R> *>(pp));
   long N=op.a.N(),M=op.b.N();
   long n = N;
   long m= M;
   map< pair<int,int>, R> *pA= new map< pair<int,int>, R>;
   map< pair<int,int>, R> & A(*pA);
   A[make_pair(n-1,m-1)] = R(); // Hack to be sure that the last term existe 
  
   for (long i=0;i<N;++i)
    for (long j=0;j<M;++j)
     { 
      R aij=op.a[i]*RNM::conj(op.b[j]);
      if (Fem2D::norm(aij)>1e-40) 
        A[make_pair(i,j)] += aij;
     } 
  delete &op;        
  return pA;
}


template<typename R>  BlockMatrix<R>::~BlockMatrix() 
{  
    if (e_Mij)
    {  cout << " del Block matrix "<< this << " " << e_Mij <<" N = " << N << " M = " << M << endl;
	for (int i=0;i<N;i++)
	{ delete [] e_Mij[i];
	    delete [] t_Mij[i];
	}
	delete [] e_Mij;
	delete [] t_Mij;
	N=0;
	M=0;
	e_Mij=0;
	t_Mij=0; }
}   

template<typename R>  RawMatrix<R>::RawMatrix(const basicAC_F0 & args,int iinit)
: init(iinit)
{
    args.SetNameParam();
    emat = args[0];
    
    const E_Array & ee= *dynamic_cast<const E_Array*>((Expression) args[1]);
    
    int N=ee.size();
    if (N==1)
    {
	C_F0 c0(ee[0]);
	coef=to<KN_<R> >(ee[0]);
	lig=0;
	col=0;
    }
    
    else if (N==3)
    {
	C_F0 c0(ee[0]),c1(ee[1]),c2(ee[2]);
	coef=to<KN_<R> >(ee[2]);
	lig=to<KN_<long> >(ee[0]);
	col=to<KN_<long> >(ee[1]);
	
    }
    
    
}
template<typename R>  BlockMatrix<R>::BlockMatrix(const basicAC_F0 & args,int iinit)
: init(iinit)
{   
    N=0;
    M=0;
    args.SetNameParam();
    emat = args[0];
    const E_Array & eMij= *dynamic_cast<const E_Array*>((Expression) args[1]);
    N=eMij.size();
    int err =0;
    for (int i=0;i<N;i++)
    {
        const E_Array* emi= dynamic_cast<const E_Array*>((Expression)  eMij[i]);
        if (!emi) err++;
        else
        { 
	    if ( i==0) 
		M = emi->size();
	    else
		if(M != emi->size()) err++;
        }
    }
    if (err) {
	CompileError(" Block matrix : [[ a, b, c], [ a,b,c ]] or Raw Matrix [a] or [ l, c, a ] ");
    }
    assert(N && M);
    e_Mij = new Expression * [N];
    t_Mij = new int * [N];
    for (int i=0;i<N;i++)
    {
	const E_Array li= *dynamic_cast<const E_Array*>((Expression)  eMij[i]);
	
	e_Mij[i] =  new Expression [M];
	t_Mij[i] = new int [M]; 
	for (int j=0; j<M;j++)  
	{
	    C_F0 c_Mij(li[j]);
	    Expression eij=c_Mij.LeftValue();
	    aType rij = c_Mij.left();
	    if ( rij == atype<long>() &&  eij->EvaluableWithOutStack() )
	    {
		long contm = GetAny<long>((*eij)(NullStack));
		/*  prev  version
		if(contm !=0) 
		CompileError(" Block matrix , Just 0 matrix");
		e_Mij[i][j]=0;
		t_Mij[i][j]=0;*/
		if(contm==0)
		{
		    e_Mij[i][j]=0;
		    t_Mij[i][j]=0;
		}
		else if ( atype<R >()->CastingFrom(rij) )
		{  		  // frev 2007 
		    e_Mij[i][j]=to<R>(c_Mij);
		    t_Mij[i][j]=7; //  just un scalaire 
		}
		else CompileError(" Block matrix , Just 0 matrix");
	    }
	    else if ( rij ==  atype<Matrice_Creuse<R> *>()) 
	    {
		e_Mij[i][j]=eij;
		t_Mij[i][j]=1;
	    } 
	    else if ( rij ==  atype<Matrice_Creuse_Transpose<R> >()) 
	    {
		e_Mij[i][j]=eij;
		t_Mij[i][j]=2;
	    } 
	    else if ( atype<KNM<R> *  >()->CastingFrom(rij) )
	      {  //  before KN_ because KNM can be cast in KN_
		  
		  e_Mij[i][j]=to<KNM<R> * >(c_Mij);
		  t_Mij[i][j]=5;
	      }	    
	    else if ( atype<KN_<R> >()->CastingFrom(rij) )
	    {  
		e_Mij[i][j]=to<KN_<R> >(c_Mij);
		t_Mij[i][j]=3;
		
	    } 
	    else if ( atype<Transpose<KN_<R> > >()->CastingFrom(rij) )
	    {  
		
		e_Mij[i][j]=to<Transpose<KN_<R> > >(c_Mij);
		t_Mij[i][j]=4;
	    } 
	    else if ( atype<Transpose< KNM<R> * > >()->CastingFrom(rij) )
	    {  
		
		e_Mij[i][j]=to<Transpose<KNM<R> *> >(c_Mij);
		t_Mij[i][j]=6;
	    }
	    else if ( atype<R >()->CastingFrom(rij) )
	    {  		  // frev 2007 
		e_Mij[i][j]=to<R>(c_Mij);
		t_Mij[i][j]=7; //  just un scalaire 
	    }
	    
	    else {  
		
		CompileError(" Block matrix ,  bad type in block matrix");
	    }
	    /*            else if   ( atype<map< pair<int,int>, R> * >()->CastingFrom(rij) ) 
	    {
		e_Mij[i][j]= to<map< pair<int,int>, R> *>(C_F0(eij,rij)).LeftValue();
		t_Mij[i][j]=10;
	    }*/
	}
	
    }
}

template<typename RR>
class  SetRawMatformMat : public OneOperator { 
public:
    typedef Matrice_Creuse<RR> *  A; // Warning  B type of  2 parameter 
    typedef Matrice_Creuse<RR> *  R;
    typedef E_Array B; //   A type of 1 parameter
    
    class CODE : public  E_F0 { public:
	Expression Mat;
	Expression lig;
	Expression col;
	Expression coef;
	bool mi;    
	    CODE(Expression a,const E_Array & tt)  
		: Mat(a),
		 mi(tt.MeshIndependent())
	    {

		    assert(&tt);
		    if(tt.size()!=3) 
			CompileError("Set raw matrix:  [ lg,col, a] = A (size !=3) ");
		    if (    aatypeknlongp->CastingFrom(tt[0].left() ) //// for  compilation error with g++ 3.2.2 (4 times)
			&&  aatypeknlongp->CastingFrom(tt[1].left() )
			&&  atype<KN<RR>* >()->CastingFrom(tt[2].left() ) )
			    {
			      lig = aatypeknlongp->CastTo(tt[0]);
			      col = aatypeknlongp->CastTo(tt[1]);
			      coef = atype<KN<RR>* >()->CastTo(tt[2]);
			    }      
			    else 
				CompileError(" we are waiting for [ lg,col,a] = A");
    }
	    
	    AnyType operator()(Stack stack)  const 
	    {
		A  a=GetAny<A>((*Mat)(stack));
		KN<long> *lg,*cl;
		KN<RR> *cc;
		lg = GetAny<KN<long>*>((*lig)(stack));
		cl = GetAny<KN<long>*>((*col)(stack));
		cc = GetAny<KN<RR>*>((*coef)(stack));
		int n=a->N(),m=a->M();
		map<pair<int,int>,RR> *M=new map<pair<int,int>,RR>;
		if (n >0 && m>0 && a->A) 
		{
		    a->A->addMatTo(RR(1.),*M);
		    // hack 
		    (*M)[make_pair(n-1,m-1)]+=RR();
		}
		int kk = M->size();
		lg->resize(kk);
		cc->resize(kk);
		cl->resize(kk);
		int k=0;
		typename map<pair<int,int>,RR>::const_iterator i;
		//if (!a->sym)
		 for (i=M->begin(); i != M->end();++i,++k)
		  {  
		    (*lg)[k]= i->first.first;
		    (*cl)[k]= i->first.second;
		    (*cc)[k]= i->second;
		  }
		    
		delete M;
		return SetAny<R>(a);
	    } 
	    bool MeshIndependent() const     {return  mi;} // 
	    ~CODE() {}
	    operator aType () const { return atype<R>();}    
    }; // end sub class CODE
    
    
public: // warning hack  A and B 
	E_F0 * code(const basicAC_F0 & args) const 
    { return  new CODE(t[1]->CastTo(args[1]),*dynamic_cast<const E_Array*>( t[0]->CastTo(args[0]).RightValue()));} 
    SetRawMatformMat():   OneOperator(atype<R>(),atype<B>(),atype<A>())  {} // warning with A and B 
    
};


R realC(Complex c) {return c.real();}
R imagC(Complex c) {return c.imag();}
template<int init>
AnyType CopyMatC2R(Stack stack,Expression emat,Expression CR2eA)
{
    typedef Complex C; 
    typedef double R; 
    using namespace Fem2D;
    Matrice_Creuse_C2R   CRMat =GetAny<Matrice_Creuse_C2R>((*CR2eA)(stack));
    Matrice_Creuse<C> *Mat=CRMat; 
    int cas = CRMat.cas; 
    Matrice_Creuse<R> * sparse_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
    MatriceMorse<C> * mr=Mat->A->toMatriceMorse(false,false);
    MatriceMorse<R> * mrr = 0;
   // cout << " CopyMatC2R:  " << init << " " <<  sparse_mat <<endl;
    if(cas==0) 
        mrr = new MatriceMorse<R>(*mr,realC);
    else if(cas==1) 
        mrr = new MatriceMorse<R>(*mr,imagC);
    else {
         cout << " cas = " << cas <<endl;
        ffassert(0); 
    }
    delete mr;
    if(!init) sparse_mat->init() ; // ???? 
    sparse_mat->A.master(mrr);
    //ffassert(0);// a faire 
    return sparse_mat;
}

template<typename R>  AnyType RawMatrix<R>::operator()(Stack stack) const
{
    MatriceMorse<R> * amorse =0; 
    KN_<R> cc(GetAny< KN_<R>  >((*coef)(stack)));
    int k= cc.N();
    int n= k;
    int m=n;
    map< pair<int,int>, R> Aij;
    bool sym=false;
    if( lig && col)
    {
	KN_<long> lg(GetAny< KN_<long>  >((*lig)(stack)));
	KN_<long> cl=(GetAny< KN_<long>  >((*col)(stack)));
	n = lg.max()+1;
	m = cl.max()+1;
	ffassert( lg.N()==k && cl.N()==k && lg.min()>=0 && lg.max()>=0);
	sym=false;
	for(int i=0;i<k;++i)
	    Aij[make_pair<int,int>(lg[i],cl[i])]+=cc[i];	
    }
    else
    {
	sym=true;
	for(int i=0;i<n;++i)
	    Aij[make_pair(i,i)]=cc[i];
    }
    //cout <<  " nxm  =" <<n<< "x" << m <<endl; 
    amorse=  new MatriceMorse<R>(n,m,Aij,sym); 
    if(verbosity)
	cout << "  -- Raw Matrix    nxm  =" <<n<< "x" << m << " nb  none zero coef. " << amorse->nbcoef << endl;
    
    Matrice_Creuse<R> * sparse_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
    if( !init) sparse_mat->init();
    // sparse_mat->pUh=0;
    // sparse_mat->pVh=0; 
    sparse_mat->A.master(amorse);
    sparse_mat->typemat=(amorse->n == amorse->m) ? TypeSolveMat(TypeSolveMat::GMRES) : TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
    
    if(verbosity>3) { cout << "  End Raw Matrix : " << endl;}
    
    return sparse_mat;   
}
template<typename R>  AnyType BlockMatrix<R>::operator()(Stack s) const
{
  typedef list<triplet<R,MatriceCreuse<R> *,bool> > * L;
   KNM<L> Bij(N,M);
   KNM<KNM_<R> * > Fij(N,M); 
   KNM<bool> cnjij(N,M); 
   KNM<R> Rij(N,M); //  to sto
  
   cnjij = false; 
   KN<long> Oi(N+1), Oj(M+1);
   if(verbosity>3) { cout << " Build Block Matrix : " << N << " x " << M << endl;}
   Bij = (L) 0;
   Oi = (long) 0;
   Oj = (long)0;
  for (int i=0;i<N;++i)
   for (int j=0;j<M;++j)
    {
      Fij(i,j)=0;
      Expression eij = e_Mij[i][j];
      int tij=t_Mij[i][j];
      if (eij) 
      {
        cnjij(i,j) = tij%2 == 0; 
        AnyType e=(*eij)(s);
        if (tij==1) Bij(i,j) = to( GetAny< Matrice_Creuse<R>* >( e)) ;
        else if  (tij==2) Bij(i,j) = to( GetAny<Matrice_Creuse_Transpose<R> >(e));
        else if (tij==3)  { KN_<R> x=GetAny< KN_<R>  >( e);  Fij(i,j) = new KNM_<R>(x,x.N(),1);}
        else if (tij==4)  { KN_<R> x=GetAny< Transpose< KN_<R> >   >( e).t ;  Fij(i,j) = new KNM_<R>(x,1,x.N());}
        else if (tij==5)  { KNM<R> * m= GetAny< KNM<R>*  >( e);  Fij(i,j) = new KNM_<R>(*m);}
        else if (tij==6)  { KNM<R> * m= GetAny< Transpose< KNM<R>* >  >( e).t;  Fij(i,j) = new KNM_<R>(m->t()); }
	else if (tij==7)   { Rij(i,j)=GetAny< R  >( e);  Fij(i,j) = new KNM_<R>(&(Rij(i,j)),1,1);}
        
   //     else if  (tij==3) {}
        else {
         cout << " Bug " << tij << endl;
         ExecError(" Type sub matrix block unknown ");
        }
      }
     }
     //  compute size of matrix
     int err=0;
    for (int i=0;i<N;++i)
     for (int j=0;j<M;++j) 
       {
        pair<long,long> nm(0,0);
        
       if (Bij(i,j)) 
         nm = get_NM( *Bij(i,j));
       else if(Fij(i,j)) 
         nm = make_pair<long,long>(Fij(i,j)->N(), Fij(i,j)->M());
          
        if (( nm.first || nm.second)  && verbosity>3)
          cout << " Block [ " << i << "," << j << " ]      =     " << nm.first << " x " << nm.second << " cnj = " << cnjij(i,j) << endl;
        if (nm.first)
          {
          if ( Oi(i+1) ==0 )  Oi(i+1)=nm.first;
          else  if(Oi(i+1) != nm.first)
            { 
                 err++;
                 cerr <<"Error Block Matrix,  size sub matrix" << i << ","<< j << " n (old) "  << Oi(i+1) 
                       << " n (new) " << nm.first << endl;

            }
          }
          if(nm.second) 
          {  
          if   ( Oj(j+1) ==0) Oj(j+1)=nm.second;
          else   if(Oj(j+1) != nm.second) 
            { 
              cerr <<"Error Block Matrix,  size sub matrix" << i << ","<< j << " m (old) "  << Oj(j+1) 
                   << " m (new) " << nm.second << endl;
              err++;}
          }
        }

    if (err)    ExecError("Error Block Matrix,  size sub matrix");
    //  gestion of zero block ????
    
    for (int j=0;j<M;++j)
    {  if(verbosity>9) cout << j << " colum size" << Oj(j+1) << endl;
        if   ( Oj(j+1) ==0) {
            Oj(j+1)=1;
            if( Oj(j+1) !=1)  err++;}
    }
    for (int i=0;i<N;++i)
    {
        if(verbosity>9) cout << i << " row size" << Oi(i+1) << endl;
        if   ( Oi(i+1) ==0) {
               Oi(i+1)=1;
               if( Oi(i+1) !=1)  err++;}
    }
    if (err)    ExecError("Error Block Matrix with  0 line or  0 colomn..");
   

   
    

//  cout << "Oi = " <<  Oi << endl;
//  cout << "Oj = " <<  Oj << endl;

    for (int i=0;i<N;++i)
      Oi(i+1) += Oi(i);
    for (int j=0;j<M;++j) // correct 10/01/2007 FH 
      Oj(j+1) += Oj(j);// correct 07/03/2010 FH
  long n=Oi(N),m=Oj(M);
  if(verbosity>3)
   {
     cout << "     Oi = " <<  Oi << endl;
     cout << "     Oj = " <<  Oj << endl;
  }
  MatriceMorse<R> * amorse =0; 
{
   map< pair<int,int>, R> Aij;
    for (int i=0;i<N;++i)
     for (int j=0;j<M;++j) 
       if (Bij(i,j)) 
         {
           if(verbosity>3)
             cout << "  Add  Block S " << i << "," << j << " =  at " << Oi(i) << " x " << Oj(j) << " conj = " << cnjij(i,j) << endl;
           BuildCombMat(Aij,*Bij(i,j),false,Oi(i),Oj(j),cnjij(i,j));
         }
       else if (Fij(i,j))
        {
           if(verbosity>3)
             cout << "  Add  Block F " << i << "," << j << " =  at " << Oi(i) << " x " << Oj(j) << endl;
           BuildCombMat(Aij,*Fij(i,j),Oi(i),Oj(j),R(1.),cnjij(i,j));// BuildCombMat
        }
        
           
  amorse=  new MatriceMorse<R>(n,m,Aij,false); 
  }
  if(verbosity)
     cout << "  -- Block Matrix NxM = " << N << "x" << M << "    nxm  =" <<n<< "x" << m << " nb  none zero coef. " << amorse->nbcoef << endl;
 
  Matrice_Creuse<R> * sparse_mat =GetAny<Matrice_Creuse<R>* >((*emat)(s));
  if(!init) sparse_mat->init();
  //sparse_mat->pUh=0;
  // sparse_mat->pVh=0; 
  sparse_mat->A.master(amorse);
  sparse_mat->typemat=(amorse->n == amorse->m) ? TypeSolveMat(TypeSolveMat::GMRES) : TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
                
     
  // cleanning    
  for (int i=0;i<N;++i)
   for (int j=0;j<M;++j)
    if(Bij(i,j)) delete Bij(i,j);
    else if(Fij(i,j))  delete Fij(i,j);  
   if(verbosity>3) { cout << "  End Build Blok Matrix : " << endl;}
   
 return sparse_mat;  

}

template<class R>
class minusMat { public:
    list<triplet<R,MatriceCreuse<R> *,bool> >  *l;
    minusMat(list<triplet<R,MatriceCreuse<R> *,bool> > *ll):
	l(new list<triplet<R,MatriceCreuse<R> *,bool> >(*ll) )
      {
	    typedef typename list<triplet<R,MatriceCreuse<R> *,bool> >::iterator lci;
	    for (lci i= l->begin();i !=l->end();++i)
		i->first*= R(-1);
      }
};

template<class R>
AnyType mM2L3 (Stack , const AnyType & pp)
{
    minusMat<R> mpp(to(GetAny<Matrice_Creuse<R> *>(pp)));
    return SetAny<minusMat<R> >(mpp);
}
/* template<class R>
AnyType mmM2L3 (Stack , const AnyType & pp)
{
    minusMat<R> &  p(GetAny<minusMat<R> >(pp));
    minusMat<R> mpp(p.l);
    delete  p.l;
    return mpp.l;
}

template<class R>
AnyType mmM2L3c (Stack , const AnyType & pp)
{
    list<triplet<R,MatriceCreuse<R> *,bool> >  *  p(GetAny<minusMat<R> >(pp))
    minusMat<R> mpp(p.l);
    delete  p.l;
    return mpp.l;
}
*/
template<class R>
class E_ForAllLoopMatrix
{  public:

    typedef R *VV;
    typedef long KK;
    
    typedef Matrice_Creuse<R> *  Tab;
    
    typedef  ForAllLoopOpBase DataL;
    const DataL *data;
    E_ForAllLoopMatrix(const DataL *t): data(t){}
    AnyType f(Stack s) const {
        Tab t= GetAny<Tab >(data->tab(s));
        KK * i   =   GetAny<KK*>(data->i(s));
        KK * j   =   GetAny<KK*>(data->j(s));
        VV  v   =   GetAny<VV >(data->v(s));
        //      cout << i << " " << j << " " << v << " " << data->epl <<     endl;
        if(verbosity>1000) {
        cout << " i " << (char*) (void *) i -  (char*)(void*) s ;
        cout << " j " << (char*) (void *) j  -  (char*)(void*) s ;
        cout << " vi " <<  (char*) (void *) v -  (char*)(void*) s ;
        cout << endl;
        }
        
        ffassert(i && v);
        MatriceCreuse<R> *m=t->A;
        MatriceMorse<R> *mm = dynamic_cast<MatriceMorse<R>*>(m);
        if(!mm) ExecError(" Matrix sparce of bad type ( not morse ) , sorry.. ");
        if(mm)
            for (long  ii=0;ii < mm->n;++ii)
            for (long  k=mm->lg[ii];k < mm->lg[ii+1];++k)
            {
                *i=ii;
                *j= mm->cl[k];
                *v =  mm->a[k];
                data->code(s);
                mm->a[k] = *v;
            }
        data->end(s);
        return Nothing  ;
    }
    
};

//  Mat real -> mat complex .... ??? FH.   april  2016 ....
struct  VirtualMatCR :public VirtualMatrice<Complex> { public:
    VirtualMatrice<double> & VM;
    typedef Complex R;
    VirtualMatCR(VirtualMatrice<double> & MM): VirtualMatrice<Complex>(MM.N,MM.M),  VM(MM) {}
    void addMatMul(const KN_<R> &  cx, KN_<R> & cy) const {
        double *px = static_cast<double*>(static_cast<void*>(cx));
        double *py = static_cast<double*>(static_cast<void*>(cy));
        KN_<double> rx(px+0,cx.N(),cx.step*2);
        KN_<double> ix(px+1,cx.N(),cx.step*2);
        KN_<double> ry(py+0,cy.N(),cy.step*2);
        KN_<double> iy(py+1,cy.N(),cy.step*2);
        VM.addMatMul(rx,ry);
        VM.addMatMul(ix,iy);
    }
    void addMatTransMul(const KN_<R> &  cx , KN_<R> & cy ) const {
        double *px = static_cast<double*>(static_cast<void*>(cx));
        double *py = static_cast<double*>(static_cast<void*>(cy));
        KN_<double> rx(px+0,cx.N(),cx.step*2);
        KN_<double> ix(px+1,cx.N(),cx.step*2);
        KN_<double> ry(py+0,cy.N(),cy.step*2);
        KN_<double> iy(py+1,cy.N(),cy.step*2);
        VM.addMatTransMul(rx,ry);
        VM.addMatTransMul(ix,iy);
    
    }
        bool WithSolver() const {return VM.WithSolver();} // by default no solver
    virtual void Solve( KN_<R> & cx ,const KN_<R> & cy) const
    { if( !VM.WithSolver()) InternalError("VirtualMatrice::solve not implemented ");
        double *px = static_cast<double*>(static_cast<void*>(cx));
        double *py = static_cast<double*>(static_cast<void*>(cy));
        KN_<double> rx(px+0,cx.N(),cx.step*2);
        KN_<double> ix(px+1,cx.N(),cx.step*2);
        KN_<double> ry(py+0,cy.N(),cy.step*2);
        KN_<double> iy(py+1,cy.N(),cy.step*2);
        VM.Solve(rx,ry);
        VM.Solve(ix,iy);
    }
    
    bool ChecknbLine  (int n) const { return VM.ChecknbLine(n); }
    bool ChecknbColumn  (int m) const  { return VM.ChecknbColumn(m); }
};


template<class R,class A,class B>    // extend (4th arg.)
class  Op2_mulvirtAvCR : public OneOperator {     //
    aType r; //  return type

  
public:
    class CODE :public  E_F0 { public:                               // extend
        Expression a0,a1;          // extend
        CODE( Expression aa0,Expression aa1) : a0(aa0), a1(aa1) {}  // extend (2th arg.)
        AnyType operator()(Stack s)  const
        {
            VirtualMatrice<Complex>  *pv = new  VirtualMatCR ((*GetAny<A>((*a0)(s))).A);
            Add2StackOfPtr2Free(s,pv);
            return SetAny<R>(R(pv,GetAny<B>((*a1)(s))));
        }
        virtual size_t nbitem() const {return a1->nbitem(); } // modif ???
        bool MeshIndependent() const  {return a0->MeshIndependent() && a1->MeshIndependent() ;}
    };

    E_F0 * code(const basicAC_F0 & args) const
    {     if ( args.named_parameter && !args.named_parameter->empty()  )
        CompileError( " They are used Named parameter ");
        
        return  new CODE( t[0]->CastTo(args[0]),
                          t[1]->CastTo(args[1]));}     // extend
    Op2_mulvirtAvCR(int preff=0):                        // 3->4
    OneOperator(map_type[typeid(R).name()],
                map_type[typeid(A).name()],
                map_type[typeid(B).name()]) {pref=preff;}
    
};


template <class R>
void AddSparseMat()
{
// aType tkrp = atype<KN<R> *>(); 
 SetMatrix_Op<R>::btype = Dcl_Type<const  SetMatrix_Op<R> * >();
 Dcl_Type<TheDiagMat<R> >();
 Dcl_Type<TheCoefMat<R> >(); // Add FH oct 2005
 Dcl_Type< map< pair<int,int>, R> * >(); // Add FH mars 2005 
 Dcl_Type<  minusMat<R>  >(); // Add FJH mars 2007
 
 basicForEachType * t_MC=atype<  Matrice_Creuse<R>* >();
    
 t_MC->SetTypeLoop(atype<  R* >(),atype<  long* >(),atype<  long* >());
    
// basicForEachType * t_MCt=atype<  Matrice_Creuse_Transpose<R> >();
// basicForEachType * t_lM=atype< list<triplet<R,MatriceCreuse<R> *,bool> > * >();
// basicForEachType * t_nM=atype<  minusMat<R> >();
 
 basicForEachType * t_MM=atype<map< pair<int,int>, R> * >();
 
TheOperators->Add("*", 
        new OneBinaryOperator<Op2_mulvirtAv<typename VirtualMatrice<R>::plusAx,Matrice_Creuse<R>*,KN_<R> > >,
        new OneBinaryOperator<Op2_mulvirtAv<typename VirtualMatrice<R>::plusAtx,Matrice_Creuse_Transpose<R>,KN_<R> > >,
        new OneBinaryOperator<Op2_mulvirtAv<typename VirtualMatrice<R>::solveAxeqb,Matrice_Creuse_inv<R>,KN_<R> > >     
        );

/*if(0)
TheOperators->Add("*", 
        new OneBinaryOperator<Op2_mulvirtAv<typename VirtualMatrice<R>::plusAx,Matrice_Creuse<R>*,KN_<R> > >( 0  ,tkrp),
        new OneBinaryOperator<Op2_mulvirtAv<typename VirtualMatrice<R>::plusAtx,Matrice_Creuse_Transpose<R>,KN_<R> > >( 0 ,tkrp),
        new OneBinaryOperator<Op2_mulvirtAv<typename VirtualMatrice<R>::solveAxeqb,Matrice_Creuse_inv<R>,KN_<R> > >( 0 ,tkrp)     
        );
*/        
TheOperators->Add("^", new OneBinaryOperatorA_inv<R>());
  
// matrix new code   FH (Houston 2004)        
 TheOperators->Add("=",
//       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation::Op*,E_F_StackF0F0>(SetMatrixInterpolation),
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const Matrix_Prod<R,R>,E_F_StackF0F0>(ProdMat<R,R,R,1>),
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,KN<R> *,E_F_StackF0F0>(DiagMat<R,1>),
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,Matrice_Creuse_Transpose<R>,E_F_StackF0F0>(CopyTrans<R,R,1>), 
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,Matrice_Creuse<R>*,E_F_StackF0F0>(CopyMat<R,R,1>) ,
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,KNM<R>*,E_F_StackF0F0>(MatFull2Sparse<R,1>) ,
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,map< pair<int,int>, R> * ,E_F_StackF0F0>(MatMap2Sparse<R,1>) ,
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,list<triplet<R,MatriceCreuse<R> *,bool> > *,E_F_StackF0F0>(CombMat<R,1>) ,
       new OneOperatorCode<BlockMatrix1<R> >()
       
       );
       
 TheOperators->Add("<-",
       new OneOperatorCode<BlockMatrix<R> >(),
//       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation::Op*,E_F_StackF0F0>(SetMatrixInterpolation),
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const Matrix_Prod<R,R>,E_F_StackF0F0>(ProdMat<R,R,R,0>),
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,KN<R> *,E_F_StackF0F0>(DiagMat<R,0>)  ,
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,Matrice_Creuse_Transpose<R>,E_F_StackF0F0>(CopyTrans<R,R,0>), 
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,Matrice_Creuse<R>*,E_F_StackF0F0>(CopyMat<R,R,0>) ,
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,KNM<R>*,E_F_StackF0F0>(MatFull2Sparse<R,0>) ,
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,map< pair<int,int>, R> * ,E_F_StackF0F0>(MatMap2Sparse<R,0>) ,
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,list<triplet<R,MatriceCreuse<R> *,bool> > *,E_F_StackF0F0>(CombMat<R,0>)
      
       
       );
TheOperators->Add("*", 
        new OneBinaryOperator<Op2_pair<Matrix_Prod<R,R>,Matrice_Creuse<R>*,Matrice_Creuse<R>*> >,
        new OneBinaryOperator<Op2_pair<Matrix_Prod<R,R>,Matrice_Creuse_Transpose<R>,Matrice_Creuse<R>* > >, 
        new OneBinaryOperator<Op2_pair<Matrix_Prod<R,R>,Matrice_Creuse_Transpose<R>,Matrice_Creuse_Transpose<R> > >,
        new OneBinaryOperator<Op2_pair<Matrix_Prod<R,R>,Matrice_Creuse<R>*,Matrice_Creuse_Transpose<R> > > ,
        new OneBinaryOperator<Op2_ListCM<R> >  , 
        new OneBinaryOperator<Op2_ListMC<R> >  ,
	new OneBinaryOperator<Op2_ListCMt<R> >  , 
        new OneBinaryOperator<Op2_ListMtC<R> >  
		  
        );
TheOperators->Add("+", 
        new OneBinaryOperator<Op2_ListCMCMadd<R> >,
        new OneBinaryOperator<Op2_ListCMMadd<R> >,
        new OneBinaryOperator<Op2_ListMCMadd<R> >,
//	new OneBinaryOperator<Op2_ListCMCMadd<R> >(t_MCt,t_lM),
//	new OneBinaryOperator<Op2_ListCMCMadd<R> >(t_MC,t_lM),
        new OneBinaryOperator<Op2_ListMMadd<R> >
       
       ); 
 TheOperators->Add("-",  
	 new OneUnaryOperator<Op1_LCMd<R> >
     );
 Add<Matrice_Creuse<R> *>("n",".",new OneOperator1<long,Matrice_Creuse<R> *>(get_mat_n<R>) );
 Add<Matrice_Creuse<R> *>("m",".",new OneOperator1<long,Matrice_Creuse<R> *>(get_mat_m<R>) );
 Add<Matrice_Creuse<R> *>("nbcoef",".",new OneOperator1<long,Matrice_Creuse<R> *>(get_mat_nbcoef<R>) );
 Add<Matrice_Creuse<R> *>("nnz",".",new OneOperator1<long,Matrice_Creuse<R> *>(get_mat_nbcoef<R>) );
 Add<Matrice_Creuse<R> *>("size",".",new OneOperator1<long,Matrice_Creuse<R> *>(get_mat_nbcoef<R>) );
    
 
 Add<Matrice_Creuse<R> *>("diag",".",new OneOperator1<TheDiagMat<R> ,Matrice_Creuse<R> *>(thediag<R>) );
 Add<Matrice_Creuse<R> *>("coef",".",new OneOperator1<TheCoefMat<R> ,Matrice_Creuse<R> *>(thecoef<R>) );

// Add<Matrice_Creuse<R> *>("setdiag",".",new OneOperator2<long,Matrice_Creuse<R> *,KN<R> *>(set_diag<R>) );
 TheOperators->Add("=", new OneOperator2<KN<R>*,KN<R>*,TheDiagMat<R> >(get_mat_daig<R>) );
 TheOperators->Add("=", new OneOperator2<TheDiagMat<R>,TheDiagMat<R>,KN<R>*>(set_mat_daig<R>) );
 
// TheOperators->Add("=", new OneOperator2<KN<R>*,KN<R>*,TheDiagMat<R> >(get_mat_daig<R>) );
// TheOperators->Add("=", new OneOperator2<TheDiagMat<R>,TheDiagMat<R>,KN<R>*>(set_mat_daig<R>) );
// ADD oct 2005
 TheOperators->Add("=", new OneOperator2<KN<R>*,KN<R>*,TheCoefMat<R> >(get_mat_coef<R>) );
 TheOperators->Add("=", new OneOperator2<TheCoefMat<R>,TheCoefMat<R>,KN<R>*>(set_mat_coef<R>) );
 
// TheOperators->Add("=", new OneOperator2<KN<R>*,KN<R>*,TheCoefMat<R> >(get_mat_coef<R>) );
// TheOperators->Add("=", new OneOperator2<TheCoefMat<R>,TheCoefMat<R>,KN<R>*>(set_mat_coef<R>) );
 
 Global.Add("set","(",new SetMatrix<R>);
 //Global.Add("psor","(",new  OneOperatorCode<Psor<R> > );
 //   a(i,i)
 atype<Matrice_Creuse<R> * >()->Add("(","",new OneOperator3_<R*,Matrice_Creuse<R> *,long,long >(1,get_elementp2mc<R>));
    
 atype<Matrice_Creuse<R> * >()->Add("[","",new OneOperator3_<R,Matrice_Creuse<R> *,long,long >(10,get_element2mc<R>));
    
 atype<KNM<R>*>()->Add("(","",new OneOperator3_<map< pair<int,int>, R> *,KNM<R>*,Inv_KN_long,Inv_KN_long >(Matrixfull2mapIJ_inv<R>));
 atype<KNM<R>*>()->Add("(","",new OneOperator3_<map< pair<int,int>, R> *,KNM<R>*,KN_<long>,KN_<long> >(Matrixfull2mapIJ<R>));
 
 atype<outProduct_KN_<R>*>()->Add("(","",new OneOperator3_<map< pair<int,int>, R> *,outProduct_KN_<R>*,Inv_KN_long,Inv_KN_long >(Matrixoutp2mapIJ_inv<R>));
 atype<outProduct_KN_<R>*>()->Add("(","",new OneOperator3_<map< pair<int,int>, R> *,outProduct_KN_<R>*,KN_<long>,KN_<long> >(Matrixoutp2mapIJ<R>));


 TheOperators->Add("=", new SetRawMatformMat<R>);



 t_MM->Add("(","",  new OneOperator3_<map< pair<int,int>, R> *,map< pair<int,int>, R> *,Inv_KN_long,Inv_KN_long >(Matrixmapp2mapIJ1<R>));
 t_MM->Add("(","",new OneOperator3_<map< pair<int,int>, R> *,map< pair<int,int>, R> *,KN_<long>,KN_<long> >(Matrixmapp2mapIJ<R>));

 t_MC->Add("(","",new OneOperator3_<map< pair<int,int>, R> *,map< pair<int,int>, R> *,Inv_KN_long,Inv_KN_long >(Matrixmapp2mapIJ1<R>,t_MC));
 t_MC->Add("(","",new OneOperator3_<map< pair<int,int>, R> *,map< pair<int,int>, R> *,KN_<long>,KN_<long> >(Matrixmapp2mapIJ<R>,t_MC));

 //atype<outProduct_KN_<R>*>()->Add("(","",new OneOperator3_<map< pair<int,int>, R> *,map< pair<int,int>, R> *,Inv_KN_long,Inv_KN_long >(Matrixmapp2mapIJ1<R>),t_lM);
 //atype<outProduct_KN_<R>*>()->Add("(","",new OneOperator3_<map< pair<int,int>, R> *,map< pair<int,int>, R> *,KN_<long>,KN_<long> >(Matrixmapp2mapIJ<R>),t_lM);

 
//map< pair<int,int>, R> * ttt=   (0);

   //   ; 
 map_type[typeid(map< pair<int,int>, R> *).name()]->AddCast(
     new E_F1_funcT<map< pair<int,int>, R> *,KNM<R>* >(Matrixfull2map<R>),
     new E_F1_funcT<map< pair<int,int>, R> *,outProduct_KN_<R>* >(Matrixoutp2map<R>),
     new E_F1_funcT<map< pair<int,int>, R> *,Matrice_Creuse<R>* >(MatriceCreuse2map<R>)

       ); 

 map_type[typeid(list<triplet<R,MatriceCreuse<R> *,bool> > *).name()]->AddCast(
     new E_F1_funcT<list<triplet<R,MatriceCreuse<R> *,bool> > *,Matrice_Creuse<R>* >(M2L3<R>),
     new E_F1_funcT<list<triplet<R,MatriceCreuse<R> *,bool> > *,Matrice_Creuse_Transpose<R> >(tM2L3<R>),						    
     new E_F1_funcT<list<triplet<R,MatriceCreuse<R> *,bool> > *,minusMat<R> >(mM2L3<R> )
     );
 /*
 map_type[typeid(minusMat<R>).name()]->AddCast( 
	new E_F1_funcT<minusMat<R>,list<triplet<R,MatriceCreuse<R> *,bool> > * >(mmM2L3<R> )
	//new E_F1_funcT<minusMat<R>,Matrice_Creuse<R>*  > * >(mmM2L3c<R> ),
	//new E_F1_funcT<minusMat<R>,Matrice_Creuse_Transpose<R> > * >(mmM2L3ct<R> )
						
      );
*/


    TheOperators->Add("{}",new ForAllLoop<E_ForAllLoopMatrix<R> >);

      
//  --- end  
}


//extern Map_type_of_map map_type_of_map ; //  to store te type 
//extern Map_type_of_map map_pair_of_type ; //  to store te type 
extern int lineno(); 
class  PrintErrorCompile : public OneOperator {
    public: 
    const char * cmm;
    E_F0 * code(const basicAC_F0 & ) const 
     { ErrorCompile(cmm,lineno());
      return 0;} 
    PrintErrorCompile(const char * cc): OneOperator(map_type[typeid(R).name()]),cmm(cc){}

};

class PrintErrorCompileIM :  public E_F0info { public:  
 typedef double  Result;
 static E_F0 *   f(const basicAC_F0 & args)  
    {   
     lgerror("\n\n *** change interplotematrix in interpole.\n  *** Bad name in previous version,\n *** sorry FH.\n\n");
     return 0;  }   
    static ArrayOfaType  typeargs() {return  ArrayOfaType(true);}
    operator aType () const { return atype<double>();} 

};

//  the 2 default sparse solver double and complex
DefSparseSolver<double>::SparseMatSolver SparseMatSolver_R ; ;
DefSparseSolver<Complex>::SparseMatSolver SparseMatSolver_C;
DefSparseSolverSym<double>::SparseMatSolver SparseMatSolverSym_R ; ;
DefSparseSolverSym<Complex>::SparseMatSolver SparseMatSolverSym_C;
// the default probleme solver 
TypeSolveMat::TSolveMat  TypeSolveMatdefaultvalue=TypeSolveMat::defaultvalue;

bool SetDefault()
{
    if(verbosity>1)
	cout << " SetDefault sparse to default" << endl;
    DefSparseSolver<double>::solver =SparseMatSolver_R;
    DefSparseSolver<Complex>::solver =SparseMatSolver_C;
    DefSparseSolverSym<double>::solver =SparseMatSolverSym_R;
    DefSparseSolverSym<Complex>::solver =SparseMatSolverSym_C;
    TypeSolveMat::defaultvalue =TypeSolveMat::SparseSolver;
    return  true;
}


template<class T>
class removeDOF_Op : public E_F0mps {
public:
    Expression A;
    Expression R;
    Expression x;
    Expression out;
    int nbarg;
    static const int n_name_param = 4;
    static basicAC_F0::name_and_type name_param[];
    Expression nargs[n_name_param];
    removeDOF_Op(const basicAC_F0&  args, Expression param1, Expression param2, Expression param3, Expression param4) : A(param1), R(param2), x(param3), out(param4),nbarg(4) {
        args.SetNameParam(n_name_param, name_param, nargs);
    }
    removeDOF_Op(const basicAC_F0&  args, Expression param2, Expression param3, Expression param4) : A(0), R(param2), x(param3), out(param4) ,nbarg(3){
        args.SetNameParam(n_name_param, name_param, nargs);
    }
    
    AnyType operator()(Stack stack) const;
};

template<class T>
basicAC_F0::name_and_type removeDOF_Op<T>::name_param[] = {
    {"symmetrize", &typeid(bool)},
    {"condensation", &typeid(KN<long>*)},
    {"interaction", &typeid(Matrice_Creuse<T>*)},
    {"eps",&typeid(double)}
};

template<class T>
class removeDOF : public OneOperator {
   const  bool withA;
public:
    removeDOF() : OneOperator(atype<long>(), atype<Matrice_Creuse<T>*>(), atype<Matrice_Creuse<double>*>(), atype<KN<T>*>(), atype<KN<T>*>()),withA(true) {}
    removeDOF(int ) : OneOperator(atype<long>(),  atype<Matrice_Creuse<double>*>(), atype<KN<T>*>(), atype<KN<T>*>()),withA(false) {}
   
    E_F0* code(const basicAC_F0& args) const {
        if(withA)
        return new removeDOF_Op<T>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]), t[3]->CastTo(args[3]));
        else
        return new removeDOF_Op<T>(args,  t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]));
    }
};
template<class T> bool cmp(const std::pair<unsigned int, T>& lhs, const std::pair<unsigned int, T>& rhs) { return lhs.first < rhs.first; }
template<class T>
AnyType removeDOF_Op<T>::operator()(Stack stack)  const {
      static const double defEPS=1e-12;
    
    Matrice_Creuse<T>* pA = A ? GetAny<Matrice_Creuse<T>* >((*A)(stack)):0;
    if(pA)
    {
    Matrice_Creuse<T>* pR = GetAny<Matrice_Creuse<T>* >((*R)(stack));
    KN<T>* pX = GetAny<KN<T>* >((*x)(stack));
    KN<T>* pOut = GetAny<KN<T>* >((*out)(stack));
    ffassert( pR && pX && pOut);
    pA->Uh = pR->Uh;
    pA->Vh = pR->Vh;
    MatriceMorse<T> *mA = static_cast<MatriceMorse<T>*>(&(*pA->A));
    MatriceMorse<T> *mR = static_cast<MatriceMorse<T>*>(&(*pR->A));
    bool symmetrize = nargs[0] ? GetAny<bool>((*nargs[0])(stack)) : false;
    double EPS=nargs[3] ? GetAny<double>((*nargs[3])(stack)) :defEPS ;
    KN<long>* condensation = nargs[1] ? GetAny<KN<long>* >((*nargs[1])(stack)) : (KN<long>*) 0;
    
    unsigned int n = condensation ? condensation->n : mR->nbcoef;
    int* lg = new int[n + 1];
    int* cl;
    T* val;
    T* b;
     if(pOut->n != n) pOut->resize(n);
    
    std::vector<signed int> tmpVec;
    if(!condensation) {
        tmpVec.resize(mA->n);
        for(unsigned int i = 0; i < n; ++i)
            tmpVec[mR->cl[i]] = i + 1;
        if(!mA->symetrique) {
            std::vector<std::pair<int, T> > tmp;
            tmp.reserve(mA->nbcoef);
            
            lg[0] = 0;
            for(unsigned int i = 0; i < n; ++i) {
                for(unsigned int j = mA->lg[mR->cl[i]]; j < mA->lg[mR->cl[i] + 1]; ++j) {
                    unsigned int col = tmpVec[mA->cl[j]];
                    if(col != 0 && abs(mA->a[j]) > EPS) {
                        if(symmetrize) {
                            if(col - 1 <= i)
                                tmp.push_back(std::make_pair(col - 1, mA->a[j]));
                        }
                        else
                            tmp.push_back(std::make_pair(col - 1, mA->a[j]));
                    }
                }
                std::sort(tmp.begin() + lg[i], tmp.end(),cmp<T> );
                // c++11 , [](const std::pair<unsigned int, T>& lhs, const std::pair<unsigned int, T>& rhs) { return lhs.first < rhs.first; });
                *(*pOut + i) = *(*pX + mR->cl[i]);
                lg[i + 1] = tmp.size();
            }
            mA->nbcoef = tmp.size();
            if(symmetrize)
                mA->symetrique = true;
            else
                mA->symetrique = false;
            
            cl = new int[tmp.size()];
            val = new T[tmp.size()];
            
            for(unsigned int i = 0; i < tmp.size(); ++i) {
                cl[i]  = tmp[i].first;
                val[i] = tmp[i].second;
            }
        }
        else {
            std::vector<std::vector<std::pair<unsigned int, T> > > tmp(n);
            for(unsigned int i = 0; i < n; ++i)
                tmp[i].reserve(mA->lg[mR->cl[i] + 1] - mA->lg[mR->cl[i]]);
            
            unsigned int nnz = 0;
            for(unsigned int i = 0; i < n; ++i) {
                for(unsigned int j = mA->lg[mR->cl[i]]; j < mA->lg[mR->cl[i] + 1]; ++j) {
                    unsigned int col = tmpVec[mA->cl[j]];
                    if(col != 0 && abs(mA->a[j]) > EPS) {
                        if(i < col - 1)
                            tmp[col - 1].push_back(make_pair(i, mA->a[j]));
                        else
                            tmp[i].push_back(make_pair(col - 1, mA->a[j]));
                        ++nnz;
                    }
                }
                *(*pOut + i) = *(*pX + mR->cl[i]);
            }
            mA->nbcoef = nnz;
            cl = new int[nnz];
            val = new T[nnz];
            nnz = 0;
            lg[0] = 0;
            for(unsigned int i = 0; i < n; ++i) {
                std::sort(tmp[i].begin(), tmp[i].end(),cmp<T>);
                // c++11, [](const std::pair<unsigned int, T>& lhs, const std::pair<unsigned int, T>& rhs) { return lhs.first < rhs.first; });
                for(typename std::vector<std::pair<unsigned int, T> >::const_iterator it = tmp[i].begin(); it != tmp[i].end(); ++it) {
                    cl[nnz] = it->first;
                    val[nnz++] = it->second;
                }
                lg[i + 1] = nnz;
            }
            
        }
        delete [] mA->cl;
        delete [] mA->lg;
        delete [] mA->a;
        mA->n = n;
        mA->m = n;
        mA->N = n;
        mA->M = n;
        mA->lg = lg;
        mA->cl = cl;
        mA->a = val;
    }
    else {
        tmpVec.reserve(mA->n);
        unsigned int i = 0, j = 1;
        for(unsigned int k = 0; k < mA->n; ++k) {
            if(k == *(*condensation + i)) {
                ++i;
                tmpVec.push_back(i);
            }
            else {
                tmpVec.push_back(-j);
                ++j;
            }
        }
        
        
        //        if(!mA->symetrique) {
        std::vector<std::pair<int, T> > tmpInterior;
        std::vector<std::pair<int, T> > tmpBoundary;
        std::vector<std::pair<int, T> > tmpInteraction;
        tmpInterior.reserve(mA->nbcoef);
        tmpBoundary.reserve(mA->nbcoef);
        tmpInteraction.reserve(mA->nbcoef);
        
        lg[0] = 0;
        for(unsigned int i = 0; i < mA->n; ++i) {
            int row = tmpVec[i];
            if(row < 0) {
                for(unsigned int j = mA->lg[i]; j < mA->lg[i + 1]; ++j) {
                    int col = tmpVec[mA->cl[j]];
                    if(col < 0)
                        tmpInterior.push_back(make_pair(-col - 1, mA->a[j]));
                    else
                        tmpInteraction.push_back(make_pair(col - 1, mA->a[j]));
                }
                
            }
            else {
                for(unsigned int j = mA->lg[i]; j < mA->lg[i + 1]; ++j) {
                    int col = tmpVec[mA->cl[j]];
                    if(col > 0)
                        tmpBoundary.push_back(make_pair(col - 1, mA->a[j]));
                }
                // std::sort(tmp.begin() + lg[i], tmp.end());
                *(*pOut + i) = *(*pX + *(*condensation + i));
                lg[i + 1] = tmpBoundary.size();
            }
        }
        cl = new int[tmpBoundary.size()];
        val = new T[tmpBoundary.size()];
        for(unsigned int i = 0; i < tmpBoundary.size(); ++i) {
            cl[i]  = tmpBoundary[i].first;
            val[i] = tmpBoundary[i].second;
        }
        //        }
        MatriceMorse<T>* m = new MatriceMorse<T>(n, n, tmpBoundary.size(), mA->symetrique, val, lg, cl, true);
        pR->typemat = TypeSolveMat(TypeSolveMat::GMRES);
        pR->A.master(m);
        m->dummy = false;
    }
    }
    else
     {
         Matrice_Creuse<T>* pR = GetAny<Matrice_Creuse<T>* >((*R)(stack));
         KN<T>* pX = GetAny<KN<T>* >((*x)(stack));
         KN<T>* pOut = GetAny<KN<T>* >((*out)(stack));
         ffassert(pR && pX && pOut);
         MatriceMorse<T> *mR = static_cast<MatriceMorse<T>*>(&(*pR->A));
        
         unsigned int n = mR->nbcoef;
         
         if(pOut->n != n) pOut->resize(n);
         
         
     
                 for(unsigned int i = 0; i < n; ++i) {
                     *(*pOut + i) = *(*pX + mR->cl[i]);
                 }
         
    }
    return 0L;
}


bool SparseDefault()
{
    return TypeSolveMat::SparseSolver== TypeSolveMat::defaultvalue;
}

bool Have_UMFPACK_=false;
bool Have_UMFPACK() { return Have_UMFPACK_;}
//OneOperator0<bool> *TheSetDefaultSolver=0; // to change the SetDefaultSolver
void  init_lgmat() 

{
  SparseMatSolver_R= DefSparseSolver<double>::solver;
  SparseMatSolver_C= DefSparseSolver<Complex>::solver;
  SparseMatSolverSym_R= DefSparseSolverSym<double>::solver;
  SparseMatSolverSym_C= DefSparseSolverSym<Complex>::solver;

  
  Dcl_Type<const  MatrixInterpolation<pfes>::Op *>(); 
  Dcl_Type<const  MatrixInterpolation<pfes3>::Op *>(); 

  map_type_of_map[make_pair(atype<Matrice_Creuse<double>* >(),atype<double*>())]=atype<Matrice_Creuse<double> *>();
  map_type_of_map[make_pair(atype<Matrice_Creuse<double>* >(),atype<Complex*>())]=atype<Matrice_Creuse<Complex> *>();
  AddSparseMat<double>();
  AddSparseMat<Complex>();
  
  Add<const MatrixInterpolation<pfes>::Op *>("<-","(", new MatrixInterpolation<pfes>);
  Add<const MatrixInterpolation<pfes>::Op *>("<-","(", new MatrixInterpolation<pfes>(1));
  Add<const MatrixInterpolation<pfes>::Op *>("<-","(", new MatrixInterpolation<pfes3>);
  Add<const MatrixInterpolation<pfes>::Op *>("<-","(", new MatrixInterpolation<pfes3>(1,1));
    
    Dcl_Type<const  RestrictArray<pfes>::Op *>();
    Dcl_Type<const  RestrictArray<pfes3>::Op *>();
  //  Add<const RestrictArray<pfes>::Op *>("<-","(", new RestrictArray<pfes>);
  //  Add<const RestrictArray<pfes3>::Op *>("<-","(", new RestrictArray<pfes3>);
  
  Global.Add("restrict","(",new RestrictArray<pfes>);// FH Jan 2014
  Global.Add("restrict","(",new RestrictArray<pfes3>);// FH Jan 2014
    
  TheOperators->Add("=",
                      new OneOperator2_<KN<long>*,KN<long>*,const RestrictArray<pfes>::Op*,E_F_StackF0F0>(SetRestrict<pfes,1>),
                      new OneOperator2_<KN<long>*,KN<long>*,const RestrictArray<pfes3>::Op*,E_F_StackF0F0>(SetRestrict<pfes3,1>)
                      );
    TheOperators->Add("<-",
                      new OneOperator2_<KN<long>*,KN<long>*,const RestrictArray<pfes>::Op*,E_F_StackF0F0>(SetRestrict<pfes,0>),
                      new OneOperator2_<KN<long>*,KN<long>*,const RestrictArray<pfes3>::Op*,E_F_StackF0F0>(SetRestrict<pfes3,0>)
                      );
  

  Global.Add("interpolate","(",new MatrixInterpolation<pfes>);
  Global.Add("interpolate","(",new MatrixInterpolation<pfes>(1));
  Global.Add("interpolate","(",new MatrixInterpolation<pfes3>);
  Global.Add("interpolate","(",new MatrixInterpolation<pfes3>(1,1));
  Global.Add("interplotematrix","(",new  OneOperatorCode<PrintErrorCompileIM>);
  zzzfff->Add("mapmatrix",atype<map< pair<int,int>, double> *>());
  zzzfff->Add("Cmapmatrix",atype<map< pair<int,int>, Complex> *>()); // a voir
  
  Global.Add("defaulttoGMRES","(",new OneOperator0<bool>(SetGMRES));
  Global.Add("defaulttoCG","(",new OneOperator0<bool>(SetCG));
  Global.New("havesparsesolver",CVariable<bool>(SparseDefault));
  Global.Add("defaultsolver","(",new OneOperator0<bool>(SetDefault));
  
  Dcl_Type< Resize<Matrice_Creuse<double> > > ();
  
  Add<Matrice_Creuse<double> *>("resize",".",new OneOperator1< Resize<Matrice_Creuse<double>  >,Matrice_Creuse<double> *>(to_Resize));
  Add<Resize<Matrice_Creuse<double> > >("(","",new OneOperator3_<Matrice_Creuse<double>  *,Resize<Matrice_Creuse<double>  > , long, long  >(resize2)); 
  // add missing in
 Dcl_Type< Resize<Matrice_Creuse<Complex> > > ();
 Add<Matrice_Creuse<Complex> *>("resize",".",new OneOperator1< Resize<Matrice_Creuse<Complex>  >,Matrice_Creuse<Complex> *>(to_Resize));
 Add<Resize<Matrice_Creuse<Complex> > >("(","",new OneOperator3_<Matrice_Creuse<Complex>  *,Resize<Matrice_Creuse<Complex>  > , long, long  >(resize2));   
 
 
 //Global.Add("defaultsolver","(",new OneOperator0<bool>(SetDefaultSolver));
 
 // pour compatibiliter 
 
 TheOperators->Add("=",
		   new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfes>::Op*,E_F_StackF0F0>(SetMatrixInterpolation<1>),
		   new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfes3>::Op*,E_F_StackF0F0>(SetMatrixInterpolation3<1>)
		   );

    
 TheOperators->Add("<-",
		   new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfes>::Op*,E_F_StackF0F0>(SetMatrixInterpolation<0>),
		   new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfes3>::Op*,E_F_StackF0F0>(SetMatrixInterpolation3<0>)
		   );
 // construction of complex matrix form a double matrix
 TheOperators->Add("=", new OneOperator2_<Matrice_Creuse<Complex>*,Matrice_Creuse<Complex>*,Matrice_Creuse<double>*,E_F_StackF0F0>(CopyMat<R,Complex,1>)
		   );
 
 TheOperators->Add("<-", new OneOperator2_<Matrice_Creuse<Complex>*,Matrice_Creuse<Complex>*,Matrice_Creuse<double>*,E_F_StackF0F0>(CopyMat<R,Complex,0>)
		   );
 // Global.Add("imag","(",new OneOperator1_<double,Complex>(Imag));
 //   Add<double>("<--","(",new OneOperator1_<double,Complex>(Real));
    Dcl_Type<Matrice_Creuse_C2R>(); 
    Add<Matrice_Creuse<Complex>*>("re",".",new OneOperator1<Matrice_Creuse_C2R ,Matrice_Creuse<Complex>* >(Build_Matrice_Creuse_C2R<0> ));
    Add<Matrice_Creuse<Complex>*>("im",".",new OneOperator1<Matrice_Creuse_C2R ,Matrice_Creuse<Complex>* >(Build_Matrice_Creuse_C2R<1> ));
    // construction of complex matrix form a double matrix
    TheOperators->Add("=", new OneOperator2_<Matrice_Creuse<Complex>*,Matrice_Creuse<double>*,Matrice_Creuse_C2R,E_F_StackF0F0>(CopyMatC2R<1>));

    TheOperators->Add("<-", new OneOperator2_<Matrice_Creuse<Complex>*,Matrice_Creuse<double>*,Matrice_Creuse_C2R,E_F_StackF0F0>(CopyMatC2R<0>));

 extern  void init_UMFPack_solver();
 init_UMFPack_solver();
 Global.Add("renumbering", "(", new removeDOF<double>);
 Global.Add("renumbering", "(", new removeDOF<Complex>);
    Global.Add("renumbering", "(", new removeDOF<double>(1));
    Global.Add("renumbering", "(", new removeDOF<Complex>(1));

    
    TheOperators->Add("*",
                     new Op2_mulvirtAvCR< VirtualMatrice<Complex>::plusAx,Matrice_Creuse<double>*,KN_<Complex> > ,
                     new Op2_mulvirtAvCR< VirtualMatrice<Complex>::plusAtx,Matrice_Creuse_Transpose<double>,KN_<Complex> > ,
                     new Op2_mulvirtAvCR< VirtualMatrice<Complex>::solveAxeqb,Matrice_Creuse_inv<R>,KN_<Complex> >
                     );
    
}

int Data_Sparse_Solver_version() { return VDATASPARSESOLVER;}
