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
#ifndef REMOVE_CODE
 // to do  F. HECHT
#ifdef __MWERKS__
#pragma optimization_level 0
#endif
#include "ff++.hpp"
#include "array_resize.hpp"

#include "AFunction_ext.hpp"
#include "CGNL.hpp"

namespace bamg { class Triangles; }
namespace Fem2D { void DrawIsoT(const R2 Pt[3],const R ff[3],const RN_ & Viso); }

extern const basicForEachType *aatypeknlongp; //// for  compilation error with g++ 3.2.2
#include "BamgFreeFem.hpp"

#include "PlotStream.hpp"

extern FILE *ThePlotStream;


// Debut FH Houston -------- avril 2004 ---------
//  class for operator on sparse matrix
//   ------------------------------------
//  pour le addition de matrice ----matrice creuse in french)
//  MatriceCreuse    real class for matrix sparse
//  Matrice_Creuse   class for matrix sparse +  poiteur on FE space def
//         to recompute matrice in case of mesh change
//  list<tuple<R,MatriceCreuse<R> *,bool> > * is liste of
//  \sum_i a_i*A_i  where a_i is a scalare and A_i is a sparse matrix
//



template<class R>
list<tuple<R,MatriceCreuse<R> *, bool > > * to(Matrice_Creuse<R> * M)
{
  list<tuple<R,MatriceCreuse<R> *,bool> >  * l=new list<tuple<R,MatriceCreuse<R> *,bool> >;
   l ->push_back(make_tuple<R,MatriceCreuse<R> *,bool>(1,M->A,false));
  return  l;
}
template<class R>
list<tuple<R,MatriceCreuse<R> *, bool > > * to(Matrice_Creuse_Transpose<R>  M)
{
  list<tuple<R,MatriceCreuse<R> *,bool> >  * l=new list<tuple<R,MatriceCreuse<R> *,bool> >;
   l ->push_back(make_tuple<R,MatriceCreuse<R> *,bool>(1,M.A->A,true));
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
struct Op2_ListCM: public binary_function<R,Matrice_Creuse<R> *,list<tuple<R,MatriceCreuse<R> *,bool> > *>
 {
   typedef tuple<R,MatriceCreuse<R> *,bool>  P;
   typedef list<P> L;
   typedef L * RR;
   typedef R  AA;
   typedef Matrice_Creuse<R> * BB;

 static  RR f(const   AA & a,const   BB & b)
  {
    RR r=  new list<P> ;
    P p(a,b->pMC(),false);
    r ->push_back(p);
    return r;}
};
template<class R> void PrintL(const char* cc, list<tuple<R,VirtualMatrix<int,R>*,bool> > const  * const lM)
{
    if(verbosity<100) return; 
    typedef typename list<tuple<R,VirtualMatrix<int,R> *,bool> >::const_iterator lconst_iterator;
    
    lconst_iterator begin=lM->begin();
    lconst_iterator end=lM->end();
    lconst_iterator i;
    cout << cc << " (" ;
    for(i=begin;i!=end;i++++)
    {
        if(std::get<1>(*i)) // M == 0 => zero matrix
        {
            VirtualMatrix<int,R>& M=*std::get<1>(*i);
            bool transpose = std::get<2>(*i) ;
            R coef=std::get<0>(*i);
            cout << "  + " << coef << "*" << &M <<"^" << transpose  ;
            
        }
    }
    
    
    cout << ") "<< endl;
}
template<class R>
struct Op2_ListMC: public binary_function<Matrice_Creuse<R> *,R,list<tuple<R,MatriceCreuse<R> *,bool> > *>
 {
   typedef tuple<R,MatriceCreuse<R> *,bool>  P;
   typedef list<P> L;
   typedef L * RR;
   typedef R  AA;
   typedef Matrice_Creuse<R> * BB;

 static  RR f(const   BB & b,const   AA & a)
  {
    RR r=  new list<P> ;
    P p(a,b->pMC(),false);
    r ->push_back(p);
    return r;}
};
//  ADD FH 16/02/2007

template<class R>
struct Op2_ListCMt: public binary_function<R,Matrice_Creuse_Transpose<R> ,list<tuple<R,MatriceCreuse<R> *,bool> > *>
{
    typedef tuple<R,MatriceCreuse<R> *,bool>  P;
    typedef list<P> L;
    typedef L * RR;
    typedef R  AA;
    typedef Matrice_Creuse_Transpose<R>  BB;

    static  RR f(const   AA & a,const   BB & b)
    {
	RR r=  new list<P> ;
        P p(a,b.A->pMC(),true);
	r ->push_back(p);
	return r;}
};

template<class R>
struct Op2_ListMtC: public binary_function<Matrice_Creuse_Transpose<R> ,R,list<tuple<R,MatriceCreuse<R> *,bool> > *>
{
    typedef tuple<R,MatriceCreuse<R> *,bool>  P;
    typedef list<P> L;
    typedef L * RR;
    typedef R  AA;
    typedef Matrice_Creuse_Transpose<R> BB;

    static  RR f(const   BB & b,const   AA & a)
    {
	RR r=  new list<P> ;
        P p(a,b.A->pMC(),true);
	r ->push_back(p);
	return r;}
};
// FIN ADD 16/02/2007



template<class R,int c=-1>
struct Op1_LCMd: public unary_function<list<tuple<R,MatriceCreuse<R> *,bool> > *,
list<tuple<R,MatriceCreuse<R> *,bool> > *  >
{  //  - ...
    typedef tuple<R,MatriceCreuse<R> *,bool>  P;
    typedef list<P> L;
    typedef L * RR;

    static   RR f(const RR & l)
    {
	typedef typename list<tuple<R,MatriceCreuse<R> *,bool> >::iterator lci;
        for(lci i= l->begin();i !=l->end();++i)
            get<0>(*i) *= R(c);
        PrintL(" - Op1_LCMd: ",l);
	return l;
    }

};

template<class R>
struct Op2_ListCMCMadd: public binary_function<list<tuple<R,MatriceCreuse<R> *,bool> > *,
                                               list<tuple<R,MatriceCreuse<R> *,bool> > *,
                                               list<tuple<R,MatriceCreuse<R> *,bool> > *  >
{  //  ... + ...
   typedef tuple<R,MatriceCreuse<R> *,bool>  P;
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
struct Op2_ListCMCMsub: public binary_function<list<tuple<R,MatriceCreuse<R> *,bool> > *,
list<tuple<R,MatriceCreuse<R> *,bool> > *,
list<tuple<R,MatriceCreuse<R> *,bool> > *  >
{  //  ... + ...
    typedef tuple<R,MatriceCreuse<R> *,bool>  P;
    typedef list<P> L;
    typedef L * RR;

    static   RR f(const RR & a,const RR & b)
    {
        Op1_LCMd<R,-1>::f(b);
        PrintL("Op2_ListCMCMsub +",a);
        PrintL(" -(-) ",b);

        a->insert(a->end(),b->begin(),b->end());
        PrintL(" =>  ",a);
        delete b;
        return a;
    }

};

template<class R,int cc=1>
struct Op2_ListMCMadd: public binary_function<Matrice_Creuse<R> *,
                                              list<tuple<R,MatriceCreuse<R> *,bool> > *,
                                               list<tuple<R,MatriceCreuse<R> *,bool> > *  >
{  //  M + ....
   typedef tuple<R,MatriceCreuse<R> *,bool> P;
   typedef list<P> L;
   typedef L * RR;
   typedef Matrice_Creuse<R> * MM;

  static   RR f(const MM & a,const RR & b)
  {
    // M  + c*L
    Op1_LCMd<R,cc>::f(b);
      PrintL("Op2_ListMCMadd M +",b);

    b->push_front(make_tuple<R,MatriceCreuse<R> *,bool>(R(1.),a->A,false));
    return b;
  }


};

template<class R,int cc=1>
struct Op2_ListCMMadd: public binary_function< list<tuple<R,MatriceCreuse<R> *,bool> > *,
                                               Matrice_Creuse<R> * ,
                                               list<tuple<R,MatriceCreuse<R> *,bool> > *>
{  //   .... + M
   typedef tuple<R,MatriceCreuse<R> *,bool> P;
   typedef list<P> L;
   typedef L * RR;
   typedef Matrice_Creuse<R> * MM;

  static   RR f(const RR & a,const MM & b)
  {
    // L + c*M
    a->push_back(make_tuple<R,MatriceCreuse<R> *,bool>(R(cc),b->A,false));
    return a;
  }


};

template<class R,int cc=1>
struct Op2_ListMMadd: public binary_function< Matrice_Creuse<R> *,
                                              Matrice_Creuse<R> * ,
                                              list<tuple<R,MatriceCreuse<R> *,bool> > *>
{  //  M + M
   typedef tuple<R,MatriceCreuse<R> *,bool> P;
   typedef list<P> L;
   typedef L * RR;
   typedef Matrice_Creuse<R> * MM;

  static   RR f(const MM & a,const MM & b)
  {
    // M + c M
    L * l=to(a);
    l->push_back(make_tuple<R,MatriceCreuse<R> *>(R(cc),b->A,false));
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

    class Op : public E_F0info { public:  // passe de l'info ..
        typedef pfes * A;
        Expression a,b,c,d;
        static const int n_name_param =0;
        Expression nargs[n_name_param];
        bool arg(int i,Stack stack,bool a) const{ return nargs[i] ? GetAny<bool>( (*nargs[i])(stack) ): a;}
        long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
        KN_<long>  arg(int i,Stack stack,KN_<long> a ) const{ return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}

    public:
        Op(const basicAC_F0 &  args,Expression aa,Expression bb,Expression cc) : a(aa),b(bb),c(cc),d(0) {
          args.SetNameParam(n_name_param,0,nargs);

        }
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
template<class pfes1,class pfes2>
class MatrixInterpolation : public OneOperator { public:

    class Op : public E_F0info { public:
       //typedef pfes * A;
       Expression a,b,c,d;
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
   MatrixInterpolation() : OneOperator(atype<const typename MatrixInterpolation<pfes1,pfes2>::Op *>(),
                                       atype<pfes1 *>(),
                                       atype<pfes2 *>()) {}
   // interpolation(Vh,xx,yy) // 2d
   MatrixInterpolation(int bidon) : OneOperator(atype<const typename MatrixInterpolation<pfes1,pfes1>::Op *>(),
                                                atype<pfes1 *>(),atype<KN_<double> >(),atype<KN_<double> >()) {}

    // interpolation(Vh,xx,yy,zz) // 3d
    MatrixInterpolation(int bidon,int bidon2) : OneOperator(atype<const typename MatrixInterpolation<pfes1,pfes1>::Op *>(),
						 atype<pfes1 *>(),atype<KN_<double> >(),atype<KN_<double> >(),atype<KN_<double> >()) {}


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

template<class pfes1,class pfes2>
basicAC_F0::name_and_type  MatrixInterpolation<pfes1,pfes2>::Op::name_param[]= {
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
          ffassert(precon);
          }

       }
       AnyType operator()(Stack stack)  const ;
    };


template<class R>
class SetMatrix : public OneOperator { public:

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

};

template<class R>
AnyType SetMatrix_Op<R>::operator()(Stack stack)  const
{
   Matrice_Creuse<R> *  A= GetAny<Matrice_Creuse<R> *>((*a)(stack));

    ffassert(A);
    Data_Sparse_Solver ds;
    bool VF=false;
    ds.factorize=0;
  // get previous value of sym ??? FH & PHT fev 2020
    int syma=-1;
    if( A->A)
    {  HashMatrix<int,R> * phm= A->pHM();
        if(phm)
        syma=phm->half;
    }

    if( !A->A) A->A.master(new MatriceMorse<R>(0,0,0,0));//  set to empty matrix .. mars 2014 FH ..
    SetEnd_Data_Sparse_Solver<R>(stack,ds,nargs,n_name_param,syma);
    if( verbosity >4) cout <<" SetMatrix_Op " << ds.sym << " "<< ds.positive << " " << syma << endl;
    VirtualMatrix<int,R> *pvm =A->pMC();
    ffassert(pvm);
    pvm->setsdp(ds.sym,ds.positive); // put the matrix in rigth format
    SetSolver<R>(stack,VF,*A->pMC(),ds);

  return Nothing;
}


/*
template<int init>
AnyType SetMatrixInterpolation(Stack,Expression ,Expression);
template<int init>
AnyType SetMatrixInterpolation3(Stack,Expression ,Expression);
template<int init>
AnyType SetMatrixInterpolationS(Stack,Expression ,Expression);
*/
template<class R>
void BuildCombMat(MatriceMorse<R> & mij,const KNM_<R> & A, int ii00=0,int jj00=0,R coef=R(1.),bool cnj=false)
{
  double eps0=numeric_limits<double>::min();
  int i,j;
  int n = A.N(),m=A.M();
  for ( i=0;i<n;i++)
   for ( j=0;j<m;j++)
          {
           R cij=coef*A(i,j);
           if (cnj)  cij = RNM::conj(cij);
             mij[ij_mat(false,ii00,jj00,i,j)] += cij;

   }

}

MatriceMorse<R> * buildInterpolationMatrix(const FESpace & Uh,const FESpace & Vh,void *data)
{  //  Uh = Vh
    MatriceMorse<R> * m=0;
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
       iU2V= idata + 4;
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
  m = new MatriceMorse<R>(n,mm,0,0);
  const  Mesh & ThU =Uh.Th; // line
  const  Mesh & ThV =Vh.Th; // colunm
  bool samemesh =  &Uh.Th == &Vh.Th;  // same Mesh
  int thecolor =0;

  int nnz =0;

  KN<int> color(ThV.nt);
  KN<int> mark(n);
  mark=0;

  color=thecolor++;
  FElement Uh0 = Uh[0];
  FElement Vh0 = Vh[0];

  FElement::aIPJ ipjU(Uh0.Pi_h_ipj());
  FElement::aR2  PtHatU(Uh0.Pi_h_R2());

  int nbdfVK= Vh0.NbDoF();
  int NVh= Vh0.N;

  int nbp= PtHatU.N();
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

   bool whatd[last_operatortype];
   for (int i=0;i<last_operatortype;i++)
     whatd[i]=false;
   whatd[op]=true; // the value of function
   KN<bool> fait(Uh.NbOfDF);
   fait=false;
    R2 Gh(1./3,1./3);
   {

	  for (int it=0;it<ThU.nt;it++)
	    {
	      thecolor++; //  change the current color
	      const Triangle & TU(ThU[it]);
	      FElement KU(Uh[it]);
	      KU.Pi_h(AipjU);
	      if (samemesh)
	        {
	          PV = PtHatU;
	          itV = it;
		  intV= false;// add July 2009 (unset varaible FH)
	        }
	      else
	       {
	          const Triangle *ts=0,*ts0=0;
                   // Search  barycenter
                   bool outside;
                   R2 G;// Add frev 2019  more cleaver F. Hecht
                   ts0=ThV.Find(TU(Gh),G,outside,ts0);// find barycenter good if imbricated mesh ....
                   if(outside) ts0=0; // not find
	          for (int i=0;i<nbp;i++)
	            {
	              ts=ThV.Find(TU(PtHatU[i]),PV[i],outside,ts0);
                        if(!outside && ts0==0) ts0=ts;
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
				      if(Abs(c)>eps)
                                          (*m)(ii,jj) += c;
				  }
			}

	          }

	      for (int df=0;df<KU.NbDoF();df++)
	          {
	           int dfu = KU(df); // le numero de df global
	           fait[dfu]=true;
	           }
	    }
    }
   // sij.clear();
    return m;
}

template<typename RdHatT1,typename RdHatT2>
void copyKNPt( KN<RdHatT1> &A, KN<RdHatT2> B)
{  A=B; };

template<>
void copyKNPt<R2,R1>( KN<R2> &A, KN<R1> B) {
for( int ik=0;ik<A.N();ik++)
    A[ik].x = B[ik].x;
}
template<>
void copyKNPt<R3,R2>( KN<R3> &A, KN<R2> B) {
    for( int ik=0;ik<A.N();ik++){
        A[ik].x = B[ik].x;
        A[ik].y = B[ik].y;
    }
}

template<>
void copyKNPt<R1,R2>( KN<R1> &A, KN<R2> B) {
for( int ik=0;ik<A.N();ik++)
    A[ik].x = B[ik].x;
}
template<>
void copyKNPt<R1,R3>( KN<R1> &A, KN<R3> B) {
    for( int ik=0;ik<A.N();ik++){
        A[ik].x = B[ik].x;
    }
}



template<typename RdHatT1,typename RdHatT2>
void copyPt( RdHatT1 &A, RdHatT2 B) {
    int d1=RdHatT1::d, d2=RdHatT2::d;
    for (int i=0 ; i< min(d1,d2) ;i++)
      A[i]=B[i];
};

template<>
void copyPt<R2,R1>( R2 &A, R1 B) {
    A.x = B.x;
}
template<>
void copyPt<R3,R2>( R3 &A, R2 B) {
    A.x = B.x;
    A.y = B.y;
}
template<>
void copyPt<R3,R1>( R3 &A, R1 B) {
    A.x = B.x;
}

// template 3D RdHat(FESpaceT1) = RdHat(FESpaceT2)
template<class FESpaceT1,class FESpaceT2>
MatriceMorse<R> * buildInterpolationMatrixT(const FESpaceT1 & Uh,const FESpaceT2 & Vh,void *data)
{
  MatriceMorse<R> * m=0;
  typedef typename FESpaceT1::Mesh Mesh1;
  typedef typename FESpaceT1::FElement FElement1;
  typedef typename Mesh1::Element Element1;
  typedef typename FESpaceT1::Rd Rd1;
  typedef typename Element1::RdHat RdHat1;
    
  typedef typename FESpaceT2::Mesh Mesh2;
  typedef typename FESpaceT2::FElement FElement2;
  typedef typename Mesh2::Element Element2;
  typedef typename FESpaceT2::Rd Rd2;
  typedef typename Element2::RdHat RdHat2;
      

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
    iU2V= idata + 4;
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
  m = new MatriceMorse<R>(n,mm,0,0);
    
  RdHat1 Gh= RdHat1::diag(1./(RdHat1::d+1));
  RdHat2 G;
    
  int n1=n+1;
  const  Mesh1 & ThU =Uh.Th; // line
  const  Mesh2 & ThV =Vh.Th; // colunm
  bool samemesh = (is_same< Mesh1, Mesh2 >::value) ? (void*)&Uh.Th == (void*)&Vh.Th : 0 ;  // same Mesh
  int thecolor =0;

  int nnz =0;

  KN<int> color(ThV.nt);
  KN<int> mark(n);
  mark=0;

  int * cl = 0;
  double *a=0;

  color=thecolor++;
  FElement1 Uh0 = Uh[0];
  FElement2 Vh0 = Vh[0];


  int nbdfVK= Vh0.NbDoF();
  int NVh= Vh0.N;

  InterpolationMatrix<RdHat1> ipmat(Uh);

  int nbp=ipmat.np; //
  KN<RdHat2> PV(nbp);   //  the PtHat in ThV mesh
  KN<int> itV(nbp); // the Triangle number
  KN<bool> intV(nbp); // ouside or not

  KNM<R> aaa(nbp,nbdfVK);


  const R eps = 1.0e-10;
  const int sfb1=Vh0.N*last_operatortype*Vh0.NbDoF();
  KN<R> kv(sfb1*nbp);
  R * v = kv;
  KN<int> ik(nbp); // the Triangle number
  op= op==3 ? op_dz : op; //   renumber op ????  dec 2010 FH.
  What_d whatd= 1<< op;
  KN<bool> fait(Uh.NbOfDF);
  fait=false;
  {

    for (int it=0;it<ThU.nt;it++)
    {
      thecolor++; //  change the current color
      const Element1 & TU(ThU[it]);
      FElement1 KU(Uh[it]);
      ipmat.set(KU);
      if (samemesh) {
        copyKNPt<RdHat2,RdHat1>(PV,ipmat.P);
        itV = it;
        intV= false;// add July 2009 (unset varaible FH)
      }
      else
      {
        const Element2 *ts=0,*ts0=0;
        bool outside;
        ts0=ThV.Find(TU(Gh),G,outside,ts0);
        if(outside) ts0=0; // bad starting tet
        for (int i=0;i<nbp;i++)
        {
          ts=ThV.Find(TU(ipmat.P[i]),PV[i],outside,ts0);
          if( ts0 ==0 && !outside) ts0=ts;
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
        int dfu = KU(ipmat.dofe[i]); // le numero de df global
        if(fait[dfu]) continue;
        int jU = ipmat.comp[i]; // la composante dans U
        int p=ipmat.p[i];  //  le point
        if (intV[p]) continue; //  ouside and inside flag => next
        R aipj = ipmat.coef[i];
        FElement2 KV(Vh[itV[p]]);
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
            if(Abs(c)>eps)
              (*m)(ii,jj) += c;
          }
        }

      }

      for (int df=0;df<KU.NbDoF();df++)
      {
        int dfu = KU(df); // le numero de df global
        fait[dfu]=true;
      }


    }

  }
  return m;
}


// ( FESpaceS/FESpaceL - FESpace)
template<class FESpaceT>
MatriceMorse<R> * funcBuildInterpolationMatrixT2(const FESpaceT & Uh,const FESpace & Vh,void *data) {
  MatriceMorse<R> * m=0;
  typedef typename FESpaceT::Mesh MeshT;
  typedef typename FESpaceT::FElement FElementT;
  typedef typename MeshT::Element ElementT;
  typedef typename FESpaceT::Rd RdT;
  typedef typename ElementT::RdHat RdHatT;
    
  typedef typename FESpace::Mesh Mesh2;
  typedef typename FESpace::FElement FElement2;
  typedef typename Mesh2::Element Element2;
  typedef typename FESpace::Rd Rd2;
  typedef typename Element2::RdHat RdHat2;
      

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
    iU2V= idata + 4;
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
  m = new MatriceMorse<R>(n,mm,0,0);
    
  RdHatT Gh= RdHatT::diag(1./(RdHatT::d+1));
  RdHat2 G;
    
  int n1=n+1;
  const  MeshT & ThU =Uh.Th; // line
  const  Mesh2 & ThV =Vh.Th; // colunm
  bool samemesh = (is_same< MeshT, Mesh2 >::value) ? (void*)&Uh.Th == (void*)&Vh.Th : 0 ;  // same Mesh
  int thecolor =0;

  int nnz =0;

  KN<int> color(ThV.nt);
  KN<int> mark(n);
  mark=0;

  int * cl = 0;
  double *a=0;

  color=thecolor++;
  FElementT Uh0 = Uh[0];
  FElement2 Vh0 = Vh[0];
  int nbdfVK= Vh0.NbDoF(), NVh= Vh0.N;

  InterpolationMatrix<RdHatT> ipmat(Uh);

  int nbp=ipmat.np; //
  KN<RdHat2> PV(nbp);   //  the PtHat in ThV mesh
  KN<int> itV(nbp); // the Triangle number
  KN<bool> intV(nbp); // ouside or not

  KNM<R> aaa(nbp,nbdfVK);


  const R eps = 1.0e-6;
  const int sfb1=Vh0.N*last_operatortype*Vh0.NbDoF();
  KN<R> kv(sfb1*nbp);
  R * v = kv;
  KN<int> ik(nbp); // the Triangle number
  
  bool whatd[last_operatortype];
  for (int i=0;i<last_operatortype;i++)
    whatd[i]=false;
  whatd[op]=true; // the value of function
    
  KN<bool> fait(Uh.NbOfDF);
  fait=false;
    double epsP=1e-6; // must be choose

    for (int it=0;it<ThU.nt;it++) {
      thecolor++; //  change the current color
      const ElementT & TU(ThU[it]);
      FElementT KU(Uh[it]);
      ipmat.set(KU);
      if (samemesh) {
        copyKNPt<RdHat2,RdHatT>(PV,ipmat.P);
        itV = it;
        intV= false;// add July 2009 (unset varaible FH)
      }
      else {
        const Element2 *ts=0,*ts0=0;
        bool outside;
        R3 P1(TU(Gh));
        R2 P12(P1.p2());
          if(abs(P1.z)>epsP) {outside=true;ts0=0;}
        else
            ts0=ThV.Find(P12,G,outside,ts0);
        if(outside) ts0=0; // bad starting tet
        for (int i=0;i<nbp;i++) {
            R3 P1(TU(ipmat.P[i]));
            R2 P12(P1.p2());
            if(abs(P1.z)>epsP) {outside=true;ts=0;}
            else ts=ThV.Find(P12,PV[i],outside,ts0);
          if( ts0 ==0 && !outside) ts0=ts;
          if(outside && verbosity>9 )
          cout << it << " " << i << " :: " << TU(ipmat.P[i]) << "  -- "<< outside << PV[i] << " " << ThV(ts) << " ->  " <<  (*ts)(PV[i]) <<endl;
          itV[i]= ThV(ts);
          intV[i]=outside && inside; //  ouside and inside flag
        }
      }

      for (int p=0;p<nbp;p++) {
        KNMK_<R> fb(v+p*sfb1,nbdfVK,NVh,last_operatortype); // valeur de fonction de base de Vh
        // ou:   fb(idf,j,0) valeur de la j composante de la fonction idf
        Vh0.tfe->FB(whatd,ThV,ThV[itV[p]],PV[p],fb);
      }

      for (int i=0;i<ipmat.ncoef;i++) { // pour tous le terme
        int dfu = KU(ipmat.dofe[i]); // le numero de df global
        if(fait[dfu]) continue;
        int jU = ipmat.comp[i]; // la composante dans U
        int p=ipmat.p[i];  //  le point
        if (intV[p]) continue; //  ouside and inside flag => next
        R aipj = ipmat.coef[i];
        FElement2 KV(Vh[itV[p]]);
        int jV=jU;
        if(iU2V) jV=iU2V[jU];

        if(jV>=0 && jV<NVh) {
          KNMK_<R> fb(v+p*sfb1,nbdfVK,NVh,last_operatortype);
          KN_<R> fbj(fb('.',jV,op));

          for (int idfv=0;idfv<nbdfVK;idfv++)
          if (Abs(fbj[idfv])>eps) {
            int dfv=KV(idfv);
            int ii=dfu, jj=dfv;
            if(transpose) Exchange(ii,jj);
            // le term dfu,dfv existe dans la matrice
            R c= fbj[idfv]*aipj;
            if(Abs(c)>eps)
              (*m)(ii,jj) += c;
          }
        }
      }
      for (int df=0;df<KU.NbDoF();df++){
        int dfu = KU(df); // le numero de df global
        fait[dfu]=true;
      }
  }
  return m;
}


// ( FESpace - FESpaceS/FESpaceL)
template<class FESpaceT>
MatriceMorse<R> * funcBuildInterpolationMatrix2T(const FESpace & Uh,const FESpaceT & Vh,void *data) {
  MatriceMorse<R> * m=0;
  typedef typename FESpace::Mesh Mesh1;
  typedef typename FESpace::FElement FElement1;
  typedef typename Mesh1::Element Element1;
  typedef typename FESpace::Rd Rd1;
  typedef typename Element1::RdHat RdHat1;
    
  typedef typename FESpaceT::Mesh MeshT;
  typedef typename FESpaceT::FElement FElementT;
  typedef typename MeshT::Element ElementT;
  typedef typename FESpaceT::Rd RdT;
  typedef typename ElementT::RdHat RdHatT;
      

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
    iU2V= idata + 4;
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
  m = new MatriceMorse<R>(n,mm,0,0);
    
  RdHat1 Gh(1./3,1./3);
  RdHatT G;
    
  int n1=n+1;
  const  Mesh1 & ThU =Uh.Th; // line
  const  MeshT & ThV =Vh.Th; // colunm
  bool samemesh = (is_same< Mesh1, MeshT >::value) ? (void*)&Uh.Th == (void*)&Vh.Th : 0 ;  // same Mesh
  int thecolor =0;

  int nnz =0;

  KN<int> color(ThV.nt);
  KN<int> mark(n);
  mark=0;

  int * cl = 0;
  double *a=0;

  color=thecolor++;
  FElement1 Uh0 = Uh[0];
  FElementT Vh0 = Vh[0];
    
  FElement1::aIPJ ipjU(Uh0.Pi_h_ipj());
  FElement1::aR2  PtHatU(Uh0.Pi_h_R2());
    
  int nbdfVK= Vh0.NbDoF(), NVh= Vh0.N;

  int nbp= PtHatU.N();
    
  KN<RdHatT> PV(nbp);   //  the PtHat in ThV mesh
  KN<int> itV(nbp); // the Triangle number
  KN<bool> intV(nbp); // ouside or not
  KN<R> AipjU(ipjU.N());
  KNM<R> aaa(nbp,nbdfVK);


  const R eps = 1.0e-6;
  const int sfb1=Vh0.N*last_operatortype*Vh0.NbDoF();
  KN<R> kv(sfb1*nbp);
  R * v = kv;
  KN<int> ik(nbp); // the Triangle number
  
  op= op==3 ? op_dz : op; //   renumber op ????  dec 2010 FH.
  What_d whatd= 1<< op;
  
  KN<bool> fait(Uh.NbOfDF);
  fait=false;
  double epsP=1e-6; // must be choose

  for (int it=0;it<ThU.nt;it++) {
    thecolor++; //  change the current color
    const Element1 & TU(ThU[it]);
    FElement1 KU(Uh[it]);
    KU.Pi_h(AipjU);
      if (samemesh) {copyKNPt(PV,PtHatU);
        //PV = PtHatU.x;   // R1 = R2
      itV = it;
      intV= false;// add July 2009 (unset varaible FH)
    }
    else {
      const ElementT *ts=0,*ts0=0;
      bool outside;
      RdHat1 P1(TU(Gh));  // R2
      RdHatT P12;  // R1
          
      copyPt<RdHat1,RdHatT>(P1,P12);
      if(abs(P1.y)>epsP) {outside=true;ts0=0;}
       else
      ts0=ThV.Find(P12,G,outside,ts0);
      if(outside) ts0=0; // bad starting tet
      
      for (int i=0;i<nbp;i++) {
        ts=ThV.Find(TU(PtHatU[i]),PV[i],outside,ts0);
        if(!outside && ts0==0) ts0=ts;
        if(outside && verbosity>9 )
          cout << it << " " << i << " :: " << TU(PtHatU[i]) << "  -- "<< outside << PV[i] << " " << ThV(ts) << " ->  " <<  (*ts)(PV[i]) <<endl;
        itV[i]= ThV(ts);
        intV[i]=outside && inside; //  ouside and inside flag
      }
    }
    for (int p=0;p<nbp;p++) {
      KNMK_<R> fb(v+p*sfb1,nbdfVK,NVh,last_operatortype); // valeur de fonction de base de Vh
      // ou:   fb(idf,j,0) valeur de la j composante de la fonction idf
      Vh0.tfe->FB(whatd,ThV,ThV[itV[p]],PV[p],fb);
    }

    for (int i=0;i<ipjU.N();i++) { // pour tous le terme
      const FElement1::IPJ &ipj_i(ipjU[i]);
      int dfu = KU(ipj_i.i); // le numero de df global
      if(fait[dfu]) continue;
      int jU = ipj_i.j; // la composante dans U
      int p=ipj_i.p;  //  le points
      if (intV[p]) continue; //  ouside and inside flag => next
      R aipj = AipjU[i];
      FElementT KV(Vh[itV[p]]);
      int jV=jU;
      if(iU2V) jV=iU2V[jU];

      if(jV>=0 && jV<NVh) {
        KNMK_<R> fb(v+p*sfb1,nbdfVK,NVh,last_operatortype);
        KN_<R> fbj(fb('.',jV,op));

        for (int idfv=0;idfv<nbdfVK;idfv++)
          if (Abs(fbj[idfv])>eps) {
            int dfv=KV(idfv);
            int ii=dfu, jj=dfv;
            if(transpose) Exchange(ii,jj);
            // le term dfu,dfv existe dans la matrice
            R c= fbj[idfv]*aipj;
            if(Abs(c)>eps)
              (*m)(ii,jj) += c;
          }
      }
    }
    for (int df=0;df<KU.NbDoF();df++){
      int dfu = KU(df); // le numero de df global
      fait[dfu]=true;
    }
  }
  return m;
}

template< >
MatriceMorse<R> * buildInterpolationMatrixT<FESpaceL,FESpace>(const FESpaceL & Uh,const FESpace & Vh,void *data)
{
    return funcBuildInterpolationMatrixT2(Uh,Vh,data);
}

template< >
MatriceMorse<R> * buildInterpolationMatrixT<FESpaceS,FESpace>(const FESpaceS & Uh,const FESpace & Vh,void *data){

     return funcBuildInterpolationMatrixT2(Uh,Vh,data);
}

template< >
MatriceMorse<R> * buildInterpolationMatrixT<FESpace,FESpaceL>(const FESpace & Uh,const FESpaceL & Vh,void *data)
{
    return funcBuildInterpolationMatrix2T(Uh,Vh,data);
}

template< >
MatriceMorse<R> * buildInterpolationMatrixT<FESpace,FESpaceS>(const FESpace & Uh,const FESpaceS & Vh,void *data){

     return funcBuildInterpolationMatrix2T(Uh,Vh,data);
}

template< >
MatriceMorse<R> * buildInterpolationMatrixT<FESpace,FESpace3>(const FESpace & Uh,const FESpace3 & Vh,void *data){

     return funcBuildInterpolationMatrix2T(Uh,Vh,data);
}


MatriceMorse<R> *  buildInterpolationMatrix1(const FESpace & Uh,const KN_<double> & xx,const KN_<double> & yy ,int *data)
{
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
   R2 Phat;
   bool outside;
    MatriceMorse<R> * m = new MatriceMorse<R>(n,mm,0,0);
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
            (*m)(i,j) += c;
      }
      }


  //  sij.clear();
   return m;
}

template<class FESpaceT>
MatriceMorse<R> *  buildInterpolationMatrixT1(const FESpaceT & Uh,const KN_<double> & xx,const KN_<double> & yy ,const KN_<double> & zz,int *data)
{
  typedef typename FESpaceT::Mesh MeshT;
  typedef typename FESpaceT::FElement  FElementT;
  typedef typename MeshT::Element  ElementT;
  typedef typename FESpaceT::Rd  RdT;
  typedef typename ElementT::RdHat  RdHatT;
    
  int op=op_id; //  value of the function
  int icomp=0;
  bool transpose=false;
  bool inside=false;
  if (data){
    transpose=data[0];
    op=data[1];
    inside=data[2];
    icomp = data[3];
    ffassert(op>=0 && op < 4);
    if(op==3) op=op_dz;//  correct missing
  }
  if(verbosity>2){
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
  const  MeshT & ThU =Uh.Th; // line

  FElementT Uh0 = Uh[0];

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
  // map< pair<int,int> , double > sij;
  MatriceMorse<R> * m = new MatriceMorse<R>(n,mm,0,0);
  RdHatT Phat;
  bool outside;

  for(int ii=0;ii<nbxx;ii++){
      if(verbosity>9) cout << " Find ThU " <<ii << ":" <<  RdT(xx[ii],yy[ii],zz[ii]) << endl;
    const ElementT *ts=ThU.Find(RdT(xx[ii],yy[ii],zz[ii]),Phat,outside);
    if(outside && !inside) continue;
    int it = ThU(ts);
    FElementT KU(Uh[it]);
    KNMK_<R> fb(v,nbdfUK,NUh,last_operatortype);
    Uh0.tfe->FB(whatd,ThU,ThU[it],Phat,fb);
    KN_<R> Fwi(fb('.',icomp,op));
    for (int idfu=0;idfu<nbdfUK;idfu++){
      int  j = ii;
      int  i = KU(idfu);
      if(transpose) Exchange(i,j);
      R c = Fwi(idfu);
      if(Abs(c)>eps)
        (*m)(i,j) += c;
    }
  }

  return m;
}


template<>
MatriceMorse<R> *  buildInterpolationMatrixT1<FESpace>(const FESpace & Uh,const KN_<double> & xx,const KN_<double> & yy ,const KN_<double> & zz,int *data)
{
  return buildInterpolationMatrix1(Uh,xx,yy,data);
}




AnyType SetMatrixInterpolation1(Stack stack,Expression emat,Expression einter,int init)
{
  using namespace Fem2D;

  Matrice_Creuse<R> * sparse_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
  const MatrixInterpolation<pfes,pfes>::Op * mi(dynamic_cast<const MatrixInterpolation<pfes,pfes>::Op *>(einter));
  ffassert(einter);
  pfes * pUh = GetAny< pfes * >((* mi->a)(stack));
  FESpace * Uh = **pUh;
  int NUh =Uh->N;
  int* data = new int[4 + NUh];
  data[0]=mi->arg(0,stack,false); // transpose not
  data[1]=mi->arg(1,stack,(long) op_id); ; // get just value
  data[2]=mi->arg(2,stack,false); ; // get just value
  data[3]=mi->arg(3,stack,0L); ; // get just value
  KN<long> U2Vc;
  U2Vc= mi->arg(4,stack,U2Vc); ;
  if( mi->c==0)
  { // old cas
  pfes * pVh = GetAny<  pfes * >((* mi->b)(stack));
  FESpace * Vh = **pVh;
  int NVh =Vh->N;

      for(int i=0;i<NUh;++i)
        data[4+i]=i;//
      for(int i=0;i<min(NUh,(int) U2Vc.size());++i)
	  data[4+i]= U2Vc[i];//
  if(verbosity>3)
	for(int i=0;i<NUh;++i)
	  {
	    cout << "The Uh componante " << i << " -> " << data[4+i] << "  Componante of Vh  " <<endl;
	  }
	  for(int i=0;i<NUh;++i)
	if(data[4+i]>=NVh)
	  {
	      cout << "The Uh componante " << i << " -> " << data[4+i] << " >= " << NVh << " number of Vh Componante " <<endl;
	      ExecError("Interpolation incompability beetween componante ");
	  }

  ffassert(Vh);
  ffassert(Uh);

  if(!init) sparse_mat->init();
      sparse_mat->typemat=0; //TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
  sparse_mat->A.master(buildInterpolationMatrix(*Uh,*Vh,data));
  }
  else
  {  // new cas mars 2006
  KN_<double>  xx = GetAny<  KN_<double>  >((* mi->b)(stack));
  KN_<double>  yy = GetAny<  KN_<double>  >((* mi->c)(stack));
  ffassert( xx.N() == yy.N());
  ffassert(Uh);

  if(!init) sparse_mat->init();
      sparse_mat->typemat=0;//TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
  sparse_mat->A.master(buildInterpolationMatrix1(*Uh,xx,yy,data));
  }
  delete [] data;
   return sparse_mat; // Warning .. no correct gestion of temp ptr ..
}



template<class pfesT1, class FESpaceT1, class pfesT2, class FESpaceT2 >
AnyType SetMatrixInterpolationT1(Stack stack,Expression emat,Expression einter,int init)
{
  using namespace Fem2D;

  Matrice_Creuse<R> * sparse_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
  const typename MatrixInterpolation<pfesT1,pfesT2>::Op * mi(dynamic_cast<const typename MatrixInterpolation<pfesT1,pfesT2>::Op *>(einter));
  ffassert(einter);
  pfesT1 * pUh = GetAny< pfesT1 * >((* mi->a)(stack));
  FESpaceT1 * Uh = **pUh;
  int NUh =Uh->N;
    
  int* data = new int[4 + NUh];
  data[0]=mi->arg(0,stack,false); // transpose not
  data[1]=mi->arg(1,stack,(long) op_id); ; // get just value
  data[2]=mi->arg(2,stack,false); ; // get just value
  data[3]=mi->arg(3,stack,0L); ; // get just value
  KN<long> U2Vc;
  U2Vc= mi->arg(4,stack,U2Vc); ;
  if( mi->c==0)
  { // old cas
    pfesT2 * pVh = GetAny<  pfesT2 * >((* mi->b)(stack));
    FESpaceT2 * Vh = **pVh;
    int NVh =Vh->N;

    for(int i=0;i<NUh;++i)
    data[4+i]=i;//
    for(int i=0;i<min(NUh,(int) U2Vc.size());++i)
    data[4+i]= U2Vc[i];//
    if(verbosity>3)
    for(int i=0;i<NUh;++i)
    {
      cout << "The Uh componante " << i << " -> " << data[4+i] << "  Componante of Vh  " <<endl;
    }
    for(int i=0;i<NUh;++i)
    if(data[4+i]>=NVh)
    {
      cout << "The Uh componante " << i << " -> " << data[4+i] << " >= " << NVh << " number of Vh Componante " <<endl;
      ExecError("Interpolation incompability beetween componante ");
    }

    ffassert(Vh);
    ffassert(Uh);
    if(!init) sparse_mat->init();
    sparse_mat->typemat=0;//(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
    sparse_mat->A.master(buildInterpolationMatrixT<FESpaceT1,FESpaceT2>(*Uh,*Vh,data));	  //  sparse_mat->A.master(new MatriceMorse<R>(*Uh,*Vh,buildInterpolationMatrix,data));
  }
  else
  {  // new cas mars 2006
    KN_<double>  xx = GetAny<  KN_<double>  >((* mi->b)(stack));
    KN_<double>  yy = GetAny<  KN_<double>  >((* mi->c)(stack));
    KN_<double>  zz = GetAny<  KN_<double>  >((* mi->d)(stack));
    ffassert( xx.N() == yy.N());
    ffassert( xx.N() == zz.N());
    ffassert(Uh);
    if(!init) sparse_mat->init();
    sparse_mat->typemat=0;//(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
    sparse_mat->A.master(buildInterpolationMatrixT1<FESpaceT1>(*Uh,xx,yy,zz,data));
  }
  delete [] data;
  return sparse_mat;
}


template<int init>
AnyType SetMatrixInterpolation(Stack stack,Expression emat,Expression einter)
{ return SetMatrixInterpolation1(stack,emat,einter,init);}
template<int init>
AnyType SetMatrixInterpolation3(Stack stack,Expression emat,Expression einter)
{ return SetMatrixInterpolationT1<pfes3,FESpace3,pfes3,FESpace3>(stack,emat,einter,init);}
template<int init>
AnyType SetMatrixInterpolationS(Stack stack,Expression emat,Expression einter)
{ return SetMatrixInterpolationT1<pfesS,FESpaceS,pfesS,FESpaceS>(stack,emat,einter,init);}
template<int init>
AnyType SetMatrixInterpolationL(Stack stack,Expression emat,Expression einter)
{ return SetMatrixInterpolationT1<pfesL,FESpaceL,pfesL,FESpaceL>(stack,emat,einter,init);}
template<int init>
AnyType SetMatrixInterpolationS3(Stack stack,Expression emat,Expression einter)
{ return SetMatrixInterpolationT1<pfesS,FESpaceS,pfes3,FESpace3>(stack,emat,einter,init);}
template<int init>
AnyType SetMatrixInterpolationL2(Stack stack,Expression emat,Expression einter)
{ return SetMatrixInterpolationT1<pfesL,FESpaceL,pfes,FESpace>(stack,emat,einter,init);}
template<int init>
AnyType SetMatrixInterpolationLS(Stack stack,Expression emat,Expression einter)
{ return SetMatrixInterpolationT1<pfesL,FESpaceL,pfesS,FESpaceS>(stack,emat,einter,init);}
template<int init>
AnyType SetMatrixInterpolationS2(Stack stack,Expression emat,Expression einter)
{ return SetMatrixInterpolationT1<pfesS,FESpaceS,pfes,FESpace>(stack,emat,einter,init);}
template<int init>
AnyType SetMatrixInterpolation2L(Stack stack,Expression emat,Expression einter)
{ return SetMatrixInterpolationT1<pfes,FESpace,pfesL,FESpaceL>(stack,emat,einter,init);}
template<int init>
AnyType SetMatrixInterpolation2S(Stack stack,Expression emat,Expression einter)
{ return SetMatrixInterpolationT1<pfes,FESpace,pfesS,FESpaceS>(stack,emat,einter,init);}
template<int init>
AnyType SetMatrixInterpolation23(Stack stack,Expression emat,Expression einter)
{ return SetMatrixInterpolationT1<pfes,FESpace,pfes3,FESpace3>(stack,emat,einter,init);}

template<class RA,class RB,class RAB,int init>
AnyType ProdMat(Stack stack,Expression emat,Expression prodmat)
{
  using namespace Fem2D;

  Matrice_Creuse<RAB> * sparse_mat =GetAny<Matrice_Creuse<RA>* >((*emat)(stack));
  const Matrix_Prod<RA,RB>  AB = GetAny<Matrix_Prod<RA,RB> >((*prodmat)(stack));

  MatriceMorse<RA> *mA= AB.A->pHM();
  MatriceMorse<RB> *mB= AB.B->pHM();
    bool ta = AB.ta, tb = AB.tb;
  if( !mA || !mB)
  {
      ExecError(" numll matrix in pod mat ");
  }
      int An= mA->n, Am =mA->m;
      int Bn= mB->n, Bm =mB->m;
      if(ta) swap(An,Am);
      if(tb) swap(Bn,Bm);

  if(  Am!= Bn) {
    cerr << "  -- Error dim ProdMat A*B : tA =" << AB.ta << " = tB " << AB.tb << endl;
    cerr << "  --MatProd " << mA->n<< " "<< mA->m << " x " << mB->n<< " "<< mB->m <<  endl;
    ExecError(" Wrong mat dim in MatProd");
  }
   MatriceMorse<RAB> *mAB=new MatriceMorse<RAB>(An,Bm,0,0);
   AddMul(*mAB,*mA,*mB,ta,tb);

  if(!init) sparse_mat->init();
    sparse_mat->typemat=0;
  sparse_mat->A.master(mAB);
  return sparse_mat;

}


template<class R,int init>
AnyType CombMat(Stack stack,Expression emat,Expression combMat)
{
  using namespace Fem2D;

  Matrice_Creuse<R> * sparse_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
  list<tuple<R,VirtualMatrix<int,R> *,bool> > *  lcB = GetAny<list<tuple<R,VirtualMatrix<int,R> *,bool> >*>((*combMat)(stack));
  HashMatrix<int,R> * AA=BuildCombMat<R>(*lcB,false,0,0);

   if(!init) sparse_mat->init();
  sparse_mat->A.master(AA);
    sparse_mat->typemat=0;
  delete lcB;
  return sparse_mat;
}
template<class R,int cc> //  July 2019 FH  A += c M +  ...
AnyType AddCombMat(Stack stack,Expression emat,Expression combMat)
{
    using namespace Fem2D;
    
    Matrice_Creuse<R> * pMCA =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
    HashMatrix<int,R> * pA=pMCA->pHM();
    ffassert(pA);
    list<tuple<R,VirtualMatrix<int,R> *,bool> > *  lcB = GetAny<list<tuple<R,VirtualMatrix<int,R> *,bool> >*>((*combMat)(stack));
    if( cc!=1) Op1_LCMd<R,cc>::f(lcB);
    BuildCombMat<R>(*pA,*lcB,false,0,0,false);
   
    delete lcB;
    return pMCA;
}

template<class R,int init>
AnyType DiagMat(Stack stack,Expression emat,Expression edia)
{
  using namespace Fem2D;
  KN<R> * diag=GetAny<KN<R>* >((*edia)(stack));
  Matrice_Creuse<R> * sparse_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
  if(!init) sparse_mat->init();
  sparse_mat->typemat=VirtualMatrix<int,R>::TS_SYM;//TypeSolveMat(TypeSolveMat::GC); //  none square matrice (morse)
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
    MatriceMorse<R> * mr=Mat->pHM();

    Matrice_Creuse<RR> * sparse_mat =GetAny<Matrice_Creuse<RR>* >((*emat)(stack));
    if(mr) {
        MatriceMorse<RR> * mrr = new MatriceMorse<RR>(mr->n,mr->m,0,0);
        *mrr = *mr;
        if(transp) mrr->dotranspose();


        if(!init) sparse_mat->init() ;
        sparse_mat->typemat=Mat->typemat; //  none square matrice (morse)
        sparse_mat->A.master(mrr);
        VirtualMatrix<int,RR> *pvm = sparse_mat->pMC();
        pvm->SetSolver(); // copy solver ???
    }
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
  if(!init) sparse_mat->init() ;
  sparse_mat->typemat=0;//(TypeSolveMat::GMRES); //  none square matrice (morse)
  sparse_mat->A.master(new MatriceMorse<R>((KNM_<R> &)*A,0.0));

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
bool set_mat_COO(Matrice_Creuse<R> * p)
{ ffassert(p ) ;
    HashMatrix<int,R> *phm=p->pHM();
    if(phm) phm->COO();
    return phm;
}
template<class R>
bool set_mat_CSR(Matrice_Creuse<R> * p)
{ ffassert(p ) ;
    HashMatrix<int,R> *phm=p->pHM();
    if(phm) phm->CSR();
    return phm;
}
template<class R>
bool set_mat_CSC(Matrice_Creuse<R> * p)
{ ffassert(p ) ;
    HashMatrix<int,R> *phm=p->pHM();
    if(phm) phm->CSC();
    return phm;
}

template<class R>
long get_mat_m(Matrice_Creuse<R> * p)
 { ffassert(p ) ;  return p->A ?p->A->m: 0  ;}

template<class R>
long get_mat_nbcoef(Matrice_Creuse<R> * p)
 { ffassert(p ) ;  return p->A ?p->A->NbCoef(): 0  ;}
template<class R>
long get_mat_half(Matrice_Creuse<R> * p)
{   return  p && p->pHM() ? p->pHM()->half: -1  ;}

template<class R>
pair<long,long> get_NM(const list<tuple<R,MatriceCreuse<R> *,bool> > & lM)
{
      typedef typename list<tuple<R,MatriceCreuse<R> *,bool> >::const_iterator lconst_iterator;

    lconst_iterator begin=lM.begin();
    lconst_iterator end=lM.end();
    lconst_iterator i;

    long n=0,m=0;
    for(i=begin;i!=end;i++++)
     {
       ffassert(get<1>(*i));
       MatriceCreuse<R> * M=get<1>(*i);
       bool t=get<2>(*i);
       int nn= M->n,mm=M->m;
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
    MatriceMorse<R> * a= ac ? ac->pHM() : 0 ;
  if(  !a || a->n <= b || c<0 || a->m <= c  )
   { cerr << " Out of bound  0 <=" << b << " < "  << (a ? a->n : 0) << ",  0 <= " << c << " < "  << (a ? a->m : 0)
           << " Matrix type = " << typeid(ac).name() << endl;
     cerr << ac << " " << a << endl;
     ExecError("Out of bound in operator Matrice_Creuse<R> (,)");}
   R *  p =a->npij(b,c);
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

// ZZZZZ
R realC(Complex c) {return c.real();}
R imagC(Complex c) {return c.imag();}


template<int cas>
newpMatrice_Creuse<double>  Build_Matrice_Creuse_C2R(Stack stack,Matrice_Creuse<Complex> * const & Mat)
{

    typedef Complex C;
    typedef double R;
    using namespace Fem2D;
    MatriceMorse<C> * mr= Mat->pHM();
    ffassert(mr);
    MatriceMorse<R> * mrr = 0;
    if(cas==0)
        mrr = new MatriceMorse<R>(*mr,realC);
    else if(cas==1)
        mrr = new MatriceMorse<R>(*mr,imagC);
    else {
        cout << " cas = " << cas <<endl;
        ffassert(0);
    }
     return newpMatrice_Creuse<double>(stack,mrr);

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
class OneBinaryOperatorAt_inv : public OneOperator { public:
    OneBinaryOperatorAt_inv() : OneOperator(atype<Matrice_Creuse_inv_trans<K> >(),atype<Matrice_Creuse_Transpose<K> >(),atype<long>()) {}
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
        return  new E_F_F0<Matrice_Creuse_inv_trans<K>,Matrice_Creuse_Transpose<K> >(Build<Matrice_Creuse_inv_trans<K>,Matrice_Creuse_Transpose<K> >,t[0]->CastTo(args[0]));
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
  void  init_get_mat_daig( KN<R> & x) {
      ffassert(A && A->A  && A->A->n == A->A->m );
         x.init(A->A->n);
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
  void   set_mat_coef(const  R & v) { ffassert(A && A->A  );
         KN<R> x(1,0,v);//  vertor constant to v
         A->A->setcoef(x);}

    R trace() { return A->A->trace(); }
 };

template<class R>
R get_trace_mat(Matrice_Creuse<R> * p)
{
    return p ? p->A->trace():0.;

}
template<class R>
bool clear_mat(Matrice_Creuse<R> * p)
{
   if(p)  p->A->clear();
   return true;
}


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
KN<R> * init_get_mat_daig(KN<R> * x,TheDiagMat<R> dm)
{
    dm.init_get_mat_daig(*x);
    return x;
}


template<class R>
TheCoefMat<R> set_mat_coef(TheCoefMat<R> dm,KN<R> * x)
{
  dm.set_mat_coef(*x);
  return dm;
}
template<class R>
TheCoefMat<R> set_mat_coef(TheCoefMat<R> dm,R  x)
{
    dm.set_mat_coef(x);
    return dm;
}
template<class R>
KN<R> * get_mat_coef(KN<R> * x,TheCoefMat<R> dm)
{
  dm.get_mat_coef(*x);
  return x;
}

template<class T> struct  Thresholding {
    Matrice_Creuse<T> *v;
    Thresholding (Matrice_Creuse<T> *vv): v(vv) {}
};

template<class R>
Matrice_Creuse<R>*thresholding2 (const Thresholding<R> &t, const double &threshold) {
    typedef HashMatrix<int,R> HMat;
    Matrice_Creuse<R> *sparse_mat = t.v;
    if (sparse_mat) {
        HMat *phm=sparse_mat->pHM() ;
        if( phm)
        {
            int n = phm->n, m = phm->m;
            int nnzo = phm->nnz;
            phm->resize(n,m,0,threshold);
            if (verbosity) {cout << "  thresholding : remove " << nnzo-phm->nnz  << " them in the matrix " << sparse_mat << " " << threshold << endl;}
        } else if (verbosity) {cout << " empty matrix " << sparse_mat << endl;}
    }
    
    return t.v;
}

template<class T>
long symmetrizeCSR (Matrice_Creuse<T> *const &sparse_mat)
{
    
    typedef HashMatrix<int,T> HMat;
    if (sparse_mat) {
        HMat *phm=sparse_mat->pHM() ;
        if( phm)
        {
            int n = phm->n, m = phm->m;
            int nnzo = phm->nnz;
            phm->resize(n,m,0,-1,true);
            if (verbosity) {cout << "  symmetrizeCSR remove " << (long) nnzo-(long) phm->nnz   << " them in the matrix " << sparse_mat << endl;}
        } else if (verbosity) {cout << " empty matrix " << sparse_mat << endl;}
    }
    
    return 1L;
}

template<class T>
Thresholding<T> to_Thresholding (Matrice_Creuse<T> *v) {return Thresholding<T>(v);}

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
newpMatrice_Creuse<R> Matrixfull2mapIJ_inv (Stack s,KNM<R>   * const & pa,const Inv_KN_long & iii,const Inv_KN_long & jjj)
{

   const  KN_<long> &ii(iii), &jj(jjj);
   const KNM<R> & a(*pa);
   int N=a.N(),M=a.M();
   long n = ii(SubArray(N)).max()+1;
   long m= jj(SubArray(M)).max()+1;


    HashMatrix<int,R> *pA= new  HashMatrix<int,R>((int) n,(int) m,0,0);
    HashMatrix<int,R> & A =*pA;

   for (long i=0;i<N;++i)
    for (long j=0;j<M;++j)
     { R aij=a(i,j);
       if(ii[i]>=0 && jj[j]>=0 && std::norm(aij)>1e-40)
         A(ii[i],jj[j]) += aij;
     }

  return newpMatrice_Creuse<R>(s,pA);
}

template<typename R>
newpMatrice_Creuse<R>  Matrixfull2mapIJ (Stack s, KNM<R>   * const & pa,const KN_<long> & ii,const  KN_<long> & jj)
{
   const KNM<R> & a(*pa);
   int N=a.N(),M=a.M();
   long n = ii.N();
   long m= jj.N();
    HashMatrix<int,R> *pA= new  HashMatrix<int,R>((int)n,(int)m,0,0);
    HashMatrix<int,R> & A =*pA;

   for (long il=0;il<n;++il)// correct juil 2017 FH N --> n
    for (long jl=0;jl<m;++jl)// correct juil 2017 FH M --> m
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
	  if (std::norm(aij)>1e-40)
           A(il,jl) += aij;
       }
     }

    return newpMatrice_Creuse<R> (s,pA);//;pA;
}

template<class R>
AnyType Matrixfull2map (Stack s , const AnyType & pp)
{
   const KNM<R> & a(*GetAny<KNM<R> *>(pp));
   int N=a.N(),M=a.M();
   int n = N;
   int m= M;
    HashMatrix<int,R> *pA= new  HashMatrix<int,R>(n,m,0,0);
    HashMatrix<int,R> & A =*pA;

   A(n-1,m-1) = R(); // Hack to be sure that the last term existe

   for (int i=0;i<N;++i)
    for (int j=0;j<M;++j)
     { R aij=a(i,j);
       if (std::norm(aij)>1e-40)
      A(i,j) += aij;
     }

  return SetAny<newpMatrice_Creuse<R> >(newpMatrice_Creuse<R> (s,pA));//
}


template<class R>
newpMatrice_Creuse<R>  Matrixoutp2mapIJ_inv (Stack s,outProduct_KN_<R>   * const & pop,const Inv_KN_long & iii,const Inv_KN_long & jjj)
{
   const KN_<long> &ii(iii), &jj(jjj);
   const outProduct_KN_<R> & op(*pop);
   long  N=op.a.N(),M=op.b.N();
   long  n = ii(SubArray(N)).max()+1;
   long m= jj(SubArray(M)).max()+1;
    HashMatrix<int,R> *pA= new  HashMatrix<int,R>((int)n,(int)m,0,0);
    HashMatrix<int,R> & A =*pA;

   for (int i=0;i<N;++i)
    for (int j=0;j<M;++j)
     {
       R aij=op.a[i]*RNM::conj(op.b[j])*op.c;
       if(ii[i]>=0 && jj[j]>=0 && std::norm(aij)>1e-40)
          A(ii[i],jj[j]) += aij;
     }
  delete pop;

 return newpMatrice_Creuse<R> (s,pA);
}


template<class R>
newpMatrice_Creuse<R>
Matrixmapp2mapIJ1 (Stack s,Matrice_Creuse<R> *const &  mcB,const Inv_KN_long & iii,const Inv_KN_long  & jjj)
{
    const KN_<long> &ii(iii), &jj(jjj);
    typedef typename map< pair<int,int>, R>::const_iterator It;

      int n=0,m=0;
     HashMatrix<int,R> *B = dynamic_cast<HashMatrix<int,R> *>(mcB->pMC());
     ffassert(B);

    int N=ii.N(),M=jj.N();
    int nn = ii.max(), mm= jj.max();
    HashMatrix<int,R> *pA= new  HashMatrix<int,R>(nn,mm,0,0);
    HashMatrix<int,R> & A =*pA;

    for (int k=0;k<B->nnz;++k)
    {
        int il =  B->i[k];
	int jl =  B->j[k];
	if ( !( 0 <= il && il < N && 0 <= jl && jl < M )  )
	{
	    cerr << " Out of Bound  in (Map)(I,J) : " << il << " " << jl << " not in " << "[0,"<<N<<"[x[0," << M << "[ \n";
	    ExecError("Out of Bound Error");
	}
	int i=ii(il);
	int j=jj(jl);
	n=max(i,n);
	m=max(j,m);
	R aij = B->aij[k];
	if(i >=0 && j>=0)
	  A(i,j) += aij;
    }
    // A[make_pair(n,m)] += R(); // Hack to be sure that the last term existe
    // delete B;
    //  FH question resize or not pA  to (n,m);
 return  newpMatrice_Creuse<R> (s,pA);//
}

template<class R>
// map< pair<int,int>, R> *
newpMatrice_Creuse<R> Matrixmapp2mapIJ (Stack s,Matrice_Creuse<R> *const &  mcB,const KN_<long> & ii,const KN_<long>  & jj)
{


 typedef typename map< pair<int,int>, R>::const_iterator It;
typedef typename multimap< int,int>::iterator  MI;
    HashMatrix<int,R> *B = dynamic_cast<HashMatrix<int,R> *>(mcB->pMC());
    ffassert(B);
    multimap< int,int > I,J;
    int N=ii.N(),M=jj.N();
    for (int i=0;i<N;++i)
	if(ii[i]>=0)
	  I.insert(make_pair(ii[i],i));
    for (int j=0;j<M;++j)
	if(jj[j]>=0)
	    J.insert(make_pair(jj[j],j));
    int n=N-1,m=M-1;// change FH  sep 2009 to have the correct size..
    HashMatrix<int,R> *pA= new  HashMatrix<int,R>(N,M,0,0);
    HashMatrix<int,R> & A =*pA;

    for (int k=0;k!=B->nnz;++k)
    {
        int il =  B->i[k];
	int jl =   B->j[k];
	R aij = B->aij[k];
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
	         A(i,j) += aij;
	    }
	}
    }
    // resier A to n,m ?????
   // A[make_pair(n,m)] += R(); // Hack to be sure that the last term existe
   // delete B;

   // return pA;
    return  newpMatrice_Creuse<R> (s,pA);//
}

template<class R>
newpMatrice_Creuse<R>
Matrixoutp2mapIJ (Stack s,outProduct_KN_<R>   * const & pop,const KN_<long> & ii,const KN_<long>  & jj)
{
   const outProduct_KN_<R> & op(*pop);
   long N=op.a.N(),M=op.b.N();
   long n=ii.N(),m=jj.N();
    R c = op.c;
    HashMatrix<int,R> *pA= new  HashMatrix<int,R>((int)n,(int)m,0,0);
    HashMatrix<int,R> & A =*pA;

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
               R aij=op.a[i]*RNM::conj(op.b[j])*c;
               if (std::norm(aij)>1e-40)
                  A(il,jl) += aij;
               }
     }
  delete pop;
  return  newpMatrice_Creuse<R> (s,pA);
}


template<class R>
AnyType Matrixoutp2map (Stack s, const AnyType & pp)
{
   const outProduct_KN_<R> & op(*GetAny<outProduct_KN_<R> *>(pp));
   long N=op.a.N(),M=op.b.N();
   long n = N;
   long m= M;
//   A[make_pair(n-1,m-1)] = R(); // Hack to be sure that the last term existe
    HashMatrix<int,R> *pA= new  HashMatrix<int,R>((int)n,(int)m,0,0);
    HashMatrix<int,R> & A =*pA;

   for (long i=0;i<N;++i)
    for (long j=0;j<M;++j)
     {
      R aij=op.a[i]*RNM::conj(op.b[j])*op.c;
      if (std::norm(aij)>1e-40)
        A(i,j) += aij;
     }
  delete &op;
 return  SetAny<newpMatrice_Creuse<R>>(newpMatrice_Creuse<R> (s,pA));//
}


template<typename R>  BlockMatrix<R>::~BlockMatrix()
{
    if (e_Mij)
    {   if(verbosity>9999) cout << " del Block matrix "<< this << " " << e_Mij <<" N = " << N << " M = " << M << endl;
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

                //V4
		A  a=GetAny<A>((*Mat)(stack));

		KN<long> *lg,*cl;
		KN<RR> *cc;
		lg = GetAny<KN<long>*>((*lig)(stack));
		cl = GetAny<KN<long>*>((*col)(stack));
		cc = GetAny<KN<RR>*>((*coef)(stack));
		int n=a->N(),m=a->M();
                HashMatrix<int,RR> *mh = a->pHM();
                mh->COO();

                int kk = mh->nnz,k1=0;
                if( mh->pij(n-1,m-1)==0 ) k1=1;
		lg->resize(kk+k1);
		cc->resize(kk+k1);
		cl->resize(kk+k1);
		int k=0;

                for ( k=0 ; k < kk;++k)
		  {
                      (*lg)[k]= mh->i[k];
		      (*cl)[k]= mh->j[k];
                      (*cc)[k]= mh->aij[k];
		  }

                if(k1)
                {
                 (*lg)[kk]= n;
                 (*cl)[kk]= m;
                 (*cc)[kk]= 0;
                }

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

template<typename R>  AnyType RawMatrix<R>::operator()(Stack stack) const
{
    MatriceMorse<R> * amorse =0;
    KN_<R> cc(GetAny< KN_<R>  >((*coef)(stack)));
    int k= cc.N();
    int n= k;
    int m=n;
    bool sym=false;
    if( lig && col)
    {
	KN_<long> lg(GetAny< KN_<long>  >((*lig)(stack)));
	KN_<long> cl=(GetAny< KN_<long>  >((*col)(stack)));
	n = lg.max()+1;
	m = cl.max()+1;
	ffassert( lg.N()==k && cl.N()==k && lg.min()>=0 && lg.max()>=0);
        amorse = new MatriceMorse<R>(n,m,k,0);
	sym=false;
	for(int i=0;i<k;++i)
	    (*amorse)[make_pair<int,int>((int)lg[i],(int)cl[i])]+=cc[i];
    }
    else
    {
        amorse = new MatriceMorse<R>(n,cc);
    }

    if(verbosity)
	cout << "  -- Raw Matrix    nxm  =" <<n<< "x" << m << " nb  none zero coef. " << amorse->nnz << endl;

    Matrice_Creuse<R> * sparse_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
    if( !init) sparse_mat->init();

    sparse_mat->A.master(amorse);
    sparse_mat->typemat=0; //(amorse->n == amorse->m) ? TypeSolveMat(TypeSolveMat::GMRES) : TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)

    if(verbosity>3) { cout << "  End Raw Matrix : " << endl;}

    return sparse_mat;
}
template<typename R>  AnyType BlockMatrix<R>::operator()(Stack s) const
{
  typedef list<tuple<R,MatriceCreuse<R> *,bool> > * L;
   KNM<L> Bij(N,M);
   KNM<KNM_<R> * > Fij(N,M);
   KNM<bool> cnjij(N,M);
   KNM<R> Rij(N,M); //  to sto

   cnjij = false;
   KN<long> Oi(N+1), Oj(M+1);
   if(verbosity>9) { cout << " Build Block Matrix : " << N << " x " << M << endl;}
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
    {  if(verbosity>99) cout << j << " colum size" << Oj(j+1) << endl;
        if   ( Oj(j+1) ==0) {
            Oj(j+1)=1;
            if( Oj(j+1) !=1)  err++;}
    }
    for (int i=0;i<N;++i)
    {
        if(verbosity>99) cout << i << " row size" << Oi(i+1) << endl;
        if   ( Oi(i+1) ==0) {
               Oi(i+1)=1;
               if( Oi(i+1) !=1)  err++;}
    }
    if (err)    ExecError("Error Block Matrix with  0 line or  0 colomn..");

    for (int i=0;i<N;++i)
      Oi(i+1) += Oi(i);
    for (int j=0;j<M;++j) // correct 10/01/2007 FH
      Oj(j+1) += Oj(j);// correct 07/03/2010 FH
  long n=Oi(N),m=Oj(M);
  if(verbosity>99)
   {
     cout << "     Oi = " <<  Oi << endl;
     cout << "     Oj = " <<  Oj << endl;
  }
  MatriceMorse<R> * amorse =0;
{
    HashMatrix<int,R>  *Aij = new  HashMatrix<int,R>( n, m,0,0);
    for (int i=0;i<N;++i)
     for (int j=0;j<M;++j)
       if (Bij(i,j))
         {
           if(verbosity>99)
             cout << "  Add  Block S " << i << "," << j << " =  at " << Oi(i) << " x " << Oj(j) << " conj = " << cnjij(i,j) << endl;
             HashMatrix<int,R> & mmij=*Aij;
             const list<tuple<R,MatriceCreuse<R>*,bool> >  &lM=*Bij(i,j);
             bool ttrans=false;
             int ii00=Oi(i);
             int jj00=Oj(j);
             bool cnj=cnjij(i,j);
            BuildCombMat(mmij,lM,ttrans,ii00,jj00,cnj);


         }
       else if (Fij(i,j))
        {
           if(verbosity>99)
             cout << "  Add  Block F " << i << "," << j << " =  at " << Oi(i) << " x " << Oj(j) << endl;
           BuildCombMat(*Aij,*Fij(i,j),Oi(i),Oj(j),R(1.),cnjij(i,j));// BuildCombMat
        }


    amorse=  Aij;
  }
  if(verbosity>9)
     cout << "  -- Block Matrix NxM = " << N << "x" << M << "    nxm  =" <<n<< "x" << m << " nb  none zero coef. " << amorse->nnz << endl;

  Matrice_Creuse<R> * sparse_mat =GetAny<Matrice_Creuse<R>* >((*emat)(s));
  if(!init) sparse_mat->init();
  sparse_mat->A.master(amorse);
    sparse_mat->typemat=0;//(amorse->n == amorse->m) ? TypeSolveMat(TypeSolveMat::GMRES) : TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)


  // cleanning
  for (int i=0;i<N;++i)
   for (int j=0;j<M;++j)
    if(Bij(i,j)) delete Bij(i,j);
    else if(Fij(i,j))  delete Fij(i,j);
   if(verbosity>9) { cout << "  End Build Blok Matrix : " << endl;}

 return sparse_mat;

}

template<class R>
class minusMat { public:
    list<tuple<R,MatriceCreuse<R> *,bool> >  *l;
    minusMat(list<tuple<R,MatriceCreuse<R> *,bool> > *ll):
	l(new list<tuple<R,MatriceCreuse<R> *,bool> >(*ll) )
      {
	    typedef typename list<tuple<R,MatriceCreuse<R> *,bool> >::iterator lci;
	    for (lci i= l->begin();i !=l->end();++i)
		get<0>(*i) *= R(-1);
      }
};

template<class R>
AnyType mM2L3 (Stack , const AnyType & pp)
{
    minusMat<R> mpp(to(GetAny<Matrice_Creuse<R> *>(pp)));
    return SetAny<minusMat<R> >(mpp);
}

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
        if(verbosity>1000) {
        cout << " i " << (char*) (void *) i -  (char*)(void*) s ;
        cout << " j " << (char*) (void *) j  -  (char*)(void*) s ;
        cout << " vi " <<  (char*) (void *) v -  (char*)(void*) s ;
        cout << endl;
        }

        ffassert(i && v);
        MatriceCreuse<R> *m=t->A;
        MatriceMorse<R> *mm = dynamic_cast<MatriceMorse<R>*>(m);
        if(!mm) ExecError(" Matrix sparse of bad type ( not HMatrix ) , sorry.. ");
        if(mm)
        for (long  kk=0;kk< mm->nnz; ++kk)
            {
                *i=mm->i[kk];
                *j= mm->j[kk];
                *v =  mm->aij[kk];
                data->code(s);
                mm->aij[kk] = *v;
            }
        return Nothing  ;
    }

};

//  Mat real -> mat complex .... ??? FH.   april  2016 ....
struct  VirtualMatCR :public RNM_VirtualMatrix<Complex>{ public:
   RNM_VirtualMatrix<double>& VM;
    typedef Complex R;
    VirtualMatCR( RNM_VirtualMatrix<double> & MM): RNM_VirtualMatrix<Complex>(MM.N,MM.M),  VM(MM) {}
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
    { if( !VM.WithSolver()) InternalError("RNM_VirtualMatrix::solve not implemented ");
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

            RNM_VirtualMatrix<Complex> *pv = new  VirtualMatCR ((*GetAny<A>((*a0)(s))).A);
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
// Norme Hashmatrix
template<class K>
double get_norme_linfty(Matrice_Creuse<K> * p){
    if(p==0) return 0.;
    HashMatrix<int,K>* ph=p->pHM();
    ffassert(ph);
    return ph->norminfty();}
template<class K>
double get_norme_l2(Matrice_Creuse<K> * p){
    if(p==0) return 0.;
    HashMatrix<int,K>* ph=p->pHM();
    ffassert(ph);
    return ph->FrobeniusNorm();}

template<class K>
double get_norme_l1(Matrice_Creuse<K> * p){
    if(p==0) return 0.;
    HashMatrix<int,K>* ph=p->pHM();
    ffassert(ph);
    return ph->norm1();}

template<class R,int Init>
Matrice_Creuse<R> * SetMatrice_Creuse(Matrice_Creuse<R> * p,newpMatrice_Creuse<R>  np)
{
    return np.set(p,Init);
}
//  Add F.H July 2019
template<class R>
Matrice_Creuse<R> * InitMatrice_Creuse_nm(Matrice_Creuse<R> * const & p,const long &n,const long &m)
{
    p->init() ;
    HashMatrix<int,R> *phm= new HashMatrix<int,R>((int) n,(int) m,0,0);
    MatriceCreuse<R> *pmc(phm);
    p->A.master(pmc);
    return p;
}
template<class R>
Matrice_Creuse<R> * InitMatrice_Creuse_n(Matrice_Creuse<R> * const & p,const long &n)
{
    p->init() ;
    HashMatrix<int,R> *phm= new HashMatrix<int,R>((int)n,(int) n,0,0);
    MatriceCreuse<R> *pmc(phm);
    p->A.master(pmc);
    return p;
}

template<class R,int c>
Matrice_Creuse<R> * AddtoMatrice_Creuse(Matrice_Creuse<R> * p,newpMatrice_Creuse<R>  np)
{
    return np.add(p,double(c));
}

template<class K, bool init>
Matrice_Creuse<K>* set_H_Eye(Matrice_Creuse<K> *pA,const  Eye eye)
{
    int n = eye.n, m=eye.m, nn= min(n,m);
    if( init) pA->init();
    pA->resize(n,m);
    HashMatrix<int,K> * pH= pA->pHM();
    ffassert(pH);
    pH->clear();
    pH->resize(n,m,nn);
    for(int i=0; i< n; ++i)
        (*pH)(i,i)=1.;
    return  pA;
}
template <class R>
void AddSparseMat()
{
 SetMatrix_Op<R>::btype = Dcl_Type<const  SetMatrix_Op<R> * >();
 Dcl_Type<TheDiagMat<R> >();
 Dcl_Type<TheCoefMat<R> >(); // Add FH oct 2005
 Dcl_Type< map< pair<int,int>, R> * >(); // Add FH mars 2005
 Dcl_Type<  minusMat<R>  >(); // Add FJH mars 2007

 basicForEachType * t_MC=atype<  Matrice_Creuse<R>* >();

 t_MC->SetTypeLoop(atype<  R* >(),atype<  long* >(),atype<  long* >());

 basicForEachType * t_MM=atype<map< pair<int,int>, R> * >();

TheOperators->Add("*",
        new OneBinaryOperator<Op2_mulvirtAv<typename RNM_VirtualMatrix<R>::plusAx,Matrice_Creuse<R>*,KN_<R> > >,
        new OneBinaryOperator<Op2_mulvirtAv<typename RNM_VirtualMatrix<R>::plusAtx,Matrice_Creuse_Transpose<R>,KN_<R> > >,
        new OneBinaryOperator<Op2_mulvirtAv<typename RNM_VirtualMatrix<R>::solveAxeqb,Matrice_Creuse_inv<R>,KN_<R> > > ,
        new OneBinaryOperator<Op2_mulvirtAv<typename RNM_VirtualMatrix<R>::solveAtxeqb,Matrice_Creuse_inv_trans<R>,KN_<R> > >
        );

TheOperators->Add("^", new OneBinaryOperatorA_inv<R>());
TheOperators->Add("^", new OneBinaryOperatorAt_inv<R>());

// matrix new code   FH (Houston 2004)
 TheOperators->Add("=",
        new OneOperator2<Matrice_Creuse<R>*,Matrice_Creuse<R>*,newpMatrice_Creuse<R> > (SetMatrice_Creuse<R,0> ),
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const Matrix_Prod<R,R>,E_F_StackF0F0>(ProdMat<R,R,R,1>),
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,KN<R> *,E_F_StackF0F0>(DiagMat<R,1>),
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,Matrice_Creuse_Transpose<R>,E_F_StackF0F0>(CopyTrans<R,R,1>),
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,Matrice_Creuse<R>*,E_F_StackF0F0>(CopyMat<R,R,1>) ,
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,KNM<R>*,E_F_StackF0F0>(MatFull2Sparse<R,1>) ,
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,list<tuple<R,MatriceCreuse<R> *,bool> > *,E_F_StackF0F0>(CombMat<R,1>) ,
       new OneOperatorCode<BlockMatrix1<R> >(),
       new OneOperator2<Matrice_Creuse<R>*,Matrice_Creuse<R>*,Eye>(set_H_Eye<R,false> )

       );
    TheOperators->Add("+=",
        new OneOperator2<Matrice_Creuse<R>*,Matrice_Creuse<R>*,newpMatrice_Creuse<R> > (AddtoMatrice_Creuse<R, 1> ),
        new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,list<tuple<R,MatriceCreuse<R> *,bool> > *,E_F_StackF0F0>(AddCombMat<R,1>));

    TheOperators->Add("-=",
                      new OneOperator2<Matrice_Creuse<R>*,Matrice_Creuse<R>*,newpMatrice_Creuse<R> > (AddtoMatrice_Creuse<R, -1> ),
                      new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,list<tuple<R,MatriceCreuse<R> *,bool> > *,E_F_StackF0F0>(AddCombMat<R,-1>));



 TheOperators->Add("<-",
       new OneOperatorCode<BlockMatrix<R> >(),
       new OneOperator2<Matrice_Creuse<R>*,Matrice_Creuse<R>*,newpMatrice_Creuse<R> > (SetMatrice_Creuse<R,1> ),
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,long > (InitMatrice_Creuse_n<R> ),
       new OneOperator3_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,long,long > (InitMatrice_Creuse_nm<R> ),
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const Matrix_Prod<R,R>,E_F_StackF0F0>(ProdMat<R,R,R,0>),
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,KN<R> *,E_F_StackF0F0>(DiagMat<R,0>)  ,
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,Matrice_Creuse_Transpose<R>,E_F_StackF0F0>(CopyTrans<R,R,0>),
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,Matrice_Creuse<R>*,E_F_StackF0F0>(CopyMat<R,R,0>) ,
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,KNM<R>*,E_F_StackF0F0>(MatFull2Sparse<R,0>) ,
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,list<tuple<R,MatriceCreuse<R> *,bool> > *,E_F_StackF0F0>(CombMat<R,0>),
       new OneOperator2<Matrice_Creuse<R>*,Matrice_Creuse<R>*,Eye>(set_H_Eye<R,true> )


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
        new OneBinaryOperator<Op2_ListMMadd<R> >

       );
    TheOperators->Add("-",
                      new OneBinaryOperator<Op2_ListCMCMsub<R> >,   //  (L) - (L)
                      new OneBinaryOperator<Op2_ListCMMadd<R,-1> >, // L - M
                      new OneBinaryOperator<Op2_ListMCMadd<R,-1> >, // M - L
                      new OneBinaryOperator<Op2_ListMMadd<R,-1> >   // M - M

                      );
 TheOperators->Add("-",
	 new OneUnaryOperator<Op1_LCMd<R,-1> >
     );
TheOperators->Add("+",
                      new OneUnaryOperator<Op1_LCMd<R,1> >
                      );

    Add<Matrice_Creuse<R> *>("COO",".",new OneOperator1<bool,Matrice_Creuse<R> *>(set_mat_COO<R>) );
    Add<Matrice_Creuse<R> *>("CSR",".",new OneOperator1<bool,Matrice_Creuse<R> *>(set_mat_CSR<R>) );
    Add<Matrice_Creuse<R> *>("CSC",".",new OneOperator1<bool,Matrice_Creuse<R> *>(set_mat_CSC<R>) );
    Add<Matrice_Creuse<R> *>("n",".",new OneOperator1<long,Matrice_Creuse<R> *>(get_mat_n<R>) );
 Add<Matrice_Creuse<R> *>("m",".",new OneOperator1<long,Matrice_Creuse<R> *>(get_mat_m<R>) );
 Add<Matrice_Creuse<R> *>("nbcoef",".",new OneOperator1<long,Matrice_Creuse<R> *>(get_mat_nbcoef<R>) );
 Add<Matrice_Creuse<R> *>("nnz",".",new OneOperator1<long,Matrice_Creuse<R> *>(get_mat_nbcoef<R>) );
 Add<Matrice_Creuse<R> *>("half",".",new OneOperator1<long,Matrice_Creuse<R> *>(get_mat_half<R>) );
 Add<Matrice_Creuse<R> *>("size",".",new OneOperator1<long,Matrice_Creuse<R> *>(get_mat_nbcoef<R>) );
 Add<Matrice_Creuse<R> *>("trace",".",new OneOperator1<R,Matrice_Creuse<R>* >(get_trace_mat<R>) );
 Add<Matrice_Creuse<R> *>("clear",".",new OneOperator1<bool,Matrice_Creuse<R>* >(clear_mat<R>) );

 Add<Matrice_Creuse<R> *>("diag",".",new OneOperator1<TheDiagMat<R> ,Matrice_Creuse<R> *>(thediag<R>) );
 Add<Matrice_Creuse<R> *>("coef",".",new OneOperator1<TheCoefMat<R> ,Matrice_Creuse<R> *>(thecoef<R>) );

 TheOperators->Add("=", new OneOperator2<KN<R>*,KN<R>*,TheDiagMat<R> >(get_mat_daig<R>) );

 TheOperators->Add("<-", new OneOperator2<KN<R>*,KN<R>*,TheDiagMat<R> >(init_get_mat_daig<R>) );
    TheOperators->Add("=", new OneOperator2<TheDiagMat<R>,TheDiagMat<R>,KN<R>*>(set_mat_daig<R>) );

// ADD oct 2005
 TheOperators->Add("=", new OneOperator2<KN<R>*,KN<R>*,TheCoefMat<R> >(get_mat_coef<R>) );
 TheOperators->Add("=", new OneOperator2<TheCoefMat<R>,TheCoefMat<R>,KN<R>*>(set_mat_coef<R>) );
 TheOperators->Add("=", new OneOperator2<TheCoefMat<R>,TheCoefMat<R>,R>(set_mat_coef<R>) );

 Global.Add("set","(",new SetMatrix<R>);
 Add<Matrice_Creuse<R> *>("linfty",".",new OneOperator1<double,Matrice_Creuse<R> *>(get_norme_linfty));
 Add<Matrice_Creuse<R> *>("l2",".",new OneOperator1<double,Matrice_Creuse<R> *>(get_norme_l2));
 Add<Matrice_Creuse<R> *>("l1",".",new OneOperator1<double,Matrice_Creuse<R> *>(get_norme_l1));
 atype<Matrice_Creuse<R> * >()->Add("(","",new OneOperator3_<R*,Matrice_Creuse<R> *,long,long >(1,get_elementp2mc<R>));

 atype<Matrice_Creuse<R> * >()->Add("[","",new OneOperator3_<R,Matrice_Creuse<R> *,long,long >(10,get_element2mc<R>));

//    typedef map< pair<int,int>, R> MAPMAT;
 typedef Matrice_Creuse<R> * MAPMATC;
 typedef   newpMatrice_Creuse<R> MAPMATN;
 atype<KNM<R>*>()->Add("(","",new OneOperator3s_<MAPMATN ,KNM<R>*,Inv_KN_long,Inv_KN_long >(Matrixfull2mapIJ_inv<R>));
 atype<KNM<R>*>()->Add("(","",new OneOperator3s_<MAPMATN ,KNM<R>*,KN_<long>,KN_<long> >(Matrixfull2mapIJ<R>));

 atype<outProduct_KN_<R>*>()->Add("(","",new OneOperator3s_<MAPMATN ,outProduct_KN_<R>*,Inv_KN_long,Inv_KN_long >(Matrixoutp2mapIJ_inv<R>));
 atype<outProduct_KN_<R>*>()->Add("(","",new OneOperator3s_<MAPMATN ,outProduct_KN_<R>*,KN_<long>,KN_<long> >(Matrixoutp2mapIJ<R>));


 TheOperators->Add("=", new SetRawMatformMat<R>);



 t_MM->Add("(","",new OneOperator3s_<MAPMATN ,MAPMATC ,Inv_KN_long,Inv_KN_long >(Matrixmapp2mapIJ1<R>));
 t_MM->Add("(","",new OneOperator3s_<MAPMATN ,MAPMATC ,KN_<long>,KN_<long> >(Matrixmapp2mapIJ<R>));

 t_MC->Add("(","",new OneOperator3s_<MAPMATN ,MAPMATC ,Inv_KN_long,Inv_KN_long >(Matrixmapp2mapIJ1<R>,t_MC));
 t_MC->Add("(","",new OneOperator3s_<MAPMATN ,MAPMATC ,KN_<long>,KN_<long> >(Matrixmapp2mapIJ<R>,t_MC));

 map_type[typeid(MAPMATN ).name()]->AddCast(
     new E_F1_funcT<MAPMATN ,KNM<R>* >(Matrixfull2map<R>),
     new E_F1_funcT<MAPMATN ,outProduct_KN_<R>* >(Matrixoutp2map<R>)

       );

 map_type[typeid(list<tuple<R,MatriceCreuse<R> *,bool> > *).name()]->AddCast(
     new E_F1_funcT<list<tuple<R,MatriceCreuse<R> *,bool> > *,Matrice_Creuse<R>* >(M2L3<R>),
     new E_F1_funcT<list<tuple<R,MatriceCreuse<R> *,bool> > *,Matrice_Creuse_Transpose<R> >(tM2L3<R>),
     new E_F1_funcT<list<tuple<R,MatriceCreuse<R> *,bool> > *,minusMat<R> >(mM2L3<R> )
     );


    TheOperators->Add("{}",new ForAllLoop<E_ForAllLoopMatrix<R> >);
   // remove 2 plugin thresholding and symmetrizeCSR
    typedef Thresholding<R> TMR;
    typedef Matrice_Creuse<R> MR;
     Dcl_Type<TMR>();
    TMR t(0);
    //thresholding2(t, 0.);
    Add<MR *>("thresholding", ".", new OneOperator1<TMR, MR *>(to_Thresholding));
    Add<TMR>("(", "", new OneOperator2_<MR *, TMR, double>(thresholding2));
    Global.Add("symmetrizeCSR", "(", new OneOperator1_<long, Matrice_Creuse<R> *>(symmetrizeCSR<R> ));
//  --- end
}

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

template<class T>
class removeDOF_Op : public E_F0mps {
public:
    Expression eA;
    Expression eR;
    Expression ex;
    Expression eout;
    bool transpose;
    static const int n_name_param = 4;
    static basicAC_F0::name_and_type name_param[];
    Expression nargs[n_name_param];
    removeDOF_Op(const basicAC_F0&  args, Expression param1, Expression param2, Expression param3, Expression param4)
    : eA(param1), eR(param2), ex(param3), eout(param4), transpose(false) {
        args.SetNameParam(n_name_param, name_param, nargs);
    }
    removeDOF_Op(const basicAC_F0&  args, Expression param2, Expression param3, Expression param4, bool t, int)
    : eA(0), eR(param2), ex(param3), eout(param4), transpose(t) {
        args.SetNameParam(n_name_param, name_param, nargs);
    }
    removeDOF_Op(const basicAC_F0&  args, Expression param1, Expression param2)
    : eA(param1), eR(param2), ex(0), eout(0), transpose(false) {
        args.SetNameParam(n_name_param, name_param, nargs);
    }

    AnyType operator()(Stack stack) const;
};

template<class T>
basicAC_F0::name_and_type removeDOF_Op<T>::name_param[] = {
    {"symmetrize", &typeid(bool)},
    {"condensation", &typeid(KN<long>*)},
    {"R", &typeid(Matrice_Creuse<double>*)},
    {"eps",&typeid(double)}
};

template<class T>
class removeDOF : public OneOperator {
   const unsigned short withA;
public:
    removeDOF() : OneOperator(atype<long>(), atype<Matrice_Creuse<T>*>(), atype<Matrice_Creuse<double>*>(), atype<KN<T>*>(), atype<KN<T>*>()),withA(1) {}
    removeDOF(int) : OneOperator(atype<long>(), atype<Matrice_Creuse<double>*>(), atype<KN<T>*>(), atype<KN<T>*>()),withA(0) {}
    removeDOF(int, int) : OneOperator(atype<long>(), atype<Matrice_Creuse<T>*>(), atype<Matrice_Creuse<double>*>()),withA(2) {}
    removeDOF(int, int, int) : OneOperator(atype<long>(), atype<Matrice_Creuse_Transpose<double>>(), atype<KN<T>*>(), atype<KN<T>*>()),withA(3) {}

    E_F0* code(const basicAC_F0& args) const {
        if(withA == 1)
        return new removeDOF_Op<T>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]), t[3]->CastTo(args[3]));
        else if(withA == 0 || withA == 3)
        return new removeDOF_Op<T>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]), withA == 3, 1);
        else
        return new removeDOF_Op<T>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]));
    }
};
template<class T> bool cmp(const std::pair<unsigned int, T>& lhs, const std::pair<unsigned int, T>& rhs) { return lhs.first < rhs.first; }
template<class T>
AnyType removeDOF_Op<T>::operator()(Stack stack)  const {

    static const double defEPS=1e-12;
    typedef double R;
    // code wri-ng no ...
    Matrice_Creuse<T>* pA = eA ? GetAny<Matrice_Creuse<T>* >((*eA)(stack)):0;
    Matrice_Creuse<R>* pR;
    if(transpose) {
        Matrice_Creuse_Transpose<R> tpR = GetAny<Matrice_Creuse_Transpose<R>>((*eR)(stack));
        pR = tpR;
    }
    else
        pR = GetAny<Matrice_Creuse<R>*>((*eR)(stack));
    KN<T>* pX = ex ? GetAny<KN<T>* >((*ex)(stack)) : 0;
    KN<T>* pOut = eout ? GetAny<KN<T>* >((*eout)(stack)) : 0;
    Matrice_Creuse<R>* pC = nargs[2] ? GetAny<Matrice_Creuse<R>*>((*nargs[2])(stack)) : 0;
    MatriceMorse<R> *mC = 0;
    if(pC && pC->A) {
        mC = static_cast<MatriceMorse<R>*>(&(*pC->A));
    }
    ffassert(pR);
    bool rhs = (pX && pOut) && (pOut->n > 0 || pX->n > 0);
    if(pA)
    {
        if(!pC)
            pC = pR;
        MatriceMorse<T> *mA = pA->pHM();
        MatriceMorse<R> *mR = pR->pHM();
        if(!mC)
            mC = mR;
        pA->Uh = pR->Uh;
        pA->Vh = pC->Vh;

        bool symmetrize = nargs[0] ? GetAny<bool>((*nargs[0])(stack)) : false;
        double EPS=nargs[3] ? GetAny<double>((*nargs[3])(stack)) :defEPS ;
        KN<long>* condensation = nargs[1] ? GetAny<KN<long>* >((*nargs[1])(stack)) : (KN<long>*) 0;
        ffassert(condensation ||( mR && mC) );
        unsigned int n = condensation ? condensation->n : mR->nnz;
        unsigned int m = condensation ? condensation->n : mC->nnz;
        KN<int> lg(n+1,0);

        if(rhs && pOut->n != n) pOut->resize(n);
        mC->COO();
        mR->COO();
        mA->CSR();
        std::vector<signed int> tmpVec;
        if(!condensation)
        {
            tmpVec.resize(mA->m);
            for(long i = 0; i < m; ++i)
                tmpVec[mC->j[i]] = i + 1;
            if(!mA->half) {
                std::vector<std::pair<int, T> > tmp;
                tmp.reserve(mA->nnz);

                lg[0] = 0;
                for(long i = 0; i < n; ++i) {
                    for(long j = mA->p[mR->j[i]]; j < mA->p[mR->j[i] + 1]; ++j) {
                        long col = tmpVec[mA->j[j]];
                        if(col != 0 && abs(mA->aij[j]) > EPS) {
                            if(symmetrize) {
                                if(col - 1 <= i)
                                    tmp.push_back(std::make_pair(col - 1, mA->aij[j]));
                            }
                            else
                                tmp.push_back(std::make_pair(col - 1, mA->aij[j]));

                        }
                    }
                    std::sort(tmp.begin() + lg[i], tmp.end(),cmp<T> );
                    // c++11 , [](const std::pair<unsigned int, T>& lhs, const std::pair<unsigned int, T>& rhs) { return lhs.first < rhs.first; });
                    if(rhs)
                        *(*pOut + i) = *(*pX + mC->j[i]);
                    lg[i + 1] = tmp.size();
                }
                mA->clear();
                mA->resize(n,m);
                MatriceMorse<T> &A = *mA;
                A.half = symmetrize;
                for(int i=0; i<n; ++i)
                {
                    for(int k= lg[i]; k < lg[i+1]; ++k)
                    {
                        int j= tmp[k].first;
                        T aij  = tmp[k].second;
                        A(i,j) =aij;
                    }
                }


            }// !mA->Half
            else {
                std::vector<std::vector<std::pair<unsigned int, T> > > tmp(n);
                for(unsigned int i = 0; i < n; ++i)
                    tmp[i].reserve(mA->p[mR->j[i] + 1] - mA->p[mR->j[i]]);

                unsigned int nnz = 0;
                for(unsigned int i = 0; i < n; ++i) {
                    for(unsigned int j = mA->p[mR->j[i]]; j < mA->p[mR->j[i] + 1]; ++j) {
                        unsigned int col = tmpVec[mA->j[j]];
                        if(col != 0 && abs(mA->aij[j]) > EPS) {
                            if(i < col - 1)
                                tmp[col - 1].push_back(make_pair(i, mA->aij[j]));
                            else
                                tmp[i].push_back(make_pair(col - 1, mA->aij[j]));
                            ++nnz;
                        }
                    }
                    if(rhs)
                        *(*pOut + i) = *(*pX + mC->j[i]);
                }
                int Half = mA->half;
                mA->clear();
                mA->resize(n,m,nnz);
                MatriceMorse<T> &A = *mA;
                A.half = Half;
                for(unsigned int i = 0; i < n; ++i) {
                    std::sort(tmp[i].begin(), tmp[i].end(),cmp<T>);
                    // c++11, [](const std::pair<unsigned int, T>& lhs, const std::pair<unsigned int, T>& rhs) { return lhs.first < rhs.first; });
                    for(typename std::vector<std::pair<unsigned int, T> >::const_iterator it = tmp[i].begin(); it != tmp[i].end(); ++it)
                          A(i,it->first) =it->second;

                }

            }
        }

        else
        {
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

            std::vector<std::pair<int, T> > tmpInterior;
            std::vector<std::pair<int, T> > tmpBoundary;
            std::vector<std::pair<int, T> > tmpInteraction;
            tmpInterior.reserve(mA->nnz);
            tmpBoundary.reserve(mA->nnz);
            tmpInteraction.reserve(mA->nnz);

            lg[0] = 0;
            for(long i = 0; i < mA->n; ++i) {
                int row = tmpVec[i];
                if(row < 0) {
                    for(unsigned int j = mA->p[i]; j < mA->p[i + 1]; ++j) {
                        int col = tmpVec[mA->j[j]];
                        if(col < 0)
                            tmpInterior.push_back(make_pair(-col - 1, mA->aij[j]));
                        else
                            tmpInteraction.push_back(make_pair(col - 1, mA->aij[j]));
                    }

                }
                else {
                    for(unsigned int j = mA->p[i]; j < mA->p[i + 1]; ++j) {
                        int col = tmpVec[mA->j[j]];
                        if(col > 0)
                            tmpBoundary.push_back(make_pair(col - 1, mA->aij[j]));
                    }
                    if(rhs)
                        *(*pOut + i) = *(*pX + *(*condensation + i));
                    lg[i + 1] = tmpBoundary.size();
                }
            }

            mA->clear();
            mA->resize(n,n);
            MatriceMorse<T> &mR = *new MatriceMorse<T>(n,m,tmpBoundary.size(),0);
            for(unsigned int i = 0; i < tmpBoundary.size(); ++i) {
                mR(i, tmpBoundary[i].first)= tmpBoundary[i].second;
            }
            ffassert(0);
            pR->typemat = 0; //TypeSolveMat(TypeSolveMat::GMRES);
            // bug ici ::: FH..            pR->A.master(&mR);
        }

    }
    else if(rhs)
    {


       MatriceMorse<R> *mR = pR->pHM();
       if(mR && mR->n && mR->m) {
        mR->COO();
        unsigned int n = mR->nnz;
        if(transpose) {
            if(pOut->n != mR->m) pOut->resize(mR->m);
            *pOut = T();
            for(unsigned int i = 0; i < n; ++i) {
                *(*pOut + mR->j[i]) = *(*pX + mR->i[i]);
            }
        }
        else {
            if(pOut->n != n) pOut->resize(n);
            for(unsigned int i = 0; i < n; ++i) {
                *(*pOut + i) = *(*pX + mR->j[i]);
            }
        }
       }
       else
           *pOut = T();
    }
    return 0L;
}


bool SparseDefault()
{
    return 1;//TypeSolveMat::SparseSolver== TypeSolveMat::defaultvalue;
}

bool Have_UMFPACK_=false;
bool Have_UMFPACK() { return Have_UMFPACK_;}

MatriceMorse<R> * removeHalf(MatriceMorse<R> & A,long half,double tol)
{

    // half < 0 => L
    // half > 0 => U
    // half = 0 => L and the result will be sym
    int nnz =0;

    if( A.half )
        return &A;//  copy
    // do alloc
    MatriceMorse<R> *r=new MatriceMorse<R>(A);
    r->RemoveHalf(half,tol);
    if(verbosity )
        cout << "  removeHalf: new nnz = "<< r->nnz << " "<< r->half << endl;

    return r;
}

newpMatrice_Creuse<R> removeHalf(Stack stack,Matrice_Creuse<R> *const & pA,long const & half,const double & tol)
{
    MatriceCreuse<R> * pa=pA->A;
    MatriceMorse<R> *pma= dynamic_cast<MatriceMorse<R>* > (pa);
    ffassert(pma);
    return newpMatrice_Creuse<R>(stack,removeHalf(*pma,half,tol));
}

bool removeHalf(Stack stack,Matrice_Creuse<R> *const & pR,Matrice_Creuse<R> *const & pA,long const & half,const double & tol)
{
    MatriceCreuse<R> * pa=pA->A;
    MatriceMorse<R> *pma= dynamic_cast<MatriceMorse<R>* > (pa);
    MatriceCreuse<R> * pr= removeHalf(*pma,half,tol);

    pR->A.master(pr);
    return true;
}
newpMatrice_Creuse<R> removeHalf(Stack stack,Matrice_Creuse<R> *const & pA,long const & half)
{
    return removeHalf(stack,pA,half,-1.);
}

template<class K>
class plotMatrix : public OneOperator {
public:

	class Op : public E_F0info {
	public:
		Expression a;

		static const int n_name_param = 1;
		static basicAC_F0::name_and_type name_param[] ;
		Expression nargs[n_name_param];
		bool arg(int i,Stack stack,bool a) const{ return nargs[i] ? GetAny<bool>( (*nargs[i])(stack) ): a;}
		long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}

	public:
		Op(const basicAC_F0 &  args,Expression aa) : a(aa) {
			args.SetNameParam(n_name_param,name_param,nargs);
		}

            AnyType operator()(Stack stack) const{

                if (mpirank == 0 && ThePlotStream) {
                    Matrice_Creuse<K>* A =GetAny<Matrice_Creuse<K>* >((*a)(stack));
                    bool wait = arg(0,stack,false);

                    PlotStream theplot(ThePlotStream);
                    theplot.SendNewPlot();
                    theplot << 3L;
                    theplot <= wait;
                    theplot.SendEndArgPlot();
                    theplot.SendMeshes();
                    theplot << 0L;
                    theplot.SendPlots();
                    theplot << 1L;
                    theplot << 31L;

                    HashMatrix<int,K>* ph=A->pHM();

                    if (!ph) {
                        theplot << 0;
                        theplot << 0;
                        theplot << 0L;
                        theplot << 0L;
                    }
                    else {
                        theplot << (int)ph->n;
                        theplot << (int)ph->m;
                        theplot << 0L;
                        theplot << (long)ph->nnz;

                        for (int i=0;i<ph->nnz;i++) {
                            theplot << ph->i[i];
                            theplot << ph->j[i];
                            theplot << 1;
                            theplot << 1;
                            theplot << 1;
                            theplot << abs(ph->aij[i]);
                        }
                    }

                    theplot.SendEndPlot();

                }

                return 0L;
            }
        };

    plotMatrix() : OneOperator(atype<long>(),atype<Matrice_Creuse<K>*>()) {}

    E_F0 * code(const basicAC_F0 & args) const
    {
        return  new Op(args,t[0]->CastTo(args[0]));
    }
};

template<class K>
basicAC_F0::name_and_type  plotMatrix<K>::Op::name_param[]= {
	{  "wait", &typeid(bool)}
};

void  init_lgmat()

{

  Dcl_Type<const  MatrixInterpolation<pfes,pfes>::Op *>();
  Dcl_Type<const  MatrixInterpolation<pfes3,pfes3>::Op *>();
  Dcl_Type<const  MatrixInterpolation<pfesS,pfesS>::Op *>();
  Dcl_Type<const  MatrixInterpolation<pfesL,pfesL>::Op *>();
  Dcl_Type<const  MatrixInterpolation<pfesS,pfes3>::Op *>();
  Dcl_Type<const  MatrixInterpolation<pfesL,pfes>::Op *>();
  Dcl_Type<const  MatrixInterpolation<pfesL,pfesS>::Op *>();
  Dcl_Type<const  MatrixInterpolation<pfesS,pfes>::Op *>();
  Dcl_Type<const  MatrixInterpolation<pfes,pfesL>::Op *>();
  Dcl_Type<const  MatrixInterpolation<pfes,pfesS>::Op *>();
  Dcl_Type<const  MatrixInterpolation<pfes,pfes3>::Op *>();
    
  map_type_of_map[make_pair(atype<Matrice_Creuse<double>* >(),atype<double*>())]=atype<Matrice_Creuse<double> *>();
  map_type_of_map[make_pair(atype<Matrice_Creuse<double>* >(),atype<Complex*>())]=atype<Matrice_Creuse<Complex> *>();
  AddSparseMat<double>();
  AddSparseMat<Complex>();

  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfes,pfes>);
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfes,pfes>(1));
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfes3,pfes3>);
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfes3,pfes3>(1,1));
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfesS,pfesS>);
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfesS,pfesS>(1,1));
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfesL,pfesL>);
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfesL,pfesL>(1,1));
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfesS,pfes3>);
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfesS,pfes3>(1,1));
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfesL,pfes>);
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfesL,pfes>(1,1));
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfesL,pfesS>);
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfesL,pfesS>(1,1));
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfesS,pfes>);
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfesS,pfes>(1,1));
 
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfes,pfesL>);
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfes,pfesL>(1,1));
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfes,pfesS>);
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfes,pfesS>(1,1));
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfes,pfes3>);
  Add<const MatrixInterpolation<pfes,pfes>::Op *>("<-","(", new MatrixInterpolation<pfes,pfes3>(1,1));
    
    
    Dcl_Type<const  RestrictArray<pfes>::Op *>();
    Dcl_Type<const  RestrictArray<pfes3>::Op *>();
    Dcl_Type<const  RestrictArray<pfesS>::Op *>();
    Dcl_Type<const  RestrictArray<pfesL>::Op *>();

  Global.Add("restrict","(",new RestrictArray<pfes>);// FH Jan 2014
  Global.Add("restrict","(",new RestrictArray<pfes3>);// FH Jan 2014
  Global.Add("restrict","(",new RestrictArray<pfesS>);// PHT Apr 2019
  Global.Add("restrict","(",new RestrictArray<pfesL>);

  TheOperators->Add("=",
                      new OneOperator2_<KN<long>*,KN<long>*,const RestrictArray<pfes>::Op*,E_F_StackF0F0>(SetRestrict<pfes,1>),
                      new OneOperator2_<KN<long>*,KN<long>*,const RestrictArray<pfes3>::Op*,E_F_StackF0F0>(SetRestrict<pfes3,1>),
                      new OneOperator2_<KN<long>*,KN<long>*,const RestrictArray<pfesS>::Op*,E_F_StackF0F0>(SetRestrict<pfesS,1>),
                      new OneOperator2_<KN<long>*,KN<long>*,const RestrictArray<pfesL>::Op*,E_F_StackF0F0>(SetRestrict<pfesL,1>)
                      );
    TheOperators->Add("<-",
                      new OneOperator2_<KN<long>*,KN<long>*,const RestrictArray<pfes>::Op*,E_F_StackF0F0>(SetRestrict<pfes,0>),
                      new OneOperator2_<KN<long>*,KN<long>*,const RestrictArray<pfes3>::Op*,E_F_StackF0F0>(SetRestrict<pfes3,0>),
                      new OneOperator2_<KN<long>*,KN<long>*,const RestrictArray<pfesS>::Op*,E_F_StackF0F0>(SetRestrict<pfesS,0>),
                      new OneOperator2_<KN<long>*,KN<long>*,const RestrictArray<pfesL>::Op*,E_F_StackF0F0>(SetRestrict<pfesL,0>)
                      );


  Global.Add("interpolate","(",new MatrixInterpolation<pfes,pfes>);
  Global.Add("interpolate","(",new MatrixInterpolation<pfes,pfes>(1));
  Global.Add("interpolate","(",new MatrixInterpolation<pfes3,pfes3>);
  Global.Add("interpolate","(",new MatrixInterpolation<pfes3,pfes3>(1,1));
  Global.Add("interpolate","(",new MatrixInterpolation<pfesS,pfesS>);
  Global.Add("interpolate","(",new MatrixInterpolation<pfesS,pfesS>(1,1));
  Global.Add("interpolate","(",new MatrixInterpolation<pfesL,pfesL>);
  Global.Add("interpolate","(",new MatrixInterpolation<pfesL,pfesL>(1,1));
    
  Global.Add("interpolate","(",new MatrixInterpolation<pfesS,pfes3>);
  Global.Add("interpolate","(",new MatrixInterpolation<pfesS,pfes3>(1,1));
  Global.Add("interpolate","(",new MatrixInterpolation<pfesL,pfes>);
  Global.Add("interpolate","(",new MatrixInterpolation<pfesL,pfes>(1,1));
  Global.Add("interpolate","(",new MatrixInterpolation<pfesL,pfesS>);
  Global.Add("interpolate","(",new MatrixInterpolation<pfesL,pfesS>(1,1));
  Global.Add("interpolate","(",new MatrixInterpolation<pfesS,pfes>);
  Global.Add("interpolate","(",new MatrixInterpolation<pfesL,pfes>(1,1));
  
  Global.Add("interpolate","(",new MatrixInterpolation<pfes,pfesL>);
  Global.Add("interpolate","(",new MatrixInterpolation<pfes,pfesL>(1,1));
  Global.Add("interpolate","(",new MatrixInterpolation<pfes,pfesS>);
  Global.Add("interpolate","(",new MatrixInterpolation<pfes,pfesS>(1,1));
  Global.Add("interpolate","(",new MatrixInterpolation<pfes,pfes3>);
  Global.Add("interpolate","(",new MatrixInterpolation<pfes,pfes3>(1,1));
    
  Global.Add("interplotematrix","(",new  OneOperatorCode<PrintErrorCompileIM>);
  zzzfff->Add("mapmatrix",atype<map< pair<int,int>, double> *>());
  zzzfff->Add("Cmapmatrix",atype<map< pair<int,int>, Complex> *>()); // a voir

  Dcl_Type< Resize<Matrice_Creuse<double> > > ();

  Add<Matrice_Creuse<double> *>("resize",".",new OneOperator1< Resize<Matrice_Creuse<double>  >,Matrice_Creuse<double> *>(to_Resize));
  Add<Resize<Matrice_Creuse<double> > >("(","",new OneOperator3_<Matrice_Creuse<double>  *,Resize<Matrice_Creuse<double>  > , long, long  >(resize2));
  // add missing in
 Dcl_Type< Resize<Matrice_Creuse<Complex> > > ();
 Add<Matrice_Creuse<Complex> *>("resize",".",new OneOperator1< Resize<Matrice_Creuse<Complex>  >,Matrice_Creuse<Complex> *>(to_Resize));
 Add<Resize<Matrice_Creuse<Complex> > >("(","",new OneOperator3_<Matrice_Creuse<Complex>  *,Resize<Matrice_Creuse<Complex>  > , long, long  >(resize2));


 // pour compatibiliter

 TheOperators->Add("=",
		   new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfes,pfes>::Op*,E_F_StackF0F0>(SetMatrixInterpolation<1>),
		   new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfes3,pfes3>::Op*,E_F_StackF0F0>(SetMatrixInterpolation3<1>),
           new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfesS,pfesS>::Op*,E_F_StackF0F0>(SetMatrixInterpolationS<1>),
           new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfesL,pfesL>::Op*,E_F_StackF0F0>(SetMatrixInterpolationL<1>),
           new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfesL,pfes>::Op*,E_F_StackF0F0>(SetMatrixInterpolationL2<1>),
           new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfesL,pfesS>::Op*,E_F_StackF0F0>(SetMatrixInterpolationLS<1>),
           new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfesS,pfes3>::Op*,E_F_StackF0F0>(SetMatrixInterpolationS3<1>),
           new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfesS,pfes>::Op*,E_F_StackF0F0>(SetMatrixInterpolationS2<1>),
           new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfes,pfesL>::Op*,E_F_StackF0F0>(SetMatrixInterpolation2L<1>),
           new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfes,pfesS>::Op*,E_F_StackF0F0>(SetMatrixInterpolation2S<1>),
           new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfes,pfes3>::Op*,E_F_StackF0F0>(SetMatrixInterpolation23<1>)
		   );


 TheOperators->Add("<-",
		   new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfes,pfes>::Op*,E_F_StackF0F0>(SetMatrixInterpolation<0>),
		   new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfes3,pfes3>::Op*,E_F_StackF0F0>(SetMatrixInterpolation3<0>),
           new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfesS,pfesS>::Op*,E_F_StackF0F0>(SetMatrixInterpolationS<0>),
		   new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfesL,pfesL>::Op*,E_F_StackF0F0>(SetMatrixInterpolationL<0>),
           new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfesL,pfesS>::Op*,E_F_StackF0F0>(SetMatrixInterpolationLS<0>),
           new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfesL,pfes>::Op*,E_F_StackF0F0>(SetMatrixInterpolationL2<0>),
           new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfesS,pfes3>::Op*,E_F_StackF0F0>(SetMatrixInterpolationS3<0>),
           new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfesS,pfes>::Op*,E_F_StackF0F0>(SetMatrixInterpolationS2<0>),
           new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfes,pfesL>::Op*,E_F_StackF0F0>(SetMatrixInterpolation2L<0>),
           new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfes,pfesS>::Op*,E_F_StackF0F0>(SetMatrixInterpolation2S<0>),
           new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation<pfes,pfes3>::Op*,E_F_StackF0F0>(SetMatrixInterpolation23<0>)
		   );
 // construction of complex matrix form a double matrix
 TheOperators->Add("=", new OneOperator2_<Matrice_Creuse<Complex>*,Matrice_Creuse<Complex>*,Matrice_Creuse<double>*,E_F_StackF0F0>(CopyMat<R,Complex,1>)
		   );

 TheOperators->Add("<-", new OneOperator2_<Matrice_Creuse<Complex>*,Matrice_Creuse<Complex>*,Matrice_Creuse<double>*,E_F_StackF0F0>(CopyMat<R,Complex,0>)
		   );
    Dcl_Type<Matrice_Creuse_C2R>();
    Add<Matrice_Creuse<Complex>*>("re",".",new OneOperator1s_<newpMatrice_Creuse<double> ,Matrice_Creuse<Complex>* >(Build_Matrice_Creuse_C2R<0> ));
    Add<Matrice_Creuse<Complex>*>("im",".",new OneOperator1s_<newpMatrice_Creuse<double>  ,Matrice_Creuse<Complex>* >(Build_Matrice_Creuse_C2R<1> ));
    // construction of complex matrix form a double matrix
 init_UMFPack_solver();
 init_HashMatrix ();

 Global.Add("renumbering", "(", new removeDOF<double>);
 Global.Add("renumbering", "(", new removeDOF<Complex>);
    Global.Add("renumbering", "(", new removeDOF<double>(1));
    Global.Add("renumbering", "(", new removeDOF<Complex>(1));
    Global.Add("renumbering", "(", new removeDOF<double>(1, 1));
    Global.Add("renumbering", "(", new removeDOF<Complex>(1, 1));
    Global.Add("renumbering", "(", new removeDOF<double>(1, 1, 1));
    Global.Add("renumbering", "(", new removeDOF<Complex>(1, 1, 1));

  Global.Add("display", "(", new plotMatrix<double>);
  Global.Add("display", "(", new plotMatrix<Complex>);

    // ZZZZZZ  ne marche pas FH....
    TheOperators->Add("*",
                     new Op2_mulvirtAvCR< RNM_VirtualMatrix<Complex>::plusAx,Matrice_Creuse<double>*,KN_<Complex> > ,
                     new Op2_mulvirtAvCR< RNM_VirtualMatrix<Complex>::plusAtx,Matrice_Creuse_Transpose<double>,KN_<Complex> > ,
                     new Op2_mulvirtAvCR< RNM_VirtualMatrix<Complex>::solveAxeqb,Matrice_Creuse_inv<R>,KN_<Complex> >,
                     new Op2_mulvirtAvCR< RNM_VirtualMatrix<Complex>::solveAtxeqb,Matrice_Creuse_inv<R>,KN_<Complex> >
                     );
     init_SparseLinearSolver();


    Global.New("DefaultSolver",CPValue<string*>(def_solver));
    Global.New("DefaultSolverSym",CPValue<string*>(def_solver_sym));
    Global.New("DefaultSolverSDP",CPValue<string*>(def_solver_sym_dp));

    Global.Add("removeHalf", "(", new OneOperator2s_<newpMatrice_Creuse<R> ,Matrice_Creuse<R> * ,long>(removeHalf));
    Global.Add("removeHalf", "(", new OneOperator3s_<newpMatrice_Creuse<R> ,Matrice_Creuse<R> * ,long,double>(removeHalf));
    Global.Add("removeHalf", "(", new OneOperator4s_<bool,Matrice_Creuse<R> * ,Matrice_Creuse<R> * ,long,double>(removeHalf));

}

int Data_Sparse_Solver_version() { return VDATASPARSESOLVER;}
#else
#include <iostream>

void  init_lgmat()
{  std::cout << "\n warning  init_lgmat EMPTY\n"<< std::endl;

}
#endif
