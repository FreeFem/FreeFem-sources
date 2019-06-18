/****************************************************************************/
/* This file is part of FreeFem++.                                          */
/*                                                                          */
/* FreeFem++ is free software: you can redistribute it and/or modify        */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFem++ is distributed in the hope that it will be useful,             */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFem++. If not, see <http://www.gnu.org/licenses/>.        */
/****************************************************************************/
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : F. Hecht
// E-MAIL  :  frederic.hecht@sorbonne-universite.fr


// *INDENT-OFF* //
//ff-c++-LIBRARY-dep: lapack blas
// *INDENT-ON* //

#include "ff++.hpp"
#include "AFunction.hpp"
#include "AFunction_ext.hpp"
#include <iostream>
#include <vector>

#ifdef __LP64__
typedef int intblas;
typedef int integer;
#else
typedef long intblas;
typedef long integer;
#endif

typedef integer logical;
typedef float LAPACK_real;
typedef double doublereal;
typedef logical (*L_fp)();
typedef integer ftnlen;
typedef complex<float> LAPACK_complex;
typedef complex<double> doublecomplex;
typedef void VOID;
#define complex LAPACK_complex
#define real LAPACK_real

#include "clapack.h"
#undef real
#undef complex
// KN<long>*CloseTo (Stack stack, double const &eps, KNM<double> *const &p, KNM<double> *const &q) {

template<class R>
long  ff_ShurComplement(Stack stack,KNM<R> *  pS,Matrice_Creuse<R> *  pmcA,KN_<long> const &  I,Data_Sparse_Solver &ds,KNM<R> *  pV=0   );

template<class R>
class ShurComplement_OP : public E_F0mps { public:
    Expression ess,eaa,eii,epV;
    
    static  aType btype;
    static const int n_name_param =NB_NAME_PARM_MAT; //  add nbiter FH 30/01/2007 11 -> 12  //add var MUMPS+autre
    static basicAC_F0::name_and_type name_param[] ;
    Expression nargs[n_name_param];
    const OneOperator * precon;
    
public:
    ShurComplement_OP(const basicAC_F0 &  args,Expression s,Expression a,Expression i,Expression ev=0) : ess(s),eaa(a),eii(i),epV(ev)  {
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
class ShurComplement : public OneOperator { public:  
    int cas;
    ShurComplement() : OneOperator(atype<long>(),atype<KNM<R> *>(),atype<Matrice_Creuse<R> *>(),atype<KN<long> *>() ) ,  cas(0) {}
    ShurComplement(int ) : OneOperator(atype<long>(),atype<KNM<R> *>(),atype<Matrice_Creuse<R> *>(),atype<KN<long> *, atype<KNM<R> *>()>() ) ,  cas(1){}

    E_F0 * code(const basicAC_F0 & args) const
    {
        if(cas==0)
        return  new ShurComplement_OP<R>(args,t[0]->CastTo(args[0]),t[1]->CastTo(args[1]),t[2]->CastTo(args[2]));
        else
        return  new ShurComplement_OP<R>(args,t[0]->CastTo(args[0]),t[1]->CastTo(args[1]),t[2]->CastTo(args[2]),t[3]->CastTo(args[3]));

    }
};

template <class R> 
basicAC_F0::name_and_type  ShurComplement_OP<R>::name_param[]= {
    LIST_NAME_PARM_MAT
    
};

template<class R>
AnyType ShurComplement_OP<R>::operator()(Stack stack)  const
{
    Matrice_Creuse<R> *  pA= GetAny<Matrice_Creuse<R> *>((*eaa)(stack));
    KNM<R> *  pE= GetAny<KNM<R> *>((*ess)(stack));
    KNM<R> *  pV=0;
    if(epV) pV=GetAny<KNM<R> *>((*epV)(stack));
    KN<long> *  pII= GetAny<KN<long> *>((*eii)(stack));
    Data_Sparse_Solver ds;
    ds.factorize=0;
    SetEnd_Data_Sparse_Solver<R>(stack,ds,nargs,n_name_param);
    return ff_ShurComplement(stack,pE,pA,*pII,ds,pV);
}
template<class R>
long  ff_ShurComplement(Stack stack,KNM<R> *  pS,Matrice_Creuse<R> *  pmcA,KN_<long> const &  I,Data_Sparse_Solver &ds,KNM<R> *  pV   )
{
    // I given numbering of Shur complement I[i] is the index in A of the i in S
    R zero(0.);
    KNM<R>   & S= *pS;
    MatriceCreuse<R> * pa=pmcA->A;
    ffassert( pa  );
    MatriceMorse<R> *pA= dynamic_cast<MatriceMorse<R>* > (pa);
    ffassert(pA);
    MatriceMorse<R> & A=*pA;
    int ni = I.N();
    int n=pA->n,m=pA->m;
    ffassert( n == m);
    S.resize(ni,ni);
    pV->resize(n,ni); //
    S=zero;
    ffassert( n>ni);
    KN<long> mark(n,-1L);
    for(int i=0; i< ni; ++i)
        mark[I[i]]=i;
    int nj=0;
    for(int i=0; i<n;++i)
        if (mark[i] <0 ) mark[i] = -2 - nj++;
    KN<int> J(nj);
    for(int i=0; i<n;++i)
    if (mark[i] <0 ) J[-mark[i]+2]= i;
    // I, J partionne in 2 set ..
    // the 4 matrix
    MatriceMorse<R> AII(ni,ni),AIJ(ni,nj), AJI(nj,ni), AJJ(nj,nj);
    int nstep = A.half;
    ffassert( nstep==0 || nstep ==1);
    for(int step = 0; step <= nstep ; ++step)
    for( int k=0; k< A.nnz;++k)
    {
        int i = A.i[k];
        int j = A.j[k];
        if( step )
        {
            if(i==j) continue;
            std::swap(i,j);
        }
        R aij = A.aij[k];
        int mi= mark[i];
        int mj = mark[j];
        int ki = mi <0 ? -mi-2 : -1;
        int kj = mj <0 ? -mj-2 : -1;
        if( mi>=0 )
        {
            if(  mj >=0 )//II
                AII(mi,mj) += aij;
            else
                AIJ(mi,kj) += aij;
        }
            else
                if(  mj >=0 )//JI
                    AJI(ki,mj) += aij;
                else// JJ
                    AJJ(ki,kj) += aij;
    }


    if(verbosity>99)
    {
    AII.COO();
    AJI.COO();
    AIJ.COO();
    AJJ.COO();

    cout << " AII "<< AII << endl;
    cout << " AIJ "<< AIJ << endl;
    cout << " AJI "<< AJI << endl;
    cout << " AJJ "<< AJJ << endl;
    }
    //  compute   AII - AIJ AJJ^-1 AJI
    // AII(:,k) = AII(:,k) - AIJ * AJJ^-1 AJI(:,k) *
    //  set solver .sur AJJ
    AJJ.setsdp(ds.sym,ds.positive); // put the matrix in rigth format
    bool VF=false;
    SetSolver<R>(stack,VF,AJJ,ds);
    // methode brutale a optimiser
    KN<R> rJ(nj),sJ(nj),sI(ni);
    AJI.CSC(); // column priority
    AII.CSC(); // column priority
    long err=0;
    for( int k=0; k<ni; ++k) // for each col
    {
        rJ=zero;
        
        for (int l = AJI.p[k]; l <AJI.p[k+1];++l)
            rJ[AJI.i[l]]= AIJ.aij[l], err+= AJI.j[l]!=k;
        AJJ.solve(sJ,rJ);
        if(pV)
        {
            for(int i=0; i< n; ++i)
            {
               int mi= mark[i];
              int ki = mi <0 ? -mi-2 : -1;
              if( mi < 0)
                  (*pV)(i,k) = sJ[ki];
              else
                  (*pV)(i,k) = R(k==mi);
            }
        }
        sI = AIJ*sJ;
        S(':',k) = -sI;
        for (int l = AII.p[k]; l <AII.p[k+1];++l)
            S(AII.i[l],k)+= AII.aij[l],   err+= AII.j[l]!=k;
    }
    ffassert(err==0);
    err=0;
    return ni;
}
template<class R>
long copy_mat(KNM<R> *  pS,Matrice_Creuse<R> *  pmcA)
{
    R zero(0.);
    KNM<R>   & S= *pS;
    MatriceCreuse<R> * pa=pmcA->A;
    ffassert( pa  );
    MatriceMorse<R> *pA= dynamic_cast<MatriceMorse<R>* > (pa);
    ffassert(pA);
    MatriceMorse<R> & A=*pA;
    int n = A.n, m= A.m;
    pS->resize(n,m);
    S = zero;
    for(long k= 0; k< A.nnz; ++k)
    {
        int i = A.i[k],j=A.j[k];
        R aij = A.aij[k];
        S(i,j) += aij;
        if( A.half && i != j) // if half take other part ..
            S(j,i)+= aij;
    }
   return 1;
}
static void Load_Init () {
    cout << " load: init ShurComplement " << endl;
    Global.Add("ShurComplement", "(", new ShurComplement<R>);
    Global.Add("ShurComplement", "(", new ShurComplement<Complex>);
    Global.Add("copy","(", new OneOperator2<long, KNM<R> *,Matrice_Creuse<R> *  >(copy_mat));
    Global.Add("copy","(", new OneOperator2<long, KNM<Complex> *,Matrice_Creuse<Complex> *  >(copy_mat));

    
}

LOADFUNC(Load_Init)
