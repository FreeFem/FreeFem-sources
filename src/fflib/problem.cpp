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
#include  <iostream>
using namespace std;

#include "rgraph.hpp"
#include "error.hpp"
#include "AFunction.hpp"

//#include "lex.hpp"
#include "HashMatrix.hpp"

#include "SparseLinearSolver.hpp"
#include "Mesh3dn.hpp"
#include "MeshPoint.hpp"
#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"
#include <set>



basicAC_F0::name_and_type  CDomainOfIntegration::name_param[]= {
    { "qft", &typeid(const Fem2D::QuadratureFormular *)},
    { "qfe", &typeid(const Fem2D::QuadratureFormular1d *)},
    { "qforder",&typeid(long)},
    { "qfnbpT",&typeid(long)},
    { "qfnbpE",&typeid(long)},
    { "optimize",&typeid(long)},
    { "binside",&typeid(double)},
    { "mortar",&typeid(bool)},
    { "qfV", &typeid(const Fem2D::GQuadratureFormular<R3> *)},
    { "levelset",&typeid(double)},
    { "mapt",&typeid(E_Array)},
    { "mapu",&typeid(E_Array)}


};

basicAC_F0::name_and_type  Problem::name_param[]= {
    {  "save",&typeid(string* )},
    {  "cadna",&typeid(KN<double>*)},
    {  "bmat",&typeid(Matrice_Creuse<R>* )},
    LIST_NAME_PARM_MAT
    /*
     {  "init", &typeid(bool)},
     {  "solver", &typeid(TypeSolveMat*)},
     {  "eps", &typeid(double) },
     {  "precon",&typeid(Polymorphic*)},
     {  "dimKrylov",&typeid(long)},
     {  "bmat",&typeid(Matrice_Creuse<R>* )},
     {  "tgv",&typeid(double )},
     {  "strategy",&typeid(long )},
     {  "save",&typeid(string* )},
     {  "cadna",&typeid(KN<double>*)},
     {  "tolpivot", &typeid(double)},
     {  "tolpivotsym", &typeid(double)},
     {  "nbiter", &typeid(long)}, // 12
     {   "paramint",&typeid(KN_<long>)}, // Add J. Morice 02/09
     {   "paramdouble",&typeid(KN_<double>)},
     {   "paramstring",&typeid(string *)},
     {   "permrow",&typeid(KN_<long>)},
     {   "permcol",&typeid(KN_<long>)},
     {   "fileparamint",&typeid(string*)}, // Add J. Morice 02/09
     {   "fileparamdouble",&typeid(string*)},
     {   "fileparamstring",&typeid(string* )},
     {   "filepermrow",&typeid(string*)},
     {   "filepermcol",&typeid(string*)} //22
     */
};

struct pair_stack_double
{
    Stack first;
    double *second;
    pair_stack_double(Stack ss,double* bb) : first(ss),second(bb) {};

};

namespace Fem2D {

    void  Expandsetoflab(Stack stack,const CDomainOfIntegration & di,set<int> & setoflab,bool &all);
    void  Expandsetoflab(Stack stack,const BC_set & bc,set<long> & setoflab);

    void Check(const Opera &Op,int N,int  M)
    {
        int err=0;
        for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++)
        {  // attention la fonction test donne la ligne
            //  et la fonction test est en second
            BilinearOperator::K ll(*l);
            pair<int,int> jj(ll.first.first),ii(ll.first.second);
            if (ii.first <0 || ii.first >= M) err++;
            if (jj.first <0 || jj.first >= N) err++;

        }
        if (err) {
            for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++)
            {  // attention la fonction test donne la ligne
                //  et la fonction test est en second
                BilinearOperator::K ll(*l);
                pair<int,int> jj(ll.first.first),ii(ll.first.second);
                cout << " +  " << jj.first << " " << jj.second << "*" << ii.first << " " << ii.second << endl;
                cout << "check ii.first ="  << ii.first << " < M " << M <<  endl;
                cout << "check jj.first ="  << jj.first << " < N " << N <<  endl;
            }
            cout << " value of N=" << N << endl;
            cout << " value of M=" << M << endl;
            ExecError("Check BilinearOperator N M");
        }
    }
    void Check(const  BC_set * bc,int N)
    {
        int err=0;
        int kk=bc->bc.size();
        for (int k=0;k<kk;k++)
        {
            pair<int,Expression> xx=bc->bc[k];
            if (xx.first >= N) {
                err++;
                cerr << " Sorry : just " << N << " componant in FE space \n"
                << "   and Boundary condition refere to " << xx.first << "+1 componant " << endl;
            }
        }
        if (err)
        ExecError("Incompatibility between  boundary condition  and FE space");
    }

    void Check(const  Ftest * fl,int N)
    {
        assert(fl);
        int err=0;
        Ftest::const_iterator kk= fl->v.end(),k;
        int ii=0;
        for (k=fl->v.begin();k<kk;k++)
        {
            ii++;
            int j=k->first.first;
            if (  j >= N) {
                err++;
                cerr << " Sorry : just " << N << " componant in FE space \n"
                << " and linear var form  refere to " << j << "+1 componant (part " << ii << ")" << endl;
            }
        }
        if (err)
        ExecError("Incompatibility between linear varf  and FE space");
    }
  template<class R>
  inline void  CheckErrorOptimisation(const R& ccc,const R& cc,const char * cmm)
    {
        if ( ccc != cc) {
             if( (abs(ccc-cc) >1e-8*(abs(cc)+abs(ccc)) ) // test for round off err
                 || (cc !=cc) ||  (ccc !=ccc) ) {// test for NaN
                cerr << cc << " != " << ccc <<  " diff "<< cc-ccc <<" => ";
                cerr << cmm << endl;
                 cerr << " remark if you add  (..  ,   optimize=2) then  you remove this check (be careful); "<< endl;
                ExecError("In Optimized version "); }}
   }
    inline void  CheckErrorOptimisation(const Complex ccc,const Complex& cc,const char * cmm)
    {
        if ( ccc != cc) {
            if( (abs(ccc.real()-ccc.real())+ abs(cc.imag()-cc.imag())  >0.5e-8*( abs(cc.real())+abs(cc.imag())+abs(ccc.real())+abs(ccc.imag()) ) )
                ||  (cc !=cc) ||  (ccc !=ccc) )
            {
                cerr << cc << " != " << ccc <<  " diff "<< cc-ccc <<" => ";
                cerr << cmm << endl;
                cerr << " remark if you add  (..  ,   optimize=2) then  you remove this check (be careful); " <<endl;

                ExecError("In Optimized version "); }}
    }

    //---------------------------------------------------------------------------------------
/*  just fo debug ...
  long bbnElementonB(Stack s)  {
    throwassert(*((long *)s));
    MeshPoint *mp = MeshPointStack(s);
    long l = 0;
    if ((mp->T) && (mp->e > -1) && (mp->d == 2))
      l = mp->Th->nTonEdge(mp->t, mp->e);
    else if (mp->d == 3 && mp->dHat == 3 && mp->T3 && (mp->f >= 0))
    {  cout << mp->t <<" .. " << mp->f << " ";
      l = mp->Th3->nElementonB(mp->t, mp->f);
    }
    else if (mp->d == 3 && mp->dHat == 2 && mp->TS && (mp->e >= 0))
      l = mp->ThS->nElementonB(mp->t, mp->e);
    return l;
  }
*/
template<class R>
    void  Element_OpVF(MatriceElementairePleine<R,FESpace3> & mat,
                       const FElement3 & Ku,const FElement3 & KKu,
                       const FElement3 & Kv,const FElement3 & KKv,
                       double * p,int ie,int iie, int label,void *bstack,R3 *B)
    {
        typedef typename FElement3::Element Element;
        ffassert(B==0);
        pair_stack_double * bs=static_cast<pair_stack_double *>(bstack);
        Stack stack= bs->first;
        double binside = *bs->second; // truc FH pour fluide de grad2 (decentrage bizard)
        ffassert(mat.onFace); //   Finite Volume or discontinuous Galerkine
        ffassert(ie>=0 && ie < 4); //  int on Face
        MeshPoint mp= *MeshPointStack(stack);
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);

        bool same = &Ku == & Kv;
        const Element & T  = Ku.T;
        const Element & TT  = KKu.T;
        bool sameT = &Ku == & KKu;
        int nTonEdge =  &Ku == &KKu ? 1 : 2;
        double cmean = 1./nTonEdge;

        throwassert(&T == &Kv.T);
        // const QuadratureFormular & FI = mat.FIT;
        const QuadratureFormular & FIb = mat.FIE;
        long npi;
        R *a=mat.a;
        R *pa=a;
        long i,j;
        long n= mat.n,m=mat.m,nx=n*m;
        assert(nx<=mat.lga);
        long N= Kv.N;
        long M= Ku.N;

        long mu=Ku.NbDoF();
        long mmu=KKu.NbDoF();
        long nv=Kv.NbDoF();
        long nnv=Kv.NbDoF();
        assert(mu==mmu && nv == nnv) ;



        const Opera &Op(*mat.bilinearform);
        bool classoptm = copt && Op.optiexpK;
        //  if (Ku.number<1 && verbosity/100 && verbosity % 10 == 2)
        if (Ku.number<1 && ( verbosity > 1 ) )
        cout << "Element_OpVF 3d P: copt = " << copt << " " << classoptm << " binside (For FH) =" << binside << " opt: " << mat.optim << endl;

         bool oldopt=1;  // juin 2007 FH ???? a voir
        int  iloop=0;
        KN<bool> unvarexp(classoptm ? Op.optiexpK->sizevar() : 1);
        if (Ku.number<1 && verbosity/100 && verbosity % 10 == 2)
        cout << "Element_Op 3d P: copt = " << copt << " " << classoptm << " opt: " << mat.optim << endl;
         //
        int lastop=0;
        lastop = 0;
        What_d Dop = Op.DiffOp(lastop);
  
        long lffv = sameT ? 0  : nv*N*lastop;
        long lffu = sameT ? 0  : mu*M*lastop;
        long loffset =  same ? 0 :  (nv+nnv)*N*lastop;

        RNMK_ fv(p,nv,N,lastop); //  the value for basic fonction in K
        RNMK_ ffv(p + lffv ,nnv,N,lastop); //  the value for basic fonction in KK
        RNMK_ fu(  (double*) fv   + loffset  ,mu,M,lastop); //  the value for basic fonction
        RNMK_ ffu( (double*) fu  + lffu  ,mmu,M,lastop); //  the value for basic fonction
        
        


        R3 NN= T.N(ie);
        double mes=NN.norme();
        NN/=mes;
 
        // compute the permutaion face ie to iie ...
        // a little tricky
        int fpe= T.facePermutation(ie);
        int fpee= TT.facePermutation(iie);
        int pr[3],ppr[3];
        SetNumPerm<3>(fpe,pr);
        SetNumPerm1<3>(fpee,ppr);
 
        for (npi=0;npi<FIb.n;npi++) // loop on the integration point
        {
            pa =a;
            GQuadraturePoint<R2> pi( FIb[npi]);
            double lpi[3];
            pi.toBary(lpi);
            double llpi[]={lpi[pr[ppr[0]]], lpi[pr[ppr[1]]], lpi[pr[ppr[2]]]};
            double coef = 0.5*mes*pi.a; // correction 0.5 050109 FH
            R3 Pt(T.PBord(ie,pi));
            R2 pii(llpi+1);
            R3 PP_t(TT.PBord(iie,pii));
            {
                static int err=0;
                R3 P(T(Pt)),PP(TT(PP_t)),D(P,PP);
               if(  D.norme2() > 1e-5)
               {
                   const Mesh3 & Th=Ku.Vh.Th;
                   cout << " pr "<< pr[0] << pr[1]<< pr[2] << " prr " << ppr[0] << ppr[1]<< ppr[2]<<" "
                   << ppr[pr[0]] <<  ppr[pr[1]] << ppr[pr[2]]  << "   " <<  pr[ppr[0]] <<  pr[ppr[1]] << pr[ppr[2]]<< endl;
                   cout << ie << " " << iie << endl;
                   cout << " T = " << Th(T[0]) << " " <<Th(T[1]) << " " << Th(T[2]) << " " << Th(T[3]) << endl;
                   cout << " TT = " << Th(TT[0]) << " " <<Th(TT[1]) << " " << Th(TT[2]) << " " << Th(TT[3]) << endl;
                   cout << fpe << " " << fpee << endl;
                   cout << P << " " << PP << " diff=D.norme2() " << D.norme2()  << endl;
                   err++;
                   ffassert(err<10);
               }
                
            }
            Ku.BF(Dop,Pt,fu);
            if(!sameT) KKu.BF(Dop,PP_t,ffu);
            if (!same) {
                Kv.BF(Dop,Pt,fv);
                if(!sameT) KKv.BF(Dop,PP_t,ffv);
            }
            if(verbosity==99999 &&( Ku.number == 3 || KKu.number == 2)) {
                cout << "intallfaces " << Ku.number << " " << fu << endl;
                cout << "intallfaces " << KKu.number << " " << ffu << endl;
            }
            MeshPointStack(stack)->set(T(Pt),Pt,Ku,label,NN,ie);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version

 
            
 


            for ( i=0;  i<n;   i++ )
            {
                int ik= mat.nik[i];
                int ikk=mat.nikk[i];

                RNM_ wi(fv(Max(ik,0),'.','.'));
                RNM_ wwi(ffv(Max(ikk,0),'.','.'));

                for ( j=0;  j<m;   j++,pa++ )
                {
                    int jk= mat.njk[j];
                    int jkk=mat.njkk[j];

                    RNM_ wj(fu(Max(jk,0),'.','.'));
                    RNM_ wwj(ffu(Max(jkk,0),'.','.'));

                    int il=0;
                    for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                    {
                        BilinearOperator::K ll(*l);
                        pair<int,int> jj(ll.first.first),ii(ll.first.second);
                        int iis = ii.second, jjs=jj.second;

                        int iicase  = iis / last_operatortype;
                        int jjcase  = jjs / last_operatortype;

                        iis %= last_operatortype;
                        jjs %= last_operatortype;
                        double w_i=0,w_j=0,ww_i=0,ww_j=0;

                        if(ik>=0) w_i =   wi(ii.first,iis );
                        if(jk>=0) w_j =   wj(jj.first,jjs );

                        if( iicase>0 && ikk>=0) ww_i =  wwi(ii.first,iis );
                        if( jjcase>0 && jkk>=0) ww_j =  wwj(jj.first,jjs );


                        if       (iicase==Code_Jump) w_i = ww_i-w_i; // jump
                        else  if (iicase==Code_Mean) {

                            w_i = cmean*  (w_i + ww_i );} // average
                        else  if (iicase==Code_OtherSide) w_i = ww_i;  // valeur de autre cote

                        if      (jjcase==Code_Jump) w_j = ww_j-w_j; // jump
                        else if (jjcase==Code_Mean) w_j = cmean*  (w_j +ww_j ); // average
                        else if (jjcase==Code_OtherSide) w_j = ww_j;  //  valeur de l'autre cote

                        // R ccc = GetAny<R>(ll.second.eval(stack));

                        R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                        if ( copt && ( mat.optim==1) && Kv.number <1)
                        {
                            R cc  =  GetAny<R>(ll.second.eval(stack));
                            CheckErrorOptimisation(cc,ccc,"Sorry error in Optimization Element_OpVF3d  (face) add:   intallface(Th,optimize=0)(...)");
                         }
                        if(verbosity==99999)
                            cout << ie << " / " << i << " " << j << "/  (" << mat.ni[i] << " "<< mat.nj[j] << ")  pi " << npi << " : " << coef << " c= " << ccc << " " << w_i << " " << w_j << " = " << coef * ccc * w_i*w_j
                            << " k/kk  " << Ku.number << " "<< KKu.number << " :: " << ik << " " << ikk  << ",  "<< jk << " " << jkk  << " ii/jj2= " << ii.second << " " << jj.second  // << " // " << bbnElementonB(stack)
                            << endl;
                        *pa += coef * ccc * w_i*w_j;
                    }
                }
            }
            // else pa += m;
        }


        pa=a;
        if ( (verbosity > 9999) ||( (verbosity > 55) && (Ku.number <=0 || KKu.number <=0 )))  {
            cout <<endl  << " face between " << Ku.number << " , " <<  KKu.number   << " =  "<<  T[0] << ", " << T[1] << ", " << T[2] << " " << nx << endl;
            cout << " K u, uu =  " << Ku.number << " " << KKu.number << " " <<  " K v, vv =  " << Kv.number << " " << KKv.number << " " <<endl;
            for (int i=0;i<n;i++)
            {
                cout << setw(2) << i << setw(4) << mat.ni[i] <<  setw(4) << mat.nik[i] << setw(4) << mat.nikk[i]  <<  " :";
                for (int j=0;j<m;j++)
                cout << setw(5)  << (*pa++) << " ";
                cout << endl;
            } }

        *MeshPointStack(stack) = mp;
        
    }

    template<class R>
    void  Element_OpVF(MatriceElementairePleine<R,FESpaceS> & mat,
                       const FElementS & Ku,const FElementS & KKu,
                       const FElementS & Kv,const FElementS & KKv,
                       double * p,int ie,int iie, int label,void *bstack,R3 *B)
    {
        ffassert(0);
    }

    template<class R>
    void  Element_OpVF(MatriceElementairePleine<R,FESpaceL> & mat,
                       const FElementL & Ku,const FElementL & KKu,
                       const FElementL & Kv,const FElementL & KKv,
                       double * p,int ie,int iie, int label,void *bstack,R3 *B)
    {
        typedef typename FElementL::E Element;
        ffassert(B==0);
        pair_stack_double * bs=static_cast<pair_stack_double *>(bstack);
        Stack stack= bs->first;
        double binside = *bs->second; // truc FH pour fluide de grad2 (decentrage bizard)
        assert(mat.onFace); //   Finite Volume or discontinuous Galerkine
        assert(ie>=0 && ie < 3); //  int on edge
        MeshPoint mp= *MeshPointStack(stack);
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);

        bool same = &Ku == & Kv;
        assert(same);
        const Element & T  = Ku.T;
        R3 NNt=T.TangenteUnitaire();
        double s = ie;
        double ss = iie;
        GQuadraturePoint<R1> pi(1.,s);
        R mes = 1.;
        R3 NN=NNt;
        if(ie==0) NN=-NN;
        R coef = 1.;
         R1 Pt(pi);
        int nTonEdge =  &Ku == &KKu ? 1 : 2;
        double cmean = 1./nTonEdge;

        throwassert(&T == &Kv.T);
        long npi=1;
        R *a=mat.a;
        R *pa=a;
        long i,j;
        long n= mat.n,m=mat.m,nx=n*m;
        assert(nx<=mat.lga);
        long N= Kv.N;
        long M= Ku.N;

        long mu=Ku.NbDoF();
        long mmu=KKu.NbDoF();
        long nv=Kv.NbDoF();
        long nnv=Kv.NbDoF();
        assert(mu==mmu && nv == nnv) ;



        const Opera &Op(*mat.bilinearform);
        bool classoptm = copt && Op.optiexpK;
        //  if (Ku.number<1 && verbosity/100 && verbosity % 10 == 2)
        if (Ku.number<1 && ( verbosity > 1 ) )
        cout << "Element_OpVF 0d P: copt = " << copt << " " << classoptm << " binside (For FH) =" << binside << " opt: " << mat.optim << endl;


        int lastop;
        lastop = 0;
        What_d Dop = Op.DiffOp(lastop);
        //assert(lastop<=3);
        int lffv = nv*N*last_operatortype;
        int lffu = mu*M*last_operatortype;
        int loffset =  same ? 0 :  (nv+nnv)*N*last_operatortype;

        RNMK_ fv(p,nv,N,lastop); //  the value for basic fonction in K
        RNMK_ ffv(p + lffv ,nnv,N,lastop); //  the value for basic fonction in KK
        RNMK_ fu(  (double*) fv   + loffset  ,mu,M,lastop); //  the value for basic fonction
        RNMK_ ffu( (double*) fu  + lffu  ,mmu,M,lastop); //  the value for basic fonction


       
        {
            pa =a;
            double coef = 1.;
            
            R1 Pt=s; //
            R1 PP_t=ss;
            assert (!binside);
            Ku.BF(Dop,Pt,fu);
            KKu.BF(Dop,PP_t,ffu);
            if (!same) { Kv.BF(Dop,Pt,fv); KKv.BF(Dop,PP_t,ffv); }
            // int label=-999999; // a passer en argument
            MeshPointStack(stack)->set(T(Pt),Pt,Kv,label,NN,NNt,ie);//   Axel
           // MeshPointStack(stack)->set(T(Pt),Pt,Kv,label, Normal,ie);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version


            for ( i=0;  i<n;   i++ )
            {
                int ik= mat.nik[i];
                int ikk=mat.nikk[i];

                RNM_ wi(fv(Max(ik,0),'.','.'));
                RNM_ wwi(ffv(Max(ikk,0),'.','.'));

                for ( j=0;  j<m;   j++,pa++ )
                {
                    int jk= mat.njk[j];
                    int jkk=mat.njkk[j];

                    RNM_ wj(fu(Max(jk,0),'.','.'));
                    RNM_ wwj(ffu(Max(jkk,0),'.','.'));

                    int il=0;
                    for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                    {
                        BilinearOperator::K ll(*l);
                        pair<int,int> jj(ll.first.first),ii(ll.first.second);
                        int iis = ii.second, jjs=jj.second;

                        int iicase  = iis / last_operatortype;
                        int jjcase  = jjs / last_operatortype;

                        iis %= last_operatortype;
                        jjs %= last_operatortype;
                        double w_i=0,w_j=0,ww_i=0,ww_j=0;

                        if(ik>=0) w_i =   wi(ii.first,iis );
                        if(jk>=0) w_j =   wj(jj.first,jjs );

                        if( iicase>0 && ikk>=0) ww_i =  wwi(ii.first,iis );
                        if( jjcase>0 && jkk>=0) ww_j =  wwj(jj.first,jjs );


                        if       (iicase==Code_Jump) w_i = ww_i-w_i; // jump
                        else  if (iicase==Code_Mean) {

                            w_i = cmean*  (w_i + ww_i );} // average
                        else  if (iicase==Code_OtherSide) w_i = ww_i;  // valeur de autre cote

                        if      (jjcase==Code_Jump) w_j = ww_j-w_j; // jump
                        else if (jjcase==Code_Mean) w_j = cmean*  (w_j +ww_j ); // average
                        else if (jjcase==Code_OtherSide) w_j = ww_j;  //  valeur de l'autre cote

                        // R ccc = GetAny<R>(ll.second.eval(stack));

                        R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                        if ( copt && ( mat.optim==1) && Kv.number <1)
                        {
                            R cc  =  GetAny<R>(ll.second.eval(stack));
                            CheckErrorOptimisation(cc,ccc,"Sorry error in Optimization Element_OpVF0d  (b) add:   int2d(Th,optimize=0)(...)");
                           }
                        *pa += coef * ccc * w_i*w_j;
                    }
                }
            }
            // else pa += m;
        }


        pa=a;
        if ( (verbosity > 9999) ||( (verbosity > 55) && (Ku.number <=0 || KKu.number <=0 )))  {
            cout <<endl  << " vertex between " << Ku.number << " , " <<  KKu.number   << " =  "<<  T[0] << ", " << T[1] << ", " << T[2] << " " << nx << endl;
            cout << " K u, uu =  " << Ku.number << " " << KKu.number << " " <<  " K v, vv =  " << Kv.number << " " << KKv.number << " " <<endl;
            for (int i=0;i<n;i++)
            {
                cout << setw(2) << i << setw(4) << mat.ni[i] <<  setw(4) << mat.nik[i] << setw(4) << mat.nikk[i]  <<  " :";
                for (int j=0;j<m;j++)
                cout << setw(5)  << (*pa++) << " ";
                cout << endl;
            } }

        *MeshPointStack(stack) = mp;
        
        
        
    }

    template<class R>
    void  Element_OpVF(MatriceElementairePleine<R,FESpace> & mat,
                       const FElement & Ku,const FElement & KKu,
                       const FElement & Kv,const FElement & KKv,
                       double * p,int ie,int iie, int label,void *bstack,R2 *B)
    {
        ffassert(B==0);
        pair_stack_double * bs=static_cast<pair_stack_double *>(bstack);
        Stack stack= bs->first;
        double binside = *bs->second; // truc FH pour fluide de grad2 (decentrage bizard)
        assert(mat.onFace); //   Finite Volume or discontinous Galerkine
        assert(ie>=0 && ie < 3); //  int on edge
        MeshPoint mp= *MeshPointStack(stack);
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);

        bool same = &Ku == & Kv;
        assert(same);
        const Triangle & T  = Ku.T;
        int nTonEdge =  &Ku == &KKu ? 1 : 2;
        double cmean = 1./nTonEdge;

        throwassert(&T == &Kv.T);
        // const QuadratureFormular & FI = mat.FIT;
        const QuadratureFormular1d & FIb = mat.FIE;
        long npi;
        R *a=mat.a;
        R *pa=a;
        long i,j;
        long n= mat.n,m=mat.m,nx=n*m;
        assert(nx<=mat.lga);
        long N= Kv.N;
        long M= Ku.N;

        long mu=Ku.NbDoF();
        long mmu=KKu.NbDoF();
        long nv=Kv.NbDoF();
        long nnv=Kv.NbDoF();
        assert(mu==mmu && nv == nnv) ;



        const Opera &Op(*mat.bilinearform);
        bool classoptm = copt && Op.optiexpK;
        //  if (Ku.number<1 && verbosity/100 && verbosity % 10 == 2)
        if (Ku.number<1 && ( verbosity > 1 ) )
        cout << "Element_OpVF P: copt = " << copt << " " << classoptm << " binside (For FH) =" << binside << " opt: " << mat.optim << endl;


        KN<bool> Dop(last_operatortype); //  sinon ca plate bizarre
        Op.DiffOp(Dop);
        int lastop=1+Dop.last([](bool x){return x;});
        //assert(lastop<=3);
        int lffv = nv*N*last_operatortype;
        int lffu = mu*M*last_operatortype;
        int loffset =  same ? 0 :  (nv+nnv)*N*last_operatortype;

        RNMK_ fv(p,nv,N,lastop); //  the value for basic fonction in K
        RNMK_ ffv(p + lffv ,nnv,N,lastop); //  the value for basic fonction in KK
        RNMK_ fu(  (double*) fv   + loffset  ,mu,M,lastop); //  the value for basic fonction
        RNMK_ ffu( (double*) fu  + lffu  ,mmu,M,lastop); //  the value for basic fonction

        R2 E=T.Edge(ie);
        double le = sqrt((E,E));
        R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
        PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]),
        PC(TriangleHat[OppositeVertex[ie]]);
        // warning the to edge are in opposite sens
        R2 PP_A(TriangleHat[VerticesOfTriangularEdge[iie][1]]),
        PP_B(TriangleHat[VerticesOfTriangularEdge[iie][0]]),
        PP_C(TriangleHat[OppositeVertex[ie]]);
        R2 Normal(E.perp()/-le);
        for (npi=0;npi<FIb.n;npi++) // loop on the integration point
        {
            pa =a;
            QuadratureFormular1dPoint pi( FIb[npi]);
            double coef = le*pi.a;
            double sa=pi.x,sb=1-sa;
            R2 Pt(PA*sa+PB*sb ); //
            R2 PP_t(PP_A*sa+PP_B*sb ); //
            if (binside) {
                Pt   = (1-binside)*Pt + binside*PC;
                PP_t  = (1-binside)*PP_t + binside*PP_C; }
            Ku.BF(Dop,Pt,fu);
            KKu.BF(Dop,PP_t,ffu);
            if (!same) { Kv.BF(Dop,Pt,fv); KKv.BF(Dop,PP_t,ffv); }
            // int label=-999999; // a passer en argument
            MeshPointStack(stack)->set(T(Pt),Pt,Kv,label, Normal,ie);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version


            for ( i=0;  i<n;   i++ )
            {
                int ik= mat.nik[i];
                int ikk=mat.nikk[i];

                RNM_ wi(fv(Max(ik,0),'.','.'));
                RNM_ wwi(ffv(Max(ikk,0),'.','.'));

                for ( j=0;  j<m;   j++,pa++ )
                {
                    int jk= mat.njk[j];
                    int jkk=mat.njkk[j];

                    RNM_ wj(fu(Max(jk,0),'.','.'));
                    RNM_ wwj(ffu(Max(jkk,0),'.','.'));

                    int il=0;
                    for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                    {
                        BilinearOperator::K ll(*l);
                        pair<int,int> jj(ll.first.first),ii(ll.first.second);
                        int iis = ii.second, jjs=jj.second;

                        int iicase  = iis / last_operatortype;
                        int jjcase  = jjs / last_operatortype;

                        iis %= last_operatortype;
                        jjs %= last_operatortype;
                        double w_i=0,w_j=0,ww_i=0,ww_j=0;

                        if(ik>=0) w_i =   wi(ii.first,iis );
                        if(jk>=0) w_j =   wj(jj.first,jjs );

                        if( iicase>0 && ikk>=0) ww_i =  wwi(ii.first,iis );
                        if( jjcase>0 && jkk>=0) ww_j =  wwj(jj.first,jjs );


                        if       (iicase==Code_Jump) w_i = ww_i-w_i; // jump
                        else  if (iicase==Code_Mean) {

                            w_i = cmean*  (w_i + ww_i );} // average
                        else  if (iicase==Code_OtherSide) w_i = ww_i;  // valeur de autre cote

                        if      (jjcase==Code_Jump) w_j = ww_j-w_j; // jump
                        else if (jjcase==Code_Mean) w_j = cmean*  (w_j +ww_j ); // average
                        else if (jjcase==Code_OtherSide) w_j = ww_j;  //  valeur de l'autre cote

                        // R ccc = GetAny<R>(ll.second.eval(stack));

                        R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                        if ( copt && ( mat.optim==1) && Kv.number <1)
                        {
                            R cc  =  GetAny<R>(ll.second.eval(stack));
                            CheckErrorOptimisation(cc,ccc,"Sorry error in Optimization Element_OpVF2d  (b) add:   int2d(Th,optimize=0)(...)");
                           /* if ( ccc != cc) {
                                cerr << cc << " != " << ccc << " => ";
                                cerr << "Sorry error in Optimization Element_OpVF2d  (b) add:  int2d(Th,optimize=0)(...)" << endl;
                                ExecError("In Optimized version "); }*/
                        }
                        *pa += coef * ccc * w_i*w_j;
                    }
                }
            }
            // else pa += m;
        }


        pa=a;
        if ( (verbosity > 9999) ||( (verbosity > 55) && (Ku.number <=0 || KKu.number <=0 )))  {
            cout <<endl  << " edge between " << Ku.number << " , " <<  KKu.number   << " =  "<<  T[0] << ", " << T[1] << ", " << T[2] << " " << nx << endl;
            cout << " K u, uu =  " << Ku.number << " " << KKu.number << " " <<  " K v, vv =  " << Kv.number << " " << KKv.number << " " <<endl;
            for (int i=0;i<n;i++)
            {
                cout << setw(2) << i << setw(4) << mat.ni[i] <<  setw(4) << mat.nik[i] << setw(4) << mat.nikk[i]  <<  " :";
                for (int j=0;j<m;j++)
                cout << setw(5)  << (*pa++) << " ";
                cout << endl;
            } }

        *MeshPointStack(stack) = mp;
    }

    //--------------------------------------------------------------------------------------




    // creating an instance of AssembleBilinearForm with MatriceCreuse
    // case 2d
    // --------- FH 120105
    template<class R>
    void AssembleBilinearForm(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                              MatriceCreuse<R>  & A, const  FormBilinear * b  )

    {
        /*FH:  case ..in 2D
         in varf ...
         standard case ..
         */
        StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
        bool sptrclean=true;
        const CDomainOfIntegration & di= *b->di;
        const Mesh * pThdi = GetAny<pmesh>( (* di.Th)(stack));
        if ( pThdi != &Th || &Uh.Th !=&Th || &Vh.Th !=&Th) {
            cout << " --Use matrix formulation ---" << endl;
            ExecError("No way to compute bilinear form with integrale of on mesh \n"
                      "  test  or unknown function  defined on another mesh! sorry to hard.   ");
        }
        SHOWVERB(cout << " FormBilinear " << endl);
        double CPU0 = CPUtime();
        MatriceElementaireSymetrique<R,FESpace> *mates =0;
        MatriceElementairePleine<R,FESpace> *matep =0;
        const int useopt=di.UseOpt(stack);
        double binside=di.binside(stack);

        //const vector<Expression>  & what(di.what);
        CDomainOfIntegration::typeofkind  kind = di.kind;
        set<int> setoflab;
        bool all=true;

        const Mesh & ThI = Th;// * GetAny<pmesh>( (* di.Th)(stack));
        bool sameMesh = &ThI == &Vh.Th &&  &ThI == &Uh.Th;

        //    const QuadratureFormular1d & FIE = di.FIE(stack);
        //    const QuadratureFormular & FIT = di.FIT(stack);
        const QuadratureFormular1d & FIEo = di.FIE(stack);
        const QuadratureFormular & FITo = di.FIT(stack);
        // const GQuadratureFormular<R3> & FIVo = di.FIV(stack);
        //  to change the quadrature on element ... may 2014 FH ..
        QuadratureFormular1d  FIE(FIEo,3);
        QuadratureFormular FIT(FITo,3);
        // GQuadratureFormular<R3>  FIV(FIVo,3);

        bool VF=b->VF();  // finite Volume or discontinuous Galerkin
        if (verbosity>2) cout << "  -- discontinuous Galerkin  =" << VF << " size of Mat =" << A.size()<< " Bytes\n";
        if (verbosity>3)
        {
            if (CDomainOfIntegration::int1d==kind) cout << "  -- boundary int border ( nQP: "<< FIE.n << ") ,"  ;
            else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all edges ( nQP: "<< FIE.n << "),"  ;
            else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF edges nQP: ("<< FIE.n << ")," ;
            else cout << "  --  int    (nQP: "<< FIT.n << " ) in "  ;
        }
        //if(di.islevelset()) InternalError("So no levelset integration type on this case (6)");
        if( di.withmap()) { ExecError(" no map  in the case (2)??");}
        if(di.islevelset() && ( (CDomainOfIntegration::int1d!=kind) && (CDomainOfIntegration::int2d!=kind) )  )
        InternalError("So no levelset integration type on no int1d case (6)");

        Expandsetoflab(stack,di, setoflab,all);
        /*
         for (size_t i=0;i<what.size();i++)
         {
         long  lab  = GetAny<long>( (*what[i])(stack));
         setoflab.insert(lab);
         if ( verbosity>3) cout << lab << " ";
         all=false;
         }*/
        if (verbosity>3) cout <<" Optimized = "<< useopt << ", ";
        const E_F0 * poptiexp0=b->b->optiexp0;

        int n_where_in_stack_opt=b->b->where_in_stack_opt.size();
        R** where_in_stack =0;
        if (n_where_in_stack_opt && useopt)
        where_in_stack = new R * [n_where_in_stack_opt];
        if (where_in_stack)
        {
            assert(b->b->v.size()==(size_t) n_where_in_stack_opt);
            for (int i=0;i<n_where_in_stack_opt;i++)
            {
                int offset=b->b->where_in_stack_opt[i];
                assert(offset>10);
                where_in_stack[i]= static_cast<R *>(static_cast<void *>((char*)stack+offset));
                *(where_in_stack[i])=0;
            }


            if(poptiexp0)
            (*poptiexp0)(stack);
            KN<bool> ok(b->b->v.size());
            {  //   remove the zero coef in the liste
                // R zero=R();
                int il=0;
                for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
                ok[il] =  ! (b->b->mesh_indep_stack_opt[il] && ( std::norm(*(where_in_stack[il])) < 1e-100 ) );
            }
            BilinearOperator b_nozer(*b->b,ok);
            if (verbosity % 10 > 3 )
            cout << "   -- nb term in bilinear form  (!0) : " << b_nozer.v.size()
            << "  total " << n_where_in_stack_opt << endl;

            if ( (verbosity/100) % 10 >= 2)
            {
                int il=0;

                for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
                cout << il << " coef (" << l->first << ") = " << *(where_in_stack[il])
                << " offset=" << b->b->where_in_stack_opt[il]
                << " dep mesh " << l->second.MeshIndependent() << b->b->mesh_indep_stack_opt[il] << endl;
            }
        }
        Stack_Ptr<R*>(stack,ElemMatPtrOffset) =where_in_stack;
        void *paramate=stack;
        pair_stack_double parammatElement_OpVF(stack,& binside);
        // parammatElement_OpVF.first = stack;
        // parammatElement_OpVF.second= & binside;

        if (verbosity >3)
        {
            if (all) cout << " all " << endl ;
            else cout << endl;
        }
        if(VF) {
            if(&Uh != &Vh || sym)
            cout << ("To Day in bilinear form with discontinuous Galerkin (2d):   \n"
                     "  test or unknown function must be  defined on the same FEspace, \n"
                     "  and the matrix is not symmetric. \n"
                     " To do other case in a future (F. Hecht) dec. 2003 ");
            if(&Uh == &Vh)
            matep= new MatriceElementairePleine<R,FESpace>(Uh,VF,FIT,FIE,useopt);
            else
            matep= new MatriceElementairePleine<R,FESpace>(Uh,Vh,VF,FIT,FIE,useopt);



            //      matep= new MatriceElementairePleine<R,FESpace>(Uh,Vh,VF,FIT,FIE);
            matep->faceelement = Element_OpVF;
            paramate= &parammatElement_OpVF;
        }
        else if (sym) {
            mates= new MatriceElementaireSymetrique<R,FESpace>(Uh,FIT,FIE,useopt);
            mates->element = Element_Op<R>;
        }
        else {
            matep= new MatriceElementairePleine<R,FESpace>(Uh,Vh,FIT,FIE,useopt);
            matep->element = Element_Op<R>;
        }
        MatriceElementaireFES<R,FESpace> & mate(*( sym? (MatriceElementaireFES<R,FESpace> *)mates : (MatriceElementaireFES<R,FESpace> *) matep));


        mate.bilinearform=b->b;

        Check(*mate.bilinearform,mate.Uh.N,mate.Vh.N);
        if(verbosity>9) cout << "  -- CPU init assemble mat " <<  CPUtime()-CPU0 << " s\n";

        if (di.kind == CDomainOfIntegration::int1d )
        {
            if(di.islevelset())
            {
                double uset = HUGE_VAL;
                R2 Q[3];
                KN<double> phi(Th.nv);phi=uset;
                double f[3];
                for(int t=0; t< Th.nt;++t)
                {
                    if ( all || setoflab.find(Th[t].lab) != setoflab.end())
                    {
                        double umx=-HUGE_VAL,umn=HUGE_VAL;
                        for(int i=0;i<3;++i)
                        {
                            int j= ThI(t,i);
                            if( phi[j]==uset)
                            {
                                MeshPointStack(stack)->setP(&ThI,t,i);
                                phi[j]= di.levelset(stack);//zzzz
                            }
                            f[i]=phi[j];
                            umx = std::max(umx,phi[j]);
                            umn = std::min(umn,phi[j]);

                        }
                        if( umn <=0 && umx >= 0)
                        {

                            int np= IsoLineK(f,Q,1e-10);
                            if(np==2)
                            {
                                /* if ( sameMesh)
                                 {

                                 Element_rhs<R>(Vh[t],*l->l,buf,stack,*B,FIE,Q[0],Q[1]);
                                 }
                                 else*/
                                //   InternalError(" No levelSet on Diff mesh :    to day  int1d of Matrix");
                                A += mate(t,10,Th[t].lab,stack,Q);
                            }
                            if(sptrclean) sptrclean=sptr->clean();
                        }
                    }
                }
            }
            else for( int e=0;e<Th.neb;e++)
            {
                if (all || setoflab.find(Th.bedges[e].lab) != setoflab.end())
                {
                    int ie,i =Th.BoundaryElement(e,ie);
                    A += mate(i,ie,Th.bedges[e].lab,stack);
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr

                }
            }
        }
        else if (di.kind == CDomainOfIntegration::intalledges)
        {
            for (int i=0;i< Th.nt; i++)
            {
                if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                for (int ie=0;ie<3;ie++)
                { // modif F.H to get the correct label in intalledges
                    int e0=VerticesOfTriangularEdge[ie][0];
                    int e1=VerticesOfTriangularEdge[ie][1];
                    int i1 = Th(Th[i][e0]),i2 = Th(Th[i][e1]);
                    BoundaryEdge * be = Th.TheBoundaryEdge(i1,i2);
                    int lab = be ? be->lab :  notalabel;

                    A += mate(i,ie,lab,paramate);
                }
                if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr

            }

        }
        else if (di.kind == CDomainOfIntegration::intallVFedges)
        {
            cerr << " a faire intallVFedges " << endl;
            ffassert(0);
            for (int i=0;i< Th.nt; i++)
            {
                if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                for (int ie=0;ie<3;ie++)
                A += mate(i,ie,Th[i].lab,paramate);
                if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr

            }

        }
        else if (di.kind == CDomainOfIntegration::int2d )
        {

            if(di.islevelset())
            {
                double uset = HUGE_VAL;
                R2 Q[2][3];
                double vol6[2];
                KN<double> phi(Th.nv);phi=uset;
                double f[3];
                for(int t=0; t< Th.nt;++t)
                {
                    if ( all || setoflab.find(Th[t].lab) != setoflab.end())
                    {
                        double umx=-HUGE_VAL,umn=HUGE_VAL;
                        for(int i=0;i<3;++i)
                        {
                            int j= ThI(t,i);
                            if( phi[j]==uset)
                            {
                                MeshPointStack(stack)->setP(&ThI,t,i);
                                phi[j]= di.levelset(stack);//zzzz
                            }
                            f[i]=phi[j];
                            umx = std::max(umx,phi[j]);
                            umn = std::min(umn,phi[j]);

                        }
                        int nt= UnderIso(f,Q, vol6,1e-14);
                        setQF<R2>(FIT,FITo,QuadratureFormular_T_1, Q,vol6,nt);
                        if(FIT.n)
                        A += mate(t,-1,Th[t].lab,stack);
                        if(sptrclean) sptrclean=sptr->clean();
                    }
                }
                FIT =FITo;
            }
            else


            for (int i=0;i< Th.nt; i++)
            {
                if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                A += mate(i,-1,Th[i].lab,stack);
                if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr

                // AA += mate;
            }
        }
        else
        InternalError(" kind of CDomainOfIntegration unknown");

        if (where_in_stack) delete [] where_in_stack;
        delete &mate;
        if(verbosity>9) cout << "  -- CPU assemble mat " <<  CPUtime()-CPU0 << " s\n";
    }

    // creating an instance of AssembleBilinearForm with MatriceCreuse
    // case 3D volume
    // --------- FH 120105
    template<class R>
    void AssembleBilinearForm(Stack stack,const FESpace3::Mesh & Th,const FESpace3 & Uh,const FESpace3 & Vh,bool sym,
                              MatriceCreuse<R>  & A, const  FormBilinear * b  )

    {
        /*FH:  case ..in 3D
         in varf ...
         standard case ..
         */


        typedef FESpace3 FESpace;
        typedef FESpace3::Mesh Mesh;
        typedef Mesh *pmesh ;
        StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
        Fem2D::MeshPoint & mp (*Fem2D::MeshPointStack(stack)), mps = mp;

        bool sptrclean=true;
        const CDomainOfIntegration & di= *b->di;
        ffassert(di.d==3);
        const Mesh * pThdi = GetAny<pmesh>( (* di.Th)(stack));
        if ( pThdi != &Th || &Uh.Th !=&Th || &Vh.Th !=&Th) {
            cout << " Use matrix formulation .... " << endl;
            ExecError("No way to compute bilinear form with integrale of on mesh \n"
                      "  test  or unknown function  defined on another mesh! sorry to hard.   ");
        }
        SHOWVERB(cout << " FormBilinear " << endl);
        MatriceElementaireSymetrique<R,FESpace> *mates =0;
        MatriceElementairePleine<R,FESpace> *matep =0;
        const int useopt=di.UseOpt(stack);
        double binside=di.binside(stack);
        if( di.withmap()) { ExecError(" no map  in the case (3)??");}

        //const vector<Expression>  & what(di.what);
        CDomainOfIntegration::typeofkind  kind = di.kind;
        set<int> setoflab;
        bool all=true;
        const QuadratureFormular1d & FIEo = di.FIE(stack);
        const QuadratureFormular & FITo = di.FIT(stack);
        const GQuadratureFormular<R3> & FIVo = di.FIV(stack);
        //  to change the quadrature on element ... may 2014 FH ..
        QuadratureFormular1d  FIE(FIEo,3);
        QuadratureFormular FIT(FITo,3);
        GQuadratureFormular<R3>  FIV(FIVo,3);

        bool VF=b->VF();  // finite Volume or discontinuous Galerkin
        if (verbosity>2) cout << "  -- discontinuous Galerkin  =" << VF << " size of Mat =" << A.size()<< " Bytes\n";
        if (verbosity>3)
        {
            if (CDomainOfIntegration::int2d==kind) cout << "  -- boundary int border ( nQP: "<< FIT.n << ") ,"  ;
            else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all edges ( nQP: "<< FIT.n << "),"  ;
            else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF edges nQP: ("<< FIE.n << ")," ;
            else cout << "  --  int3d   (nQP: "<< FIV.n << " ) in "  ;
            if(di.islevelset()) cout << " ( int on Levelset) " << endl;

        }
        if(di.islevelset() && (CDomainOfIntegration::int2d!=kind) && (CDomainOfIntegration::int3d!=kind))
        InternalError("Sorry no levelset integration type on no int[2|3]d case");

        Expandsetoflab(stack,di, setoflab,all);
        /*
         for (size_t i=0;i<what.size();i++)
         {
         long  lab  = GetAny<long>( (*what[i])(stack));
         setoflab.insert(lab);
         if ( verbosity>3) cout << lab << " ";
         all=false;
         }*/
        if (verbosity>3) cout <<" Optimized = "<< useopt << ", ";
        const E_F0 * poptiexp0=b->b->optiexp0;

        int n_where_in_stack_opt=b->b->where_in_stack_opt.size();
        R** where_in_stack =0;
        if (n_where_in_stack_opt && useopt)
        where_in_stack = new R * [n_where_in_stack_opt];
        if (where_in_stack)
        {
            assert(b->b->v.size()==(size_t) n_where_in_stack_opt);
            for (int i=0;i<n_where_in_stack_opt;i++)
            {
                int offset=b->b->where_in_stack_opt[i];
                assert(offset>10);
                where_in_stack[i]= static_cast<R *>(static_cast<void *>((char*)stack+offset));
                *(where_in_stack[i])=0;
            }


            if(poptiexp0)
            (*poptiexp0)(stack);
            KN<bool> ok(b->b->v.size());
            {  //   remove the zero coef in the liste
                // R zero=R();
                int il=0;
                for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
                ok[il] =  ! (b->b->mesh_indep_stack_opt[il] && ( std::norm(*(where_in_stack[il])) < 1e-100 ) );
            }
            BilinearOperator b_nozer(*b->b,ok);
            if (verbosity % 10 > 3 )
            cout << "   -- nb term in bilinear form  (!0) : " << b_nozer.v.size()
            << "  total " << n_where_in_stack_opt << endl;

            if ( (verbosity/100) % 10 >= 2)
            {
                int il=0;

                for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
                cout << il << " coef (" << l->first << ") = " << *(where_in_stack[il])
                << " offset=" << b->b->where_in_stack_opt[il]
                << " dep mesh " << l->second.MeshIndependent() << b->b->mesh_indep_stack_opt[il] << endl;
            }
        }
        Stack_Ptr<R*>(stack,ElemMatPtrOffset) =where_in_stack;
        void *paramate=stack;
        pair_stack_double parammatElement_OpVF(stack, & binside);
        parammatElement_OpVF.first = stack;
        parammatElement_OpVF.second= & binside;

        if (verbosity >3)
        {
            if (all) cout << " all " << endl ;
            else cout << endl;
        }
        if(VF) {
            if(&Uh != &Vh || sym)
            cout <<  ("To Day in bilinear form with discontinuous Galerkin (3d):   \n"
                      "  test or unknown function must be  defined on the same FEspace, \n"
                      "  and the matrix is not symmetric. \n"
                      " To do other case in a future (F. Hecht) dec. 2014 ");
            if(&Uh == &Vh)
            matep= new MatriceElementairePleine<R,FESpace>(Uh,VF,FIV,FIT,useopt);
            else
            matep= new MatriceElementairePleine<R,FESpace>(Uh,Vh,VF,FIV,FIT,useopt);
            matep->faceelement = Element_OpVF;
            paramate= &parammatElement_OpVF;
        }
        else if (sym) {
            mates= new MatriceElementaireSymetrique<R,FESpace>(Uh,FIV,FIT,useopt);
            mates->element = Element_Op<R>;
        }
        else {
            matep= new MatriceElementairePleine<R,FESpace>(Uh,Vh,FIV,FIT,useopt);
            matep->element = Element_Op<R>;
        }
        MatriceElementaireFES<R,FESpace> & mate(*( sym? (MatriceElementaireFES<R,FESpace> *)mates : (MatriceElementaireFES<R,FESpace> *) matep));


        mate.bilinearform=b->b;

        Check(*mate.bilinearform,mate.Uh.N,mate.Vh.N);

        if (di.kind == CDomainOfIntegration::int2d )
        {

            if(di.islevelset())
            {
                if(verbosity>99) cout << " int2d on levelset in 3d " << endl;
                double uset = HUGE_VAL;
                R3 Q[4];
                KN<double> phi(Th.nv);phi=uset;
                double f[4];
                for(int t=0; t< Th.nt;++t)
                {
                    if ( all || setoflab.find(Th[t].lab) != setoflab.end())
                    {
                        double umx=-HUGE_VAL,umn=HUGE_VAL;
                        for(int i=0;i<4;++i)
                        {
                            int j= Th(t,i);
                            if( phi[j]==uset)
                            {
                                MeshPointStack(stack)->setP(&Th,t,i);
                                phi[j]= di.levelset(stack);//zzzz
                            }
                            f[i]=phi[j];
                            umx = std::max(umx,phi[j]);
                            umn = std::min(umn,phi[j]);

                        }
                        if( umn <=0 && umx >= 0)
                        {
                            int np= IsoLineK(f,Q,1e-10);// ca code ...
                            //  cout <<umn << " " << umx << " " << np << endl;

                            if(np>2 )
                            {
                                if( verbosity > 999 ) cout << " -- int " << np << " on:  " << Q[0] << " " << Q[1] << " " << Q[2] << " " << Q[3] << endl;
                                A += mate(t,10+np,Th[t].lab,stack,Q);
                            }
                            if(sptrclean) sptrclean=sptr->clean();
                        }
                    }}

            }
            else
            for( int e=0;e<Th.nbe;e++)
            {
                if (all || setoflab.find(Th.be(e).lab) != setoflab.end())
                {
                    int ie,i =Th.BoundaryElement(e,ie);
                    A += mate(i,ie,Th.be(e).lab,stack);
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr

                }
            }
        }
        else if (di.kind == CDomainOfIntegration::intallfaces  )
        {
            for (int i=0;i< Th.nt; i++)
            {
                if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                for (int ie=0;ie<4;ie++)
                A += mate(i,ie,Th[i].lab,paramate);
                if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr

            }

        }
        else if (di.kind == CDomainOfIntegration::intallVFedges)
        {
            cerr << " a faire intallVFedges " << endl;
            ffassert(0);
            for (int i=0;i< Th.nt; i++)
            {
                if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                for (int ie=0;ie<3;ie++)
                A += mate(i,ie,Th[i].lab,paramate);
                if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr

            }

        }
        else if (di.kind == CDomainOfIntegration::int3d )
        {
            if(di.islevelset())  //  may 2014 FH ...
            {   // int3d levelset < 0
                double llevelset = 0;
                const double uset = std::numeric_limits<double>::max();
                // cout << " uset ="<<uset << endl;
                R3 Q[3][4];
                double vol6[3];
                KN<double> phi(Th.nv);
                phi=uset;
                double f[4];

                for (int t=0;t< Th.nt; t++)
                {

                    const Mesh3::Element & K(Th[t]);
                    if ( all || setoflab.find(Th[t].lab) != setoflab.end())

                    {
                        double umx=std::numeric_limits<double>::lowest(),umn=std::numeric_limits<double>::max();
                        for(int i=0;i<4;++i)
                        {
                            int j= Th(t,i);
                            if( phi[j]==uset)
                            {
                                MeshPointStack(stack)->setP(&Th,t,i);
                                phi[j]= di.levelset(stack);//zzzz
                            }
                            f[i]=phi[j];
                        }
                        int ntets= UnderIso(f,Q, vol6,1e-14);
                        setQF<R3>(FIV,FIVo,QuadratureFormular_Tet_1, Q,vol6,ntets);
                        if(FIV.n)
                        {
                            A += mate(t,-1,Th[t].lab,stack);
                            if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr

                        }
                    }
                }
                FIV=FIVo;
            }
            else
            for (int i=0;i< Th.nt; i++)
            {
                if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                A += mate(i,-1,Th[i].lab,stack);
                if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr

                // AA += mate;
            }

        }
        else
        {
            cerr << " kind of CDomainOfIntegration unknown ?? " << di.kind << endl;
            InternalError(" kind of CDomainOfIntegration unknown");
        }

        if (where_in_stack) delete [] where_in_stack;
        delete &mate;
        mp = mps;// restore data x,yz
    }



    // template struct to obtain the original type mesh - particular case for meshS FEM
  /*  template<class FE>    struct Trait_MESHO {
        // By default Mesh == MeshO
        typedef typename FE::Mesh MeshO;//   Mesh origin for Mesh S
        typedef typename FE::Mesh Mesh; //   Mesh Mesh S
        static Mesh * topmesh(MeshO *p) {return p;}
    };*/
   /* template<>    struct Trait_MESHO<FESpaceS> {
        // By default Mesh == MeshS and MeshO == Mesh3
        typedef  Mesh3 MeshO;
        typedef typename FESpaceS::Mesh Mesh;
        static Mesh * topmesh(MeshO *p) {return p->getMeshS();}
    };*/








    // creating an instance of AssembleBilinearForm with MatriceCreuse
    // case 3D surface
    template<class R>
    void AssembleBilinearForm(Stack stack,const MeshS & Th,const FESpaceS & Uh,const FESpaceS & Vh,bool sym,
                              MatriceCreuse<R>  & A, const  FormBilinear * b  )

    {
        /*FH:  case ..in 2D
         in varf ...
         standard case ..
         */

        typedef typename  Mesh::RdHat RdHat;
        StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
        bool sptrclean=true;
        const CDomainOfIntegration & di= *b->di;
         ffassert(di.d==3);

        SHOWVERB(cout << " FormBilinear " << endl);
        double CPU0 = CPUtime();
        MatriceElementaireSymetrique<R,FESpaceS> *mates =0;
        MatriceElementairePleine<R,FESpaceS> *matep =0;
        const int useopt=di.UseOpt(stack);
        double binside=di.binside(stack);

        //const vector<Expression>  & what(di.what);
        CDomainOfIntegration::typeofkind  kind = di.kind;
        set<int> setoflab;
        bool all=true;

        const MeshS & ThI = Th;// * GetAny<pmesh>( (* di.Th)(stack));
        bool sameMesh = &ThI == &Vh.Th &&  &ThI == &Uh.Th;

        //    const QuadratureFormular1d & FIE = di.FIE(stack);
        //    const QuadratureFormular & FIT = di.FIT(stack);
        const QuadratureFormular1d & FIEo = di.FIE(stack);
        const QuadratureFormular & FITo = di.FIT(stack);
        // const GQuadratureFormular<R3> & FIVo = di.FIV(stack);
        //  to change the quadrature on element ... may 2014 FH ..
        QuadratureFormular1d  FIE(FIEo,3);
        QuadratureFormular FIT(FITo,3);
        // GQuadratureFormular<R3>  FIV(FIVo,3);

        bool VF=b->VF();  // finite Volume or discontinuous Galerkin
        if (verbosity>2) cout << "  -- discontinuous Galerkin  =" << VF << " size of Mat =" << A.size()<< " Bytes\n";
        if (verbosity>3)
        {
            if (CDomainOfIntegration::int1d==kind) cout << "  -- boundary int border ( nQP: "<< FIE.n << ") ,"  ;
            else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all edges ( nQP: "<< FIE.n << "),"  ;
            else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF edges nQP: ("<< FIE.n << ")," ;
            else cout << "  --  int    (nQP: "<< FIT.n << " ) in "  ;
        }
        //if(di.islevelset()) InternalError("So no levelset integration type on this case (6)");
        if( di.withmap()) { ExecError(" no map  in the case (2)??");}
        if(di.islevelset() && ( (CDomainOfIntegration::int1d!=kind) && (CDomainOfIntegration::int2d!=kind) )  )
        InternalError("So no levelset integration type on no int1d case (6)");

        Expandsetoflab(stack,di, setoflab,all);
        /*
         for (size_t i=0;i<what.size();i++)
         {
         long  lab  = GetAny<long>( (*what[i])(stack));
         setoflab.insert(lab);
         if ( verbosity>3) cout << lab << " ";
         all=false;
         }*/
        if (verbosity>3) cout <<" Optimized = "<< useopt << ", ";
        const E_F0 * poptiexp0=b->b->optiexp0;

        int n_where_in_stack_opt=b->b->where_in_stack_opt.size();
        R** where_in_stack =0;
        if (n_where_in_stack_opt && useopt)
        where_in_stack = new R * [n_where_in_stack_opt];
        if (where_in_stack)
        {
            assert(b->b->v.size()==(size_t) n_where_in_stack_opt);
            for (int i=0;i<n_where_in_stack_opt;i++)
            {
                int offset=b->b->where_in_stack_opt[i];
                assert(offset>10);
                where_in_stack[i]= static_cast<R *>(static_cast<void *>((char*)stack+offset));
                *(where_in_stack[i])=0;
            }


            if(poptiexp0)
            (*poptiexp0)(stack);
            KN<bool> ok(b->b->v.size());
            {  //   remove the zero coef in the liste
                // R zero=R();
                int il=0;
                for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
                ok[il] =  ! (b->b->mesh_indep_stack_opt[il] && ( std::norm(*(where_in_stack[il])) < 1e-100 ) );
            }
            BilinearOperator b_nozer(*b->b,ok);
            if (verbosity % 10 > 3 )
            cout << "   -- nb term in bilinear form  (!0) : " << b_nozer.v.size()
            << "  total " << n_where_in_stack_opt << endl;

            if ( (verbosity/100) % 10 >= 2)
            {
                int il=0;

                for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
                cout << il << " coef (" << l->first << ") = " << *(where_in_stack[il])
                << " offset=" << b->b->where_in_stack_opt[il]
                << " dep mesh " << l->second.MeshIndependent() << b->b->mesh_indep_stack_opt[il] << endl;
            }
        }
        Stack_Ptr<R*>(stack,ElemMatPtrOffset) =where_in_stack;
        void *paramate=stack;
        pair_stack_double parammatElement_OpVF(stack,& binside);
        // parammatElement_OpVF.first = stack;
        // parammatElement_OpVF.second= & binside;

        if (verbosity >3)
        {
            if (all) cout << " all " << endl ;
            else cout << endl;
        }
        if(VF) {
            if(&Uh != &Vh || sym)
            cout << ("To Day in bilinear form with discontinuous Galerkin (2d):   \n"
                     "  test or unknown function must be  defined on the same FEspace, \n"
                     "  and the matrix is not symmetric. \n"
                     " To do other case in a future (F. Hecht) dec. 2003 ");
            if(&Uh == &Vh)
            matep= new MatriceElementairePleine<R,FESpaceS>(Uh,VF,FIT,FIE,useopt);
            else
            matep= new MatriceElementairePleine<R,FESpaceS>(Uh,Vh,VF,FIT,FIE,useopt);


            //      matep= new MatriceElementairePleine<R,FESpace>(Uh,Vh,VF,FIT,FIE);
            matep->faceelement = Element_OpVF;
            paramate= &parammatElement_OpVF;
        }
        else if (sym) {
            mates= new MatriceElementaireSymetrique<R,FESpaceS>(Uh,FIT,FIE,useopt);
            mates->element = Element_Op<R>;
        }
        else {
            matep= new MatriceElementairePleine<R,FESpaceS>(Uh,Vh,FIT,FIE,useopt);
            matep->element = Element_Op<R>;
        }
        MatriceElementaireFES<R,FESpaceS> & mate(*( sym? (MatriceElementaireFES<R,FESpaceS> *)mates : (MatriceElementaireFES<R,FESpaceS> *) matep));


        mate.bilinearform=b->b;

        Check(*mate.bilinearform,mate.Uh.N,mate.Vh.N);
        if(verbosity>9) cout << "  -- CPU init assemble mat " <<  CPUtime()-CPU0 << " s\n";

        if (di.kind == CDomainOfIntegration::int1d )
        {
            if(di.islevelset())
            {
                double uset = HUGE_VAL;
                R2 Q[3];
                KN<double> phi(Th.nv);phi=uset;
                double f[3];
                for(int t=0; t< Th.nt;++t)
                {
                    if ( all || setoflab.find(Th[t].lab) != setoflab.end())
                    {
                        double umx=-HUGE_VAL,umn=HUGE_VAL;
                        for(int i=0;i<3;++i)
                        {
                            int j= ThI(t,i);
                            if( phi[j]==uset)
                            {
                                MeshPointStack(stack)->setP(&ThI,t,i);
                                phi[j]= di.levelset(stack);//zzzz
                            }
                            f[i]=phi[j];
                            umx = std::max(umx,phi[j]);
                            umn = std::min(umn,phi[j]);

                        }
                        if( umn <=0 && umx >= 0)
                        {

                            int np= IsoLineK(f,Q,1e-10);
                            if(np==2)
                            {
                                /* if ( sameMesh)
                                 {

                                 Element_rhs<R>(Vh[t],*l->l,buf,stack,*B,FIE,Q[0],Q[1]);
                                 }
                                 else*/
                                //   InternalError(" No levelSet on Diff mesh :    to day  int1d of Matrix");
                                A += mate(t,10,Th[t].lab,stack,Q);
                            }
                            if(sptrclean) sptrclean=sptr->clean();
                        }
                    }
                }
            }
            else for( int e=0;e<Th.nbe;e++)
            {
                if (all || setoflab.find(Th.be(e).lab) != setoflab.end())
                {
                    int ie,i =Th.BoundaryElement(e,ie);
                    A += mate(i,ie,Th.be(e).lab,stack);
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr

                }
            }
        }
        else if (di.kind == CDomainOfIntegration::intalledges)
        {
            ffassert(0);
           /* for (int i=0;i< Th.nt; i++)
            {
                if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                for (int ie=0;ie<3;ie++)
                { // modif F.H to get the correct label in intalledges
                    int e0=VerticesOfTriangularEdge[ie][0];
                    int e1=VerticesOfTriangularEdge[ie][1];
                    int i1 = Th(Th[i][e0]),i2 = Th(Th[i][e1]);
                    BoundaryEdge * be = Th.TheBoundaryEdge(i1,i2);

                    int lab = be ? be->lab :  notalabel;

                    A += mate(i,ie,lab,paramate);
                }
                if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr

            }
           */
        }
       else if (di.kind == CDomainOfIntegration::intallVFedges)
        {
            ffassert(0);
           /* cerr << " a faire intallVFedges " << endl;
            ffassert(0);
            for (int i=0;i< Th.nt; i++)
            {
                if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                for (int ie=0;ie<3;ie++)
                A += mate(i,ie,Th[i].lab,paramate);
                if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr

            }
           */
        }
        else if (di.kind == CDomainOfIntegration::int2d )
        {
            if(di.islevelset())
            {
                double uset = HUGE_VAL;
                R2 Q[2][3];
                double vol6[2];
                KN<double> phi(Th.nv);phi=uset;
                double f[3];
                for(int t=0; t< Th.nt;++t)
                {
                    if ( all || setoflab.find(Th[t].lab) != setoflab.end())
                    {
                        double umx=-HUGE_VAL,umn=HUGE_VAL;
                        for(int i=0;i<3;++i)
                        {
                            int j= ThI(t,i);
                            if( phi[j]==uset)
                            {
                                MeshPointStack(stack)->setP(&ThI,t,i);
                                phi[j]= di.levelset(stack);//zzzz
                            }
                            f[i]=phi[j];
                            umx = std::max(umx,phi[j]);
                            umn = std::min(umn,phi[j]);

                        }
                        int nt= UnderIso(f,Q, vol6,1e-14);
                        setQF<R2>(FIT,FITo,QuadratureFormular_T_1, Q,vol6,nt);
                        if(FIT.n)
                        A += mate(t,-1,Th[t].lab,stack);
                        if(sptrclean) sptrclean=sptr->clean();
                    }
                }
                FIT =FITo;
            }
            else


            for (int i=0;i< Th.nt; i++)
            {
                if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                A += mate(i,-1,Th[i].lab,stack);
                if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr

                // AA += mate;
            }
        }
        else
        InternalError(" kind of CDomainOfIntegration unknown");

        if (where_in_stack) delete [] where_in_stack;
        delete &mate;
        if(verbosity>9) cout << "  -- CPU assemble mat " <<  CPUtime()-CPU0 << " s\n";
    }


    // creating an instance of AssembleBilinearForm with MatriceCreuse
    // case 3D curve
    template<class R>
    void AssembleBilinearForm(Stack stack,const MeshL & Th,const FESpaceL & Uh,const FESpaceL & Vh,bool sym,
                              MatriceCreuse<R>  & A, const  FormBilinear * b  )

    {

        typedef typename  Mesh::RdHat RdHat;
        StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
        bool sptrclean=true;
        const CDomainOfIntegration & di= *b->di;
        ffassert(di.d==3);

        SHOWVERB(cout << " FormBilinear " << endl);
        double CPU0 = CPUtime();
        MatriceElementaireSymetrique<R,FESpaceL> *mates =0;
        MatriceElementairePleine<R,FESpaceL> *matep =0;
        const int useopt=di.UseOpt(stack);
        double binside=di.binside(stack);

        //const vector<Expression>  & what(di.what);
        CDomainOfIntegration::typeofkind  kind = di.kind;
        set<int> setoflab;
        bool all=true;

        const MeshL & ThI = Th;// * GetAny<pmesh>( (* di.Th)(stack));
        bool sameMesh = &ThI == &Vh.Th &&  &ThI == &Uh.Th;

        const GQuadratureFormular<R1> & FITo = di.FIE(stack);
        GQuadratureFormular<R1>  FIT(FITo,3);

        bool VF=b->VF();  // finite Volume or discontinuous Galerkin
        if (verbosity>2) cout << "  -- discontinuous Galerkin  =" << VF << " size of Mat =" << A.size()<< " Bytes\n";

        if( di.withmap()) { ExecError(" no map  in the case (2)??");}
        if(di.islevelset() && ( (CDomainOfIntegration::int1d!=kind) && (CDomainOfIntegration::int2d!=kind) )  )
            InternalError("So no levelset integration type on no int1d case (6)");

        Expandsetoflab(stack,di, setoflab,all);

        if (verbosity>3) cout <<" Optimized = "<< useopt << ", ";
        const E_F0 * poptiexp0=b->b->optiexp0;

        int n_where_in_stack_opt=b->b->where_in_stack_opt.size();
        R** where_in_stack =0;
        if (n_where_in_stack_opt && useopt)
            where_in_stack = new R * [n_where_in_stack_opt];
        if (where_in_stack)
        {
            assert(b->b->v.size()==(size_t) n_where_in_stack_opt);
            for (int i=0;i<n_where_in_stack_opt;i++)
            {
                int offset=b->b->where_in_stack_opt[i];
                assert(offset>10);
                where_in_stack[i]= static_cast<R *>(static_cast<void *>((char*)stack+offset));
                *(where_in_stack[i])=0;
            }


            if(poptiexp0)
                (*poptiexp0)(stack);
            KN<bool> ok(b->b->v.size());
            {  //   remove the zero coef in the liste
                // R zero=R();
                int il=0;
                for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
                    ok[il] =  ! (b->b->mesh_indep_stack_opt[il] && ( std::norm(*(where_in_stack[il])) < 1e-100 ) );
            }
            BilinearOperator b_nozer(*b->b,ok);
            if (verbosity % 10 > 3 )
                cout << "   -- nb term in bilinear form  (!0) : " << b_nozer.v.size()
                << "  total " << n_where_in_stack_opt << endl;

            if ( (verbosity/100) % 10 >= 2)
            {
                int il=0;

                for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
                    cout << il << " coef (" << l->first << ") = " << *(where_in_stack[il])
                    << " offset=" << b->b->where_in_stack_opt[il]
                    << " dep mesh " << l->second.MeshIndependent() << b->b->mesh_indep_stack_opt[il] << endl;
            }
        }
        Stack_Ptr<R*>(stack,ElemMatPtrOffset) =where_in_stack;
        void *paramate=stack;
        pair_stack_double parammatElement_OpVF(stack,& binside);
        // parammatElement_OpVF.first = stack;
        // parammatElement_OpVF.second= & binside;

        if (verbosity >3)
        {
            if (all) cout << " all " << endl ;
            else cout << endl;
        }
        if(VF) {
            if(&Uh != &Vh || sym)
                cout << ("To Day in bilinear form with discontinuous Galerkin (2d):   \n"
                         "  test or unknown function must be  defined on the same FEspace, \n"
                         "  and the matrix is not symmetric. \n"
                         " To do other case in a future (F. Hecht) dec. 2003 ");
            if(&Uh == &Vh)
                matep= new MatriceElementairePleine<R,FESpaceL>(Uh,VF,FIT,0,useopt);
            else
                matep= new MatriceElementairePleine<R,FESpaceL>(Uh,Vh,VF,FIT,0,useopt);

            matep->faceelement = Element_OpVF;
            paramate= &parammatElement_OpVF;
        }
        else if (sym) {
            mates= new MatriceElementaireSymetrique<R,FESpaceL>(Uh,FIT,0,useopt);
            mates->element = Element_Op<R>;
        }
        else {
            matep= new MatriceElementairePleine<R,FESpaceL>(Uh,Vh,FIT,0,useopt);
            matep->element = Element_Op<R>;
        }
        MatriceElementaireFES<R,FESpaceL> & mate(*( sym? (MatriceElementaireFES<R,FESpaceL> *)mates : (MatriceElementaireFES<R,FESpaceL> *) matep));


        mate.bilinearform=b->b;

        Check(*mate.bilinearform,mate.Uh.N,mate.Vh.N);
        if(verbosity>9) cout << "  -- CPU init assemble mat " <<  CPUtime()-CPU0 << " s\n";

        if (di.kind == CDomainOfIntegration::int1d ) {
            for (int i=0;i< Th.nt; i++) {
                    if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                        A += mate(i,-1,Th[i].lab,stack);
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
            }
        }

        else if (di.kind == CDomainOfIntegration::int0d ) {
            for( int e=0;e<Th.nbe;e++) {
                if (all || setoflab.find(Th.be(e).lab) != setoflab.end()) {
                    int ie,i =Th.BoundaryElement(e,ie);
                    A += mate(i,ie,Th.be(e).lab,stack);
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                }
            }
        }
        else if (di.kind == CDomainOfIntegration::intall0d ) {// add FH juin 2021 
            for( int k=0;k<Th.nt;k++) {
                if (all || setoflab.find(Th[k].lab) != setoflab.end()) {
                    for( int ie=0;ie<2;ie++)
                    A += mate(k,ie,Th[k].lab,&parammatElement_OpVF);
                    if(sptrclean) sptrclean=sptr->clean();
                }
            }
        }

        
        else
            InternalError(" kind of CDomainOfIntegration unknown");

        if (where_in_stack) delete [] where_in_stack;
        delete &mate;
        if(verbosity>9) cout << "  -- CPU assemble mat " <<  CPUtime()-CPU0 << " s\n";
    }

    // 3D curve / 2D on meshL
    template<class R>
    void AssembleBilinearForm(Stack stack,const MeshL & Th,const FESpaceL & Uh,const FESpace & Vh,bool sym,
                          MatriceCreuse<R>  & A, const  FormBilinear * b  )

    {
      ffassert(0);
    }

    // 2D / 3D curve on meshL
    template<class R>
    void AssembleBilinearForm(Stack stack,const MeshL & Th,const FESpace & Uh,const FESpaceL & Vh,bool sym,
                         MatriceCreuse<R>  & A, const  FormBilinear * b  )
    {
      ffassert(0);
    }

    // 3D Surf / 3D volume on meshS
    template<class R>
    void AssembleBilinearForm(Stack stack,const MeshS & Th,const FESpaceS & Uh,const FESpace3 & Vh,bool sym,
                         MatriceCreuse<R>  & A, const  FormBilinear * b  )
    {
      ffassert(0);
    }
   // 3D volume / 3D Surf on meshS
   template<class R>
   void AssembleBilinearForm(Stack stack,const MeshS & Th,const FESpace3 & Uh,const FESpaceS & Vh,bool sym,
                        MatriceCreuse<R>  & A, const  FormBilinear * b  )
   {
     ffassert(0);
   }
   // 3D Surf / 3D curve on meshL
   template<class R>
   void AssembleBilinearForm(Stack stack,const MeshL & Th,const FESpaceS & Uh,const FESpaceL & Vh,bool sym,
                        MatriceCreuse<R>  & A, const  FormBilinear * b  )
   {
     ffassert(0);
   }
   // 3D curve / 3D Surf on meshL
   template<class R>
   void AssembleBilinearForm(Stack stack,const MeshL & Th,const FESpaceL & Uh,const FESpaceS & Vh,bool sym,
                        MatriceCreuse<R>  & A, const  FormBilinear * b  )
   {
     ffassert(0);
   }

// end 3d



    ////////////////////////////////////////////////
    // AddMatElem
    ////////////////////////////////////////////////

    // case 2d
    // --------- FH 170605

    template<class R>
    void  AddMatElem(MatriceMap<R> & A,const Mesh & Th,const BilinearOperator & Op,bool sym,int it,  int ie,int label,
                     const FESpace & Uh,const FESpace & Vh,
                     const QuadratureFormular & FI,
                     const QuadratureFormular1d & FIb,
                     double *p,   void *vstack, bool intmortar=false,R2 *Q=0)
    {
        //cout << "AddMatElem" << Q << " "  << ie << endl;
        Stack stack=pvoid2Stack(vstack);
        MeshPoint mp= *MeshPointStack(stack);
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        const Mesh & Thu(Uh.Th);
        const Mesh & Thv(Vh.Th);

        bool same = &Uh == & Vh;
        const Triangle & T  = Th[it];
        long npi;
        long i,j;
        bool classoptm = copt && Op.optiexpK;
        assert(Op.MaxOp() <last_operatortype);
        //


        KN<bool> Dop(last_operatortype);
        Op.DiffOp(Dop);
        int lastop=1+Dop.last([](bool x){return x;});
        //assert(lastop<=3);

        if (ie<0)
        {
            for (npi=0;npi<FI.n;npi++) // loop on the integration point
            {
                QuadraturePoint pi(FI[npi]);
                double coef = T.area*pi.a;
                R2 Pt(pi),Ptu,Ptv;
                R2 P(T(Pt));
                bool outsideu,outsidev;
                // ici trouve le T
                int iut=0,ivt=0;
                const Triangle * tu,*tv;
                if(&Th == & Thu )
                {
                    tu =&T;
                    Ptu=Pt;
                }
                else
                {
                    tu= Thu.Find(P,Ptu,outsideu);
                    if( !tu ||  outsideu) {
                        if(verbosity>100) cout << " On a pas trouver (u) " << P << " " << endl;
                        continue;}}
                if(same)
                {
                    tv=tu;
                    outsidev=outsideu;
                    Ptv=Ptu;
                }
                else
                {
                    if(&Th == & Thv )
                    {
                        tv =&T;
                        Ptv=Pt;
                    }
                    else
                    {
                        tv= Thv.Find(P,Ptv,outsidev);
                        if( !tv || outsidev) {
                            if(verbosity>100) cout << " On a pas trouver (v) " << P << " " << endl;
                            continue;
                        }}
                }
                iut = Thu(tu);
                ivt = Thv(tv);
                if( verbosity>1000) cout << " T " << it  << "  iut " << iut << " ivt " << ivt  <<  endl ;
                FElement Ku(Uh[iut]);
                FElement Kv(Vh[ivt]);
                long n= Kv.NbDoF() ,m=Ku.NbDoF();
                long N= Kv.N;
                long M= Ku.N;
                RNMK_ fv(p,n,N,lastop); //  the value for basic fonction
                RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction


                Ku.BF(Dop,Ptu,fu);
                MeshPointStack(stack)->set(Th,P,Pt,T,label);
                if (classoptm) (*Op.optiexpK)(stack); // call optim version
                if (!same) Kv.BF(Dop,Ptv,fv);
                for ( i=0;  i<n;   i++ )
                {

                    // attention la fonction test donne la ligne
                    //  et la fonction test est en second
                    int ig = Kv(i);
                    RNM_ wi(fv(i,'.','.'));
                    for ( j=0;  j<m;   j++ )
                    {
                        RNM_ wj(fu(j,'.','.'));
                        int il=0;
                        int jg(Ku(j));
                        if ( !sym ||  ig <= jg )
                        for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                        {  // attention la fonction test donne la ligne
                            //  et la fonction test est en second
                            BilinearOperator::K ll(*l);
                            pair<int,int> jj(ll.first.first),ii(ll.first.second);
                            double w_i =  wi(ii.first,ii.second);
                            double w_j =  wj(jj.first,jj.second);
                            R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack))   ;
                            if( verbosity>1000) cout << ig << " " << jg << " "  <<  " " << ccc << " " <<  coef * ccc * w_i*w_j << " on T \n"   ;
                            double wij =  w_i*w_j;
                            if (abs(wij)>= 1e-10)
                            A[make_pair(ig,jg)] += coef * ccc * wij;
                        }
                    }
                }
            }
        }
        else // int on edge ie
        {
            R2 PA,PB,E;
            if(Q)
            {
                PA=Q[0];
                PB=Q[1];
                E=T(PB)-T(PA);
                // cout << " AddMAtElem " <<  PA <<  " " << PB << " "<< sqrt((E,E))<< endl;
            }
            else
            {
                PA=TriangleHat[VerticesOfTriangularEdge[ie][0]];
                PB=TriangleHat[VerticesOfTriangularEdge[ie][1]];
                E=T.Edge(ie);
            }
            double le = sqrt((E,E));

            for (npi=0;npi<FIb.n;npi++) // loop on the integration point
            {
                QuadratureFormular1dPoint pi( FIb[npi]);
                double sa=pi.x,sb=1-sa;
                double coef = le*pi.a;

                R2 Pt(PA*sa+PB*sb ); //

                R2 Ptu,Ptv;
                R2 P(T(Pt));
                bool outsideu,outsidev;
                // ici trouve le T
                int iut=0,ivt=0;
                const Triangle * tu, *tv;
                if(&Th == & Thu )
                {
                    tu =&T;
                    Ptu=Pt;
                }
                else
                {
                    tu= Thu.Find(P,Ptu,outsideu);
                    if( !tu ||  (outsideu && !intmortar) )  {
                        //R dd=-1;
                        //if(tu) { R2 PP((*tu)(Ptu)),PPP(P,PP) ; cout << PP << " " << sqrt( (PPP,PPP) ) <<"    "; }
                        if(verbosity>100) cout << " On a pas trouver (u) " << P << " " <<Ptu << " " << tu <<   endl;
                        continue;}}
                iut = Thu(tu);
                if(same)
                {
                    tv=tu;
                    outsidev=outsideu;
                    Ptv=Ptu;
                    ivt=iut;
                }
                else
                {
                    if(&Th == & Thv )
                    {
                        tv =&T;
                        Ptv=Pt;
                    }
                    else {
                        tv= Thv.Find(P,Ptv,outsidev);
                        if( !tv || (outsidev&& !intmortar))  {
                            if(verbosity>100) cout << " On a pas trouver (v) " << P << " " << endl;
                            continue;}}
                    ivt = Thv(tv);
                }
                FElement Ku(Uh[iut]);
                FElement Kv(Vh[ivt]);
                long n= Kv.NbDoF() ,m=Ku.NbDoF();
                long N= Kv.N;
                long M= Ku.N;
                //  cout << P << " " <<  Pt << " " <<  iut << " " << ivt  << "  Ptu : " << Ptu << " Ptv: " << Ptv << " n:" << n << " m:" << m << endl;
                RNMK_ fv(p,n,N,lastop); //  the value for basic fonction
                RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction

                Ku.BF(Dop,Ptu,fu);
                if( !same)
                Kv.BF(Dop,Ptv,fv);


                // int label=-999999; // a passer en argument
                MeshPointStack(stack)->set(Th,P,Pt,T,label,R2(E.y,-E.x)/le,ie);
                if (classoptm) (*Op.optiexpK)(stack); // call optim version


                for ( i=0;  i<n;   i++ )
                // if (onWhatIsEdge[ie][Kv.DFOnWhat(i)]) // juste the df on edge bofbof generaly wrong FH dec 2003
                {
                    RNM_ wi(fv(i,'.','.'));
                    int ig=Kv(i);
                    for ( j=0;  j<m;   j++ )
                    {
                        RNM_ wj(fu(j,'.','.'));
                        int il=0;
                        int jg=Ku(j);
                        if( ! sym || ig <= jg )
                        for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                        {
                            BilinearOperator::K ll(*l);
                            pair<int,int> jj(ll.first.first),ii(ll.first.second);
                            double w_i =  wi(ii.first,ii.second);
                            double w_j =  wj(jj.first,jj.second);
                            // R ccc = GetAny<R>(ll.second.eval(stack));

                            R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                            double wij =  w_i*w_j;
                            if (abs(wij)>= 1e-10&& (verbosity>1000))
                            cout << " \t\t\t" << ig << " " << jg << " "  <<  ccc <<  " " <<  coef * ccc * w_i*w_j << " on edge \n" ;
                            if (abs(wij)>= 1e-10)
                            A[make_pair(ig,jg)] += wij*coef*ccc ;
                        }
                    }
                }
            }
        }

        *MeshPointStack(stack) = mp;
    }



    template<class R>
    void  AddMatElem(Expression const *const  mapu,Expression const * const mapt, MatriceMap<R> & A,const Mesh & Th,const BilinearOperator & Op,bool sym,int it,  int ie,int label,
                     const FESpace & Uh,const FESpace & Vh,
                     const QuadratureFormular & FI,
                     const QuadratureFormular1d & FIb,
                     double *p,   void *vstack, bool intmortar=false,R2 *Q=0)
    {
        //cout << "AddMatElem" << Q << " "  << ie << endl;
        Stack stack=pvoid2Stack(vstack);
        MeshPoint mp= *MeshPointStack(stack);
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        const Mesh & Thu(Uh.Th);
        const Mesh & Thv(Vh.Th);

        bool same = (&Uh == & Vh) && !mapu && !mapt;
        bool sameu =  &Th == & Thu && !mapu ;
        bool samev =  &Th == & Thv && !mapt ;
        const Triangle & T  = Th[it];
        long npi;
        long i,j;
        bool classoptm = copt && Op.optiexpK;
        assert(Op.MaxOp() <last_operatortype);
        //


        KN<bool> Dop(last_operatortype);
        Op.DiffOp(Dop);
        int lastop=1+Dop.last([](bool x){return x;});
        //assert(lastop<=3);

        if (ie<0)
        {
            for (npi=0;npi<FI.n;npi++) // loop on the integration point
            {
                QuadraturePoint pi(FI[npi]);
                double coef = T.area*pi.a;
                R2 Pt(pi),Ptu,Ptv;
                R2 P(T(Pt)),Pu(P),Pv(P);
                MeshPointStack(stack)->set(Th,P,Pt,T,label);

                if(mapu)
                Pu = R2( GetAny<double>((*mapu[0])(vstack)), GetAny<double>((*mapu[1])(vstack)));
                if(mapt)
                Pv = R2( GetAny<double>((*mapt[0])(vstack)), GetAny<double>((*mapt[1])(vstack)));
                if(verbosity>9999 && (mapu || mapt) )
                cout << " mapinng: " << P << " AddMatElem + map  -> (u) " << Pu << "  (t) ->"<< Pv << endl;
                bool outsideu,outsidev;
                // ici trouve le T
                int iut=0,ivt=0;
                const Triangle * tu,*tv;
                if(sameu )
                {
                    tu =&T;
                    Ptu=Pt;
                }
                else
                {
                    tu= Thu.Find(Pu,Ptu,outsideu);
                    if( !tu ||  outsideu) {
                        if(verbosity>100) cout << " On a pas trouver (u) " << P << " " << endl;
                        continue;}}
                if(same )
                {
                    tv=tu;
                    outsidev=outsideu;
                    Ptv=Ptu;
                }
                else if(samev)
                {
                    tv =&T;
                    Ptv=Pt;

                }
                else
                {
                    tv= Thv.Find(Pv,Ptv,outsidev);
                    if( !tv || outsidev) {
                        if(verbosity>100) cout << " On a pas trouver (v) " << P << " " << endl;
                        continue;
                    }
                }
                iut = Thu(tu);
                ivt = Thv(tv);
                if( verbosity>1000) cout << " T " << it  << "  iut " << iut << " ivt " << ivt  <<  endl ;
                FElement Ku(Uh[iut]);
                FElement Kv(Vh[ivt]);
                long n= Kv.NbDoF() ,m=Ku.NbDoF();
                long N= Kv.N;
                long M= Ku.N;
                RNMK_ fv(p,n,N,lastop); //  the value for basic fonction
                RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction


                Ku.BF(Dop,Ptu,fu);
                if (classoptm) (*Op.optiexpK)(stack); // call optim version
                if (!same) Kv.BF(Dop,Ptv,fv);
                for ( i=0;  i<n;   i++ )
                {

                    // attention la fonction test donne la ligne
                    //  et la fonction test est en second
                    int ig = Kv(i);
                    RNM_ wi(fv(i,'.','.'));
                    for ( j=0;  j<m;   j++ )
                    {
                        RNM_ wj(fu(j,'.','.'));
                        int il=0;
                        int jg(Ku(j));
                        if ( !sym ||  ig <= jg )
                        for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                        {  // attention la fonction test donne la ligne
                            //  et la fonction test est en second
                            BilinearOperator::K ll(*l);
                            pair<int,int> jj(ll.first.first),ii(ll.first.second);
                            double w_i =  wi(ii.first,ii.second);
                            double w_j =  wj(jj.first,jj.second);
                            R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack))   ;
                            if( verbosity>1000) cout << ig << " " << jg << " "  <<  " " << ccc << " " <<  coef * ccc * w_i*w_j << " on T \n"   ;
                            double wij =  w_i*w_j;
                            if (abs(wij)>= 1e-10)
                            A[make_pair(ig,jg)] += coef * ccc * wij;
                        }
                    }
                }
            }
        }
        else // int on edge ie
        {
            R2 PA,PB,E;
            if(Q)
            {
                PA=Q[0];
                PB=Q[1];
                E=T(PB)-T(PA);
                // cout << " AddMAtElem " <<  PA <<  " " << PB << " "<< sqrt((E,E))<< endl;
            }
            else
            {
                PA=TriangleHat[VerticesOfTriangularEdge[ie][0]];
                PB=TriangleHat[VerticesOfTriangularEdge[ie][1]];
                E=T.Edge(ie);
            }
            double le = sqrt((E,E));

            for (npi=0;npi<FIb.n;npi++) // loop on the integration point
            {
                QuadratureFormular1dPoint pi( FIb[npi]);
                double sa=pi.x,sb=1-sa;
                double coef = le*pi.a;

                R2 Pt(PA*sa+PB*sb ); //

                R2 Ptu,Ptv;
                R2 P(T(Pt)),Pu(P),Pv(P);
                MeshPointStack(stack)->set(Th,P,Pt,T,label,R2(E.y,-E.x)/le,ie);
                if(mapu)
                Pu = R2( GetAny<double>((*mapu[0])(vstack)), GetAny<double>((*mapu[1])(vstack)));
                if(mapt)
                Pv = R2( GetAny<double>((*mapt[0])(vstack)), GetAny<double>((*mapt[1])(vstack)));

                bool outsideu,outsidev;
                // ici trouve le T
                int iut=0,ivt=0;
                const Triangle * tu, *tv;
                if(sameu )
                {
                    tu =&T;
                    Ptu=Pt;
                }
                else
                {
                    tu= Thu.Find(Pu,Ptu,outsideu);
                    if( !tu ||  (outsideu && !intmortar) )  {
                        //R dd=-1;
                        //if(tu) { R2 PP((*tu)(Ptu)),PPP(P,PP) ; cout << PP << " " << sqrt( (PPP,PPP) ) <<"    "; }
                        if(verbosity>100) cout << " On a pas trouver (u) " << P << " " <<Ptu << " " << tu <<   endl;
                        continue;}}
                iut = Thu(tu);


                if(samev)
                {
                    tv =&T;
                    Ptv=Pt;
                }
                else {
                    tv= Thv.Find(Pv,Ptv,outsidev);
                    if( !tv || (outsidev&& !intmortar))  {
                        if(verbosity>100) cout << " On a pas trouver (v) " << P << " " << endl;
                        continue;}}
                ivt = Thv(tv);

                FElement Ku(Uh[iut]);
                FElement Kv(Vh[ivt]);
                long n= Kv.NbDoF() ,m=Ku.NbDoF();
                long N= Kv.N;
                long M= Ku.N;
                //  cout << P << " " <<  Pt << " " <<  iut << " " << ivt  << "  Ptu : " << Ptu << " Ptv: " << Ptv << " n:" << n << " m:" << m << endl;
                RNMK_ fv(p,n,N,lastop); //  the value for basic fonction
                RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction

                Ku.BF(Dop,Ptu,fu);
                if( !same)
                Kv.BF(Dop,Ptv,fv);


                // int label=-999999; // a passer en argument

                if (classoptm) (*Op.optiexpK)(stack); // call optim version


                for ( i=0;  i<n;   i++ )
                // if (onWhatIsEdge[ie][Kv.DFOnWhat(i)]) // juste the df on edge bofbof generaly wrong FH dec 2003
                {
                    RNM_ wi(fv(i,'.','.'));
                    int ig=Kv(i);
                    for ( j=0;  j<m;   j++ )
                    {
                        RNM_ wj(fu(j,'.','.'));
                        int il=0;
                        int jg=Ku(j);
                        if( ! sym || ig <= jg )
                        for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                        {
                            BilinearOperator::K ll(*l);
                            pair<int,int> jj(ll.first.first),ii(ll.first.second);
                            double w_i =  wi(ii.first,ii.second);
                            double w_j =  wj(jj.first,jj.second);
                            // R ccc = GetAny<R>(ll.second.eval(stack));

                            R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                            double wij =  w_i*w_j;
                            if (abs(wij)>= 1e-10&& (verbosity>1000))
                            cout << " \t\t\t" << ig << " " << jg << " "  <<  ccc <<  " " <<  coef * ccc * w_i*w_j << " on edge \n" ;
                            if (abs(wij)>= 1e-10)
                            A[make_pair(ig,jg)] += wij*coef*ccc ;
                        }
                    }
                }
            }
        }

        *MeshPointStack(stack) = mp;
    }


    //3D volume
    template<class R>
    void  AddMatElem(MatriceMap<R> & A,const Mesh3 & Th,const BilinearOperator & Op,bool sym,int it,  int ie,int label,
                     const FESpace3 & Uh,const FESpace3 & Vh,
                     const Fem2D::GQuadratureFormular<R3>  & FI,
                     const  QuadratureFormular & FIb,
                     double *p,   void *vstack, bool intmortar=false)
    {

        Stack stack=pvoid2Stack(vstack);
        MeshPoint mp= *MeshPointStack(stack);
        static int count =0; // non test FH .........................
        if(count++ < 1) {
            cout << " Warning : Assemble Matrix with incompatible 3d meshes in test (FH)  " << endl;
            cout << " ------------------------------------------------------------- " << endl;
        }
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        const Mesh3 & Thu(Uh.Th);
        const Mesh3 & Thv(Vh.Th);

        bool same = &Uh == & Vh;
        const Tet & T  = Th[it];
        long npi;
        long i,j;
        bool classoptm = copt && Op.optiexpK;
        assert(Op.MaxOp() <last_operatortype);
        //
        int lastop=0;
        lastop = 0;
        What_d Dop = Op.DiffOp(lastop);


        //assert(lastop<=3);

        if (ie<0)
        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            GQuadraturePoint<R3> pi(FI[npi]);
            double coef = T.mesure()*pi.a;
            R3 Pt(pi),Ptu,Ptv;
            R3 P(T(Pt));
            bool outsideu,outsidev;
            // ici trouve le T
            int iut=0,ivt=0;
            const Tet * tu,*tv;
            if(&Th == & Thu )
            {
                tu =&T;
                Ptu=Pt;
            }
            else
            {
                tu= Thu.Find(P,Ptu,outsideu);
                if( !tu ||  outsideu) {
                    if(verbosity>100) cout << " On a pas trouver (u) " << P << " " << endl;
                    continue;}}
            if(same)
            {
                tv=tu;
                outsidev=outsideu;
                Ptv=Ptu;
            }
            else
            {
                if(&Th == & Thv )
                {
                    tv =&T;
                    Ptv=Pt;
                }
                else
                {
                    tv= Thv.Find(P,Ptv,outsidev);
                    if( !tv || outsidev) {
                        if(verbosity>100) cout << " On a pas trouver (v) " << P << " " << endl;
                        continue;
                    }}
            }
            iut = Thu(tu);
            ivt = Thv(tv);
            if( verbosity>1000) cout << " T " << it  << "  iut " << iut << " ivt " << ivt  <<  endl ;
            FElement3 Ku(Uh[iut]);
            FElement3 Kv(Vh[ivt]);
            long n= Kv.NbDoF() ,m=Ku.NbDoF();
            long N= Kv.N;
            long M= Ku.N;
            RNMK_ fv(p,n,N,(long) lastop); //  the value for basic fonction
            RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,(long) lastop); //  the value for basic fonction


            Ku.BF(Dop,Ptu,fu);
            MeshPointStack(stack)->set(Th,P,Pt,T,label);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            if (!same) Kv.BF(Dop,Ptv,fv);
            for ( i=0;  i<n;   i++ )
            {

                // attention la fonction test donne la ligne
                //  et la fonction test est en second
                int ig = Kv(i);
                RNM_ wi(fv(i,'.','.'));
                for ( j=0;  j<m;   j++ )
                {
                    RNM_ wj(fu(j,'.','.'));
                    int il=0;
                    int jg(Ku(j));
                    if ( !sym ||  ig <= jg )
                    for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                    {  // attention la fonction test donne la ligne
                        //  et la fonction test est en second
                        BilinearOperator::K ll(*l);
                        pair<int,int> jj(ll.first.first),ii(ll.first.second);
                        double w_i =  wi(ii.first,ii.second);
                        double w_j =  wj(jj.first,jj.second);
                        R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack))   ;
                        if( verbosity>1000) cout << ig << " " << jg << " "  <<  " " << ccc << " " <<  coef * ccc * w_i*w_j << " on T \n"   ;
                        double wij =  w_i*w_j;
                        if (abs(wij)>= 1e-10)
                        A[make_pair(ig,jg)] += coef * ccc * wij;
                    }
                }
            }
        }
        else // int on edge ie
        for (npi=0;npi<FIb.n;npi++) // loop on the integration point
        {

            GQuadraturePoint<R2> pi( FIb[npi]);
            R3 NN= T.N(ie);
            double mes=NN.norme();
            NN/=mes;
            double coef = 0.5*mes*pi.a; // correction 0.5 050109 FH
            R3 Pt(T.PBord(ie,pi));
            //Ku.BF(Dop,Pt,fu);



            R3 Ptu,Ptv;
            R3 P(T(Pt));
            bool outsideu,outsidev;
            // ici trouve le T
            int iut=0,ivt=0;
            const Tet * tu, *tv;
            if(&Th == & Thu )
            {
                tu =&T;
                Ptu=Pt;
            }
            else
            {
                tu= Thu.Find(P,Ptu,outsideu);
                if( !tu ||  (outsideu && !intmortar) )  {
                    //R dd=-1;
                    //if(tu) { R2 PP((*tu)(Ptu)),PPP(P,PP) ; cout << PP << " " << sqrt( (PPP,PPP) ) <<"    "; }
                    if(verbosity>100) cout << " On a pas trouver (u) " << P << " " <<Ptu << " " << tu <<   endl;
                    continue;}}
            iut = Thu(tu);
            if(same)
            {
                tv=tu;
                outsidev=outsideu;
                Ptv=Ptu;
                ivt=iut;
            }
            else
            {
                if(&Th == & Thv )
                {
                    tv =&T;
                    Ptv=Pt;
                }
                else {
                    tv= Thv.Find(P,Ptv,outsidev);
                    if( !tv || (outsidev&& !intmortar))  {
                        if(verbosity>100) cout << " On a pas trouver (v) " << P << " " << endl;
                        continue;}}
                ivt = Thv(tv);
            }
            FElement3 Ku(Uh[iut]);
            FElement3 Kv(Vh[ivt]);
            long n= Kv.NbDoF() ,m=Ku.NbDoF();
            long N= Kv.N;
            long M= Ku.N;
            //  cout << P << " " <<  Pt << " " <<  iut << " " << ivt  << "  Ptu : " << Ptu << " Ptv: " << Ptv << " n:" << n << " m:" << m << endl;
            RNMK_ fv(p,n,N,lastop); //  the value for basic fonction
            RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction

            Ku.BF(Dop,Ptu,fu);
            if( !same)
            Kv.BF(Dop,Ptv,fv);


            // int label=-999999; // a passer en argument
            MeshPointStack(stack)->set(Th,P,Pt,T,label,NN,ie);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version


            for ( i=0;  i<n;   i++ )
            // if (onWhatIsEdge[ie][Kv.DFOnWhat(i)]) // juste the df on edge bofbof generaly wrong FH dec 2003
            {
                RNM_ wi(fv(i,'.','.'));
                int ig=Kv(i);
                for ( j=0;  j<m;   j++ )
                {
                    RNM_ wj(fu(j,'.','.'));
                    int il=0;
                    int jg=Ku(j);
                    if( ! sym || ig <= jg )
                    for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                    {
                        BilinearOperator::K ll(*l);
                        pair<int,int> jj(ll.first.first),ii(ll.first.second);
                        double w_i =  wi(ii.first,ii.second);
                        double w_j =  wj(jj.first,jj.second);
                        // R ccc = GetAny<R>(ll.second.eval(stack));

                        R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                        double wij =  w_i*w_j;
                        if (abs(wij)>= 1e-10&& (verbosity>1000))
                        cout << " \t\t\t" << ig << " " << jg << " "  <<  ccc <<  " " <<  coef * ccc * w_i*w_j << " on edge \n" ;
                        if (abs(wij)>= 1e-10)
                        A[make_pair(ig,jg)] += wij*coef*ccc ;
                    }
                }
            }
        }


        *MeshPointStack(stack) = mp;
    }


 // 3D surface case
    template<class R>
    void  AddMatElem(MatriceMap<R> & A,const MeshS & Th,const BilinearOperator & Op,bool sym,int it,  int ie,int label,
                     const FESpaceS & Uh,const FESpaceS & Vh,
                     const QuadratureFormular & FI,
                     const QuadratureFormular1d & FIb,
                     double *p,   void *vstack, bool intmortar=false)
    {
        Stack stack=pvoid2Stack(vstack);
        MeshPoint mp= *MeshPointStack(stack);
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        const MeshS & Thu(Uh.Th);
        const MeshS & Thv(Vh.Th);

        bool same = &Uh == & Vh;
        const TriangleS & T  = Th[it];
        long npi;
        long i,j;
        bool classoptm = copt && Op.optiexpK;
        assert(Op.MaxOp() <last_operatortype);
        //

        //KN<bool> Dop(last_operatortype);
        //Op.DiffOp(Dop);
        //int lastop=1+Dop.last(binder1st<equal_to<bool> >(equal_to<bool>(),true));
        //assert(lastop<=3);

        int lastop=0;
        What_d Dop = Op.DiffOp(lastop);

        if (ie<0)
        {
            for (npi=0;npi<FI.n;npi++) // loop on the integration point
            {
                QuadraturePoint pi(FI[npi]);
                double coef = T.mesure()*pi.a;
                R2 Pt(pi),Ptu,Ptv;
                R3 P(T(Pt));
                bool outsideu,outsidev;
                // ici trouve le T
                int iut=0,ivt=0;
                const TriangleS * tu,*tv;
                if(&Th == & Thu )
                {
                    tu =&T;
                    Ptu=Pt;
                }
                else
                {
                    tu= Thu.Find(P,Ptu,outsideu);
                    if( !tu ||  outsideu) {
                        if(verbosity>100) cout << " On a pas trouver (u) " << P << " " << endl;
                        continue;}}
                if(same)
                {
                    tv=tu;
                    outsidev=outsideu;
                    Ptv=Ptu;
                }
                else
                {
                    if(&Th == & Thv )
                    {
                        tv =&T;
                        Ptv=Pt;
                    }
                    else
                    {
                        tv= Thv.Find(P,Ptv,outsidev);
                        if( !tv || outsidev) {
                            if(verbosity>100) cout << " On a pas trouver (v) " << P << " " << endl;
                            continue;
                        }}
                }
                iut = Thu(tu);
                ivt = Thv(tv);
                if( verbosity>1000) cout << " T " << it  << "  iut " << iut << " ivt " << ivt  <<  endl ;
                FElementS Ku(Uh[iut]);
                FElementS Kv(Vh[ivt]);
                long n= Kv.NbDoF() ,m=Ku.NbDoF();
                long N= Kv.N;
                long M= Ku.N;
                RNMK_ fv(p,n,N,lastop); //  the value for basic fonction
                RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction


                Ku.BF(Dop,Ptu,fu);
                MeshPointStack(stack)->set(Th,P,Pt,T,label);
                if (classoptm) (*Op.optiexpK)(stack); // call optim version
                if (!same) Kv.BF(Dop,Ptv,fv);
                for ( i=0;  i<n;   i++ )
                {

                    // attention la fonction test donne la ligne
                    //  et la fonction test est en second
                    int ig = Kv(i);
                    RNM_ wi(fv(i,'.','.'));
                    for ( j=0;  j<m;   j++ )
                    {
                        RNM_ wj(fu(j,'.','.'));
                        int il=0;
                        int jg(Ku(j));
                        if ( !sym ||  ig <= jg )
                            for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                            {  // attention la fonction test donne la ligne
                                //  et la fonction test est en second
                                BilinearOperator::K ll(*l);
                                pair<int,int> jj(ll.first.first),ii(ll.first.second);
                                double w_i =  wi(ii.first,ii.second);
                                double w_j =  wj(jj.first,jj.second);
                                R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack))   ;
                                if( verbosity>1000) cout << ig << " " << jg << " "  <<  " " << ccc << " " <<  coef * ccc * w_i*w_j << " T \n"   ;
                                double wij =  w_i*w_j;
                                if (abs(wij)>= 1e-10)
                                A[make_pair(ig,jg)] += coef * ccc * wij;
                            }
                    }
                }
            }
        }
        else // int on edge ie
        {
          //  ffassert(0);
            R2 PA,PB;
            PA=TriangleHat[VerticesOfTriangularEdge[ie][0]];
            PB=TriangleHat[VerticesOfTriangularEdge[ie][1]];
            R3 E=T.Edge(ie);
            double le = sqrt((E,E));
        
            for (npi=0;npi<FIb.n;npi++) // loop on the integration point
            {
                QuadratureFormular1dPoint pi( FIb[npi]);
                double sa=pi.x,sb=1-sa;
                double coef = le*pi.a;
                R2 Pt(PA*sa+PB*sb );
                R2 Ptu,Ptv;
                R3 P(T(Pt));

                bool outsideu,outsidev;
                // ici trouve le T
                int iut=0,ivt=0;
                const TriangleS* tu, *tv;
                if(&Th == & Thu )
                {
                    tu =&T;
                    Ptu=Pt;
                }
                else
                {
                    tu= Thu.Find(P,Ptu,outsideu);        ////// probleme
                    if( !tu ||  (outsideu && !intmortar) )  {
                        //R dd=-1;
                        //if(tu) { R2 PP((*tu)(Ptu)),PPP(P,PP) ; cout << PP << " " << sqrt( (PPP,PPP) ) <<"    "; }
                        if(verbosity>100) cout << " On a pas trouver (u) " << P << " " <<Ptu << " " << tu <<   endl;
                        continue;}}

                iut = Thu(tu);

                if(same)
                {
                    tv=tu;
                    outsidev=outsideu;
                    Ptv=Ptu;
                    ivt=iut;
                }
                else
                {
                    if(&Th == & Thv )
                    {
                        tv =&T;
                        Ptv=Pt;
                    }
                    else {
                        tv= Thv.Find(P,Ptv,outsidev);
                        if( !tv || (outsidev&& !intmortar))  {
                            if(verbosity>100) cout << " On a pas trouver (v) " << P << " " << endl;
                            continue;}}
                    ivt = Thv(tv);
                }

                FElementS Ku(Uh[iut]);
                FElementS Kv(Vh[ivt]);
                long n= Kv.NbDoF() ,m=Ku.NbDoF();
                long N= Kv.N;
                long M= Ku.N;
                //  cout << P << " " <<  Pt << " " <<  iut << " " << ivt  << "  Ptu : " << Ptu << " Ptv: " << Ptv << " n:" << n << " m:" << m << endl;
                RNMK_ fv(p,n,N,lastop); //  the value for basic fonction
                RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction

                Ku.BF(Dop,Ptu,fu);
                if( !same)
                    Kv.BF(Dop,Ptv,fv);

                R3 NN= T.N(ie); //dHat=2
                // int label=-999999; // a passer en argument
                MeshPointStack(stack)->set(Th,P,Pt,T,label,NN,ie);
                if (classoptm) (*Op.optiexpK)(stack); // call optim version


                for ( i=0;  i<n;   i++ )
                    // if (onWhatIsEdge[ie][Kv.DFOnWhat(i)]) // juste the df on edge bofbof generaly wrong FH dec 2003
                {
                    RNM_ wi(fv(i,'.','.'));
                    int ig=Kv(i);
                    for ( j=0;  j<m;   j++ )
                    {
                        RNM_ wj(fu(j,'.','.'));
                        int il=0;
                        int jg=Ku(j);
                        if( ! sym || ig <= jg )
                            for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                            {
                                BilinearOperator::K ll(*l);
                                pair<int,int> jj(ll.first.first),ii(ll.first.second);
                                double w_i =  wi(ii.first,ii.second);
                                double w_j =  wj(jj.first,jj.second);
                                // R ccc = GetAny<R>(ll.second.eval(stack));

                                R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                                double wij =  w_i*w_j;
                                if (abs(wij)>= 1e-10&& (verbosity>1000))
                                    cout << " \t\t\t" << ig << " " << jg << " "  <<  ccc <<  " " <<  coef * ccc * w_i*w_j << " on edge \n" ;
                                if (abs(wij)>= 1e-10)
                                    A[make_pair(ig,jg)] += wij*coef*ccc ;
                            }
                    }
                }}
            }


        *MeshPointStack(stack) = mp;
    }




    // 3D curve case
    template<class R>
    void  AddMatElem(MatriceMap<R> & A,const MeshL & Th,const BilinearOperator & Op,bool sym,int it,  int ie,int label,
                     const FESpaceL & Uh,const FESpaceL & Vh,
                     const GQuadratureFormular<R1> & FI,
                     const QuadratureFormular1d & FIb,
                     double *p,   void *vstack, bool intmortar=false)
    {
        Stack stack=pvoid2Stack(vstack);
        MeshPoint mp= *MeshPointStack(stack);
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        const MeshL & Thu(Uh.Th);
        const MeshL & Thv(Vh.Th);

        bool same = &Uh == & Vh;
        const EdgeL & T  = Th[it];
        long npi;
        long i,j;
        bool classoptm = copt && Op.optiexpK;
        assert(Op.MaxOp() <last_operatortype);

        int lastop=0;
        What_d Dop = Op.DiffOp(lastop);

        if (ie<0)
        {
            for (npi=0;npi<FI.n;npi++) // loop on the integration point
            {
                GQuadraturePoint<R1> pi(FI[npi]);
                double coef = T.mesure()*pi.a;
                R1 Pt(pi),Ptu,Ptv;
                R3 P(T(Pt));
                bool outsideu,outsidev;
                // ici trouve le T
                int iut=0,ivt=0;
                const EdgeL * tu,*tv;
                if(&Th == & Thu ) {
                    tu =&T;
                    Ptu=Pt;
                }
                else {
                    tu= Thu.Find(P,Ptu,outsideu);
                    if( !tu ||  outsideu) {
                        if(verbosity>100) cout << " On a pas trouver (u) " << P << " " << endl;
                        continue;}
                }
                if(same) {
                    tv=tu;
                    outsidev=outsideu;
                    Ptv=Ptu;
                }
                else {
                    if(&Th == & Thv ) {
                        tv =&T;
                        Ptv=Pt;
                    }
                    else {
                        tv= Thv.Find(P,Ptv,outsidev);
                        if( !tv || outsidev) {
                            if(verbosity>100) cout << " On a pas trouver (v) " << P << " " << endl;
                            continue;
                        }
                    }
                }
                iut = Thu(tu);
                ivt = Thv(tv);
                if( verbosity>1000) cout << " T " << it  << "  iut " << iut << " ivt " << ivt  <<  endl ;
                FElementL Ku(Uh[iut]), Kv(Vh[ivt]);
                long n= Kv.NbDoF() ,m=Ku.NbDoF();
                long N= Kv.N;
                long M= Ku.N;
                RNMK_ fv(p,n,N,lastop); //  the value for basic fonction
                RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction


                Ku.BF(Dop,Ptu,fu);
                R3 NNt=T.TangenteUnitaire();
                MeshPointStack(stack)->set(Th,P,Pt,T,NNt,label);
                if (classoptm) (*Op.optiexpK)(stack); // call optim version
                if (!same) Kv.BF(Dop,Ptv,fv);
                for ( i=0;  i<n;   i++ ) {
                    // attention la fonction test donne la ligne
                    //  et la fonction test est en second
                    int ig = Kv(i);
                    RNM_ wi(fv(i,'.','.'));
                    for ( j=0;  j<m;   j++ ) {
                        RNM_ wj(fu(j,'.','.'));
                        int il=0;
                        int jg(Ku(j));
                        if ( !sym ||  ig <= jg )
                            for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++) {
                                // attention la fonction test donne la ligne
                                //  et la fonction test est en second
                                BilinearOperator::K ll(*l);
                                pair<int,int> jj(ll.first.first),ii(ll.first.second);
                                double w_i =  wi(ii.first,ii.second);
                                double w_j =  wj(jj.first,jj.second);
                                R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack))   ;
                                if( verbosity>1000) cout << ig << " " << jg << " "  <<  " " << ccc << " " <<  coef * ccc * w_i*w_j << " T \n"   ;
                                double wij =  w_i*w_j;
                                if (abs(wij)>= 1e-10)
                                    A[make_pair(ig,jg)] += coef * ccc * wij;
                            }
                    }
                }
            }
        }
        else // int on point ie
        //{
            ffassert(0);
            /*R2 PA,PB;
            PA=TriangleHat[VerticesOfTriangularEdge[ie][0]];
            PB=TriangleHat[VerticesOfTriangularEdge[ie][1]];
            R3 E=T.Edge(ie);
            double le = sqrt((E,E));
            for (npi=0;npi<FIb.n;npi++) // loop on the integration point
            {
                QuadratureFormular1dPoint pi( FIb[npi]);
                double sa=pi.x,sb=1-sa;
                double coef = le*pi.a;
                R2 Pt(PA*sa+PB*sb );
                R2 Ptu,Ptv;
                R3 P(T(Pt));

                bool outsideu,outsidev;
                // ici trouve le T
                int iut=0,ivt=0;
                const TriangleS* tu, *tv;
                if(&Th == & Thu )
                {
                    tu =&T;
                    Ptu=Pt;
                }
                else
                {
                    tu= Thu.Find(P,Ptu,outsideu);        ////// probleme
                    if( !tu ||  (outsideu && !intmortar) )  {
                        //R dd=-1;
                        //if(tu) { R2 PP((*tu)(Ptu)),PPP(P,PP) ; cout << PP << " " << sqrt( (PPP,PPP) ) <<"    "; }
                        if(verbosity>100) cout << " On a pas trouver (u) " << P << " " <<Ptu << " " << tu <<   endl;
                        continue;}}

                iut = Thu(tu);

                if(same)
                {
                    tv=tu;
                    outsidev=outsideu;
                    Ptv=Ptu;
                    ivt=iut;
                }
                else
                {
                    if(&Th == & Thv )
                    {
                        tv =&T;
                        Ptv=Pt;
                    }
                    else {
                        tv= Thv.Find(P,Ptv,outsidev);
                        if( !tv || (outsidev&& !intmortar))  {
                            if(verbosity>100) cout << " On a pas trouver (v) " << P << " " << endl;
                            continue;}}
                    ivt = Thv(tv);
                }

                FElementS Ku(Uh[iut]);
                FElementS Kv(Vh[ivt]);
                long n= Kv.NbDoF() ,m=Ku.NbDoF();
                long N= Kv.N;
                long M= Ku.N;
                //  cout << P << " " <<  Pt << " " <<  iut << " " << ivt  << "  Ptu : " << Ptu << " Ptv: " << Ptv << " n:" << n << " m:" << m << endl;
                RNMK_ fv(p,n,N,lastop); //  the value for basic fonction
                RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction

                Ku.BF(Dop,Ptu,fu);
                if( !same)
                    Kv.BF(Dop,Ptv,fv);

                R3 NN= T.N(ie); //dHat=2
                // int label=-999999; // a passer en argument
                MeshPointStack(stack)->set(Th,P,Pt,T,label,NN,ie);
                if (classoptm) (*Op.optiexpK)(stack); // call optim version


                for ( i=0;  i<n;   i++ )
                    // if (onWhatIsEdge[ie][Kv.DFOnWhat(i)]) // juste the df on edge bofbof generaly wrong FH dec 2003
                {
                    RNM_ wi(fv(i,'.','.'));
                    int ig=Kv(i);
                    for ( j=0;  j<m;   j++ )
                    {
                        RNM_ wj(fu(j,'.','.'));
                        int il=0;
                        int jg=Ku(j);
                        if( ! sym || ig <= jg )
                            for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                            {
                                BilinearOperator::K ll(*l);
                                pair<int,int> jj(ll.first.first),ii(ll.first.second);
                                double w_i =  wi(ii.first,ii.second);
                                double w_j =  wj(jj.first,jj.second);
                                // R ccc = GetAny<R>(ll.second.eval(stack));

                                R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                                double wij =  w_i*w_j;
                                if (abs(wij)>= 1e-10&& (verbosity>1000))
                                    cout << " \t\t\t" << ig << " " << jg << " "  <<  ccc <<  " " <<  coef * ccc * w_i*w_j << " on edge \n" ;
                                if (abs(wij)>= 1e-10)
                                    A[make_pair(ig,jg)] += wij*coef*ccc ;
                            }
                    }
                }}
        }*/


        *MeshPointStack(stack) = mp;
    }


// 3D curve / 2D on meshL
template<class R>
void  AddMatElem(MatriceMap<R> & A,const MeshL & Th,const BilinearOperator & Op,bool sym,int it,int ie,int label,
                 const FESpaceL & Uh,const FESpace & Vh,const GQuadratureFormular<R1> & FI,const QuadratureFormular1d & FIb,
                 double *p,void *vstack, bool intmortar=false)
{

    Stack stack=pvoid2Stack(vstack);
    MeshPoint mp= *MeshPointStack(stack);
    R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
    const MeshL & Thu(Uh.Th);
    const Mesh & Thv(Vh.Th);

    bool same = false;
    const EdgeL &T = Th[it];
    long npi;
    long i,j;
    bool classoptm = copt && Op.optiexpK;
    assert(Op.MaxOp() <last_operatortype);

    int lastop=0;
    What_d Dop = Op.DiffOp(lastop);
    KN<bool> Dop2(last_operatortype);
    Op.DiffOp(Dop2);
    int lastop2=lastop;//1+Dop2.last(binder1st<equal_to<bool> >(equal_to<bool>(),true));
    double epsP=1e-6; // must be choose
    
    if (ie<0)
    {
        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            GQuadraturePoint<R1> pi(FI[npi]);
            double coef = T.mesure()*pi.a;
            R1 Pt(pi),Ptu;
            R2 Ptv;
            R3 P(T(Pt));
            R2 PP(P.p2());
            bool outsideu,outsidev;
            // ici trouve le T
            int iut=0,ivt=0;
            const EdgeL * tu;
            const Triangle *tv;
            if(&Th == & Thu) {
                tu =&T;
                Ptu=Pt;
            }
            else {
                tu= Thu.Find(P,Ptu,outsideu);
                if( !tu ||  outsideu) {
                    if(verbosity>100) cout << " On a pas trouver (u) " << P << " " << endl;
                    continue;
                }
            }
            if(abs(P.z)>epsP) {outsidev=true;tv=0;}
            else tv= Thv.Find(PP,Ptv,outsidev);
            if( !tv || outsidev) {
                if(verbosity>100) cout << " On a pas trouver (v) " << P << " " << endl;
                continue;
            }
            iut = Thu(tu);
            ivt = Thv(tv);
            if( verbosity>1000) cout << " T " << it  << "  iut " << iut << " ivt " << ivt  <<  endl ;
            FElementL Ku(Uh[iut]);
            FElement Kv(Vh[ivt]);
            long n= Kv.NbDoF() ,m=Ku.NbDoF();
            long N= Kv.N, M= Ku.N;
            RNMK_ fv(p,n,N,lastop2); //  the value for basic fonction       // AXEL  lastop2
            RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction     // lastop
          //  cout << " same ?0:n*N*lastop " << same ?0:n*N*lastop << endl;

            Ku.BF(Dop,Ptu,fu);
            MeshPointStack(stack)->set(Th,P,Pt,T,label);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            if (!same) Kv.BF(Dop2,Ptv,fv);
            for ( i=0;  i<n;   i++ ) {
                // attention la fonction test donne la ligne
                //  et la fonction test est en second
                int ig = Kv(i);
                RNM_ wi(fv(i,'.','.'));
                for ( j=0;  j<m;   j++ ) {
                    RNM_ wj(fu(j,'.','.'));
                    int il=0;
                    int jg(Ku(j));
                    if ( !sym ||  ig <= jg )
                        for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++) {
                            // attention la fonction test donne la ligne
                            //  et la fonction test est en second
                            BilinearOperator::K ll(*l);
                            pair<int,int> jj(ll.first.first),ii(ll.first.second);
                            double w_i =  wi(ii.first,ii.second);
                            double w_j =  wj(jj.first,jj.second);
                            R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                            if( verbosity>1000) cout << ig << " " << jg << " "  <<  " " << ccc << " " <<  coef * ccc * w_i*w_j << " T \n"   ;
                            double wij =  w_i*w_j;
                            if (abs(wij)>= 1e-10)
                                A[make_pair(ig,jg)] += coef * ccc * wij;
                        }
                }
            }
        }
    }
    else // int on point ie
        ffassert(0);
 
    *MeshPointStack(stack) = mp;
}

// 2D / 3D curve on meshL
template<class R>
void  AddMatElem(MatriceMap<R> & A,const MeshL & Th,const BilinearOperator & Op,bool sym,int it,int ie,int label,
                 const FESpace & Uh,const FESpaceL & Vh,const GQuadratureFormular<R1> & FI,const QuadratureFormular1d & FIb,
                 double *p,void *vstack, bool intmortar=false)
{

    Stack stack=pvoid2Stack(vstack);
    MeshPoint mp= *MeshPointStack(stack);
    R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
    const Mesh & Thu(Uh.Th);
    const MeshL & Thv(Vh.Th);

    bool same = false;
    const EdgeL &T = Th[it];
    long npi;
    long i,j;
    bool classoptm = copt && Op.optiexpK;
    assert(Op.MaxOp() <last_operatortype);

    int lastop=0;
    What_d Dop = Op.DiffOp(lastop);
    KN<bool> Dop2(last_operatortype);
    Op.DiffOp(Dop2);
    int lastop2=lastop;//1+Dop2.last(binder1st<equal_to<bool> >(equal_to<bool>(),true));
    double epsP=1e-6; // must be choose
    
    if (ie<0) {
        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            GQuadraturePoint<R1> pi(FI[npi]);
            double coef = T.mesure()*pi.a;
            R1 Pt(pi),Ptv;
            R2 Ptu;
            R3 P(T(Pt));
            R2 PP(P.p2());
            bool outsideu,outsidev;
            // ici trouve le T
            int iut=0,ivt=0;
            const Triangle * tu;
            const EdgeL *tv;
            if(abs(P.z)>epsP) {outsidev=true;tu=0;}
            else tu= Thu.Find(PP,Ptu,outsideu);
              if( !tu ||  outsideu) {
                if(verbosity>100) cout << " On a pas trouver (u) " << P << " " << endl;
                continue;
              }
                        
            if(&Th == & Thv ) {
              tv =&T;
              Ptv=Pt;
            }
            else {
              tv= Thv.Find(P,Ptv,outsideu);
              if( !tv ||  outsideu) {
                if(verbosity>100) cout << " On a pas trouver (u) " << P << " " << endl;
                continue;}
            }
            
            iut = Thu(tu);
            ivt = Thv(tv);
            
            if( verbosity>1000) cout << " T " << it  << "  iut " << iut << " ivt " << ivt  <<  endl ;
            FElement Ku(Uh[iut]);
            FElementL Kv(Vh[ivt]);
            long n= Kv.NbDoF() ,m=Ku.NbDoF();
            long N= Kv.N, M= Ku.N;
            RNMK_ fv(p,n,N,lastop2); //  the value for basic fonction       // AXEL  lastop2
            RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction     // lastop
          //  cout << " same ?0:n*N*lastop " << same ?0:n*N*lastop << endl;

            Ku.BF(Dop2,Ptu,fu);
            MeshPointStack(stack)->set(Th,P,Pt,T,label);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            if (!same) Kv.BF(Dop,Ptv,fv);
            for ( i=0;  i<n;   i++ ) {
                // attention la fonction test donne la ligne
                //  et la fonction test est en second
                int ig = Kv(i);
                RNM_ wi(fv(i,'.','.'));
                for ( j=0;  j<m;   j++ ) {
                    RNM_ wj(fu(j,'.','.'));
                    int il=0;
                    int jg(Ku(j));
                    if ( !sym ||  ig <= jg )
                        for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++) {
                            // attention la fonction test donne la ligne
                            //  et la fonction test est en second
                            BilinearOperator::K ll(*l);
                            pair<int,int> jj(ll.first.first),ii(ll.first.second);
                            double w_i =  wi(ii.first,ii.second);
                            double w_j =  wj(jj.first,jj.second);
                            R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                            if( verbosity>1000) cout << ig << " " << jg << " "  <<  " " << ccc << " " <<  coef * ccc * w_i*w_j << " T \n"   ;
                            double wij =  w_i*w_j;
                            if (abs(wij)>= 1e-10)
                                A[make_pair(ig,jg)] += coef * ccc * wij;
                        }
                }
            }
        }
    }
    else // int on point ie
        ffassert(0);
 
    *MeshPointStack(stack) = mp;
}

// 3D Surf / 3D volume on meshS
template<class R>
void  AddMatElem(MatriceMap<R> & A,const MeshS & Th,const BilinearOperator & Op,bool sym,int it,int ie,int label,
                 const FESpaceS & Uh,const FESpace3 & Vh,const QuadratureFormular & FI,const QuadratureFormular1d & FIb,
                 double *p,void *vstack, bool intmortar=false)
{

    Stack stack=pvoid2Stack(vstack);
    MeshPoint mp= *MeshPointStack(stack);
    R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
    const MeshS & Thu(Uh.Th);
    const Mesh3 & Thv(Vh.Th);

    bool same = false;
    const TriangleS &T = Th[it];
    long npi;
    long i,j;
    bool classoptm = copt && Op.optiexpK;
    assert(Op.MaxOp() <last_operatortype);

    int lastop=0;
    What_d Dop = Op.DiffOp(lastop);
    double epsP=1e-6; // must be choose
    
    if (ie<0)
    {
        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            QuadraturePoint pi(FI[npi]);
            double coef = T.mesure()*pi.a;
            R2 Pt(pi),Ptu;
            R3 Ptv;
            R3 P(T(Pt));
            R3 PP(P); // correction FH. Sep 2022.. 
            bool outsideu,outsidev;
            // ici trouve le T
            int iut=0,ivt=0;
            const TriangleS * tu;
            const Tet *tv;
            if(&Th == & Thu) {
                tu =&T;
                Ptu=Pt;
            }
            else {
                tu= Thu.Find(P,Ptu,outsideu);
                if( !tu ||  outsideu) {
                    if(verbosity>100) cout << " On a pas trouver (u) " << P << " " << endl;
                    continue;
                }
            }
            tv= Thv.Find(PP,Ptv,outsidev);
            if( !tv || outsidev) {
                if(verbosity>100) cout << " On a pas trouver (v) " << P << " " << endl;
                continue;
            }
            iut = Thu(tu);
            ivt = Thv(tv);
            if( verbosity>1000) cout << " T " << it  << "  iut " << iut << " ivt " << ivt  <<  endl ;
            FElementS Ku(Uh[iut]);
            FElement3 Kv(Vh[ivt]);
            long n= Kv.NbDoF() ,m=Ku.NbDoF();
            long N= Kv.N, M= Ku.N;
            RNMK_ fv(p,n,N,lastop); //  the value for basic fonction
            RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction     // lastop
          //  cout << " same ?0:n*N*lastop " << same ?0:n*N*lastop << endl;

            Ku.BF(Dop,Ptu,fu);
            MeshPointStack(stack)->set(Th,P,Pt,T,label);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            if (!same) Kv.BF(Dop,Ptv,fv);
            for ( i=0;  i<n;   i++ ) {
                // attention la fonction test donne la ligne
                //  et la fonction test est en second
                int ig = Kv(i);
                RNM_ wi(fv(i,'.','.'));
                for ( j=0;  j<m;   j++ ) {
                    RNM_ wj(fu(j,'.','.'));
                    int il=0;
                    int jg(Ku(j));
                    if ( !sym ||  ig <= jg )
                        for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++) {
                            // attention la fonction test donne la ligne
                            //  et la fonction test est en second
                            BilinearOperator::K ll(*l);
                            pair<int,int> jj(ll.first.first),ii(ll.first.second);
                            double w_i =  wi(ii.first,ii.second);
                            double w_j =  wj(jj.first,jj.second);
                            R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                            if( verbosity>1000) cout << ig << " " << jg << " "  <<  " " << ccc << " " <<  coef * ccc * w_i*w_j << " T \n"   ;
                            double wij =  w_i*w_j;
                            if (abs(wij)>= 1e-10)
                                A[make_pair(ig,jg)] += coef * ccc * wij;
                        }
                }
            }
        }
    }
    else // int on point ie
        ffassert(0);
 
    *MeshPointStack(stack) = mp;
}
// 3D volume / 3D Surf on meshS
template<class R>
void  AddMatElem(MatriceMap<R> & A,const MeshS & Th,const BilinearOperator & Op,bool sym,int it,int ie,int label,
                 const FESpace3 & Uh,const FESpaceS & Vh,const QuadratureFormular & FI,const QuadratureFormular1d & FIb,
                 double *p,void *vstack, bool intmortar=false)
{

    Stack stack=pvoid2Stack(vstack);
    MeshPoint mp= *MeshPointStack(stack);
    R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
    const Mesh3 & Thu(Uh.Th);
    const MeshS & Thv(Vh.Th);

    bool same = false;
    const TriangleS &T = Th[it];
    long npi;
    long i,j;
    bool classoptm = copt && Op.optiexpK;
    assert(Op.MaxOp() <last_operatortype);

    int lastop=0;
    What_d Dop = Op.DiffOp(lastop);
    
    if (ie<0)
    {
        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            QuadraturePoint pi(FI[npi]);
            double coef = T.mesure()*pi.a;
            R2 Pt(pi),Ptv;
            R3 Ptu;
            R3 P(T(Pt));
            R2 PP(P.p2());
            bool outsideu,outsidev;
            // ici trouve le T
            int iut=0,ivt=0;
            const Tet * tu;
            const TriangleS *tv;
            tu= Thu.Find(P,Ptu,outsideu);
              if( !tu ||  outsideu) {
                if(verbosity>100) cout << " On a pas trouver (u) " << P << " " << endl;
                continue;
              }
                        
            if(&Th == & Thv ) {
              tv =&T;
              Ptv=Pt;
            }
            else {
              tv= Thv.Find(PP,Ptv,outsidev);
              if( !tv ||  outsidev) {
                if(verbosity>100) cout << " On a pas trouver (u) " << P << " " << endl;
                continue;}
            }
            iut = Thu(tu);
            ivt = Thv(tv);
            if( verbosity>1000) cout << " T " << it  << "  iut " << iut << " ivt " << ivt  <<  endl ;
            FElement3 Ku(Uh[iut]);
            FElementS Kv(Vh[ivt]);
            long n= Kv.NbDoF() ,m=Ku.NbDoF();
            long N= Kv.N, M= Ku.N;
            RNMK_ fv(p,n,N,lastop); //  the value for basic fonction
            RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction     // lastop

            Ku.BF(Dop,Ptu,fu);
            MeshPointStack(stack)->set(Th,P,Pt,T,label);
            if (classoptm) (*Op.optiexpK)(stack);  // call optim version
            if (!same) Kv.BF(Dop,Ptv,fv);
            for ( i=0;  i<n;   i++ ) {
                // attention la fonction test donne la ligne
                //  et la fonction test est en second
                int ig = Kv(i);
                RNM_ wi(fv(i,'.','.'));
                for ( j=0;  j<m;   j++ ) {
                    RNM_ wj(fu(j,'.','.'));
                    int il=0;
                    int jg(Ku(j));
                    if ( !sym ||  ig <= jg )
                        for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++) {
                            // attention la fonction test donne la ligne
                            //  et la fonction test est en second
                            BilinearOperator::K ll(*l);
                            pair<int,int> jj(ll.first.first),ii(ll.first.second);
                            double w_i =  wi(ii.first,ii.second);
                            double w_j =  wj(jj.first,jj.second);
                            R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                            if( verbosity>1000) cout << ig << " " << jg << " "  <<  " " << ccc << " " <<  coef * ccc * w_i*w_j << " T \n"   ;
                            double wij =  w_i*w_j;
                            if (abs(wij)>= 1e-10)
                                A[make_pair(ig,jg)] += coef * ccc * wij;
                        }
                }
            }
        }
    }
    else // int on point ie
        ffassert(0);
 
    *MeshPointStack(stack) = mp;
}


// 3D curve / 3D Surf on meshL
template<class R>
void  AddMatElem(MatriceMap<R> & A,const MeshL & Th,const BilinearOperator & Op,bool sym,int it,int ie,int label,
                 const FESpaceL & Uh,const FESpaceS & Vh,const GQuadratureFormular<R1> & FI,const QuadratureFormular1d & FIb,
                 double *p,void *vstack, bool intmortar=false)
{

    Stack stack=pvoid2Stack(vstack);
    MeshPoint mp= *MeshPointStack(stack);
    R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
    const MeshL & Thu(Uh.Th);
    const MeshS & Thv(Vh.Th);

    bool same = false;
    const EdgeL &T = Th[it];
    long npi;
    long i,j;
    bool classoptm = copt && Op.optiexpK;
    assert(Op.MaxOp() <last_operatortype);

    int lastop=0;
    What_d Dop = Op.DiffOp(lastop);
    if (ie<0) {
        for (npi=0;npi<FI.n;npi++) {
           GQuadraturePoint<R1> pi(FI[npi]);
            double coef = T.mesure()*pi.a;
            R1 Pt(pi),Ptu;
            R2 Ptv;
            R3 P(T(Pt));
            bool outsideu,outsidev;
            // ici trouve le T
            int iut=0,ivt=0;
            const EdgeL * tu;
            const TriangleS *tv;
            if(&Th == & Thu) {
                tu =&T;
                Ptu=Pt;
            }
            else {
                tu= Thu.Find(P,Ptu,outsideu);
                if( !tu ||  outsideu) {
                    if(verbosity>100) cout << " On a pas trouver (u) " << P << " " << endl;
                    continue;
                }
            }
            
           tv= Thv.Find(P,Ptv,outsidev);
           if( !tv || outsidev) {
                if(verbosity>100) cout << " On a pas trouver (v) " << P << " " << endl;
                continue;
            }
            iut = Thu(tu);
            ivt = Thv(tv);
            if( verbosity>1000) cout << " T " << it  << "  iut " << iut << " ivt " << ivt  <<  endl ;
            FElementL Ku(Uh[iut]);
            FElementS Kv(Vh[ivt]);
            long n= Kv.NbDoF() ,m=Ku.NbDoF();
            long N= Kv.N, M= Ku.N;
            RNMK_ fv(p,n,N,lastop); //  the value for basic fonction
            RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction     // lastop

            Ku.BF(Dop,Ptu,fu);
            MeshPointStack(stack)->set(Th,P,Pt,T,label);
            
            if (classoptm) (*Op.optiexpK)(stack);   // call optim version
            if (!same) Kv.BF(Dop,Ptv,fv);
            for ( i=0;  i<n;   i++ ) {
                // attention la fonction test donne la ligne
                //  et la fonction test est en second
                int ig = Kv(i);
                RNM_ wi(fv(i,'.','.'));
                for ( j=0;  j<m;   j++ ) {
                    RNM_ wj(fu(j,'.','.'));
                    int il=0;
                    int jg(Ku(j));
                    if ( !sym ||  ig <= jg )
                        for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++) {
                            // attention la fonction test donne la ligne
                            //  et la fonction test est en second
                            BilinearOperator::K ll(*l);
                            pair<int,int> jj(ll.first.first),ii(ll.first.second);
                            double w_i =  wi(ii.first,ii.second);
                            double w_j =  wj(jj.first,jj.second);
                            R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                            if( verbosity>1000) cout << ig << " " << jg << " "  <<  " " << ccc << " " <<  coef * ccc * w_i*w_j << " T \n"   ;
                            double wij =  w_i*w_j;
                            if (abs(wij)>= 1e-10)
                                A[make_pair(ig,jg)] += coef * ccc * wij;
                        }
                }
            }
        }
    }
    else // int on point ie
        ffassert(0);
 
    *MeshPointStack(stack) = mp;
}

// 3D Surf / 3D curve on meshL
template<class R>
void  AddMatElem(MatriceMap<R> & A,const MeshL & Th,const BilinearOperator & Op,bool sym,int it,int ie,int label,
                 const FESpaceS & Uh,const FESpaceL & Vh,const GQuadratureFormular<R1> & FI,const QuadratureFormular1d & FIb,
                 double *p,void *vstack, bool intmortar=false)
{

    Stack stack=pvoid2Stack(vstack);
    MeshPoint mp= *MeshPointStack(stack);
    R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
    const MeshS & Thu(Uh.Th);
    const MeshL & Thv(Vh.Th);

    bool same = false;
    const EdgeL &T = Th[it];
    long npi;
    long i,j;
    bool classoptm = copt && Op.optiexpK;
    assert(Op.MaxOp() <last_operatortype);

    int lastop=0;
    What_d Dop = Op.DiffOp(lastop);
    double epsP=1e-6; // must be choose
    
    if (ie<0) {
        for (npi=0;npi<FI.n;npi++) {
            GQuadraturePoint<R1> pi(FI[npi]);
            double coef = T.mesure()*pi.a;
            R1 Pt(pi),Ptv;
            R2 Ptu;
            R3 P(T(Pt));
            R2 PP(P.p2());
            bool outsideu,outsidev;
            // ici trouve le T
            int iut=0,ivt=0;
            const TriangleS * tu;
            const EdgeL *tv;
            if(abs(P.z)>epsP) {outsidev=true;tu=0;}
            else tu= Thu.Find(PP,Ptu,outsideu);
              if( !tu ||  outsideu) {
                if(verbosity>100) cout << " On a pas trouver (u) " << P << " " << endl;
                continue;
              }
            if(&Th == & Thv ) {
              tv =&T;
              Ptv=Pt;
            }
            else {
              tv= Thv.Find(P,Ptv,outsideu);
              if( !tv ||  outsideu) {
                if(verbosity>100) cout << " On a pas trouver (u) " << P << " " << endl;
                continue;}
            }
            
            iut = Thu(tu);
            ivt = Thv(tv);
            if( verbosity>1000) cout << " T " << it  << "  iut " << iut << " ivt " << ivt  <<  endl ;
            FElementS Ku(Uh[iut]);
            FElementL Kv(Vh[ivt]);
            long n= Kv.NbDoF() ,m=Ku.NbDoF();
            long N= Kv.N, M= Ku.N;
            RNMK_ fv(p,n,N,lastop); //  the value for basic fonction
            RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction     // lastop

            Ku.BF(Dop,Ptu,fu);
            MeshPointStack(stack)->set(Th,P,Pt,T,label);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            if (!same) Kv.BF(Dop,Ptv,fv);
            for ( i=0;  i<n;   i++ ) {
                // attention la fonction test donne la ligne
                //  et la fonction test est en second
                int ig = Kv(i);
                RNM_ wi(fv(i,'.','.'));
                for ( j=0;  j<m;   j++ ) {
                    RNM_ wj(fu(j,'.','.'));
                    int il=0;
                    int jg(Ku(j));
                    if ( !sym ||  ig <= jg )
                        for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++) {
                            // attention la fonction test donne la ligne
                            //  et la fonction test est en second
                            BilinearOperator::K ll(*l);
                            pair<int,int> jj(ll.first.first),ii(ll.first.second);
                            double w_i =  wi(ii.first,ii.second);
                            double w_j =  wj(jj.first,jj.second);
                            R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                            if( verbosity>1000) cout << ig << " " << jg << " "  <<  " " << ccc << " " <<  coef * ccc * w_i*w_j << " T \n"   ;
                            double wij =  w_i*w_j;
                            if (abs(wij)>= 1e-10)
                                A[make_pair(ig,jg)] += coef * ccc * wij;
                        }
                }
            }
        }
    }
    else // int on point ie
        ffassert(0);
 
    *MeshPointStack(stack) = mp;
}

   ////////////////////////////////////////////////
   // AssembleBilinearForm
   ////////////////////////////////////////////////


    // creating an instance of AssembleBilinearForm with map
    // case 2d
    template<class R>
    void AssembleBilinearForm(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                              MatriceMap<R>  & A, const  FormBilinear * b  )

    {
        /*FH:  case ..in 2D
         in varf ...
         all mesh can can be different ....
         */
        StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
        bool sptrclean=true;
        //     sptr->clean(); // modif FH mars 2006  clean Ptr

        const CDomainOfIntegration & di= *b->di;
        const Mesh * pThdi = GetAny<pmesh>( (* di.Th)(stack));
        SHOWVERB(cout << " FormBilinear () " << endl);
        //MatriceElementaireSymetrique<R> *mates =0;
        // MatriceElementairePleine<R> *matep =0;
        const int useopt=di.UseOpt(stack);
        //double binside=di.binside(stack);
        const bool intmortar=di.intmortar(stack);

        if ( verbosity >1)
        {
            cout << " Integral(1)   on Th "<< &Th << " nv :  " << Th.nv << " nt : " << Th.nt << endl;
            cout << "        Th/ u "<< &Uh.Th << " nv : " << Uh.Th.nv << "   nt : " << Uh.Th.nt << endl;
            cout << "        Th/ v "<< &Vh.Th << " nv : " << Vh.Th.nv << "   nt : " << Vh.Th.nt << endl;
            cout << "        suppose in mortar " << intmortar << "   levelset=  " << di.islevelset() << " withmap: " << di.withmap() << endl;
        }
        Expression  const * const mapt=*di.mapt?di.mapt:0 ;
        Expression  const * const mapu=*di.mapu?di.mapu:0 ;
        bool withmap =di.withmap();
        //   ExecError(" no map  in the case (4) ??");}
        ffassert(pThdi == & Th);
        //const vector<Expression>  & what(di.what);
        CDomainOfIntegration::typeofkind  kind = di.kind;
        set<int> setoflab;
        bool all=true;
        const QuadratureFormular1d & FIE = di.FIE(stack);
        const QuadratureFormular & FITo = di.FIT(stack);
        QuadratureFormular FIT(FITo,3);
        bool VF=b->VF();  // finite Volume or discontinuous Galerkin
        if (verbosity>2) cout << "  -- discontinuous Galerkin  =" << VF << " size of Mat =" << A.size()<< " Bytes\n";
        if (verbosity>3)
        {
            if (CDomainOfIntegration::int1d==kind) cout << "  -- boundary int border ( nQP: "<< FIE.n << ") ,"  ;
            else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all edges ( nQP: "<< FIE.n << "),"  ;
            else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF edges nQP: ("<< FIE.n << ")," ;
            else cout << "  --  int 2d   (nQP: "<< FIT.n << " ) in "  ;
        }
        // if(di.islevelset()) InternalError("Sorry no levelset integration type on this case (1)");
        if(di.islevelset() && (CDomainOfIntegration::int1d!=kind) &&  (CDomainOfIntegration::int2d!=kind) )
        InternalError("Sorry no levelset integration type on no int1d case");

        /*
         if (verbosity>3)
         if (CDomainOfIntegration::int1d==kind) cout << "  -- boundary int border  " ;
         else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all edges, "   ;
         else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF edges, "   ;
         else cout << "  --  int  in  " ; */
        Expandsetoflab(stack,di, setoflab,all);
        /*
         for (size_t i=0;i<what.size();i++)
         {long  lab  = GetAny<long>( (*what[i])(stack));
         setoflab.insert(lab);
         if ( verbosity>3) cout << lab << " ";
         all=false;
         }*/
        if (verbosity>3) cout <<" Optimized = "<< useopt << ", ";
        const E_F0 * poptiexp0=b->b->optiexp0;
        // const E_F0 & optiexpK=*b->b->optiexpK;
        int n_where_in_stack_opt=b->b->where_in_stack_opt.size();
        R** where_in_stack =0;
        if (n_where_in_stack_opt && useopt)
        where_in_stack = new R * [n_where_in_stack_opt];
        if (where_in_stack)
        {
            assert(b->b->v.size()==(size_t) n_where_in_stack_opt);
            for (int i=0;i<n_where_in_stack_opt;i++)
            {
                int offset=b->b->where_in_stack_opt[i];
                assert(offset>10);
                where_in_stack[i]= static_cast<R *>(static_cast<void *>((char*)stack+offset));
                *(where_in_stack[i])=0;
            }


            if(poptiexp0)
            (*poptiexp0)(stack);
            KN<bool> ok(b->b->v.size());
            {  //   remove the zero coef in the liste
                // R zero=R();
                int il=0;
                for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
                ok[il] =  ! (b->b->mesh_indep_stack_opt[il] && ( std::norm(*(where_in_stack[il])) < 1e-100 ) );
            }
            BilinearOperator b_nozer(*b->b,ok);
            if (verbosity % 10 > 3 )
            cout << "   -- nb term in bilinear form  (!0) : " << b_nozer.v.size()
            << "  total " << n_where_in_stack_opt << endl;

            if ( (verbosity/100) % 10 >= 2)
            {
                int il=0;

                for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
                cout << il << " coef (" << l->first << ") = " << *(where_in_stack[il])
                << " offset=" << b->b->where_in_stack_opt[il]
                << " dep mesh " << l->second.MeshIndependent() << b->b->mesh_indep_stack_opt[il] << endl;
            }
        }
        Stack_Ptr<R*>(stack,ElemMatPtrOffset) =where_in_stack;

        KN<double>  p(Vh.esize()+ Uh.esize() );


        if (verbosity >3)
        {
            if (all) cout << " all " << endl ;
            else cout << endl;
        }

        if (di.kind == CDomainOfIntegration::int1d )
        {

            if(di.islevelset())
            {
                double uset = HUGE_VAL;
                R2 Q[2];
                double vol6[2];
                KN<double> phi(Th.nv);phi=uset;
                double f[3], ll=0;
                for(int t=0; t< Th.nt;++t)
                {
                    if ( all || setoflab.find(Th[t].lab) != setoflab.end())
                    {
                        double umx=-HUGE_VAL,umn=HUGE_VAL;
                        for(int i=0;i<3;++i)
                        {
                            int j= Th(t,i);
                            if( phi[j]==uset)
                            {
                                MeshPointStack(stack)->setP(&Th,t,i);
                                phi[j]= di.levelset(stack);//zzzz
                            }
                            f[i]=phi[j];
                            umx = std::max(umx,phi[j]);
                            umn = std::min(umn,phi[j]);

                        }
                        int ntp= IsoLineK(f,Q,1e-10);
                        if(verbosity>999 && ntp==2)
                        {
                            const Triangle &T = Th[t];
                            R2 E(T(Q[0]),T(Q[1]));
                            double le=sqrt((E,E));
                            ll += le;
                            cout << "\t\t" << ntp <<" :  " << Q[0] << " " << Q[1] << " ;  "
                            << f[0] << " " << f[1] << " " << f[2] << "  " << le << " / " << ll<<endl;
                        }
                        if( ntp==2)
                        { if( withmap)
                            AddMatElem(mapu,mapt,A,Th,*b->b,sym,t,10,Th[t].lab,Uh,Vh,FIT,FIE,p,stack,intmortar,Q);
                            else
                            AddMatElem(A,Th,*b->b,sym,t,10,Th[t].lab,Uh,Vh,FIT,FIE,p,stack,intmortar,Q);
                            if(sptrclean) sptrclean=sptr->clean();
                        }
                    }
                }
                FIT =FITo;
            }


            else
            {
                for( int e=0;e<Th.neb;e++)
                {
                    if (all || setoflab.find(Th.bedges[e].lab) != setoflab.end())
                    {
                        int ie,i =Th.BoundaryElement(e,ie);
                        if( withmap)
                        AddMatElem(mapu,mapt,A,Th,*b->b,sym,i,ie,Th.bedges[e].lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
                        else
                        AddMatElem(A,Th,*b->b,sym,i,ie,Th.bedges[e].lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
                        if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                    }
                }
            }}
        else if (di.kind == CDomainOfIntegration::intalledges)
        {
            cerr << " Sorry no implement to hard  "<< endl;
            ExecError("FH: no intalledges on diff mesh ???");
            ffassert(0); // a faire
            if(withmap)
            for (int i=0;i< Th.nt; i++)
            {
                if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                for (int ie=0;ie<3;ie++)
                AddMatElem(mapu,mapt,A,Th,*b->b,sym,i,ie,Th[i].lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
                if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr


            }

            else
            for (int i=0;i< Th.nt; i++)
            {
                if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                for (int ie=0;ie<3;ie++)
                AddMatElem(A,Th,*b->b,sym,i,ie,Th[i].lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
                if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr


            }

        }
        else if (di.kind == CDomainOfIntegration::intallVFedges)
        {

            cerr << " a faire intallVFedges " << endl;
            ffassert(0);

        }
        else if (di.kind == CDomainOfIntegration::int2d )
        {
            // cerr << " a faire CDomainOfIntegration::int2d  " << endl;
            if(di.islevelset())
            {
                double uset = HUGE_VAL;
                R2 Q[2][3];
                double vol6[2];
                KN<double> phi(Th.nv);phi=uset;
                double f[3];
                for(int t=0; t< Th.nt;++t)
                {
                    if ( all || setoflab.find(Th[t].lab) != setoflab.end())
                    {
                        double umx=-HUGE_VAL,umn=HUGE_VAL;
                        for(int i=0;i<3;++i)
                        {
                            int j= Th(t,i);
                            if( phi[j]==uset)
                            {
                                MeshPointStack(stack)->setP(&Th,t,i);
                                phi[j]= di.levelset(stack);//zzzz
                            }
                            f[i]=phi[j];
                            umx = std::max(umx,phi[j]);
                            umn = std::min(umn,phi[j]);

                        }
                        int nt= UnderIso(f,Q, vol6,1e-14);
                        setQF<R2>(FIT,FITo,QuadratureFormular_T_1, Q,vol6,nt);
                        if(FIT.n)
                        {
                            if(withmap)
                            AddMatElem(mapu,mapt,A,Th,*b->b,sym,t,-1,Th[t].lab,Uh,Vh,FIT,FIE,p,stack);
                            else
                            AddMatElem(A,Th,*b->b,sym,t,-1,Th[t].lab,Uh,Vh,FIT,FIE,p,stack);
                        }
                        if(sptrclean) sptrclean=sptr->clean();
                    }
                }
                FIT =FITo;
            }
            else

            {
                if(withmap)
                for (int i=0;i< Th.nt; i++)
                {
                    if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                    AddMatElem(mapu,mapt,A,Th,*b->b,sym,i,-1,Th[i].lab,Uh,Vh,FIT,FIE,p,stack);
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                }

                else
                for (int i=0;i< Th.nt; i++)
                {
                    if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                    AddMatElem(A,Th,*b->b,sym,i,-1,Th[i].lab,Uh,Vh,FIT,FIE,p,stack);
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                }

            }}
        else
        InternalError(" kind of CDomainOfIntegration unknown");

        if (where_in_stack) delete [] where_in_stack;
    }

    // creating an instance of AssembleBilinearForm with map
    // case 3D volume
    template<class R>
    void AssembleBilinearForm(Stack stack,const Mesh3 & Th,const FESpace3 & Uh,const FESpace3 & Vh,bool sym,
                              MatriceMap<R>  & A, const  FormBilinear * b  )

    {
        /*FH:  case ..in 3D
         in varf ...
         all mesh can can be different ....
         */

        StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
        bool sptrclean=true;
        //     sptr->clean(); // modif FH mars 2006  clean Ptr
        Fem2D::MeshPoint & mp (*Fem2D::MeshPointStack(stack)), mps = mp;

        const CDomainOfIntegration & di= *b->di;
        const Mesh3 * pThdi = GetAny<pmesh3>( (* di.Th)(stack));
        SHOWVERB(cout << " FormBilinear () " << endl);
        //MatriceElementaireSymetrique<R> *mates =0;
        // MatriceElementairePleine<R> *matep =0;
        const int useopt=di.UseOpt(stack);
        //double binside=di.binside(stack);
        const bool intmortar=di.intmortar(stack);
        if ( verbosity >1)
        {
            cout << " Integral(2)   on Th "<< &Th << " nv :  " << Th.nv << " nt : " << Th.nt << endl;
            cout << "        Th/ u "<< &Uh.Th << " nv : " << Uh.Th.nv << "   nt : " << Uh.Th.nt << endl;
            cout << "        Th/ v "<< &Vh.Th << " nv : " << Vh.Th.nv << "   nt : " << Vh.Th.nt << endl;
            cout << "        suppose in mortar " << intmortar << endl;
        }
        assert(pThdi == & Th);
        //const vector<Expression>  & what(di.what);
        CDomainOfIntegration::typeofkind  kind = di.kind;

        set<int> setoflab;
        bool all=true;
        // const QuadratureFormular1d & FIEo = di.FIE(stack);
        const QuadratureFormular & FITo = di.FIT(stack);
        const GQuadratureFormular<R3> & FIVo = di.FIV(stack);
        //  to change the quadrature on element ... may 2014 FH ..
        // QuadratureFormular1d  FIE(FIEo,3);
        QuadratureFormular FIT(FITo,3);
        GQuadratureFormular<R3>  FIV(FIVo,3);



        //    const QuadratureFormular & FIT = di.FIT(stack);
        //    const Fem2D::GQuadratureFormular<R3> & FIV = di.FIV(stack);
        bool VF=b->VF();  // finite Volume or discontinuous Galerkin
        if (verbosity>2) cout << "  -- discontinuous Galerkin  =" << VF << " size of Mat =" << A.size()<< " Bytes\n";
        if (verbosity>3)
        {
            if (CDomainOfIntegration::int2d==kind) cout << "  -- boundary int border ( nQP: "<< FIT.n << ") ,"  ;
            else  if (CDomainOfIntegration::intallfaces==kind) cout << "  -- boundary int all edges ( nQP: "<< FIT.n << "),"  ;
            //else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF edges nQP: ("<< FIT.n << ")," ;
            else cout << "  --  int 3d   (nQP: "<< FIV.n << " ) in "  ;
        }
        if(di.islevelset()) InternalError("Sorry no levelset integration type on this case (2)");
        if(di.islevelset() && (CDomainOfIntegration::int2d!=kind) && (CDomainOfIntegration::int3d!=kind) ) InternalError("Sorry no levelset integration type on no int2d case");

        Expandsetoflab(stack,di, setoflab,all);
        /*
         for (size_t i=0;i<what.size();i++)
         {long  lab  = GetAny<long>( (*what[i])(stack));
         setoflab.insert(lab);
         if ( verbosity>3) cout << lab << " ";
         all=false;
         }*/
        if (verbosity>3) cout <<" Optimized = "<< useopt << ", ";
        const E_F0 *poptiexp0=b->b->optiexp0;
        // const E_F0 & optiexpK=*b->b->optiexpK;
        int n_where_in_stack_opt=b->b->where_in_stack_opt.size();
        R** where_in_stack =0;
        if (n_where_in_stack_opt && useopt)
        where_in_stack = new R * [n_where_in_stack_opt];
        if (where_in_stack)
        {
            assert(b->b->v.size()==(size_t) n_where_in_stack_opt);
            for (int i=0;i<n_where_in_stack_opt;i++)
            {
                int offset=b->b->where_in_stack_opt[i];
                assert(offset>10);
                where_in_stack[i]= static_cast<R *>(static_cast<void *>((char*)stack+offset));
                *(where_in_stack[i])=0;
            }


            if(poptiexp0)
            (*poptiexp0)(stack);
            KN<bool> ok(b->b->v.size());
            {  //   remove the zero coef in the liste
                // R zero=R();
                int il=0;
                for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
                ok[il] =  ! (b->b->mesh_indep_stack_opt[il] && ( std::norm(*(where_in_stack[il])) < 1e-100 ) );
            }
            BilinearOperator b_nozer(*b->b,ok);
            if (verbosity % 10 > 3 )
            cout << "   -- nb term in bilinear form  (!0) : " << b_nozer.v.size()
            << "  total " << n_where_in_stack_opt << endl;

            if ( (verbosity/100) % 10 >= 2)
            {
                int il=0;

                for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
                cout << il << " coef (" << l->first << ") = " << *(where_in_stack[il])
                << " offset=" << b->b->where_in_stack_opt[il]
                << " dep mesh " << l->second.MeshIndependent() << b->b->mesh_indep_stack_opt[il] << endl;
            }
        }
        Stack_Ptr<R*>(stack,ElemMatPtrOffset) =where_in_stack;

        KN<double>  p(Vh.esize()+ Uh.esize() );


        if (verbosity >3)
        {
            if (all) cout << " all " << endl ;
            else cout << endl;
        }

        if (di.kind == CDomainOfIntegration::int2d )
        {
            for( int e=0;e<Th.nbe;e++)
            {
                if (all || setoflab.find(Th.be(e).lab) != setoflab.end())
                {
                    int ie,i =Th.BoundaryElement(e,ie);
                    AddMatElem(A,Th,*b->b,sym,i,ie,Th.be(e).lab,Uh,Vh,FIV,FIT,p,stack,intmortar);
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                }
            }
        }
        else if (di.kind == CDomainOfIntegration::intallfaces)
        {
            ffassert(0); // a faire

            for (int i=0;i< Th.nt; i++)
            {
                if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                for (int ie=0;ie<3;ie++)
                AddMatElem(A,Th,*b->b,sym,i,ie,Th[i].lab,Uh,Vh,FIV,FIT,p,stack,intmortar);

                if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr


            }

        }
        /* else if (di.kind == CDomainOfIntegration::intallVFedges)
         {

         cerr << " a faire intallVFedges " << endl;
         ffassert(0);

         }  */
        else if (di.kind == CDomainOfIntegration::int3d )
        {
            if(di.islevelset())  //  may 2014 FH ...
            {   // int3d levelset < 0
                double llevelset = 0;
                const double uset = std::numeric_limits<double>::max();
                // cout << " uset ="<<uset << endl;
                R3 Q[3][4];
                double vol6[3];
                KN<double> phi(Th.nv);
                phi=uset;
                double f[4];

                for (int t=0;t< Th.nt; t++)
                {

                    const Mesh3::Element & K(Th[t]);
                    if (all || setoflab.find(Th[t].lab) != setoflab.end())

                    {
                        double umx=std::numeric_limits<double>::lowest(),umn=std::numeric_limits<double>::max();
                        for(int i=0;i<4;++i)
                        {
                            int j= Th(t,i);
                            if( phi[j]==uset)
                            {
                                MeshPointStack(stack)->setP(&Th,t,i);
                                phi[j]= di.levelset(stack);//zzzz
                            }
                            f[i]=phi[j];
                        }
                        int ntets= UnderIso(f,Q, vol6,1e-14);
                        setQF<R3>(FIV,FIVo,QuadratureFormular_Tet_1, Q,vol6,ntets);
                        if(FIV.n)
                        {
                            AddMatElem(A,Th,*b->b,sym,t,-1,Th[t].lab,Uh,Vh,FIV,FIT,p,stack);
                            if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr

                        }

                    }
                }
                FIV = FIVo;

            }
            else

            {
                // cerr << " a faire CDomainOfIntegration::int3d  " << endl;
                for (int i=0;i< Th.nt; i++)
                {
                    if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                    AddMatElem(A,Th,*b->b,sym,i,-1,Th[i].lab,Uh,Vh,FIV,FIT,p,stack);
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                }

            } }
        else
        InternalError(" kind of CDomainOfIntegration unknown");

        if (where_in_stack) delete [] where_in_stack;
        mp=mps;// restore x,y,z

    }

    // creating an instance of AssembleBilinearForm with map
    // case 3D surface
    template<class R>
    void AssembleBilinearForm(Stack stack,const MeshS & Th,const FESpaceS & Uh,const FESpaceS & Vh,bool sym,
                              MatriceMap<R>  & A, const  FormBilinear * b  )

    {
        StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
        bool sptrclean=true;
        //     sptr->clean(); // modif FH mars 2006  clean Ptr

        const CDomainOfIntegration & di= *b->di;
        ///typedef typename Trait_MESHO<FESpaceS>::MeshO * pmeshO;
       // pmeshO  ThbfO = GetAny<pmeshO>((*b->di->Th)(stack)); // case 3D surface ThbfO =
        pmeshS  pThdi = GetAny<pmeshS>((*b->di->Th)(stack)); //Trait_MESHO<FESpaceS>::topmesh(ThbfO);  //

        SHOWVERB(cout << " FormBilinear () " << endl);
        //MatriceElementaireSymetrique<R> *mates =0;
        // MatriceElementairePleine<R> *matep =0;
        const int useopt=di.UseOpt(stack);
        //double binside=di.binside(stack);
        const bool intmortar=di.intmortar(stack);
        if ( verbosity >1)
        {
            cout << " Integral(3)   on Th "<< &Th << " nv :  " << Th.nv << " nt : " << Th.nt << endl;
            cout << "        Th/ u "<< &Uh.Th << " nv : " << Uh.Th.nv << "   nt : " << Uh.Th.nt << endl;
            cout << "        Th/ v "<< &Vh.Th << " nv : " << Vh.Th.nv << "   nt : " << Vh.Th.nt << endl;
            cout << "        suppose in mortar " << intmortar << "   levelset=  " << di.islevelset() << " withmap: " << di.withmap() << endl;
        }
        Expression  const * const mapt=*di.mapt?di.mapt:0 ;
        Expression  const * const mapu=*di.mapu?di.mapu:0 ;
        bool withmap =di.withmap();
        //   ExecError(" no map  in the case (4) ??");}
        assert(pThdi == & Th);
        //const vector<Expression>  & what(di.what);
        CDomainOfIntegration::typeofkind  kind = di.kind;
        set<int> setoflab;
        bool all=true;
        const QuadratureFormular1d & FIE = di.FIE(stack);
        const QuadratureFormular & FITo = di.FIT(stack);
        QuadratureFormular FIT(FITo,3);
        bool VF=b->VF();  // finite Volume or discontinuous Galerkin
        if (verbosity>2) cout << "  -- discontinuous Galerkin  =" << VF << " size of Mat =" << A.size()<< " Bytes\n";
        if (verbosity>3)
        {
            if (CDomainOfIntegration::int1d==kind) cout << "  -- boundary int border ( nQP: "<< FIE.n << ") ,"  ;
            else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all edges ( nQP: "<< FIE.n << "),"  ;
            else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF edges nQP: ("<< FIE.n << ")," ;
            else cout << "  --  int 2d   (nQP: "<< FIT.n << " ) in "  ;
        }
        // if(di.islevelset()) InternalError("Sorry no levelset integration type on this case (1)");
        if(di.islevelset() && (CDomainOfIntegration::int1d!=kind) &&  (CDomainOfIntegration::int2d!=kind) )
            InternalError("Sorry no levelset integration type on no int1d case");

        /*
         if (verbosity>3)
         if (CDomainOfIntegration::int1d==kind) cout << "  -- boundary int border  " ;
         else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all edges, "   ;
         else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF edges, "   ;
         else cout << "  --  int  in  " ; */
        Expandsetoflab(stack,di, setoflab,all);
        /*
         for (size_t i=0;i<what.size();i++)
         {long  lab  = GetAny<long>( (*what[i])(stack));
         setoflab.insert(lab);
         if ( verbosity>3) cout << lab << " ";
         all=false;
         }*/
        if (verbosity>3) cout <<" Optimized = "<< useopt << ", ";
        const E_F0 * poptiexp0=b->b->optiexp0;
        // const E_F0 & optiexpK=*b->b->optiexpK;
        int n_where_in_stack_opt=b->b->where_in_stack_opt.size();
        R** where_in_stack =0;
        if (n_where_in_stack_opt && useopt)
            where_in_stack = new R * [n_where_in_stack_opt];
        if (where_in_stack)
        {
            assert(b->b->v.size()==(size_t) n_where_in_stack_opt);
            for (int i=0;i<n_where_in_stack_opt;i++)
            {
                int offset=b->b->where_in_stack_opt[i];
                assert(offset>10);
                where_in_stack[i]= static_cast<R *>(static_cast<void *>((char*)stack+offset));
                *(where_in_stack[i])=0;
            }


            if(poptiexp0)
                (*poptiexp0)(stack);
            KN<bool> ok(b->b->v.size());
            {  //   remove the zero coef in the liste
                // R zero=R();
                int il=0;
                for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
                    ok[il] =  ! (b->b->mesh_indep_stack_opt[il] && ( std::norm(*(where_in_stack[il])) < 1e-100 ) );
            }
            BilinearOperator b_nozer(*b->b,ok);
            if (verbosity % 10 > 3 )
                cout << "   -- nb term in bilinear form  (!0) : " << b_nozer.v.size()
                << "  total " << n_where_in_stack_opt << endl;

            if ( (verbosity/100) % 10 >= 2)
            {
                int il=0;

                for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
                    cout << il << " coef (" << l->first << ") = " << *(where_in_stack[il])
                    << " offset=" << b->b->where_in_stack_opt[il]
                    << " dep mesh " << l->second.MeshIndependent() << b->b->mesh_indep_stack_opt[il] << endl;
            }
        }
        Stack_Ptr<R*>(stack,ElemMatPtrOffset) =where_in_stack;

        KN<double>  p(Vh.esize()+ Uh.esize() );


        if (verbosity >3)
        {
            if (all) cout << " all " << endl ;
            else cout << endl;
        }

        if (di.kind == CDomainOfIntegration::int1d )
        {

            if(di.islevelset())
            {
                double uset = HUGE_VAL;
                R2 Q[2];
                double vol6[2];
                KN<double> phi(Th.nv);phi=uset;
                double f[3], ll=0;
                for(int t=0; t< Th.nt;++t)
                {
                    if ( all || setoflab.find(Th[t].lab) != setoflab.end())
                    {
                        double umx=-HUGE_VAL,umn=HUGE_VAL;
                        for(int i=0;i<3;++i)
                        {
                            int j= Th(t,i);
                            if( phi[j]==uset)
                            {
                                MeshPointStack(stack)->setP(&Th,t,i);
                                phi[j]= di.levelset(stack);//zzzz
                            }
                            f[i]=phi[j];
                            umx = std::max(umx,phi[j]);
                            umn = std::min(umn,phi[j]);

                        }
                        int ntp= IsoLineK(f,Q,1e-10);
                        if(verbosity>999 && ntp==2)
                        {
                            const TriangleS &T = Th[t];
                            R3 E(T(Q[0]),T(Q[1]));
                            double le=sqrt((E,E));
                            ll += le;
                            cout << "\t\t" << ntp <<" :  " << Q[0] << " " << Q[1] << " ;  "
                            << f[0] << " " << f[1] << " " << f[2] << "  " << le << " / " << ll<<endl;
                        }
                        if( ntp==2)
                        { //if( withmap)
                           // AddMatElem(mapu,mapt,A,Th,*b->b,sym,t,10,Th[t].lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
                        //else
                            AddMatElem(A,Th,*b->b,sym,t,10,Th[t].lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
                            if(sptrclean) sptrclean=sptr->clean();
                        }
                    }
                }
                FIT =FITo;
            }


            else
            {
                for( int e=0;e<Th.nbe;e++)
                {
                    if (all || setoflab.find(Th.be(e).lab) != setoflab.end())
                    {
                        int ie,i =Th.BoundaryElement(e,ie);
                        //if( withmap)
                        //    AddMatElem(mapu,mapt,A,Th,*b->b,sym,i,ie,Th.be(e).lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
                        //else
                            AddMatElem(A,Th,*b->b,sym,i,ie,Th.be(e).lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
                        if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                    }
                }
            }}
        /*else if (di.kind == CDomainOfIntegration::intalledges)
        {
            cerr << " Sorry no implement to hard  "<< endl;
            ExecError("FH: no intalledges on diff mesh ???");
            ffassert(0); // a faire
            if(withmap)
                for (int i=0;i< Th.nt; i++)
                {
                    if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                        for (int ie=0;ie<3;ie++)
                            AddMatElem(mapu,mapt,A,Th,*b->b,sym,i,ie,Th[i].lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr


                }

            else
                for (int i=0;i< Th.nt; i++)
                {
                    if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                        for (int ie=0;ie<3;ie++)
                            AddMatElem(A,Th,*b->b,sym,i,ie,Th[i].lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr


                }

        }*/
        else if (di.kind == CDomainOfIntegration::intallVFedges)
        {

            cerr << " a faire intallVFedges " << endl;
            ffassert(0);

        }
        else if (di.kind == CDomainOfIntegration::int2d ) {
            // cerr << " a faire CDomainOfIntegration::int2d  " << endl;
            if(di.islevelset())
            {
                double uset = HUGE_VAL;
                R2 Q[2][3];
                double vol6[2];
                KN<double> phi(Th.nv);phi=uset;
                double f[3];
                for(int t=0; t< Th.nt;++t)
                {
                    if ( all || setoflab.find(Th[t].lab) != setoflab.end())
                    {
                        double umx=-HUGE_VAL,umn=HUGE_VAL;
                        for(int i=0;i<3;++i)
                        {
                            int j= Th(t,i);
                            if( phi[j]==uset)
                            {
                                MeshPointStack(stack)->setP(&Th,t,i);
                                phi[j]= di.levelset(stack);//zzzz
                            }
                            f[i]=phi[j];
                            umx = std::max(umx,phi[j]);
                            umn = std::min(umn,phi[j]);

                        }
                        int nt= UnderIso(f,Q, vol6,1e-14);
                        setQF<R2>(FIT,FITo,QuadratureFormular_T_1, Q,vol6,nt);
                        if(FIT.n)
                        {
                            //if(withmap)
                             //   AddMatElem(mapu,mapt,A,Th,*b->b,sym,t,-1,Th[t].lab,Uh,Vh,FIT,FIE,p,stack);
                            //else
                                AddMatElem(A,Th,*b->b,sym,t,-1,Th[t].lab,Uh,Vh,FIT,FIE,p,stack);
                        }
                        if(sptrclean) sptrclean=sptr->clean();
                    }
                }
                FIT =FITo;
            }
            else

            {
                //if(withmap)
                 //   for (int i=0;i< Th.nt; i++)
                  //  {
                   //     if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                    //        AddMatElem(mapu,mapt,A,Th,*b->b,sym,i,-1,Th[i].lab,Uh,Vh,FIT,FIE,p,stack);
                    //    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                   // }

                //else
                    for (int i=0;i< Th.nt; i++)
                    {
                        if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                            AddMatElem(A,Th,*b->b,sym,i,-1,Th[i].lab,Uh,Vh,FIT,FIE,p,stack);
                        if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                    }

            }

        }
        else
            InternalError(" kind of CDomainOfIntegration unknown");

        if (where_in_stack) delete [] where_in_stack;
            }



    // creating an instance of AssembleBilinearForm with map
    // case 3D curve
    template<class R>
    void AssembleBilinearForm(Stack stack,const MeshL & Th,const FESpaceL & Uh,const FESpaceL & Vh,bool sym,
                              MatriceMap<R>  & A, const  FormBilinear * b  )

    {
        StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
        bool sptrclean=true;
        //     sptr->clean(); // modif FH mars 2006  clean Ptr

        const CDomainOfIntegration & di= *b->di;
        pmeshL  pThdi = GetAny<pmeshL>((*b->di->Th)(stack));

        SHOWVERB(cout << " FormBilinear () " << endl);

        const int useopt=di.UseOpt(stack);
        //double binside=di.binside(stack);
        const bool intmortar=di.intmortar(stack);
        if ( verbosity >1)
        {
            cout << " Integral(4)   on Th "<< &Th << " nv :  " << Th.nv << " nt : " << Th.nt << endl;
            cout << "        Th/ u "<< &Uh.Th << " nv : " << Uh.Th.nv << "   nt : " << Uh.Th.nt << endl;
            cout << "        Th/ v "<< &Vh.Th << " nv : " << Vh.Th.nv << "   nt : " << Vh.Th.nt << endl;
            cout << "        suppose in mortar " << intmortar << "   levelset=  " << di.islevelset() << " withmap: " << di.withmap() << endl;
        }
        Expression  const * const mapt=*di.mapt?di.mapt:0 ;
        Expression  const * const mapu=*di.mapu?di.mapu:0 ;
        bool withmap =di.withmap();
        assert(pThdi == & Th);
        CDomainOfIntegration::typeofkind  kind = di.kind;
        set<int> setoflab;
        bool all=true;

        const GQuadratureFormular<R1> & FITo = di.FIE(stack);
        GQuadratureFormular<R1> FIT(FITo,3);

        bool VF=b->VF();  // finite Volume or discontinuous Galerkin
        if (verbosity>2) cout << "  -- discontinuous Galerkin  =" << VF << " size of Mat =" << A.size()<< " Bytes\n";

        // if(di.islevelset()) InternalError("Sorry no levelset integration type on this case (1)");
        if(di.islevelset() && (CDomainOfIntegration::int1d!=kind))
            InternalError("Sorry no levelset integration type on no int1d case");

        Expandsetoflab(stack,di, setoflab,all);

        if (verbosity>3) cout <<" Optimized = "<< useopt << ", ";
        const E_F0 * poptiexp0=b->b->optiexp0;
        // const E_F0 & optiexpK=*b->b->optiexpK;
        int n_where_in_stack_opt=b->b->where_in_stack_opt.size();
        R** where_in_stack =0;
        if (n_where_in_stack_opt && useopt)
            where_in_stack = new R * [n_where_in_stack_opt];
        if (where_in_stack) {
            assert(b->b->v.size()==(size_t) n_where_in_stack_opt);
            for (int i=0;i<n_where_in_stack_opt;i++) {
                int offset=b->b->where_in_stack_opt[i];
                assert(offset>10);
                where_in_stack[i]= static_cast<R *>(static_cast<void *>((char*)stack+offset));
                *(where_in_stack[i])=0;
            }


            if(poptiexp0)
                (*poptiexp0)(stack);
            KN<bool> ok(b->b->v.size());
            {  //   remove the zero coef in the liste
                // R zero=R();
                int il=0;
                for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
                    ok[il] =  ! (b->b->mesh_indep_stack_opt[il] && ( std::norm(*(where_in_stack[il])) < 1e-100 ) );
            }
            BilinearOperator b_nozer(*b->b,ok);
            if (verbosity % 10 > 3 )
                cout << "   -- nb term in bilinear form  (!0) : " << b_nozer.v.size()
                << "  total " << n_where_in_stack_opt << endl;

            if ( (verbosity/100) % 10 >= 2) {
                int il=0;

                for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
                    cout << il << " coef (" << l->first << ") = " << *(where_in_stack[il])
                    << " offset=" << b->b->where_in_stack_opt[il]
                    << " dep mesh " << l->second.MeshIndependent() << b->b->mesh_indep_stack_opt[il] << endl;
            }
        }
        Stack_Ptr<R*>(stack,ElemMatPtrOffset) =where_in_stack;

        KN<double>  p(Vh.esize()+ Uh.esize() );


        if (verbosity >3) {
            if (all) cout << " all " << endl ;
            else cout << endl;
        }

        if (di.kind == CDomainOfIntegration::int1d ) {
            if(di.islevelset())   ////// must be check
                ffassert(0);
            else {
                 for (int i=0;i< Th.nt; i++) {
                    if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                        AddMatElem(A,Th,*b->b,sym,i,-1,Th[i].lab,Uh,Vh,FIT,0,p,stack);
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                }
              }
        }
        else { cout << " di.kind " << di.kind << endl;
            InternalError(" kind of CDomainOfIntegration unknown");
        }
        if (where_in_stack) delete [] where_in_stack;
    }

   // creating an instance of AssembleBilinearForm with map
   // case 3D curve / 2D on meshL
   template<class R>
   void AssembleBilinearForm(Stack stack,const MeshL & Th,const FESpaceL & Uh,const FESpace & Vh,bool sym,
                             MatriceMap<R>  & A, const  FormBilinear * b  )

   {
       StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
       bool sptrclean=true;

       const CDomainOfIntegration & di= *b->di;
       pmeshL  pThdi = GetAny<pmeshL>((*b->di->Th)(stack));

       SHOWVERB(cout << " FormBilinear () " << endl);

       const int useopt=di.UseOpt(stack);
       //double binside=di.binside(stack);
       const bool intmortar=di.intmortar(stack);
       if ( verbosity >1)
       {
           cout << " Integral(5)   on Th "<< &Th << " nv :  " << Th.nv << " nt : " << Th.nt << endl;
           cout << "        Th/ u "<< &Uh.Th << " nv : " << Uh.Th.nv << "   nt : " << Uh.Th.nt << endl;
           cout << "        Th/ v "<< &Vh.Th << " nv : " << Vh.Th.nv << "   nt : " << Vh.Th.nt << endl;
           cout << "        suppose in mortar " << intmortar << "   levelset=  " << di.islevelset() << " withmap: " << di.withmap() << endl;
       }
       Expression  const * const mapt=*di.mapt?di.mapt:0 ;
       Expression  const * const mapu=*di.mapu?di.mapu:0 ;
       bool withmap =di.withmap();
       assert(pThdi == & Th);
       CDomainOfIntegration::typeofkind  kind = di.kind;
       set<int> setoflab;
       bool all=true;

       const GQuadratureFormular<R1> & FITo = di.FIE(stack);
       GQuadratureFormular<R1> FIT(FITo,3);

       bool VF=b->VF();  // finite Volume or discontinuous Galerkin
       if (verbosity>2) cout << "  -- discontinuous Galerkin  =" << VF << " size of Mat =" << A.size()<< " Bytes\n";

       // if(di.islevelset()) InternalError("Sorry no levelset integration type on this case (1)");
       if(di.islevelset() && (CDomainOfIntegration::int1d!=kind))
           InternalError("Sorry no levelset integration type on no int1d case");

       Expandsetoflab(stack,di, setoflab,all);

       if (verbosity>3) cout <<" Optimized = "<< useopt << ", ";
       const E_F0 * poptiexp0=b->b->optiexp0;
       // const E_F0 & optiexpK=*b->b->optiexpK;
       int n_where_in_stack_opt=b->b->where_in_stack_opt.size();
       R** where_in_stack =0;
       if (n_where_in_stack_opt && useopt)
           where_in_stack = new R * [n_where_in_stack_opt];
       if (where_in_stack) {
           assert(b->b->v.size()==(size_t) n_where_in_stack_opt);
           for (int i=0;i<n_where_in_stack_opt;i++) {
               int offset=b->b->where_in_stack_opt[i];
               assert(offset>10);
               where_in_stack[i]= static_cast<R *>(static_cast<void *>((char*)stack+offset));
               *(where_in_stack[i])=0;
           }


           if(poptiexp0)
               (*poptiexp0)(stack);
           KN<bool> ok(b->b->v.size());
           int il=0;
           for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
             ok[il] =  ! (b->b->mesh_indep_stack_opt[il] && ( std::norm(*(where_in_stack[il])) < 1e-100 ) );
           
           BilinearOperator b_nozer(*b->b,ok);
           if (verbosity % 10 > 3 )
               cout << "   -- nb term in bilinear form  (!0) : " << b_nozer.v.size()
               << "  total " << n_where_in_stack_opt << endl;

           if ( (verbosity/100) % 10 >= 2) {
               int il=0;

               for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
                   cout << il << " coef (" << l->first << ") = " << *(where_in_stack[il])
                   << " offset=" << b->b->where_in_stack_opt[il]
                   << " dep mesh " << l->second.MeshIndependent() << b->b->mesh_indep_stack_opt[il] << endl;
           }
       }
       Stack_Ptr<R*>(stack,ElemMatPtrOffset) =where_in_stack;

       KN<double>  p(Vh.esize()+ Uh.esize() );


       if (verbosity >3) {
           if (all) cout << " all " << endl ;
           else cout << endl;
       }

       if (di.kind == CDomainOfIntegration::int1d ) {
           if(di.islevelset())
               ffassert(0);
           else {
                for (int i=0;i< Th.nt; i++) {
                   if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                       AddMatElem(A,Th,*b->b,sym,i,-1,Th[i].lab,Uh,Vh,FIT,0,p,stack);
                   if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
               }
             }
       }
       else { cout << " di.kind " << di.kind << endl;
           InternalError(" kind of CDomainOfIntegration unknown");
       }
       if (where_in_stack) delete [] where_in_stack;
   }

   // case 2D / 3D curve on meshL
   template<class R>
   void AssembleBilinearForm(Stack stack,const MeshL & Th,const FESpace & Uh,const FESpaceL & Vh,bool sym,
                          MatriceMap<R>  & A, const  FormBilinear * b  )

   {
       StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
       bool sptrclean=true;

       const CDomainOfIntegration & di= *b->di;
       pmeshL  pThdi = GetAny<pmeshL>((*b->di->Th)(stack));

       SHOWVERB(cout << " FormBilinear () " << endl);

       const int useopt=di.UseOpt(stack);
       //double binside=di.binside(stack);
       const bool intmortar=di.intmortar(stack);
       if ( verbosity >1) {
        cout << " Integral(6)   on Th "<< &Th << " nv :  " << Th.nv << " nt : " << Th.nt << endl;
        cout << "        Th/ u "<< &Uh.Th << " nv : " << Uh.Th.nv << "   nt : " << Uh.Th.nt << endl;
        cout << "        Th/ v "<< &Vh.Th << " nv : " << Vh.Th.nv << "   nt : " << Vh.Th.nt << endl;
        cout << "        suppose in mortar " << intmortar << "   levelset=  " << di.islevelset() << " withmap: " << di.withmap() << endl;
       }
       Expression  const * const mapt=*di.mapt?di.mapt:0 ;
       Expression  const * const mapu=*di.mapu?di.mapu:0 ;
       bool withmap =di.withmap();
       assert(pThdi == & Th);
       CDomainOfIntegration::typeofkind  kind = di.kind;
       set<int> setoflab;
       bool all=true;

       const GQuadratureFormular<R1> & FITo = di.FIE(stack);
       GQuadratureFormular<R1> FIT(FITo,3);

       bool VF=b->VF();  // finite Volume or discontinuous Galerkin
       if (verbosity>2) cout << "  -- discontinuous Galerkin  =" << VF << " size of Mat =" << A.size()<< " Bytes\n";

       // if(di.islevelset()) InternalError("Sorry no levelset integration type on this case (1)");
       if(di.islevelset() && (CDomainOfIntegration::int1d!=kind))
         InternalError("Sorry no levelset integration type on no int1d case");

       Expandsetoflab(stack,di, setoflab,all);

       if (verbosity>3) cout <<" Optimized = "<< useopt << ", ";
       const E_F0 * poptiexp0=b->b->optiexp0;
    // const E_F0 & optiexpK=*b->b->optiexpK;
       int n_where_in_stack_opt=b->b->where_in_stack_opt.size();
       R** where_in_stack =0;
       if (n_where_in_stack_opt && useopt)
         where_in_stack = new R * [n_where_in_stack_opt];
       if (where_in_stack) {
         assert(b->b->v.size()==(size_t) n_where_in_stack_opt);
         for (int i=0;i<n_where_in_stack_opt;i++) {
           int offset=b->b->where_in_stack_opt[i];
           assert(offset>10);
           where_in_stack[i]= static_cast<R *>(static_cast<void *>((char*)stack+offset));
           *(where_in_stack[i])=0;
         }
         if(poptiexp0)
           (*poptiexp0)(stack);
         KN<bool> ok(b->b->v.size());
         int il=0;
         for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
           ok[il] =  ! (b->b->mesh_indep_stack_opt[il] && ( std::norm(*(where_in_stack[il])) < 1e-100 ) );
        
         BilinearOperator b_nozer(*b->b,ok);
         if (verbosity % 10 > 3 ) cout << "   -- nb term in bilinear form  (!0) : " << b_nozer.v.size() << "  total " << n_where_in_stack_opt << endl;

         if ( (verbosity/100) % 10 >= 2) {
           int il=0;
           for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
              cout << il << " coef (" << l->first << ") = " << *(where_in_stack[il]) << " offset=" << b->b->where_in_stack_opt[il]
              << " dep mesh " << l->second.MeshIndependent() << b->b->mesh_indep_stack_opt[il] << endl;
         }
       }
       Stack_Ptr<R*>(stack,ElemMatPtrOffset) =where_in_stack;
       KN<double>  p(Vh.esize()+ Uh.esize() );

       if (verbosity >3) {
         if (all) cout << " all " << endl ;
         else cout << endl;
       }

       if (di.kind == CDomainOfIntegration::int1d ) {
         if(di.islevelset())
           ffassert(0);
         else {
           for (int i=0;i< Th.nt; i++) {
             if ( all || setoflab.find(Th[i].lab) != setoflab.end())
               AddMatElem(A,Th,*b->b,sym,i,-1,Th[i].lab,Uh,Vh,FIT,0,p,stack);
             if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
            }
          }
       }
       else { cout << " di.kind " << di.kind << endl;
         InternalError(" kind of CDomainOfIntegration unknown");
       }
       if (where_in_stack) delete [] where_in_stack;
   }

// case 3D Surf / 3D volume on meshS
 template<class R>
 void AssembleBilinearForm(Stack stack,const MeshS & Th,const FESpaceS & Uh,const FESpace3 & Vh,bool sym,
                           MatriceMap<R>  & A, const  FormBilinear * b  )

{
 StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
 bool sptrclean=true;
 //     sptr->clean(); // modif FH mars 2006  clean Ptr

 const CDomainOfIntegration & di= *b->di;
 ///typedef typename Trait_MESHO<FESpaceS>::MeshO * pmeshO;
// pmeshO  ThbfO = GetAny<pmeshO>((*b->di->Th)(stack)); // case 3D surface ThbfO =
 pmeshS  pThdi = GetAny<pmeshS>((*b->di->Th)(stack)); //Trait_MESHO<FESpaceS>::topmesh(ThbfO);  //

 SHOWVERB(cout << " FormBilinear () " << endl);
 //MatriceElementaireSymetrique<R> *mates =0;
 // MatriceElementairePleine<R> *matep =0;
 const int useopt=di.UseOpt(stack);
 //double binside=di.binside(stack);
 const bool intmortar=di.intmortar(stack);
 if ( verbosity >1)
 {
     cout << " Integral(7)   on Th "<< &Th << " nv :  " << Th.nv << " nt : " << Th.nt << endl;
     cout << "        Th/ u "<< &Uh.Th << " nv : " << Uh.Th.nv << "   nt : " << Uh.Th.nt << endl;
     cout << "        Th/ v "<< &Vh.Th << " nv : " << Vh.Th.nv << "   nt : " << Vh.Th.nt << endl;
     cout << "        suppose in mortar " << intmortar << "   levelset=  " << di.islevelset() << " withmap: " << di.withmap() << endl;
 }
 Expression  const * const mapt=*di.mapt?di.mapt:0 ;
 Expression  const * const mapu=*di.mapu?di.mapu:0 ;
 bool withmap =di.withmap();
 //   ExecError(" no map  in the case (4) ??");}
 assert(pThdi == & Th);
 //const vector<Expression>  & what(di.what);
 CDomainOfIntegration::typeofkind  kind = di.kind;
 set<int> setoflab;
 bool all=true;
 const QuadratureFormular1d & FIE = di.FIE(stack);
 const QuadratureFormular & FITo = di.FIT(stack);
 QuadratureFormular FIT(FITo,3);
 bool VF=b->VF();  // finite Volume or discontinuous Galerkin
 if (verbosity>2) cout << "  -- discontinuous Galerkin  =" << VF << " size of Mat =" << A.size()<< " Bytes\n";
 if (verbosity>3)
 {
     if (CDomainOfIntegration::int1d==kind) cout << "  -- boundary int border ( nQP: "<< FIE.n << ") ,"  ;
     else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all edges ( nQP: "<< FIE.n << "),"  ;
     else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF edges nQP: ("<< FIE.n << ")," ;
     else cout << "  --  int 2d   (nQP: "<< FIT.n << " ) in "  ;
 }
 // if(di.islevelset()) InternalError("Sorry no levelset integration type on this case (1)");
 if(di.islevelset() && (CDomainOfIntegration::int1d!=kind) &&  (CDomainOfIntegration::int2d!=kind) )
     InternalError("Sorry no levelset integration type on no int1d case");

 /*
  if (verbosity>3)
  if (CDomainOfIntegration::int1d==kind) cout << "  -- boundary int border  " ;
  else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all edges, "   ;
  else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF edges, "   ;
  else cout << "  --  int  in  " ; */
 Expandsetoflab(stack,di, setoflab,all);
 /*
  for (size_t i=0;i<what.size();i++)
  {long  lab  = GetAny<long>( (*what[i])(stack));
  setoflab.insert(lab);
  if ( verbosity>3) cout << lab << " ";
  all=false;
  }*/
 if (verbosity>3) cout <<" Optimized = "<< useopt << ", ";
 const E_F0 * poptiexp0=b->b->optiexp0;
 // const E_F0 & optiexpK=*b->b->optiexpK;
 int n_where_in_stack_opt=b->b->where_in_stack_opt.size();
 R** where_in_stack =0;
 if (n_where_in_stack_opt && useopt)
     where_in_stack = new R * [n_where_in_stack_opt];
 if (where_in_stack)
 {
     assert(b->b->v.size()==(size_t) n_where_in_stack_opt);
     for (int i=0;i<n_where_in_stack_opt;i++)
     {
         int offset=b->b->where_in_stack_opt[i];
         assert(offset>10);
         where_in_stack[i]= static_cast<R *>(static_cast<void *>((char*)stack+offset));
         *(where_in_stack[i])=0;
     }


     if(poptiexp0)
         (*poptiexp0)(stack);
     KN<bool> ok(b->b->v.size());
     {  //   remove the zero coef in the liste
         // R zero=R();
         int il=0;
         for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
             ok[il] =  ! (b->b->mesh_indep_stack_opt[il] && ( std::norm(*(where_in_stack[il])) < 1e-100 ) );
     }
     BilinearOperator b_nozer(*b->b,ok);
     if (verbosity % 10 > 3 )
         cout << "   -- nb term in bilinear form  (!0) : " << b_nozer.v.size()
         << "  total " << n_where_in_stack_opt << endl;

     if ( (verbosity/100) % 10 >= 2)
     {
         int il=0;

         for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
             cout << il << " coef (" << l->first << ") = " << *(where_in_stack[il])
             << " offset=" << b->b->where_in_stack_opt[il]
             << " dep mesh " << l->second.MeshIndependent() << b->b->mesh_indep_stack_opt[il] << endl;
     }
 }
 Stack_Ptr<R*>(stack,ElemMatPtrOffset) =where_in_stack;

 KN<double>  p(Vh.esize()+ Uh.esize() );


 if (verbosity >3)
 {
     if (all) cout << " all " << endl ;
     else cout << endl;
 }

 if (di.kind == CDomainOfIntegration::int1d )
 {

     if(di.islevelset())
     {
         double uset = HUGE_VAL;
         R2 Q[2];
         double vol6[2];
         KN<double> phi(Th.nv);phi=uset;
         double f[3], ll=0;
         for(int t=0; t< Th.nt;++t)
         {
             if ( all || setoflab.find(Th[t].lab) != setoflab.end())
             {
                 double umx=-HUGE_VAL,umn=HUGE_VAL;
                 for(int i=0;i<3;++i)
                 {
                     int j= Th(t,i);
                     if( phi[j]==uset)
                     {
                         MeshPointStack(stack)->setP(&Th,t,i);
                         phi[j]= di.levelset(stack);//zzzz
                     }
                     f[i]=phi[j];
                     umx = std::max(umx,phi[j]);
                     umn = std::min(umn,phi[j]);

                 }
                 int ntp= IsoLineK(f,Q,1e-10);
                 if(verbosity>999 && ntp==2)
                 {
                     const TriangleS &T = Th[t];
                     R3 E(T(Q[0]),T(Q[1]));
                     double le=sqrt((E,E));
                     ll += le;
                     cout << "\t\t" << ntp <<" :  " << Q[0] << " " << Q[1] << " ;  "
                     << f[0] << " " << f[1] << " " << f[2] << "  " << le << " / " << ll<<endl;
                 }
                 if( ntp==2)
                 { //if( withmap)
                    // AddMatElem(mapu,mapt,A,Th,*b->b,sym,t,10,Th[t].lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
                 //else
                     AddMatElem(A,Th,*b->b,sym,t,10,Th[t].lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
                     if(sptrclean) sptrclean=sptr->clean();
                 }
             }
         }
         FIT =FITo;
     }


     else
     {
         for( int e=0;e<Th.nbe;e++)
         {
             if (all || setoflab.find(Th.be(e).lab) != setoflab.end())
             {
                 int ie,i =Th.BoundaryElement(e,ie);
                 //if( withmap)
                 //    AddMatElem(mapu,mapt,A,Th,*b->b,sym,i,ie,Th.be(e).lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
                 //else
                     AddMatElem(A,Th,*b->b,sym,i,ie,Th.be(e).lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
                 if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
             }
         }
     }}
 /*else if (di.kind == CDomainOfIntegration::intalledges)
 {
     cerr << " Sorry no implement to hard  "<< endl;
     ExecError("FH: no intalledges on diff mesh ???");
     ffassert(0); // a faire
     if(withmap)
         for (int i=0;i< Th.nt; i++)
         {
             if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                 for (int ie=0;ie<3;ie++)
                     AddMatElem(mapu,mapt,A,Th,*b->b,sym,i,ie,Th[i].lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
             if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr


         }

     else
         for (int i=0;i< Th.nt; i++)
         {
             if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                 for (int ie=0;ie<3;ie++)
                     AddMatElem(A,Th,*b->b,sym,i,ie,Th[i].lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
             if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr


         }

 }*/
 else if (di.kind == CDomainOfIntegration::intallVFedges)
 {

     cerr << " a faire intallVFedges " << endl;
     ffassert(0);

 }
 else if (di.kind == CDomainOfIntegration::int2d ) {
     // cerr << " a faire CDomainOfIntegration::int2d  " << endl;
     if(di.islevelset())
     {
         double uset = HUGE_VAL;
         R2 Q[2][3];
         double vol6[2];
         KN<double> phi(Th.nv);phi=uset;
         double f[3];
         for(int t=0; t< Th.nt;++t)
         {
             if ( all || setoflab.find(Th[t].lab) != setoflab.end())
             {
                 double umx=-HUGE_VAL,umn=HUGE_VAL;
                 for(int i=0;i<3;++i)
                 {
                     int j= Th(t,i);
                     if( phi[j]==uset)
                     {
                         MeshPointStack(stack)->setP(&Th,t,i);
                         phi[j]= di.levelset(stack);//zzzz
                     }
                     f[i]=phi[j];
                     umx = std::max(umx,phi[j]);
                     umn = std::min(umn,phi[j]);

                 }
                 int nt= UnderIso(f,Q, vol6,1e-14);
                 setQF<R2>(FIT,FITo,QuadratureFormular_T_1, Q,vol6,nt);
                 if(FIT.n)
                 {
                     //if(withmap)
                      //   AddMatElem(mapu,mapt,A,Th,*b->b,sym,t,-1,Th[t].lab,Uh,Vh,FIT,FIE,p,stack);
                     //else
                         AddMatElem(A,Th,*b->b,sym,t,-1,Th[t].lab,Uh,Vh,FIT,FIE,p,stack);
                 }
                 if(sptrclean) sptrclean=sptr->clean();
             }
         }
         FIT =FITo;
     }
     else

     {
         //if(withmap)
          //   for (int i=0;i< Th.nt; i++)
           //  {
            //     if ( all || setoflab.find(Th[i].lab) != setoflab.end())
             //        AddMatElem(mapu,mapt,A,Th,*b->b,sym,i,-1,Th[i].lab,Uh,Vh,FIT,FIE,p,stack);
             //    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
            // }

         //else
             for (int i=0;i< Th.nt; i++)
             {
                 if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                     AddMatElem(A,Th,*b->b,sym,i,-1,Th[i].lab,Uh,Vh,FIT,FIE,p,stack);
                 if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
             }

     }

 }
 else
     InternalError(" kind of CDomainOfIntegration unknown");

 if (where_in_stack) delete [] where_in_stack;
     }





// case 3D volume / 3D Surf on meshS
 template<class R>
 void AssembleBilinearForm(Stack stack,const MeshS & Th,const FESpace3 & Uh,const FESpaceS & Vh,bool sym,
                           MatriceMap<R>  & A, const  FormBilinear * b  )

{
 StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
 bool sptrclean=true;
 //     sptr->clean(); // modif FH mars 2006  clean Ptr

 const CDomainOfIntegration & di= *b->di;
pmeshS  pThdi = GetAny<pmeshS>((*b->di->Th)(stack));

 SHOWVERB(cout << " FormBilinear () " << endl);
 const int useopt=di.UseOpt(stack);
 //double binside=di.binside(stack);
 const bool intmortar=di.intmortar(stack);
 if ( verbosity >1)
 {
     cout << " Integral(8)   on Th "<< &Th << " nv :  " << Th.nv << " nt : " << Th.nt << endl;
     cout << "        Th/ u "<< &Uh.Th << " nv : " << Uh.Th.nv << "   nt : " << Uh.Th.nt << endl;
     cout << "        Th/ v "<< &Vh.Th << " nv : " << Vh.Th.nv << "   nt : " << Vh.Th.nt << endl;
     cout << "        suppose in mortar " << intmortar << "   levelset=  " << di.islevelset() << " withmap: " << di.withmap() << endl;
 }
 Expression  const * const mapt=*di.mapt?di.mapt:0 ;
 Expression  const * const mapu=*di.mapu?di.mapu:0 ;
 bool withmap =di.withmap();
 //   ExecError(" no map  in the case (4) ??");}
 assert(pThdi == & Th);
 //const vector<Expression>  & what(di.what);
 CDomainOfIntegration::typeofkind  kind = di.kind;
 set<int> setoflab;
 bool all=true;
 const QuadratureFormular1d & FIE = di.FIE(stack);
 const QuadratureFormular & FITo = di.FIT(stack);
 QuadratureFormular FIT(FITo,3);
 bool VF=b->VF();  // finite Volume or discontinuous Galerkin
 if (verbosity>2) cout << "  -- discontinuous Galerkin  =" << VF << " size of Mat =" << A.size()<< " Bytes\n";
 if (verbosity>3)
 {
     if (CDomainOfIntegration::int1d==kind) cout << "  -- boundary int border ( nQP: "<< FIE.n << ") ,"  ;
     else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all edges ( nQP: "<< FIE.n << "),"  ;
     else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF edges nQP: ("<< FIE.n << ")," ;
     else cout << "  --  int 2d   (nQP: "<< FIT.n << " ) in "  ;
 }
 // if(di.islevelset()) InternalError("Sorry no levelset integration type on this case (1)");
 if(di.islevelset() && (CDomainOfIntegration::int1d!=kind) &&  (CDomainOfIntegration::int2d!=kind) )
     InternalError("Sorry no levelset integration type on no int1d case");

 /*
  if (verbosity>3)
  if (CDomainOfIntegration::int1d==kind) cout << "  -- boundary int border  " ;
  else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all edges, "   ;
  else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF edges, "   ;
  else cout << "  --  int  in  " ; */
 Expandsetoflab(stack,di, setoflab,all);
 /*
  for (size_t i=0;i<what.size();i++)
  {long  lab  = GetAny<long>( (*what[i])(stack));
  setoflab.insert(lab);
  if ( verbosity>3) cout << lab << " ";
  all=false;
  }*/
 if (verbosity>3) cout <<" Optimized = "<< useopt << ", ";
 const E_F0 * poptiexp0=b->b->optiexp0;
 // const E_F0 & optiexpK=*b->b->optiexpK;
 int n_where_in_stack_opt=b->b->where_in_stack_opt.size();
 R** where_in_stack =0;
 if (n_where_in_stack_opt && useopt)
     where_in_stack = new R * [n_where_in_stack_opt];
 if (where_in_stack)
 {
     assert(b->b->v.size()==(size_t) n_where_in_stack_opt);
     for (int i=0;i<n_where_in_stack_opt;i++)
     {
         int offset=b->b->where_in_stack_opt[i];
         assert(offset>10);
         where_in_stack[i]= static_cast<R *>(static_cast<void *>((char*)stack+offset));
         *(where_in_stack[i])=0;
     }


     if(poptiexp0)
         (*poptiexp0)(stack);
     KN<bool> ok(b->b->v.size());
     {  //   remove the zero coef in the liste
         // R zero=R();
         int il=0;
         for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
             ok[il] =  ! (b->b->mesh_indep_stack_opt[il] && ( std::norm(*(where_in_stack[il])) < 1e-100 ) );
     }
     BilinearOperator b_nozer(*b->b,ok);
     if (verbosity % 10 > 3 )
         cout << "   -- nb term in bilinear form  (!0) : " << b_nozer.v.size()
         << "  total " << n_where_in_stack_opt << endl;

     if ( (verbosity/100) % 10 >= 2)
     {
         int il=0;

         for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
             cout << il << " coef (" << l->first << ") = " << *(where_in_stack[il])
             << " offset=" << b->b->where_in_stack_opt[il]
             << " dep mesh " << l->second.MeshIndependent() << b->b->mesh_indep_stack_opt[il] << endl;
     }
 }
 Stack_Ptr<R*>(stack,ElemMatPtrOffset) =where_in_stack;

 KN<double>  p(Vh.esize()+ Uh.esize() );


 if (verbosity >3)
 {
     if (all) cout << " all " << endl ;
     else cout << endl;
 }

 if (di.kind == CDomainOfIntegration::int1d )
 {

     if(di.islevelset())
     {
         double uset = HUGE_VAL;
         R2 Q[2];
         double vol6[2];
         KN<double> phi(Th.nv);phi=uset;
         double f[3], ll=0;
         for(int t=0; t< Th.nt;++t)
         {
             if ( all || setoflab.find(Th[t].lab) != setoflab.end())
             {
                 double umx=-HUGE_VAL,umn=HUGE_VAL;
                 for(int i=0;i<3;++i)
                 {
                     int j= Th(t,i);
                     if( phi[j]==uset)
                     {
                         MeshPointStack(stack)->setP(&Th,t,i);
                         phi[j]= di.levelset(stack);//zzzz
                     }
                     f[i]=phi[j];
                     umx = std::max(umx,phi[j]);
                     umn = std::min(umn,phi[j]);

                 }
                 int ntp= IsoLineK(f,Q,1e-10);
                 if(verbosity>999 && ntp==2)
                 {
                     const TriangleS &T = Th[t];
                     R3 E(T(Q[0]),T(Q[1]));
                     double le=sqrt((E,E));
                     ll += le;
                     cout << "\t\t" << ntp <<" :  " << Q[0] << " " << Q[1] << " ;  "
                     << f[0] << " " << f[1] << " " << f[2] << "  " << le << " / " << ll<<endl;
                 }
                 if( ntp==2)
                 { //if( withmap)
                    // AddMatElem(mapu,mapt,A,Th,*b->b,sym,t,10,Th[t].lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
                 //else
                     AddMatElem(A,Th,*b->b,sym,t,10,Th[t].lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
                     if(sptrclean) sptrclean=sptr->clean();
                 }
             }
         }
         FIT =FITo;
     }


     else
     {
         for( int e=0;e<Th.nbe;e++)
         {
             if (all || setoflab.find(Th.be(e).lab) != setoflab.end())
             {
                 int ie,i =Th.BoundaryElement(e,ie);
                 //if( withmap)
                 //    AddMatElem(mapu,mapt,A,Th,*b->b,sym,i,ie,Th.be(e).lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
                 //else
                     AddMatElem(A,Th,*b->b,sym,i,ie,Th.be(e).lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
                 if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
             }
         }
     }}
 /*else if (di.kind == CDomainOfIntegration::intalledges)
 {
     cerr << " Sorry no implement to hard  "<< endl;
     ExecError("FH: no intalledges on diff mesh ???");
     ffassert(0); // a faire
     if(withmap)
         for (int i=0;i< Th.nt; i++)
         {
             if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                 for (int ie=0;ie<3;ie++)
                     AddMatElem(mapu,mapt,A,Th,*b->b,sym,i,ie,Th[i].lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
             if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr


         }

     else
         for (int i=0;i< Th.nt; i++)
         {
             if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                 for (int ie=0;ie<3;ie++)
                     AddMatElem(A,Th,*b->b,sym,i,ie,Th[i].lab,Uh,Vh,FIT,FIE,p,stack,intmortar);
             if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr


         }

 }*/
 else if (di.kind == CDomainOfIntegration::intallVFedges)
 {

     cerr << " a faire intallVFedges " << endl;
     ffassert(0);

 }
 else if (di.kind == CDomainOfIntegration::int2d ) {
     // cerr << " a faire CDomainOfIntegration::int2d  " << endl;
     if(di.islevelset())
     {
         double uset = HUGE_VAL;
         R2 Q[2][3];
         double vol6[2];
         KN<double> phi(Th.nv);phi=uset;
         double f[3];
         for(int t=0; t< Th.nt;++t)
         {
             if ( all || setoflab.find(Th[t].lab) != setoflab.end())
             {
                 double umx=-HUGE_VAL,umn=HUGE_VAL;
                 for(int i=0;i<3;++i)
                 {
                     int j= Th(t,i);
                     if( phi[j]==uset)
                     {
                         MeshPointStack(stack)->setP(&Th,t,i);
                         phi[j]= di.levelset(stack);//zzzz
                     }
                     f[i]=phi[j];
                     umx = std::max(umx,phi[j]);
                     umn = std::min(umn,phi[j]);

                 }
                 int nt= UnderIso(f,Q, vol6,1e-14);
                 setQF<R2>(FIT,FITo,QuadratureFormular_T_1, Q,vol6,nt);
                 if(FIT.n)
                 {
                     //if(withmap)
                      //   AddMatElem(mapu,mapt,A,Th,*b->b,sym,t,-1,Th[t].lab,Uh,Vh,FIT,FIE,p,stack);
                     //else
                         AddMatElem(A,Th,*b->b,sym,t,-1,Th[t].lab,Uh,Vh,FIT,FIE,p,stack);
                 }
                 if(sptrclean) sptrclean=sptr->clean();
             }
         }
         FIT =FITo;
     }
     else

     {
         //if(withmap)
          //   for (int i=0;i< Th.nt; i++)
           //  {
            //     if ( all || setoflab.find(Th[i].lab) != setoflab.end())
             //        AddMatElem(mapu,mapt,A,Th,*b->b,sym,i,-1,Th[i].lab,Uh,Vh,FIT,FIE,p,stack);
             //    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
            // }

         //else
             for (int i=0;i< Th.nt; i++)
             {
                 if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                     AddMatElem(A,Th,*b->b,sym,i,-1,Th[i].lab,Uh,Vh,FIT,FIE,p,stack);
                 if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
             }

     }

 }
 else
     InternalError(" kind of CDomainOfIntegration unknown");

 if (where_in_stack) delete [] where_in_stack;
     }


// 3D Surf / 3D curve on meshL
 template<class R>
 void AssembleBilinearForm(Stack stack,const MeshL & Th,const FESpaceS & Uh,const FESpaceL & Vh,bool sym,
                           MatriceMap<R>  & A, const  FormBilinear * b  )

{
    StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
    bool sptrclean=true;

    const CDomainOfIntegration & di= *b->di;
    pmeshL  pThdi = GetAny<pmeshL>((*b->di->Th)(stack));

    SHOWVERB(cout << " FormBilinear () " << endl);

    const int useopt=di.UseOpt(stack);
    //double binside=di.binside(stack);
    const bool intmortar=di.intmortar(stack);
    if ( verbosity >1) {
     cout << " Integral(9)   on Th "<< &Th << " nv :  " << Th.nv << " nt : " << Th.nt << endl;
     cout << "        Th/ u "<< &Uh.Th << " nv : " << Uh.Th.nv << "   nt : " << Uh.Th.nt << endl;
     cout << "        Th/ v "<< &Vh.Th << " nv : " << Vh.Th.nv << "   nt : " << Vh.Th.nt << endl;
     cout << "        suppose in mortar " << intmortar << "   levelset=  " << di.islevelset() << " withmap: " << di.withmap() << endl;
    }
    Expression  const * const mapt=*di.mapt?di.mapt:0 ;
    Expression  const * const mapu=*di.mapu?di.mapu:0 ;
    bool withmap =di.withmap();
    assert(pThdi == & Th);
    CDomainOfIntegration::typeofkind  kind = di.kind;
    set<int> setoflab;
    bool all=true;

    const GQuadratureFormular<R1> & FITo = di.FIE(stack);
    GQuadratureFormular<R1> FIT(FITo,3);

    bool VF=b->VF();  // finite Volume or discontinuous Galerkin
    if (verbosity>2) cout << "  -- discontinuous Galerkin  =" << VF << " size of Mat =" << A.size()<< " Bytes\n";

    // if(di.islevelset()) InternalError("Sorry no levelset integration type on this case (1)");
    if(di.islevelset() && (CDomainOfIntegration::int1d!=kind))
      InternalError("Sorry no levelset integration type on no int1d case");

    Expandsetoflab(stack,di, setoflab,all);

    if (verbosity>3) cout <<" Optimized = "<< useopt << ", ";
    const E_F0 * poptiexp0=b->b->optiexp0;
 // const E_F0 & optiexpK=*b->b->optiexpK;
    int n_where_in_stack_opt=b->b->where_in_stack_opt.size();
    R** where_in_stack =0;
    if (n_where_in_stack_opt && useopt)
      where_in_stack = new R * [n_where_in_stack_opt];
    if (where_in_stack) {
      assert(b->b->v.size()==(size_t) n_where_in_stack_opt);
      for (int i=0;i<n_where_in_stack_opt;i++) {
        int offset=b->b->where_in_stack_opt[i];
        assert(offset>10);
        where_in_stack[i]= static_cast<R *>(static_cast<void *>((char*)stack+offset));
        *(where_in_stack[i])=0;
      }
      if(poptiexp0)
        (*poptiexp0)(stack);
      KN<bool> ok(b->b->v.size());
      int il=0;
      for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
        ok[il] =  ! (b->b->mesh_indep_stack_opt[il] && ( std::norm(*(where_in_stack[il])) < 1e-100 ) );
     
      BilinearOperator b_nozer(*b->b,ok);
      if (verbosity % 10 > 3 ) cout << "   -- nb term in bilinear form  (!0) : " << b_nozer.v.size() << "  total " << n_where_in_stack_opt << endl;

      if ( (verbosity/100) % 10 >= 2) {
        int il=0;
        for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
           cout << il << " coef (" << l->first << ") = " << *(where_in_stack[il]) << " offset=" << b->b->where_in_stack_opt[il]
           << " dep mesh " << l->second.MeshIndependent() << b->b->mesh_indep_stack_opt[il] << endl;
      }
    }
    Stack_Ptr<R*>(stack,ElemMatPtrOffset) =where_in_stack;
    KN<double>  p(Vh.esize()+ Uh.esize() );

    if (verbosity >3) {
      if (all) cout << " all " << endl ;
      else cout << endl;
    }

    if (di.kind == CDomainOfIntegration::int1d ) {
      if(di.islevelset())
        ffassert(0);
      else {
        for (int i=0;i< Th.nt; i++) {
          if ( all || setoflab.find(Th[i].lab) != setoflab.end())
            AddMatElem(A,Th,*b->b,sym,i,-1,Th[i].lab,Uh,Vh,FIT,0,p,stack);
          if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
         }
       }
    }
    else { cout << " di.kind " << di.kind << endl;
      InternalError(" kind of CDomainOfIntegration unknown");
    }
    if (where_in_stack) delete [] where_in_stack;
}





// 3D curve / 3D Surf on meshL
  template<class R>
  void AssembleBilinearForm(Stack stack,const MeshL & Th,const FESpaceL & Uh,const FESpaceS & Vh,bool sym,
                            MatriceMap<R>  & A, const  FormBilinear * b  )

  {
      StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
      bool sptrclean=true;

      const CDomainOfIntegration & di= *b->di;
      pmeshL  pThdi = GetAny<pmeshL>((*b->di->Th)(stack));

      SHOWVERB(cout << " FormBilinear () " << endl);

      const int useopt=di.UseOpt(stack);
      //double binside=di.binside(stack);
      const bool intmortar=di.intmortar(stack);
      if ( verbosity >1)
      {
          cout << " Integral(10)   on Th "<< &Th << " nv :  " << Th.nv << " nt : " << Th.nt << endl;
          cout << "        Th/ u "<< &Uh.Th << " nv : " << Uh.Th.nv << "   nt : " << Uh.Th.nt << endl;
          cout << "        Th/ v "<< &Vh.Th << " nv : " << Vh.Th.nv << "   nt : " << Vh.Th.nt << endl;
          cout << "        suppose in mortar " << intmortar << "   levelset=  " << di.islevelset() << " withmap: " << di.withmap() << endl;
      }
      Expression  const * const mapt=*di.mapt?di.mapt:0 ;
      Expression  const * const mapu=*di.mapu?di.mapu:0 ;
      bool withmap =di.withmap();
      assert(pThdi == & Th);
      CDomainOfIntegration::typeofkind  kind = di.kind;
      set<int> setoflab;
      bool all=true;

      const GQuadratureFormular<R1> & FITo = di.FIE(stack);
      GQuadratureFormular<R1> FIT(FITo,3);

      bool VF=b->VF();  // finite Volume or discontinuous Galerkin
      if (verbosity>2) cout << "  -- discontinuous Galerkin  =" << VF << " size of Mat =" << A.size()<< " Bytes\n";

      // if(di.islevelset()) InternalError("Sorry no levelset integration type on this case (1)");
      if(di.islevelset() && (CDomainOfIntegration::int1d!=kind))
          InternalError("Sorry no levelset integration type on no int1d case");

      Expandsetoflab(stack,di, setoflab,all);

      if (verbosity>3) cout <<" Optimized = "<< useopt << ", ";
      const E_F0 * poptiexp0=b->b->optiexp0;
      // const E_F0 & optiexpK=*b->b->optiexpK;
      int n_where_in_stack_opt=b->b->where_in_stack_opt.size();
      R** where_in_stack =0;
      if (n_where_in_stack_opt && useopt)
          where_in_stack = new R * [n_where_in_stack_opt];
      if (where_in_stack) {
          assert(b->b->v.size()==(size_t) n_where_in_stack_opt);
          for (int i=0;i<n_where_in_stack_opt;i++) {
              int offset=b->b->where_in_stack_opt[i];
              assert(offset>10);
              where_in_stack[i]= static_cast<R *>(static_cast<void *>((char*)stack+offset));
              *(where_in_stack[i])=0;
          }


          if(poptiexp0)
              (*poptiexp0)(stack);
          KN<bool> ok(b->b->v.size());
          int il=0;
          for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
            ok[il] =  ! (b->b->mesh_indep_stack_opt[il] && ( std::norm(*(where_in_stack[il])) < 1e-100 ) );
          
          BilinearOperator b_nozer(*b->b,ok);
          if (verbosity % 10 > 3 )
              cout << "   -- nb term in bilinear form  (!0) : " << b_nozer.v.size()
              << "  total " << n_where_in_stack_opt << endl;

          if ( (verbosity/100) % 10 >= 2) {
              int il=0;

              for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
                  cout << il << " coef (" << l->first << ") = " << *(where_in_stack[il])
                  << " offset=" << b->b->where_in_stack_opt[il]
                  << " dep mesh " << l->second.MeshIndependent() << b->b->mesh_indep_stack_opt[il] << endl;
          }
      }
      Stack_Ptr<R*>(stack,ElemMatPtrOffset) =where_in_stack;

      KN<double>  p(Vh.esize()+ Uh.esize() );


      if (verbosity >3) {
          if (all) cout << " all " << endl ;
          else cout << endl;
      }

      if (di.kind == CDomainOfIntegration::int1d ) {
          if(di.islevelset())
              ffassert(0);
          else {
               for (int i=0;i< Th.nt; i++) {
                  if ( all || setoflab.find(Th[i].lab) != setoflab.end())
                      AddMatElem(A,Th,*b->b,sym,i,-1,Th[i].lab,Uh,Vh,FIT,0,p,stack);
                  if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
              }
            }
      }
      else { cout << " di.kind " << di.kind << endl;
          InternalError(" kind of CDomainOfIntegration unknown");
      }
      if (where_in_stack) delete [] where_in_stack;
  }



    // --------- FH 170605
    ////////////////////////////////////////////////
    // Element_Op for MatriceElementairePleine
    ////////////////////////////////////////////////


    // xxxxxxxxxxxxxxxxx  modif a faire
    // creating an instance of Element_Op with MatriceElementairePleine
    // case 2d
    template<class R>
    void  Element_Op(MatriceElementairePleine<R,FESpace> & mat,const FElement & Ku,const FElement & Kv,double * p,int ie,int label,void *vstack,R2 *B)
    {
        Stack stack=pvoid2Stack(vstack);
        typedef  FElement::Element Element;
        MeshPoint mp= *MeshPointStack(stack);
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);

        bool same = &Ku == & Kv;
        const Element & T  = Ku.T;
        throwassert(&T == &Kv.T);
        const QuadratureFormular & FI = mat.FIT;
        const QuadratureFormular1d & FIb = mat.FIE;
        long npi;
        R *a=mat.a;
        R *pa=a;
        long i,j;
        long n= mat.n,m=mat.m,nx=n*m;
        long N= Kv.N;
        long M= Ku.N;





        const Opera &Op(*mat.bilinearform);
        bool classoptm = copt && Op.optiexpK;
        bool oldopt=1;  // juin 2007 FH ???? a voir
        int  iloop=0;
        KN<bool> unvarexp(classoptm ? Op.optiexpK->sizevar() : 1);
        if (Ku.number<1 && verbosity/100 && verbosity % 10 == 2)
        cout << "Element_Op P: copt = " << copt << " " << classoptm << " opt: " << mat.optim << endl;
        assert(Op.MaxOp() <last_operatortype);


        KN<bool> Dop(last_operatortype);
        Op.DiffOp(Dop);
  //      int lastop=1+Dop.last(binder1st<equal_to<bool> >(equal_to<bool>(),true);
        int lastop=1+Dop.last([](bool x){return x;});
        //assert(lastop<=3);
        RNMK_ fv(p,n,N,lastop); //  the value for basic fonction
        RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction

        for (i=0;i< nx;i++)
        *pa++ = 0.;
        if (ie<0 )//&& B==0)
        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            QuadraturePoint pi(FI[npi]);
            R mes = B ? B->x : T.area;
            R coef = mes *pi.a;
            R2 Pt(pi);
            pa =a;
            Ku.BF(Dop,Pt,fu);
            MeshPointStack(stack)->set(T(Pt),Pt,Kv);
            if (classoptm) {
                if( oldopt) (*Op.optiexpK)(stack); // call old optim version
                else Op.optiexpK->eval(stack,iloop++,unvarexp); // new optim version
            }
            if (!same) Kv.BF(Dop,Pt,fv);
            int il=0;
            for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
            {  // attention la fonction test donne la ligne
                //  et la fonction test est en second
                BilinearOperator::K ll(*l);
                //          pair<int,int> jj(ll.first.first),ii(ll.first.second);
                long jcomp= ll.first.first.first,jop=ll.first.first.second;
                long icomp= ll.first.second.first,iop=ll.first.second.second;

                R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                if ( copt && ( mat.optim==1) && Kv.number <1)
                {
                    R cc  =  GetAny<R>(ll.second.eval(stack));
                    //cout << *(copt[il]) << " == " <<  cc << endl;
                    CheckErrorOptimisation(cc,ccc,"Sorry error in Optimization (e) add:  int2d(Th,optimize=0)(...)");
                   /* if ( ccc != cc) {
                        cerr << cc << " != " << ccc << " => ";
                        cerr << "Sorry error in Optimization (e) add:  int2d(Th,optimize=0)(...)" << endl;
                        ExecError("In Optimized version "); }*/
                }
                int fi=Kv.dfcbegin(icomp);
                int li=Kv.dfcend(icomp);
                int fj=Ku.dfcbegin(jcomp);
                int lj=Ku.dfcend(jcomp);
                ccc *= coef;

                // attention la fonction test donne la ligne
                //  et la fonction test est en second
                for ( i=fi;  i<li;   i++ )
                {
                    for ( j=fj;  j<lj;   j++ )
                    {
                        R w_i =  fv(i,icomp,iop);
                        R w_j =  fu(j,jcomp,jop);
                        mat(i,j) += ccc * w_i*w_j;
                    }
                }
            }
        }
        else if(B)
        {  // int on isovalue ...
            R2 PA(B[0]),PB(B[1]);
            R2 A=T(PA),B=T(PB);
            R2 E(A,B);
            double le = sqrt((E,E));
            //  cout << " xxxx "<< le << " "<< A << " " << B << endl;
            if(le > 1e-15) // bofbof ????
            for (npi=0;npi<FIb.n;npi++) // loop on the integration point
            {
                pa =a;
                QuadratureFormular1dPoint pi( FIb[npi]);
                double coef = le*pi.a;
                double sa=pi.x,sb=1-sa;
                R2 Pt(PA*sa+PB*sb ); //
                Ku.BF(Dop,Pt,fu);
                if (!same) Kv.BF(Dop,Pt,fv);
                // int label=-999999; // a passer en argument
                MeshPointStack(stack)->set(T(Pt),Pt,Kv,-1,R2(E.y,-E.x)/le,-1);
                if (classoptm) (*Op.optiexpK)(stack); // call optim version
                int il=0;
                for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                {  // attention la fonction test donne la ligne
                    //  et la fonction test est en second
                    BilinearOperator::K ll(*l);
                    //          pair<int,int> jj(ll.first.first),ii(ll.first.second);
                    long jcomp= ll.first.first.first,jop=ll.first.first.second;
                    long icomp= ll.first.second.first,iop=ll.first.second.second;


                    R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                    if ( copt && ( mat.optim==1) && Kv.number <1)
                    {
                        R cc  =  GetAny<R>(ll.second.eval(stack));
                        //cout << *(copt[il]) << " == " <<  cc << endl;
                        CheckErrorOptimisation(cc,ccc,"Sorry error in Optimization (f) add:  int2d(Th,optimize=0)(...)");
                       /* if ( ccc != cc) {
                            cerr << cc << " != " << ccc << " => ";
                            cerr << "Sorry error in Optimization (f) add:  int2d(Th,optimize=0)(...)" << endl;
                            ExecError("In Optimized version "); }*/
                    }
                    int fi=Kv.dfcbegin(icomp);
                    int li=Kv.dfcend(icomp);
                    int fj=Ku.dfcbegin(jcomp);
                    int lj=Ku.dfcend(jcomp);
                    ccc *= coef;

                    // attention la fonction test donne la ligne
                    //  et la fonction test est en second

                    for ( i=fi;  i<li;   i++ )
                    {
                        for ( j=fj;  j<lj;   j++ )
                        {
                            R w_i =  fv(i,icomp,iop);
                            R w_j =  fu(j,jcomp,jop);
                            mat(i,j) += ccc * w_i*w_j;
                        }
                    }
                }
            }
        }
        else // int on edge ie
        for (npi=0;npi<FIb.n;npi++) // loop on the integration point
        {
            pa =a;
            QuadratureFormular1dPoint pi( FIb[npi]);
            R2 E=T.Edge(ie);
            double le = sqrt((E,E));
            double coef = le*pi.a;
            double sa=pi.x,sb=1-sa;
            R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
            PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]);
            R2 Pt(PA*sa+PB*sb ); //
            Ku.BF(Dop,Pt,fu);
            if (!same) Kv.BF(Dop,Pt,fv);
            // int label=-999999; // a passer en argument
            MeshPointStack(stack)->set(T(Pt),Pt,Kv,label,R2(E.y,-E.x)/le,ie);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            int il=0;
            for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
            {  // attention la fonction test donne la ligne
                //  et la fonction test est en second
                BilinearOperator::K ll(*l);
                //          pair<int,int> jj(ll.first.first),ii(ll.first.second);
                long jcomp= ll.first.first.first,jop=ll.first.first.second;
                long icomp= ll.first.second.first,iop=ll.first.second.second;


                R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                if ( copt && ( mat.optim==1) && Kv.number <1)
                {
                    R cc  =  GetAny<R>(ll.second.eval(stack));
                    //cout << *(copt[il]) << " == " <<  cc << endl;
                    CheckErrorOptimisation(cc,ccc,"Sorry error in Optimization (g) add:  int2d(Th,optimize=0)(...)");
                   /* if ( ccc != cc) {
                        cerr << cc << " != " << ccc << " => ";
                        cerr << "Sorry error in Optimization (g) add:  int2d(Th,optimize=0)(...)" << endl;
                        ExecError("In Optimized version "); }*/
                }
                int fi=Kv.dfcbegin(icomp);
                int li=Kv.dfcend(icomp);
                int fj=Ku.dfcbegin(jcomp);
                int lj=Ku.dfcend(jcomp);
                ccc *= coef;

                // attention la fonction test donne la ligne
                //  et la fonction test est en second

                for ( i=fi;  i<li;   i++ )
                {
                    for ( j=fj;  j<lj;   j++ )
                    {
                        R w_i =  fv(i,icomp,iop);
                        R w_j =  fu(j,jcomp,jop);
                        mat(i,j) += ccc * w_i*w_j;
                    }
                }
            }

            /*
             for ( i=0;  i<n;   i++ )
             // if (onWhatIsEdge[ie][Kv.DFOnWhat(i)]) // juste the df on edge bofbof generaly wrong FH dec 2003
             {
             RNM_ wi(fv(i,'.','.'));
             for ( j=0;  j<m;   j++,pa++ )
             {
             RNM_ wj(fu(j,'.','.'));
             int il=0;
             for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
             // if (onWhatIsEdge[ie][Kv.DFOnWhat(j)]) // juste the df on edge bofbof generaly wrong FH dec 2003
             {
             BilinearOperator::K ll(*l);
             pair<int,int> jj(ll.first.first),ii(ll.first.second);

             double w_i =  wi(ii.first,ii.second);
             double w_j =  wj(jj.first,jj.second);
             // R ccc = GetAny<R>(ll.second.eval(stack));

             R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
             if ( copt && ( mat.optim==1) && Kv.number <1)
             {
             R cc  =  GetAny<R>(ll.second.eval(stack));
             if ( ccc != cc) {
             cerr << cc << " != " << ccc << " => ";
             cerr << "Sorry error in Optimization (h) add:  int2d(Th,optimize=0)(...)" << endl;
             ExecError("In Optimized version "); }
             }
             *pa += coef * ccc * w_i*w_j;
             }
             }
             }
             // else pa += m;  FH dec 2003
             */
        }


        /*  pa=a;
         if (Ku.Vh.Th(T) >=0 ) {
         cout <<endl  << " Triangle " << Ku.Vh.Th(T) << " =  "<<  T[0] << ", " << T[1] << ", " << T[2] << " " << nx << endl;
         for (int i=0;i<n;i++)
         {
         cout << setw(2) << i << setw(4) << mat.ni[i] << " :";
         for (int j=0;j<m;j++)
         cout << setw(5)  << (*pa++) << " ";
         cout << endl;
         } }
         */
        *MeshPointStack(stack) = mp;
    }



    // creating an instance of Element_Op with MatriceElementairePleine
    // case 3D volume
    template<class R>
    void  Element_Op(MatriceElementairePleine<R,FESpace3> & mat,const FElement3 & Ku,const FElement3 & Kv,double * p,int ie,int label,void *vstack,R3 *B)
    {
        //  ffassert(B==0);
        Stack stack=pvoid2Stack(vstack);
        //    ffassert(0);
        typedef  FElement3::Element Element;
        MeshPoint mp= *MeshPointStack(stack);
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);

        bool same = &Ku == & Kv;
        const Element & T  = Ku.T;
        throwassert(&T == &Kv.T);
        const GQuadratureFormular<R3> & FI = mat.FIT;
        const GQuadratureFormular<R2> & FIb = mat.FIE;
        long npi;
        R *a=mat.a;
        R *pa=a;
        long i,j;
        long n= mat.n,m=mat.m,nx=n*m;
        long N= Kv.N;
        long M= Ku.N;





        const Opera &Op(*mat.bilinearform);
        bool classoptm = copt && Op.optiexpK;
        bool oldopt=1;  // juin 2007 FH ???? a voir
        int  iloop=0;
        KN<bool> unvarexp(classoptm ? Op.optiexpK->sizevar() : 1);
        if (Ku.number<1 && verbosity/100 && verbosity % 10 == 2)
        cout << "Element_Op 3d P: copt = " << copt << " " << classoptm << " opt: " << mat.optim << endl;
        assert(Op.MaxOp() <last_operatortype);
        //
        int lastop;
        lastop = 0;
        What_d Dop = Op.DiffOp(lastop);
        //KN<bool> Dop(last_operatortype);
        //p.DiffOp(Dop);
        //int lastop=1+Dop.last(binder1st<equal_to<bool> >(equal_to<bool>(),true));
        //assert(lastop<=3);
        RNMK_ fv(p,n,N,lastop); //  the value for basic fonction
        RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction

        for (i=0;i< nx;i++)
        *pa++ = 0.;
        if (ie<0)
        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            GQuadraturePoint<R3> pi(FI[npi]);
            R coef = T.mesure()*pi.a;
            R3 Pt(pi);
            pa =a;
            Ku.BF(Dop,Pt,fu);
            MeshPointStack(stack)->set(T(Pt),Pt,Kv);
            if (classoptm) {
                if( oldopt) (*Op.optiexpK)(stack); // call old optim version
                else Op.optiexpK->eval(stack,iloop++,unvarexp); // new optim version
            }
            if (!same) Kv.BF(Dop,Pt,fv);
            int il=0;
            for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
            {  // attention la fonction test donne la ligne
                //  et la fonction test est en second
                BilinearOperator::K ll(*l);
                //          pair<int,int> jj(ll.first.first),ii(ll.first.second);
                long jcomp= ll.first.first.first,jop=ll.first.first.second;
                long icomp= ll.first.second.first,iop=ll.first.second.second;

                R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                if ( copt && ( mat.optim==1) && Kv.number <1)
                {
                    R cc  =  GetAny<R>(ll.second.eval(stack));
                    //cout << *(copt[il]) << " == " <<  cc << endl;
                    CheckErrorOptimisation(cc,ccc,"Sorry error in Optimization Element_Op plein 3d (a) add:  int2d(Th,optimize=0)(...)");
                     /*   if ( ccc != cc) {
                            cerr << cc << " != " << ccc << " => ";
                            cerr << "Sorry error in Optimization Element_Op plein 3d (a) add:  int2d(Th,optimize=0)(...)" << endl;
                            ExecError("In Optimized version "); }*/
                    }
                int fi=Kv.dfcbegin(icomp);
                int li=Kv.dfcend(icomp);
                int fj=Ku.dfcbegin(jcomp);
                int lj=Ku.dfcend(jcomp);
                ccc *= coef;

                // attention la fonction test donne la ligne
                //  et la fonction test est en second

                for ( i=fi;  i<li;   i++ )
                {
                    for ( j=fj;  j<lj;   j++ )
                    {
                        R w_i =  fv(i,icomp,iop);
                        R w_j =  fu(j,jcomp,jop);
                        mat(i,j) += ccc * w_i*w_j;
                    }
                }
            }
        }
        else if(B)
        {  // int on leveset
            int np = ie-10; //= (B[0].x == B[3].x ) && (B[0].y == B[3].y ) && (B[0].z == B[3].z ) ? 3 : 4;
            if(verbosity>999) cout << "    Ass mat pleine /"<< np << endl;
            assert( np==3 || np==4);
            // XXXXXXX
            double epsmes3=T.mesure()*T.mesure()*1e-18;
            R3 PP[4];
            double l[3];
            for(int i=0; i< np; ++i)
            PP[i]= T(B[i]);

            for( int i =0; i+1 < np; i+=2)
            { // 0,1,, a and 2,3,0.
                int i0=i,i1=i+1,i2=(i+2)%np;
                R3 NN= R3(PP[i0],PP[i1])^R3(PP[i0],PP[i2]);
                double mes2 = (NN,NN);
                double mes = sqrt(mes2);

                if(mes2*mes <epsmes3) continue; //  too small
                NN /= mes;
                mes *= 0.5;
                if(verbosity>999)
                cout << " --int on leveset3d " << np << " " << mes << " " << i0<<i1<<i2 <<endl;
                double asum=0;
                for (npi=0;npi<FIb.n;npi++) // loop on the integration point
                {
                    GQuadraturePoint<R2>  pi( FIb[npi]);
                    // cout << " %% " << npi << " " << pi.a << " " << pi.x << " " << pi.y << endl;
                    asum+= pi.a;
                    pi.toBary(l);
                    R3 Pt( l[0]*B[i0]+l[1]*B[i1]+l[2]*B[i2]); //
                    double coef = mes*pi.a; // correction 0.5 050109 FH
                    Ku.BF(Dop,Pt,fu);
                    if (!same) Kv.BF(Dop,Pt,fv);
                    MeshPointStack(stack)->set(T(Pt),Pt,Ku,label,NN,ie);
                    if (classoptm) (*Op.optiexpK)(stack); // call optim version

                    pa=a;
                    for (int i=0;  i<n;   i++ )
                    {
                        RNM_ wi(fv(i,'.','.'));
                        for (int  j=0;  j<m;   j++,pa++ )
                        {
                            RNM_ wj(fu(j,'.','.'));
                            int il=0;
                            for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                            {
                                BilinearOperator::K ll(*l);
                                pair<int,int> jj(ll.first.first),ii(ll.first.second);

                                double w_i =  wi(ii.first,ii.second);
                                double w_j =  wj(jj.first,jj.second);

                                R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                                if ( copt && ( mat.optim==1) && Kv.number <1)
                                {
                                    R cc  =  GetAny<R>(ll.second.eval(stack));
                                    CheckErrorOptimisation(cc,ccc,"Sorry error in Optimization  Element_Op plein 3d (b) add:  int2d(Th,optimize=0)(...)");
                                   /* if ( ccc != cc) {
                                        cerr << cc << " != " << ccc << " => ";
                                        cerr << "Sorry error in Optimization  Element_Op plein 3d (b) add:  int2d(Th,optimize=0)(...)" << endl;
                                        ExecError("In Optimized version "); }*/
                                }
                                if(verbosity>999)
                                cout << " -- int on leveset3d  aij = "<< pi.a* ccc * w_i*w_j <<" " << ccc << " " << w_i*w_j <<endl;
                                *pa += coef * ccc * w_i*w_j;
                            }
                        }
                    }
                    if(verbosity>999) cout << " ++\n";
                }

                if(verbosity>999) cout << " @@ "<< asum << endl;;

            }

        }// end int level set ...
        else // int on edge ie
        for (npi=0;npi<FIb.n;npi++) // loop on the integration point
        {
            pa =a;
            GQuadraturePoint<R2> pi( FIb[npi]);
            R3 NN= T.N(ie);
            double mes=NN.norme();
            NN/=mes;
            double coef = 0.5*mes*pi.a; // correction 0.5 050109 FH
            R3 Pt(T.PBord(ie,pi));
            Ku.BF(Dop,Pt,fu);
            if (!same) Kv.BF(Dop,Pt,fv);
            MeshPointStack(stack)->set(T(Pt),Pt,Ku,label,NN,ie);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version


            for ( i=0;  i<n;   i++ )
            {
                RNM_ wi(fv(i,'.','.'));
                for ( j=0;  j<m;   j++,pa++ )
                {
                    RNM_ wj(fu(j,'.','.'));
                    int il=0;
                    for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                    {
                        BilinearOperator::K ll(*l);
                        pair<int,int> jj(ll.first.first),ii(ll.first.second);

                        double w_i =  wi(ii.first,ii.second);
                        double w_j =  wj(jj.first,jj.second);

                        R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                        if ( copt && ( mat.optim==1) && Kv.number <1)
                        {
                            R cc  =  GetAny<R>(ll.second.eval(stack));
                            CheckErrorOptimisation(cc,ccc,"Sorry error in Optimization Element_Op plein 3d (c) add:  int2d(Th,optimize=0)(...)");
                            /*
                            if ( ccc != cc) {
                                cerr << cc << " != " << ccc << " => ";
                                cerr << "Sorry error in Optimization  Element_Op plein 3d (c) add:  int2d(Th,optimize=0)(...)" << endl;
                                ExecError("In Optimized version "); }*/
                        }
                        *pa += coef * ccc * w_i*w_j;
                    }
                }
            }
        }


        if (Ku.Vh.Th(T) <1 && verbosity>100) {
            pa=mat.a;
            cout <<endl  << " Tet " << Ku.Vh.Th(T) << " =  " << T  << " " << nx << endl;
            for (int i=0;i<n;i++)
            {
                cout << setw(2) << i << setw(4) << mat.ni[i] << " :";
                for (int j=0;j<m;j++)
                cout << setw(5)  << (*pa++) << " ";
                cout << endl;
            } }


    }


    // creating an instance of Element_Op with MatriceElementairePleine
    // case 3D surface
    template<class R>
    void  Element_Op(MatriceElementairePleine<R,FESpaceS> & mat,const FElementS & Ku,const FElementS & Kv,double * p,int ie,int label,void *vstack,R3 *B)
    {
        Stack stack=pvoid2Stack(vstack);
        typedef  FElementS::Element Element;

        MeshPoint mp= *MeshPointStack(stack);
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);

        bool same = &Ku == & Kv;
        const Element & T  = Ku.T;
        throwassert(&T == &Kv.T);
        const QuadratureFormular & FI = mat.FIT;
        const QuadratureFormular1d & FIb = mat.FIE;
        long npi;
        R *a=mat.a;
        R *pa=a;
        long i,j;
        long n= mat.n,m=mat.m,nx=n*m;
        long N= Kv.N;
        long M= Ku.N;


        const Opera &Op(*mat.bilinearform);
        bool classoptm = copt && Op.optiexpK;
        bool oldopt=1;  // juin 2007 FH ???? a voir
        int  iloop=0;
        KN<bool> unvarexp(classoptm ? Op.optiexpK->sizevar() : 1);
        if (Ku.number<1 && verbosity/100 && verbosity % 10 == 2)
        cout << "Element_Op 3d P: copt = " << copt << " " << classoptm << " opt: " << mat.optim << endl;
        assert(Op.MaxOp() <last_operatortype);

        int lastop;
        lastop = 0;
        What_d Dop = Op.DiffOp(lastop);

        //assert(lastop<=3);
        RNMK_ fv(p,n,N,lastop); //  the value for basic fonction
        RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction

        for (i=0;i< nx;i++)
        *pa++ = 0.;
        if (ie<0 )//&& B==0)
        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            QuadraturePoint pi(FI[npi]);
            R mes = B ? B->x : T.mesure();
            R coef = mes *pi.a;
            R2 Pt(pi);
            pa =a;
            Ku.BF(Dop,Pt,fu);
            MeshPointStack(stack)->set(T(Pt),Pt,Kv);
            if (classoptm) {
                if( oldopt) (*Op.optiexpK)(stack); // call old optim version
                else Op.optiexpK->eval(stack,iloop++,unvarexp); // new optim version
            }
            if (!same) Kv.BF(Dop,Pt,fv);
            int il=0;
            for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
            {  // attention la fonction test donne la ligne
                //  et la fonction test est en second
                BilinearOperator::K ll(*l);
                //          pair<int,int> jj(ll.first.first),ii(ll.first.second);
                long jcomp= ll.first.first.first,jop=ll.first.first.second;
                long icomp= ll.first.second.first,iop=ll.first.second.second;

                R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                if ( copt && ( mat.optim==1) && Kv.number <1)
                {
                    R cc  =  GetAny<R>(ll.second.eval(stack));
                    //cout << *(copt[il]) << " == " <<  cc << endl;
                    CheckErrorOptimisation(cc,ccc,"Sorry error in Optimization (e) add:  int2d(Th,optimize=0)(...)");
                    /*if ( ccc != cc) {
                        cerr << cc << " != " << ccc << " => ";
                        cerr << "Sorry error in Optimization (e) add:  int2d(Th,optimize=0)(...)" << endl;
                        ExecError("In Optimized version "); }*/
                }
                int fi=Kv.dfcbegin(icomp);
                int li=Kv.dfcend(icomp);
                int fj=Ku.dfcbegin(jcomp);
                int lj=Ku.dfcend(jcomp);
                ccc *= coef;

                // attention la fonction test donne la ligne
                //  et la fonction test est en second
                for ( i=fi;  i<li;   i++ )
                {
                    for ( j=fj;  j<lj;   j++ )
                    {
                        R w_i =  fv(i,icomp,iop);
                        R w_j =  fu(j,jcomp,jop);
                        mat(i,j) += ccc * w_i*w_j;
                    }
                }
            }
        }
        else if(B)
         ffassert(0);
        else // int on edge ie
        for (npi=0;npi<FIb.n;npi++) // loop on the integration point
        {
            pa =a;
            QuadratureFormular1dPoint pi( FIb[npi]);
            R3 E=T.Edge(ie);
            double le = sqrt((E,E));
            double coef = le*pi.a;
            double sa=pi.x,sb=1-sa;
            R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
            PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]);
            R2 Pt(PA*sa+PB*sb ); //
            Ku.BF(Dop,Pt,fu);
            // surface normal
            R3 NNt=T.NormalTUnitaire();
            // exterior normal (flux)
            R3 NN=T.N(ie);
            NN /= NN.norme();
            if (!same) Kv.BF(Dop,Pt,fv);
            // int label=-999999; // a passer en argument
            MeshPointStack(stack)->set(T(Pt),Pt,Kv,label,NN,NNt,ie);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            int il=0;
            for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
            {  // attention la fonction test donne la ligne
                //  et la fonction test est en second
                BilinearOperator::K ll(*l);
                //          pair<int,int> jj(ll.first.first),ii(ll.first.second);
                long jcomp= ll.first.first.first,jop=ll.first.first.second;
                long icomp= ll.first.second.first,iop=ll.first.second.second;


                R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                if ( copt && ( mat.optim==1) && Kv.number <1)
                {
                    R cc  =  GetAny<R>(ll.second.eval(stack));
                    //cout << *(copt[il]) << " == " <<  cc << endl;
                    CheckErrorOptimisation(cc,ccc,"Sorry error in Optimization (g) add:  int2d(Th,optimize=0)(...)");
                    /*
                    if ( ccc != cc) {
                        cerr << cc << " != " << ccc << " => ";
                        cerr << "Sorry error in Optimization (g) add:  int2d(Th,optimize=0)(...)" << endl;
                        ExecError("In Optimized version "); }*/
                }
                int fi=Kv.dfcbegin(icomp);
                int li=Kv.dfcend(icomp);
                int fj=Ku.dfcbegin(jcomp);
                int lj=Ku.dfcend(jcomp);
                ccc *= coef;

                // attention la fonction test donne la ligne
                //  et la fonction test est en second

                for ( i=fi;  i<li;   i++ )
                {
                    for ( j=fj;  j<lj;   j++ )
                    {
                        R w_i =  fv(i,icomp,iop);
                        R w_j =  fu(j,jcomp,jop);
                        mat(i,j) += ccc * w_i*w_j;
                    }
                }
            }

            /*
             for ( i=0;  i<n;   i++ )
             // if (onWhatIsEdge[ie][Kv.DFOnWhat(i)]) // juste the df on edge bofbof generaly wrong FH dec 2003
             {
             RNM_ wi(fv(i,'.','.'));
             for ( j=0;  j<m;   j++,pa++ )
             {
             RNM_ wj(fu(j,'.','.'));
             int il=0;
             for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
             // if (onWhatIsEdge[ie][Kv.DFOnWhat(j)]) // juste the df on edge bofbof generaly wrong FH dec 2003
             {
             BilinearOperator::K ll(*l);
             pair<int,int> jj(ll.first.first),ii(ll.first.second);

             double w_i =  wi(ii.first,ii.second);
             double w_j =  wj(jj.first,jj.second);
             // R ccc = GetAny<R>(ll.second.eval(stack));

             R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
             if ( copt && ( mat.optim==1) && Kv.number <1)
             {
             R cc  =  GetAny<R>(ll.second.eval(stack));
             if ( ccc != cc) {
             cerr << cc << " != " << ccc << " => ";
             cerr << "Sorry error in Optimization (h) add:  int2d(Th,optimize=0)(...)" << endl;
             ExecError("In Optimized version "); }
             }
             *pa += coef * ccc * w_i*w_j;
             }
             }
             }
             // else pa += m;  FH dec 2003
             */
        }


        /*  pa=a;
         if (Ku.Vh.Th(T) >=0 ) {
         cout <<endl  << " Triangle " << Ku.Vh.Th(T) << " =  "<<  T[0] << ", " << T[1] << ", " << T[2] << " " << nx << endl;
         for (int i=0;i<n;i++)
         {
         cout << setw(2) << i << setw(4) << mat.ni[i] << " :";
         for (int j=0;j<m;j++)
         cout << setw(5)  << (*pa++) << " ";
         cout << endl;
         } }
         */
        *MeshPointStack(stack) = mp;
    }

    // creating an instance of Element_Op with MatriceElementairePleine
    // case 3D curve
    template<class R>
    void  Element_Op(MatriceElementairePleine<R,FESpaceL> & mat,const FElementL & Ku,const FElementL & Kv,double * p,int ie,int label,void *vstack,R3 *B)
    {
        Stack stack=pvoid2Stack(vstack);
        typedef  FElementL::Element Element;

        MeshPoint mp= *MeshPointStack(stack);
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);

        bool same = &Ku == & Kv;
        const Element & T  = Ku.T;
        throwassert(&T == &Kv.T);
        //const QuadratureFormular & FI = mat.FIT;
        const GQuadratureFormular<R1> & FI = mat.FIT;
        long npi;
        R *a=mat.a;
        R *pa=a;
        long i,j;
        long n= mat.n,m=mat.m,nx=n*m;
        long N= Kv.N;
        long M= Ku.N;


        const Opera &Op(*mat.bilinearform);
        bool classoptm = copt && Op.optiexpK;
        bool oldopt=1;  // juin 2007 FH ???? a voir
        int  iloop=0;
        KN<bool> unvarexp(classoptm ? Op.optiexpK->sizevar() : 1);
        if (Ku.number<1 && verbosity/100 && verbosity % 10 == 2)
            cout << "Element_Op 3d P: copt = " << copt << " " << classoptm << " opt: " << mat.optim << endl;
        assert(Op.MaxOp() <last_operatortype);

        int lastop;
        lastop = 0;
        What_d Dop = Op.DiffOp(lastop);

        //assert(lastop<=3);
        RNMK_ fv(p,n,N,lastop); //  the value for basic fonction
        RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction

        for (i=0;i< nx;i++)
            *pa++ = 0.;
        int ll=-1; //bof bof
         R3 NNt=T.TangenteUnitaire();
        R3 NN;
        if (ie<0 )//&& B==0)
            for (npi=0;npi<FI.n;npi++) // loop on the integration point
            {
                GQuadraturePoint<R1> pi(FI[npi]);
                R mes = B ? B->x : T.mesure();
                 R coef = mes *pi.a;
                R1 Pt(pi);
                pa =a;
                Ku.BF(Dop,Pt,fu);
                MeshPointStack(stack)->set(T(Pt),Pt,Kv,-1,NN,NNt,-1);// non on boundary ,NNt,ll);//,label,NN,NNt,ie);   Axel
                if (classoptm) {
                    if( oldopt) (*Op.optiexpK)(stack); // call old optim version
                    else Op.optiexpK->eval(stack,iloop++,unvarexp); // new optim version
                }
                if (!same) Kv.BF(Dop,Pt,fv);
                int il=0;
                for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                {  // attention la fonction test donne la ligne
                    //  et la fonction test est en second
                    BilinearOperator::K ll(*l);
                    //          pair<int,int> jj(ll.first.first),ii(ll.first.second);
                    long jcomp= ll.first.first.first,jop=ll.first.first.second;
                    long icomp= ll.first.second.first,iop=ll.first.second.second;

                    R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                    if ( copt && ( mat.optim==1) && Kv.number <1)
                    {
                        R cc  =  GetAny<R>(ll.second.eval(stack));
                        //cout << *(copt[il]) << " == " <<  cc << endl;
                        CheckErrorOptimisation(cc,ccc,"Sorry error in Optimization (e) add:  int2d(Th,optimize=0)(...)");
                        /*if ( ccc != cc) {
                         cerr << cc << " != " << ccc << " => ";
                         cerr << "Sorry error in Optimization (e) add:  int2d(Th,optimize=0)(...)" << endl;
                         ExecError("In Optimized version "); }*/
                    }
                    int fi=Kv.dfcbegin(icomp);
                    int li=Kv.dfcend(icomp);
                    int fj=Ku.dfcbegin(jcomp);
                    int lj=Ku.dfcend(jcomp);
                    ccc *= coef;

                    // attention la fonction test donne la ligne
                    //  et la fonction test est en second
                    for ( i=fi;  i<li;   i++ )
                    {
                        for ( j=fj;  j<lj;   j++ )
                        {
                            R w_i =  fv(i,icomp,iop);
                            R w_j =  fu(j,jcomp,jop);
                            mat(i,j) += ccc * w_i*w_j;
                        }
                    }
                }
            }
        else if(B)
            ffassert(0);
        else // int on vetex  ie
  { // Add F.H sep 2020 for int0d(Th)( u*v )
      double s = ie;
      GQuadraturePoint<R1> pi(1.,s);
      R mes = 1.;
      R3 NN=NNt;
      if(ie==0) NN=-NN;
      R coef = 1.;
       R1 Pt(pi);
       pa =a;
       Ku.BF(Dop,Pt,fu);
       MeshPointStack(stack)->set(T(Pt),Pt,Kv,label,NN,NNt,ie);//   Axel
       if (classoptm) {
           if( oldopt) (*Op.optiexpK)(stack); // call old optim version
           else Op.optiexpK->eval(stack,iloop++,unvarexp); // new optim version
       }
       if (!same) Kv.BF(Dop,Pt,fv);
       int il=0;
       for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
       {  // attention la fonction test donne la ligne
           //  et la fonction test est en second
           BilinearOperator::K ll(*l);
           //          pair<int,int> jj(ll.first.first),ii(ll.first.second);
           long jcomp= ll.first.first.first,jop=ll.first.first.second;
           long icomp= ll.first.second.first,iop=ll.first.second.second;

           R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
           if ( copt && ( mat.optim==1) && Kv.number <1)
           {
               R cc  =  GetAny<R>(ll.second.eval(stack));
               //cout << *(copt[il]) << " == " <<  cc << endl;
               CheckErrorOptimisation(cc,ccc,"Sorry error in Optimization (e) add:  int2d(Th,optimize=0)(...)");
               /*if ( ccc != cc) {
                cerr << cc << " != " << ccc << " => ";
                cerr << "Sorry error in Optimization (e) add:  int2d(Th,optimize=0)(...)" << endl;
                ExecError("In Optimized version "); }*/
           }
           int fi=Kv.dfcbegin(icomp);
           int li=Kv.dfcend(icomp);
           int fj=Ku.dfcbegin(jcomp);
           int lj=Ku.dfcend(jcomp);
           ccc *= coef;

           // attention la fonction test donne la ligne
           //  et la fonction test est en second
           for ( i=fi;  i<li;   i++ )
           {
               for ( j=fj;  j<lj;   j++ )
               {
                   R w_i =  fv(i,icomp,iop);
                   R w_j =  fu(j,jcomp,jop);
                   mat(i,j) += ccc * w_i*w_j;
               }
           }
       }
   }
            

        *MeshPointStack(stack) = mp;
    }



    ////////////////////////////////////////////////
    // Element_Op for MatriceElementaireSymetrique
    ////////////////////////////////////////////////
    // using to define new solver

    // creating an instance of Element_Op with MatriceElementaireSymetrique
    // case 2d
    // xxxxxxxxxxxxxxxxx  modif a faire
    template<class R>
    void  Element_Op(MatriceElementaireSymetrique<R,FESpace> & mat,const FElement & Ku,double * p,int ie,int label, void * vstack,R2*B)
    {
        Stack stack=pvoid2Stack(vstack);
        MeshPoint mp= *MeshPointStack(stack);
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        const Triangle & T  = Ku.T;
        //  const QuadratureFormular & FI = QuadratureFormular_T_2;
        //  const QuadratureFormular1d & FIb = QF_GaussLegendre2;
        const QuadratureFormular & FI = mat.FIT;
        const QuadratureFormular1d & FIb = mat.FIE;
        long npi;
        R *a=mat.a;
        R *pa=a;
        long i,j;
        long n= mat.n,m=mat.m,nx=n*(m+1)/2;
        long N= Ku.N;
        //long M=N;
        // bool show = Ku.Vh.Th(T)==0;
        //    char * xxx[] ={" u"," v"," p"," q"," r"};
        //char * xxxx[] ={" u'"," v'"," p'"," q'"," r'"};
        //char * yyy[] ={" ","_x ","_y "};


        throwassert(mat.bilinearform);

        const Opera &Op(*mat.bilinearform);
        bool classoptm = copt && Op.optiexpK;
        // assert(  (copt !=0) ||  (Op.where_in_stack_opt.size() !=0) );
        if (Ku.number<1  && verbosity/100 && verbosity % 10 == 2 )
        cout << "Element_Op S: copt = " << copt << " " << classoptm << " opt "<< mat.optim << endl;
        assert(Op.MaxOp() <last_operatortype);


        KN<bool> Dop(last_operatortype);
        Op.DiffOp(Dop);
        int lastop=1+Dop.last([](bool x){return x;});
        // assert(lastop<=3);

        RNMK_ fu(p,n,N,lastop); //  the value for basic fonction

        pa =a;
        for (i=0;i< nx;i++)
        *pa++ = 0.;

        if (ie<0)
        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            QuadraturePoint pi(FI[npi]);
            double mes= B ? B->x :T.area;
            double coef = mes*pi.a;
            R2 Pt(pi);
            pa =a;
            Ku.BF(Dop,Pt,fu);
            MeshPointStack(stack)->set(T(pi),pi,Ku);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            int il=0;
            for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
            {  // attention la fonction test donne la ligne
                //  et la fonction test est en second
                BilinearOperator::K ll(*l);
                //          pair<int,int> jj(ll.first.first),ii(ll.first.second);
                long jcomp= ll.first.first.first,jop=ll.first.first.second;
                long icomp= ll.first.second.first,iop=ll.first.second.second;

                R c = copt ? *(copt[il]): GetAny<R>(ll.second.eval(stack));
                if ( copt && Ku.number <1)
                {
                    R cc  =  GetAny<R>(ll.second.eval(stack));
                    // cout << *(copt[il]) << " == " <<  cc << endl;
                    CheckErrorOptimisation(c,cc,"Sorry error in Optimization (l) add:  int2d(Th,optimize=0)(...)");
                    /*
                    if ( c != cc) {
                        cerr << c << " != " << cc << " => ";
                        cerr << "Sorry error in Optimization (l) add:  int2d(Th,optimize=0)(...)" << endl;
                        ExecError("In Optimized version "); }*/
                }
                c *= coef ;
                long fi=Ku.dfcbegin(icomp);
                long li=Ku.dfcend(icomp);
                long fj=Ku.dfcbegin(jcomp);
                long lj=Ku.dfcend(jcomp);
                if (verbosity>10 && Ku.Vh.Th(T) < 1 && npi < 1)
                cout << " ic "<< icomp << fi<< " "<< lj << " "<< " c "<< jcomp << " " <<fj << " "<< lj << endl;
                for ( i=fi;  i<li;   i++ )
                for ( j=fj;  j<min(lj,i+1);  j++ ) //
                {
                    R w_i =  fu(i,icomp,iop);
                    R w_j =  fu(j,jcomp,jop);

                    mat(i,j)  +=  c * w_i*w_j;

                }

            }

            /*
             for ( i=0;  i<n;   i++ )
             {
             RNM_ wi(fu(i,'.','.'));
             //    if (Ku.Vh.Th(T) < 1) cout << i <<" " <<Pt<< "wi =" << wi ;
             for ( j=0;  j<=i;  j++,pa++ ) //
             {

             RNM_ wj(fu(j,'.','.'));
             //   if (Ku.Vh.Th(T) < 1) cout << j <<" " <<Pt<< "wj =" << wj ;
             int il=0;
             for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
             {
             const  BilinearOperator::K & ll(*l);
             pair<int,int> ii(ll.first.first),jj(ll.first.second);
             double w_i =  wi(ii.first,ii.second);
             double w_j =  wj(jj.first,jj.second);

             R c = copt ? *(copt[il]): GetAny<R>(ll.second.eval(stack));
             if ( copt && Ku.number <1)
             {
             R cc  =  GetAny<R>(ll.second.eval(stack));
             // cout << *(copt[il]) << " == " <<  cc << endl;
             if ( c != cc) {
             cerr << c << " != " << cc << " => ";
             cerr << "Sorry error in Optimization (m) add:  int2d(Th,optimize=0)(...)" << endl;
             ExecError("In Optimized version "); }
             }

             *pa += coef * c * w_i*w_j;
             }
             }

             }*/

        }
        else if(B)
        {
            R2 PA(B[0]),PB(B[1]);
            R2 A=T(PA),B=T(PB);
            R2 E(A,B);
            double le = sqrt((E,E));
            if(le > 1e-15)
            for (npi=0;npi<FIb.n;npi++) // loop on the integration point
            {

                pa =a;
                QuadratureFormular1dPoint pi( FIb[npi]);

                double coef = le*pi.a;
                double sa=pi.x,sb=1-sa;
                R2 Pt(PA*sa+PB*sb ); //
                Ku.BF(Dop,Pt,fu);
                // int label=-999999; // a passer en argument
                MeshPointStack(stack)->set(T(Pt),Pt,Ku,0,R2(E.y,-E.x)/le,0);
                if (classoptm) (*Op.optiexpK)(stack); // call optim version

                int il=0;
                for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                {  // attention la fonction test donne la ligne
                    //  et la fonction test est en second
                    BilinearOperator::K ll(*l);
                    //          pair<int,int> jj(ll.first.first),ii(ll.first.second);
                    long jcomp= ll.first.first.first,jop=ll.first.first.second;
                    long icomp= ll.first.second.first,iop=ll.first.second.second;

                    R c = copt ? *(copt[il]): GetAny<R>(ll.second.eval(stack));
                    if ( copt && Ku.number <1)
                    {
                        R cc  =  GetAny<R>(ll.second.eval(stack));
                        // cout << *(copt[il]) << " == " <<  cc << endl;
                        CheckErrorOptimisation(c,cc,"Sorry error in Optimization (n) add:  int2d(Th,optimize=0)(...)");
                        /*
                        if ( c != cc) {
                            cerr << c << " != " << cc << " => ";
                            cerr << "Sorry error in Optimization (n) add:  int2d(Th,optimize=0)(...)" << endl;
                            ExecError("In Optimized version "); }*/
                    }
                    c *= coef ;
                    long fi=Ku.dfcbegin(icomp);
                    long li=Ku.dfcend(icomp);
                    long fj=Ku.dfcbegin(jcomp);
                    long lj=Ku.dfcend(jcomp);

                    for ( i=fi;  i<li;   i++ )
                    for ( j=fj;  j<min(lj,i+1);  j++,pa++ ) //
                    {
                        R w_i =  fu(i,icomp,iop);
                        R w_j =  fu(j,jcomp,jop);

                        mat(i,j)  +=  c * w_i*w_j;
                        /*
                         if (Ku.Vh.Th(T) < 1 && npi < 1 && i < 1 && j < 1 )
                         cout <<" + " << c << " (" <<coef << " " << w_i << " " << w_j << " " << jj.first << " " << jj.second << ") " ;
                         */
                    }

                }
            }
        }
        else    // int on edge ie
        for (npi=0;npi<FIb.n;npi++) // loop on the integration point
        {

            pa =a;
            QuadratureFormular1dPoint pi( FIb[npi]);
            R2 E=T.Edge(ie);
            double le = sqrt((E,E));
            double coef = le*pi.a;
            double sa=pi.x,sb=1-sa;
            R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
            PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]);
            R2 Pt(PA*sa+PB*sb ); //
            Ku.BF(Dop,Pt,fu);
            // int label=-999999; // a passer en argument
            MeshPointStack(stack)->set(T(Pt),Pt,Ku,label,R2(E.y,-E.x)/le,ie);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version

            int il=0;
            for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
            {  // attention la fonction test donne la ligne
                //  et la fonction test est en second
                BilinearOperator::K ll(*l);
                //          pair<int,int> jj(ll.first.first),ii(ll.first.second);
                long jcomp= ll.first.first.first,jop=ll.first.first.second;
                long icomp= ll.first.second.first,iop=ll.first.second.second;

                R c = copt ? *(copt[il]): GetAny<R>(ll.second.eval(stack));
                if ( copt && Ku.number <1)
                {
                    R cc  =  GetAny<R>(ll.second.eval(stack));
                    // cout << *(copt[il]) << " == " <<  cc << endl;
                    CheckErrorOptimisation(c,cc,"Sorry error in Optimization (o) add:  int2d(Th,optimize=0)(...)");
                    /*
                    if ( c != cc) {
                        cerr << c << " != " << cc << " => ";
                        cerr << "Sorry error in Optimization (o) add:  int2d(Th,optimize=0)(...)" << endl;
                        ExecError("In Optimized version "); }*/
                }
                c *= coef ;
                long fi=Ku.dfcbegin(icomp);
                long li=Ku.dfcend(icomp);
                long fj=Ku.dfcbegin(jcomp);
                long lj=Ku.dfcend(jcomp);

                for ( i=fi;  i<li;   i++ )
                for ( j=fj;  j<min(lj,i+1);  j++,pa++ ) //
                {
                    R w_i =  fu(i,icomp,iop);
                    R w_j =  fu(j,jcomp,jop);

                    mat(i,j)  +=  c * w_i*w_j;
                    /*
                     if (Ku.Vh.Th(T) < 1 && npi < 1 && i < 1 && j < 1 )
                     cout <<" + " << c << " (" <<coef << " " << w_i << " " << w_j << " " << jj.first << " " << jj.second << ") " ;
                     */
                }

            }

            /*
             for ( i=0;  i<n;   i++ )
             // if ( onWhatIsEdge[ie][Ku.DFOnWhat(i)]) // generaly wrong FH dec 2003
             {
             RNM_ wi(fu(i,'.','.'));
             for ( j=0;  j<=i;   j++,pa++ )
             {
             RNM_ wj(fu(j,'.','.'));
             int il=0;
             for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
             // if (onWhatIsEdge[ie][Ku.DFOnWhat(j)]) // generaly wrong FH dec 2003
             {
             BilinearOperator::K ll(*l);
             pair<int,int> ii(ll.first.first),jj(ll.first.second);
             double w_i =  wi(ii.first,ii.second);
             double w_j =  wj(jj.first,jj.second);
             // R ccc = GetAny<R>(ll.second.eval(stack));
             R ccc = copt ? *(copt[il]): GetAny<R>(ll.second.eval(stack));
             if ( copt && Ku.number <1)
             {
             R cc  =  GetAny<R>(ll.second.eval(stack));
             if ( ccc != cc) {
             cerr << ccc << " != " << cc << ", xy = "<< T(Pt) << " => ";
             cerr << "Sorry error in Optimization (d)  add:  int2d(Th,optimize=0)(...)" << endl;
             ExecError("In Optimized version "); }
             }

             *pa += coef * ccc * w_i*w_j;
             }
             }
             } //else pa+= i+1;
             */
        }

        /*
         pa=a;
         if (Ku.Vh.Th(T) <=0 ) {
         cout <<endl  << " Triangle " << Ku.Vh.Th(T) << " =  "<<  T[0] << ", " << T[1] << ", " << T[2] << " " << nx << endl;
         for (int i=0;i<n;i++)
         {
         cout << setw(2) << i << setw(4) << mat.ni[i] << " :";
         for (int j=0;j<=i;j++)
         cout << setw(5)  << (*pa++) << " ";
         cout << endl;
         } }
         pa=a;
         for (int i=0;i<n;i++)
         cout << mat.ni[i] << " " ;
         for (int i=0;i<n;i++)
         for (int j=0;j<n;j++,pa++)
         if ( mat.ni[i]==150 && mat.nj[j] == 150)
         cout << "a_150,150 = "<< *pa ;
         cout << endl;
         */

        *MeshPointStack(stack) = mp;

    }



    // creating an instance of Element_Op MatriceElementaireSymetrique
    // case 3D volume
    template<class R>
    void  Element_Op(MatriceElementaireSymetrique<R,FESpace3> & mat,const FElement3 & Ku,double * p,int ie,int label, void * vstack,R3 *B)
    {
        //    ffassert(B==0);
        Stack stack=pvoid2Stack(vstack);
        typedef FESpace3 FESpace;
        typedef typename FESpace3::Mesh Mesh;
        typedef Mesh *pmesh ;
        typedef typename Mesh::Element Element;
        MeshPoint mp= *MeshPointStack(stack);
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        const Element & T  = Ku.T;

        const GQuadratureFormular<R3> & FI = mat.FIT;
        const GQuadratureFormular<R2> & FIb = mat.FIE;

        long npi;
        R *a=mat.a;
        R *pa=a;
        long i,j;
        long n= mat.n,m=mat.m,nx=n*(m+1)/2;
        long N= Ku.N;

        assert(mat.bilinearform);

        const Opera &Op(*mat.bilinearform);
        bool classoptm = copt && Op.optiexpK;
        // assert(  (copt !=0) ||  (Op.where_in_stack_opt.size() !=0) );
        int lastop;
        What_d Dop = Op.DiffOp(lastop);

        if (Ku.number<1  && verbosity/100 && verbosity % 10 == 2 )
        cout << "Element_Op S 3d: copt = " << copt << " " << classoptm << " lastop = "<< lastop << " Dop " << Dop << " opt: " << mat.optim << endl;
        assert(Op.MaxOp() <last_operatortype);

        RNMK_ fu(p,n,N,lastop); //  the value for basic fonction


        pa =a;
        for (i=0;i< nx;i++)
        *pa++ = 0.;

        if (ie<0)
        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            GQuadraturePoint<R3> pi(FI[npi]);
            double coef = T.mesure()*pi.a;
            //R3 Pt(pi);
            pa =a;
            Ku.BF(Dop,pi,fu);
            MeshPointStack(stack)->set(T(pi),pi,Ku);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            int il=0;
            for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
            {  // attention la fonction test donne la ligne
                //  et la fonction test est en second
                BilinearOperator::K ll(*l);
                //          pair<int,int> jj(ll.first.first),ii(ll.first.second);
                long jcomp= ll.first.first.first,jop=ll.first.first.second;
                long icomp= ll.first.second.first,iop=ll.first.second.second;


                R c = copt ? *(copt[il]): GetAny<R>(ll.second.eval(stack));
                if ( copt && Ku.number <1)
                {
                    R cc  =  GetAny<R>(ll.second.eval(stack));
                    // cout << *(copt[il]) << " == " <<  cc << endl;
                    CheckErrorOptimisation(c,cc,"Sorry error in Optimization (i) add:  int2d(Th,optimize=0)(...)");
                    /*
                    if ( c != cc) {
                        cerr << c << " != " << cc << " => ";
                        cerr << "Sorry error in Optimization (i) add:  int2d(Th,optimize=0)(...)" << endl;
                        ExecError("In Optimized version "); }*/
                }
                c *= coef ;
                long fi=Ku.dfcbegin(icomp);
                long li=Ku.dfcend(icomp);
                long fj=Ku.dfcbegin(jcomp);
                long lj=Ku.dfcend(jcomp);

                for ( i=fi;  i<li;   i++ )
                for ( j=fj;  j<min(lj,i+1);  j++,pa++ ) //
                {
                    R w_i =  fu(i,icomp,iop);
                    R w_j =  fu(j,jcomp,jop);

                    mat(i,j)  +=  c * w_i*w_j;

                    /*
                     if (Ku.Vh.Th(T) < 1 && npi < 1 && i < 1 && j < 1 )
                     cout <<" + " << c << " (" <<coef << " " << w_i << " " << w_j << " " << jj.first << " " << jj.second << ") " ;
                     */
                }

            }

        }
        else if(B)
        {  // int on leveset
            int np = ie-10; //= (B[0].x == B[3].x ) && (B[0].y == B[3].y ) && (B[0].z == B[3].z ) ? 3 : 4;
            if(verbosity>999) cout << "    Ass mat pleine /"<< np << endl;
            assert( np==3 || np==4);
            // XXXXXXX
            double epsmes3=T.mesure()*T.mesure()*1e-18;
            R3 PP[4];
            double l[3];
            for(int i=0; i< np; ++i)
            PP[i]= T(B[i]);

            for( int i =0; i+1 < np; i+=2)
            { // 0,1,, a and 2,3,0.
                int i0=i,i1=i+1,i2=(i+2)%np;
                R3 NN= R3(PP[i0],PP[i1])^R3(PP[i0],PP[i2]);
                double mes2 = (NN,NN);
                double mes = sqrt(mes2);

                if(mes2*mes <epsmes3) continue; //  too small
                NN /= mes;
                mes *= 0.5;
                if(verbosity>999)
                cout << " --int on leveset3d " << np << " " << mes << " " << i0<<i1<<i2 <<endl;
                double asum=0;
                for (npi=0;npi<FIb.n;npi++) // loop on the integration point
                {
                    GQuadraturePoint<R2>  pi( FIb[npi]);
                    // cout << " %% " << npi << " " << pi.a << " " << pi.x << " " << pi.y << endl;
                    asum+= pi.a;
                    pi.toBary(l);
                    R3 Pt( l[0]*B[i0]+l[1]*B[i1]+l[2]*B[i2]); //
                    double coef = mes*pi.a; // correction 0.5 050109 FH
                    Ku.BF(Dop,Pt,fu);
                    MeshPointStack(stack)->set(T(Pt),Pt,Ku,label,NN,ie);
                    if (classoptm) (*Op.optiexpK)(stack); // call optim version

                    pa=a;
                    int il=0;
                    for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                    {  // attention la fonction test donne la ligne
                        //  et la fonction test est en second
                        BilinearOperator::K ll(*l);
                        //          pair<int,int> jj(ll.first.first),ii(ll.first.second);
                        long jcomp= ll.first.first.first,jop=ll.first.first.second;
                        long icomp= ll.first.second.first,iop=ll.first.second.second;

                        R c = copt ? *(copt[il]): GetAny<R>(ll.second.eval(stack));
                        if ( copt && Ku.number <1)
                        {
                            R cc  =  GetAny<R>(ll.second.eval(stack));
                            // cout << *(copt[il]) << " == " <<  cc << endl;
                            CheckErrorOptimisation(c,cc,"Sorry error in Optimization (j) add:  int2d(Th,optimize=0)(...)");
                            /*
                            if ( c != cc) {
                                cerr << c << " != " << cc << " => ";
                                cerr << "Sorry error in Optimization (j) add:  int2d(Th,optimize=0)(...)" << endl;
                                ExecError("In Optimized version "); }*/
                        }
                        c *= coef ;
                        long fi=Ku.dfcbegin(icomp);
                        long li=Ku.dfcend(icomp);
                        long  fj=Ku.dfcbegin(jcomp);
                        long  lj=Ku.dfcend(jcomp);

                        for (long i=fi;  i<li;   i++ )
                        for (long j=fj;  j<min(lj,i+1);  j++,pa++ ) //
                        {
                            R w_i =  fu(i,icomp,iop);
                            R w_j =  fu(j,jcomp,jop);

                            mat(i,j)  +=  c * w_i*w_j;

                            /*
                             if (Ku.Vh.Th(T) < 1 && npi < 1 && i < 1 && j < 1 )
                             cout <<" + " << c << " (" <<coef << " " << w_i << " " << w_j << " " << jj.first << " " << jj.second << ") " ;
                             */
                        }

                    }



                }


            }

        }// end int level set ...
        else
        // int on edge ie
        for (npi=0;npi<FIb.n;npi++) // loop on the integration point
        {

            pa =a;
            GQuadraturePoint<R2> pi( FIb[npi]);
            R3 NN= T.N(ie);
            double mes=NN.norme();
            NN/=mes;
            mes *=0.5;
            double coef = mes*pi.a; // correction 0.5 050109 FH
            R3 Pt(T.PBord(ie,pi));
            Ku.BF(Dop,Pt,fu);
            // int label=-999999; // a passer en argument
            MeshPointStack(stack)->set(T(Pt),Pt,Ku,label,NN,ie);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            int il=0;
            for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
            {  // attention la fonction test donne la ligne
                //  et la fonction test est en second
                BilinearOperator::K ll(*l);
                //          pair<int,int> jj(ll.first.first),ii(ll.first.second);
                long jcomp= ll.first.first.first,jop=ll.first.first.second;
                long icomp= ll.first.second.first,iop=ll.first.second.second;

                R c = copt ? *(copt[il]): GetAny<R>(ll.second.eval(stack));
                if ( copt && Ku.number <1)
                {
                    R cc  =  GetAny<R>(ll.second.eval(stack));
                    // cout << *(copt[il]) << " == " <<  cc << endl;
                    CheckErrorOptimisation(c,cc,"Sorry error in Optimization (k) add:  int2d(Th,optimize=0)(...)");
                    /*
                    if ( c != cc) {
                        cerr << c << " != " << cc << " => ";
                        cerr << "Sorry error in Optimization (k) add:  int2d(Th,optimize=0)(...)" << endl;
                        ExecError("In Optimized version "); }*/
                }
                c *= coef ;
                long fi=Ku.dfcbegin(icomp);
                long li=Ku.dfcend(icomp);
                long  fj=Ku.dfcbegin(jcomp);
                long  lj=Ku.dfcend(jcomp);

                for (long i=fi;  i<li;   i++ )
                for (long j=fj;  j<min(lj,i+1);  j++,pa++ ) //
                {
                    R w_i =  fu(i,icomp,iop);
                    R w_j =  fu(j,jcomp,jop);

                    mat(i,j)  +=  c * w_i*w_j;

                    /*
                     if (Ku.Vh.Th(T) < 1 && npi < 1 && i < 1 && j < 1 )
                     cout <<" + " << c << " (" <<coef << " " << w_i << " " << w_j << " " << jj.first << " " << jj.second << ") " ;
                     */
                }

            }


        }


        pa=a;
        if (Ku.Vh.Th(T) <0 & verbosity>100) {
            cout <<endl  << " Tet " << Ku.Vh.Th(T) << " =  "<<  T << "  nx= " << nx << endl;
            for (int i=0;i<n;i++)
            {
                cout << setw(2) << i << setw(4) << mat.ni[i] << " :";
                for (int j=0;j<=i;j++)
                cout << setw(5)  << (*pa++) << " ";
                cout << endl;
            } }
        /*
         pa=a;
         for (int i=0;i<n;i++)
         cout << mat.ni[i] << " " ;
         for (int i=0;i<n;i++)
         for (int j=0;j<n;j++,pa++)
         if ( mat.ni[i]==150 && mat.nj[j] == 150)
         cout << "a_150,150 = "<< *pa ;
         cout << endl;
         */

        *MeshPointStack(stack) = mp;

    }

    // creating an instance of Element_Op with MatriceElementaireSymetrique
    // case 3D surface
    // xxxxxxxxxxxxxxxxx  modif a faire
    template<class R>
    void  Element_Op(MatriceElementaireSymetrique<R,FESpaceS> & mat,const FElementS & Ku,double * p,int ie,int label, void * vstack,R3 *B)
    {
        Stack stack=pvoid2Stack(vstack);
        MeshPoint mp= *MeshPointStack(stack);
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        const TriangleS & T  = Ku.T;
        //  const QuadratureFormular & FI = QuadratureFormular_T_2;
        //  const QuadratureFormular1d & FIb = QF_GaussLegendre2;
        const QuadratureFormular & FI = mat.FIT;
        const QuadratureFormular1d & FIb = mat.FIE;
        long npi;
        R *a=mat.a;
        R *pa=a;
        long i,j;
        long n= mat.n,m=mat.m,nx=n*(m+1)/2;
        long N= Ku.N;
        //long M=N;
        // bool show = Ku.Vh.Th(T)==0;
        //    char * xxx[] ={" u"," v"," p"," q"," r"};
        //char * xxxx[] ={" u'"," v'"," p'"," q'"," r'"};
        //char * yyy[] ={" ","_x ","_y "};


        throwassert(mat.bilinearform);

        const Opera &Op(*mat.bilinearform);
        bool classoptm = copt && Op.optiexpK;
        // assert(  (copt !=0) ||  (Op.where_in_stack_opt.size() !=0) );
        if (Ku.number<1  && verbosity/100 && verbosity % 10 == 2 )
            cout << "Element_Op S: copt = " << copt << " " << classoptm << " opt "<< mat.optim << endl;
        assert(Op.MaxOp() <last_operatortype);

        int lastop=0;
        What_d Dop = Op.DiffOp(lastop);

        // assert(lastop<=3);

        RNMK_ fu(p,n,N,lastop); //  the value for basic fonction

        pa =a;
        for (i=0;i< nx;i++)
            *pa++ = 0.;

        if (ie<0)
            for (npi=0;npi<FI.n;npi++) // loop on the integration point
            {
                QuadraturePoint pi(FI[npi]);
                double mes= B ? B->x :T.mesure();
                double coef = mes*pi.a;
                R2 Pt(pi);
                pa =a;
                Ku.BF(Dop,Pt,fu);
                MeshPointStack(stack)->set(T(pi),pi,Ku);
                if (classoptm) (*Op.optiexpK)(stack); // call optim version
                int il=0;
                for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                {  // attention la fonction test donne la ligne
                    //  et la fonction test est en second
                    BilinearOperator::K ll(*l);
                    //          pair<int,int> jj(ll.first.first),ii(ll.first.second);
                    long jcomp= ll.first.first.first,jop=ll.first.first.second;
                    long icomp= ll.first.second.first,iop=ll.first.second.second;

                    R c = copt ? *(copt[il]): GetAny<R>(ll.second.eval(stack));
                    if ( copt && Ku.number <1)
                    {
                        R cc  =  GetAny<R>(ll.second.eval(stack));
                        // cout << *(copt[il]) << " == " <<  cc << endl;
                        CheckErrorOptimisation(c,cc,"Sorry error in Optimization (l) add:  int2d(Th,optimize=0)(...)");
                        /*
                        if ( c != cc) {
                            cerr << c << " != " << cc << " => ";
                            cerr << "Sorry error in Optimization (l) add:  int2d(Th,optimize=0)(...)" << endl;
                            ExecError("In Optimized version "); }*/
                    }
                    c *= coef ;
                    long fi=Ku.dfcbegin(icomp);
                    long li=Ku.dfcend(icomp);
                    long fj=Ku.dfcbegin(jcomp);
                    long lj=Ku.dfcend(jcomp);
                    if (verbosity>10 && Ku.Vh.Th(T) < 1 && npi < 1)
                        cout << " ic "<< icomp << fi<< " "<< lj << " "<< " c "<< jcomp << " " <<fj << " "<< lj << endl;
                    for ( i=fi;  i<li;   i++ )
                        for ( j=fj;  j<min(lj,i+1);  j++ ) //
                        {
                            R w_i =  fu(i,icomp,iop);
                            R w_j =  fu(j,jcomp,jop);

                            mat(i,j)  +=  c * w_i*w_j;

                        }

                }

                /*
                 for ( i=0;  i<n;   i++ )
                 {
                 RNM_ wi(fu(i,'.','.'));
                 //    if (Ku.Vh.Th(T) < 1) cout << i <<" " <<Pt<< "wi =" << wi ;
                 for ( j=0;  j<=i;  j++,pa++ ) //
                 {

                 RNM_ wj(fu(j,'.','.'));
                 //   if (Ku.Vh.Th(T) < 1) cout << j <<" " <<Pt<< "wj =" << wj ;
                 int il=0;
                 for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                 {
                 const  BilinearOperator::K & ll(*l);
                 pair<int,int> ii(ll.first.first),jj(ll.first.second);
                 double w_i =  wi(ii.first,ii.second);
                 double w_j =  wj(jj.first,jj.second);

                 R c = copt ? *(copt[il]): GetAny<R>(ll.second.eval(stack));
                 if ( copt && Ku.number <1)
                 {
                 R cc  =  GetAny<R>(ll.second.eval(stack));
                 // cout << *(copt[il]) << " == " <<  cc << endl;
                 if ( c != cc) {
                 cerr << c << " != " << cc << " => ";
                 cerr << "Sorry error in Optimization (m) add:  int2d(Th,optimize=0)(...)" << endl;
                 ExecError("In Optimized version "); }
                 }

                 *pa += coef * c * w_i*w_j;
                 }
                 }

                 }*/

            }
        else if(B)

            ffassert(0);

        else    // int on edge ie
            for (npi=0;npi<FIb.n;npi++) // loop on the integration point
            {

                pa =a;
                QuadratureFormular1dPoint pi( FIb[npi]);
                R3 E=T.Edge(ie);
                double le = sqrt((E,E));
                double coef = le*pi.a;
                double sa=pi.x,sb=1-sa;
                R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
                PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]);
                R2 Pt(PA*sa+PB*sb ); //
                Ku.BF(Dop,Pt,fu);
                // int label=-999999; // a passer en argument
                MeshPointStack(stack)->set(T(Pt),Pt,Ku,label,R2(E.y,-E.x)/le,ie);
                if (classoptm) (*Op.optiexpK)(stack); // call optim version

                int il=0;
                for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                {  // attention la fonction test donne la ligne
                    //  et la fonction test est en second
                    BilinearOperator::K ll(*l);
                    //          pair<int,int> jj(ll.first.first),ii(ll.first.second);
                    long jcomp= ll.first.first.first,jop=ll.first.first.second;
                    long icomp= ll.first.second.first,iop=ll.first.second.second;

                    R c = copt ? *(copt[il]): GetAny<R>(ll.second.eval(stack));
                    if ( copt && Ku.number <1)
                    {
                        R cc  =  GetAny<R>(ll.second.eval(stack));
                        // cout << *(copt[il]) << " == " <<  cc << endl;
                        CheckErrorOptimisation(c,cc,"Sorry error in Optimization (o) add:  int2d(Th,optimize=0)(...)");
                        /*
                        if ( c != cc) {
                            cerr << c << " != " << cc << " => ";
                            cerr << "Sorry error in Optimization (o) add:  int2d(Th,optimize=0)(...)" << endl;
                            ExecError("In Optimized version "); }*/
                    }
                    c *= coef ;
                    long fi=Ku.dfcbegin(icomp);
                    long li=Ku.dfcend(icomp);
                    long fj=Ku.dfcbegin(jcomp);
                    long lj=Ku.dfcend(jcomp);

                    for ( i=fi;  i<li;   i++ )
                        for ( j=fj;  j<min(lj,i+1);  j++,pa++ ) //
                        {
                            R w_i =  fu(i,icomp,iop);
                            R w_j =  fu(j,jcomp,jop);

                            mat(i,j)  +=  c * w_i*w_j;
                            /*
                             if (Ku.Vh.Th(T) < 1 && npi < 1 && i < 1 && j < 1 )
                             cout <<" + " << c << " (" <<coef << " " << w_i << " " << w_j << " " << jj.first << " " << jj.second << ") " ;
                             */
                        }

                }

                /*
                 for ( i=0;  i<n;   i++ )
                 // if ( onWhatIsEdge[ie][Ku.DFOnWhat(i)]) // generaly wrong FH dec 2003
                 {
                 RNM_ wi(fu(i,'.','.'));
                 for ( j=0;  j<=i;   j++,pa++ )
                 {
                 RNM_ wj(fu(j,'.','.'));
                 int il=0;
                 for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                 // if (onWhatIsEdge[ie][Ku.DFOnWhat(j)]) // generaly wrong FH dec 2003
                 {
                 BilinearOperator::K ll(*l);
                 pair<int,int> ii(ll.first.first),jj(ll.first.second);
                 double w_i =  wi(ii.first,ii.second);
                 double w_j =  wj(jj.first,jj.second);
                 // R ccc = GetAny<R>(ll.second.eval(stack));
                 R ccc = copt ? *(copt[il]): GetAny<R>(ll.second.eval(stack));
                 if ( copt && Ku.number <1)
                 {
                 R cc  =  GetAny<R>(ll.second.eval(stack));
                 if ( ccc != cc) {
                 cerr << ccc << " != " << cc << ", xy = "<< T(Pt) << " => ";
                 cerr << "Sorry error in Optimization (d)  add:  int2d(Th,optimize=0)(...)" << endl;
                 ExecError("In Optimized version "); }
                 }

                 *pa += coef * ccc * w_i*w_j;
                 }
                 }
                 } //else pa+= i+1;
                 */
            }

        /*
         pa=a;
         if (Ku.Vh.Th(T) <=0 ) {
         cout <<endl  << " Triangle " << Ku.Vh.Th(T) << " =  "<<  T[0] << ", " << T[1] << ", " << T[2] << " " << nx << endl;
         for (int i=0;i<n;i++)
         {
         cout << setw(2) << i << setw(4) << mat.ni[i] << " :";
         for (int j=0;j<=i;j++)
         cout << setw(5)  << (*pa++) << " ";
         cout << endl;
         } }
         pa=a;
         for (int i=0;i<n;i++)
         cout << mat.ni[i] << " " ;
         for (int i=0;i<n;i++)
         for (int j=0;j<n;j++,pa++)
         if ( mat.ni[i]==150 && mat.nj[j] == 150)
         cout << "a_150,150 = "<< *pa ;
         cout << endl;
         */

        *MeshPointStack(stack) = mp;


    }

    // creating an instance of Element_Op with MatriceElementaireSymetrique
    // case 3D surface
    // xxxxxxxxxxxxxxxxx  modif a faire
    template<class R>
    void  Element_Op(MatriceElementaireSymetrique<R,FESpaceL> & mat,const FElementL & Ku,double * p,int ie,int label, void * vstack,R3 *B)
    {
        ffassert(0);
    }


   //////////////////////////////////
   // Element_rhs
   //////////////////////////////////

    // #pragma optimization_level 0
    // creating an instance of Element_rhs
    // case 2d
    template<class R>
    void  Element_rhs(const FElement & Kv,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular & FI = QuadratureFormular_T_2,int optim=1)
    {
        Stack stack=pvoid2Stack(vstack);
        MeshPoint mp=*MeshPointStack(stack) ;
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        const Triangle & T  = Kv.T;
        //  const QuadratureFormular & FI = QuadratureFormular_T_2;
        //  const QuadratureFormular & FI = QuadratureFormular_T_2;
        long npi;
        long i,n=Kv.NbDoF(),N=Kv.N;

        //  bool show = Kv.Vh.Th(T)==0;
        //  char * xxx[] ={" u"," v,"," p"," q"," r"};
        // char * xxxx[] ={" u'"," v',"," p'"," q'"," r'"};
        // char * yyy[] ={" ","_x ","_y "};

        bool classoptm = copt && Op.optiexpK;
        // assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
        if (Kv.number<1  && verbosity/100 && verbosity % 10 == 2)
        cout << "Element_rhs S0: copt = " << copt << " " << classoptm << " opt " << optim << endl;


        KN<bool> Dop(last_operatortype);
        Op.DiffOp(Dop);
        int lastop=1+Dop.last([](bool x){return x;});
        assert(Op.MaxOp() <last_operatortype);

        //  assert(lastop<=3);


        RNMK_ fu(p,n,N,lastop); //  the value for basic fonction

        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            QuadraturePoint pi(FI[npi]);
            double coef = T.area*pi.a;
            R2 Pt(pi);
            Kv.BF(Dop,Pt,fu);
            MeshPointStack(stack)->set(T(Pt),Pt,Kv);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            for ( i=0;  i<n;   i++ )
            {
                RNM_ wi(fu(i,'.','.'));
                int il=0;
                for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                {
                    LOperaD::K ll(*l);
                    pair<int,int> ii(ll.first);
                    double w_i =  wi(ii.first,ii.second);
                    //copt=0;
                    R c = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack)); //GetAny<double>(ll.second.eval(stack));
                    if ( copt && ( optim==1) && Kv.number <1)
                    {
                        R cc  =  GetAny<R>(ll.second.eval(stack));
                        CheckErrorOptimisation(c,cc, "Sorry error in Optimization add:  (p) int2d(Th,optimize=0)(...)" );
                        /*
                        if ( c != cc) {
                            cerr << c << " != " << cc << " => ";
                            cerr << "Sorry error in Optimization add:  (p) int2d(Th,optimize=0)(...)" << endl;
                            ExecError("In Optimized version "); }*/
                    }
                    //if (Kv.number<5) cout << il<< " " << i << "  c== " <<  c << endl;
                    R a = coef * c * w_i;
                    B[Kv(i)] += a;
                }
            }


        }
        *MeshPointStack(stack) = mp;


    }

    // creating an instance of Element_rhs
    // case 3D volume
    template<class R>
    void  Element_rhs(const FElement3 & Kv,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const GQuadratureFormular<R3> & FI = QuadratureFormular_Tet_2,int optim=1)
    {
        Stack stack=pvoid2Stack(vstack);
        typedef  FElement3::Element Element;
        MeshPoint mp=*MeshPointStack(stack) ;
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        const Element & T  = Kv.T;
        //  const QuadratureFormular & FI = QuadratureFormular_T_2;
        //  const QuadratureFormular & FI = QuadratureFormular_T_2;
        long npi;
        long i,n=Kv.NbDoF(),N=Kv.N;


        bool classoptm = copt && Op.optiexpK;
        // assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
        if (Kv.number<1  && verbosity/100 && verbosity % 10 == 2)
        cout << "Element_rhs S0: copt = " << copt << " " << classoptm << " opt: " << optim << endl;


        int lastop;
        What_d Dop = Op.DiffOp(lastop);
        assert(Op.MaxOp() <last_operatortype);

        //  assert(lastop<=3);


        RNMK_ fu(p,n,N,lastop); //  the value for basic fonction

        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            GQuadraturePoint<R3> pi(FI[npi]);
            double coef = T.mesure()*pi.a;
            R3 Pt(pi);
            Kv.BF(Dop,Pt,fu);
            MeshPointStack(stack)->set(T(Pt),Pt,Kv);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            for ( i=0;  i<n;   i++ )
            {
                RNM_ wi(fu(i,'.','.'));
                int il=0;
                for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                {
                    LOperaD::K ll(*l);
                    pair<int,int> ii(ll.first);
                    double w_i =  wi(ii.first,ii.second);
                    //copt=0;
                    R c = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack)); //GetAny<double>(ll.second.eval(stack));
                    if ( copt && ( optim==1) && Kv.number <1)
                    {
                        R cc  =  GetAny<R>(ll.second.eval(stack));
                        CheckErrorOptimisation(c,cc,"Sorry error in Optimization (q) add:  int2d(Th,optimize=0)(...)");
                        /*
                        if ( c != cc) {
                            cerr << c << " != " << cc << " => ";
                            cerr << "Sorry error in Optimization (q) add:  int2d(Th,optimize=0)(...)" << endl;
                            ExecError("In Optimized version "); }*/
                    }
                    //if (Kv.number<5) cout << il<< " " << i << "  c== " <<  c << endl;
                    R a = coef * c * w_i;
                    B[Kv(i)] += a;
                }
            }


        }
        *MeshPointStack(stack) = mp;

    }


    // creating an instance of Element_rhs
    // case 3D surface
    template<class R>
    void  Element_rhs(const FElementS & Kv,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular & FI = QuadratureFormular_T_2,int optim=1)
    {
        Stack stack=pvoid2Stack(vstack);
        MeshPoint mp=*MeshPointStack(stack) ;
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        const TriangleS & T  = Kv.T;
        //  const QuadratureFormular & FI = QuadratureFormular_T_2;
        //  const QuadratureFormular & FI = QuadratureFormular_T_2;
        long npi;
        long i,n=Kv.NbDoF(),N=Kv.N;

        //  bool show = Kv.Vh.Th(T)==0;
        //  char * xxx[] ={" u"," v,"," p"," q"," r"};
        // char * xxxx[] ={" u'"," v',"," p'"," q'"," r'"};
        // char * yyy[] ={" ","_x ","_y "};

        bool classoptm = copt && Op.optiexpK;
        // assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
        if (Kv.number<1  && verbosity/100 && verbosity % 10 == 2)
            cout << "Element_rhs S0: copt = " << copt << " " << classoptm << " opt " << optim << endl;

        int lastop=0;
        What_d Dop = Op.DiffOp(lastop);

        //  assert(lastop<=3);


        RNMK_ fu(p,n,N,lastop); //  the value for basic fonction

        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            QuadraturePoint pi(FI[npi]);
            double coef = T.mesure()*pi.a;
            R2 Pt(pi);
            Kv.BF(Dop,Pt,fu);
            MeshPointStack(stack)->set(T(Pt),Pt,Kv);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            for ( i=0;  i<n;   i++ )
            {
                RNM_ wi(fu(i,'.','.'));
                int il=0;
                for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                {
                    LOperaD::K ll(*l);
                    pair<int,int> ii(ll.first);
                    double w_i =  wi(ii.first,ii.second);
                    //copt=0;
                    R c = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack)); //GetAny<double>(ll.second.eval(stack));
                    if ( copt && ( optim==1) && Kv.number <1)
                    {
                        R cc  =  GetAny<R>(ll.second.eval(stack));
                        CheckErrorOptimisation(c,cc,"Sorry error in Optimization Element_OpVF2d  (b) add:  int2d(Th,optimize=0)(...)");
                        /*
                        if ( c != cc) {
                            cerr << c << " != " << cc << " => ";
                            cerr << "Sorry error in Optimization add:  (p) int2d(Th,optimize=0)(...)" << endl;
                            ExecError("In Optimized version "); }*/
                    }
                    //if (Kv.number<5) cout << il<< " " << i << "  c== " <<  c << endl;
                    R a = coef * c * w_i;
                    B[Kv(i)] += a;
                }
            }


        }
        *MeshPointStack(stack) = mp;
    }

    // case 3D curve YYYY

template<class R>
void  Element_rhsVF(const FElementL & Kv,const FElementL & KKv,int ie,int iie,int label,const LOperaD &Op,double * p,int *ip,void  * bstack,KN_<R> & B,
                    int optim=1)
{
    
    
    typedef typename FElementL::E Element;
    pair_stack_double * bs=static_cast<pair_stack_double *>(bstack);
    Stack stack= bs->first;
    double binside = *bs->second; // truc FH pour fluide de grad2 (decentrage bizard)
    assert(ie>=0 && ie < 2); //  int verter
    MeshPoint mp= *MeshPointStack(stack);
    R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
    
    const Element & T  = Kv.T;
    R3 NNt=T.TangenteUnitaire();
    double s = ie;
    double ss = iie;
    GQuadraturePoint<R1> pi(1.,s);
    R mes = 1.;
    R3 NN=NNt;
    if(ie==0) NN=-NN;
    R coef = 1.;
  
    int nTonEdge =  &Kv == &KKv ? 1 : 2;
    double cmean = 1./nTonEdge;
    
    throwassert(&T == &Kv.T);
    long npi=1;
    long i,j;
    long N= Kv.N;
    
    long nv=Kv.NbDoF();
    long nnv=KKv.NbDoF();
    assert(nv == nnv) ;
    
    int lp =nv*2;
    KN_<int> pp(ip,lp),pk(ip+lp,lp),pkk(ip+2*lp,lp);
    int n = BuildMEK_KK(lp,pp,pk,pkk,&Kv,&KKv);
    
    
    bool classoptm = copt && Op.optiexpK;
    //  if (Ku.number<1 && verbosity/100 && verbosity % 10 == 2)
    if (Kv.number<1 && ( verbosity > 1 ) )
        cout << "Element_rhsVF 0d P: copt = " << copt << " " << classoptm << " binside (For FH) =" << binside  << endl;
    
    
    int lastop;
    lastop = 0;
    What_d Dop = Op.DiffOp(lastop);
    //assert(lastop<=3);
    int lffv = nv*N*last_operatortype;
    int loffset =  0 ;
    
    RNMK_ fv(p,nv,N,lastop); //  the value for basic fonction in K
    RNMK_ ffv(p + lffv ,nnv,N,lastop); //  the value for basic fonction in KK
    
    
    
    
 
    
    R1 Pt=s; //
    R1 PP_t=ss;
    assert (!binside);
    Kv.BF(Dop,Pt,fv);
    KKv.BF(Dop,PP_t,ffv);
    MeshPointStack(stack)->set(T(Pt),Pt,Kv,label,NN,NNt,ie);//   Axel
    if (classoptm) (*Op.optiexpK)(stack); // call optim version
    
    
    for ( i=0;  i<nv;   i++ )
    {
        
        int ik= pk[i];
        int ikk=pkk[i];
        int dofik=ik>=0? Kv(ik):-1;
        int dofikk=ikk>=0? KKv(ikk):-1;
        
        RNM_ wi(fv(Max(ik,0),'.','.'));
        RNM_ wwi(ffv(Max(ikk,0),'.','.'));
        
        
        int il=0;
        if(dofik >=0 || dofikk>=0 )
            for (LinearOperatorD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
        {
            LOperaD::K ll(*l);
            pair<int,int> ii(ll.first);
            int iis = ii.second;
            int iicase  = iis / last_operatortype;
            iis %= last_operatortype;
            double w_i=0,ww_i=0;
            if(ik>=0) w_i =   wi(ii.first,iis );
            if( iicase>0 )
            {
                if( ikk>=0) ww_i =  wwi(ii.first,iis );
                if       (iicase==Code_Jump)      w_i = -w_i; ///(w_i = ww_i-w_i); // jump
                else  if (iicase==Code_Mean)      ww_i=w_i = cmean*  (w_i + ww_i ); // average
                else  if (iicase==Code_OtherSide) std::swap(w_i,ww_i);  // valeur de autre cote
                else ffassert(0);
            }
            R c =copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
            // FFCS - removing what is probably a small glitch
            if ( copt && ( optim==1) && Kv.number<1)
            {
                R cc  =  GetAny<R>(ll.second.eval(stack));
                if ( c != cc) {
                    cerr << c << " =! " << cc << endl;
                    cerr << "Sorry error in Optimization (x) add:  int0d(Th,optimize=0)(...)" << endl;
                    ExecError("In Optimized version "); }
            }
            
            
            if(dofik>=0) B[dofik] += coef * c * w_i;
            if(dofikk>=0) B[dofikk] += coef * c * ww_i;
            
        }
    }
    
    
    
    
    *MeshPointStack(stack) = mp;
    
}
    template<class R>
    void  Element_rhs(const FElementL & Kv,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const GQuadratureFormular<R1> & FI = QF_GaussLegendre2,int optim=1)
    {
        Stack stack=pvoid2Stack(vstack);
        MeshPoint mp=*MeshPointStack(stack) ;
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        const EdgeL & T  = Kv.T;
        long npi;
        long i,n=Kv.NbDoF(),N=Kv.N;

        bool classoptm = copt && Op.optiexpK;
        // assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
        if (Kv.number<1  && verbosity/100 && verbosity % 10 == 2)
            cout << "Element_rhs S0: copt = " << copt << " " << classoptm << " opt " << optim << endl;

        int lastop=0;
        What_d Dop = Op.DiffOp(lastop);

        RNMK_ fu(p,n,N,lastop); //  the value for basic fonction

        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            GQuadraturePoint<R1> pi(FI[npi]);
            double coef = T.mesure()*pi.a;
            R1 Pt(pi);
            Kv.BF(Dop,Pt,fu);
            MeshPointStack(stack)->set(T(Pt),Pt,Kv);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            for ( i=0;  i<n;   i++ )
            {
                RNM_ wi(fu(i,'.','.'));
                int il=0;
                for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                {
                    LOperaD::K ll(*l);
                    pair<int,int> ii(ll.first);
                    double w_i =  wi(ii.first,ii.second);
                    //copt=0;
                    R c = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack)); //GetAny<double>(ll.second.eval(stack));
                    if ( copt && ( optim==1) && Kv.number <1)
                    {
                        R cc  =  GetAny<R>(ll.second.eval(stack));
                        CheckErrorOptimisation(c,cc,"Sorry error in Optimization Element_OpVF3dcurve  (b) add:  int2d(Th,optimize=0)(...)");
                    }
                    R a = coef * c * w_i;
                    B[Kv(i)] += a;
                }
            }


        }
        *MeshPointStack(stack) = mp;
    }

    // #pragma optimization_level 0

    // creating an instance of Element_rhs
    // case 2d
    template<class R>
    void  Element_rhs(const  Mesh & ThI,const Triangle & KI,
                      const FESpace & Vh,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular & FI = QuadratureFormular_T_2,int optim=1)
    {
        Stack stack=pvoid2Stack(vstack);
        MeshPoint mp=*MeshPointStack(stack) ;
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        //    int maxd = Op.MaxOp();
        //    assert(maxd<last_operatortype);
        const Triangle * Kp=0;

        bool classoptm = copt && Op.optiexpK;
        // assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
        if (ThI(KI)<1 && verbosity/100 && verbosity % 10 == 2)

        cout << "Element_rhs 3: copt = " << copt << " " << classoptm <<" opt " << optim<< endl;

        KN<bool> Dop(last_operatortype);
        Op.DiffOp(Dop);
        int lastop=1+Dop.last([](bool x){return x;});
        assert(Op.MaxOp() <last_operatortype);

        // assert(lastop<=3);

        for (long npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            QuadraturePoint pi(FI[npi]);
            R2 PI(KI(pi));
            double coef = KI.area*pi.a;
            MeshPointStack(stack)->set(ThI,PI,pi,KI,KI.lab);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            bool outside;
            R2 Pt;
            const Triangle & K  = *Vh.Th.Find(PI,Pt,outside,Kp);
            if ( ! outside)
            {
                const  FElement  Kv= Vh[K];
                long i,n=Kv.NbDoF(),N=Kv.N;
                RNMK_ fu(p,n,N,lastop); //  the value for basic fonction
                Kv.BF(Dop,Pt,fu);

                for ( i=0;  i<n;   i++ )
                {
                    RNM_ wi(fu(i,'.','.'));
                    int il=0;
                    for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                    {
                        LOperaD::K ll(*l);
                        pair<int,int> ii(ll.first);

                        double w_i =  wi(ii.first,ii.second);

                        R c = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));;//GetAny<double>(ll.second.eval(stack));
                        if ( copt && ThI(KI) <1)
                        {
                            R cc  =  GetAny<R>(ll.second.eval(stack));
                            CheckErrorOptimisation(c,cc,"Sorry error in Optimization (s) add:  int2d(Th,optimize=0)(...)");
                            /*
                            if ( c != cc) {
                                cerr << c << " != " << cc << " => ";
                                cerr << "Sorry error in Optimization (s) add:  int2d(Th,optimize=0)(...)" << endl;
                                ExecError("In Optimized version "); }*/
                        }

                        R a = coef * c * w_i;
                        B[Kv(i)] += a;
                    }
                }
            }
            Kp = & K;
        }
        *MeshPointStack(stack) = mp;


    }


    // creating an instance of Element_rhs
    // case 3D volume
    template<class R>
    void  Element_rhs(const  Mesh3 & ThI,const Mesh3::Element & KI,
                      const FESpace3 & Vh,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const GQuadratureFormular<R3> & FI,int optim)
    {
        Stack stack=pvoid2Stack(vstack);
        // AFAIRE("Element_rhs 3d diff meshes");
        static int count=0;
        if(count++<1)
        {
            cout << "Warning:  Element_rhs 3 3d diff meshes in test (FH) " << endl;
            cout << "--------------------------------------------------- " << endl;
        }
        MeshPoint mp=*MeshPointStack(stack) ;
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        //    int maxd = Op.MaxOp();
        //    assert(maxd<last_operatortype);
        const Tet * Kp=0;

        bool classoptm = copt && Op.optiexpK;
        // assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
        if (ThI(KI)<1 && verbosity/100 && verbosity % 10 == 2)

        cout << "Element_rhs 3d  3: copt = " << copt << " " << classoptm <<" opt " <<optim << endl;

        assert(Op.MaxOp() <last_operatortype);
        //
        int lastop=0;
        lastop = 0;
        What_d Dop = Op.DiffOp(lastop);
        assert(Op.MaxOp() <last_operatortype);

        // assert(lastop<=3);

        for (long npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            GQuadraturePoint<R3> pi(FI[npi]);
            R3 PI(KI(pi));
            double coef = KI.mesure()*pi.a;
            MeshPointStack(stack)->set(ThI,PI,pi,KI,KI.lab);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            bool outside;
            R3 Pt;
            const Tet & K  = *Vh.Th.Find(PI,Pt,outside,Kp);
            if ( ! outside)
            {
                const  FElement3  Kv= Vh[K];
                long i,n=Kv.NbDoF(),N=Kv.N;
                RNMK_ fu(p,n,N,lastop); //  the value for basic fonction
                Kv.BF(Dop,Pt,fu);

                for ( i=0;  i<n;   i++ )
                {
                    RNM_ wi(fu(i,'.','.'));
                    int il=0;
                    for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                    {
                        LOperaD::K ll(*l);
                        pair<int,int> ii(ll.first);

                        double w_i =  wi(ii.first,ii.second);

                        R c = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));;//GetAny<double>(ll.second.eval(stack));
                        if ( copt && ThI(KI) <1 && optim==1)
                        {
                            R cc  =  GetAny<R>(ll.second.eval(stack));
                            if ( c != cc) {
                                cerr << c << " != " << cc << " => ";
                                cerr << "Sorry error in Optimization (r) add:  int2d(Th,optimize=0)(...)" << endl;
                                ExecError("In Optimized version "); }
                        }

                        R a = coef * c * w_i;
                        B[Kv(i)] += a;
                    }
                }
            }
            Kp = & K;
        }
        *MeshPointStack(stack) = mp;


    }

    // creating an instance of Element_rhs
    // case 3D surface
    template<class R>
    void  Element_rhs(const  MeshS & ThI,const TriangleS & KI,
                      const FESpaceS & Vh,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular & FI = QuadratureFormular_T_2,int optim=1)
    {
        ffassert(0);
    }

    // creating an instance of Element_rhs
    // case 3D curve
    template<class R>
    void  Element_rhs(const  MeshL & ThI,const EdgeL & KI,
                      const FESpaceL & Vh,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const GQuadratureFormular<R1> & FI = QF_GaussLegendre2,int optim=1)
    {
        ffassert(0);
    }

    // creating an instance of Element_rhs
    // case 2d
    template<class R>
    void  Element_rhs(Expression const * const mapt,const  Mesh & ThI,const Triangle & KI,
                      const FESpace & Vh,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular & FI = QuadratureFormular_T_2,int optim=1)
    {

        Stack stack=pvoid2Stack(vstack);
        MeshPoint mp=*MeshPointStack(stack) ;
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        //    int maxd = Op.MaxOp();
        //    assert(maxd<last_operatortype);
        const Triangle * Kp=0;

        bool classoptm = copt && Op.optiexpK;
        // assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
        if (ThI(KI)<1 && verbosity/100 && verbosity % 10 == 2)

        cout << "Element_rhs 3: copt = " << copt << " " << classoptm <<" opt " << optim<< endl;

        KN<bool> Dop(last_operatortype);
        Op.DiffOp(Dop);
        int lastop=1+Dop.last([](bool x){return x;});
        assert(Op.MaxOp() <last_operatortype);

        // assert(lastop<=3);

        for (long npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            QuadraturePoint pi(FI[npi]);
            R2 PIo(KI(pi)),PI(PIo);
            double coef = KI.area*pi.a;
            MeshPointStack(stack)->set(ThI,PI,pi,KI,KI.lab);
            if(mapt)
            { // move poit
                PI= R2( GetAny<double>((*mapt[0])(vstack)), GetAny<double>((*mapt[1])(vstack)));
                if(verbosity>9999) cout << "  Element_rhs mapt =" << PIo << " -> " <<PI << endl;
            }
            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            bool outside;
            R2 Pt;
            const Triangle & K  = *Vh.Th.Find(PI,Pt,outside,Kp);
            if ( ! outside)
            {
                const  FElement  Kv= Vh[K];
                long i,n=Kv.NbDoF(),N=Kv.N;
                RNMK_ fu(p,n,N,lastop); //  the value for basic fonction
                Kv.BF(Dop,Pt,fu);

                for ( i=0;  i<n;   i++ )
                {
                    RNM_ wi(fu(i,'.','.'));
                    int il=0;
                    for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                    {
                        LOperaD::K ll(*l);
                        pair<int,int> ii(ll.first);

                        double w_i =  wi(ii.first,ii.second);

                        R c = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));;//GetAny<double>(ll.second.eval(stack));
                        if ( copt && ThI(KI) <1)
                        {
                            R cc  =  GetAny<R>(ll.second.eval(stack));
                            if ( c != cc) {
                                cerr << c << " != " << cc << " => ";
                                cerr << "Sorry error in Optimization (s) add:  int2d(Th,optimize=0)(...)" << endl;
                                ExecError("In Optimized version "); }
                        }

                        R a = coef * c * w_i;
                        B[Kv(i)] += a;
                    }
                }
            }
            Kp = & K;
        }
        *MeshPointStack(stack) = mp;


    }
    // creating an instance of Element_rhs
    // case 3D volume
    template<class R>
    void  Element_rhs(const FElement3 & Kv,int ie,int label,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular & FI ,bool alledges=false,int optim=1)
    {
        //   AFAIRE("Element_rhs on border");
        Stack stack=pvoid2Stack(vstack);
        typedef  FElement3::Element Element;
        MeshPoint mp=*MeshPointStack(stack) ;
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        const Element & T  = Kv.T;
        long npi;
        long i,n=Kv.NbDoF(),N=Kv.N;

        bool classoptm = copt && Op.optiexpK;
        // assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
        if (Kv.number<1 && verbosity/100 && verbosity % 10 == 2)
        cout << "Element_rhs 3d S: copt = " << copt << " " << classoptm <<" opt " << optim << endl;
        int lastop;
        What_d Dop = Op.DiffOp(lastop);

        assert(Op.MaxOp() <last_operatortype);
        // assert(lastop<=3);

        RNMK_ fu(p,n,N,lastop); //  the value for basic fonction

        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            GQuadraturePoint<R2> pi( FI[npi]);
            R3 NN=T.N(ie);
            double le= NN.norme();
            NN /= le;
            double coef = le*pi.a*0.5;// correction 050109 FH
            R3 Pt(T.PBord(ie,pi));
            //
            Kv.BF(Dop,Pt,fu);
            MeshPointStack(stack)->set(T(Pt),Pt,Kv,label,NN,ie);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version

            for ( i=0;  i<n;   i++ )
            // if (alledges || onWhatIsEdge[ie][Kv.DFOnWhat(i)]) // bofbof faux si il y a des derives ..
            {
                RNM_ wi(fu(i,'.','.'));
                int il=0;
                for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                {
                    LOperaD::K ll(*l);
                    pair<int,int> ii(ll.first);
                    double w_i =  wi(ii.first,ii.second);
                    R c =copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                    // FFCS - removing what is probably a small glitch
                    if ( copt && ( optim==1) && Kv.number<1)
                    {
                        R cc  =  GetAny<R>(ll.second.eval(stack));
                        if ( c != cc) {
                            cerr << c << " =! " << cc << endl;
                            cerr << "Sorry error in Optimization (t) add:  int2d(Th,optimize=0)(...)" << endl;
                            ExecError("In Optimized version "); }
                    }


                    //= GetAny<double>(ll.second.eval(stack));
                    if(verbosity>100000)
                        cout << "Element_rhs 3d S: " <<coef << " " << c << " "<< w_i << " = "<< coef * c * w_i << " += " << Kv(i) <<endl;
                    B[Kv(i)] += coef * c * w_i;
                }
            }


        }
        *MeshPointStack(stack) = mp;

    }

    // creating an instance of Element_rhs
    // case 3D surface
    template<class R>
    void  Element_rhs(Expression const * const mapt,const  MeshS & ThI,const TriangleS & KI,
                      const FESpaceS & Vh,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular & FI = QuadratureFormular_T_2,int optim=1)
    {
        ffassert(0);
    }
    // creating an instance of Element_rhs
    // case 3D curve
    template<class R>
    void  Element_rhs(Expression const * const mapt,const  MeshL & ThI,const EdgeL & KI,
                      const FESpaceL & Vh,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular & FI = QuadratureFormular_T_2,int optim=1)
    {
        ffassert(0);
    }


    // creating an instance of Element_rhs
    // case 2d
    template<class R>
    void  Element_rhs(const FElement & Kv,int ie,int label,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular1d & FI = QF_GaussLegendre2,bool alledges=false,int optim=1)
    {
        Stack stack=pvoid2Stack(vstack);
        MeshPoint mp=*MeshPointStack(stack) ;
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        const Triangle & T  = Kv.T;
        // const QuadratureFormular1d & FI = QF_GaussLegendre2;
        long npi;
        long i,n=Kv.NbDoF(),N=Kv.N;

        //  bool show = Kv.Vh.Th(T)==0;
        // char * xxx[] ={" u"," v,"," p"," q"," r"};
        // char * xxxx[] ={" u'"," v',"," p'"," q'"," r'"};
        // char * yyy[] ={" ","_x ","_y "};

        bool classoptm = copt && Op.optiexpK;
        // assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
        if (Kv.number<1 && verbosity/100 && verbosity % 10 == 2)
        cout << "Element_rhs S: copt = " << copt << " " << classoptm << "opt " << optim << endl;
        KN<bool> Dop(last_operatortype);
        Op.DiffOp(Dop);
        int lastop=1+Dop.last([](bool x){return x;});
        assert(Op.MaxOp() <last_operatortype);
        // assert(lastop<=3);

        RNMK_ fu(p,n,N,lastop); //  the value for basic fonction

        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            QuadratureFormular1dPoint pi( FI[npi]);
            R2 E=T.Edge(ie);
            double le = sqrt((E,E));
            double coef = le*pi.a;
            double sa=pi.x,sb=1-sa;
            R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
            PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]);
            R2 Pt(PA*sa+PB*sb ); //
            Kv.BF(Dop,Pt,fu);
            MeshPointStack(stack)->set(T(Pt),Pt,Kv,label,R2(E.y,-E.x)/le,ie);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version

            for ( i=0;  i<n;   i++ )
            // if (alledges || onWhatIsEdge[ie][Kv.DFOnWhat(i)]) // bofbof faux si il y a des derives ..
            {
                RNM_ wi(fu(i,'.','.'));
                int il=0;
                for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                {
                    LOperaD::K ll(*l);
                    pair<int,int> ii(ll.first);
                    double w_i =  wi(ii.first,ii.second);
                    R c =copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                    // FFCS - removing what is probably a small glitch
                    if ( copt && ( optim==1) && Kv.number<1)
                    {
                        R cc  =  GetAny<R>(ll.second.eval(stack));
                        if ( c != cc) {
                            cerr << c << " =! " << cc << endl;
                            cerr << "Sorry error in Optimization (v) add:  int2d(Th,optimize=0)(...)" << endl;
                            ExecError("In Optimized version "); }
                    }


                    //= GetAny<double>(ll.second.eval(stack));

                    B[Kv(i)] += coef * c * w_i;
                }
            }


        }
        *MeshPointStack(stack) = mp;

    }


    // creating an instance of Element_rhs
    // case 3D volume isoline ... levelset ...
    template<class R>
    void  Element_rhs(const FElement3 & Kv,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular & FI ,int np, R3 *Q,int optim)
    {
        //   AFAIRE("Element_rhs on border");
        Stack stack=pvoid2Stack(vstack);
        typedef  FElement3::Element Element;

        MeshPoint mp=*MeshPointStack(stack) ;
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        const Element & K  = Kv.T;
        const Mesh3 & Th= Kv.Vh.Th;
        double epsmes3=K.mesure()*K.mesure()*1e-18;
        long npi;
        long n=Kv.NbDoF(),N=Kv.N;
        double l[3];

        bool classoptm = copt && Op.optiexpK;
        // assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
        if (Kv.number<1 && verbosity/100 && verbosity % 10 == 2)
        cout << "Element_rhs 3d S(levelset): copt = " << copt << " " << classoptm << " opt " << optim << endl;
        int lastop;
        What_d Dop = Op.DiffOp(lastop);

        assert(Op.MaxOp() <last_operatortype);
        // assert(lastop<=3);

        RNMK_ fu(p,n,N,lastop); //  the value for basic fonction
        R3 PP[4];
        for(int i=0; i< np; ++i)
        PP[i]= K(Q[i]);

        for( int iii =0; iii+1 < np; iii+=2)
        { // 0,1,, a and 2,3,0.
            int i0=iii,i1=iii+1,i2=(iii+2)%np;
            R3 NN= R3(PP[i0],PP[i1])^R3(PP[i0],PP[i2]);
            double mes2 = (NN,NN);
            double mes = sqrt(mes2);
            if(mes2*mes <epsmes3) continue; //  too small
            NN /= mes;
            mes *= 0.5;
            // cout << " Element_rhs::mes " << mes << " " << iii << endl;

            for (npi=0;npi<FI.n;npi++) // loop on the integration point
            {
                GQuadraturePoint<R2>  pi( FI[npi]);
                pi.toBary(l);
                R3 Pt( l[0]*Q[i0]+l[1]*Q[i1]+l[2]*Q[i2]); //
                MeshPointStack(stack)->set(Th,K(Pt),Pt,K,-1,NN,-1);
                //
                Kv.BF(Dop,Pt,fu);
                //        MeshPointStack(stack)->set(K(Pt),Pt,Kv,label,NN,ie);
                if (classoptm) (*Op.optiexpK)(stack); // call optim version
                double coef = mes*pi.a;
                for (int  i=0;  i<n;   i++ )
                // if (alledges || onWhatIsEdge[ie][Kv.DFOnWhat(i)]) // bofbof faux si il y a des derives ..
                {
                    RNM_ wi(fu(i,'.','.'));
                    int il=0;
                    for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                    {
                        LOperaD::K ll(*l);
                        pair<int,int> ii(ll.first);
                        double w_i =  wi(ii.first,ii.second);
                        R c =copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                        // FFCS - removing what is probably a small glitch
                        if ( copt && ( optim==1) && Kv.number<1)
                        {
                            R cc  =  GetAny<R>(ll.second.eval(stack));
                            if ( c != cc) {
                                cerr << c << " =! " << cc << endl;
                                cerr << "Sorry error in Optimization (u) add:  int2d(Th,optimize=0)(...)" << endl;
                                ExecError("In Optimized version "); }
                        }


                        //= GetAny<double>(ll.second.eval(stack));

                        B[Kv(i)] += coef * c * w_i;
                    }
                }


            }
        }
        *MeshPointStack(stack) = mp;

    }

    // creating an instance of Element_rhs
    // case 3D surface
    template<class R>
    void  Element_rhs(const FElementS & Kv,int ie,int label,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular1d & FI = QF_GaussLegendre2,bool alledges=false,int optim=1)
    {
        Stack stack=pvoid2Stack(vstack);
        MeshPoint mp=*MeshPointStack(stack) ;
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        const TriangleS & T  = Kv.T;
        // const QuadratureFormular1d & FI = QF_GaussLegendre2;
        long npi;
        long i,n=Kv.NbDoF(),N=Kv.N;

        //  bool show = Kv.Vh.Th(T)==0;
        // char * xxx[] ={" u"," v,"," p"," q"," r"};
        // char * xxxx[] ={" u'"," v',"," p'"," q'"," r'"};
        // char * yyy[] ={" ","_x ","_y "};

        bool classoptm = copt && Op.optiexpK;
        // assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
        if (Kv.number<1 && verbosity/100 && verbosity % 10 == 2)
            cout << "Element_rhs S: copt = " << copt << " " << classoptm << "opt " << optim << endl;

        int lastop=0;
        What_d Dop = Op.DiffOp(lastop);
        assert(Op.MaxOp() <last_operatortype);

        // assert(lastop<=3);

        RNMK_ fu(p,n,N,lastop); //  the value for basic fonction
        // surface normal
         R3 NNt=T.NormalTUnitaire();
         // exterior normal (flux)
         R3 NN=T.N(ie);
         NN /= NN.norme();
        //cout << " NN= " << NN << endl; 
        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            QuadratureFormular1dPoint pi( FI[npi]);
            R3 E=T.Edge(ie);
            double le = sqrt((E,E));
            double coef = le*pi.a;
            double sa=pi.x,sb=1-sa;
            R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
            PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]);
            R2 Pt(PA*sa+PB*sb );
            Kv.BF(Dop,Pt,fu);
             MeshPointStack(stack)->set(T(Pt),Pt,Kv,label,NN,NNt,ie);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version

            for ( i=0;  i<n;   i++ )
                // if (alledges || onWhatIsEdge[ie][Kv.DFOnWhat(i)]) // bofbof faux si il y a des derives ..
            {
                RNM_ wi(fu(i,'.','.'));
                int il=0;
                for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                {
                    LOperaD::K ll(*l);
                    pair<int,int> ii(ll.first);
                    double w_i =  wi(ii.first,ii.second);
                    R c =copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                    // FFCS - removing what is probably a small glitch
                    if ( copt && ( optim==1) && Kv.number<1)
                    {
                        R cc  =  GetAny<R>(ll.second.eval(stack));
                        if ( c != cc) {
                            cerr << c << " =! " << cc << endl;
                            cerr << "Sorry error in Optimization (v) add:  int2d(Th,optimize=0)(...)" << endl;
                            ExecError("In Optimized version "); }
                    }


                    //= GetAny<double>(ll.second.eval(stack));

                    B[Kv(i)] += coef * c * w_i;
                }
            }


        }
        *MeshPointStack(stack) = mp;
    }

    // creating an instance of Element_rhs
    // case 3D curve
    template<class R>
    void  Element_rhs(const FElementL & Kv,int ie,int label,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular1d & FI = QF_GaussLegendre2,
                      bool alledges=false,int optim=1)
    {
       // cout <<  "MMMMMMMMM " <<ie << " " << label << " " <<  optim<< endl;
        
           Stack stack=pvoid2Stack(vstack);
            MeshPoint mp=*MeshPointStack(stack) ;
            R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
            const EdgeL & T  = Kv.T;
            // const QuadratureFormular1d & FI = QF_GaussLegendre2;
            long npi;
            long i,n=Kv.NbDoF(),N=Kv.N;

            //  bool show = Kv.Vh.Th(T)==0;
            // char * xxx[] ={" u"," v,"," p"," q"," r"};
            // char * xxxx[] ={" u'"," v',"," p'"," q'"," r'"};
            // char * yyy[] ={" ","_x ","_y "};

            bool classoptm = copt && Op.optiexpK;
            // assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
            if (Kv.number<1 && verbosity/100 && verbosity % 10 == 2)
            cout << "Element_rhs L: copt = " << copt << " " << classoptm << "opt " << optim << " int0d " << endl;
            int lastop=0;
            What_d Dop = Op.DiffOp(lastop);
             
            // assert(lastop<=3);

            RNMK_ fu(p,n,N,lastop); //  the value for basic fonction

            
            {
                double s = ie;
                QuadratureFormular1dPoint pi(1.,s);
               R1 Pt(pi);
               Kv.BF(Dop,Pt,fu);
                //  calcul de N ...
                R3 Nt=T.TangenteUnitaire(),NN=Nt;
                if(ie==0) NN=-NN;
                MeshPointStack(stack)->set(T(Pt),Pt,Kv,label,NN,Nt,ie);//;set(T(Pt),Pt,Kv);//  pas bon ...
                if (classoptm) (*Op.optiexpK)(stack); // call optim version
               
                for ( i=0;  i<n;   i++ )
                 {
                    RNM_ wi(fu(i,'.','.'));
                    int il=0;
                    for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                    {
                        LOperaD::K ll(*l);
                        pair<int,int> ii(ll.first);
                        double w_i =  wi(ii.first,ii.second);
                        R c =copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
 
                        B[Kv(i)] +=  c * w_i;
                    }
                }


            }
            *MeshPointStack(stack) = mp;

        
       }


    // end 3d


    // creating an instance of Element_rhs
    // case 2d
    template<class R>
    void  Element_rhs(const FElement & Kv,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular1d & FI ,const R2 & PPA,const R2 &PPB,int optim)
    {
        Stack stack=pvoid2Stack(vstack);
        MeshPoint mp=*MeshPointStack(stack) ;
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        const Triangle & T  = Kv.T;
        R2 PA=T(PPA),PB=T(PPB);
        // const QuadratureFormular1d & FI = QF_GaussLegendre2;
        long npi;
        long i,n=Kv.NbDoF(),N=Kv.N;

        //  bool show = Kv.Vh.Th(T)==0;
        // char * xxx[] ={" u"," v,"," p"," q"," r"};
        // char * xxxx[] ={" u'"," v',"," p'"," q'"," r'"};
        // char * yyy[] ={" ","_x ","_y "};

        bool classoptm = copt && Op.optiexpK;
        // assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
        if (Kv.number<1 && verbosity/100 && verbosity % 10 == 2)
        cout << "Element_rhs(levelset) S: copt = " << copt << " " << classoptm <<" opt " << optim << endl;
        KN<bool> Dop(last_operatortype);
        Op.DiffOp(Dop);
        int lastop=1+Dop.last([](bool x){return x;});
        assert(Op.MaxOp() <last_operatortype);
        // assert(lastop<=3);

        RNMK_ fu(p,n,N,lastop); //  the value for basic fonction

        R2 E(PA,PB);
        double le = sqrt((E,E));

        //cout << " Element_rhs 2d " << PA << " " << PB << " " << le << " " << Kv.number <<  endl;
        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            QuadratureFormular1dPoint pi( FI[npi]);
            double coef = le*pi.a;
            double sa=pi.x,sb=1-sa;
            R2 Pt(PPA*sa+PPB*sb ); //
            Kv.BF(Dop,Pt,fu);
            MeshPointStack(stack)->set(T(Pt),Pt,Kv,0,R2(E.y,-E.x)/le,0);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version

            for ( i=0;  i<n;   i++ )
            // if (alledges || onWhatIsEdge[ie][Kv.DFOnWhat(i)]) // bofbof faux si il y a des derives ..
            {
                RNM_ wi(fu(i,'.','.'));
                int il=0;
                for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                {
                    LOperaD::K ll(*l);
                    pair<int,int> ii(ll.first);
                    double w_i =  wi(ii.first,ii.second);
                    R c =copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                    // FFCS - removing what is probably a small glitch
                    if ( copt && ( optim==1) && Kv.number<1)
                    {
                        R cc  =  GetAny<R>(ll.second.eval(stack));
                        if ( c != cc) {
                            cerr << c << " =! " << cc << endl;
                            cerr << "Sorry error in Optimization (w)  add:  int2d(Th,optimize=0)(...)" << endl;
                            ExecError("In Optimized version "); }
                    }


                    //= GetAny<double>(ll.second.eval(stack));
                  //  cout << "         " << coef<< " " << c << " " << w_i << " " << Kv(i) << " | " << Pt <<  endl;
                    B[Kv(i)] += coef * c * w_i;
                }
            }


        }
        *MeshPointStack(stack) = mp;

    }


    // case 3D surface
    template<class R>
    void  Element_rhs(const FElementS & Kv,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular1d & FI ,const R3 & PPA,const R3 &PPB,int optim)
    {
        ffassert(0);
    }

    // case 3D curve
    template<class R>
    void  Element_rhs(const FElementL & Kv,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular1d & FI ,const R3 & PPA,const R3 &PPB,int optim)
    {
        ffassert(0);
    }



    // creating an instance of Element_rhsVF
    // case 2d
    template<class R>
    void  Element_rhsVF(const FElement & Kv,const FElement & KKv,int ie,int iie,int label,const LOperaD &Op,double * p,int *ip,void  * bstack,KN_<R> & B,
                        const QuadratureFormular1d & FI = QF_GaussLegendre2,int optim=1)
    // sier of ip
    //  version correct the  29 april 2015 by. FH
    //  missing before in case of jump, mean , .. in test functions
    //  Thank to Lucas Franceschini <lucas.franceschini@ensta-paristech.fr>
    {
        pair_stack_double * bs=static_cast<pair_stack_double *>(bstack);
        Stack stack= bs->first;
        double binside = *bs->second; // truc FH pour fluide de grad2 (decentrage bizard)
        bool onborder= &Kv.T == &KKv.T;
        const FElement *pKKv= !onborder ?  & KKv : 0;
        MeshPoint mp=*MeshPointStack(stack) ;
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);
        const Triangle & T  = Kv.T;
        // const QuadratureFormular1d & FI = QF_GaussLegendre2;
        long npi;
        long i,nv=Kv.NbDoF(),N=Kv.N;
        long nnv=KKv.NbDoF();
        assert(nv==nnv);
        //  bool show = Kv.Vh.Th(T)==0;
        // char * xxx[] ={" u"," v,"," p"," q"," r"};
        // char * xxxx[] ={" u'"," v',"," p'"," q'"," r'"};
        // char * yyy[] ={" ","_x ","_y "};

        bool classoptm = copt && Op.optiexpK;
        // assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
        if (Kv.number<1 && verbosity/100 && verbosity % 10 == 2)
        cout << "Element_rhs S: copt = " << copt << " " << classoptm << " opt " << optim << endl;
        KN<bool> Dop(last_operatortype);
        Op.DiffOp(Dop);
        int lastop=1+Dop.last([](bool x){return x;});
        //assert(Op.MaxOp() <last_operatortype);
        // assert(lastop<=3);
        int lffv = nv*N*last_operatortype;
        int lp =nv*2;
        KN_<int> pp(ip,lp),pk(ip+lp,lp),pkk(ip+2*lp,lp);
        int n = BuildMEK_KK(lp,pp,pk,pkk,&Kv,pKKv);
        RNMK_ fu(p,nv,N,lastop); //  the value for basic fonction
        RNMK_ ffu( (double*) p  + lffv  ,nv,N,lastop); //  the value for basic fonction

        R2 E=T.Edge(ie);
        double le = sqrt((E,E));
        R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
        PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]),
        PC(TriangleHat[OppositeVertex[ie]]);
        // warning the to edge are in opposite sens
        R2 PP_A(TriangleHat[VerticesOfTriangularEdge[iie][1]]),
        PP_B(TriangleHat[VerticesOfTriangularEdge[iie][0]]),
        PP_C(TriangleHat[OppositeVertex[ie]]);
        R2 Normal(E.perp()/-le);
        double cmean = onborder ? 1. : 0.5;
        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            QuadratureFormular1dPoint pi( FI[npi]);
            R2 E=T.Edge(ie);
            double le = sqrt((E,E));
            double coef = le*pi.a;
            double sa=pi.x,sb=1-sa;
            R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
            PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]);
            R2 Pt(PA*sa+PB*sb ); //
            R2 PP_t(PP_A*sa+PP_B*sb ); //
            if (binside) {
                Pt   = (1-binside)*Pt + binside*PC;
                PP_t  = (1-binside)*PP_t + binside*PP_C; }
            Kv.BF(Dop,Pt,fu);
            if(onborder)
            ffu=0;
            else
            KKv.BF(Dop,PP_t,ffu);

            MeshPointStack(stack)->set(T(Pt),Pt,Kv,label,R2(E.y,-E.x)/le,ie);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version

            for ( i=0;  i<n;   i++ )
            // if (alledges || onWhatIsEdge[ie][Kv.DFOnWhat(i)]) // bofbof faux si il y a des derives ..
            {
                int ik= pk[i];
                int ikk=pkk[i];

                RNM_ wi(fu(Max(ik,0),'.','.'));
                RNM_ wwi(ffu(Max(ikk,0),'.','.'));
                int il=0;
                int dofik=ik>=0? Kv(ik):-1;
                int dofikk=ikk>=0? KKv(ikk):-1;

                for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                {


                    LOperaD::K ll(*l);
                    pair<int,int> ii(ll.first);
                    int iis = ii.second;
                    int iicase  = iis / last_operatortype;
                    iis %= last_operatortype;
                    double w_i=0,ww_i=0;
                    if(ik>=0) w_i =   wi(ii.first,iis );
                    if( iicase>0 )
                    {
                        if( ikk>=0) ww_i =  wwi(ii.first,iis );
                        if       (iicase==Code_Jump)      w_i = -w_i; ///(w_i = ww_i-w_i); // jump
                        else  if (iicase==Code_Mean)      ww_i=w_i = cmean*  (w_i + ww_i ); // average
                        else  if (iicase==Code_OtherSide) std::swap(w_i,ww_i);  // valeur de autre cote
                        else ffassert(0);
                    }
                    R c =copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                    // FFCS - removing what is probably a small glitch
                    if ( copt && ( optim==1) && Kv.number<1)
                    {
                        R cc  =  GetAny<R>(ll.second.eval(stack));
                        if ( c != cc) {
                            cerr << c << " =! " << cc << endl;
                            cerr << "Sorry error in Optimization (x) add:  int2d(Th,optimize=0)(...)" << endl;
                            ExecError("In Optimized version "); }
                    }


                    //= GetAny<double>(ll.second.eval(stack));

                    if(dofik>=0) B[dofik] += coef * c * w_i;
                    if(dofikk>=0) B[dofikk] += coef * c * ww_i;

                }
            }


        }
        *MeshPointStack(stack) = mp;

    }

// creating an instance of Element_rhs
// case 3D volume


//  VFVF
template<class R>
void  Element_rhsVF(const FElement3 & Kv,const FElement3 & KKv,int ie,int iie,int label,const LOperaD &Op,double * p,int *ip,void  * bstack,KN_<R> & B,
                    const QuadratureFormular & FIb,
                    int optim=1)
{
   
      int intmortar=0;
  //     AFAIRE("Element_rhsVF 3d  "); A FAit 23 sep 2022 FH
    typedef typename FElement3::Element Element;
    ffassert(B!=0);
    pair_stack_double * bs=static_cast<pair_stack_double *>(bstack);
    Stack stack= bs->first;
    double binside = *bs->second; // truc FH pour fluide de grad2 (decentrage bizard)
    ffassert(ie>=0 && ie < 4); //  int on Face
    MeshPoint mp= *MeshPointStack(stack);
    R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);

    const Element & T  = Kv.T;
    const Element & TT  = KKv.T;
    bool sameT = &T == & TT;
    int nTonEdge =  &T == &TT ? 1 : 2;
    double cmean = 1./nTonEdge;
    bool onborder= &T == &TT;
    const FElement3 *pKKv= !onborder ?  & KKv : 0;

    long npi;
    long i;
    long N= Kv.N;
    long nv=Kv.NbDoF();
    long nnv=Kv.NbDoF();
    assert(nv == nnv) ;
    int lp =nv*2;
    KN_<int> pp(ip,lp),pk(ip+lp,lp),pkk(ip+2*lp,lp);
    int n = BuildMEK_KK(lp,pp,pk,pkk,&Kv,pKKv);

    bool classoptm = copt && Op.optiexpK;// a voir !!!!
    if (Kv.number<1 && ( verbosity > 1 ) )
    cout << "Element_rhsVF 3d P: copt = " << copt << " " << classoptm << " binside (For FH) =" << binside << " opt: " << optim << endl;

     bool oldopt=1;
    int  iloop=0;
    KN<bool> unvarexp(classoptm ? Op.optiexpK->sizevar() : 1);
    if (Kv.number<1 && verbosity/100 && verbosity % 10 == 2)
    cout << "Element_rhsVF 3d P: copt = " << copt << " " << classoptm << " opt: " << optim << endl;
     //
    int lastop=0;
    lastop = 0;
    What_d Dop = Op.DiffOp(lastop);

    long lffv = sameT ? 0  : nv*N*lastop;

    RNMK_ fv(p,nv,N,lastop); //  the value for basic fonction in K
    RNMK_ ffv(p + lffv ,nnv,N,lastop); //  the value for basic fonction in KK
     
    


    R3 NN= T.N(ie);
    double mes=NN.norme();
    NN/=mes;
//  PB for find common dof
    // compute the permutaion face ie to iie ...
    // a little tricky
    int fpe= T.facePermutation(ie);
    int fpee= TT.facePermutation(iie);
    int pr[3],ppr[3];
    SetNumPerm<3>(fpe,pr);
    SetNumPerm1<3>(fpee,ppr);

    for (npi=0;npi<FIb.n;npi++) // loop on the integration point
    {
  
        GQuadraturePoint<R2> pi( FIb[npi]);
        double lpi[3];
        pi.toBary(lpi);
        double llpi[]={lpi[pr[ppr[0]]], lpi[pr[ppr[1]]], lpi[pr[ppr[2]]]};
        double coef = 0.5*mes*pi.a; // correction 0.5 050109 FH
        R3 Pt(T.PBord(ie,pi));
        R2 pii(llpi+1);
        R3 PP_t(TT.PBord(iie,pii));
        {
            static int err=0;
            R3 P(T(Pt)),PP(TT(PP_t)),D(P,PP);
           if(  D.norme2() > 1e-5)
           {
               const Mesh3 & Th=Kv.Vh.Th;
               cout << " pr "<< pr[0] << pr[1]<< pr[2] << " prr " << ppr[0] << ppr[1]<< ppr[2]<<" "
               << ppr[pr[0]] <<  ppr[pr[1]] << ppr[pr[2]]  << "   " <<  pr[ppr[0]] <<  pr[ppr[1]] << pr[ppr[2]]<< endl;
               cout << ie << " " << iie << endl;
               cout << " T = " << Th(T[0]) << " " <<Th(T[1]) << " " << Th(T[2]) << " " << Th(T[3]) << endl;
               cout << " TT = " << Th(TT[0]) << " " <<Th(TT[1]) << " " << Th(TT[2]) << " " << Th(TT[3]) << endl;
               cout << fpe << " " << fpee << endl;
               cout << P << " " << PP << " diff=D.norme2() " << D.norme2()  << endl;
               err++;
               ffassert(err<10);
           }
            
        }
        Kv.BF(Dop,Pt,fv);
        if(!sameT) KKv.BF(Dop,PP_t,ffv);
        if(verbosity==99999 &&( Kv.number == 3 || KKv.number == 2)) {
            cout << "intallfaces " << Kv.number << " " << fv << endl;
            cout << "intallfaces " << KKv.number << " " << ffv << endl;
        }
        MeshPointStack(stack)->set(T(Pt),Pt,Kv,label,NN,ie);
        if (classoptm) (*Op.optiexpK)(stack); // call optim version


        



        for ( i=0;  i<n;   i++ )
        {
            int ik= pk[i];
            int ikk=pkk[i];
            int dofik=ik>=0? Kv(ik):-1;
            int dofikk=ikk>=0? KKv(ikk):-1;

            RNM_ wi(fv(Max(ik,0),'.','.'));
            RNM_ wwi(ffv(Max(ikk,0),'.','.'));


                int il=0;
                for (auto l=Op.v.begin();l!=Op.v.end();l++,il++)
                {
                    auto ll(*l);
                    pair<int,int> ii(ll.first);
                    int iis = ii.second;

                    int iicase  = iis / last_operatortype;

                    iis %= last_operatortype;
              /*
                    double w_i=0,w_j=0,ww_i=0;

                    if(ik>=0) w_i =   wi(ii.first,iis );
                

                    if( iicase>0 && ikk>=0) ww_i =  wwi(ii.first,iis );
                    


                    if       (iicase==Code_Jump) ww_i=-(w_i = ww_i-w_i); // jump
                    else  if (iicase==Code_Mean) {

                        ww_i=w_i = cmean*  (w_i + ww_i );} // average
                    else  if (iicase==Code_OtherSide) swap(w_i,ww_i);  // valeur de autre cote

 
                    // R ccc = GetAny<R>(ll.second.eval(stack));

                    R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                    if ( copt && ( optim==1) && Kv.number <1)
                    {
                        R cc  =  GetAny<R>(ll.second.eval(stack));
                        CheckErrorOptimisation(cc,ccc,"Sorry error in Optimization Element_OpVF3d  (face) add:   intallface(Th,optimize=0)(...)");
                     }
                    if(verbosity==99999)
                        cout << ie << " / " << i <<  "/  "<<"  pi " << npi << " : " << coef << " c= " << ccc << " " << w_i  << " = " << coef * ccc * w_i
                        << " k/kk  " << Kv.number << " "<< KKv.number << " :: " << ik << " " << ikk  << ",  " << " ii" << ii.second << " w= "   // << " // " << bbnElementonB(stack)
                        << w_i << " " << ww_i << endl;
                    if(dofik>=0) B[dofik] += coef * ccc * w_i;
                    if(dofikk>=0) B[dofikk] += coef * ccc * ww_i;
                    */
                    
                    double w_i=0,ww_i=0;
                    if(ik>=0) w_i =   wi(ii.first,iis );
                    if( iicase>0 )
                    {
                        if( ikk>=0) ww_i =  wwi(ii.first,iis );
                        if       (iicase==  Code_Jump)      w_i = -w_i; ///(w_i = ww_i-w_i); // jump
                        else  if (iicase==Code_Mean)      ww_i=w_i = cmean*  (w_i + ww_i ); // average
                        else  if (iicase==Code_OtherSide) std::swap(w_i,ww_i);  // valeur de autre cote
                        else ffassert(0);
                    }
                    R c =copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                    // FFCS - removing what is probably a small glitch
                    if ( copt && ( optim==1) && Kv.number<1)
                    {
                        R cc  =  GetAny<R>(ll.second.eval(stack));
                        if ( c != cc) {
                            cerr << c << " =! " << cc << endl;
                            cerr << "Sorry error in Optimization (x) add:  int2d(Th,optimize=0)(...)" << endl;
                            ExecError("In Optimized version "); }
                    }

                    if(verbosity>999)
                    cout<< iicase << " : " << dofik << " " << dofikk << " :: " << coef << " "<< c << " " << w_i << " " << ww_i << " K  = "<< Kv.number<< " " << KKv.number << " ?? ";
 
                    //= GetAny<double>(ll.second.eval(stack));

                    if(dofik>=0) B[dofik] += coef * c * w_i;
                    if(dofikk>=0) B[dofikk] += coef * c * ww_i;
                    if(verbosity>999)
                    {
                        if(dofik>=0)   cout << " + "<< dofik << " " << coef * c * w_i << " b " << B[dofik]   << " , cm  " << cmean;
                        if(dofikk>=0)  cout << " + "<< dofikk <<" " <<  coef * c * ww_i << " b " << B[dofikk] << " , cm  " << cmean;
                        cout << endl;
                    }

               }
            }
        }
    


    *MeshPointStack(stack) = mp;

    }
// creating an instance of Element_rhs
// case 3D volume
template<class R>
void  Element_rhs(const  Mesh3 & ThI,const Mesh3::Element & KI, const FESpace3 & Vh,
                  int ie,int label,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                  const QuadratureFormular & FI,bool alledges=false,int optim=1)
{
    Stack stack=pvoid2Stack(vstack);
    int intmortar=0;
    //  AFAIRE("Element_rhs 3d on surface  2 diff mesh ");
    static int count =0;
    if(count++<1)
    {
        cout << " Element_rhs 3d on surface  2 diff mesh int test (FH)" << endl;
        cout << " -----------------------------------------------------" << endl;
    }
    // integration 1d on 2 diff mesh


    MeshPoint mp=*MeshPointStack(stack) ;
    R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);


    bool classoptm = copt && Op.optiexpK;
    //assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
    if (ThI(KI)<1 && verbosity/100 && verbosity % 10 == 2)
    cout << "Element_rhs S: copt = " << copt << " " << classoptm << " opt "<< optim <<endl;
    assert(Op.MaxOp() <last_operatortype);
    //
    int lastop=0;
    lastop = 0;
    What_d Dop = Op.DiffOp(lastop);
    // assert(lastop<=3);
    const Tet & T  = KI;
    long npi;

    const Tet * Kp=0;

    for (npi=0;npi<FI.n;npi++) // loop on the integration point
    {

        GQuadraturePoint<R2> pi(FI[npi]);
        R3 NN= T.N(ie);
        double mes=NN.norme();
        NN/=mes;
        double coef = 0.5*mes*pi.a; //
        R3 Pt(T.PBord(ie,pi)),PI(T(Pt));



        MeshPointStack(stack)->set(ThI,PI,Pt,KI,label,NN,ie);
        if (classoptm) (*Op.optiexpK)(stack); // call optim version
        bool outside;
        R3 PIt;
        const Tet & K  = *Vh.Th.Find(PI,PIt,outside,Kp);
        if ( ! outside || intmortar) //  FH march 2009 ???
        {
            const  FElement3  Kv= Vh[K];
            long i,n=Kv.NbDoF(),N=Kv.N;
            RNMK_ fu(p,n,N,lastop); //  the value for basic fonction
            Kv.BF(Dop,PIt,fu);

            for ( i=0;  i<n;   i++ )
            //   if (alledges || onWhatIsFace[ie][Kv.DFOnWhat(i)]) // bofbof faux si il y a des derives ..
            {
                RNM_ wi(fu(i,'.','.'));
                int il=0;
                for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                {
                    LOperaD::K ll(*l);
                    pair<int,int> ii(ll.first);
                    double w_i =  wi(ii.first,ii.second);
                    R c =copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                    // FFCS - removing what is probably a small glitch
                    if ( copt && ( optim==1) && Kv.number<1)
                    {
                        R cc  =  GetAny<R>(ll.second.eval(stack));
                        if ( c != cc) {
                            cerr << c << " =! " << cc << endl;
                            cerr << "Sorry error in Optimization (y) add:  int1d(Th,optimize=0)(...)" << endl;
                            ExecError("In Optimized version "); }
                    }


                    //= GetAny<double>(ll.second.eval(stack));

                    B[Kv(i)] += coef * c * w_i;
                }
            }

        }
    }
    *MeshPointStack(stack) = mp;

}
    // creating an instance of Element_rhs
    // case 3D surface
    template<class R>
    void  Element_rhs(const  MeshS & ThI,const TriangleS & KI, const FESpaceS & Vh,
                      int ie,int label,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular1d & FI = QF_GaussLegendre2,bool alledges=false,bool intmortar=false,
                      R2 *Q=0,int optim=1)
    {
        ffassert(0);
    }

    // creating an instance of Element_rhs
    // case 3D curve
    template<class R>
    void  Element_rhs(const  MeshL & ThI,const EdgeL & KI, const FESpaceL & Vh,
                      int ie,int label,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular1d & FI = QF_GaussLegendre2,bool alledges=false,bool intmortar=false,
                      R2 *Q=0,int optim=1)
    {
        ffassert(0);
    }

    // creating an instance of Element_rhs
    // case 2d
    template<class R>
    void  Element_rhs(const  Mesh & ThI,const Triangle & KI, const FESpace & Vh,
                      int ie,int label,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular1d & FI = QF_GaussLegendre2,bool alledges=false,bool intmortar=false,
                      R2 *Q=0,int optim=1)
    {
        // integration 1d on 2 diff mesh

        Stack stack=pvoid2Stack(vstack);
        MeshPoint mp=*MeshPointStack(stack) ;
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);


        bool classoptm = copt && Op.optiexpK;
        //assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
        if (ThI.number(KI)<1 && verbosity/100 && verbosity % 10 == 2)
        cout << "Element_rhs S: copt = " << copt << " " << classoptm << " opt " << optim << endl;
        KN<bool> Dop(last_operatortype);
        Op.DiffOp(Dop);
        int lastop=1+Dop.last([](bool x){return x;});
        assert(Op.MaxOp() <last_operatortype);
        // assert(lastop<=3);
        const Triangle & T  = KI;
        long npi;

        const Triangle * Kp=0;
        R2 PA,PB,E;
        if( Q==0)
        {
            PA=TriangleHat[VerticesOfTriangularEdge[ie][0]];
            PB=TriangleHat[VerticesOfTriangularEdge[ie][1]];
            E=T.Edge(ie);
        }
        else
        {
            PA=Q[0];
            PB=Q[1];
            E=T(PB)-T(PA);
        }
        double le = sqrt((E,E));

        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            QuadratureFormular1dPoint pi( FI[npi]);


            double coef = le*pi.a;
            double sa=pi.x,sb=1-sa;
            R2 Pt(PA*sa+PB*sb ); //
            R2 PI(KI(Pt));
            //   Kv.BF(Dop,Pt,fu);
            MeshPointStack(stack)->set(ThI,PI,Pt,KI,label,R2(E.y,-E.x)/le,ie);
            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            bool outside;
            R2 PIt;
            const Triangle & K  = *Vh.Th.Find(PI,PIt,outside,Kp);
            if ( ! outside || intmortar) //  FH march 2009 ???
            {
                const  FElement  Kv= Vh[K];
                long i,n=Kv.NbDoF(),N=Kv.N;
                RNMK_ fu(p,n,N,lastop); //  the value for basic fonction
                Kv.BF(Dop,PIt,fu);

                for ( i=0;  i<n;   i++ )
                // if (alledges || onWhatIsEdge[ie][Kv.DFOnWhat(i)]) // bofbof faux si il y a des derives ..
                {
                    RNM_ wi(fu(i,'.','.'));
                    int il=0;
                    for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                    {
                        LOperaD::K ll(*l);
                        pair<int,int> ii(ll.first);
                        double w_i =  wi(ii.first,ii.second);
                        R c =copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                        // FFCS - removing what is probably a small glitch
                        if ( copt && ( optim==1) && Kv.number<1)
                        {
                            R cc  =  GetAny<R>(ll.second.eval(stack));
                            if ( c != cc) {
                                cerr << c << " =! " << cc << endl;
                                cerr << "Sorry error in Optimization (z) add:  int1d(Th,optimize=0)(...)" << endl;
                                ExecError("In Optimized version "); }
                        }


                        //= GetAny<double>(ll.second.eval(stack));

                        B[Kv(i)] += coef * c * w_i;
                    }
                }

            }
        }
        *MeshPointStack(stack) = mp;

    }
    // creating an instance of Element_rhs
    // case 2d
    template<class R>
    void  Element_rhs(Expression const * const mapt,const  Mesh & ThI,const Triangle & KI, const FESpace & Vh,
                      int ie,int label,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular1d & FI = QF_GaussLegendre2,bool alledges=false,bool intmortar=false,
                      R2 *Q=0,int optim=1)
    {
        // integration 1d on 2 diff mesh
        //  ffassert(0);
        Stack stack=pvoid2Stack(vstack);
        MeshPoint mp=*MeshPointStack(stack) ;
        R ** copt = Stack_Ptr<R*>(stack,ElemMatPtrOffset);


        bool classoptm = copt && Op.optiexpK;
        //assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
        if (ThI.number(KI)<1 && verbosity/100 && verbosity % 10 == 2)
        cout << "Element_rhs S: copt = " << copt << " " << classoptm << " opt " << optim << endl;
        KN<bool> Dop(last_operatortype);
        Op.DiffOp(Dop);
        int lastop=1+Dop.last([](bool x){return x;});
        assert(Op.MaxOp() <last_operatortype);
        // assert(lastop<=3);
        const Triangle & T  = KI;
        long npi;

        const Triangle * Kp=0;
        R2 PA,PB,E;
        if( Q==0)
        {
            PA=TriangleHat[VerticesOfTriangularEdge[ie][0]];
            PB=TriangleHat[VerticesOfTriangularEdge[ie][1]];
            E=T.Edge(ie);
        }
        else
        {
            PA=Q[0];
            PB=Q[1];
            E=T(PB)-T(PA);
        }
        double le = sqrt((E,E));

        for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
            QuadratureFormular1dPoint pi( FI[npi]);


            double coef = le*pi.a;
            double sa=pi.x,sb=1-sa;
            R2 Pt(PA*sa+PB*sb ); //
            R2 PIo(KI(Pt)),PI(PIo);
            //   Kv.BF(Dop,Pt,fu);
            MeshPointStack(stack)->set(ThI,PI,Pt,KI,label,R2(E.y,-E.x)/le,ie);
            if(mapt)
            { // move poit
                PI= R2( GetAny<double>((*mapt[0])(vstack)), GetAny<double>((*mapt[1])(vstack)));
                if(verbosity>9999) cout << "  Element_rhs(2) mapt =" << PIo << " -> " <<PI << endl;
            }

            if (classoptm) (*Op.optiexpK)(stack); // call optim version
            bool outside;
            R2 PIt;
            const Triangle & K  = *Vh.Th.Find(PI,PIt,outside,Kp);
            if ( ! outside || intmortar) //  FH march 2009 ???
            {
                const  FElement  Kv= Vh[K];
                long i,n=Kv.NbDoF(),N=Kv.N;
                RNMK_ fu(p,n,N,lastop); //  the value for basic fonction
                Kv.BF(Dop,PIt,fu);

                for ( i=0;  i<n;   i++ )
                // if (alledges || onWhatIsEdge[ie][Kv.DFOnWhat(i)]) // bofbof faux si il y a des derives ..
                {
                    RNM_ wi(fu(i,'.','.'));
                    int il=0;
                    for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                    {
                        LOperaD::K ll(*l);
                        pair<int,int> ii(ll.first);
                        double w_i =  wi(ii.first,ii.second);
                        R c =copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                        // FFCS - removing what is probably a small glitch
                        if ( copt && ( optim==1) && Kv.number<1)
                        {
                            R cc  =  GetAny<R>(ll.second.eval(stack));
                            if ( c != cc) {
                                cerr << c << " =! " << cc << endl;
                                cerr << "Sorry error in Optimization (z) add:  int1d(Th,optimize=0)(...)" << endl;
                                ExecError("In Optimized version "); }
                        }


                        //= GetAny<double>(ll.second.eval(stack));

                        B[Kv(i)] += coef * c * w_i;
                    }
                }

            }
        }
        *MeshPointStack(stack) = mp;

    }


    // creating an instance of Element_rhs
    // case 3D surface
    template<class R>
    void  Element_rhs(Expression const * const mapt,const  MeshS & ThI,const TriangleS & KI, const FESpaceS & Vh,
                      int ie,int label,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular1d & FI = QF_GaussLegendre2,bool alledges=false,bool intmortar=false,
                      R2 *Q=0,int optim=1)
    {
     ffassert(0);
    }
    // creating an instance of Element_rhs
    // case 3D curve
    template<class R>
    void  Element_rhs(Expression const * const mapt,const  MeshL & ThI,const EdgeL & KI, const FESpaceL & Vh,
                      int ie,int label,const LOperaD &Op,double * p,void * vstack,KN_<R> & B,
                      const QuadratureFormular1d & FI = QF_GaussLegendre2,bool alledges=false,bool intmortar=false,
                      R1 *Q=0,int optim=1)
    {
        ffassert(0);
    }

// generic template for AssembleVarForm
template<class R,typename MC,class MMesh,class FESpace1,class FESpace2>
bool AssembleVarForm(Stack stack,const MMesh & Th,const FESpace1 & Uh,const FESpace2 & Vh,bool sym,
                     MC  * A,KN_<R> * B,const list<C_F0> &largs)
{ // return true if BC
    
    typedef MMesh * pmesh; // integration mesh type
    bool ret=false;
    typedef DotStar_KN_<R> DotStar;
    typedef DotSlash_KN_<R> DotSlash;
    list<C_F0>::const_iterator ii,ib=largs.begin(),
    ie=largs.end();
    using namespace FreeFempp;
    TypeVarForm<R> *tvf=TypeVarForm<R>::Global;
    assert( tvf);
    for (ii=ib;ii != ie;ii++)
    {
        Expression e=ii->LeftValue();
        aType r = ii->left();
        //    if(verbosity > 99)   cout <<  << "AssembleVarForm " <<  << " " <<  (*A)(0,3) << endl;
        if (r==  tvf->tFB)
        {
            if (A)
            {
                const  FormBilinear * bf =dynamic_cast<const  FormBilinear *>(e);
                if((bf->di->d != MMesh::Rd::d) || (bf->di->dHat != MMesh::RdHat::d) )
                {
                    cout << " Errer: int "<< bf->di->dHat << "d case ( Bilinear Form ) on Mesh"<< MMesh::RdHat::d <<" Bizarre  !!!!! "<< endl;
                    cout << "    dim coord (template) "<<MMesh::Rd::d << " mesh : " << bf->di->d << endl;

                  ffassert(0);        
          
                }
                else {
                    pmesh  Thbf= GetAny<pmesh>((*bf->di->Th)(stack));
                    if(Thbf) AssembleBilinearForm<R>( stack,*Thbf,Uh,Vh,sym,*A,bf);
                }
            }
        }
        else if (r==tvf->tMat)
        {
            if (A)
                InternalError(" Add sparse matrice; to do, sorry");
        }
        else if (r==tvf->tFL)
        {
            if (B) {
                const  FormLinear * bf =dynamic_cast<const  FormLinear *>(e);
                if(bf->di->d != MMesh::Rd::d )
                {
                    
                    if( bf->di->dHat==2)
                    {
                        cout << " int on MeshS toDo  ( Linear Form )" << endl;
                        ffassert(0);
                        
                    }
                    else  if( bf->di->dHat==1)
                    {
                        cout << " int on MeshL toDo  ( Linear Form )" << endl;
                        ffassert(0);
                        
                    }
                    else if(bf->di->d != MMesh::Rd::d ){
                        cout << " int 2d case  on Mesh ( Linear Form )"<< MMesh::Rd::d <<" debile !!!!! "<< endl;
                        ffassert(0);}
                }
                
                else
                {
                    pmesh  Thbf= GetAny<pmesh>((*bf->di->Th)(stack));
                    if(Thbf) AssembleLinearForm<R>( stack,*Thbf, Vh, B,bf);
                }}
        }
        else if (r==tvf->tTab)
        {
            if ( B)
                *B += *GetAny<KN<R> *>( (*e)(stack) );
        }
        else if (r==tvf->tDotStar)
        {
            if ( B)
            {
                DotStar ab=GetAny<DotStar>( (*e)(stack) );
                *B += ab;
            }
        }
        else if (r==tvf->tMatX)
        {
            if ( B)
            {
                *B += GetAny<typename RNM_VirtualMatrix<R>::plusAx >( (*e)(stack) )  ;
            }
        }
        else if (r==tvf->tMatTX)
        {
            if ( B)
            {
                *B += GetAny<typename RNM_VirtualMatrix<R>::plusAtx >( (*e)(stack) )  ;
            }
        }
        else if (r== tvf->tBC)
            ret=true;
        else
        {
            cerr << "AssembleVarForm  invalid type : " << * r <<  endl;
            throw(ErrorExec("AssembleVarForm invalid type in varf",1));
        }
    }
    return ret;
}

    template<class R,class MMesh,class FESpace1,class FESpace2>
    void AssembleBC(Stack stack,const MMesh & Th,const FESpace1 & Uh,const FESpace2 & Vh,bool sym,
                    MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X, const list<C_F0> &largs , double tgv  )
    {
        list<C_F0>::const_iterator ii,ib=largs.begin(),
        ie=largs.end();
        aType tBC( atype<const  BC_set  *>()) ;
        for (ii=ib;ii != ie;ii++)
        {
            Expression e=ii->LeftValue();
            aType r = ii->left();
            if (r==tBC)
            AssembleBC(stack,Th,Uh,Vh,sym,A,B,X, dynamic_cast<const  BC_set *>(e),tgv);
        }

    }

   //////////////////////////////////
   // AssembleBC
   //////////////////////////////////

  // creating an instance of AssembleBC
    // case 2d
    template<class R>
    void AssembleBC(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                    MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X, const  BC_set * bc, double tgv  )

    {
        MeshPoint *mps= MeshPointStack(stack),mp=*mps;
        StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
        bool sptrclean=true;
        //     sptr->clean(); // modif FH mars 2006  clean Ptr

        int ktbc=0, nbon =0;
        bool Aii = A && A->n == A->m;

        int Nbcomp=Vh.N;
        Check(bc,Nbcomp);
        ffassert(Vh.N == Uh.N);
        TabFuncArg tabexp(stack,Vh.N);
        KN<double> buf((long int)Vh.MaximalNbOfDF() * (long int)last_operatortype * (long int)Vh.N);
        int ndofBC = Aii ?  A->n : 1;
        KN<char> onBC(ndofBC);
        onBC= '\0';

        KN<R> gg(buf);
        if ( B && B->N() != Vh.NbOfDF) ExecError("AssembleBC size rhs and nb of DF of Vh");
        if(verbosity>99) cout << " Problem : BC_set "<< typeid(R).name() << " " ;
        nbon =bc->on.size();
        set<long> on;
        Expandsetoflab(stack,*bc, on);
        /*
         for (int i=0;i<nbon;i++)
         {
         long  lab  = GetAny<long>( (*bc->on[i])(stack));
         if(verbosity>99) cout << lab << " " ;
         on.insert(lab);
         }
         if(verbosity>99)
         cout << endl;
         */
        int kk=bc->bc.size();

        const int dim=Vh.N;
        FElement::aIPJ ipj(Vh[0].Pi_h_ipj());
        FElement::aR2  PtHat(Vh[0].Pi_h_R2());

        KN<int> PtonB(PtHat.N());

        KN<double>   Aipj(ipj.N());
        KNM<R>  Vp(dim,PtHat.N());
        double tgv1=tgv <0? 1: tgv; // change 21 dec 2010 FH (Hack of ILU)
        for (int ib=0;ib<Th.neb;ib++)
        {
            int ie;
            int it = Th.BoundaryElement(ib,ie);
            int r =Th.bedges[ib].lab;
            if (on.find(r) != on.end() )
            {
                const FElement K(Uh[it]);
                R2 E=K.T.Edge(ie);
                double le = sqrt((E,E));

                ktbc++;
                if(verbosity>99)   cout << "BC " << it << " " << ie << " lab=" << r <<  ":\t"
                << K.T[VerticesOfTriangularEdge[ie][0]] << "; "
                << K.T[VerticesOfTriangularEdge[ie][1]] << " E=" << K.T.Edge(ie) << endl;

                for (int k=0;k<kk;k++)
                {
                    gg=R();
                    pair<int,Expression> xx=bc->bc[k];
                    tabexp=0;
                    int comp = xx.first;
                    tabexp[comp]=xx.second;
                    // while  (comp+1 <Nbcomp && which_uh[comp+1] == which_uh[comp])
                    while  (comp+1 <Nbcomp && Uh.dim_which_sub_fem[comp+1] == Uh.dim_which_sub_fem[comp])
                    {  // the right
                        k++; // NEXT COMP
                        comp++;
                        if (k<kk && (comp == bc->bc[k].first) )
                        tabexp[comp]=bc->bc[k].second;
                        else
                        CompileError("In Boundary condition the vector FESpace , we must have:"
                                     " all componant, in the right order");

                    }
                    // cout << " k "<< k << " " << comp << " " << " Nbcomp=" << Nbcomp << " " << Uh.dim_which_sub_fem[comp] << " " << Uh.dim_which_sub_fem[comp+1] <<  endl;
#ifdef OLDPih
                    K.Pi_h(gg,F_Pi_h,buf,&tabexp);

#else
                    K.Pi_h(Aipj);
                    PtonB = 0;
                    for (int i=0;i<Aipj.N();i++)
                    PtonB[ipj[i].p] += onWhatIsEdge[ie][K.DFOnWhat(ipj[i].i)] ;
                    // cout << "   bc->complextype:  " << bc->complextype << endl;
                    for (int p=0;p<PtHat.N();p++)
                    if (PtonB[p]) // in on boundary
                    {
                        mps->set(K.T(PtHat[p]),PtHat[p],K,r,R2(E.y,-E.x)/le,ie); // la normal bofbof ?
                        KN_<R> Vpp(Vp('.',p));
                        Vpp=R();
                        for (int j=0;j<dim;j++)
                        if (tabexp[j])
                        {
                            if(bc->complextype) // FH may 2007  MatriceCreuse
                            Vpp[j]=GetAny<R>( (*tabexp[j])(stack) );
                            else
                            Vpp[j]=GetAny<double>( (*tabexp[j])(stack) );
                        }
                        else Vpp[j]=0.;
                    }
                    //cout << " ..... Vp " << Vp << " " << bc->complextype << " " << bc << endl;
                    for (int i=0;i<Aipj.N();i++)
                    {
                        const FElement::IPJ &ipj_i(ipj[i]);
                        gg[ipj_i.i] += Aipj[i]*Vp(ipj_i.j,ipj_i.p);
                    }
#endif
                    int nbdf = K.NbDoF();
                    for (int df=0;df<nbdf;df++)
                    // if (K.FromFE(df)==which_uh[xx.first] && onWhatIsEdge[ie][K.DFOnWhat(df)] )
                    {
                        //  cout << df << " from = " << K.FromFE(df) << "   dim .. " << Uh.dim_which_sub_fem[xx.first] << "  first " << xx.first << " " << onWhatIsEdge[ie][K.DFOnWhat(df)] << endl;
                        if (K.FromASubFE(df)==Uh.dim_which_sub_fem[xx.first] && onWhatIsEdge[ie][K.DFOnWhat(df)] )
                        {
                            // cout << k << " df=" << df <<  " g= " << gg[df] <<" " << gg(FromTo(0,2)) << endl;
                            int ddf=K(df);
                            // AA(ddf,ddf) =tgv;
                            if (Aii)  onBC[ddf]='1'; ;//A->SetBC(ddf, tgv);// change 21 dec 2010 FH (Hack of ILU)
                            if (B) (*B)[ddf]=  tgv1*gg[df];
                            if (X) (*X)[ddf]=gg[df];
                        }
                    }
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                }
            }
        }
        if( Aii) A->SetBC(onBC,tgv);
        if (! ktbc  && nbon && verbosity )
        {
            cout << " Warning: -- Your set of boundary condition is incompatible with the mesh label." << endl;
        }
        *mps =mp;
    }

    // creating an instance of AssembleBC
    // case 3D volume
    template<class R>
    void AssembleBC(Stack stack,const Mesh3 & Th,const FESpace3 & Uh,const FESpace3 & Vh,bool sym,
                    MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X, const  BC_set * bc, double tgv  )

    {
        typedef Mesh3 Mesh;
        typedef typename FESpace3::FElement FElement;
        typedef typename Mesh::BorderElement BorderElement;
        typedef typename Mesh::Rd Rd;
        typedef typename Mesh::Element Element;
        typedef typename Mesh::RdHat RdHat;

        MeshPoint *mps= MeshPointStack(stack),mp=*mps;
        StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
        bool sptrclean=true;
        //     sptr->clean(); // modif FH mars 2006  clean Ptr

        int ktbc=0, nbon =0;
        bool Aii = A && A->n == A->m;
        int ndofBC = Aii ?  A->n : 1;
        KN<char> onBC(ndofBC);
        onBC= '\0';

        int Nbcomp=Vh.N;
        Check(bc,Nbcomp);
        assert(Vh.N == Uh.N);
        TabFuncArg tabexp(stack,Vh.N);
        KN<double> buf((long int)Vh.MaximalNbOfDF() * (long int)last_operatortype * (long int)Vh.N);
        KN<R> gg(buf);
        if ( B && B->N() != Vh.NbOfDF) ExecError("AssembleBC size rhs and nb of DF of Vh");
        if(verbosity>99) cout << " Problem : BC_set "<< typeid(R).name() << " " ;
        nbon =bc->on.size();
        set<long> on;
        Expandsetoflab(stack,*bc, on);
       /*for (int i=0;i<nbon;i++)
         {
         long  lab  = GetAny<long>( (*bc->on[i])(stack));
         if(verbosity>99) cout << lab << " " ;
         on.insert(lab);
         }
         if(verbosity>99)
         cout << endl;*/
        int kk=bc->bc.size();

        const int dim=Vh.N;

        InterpolationMatrix<RdHat> ipmat(Vh);
        int npPh = Vh.maxNbPtforInterpolation;
        KN<int> PtonB(npPh);
        KNM<R>   Vp(npPh,dim);
        Vp=R();
        KN<R>  Vdf(Vh.MaxNbDFPerElement);
        double tgv1=tgv <0? 1: tgv;
        map<int,int> lll;
        for (int ib=0;ib<Th.nbe;ib++)
        {
            int ie;
            int it = Th.BoundaryElement(ib,ie);

            //const BorderElement &be=Th.be(ib);
            int r =Th.be(ib).lab;
            lll[r]++;
            if (on.find(r) != on.end() )
            {
                const FElement K(Uh[it]);
                ipmat.set(K);

                //R2 E=K.T.Edge(ie);
                //double le = be.mesure();

                ktbc++;
                /*
                 if(verbosity>99)   cout << "BC " << it << " " << ie << " lab=" << r <<  ":\t"
                 << K.T[VerticesOfTriangularEdge[ie][0]] << "; "
                 << K.T[VerticesOfTriangularEdge[ie][1]] << " E=" << K.T.Edge(ie) << endl;
                */
                for (int k=0;k<kk;k++)
                {
                    gg=R();
                    pair<int,Expression> xx=bc->bc[k];
                    tabexp=0;
                    int comp = xx.first;
                    tabexp[comp]=xx.second;
                    // while  (comp+1 <Nbcomp && which_uh[comp+1] == which_uh[comp])
                    while  (comp+1 <Nbcomp && Uh.dim_which_sub_fem[comp+1] == Uh.dim_which_sub_fem[comp])
                    {  // the right
                        k++; // NEXT COMP
                        comp++;
                        if (k<kk && (comp == bc->bc[k].first) )
                        tabexp[comp]=bc->bc[k].second;
                        else
                        CompileError("In Boundary condition the vector FESpace , we must have:"
                                     " all componant, in the right order");

                    }
                    int nbdf=K.NbDoF() ;
                    //ipmat.set(it);
                    PtonB = 0;
                    Rd NN=K.T.N(ie);
                    NN /= NN.norme();
                    for (int i=0;i<ipmat.ncoef;i++)
                    PtonB[ipmat.p[i]] +=  Element::onWhatBorder[ie][K.DFOnWhat(ipmat.dofe[i])] ;


                    for (int p=0;p<ipmat.np;p++)
                    if (PtonB[p]) // in on boundary
                    {
                        const RdHat & PtHat(ipmat.P[p]);
                        mps->set(K.T(PtHat),PtHat,K,r,NN,ie); // la normal bofbof ?
                        KN_<R> Vpp(Vp(p,'.'));
                        for (int j=0;j<dim;j++)
                        if (tabexp[j])
                        if(bc->complextype) // FH may 2007
                        Vpp[j]=GetAny<R>( (*tabexp[j])(stack) );
                        else
                        Vpp[j]=GetAny<double>( (*tabexp[j])(stack) );

                        else Vpp[j]=0.;
                    }
                    // cout << " Vp:  " << Vp << endl;
                    K.Pi_h(Vp,Vdf,ipmat);
                    for (int df=0;df<nbdf;df++)
                    {
                        if (K.FromASubFE(df)==Uh.dim_which_sub_fem[xx.first] && Element::onWhatBorder[ie][K.DFOnWhat(df)] )
                        {
                            int ddf=K(df);
                            // cout << ddf << " " << df << " " << Vdf[df] << " " << it << " ib = " << ib  << " == " << Th(Th[it][df]) <<  endl;
                            //if (Aii)  A->SetBC(ddf,tgv);// change 21 dec 2010 FH (Hack of ILU)
                            if (Aii)  onBC[ddf]='1'; ;//   april 2018 FH
                            if (B) (*B)[ddf]=tgv1*Vdf[df];
                            if (X) (*X)[ddf]=Vdf[df];
                        }
                    }
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                }
            }
        }
        if( Aii) A->SetBC(onBC,tgv);
        if (! ktbc  && nbon && verbosity>1 )
        {
            cout << " Warning: -- Your set of boundary condition is incompatible with the mesh label." << endl;
            if(verbosity>4)
            for (map<int,int>::const_iterator i=lll.begin();i!=lll.end();i++)
                if( on.find(i->first) != on.end() )
                    cout << " on: missing lab " << i-> first << "  nb " << i->second  << endl;
        }
        *mps =mp;
    }

    // creating an instance of AssembleBC
    // case 3D surface
    template<class R>
    void AssembleBC(Stack stack,const MeshS & Th,const FESpaceS & Uh,const FESpaceS & Vh,bool sym,
                    MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X, const  BC_set * bc, double tgv  )
    {
        typedef MeshS Mesh;
        typedef typename FESpaceS::FElement FElement;
        typedef typename Mesh::BorderElement BorderElement;
        typedef typename Mesh::Rd Rd;
        typedef typename Mesh::Element Element;
        typedef typename Mesh::RdHat RdHat;

        MeshPoint *mps= MeshPointStack(stack),mp=*mps;
        StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
        bool sptrclean=true;
        //     sptr->clean(); // modif FH mars 2006  clean Ptr

        int ktbc=0, nbon =0;
        bool Aii = A && A->n == A->m;
        int ndofBC = Aii ?  A->n : 1;
        KN<char> onBC(ndofBC);
        onBC= '\0';

        int Nbcomp=Vh.N;
        Check(bc,Nbcomp);
        assert(Vh.N == Uh.N);
        TabFuncArg tabexp(stack,Vh.N);
        KN<double> buf((long int)Vh.MaximalNbOfDF() * (long int)last_operatortype * (long int)Vh.N);
        KN<R> gg(buf);
        if ( B && B->N() != Vh.NbOfDF) ExecError("AssembleBC size rhs and nb of DF of Vh");
        if(verbosity>99) cout << " Problem : BC_set "<< typeid(R).name() << " " ;
        nbon =bc->on.size();
        set<long> on;
        Expandsetoflab(stack,*bc, on);
        /*
         for (int i=0;i<nbon;i++)
         {
         long  lab  = GetAny<long>( (*bc->on[i])(stack));
         if(verbosity>99) cout << lab << " " ;
         on.insert(lab);
         }
         if(verbosity>99)
         cout << endl;*/
        int kk=bc->bc.size();

        const int dim=Vh.N;

        InterpolationMatrix<RdHat> ipmat(Vh);
        int npPh = Vh.maxNbPtforInterpolation;
        KN<int> PtonB(npPh);
        KNM<R>   Vp(npPh,dim);
        Vp=R();
        KN<R>  Vdf(Vh.MaxNbDFPerElement);
        double tgv1=tgv <0? 1: tgv;
        map<int,int> lll;
        for (int ib=0;ib<Th.nbe;ib++)
        {
            int ie;
            int it = Th.BoundaryElement(ib,ie);

            //const BorderElement &be=Th.be(ib);
            int r =Th.be(ib).lab;
            lll[r]++;
            if (on.find(r) != on.end() )
            {
                const FElement K(Uh[it]);
                ipmat.set(K);

                //R2 E=K.T.Edge(ie);
                //double le = be.mesure();

                ktbc++;
                /*
                 if(verbosity>99)   cout << "BC " << it << " " << ie << " lab=" << r <<  ":\t"
                 << K.T[VerticesOfTriangularEdge[ie][0]] << "; "
                 << K.T[VerticesOfTriangularEdge[ie][1]] << " E=" << K.T.Edge(ie) << endl;
                 */
                for (int k=0;k<kk;k++)
                {
                    gg=R();
                    pair<int,Expression> xx=bc->bc[k];
                    tabexp=0;
                    int comp = xx.first;
                    tabexp[comp]=xx.second;
                    // while  (comp+1 <Nbcomp && which_uh[comp+1] == which_uh[comp])
                    while  (comp+1 <Nbcomp && Uh.dim_which_sub_fem[comp+1] == Uh.dim_which_sub_fem[comp])
                    {  // the right
                        k++; // NEXT COMP
                        comp++;
                        if (k<kk && (comp == bc->bc[k].first) )
                            tabexp[comp]=bc->bc[k].second;
                        else
                            CompileError("In Boundary condition the vector FESpace , we must have:"
                                         " all componant, in the right order");

                    }
                    int nbdf=K.NbDoF() ;
                    //ipmat.set(it);
                    PtonB = 0;

                    R3 E=K.T.Edge(ie);
                    double le = sqrt((E,E));
                    // surface normal
                    Rd NNt=K.T.NormalTUnitaire();
                    // exterior normal (flux)
                    Rd NN=K.T.N(ie);
                    NN /= NN.norme();
                    
                    for (int i=0;i<ipmat.ncoef;i++)
                        PtonB[ipmat.p[i]] +=  Element::onWhatBorder[ie][K.DFOnWhat(ipmat.dofe[i])] ;


                    for (int p=0;p<ipmat.np;p++)
                        if (PtonB[p]) // in on boundary
                        {
                            const RdHat & PtHat(ipmat.P[p]);
                            mps->set(K.T(PtHat),PtHat,K,r,NN,NNt,ie); // la normal bofbof ?
                            KN_<R> Vpp(Vp(p,'.'));
                            for (int j=0;j<dim;j++)
                                if (tabexp[j])
                                    if(bc->complextype) // FH may 2007
                                        Vpp[j]=GetAny<R>( (*tabexp[j])(stack) );
                                    else
                                        Vpp[j]=GetAny<double>( (*tabexp[j])(stack) );

                                    else Vpp[j]=0.;
                        }
                    // cout << " Vp:  " << Vp << endl;
                    K.Pi_h(Vp,Vdf,ipmat);
                    for (int df=0;df<nbdf;df++)
                    {
                        if (K.FromASubFE(df)==Uh.dim_which_sub_fem[xx.first] && Element::onWhatBorder[ie][K.DFOnWhat(df)] )
                        {
                            int ddf=K(df);
                            // cout << ddf << " " << df << " " << Vdf[df] << " " << it << " ib = " << ib  << " == " << Th(Th[it][df]) <<  endl;
                            //if (Aii)  A->SetBC(ddf,tgv);// change 21 dec 2010 FH (Hack of ILU)
                            if (Aii)  onBC[ddf]='1'; ;//   april 2018 FH
                            if (B) (*B)[ddf]=tgv1*Vdf[df];
                            if (X) (*X)[ddf]=Vdf[df];
                        }
                    }
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                }
            }
        }
        if( Aii) A->SetBC(onBC,tgv);
        if (! ktbc  && nbon && verbosity>1 )
        {
            cout << " Warning: -- Your set of boundary condition is incompatible with the mesh label." << endl;
            if(verbosity>9)
            for (map<int,int>::const_iterator i=lll.begin();i!=lll.end();i++)
              if( on.find(i->first) != on.end() )
                  cout << " on: missing lab " << i-> first << "  nb " << i->second  << endl;
        }
        *mps =mp;
    }

    // creating an instance of AssembleBC
    // case 3D curve
    template<class R>
    void AssembleBC(Stack stack,const MeshL & Th,const FESpaceL & Uh,const FESpaceL & Vh,bool sym,
                    MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X, const  BC_set * bc, double tgv  )
    {
        typedef MeshL Mesh;
        typedef typename FESpaceL::FElement FElement;
        typedef typename Mesh::BorderElement BorderElement;
        typedef typename Mesh::Rd Rd;
        typedef typename Mesh::Element Element;
        typedef typename Mesh::RdHat RdHat;

        MeshPoint *mps= MeshPointStack(stack),mp=*mps;
        StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
        bool sptrclean=true;
        //     sptr->clean(); // modif FH mars 2006  clean Ptr

        int ktbc=0, nbon =0;
        bool Aii = A && A->n == A->m;
        int ndofBC = Aii ?  A->n : 1;
        KN<char> onBC(ndofBC);
        onBC= '\0';

        int Nbcomp=Vh.N;
        Check(bc,Nbcomp);
        assert(Vh.N == Uh.N);
        TabFuncArg tabexp(stack,Vh.N);
        KN<double> buf((long int)Vh.MaximalNbOfDF() * (long int)last_operatortype * (long int)Vh.N);
        KN<R> gg(buf);
        if ( B && B->N() != Vh.NbOfDF) ExecError("AssembleBC size rhs and nb of DF of Vh");
        if(verbosity>99) cout << " Problem : BC_set "<< typeid(R).name() << " " ;
        nbon =bc->on.size();
        set<long> on;
        Expandsetoflab(stack,*bc, on);

        int kk=bc->bc.size();

        const int dim=Vh.N;

        InterpolationMatrix<RdHat> ipmat(Vh);
        int npPh = Vh.maxNbPtforInterpolation;
        KN<int> PtonB(npPh);
        KNM<R>   Vp(npPh,dim);
        Vp=R();
        KN<R>  Vdf(Vh.MaxNbDFPerElement);
        double tgv1=tgv <0? 1: tgv;
        map<int,int> lll;
        for (int ib=0;ib<Th.nbe;ib++) {
            int ie;
            int it = Th.BoundaryElement(ib,ie);

            int r =Th.be(ib).lab;
            lll[r]++;
            if (on.find(r) != on.end() ) {
                const FElement K(Uh[it]);
                ipmat.set(K);
                 ktbc++;

                for (int k=0;k<kk;k++) {
                    gg=R();
                    pair<int,Expression> xx=bc->bc[k];
                    tabexp=0;
                    int comp = xx.first;
                    tabexp[comp]=xx.second;
                    // while  (comp+1 <Nbcomp && which_uh[comp+1] == which_uh[comp])
                    while  (comp+1 <Nbcomp && Uh.dim_which_sub_fem[comp+1] == Uh.dim_which_sub_fem[comp])
                    {  // the right
                        k++; // NEXT COMP
                        comp++;
                        if (k<kk && (comp == bc->bc[k].first) )
                            tabexp[comp]=bc->bc[k].second;
                        else
                            CompileError("In Boundary condition the vector FESpace , we must have:"
                                         " all componant, in the right order");

                    }
                    int nbdf=K.NbDoF() ;
                    //ipmat.set(it);
                    PtonB = 0;
                    R3 NNt=K.T.NormalTUnitaire();
                    // exterior normal (flux)
                    Rd NN=K.T.N(ie);
                    NN /= NN.norme();
                    
                    for (int i=0;i<ipmat.ncoef;i++)
                        PtonB[ipmat.p[i]] +=  Element::onWhatBorder[ie][K.DFOnWhat(ipmat.dofe[i])] ;


                    for (int p=0;p<ipmat.np;p++)
                        if (PtonB[p]) // in on boundary
                        {
                            const RdHat & PtHat(ipmat.P[p]);
                            mps->set(K.T(PtHat),PtHat,K,r,NN,NNt,ie);
                            KN_<R> Vpp(Vp(p,'.'));
                            for (int j=0;j<dim;j++)
                                if (tabexp[j])
                                    if(bc->complextype) // FH may 2007
                                        Vpp[j]=GetAny<R>( (*tabexp[j])(stack) );
                                    else
                                        Vpp[j]=GetAny<double>( (*tabexp[j])(stack) );

                                    else Vpp[j]=0.;
                        }
                    K.Pi_h(Vp,Vdf,ipmat);
                    for (int df=0;df<nbdf;df++)
                    {
                        if (K.FromASubFE(df)==Uh.dim_which_sub_fem[xx.first] && Element::onWhatBorder[ie][K.DFOnWhat(df)] )
                        {
                            int ddf=K(df);
                            if (Aii)  onBC[ddf]='1'; ;//   april 2018 FH
                            if (B) (*B)[ddf]=tgv1*Vdf[df];
                            if (X) (*X)[ddf]=Vdf[df];
                        }
                    }
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                }
            }
        }
        if( Aii) A->SetBC(onBC,tgv);
        if (! ktbc  && nbon && verbosity )
        {
            cout << " Warning: -- Your set of boundary condition is incompatible with the mesh label." << endl;
            for (map<int,int>::const_iterator i=lll.begin();i!=lll.end();i++)
                cout << " lab " << i-> first << "  nb " << i->second  << endl;
        }
        *mps =mp;
    }

   // case 3D curve / 2D on meshL
   template<class R>
   void AssembleBC(Stack stack,const MeshL & Th,const FESpaceL & Uh,const FESpace & Vh,bool sym,
                   MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X, const  BC_set * bc, double tgv  )
   {
       ffassert(0);
   }

   // case 2D / 3D curve on meshL
   template<class R>
   void AssembleBC(Stack stack,const MeshL & Th,const FESpace & Uh,const FESpaceL & Vh,bool sym,
                   MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X, const  BC_set * bc, double tgv  )
   {
       ffassert(0);
   }
    // case 3D Surf / 3D volume on meshS
    template<class R>
    void AssembleBC(Stack stack,const MeshS & Th,const FESpaceS & Uh,const FESpace3 & Vh,bool sym,
                   MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X, const  BC_set * bc, double tgv  )
    {
       ffassert(0);
    }
    // case 3D volume / 3D Surf on meshS
    template<class R>
    void AssembleBC(Stack stack,const MeshS & Th,const FESpace3 & Uh,const FESpaceS & Vh,bool sym,
                   MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X, const  BC_set * bc, double tgv  )
    {
       ffassert(0);
    }
    // case 3D curve / 3D Surf on meshL
    template<class R>
    void AssembleBC(Stack stack,const MeshL & Th,const FESpaceL & Uh,const FESpaceS & Vh,bool sym,
                   MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X, const  BC_set * bc, double tgv  )
    {
       ffassert(0);
    }
    // case 3D Surf / 3D curve on meshL
    template<class R>
    void AssembleBC(Stack stack,const MeshL & Th,const FESpaceS & Uh,const FESpaceL & Vh,bool sym,
                   MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X, const  BC_set * bc, double tgv  )
    {
       ffassert(0);
    }

    void  Expandsetoflab(Stack stack,const BC_set & bc,set<long> & setoflab);
    void  Expandsetoflab(Stack stack,const CDomainOfIntegration & di,set<int> & setoflab,bool &all);

    //////////////////////////////////
    // AssembleLinearForm
    //////////////////////////////////

    // creating an instance of AssembleLinearForm
    // case 2d
    template<class R>
    void AssembleLinearForm(Stack stack,const Mesh & Th,const FESpace & Vh,KN_<R> * B,const  FormLinear * l )
    {
        StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
        bool sptrclean=true;
        //     sptr->clean(); // modif FH mars 2006  clean Ptr
        Check(l->l,Vh.N);
        if ( B && B->N() != Vh.NbOfDF) ExecError("AssembleLinearForm size rhs and nb of DF of Vh");
        // if ( & Th != &Vh.Th ) ExecError("AssembleLinearForm on different meshes  ( not implemented FH).");
        KN<double> buf(Vh.MaximalNbOfDF()*last_operatortype*Vh.N*2);

        //            const  FormLinear * l=dynamic_cast<const  FormLinear *>(e);
        const CDomainOfIntegration & di= *l->di;
        const Mesh & ThI = Th;// * GetAny<pmesh>( (* di.Th)(stack));
        bool sameMesh = &ThI == &Vh.Th;
        const bool intmortar=di.intmortar(stack);

        SHOWVERB(cout << " FormLinear " << endl);
        // const vector<Expression>  & what(di.what);

        CDomainOfIntegration::typeofkind  kind = di.kind;
        const QuadratureFormular1d & FIE = di.FIE(stack);
        const QuadratureFormular & FIT = di.FIT(stack);
        const int useopt=di.UseOpt(stack);
        double binside=di.binside(stack);  // truc FH pour fluide de grad2 (decentrage bizard)
        //  cout << "AssembleLinearForm " << l->l->v.size() << endl;
        set<int> setoflab;
        bool all=true;
        bool VF=l->VF();  // finite Volume or discontinuous Galerkin
        if (verbosity>2) cout << "  -- AssembleLinearForm 2, discontinuous Galerkin  =" << VF << " binside = "<< binside
        << " levelset integration " <<di.islevelset()<< " withmap: "<<  di.withmap() << "\n";
        //  if( di.withmap()) { ExecError(" no map  in the case (6)??");}
        Expression  const * const mapt=di.mapt[0] ? di.mapt:0;
        sameMesh = sameMesh && !mapt; //

        if (verbosity>3)
        {

            if (CDomainOfIntegration::int1d==kind) cout << "  -- boundary int border ( nQP: "<< FIE.n << ") ";
            else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all edges ( nQP: "<< FIE.n << "),"  ;
            else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF edges nQP: ("<< FIE.n << ")," ;
            else cout << "  --  int 2d  (nQP: "<< FIT.n << " ) in "  ;
            cout << ", samemesh :"<< sameMesh<< " int mortar: " << intmortar ;
        }
        /*
         if ( verbosity>3)
         if (kind==CDomainOfIntegration::int1d) cout << "  -- boundary int border " ;
         else if (kind==CDomainOfIntegration::intalledges) cout << "  -- boundary int all edges " ;
         else if (kind==CDomainOfIntegration::intallVFedges) cout << "  -- boundary int all edges " ;
         else cout << "  -- boundary int  " ;
         */
        if(di.islevelset() && ( (CDomainOfIntegration::int1d!=kind) && (CDomainOfIntegration::int2d!=kind) )  )
        InternalError("So no levelset integration type on no int1d/int2d case (4)");
        Expandsetoflab(stack,di, setoflab,all);
        /*
         for (size_t i=0;i<what.size();i++)
         {long  lab  = GetAny<long>( (*what[i])(stack));
         setoflab.insert(lab);
         if ( verbosity>3) cout << lab << " ";
         all=false;
         } */
        if (verbosity>3) cout << " Optimized = "<< useopt << ", ";

        const E_F0 * poptiexp0=l->l->optiexp0;
        // const E_F0 & optiexpK=*l->l->optiexpK;
        int n_where_in_stack_opt=l->l->where_in_stack_opt.size();
        R** where_in_stack =0;
        if (n_where_in_stack_opt && useopt)
        where_in_stack = new R * [n_where_in_stack_opt];
        if (where_in_stack)
        {
            assert(l->l->v.size()==(size_t) n_where_in_stack_opt);
            for (int i=0;i<n_where_in_stack_opt;i++)
            {
                int offset=l->l->where_in_stack_opt[i];
                assert(offset>10);
                where_in_stack[i]= static_cast<R *>(static_cast<void *>((char*)stack+offset));
                *(where_in_stack[i])=0;
            }
            if(poptiexp0) (*poptiexp0)(stack);

            if( (verbosity/100) && verbosity % 10 == 2)
            {
                int il=0;

                for (LinearOperatorD::const_iterator ll=l->l->v.begin();ll!=l->l->v.end();ll++,il++)
                cout << il << " coef (" << ll->first << ") = " << *(where_in_stack[il]) << " offset=" << l->l->where_in_stack_opt[il] <<endl;

                for (int i=0;i<n_where_in_stack_opt;i++)
                cout << "const coef " << i << " = " << *(where_in_stack[i]) << endl;
            }
        }
        Stack_Ptr<R*>(stack,ElemMatPtrOffset) =where_in_stack;

        KN<int>   ip(Vh.MaxNbDFPerElement*6);
        if (verbosity >3)
        {
            if (all) cout << " all " << endl ;
            else cout << endl;
        }
        if(di.islevelset() && (kind !=CDomainOfIntegration::int1d)&& (kind !=CDomainOfIntegration::int2d))
        InternalError(" Sorry No levelSet integral for is case ..(5)");


        if (kind==CDomainOfIntegration::int1d)
        {


            if(VF) InternalError(" no jump or average in int1d of RHS");
            if(di.islevelset())
            {
                double uset = HUGE_VAL;
                R2 Q[3];
                KN<double> phi(ThI.nv);phi=uset;
                double f[3];
                for(int t=0; t< ThI.nt;++t)
                {
                    double umx=-HUGE_VAL,umn=HUGE_VAL;
                    for(int i=0;i<3;++i)
                    {
                        int j= ThI(t,i);
                        if( phi[j]==uset)
                        {
                            MeshPointStack(stack)->setP(&ThI,t,i);
                            phi[j]= di.levelset(stack);//zzzz
                        }
                        f[i]=phi[j];
                        umx = std::max(umx,phi[j]);
                        umn = std::min(umn,phi[j]);

                    }
                    if( umn <=0 && umx >= 0)
                    {

                        int np= IsoLineK(f,Q,1e-10);
                        if(np==2)
                        {
                            if ( sameMesh )
                            {/*
                              void  Element_rhs(const FElement & Kv,const LOperaD &Op,double * p,void * stack,KN_<R> & B,
                              const QuadratureFormular1d & FI ,const R2 & PA,const R2 &PB)

                              */
                                Element_rhs<R>(Vh[t],*l->l,buf,stack,*B,FIE,Q[0],Q[1],useopt);
                            }
                            else if(!mapt)
                            Element_rhs<R>(ThI,ThI[t],Vh,0,ThI[t].lab,*l->l,buf,stack,*B,FIE,false,intmortar,Q,useopt);
                            else
                            Element_rhs<R>(mapt,ThI,ThI[t],Vh,0,ThI[t].lab,*l->l,buf,stack,*B,FIE,false,intmortar,Q,useopt);

                            //InternalError(" No levelSet on Diff mesh :    to day  int1d of RHS");
                        }
                        if(sptrclean) sptrclean=sptr->clean();
                    }
                }

            }
            else
            for( int e=0;e<ThI.neb;e++)
            {
                if (all || setoflab.find(ThI.bedges[e].lab) != setoflab.end())
                {
                    int ie,i =ThI.BoundaryElement(e,ie);
                    if ( sameMesh )
                    Element_rhs<R>(Vh[i],ie,Th.bedges[e].lab,*l->l,buf,stack,*B,FIE,false,useopt);
                    else if(!mapt)
                    Element_rhs<R>(ThI,ThI[i],Vh,ie,Th.bedges[e].lab,*l->l,buf,stack,*B,FIE,false,intmortar,0,useopt);
                    else
                    Element_rhs<R>(mapt,ThI,ThI[i],Vh,ie,Th.bedges[e].lab,*l->l,buf,stack,*B,FIE,false,intmortar,0,useopt);
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                }
            }
        }
        else if (kind==CDomainOfIntegration::intalledges)
        {
            ffassert(mapt==0);
            if(VF)
            {
                pair_stack_double bstack(stack,& binside);

                //bstack.first = stack;
                //bstack.second= & binside;

                //InternalError(" Today no jump or average in intalledges of RHS ");
                for (int i=0;i< ThI.nt; i++)
                if (all || setoflab.find(ThI[i].lab) != setoflab.end())
                {

                    for (int ie=0;ie<3;ie++)
                    if ( sameMesh)
                    {
                        int iie=ie,ii=Th.ElementAdj(i,iie);
                        if(ii<0) ii=i;//  sur le bord
                        const Triangle & K(ThI[i]);
                        int e0=VerticesOfTriangularEdge[ie][0];
                        int e1=VerticesOfTriangularEdge[ie][1];
                        int i1 = ThI(K[e0]),i2 = ThI(K[e1]);
                        BoundaryEdge * be = ThI.TheBoundaryEdge(i1,i2);
                        int lab = be ? be->lab :  notalabel;

                        Element_rhsVF<R>(Vh[i],Vh[ii],ie,iie,lab,*l->l,buf,ip,&bstack,*B,FIE,useopt);
                    }
                    else
                    InternalError("To Do") ;
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                }

            }
            else
            for (int i=0;i< ThI.nt; i++)
            if (all || setoflab.find(ThI[i].lab) != setoflab.end())
            {
                for (int ie=0;ie<3;ie++)
                if ( sameMesh)
                Element_rhs<R>(Vh[i],ie,Th[i].lab,*l->l,buf,stack,*B,FIE,true,useopt);
                else
                InternalError("To Do") ;
                if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
            }
        }
        else if (kind==CDomainOfIntegration::intallVFedges)
        {
            cerr << " intallVFedges a faire" << endl;

            InternalError(" intallVFedges a faire ");

            ffassert(0);
            for (int i=0;i< ThI.nt; i++)
            {
                if (all || setoflab.find(ThI[i].lab) != setoflab.end())
                for (int ie=0;ie<3;ie++)
                {
                    if ( sameMesh)
                    Element_rhs<R>(Vh[i],ie,Th[i].lab,*l->l,buf,stack,*B,FIE,true,useopt);
                    else
                    InternalError("To Do") ;
                }
                if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
            }
        }

        else if (kind==CDomainOfIntegration::int2d){
            if(di.islevelset())
            {
                QuadratureFormular FITM(FIT);
                double uset = HUGE_VAL;
                R2 Q[4];
                KN<double> phi(Th.nv);phi=uset;
                double f[3];
                for(int t=0; t< Th.nt;++t)
                {
                    if ( all || setoflab.find(ThI[t].lab) != setoflab.end())
                    {
                        double umx=-HUGE_VAL,umn=HUGE_VAL;
                        for(int i=0;i<3;++i)
                        {
                            int j= ThI(t,i);
                            if( phi[j]==uset)
                            {
                                MeshPointStack(stack)->setP(&ThI,t,i);
                                phi[j]= di.levelset(stack);//zzzz
                            }
                            f[i]=phi[j];
                            umx = std::max(umx,phi[j]);
                            umn = std::min(umn,phi[j]);

                        }
                        if( umx <=0 )
                        Element_rhs<R>(Vh[t],*l->l,buf,stack,*B,FIT,useopt);
                        else if( umn <0 )
                        { // coupe ..
                            int i0 = 0, i1 = 1, i2 =2;

                            if( f[i0] > f[i1] ) swap(i0,i1) ;
                            if( f[i0] > f[i2] ) swap(i0,i2) ;
                            if( f[i1] > f[i2] ) swap(i1,i2) ;

                            double c = (f[i2]-f[i1])/(f[i2]-f[i0]); // coef Up Traing
                            if( f[i1] < 0 ) {double y=f[i2]/(f[i2]-f[i1]); c *=y*y; }
                            else {double y=f[i0]/(f[i0]-f[i1]) ; c = 1.- (1.-c)*y*y; };
                            assert( c > 0 && c < 1);
                            double arean = (1-c)*Th[t].area;
                            FITM=FIT;
                            FITM*=1-c;
                            ffassert(mapt==0);
                            Element_rhs<R>(Vh[t],*l->l,buf,stack,*B,FITM,useopt);
                        }
                        if(sptrclean) sptrclean=sptr->clean();
                    }
                }
            }
            else
            for (int i=0;i< ThI.nt; i++)
            if (all || setoflab.find(ThI[i].lab) != setoflab.end())
            {
                if ( sameMesh )
                Element_rhs<R>(Vh[i],*l->l,buf,stack,*B,FIT,useopt);
                else if(!mapt)
                Element_rhs<R>(ThI,ThI[i],Vh,*l->l,buf,stack,*B,FIT,useopt);
                else
                Element_rhs<R>(mapt,ThI,ThI[i],Vh,*l->l,buf,stack,*B,FIT,useopt);
                if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
            }
        }

        if (n_where_in_stack_opt) delete [] where_in_stack;

    }


    // creating an instance of AssembleLinearForm
    // case 3D volume
    template<class R>
    void AssembleLinearForm(Stack stack,const Mesh3 & Th,const FESpace3 & Vh,KN_<R> * B,const  FormLinear * l )
    {
        typedef FESpace3 FESpace;
        typedef FESpace3::Mesh Mesh;
        typedef Mesh *pmesh ;
        typedef Mesh::BorderElement BorderElement;
        StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
        bool sptrclean=true;
        //     sptr->clean(); // modif FH mars 2006  clean Ptr
        Check(l->l,Vh.N);
        if ( B && B->N() != Vh.NbOfDF) ExecError("AssembleLinearForm size rhs and nb of DF of Vh");
        // if ( & Th != &Vh.Th ) ExecError("AssembleLinearForm on different meshes  ( not implemented FH).");
        KN<double> buf(Vh.MaximalNbOfDF()*last_operatortype*Vh.N*2);

        //            const  FormLinear * l=dynamic_cast<const  FormLinear *>(e);
        const CDomainOfIntegration & di= *l->di;
        ffassert(di.d==3);
        // const Mesh * pThdi = GetAny<pmesh>( (* di.Th)(stack));

        const Mesh & ThI = Th;// * GetAny<pmesh>( (* di.Th)(stack));
        bool sameMesh = &ThI == &Vh.Th;

        SHOWVERB(cout << " FormLinear " << endl);
        //const vector<Expression>  & what(di.what);

        CDomainOfIntegration::typeofkind  kind = di.kind;
        //const QuadratureFormular1d & FIE = di.FIE(stack);
        //  const QuadratureFormular & FIT = di.FIT(stack);
        // const GQuadratureFormular<R3> & FIV = di.FIV(stack);

        // const QuadratureFormular1d & FIEo = di.FIE(stack);
        const QuadratureFormular & FITo = di.FIT(stack);
        const GQuadratureFormular<R3> & FIVo = di.FIV(stack);
        //  to change the quadrature on element ... may 2014 FH ..
        // QuadratureFormular1d  FIE(FIEo,3);
        QuadratureFormular FIT(FITo,3);
        GQuadratureFormular<R3>  FIV(FIVo,3);

        const int useopt=di.UseOpt(stack);
        double binside=di.binside(stack);  // truc FH pour fluide de grad2 (decentrage bizard)
        //  cout << "AssembleLinearForm " << l->l->v.size() << endl;
        set<int> setoflab;
        bool all=true;
        bool VF=l->VF();  // finite Volume or discontinuous Galerkin
        if (verbosity>2) cout << "  -- AssembleLinearForm 1,  discontinuous Galerkin  =" << VF << " binside = "<< binside <<"\n";

        if (verbosity>3)
        {
            if (CDomainOfIntegration::int2d==kind) cout << "  -- boundary int border ( nQP: "<< FIT.n << ") , samemesh: " << sameMesh << " "   ;
            else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all edges ( nQP: "<< FIT.n << "),"  ;
            else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF edges nQP: ("<< FIT.n << ")," ;
            else cout << "  --  int 3d   (nQP: "<< FIV.n << " ) in "  ;
        }
        if( di.withmap()) { ExecError(" no map  in the case (5)??");}

        //  if(di.islevelset()) InternalError("So no levelset integration type on this case (3)");
        if(di.islevelset() && (CDomainOfIntegration::int2d!=kind) && (CDomainOfIntegration::int3d!=kind) )
        InternalError("So no levelset intgeration type on no int2d/3d case");
        /*
         if ( verbosity>3)
         if (kind==CDomainOfIntegration::int1d) cout << "  -- boundary int border " ;
         else if (kind==CDomainOfIntegration::intalledges) cout << "  -- boundary int all edges " ;
         else if (kind==CDomainOfIntegration::intallVFedges) cout << "  -- boundary int all edges " ;
         else cout << "  -- boundary int  " ;
         */

        Expandsetoflab(stack,di, setoflab,all);
        /*
         for (size_t i=0;i<what.size();i++)
         if(di.whatis[i] ==0)
         {
         long  lab  = GetAny<long>( (*what[i])(stack));
         setoflab.insert(lab);
         if ( verbosity>3) cout << lab << " ";
         all=false;
         }
         else
         {
         KN<long>  labs( GetAny<KN_<long> >( (*what[i])(stack)));
         for (long j=0; j<labs.N(); ++j) {
         setoflab.insert(labs[j]);
         if ( verbosity>3) cout << labs[j] << " ";
         }
         all=false;
         }*/

        if (verbosity>3) cout << " Optimized = "<< useopt << ", ";

        const E_F0 * poptiexp0=l->l->optiexp0;
        // const E_F0 & optiexpK=*l->l->optiexpK;
        int n_where_in_stack_opt=l->l->where_in_stack_opt.size();
        R** where_in_stack =   0;
        if (n_where_in_stack_opt && useopt)
        where_in_stack = new R * [n_where_in_stack_opt];
        if (where_in_stack)
        {
            assert(l->l->v.size()==(size_t) n_where_in_stack_opt);
            for (int i=0;i<n_where_in_stack_opt;i++)
            {
                int offset=l->l->where_in_stack_opt[i];
                assert(offset>10);
                where_in_stack[i]= static_cast<R *>(static_cast<void *>((char*)stack+offset));
                *(where_in_stack[i])=0;
            }
            if(poptiexp0) (*poptiexp0)(stack);

            if( (verbosity/100) && verbosity % 10 == 2)
            {
                int il=0;

                for (LinearOperatorD::const_iterator ll=l->l->v.begin();ll!=l->l->v.end();ll++,il++)
                cout << il << " coef (" << ll->first << ") = " << *(where_in_stack[il]) << " offset=" << l->l->where_in_stack_opt[il] <<endl;

                for (int i=0;i<n_where_in_stack_opt;i++)
                cout << "const coef " << i << " = " << *(where_in_stack[i]) << endl;
            }
        }
        Stack_Ptr<R*>(stack,ElemMatPtrOffset) =where_in_stack;

        KN<int>   ip(Vh.MaxNbDFPerElement*6);
        if (verbosity >3)
        {
            if (all) cout << " all " << endl ;
            else cout << endl;
        }
        if (kind==CDomainOfIntegration::int2d)
        { //AFAIRE("3D Element RHS CDomainOfIntegration::int2d");
            double  ss =0;
            if(VF) InternalError(" no jump or average in int1d of RHS");
            if(di.islevelset()) // init on level set (of RHS)
            {
                double uset = HUGE_VAL;
                R3 Q[4];
                KN<double> phi(ThI.nv);phi=uset;
                double f[4];
                for(int t=0; t< ThI.nt;++t)
                {

                    double umx=-HUGE_VAL,umn=HUGE_VAL;
                    for(int i=0;i<4;++i)
                    {
                        int j= ThI(t,i);
                        if( phi[j]==uset)
                        {
                            MeshPointStack(stack)->setP(&ThI,t,i);
                            phi[j]= di.levelset(stack);//zzzz
                        }
                        f[i]=phi[j];
                        umx = std::max(umx,phi[j]);
                        umn = std::min(umn,phi[j]);

                    }
                    if( umn <=0 && umx >= 0)
                    {

                        int np= IsoLineK(f,Q,1e-10);// ca code ...
                        if(np==3 || np==4)
                        {  //  if(np==3) Q[3]=Q[0]; // same 0 == 3 bofbof ??? FH
                            //   cout << " Q[0]" << Q[0] << endl;
                            if( verbosity> 99)
                            {
                                R3 PP[4];
                                const Tet  &K(ThI[t]);
                                for(int i=0; i< np; ++i)
                                PP[i]= K(Q[i]);
                                for( int i =0; i+1 < np; i+=2)
                                {
                                    int i0=i,i1=i+1,i2=(i+2)%np;
                                    R3 NN= R3(PP[i0],PP[i1])^R3(PP[i0],PP[i2]);
                                    double mes2 = (NN,NN);
                                    double mes = sqrt(mes2)/2;
                                    ss+= mes;
                                    //cout << "mes " << mes << " " << i << " , ";
                                }
                            }

                            if ( sameMesh)
                            Element_rhs<R>(Vh[t],*l->l,buf,stack,*B,FIT,np,Q,useopt);
                            else
                            //    else
                            InternalError(" No levelSet on Diff mesh3 :    to day  int2d of RHS");
                            //    Element_rhs<R>(ThI,ThI[t],Vh,-1,lab,*l->l,buf,stack,*B,FIT,false);
                        }
                        if(sptrclean) sptrclean=sptr->clean();
                    }
                }
                if( verbosity> 99)
                cout << "          surf levelset = " << ss << endl;

            }
            else
            for( int e=0;e<ThI.nbe;e++)
            {
                if (all || setoflab.find(ThI.be(e).lab) != setoflab.end())
                {
                    int ie,i =ThI.BoundaryElement(e,ie);
                    if ( sameMesh)
                    Element_rhs<R>(Vh[i],ie,Th.be(e).lab,*l->l,buf,stack,*B,FIT,false,useopt);
                    else
                    Element_rhs<R>(ThI,ThI[i],Vh,ie,Th.be(e).lab,*l->l,buf,stack,*B,FIT,false,useopt);
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                }
            }
        }
        else if (kind==CDomainOfIntegration::intalledges)
        {     InternalError("3D Element RHS CDomainOfIntegration::intalledges :  stupide !!!");
        }
        else if (kind==CDomainOfIntegration::intallVFedges)
        {
            cerr << " intallVFedges a faire" << endl;

            InternalError(" intallVFedges a stupide!!! ");

            ffassert(0);
        }

        else if(kind==CDomainOfIntegration::int3d) {
            if(di.islevelset())  //  may 2014 FH ...
            {   // int3d levelset < 0
                double llevelset = 0;
                const double uset = std::numeric_limits<double>::max();
                // cout << " uset ="<<uset << endl;
                R3 Q[3][4];
                double vol6[3];
                KN<double> phi(Th.nv);
                phi=uset;
                double f[4];

                for (int t=0;t< Th.nt; t++)
                {

                    const Mesh3::Element & K(ThI[t]);
                    if (all || setoflab.find(ThI[t].lab) != setoflab.end())

                    {
                        double umx=std::numeric_limits<double>::lowest(),umn=std::numeric_limits<double>::max();
                        for(int i=0;i<4;++i)
                        {
                            int j= ThI(t,i);
                            if( phi[j]==uset)
                            {
                                MeshPointStack(stack)->setP(&ThI,t,i);
                                phi[j]= di.levelset(stack);//zzzz
                            }
                            f[i]=phi[j];
                        }
                        int ntets= UnderIso(f,Q, vol6,1e-14);
                        setQF<R3>(FIV,FIVo,QuadratureFormular_Tet_1, Q,vol6,ntets);
                        if(FIV.n)
                        {
                            if ( sameMesh )
                            Element_rhs<R>(Vh[t],*l->l,buf,stack,*B,FIV,useopt);
                            else
                            Element_rhs<R>(ThI,ThI[t],Vh,*l->l,buf,stack,*B,FIV,useopt);
                            if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr

                        }
                    }
                }
                FIV=FIVo;
            }
            else
            {

                for (int i=0;i< ThI.nt; i++)
                if (all || setoflab.find(ThI[i].lab) != setoflab.end())
                {
                    if ( sameMesh )
                    Element_rhs<R>(Vh[i],*l->l,buf,stack,*B,FIV,useopt);
                    else
                    Element_rhs<R>(ThI,ThI[i],Vh,*l->l,buf,stack,*B,FIV,useopt);
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                }}
        }
        else  if(kind==CDomainOfIntegration::intallfaces    ) {

            if(VF)
            {
                if ( !sameMesh) InternalError(" no jump or average in intallfaces for RHS in not samemesh");
             //   InternalError(" no jump or average in intallfaces for RHS");
                pair_stack_double bstack(stack,& binside);

                for (int i=0;i< ThI.nt; i++)
                if (all || setoflab.find(ThI[i].lab) != setoflab.end())
                {

                    for (int ie=0;ie<4;ie++)
                   
                    {
                        int iie=ie,ii=Th.ElementAdj(i,iie);
                        BorderElement * be = 0; // FIND BOUNDARY ELEMENT !!!
                        if(ii<0) ii=i;//  sur le bord
                        const Tet & K(ThI[i]);
                        
                        int lab = be ? be->lab :  notalabel;
                        Element_rhsVF<R>(Vh[i],Vh[ii],ie,iie,lab,*l->l,buf,ip,&bstack,*B,FIT,useopt);
                    }
                    
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
    
                }
                
            }
            else
            {
            for(int i=0;i<ThI.nt; i++)
              for(int ie=0;ie<Mesh3::nea; ie++)
               {
                int lab=0;
                // if face on bord get the lab ???
                if ( sameMesh)
                Element_rhs<R>(Vh[i],ie,lab,*l->l,buf,stack,*B,FIT,false,useopt);
                else
                Element_rhs<R>(ThI,ThI[i],Vh,ie,lab,*l->l,buf,stack,*B,FIT,false,useopt);
                if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr

               }
            }
        }
        else
        {
            cout << " Strange (unknows) kind = " << kind << endl;
            ffassert(0);
        }
        if (n_where_in_stack_opt) delete [] where_in_stack;

    }







// creating an instance of AssembleLinearForm
// case surface 3d
template<class R>
void AssembleLinearForm(Stack stack,const MeshS & Th,const FESpaceS & Vh,KN_<R> * B,const  FormLinear * l )
    {

        StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
        bool sptrclean=true;
        //     sptr->clean(); // modif FH mars 2006  clean Ptr
        Check(l->l,Vh.N);
        if ( B && B->N() != Vh.NbOfDF) ExecError("AssembleLinearForm size rhs and nb of DF of Vh");
        // if ( & Th != &Vh.Th ) ExecError("AssembleLinearForm on different meshes  ( not implemented FH).");
        KN<double> buf(Vh.MaximalNbOfDF()*last_operatortype*Vh.N*2);

        //            const  FormLinear * l=dynamic_cast<const  FormLinear *>(e);
        const CDomainOfIntegration & di= *l->di;
        ffassert(di.d==3);

        const MeshS & ThI = Th;// * GetAny<pmesh>( (* di.Th)(stack));
        bool sameMesh = &ThI == &Vh.Th;

        const bool intmortar=di.intmortar(stack);

        SHOWVERB(cout << " FormLinear " << endl);
        // const vector<Expression>  & what(di.what);

        CDomainOfIntegration::typeofkind  kind = di.kind;
        const QuadratureFormular1d & FIE = di.FIE(stack);
        const QuadratureFormular & FIT = di.FIT(stack);
        const int useopt=di.UseOpt(stack);
        double binside=di.binside(stack);  // truc FH pour fluide de grad2 (decentrage bizard)
        //  cout << "AssembleLinearForm " << l->l->v.size() << endl;
        set<int> setoflab;
        bool all=true;
        bool VF=l->VF();  // finite Volume or discontinuous Galerkin

        if (verbosity>2) cout << "  -- AssembleLinearForm S, discontinuous Galerkin  =" << VF << " binside = "<< binside
            << " levelset integration " <<di.islevelset()<< " withmap: "<<  di.withmap() << "\n";
        //  if( di.withmap()) { ExecError(" no map  in the case (6)??");}
        Expression  const * const mapt=di.mapt[0] ? di.mapt:0;
        sameMesh = sameMesh && !mapt; //

        if (verbosity>3)
        {

            if (CDomainOfIntegration::int1d==kind) cout << "  -- boundary int border ( nQP: "<< FIE.n << ") ";
            else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all edges ( nQP: "<< FIE.n << "),"  ;
            else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF edges nQP: ("<< FIE.n << ")," ;
            else cout << "  --  int 2d  (nQP: "<< FIT.n << " ) in "  ;
            cout << ", samemesh :"<< sameMesh<< " int mortar: " << intmortar ;
        }
         if(di.islevelset() && ( (CDomainOfIntegration::int1d!=kind) && (CDomainOfIntegration::int2d!=kind) )  )
            InternalError("So no levelset integration type on no int1d/int2d case (4)");
        Expandsetoflab(stack,di, setoflab,all);
         if (verbosity>3) cout << " Optimized = "<< useopt << ", ";
        const E_F0 * poptiexp0=l->l->optiexp0;
        // const E_F0 & optiexpK=*l->l->optiexpK;
        int n_where_in_stack_opt=l->l->where_in_stack_opt.size();
        R** where_in_stack =0;
        if (n_where_in_stack_opt && useopt)
            where_in_stack = new R * [n_where_in_stack_opt];
        if (where_in_stack)
        {
            assert(l->l->v.size()==(size_t) n_where_in_stack_opt);
            for (int i=0;i<n_where_in_stack_opt;i++)
            {
                int offset=l->l->where_in_stack_opt[i];
                assert(offset>10);
                where_in_stack[i]= static_cast<R *>(static_cast<void *>((char*)stack+offset));
                *(where_in_stack[i])=0;
            }
            if(poptiexp0) (*poptiexp0)(stack);

            if( (verbosity/100) && verbosity % 10 == 2)
            {
                int il=0;

                for (LinearOperatorD::const_iterator ll=l->l->v.begin();ll!=l->l->v.end();ll++,il++)
                    cout << il << " coef (" << ll->first << ") = " << *(where_in_stack[il]) << " offset=" << l->l->where_in_stack_opt[il] <<endl;

                for (int i=0;i<n_where_in_stack_opt;i++)
                    cout << "const coef " << i << " = " << *(where_in_stack[i]) << endl;
            }
        }
        Stack_Ptr<R*>(stack,ElemMatPtrOffset) =where_in_stack;

        KN<int>   ip(Vh.MaxNbDFPerElement*6);
        if (verbosity >3)
        {
            if (all) cout << " all " << endl ;
            else cout << endl;
        }
        if(di.islevelset() && (kind !=CDomainOfIntegration::int1d)&& (kind !=CDomainOfIntegration::int2d))
            InternalError(" Sorry No levelSet integral for is case ..(5)");


        if (kind==CDomainOfIntegration::int1d)
        {


            if(VF) InternalError(" no jump or average in int1d of RHS");
            if(di.islevelset())
            {
                double uset = HUGE_VAL;
                R2 Q[3];
                KN<double> phi(ThI.nv);phi=uset;
                double f[3];
                for(int t=0; t< ThI.nt;++t)
                {
                    double umx=-HUGE_VAL,umn=HUGE_VAL;
                    for(int i=0;i<3;++i)
                    {
                        int j= ThI(t,i);
                        if( phi[j]==uset)
                        {
                            MeshPointStack(stack)->setP(&ThI,t,i);
                            phi[j]= di.levelset(stack);//zzzz
                        }
                        f[i]=phi[j];
                        umx = std::max(umx,phi[j]);
                        umn = std::min(umn,phi[j]);

                    }
                    if( umn <=0 && umx >= 0)
                    {

                        int np= IsoLineK(f,Q,1e-10);
                        if(np==2)
                        {
                            if ( sameMesh )
                            {/*
                              void  Element_rhs(const FElement & Kv,const LOperaD &Op,double * p,void * stack,KN_<R> & B,
                              const QuadratureFormular1d & FI ,const R2 & PA,const R2 &PB)

                              */
                                Element_rhs<R>(Vh[t],*l->l,buf,stack,*B,FIE,Q[0],Q[1],useopt);
                            }
                            else if(!mapt)
                            Element_rhs<R>(ThI,ThI[t],Vh,0,ThI[t].lab,*l->l,buf,stack,*B,FIE,false,intmortar,Q,useopt);
                            else
                            ffassert(0); //Element_rhs<R>(mapt,ThI,ThI[t],Vh,0,ThI[t].lab,*l->l,buf,stack,*B,FIE,false,intmortar,Q,useopt);

                            //InternalError(" No levelSet on Diff mesh :    to day  int1d of RHS");
                        }
                        if(sptrclean) sptrclean=sptr->clean();
                    }
                }

            }
            else
                for( int e=0;e<ThI.nbe;e++)
                {
                    if (all || setoflab.find(ThI.be(e).lab) != setoflab.end())
                    {
                        int ie,i =ThI.BoundaryElement(e,ie);
                        if ( sameMesh )
                        Element_rhs<R>(Vh[i],ie,Th.be(e).lab,*l->l,buf,stack,*B,FIE,false,useopt);
                        else if(!mapt)
                        Element_rhs<R>(ThI,ThI[i],Vh,ie,Th.be(e).lab,*l->l,buf,stack,*B,FIE,false,intmortar,0,useopt);
                        else
                            ffassert(0);//Element_rhs<R>(mapt,ThI,ThI[i],Vh,ie,Th.be(e).lab,*l->l,buf,stack,*B,FIE,false,intmortar,0,useopt);
                        if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                    }
                }
        }
        else if (kind==CDomainOfIntegration::intalledges)
        {
            ffassert(mapt==0);
            cerr << " intalledges on fespaceS (to do) " << endl;
            InternalError(" intalledges (to do) ");
            ffassert(0);

            /*if(VF)   code 2d .... ici
            {
                pair_stack_double bstack(stack,& binside);

                //bstack.first = stack;
                //bstack.second= & binside;

                //InternalError(" Today no jump or average in intalledges of RHS ");
                for (int i=0;i< ThI.nt; i++)
                    if (all || setoflab.find(ThI[i].lab) != setoflab.end())
                    {

                        for (int ie=0;ie<3;ie++)
                            if ( sameMesh)
                            {
                                int iie=ie,ii=Th.ElementAdj(i,iie);
                                if(ii<0) ii=i;//  sur le bord
                                const TriangleS & K(ThI[i]);
                                int e0=VerticesOfTriangularEdge[ie][0];
                                int e1=VerticesOfTriangularEdge[ie][1];
                                int i1 = ThI(K[e0]),i2 = ThI(K[e1]);
                                BoundaryEdge * be = ThI.TheBoundaryEdge(i1,i2);
                                int lab = be ? be->lab :  notalabel;

                                Element_rhsVF<R>(Vh[i],Vh[ii],ie,iie,lab,*l->l,buf,ip,&bstack,*B,FIE,useopt);
                            }
                            else
                                InternalError("To Do") ;
                        if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                }

            }
            else
                for (int i=0;i< ThI.nt; i++)
                    if (all || setoflab.find(ThI[i].lab) != setoflab.end())
                    {
                        for (int ie=0;ie<3;ie++)
                            if ( sameMesh)
                            Element_rhs<R>(Vh[i],ie,Th[i].lab,*l->l,buf,stack,*B,FIE,true,useopt);
                            else
                                InternalError("To Do") ;
                        if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                    }*/
        }
        else if (kind==CDomainOfIntegration::intallVFedges)
        {
            cerr << " intallVFedges a faire" << endl;

            InternalError(" intallVFedges a faire ");

            ffassert(0);
            for (int i=0;i< ThI.nt; i++)
            {
                if (all || setoflab.find(ThI[i].lab) != setoflab.end())
                    for (int ie=0;ie<3;ie++)
                    {
                        if ( sameMesh)
                        Element_rhs<R>(Vh[i],ie,Th[i].lab,*l->l,buf,stack,*B,FIE,true,useopt);
                        else
                            InternalError("To Do") ;
                    }
                if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
            }
        }

        else if (kind==CDomainOfIntegration::int2d){
            if(di.islevelset())
            {
                QuadratureFormular FITM(FIT);
                double uset = HUGE_VAL;
                R2 Q[4];
                KN<double> phi(Th.nv);phi=uset;
                double f[3];
                for(int t=0; t< Th.nt;++t)
                {
                    if ( all || setoflab.find(ThI[t].lab) != setoflab.end())
                    {
                        double umx=-HUGE_VAL,umn=HUGE_VAL;
                        for(int i=0;i<3;++i)
                        {
                            int j= ThI(t,i);
                            if( phi[j]==uset)
                            {
                                MeshPointStack(stack)->setP(&ThI,t,i);
                                phi[j]= di.levelset(stack);//zzzz
                            }
                            f[i]=phi[j];
                            umx = std::max(umx,phi[j]);
                            umn = std::min(umn,phi[j]);

                        }
                        if( umx <=0 )
                        {Element_rhs<R>(Vh[t],*l->l,buf,stack,*B,FIT,useopt);}
                        else if( umn <0 )
                        { // coupe ..
                            int i0 = 0, i1 = 1, i2 =2;

                            if( f[i0] > f[i1] ) swap(i0,i1) ;
                            if( f[i0] > f[i2] ) swap(i0,i2) ;
                            if( f[i1] > f[i2] ) swap(i1,i2) ;

                            double c = (f[i2]-f[i1])/(f[i2]-f[i0]); // coef Up Traing
                            if( f[i1] < 0 ) {double y=f[i2]/(f[i2]-f[i1]); c *=y*y; }
                            else {double y=f[i0]/(f[i0]-f[i1]) ; c = 1.- (1.-c)*y*y; };
                            assert( c > 0 && c < 1);
                            double arean = (1-c)*Th[t].mesure();
                            FITM=FIT;
                            FITM*=1-c;
                            ffassert(mapt==0);
                            Element_rhs<R>(Vh[t],*l->l,buf,stack,*B,FITM,useopt);
                        }
                        if(sptrclean) sptrclean=sptr->clean();
                    }
                }
            }
            else {
                for (int i=0;i< ThI.nt; i++)
                    if (all || setoflab.find(ThI[i].lab) != setoflab.end())
                    {
                        if ( sameMesh )
                        Element_rhs<R>(Vh[i],*l->l,buf,stack,*B,FIT,useopt);
                        else
                        Element_rhs<R>(ThI,ThI[i],Vh,*l->l,buf,stack,*B,FIT,useopt);

                        if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                    } }
        }

        if (n_where_in_stack_opt) delete [] where_in_stack;

    }

    // creating an instance of AssembleLinearForm
    // case 3d curve
    template<class R>
    void AssembleLinearForm(Stack stack,const MeshL & Th,const FESpaceL & Vh,KN_<R> * B,const  FormLinear * l )
    {
        typedef typename MeshL::Element Element;
        typedef typename MeshL::BorderElement BorderElement;

        StackOfPtr2Free * sptr = WhereStackOfPtr2Free(stack);
        bool sptrclean=true;
        //     sptr->clean(); // modif FH mars 2006  clean Ptr
        Check(l->l,Vh.N);
        if ( B && B->N() != Vh.NbOfDF) ExecError("AssembleLinearForm size rhs and nb of DF of Vh");
        // if ( & Th != &Vh.Th ) ExecError("AssembleLinearForm on different meshes  ( not implemented FH).");
        KN<double> buf(Vh.MaximalNbOfDF()*last_operatortype*Vh.N*2);
        const CDomainOfIntegration & di= *l->di;
        ffassert(di.d==3);

        const MeshL & ThI = Th;
        bool sameMesh = &ThI == &Vh.Th;

        const bool intmortar=di.intmortar(stack);

        SHOWVERB(cout << " FormLinear " << endl);
        // const vector<Expression>  & what(di.what);

        CDomainOfIntegration::typeofkind  kind = di.kind;
        const GQuadratureFormular<R1> & FIT = di.FIE(stack);
        //const QuadratureFormular & FIT = di.FIT(stack);
        const int useopt=di.UseOpt(stack);
        double binside=di.binside(stack);  // truc FH pour fluide de grad2 (decentrage bizard)
        //  cout << "AssembleLinearForm " << l->l->v.size() << endl;
        set<int> setoflab;
        bool all=true;
        bool VF=l->VF();  // finite Volume or discontinuous Galerkin

        if (verbosity>2) cout << "  -- AssembleLinearForm L, discontinuous Galerkin  =" << VF << " binside = "<< binside
            << " levelset integration " <<di.islevelset()<< " withmap: "<<  di.withmap() << " kind int: " << kind<< "\n";
        //  if( di.withmap()) { ExecError(" no map  in the case (6)??");}
        Expression  const * const mapt=di.mapt[0] ? di.mapt:0;
        sameMesh = sameMesh && !mapt; //

        if(di.islevelset() && CDomainOfIntegration::int1d!=kind   )
            InternalError("So no levelset integration type on no int1d/int2d case (4)");
        Expandsetoflab(stack,di, setoflab,all);
        if (verbosity>3) cout << " Optimized = "<< useopt << ", ";
        const E_F0 * poptiexp0=l->l->optiexp0;
        // const E_F0 & optiexpK=*l->l->optiexpK;
        int n_where_in_stack_opt=l->l->where_in_stack_opt.size();
        R** where_in_stack =0;
        if (n_where_in_stack_opt && useopt)
            where_in_stack = new R * [n_where_in_stack_opt];
        if (where_in_stack)
        {
            assert(l->l->v.size()==(size_t) n_where_in_stack_opt);
            for (int i=0;i<n_where_in_stack_opt;i++)
            {
                int offset=l->l->where_in_stack_opt[i];
                assert(offset>10);
                where_in_stack[i]= static_cast<R *>(static_cast<void *>((char*)stack+offset));
                *(where_in_stack[i])=0;
            }
            if(poptiexp0) (*poptiexp0)(stack);

            if( (verbosity/100) && verbosity % 10 == 2)
            {
                int il=0;

                for (LinearOperatorD::const_iterator ll=l->l->v.begin();ll!=l->l->v.end();ll++,il++)
                    cout << il << " coef (" << ll->first << ") = " << *(where_in_stack[il]) << " offset=" << l->l->where_in_stack_opt[il] <<endl;

                for (int i=0;i<n_where_in_stack_opt;i++)
                    cout << "const coef " << i << " = " << *(where_in_stack[i]) << endl;
            }
        }
        Stack_Ptr<R*>(stack,ElemMatPtrOffset) =where_in_stack;

        KN<int>   ip(Vh.MaxNbDFPerElement*6);
        if (verbosity >3)
        {
            if (all) cout << " all " << endl ;
            else cout << endl;
        }
        if(di.islevelset() && (kind !=CDomainOfIntegration::int1d))
            InternalError(" Sorry No levelSet integral for is case ..(5)");

        if (kind==CDomainOfIntegration::int1d){

            if(VF) InternalError(" no jump or average in int1d of RHS");
            if(di.islevelset())
            {cout << " ok di.islevelset() " << di.islevelset() << endl;}
    
            for (int i=0;i< ThI.nt; i++)
                if (all || setoflab.find(ThI[i].lab) != setoflab.end())
                {
                    if ( sameMesh )
                        Element_rhs<R>(Vh[i],*l->l,buf,stack,*B,FIT,useopt);
                    else
                        Element_rhs<R>(ThI,ThI[i],Vh,*l->l,buf,stack,*B,FIT,useopt);

                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                }


        }
        else if(kind==CDomainOfIntegration::int0d){ // to

            for( int e=0;e<ThI.nbe;e++)
            {
                if (all || setoflab.find(ThI.be(e).lab) != setoflab.end())
                {
                    int ie,i =ThI.BoundaryElement(e,ie);
                    if ( sameMesh )
                      Element_rhs<R>(Vh[i],ie,Th.be(e).lab,*l->l,buf,stack,*B,FIT,false,useopt);
                    else if(!mapt)
                      Element_rhs<R>(ThI,ThI[i],Vh,ie,Th.be(e).lab,*l->l,buf,stack,*B,FIT,false,intmortar,0,useopt);
                    else
                        ffassert(0);//Element_rhs<R>(mapt,ThI,ThI[i],Vh,ie,Th.be(e).lab,*l->l,buf,stack,*B,FIE,false,intmortar,0,useopt);
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                }
            }

        }
        else if(kind==CDomainOfIntegration::intall0d){ // to
            // wrong to day ...
            
            ffassert(mapt==0);
            if(VF)
            { // Add juin 2021 ...
                pair_stack_double bstack(stack,& binside);

                 for (int i=0;i< ThI.nt; i++)
                if (all || setoflab.find(ThI[i].lab) != setoflab.end())
                {

                    for (int ie=0;ie<2;ie++)
                    if ( sameMesh)
                    {
                        int iie=ie,ii=Th.ElementAdj(i,iie);
                        if(ii<0) ii=i;//  sur le bord
                        const Element & K(ThI[i]);
                        BorderElement * be=0; // to def
                        int lab = be ? be->lab :  notalabel;

                        Element_rhsVF<R>(Vh[i],Vh[ii],ie,iie,lab,*l->l,buf,ip,&bstack,*B,useopt);
                    }
                    else
                    InternalError("To Do") ;
                    if(sptrclean) sptrclean=sptr->clean();
                }

            }
            else
                
            for( int i=0;i<ThI.nt;i++)
            {
                if (all || setoflab.find(ThI[i].lab) != setoflab.end())
                {
                    for(int ie=0; ie<2; ++ie)
                    {
                  
                    if ( sameMesh )
                      Element_rhs<R>(Vh[i],ie,Th[i].lab,*l->l,buf,stack,*B,FIT,false,useopt);
                    else if(!mapt)
                      Element_rhs<R>(ThI,ThI[i],Vh,ie,Th[i].lab,*l->l,buf,stack,*B,FIT,false,intmortar,0,useopt);
                    else
                        ffassert(0);//Element_rhs<R>(mapt,ThI,ThI[i],Vh,ie,Th.be(e).lab,*l->l,buf,stack,*B,FIE,false,intmortar,0,useopt);
                    if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
                    }
                }
            }

        }

        else
        {
            cerr << " Error integration on MeshL unknown  kind " << kind << endl;
            ffassert(0);
            
        }
        if (n_where_in_stack_opt) delete [] where_in_stack;

    }


   // creating an instance of AssembleLinearForm
   // 3D curve / 2D on meshL
   template<class R>
   void AssembleLinearForm(Stack stack,const MeshL & Th,const FESpace & Vh,KN_<R> * B,const  FormLinear * l )
   {
       ffassert(0);
   }
   // creating an instance of AssembleLinearForm
   // 3D Surf / 3D volume on meshS
   template<class R>
   void AssembleLinearForm(Stack stack,const MeshS & Th,const FESpace3 & Vh,KN_<R> * B,const  FormLinear * l )
   {
       ffassert(0);
   }
   // 3D curve / 3D Surf on meshL
   template<class R>
   void AssembleLinearForm(Stack stack,const MeshL & Th,const FESpaceS & Vh,KN_<R> * B,const  FormLinear * l )
   {
       ffassert(0);
   }

}// END of NameSpace Fem2D


bool isVF(const list<C_F0> & largs)  // true => VF type of Matrix
{
    list<C_F0>::const_iterator ii,ib=largs.begin(),
    ie=largs.end();

    bool VVF =false;
    int kk=0,err=0;
    
    for (ii=ib;ii != ie;ii++)
    {
        kk++;
        Expression e=ii->LeftValue();
        aType r = ii->left();
        if (r==atype<const  FormBilinear *>())
        {
            const  FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
            bool vvf  = bb->VF();
            if( vvf && (bb->di->kind != CDomainOfIntegration::intalledges)
                    && (bb->di->kind != CDomainOfIntegration::intallVFedges)
                    && (bb->di->kind != CDomainOfIntegration::intallfaces)
                    && (bb->di->kind != CDomainOfIntegration::intall0d)) // Add 15 juin 2021 FH?
            {
                if(err==0) cerr << "\n\n"; 
                cerr << " ** Fatal error in term "<< kk << " of the varf form (integral, on , ... ) " << endl;
                err++;
            }
            VVF = vvf || VVF;
        }
    }
    if(err)
    {
        cerr << " ** number  " << err << " of  error the varf form, with " << kk << " terms "<< endl;
        CompileError("Sorry, no  jump, mean, otherside in bilinear term must be in integral of type  intalledges,  intallVFedges or intallfaces");
    }
    return VVF;
}


bool isSameMesh(const list<C_F0> & largs,const void * Thu,const void * Thv,Stack stack)  // true => VF type of Matrix
{
    if( Thv != Thu ) return false;
    list<C_F0>::const_iterator ii,ib=largs.begin(),
    ie=largs.end();

    // bool VVF =false;
    for (ii=ib;ii != ie;ii++)
    {
        Expression e=ii->LeftValue();
        aType r = ii->left();
        if (r==atype<const  FormBilinear *>())
        {
            const  FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
            const void *  Thbf = GetAny<const void *>((*bb->di->Th)(stack));
            if (Thbf != Thu) return false;
        }
        else if (r==atype<const  FormLinear *>())
        {
            const  FormLinear * bb=dynamic_cast<const  FormLinear *>(e);
            const void * Thbf = GetAny<const void *>((*bb->di->Th)(stack));
            if (Thbf != Thu) return false;
        }
    }
    return true;
}

template<class R,class FESpace,class v_fes>
void InitProblem( int Nb, const FESpace & Uh,
                 const FESpace & Vh,
                 KN<R> *&B,KN<R> *&X,vector<  pair< FEbase<R,v_fes> * ,int> > &u_hh,
                 Data_Sparse_Solver    *ds ,//    *typemat ,
                 vector<  FEbase<R,v_fes> *  > & u_h,const FESpace ** LL, bool initx )
{
    typedef typename  FESpace::Mesh Mesh;
    typedef typename  FESpace::FElement FElement;
    typedef typename  Mesh::Element Element;
    typedef typename  Mesh::Vertex Vertex;
    typedef typename  Mesh::RdHat RdHat;
    typedef typename  Mesh::Rd Rd;

    *B=R();

    //  bool initx = typemat->t==TypeSolveMat::GC;

    const  Mesh & Th(Uh.Th);

    if (initx)
    {
        if (!X || (X =B) )
        X=new KN<R>(B->N());
        const FEbase<R,v_fes> & u_h0 = *(u_h[0]);
        const FESpace  * u_Vh = u_h0.Vh ;

        if ( u_Vh==0  || &((u_h[0])->Vh->Th) != &Th )
        {
            *X=R();
            if(verbosity>1)
            cout << "   -- Change of Mesh " << (u_Vh ? & (*(u_h[0])).Vh->Th: 0 )
            << "  " << &Th <<  endl;
        }
        else
        { //  copy the previous soluton to initialize CG, GMRES, etc ...
            if (Nb==1)
            {  // modif  FH 0701/2005 + april 2006
                if(u_h[0]->x() && u_h[0]->x()->N() != X->N() )// correction juin FH 2021 .. 
                cout << " bug ???? " << endl;
                if (u_h[0]->x() && u_h[0]->x()->N() == X->N() )
                *X= * u_h[0]->x();
                else
                *X=R();
            }
            else { // dispatch the solution
                const FElement ** sK= new const FElement * [Nb];
                KN<R> ** sol= new KN<R> * [Nb];
                for (int i=0;i<Nb;i++) {

                    sol[i] = (*(u_h[i])).x() ;
                }

                for (int it=0;it<Th.nt;it++)
                {
                    const FElement K(Uh[it]);
                    const int nbdf=K.NbDoF();
                    for (int i=0;i<Nb;i++)
                    sK[i]= new FElement( (*LL[i])[it]) ;
                    for (int df=0;df< nbdf;df++)
                    {  int kfe=K.FromFE(df);
                        int kdf=K.FromDF(df);
                        if (sol[kfe]) {
                            const FElement & SK(*sK[kfe]);
                            (*X)[K(df)]= (*sol[kfe])[SK(kdf)] ;
                        }
                        else (*X)[K(df)]= R();
                    }
                    for (int i=0;i<Nb;i++)
                    delete sK[i];
                }
                delete [] sol;
                delete [] sK;
            }}
    }


}



template<class R>
 MatriceCreuse<typename CadnaType<R>::Scalaire> * DefSolverCadna(
  Stack stack,
  MatriceCreuse<R>  & A,
  Data_Sparse_Solver & ds
/*  long NbSpace ,
  long itmax,
  double & eps,
  bool initmat,
  int strategy,
  const OneOperator *precon,
  double tgv,
  double tol_pivot, double tol_pivot_sym
*/
)
{
   typedef typename CadnaType<R>::Scalaire R_st;
    /*
 //  MatriceCreuse<R_st> *CadnaMat;
    if (ds.typemat->profile)
      {
        if(verbosity>5) cout << " Matrix skyline type:" << ds.typemat->t <<endl;
        MatriceProfile<R> & AAA(dynamic_cast<MatriceProfile<R> &>(A));
        MatriceProfile<R_st> &AA(*new MatriceProfile<R_st>(AAA)); //

        throwassert(&AA);
        double tol_pivot1= (ds.tol_pivot>0) ? ds.tol_pivot : EPSILON/8.;
       // cout << " tol_pivot1 " <<tol_pivot1 <<  endl;
        switch (ds.typemat->t) {
        case TypeSolveMat::LU       : AA.LU(tol_pivot1); break;
        case TypeSolveMat::CROUT    : AA.crout(tol_pivot1); break;
        case TypeSolveMat::CHOLESKY : AA.cholesky(tol_pivot1); break;
        default:
          cerr << " type resolution " << ds.typemat->t << endl;
          CompileError("type resolution profile inconnue"); break;
        }
        return &AA;
      }
    else */
      {
         ExecError("matrix HMAT & CADNA are incompatible today, sorry!");
         return 0;
      }
   return 0;
  }

template<class R,class FESpace,class v_fes>
void   DispatchSolution(const typename FESpace::Mesh & Th,int Nb, vector<  FEbase<R,v_fes> * > & u_h,KN<R> * X,KN<R> * B,const FESpace **  LL,const FESpace &  Uh)
{
    typedef typename  FESpace::Mesh Mesh;
    typedef typename  FESpace::FElement FElement;
    typedef typename  Mesh::Element Element;
    typedef typename  Mesh::Vertex Vertex;
    typedef typename  Mesh::RdHat RdHat;
    typedef typename  Mesh::Rd Rd;

    // dispatch the solution
    if (Nb==1)  {
        *(u_h[0])=X;
        if (X != B ) delete B;  }
    else {
        const FElement ** sK= new const FElement * [Nb];

        KN<R> ** sol= new KN<R> * [Nb];
        for (int i=0;i<Nb;i++) {
            sol[i]= new KN<R>( LL[i]->NbOfDF) ;
            *(u_h[i]) = sol[i];
        }

        for (int it=0;it<Th.nt;it++)
        {
            const FElement K(Uh[it]);
            const int nbdf=K.NbDoF();
            for (int i=0;i<Nb;i++)
            sK[i]= new FElement( (*LL[i])[it]) ;
            for (int df=0;df< nbdf;df++)
            {  int kfe=K.FromFE(df);
                int kdf=K.FromDF(df);
                const FElement & SK(*sK[kfe]);
                (*sol[kfe])[SK(kdf)] = (*X)[K(df)];
            }
            for (int i=0;i<Nb;i++)
            delete sK[i];

        }

        delete [] sK;
        delete [] sol;
        if (X != B && X ) delete X;
        delete B;
    }
}
/*
#ifdef HAVE_LIBUMFPACK
TypeSolveMat::TSolveMat  TypeSolveMat::defaultvalue=TypeSolveMat::SparseSolver;
#else
TypeSolveMat::TSolveMat  TypeSolveMat::defaultvalue=TypeSolveMat::LU;
#endif
*/

template<class R,class FESpace,class v_fes>    // TODO if coupling FE wit problem
AnyType Problem::eval(Stack stack,Data<FESpace> * data,CountPointer<MatriceCreuse<R> > & dataA,
                      MatriceCreuse< typename CadnaType<R>::Scalaire >   * & cadnamat ) const
{
    typedef typename  FESpace::Mesh MeshT;
    typedef typename  FESpace::FElement FElement;
    typedef typename  MeshT::Element Element;
    typedef typename  MeshT::Vertex Vertex;
    typedef typename  MeshT::RdHat RdHat;
    typedef typename  MeshT::Rd Rd;

    using namespace Fem2D;
    typedef typename CadnaType<R>::Scalaire R_st;
    MeshPoint *mps= MeshPointStack(stack),mp=*mps;
    Data_Sparse_Solver ds;
    /* long NbSpace = 50;
     long itmax=0;
     double epsilon=1e-6;*/
    string save;

     KN<double>* cadna=0;

    if (nargs[0]) save = *GetAny<string*>((*nargs[0])(stack));
    if (nargs[1]) cadna= GetAny<KN<double>* >((*nargs[1])(stack));

    SetEnd_Data_Sparse_Solver<R>(stack,ds,nargs,n_name_param);


    //  for the gestion of the PTR.
    WhereStackOfPtr2Free(stack)=new StackOfPtr2Free(stack);// FH aout 2007

    bool sym = ds.sym;

    list<C_F0>::const_iterator ii,ib=op->largs.begin(),
    ie=op->largs.end();
    int Nbcomp2=var.size(),Nbcomp=Nbcomp2/2; // nb de composante
    throwassert(Nbcomp2==2*Nbcomp);
    //  Data *data= dataptr(stack);
    //   data->init();
    KN<int>  which_comp(Nbcomp2),which_uh(Nbcomp2);

    TabFuncArg tabexp(stack,Nbcomp);
    typedef pair< FEbase<R,v_fes> *,int> pfer;
    vector< pair< FEbase<R,v_fes> *,int> > u_hh(Nbcomp2);
    // u_hh.first --> FEbase *
    // u_hh.second --> numero de la composante dans le cas d'un vectorial FESpace (ex: pfes*_tefk)

    for (size_t i=0;i<var.size();i++)
    u_hh[i] = GetAny< pfer  >( (*(var[i]))(stack));
    for (size_t i=0;i<var.size();i++)
    u_hh[i].first->newVh();
    //   compression pour les cas vectoriel
    int kkk=0;
    for (int i=0;i<Nbcomp2;i++)
    {
        if ( u_hh[i].second==0)
        kkk++;
        else {
            throwassert(u_hh[i].second==(u_hh[i-1].second+1));} // verificationd que les composantes des vectorial FESpace sont dans le bon ordre
        which_uh[i]=kkk-1;               // set the numero of the FESpace
        which_comp[i]=u_hh[i].second;    // set the component of the FESpace (=0 for scalar FESpace by definition)
    }

    vector<  FEbase<R,v_fes> * > u_h(kkk); // the list of the true pointer to FEbase
    kkk= 0;
    for (int i=0;i<Nbcomp2;i++)
    if ( u_hh[i].second==0) u_h[kkk++]=u_hh[i].first;
    const int  Nb2 = kkk, Nb=Nb2/2; // nb of FESpace
    throwassert(Nb2==2*Nb);

    //const FESpace ** LL = new  const FESpace *[var.size()];
    KN<const FESpace *> LL(var.size());    // creation de la liste des FESpace" Nb2 <= var.size() "
    for (int i=0;i<Nb2;i++)
    LL[i]= (*(u_h[i])).newVh();
    SHOWVERB(cout << "Problem  " << Nb << endl);

    //   const de

    //  const FESpace * Uhh , *Vhh;
    // check we have same Th for each FESpace
    const MeshT * pTh= &LL[0]->Th;
    for (int i=0;i<Nb2;i++)
    if ( &LL[i]->Th != pTh)
    ExecError("all the finites elements spaces must be defined on the same mesh in solve");

    if ( pTh != data->pTh )
    {
        // set value of the struct data 
        ds.initmat = true;
        data->pTh=pTh;
        if (Nb==1)
        { //  cas scalaire
            data->Uh=LL[0];
            data->Vh=LL[1]; }
        else
        { //  cas vectoriel 
            // check that we have the same finconue space and ftest space
            bool same=true;
            for (int i=0;i<Nb;i++)
            if ( LL[i] != LL[Nb+i] )
            {
                same = false;
                break;
            }
            if(!same)
            InternalError("Methode de Galerkine (a faire)");
            else
            {
                // check if we have only one FESpace
                bool unique=true;
                for (int i=1;i<Nb;i++)
                if ( LL[0] != LL[i])
                {
                    unique = false;
                    break;
                }
                else if(LL[i]->FirstDfOfNodeData)// Correct jan 2015 not FE product this case
                {
                    unique = false;
                    break;
                }

                if (unique)
                data->Uh.master( new FESpace(*LL[0],Nb));
                else
                data->Uh.master(new FESpace(LL,Nb));
                data->Vh=data->Uh;
            }

        }
    }

    const FESpace & Uh(*data->Uh);
    const FESpace & Vh(*data->Vh);
    throwassert(Nbcomp==Uh.N && Nbcomp==Vh.N);
    KN<R> *B=new KN<R>(Vh.NbOfDF);
    KN<R> *X=B; //
    const  MeshT & Th(Uh.Th);
    bool initx = true; //typemat->t==TypeSolveMat::GC ; //  make x and b different in all case
    // more safe for the future ( 4 days lose with is optimization FH )

    InitProblem<R,FESpace,v_fes>(  Nb,  Uh, Vh, B, X,u_hh,&ds , u_h,  LL,  initx);

    if(verbosity>2) cout << "   Problem(): initmat " << ds.initmat << " VF (discontinuous Galerkin) = " << VF << endl;



    if (ds.initmat)
     {
       {
          if ( &Uh == & Vh )
            dataA.master(new MatriceMorse<R>( Vh.NbOfDF,Vh.NbOfDF,sym));
          else
            dataA.master(new MatriceMorse<R>(Vh.NbOfDF,Uh.NbOfDF,false));
        }
        MatriceCreuse<R>  & AA(dataA);
       if(verbosity>1) cout <<  "   -- size of Matrix " << AA.size()<< " Bytes" <</* " skyline =" <<ds.typemat->profile <<*/ endl;
      }
    MatriceCreuse<R>  & A(dataA);
    if  (AssembleVarForm( stack,Th,Uh,Vh,sym, ds.initmat ? &A:0 , B, op->largs))
    {
        *B = - *B;
        // hach FH
        for (int i=0, n= B->N(); i< n; i++)
        if( abs((*B)[i]) < 1.e-60 ) (*B)[i]=0;

        AssembleBC<R,MeshT,FESpace,FESpace>     ( stack,Th,Uh,Vh,sym, ds.initmat ? &A:0 , B, initx ? X:0,  op->largs, ds.tgv );   // TODO with problem
    }
    else
    *B = - *B;

    dynamic_cast<HashMatrix<int,R> *>(&A)->half = ds.sym;

    MatriceCreuse<R_st>  * ACadna = 0;


    try {

        if (ds.initmat)
        {
         //   if(cadna)
         //   ACadna = DefSolverCadna( stack,A, ds);
         //   else
            DefSolver(stack,  A, ds);
        }



        // if(verbosity>3) cout << "   B  min " << B->min() << " ,  max = " << B->max() << endl;
        if( save.length() )
        {
            string savem=save+".matrix";
            string saveb=save+".b";
            {
                ofstream outmtx( savem.c_str());
                A.dump(outmtx)  << endl;
            }
            {
                ofstream outb(saveb.c_str());
                outb<< *B << endl;
            }

        }
        if (verbosity>99)
        {
            cout << " X= " << *X << endl;
            cout << " B= " << *B << endl;
        }

        if(ACadna)
        {
            KN<R_st> XX(*X);
            KN<R_st> BB(*B);
            ACadna->Solve(XX,BB);
            *X=XX;
            *cadna =-1.;

#ifdef HAVE_CADNA
            R_st xxmin = XX.min();
            R_st xxmax = XX.max();
            cout  << "    cadna:      min " <<  xxmin << "/ nd " << cestac(xxmin)
            << " ,   max " << xxmax << " / nd " << cestac(xxmax)   << endl ;
            int nn= XX.N();
            if ( cadna->N() == nn )
            for (int i=0;i<nn;++i)
            (*cadna)[i] = cestac(XX[i]);
            else
            cerr << "Warning: Sorry array is incorrect size to store cestac "
            << nn << " != " << cadna->N() << endl;
#endif
        }
        else

        A.Solve(*X,*B);

        if (verbosity>99)
        {
            cout << " X= " << *X << endl;
        }
    }
    catch (...)
    {
        if(verbosity) cout << " catch an erreur in  solve  =>  set  sol = 0 !!!!!!! "   <<  endl;
        *X=R(); // erreur set the sol of zero ????
        DispatchSolution<R,FESpace,v_fes>(Th,Nb,u_h,X,B,LL,Uh);
        throw ;
    }
    DispatchSolution<R,FESpace,v_fes>(Th,Nb,u_h,X,B,LL,Uh);


    if (verbosity)
    {cout << "  -- Solve : \n" ;
        for (int i=0;i<Nb;i++)
        cout  << "          min " << (u_h[i])->x()->min() << "  max " << (u_h[i])->x()->max() << endl ;
    }

    // delete [] LL;
    // if (save) delete save; // clean memory
    *mps=mp;
    return SetAny<const Problem *>(this);
}


// dimProblem read the number of arguments of problem ex: problem a(u,v) or a([u1,u2], [v1,v2])
int dimProblem(const ListOfId &l)
{
    int dim=0;
    int nb=l.size();//,nbarray=0;//,n=0,
    //const UnId *p1;
    for(int i=0; i<nb; ++i)
    {
        if(l[i].e ==0)// to miss name parameter solver=ddd
        {
        if (l[i].array)
        {
            ListOfId * array=l[i].array;
            for(int j=0; j<array->size(); ++j)
            {
                const UnId & idi( (*array)[j]);
                if (idi.r == 0 && idi.re  == 0 && idi.array==0 )
                {
                    C_F0 c=::Find( idi.id);
                    if(BCastTo<pfec>(c) ) ffassert(dim==0 || dim==2),dim=2;
                    if(BCastTo<pfer>(c) ) ffassert(dim==0 || dim==2),dim=2;
                    if(BCastTo<pf3c>(c) ) ffassert(dim==0 || dim==3),dim=3;
                    if(BCastTo<pf3r>(c) ) ffassert(dim==0 || dim==3),dim=3;
                    if(BCastTo<pfSr>(c) ) ffassert(dim==0 || dim==4),dim=4;
                    if(BCastTo<pfSc>(c) ) ffassert(dim==0 || dim==4),dim=4;
                    if(BCastTo<pfLr>(c) ) ffassert(dim==0 || dim==5),dim=5;
                    if(BCastTo<pfLc>(c) ) ffassert(dim==0 || dim==5),dim=5;
                }
            }

        }
        else
        {
            C_F0 c=::Find(l[i].id);
            if(BCastTo<pfec>(c) ) ffassert(dim==0 || dim==2),dim=2;
            if(BCastTo<pfer>(c) ) ffassert(dim==0 || dim==2),dim=2;
            if(BCastTo<pf3c>(c) ) ffassert(dim==0 || dim==3),dim=3;
            if(BCastTo<pf3r>(c) ) ffassert(dim==0 || dim==3),dim=3;
            if(BCastTo<pfSr>(c) ) ffassert(dim==0 || dim==4),dim=4;
            if(BCastTo<pfSc>(c) ) ffassert(dim==0 || dim==4),dim=4;
            if(BCastTo<pfLr>(c) ) ffassert(dim==0 || dim==5),dim=5;
            if(BCastTo<pfLc>(c) ) ffassert(dim==0 || dim==5),dim=5;
        }
        }
    }
    ffassert(dim);
    return dim;

}

AnyType Problem::operator()(Stack stack) const
{
    if(dim==2) {
        Data<FESpace> *data= dataptr(stack);
        if (complextype)
            return eval<Complex,FESpace,v_fes>(stack,data,data->AC,data->AcadnaC);
        else
            return eval<double,FESpace,v_fes>(stack,data,data->AR,data->AcadnaR);
    }
    else if(dim==3) {
        Data<FESpace3> *data= dataptr3(stack);
        if (complextype)
            return eval<Complex,FESpace3,v_fes3>(stack,data,data->AC,data->AcadnaC);
        else
            return eval<double,FESpace3,v_fes3>(stack,data,data->AR,data->AcadnaR);
    }
    else if(dim==4) {
        Data<FESpaceS> *data= dataptrS(stack);
        if (complextype)
            return eval<Complex,FESpaceS,v_fesS>(stack,data,data->AC,data->AcadnaC);
        else
            return eval<double,FESpaceS,v_fesS>(stack,data,data->AR,data->AcadnaR);
    }
    else if(dim==5) {
        Data<FESpaceL> *data= dataptrL(stack);
        if (complextype)
            return eval<Complex,FESpaceL,v_fesL>(stack,data,data->AC,data->AcadnaC);
        else
            return eval<double,FESpaceL,v_fesL>(stack,data,data->AR,data->AcadnaR);
    }

    else ffassert(0);
}

template<class pfer,class pfec>
bool GetBilinearParam(const ListOfId &l,basicAC_F0::name_and_type *name_param,int n_name_param,
                      Expression *nargs,int & N,int & M,  vector<Expression> & var )
{
    bool unset=true,complextype=false;

    for (int i=0;i<n_name_param;i++)
    nargs[i]=0;
    int nb=l.size(),n=0,nbarray=0;
    ListOfId * array[2];
    for (int i=0;i<nb;i++)
    if (l[i].r == 0 && l[i].re  == 0 && l[i].array == 0)
    n++;
    else if (l[i].array) array[Min(nbarray++,1)] = l[i].array;
    else
    {
        bool ok=false;
        for (int j=0;j<n_name_param;j++)
        if (!strcmp(l[i].id,name_param[j].name))
        {
            ok = !nargs[j];
            nargs[j]= map_type[name_param[j].type->name()]->CastTo(C_F0(l[i].e,l[i].re));
            break;
        }
        if (!ok)
        {
            cerr << " Error name argument " << l[i].id << " the kown arg : ";
            for (int k=0;k<n_name_param;k++)
            cerr << name_param[k].name << " ";
            cerr << endl;
            CompileError("Unknown name argument or two times same name argument ");
        }
    }

    if (nbarray)
    { // new version ok
        if(nbarray!=2)
        CompileError(" Must have 2 array, one for unknown functions, one for test functions");
        N = array[0]->size();
        M = array[1]->size();
        var.resize(N+M);
        for (size_t k=0,j=0;k<2;k++)
        for  (size_t i=0;i<array[k]->size();i++)
        {
            const UnId & idi((*array[k])[i]);
            if (idi.r == 0 && idi.re  == 0 && idi.array==0 )
            { C_F0 c=::Find( idi.id);
                if (unset)
                complextype =  BCastTo<pfec>(c) , unset=false;

                if(complextype)
                var[j++]=CastTo<pfec>(c);
                else
                var[j++]=CastTo<pfer>(c);
            }
            else
            CompileError(" Just Variable in array parameter ");
        }
    }
    else
    { // old version
        assert(n%2==0);
        N=n/2;
        M=N;
        var.resize(N+M);
        for  (size_t i=0,j=0;i<l.size();i++)
        if (l[i].r == 0 && l[i].re  == 0 && l[i].array==0 )
        {
            C_F0 c=::Find(l[i].id);
            if (unset)
            complextype =  BCastTo<pfec>(c) , unset=false;
            if(complextype)
            var[j++]=CastTo<pfec>(c);
            else
            var[j++]=CastTo<pfer>(c);
        }

    }
    return complextype;
}


/*
 int DimForm( list<C_F0> & largs)
 {
 int dim=0;
 list<C_F0>::iterator ii,ib=largs.begin(),
 ie=largs.end();
 for (ii=ib;ii != ie;ii++)
 {
 Expression e=ii->LeftValue();
 aType r = ii->left();
 if (r==atype<const  FormBilinear *>())
 {
 const  FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
 if(dim) ffassert(bb->d==dim);
 else
 dim=bb->d;
 }
 else if (r==atype<const  FormLinear *>())
 {
 const  FormLinear * ll=dynamic_cast<const  FormLinear *>(e);
 if(dim) ffassert(bb->d==dim);
 else
 dim=bb->d;
 }
 else if (r == atype<const  BC_set *>())
 {
 const  BC_set * bc=dynamic_cast<const  BC_set *>(e);
 if (bc->complextype)  complextype=true;
 }
 }
 }*/
bool CheckSizeOfForm( list<C_F0> & largs ,int N,int M)
{
    list<C_F0>::iterator ii,ib=largs.begin(),
    ie=largs.end();
    for (ii=ib;ii != ie;ii++)
    {
        Expression e=ii->LeftValue();
        aType r = ii->left();
        if (r==atype<const  FormBilinear *>())
        {
            const  FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
            const Foperator * b=const_cast<  Foperator *>(bb->b);
            //  verif
         }
        else if (r==atype<const  FormLinear *>())
        {
            const FormLinear * ll=dynamic_cast<const  FormLinear *>(e);
            const Ftest * l= const_cast<Ftest *>(ll->l);
        }
        else if (r==atype<const  BC_set *>())
        {

            const BC_set * bc= dynamic_cast<const  BC_set *>(e);
         }

    }
    return true;
}

bool FieldOfForm( list<C_F0> & largs ,bool complextype)  // true => complex problem
{
    //  bool   iscomplextype=complextype;
    list<C_F0>::iterator ii,ib=largs.begin(),
    ie=largs.end();
    // bool complextype =false;
    for (ii=ib;ii != ie;ii++)
    {
        Expression e=ii->LeftValue();
        aType r = ii->left();
        if (r==atype<const  FormBilinear *>())
        {
            const  FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
            if (! bb->b->mappable(BCastToR))
            complextype=true;
        }
        else if (r==atype<const  FormLinear *>())
        {
            const  FormLinear * ll=dynamic_cast<const  FormLinear *>(e);
            if (! ll->l->mappable(BCastToR))
            complextype=true;
        }
        else if (r == atype<const  BC_set *>())
        {
            const  BC_set * bc=dynamic_cast<const  BC_set *>(e);
            if (bc->complextype)  complextype=true;

        }
    }

    for (ii=ib;ii != ie;ii++)
    {
        Expression e=ii->LeftValue();
        aType r = ii->left();
        if (r==atype<const  FormBilinear *>())
        {
            FormBilinear * bb=new FormBilinear(*dynamic_cast<const FormBilinear *>(e));
            Foperator * b=const_cast<  Foperator *>(bb->b);
            // const Foperator * b=bb->b;
            //cout << b <<  " bb->b " <<  bb->b << " " <<  bb->b <<  " " << bb->b->isoptimize <<endl;
            assert(b->isoptimize==false);
            if (complextype)  b->mapping(&CCastToC);
            else b->mapping(&CCastToR) ;
            Foperator * bn = b->Optimize(currentblock);
            *bb->b = *bn;
            *ii=C_F0(bb,r);
        }
        else if (r==atype<const  FormLinear *>())
        {
            FormLinear * ll=new FormLinear(*dynamic_cast<const  FormLinear *>(e));
            Ftest * l= const_cast<Ftest *>(ll->l);
            if (complextype)  l->mapping(&CCastToC) ;
            else l->mapping(&CCastToR) ;
            Ftest * ln = l->Optimize(currentblock);
            *ll->l=*ln;
            *ii=C_F0(ll,r);
            //cout << l <<   " ll->l " <<  ll->l << " " << ll->l->isoptimize <<endl;
        }
        else if (r==atype<const  BC_set *>())
        {// modif FH  mai 2007  A FAIRE il y a un bug ici XXXXXXXXXXXXX

            BC_set * bc= new BC_set(*dynamic_cast<const  BC_set *>(e));
            if (complextype && !bc->complextype) {
                bc->CastToK<Complex>() ;
                if(verbosity > 10) cout << " Bc to complex " << endl;
            }
            //else bc->mapping(&CCastToR) ;
            //cout << l <<   " ll->l " <<  ll->l << " " << ll->l->isoptimize <<endl;
            *ii=C_F0(bc,r);
        }

    }
    return complextype;
}


Problem::Problem(const C_args * ca,const ListOfId &l,size_t & top) :
op(new C_args(*ca)),
var(l.size()),
VF(false),
offset(align8(top)),
dim(dimProblem(l))
{
    if( verbosity > 999)  cout << "Problem : ----------------------------- " << top << " dim = " << dim<<" " << nargs <<  endl;
    top = offset + max(sizeof(Data<FESpace>),sizeof(Data<FESpace>));

    bool iscomplex;
    if(dim==2)
    iscomplex=GetBilinearParam<pfer,pfec>(l,name_param,n_name_param,nargs, Nitem,Mitem,var);
    else if (dim==3)
    iscomplex=GetBilinearParam<pf3r,pf3c>(l,name_param,n_name_param,nargs, Nitem,Mitem,var);
    else if (dim==4)  // dim = 4 for a 3D surface problem
    iscomplex=GetBilinearParam<pfSr,pfSc>(l,name_param,n_name_param,nargs, Nitem,Mitem,var);
    else if (dim==5)  // dim = 5 for a 3D curve problem
        iscomplex=GetBilinearParam<pfLr,pfLc>(l,name_param,n_name_param,nargs, Nitem,Mitem,var);
    else ffassert(0); // bug

    precon = 0; //  a changer
    if ( nargs[3+3])
    {
        const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(nargs[3+3]);
        assert(op);
        precon = op->Find("(",ArrayOfaType(atype<KN<R>* >(),false));
        ffassert(precon);
    }

    VF=isVF(op->largs);
    // cout << " Problem ) VF = " << VF << endl;
    complextype =  FieldOfForm(op->largs,iscomplex)  ;  // Warning do the casting of all expression in double or complex
    if( complextype && !iscomplex)
    CompileError("Error: Problem  a complex problem with no complex FE function ");
    if( verbosity > 1)
    cout << "  -- Problem type  ( complex : " << complextype << " )  "  <<endl;
}

Expression IsFebaseArray(Expression f)
{
    assert(f);
    size_t N=f->nbitem();
    E_Array * vvi(dynamic_cast< E_Array *>(f));
    if ( ! vvi) return 0;
    E_Array & vi(*vvi);
    Expression febase=0;
    for (size_t i=0;i<N;i++)
    {
        assert(vi[i].left() == atype<pfer>() );
        const E_FEcomp<R,v_fes> * comp=dynamic_cast<const E_FEcomp<R,v_fes> *>( vi[i].LeftValue()) ;
        if (!(comp && comp->comp == (int) i  && comp->N == (int) N)) return 0;
        if (!febase) febase = comp->a0;
        else if(comp->a0 != febase) return 0;
    }
    return febase;
}
template<class VFES1, class VFES2>
Call_FormBilinear< VFES1, VFES2>::Call_FormBilinear(Expression * na,Expression  BB,Expression fi, Expression fj)
: nargs(na),largs(),N(fi->nbitem()),M(fj->nbitem()),
euh(fi), evh(fj)
{
    assert(nargs );
    const C_args * LLL=dynamic_cast<const C_args *>(BB);
    if (!LLL)
    CompileError("Sorry the variationnal form (varf)  is not a the variationnal form (type const C_args *)");
    largs=LLL->largs;
}
template<class VFES>
Call_FormLinear<VFES>::Call_FormLinear(Expression *na,Expression  LL, Expression ft)
:largs(),nargs(na),N(ft->nbitem()),
ppfes(ft)//IsFebaseArray(ft))
{
    const C_args * LLL=dynamic_cast<const C_args *>(LL);
    if ( !LLL) CompileError("The parameter of a LinearForm must be a array of all componate of FE function");
    largs=LLL->largs;
}
bool C_args::IsLinearOperator() const {
    //  int n=largs.size();
    aType tRn =atype<KN<R>* >();
    aType tCn =atype<KN<Complex>* >();
    for (const_iterator i=largs.begin(); i != largs.end();i++)
    {
        C_F0  c= *i;
        // Expression e=c;
        aType r=c.left();
        if (     ( r != atype<const  FormLinear *>() )
 	       &&  ( r != atype<const  BC_set *>() )
 	       &&  ( r != atype<RNM_VirtualMatrix<R>::plusAx >() )
 	       &&  ( r != atype<RNM_VirtualMatrix<R>::plusAtx >() )
 	       &&  ( r != atype<RNM_VirtualMatrix<Complex>::plusAx >() )
 	       &&  ( r != atype<RNM_VirtualMatrix<Complex>::plusAtx >() )
 	       &&  ( r != tRn)
 	       &&  ( r != tCn)
            ) return false;
    }
    return true;}

bool C_args::IsBilinearOperator() const {
    //int n=largs.size();
    aType tRn =atype<Matrice_Creuse<R>* >();
    aType tCn =atype<Matrice_Creuse<Complex>* >();
    for (const_iterator i=largs.begin(); i != largs.end();i++)
    {
        C_F0  c= *i;
        //Expression e=c;
        aType r=c.left();
        if (     ( r!= atype<const  FormBilinear *>() )
            &&  ( r != atype<const  BC_set *>() )
            &&  ( r != tRn)
            &&  ( r != tCn)
            ) return false;
    }
    return true;}


void SetArgsFormLinear(const ListOfId *lid,int ordre)
{
    //  the local parameter are
    //  ordre ==2 => bilinear form  unknown (newU_) and test function (newV_)
    //  ordre ==1 =>   linear form just  test function (newV_)
    // ---------------------
    throwassert(ordre >0 && ordre <=2 && (lid || lid->size()>0 ) );
    const ListOfId & l(*lid);
    int nb=l.size();
    int n=0;
    C_F0 type,init;
    int nbarray=0;
    ListOfId * array[2];
    aType uh=atype<const finconnue*>(),vh=atype<const ftest*>();

    for (int i=0;i<nb;i++)
    if (l[i].r == 0 &&  l[i].re == 0 && l[i].id  ) n++;
    else if (l[i].array)
    array[Min(nbarray++,2)] = l[i].array;
    if (nbarray && n==0)
    {  //

        if(nbarray!=ordre)
        { cerr << " form " << ordre << " == " << nbarray << " Nb of Array "<<endl;
            CompileError(" Must have 1 or 2 array, one for unknown functions, one for test functions");
        }
        for (int k=0;k<ordre;k++)
        for  (int i=0,iend=array[k]->size();i<iend;i++)
        {
            const UnId & idi((*array[k])[i].id);
            if (idi.r == 0 && idi.re  == 0 && idi.array==0 )
            {
                if (k==ordre-2)  //  unknow function just in case of bilinear form
                currentblock->NewID(uh,idi.id,C_F0(newU_(i),uh));
                else   //  test function
                currentblock->NewID(vh,idi.id,C_F0(newV_(i),vh));
            }
            else
            CompileError(" Just Variable in array parameter ");
        }
    }
    else if (nbarray==0)
    {    // a supprimer  to remove   in case of bilinear

        SHOWVERB(cout << "SetArgs:: form  set parameter " << endl);
        if( ! ( ordre==1 || n%2==0) )
        CompileError(" Error in test or unknown function (odd number of function) ");
        ffassert( ordre==1 || n%2==0);
        int nn=ordre==1 ? 0 : n/2; // order == 1 => no unknown function just test function

        for (int i=0,j=0;i<nb;i++)
        if (l[i].r == 0 && l[i].re  == 0 && l[i].array==0)
        {
            SHOWVERB(cout <<"  " <<  l[i].id  << " " << (j<nn) << endl);
            if (j<nn)
            currentblock->NewID(uh,l[i].id,C_F0(newU_(j%nn),uh));
            else
            currentblock->NewID(vh,l[i].id,C_F0(newV_(j%nn),vh));
            j++;
        }
    }
    else
    {
        CompileError(" Sorry you mixte formulation with and without array ");
    }
}

const Fem2D::GQuadratureFormular<R3> & CDomainOfIntegration::FIV(Stack stack) const
{
    using namespace Fem2D;
    if (nargs[8]) return  *GetAny<const Fem2D::GQuadratureFormular<R3> *>((*nargs[8])(stack));
    int exact = 5;
    if (nargs[2]) exact=  GetAny<long>((*nargs[2])(stack))-1;
    GQuadratureFormular<R3> *qf=QF_Simplex<R3>(exact);//QF_Tria_exact(exact);
    if(verbosity>99 && qf ) cout << "   QF Tet  n:" << qf->n << " exact = " << exact <<  endl;
    if(qf) return *qf;
    /*
     if( QuadratureFormular_T_1.exact >= exact ) return QuadratureFormular_T_1;
     if( QuadratureFormular_T_2.exact >= exact ) return QuadratureFormular_T_2;
     if( QuadratureFormular_T_5.exact >= exact ) return QuadratureFormular_T_5;
     if( QuadratureFormular_T_7.exact >= exact ) return QuadratureFormular_T_7;
     if( QuadratureFormular_T_9.exact >= exact ) return QuadratureFormular_T_9;
     */
    static long count = 0;
    if(verbosity > 1 && count++ < 5)
    cerr << "Warning :  Max Order of the Quadrature Formular Tet  is 6 and expect: " << exact+1
    <<  endl;
    //  ExecError(" We find  no Quadrature Formular on Tet for this  order: too high");
    return QuadratureFormular_Tet_5;
}

const Fem2D::QuadratureFormular & CDomainOfIntegration::FIT(Stack stack) const
{
    using namespace Fem2D;
    if (nargs[0]) return  *GetAny<const Fem2D::QuadratureFormular *>((*nargs[0])(stack));
    int exact = 5;
    if (nargs[2]) exact=  GetAny<long>((*nargs[2])(stack))-1;
    QuadratureFormular *qf=QF_Simplex<R2>(exact);//QF_Tria_exact(exact);
    if(verbosity>99 && qf ) cout << "   QF Tria  n:" << qf->n << " exact = " << exact <<  endl;
    if(qf) return *qf;
    /*
     if( QuadratureFormular_T_1.exact >= exact ) return QuadratureFormular_T_1;
     if( QuadratureFormular_T_2.exact >= exact ) return QuadratureFormular_T_2;
     if( QuadratureFormular_T_5.exact >= exact ) return QuadratureFormular_T_5;
     if( QuadratureFormular_T_7.exact >= exact ) return QuadratureFormular_T_7;
     if( QuadratureFormular_T_9.exact >= exact ) return QuadratureFormular_T_9;
     */
    cerr << " Order of the Quadature Formular: order = " << exact+1 << " exact = " << exact << endl;
    ExecError("Sorry,  we find  no Quadrature Formular on Triangle for this  order: too high.");
    return QuadratureFormular_T_1;
}
const Fem2D::QuadratureFormular1d & CDomainOfIntegration::FIE(Stack stack) const
{
    using namespace Fem2D;
    if (nargs[1]) return  *GetAny<const Fem2D::QuadratureFormular1d *>((*nargs[1])(stack));
    int exact = 5;
    if (nargs[2]) exact=  GetAny<long>((*nargs[2])(stack))-1;
    QuadratureFormular1d *qf=QF_Simplex<R1>(exact);//QF_1d_exact(exact);
    if(verbosity>99 && qf ) cout << "   QF 1d  n:" << qf->n << " exact = " << exact <<  endl;
    if(qf) return *qf;
    /*
     if( 1 >= exact ) return QF_GaussLegendre1;
     if( 3 >= exact ) return QF_GaussLegendre2;
     if( 5 >= exact ) return QF_GaussLegendre3;
     if( 7 >= exact ) return QF_GaussLegendre4;
     if( 9 >= exact ) return QF_GaussLegendre5;
     */
    cerr << " Ordre of the Integration Formular on Edge, order = " << exact+1 << " exact = " << exact << endl;
    ExecError(" We find  no Quadrature Formular on Edge  for this  order:  too high.");
    return QF_GaussLegendre1;
}


namespace Fem2D {


    // general template
    template  void AssembleLinearForm<double>(Stack stack,const Mesh & Th,const FESpace & Vh,KN_<double> * B,const  FormLinear * const l);

    template   void AssembleBilinearForm<double>(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                                                 MatriceCreuse<double>  & A, const  FormBilinear * b  );

    template   void AssembleBilinearForm<double>(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                                                 MatriceMap<double> & A, const  FormBilinear * b  );

    template  void AssembleLinearForm<Complex>(Stack stack,const Mesh & Th,const FESpace & Vh,KN_<Complex> * B,const  FormLinear * const l);

    template   void AssembleBilinearForm<Complex>(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                                                  MatriceCreuse<Complex>  & A, const  FormBilinear * b  );

    template   void AssembleBilinearForm<Complex>(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                                                  MatriceMap<Complex> & A, const  FormBilinear * b  );



    /////// 2d case
    // instantation for type double
    template  bool AssembleVarForm<double,MatriceCreuse<double>,Mesh,FESpace,FESpace >(Stack stack,const Mesh & Th,
                                                                          const FESpace & Uh,const FESpace & Vh,bool sym,
                                                                          MatriceCreuse<double>  * A,KN_<double> * B,const list<C_F0> &largs );
    template  bool AssembleVarForm<double,MatriceMap<double>,Mesh,FESpace,FESpace>(Stack stack,const Mesh & Th,
                                                                                const FESpace & Uh,const FESpace & Vh,bool sym,
                                                                                MatriceMap<double>  * A,KN_<double> * B,const list<C_F0> &largs );
    template   void AssembleBC<double,Mesh,FESpace,FESpace>(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                                               MatriceCreuse<double>  * A,KN_<double> * B,KN_<double> * X, const list<C_F0> &largs , double tgv  );
    // instantation for type complex
    template  bool AssembleVarForm<Complex,MatriceCreuse<Complex>,Mesh,FESpace,FESpace>(Stack stack,const Mesh & Th,
                                                                            const FESpace & Uh,const FESpace & Vh,bool sym,
                                                                            MatriceCreuse<Complex>  * A,KN_<Complex> * B,const list<C_F0> &largs );

    template  bool AssembleVarForm<Complex,MatriceMap<Complex>,Mesh,FESpace,FESpace >(Stack stack,const Mesh & Th,
                                                                                  const FESpace & Uh,const FESpace & Vh,bool sym,
                                                                                  MatriceMap<Complex> * A,KN_<Complex> * B,const list<C_F0> &largs );

    template   void AssembleBC<Complex,Mesh,FESpace,FESpace>(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                                                MatriceCreuse<Complex>  * A,KN_<Complex> * B,KN_<Complex> * X, const list<C_F0> &largs , double tgv  );


    /////// 3D volume case
    // instantation for type double
    template  bool AssembleVarForm<double,MatriceCreuse<double>,Mesh3,FESpace3,FESpace3>(Stack stack,const Mesh3 & Th,
                                                                           const FESpace3 & Uh,const FESpace3 & Vh,bool sym,
                                                                           MatriceCreuse<double>  * A,KN_<double> * B,const list<C_F0> &largs );
    template  bool AssembleVarForm<double,MatriceMap<double>,Mesh3,FESpace3,FESpace3 >(Stack stack,const Mesh3 & Th,
                                                                                 const FESpace3 & Uh,const FESpace3 & Vh,bool sym,
                                                                                 MatriceMap<double>  * A,KN_<double> * B,const list<C_F0> &largs );
    template   void AssembleBC<double,Mesh3,FESpace3,FESpace3>(Stack stack,const Mesh3 & Th,const FESpace3 & Uh,const FESpace3 & Vh,bool sym,
                                                MatriceCreuse<double>  * A,KN_<double> * B,KN_<double> * X, const list<C_F0> &largs , double tgv  );

    // instantation for type complex
    template  bool AssembleVarForm<Complex,MatriceCreuse<Complex>,Mesh3,FESpace3,FESpace3>(Stack stack,const Mesh3 & Th,
                                                                             const FESpace3 & Uh,const FESpace3 & Vh,bool sym,
                                                                             MatriceCreuse<Complex>  * A,KN_<Complex> * B,const list<C_F0> &largs );
    template  bool AssembleVarForm<Complex,MatriceMap<Complex>,Mesh3,FESpace3,FESpace3>(Stack stack,const Mesh3 & Th,
                                                                                   const FESpace3 & Uh,const FESpace3 & Vh,bool sym,
                                                                                   MatriceMap<Complex> * A,KN_<Complex> * B,const list<C_F0> &largs );
    template   void AssembleBC<Complex,Mesh3,FESpace3,FESpace3>(Stack stack,const Mesh3 & Th,const FESpace3 & Uh,const FESpace3 & Vh,bool sym,
                                                 MatriceCreuse<Complex>  * A,KN_<Complex> * B,KN_<Complex> * X, const list<C_F0> &largs , double tgv  );






    /////// 3D surface case
    // instantation for type double

    template  bool AssembleVarForm<double,MatriceCreuse<double>,MeshS,FESpaceS,FESpaceS>(Stack stack,const MeshS & Th,
                                                                           const FESpaceS & Uh,const FESpaceS & Vh,bool sym,
                                                                           MatriceCreuse<double>  * A,KN_<double> * B,const list<C_F0> &largs );
    template  bool AssembleVarForm<double,MatriceMap<double>,MeshS,FESpaceS,FESpaceS>(Stack stack,const MeshS & Th,
                                                                                 const FESpaceS & Uh,const FESpaceS & Vh,bool sym,
                                                                                 MatriceMap<double>  * A,KN_<double> * B,const list<C_F0> &largs );
    template   void AssembleBC<double,MeshS,FESpaceS,FESpaceS>(Stack stack,const MeshS & Th,const FESpaceS & Uh,const FESpaceS & Vh,bool sym,
                                                MatriceCreuse<double>  * A,KN_<double> * B,KN_<double> * X, const list<C_F0> &largs , double tgv  );

    // instantation for type complex
    template  bool AssembleVarForm<Complex,MatriceCreuse<Complex>,MeshS,FESpaceS,FESpaceS>(Stack stack,const MeshS & Th,
                                                                             const FESpaceS & Uh,const FESpaceS & Vh,bool sym,
                                                                             MatriceCreuse<Complex>  * A,KN_<Complex> * B,const list<C_F0> &largs );
    template  bool AssembleVarForm<Complex,MatriceMap<Complex>,MeshS,FESpaceS,FESpaceS>(Stack stack,const MeshS & Th,
                                                                                   const FESpaceS & Uh,const FESpaceS & Vh,bool sym,
                                                                                   MatriceMap<Complex> * A,KN_<Complex> * B,const list<C_F0> &largs );
    template   void AssembleBC<Complex,MeshS,FESpaceS,FESpaceS>(Stack stack,const MeshS & Th,const FESpaceS & Uh,const FESpaceS & Vh,bool sym,
                                                MatriceCreuse<Complex>  * A,KN_<Complex> * B,KN_<Complex> * X, const list<C_F0> &largs , double tgv
                                                );


    /////// 3D  curve
    // instantation for type double

    template  bool AssembleVarForm<double,MatriceCreuse<double>,MeshL,FESpaceL,FESpaceL>(Stack stack,const MeshL & Th,
                                                                           const FESpaceL & Uh,const FESpaceL & Vh,bool sym,
                                                                           MatriceCreuse<double>  * A,KN_<double> * B,const list<C_F0> &largs );
    template  bool AssembleVarForm<double,MatriceMap<double>,MeshL,FESpaceL,FESpaceL>(Stack stack,const  FESpaceL::Mesh & Th,
                                                                        const FESpaceL & Uh,const FESpaceL & Vh,bool sym,
                                                                        MatriceMap<double>  * A,KN_<double> * B,const list<C_F0> &largs );
    template   void AssembleBC<double,MeshL,FESpaceL,FESpaceL>(Stack stack,const MeshL & Th,const FESpaceL & Uh,const FESpaceL & Vh,bool sym,
                                                MatriceCreuse<double>  * A,KN_<double> * B,KN_<double> * X, const list<C_F0> &largs , double tgv  );

    // instantation for type complex
    template  bool AssembleVarForm<Complex,MatriceCreuse<Complex>,MeshL,FESpaceL,FESpaceL>(Stack stack,const MeshL & Th,
                                                                             const FESpaceL & Uh,const FESpaceL & Vh,bool sym,
                                                                             MatriceCreuse<Complex>  * A,KN_<Complex> * B,const list<C_F0> &largs );
    template  bool AssembleVarForm<Complex,MatriceMap<Complex>,MeshL,FESpaceL,FESpaceL>(Stack stack,const MeshL & Th,
                                                                          const FESpaceL & Uh,const FESpaceL & Vh,bool sym,
                                                                          MatriceMap<Complex> * A,KN_<Complex> * B,const list<C_F0> &largs );
    template   void AssembleBC<Complex,MeshL,FESpaceL,FESpaceL>(Stack stack,const MeshL & Th,const FESpaceL & Uh,const FESpaceL & Vh,bool sym,
                                                 MatriceCreuse<Complex>  * A,KN_<Complex> * B,KN_<Complex> * X, const list<C_F0> &largs , double tgv
                                                 );

    /////// 3D  curve / 2D on meshL
    // instantation for type double

    template bool AssembleVarForm<double,MatriceCreuse<double>,MeshL,FESpaceL,FESpace>(Stack stack,const MeshL & Th,
                                                                           const FESpaceL & Uh,const FESpace & Vh,bool sym,
                                                                           MatriceCreuse<double>  * A,KN_<double> * B,const list<C_F0> &largs );
    template bool AssembleVarForm<double,MatriceMap<double>,MeshL,FESpaceL,FESpace>(Stack stack,const  FESpaceL::Mesh & Th,
                                                                        const FESpaceL & Uh,const FESpace & Vh,bool sym,
                                                                        MatriceMap<double>  * A,KN_<double> * B,const list<C_F0> &largs );
    template void AssembleBC<double,MeshL,FESpaceL,FESpace>(Stack stack,const MeshL & Th,const FESpaceL & Uh,const FESpace & Vh,bool sym,
                                                MatriceCreuse<double>  * A,KN_<double> * B,KN_<double> * X, const list<C_F0> &largs , double tgv  );

    // instantation for type complex
    template bool AssembleVarForm<Complex,MatriceCreuse<Complex>,MeshL,FESpaceL,FESpace>(Stack stack,const MeshL & Th,
                                                                             const FESpaceL & Uh,const FESpace & Vh,bool sym,
                                                                             MatriceCreuse<Complex>  * A,KN_<Complex> * B,const list<C_F0> &largs );
    template bool AssembleVarForm<Complex,MatriceMap<Complex>,MeshL,FESpaceL,FESpace>(Stack stack,const MeshL & Th,
                                                                          const FESpaceL & Uh,const FESpace & Vh,bool sym,
                                                                          MatriceMap<Complex> * A,KN_<Complex> * B,const list<C_F0> &largs );
    template void AssembleBC<Complex,MeshL,FESpaceL,FESpace>(Stack stack,const MeshL & Th,const FESpaceL & Uh,const FESpace & Vh,bool sym,
                                                 MatriceCreuse<Complex>  * A,KN_<Complex> * B,KN_<Complex> * X, const list<C_F0> &largs , double tgv
                                                 );
    /////// 2D / 3D  curve on meshL
    // instantation for type double

    template bool AssembleVarForm<double,MatriceCreuse<double>,MeshL,FESpace,FESpaceL>(Stack stack,const MeshL & Th,
                                                                            const FESpace & Uh,const FESpaceL & Vh,bool sym,
                                                                            MatriceCreuse<double>  * A,KN_<double> * B,const list<C_F0> &largs );
    template bool AssembleVarForm<double,MatriceMap<double>,MeshL,FESpace,FESpaceL>(Stack stack,const MeshL & Th,
                                                                          const FESpace & Uh,const FESpaceL & Vh,bool sym,
                                                                          MatriceMap<double>  * A,KN_<double> * B,const list<C_F0> &largs );
    template void AssembleBC<double,MeshL,FESpace,FESpaceL>(Stack stack,const MeshL & Th,const FESpace & Uh,const FESpaceL & Vh,bool sym,
                                                 MatriceCreuse<double>  * A,KN_<double> * B,KN_<double> * X, const list<C_F0> &largs , double tgv  );

    // instantation for type complex
    template bool AssembleVarForm<Complex,MatriceCreuse<Complex>,MeshL,FESpace,FESpaceL>(Stack stack,const MeshL & Th,
                                                                             const FESpace & Uh,const FESpaceL & Vh,bool sym,
                                                                             MatriceCreuse<Complex>  * A,KN_<Complex> * B,const list<C_F0> &largs );
    template bool AssembleVarForm<Complex,MatriceMap<Complex>,MeshL,FESpace,FESpaceL>(Stack stack,const MeshL & Th,
                                                                        const FESpace & Uh,const FESpaceL & Vh,bool sym,
                                                                        MatriceMap<Complex> * A,KN_<Complex> * B,const list<C_F0> &largs );
    template void AssembleBC<Complex,MeshL,FESpace,FESpaceL>(Stack stack,const MeshL & Th,const FESpace & Uh,const FESpaceL & Vh,bool sym,
                                               MatriceCreuse<Complex>  * A,KN_<Complex> * B,KN_<Complex> * X, const list<C_F0> &largs , double tgv );
   /////// 3D Surf / 3D volume on meshS
   // instantation for type double

   template bool AssembleVarForm<double,MatriceCreuse<double>,MeshS,FESpaceS,FESpace3>(Stack stack,const MeshS & Th,
                                                                           const FESpaceS & Uh,const FESpace3 & Vh,bool sym,
                                                                           MatriceCreuse<double>  * A,KN_<double> * B,const list<C_F0> &largs );
   template bool AssembleVarForm<double,MatriceMap<double>,MeshS,FESpaceS,FESpace3>(Stack stack,const MeshS & Th,
                                                                         const FESpaceS & Uh,const FESpace3 & Vh,bool sym,
                                                                         MatriceMap<double>  * A,KN_<double> * B,const list<C_F0> &largs );
   template void AssembleBC<double,MeshS,FESpaceS,FESpace3>(Stack stack,const MeshS & Th,const FESpaceS & Uh,const FESpace3 & Vh,bool sym,
                                                MatriceCreuse<double>  * A,KN_<double> * B,KN_<double> * X, const list<C_F0> &largs , double tgv  );

   // instantation for type complex
   template bool AssembleVarForm<Complex,MatriceCreuse<Complex>,MeshS,FESpaceS,FESpace3>(Stack stack,const MeshS & Th,
                                                                            const FESpaceS & Uh,const FESpace3 & Vh,bool sym,
                                                                            MatriceCreuse<Complex>  * A,KN_<Complex> * B,const list<C_F0> &largs );
   template bool AssembleVarForm<Complex,MatriceMap<Complex>,MeshS,FESpaceS,FESpace3>(Stack stack,const MeshS & Th,
                                                                       const FESpaceS & Uh,const FESpace3 & Vh,bool sym,
                                                                       MatriceMap<Complex> * A,KN_<Complex> * B,const list<C_F0> &largs );
   template void AssembleBC<Complex,MeshS,FESpaceS,FESpace3>(Stack stack,const MeshS & Th,const FESpaceS & Uh,const FESpace3 & Vh,bool sym,
                                              MatriceCreuse<Complex>  * A,KN_<Complex> * B,KN_<Complex> * X, const list<C_F0> &largs , double tgv );
   /////// 3D volume / 3D Surf on meshS
   // instantation for type double

   template bool AssembleVarForm<double,MatriceCreuse<double>,MeshS,FESpace3,FESpaceS>(Stack stack,const MeshS & Th,
                                                                           const FESpace3 & Uh,const FESpaceS & Vh,bool sym,
                                                                           MatriceCreuse<double>  * A,KN_<double> * B,const list<C_F0> &largs );
   template bool AssembleVarForm<double,MatriceMap<double>,MeshS,FESpace3,FESpaceS>(Stack stack,const MeshS & Th,
                                                                         const FESpace3 & Uh,const FESpaceS & Vh,bool sym,
                                                                         MatriceMap<double>  * A,KN_<double> * B,const list<C_F0> &largs );
   template void AssembleBC<double,MeshS,FESpace3,FESpaceS>(Stack stack,const MeshS & Th,const FESpace3 & Uh,const FESpaceS & Vh,bool sym,
                                                MatriceCreuse<double>  * A,KN_<double> * B,KN_<double> * X, const list<C_F0> &largs , double tgv  );

   // instantation for type complex
   template bool AssembleVarForm<Complex,MatriceCreuse<Complex>,MeshS,FESpace3,FESpaceS>(Stack stack,const MeshS & Th,
                                                                            const FESpace3 & Uh,const FESpaceS & Vh,bool sym,
                                                                            MatriceCreuse<Complex>  * A,KN_<Complex> * B,const list<C_F0> &largs );
   template bool AssembleVarForm<Complex,MatriceMap<Complex>,MeshS,FESpace3,FESpaceS>(Stack stack,const MeshS & Th,
                                                                       const FESpace3 & Uh,const FESpaceS & Vh,bool sym,
                                                                       MatriceMap<Complex> * A,KN_<Complex> * B,const list<C_F0> &largs );
   template void AssembleBC<Complex,MeshS,FESpace3,FESpaceS>(Stack stack,const MeshS & Th,const FESpace3 & Uh,const FESpaceS & Vh,bool sym,
                                              MatriceCreuse<Complex>  * A,KN_<Complex> * B,KN_<Complex> * X, const list<C_F0> &largs , double tgv );

   /////// 3D  curve / 3D Surf on meshL
   // instantation for type double

   template bool AssembleVarForm<double,MatriceCreuse<double>,MeshL,FESpaceL,FESpaceS>(Stack stack,const MeshL & Th,
                                                                          const FESpaceL & Uh,const FESpaceS & Vh,bool sym,
                                                                          MatriceCreuse<double>  * A,KN_<double> * B,const list<C_F0> &largs );
   template bool AssembleVarForm<double,MatriceMap<double>,MeshL,FESpaceL,FESpaceS>(Stack stack,const MeshL & Th,
                                                                       const FESpaceL & Uh,const FESpaceS & Vh,bool sym,
                                                                       MatriceMap<double>  * A,KN_<double> * B,const list<C_F0> &largs );
   template void AssembleBC<double,MeshL,FESpaceL,FESpaceS>(Stack stack,const MeshL & Th,const FESpaceL & Uh,const FESpaceS & Vh,bool sym,
                                               MatriceCreuse<double>  * A,KN_<double> * B,KN_<double> * X, const list<C_F0> &largs , double tgv  );

   // instantation for type complex
   template bool AssembleVarForm<Complex,MatriceCreuse<Complex>,MeshL,FESpaceL,FESpaceS>(Stack stack,const MeshL & Th,
                                                                            const FESpaceL & Uh,const FESpaceS & Vh,bool sym,
                                                                            MatriceCreuse<Complex>  * A,KN_<Complex> * B,const list<C_F0> &largs );
   template bool AssembleVarForm<Complex,MatriceMap<Complex>,MeshL,FESpaceL,FESpaceS>(Stack stack,const MeshL & Th,
                                                                         const FESpaceL & Uh,const FESpaceS & Vh,bool sym,
                                                                         MatriceMap<Complex> * A,KN_<Complex> * B,const list<C_F0> &largs );
   template void AssembleBC<Complex,MeshL,FESpaceL,FESpaceS>(Stack stack,const MeshL & Th,const FESpaceL & Uh,const FESpaceS & Vh,bool sym,
                                                MatriceCreuse<Complex>  * A,KN_<Complex> * B,KN_<Complex> * X, const list<C_F0> &largs , double tgv
                                                );
   /////// 3D Surf / 3D  curve on meshL
   // instantation for type double

   template bool AssembleVarForm<double,MatriceCreuse<double>,MeshL,FESpaceS,FESpaceL>(Stack stack,const MeshL & Th,
                                                                           const FESpaceS & Uh,const FESpaceL & Vh,bool sym,
                                                                           MatriceCreuse<double>  * A,KN_<double> * B,const list<C_F0> &largs );
   template bool AssembleVarForm<double,MatriceMap<double>,MeshL,FESpaceS,FESpaceL>(Stack stack,const MeshL & Th,
                                                                         const FESpaceS & Uh,const FESpaceL & Vh,bool sym,
                                                                         MatriceMap<double>  * A,KN_<double> * B,const list<C_F0> &largs );
   template void AssembleBC<double,MeshL,FESpaceS,FESpaceL>(Stack stack,const MeshL & Th,const FESpaceS & Uh,const FESpaceL & Vh,bool sym,
                                                MatriceCreuse<double>  * A,KN_<double> * B,KN_<double> * X, const list<C_F0> &largs , double tgv  );

   // instantation for type complex
   template bool AssembleVarForm<Complex,MatriceCreuse<Complex>,MeshL,FESpaceS,FESpaceL>(Stack stack,const MeshL & Th,
                                                                            const FESpaceS & Uh,const FESpaceL & Vh,bool sym,
                                                                            MatriceCreuse<Complex>  * A,KN_<Complex> * B,const list<C_F0> &largs );
   template bool AssembleVarForm<Complex,MatriceMap<Complex>,MeshL,FESpaceS,FESpaceL>(Stack stack,const MeshL & Th,
                                                                       const FESpaceS & Uh,const FESpaceL & Vh,bool sym,
                                                                       MatriceMap<Complex> * A,KN_<Complex> * B,const list<C_F0> &largs );
   template void AssembleBC<Complex,MeshL,FESpaceS,FESpaceL>(Stack stack,const MeshL & Th,const FESpaceS & Uh,const FESpaceL & Vh,bool sym,
                                              MatriceCreuse<Complex>  * A,KN_<Complex> * B,KN_<Complex> * X, const list<C_F0> &largs , double tgv );
  
}

template class Call_FormLinear<v_fes>;
template class Call_FormLinear<v_fes3>;
template class Call_FormLinear<v_fesS>;
template class Call_FormLinear<v_fesL>;
template class Call_FormLinear<vect_generic_v_fes>; // Morice: added vector FESpace (composite FESpace)

template class Call_FormBilinear<v_fes,v_fes>;
template class Call_FormBilinear<v_fes3,v_fes3>;
template class Call_FormBilinear<v_fesS,v_fesS>;
template class Call_FormBilinear<v_fesL,v_fesL>;

template class Call_FormBilinear<v_fesL, v_fesS>; //  3D curve / 3D Surf on meshL and bem
template class Call_FormBilinear<v_fesS, v_fesL>; //  3D Surf / 3D curve on meshL
template class Call_FormBilinear<v_fesL, v_fes>;  //  3D curve / 2D on meshL
template class Call_FormBilinear<v_fes, v_fesL>;  //  2D / 3D curve on meshL
template class Call_FormBilinear<v_fesS, v_fes3>;  //  3D Surf / 3D volume on meshS
template class Call_FormBilinear<v_fes3, v_fesS>;  //  3D volume / 3D Surf on meshS
template class Call_FormBilinear<v_fesS, v_fes>;

template class Call_FormBilinear<vect_generic_v_fes, vect_generic_v_fes>; // Morice: added vector FESpace (composite FESpace)
/*

#ifndef FFLANG
#ifdef PARALLELE

#define BOOST_NO_CXX17_IF_CONSTEXPR
#include <ff++.hpp>
#include <AFunction_ext.hpp>
#include <lgfem.hpp>
#include <R3.hpp>

#include <htool/htool.hpp>

// include the bemtool library .... path define in where library
//#include <bemtool/operator/block_op.hpp>
#include <bemtool/tools.hpp>
#include <bemtool/fem/dof.hpp>
#include <bemtool/operator/operator.hpp>
#include <bemtool/miscellaneous/htool_wrap.hpp>
// #include "PlotStream.hpp"

#include "common.hpp"

// extern FILE *ThePlotStream;

using namespace std;
//using namespace htool;
//using namespace bemtool;

#include <type_traits>
/*
typedef   LinearComb<MGauche,C_F0> Finconnue;
typedef   LinearComb<MDroit,C_F0> Ftest;
typedef  const Finconnue  finconnue;
typedef  const Ftest ftest;
class CDomainOfIntegration;
class FormBilinear;
*/
/*
#include "bem.hpp"

#endif
#endif
*/
/*
template<class R>
AnyType OpMatrixtoBilinearFormVG<R>::Op::operator()(Stack stack)  const
{
  assert(b && b->nargs);

  pvectgenericfes  * pUh= GetAny<pvectgenericfes *>((*b->euh)(stack));
  pvectgenericfes  * pVh= GetAny<pvectgenericfes *>((*b->evh)(stack));
  ffassert( *pUh && *pVh ); 

  int NpUh = (*pUh)->N; // number of fespace in pUh
  int NpVh = (*pVh)->N; // number of fespace in pVh

  KN<int> UhNbOfDf = (*pUh)->vectOfNbOfDF();
  KN<int> VhNbOfDf = (*pVh)->vectOfNbOfDF();

  KN<int> UhNbItem = (*pUh)->vectOfNbitem();
  KN<int> VhNbItem = (*pVh)->vectOfNbitem();

  KN<int> beginBlockUh(NpUh); // index of the first elment of a block
  KN<int> beginBlockVh(NpVh);

  // loop over the index 
  int UhtotalNbItem=0;
  cout << "finc: FESpace" << endl;
  for(int i=0; i< NpUh; i++){
    cout << "component i=" << i << ", NbItem["<<i<<"]=" <<  UhNbItem[i] << ", NbDof["<<i<<"]=" << UhNbOfDf[i] << endl;
    beginBlockUh[i] = UhtotalNbItem;
    UhtotalNbItem += UhNbItem[i];
    ffassert( UhNbItem[i] > 0 && UhNbOfDf[i] > 0);
  }

  int VhtotalNbItem=0;
  cout << "ftest: FESpace" << endl;
  for(int i=0; i< NpVh; i++){
    cout << "component i=" << i << ", NbItem["<<i<<"]=" <<  VhNbItem[i] << ", NbDof["<<i<<"]=" << VhNbOfDf[i] << endl;
    beginBlockVh[i] = VhtotalNbItem;
    VhtotalNbItem += VhNbItem[i];
    ffassert( VhNbItem[i] > 0 && VhNbOfDf[i] > 0);
  }
  
  // index for the construction of the block
  KN<int> indexBlockUh(UhtotalNbItem);
  KN<int> localIndexInTheBlockUh(UhtotalNbItem);
  { 
    // ========================
    //
    // varf([u0,u1,...,u4], ... ) 
    // varf([ [u0_blk1,u1_blk1],[u0_blk2,u1_blk2,u2_blk2] ], ... ) 

    // For u4, on a :: current_index = 4
    //              :: indexBlockUh = 2
    //              :: localIndexInThBlock = 3
    int current_index=0;
    for(int i=0; i<NpUh; i++){
      for(int j=0; j<UhNbItem[i]; j++){
        indexBlockUh[current_index] = i;
        localIndexInTheBlockUh[current_index] = j;
        current_index++;
      }
    }
    ffassert(current_index==UhtotalNbItem);
  }

  KN<int> indexBlockVh(VhtotalNbItem);
  KN<int> localIndexInTheBlockVh(VhtotalNbItem);
  { 
    int current_index=0;
    for(int i=0; i<NpVh; i++){
      for(int j=0; j<VhNbItem[i]; j++){
        indexBlockVh[current_index] = i;
        localIndexInTheBlockVh[current_index] = j;
        current_index++;
      }
    }
    ffassert(current_index==VhtotalNbItem);
  }
  cout <<"========================================================" << endl;
  cout <<"=                                                      =" << endl;
  cout <<"= indexBlockUh=                                        =" << endl;
  cout << indexBlockUh << endl;
  cout <<"=                                                      =" << endl;
  cout <<"= localIndexInTheBlockUh=                              =" << endl;
  cout << localIndexInTheBlockUh << endl;

  cout <<"========================================================" << endl;
  cout <<"=                                                      =" << endl;
  cout <<"= indexBlockVh=                                        =" << endl;
  cout << indexBlockVh << endl;
  cout <<"=                                                      =" << endl;
  cout <<"= localIndexInTheBlockVh=                              =" << endl;
  cout << localIndexInTheBlockVh << endl;


  #ifndef FFLANG
  #ifdef PARALLELE
  cout << "====    define parallele  =====" << endl;
  exit(0);
  #endif
  #endif
  


  //
  const list<C_F0> & largs=b->largs; 

  KNM< list<C_F0> > block_largs( (long)NpUh, (long)NpVh );

  // impression des information de la composition largs
  list<C_F0>::const_iterator ii,ib=largs.begin(),ie=largs.end(); 


  // necessaire :: UhtotalNbItem, indexBlockUh

  // Loop to put each term of the varf in each correct block

  // loop over largs information 
  cout << "loop over the integral" << endl;

  int count_integral = 0;
  for (ii=ib;ii != ie;ii++) {
    count_integral++;
    cout <<"========================================================" << endl;
    cout <<"=                                                      =" << endl;
    cout << "reading the " << count_integral << "-th term of the variational form used to define the matrix" << endl;
    Expression e=ii->LeftValue();
    aType r = ii->left();
    cout << "e=" << e << ", " << "r=" << r << endl;
    cout <<"=                                                      =" << endl;

    // ***************************************
    // Case FormBillinear
    // ***************************************
    if (r==atype<const  FormBilinear *>() ){
      const FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
      const CDomainOfIntegration & di= *bb->di;

      cout << "di.kind=" << di.kind << endl;
      cout << "di.dHat=" << di.dHat << endl;
      cout << "di.d=" << di.d << endl;
      cout << "di.Th=" << di.Th << endl;

      int    d = di.d;
      int dHat = di.dHat;

      // Sert a verifier que "*bb->di->Th" est du bon type ==> A enlever
      
      // Recuperation du pointeur sur le maillage de l'integrale
      //if(d==2){ // 3d
      //  pmesh Thtest=GetAny<pmesh >((*bb->di->Th)(stack));
      //  cout << "pointeur du Th de l'integrale =" << Thtest << endl;
      //}
      //else if(d==3 && dHat==3){
      //  pmesh3 Thtest=GetAny<pmesh3 >((*bb->di->Th)(stack));
      //  cout << "pointeur du Th de l'integrale =" << Thtest << endl;
      //}
      //else if(d==3 && dHat==2){
      //  pmeshS Thtest=GetAny<pmeshS >((*bb->di->Th)(stack));
      //  cout << "pointeur du Th de l'integrale =" << Thtest << endl;
      //}
      //else if(d==3 && dHat==1){
      //  pmeshL Thtest=GetAny<pmeshL >((*bb->di->Th)(stack));
      //  cout << "pointeur du Th de l'integrale =" << Thtest << endl;
      //}
      //else{ ffassert(0); }// a faire
      
      BilinearOperator * Op=const_cast<  BilinearOperator *>(bb->b);
      if (Op == NULL) {
        if(mpirank == 0) cout << "dynamic_cast error" << endl; 
        ffassert(0);
      }
      
      size_t Opsize= Op->v.size();
      cout << " loop over the term inside the integral" << endl;
      cout << " Number of term in the integral:: Op->v.size()=" << Op->v.size() << endl;

    
      // index to check if a integral is defined on multi block
      int indexOfBlockUh = -1; // A changer de nom
      int indexOfBlockVh = -1; // A changer de nom
      for(size_t jj=0; jj<Opsize; jj++){
        // attention la fonction test donne la ligne
        //  et la fonction test est en second
        BilinearOperator::K ll = Op->v[jj];
        pair<int,int> finc(ll.first.first), ftest(ll.first.second);
        cout << " operateur jj= " << jj << endl;
        cout << " FormBilinear: number of unknown finc=" <<  finc.first << " ,ftest= " << ftest.first << endl;
        cout << " FormBilinear: operator order finc   =" << finc.second << " ,ftest= " << ftest.second << endl; // ordre   only op_id=0
        
        // Fred fait peut être un message après ????
        // verification que la taille des tableaux des fonctions tests et de la fonction inconnue``
        // sont correctes.  
        ffassert( -1  < finc.first  && finc.first < UhtotalNbItem);
        ffassert( -1  < ftest.first && ftest.first < VhtotalNbItem);

        // finc.first : index de component de la fonction inconnue
        // ftest.first: index de component de la fonction test
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // finc.second : renvoie l'index du type de l'operateur: Id, dx(), dy(), dz(), dxx(), dxy(), ...
        //
        // la liste des index des operateurs est definis dans [[????]].
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // exemple vf([u1,u2,..,u30],[v1,v2,...,v30]) = int2d(Th)(dx(u20)*v15)
        //      finc.first  = 20 , ftest.first = 15
        //      finc.second = 1 , ftest.second = 0

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if( jj== 0 ){
          indexOfBlockUh = indexBlockUh[finc.first];
          indexOfBlockVh = indexBlockVh[ftest.first];
        }
        else if( indexOfBlockUh != indexBlockUh[finc.first] ){
          cerr << "The " << count_integral <<"-th integral(s) contains the constribution of two different blocks:" << endl;
          cerr << "the first term correspond to block (" << indexOfBlockUh << " , " <<  indexOfBlockVh << ")" << endl;
          cerr << "the "<<jj <<"-th term correspond to block (" << indexBlockUh(finc.first) << " " <<  indexBlockVh(ftest.first) << ")" << endl;
          cerr << "You need to separate the integral in individual part." << endl;
          cerr << "Remark: scalar product in N dimension correspond to N terms inside the integral" << endl;
          cerr << "A ameliorer Jacques." << endl;
          ffassert(0);
        }
        else if( indexOfBlockVh != indexBlockVh[ftest.first] ){
          cerr << "The " << count_integral <<"-th integral(s) contains the constribution of two different blocks:" << endl;
          cerr << "the first term correspond to block (" << indexOfBlockUh << " , " <<  indexOfBlockVh << ")" <<endl;
          cerr << "the "<<jj <<"-th term correspond to block (" << indexBlockUh(finc.first) << ", " <<  indexBlockVh(ftest.first) << ")" << endl;
          cerr << "You need to separate the integral in individual part." << endl;
          cerr << "Remark: scalar product in N dimension correspond to N terms inside the integral" << endl;
          cerr << "A ameliorer Jacques." << endl;
          ffassert(0);
        }

        ffassert( indexOfBlockUh == indexBlockUh(finc.first) );
        ffassert( indexOfBlockVh == indexBlockVh(ftest.first) );
      }

    ffassert( indexOfBlockUh >= 0 && indexOfBlockVh >= 0);
    
    // A faire :: recuperation des éléments pour chacun des blocs
    // Actuellement, on associe une intégrale par block ==> 
    
    // change the index of the block

    //changeIndexFunctionInconnue(*Op, index_operator_finc, new_index_funct_finc );
  
    //changeIndexFunctionTest(*Op, index_operator_ftest, new_index_funct_ftest  );
    
    block_largs(indexOfBlockUh,indexOfBlockVh).push_back(*ii);

    cout << "The " << count_integral <<"-th integral(s) is added to the block (" << indexOfBlockUh << " , " <<  indexOfBlockVh << ")" <<endl;

    }
#ifndef FFLANG
#ifdef PARALLELE
    // ******************************************
    // Case BemKFormBilinear (KERNEL FORM ONLY)
    // ******************************************
    else if (r==atype<const BemFormBilinear *>() ){
      BemFormBilinear * bbtmp= dynamic_cast< BemFormBilinear *>(e);
      int VVFBEM = bbtmp->type;

      if(VVFBEM ==1){
        BemKFormBilinear * bb=new BemKFormBilinear(*dynamic_cast<const BemKFormBilinear *>(e));
        FoperatorKBEM * b=const_cast<  FoperatorKBEM *>(bb->b);
        if (b == NULL) { if(mpirank == 0) cout << "dynamic_cast error" << endl; exit(0);}

        int indexOfBlockUh = -1; // A changer de nom
        int indexOfBlockVh = -1; // A changer de nom


        // loop over the index of finconnue
        LOperaG * OpG = const_cast<LOperaG *>(b->fi);
        ffassert( OpG->v.size() == 1);
        size_t jj =0;
        for (LOperaG::const_iterator lop=OpG->v.begin();lop!=OpG->v.end();lop++){

          LOperaG::K lf(*lop);
          pair<int,int> finc(lf.first);
          cout << " operateur jj= " << jj << endl;
          cout << " BemFormLinear: number of unknown finc= " << finc.first << endl;
          ffassert( -1  < finc.first && finc.first < UhtotalNbItem);     // check the index 

          if( jj== 0 ){
            indexOfBlockUh = indexBlockUh[finc.first];
          }
          else if( indexOfBlockUh != indexBlockUh[finc.first] ){
            cerr << "The " << count_integral <<"-th term of the varitional form contains the constribution of two different FESpace:" << endl;
            cerr << "This terms correspond to a BEM integral terms" << endl;
            cerr << "the first term correspond to element " << indexOfBlockUh << " of the Composite FESpace (Finconnu)." << endl;
            cerr << "the "<< jj <<"-th term correspond to element " << indexBlockUh(finc.first) << endl;
            cerr << "In a composite FESpace, you need to define a BEM integral for each FESpace individually." << endl;
            cerr << "A ameliorer Jacques." << endl;
            ffassert(0);
          }
          jj+=1;
        }


        // Loop over the index of ftest
        LOperaD * OpD = const_cast<LOperaD *>(b->ft);
        ffassert( OpD->v.size() == 1);
        jj =0; // reinitialisation ton zero
        for (LOperaD::const_iterator lop=OpD->v.begin();lop!=OpD->v.end();lop++){
          
          LOperaD::K lf(*lop);
          pair<int,int> ftest(lf.first);
          cout << " operateur jj= " << jj << endl;
          cout << " BemFormLinear: number of unknown ftest= " << ftest.first << endl;
          ffassert( -1  < ftest.first && ftest.first < VhtotalNbItem);    // check the index 

          if( jj== 0 ){
            indexOfBlockVh = indexBlockVh[ftest.first];
          }
          else if( indexOfBlockVh != indexBlockVh[ftest.first] ){
            cerr << "The " << count_integral <<"-th term of the varitional form contains the constribution of two different FESpace:" << endl;
            cerr << "This terms correspond to a BEM integral terms" << endl;
            cerr << "the first term correspond to element " << indexOfBlockVh << " of the Composite FESpace (Ftest)." << endl;
            cerr << "the "<< jj <<"-th term correspond to element " << indexBlockVh(ftest.first) << endl;
            cerr << "In a composite FESpace, you need to define a BEM integral term for each FESpace individually." << endl;
            cerr << "A ameliorer Jacques." << endl;
            ffassert(0);
          }
          jj+=1;
        }
        block_largs(indexOfBlockUh,indexOfBlockVh).push_back(*ii); 

      }else if(VVFBEM == 2){
        BemPFormBilinear * bb=new BemPFormBilinear(*dynamic_cast<const BemPFormBilinear *>(e));
        FoperatorPBEM * b=const_cast<  FoperatorPBEM *>(bb->b);
        if (b == NULL) { if(mpirank == 0) cout << "dynamic_cast error" << endl; }

        cerr << " BEM Potential in composite FESpace in construction " << endl;
        ffassert(0);
      }

    }
#endif
#endif
    else if(r == atype<const  BC_set  *>()){
      cout << " BC in variational form " << endl;
      
      const BC_set * bc=dynamic_cast<const  BC_set *>(e);
      
      // index to check if a integral is defined on multi block
      int indexOfBlockUh = -1;
    
      int kk=bc->bc.size();
      for (int k=0;k<kk;k++)
      {
          pair<int,Expression> xx=bc->bc[k];
          ffassert( -1  < xx.first  && xx.first < UhtotalNbItem); // check the value of index of the component of the varf
        
          if( k == 0) indexOfBlockUh = indexBlockUh[xx.first]; // index of the block Uh
          else if( indexOfBlockUh != indexBlockUh[xx.first] ){
            cerr << "The " << count_integral <<"-th term of the varitional form contains the constribution of two different FESpace:" << endl;
            cerr << "This terms correspond to Boundary condition" << endl;
            cerr << "the first term correspond to element " << indexOfBlockUh << " of the Composite FESpace " << endl;
            cerr << "the "<<kk <<"-th term correspond to element " << indexBlockUh(xx.first) << endl;
            cerr << "In a composite FESpace, you need to define a BC for each FESpace individually." << endl;
            cerr << "A ameliorer Jacques." << endl;
            ffassert(0);
          }
      }
      // Added the boundary condition in the largs block
      block_largs(indexOfBlockUh,indexOfBlockUh).push_back(*ii); 
        
      //ffassert(0);
    }
    else{
      #ifndef FFLANG
      #ifdef PARALLELE
      cout << "r=" << r << ", " << atype<const BemFormBilinear *>() << endl;
      #endif
      #endif
      cerr << "Composite FESpace :: bilinear form only " << endl;
      cerr << "      uniquement terme bilineaire + BC " << endl;
      ffassert(0);
    }
  }

  // check the list of <<largs>> for each block
  cout <<"========================================================" << endl;
  cout <<"=                                                      =" << endl;
  cout <<"= check the list of largs of each block                =" << endl;
  cout <<"=                                                      =" << endl;
  for( int i=0; i<NpUh; i++){
    for( int j=0; j<NpVh; j++){
      cout<< " block ( "<< i << " , " << j <<  " ) " ;
      const list<C_F0> & b_largs=block_largs(i,j); 
      cout<< ", size of the list=" << b_largs.size() << endl;
      // impression des information de la composition largs
      list<C_F0>::const_iterator b_ii,b_ib=b_largs.begin(),b_ie=b_largs.end(); 
      for (b_ii=b_ib;b_ii != b_ie;b_ii++){
        Expression e=b_ii->LeftValue();
        aType r = b_ii->left();
        cout << "e=" << e << ", r=" << r << endl;
      }
    }
  }
  
  // Info necessaire :: " block_largs, localIndexInTheBlockUh, localIndexInTheBlockVh, NpUh, NpVh  

  // put the right number of each component of each block
  for( int i=0; i<NpUh; i++){
      for( int j=0; j<NpVh; j++){
        
        const list<C_F0> *b_largs=&block_largs(i,j); 
        list<C_F0>::const_iterator b_ii,b_ib=b_largs->begin(),b_ie=b_largs->end(); 
        for (b_ii=b_ib;b_ii != b_ie;b_ii++){
          Expression e=b_ii->LeftValue();
          aType r = b_ii->left();
          // Case FormBilinear
          if (r==atype<const  FormBilinear *>() ){
            const FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
      
            BilinearOperator * Op=const_cast<  BilinearOperator *>(bb->b);
            if (Op == NULL) {
              if(mpirank == 0) cout << "dynamic_cast error" << endl; 
            ffassert(0);
            }
      
            size_t Opsize= Op->v.size();
            cout << " loop over the term inside the integral" << endl;
            cout << " Number of term in the integral:: Op->v.size()=" << Op->v.size() << endl;

        
            KN<size_t> index_operator_finc(Opsize);
            KN<int>    new_index_funct_finc(Opsize);

            KN<size_t> index_operator_ftest(Opsize);
            KN<int>    new_index_funct_ftest(Opsize);

            // index to check if a integral is defined on multi block
            for(size_t jj=0; jj<Opsize; jj++){
              // attention la fonction test donne la ligne
              //  et la fonction test est en second
              BilinearOperator::K ll = Op->v[jj];
              pair<int,int> finc(ll.first.first), ftest(ll.first.second);

              long jj2= jj;

              index_operator_finc[ jj2] = jj;
              new_index_funct_finc[ jj2] = localIndexInTheBlockUh(finc.first);
            
              index_operator_ftest[ jj2]  = jj;
              new_index_funct_ftest[ jj2] = localIndexInTheBlockVh(ftest.first);

            }
            changeIndexFunctionInconnue(*Op, index_operator_finc, new_index_funct_finc );
      
            changeIndexFunctionTest(*Op, index_operator_ftest, new_index_funct_ftest  );          
          }  

          #ifndef FFLANG
          #ifdef PARALLELE
          // ******************************************
          // Case BemKFormBilinear (KERNEL FORM ONLY)
          // ******************************************
          else if (r==atype<const BemFormBilinear *>() ){
            BemFormBilinear * bbtmp= dynamic_cast< BemFormBilinear *>(e);
            int VVFBEM = bbtmp->type;

            if(VVFBEM ==1){
              BemKFormBilinear * bb=new BemKFormBilinear(*dynamic_cast<const BemKFormBilinear *>(e));
              FoperatorKBEM * b=const_cast<  FoperatorKBEM *>(bb->b);
              if (b == NULL) { if(mpirank == 0) cout << "dynamic_cast error" << endl; exit(0);}

              // loop over the index of finconnue
              LOperaG * OpG = const_cast<LOperaG *>(b->fi);
              ffassert( OpG->v.size() == 1);
              size_t Opsize= OpG->v.size();

              for(size_t jjj=0; jjj<Opsize; jjj++){
                LOperaG::K *lf=&(OpG->v[jjj]);
                OpG->v[jjj].first.first = localIndexInTheBlockUh( OpG->v[jjj].first.first );
              
                pair<int,int> finc(lf->first);
                cout << " new value :: block i,j=" << i << ","<< j << ", operateur jj= " << jjj << endl;
                cout << " BemormBilinear: number of unknown finc = " << finc.first << endl;
                cout << " BemFormBilinear: operator order   finc = " << finc.second << endl; 
              }


              // Loop over the index of ftest
              LOperaD * OpD = const_cast<LOperaD *>(b->ft);
              ffassert( OpD->v.size() == 1);
              Opsize= OpD->v.size();
              for(size_t jjj=0; jjj<Opsize; jjj++){
                LOperaD::K *lf=&(OpD->v[jjj]);
                OpD->v[jjj].first.first = localIndexInTheBlockVh( OpD->v[jjj].first.first );
              
                pair<int,int> ftest(lf->first);
                cout << " new value :: block i,j=" << i << ","<< j << ", operateur jj= " << jjj << endl;
                cout << " BemormBilinear: number of unknown ftest = " << ftest.first << endl;
                cout << " BemFormBilinear: operator order   ftest = " << ftest.second << endl; 
              }
            }
            else if(VVFBEM == 2){
              BemPFormBilinear * bb=new BemPFormBilinear(*dynamic_cast<const BemPFormBilinear *>(e));
              FoperatorPBEM * b=const_cast<  FoperatorPBEM *>(bb->b);
              if (b == NULL) { if(mpirank == 0) cout << "dynamic_cast error" << endl; }

              cerr << " BEM Potential in composite FESpace in construction " << endl;
              ffassert(0);
            }

          }
          #endif
          #endif

          // case BC_set 
          else if(r == atype<const  BC_set  *>()){
            ffassert( i == j ); // diagonal block 
            BC_set * bc=dynamic_cast< BC_set *>(e); // on ne peut pas utiliser " const BC_set * " ou autrement erreur ce ompilation:  Morice

            //KN<int>  new_index_funct_finc( bc.size() );
            int kk=bc->bc.size();
            //pair<int,Expression>  &bc_ib(bc->bc.begin());
            
            for (int k=0;k<kk;k++)
            {
              pair<int,Expression> &xx2= bc->bc[k];
              //new_index_funct_finc[k] = localIndexInTheBlockUh(bc[k].first);
              // change the index of the component to correspond to the index in the block
              xx2.first = localIndexInTheBlockUh(xx2.first);
              //bc->changeNumberOfComponent(k,localIndexInTheBlockUh(xx.first));
              //bc->bc[k].first = localIndexInTheBlockUh( bc->bc[k].first );
            }
          }
        }
        // listOfComponentBilinearForm(*b_largs);
      }
  }
  
  cout <<"========================================================" << endl;
  cout <<"=                                                      =" << endl;
  cout <<"= check the component of Bilinear Form                 =" << endl;
  cout <<"=                                                      =" << endl;
  listOfComponentBilinearForm(largs);
  cout <<"=                                                      =" << endl;
  cout <<"========================================================" << endl;

  //===   Information of the global matrix    ===// 

  // check if we have a square matrix
  bool A_is_square= (void*)pUh == (void*)pVh || ((*pUh)->totalNbOfDF()) == ( (*pVh)->totalNbOfDF()) ;
  cout << "A_is_square=" << A_is_square << endl;

  // === simple check if A is symetrical === // 
  // voir avec les autres.
  bool A_is_maybe_sym = (void*)pUh == (void*)pVh; 

  // VF == true => VF type of Matrix
  bool VF=isVF(b->largs);    //=== used to set the solver ??? block matrix ??? ===/

  // set parameteer of the matrix :: 
  Data_Sparse_Solver ds;
  ds.factorize=0;
  ds.initmat=true;
  int np = OpCall_FormBilinear_np::n_name_param - NB_NAME_PARM_HMAT;
  SetEnd_Data_Sparse_Solver<R>(stack,ds, b->nargs,np);

  // set ds.sym = 0 
  ds.sym = 0;
  if(verbosity)
    cout << " we consider the block matrix as a non symetric matrix " << endl; 

  // J'ai repris ce qu'il y avait. 
  // PAC(e)     :: Attention peut être pas compatible avec les matrices bloques.
  // A repenser :: surtout pour le parametre symetrique? on le met ce parametre à zéro pour l'instant.
  // set ds.sym = 0 

  ds.sym = 0;
  if(verbosity)
    cout << " === we consider the block matrix as a non symetric matrix === (to be change in the future)" << endl; 

  if (! A_is_square )
   {
     if(verbosity>3) cout << " -- the solver  is un set  on rectangular matrix  " << endl;
    }

  // A quoi cela correspond?? Gestion du stack + autre
  WhereStackOfPtr2Free(stack)=new StackOfPtr2Free(stack);// FH aout 2007

  Matrice_Creuse<R> & A( * GetAny<Matrice_Creuse<R>*>((*a)(stack)));
  if(init) A.init();                                            //
  cout << " A.N=" <<  A.N() << endl;
  cout << " A.M=" <<  A.M() << endl;
  if( ! pUh || ! pVh) return SetAny<Matrice_Creuse<R>  *>(&A);  //

  // need to define the size of the entire matrix here
  A.resize( (*pVh)->totalNbOfDF(), (*pUh)->totalNbOfDF() ); 
  // test function (Vh) are the line
  // inconnu function (Uh) are the column

  
  // Assemble the variationnal form
  int maxJVh=NpVh;
  
  int offsetMatrixUh = 0;
  // loop over the block
  for( int i=0; i<NpUh; i++){
    int offsetMatrixVh = 0;
    if( ds.sym > 0 ){ maxJVh=(i+1); ffassert(maxJVh<NpVh);}
    for( int j=0; j<maxJVh; j++){
      cout << "offsetMatrixUh= " << offsetMatrixUh << ", offsetMatrixVh= " << offsetMatrixVh << endl;
      
      // construction du block (i,j)
      const list<C_F0> & b_largs=block_largs(i,j); 

      //const void * PUh = (void *) (*pUh)->vect[i]->getpVh();
      //const void * PVh = (void *) (*pVh)->vect[j]->getpVh();

      // size of the block
      int N_block = UhNbOfDf[i];
      int M_block = VhNbOfDf[j];
      
      Matrice_Creuse<R> *CCC = new Matrice_Creuse<R>() ;
      CCC->resize(M_block,N_block); // test function (Vh) are the line and inconnu function (Uh) are the column
      cout << "block:  i=" << i << "j=" << j <<  " (N,M)=" << M_block << " " << N_block << endl;
      //cout << "CCC=" << CCC<< endl;
      Matrice_Creuse<R> &BBB(*CCC); 
      //cout << "BBB=" << BBB<< endl;       
      // Faire les cast ici + tester ce soir ou demain

      // cas ::  Mesh, v_fes, v_fes
      if( (*pUh)->typeFE[i] == 2 && (*pVh)->typeFE[j] == 2 ){

        // ==== FESpace 2d : inconnue et test  ===
        const FESpace * PUh = (FESpace *) (*pUh)->vect[i]->getpVh();
        const FESpace * PVh = (FESpace *) (*pVh)->vect[j]->getpVh();
        creationBlockOfMatrixToBilinearForm< R, Mesh, FESpace,FESpace>( PUh, PVh, ds.sym, ds.tgv, b_largs, stack, BBB);
      }
      // cas ::  Mesh3, v_fes3, v_fes3 
      else if( (*pUh)->typeFE[i] == 3 && (*pVh)->typeFE[j] == 3 ){

        // ==== FESpace 3d : inconnue et test ===
        const FESpace3 * PUh = (FESpace3 *) (*pUh)->vect[i]->getpVh();
        const FESpace3 * PVh = (FESpace3 *) (*pVh)->vect[j]->getpVh();
        creationBlockOfMatrixToBilinearForm<R,Mesh3,FESpace3,FESpace3>( PUh, PVh, ds.sym, ds.tgv, b_largs, stack, BBB);
      }
      // cas :: MeshS, v_fesS, v_fesS 
      else if( (*pUh)->typeFE[i] == 4 && (*pVh)->typeFE[j] == 4 ){

        // ==== FESpace 3d Surf: inconnue et test ===
        const FESpaceS * PUh = (FESpaceS *) (*pUh)->vect[i]->getpVh();
        const FESpaceS * PVh = (FESpaceS *) (*pVh)->vect[j]->getpVh();
        creationBlockOfMatrixToBilinearForm<R,MeshS,FESpaceS,FESpaceS>( PUh, PVh, ds.sym, ds.tgv, b_largs, stack, BBB);
      }
      // cas :: MeshL, v_fesL, v_fesL
      else if( (*pUh)->typeFE[i] == 5 && (*pVh)->typeFE[j] == 5 ){

        // ==== FESpace 3d Curve: inconnue et test ===
        const FESpaceL * PUh = (FESpaceL *) (*pUh)->vect[i]->getpVh();
        const FESpaceL * PVh = (FESpaceL *) (*pVh)->vect[j]->getpVh();
        creationBlockOfMatrixToBilinearForm<R,MeshL,FESpaceL,FESpaceL>( PUh, PVh, ds.sym, ds.tgv, b_largs, stack, BBB);
      }
      // cas :: MeshL, v_fesL, v_fes
      else if( (*pUh)->typeFE[i] == 5 && (*pVh)->typeFE[j] == 2 ){

        // ==== FESpace 3d Curve: inconnue et 2d : test ===
        const FESpaceL * PUh = (FESpaceL *) (*pUh)->vect[i]->getpVh();
        const FESpace * PVh = (FESpace *) (*pVh)->vect[j]->getpVh();
        creationBlockOfMatrixToBilinearForm<R,MeshL,FESpaceL,FESpace>( PUh, PVh, ds.sym, ds.tgv, b_largs, stack, BBB);
      }
      // cas :: MeshL, v_fes, v_fesL
      else if( (*pUh)->typeFE[i] == 2 && (*pVh)->typeFE[j] == 5 ){

        // ==== FESpace 2d: inconnue et 3d Curve: test ===
        const FESpace * PUh = (FESpace *) (*pUh)->vect[i]->getpVh();
        const FESpaceL * PVh = (FESpaceL *) (*pVh)->vect[j]->getpVh();
        creationBlockOfMatrixToBilinearForm<R,MeshL,FESpace,FESpaceL>( PUh, PVh, ds.sym, ds.tgv, b_largs, stack, BBB);
      }
      // cas :: new OpMatrixtoBilinearForm< double, MeshS, v_fesS, v_fes3 >,      // 3D Surf / 3D volume on meshS
      else if( (*pUh)->typeFE[i] == 4 && (*pVh)->typeFE[j] == 3 ){

        // ==== FESpace 3d Surf: inconnue et 3d : test ===
        const FESpaceS * PUh = (FESpaceS *) (*pUh)->vect[i]->getpVh();
        const FESpace3 * PVh = (FESpace3 *) (*pVh)->vect[j]->getpVh();
        creationBlockOfMatrixToBilinearForm<R,MeshS,FESpaceS,FESpace3>( PUh, PVh, ds.sym, ds.tgv, b_largs, stack, BBB);
      } 
      // cas :: new OpMatrixtoBilinearForm< double, MeshS, v_fes3, v_fesS >,     // 3D volume / 3D Surf on meshS
      else if( (*pUh)->typeFE[i] == 3 && (*pVh)->typeFE[j] == 4 ){

        // ==== FESpace 3d : inconnue et 3d Surf : test ===
        const FESpace3 * PUh = (FESpace3 *) (*pUh)->vect[i]->getpVh();
        const FESpaceS * PVh = (FESpaceS *) (*pVh)->vect[j]->getpVh();
        creationBlockOfMatrixToBilinearForm<R,MeshS,FESpace3,FESpaceS>( PUh, PVh, ds.sym, ds.tgv, b_largs, stack, BBB);
      } 
      // cas :: new OpMatrixtoBilinearForm< double, MeshL, v_fesL, v_fesS >,       // 3D curve / 3D Surf on meshL
      else if( (*pUh)->typeFE[i] == 5 && (*pVh)->typeFE[j] == 4 ){

        // ====  FESpace 3d Curve : inconnue et 3d Surf : test ===
        const FESpaceL * PUh = (FESpaceL *) (*pUh)->vect[i]->getpVh();
        const FESpaceS * PVh = (FESpaceS *) (*pVh)->vect[j]->getpVh();
        creationBlockOfMatrixToBilinearForm<R,MeshL,FESpaceL,FESpaceS>( PUh, PVh, ds.sym, ds.tgv, b_largs, stack, BBB);
      }
      // cas :: new OpMatrixtoBilinearForm< double, MeshL, v_fesS, v_fesL >);       // 3D Surf / 3D curve on meshL
        else if( (*pUh)->typeFE[i] == 4 && (*pVh)->typeFE[j] == 5 ){

        // ====  FESpace 3d Surf : inconnue et 3d Curve : test ===
        const FESpaceS * PUh = (FESpaceS *) (*pUh)->vect[i]->getpVh();
        const FESpaceL * PVh = (FESpaceL *) (*pVh)->vect[j]->getpVh();
        creationBlockOfMatrixToBilinearForm<R,MeshL,FESpaceS,FESpaceL>( PUh, PVh, ds.sym, ds.tgv, b_largs, stack, BBB);
      }
      else{
        cerr << " =: Pas prise en compte des FESpace inconnue de type := "<< typeFEtoString( (*pUh)->typeFE[i] ) << endl;
        cerr << " =:                 avec des FESpace test de type    := "<< typeFEtoString( (*pVh)->typeFE[j] ) << endl;
        ffassert(0);
      }
      //cout << "BBB=" << BBB<< endl;
      A.pHM()->Add( BBB.pHM(), R(1), false, offsetMatrixVh, offsetMatrixUh ); // test function (Vh) are the line and inconnu function (Uh) are the column
      //cout << "BBB=" << BBB<< endl;
      //cout << "A=" << A << endl;
      delete CCC;
    offsetMatrixVh += VhNbOfDf[j];
    }
    offsetMatrixUh += UhNbOfDf[i];
  }
  
  A.pHM()->half = ds.sym;
  if (A_is_square)
    SetSolver(stack,VF,*A.A,ds);

  // === re-szet the original value of number of the component in the 
  //          - BilinearForm 
  //          - BC_set. 
  // ===

  // bilinear Form

  
  // Info necessaire :: " block_largs, beginBlockUh, beginBlockVh, NpUh, NpVh  

  // put the right number of component of each block
  for( int i=0; i<NpUh; i++){
      for( int j=0; j<NpVh; j++){
        
        const list<C_F0> *b_largs=&block_largs(i,j); 
        list<C_F0>::const_iterator b_ii,b_ib=b_largs->begin(),b_ie=b_largs->end(); 
        for (b_ii=b_ib;b_ii != b_ie;b_ii++){
          Expression e=b_ii->LeftValue();
          aType r = b_ii->left();

          // bilinear case
          if (r==atype<const  FormBilinear *>() ){
            const FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
      
            BilinearOperator * Op=const_cast<  BilinearOperator *>(bb->b);
            if (Op == NULL) {
              if(mpirank == 0) cout << "dynamic_cast error" << endl; 
            ffassert(0);
            }
      
            size_t Opsize= Op->v.size();
      
            KN<size_t> index_operator_finc(Opsize);
            KN<int>    new_index_funct_finc(Opsize);

            KN<size_t> index_operator_ftest(Opsize);
            KN<int>    new_index_funct_ftest(Opsize);

            // index to check if a integral is defined on multi block
            for(size_t jj=0; jj<Opsize; jj++){
              // attention la fonction test donne la ligne
              //  et la fonction test est en second
              BilinearOperator::K ll = Op->v[jj];
              pair<int,int> finc(ll.first.first), ftest(ll.first.second);

              long jj2= jj;

              index_operator_finc[ jj2] = jj;
              new_index_funct_finc[ jj2] = beginBlockUh[i]+finc.first;
            
              index_operator_ftest[ jj2]  = jj;
              new_index_funct_ftest[ jj2] = beginBlockVh[i]+ftest.first; 
            }
            changeIndexFunctionInconnue(*Op, index_operator_finc, new_index_funct_finc );
      
            changeIndexFunctionTest(*Op, index_operator_ftest, new_index_funct_ftest  );          
          }  
      
#ifndef FFLANG
#ifdef PARALLELE
          // ******************************************
          // Case BemKFormBilinear (KERNEL FORM ONLY)
          // ******************************************
          else if (r==atype<const BemFormBilinear *>() ){
            BemFormBilinear * bbtmp= dynamic_cast< BemFormBilinear *>(e);
            int VVFBEM = bbtmp->type;

            if(VVFBEM ==1){
              BemKFormBilinear * bb=new BemKFormBilinear(*dynamic_cast<const BemKFormBilinear *>(e));
              FoperatorKBEM * b=const_cast<  FoperatorKBEM *>(bb->b);
              if (b == NULL) { if(mpirank == 0) cout << "dynamic_cast error" << endl; exit(0);}

              // loop over the index of finconnue
              LOperaG * OpG = const_cast<LOperaG *>(b->fi);
              ffassert( OpG->v.size() == 1);
              size_t Opsize= OpG->v.size();

              for(size_t jjj=0; jjj<Opsize; jjj++){
                LOperaG::K *lf=&(OpG->v[jjj]);
                OpG->v[jjj].first.first += beginBlockUh[i];
                
                //pair<int,int> finc(lf->first);
                //cout << " new value :: block i,j=" << i << ","<< j << ", operateur jj= " << jjj << endl;
                //cout << " BemormBilinear: number of unknown finc = " << finc.first << endl;
                //cout << " BemFormBilinear: operator order   finc = " << finc.second << endl; 
                
              }


              // Loop over the index of ftest
              LOperaD * OpD = const_cast<LOperaD *>(b->ft);
              ffassert( OpD->v.size() == 1);
              Opsize= OpD->v.size();
              for(size_t jjj=0; jjj<Opsize; jjj++){
                LOperaD::K *lf=&(OpD->v[jjj]);
                OpD->v[jjj].first.first += beginBlockVh[j]; 
                
                //pair<int,int> finc(lf->first);
                //cout << " new value :: block i,j=" << i << ","<< j << ", operateur jj= " << jjj << endl;
                //cout << " BemormBilinear: number of unknown ftest = " << ftest.first << endl;
                //cout << " BemFormBilinear: operator order   ftest = " << ftest.second << endl;
                
              }
            }
            else if(VVFBEM == 2){
              BemPFormBilinear * bb=new BemPFormBilinear(*dynamic_cast<const BemPFormBilinear *>(e));
              FoperatorPBEM * b=const_cast<  FoperatorPBEM *>(bb->b);
              if (b == NULL) { if(mpirank == 0) cout << "dynamic_cast error" << endl; }

              cerr << " BEM Potential in composite FESpace in construction " << endl;
              ffassert(0);
            }

          }
#endif
#endif  
          // BC_set
          // case BC_set 
          else if(r == atype<const  BC_set  *>()){
            ffassert( i == j ); // diagonal block 
            BC_set * bc=dynamic_cast<BC_set *>(e);
        
            //KN<int>  new_index_funct_finc( bc.size() );
            int kk=bc->bc.size();
            for (int k=0;k<kk;k++)
            {
              //bc->bc[k].first += beginBlockUh[i];
              pair<int,Expression> &xx=bc->bc[k];
              xx.first += beginBlockUh[i];
            }
          }
        }
      }
  }


  cout <<"========================================================" << endl;
  cout <<"=                                                      =" << endl;
  cout <<"= check the component of Bilinear Form                 =" << endl;
  cout <<"=                                                      =" << endl;
  listOfComponentBilinearForm(largs);
  cout <<"=                                                      =" << endl;
  cout <<"========================================================" << endl;
  
  return SetAny<Matrice_Creuse<R>  *>(&A);
}

*/