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
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include "ff++.hpp"
#include "AFunction_ext.hpp"
#include <vector>


long ichol(MatriceMorse<R> & A,MatriceMorse<R> &  L,double tgv)
{
    // cf https://en.wikipedia.org/wiki/Incomplete_Cholesky_factorization
    cout << " tgv " << tgv << endl;
    ffassert( A.half && L.half);
    ffassert( A.n == L.n);
    int n =A.n,i,j,k,kk;
    double tgv1= sqrt(tgv);
    double tgve =tgv*0.99999999;
    if(tgve < 1) tgve=1e200;
    double nan=sqrt(-1.);
    for(int k=0; k<L.nnz;++k)
        L.aij[k]=nan;
    int BC = 0;
    long err=0;
    A.CSR();
    L.CSR();
    for(int i=0; i< n; ++i)
    {
        int ai1=A.p[i+1]-1;
        int ai0=A.p[i];
        int li1=L.p[i+1]-1;
        int li0=L.p[i];
        double  Aii=A.aij[ai1];
        if (Aii > tgve)
        { // B.C
            for (kk=li0;kk<li1;kk++)
                L.aij[kk]=0; //  remove row and col
            L.aij[li1]=1;
            BC++;
        }
        else
        {
            for (kk=li0;kk<li1;kk++) //  Build Lij is existe j<i
            {
                int j = L.j[kk]; // j < i
                ffassert(j<i);
                int lj1=L.p[j+1]-1;
                int lj0=L.p[j];
                
                double *pAij = A.pij(i,j) ;
                double Lij = pAij ? *pAij: 0.,Aij=Lij;
                for(int kkk= lj0; kkk<lj1; ++kkk)// loop  row j
                {  //cout << " ?? " << kkk << " "<< lj0 << " " <<lj1 <<endl;
                    int k = L.j[kkk];
                    //cout << " @@@" << i << " " << j << " ( " << lj0 << " " << lj1 << ")  " << k << " // " << kkk <<  " " << endl;
                    ffassert(k >=0 && k < j);
                    double Ljk = L.aij[kkk], *pLik=L.pij(i,k), Lik = pLik ? *pLik : 0.;
                    //cout << " *** " << k << " " <<Lik *Ljk <<endl;
                    Lij -= Lik *Ljk;
                }
                Lij /=  L(j,j);
                L.aij[kk] =Lij;
                //cout <<kk << " " << j << " " << Lij << " "<< Aij << " , ";
                
            }
            // cout << " **" << endl;
            for(int k= li0; k<li1; ++k)
                Aii -= L.aij[k]*L.aij[k];
            if( Aii <=1e-30) {
                if(err<10 && verbosity )
                    cout << "   ichol neg pivot:"<< i << " " << Aii << " " << A.aij[ai1] << endl;
                Aii=1;// Bof Bof !!!
                err++;
            }
            double Lii = sqrt(Aii);
             L.aij[li1] =Lii;
            
        }
    }
    if(verbosity>2) cout << "  -- ichol:  N BC = " << BC << " nberr "<< err <<endl;
    return err;
}

inline R pscal(R*L,int *cl,int kl,int kl1,int i, MatriceMorse<R> &  Ut,int j )
{
    int ku = Ut.p[j],ku1=Ut.p[j]-1;
    int k= min(i,j); //  common  part
    R r =0;
    //cout << " pscal: "<<  i << " " << j << "  min: " << k << endl;
    for(int l=kl;l<kl1;++l)
    {
        
        int jl = cl[l];
        //cout << "     ##" <<i << " " << jl << " " << (jl > k) << endl;
        if( jl >= k) break;
        R Lijl = L[l];
        R * pUtjjl = Ut.pij(j,jl);
        if(pUtjjl) { r += Lijl* *pUtjjl;
            //cout <<   "   **  "<< Lijl << " " << *pUtjjl << " " << jl <<  " " << r << endl;
        }
    }
    ffassert ( r==r );
    return r;
}
long iLU(MatriceMorse<R> & A,MatriceMorse<R> &  L,MatriceMorse<R> &  Ut,double tgv)
{
    /*  Algo LU :
 L = L + I, U = U+D
 for(int i=0;i<n; ++i)
 {
 for(int j=0;j<i;++j) L(i,j) = (A(i,j) - (L(i,':'),U(':',j)))/ D(j,j);
 for(int j=0;j<i;++j) U(j,i) = (A(j,i) - (L(j,':'),U(':',i))) ;
 D(i,i) = A(i,i) - (L(i,':'),U(':',i));
 }

 */
      A.CSR();
      L.CSR();
      Ut.dotranspose();
      Ut.CSR();
    if( verbosity > 4)
        cout << "   - ILU fact:   tgv " << tgv << endl;
    ffassert( A.n == L.n);
    ffassert( A.n == Ut.n);
    int sym = A.half    ;
    int n =A.n,i,j,k,kk;
    double tgve =tgv*0.999;
    if(tgve < 1) tgve=1e200;
    double NaN=sqrt(-1.);
    fill(L.aij,L.aij+L.nnz,NaN);
    fill(Ut.aij,Ut.aij+Ut.nnz,NaN);
    int BC = 0;
    KN<int> wbc(n);
    long err=0;
    double mUii=1e200;
    for(int i=0; i< n; ++i)
    {
        int ai1=A.p[i+1]-1;
        int ai0=A.p[i];
        int li1=L.p[i+1]-1;
        int li0=L.p[i];
        int ui1=Ut.p[i+1]-1;
        int ui0=Ut.p[i];
        err += Ut.j[ui1] != i;
        err += L.j[li1] != i;
        ffassert( L.j[li1]==i && Ut.j[ui1]==i);
        double  Aii=A(i,i),Uii;
        
        int BCi ;
         wbc[i]  = BCi= (Aii > tgve);
        if (BCi)
        { // B.C
            fill(L.aij+li0,L.aij+li1,0.);
            fill(Ut.aij+ui0,Ut.aij+ui1,0.);
            L.aij[li1]=1.;
            Ut.aij[ui1]=Aii;
            BC++;
        }
       else
        {
            for(int l=li0;l<li1;++l) // coef of  L non zero
            {
                int j   = L.j[l];
                R *pAij = A.pij(i,j), Aij = pAij ? *pAij : 0.;
                
                R Ujj=Ut(j,j);
                ffassert(j<i);
                L.aij[l] = (Aij - pscal(L.aij,L.j,li0,li1,i, Ut,j)) / Ujj;
            }
             for(int u=ui0;u<ui1;++u) // coef of  Ut  non zero
            {
                int j   = Ut.j[u];// Ut(j,i) == U(j,i)
                R *pAji = sym ? A.pij(i,j) : A.pij(j,i), Aji = pAji ? *pAji : 0.;
                if( wbc[j] )  Aji=0; // remove row term  if BC. on j  ...
                ffassert(j<i);// transpose
                Ut.aij[u] = (Aji - pscal(Ut.aij,Ut.j,ui0,ui1,i, L,j));
            }
            Uii= Aii - pscal(Ut.aij,Ut.j,ui0,ui1,i, L,i);
            L(i,i) =1.;

            mUii = min(mUii,abs(Uii));
            
            if (abs(Uii)< 1e-30)
            {
                if(verbosity && err<10)
                    cerr << "    error: ILU nul pivot "<< i << " " << Uii << endl;
                Uii =1;
                err++;
            }
            Ut(i,i) = Uii;
        }
    }
    if(verbosity>2 || err ) cout << "   - ILU: Nb BC = " << BC << "nb err =" << err << " main Uii " << mUii << endl;
    Ut.dotranspose();
    Ut.CSC();
    L.CSR();
    return  err;
}


long ff_ilu (Matrice_Creuse<R> * const & pcA,Matrice_Creuse<R> * const & pcL,Matrice_Creuse<R> * const & pcU,double const & tgv)
{
    MatriceCreuse<R> * pa=pcA->A;
    MatriceCreuse<R> * pl=pcL->A;
    MatriceCreuse<R> * pu=pcU->A;
    ffassert( pa  && pl && pu);
    MatriceMorse<R> *pA= dynamic_cast<MatriceMorse<R>* > (pa);
    MatriceMorse<R> *pL = dynamic_cast<MatriceMorse<R>* > (pl);
    MatriceMorse<R> *pU = dynamic_cast<MatriceMorse<R>* > (pu);
    ffassert(pL && pA && pU);
    
    return iLU(*pA,*pL,*pU,tgv);
}

long ff_ichol (Matrice_Creuse<R> * const & pcA,Matrice_Creuse<R> * const & pcL,double const & tgv)
{
    MatriceCreuse<R> * pa=pcA->A;
    MatriceCreuse<R> * pl=pcL->A;
    ffassert( pa  && pl );
    MatriceMorse<R> *pA= dynamic_cast<MatriceMorse<R>* > (pa);
    MatriceMorse<R> *pL = dynamic_cast<MatriceMorse<R>* > (pl);
    ffassert(pL && pA);
    
    return ichol(*pA,*pL,tgv);
}
long ff_ilu (Matrice_Creuse<R> *  const & pcA,Matrice_Creuse<R> *  const & pcL,Matrice_Creuse<R> * const & pcU)
{
    return ff_ilu(pcA,pcL,pcU, ff_tgv);
}
long ff_ichol (Matrice_Creuse<R> *  pcA,Matrice_Creuse<R> *  pcL)
{
    return ff_ichol(pcA,pcL,ff_tgv);
}
void LU_solve(MatriceMorse<R> &T,int cas,KN<double> & b,bool trans)
{
    int n =T.n,i,j,k,k1,k0;
    if(cas<0) T.CSR();
    else if  (cas>0) T.CSC();
    int  *ij= cas<0? T.j : T.i;
    ffassert(cas!=0);
    ffassert( n == b.N());
    if( trans == (cas <0) ) // (cas <0 et trans) or (not et  cas >0)
    { // U = L'
        if(verbosity>9)
            cout  << " LU_solve:: Remonte:  " << cas << " " << trans << endl;
        for(int i=n-1; i>=0; --i)
        {
            k0 = T.p[i];
            k1 = T.p[i+1]-1;
            b[i] /= T.aij[k1];
            
            for (k=k0;k<k1;k++)
            {
                int j = ij[k];
                b[j] -= b[i]*T.aij[k];
            }
            
            assert(ij[k] == i);
        }
    }
    else // cas >0 et not trans
    {
        if(verbosity>9)
        cout  << " LU_solve:: Descente:  " << cas << " " << trans << endl;
        for(int i=0; i< n; ++i)
        {
            R bi= b[i];
            for (k=T.p[i];k<T.p[i+1]-1;k++)
            {
                int j = ij[k];
                bi -= b[j]*T.aij[k];
            }
            b[i] = bi/ T.aij[k];
            assert(ij[k] == i);
        }
        
    }
    
    
    
}
bool ff_ichol_solve(Matrice_Creuse<R> * pcL,KN<double> * b)
{
    // L L' u = b =>  L uu = b;  L' u = uu;
    MatriceCreuse<R> * pl=pcL->A;
    ffassert(pl );
    MatriceMorse<R> *pL = dynamic_cast<MatriceMorse<R>* > (pl);
    ffassert(pL );
    LU_solve(*pL,-1,*b,0);
    LU_solve(*pL,-1,*b,1);
    
    return true;
}
bool ff_ilu_solve(Matrice_Creuse<R> * const & pcL,Matrice_Creuse<R> *const &  pcU,KN<double> * const & b)
{
    // L Ut u = b =>  L uu = b;  L' u = uu;
    MatriceCreuse<R> * pl=pcL->A;
    ffassert(pl );
    MatriceMorse<R> *pL = dynamic_cast<MatriceMorse<R>* > (pl);
    ffassert(pL );
    MatriceCreuse<R> * pu=pcU->A;
    ffassert(pu );
    MatriceMorse<R> *pU = dynamic_cast<MatriceMorse<R>* > (pu);
    ffassert(pl );
    LU_solve(*pL,-1,*b,0);
    LU_solve(*pU, 1,*b,0);
    
    return true;
}


static void Load_Init () {
    cout << " load: init Incomplete Cholesky " << endl;
    Global.Add("ichol", "(", new OneOperator2<long,Matrice_Creuse<R> * ,Matrice_Creuse<R> * >(ff_ichol));
    Global.Add("ichol", "(", new OneOperator3_<long,Matrice_Creuse<R> * ,Matrice_Creuse<R> * ,double >(ff_ichol));
    Global.Add("iLU", "(", new OneOperator4_<long,Matrice_Creuse<R> * ,Matrice_Creuse<R> * ,Matrice_Creuse<R> *,double >(ff_ilu));
    Global.Add("iLU", "(", new OneOperator3_<long,Matrice_Creuse<R> * ,Matrice_Creuse<R> * ,Matrice_Creuse<R> * >(ff_ilu));
    Global.Add("iluSolve", "(", new OneOperator3_<bool ,Matrice_Creuse<R> * ,Matrice_Creuse<R> * , KN<R> *>(ff_ilu_solve));
    Global.Add("icholSolve", "(", new OneOperator2<bool ,Matrice_Creuse<R> * , KN<R> *>(ff_ichol_solve));


}

LOADFUNC(Load_Init)
