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
// AUTHORS : ...
// E-MAIL  : ...


// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include "ff++.hpp"
#include <vector>

void ichol(MatriceMorse<R> & A,MatriceMorse<R> &  L,double tgv)
{
    // cf https://en.wikipedia.org/wiki/Incomplete_Cholesky_factorization
    cout << " tgv " << tgv << endl;
    ffassert( A.symetrique && L.symetrique);
    ffassert( A.n == L.n);
    int n =A.n,i,j,k,kk;
    double tgve =tgv*0.99999999;
    if(tgve < 1) tgve=1e200;
    double nan=sqrt(-1.);
    for(int k=0; k<L.nbcoef;++k)
        L.a[k]=nan;
    int BC = 0;
    for(int i=0; i< n; ++i)
    {
        int ai1=A.lg[i+1]-1;
        int ai0=A.lg[i];
        int li1=L.lg[i+1]-1;
        int li0=L.lg[i];
        double  Aii=A.a[ai1];
        if (Aii > tgve)
        { // B.C
            for (kk=li0;kk<li1;kk++)
                L.a[kk]=0; //  remove row and col
            L.a[li1]=1.;
            BC++;
        }
        else
        {
            for (kk=li0;kk<li1;kk++) //  Build Lij is existe j<i
            {
                int j = L.cl[kk]; // j < i
                ffassert(j<i);
                int lj1=L.lg[j+1]-1;
                int lj0=L.lg[j];
                
                double *pAij = A.pij(i,j) ;
                double Lij = pAij ? *pAij: 0.,Aij=Lij;
                for(int kkk= lj0; kkk<lj1; ++kkk)// loop  row j
                {  //cout << " ?? " << kkk << " "<< lj0 << " " <<lj1 <<endl;
                    int k = L.cl[kkk];
                    //cout << " @@@" << i << " " << j << " ( " << lj0 << " " << lj1 << ")  " << k << " // " << kkk <<  " " << endl;
                    ffassert(k >=0 && k < j);
                    double Ljk = L.a[kkk], *pLik=L.pij(i,k), Lik = pLik ? *pLik : 0.;
                    //cout << " *** " << k << " " <<Lik *Ljk <<endl;
                    Lij -= Lik *Ljk;
                }
                Lij /=  L(j,j);
                L.a[kk] =Lij;
                //cout <<kk << " " << j << " " << Lij << " "<< Aij << " , ";
                
            }
            // cout << " **" << endl;
            for(int k= li0; k<li1; ++k)
                Aii -= L.a[k]*L.a[k];
            if( Aii <=0.) { cout << i << " " << Aii << " " << A.a[ai1] << endl; Aii=1;}
            double Lii = sqrt(Aii);
            //  cout << " L " << i << " " << Lii << " " << A.a[ai1]  << endl;
            L.a[li1] =Lii;
            
        }
    }
    cout << " N BC = " << BC <<endl;
}



bool ff_ichol (Matrice_Creuse<R> * const & pcA,Matrice_Creuse<R> * const & pcL,double const & tgv)
{
    MatriceCreuse<R> * pa=pcA->A;
    MatriceCreuse<R> * pl=pcL->A;
    ffassert( pa  && pl );
    MatriceMorse<R> *pA= dynamic_cast<MatriceMorse<R>* > (pa);
    MatriceMorse<R> *pL = dynamic_cast<MatriceMorse<R>* > (pl);
    ffassert(pL && pA);
    ichol(*pA,*pL,tgv);
    return true;
}
bool ff_ichol (Matrice_Creuse<R> *  pcA,Matrice_Creuse<R> *  pcL)
{
    return ff_ichol(pcA,pcL,ff_tgv);
}
void ichol_solve(MatriceMorse<R> &L,KN<double> & b,bool trans)
{
    int n =L.n,i,j,k,k1,k0;
    ffassert(L.symetrique);
    ffassert( L.n == b.N());
    if(trans)
    {
        for(int i=n-1; i>=0; --i)
        {
            k0 = L.lg[i];
            k1 = L.lg[i+1]-1;
            b[i] /= L.a[k1];
            
            for (k=k0;k<k1;k++)
            {
                int j = L.cl[k];
                b[j] -= b[i]*L.a[k];
            }
            
            assert(L.cl[k] == i);
        }
    }
    else
    {
        for(int i=0; i< n; ++i)
        {
            R bi= b[i];
            for (k=L.lg[i];k<L.lg[i+1]-1;k++)
            {
                int j = L.cl[k];
                bi -= b[j]*L.a[k];
            }
            b[i] = bi/ L.a[k];
            assert(L.cl[k] == i);
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
    ichol_solve(*pL,*b,0);
    ichol_solve(*pL,*b,1);
    
    return true;
}


static void Load_Init () {
    cout << " lood: init Incomplete Cholesky " << endl;
    Global.Add("ichol", "(", new OneOperator2<bool,Matrice_Creuse<R> * ,Matrice_Creuse<R> * >(ff_ichol));
    Global.Add("ichol", "(", new OneOperator3_<bool,Matrice_Creuse<R> * ,Matrice_Creuse<R> * ,double >(ff_ichol));
    Global.Add("icholSolve", "(", new OneOperator2<bool ,Matrice_Creuse<R> * , KN<R> *>(ff_ichol_solve));
    
}

LOADFUNC(Load_Init)
