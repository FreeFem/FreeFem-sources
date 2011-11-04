// ORIG-DATE:   September 2010
// -*- Mode : c++ -*%
//
// SUMMARY  : seuillage des matrices EF de freefem++   
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
// AUTHOR   : Jacques Morice
// E-MAIL   : jacques.morice@ann.jussieu.fr
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

 Thank to the ARN ()  FF2A3 grant
 ref:ANR-07-CIS7-002-01 
 */
   
#include "ff++.hpp"

//using namespace Fem2D;

template<class T> struct  Thresholding{ 
  Matrice_Creuse<T> *v;
  Thresholding( Matrice_Creuse<T> * vv) : v(vv) {}
 }; 

template<class R> 
Matrice_Creuse<R> *thresholding2(const Thresholding<R> & t,const double &threshold)
 {  
     Matrice_Creuse<R> * sparse_mat =t.v;
     if(sparse_mat) 
       {
	 int n=sparse_mat->N(),m=sparse_mat->M();
	 map<pair<int,int>,R> M;
	 if (n >0 && m>0 && sparse_mat->A) 
	   {
	     int nrt = sparse_mat->A->NbCoef();
	     sparse_mat->A->addMatTo(R(1.),M,false,0,0,false,threshold);	    
	     (M)[make_pair(n-1,m-1)]+=R();
	     bool sym=false; // bof bof  
	     sparse_mat->typemat=TypeSolveMat(TypeSolveMat::GMRES); //  none square matrice (morse)
	     sparse_mat->A.master(new MatriceMorse<R>(n,m,M,sym));
	     nrt-=sparse_mat->A->NbCoef();
	     if(verbosity) cout << "  thresholding= remove " << nrt << " them in the matrix " << sparse_mat << " " << threshold << endl; 
	   }
	 else if(verbosity) cout << " empty matrix " <<sparse_mat << endl;
       }
  return  t.v;
 }

template<class T> 
Thresholding<T> to_Thresholding( Matrice_Creuse<T> *v){ return Thresholding<T>(v);}


class Init1 { public:
  Init1();
};
 
LOADINIT(Init1)  //  une variable globale qui serat construite  au chargement dynamique 

Init1::Init1()
{  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 
    typedef Thresholding<double>  TMR ;
    typedef Thresholding<Complex>  TMC ;
    typedef Matrice_Creuse<double>  MR ;
    typedef Matrice_Creuse<Complex>  MC ;
    Dcl_Type< TMR > ();
    Dcl_Type< TMC > ();
    TMR t(0);
    thresholding2(t,0.);
    Add<MR *>("thresholding",".",new OneOperator1< TMR ,MR *>(to_Thresholding));
    Add<TMR   >("(","",new OneOperator2_<MR  *, TMR , double >(thresholding2)); 
    Add<MC *>("thresholding",".",new OneOperator1< TMC ,MC *>(to_Thresholding));
    Add<TMC >("(","",new OneOperator2_<MC *, TMC, double >(thresholding2)); 
    
}
