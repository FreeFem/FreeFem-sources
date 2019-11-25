/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Jacques Morice
// E-MAIL  : jacques.morice@ann.jussieu.fr

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

/*
 * Thank to the ARN ()  FF2A3 grant
 * ref:ANR-07-CIS7-002-01
 */

#include "ff++.hpp"
/*
// using namespace Fem2D;

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
                  if (verbosity) {cout << "  thresholding : remove " << nnzo-phm->nnz  << " them in
the matrix " << sparse_mat << " " << threshold << endl;} } else if (verbosity) {cout << " empty
matrix " << sparse_mat << endl;}
        }

        return t.v;
}

template<class T>
Thresholding<T> to_Thresholding (Matrice_Creuse<T> *v) {return Thresholding<T>(v);}
*/
static void Load_Init( ) {    // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
  if (mpirank == 0)
    cerr << " ++ WARNING Obsolete plugin: thresholding  is now in FreeFem++ core " << endl;
  /*
      typedef Thresholding<double> TMR;
      typedef Thresholding<Complex> TMC;
      typedef Matrice_Creuse<double> MR;
      typedef Matrice_Creuse<Complex> MC;
      Dcl_Type<TMR>();
      Dcl_Type<TMC>();
      //TMR t(0);
      //thresholding2(t, 0.);
      Add<MR *>("thresholding", ".", new OneOperator1<TMR, MR *>(to_Thresholding));
      Add<TMR>("(", "", new OneOperator2_<MR *, TMR, double>(thresholding2));
      Add<MC *>("thresholding", ".", new OneOperator1<TMC, MC *>(to_Thresholding));
      Add<TMC>("(", "", new OneOperator2_<MC *, TMC, double>(thresholding2));
  */
}

LOADFUNC(Load_Init)
