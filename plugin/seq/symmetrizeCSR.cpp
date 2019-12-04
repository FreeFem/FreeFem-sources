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
// SUMMARY : remove the upper part of a CSR if the supplied matrix is not symmetric
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Pierre Jolivet
// E-MAIL  : pierre.jolivet@enseeiht.fr

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

#include "ff++.hpp"
/*
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
                  if (verbosity) {cout << "  symmetrizeCSR remove " << (long) nnzo-(long) phm->nnz
<< " them in the matrix " << sparse_mat << endl;} } else if (verbosity) {cout << " empty matrix " <<
sparse_mat << endl;}
        }

        return 1L;
}
*/
static void Load_Init( ) {
  if (mpirank == 0)
    cerr << " ++ WARNING Obsolete plugin: symmetrizeCSR  is now in FreeFem++ core " << endl;
  /*
          Global.Add("symmetrizeCSR", "(", new OneOperator1_<long, Matrice_Creuse<double>
     *>(symmetrizeCSR<double> )); Global.Add("symmetrizeCSR", "(", new OneOperator1_<long,
     Matrice_Creuse<std::complex<double> > *>(symmetrizeCSR<std::complex<double> > ));
   */
}

LOADFUNC(Load_Init)
