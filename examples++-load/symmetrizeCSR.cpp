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
// SUMMARY : remove the upper part of a CSR if the supplied matrix is not symmetric
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Pierre Jolivet
// E-MAIL  : pierre.jolivet@enseeiht.fr

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include "ff++.hpp"

template<class T>
long symmetrizeCSR (Matrice_Creuse<T> *const &A) {
	MatriceMorse<T> *mA = static_cast<MatriceMorse<T> *>(&(*A->A));
	if (!mA->symetrique) {
		mA->symetrique = true;
		std::vector<int> cl;
		std::vector<T> a;
		a.reserve(mA->nbcoef);
		cl.reserve(mA->nbcoef);
		unsigned int save = mA->lg[0];

		for (unsigned int i = 0; i < mA->n; ++i) {
			for (unsigned int j = save; j < mA->lg[i + 1]; ++j) {
				int col = mA->cl[j];
				if (col <= i) {
					T val = mA->a[j];
					if (abs(val) > 1e-14) {
						a.push_back(val);
						cl.push_back(col);
					}
				} else {
					break;
				}
			}

			save = mA->lg[i + 1];
			mA->lg[i + 1] = cl.size();
		}

		delete [] mA->cl;
		delete [] mA->a;
		int *col = new int[cl.size()];
		T *val = new T[cl.size()];

		for (unsigned int i = 0; i < cl.size(); ++i) {
			col[i] = cl[i];
			val[i] = a[i];
		}

		mA->cl = col;
		mA->a = val;
		mA->nbcoef = cl.size();
	}

	return 1L;
}

static void Load_Init () {
	Global.Add("symmetrizeCSR", "(", new OneOperator1_<long, Matrice_Creuse<double> *>(symmetrizeCSR<double> ));
	Global.Add("symmetrizeCSR", "(", new OneOperator1_<long, Matrice_Creuse<std::complex<double> > *>(symmetrizeCSR<std::complex<double> > ));
}

LOADFUNC(Load_Init)
