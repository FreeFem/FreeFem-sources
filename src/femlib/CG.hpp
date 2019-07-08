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
// AUTHORS : ...
// E-MAIL  : ...

#ifndef _CG_HPP_
#define _CG_HPP_

#include <algorithm>
#include <ctime>

extern long verbosity;

template<class TypeIndex=int, class TypeScalar=double>
struct CGMatVirt {
public:
  typedef TypeIndex I;
  typedef TypeScalar R;

  I n, m;
  mutable int it;
  mutable double cpu;
  virtual R *addmatmul(R *x, R *Ax) const = 0;
  virtual ~CGMatVirt() {
    if (verbosity > 5)
    std::cout << " cpu CGMatVirt " << cpu << " s / " << it << " nb mul " << std::endl;
  }

  virtual R *matmul(R *x, R *Ax) const {
    it++;
    std::fill(Ax, Ax+n, 0.);
    double t0 = ((double) clock())/CLOCKS_PER_SEC;
    R*p = addmatmul(x, Ax);
    cpu += ((double)clock())/CLOCKS_PER_SEC - t0;
    return p;
  }
  virtual void SetInitWithBC(R*rhs, R *x) const {} // do nothing by default ..
  CGMatVirt(int nn, int mm=-1) : n(nn), m(mm < 0 ? nn : mm), cpu(0.), it(0) {}
  virtual int *pwcl() const { return 0; } // array know if node with BC (TGV)
};

template<class TypeIndex=int, class TypeScalar=double>
int ConjugueGradient(CGMatVirt<TypeIndex, TypeScalar> &A, // fonction et pointeur data pour A
                     CGMatVirt<TypeIndex,TypeScalar> &C, // fonction et pointeur data pour C
                     TypeScalar *b, // second membre
                     TypeScalar *x, // solution qui contient une initialisation
                     int nbitermax,
                     double eps,
                     int niveauimpression);

template<typename K, typename Z>
bool fgmres(CGMatVirt<Z, K> &A, // fonction et pointeur data pour A
            CGMatVirt<Z, K> &C, int leftC,
            K *y,
            K *x,
            double tol,
            int maxits,
            int restart=50,
            int verbo=3,
            int *perm=0);

template<class I, class K>
K *myscopy(I n, const K *x, K *y);
template<class I, class K>
K *myscal(I n, K a, K *x);
template<class TypeIndex=int, class TypeScalar=double>
inline double *ProduitMatVec(const CGMatVirt<TypeIndex, TypeScalar> *A, TypeScalar *x, TypeScalar *Ax) { return A->matmul(x, Ax); }
template<class TypeIndex=int,class TypeScalar=double>
inline double *ProduitMatVec(const CGMatVirt<TypeIndex, TypeScalar> &A, TypeScalar *x, TypeScalar *Ax) { return A.matmul(x, Ax); }
template<class I, class K> K *mysaxpy(I n, K a, const K *x, K *y);

template<class Z=int, class R=double>
struct CGMatVirtId : public CGMatVirt<Z, R> {
  CGMatVirtId(Z nn): CGMatVirt<Z, R> (nn, nn) {}
  R *matmul(R *x, R *Ax) const { myscopy(this->n, x, Ax); return Ax; }
  R *addmatmul(R *x, R *Ax) const { mysaxpy(this->n,R(1.), x, Ax); return Ax; }

};

#endif // _CG_HPP_
