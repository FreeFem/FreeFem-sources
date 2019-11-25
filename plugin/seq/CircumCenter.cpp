#include "ff++.hpp"

/*!
 * \brief CircumCenter
 * \param pTh Fem2D::Mesh const *const &
 * \param pcx KN<double> * const &
 * \param pcy KN<double> * const
 * \return 0L
 */
long CircumCenter(Fem2D::Mesh const *const &pTh, KN< double > *const &pcx,
                  KN< double > *const &pcy) {
  KN< double > &cx = *pcx;
  KN< double > &cy = *pcy;
  const Mesh &Th = *pTh;

  ffassert(Th.nt == cx.N( ));
  ffassert(Th.nt == cy.N( ));
  for (int k = 0; k < Th.nt; ++k) {
    const Mesh::Element &K = Th[k];
    R2 A = K[0], B = K[1], C = K[2];
    R2 a = (B + C) / 2, aa(a + R2(B, C).perp( ));
    R2 b = (A + C) / 2, bb(b + R2(A, C).perp( ));
    // intersection ddroite [a,aa] and [b,bb]
    double ab = -det(a, aa, b);
    double abb = det(a, aa, bb);
    double s = ab + abb;
    R2 P = bb * ab / s + b * abb / s;
    cx[k] = P.x;
    cy[k] = P.y;
    //  verif
  }
  return 0L;
}

static void Load_Init( ) {    // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
  if (verbosity) {
    cout << " lood: CircumCenter  " << endl;
  }

  Global.Add("CircumCenter", "(",
             new OneOperator3_< long, pmesh, KN< double > *, KN< double > * >(CircumCenter));
}
LOADFUNC(Load_Init)
