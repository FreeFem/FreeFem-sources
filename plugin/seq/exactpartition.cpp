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
// SUMMARY : Remove the upper part of a CSR if the supplied matrix is not symmetric
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Pierre Jolivet
// E-MAIL  : pierre.jolivet@enseeiht.fr

// *INDENT-OFF* //
// ff-c++-LIBRARY-dep:
// ff-c++-cpp-dep:
// *INDENT-ON* //

#include "ff++.hpp"
long exactpartition(int n, int m, double **v, long *J) {
  int kkerr = 0;
  long N = (1 << 30);

  cout << " exactpartition " << n << " " << m << " N =" << N << endl;

  for (int i = 0; i < m; ++i) {
    long s = 0, j0 = N;

    for (int j = 0; j < n; ++j) {
      if (v[j]) {
        int jj = J[j];
        double w = v[j][i];
        long lw = lrint(w * N);
        if (lw && jj < j0) {
          j0 = j;    // prco le plus petit non nulle
        }

        v[j][i] = (double)lw;
        s += lw;
      }
    }

    ffassert(s && j0 < N);
    double ss = 0;

    for (int j = 0; j < n; ++j) {
      if (v[j]) {
        double w = v[j][i] / s;
        long lw = lrint(w * N);
        double we = (double)lw / (double)N;
        v[j][i] = we;
        ss += we;
      }
    }

    double err = ss - 1.;
    assert(fabs(err) * N < 10);    // bofbof ???? FH.
    v[j0][i] -= err;
    // verif ..
    ss = 0;

    for (int j = 0; j < n; ++j) {
      if (v[j]) {
        ss += v[j][i];
      }
    }

    kkerr += (ss != 1.);
  }

  ffassert(kkerr == 0);
  return 0;
}

long exactpartition(FEbaseArrayKn< double > *const &p, KN< long > *const &pj) {
  int n = p->N, m = 0;
  double **v = new double *[n];
  int kerr = 0;

  for (int i = 0; i < n; ++i) {
    KN< double > *vi = p->get(i);
    int mi = vi ? vi->N( ) : 0;
    if (m == 0) {
      m = mi;
    } else if (m != mi) {
      kerr++;
    }

    // v[i] = *vi;
    v[i] = vi ? *vi : NULL;
  }

  ffassert(kerr == 0);
  ffassert(pj->N( ) >= n);
  exactpartition(n, m, v, *pj);
  delete[] v;
  return 0;
}

long exactpartition(KN< KN< double > > *const &p, KN< long > *const &pj) {
  int n = p->N( ), m = 0;
  double **v = new double *[n];
  int kerr = 0;

  for (int i = 0; i < n; ++i) {
    KN_< double > vi = (*p)(i);
    int mi = vi ? vi.N( ) : 0;
    if (m == 0) {
      m = mi;
    } else if (m != mi) {
      kerr++;
    }

    if (mi == 0) {
      v[i] = 0;
    } else {
      v[i] = vi;
    }
  }

  ffassert(kerr == 0);
  ffassert(pj->N( ) >= n);
  exactpartition(n, m, v, *pj);
  delete[] v;
  return 0;
}

static void Load_Init( ) {
  // to be sure  to have unique  add
  if (!Global.Find("exactpartition").NotNull( )) {
    Global.Add("exactpartition", "(",
               new OneOperator2_< long, FEbaseArrayKn< double > *, KN< long > * >(exactpartition));
    // KN<KN<double> >
    Global.Add("exactpartition", "(",
               new OneOperator2_< long, KN< KN< double > > *, KN< long > * >(exactpartition));
  }
}

LOADFUNC(Load_Init)
