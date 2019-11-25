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

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <cmath>
using namespace std;
#include "RNM.hpp"
#include "bmo.hpp"

int irand_(int i) {
#ifdef _WIN32
  srand(i);
  return rand( );
#else
  srandom(i);
  return random( );
#endif
}

static /* Subroutine */ double xrandme(integer ii) {
#ifdef _WIN32
  static double xrd2, xrd3, xrd4, xrd5, xrd6;
  /* System generated locals */
  int i__1, i__2, i__3, i__4;
  double d__1;
  double xilim = 2147483648.;
  xrd2 = (d__1 = irand_(ii) / xilim, abs(d__1));
  i__1 = irand_(ii);
  xrd3 = (d__1 = irand_(i__1) / xilim, abs(d__1));
  i__2 = irand_(ii);
  i__1 = irand_(i__2);
  xrd4 = (d__1 = irand_(i__1) / xilim, abs(d__1));
  i__3 = irand_(ii);
  i__2 = irand_(i__3);
  i__1 = irand_(i__2);
  xrd5 = (d__1 = irand_(i__1) / xilim, abs(d__1));
  i__4 = irand_(ii);
  i__3 = irand_(i__4);
  i__2 = irand_(i__3);
  i__1 = irand_(i__2);
  xrd6 = (d__1 = irand_(i__1) / xilim, abs(d__1));
  double xx = (xrd2 + xrd3 + xrd4 + xrd5 + xrd6) / 5.;
#else
  // srandom(ii);
  long r = random( );
  // r = random();	//TODO why twice?
  const unsigned long rmax = (1UL << 31) - 1;
  double xx = (double)r / (double)rmax;
#endif
  // cout << " \t\t\trand = " << xx << " " << ii << endl;
  return xx;
}

static istream &Eat2LN(istream &f) {
  int c;

  while ((c = f.get( )) != '\n') {
    cout << char(c);
    assert(f.good( ));
  }

  cout << endl;
  return f;
}

/*!
 * \brief main
 * \param xx Vect &
 * \param xxmin Vect &
 * \param xxmax Vect &
 */
double BijanMO::main(Vect &xx, Vect &xxmin, Vect &xxmax) {
  /* Local variables */
  double costsave;
  integer irestart;
  double f;
  Vect v(ndim), v0(ndim), x1(ndim), hgc(ndim), fpx(ndim), fpx0(ndim), temp(ndim), xmin(ndim),
    xmax(ndim), xsave(ndim), vinit(ndim), xsave0(ndim);
  double rho;
  double rho0;
  double rho00;
  integer iter1;
  double gnorm;
  integer iterbvp, itersom;
  double irestart2;

  ncstr = 0;
  nbeval = 0;
  nbevalp = 0;
  // init ..
  ffassert(ndim == xx.N( ));
  ffassert(ndim == xxmin.N( ));
  ffassert(ndim == xxmax.N( ));
  vinit = xx;
  xmin = xxmin;
  xmax = xxmax;

  finit = func(vinit);
  if (debug) {
    cout << " ndim = " << ndim << endl;
  }

  f = finit;
  xopt1 = xoptg = vinit;
  if (finit < epsij) {
    cstropt = cstr;
    fseulopt = fseul;
    goto L9101;
  }

  epsij *= finit;
  if (debug) {
    cout << " F = " << finit << endl;
  }

  if (ncstr > 0) {
    if (debug) {
      cout << " CSTR = ";

      for (int i = 0; i < ncstr; ++i) {
        cout << cstr[i] << " ";
      }

      cout << endl;
    }
  }

  if (debug) {
    cout << finit << " " << 1. << " " << xoptg[0] << " " << xoptg[1] << " /J/  " << endl;
  }

  itersom = 0;
  costsaveming = finit;
  irestart2 = 1;

  for (irestart = 1; irestart <= nbrestart; ++irestart) {
    xsave = vinit;
    irestart2 *= 2;
    rho00 = rho000 / irestart2;

    double iter2 = 1;

    for (iter1 = 1; iter1 <= nbext1; ++iter1) {
      costsavemin = finit;
      v = xsave;
      iter2 *= 2;
      rho0 = rho00 / iter2;
      double iterbvp2 = 1;

      for (iterbvp = 1; iterbvp <= nbbvp; ++iterbvp) {
        iterbvp2 *= 2;
        ++itersom;
        x1 = v;
        rho = rho0 / iterbvp2;
        if (debug > 4) {
          cout << "MM " << irestart << " " << iter1 << " " << iterbvp << " " << rho
               << " ------------------------------ \n";
        }

        gradopt(x1, fpx, temp, rho, f, gnorm, fpx0, hgc);

        if (costsaveming < epsij) {
          break;
        }

        if (iterbvp >= 2) {
          tir(v, fpx);
        } else {
          rand(v);
        }

        v0 = v;
      }

      if (debug) {
        cout.precision(15);
        cout << " F = " << costsavemin << " FM = " << costsaveming << endl;
      }

      if (costsaveming < epsij) {
        goto L9101;
      }

      costsave = f;

      if (iter1 >= 2) {
        tir(v, fpx);
      } else {
        rand(v);
      }

      xsave0 = xsave;
    }
  }

L9101:

  result(xoptg, vinit);
  cout << "-------------------------------------------\n";
  cout.precision(15);
  cout << " FM = " << costsaveming << " nb eval J : " << nbeval << " nbevalp : " << nbevalp << endl;
  if (ncstr > 0) {
    cout << "-------------------------------------------\n";
    cout << "F seul = " << fseulopt << endl;
    cout << "-------------------------------------------\n";
    cout << " CSTR = ";

    for (int i = 0; i < ncstr; ++i) {
      cout << cstropt[i] << " ";
    }

    cout << endl;
    cout << "-------------------------------------------\n";
    if (ndim < 20) {
      cout << " x = ";

      for (int i = 0; i < ndim; ++i) {
        cout << xoptg[i] << " ";
      }
    }

    cout << "-------------------------------------------\n";
  }

  xx = xoptg;
  return fseulopt;
}

/*!
 * \brief tir
 * \param v Vect &
 * \param fpx Vect &
 */
void BijanMO::tir(Vect &v, Vect &fpx) {
  for (int i = 0; i < ndim; ++i) {
    double vi = v[i], x0 = xmin[i], x1 = xmax[i], fpxi = -fpx[i];
    fpxi = min(fpxi, (x1 - vi) * 0.95);
    fpxi = max(fpxi, (x0 - vi) * 0.95);
    vi = max(min(vi + fpxi, x1), x0);
    v[i] = vi;
    fpx[i] = fpxi;
  }
}

/*!
 * \brief rand
 * \param v Vect &
 */
void BijanMO::rand(Vect &v) {
  if (diagrand) {
    double xrdran = xrandme(nbeval + nbevalp);

    for (int ii = 0; ii < ndim; ++ii) {
      v(ii) = xmin(ii) + xrdran * (xmax(ii) - xmin(ii));
      v(ii) = max(min(v(ii), xmax(ii)), xmin(ii));
    }
  } else {
    for (int ii = 0; ii < ndim; ++ii) {
      double xrdran = xrandme(nbeval + nbevalp);
      v(ii) = xmin(ii) + xrdran * (xmax(ii) - xmin(ii));
      v(ii) = max(min(v(ii), xmax(ii)), xmin(ii));
    }
  }
}

/*!
 * \brief gradopt
 * \param x1 Vect &
 * \param fpx Vect &
 * \param temp Vect &
 * \param rho double
 * \param f double
 * \param gnorm double
 * \param fpx0 Vect &
 * \param hgc Vect &
 */
int BijanMO::gradopt(Vect &x1, Vect &fpx, Vect &temp, double &rho, double &f, double &gnorm,
                     Vect &fpx0, Vect &hgc) {
  /* Local variables */
  integer ii;
  integer igc, igr;
  double xmod, gamgc;
  double xmodd, xnorm, gnorm0 = 0.;

  /* Function Body */
  f = func(x1);

  igc = typealgo;    // 1 =>   CG  , other : descent
  hgc = 0.;

  for (igr = 1; igr <= nbgrad; ++igr) {
    xnorm = 0.;
    fpx0 = fpx;
    xnorm = fpx.norm( );
    nbeval = -nbeval;
    funcp(x1, fpx, f);
    nbeval = -nbeval;

    gamgc = 0.;

    if (igc == 1 && igr >= 2 && xnorm > 1e-10) {
      for (ii = 0; ii < ndim; ++ii) {
        gamgc += (fpx[ii] - fpx0[ii]) * fpx[ii] / xnorm;
      }
    }

    for (ii = 0; ii < ndim; ++ii) {
      hgc(ii) = fpx(ii) + gamgc * hgc(ii);
    }

    if (debug > 5) {
      cout << "\t\t\t" << rho << " " << hgc(0) << " " << hgc(1) << "\n";
    }

    f = ropt_dicho(x1, temp, rho, hgc, f);

    xmod = 0.;
    xmodd = 0.;

    for (ii = 0; ii < ndim; ++ii) {
      double x0 = x1(ii);
      x1(ii) = x1(ii) - rho * hgc(ii);
      x1(ii) = min(x1(ii), xmax(ii));
      x1(ii) = max(x1(ii), xmin(ii));
      xmod = xmod + abs(x1(ii) - x0);
      xmodd = xmodd + abs(x1(ii));
    }

    f = func(x1);

    gnorm = fpx.l2( );
    if (igr == 1) {
      gnorm0 = gnorm;
    }

    if (gnorm0 < 1e-6) {
      return 0;
    }

    gnorm /= gnorm0;
    if (histpath) {
      ofstream fhist(histpath->c_str( ), ios::app);
      fhist.precision(16);
      fhist << f << " " << gnorm * gnorm0 << " ";
      int n1 = min(ndim, 10);

      for (int i = 0; i < n1; ++i) {
        fhist << x1[i] << " ";
      }
    }

    if (histcpath) {
      ofstream fhist(histcpath->c_str( ), ios::app);
      fhist.precision(16);
      fhist << fseul << endl;

      for (int i = 0; i < ncstr; ++i) {
        fhist << cstr[i] << (i % 4 ? '\n' : '\t') << '\t';
      }
    }

    if (debug > 2) {
      cout << "\t\t\t " << f << " " << gnorm * gnorm0 << " " << x1[0] << " " << x1[1] << " /J/ "
           << endl;
    }

    if (f < costsaveming) {
      costsaveming = f;
      gnormsave = gnorm;
      cstropt = cstr;
      fseulopt = fseul;
      xoptg = x1;
    }

    if (f < costsavemin) {
      costsavemin = f;
      xopt1 = x1;
    }

    if (f < epsij) {
      break;
    }

    if (gnorm < 1e-6 || gnorm * gnorm0 < 1e-6) {
      break;
    }
  }

  if (debug > 3) {    //
    cout << "\t\t\t opt: rho = " << rho << " F = " << f << endl;
  }

  return 0;
}

/*!
 * \brief ropt_dicho
 * \param x Vect
 * \param temp Vect
 * \param ro double &
 * \param rho double
 * \param g Vect
 * \param ccout double
 * \return double
 */
double BijanMO::ropt_dicho(Vect x, Vect temp, double &ro, Vect g, double ccout) {
  integer j, l;
  double s, fm, sd, sn, pr;
  static double fmin[3] = {0, 0, 0};
  integer numi;
  double romin[3];
  integer numimax;

  /* Function Body */
  numi = 0;
  numimax = 5;
L240:
  romin[0] = ro * .5;
  romin[1] = ro;
  romin[2] = ro * 2.;
  l = 0;
L300:
  fmin[l] = fun(x, temp, g, romin[l]);
  ++l;
  ++numi;
  if (l == 1 & fmin[0] > ccout) {
    ro *= .5;
    /* ******  test d'arret */
    if (abs(ro) < 1e-5 || numi > numimax) {
      goto L500;
    }

    goto L240;
  }

L360:
  if (l < 2) {
    goto L300;
  }

  if (fmin[0] < fmin[1]) {
    goto L380;
  }

L370:
  if (l < 3) {
    goto L300;
  }

  if (fmin[1] <= fmin[2]) {
    goto L450;
  } else {
    goto L420;
  }

L380:
  l = 3;
  romin[2] = romin[1];
  fmin[2] = fmin[1];
  romin[1] = romin[0];
  fmin[1] = fmin[0];
  romin[0] *= .5;
  ++numi;
  fmin[0] = fun(x, temp, g, romin[0]);
  goto L360;
L420:
  romin[0] = romin[1];
  fmin[0] = fmin[1];
  romin[1] = romin[2];
  fmin[1] = fmin[2];
  romin[2] *= 2.;
  ++numi;
  fmin[2] = fun(x, temp, g, romin[2]);
  goto L370;
L450:
  ro = romin[1];
  if (abs(fmin[1] - fmin[2]) * 2 / (fmin[1] + fmin[2]) < 1e-4 || numi > numimax) {
    goto L500;
  }

  /* ****** calcul de ro interpole */
  sn = 0.;
  sd = 0.;

  for (int i = 0; i < 3; ++i) {
    s = 0.;
    pr = 1.;

    for (j = 0; j < 3; ++j) {
      if (i != j) {
        s += romin[j];
        pr *= romin[i] - romin[j];
      }
    }

    sn += fmin[i] * s / pr;
    sd += fmin[i] / pr;
  }

  ro = sn / sd / 2.;
  if (debug > 5) {
    cout << "\t\t\t\tro int  = " << ro << " " << l << endl;
  }

L500:
  fm = fun(x, temp, g, ro);
  ccout = fm;
  if (fm > fmin[1]) {
    ro = romin[1];
    ccout = fmin[1];
  }

  if (debug > 4) {
    cout << "\t\t\t\tdicho : " << ro << " " << ccout << " " << l << endl;
  }

  return ccout;
} /* ropt_dicho__ */

/*!
 * \brief fun
 * \param x Vect &
 * \param temp Vect &
 * \param ro double
 * \return double
 */
double BijanMO::fun(Vect &x, Vect &temp, Vect &g, double ro) {
  for (int ii = 0; ii < ndim; ++ii) {
    temp[ii] = x[ii] - ro * g[ii];
    temp[ii] = max(min(temp[ii], xmax[ii]), xmin[ii]);
  }

  if (debug > 5) {
    cout << "                ro = " << ro << endl;
  }

  return func(temp);
}

/*!
 * \brief funcp
 * \param x Vect &
 * \param fpx Vect &
 * \param double f
 */
void BijanMO::funcp(Vect &x, Vect &fpx, double f) {
  /* Local variables */
  double fp, x00;

  /* Function Body */
  nbevalp = nbevalp + 1;
  double *ok = DJ(x, fpx);
  if (!ok) {
    for (int ii = 0; ii < ndim; ++ii) {
      x00 = x[ii];
      double epsifd = max(min(abs(x00) * epsfd, 100 * epsfd), epsfd / 100.);
      if (x00 + epsifd <= xmax[ii]) {
        x[ii] = x00 + epsifd;
        fp = func(x);
      } else {
        x[ii] = x00 - epsifd;
        fp = func(x);
        epsifd = -epsifd;
      }

      fpx[ii] = (fp - f) / epsifd;
      x[ii] = x00;
    }
  }
}

/*!
 * \brief funcapp
 * \param x Vect &
 * \param fpx Vect &
 * \return double
 */
double BijanMO::funcapp(Vect &x, Vect &fpx) {
  /* Local variables */
  integer kk;
  double diffucarte;
  integer nbevalsave;
  double diffucarte0;
  integer kkk;
  double xcoef, fapp = 0.;

  /* Function Body */
  nbevalsave = min(nbeval, nbsol);

  diffucarte0 = 100.;
  diffucarte = diffucarte0;
  double itest2 = 1.;

  for (int itest = 0; itest <= 5; ++itest) {
    itest2 *= 2;
    fapp = 0.;
    fpx = 0.;
    xcoef = 0.;

    for (kk = 0; kk < nbevalsave; ++kk) {
      double d = 0.;

      for (kkk = 0; kkk < ndim; ++kkk) {
        double dd = (x[kkk] - xfeval(kk, kkk)) / (xmax[kkk] - xmin[kkk]);
        d += dd * dd;
      }

      double vloc = exp(-d * diffucarte);
      fapp += feval[kk] * vloc;

      for (kkk = 0; kkk < ndim; ++kkk) {
        double xcc = (x[kkk] - xfeval(kk, kkk)) / (xmax[kkk] - xmin[kkk]);
        fpx[kkk] -= diffucarte * 2 * xcc * vloc;
      }

      xcoef += vloc;
    }

    if (xcoef > 1e-6) {
      fapp /= xcoef;
      fpx /= xcoef;
      break;
    } else {
      diffucarte = diffucarte0 / itest2;
    }
  }

  if (debug > 3) {
    cout << "                fapp = " << fapp << " " << nbeval << x[0] << " " << x[1] << endl;
  }

  return fapp;
}
