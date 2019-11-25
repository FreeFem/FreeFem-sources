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

// compile and link with ff-c++ metric_Pk.cpp

// *INDENT-OFF* //
// ff-c++-LIBRARY-dep:
// ff-c++-cpp-dep:
// *INDENT-ON* //

#include <iostream>
#include <cfloat>
#include <cmath>
using namespace std;
#include "ff++.hpp"
using namespace Fem2D;

#include <vector>

#include "TensorK.hpp"
// the main class

class MetricPk : public E_F0mps {
 public:
  static basicAC_F0::name_and_type name_param[];
  static const int n_name_param = 10;
  Expression nargs[n_name_param];    // stocker les argunments nommes

  typedef KN_< double > Result;
  Expression expTh;
  Expression expu;

  MetricPk(const basicAC_F0 &args) {
    args.SetNameParam(n_name_param, name_param, nargs);    // les arguments nommes
    expTh = to< pmesh >(args[0]);                          // a the expression to get the mesh
    expu = to< double >(args[1]);                          // a the expression to get the mesh
  }

  double arg(int i, Stack stack, double a) const {
    return nargs[i] ? GetAny< double >((*nargs[i])(stack)) : a;
  }

  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

  bool arg(int i, Stack stack, bool a) const {
    return nargs[i] ? GetAny< bool >((*nargs[i])(stack)) : a;
  }

  KN< double > *arg(int i, Stack stack, KN< double > *a) const {
    return nargs[i] ? GetAny< KN< double > * >((*nargs[i])(stack)) : a;
  }

  ~MetricPk( ) {}

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< pmesh >( ), atype< double >( ));
    ;
  }

  static E_F0 *f(const basicAC_F0 &args) { return new MetricPk(args); }

  AnyType operator( )(Stack s) const;    // la vraie fonction qui fait faire le boulot
};

basicAC_F0::name_and_type MetricPk::name_param[MetricPk::n_name_param] = {
  {"kDeg", &typeid(long)},
  {"rDeg", &typeid(long)},
  {"iterJacobiDeriv", &typeid(long)},
  {"iterJacobiMetric", &typeid(long)},
  {"Derivatives", &typeid(KN< double > *)},
  {"rmax", &typeid(double)},
  {"mass", &typeid(double)},
  {"TriangulationType", &typeid(long)},
  {"MetricType", &typeid(long)},
  {"pExp", &typeid(double)}};
AnyType MetricPk::operator( )(Stack stack) const {
  /************* r�cup�ration des arguments ****************/

  const long k_deg =
    arg(0, stack, 2L);             // Finite element of degree k_deg will be used for approximation.
  const long m_deg = k_deg + 1;    // Derivatives of degree m_deg = k_deg+1 will be estimated.
  const long m_dim =
    m_deg + 1;    // The description of these derivatives requires m_dim = m_deg+1 coefficients.
  const long nDOFt =
    ((k_deg + 1) * (k_deg + 2)) / 2;    // Number of Lagrange points on each triangle.
  const long r_deg =
    arg(1, stack, 1L);    // The function is approximated in the W^r,p Sobolev semi-norm.
  const double p_exp = arg(9, stack, 2.);
  const long iterJacobiDeriv =
    arg(2, stack, 3L);    // The derivatives are slightly smoothed before use.
  const long iterJacobiMetric =
    arg(3, stack, 3L);    // The riemannian metric is slightly smoothed before being returned.
  const double rmax = arg(5, stack, 1.);    // Not used yet. (Lower bound for the metric)
  const double mass =
    arg(6, stack, 1000.);    // Mass of the metric returned, i.e. mass = int sqrt(det M).
  // In practice, bamg produces a mesh with nt=2*mass elements.
  const TensorK::triangulation_type ttype = TensorK::triangulation_type(
    arg(7, stack, long(0)));    // Type of triangulation on which approx will be done.
  const TensorK::which_matrix wmat =
    TensorK::which_matrix(arg(8, stack, long(1)));    // Type of metric. Do not change.
  TensorK tk(m_deg, r_deg, ttype, wmat, p_exp);

  cout << "Approximation of " << r_deg << "th derivatives using finite elements of degree " << k_deg
       << ", in the L^" << p_exp << " norm.\n";
  cout
    << "Triangulation type : " << ttype
    << "; Graded=0, Quasi_Acute(refined)=1, Quasi_Acute_Unrefined=2, Quasi_Acute_Proved(refined)=3"
    << endl;

  const Mesh *pTh = GetAny< pmesh >((*expTh)(stack));
  ffassert(pTh);
  const Mesh &Th = *pTh;

  /************ initialisations ***************/

  // the metric
  const int nv = Th.nv;
  KN< double > *pMetric = new KN< double >(nv * 3);
  KN< double > &metric = *pMetric;
  metric = 0.;

  if (!tk.is_valid) {
    cout << "Error : Unsupported parameters for MetricPk!\n";
    Add2StackOfPtr2Free(stack, pMetric);
    return SetAny< KN< double > >(
      metric);    // identically zero metric is returned in case of error
  }

  // Lagrange points.
  const R2 QLagrange[5][15] = {
    {R2(1. / 3., 1. / 3.)},    // k_deg=0, barycenter.
    {R2(0, 0), R2(1, 0), R2(0, 1)},
    {R2(0, 0), R2(1, 0), R2(0, 1), R2(0.5, 0.5), R2(0, 0.5), R2(0.5, 0)},
    {R2(0, 0), R2(1, 0), R2(0, 1), R2(2. / 3., 1. / 3.), R2(1. / 3., 2. / 3.), R2(0, 2. / 3.),
     R2(0, 1. / 3.), R2(1. / 3., 0), R2(2. / 3., 0), R2(1 / 3., 1 / 3.)},
    {R2(0, 0), R2(1, 0), R2(0, 1), R2(3. / 4., 1. / 4.), R2(1. / 2., 1. / 2.), R2(1. / 4., 3. / 4.),
     R2(0., 3. / 4.), R2(0., 1. / 2.), R2(0., 1. / 4.), R2(1. / 4., 0), R2(1. / 2., 0),
     R2(3. / 4., 0), R2(1. / 4., 1. / 2.), R2(1. / 2., 1. / 4.), R2(1. / 4., 1. / 4.)}};
  std::vector< R > aires;
  aires.resize(nv);    // area of the "cell" surrounding a point
  std::vector< R > Deriv;
  Deriv.resize(m_dim * nv);    // estimate of the derivatives at a point
  std::vector< R > DOFt;
  DOFt.resize(nDOFt);    // degrees of freedom on a triangle

  // le bord
  // Th.ElementAdj(k,ie); renvoie k', �crase ie par ie';
  // is frontiere k' <0 ou k'==k
  // Seg
  // BoundaryElement (+r�cent)
  // Th.nbe ou Th.neb
  // Th.be(i)  Th.be(i)[0] ou Th.be(i)[1]
  std::vector< int > nextv;
  nextv.resize(nv);
  fill(nextv.begin( ), nextv.end( ), -1);

  for (int i = 0; i < Th.nt; ++i) {
    for (int ie = 0; ie < 3; ++ie) {
      int iee = ie;
      const int j = Th.ElementAdj(i, iee);
      if (j == i || j < 0) {
        nextv[Th(i, (ie + 1) % 3)] =
          Th(i, (ie + 2) % 3);    // nextv[Th(i,(ie+1)%3)]=Th(i,(ie+2)%3);
      }
    }
  }

  /********* estimation des d�riv�es ***********/
  for (int i = 0; i < Th.nt; ++i) {
    const Triangle &K = Th[i];
    const R2 sommets[3] = {K(R2(0, 0)), K(R2(1, 0)), K(R2(0, 1))};
    const R aire = det(sommets[0], sommets[1], sommets[2]) / 2;
    const R2 invHauteur[3] = {perp(sommets[1] - sommets[2]) / (2 * aire),
                              perp(sommets[2] - sommets[0]) / (2 * aire),
                              perp(sommets[0] - sommets[1]) / (2 * aire)};

    for (int dof = 0; dof < nDOFt; ++dof) {    // Degrees of freedom on the triangle
      const R2 Point = K(QLagrange[k_deg][dof]);
      MeshPointStack(stack)->set(Point.x, Point.y);
      DOFt[dof] = GetAny< double >((*expu)(stack));
    }

    double f[m_deg];    // contains the derivatives of order

    switch (m_deg) {
      case 2: {    // accolades n�cessaires pour d�clarer des variables dans un case
        double f[2];
        tk.getDerivatives(DOFt, invHauteur, f);    // f={fx,fy}

        for (int j = 0; j < 3; ++j) {
          const int s = Th(i, j);    // le sommet j du triangle i
          Deriv[m_dim * s + 0] +=
            aire * f[0] * invHauteur[j].x;    // contribution � l'estimation des d�riv�es secondes.
          Deriv[m_dim * s + 1] +=
            aire * (f[0] * invHauteur[j].y / 2. + f[1] * invHauteur[j].x / 2.);
          Deriv[m_dim * s + 2] += aire * f[1] * invHauteur[j].y;
          aires[s] += aire;
        }

        break;
      }
      case 3: {
        double f[3];
        tk.getDerivatives(DOFt, invHauteur, f);    // f={fxx,fxy,fyy}

        for (int j = 0; j < 3; ++j) {
          const int s = Th(i, j);
          Deriv[m_dim * s + 0] +=
            aire * f[0] *
            invHauteur[j].x;    // contribution � l'estimation des d�riv�es troisi�mes.
          Deriv[m_dim * s + 1] +=
            aire * (f[0] * invHauteur[j].y / 3. + f[1] * invHauteur[j].x * 2. / 3.);
          Deriv[m_dim * s + 2] +=
            aire * (f[2] * invHauteur[j].x / 3. + f[1] * invHauteur[j].y * 2. / 3.);
          Deriv[m_dim * s + 3] += aire * f[2] * invHauteur[j].y;
          aires[s] += aire;
        }

        break;
      }
      case 4: {
        double f[4];
        tk.getDerivatives(DOFt, invHauteur, f);    // f={fxxx,fxxy,fxyy,fyyy}

        for (int j = 0; j < 3; ++j) {
          const int s = Th(i, j);
          Deriv[m_dim * s + 0] += f[0] * invHauteur[j].x *
                                  aire;    // contribution � l'estimation des d�riv�es quatri�mes.
          Deriv[m_dim * s + 1] +=
            (f[0] * invHauteur[j].y / 4. + f[1] * invHauteur[j].x * 3. / 4.) * aire;
          Deriv[m_dim * s + 2] +=
            (f[1] * invHauteur[j].y / 2. + f[2] * invHauteur[j].x / 2.) * aire;
          Deriv[m_dim * s + 3] +=
            (f[2] * invHauteur[j].y * 3. / 4. + f[3] * invHauteur[j].x / 4.) * aire;
          Deriv[m_dim * s + 4] += f[3] * invHauteur[j].y * aire;
          aires[s] += aire;
        }

        break;
      }
      case 5: {
        double f[5];
        tk.getDerivatives(DOFt, invHauteur, f);    // f={fxxxx,fxxxy,fxxyy,fxyyy,fyyyy}

        for (int j = 0; j < 3; ++j) {
          const int s = Th(i, j);
          Deriv[m_dim * s + 0] += f[0] * invHauteur[j].x *
                                  aire;    // contribution � l'estimation des d�riv�es quatri�mes.
          Deriv[m_dim * s + 1] +=
            (f[0] * invHauteur[j].y / 5. + f[1] * invHauteur[j].x * 4. / 5.) * aire;
          Deriv[m_dim * s + 2] +=
            (f[1] * invHauteur[j].y * 2. / 5. + f[2] * invHauteur[j].x * 3. / 5.) * aire;
          Deriv[m_dim * s + 3] +=
            (f[2] * invHauteur[j].y * 3. / 5. + f[3] * invHauteur[j].x * 2. / 5.) * aire;
          Deriv[m_dim * s + 4] +=
            (f[3] * invHauteur[j].y * 4. / 5. + f[4] * invHauteur[j].x / 5.) * aire;
          Deriv[m_dim * s + 5] += f[4] * invHauteur[j].y * aire;
          aires[s] += aire;
        }
      }    // case m_deg==5
    }      // switch m_deg
  }        // for i triangle

  for (int i = 0; i < nv; ++i) {
    for (int j = 0; j < m_dim; ++j) {
      Deriv[m_dim * i + j] *= 3 / aires[i];
    }
  }

  // Estimating derivatives on the boundary by averaging neighboring estimates in the interior.
  // First the graph distance to the interior is estimated, using Dijkstra's algorithm.
  {
    multimap< int, int > connectivity;    // this is probably already computed by FreeFem. Where ?
    vector< pair< int, int > > dist;      // (distance,number)
    set< int > computed;

    for (int i = 0; i < Th.nt; ++i) {    // Obtaining the connectivity
      const int u = Th(i, 0), v = Th(i, 1), w = Th(i, 2);
      if (nextv[u] == -1 && nextv[v] == -1 && nextv[w] == -1) {
        continue;    // attention is restricted to the boundary.
      }

      connectivity.insert(pair< int, int >(u, v));
      connectivity.insert(pair< int, int >(v, u));
      connectivity.insert(pair< int, int >(v, w));
      connectivity.insert(pair< int, int >(w, v));
      connectivity.insert(pair< int, int >(w, u));
      connectivity.insert(pair< int, int >(u, w));
      if (nextv[u] == -1 && computed.insert(u).second == true) {
        dist.push_back(pair< int, int >(u, 0));
      }

      if (nextv[v] == -1 && computed.insert(v).second == true) {
        dist.push_back(pair< int, int >(v, 0));
      }

      if (nextv[w] == -1 && computed.insert(w).second == true) {
        dist.push_back(pair< int, int >(w, 0));
      }
    }

    for (int i = 0; i < dist.size( ); ++i) {    // Dijkstra's algorithm.
      const int u = dist[i].first;
      const int du = dist[i].second;
      computed.insert(u);
      const pair< multimap< int, int >::iterator, multimap< int, int >::iterator > ret =
        connectivity.equal_range(u);

      for (multimap< int, int >::iterator it = ret.first; it != ret.second; ++it) {
        const int v = it->second;
        if (computed.insert(v).second) {
          dist.push_back(pair< int, int >(v, du + 1));
        }
      }
    }

    map< int, int > dist_sorted;
    dist_sorted.insert(dist.begin( ), dist.end( ));

    for (int i = 0; i < dist.size( ); ++i) {
      const int u = dist[i].first;
      const int du = dist[i].second;
      if (du == 0) {
        continue;
      }

      const pair< multimap< int, int >::iterator, multimap< int, int >::iterator > neighbors =
        connectivity.equal_range(u);
      int closer_neighbors = 0;

      for (multimap< int, int >::iterator it = neighbors.first; it != neighbors.second; ++it) {
        const int v = it->second;
        const int dv = dist_sorted.find(v)->second;
        if (du != dv + 1) {
          continue;
        }

        closer_neighbors++;

        for (int k = 0; k < m_dim; ++k) {
          Deriv[m_dim * u + k] += Deriv[m_dim * v + k];
        }
      }

      for (int k = 0; k < m_dim; ++k) {
        Deriv[m_dim * u + k] /= closer_neighbors;
      }
    }
  }

  // Averaging derivatives
  std::vector< int > cardNeighbors;
  cardNeighbors.resize(nv);

  for (int i = 0; i < Th.nt; ++i) {
    for (int j = 0; j < 3; ++j) {
      cardNeighbors[Th(i, j)] += 3;
    }
  }

  {
    std::vector< R > DerivNew;
    DerivNew.resize(m_dim * nv);

    for (int r = 0; r < iterJacobiDeriv; ++r) {
      fill(DerivNew.begin( ), DerivNew.end( ), 0);

      for (int i = 0; i < Th.nt; ++i) {
        const int u = Th(i, 0);
        const int v = Th(i, 1);
        const int w = Th(i, 2);

        for (int k = 0; k < m_dim; ++k) {
          const R sum = Deriv[m_dim * u + k] + Deriv[m_dim * v + k] + Deriv[m_dim * w + k];
          DerivNew[m_dim * u + k] += sum;
          DerivNew[m_dim * v + k] += sum;
          DerivNew[m_dim * w + k] += sum;
        }    // for k
      }      // for i triangles

      for (int i = 0; i < nv; ++i) {
        for (int k = 0; k < m_dim; ++k) {
          Deriv[m_dim * i + k] = DerivNew[m_dim * i + k] / cardNeighbors[i];
        }
      }
    }
  }

  // Exporting the derivatives, if required
  KN< double > *pDerivRes = arg(4, stack, (KN< double > *)0);
  if (pDerivRes) {
    ffassert(pDerivRes->N( ) == m_dim * nv);
    KN< double > &DerivRes = *pDerivRes;

    for (int i = 0; i < m_dim * nv; ++i) {
      DerivRes[i] = Deriv[i];
    }
  }

  /**************** Computing the "L^infinity homogeneous" metric ************/

  std::vector< R > ih_metric;
  ih_metric.resize(3 * nv);

  for (int i = 0; i < nv; ++i) {
    tk.getM(&Deriv[m_dim * i], &ih_metric[3 * i]);
  }

  // Jacobi iterations
  // the power -1/2 of the metrics, which is homogenous to a distance, is averaged.
  // for(int i=0; i<nv; ++i) if(ih_metric[3*i]==0 && ih_metric[3*i+1]==0 && ih_metric[3*i+2]==0) {
  // cout << Th(i).x << "," << Th(i).y << endl;
  // for(int j=0; j<m_dim; ++j) cout << Deriv[m_dim*i+j] << ",";
  // cout << "jacobi deriv : " << iterJacobiDeriv << endl;
  // }

  // Averaging the metric
  {
    for (int i = 0; i < nv; ++i) {
      TensorK::PowSym(&ih_metric[3 * i], -0.5);
    }

    std::vector< double > ih_metricNew;
    ih_metricNew.resize(3 * nv);

    for (int r = 0; r < iterJacobiMetric; ++r) {
      fill(ih_metricNew.begin( ), ih_metricNew.end( ), 0);

      for (int i = 0; i < Th.nt; ++i) {
        const int u = Th(i, 0);
        const int v = Th(i, 1);
        const int w = Th(i, 2);

        for (int k = 0; k < 3; ++k) {
          const R sum = ih_metric[3 * u + k] + ih_metric[3 * v + k] + ih_metric[3 * w + k];
          ih_metricNew[3 * u + k] += sum;
          ih_metricNew[3 * v + k] += sum;
          ih_metricNew[3 * w + k] += sum;
        }    // for k
      }      // for i triangles

      for (int i = 0; i < nv; ++i) {
        for (int k = 0; k < 3; ++k) {
          ih_metric[3 * i + k] = ih_metricNew[3 * i + k] / cardNeighbors[i];
        }
      }
    }

    for (int i = 0; i < nv; ++i) {
      TensorK::PowSym(&ih_metric[3 * i], -2);
    }
  }

  /**************** Multiplicator to obtain prescribed mass and to balance errors ****************/
  // Note : future versions could include hmin and hmax parameters.
  {
    const long d_dim = 2;    // Space Dimension : 2

    // M : Metric for the L^p norm, M_0 metric for the L^infty norm
    const R alpha = -1. / ((m_deg - r_deg) * p_exp + d_dim);    // M = (det M_0)^alpha M_0
    const R beta =
      (m_deg - r_deg) * p_exp /
      (((m_deg - r_deg) * p_exp + d_dim) *
       2.);    // integrate sqrt(det M) = integrate (det M_0)^beta; beta = (alpha*d_dim+1)/2.

    for (int i = 0; i < nv; ++i) {
      aires[i] /= 3;
    }

    R totalArea = 0;

    for (int i = 0; i < nv; ++i) {
      totalArea += aires[i];
    }

    R totalMass = 0;

    for (int i = 0; i < nv; ++i) {
      totalMass +=
        aires[i] * pow(TensorK::det(&ih_metric[3 * i]), beta);    // Eigen[i][0]*Eigen[i][1]
    }

    const R gamma =
      pow(mass / totalMass,
          1. / (beta * d_dim));    // gamma satisfies integrate (det (gamma M_0) )^beta = mass

    for (int i = 0; i < nv; ++i) {
      const R lambda = gamma * pow(square(gamma) * TensorK::det(&ih_metric[3 * i]),
                                   alpha);    // Eigen[i][0]*Eigen[i][1]

      for (int j = 0; j < 3; ++j) {
        metric[3 * i + j] = lambda * ih_metric[3 * i + j];
      }
    }

    R obtainedMass = 0;

    for (int i = 0; i < nv; ++i) {
      obtainedMass += aires[i] * sqrt(TensorK::det(metric + 3 * i));
    }

    cout << "Desired Mass : " << mass << "; obtained mass " << obtainedMass << endl;
  }

  Add2StackOfPtr2Free(stack, pMetric);
  return SetAny< KN< double > >(metric);

  /* //example by F.Hecht
   * for( int k=0;k<Th.nt; ++k)
   * {
   *  const Triangle &K=Th[k];
   *  for (int i=0;i<6;++i)
   *  {
   *      R2 Pi=K(Q2[i]);
   *      MeshPointStack(stack)->set(Pi.x,Pi.y);
   *      R uP= GetAny<double>((*expu)(stack));
   *      cout << Pi << " " << uP << endl;
   *      if(i<3)
   *      { int s= Th(k,i); // le sommet i du triangle k.
   *          metric[3*s+0] =uP;
   *          metric[3*s+1] =uP+1;
   *          metric[3*s+2] =uP+3;
   *      }
   *
   *  }
   * }
   */
}

static void Load_Init( ) {
  cout << "\n  -- lood: init MetricPk\n";
  Global.Add("MetricPk", "(", new OneOperatorCode< MetricPk >( ));
}

LOADFUNC(Load_Init)
