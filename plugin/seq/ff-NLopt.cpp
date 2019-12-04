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
// AUTHOR  : Sylvain Auliac
// E-MAIL  : auliac@ann.jussieu.fr

/* clang-format off */
//ff-c++-LIBRARY-dep: nlopt
//ff-c++-cpp-dep:
/* clang-format on */

#include <iostream>
#include <stack>
#include <vector>
using namespace std;
#include "ff++.hpp"

#include <nlopt.hpp>

extern Block *currentblock;

typedef double R;
typedef double (*NLoptFuncType)(unsigned, const double *, double *, void *);

template< nlopt::algorithm ALGO >
struct Info {
  static const bool DF = true, SA = false;
  static const char *name;
};    // DF=Derivative Free , SA=need a Sub Algorithm
template<>
struct Info< nlopt::LD_LBFGS_NOCEDAL > {
  static const bool DF = false, SA = false;
  static const char *name;
};
template<>
struct Info< nlopt::LD_LBFGS > {
  static const bool DF = false, SA = false;
  static const char *name;
};
template<>
struct Info< nlopt::LD_VAR1 > {
  static const bool DF = false, SA = false;
  static const char *name;
};
template<>
struct Info< nlopt::LD_VAR2 > {
  static const bool DF = false, SA = false;
  static const char *name;
};
template<>
struct Info< nlopt::LD_TNEWTON > {
  static const bool DF = false, SA = false;
  static const char *name;
};
template<>
struct Info< nlopt::LD_TNEWTON_RESTART > {
  static const bool DF = false, SA = false;
  static const char *name;
};
template<>
struct Info< nlopt::LD_TNEWTON_PRECOND > {
  static const bool DF = false, SA = false;
  static const char *name;
};
template<>
struct Info< nlopt::LD_TNEWTON_PRECOND_RESTART > {
  static const bool DF = false, SA = false;
  static const char *name;
};
template<>
struct Info< nlopt::LD_MMA > {
  static const bool DF = false, SA = false;
  static const char *name;
};
template<>
struct Info< nlopt::LD_AUGLAG > {
  static const bool DF = false, SA = true;
  static const char *name;
};
template<>
struct Info< nlopt::LD_AUGLAG_EQ > {
  static const bool DF = false, SA = true;
  static const char *name;
};
template<>
struct Info< nlopt::LD_SLSQP > {
  static const bool DF = false, SA = false;
  static const char *name;
};
template<>
struct Info< nlopt::GD_STOGO > {
  static const bool DF = false, SA = false;
  static const char *name;
};
template<>
struct Info< nlopt::GD_STOGO_RAND > {
  static const bool DF = false, SA = false;
  static const char *name;
};
template<>
struct Info< nlopt::GD_MLSL > {
  static const bool DF = false, SA = true;
  static const char *name;
};
template<>
struct Info< nlopt::GD_MLSL_LDS > {
  static const bool DF = false, SA = true;
  static const char *name;
};
template<>
struct Info< nlopt::GN_MLSL > {
  static const bool DF = true, SA = true;
  static const char *name;
};
template<>
struct Info< nlopt::GN_MLSL_LDS > {
  static const bool DF = true, SA = true;
  static const char *name;
};
template<>
struct Info< nlopt::LN_AUGLAG > {
  static const bool DF = true, SA = true;
  static const char *name;
};
template<>
struct Info< nlopt::LN_AUGLAG_EQ > {
  static const bool DF = true, SA = true;
  static const char *name;
};
template<>
struct Info< nlopt::G_MLSL > {
  static const bool DF = true, SA = true;
  static const char *name;
};
template<>
struct Info< nlopt::G_MLSL_LDS > {
  static const bool DF = true, SA = true;
  static const char *name;
};
template<>
struct Info< nlopt::AUGLAG > {
  static const bool DF = true, SA = true;
  static const char *name;
};
template<>
struct Info< nlopt::AUGLAG_EQ > {
  static const bool DF = true, SA = true;
  static const char *name;
};

template< nlopt::algorithm ALGO >
const char *Info< ALGO >::name = "ALGORITHM";
template<>
const char *Info< nlopt::GN_DIRECT >::name = "Dividing Rectangles";
template<>
const char *Info< nlopt::GN_DIRECT_L >::name = "Locally Biased Dividing Rectangles";
template<>
const char *Info< nlopt::GN_DIRECT_L_RAND >::name = "Randomized Locally Biased Dividing Rectangles";
template<>
const char *Info< nlopt::GN_DIRECT_NOSCAL >::name = "Dividing Rectangles (no scaling)";
template<>
const char *Info< nlopt::GN_DIRECT_L_NOSCAL >::name =
  "Locally Biased Dividing Rectangles (no scaling)";
template<>
const char *Info< nlopt::GN_DIRECT_L_RAND_NOSCAL >::name =
  "Randomized Locally Biased Dividing Rectangles (no scaling)";
template<>
const char *Info< nlopt::GN_ORIG_DIRECT >::name = "Original Glabonsky's Dividing Rectangles";
template<>
const char *Info< nlopt::GN_ORIG_DIRECT_L >::name =
  "Original Glabonsky's Locally Biased Dividing Rectangles";
const char *Info< nlopt::GD_STOGO >::name = "StoGO";
const char *Info< nlopt::GD_STOGO_RAND >::name = "Randomized StoGO";
const char *Info< nlopt::LD_LBFGS_NOCEDAL >::name = "Nocedal's Low-Storage BFGS";
const char *Info< nlopt::LD_LBFGS >::name = "Low-Storage BFGS";
template<>
const char *Info< nlopt::LN_PRAXIS >::name = "Principal Axis";
const char *Info< nlopt::LD_VAR1 >::name = "Rank-1 Shifted Limited Memory Variable Metric";
const char *Info< nlopt::LD_VAR2 >::name = "Rank-2 Shifted Limited Memory Variable Metric";
const char *Info< nlopt::LD_TNEWTON >::name = "Truncated Newton";
const char *Info< nlopt::LD_TNEWTON_RESTART >::name =
  "Steepest Descent Restarting Truncated Newton";
const char *Info< nlopt::LD_TNEWTON_PRECOND >::name = "BFGS Preconditionned Truncated Newton";
const char *Info< nlopt::LD_TNEWTON_PRECOND_RESTART >::name =
  "BFGS Precondionned Truncated Newton with Steepest Descent Resrtarting";
template<>
const char *Info< nlopt::GN_CRS2_LM >::name = "Controlled Random Search with Local Mutation";
const char *Info< nlopt::GN_MLSL >::name = "Multi-Level Single-Linkage (derivative free)";
const char *Info< nlopt::GD_MLSL >::name =
  "Multi-Level Single-Linkage (with gradient-based local search)";
const char *Info< nlopt::GN_MLSL_LDS >::name =
  "Low Discrepancy Sequence Multi-Level Single-Linkage (derivative free)";
const char *Info< nlopt::GD_MLSL_LDS >::name =
  "Low Discrepancy Sequence Multi-Level Single-Linkage (with gradient-based local search)";
const char *Info< nlopt::LD_MMA >::name = "Method of Moving Asymptotes";
template<>
const char *Info< nlopt::LN_COBYLA >::name = "Constrained Optimization by Linear Approximations";
template<>
const char *Info< nlopt::LN_NEWUOA >::name = "NEWUOA";
template<>
const char *Info< nlopt::LN_NEWUOA_BOUND >::name = "NEWUOA for bounded optimization";
template<>
const char *Info< nlopt::LN_NELDERMEAD >::name = "Nelder-Mead Simplex";
template<>
const char *Info< nlopt::LN_SBPLX >::name = "Subplex";
const char *Info< nlopt::LN_AUGLAG >::name =
  "Inequality/Equality Constraints Augmented Lagrangian (derivative free)";
const char *Info< nlopt::LD_AUGLAG >::name =
  "Inequality/Equality Constraints Augmented Lagrangian (with gradient-based subsidiary search)";
const char *Info< nlopt::LN_AUGLAG_EQ >::name =
  "Equality Constraints Augmented Lagrangian (derivative free)";
const char *Info< nlopt::LD_AUGLAG_EQ >::name =
  "Equality Constraints Augmented Lagrangian (with gradient-based subsidiary search)";
template<>
const char *Info< nlopt::LN_BOBYQA >::name = "BOBYQA";
template<>
const char *Info< nlopt::GN_ISRES >::name = "Improved Stochastic Ranking Evolution Strategy";
const char *Info< nlopt::LD_SLSQP >::name = "Sequential Least-Squares Quadratic Programming";
const char *Info< nlopt::G_MLSL >::name = "Multi-Level Single-Linkage";
const char *Info< nlopt::G_MLSL_LDS >::name = "Low Discrepancy Multi-level Single-Linkage";
const char *Info< nlopt::AUGLAG >::name = "Inequality/Equality Constraints Augmented Lagrangian";
const char *Info< nlopt::AUGLAG_EQ >::name = "Equality Constraints Augmented Lagrangian";
inline void Sonde(int i) { cout << "sonde " << i << endl; }

typedef KN_< R > Rn_;
typedef KN< R > Rn;
typedef KNM_< R > Rnm_;
typedef KNM< R > Rnm;

/*template<class T> inline std::vector<T> KnToStdVect(const KN<T> &V)
 * {
 *      std::vector<T> v(V.n);
 *      for(int i=0;i<v.size();++i) v[i] = V[i];
 *      return v;
 * }*/

template< class T >
inline std::vector< T > KnToStdVect(const KN_< T > &V) {
  std::vector< T > v(V.n);

  for (int i = 0; i < v.size( ); ++i) {
    v[i] = V[i];
  }

  return v;
}

template< class K >
class ffcalfunc    // to call the freefem function .. J, constraints, and derivatives
{
 public:
  Stack stack;
  Expression JJ, theparame;

  ffcalfunc(const ffcalfunc &f) : stack(f.stack), JJ(f.JJ), theparame(f.theparame) {}

  ffcalfunc(Stack s, Expression JJJ, Expression epar) : stack(s), JJ(JJJ), theparame(epar) {}

  K J(Rn_ x) const {
    KN< double > *p = GetAny< KN< double > * >((*theparame)(stack));
    *p = x;
    K ret = GetAny< K >((*JJ)(stack));
    WhereStackOfPtr2Free(stack)->clean( );
    return ret;
  }
};

typedef ffcalfunc< double > *ScalarFunc;
typedef ffcalfunc< Rn > *VectorFunc;
typedef ffcalfunc< Rnm > *MatrixFunc;

class GenericOptimizer {
 public:
  GenericOptimizer(nlopt::algorithm ALGO, int dim = 0)
    : opt(ALGO, dim), x(0), econsttol(0), iconsttol(0), econstrained(false), iconstrained(false),
      fit(0), d_fit(0), equaconst(0), d_equaconst(0), ineqconst(0), d_ineqconst(0), subopt(0) {}

  GenericOptimizer(nlopt::algorithm ALGO, const ffcalfunc< R > &_ff, Rn &xstart)
    : opt(ALGO, xstart.n), x(&xstart), econsttol(0), iconsttol(0), econstrained(false),
      iconstrained(false), fit(new ffcalfunc< R >(_ff)), d_fit(0), equaconst(0), d_equaconst(0),
      ineqconst(0), d_ineqconst(0), subopt(0) {
    opt.set_min_objective(NLoptFunc, static_cast< void * >(this));
  }

  virtual ~GenericOptimizer( ) {
    Clean(fit);
    Clean(d_fit);
    Clean(equaconst);
    Clean(d_equaconst);
    Clean(ineqconst);
    Clean(d_ineqconst);
    Clean(subopt);
  }

  double operator( )( ) {
    double minf;

    vector< double > vv(x->n);

    for (int i = 0; i < vv.size( ); ++i) {
      vv[i] = (*x)[i];
    }

    opt.optimize(vv, minf);

    for (int i = 0; i < vv.size( ); ++i) {
      (*x)[i] = vv[i];
    }

    return minf;
  }

  virtual bool DF( ) const { return true; }

  virtual bool SA( ) const { return false; }

  virtual const char *Name( ) const { return "Generic Algorithm"; }

  virtual nlopt::algorithm Tag( ) const = 0;
  GenericOptimizer &SetEqualityConstraintsTolerance(const Rn_ &val) {
    econsttol = val;
    return *this;
  }

  GenericOptimizer &SetInequalityConstraintsTolerance(const Rn_ &val) {
    iconsttol = val;
    return *this;
  }

  GenericOptimizer &SetSCXRelativeTolerance(const double val) {
    opt.set_xtol_rel(val);
    return *this;
  }    // SC = stopping criteria

  GenericOptimizer &SetSCXAbsoluteTolerance(const Rn_ &val) {
    opt.set_xtol_abs(KnToStdVect(val));
    return *this;
  }

  GenericOptimizer &SetLowerBounds(const Rn_ &lb) {
    opt.set_lower_bounds(KnToStdVect(lb));
    return *this;
  }

  GenericOptimizer &SetUpperBounds(const Rn_ &ub) {
    opt.set_upper_bounds(KnToStdVect(ub));
    return *this;
  }

  GenericOptimizer &SetSCStopFunctionValue(const double val) {
    opt.set_stopval(val);
    return *this;
  }

  GenericOptimizer &SetSCRelativeFunctionTolerance(const double val) {
    opt.set_ftol_rel(val);
    return *this;
  }

  GenericOptimizer &SetSCAbsoluteFunctionTolerance(const double val) {
    opt.set_ftol_abs(val);
    return *this;
  }

  GenericOptimizer &SetSCMaxFunctionEvaluations(const long val) {
    opt.set_maxeval(val);
    return *this;
  }

  GenericOptimizer &SetSCEllapsedTime(const double val) {
    opt.set_maxtime(val);
    return *this;
  }

  GenericOptimizer &SetPopulationSize(const int val) {
    opt.set_population(static_cast< unsigned >(val));
    return *this;
  }

  virtual GenericOptimizer &SetVectorStorage(const int val) {
    opt.set_vector_storage(static_cast< unsigned >(val));
    return *this;
  }

  GenericOptimizer &SetObjectiveFunctionGradient(const ffcalfunc< Rn > &g) {
    Clean(d_fit) = new ffcalfunc< Rn >(g);
    return *this;
  }

  GenericOptimizer &SetEqualityConstraintFunction(const ffcalfunc< Rn > &f) {
    Clean(equaconst) = new ffcalfunc< Rn >(f);
    return *this;
  }

  GenericOptimizer &SetEqualityConstraintGradient(const ffcalfunc< Rnm > &g) {
    Clean(d_equaconst) = new ffcalfunc< Rnm >(g);
    return *this;
  }

  GenericOptimizer &SetInequalityConstraintFunction(const ffcalfunc< Rn > &f) {
    Clean(ineqconst) = new ffcalfunc< Rn >(f);
    return *this;
  }

  GenericOptimizer &SetInequalityConstraintGradient(const ffcalfunc< Rnm > &g) {
    Clean(d_ineqconst) = new ffcalfunc< Rnm >(g);
    return *this;
  }

  GenericOptimizer &SetEqualityConstraints( ) {
    if (econstrained) {
      opt.remove_equality_constraints( );
    }

    Rn etestv = equaconst->J(*x);
    if (econsttol.n == 0) {
      econsttol.resize(etestv.n);
      econsttol = 1.e-12;
    } else {
      assert(econsttol.n == etestv.n);
    }

    opt.add_equality_mconstraint(NLoptECDF, static_cast< void * >(this), KnToStdVect(econsttol));
    econstrained = true;
    return *this;
  }

  GenericOptimizer &SetInequalityConstraints( ) {
    if (iconstrained) {
      opt.remove_inequality_constraints( );
    }

    Rn itestv = ineqconst->J(*x);
    if (iconsttol.n == 0) {
      iconsttol.resize(itestv.n);
      iconsttol = 1.e-12;
    } else {
      assert(iconsttol.n == itestv.n);
    }

    opt.add_inequality_mconstraint(NLoptICDF, static_cast< void * >(this), KnToStdVect(iconsttol));
    iconstrained = true;
    return *this;
  }

  static double NLoptFunc(const std::vector< double > &xx, std::vector< double > &grad,
                          void *data) {
    GenericOptimizer *pthis = static_cast< GenericOptimizer * >(data);
    int n = xx.size( );
    Rn X(n);

    for (int i = 0; i < n; ++i) {
      X[i] = xx[i];
    }

    if (grad.size( ) && pthis->d_fit) {
      Rn dJ = pthis->d_fit->J(X);

      for (int i = 0; i < n; ++i) {
        grad[i] = dJ[i];
      }
    }

    return pthis->fit->J(X);
  }

  static void NLoptECDF(unsigned m, double *result, unsigned n, const double *xx, double *grad,
                        void *data) {
    GenericOptimizer *pthis = static_cast< GenericOptimizer * >(data);

    if (pthis->equaconst) {
      Rn X(n);

      for (int k = 0; k < n; ++k) {
        X[k] = xx[k];
      }

      Rn Y = pthis->equaconst->J(X);
      assert(Y.n == m);

      for (int i = 0; i < m; ++i) {
        result[i] = Y[i];
      }

      if (grad) {
        assert(pthis->d_equaconst);
        Rnm dconst = pthis->d_equaconst->J(X);
        assert(dconst.N( ) == m && dconst.M( ) == n);

        for (int i = 0; i < m; ++i) {
          for (int j = 0; j < n; ++j) {
            grad[i * n + j] = dconst(i, j);
          }
        }
      }
    } else {
      for (int i = 0; i < m; ++i) {
        result[i] = 0;
        if (grad) {
          for (int j = 0; j < n; ++j) {
            grad[i * n + j] = 0;
          }
        }
      }
    }
  }

  static void NLoptICDF(unsigned m, double *result, unsigned n, const double *xx, double *grad,
                        void *data) {
    GenericOptimizer *pthis = static_cast< GenericOptimizer * >(data);

    if (pthis->ineqconst) {
      Rn X(n);

      for (int k = 0; k < n; ++k) {
        X[k] = xx[k];
      }

      Rn Y = pthis->ineqconst->J(X);
      assert(Y.N( ) == m);

      for (int i = 0; i < m; ++i) {
        result[i] = Y[i];
      }

      if (grad) {
        assert(pthis->d_ineqconst);
        Rnm dconst = pthis->d_ineqconst->J(X);
        assert(dconst.N( ) == m && dconst.M( ) == n);

        for (int i = 0; i < m; ++i) {
          for (int j = 0; j < n; ++j) {
            grad[i * n + j] = dconst(i, j);
          }
        }
      }
    } else {
      for (int i = 0; i < m; ++i) {
        result[i] = 0;
        if (grad) {
          for (int j = 0; j < n; ++j) {
            grad[i * n + j] = 0;
          }
        }
      }
    }
  }

  nlopt::opt opt, *subopt;
  Rn *x, econsttol, iconsttol;
  bool iconstrained, econstrained;
  ScalarFunc fit;
  VectorFunc d_fit, equaconst, ineqconst;
  MatrixFunc d_equaconst, d_ineqconst;

 private:
  GenericOptimizer( );
  template< class T >
  static T *&Clean(T *&p) {
    if (p) {
      delete p;
    }

    p = 0;
    return p;
  }
};

template< nlopt::algorithm ALGO >
class Optimizer : public GenericOptimizer {
 public:
  Optimizer(int dim = 0) : GenericOptimizer(ALGO, dim) {}

  Optimizer(const ffcalfunc< R > &_ff, Rn &xstart) : GenericOptimizer(ALGO, _ff, xstart) {}

  ~Optimizer( ) {}

  bool DF( ) const { return Info< ALGO >::DF; }

  bool SA( ) const { return Info< ALGO >::SA; }

  nlopt::algorithm Tag( ) const { return ALGO; }

  const char *Name( ) const { return Info< ALGO >::name; }

 private:
};

template< bool a >
struct MyCheck {
  MyCheck( ) {}
};

template<>
struct MyCheck< false > {
 private:
  MyCheck( ) {}
};

template< nlopt::algorithm ALGO >
class SAOptimizer : public GenericOptimizer {
 public:
  SAOptimizer(int dim = 0) : GenericOptimizer(ALGO, dim), subopt(0) {
    MyCheck< Info< ALGO >::SA >( );
  }

  SAOptimizer(const ffcalfunc< R > &_ff, Rn &xstart)
    : GenericOptimizer(ALGO, _ff, xstart), subopt(0) {
    MyCheck< Info< ALGO >::SA >( );
  }

  ~SAOptimizer( ) {
    if (subopt) {
      delete subopt;
    }

    subopt = 0;
  }

  GenericOptimizer *subopt;

  bool DF( ) const { return Info< ALGO >::DF && (subopt ? subopt->DF( ) : false); }

  bool SA( ) const { return true; }

  nlopt::algorithm Tag( ) const { return ALGO; }

  const char *Name( ) const { return Info< ALGO >::name; }

  GenericOptimizer &SetSubOptimizer(const string &name = string( ), bool save = 1);
  GenericOptimizer &SetSASCXRelativeTolerance(const double val) {
    if (subopt) {
      subopt->opt.set_xtol_rel(val);
    }

    return *this;
  }    // SC = stopping criteria

  GenericOptimizer &SetSASCXAbsoluteTolerance(const Rn_ &val) {
    if (subopt) {
      subopt->opt.set_xtol_abs(KnToStdVect(val));
    }

    return *this;
  }

  GenericOptimizer &SetSASCStopFunctionValue(const double val) {
    if (subopt) {
      subopt->opt.set_stopval(val);
    }

    return *this;
  }

  GenericOptimizer &SetSASCRelativeFunctionTolerance(const double val) {
    if (subopt) {
      subopt->opt.set_ftol_rel(val);
    }

    return *this;
  }

  GenericOptimizer &SetSASCAbsoluteFunctionTolerance(const double val) {
    if (subopt) {
      subopt->opt.set_ftol_abs(val);
    }

    return *this;
  }

  GenericOptimizer &SetSASCMaxFunctionEvaluations(const long val) {
    if (subopt) {
      subopt->opt.set_maxeval(val);
    }

    return *this;
  }

  GenericOptimizer &SetSASCEllapsedTime(const double val) {
    if (subopt) {
      subopt->opt.set_maxtime(val);
    }

    return *this;
  }

  GenericOptimizer &SetSAPopulationSize(const int val) {
    if (subopt) {
      subopt->opt.set_population(static_cast< unsigned >(val));
    }

    return *this;
  }

  GenericOptimizer &SetVectorStorage(const int val) {
    if (subopt) {
      subopt->opt.set_vector_storage(static_cast< unsigned >(val));
    }

    return *this;
  }
};

template< nlopt::algorithm ALGO >
GenericOptimizer &SAOptimizer< ALGO >::SetSubOptimizer(const string &name, bool save) {
  if (!subopt) {
    if (name == "DIRECT") {
      subopt = new Optimizer< nlopt::GN_DIRECT >(x->n);
    } else if (name == "DIRECTL") {
      subopt = new Optimizer< nlopt::GN_DIRECT_L >(x->n);
    } else if (name == "DIRECTLRand") {
      subopt = new Optimizer< nlopt::GN_DIRECT_L_RAND >(x->n);
    } else if (name == "DIRECTNoScal") {
      subopt = new Optimizer< nlopt::GN_DIRECT_NOSCAL >(x->n);
    } else if (name == "DIRECTLNoScal") {
      subopt = new Optimizer< nlopt::GN_DIRECT_L_NOSCAL >(x->n);
    } else if (name == "DIRECTLRandNoScal") {
      subopt = new Optimizer< nlopt::GN_DIRECT_L_RAND_NOSCAL >(x->n);
    } else if (name == "OrigDIRECT") {
      subopt = new Optimizer< nlopt::GN_ORIG_DIRECT >(x->n);
    } else if (name == "OrigDIRECTL") {
      subopt = new Optimizer< nlopt::GN_ORIG_DIRECT_L >(x->n);
    } else if (name == "StoGO") {
      subopt = new Optimizer< nlopt::GD_STOGO >(x->n);
    } else if (name == "StoGORand") {
      subopt = new Optimizer< nlopt::GD_STOGO_RAND >(x->n);
    } else if (name == "LBFGS") {
      subopt = new Optimizer< nlopt::LD_LBFGS >(x->n);
    } else if (name == "PRAXIS") {
      subopt = new Optimizer< nlopt::LN_PRAXIS >(x->n);
    } else if (name == "Var1") {
      subopt = new Optimizer< nlopt::LD_VAR1 >(x->n);
    } else if (name == "Var2") {
      subopt = new Optimizer< nlopt::LD_VAR2 >(x->n);
    } else if (name == "TNewton") {
      subopt = new Optimizer< nlopt::LD_TNEWTON >(x->n);
    } else if (name == "TNewtonRestart") {
      subopt = new Optimizer< nlopt::LD_TNEWTON_RESTART >(x->n);
    } else if (name == "TNewtonPrecond") {
      subopt = new Optimizer< nlopt::LD_TNEWTON_PRECOND >(x->n);
    } else if (name == "TNewtonPrecondRestart") {
      subopt = new Optimizer< nlopt::LD_TNEWTON_PRECOND_RESTART >(x->n);
    } else if (name == "CRS2") {
      subopt = new Optimizer< nlopt::GN_CRS2_LM >(x->n);
    } else if (name == "MMA") {
      subopt = new Optimizer< nlopt::LD_MMA >(x->n);
    } else if (name == "COBYLA") {
      subopt = new Optimizer< nlopt::LN_COBYLA >(x->n);
    } else if (name == "NEWUOA") {
      subopt = new Optimizer< nlopt::LN_NEWUOA >(x->n);
    } else if (name == "NEWUOABound") {
      subopt = new Optimizer< nlopt::LN_NEWUOA_BOUND >(x->n);
    } else if (name == "NelderMead") {
      subopt = new Optimizer< nlopt::LN_NELDERMEAD >(x->n);
    } else if (name == "Sbplx") {
      subopt = new Optimizer< nlopt::LN_SBPLX >(x->n);
    } else if (name == "BOBYQA") {
      subopt = new Optimizer< nlopt::LN_BOBYQA >(x->n);
    } else if (name == "ISRES") {
      subopt = new Optimizer< nlopt::GN_ISRES >(x->n);
    } else if (name == "SLSQP") {
      subopt = new Optimizer< nlopt::LD_SLSQP >(x->n);
    } else {
      cout << "Warning: unknown or unauthorized optimizer name passed as sub algorithm to "
           << Info< ALGO >::name << endl;
    }
  }

  if (subopt && save) {
    opt.set_local_optimizer(subopt->opt);
  }

  return *this;
}

template< nlopt::algorithm ALGO, bool SA = Info< ALGO >::SA >
class OptimNLopt : public OneOperator {
 public:
  const int cas;

  class E_NLopt : public E_F0mps {
   public:
    const int cas;
    static basicAC_F0::name_and_type name_param[];
    static const int n_name_param = 18;
    Expression nargs[n_name_param];
    Expression X;
    C_F0 inittheparam, theparam, closetheparam;
    Expression JJ;
    Expression GradJ, EIConst, EGradIConst, EEConst, EGradEConst;
    long arg(int i, Stack stack, long a) const {
      return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
    }

    R arg(int i, Stack stack, R a) const { return nargs[i] ? GetAny< R >((*nargs[i])(stack)) : a; }

    Rn_ arg(int i, Stack stack, Rn_ a) const {
      return nargs[i] ? GetAny< Rn_ >((*nargs[2])(stack)) : a;
    }

    template< typename T >
    T Arg(int i, Stack s) const {
      return GetAny< T >((*nargs[i])(s));
    }

    E_NLopt(const basicAC_F0 &args, int cc) : cas(cc) {
      int nbj = args.size( ) - 1;

      Block::open(currentblock);    // make a new block to
      X = to< Rn * >(args[nbj]);
      C_F0 X_n(args[nbj], "n");
      // the expression to init the theparam of all
      inittheparam =
        currentblock->NewVar< LocalVariable >("the parameter", atype< KN< R > * >( ), X_n);
      theparam = currentblock->Find("the parameter");    // the expression for the parameter
      args.SetNameParam(n_name_param, name_param, nargs);
      const Polymorphic *opJ = 0;
      if (nbj > 0) {
        opJ = dynamic_cast< const Polymorphic * >(args[0].LeftValue( ));
        assert(opJ);
      }

      JJ = to< R >(C_F0(opJ, "(", theparam));
      const Polymorphic *gradient = nargs[0] ? dynamic_cast< const Polymorphic * >(nargs[0]) : 0,
                        *iconst = nargs[1] ? dynamic_cast< const Polymorphic * >(nargs[1]) : 0,
                        *gradiconst = nargs[2] ? dynamic_cast< const Polymorphic * >(nargs[2]) : 0,
                        *econst = nargs[3] ? dynamic_cast< const Polymorphic * >(nargs[3]) : 0,
                        *gradeconst = nargs[4] ? dynamic_cast< const Polymorphic * >(nargs[4]) : 0;
      if (gradient) {
        GradJ = to< Rn_ >(C_F0(gradient, "(", theparam));
      }

      if (iconst) {
        EIConst = to< Rn_ >(C_F0(iconst, "(", theparam));
      }

      if (gradiconst) {
        EGradIConst = to< Rnm_ >(C_F0(gradiconst, "(", theparam));
      }

      if (econst) {
        EEConst = to< Rn_ >(C_F0(econst, "(", theparam));
      }

      if (gradeconst) {
        EGradEConst = to< Rnm_ >(C_F0(gradeconst, "(", theparam));
      }

      // closetheparam=currentblock->close(currentblock);   // the cleanning block expression
      closetheparam = C_F0((Expression)Block::snewclose(currentblock), atype< void >( ));
    }

    virtual AnyType operator( )(Stack stack) const {
      double cost = 1e100;

      WhereStackOfPtr2Free(stack) = new StackOfPtr2Free(stack);    // FH mars 2005
      Rn &x = *GetAny< Rn * >((*X)(stack));
      long n = x.N( );
      const bool gradient = nargs[0] ? dynamic_cast< const Polymorphic * >(nargs[0]) : 0,
                 iconst = nargs[1] ? dynamic_cast< const Polymorphic * >(nargs[1]) : 0,
                 gradiconst = nargs[2] ? dynamic_cast< const Polymorphic * >(nargs[2]) : 0,
                 econst = nargs[3] ? dynamic_cast< const Polymorphic * >(nargs[3]) : 0,
                 gradeconst = nargs[4] ? dynamic_cast< const Polymorphic * >(nargs[4]) : 0;
      long iprint = verbosity;
      ffcalfunc< double > ffJ(stack, JJ, theparam);

      Optimizer< ALGO > optim(ffJ, x);
      if (nargs[5]) {
        optim.SetLowerBounds(Arg< Rn_ >(5, stack));
      }

      if (nargs[6]) {
        optim.SetUpperBounds(Arg< Rn_ >(6, stack));
      }

      if (nargs[7]) {
        optim.SetSCStopFunctionValue(Arg< R >(7, stack));
      }

      if (nargs[8]) {
        optim.SetEqualityConstraintsTolerance(Arg< Rn_ >(8, stack));
      }

      if (nargs[9]) {
        optim.SetSCXRelativeTolerance(Arg< R >(9, stack));
      }

      if (nargs[10]) {
        optim.SetSCXAbsoluteTolerance(Arg< Rn_ >(10, stack));
      }

      if (nargs[11]) {
        optim.SetSCRelativeFunctionTolerance(Arg< R >(11, stack));
      }

      if (nargs[12]) {
        optim.SetSCAbsoluteFunctionTolerance(Arg< R >(12, stack));
      }

      if (nargs[13]) {
        optim.SetSCMaxFunctionEvaluations(Arg< long >(13, stack));
      }

      if (nargs[14]) {
        optim.SetSCEllapsedTime(Arg< R >(14, stack));
      }

      if (nargs[15]) {
        optim.SetInequalityConstraintsTolerance(Arg< Rn_ >(15, stack));
      }

      if (nargs[16]) {
        optim.SetPopulationSize(static_cast< int >(Arg< long >(16, stack)));
      }

      if (nargs[17]) {
        optim.SetVectorStorage(static_cast< int >(Arg< long >(17, stack)));
        if (optim.DF( )) {
          cout << "Warning: in " << optim.Name( )
               << " algorithm - using nGradStored is pointless (no gradient to store in a "
                  "derivative free context)."
               << endl;
        } else if (ALGO == nlopt::LD_SLSQP || ALGO == nlopt::LD_MMA) {
          cout << "Warning: nGradStored can't be used with " << optim.Name( )
               << ", parameter will be ignored." << endl;
        }
      }

      if (econst) {
        optim.SetEqualityConstraintFunction(ffcalfunc< Rn >(stack, EEConst, theparam));
      }

      if (iconst) {
        optim.SetInequalityConstraintFunction(ffcalfunc< Rn >(stack, EIConst, theparam));
      }

      if (optim.DF( )) {
        if (gradient) {
          cout
            << "Warning: in " << optim.Name( )
            << " algorithm - derivative free algorithm will ignore the objective function gradient."
            << endl;
        }

        if (gradiconst) {
          cout << "Warning: in " << optim.Name( )
               << " algorithm - derivative free algorithm will ignore the inequality constraints "
                  "gradient."
               << endl;
          if (!iconst) {
            cout << "Also note that this gradient has been provided for an inexisting set of "
                    "inequality constraints!"
                 << endl;
          }
        }

        if (gradeconst) {
          cout << "Warning: in " << optim.Name( )
               << " algorithm - derivative free algorithm will ignore the equality constraints "
                  "gradient."
               << endl;
          if (!econst) {
            cout << "Also note that this gradient has been provided for an inexisting set of "
                    "equality constraints!"
                 << endl;
          }
        }
      } else {
        if (gradient) {
          optim.SetObjectiveFunctionGradient(ffcalfunc< Rn >(stack, GradJ, theparam));
        } else {
          cout << "Warning: in " << optim.Name( )
               << " algorithm - no objective function gradient has been provided (choose a "
                  "derivative free algorithm if it is not available)."
               << endl;
        }

        if (econst) {
          if (gradeconst) {
            optim.SetEqualityConstraintGradient(ffcalfunc< Rnm >(stack, EGradEConst, theparam));
          } else {
            cout << "Warning: in " << optim.Name( )
                 << " algorithm - no equality constraints gradients has been provided." << endl;
          }
        } else if (gradeconst) {
          cout << "Warning: in " << optim.Name( )
               << " algorithm - gradients have been provided for an inexisting set of equality "
                  "constraints."
               << endl;
        }

        if (iconst) {
          if (gradiconst) {
            optim.SetInequalityConstraintGradient(ffcalfunc< Rnm >(stack, EGradIConst, theparam));
          } else {
            cout << "Warning: in " << optim.Name( )
                 << " algorithm - no inequality constraints gradients has been provided." << endl;
          }
        } else if (gradiconst) {
          cout << "Warning: in " << optim.Name( )
               << " algorithm - gradients have been provided for an inexisting set of inequality "
                  "constraints."
               << endl;
        }
      }

      if (econst) {
        optim.SetEqualityConstraints( );
      }

      if (iconst) {
        optim.SetInequalityConstraints( );
      }

      if (verbosity > 1) {
        cout << Info< ALGO >::name << " starting..." << endl;
      }

      try {
        cost = optim( );
      } catch (const nlopt::roundoff_limited &) {
        cout << " nlopt roundoff limited" << endl;
      } catch (const nlopt::forced_stop &) {
        cout << " nlopt forced stop" << endl;
      } catch (const std::runtime_error &) {
        cout << "runtime error" << endl;
      } catch (const std::invalid_argument &) {
        cout << "invalid argument" << endl;
      } catch (const std::bad_alloc &) {
        cout << "bad alloc" << endl;
      }

      closetheparam.eval(stack);                // clean memory
      WhereStackOfPtr2Free(stack)->clean( );    // FH mars 2005
      return cost;                              // SetAny<long>(0);  Modif FH  july 2005
    }

    operator aType( ) const { return atype< double >( ); }
  };

  E_F0 *code(const basicAC_F0 &args) const { return new E_NLopt(args, cas); }

  OptimNLopt(int c)
    : OneOperator(atype< double >( ), atype< Polymorphic * >( ), atype< KN< R > * >( )), cas(c) {}
};

template< nlopt::algorithm ALGO, bool SA >
basicAC_F0::name_and_type OptimNLopt< ALGO, SA >::E_NLopt::name_param[] = {
  {"grad", &typeid(Polymorphic *)},
  {"IConst", &typeid(Polymorphic *)},
  {"gradIConst", &typeid(Polymorphic *)},
  {"EConst", &typeid(Polymorphic *)},
  {"gradEConst", &typeid(Polymorphic *)},
  {"lb", &typeid(KN_< double >)},
  {"ub", &typeid(KN_< double >)},
  {"stopFuncValue", &typeid(double)},
  {"tolEConst", &typeid(KN_< double >)},
  {"stopRelXTol", &typeid(double)},
  {"stopAbsXTol", &typeid(KN_< double >)},
  {"stopRelFTol", &typeid(double)},
  {"stopAbsFTol", &typeid(double)},
  {"stopMaxFEval", &typeid(long)},
  {"stopTime", &typeid(double)},
  {"tolIConst", &typeid(KN_< double >)},
  {"popSize", &typeid(long)},
  {"nGradStored", &typeid(long)}};

template< nlopt::algorithm ALGO >
class OptimNLopt< ALGO, true > : public OneOperator {
 public:
  const int cas;

  class E_NLopt : public E_F0mps {
   public:
    const int cas;
    static basicAC_F0::name_and_type name_param[];
    static const int n_name_param = 27;
    Expression nargs[n_name_param];
    Expression X;
    C_F0 inittheparam, theparam, closetheparam;
    Expression JJ;
    Expression GradJ, EIConst, EGradIConst, EEConst, EGradEConst;
    long arg(int i, Stack stack, long a) const {
      return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
    }

    R arg(int i, Stack stack, R a) const { return nargs[i] ? GetAny< R >((*nargs[i])(stack)) : a; }

    Rn_ arg(int i, Stack stack, Rn_ a) const {
      return nargs[i] ? GetAny< Rn_ >((*nargs[2])(stack)) : a;
    }

    template< typename T >
    T Arg(int i, Stack s) const {
      return GetAny< T >((*nargs[i])(s));
    }

    E_NLopt(const basicAC_F0 &args, int cc) : cas(cc) {
      int nbj = args.size( ) - 1;

      Block::open(currentblock);    // make a new block to
      X = to< Rn * >(args[nbj]);
      C_F0 X_n(args[nbj], "n");
      // the expression to init the theparam of all
      inittheparam =
        currentblock->NewVar< LocalVariable >("the parameter", atype< KN< R > * >( ), X_n);
      theparam = currentblock->Find("the parameter");    // the expression for the parameter
      args.SetNameParam(n_name_param, name_param, nargs);
      const Polymorphic *opJ = 0;
      if (nbj > 0) {
        opJ = dynamic_cast< const Polymorphic * >(args[0].LeftValue( ));
        assert(opJ);
      }

      JJ = to< R >(C_F0(opJ, "(", theparam));
      const Polymorphic *gradient = nargs[0] ? dynamic_cast< const Polymorphic * >(nargs[0]) : 0,
                        *iconst = nargs[1] ? dynamic_cast< const Polymorphic * >(nargs[1]) : 0,
                        *gradiconst = nargs[2] ? dynamic_cast< const Polymorphic * >(nargs[2]) : 0,
                        *econst = nargs[3] ? dynamic_cast< const Polymorphic * >(nargs[3]) : 0,
                        *gradeconst = nargs[4] ? dynamic_cast< const Polymorphic * >(nargs[4]) : 0;
      if (gradient) {
        GradJ = to< Rn_ >(C_F0(gradient, "(", theparam));
      }

      if (iconst) {
        EIConst = to< Rn_ >(C_F0(iconst, "(", theparam));
      }

      if (gradiconst) {
        EGradIConst = to< Rnm_ >(C_F0(gradiconst, "(", theparam));
      }

      if (econst) {
        EEConst = to< Rn_ >(C_F0(econst, "(", theparam));
      }

      if (gradeconst) {
        EGradEConst = to< Rnm_ >(C_F0(gradeconst, "(", theparam));
      }

      closetheparam = C_F0((Expression)Block::snewclose(currentblock), atype< void >( ));
    }

    virtual AnyType operator( )(Stack stack) const {
      double cost = 1e100;

      WhereStackOfPtr2Free(stack) = new StackOfPtr2Free(stack);    // FH mars 2005
      Rn &x = *GetAny< Rn * >((*X)(stack));
      long n = x.N( );
      const bool gradient = nargs[0] ? dynamic_cast< const Polymorphic * >(nargs[0]) : 0,
                 iconst = nargs[1] ? dynamic_cast< const Polymorphic * >(nargs[1]) : 0,
                 gradiconst = nargs[2] ? dynamic_cast< const Polymorphic * >(nargs[2]) : 0,
                 econst = nargs[3] ? dynamic_cast< const Polymorphic * >(nargs[3]) : 0,
                 gradeconst = nargs[4] ? dynamic_cast< const Polymorphic * >(nargs[4]) : 0;
      long iprint = verbosity;
      ffcalfunc< double > ffJ(stack, JJ, theparam);

      SAOptimizer< ALGO > optim(ffJ, x);
      if (nargs[5]) {
        optim.SetLowerBounds(Arg< Rn_ >(5, stack));
      }

      if (nargs[6]) {
        optim.SetUpperBounds(Arg< Rn_ >(6, stack));
      }

      if (nargs[7]) {
        optim.SetSCStopFunctionValue(Arg< R >(7, stack));
      }

      if (nargs[8]) {
        optim.SetEqualityConstraintsTolerance(Arg< Rn_ >(8, stack));
      }

      if (nargs[9]) {
        optim.SetSCXRelativeTolerance(Arg< R >(9, stack));
      }

      if (nargs[10]) {
        optim.SetSCXAbsoluteTolerance(Arg< Rn_ >(10, stack));
      }

      if (nargs[11]) {
        optim.SetSCRelativeFunctionTolerance(Arg< R >(11, stack));
      }

      if (nargs[12]) {
        optim.SetSCAbsoluteFunctionTolerance(Arg< R >(12, stack));
      }

      if (nargs[13]) {
        optim.SetSCMaxFunctionEvaluations(Arg< long >(13, stack));
      }

      if (nargs[14]) {
        optim.SetSCEllapsedTime(Arg< R >(14, stack));
      }

      if (nargs[15]) {
        optim.SetInequalityConstraintsTolerance(Arg< Rn_ >(15, stack));
      }

      if (nargs[16]) {
        optim.SetPopulationSize(static_cast< int >(Arg< long >(16, stack)));
      }

      if (nargs[17]) {
        optim.SetSubOptimizer(*Arg< string * >(17, stack), 0);
      } else {
        cout << "Warning: in " << optim.Name( )
             << " algorithm - you have to specify a local optimizer, aboarting optimization (use "
                "the subOpt named parameter)."
             << endl;
      }

      if (nargs[18]) {
        optim.SetSASCStopFunctionValue(Arg< R >(18, stack));
      }

      if (nargs[19]) {
        optim.SetSASCXRelativeTolerance(Arg< R >(19, stack));
      }

      if (nargs[20]) {
        optim.SetSASCXAbsoluteTolerance(Arg< Rn_ >(20, stack));
      }

      if (nargs[21]) {
        optim.SetSASCRelativeFunctionTolerance(Arg< R >(21, stack));
      }

      if (nargs[22]) {
        optim.SetSASCAbsoluteFunctionTolerance(Arg< R >(22, stack));
      }

      if (nargs[23]) {
        optim.SetSASCMaxFunctionEvaluations(Arg< long >(23, stack));
      }

      if (nargs[24]) {
        optim.SetSASCEllapsedTime(Arg< R >(24, stack));
      }

      if (nargs[25]) {
        optim.SetSAPopulationSize(static_cast< int >(Arg< long >(25, stack)));
      }

      if (nargs[26]) {
        optim.SetVectorStorage(static_cast< int >(Arg< long >(26, stack)));
        if (optim.DF( )) {
          cout << "Warning: in " << optim.subopt->Name( )
               << " algorithm - using nGradStored is pointless (no gradient to store in a "
                  "derivative free context)."
               << endl;
        } else if (optim.subopt->Tag( ) == nlopt::LD_SLSQP ||
                   optim.subopt->Tag( ) == nlopt::LD_MMA) {
          cout << "Warning: nGradStored can't be used with " << optim.Name( )
               << ", parameter will be ignored." << endl;
        }
      }

      optim.SetSubOptimizer( );

      if (econst) {
        optim.SetEqualityConstraintFunction(ffcalfunc< Rn >(stack, EEConst, theparam));
      }

      if (iconst) {
        optim.SetInequalityConstraintFunction(ffcalfunc< Rn >(stack, EIConst, theparam));
      }

      if (optim.subopt) {
        if (optim.subopt->DF( )) {
          if (gradient) {
            cout << "Warning: in " << optim.Name( )
                 << " algorithm - derivative free sub-algorithm will ignore the objective function "
                    "gradient."
                 << endl;
          }

          if (gradiconst) {
            cout << "Warning: in " << optim.Name( )
                 << " algorithm - derivative free sub-algorithm will ignore the inequality "
                    "constraints gradient."
                 << endl;
            if (!iconst) {
              cout << "Also note that this gradient has been provided for an inexisting set of "
                      "inequality constraints!"
                   << endl;
            }
          }

          if (gradeconst) {
            cout << "Warning: in " << optim.Name( )
                 << " algorithm - derivative free sub-algorithm will ignore the equality "
                    "constraints gradient."
                 << endl;
            if (!econst) {
              cout << "Also note that this gradient has been provided for an inexisting set of "
                      "equality constraints!"
                   << endl;
            }
          }
        } else {
          if (gradient) {
            optim.SetObjectiveFunctionGradient(ffcalfunc< Rn >(stack, GradJ, theparam));
          } else {
            cout << "Warning: in " << optim.Name( )
                 << " algorithm - no objective function gradient has been provided (choose a "
                    "derivative free local search if it is not available)."
                 << endl;
          }

          if (econst) {
            if (gradeconst) {
              optim.SetEqualityConstraintGradient(ffcalfunc< Rnm >(stack, EGradEConst, theparam));
            } else {
              cout << "Warning: in " << optim.Name( )
                   << " algorithm - no equality constraints gradients has been provided." << endl;
            }
          } else if (gradeconst) {
            cout << "Warning: in " << optim.Name( )
                 << " algorithm - gradients have been provided for an inexisting set of equality "
                    "constraints."
                 << endl;
          }

          if (iconst) {
            if (gradiconst) {
              optim.SetInequalityConstraintGradient(ffcalfunc< Rnm >(stack, EGradIConst, theparam));
            } else {
              cout << "Warning: in " << optim.Name( )
                   << " algorithm - no inequality constraints gradients has been provided." << endl;
            }
          } else if (gradiconst) {
            cout << "Warning: in " << optim.Name( )
                 << " algorithm - gradients have been provided for an inexisting set of inequality "
                    "constraints."
                 << endl;
          }
        }

        if (econst) {
          optim.SetEqualityConstraints( );
        }

        if (iconst) {
          optim.SetInequalityConstraints( );
        }

        if (verbosity > 1) {
          cout << Info< ALGO >::name << " starting..." << endl;
        }

        try {
          cost = optim( );
        } catch (const nlopt::roundoff_limited &) {
          cout << " nlopt roundoff limited" << endl;
        } catch (const nlopt::forced_stop &) {
          cout << " nlopt forced stop" << endl;
        } catch (const std::runtime_error &) {
          cout << "runtime error" << endl;
        } catch (const std::invalid_argument &) {
          cout << "invalid argument" << endl;
        } catch (const std::bad_alloc &) {
          cout << "bad alloc" << endl;
        }
      }

      closetheparam.eval(stack);                // clean memory
      WhereStackOfPtr2Free(stack)->clean( );    // FH mars 2005
      return cost;                              // SetAny<long>(0);  Modif FH  july 2005
    }

    operator aType( ) const { return atype< double >( ); }
  };

  E_F0 *code(const basicAC_F0 &args) const { return new E_NLopt(args, cas); }

  OptimNLopt(int c)
    : OneOperator(atype< double >( ), atype< Polymorphic * >( ), atype< KN< R > * >( )), cas(c) {}
};

template< nlopt::algorithm ALGO >
basicAC_F0::name_and_type OptimNLopt< ALGO, true >::E_NLopt::name_param[] = {
  {"grad", &typeid(Polymorphic *)},
  {"IConst", &typeid(Polymorphic *)},
  {"gradIConst", &typeid(Polymorphic *)},
  {"EConst", &typeid(Polymorphic *)},
  {"gradEConst", &typeid(Polymorphic *)},
  {"lb", &typeid(KN_< double >)},
  {"ub", &typeid(KN_< double >)},
  {"stopFuncValue", &typeid(double)},
  {"tolEConst", &typeid(KN_< double >)},
  {"stopRelXTol", &typeid(double)},
  {"stopAbsXTol", &typeid(KN_< double >)},
  {"stopRelFTol", &typeid(double)},
  {"stopAbsFTol", &typeid(double)},
  {"stopMaxFEval", &typeid(long)},
  {"stopTime", &typeid(double)},
  {"tolIConst", &typeid(KN_< double >)},
  {"popSize", &typeid(long)},    // 16
  {"subOpt", &typeid(string *)},
  {"SOStopFuncValue", &typeid(double)},
  {"SOStopRelXTol", &typeid(double)},
  {"SOStopAbsXTol", &typeid(KN_< double >)},
  {"SOStopRelFTol", &typeid(double)},
  {"SOStopAbsFTol", &typeid(double)},
  {"SOStopMaxFEval", &typeid(long)},
  {"SOStopTime", &typeid(double)},
  {"SOPopSize", &typeid(long)},
  {"nGradStored", &typeid(long)}};

/*  class Init { public:
 * Init();
 * };
 *
 * $1 */

static void Load_Init( ) {
  Global.Add("nloptDIRECT", "(", new OptimNLopt< nlopt::GN_DIRECT >(1));
  Global.Add("nloptDIRECTL", "(", new OptimNLopt< nlopt::GN_DIRECT_L >(1));
  Global.Add("nloptDIRECTLRand", "(", new OptimNLopt< nlopt::GN_DIRECT_L_RAND >(1));
  Global.Add("nloptDIRECTNoScal", "(", new OptimNLopt< nlopt::GN_DIRECT_NOSCAL >(1));
  Global.Add("nloptDIRECTLNoScal", "(", new OptimNLopt< nlopt::GN_DIRECT_L_NOSCAL >(1));
  Global.Add("nloptDIRECTLRandNoScal", "(", new OptimNLopt< nlopt::GN_DIRECT_L_RAND_NOSCAL >(1));
  Global.Add("nloptOrigDIRECT", "(", new OptimNLopt< nlopt::GN_ORIG_DIRECT >(1));
  Global.Add("nloptOrigDIRECTL", "(", new OptimNLopt< nlopt::GN_ORIG_DIRECT_L >(1));
  Global.Add("nloptStoGO", "(", new OptimNLopt< nlopt::GD_STOGO >(1));
  Global.Add("nloptStoGORand", "(", new OptimNLopt< nlopt::GD_STOGO_RAND >(1));
  // Global.Add("nloptLBFGSNocedal",					"(",new
  // OptimNLopt<nlopt::LD_LBFGS_NOCEDAL>(1));
  // //Invalid argument
  Global.Add("nloptLBFGS", "(", new OptimNLopt< nlopt::LD_LBFGS >(1));
  Global.Add("nloptPRAXIS", "(", new OptimNLopt< nlopt::LN_PRAXIS >(1));
  Global.Add("nloptVar1", "(", new OptimNLopt< nlopt::LD_VAR1 >(1));
  Global.Add("nloptVar2", "(", new OptimNLopt< nlopt::LD_VAR2 >(1));
  Global.Add("nloptTNewton", "(", new OptimNLopt< nlopt::LD_TNEWTON >(1));
  Global.Add("nloptTNewtonRestart", "(", new OptimNLopt< nlopt::LD_TNEWTON_RESTART >(1));
  Global.Add("nloptTNewtonPrecond", "(", new OptimNLopt< nlopt::LD_TNEWTON_PRECOND >(1));
  Global.Add("nloptTNewtonPrecondRestart", "(",
             new OptimNLopt< nlopt::LD_TNEWTON_PRECOND_RESTART >(1));
  Global.Add("nloptCRS2", "(", new OptimNLopt< nlopt::GN_CRS2_LM >(1));
  Global.Add("nloptMMA", "(", new OptimNLopt< nlopt::LD_MMA >(1));
  Global.Add("nloptCOBYLA", "(", new OptimNLopt< nlopt::LN_COBYLA >(1));
  Global.Add("nloptNEWUOA", "(", new OptimNLopt< nlopt::LN_NEWUOA >(1));
  Global.Add("nloptNEWUOABound", "(", new OptimNLopt< nlopt::LN_NEWUOA_BOUND >(1));
  Global.Add("nloptNelderMead", "(", new OptimNLopt< nlopt::LN_NELDERMEAD >(1));
  Global.Add("nloptSbplx", "(", new OptimNLopt< nlopt::LN_SBPLX >(1));
  Global.Add("nloptBOBYQA", "(", new OptimNLopt< nlopt::LN_BOBYQA >(1));
  Global.Add("nloptISRES", "(", new OptimNLopt< nlopt::GN_ISRES >(1));
  Global.Add("nloptSLSQP", "(", new OptimNLopt< nlopt::LD_SLSQP >(1));
  Global.Add("nloptMLSL", "(", new OptimNLopt< nlopt::G_MLSL >(1));
  Global.Add("nloptMLSLLDS", "(", new OptimNLopt< nlopt::G_MLSL_LDS >(1));
  Global.Add("nloptAUGLAG", "(", new OptimNLopt< nlopt::AUGLAG >(1));
  Global.Add("nloptAUGLAGEQ", "(", new OptimNLopt< nlopt::AUGLAG_EQ >(1));
}

LOADFUNC(Load_Init)
