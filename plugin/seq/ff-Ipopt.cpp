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
// AUTHORS : Sylvain Auliac
// E-MAIL  : auliac@ann.jussieu.fr

// *INDENT-OFF* //
// ff-c++-LIBRARY-dep: ipopt mumps_seq blas libseq fc
// ff-c++-cpp-dep:
// *INDENT-ON* //

// using namespace std;

// TODO: remove this block as soon as autoconf is removed from FreeFem++
#ifndef CMAKE
#include "../../config.h"
#endif

#include "coin/IpTNLP.hpp"
#include "coin/IpIpoptApplication.hpp"
#include "ff++.hpp"

extern Block *currentblock;

typedef double R;
typedef KN_< R > Rn_;
typedef KN< R > Rn;
typedef KNM_< R > Rnm_;
typedef KNM< R > Rnm;

/*****************************************************************************************************************************
 *	Some misc. function usefull later...
 *****************************************************************************************************************************/

// A variadic function to add an undefinite number of elements to a set of short int
// This is used to define the set of named parameter which are not used when some assumptions
// upon the optimization poblem functions are met
void AddElements(std::set< unsigned short > &_set, int amount, int first, ...) {
  int elem = 0;
  va_list vl;

  va_start(vl, first);
  _set.insert(first);

  for (int i = 1; i < amount; i++) {
    elem = va_arg(vl, int);
    _set.insert(elem);
  }

  va_end(vl);
}

// A raw pointer cleaner
template< class T >
inline void clean(T *p) {
  if (p) {
    delete p;
    p = 0;
  }
}

// Pair compare (certainly already implemented in the STL with KeyLess...)
inline bool operator<=(const std::pair< int, int > &l, const std::pair< int, int > &r) {
  return (l.first < r.first) || (l.first == r.first && l.second <= r.second);
}

// Some logical operators (exclussive or and its negation)
inline bool XOR(bool a, bool b) { return (!a && b) || (a && !b); }

inline bool NXOR(bool a, bool b) { return !XOR(a, b); }

// A debug tool
#ifdef DEBUG
inline void SONDE( ) {
  static int i = 1;
  cout << "SONDE " << i << endl;
  ++i;
}

#else
inline void SONDE( ) {}

#endif

/*****************************************************************************************************************************
 *	FreeFem function callers
 *  ffcalfunc : template abstract mother class with a pointer to the freefem stack and the J virtual
 *method which computes the function
 *****************************************************************************************************************************/
template< class K >
class ffcalfunc {
 public:
  Stack stack;
  ffcalfunc(const ffcalfunc &f) : stack(f.stack) {}

  ffcalfunc(Stack _stack) : stack(_stack) {}

  virtual K J(Rn_) const = 0;
  virtual ~ffcalfunc( ) {}
};

/*****************************************************************************************************************************
 *	GeneralFunc : Most general case (specialized for sparse matrix returning functions, because
 *IPOPT need the hessian func to take some additional parameters).
 *    @theparame: ff expression of the parameter of the ff function, computing J(x) need the
 *associated KN to be set to the values of x
 *    @JJ       : ff expression of the function
 *****************************************************************************************************************************/
template< class K >
class GeneralFunc : public ffcalfunc< K > {
 public:
  Expression JJ, theparame;
  GeneralFunc(const GeneralFunc &f) : ffcalfunc< K >(f), JJ(f.JJ), theparame(f.theparame) {}

  GeneralFunc(Stack s, Expression JJJ, Expression epar)
    : ffcalfunc< K >(s), JJ(JJJ), theparame(epar) {}

  K J(Rn_ x) const {
    KN< double > *p = GetAny< KN< double > * >((*theparame)(this->stack));
    *p = x;
    K ret = GetAny< K >((*JJ)(this->stack));
    WhereStackOfPtr2Free(this->stack)->clean( );
    return ret;
  }
};

/*****************************************************************************************************************************
 *	P2ScalarFunc: encapsulate a function which is the sum of a bilinear and a linear form (no
 *constant part since it will be used as fitness function). It also handles the case of pure
 *quadratic or linear forms.
 *    @vf       : If true J will compute 0.5xMx - bx (x is the solution of Mx = b in the
 *unconstrained optimization process) if false xMx + bx is returned
 *    @M        : FF expression of the matrix of the bilinear form (null pointer for linear form
 *case)
 *    @b        : FF expression of the vector representation of the linear part (null for pure
 *quadratic case)
 *****************************************************************************************************************************/
class P2ScalarFunc : public ffcalfunc< R > {
 public:
  const bool vf;
  Expression M, b;    // Matrix of the quadratic part, vector of the linear part
  P2ScalarFunc(const P2ScalarFunc &f) : ffcalfunc< R >(f), M(f.M), b(f.b), vf(f.vf) {}

  P2ScalarFunc(Stack s, Expression _M, Expression _b, bool _vf = false)
    : ffcalfunc< R >(s), M(_M), b(_b), vf(_vf) {}

  R J(Rn_ x) const {
    Rn tmp(x.N( ), 0.);

    if (M) {
      Matrice_Creuse< R > *a = GetAny< Matrice_Creuse< R > * >((*M)(stack));
      MatriceMorse< R > *A = dynamic_cast< MatriceMorse< R > * >(&(*a->A));
      assert(A);
      tmp = (*A) * x;
      if (vf) {
        tmp /= 2.;
      }
    }

    if (b) {
      Rn *B = GetAny< Rn * >((*b)(stack));
      tmp += *B;
    }

    R res = 0.;

    for (int i = 0; i < x.N( ); ++i) {
      res += x[i] * tmp[i];
    }

    return res;
  }
};

/*****************************************************************************************************************************
 *	P1VectorFunc: encapsulate a function which is the sum of a linear part and a constant,
 *mostly used for affine/linear constraints, or for P2 fitness function gradient
 *    @vf       : Set to true if this is expected the gradient of a P2 scalar function associated to
 *Ax=b linear system J will then return Ax - b. Otherwize Ax+b is returned.
 *    @M        : FF expression of the matrix of the linear part
 *    @b        : FF expression of the vector representation of the constant part
 *****************************************************************************************************************************/
class P1VectorFunc : public ffcalfunc< Rn > {
 public:
  const bool vf;
  Expression M, b;
  P1VectorFunc(const P1VectorFunc &f) : ffcalfunc< Rn >(f), M(f.M), b(f.b), vf(f.vf) {}

  P1VectorFunc(Stack s, Expression _M, Expression _b, bool _vf = false)
    : ffcalfunc< Rn >(s), M(_M), b(_b), vf(_vf) {}

  Rn J(Rn_ x) const {
    Rn tmp(0);

    if (M) {
      Matrice_Creuse< R > *a = GetAny< Matrice_Creuse< R > * >((*M)(stack));
      MatriceMorse< R > *A = dynamic_cast< MatriceMorse< R > * >(&(*a->A));
      assert(A);
      if (tmp.N( ) != A->n) {
        tmp.resize(A->n);
        tmp = 0.;
      }

      tmp = (*A) * x;
    }

    if (b) {
      Rn *B = GetAny< Rn * >((*b)(stack));
      if (tmp.N( ) != B->N( )) {
        tmp.resize(B->N( ));
        tmp = 0.;
      }

      tmp += *B;
    }

    return tmp;
  }
};

/*****************************************************************************************************************************
 *	ffcalfunc<Matrice_Creuse>R>*>: specialization for sparse matrix returning function. When it
 *encapsulates the hessian function of the lagragian, non-linear constraints will need the
 *additional obj_factor and lagrange multiplier parameters. The one parameter version of J is called
 *if there is no non-linear constraints or if the objects represents the jacobian of the
 *constraints.
 *****************************************************************************************************************************/
template<>
class ffcalfunc< Matrice_Creuse< R > * > {
 public:
  typedef Matrice_Creuse< R > *K;
  Stack stack;
  ffcalfunc(const ffcalfunc &f) : stack(f.stack) {}

  ffcalfunc(Stack s) : stack(s) {}

  virtual K J(Rn_) const = 0;
  virtual K J(Rn_, double, Rn_) const = 0;
  virtual bool NLCHPEnabled( ) const = 0;    // Non Linear Constraints Hessian Prototype
  virtual ~ffcalfunc( ) {}
};

/*****************************************************************************************************************************
 *	GeneralSparseMatFunc: general case of sparse matrix returning function. Members datas added
 *are ff expression of the scalar objective factor and vectorial lagrange multipliers.
 *****************************************************************************************************************************/
class GeneralSparseMatFunc : public ffcalfunc< Matrice_Creuse< R > * > {
 private:
  typedef ffcalfunc< Matrice_Creuse< R > * > FFF;

 public:
  Expression JJ, param, paramlm, paramof;
  GeneralSparseMatFunc(const GeneralSparseMatFunc &f)
    : FFF(f), JJ(f.JJ), param(f.param), paramlm(f.paramlm), paramof(f.paramof){};
  GeneralSparseMatFunc(Stack s, Expression JJJ, Expression epar, Expression eparof = 0,
                       Expression eparlm = 0)
    : FFF(s), JJ(JJJ), param(epar), paramlm(eparlm), paramof(eparof) {
    ffassert(NXOR(paramlm, paramof));
  }

  bool NLCHPEnabled( ) const { return paramlm && paramof; }

  K J(Rn_ x) const {
    KN< double > *p = GetAny< KN< double > * >((*param)(stack));
    *p = x;
    K ret = GetAny< K >((*JJ)(stack));
    // cout << "call to ffcalfunc.J with " << *p << " and ret=" << ret << endl;
    WhereStackOfPtr2Free(stack)->clean( );
    return ret;
  }

  K J(Rn_ x, double of, Rn_ lm) const {
    if (paramlm && paramof) {
      KN< double > *p = GetAny< KN< double > * >((*param)(stack));
      double *pof = GetAny< double * >((*paramof)(stack));
      KN< double > *plm = GetAny< KN< double > * >((*paramlm)(stack));
      *p = x;
      *pof = of;
      int m = lm.N( ), mm = plm->N( );
      if ((m != mm) && mm) {
        cout << " ff-ipopt H : big bug int size ???" << m << " != " << mm << endl;
        abort( );
      }

      ;
      *plm = lm;
      K ret = GetAny< K >((*JJ)(stack));
      // cout << "call to ffcalfunc.J with " << *p << " and ret=" << ret << endl;
      WhereStackOfPtr2Free(stack)->clean( );
      return ret;
    } else {
      return J(x);
    }
  }
};

/*****************************************************************************************************************************
 *	ConstantSparseMatFunc: Encapsulate a constant matrix returning function. Just contains the
 *ff expression of the matrix (and stack inherited from mother class), this matrix is returned
 *regardless of x.
 *****************************************************************************************************************************/
class ConstantSparseMatFunc : public ffcalfunc< Matrice_Creuse< R > * > {
 private:
  typedef ffcalfunc< Matrice_Creuse< R > * > FFF;

 public:
  Expression M;    // Expression of the matrix
  ConstantSparseMatFunc(const ConstantSparseMatFunc &f) : FFF(f), M(f.M) {}

  ConstantSparseMatFunc(Stack s, Expression _M) : FFF(s), M(_M) {}

  bool NLCHPEnabled( ) const { return false; }

  K J(Rn_) const {
    K ret = M ? GetAny< K >((*M)(stack)) : 0;

    WhereStackOfPtr2Free(stack)->clean( );
    return ret;
  }

  K J(Rn_ x, double, Rn_) const { return J(x); }
};

typedef ffcalfunc< double > ScalarFunc;
typedef ffcalfunc< Rn > VectorFunc;
typedef ffcalfunc< Rnm > FullMatrixFunc;
typedef ffcalfunc< Matrice_Creuse< R > * > SparseMatFunc;

/*****************************************************************************************************************************
 *	SparseMatStructure: a class for sparse matrix structure management (mostly merging). The
 *most interesting methods in this class are : AddMatrix : merge the structure of the given matrix
 *to the structure of current object AddArrays : merge structure in arrays form to the current
 *object ToKn      : allocate the raws and cols pointers and fill them with the std::set<Z2> form of
 *the structure structure is then emptied if this method is passed a true value
 * ==> update 28/03/2012, autostruct proved useless since the structure merging can be done with
 *operator + (I did not no whether nullify coefficients where removed from the result but it
 *actually doesn't so the structure of the lagrangian hessian can be guessed exactly by evaluating
 *on a point yeilding the biggest fitness function hessian along with a dual vector filled with 1).
 *****************************************************************************************************************************/
class SparseMatStructure {
 public:
  typedef std::pair< int, int > Z2;
  typedef std::set< Z2 > Structure;
  typedef std::pair< KN< int >, KN< int > > Zn2;
  typedef Structure::const_iterator const_iterator;
  typedef Structure::iterator iterator;

  SparseMatStructure(bool _sym = 0) : structure( ), sym(_sym), n(0), m(0), raws(0), cols(0) {}

  SparseMatStructure(Matrice_Creuse< R > *M, bool _sym = 0)
    : structure( ), sym(_sym), n(M->N( )), m(M->M( )), raws(0), cols(0) {
    this->AddMatrix(M);
  }

  template< class INT >
  SparseMatStructure(const KN< INT > &I, const KN< INT > &J, bool _sym = 0)
    : structure( ), sym(_sym), n(I.max( )), m(J.max( )), raws(0), cols(0) {
    this->AddArrays(I, J);
  }

  ~SparseMatStructure( ) {
    if (raws) {
      delete raws;
    }

    if (cols) {
      delete cols;
    }
  }

  const_iterator begin( ) const { return structure.begin( ); }

  iterator begin( ) { return structure.begin( ); }

  const_iterator end( ) const { return structure.end( ); }

  iterator end( ) { return structure.end( ); }

  // Structure& operator()() {return structure;}
  // const Structure& operator()() const {return structure;}
  bool empty( ) const { return structure.empty( ) && !raws && !cols; }

  int N( ) const { return n; }

  int M( ) const { return m; }

  SparseMatStructure &clear( ) {
    structure.clear( );
    if (raws) {
      delete raws;
    }

    if (cols) {
      delete cols;
    }

    sym = false;
    n = 0;
    m = 0;
    return *this;
  }

  int size( ) const { return structure.size( ) ? structure.size( ) : (raws ? raws->N( ) : 0); }

  SparseMatStructure &AddMatrix(Matrice_Creuse< R > *);
  template< class INT >
  SparseMatStructure &AddArrays(const KN< INT > &, const KN< INT > &);
  SparseMatStructure &ToKn(bool emptystruct = false);

  KN< int > &Raws( ) { return *raws; }

  KN< int > const &Raws( ) const { return *raws; }

  KN< int > &Cols( ) { return *cols; }

  KN< int > const &Cols( ) const { return *cols; }

 private:
  int n, m;
  Structure structure;
  bool sym;
  KN< int > *raws, *cols;
};

SparseMatStructure &SparseMatStructure::ToKn(bool emptystruct) {
  if (raws) {
    delete raws;
  }

  if (cols) {
    delete cols;
  }

  raws = new KN< int >(structure.size( ));
  cols = new KN< int >(structure.size( ));
  int k = 0;

  for (const_iterator i = begin( ); i != end( ); ++i) {
    (*raws)[k] = i->first;
    (*cols)[k] = i->second;
    ++k;
  }

  if (emptystruct) {
    structure.clear( );
  }

  return *this;
}

SparseMatStructure &SparseMatStructure::AddMatrix(Matrice_Creuse< R > *const _M) {
  n = n > _M->N( ) ? n : _M->N( );
  m = m > _M->M( ) ? m : _M->M( );
  MatriceMorse< R > *M = _M->pHM( );
  if (!M) {
    cerr << " Err= "
         << " Matrix is not morse or CSR " << &(*_M->A) << endl;
    ffassert(M);
  }
  M->CSR( );
  {
    if (!sym || (sym && M->half)) {
      for (int i = 0; i < M->N; ++i) {
        for (int k = M->p[i]; k < M->p[i + 1]; ++k) {
          structure.insert(Z2(i, M->j[k]));
        }
      }
    } else {    // sym && !M->symetrique
      for (int i = 0; i < M->N; ++i) {
        for (int k = M->p[i]; k < M->p[i + 1]; ++k) {
          if (i >= M->j[k]) {
            structure.insert(Z2(i, M->j[k]));
          }
        }
      }
    }
  }
  return *this;
}

template< class INT >
SparseMatStructure &SparseMatStructure::AddArrays(const KN< INT > &I, const KN< INT > &J) {
  ffassert(I.N( ) == J.N( ));
  n = n > I.max( ) + 1 ? n : I.max( ) + 1;
  m = m > J.max( ) + 1 ? m : J.max( ) + 1;
  if (!sym) {
    for (int k = 0; k < I.N( ); ++k) {
      structure.insert(Z2(I[k], J[k]));
    }
  } else {
    for (int k = 0; k < I.N( ); ++k) {
      if (I[k] >= J[k]) {
        structure.insert(Z2(I[k], J[k]));
      }
    }
  }

  return *this;
}

/*****************************************************************************************************************************
 *	ffNLP : Derived from the TNLP non-linear problem wrapper class of Ipopt. Virtual methods are
 *defined as explain in the IPOPT documentation. Some of them are tricky because the sparse matrix
 *format in freefem is CRS, whereas IPOPT use COO storage. It is even more tricky because most of
 *time, FreeFem will remove null coefficient from the structure, leading to non constant indexing of
 *the coefficient through the algorithm in case of very non linear functions. As IPOPT need a
 *constant structure, a FindIndex method involving a dichotomic search has been implemented to
 *prevent the errors related to that.
 *****************************************************************************************************************************/
using namespace Ipopt;

class ffNLP : public TNLP {
 public:
  ffNLP( ) : xstart(0) {}

  ffNLP(Rn &, const Rn &, const Rn &, const Rn &, const Rn &, ScalarFunc *, VectorFunc *,
        SparseMatFunc *, VectorFunc *, SparseMatFunc *);
  ffNLP(Rn &, const Rn &, const Rn &, const Rn &, const Rn &, ScalarFunc *, VectorFunc *,
        SparseMatFunc *, VectorFunc *, SparseMatFunc *, int, int, int);
  virtual ~ffNLP( );

  bool get_nlp_info(Index &, Index &, Index &, Index &, IndexStyleEnum &);    // the IPOPT methods
  bool get_bounds_info(Index, Number *, Number *, Index, Number *, Number *);
  bool get_starting_point(Index, bool, Number *, bool, Number *, Number *, Index, bool, Number *);
  bool eval_f(Index, const Number *, bool, Number &);
  bool eval_grad_f(Index, const Number *, bool, Number *);
  bool eval_g(Index, const Number *, bool, Index, Number *);
  bool eval_jac_g(Index, const Number *, bool, Index, Index, Index *, Index *, Number *);
  bool eval_h(Index, const Number *, bool, Number, Index, const Number *, bool, Index, Index *,
              Index *, Number *);
  void finalize_solution(SolverReturn, Index, const Number *, const Number *, const Number *, Index,
                         const Number *, const Number *, Number, const IpoptData *ip_data,
                         IpoptCalculatedQuantities *ip_cq);

  template< class INT >
  ffNLP &SetHessianStructure(const KN< INT > &, const KN< INT > &, bool reset = 0);
  template< class INT >
  ffNLP &SetJacobianStructure(const KN< INT > &, const KN< INT > &, bool reset = 0);
  enum Level { do_nothing, user_defined, one_evaluation, basis_analysis };
  ffNLP &BuildMatrixStructures(Level, Level, int);
  ffNLP &EnableCheckStruct( ) {
    checkstruct = true;
    return *this;
  }

  ffNLP &DisableCheckStruct( ) {
    checkstruct = false;
    return *this;
  }

  Rn lambda_start, x_start, uz_start, lz_start;
  double sigma_start;
  double final_value;

 private:
  // algorithm datas
  Rn *xstart, xl, xu, gl, gu;
  ScalarFunc *fitness;    // Pointers to functions wrappers
  VectorFunc *dfitness, *constraints;
  SparseMatFunc *hessian, *dconstraints;
  int mm, nnz_jac, nnz_h;    // duplicated datas? did not seems to be reachable in the base class
  // bool sym;
  bool checkstruct;
  SparseMatStructure HesStruct, JacStruct;

  // some static functions...
  template< class A, class B >
  static void KnToPtr(const KN< A > &a, B *b) {
    for (int i = 0; i < a.N( ); ++i) {
      b[i] = a[i];
    }
  }    // Fill a pointer with a KN

  template< class A, class B >
  static void KnFromPtr(KN< A > &a, B const *b) {
    for (int i = 0; i < a.N( ); ++i) {
      a[i] = b[i];
    }
  }    // Fill a KN with a pointer <-- to avoid the use of const_cast

  static int FindIndex(const KN< int > &irow, const KN< int > &jrow, int i, int j, int kmin,
                       int kmax);
};

ffNLP::ffNLP(Rn &x, const Rn &_xl, const Rn &_xu, const Rn &_gl, const Rn &_gu,
             ScalarFunc *_fitness, VectorFunc *_dfitness, SparseMatFunc *_hessian,
             VectorFunc *_constraints, SparseMatFunc *_dconstraints)
  : xstart(&x), xl(_xl), xu(_xu), gl(_gl), gu(_gu),
    final_value(299792458.),    // sym(0),unsymind(),
    fitness(_fitness), dfitness(_dfitness), constraints(_constraints), uz_start( ), lz_start( ),
    hessian(_hessian), dconstraints(_dconstraints), mm(-1), nnz_jac(-1), nnz_h(-1), HesStruct(true),
    JacStruct(false), sigma_start(1.), lambda_start( ), x_start(x), checkstruct(1) {}

ffNLP::ffNLP(Rn &x, const Rn &_xl, const Rn &_xu, const Rn &_gl, const Rn &_gu,
             ScalarFunc *_fitness, VectorFunc *_dfitness, SparseMatFunc *_hessian,
             VectorFunc *_constraints, SparseMatFunc *_dconstraints, int _mm, int _nnz_jac,
             int _nnz_h)
  : xstart(&x), xl(_xl), xu(_xu), gl(_gl), gu(_gu), hessian(_hessian),
    final_value(299792458.),    // sym(0),unsymind(),
    fitness(_fitness), dfitness(_dfitness), constraints(_constraints), dconstraints(_dconstraints),
    uz_start( ), lz_start( ), mm(_mm), nnz_jac(_nnz_jac), nnz_h(_nnz_h), HesStruct(true),
    JacStruct(false), sigma_start(1.), lambda_start( ), x_start(x), checkstruct(1) {}

ffNLP::~ffNLP( ) {
  /*
   * clean(fitness);
   * clean(dfitness);
   * clean(constraints);
   * clean(hessian);
   * clean(dconstraints);
   */
}

template< class INT >
ffNLP &ffNLP::SetHessianStructure(const KN< INT > &I, const KN< INT > &J, bool reset) {
  if (reset) {
    HesStruct.clear( );
  }

  HesStruct.AddArrays(I, J);
  return *this;
}

template< class INT >
ffNLP &ffNLP::SetJacobianStructure(const KN< INT > &I, const KN< INT > &J, bool reset) {
  if (reset) {
    JacStruct.clear( );
  }

  JacStruct.AddArrays(I, J);
  return *this;
}

ffNLP &ffNLP::BuildMatrixStructures(Level hlvl, Level jlvl, int _mm) {
  if (jlvl != do_nothing && dconstraints) {
    if (jlvl == user_defined) {
      ffassert(JacStruct.size( ));
    } else if ((jlvl == one_evaluation || jlvl == basis_analysis) && dconstraints) {
      JacStruct.AddMatrix(dconstraints->J(x_start));
    }
  }

  if (hlvl != do_nothing && hessian) {
    if (hlvl == user_defined) {
      ffassert(HesStruct.size( ));
    } else if (hlvl == one_evaluation || !hessian->NLCHPEnabled( )) {
      Rn lms = lambda_start;
      lms = 1.;
      HesStruct.AddMatrix(hessian->J(x_start, sigma_start, lms));
    } else if (hlvl == basis_analysis) {
      {
        Rn lambda(_mm, 0.);
        HesStruct.AddMatrix(hessian->J(x_start, 1., lambda));
      }

      for (int i = 0; i < _mm; ++i) {
        Rn lambda(_mm, 0.);
        lambda[i] = 1.;
        HesStruct.AddMatrix(hessian->J(x_start, 0., lambda));
        lambda[i] = 0.;
      }
    }
  }

  JacStruct.ToKn( );
  HesStruct.ToKn( );
  return *this;
}

int ffNLP::FindIndex(const KN< int > &irow, const KN< int > &jcol, int i, int j, int kmin,
                     int kmax) {
  typedef std::pair< int, int > Z2;
  Z2 ij(i, j), ijmin(irow[kmin], jcol[kmin]), ijmax(irow[kmax], jcol[kmax]);
  if (abs(kmin - kmax) <= 1) {
    if (ij == ijmin) {
      return kmin;
    } else if (ij == ijmax) {
      return kmax;
    } else {
      return -1;
    }
  } else {
    int knew = (kmin + kmax) / 2;
    Z2 ijnew(irow[knew], jcol[knew]);
    if (ij <= ijnew) {
      return FindIndex(irow, jcol, i, j, kmin, knew);
    } else {
      return FindIndex(irow, jcol, i, j, knew, kmax);
    }
  }
}

bool ffNLP::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag,
                         IndexStyleEnum &index_style) {
  bool ret = true;

  n = xstart ? xstart->N( ) : (ret = 0);
  mm = m = constraints ? JacStruct.N( ) : 0;
  nnz_jac = nnz_jac_g = constraints ? JacStruct.size( ) : 0;
  nnz_h = nnz_h_lag = HesStruct.size( );
  index_style = TNLP::C_STYLE;
  return ret;
}

bool ffNLP::get_bounds_info(Index n, Number *x_l, Number *x_u, Index m, Number *g_l, Number *g_u) {
  KnToPtr(xl, x_l);
  KnToPtr(xu, x_u);
  if (mm) {
    KnToPtr(gl, g_l);
  }

  if (mm) {
    KnToPtr(gu, g_u);
  }

  /* DEBUG
   * cout << "constraints lower bound = (";
   * for(int i=0;i<m;++i) cout << g_l[i] <<  (i<m-1 ? ',':')');
   * cout << endl << "constraints upper bound = (";
   * for(int i=0;i<m;++i) cout << g_u[i] <<  (i<m-1 ? ',':')');
   * cout << endl;*/
  return true;
}

bool ffNLP::get_starting_point(Index n, bool init_x, Number *x, bool init_z, Number *z_L,
                               Number *z_U, Index m, bool init_lambda, Number *lambda) {
  assert(init_x == true);
  assert(xstart->N( ) == n);
  KnToPtr(*xstart, x);
  if (init_z) {
    if (uz_start.N( ) != n) {
      if (xu.min( ) < 1.e19) {
        cout << "ff-IPOPT warm start : upper simple bounds start multipliers array doesn't have "
                "the expected size ("
             << uz_start.N( ) << "!=" << n << ")." << endl;
        cout << "                   ";
        if (uz_start.N( ) == 0) {
          cout << "maybe because no upper bounds multiplier has been given. " << endl;
        }

        cout << " Initializing them to 1..." << endl;
      }

      uz_start.resize(n);
      uz_start = 1.;
    }

    if (lz_start.N( ) != n) {
      if (xl.max( ) > -1e19) {
        cout << "ff-IPOPT warm start : lower simple bounds start multipliers array doesn't have "
                "the expected size ("
             << lz_start.N( ) << "!=" << n << ")." << endl;
        cout << "                   ";
        if (lz_start.N( ) == 0) {
          cout << "maybe because no lower bounds multiplier has been given. " << endl;
        }

        cout << " Initializing them to 1..." << endl;
      }

      lz_start.resize(n);
      lz_start = 1.;
    }

    KnToPtr(uz_start, z_U);
    KnToPtr(lz_start, z_L);
  }

  if (init_lambda) {
    if (lambda_start.N( ) != m) {
      cout << "ff-IPOPT warm start : constraints start multipliers array doesn't have the expected "
              "size ("
           << lambda_start.N( ) << "!=" << m << ")." << endl;
      cout << "                   ";
      if (lambda_start.N( ) == 0) {
        cout << "maybe because no constraints multiplier has been given. " << endl;
      }

      cout << " Initializing them to 1..." << endl;
      lambda_start.resize(m);
      lambda_start = 1.;
    }

    KnToPtr(lambda_start, lambda);
  }

  return true;
}

bool ffNLP::eval_f(Index n, const Number *x, bool new_x, Number &obj_value) {
  assert(n == xstart->N( ));
  Rn X(n);
  KnFromPtr(X, x);
  obj_value = fitness->J(X);
  return true;
}

bool ffNLP::eval_grad_f(Index n, const Number *x, bool new_x, Number *grad_f) {
  assert(n == xstart->N( ));
  Rn X(n);
  KnFromPtr(X, x);
  Rn _grad_f = dfitness->J(X);
  KnToPtr(_grad_f, grad_f);
  return true;
}

bool ffNLP::eval_g(Index n, const Number *x, bool new_x, Index m, Number *g) {
  Rn X(n);

  KnFromPtr(X, x);
  if (constraints) {
    Rn _g = constraints->J(X);
    KnToPtr(_g, g);
  }

  return true;
}

bool ffNLP::eval_jac_g(Index n, const Number *x, bool new_x, Index m, Index nele_jac, Index *iRow,
                       Index *jCol, Number *values) {
  assert(n == xstart->N( ));
  Rn X(n);
  if (x) {
    KnFromPtr(X, x);
  } else {
    X = *xstart;
  }

  if (values == 0) {
    int k = 0;

    for (SparseMatStructure::const_iterator i = JacStruct.begin( ); i != JacStruct.end( ); ++i) {
      iRow[k] = i->first;
      jCol[k] = i->second;
      ++k;
    }
  } else if (dconstraints) {
    Matrice_Creuse< R > *M = dconstraints->J(X);
    MatriceMorse< R > *MM = dynamic_cast< MatriceMorse< R > * >(&(*M->A));    // ugly!
    MM->CSR( );
    for (int i = 0; i < MM->N; ++i) {
      for (int k = MM->p[i]; k < MM->p[i + 1]; ++k) {
        if (checkstruct) {
          int kipopt =
            FindIndex(JacStruct.Raws( ), JacStruct.Cols( ), i, MM->j[k], 0, nele_jac - 1);
          if (kipopt >= 0) {
            values[kipopt] = MM->aij[k];
          }
        } else {
          values[k] = MM->aij[k];
        }
      }
    }
  }

  return true;
}

bool ffNLP::eval_h(Index n, const Number *x, bool new_x, Number obj_factor, Index m,
                   const Number *lambda, bool new_lambda, Index nele_hess, Index *iRow, Index *jCol,
                   Number *values) {
  Rn X(n), L(m);

  if (x) {
    KnFromPtr(X, x);
  } else {
    X = *xstart;
  }

  if (lambda) {
    KnFromPtr(L, lambda);
  } else {
    L = 0.;
  }

  bool NLCHPE = hessian->NLCHPEnabled( );
  Number _obj_factor = NLCHPE ? 1. : obj_factor;
  if (values == 0) {
    int k = 0;

    for (SparseMatStructure::const_iterator i = HesStruct.begin( ); i != HesStruct.end( ); ++i) {
      iRow[k] = i->first;
      jCol[k] = i->second;
      ++k;
    }
  } else {
    Matrice_Creuse< R > *M = 0;
    if (NLCHPE) {
      M = hessian->J(X, obj_factor, L);
    } else {
      M = hessian->J(X);
    }

    MatriceMorse< R > *MM = dynamic_cast< MatriceMorse< R > * >(&(*M->A));    // ugly!
    MM->CSR( );
    if (MM) {
      if (checkstruct) {
        for (int i = 0; i < MM->N; ++i) {
          for (int k = MM->p[i]; k < MM->p[i + 1]; ++k) {
            int kipopt =
              FindIndex(HesStruct.Raws( ), HesStruct.Cols( ), i, MM->j[k], 0, nele_hess - 1);
            if (kipopt >= 0) {
              values[kipopt] = _obj_factor * (MM->aij[k]);
            }
          }
        }
      } else if (!MM->half) {
        for (int i = 0, kipopt = 0; i < MM->N; ++i) {
          for (int k = MM->p[i]; k < MM->p[i + 1]; ++k) {
            if (i >= MM->j[k]) {
              values[kipopt] = _obj_factor * (MM->aij[k]);
              ++kipopt;
            }
          }
        }
      } else {
        for (int i = 0; i < MM->N; ++i) {
          for (int k = MM->p[i]; k < MM->p[i + 1]; ++k) {
            values[k] = _obj_factor * (MM->aij[k]);
          }
        }
      }
    }
  }

  return true;
}

void ffNLP::finalize_solution(SolverReturn status, Index n, const Number *x, const Number *z_L,
                              const Number *z_U, Index m, const Number *g, const Number *lambda,
                              Number obj_value, const IpoptData *ip_data,
                              IpoptCalculatedQuantities *ip_cq) {
  KnFromPtr(*xstart, x);
  KnFromPtr(lambda_start, lambda);
  KnFromPtr(lz_start, z_L);
  KnFromPtr(uz_start, z_U);
  final_value = obj_value;
}

/*****************************************************************************************************************************
 *	Assumptions : these are tags used as template parameters for case specific function wrapping
 *or warning message in the interface. Some case can be added here (but some class has to be
 *specialized for the new cases) AssumptionF : undeff              --> undefined case (not used)
 *                no_assumption_f     --> most general case when the fitness function and all its
 *derivative are coded in the freefem script with the func keyword (type Polymorphic in c++). These
 *functions are then wrapped in GeneralFunc objects. P2_f                --> no longer used (it was
 *used for fitness and its gradient defined as func in freefem script, while the hessian is a
 *constant matrix directly given to the interface, but it leads to ambiguities). unavailable_hessian
 *--> fitness function and its gradients coded with func in the ff script, wrapped into GeneralFunc
 *objects, without second order derivative function. Enables the BFGS option of IPOPT. mv_P2_f -->
 *fitness function is a P2 function which will be defined by a [matrix,vector] array. The functions
 *are passed to the ffNLP object as P2ScalarFunc, P1VectorFunc and ConstantSparseMatFunc
 *respectively for the fitness function, its gradient and its hessian (with all vf=1). quadratic_f
 *--> f is a pure quadratic fonction, defined by a single matrix. Same type as mv_P2_f for function
 *wrappers with a vf=0 tag. linear_f            --> f is a linear form, defined by a single vector.
 *Same type as mv_P2_f for function wrappers with a vf=0 tag. AssumptionG : undeff              -->
 *undefined case (not used) no_assumption_f     --> most general case when the constraint functions
 *and all its derivative are coded in the freefem script with the func keyword (type Polymorphic in
 *c++). These functions are then wrapped in GeneralFunc objects. P1_g                --> no longer
 *used (it was used for constraints defined as func in freefem script , while the jacobian is a
 *constant matrix directly given to the interface, but it leads to ambiguities). mv_P1_g -->
 *Constraints function is a P1 function which will be defined by a [matrix,vector] array. The
 *functions are passed to the ffNLP object as P1VectorFunc and ConstantSparseMatFunc respectively
 *for the constraints and its jacobian (with all vf=0). linear_g            --> Constraints are
 *linear, defined by a single matrix. Same type as mv_P1_g for function wrappers with a vf=0 tag.
 *  Case : templatized with a pair of AssumptionF and AssumptionG, is used to build different
 *constructor for the interface class in order to overload the freefem function which will call
 *IPOPT
 *****************************************************************************************************************************/

enum AssumptionF {
  undeff,
  no_assumption_f,
  P2_f,
  unavailable_hessian,
  mv_P2_f,
  quadratic_f,
  linear_f
};
enum AssumptionG { undefg, without_constraints, no_assumption_g, P1_g, mv_P1_g, linear_g };

template< AssumptionF AF, AssumptionG AG >
struct Case {
  Case( ) {}

  static const AssumptionF af = AF;
  static const AssumptionG ag = AG;
};

/*****************************************************************************************************************************
 *	CheckMatrixVectorPair : Small function taking an E_Array and check whether the type of the 2
 *objects contained in the array are matrix and vector. Returns false if types are not
 *matrix/vector. order is modified to know whether the matrix is in first position or not.
 *****************************************************************************************************************************/
bool CheckMatrixVectorPair(const E_Array *mv, bool &order) {
  const aType t1 = (*mv)[0].left( ), t2 = (*mv)[1].left( );

  if (NXOR(t1 == atype< Matrice_Creuse< R > * >( ), t2 == atype< Matrice_Creuse< R > * >( ))) {
    return false;
  } else if (NXOR(t1 == atype< Rn * >( ), t2 == atype< Rn * >( ))) {
    return false;
  } else {
    order = (t1 == atype< Matrice_Creuse< R > * >( ));
    return true;
  }
}

/*****************************************************************************************************************************
 *	The following class offers a polymorphic way to build the function wrappers to pass to the
 *ffNLP object Each element of the assumption enum define a "FunctionDatas" class in which the
 *constructor and the operator() makes case specific task. If some new value in the Assumption enums
 *are to be added, the FitnessFunctionDatas and/or ConstraintFunctionDatas with the new value as
 *template parameter has to be specialized. What should the method do is (exemple at the end of the
 *file with already coded cases): Constructor : define the Expression members using the arguments
 *passed to the IPOPT function in the script operator()  : allocate with appropriate dynamic type
 *the ScalarFunc, VectorFunc, SparseMatFunc pointers, and display some case dependant errors or
 *warnings (note that there is no ScalarFunc ptr to allocate for constraints)
 *****************************************************************************************************************************/
class GenericFitnessFunctionDatas {
 public:
  static GenericFitnessFunctionDatas *New(AssumptionF, const basicAC_F0 &, Expression const *,
                                          const C_F0 &, const C_F0 &, const C_F0 &);
  bool CompletelyNonLinearConstraints;
  Expression JJ, GradJ, Hessian;
  GenericFitnessFunctionDatas( )
    : CompletelyNonLinearConstraints(true), JJ(0), GradJ(0), Hessian(0) {}

  virtual const AssumptionF A( ) const { return undeff; }

  virtual void operator( )(Stack, const C_F0 &, const C_F0 &, const C_F0 &, Expression const *,
                           ScalarFunc *&, VectorFunc *&, SparseMatFunc *&,
                           bool) const = 0;    // Build the functions
  virtual ~GenericFitnessFunctionDatas( ) {}
};

template< AssumptionF AF >
class FitnessFunctionDatas
  : public GenericFitnessFunctionDatas    // not really a template, since most of the methods of all
                                          // cases have to be specialized
{
 public:
  FitnessFunctionDatas(const basicAC_F0 &, Expression const *, const C_F0 &, const C_F0 &,
                       const C_F0 &);
  const AssumptionF A( ) const { return AF; }

  void operator( )(Stack, const C_F0 &, const C_F0 &, const C_F0 &, Expression const *,
                   ScalarFunc *&, VectorFunc *&, SparseMatFunc *&, bool) const;
};

class GenericConstraintFunctionDatas {
 public:
  static GenericConstraintFunctionDatas *New(AssumptionG, const basicAC_F0 &, Expression const *,
                                             const C_F0 &);
  Expression Constraints, GradConstraints;
  GenericConstraintFunctionDatas( ) : Constraints(0), GradConstraints(0) {}

  virtual const AssumptionG A( ) const { return undefg; }

  virtual const bool WC( ) const = 0;    // with constraints
  virtual void operator( )(Stack, const C_F0 &, Expression const *, VectorFunc *&, SparseMatFunc *&,
                           bool) const = 0;    // build the functions`
  virtual ~GenericConstraintFunctionDatas( ) {}
};

template< AssumptionG AG >
class ConstraintFunctionDatas : public GenericConstraintFunctionDatas {
 public:
  ConstraintFunctionDatas(const basicAC_F0 &, Expression const *, const C_F0 &);
  const AssumptionG A( ) const { return AG; }

  const bool WC( ) const { return AG != without_constraints; }

  void operator( )(Stack, const C_F0 &, Expression const *, VectorFunc *&, SparseMatFunc *&,
                   bool) const;
};

/*****************************************************************************************************************************
 *	OptimIpopt & OptimIpopt::E_Ipopt - The interface class
 *  Do the link beetween freefem and Ipopt
 *****************************************************************************************************************************/
class OptimIpopt : public OneOperator {
 public:
  const AssumptionF AF;
  const AssumptionG AG;
  class E_Ipopt : public E_F0mps {
   private:
    bool spurious_cases;

   public:
    const AssumptionF AF;
    const AssumptionG AG;
    const bool WC;
    std::set< unsigned short > unused_name_param;    // In some case, some parameter are usless,
                                                     // this is the list of their index in nargs
    void InitUNP( );    // Initialize unusued_name_param at freefem compile time
    static basicAC_F0::name_and_type name_param[];
    static const int n_name_param = 29;
    Expression nargs[n_name_param];
    Expression X;
    mutable Rn lm;
    C_F0 L_m;
    C_F0 inittheparam, theparam, closetheparam;
    C_F0 initobjfact, objfact;
    GenericFitnessFunctionDatas *fitness_datas;
    GenericConstraintFunctionDatas *constraints_datas;
    bool arg(int i, Stack stack, bool a) const {
      return nargs[i] ? GetAny< bool >((*nargs[i])(stack)) : a;
    }

    long arg(int i, Stack stack, long a) const {
      return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
    }

    R arg(int i, Stack stack, R a) const { return nargs[i] ? GetAny< R >((*nargs[i])(stack)) : a; }

    Rn_ arg(int i, Stack stack, Rn_ a) const {
      return nargs[i] ? GetAny< Rn_ >((*nargs[i])(stack)) : a;
    }

    template< typename T >
    T Arg(int i, Stack s) const {
      return GetAny< T >((*nargs[i])(s));
    }

    E_Ipopt(const basicAC_F0 &args, AssumptionF af, AssumptionG ag)
      : lm( ), L_m(CPValue(lm)), AF(af), AG(ag), WC(ag != without_constraints),
        unused_name_param( ), spurious_cases(false), fitness_datas(0), constraints_datas(0) {
      InitUNP( );
      int nbj = args.size( ) - 1;
      Block::open(currentblock);    // make a new block to
      X = to< Rn * >(args[nbj]);
      C_F0 X_n(args[nbj], "n");
      // the expression to init the theparam of all
      inittheparam =
        currentblock->NewVar< LocalVariable >("the parameter", atype< KN< R > * >( ), X_n);
      initobjfact = currentblock->NewVar< LocalVariable >("objective factor", atype< double * >( ));
      theparam = currentblock->Find("the parameter");    // the expression for the parameter
      objfact = currentblock->Find("objective factor");
      args.SetNameParam(n_name_param, name_param, nargs);
      fitness_datas = GenericFitnessFunctionDatas::New(
        AF, args, nargs, theparam, objfact, L_m);    // Creates links to the freefem objects
      constraints_datas =
        GenericConstraintFunctionDatas::New(AG, args, nargs, theparam);    // defining the functions
      spurious_cases = AG == no_assumption_g &&
                       (AF == P2_f || AF == mv_P2_f || AF == quadratic_f || AF == linear_f);
      closetheparam = C_F0((Expression)Block::snewclose(currentblock), atype< void >( ));
    }

    ~E_Ipopt( ) {
      if (fitness_datas) {
        delete fitness_datas;
      }

      if (constraints_datas) {
        delete constraints_datas;
      }
    }

    virtual AnyType operator( )(Stack stack) const {
      double cost = nan("");

      WhereStackOfPtr2Free(stack) = new StackOfPtr2Free(stack);    // FH mars 2005
      Rn &x = *GetAny< Rn * >((*X)(stack));
      {
        Expression test(theparam);    // in some case the KN object associated to the param is never
                                      // initialized, leading to failed assertion in KN::destroy
        Rn *tt = GetAny< Rn * >((*test)(stack));    // this lines prevent this to happen
        *tt = x;
      }
      long n = x.N( );
      bool warned = false;
      cout << endl;
      if (spurious_cases) {
        cout << "ff-IPOPT Spurious case detected : the hessian is defined as a constant matrix but "
                "constraints are given in function form."
             << endl;
        cout << "If they are not affine, the optimization is likely to fail. In this case, try one "
                "of the following suggestions:"
             << endl;
        cout << "  - if constraints have computable hessians, use function form for the fitness "
                "function and all its derivatives"
             << endl;
        cout
          << "    and check the documentation to know how to express the whole lagrangian hessian."
          << endl;
        cout << "  - if constraints hessians are difficult to obtain, force the BFGS mode using "
                "named parameter "
             << name_param[12].name << '.' << endl;
        cout << "Do not worry about this message if you know all your constraints has a constant "
                "null hessian."
             << endl
             << endl;
      }

      if (nargs[7]) {
        cout << "ff-IPOPT : the named parameter autostruct is no longer used in this version of "
                "the interface."
             << endl;
      }

      // Detection of mixed case dependant warnings or error
      for (int i = 0; i < n_name_param; ++i) {
        if (nargs[i] && unused_name_param.find(i) != unused_name_param.end( )) {
          cout << "ff-IPOPT Warning: named parameter " << name_param[i].name
               << " is useless for the problem you have set." << endl;
          warned = true;
        }
      }

      if (nargs[4] && nargs[5] && nargs[7]) {
        cout << "ff-IPOPT Warning: both " << name_param[4].name << " and " << name_param[5].name
             << " has been defined, so " << name_param[7].name;
        cout << " will be ignored." << endl;
      }

      if (warned) {
        if (!WC && AF == unavailable_hessian && nargs[8]) {
          cout << "  ==> " << name_param[8].name
               << " is useless because there should not be any function returning matrix in your "
                  "problem,"
               << endl;
          cout << "      (2 functions can only be J and dJ). You may as well have forgotten one "
                  "function (IPOPT will certainly crash if so)."
               << endl;
        }

        if (AF != no_assumption_f && AF != unavailable_hessian && AG != no_assumption_g &&
            nargs[5]) {
          cout << "  ==> your lagrangian hessian is a constant matrix, so there is no need to "
                  "specify its structure with "
               << name_param[5].name << endl;
          cout << "      since it is contained in the matrix object." << endl;
        }

        if (AF != no_assumption_f && AF != unavailable_hessian && AG != no_assumption_g &&
            nargs[7]) {
          cout << "  ==> " << name_param[7].name
               << " will be ignored since all matrices are constants and constraints do not"
               << endl;
          cout << "      contribute to the hessian, matrix structure determination is trivial."
               << endl;
        }

        if (AF == unavailable_hessian && AG != no_assumption_g && (nargs[7] || nargs[8])) {
          cout << "  ==> " << name_param[7].name << " or " << name_param[8].name
               << " will be ignored since the only matrix you have passed is constant. " << endl;
          cout << "      Or maybe did you forget to pass a function (IPOPT will certainly crash if "
                  "so)."
               << endl;
        }

        if (AF != no_assumption_f && AF != unavailable_hessian && AG != no_assumption_g &&
            nargs[8]) {
          cout << "  ==> no need to use " << name_param[8].name
               << " since all matrices are constant, structures won't change through the algorithm,"
               << endl;
          cout << "      it is automatically set to the default disabling value." << endl;
        }
      }

      long iprint = verbosity;
      ScalarFunc *ffJ = 0;
      VectorFunc *ffdJ = 0, *ffC = 0;
      SparseMatFunc *ffH = 0, *ffdC = 0;

      (*fitness_datas)(stack, theparam, objfact, L_m, nargs, ffJ, ffdJ, ffH,
                       warned);    // Fill the functions
      (*constraints_datas)(stack, theparam, nargs, ffC, ffdC, warned);

      Rn xl(n), xu(n), gl(nargs[2] ? Arg< Rn_ >(2, stack).N( ) : 0),
        gu(nargs[3] ? Arg< Rn_ >(3, stack).N( ) : 0);
      int mmm = 0;
      if (WC && (gl.N( ) + gu.N( )) == 0) {
        cout << "IPOPT Warning : constrained problem without constraints bounds." << endl;
        mmm = ffC->J(x).N( );
      } else {
        mmm = gl.N( ) > gu.N( ) ? gl.N( ) : gu.N( );
      }

      Rn_ *lag_mul = 0, *l_z = 0, *u_z = 0;    // Rn(mmm,1.);
      // int niter=arg(6,stack,100L);
      int autostructmode = ffNLP::one_evaluation;
      bool checkindex =
             (AF != no_assumption_f && AG != no_assumption_g) ? false : arg(8, stack, true),
           cberror = false;

      if (nargs[0]) {
        xl = Arg< Rn_ >(0, stack);
      } else {
        xl = -1.e19;
      }

      if (nargs[1]) {
        xu = Arg< Rn_ >(1, stack);
      } else {
        xu = 1.e19;
      }

      if (nargs[2]) {
        gl = Arg< Rn_ >(2, stack);
      } else {
        gl.resize(mmm);
        gl = -1.e19;
      }

      if (nargs[3]) {
        gu = Arg< Rn_ >(3, stack);
      } else {
        gu.resize(mmm);
        gu = 1.e19;
      }

      const E_Array *ejacstruct = (WC && AF == no_assumption_f && AG == no_assumption_g && nargs[4])
                                    ? dynamic_cast< const E_Array * >(nargs[4])
                                    : 0,
                    *ehesstruct = (AF == no_assumption_f && nargs[5])
                                    ? dynamic_cast< const E_Array * >(nargs[5])
                                    : 0;

      if (nargs[6] && WC) {
        lag_mul = new Rn_(GetAny< Rn_ >((*nargs[6])(stack)));
      }

      if (nargs[21]) {
        l_z = new Rn_(GetAny< Rn_ >((*nargs[21])(stack)));
      }

      if (nargs[20]) {
        u_z = new Rn_(GetAny< Rn_ >((*nargs[20])(stack)));
      }

      SmartPtr< TNLP > optim = new ffNLP(x, xl, xu, gl, gu, ffJ, ffdJ, ffH, ffC, ffdC);
      ffNLP *_optim = dynamic_cast< ffNLP * >(&(*optim));
      assert(_optim);
      if (WC && nargs[6]) {
        _optim->lambda_start = *lag_mul;
      } else if (WC) {
        _optim->lambda_start.resize(mmm);
        _optim->lambda_start = 1.;
      }

      _optim->sigma_start = 1.;
      if (nargs[21] && nargs[0]) {
        _optim->lz_start = *l_z;
      } else if (nargs[0]) {
        _optim->lz_start.resize(xl.N( ));
        _optim->lz_start = 1.;
      }

      if (nargs[20] && nargs[1]) {
        _optim->uz_start = *u_z;
      } else if (nargs[1]) {
        _optim->uz_start.resize(xu.N( ));
        _optim->uz_start = 1.;
      }

      if (ejacstruct) {
        if (ejacstruct->nbitem( ) != 2) {
          ExecError(
            "\nSorry, we were expecting an array with two componants in structjac=[iraw,jcol]");
        }

        if ((*ejacstruct)[0].left( ) != atype< KN< long > * >( )) {
          CompileError("Sorry, array componants in structjac=[iraw,jcol] must be integer arrays");
        }

        if ((*ejacstruct)[1].left( ) != atype< KN< long > * >( )) {
          CompileError("Sorry, array componants in structjac=[iraw,jcol] must be integer arrays");
        }

        Expression raws = (*ejacstruct)[0], cols = (*ejacstruct)[1];
        _optim->SetJacobianStructure(*GetAny< KN< long > * >((*raws)(stack)),
                                     *GetAny< KN< long > * >((*cols)(stack)), true);
      }

      if (ehesstruct) {
        if (ehesstruct->nbitem( ) != 2) {
          ExecError(
            "\nSorry, we were expecting an array with two componants in structhess=[iraw,jcol]");
        }

        if ((*ehesstruct)[0].left( ) != atype< KN< long > * >( )) {
          CompileError("Sorry, array componants in structhess=[iraw,jcol] must be integer arrays");
        }

        if ((*ehesstruct)[1].left( ) != atype< KN< long > * >( )) {
          CompileError("Sorry, array componants in structhess=[iraw,jcol] must be integer arrays");
        }

        Expression raws = (*ehesstruct)[0], cols = (*ehesstruct)[1];
        _optim->SetHessianStructure(*GetAny< KN< long > * >((*raws)(stack)),
                                    *GetAny< KN< long > * >((*cols)(stack)), true);
      }

      ffNLP::Level lh = ehesstruct ? ffNLP::user_defined : ffNLP::Level(autostructmode),
                   lj = ejacstruct ? ffNLP::user_defined : ffNLP::Level(autostructmode);
      if (AF == unavailable_hessian) {
        lh = ffNLP::do_nothing;
      }

      _optim->BuildMatrixStructures(lh, lj, mmm);
      if (checkindex) {
        _optim->EnableCheckStruct( );
      }

      SmartPtr< IpoptApplication > app = new IpoptApplication( );

      if (nargs[9]) {
        app->Options( )->SetNumericValue("tol", GetAny< double >((*nargs[9])(stack)));
      }

      if (nargs[10]) {
        app->Options( )->SetIntegerValue("max_iter", GetAny< long >((*nargs[10])(stack)));
      }

      if (nargs[11]) {
        app->Options( )->SetNumericValue("max_cpu_time", GetAny< double >((*nargs[11])(stack)));
      }

      bool bfgs = nargs[12] ? GetAny< bool >((*nargs[12])(stack)) : false;
      if (AF == unavailable_hessian || bfgs) {
        if (AF == unavailable_hessian && !bfgs) {
          cout << "IPOPT Note : No hessian given ==> LBFGS hessian approximation enabled" << endl;
        }

        app->Options( )->SetStringValue("hessian_approximation", "limited-memory");
      }

      if (nargs[13]) {
        string derivative_test = *GetAny< string * >((*nargs[13])(stack));
        app->Options( )->SetStringValue("derivative_test", derivative_test.c_str( ));
      }

      if (nargs[14]) {
        string options_file = *GetAny< string * >((*nargs[14])(stack));
        app->Options( )->SetStringValue("option_file_name", options_file.c_str( ));
      }

      if (nargs[15]) {
        app->Options( )->SetIntegerValue("print_level", GetAny< long >((*nargs[15])(stack)));
      }

      if (AG == without_constraints || AG == mv_P1_g || AG == linear_g) {
        app->Options( )->SetStringValue("jac_c_constant", "yes");
        app->Options( )->SetStringValue("jac_d_constant", "yes");
      }

      if (AF == mv_P2_f || AF == quadratic_f || AF == linear_f) {
        app->Options( )->SetStringValue("hessian_constant", "yes");
      }

      if (nargs[16]) {
        app->Options( )->SetNumericValue("derivative_test_perturbation",
                                         GetAny< double >((*nargs[16])(stack)));
      }

      if (nargs[17]) {
        app->Options( )->SetNumericValue("derivative_test_tol",
                                         GetAny< double >((*nargs[16])(stack)));
      }

      if (nargs[18]) {
        app->Options( )->SetStringValue("fixed_variable_treatment",
                                        GetAny< string * >((*nargs[18])(stack))->c_str( ));
      }

      if (nargs[19]) {
        app->Options( )->SetStringValue("warm_start_init_point", "yes");
        if (WC && !nargs[6]) {
          cout << "ff-IPOPT Warning : warm start for constrained problem without initial "
                  "constraints dual variables ("
               << name_param[6].name << " parameter)." << endl;
          cout << "                   ==> Starting with " << name_param[6].name << "=(1,1,...,1)."
               << endl;
        }

        if (nargs[0] && !nargs[21]) {
          cout << "ff-IPOPT Warning : warm start with simple lower bounds without initial lower "
                  "bounds dual variables ("
               << name_param[21].name << " parameter)." << endl;
          cout << "                   ==> Starting with " << name_param[21].name << "=(1,1,...,1)."
               << endl;
        }

        if (nargs[1] && !nargs[20]) {
          cout << "ff-IPOPT Warning : warm start with simple upper bounds without initial upper "
                  "bounds dual variables ("
               << name_param[20].name << " parameter)." << endl;
          cout << "                   ==> Starting with " << name_param[20].name << "=(1,1,...,1)."
               << endl;
        }

        if (l_z) {
          _optim->lz_start = *l_z;
        }

        if (u_z) {
          _optim->uz_start = *u_z;
        }

        if (lag_mul) {
          _optim->lambda_start = *lag_mul;
        }
      }

      if (nargs[22]) {
        app->Options( )->SetNumericValue("mu_init", GetAny< double >((*nargs[22])(stack)));
      } else {
        app->Options( )->SetStringValue("mu_strategy", "adaptive");
      }

      if (nargs[23]) {
        app->Options( )->SetNumericValue("mumps_pivtol", GetAny< double >((*nargs[23])(stack)));
      }

      if (nargs[24]) {
        app->Options( )->SetNumericValue("bound_relax_factor",
                                         GetAny< double >((*nargs[24])(stack)));
      }

      if (nargs[25]) {
        app->Options( )->SetStringValue("mu_strategy",
                                        GetAny< string * >((*nargs[25])(stack))->c_str( ));
      }

      if (nargs[27]) {
        app->Options( )->SetNumericValue("mu_min", GetAny< double >((*nargs[27])(stack)));
      }

      if (nargs[28]) {
        if (!GetAny< bool >((*nargs[28])(stack))) {
          app->Options( )->SetStringValue("accept_every_trial_step", "yes");
        }
      }

      if (verbosity > 1) {
        app->Options( )->SetStringValue("print_user_options", "yes");
      }

      app->Options( )->SetStringValue("output_file", "ipopt.out");
      if (AF != no_assumption_f && AF != unavailable_hessian && AG != no_assumption_g) {
        app->Options( )->SetStringValue("mehrotra_algorithm", "yes");
      }

      ApplicationReturnStatus status;
      app->Initialize( );

      // Ask Ipopt to solve the problem
      status = app->OptimizeTNLP(optim);

      if (lag_mul) {
        *lag_mul = _optim->lambda_start;
      }

      if (l_z) {
        *l_z = _optim->lz_start;
      }

      if (u_z) {
        *u_z = _optim->uz_start;
      }

      cost = _optim->final_value;

      if (nargs[26]) {
        double *pfv = GetAny< double * >((*nargs[26])(stack));
        *pfv = cost;
      }

      if (verbosity) {
        if (status == Solve_Succeeded) {
          printf("\n\n*** Ipopt succeeded \n");
        } else if (static_cast< int >(status) < 0) {
          printf("\n\n*** Ipopt failure!\n");
        } else {
          printf("\n\n*** Ipopt mixed results.\n");
        }
      }

      clean(lag_mul);
      clean(l_z);
      clean(u_z);
      clean(ffJ);
      clean(ffdJ);
      clean(ffH);
      clean(ffC);
      clean(ffdC);
      if (lm) {
        lm.destroy( );    // clean memory of LM
      }

      closetheparam.eval(stack);                // clean memory
      WhereStackOfPtr2Free(stack)->clean( );    // FH mars 2005
      return SetAny< long >(static_cast< long >(
        static_cast< int >(status)));    // SetAny<long>(0);  Modif FH  july 2005
    }

    operator aType( ) const { return atype< long >( ); }
  };

  E_F0 *code(const basicAC_F0 &args) const { return new E_Ipopt(args, AF, AG); }

  // Constructors - they define the different prototype of the overloaded IPOPT function reachable
  // in freefem scripts

  OptimIpopt(Case< no_assumption_f, no_assumption_g >)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< Polymorphic * >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< KN< R > * >( )),
      AF(no_assumption_f), AG(no_assumption_g) {}

  OptimIpopt(Case< no_assumption_f, without_constraints >)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< Polymorphic * >( ), atype< KN< R > * >( )),
      AF(no_assumption_f), AG(without_constraints) {}

  OptimIpopt(Case< no_assumption_f, P1_g >)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< Matrice_Creuse< R > * >( ), atype< KN< R > * >( )),
      AF(no_assumption_f), AG(P1_g) {}

  OptimIpopt(Case< no_assumption_f, mv_P1_g >)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< Polymorphic * >( ), atype< E_Array >( ), atype< KN< R > * >( )),
      AF(no_assumption_f), AG(mv_P1_g) {}

  OptimIpopt(Case< no_assumption_f, linear_g >)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< Polymorphic * >( ), atype< Matrice_Creuse< R > * >( ),
                  atype< KN< R > * >( )),
      AF(no_assumption_f), AG(linear_g) {}

  OptimIpopt(Case< P2_f, P1_g >)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< Matrice_Creuse< R > * >( ), atype< Polymorphic * >( ),
                  atype< Matrice_Creuse< R > * >( ), atype< KN< R > * >( )),
      AF(P2_f), AG(P1_g) {}

  OptimIpopt(Case< P2_f, without_constraints >)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< Matrice_Creuse< R > * >( ), atype< KN< R > * >( )),
      AF(P2_f), AG(without_constraints) {}

  OptimIpopt(Case< P2_f, no_assumption_g >)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< Matrice_Creuse< R > * >( ), atype< Polymorphic * >( ),
                  atype< Polymorphic * >( ), atype< KN< R > * >( )),
      AF(P2_f), AG(no_assumption_g) {}

  OptimIpopt(Case< P2_f, mv_P1_g >)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< Matrice_Creuse< R > * >( ), atype< E_Array >( ), atype< KN< R > * >( )),
      AF(P2_f), AG(mv_P1_g) {}

  OptimIpopt(Case< P2_f, linear_g >)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< Matrice_Creuse< R > * >( ), atype< Matrice_Creuse< R > * >( ),
                  atype< KN< R > * >( )),
      AF(P2_f), AG(linear_g) {}

  OptimIpopt(Case< unavailable_hessian, no_assumption_g >)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< Polymorphic * >( ), atype< Polymorphic * >( ), atype< KN< R > * >( )),
      AF(unavailable_hessian), AG(no_assumption_g) {}

  OptimIpopt(Case< unavailable_hessian, without_constraints >)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< KN< R > * >( )),
      AF(unavailable_hessian), AG(without_constraints) {}

  OptimIpopt(Case< unavailable_hessian, P1_g >)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< Polymorphic * >( ), atype< Matrice_Creuse< R > * >( ),
                  atype< KN< R > * >( )),
      AF(unavailable_hessian), AG(P1_g) {}

  OptimIpopt(Case< unavailable_hessian, mv_P1_g >)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< E_Array >( ), atype< KN< R > * >( )),
      AF(unavailable_hessian), AG(mv_P1_g) {}

  OptimIpopt(Case< unavailable_hessian, linear_g >)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< Matrice_Creuse< R > * >( ), atype< KN< R > * >( )),
      AF(unavailable_hessian), AG(linear_g) {}

  OptimIpopt(Case< mv_P2_f, no_assumption_g >)
    : OneOperator(atype< long >( ), atype< E_Array >( ), atype< Polymorphic * >( ),
                  atype< Polymorphic * >( ), atype< KN< R > * >( )),
      AF(mv_P2_f), AG(no_assumption_g) {}

  OptimIpopt(Case< mv_P2_f, without_constraints >)
    : OneOperator(atype< long >( ), atype< E_Array >( ), atype< KN< R > * >( )), AF(mv_P2_f),
      AG(without_constraints) {}

  OptimIpopt(Case< mv_P2_f, P1_g >)
    : OneOperator(atype< long >( ), atype< E_Array >( ), atype< Polymorphic * >( ),
                  atype< Matrice_Creuse< R > * >( ), atype< KN< R > * >( )),
      AF(mv_P2_f), AG(P1_g) {}

  OptimIpopt(Case< mv_P2_f, mv_P1_g >)
    : OneOperator(atype< long >( ), atype< E_Array >( ), atype< E_Array >( ),
                  atype< KN< R > * >( )),
      AF(mv_P2_f), AG(mv_P1_g) {}

  OptimIpopt(Case< mv_P2_f, linear_g >)
    : OneOperator(atype< long >( ), atype< E_Array >( ), atype< Matrice_Creuse< R > * >( ),
                  atype< KN< R > * >( )),
      AF(mv_P2_f), AG(linear_g) {}

  OptimIpopt(Case< quadratic_f, no_assumption_g >)
    : OneOperator(atype< long >( ), atype< Matrice_Creuse< R > * >( ), atype< Polymorphic * >( ),
                  atype< Polymorphic * >( ), atype< KN< R > * >( )),
      AF(quadratic_f), AG(no_assumption_g) {}

  OptimIpopt(Case< quadratic_f, without_constraints >)
    : OneOperator(atype< long >( ), atype< Matrice_Creuse< R > * >( ), atype< KN< R > * >( )),
      AF(quadratic_f), AG(without_constraints) {}

  OptimIpopt(Case< quadratic_f, P1_g >)
    : OneOperator(atype< long >( ), atype< Matrice_Creuse< R > * >( ), atype< Polymorphic * >( ),
                  atype< Matrice_Creuse< R > * >( ), atype< KN< R > * >( )),
      AF(quadratic_f), AG(P1_g) {}

  OptimIpopt(Case< quadratic_f, mv_P1_g >)
    : OneOperator(atype< long >( ), atype< Matrice_Creuse< R > * >( ), atype< E_Array >( ),
                  atype< KN< R > * >( )),
      AF(quadratic_f), AG(mv_P1_g) {}

  OptimIpopt(Case< quadratic_f, linear_g >)
    : OneOperator(atype< long >( ), atype< Matrice_Creuse< R > * >( ),
                  atype< Matrice_Creuse< R > * >( ), atype< KN< R > * >( )),
      AF(quadratic_f), AG(linear_g) {}

  OptimIpopt(Case< linear_f, no_assumption_g >)
    : OneOperator(atype< long >( ), atype< KN< R > * >( ), atype< Polymorphic * >( ),
                  atype< Polymorphic * >( ), atype< KN< R > * >( )),
      AF(linear_f), AG(no_assumption_g) {}

  OptimIpopt(Case< linear_f, without_constraints >)
    : OneOperator(atype< long >( ), atype< KN< R > * >( ), atype< KN< R > * >( )), AF(linear_f),
      AG(without_constraints) {}

  OptimIpopt(Case< linear_f, P1_g >)
    : OneOperator(atype< long >( ), atype< KN< R > * >( ), atype< Polymorphic * >( ),
                  atype< Matrice_Creuse< R > * >( ), atype< KN< R > * >( )),
      AF(linear_f), AG(P1_g) {}

  OptimIpopt(Case< linear_f, mv_P1_g >)
    : OneOperator(atype< long >( ), atype< KN< R > * >( ), atype< E_Array >( ),
                  atype< KN< R > * >( )),
      AF(linear_f), AG(mv_P1_g) {}

  OptimIpopt(Case< linear_f, linear_g >)
    : OneOperator(atype< long >( ), atype< KN< R > * >( ), atype< Matrice_Creuse< R > * >( ),
                  atype< KN< R > * >( )),
      AF(linear_f), AG(linear_g) {}
};

/*
 * enum AssumptionF {no_assumption_f, P2_f, unavailable_hessian, mv_P2_f, quadratic_f,linear_f};
 * enum AssumptionG {without_constraints, no_assumption_g, P1_g, mv_P1_g, linear_g};
 */

// Case dependant initialization of unused_name_param
// exemple : AF==no_assumption_f && AG==without_constraints ==> no constraint related named
// parameter should be used (index 2,3,4,6 - first integer is how many are not used)
void OptimIpopt::E_Ipopt::InitUNP( ) {
  if (AF == no_assumption_f && AG == no_assumption_g) {
  }    // no unused named parameter

  if (AF == no_assumption_f && AG == without_constraints) {
    AddElements(unused_name_param, 4, 2, 3, 4, 6);
  }

  if (AF == no_assumption_f && AG == P1_g) {
    AddElements(unused_name_param, 1, 4);
  }

  if (AF == no_assumption_f && AG == mv_P1_g) {
    AddElements(unused_name_param, 1, 4);
  }

  if (AF == no_assumption_f && AG == linear_g) {
    AddElements(unused_name_param, 1, 4);
  }

  if (AF == P2_f && AG == P1_g) {
    AddElements(unused_name_param, 5, 4, 5, 7, 8, 12);
  }

  if (AF == P2_f && AG == without_constraints) {
    AddElements(unused_name_param, 8, 2, 3, 4, 5, 6, 7, 8, 12);
  }

  if (AF == P2_f && AG == no_assumption_g) {
    AddElements(unused_name_param, 1, 5);
  }

  if (AF == P2_f && AG == mv_P1_g) {
    AddElements(unused_name_param, 5, 4, 5, 7, 8, 12);
  }

  if (AF == P2_f && AG == linear_g) {
    AddElements(unused_name_param, 5, 4, 5, 7, 8, 12);
  }

  if (AF == unavailable_hessian && AG == no_assumption_g) {
    AddElements(unused_name_param, 1, 5);
  }

  if (AF == unavailable_hessian && AG == without_constraints) {
    AddElements(unused_name_param, 7, 2, 3, 4, 5, 6, 7, 8);
  }

  if (AF == unavailable_hessian && AG == P1_g) {
    AddElements(unused_name_param, 4, 4, 5, 7, 8);
  }

  if (AF == unavailable_hessian && AG == mv_P1_g) {
    AddElements(unused_name_param, 4, 4, 5, 7, 8);
  }

  if (AF == unavailable_hessian && AG == linear_g) {
    AddElements(unused_name_param, 4, 4, 5, 7, 8);
  }

  if (AF == mv_P2_f && AG == without_constraints) {
    AddElements(unused_name_param, 8, 2, 3, 4, 5, 6, 7, 8, 12);
  }

  if (AF == mv_P2_f && AG == no_assumption_g) {
    AddElements(unused_name_param, 1, 5);
  }

  if (AF == mv_P2_f && AG == P1_g) {
    AddElements(unused_name_param, 5, 4, 5, 7, 8, 12);
  }

  if (AF == mv_P2_f && AG == mv_P1_g) {
    AddElements(unused_name_param, 5, 4, 5, 7, 8, 12);
  }

  if (AF == mv_P2_f && AG == linear_g) {
    AddElements(unused_name_param, 5, 4, 5, 7, 8, 12);
  }

  if (AF == quadratic_f && AG == without_constraints) {
    AddElements(unused_name_param, 8, 2, 3, 4, 5, 6, 7, 8, 12);
  }

  if (AF == quadratic_f && AG == no_assumption_g) {
    AddElements(unused_name_param, 1, 5);
  }

  if (AF == quadratic_f && AG == P1_g) {
    AddElements(unused_name_param, 5, 4, 5, 7, 8, 12);
  }

  if (AF == quadratic_f && AG == mv_P1_g) {
    AddElements(unused_name_param, 5, 4, 5, 7, 8, 12);
  }

  if (AF == quadratic_f && AG == linear_g) {
    AddElements(unused_name_param, 5, 4, 5, 7, 8, 12);
  }

  if (AF == linear_f && AG == without_constraints) {
    AddElements(unused_name_param, 8, 2, 3, 4, 5, 6, 7, 8, 12);
  }

  if (AF == linear_f && AG == no_assumption_g) {
    AddElements(unused_name_param, 1, 5);
  }

  if (AF == linear_f && AG == P1_g) {
    AddElements(unused_name_param, 5, 4, 5, 7, 8, 12);
  }

  if (AF == linear_f && AG == mv_P1_g) {
    AddElements(unused_name_param, 5, 4, 5, 7, 8, 12);
  }

  if (AF == linear_f && AG == linear_g) {
    AddElements(unused_name_param, 5, 4, 5, 7, 8, 12);
  }
}

basicAC_F0::name_and_type OptimIpopt::E_Ipopt::name_param[] = {
  // DONT CHANGE THE ORDER!!!! If some parameters need to be added, add them after the last one
  // otherwize warning message of this interface will be a mess
  {"lb", &typeid(KN_< double >)},      // 0  -  lower bound on optimization parameter X
  {"ub", &typeid(KN_< double >)},      // 1  -  upper bound on optimization parameter X
  {"clb", &typeid(KN_< double >)},     // 2  -  constraints lower bounds
  {"cub", &typeid(KN_< double >)},     // 3  -  constraints upper bounds
  {"structjacc", &typeid(E_Array)},    // 4  -  constraints jacobian structure
  {"structhess", &typeid(E_Array)},    // 5  -  lagrangian hessian structure
  {"lm", &typeid(KN_< double >)},      // 6  -  lagrange multiplier (for autostruct or to get their
                                       // value at the end of the algorithm)
  {"autostruct", &typeid(long)},       // 7  -  automatic structure determination
  {"checkindex", &typeid(bool)},       // 8  -  whether to use the FindIndex function or not
  {"tol", &typeid(double)},            // 9  -  stopping criteria tol
  {"maxiter", &typeid(long)},          // 10 -  stopping criteria : maximum number of iterations
  {"maxcputime", &typeid(double)},     // 11 -  stopping criteria : maximum CPU time
  {"bfgs", &typeid(bool)},             // 12 -  force the bfgs hessian approximation
  {"derivativetest", &typeid(string *)},    // 13 -  call the derivative checker
  {"optfile", &typeid(string *)},    // 14 -  set the ipopt option file name (default is ipopt.opt)
  {"printlevel", &typeid(long)},     // 15 -  controls IPOPT print level output
  {"dth", &typeid(double)},          // 16 -  perturbation for finite difference derivative test
  {"dttol", &typeid(double)},    // 17 -  relative tolerence for the derivative test error detection
  {"fixedvar", &typeid(string *)},      // 18 -  remove the equality simple bounds from problem
  {"warmstart", &typeid(bool)},         // 19 -  do we initialize multipliers with given values
  {"uz", &typeid(KN_< double >)},       // 20 -  simple upper bounds dual variable
  {"lz", &typeid(KN_< double >)},       // 21 -  simple lower bounds dual variable
  {"muinit", &typeid(double)},          // 22 -  barrier parameter initialization
  {"pivtol", &typeid(double)},          // 23 -  pivot tolerance for the linear solver
  {"brf", &typeid(double)},             // 24 -  bounds relax factor
  {"mustrategy", &typeid(string *)},    // 25 -  strategy for barrier parameter update
  {"objvalue", &typeid(double *)},      // 26 -  to get the last objective function value
  {"mumin", &typeid(double)},           // 27 -  minimal value for the barrier parameter
  {"linesearch",
   &typeid(bool)}    // 28 -  use the line search or not (if no, the usual Newton step is kept)
                     // {"osf",				&typeid(double) }
                     // //26 -  objective function scalling factor
};

static void Load_Init( ) {
  Global.Add("IPOPT", "(", new OptimIpopt(Case< no_assumption_f, no_assumption_g >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< no_assumption_f, without_constraints >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< no_assumption_f, mv_P1_g >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< no_assumption_f, linear_g >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< unavailable_hessian, no_assumption_g >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< unavailable_hessian, without_constraints >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< unavailable_hessian, mv_P1_g >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< unavailable_hessian, linear_g >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< mv_P2_f, no_assumption_g >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< mv_P2_f, without_constraints >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< mv_P2_f, mv_P1_g >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< mv_P2_f, linear_g >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< quadratic_f, no_assumption_g >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< quadratic_f, without_constraints >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< quadratic_f, mv_P1_g >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< quadratic_f, linear_g >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< linear_f, no_assumption_g >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< linear_f, without_constraints >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< linear_f, mv_P1_g >( )));
  Global.Add("IPOPT", "(", new OptimIpopt(Case< linear_f, linear_g >( )));
}

/*****************************************************************************************************************************
 *	Specialization of functions builders' constructor and operator()
 *****************************************************************************************************************************/

template<>
FitnessFunctionDatas< no_assumption_f >::FitnessFunctionDatas(const basicAC_F0 &args,
                                                              Expression const *nargs,
                                                              const C_F0 &theparam,
                                                              const C_F0 &objfact, const C_F0 &L_m)
  : GenericFitnessFunctionDatas( ) {
  const Polymorphic *opJ = dynamic_cast< const Polymorphic * >(args[0].LeftValue( )),
                    *opdJ = dynamic_cast< const Polymorphic * >(args[1].LeftValue( )),
                    *opH = dynamic_cast< const Polymorphic * >(args[2].LeftValue( ));
  ArrayOfaType hprototype2(atype< KN< R > * >( ), atype< double >( ), atype< KN< R > * >( )),
    hprototype1(atype< KN< R > * >( ));

  JJ = to< R >(C_F0(opJ, "(", theparam));
  GradJ = to< Rn_ >(C_F0(opdJ, "(", theparam));
  if (opH->Find("(", hprototype2)) {
    CompletelyNonLinearConstraints = true;
    Hessian = to< Matrice_Creuse< R > * >(C_F0(opH, "(", theparam, objfact, L_m));
  } else if (opH->Find("(", hprototype1)) {
    CompletelyNonLinearConstraints =
      false;    // When constraints are affine, lagrange multipliers are not used in the hessian,
                // obj_factor is also hidden to the user
    Hessian = to< Matrice_Creuse< R > * >(C_F0(opH, "(", theparam));
  } else {
    CompileError(
      "Error, wrong hessian function prototype. Must be either (real[int] &) or (real[int] "
      "&,real,real[int] &)");
  }
}

template<>
void FitnessFunctionDatas< no_assumption_f >::operator( )(Stack stack, const C_F0 &theparam,
                                                          const C_F0 &objfact, const C_F0 &L_m,
                                                          Expression const *nargs, ScalarFunc *&ffJ,
                                                          VectorFunc *&ffdJ, SparseMatFunc *&ffH,
                                                          bool warned) const {
  ffJ = new GeneralFunc< R >(stack, JJ, theparam);
  ffdJ = new GeneralFunc< Rn >(stack, GradJ, theparam);
  if (CompletelyNonLinearConstraints) {
    ffH = new GeneralSparseMatFunc(stack, Hessian, theparam, objfact, L_m);
  } else {
    ffH = new GeneralSparseMatFunc(stack, Hessian, theparam);
  }
}

template<>
FitnessFunctionDatas< P2_f >::FitnessFunctionDatas(const basicAC_F0 &args, Expression const *nargs,
                                                   const C_F0 &theparam, const C_F0 &objfact,
                                                   const C_F0 &L_m)
  : GenericFitnessFunctionDatas( ) {
  CompletelyNonLinearConstraints = false;
  const Polymorphic *opJ = dynamic_cast< const Polymorphic * >(args[0].LeftValue( )),
                    *opdJ = dynamic_cast< const Polymorphic * >(args[1].LeftValue( ));
  JJ = to< R >(C_F0(opJ, "(", theparam));
  GradJ = to< Rn_ >(C_F0(opdJ, "(", theparam));
  Hessian = to< Matrice_Creuse< R > * >(args[2]);
}

template<>
void FitnessFunctionDatas< P2_f >::operator( )(Stack stack, const C_F0 &theparam,
                                               const C_F0 &objfact, const C_F0 &L_m,
                                               Expression const *nargs, ScalarFunc *&ffJ,
                                               VectorFunc *&ffdJ, SparseMatFunc *&ffH,
                                               bool warned) const {
  if (warned && nargs[5]) {
    cout << "  ==> your lagrangian hessian is a constant matrix, so there is no need to specify "
            "its structure with ";
    cout << OptimIpopt::E_Ipopt::name_param[5].name << endl;
    cout << "      since it is contained in the matrix object." << endl;
  }

  ffJ = new GeneralFunc< R >(stack, JJ, theparam);
  ffdJ = new GeneralFunc< Rn >(stack, GradJ, theparam);
  ffH = new ConstantSparseMatFunc(stack, Hessian);
}

template<>
FitnessFunctionDatas< unavailable_hessian >::FitnessFunctionDatas(const basicAC_F0 &args,
                                                                  Expression const *nargs,
                                                                  const C_F0 &theparam,
                                                                  const C_F0 &objfact,
                                                                  const C_F0 &L_m)
  : GenericFitnessFunctionDatas( ) {
  CompletelyNonLinearConstraints = false;
  const Polymorphic *opJ = dynamic_cast< const Polymorphic * >(args[0].LeftValue( )),
                    *opdJ = dynamic_cast< const Polymorphic * >(args[1].LeftValue( ));
  JJ = to< R >(C_F0(opJ, "(", theparam));
  GradJ = to< Rn_ >(C_F0(opdJ, "(", theparam));
}

template<>
void FitnessFunctionDatas< unavailable_hessian >::operator( )(
  Stack stack, const C_F0 &theparam, const C_F0 &objfact, const C_F0 &L_m, Expression const *nargs,
  ScalarFunc *&ffJ, VectorFunc *&ffdJ, SparseMatFunc *&ffH, bool warned) const {
  if (warned && nargs[5]) {
    cout << "  ==> no hessian has been given, the LBFGS mode has been enabled, thus making ";
    cout << OptimIpopt::E_Ipopt::name_param[5].name << " useless. You may also" << endl
         << "      have forgoten a function (IPOPT will certainly crash if so)." << endl;
  }

  ffJ = new GeneralFunc< R >(stack, JJ, theparam);
  ffdJ = new GeneralFunc< Rn >(stack, GradJ, theparam);
  ffH = 0;
}

template<>
FitnessFunctionDatas< mv_P2_f >::FitnessFunctionDatas(const basicAC_F0 &args,
                                                      Expression const *nargs, const C_F0 &theparam,
                                                      const C_F0 &objfact, const C_F0 &L_m)
  : GenericFitnessFunctionDatas( ) {
  const E_Array *M_b = dynamic_cast< const E_Array * >(args[0].LeftValue( ));

  if (M_b->nbitem( ) != 2) {
    CompileError(
      "\nSorry, we were expecting an array with two componants, either [M,b] or [b,M] for the "
      "affine constraints expression.");
  }

  bool order = true;
  if (CheckMatrixVectorPair(M_b, order)) {
    Hessian = to< Matrice_Creuse< R > * >((*M_b)[order ? 0 : 1]);
    GradJ = to< Rn * >((*M_b)[order ? 1 : 0]);    // This is gradJ evaluated on x=0
  }
}

template<>
void FitnessFunctionDatas< mv_P2_f >::operator( )(Stack stack, const C_F0 &theparam,
                                                  const C_F0 &objfact, const C_F0 &L_m,
                                                  Expression const *nargs, ScalarFunc *&ffJ,
                                                  VectorFunc *&ffdJ, SparseMatFunc *&ffH,
                                                  bool warned) const {
  if (warned && nargs[5]) {
    cout << "  ==> your lagrangian hessian is a constant matrix, so there is no need to specify "
            "its structure with ";
    cout << OptimIpopt::E_Ipopt::name_param[5].name << endl;
    cout << "      since it is contained in the matrix object." << endl;
  }

  ffJ = new P2ScalarFunc(stack, Hessian, GradJ, true);
  ffdJ = new P1VectorFunc(stack, Hessian, GradJ, true);
  ffH = new ConstantSparseMatFunc(stack, Hessian);
}

template<>
FitnessFunctionDatas< quadratic_f >::FitnessFunctionDatas(const basicAC_F0 &args,
                                                          Expression const *nargs,
                                                          const C_F0 &theparam, const C_F0 &objfact,
                                                          const C_F0 &L_m)
  : GenericFitnessFunctionDatas( ) {
  Hessian = to< Matrice_Creuse< R > * >(args[0]);
}

template<>
void FitnessFunctionDatas< quadratic_f >::operator( )(Stack stack, const C_F0 &theparam,
                                                      const C_F0 &objfact, const C_F0 &L_m,
                                                      Expression const *nargs, ScalarFunc *&ffJ,
                                                      VectorFunc *&ffdJ, SparseMatFunc *&ffH,
                                                      bool warned) const {
  if (warned && nargs[5]) {
    cout << "  ==> your lagrangian hessian is a constant matrix, so there is no need to specify "
            "its structure with ";
    cout << OptimIpopt::E_Ipopt::name_param[5].name << endl;
    cout << "      since it is contained in the matrix object." << endl;
  }

  ffJ = new P2ScalarFunc(stack, Hessian, 0, true);
  ffdJ = new P1VectorFunc(stack, Hessian, 0, true);
  ffH = new ConstantSparseMatFunc(stack, Hessian);
}

template<>
FitnessFunctionDatas< linear_f >::FitnessFunctionDatas(const basicAC_F0 &args,
                                                       Expression const *nargs,
                                                       const C_F0 &theparam, const C_F0 &objfact,
                                                       const C_F0 &L_m)
  : GenericFitnessFunctionDatas( ) {
  GradJ = to< Rn * >(args[0]);
}

template<>
void FitnessFunctionDatas< linear_f >::operator( )(Stack stack, const C_F0 &theparam,
                                                   const C_F0 &objfact, const C_F0 &L_m,
                                                   Expression const *nargs, ScalarFunc *&ffJ,
                                                   VectorFunc *&ffdJ, SparseMatFunc *&ffH,
                                                   bool warned) const {
  if (warned && nargs[5]) {
    cout << "  ==> your lagrangian hessian is a null matrix, so there is no need to specify its "
            "structure with ";
    cout << OptimIpopt::E_Ipopt::name_param[5].name << endl;
    cout << "      since it is empty." << endl;
  }

  ffJ = new P2ScalarFunc(stack, 0, GradJ);
  ffdJ = new P1VectorFunc(stack, 0, GradJ);
  ffH = 0;
}

template<>
ConstraintFunctionDatas< without_constraints >::ConstraintFunctionDatas(const basicAC_F0 &args,
                                                                        Expression const *nargs,
                                                                        const C_F0 &theparam)
  : GenericConstraintFunctionDatas( ) {}

template<>
void ConstraintFunctionDatas< without_constraints >::operator( )(Stack stack, const C_F0 &theparam,
                                                                 Expression const *nargs,
                                                                 VectorFunc *&ffC,
                                                                 SparseMatFunc *&ffdC,
                                                                 bool warned) const {
  if (warned) {
    if (nargs[2] || nargs[3]) {
      cout << "  ==> Some constraints bounds have been defined while no constraints function has "
              "been passed."
           << endl;
    }

    if (nargs[4]) {
      cout << "  ==> A structure has been provided for the constraints jacobian but there is no "
              "constraint function."
           << endl;
    }

    if (nargs[6]) {
      cout << "  ==> Unconstrained problem make the use of "
           << OptimIpopt::E_Ipopt::name_param[6].name
           << " pointless (see the documentation for more details)." << endl;
    }
  }

  ffC = 0;
  ffdC = 0;
}

template<>
ConstraintFunctionDatas< no_assumption_g >::ConstraintFunctionDatas(const basicAC_F0 &args,
                                                                    Expression const *nargs,
                                                                    const C_F0 &theparam)
  : GenericConstraintFunctionDatas( ) {
  int nbj = args.size( ) - 1;
  const Polymorphic *opG = dynamic_cast< const Polymorphic * >(args[nbj - 2].LeftValue( )),
                    *opjG = dynamic_cast< const Polymorphic * >(args[nbj - 1].LeftValue( ));

  Constraints = to< Rn_ >(C_F0(opG, "(", theparam));
  GradConstraints = to< Matrice_Creuse< R > * >(C_F0(opjG, "(", theparam));
}

template<>
void ConstraintFunctionDatas< no_assumption_g >::operator( )(Stack stack, const C_F0 &theparam,
                                                             Expression const *nargs,
                                                             VectorFunc *&ffC, SparseMatFunc *&ffdC,
                                                             bool) const {
  ffC = new GeneralFunc< Rn >(stack, Constraints, theparam);
  ffdC = new GeneralSparseMatFunc(stack, GradConstraints, theparam);
}

template<>
ConstraintFunctionDatas< P1_g >::ConstraintFunctionDatas(const basicAC_F0 &args,
                                                         Expression const *nargs,
                                                         const C_F0 &theparam)
  : GenericConstraintFunctionDatas( ) {
  int nbj = args.size( ) - 1;
  const Polymorphic *opG = dynamic_cast< const Polymorphic * >(args[nbj - 2].LeftValue( ));

  Constraints = to< Rn_ >(C_F0(opG, "(", theparam));
  GradConstraints = to< Matrice_Creuse< R > * >(args[nbj - 1]);
}

template<>
void ConstraintFunctionDatas< P1_g >::operator( )(Stack stack, const C_F0 &theparam,
                                                  Expression const *nargs, VectorFunc *&ffC,
                                                  SparseMatFunc *&ffdC, bool warned) const {
  if (warned && nargs[4]) {
    cout << "  ==> your constraints jacobian is a constant matrix, there is no need to specify its "
            "structure with "
         << OptimIpopt::E_Ipopt::name_param[4].name << endl;
    cout << "      since it is contained in the matrix object." << endl;
  }

  ffC = new GeneralFunc< Rn >(stack, Constraints, theparam);
  ffdC = new ConstantSparseMatFunc(stack, GradConstraints);
}

template<>
ConstraintFunctionDatas< mv_P1_g >::ConstraintFunctionDatas(const basicAC_F0 &args,
                                                            Expression const *nargs,
                                                            const C_F0 &theparam)
  : GenericConstraintFunctionDatas( ) {
  int nbj = args.size( ) - 1;
  const E_Array *M_b = dynamic_cast< const E_Array * >(args[nbj - 1].LeftValue( ));

  if (M_b->nbitem( ) != 2) {
    CompileError(
      "\nSorry, we were expecting an array with two componants, either [M,b] or [b,M] for the "
      "affine constraints expression.");
  }

  bool order = true;
  if (CheckMatrixVectorPair(M_b, order)) {
    GradConstraints = to< Matrice_Creuse< R > * >((*M_b)[order ? 0 : 1]);
    Constraints = to< Rn * >((*M_b)[order ? 1 : 0]);    // Constraint on x=0
  } else {
    CompileError(
      "\nWrong types in the constraints [matrix,vector] pair, expecting a sparse matrix and "
      "real[int].");
  }
}

template<>
void ConstraintFunctionDatas< mv_P1_g >::operator( )(Stack stack, const C_F0 &theparam,
                                                     Expression const *nargs, VectorFunc *&ffC,
                                                     SparseMatFunc *&ffdC, bool warned) const {
  if (warned && nargs[4]) {
    cout << "  ==> your constraints jacobian is a constant matrix, there is no need to specify its "
            "structure with "
         << OptimIpopt::E_Ipopt::name_param[4].name << endl;
    cout << "      since it is contained in the matrix object." << endl;
  }

  ffC = new P1VectorFunc(stack, GradConstraints, Constraints);
  ffdC = new ConstantSparseMatFunc(stack, GradConstraints);
}

template<>
ConstraintFunctionDatas< linear_g >::ConstraintFunctionDatas(const basicAC_F0 &args,
                                                             Expression const *nargs,
                                                             const C_F0 &theparam)
  : GenericConstraintFunctionDatas( ) {
  int nbj = args.size( ) - 1;

  GradConstraints = to< Matrice_Creuse< R > * >(args[nbj - 1]);
}

template<>
void ConstraintFunctionDatas< linear_g >::operator( )(Stack stack, const C_F0 &theparam,
                                                      Expression const *nargs, VectorFunc *&ffC,
                                                      SparseMatFunc *&ffdC, bool warned) const {
  if (warned && nargs[4]) {
    cout << "  ==> your constraints jacobian is a constant matrix, there is no need to specify its "
            "structure with "
         << OptimIpopt::E_Ipopt::name_param[4].name << endl;
    cout << "      since it is contained in the matrix object." << endl;
  }

  ffC = new P1VectorFunc(stack, GradConstraints, 0);
  ffdC = new ConstantSparseMatFunc(stack, GradConstraints);
}

GenericFitnessFunctionDatas *GenericFitnessFunctionDatas::New(AssumptionF AF,
                                                              const basicAC_F0 &args,
                                                              Expression const *nargs,
                                                              const C_F0 &theparam,
                                                              const C_F0 &objfact, const C_F0 &lm) {
  switch (AF) {
    case no_assumption_f:
      return new FitnessFunctionDatas< no_assumption_f >(args, nargs, theparam, objfact, lm);
      break;
    case P2_f:
      return new FitnessFunctionDatas< P2_f >(args, nargs, theparam, objfact, lm);
      break;
    case unavailable_hessian:
      return new FitnessFunctionDatas< unavailable_hessian >(args, nargs, theparam, objfact, lm);
      break;
    case mv_P2_f:
      return new FitnessFunctionDatas< mv_P2_f >(args, nargs, theparam, objfact, lm);
      break;
    case quadratic_f:
      return new FitnessFunctionDatas< quadratic_f >(args, nargs, theparam, objfact, lm);
      break;
    case linear_f:
      return new FitnessFunctionDatas< linear_f >(args, nargs, theparam, objfact, lm);
      break;
    default:
      return 0;
      break;
  }
}

GenericConstraintFunctionDatas *GenericConstraintFunctionDatas::New(AssumptionG AG,
                                                                    const basicAC_F0 &args,
                                                                    Expression const *nargs,
                                                                    const C_F0 &theparam) {
  switch (AG) {
    case no_assumption_g:
      return new ConstraintFunctionDatas< no_assumption_g >(args, nargs, theparam);
      break;
    case without_constraints:
      return new ConstraintFunctionDatas< without_constraints >(args, nargs, theparam);
      break;
    case P1_g:
      return new ConstraintFunctionDatas< P1_g >(args, nargs, theparam);
      break;
    case mv_P1_g:
      return new ConstraintFunctionDatas< mv_P1_g >(args, nargs, theparam);
      break;
    case linear_g:
      return new ConstraintFunctionDatas< linear_g >(args, nargs, theparam);
      break;
    default:
      return 0;
      break;
  }
}

LOADFUNC(Load_Init)
