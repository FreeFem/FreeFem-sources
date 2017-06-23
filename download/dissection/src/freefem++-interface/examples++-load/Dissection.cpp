//   for automatic  compilation with ff-c++
//ff-c++-LIBRARY-dep:  dissection blas pthread
//ff-c++-cpp-dep: 
//  
//  file to add Dissection solver with dynamic load.
#include  <iostream>
#include <list>
using namespace std;

#include "rgraph.hpp"
#include "error.hpp"
#include "AFunction.hpp"


#include "MatriceCreuse_tpl.hpp"

#include <wchar.h>

#include "Dissection.hpp"

template<typename T>
int generate_CSR(list<int>* ind_cols_tmp, list<T>* val_tmp, 
		 int nrow, int nrow1, int *old2new, int *new2old,
		 int *ptrow, int *ind_col, T* val)
{
  int nnz = 0;
  //  ind_cols_tmp = new list<int>[nrow];
  //  val_tmp = new list<T>[nrow];
  for (int i = 0; i < nrow1; i++) {
    const int ii = new2old[i];
    bool diag_flag = false;
    for (int k = ptrow[ii]; k < ptrow[ii + 1]; k++) {
      const int j = old2new[ind_col[k]];
      if (i == j) {
	diag_flag = true;
      }
    //    fprintf(stderr, "%d %d -> %d %d \n", i0, j0, ii, jj);
      if (ind_cols_tmp[i].empty()) {
	ind_cols_tmp[i].push_back(j);
	val_tmp[i].push_back(val[k]);
	nnz++;
      }
      else {
	if (ind_cols_tmp[i].back() < j) {
	  ind_cols_tmp[i].push_back(j);
	  val_tmp[i].push_back(val[k]);
	  nnz++;
	}
	else {
	  typename list<T>::iterator iv = val_tmp[i].begin();
	  list<int>::iterator it = ind_cols_tmp[i].begin();
	  for ( ; it != ind_cols_tmp[i].end(); ++it, ++iv) {
	    if (*it == j) {
	      fprintf(stderr, "already exits? (%d %d)\n", ii, j);
	      break;
	    }
	    if (*it > j) {
	      ind_cols_tmp[i].insert(it, j);
	      val_tmp[i].insert(iv, val[k]);
	      nnz++;
	      break;
	    }
	  }
	}
      }
    } // loop : k
    if (!diag_flag) {
      fprintf(stderr, "%s %d : adding zero-entry %d\n",
	      __FILE__, __LINE__, i);
      typename list<T>::iterator iv = val_tmp[i].begin();
      list<int>::iterator it = ind_cols_tmp[i].begin();
      for ( ; it != ind_cols_tmp[i].end(); ++it, ++iv) {
	if ((*it) > i) {
	  ind_cols_tmp[i].insert(it, i);
	  val_tmp[i].insert(iv, 0.0);
	  nnz++;
	  break;
	}
      } // loop : iv
    } // if (!diag_flag)
  } // loop : i
  return nnz;
}

template
int generate_CSR<double>(list<int>* ind_cols_tmp, list<double>* val_tmp, 
			 int nrow, int nrow1, int *old2new, int *new2old,
			 int *ptrow, int *ind_col, double * val);

template
int generate_CSR<Complex>(list<int>* ind_cols_tmp,
			  list<Complex>* val_tmp, 
			  int nrow, int nrow1, int *old2new,
			  int *new2old,
			  int *ptrow, int *ind_col,
			  Complex * val);

template<typename T>
int copy_CSR(int *ptrows, int *indcols, T* coefs, int nrow, 
	     list<int>* ind_cols_tmp, list<T>* val_tmp,
	     int lower, int upper)
{
  ptrows[0] = 0;
  int k = 0;
  for (int i = 0; i < nrow; i++) {
    list<int>::iterator it = ind_cols_tmp[i].begin();
    typename list<T>::iterator iv = val_tmp[i].begin();
    for ( ; it != ind_cols_tmp[i].end(); ++it, ++iv) {
      if ((*it) >= lower && (*it) < upper) {
	indcols[k] = *it;
	coefs[k] = *iv;
	k++;
      }
    }
    ptrows[i + 1] = k;
  } // loop : i
  return k;
}

template
int copy_CSR<double>(int *ptrows, int *indcols, double* coefs, int nrow, 
		     list<int>* ind_cols_tmp, list<double>* val_tmp,
		     int lower, int upper);
template
int copy_CSR<Complex>(int *ptrows, int *indcols,
		      Complex* coefs, int nrow, 
		      list<int>* ind_cols_tmp,
		      list<Complex>* val_tmp,
		      int lower, int upper);

template<class R>
class SolveDissection:   public MatriceMorse<R>::VirtualSolver  {
  double _tgv;
  double _pivot;
  int _dim0;
  int *_new2old;
  int *_ptrows0;
  int *_indcols0;
  double *_coefs0;
  int *_ptrows1;
  int *_indcols1;
  double *_coefs1;
  double *_xtmp;
  uint64_t *_dslv;
public:

  SolveDissection(const MatriceMorse<R> &A,
		  int strategy_, double ttgv, 
		  double pivot) :
    _tgv(ttgv),
    _pivot(pivot)
  {
    const int dim = A.n;
   int num_threads = 1;
   if (getenv("DISSECTION_NUM_THREADS")) {
      sscanf(getenv("DISSECTION_NUM_THREADS"), "%d", &num_threads);
    fprintf(stderr,
	    "environmental variable DISSECTION_NUM_THREADS = %d\n",
	    num_threads);
   }
    _dslv = new uint64_t;
    diss_init(*_dslv, 0, 1, 0, num_threads, ((verbosity > 0) ? 1 : 0)); 
    // real matrix, with double precision factorization, # of threads
    const int decomposer = (strategy_ % 100) / 10;
    const int scaling = strategy_ == 0 ? 2 : strategy_ % 10;
    // sym + lower + isWhole = 1 + 2 + 4
    const int sym = (strategy_ / 100) ? 5 : 0;

    int *old2new = new int[dim];
    _new2old = new int[dim];
    {
      int m = 0;
      for (int i = 0; i < dim; i++) {
	for (int k = (A.lg)[i]; k < (A.lg)[i + 1]; k++) {
	  const int j = (A.cl)[k];
	  if (i == j) {
	    if ((A.a)[k] != ttgv) {
	      _new2old[m] = i;
	      m++;
	    }
	    break;
	  }
	}
      }
      _dim0 = m;
      for (int i = 0; i < dim; i++) {
	for (int k = (A.lg)[i]; k < (A.lg)[i + 1]; k++) {
	  const int j = (A.cl)[k];
	  if (i == j) {
	    if ((A.a)[k] == ttgv) {
	      _new2old[m] = i;
	      m++;
	    }
	    break;
	  }
	}
      }
      if (verbosity > 10) {
	cout << "m = " << m << " dim = " << dim << endl;
      }
    }
    for (int i = 0; i < dim; i++) {
      old2new[_new2old[i]] = i;
    }
    //    
    _xtmp = new double[_dim0];
    int nnz;
    list<int> *indcols_tmp = new list<int>[_dim0];
    list<double> *coefs_tmp = new list<double>[_dim0];
    nnz = generate_CSR<double>(indcols_tmp, coefs_tmp, dim, _dim0, old2new,
			       _new2old,
			       (int *)A.lg, (int *)A.cl, (double *)A.a);
    delete [] old2new;
    if (verbosity > 10) {
      cout << "nnz " << nnz << endl;
    }
    _ptrows0 = new int[_dim0 + 1];
    _ptrows1 = new int[_dim0 + 1];
    int *indcol_tmp = new int[nnz];
    double *coef_tmp = new double[nnz];
    int nnz1;
    nnz1 = copy_CSR<double>(_ptrows0, indcol_tmp, coef_tmp, _dim0, indcols_tmp, 
			    coefs_tmp, 0, _dim0);
    if (verbosity > 0) {
      cout << "nnz1 " << nnz1 << endl;
    }
    _indcols0 = new int[nnz1];
    _coefs0 = new double[nnz1];
    for (int i = 0; i < nnz1; i++) {
      _indcols0[i] = indcol_tmp[i];
      _coefs0[i] = coef_tmp[i];
    }
    nnz1 = copy_CSR<double>(_ptrows1, indcol_tmp, coef_tmp, _dim0, 
			    indcols_tmp, coefs_tmp, _dim0, dim);
    if (verbosity > 10) {
      cout << "nnz1 " << nnz1 << endl;
    }
    _indcols1 = new int[nnz1];
    _coefs1 = new double[nnz1];
    for (int i = 0; i < nnz1; i++) {
      _indcols1[i] = indcol_tmp[i];
      _coefs1[i] = coef_tmp[i];
    }

    delete [] indcols_tmp;
    delete [] coefs_tmp;
    delete [] indcol_tmp;
    delete [] coef_tmp;
    diss_s_fact(*_dslv, _dim0, _ptrows0, _indcols0, sym, decomposer);
    const int indefinite_flag = 1;
    const double eps_pivot = _pivot == (-1.0) ? 1.0e-2 : _pivot;
    diss_n_fact(*_dslv, _coefs0, scaling, eps_pivot, 
		indefinite_flag);
    int n0;
    diss_get_kern_dim(*_dslv, &n0);
    if (n0 > 0) {
      cout << "the matrix with size = " << _dim0 << " is singular with "
	   << n0 << " dimenisonal kernel." 
	   << endl;
    }
  }

  void Solver(const MatriceMorse<R> &A,KN_<R> &x,const KN_<R> &b) const {
    ffassert (&x[0] != &b[0]);
    int dim = A.n;
    // #x_1[] = _dim0; x_1 = b_1 - A_12 x_2 
    for (int i = 0; i < _dim0; i++) {
      _xtmp[i] = b[_new2old[i]];
    }
    for (int i = _dim0; i < dim; i++) {
      const int ii = _new2old[i];
      x[ii] = b[ii] / _tgv;
    }
    for (int i = 0; i < _dim0; i++) {
      for (int k = _ptrows1[i]; k < _ptrows1[i + 1]; k++) {
	_xtmp[i] -= _coefs1[k] * x[_new2old[_indcols1[k]]];
      }
    }
    const int projection = 1;
    const int transpose = 0;
    diss_solve_1(*_dslv, _xtmp, projection, transpose);
    for (int i = 0; i < _dim0; i++) {
      x[_new2old[i]] = _xtmp[i];
    }
  }

  ~SolveDissection() {
    diss_free(*_dslv);
    delete _dslv;
    delete [] _new2old;
    delete [] _ptrows0;
    delete [] _indcols0;
    delete [] _coefs0;
    delete [] _ptrows1;
    delete [] _indcols1;
    delete [] _coefs1;
    delete [] _xtmp;
  }

  void addMatMul(const KN_<R> & x, KN_<R> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax += (const MatriceMorse<R> &) (*this) * x; 
  }

}; 


template<>
class SolveDissection<Complex> : public MatriceMorse<Complex>::VirtualSolver {
  // double eps;
  //  mutable double  epsr;
  double _tgv;
  double _pivot; //, _pivot_sym;
  int _dim0;
  int *_new2old;
  int *_ptrows0;
  int *_indcols0;
  Complex *_coefs0;
  int *_ptrows1;
  int *_indcols1;
  Complex *_coefs1;
  Complex *_xtmp;
  uint64_t *_dslv;
public:
  SolveDissection(const MatriceMorse<Complex> &A,
		  int strategy_, double ttgv, 
		  double pivot) : 
    _tgv(ttgv),
    _pivot(pivot)
  {
    const int dim = A.n;
    int num_threads = 1;
    if (getenv("DISSECTION_NUM_THREADS")) {
      sscanf(getenv("DISSECTION_NUM_THREADS"), "%d", &num_threads);
      fprintf(stderr,
	      "environmental variable DISSECTION_NUM_THREADS = %d\n",
	      num_threads);
    }

    _dslv = new uint64_t;
    diss_init(*_dslv, 0, 2, 0, num_threads, ((verbosity > 0) ? 1 : 0)); // 
    // complex matrix, with double precision factorization, # of threads
    const int decomposer = (strategy_ % 100) / 10;
    const int scaling = strategy_ == 0 ? 2 : strategy_ % 10;
    // sym + lower + isWhole = 1 + 2 + 4
    const int sym = (strategy_ / 100) ? 5 : 0;

    int *old2new = new int[dim];
    _new2old = new int[dim];
    {
      int m = 0;
      for (int i = 0; i < dim; i++) {
	for (int k = (A.lg)[i]; k < (A.lg)[i + 1]; k++) {
	  const int j = (A.cl)[k];
	  if (i == j) {
	    if ((A.a)[k] != ttgv) {
	      _new2old[m] = i;
	      m++;
	    }
	    break;
	  }
	}
      }
      _dim0 = m;
      for (int i = 0; i < dim; i++) {
	for (int k = (A.lg)[i]; k < (A.lg)[i + 1]; k++) {
	  const int j = (A.cl)[k];
	  if (i == j) {
	    if ((A.a)[k] == ttgv) {
	      _new2old[m] = i;
	      m++;
	    }
	    break;
	  }
	}
      }
      if (verbosity > 10) {
	cout << "m = " << m << " dim = " << dim << endl;
      }
    }
    for (int i = 0; i < dim; i++) {
      old2new[_new2old[i]] = i;
    }
    //    
    _xtmp = new Complex[_dim0];
    int nnz;
    list<int> *indcols_tmp = new list<int>[_dim0];
    list<Complex> *coefs_tmp = new list<Complex>[_dim0];
    nnz = generate_CSR<Complex>(indcols_tmp, coefs_tmp, dim, _dim0, old2new,
				_new2old,
				(int *)A.lg, (int *)A.cl, (Complex *)A.a);
    delete [] old2new;
    if (verbosity > 10) {
      cout << "nnz " << nnz << endl;
    }
    _ptrows0 = new int[_dim0 + 1];
    _ptrows1 = new int[_dim0 + 1];
    int *indcol_tmp = new int[nnz];
    Complex *coef_tmp = new Complex[nnz];
    int nnz1;
    nnz1 = copy_CSR<Complex>(_ptrows0, indcol_tmp, coef_tmp, _dim0,
			     indcols_tmp, coefs_tmp, 0, _dim0);
    if (verbosity > 10) {
      cout << "nnz1 " << nnz1 << endl;
    }
    _indcols0 = new int[nnz1];
    _coefs0 = new Complex[nnz1];
    for (int i = 0; i < nnz1; i++) {
      _indcols0[i] = indcol_tmp[i];
      _coefs0[i] = coef_tmp[i];
    }
    nnz1 = copy_CSR<Complex>(_ptrows1, indcol_tmp, coef_tmp, _dim0, 
			     indcols_tmp, coefs_tmp, _dim0, dim);
    if (verbosity > 10) {
      cout << "nnz1 " << nnz1 << endl;
    }
    _indcols1 = new int[nnz1];
    _coefs1 = new Complex[nnz1];
    for (int i = 0; i < nnz1; i++) {
      _indcols1[i] = indcol_tmp[i];
      _coefs1[i] = coef_tmp[i];
    }
    delete [] indcols_tmp;
    delete [] coefs_tmp;
    delete [] indcol_tmp;
    delete [] coef_tmp;
    diss_s_fact(*_dslv, _dim0, _ptrows0, _indcols0, sym, decomposer);
    const int indefinite_flag = 1;
    const double eps_pivot = _pivot == (-1.0) ? 1.0e-2 : _pivot;
    diss_n_fact(*_dslv, (double *)_coefs0, scaling, eps_pivot, 
		indefinite_flag);
    int n0;
    diss_get_kern_dim(*_dslv, &n0);
    if (n0 > 0) {
      cout << "the matrix with size = " << _dim0 << " is singular with "
	   << n0 << " dimenisonal kernel." 
	   << endl;
    }
  }

  void Solver(const MatriceMorse<Complex> &A,KN_<Complex> &x,const KN_<Complex> &b) const  { 
    ffassert (&x[0] != &b[0]);
    int dim = A.n;
    // #x_1[] = _dim0; x_1 = b_1 - A_12 x_2 
    for (int i = 0; i < _dim0; i++) {
      _xtmp[i] = b[_new2old[i]];
    }
    for (int i = _dim0; i < dim; i++) {
      const int ii = _new2old[i];
      x[ii] = b[ii] / _tgv;
    }
    for (int i = 0; i < _dim0; i++) {
      for (int k = _ptrows1[i]; k < _ptrows1[i + 1]; k++) {
	_xtmp[i] -= _coefs1[k] * x[_new2old[_indcols1[k]]];
      }
    }
    const int projection = 1;
    const int transpose = 0;
    diss_solve_1(*_dslv, (double *)_xtmp, projection, transpose);
    for (int i = 0; i < _dim0; i++) {
      x[_new2old[i]] = _xtmp[i];
    }
  }

  ~SolveDissection() {
    diss_free(*_dslv);
    delete _dslv;
  }

  void addMatMul(const KN_<Complex> & x, KN_<Complex> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax += (const MatriceMorse<Complex> &) (*this) * x; 
  }

}; 

inline MatriceMorse<double>::VirtualSolver *
BuildSolverIDissection(DCL_ARG_SPARSE_SOLVER(double,A))
{
  if( verbosity>9)
    cout << " BuildSolverDissection<double>" << endl;
  return new SolveDissection<double>(*A,
				     ds.strategy,ds.tgv,ds.tol_pivot);
}

inline MatriceMorse<Complex>::VirtualSolver *
BuildSolverIDissection(DCL_ARG_SPARSE_SOLVER(Complex,A))
{
  if( verbosity>9)
    cout << " BuildSolverDissection<Complex>" << endl;
  return new SolveDissection<Complex>(*A,
				      ds.strategy,
				      ds.tgv,ds.tol_pivot);
}


//  the 2 default sparse solver double and complex
DefSparseSolver<double>::SparseMatSolver SparseMatSolver_R ;
DefSparseSolver<Complex>::SparseMatSolver SparseMatSolver_C;
DefSparseSolverSym<double>::SparseMatSolver SparseMatSolverSym_R ;
DefSparseSolverSym<Complex>::SparseMatSolver SparseMatSolverSym_C;
// the default probleme solver 
TypeSolveMat::TSolveMat  TypeSolveMatdefaultvalue=TypeSolveMat::defaultvalue;

bool SetDissection()
{
    if(verbosity>1)
	cout << " SetDefault sparse solver to Dissection" << endl;
    DefSparseSolver<double>::solver  =BuildSolverIDissection;
    DefSparseSolver<Complex>::solver =BuildSolverIDissection;    
    DefSparseSolverSym<double>::solver  =BuildSolverIDissection;
    DefSparseSolverSym<Complex>::solver =BuildSolverIDissection;    
    TypeSolveMat::defaultvalue =TypeSolveMatdefaultvalue;
    return  true;
}


void init22()
{    
  SparseMatSolver_R= DefSparseSolver<double>::solver;
  SparseMatSolver_C= DefSparseSolver<Complex>::solver;
  SparseMatSolverSym_R= DefSparseSolverSym<double>::solver;
  SparseMatSolverSym_C= DefSparseSolverSym<Complex>::solver;

  if(verbosity>1)
    cout << "\n Add: Dissection:  defaultsolver defaultsolverDissection" << endl;
  TypeSolveMat::defaultvalue=TypeSolveMat::SparseSolver; 
  DefSparseSolver<double>::solver =BuildSolverIDissection;
  DefSparseSolver<Complex>::solver =BuildSolverIDissection;

  if(! Global.Find("defaulttoDissection").NotNull() )
    Global.Add("defaulttoDissection","(",new OneOperator0<bool>(SetDissection));  
}


LOADFUNC(init22);
