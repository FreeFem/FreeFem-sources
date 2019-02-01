 /*! \file DissectionSolver.cpp
    \brief task mangemanet of dissection algorithm
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Mar. 30th 2012
    \date   Jul. 12th 2015
    \date   Nov. 30th 2016
*/

// This file is part of Dissection
// 
// Dissection is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Linking Dissection statically or dynamically with other modules is making
// a combined work based on Disssection. Thus, the terms and conditions of 
// the GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of Dissection give you 
// permission to combine Dissection program with free software programs or 
// libraries that are released under the GNU LGPL and with independent modules 
// that communicate with Dissection solely through the Dissection-fortran 
// interface. You may copy and distribute such a system following the terms of 
// the GNU GPL for Dissection and the licenses of the other code concerned, 
// provided that you include the source code of that other code when and as
// the GNU GPL requires distribution of source code and provided that you do 
// not modify the Dissection-fortran interface.
//
// Note that people who make modified versions of Dissection are not obligated 
// to grant this special exception for their modified versions; it is their
// choice whether to do so. The GNU General Public License gives permission to 
// release a modified version without this exception; this exception also makes
// it possible to release a modified version which carries forward this
// exception. If you modify the Dissection-fortran interface, this exception 
// does not apply to your modified version of Dissection, and you must remove 
// this exception when you distribute your modified version.
//
// This exception is an additional permission under section 7 of the GNU 
// General Public License, version 3 ("GPLv3")
//
// Dissection is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Dissection.  If not, see <http://www.gnu.org/licenses/>.
//

#include <float.h>
#include <vector>
#include <algorithm>

#include "Driver/DissectionSolver.hpp"
#include "Driver/C_BlasRoutines.hpp"
#include "Driver/C_KernDetect.hpp"
#include "Driver/CopyMatrix.hpp"
#include "Driver/TridiagBlockMatrix.hpp"
#include "Splitters/MetisSplitter.hpp"
#include "Algebra/SparseRenumbering.hpp"
#include "Algebra/VectorArray.hpp"
#include "Compiler/DissectionIO.hpp"

#define NORMALIZE_KERNEL_BASIS

template<typename T, typename U, typename W, typename Z>
const T DissectionSolver<T, U, W, Z>::_one = T(1.0);
template<typename T, typename U, typename W, typename Z>
const T DissectionSolver<T, U, W, Z>::_none = T(-1.0);
template<typename T, typename U, typename W, typename Z>
const T DissectionSolver<T, U, W, Z>::_zero = T(0.0);

template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
Destroy(void) 
{
  for (int m = 0; m < _graph_colors;m++) {
    if (!_tridiagQueue[m]->tridiagSolver()) {
      delete _dissectionQueue[m];
      for (typename vector<DissectionMatrix<T, U>* >::const_iterator it =
	     _dissectionMatrix[m].begin(); 
	   it != _dissectionMatrix[m].end(); ++it) {
	delete (*it);
      }  
    }
    else {
      delete _tridiagMatrix[m];
    }
  }
  delete [] _dissectionQueue;
  delete [] _dissectionMatrix;
  delete [] _tridiagMatrix;

  delete _ptDA;
  delete[] _precDiag; 

  NumericFree();
  delete [] _Schur;
  delete [] _kernel;
  delete [] _singIdx;       // added 01 Oct.2013 Atsushi
  for (int m = 0; m < _graph_colors;m++) {
    if (_with_btree && !_tridiagQueue[m]->tridiagSolver()) {
      delete _btree[m];               // added 01 Oct.2013 Atsushi
    }
    delete _tridiagQueue[m];
  }
  delete [] _tridiagQueue;
  delete [] _btree;
  _index_isolated.clear();
} 

template
void DissectionSolver<double>::Destroy(void);

template
void DissectionSolver<quadruple>::Destroy(void);

template
void DissectionSolver<quadruple, quadruple, double, double>::Destroy(void);

template
void DissectionSolver<double, double, quadruple, quadruple>::Destroy(void);

template
void DissectionSolver<complex<double>, double>::Destroy(void);

template
void DissectionSolver<complex<quadruple>, quadruple>::Destroy(void);

template
void DissectionSolver<complex<quadruple>, quadruple,
		      complex<double>, double>::Destroy(void);

template
void DissectionSolver<complex<double>, double,
		      complex<quadruple>, quadruple>::Destroy(void);

template
void DissectionSolver<float>::Destroy(void);

template
void DissectionSolver<complex<float>, float>::Destroy(void);

//

template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
NumericFree(void) 
{
  diss_printf(_verbose, _fp,
	      "%s %d : void DissectionSolver::NumericFree()",
	      __FILE__, __LINE__ );
  if (_status_factorized) {
    diss_printf(_verbose, _fp, "start ");
    for (int m = 0; m < _graph_colors; m++) {
      _Schur[m].free();
      _kernel[m].free();
      _singIdx[m].clear(); // added 01 Oct.2013 Atsushi
    }
    diss_printf(_verbose, _fp, "end");
  }
  _status_factorized = false;
  diss_printf(_verbose, _fp, ".\n");
} 

template
void DissectionSolver<double>::NumericFree(void);

template
void DissectionSolver<quadruple>::NumericFree(void);

template
void DissectionSolver<quadruple, quadruple, double, double>::NumericFree(void);

template
void DissectionSolver<double, double, quadruple, quadruple>::NumericFree(void);

template
void DissectionSolver<complex<double>, double>::NumericFree(void);

template
void DissectionSolver<complex<quadruple>, quadruple>::NumericFree(void);
//

template
void DissectionSolver<complex<quadruple>, quadruple, complex<double>, double>::NumericFree(void);

template
void DissectionSolver<complex<double>, double, complex<quadruple>, quadruple>::NumericFree(void);

template
void DissectionSolver<float>::NumericFree(void);

template
void DissectionSolver<complex<float>, float>::NumericFree(void);

//

template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
SaveCSRMatrix(const int called,	const T *coefs_)
{
  diss_printf(true, stderr,
	      "%s %d : specialized template is not yet defined.\n",
	      __FILE__, __LINE__);
}

template
void DissectionSolver<double>::SaveCSRMatrix(const int called,
					     const double *coefs);

template
void DissectionSolver<complex<double>, double>::
SaveCSRMatrix(const int called,
	      const complex<double> *coefs);

template
void DissectionSolver<quadruple>::
SaveCSRMatrix(const int called,
	      const quadruple *coefs);

template
void DissectionSolver<complex<quadruple>, quadruple>::
SaveCSRMatrix(const int called,
	      const complex<quadruple> *coefs);

template
void DissectionSolver<float>::SaveCSRMatrix(const int called,
					     const float *coefs);

template
void DissectionSolver<complex<float>, float>::
SaveCSRMatrix(const int called,
	      const complex<float> *coefs);
//

template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
SaveMMMatrix(const int called, const T *coefs_)
{
  diss_printf(true, stderr,
	      "%s %d : specialized template is not yet defined.\n",
	      __FILE__, __LINE__);
}


template<>
void DissectionSolver<double>::SaveMMMatrix(const int called,
					    const double *coefs_)
{
  SaveMMMatrix_(_dim,
		_ptDA->nnz(),
		_ptDA->isSymmetric(),
		_ptDA->isUpper(),
		_ptDA->getRows(),
		_ptDA->getIndCols(),
		called, coefs_);
  
}

template<>
void DissectionSolver<complex<double>, double>::
SaveMMMatrix(const int called,
	      const complex<double> *coefs_)
{
  SaveMMMatrix_(_dim,
		_ptDA->nnz(),
		_ptDA->isWhole() ? false : _ptDA->isSymmetric(),
		_ptDA->isUpper(),
		_ptDA->getRows(),
		_ptDA->getIndCols(),
		called, coefs_);
}

template
void DissectionSolver<quadruple>::
SaveMMMatrix(const int called,
	     const quadruple *coefs_);
template
void DissectionSolver<complex<quadruple>, quadruple>::
SaveMMMatrix(const int called,
	     const complex<quadruple> *coefs_);
//

void SaveMMMatrix_(const int dim,
		   const int nnz,
		   const bool isSymmetric,
		   const bool isUpper,
		   const int *ptrows,
		   const int *indcols,
		   const int called,
		   const double *coefs_)
  
{
 // MatrixMarket format : one-based index, lower part stored for symmetric
  char fname[256];
  int pid = get_process_id();
  FILE *fp;
  sprintf(fname, "matrix.%04d.%06d.data", called, pid);
  if ((fp = fopen(fname, "w")) != NULL) {
    fprintf(fp, "%%%%MatrixMarket matrix coordinate real %s\n%d %d %d\n",
	    
	    isSymmetric? "symmetric" : "general",
	    dim, dim, nnz);
    if (isUpper) {
      for (int i = 0; i < dim; i++) {
	for (int k = ptrows[i]; k < ptrows[i + 1]; k++) {
	  const int jj = indcols[k] + 1;
	  fprintf(fp, "%d %d %.16e\n", jj, (i + 1), coefs_[k]);
	}
      }
    }
    else {
      for (int i = 0; i < dim; i++) {
	for (int k = ptrows[i]; k < ptrows[i + 1]; k++) {
	  const int jj = indcols[k] + 1;
	  fprintf(fp, "%d %d %.16e\n", (i + 1), jj, coefs_[k]);
	}
      }
    }
    fclose(fp);
  }
  else {
    fprintf(stderr,
	    "%s %d : fail to open %s\n",
	    __FILE__, __LINE__, fname);
    exit(-1);
  }
}


void SaveMMMatrix_(const int dim,
		   const int nnz,
		   const bool isSymmetric,
		   const bool isUpper,
		   const int *ptrows,
		   const int *indcols,
		   const int called,
		   const complex<double> *coefs_)
{
 // MatrixMarket format : one-based index, lower part stored for symmetric
  char fname[256];
  int pid = get_process_id();
  FILE *fp;
  sprintf(fname, "matrix.%04d.%06d.data", called, pid);
  if ((fp = fopen(fname, "w")) != NULL) {
    fprintf(fp, "%%%%MatrixMarket matrix coordinate complex %s\n%d %d %d\n",
	    
	    isSymmetric? "symmetric" : "general",
	    dim, dim, nnz);
    if (isUpper) {
      for (int i = 0; i < dim; i++) {
	for (int k = ptrows[i]; k < ptrows[i + 1]; k++) {
	  const int jj = indcols[k] + 1;
	  fprintf(fp, "%d %d %.16e %.16e\n", jj, (i + 1), 
		  coefs_[k].real(), coefs_[k].imag());
	}
      }
    }
    else {
      for (int i = 0; i < dim; i++) {
	for (int k = ptrows[i]; k < ptrows[i + 1]; k++) {
	  const int jj = indcols[k] + 1;
	  	  fprintf(fp, "%d %d %.16e %.16e\n", (i + 1), jj, 
		  coefs_[k].real(), coefs_[k].imag());
	}
      }
    }
    fclose(fp);
  }
  else {
    fprintf(stderr,
	    "%s %d : fail to open %s\n",
	    __FILE__, __LINE__, fname);
    exit(-1);
  }
}

// copy from higher precision than T and U : quadruple <- T = double
template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
CopyQueueFwBw(DissectionSolver<W, Z, T, U> &qdslv)
{
  bool isSym = qdslv.ptDA()->isSymmetric();
  bool isUpper = qdslv.ptDA()->isUpper();
  bool isWhole = qdslv.ptDA()->isWhole();
  _dim = qdslv.dimension();
  _graph_colors = qdslv.graph_colors();
  _scaling = qdslv.scaling();
  _verbose = qdslv.verbose();
  _num_threads = qdslv.num_threads();
  _nsing = qdslv.nsing();
  _fp = qdslv.get_filedescriptor();
  bool verbose = qdslv.verbose();
  //  T *coefs;
  _dissectionQueue = new DissectionQueue<T, U>*[_graph_colors];  
  _tridiagQueue = new TridiagQueue<T, U>*[_graph_colors];
  _dissectionMatrix =
    new vector<DissectionMatrix<T, U>* >[_graph_colors];
  _tridiagMatrix = new TridiagBlockMatrix<T, U>*[_graph_colors];
  
  _Schur = new SchurMatrix<T>[_graph_colors];
  _kernel = new KernelMatrix<T>[_graph_colors];
  _singIdx = new vector<int>[_graph_colors];
  for (int m = 0; m < _graph_colors; m++) {
    _Schur[m].getAcol() = new SparseMatrix<T>(); // dummy allocation
    _Schur[m].getArow() = new SparseMatrix<T>(); // dummy allocation
  }
  DissectionQueue<W, Z>** dQ = qdslv.getDissectionQueue();
  TridiagQueue<W, Z>** tQ = qdslv.getTridiagQueue();
  _precDiag = new U[_dim];
  for (int i = 0; i < _dim; i++) {
    _precDiag[i] = conv_prec<U, Z>(qdslv.addrPrecDiag()[i]);
  }
  
  _ptDA = new SparseMatrix<T>(isSym, isUpper, isWhole);
  
  CopySparseMatrix(_ptDA, qdslv.ptDA());
  //  coefs = _ptDA->getCoef();
  for (int m = 0; m < _graph_colors; m++) {
    if (tQ[m]->tridiagSolver()) {
      const int dim = tQ[m]->dimension();
      _tridiagMatrix[m] = new TridiagBlockMatrix<T, U>(dim,
						       SIZE_B1,
						       isSym,
						       0, // no other Tridiag
						       verbose,
						       _fp);
      CopyTridiagBlockMatrix(*_tridiagMatrix[m],
			     *qdslv.getTridiagBlockMatrix()[m],
			     _ptDA->getCoef());

      _tridiagQueue[m] = new TridiagQueue<T, U>(true, verbose, _fp);
      _tridiagQueue[m]->generate_queue(_tridiagMatrix[m],
				       tQ[m]->dimension(),
				       tQ[m]->nnz(),
				       tQ[m]->isMapped(),
				       tQ[m]->remap_eqn(),
				       tQ[m]->ptRows(),
				       tQ[m]->indCols(),
				       tQ[m]->indVals(),
				       _ptDA->getCoef());
    }
    else { //     if (tQ[m]->tridiagSolver()) 
      Dissection::Tree* btree = qdslv.btree()[m];
      _tridiagQueue[m] = new TridiagQueue<T, U>(false, verbose, _fp);
      const int nb_doms = btree->NumberOfSubdomains();
      _dissectionMatrix[m].resize(nb_doms);
      _dissectionQueue[m] = new DissectionQueue<T, U>(btree,
						      _dissectionMatrix[m],
						      _num_threads,
						      isSym,
						      verbose,
						      _fp);
      
      typename vector<DissectionMatrix<W, Z>* >::const_iterator it = qdslv.getDissectionMatrix()[m].begin();
      typename vector<DissectionMatrix<T, U>* >::const_iterator jt = _dissectionMatrix[m].begin();
      for ( ; it != qdslv.getDissectionMatrix()[m].end(); ++it, ++jt) {
	if ((*it)->islast()) {
	  int color = (*it)->ColorTridiagBlockMatrix();
	  for (int i = 0; i < color; i++) {
	    CopyTridiagBlockMatrix(*(*jt)->addrtridiagBlock()[i],
				   *(*it)->addrtridiagBlock()[i],
				   _ptDA->getCoef());
	  } // loop : i
	}   // if ((*it)->islast())
	else {
	  CopySquareBlockMatrix((*jt)->diagBlock(), (*it)->diagBlock());
	}
	SquareBlockMatrix<T>* diag = (*jt)->addrdiagBlock();
	RectBlockMatrix<T>* lower = (*jt)->addrlowerBlock();
	if (!isSym) {
	  CopyRectBlockMatrix(*lower, (*it)->lowerBlock());
	}
	RectBlockMatrix<T>* upper = (*jt)->addrupperBlock();
	CopyRectBlockMatrix(*upper, (*it)->upperBlock());
	
	CopyDissectionMatrix((*jt),
			     (*it),
			     diag,
			     lower,
			     upper);
      } // loop : it, jt

      _dissectionQueue[m]->generate_queue_fwbw(_dissectionMatrix[m],
					       dQ[m]->dimension(),
					       dQ[m]->nnz(), _ptDA->getCoef());
    } // if (tQ[m]->tridiagSolver())
    CopyKernelMatrix(_kernel[m], qdslv.getKernelMatrix()[m]);
    _singIdx[m] = qdslv.getSingVal()[m];
  } // loop : m
  //  qdslv.GetMatrixScaling(_precDiag);
  _index_isolated = qdslv.getIndexIsolated();
  _status_factorized = true;
  _btree = new Dissection::Tree*[_graph_colors]; // dummy allocation
  _with_btree = false;
}

template
void DissectionSolver<double, double, quadruple, quadruple>::
CopyQueueFwBw(DissectionSolver<quadruple, quadruple, double, double> &qdslv);

template
void DissectionSolver<complex<double>, double, complex<quadruple>, quadruple>::
CopyQueueFwBw(DissectionSolver<complex<quadruple>, quadruple,
	      complex<double>, double> &qdslv);

template
void DissectionSolver<quadruple, quadruple, double, double>::
CopyQueueFwBw(DissectionSolver<double, double, quadruple, quadruple> &qdslv);

template
void DissectionSolver<complex<quadruple>, quadruple, complex<double>, double>::
CopyQueueFwBw(DissectionSolver<complex<double>, double, complex<quadruple>, quadruple> &qdslv);

template
void DissectionSolver<float, float, double, double>::
CopyQueueFwBw(DissectionSolver<double, double, float, float> &qdslv);

template
void DissectionSolver<complex<float>, float, complex<double>, double>::
CopyQueueFwBw(DissectionSolver<complex<double>, double, complex<float>, float> &qdslv);


template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
SymbolicFact(const int dim_,
	     const int *ptRows,
	     const int *indCols,
	     const bool isSym,
	     const bool isUpper,
	     const bool isWhole,
	     const int decomposer_, 
	     const int nbLevels_,
	     const int minNodes)
{
  SymbolicFact_(dim_,
		ptRows[dim_],
		false,
		ptRows,
		indCols,
		(long long int *)NULL, //  const long long int *ptRows,
		(long long int *)NULL, //  const long long int *indCols,
		isSym,
		isUpper,
		isWhole,
		decomposer_, 
		nbLevels_,
		minNodes);
}

template
void DissectionSolver<double>::SymbolicFact(const int dim_,
					    const int *ptRows,
					    const int *indCols,
					    const bool isSym,
					    const bool isUpper,
					    const bool isWhole,
					    const int decomposer,
					    const int nbLevels_,
					    const int minNodes);

template
void DissectionSolver<quadruple>::SymbolicFact(const int dim_,
					       const int *ptRows,
					       const int *indCols,
					       const bool isSym,
					       const bool isUpper,
					       const bool isWhole,
					       const int decomposer, 
					       const int nbLevels_,
					       const int minNodes);

template
void DissectionSolver<complex<double>, double>::
SymbolicFact(const int dim_,
	     const int *ptRows,
	     const int *indCols,
	     const bool isSym,
	     const bool isUpper,
	     const bool isWhole,
	     const int decomposer,
	     const int nbLevels_,
	     const int minNodes);

template
void DissectionSolver<complex<quadruple>, quadruple>::
SymbolicFact(const int dim_,
	     const int *ptRows,
	     const int *indCols,
	     const bool isSym,
	     const bool isUpper,
	     const bool isWhole,
	     const int decomposer,
	     const int nbLevels_,
	     const int minNodes);

template
void DissectionSolver<double, double,
		      quadruple, quadruple>::SymbolicFact(const int dim_,
							  const int *ptRows,
							  const int *indCols,
							  const bool isSym,
							  const bool isUpper,
							  const bool isWhole,
							  const int decomposer,
							  const int nbLevels_,
							  const int minNodes);
template
void DissectionSolver<quadruple, quadruple,
		      double, double>::SymbolicFact(const int dim_,
						    const int *ptRows,
						    const int *indCols,
						    const bool isSym,
						    const bool isUpper,
						    const bool isWhole,
						    const int decomposer, 
						    const int nbLevels_,
						    const int minNodes);

template
void DissectionSolver<complex<double>, double, complex<quadruple>, quadruple>::
SymbolicFact(const int dim_,
	     const int *ptRows,
	     const int *indCols,
	     const bool isSym,
	     const bool isUpper,
	     const bool isWhole,
	     const int decomposer,
	     const int nbLevels_,
	     const int minNodes);

template
void DissectionSolver<complex<quadruple>, quadruple, complex<double>, double>::
SymbolicFact(const int dim_,
	     const int *ptRows,
	     const int *indCols,
	     const bool isSym,
	     const bool isUpper,
	     const bool isWhole,
	     const int decomposer,
	     const int nbLevels_,
	     const int minNodes);

template
void DissectionSolver<float>::SymbolicFact(const int dim_,
					   const int *ptRows,
					   const int *indCols,
					   const bool isSym,
					   const bool isUpper,
					   const bool isWhole,
					   const int decomposer,
					   const int nbLevels_,
					   const int minNodes);

template
void DissectionSolver<complex<float>, float>::SymbolicFact(const int dim_,
							   const int *ptRows,
							   const int *indCols,
							   const bool isSym,
							   const bool isUpper,
							   const bool isWhole,
							   const int decomposer,
							   const int nbLevels_,
							   const int minNodes);
//

template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
SymbolicFact(const int dim_,
	     const long long int *ptRows,
	     const long long int *indCols,
	     const bool isSym,
	     const bool isUpper,
	     const bool isWhole,
	     const int decomposer_, 
	     const int nbLevels_,
	     const int minNodes)
{
  SymbolicFact_(dim_,
		(int)ptRows[dim_],
		true,
		(int *)NULL,         //  const int *ptRows,
		(int *)NULL,         //  const int *indCols,
		ptRows,
		indCols,
		isSym,
		isUpper,
		isWhole,
		decomposer_, 
		nbLevels_,
		minNodes);
}

template
void DissectionSolver<double>::SymbolicFact(const int dim_,
					    const long long int *ptRows,
					    const long long int *indCols,
					    const bool isSym,
					    const bool isUpper,
					    const bool isWhole,
					    const int decomposer,
					    const int nbLevels_,
					    const int minNodes);
template
void DissectionSolver<quadruple>::SymbolicFact(const int dim_,
					       const long long int *ptRows,
					       const long long int *indCols,
					       const bool isSym,
					       const bool isUpper,
					       const bool isWhole,
					       const int decomposer, 
					       const int nbLevels_,
					       const int minNodes);

template
void DissectionSolver<complex<double>, double>::
SymbolicFact(const int dim_,
	     const long long int *ptRows,
	     const long long int *indCols,
	     const bool isSym,
	     const bool isUpper,
	     const bool isWhole,
	     const int decomposer,
	     const int nbLevels_,
	     const int minNodes);

template
void DissectionSolver<complex<quadruple>, quadruple>::
SymbolicFact(const int dim_,
	     const long long int *ptRows,
	     const long long int *indCols,
	     const bool isSym,
	     const bool isUpper,
	     const bool isWhole,
	     const int decomposer,
	     const int nbLevels_,
	     const int minNodes);

template
void DissectionSolver<float>::SymbolicFact(const int dim_,
					   const long long int *ptRows,
					   const long long int *indCols,
					   const bool isSym,
					   const bool isUpper,
					   const bool isWhole,
					   const int decomposer,
					   const int nbLevels_,
					   const int minNodes);

template
void DissectionSolver<complex<float>, float>::SymbolicFact(const int dim_,
							   const long long int *ptRows,
							   const long long int *indCols,
							   const bool isSym,
							   const bool isUpper,
							   const bool isWhole,
							   const int decomposer,
							   const int nbLevels_,
							   const int minNodes);

//
template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
SymbolicFact_(const int dim_,
	      const int nz_,
	      const bool flagint64,
	      const int *ptRows,
	      const int *indCols,
	      const long long int *ptRows64,
	      const long long int *indCols64,
	      const bool isSym,
	      const bool isUpper,
	      const bool isWhole,
	      const int decomposer_, 
	      const int nbLevels_,
	      const int minNodes)
{
  // decomposer = 0: SCOTCH, 1 : METIS, 2 : TRIDIAG(Cuthill-McKee)
  int dim = dim_;
  clock_t t0_cpu, t1_cpu, t2_cpu;
  elapsed_t t0_elapsed, t1_elapsed, t2_elapsed;

  _dim = dim_;
  const int nz = nz_;
  int nbLevels, decomposer;
  bool flag, berr;
  int *map_eqn, *remap_eqn;
  int num_threads;
  
  t0_cpu = clock();
  get_realtime(&t0_elapsed);

  if (flagint64) {
      _ptDA = new SparseMatrix<T>(dim, nz, ptRows64, indCols64,
				  isSym, isUpper, isWhole);
  }
  else {
    _ptDA = new SparseMatrix<T>(dim, nz, ptRows,  indCols,
				isSym, isUpper, isWhole);
  }

  diss_printf(_verbose, _fp,
	      "%s %d : isSym = %s isUpper = %s isWhole = %s decomposer = %d\n",
	      __FILE__, __LINE__,
	      isSym ? "ture" : "false",
	      isUpper ? "ture" : "false",
	      isWhole ? "ture" : "false",
	      decomposer_);
  switch(nbLevels_) {
  case 1 :
    diss_printf(_verbose, _fp, 
		"%s %d : Dissection :: not activated, switch to tridiag\n",
		__FILE__, __LINE__);
    decomposer = TRIDIAG_DECOMPOSER;
    break;
  case (-1):
      nbLevels = (int)log2((double)dim_ / (double)minNodes);
      if (nbLevels < 2) {
        decomposer = TRIDIAG_DECOMPOSER;
      }
      else {
        decomposer = decomposer_;
      }
      break;
  default:
    nbLevels = nbLevels_;
    decomposer = decomposer_;
    break;
  } 

  map_eqn = new int[dim];   // delete [] in DissectionSolver::Destroy(void) 
  remap_eqn = new int[dim]; // delete [] in DissectionSolver::Destroy(void) 
  for (int i = 0; i < dim; i++) {
    map_eqn[i] = i;
    remap_eqn[i] = i;
  }

  CSR_indirect *unsym_csr = new CSR_indirect;
  bool unsym_csr_alloc = false;
  //  int *ptUnsymRows, *indUnsymCols, *indVals, *indSymVals;

  if (isSym) {
    if (isWhole) {
      unsym_csr->n = dim;
      unsym_csr->nnz = nz;
      unsym_csr->ptRows = _ptDA->getRows(); // (int *)ptRows;
      unsym_csr->indCols = _ptDA->getIndCols(); //(int *)indCols;
      unsym_csr->indVals = new int[nz];
      unsym_csr->indVals_unsym = new int[nz]; // for safe use of delete []
    }
    else {
      const int nnz0 = nz * 2 - dim;
      unsym_csr->n = dim;
      unsym_csr->nnz = nnz0;
      unsym_csr->ptRows = new int[dim + 1];
      unsym_csr->indCols = new int[nnz0];
      unsym_csr_alloc = true;
      unsym_csr->indVals = new int[nnz0];
      unsym_csr->indVals_unsym = (int *)NULL; // for safe use of delete []
      // only used as extend symmetric symbolic structure to unsymmetric
      int nnz1;
      nnz1 = CSR_sym2unsym(unsym_csr, 
			   _ptDA->getRows(), _ptDA->getIndCols(), 
			   map_eqn, remap_eqn, //
			   dim, isUpper);
      if (nnz1 != nnz0) {
	fprintf(stderr, 
		"%s %d : symmetric matrix has no diagonal entry %d != %d\n",
		__FILE__, __LINE__, nnz0, nnz1);
	exit(-1);
      }
    }
  }
  else {
    unsym_csr->n = dim;
    unsym_csr->nnz = nz;
    unsym_csr->indVals = new int[nz];
    unsym_csr->indVals_unsym = new int[nz];
    unsym_csr->ptRows = _ptDA->getRows(); //(int *)ptRows;
    unsym_csr->indCols = _ptDA->getIndCols(); //(int *)indCols;
  }
  int* graph_mask = new int[dim];

  {
    _graph_colors = getColorMaskCSR(graph_mask, unsym_csr, _verbose, _fp);
    int count = 0;
    for (int i = 0; i < dim; i++) {
      if (graph_mask[i] == 0) {
	count++;
      }
    }
    _index_isolated.resize(count);
    count = 0;
    for (int i = 0; i < dim; i++) {
      if (graph_mask[i] == 0) {
	_index_isolated[count] = i;
	count++;
      }
    }
  }
  _btree = new Dissection::Tree*[_graph_colors];
  _with_btree = true;
  _dissectionMatrix = new vector<DissectionMatrix<T, U>* >[_graph_colors];
  _dissectionQueue = new DissectionQueue<T, U>*[_graph_colors];
  _tridiagMatrix = new TridiagBlockMatrix<T, U>*[_graph_colors];
  _tridiagQueue = new TridiagQueue<T, U>*[_graph_colors];
  // loop over each connected graph
  for (int m = 0; m < _graph_colors; m++) {
    bool flag_dissection = false;
    int i0, i1, i2;
    const int color = m + 1;
    int count = 0;
    int dim1, dim2;
    dim1 = 0;
    dim2 = 0;
    for (int i = 0; i < dim; i++) {
      if (graph_mask[i] == color) {
	dim1++;
      }
      else if (graph_mask[i] == (-color)) {
	dim2++;
      }
    }
    count = dim1 + dim2;
    i0 = 0;
    i1 = dim1;
    i2 = count;
    for (int i = 0; i < dim; i++) {
      if (graph_mask[i] == color) {
	remap_eqn[i0++] = i;     // remap : new(mapped) to old(original) index
      }
      else if (graph_mask[i] == (-color)) {
	remap_eqn[i1++] = i;     // remap : new(mapped) to old(original) index
      }
      else {
	remap_eqn[i2++] = i;
      }
    }
    for (int i = 0; i < dim; i++) {
      map_eqn[remap_eqn[i]] = i; // map : old(original) to new (mapped) index
    }
    if (isSym) {
      if (isWhole) {
	if (m == 0) {  //reconnect to newly allocated array
	  unsym_csr->ptRows = new int[dim + 1]; 
	  unsym_csr->indCols = new int[nz];
	  unsym_csr_alloc = true;
	}
	CSR_unsym2unsym(unsym_csr,
			_ptDA->getRows(), _ptDA->getIndCols(),
			map_eqn, remap_eqn, //
			dim, _verbose, _fp);
      }
      else {
	CSR_sym2unsym(unsym_csr,
		      _ptDA->getRows(), _ptDA->getIndCols(),
		      map_eqn, remap_eqn, //
		      dim, isUpper, _verbose, _fp);
      }
    }
    else {
      if (m == 0) {  //reconnect to newly allocated array
	unsym_csr->ptRows = new int[dim + 1];
	unsym_csr->indCols = new int[nz];
	unsym_csr_alloc = true;
      }
      CSR_unsym2unsym(unsym_csr,
		      _ptDA->getRows(), _ptDA->getIndCols(),
		      map_eqn, remap_eqn, //
		      dim, _verbose, _fp);
    }

    flag = false;
    
    if (_graph_colors > 1) {
      const int level_tmp = (int)log2((double)count / (double)minNodes);
      nbLevels = nbLevels < level_tmp ? nbLevels : level_tmp;
      nbLevels = nbLevels < 2 ? 2 : nbLevels; //
    }
    diss_printf(_verbose, _fp,
		"%s %d : count = %d decomposer = %d\n", __FILE__, __LINE__,
		count, decomposer);

    if ((count > SIZE_TRIDIAG) && ((decomposer == SCOTCH_DECOMPOSER) ||
				   (decomposer == METIS_DECOMPOSER))) {
      flag_dissection = true;
      diss_printf(_verbose, _fp,
		  "%s %d : %s applied to color %d neq = %d with %d levels\n",
		  __FILE__, __LINE__,
		  (decomposer == METIS_DECOMPOSER) ? "METIS" : "SCOTCH",
		  color, count, 
		  nbLevels);
      while (flag == false) { 
	_btree[m] = 
	  new Dissection::Tree(_fp, 
			       berr, 
			       (unsigned)count,
			       unsym_csr,
			       isSym,
			       remap_eqn,
			       nbLevels, 
			       minNodes,
#ifdef NO_METIS
			       NULL,
#else
			       ((decomposer == METIS_DECOMPOSER) ? MetisSplitter : NULL),
#endif
			       true, _verbose);
	diss_printf(_verbose, _fp, 
		    "%s %d : Dissection::Tree %s decomposition\n",
		    __FILE__, __LINE__,
		    (decomposer == METIS_DECOMPOSER) ? "METIS" : "SCOTCH");
	flag = true;
	for (int d = 1; d <= _btree[m]->NumberOfSubdomains(); d++) {
	  //	  if (_btree[m]->sizeOfDomain(d) == 0) {
	  const int isd = _btree[m]->sizeOfDomain(d);
	  if (isd <= DIM_AUG_KERN) {
	    if (d > 1) {
	      diss_printf(_verbose, _fp, 
			  "%s %d : Dissection:: %d-th too small nrow = %d\n",
			  __FILE__, __LINE__, d, isd);
	    }
	    else {
	      fprintf(stderr,
		      "%s %d : too small first bisector %d not yet accepted\n",
		      __FILE__, __LINE__, isd);
	      exit(-1);
	    }
	  } // if (isd == 0)
	} // loop : d
	if (flag) {
	  break;
	}
	else {
	  if (berr == true) {
	    delete _btree[m];  // for retry of graph partitioning
	  }
	}
      } // while
      //      _nbLevels = _btree[m]->NumberOfLevels();
      diss_printf(_verbose, _fp, 
		  "%s %d : Tree strategy = 0 : nblevels = %d -> %d\n",
		__FILE__, __LINE__,
		nbLevels, _btree[m]->NumberOfLevels());

      if ((_btree[m]->NumberOfLevels() == 1) 
	  || ((nbLevels - _btree[m]->NumberOfLevels()) > 3) || berr == false) {
#ifdef NO_METIS
	flag_dissection = false;
#else
	flag_dissection = true;
	diss_printf(_verbose, _fp, 
		    "%s %d :: retry partitioning by METIS : %d\n",
		    __FILE__, __LINE__, nbLevels);
	nbLevels = nbLevels > 4 ? (nbLevels - 3) : 2;
	flag = false;
	while (flag == false) { 
	  _btree[m] = 
	    new Dissection::Tree(_fp, berr, (unsigned)count,
				 unsym_csr,
				 isSym,
				 remap_eqn,
				 nbLevels, minNodes, 
				 MetisSplitter,
				 true, _verbose);
	  diss_printf(_verbose, _fp, 
		      "%s %d :: Tree METIS decomposition **\n",
		      __FILE__, __LINE__ );
	  flag = true;
	  for (int d = 1; d <= _btree[m]->NumberOfSubdomains(); d++) {
	    if (_btree[m]->sizeOfDomain(d) == 0) {
	      flag = false;
	    }
	  }
	  if (flag) {
	    break;
	  }
	  else {
	    delete _btree[m];
	    nbLevels--;
	  }
	} // while (flag == flase)
#endif
      } // if
    } // if ((count > SIZE_TRIDIAG) && ((decomposer == SCOTCH_DECOMPOSER)...

    if (flag_dissection) {
       bool verify_root = (_btree[m]->sizeOfDomain(1) <= SIZE_B1);
      int nb_doms = _btree[m]->NumberOfSubdomains();
      _dissectionMatrix[m].resize(nb_doms);
      diss_printf(_verbose, _fp,
		"%s %d : ", __FILE__, __LINE__);
      diss_printf(_verbose, _fp, 
		  "num_threads = %d dim = %d dim_diss. = %d nbLevels = %d\n", 
		  _num_threads, _dim, count, _btree[m]->NumberOfLevels());
      for (int d = 1; d <= _btree[m]->NumberOfSubdomains(); d++) {
	diss_printf(_verbose, _fp,
		    "( %d % d) ", d, _btree[m]->sizeOfDomain(d));
      }
      diss_printf(_verbose, _fp, "\n");
      const int level_last = _btree[m]->NumberOfLevels();
      const int nb_doms_dense0 = (1U << (level_last - 1));
      if (_num_threads > nb_doms_dense0 || verify_root) {
	diss_printf(_verbose, _fp,
		     "%s %d : ", __FILE__, __LINE__);
	diss_printf(_verbose, _fp, 
		    "thread number reduced : %d -> %d\n",
		    _num_threads, 1);
	diss_printf(_verbose, _fp,
		    "which is used as number of threads for computation\n");
	//	num_threads = nb_doms_dense0;
	num_threads = 1;
      }
      else {
	num_threads = _num_threads;
      }
      _dissectionQueue[m] = new DissectionQueue<T, U>(_btree[m], 
						      _dissectionMatrix[m],
						      num_threads,
						      isSym,
						      _verbose,
						      _fp);
      // generation of queue for numerical factorization is treated as a part of
      // symbolic factorization : _queue_symb, {_queue_static, _queue_dynamic}
      _dissectionQueue[m]->generate_queue(_dissectionMatrix[m], 
					  _ptDA->nnz(), _ptDA->getCoef());
      _dissectionQueue[m]->generate_queue_fwbw(_dissectionMatrix[m], 
					       _ptDA->dimension(),
					       _ptDA->nnz(), _ptDA->getCoef());
      _tridiagQueue[m] = new TridiagQueue<T, U>(false, _verbose, _fp);
    }  //   if (flag_dissection) 
    else {
	diss_printf(_verbose, _fp,
		    "%s %d : TridiagSolver : m = %d count = %d\n",
		    __FILE__, __LINE__, m, count);
      const int nnz = unsym_csr->ptRows[count];
      bool isMapped = ((_graph_colors > 1) || (_index_isolated.size() > 0));
      _tridiagMatrix[m] = new TridiagBlockMatrix<T, U>(count, SIZE_B1, isSym,
						       0, // no other Tridiag
						       _verbose, _fp);
      _tridiagQueue[m] = new TridiagQueue<T, U>(true, _verbose, _fp);
      _tridiagQueue[m]->generate_queue(_tridiagMatrix[m], count, nnz,
				       isMapped, remap_eqn, 
				       unsym_csr->ptRows,
				       unsym_csr->indCols,
				       unsym_csr->indVals,
				       _ptDA->getCoef());
    }
  }  // loop : m      
  if (unsym_csr_alloc) {
    delete [] unsym_csr->ptRows;
    delete [] unsym_csr->indCols;
  }
  delete [] unsym_csr->indVals;
  delete [] unsym_csr->indVals_unsym;   // isSym == (int *)NULL
  delete unsym_csr;
  delete [] map_eqn;
  delete [] remap_eqn;

  delete [] graph_mask;      // added 01 Oct.2013 Atsushi
  t1_cpu = clock();
  get_realtime(&t1_elapsed);

  for (int m = 0; m < _graph_colors; m++) {
    if (_tridiagQueue[m]->tridiagSolver()) {
      _tridiagQueue[m]->exec_symb_fact();
    }
    else {
      _dissectionQueue[m]->exec_symb_fact();
    }
  }

  _precDiag = new U[_dim]; // preparation for numerical fact.

  _Schur = new SchurMatrix<T>[_graph_colors];
  _kernel = new KernelMatrix<T>[_graph_colors];
  _singIdx = new vector<int>[_graph_colors];
  _status_factorized = false;
  t2_cpu = clock();
  get_realtime(&t2_elapsed);
  diss_printf(_verbose, _fp,
		"graph paritioner. : cpu time = %.4e elapsed time = %.4e\n", 
	      (double)(t1_cpu - t0_cpu) / (double)CLOCKS_PER_SEC,
	      convert_time(t1_elapsed, t0_elapsed));

  diss_printf(_verbose, _fp,
	      "queue symb. fact. : cpu time = %.4e elapsed time = %.4e\n", 
	      (double)(t2_cpu - t1_cpu) / (double)CLOCKS_PER_SEC,
	      convert_time(t2_elapsed, t1_elapsed));
}


template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
NumericFact(const int called,
	    T *coefs,
	    const int scaling,
	    const double eps_pivot,
	    const bool kernel_detection_all,
	    const int dim_augkern,
	    const double machine_eps_,
	    const bool assume_invertible,
	    const bool higher_precision)
{
  clock_t t0_cpu, t1_cpu, t2_cpu, t3_cpu;
  elapsed_t t0_elapsed, t1_elapsed, t2_elapsed, t3_elapsed;
  const U eps_machine = machine_eps_ < 0.0 ? machine_epsilon<U, Z>() : U(machine_eps_);
  t0_cpu = clock();
  _assume_invertible = assume_invertible;
  _dim_augkern = dim_augkern;
  get_realtime(&t0_elapsed);

  diss_printf(_verbose, _fp, "%s %d : eps_machine = %s scaling = %d\n",
	      __FILE__, __LINE__,
	      tostring<U>(eps_machine).c_str(), scaling);
  
  _scaling = scaling;

  //  _ptDA->normalize(_scaling, coefs, _precDiag);
  normalize<T, U>(_scaling, coefs, _ptDA, _precDiag);
  t1_cpu = clock();
  get_realtime(&t1_elapsed);
  for (int m = 0; m < _graph_colors; m++) {
    if (_tridiagQueue[m]->tridiagSolver()) {
      _tridiagQueue[m]->exec_num_fact(called, 
				      eps_pivot, 
				      true, // matrix will not be decomposed
				      dim_augkern,
				      eps_machine,
				      higher_precision);
    }
    else {
      _dissectionQueue[m]->exec_num_fact(called, 
					 eps_pivot, 
					 kernel_detection_all, 
					 dim_augkern,
					 eps_machine,
					 higher_precision);
    }
  }
  t2_cpu = clock();
  get_realtime(&t2_elapsed);

  // dense {0 [1 2] [3 4 5 6] .... } + sparse {[2^(_nbLevels-1)-1... ]}
  int sz_total = 0;
  VectorArray<T> work(_dim);
  for (int m = 0; m < _graph_colors; m++) {
    if (_tridiagQueue[m]->tridiagSolver()) {
      const int nrow = _tridiagMatrix[m]->nrow();
      const int nn0 = nrow - _tridiagMatrix[m]->rank();
      if (nn0 > 0) {
	_tridiagMatrix[m]->SingularNode(_kernel[m].getKernListEq());
	//					_ptDA_kern_list_eq[m]);
	_kernel[m].set_dimension(nn0);
	_kernel[m].getKernProj().init(nn0);
	//_ptDA_kern_proj[m].init(nn0); 
	_tridiagMatrix[m]->KernelBasis(false, _kernel[m].getKernBasis());
				       //_ptDA_kern_basis[m]); // isTrans
	if (_tridiagQueue[m]->isMapped()) {
	  int *remap_eqn = _tridiagQueue[m]->remap_eqn();
	  ColumnMatrix<T> kernel_tmp;
	  kernel_tmp.init(_dim, nn0);
	  for (int n = 0; n < nn0; n++) {
	    for (int i = 0; i < _dim; i++) {
	      kernel_tmp(i, n) = _zero;
	    }
	    for (int i = 0; i < nrow; i++) {
	      kernel_tmp(remap_eqn[i], n) = _kernel[m].getKernBasis()(i, n);
	    }
	  }
	  _kernel[m].getKernBasis().free();
	  _kernel[m].getKernBasis().init(_dim, nn0);
	  for (int n = 0; n < nn0; n++) {
	    for (int i = 0; i < _dim; i++) {
	      _kernel[m].getKernBasis()(i, n) =  kernel_tmp(i, n);
	    }
	  }
	  kernel_tmp.free();
	} // 	if (_tridiagQueue[m]->isMapped()) 
	T *kern_basis = _kernel[m].getKernBasis().addrCoefs();
	              //_ptDA_kern_basis[m].addrCoefs();
	diss_printf(_verbose, _fp,
		    "%s %d : residual of kernel vertors : %d\n",
		    __FILE__, __LINE__, nn0);
      for (int i = 0; i < nn0; i++) {
	_ptDA->prod((kern_basis + (i * _dim)), work.addrCoefs());
	U stmp = blas_l2norm<T, U>(_dim, work.addrCoefs(), 1);
	diss_printf(_verbose, _fp, "%d %s\n", i, tostring<U>(stmp).c_str());
      }
      if(_scaling) {
	for (int j = 0; j < nn0; j++) {
	  const int jtmp = j * _dim;
	  for (int i = 0; i < _dim; i++) {
	    kern_basis[i + jtmp] *= _precDiag[i];
	  }
	}
      }
      _kernel[m].getTKernBasis().init(_dim, nn0);
      
      // normalize each kernel_basis
      for (int j = 0; j < nn0; j++) {
	U stmp(0.0);
	U one(1.0);
	stmp = blas_l2norm<T, U>(_dim, kern_basis + (j * _dim), 1);
	stmp = one / stmp;
	blas_scal2<T, U>(_dim, stmp, (kern_basis + (j * _dim)), 1);
      }
      
      for (int i = 0; i < nn0; i++) {
	for(int j = i; j < nn0; j++) {
	  //	    _ptDA_kern_proj[m](i, j) =
	  _kernel[m].getKernProj()(i,j) =
	    blas_dot<T>(_dim,
			(kern_basis + (i * _dim)), 1,
			(kern_basis + (j * _dim)), 1);
	}
	// symmetrize
	  for (int j = (i + 1); j < nn0; j++) {
	    _kernel[m].getKernProj()(j,i) = _kernel[m].getKernProj()(i,j);
	    //	    _ptDA_kern_proj[m](j, i) = _ptDA_kern_proj[m](i, j);
	  }
	}
	full_ldlh<T>(nn0, _kernel[m].getKernProj().addrCoefs(), nn0);
	if (!_tridiagMatrix[m]->isSym()) {  
	  T *kernt_basis = _kernel[m].getTKernBasis().addrCoefs();
	  _tridiagMatrix[m]->KernelBasis(true, _kernel[m].getTKernBasis());
	  if (_tridiagQueue[m]->isMapped()) {
	    int *remap_eqn = _tridiagQueue[m]->remap_eqn();
	    ColumnMatrix<T> kernel_tmp;
	    kernel_tmp.init(_dim, nn0);
	    for (int n = 0; n < nn0; n++) {
	      for (int i = 0; i < _dim; i++) {
		kernel_tmp(i, n) = _zero;
	      }
	      for (int i = 0; i < nrow; i++) {
		kernel_tmp(remap_eqn[i], n) = _kernel[m].getTKernBasis()(i, n);
	      }
	    }
	    _kernel[m].getTKernBasis().free();
	    _kernel[m].getTKernBasis().init(_dim, nn0);
	    for (int n = 0; n < nn0; n++) {
	      for (int i = 0; i < _dim; i++) {
		_kernel[m].getTKernBasis()(i, n) =  kernel_tmp(i, n);
	      }
	    }
	    kernel_tmp.free();
	  } // 	if (_tridiagQueue[m]->isMapped()) 
	  diss_printf(_verbose, _fp,
		      "%s %d : residual of transposed kernel vertors : %d\n",
		      __FILE__, __LINE__, nn0);
	  for (int i = 0; i < nn0; i++) {
	    _ptDA->prodt((kernt_basis + (i * _dim)), work.addrCoefs());
	    U stmp = blas_l2norm<T, U>(_dim, work.addrCoefs(), 1);
	    diss_printf(_verbose, _fp, "%d %s\n", i, tostring<U>(stmp).c_str());
	  }
	  if(_scaling) {
	    for (int j = 0; j < nn0; j++) {
	      const int jtmp = j * _dim;
	      for (int i = 0; i < _dim; i++) {
		kernt_basis[i + jtmp] *= _precDiag[i];
	      }
	    }
	  }
	} //  if (!_tridiagMatrix[m]->isSym())
      } // if (nn0 > 0)
      else {
	_kernel[m].set_dimension(0);
	_kernel[m].getKernProj().init(0);
	//	_ptDA_kern_proj[m].init(0); // = new SquareMatrix<double>();
	_kernel[m].getKernBasis().init(0, 0); //  = new ColumnMatrix<double>();
	_kernel[m].getKernListEq().clear();
	_kernel[m].getKernListEqLeft().clear();
      }
      _singIdx[m].resize(0);
      _Schur[m].getAcol() = new SparseMatrix<T>(); // dummy allocation
      _Schur[m].getArow() = new SparseMatrix<T>(); // dummy allocation
      //   if (_status_factorized) {  // for safty of second call
	_Schur[m].getSlduList().clear(); //
	_Schur[m].getSlduListLeft().clear(); //
	_Schur[m].getSldu().free();
    } //  if (_tridiagQueue[m]->tridiagSolver())
    else {
      int sz = 0, sz1 = 0;
      bool kernel_detected = true;
      vector<int> kernel_found;
      //      int id_dom = 1;
      for (int j = 0; j < _dissectionMatrix[m].size(); j++ ) {
	vector<int>& singIdx = _dissectionMatrix[m][j]->singIdxPermute();
	if (singIdx.size() > 0) { 
	  diss_printf(_verbose, _fp,
		      "%s %d : ", __FILE__, __LINE__);
	  diss_printf(_verbose, _fp,
		      "%d : %s : %d :: %d [ ", (j + 1),  // selfIndex() == j
		      _dissectionMatrix[m][j]->KernelDetected() ? "true" : "false",
		      _dissectionMatrix[m][j]->nrow(),
		      (int)singIdx.size());
	  vector<int>::const_iterator it =  singIdx.begin();
	  for ( ;it !=  singIdx.end(); ++it) {   
	    diss_printf(_verbose, _fp, "%d ", *it);
	  }
	  diss_printf(_verbose, _fp, "]\n");
	} // if (singIdx.size() > 0) 
	if (singIdx.size() > 0) {
	  //	  id_dom = j + 1;
	  kernel_found.push_back(j);
	  sz += singIdx.size(); 
	  if (!_dissectionMatrix[m][j]->KernelDetected()) {
	    //if (true) { // 29 Oct.2014 debug
	    sz1 += singIdx.size();
	    diss_printf(_verbose, _fp,
		    "%s %d : ", __FILE__, __LINE__);
	    diss_printf(_verbose, _fp, "domain %d kernel is not detected\n", j);
	    kernel_detected = false;
	  }
	}
      } // loop : j	
      sz_total += sz;
      if (sz > 0) {
	if (kernel_detected && (!_assume_invertible)) {
	  _singIdx[m].resize(sz); 
	  int k = 0;
	  for (int j = 0; j < _dissectionMatrix[m].size(); j++ ) {
	    const int *loc2glob = _btree[m]->getDiagLoc2Glob(j + 1);
	    for (vector<int>::const_iterator it = 
		   _dissectionMatrix[m][j]->singIdxPermute().begin();
		 it !=  _dissectionMatrix[m][j]->singIdxPermute().end();
		 ++it, k++) {
	      _singIdx[m][k] = loc2glob[*it];
	    } // loop : d
	  }
#if 0
	  if (_verbose)  {
	    std::sort(_singIdx[m].begin(), _singIdx[m].end());
	    char fname[256];
	    FILE *fp;
	    int pid = get_process_id();
	    sprintf(fname, "singularindex.%06d.data", pid);
	    if ((fp = fopen(fname, "w")) != NULL) {
	      fprintf(fp, "# %d \n",  (int)_singIdx[m].size());
	      for (int i = 0; i < _singIdx[m].size(); i++) {
		fprintf(fp, "%d\n", _singIdx[m][i]);
	      }
	      fclose(fp);
	    }
	    else {
	      fprintf(stderr,
		      "%s %d : fail to open %s\n",
		      __FILE__, __LINE__, fname);
	      exit(-1);
	    }
	  } // if (_verbose)

	  BuildKernels(_singIdx[m], 
		       sz,
		       _Schur[m], _kernel[m]);
#else
	  vector<SquareBlockMatrix<T>* > diags; // dummy
	  vector<vector<int> >augkern_indexes;  // dummy

	  BuildKernelsDetection(sz, 
				_singIdx[m],
				augkern_indexes,
				diags,
				eps_pivot,
				dim_augkern,
				_Schur[m],
				_kernel[m], false);
#endif
	  // need to copy kern basis projection matrix from global to m-th array
	} // if (kernel_detected)
	else {
	  diss_printf(_verbose, _fp, 
		      "%s %d : candidates of null pivots = %d, to be fixed\n",
		      __FILE__, __LINE__, sz);
	  _singIdx[m].resize(sz + dim_augkern); //= vector<int>
	  int k = 0;
	  for (int d = 1; d <= _btree[m]->NumberOfSubdomains(); d++) {
	    const int j = _btree[m]->selfIndex(d);
	    if (_dissectionMatrix[m][j]->KernelDetected() &&
		(_dissectionMatrix[m][j]->diagBlock().dim_kern() > 0)) {
	      _dissectionMatrix[m][j]->diagBlock().set_KernelDetected(false);
	    }
	    const int *loc2glob = _btree[m]->getDiagLoc2Glob(d);
	    for (vector<int>::const_iterator it = 
		   _dissectionMatrix[m][j]->singIdxPermute().begin();
		 it !=  _dissectionMatrix[m][j]->singIdxPermute().end(); 
		 ++it, k++) {
	      _singIdx[m][k] = loc2glob[*it];
	    }
	  } // loop : d
      // selecting dim_augkern indices for comparison in the kernel detection
	  diss_printf(_verbose, _fp,
		      "%s %d : size = %d : ", __FILE__, __LINE__,
		      (int)_singIdx[m].size());
	  for (int i = 0; i < _singIdx[m].size(); i++) {
	    diss_printf(_verbose, _fp, "%d ", _singIdx[m][i]);
	  }
	  diss_printf(_verbose, _fp, "\n");
	  vector<SquareBlockMatrix<T>* > diags;
	  vector<vector<int> >augkern_indexes;
	  int nn = 1; // level of dense part with 1 indexing
	  k = 0;
	  while (nn <= dim_augkern) {
	    SquareBlockMatrix<T>* diag =
		   &_dissectionMatrix[m][_btree[m]->selfIndex(nn)]->diagBlock();
	    const int *loc2glob = _btree[m]->getDiagLoc2Glob(nn);
	    vector<int> &permute = diag->getPermute();
	    vector<int> augkern_index;
	    //	    augkern_index.resize(dim_augkern);
	    int n = diag->dimension() - 1;  //diag.dimension() == permute.size()
	    int k0 = 0;
	    while ((k < dim_augkern) && (n >= 0)) {
	      const int itmp = loc2glob[permute[n]];
	      bool flag = true;
	      for (int i = 0; i < sz; i++) {
		if(_singIdx[m][i] == itmp) {
		  flag = false;
		  break;
		}
	      }
	      if (flag) {
		_singIdx[m][sz + k0] = itmp;
		augkern_index.push_back(n);
		k0++;
		k++;
	      }
	      n--;
	    } // while
	    diags.push_back(diag);
	    augkern_indexes.push_back(augkern_index);
	    if (k == dim_augkern) {
	      break;
	    }
	    nn++; // 
	  } // while (nn <= dim_augkern)
	  diss_printf(_verbose, _fp,
		      "%s %d : adding local node using %d diag matrices ",
		      __FILE__, __LINE__, (int)augkern_indexes.size());
	  for (vector<vector<int> >::const_iterator it = augkern_indexes.begin();
	       it != augkern_indexes.end(); ++it) {
	    for (vector<int>::const_iterator jt = (*it).begin();
		 jt != (*it).end();
		 ++jt) {
	      diss_printf(_verbose, _fp, "%d ", (*jt));
	    }
	  }
	  diss_printf(_verbose, _fp, "to the last level.\n");
	  diss_printf(_verbose, _fp,
		      "%s %d : singIdx = %d : ",
		      __FILE__, __LINE__, (int)_singIdx[m].size());
	  for (int i = 0; i < _singIdx[m].size(); i++) {
	    diss_printf(_verbose, _fp, "%d ", _singIdx[m][i]);
	  }
	  diss_printf(_verbose, _fp, "\n");
	  std::sort(_singIdx[m].begin(), _singIdx[m].end());
	  diss_printf(_verbose, _fp,
		      "%s %d : singIdx sorted = %d : ",
		      __FILE__, __LINE__, (int)_singIdx[m].size());
	  for (int i = 0; i < _singIdx[m].size(); i++) {
	    diss_printf(_verbose,_fp, "%d ", _singIdx[m][i]);
	  }
	  diss_printf(_verbose, _fp, "\n");
	  int kern_dim = 0;
	  BuildKernelsDetection(kern_dim, 
				_singIdx[m],
				augkern_indexes, 
				diags,
				eps_pivot,
				dim_augkern,
				_Schur[m],
				_kernel[m],
				true);
	  if (kern_dim == (-1)) {
	    int pid = get_process_id();
	    SaveCSRMatrix(_called, coefs);
	    fprintf(stderr, 
		    "%s %d : DissectionSolver CSR matrix dumped : %d %d\n", 
		    __FILE__, __LINE__, _called, pid);
	    exit(-1);
	  }
	}  // if (kernel_detected)
      } // if (sz > 0) 
      else {
	_singIdx[m].resize(0);
	_Schur[m].getAcol() = new SparseMatrix<T>(); // dummy allocation
	_Schur[m].getArow() = new SparseMatrix<T>(); // dummy allocation
	  //	_ptDA_arow[m] = new SparseMatrix<T, W, Z>(true);
	_kernel[m].set_dimension(0);
	_kernel[m].getKernProj().init(0);   //  _ptDA_kern_proj[m].init(0); //
	_kernel[m].getKernBasis().init(0, 0); //
	_kernel[m].getKernListEq().clear();
	_kernel[m].getKernListEqLeft().clear();
	//	if (_status_factorized) {          // for safty of second call
	  _Schur[m].getSlduList().clear(); // _ptDA_sldu_list[m].clear();
	  _Schur[m].getSlduListLeft().clear(); // _ptDA_sldu_list[m].clear(); 
	  _Schur[m].getSldu().free();      // 	  _ptDA_sldu[m].free();
      }
    } //  if (_tridiagQueue[m]->tridiagSolver()) 
  }   // loop : _graph_colors
  work.free();
  {
    int itmp = 0;
    for (int m = 0; m < _graph_colors; m++) {
      itmp += _singIdx[m].size();
    }
    _nsing = itmp; // keep total dimension of kernel of singluar blocks
  }
  t3_cpu = clock();
  get_realtime(&t3_elapsed);

  //clock_gettime(CLOCK_REALTIME, &ts3);
  diss_printf(_verbose, _fp,
	      "queue num. fact. : cpu time = %.4e elapsed time = %.4e\n", 
	      (double)(t2_cpu - t1_cpu) / (double)CLOCKS_PER_SEC,
	      convert_time(t2_elapsed, t1_elapsed));
  diss_printf(_verbose, _fp,
	      "total num. fact. : cpu time = %.4e elapsed time = %.4e\n", 
	      (double)(t3_cpu - t0_cpu) / (double)CLOCKS_PER_SEC,
	      convert_time(t2_elapsed, t0_elapsed));
  diss_printf(_verbose, _fp,
	      "scaling matrix   : cpu time = %.4e elapsed time = %.4e\n", 
	      (double)(t1_cpu - t0_cpu) / (double)CLOCKS_PER_SEC,
	      convert_time(t1_elapsed, t0_elapsed));
  if (sz_total > 0) {
    diss_printf(_verbose, _fp,
		"kernel vec. gen  : cpu time = %.4e elapsed time = %.4e\n", 
	      (double)(t3_cpu - t2_cpu) / (double)CLOCKS_PER_SEC,
	      convert_time(t3_elapsed, t2_elapsed));
  } // (sz_total > 0) {
  // check error
  // this should be a function with specialized template
  if (_verbose) {
    int count = 0;
    for (int m = 0; m < _graph_colors; m++) {
      diss_printf(_verbose, _fp,
	      "%s %d :negative diagnol entries : color = %d\n",
	      __FILE__, __LINE__, m);
      if (_tridiagQueue[m]->tridiagSolver()) {
	int count0;
	count0 = _tridiagMatrix[m]->NumNegativeDiags();
	if (count0 > 0) {
	  diss_printf(_verbose, _fp,
		      "%d / %d\n", count0, _tridiagMatrix[m]->nrow());
	}
	count += count0;
      }
      else {
	int level;
	level = (_btree[m]->NumberOfLevels() - 1);
	const unsigned begdom = 1U << level;
	const unsigned enddom = begdom * 2;
	
	for (int d = begdom; d < enddom; d++) {   
	  const int d0 = _btree[m]->selfIndex(d);
	  //	const void *diagsparse = _dissectionMatrix[m][d0]->diagSparse();
	  TridiagBlockMatrix<T, U> **tridiag = _dissectionMatrix[m][d0]->addrtridiagBlock();
	  int colors = _dissectionMatrix[m][d0]->ColorTridiagBlockMatrix();
	  int count0 = 0;
	  for (int i = 0; i < colors; i++) {
	    count0 += tridiag[i]->NumNegativeDiags();
	  }
	  if (count0 > 0) {
	    diss_printf(_verbose, _fp, "%d : %d / %d\n",
		    d0, count0, _dissectionMatrix[m][d0]->nrow());
	  }
	  count += count0;
	}
	for (int level = (_btree[m]->NumberOfLevels() - 2); level >= 0; level--) {
	  const unsigned begdom = 1U << level;
	  const unsigned enddom = begdom * 2;
	  
	  for (int d = begdom; d < enddom; d++) {   
	    const int d0 = _btree[m]->selfIndex(d);
	    int count0;
	    SquareBlockMatrix<T>& diag = _dissectionMatrix[m][d0]->diagBlock();
	    count0 = count_diag_negative<T>(diag);
	    if (count0 > 0) {
	      diss_printf(_verbose, _fp, "%d : %d / %d\n",
		      d0, count0, _dissectionMatrix[m][d0]->nrow());
	    }
	    count += count0;
	  }
	}
      } // (_tridiagQueue[m]->tridiagSolver()) {
      {
	int count0;
	SubSquareMatrix<T>& diag = _Schur[m].getSldu();
	count0 = count_diag_negative<T>(diag);
	if (count0 > 0) {
	  diss_printf(_verbose, _fp,
		      "-1 : %d / %d\n", count0, diag.dimension());
	}
	count += count0;
      }
    } // loop : m
    if (count > 0) {
      diss_printf(_verbose, _fp,
		  "%s %d :negative diagonal entries = %d / %d\n", 
		  __FILE__, __LINE__, count, _dim);
    }
  } // if (_verbose)
    //
  bool kernel_flag = false;
  for (int m = 0; m < _graph_colors; m++) {
    if (_kernel[m].dimension() != 0) {
      kernel_flag = true;
      break;
    }
  }
#if 0
  if (kernel_flag && (!_ptDA->isSymmetric())) {
    ComputeTransposedKernels();
  }
#endif
  if(_verbose) {
    int n1 = 0;
    for (int m = 0; m < _graph_colors; m++) {
      n1 += _Schur[m].getSlduList().size();
    }
    VectorArray<T> xx(_dim);
    VectorArray<T> yy(_dim);
    VectorArray<T> ww(_dim);
    for (int i = 0; i < _dim; i++) {
      ww[i] = T(2.0 * ((double)rand() / (double)RAND_MAX) - 1.0);
//    ww[i] = T((double)(i % 11));
    }
    VectorArray<T> bb(_dim);
    bb.ZeroClear();
    diss_printf(_verbose, _fp, "%s %d : compute error in (Ker A)^\\perp\n",
		__FILE__, __LINE__);
    SpMV(ww.addrCoefs(), xx.addrCoefs(), true); //  isScaling = true
    if (kernel_flag) {
      ProjectionKernelOrthSingle(xx.addrCoefs(), "creating solution", false);
    }
    SpMV(xx.addrCoefs(), bb.addrCoefs(), true);  //  isScaling = true
    if (kernel_flag) {
      ProjectionKernelOrthSingle(bb.addrCoefs(), "RHS", true);
      //      ProjectionImageSingle(bb.addrCoefs(), "RHS", false);
    }
    blas_copy<T>(_dim, bb.addrCoefs(), 1, yy.addrCoefs(), 1);
    SolveSingle(yy.addrCoefs(), false, false, true); // isScaling = true
    if (kernel_flag) {
      ProjectionKernelOrthSingle(yy.addrCoefs(), "computed solution", false);
    }
    SpMV(yy.addrCoefs(), ww.addrCoefs(), true);
    U norm0, norm1;
    norm0 = blas_l2norm<T, U>(_dim, xx.addrCoefs(), 1); 
    blas_axpy<T>(_dim, _none, xx.addrCoefs(), 1, yy.addrCoefs(), 1);
    norm1 = blas_l2norm<T, U>(_dim, yy.addrCoefs(), 1);
    diss_printf(_verbose, _fp, "%s %d : error = %s / %s = %s\n",
	    __FILE__, __LINE__,
	    tostring<U>(norm1).c_str(), tostring<U>(norm0).c_str(),
	    tostring<U>(norm1 / norm0).c_str());
    if (todouble<U>(norm1 / norm0) > 1.0e-5) {
      int pid = get_process_id();
      int nnz = _ptDA->nnz();
      diss_printf(_verbose, stderr, "%s %d : too large error = %s / %s = %s\n",
		  __FILE__, __LINE__,
		  tostring<U>(norm1).c_str(), tostring<U>(norm0).c_str(),
		  tostring<U>(norm1 / norm0).c_str());
      //      SaveMMMatrix(_called, coefs);
      //      _status_factorized = false;
      //      return;
    }
    norm0 = blas_l2norm<T, U>(_dim, bb.addrCoefs(), 1); 
    blas_axpy<T>(_dim, _none, bb.addrCoefs(), 1, ww.addrCoefs(), 1);
    norm1 = blas_l2norm<T, U>(_dim, ww.addrCoefs(), 1);
    diss_printf(_verbose, _fp, "%s %d : residual = %s / %s = %s\n",
	    __FILE__, __LINE__,
	    tostring<U>(norm1).c_str(), tostring<U>(norm0).c_str(),
	    tostring<U>(norm1 / norm0).c_str());

  } //  if (_verbose)
  _status_factorized = true;
}

template
void DissectionSolver<double>::NumericFact(const int called,
					   double *coefs,
					   const int scaling,
					   const double eps_pivot,
					   const bool kernel_detection_all,
					   const int dim_augkern,
					   const double machine_eps_,
					   const bool assume_invertible,
					   const bool higher_precision);

template
void DissectionSolver<complex<double>, double>::
NumericFact(const int called,
	    complex<double> *coefs,
	    const int scaling,
	    const double eps_pivot,
	    const bool kernel_detection_all,
	    const int dim_augkern,
	    const double machine_eps_,
	    const bool assume_invertible,
	    const bool higher_precision);

template
void DissectionSolver<quadruple>::NumericFact(const int called,
					      quadruple *coefs,
					      const int scaling,
					      const double eps_pivot,
					      const bool kernel_detection_all,
					      const int dim_augkern,
					      const double machine_eps_,
					      const bool assume_invertible,
					      const bool higher_precision);
template
void DissectionSolver<complex<quadruple>, quadruple>::
NumericFact(const int called,
	    complex<quadruple> *coefs,
	    const int scaling,
	    const double eps_pivot,
	    const bool kernel_detection_all,
	    const int dim_augkern,
	    const double machine_eps_,
	    const bool assume_invertible,
	    const bool higher_precision);

template
void DissectionSolver<double, double, quadruple, quadruple>::
NumericFact(const int called,
	    double *coefs,
	    const int scaling,
	    const double eps_pivot,
	    const bool kernel_detection_all,
	    const int dim_augkern,
	    const double machine_eps_,
	    const bool assume_invertible,
	    const bool higher_precision);

template
void DissectionSolver<complex<double>, double, complex<quadruple>, quadruple>::
NumericFact(const int called,
	    complex<double> *coefs,
	    const int scaling,
	    const double eps_pivot,
	    const bool kernel_detection_all,
	    const int dim_augkern,
	    const double machine_eps_,
	    const bool assume_invertible,
	    const bool higher_precision);

template
void DissectionSolver<quadruple, quadruple, double, double>::
NumericFact(const int called,
	    quadruple *coefs,
	    const int scaling,
	    const double eps_pivot,
	    const bool kernel_detection_all,
	    const int dim_augkern,
	    const double machine_eps_,
	    const bool assume_invertible,
	    const bool higher_precision);
template
void DissectionSolver<complex<quadruple>, quadruple, complex<double>, double>::
NumericFact(const int called,
	    complex<quadruple> *coefs,
	    const int scaling,
	    const double eps_pivot,
	    const bool kernel_detection_all,
	    const int dim_augkern,
	    const double machine_eps_,
	    const bool assume_invertible,
	    const bool higher_precision);

template
void DissectionSolver<float>::NumericFact(const int called,
					  float *coefs,
					  const int scaling,
					  const double eps_pivot,
					  const bool kernel_detection_all,
					  const int dim_augkern,
					  const double machine_eps_,
					  const bool assume_invertible,
					  const bool higher_precision);
template
void DissectionSolver<complex<float>, float>::
NumericFact(const int called,
	    complex<float> *coefs,
	    const int scaling,
	    const double eps_pivot,
	    const bool kernel_detection_all,
	    const int dim_augkern,
	    const double machine_eps_,
	    const bool assume_invertible,
	    const bool higher_precision);

//

template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
GetNullPivotIndices(int *pivots)
{
  int i0 = 0;
  for (int m = 0; m < _graph_colors; m++) {
    const int n0 = _kernel[m].dimension(); //_ptDA_kern_proj[m].dimension();
  // should be optimized : replaced by memcopy()?
    for (int i = 0; i < n0; i++, i0++) {
      pivots[i0] = _kernel[m].getKernListEq()[i]; //_ptDA_kern_list_eq[m][i];
    }
  } // loop : _graph_colors : m
}

template
void DissectionSolver<double>::GetNullPivotIndices(int *pivots);

template
void DissectionSolver<quadruple>::GetNullPivotIndices(int *pivots);

template
void DissectionSolver<complex<double>, double>::GetNullPivotIndices(int *pivots);

template
void DissectionSolver<complex<quadruple>, quadruple>::GetNullPivotIndices(int *pivots);

template
void DissectionSolver<double, double, quadruple, quadruple>::GetNullPivotIndices(int *pivots);

template
void DissectionSolver<quadruple, quadruple, double, double>::GetNullPivotIndices(int *pivots);

template
void DissectionSolver<float>::GetNullPivotIndices(int *pivots);

template
void DissectionSolver<complex<float>, float>::GetNullPivotIndices(int *pivots);

template
void DissectionSolver<float, float, double, double>::GetNullPivotIndices(int *pivots);

//

template<typename T, typename U, typename W, typename Z>
int DissectionSolver<T, U, W, Z>::
GetMaxColors()
{
  int nn = 0;
  for (int m = 0; m < _graph_colors; m++) {
    if (!_tridiagQueue[m]->tridiagSolver()) {
      nn++;
    }
  }
  return nn;
}

template
int DissectionSolver<double>::GetMaxColors();

template
int DissectionSolver<quadruple>::GetMaxColors();

template
int DissectionSolver<complex<double>, double>::GetMaxColors();

template
int DissectionSolver<complex<quadruple>, quadruple>::GetMaxColors();

template
int DissectionSolver<double, double, quadruple, quadruple>::GetMaxColors();

template
int DissectionSolver<quadruple, quadruple, double, double>::GetMaxColors();

template
int DissectionSolver<float>::GetMaxColors();

template
int DissectionSolver<complex<float>, float>::GetMaxColors();
//

template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
GetSmallestPivotIndices(const int n, int *pivots)
{
  int nn = 0;
  for (int m = 0; m < _graph_colors; m++) {
    if (!_tridiagQueue[m]->tridiagSolver()) {
      SquareBlockMatrix<T>* diag =
	&_dissectionMatrix[m][_btree[m]->selfIndex(1)]->diagBlock();
      const int *loc2glob = _btree[m]->getDiagLoc2Glob(1);
      vector<int> &permute = diag->getPermute();
      int ndiag = diag->dimension();
      if (n > ndiag) {
	diss_printf(_verbose, _fp, "%s %d : GetSmallestPivotIndices %d > %d\n",
		__FILE__, __LINE__, n, ndiag);
      }
      for (int k = (ndiag - 1); k >= (ndiag - n); k--, nn++) {
	pivots[nn] = loc2glob[permute[k]];
      }
    }
  }  // loop : m
}

template
void DissectionSolver<double>::GetSmallestPivotIndices(const int n,
						       int *pivots);

template
void DissectionSolver<quadruple>::GetSmallestPivotIndices(const int n, int *pivots);

template
void DissectionSolver<complex<double>, double>::GetSmallestPivotIndices(const int n, int *pivots);

template
void DissectionSolver<complex<quadruple>, quadruple>::GetSmallestPivotIndices(const int n, int *pivots);

template
void DissectionSolver<double, double, quadruple, quadruple>::GetSmallestPivotIndices(const int n, int *pivots);

template
void DissectionSolver<quadruple, quadruple, double, double>::GetSmallestPivotIndices(const int n, int *pivots);

template
void DissectionSolver<float>::GetSmallestPivotIndices(const int n,
						       int *pivots);

template
void DissectionSolver<complex<float>, float>::GetSmallestPivotIndices(const int n, int *pivots);

//
template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
GetKernelVectors(T *kern_basis)
{
  int ii = 0;
  for (int m = 0; m < _graph_colors; m++) {
    const int n0 = _kernel[m].dimension(); //_ptDA_kern_proj[m].dimension();
    // should be optimized : replaced by memcopy()?
    for (int j = 0; j < n0; j++) {
      for (int i = 0; i < _dim; i++, ii++) {
	kern_basis[ii] = _kernel[m].getKernBasis().addrCoefs()[ii];
      }
    }
  } // loop : _graph_colors : m
}

template
void DissectionSolver<double>::GetKernelVectors(double *kern_basis);

template
void DissectionSolver<quadruple>::GetKernelVectors(quadruple *kern_basis);

template
void DissectionSolver<complex<double>, double>::GetKernelVectors(complex<double> *kern_basis);

template
void DissectionSolver<complex<quadruple>, quadruple>::GetKernelVectors(complex<quadruple> *kern_basis);

template
void DissectionSolver<double, double, quadruple, quadruple>::GetKernelVectors(double *kern_basis);

template
void DissectionSolver<quadruple, quadruple, double, double>::GetKernelVectors(quadruple *kern_basis);

template
void DissectionSolver<float>::GetKernelVectors(float *kern_basis);

template
void DissectionSolver<complex<float>, float>::GetKernelVectors(complex<float> *kern_basis);
//

template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
GetTransKernelVectors(T *kernt_basis)
{
  int ii = 0;
  for (int m = 0; m < _graph_colors; m++) {
    const int n0 = _kernel[m].dimension(); //_ptDA_kern_proj[m].dimension();
    // should be optimized : replaced by memcopy()?
    for (int j = 0; j < n0; j++) {
      for (int i = 0; i < _dim; i++, ii++) {
	kernt_basis[ii] = _kernel[m].getTKernBasis().addrCoefs()[ii];
      }
    }
  } // loop : _graph_colors : m
}

template
void DissectionSolver<double>::GetTransKernelVectors(double *kernt_basis);


template
void DissectionSolver<quadruple>::GetTransKernelVectors(quadruple *kernt_basis);

template
void DissectionSolver<complex<double>, double>::GetTransKernelVectors(complex<double> *kernt_basis);

template
void DissectionSolver<complex<quadruple>, quadruple>::GetTransKernelVectors(complex<quadruple> *kernt_basis);

template
void DissectionSolver<double, double, quadruple, quadruple>::GetTransKernelVectors(double *kernt_basis);

template
void DissectionSolver<quadruple, quadruple, double, double>::GetTransKernelVectors(quadruple *kernt_basis);

template
void DissectionSolver<float>::GetTransKernelVectors(float *kernt_basis);

template
void DissectionSolver<complex<float>, float>::GetTransKernelVectors(complex<float> *kernt_basis);

//

template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
GetMatrixScaling(Z *weight)
{
  for (int i = 0; i < _dim; i++) {
    weight[i] = _precDiag[i];
  }
}

template
void DissectionSolver<double>::GetMatrixScaling(double *weight);

template
void DissectionSolver<complex<double>, double>::
GetMatrixScaling(double *weight);

template
void DissectionSolver<quadruple>::GetMatrixScaling(quadruple *weight);

template
void DissectionSolver<complex<quadruple>, quadruple>::
GetMatrixScaling(quadruple *weight);

template
void DissectionSolver<float>::GetMatrixScaling(float *weight);

template
void DissectionSolver<complex<float>, float>::
GetMatrixScaling(float *weight);
//

template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
ProjectionImageSingle(T *x, string name)
{
  for (int m = 0; m < _graph_colors; m++) {
    T *kern_basis0 = _kernel[m].getKernBasis().addrCoefs();
    T *kern_basis1 = (_ptDA->isSymmetric() ? 
		      _kernel[m].getKernBasis().addrCoefs() : 
		      _kernel[m].getTKernBasis().addrCoefs());
    T *kern_proj = (_ptDA->isSymmetric() ? 
		    _kernel[m].getKernProj().addrCoefs() : 
		    _kernel[m].getNTKernProj().addrCoefs());
    
    const int n0 = _kernel[m].dimension();

    VectorArray<T> xx(n0);
    for (int j = 0; j < n0; j++) {
      xx[j] = blas_dot<T>(_dim,
			  (kern_basis1 + j * _dim),
			  1,
			  x, 1);
    }
    diss_printf(_verbose, _fp,
		"%s %d : %s : check orthogonality of the given vector \n",
		__FILE__, __LINE__,
		name.c_str());
    for (int j = 0; j < n0; j++) {
      diss_printf(_verbose, _fp, "%.6e ", blas_abs<T, double>(xx[j]));
    }
    diss_printf(_verbose, _fp, "\n");
    blas_trsv<T>(CblasLower, CblasNoTrans, CblasUnit,
		 n0, kern_proj,
		 n0, xx.addrCoefs(), 1);
    {
      int itmp = 0;
      for (int j = 0; j < n0; j++, itmp += (n0 + 1)) {
	xx[j] *= kern_proj[itmp];
      }
    }
    
    blas_trsv<T>(CblasUpper, CblasNoTrans, CblasUnit,
		 n0, kern_proj,
		 n0, xx.addrCoefs(), 1);
    diss_printf(_verbose, _fp,
		"%s %d : %s : solution of kernel adjustment : ",
		__FILE__, __LINE__, name.c_str());
    for (int j = 0; j < n0; j++) {
      diss_printf(_verbose, _fp, "%.6e ", blas_abs<T, double>(xx[j]));
    } 
    diss_printf(_verbose, _fp, "\n");

    for (int i = 0; i < _dim; i++) {
      for (int j = 0; j < n0; j++) {
	x[i] -= xx[j] * kern_basis0[i + j * _dim];
      }
    }
    
    for (int j = 0; j < n0; j++) {
      xx[j] = blas_dot<T>(_dim,
			  (kern_basis1 + j * _dim), 1,
			  x, 1);
    }
    diss_printf(_verbose, _fp,
		"%s %d : %s : after projection, component of each kernel \n",
		__FILE__, __LINE__, name.c_str());
    for (int j = 0; j < n0; j++) {
      diss_printf(_verbose, _fp, "%.6e ", blas_abs<T, double>(xx[j]));
    }
    diss_printf(_verbose, _fp, "\n");
    xx.free();
  }  // loop : _graph_colors : m
}

template
void DissectionSolver<double>::ProjectionImageSingle(double *x, string name);

template
void DissectionSolver<quadruple>::ProjectionImageSingle(quadruple *x,
							string name);
template
void DissectionSolver<complex<double>, double>::
ProjectionImageSingle(complex<double> *x, string name);

template
void DissectionSolver<complex<quadruple>, quadruple>::
ProjectionImageSingle(complex<quadruple> *x, string name);

template
void DissectionSolver<float>::ProjectionImageSingle(float *x, string name);

template
void DissectionSolver<complex<float>, float>::
ProjectionImageSingle(complex<float> *x, string name);

template
void DissectionSolver<double, double, quadruple, quadruple>::ProjectionImageSingle(double *x, string name);

template
void DissectionSolver<quadruple, quadruple, double, double>::ProjectionImageSingle(quadruple *x, string name);

//
template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
ProjectionKernelOrthSingle(T *x, string name, bool isTrans)
{
  for (int m = 0; m < _graph_colors; m++) {
    const bool flag_trans = isTrans && (!_ptDA->isSymmetric());
    
    T *kern_basis = flag_trans ? _kernel[m].getTKernBasis().addrCoefs() :
                                   _kernel[m].getKernBasis().addrCoefs();
    T *kern_proj = flag_trans ? _kernel[m].getTKernProj().addrCoefs() : 
                                  _kernel[m].getKernProj().addrCoefs();
    const int n0 = _kernel[m].dimension();
  
    VectorArray<T> xx(n0);
    for (int j = 0; j < n0; j++) {
      xx[j] = blas_dot<T>(_dim,
			  (kern_basis + j * _dim),
			  1,
			  x, 1);
    }
    diss_printf(_verbose, _fp,
		"%s %d : %s : check orthogonality of the given vector\n",
		__FILE__, __LINE__, name.c_str());
    for (int j = 0; j < n0; j++) {
      diss_printf(_verbose, _fp, "%.6e ", blas_abs<T, double>(xx[j]));
    } 
    diss_printf(_verbose, _fp, "\n");

    blas_trsv<T>(CblasLower, CblasNoTrans, CblasUnit,
		 n0, kern_proj,
		 n0, xx.addrCoefs(), 1);
    {
      int itmp = 0;
      for (int j = 0; j < n0; j++, itmp += (n0 + 1)) {
	xx[j] *= kern_proj[itmp];
      }
    }
    
    blas_trsv<T>(CblasUpper, CblasNoTrans, CblasUnit,
		 n0, kern_proj,
		 n0, xx.addrCoefs(), 1);
    diss_printf(_verbose, _fp,
		"%s %d : %s : solution of kernel adjustment\n",
		__FILE__, __LINE__, name.c_str());
    for (int j = 0; j < n0; j++) {
      diss_printf(_verbose, _fp, "%.6e ", blas_abs<T, double>(xx[j]));
    } 
    diss_printf(_verbose, _fp, "\n");

    for (int i = 0; i < _dim; i++) {
      for (int j = 0; j < n0; j++) {
	x[i] -= xx[j] * kern_basis[i + j * _dim];
      }
    }
    
    for (int j = 0; j < n0; j++) {
      xx[j] = blas_dot<T>(_dim,
			  (kern_basis + j * _dim), 1,
			  x, 1);
    }
    diss_printf(_verbose, _fp,
		"%s %d : %s : after projection, component of each kernel\n",
		__FILE__, __LINE__, name.c_str());
    for (int j = 0; j < n0; j++) {
      diss_printf(_verbose, _fp, "%.6e ", blas_abs<T, double>(xx[j]));
    } 
    diss_printf(_verbose, _fp, "\n");
    xx.free();
  } // loop : _graph_colors : m
}

template
void DissectionSolver<double>::ProjectionKernelOrthSingle(double *x,
							  string name,
							  bool isTrans);
template
void DissectionSolver<quadruple>::ProjectionKernelOrthSingle(quadruple *x,
							     string name,
							     bool isTrans);

template
void DissectionSolver<complex<double>, double>::
ProjectionKernelOrthSingle(complex<double> *x, string name, bool isTrans);

template
void DissectionSolver<complex<quadruple>, quadruple>::
ProjectionKernelOrthSingle(complex<quadruple> *x,  string name, bool isTrans);

template
void DissectionSolver<float>::ProjectionKernelOrthSingle(float *x,
							 string name,
							 bool isTrans);

template
void DissectionSolver<complex<float>, float>::
ProjectionKernelOrthSingle(complex<float> *x, string name, bool isTrans);

//
//

template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
SpMV(const T *x, T *y, bool scaling_flag)
{
  if (_scaling && scaling_flag) {
    VectorArray<T> w(_dim);
    for (int i = 0; i < _dim; i++) {
      //      w[i] = x[i] / fromreal(_precDiag[i]); // conversion if T == complex<U>
      w[i] = x[i] / _precDiag[i]; // without conversion when complex
    }
    _ptDA->prod(w.addrCoefs(), y); 
    for (int i = 0; i < _dim; i++) {
      //      y[i] /= fromreal(_precDiag[i]);
      y[i] /= _precDiag[i];
    }
    w.free();
  }
  else {
    _ptDA->prod(x, y); 
  }
}

template
void DissectionSolver<double>::SpMV(const double *x, double *y, bool scaling_flag);

template
void DissectionSolver<quadruple>::SpMV(const quadruple *x, quadruple *y, bool scaling_flag);

template
void DissectionSolver<complex<double>, double>::
SpMV(const complex<double> *x, complex<double> *y, bool scaling_flag);

template
void DissectionSolver<complex<quadruple>, quadruple>::
SpMV(const complex<quadruple> *x, complex<quadruple> *y, bool scaling_flag);

template
void DissectionSolver<double, double, quadruple, quadruple>::SpMV(const double *x, double *y, bool scaling_flag);

template
void DissectionSolver<quadruple, quadruple, double, double>::SpMV(const quadruple *x, quadruple *y, bool scaling_flag);

template
void DissectionSolver<complex<double>, double, complex<quadruple>, quadruple>::
SpMV(const complex<double> *x, complex<double> *y, bool scaling_flag);

template
void DissectionSolver<complex<quadruple>, quadruple, complex<double>, double>::
SpMV(const complex<quadruple> *x, complex<quadruple> *y, bool scaling_flag);

template
void DissectionSolver<float>::SpMV(const float *x, float *y, bool scaling_flag);

template
void DissectionSolver<complex<float>, float>::
SpMV(const complex<float> *x, complex<float> *y, bool scaling_flag);
//

template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
SpMtV(const T *x, T *y, bool scaling_flag)
{
  if (_scaling && scaling_flag) {
    VectorArray<T> w(_dim);
    for (int i = 0; i < _dim; i++) {
      //      w[i] = x[i] / fromreal(_precDiag[i]);
      w[i] = x[i] / _precDiag[i];
    }
    _ptDA->prodt(w.addrCoefs(), y); 
    for (int i = 0; i < _dim; i++) {
      //      y[i] /= fromreal(_precDiag[i]);
      y[i] /= _precDiag[i];
    }
    w.free();
  }
  else {
    _ptDA->prodt(x, y); 
  }
}

template
void DissectionSolver<double>::SpMtV(const double *x, double *y, bool scaling_flag);

template
void DissectionSolver<quadruple>::SpMtV(const quadruple *x, quadruple *y, bool scaling_flag);

template
void DissectionSolver<complex<double>, double>::
SpMtV(const complex<double> *x, complex<double> *y, bool scaling_flag);

template
void DissectionSolver<complex<quadruple>, quadruple>::
SpMtV(const complex<quadruple> *x, complex<quadruple> *y, bool scaling_flag);

template
void DissectionSolver<float>::SpMtV(const float *x, float *y, bool scaling_flag);

template
void DissectionSolver<complex<float>, float>::
SpMtV(const complex<float> *x, complex<float> *y, bool scaling_flag);

//

template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
SolveSingle(T *x, bool projection,
	    bool isTrans,
	    bool isScaling,
	    const int nexcls)
{
  int n0, n1;
  n0 = 0;
  for (int m = 0; m < _graph_colors; m++) {
    n0 += _kernel[m].dimension();
  }
  n1 = 0;
  for (int m = 0; m < _graph_colors; m++) {
    n1 += _Schur[m].getSlduList().size();
  }
  diss_printf(_verbose, _fp, "%s %d : colors = %d n1 = %d nexcls = %d\n",
	      __FILE__, __LINE__, _graph_colors, n1, nexcls);

  VectorArray<T> yy(n1); // allocation with size 0 is defined in PlainMatirx.hpp
  elapsed_t t0_elapsed, t1_elapsed;
  get_realtime(&t0_elapsed);

  if (n0 > 0 && projection) {
    diss_printf(_verbose, _fp, 
		"%s %d : projection of the scaled RHS onto the image\n", 
		__FILE__, __LINE__);
    ProjectionKernelOrthSingle(x, "scaled RHS", true); // orthogonal to ker A^T
  }
  if (isScaling) {
    for (int i = 0; i < _dim; i++) {
      x[i] *= _precDiag[i];
    }
  }
  if (isTrans) {
    if (n1 > 0) {
      int i0 = 0;
      for (int m = 0; m < _graph_colors; m++) {
	int nn1 = _Schur[m].getSlduList().size();
	for (int i = 0; i < nn1; i++, i0++) {
	  const int ii = _Schur[m].getSldu().loc2glob()[i];
	  yy[i0] = x[ii];
	  x[ii] = _zero;  // the matrix is factorized only regular node, 
	  // elements for suspicious pivots are set to be 0
	} // loop : i
      } // loop : m
      diss_printf(_verbose, _fp, 
		  "%s %d : _ptDA_sldu[] are solved with n1=%d\n", 
		  __FILE__, __LINE__, n1);
      int mm0 = 0;
      for (int m = 0; m < _graph_colors; m++) {
	const int nnn1 = _Schur[m].getSlduList().size();
	if (nnn1 > nexcls) {
	  const int nn1 = nnn1 - nexcls;
	// y2 = A32^T - (A_11^-1 A_12)^T A_31^T
	  blas_gemv<T>(CblasTrans, _dim, nn1, _none,
		       _Schur[m].getScol().addrCoefs(), _dim,
		       x, 1, _one, yy.addrCoefs() + mm0, 1);
	  // S22 x2 = y2
	  blas_trsv<T>(CblasUpper, CblasTrans, CblasUnit,
		       nn1, 
		       _Schur[m].getSldu().addrCoefs(),
		       nn1, 
		       yy.addrCoefs() + mm0, 1);
	  for (int i = 0; i < nn1; i++) {
	    yy[mm0 + i] *= _Schur[m].getSldu()(i, i); // should be optimized
	  }
	  blas_trsv<T>(CblasLower, CblasTrans, CblasUnit,
		     nn1,
		       _Schur[m].getSldu().addrCoefs(),
		       nn1, 
		       yy.addrCoefs() + mm0, 1);
	  // A_11^T x_1 = A_31^T - A_21^T x_2
	  for (int i = 0; i < nn1; i++) {
	    const int ii = _Schur[m].getSlduList()[i];
	    const SparseMatrix<T> * const arow = _Schur[m].getArow();
	    for (int k = arow->ptRow(ii); k < arow->ptRow(ii + 1); k++) {
	      const int j = arow->indCol(k);
	      x[j] -= arow->Coef(k) * yy[mm0 + i];
	    }
	  } // loop : i
	}
	else {
	  // change singIdx0
	}
	mm0 += nnn1; // move to the next block
      } // loop : _graph_colors : m
    } // if (n1 > 0)
    SolveScaled(x, 1, true);
    if (n1 > 0) {
      int mm0 = 0;
      for (int m = 0; m < _graph_colors; m++) {
	const int nn1 = _Schur[m].getSlduList().size();
    	for (int i = 0; i < nn1; i++) {
	  const int ii = _Schur[m].getSldu().loc2glob()[i];
	  x[ii] = yy[mm0 + i];
	}
	mm0 += nn1; // move to the next block
      }  // loop : m
    }    // if (n1 > 0)
  }      // if (isTrans)
  else {
    if (n1 > 0) {
      int i0 = 0;
      for (int m = 0; m < _graph_colors; m++) {
	int nn1 = _Schur[m].getSlduList().size();
	for (int i = 0; i < nn1; i++, i0++) {
	  const int ii = _Schur[m].getSldu().loc2glob_left()[i];
	  yy[i0] = x[ii];
	  x[ii] = _zero;  // the matrix is factorized only regular node, 
	  // elements for suspicious pivots are set to be 0
	} // loop : i
	for (int i = 0; i < _kernel[m].getKernListEqLeft().size(); i++) {
	  const int ii = _kernel[m].getKernListEqLeft()[i];
	  x[ii] = _zero;
	}
      } // loop : m
    }
    for (int m = 0; m < _graph_colors; m++) {
      const int nnn1 = _Schur[m].getSlduList().size();
      diss_printf(_verbose, _fp, "%s %d : color = %d %d\n",
		  __FILE__, __LINE__, m, nnn1);
      if (nnn1 < nexcls) {
	if (!_tridiagQueue[m]->tridiagSolver()) {
	  SquareBlockMatrix<T>* diag =
	    &_dissectionMatrix[m][_btree[m]->selfIndex(1)]->diagBlock();
	  vector<int> &singIdx0 = diag->getSingIdx0();
	  int pseudosing;
	  if (singIdx0.size() > 0) {
	    pseudosing = std::min<int>((int)singIdx0.front(),
				       (int)singIdx0.back()) - 1;
	  }
	  else {
	    pseudosing = diag->dimension() - 1;
	  }
	  for (int i = nnn1; i < nexcls; i++, pseudosing--) {
	    singIdx0.push_back(pseudosing);
	  }
	}
      }
    }
    for (int m = 0; m < _graph_colors; m++) {
      if (!_tridiagQueue[m]->tridiagSolver()) {
	SquareBlockMatrix<T>* diag =
	    &_dissectionMatrix[m][_btree[m]->selfIndex(1)]->diagBlock();
	vector<int> &singIdx0 = diag->getSingIdx0();
	if (singIdx0.size() > 0) {
	  diss_printf(_verbose, _fp, "%s %d : color = %d modify dim = %d : ",
		  __FILE__, __LINE__, m, diag->dimension());
	  for (vector<int>::const_iterator it = singIdx0.begin();
	       it != singIdx0.end(); ++it) {
	    diss_printf(_verbose, _fp, " %d", *it);
	  }
	  diss_printf(_verbose, _fp, "\n");
	} //if (_verbose && (singIdx0.size() > 0))
      }
    }
    SolveScaled(x, 1, false);
    if (n1 > 0) {
      diss_printf(_verbose, _fp, 
		  "%s %d : _ptDA_sldu[] are solved with n1=%d\n", 
		  __FILE__, __LINE__, n1);
      int mm0 = 0;
      for (int m = 0; m < _graph_colors; m++) {
	const int nnn1 = _Schur[m].getSlduList().size();
	if (nnn1 > nexcls) {
	  const int nn1 = nnn1 - nexcls;
	  diss_printf(_verbose, _fp, "%s %d : sldu %d\n", __FILE__, __LINE__,
		      nn1);
	  for (int i = 0; i < nn1; i++) {
	    T stmp(0.0);
	    const int ii = _Schur[m].getSlduListLeft()[i];
	    const SparseMatrix<T> * const arow = _Schur[m].getArow();
	    for (int k = arow->ptRow(ii); k < arow->ptRow(ii + 1); k++) {
	      const int j = arow->indCol(k);
	      stmp += arow->Coef(k) * x[j];
	    }
	    yy[mm0 + i] -= stmp;
	  }
	  blas_trsv<T>(CblasLower, CblasNoTrans, CblasUnit,
		     nn1, _Schur[m].getSldu().addrCoefs(), nn1, 
		       yy.addrCoefs() + mm0, 1);
	  for (int i = 0; i < nn1; i++) {
	    yy[mm0 + i] *= _Schur[m].getSldu()(i, i); // should be optimized
	  }
	  blas_trsv<T>(CblasUpper, CblasNoTrans, CblasUnit,
		       nn1, _Schur[m].getSldu().addrCoefs(), nn1, 
		       yy.addrCoefs() + mm0, 1);
	  for (int i = 0; i < nn1; i++) {
	    const int ii = _Schur[m].getSldu().loc2glob()[i]; 
	    x[ii] = yy[mm0 + i];
	  }
	}
      } // loop : m
      {
	int mm0 = 0;
	for (int m = 0; m < _graph_colors; m++) {
	  const int nnn1 = _Schur[m].getSlduList().size();
	  if (nnn1 > nexcls) {
	    const int nn1 = nnn1 - nexcls;
	    blas_gemv<T>(CblasNoTrans, _dim, nn1, _none,
			 _Schur[m].getScol().addrCoefs(), _dim,
			 yy.addrCoefs() + mm0, 1, _one, x, 1);
	    mm0 += nnn1; // move to the next block
	  }
	  for (int i = 0; i < _kernel[m].getKernListEq().size(); i++) {
	    const int ii = _kernel[m].getKernListEq()[i];
	    x[ii] = _zero;
	  }
	} // loop : _graph_colors : m
      } // scope of mm0
    }   // if (n1 > 0)
    else {
      for (int m = 0; m < _graph_colors; m++) {
	const int nnn1 = _Schur[m].getSlduList().size();
	if (nnn1 < nexcls) {
	  if (!_tridiagQueue[m]->tridiagSolver()) {
	    SquareBlockMatrix<T>* diag =
	      &_dissectionMatrix[m][_btree[m]->selfIndex(1)]->diagBlock();
	    vector<int> &singIdx0 = diag->getSingIdx0();
	    if (singIdx0.size() > 0) {
	      diss_printf(_verbose, _fp, "%s %d : restore ",
			  __FILE__, __LINE__);
	      for (vector<int>::const_iterator it = singIdx0.begin();
		   it != singIdx0.end(); ++it) {
		diss_printf(_verbose, _fp, " %d", *it);
	      }
	      diss_printf(_verbose, _fp, "\n");
	      for (int i = nnn1; i < nexcls; i++) {
		singIdx0.pop_back();
	      }
	    }
	  } 
	} // if (nn1 < nexcls)
      }   // loop : m
    }     // if (n1 > 0)
  }       // if (isTrans)
  // adjust the kernel  : x = x - N (N{^T} N)^-1 N^{T} x
  if (isScaling) {
    for (int i = 0; i < _dim; i++) {
      x[i] *= _precDiag[i];
    }
  }
  //  
  if (_index_isolated.size() > 0) {
    for (vector<int>::const_iterator it = _index_isolated.begin();
	 it != _index_isolated.end(); ++it) {
      T atmp(1.0);
      for (int k = _ptDA->ptRow((*it)); k < _ptDA->ptRow((*it) + 1); k++) {
	if (_ptDA->indCol(k) == (*it)) {
	  atmp = _one / _ptDA->Coef(k);
	  break;
	}
      }
      x[(*it)] *= atmp;
    } // loop : it
  }   //  if (_index_isolated.size() > 0)
  if (n0 > 0 && projection) {
    diss_printf(_verbose, _fp, 
		"%s %d : projection orthogonal to the kernel\n", 
		__FILE__, __LINE__);
    ProjectionKernelOrthSingle(x, "scaled RHS", false); // orthogonal to ker A
  }
  get_realtime(&t1_elapsed);
}

template
void DissectionSolver<double>::SolveSingle(double *x, bool projection,
					   bool isTrans, bool isScaling,
					   const int nexcls);

template
void DissectionSolver<quadruple>::SolveSingle(quadruple *x, bool projection,
					      bool isTrans, bool isScaling,
					      const int nexcls);

template
void DissectionSolver<complex<double>, double>::
SolveSingle(complex<double> *x, bool projection, bool isTrans, bool isScaling,
	    const int nexcls);

template
void DissectionSolver<complex<quadruple>, quadruple>::
SolveSingle(complex<quadruple> *x, bool projection, bool isTrans,
	    bool isScaling, const int nexcls);
//
template
void DissectionSolver<double, double, quadruple, quadruple>::SolveSingle(double *x, bool projection,
					   bool isTrans, bool isScaling,
					   const int nexcls);

template
void DissectionSolver<quadruple, quadruple, double, double>::SolveSingle(quadruple *x, bool projection,
					      bool isTrans, bool isScaling,
					      const int nexcls);

template
void DissectionSolver<complex<double>, double, complex<quadruple>, quadruple>::
SolveSingle(complex<double> *x, bool projection, bool isTrans, bool isScaling,
	    const int nexcls);

template
void DissectionSolver<complex<quadruple>, quadruple, complex<double>, double>::
SolveSingle(complex<quadruple> *x, bool projection, bool isTrans,
	    bool isScaling, const int nexcls);

template
void DissectionSolver<float>::SolveSingle(float *x, bool projection,
					   bool isTrans, bool isScaling,
					   const int nexcls);

template
void DissectionSolver<complex<float>, float>::
SolveSingle(complex<float> *x, bool projection, bool isTrans, bool isScaling,
	    const int nexcls);

//
template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
SolveScaled(T *x, int nrhs, bool isTrans)
{
  for (int m = 0; m < _graph_colors; m++) {
    if (_tridiagQueue[m]->tridiagSolver()) {
      _tridiagQueue[m]->exec_fwbw(x, nrhs, isTrans);
    }
    else {
      //if (false) {
      if (nrhs == 1) {
        _dissectionQueue[m]->exec_fwbw_seq(x, nrhs, isTrans);
      }
      else {
	_dissectionQueue[m]->exec_fwbw(x, nrhs, isTrans);
      }
    }
  } // loop : m
}

template
void DissectionSolver<double>::SolveScaled(double *x, int nrhs, bool isTrans);

template
void DissectionSolver<complex<double>, double>::
SolveScaled(complex<double> *x, 
	    int nrhs, bool isTrans);

template
void DissectionSolver<quadruple>::SolveScaled(quadruple *x, int nrhs,
					      bool isTrans);

template
void DissectionSolver<complex<quadruple>, quadruple>::
SolveScaled(complex<quadruple> *x, 
	    int nrhs, bool isTrans);

template
void DissectionSolver<float>::SolveScaled(float *x, int nrhs, bool isTrans);


template
void DissectionSolver<complex<float>, float>::
SolveScaled(complex<float> *x, 
	    int nrhs, bool isTrans);
//

template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
SolveMulti(T *x, int nrhs, bool projection, 
	   bool isTrans, bool isScaling, const int nexcls)
{
  const int nsing = _nsing; //_singIdx.size();
  VectorArray<T> vtmp(nrhs);
  vector<int> nn1(_graph_colors);
    //  nn1.resize(_graph_colors);
  int i0, mm0, n0, n1;

  int itmp;
  ColumnMatrix<T> xx;
  ColumnMatrix<T> yy;

  elapsed_t t0_elapsed, t1_elapsed;

  n0 = 0;
  for (int m = 0; m < _graph_colors; m++) {
    n0 += _kernel[m].dimension();
  }

  n1 = 0;
  if (nsing > 0) {
    for (int m = 0; m < _graph_colors; m++) {
      nn1[m] = _Schur[m].getSlduList().size();
      n1 += nn1[m];
    }
  }
  get_realtime(&t0_elapsed);

  itmp = 0;
  xx.init(_dim, nrhs, x, false);
  if (isScaling) {
    for (int i = 0; i < _dim; i++) {
      for (int n = 0; n < nrhs; n++) {
	xx(i, n) *= _precDiag[i];
      }
    } // loop : n
  } 
  if (n0 > 0 && projection) {
    diss_printf(_verbose, _fp, 
		"%s %d : projection of the scaled RHS onto the image\n",
		__FILE__, __LINE__);
    for (int n = 0; n < nrhs; n++) {
      ProjectionKernelOrthSingle(xx.addrCoefs() + (n * _dim),
				 "scaled RHS", true); // orthogonal to ker A^T
    }
  }
  if (n1 > 0) {
    yy.init(n1, nrhs);
  }           // if (n1 > 0)
  if (isTrans) {
    if (n1 > 0) {    //    if (nsing > 0) {
      i0 = 0;
      for (int m = 0; m < _graph_colors; m++) {
	for (int i = 0; i < nn1[m]; i++, i0++) {
	  for (int n = 0; n < nrhs; n++) {
	    const int ii = _Schur[m].getSldu().loc2glob()[i];
	    yy(i0, n) = xx(ii, n);     //	  yy[n][i0] = xx[n][ii];
	    xx(ii, n) = _zero; 	 
	  }     // loop : n
	}       // loop : i
      }         // loop : m
      mm0 = 0;
      for (int m = 0; m < _graph_colors; m++) {
	if (nn1[m] > 0) {
	  blas_gemm<T>(CblasTrans, CblasNoTrans, nn1[m], nrhs, _dim, _none,
		       _Schur[m].getScol().addrCoefs(), _dim,
		       xx.addrCoefs(), _dim, _one, yy.addrCoefs() + mm0, n1);
	  mm0 += nn1[m];
	}
      }
      mm0 = 0;
      i0 = 0;
      for (int m = 0; m < _graph_colors; m++) {
	if (nn1[m] > 0) {
	  blas_trsm<T>(CblasLeft,
		       CblasUpper, CblasTrans, CblasUnit,
		       nn1[m], nrhs, _one, // alpha
		       _Schur[m].getSldu().addrCoefs(), nn1[m], 
		       yy.addrCoefs() + mm0, n1); // yy[i] = &y[i * n1]
	  for (int i = 0; i < nn1[m]; i++, i0++) {
	    for (int n = 0; n < nrhs; n++) {
	      //yy[n][i0] *= _Schur[m].getSldu()(i, i); // should be optimized
	      yy(i0, n) *= _Schur[m].getSldu()(i, i); // should be optimized
	    }
	  }
	  blas_trsm<T>(CblasLeft,
		       CblasLower, CblasTrans, CblasUnit,
		       nn1[m], nrhs, _one, // alpha
		       _Schur[m].getSldu().addrCoefs(), nn1[m], 
		       yy.addrCoefs() + mm0, n1);
	  mm0 += _Schur[m].getSlduList().size();
	}
      }
      i0 = 0;
      for (int m = 0; m < _graph_colors; m++) {
	for (int i = 0; i < nn1[m]; i++, i0++) {
	  const int ii = _Schur[m].getSlduList()[i];
	  const SparseMatrix<T> * const arow = _Schur[m].getArow();
	  for (int k = arow->ptRow(ii); k < arow->ptRow(ii + 1); k++) {
	    const int j = arow->indCol(k);
	    const T xtmp = arow->Coef(k);
	    for (int n = 0; n < nrhs; n++) {
	      xx(j, n) -= xtmp * yy(i0, n);  	   //	   xx[n][j] -= xtmp * yy[n][i0];
	    }
	  }   // loop : n
	}     // loop : i     
      }       // loop : m
    }         // if (nsing > 0)
    SolveScaled(xx.addrCoefs(), nrhs, true);
    if (nsing > 0) {
      i0 = 0;
      for (int m = 0; m < _graph_colors; m++) {
	for (int i = 0; i < nn1[m]; i++, i0++) {
	  for (int n = 0; n < nrhs; n++) {
	    const int ii = _Schur[m].getSldu().loc2glob()[i];
	    xx(ii, n) = yy(i0, n);   	    //	    xx[n][ii] = yy[n][i0];
	  }
	} // loop : i
      }   // loop : m
    }     // if (nsing > 0)
  }
  else {
    if (n1 > 0) {    //    if (nsing > 0) {
      i0 = 0;
      for (int m = 0; m < _graph_colors; m++) {
	for (int i = 0; i < nn1[m]; i++, i0++) {
	  for (int n = 0; n < nrhs; n++) {
	    const int ii = _Schur[m].getSldu().loc2glob_left()[i];
	    yy(i0, n) = xx(ii, n);     //	  yy[n][i0] = xx[n][ii];
	    xx(ii, n) = _zero; 	 
	  }     // loop : n
	}       // loop : i
      }
    }
    for (int m = 0; m < _graph_colors; m++) {
      const int nnn1 = _Schur[m].getSlduList().size();
      diss_printf(_verbose, _fp, "%s %d : color = %d %d\n",
		  __FILE__, __LINE__, m, nnn1);
      if (nnn1 < nexcls) {
	if (!_tridiagQueue[m]->tridiagSolver()) {
	  SquareBlockMatrix<T>* diag =
	    &_dissectionMatrix[m][_btree[m]->selfIndex(1)]->diagBlock();
	  vector<int> &singIdx0 = diag->getSingIdx0();
	  int pseudosing;
	  if (singIdx0.size() > 0) {
	    pseudosing = std::min<int>((int)singIdx0.front(),
				       (int)singIdx0.back()) - 1;
	  }
	  else {
	    pseudosing = diag->dimension() - 1;
	  }
	  for (int i = nnn1; i < nexcls; i++, pseudosing--) {
	    singIdx0.push_back(pseudosing);
	  }
	}
      }
    }
    for (int m = 0; m < _graph_colors; m++) {
      if (!_tridiagQueue[m]->tridiagSolver()) {
	SquareBlockMatrix<T>* diag =
	    &_dissectionMatrix[m][_btree[m]->selfIndex(1)]->diagBlock();
	vector<int> &singIdx0 = diag->getSingIdx0();
	if (singIdx0.size() > 0) {
	  diss_printf(_verbose, _fp, "%s %d : color = %d modify dim = %d : ",
		      __FILE__, __LINE__, m, diag->dimension());
	  for (vector<int>::const_iterator it = singIdx0.begin();
	       it != singIdx0.end(); ++it) {
	    diss_printf(_verbose, _fp, " %d", *it);
	  }
	  diss_printf(_verbose, _fp, "\n");
	} //if (singIdx0.size() > 0) 
      }
    }
    SolveScaled(xx.addrCoefs(), nrhs, false);
    if (nsing > 0) {
      i0 = 0;
      for (int m = 0; m < _graph_colors; m++) {
	for (int i = 0; i < nn1[m]; i++, i0++) {
	  vtmp.ZeroClear();
	  const int ii = _Schur[m].getSlduListLeft()[i];
	  const SparseMatrix<T> * const arow = _Schur[m].getArow();
	  for (int k = arow->ptRow(ii); k < arow->ptRow(ii + 1); k++) {
	    const int j = arow->indCol(k);
	    const T xtmp = arow->Coef(k);
	    for (int n = 0; n < nrhs; n++) {
	      vtmp[n] += xtmp * xx(j, n); //xx[n][j];
	    }
	  }   // loop : n
	  for (int n = 0; n < nrhs; n++) {
	    yy(i0, n) -= vtmp[n];  	    //	    yy[n][i0] -= vtmp[n];
	  }   // loop : n
	}     // loop : i     
      }       // loop : m
      mm0 = 0;
      i0 = 0;
      for (int m = 0; m < _graph_colors; m++) {
	if (nn1[m] > 0) {
	  blas_trsm<T>(CblasLeft,
		       CblasLower, CblasNoTrans, CblasUnit,
		       nn1[m], nrhs, _one, // alpha
		       _Schur[m].getSldu().addrCoefs(), nn1[m], 
		       yy.addrCoefs() + mm0, n1); // yy[i] = &y[i * n1]
	  for (int i = 0; i < nn1[m]; i++, i0++) {
	    for (int n = 0; n < nrhs; n++) {
	      //	      yy[n][i0] *= _Schur[m].getSldu()(i, i); // should be optimized
	      yy(i0, n) *= _Schur[m].getSldu()(i, i); // should be optimized
	    }
	  }
	  blas_trsm<T>(CblasLeft,
		       CblasUpper, CblasNoTrans, CblasUnit,
		       nn1[m], nrhs, _one, _Schur[m].getSldu().addrCoefs(),
		       nn1[m], 
		       yy.addrCoefs() + mm0, n1);
	  mm0 += _Schur[m].getSlduList().size();
	} // if (nn1[m] > 0)
      } // loop : m
      i0 = 0;
      for (int m = 0; m < _graph_colors; m++) {
	for (int i = 0; i < nn1[m]; i++, i0++) {
	  for (int n = 0; n < nrhs; n++) {
	    const int ii = _Schur[m].getSldu().loc2glob()[i];
	    xx(ii, n)= yy(i0, n);   	    //	    xx[n][ii] = yy[n][i0];
	  }
	}
      } // loop : m
      mm0 = 0;
      for (int m = 0; m < _graph_colors; m++) {
	if (nn1[m] > 0) {
	  blas_gemm<T>(CblasNoTrans, CblasNoTrans, _dim, nrhs, nn1[m], _none,
		       _Schur[m].getScol().addrCoefs(), _dim,
		       yy.addrCoefs() + mm0, n1, _one, xx.addrCoefs(), _dim);
	  mm0 += nn1[m];
	}
      }
    }
  } // if (isTrans) 
    // adjust the kernel  : x = x - N (N{^T} N)^-1 N^{T} x
  if (isScaling) {
    for (int n = 0; n < nrhs; n++) {
      for (int i = 0; i < _dim; i++) {
	xx(i, n) *= _precDiag[i];  	//	xx[n][i] *= _precDiag[i];
      }
    }
  }
  else { // inversion of isolated diagonal entries
    for (vector<int>::const_iterator it = _index_isolated.begin();
	 it != _index_isolated.end(); ++it) {
      T atmp(1.0);
      for (int k = _ptDA->ptRow((*it)); k < _ptDA->ptRow((*it) + 1); k++) {
	if (_ptDA->indCol(k) == (*it)) {
	  atmp = _one / _ptDA->Coef(k);
	  break;
	}
      }
      for (int n = 0; n < nrhs; n++) {
	xx((*it), n) *= atmp;  	//	xx[n][(*it)] *= atmp;
      }
    } // loop : it 
  }     // if (nsing > 0)
  if (n0 > 0 && projection) {
    diss_printf(_verbose, _fp,
		"%s %d : projection of the scaled solution onto the image\n",
	      __FILE__, __LINE__);
    for (int n = 0; n < nrhs; n++) {
      ProjectionKernelOrthSingle(xx.addrCoefs() + (n * _dim),
				 "scaled solution",
				 false); // orthogonal to ker A
    }
  }
  //  if (n1 > 0) {
    //    delete [] y;
    //yy.free();
 //  }
//  delete [] yy;
  //  delete [] xx;
  xx.free();
  yy.free();
  //  delete [] vtmp;
  //  delete [] nn1;
  get_realtime(&t1_elapsed);
}

template
void DissectionSolver<double>::SolveMulti(double *x, int nrhs,
					  bool projection, 
					  bool isTrans, bool isScaling,
					  const int nexcls);

template
void DissectionSolver<quadruple>::SolveMulti(quadruple *x, int nrhs,
					     bool projection, 
					     bool isTrans,
					     bool isScaling,
					     const int nexcls);

template
void DissectionSolver<complex<double>, double>::
SolveMulti(complex<double> *x, int nrhs, bool projection, 
	   bool isTrans, bool isScaling,
	   const int nexcls);

template
void DissectionSolver<complex<quadruple>, quadruple>::
SolveMulti(complex<quadruple> *x, int nrhs, bool projection, 
	   bool isTrans, bool isScaling,
	   const int nexcls);

template
void DissectionSolver<double, double, quadruple, quadruple>::
SolveMulti(double *x, int nrhs,
	   bool projection, 
	   bool isTrans, bool isScaling,
	   const int nexcls);

template
void DissectionSolver<complex<double>, double,
		      complex<quadruple>, quadruple>::
SolveMulti(complex<double> *x, int nrhs,
	   bool projection, 
	   bool isTrans,
	   bool isScaling,
	   const int nexcls);

template
void DissectionSolver<float>::SolveMulti(float *x, int nrhs,
					  bool projection, 
					  bool isTrans, bool isScaling,
					  const int nexcls);
template
void DissectionSolver<complex<float>, float>::
SolveMulti(complex<float> *x, int nrhs, bool projection, 
	   bool isTrans, bool isScaling,
	   const int nexcls);
//

template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
BuildSingCoefs(T *DSsingCoefs, 
	       SparseMatrix<T> *DCsCoefs,
	       T *DBsCoefs,
	       vector<int>& singIdx,
	       const bool isTrans)
{
  const bool flag_trans = isTrans && (!_ptDA->isSymmetric());
  T *DSsingCoefs0;
  const int nsing = singIdx.size();

  if (flag_trans) {
    DSsingCoefs0 = new T [nsing * nsing];
  }
  else {
    DSsingCoefs0 = DSsingCoefs;
  }
  _ptDA->extractSquareMatrix(DSsingCoefs0, singIdx);
  if (flag_trans) {
    for (int i = 0; i < nsing; i++) {
      for (int j = 0; j < nsing;j++) {
	DSsingCoefs[i + j * nsing] = DSsingCoefs0[j + i * nsing];
      }
    }
  }
  
  if ( _ptDA->isSymmetric() ) {
    // singIdx[] is sorted in increasing order
    // col_ind[] of symmetric matrix with CSR also in increasing order

    // Ssing = Asing - Cs*Am-1*(Bs)
    // ?? (*it == 0) where it.column() == singIdx[k] ??
    for (int i = 0; i < nsing; i++) {
      for (int j = i; j < nsing;j++) {
	//  [dense matrix As] - [sparse matrix Cs] * [dense matrix Am^-1 (Bs)]
	T stmp(0.0);
	for (int k = DCsCoefs->ptRow(i); k < DCsCoefs->ptRow(i + 1); k++) {
	  const int icol = DCsCoefs->indCol(k);
	  stmp += DCsCoefs->Coef(k) * DBsCoefs[icol + j * _dim];
	}
	DSsingCoefs[i + j * nsing] -= stmp;
      }
      // symmetrize
      for (int j = i + 1; j < nsing;j++) {
	DSsingCoefs[j + i * nsing] = DSsingCoefs[i + j * nsing];
      }
    }
  } //   if ( _ptDA->isSymmetric() ) {
  else {
    for (int i = 0; i < nsing; i++) {
      for (int j = 0; j < nsing;j++) {
	//  [dense matrix As] - [sparse matrix Cs] * [dense matrix Am^-1 (Bs)]
	T stmp(0.0);
	for (int k = DCsCoefs->ptRow(i); k < DCsCoefs->ptRow(i + 1); k++) {
	  const int icol = DCsCoefs->indCol(k);
	  stmp += DCsCoefs->Coef(k) * DBsCoefs[icol + j * _dim];
	}
	DSsingCoefs[i + j * nsing] -= stmp;
      }
    }
  } //   if ( _ptDA->isSymmetric() ) {
  diss_printf(_verbose, _fp,
	      "%s %d : Schur complement on the singular dofs set %d x %d\n",
	    __FILE__, __LINE__,
	    nsing, nsing);
  if (flag_trans) {
    delete [] DSsingCoefs0;
  }
}

template
void DissectionSolver<double>::
BuildSingCoefs(double *DSsingCoefs, 
	       SparseMatrix<double> *DCsCoefs,
	       double *DBsCoefs,
	       vector<int>& singIdx,
	       const bool isTrans);
template
void DissectionSolver<complex<double>, double>::
BuildSingCoefs(complex<double> *DSsingCoefs, 
	       SparseMatrix<complex<double> > *DCsCoefs,
	       complex<double> *DBsCoefs,
	       vector<int>& singIdx,
	       const bool isTrans);

template
void DissectionSolver<quadruple>::
BuildSingCoefs(quadruple *DSsingCoefs, 
	       SparseMatrix<quadruple> *DCsCoefs,
	       quadruple *DBsCoefs,
	       vector<int>& singIdx,
	       const bool isTrans);
template
void DissectionSolver<complex<quadruple>, quadruple>::
BuildSingCoefs(complex<quadruple> *DSsingCoefs, 
	       SparseMatrix<complex<quadruple> > *DCsCoefs,
	       complex<quadruple> *DBsCoefs,
	       vector<int>& singIdx,
	       const bool isTrans);

template
void DissectionSolver<float>::
BuildSingCoefs(float *DSsingCoefs, 
	       SparseMatrix<float> *DCsCoefs,
	       float *DBsCoefs,
	       vector<int>& singIdx,
	       const bool isTrans);
template
void DissectionSolver<complex<float>, float>::
BuildSingCoefs(complex<float> *DSsingCoefs, 
	       SparseMatrix<complex<float> > *DCsCoefs,
	       complex<float> *DBsCoefs,
	       vector<int>& singIdx,
	       const bool isTrans);
//

template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
BuildKernels(vector<int> &singIdx_, 
	     int n2,
	     SchurMatrix<T> &Schur,
	     KernelMatrix<T> &kernel)
{
  elapsed_t t0_elapsed, t1_elapsed;
  //  struct timespec ts0, ts1;
  const int n0 = singIdx_.size();
  const int n1 = n0 - n2; // regular nodes
  const int nsing = n2; // singIdx.size();
  ColumnMatrix<T> DBsCoefs(_dim, nsing);
  diss_printf(_verbose, _fp,
	      "%s %d : BuildKernels n0 = %d n1 = %d n2 = %d\n",
	      __FILE__, __LINE__, n0, n1, n2);
  if (n1 > 0) {
    //IndiceArray regularVal(n1);
    vector<int> regularVal;
    regularVal.resize(n1);
    //    SparseMatrix<T>* ptDA_acol;
    for (int i = 0; i < n1; i++) {
      regularVal[i] = singIdx_[i]; 
    }
    // to sort by increasing-order
    std::sort(regularVal.begin(), regularVal.end()); //regularVal.renumber();
    Schur.getAcol() = _ptDA->PartialCopyCSR(regularVal, n1, true);
    Schur.getArow() = _ptDA->PartialCopyCSR(regularVal, n1, false);
    // nullify [ A_12 / A_22 / A_23 ] -> [ A_12 / 0 / 0 ]
    for (int i = 0; i < n1; i++) {
      for (int k = Schur.getAcol()->ptRow(i);
	   k < Schur.getAcol()->ptRow(i + 1); k++) {
	// elimination from the CSR foramt is best for performance
	for (int m = 0; m < n0; m++) {
	  if (Schur.getAcol()->indCol(k) == singIdx_[m]) {
	    Schur.getAcol()->Coef(k) = _zero;
	  }
	}
      }
    }

    Schur.getScol().init(_dim, n1);
    for (int j = 0; j < n1; j++) {
      for (int i = 0; i < _dim; i++) {
	Schur.getScol()(i, j) = _zero;    // ZeroClear()
      }
      for (int k = Schur.getAcol()->ptRow(j);
	   k < Schur.getAcol()->ptRow(j + 1); k++) {
	const int icol = Schur.getAcol()->indCol(k);
	Schur.getScol()(icol, j) = Schur.getAcol()->Coef(k);
      }
    }
    //IndiceArray s_list_eq(n1); // local -> global mapping
    vector<int> s_list_eq;
    s_list_eq.resize(n1);
    for (int i = 0; i < n1; i++) {
      s_list_eq[i] = regularVal[i]; 
    }
    Schur.getSldu().init(s_list_eq);    //    ptDA_sldu.init
    // _Schur.getScol() = [ A_11^-1 A_12 / 0 / 0 ]
    SolveScaled(Schur.getScol().addrCoefs(), n1, false);
    //Schur.getScol().addrCoefs(), n1, false);

    BuildSingCoefs(Schur.getSldu().addrCoefs(), //ptDA_sldu.addrCoefs(),
		   Schur.getArow(),
		   Schur.getScol().addrCoefs(), //Schur.getScol().addrCoefs(),
		   regularVal);
    int *pivot_width = &(Schur.getSldu().getPivotWidth()[0]); //&ptDA_sldu.getPivotWidth()[0];
    T *d1 = Schur.getSldu().addr2x2();
    int *permute = &Schur.getSldu().getPermute()[0];
#if 0
    full_sym_2x2BK<T>(n1, Schur.getSldu().addrCoefs(), d1, pivot_width,
		      permute);
    //    fprintf(_fp, "Bunch-Kaufman permutation : ");
    diss_printf(_verbose, _fp, "%s %d : permutation for Schur complement: ",
		__FILE__, __LINE__);
    for (int i = 0; i < n1; i++) {
      diss_printf(_verbose, _fp, "%3d ", pivot_width[i]);
    }
    diss_printf(_verbose, _fp, "\n");
    for (int i = 0; i < n1; i++) {
      diss_printf(_verbose, _fp, "%3d ", permute[i]);
    }
    diss_printf(_verbose, _fp, "\n");
    Schur.getSlduList().resize(n1);
    for (int i = 0; i < n1; i++) {
      Schur.getSlduList()[i] = permute[i]; // ?? 20 Apr.2016 ??
    }
    int itmp2x2 = 0;
    for (int i = 0; i < n1; i++) {
      if (pivot_width[i] == 20) {
	itmp2x2++;
      }
    }
    Schur.getSldu().getPivot2x2().resize(itmp2x2);
    itmp2x2 = 0;
    for (int i = 0; i < n1; i++) {
      if (pivot_width[i] == 20) {
	Schur.getSldu().getPivot2x2()[itmp2x2] = i;
	itmp2x2++;
      }
    }
    // 1x1 + 2x2
    //    delete ptDA_acol;
#else
  {
    int nn0, n0 = 0;
    double fop, pivot_ref = 1.0; // should be automatically defined in ldlt/ldu
    double eps_piv = machine_epsilon<double, double>();
    if (_ptDA->isSymmetric()) {
      full_ldlt_permute<T, U>(&nn0, n0, n1,
			      Schur.getSldu().addrCoefs(), n1,
			      &pivot_ref,
			      permute, eps_piv, &fop);
    }
    else {
      full_ldu_permute<T, U>(&nn0, n0, n1,
			     Schur.getSldu().addrCoefs(), n1,
			     &pivot_ref,
			     permute, eps_piv, &fop);
    }
    diss_printf(_verbose, _fp, 
		"%s %d : factorization with eps_piv = %g : dim kern = %d\n",
		__FILE__, __LINE__, eps_piv, nn0);
    diss_printf(_verbose, _fp, "permute[] = ");
    for (int i = 0; i < n1; i++) {
      diss_printf(_verbose, _fp, "%d ", permute[i]);
    }
    diss_printf(_verbose, _fp, "\n");
  }
  Schur.getSlduList().resize(n1);
  for (int i = 0; i < n1; i++) {
    Schur.getSlduList()[i] = permute[i]; 
  }
#endif
  }  // if (n1 > 0)
  else { // if (n1 == 0)
    Schur.getAcol() = new SparseMatrix<T>(); // dummy allocation
    Schur.getArow() = new SparseMatrix<T>(); // dummy allocation
    Schur.getSlduList().clear();
  } // if (n1 > 0)
  
  if (n2 == 0) {
    kernel.set_dimension(0);
    kernel.getKernProj().init(0);
    kernel.getKernBasis().init(0, 0);
  }
  else {
    SparseMatrix<T> *DCsCoefs;
    vector<int> singIdx;
    singIdx.resize(n2);
    for (int i = 0; i < n2; i++) {
      singIdx[i] = singIdx_[n0 - n2 + i]; 
    }
    // to sort by increasing-order
    std::sort(singIdx.begin(), singIdx.end()); //    singIdx.renumber();
    //DCsCoefs = PartialCopyCSR(_ptDA, singIdx, nsing, true);
    DCsCoefs = _ptDA->PartialCopyCSR(singIdx, nsing, true); // 28 Dec.2015
    // clear element of DBsCoefs belonging to diagnal block
    for (int i = 0; i < nsing; i++) {
      for (int k = DCsCoefs->ptRow(i); k < DCsCoefs->ptRow(i + 1); k++) {
	// elimination from the CSR foramt is best for performance
	for (int m = 0; m < nsing; m++) {
	  if (DCsCoefs->indCol(k) == singIdx[m]) {
	    DCsCoefs->Coef(k) = _zero;
	  }
	}
      }
    }
    diss_printf(_verbose, _fp,
		"%s %d : Cs matrix: row = %d\n",
		__FILE__, __LINE__, nsing);
    // generate DBsCoefs[]
    // copy from lower part (sparse) with transporse to upper part (dense)
    DBsCoefs.ZeroClear();
      for (int j = 0; j < nsing; j++) {
      for (int k = DCsCoefs->ptRow(j); k < DCsCoefs->ptRow(j + 1); k++) {
	const int icol = DCsCoefs->indCol(k);
	DBsCoefs(icol, j) = DCsCoefs->Coef(k);
      }
    }
    // access with stride size _dim
    for (int i = 0; i < nsing; i++) {
      for (int j = 0; j < nsing; j++) {
	DBsCoefs(singIdx[i], j) = _zero;
      }
    }
    diss_printf(_verbose, _fp,
		"%s %d Bs matrix: singular columns of global matrix %d x %d\n",
		__FILE__, __LINE__,
		_dim, nsing);
  // Compute Am-1*Bs
  // scaling flag is off : internal solver A_11^-1 [A^12 A^13]
    if (n1 > 0) {
      ColumnMatrix<T> yy(n1, nsing);
      ColumnMatrix<T> yy1(n1, nsing);
      for (int m = 0; m < nsing; m++) {
	for (int i = 0; i < n1; i++) {
	  int ii = Schur.getSldu().loc2glob()[i];
	  // yy[i + m * n1] = DBsCoefs[ii + m * _dim]; save data of singular nds
	  // DBsCoefs[ii + m * _dim] = zero; SolveScaled works for regular nodes
	  yy(i, m) = DBsCoefs(ii, m);
	  DBsCoefs(ii, m) = _zero;
	}
      } // loop : j
      get_realtime(&t0_elapsed);
      SolveScaled(DBsCoefs.addrCoefs(), nsing, false);
      get_realtime(&t1_elapsed);
      diss_printf(_verbose, _fp,
		  "%s %d : dissection solve : %d RHS (sec.) = %.6f\n", 
		  __FILE__, __LINE__, nsing,
		  convert_time(t1_elapsed, t0_elapsed));
      for (int m = 0; m < nsing; m++) {    // Block SpMV
	for (int i = 0; i < n1; i++) {
	  T stmp(0.0); // = _zero;
	  const int ii = Schur.getSlduList()[i];
	  for (int k = Schur.getArow()->ptRow(ii);
	       k < Schur.getArow()->ptRow(ii + 1); 
	       k++) {
	    const int j = Schur.getArow()->indCol(k);
	    stmp += Schur.getArow()->Coef(k) * DBsCoefs(j, m);
	  }
	  yy(i, m) -= stmp;
	}
      } // loop : m
      blas_trsm<T>(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
		   n1, nsing, _one, // alpha
		   Schur.getSldu().addrCoefs(), n1, 
		   yy.addrCoefs(), n1);
      
      for (int i = 0; i < n1; i++) {
	T stmp = Schur.getSldu()(i, i);
	for (int m = 0; m < nsing; m++) {
	  yy1(i, m) = yy(i, m) * stmp; // should be optimized
	} // loop : m
      } 
      for (int i = 0; i < Schur.getSldu().getPivot2x2().size(); i++) {
	const int ii = Schur.getSldu().getPivot2x2()[i];
	T *d1 = Schur.getSldu().addr2x2();
	for (int m = 0; m < nsing; m++) {
	  yy1(ii, m) += d1[m] * yy((i + 1), m);
	  yy1((ii + 1), m) += d1[m + 1] * yy(i, m);
	}
      }
      yy.copy(yy1);
      blas_trsm<T>(CblasLeft, CblasUpper, CblasNoTrans, CblasUnit,
		   n1, nsing, _one,
		   Schur.getSldu().addrCoefs(), n1, 
		   yy.addrCoefs(), n1);
      for (int i = 0; i < n1; i++) {
	const int ii = Schur.getSldu().loc2glob()[i];
	for (int m = 0; m < nsing; m++) {
	  DBsCoefs(ii,  m) = yy(i, m);
	} // loop : m
      } 
      blas_gemm<T>(CblasNoTrans, CblasNoTrans,
		   _dim, nsing, n1, 
		   _none, // alpha
		   Schur.getScol().addrCoefs(), _dim,
		   yy.addrCoefs(), n1,
		   _one,  // beta
		   DBsCoefs.addrCoefs(), _dim);
    } // if (n1 > 0)
    else {
      get_realtime(&t0_elapsed);
      SolveScaled(DBsCoefs.addrCoefs(), nsing, false);
      get_realtime(&t1_elapsed);
      diss_printf(_verbose, _fp,
		  "%s %d : dissection solve : %d RHS (sec.) = %.6f\n", 
		  __FILE__, __LINE__, nsing,
		  convert_time(t1_elapsed, t0_elapsed));
    }
      // clear diagonal blocks of Bs in the place of Schur complement
    for (int i = 0; i < nsing; i++) {
      DBsCoefs(singIdx[i], i) = _none;
    }
    
    kernel.getKernListEq().resize(nsing); // = IndiceArray(nsing);
    for (int i = 0; i < nsing; i++) {
      kernel.getKernListEq()[i]= singIdx[i];
    }
    diss_printf(_verbose, _fp,
		 "%s %d : kern_list_eq[] = ",
		 __FILE__, __LINE__);
    for (int i = 0; i < kernel.getKernListEq().size(); i++) {
      diss_printf(_verbose, _fp, "%d ", kernel.getKernListEq()[i]);
    }
    diss_printf(_verbose, _fp, "\n");

    kernel.getKernBasis().init(_dim, nsing);
    T *kern_basis = kernel.getKernBasis().addrCoefs();
    kernel.getKernBasis().copy(DBsCoefs);
  // normalize each kernel_basis
    for (int j = 0; j < nsing; j++) {
      U stmp(0.0);
      const U one(1.0);
      stmp = one / blas_l2norm<T, U>(_dim, (kern_basis + (j * _dim)), 1);
      blas_scal2<T, U>(_dim, stmp, (kern_basis + (j * _dim)), 1);
    }
    diss_printf(_verbose, _fp,
		"%s %d : check kern_basis belonging to the kernel of A\n",
		__FILE__, __LINE__);
    VectorArray<T> v(_dim);
    diss_printf(_verbose, _fp,
		"%s %d : orthogonality of the kernel vectors\n",
	      __FILE__, __LINE__);
    for (int j = 0; j < nsing; j++) {
      _ptDA->prod((kern_basis + (j * _dim)), v.addrCoefs());
      double norm_l2, norm_infty;
      calc_relative_norm<T>(&norm_l2, &norm_infty, v.addrCoefs(),
			    kern_basis + (j * _dim), _dim);
      diss_printf(_verbose, _fp,
		  "%d-th kernel scaled : l2_norm = %.6e, infty_norm = %.6e / ",
		  j, norm_l2, norm_infty);
      if(_scaling) {
	calc_relative_normscaled<T, U>(&norm_l2, &norm_infty,
				       v.addrCoefs(),
				       kern_basis + (j * _dim),
				       &_precDiag[0], _dim);
	diss_printf(_verbose, _fp,
		    "original : l2_norm = %.6e, infty_norm = %.6e\n",
		    norm_l2, norm_infty);
      } // if (_scaling)
      diss_printf(_verbose, _fp, "\n");
    }
#ifdef DEBUG_DATA
    fout.close();
#endif

    if(_scaling) {
      for (int j = 0; j < nsing; j++) {
	const int jtmp = j * _dim;
	for (int i = 0; i < _dim; i++) {
	  kern_basis[i + jtmp] *= _precDiag[i];
	}
      }
    }
#ifdef NORMALIZE_KERNEL_BASIS
  // normalize each kernel_basis
    for (int j = 0; j < nsing; j++) {
      U stmp(0.0);
      const U one(1.0);
      stmp = one / blas_l2norm<T, U>(_dim, (kern_basis + (j * _dim)), 1);
      blas_scal2<T, U>(_dim, stmp, (kern_basis + (j * _dim)), 1);
    }
#endif
    kernel.set_dimension(nsing);
    kernel.getKernProj().init(nsing);
    
    for (int i = 0; i < nsing; i++) {
      for(int j = i; j < nsing; j++) {
	kernel.getKernProj()(i, j) = blas_dot<T>(_dim,
						 (kern_basis + (i * _dim)), 1,
						 (kern_basis + (j * _dim)), 1);
      }
      // symmetrize
      for (int j = (i + 1); j < nsing; j++) {
	kernel.getKernProj()(j, i) = kernel.getKernProj()(i, j);
      }
    }
    full_ldlh<T>(nsing, kernel.getKernProj().addrCoefs(), nsing);
  // inverse of diagonal part is also storead in the factorized matrix
    delete DCsCoefs;
  } //  if (n2 == 0) 
  //  delete [] DBsCoefs;
}

template
void DissectionSolver<double>::
BuildKernels(vector<int> &singIdx_, 
	     int n2,
	     SchurMatrix<double> &Schur,
	     KernelMatrix<double> &kernel);

template
void DissectionSolver<complex<double>, double>::
BuildKernels(vector<int> &singIdx_, 
	     int n2,
	     SchurMatrix<complex<double> > &Schur,
	     KernelMatrix<complex<double> > &kernel);

template
void DissectionSolver<quadruple>::
BuildKernels(vector<int> &singIdx_, 
	     int n2,
	     SchurMatrix<quadruple> &Schur,
	     KernelMatrix<quadruple> &kernel);

template
void DissectionSolver<complex<quadruple>, quadruple>::
BuildKernels(vector<int> &singIdx_, 
	     int n2,
	     SchurMatrix<complex<quadruple> > &Schur,
	     KernelMatrix<complex<quadruple> > &kernel);

template
void DissectionSolver<float>::
BuildKernels(vector<int> &singIdx_, 
	     int n2,
	     SchurMatrix<float> &Schur,
	     KernelMatrix<float> &kernel);

template
void DissectionSolver<complex<float>, float>::
BuildKernels(vector<int> &singIdx_, 
	     int n2,
	     SchurMatrix<complex<float> > &Schur,
	     KernelMatrix<complex<float> > &kernel);

//

template<typename T, typename U, typename W, typename Z>
void DissectionSolver<T, U, W, Z>::
BuildKernelsDetection(int &n0,
		      vector<int> &singIdx, 
		      vector<vector<int> >& augkern_indexes,
		      vector<SquareBlockMatrix<T>* >& diags, 
		      const double eps_piv,
		      const int dim_augkern,
		      SchurMatrix<T> &Schur,
		      KernelMatrix<T> &kernel,
		      const bool enableDetection)
{
  elapsed_t t0_elapsed, t1_elapsed, t2_elapsed, t3_elapsed;
  //  struct timespec ts0, ts1, ts2, ts3;

  //  vector<int>& singIdx = singIdx; // singIdx is inherited from older version
  const int nsing = singIdx.size();
  ColumnMatrix<T> DBsCoefs(_dim, nsing);
  ColumnMatrix<T> DCsCoefs(_dim, nsing);
  ColumnMatrix<T> DSsingCoefs(nsing, nsing);
  ColumnMatrix<T> DSsingCoefs0(nsing, nsing);
  ColumnMatrix<T> DSsingCoefs1;
  ColumnMatrix<T> DSsingCoefs2(nsing, nsing);
  ColumnMatrix<T> DSsingCoefs3(nsing, nsing);
  ColumnMatrix<T> DS0(nsing, nsing);
  ColumnMatrix<T> DS2(nsing, nsing);
  double pivot_ref;
  int *permute0 = new int[nsing];
  VectorArray<T> u(_dim);
  VectorArray<T> v(_dim);
  //
  VectorArray<T> save_diag(dim_augkern);
  {
    int i = 0;
    diss_printf(_verbose, _fp,
		"%s %d : nullifying diagonal entries\n",
		__FILE__, __LINE__);
    typename vector<SquareBlockMatrix<T>* >::const_iterator kt = diags.begin();
    for (vector<vector<int> >::const_iterator it = augkern_indexes.begin();
	 it != augkern_indexes.end(); ++it, ++kt) {
      SquareBlockMatrix<T> &diag = *(*kt);
      const int n0_diag = diag.rank();
      diss_printf(_verbose, _fp, "diagonal entries modifed : %d -> %d\n",
		  n0_diag, n0_diag - (int)(*it).size());
      for (vector<int>::const_iterator jt = (*it).begin(); jt != (*it).end();
	   ++jt) {
	diss_printf(_verbose, _fp, "%d : %s\n",
		    *jt, tostring<T>(diag.diag(*jt)).c_str());
	save_diag[i++] = diag.diag((*jt));
	diag.diag((*jt)) = _zero;
      }
      diag.set_rank(n0_diag - (*it).size());
    }
  }
  // singIdx is sorted with increasing-order
  Schur.getAcol() = _ptDA->PartialCopyCSR(singIdx, nsing, true);  // 28 Dec.2015
  Schur.getArow() = _ptDA->PartialCopyCSR(singIdx, nsing, false); // 28 Dec.2015
  //  nullify diagonal blocks of singular nodes
  //  [A_21 A_22 A_23 / A_31 A_32 A_33 ] -> [A_21 0  0 / A_31 0  0 ]
  // elimination of CSR entries is best for performance, but now replacing by 0
  {   // scope of acol
    SparseMatrix<T> *acol = Schur.getAcol();
    SparseMatrix<T> *arow = Schur.getArow();
    for (int i = 0; i < nsing; i++) {
      for (int k = acol->ptRow(i); k < acol->ptRow(i + 1); k++) {
	for (int n = 0; n < nsing; n++) {
	  if (acol->indCol(k) == singIdx[n]) {
	    acol->Coef(k) = _zero;
	  }
	}
      }
      for (int k = arow->ptRow(i); k < arow->ptRow(i + 1); k++) {
	for (int n = 0; n < nsing; n++) {
	  if (arow->indCol(k) == singIdx[n]) {
	    arow->Coef(k) = _zero;
	  }
	}
      }
    }
    diss_printf(_verbose, _fp, 
		"%s %d : Cs matrix: from global sparse matrix : row = %d\n",
		__FILE__, __LINE__, nsing);
    // generate DBsCoefs[]
    // copy from lower part (sparse) with transporse to upper part (dense)
    DBsCoefs.ZeroClear();
    for (int j = 0; j < nsing; j++) {
      for (int k = acol->ptRow(j); k < acol->ptRow(j + 1); k++) {
	const int icol = acol->indCol(k);
	DBsCoefs(icol, j) = acol->Coef(k);
      }
    }
    DCsCoefs.ZeroClear();
    for (int j = 0; j < nsing; j++) {
      for (int k = arow->ptRow(j); k < arow->ptRow(j + 1); k++) {
	const int icol = arow->indCol(k);
	DCsCoefs(icol, j) = arow->Coef(k);
      }
    }
  }  // scope of acol
  diss_printf(_verbose, _fp,
	      "%s %d : Bs matrix: from global sparse matrix %d x %d\n",
	      __FILE__, __LINE__, _dim, nsing);

  // Compute Am-1*Bs
  // scaling flag is off : internal solver A_11^-1 [A^12 A^13]
  get_realtime(&t0_elapsed);
  //  clock_gettime(CLOCK_REALdoubleIME, &ts0);
  SolveScaled(DBsCoefs.addrCoefs(), nsing, false);
  if (!_ptDA->isSymmetric()) {
    SolveScaled(DCsCoefs.addrCoefs(), nsing, true);
  }
  get_realtime(&t1_elapsed);
  //  clock_gettime(CLOCK_REALTIME, &ts1);
  diss_printf(_verbose, _fp,
	      "%s %d : dissection solve : %d RHS (sec.) = %.6f\n", 
	      __FILE__, __LINE__,
	      nsing, convert_time(t1_elapsed, t0_elapsed));

  BuildSingCoefs(DSsingCoefs.addrCoefs(),
		 Schur.getArow(), DBsCoefs.addrCoefs(), singIdx);
  if (!_ptDA->isSymmetric()) {
    BuildSingCoefs(DSsingCoefs2.addrCoefs(),
		   Schur.getAcol(), DCsCoefs.addrCoefs(), singIdx, true);
  }
  // copy DSsingCoefs
  DS0.copy(DSsingCoefs);
  DS2.copy(DSsingCoefs2);
  // set candidate of kernel vectors before applying Schur complement on 
  // regular part of suspicious pivots
  for (int i = 0; i < nsing; i++) {
    DBsCoefs(singIdx[i], i) = _none;
    DCsCoefs(singIdx[i], i) = _none;
  }
  get_realtime(&t2_elapsed);
  DSsingCoefs0.copy(DSsingCoefs);

  diss_printf(_verbose, _fp,
	      "%s %d : DSsingCoefs0[]\n",
	      __FILE__, __LINE__);
  for (int i = 0; i < nsing; i++) {
    for (int j = 0; j < nsing; j++) {
      diss_printf(_verbose, _fp, "%s ",
		  tostring<T>(DSsingCoefs0(i, j)).c_str());
    }
    diss_printf(_verbose, _fp, "\n");
  }
  if (!_ptDA->isSymmetric()) {
    diss_printf(_verbose, _fp,
		"%s %d : transposed DSsingCoefs2[]\n",
		__FILE__, __LINE__);
    for (int i = 0; i < nsing; i++) {
      for (int j = 0; j < nsing; j++) {
	diss_printf(_verbose, _fp, "%s ",
		    tostring<T>(DSsingCoefs2(j, i)).c_str());
      }
      diss_printf(_verbose, _fp, "\n");
    }
    diss_printf(_verbose, _fp,
		"%s %d : unsymmetry\n",
		__FILE__, __LINE__);
    for (int i = 0; i < nsing; i++) {
      for (int j = 0; j < nsing; j++) {
	diss_printf(_verbose, _fp, "%s ",
		    tostring<T>(DSsingCoefs0(i, j) - DSsingCoefs2(j, i)).c_str());
      }
      diss_printf(_verbose, _fp, "\n");
    }
  }
  bool flag_unsym_permute = false;
  if (enableDetection) {
  pivot_ref = 0.0;
  for (int i = 0; i < nsing; i++) { // find maximum of abs value of diagonal
    const double stmp = blas_abs<T, double>(DSsingCoefs0(i, i));
    pivot_ref = stmp > pivot_ref ? stmp : pivot_ref;
  }
  
  diss_printf(_verbose, _fp,
	      "%s %d : pviot_ref = %.8e\n",
	    __FILE__, __LINE__, pivot_ref);

  // permute0[] is initialized in ddfull_sym_gauss_part
  int nn = 0;
  {
    int nn0;
    double fop;
    if (_ptDA->isSymmetric()) {
      full_ldlt_permute<T, U>(&nn0, nn, nsing,
			      DSsingCoefs0.addrCoefs(), nsing,
			      &pivot_ref,
			      permute0, eps_piv, &fop);
    }
    else {
      full_ldu_permute<T, U>(&nn0, nn, nsing,
			     DSsingCoefs0.addrCoefs(), nsing,
			     &pivot_ref,
			     permute0, eps_piv, &fop);
    }
    nn = nn0;
  }
  diss_printf(_verbose, _fp, 
	      "%s %d : factorization with eps_piv = %g : dim kern = %d\n",
	      __FILE__, __LINE__, eps_piv, nn);
  diss_printf(_verbose, _fp, "permute[] = ");
  for (int i = 0; i < nsing; i++) {
    diss_printf(_verbose, _fp, "%d ", permute0[i]);
  }
  diss_printf(_verbose, _fp, "\n");
#if 0
  diss_printf(_verbose, _fp,
	      "%s %d : DSsingCoefs0[]\n",
	      __FILE__, __LINE__);
  for (int i = 0; i < nsing; i++) {
    for (int j = 0; j < nsing; j++) {
      diss_printf(_verbose, _fp, "%s ",
		  tostring<T>(DSsingCoefs0(i, j)).c_str());
    }
    diss_printf(_verbose, _fp, "\n");
  }
#endif
  if (nn == 0) {
    n0 = 0;
  }
  else {
    const int nsing1 = nn + dim_augkern;
    const int nsing0 = nsing - nsing1;
    DSsingCoefs1.init(nsing1, nsing1); 
    if (nsing0 == 0) {
      DSsingCoefs1.copy(DSsingCoefs);
    }
    else{ // recompute smaller Schur complement with size nsing1

      // nullify rows and columns
      for (int j = nsing0; j < nsing; j++) {
	for (int i = 0; i < nsing; i++) {
	  DSsingCoefs0(i, j) = _zero;
	  DSsingCoefs0(j, i) = _zero;
	}
      }
    // copy upper block
      for (int j = nsing0; j < nsing; j++) {
	for (int i = 0; i < nsing0; i++) {
	  DSsingCoefs0(i, j) = DSsingCoefs(permute0[i], permute0[j]);
	}
      }
      if (_ptDA->isSymmetric()) {
	int ii = nsing0;
	for (int i = 0; i < nsing1; i++, ii++) { 
	  int jj = nsing0;
	  for (int j = 0; j <= i; j++, jj++) {
	    DSsingCoefs1(i, j) = DSsingCoefs(permute0[ii], permute0[jj]);
	  }
	  // symmetrize
	  for (int j = 0; j < i; j++) {
	    DSsingCoefs1(j, i) = DSsingCoefs1(i, j);
	  }
	}
      }
      else {
	int ii = nsing0;
	for (int i = 0; i < nsing1; i++, ii++) { 
	  int jj = nsing0;
	  for (int j = 0; j < nsing1; j++, jj++) {
	    DSsingCoefs1(i, j) =  DSsingCoefs(permute0[ii], permute0[jj]);
	  }
	}
	// copy lower block
	for (int j = 0; j < nsing0; j++) {
	  for (int i = nsing0; i < nsing; i++) {
	    DSsingCoefs0(i, j) = DSsingCoefs(permute0[i], permute0[j]);
	  }
	}
      }
      blas_trsm<T>(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
		   nsing0, nsing1, 
		   _one, // alpha
		   DSsingCoefs0.addrCoefs(), nsing, 
		   DSsingCoefs0.addrCoefs() + (nsing0 * nsing), nsing);
      if (_ptDA->isSymmetric()) {
	for (int j = nsing0; j < nsing; j++) {
	  for (int i = 0; i < nsing0; i++) {
	    DSsingCoefs0(j, i) = DSsingCoefs0(i, j);
	  }
	}
      }
      else {
	blas_trsm<T>(CblasRight, CblasUpper, CblasNoTrans, CblasUnit,
		     nsing1, nsing0, 
		     _one, // alpha
		     DSsingCoefs0.addrCoefs(), nsing, 
		     DSsingCoefs0.addrCoefs() + nsing0, nsing);
      } //  if (_ptDA->isSymmetric()) 
      // scaling upper block by diagonal
      for (int j = nsing0; j < nsing; j++) {
	for (int i = 0; i < nsing0; i++) {
	  DSsingCoefs0(i, j) *= DSsingCoefs0(i, i);
	}
      }
      blas_gemm<T>(CblasNoTrans, CblasNoTrans,
		   nsing1, nsing1, nsing0, 
		   _none, // alpha
		   DSsingCoefs0.addrCoefs() + nsing0, nsing, 
		   DSsingCoefs0.addrCoefs() + (nsing0 * nsing), nsing, 
		   _one, // beta
		   DSsingCoefs1.addrCoefs(), nsing1);
      if (_ptDA->isSymmetric()) {         // symmetrize
	for (int j = 0; j < nsing1; j++) {
	  for (int i = 0; i < j; i++) {
	    //	    DSsingCoefs1[j + i * nsing1] = DSsingCoefs1[i + j * nsing1];
	    DSsingCoefs1(j, i) = DSsingCoefs1(i, j);
	  }
	} 
      }
      diss_printf(_verbose, _fp,
		  "%s %d : DSsingCoefs1[]\n",
		  __FILE__, __LINE__);
      for (int i = 0; i < nsing1; i++) {
	for (int j = 0; j < nsing1; j++) {
	  diss_printf(_verbose,_fp, "%s ",
		      tostring<T>(DSsingCoefs1(i, j)).c_str());
	}
	diss_printf(_verbose, _fp, "\n");
      }
    }  //   if (nsing0 <= 0)
    bool flag;
    const U eps_machine = machine_epsilon<U, Z>();
    DSsingCoefs3.copy(DSsingCoefs2);
    int n0n, n0t;
    flag = ComputeDimKernel<T, U>(&n0n, &flag_unsym_permute,
				  DSsingCoefs1.addrCoefs(), nsing1,
				  _ptDA->isSymmetric(),
				  dim_augkern,
				  eps_machine,
				  eps_piv,
				  _verbose,
				  _fp);

    flag = ComputeDimKernel<T, U>(&n0t, &flag_unsym_permute,
				  DSsingCoefs3.addrCoefs(), nsing1,
				  _ptDA->isSymmetric(),
				  dim_augkern,
				  eps_machine,
				  eps_piv,
				  _verbose,
				  _fp);
    
    n0 = n0t > n0n ? n0t : n0n;
    // n0 = n0n;
    if (flag == false) {
      fprintf(stderr,
	      "%s %d : ERROR: kernel detection failed!\n",
	      __FILE__, __LINE__);
      //exit(-1);
    }
    if (_assume_invertible) {
      n0 = 0;
    }
  } // if (nn > 0)

  delete [] permute0;
  get_realtime(&t3_elapsed);
  //  clock_gettime(CLOCK_REALTIME, &ts3);
  diss_printf(_verbose, _fp, 
	      "%s %d : detection of the _dim. of the kernel (sec.) = %.6f\n",
	      __FILE__, __LINE__,
	      convert_time(t3_elapsed, t2_elapsed));
  } // if (enableDetection)

  // LDLt factorization of DS0(i,j)
  vector<int> permute, permute_left, permute_tright, permute_tleft;
  permute.resize(nsing);


  const double machine_eps_double = machine_epsilon<double, double>();
  const int n1 = nsing - n0;
  
  //  flag_unsym_permute = false;
  Schur.setFullPivoting(flag_unsym_permute);

  if (flag_unsym_permute) {
    int nn0;
    double fop;
    double last_pivot = 1.0;
    permute_left.resize(nsing);
    ldu_full_permute<T, U>(&nn0, n0, nsing, 
			   DSsingCoefs.addrCoefs(), nsing, &last_pivot,
			   &permute[0], &permute_left[0],
			   machine_eps_double, &fop);
#if 0
    last_pivot = 1.0;
    permute_tright.resize(nsing);
    permute_tleft.resize(nsing);
    ldu_full_permute<T, U>(&nn0, n0, nsing, 
			   DSsingCoefs2.addrCoefs(), nsing, &last_pivot,
			   &permute_tright[0], &permute_tleft[0],
			   machine_eps_double, &fop);
#endif
    n0 = nn0;
  }
  else {
    int nn0;
    double fop;
    double last_pivot = 1.0;
    if (_ptDA->isSymmetric()) {
      full_ldlt_permute<T, U>(&nn0, n0, nsing,
			      DSsingCoefs.addrCoefs(), nsing, &last_pivot,
			       &permute[0], machine_eps_double, &fop);
    }
    else {
      full_ldu_permute<T, U>(&nn0, n0, nsing,
			     DSsingCoefs.addrCoefs(), nsing, &last_pivot,
			     &permute[0], machine_eps_double, &fop);
    }
    n0 = nn0;
    permute_left = permute;
  }
  
  if (n0 != (nsing - n1)) {
    diss_printf(_verbose, _fp,
		"%s %d : factroization of suspicious pivots fails : %d -> %d\n",
		__FILE__, __LINE__, (nsing - n1), n0);
  }
  else {
    diss_printf(_verbose, _fp,
		"%s %d : pivot = %s after n0 = %d\n",
		__FILE__, __LINE__,
		flag_unsym_permute ? "full" : "symmetric", n0);
  } 
// keep offdiagonal part corresponding to the last Schur complement
  Schur.getScol().init(_dim, n1);
  for (int j = 0; j < n1; j++) {
    blas_copy<T>(_dim, DBsCoefs.addrCoefs() + (permute[j] * _dim), 1,
		 Schur.getScol().addrCoefs() + (j * _dim), 1);
  }
  if (n0 > 0) {
    Schur.getSchur().init(n1, n1); // copy with permutation for debugging 
    for (int j = 0; j < n1; j++) {
      for (int i = 0; i < n1; i++) {
	Schur.getSchur()(i, j) = DS0(permute_left[i], permute[j]);
      }
    }
  // d23 = (S_22)^{-1} (-S_23)
    ColumnMatrix<T> d23(nsing, n0);
    ColumnMatrix<T> d32(nsing, n0);

    d23.ZeroClear();   
    for (int j = 0; j < n0; j++) {
      for (int i = 0; i < n1; i++) {
	d23(i, j) = DS0(permute_left[i], permute[j + n1]);
      }
    }
    d32.ZeroClear();
    for (int j = 0; j < n0; j++) {
      for (int i = 0; i < n1; i++) {
	d32(i, j) = DS0(permute_left[j + n1], permute[i]);
      //d32(i, j) = DS2(permute_tleft[i], permute_tright[j + n1]);
      }
    }
    //    alpha = 1.0;
    blas_trsm<T>(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
		 n1, n0, _one, DSsingCoefs.addrCoefs(), nsing,
		 d23.addrCoefs(), nsing);
    {
      for (int i = 0; i < n1; i++) {
	const T stmp = DSsingCoefs(i, i);
	for (int j = 0; j < n0; j++) {
	  d23(i, j) *= stmp;
	}
      }
    }
    blas_trsm<T>(CblasLeft, CblasUpper, CblasNoTrans, CblasUnit,
		 n1, n0, _one, DSsingCoefs.addrCoefs(), nsing,
		 d23.addrCoefs(), nsing);
    for (int i = 0; i < n0; i++) {
      const int ii = i + n1;
      d23(ii, i) = _none;
    }
    blas_trsm<T>(CblasLeft, CblasUpper, CblasTrans, CblasUnit,
    		 n1, n0, _one, DSsingCoefs.addrCoefs(), nsing,
		 //blas_trsm<T>(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
		 //n1, n0, _one, DSsingCoefs2.addrCoefs(), nsing,
		 d32.addrCoefs(), nsing);
    {
      for (int i = 0; i < n1; i++) {
	const T stmp = DSsingCoefs(i, i);
	//const T stmp = DSsingCoefs2(i, i);
	for (int j = 0; j < n0; j++) {
	  d32(i, j) *= stmp;
	}
      }
    }
    blas_trsm<T>(CblasLeft, CblasLower, CblasTrans, CblasUnit,
    		 n1, n0, _one, DSsingCoefs.addrCoefs(), nsing,
		 // blas_trsm<T>(CblasLeft, CblasUpper, CblasNoTrans, CblasUnit,
		 // n1, n0, _one, DSsingCoefs2.addrCoefss(), nsing,
		 d32.addrCoefs(), nsing);
    
    for (int i = 0; i < n0; i++) {
      const int ii = i + n1;
      d32(ii, i) = _none;
    }
  // generating the kernel vectors
  // keep the vectors from suspicious pivots
    kernel.getKernBasis().free(); // the 2nd allocation _dim * nsing -> _dim * n
    kernel.getKernBasis().init(_dim, n0);
    kernel.getTKernBasis().free(); //
    kernel.getTKernBasis().init(_dim, n0);
    T *kern_basis = kernel.getKernBasis().addrCoefs();
    T *kernt_basis = kernel.getTKernBasis().addrCoefs();
    ColumnMatrix<T> d23_permuted(nsing, n0);
#if 1
      ColumnMatrix<T> dtest(nsing, n0);
      // for debugging
      d23_permuted.ZeroClear();
      for (int j = 0; j < n0; j++) {
	for (int i = 0; i < nsing; i++) {
	  d23_permuted(permute[i], j) = d23(i, j);
	}
      }
      blas_gemm<T>(CblasNoTrans, CblasNoTrans,
		   nsing, n0, nsing, _one,
		   DS0.addrCoefs(), nsing,
		   d23_permuted.addrCoefs(), nsing, _zero,
		   dtest.addrCoefs(), nsing);
      diss_printf(_verbose, _fp, "%s %d verify kernel : %d\n",
		  __FILE__, __LINE__, n0);
      for (int j = 0; j < n0; j++) {
	U stmp = blas_l2norm<T, U>(nsing, dtest.addrCoefs() + j * nsing, 1);
	for (int i = 0; i < nsing; i++) {
	  diss_printf(_verbose, _fp, "%d : %d %s %s\n",
		      i, permute[i],
		      tostring<T>(d23_permuted(i, j)).c_str(),
		      tostring<T>(dtest[i + j * nsing]).c_str());
	}
	diss_printf(_verbose, _fp, "%d %s\n", j, tostring<U>(stmp).c_str());
      }
      if (!_ptDA->isSymmetric()) {  
	d23_permuted.ZeroClear();
	for (int j = 0; j < n0; j++) {
	  for (int i = 0; i < nsing; i++) {
	    d23_permuted(permute_left[i], j) = d32(i, j);
	    //d23_permuted(permute_tright[i], j) = d32(i, j);
	  }
	}
	blas_gemm<T>(CblasTrans, CblasNoTrans,
		     nsing, n0, nsing, _one,
		     DS0.addrCoefs(), nsing,
		     //DS2.addrCoefs(), nsing, 
		     d23_permuted.addrCoefs(), nsing, _zero,
		     dtest.addrCoefs(), nsing);
	diss_printf(_verbose, _fp, "%s %d verify transposed kernel : %d\n",
		    __FILE__, __LINE__, n0);
	for (int j = 0; j < n0; j++) {
	  U stmp = blas_l2norm<T, U>(nsing, dtest.addrCoefs() + j * nsing, 1);
	  for (int i = 0; i < nsing; i++) {
	    diss_printf(_verbose, _fp, "%d : %d %s %s\n",
			i, permute_left[i],
			//i, permute_tright[i],
			tostring<T>(d23_permuted(i, j)).c_str(),
			tostring<T>(dtest[i + j * nsing]).c_str());
	  }
	  diss_printf(_verbose, _fp, "%d %s\n", j, tostring<U>(stmp).c_str());
	}
      }
#endif
      //
   d23_permuted.ZeroClear();
   if (flag_unsym_permute) {
      for (int j = 0; j < n0; j++) {
	for (int i = 0; i < nsing; i++) {
	  d23_permuted(permute[i], j) = d23(i, j);
	}
      }
      //  A11^-1 [ A_13 A_12] permute_right [S22^-1 S_23] - p_r [S22^-1 S_23]
      //                                    [     -I    ]       [     -I    ]
      kernel.getKernBasis().ZeroClear();
      for (int i = 0; i < nsing; i++) {
	DBsCoefs(singIdx[i], i) = _zero;
      }

      blas_gemm<T>(CblasNoTrans, CblasNoTrans,
		   _dim, n0, nsing, _one, // alpha
		   DBsCoefs.addrCoefs(), _dim, //
		   d23_permuted.addrCoefs(), nsing, _zero, // beta
		   kern_basis, _dim);
      for (int j = 0; j < n0; j++) {
	for (int i = 0; i < nsing; i++) {
	  kernel.getKernBasis()(singIdx[i], j) = -d23_permuted(i, j);
	}
      }
      // compute transposed kernel
      d23_permuted.ZeroClear();
      for (int j = 0; j < n0; j++) {
	for (int i = 0; i < nsing; i++) {
	  d23_permuted(permute_left[i], j) = d32(i, j);
	  //d23_permuted(permute_tright[i], j) = d32(i, j);
	}
      }
      // A11^-T [ A_31^T A_12^T] permute_l [S22^-T S_32^T] - p_l [S22^-T S_23^T]
      //                                   [     -I      ]       [     -I      ]
      kernel.getTKernBasis().ZeroClear();
      for (int i = 0; i < nsing; i++) {
	DCsCoefs(singIdx[i], i) = _zero;
      }
      blas_gemm<T>(CblasNoTrans, CblasNoTrans,
		   _dim, n0, nsing, _one, // alpha
		   DCsCoefs.addrCoefs(), _dim, //
		   d23_permuted.addrCoefs(), nsing, _zero, // beta
		   kernt_basis, _dim);
      for (int j = 0; j < n0; j++) {
	for (int i = 0; i < nsing; i++) {
	  kernel.getTKernBasis()(singIdx[i], j) = -d23_permuted(i, j);
	}
      }
      //      d23_permuted.free();
   }
   else { //  if (flag_unsym_permute) 
     for (int i = 0; i < nsing; i++) {
       DBsCoefs(singIdx[i], i) = _none;
      }
      for (int j = 0; j < n0; j++) {
	for (int i = 0; i < _dim; i++) {
	  kernel.getKernBasis()(i, j) = DBsCoefs(i, permute[j + n1]);
	}
      }
      // update the kernel using singular part of the Schur complement, d23
      // kern_basis = - kern_basis - scol * d23
      //              - [ A11^-1 A_13] + [A_11^-1 A_12] [S22^-1 S_23]
      //                [      0     ]   [     -I     ]
      //                [     -I     ]   [      0     ]
      blas_gemm<T>(CblasNoTrans, CblasNoTrans,
		   _dim, n0, n1, _one, // alpha
		   Schur.getScol().addrCoefs(), _dim, // none
		   d23.addrCoefs(), nsing, _none, // beta
		   kern_basis, _dim);

      // kernel of the transposed matrix
      ColumnMatrix<T> DCsCoefs12(_dim, n1);
      for (int i = 0; i < nsing; i++) {
	DCsCoefs(singIdx[i], i) = _none;
      }
      for (int j = 0; j < n1; j++) {
	blas_copy<T>(_dim, DCsCoefs.addrCoefs() + (permute[j] * _dim), 1,
		     DCsCoefs12.addrCoefs() + (j * _dim), 1);
      }
      for (int j = 0; j < n0; j++) {
	for (int i = 0; i < _dim; i++) {
	  kernel.getTKernBasis()(i, j) = DCsCoefs(i, permute[j + n1]);
	}
      }
      blas_gemm<T>(CblasNoTrans, CblasNoTrans,
		   _dim, n0, n1, _one, // alpha
		   DCsCoefs12.addrCoefs(), _dim, // none
		   d32.addrCoefs(), nsing, _none, // beta
		   kernt_basis, _dim);
      DCsCoefs12.free();
   } // if (flag_unsym_permute)
    // normalize each kernel_basis
    for (int j = 0; j < n0; j++) {
      U stmp(0.0);
      const U one(1.0);
      stmp = one / blas_l2norm<T, U>(_dim, (kern_basis + (j * _dim)), 1);
      blas_scal2<T, U>(_dim, stmp, (kern_basis + (j * _dim)), 1);
    }
    diss_printf(_verbose, _fp,
		"%s %d : == ptDA * kernel\n", __FILE__, __LINE__);

    for (int j = 0; j < n0; j++) {
      _ptDA->prod((kern_basis + (j * _dim)), v.addrCoefs());

      double norm_l2, norm_infty;
      calc_relative_norm<T>(&norm_l2, &norm_infty,
			    v.addrCoefs(), (kern_basis + (j * _dim)), _dim);
      diss_printf(_verbose, _fp,
		  "%d -th kernel scaled : norm_l2 = %.6e, norm_infty = %.6e / ",
		  j, norm_l2, norm_infty);
      if(_scaling) {
	calc_relative_normscaled<T, U>(&norm_l2, &norm_infty,
				       v.addrCoefs(), kern_basis + (j * _dim),
				       &_precDiag[0], _dim);
	diss_printf(_verbose, _fp,
		    "original : norm_l2 = %.6e, norm_infty = %.6e\n",
		    norm_l2, norm_infty);
      }
      else {
	diss_printf(_verbose, _fp, "\n");
      }// if (_scaling)
      for (int i = 0; i < nsing; i++) {
	diss_printf(_verbose, _fp, "%d %d %s\n",
		    i, singIdx[i],
		    tostring<T>(v(singIdx[i], j)).c_str());
      }
    }
    if (!_ptDA->isSymmetric()) {  
      for (int j = 0; j < n0; j++) {
	U stmp(0.0);
	const U one(1.0);
	stmp = one / blas_l2norm<T, U>(_dim, (kernt_basis + (j * _dim)), 1);
	blas_scal2<T, U>(_dim, stmp, (kernt_basis + (j * _dim)), 1);
      }
      diss_printf(_verbose, _fp,
		  "%s %d : == ptDA^T * transposed kernel\n",
		  __FILE__, __LINE__);
	
      for (int j = 0; j < n0; j++) {
	_ptDA->prodt((kernt_basis + (j * _dim)), v.addrCoefs());

	double norm_l2, norm_infty;
	calc_relative_norm<T>(&norm_l2, &norm_infty,
			      v.addrCoefs(), (kernt_basis + (j * _dim)), _dim);
	diss_printf(_verbose, _fp,
		    "%d -th  scaled : norm_l2 = %.6e, norm_infty = %.6e / ",
		    j, norm_l2, norm_infty);
	if(_scaling) {
	  calc_relative_normscaled<T, U>(&norm_l2, &norm_infty,
					 v.addrCoefs(),
					 kernt_basis + (j * _dim),
					 &_precDiag[0], _dim);
	  diss_printf(_verbose, _fp,
		      "original : norm_l2 = %.6e, norm_infty = %.6e\n",
		      norm_l2, norm_infty);
	}
	else {
	  diss_printf(_verbose, _fp, "\n");
	}
	diss_printf(_verbose, _fp, "%s %d\n", __FILE__, __LINE__);
  	for (int i = 0; i < nsing; i++) {
	  diss_printf(_verbose, _fp, "%d %d %s\n",
		      i, singIdx[i],
		      tostring<T>(v(singIdx[i], j)).c_str());
	}
      }
    } //  if (!_ptDA->isSymmetric()) 
    if(_scaling) {
      for (int j = 0; j < n0; j++) {
	const int jtmp = j * _dim;
	for (int i = 0; i < _dim; i++) {
	  kern_basis[i + jtmp] *= _precDiag[i];
	}
      }
      for (int j = 0; j < n0; j++) {
	const int jtmp = j * _dim;
	for (int i = 0; i < _dim; i++) {
	  kernt_basis[i + jtmp] *= _precDiag[i];
	}
      }
    }
#ifdef NORMALIZE_KERNEL_BASIS
  // normalize each kernel_basis
    for (int j = 0; j < n0; j++) {
      U stmp(0.0);
      const U one(1.0);
      stmp = one / blas_l2norm<T, U>(_dim, (kern_basis + (j * _dim)), 1);
      blas_scal2<T, U>(_dim, stmp, (kern_basis + (j * _dim)), 1);
      stmp = one / blas_l2norm<T, U>(_dim, (kernt_basis + (j * _dim)), 1);
      blas_scal2<T, U>(_dim, stmp, (kernt_basis + (j * _dim)), 1);
    
    }
#endif
    kernel.getKernProj().free();     // second alloctation : nsing -> n0
    kernel.set_dimension(n0);
    kernel.getKernProj().init(n0);
    kernel.getTKernProj().free();     // second alloctation : nsing -> n0
    if (!_ptDA->isSymmetric()) {  
      kernel.getTKernProj().init(n0);
      kernel.getNTKernProj().init(n0);
    }
    for (int i = 0; i < n0; i++) {
      for(int j = i; j < n0; j++) {
	kernel.getKernProj()(i, j) = blas_dot<T>(_dim,
						 (kern_basis + (i * _dim)), 1,
						 (kern_basis + (j * _dim)), 1);
      }
      // symmetrize
      for (int j = (i + 1); j < n0; j++) {
	kernel.getKernProj()(j, i) = kernel.getKernProj()(i, j);
      }
    }
    full_ldlh<T>(n0, kernel.getKernProj().addrCoefs(), n0);
    // inverse of diagonal part is also storead in the factorized matrix
    if (!_ptDA->isSymmetric()) {      
      for (int i = 0; i < n0; i++) {
	for(int j = i; j < n0; j++) {
	  kernel.getTKernProj()(i, j) = 
	    blas_dot<T>(_dim, 
			(kernt_basis + (i * _dim)), 1,
			(kernt_basis + (j * _dim)), 1);
	}
	// symmetrize
	for (int j = (i + 1); j < n0; j++) {
	  kernel.getTKernProj()(j, i) = kernel.getTKernProj()(i, j);
	}
      }
      full_ldlh<T>(n0, kernel.getTKernProj().addrCoefs(), n0);
      // for oblique projection
      for (int i = 0; i < n0; i++) {
	for(int j = 0; j < n0; j++) {
	  kernel.getNTKernProj()(i, j) = 
	    blas_dot<T>(_dim, 
			(kernt_basis + (i * _dim)), 1,
			(kern_basis + (j * _dim)), 1);
	}
      }
      diss_printf(_verbose, _fp,
		  "%s %d : orthogonality matrix kernels of A and A^T\n",
		  __FILE__, __LINE__);
      for (int i = 0; i < n0; i++) {
	diss_printf(_verbose, _fp, "%d : ", i);
	for(int j = 0; j < n0; j++) {
	  diss_printf(_verbose, _fp, "%.16e ",
		      blas_abs<T, double>(kernel.getNTKernProj()(i, j)));
	}
	diss_printf(_verbose, _fp, "\n");
      }
      full_ldu<T>(n0, kernel.getNTKernProj().addrCoefs(), n0);
    }
  } // if (n0 > 0)
  else {
    kernel.set_dimension(0);
  }
  //
  diss_printf(_verbose, _fp,
	    "%s %d : nsing = %d dim_augkern = %d n1 = %d n0 = %d\n",
	    __FILE__, __LINE__, nsing, dim_augkern, n1, n0);

  // for some applications which need index of singular entries
  kernel.setFullPivoting(flag_unsym_permute);
  kernel.getKernListEq().resize(n0);
  kernel.getKernListEqLeft().resize(n0);
  for (int i = 0; i < n0; i++) {
    kernel.getKernListEq()[i]= singIdx[permute[i + n1]];
    kernel.getKernListEqLeft()[i]= singIdx[permute_left[i + n1]];
  }
  diss_printf(_verbose, _fp, "%s %d : kern_list_eq[] = ", __FILE__, __LINE__);
  for (int i = 0; i < kernel.getKernListEq().size(); i++) {
    diss_printf(_verbose,_fp, "%d ", kernel.getKernListEq()[i]);
  }
  diss_printf(_verbose, _fp, "\n");
  if (flag_unsym_permute) {
    diss_printf(_verbose, _fp, "%s %d : kern_list_eq_left[] = ",
		__FILE__, __LINE__);
    for (int i = 0; i < kernel.getKernListEqLeft().size(); i++) {
      diss_printf(_verbose,_fp, "%d ", kernel.getKernListEqLeft()[i]);
    }
    diss_printf(_verbose, _fp, "\n");
  }
  if (n1 == dim_augkern && (!flag_unsym_permute)) {
    diss_printf(_verbose, _fp,
		"%s %d dimension of detected kernel == suspicous one = %d\n",
		__FILE__, __LINE__, n0);
    {
      int i = 0;
      typename vector<SquareBlockMatrix<T>* >::const_iterator kt = diags.begin();
      for (vector<vector<int> >::const_iterator it = augkern_indexes.begin();
	   it != augkern_indexes.end(); ++it, ++kt) {
	SquareBlockMatrix<T> &diag = *(*kt);
	const int n0_diag = diag.rank();
	for (vector<int>::const_iterator jt = (*it).begin(); jt != (*it).end();
	     ++jt, i++) {
	  diag.diag((*jt)) = save_diag[i]; // i++
	}
	diag.set_rank(n0_diag + (*it).size()); // restore
      }
    }
    Schur.getSlduList().resize(0);
  }
  else {  //  if (n1 == dim_augkern) {
    Schur.getSlduList().resize(n1);
    Schur.getSlduListLeft().resize(n1);
    for (int i = 0; i < n1; i++) {
      Schur.getSlduList()[i] = permute[i];
      Schur.getSlduListLeft()[i] = permute_left[i];
    }
    vector<int> s_list_eq, s_list_eq_left;
    s_list_eq.resize(n1);
    s_list_eq_left.resize(n1);
    // local2global index for the last Schur complement
    for (int i = 0; i < n1; i++) {
      s_list_eq[i] = singIdx[permute[i]];
      s_list_eq_left[i] = singIdx[permute_left[i]];
    }
    diss_printf(_verbose, _fp, "%s %d : s_list_eq[] = ",  __FILE__, __LINE__);
    for (int i = 0; i < s_list_eq.size(); i++) {
      diss_printf(_verbose, _fp, "%d ", s_list_eq[i]);
    }
    diss_printf(_verbose, _fp, "\n");
    if (flag_unsym_permute) {
      diss_printf(_verbose, _fp, "%s %d : s_list_eq_left[] = ",
		  __FILE__, __LINE__);
      for (int i = 0; i < s_list_eq_left.size(); i++) {
	diss_printf(_verbose, _fp, "%d ", s_list_eq_left[i]);
      }
      diss_printf(_verbose, _fp, "\n");
    }
    Schur.getSldu().init(flag_unsym_permute, s_list_eq, s_list_eq_left);
    
  // copy regular part n1*n1 from nsing*nsing matrix DSsingCoefs[]
    for (int j = 0; j < n1; j++) {
      blas_copy<T>(n1, DSsingCoefs.addrCoefs() + j * nsing, 1,
		   Schur.getSldu().addrCoefs() + j * n1, 1);
    }
    //  <---- n1 ---->   : given by pemrute[j]
    // [ A_11^-1 A_12 ] -> [ A_11^-1 A_12  ]
    // [  -I          ]    [   0           ] : singIdx[permute[i]] 
    // [   0          ]    [   0           ] :  0<= i < nsing
    for (int i = 0; i < n1; i++) {
	Schur.getScol()(s_list_eq[i], i) = _zero; 
    }
  } //  if (n1 == dim_augkern) {

  //  delete ptDA_acol; // 28 Dec.2015
  //  delete [] save_diag;
  //  delete [] DBsCoefs;
  //  delete [] DBsCoefs0;
  //  delete [] DSsingCoefs;
  //  delete [] DSsingCoefs0;
  //  delete [] DS0;
  //  delete [] DS2;
}

template
void DissectionSolver<double>::
BuildKernelsDetection(int &n0,
		      vector<int> &singIdx, 
		      vector<vector<int> > &augkern_indexes, 
		      vector<SquareBlockMatrix<double>* >& diags, 
		      const double eps_piv,
		      const int dim_augkern,
		      SchurMatrix<double> &Schur,
		      KernelMatrix<double> &kernel,
		      const bool enableDetection);

template
void DissectionSolver<complex<double>, double>::
BuildKernelsDetection(int &n0,
		      vector<int> &singIdx, 
		      vector<vector<int> >& augkern_indexes, 
		      vector<SquareBlockMatrix<complex<double> >* >& diags, 
		      const double eps_piv,
		      const int dim_augkern,
		      SchurMatrix<complex<double> > &Schur,
		      KernelMatrix<complex<double> > &kernel,
		      const bool enableDetection);

template
void DissectionSolver<quadruple>::
BuildKernelsDetection(int &n0,
		      vector<int> &singIdx, 
		      vector<vector<int> > &augkern_indexes, 
		      vector<SquareBlockMatrix<quadruple>* >& diags, 
		      const double eps_piv,
		      const int dim_augkern,
		      SchurMatrix<quadruple> &Schur,
		      KernelMatrix<quadruple> &kernel,
		      const bool enableDetection);

template
void DissectionSolver<complex<quadruple>, quadruple>::
BuildKernelsDetection(int &n0,
		      vector<int> &singIdx, 
		      vector<vector<int> >& augkern_indexes, 
		      vector<SquareBlockMatrix<complex<quadruple> >* >& diags, 
		      const double eps_piv,
		      const int dim_augkern,
		      SchurMatrix<complex<quadruple> > &Schur,
		      KernelMatrix<complex<quadruple> > &kernel,
		      const bool enableDetection);
template
void DissectionSolver<float>::
BuildKernelsDetection(int &n0,
		      vector<int> &singIdx, 
		      vector<vector<int> > &augkern_indexes, 
		      vector<SquareBlockMatrix<float>* >& diags, 
		      const double eps_piv,
		      const int dim_augkern,
		      SchurMatrix<float> &Schur,
		      KernelMatrix<float> &kernel,
		      const bool enableDetection);

template
void DissectionSolver<complex<float>, float>::
BuildKernelsDetection(int &n0,
		      vector<int> &singIdx, 
		      vector<vector<int> >& augkern_indexes, 
		      vector<SquareBlockMatrix<complex<float> >* >& diags, 
		      const double eps_piv,
		      const int dim_augkern,
		      SchurMatrix<complex<float> > &Schur,
		      KernelMatrix<complex<float> > &kernel,
		      const bool enableDetection);

//

template<typename T, typename U, typename W, typename Z>
int DissectionSolver<T, U, W, Z>::
kern_dimension(void) 
{ 
  int itmp = 0;
  for (int m = 0; m < _graph_colors; m++) {
    itmp += _kernel[m].dimension();
  }
  return itmp; 
}

//
template
int DissectionSolver<double>::kern_dimension(void);

template
int DissectionSolver<quadruple>::kern_dimension(void);

template
int DissectionSolver<complex<double>, double>::kern_dimension(void);

template
int DissectionSolver<complex<quadruple>, quadruple>::kern_dimension(void);

template
int DissectionSolver<double, double, quadruple, quadruple>::kern_dimension(void);

template
int DissectionSolver<quadruple, quadruple, double, double>::kern_dimension(void);

template
int DissectionSolver<float>::kern_dimension(void);

template
int DissectionSolver<complex<float>, float>::kern_dimension(void);

//

//
template<typename T, typename U, typename W, typename Z>
int DissectionSolver<T, U, W, Z>::
postponed_pivots(void) 
{ 
  int itmp = 0;
  for (int m = 0; m < _graph_colors; m++) {
    itmp += _singIdx[m].size();
  }
  return itmp; 
}

//
template
int DissectionSolver<double>::postponed_pivots(void);

template
int DissectionSolver<quadruple>::postponed_pivots(void);

template
int DissectionSolver<complex<double>, double>::postponed_pivots(void);

template
int DissectionSolver<complex<quadruple>, quadruple>::postponed_pivots(void);

template
int DissectionSolver<double, double, quadruple, quadruple>::postponed_pivots(void);

template
int DissectionSolver<quadruple, quadruple, double, double>::postponed_pivots(void);

template
int DissectionSolver<float>::postponed_pivots(void);

template
int DissectionSolver<complex<float>, float>::postponed_pivots(void);

//

template<typename T, typename U, typename W, typename Z>
int DissectionSolver<T, U, W, Z>::
ComputeTransposedKernels(void)
{
  //  int n0_max = 0;
  if (_ptDA->isSymmetric()) {
    diss_printf(_verbose, _fp,
	    "%s %d : ComputeTransposedKernels : matrix is symmetric\n",
	    __FILE__, __LINE__);
    return (-1);
  }

  for (int m = 0; m < _graph_colors; m++) {
    int n0 = _kernel[m].dimension();
#if 0  
    if (!_tridiagQueue[m]->tridiagSolver()) {
      _kernel[m].getTKernBasis().init(_dim, n0);
    }
#endif
    _kernel[m].getTKernProj().init(n0);
    _kernel[m].getNTKernProj().init(n0);
  }
  for (int m = 0; m < _graph_colors; m++) {
    int n0 = _kernel[m].dimension();
    T *kernt_basis = _kernel[m].getTKernBasis().addrCoefs();
    T *kernn_basis = _kernel[m].getKernBasis().addrCoefs();
#if 0
    if (!_tridiagQueue[m]->tridiagSolver()) {
      _kernel[m].getTKernBasis().ZeroClear();
      for (int i = 0; i < n0; i++) {
	const int ii = _kernel[m].getKernListEq()[i];
	for (int k = _ptDA->ptRow(ii); k < _ptDA->ptRow(ii + 1); k++) {
	  _kernel[m].getTKernBasis()(_ptDA->indCol(k), i) = _ptDA->Coef(k);
	}
      }
      for (int j = 0; j < n0; j++) {
	for (int i = 0; i < n0; i++) {
	  _kernel[m].getTKernBasis()(_kernel[m].getKernListEq()[i], j) = _zero;
	}
      }
      if (n0 > 1) {
	SolveMulti(kernt_basis, n0, false, true, false);
      }
      else {
	                 // porjection, isTrans, isScaling
	SolveSingle(kernt_basis, false, true, false); 
      }
      // nullify singular block 
      for (int j = 0; j < n0; j++) {
	for (int i = 0; i < n0; i++) {
	  _kernel[m].getTKernBasis()(_kernel[m].getKernListEqLeft()[i], j) = _zero;
	}
      }
      // set -I on the diagonal entries of singular block
      for (int i = 0; i < n0; i++) {
	_kernel[m].getTKernBasis()(_kernel[m].getKernListEqLeft()[i], i) = _none;
      }
      if(_scaling) {
	for (int j = 0; j < n0; j++) {
	  const int jtmp = j * _dim;
	  for (int i = 0; i < _dim; i++) {
	    kernt_basis[i + jtmp] *= _precDiag[i];
	  }
	}
      }
    } //  if (_tridiagQueue[m]->tridiagSolver())
    else {
      diss_printf(_verbose, _fp,
		  "%s %d : transposed kernel from tridiag_kernelt_basis()\n",
		  __FILE__, __LINE__);
    }
#endif
    VectorArray<T> vn(_dim);
    VectorArray<T> vt(_dim);
    diss_printf(_verbose, _fp,
		"%s %d : check kern_basis belonging to the kernel of A\n",
	      __FILE__, __LINE__);
    for (int j = 0; j < n0; j++) {
      SpMV(&kernn_basis[j * _dim], vn.addrCoefs()); //&vn[j][0]);
      double norm_l2, norm_infty;
      calc_relative_norm<T>(&norm_l2, &norm_infty, vn.addrCoefs(), //&vn[j][0],
			    &kernn_basis[j * _dim], _dim);
      diss_printf(_verbose, _fp, "%s %d : %d-th : l2_norm = %.6e, infty_norm = %.6e\n",
		__FILE__, __LINE__, j, norm_l2, norm_infty);
    }
    diss_printf(_verbose, _fp,
		"%s %d : check kern_basis belonging to the kernel of A^T\n",
		__FILE__, __LINE__);
    for (int j = 0; j < n0; j++) {
      SpMtV(&kernt_basis[j * _dim], vt.addrCoefs()); 
      double norm_l2, norm_infty;
      calc_relative_norm<T>(&norm_l2, &norm_infty, vt.addrCoefs(),  //&vt[j][0],
			    &kernt_basis[j * _dim], _dim);
      diss_printf(_verbose, _fp,
		  "%s %d : %d-th : l2_norm = %.6e, infty_norm = %.6e\n",
		  __FILE__, __LINE__, j, norm_l2, norm_infty);
    }  // loop : j
    vn.free();
    vt.free();
#ifdef NORMALIZE_KERNEL_BASIS2
      // normalize each kernel_basis
    for (int j = 0; j < n0; j++) {
      U stmp(0.0);
      const U one(1.0);
      stmp = blas_l2norm<T, U>(_dim, (kernt_basis + (j * _dim)), 1);
      stmp = one / stmp;
      blas_scal2<T, U>(_dim, stmp, (kernt_basis + (j * _dim)), 1);
    }
#endif
    for (int i = 0; i < n0; i++) {
      for(int j = i; j < n0; j++) {
	_kernel[m].getTKernProj()(i, j) = 
	  blas_dot<T>(_dim, 
		      (kernt_basis + (i * _dim)), 1,
		      (kernt_basis + (j * _dim)), 1);
      }
      // symmetrize
      for (int j = (i + 1); j < n0; j++) {
	_kernel[m].getTKernProj()(j, i) = _kernel[m].getTKernProj()(i, j);
      }
    }
    full_ldlh<T>(n0, _kernel[m].getTKernProj().addrCoefs(), n0);
// for oblique projection
    for (int i = 0; i < n0; i++) {
      for(int j = 0; j < n0; j++) {
	_kernel[m].getNTKernProj()(i, j) = 
	  blas_dot<T>(_dim, 
		      (kernt_basis + (i * _dim)), 1,
		      (kernn_basis + (j * _dim)), 1);
      }
    }
    diss_printf(_verbose, _fp,
	      "%s %d : orthogonality matrix kernels of A and A^T : %d\n",
		__FILE__, __LINE__, m);
    for (int i = 0; i < n0; i++) {
      diss_printf(_verbose, _fp, "%d : ", i);
      for(int j = 0; j < n0; j++) {
	diss_printf(_verbose, _fp, "%.16e ",
		    blas_abs<T, double>(_kernel[m].getNTKernProj()(i, j)));
      }
      diss_printf(_verbose, _fp, "\n");
    }
    full_ldu<T>(n0, _kernel[m].getNTKernProj().addrCoefs(), n0);
  } // loop : m
  return 1;
}

template
int DissectionSolver<double>::ComputeTransposedKernels(void);

template
int DissectionSolver<complex<double>, double>::ComputeTransposedKernels(void);


template
int DissectionSolver<quadruple>::ComputeTransposedKernels(void);

template
int DissectionSolver<complex<quadruple>, quadruple>::
ComputeTransposedKernels(void);

template
int DissectionSolver<float>::ComputeTransposedKernels(void);

template
int DissectionSolver<complex<float>, float>::ComputeTransposedKernels(void);

