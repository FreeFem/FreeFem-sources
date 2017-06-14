/*! \file SparseRenumbering.cpp
    \brief tridiagonal factorization algorithm with Cuthill-McKee
    \author François-Xavier Roux, ONERA, Laboratoire Jacques-Louis Lions
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
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

#include <cstdio>
#include <cstdlib>
#include <list>
#include "Algebra/SparseRenumbering.hpp"
#include "Driver/DissectionDefault.hpp"

using std::list;

void CMK_number(const int dim, const int *ptrows, const int *indcols,
		vector<int> &new2old, const bool verbose, FILE *fp)
{
  vector<int> list, indic, connect;
  vector<bool> mask;
  int i1_prev, i2_prev, i0, i1;
  double profile, min_profile;

  connect.resize(dim);
  list.resize(dim);
  indic.resize(dim);
  mask.resize(dim);
  
  // compute initial profile
  i1_prev = 0;
  i2_prev = indcols[ptrows[0]];
  for (int k = (ptrows[0] + 1); k < ptrows[1]; k++) {
    i2_prev = indcols[k] > i2_prev ? indcols[k] : i2_prev;
  }
  min_profile = ((double)(i2_prev - i1_prev + 2)
		 * (double)(i2_prev - i1_prev + 1)) * 0.5;
  while (i2_prev < (dim - 1)) {
    int i1 = i2_prev + 1;
    int i2 = indcols[ptrows[i1_prev]];
    for (int k = (ptrows[i1_prev] + 1); k < ptrows[i2_prev + 1]; k++) {
      i2 = indcols[k] > i2 ? indcols[k] : i2;
    }
    min_profile += (double)(i2 - i1 + 2) * (double)(i2 - i1 + 1) * 0.5;
    i1_prev = i1;
    i2_prev = i2;
  }
  if (verbose) {
    fprintf(fp, "%s %d : Cuthill McKee renumbering : init block : %g :: ",
	    __FILE__, __LINE__, min_profile);
  }

  // initialize frontal renumbering from an excentric point
  for (int i = 0; i < dim; i++) {
    connect[i] = ptrows[i + 1] - ptrows[i];
    mask[i] = false;
  }
  i0 = 0; // initialization : critical when connent[] == dim
  {
    int imin = dim;    // find index that attains minimun of connect[]
    for (int i = 0; i < dim; i++) {
      if (imin > connect[i]) {
	imin = connect[i];
	i0 = i;
      }
    }
  }
  mask[i0] = true;
  profile = frontal_numb(dim, ptrows, indcols, i0, list, indic, connect);
  i1 = list[dim - 1];

  if (verbose) {
    fprintf(fp, "- > node %d : %g :: ", i0, profile);
  }
  if (profile < min_profile) {
    min_profile = profile;
  }
  else {
    for (int i = 0; i < dim; i++) {
      list[i] = i;
    }
  }
  profile = min_profile;
  min_profile = min_profile + 1.0;
  
  // repeat frontal numbering as long as profile improves
  while (profile < min_profile) {
    min_profile = profile;
    for (int i = 0; i < dim; i++) {
      new2old[i] = list[i];
    }
    for (int i = 0; i < dim; i++) {
      connect[i] = ptrows[i + 1] - ptrows[i];
    }
    {
      i0 = 0; // initialization : critical when connent[] == dim
      int imin = dim;
      for (int k = ptrows[i1]; k < ptrows[i1 + 1]; k++) {
	const int ii = indcols[k];
	if (imin > connect[ii]) {
	  imin = connect[ii];
	  i0 = ii;
	}
      }
      if (mask[i0]) {
	for (int i = (dim - 1); i >=0; i--) {
	  if (!mask[list[i]]) {
	    i0 = list[i];
	    break;
	  }
	}
      }
    }
    mask[i0] = true;
    profile = frontal_numb(dim, ptrows, indcols, i0, list, indic, connect);
    i1 = list[dim - 1];
    if (verbose) {
      fprintf(fp, "- > node: %d : %g :: ", i0, profile);
    }
  } // while

  if (verbose) {
    fprintf(fp, "-> optimazied : %g\n", min_profile);
  }
}

double frontal_numb(const int dim, const int *ptrows, const int *indcols,
		    const int i0, 
		    vector<int> &list,
		    vector<int> &indic, vector<int> &connect)
{
  int n, nn, i1, i2, i1_prev, i2_prev, k1, k2;
  double profile;
  // compute nodal connectivity degree
  for (int i = 0; i < dim; i++) {
    connect[i] = ptrows[i + 1] - ptrows[i];
  }
#if 1
  for (int i = 0; i < dim; i++) { 
    indic[i] = 0;
  }
  n = 0;
  nn = (-1);
  list[0] = i0;
  indic[i0] = 1;
  while (n < (dim - 1)) {
    nn++;
    int j = list[nn];
    // renumber un-numbered neigbours of j, in increasing connectivity order
    int k = 0;
    for (int l = ptrows[j]; l < ptrows[j + 1]; l++) {
      const int i = indcols[l];
      if (indic[i] == 0) {
	int ki = 1;
	for ( ; ki <= k; ki++) {
	  if (connect[i] < connect[list[n + ki]]) {
	    break;
	  }
	}
	k++;
	for (int kk = k; kk >= (ki + 1); kk--) {
	  list[n + kk] = list[n + kk - 1];
	}
	list[n + ki] = i;
      } // if (indic[i] == 0)    
    }    // update connectivity degree of un-numbered nodes
    for (int kk = (n + 1); kk <= (n + k); kk++) {
      const int ii = list[kk];
      for (int m = ptrows[ii]; m < ptrows[ii + 1]; m++) {
	connect[indcols[m]]--;
      }
    }
    for (int kk = (n + 1); kk <= (n + k); kk++) {
      indic[list[kk]] = 1;
    }
    n += k;
  } // while (n < (dim - 1))
    // build old to new correspondance
  for (int i = 0; i < dim; i++) {
    indic[list[i]] = i;
  }
#else
  std::list<int> idexlist;
  vector<bool> mask;
  mask.resize(dim, false);
  idexlist.push_back(i0);
  std::list<int>::iterator it = idexlist.begin();
  mask[(*it)] = true;
  n = 0;
  std::list<int>::iterator nt = idexlist.begin();
  while (n < (dim - 1)) {
    int k = 0;
    for (int l = ptrows[(*it)]; l < ptrows[(*it) + 1]; l++) {
      const int ii = indcols[l];
      if (!mask[ii]) {
	bool flag = false;
	std::list<int>::iterator kt = nt;
	++kt;
	for (int kk = 0 ; kk < k; kk++, ++kt) {
	  if (connect[ii] < connect[(*kt)]) {
	    idexlist.insert(kt, ii);
	    flag = true;
	    break;
	  }
	}
	if (!flag) {
	  idexlist.push_back(ii);
	}
	k++;
      }
    }
    ++it;
    std::list<int>::iterator kt = nt;
    int kk = 0;
    for ( ; kk < k; ++kt, kk++) {
      for (int m = ptrows[(*kt)]; m < ptrows[(*kt) + 1]; m++) {
	connect[indcols[m]]--;
      }
      mask[(*kt)] = true;
      ++nt;
    }
    n += k;
  }
  {
    int i = 0;
    std::list<int>::iterator it = idexlist.begin();
    for (; it != idexlist.end(); ++it, i++) {
      list[i] = (*it);
      indic[(*it)] = i;
    }
  }
  idexlist.clear();
#endif
  
  // compute profile
  i1_prev = 0;
  k1 = ptrows[list[i1_prev]];
  k2 = ptrows[list[i1_prev] + 1];
  i2_prev = indic[indcols[k1]];
  for (int k = (k1 + 1); k < k2; k++) {
    i2_prev = indic[indcols[k]] > i2_prev ? indic[indcols[k]] : i2_prev;
  }
  profile = ((double)(i2_prev - i1_prev + 1) *
	     (double)(i2_prev - i1_prev + 2) * 0.5);
  while (i2_prev < (dim - 1)) {
    i1 = i2_prev + 1;
    i2 = i1;
    for (int i = i1_prev; i <= i2_prev; i++) {
      k1 = ptrows[list[i]];
      k2 = ptrows[list[i] + 1];
      i2 = std::max(i2, indic[indcols[k1]]);
      for (int k = (k1 + 1); k < k2; k++) {
	  i2 = indic[indcols[k]] > i2 ? indic[indcols[k]] : i2;
      }
    } 
    profile += ((double)(i2 - i1 + 1) * (double)(i2 - i1 + 2) * 0.5);
    i1_prev = i1;
    i2_prev = i2;
  }
  return profile;
}
    
int point_front(const int dim, const int *ptrows, const int *indcols,
		vector<int> &new2old, vector<int> &p_front)
{
  int nfront, ilast;
  vector<int> old2new, p_front1;
  vector<bool> mask;
  const int size_thrs = 24;
  old2new.resize(dim);
  for (int i = 0; i < dim; i++) {
    old2new[new2old[i]] = i;
  }
  p_front.resize(dim + 1); // for safty
  
  nfront = 0;
  p_front[0] = 0;
  ilast = (-1);
  while (ilast < (dim - 1)) { // loop ilast
    ilast++;
    bool flag_cont = false;
    bool flag_break = false;
    const int ii = new2old[ilast];
    nfront++;
    int isup = (-1);
    for (int k = ptrows[ii]; k < ptrows[ii + 1]; k++) {
      const int jj = old2new[indcols[k]];
      isup = jj > isup ? jj : isup;
    }
    p_front[nfront] = isup + 1;
    while (isup > ilast) {
      ilast = isup;
      for (int i = p_front[nfront - 1]; i < p_front[nfront]; i++) {
	const int ii = new2old[i];
	for (int k = ptrows[ii]; k < ptrows[ii + 1]; k++) {
	  const int jj = old2new[indcols[k]];
	  isup = jj > isup ? jj : isup;
	}
      } // loop : i
      if (isup == ilast) {
	flag_cont = true;
	break; // continue loop_ilast
      }
      nfront++;
      p_front[nfront] = isup + 1;
      if (isup == (dim - 1)) {
	flag_break = true;
	break; // exit loop_ilast
      }
    }
    if (flag_cont) {
      continue;
    }
    if (flag_break) {
      break;
    }
  } // loop ilast
  old2new.clear();
  mask.resize(nfront + 1);
  p_front1.resize(nfront + 1);
  for (int k = 0; k < nfront; k++) {
    const int itmp = p_front[k + 1] - p_front[k];
    if (itmp >= size_thrs) {
      mask[k] = true;
    }
    else {
      mask[k] = false;
    }
  }
  mask[nfront] = true;
  int n = 0;
  {
    int k = 0;
    p_front1[0] = p_front[0];
    int k0 = k;
    while (k < nfront) {
      int flag = 0;
      if (!mask[k0]) {
	if (!mask[k + 1]) {
	  if ((p_front[k + 1] - p_front[k0]) < size_thrs) {
	    flag = 2;
	  }
	  else {
	    flag = (-2);
	  }
	}
	else {
	  if ((k == 0) && mask[1]) {
	    flag = 2;
	  }
	  else {
	    if ((k0 == k) || (k == (nfront - 1))) {
	      flag = 0;
	    }
	    else {
	      flag = (-1);
	    }
	  }
	} //  if (!mask[k + 1]) 
      } // if (!mask[k0])
      else {
	if (k == (nfront - 2) && !mask[nfront - 1]) {
	  flag = 1;
	}	
      }
      switch (flag) {
      case 1:
	p_front1[n + 1] = p_front[k + 2];
	n++;
	k += 2;
	k0 = k;
	break;
      case (-2):
	p_front1[n + 1] = p_front[k];
	n++;
	k0 = k;
	break;
      case (-1):
	p_front1[n + 1] = p_front[k];
	p_front1[n + 2] = p_front[k + 1];
	n += 2;
	k++;
	k0 = k;
	break;
      case (0):
	p_front1[n + 1] = p_front[k + 1];
	n++;
	k++;
	k0 = k;
	break;
      case 2:
	k++;
	break;
      }
    } // while (k < nfront)
  } // scope for k, n etc.
  if ((p_front1[0] != 0)|| (p_front1[n] != p_front[nfront])) {
    fprintf(stderr, "%s %d : error %d -> %d\n", __FILE__, __LINE__,
	    p_front[nfront], p_front1[n]);
    exit(-1);
  }
  nfront = n;
  p_front.resize(nfront + 1);
  for (int k = 0; k <= n; k++) {
    p_front[k] = p_front1[k];
  }
  return nfront;
}

int getColorMaskCSR(int *color_mask, const CSR_indirect *csr, 
		    const bool verbose, FILE *fp)
{
  const int dim = csr->n;
  int *prow = csr->ptRows;
  int *indcols = csr->indCols;

  for (int i = 0; i < dim; i++) {
    color_mask[i] = (-1);
  }

  list<int> graph_start;
  int num_isolated = 0;
  for (int i = 0 ; i < dim; i++) {
    if ((prow[i] + 1) == prow[i + 1]) {
      num_isolated++;
      color_mask[i] = 0;
    }
  }
  if (verbose && (num_isolated > 0)) {
    fprintf(fp, "%s %d : isolated = %d\n", __FILE__, __LINE__, num_isolated);
  }
  int color = 0;
  for (int i = 0; i < dim; i++) {
    if (color_mask[i] == (-1)) {
      color++;
      color_mask[i] = color;
      graph_start.clear();
      graph_start.push_back(i);
      for (list<int>::iterator it = graph_start.begin(); 
	   it != graph_start.end(); ++it) {
	for (int k = prow[(*it)]; k < prow[(*it) + 1]; k++) {
	  if (color_mask[indcols[k]] == (-1)) {
	    color_mask[indcols[k]] = color;
	    graph_start.push_back(indcols[k]);
	  }
	}
      } // loop : it
    }
  } // loop : i
  // post process merge color with small size into previous color
  int reduced = 0;
  if (verbose) {
    fprintf(fp, "%s %d : before color = %d\n",
	    __FILE__, __LINE__, color);
  }
  int m = 1;
  while (m <= color) {
    int count = 0;
    for (int i = 0; i < dim; i++) {
      if (color_mask[i] == m) {
	count++;
      }
    }
    if ((count <= DIM_AUG_KERN) && (color > 1)){ // to treat very small matrix
      reduced++;
      color--;
      for (int i = 0; i < dim; i++) {
	if (color_mask[i] == m) {
	  color_mask[i] = (-1);
	}
	else if (color_mask[i] > m) {
	  color_mask[i]--;
	}
      }
    }
    else {
      m++;
    }
  }
  if (verbose) {
    fprintf(fp, "%s %d : after fused = %d color = %d\n",
	    __FILE__, __LINE__, reduced, color);
  }
  return color;
}

