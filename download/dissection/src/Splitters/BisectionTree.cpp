/*! \file   BisectionTree.cpp
    \brief  ordeing by graph decomposer, SCOTCH or METIS
    \author Xavier Juvigny, ONERA
    \date   Jul.  2nd 2012
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

// ==============================================================
// ==     Bisection tree : Implementation of some methods     ==
// ==============================================================
# include <cstdlib>
# include <cstring>
# include <stdexcept>
# include <algorithm>
# include <map>
# include <time.h>
# include "BisectionTree.hpp"
# include "ScotchSplitter.hpp"
# include "MetisSplitter.hpp"
# include "BitTools/BitManipulations.hpp"
# include "Compiler/DebugUtils.hpp"

#include <vector>
#include <list>

using namespace Dissection;
using std::vector;
using std::list;

unsigned*
Tree::compLoc2Dom(unsigned dim) const
{
  /** The splitters sort the domains to be compatible
      with the domain numerotation cared in this
      class...
  */
  unsigned nbDoms = (1<<_nbLevels) - 1;
  unsigned* loc2Dom = new unsigned[dim];
  for (unsigned iDom = 0; iDom < nbDoms; iDom++)
  {
    for (unsigned iNode = _ptOnDomains[iDom];
	 iNode<_ptOnDomains[iDom]+_sizeOfDomains[iDom];
	 iNode++) {
      assert(iNode<dim);
      loc2Dom[iNode] = iDom+1;	
    }
  }
  return loc2Dom;
}
// ==============================================================
void
Tree::renumberingInterface(unsigned nperm, unsigned* perm,
			   const SetOfStrips& paral, 
			   const SetOfStrips& parar,
			   const SetOfStrips& seq,
			   unsigned layer, unsigned iDom ) const
{
  if (nperm<3) return;
  long endPerm = long(nperm-1);
  unsigned endParal= 0;
  unsigned endParar= 0;
  unsigned endSeq  = 0;
  // Begin sorting perm indices which are in paral :
  if (paral.numberOfStrips() > 0) {
    unsigned iperm = 0;
    while (long(iperm) <= endPerm) {
      if (paral.contains(perm[iperm]))
	iperm++;
      else {
	std::swap(perm[iperm],perm[endPerm]);
	endPerm -= 1;
      }
    }
    endParal = iperm;
    if (endParal > 1)
      std::sort(perm,perm+endParal);
  }
  else endParal = 0;
  // Sorting sequential indices given
  if (seq.numberOfStrips() > 0) {
    unsigned iperm = endParal;
    endPerm        = nperm-1;
    while ((long)iperm <= endPerm) {
      if (seq.contains(perm[iperm]))
	++iperm;
      else {
	std::swap(perm[iperm],perm[endPerm]);
	--endPerm;
      }
    }
    endSeq = iperm;
    if (endSeq-endParal > 1)
      std::sort(perm+endParal,perm+endSeq);
  }
  else endSeq = endParal;
  // Sorting parallel indices from right interface :
  if (parar.numberOfStrips() > 0) {
    unsigned iperm = endSeq;
    endPerm        = nperm-1;
    while (long(iperm) <= endPerm) {
      if (parar.contains(perm[iperm]))
	iperm++;
      else {
	std::swap(perm[iperm],perm[endPerm]);
	endPerm -= 1;
      }
    }
    endParar = iperm;
    if (endParar-endSeq > 1)
      std::sort(perm+endSeq,perm+endParar);
  }
  else endParar = endSeq;
  if (endParar < nperm-1)
    std::sort(perm+endParar,perm+nperm);

  // Recursive calling :
  if (nodeLayer(iDom) < _nbLevels-2) {
    SetOfStrips paral2, parar2, seq2;
    iDom *= 2;
    unsigned leftChild  = (iDom<<1);
    unsigned rghtChild  = leftChild+1;
    SetOfStrips& lStrip = _interco[leftChild-1][layer];
    SetOfStrips& rStrip = _interco[rghtChild-1][layer];
    lStrip.para_union(rStrip, paral2, parar2,seq2);
    if (endParal > 0) {
      renumberingInterface(endParal, perm, paral2, parar2,
			   seq2, layer, iDom );	
    }
    // For sequential part, given priority to the left
    // interface if right interface has parallel indices
    if (endParar > endSeq) {
      renumberingInterface(endSeq-endParal, perm+endParal, 
			   paral2, parar2, seq2, layer, iDom );
      iDom += 1;
      unsigned leftChild  = (iDom<<1);
      unsigned rghtChild  = leftChild+1;
      SetOfStrips& lStrip = _interco[leftChild-1][layer];
      SetOfStrips& rStrip = _interco[rghtChild-1][layer];
      lStrip.para_union(rStrip, paral2, parar2,seq2);
      renumberingInterface(endParar-endSeq, perm+endSeq, 
			   paral2, parar2, seq2, layer, iDom );
    } else {
      iDom += 1;
      unsigned leftChild  = (iDom<<1);
      unsigned rghtChild  = leftChild+1;
      SetOfStrips& lStrip = _interco[leftChild-1][layer];
      SetOfStrips& rStrip = _interco[rghtChild-1][layer];
      lStrip.para_union(rStrip, paral2, parar2,seq2);	
      renumberingInterface(endSeq-endParal, perm+endParal, 
			   paral2, parar2, seq2, layer, iDom );
    }
  }
}
// ==============================================================
Tree::Tree(FILE *fp,
	   bool &berr, unsigned dim,
	   const CSR_indirect *csr,
	   const bool isSym, //const bool isLower,
	   const int *remap_eqn, //const int* map_indcols,
	   unsigned nbMaxLevel, unsigned minSize, 
	   splitter spltFct,
	   bool checkData, const bool verbose) :
  _verbose(verbose), _glob2dom(NULL), _isSym(isSym)
{
  const int *ptRows = csr->ptRows;
  const int *indCols = csr->indCols;
  const int *indVals = csr->indVals;
  const int *unsym_upper2lower = csr->indVals_unsym;
  bool ok;
  if (!spltFct) spltFct = &ScotchSplitter;
  // 0. Remove loops for diagonal elements (Metis and Scotch
  //    don't like loops into graph)
  int * ptRows2, * indCols2;
  int *indCols_sbdmn, *indCols_idxStrip, *indCols_idxSbdmn;
  unsigned *l2d;
  //  int *ptRow_diag; // 01 Feb 2013 : Atsushi
#if 0
  struct timespec tt0, tt1, tt2, tt3, tt4, tt5;
  clock_gettime(CLOCK_REALTIME, &tt0);
#endif
  //  TRACE("Removing loops\n");
  removeLoops(dim, ptRows, indCols, ptRows2, indCols2);
  // 1. Use splitter utility (Scotch, Metis or other tool)
  //    which convert CSR format into internal format
  //    (by example, removing diagonal part in the graph)
  //  TRACE("Call Scotch/Metis splitter\n");
  _local2global = new int[dim];
  _global2local = new int[dim];
#if 0
  clock_gettime(CLOCK_REALTIME, &tt1);
#endif
  int  nbDoms;
#if 0
  fich = fopen("splitter.data", "r");
  fscanf(fich, "# %d %d", &itmp, &nbDoms);
  if (itmp != dim) {
    fprintf(stderr, "splitter.data error : %d != %d\n", itmp, dim);
  }
  for (int i = 0 ; i < dim; i++) {
    fscanf(fich, "%d %d", &_global2local[i], &_local2global[i]);
  }
  for (int i = 0; i < (nbDoms + 1); i++) {
    fscanf(fich, "%d", &_ptOnDomains[i]);
  }
  for (int i = 0; i < nbDoms; i++) {
    fscanf(fich, "%d", &_sizeOfDomains[i]);
  }
  fclose(fich);
#else
  ok = (*spltFct)(dim, ptRows2,indCols2,nbMaxLevel,minSize,
		  _global2local,_local2global,
		  nbDoms,_ptOnDomains,_sizeOfDomains, 
		  checkData, verbose, fp);
#if 0
  fich = fopen("splitter.data", "w");
  fprintf(fich, "# %d %d\n", dim, nbDoms);
  for (int i = 0 ; i < dim; i++) {
    fprintf(fich, "%d %d\n", _global2local[i], _local2global[i]);
  }
  for (int i = 0; i < (nbDoms + 1); i++) {
    fprintf(fich, "%d\n", _ptOnDomains[i]);
  }
  for (int i = 0; i < nbDoms; i++) {
    fprintf(fich, "%d\n", _sizeOfDomains[i]);
  }
  fclose(fich);
#endif
#endif
#if 0
  if (!ok) {
    berr = false;
    std::cerr << __FILE__ << " : " << __LINE__ << std::endl;
    delete [] _sizeOfDomains;
    delete [] _ptOnDomains;
    delete [] _local2global;
    delete [] _global2local;
    if (_interco) delete [] _interco;
    return;
    // 4 Jul.2013 Atsushi
    throw std::runtime_error("Failed splitting matrix graph");
  }
  clock_gettime(CLOCK_REALTIME, &tt2);
#endif
  _nbLevels = nbMaxLevel;
  assert(nbDoms == (1<<_nbLevels)-1);

  l2d = compLoc2Dom(dim);

# if defined(DISSECTION_DEBUG)
  for (unsigned ii = 0; ii < dim; ii++)
      assert(_global2local[_local2global[ii]] == ii);
  // Save profile of the skeleton of the matrix with new ordering
  // into a file.
  fich = fopen("matrixProfile.dat", "w");
  fprintf(fich,"%d\n",dim);
  int nz = 0;
  fprintf(fich,"%d ",nz);
  for (int iLoc = 0; iLoc < dim; iLoc++) {      
      int iGlob = _local2global[iLoc];
      int nbCol = ptRows[iGlob+1]-ptRows[iGlob];
      nz += nbCol;
      fprintf(fich,"%d ",nz);
  }
  CHECK(nz == ptRows[dim], "Wrong number of non zero coefficients?");
  fprintf(fich,"\n");
  for (int iLoc = 0; iLoc < dim; iLoc++) {      
      int iGlob = _local2global[iLoc];
      for (int jCol = ptRows[iGlob]; jCol<ptRows[iGlob+1];jCol++) 
      {
	int jGlob = indCols[jCol];
	int jLoc  = _global2local[jGlob];
	fprintf(fich,"%d ",jLoc);
      }
  }
  fprintf(fich,"\n");
  fclose(fich);
  // l2d keeps the index of subdomain where node is included
  // ptDom keeps begining address of subdomain with sorting : nb_sbdm ... 1
  // l2d = compLoc2Dom(dim);
  unsigned* ptDom = new unsigned[(1<<_nbLevels)-1];
  fich = fopen("sparseMatrix3.dat", "w");
  fprintf(fich,"%d\n",dim);
  nz = 0;
  unsigned d = 0;
  fprintf(fich,"%d ",nz);
  for (long iLvl = _nbLevels-1; iLvl >=0; iLvl--) {
    unsigned begDom = (1<<iLvl);
    unsigned endDom = 2*begDom;
    for (unsigned iDom = begDom; iDom < endDom; iDom++) {
      ptDom[iDom-1] = d;
      d += _sizeOfDomains[iDom-1];
      for (unsigned iLoc = _ptOnDomains[iDom-1]; 
	   iLoc < _ptOnDomains[iDom-1]+_sizeOfDomains[iDom-1]; 
	   iLoc++)
      {	  
	  int iGlob = _local2global[iLoc];
	  int nbCol = ptRows[iGlob+1]-ptRows[iGlob];
	  nz += nbCol;
	  fprintf(fich,"%d ",nz);
      }
    }
  }
  assert(nz == ptRows[dim]);
  CHECK(nz == ptRows[dim], "Wrong number of non zero coefficients?");
  fprintf(fich,"\n");
  for (long iLvl = _nbLevels-1; iLvl >=0; iLvl--) {
    unsigned begDom = (1<<iLvl);
    unsigned endDom = 2*begDom;
    for (unsigned iDom = begDom; iDom < endDom; iDom++) {
      for (unsigned iLoc = _ptOnDomains[iDom-1]; 
	   iLoc < _ptOnDomains[iDom-1]+_sizeOfDomains[iDom-1]; 
	   iLoc++)
      {
	int iGlob = _local2global[iLoc];
	for (int jCol = ptRows[iGlob]; jCol<ptRows[iGlob+1];jCol++) {
	  int jGlob  = indCols[jCol];
	  int jLoc1  = _global2local[jGlob];
	  unsigned dl = l2d[jLoc1];
	  int jGlob1 = _local2global[jLoc1];
	  jLoc1 -= _ptOnDomains[dl-1];
	  assert(jLoc1 < _sizeOfDomains[dl-1]);
	  jLoc1 += ptDom[dl-1];
	  assert(jLoc1 < dim);
	  fprintf(fich,"%d ",jLoc1);
	}
      }
    }
  }
  fprintf(fich,"\n");
  fclose(fich);  
  delete [] ptDom;
# endif
  // 2. Compute local indices to domains containing this indice
  //    to speed-up symbolic computation
  unsigned* loc2dom = compLoc2Dom(dim);
  // 3. Allocation of lists of connection per levels of ancestors
  //    for each node of the bisection tree
  // NB : Some optimization to do here. Avoiding uniqueness method
  //      checking for each connected nodes if it's in the
  //      connection list, thanks to a logical array with size
  //      equal to indices ?
  std::list<unsigned>** connections = new std::list<unsigned>*[nbDoms];
  assert(connections != NULL);
  connections[0] = NULL;
  for (unsigned iDom = 2; iDom <= nbDoms; iDom++) {
    unsigned layer = nodeLayer(iDom);
    assert(layer < _nbLevels);
    assert(layer > 0);
    connections[iDom-1] = new std::list<unsigned>[layer];
    // a. Filling interconnection between iDom and his ancestors
    //
    //   i. For each node in iDom :
    for (unsigned idxVert = _ptOnDomains[iDom-1];
	 idxVert<_ptOnDomains[iDom-1]+_sizeOfDomains[iDom-1]; idxVert++)
    {
      unsigned indGlob = _local2global[idxVert];
      // ii. For global indice node indGlob, search connected
      //     nodes from the matrix skeleton :
      for (unsigned ptRow = ptRows[indGlob]; ptRow < ptRows[indGlob+1];
	   ptRow++)
      {
	unsigned connectGlobInd = indCols[ptRow];
	unsigned connectLocInd  = _global2local[connectGlobInd];
	unsigned ancestor = loc2dom[connectLocInd];
	connectLocInd -= _ptOnDomains[ancestor-1];

	if (ancestor < iDom) // If ancestor is an ancestor (and not itself or a child)
	{
	  unsigned lvlLocDom = nodeLayer(ancestor);
	  connections[iDom-1][lvlLocDom].push_back(connectLocInd);
	}
      }// End for (unsigned ptRow
    }// End for (unsigned idxVert
    for (unsigned iLayer = 0; iLayer < layer; iLayer++)
    {
      connections[iDom-1][iLayer].sort();
      connections[iDom-1][iLayer].unique();
    }
  }// End for (iDom
  // 4. Convert the nodes connections lists into strips 
  //    connections lists
  assert(connections != NULL);
  SetOfStrips** stripConnections =
      new SetOfStrips*[nbDoms];
  assert(stripConnections != NULL);
  stripConnections[0] = NULL;
  for (unsigned iDom = 2; iDom <= nbDoms; iDom++) {
    unsigned layer = nodeLayer(iDom);
    assert(layer < _nbLevels);
    stripConnections[iDom-1] = new SetOfStrips[layer];
    for (unsigned iLayer = 0; iLayer < layer; iLayer++)
    {
      assert(connections != NULL);
      stripConnections[iDom-1][iLayer]=
	  SetOfStrips(connections[iDom-1][iLayer]);
    }
  }
  // OK, initial connection is decribed as strip 
  // with local indices per domains.
  // 5. Begin symbolic factorization :
  _interco = symbolicFactorization(stripConnections);
  // 6. Renumbering father--children to have one parallel strip 
  //    and one sequential strip for interconnection between 
  //    children and father.
  //    Renumbering :
  //    Parallel indices from left child | Sequential indices from
  //    left and right children | Parallel indices from right children
  for (unsigned iDom = 1; iDom < (1<<(_nbLevels-1)); iDom++)
  {
    SetOfStrips paral, parar, seq;
    SetOfStrips::iterator itS;
    unsigned *perm = new unsigned[sizeOfDomain(iDom)];
    for (unsigned i = 0; i < sizeOfDomain(iDom); i++) 
      perm[i] = i;
    unsigned layer = nodeLayer(iDom);
    assert(layer < _nbLevels-1);
    unsigned leftChild  = (iDom<<1);
    unsigned rghtChild  = leftChild+1;
    SetOfStrips& lStrip = _interco[leftChild-1][layer];
    SetOfStrips& rStrip = _interco[rghtChild-1][layer];
    lStrip.para_union(rStrip, paral, parar,seq);
    renumberingInterface(sizeOfDomain(iDom), perm, paral, parar,
			 seq, layer, iDom );
    unsigned* invperm = new unsigned[sizeOfDomain(iDom)];
    for (unsigned ii = 0; ii < sizeOfDomain(iDom); ii++)
      invperm[perm[ii]] =  ii;

    unsigned curDom = iDom;
    unsigned nbDoms = 1;
    for (unsigned iLayer = layer+1; iLayer < _nbLevels; iLayer++)
    {
      curDom *= 2;
      nbDoms *= 2;
      for (unsigned jDom = curDom; jDom < curDom+nbDoms; jDom++)
      {
	// Retrieve connection between jDom and iDom
	// in old renumbering
	SetOfStrips& strips = _interco[jDom-1][layer];
	// and convert into list 
	std::list<unsigned> locConnect = strips;
	// Convert old renumbering into new renumbering in the list
	for (std::list<unsigned>::iterator itL = locConnect.begin();
	     itL != locConnect.end(); itL++) {
	  (*itL) = invperm[(*itL)];
	}
	locConnect.sort();
	// Build new strips coming from new renumbering
	strips = SetOfStrips(locConnect);
      }// for (jDom
    }// for (iLayer
    // Renumbering the local indices for subdomain iDom :
    // ------------------------------------------------
    int* tempInd = new int[sizeOfDomain(iDom)];
    memcpy(tempInd, _local2global+ _ptOnDomains[iDom-1],
	   sizeof(int)*sizeOfDomain(iDom));
    for (unsigned i = 0; i < sizeOfDomain(iDom); i++) {
	_local2global[_ptOnDomains[iDom-1]+i] =
	  tempInd[perm[i]];
    }
    delete [] tempInd; // 01 Feb.2013 Atsushi
    delete [] invperm;
    delete [] perm;

    for (unsigned i = 0; i < sizeOfDomain(iDom); i++) {
      int l2g = _local2global[_ptOnDomains[iDom-1]+i];
      _global2local[l2g] = _ptOnDomains[iDom-1]+i;
    }
  }// for (iDom
  // At last, create interconnections between a node set of strips
  // and his ancestors (out of diagonal and diagonals blocks).
  // Nb : For layer 0, no ancestors, so...
  for (unsigned iLayer = 1; iLayer < _nbLevels; iLayer++) {
    // For each domain contained in layer iLayer :
    for (unsigned iDom = (1<<iLayer); iDom < (1<<(iLayer+1)); iDom++) {
      FathersStrips& connection_i = _interco[iDom-1];
      unsigned jDom = iDom;
      // For each ancestor of our domain iDom (beginning from the father to root) :
      for (unsigned jLayer = iLayer; jLayer > 0; jLayer--) {
	// Compute index of the ancestor in layer jLayer :
	jDom = jDom/2;
	FathersStrips& connection_j = _interco[jDom-1];
	// Retrieve the set of strips describing the direct connection between the
	// strips 
	for (unsigned kLayer = jLayer; kLayer > 0; kLayer --) {
	  SetOfStrips& stripForRowBlock_jk = connection_i.getRowSetOfStrips(jLayer-1,kLayer-1);
	  SetOfStrips& stripForColBlock_jk = connection_i.getColSetOfStrips(jLayer-1,kLayer-1);
	  // Most of data are copy of out of diagonal strips :
	  stripForRowBlock_jk = connection_i[jLayer-1];
	  stripForColBlock_jk = connection_i[kLayer-1];
	  if (kLayer != jLayer) {// In this case, must change the destination index :
	    SetOfStrips& outOfDiag_strips_jk = connection_j[kLayer-1];
	    SetOfStrips::iterator it_DstBjk = stripForColBlock_jk.begin();
	    for (SetOfStrips::iterator it_Bjk = outOfDiag_strips_jk.begin();
		 it_Bjk != outOfDiag_strips_jk.end(); it_Bjk++) {
	      unsigned decal;
	      while ((it_DstBjk != stripForColBlock_jk.end()) &&
		     (*it_DstBjk).inside(*it_Bjk,decal)) {
		(*it_DstBjk).begin_dst = (*it_Bjk).begin_src + decal;
		it_DstBjk ++;
	      }
	    }
	  }
	}
      }
    }
  }
# if defined(DISSECTION_DEBUG)
  {
    unsigned iDom;
  // DEBUGGING : Save filling of the matrix
    fich = fopen("filledMatrix2.dat", "w");
    fprintf(fich,"%d\n",_nbLevels);
    for (iDom = 1; iDom <= nbDoms; iDom++)
      fprintf(fich,"%d %d\n",_ptOnDomains[iDom - 1],
	      _sizeOfDomains[iDom - 1]);
    // ancestors connection :
    // added by Atsush :
    // iDom ==1 does not have father, hence it is not shown here
    // strips are shown in order of Idom : {2, 3} {4, 5, 6, 7}
    // NB : index of _interco[] starts at 0, (iDom - 1) is obligation!!
    for (unsigned iLevel = 1; iLevel < _nbLevels; iLevel++)
      {
	// First domain indice for this level:
	unsigned begDom = (1<<iLevel);
	// Last domain indice for this level + 1:
	unsigned endDom = 2*begDom; // (1 << (iLevel + 1)); 
	for (iDom = begDom; iDom < endDom; iDom++) {
	  fprintf(fich,"# iDom = %d\n", iDom);
	  for (unsigned iLayer = 0; iLayer < iLevel; iLayer++) {
	    fprintf(fich,"# iLayer = %d\n", iLayer);
	    fprintf(fich,"%d\t",_interco[iDom-1][iLayer].numberOfStrips());
	    for (SetOfStrips::iterator itS = _interco[iDom-1][iLayer].begin();
		 itS != _interco[iDom-1][iLayer].end(); itS++) {
	      fprintf(fich, "(%d %d) %d\t",
		      (*itS).begin_src, (*itS).begin_dst, (*itS).width);
	    }
	    fprintf(fich,"\n");
	  }// for iLayer
	}// for iDom
      }// for iLevel
    fclose(fich);
  }
#endif
  // Atsushi : 03 Mar.2012
#if 0
  clock_gettime(CLOCK_REALTIME, &tt3);
#endif
  //
  std::map<int, int>* indCols_bylvl = new std::map<int, int>[_nbLevels];

  int *ptRow_diag = new int[dim];
  int nnzh = (ptRows[dim] + dim) / 2;
  int *indVals_tmp = new int [nnzh];
    
  indCols_sbdmn = new int[ptRows[dim]];
  indCols_idxStrip = new int[ptRows[dim]];
  indCols_idxSbdmn = new int[ptRows[dim]];

#if 0
  fich = fopen("debug.dat", "w");
#endif
  for (int n = (_nbLevels - 1); n >=0; n--) {
    unsigned begDom = (1 << n);
    unsigned endDom = 2 * begDom; // (1 << (n + 1)); 
    unsigned dd;
    for (unsigned d = begDom; d < endDom; d++) {
      const unsigned d1 = d - 1;  // to access C array
      for (unsigned i = _ptOnDomains[d1]; 
	   i < (_ptOnDomains[d1] + _sizeOfDomains[d1]); 
	   i++) {
	for (int m = 0; m <= n; m++) {
	  indCols_bylvl[m].clear();
	}
	const int ii = _local2global[i];
	for (int k = ptRows[ii]; k < ptRows[ii + 1]; k++) {
	  const int j = indCols[k];
	  const int jj = _global2local[j];
	  dd = l2d[jj];
	  if (dd <= d) {
	    unsigned ndl = nodeLayer(dd);
	    indCols_bylvl[ndl].insert(std::map<int,int>::value_type(jj, k));
	  }
	} // loop : k
#if 0
	dd = d;
	for (int m = n; m >= 0; m--) {
	  const unsigned dd1 = dd - 1;
	  if (indCols_bylvl[m].size() > 0) {
	    fprintf(fich, 
		    "i = %d (global = %d) sbdmn = %d level = %d size = %d\n", 
		    i, ii, l2d[i], m, indCols_bylvl[m].size());
	    for (std::map<int, int>::iterator it = indCols_bylvl[m].begin();
		 it != indCols_bylvl[m].end();
		 ++it) {
	      fprintf(fich, "( %d %d ) ", (*it).first, (*it).second);
	    }
	    fprintf(fich, "\n");
	    if (m == n) {
	      fprintf(fich, "dd = %d m = %d : %d %d\n", 
		      dd, m, _ptOnDomains[dd1],
		      _sizeOfDomains[dd1]);
	    }
	    else {
	      fprintf(fich, "dd = %d m = %d : %d ", 
		      dd, m, _ptOnDomains[dd1]);
	      for (SetOfStrips::iterator itQ = _interco[d1][m].begin();
		   itQ != _interco[d1][m].end();
		   ++itQ) {
		fprintf(fich, "(%d %d %d) ", 
			(*itQ).begin_src,
			(*itQ).begin_dst,
			(*itQ).width);
	      }
	      fprintf(fich, "\n");
	    }
	  }
	  dd /= 2;
	}
#endif
	// m == n : diagonal
	dd = d;
	for (std::map<int, int>::iterator it = indCols_bylvl[n].begin();
	     it != indCols_bylvl[n].end();
	     ++it) {
	  indCols_sbdmn[(*it).second] = dd;
	  indCols_idxStrip[(*it).second] = // within one strip
	  indCols_idxSbdmn[(*it).second] = 
	    ((*it).first - _ptOnDomains[dd - 1]);
	}
      	for (int m = (n - 1); m >= 0; m--) {
	  dd /= 2;
	  const unsigned dd1 = dd - 1;
	  SetOfStrips::iterator is = _interco[d1][m].begin();
	  int width_strips = 0;
	  for (std::map<int, int>::iterator it = indCols_bylvl[m].begin();
	       it != indCols_bylvl[m].end();
	       ++it) {
	    int endStrip = ((*is).begin_dst + _ptOnDomains[dd1] + 
			    (*is).width);
	    while (endStrip < (*it).first) {
	      width_strips += (*is).width;
	      ++is;
	      if (is == _interco[d1][m].end()) {
#if 0
		  berr = false;
		  delete [] _sizeOfDomains;
		  delete [] _ptOnDomains;
		  delete [] _local2global;
		  delete [] _global2local;
		  if (_interco) delete [] _interco;
		  return; 
		// 4 Jul.2013 Atsushi
		  throw std::runtime_error("Failed map CRS in strip");
#endif
	      }
	      else {
		endStrip = ((*is).begin_dst + _ptOnDomains[dd1] + 
			    (*is).width);
	      }
	    }
	    const int begStrip = (*is).begin_dst + _ptOnDomains[dd1];
	    if ((begStrip <= (*it).first) && ((*it).first < endStrip)) {
#if 0
	      fprintf(fich, "(%d %d : %d [%d %d])\n", 
		      (*it).first, (*it).second, dd, begStrip, (*is).width);
#endif
	      indCols_sbdmn[(*it).second] = dd;
	      indCols_idxStrip[(*it).second] = 
		(((*it).first - begStrip) + width_strips);
	      indCols_idxSbdmn[(*it).second] = 
		((*it).first - _ptOnDomains[dd1]);
	    }
	  } // loop : it
	} // loop : m
      }  // loop : i
    }
  }
#if 0
  fclose(fich);
#endif

#if 0
  fich = fopen("sparseMatrix2.dat", "w");
  fprintf(fich, "dim = %d\n", dim);

  for (int n = (_nbLevels - 1); n >=0; n--) {
    unsigned begDom = (1 << n);
    unsigned endDom = 2 * begDom; // (1 << (n + 1)); 
    for (unsigned d = begDom; d < endDom; d++) {
      for (int ii = _ptOnDomains[d - 1]; 
	   ii < (_ptOnDomains[d - 1] + _sizeOfDomains[d - 1]); 
	   ii++) {
	const int i = _local2global[ii];
	fprintf(fich, "i = %d (%d) ptRows[] = %d width = %d\n", 
		i, 
		_global2local[i],
		ptRows[i], (ptRows[i + 1] - ptRows[i]));
	for (int k = ptRows[i]; k < ptRows[i + 1]; k++) {
	  int jj = _global2local[indCols[k]];
	  int dd = l2d[jj];
	  if (dd <= d) {
	  //	  fprintf(fich, "%d %d [%d : %d %d]\n", 
	    if ((indCols_sbdmn[k] >= 1) &&
		(indCols_sbdmn[k] < (1 << (_nbLevels + 1)))) {
	      const int itmp = (indCols_idxSbdmn[k] +
				_ptOnDomains[indCols_sbdmn[k] - 1]);
	      fprintf(fich, "[%d : %d] %d(%d) %d [%d->%d :: %d : %d %d]", 
		      d, dd,
		      indCols[k],
		      _global2local[indCols[k]],
		      itmp,
		      k,
		      toSym[k],
		      indCols_sbdmn[k],
		      indCols_idxStrip[k],
		      indCols_idxSbdmn[k]);
	      if (_global2local[indCols[k]] != itmp) {
		fprintf(fich, " mismatch\n");
	      }
	      else {
		fprintf(fich, "\n");
	      }
	    }
	    else {
	      fprintf(fich, 
		      "err [%d : %d] i=%d(%d) k=%d j=%d(%d) [%d : %d %d] ", 
		      d, dd,
		      i,           // unsymmetirized CSR
		      ii,          // permuted by METIS 
		      k,
		      indCols[k],  // unsymmetirized CSR
		      jj,          // permuted by METIS
		      indCols_sbdmn[k],
		      indCols_idxStrip[k],
		      indCols_idxSbdmn[k]);
	      unsigned m = nodeLayer(dd);
	      fprintf(fich, "%d (%d th father = %d)\n",
		      _ptOnDomains[nthfatherIndex(d, (n - m)) - 1],
		      (n - m),
		      nthfatherIndex(d, (n - m)));
	      for (SetOfStrips::iterator itQ = _interco[d - 1][m].begin();
		  itQ != _interco[d - 1][m].end();
		   ++itQ) {
		fprintf(fich, "(%d %d %d) ", 
			(*itQ).begin_src,
			(*itQ).begin_dst,
			(*itQ).width);
	      }
	      fprintf(fich, "\n");
	    }
	  } // if (dd <= d)
	} // loop : j
      }   // loop : i
    }     // loop : d
  }       // loop : n
  fclose(fich);
#endif
  //
#if 0
  for (int i = 0; i < dim; i++) {
    bool flag = false;
    for (int k = ptRows[i]; k < ptRows[i + 1]; k++) {
      if (i == indCols[k]) {
	ptRow_diag[i] = k;
	flag = true;
	break;
      }
    }
    assert(flag);
    CHECK(flag, "no diagonal entry in unsymmeterized CSR\n");
  }
#endif
#if 0
  clock_gettime(CLOCK_REALTIME, &tt4);
#endif
  // gnereating map from diagonal matrix / offdiagonal strips to matrix value 
  // in original CSR format
  int itmp = 0;
  for (int d = 0; d < nbDoms; d++) {
    if (itmp < _sizeOfDomains[d]) {
      itmp =  _sizeOfDomains[d];
    }
  }
  if (verbose) {
    fprintf(fp, "max of size of subdomains = %d\n", itmp);
  }
  std::map<int, int> *sparse_row_diag = new std::map<int, int>[itmp];
  std::map<int, int> *sparse_row_offdiag = new std::map<int, int>[itmp];
  
  _csr_diag = new CSR_indirect[nbDoms];
  _csr_offdiag = new CSR_indirect[nbDoms];
  _sizeOfFathersStrips = new int[nbDoms];

  int *width_strips_sbdmn = new int[_nbLevels];

  for (int n = (_nbLevels - 1); n >=0; n--) {
    unsigned begDom = (1 << n);
    unsigned endDom = 2 * begDom; // (1 << (n + 1)); 
    for (unsigned d = begDom; d < endDom; d++) {
      const int d1 = d - 1;
      for (int m = nodeLayer(d) - 1; m >= 0; m--) {
	width_strips_sbdmn[m] = 0;
      }
      for (int m = (nodeLayer(d) - 1); m > 0; m--) {
	width_strips_sbdmn[m - 1] = width_strips_sbdmn[m]; 
	for (SetOfStrips::iterator it = _interco[d1][m].begin();
	       it != _interco[d1][m].end(); ++it) {
	  width_strips_sbdmn[m - 1] += (*it).width;
	}
      }
      for (int i = 0; i < _sizeOfDomains[d1]; i++) {
	sparse_row_diag[i].clear();
	sparse_row_offdiag[i].clear();
      }
      int ii0, ii;
      for (ii0 = 0, ii = _ptOnDomains[d1];
	   ii0 < _sizeOfDomains[d1]; 
	   ii++, ii0++) {
	const int i = _local2global[ii];
	
	for (int k = ptRows[i]; k < ptRows[i + 1]; k++) {
	  int jj = _global2local[indCols[k]];
	   unsigned dd = l2d[jj];
	  if (dd == d) {
#ifdef DEBUG_MAPPING_CSR
	    sparse_row_diag[ii0].insert(std::map<int, int>::value_type(indCols_idxStrip[k], k));
#else
	    sparse_row_diag[ii0].insert(std::map<int, int>::value_type(indCols_idxStrip[k], k)); //(isLower ? toSym[k] : k)));
#endif
	  }
	  else {
	    if (dd < d) {
#ifdef DEBUG_MAPPING_CSR
	      sparse_row_offdiag[ii0].insert(std::map<int, int>::value_type(indCols_idxStrip[k] + width_strips_sbdmn[nodeLayer(dd)], k));
#else
	      sparse_row_offdiag[ii0].insert(std::map<int, int>::value_type(indCols_idxStrip[k] + width_strips_sbdmn[nodeLayer(dd)], k)); //(isLower ? toSym[k] : k)));
#endif
	    }
	  }
	} // loop : k
      }   // loop : ii, ii0
      _csr_diag[d1].n = _sizeOfDomains[d1];
      _csr_offdiag[d1].n = _sizeOfDomains[d1];
      _csr_diag[d1].ptRows = new int [_csr_diag[d1].n + 1];
      _csr_offdiag[d1].ptRows = new int [_csr_offdiag[d1].n + 1];
      _csr_diag[d1].ptRows[0] = 0;
      _csr_offdiag[d1].ptRows[0] = 0;
      for (int i = 0; i < _csr_diag[d1].n; i++) {
	_csr_diag[d1].ptRows[i + 1] = (_csr_diag[d1].ptRows[i] + 
				      sparse_row_diag[i].size());
	_csr_offdiag[d1].ptRows[i + 1] = (_csr_offdiag[d1].ptRows[i] + 
					 sparse_row_offdiag[i].size());
      }
      _csr_diag[d1].nnz = _csr_diag[d1].ptRows[_csr_diag[d1].n];
      _csr_offdiag[d1].nnz = _csr_offdiag[d1].ptRows[_csr_offdiag[d1].n];
      _csr_diag[d1].indCols = new int [_csr_diag[d1].nnz];
      _csr_offdiag[d1].indCols = new int [_csr_offdiag[d1].nnz];
      _csr_diag[d1].indVals = new int [_csr_diag[d1].nnz];
      _csr_offdiag[d1].indVals = new int [_csr_offdiag[d1].nnz];
      _csr_diag[d1].indVals0 = new int [_csr_diag[d1].nnz];
      _csr_offdiag[d1].indVals0 = new int [_csr_offdiag[d1].nnz];
      if(!isSym) {
	_csr_offdiag[d1].indVals_unsym = new int [_csr_offdiag[d1].nnz];
      }
#ifdef DEBUG_MAPPING_CSR
      _csr_diag[d1].indVals2 = new int [_csr_diag[d1].nnz];
      _csr_offdiag[d1].indVals2 = new int [_csr_offdiag[d1].nnz];
#endif      
      int k0, k1;
      k0 = 0; 
      k1 = 0;

      for (int i = 0; i < _csr_diag[d1].n; i++){
	for (std::map<int, int>::iterator it = sparse_row_diag[i].begin();
	     it != sparse_row_diag[i].end(); ++it, k0++) {
	  _csr_diag[d1].indCols[k0] = (*it).first;
#ifdef DEBUG_MAPPING_CSR
	  _csr_diag[d1].indVals2[k0] = (*it).second;
	  _csr_diag[d1].indVals[k0] = toSym[(*it).second];
#else
	  _csr_diag[d1].indVals0[k0] = (*it).second; //map_indcols[(*it).second];
	  _csr_diag[d1].indVals[k0] = indVals[(*it).second]; //remap_indcols[(*it).second];
#endif
	}
	for (std::map<int, int>::iterator it = sparse_row_offdiag[i].begin();
	     it != sparse_row_offdiag[i].end(); ++it, k1++) {
	  _csr_offdiag[d1].indCols[k1] = (*it).first;
#ifdef DEBUG_MAPPING_CSR
	  _csr_offdiag[d1].indVals2[k1] = (*it).second;
	  _csr_offdiag[d1].indVals[k1] = toSym[(*it).second];
#else
	  _csr_offdiag[d1].indVals0[k1] = (*it).second; //remap_indcols[(*it).second];
	  _csr_offdiag[d1].indVals[k1] = indVals[(*it).second]; //remap_indcols[(*it).second];
#endif
	}
      }
      if (!isSym) {
	int k = 0;
	for (int i = 0; i < _csr_diag[d1].n; i++){
	  for (std::map<int, int>::iterator it = sparse_row_offdiag[i].begin();
	       it != sparse_row_offdiag[i].end(); ++it, k++) {
	    _csr_offdiag[d1].indVals_unsym[k] = unsym_upper2lower[(*it).second];
	    //	      unsym_upper2lower[(*it).second];
	  }
	}
      } // if (!isSym)
    } // loop : d
  } // loop : n 


  delete [] sparse_row_diag;
  delete [] sparse_row_offdiag;
  //  delete [] l2d;
#if 0
  fich = fopen("sparseMatrix4.dat", "w");
  fprintf(fich, "dim = %d\n", dim);
  fprintf(fich, "# of subdomais = %d\n", nbDoms);
  for (int d = 0; d < nbDoms; d++) {
    fprintf(fich, "subdomain = %d : %d diagonal n = %d nnz = %d\n", 
	    d, _sizeOfDomains[d],
	    _csr_diag[d].n, _csr_diag[d].nnz);
    for (int i = 0; i < (_csr_diag[d].n + 1); i++) {
      fprintf(fich, "%d ", _csr_diag[d].ptRows[i]);
    }
    fprintf(fich, "\n");
    for (int i = 0; i < _csr_diag[d].n; i++) {
      for (int k = _csr_diag[d].ptRows[i]; 
	   k < _csr_diag[d].ptRows[i + 1]; 
	   k++) {
	fprintf(fich, "(%d %d) ",
		_csr_diag[d].indCols[k], _csr_diag[d].indVals0[k]);
      }
    }
    fprintf(fich, "\n");
    fprintf(fich, "subdomain = %d : %d offdiagonal n = %d nnz = %d\n", 
	    d, _sizeOfDomains[d],
	    _csr_offdiag[d].n, _csr_offdiag[d].nnz);
    for (int i = 0; i < (_csr_offdiag[d].n + 1); i++) {
      fprintf(fich, "%d ", _csr_offdiag[d].ptRows[i]);
    }
    fprintf(fich, "\n");
    for (int i = 0; i < _csr_offdiag[d].n; i++) {
      for (int k = _csr_offdiag[d].ptRows[i]; 
	   k < _csr_offdiag[d].ptRows[i + 1]; 
	   k++) {
	fprintf(fich, "(%d %d) ",
		_csr_offdiag[d].indCols[k], _csr_offdiag[d].indVals0[k]);
      }
    }
    fprintf(fich, "\n");
  } // loop : d
#endif
#if 0
  if (isSym) {
    for (int i = 0; i < dim; i++) {
      bool flag = false;
      for (int k = ptRows[i]; k < ptRows[i + 1]; k++) {
	if (i == indCols[k]) {
	  ptRow_diag[i] = k;
	flag = true;
	break;
	}
      }
      assert(flag);
      CHECK(flag, "no diagonal entry in unsymmeterized CSR\n");
    }
    // prepare CSR for upper symmetric format
    for (int i = 0; i < nnzh; i++) {
      indVals_tmp[i] = 0;
    }
    //    fprintf(stderr, "half of nnz = %d\n", nnzh);
    for (int d = 0; d < nbDoms; d++) {
      for (int i = 0; i < _csr_diag[d].n; i++) {
	for (int k = _csr_diag[d].ptRows[i]; k < _csr_diag[d].ptRows[i + 1]; 
	     k++) {
	  if (_csr_diag[d].indVals0[k] >= nnzh) {
	    fprintf(stderr, 
		    "too large element in upper symmetric %d : %d %d -> %d\n", 
		    d, i, k, _csr_diag[d].indVals0[k]);
	  }
	  else {
	    indVals_tmp[_csr_diag[d].indVals0[k]]++;
	  }
	}
	for (int k = _csr_offdiag[d].ptRows[i]; 
	     k < _csr_offdiag[d].ptRows[i + 1];
	     k++) {
	  if (_csr_offdiag[d].indVals0[k] >= nnzh) {
	    fprintf(stderr, 
		    "too large element in upper symmetric %d : %d %d -> %d\n", 
		    d, i, k, _csr_offdiag[d].indVals0[k]);
	  }
	  else {
	    indVals_tmp[_csr_offdiag[d].indVals0[k]]++;
	  }
	}
      }
    }
    for (int k = 0; k < nnzh; k++) {
      if (indVals_tmp[k] == 0) {
	fprintf(stderr, "non touched element in upper symmetric %d\n", k);
      }
    }
  }
#endif
#if 0
  clock_gettime(CLOCK_REALTIME, &tt5);
#endif
  _nbDoms = nbDoms;
  for (int d = 0; d < nbDoms; d++) {
    _sizeOfFathersStrips[d] = 0;
    int level = nodeLayer(d + 1);
    for (int m = level - 1; m >= 0; m--) {
      SetOfStrips &strps = _interco[d][level - 1 - m];
      _sizeOfFathersStrips[d] += strps.numberOfIndices();
    }
  }
    
  // Create glob 2 dom links :
  _glob2dom = new int[dim];
  for (unsigned iDom = 0; iDom < (1<<_nbLevels)-1; iDom++) {
    for (unsigned iNode = _ptOnDomains[iDom];
	 iNode < _ptOnDomains[iDom+1]; iNode++) {
      _glob2dom[_local2global[iNode]] = iDom + 1;
    }
  }

  _loc2glob_diag = new int*[nbDoms];
  _loc2glob_offdiag = new int *[nbDoms];
  for (int d = 0; d < nbDoms; d++) { 
    _loc2glob_diag[d] = new int[_sizeOfDomains[d]];
    int ii = _ptOnDomains[d];
    for (int i = 0; i < _sizeOfDomains[d]; i++, ii++) {
//    _loc2glob_diag[d][i] = _local2global[ii];  
//    16 Jul.2013 Atsushi getDiagLoc2Glob(nd)[] = loc2glob_diag[nd][]
      _loc2glob_diag[d][i] = remap_eqn[_local2global[ii]]; // _local2global[ii]; //
                                  // 19 Dec. 
    }
    int itmp = 0;
    FathersStrips &fstrps = _interco[d];
    const int layer = nodeLayer(d + 1);
    for (int m = (layer - 1); m >= 0; m--) {
      itmp += fstrps[m].numberOfIndices();
    }
    _loc2glob_offdiag[d] = new int[_sizeOfFathersStrips[d]];
    int i0 = 0;  // index inside of strips which are stored contiguously
    for (int m = (layer - 1); m >= 0; m--) {
      for (SetOfStrips::iterator it = fstrps[m].begin();
	   it != fstrps[m].end(); ++it) {
	const int n = nthfatherIndex((d + 1), (layer - m));
	ii = (_ptOnDomains[n - 1] + (*it).begin_dst);
	for (int j = 0; j < (*it).width; j++, i0++, ii++) {
//	  _loc2glob_offdiag[d][i0] = _local2global[ii];
//    16 Jul.2013 Atsushi getOffdiagLoc2Glob(nd)[] = loc2glob_offdiag[nd][]
	  _loc2glob_offdiag[d][i0] = remap_eqn[_local2global[ii]]; // _local2global[ii]; 
	                          // 19 Dec. 
	}
      }
    }
  }

#if 0
  draw_csr("matrix-debug.ps", dim, 
	   ptRows, indCols, ptRow_diag, indVals, l2d);
  exit(-1);
#endif

  // Destroy temporary data
  for (int d = 0; d < nbDoms; d++) {
    delete [] _csr_diag[d].indVals0;
    delete [] _csr_offdiag[d].indVals0;
    _csr_diag[d].indVals0 = NULL;
    _csr_offdiag[d].indVals0 = NULL;
  }
  for (unsigned iDom = 2; iDom <= nbDoms; iDom++) {
    delete [] stripConnections[iDom-1];
    delete [] connections[iDom-1];
  }
  delete [] stripConnections;
  delete [] connections;
  delete [] loc2dom;
  delete [] indCols2;
  delete [] ptRows2 ;
  delete [] l2d;
  // 01 Feb.2013 : Atsushi -- begin
  delete [] indVals_tmp;
  delete [] ptRow_diag;
  delete [] indCols_bylvl;
  delete [] indCols_sbdmn;
  delete [] indCols_idxStrip;
  delete [] indCols_idxSbdmn;
  delete [] width_strips_sbdmn;
  // 01 Feb.2013 : Atsushi -- end
  berr = true;
#if 0
  fprintf(stderr, "elapsed time (sec.): \n");
  fprintf(stderr, "preprocess:         %g\n", 
	  ((tt1.tv_sec - tt0.tv_sec) + (tt1.tv_nsec - tt0.tv_nsec) * 1.0e-9));
  fprintf(stderr, "%6s:             %g\n", 
	  spltFct == &ScotchSplitter ? "SCOTCH" : "METIS ",
	  ((tt2.tv_sec - tt1.tv_sec) + (tt2.tv_nsec - tt1.tv_nsec) * 1.0e-9));
  fprintf(stderr, "generating strips:  %g\n",
	  ((tt3.tv_sec - tt2.tv_sec) + (tt3.tv_nsec - tt2.tv_nsec) * 1.0e-9));
  fprintf(stderr, "+ CSR -> strip:    %g\n",
  	  ((tt4.tv_sec - tt3.tv_sec) + (tt4.tv_nsec - tt3.tv_nsec) * 1.0e-9));
  fprintf(stderr, "+ local CSR index: %g\n",
	  ((tt5.tv_sec - tt4.tv_sec) + (tt5.tv_nsec - tt4.tv_nsec) * 1.0e-9));
  fprintf(stderr, "postprocess index:  %g\n",
	  (((tt4.tv_sec - tt3.tv_sec) + (tt4.tv_nsec - tt3.tv_nsec) * 1.0e-9)+
	   ((tt5.tv_sec - tt4.tv_sec) + (tt5.tv_nsec - tt4.tv_nsec) * 1.0e-9)));
#endif
}
// --------------------------------------------------------------
Tree::~Tree()
{
  delete [] _sizeOfDomains;
  delete [] _ptOnDomains;
  delete [] _local2global;
  delete [] _global2local;
  if (_interco) delete [] _interco;
  //unsigned nDoms = (1<<_nbLevels)-1;
 // 01 Feb.2013 : Atsushi -- begin
  delete [] _glob2dom;
  delete [] _sizeOfFathersStrips;
  for (int i = 0; i < _nbDoms; i++) {
    delete [] _loc2glob_diag[i];
    delete [] _loc2glob_offdiag[i];
    delete [] _csr_diag[i].ptRows;
    delete [] _csr_diag[i].indCols;
    delete [] _csr_diag[i].indVals;
    //    delete [] _csr_diag[i].indVals0;
    delete [] _csr_offdiag[i].ptRows;
    delete [] _csr_offdiag[i].indCols;
    delete [] _csr_offdiag[i].indVals;
    //    delete [] _csr_offdiag[i].indVals0;
  }
  if (!_isSym) {
    for (int i = 0; i < _nbDoms; i++) {
      delete [] _csr_offdiag[i].indVals_unsym;
    }
  }
  delete [] _loc2glob_diag;
  delete [] _loc2glob_offdiag;
  delete [] _csr_diag;
  delete [] _csr_offdiag;
 // 01 Feb.2013 : Atsushi -- end
}
// ==============================================================
bool
Tree::save(FILE* stream) const
{
    // To do
    return true;
}
// --------------------------------------------------------------
bool
Tree::load(FILE* stream)
{
    // To do
    return true;
}
// --------------------------------------------------------------
void
Tree::printInfo(FILE* stream, unsigned verboseLevel) const
{
  unsigned l,n,p;
  fprintf(stream,"Number of bisection level : %u\n",_nbLevels);
  for (l=0; l < _nbLevels; l++) {
    double meanSize = 0., varSize = 0.;
    for (n = 1<<l; n < 1<<(l+1); n++) {
	unsigned szDom = sizeOfDomain(n);
	meanSize += szDom;
	varSize  += szDom*szDom;
    }
    meanSize /= (1<<l);
    varSize = varSize/(1<<l) - meanSize*meanSize;
    fprintf(stream, "Mean number of nodes at level %u : %lg\n",l,
	    meanSize);
    fprintf(stream,
	    "Variance for number of nodes at level %u : %lg\n",
	    l, varSize);
  }
  if (verboseLevel <= 1) return;
  for (n=1; n < (1<<_nbLevels); n++) {
    unsigned idxStart = _ptOnDomains[n-1];
    unsigned idxEnd   = idxStart + sizeOfDomain(n);
    fprintf(stream, 
	    "Global index of %d nodes in subdomain %u :\n",
	    _sizeOfDomains[n-1],n);
    for (p=idxStart; p < idxEnd; p++) {
      fprintf(stream,"%d\t",_local2global[p]);	  
    }
    fprintf(stream, "\n");
  }
}
// ==============================================================
void
Tree::removeLoops(int dim,
		  const int* ptRows, const int* indCols,
		  int*& ptRows2, int*& indCols2)
{
    // Optimize with bsearch or equivalent ?
    int nz = ptRows[dim];
    ptRows2 = new int[dim+1];
    indCols2= new int[nz-dim];
    ptRows2[dim] = nz-dim;
    for (int iRow = 0; iRow < dim; iRow++) {
	ptRows2[iRow] = ptRows[iRow]-iRow;
	int jColInd2 = ptRows2[iRow];
	for (int jColInd = ptRows[iRow];
	     jColInd < ptRows[iRow+1];
	     jColInd++)
	{
	    if (indCols[jColInd] != iRow) {
		indCols2[jColInd2] = indCols[jColInd];
		jColInd2++;
		CHECK(jColInd2<=nz-dim,
		      "jColInd2 is greater than number of non zero elements !"); 
	    }
	}
    }
}
// ==============================================================
FathersStrips*
Tree::symbolicFactorization(SetOfStrips** connections)
{
  bool verbose = false;
  unsigned iLevel, iDom, layer, iLayer, jLayer;
  unsigned nbDoms = (1U<<_nbLevels)-1;
  if (verbose) {
    fprintf(stderr, "nbDoms = %d\n", nbDoms);
  }
  // Initialize initial filling for each block :
  FathersStrips* factConnect = new FathersStrips[nbDoms];
  assert(factConnect != NULL);
  for (iDom = 2; iDom <= nbDoms; iDom++)
  {// For iDom==1, nothing to do...
    layer = nodeLayer(iDom);
    if (verbose) {
      fprintf(stderr, "iDom  %d at level %d\n", iDom, layer);
    }
    factConnect[iDom-1].alloc(layer);
    for (iLayer = 0; iLayer < layer; iLayer++)
    {
      if (verbose) {
	fprintf(stderr, "initial connection for %d with layer %d\n",
		iDom, iLayer);
      }
      // Initial connection block with his ancestors :
      factConnect[iDom-1][iLayer] = connections[iDom-1][iLayer];
    }
  }
  // Start symbolic factorization :
  for (iLevel = _nbLevels-1; iLevel > 0; iLevel--)
  {
    // First domain indice for this level:
    unsigned begDom = (1<<iLevel);
    // Last domain indice for this level + 1:
    unsigned endDom = 2*begDom;
    // Gauss symbolic elimination for each domain iDom at level iLevel
    for (iDom = begDom; iDom < endDom; iDom++)
    {
      if (verbose) {
	fprintf(stderr, "Eliminating dom %d \n", iDom);
      }
      unsigned ancestDom = iDom;
      for (iLayer = iLevel-1; iLayer > 0; iLayer--) {
	ancestDom = ancestDom/2;// Compute indice (starting to 1)
	                        // of next ancestor :
	if (verbose) {
	  fprintf(stderr, "ancestDom = %d\niLayer = %d\n", ancestDom, iLayer);
	}
	for (jLayer = iLayer; jLayer > 0; jLayer--) {
	  if (verbose) {
	    fprintf(stderr, "jLayer = %d\n", (jLayer - 1));
	  }
  	  factConnect[ancestDom-1][jLayer-1] += 
	                        factConnect[iDom-1][jLayer-1];
	}
      }
    }
  }
  // end symbolic factorization
  // ............................................................
# if defined(DISSECTION_DEBUG)
  // DEBUGGING : Save filling of the matrix
  FILE* fich = fopen("filledMatrix.dat", "w");
  fprintf(fich,"%d\n",_nbLevels);
  for (iDom = 0; iDom < nbDoms; iDom++)
    fprintf(fich,"%d %d\n",_ptOnDomains[iDom],
	    _sizeOfDomains[iDom]);
  // ancestors connection :
  for (iLevel = 1; iLevel < _nbLevels; iLevel++)
  {
    // First domain indice for this level:
    unsigned begDom = (1<<iLevel);
    // Last domain indice for this level + 1:
    unsigned endDom = 2*begDom;
    for (iDom = begDom; iDom < 2*begDom; iDom++) {
      for (iLayer = 0; iLayer < iLevel; iLayer++) {
	fprintf(fich,"%d\t",factConnect[iDom-1][iLayer].numberOfStrips());
	for (SetOfStrips::iterator itS = factConnect[iDom-1][iLayer].begin();
	     itS != factConnect[iDom-1][iLayer].end(); itS++) {
	  fprintf(fich, "%d %d\t",(*itS).begin_dst, (*itS).width);
	}
	fprintf(fich,"\n");
      }// for iLayer
    }// for iDom
  }// for iLevel
  fclose(fich);  
# endif  
  return factConnect;
}

void Tree::draw_csr(char *filename, 
		    int dim,
		    const int *ptRows,
		    const int *indCols,
		    const int *ptRow_diag,
		    const int *indVals,
		    const unsigned *l2d)
{
  FILE *fp;
  char color[256]; 
  double xscale, yscale;
  double xpos, ypos;
  int *posDom = new int[_nbDoms];
  xscale = yscale =  380.0 / (double)(dim + 2);
  xpos = 10.0;
  ypos = 390.0;
  if((fp = fopen(filename, "w")) == NULL) {
    exit(-1);
  }
  fprintf(fp, "%%!PS-Adobe-3.0 EPSF-3.0\n%%%%BoundingBox: 5 5 395 395\n");
  fprintf(fp, "/rr { %g } def\n", 
	  0.05 < (xscale * 0.45) ? 0.05 : (xscale *0.45));
  //  fprintf(fp, "/rrr { %g } def\n", 0.5);
  fprintf(fp, "/n { newpath } def\n");
  fprintf(fp, "/c { closepath } def\n");
  fprintf(fp, "/s { 0.01 setlinewidth stroke } def\n");
  fprintf(fp, "/f { fill } def\n");
  fprintf(fp, "/sfill { 0.01 setlinewidth fill } def\n");
  fprintf(fp, "/gs { gsave } def\n");
  fprintf(fp, "/gr { grestore } def\n");
  fprintf(fp, "/rl { rlineto } def\n");
  fprintf(fp, "/m { moveto } def\n");
  fprintf(fp, "/black { 0 0 0 setrgbcolor } def\n");
  fprintf(fp, "/red { 1 0 0 setrgbcolor } def\n");
  fprintf(fp, "/green { 0 1 0 setrgbcolor } def\n");
  fprintf(fp, "/blue { 0 0 1 setrgbcolor } def\n");
  fprintf(fp, "/yellow { 1 1 0 setrgbcolor } def\n");
  fprintf(fp, "/magenta { 1 0 1 setrgbcolor } def\n");
  fprintf(fp, "/cyan { 0 1 1 setrgbcolor } def\n");
  fprintf(fp, "/darkgreen { 0 0.5 0 setrgbcolor } def\n");  
  fprintf(fp, "/gray95 { 0.95  setgray } def\n");
  fprintf(fp, "/graycyan { 0.2 0.6 0.6 setrgbcolor } def\n");
  fprintf(fp, "n 10 10 m 380 0 rl 0 380 rl -380 0 rl 0 -380 rl c gray95 fill\n");
  int pos = 0;
  int dense = 0;
  posDom[0] = 0;
  for (int n = (_nbLevels - 1); n >=0; n--) {
    unsigned begDom = (1 << n);
    unsigned endDom = 2 * begDom; // (1 << (n + 1)); 
    for (unsigned d = begDom; d < endDom; d++) {
      const int d1 = d - 1;
      posDom[d1] = pos;
      pos += _sizeOfDomains[d1];
    }
  }
  // drawing position of non-zero elemnet laying over domain decompositions
  //strcpy(color, "yellow");
  for (int n = (_nbLevels - 1); n >=0; n--) {
    unsigned begDom = (1 << n);
    unsigned endDom = 2 * begDom; // (1 << (n + 1)); 
    for (unsigned d = begDom; d < endDom; d++) {
      const int d1 = d - 1;
      int dd = d;
      if (n == (_nbLevels - 1)) {
	strcpy(color, "blue");
      }
      else {
	strcpy(color, "red");
	dense += _sizeOfDomains[d1] * _sizeOfDomains[d1];
      }
      fprintf(fp, "n %g %g m %g 0 rl 0 %g rl %g 0 rl c gs %s f gr %s s\n",
	      xpos + (double)posDom[d1] * xscale,
	      ypos - (double)(posDom[d1] + _sizeOfDomains[d1]) * yscale,
	      _sizeOfDomains[d1] * xscale,
	      (double)_sizeOfDomains[d1] * yscale,
	      (-1.0) * (double)_sizeOfDomains[d1] * xscale,
	      color, color);
  
      dd = d;
      if (n == (_nbLevels - 1)) {
	strcpy(color, "darkgreen");
      }
      else {
	strcpy(color, "graycyan");
      }
      for (int m = (n - 1); m >= 0; m--) {
	dd /= 2;
	const unsigned dd1 = dd - 1;
	for (SetOfStrips::iterator it = _interco[d1][m].begin();
	     it != _interco[d1][m].end(); ++it) {
	  fprintf(fp, "n %g %g m %g 0 rl 0 %g rl %g 0 rl c gs %s f gr %s s\n",
		  xpos + (double)(posDom[dd1] + (*it).begin_dst) * xscale,
		  ypos - (double)(posDom[d1] + _sizeOfDomains[d1]) * yscale,
		  (double)(*it).width * xscale,
		  (double)_sizeOfDomains[d1] * yscale,
		  (-1.0) * (double)(*it).width * xscale,
		  color, color);
	  if (n < (_nbLevels - 1)) {
	    dense += _sizeOfDomains[d1] * (*it).width;
	  }
	}
      }
    }
  } // loop : n
#if 1    
  int k = 0;
  for (int i = 0; i < dim; i++) {
    for (int kk = ptRow_diag[i]; kk < ptRows[i + 1]; kk++, k++) {
      switch (indVals[k]) {
      case 0 :
	strcpy(color, "red");
	break;
      case 1:
	strcpy(color, "black");
	break;
      case 2:
	strcpy(color, "cyan");
	break;
      default:
	strcpy(color, "green");
	break;
      }
      int j = indCols[kk];
#if 1
      int ii = _global2local[i];
      int jj = _global2local[j];
      int di = l2d[ii];
      int dj = l2d[jj];
#else
      unsigned ii, jj, di, dj;
      getGlob2Loc_dom(i, ii, di);
      getGlob2Loc_dom(j, jj, dj);
#endif
      if (dj > di) {
	int itmp;
	itmp = ii;
	ii = jj;
	jj = itmp;
	itmp = di;
	di = dj;
	dj = itmp;
	if (indVals[k] == 0) {
	  strcpy(color, "magenta");
	}
      }

      fprintf(fp,"n %g %g rr 0 360 arc %s sfill\n", 
	      xpos + 
	      ((double)(jj - _ptOnDomains[dj - 1] + posDom[dj - 1]) + 0.5) * xscale,
	      ypos - ((double)(ii - _ptOnDomains[di - 1] + posDom[di - 1]) + 0.5) * yscale,
	      color);
      if (dj == di) {
	strcpy(color, "graycyan");
	fprintf(fp,"n %g %g rr 0 360 arc %s sfill\n", 
		xpos + ((double)(ii - _ptOnDomains[di - 1] + posDom[di - 1]) + 0.5) * xscale,
		ypos - ((double)(jj - _ptOnDomains[di - 1] + posDom[di - 1]) + 0.5) * yscale,
		color);
      }
    }
  }
#endif
  fprintf(fp, "%%dense=%d\n", dense);
  fprintf(fp, "showpage\n");
  fclose(fp);
}
