/*! \file   MetisSplitter.cpp
    \brief  to call grpah decomposer : METIS
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

// ============================================================
// ==     Implementation of splitter using Metis library     ==
// ============================================================
#include <stdexcept>
#include <cstdio>
#include <cstring>
#include <algorithm>
#ifndef NO_METIS
#include "metis.h"
#endif

#include "Compiler/DebugUtils.hpp"
#include "BitTools/BitManipulations.hpp"
#include "MetisSplitter.hpp"

static unsigned compBegOfDomains(unsigned invLevel,
		 	         unsigned&begDom,
			         unsigned indDom,
			         const int* sizeOfDomains,
			         int* ptOnDomains)
{
    if (invLevel==1) {
      ptOnDomains[indDom-1] = begDom;
      return sizeOfDomains[indDom-1];
    }
    else {
	begDom += compBegOfDomains(invLevel-1, begDom,
				   2*indDom, sizeOfDomains,
				   ptOnDomains);
	begDom += compBegOfDomains(invLevel-1, begDom,
				   2*indDom+1, sizeOfDomains,
				   ptOnDomains);
	ptOnDomains[indDom-1] = begDom;
	return sizeOfDomains[indDom-1];
    }
}
			     

bool
MetisSplitter(unsigned dim, 
	      const int* ptRows, const int* indCols, 
	      unsigned& nbMaxLevels, unsigned minSize,
	      int* loc2glob, int* glob2loc, int& nbDoms, 
	      int*& ptOnDomains, int*& sizeOfDomains,
	      bool checkData, const bool verbose, FILE *fp)
{
#ifdef NO_METIS
  if (verbose) {
    fprintf(stderr, "%s %d : Metis is not linked\n", __FILE__, __LINE__);
  }
  return false;
#else
  /** Check inputs */
  CHECK(ptRows  !=NULL, "Null pointer for ptRows  !");
  CHECK(indCols !=NULL, "Null pointer for indCols !");
  CHECK(loc2glob!=NULL, "Null pointer for loc2glob!");
  CHECK(glob2loc!=NULL, "Null pointer for glob2loc!");
  int ierr;
  // Compute nbMaxLevels according to the minSize parameter
  int maxBlocks = dim/minSize;
  unsigned nbLevels = highestbit(maxBlocks)-1;
  nbMaxLevels = (nbMaxLevels<nbLevels ? nbMaxLevels : nbLevels);
  nbDoms = 1<<(nbMaxLevels-1);
  if (dim==0) return true;
  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_PTYPE    ] = METIS_PTYPE_RB;
  options[METIS_OPTION_OBJTYPE  ] = METIS_OBJTYPE_NODE;
  options[METIS_OPTION_IPTYPE   ] = METIS_IPTYPE_NODE;
  options[METIS_OPTION_CTYPE    ] = METIS_CTYPE_RM;
  options[METIS_OPTION_RTYPE    ] = METIS_RTYPE_SEP1SIDED;
  options[METIS_OPTION_NSEPS    ] = 1;/*Default value */
  options[METIS_OPTION_NITER    ] = 10;/* Default value */
  options[METIS_OPTION_UFACTOR  ] = 1;
  options[METIS_OPTION_COMPRESS ] = 0;/* No compress */
  options[METIS_OPTION_CCORDER  ] = 0;/* No connect component detection */
  options[METIS_OPTION_SEED     ] = -1;/* Seed of random algo */
  options[METIS_OPTION_PFACTOR  ] = 0;/* No dense nodes removal */
  options[METIS_OPTION_NUMBERING] = (idx_t)0; //offset;
# ifdef DEBUG
  options[METIS_OPTION_DBGLVL   ] = 255;/* Full debug */  
# else
  options[METIS_OPTION_DBGLVL   ] = METIS_DBG_TIME; /* only time profiling */  
# endif
//  Opt [0] = 0 ;	/* use defaults */
//  Opt [1] = 3 ;	/* matching type */
//  Opt [2] = 1 ;	/* init. partitioning algo*/
//  Opt [3] = 2 ;	/* refinement algorithm */
//  Opt [4] = 0 ;	/* no debug */
//  Opt [5] = 0 ;	/* unused */
//  Opt [6] = 0 ;	/* unused */
//  Opt [7] = -1 ;	/* random seed */
  // 
  idx_t nvtxs = dim;
  idx_t *xadj, *adjncy, *perm, *iperm;
  if (sizeof(idx_t)==sizeof(int)) {// Ok, simple pointer may be ok
      xadj    = (idx_t*)ptRows;
      adjncy  = (idx_t*)indCols;
      perm    = (idx_t*)loc2glob;
      iperm   = (idx_t*)glob2loc;
  }
  else {// Arf, some copy to do :(
      xadj = new idx_t[dim+1];
      for (unsigned ix = 0; ix <= dim; ix++)
	  xadj[ix] = ptRows[ix];
      adjncy = new idx_t[ptRows[dim]];
      for (unsigned ix = 0; ix < ptRows[dim]; ++ix)
	  adjncy[ix] = indCols[ix];
      perm  = new idx_t[dim];
      iperm = new idx_t[dim];
  }
  if (verbose) {
    fprintf(fp, "%s %d : nbDoms = %d\n",
	    __FILE__, __LINE__, nbDoms);
  }
  idx_t* sizes = new idx_t[2*nbDoms-1];
  /** This metis subroutine make dissection with some history
      in sizes array
  */
#if 1
  int err = METIS_NodeNDP(nvtxs, xadj, adjncy, NULL, 
			  idx_t(nbDoms), options, 
			  iperm, perm, sizes);
  if (verbose && (err != METIS_OK)) {
      fprintf(fp, "Failed to complete Metis bisection !\n");
      if (err == METIS_ERROR_INPUT)
	fprintf(fp, "\t Wrong arguments given to Metis\n");
      else if (err == METIS_ERROR_MEMORY)
	fprintf(fp, "\t Metis fails allocate some memory\n");
      else
	fprintf(fp, "\t Other error in Metis\n");
  }
#if 0
  FILE *fpp;
  fpp = fopen("metis.data", "w");
  fprintf(fpp, "%d %d\n", dim, nbDoms);
  for (int i = 0; i < dim; i++) {
    fprintf(fpp, "%d\n", perm[i]);
  }
  for (int i = 0; i < dim; i++) {
    fprintf(fpp, "%d\n", iperm[i]);
  }
  for (int i = 0; i < (2 * nbDoms - 1) ; i++) {
    fprintf(fpp, "%d\n", sizes[i]);
  }
  fclose(fpp);
#endif
#else
  FILE *fpp;
  int itmp, jtmp;
  fpp = fopen("metis.data", "r");
  fscanf(fpp, "%d %d", &itmp, &jtmp);
  for (int i = 0; i < dim; i++) {
    fscanf(fpp, "%d", &perm[i]);
  }
  for (int i = 0; i < dim; i++) {
    fscanf(fpp, "%d", &iperm[i]);
  }
  for (int i = 0; i < (2 * nbDoms - 1) ; i++) {
    fscanf(fpp, "%d", &sizes[i]);
  }
  fclose(fpp);
#endif
  //  TRACE("Finish to split ;)\n");

  nbDoms = nbDoms*2-1;
  // Compute where start each domain per bisection level
  // and compute the size of each subdomain
  ptOnDomains   = new int[nbDoms+1];
  sizeOfDomains = new int[nbDoms];
  unsigned ptSize = 0;
  unsigned begDom = 0;
  // Metis store size of each domains and interface
  // per level of tree, beginning with the leaves of the tree
  for (unsigned iLvl = nbMaxLevels; iLvl > 0; iLvl --)
  {
    for (unsigned iDom = (1<<(iLvl-1)); iDom < (1<<iLvl);
	 ++iDom)
    {
      assert(ptSize < nbDoms);
      sizeOfDomains[iDom-1] = sizes[ptSize];
      ptSize++;
    }
  }
  delete [] sizes;
  // But renumbering is done as left child - right child - father
  // way !
  compBegOfDomains(nbMaxLevels, begDom, 1, sizeOfDomains,ptOnDomains);
  ptOnDomains[nbDoms] = dim;  
# ifdef DEBUG
  for (unsigned idom = 0; idom < nbDoms; idom++)
      std::cerr << "dom " << idom+1 << " => "
		<< "start : " << ptOnDomains[idom] 
		<< " size : " << sizeOfDomains[idom] 
		<< " next : " << ptOnDomains[idom]+sizeOfDomains[idom]
		<< std::endl;
# endif
  //assert(begDom==dim);
  assert(ptOnDomains[0]+sizeOfDomains[0]==dim);
  if (sizeof(idx_t)!=sizeof(int)) {
      for (unsigned i = 0; i < dim; i++)
      {
	loc2glob[i] =  perm[i];
	glob2loc[i] = iperm[i];
      }
      delete [] iperm;
      delete [] perm;
      delete [] adjncy;
      delete [] xadj;
  }
  return true;
#endif
}
