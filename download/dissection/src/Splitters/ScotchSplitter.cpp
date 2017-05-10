/*! \file   ScotchSplitter.cpp
    \brief  to call grpah decomposer : SCOTCH
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
// ==     Implementation of splitter using Scotch library      ==
// ==============================================================
#include <stdexcept>
#include <cstdio>
#include <cstring>
#include <algorithm>
extern"C" {
#include <scotch.h>
}

#include "Compiler/DebugUtils.hpp"
#include "BitTools/BitManipulations.hpp"
#include "ScotchSplitter.hpp"

#define TOO_SMALL 16

static unsigned
compLevelOfDoms(int nbDoms, int maxLevels,
		const int* treetab, int* levels,
		int*& nbDomsPerLevel)
{
  int iDebScan, iEndScan, iCurScan;
  //  DBG_PRINT("Number of levels : %d\n",maxLevels);
  //  DBG_PRINT("Number of doms   : %d\n",nbDoms   );
  // Search the effective number of levels given by Scotch :
  int nbLvls = 0;
  for (unsigned i = 0; i < nbDoms; i++) {
    int locLvl = 1;
    int indDom = i;
    while (indDom != nbDoms-1) {
      locLvl ++;
      indDom = treetab[indDom];
    }
    nbLvls = std::max(nbLvls, locLvl);
  }
  //  DBG_PRINT("Effective level : %d\n",nbLvls);
  nbDomsPerLevel = new int[nbLvls];
  for (int i = 0; i < nbDoms-1; i++)
      levels[i] = -1;
  assert(treetab[nbDoms-1]==-1);
  levels[nbDoms-1] = 0;
  iDebScan = 0;
  while (iDebScan < nbDoms) {
    // Search next domains without level of bissection computation :
    for ( ;(iDebScan<nbDoms)&&(levels[iDebScan] != -1);iDebScan++);
    if (iDebScan<nbDoms) {
      assert(iDebScan < nbDoms);
      assert(iDebScan >=0);
      // Run child to father into the tree to set levels on some domains :
      // We assume than this domain is on a leaf of the tree (uppest layer)
      levels[iDebScan] = nbLvls-1;//maxLevels-1;
      // From this domain, we go through all ancestors until
      // we find a ancestor with defined level (we known here than
      // the root interface has his level defined at 0)
      iCurScan = iDebScan;
      iEndScan = treetab[iDebScan];
      while (-1==levels[iEndScan]) {
	assert(iEndScan<nbDoms);
        levels[iEndScan] = levels[iCurScan]-1;
        iCurScan = iEndScan;// Current is current ancestor now
        iEndScan = treetab[iCurScan];// <-- next ancestor
      }
      if (levels[iCurScan]!=levels[iEndScan]+1) { // Ooops, wrong level for the first domain
	int diff = levels[iCurScan]-levels[iEndScan]-1;// Correction to do
#if 0 // verbose
	if (diff < 0) {
	  std::cerr << "iCurScan = " << iCurScan << ", iEndScan = " << iEndScan << std::endl;
	  std::cerr << "levels[" << iCurScan << "] = " << levels[iCurScan]
		    << ", levels[" << iEndScan << "] = " << levels[iEndScan] << std::endl;
	}
#endif
	assert(diff > 0);
	iCurScan = iDebScan;
	while (iCurScan != iEndScan) {
	    levels[iCurScan] -= diff;
	    iCurScan = treetab[iCurScan];
	}
      }
    }// if iDebScan<nbDoms
    // ok, next domain, please ?
  }
  // Count number of domains per levels :
  memset(nbDomsPerLevel,0,nbLvls*sizeof(int));
  for (int i = 0; i < nbDoms; i++) {
    assert(levels[i] >= 0);
#if 0 // verbose
    if (levels[i]>=nbLvls) std::cerr << "levels[" << i << "] = " << levels[i]
					<< std::endl;
#endif
    CHECK((levels[i]<nbLvls),"Wrong levels value ?");
    CHECK((i<nbDoms), "Wrong value for i");
    nbDomsPerLevel[levels[i]] += 1;
  }
  return nbLvls;
}

bool
ScotchSplitter(unsigned dim, 
	       const int* ptRows, const int* indCols, 
	       unsigned& nbMaxLevels, unsigned minSize,
	       int* loc2glob, int* glob2loc, int& nbDoms, 
	       int*& ptOnDomains, int*& sizeOfDomains,
	       bool checkData, const bool verbose, FILE *fp)
{
  int ierr;
  // check consistency between Scotch and Dissection library
  CHECK(sizeof(int) == sizeof(SCOTCH_Num),
	"Incompatible integer representation between Scotch and Dissection");
  if (verbose) {
    int vers, rela, patc;
    SCOTCH_version(&vers, &rela, &patc);
    fprintf(fp, "%s %d : Soctch version : %d.%d.%d\n",
	    __FILE__, __LINE__, vers, rela, patc);
    //      fprintf(fp, "%s %d : 0 : ptRows = %p indCols = %p\n", 
    //	    __FILE__, __LINE__, (void *)ptRows, (void *)indCols);
  }
  // Allocate and initialize the Scotch graph
  SCOTCH_Graph ptGraph;
  ierr = SCOTCH_graphInit(&ptGraph);
  CHECK(ierr==0,"Fail initializing Scotch graph");
  //  TRACE("Building Scotch graph\n");
  ierr = SCOTCH_graphBuild(&ptGraph,
			   (SCOTCH_Num)0, // offset,
			   (SCOTCH_Num)dim,
			   (const SCOTCH_Num*)ptRows,
			   (const SCOTCH_Num*)ptRows+1,
			   NULL, NULL, (SCOTCH_Num)ptRows[dim],
			   (const SCOTCH_Num*)indCols, NULL);
  //  if (verbose) {
  // fprintf(fp, "%s %d : 1 : ptRows = %p indCols = %p\n", 
  //	    __FILE__, __LINE__, (void *)ptRows, (void *)indCols);
  // }
  CHECK(ierr==0,"Scotch graph building failed !");
  if (checkData) {
    //TRACE("Check Scotch graph\n");
    ierr = SCOTCH_graphCheck(&ptGraph);
    if (ierr) {
      if (verbose) {
	fprintf(fp, "Failed the checking of the graph : bad data ?\n");
      }
      SCOTCH_graphFree(&ptGraph);
      return false;
    }
  }
  // Allocate and initialize the strategy wanted for Scotch
  // TRACE("Initialize strategy for splitting\n");
  SCOTCH_Strat ptStrat;
  ierr = SCOTCH_stratInit(&ptStrat);
  CHECK(ierr==0,
	"Failed initializing Scotch strategy structure");

  char *str_Strat = new char[1024];
  int nbLvls = std::min(unsigned(nbMaxLevels),
			highestbit(unsigned(dim/minSize)));
  if (verbose) {
    fprintf(fp, "nbLevels = %d\n", nbLvls);
  }
  sprintf(str_Strat,"c{rat=0.7,cpr=n{sep=/((levl<%d)|(vert>%d))?m{type=h,rat=0.7,vert=100,low=h{pass=10},asc=b{width=3,bnd=f{bal=0.2},org=h{pass=10}f{bal=0.2}}}|m{type=h,rat=0.7,vert=100,low=h{pass=10},asc=b{width=3,bnd=f{bal=0.2},org=h{pass=10}f{bal=0.2}}};,ole=f{cmin=%d,cmax=%d,frat=0.05},ose=s},unc=n{sep=/(levl<%d)?(m{type=h,rat=0.7,vert=100,low=h{pass=10},asc=b{width=3,bnd=f{bal=0.2},org=h{pass=10}f{bal=0.2}}})|m{type=h,rat=0.7,vert=100,low=h{pass=10},asc=b{width=3,bnd=f{bal=0.2},org=h{pass=10}f{bal=0.2}}};,ole=f{cmin=%d,cmax=%d,frat=0.05},ose=s}}",
	  nbLvls-1,2*minSize-1,minSize,dim,nbLvls-1,minSize,dim);
  //  DBG_PRINT("Strategy string : %s\n", str_Strat);
  ierr = SCOTCH_stratGraphOrder(&ptStrat, str_Strat);
  delete [] str_Strat;
  CHECK(ierr==0,
	"Failed build graph ordering strategy for Scotch");  
  // Ordering with nested bisection :
  // TRACE("Split the graph\n");
  int* rangtab = new int[dim+1];
  int* treetab = new int[dim];
  int nbSplitDoms;
  bool repeat = true;
  int lastCompleteLevel;
  int *levels, *nbDomsPerLevels;
  SCOTCH_randomReset();
  while (repeat) {
    ierr = SCOTCH_graphOrder(&ptGraph, &ptStrat,
			     (SCOTCH_Num*)loc2glob, 
			     (SCOTCH_Num*)glob2loc, 
			     (SCOTCH_Num*)&nbSplitDoms,
			     (SCOTCH_Num*)rangtab,
			     (SCOTCH_Num*)treetab);
    if (ierr) {
      fprintf(fp, "Failed reordering sparse matrix graph !\n");
      SCOTCH_stratExit(&ptStrat);
      SCOTCH_graphFree(&ptGraph);
      return false;
    }

    levels = new int[nbSplitDoms];
//    int *nbDomsPerLevels;// = new int[nbLvls];
    unsigned nbLvlsScotch= compLevelOfDoms(nbSplitDoms, nbLvls, treetab, levels,
					   nbDomsPerLevels);
# if 0 //defined(DISSECTION_DEBUG)
    printf("Level of each domains : \n");
    for (int i = 0; i < nbSplitDoms; i++)
      printf("Domain % d : level %d\n",i+1,levels[i]);
    printf("Domains per level :\n");
    for (int i = 0; i < nbLvls; i++) {
      printf("Level %d : %d domains\n",i,nbDomsPerLevels[i]);
    }
# endif
    /** Search last level where number of domains is a power of two
     */
    lastCompleteLevel = 0;
    while ((lastCompleteLevel<nbLvlsScotch) && 
	   ((1<<lastCompleteLevel) == nbDomsPerLevels[lastCompleteLevel]) ) {
      lastCompleteLevel ++;
    }
    lastCompleteLevel = std::min(lastCompleteLevel, nbLvls);

    nbDoms = (1<<lastCompleteLevel)-1;

    // Search where start each domain per bisection level 
    // and compute the size of each subdomain :
    //  int indDom = 0;
    bool flag_size_check = false;
    for (int i = 0; i < nbSplitDoms; i++) {
      int sz = 0;
      while (levels[i]>=lastCompleteLevel) {
	sz += rangtab[i+1]-rangtab[i];
	i++;
      }
      if (sz+rangtab[i+1]-rangtab[i] <= TOO_SMALL) {
	flag_size_check = true;
	break;
      }
      //DBG_PRINT("Domain %d begin at %d\n",indDom+1,begDom);
    } // loop : i
    if (!flag_size_check) {
      repeat = false;
      break;
    }
    else {
      delete [] levels;
      delete [] nbDomsPerLevels;
    }
  } // while (repeat)
  ptOnDomains   = new int[nbDoms+1];
  sizeOfDomains = new int[nbDoms];
  memset(sizeOfDomains, 0, nbDoms*sizeof(int));
  int* indDomPerLevel = new int[lastCompleteLevel+1];
  memset(indDomPerLevel,0,(lastCompleteLevel+1)*sizeof(int));
  int begDom = 0;
  for (int i = 0; i < nbSplitDoms; i++) {
    int sz = 0;
    while (levels[i]>=lastCompleteLevel) {
      sz += rangtab[i+1]-rangtab[i];
      i++;
    }
    int indDom = (1<<levels[i])-1+indDomPerLevel[levels[i]];
    // DBG_PRINT("level %d : current dom = %d\n", levels[i],indDom+1);
    ptOnDomains[indDom] = begDom;
    //DBG_PRINT("Domain %d begin at %d\n",indDom+1,begDom);
    sizeOfDomains[indDom] = sz+rangtab[i+1]-rangtab[i];
    //DBG_PRINT("Domain %d size of %d\n",indDom+1,sizeOfDomains[indDom]);
    begDom += sizeOfDomains[indDom];
    indDomPerLevel[levels[i]] += 1;
  }
  if (verbose) {
    fprintf(fp, "%s %d : indDomPerlevel[lastCompleteLevel] = %d\n",
	    __FILE__, __LINE__, indDomPerLevel[lastCompleteLevel]);
  }
  ptOnDomains[nbDoms] = dim;
  nbMaxLevels = lastCompleteLevel;
  //  if (verbose) {
  //    fprintf(fp, "%s %d : 5 : ptRows = %p indCols = %p\n", 
  //	    __FILE__, __LINE__, (void *)ptRows, (void *)indCols);
  //}
  // Cleaning temporary arrays :
  delete [] indDomPerLevel;
  delete [] nbDomsPerLevels;
  delete [] levels;
  delete [] treetab;
  delete [] rangtab;
  // Cleaning all Scotch structures :
  SCOTCH_stratExit(&ptStrat);
  SCOTCH_graphFree(&ptGraph);
  return true;
}
