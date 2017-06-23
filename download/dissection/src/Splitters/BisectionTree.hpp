/*! \file   BisectionTree.hpp
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
// ==   Bisection tree : stored by components in 1D arrays    ==
// ==============================================================
#ifndef _DISSECTION_SPLITTERS_BISECTIONTREE_HPP_
#define _DISSECTION_SPLITTERS_BISECTIONTREE_HPP_
#include <cstdio>
// #include "BitTools/BitManipulations.hpp"
#include "BitManipulations.hpp"
#include "BisectionInterConnection.hpp"
#include "Splitter.hpp"
#include "Algebra/CSR_matrix.hpp"

namespace Dissection
{
  /** @brief Nested Bisection tree

      The nested bisection tree is the core of the
      symbolic factorization.
      Nodes are numbered from root,starting at 1, 
      to leaves by layer order.

      By example :      (layer level)
      -----------  1          0
                  / \
                 2   3        1
                /\   /\
               4  5 6  7      2

      Data are stored components by components, ordered as above :
                       [D1,D2,D3,D4,D5,D6,D7]

      Each layer l data begin at 2^l (0 <= l < number of layers)

      A node n is in layer highest_one_idx(n)

      The brother of a node n (node sharing same father) is n^1
      (^ means xor operation as C convention).

      Ancestors are n>>p where 
      0 < p < number of layers - level of domain
   */
  class Tree
  {
  public:
      /**@name Constructors and destructor */
      //@{
      /** @brief Build bisection tree

	  From the sparse matrix CSR structure,
	  build the bisection tree for symbolic factorization
	  computing some block permutations to speed-up data
	  transfert between blocks.

	  @param dim     Dimension of the sparse matrix
	  @param ptRows  Indices pointing on the beginning of 
	                 each row of the sparse matrix.
			 The last value ptRows[dim] provides the
			 number of non zero coefficients of the
			 sparse matrix.
	  @param indCols Column indice of each non zero coefficients
	                 of the sparse matrix. Values must be stored
			 and sorted (increasing) by row
	  @param nbMaxLevel Maximal number of iteration for the
	                    bisection method.
	  @param minSize    Minimal number of nodes per subdomains
	                    in the leaves of the tree.
          @param spltFct The nested bisection library function to use
	                 By default, Scotch is used to do nested
			 bisection
          @param checkData  If true, verify some consistencies in the given sparse
	                    matrix structure. If false, assume than
			    provided sparse matrix structure is
			    consistency.
       */
    Tree(FILE *fp, bool &berr, unsigned dim,
	 const CSR_indirect *csr,
	 const bool isSym, //const bool isLower,
	 const int *remap_eqn, //const int *remap_indcols,
	 unsigned nbMaxLevel=8, unsigned minSize = 120, 
	 splitter spltFct = NULL,
	 bool checkData = false, bool verbose = true);
    
      /** @brief Destructor
       */
    ~Tree();
      //@}

      /** @name Getters and setters */
      //@{
      /// Return index of the brother (i.e having same father) node
      /// of a node nd.
    static int selfIndex(unsigned nd) {
      return (int)(nd - 1);
    }
    static int Index2Node(unsigned nd) {
      return (int)(nd + 1);
    }
    static unsigned brotherIndex(unsigned nd) {
      return (nd ^ 1);
    }
    /// Return index of the father node index of a node nd
    static unsigned fatherIndex(unsigned nd) {
      return (nd / 2);
    }
    /// Return index of nth forerunner of domain nd :
    static unsigned nthfatherIndex(unsigned nd, unsigned nth) {
      return (nd >> nth);
    }
    /// Return index of the first child index for node nd :
    static unsigned childIndex(unsigned nd) {
      return (nd * 2);
    }
      /// Return layer number where lies the domain nd :
    static unsigned nodeLayer(unsigned nd) {
      return highest_one_idx(nd);
    }
    /// Return 
    unsigned NumberOfLevels() { 
      return _nbLevels;
    }
    unsigned NumberOfSubdomains() {
      return _nbDoms;
    }
    // 
    unsigned NumberOfSubdomains(unsigned level) {
      // (_sbLevels - 1) == (void)
      if (level < _nbLevels) {
	return ((1 << (level + 1)) - 1);  // including own level
      }
      else {
	return 0U;
      }
    }

      /// Return number of internal nodes of a domain nd (>=1)
    unsigned sizeOfDomain(unsigned nd) const {
      return _sizeOfDomains[selfIndex(nd)];
    }

    unsigned sizeOfFathersStrips(unsigned nd) const {
      return _sizeOfFathersStrips[selfIndex(nd)];
    }

    /// Return connection between a domain nd (nd >= 1) and his ancestors
    const FathersStrips& getFathersStrips(unsigned nd) {
      return _interco[selfIndex(nd)];
    }
    const CSR_indirect& getDiagCSR(unsigned nd) {
      return _csr_diag[selfIndex(nd)];
    }
    const CSR_indirect& getOffdiagCSR(unsigned nd) {
      return _csr_offdiag[selfIndex(nd)];
    }
    int *getDiagLoc2Glob(unsigned nd) {
      return _loc2glob_diag[selfIndex(nd)];
    }
    int *getOffdiagLoc2Glob(unsigned nd) {
      return _loc2glob_offdiag[selfIndex(nd)];
    }

    bool save(FILE* stream) const;
      /** @brief Load bisection tree contents */
    bool load(FILE* stream)      ;
      /** @brief Print information on bisection tree 

	  A verbose level can be provided. Depending of
	  his value, informations printed are :
	  
	  >=1 : Some statistics data (number of domains,
	        mean number of nodes per domains per layer,
		variance of number of nodes per domains,...)
	  >=2 : Global indices of nodes contained per domains
	  
	  >=3 : Block permutation from a domain and his forerunnings

	  @param verboseLevel Choose the level of verbose
	                      to print some information for
			      dissection tree.
       */
    void printInfo(FILE* stream, 
		   unsigned verboseLevel = 1) const;
    
      //@}
  private:
    bool _verbose;
      /* Allocate and fill an array with indices of domain 
	 containing each local index.
	 The caller of this method must free the returned array.
      */
    unsigned* compLoc2Dom(unsigned dim) const;
      /* Remove self references from the skeleton graph of the sparse matrix
	 provided to the splitter
      */
    void removeLoops(int dim, 
		     const int* ptRows, const int* indCols,
		     int*& ptRows2, int*& indCols2);
      /* Do symbolic factorization */
    FathersStrips* symbolicFactorization(SetOfStrips** connections);
      /* Renumbering interfaces to optimize number of strips */
    void renumberingInterface(unsigned nperm, unsigned* perm,
			      const SetOfStrips& paral, 
			      const SetOfStrips& parar,
			      const SetOfStrips& seq,
			      unsigned layer, unsigned iDom) const;
    void draw_csr(char *filename, 
		  int dim,
		  const int *ptRows,
		  const int *indCols,
		  const int *ptRow_diag,
		  const int *indVals,
		  const unsigned *l2d);

      /// Number of levels done in nested bisection algorithm

    unsigned _nbDoms;
    unsigned _nbLevels;
      /// Number of nodes (internal only for each leaf)
    int *_sizeOfDomains;
    int *_sizeOfFathersStrips;
      /// Indice (starting at 0) of the beginning of each domain 
      /// in local 2 global correspondance.
    int *_ptOnDomains;
      /// Local 2 global array for each subdomain (partition)
    int *_local2global;
      /// global 2 local array (partition)
    int *_global2local;
      ///
    FathersStrips* _interco;
      /// For global indices, store the domain containing it :
    int* _glob2dom;
    int **_loc2glob_diag;
    int **_loc2glob_offdiag;
    bool _isSym;
    CSR_indirect *_csr_diag;
    CSR_indirect *_csr_offdiag;
  };
}

#endif
