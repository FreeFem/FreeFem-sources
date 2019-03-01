/*! \file   BisectionInterConnetction.hpp
    \brief  Inter connection between off-diag strips
    \author Xavier Juvigny, ONERA
    \date   Jul.  2nd 2012
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

// ==============================================================
// ==     Bisection interconnection : store the indices of the==
// == nodes of an interface connected with nodes of a child    ==
// == interface (or domain). The indices are stored as local of==
// == the interface and global.                                ==
// ==============================================================
#ifndef _DISSECTION_SPLITTERS_BISECTIONINTERCONNECTION_HPP_
#define _DISSECTION_SPLITTERS_BISECTIONINTERCONNECTION_HPP_
# include <list>

namespace Dissection
{
  struct Strip {
    unsigned begin_src;// Starting indice in source block
    unsigned begin_dst;// Starting indice in target block
    unsigned width;

    ///@name Constructors and destructor
    ///@{
    /** Default constructor. Initialize with unvalid indices
     */
      Strip() :
	begin_dst(unsigned(-1)), width(unsigned(-1))
	  {}
    /** Simple constructor initializing directly members of the instance
     */
    Strip(unsigned first, unsigned sz) :
	begin_src(0), begin_dst(first), width(sz)
	  {}
    /** Simple constructor which initialize members with starting
	  indices for source block
    */
    Strip(unsigned firstSnd, unsigned firstRcv, unsigned sz) :
	begin_src(firstSnd), begin_dst(firstRcv), width(sz)
	  {}
    /** Copy constructor (default behaviour) */
    Strip(const Strip& strp) :
	begin_src(strp.begin_src), begin_dst(strp.begin_dst),
	width(strp.width)
	  {}
    /** From a list of indices, compute the biggest valid strip from
	the first indices of the list.

	  @return The indices include in the strip are removed from the list.
    */
    Strip(std::list<unsigned>& indices);
    /** Return first indice contiguous to the last indice of the strip
     */
    unsigned nextIndice() const
    {
      return begin_dst + width;
    }
    /** Check if a strip S is contained inside another strip T.
	If true, initialize decal with the difference between
	the first indice of T and S.

	If false, decal is undefined
    */
    bool inside(const Strip& T, unsigned& decal) const;
    /** Check if a strip S contains another strip T.
	If true, initialize decal with the difference between
	the first indice of S and T.
	
	If false, decal is undefined
    */
    bool contains(const Strip& T, unsigned& decal) const
    {
      return T.inside(*this,decal);
    }
    /** Check if an indice is contained by the strip */
    bool contains(unsigned ind) const
    {
      return ((ind>=begin_dst)&&(ind<nextIndice()));
    }
  };
  // --------------------------------------------------------------
  class SetOfStrips
  {
  public:
    typedef std::list<Strip>::iterator iterator;
    typedef std::list<Strip>::const_iterator const_iterator;
    ///@name Constructors and destructor
    //@{
    /** @brief Default constructor 

	Build a empty set of strips
    */
    SetOfStrips() : _lst_strips() {}
    /** @brief Build a list of strips from a set of indices 

	At the exit of the constructor, indices will be an empty
	list.
     */
    SetOfStrips(std::list<unsigned>& indices);
    //@}

    ///@name Getters/Setters
    //@{
    /// @brief Iterator on the first strip
    iterator begin() { return _lst_strips.begin(); }
    /// @brief Const iterator on the first strip
    const_iterator begin() const { return _lst_strips.begin(); }
    /// @brief Iterator on the end of the set 
    iterator end() { return _lst_strips.end(); }
    /// @brief Const iterator on the end of the set 
    const_iterator end() const { return _lst_strips.end(); }
    /// @brief Number of indices included into the set of strips
    unsigned numberOfIndices() const;
    /// @brief Return the number of strips
    unsigned numberOfStrips() const {
      return unsigned(_lst_strips.size());
    }
    //@}

    ///@name Operations on SetOfStrips
    //@{
    /** Check if an indice is contained by the strip */
    bool contains(unsigned ind) const
    {
      bool ok = false;
      for (const_iterator it = begin(); it != end(); it++)
	ok |= (*it).contains(ind);
      return ok;
    }
    /** Add a new strip in the SetOfStrips */
    void push_back(const Strip& strp) {
      _lst_strips.push_back(strp);
    }
    /** Add a new SetOfStrips in the SetOfStrips */
    void push_back(const SetOfStrips& strps);

    /** Convert set of strips as list of indices */
    operator std::list<unsigned>();
    /** Compute the union of current strips indices 
	and strips coming from strips parameter.
    */
    SetOfStrips& operator += (const SetOfStrips& strips);
    /** \brief Compute the union of current strips indices
	       with indices of another strips as union of
	       intersections and differences (for parallel
	       managing).

	The first indice in the source is not computed here,
	the main reason is than we renumber after this computation.
    */
    void para_union(const SetOfStrips& strips,
		    SetOfStrips& paraStrips1,
		    SetOfStrips& paraStrips2,
		    SetOfStrips& seqStrips) const;
    //@}
  private:
    std::list<Strip> _lst_strips;
  };
  // --------------------------------------------------------------
  /** @brief Interconnection between a block B and his ancestors
	     in the bisection tree at the time when the block B 
	     must be eliminated by the factorization algorithm.

      During numerical factorization, Schur complements must
	be computed by the LDU factorization. Schur complements
	are computed for ancestors of a domain and added
	on a part of the condensed problem on these ancestors.
	This class provides some helpers to keep which nodes
	indices in each ancestor are changed by the Schur
	complement.

      NB : For the direct father, only one strip is stored
      ---  for the reason that all nodes of the father are connected
           with the domain.
  */
  class FathersStrips
  {
  public:
    /// @name Constructors and destructor
    //@{
    /** @brief Initialize interconnection with empty connections
     */
    FathersStrips();
    /** @brief Initialize interconnection with his 
	       Schur complement for a node at level lvl (>=1).
    */
    FathersStrips(unsigned lvl);
    /** @brief Copy constructor

	Steal pointer of the copied object which is destroyed !
    */
    FathersStrips(FathersStrips& connections);
    /// @brief Destructor
    ~FathersStrips();
    //@}

    void alloc(unsigned lvl);

    /// @name Getters/Setters
    //@{
    /** @brief Return the set of strips for the out of
	       diagonal block coming from ancestor at
	       level lvl.
    */
    SetOfStrips& operator [] (unsigned lvl);
    /** @brief Return the set of strips for the out of
	       diagonal block coming from ancestor at
	       level lvl.
    */
    const SetOfStrips& operator [] (unsigned lvl) const;
    /** @brief Return the set of strips for row indices
	       for interaction block between
	       the ancestor at level l1 and ancestor at
	       level l2 (0 <= l1,l2<lvl)
    */
    SetOfStrips& getRowSetOfStrips(unsigned l1, unsigned l2);
    /** @brief Return the set of strips for row indices
	       for interaction block between
	       the ancestor at level l1 and ancestor at
	       level l2 (0 <= l1,l2<lvl)
    */
    const SetOfStrips& 
    getRowSetOfStrips(unsigned l1, unsigned l2) const;
    /** @brief Return the set of strips for column indices
	       for interaction block between
	       the ancestor at level l1 and ancestor at
	       level l2 (0 <= l1,l2<lvl)
    */
    SetOfStrips& getColSetOfStrips(unsigned l1, unsigned l2);
    /** @brief Return the set of strips for column indices
	       for interaction block between
	       the ancestor at level l1 and ancestor at
	       level l2 (0 <= l1,l2<lvl)
    */
    const SetOfStrips& 
    getColSetOfStrips(unsigned l1, unsigned l2) const;
      
    //@}
  private:
    unsigned _level; // Level of the node
    // Local indices (referring to each ancestors) for out
    // of diagonal block for current domain stored per ancestor
    SetOfStrips* _outOfDiags;
    // Set of strips for the out of diagonal blocks coming from direct
    // father with another ancestor. This strips are splitted
    // into sequential strips (strips within indices are cared for
    // both children) and parallel strips (one only child contains
    // these indices).
    // NB : For the diagonal block coming from the direct parent,
    // ---  none strip is necessary...
    SetOfStrips* _seq_strips,* _par_strips;

    // Set of strips for row and columns on upper triangular part
    // of the Schur complement (swap row and col for lower part)
    // which not contains a direct parent.
    SetOfStrips* _row_strips;
    SetOfStrips* _col_strips;
  };
}
// std::ostream& operator << (std::ostream& out, 
//			   const Dissection::Strip& strp);
//std::ostream& operator << (std::ostream& out, 
//			   const Dissection::SetOfStrips& set);
#endif
