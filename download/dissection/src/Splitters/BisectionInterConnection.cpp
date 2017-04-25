/*! \file   BisectionInterConnetction.cpp
    \brief  Inter connection between off-diag strips
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
// ==      Implementation of the Interconnection methods       ==
// ==============================================================
# include <cstdlib>
# include <list>
# include <algorithm>
# include <cassert>
// # include <iostream>
# include "BitTools/BitManipulations.hpp"
# include "Compiler/DebugUtils.hpp"
# include "BisectionInterConnection.hpp"
using namespace Dissection;

// ==============================================================
// ==                       Strip methods                      ==
// ==============================================================
Strip::Strip(std::list<unsigned>& indices)
{
  assert(!indices.empty());
  /* Next contiguous index. */
  unsigned nextContInd;
  /* First local index is the first element of the list */
  begin_dst = indices.front();
  /* Which one remove from the list */
  indices.pop_front();
  nextContInd = begin_dst + 1;
  while ((!indices.empty()) && (nextContInd == indices.front()))
  {
    nextContInd += 1;
    indices.pop_front();
  }
  width = nextContInd - begin_dst;
}
// ..............................................................
bool
Strip::inside(const Strip& T, unsigned& decal) const
{
  bool flg = (begin_dst>=T.begin_dst) && 
             (nextIndice()<=T.nextIndice());
  decal = unsigned(T.begin_dst - begin_dst);
  return flg;
}
// --------------------------------------------------------------
#if 0
std::ostream& 
operator << (std::ostream& out, const Dissection::Strip& strp)
{
    out << "Strip(first " << strp.begin_dst << ", size " 
	<< strp.width << " )";
    return out;
}
#endif
// ==============================================================
// ==                     SetOfStrips methods                  ==
// ==============================================================
SetOfStrips::SetOfStrips(std::list<unsigned>& indices) :
    _lst_strips()
{
  unsigned start = 0;
  while (!indices.empty()) {
    _lst_strips.push_back(Strip(indices));
    Strip& strp = _lst_strips.back();
    strp.begin_src = start;
    start += strp.width;
  }
  assert(indices.empty());
}
// ..............................................................
void 
SetOfStrips::push_back(const SetOfStrips& strps)
{
  for (const_iterator it = strps.begin();
       it != strps.end(); ++it) {
    _lst_strips.push_back(*it);
  }
}
// ..............................................................
SetOfStrips&
SetOfStrips::operator += (const SetOfStrips& strips)
{
  std::list<Strip> unionOfStrips;
  const std::list<Strip>&          lstS2= strips._lst_strips;
  std::list<Strip>::iterator       itS1 = _lst_strips.begin();
  std::list<Strip>::const_iterator itS2 = lstS2.begin();
  // Easy cases when current or parameter set of strips
  // is empty
  // --------------------------------------------------
  if (itS2 == lstS2.end()) return *this;
  if (itS1 == _lst_strips.end())
  {
    _lst_strips = lstS2;
    return *this;
  }
  // General case, when neither set is empty :
  // ----------------------------------------
  Strip curStrip;
  unsigned nextCurStripInd;
  if ((*itS1).begin_dst < (*itS2).begin_dst)
  {
    curStrip = *itS1;
    itS1++;
  }
  else
  {
    curStrip = *itS2;
    itS2++;
  }
  nextCurStripInd = curStrip.nextIndice();
  while ((itS1 != _lst_strips.end()) &&
	 (itS2 != lstS2.end()       ))
  {
    // Cases where we bind strips : Update current strip
    if ((*itS1).begin_dst<=nextCurStripInd)
    {
      unsigned nextS1StripInd = (*itS1).nextIndice();
      if (nextS1StripInd>nextCurStripInd) {
	curStrip.width   = nextS1StripInd-curStrip.begin_dst;
	nextCurStripInd = nextS1StripInd;
      }
      itS1++;
    }
    else if ((*itS2).begin_dst<=nextCurStripInd)
    {
      unsigned nextS2StripInd = (*itS2).nextIndice();
      if (nextS2StripInd>nextCurStripInd) {
	curStrip.width   = nextS2StripInd-curStrip.begin_dst;
	nextCurStripInd = nextS2StripInd;
      }
      itS2++;
    }
    else
    { // Case where we can't bind strips. 
      // --------------------------------
      // 1. Register current strip :
      unionOfStrips.push_back(curStrip);
      // Create new strip :
      if ((*itS1).begin_dst < (*itS2).begin_dst)
      {
	curStrip = *itS1;
	itS1++;
      }
      else
      {
	curStrip = *itS2;
	itS2++;
      }
      nextCurStripInd = curStrip.nextIndice();
    }
  }// While (((itS1 != _lst_strips.end())&&(itS2 != lstS2.end()))

  // One of set of strips is cared in the union.
  // Finish the union computation with the remained strips in
  // one of set of strips :
  if (itS1 != _lst_strips.end())
  {
    while ((itS1 != _lst_strips.end()) && ((*itS1).begin_dst<=nextCurStripInd))
    {
      unsigned nextS1StripInd = (*itS1).nextIndice();
      if (nextS1StripInd>nextCurStripInd) {
	curStrip.width   = nextS1StripInd-curStrip.begin_dst;
	nextCurStripInd = nextS1StripInd;
      }
      itS1++;
    }
    unionOfStrips.push_back(curStrip);
    while (itS1 != _lst_strips.end()){
	unionOfStrips.push_back(*itS1);
	itS1++;
    }
  }
  if (itS2 != lstS2.end())
  {
    while ((itS2 != lstS2.end()) &&
	   ((*itS2).begin_dst<=nextCurStripInd))
	   
    {
      unsigned nextS2StripInd = (*itS2).begin_dst+(*itS2).width;
      if (nextS2StripInd>nextCurStripInd) {
	curStrip.width   = nextS2StripInd-curStrip.begin_dst;
	nextCurStripInd = nextS2StripInd;
      }
      itS2++;
    }
    unionOfStrips.push_back(curStrip);
    while (itS2 != lstS2.end()){
      unionOfStrips.push_back(*itS2);
      itS2++;
    }
  }
  // Update src indices in union of strips :
  unsigned start = 0;
  for (itS1 = unionOfStrips.begin(); itS1 != unionOfStrips.end();
       itS1++) {
    (*itS1).begin_src = start;
    start += (*itS1).width;
  }
  // Copy new list of strips in current strip :
  _lst_strips = unionOfStrips;
  return *this;
}
// ............................................................
unsigned
SetOfStrips::numberOfIndices() const
{
  unsigned width = 0U;
  for (const_iterator it = begin(); it != end(); it++)
  {
    width += (*it).width;
  }
  return width;
}
// ............................................................
SetOfStrips::operator std::list<unsigned>()
{
  std::list<unsigned> lst;
  for (SetOfStrips::iterator itS = begin(); itS != end(); itS++)
  {
    for (unsigned l = (*itS).begin_dst; l < (*itS).nextIndice();l++)
      lst.push_back(l);
  }
  return lst;
}
// ............................................................
void
SetOfStrips::para_union(const SetOfStrips& strips,
			SetOfStrips& paraStrips1,
			SetOfStrips& paraStrips2,
			SetOfStrips& seqStrips) const
{
  const std::list<Strip>&          lstS1= _lst_strips;
  const std::list<Strip>&          lstS2= strips._lst_strips;
  std::list<Strip>::const_iterator itS1 = lstS1.begin();
  std::list<Strip>::const_iterator itS2 = lstS2.begin();
  std::list<Strip>&                lstP1= paraStrips1._lst_strips;
  std::list<Strip>&                lstP2= paraStrips2._lst_strips;
  std::list<Strip>&                lstS = seqStrips._lst_strips;
  // Easy cases when current or parameter set of strips
  // is empty
  // --------------------------------------------------
  if (itS2 == lstS2.end()) {
    paraStrips1 = *this;
    return;
  }
  if (itS1 == lstS1.end())
  {
    paraStrips2 = strips;
    return;
  }
  // General case : neither set is empty
  // -----------------------------------
  bool  curStripSeq;// Parallel(false) or sequential(true) strip ?
  unsigned parFromStrip; // Number from which is build the strip
  unsigned begCurStrip, szCurStrip, endCurStrip;

  begCurStrip=std::min((*itS1).begin_dst,(*itS2).begin_dst);
  if ((*itS1).begin_dst == (*itS2).begin_dst)
  {// First strip is sequential in this case :
    curStripSeq = true;    
    if ((*itS1).width < (*itS2).width) {
      szCurStrip  = (*itS1).width;
      itS1++;
    }
    else if ((*itS2).width < (*itS1).width) {
      szCurStrip  = (*itS2).width;
      itS2++;
    }
    else {
      szCurStrip  = (*itS1).width;
      itS1++;
      itS2++;
    }
    lstS.push_back(Strip(begCurStrip,szCurStrip));
  }
  else {// First strip is parallel
    curStripSeq = false;
    if ((*itS1).begin_dst < (*itS2).begin_dst)
    {
      unsigned endS1 = (*itS1).nextIndice()-1;
      parFromStrip   = 1;
      if (endS1 < (*itS2).begin_dst) {
	szCurStrip = (*itS1).width;
	itS1++;
      }
      else
	szCurStrip = (*itS2).begin_dst - begCurStrip;
    }
    else
    {
      unsigned endS2 = (*itS2).nextIndice()-1;
      parFromStrip   = 2;
      if (endS2 < (*itS1).begin_dst) {
	szCurStrip = (*itS2).width;
	itS2++;
      }
      else
	szCurStrip = (*itS1).begin_dst - begCurStrip;
    }
    if (1 == parFromStrip)
      lstP1.push_back(Strip(begCurStrip,szCurStrip));
    else {
      assert(2 == parFromStrip);
      lstP2.push_back(Strip(begCurStrip,szCurStrip));
    }
  }
  endCurStrip = begCurStrip + szCurStrip - 1;
  while ((itS1 != lstS1.end()) && (itS2 != lstS2.end()))
  {
    // Search first indice for the next strip and if it's
    // sequential or not
    if ( ((*itS1).begin_dst <= endCurStrip) ||
	 ((*itS2).begin_dst <= endCurStrip) ) {
      // We are in this case :
      //  |---------|        <---- itS1 or itS2
      //       |-----------| <---- itS2 or itS1
      //  |----|             <---- The current union strip
      // or this case :
      //  |---------|        <---- itS1 or itS2
      //  |----|             <---- itS2 or itS1
      //  |----|             <---- The current union strip
      begCurStrip = endCurStrip+1;
      // In this case, we switch between parallel and sequential
      // strips
      curStripSeq = !curStripSeq;
    }
    else
    {
      // We are in this case :
      //      |-----|        <---- itS1 or itS2
      //         |---------| <---- itS2 or itS1
      // |--|                <---- The current union strip

      begCurStrip=
	  std::min((*itS1).begin_dst,(*itS2).begin_dst);
      // In this case, the next strip is sequential only
      // if both strips to union begin with same indices
      curStripSeq=((*itS1).begin_dst==(*itS2).begin_dst);
    }
    // Search size of the next strip :
    if (curStripSeq)
    {// Strips S1 and S2 have common indices
      unsigned endS1 = (*itS1).nextIndice()-1;
      unsigned endS2 = (*itS2).nextIndice()-1;
      if (endS1 < endS2) {
	szCurStrip = endS1-begCurStrip+1;
	itS1++;
      }
      else if (endS2 < endS1) {
	szCurStrip  = endS2-begCurStrip+1;
	itS2++;
      }
      else {
	szCurStrip  = endS1-begCurStrip+1;
	itS1++;
	itS2++;
      }
      lstS.push_back(Strip(begCurStrip,szCurStrip));
    }
    else {// We care about some part of S1 or S2
	  // having exclusive indices
      if ((*itS1).begin_dst <= begCurStrip)
      {// S1 has some exclusive indices
        unsigned endS1 = (*itS1).nextIndice()-1;
	parFromStrip   = 1;
        if (endS1 < (*itS2).begin_dst) {
	  szCurStrip = endS1 - begCurStrip + 1;
	  itS1++;
        }
        else
	  szCurStrip = (*itS2).begin_dst - begCurStrip;
      }
      else {
	// S2 has some exclusive indices
        unsigned endS2 = (*itS2).nextIndice()-1;
	parFromStrip   = 2;
        if (endS2 < (*itS1).begin_dst) {
	  szCurStrip = endS2 - begCurStrip + 1;
	  itS2++;
        }
        else
	  szCurStrip = (*itS1).begin_dst - begCurStrip;
      }
      if (1 == parFromStrip)
        lstP1.push_back(Strip(begCurStrip,szCurStrip));
      else {
	assert(2 == parFromStrip);
	lstP2.push_back(Strip(begCurStrip,szCurStrip));
      }
    }
    endCurStrip = begCurStrip + szCurStrip - 1;
  }// end while
  // Now, one of S1 or (and) S2 has all indices in union
  if (itS1==lstS1.end())
  {
    if (itS2!=lstS2.end())
    {
     if ((*itS2).begin_dst<endCurStrip) {
      unsigned endS2 = (*itS2).nextIndice()-1;
      lstP2.push_back(Strip(endCurStrip+1,endS2-endCurStrip));
      itS2++;
     }
     while (itS2!=lstS2.end())
       lstP2.push_back(*itS2++);     
    }
  }
  else {
    if ((*itS1).begin_dst<endCurStrip) {
      unsigned endS1 = (*itS1).nextIndice()-1;
      lstP1.push_back(Strip(endCurStrip+1,endS1-endCurStrip));
      itS1++;
    }
    while (itS1!=lstS1.end())
      lstP1.push_back(*itS1++);
  }
  // End of the algorithm....
}
// ............................................................
#if 0
std::ostream& operator << (std::ostream& out, 
			   const Dissection::SetOfStrips& set)
{
    out << "[ ";
    for (Dissection::SetOfStrips::const_iterator itS = set.begin();
	 itS != set.end(); )
    {
	out << (*itS);
	itS++;
	if (itS != set.end()) out << ", ";
    }
    out << " ]";
    return out;
}
#endif
// ==============================================================
// ==          FathersStrips methods implementation          ==
// ==============================================================
FathersStrips::FathersStrips() :
    _level(0),
    _outOfDiags(NULL),
    _seq_strips(NULL),
    _par_strips(NULL),
    _row_strips(NULL),
    _col_strips(NULL)
{}
// ------------------------------------------------------------
FathersStrips::FathersStrips(unsigned lvl) :
    _level(lvl),
    _outOfDiags(NULL),
    _seq_strips(NULL),
    _par_strips(NULL),
    _row_strips(NULL),
    _col_strips(NULL)
{
  if (lvl>0)
  {
    _outOfDiags = new SetOfStrips[lvl];
    _row_strips = new SetOfStrips[((lvl+1)*lvl)/2];
    _col_strips = new SetOfStrips[((lvl+1)*lvl)/2];
    assert(_row_strips!=NULL);
    assert(_col_strips!=NULL);
  }
}
// ------------------------------------------------------------
FathersStrips::FathersStrips(FathersStrips& inter) :
    _level(inter._level),
    _outOfDiags(inter._outOfDiags),
    _seq_strips(inter._seq_strips),
    _par_strips(inter._par_strips),
    _row_strips(inter._row_strips),
    _col_strips(inter._col_strips)
{
  inter._level = 0;
  inter._outOfDiags = NULL;
  inter._seq_strips = NULL;
  inter._par_strips = NULL;
  inter._row_strips = NULL;
  inter._col_strips = NULL;
}
// ------------------------------------------------------------
FathersStrips::~FathersStrips()
{
  if (_outOfDiags) delete [] _outOfDiags;
  if (_seq_strips) delete [] _seq_strips;
  if (_par_strips) delete [] _par_strips;
  if (_row_strips) delete [] _row_strips;
  if (_col_strips) delete [] _col_strips;
}
// ------------------------------------------------------------
void
FathersStrips::alloc(unsigned lvl)
{
  if (_outOfDiags) delete [] _outOfDiags;
  if (_seq_strips) delete [] _seq_strips;
  if (_par_strips) delete [] _par_strips;
  if (_row_strips) delete [] _row_strips;
  if (_col_strips) delete [] _col_strips;
  _level = lvl;
  if (lvl>0)
  {
    _outOfDiags = new SetOfStrips[lvl];
    _row_strips = new SetOfStrips[((lvl+1)*lvl)/2];
    _col_strips = new SetOfStrips[((lvl+1)*lvl)/2];
    assert(_row_strips!=NULL);
    assert(_col_strips!=NULL);
  }
}
// ------------------------------------------------------------
SetOfStrips&
FathersStrips::operator [] (unsigned lvl)
{
  assert(lvl<_level);
  assert(_outOfDiags!=NULL);
  return _outOfDiags[lvl];
}
// ------------------------------------------------------------
const SetOfStrips&
FathersStrips::operator [] (unsigned lvl) const
{
  assert(lvl<_level);
  assert(_outOfDiags!=NULL);
  return _outOfDiags[lvl];
}
// ____________________________________________________________
SetOfStrips& 
FathersStrips::getRowSetOfStrips(unsigned l1, unsigned l2)
{
  unsigned ind;

  assert(l1<_level);
  assert(l2<_level);

  if (l1 > l2)
  {
    ind = l2+((l1+1)*l1)/2;
    return _col_strips[ind];
  }
  ind = l1+((l2+1)*l2)/2;
  return _row_strips[ind];
}
// ------------------------------------------------------------
const SetOfStrips& 
FathersStrips::getRowSetOfStrips(unsigned l1, unsigned l2) const
{
  unsigned ind;

  assert(l1<_level);
  assert(l2<_level);

  if (l1 > l2)
  {
    ind = l2+((l1+1)*l1)/2;
    return _col_strips[ind];
  }
  ind = l1+((l2+1)*l2)/2;
  return _row_strips[ind];
}
// ------------------------------------------------------------
SetOfStrips& 
FathersStrips::getColSetOfStrips(unsigned l1, unsigned l2)
{
  unsigned ind;

  assert(l1<_level);
  assert(l2<_level);

  if (l1 > l2)
  {
    ind = l2+((l1+1)*l1)/2;
    return _row_strips[ind];
  }
  ind = l1+((l2+1)*l2)/2;
  return _col_strips[ind];
}
// ------------------------------------------------------------
const SetOfStrips& 
FathersStrips::getColSetOfStrips(unsigned l1, unsigned l2) const
{
  unsigned ind;

  assert(l1<_level);
  assert(l2<_level);

  if (l1 > l2)
  {
    ind = l2+((l1+1)*l1)/2;
    return _row_strips[ind];
  }
  ind = l1+((l2+1)*l2)/2;
  return _col_strips[ind];
}
// --------------------------------------------------------------
