/*! \file   BitManipulations.hpp
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

// ==============================================================
// ==   Some function to manipulate easily bits on integers    ==
// ==============================================================

#ifndef _DISSECTION_BITTOOLS_BITMANIPULATIONS_HPP_
# define _DISSECTION_BITTOOLS_BITMANIPULATIONS_HPP_

/** @brief Return the highest power of 2 which 
           is lesser or equal to x
*/
inline unsigned highestbit(unsigned x) 
{
  if (0==x) return 0;
# if defined(USE_X86_ASM)
  __asm__("bsr{l}\t%0, %0\n\t" : "=r" (x) : "0" (x));
  //x = asm("bsrl %0 %0" : "=r" (x) : "0" (x));
  return 1U<<(x-1);
# else
  x |= x>>1;
  x |= x>>2;
  x |= x>>4;
  x |= x>>8;
  x |= x>>16;
  return x ^ (x>>1); 
# endif
}
// --------------------------------------------------------------
// @brief Return the position of the highest bit of x
inline unsigned highest_one_idx(unsigned x) 
{
# if defined(USE_X86_ASM)
  __asm__("bsr{l}\t%0, %0\n\t" : "=r" (x) : "0" (x));
  return x;
# else
  unsigned r = 0;
  if (x & 0xffff0000U) { x >>= 16; r += 16; }
  if (x & 0x0000ff00U) { x >>= 8 ; r +=  8; }
  if (x & 0x000000f0U) { x >>= 4 ; r +=  4; }
  if (x & 0x0000000cU) { x >>= 2 ; r +=  2; }
  if (x & 0x00000002U) { x >>= 1 ; r +=  1; }
  return r;
# endif
}

#endif
