/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : The PCM class
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Dion Crannitch
// Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-universite.fr

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

/*
 * Author: Dion Crannitch (adapted from Dave Mason's complexMap class)
 * Modified : 22 June 1999
 * Last Modified: 20 Oct. 2009  (F. Hecht,  Frederic.Hecht@upmc.fr)   for g++ v4
 *
 * Declarations for PCM class, a data structure containing a 2D vector
 * image, with encapsulated operations for reading and writing PCM files.
 *
 */

#ifndef __PCM__
#define __PCM__

#include <cstdio>
using namespace std;
typedef struct {
  float r, i;
} pcm_complex;

class PCM {
  void CalcMax( );

 public:
  int width, height;
  unsigned long pixels;
  float max;
  pcm_complex *image;
  PCM(int w, int h);
  PCM(const char *filename);
  pcm_complex *Get(int x, int y);
  void Set(int x, int y, pcm_complex c);
  void Load(const char *filename);
  void Save(const char *filename);
  ~PCM( );

 private:
  PCM(const PCM &);
  void operator=(const PCM &);
};

#endif
