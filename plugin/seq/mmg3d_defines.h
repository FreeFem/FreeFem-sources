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
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

#define EPS 1.e-06
#define EPS1 1.e-9
#define EPS2 1.e-12
#define EPSOK 1.e-18
#define EPS30 1.e-30

#define ALPHAC 0.20412415 /* sqrt(6)/12 */
#define ALPHAD 0.04811252 /* 1.0/(12*sqrt(3)) */
#define BETAC 0.03928371  /* sqrt(2)/36 */

#define LLONG 1.3
#define LSHORT 0.72
#define LFILT 0.7
#define QDEGRAD 2.45

#define LONMAX 4096
#define NPMAX 500000
#define NTMAX 1000000
#define NEMAX 3000000

#define PRECI 1
#define BUCKSIZ 64

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) < (b) ? (b) : (a))

#define M_NOTAG (0)
#define M_UNUSED (1 << 0)
#define M_BDRY (1 << 1)
#define M_MOVE (1 << 2)
#define M_CAVITY (1 << 3)
#define ALL_BDRY 63

#ifdef INT_MAX
#undef INT_MAX
#undef SHORT_MAX
#endif
#define INT_MAX 0x7fffffff
#define SHORT_MAX 0x7fff

extern unsigned char MMG_idir[4][3];
extern unsigned char MMG_inxt[7];
extern unsigned char MMG_iarf[4][3];
extern unsigned char MMG_iare[6][2];
extern unsigned char MMG_ifar[6][2];
extern unsigned char MMG_isar[6][2];
extern unsigned char MMG_arpt[4][3];
