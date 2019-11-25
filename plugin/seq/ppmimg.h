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

#ifndef _PPMIMG_H_
#define _PPMIMG_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef unsigned char ubyte;

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

enum imgtyp { DEFAULT = 0, P2, P3, P4, P5, P6, GREY, RGB, RED, GREEN, BLUE, COLOR };

typedef struct {
  short sizeX, sizeY;
  ubyte type;
  ubyte *data;
} PPMimage;
typedef PPMimage *pPPMimage;

/* prototypes */
PPMimage *loadPPM(const char *imgname, ubyte quiet);
int savePPM(const char *imgname, pPPMimage img);
pPPMimage diffImg(pPPMimage bits, pPPMimage img);
void freePPMimage(pPPMimage &image);

pRnm PPMimage2Rnm(pPPMimage const &img);
pPPMimage Rnm2PPMimage(pRnm const &array);

#ifdef __cplusplus
}
#endif

#endif /*_PPMIMG_H_*/
