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
/* SUMMARY : ...                                                            */
/* LICENSE : LGPLv3                                                         */
/* ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE           */
/* AUTHORS : Pascal Frey                                                    */
/* E-MAIL  : pascal.frey@sorbonne-universite.fr                             */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef _IMAGE_H_
#define _IMAGE_H_

enum imgtyp { DEFAULT = 0, P2, P3, P5, P6, PS, GREY, RGB, RED, GREEN, BLUE, COLOR };

typedef struct {
  int sizeX, sizeY;
  GLubyte *data;
} PPMimage;
typedef PPMimage *pPPMimage;
typedef struct {
  unsigned char idfield_len;
  unsigned char cmap_type;
  unsigned char image_type;
  unsigned char cmap_spec[5];
  unsigned char x_orig[2];
  unsigned char y_orig[2];
  unsigned char width[2];
  unsigned char height[2];
  unsigned char pixel_size;
  unsigned char image_desc;
} TGAheader;

#endif /* _IMAGE_H_ */

#ifdef __cplusplus
}
#endif
