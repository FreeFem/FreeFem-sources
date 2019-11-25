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
// SUMMARY : Methods for the PCM class
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Dion Crannitch
// Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-universite.fr

// *INDENT-OFF* //
// ff-c++-LIBRARY-dep:
// ff-c++-cpp-dep:
// *INDENT-ON* //

/*
 * Author: Dion Crannitch (adapted from Dave Mason's complexMap class)
 * Last Modified : 22 June 1999
 * Last Modified: 20 Oct. 2009  (F. Hecht,  Frederic.Hecht@upmc.fr)   for g++ v4
 *  change ENDIANESS interface
 *         add  missing read max value in file
 *         add const in char * declaration  for g++-4.
 *         change 32 in 4  for the size of float .
 */

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include "pcm.hpp"
#include <iostream>
using namespace std;
#define BINARY_IN ios::binary | ios::in
#define BINARY_OUT ios::binary | ios::out
#define TOKEN_SIZE 100

PCM::PCM(int w, int h) {
  width = w;
  height = h;
  pixels = (long)w * h;
  max = 0.0;
  image = new pcm_complex[pixels];
}

PCM::PCM(const char *filename) : image(0) { Load(filename); }

void PCM::CalcMax( ) {
  unsigned long p;
  float i, r, m;

  max = m = 0.0;

  for (p = 0; p < pixels; p++) {
    r = image[p].r;
    i = image[p].i;
    m = r * r + i * i;
    if (m > max) {
      max = m;
    }
  }

  max = sqrt(max);
}

void PCM::Set(int x, int y, pcm_complex c) {
  if ((x >= 0) && (y >= 0) && (x < width) && (y < height)) {
    image[(y * width) + x] = c;
  }
}

pcm_complex *PCM::Get(int x, int y) {
  if ((x >= 0) && (y >= 0) && (x < width) && (y < height)) {
    return &image[(y * width) + x];
  } else {
    return NULL;
  }
}

PCM::~PCM( ) { delete[] image; }

// fatal error function - print error message then exit
void fatal_error(const char *mesg) {
  fprintf(stderr, "%s\nFatal error - exitting\n", mesg);
  exit(0);
}

/* For compatibility between systems, the floats in PCM files are
 * always saved as little endian. Therefore on big endian platforms,
 * they must be converted during load and save procedures.
 */

// change float from BIG_ENDIAN to LITTLE_ENDIAN and vice versa
void swap_float_endian(float *number) {
  if (sizeof(float) != 4) {
    fatal_error("PCM -> sizeof(float) != 4.\nProgram cannot continue - exiting.");
  }

  float in = *number;
  float out = 0.;
  ((char *)&out)[0] = ((char *)&in)[3];
  ((char *)&out)[1] = ((char *)&in)[2];
  ((char *)&out)[2] = ((char *)&in)[1];
  ((char *)&out)[3] = ((char *)&in)[0];
  *number = out;
}

// dummy function
void do_nothing(float *number) {
  // do nothing
}

// get a token from the input stream
void extract_token(ifstream &input_stream, char *token, int max_size) {
  int count = 0;
  char ch;

  // skip any whitespace or comments
  do {
    input_stream.read((char *)&ch, sizeof(char));
    // if we've hit a comment char then read till the next "\n"
    if (ch == '#') {
      while (ch != '\n') {
        input_stream.read((char *)&ch, sizeof(char *));
      }
    }
  } while (ch == ' ' || ch == '\t' || ch == '\n');

  // copy data into 'token'
  do {
    if (count >= max_size - 1) {
      fatal_error("extract_token -> token too large");
    }

    token[count++] = ch;
    input_stream.read((char *)&ch, sizeof(char));
  } while (ch != ' ' && ch != '\t' && ch != '\n' && ch != '.');

  input_stream.putback(ch);
  token[count] = '\0';
}

void PCM::Load(const char *filename) {
  char token[TOKEN_SIZE];
  ifstream input_stream(filename, BINARY_IN);

  if (!input_stream) {
    fatal_error("PCM::Load -> file not found.");
  }

  // check magic number
  extract_token(input_stream, token, TOKEN_SIZE);
  if (strcmp(token, "PC") != 0) {
    fprintf(stderr, "Magic number \"%s\" != PC\n", token);
    fatal_error("PCM::Load -> bad magic number");
  }

  // get dimensions
  extract_token(input_stream, token, TOKEN_SIZE);
  width = atoi(token);
  extract_token(input_stream, token, TOKEN_SIZE);
  height = atoi(token);
  extract_token(input_stream, token, TOKEN_SIZE);
  max = atof(token);
  cout << " pcm : " << width << "x" << height << "  max :" << max << endl;
  // Reallocate memory, if necessary.
  unsigned long p = (long)(width * height);
  if (p != pixels) {
    pixels = p;
    if (image) {
      delete[] image;
      image = NULL;
    }
  }

  if (!image) {
    image = new pcm_complex[pixels];
  }

  // skip max cols
  extract_token(input_stream, token, TOKEN_SIZE);
  char ch;
  input_stream.read((char *)&ch, sizeof(char));

  // if machine is big endian then we need to convert floats from little endian
  void (*endian_filter)(float *);
  int i4 = 1;
  char zero_is_little_endian = *((static_cast< const char * >((void *)&i4)) + 3);

  if (zero_is_little_endian == 0) {
    endian_filter = (void (*)(float *))do_nothing;
  } else {
    endian_filter = (void (*)(float *))swap_float_endian;
  }

  // read data
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      pcm_complex c;
      input_stream.read((char *)(void *)&c.r, sizeof(c.r));
      input_stream.read((char *)(void *)&c.i, sizeof(c.i));
      if (x < 0 && y < 0) {
        cout << x << y << "   " << c.r << " " << c.i << endl;
      }

      endian_filter(&c.r);
      endian_filter(&c.i);
      Set(x, y, c);
    }
  }

  input_stream.close( );
}

void PCM::Save(const char *filename) {
  ofstream output_stream(filename, BINARY_OUT);

  if (!output_stream) {
    fatal_error("PCM::Save -> error creating file.");
  }

  // find the maximum vector magnitude
  CalcMax( );

  char header[100];
  sprintf(header, "PC\n%d %d\n%f\n", width, height, max);
  output_stream.write(header, strlen(header));

  void (*endian_filter)(float *);
  int i4 = 1;
  char zero_is_little_endian = *((static_cast< const char * >((void *)&i4)) + 3);

  if (zero_is_little_endian == 0) {
    endian_filter = (void (*)(float *))do_nothing;
  } else {
    endian_filter = (void (*)(float *))swap_float_endian;
  }

  // write data
  pcm_complex *c;

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      c = Get(x, y);
      if (c) {
        endian_filter(&(c->r));
        endian_filter(&(c->i));
        output_stream.write((const char *)(const void *)&(c->r), sizeof(c->r));
        output_stream.write((const char *)(const void *)&(c->i), sizeof(c->i));
      }
    }
  }

  output_stream.close( );
}
