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
// SUMMARY : READ/WRITE MESH AND SOLUTION IN FORMAT VTK
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Jacques Morice
// E-MAIL  : jacques.morice@ann.jussieu.fr

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

/*
 * Thank to the ARN ()  FF2A3 grant
 * ref:ANR-07-CIS7-002-01
 */

// FH   July 2009
// comment all
// Th3_t->BuildBound();
// Th3_t->BuildAdj();
// Th3_t->Buildbnormalv();
// Th3_t->BuildjElementConteningVertex();
// is now in the constructor of Mesh3 to be consistante.
// FH  nev 2009
// correct  gestion of nameofuser    variable
// JM :: VTU extension file

#include <iostream>
#include <cfloat>
#include <cmath>
#include <complex>
using namespace std;
/*
 using namespace std;
 #include "error.hpp"
 #include "AFunction.hpp"
 using namespace std;
 #include "rgraph.hpp"
 #include "RNM.hpp"
 #include "fem.hpp"

 #include "FESpacen.hpp"
 #include "FESpace.hpp"

 #include "MatriceCreuse_tpl.hpp"
 #include "MeshPoint.hpp"
 #include "Operator.hpp"
 #include "lex.hpp"

 #include "lgfem.hpp"
 #include "lgmesh3.hpp"
 #include "lgsolver.hpp"
 #include "problem.hpp"
 // #include "LayerMesh.hpp"
 // #include "TransfoMesh_v2.hpp"
 #include "msh3.hpp"
 // #include "GQuadTree.hpp"
 // #include "lex.hpp"
 */
#include "ff++.hpp"
#include <set>
#include <vector>
#include <list>
#include <fstream>

using namespace Fem2D;

static const char *EncodeB64_LoopByte =
  "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
  "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
  "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
  "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/ ";
static const char *EncodeB64_Mul4Byte =
  "AAAABBBBCCCCDDDDEEEEFFFFGGGGHHHHIIIIJJJJKKKKLLLLMMMMNNNNOOOOPPPPQQQQRRRRSSSSTTTTUUUUVVVVWWWWXXXX"
  "YYYYZZZZaaaabbbbccccddddeeeeffffgggghhhhiiiijjjjkkkkllllmmmmnnnnooooppppqqqqrrrrssssttttuuuuvvvv"
  "wwwwxxxxyyyyzzzz0000111122223333444455556666777788889999++++//// ";

// =====================
// Test for BigEndian
// =====================
bool isBigEndian( ) {
  const int i = 1;
  char *c = (char *)&i;

  if (c[sizeof(int) - 1] == 0) {
    return false;
  } else {
    return true;
  }
}

void encodeB64_3Bytes(unsigned char *in3Bytes, unsigned char *out4Bytes) {
  if (in3Bytes == NULL || out4Bytes == NULL) {
    return;
  }

  out4Bytes[0] = EncodeB64_Mul4Byte[in3Bytes[0]];
  out4Bytes[1] = EncodeB64_LoopByte[((in3Bytes[0] & 0x03) << 4) + ((in3Bytes[1] & 0xf0) >> 4)];
  out4Bytes[2] = EncodeB64_LoopByte[((in3Bytes[1] & 0x0f) << 2) + ((in3Bytes[2] & 0xc0) >> 6)];
  out4Bytes[3] = EncodeB64_LoopByte[in3Bytes[2]];
}

// ===========================
//
// ===========================
int encodeB64(int n, unsigned char *inBytes, unsigned char *outBytes) {
  if (inBytes == NULL || outBytes == NULL || n <= 0) {
    return 0;
  }

  int m = n - (n % 3);
  int ii, jj;
  jj = 0;

  for (ii = 0; ii < m; ii += 3) {
    encodeB64_3Bytes(inBytes + ii, outBytes + jj);
    jj = jj + 4;
  }

  if (m != n) {
    unsigned char lastBytes[3] = {0, 0, 0};
    // give last (n-m) chars
    lastBytes[0] = inBytes[ii];
    if ((n - m) == 2) {
      lastBytes[1] = inBytes[ii + 1];
    }

    // encoding chars
    encodeB64_3Bytes(lastBytes, outBytes + jj);
    // definition of bytes which is not encoding
    outBytes[jj + 3] = '=';
    if ((n - m) == 1) {
      outBytes[jj + 2] = '=';
    }

    jj = jj + 4;
  }

  return jj;
}

int runEncodeB64(int n, unsigned char *inBytes, unsigned char *outBytes) {
  static int nbcached = 0;
  static unsigned char charCache[3];
  int l = 0;
  int nn = n;
  int nbcachedInBytes = 0;

  if (n == 0) {
    // no bytes are given :: the bytes cached is encoding
    if (nbcached > 0) {
      l = encodeB64(nbcached, charCache, outBytes);
      nbcached = 0;
    }
  } else {
    // bytes are given ::
    if (nbcached > 0) {
      // give chars to the Cache
      if (nn > 0) {
        charCache[nbcached++] = inBytes[nbcachedInBytes];
        nbcachedInBytes++;
        nn--;
      }

      if (nbcached < 3 && nn > 0) {
        charCache[nbcached++] = inBytes[nbcachedInBytes];
        nbcachedInBytes++;
        nn--;
      }

      if (nbcached == 3) {
        // the cache is filled :: encode bytes
        l = encodeB64(nbcached, charCache, outBytes);
        outBytes += l;
        nbcached = 0;
      }
    }

    if (nn == 0) {
      return l;
    }

    unsigned char *newInBytes = &inBytes[nbcachedInBytes];
    // unsigned char *newOutBytes = &outBytes[l];
    int m = nn - nn % 3;
    if (nn == m) {
      l += encodeB64(nn, newInBytes, outBytes);
      return l;
    }

    // cache left overs
    charCache[nbcached++] = newInBytes[m];
    if (m + 1 < nn) {
      charCache[nbcached++] = newInBytes[m + 1];
    }

    l += encodeB64(m, newInBytes, outBytes);
  }

  return l;
}

char *newcopy(const char *s) {
  char *r(new char[strlen(s) + 1]);

  strcpy(r, s);
  return r;
}

char *newcopy(const string *s) {
  char *r(new char[s->size( ) + 1]);

  strcpy(r, s->c_str( ));
  return r;
}

// Tables of element type of VTK considered in Freefem++
static const int nvElemVTK[25] = {1, 0, 2, 0, 3, 0, 0, 0, 0, 4, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

// we considerer only vertex, edges, triangles and tetrahedrons in Freefem++
// 1  :: Vertex Corner
// 3  :: Edge/line
// 5  :: triangles
// 10  :: tetrahedrons
enum FFppCells { VTK_VERTEX = 1, VTK_EDGE = 3, VTK_TRI = 5, VTK_TET = 10 };
static const int NbColorTable = 30;
// Table of colors for Labels of elements :: RGB
static const float ColorTable[30][3] = {
  {1.0, 0.0, 0.0}, /* red    */
  {1.0, 1.0, 0.0}, /* yellow */
  {1.0, 0.0, 1.0}, /* ????   */
  {0.0, 1.0, 0.0}, /* green  */
  {0.0, 1.0, 1.0}, /* cyan   */
  {0.0, 0.0, 1.0}, /* blue   */
  {1.0, 0.5, 0.0}, /* orange */
  {0.5, 0.0, 1.0}, /* violet */
  {1.0, 0.0, 0.5}, /* ???  */
  {0.5, 1.0, 0.0}, /* ???  */
  {0.0, 1.0, 0.5}, /* ???  */
  {0.0, 0.5, 1.0}, /* ???  */
  {1.0, 0.5, 0.5}, /*  ???    */
  {1.0, 1.0, 0.5}, /*  ???    */
  {1.0, 0.5, 1.0}, /*  ???    */
  {0.5, 1.0, 0.5}, /*  ???    */
  {0.5, 1.0, 1.0}, /*  ???   */
  {0.5, 0.5, 1.0}, /*  ???   */
  {0.4, 0.0, 0.0}, /* dark blue    */
  {0.4, 0.4, 0.0}, /* dark yellow */
  {0.4, 0.0, 0.4}, /* dark ????   */
  {0.0, 0.4, 0.0}, /* dark green  */
  {0.0, 0.4, 0.4}, /* dark cyan   */
  {0.0, 0.0, 0.4}, /* dark blue   */
  {0.8, 0.0, 0.0}, /*  ???    */
  {0.8, 0.8, 0.0}, /*  ???    */
  {0.8, 0.0, 0.8}, /*  ???    */
  {0.0, 0.8, 0.0}, /*  ???    */
  {0.0, 0.8, 0.8}, /*  ???    */
  {0.0, 0.0, 0.8}, /*  ???    */
};                 // a voir F.Hecht

namespace FreeFEM {
void SwapBytes(char *array, int size, int n) {
  char *x = new char[size];

  for (int i = 0; i < n; i++) {
    char *a = &array[i * size];
    memcpy(x, a, size);

    for (int c = 0; c < size; c++) {
      a[size - 1 - c] = x[c];
    }
  }

  delete[] x;
}
}

//= =============================================
// Fichier Format .vtu
//= =============================================
// general functions

void VTU_BEGIN(FILE *fp) {
  string version = "1.0";

  fprintf(fp, "<?xml version=\"%s\"?>\n", version.c_str( ));
}

void VTU_VTKFILE(FILE *fp, bool bigEndian) {
  string type("UnstructuredGrid");
  string byte_big("BigEndian");
  string byte_little("LittleEndian");
  string version("0.1");

  fprintf(fp, "<VTKFile type=\"%s\"", type.c_str( ));
  fprintf(fp, " version=\"%s\"", version.c_str( ));
  if (bigEndian) {
    fprintf(fp, " byte_order=\"%s\">\n", byte_big.c_str( ));
  } else {
    fprintf(fp, " byte_order=\"%s\">\n", byte_little.c_str( ));
  }

  // fprintf(fp,"<%s>\n",type.c_str());
}

void VTU_PIECE(FILE *fp, const int &nv, const int &nc) {
  fprintf(fp, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nv, nc);
}

void VTU_DATA_ARRAY(FILE *fp, const string &type, const string &name, bool binary) {
  fprintf(fp, "<DataArray type=\"%s\"", type.c_str( ));
  fprintf(fp, " Name=\"%s\"", name.c_str( ));

  if (binary) {
    fprintf(fp, " format=\"binary\"");
  } else {
    fprintf(fp, " format=\"ascii\"");
  }

  fprintf(fp, ">\n");
}

void VTU_DATA_ARRAY(FILE *fp, const string &type, const string &name, const long &noc,
                    bool binary) {
  fprintf(fp, "<DataArray type=\"%s\"", type.c_str( ));
  fprintf(fp, " Name=\"%s\"", name.c_str( ));
  fprintf(fp, " NumberOfComponents=\"%ld\"", noc);

  if (binary) {
    fprintf(fp, " format=\"binary\"");
  } else {
    fprintf(fp, " format=\"ascii\"");
  }

  fprintf(fp, ">\n");
}

void BEGINTYPE_VTU(FILE *fp, string begintype) { fprintf(fp, "<%s>\n", begintype.c_str( )); }

void ENDTYPE_VTU(FILE *fp, string endtype) { fprintf(fp, "</%s>\n", endtype.c_str( )); }
void writebin64(FILE *fp, int nc) {
  unsigned char ElementChars[256];
  unsigned nbytes = nc;
  int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
  ElementChars[l] = 0;
  fwrite(&ElementChars, l, 1, fp);
}
void writebin64flush(FILE *fp) {
  unsigned char ElementChars[256];
  int l = runEncodeB64(0, NULL, ElementChars);
  ElementChars[l] = 0;
  fwrite(&ElementChars, l, 1, fp);
}

void VTU_WRITE_MESH(FILE *fp, const Mesh &Th, bool binary, int datasize, bool surface) {
  int nc, nv, nconnex;

  if (surface) {
    nc = Th.nt + Th.neb;
  } else {
    nc = Th.nt;
  }

  if (surface) {
    nconnex = 3 * Th.nt + 2 * Th.neb;
  } else {
    nconnex = 3 * Th.nt;
  }

  unsigned char ElementChars[256];
  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" ");
  if (isBigEndian( )) {
    fprintf(fp, "byte_order=\"BigEndian\">\n");
  } else {
    fprintf(fp, " byte_order=\"LittleEndian\">\n");
  }

  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\" %d\">\n", Th.nv, nc);

  fprintf(fp, "<Points>\n");
  fprintf(fp, "<DataArray type=\"Float32\" NumberOfComponents=\"3\"");

  // A definir coord

  float *coord = new float[3 * Th.nv];

  for (int ii = 0; ii < Th.nv; ii++) {
    coord[3 * ii] = Th.vertices[ii].x;
    coord[3 * ii + 1] = Th.vertices[ii].y;
    coord[3 * ii + 2] = 0;    // Th.vertices[ii].z;
  }

  if (binary) {
    fprintf(fp, " format=\"binary\">\n	");
    {
      unsigned nbytes = 3 * Th.nv * sizeof(float);
      int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      for (long i = 0; i < Th.nv; i++) {
        l = runEncodeB64(3 * sizeof(float), (unsigned char *)(coord + 3 * i), ElementChars);
        ElementChars[l] = 0;
        fwrite(&ElementChars, l, 1, fp);
      }

      // flush buffer
      l = runEncodeB64(0, NULL, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);
    }
  } else {
    fprintf(fp, " format=\"ascii\">\n");

    for (long i = 0; i < Th.nv; i++) {
      fprintf(fp, "%f %f %f ", coord[i * 3 + 0], coord[i * 3 + 1], coord[i * 3 + 2]);
    }
  }

  fprintf(fp, "\n</DataArray>\n");
  fprintf(fp, "</Points>\n");
  fprintf(fp, "<Cells>\n");

  delete[] coord;
  long *ien = new long[nconnex];

  for (long it = 0; it < Th.nt; it++) {
    const Mesh::Triangle &K(Th.t(it));

    for (int ii = 0; ii < 3; ii++) {
      ien[3 * it + ii] = Th.operator( )(K[ii]);
    }
  }

  if (surface) {
    for (long ibe = 0; ibe < Th.neb; ibe++) {
      const Mesh::BorderElement &K(Th.be(ibe));

      for (int ii = 0; ii < 2; ii++) {
        ien[3 * Th.nt + 2 * ibe + ii] = Th.operator( )(K[ii]);
      }
    }
  }

  fprintf(fp, "<DataArray type=\"Int32\" Name=\"connectivity\" ");

  // need ien vector
  if (binary) {
    fprintf(fp, "format=\"binary\">\n	");
    unsigned nbytes = nconnex * sizeof(int);
    int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
    ElementChars[l] = 0;
    fwrite(&ElementChars, l, 1, fp);

    for (long i = 0; i < Th.nt; i++) {
      l = runEncodeB64(sizeof(int), (unsigned char *)(ien + i * 3), ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      l = runEncodeB64(sizeof(int), (unsigned char *)(ien + i * 3 + 1), ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      l = runEncodeB64(sizeof(int), (unsigned char *)(ien + i * 3 + 2), ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);
    }

    if (surface) {
      for (long i = 0; i < Th.neb; i++) {
        l = runEncodeB64(sizeof(int), (unsigned char *)(ien + 3 * Th.nt + i * 2), ElementChars);
        ElementChars[l] = 0;
        fwrite(&ElementChars, l, 1, fp);

        l = runEncodeB64(sizeof(int), (unsigned char *)(ien + 3 * Th.nt + i * 2 + 1), ElementChars);
        ElementChars[l] = 0;
        fwrite(&ElementChars, l, 1, fp);
      }
    }

    // flush buffer
    l = runEncodeB64(0, NULL, ElementChars);
    ElementChars[l] = 0;
    fwrite(&ElementChars, l, 1, fp);
  } else {
    fprintf(fp, " format=\"ascii\">\n");

    for (long i = 0; i < Th.nt; i++) {
      fprintf(fp, "%ld %ld %ld ", ien[i * 3 + 0], ien[i * 3 + 1], ien[i * 3 + 2]);
    }

    if (surface) {
      for (long i = 0; i < Th.neb; i++) {
        fprintf(fp, "%ld %ld ", ien[i * 2 + 3 * Th.nt], ien[i * 2 + 3 * Th.nt + 1]);
      }
    }
  }

  fprintf(fp, "\n</DataArray>\n");
  delete[] ien;

  fprintf(fp, "<DataArray type=\"Int32\" Name=\"offsets\" ");
  if (binary) {
    fprintf(fp, "format=\"binary\">\n    ");
    {
      unsigned nbytes = nc * sizeof(int);
      int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);
      long nelem = 3;

      for (long i = nelem; i <= nelem * Th.nt; i += nelem) {
        l = runEncodeB64(sizeof(int), (unsigned char *)&i, ElementChars);
        ElementChars[l] = 0;
        fwrite(&ElementChars, l, 1, fp);
      }

      if (surface) {
        nelem = 2;

        for (long i = nelem + 3 * Th.nt; i <= nelem * Th.neb + 3 * Th.nt; i += nelem) {
          l = runEncodeB64(sizeof(int), (unsigned char *)&i, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      }

      // flush buffer
      l = runEncodeB64(0, NULL, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);
    }
  } else {
    fprintf(fp, "format=\"ascii\" >\n");
    long nelem = 3;

    for (long i = nelem; i <= nelem * Th.nt; i += nelem) {
      fprintf(fp, "%ld ", i);
    }

    if (surface) {
      nelem = 2;

      for (long i = nelem; i <= nelem * Th.neb; i += nelem) {
        fprintf(fp, "%ld ", i + 3 * Th.nt);
      }
    }
  }

  fprintf(fp, "\n</DataArray>\n");

  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" ");
  if (binary) {
    fprintf(fp, "format=\"binary\" >\n    ");
    {
      unsigned nbytes = nc;
      int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      for (long i = 0; i < Th.nt; i++) {
        unsigned char types = 5;
        l = runEncodeB64(1, &types, ElementChars);
        ElementChars[l] = 0;
        fwrite(&ElementChars, l, 1, fp);
      }

      if (surface) {
        for (long i = 0; i < Th.neb; i++) {
          unsigned char types = 3;
          l = runEncodeB64(1, &types, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      }

      // flush buffer
      l = runEncodeB64(0, NULL, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);
    }
  } else {
    fprintf(fp, "format=\"ascii\" >\n");

    for (long i = 0; i < Th.nt; i++) {
      unsigned int types = 5;
      fprintf(fp, "%d ", types);
    }

    if (surface) {
      for (long i = 0; i < Th.neb; i++) {
        unsigned int types = 3;
        fprintf(fp, "%d ", types);
      }
    }
  }

  fprintf(fp, "\n</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  //---------------------------------- LABELS WITH VTU -------------------------------------------//
  fprintf(fp, "<CellData Scalars=\"Label\">\n");
  if (binary) {
    fprintf(fp, "<DataArray type=\"Int32\" Name=\"Label\" format=\"binary\">\n");
    int label;
    int nl = Th.nt;
    if (surface) nl += Th.neb;
    nl *= sizeof(int);
    //  cout << " nl =="<< nl << endl;
    writebin64(fp, nl);
    for (int it = 0; it < Th.nt; it++) {
      const Mesh::Triangle &K(Th.t(it));
      label = K.lab;
      writebin64(fp, label);
    }

    if (surface) {
      for (int ibe = 0; ibe < Th.neb; ibe++) {
        const Mesh::BorderElement &K(Th.be(ibe));
        label = K.lab;
        writebin64(fp, label);
      }
    }
    writebin64flush(fp);
  } else {
    fprintf(fp, "<DataArray type=\"Int32\" Name=\"Label\" format=\"ascii\">\n");
    int label;

    for (int it = 0; it < Th.nt; it++) {
      const Mesh::Triangle &K(Th.t(it));
      label = K.lab;
      fprintf(fp, "%d ", label);
    }

    if (surface) {
      for (int ibe = 0; ibe < Th.neb; ibe++) {
        const Mesh::BorderElement &K(Th.be(ibe));
        label = K.lab;
        fprintf(fp, "%d ", label);
      }
    }
  }
  fprintf(fp, "\n</DataArray>\n");
  //---------------------------------- LABELS WITH VTU -------------------------------------------//
}

void VTU_WRITE_MESH(FILE *fp, const Mesh3 &Th, bool binary, int datasize, bool surface) {
  int nc, nv, nconnex;

  if (surface) {
    nc = Th.nt + Th.nbe;
  } else {
    nc = Th.nt;
  }

  if (surface) {
    nconnex = 4 * Th.nt + 3 * Th.nbe;
  } else {
    nconnex = 4 * Th.nt;
  }

  unsigned char ElementChars[256];
  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" ");
  if (isBigEndian( )) {
    fprintf(fp, "byte_order=\"BigEndian\">\n");
  } else {
    fprintf(fp, " byte_order=\"LittleEndian\">\n");
  }

  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\" %d\">\n", Th.nv, nc);

  fprintf(fp, "<Points>\n");
  fprintf(fp, "<DataArray type=\"Float32\" NumberOfComponents=\"3\"");

  // A definir coord

  float *coord = new float[3 * Th.nv];

  for (int ii = 0; ii < Th.nv; ii++) {
    coord[3 * ii] = Th.vertices[ii].x;
    coord[3 * ii + 1] = Th.vertices[ii].y;
    coord[3 * ii + 2] = Th.vertices[ii].z;
  }

  if (binary) {
    fprintf(fp, " format=\"binary\">\n	");
    {
      unsigned nbytes = 3 * Th.nv * sizeof(float);
      int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      for (long i = 0; i < Th.nv; i++) {
        l = runEncodeB64(3 * sizeof(float), (unsigned char *)(coord + 3 * i), ElementChars);
        ElementChars[l] = 0;
        fwrite(&ElementChars, l, 1, fp);
      }

      // flush buffer
      l = runEncodeB64(0, NULL, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);
    }
  } else {
    fprintf(fp, " format=\"ascii\">\n");

    for (long i = 0; i < Th.nv; i++) {
      fprintf(fp, "%f %f %f ", coord[i * 3 + 0], coord[i * 3 + 1], coord[i * 3 + 2]);
    }
  }

  fprintf(fp, "\n</DataArray>\n");
  fprintf(fp, "</Points>\n");
  fprintf(fp, "<Cells>\n");

  delete[] coord;
  long *ien = new long[nconnex];

  for (long it = 0; it < Th.nt; it++) {
    const Tet &K(Th.elements[it]);

    for (int ii = 0; ii < 4; ii++) {
      ien[4 * it + ii] = Th.operator( )(K[ii]);
    }
  }

  if (surface) {
    for (long ibe = 0; ibe < Th.nbe; ibe++) {
      const Triangle3 &K(Th.be(ibe));

      for (int ii = 0; ii < 3; ii++) {
        ien[4 * Th.nt + 3 * ibe + ii] = Th.operator( )(K[ii]);
      }
    }
  }

  fprintf(fp, "<DataArray type=\"Int32\" Name=\"connectivity\" ");

  // need ien vector
  if (binary) {
    fprintf(fp, "format=\"binary\">\n	");
    unsigned nbytes = nconnex * sizeof(int);
    int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
    ElementChars[l] = 0;
    fwrite(&ElementChars, l, 1, fp);

    for (long i = 0; i < Th.nt; i++) {
      long nelem = 4;
      l = runEncodeB64(sizeof(int), (unsigned char *)(ien + i * nelem), ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      l = runEncodeB64(sizeof(int), (unsigned char *)(ien + i * nelem + 1), ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      l = runEncodeB64(sizeof(int), (unsigned char *)(ien + i * nelem + 2), ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      l = runEncodeB64(sizeof(int), (unsigned char *)(ien + i * nelem + 3), ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);
    }

    if (surface) {
      for (long i = 0; i < Th.nbe; i++) {
        long nelem = 3;
        l = runEncodeB64(sizeof(int), (unsigned char *)(ien + 4 * Th.nt + i * nelem), ElementChars);
        ElementChars[l] = 0;
        fwrite(&ElementChars, l, 1, fp);

        l = runEncodeB64(sizeof(int), (unsigned char *)(ien + 4 * Th.nt + i * nelem + 1),
                         ElementChars);
        ElementChars[l] = 0;
        fwrite(&ElementChars, l, 1, fp);

        l = runEncodeB64(sizeof(int), (unsigned char *)(ien + 4 * Th.nt + i * nelem + 2),
                         ElementChars);
        ElementChars[l] = 0;
        fwrite(&ElementChars, l, 1, fp);
      }
    }

    // flush buffer
    l = runEncodeB64(0, NULL, ElementChars);
    ElementChars[l] = 0;
    fwrite(&ElementChars, l, 1, fp);
  } else {
    fprintf(fp, " format=\"ascii\">\n");

    for (long i = 0; i < Th.nt; i++) {
      fprintf(fp, "%ld %ld %ld %ld ", ien[i * 4 + 0], ien[i * 4 + 1], ien[i * 4 + 2],
              ien[i * 4 + 3]);    // J.Morice 01/11
    }

    if (surface) {
      for (long i = 0; i < Th.nbe; i++) {
        fprintf(fp, "%ld %ld %ld ", ien[i * 3 + 4 * Th.nt], ien[i * 3 + 4 * Th.nt + 1],
                ien[i * 3 + 4 * Th.nt + 2]);
      }
    }
  }

  fprintf(fp, "\n</DataArray>\n");
  delete[] ien;

  fprintf(fp, "<DataArray type=\"Int32\" Name=\"offsets\" ");
  if (binary) {
    fprintf(fp, "format=\"binary\">\n	");
    {
      unsigned nbytes = nc * sizeof(int);
      int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);
      long nelem = 4;

      for (long i = nelem; i <= nelem * Th.nt; i += nelem) {
        l = runEncodeB64(sizeof(int), (unsigned char *)&i, ElementChars);
        ElementChars[l] = 0;
        fwrite(&ElementChars, l, 1, fp);
      }

      if (surface) {
        nelem = 3;

        for (long i = nelem + 4 * Th.nt; i <= nelem * Th.nbe + 4 * Th.nt; i += nelem) {
          l = runEncodeB64(sizeof(int), (unsigned char *)&i, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      }

      // flush buffer
      l = runEncodeB64(0, NULL, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);
    }
  } else {
    fprintf(fp, "format=\"ascii\" >\n");
    long nelem = 4;

    for (long i = nelem; i <= nelem * Th.nt; i += nelem) {
      fprintf(fp, "%ld ", i);
    }

    if (surface) {
      nelem = 3;

      for (long i = nelem + 4 * Th.nt; i <= nelem * Th.nbe + 4 * Th.nt; i += nelem) {
        fprintf(fp, "%ld ", i);
      }
    }
  }

  fprintf(fp, "\n</DataArray>\n");

  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" ");
  if (binary) {
    fprintf(fp, "format=\"binary\" >\n	");
    {
      unsigned nbytes = nc;
      int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      for (long i = 0; i < Th.nt; i++) {
        unsigned char types = 10;
        l = runEncodeB64(1, &types, ElementChars);
        ElementChars[l] = 0;
        fwrite(&ElementChars, l, 1, fp);
      }

      if (surface) {
        for (long i = 0; i < Th.nbe; i++) {
          unsigned char types = 5;
          l = runEncodeB64(1, &types, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      }

      // flush buffer
      l = runEncodeB64(0, NULL, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);
    }
  } else {
    fprintf(fp, "format=\"ascii\" >\n");

    for (long i = 0; i < Th.nt; i++) {
      unsigned int types = 10;
      fprintf(fp, "%d ", types);
    }

    if (surface) {
      for (long i = 0; i < Th.nbe; i++) {
        unsigned int types = 5;
        fprintf(fp, "%d ", types);
      }
    }
  }

  fprintf(fp, "\n</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  //---------------------------------- LABELS WITH VTU -------------------------------------------//
  fprintf(fp, "<CellData Scalars=\"Label\">\n");
  if (binary) {
    fprintf(fp, "<DataArray type=\"Int32\" Name=\"Label\" format=\"binary\">\n");
    int label;
    int nl = Th.nt;
    if (surface) nl += Th.nbe;
    nl *= sizeof(int);
    //  cout << " nl =="<< nl << endl;
    writebin64(fp, nl);
    for (int it = 0; it < Th.nt; it++) {
      label = Th[it].lab;
      writebin64(fp, label);
    }

    if (surface) {
      for (int ibe = 0; ibe < Th.nbe; ibe++) {
        label = Th.be(ibe).lab;
        writebin64(fp, label);
      }
    }
    writebin64flush(fp);
  } else {
    fprintf(fp, "<DataArray type=\"Int32\" Name=\"Label\" format=\"ascii\">\n");
    int label;

    for (int it = 0; it < Th.nt; it++) {
      const Tet &K(Th.elements[it]);
      label = K.lab;
      fprintf(fp, "%d\n", label);
    }

    if (surface) {
      for (int ibe = 0; ibe < Th.nbe; ibe++) {
        const Triangle3 &K(Th.be(ibe));
        label = K.lab;
        fprintf(fp, "%d\n", label);
      }
    }
  }
  fprintf(fp, "\n</DataArray>\n");
  //---------------------------------- LABELS WITH VTU -------------------------------------------//
}

template< class MMesh >
void VTU_WRITE_MESHT(FILE *fp, const MMesh &Th, bool binary, int datasize, bool surface) {
  typedef typename MMesh::Element T;
  typedef typename MMesh::BorderElement B;

  int nv, nconnex;
  int nc = surface ? Th.nt + Th.nbe : Th.nt;
  if (std::is_same< MMesh, MeshS >::value)
    nconnex = surface ? nconnex = 3 * Th.nt + 2 * Th.nbe : 3 * Th.nt;
  else if (std::is_same< MMesh, MeshL >::value)
    nconnex = surface ? nconnex = 2 * Th.nt + Th.nbe : 2 * Th.nt;

  unsigned char ElementChars[256];
  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" ");
  if (isBigEndian( ))
    fprintf(fp, "byte_order=\"BigEndian\">\n");
  else
    fprintf(fp, " byte_order=\"LittleEndian\">\n");

  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\" %d\">\n", Th.nv, nc);

  fprintf(fp, "<Points>\n");
  fprintf(fp, "<DataArray type=\"Float32\" NumberOfComponents=\"3\"");

  // A definir coord

  float *coord = new float[3 * Th.nv];

  for (int ii = 0; ii < Th.nv; ii++) {
    coord[3 * ii] = Th.vertices[ii].x;
    coord[3 * ii + 1] = Th.vertices[ii].y;
    coord[3 * ii + 2] = Th.vertices[ii].z;
  }

  if (binary) {
    fprintf(fp, " format=\"binary\">\n    ");
    {
      unsigned nbytes = 3 * Th.nv * sizeof(float);
      int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      for (long i = 0; i < Th.nv; i++) {
        l = runEncodeB64(3 * sizeof(float), (unsigned char *)(coord + 3 * i), ElementChars);
        ElementChars[l] = 0;
        fwrite(&ElementChars, l, 1, fp);
      }

      // flush buffer
      l = runEncodeB64(0, NULL, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);
    }
  } else {
    fprintf(fp, " format=\"ascii\">\n");

    for (long i = 0; i < Th.nv; i++) {
      fprintf(fp, "%f %f %f ", coord[i * 3 + 0], coord[i * 3 + 1], coord[i * 3 + 2]);
    }
  }

  fprintf(fp, "\n</DataArray>\n");
  fprintf(fp, "</Points>\n");
  fprintf(fp, "<Cells>\n");

  delete[] coord;
  long *ien = new long[nconnex];

  for (long it = 0; it < Th.nt; it++) {
    const T &K(Th.t(it));
    for (int ii = 0; ii < T::nv; ii++) ien[(T::nv)*it + ii] = Th.operator( )(K[ii]);
  }

  if (surface) {
    for (long ibe = 0; ibe < Th.nbe; ibe++) {
      const B &K(Th.be(ibe));
      for (int ii = 0; ii < B::nv; ii++)
        ien[(T::nv)*Th.nt + (B::nv)*ibe + ii] = Th.operator( )(K[ii]);
    }
  }

  fprintf(fp, "<DataArray type=\"Int32\" Name=\"connectivity\" ");

  // need ien vector
  if (binary) {
    fprintf(fp, "format=\"binary\">\n    ");
    unsigned nbytes = nconnex * sizeof(int);
    int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
    ElementChars[l] = 0;
    fwrite(&ElementChars, l, 1, fp);

    for (long i = 0; i < Th.nt; i++) {
      l = runEncodeB64(sizeof(int), (unsigned char *)(ien + i * (T::nv)), ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      l = runEncodeB64(sizeof(int), (unsigned char *)(ien + i * (T::nv) + 1), ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      if (std::is_same< MMesh, MeshS >::value) {
        l = runEncodeB64(sizeof(int), (unsigned char *)(ien + i * (T::nv) + 2), ElementChars);
        ElementChars[l] = 0;
        fwrite(&ElementChars, l, 1, fp);
      }
    }

    if (surface) {
      for (long i = 0; i < Th.nbe; i++) {
        l = runEncodeB64(sizeof(int), (unsigned char *)(ien + (T::nv)*Th.nt + i * (B::nv)),
                         ElementChars);
        ElementChars[l] = 0;
        fwrite(&ElementChars, l, 1, fp);
        if (std::is_same< MMesh, MeshS >::value) {
          l = runEncodeB64(sizeof(int), (unsigned char *)(ien + (T::nv)*Th.nt + i * (B::nv) + 1),
                           ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      }
    }

    // flush buffer
    l = runEncodeB64(0, NULL, ElementChars);
    ElementChars[l] = 0;
    fwrite(&ElementChars, l, 1, fp);
  } else {
    fprintf(fp, " format=\"ascii\">\n");

    for (long i = 0; i < Th.nt; i++)
      if (std::is_same< MMesh, MeshS >::value)
        fprintf(fp, "%ld %ld %ld ", ien[i * 3 + 0], ien[i * 3 + 1], ien[i * 3 + 2]);
      else if (std::is_same< MMesh, MeshL >::value)
        fprintf(fp, "%ld %ld ", ien[i * 2 + 0], ien[i * 2 + 1]);

    if (surface) {
      for (long i = 0; i < Th.nbe; i++)
        if (std::is_same< MMesh, MeshS >::value)
          fprintf(fp, "%ld %ld ", ien[i * 2 + 3 * Th.nt], ien[i * 2 + 3 * Th.nt + 1]);
        else if (std::is_same< MMesh, MeshL >::value)
          fprintf(fp, "%ld ", ien[i + 2 * Th.nt]);
    }
  }

  fprintf(fp, "\n</DataArray>\n");
  delete[] ien;

  fprintf(fp, "<DataArray type=\"Int32\" Name=\"offsets\" ");
  // long nelem;
  if (binary) {
    fprintf(fp, "format=\"binary\">\n    ");
    {
      unsigned nbytes = nc * sizeof(int);
      int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      for (long i = (T::nv); i <= (T::nv)*Th.nt; i += (T::nv)) {
        l = runEncodeB64(sizeof(int), (unsigned char *)&i, ElementChars);
        ElementChars[l] = 0;
        fwrite(&ElementChars, l, 1, fp);
      }

      if (surface) {
        for (long i = (B::nv) + (T::nv)*Th.nt; i <= (B::nv)*Th.nbe + (T::nv)*Th.nt; i += (B::nv)) {
          l = runEncodeB64(sizeof(int), (unsigned char *)&i, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      }

      // flush buffer
      l = runEncodeB64(0, NULL, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);
    }
  } else {
    fprintf(fp, "format=\"ascii\" >\n");
    for (long i = (T::nv); i <= (T::nv)*Th.nt; i += (T::nv)) fprintf(fp, "%ld ", i);

    if (surface) {
      for (long i = (B::nv); i <= (B::nv)*Th.nbe; i += (B::nv))
        fprintf(fp, "%ld ", i + (T::nv)*Th.nt);
    }
  }

  fprintf(fp, "\n</DataArray>\n");

  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" ");
  unsigned char types;
  if (binary) {
    fprintf(fp, "format=\"binary\" >\n    ");
    {
      unsigned nbytes = nc;
      int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      for (long i = 0; i < Th.nt; i++) {
        if (std::is_same< MMesh, MeshS >::value)
          types = 5;
        else if (std::is_same< MMesh, MeshL >::value)
          types = 3;
        l = runEncodeB64(1, &types, ElementChars);
        ElementChars[l] = 0;
        fwrite(&ElementChars, l, 1, fp);
      }

      if (surface) {
        for (long i = 0; i < Th.nbe; i++) {
          if (std::is_same< MMesh, MeshS >::value)
            types = 3;
          else if (std::is_same< MMesh, MeshL >::value)
            types = 1;
          l = runEncodeB64(1, &types, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      }

      // flush buffer
      l = runEncodeB64(0, NULL, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);
    }
  } else {
    fprintf(fp, "format=\"ascii\" >\n");

    for (long i = 0; i < Th.nt; i++) {
      if (std::is_same< MMesh, MeshS >::value)
        types = 5;
      else if (std::is_same< MMesh, MeshL >::value)
        types = 3;
      fprintf(fp, "%d ", types);
    }

    if (surface) {
      for (long i = 0; i < Th.nbe; i++) {
        if (std::is_same< MMesh, MeshS >::value)
          types = 3;
        else if (std::is_same< MMesh, MeshL >::value)
          types = 1;
        fprintf(fp, "%d ", types);
      }
    }
  }

  fprintf(fp, "\n</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  //---------------------------------- LABELS WITH VTU -------------------------------------------//
  fprintf(fp, "<CellData Scalars=\"Label\">\n");
  int label;
  if (binary) {
    fprintf(fp, "<DataArray type=\"Int32\" Name=\"Label\" format=\"binary\">\n");

    int nl = Th.nt;
    if (surface) nl += Th.nbe;
    nl *= sizeof(int);
    writebin64(fp, nl);
    for (int it = 0; it < Th.nt; it++) {
      label = Th[it].lab;
      writebin64(fp, label);
    }
    if (surface) {
      for (int ibe = 0; ibe < Th.nbe; ibe++) {
        label = Th.be(ibe).lab;
        writebin64(fp, label);
      }
    }
    writebin64flush(fp);
  } else {
    fprintf(fp, "<DataArray type=\"Int32\" Name=\"Label\" format=\"ascii\">\n");
    for (int it = 0; it < Th.nt; it++) {
      const T &K(Th.elements[it]);
      label = K.lab;
      fprintf(fp, "%d\n", label);
    }
    if (surface) {
      for (int ibe = 0; ibe < Th.nbe; ibe++) {
        const B &K(Th.be(ibe));
        label = K.lab;
        fprintf(fp, "%d\n", label);
      }
    }
  }
  fprintf(fp, "\n</DataArray>\n");
  fprintf(fp, "</CellData>\n");
  //---------------------------------- LABELS WITH VTU -------------------------------------------//
}

//= =============================================
// LOAD DE FICHIER .vtk (2D)
//= =============================================

class VTK_LoadMesh_Op : public E_F0mps {
 public:
  Expression filename;
  static const int n_name_param = 5;    //
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
   
  KN<KN<double> >* arg(int i, Stack stack,KN<KN<double> >*p) const {
      return nargs[i] ? GetAny< KN<KN<double> >* >((*nargs[i])(stack)) : p;
  }
    
 public:
  VTK_LoadMesh_Op(const basicAC_F0 &args, Expression ffname) : filename(ffname) {
    if (verbosity) {
      cout << "Load mesh given by VTK " << endl;
    }

    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type VTK_LoadMesh_Op::name_param[] = {{"reft", &typeid(long)},
                                                           {"swap", &typeid(bool)},
                                                           {"refe", &typeid(long)},
                                                           {"namelabel", &typeid(string)},
                                                           {"fields", &typeid( KN<KN<double> >*)}};

class VTK_LoadMesh : public OneOperator {
 public:
  VTK_LoadMesh( ) : OneOperator(atype< pmesh >( ), atype< string * >( )) {}

  E_F0 *code(const basicAC_F0 &args) const {
    return new VTK_LoadMesh_Op(args, t[0]->CastTo(args[0]));
  }
};

Mesh *VTK_Load(const string &filename, bool bigEndian, KN<KN<double> >* pfields) {
  // swap = bigEndian or not bigEndian
  // variable freefem++
  int nv, nt = 0, nbe = 0;
  int nerr = 0;
  char *res;
  map< int, int > mapnumv;

  // Reading Mesh in vtk formats

  FILE *fp = fopen(filename.c_str( ), "rb");
  if (!fp) {
    cerr << "Unable to open file " << filename.c_str( ) << endl;
    exit(1);
  }

  char buffer[256], buffer2[256];

  res = fgets(buffer, sizeof(buffer), fp);    // version line
  res = fgets(buffer, sizeof(buffer), fp);    // title

  fscanf(fp, "%s", buffer);    // ASCII or BINARY
  bool binary = false;
  if (!strncmp(buffer, "BINARY", 6)) {
    binary = true;
  }

  if (fscanf(fp, "%s %s", buffer, buffer2) != 2) {
    cerr << "error in reading vtk files" << filename << endl;
    ExecError("load vtk files");
  }

  if (strcmp(buffer, "DATASET") || strcmp(buffer2, "UNSTRUCTURED_GRID")) {
    cout << "VTK reader can only read unstructured datasets" << endl;
    ExecError("load  vtk files");
    exit(1);
  }

  // read mesh vertices
  if (fscanf(fp, "%s %d %s\n", buffer, &nv, buffer2) != 3) {
    cout << "error in reading vtk files" << endl;
    ExecError("load  vtk files reah vertices");
    exit(1);
  }

  if (strcmp(buffer, "POINTS") || !nv) {
    cerr << "No points in dataset" << endl;
    ExecError("load  vtk files: No points in dataset");
    exit(1);
  }

  int datasize;
  if (!strncmp(buffer2, "double", 6)) {
    datasize = sizeof(double);
  } else if (!strncmp(buffer2, "float", 5)) {
    datasize = sizeof(float);
  } else {
    cout << "VTK reader only accepts float or double datasets" << endl;
    ExecError("load  vtk files VTK reader only accepts float or double datasets");
    exit(1);
  }

  if (verbosity > 1) {
    cout << "   vtkio: Reading points " << nv << endl;
  }

  Mesh::Vertex *vff = new Mesh::Vertex[nv];

  for (int i = 0; i < nv; i++) {
    double xyz[3];
    if (binary) {
      if (datasize == sizeof(float)) {
        float f[3];
        if (fread(f, sizeof(float), 3, fp) != 3) {
          ExecError("load  vtk files VTK reader only accepts float or double datasets");
          ExecError("error in reading vtk file");
        }

        if (!bigEndian) {
          FreeFEM::SwapBytes((char *)f, sizeof(float), 3);
        }

        for (int j = 0; j < 3; j++) {
          xyz[j] = f[j];
        }
      } else {
        if (fread(xyz, sizeof(double), 3, fp) != 3) {
          cout << "error in reading vtk files" << endl;
          ExecError("error in reading vtk file");
        }

        if (!bigEndian) {
          FreeFEM::SwapBytes((char *)xyz, sizeof(double), 3);
        }
      }
    } else {
      if (fscanf(fp, "%lf %lf %lf", &xyz[0], &xyz[1], &xyz[2]) != 3) {
        cout << "error in reading vtk files" << endl;
        ExecError("error in reading vtk file");
      }
    }

    vff[i].x = xyz[0];
    vff[i].y = xyz[1];
    if (abs(xyz[2]) > 1.e-7) {
      cout << "we are plotted a two dimensional mesh: z coordinate must be 0" << endl;
      ExecError(
        "error in reading vtk file,we are plotted a two dimensional mesh: z coordinate must be 0");
    }

    vff[i].lab = 1;
  }

  // read mesh elements
  int numElements, numElements2, totalNumInt;
  if (fscanf(fp, "%s %d %d\n", buffer, &numElements, &totalNumInt) != 3) {
    cout << "error in reading vtk files" << endl;
    ExecError("error in reading vtk file");
  }

  if (strcmp(buffer, "CELLS") || !numElements) {
    cout << "No cells in dataset" << endl;
    ExecError("error in reading vtk file");
  }

  if(verbosity>2) cout << "Reading cells" << numElements << endl;

  int *IntCells = new int[totalNumInt - numElements];
  int *firstCell = new int[numElements + 1];
  int *TypeCells = new int[numElements];
  int numIntCells = 0;

  for (unsigned int i = 0; i < numElements; i++) {
    int numVerts, n[100];

    for (int ii = 0; ii < 100; ii++) {
      n[ii] = -1;
    }

    if (binary) {
      if (fread(&numVerts, sizeof(int), 1, fp) != 1) {
        cout << "error in reading VTK files " << endl;
        ExecError("error in reading vtk file");
      }

      if (!bigEndian) {
        FreeFEM::SwapBytes((char *)&numVerts, sizeof(int), 1);
      }

      if ((int)fread(n, sizeof(int), numVerts, fp) != numVerts) {
        cout << "error in reading VTK files " << endl;
        ExecError("error in reading vtk file");
      }

      if (!bigEndian) {
        FreeFEM::SwapBytes((char *)n, sizeof(int), numVerts);
      }
    } else {
      if (fscanf(fp, "%d", &numVerts) != 1) {
        cout << "error in reading VTK files " << endl;
        ExecError("error in reading vtk file");
      }

      for (int j = 0; j < numVerts; j++) {
        if (fscanf(fp, "%d", &n[j]) != 1) {
          cout << "error in reading VTK files " << endl;
          ExecError("error in reading vtk file");
        }
      }
    }

    firstCell[i] = numIntCells;

    for (int j = 0; j < numVerts; j++) {
      if (n[j] >= 0 && n[j] < nv) {
        IntCells[numIntCells] = n[j];
        numIntCells++;
      } else {
        cout << "Bad vertex index" << endl;
        ExecError("error in reading vtk file");
      }
    }
  }

  firstCell[numElements] = totalNumInt - numElements;

  if (fscanf(fp, "%s %d\n", buffer, &numElements2) != 2) {
    cout << " Error in reading CELL_TYPES ARGUMENT " << endl;
    ExecError("error in reading vtk file");
  }

  if (strcmp(buffer, "CELL_TYPES") || numElements2 != (int)numElements) {
    cout << "No or invalid number of cells types" << endl;
    ExecError("error in reading vtk file");
  }

  // 2D

  for (unsigned int i = 0; i < numElements; i++) {
    int type;
    if (binary) {
      if (fread(&type, sizeof(int), 1, fp) != 1) {
        cout << "bug in readings cell types" << endl;
        ExecError("error in reading vtk file");
      }

      if (!bigEndian) {
        FreeFEM::SwapBytes((char *)&type, sizeof(int), 1);
      }
    } else {
      if (fscanf(fp, "%d", &type) != 1) {
        cout << "bug in readings cell types" << endl;
        ExecError("error in reading vtk file");
      }
    }

    TypeCells[i] = type;

    switch (type) {
      case 1:    // Vertex
        if (nerr++ < 10 && verbosity) {
          cout << "this type of cell (vertex) is not taking account in Freefem++ " << type << endl;
        }

        break;
      case 3:     // Edge/line
        nbe++;    // 2D
        break;
      case 5:    // Triangle
        nt++;    // 2D
        break;
      case 10:    // Tetrah�dre
        cout << "We are loading a three dimensional mesh. Three is no tetrahedron." << endl;
        ExecError("error in reading vtk file");
        break;
      default:
        cout << "Error :: This type of cell is not considered in Freefem++ " << type << endl;
        ExecError("error in reading vtk file");
        break;
    }
  }

  if(pfields) {
    if(verbosity>1)   cout << " try  reading POINT_DATA  " << endl;
         
    int nbp=0,nbf=0,err=0;
    if (fscanf(fp, "%s %d", buffer, &nbp) != 2) {
      if (verbosity) cout << "error in reading vtk files pfields" << endl;
      err++;
    }

    if (strcmp(buffer, "POINT_DATA"))  {
      if (verbosity) cout << "VTK reader can only read POINT_DATA datasets:  not " << buffer<<  endl;
      err++;
    }
    if ((!err) &&(fscanf(fp, "%s %s %d", buffer, buffer2,&nbf) != 3)) {
      if (verbosity) cout << "error in reading vtk files FIELD FieldData" << endl;
      err++;
    }

    if(err) nbf=0;
      else pfields->resize(nbf);
      for(int nf=0 ; nf < nbf; nf++) {

      int m,nv;
      // read mesh vertices
      if (fscanf(fp, "%s %d %d %s\n", buffer,&m, &nv, buffer2) != 4) {
        if (verbosity) cout << "error in reading vtk files " << endl;
        err++;
        break;
      }
      int n = m*nv;
      if(verbosity) cout << " reading "<< buffer << " "<< m << " " << nv << " "<< buffer2 << endl;
      int datasize;
      if (!strncmp(buffer2, "double", 6))
        datasize = sizeof(double);
      else if (!strncmp(buffer2, "float", 5))
        datasize = sizeof(float);
      else {
          if (verbosity) cout << "VTK reader only accepts float or double datasets" << endl;
          err++;
          break;
      }
      //  read data ..
      if(err) break;
      (*pfields)[nf].resize(n);
      double* pv=&(*pfields)[nf][0];
      for(int i=0; i<n; ++i) {
        if (binary) {
          if (datasize == sizeof(float)) {
            float f[1];
            if (fread(f, sizeof(float), 1, fp) != 1) {
              if (verbosity) cout << "error in reading  vtk fields float at " << nf << " / " << i << endl;
              err++;
              break;
            }
            if (!bigEndian)
              FreeFEM::SwapBytes((char *)f, sizeof(float), 1);
              *pv++=*f;
          }
          else {
            if (fread(pv++, sizeof(double), 1, fp) != 1) {
              if (verbosity) cout << "error in reading  vtk fields double at " << nf << " / " << i << endl;
              err++;
              break;
            }
            if (!bigEndian)
              FreeFEM::SwapBytes((char *)pv, sizeof(double), 1);
          }
        }
        else {
          if (fscanf(fp, "%lf", pv++) != 1) {
            if (verbosity) cout << "error in reading vtk files ascii fields" << nf << " / " << i <<endl;
            err++;
          }
        }
        if(err) break;
      }
       if(err) break;
    }
    if( err ) (*pfields).resize(0);
  } // end if(pfields)
    
  fclose(fp);

  // 2D Versions
  Mesh::Triangle *tff;
  if (nt > 0) {
    tff = new Mesh::Triangle[nt];
  }

  Mesh::Triangle *ttff = tff;
  Mesh::BorderElement *bff;
  if (nbe > 0) {
    bff = new Mesh::BorderElement[nbe];
  }

  Mesh::BorderElement *bbff = bff;

  for (unsigned int i = 0; i < numElements; i++) {
    int type = TypeCells[i];
    int iv[3];
    int label = 1;

    switch (type) {
      case 1:    // Vertex

        if (nerr++ < 10 && verbosity) {
          cout << "this type of cell (vertex) is not taking account in Freefem++ " << type << " "
               << endl;
        }

        break;
      case 3:    // Edge/line
        assert((firstCell[i + 1] - firstCell[i]) == 2);

        for (int j = firstCell[i]; j < firstCell[i + 1]; j++) {
          iv[j - firstCell[i]] = IntCells[j];
        }

        (bbff++)->set(vff, iv[0], iv[1], label);
        break;
      case 5:    // Triangle
        assert((firstCell[i + 1] - firstCell[i]) == 3);

        for (int j = firstCell[i]; j < firstCell[i + 1]; j++) {
          iv[j - firstCell[i]] = IntCells[j];
        }

        (ttff++)->set(vff, iv[0], iv[1], iv[2], label);
        break;
      default:
        break;
    }
  }

  delete[] IntCells;
  delete[] firstCell;
  delete[] TypeCells;

  Mesh *pTh = new Mesh(nv, nt, nbe, vff, tff, bff);
  return pTh;
}

AnyType VTK_LoadMesh_Op::operator( )(Stack stack) const {
  string *pffname = GetAny< string * >((*filename)(stack));
  bool swap = false;
  int reftri = 1;
  int refedges = 1;

  if (nargs[0]) {
    reftri = GetAny< long >((*nargs[0])(stack));
  }

  if (nargs[1]) {
    swap = GetAny< bool >((*nargs[1])(stack));
  }

  if (nargs[2]) {
    refedges = GetAny< long >((*nargs[2])(stack));
  }

  string *DataLabel;
  if (nargs[3]) {
    DataLabel = GetAny< string * >((*nargs[3])(stack));
  }
  KN<KN<double> > * pfields=0;
  pfields=arg(4, stack,pfields);
  Mesh *Th = VTK_Load(*pffname, swap, pfields);

  // A faire fonction pour changer le label

  Add2StackOfPtr2FreeRC(stack, Th);

  return Th;
}

//= =============================================
// ECRITURE DE FICHIER .vtk (2D)
//= =============================================

class VTK_WriteMesh_Op : public E_F0mps {
 public:
  typedef long Result;
  Expression eTh;
  Expression filename;
  struct Expression2 {
    string name;
    long what;       // 1 scalar, 2 vector, 3 symtensor
    long nbfloat;    // 1 scalar, 2 vector (3D), 3 symtensor(3D)
    Expression e[3];
    Expression2( ) {
      e[0] = 0;
      e[1] = 0;
      e[2] = 0;
      what = 0;
      nbfloat = 0;
    };
    Expression &operator[](int i) { return e[i]; }

    double eval(int i, Stack stack) const {
      if (e[i]) {
        return GetAny< double >((*e[i])(stack));
      } else {
        return 0;
      }
    }

    void writesolutionP0_float_binary(FILE *fp, const Mesh &Th, Stack stack, bool surface,
                                      bool bigEndian) const {
      MeshPoint *mp3(MeshPointStack(stack));
      R2 Cdg_hat = R2(1. / 3., 1. / 3.);

      for (int it = 0; it < Th.nt; it++) {
        const Mesh::Triangle &K(Th.t(it));
        mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

        for (int j = 0; j < nbfloat; j++) {
          float value = eval(j, stack);

          if (!bigEndian) {
            FreeFEM::SwapBytes((char *)&value, sizeof(float), 1);
          }

          fwrite(&value, sizeof(float), 1, fp);
        }
      }

      if (surface) {
        for (int ibe = 0; ibe < Th.neb; ibe++) {
          // determination du triangle contenant cette edge
          int ie;
          int it = Th.BoundaryElement(ibe, ie);
          const Mesh::Triangle &K(Th.t(it));
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < nbfloat; j++) {
            float value = eval(j, stack);

            if (!bigEndian) {
              FreeFEM::SwapBytes((char *)&value, sizeof(float), 1);
            }

            fwrite(&value, sizeof(float), 1, fp);
          }
        }
      }

      fprintf(fp, "\n");
    }

    void writesolutionP0_float(FILE *fp, const Mesh &Th, Stack stack, bool surface) const {
      MeshPoint *mp3(MeshPointStack(stack));
      R2 Cdg_hat = R2(1. / 3., 1. / 3.);

      for (int it = 0; it < Th.nt; it++) {
        const Mesh::Triangle &K(Th.t(it));
        mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

        for (int j = 0; j < nbfloat; j++) {
          float value = eval(j, stack);

          fprintf(fp, "%.8e ", value);
        }
      }

      if (surface) {
        for (int ibe = 0; ibe < Th.neb; ibe++) {
          // determination du triangle contenant cette edge
          int ie;
          int it = Th.BoundaryElement(ibe, ie);
          const Mesh::Triangle &K(Th.t(it));
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < nbfloat; j++) {
            float value = eval(j, stack);

            fprintf(fp, "%.8e ", value);
          }
        }
      }

      fprintf(fp, "\n");
    }

    void writesolutionP0_float_XML(FILE *fp, const Mesh &Th, Stack stack, bool surface) const {
      MeshPoint *mp3(MeshPointStack(stack));
      R2 Cdg_hat = R2(1. / 3., 1. / 3.);
      unsigned char ElementChars[256];
      long nc = Th.nt;

      if (surface) {
        nc = nc + Th.neb;
      }

      unsigned nbytes = nc * sizeof(float) * (*this).nbfloat;
      int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      for (int it = 0; it < Th.nt; it++) {
        const Mesh::Triangle &K(Th.t(it));
        mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

        for (int j = 0; j < nbfloat; j++) {
          float value = eval(j, stack);

          l = runEncodeB64(sizeof(float), (unsigned char *)&value, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      }

      if (surface) {
        for (int ibe = 0; ibe < Th.neb; ibe++) {
          // determination du triangle contenant cette edge
          int ie;
          int it = Th.BoundaryElement(ibe, ie);
          const Mesh::Triangle &K(Th.t(it));
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < nbfloat; j++) {
            float value = eval(j, stack);

            l = runEncodeB64(sizeof(float), (unsigned char *)&value, ElementChars);
            ElementChars[l] = 0;
            fwrite(&ElementChars, l, 1, fp);
          }
        }
      }

      // flush buffer
      l = runEncodeB64(0, NULL, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);
      fprintf(fp, "\n");
    }

    void writesolutionP0_float(FILE *fp, const Mesh &Th, Stack stack, bool surface, bool binary,
                               bool bigEndian, bool XML = false) const {
      if (binary) {
        if (!XML) {
          (*this).writesolutionP0_float_binary(fp, Th, stack, surface, bigEndian);
        } else {
          (*this).writesolutionP0_float_XML(fp, Th, stack, surface);
        }
      } else {
        (*this).writesolutionP0_float(fp, Th, stack, surface);
      }
    }

    void writesolutionP0_double_binary(FILE *fp, const Mesh &Th, Stack stack, bool surface,
                                       bool bigEndian) const {
      MeshPoint *mp3(MeshPointStack(stack));
      R2 Cdg_hat = R2(1. / 3., 1. / 3.);

      for (int it = 0; it < Th.nt; it++) {
        const Mesh::Triangle &K(Th.t(it));
        mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

        for (int j = 0; j < (*this).nbfloat; j++) {
          double value = (*this).eval(j, stack);

          if (!bigEndian) {
            FreeFEM::SwapBytes((char *)&value, sizeof(double), 1);
          }

          fwrite(&value, sizeof(double), 1, fp);
        }
      }

      if (surface) {
        for (int ibe = 0; ibe < Th.neb; ibe++) {
          // determination du triangle contenant cette edge
          int ie;
          int it = Th.BoundaryElement(ibe, ie);
          const Mesh::Triangle &K(Th.t(it));
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < (*this).nbfloat; j++) {
            double value = (*this).eval(j, stack);

            if (!bigEndian) {
              FreeFEM::SwapBytes((char *)&value, sizeof(double), 1);
            }

            fwrite(&value, sizeof(double), 1, fp);
          }
        }
      }

      fprintf(fp, "\n");
    }

    void writesolutionP0_double(FILE *fp, const Mesh &Th, Stack stack, bool surface) const {
      MeshPoint *mp3(MeshPointStack(stack));
      R2 Cdg_hat = R2(1. / 3., 1. / 3.);
      unsigned char ElementChars[256];

      for (int it = 0; it < Th.nt; it++) {
        const Mesh::Triangle &K(Th.t(it));
        mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

        for (int j = 0; j < (*this).nbfloat; j++) {
          double value = (*this).eval(j, stack);

          fprintf(fp, "%.16e ", value);
        }
      }

      if (surface) {
        for (int ibe = 0; ibe < Th.neb; ibe++) {
          // determination du triangle contenant cette edge
          int ie;
          int it = Th.BoundaryElement(ibe, ie);
          const Mesh::Triangle &K(Th.t(it));
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < (*this).nbfloat; j++) {
            double value = (*this).eval(j, stack);
            fprintf(fp, "%.16e ", value);
          }
        }
      }

      fprintf(fp, "\n");
    }

    void writesolutionP0_double_XML(FILE *fp, const Mesh &Th, Stack stack, bool surface) const {
      unsigned char ElementChars[256];
      MeshPoint *mp3(MeshPointStack(stack));
      R2 Cdg_hat = R2(1. / 3., 1. / 3.);
      long nc = Th.nt;

      if (surface) {
        nc = nc + Th.neb;
      }

      unsigned nbytes = nc * sizeof(double) * (*this).nbfloat;
      int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      for (int it = 0; it < Th.nt; it++) {
        const Mesh::Triangle &K(Th.t(it));
        mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

        for (int j = 0; j < (*this).nbfloat; j++) {
          double value = (*this).eval(j, stack);

          l = runEncodeB64(sizeof(double), (unsigned char *)&value, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      }

      if (surface) {
        for (int ibe = 0; ibe < Th.neb; ibe++) {
          // determination du triangle contenant cette edge
          int ie;
          int it = Th.BoundaryElement(ibe, ie);
          const Mesh::Triangle &K(Th.t(it));
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < (*this).nbfloat; j++) {
            double value = (*this).eval(j, stack);

            l = runEncodeB64(sizeof(double), (unsigned char *)&value, ElementChars);
            ElementChars[l] = 0;
            fwrite(&ElementChars, l, 1, fp);
          }
        }
      }

      // flush buffer
      l = runEncodeB64(0, NULL, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      fprintf(fp, "\n");
    }

    void writesolutionP0_double(FILE *fp, const Mesh &Th, Stack stack, bool surface, bool binary,
                                bool bigEndian, bool XML = false) const {
      if (binary) {
        if (!XML) {
          (*this).writesolutionP0_double_binary(fp, Th, stack, surface, bigEndian);
        } else {
          (*this).writesolutionP0_double_XML(fp, Th, stack, surface);
        }
      } else {
        (*this).writesolutionP0_double(fp, Th, stack, surface);
      }
    }

    void writesolutionP1_float(FILE *fp, const Mesh &Th, Stack stack, bool binary, bool bigEndian,
                               bool XML = false) const {
      unsigned char ElementChars[256];
      MeshPoint *mp3(MeshPointStack(stack));

      KN< double > valsol(Th.nv * (*this).nbfloat);
      KN< int > takemesh(Th.nv);
      takemesh = 0;
      valsol = 0.;

      for (int it = 0; it < Th.nt; it++) {
        for (int iv = 0; iv < 3; iv++) {
          int i = Th(it, iv);
          mp3->setP(&Th, it, iv);

          for (int j = 0; j < (*this).nbfloat; j++) {
            valsol[i * (*this).nbfloat + j] =
              valsol[i * (*this).nbfloat + j] + (*this).eval(j, stack);
          }

          takemesh[i] = takemesh[i] + 1;
        }
      }

      if (binary) {
        if (!XML) {
          if (!bigEndian) {
            for (int iv = 0; iv < Th.nv; iv++) {
              for (int j = 0; j < (*this).nbfloat; j++) {
                valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
                float value = valsol[iv * (*this).nbfloat + j];

                FreeFEM::SwapBytes((char *)&value, sizeof(float), 1);
                fwrite(&value, sizeof(float), 1, fp);
              }
            }
          } else {
            for (int iv = 0; iv < Th.nv; iv++) {
              for (int j = 0; j < (*this).nbfloat; j++) {
                valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
                float value = valsol[iv * (*this).nbfloat + j];

                fwrite(&value, sizeof(float), 1, fp);
              }
            }
          }
        } else {
          long nc = Th.nv;
          unsigned nbytes = nc * sizeof(float) * (*this).nbfloat;
          int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);

          for (int iv = 0; iv < Th.nv; iv++) {
            for (int j = 0; j < (*this).nbfloat; j++) {
              valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
              float value = valsol[iv * (*this).nbfloat + j];

              l = runEncodeB64(sizeof(float), (unsigned char *)&value, ElementChars);
              ElementChars[l] = 0;
              fwrite(&ElementChars, l, 1, fp);
            }
          }

          // flush buffer
          l = runEncodeB64(0, NULL, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      } else {
        for (int iv = 0; iv < Th.nv; iv++) {
          for (int j = 0; j < (*this).nbfloat; j++) {
            valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
            float value = valsol[iv * (*this).nbfloat + j];

            fprintf(fp, "%.8e ", value);
          }
        }
      }

      fprintf(fp, "\n");
    }

    void writesolutionP1_double(FILE *fp, const Mesh &Th, Stack stack, bool binary, bool bigEndian,
                                bool XML = false) const {
      unsigned char ElementChars[256];
      MeshPoint *mp3(MeshPointStack(stack));

      KN< double > valsol(Th.nv * (*this).nbfloat);
      KN< int > takemesh(Th.nv);
      takemesh = 0;
      valsol = 0.;

      for (int it = 0; it < Th.nt; it++) {
        for (int iv = 0; iv < 3; iv++) {
          int i = Th(it, iv);
          mp3->setP(&Th, it, iv);

          for (int j = 0; j < (*this).nbfloat; j++) {
            valsol[i * (*this).nbfloat + j] =
              valsol[i * (*this).nbfloat + j] + (*this).eval(j, stack);
          }

          takemesh[i] = takemesh[i] + 1;
        }
      }

      if (binary) {
        if (!XML) {
          if (!bigEndian) {
            for (int iv = 0; iv < Th.nv; iv++) {
              for (int j = 0; j < (*this).nbfloat; j++) {
                valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
                double value = valsol[iv * (*this).nbfloat + j];

                FreeFEM::SwapBytes((char *)&value, sizeof(double), 1);
                fwrite(&value, sizeof(double), 1, fp);
              }
            }
          } else {
            for (int iv = 0; iv < Th.nv; iv++) {
              for (int j = 0; j < (*this).nbfloat; j++) {
                valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
                double value = valsol[iv * (*this).nbfloat + j];

                fwrite(&value, sizeof(double), 1, fp);
              }
            }
          }
        } else {
          long nc = Th.nv;
          unsigned nbytes = nc * sizeof(double) * (*this).nbfloat;
          int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);

          for (int iv = 0; iv < Th.nv; iv++) {
            for (int j = 0; j < (*this).nbfloat; j++) {
              valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
              double value = valsol[iv * (*this).nbfloat + j];

              l = runEncodeB64(sizeof(double), (unsigned char *)&value, ElementChars);
              ElementChars[l] = 0;
              fwrite(&ElementChars, l, 1, fp);
            }
          }

          // flush buffer
          l = runEncodeB64(0, NULL, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      } else {
        for (int iv = 0; iv < Th.nv; iv++) {
          for (int j = 0; j < (*this).nbfloat; j++) {
            valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
            double value = valsol[iv * (*this).nbfloat + j];

            fprintf(fp, "%.16e ", value);
          }
        }
      }

      fprintf(fp, "\n");
    }
  };
  vector< Expression2 > l;
#ifndef COMMON_HPDDM_PARALLEL_IO
  static const int n_name_param = 8;
#else
  static const int n_name_param = 9;
#endif
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

 public:
  VTK_WriteMesh_Op(const basicAC_F0 &args) : l(args.size( ) - 2) {
    int nbofsol;
    int ddim = 2;
    int stsize = 3;
    int sca = 0, vec = 0, ten = 0;
    string scas("scalaire");
    string vecs("vector");
    string tens("tensor");

    if (verbosity > 2) {
      cout << "Write Mesh and Solutions in VTK Formats" << endl;
    }

    args.SetNameParam(n_name_param, name_param, nargs);

    if (BCastTo< string * >(args[0])) {
      filename = CastTo< string * >(args[0]);
    }

    if (BCastTo< pmesh >(args[1])) {
      eTh = CastTo< pmesh >(args[1]);
    }

    nbofsol = l.size( );

    for (size_t i = 2; i < args.size( ); i++) {
      size_t jj = i - 2;

      if (BCastTo< double >(args[i])) {
        l[jj].what = 1;
        l[jj].nbfloat = 1;
        l[jj][0] = to< double >(args[i]);

        char number[16];
        sprintf(number, "%li", jj + 1);
        l[jj].name = scas;
        l[jj].name += number;
        sca++;
      } else if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a0 = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        // cout << "taille" << a0->size() << endl;
        if (a0->size( ) != ddim && a0->size( ) != stsize) {
          CompileError(
            "savesol in 2D: vector solution is 2 composant, tensor solution is 3 composant");
        }

        if (a0->size( ) == ddim) {
          // vector solution
          l[jj].what = 2;
          l[jj].nbfloat = ddim;

          for (int j = 0; j < ddim; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }

          char number[16];
          sprintf(number, "%li", jj + 1);
          l[jj].name = vecs;
          l[jj].name += number;
          vec++;
        } else if (a0->size( ) == stsize) {
          // symmetric tensor solution
          l[jj].what = 3;
          l[jj].nbfloat = stsize;

          for (int j = 0; j < stsize; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }

          char number[16];
          sprintf(number, "%li", jj + 1);
          l[jj].name = tens;
          l[jj].name += number;
          ten++;
        }
      } else {
        cout << " arg " << i << " " << args[i].left( ) << endl;
        CompileError("save solution in 2D in format VTK: Sorry no way to save this kind of data");
      }
    }
  }

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< string * >( ), atype< pmesh >( ), true);
  }    // all type

  static E_F0 *f(const basicAC_F0 &args) { return new VTK_WriteMesh_Op(args); }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type VTK_WriteMesh_Op::name_param[] = {{"dataname", &typeid(string *)},
                                                            {"withsurfacemesh", &typeid(bool)},
                                                            {"order", &typeid(KN_< long >)},
                                                            // A rajouter dans le 3D
                                                            {"floatmesh", &typeid(bool)},
                                                            {"floatsol", &typeid(bool)},
                                                            {"bin", &typeid(bool)},
                                                            {"swap", &typeid(bool)},
#ifdef COMMON_HPDDM_PARALLEL_IO
                                                            {"communicator", &typeid(pcommworld)},
#endif
                                                            {"append", &typeid(bool)}};

void VTK_WRITE_MESH(const string &filename, FILE *fp, const Mesh &Th, bool binary, int datasize,
                    bool surface, bool bigEndian) {
  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "%s, Created by Freefem++ \n", filename.c_str( ));
  if (binary) {
    fprintf(fp, "BINARY\n");
  } else {
    fprintf(fp, "ASCII\n");
  }

  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
  // get all the entities in the model

  // write mesh vertices

  if (datasize == sizeof(float)) {
    fprintf(fp, "POINTS %d float\n", Th.nv);

    for (unsigned int i = 0; i < Th.nv; i++) {
      const Mesh::Vertex &P = Th.vertices[i];
      float f[3];
      f[0] = P.x;
      f[1] = P.y;
      f[2] = 0.;    // P.z; 3D case
      if (binary) {
        if (!bigEndian) {
          FreeFEM::SwapBytes((char *)&f, sizeof(float), 3);
        }

        fwrite(&f, sizeof(float), 3, fp);
      } else {
        fprintf(fp, "%.8g %.8g %.8g\n", P.x, P.y, 0.);
      }
    }
  } else if (datasize == sizeof(double)) {
    fprintf(fp, "POINTS %d double\n", Th.nv);

    for (unsigned int i = 0; i < Th.nv; i++) {
      const Mesh::Vertex &P = Th.vertices[i];
      double f[3];
      f[0] = P.x;
      f[1] = P.y;
      f[2] = 0.;    // P.z; 3D case
      if (binary) {
        if (!bigEndian) {
          FreeFEM::SwapBytes((char *)&f, sizeof(double), 3);
        }

        fwrite((unsigned char *)&f, sizeof(double), 3, fp);
      } else {
        fprintf(fp, "%.15lg %.15lg %.15lg\n", f[0], f[1], f[2]);
      }
    }
  }

  fprintf(fp, "\n");
  if (verbosity > 1) {
    printf("writing vertices is finish\n");
  }

  if (verbosity > 1) {
    printf("writing elements now\n");
  }

  //= ===============
  // CELL
  //= ===============
  // loop over all elements we need to save and count vertices
  int numElements, totalNumInt;
  if (surface) {
    numElements = Th.nt + Th.neb;
    totalNumInt = Th.nt * 3 + Th.neb * 2 + numElements;
  } else {
    numElements = Th.nt;
    totalNumInt = Th.nt * 3 + numElements;
  }

  if (verbosity > 1) {
    printf("writing cells \n");
  }

  // print vertex indices in ascii or binary
  fprintf(fp, "CELLS %d %d\n", numElements, totalNumInt);

  if (binary) {
    int IntType = 3;
    if (verbosity > 1) {
      printf("writing triangle elements \n");
    }

    for (int it = 0; it < Th.nt; it++) {
      const Mesh::Triangle &K(Th.t(it));
      int iv[IntType + 1];

      iv[0] = IntType;

      for (int ii = 0; ii < IntType; ii++) {
        iv[ii + 1] = Th.operator( )(K[ii]);
      }

      if (!bigEndian) {
        FreeFEM::SwapBytes((char *)&iv, sizeof(int), IntType + 1);
      }

      fwrite(&iv, sizeof(int), IntType + 1, fp);
    }

    if (surface) {
      if (verbosity > 1) {
        printf("writing edge elements \n");
      }

      IntType = 2;

      for (int ibe = 0; ibe < Th.neb; ibe++) {
        const Mesh::BorderElement &K(Th.be(ibe));
        int iv[IntType + 1];
        iv[0] = IntType;

        for (int ii = 0; ii < IntType; ii++) {
          iv[ii + 1] = Th.operator( )(K[ii]);
        }

        if (!bigEndian) {
          FreeFEM::SwapBytes((char *)&iv, sizeof(int), IntType + 1);
        }

        fwrite(&iv, sizeof(int), IntType + 1, fp);
      }
    }
  } else {
    int IntType = 3;
    if (verbosity > 1) {
      printf("writing triangle elements \n");
    }

    for (int it = 0; it < Th.nt; it++) {
      const Mesh::Triangle &K(Th.t(it));
      int iv[IntType + 1];
      iv[0] = IntType;

      for (int ii = 0; ii < IntType; ii++) {
        iv[ii + 1] = Th.operator( )(K[ii]);
      }

      fprintf(fp, "%d %d %d %d\n", iv[0], iv[1], iv[2], iv[3]);
    }

    if (surface) {
      if (verbosity > 1) {
        printf("writing edge elements \n");
      }

      IntType = 2;

      for (int ibe = 0; ibe < Th.neb; ibe++) {
        const Mesh::BorderElement &K(Th.be(ibe));
        int iv[IntType + 1];
        iv[0] = IntType;

        for (int ii = 0; ii < IntType; ii++) {
          iv[ii + 1] = Th.operator( )(K[ii]);
        }

        fprintf(fp, "%d %d %d\n", iv[0], iv[1], iv[2]);
      }
    }
  }

  fprintf(fp, "\n");

  // CELL_TYPE
  // print element types in ascii or binary

  fprintf(fp, "CELL_TYPES %d\n", numElements);
  if (binary) {
    int type;

    for (int it = 0; it < Th.nt; it++) {
      type = VTK_TRI;
      if (!bigEndian) {
        FreeFEM::SwapBytes((char *)&type, sizeof(int), 1);
      }

      fwrite(&type, sizeof(int), 1, fp);
    }

    if (surface) {
      for (int ibe = 0; ibe < Th.neb; ibe++) {
        type = VTK_EDGE;
        if (!bigEndian) {
          FreeFEM::SwapBytes((char *)&type, sizeof(int), 1);
        }

        fwrite(&type, sizeof(int), 1, fp);
      }
    }
  } else {
    int type;
    type = VTK_TRI;

    for (int it = 0; it < Th.nt; it++) {
      fprintf(fp, "%d ", type);
    }

    if (surface) {
      type = VTK_EDGE;

      for (int ibe = 0; ibe < Th.neb; ibe++) {
        // fprintf(fp,"%d%c",type,(ibe%10==9)? '\n' : ' ');
        fprintf(fp, "%d ", type);
      }
    }
  }

  fprintf(fp, "\n");

  //= ================================
  // WRITE SOLUTION IN FORMAT VTK
  // LABEL OF ELEMENTS
  //= ================================

  list< int > list_label_Elem;
  // list<int> list_label_Border_Elem;
  {
    for (int it = 0; it < Th.nt; it++) {
      const Mesh::Triangle &K(Th.t(it));
      list< int >::const_iterator ilist;
      int labOk = 0;

      for (ilist = list_label_Elem.begin( ); ilist != list_label_Elem.end( ); ilist++) {
        if (*ilist == K.lab) {
          labOk = 1;
          break;
        }
      }

      if (labOk == 0) {
        list_label_Elem.push_back(K.lab);
      }
    }

    if (surface) {
      for (int ibe = 0; ibe < Th.neb; ibe++) {
        const Mesh::BorderElement &K(Th.be(ibe));
        list< int >::const_iterator ilist;
        int labOk = 0;

        for (ilist = list_label_Elem.begin( ); ilist != list_label_Elem.end( ); ilist++) {
          if (*ilist == K.lab) {
            labOk = 1;
            break;
          }
        }

        if (labOk == 0) {
          list_label_Elem.push_back(K.lab);
        }
      }
    }
  }
  list_label_Elem.sort( );

  //= ================================
  //= ================================

  fprintf(fp, "CELL_DATA %d\n", numElements);
  int cell_fd = 1;
  int cell_lab = 1;

  fprintf(fp, "Scalars  Label int %d\n", cell_fd);
  fprintf(fp, "LOOKUP_TABLE FreeFempp_table\n");
  // Determination des labels
  if (binary) {
    int label;

    for (int it = 0; it < Th.nt; it++) {
      const Mesh::Triangle &K(Th.t(it));
      label = K.lab;
      if (!bigEndian) {
        FreeFEM::SwapBytes((char *)&label, sizeof(int), 1);
      }

      fwrite(&label, sizeof(int), 1, fp);
    }

    if (surface) {
      for (int ibe = 0; ibe < Th.neb; ibe++) {
        const Mesh::BorderElement &K(Th.be(ibe));
        label = K.lab;
        if (!bigEndian) {
          FreeFEM::SwapBytes((char *)&label, sizeof(int), 1);
        }

        fwrite(&label, sizeof(int), 1, fp);
      }
    }
  } else {
    int label;

    for (int it = 0; it < Th.nt; it++) {
      const Mesh::Triangle &K(Th.t(it));
      label = K.lab;
      fprintf(fp, "%d\n", label);
    }

    if (surface) {
      for (int ibe = 0; ibe < Th.neb; ibe++) {
        const Mesh::BorderElement &K(Th.be(ibe));
        label = K.lab;
        fprintf(fp, "%d\n", label);
      }
    }
  }

  fprintf(fp, "\n");

  int size_list = 0;
  list< int >::const_iterator ilist;

  for (ilist = list_label_Elem.begin( ); ilist != list_label_Elem.end( ); ilist++) {
    size_list++;
  }

  fprintf(fp, "LOOKUP_TABLE FreeFempp_table %d\n", size_list);
  {
    list< int >::const_iterator ilist;

    for (ilist = list_label_Elem.begin( ); ilist != list_label_Elem.end( ); ilist++) {
      if (binary) {
        int tab[4];
        tab[0] = (int)ColorTable[abs(*ilist) % NbColorTable][0] * 255;
        tab[1] = (int)ColorTable[abs(*ilist) % NbColorTable][1] * 255;
        tab[2] = (int)ColorTable[abs(*ilist) % NbColorTable][2] * 255;
        tab[3] = 255;

        for (int itab = 0; itab < 4; itab++) {
          char newvalue[sizeof(int)];
          int bid0 = sprintf(newvalue, "%s", (char *)&tab[itab]);
          fwrite(&newvalue, sizeof(unsigned char), 1, fp);
        }
      } else {
        float tab[4];
        tab[0] = ColorTable[abs(*ilist) % NbColorTable][0];
        tab[1] = ColorTable[abs(*ilist) % NbColorTable][1];
        tab[2] = ColorTable[abs(*ilist) % NbColorTable][2];
        tab[3] = 1.;
        fprintf(fp, "%.8f %.8f %.8f %.8f\n", tab[0], tab[1], tab[2], tab[3]);
      }
    }
  }
  fprintf(fp, "\n");
}

AnyType VTK_WriteMesh_Op::operator( )(Stack stack) const {
  string *pffname = GetAny< string * >((*filename)(stack));
  Mesh *pTh = GetAny< Mesh * >((*eTh)(stack));

  ffassert(pTh);
  Mesh &Th = *pTh;
  bool bigEndian = isBigEndian( );
  bool swap = bigEndian;
  bool binary = true;
  bool surface = true;
  bool floatmesh = false;
  bool floatsol = false;
  string *dataname;
  int nbofsol = l.size( );
  KN< int > order(nbofsol);

  char *nameofuser[nbofsol];

  for (int ii = 0; ii < nbofsol; ii++) {
    order[ii] = 0;
  }

  if (nargs[0]) {
    dataname = GetAny< string * >((*nargs[0])(stack));
  }

  if (nargs[1]) {
    surface = GetAny< bool >((*nargs[1])(stack));
  }

  if (nargs[2]) {
    order = GetAny< KN_< long > >((*nargs[2])(stack));
  }

  if (nargs[3]) {
    floatmesh = GetAny< bool >((*nargs[3])(stack));
  }

  if (nargs[4]) {
    floatsol = GetAny< bool >((*nargs[4])(stack));
  }

  if (nargs[5]) {
    binary = GetAny< bool >((*nargs[5])(stack));
  }

  if (nargs[6]) {
    swap = GetAny< bool >((*nargs[6])(stack));
  }
#ifdef COMMON_HPDDM_PARALLEL_IO
  parallelIO(pffname, nargs[7] ? (MPI_Comm *)GetAny< pcommworld >((*nargs[7])(stack)) : 0,
             nargs[8] && GetAny< bool >((*nargs[8])(stack)));
#endif

  int datasize = floatmesh ? sizeof(float) : sizeof(double);

  int datasizeSol = floatsol ? sizeof(float) : sizeof(double);

  int iii = 0;
  if (nargs[0]) {
    // char *data = newcopy(dataname->c_str());
    char *data = newcopy(dataname);
    char *name = strtok(data, " \t\n");

    nameofuser[iii] = newcopy(name);
    if (verbosity > 5) {
      cout << "   iovtk : value of iii  =" << iii << "  \"" << nameofuser[iii] << "\"\n";
    }

    iii++;
    {
      while ((name = strtok(NULL, " \t\n\0"))) {
        if (iii >= nbofsol) {
          if (verbosity > 5) {
            cout << "   iovtk : The number of data name is too large " << endl;
          }

          break;
        }

        nameofuser[iii] = newcopy(name);
        if (verbosity > 5) {
          cout << "   iovtk : value of iii  =" << iii << "  \"" << nameofuser[iii] << "\"\n";
        }

        iii++;
      }

      if (iii < nbofsol) {
        if (verbosity > 6) {
          cout << "   iovtk:  The number of data name is too small, we give default name " << endl;
        }
      }

      delete[] data;
    }
  }

  if (iii < nbofsol) {
    for (int iiii = iii; iiii < nbofsol; iiii++) {
      // char *dataff = new char[l[iii].name.size()+1];
      // strcpy(dataff, l[iii].name.c_str());
      nameofuser[iiii] = newcopy(l[iiii].name.c_str( ));    // dataff;
    }
  }

  FILE *fp = fopen((*pffname).c_str( ), "wb");
  // FILE *fp = fopen( newcopy(pffname), "wb");
  if (!fp) {
    cerr << "Unable to open file " << (*pffname).c_str( ) << endl;
    ExecError("error in reading vtk file");
  }

  // determination of number of order 0 et 1.
  int Norder0 = 0;

  for (int ii = 0; ii < nbofsol; ii++) {
    if (order[ii] == 0) {
      Norder0++;
    }
  }

  char *pch = newcopy(pffname);
  int VTK_FILE = 0;
  int ls = 0;
  int lll = strlen(pch);
  if (!strcmp(pch + lll - (ls = 4), ".vtk")) {
    VTK_FILE = 1;
  } else if (!strcmp(pch + lll - (ls = 4), ".vtu")) {
    VTK_FILE = 2;
  }

  if (verbosity) {
    cout << " " << pffname << " VTK_FILE " << VTK_FILE << endl;
  }

  if (VTK_FILE == 1) {
    // CAS VTK
    VTK_WRITE_MESH(*pffname, fp, Th, binary, datasize, surface, swap);

    // WRITE SOLUTIONS
    if (datasizeSol == sizeof(float)) {
      if (Norder0 > 0) {
        fprintf(fp, "FIELD FieldData %d\n", Norder0);

        for (int ii = 0; ii < nbofsol; ii++) {
          if (order[ii] == 0) {
            int nsol;
            if (surface) {
              nsol = Th.nt + Th.neb;
            } else {
              nsol = Th.nt;
            }

            fprintf(fp, "%s %ld %d float\n", nameofuser[ii], l[ii].nbfloat, nsol);

            // changement ecriture solution
            l[ii].writesolutionP0_float(fp, Th, stack, surface, binary, swap);
          }
        }
      }

      if (Norder0 < nbofsol) {
        fprintf(fp, "POINT_DATA %d\n", Th.nv);
        fprintf(fp, "FIELD FieldData %d\n", nbofsol - Norder0);

        for (int ii = 0; ii < nbofsol; ii++) {
          if (order[ii] == 1) {
            // fprintf(fp,"%s %d %d float\n",l[ii].name.c_str(),l[ii].nbfloat,Th.nv);
            fprintf(fp, "%s %ld %d float\n", nameofuser[ii], l[ii].nbfloat, Th.nv);
            if (verbosity > 5) {
              cout << "name of data(" << ii << ")=" << nameofuser[ii] << " " << l[ii].name << endl;
            }

            // changement ecriture solution
            l[ii].writesolutionP1_float(fp, Th, stack, binary, swap);
          }
        }
      }
    }

    if (datasizeSol == sizeof(double)) {
      if (Norder0 > 0) {
        fprintf(fp, "FIELD FieldData %d\n", Norder0);

        for (int ii = 0; ii < nbofsol; ii++) {
          if (order[ii] == 0) {
            int nsol;
            if (surface) {
              nsol = Th.nt + Th.neb;
            } else {
              nsol = Th.nt;
            }

            fprintf(fp, "%s %ld %d double\n", nameofuser[ii], l[ii].nbfloat, nsol);
            // changement ecriture solution
            l[ii].writesolutionP0_double(fp, Th, stack, surface, binary, swap);
          }
        }
      }

      if (Norder0 < nbofsol) {
        fprintf(fp, "POINT_DATA %d\n", Th.nv);
        fprintf(fp, "FIELD FieldData %d\n", nbofsol - Norder0);

        for (int ii = 0; ii < nbofsol; ii++) {
          if (order[ii] == 1) {
            fprintf(fp, "%s %ld %d double\n", nameofuser[ii], l[ii].nbfloat, Th.nv);
            if (verbosity > 5) {
              cout << "name of data(" << ii << ")=" << nameofuser[ii] << endl;
            }

            // changement ecriture solution
            l[ii].writesolutionP1_double(fp, Th, stack, binary, swap);
          }
        }
      }
    }
  } else if (VTK_FILE == 2) {
    /*
     * // Solution Order
     * // order 1
     * if( Norder0 != nbofsol){
     *   BEGINTYPE_VTU( fp, "PointData");
     *   for(int ii=0; ii< nbofsol; ii++){
     *     if(order[ii] == 1){
     *
     *       if(datasize == sizeof(float)){
     *         VTU_DATA_ARRAY( fp, "Float32", nameofuser[ii] , l[ii].nbfloat, binary);
     *         l[ii].writesolutionP1_float( fp, Th, stack, binary, swap);
     *       }
     *       else if(datasize == sizeof(double)) {
     *         VTU_DATA_ARRAY( fp, "Float64", nameofuser[ii] , l[ii].nbfloat, binary);
     *         l[ii].writesolutionP1_double( fp, Th, stack, binary, swap);
     *       }
     *
     *       ENDTYPE_VTU( fp, "DataArray");
     *     }
     *   }
     *   ENDTYPE_VTU( fp, "PointData");
     * }
     * // order 0
     * if( Norder0 > 0 ){
     *  BEGINTYPE_VTU( fp, "CellData");
     *  for(int ii=0; ii< nbofsol; ii++){
     *    if(order[ii] == 0){
     *
     *      if(datasize == sizeof(float)){
     *         VTU_DATA_ARRAY( fp, "Float32", nameofuser[ii] , l[ii].nbfloat, binary);
     *         l[ii].writesolutionP0_float( fp, Th, stack, surface, binary, swap);
     *       }
     *       else if(datasize == sizeof(double)) {
     *         VTU_DATA_ARRAY( fp, "Float64", nameofuser[ii] , l[ii].nbfloat, binary);
     *         l[ii].writesolutionP0_double( fp, Th, stack, surface, binary, swap);
     *       }
     *      ENDTYPE_VTU( fp, "DataArray");
     *    }
     *  }
     *  ENDTYPE_VTU( fp, "CellData");
     * }
     * long offsetsol=0;
     * bool encode64=0;
     */
    VTU_WRITE_MESH(fp, Th, binary, datasize, surface);
    // Solution Order
    // order 0
    if (Norder0 > 0) {
      for (int ii = 0; ii < nbofsol; ii++) {
        if (order[ii] == 0) {
          if (datasize == sizeof(float)) {
            VTU_DATA_ARRAY(fp, "Float32", nameofuser[ii], l[ii].nbfloat, binary);
            l[ii].writesolutionP0_float(fp, Th, stack, surface, binary, swap, 1);
          } else if (datasize == sizeof(double)) {
            VTU_DATA_ARRAY(fp, "Float64", nameofuser[ii], l[ii].nbfloat, binary);
            l[ii].writesolutionP0_double(fp, Th, stack, surface, binary, swap, 1);
          }

          ENDTYPE_VTU(fp, "DataArray");
        }
      }
    }
    ENDTYPE_VTU(fp, "CellData");
    // order 1
    if (Norder0 != nbofsol) {
      BEGINTYPE_VTU(fp, "PointData");

      for (int ii = 0; ii < nbofsol; ii++) {
        if (order[ii] == 1) {
          if (datasize == sizeof(float)) {
            VTU_DATA_ARRAY(fp, "Float32", nameofuser[ii], l[ii].nbfloat, binary);
            l[ii].writesolutionP1_float(fp, Th, stack, binary, swap, 1);
          } else if (datasize == sizeof(double)) {
            VTU_DATA_ARRAY(fp, "Float64", nameofuser[ii], l[ii].nbfloat, binary);
            l[ii].writesolutionP1_double(fp, Th, stack, binary, swap, 1);
          }

          ENDTYPE_VTU(fp, "DataArray");
        }
      }

      ENDTYPE_VTU(fp, "PointData");
    }

    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
  } else {
    cout << " iovtk extension file is not correct (" << VTK_FILE << " != 1 or 2 ) " << endl;
    ExecError(" iovtk : extension file");
  }

  // close file fp
  fclose(fp);

  for (int iiii = 0; iiii < nbofsol; iiii++) {
    delete[] nameofuser[iiii];
  }

  delete[] pch;
  return (Mesh *)NULL;
}

//= =============================================
// FIN ECRITURE DE FICHIER .vtk (2D)
//= =============================================

//= ======================
// FIN 2D Fichier .vtk
//= ======================

//= =============================================
// LOAD DE FICHIER .vtk (3D) FOR MESH3
//= =============================================

class VTK_LoadMesh3_Op : public E_F0mps {
 public:
  Expression filename;
  static const int n_name_param = 8;    //
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  int arg(int i, Stack stack, int a) const {
    return nargs[i] ? GetAny< int >((*nargs[i])(stack)) : a;
  }
  bool arg(int i, Stack stack, bool a) const {
    return nargs[i] ? GetAny< bool >((*nargs[i])(stack)) : a;
  }
  double arg(int i, Stack stack, double a) const {
    return nargs[i] ? GetAny< double >((*nargs[i])(stack)) : a;
  }
   KN<KN<double> >* arg(int i, Stack stack,KN<KN<double> >*p) const {
        return nargs[i] ? GetAny< KN<KN<double> >* >((*nargs[i])(stack)) : p;
    }

 public:
  VTK_LoadMesh3_Op(const basicAC_F0 &args, Expression ffname) : filename(ffname) {
    if (verbosity) {
      cout << "Load mesh given by VTK " << endl;
    }

    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type VTK_LoadMesh3_Op::name_param[] = {
  {"reftet", &typeid(long)},         {"swap", &typeid(bool)},
  {"refface", &typeid(long)},        {"namelabel", &typeid(string)},
  {"cleanmesh", &typeid(bool)},      {"removeduplicate", &typeid(bool)},
  {"precisvertice", &typeid(double)},
  {"fields", &typeid( KN<KN<double> >*)}
    
};

class VTK_LoadMesh3 : public OneOperator {
 public:
  VTK_LoadMesh3( ) : OneOperator(atype< pmesh3 >( ), atype< string * >( )) {}

  E_F0 *code(const basicAC_F0 &args) const {
    return new VTK_LoadMesh3_Op(args, t[0]->CastTo(args[0]));
  }
};

Mesh3 *VTK_Load3(const string &filename, bool bigEndian, bool cleanmesh, bool removeduplicate,
                 double precisvertice,KN<KN<double> >* pfields) {
  // swap = bigEndian or not bigEndian
  // variable freefem++
  int nv, nt = 0, nbe = 0;
  int nerr = 0;
  char *res;
  // Reading Mesh in vtk formats
  FILE *fp = fopen(filename.c_str( ), "rb");

  if (!fp) {
    cerr << "Unable to open file " << filename.c_str( ) << endl;
    ExecError("error in reading vtk file");
  }

  char buffer[256], buffer2[256],buffer3[256];

  res = fgets(buffer, sizeof(buffer), fp);    // version line
  res = fgets(buffer, sizeof(buffer), fp);    // title

  fscanf(fp, "%s", buffer);    // ASCII or BINARY
  bool binary = false;
  if (!strcmp(buffer, "BINARY")) {
    binary = true;
  }
    if(verbosity>9)  cout << " binary = " <<binary << " bigEndian: " << bigEndian  <<endl;
  if (fscanf(fp, "%s %s", buffer, buffer2) != 2) {
    cout << "error in reading vtk files" << endl;
    ExecError("error in reading vtk file");
  }

  if (strcmp(buffer, "DATASET") || strcmp(buffer2, "UNSTRUCTURED_GRID")) {
    cout << "VTK reader can only read unstructured datasets" << endl;
    ExecError("error in reading vtk file");
  }

  // read mesh vertices
  if (fscanf(fp, "%s %d %s\n", buffer, &nv, buffer2) != 3) {
    cout << "error in reading vtk files" << endl;
    ExecError("error in reading vtk file");
  }

  if (strcmp(buffer, "POINTS") || !nv) {
    cerr << "No points in dataset" << endl;
    ExecError("error in reading vtk file");
  }

  int datasize;
  if (!strncmp(buffer2, "double", 6)) {
    datasize = sizeof(double);
  } else if (!strncmp(buffer2, "float", 5)) {
    datasize = sizeof(float);
  } else {
    cout << "VTK reader only accepts float or double datasets" << endl;
    ExecError("error in reading vtk file");
  }

  if (verbosity > 3) {
    cout << "Reading points" << nv << " buffer2 " << buffer2 << " binary " << binary << " datasize "
         << datasize << " " << sizeof(float) << endl;
  }

  Vertex3 *vff = new Vertex3[nv];
  
  for (int i = 0; i < nv; i++) {
    if (verbosity > 19) {
      cout << " i=" << i << endl;
    }

    double xyz[3];
    if (binary) {
      if (datasize == sizeof(float)) {
        float f[3];
        if (fread(f, sizeof(float), 3, fp) != 3) {
            cout << "error in reading vtk files item: " << i << endl;
          ExecError("error in reading vtk file");
        }

        if (!bigEndian) {
          FreeFEM::SwapBytes((char *)f, sizeof(float), 3);
        }

        for (int j = 0; j < 3; j++) {
          xyz[j] = f[j];
        }
      } else {
        if (fread(xyz, sizeof(double), 3, fp) != 3) {
          cout << "error in reading vtk files item " <<i<<  endl;
          ExecError("error in reading vtk file");
        }
        //     cout <<  i << " : " << xyz[0] << " " <<  xyz[1] << " " <<  xyz[2] << endl;
        if (!bigEndian) {
          FreeFEM::SwapBytes((char *)xyz, sizeof(double), 3);
        //     cout <<  i << " : " << xyz[0] << " " <<  xyz[1] << " " <<  xyz[2] << endl;
        }
      }
    } else {
      if (verbosity > 10) cout <<"  iovtk: datasize " << datasize << " == " << sizeof(float) << endl;
      if (datasize == sizeof(float)) {
        if (fscanf(fp, "%lf %lf %lf", &xyz[0], &xyz[1], &xyz[2]) != 3) {
          cout << "error in reading vtk files (float)" << endl;
          ExecError("error in reading vtk file");
        }
      } else {
        if (fscanf(fp, "%lf %lf %lf", &xyz[0], &xyz[1], &xyz[2]) != 3) {
          cout << "error in reading vtk files (double)" << endl;
          ExecError("error in reading vtk file");
        }
      }
    }

    vff[i].x = xyz[0];
    vff[i].y = xyz[1];
    vff[i].z = xyz[2];
    vff[i].lab = 1;
    if (verbosity > 9) {
        cout <<  i << " : " << xyz[0] << " " <<  xyz[1] << " " <<  xyz[2] << endl;
    }
  }

  // read mesh elements
  int numElements, numElements2, totalNumInt,kk,fpos=ftell(fp);
  if ((kk=fscanf(fp, "%s %d %d\n", buffer, &numElements, &totalNumInt)) != 3) {
    cout << "error in " << fpos<< " " <<  buffer << " reading vtk files" << numElements << " " << totalNumInt << " " << kk<< endl;
    ExecError("error in reading vtk file");
  }

  if (verbosity > 3) {
      printf("reading %d parameter %s %d %d\n",fpos, buffer, numElements, totalNumInt);
  }

  if (strncmp(buffer, "CELLS", 5) || !numElements) {
    cout << "No cells in dataset" << endl;
    ExecError("error in reading vtk file");
  }

  if (verbosity > 3) {
    cout << "Reading cells" << numElements << endl;
  }

  int *IntCells = new int[totalNumInt - numElements];
  int *firstCell = new int[numElements + 1];
  int *TypeCells = new int[numElements];
  int numIntCells = 0;

  for (unsigned int i = 0; i < numElements; i++) {
    int numVerts, n[100];
    if (verbosity > 9) {
      cout << "i=" << i << " " << numElements << endl;
    }

    for (int ii = 0; ii < 100; ii++) {
      n[ii] = -1;
    }

    if (binary) {
      if (fread(&numVerts, sizeof(int), 1, fp) != 1) {
        cout << "error in reading VTK files " << endl;
        ExecError("error in reading vtk file");
      }

      if (!bigEndian) {
        FreeFEM::SwapBytes((char *)&numVerts, sizeof(int), 1);
      }

      if ((int)fread(n, sizeof(int), numVerts, fp) != numVerts) {
        cout << "error in reading VTK files " << endl;
        ExecError("error in reading vtk file");
      }

      if (!bigEndian) {
        FreeFEM::SwapBytes((char *)n, sizeof(int), numVerts);
      }
    } else {
      if (fscanf(fp, "%d", &numVerts) != 1) {
        cout << "error in reading VTK files " << endl;
        ExecError("error in reading vtk file");
      }

     // cout << "numVerts" << numVerts << endl;

      for (int j = 0; j < numVerts; j++) {
        if (fscanf(fp, "%d", &n[j]) != 1) {
          cout << "error in reading VTK files " << endl;
          ExecError("error in reading vtk file");
        }

        if (verbosity > 9) {
          cout << "  n[j]" << n[j] << endl;
        }
      }
    }

    firstCell[i] = numIntCells;

    for (int j = 0; j < numVerts; j++) {
      if (n[j] >= 0 && n[j] < nv) {
        IntCells[numIntCells] = n[j];
        numIntCells++;
      } else {
        cout << "Bad vertex index" << endl;
        ExecError("error in reading vtk file");
      }
    }
  }

  firstCell[numElements] = totalNumInt - numElements;

  if (fscanf(fp, "%s %d\n", buffer, &numElements2) != 2) {
    cout << " Error in reading CELL_TYPES ARGUMENT " << endl;
    ExecError("error in reading vtk file");
  }

  if (strcmp(buffer, "CELL_TYPES") || numElements2 != (int)numElements) {
    cout << "No or invalid number of cells types" << endl;
    ExecError("error in reading vtk file");
  }

  if (verbosity > 3) {
    printf("reading parameter %s %d\n", buffer, numElements2);
  }

  // 3D

  for (unsigned int i = 0; i < numElements; i++) {
    int type;
    if (binary) {
      if (fread(&type, sizeof(int), 1, fp) != 1) {
        cout << "bug in readings cell types" << endl;
        ExecError("error in reading vtk file");
      }

      if (!bigEndian) {
        FreeFEM::SwapBytes((char *)&type, sizeof(int), 1);
      }
    } else {
      if (fscanf(fp, "%d", &type) != 1) {
        cout << "bug in readings cell types" << endl;
        ExecError("error in reading vtk file");
      }
    }

    TypeCells[i] = type;

    switch (type) {
      case 1:    // Vertex
        if (nerr++ < 3 && verbosity) {
         if (verbosity) cout << "this type of cell (vertex) is not taking account in Freefem++ " << endl;
        }

        break;
      case 3:    // Edge/line
       if (verbosity) cout << "this type of cell is not taking account in Freefem++ for a two dimensional mesh"
             << endl;    // 3D
        break;
      case 5:     // Triangle
        nbe++;    // 3D
        break;
      case 10:    // Tetrah�dre
        nt++;
        break;
      default:
        cout << "Error :: This type of cell is not considered in Freefem++" << endl;
        ExecError("error in reading vtk file");
        break;
    }
  }

    int nbp=0,nbf=0, err=0;
    if (fscanf(fp, "%s %d", buffer, &nbp) != 2)
        {   cout << "error in reading vtk files pfields" << endl;
            err++;}
    int startdatapoint=0;
    if(err==0)
    {
        int nf=-1;
        /*
         CELL_DATA 209726
         Scalars  Label int 1
         LOOKUP_TABLE FreeFempp_table
         ....
         LOOKUP_TABLE FreeFempp_table 7
         4*7 value
         */
        
        if (strcmp(buffer, "CELL_DATA"))  { //  read region number if exist
            if (strcmp(buffer, "POINT_DATA"))  {
                if (verbosity)cout << "VTK reader can only read CELL_DATA or POINT_DATA datasets:  not " << buffer<< " " << nbp<<  endl;
                err=1;
            }
            else startdatapoint=1;
        }
        else {
            if ((!err) &&(fscanf(fp, "%s %s %s %d\n", buffer, buffer2,buffer3,&nbf) != 4)) {
                if (verbosity)cout << "error in reading vtk files FIELD FieldData" << endl;
                err++;
            }}
        if(startdatapoint==0)
        {
        if( strcmp(buffer3, "int") !=0)// not integer
            err++;
        if ((!err) &&(fscanf(fp, "%s %s\n", buffer, buffer2) != 2))
            err++;
        // read nbf
        if (verbosity)cout << " err= " << err << " read nbp "<< nbp << endl;
        if(err==0)
            for( nf=0 ; nf < nbp; nf++)
            {
                int ii[1];
                if (binary)
                { if (fread(ii, sizeof(int), 1, fp) != 1) err++;}
                else
                {      if (fscanf(fp, "%d", ii) != 1) err++;}
                if(err) break;
            }
        if(err) cout << " err reading CELL_DATA  at " << nf << endl;
        
        if ((!err) &&(fscanf(fp, "%s %s %d\n", buffer, buffer2,&nbf) != 3 ) ) err++;
        nf =-1;
        if(err==0)
            for( nf=0 ; nf < nbf; nf++)
            {
                float f[4];
                char cc[4];
                if (binary)
                { if (fread(cc, sizeof(char), 4, fp) != 4) err++;}
                else
                {      if (fscanf(fp, "%f %f %f %fa",f+0,f+1,f+2,f+3) != 4) err++;}
                if(err) break;
            }
        if(err&& nf>=0 && verbosity) cout << " err LOOKUP_TABLE FreeFempp_table at " << nf << " " << err << endl;
        
        startdatapoint=0;
        }
        }
           
           

  if(pfields && err==0) {
    if(verbosity>1)   cout << " try  reading POINT_DATA  " << startdatapoint << endl;
    /*
    POINT_DATA 32436
    FIELD FieldData 2
    Velocity 3 32436 float
      
    */
    nbp =0;//  no POINT_DATA
    nbf=0;
    if (startdatapoint==0)
    {
      if (fscanf(fp, "%s", buffer) != 1) {
        cout << "error in reading vtk files pfields" << endl;
        err++;
      }
        else {if (verbosity) cout << " buff: "<< buffer << nbp << endl;}
    }
      if (strcmp(buffer, "POINT_DATA")==0) {
        if(startdatapoint==0)
         if (fscanf(fp, "%d", &nbp) != 1) err++;
        if (fscanf(fp, "%s", buffer) != 1) err++;
      }
      if(err == 0 && strcmp(buffer, "FIELD")!=0) {
        cout << "VTK reader can only read FIELD/POINT_DATA datasets:  not " << buffer<<  endl;
        err++;
      }
        
      if ((!err) &&(fscanf(fp, "%s %d", buffer2,&nbf) != 2)) {
        cout << "error in reading vtk files FIELD FieldData" << endl;
        err++;
      }
        
      if( err) nbf=0;
      else pfields->resize(nbf);
      for(int nf=0 ; nf < nbf; nf++) {
                
        int m,nv;
        // read mesh vertices
        if (fscanf(fp, "%s %d %d %s\n", buffer,&m, &nv, buffer2) != 4) {
          cout << "error in reading vtk files " << endl;
          err++;
          break;
        }
        int n = m*nv;
        if(verbosity) cout << " reading "<< buffer << " "<< m << " " << nv << " "<< buffer2 << endl;
        int datasize;
        if (!strncmp(buffer2, "double", 6))
          datasize = sizeof(double);
        else if (!strncmp(buffer2, "float", 5))
           datasize = sizeof(float);
        else {
         if (verbosity) cout << "VTK reader only accepts float or double datasets" << endl;
          err++;
          break;
        }
        //  read data ..
        if(err) break;
        (*pfields)[nf].resize(n);
        double* pv=&(*pfields)[nf][0];
        for(int i=0; i<n; ++i) {
          if (binary) {
            double f[1];
            if (datasize == sizeof(float)) {
              if (fread(f, sizeof(float), 1, fp) != 1) {
                if (verbosity)cout << "error in reading  vtk fields float at " << nf << " / " << i << endl;
                err++;
                break;
              }
            }
            else {
              if (fread(f, sizeof(double), 1, fp) != 1) {
                if (verbosity)cout << "error in reading  vtk fields double at " << nf << " / " << i << endl;
                err++;
                break;
              }
            }
            if (!bigEndian)
              FreeFEM::SwapBytes((char *)f, sizeof(float), 1);
            *pv++=*f;
          }
          else
            if (fscanf(fp, "%lf", pv++) != 1) {
              if (verbosity)cout << "error in reading vtk files ascii fields" << nf << " / " << i <<endl;
              err++;
            }
          if(err) break;
        }
        if(err) break;
      }
    if( err ) (*pfields).resize(0);
  } // end if(pfields)
  fclose(fp);

  // 3D versions

  Tet *tff;
  if (nt > 0)
    tff = new Tet[nt];
  else
    ExecError(
      "error in reading vtk file: not Tetrahedrons in vtk file, if it's a surface mesh use "
      "vtkloadS function");
  Tet *ttff = tff;
  Triangle3 *bff = new Triangle3[nbe];
  Triangle3 *bbff = bff;
  int   badorient=0;
  for (unsigned int i = 0; i < numElements; i++) {
    int type = TypeCells[i];
    int ivb[3], ivt[4];
    int label = 1;

    switch (type) {
      case 1:    // Vertex
        if (nerr++ < 3 && verbosity) {
          cout << "this type of cell is not taking account in Freefem++ " << endl;
        }

        break;
      case 3:    // Edge/line
        break;
      case 5:    // Triangle
        if(verbosity>9) cout << i << " " << firstCell[i + 1] << " " << firstCell[i] << endl;
        assert((firstCell[i + 1] - firstCell[i]) == 3);

        for (int j = firstCell[i]; j < firstCell[i + 1]; j++) {
          ivb[j - firstCell[i]] = IntCells[j];
        }

        (bbff++)->set(vff, ivb, label);
        break;
      case 10:    // Tetrah�dre
        assert((firstCell[i + 1] - firstCell[i]) == 4);

        for (int j = firstCell[i]; j < firstCell[i + 1]; j++) {
          ivt[j - firstCell[i]] = IntCells[j];
        }

        (ttff)->set(vff, ivt, label);
            if(ttff->mesure() < 0) {
                badorient ++;
               
                if( verbosity && badorient <11 )
                    cout << "   ** " <<  badorient << "  bad Tet "<< (ttff-tff) << " : "
                    << ivt[0]<< " " <<ivt[1] << " " <<ivt[2]<<" " << ivt[3] << " Vol="<<ttff->mesure()<< endl;
             (ttff)->set(vff, ivt, label);
            }
        ++ttff;
        break;
      default:
        break;
    }
  }
    if( badorient ) {
        cout << "  iovtk:loadmesh3:  Fatal error some  Tet with bad oriantation (vol <0) :   nb = " <<badorient << endl;
        ffassert(0);
    }
  delete[] IntCells;
  delete[] firstCell;
  delete[] TypeCells;
/*
 Mesh3::Mesh3(int nnv, int nnt, int nnbe, Vertex3 *vv, Tet *tt, Triangle3 *bb, bool cleanmesh, bool removeduplicate, bool rebuildboundary, int orientation, double precis_mesh)
 
 */
    
  Mesh3 *pTh = new Mesh3(nv, nt, nbe, vff, tff, bff, cleanmesh || (nbe==0), removeduplicate, (nbe==0) , 0, precisvertice);
  return pTh;
}

AnyType VTK_LoadMesh3_Op::operator( )(Stack stack) const {
  string *pffname = GetAny< string * >((*filename)(stack));

  int reftetra(arg(0, stack, 1));
  bool swap(arg(1, stack, false));
  int reftri(arg(2, stack, 1));

  string *DataLabel;
  if (nargs[3]) {
    DataLabel = GetAny< string * >((*nargs[3])(stack));
  }

  bool cleanmesh(arg(4, stack, false));
  bool removeduplicate(arg(5, stack, false));
  double precisvertice(arg(6, stack, 1e-6));
    KN<KN<double> > * pfields=0;
    pfields=arg(7, stack,pfields);
  Mesh3 *Th = VTK_Load3(*pffname, swap, cleanmesh, removeduplicate, precisvertice,pfields);

  // A faire fonction pour changer le label

  Add2StackOfPtr2FreeRC(stack, Th);

  return Th;
}

//= =============================================
// ECRITURE DE FICHIER .vtk (3D) FOR MESH3
//= =============================================

class VTK_WriteMesh3_Op : public E_F0mps {
 public:
  typedef long Result;
  Expression eTh;
  Expression filename;
  struct Expression2 {
    string name;
    long what;       // 1 scalar, 2 vector, 3 symtensor
    long nbfloat;    // 1 scalar, 2 vector (3D), 3 symtensor(3D)
    Expression e[6];
    Expression2( ) {
      e[0] = 0;
      e[1] = 0;
      e[2] = 0;
      e[3] = 0;
      e[4] = 0;
      e[5] = 0;
      what = 0;
      nbfloat = 0;
    };
    Expression &operator[](int i) { return e[i]; }

    double eval(int i, Stack stack) const {
      if (e[i]) {
        return GetAny< double >((*e[i])(stack));
      } else {
        return 0;
      }
    }

    void writesolutionP0_float_XML(FILE *fp, const Mesh3 &Th, Stack stack, bool surface) const {
      unsigned char ElementChars[256];
      MeshPoint *mp3(MeshPointStack(stack));
      R3 Cdg_hat = R3(1. / 4., 1. / 4., 1 / 4.);
      long nc = Th.nt;

      if (surface) {
        nc = nc + Th.nbe;
      }

      unsigned nbytes = nc * sizeof(float) * (*this).nbfloat;
      int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      for (int it = 0; it < Th.nt; it++) {
        const Tet &K(Th.elements[it]);
        mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

        for (int j = 0; j < nbfloat; j++) {
          float value = eval(j, stack);

          l = runEncodeB64(sizeof(float), (unsigned char *)&value, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      }

      if (surface) {
        for (int ibe = 0; ibe < Th.nbe; ibe++) {
          // determination du triangle contenant cette edge
          int ie;
          int it = Th.BoundaryElement(ibe, ie);
          const Tet &K(Th.elements[it]);
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < nbfloat; j++) {
            float value = eval(j, stack);

            l = runEncodeB64(sizeof(float), (unsigned char *)&value, ElementChars);
            ElementChars[l] = 0;
            fwrite(&ElementChars, l, 1, fp);
          }
        }
      }

      // flush buffer
      l = runEncodeB64(0, NULL, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);
      fprintf(fp, "\n");
    }

    void writesolutionP0_float_binary(FILE *fp, const Mesh3 &Th, Stack stack, bool surface,
                                      bool bigEndian) const {
      MeshPoint *mp3(MeshPointStack(stack));
      R3 Cdg_hat = R3(1. / 4., 1. / 4., 1 / 4.);

      if (!bigEndian) {
        for (int it = 0; it < Th.nt; it++) {
          const Tet &K(Th.elements[it]);
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < nbfloat; j++) {
            float value = eval(j, stack);

            FreeFEM::SwapBytes((char *)&value, sizeof(float), 1);
            fwrite(&value, sizeof(float), 1, fp);
          }
        }

        if (surface) {
          for (int ibe = 0; ibe < Th.nbe; ibe++) {
            // determination du triangle contenant cette edge
            int ie;
            int it = Th.BoundaryElement(ibe, ie);
            const Tet &K(Th.elements[it]);
            mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

            for (int j = 0; j < nbfloat; j++) {
              float value = eval(j, stack);

              FreeFEM::SwapBytes((char *)&value, sizeof(float), 1);
              fwrite(&value, sizeof(float), 1, fp);
            }
          }
        }
      } else {
        for (int it = 0; it < Th.nt; it++) {
          const Tet &K(Th.elements[it]);
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < nbfloat; j++) {
            float value = eval(j, stack);

            fwrite(&value, sizeof(float), 1, fp);
          }
        }

        if (surface) {
          for (int ibe = 0; ibe < Th.nbe; ibe++) {
            // determination du triangle contenant cette edge
            int ie;
            int it = Th.BoundaryElement(ibe, ie);
            const Tet &K(Th.elements[it]);
            mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

            for (int j = 0; j < nbfloat; j++) {
              float value = eval(j, stack);

              fwrite(&value, sizeof(float), 1, fp);
            }
          }
        }
      }

      fprintf(fp, "\n");
    }

    void writesolutionP0_float(FILE *fp, const Mesh3 &Th, Stack stack, bool surface) const {
      MeshPoint *mp3(MeshPointStack(stack));
      R3 Cdg_hat = R3(1. / 4., 1. / 4., 1 / 4.);

      for (int it = 0; it < Th.nt; it++) {
        const Tet &K(Th.elements[it]);
        mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

        for (int j = 0; j < nbfloat; j++) {
          float value = eval(j, stack);

          fprintf(fp, "%.8e ", value);
        }
      }

      if (surface) {
        for (int ibe = 0; ibe < Th.nbe; ibe++) {
          // determination du triangle contenant cette edge
          int ie;
          int it = Th.BoundaryElement(ibe, ie);
          const Tet &K(Th.elements[it]);
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < nbfloat; j++) {
            float value = eval(j, stack);

            fprintf(fp, "%.8e ", value);
          }
        }
      }

      fprintf(fp, "\n");
    }

    void writesolutionP0_float(FILE *fp, const Mesh3 &Th, Stack stack, bool surface, bool binary,
                               bool bigEndian, bool XML = false) const {
      if (binary) {
        if (!XML) {
          (*this).writesolutionP0_float_binary(fp, Th, stack, surface, bigEndian);
        } else {
          (*this).writesolutionP0_float_XML(fp, Th, stack, surface);
        }
      } else {
        (*this).writesolutionP0_float(fp, Th, stack, surface);
      }
    }

    void writesolutionP0_double_binary(FILE *fp, const Mesh3 &Th, Stack stack, bool surface,
                                       bool bigEndian) const {
      MeshPoint *mp3(MeshPointStack(stack));
      R3 Cdg_hat = R3(1. / 4., 1. / 4., 1. / 4.);

      if (!bigEndian) {
        for (int it = 0; it < Th.nt; it++) {
          const Tet &K(Th.elements[it]);
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < (*this).nbfloat; j++) {
            double value = (*this).eval(j, stack);

            FreeFEM::SwapBytes((char *)&value, sizeof(double), 1);
            fwrite(&value, sizeof(double), 1, fp);
          }
        }

        if (surface) {
          for (int ibe = 0; ibe < Th.nbe; ibe++) {
            // determination du triangle contenant cette edge
            int ie;
            int it = Th.BoundaryElement(ibe, ie);
            const Tet &K(Th.elements[it]);
            mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

            for (int j = 0; j < (*this).nbfloat; j++) {
              double value = (*this).eval(j, stack);

              FreeFEM::SwapBytes((char *)&value, sizeof(double), 1);
              fwrite(&value, sizeof(double), 1, fp);
            }
          }
        }
      } else {
        for (int it = 0; it < Th.nt; it++) {
          const Tet &K(Th.elements[it]);
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < (*this).nbfloat; j++) {
            double value = (*this).eval(j, stack);

            fwrite(&value, sizeof(double), 1, fp);
          }
        }

        if (surface) {
          for (int ibe = 0; ibe < Th.nbe; ibe++) {
            // determination du triangle contenant cette edge
            int ie;
            int it = Th.BoundaryElement(ibe, ie);
            const Tet &K(Th.elements[it]);
            mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

            for (int j = 0; j < (*this).nbfloat; j++) {
              double value = (*this).eval(j, stack);

              fwrite(&value, sizeof(double), 1, fp);
            }
          }
        }
      }

      fprintf(fp, "\n");
    }

    void writesolutionP0_double_XML(FILE *fp, const Mesh3 &Th, Stack stack, bool surface) const {
      unsigned char ElementChars[256];
      MeshPoint *mp3(MeshPointStack(stack));
      R3 Cdg_hat = R3(1. / 4., 1. / 4., 1. / 4.);
      long nc = Th.nt;

      if (surface) {
        nc = nc + Th.nbe;
      }

      unsigned nbytes = nc * sizeof(double) * (*this).nbfloat;
      int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      for (int it = 0; it < Th.nt; it++) {
        const Tet &K(Th.elements[it]);
        mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

        for (int j = 0; j < (*this).nbfloat; j++) {
          double value = (*this).eval(j, stack);

          l = runEncodeB64(sizeof(double), (unsigned char *)&value, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      }

      if (surface) {
        for (int ibe = 0; ibe < Th.nbe; ibe++) {
          // determination du triangle contenant cette edge
          int ie;
          int it = Th.BoundaryElement(ibe, ie);
          const Tet &K(Th.elements[it]);
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < (*this).nbfloat; j++) {
            double value = (*this).eval(j, stack);

            l = runEncodeB64(sizeof(double), (unsigned char *)&value, ElementChars);
            ElementChars[l] = 0;
            fwrite(&ElementChars, l, 1, fp);
          }
        }
      }

      // flush buffer
      l = runEncodeB64(0, NULL, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);
      fprintf(fp, "\n");
    }

    void writesolutionP0_double(FILE *fp, const Mesh3 &Th, Stack stack, bool surface) const {
      MeshPoint *mp3(MeshPointStack(stack));
      R3 Cdg_hat = R3(1. / 4., 1. / 4., 1. / 4.);

      for (int it = 0; it < Th.nt; it++) {
        const Tet &K(Th.elements[it]);
        mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

        for (int j = 0; j < (*this).nbfloat; j++) {
          double value = (*this).eval(j, stack);
          fprintf(fp, "%.16e ", value);
        }
      }

      if (surface) {
        for (int ibe = 0; ibe < Th.nbe; ibe++) {
          // determination du triangle contenant cette edge
          int ie;
          int it = Th.BoundaryElement(ibe, ie);
          const Tet &K(Th.elements[it]);
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < (*this).nbfloat; j++) {
            double value = (*this).eval(j, stack);

            fprintf(fp, "%.16e ", value);
          }
        }
      }

      // fprintf(fp,"\n");
    }

    void writesolutionP0_double(FILE *fp, const Mesh3 &Th, Stack stack, bool surface, bool binary,
                                bool bigEndian, bool XML = false) const {
      if (binary) {
        if (!XML) {
          (*this).writesolutionP0_double_binary(fp, Th, stack, surface, bigEndian);
        } else {
          (*this).writesolutionP0_double_XML(fp, Th, stack, surface);
        }
      } else {
        (*this).writesolutionP0_double(fp, Th, stack, surface);
      }
    }

    void writesolutionP1_float(FILE *fp, const Mesh3 &Th, Stack stack, bool binary, bool bigEndian,
                               bool XML = false) const {
      unsigned char ElementChars[256];
      MeshPoint *mp3(MeshPointStack(stack));

      KN< double > valsol(Th.nv * (*this).nbfloat);
      KN< int > takemesh(Th.nv);
      takemesh = 0;
      valsol = 0.;

      for (int it = 0; it < Th.nt; it++) {
        for (int iv = 0; iv < 4; iv++) {
          int i = Th(it, iv);
          mp3->setP(&Th, it, iv);

          for (int j = 0; j < (*this).nbfloat; j++) {
            valsol[i * (*this).nbfloat + j] =
              valsol[i * (*this).nbfloat + j] + (*this).eval(j, stack);
          }

          takemesh[i] = takemesh[i] + 1;
        }
      }

      if (binary) {
        if (!XML) {
          if (!bigEndian) {
            for (int iv = 0; iv < Th.nv; iv++) {
              for (int j = 0; j < (*this).nbfloat; j++) {
                valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
                float value = valsol[iv * (*this).nbfloat + j];

                if (!bigEndian) {
                  FreeFEM::SwapBytes((char *)&value, sizeof(float), 1);
                }

                fwrite(&value, sizeof(float), 1, fp);
              }
            }
          } else {
            for (int iv = 0; iv < Th.nv; iv++) {
              for (int j = 0; j < (*this).nbfloat; j++) {
                valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
                float value = valsol[iv * (*this).nbfloat + j];

                fwrite(&value, sizeof(float), 1, fp);
              }
            }
          }
        } else {
          long nc = Th.nv;
          unsigned nbytes = nc * sizeof(float) * (*this).nbfloat;
          int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);

          for (int iv = 0; iv < Th.nv; iv++) {
            for (int j = 0; j < (*this).nbfloat; j++) {
              valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
              float value = valsol[iv * (*this).nbfloat + j];

              l = runEncodeB64(sizeof(float), (unsigned char *)&value, ElementChars);
              ElementChars[l] = 0;
              fwrite(&ElementChars, l, 1, fp);
            }
          }

          // flush buffer
          l = runEncodeB64(0, NULL, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      } else {
        for (int iv = 0; iv < Th.nv; iv++) {
          for (int j = 0; j < (*this).nbfloat; j++) {
            valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
            float value = valsol[iv * (*this).nbfloat + j];

            fprintf(fp, " %.8e\n", value);
          }
        }
      }

      fprintf(fp, "\n");
    }

    void writesolutionP1_double(FILE *fp, const Mesh3 &Th, Stack stack, bool binary, bool bigEndian,
                                bool XML = false) const {
      unsigned char ElementChars[256];
      MeshPoint *mp3(MeshPointStack(stack));

      KN< double > valsol(Th.nv * (*this).nbfloat);
      KN< int > takemesh(Th.nv);
      takemesh = 0;
      valsol = 0.;

      for (int it = 0; it < Th.nt; it++) {
        for (int iv = 0; iv < 4; iv++) {
          int i = Th(it, iv);
          mp3->setP(&Th, it, iv);

          for (int j = 0; j < (*this).nbfloat; j++) {
            valsol[i * (*this).nbfloat + j] =
              valsol[i * (*this).nbfloat + j] + (*this).eval(j, stack);
          }

          takemesh[i] = takemesh[i] + 1;
        }
      }

      if (binary) {
        if (!XML) {
          if (!bigEndian) {
            for (int iv = 0; iv < Th.nv; iv++) {
              for (int j = 0; j < (*this).nbfloat; j++) {
                valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
                double value = valsol[iv * (*this).nbfloat + j];

                FreeFEM::SwapBytes((char *)&value, sizeof(double), 1);
                fwrite(&value, sizeof(double), 1, fp);
              }
            }
          } else {
            for (int iv = 0; iv < Th.nv; iv++) {
              for (int j = 0; j < (*this).nbfloat; j++) {
                valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
                double value = valsol[iv * (*this).nbfloat + j];

                fwrite(&value, sizeof(double), 1, fp);
              }
            }
          }
        } else {
          long nc = Th.nv;
          unsigned nbytes = nc * sizeof(double) * (*this).nbfloat;
          int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);

          for (int iv = 0; iv < Th.nv; iv++) {
            for (int j = 0; j < (*this).nbfloat; j++) {
              valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
              double value = valsol[iv * (*this).nbfloat + j];

              l = runEncodeB64(sizeof(double), (unsigned char *)&value, ElementChars);
              ElementChars[l] = 0;
              fwrite(&ElementChars, l, 1, fp);
            }
          }

          // flush buffer
          l = runEncodeB64(0, NULL, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      } else {
        for (int iv = 0; iv < Th.nv; iv++) {
          for (int j = 0; j < (*this).nbfloat; j++) {
            valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
            double value = valsol[iv * (*this).nbfloat + j];

            fprintf(fp, "%.16e ", value);
          }
        }
      }

      fprintf(fp, "\n");
    }
  };
  vector< Expression2 > l;
#ifndef COMMON_HPDDM_PARALLEL_IO
  static const int n_name_param = 8;
#else
  static const int n_name_param = 9;
#endif
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

 public:
  VTK_WriteMesh3_Op(const basicAC_F0 &args) : l(args.size( ) - 2) {
    int nbofsol;
    int ddim = 3;
    int stsize = 6;
    int sca = 0, vec = 0, ten = 0;
    string scas("scalaire");
    string vecs("vector");
    string tens("tensor");

    if (verbosity > 2) {
      cout << "Write Mesh and Solutions in VTK Formats" << endl;
    }

    args.SetNameParam(n_name_param, name_param, nargs);

    if (BCastTo< string * >(args[0])) {
      filename = CastTo< string * >(args[0]);
    }

    if (BCastTo< pmesh3 >(args[1])) {
      eTh = CastTo< pmesh3 >(args[1]);
    }

    nbofsol = l.size( );

    for (size_t i = 2; i < args.size( ); i++) {
      size_t jj = i - 2;

      if (BCastTo< double >(args[i])) {
        l[jj].what = 1;
        l[jj].nbfloat = 1;
        l[jj][0] = to< double >(args[i]);

        char number[16];
        sprintf(number, "%li", jj + 1);
        l[jj].name = scas;
        l[jj].name += number;
        sca++;
      } else if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a0 = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        // cout << "taille" << a0->size() << endl;
        if (a0->size( ) != ddim && a0->size( ) != stsize) {
          CompileError(
            "savesol in 3D: vector solution is 3 composant, tensor solution is 6 composant");
        }

        if (a0->size( ) == ddim) {
          // vector solution
          l[jj].what = 2;
          l[jj].nbfloat = ddim;

          for (int j = 0; j < ddim; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }

          char number[16];
          sprintf(number, "%li", jj + 1);
          l[jj].name = vecs;
          l[jj].name += number;
          vec++;
        } else if (a0->size( ) == stsize) {
          // symmetric tensor solution
          l[jj].what = 3;
          l[jj].nbfloat = stsize;

          for (int j = 0; j < stsize; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }

          char number[16];
          sprintf(number, "%li", jj + 1);
          l[jj].name = tens;
          l[jj].name += number;
          ten++;
        }
      } else {
        cout << " arg " << i << " " << args[i].left( ) << endl;
        CompileError("savesol in 2D: Sorry no way to save this kind of data");
      }
    }
  }

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< string * >( ), atype< pmesh3 >( ), true);
  }    // all type

  static E_F0 *f(const basicAC_F0 &args) { return new VTK_WriteMesh3_Op(args); }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type VTK_WriteMesh3_Op::name_param[] = {{"dataname", &typeid(string *)},
                                                             {"withsurfacemesh", &typeid(bool)},
                                                             {"order", &typeid(KN_< long >)},
                                                             // A rajouter dans le 3D
                                                             {"floatmesh", &typeid(bool)},
                                                             {"floatsol", &typeid(bool)},
                                                             {"bin", &typeid(bool)},
                                                             {"swap", &typeid(bool)},
#ifdef COMMON_HPDDM_PARALLEL_IO
                                                             {"communicator", &typeid(pcommworld)},
#endif
                                                             {"append", &typeid(bool)}};

void VTK_WRITE_MESH3(const string &filename, FILE *fp, const Mesh3 &Th, bool binary, int datasize,
                     bool surface, bool bigEndian) {
  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "%s, Created by Freefem++ \n", filename.c_str( ));
  if (binary) {
    fprintf(fp, "BINARY\n");
  } else {
    fprintf(fp, "ASCII\n");
  }

  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
  // get all the entities in the model

  // write mesh vertices
  if (datasize == sizeof(float)) {
    fprintf(fp, "POINTS %d float\n", Th.nv);
  }
  else if (datasize == sizeof(double)) {
    fprintf(fp, "POINTS %d double\n", Th.nv);
  }
  else
      ffassert(0); 
  if (datasize == sizeof(float)) {
    for (unsigned int i = 0; i < Th.nv; i++) {
      const Vertex3 &P = Th.vertices[i];
      float f[3];
      f[0] = P.x;
      f[1] = P.y;
      f[2] = P.z;
      if (binary) {
        if (!bigEndian) {
          FreeFEM::SwapBytes((char *)&f, sizeof(float), 3);
        }

        fwrite(&f, sizeof(float), 3, fp);
      } else {
        fprintf(fp, "%.8f %.8f %.8f\n", f[0], f[1], f[2]);
      }
    }
  } else if (datasize == sizeof(double)) {
    for (unsigned int i = 0; i < Th.nv; i++) {
      const Vertex3 &P = Th.vertices[i];
      double f[3];
      f[0] = P.x;
      f[1] = P.y;
      f[2] = P.z;    // 3D case
       if (binary) {
        if (!bigEndian) {
          FreeFEM::SwapBytes((char *)&f, sizeof(double), 3);
        }
        fwrite(&f, sizeof(double), 3, fp);
      } else {
        fprintf(fp, "%lf %lf %lf\n", f[0], f[1], f[2]);
      }
    }
  }

  fprintf(fp, "\n");
  if (verbosity > 1) {
    printf("writing vertices is finish %ld \n",ftell(fp));
  }

  if (verbosity > 1) {
    printf("writing elements now\n");
  }

  //= ===============
  // CELL
  //= ===============
  // loop over all elements we need to save and count vertices
  int numElements, totalNumInt;
  if (surface) {
    numElements = Th.nt + Th.nbe;
    totalNumInt = Th.nt * 4 + Th.nbe * 3 + numElements;
  } else {
    numElements = Th.nt;
    totalNumInt = Th.nt * 4 + numElements;
  }

  if (verbosity > 1) {
    printf("writing cells \n");
  }

  // print vertex indices in ascii or binary
  fprintf(fp, "CELLS %d %d\n", numElements, totalNumInt);

  if (binary) {
    int IntType = 4;
    if (verbosity > 1) {
      printf("writing tetrahedron elements \n");
    }

    for (int it = 0; it < Th.nt; it++) {
      const Tet &K(Th.elements[it]);
      int iv[IntType + 1];

      iv[0] = IntType;

      for (int ii = 0; ii < IntType; ii++) {
        iv[ii + 1] = Th.operator( )(K[ii]);
      }

      if (!bigEndian) {
        FreeFEM::SwapBytes((char *)&iv, sizeof(int), IntType + 1);
      }

      fwrite(&iv, sizeof(int), IntType + 1, fp);
    }

    if (surface) {
      if (verbosity > 1) {
        printf("writing triangle elements \n");
      }

      IntType = 3;

      for (int ibe = 0; ibe < Th.nbe; ibe++) {
        const Triangle3 &K(Th.be(ibe));
        int iv[IntType + 1];
        iv[0] = IntType;

        for (int ii = 0; ii < IntType; ii++) {
          iv[ii + 1] = Th.operator( )(K[ii]);
        }

        if (!bigEndian) {
          FreeFEM::SwapBytes((char *)&iv, sizeof(int), IntType + 1);
        }

        fwrite(&iv, sizeof(int), IntType + 1, fp);
      }
    }
  } else {
    int IntType = 4;
    if (verbosity > 1) {
      printf("writing tetrahedron elements \n");
    }

    for (int it = 0; it < Th.nt; it++) {
      const Tet &K(Th.elements[it]);
      int iv[IntType + 1];
      iv[0] = IntType;

      for (int ii = 0; ii < IntType; ii++) {
        iv[ii + 1] = Th.operator( )(K[ii]);
      }

      fprintf(fp, "%d %d %d %d %d\n", iv[0], iv[1], iv[2], iv[3], iv[4]);
    }

    if (surface) {
      if (verbosity > 1) {
        printf("writing triangle elements \n");
      }

      IntType = 3;

      for (int ibe = 0; ibe < Th.nbe; ibe++) {
        const Triangle3 &K(Th.be(ibe));
        int iv[IntType + 1];
        iv[0] = IntType;

        for (int ii = 0; ii < IntType; ii++) {
          iv[ii + 1] = Th.operator( )(K[ii]);
        }

        fprintf(fp, "%d %d %d %d\n", iv[0], iv[1], iv[2], iv[3]);
      }
    }
  }

  fprintf(fp, "\n");

  // CELL_TYPE
  // print element types in ascii or binary

  fprintf(fp, "CELL_TYPES %d\n", numElements);
  if (binary) {
    unsigned int type;

    for (int it = 0; it < Th.nt; it++) {
      type = VTK_TET;
      if (!bigEndian) {
        FreeFEM::SwapBytes((char *)&type, sizeof(int), 1);
      }

      fwrite(&type, sizeof(int), 1, fp);
    }

    if (surface) {
      for (int ibe = 0; ibe < Th.nbe; ibe++) {
        type = VTK_TRI;
        if (!bigEndian) {
          FreeFEM::SwapBytes((char *)&type, sizeof(int), 1);
        }

        fwrite(&type, sizeof(int), 1, fp);
      }
    }
  } else {
    unsigned int type;

    for (int it = 0; it < Th.nt; it++) {
      type = VTK_TET;
      fprintf(fp, "%d ", type);
    }

    if (surface) {
      for (int ibe = 0; ibe < Th.nbe; ibe++) {
        type = VTK_TRI;
        fprintf(fp, "%d ", type);
      }
    }
  }

  fprintf(fp, "\n");

  //= ================================
  // WRITE SOLUTION IN FORMAT VTK
  // LABEL OF ELEMENTS
  //= ================================

  list< int > list_label_Elem;
  // list<int> list_label_Border_Elem;
  {
    for (int it = 0; it < Th.nt; it++) {
      const Tet &K(Th.elements[it]);
      list< int >::const_iterator ilist;
      int labOk = 0;

      for (ilist = list_label_Elem.begin( ); ilist != list_label_Elem.end( ); ilist++) {
        if (*ilist == K.lab) {
          labOk = 1;
          break;
        }
      }

      if (labOk == 0) {
        list_label_Elem.push_back(K.lab);
      }
    }

    if (surface) {
      for (int ibe = 0; ibe < Th.nbe; ibe++) {
        const Triangle3 &K(Th.be(ibe));
        list< int >::const_iterator ilist;
        int labOk = 0;

        for (ilist = list_label_Elem.begin( ); ilist != list_label_Elem.end( ); ilist++) {
          if (*ilist == K.lab) {
            labOk = 1;
            break;
          }
        }

        if (labOk == 0) {
          list_label_Elem.push_back(K.lab);
        }
      }
    }
  }
  list_label_Elem.sort( );

  //= ================================
  //= ================================

  fprintf(fp, "CELL_DATA %d\n", numElements);
  int cell_fd = 1;
  fprintf(fp, "Scalars  Label int %d\n", cell_fd);
  fprintf(fp, "LOOKUP_TABLE FreeFempp_table\n");
  // Determination des labels
  if (binary) {
    int label;

    for (int it = 0; it < Th.nt; it++) {
      const Tet &K(Th.elements[it]);
      label = K.lab;
      if (!bigEndian) {
        FreeFEM::SwapBytes((char *)&label, sizeof(int), 1);
      }

      fwrite(&label, sizeof(int), 1, fp);
    }

    if (surface) {
      for (int ibe = 0; ibe < Th.nbe; ibe++) {
        const Triangle3 &K(Th.be(ibe));
        label = K.lab;
        if (!bigEndian) {
          FreeFEM::SwapBytes((char *)&label, sizeof(int), 1);
        }

        fwrite(&label, sizeof(int), 1, fp);
      }
    }
  } else {
    int label;

    for (int it = 0; it < Th.nt; it++) {
      const Tet &K(Th.elements[it]);
      label = K.lab;
      fprintf(fp, "%d\n", label);
    }

    if (surface) {
      for (int ibe = 0; ibe < Th.nbe; ibe++) {
        const Triangle3 &K(Th.be(ibe));
        label = K.lab;
        fprintf(fp, "%d\n", label);
      }
    }
  }

  fprintf(fp, "\n");

  int size_list = 0;
  list< int >::const_iterator ilist;

  for (ilist = list_label_Elem.begin( ); ilist != list_label_Elem.end( ); ilist++) {
    size_list++;
  }

  fprintf(fp, "LOOKUP_TABLE FreeFempp_table %d\n", size_list);
  {
    list< int >::const_iterator ilist;

    for (ilist = list_label_Elem.begin( ); ilist != list_label_Elem.end( ); ilist++) {
      if (binary) {
        int tab[4];
        tab[0] = (int)ColorTable[abs(*ilist) % NbColorTable][0] * 255;
        tab[1] = (int)ColorTable[abs(*ilist) % NbColorTable][1] * 255;
        tab[2] = (int)ColorTable[abs(*ilist) % NbColorTable][2] * 255;
        tab[3] = 255;

        for (int itab = 0; itab < 4; itab++) {
          char newvalue[sizeof(int)];
          int bid0 = sprintf(newvalue, "%s", (char *)&tab[itab]);
          fwrite(&newvalue, sizeof(unsigned char), 1, fp);
        }
      } else {
        float tab[4];
        tab[0] = ColorTable[abs(*ilist) % NbColorTable][0];
        tab[1] = ColorTable[abs(*ilist) % NbColorTable][1];
        tab[2] = ColorTable[abs(*ilist) % NbColorTable][2];
        tab[3] = 1.0;

        fprintf(fp, "%.8f %.8f %.8f %.8f\n", tab[0], tab[1], tab[2], tab[3]);
      }
    }
  }
  fprintf(fp, "\n");
}

AnyType VTK_WriteMesh3_Op::operator( )(Stack stack) const {
  string *pffname = GetAny< string * >((*filename)(stack));
  Mesh3 *pTh = GetAny< Mesh3 * >((*eTh)(stack));

  ffassert(pTh);
  Mesh3 &Th = *pTh;
  bool bigEndian = isBigEndian( );
  bool swap = bigEndian;
  bool binary = true;
  bool surface = true;
  bool floatmesh = false;
  bool floatsol = false;

  string *dataname;
  int nbofsol = l.size( );
  KN< int > order(nbofsol);

  char *nameofuser[nbofsol];

  for (int ii = 0; ii < nbofsol; ii++) {
    order[ii] = 0;
  }

  if (nargs[0]) {
    dataname = GetAny< string * >((*nargs[0])(stack));
  }

  if (nargs[1]) {
    surface = GetAny< bool >((*nargs[1])(stack));
  }

  if (nargs[2]) {
    order = GetAny< KN_< long > >((*nargs[2])(stack));
  }

  if (nargs[3]) {
    floatmesh = GetAny< bool >((*nargs[3])(stack));
  }

  if (nargs[4]) {
    floatsol = GetAny< bool >((*nargs[4])(stack));
  }

  if (nargs[5]) {
    binary = GetAny< bool >((*nargs[5])(stack));
  }

  if (nargs[6]) {
    swap = GetAny< bool >((*nargs[6])(stack));
  }
#ifdef COMMON_HPDDM_PARALLEL_IO
  parallelIO(pffname, nargs[7] ? (MPI_Comm *)GetAny< pcommworld >((*nargs[7])(stack)) : 0,
             nargs[8] && GetAny< bool >((*nargs[8])(stack)));
#endif

  int datasize = floatmesh ? sizeof(float) : sizeof(double);

  int datasizeSol = floatsol ? sizeof(float) : sizeof(double);

  int iii = 0;
  if (nargs[0]) {
    char *data = newcopy(dataname);
    if (verbosity > 5) {
      cout << "   iovtk writeMesh3: names  \"" << data << "\"" << endl;
    }

    char *name = strtok(data, " \n\0\t");
    nameofuser[iii] = newcopy(name);
    if (verbosity > 5) {
      cout << "   iovtk writeMesh3:value of iii=" << iii << " " << nameofuser[iii] << endl;
    }

    iii++;
    {
      while ((name = strtok(NULL, " \n\0\t"))) {
        if (iii >= nbofsol) {
          if (verbosity) {
            cout << "   iovtk writeMesh3: The number of data name is too large " << endl;
          }

          break;
        }

        nameofuser[iii] = newcopy(name);
        if (verbosity > 5) {
          cout << "   iovtk writeMesh3:value of iii=" << iii << " " << nameofuser[iii] << endl;
        }

        iii++;
      }

      if (iii < nbofsol) {
        if (verbosity) {
          cout << "   iovtk writeMesh3: The number of data name is too small, we give default name "
               << endl;
        }
      }
    }
  }

  if (iii < nbofsol) {
    for (int iiii = iii; iiii < nbofsol; iiii++) {
      nameofuser[iiii] = newcopy(l[iiii].name.c_str( ));
    }
  }

  // determination of number of order 0 et 1.
  int Norder0 = 0;

  for (int ii = 0; ii < nbofsol; ii++) {
    if (order[ii] == 0) {
      Norder0++;
    }
  }

  // lecture du nom des variables

  FILE *fp = fopen((*pffname).c_str( ), "wb");
  if (!fp) {
    cerr << "Unable to open file " << (*pffname).c_str( ) << endl;
    ExecError("error in reading vtk file");
  }

  // type of VTK FILE
  char *pch = newcopy(pffname);
  int VTK_FILE = 0;
  int ls = 0;
  int lll = strlen(pch);
  if (!strcmp(pch + lll - (ls = 4), ".vtk")) {
    VTK_FILE = 1;
  } else if (!strcmp(pch + lll - (ls = 4), ".vtu")) {
    VTK_FILE = 2;
  }

  if (VTK_FILE == 1) {
    VTK_WRITE_MESH3(*pffname, fp, Th, binary, datasize, surface, swap);

    if (datasizeSol == sizeof(float)) {
      if (Norder0 > 0) {
        fprintf(fp, "FIELD FieldData %d\n", Norder0);

        for (int ii = 0; ii < nbofsol; ii++) {
          if (order[ii] == 0) {
            int nsol;
            if (surface) {
              nsol = Th.nt + Th.nbe;
            } else {
              nsol = Th.nt;
            }

            fprintf(fp, "%s %ld %d float\n", nameofuser[ii], l[ii].nbfloat, nsol);
            if (verbosity > 5) {
              cout << "   iovtk writeMesh3: name of data(" << ii << ")=" << nameofuser[ii] << endl;
            }

            // changement ecriture solution
            l[ii].writesolutionP0_float(fp, Th, stack, surface, binary, swap);
            // fprintf(fp,"\n");
          }
        }
      }

      if (Norder0 < nbofsol) {
        fprintf(fp, "POINT_DATA %d\n", Th.nv);
        fprintf(fp, "FIELD FieldData %d\n", nbofsol - Norder0);

        for (int ii = 0; ii < nbofsol; ii++) {
          if (order[ii] == 1) {
            fprintf(fp, "%s %ld %d float\n", nameofuser[ii], l[ii].nbfloat, Th.nv);
            if (verbosity > 5) {
              cout << "   iovtk writeMesh3:name of data(" << ii << ")=" << nameofuser[ii] << endl;
            }

            // changement ecriture solution
            l[ii].writesolutionP1_float(fp, Th, stack, binary, swap);
          }
        }
      }
    }

    if (datasizeSol == sizeof(double)) {
      if (Norder0 > 0) {
        fprintf(fp, "FIELD FieldData %d\n", Norder0);

        for (int ii = 0; ii < nbofsol; ii++) {
          if (order[ii] == 0) {
            int nsol;
            if (surface) {
              nsol = Th.nt + Th.nbe;
            } else {
              nsol = Th.nt;
            }

            fprintf(fp, "%s %ld %d double\n", nameofuser[ii], l[ii].nbfloat, nsol);
            if (verbosity > 5) {
              cout << "   iovtk writeMesh3:name of data(" << ii << ")=" << nameofuser[ii] << endl;
            }

            // changement ecriture solution
            l[ii].writesolutionP0_double(fp, Th, stack, surface, binary, swap);
          }
        }
      }

      if (Norder0 < nbofsol) {
        fprintf(fp, "POINT_DATA %d\n", Th.nv);
        fprintf(fp, "FIELD FieldData %d\n", nbofsol - Norder0);

        for (int ii = 0; ii < nbofsol; ii++) {
          if (order[ii] == 1) {
            fprintf(fp, "%s %ld %d double\n", nameofuser[ii], l[ii].nbfloat, Th.nv);
            if (verbosity > 5) {
              cout << "   iovtk writeMesh3:name of data(" << ii << ")=" << nameofuser[ii] << endl;
            }

            // changement ecriture solution
            l[ii].writesolutionP1_double(fp, Th, stack, binary, swap);
          }
        }
      }
    }
  } else if (VTK_FILE == 2) {
    VTU_WRITE_MESH(fp, Th, binary, datasize, surface);
    // Solution Order
    // order 0
    if (Norder0 > 0) {
      for (int ii = 0; ii < nbofsol; ii++) {
        if (order[ii] == 0) {
          if (datasize == sizeof(float)) {
            VTU_DATA_ARRAY(fp, "Float32", nameofuser[ii], l[ii].nbfloat, binary);
            l[ii].writesolutionP0_float(fp, Th, stack, surface, binary, swap, 1);
          } else if (datasize == sizeof(double)) {
            VTU_DATA_ARRAY(fp, "Float64", nameofuser[ii], l[ii].nbfloat, binary);
            l[ii].writesolutionP0_double(fp, Th, stack, surface, binary, swap, 1);
          }

          ENDTYPE_VTU(fp, "DataArray");
        }
      }
    }
    ENDTYPE_VTU(fp, "CellData");
    // order 1
    if (Norder0 != nbofsol) {
      BEGINTYPE_VTU(fp, "PointData");

      for (int ii = 0; ii < nbofsol; ii++) {
        if (order[ii] == 1) {
          if (datasize == sizeof(float)) {
            VTU_DATA_ARRAY(fp, "Float32", nameofuser[ii], l[ii].nbfloat, binary);
            l[ii].writesolutionP1_float(fp, Th, stack, binary, swap, 1);
          } else if (datasize == sizeof(double)) {
            VTU_DATA_ARRAY(fp, "Float64", nameofuser[ii], l[ii].nbfloat, binary);
            l[ii].writesolutionP1_double(fp, Th, stack, binary, swap, 1);
          }

          ENDTYPE_VTU(fp, "DataArray");
        }
      }
      ENDTYPE_VTU(fp, "PointData");
    }

    ENDTYPE_VTU(fp, "Piece");
    ENDTYPE_VTU(fp, "UnstructuredGrid");
    ENDTYPE_VTU(fp, "VTKFile");
  } else {
    cout << "extension file of VTK is not correct" << endl;
    exit(1);
  }

  fclose(fp);

  for (int iiii = 0; iiii < nbofsol; iiii++) {
    delete[] nameofuser[iiii];
  }

  delete[] pch;
  return (Mesh3 *)NULL;
}

//= =============================================
// FIN ECRITURE DE FICHIER .vtk (3D VOLUME)
//= =============================================

//= ======================
// FIN 3D VOLUME Fichier .vtk
//= ======================

//= =============================================
// ECRITURE DE FICHIER .vtk
//= =============================================

template< class MMesh >
class VTK_WriteMeshT_Op : public E_F0mps {
 public:
  typedef long Result;
  Expression eTh;
  typedef const MMesh *ppmesh;
  /*typedef typename MMesh::Element T;
   typedef typename MMesh::Element::RdHat TRdHat;

  typedef typename MMesh::BorderElement B;
  typedef typename MMesh::Vertex V;

  typedef typename MMesh::BorderElement::RdHat BRdHat;*/

  Expression filename;
  struct Expression2 {
    string name;
    long what;       // 1 scalar, 2 vector, 3 symtensor
    long nbfloat;    // 1 scalar, 2 vector (3D), 3 symtensor(3D)
    Expression e[3];
    Expression2( ) {
      e[0] = 0;
      e[1] = 0;
      e[2] = 0;
      what = 0;
      nbfloat = 0;
    };
    Expression &operator[](int i) { return e[i]; }

    double eval(int i, Stack stack) const {
      if (e[i])
        return GetAny< double >((*e[i])(stack));
      else
        return 0;
    }

    void writesolutionP0_float_binary(FILE *fp, const MMesh &Th, Stack stack, bool surface,
                                      bool bigEndian) const {
      MeshPoint *mp3(MeshPointStack(stack));
      typedef typename MMesh::Element T;
      typedef typename MMesh::Element::RdHat TRdHat;
      double k = T::nv;
      TRdHat Cdg_hat = TRdHat::diag(1. / k);

      if (!bigEndian) {
        for (int it = 0; it < Th.nt; it++) {
          const T &K(Th.elements[it]);
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < nbfloat; j++) {
            float value = eval(j, stack);
            FreeFEM::SwapBytes((char *)&value, sizeof(float), 1);
            fwrite(&value, sizeof(float), 1, fp);
          }
        }

        if (surface) {
          for (int ibe = 0; ibe < Th.nbe; ibe++) {
            // determine the element contening this border
            int ie;
            int it = Th.BoundaryElement(ibe, ie);
            const T &K(Th.elements[it]);
            mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

            for (int j = 0; j < nbfloat; j++) {
              float value = eval(j, stack);

              FreeFEM::SwapBytes((char *)&value, sizeof(float), 1);
              fwrite(&value, sizeof(float), 1, fp);
            }
          }
        }
      } else {
        for (int it = 0; it < Th.nt; it++) {
          const T &K(Th.elements[it]);
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < nbfloat; j++) {
            float value = eval(j, stack);

            fwrite(&value, sizeof(float), 1, fp);
          }
        }

        if (surface) {
          for (int ibe = 0; ibe < Th.nbe; ibe++) {
            // determine the element contening this border
            int ie;
            int it = Th.BoundaryElement(ibe, ie);
            const T &K(Th.elements[it]);
            mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

            for (int j = 0; j < nbfloat; j++) {
              float value = eval(j, stack);

              fwrite(&value, sizeof(float), 1, fp);
            }
          }
        }
      }

      fprintf(fp, "\n");
    }

    void writesolutionP0_float(FILE *fp, const MMesh &Th, Stack stack, bool surface) const {
      MeshPoint *mp3(MeshPointStack(stack));
      typedef typename MMesh::Element T;
      typedef typename MMesh::Element::RdHat TRdHat;
      double k = T::nv;
      TRdHat Cdg_hat = TRdHat::diag(1. / k);

      for (int it = 0; it < Th.nt; it++) {
        const T &K(Th.t(it));
        mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

        for (int j = 0; j < nbfloat; j++) {
          float value = eval(j, stack);

          fprintf(fp, "%.8e ", value);
        }
      }

      if (surface) {
        for (int ibe = 0; ibe < Th.nbe; ibe++) {
          // determination du triangle contenant cette edge
          int ie;
          int it = Th.BoundaryElement(ibe, ie);
          const T &K(Th.t(it));
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < nbfloat; j++) {
            float value = eval(j, stack);

            fprintf(fp, "%.8e ", value);
          }
        }
      }

      fprintf(fp, "\n");
    }

    void writesolutionP0_float_XML(FILE *fp, const MMesh &Th, Stack stack, bool surface) const {
      MeshPoint *mp3(MeshPointStack(stack));
      typedef typename MMesh::Element T;
      typedef typename MMesh::Element::RdHat TRdHat;
      double k = T::nv;
      TRdHat Cdg_hat = TRdHat::diag(1. / k);
      unsigned char ElementChars[256];
      long nc = Th.nt;

      if (surface) nc = nc + Th.nbe;

      unsigned nbytes = nc * sizeof(float) * (*this).nbfloat;
      int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      for (int it = 0; it < Th.nt; it++) {
        const T &K(Th.t(it));
        mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

        for (int j = 0; j < nbfloat; j++) {
          float value = eval(j, stack);

          l = runEncodeB64(sizeof(float), (unsigned char *)&value, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      }

      if (surface) {
        for (int ibe = 0; ibe < Th.nbe; ibe++) {
          // determination du triangle contenant cette edge
          int ie;
          int it = Th.BoundaryElement(ibe, ie);
          const T &K(Th.t(it));
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < nbfloat; j++) {
            float value = eval(j, stack);

            l = runEncodeB64(sizeof(float), (unsigned char *)&value, ElementChars);
            ElementChars[l] = 0;
            fwrite(&ElementChars, l, 1, fp);
          }
        }
      }

      // flush buffer
      l = runEncodeB64(0, NULL, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);
      fprintf(fp, "\n");
    }

    void writesolutionP0_float(FILE *fp, const MMesh &Th, Stack stack, bool surface, bool binary,
                               bool bigEndian, bool XML = false) const {
      if (binary) {
        if (!XML) {
          (*this).writesolutionP0_float_binary(fp, Th, stack, surface, bigEndian);
        } else {
          (*this).writesolutionP0_float_XML(fp, Th, stack, surface);
        }
      } else {
        (*this).writesolutionP0_float(fp, Th, stack, surface);
      }
    }

    void writesolutionP0_double_binary(FILE *fp, const MMesh &Th, Stack stack, bool surface,
                                       bool bigEndian) const {
      MeshPoint *mp3(MeshPointStack(stack));
      typedef typename MMesh::Element T;
      typedef typename MMesh::Element::RdHat TRdHat;
      double k = T::nv;
      TRdHat Cdg_hat = TRdHat::diag(1. / k);

      if (!bigEndian) {
        for (int it = 0; it < Th.nt; it++) {
          const T &K(Th.elements[it]);
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < (*this).nbfloat; j++) {
            double value = (*this).eval(j, stack);

            FreeFEM::SwapBytes((char *)&value, sizeof(double), 1);
            fwrite(&value, sizeof(double), 1, fp);
          }
        }

        if (surface) {
          for (int ibe = 0; ibe < Th.nbe; ibe++) {
            // determination du triangle contenant cette edge
            int ie;
            int it = Th.BoundaryElement(ibe, ie);
            const T &K(Th.elements[it]);
            mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

            for (int j = 0; j < (*this).nbfloat; j++) {
              double value = (*this).eval(j, stack);

              FreeFEM::SwapBytes((char *)&value, sizeof(double), 1);
              fwrite(&value, sizeof(double), 1, fp);
            }
          }
        }
      } else {
        for (int it = 0; it < Th.nt; it++) {
          const T &K(Th.elements[it]);
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < (*this).nbfloat; j++) {
            double value = (*this).eval(j, stack);

            fwrite(&value, sizeof(double), 1, fp);
          }
        }

        if (surface) {
          for (int ibe = 0; ibe < Th.nbe; ibe++) {
            // determination du triangle contenant cette edge
            int ie;
            int it = Th.BoundaryElement(ibe, ie);
            const T &K(Th.elements[it]);
            mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

            for (int j = 0; j < (*this).nbfloat; j++) {
              double value = (*this).eval(j, stack);

              fwrite(&value, sizeof(double), 1, fp);
            }
          }
        }
      }

      fprintf(fp, "\n");
    }

    void writesolutionP0_double(FILE *fp, const MMesh &Th, Stack stack, bool surface) const {
      MeshPoint *mp3(MeshPointStack(stack));
      typedef typename MMesh::Element T;
      typedef typename MMesh::Element::RdHat TRdHat;
      double k = T::nv;
      TRdHat Cdg_hat = TRdHat::diag(1. / k);

      for (int it = 0; it < Th.nt; it++) {
        const T &K(Th.t(it));
        mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

        for (int j = 0; j < (*this).nbfloat; j++) {
          double value = (*this).eval(j, stack);

          fprintf(fp, "%.16e ", value);
        }
      }

      if (surface) {
        for (int ibe = 0; ibe < Th.nbe; ibe++) {
          // determination du triangle contenant cette edge
          int ie;
          int it = Th.BoundaryElement(ibe, ie);
          const T &K(Th.t(it));
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < (*this).nbfloat; j++) {
            double value = (*this).eval(j, stack);
            fprintf(fp, "%.16e ", value);
          }
        }
      }

      fprintf(fp, "\n");
    }

    void writesolutionP0_double_XML(FILE *fp, const MMesh &Th, Stack stack, bool surface) const {
      unsigned char ElementChars[256];
      MeshPoint *mp3(MeshPointStack(stack));
      typedef typename MMesh::Element T;
      typedef typename MMesh::Element::RdHat TRdHat;
      double k = T::nv;
      TRdHat Cdg_hat = TRdHat::diag(1. / k);
      long nc = Th.nt;

      if (surface) nc = nc + Th.nbe;

      unsigned nbytes = nc * sizeof(double) * (*this).nbfloat;
      int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);

      for (int it = 0; it < Th.nt; it++) {
        const T &K(Th.t(it));
        mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

        for (int j = 0; j < (*this).nbfloat; j++) {
          double value = (*this).eval(j, stack);

          l = runEncodeB64(sizeof(double), (unsigned char *)&value, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      }

      if (surface) {
        for (int ibe = 0; ibe < Th.nbe; ibe++) {
          // determination du triangle contenant cette edge
          int ie;
          int it = Th.BoundaryElement(ibe, ie);
          const T &K(Th.t(it));
          mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

          for (int j = 0; j < (*this).nbfloat; j++) {
            double value = (*this).eval(j, stack);

            l = runEncodeB64(sizeof(double), (unsigned char *)&value, ElementChars);
            ElementChars[l] = 0;
            fwrite(&ElementChars, l, 1, fp);
          }
        }
      }

      // flush buffer
      l = runEncodeB64(0, NULL, ElementChars);
      ElementChars[l] = 0;
      fwrite(&ElementChars, l, 1, fp);
      fprintf(fp, "\n");
    }

    void writesolutionP0_double(FILE *fp, const MMesh &Th, Stack stack, bool surface, bool binary,
                                bool bigEndian, bool XML = false) const {
      if (binary)
        if (!XML)
          (*this).writesolutionP0_double_binary(fp, Th, stack, surface, bigEndian);
        else
          (*this).writesolutionP0_double_XML(fp, Th, stack, surface);

      else
        (*this).writesolutionP0_double(fp, Th, stack, surface);
    }

    void writesolutionP1_float(FILE *fp, const MMesh &Th, Stack stack, bool binary, bool bigEndian,
                               bool XML = false) const {
      unsigned char ElementChars[256];
      MeshPoint *mp3(MeshPointStack(stack));

      KN< double > valsol(Th.nv * (*this).nbfloat);
      KN< int > takemesh(Th.nv);
      takemesh = 0;
      valsol = 0.;

      for (int it = 0; it < Th.nt; it++)
        for (int iv = 0; iv < 3; iv++) {
          int i = Th(it, iv);
          mp3->setP(&Th, it, iv);

          for (int j = 0; j < (*this).nbfloat; j++)
            valsol[i * (*this).nbfloat + j] =
              valsol[i * (*this).nbfloat + j] + (*this).eval(j, stack);
          takemesh[i] = takemesh[i] + 1;
        }

      if (binary) {
        if (!XML) {
          if (!bigEndian) {
            for (int iv = 0; iv < Th.nv; iv++) {
              for (int j = 0; j < (*this).nbfloat; j++) {
                valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
                float value = valsol[iv * (*this).nbfloat + j];
                FreeFEM::SwapBytes((char *)&value, sizeof(float), 1);
                fwrite(&value, sizeof(float), 1, fp);
              }
            }
          } else {
            for (int iv = 0; iv < Th.nv; iv++) {
              for (int j = 0; j < (*this).nbfloat; j++) {
                valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
                float value = valsol[iv * (*this).nbfloat + j];

                fwrite(&value, sizeof(float), 1, fp);
              }
            }
          }
        } else {
          long nc = Th.nv;
          unsigned nbytes = nc * sizeof(float) * (*this).nbfloat;
          int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);

          for (int iv = 0; iv < Th.nv; iv++) {
            for (int j = 0; j < (*this).nbfloat; j++) {
              valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
              float value = valsol[iv * (*this).nbfloat + j];

              l = runEncodeB64(sizeof(float), (unsigned char *)&value, ElementChars);
              ElementChars[l] = 0;
              fwrite(&ElementChars, l, 1, fp);
            }
          }

          // flush buffer
          l = runEncodeB64(0, NULL, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      } else {
        for (int iv = 0; iv < Th.nv; iv++)
          for (int j = 0; j < (*this).nbfloat; j++) {
            valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
            float value = valsol[iv * (*this).nbfloat + j];

            fprintf(fp, "%.8e ", value);
          }
      }
      fprintf(fp, "\n");
    }

    void writesolutionP1_double(FILE *fp, const MMesh &Th, Stack stack, bool binary, bool bigEndian,
                                bool XML = false) const {
      unsigned char ElementChars[256];
      MeshPoint *mp3(MeshPointStack(stack));

      KN< double > valsol(Th.nv * (*this).nbfloat);
      KN< int > takemesh(Th.nv);
      takemesh = 0;
      valsol = 0.;

      for (int it = 0; it < Th.nt; it++)
        for (int iv = 0; iv < 3; iv++) {
          int i = Th(it, iv);
          mp3->setP(&Th, it, iv);

          for (int j = 0; j < (*this).nbfloat; j++)
            valsol[i * (*this).nbfloat + j] =
              valsol[i * (*this).nbfloat + j] + (*this).eval(j, stack);
          takemesh[i] = takemesh[i] + 1;
        }

      if (binary) {
        if (!XML) {
          if (!bigEndian)
            for (int iv = 0; iv < Th.nv; iv++)
              for (int j = 0; j < (*this).nbfloat; j++) {
                valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
                double value = valsol[iv * (*this).nbfloat + j];

                FreeFEM::SwapBytes((char *)&value, sizeof(double), 1);
                fwrite(&value, sizeof(double), 1, fp);
              }

          else
            for (int iv = 0; iv < Th.nv; iv++)
              for (int j = 0; j < (*this).nbfloat; j++) {
                valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
                double value = valsol[iv * (*this).nbfloat + j];

                fwrite(&value, sizeof(double), 1, fp);
              }
        } else {
          long nc = Th.nv;
          unsigned nbytes = nc * sizeof(double) * (*this).nbfloat;
          int l = runEncodeB64(sizeof(int), (unsigned char *)&nbytes, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);

          for (int iv = 0; iv < Th.nv; iv++)
            for (int j = 0; j < (*this).nbfloat; j++) {
              valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
              double value = valsol[iv * (*this).nbfloat + j];

              l = runEncodeB64(sizeof(double), (unsigned char *)&value, ElementChars);
              ElementChars[l] = 0;
              fwrite(&ElementChars, l, 1, fp);
            }
          // flush buffer
          l = runEncodeB64(0, NULL, ElementChars);
          ElementChars[l] = 0;
          fwrite(&ElementChars, l, 1, fp);
        }
      } else {
        for (int iv = 0; iv < Th.nv; iv++)
          for (int j = 0; j < (*this).nbfloat; j++) {
            valsol[iv * (*this).nbfloat + j] = valsol[iv * (*this).nbfloat + j] / takemesh[iv];
            double value = valsol[iv * (*this).nbfloat + j];

            fprintf(fp, "%.16e ", value);
          }
      }

      fprintf(fp, "\n");
    }
  };

  vector< Expression2 > l;
  static const int n_name_param = 7;
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }
  int arg(int i, Stack stack, int a) const {
    return nargs[i] ? GetAny< int >((*nargs[i])(stack)) : a;
  }
  bool arg(int i, Stack stack, bool a) const {
    return nargs[i] ? GetAny< bool >((*nargs[i])(stack)) : a;
  }
  double arg(int i, Stack stack, double a) const {
    return nargs[i] ? GetAny< double >((*nargs[i])(stack)) : a;
  }

 public:
  VTK_WriteMeshT_Op(const basicAC_F0 &args) : l(args.size( ) - 2) {
    int nbofsol;
    int ddim = 3;
    int stsize = 3;
    int sca = 0, vec = 0, ten = 0;
    string scas("scalaire");
    string vecs("vector");
    string tens("tensor");

    if ((std::is_same< MMesh, MeshS >::value) && (verbosity > 2))
      cout << "Write MeshS and Solutions in VTK Formats" << endl;
    else if ((std::is_same< MMesh, MeshL >::value) && (verbosity > 2))
      cout << "Write MeshL and Solutions in VTK Formats" << endl;

    args.SetNameParam(n_name_param, name_param, nargs);

    if (BCastTo< string * >(args[0])) filename = CastTo< string * >(args[0]);

    if (BCastTo< ppmesh >(args[1])) eTh = CastTo< ppmesh >(args[1]);

    nbofsol = l.size( );

    for (size_t i = 2; i < args.size( ); i++) {
      size_t jj = i - 2;

      if (BCastTo< double >(args[i])) {
        l[jj].what = 1;
        l[jj].nbfloat = 1;
        l[jj][0] = to< double >(args[i]);

        char number[16];
        sprintf(number, "%li", jj + 1);
        l[jj].name = scas;
        l[jj].name += number;
        sca++;
      } else if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a0 = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        // cout << "taille" << a0->size() << endl;
        if (a0->size( ) != ddim && a0->size( ) != stsize) {
          CompileError(
            "savesol in 3D: vector solution is 3 composant, tensor solution is 3 composant");
        }

        if (a0->size( ) == ddim) {
          // vector solution
          l[jj].what = 2;
          l[jj].nbfloat = ddim;

          for (int j = 0; j < ddim; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }

          char number[16];
          sprintf(number, "%li", jj + 1);
          l[jj].name = vecs;
          l[jj].name += number;
          vec++;
        } else if (a0->size( ) == stsize) {
          // symmetric tensor solution
          l[jj].what = 3;
          l[jj].nbfloat = stsize;

          for (int j = 0; j < stsize; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }

          char number[16];
          sprintf(number, "%li", jj + 1);
          l[jj].name = tens;
          l[jj].name += number;
          ten++;
        }
      } else {
        cout << " arg " << i << " " << args[i].left( ) << endl;
        CompileError("savesol in 3D: Sorry no way to save this kind of data");
      }
    }
  }

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< string * >( ), atype< ppmesh >( ), true);
  }    // all type

  static E_F0 *f(const basicAC_F0 &args) { return new VTK_WriteMeshT_Op< MMesh >(args); }

  AnyType operator( )(Stack stack) const;
};

template<>
basicAC_F0::name_and_type VTK_WriteMeshT_Op< MeshS >::name_param[] = {
  {"dataname", &typeid(string *)},
  {"withsurfacemesh", &typeid(bool)},
  {"order", &typeid(KN_< long >)},
  // A rajouter dans le 3D
  {"floatmesh", &typeid(bool)},
  {"floatsol", &typeid(bool)},
  {"bin", &typeid(bool)},
  {"swap", &typeid(bool)}};

template<>
basicAC_F0::name_and_type VTK_WriteMeshT_Op< MeshL >::name_param[] = {
  {"dataname", &typeid(string *)},
  {"withsurfacemesh", &typeid(bool)},
  {"order", &typeid(KN_< long >)},
  // A rajouter dans le 3D
  {"floatmesh", &typeid(bool)},
  {"floatsol", &typeid(bool)},
  {"bin", &typeid(bool)},
  {"swap", &typeid(bool)}};

template< class MMesh >
void VTK_WRITE_MESHT(const string &filename, FILE *fp, const MMesh &Th, bool binary, int datasize,
                     bool surface, bool bigEndian) {
  typedef typename MMesh::Element T;
  typedef typename MMesh::BorderElement B;
  typedef typename MMesh::Vertex V;

  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "%s, Created by Freefem++ \n", filename.c_str( ));
  if (binary)
    fprintf(fp, "BINARY\n");
  else
    fprintf(fp, "ASCII\n");

  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
  // get all the entities in the model

  // write mesh vertices

  if (datasize == sizeof(float)) {
    fprintf(fp, "POINTS %d float\n", Th.nv);

    for (unsigned int i = 0; i < Th.nv; i++) {
      const V &P = Th.vertices[i];
      float f[3];
      f[0] = P.x;
      f[1] = P.y;
      f[2] = P.z;
      if (binary) {
        if (!bigEndian) FreeFEM::SwapBytes((char *)&f, sizeof(float), 3);
        fwrite(&f, sizeof(float), 3, fp);
      } else
        fprintf(fp, "%.8g %.8g %.8g\n", P.x, P.y, P.z);
    }
  } else if (datasize == sizeof(double)) {
    fprintf(fp, "POINTS %d double\n", Th.nv);

    for (unsigned int i = 0; i < Th.nv; i++) {
      const V &P = Th.vertices[i];
      double f[3];
      f[0] = P.x;
      f[1] = P.y;
      f[2] = P.z;
      if (binary) {
        if (!bigEndian) FreeFEM::SwapBytes((char *)&f, sizeof(double), 3);
        fwrite((unsigned char *)&f, sizeof(double), 3, fp);
      } else
        fprintf(fp, "%.15lg %.15lg %.15lg\n", f[0], f[1], f[2]);
    }
  }
  fprintf(fp, "\n");
  if (verbosity > 1) printf("writing vertices is finish, writing elements now\n");

  //= ===============
  // CELL
  //= ===============
  // loop over all elements we need to save and count vertices
  int numElements = surface ? Th.nt + Th.nbe : Th.nt;
  int totalNumInt =
    surface ? Th.nt * (T::nv) + Th.nbe * (B::nv) + numElements : Th.nt * (T::nv) + numElements;

  if (verbosity > 1) printf("writing cells \n");

  // print vertex indices in ascii or binary
  fprintf(fp, "CELLS %d %d\n", numElements, totalNumInt);
  if (binary) {
    if (verbosity > 1) {
      printf("writing elements \n");
    }
    for (int it = 0; it < Th.nt; it++) {
      const T &K(Th.t(it));
      int iv[(T::nv) + 1];
      iv[0] = (T::nv);

      for (int ii = 0; ii < (T::nv); ii++) iv[ii + 1] = Th.operator( )(K[ii]);
      if (!bigEndian) FreeFEM::SwapBytes((char *)&iv, sizeof(int), (T::nv) + 1);
      fwrite(&iv, sizeof(int), (T::nv) + 1, fp);
    }

    if (surface) {
      if (verbosity > 1) {
        printf("writing border elements \n");
      }
      for (int ibe = 0; ibe < Th.nbe; ibe++) {
        const B &K(Th.be(ibe));
        int iv[(B::nv) + 1];
        iv[0] = (B::nv);

        for (int ii = 0; ii < (B::nv); ii++) iv[ii + 1] = Th.operator( )(K[ii]);
        if (!bigEndian) FreeFEM::SwapBytes((char *)&iv, sizeof(int), (B::nv) + 1);
        fwrite(&iv, sizeof(int), (B::nv) + 1, fp);
      }
    }
  } else {
    if (verbosity > 1) printf("writing  elements \n");
    for (int it = 0; it < Th.nt; it++) {
      const T &K(Th.t(it));
      int iv[(T::nv) + 2];
      iv[0] = (T::nv);
      for (int ii = 0; ii < (T::nv); ii++) iv[ii + 1] = Th.operator( )(K[ii]);
        if (T::nv==3)
        fprintf(fp, "%d %d %d %d\n", iv[0], iv[1], iv[2], iv[3]);
      else if (T::nv==2)
        fprintf(fp, "%d %d %d\n", iv[0], iv[1], iv[2]);
    }
    if (surface) {
      if (verbosity > 1) printf("writing border elements \n");
      for (int ibe = 0; ibe < Th.nbe; ibe++) {
        const B &K(Th.be(ibe));
        int iv[(B::nv) + 2];
        iv[0] = (B::nv);
        for (int ii = 0; ii < (B::nv); ii++) iv[ii + 1] = Th.operator( )(K[ii]);
        if (B::nv==2)
          fprintf(fp, "%d %d %d\n", iv[0], iv[1], iv[2]);
        else if (B::nv==1)
          fprintf(fp, "%d %d\n", iv[0], iv[1]);
      }
    }
  }

  fprintf(fp, "\n");

  // CELL_TYPE
  // print element types in ascii or binary
  int type;
  fprintf(fp, "CELL_TYPES %d\n", numElements);
  if (binary) {

    for (int it = 0; it < Th.nt; it++) {
      if (std::is_same< MMesh, MeshS >::value)
        type = VTK_TRI;
      else if (std::is_same< MMesh, MeshL >::value)
        type = VTK_EDGE;
      if (!bigEndian) FreeFEM::SwapBytes((char *)&type, sizeof(int), 1);
      fwrite(&type, sizeof(int), 1, fp);
    }

    if (surface) {
      for (int ibe = 0; ibe < Th.nbe; ibe++) {
        if (std::is_same< MMesh, MeshS >::value)
          type = VTK_EDGE;
        else if (std::is_same< MMesh, MeshL >::value)
          type = VTK_VERTEX;
        if (!bigEndian) FreeFEM::SwapBytes((char *)&type, sizeof(int), 1);
        fwrite(&type, sizeof(int), 1, fp);
      }
    }
  } else {
    if (std::is_same< MMesh, MeshS >::value)
      type = VTK_TRI;
    else if (std::is_same< MMesh, MeshL >::value)
      type = VTK_EDGE;
    for (int it = 0; it < Th.nt; it++) {
      fprintf(fp, "%d ", type);
    }

    if (surface) {
      if (std::is_same< MMesh, MeshS >::value)
        type = VTK_EDGE;
      else if (std::is_same< MMesh, MeshL >::value)
        type = VTK_VERTEX;
      for (int ibe = 0; ibe < Th.nbe; ibe++) {
        fprintf(fp, "%d ", type);
      }
    }
  }

  fprintf(fp, "\n");

  //= ================================
  // WRITE SOLUTION IN FORMAT VTK
  // LABEL OF ELEMENTS
  //= ================================

  list< int > list_label_Elem;
  for (int it = 0; it < Th.nt; it++) {
    const T &K(Th.t(it));
    list< int >::const_iterator ilist;
    int labOk = 0;

    for (ilist = list_label_Elem.begin( ); ilist != list_label_Elem.end( ); ilist++)
      if (*ilist == K.lab) {
        labOk = 1;
        break;
      }
    if (labOk == 0) list_label_Elem.push_back(K.lab);
  }

  if (surface) {
    for (int ibe = 0; ibe < Th.nbe; ibe++) {
      const B &K(Th.be(ibe));
      list< int >::const_iterator ilist;
      int labOk = 0;

      for (ilist = list_label_Elem.begin( ); ilist != list_label_Elem.end( ); ilist++)
        if (*ilist == K.lab) {
          labOk = 1;
          break;
        }
      if (labOk == 0) list_label_Elem.push_back(K.lab);
    }
  }
  list_label_Elem.sort( );

  //= ================================
  //= ================================

  fprintf(fp, "CELL_DATA %d\n", numElements);
  int cell_fd = 1;
  int cell_lab = 1;
  fprintf(fp, "Scalars  Label int %d\n", cell_fd);
  fprintf(fp, "LOOKUP_TABLE FreeFempp_table\n");
  int label;
  if (binary) {
    for (int it = 0; it < Th.nt; it++) {
      const T &K(Th.t(it));
      label = K.lab;
      if (!bigEndian) FreeFEM::SwapBytes((char *)&label, sizeof(int), 1);
      fwrite(&label, sizeof(int), 1, fp);
    }
    if (surface) {
      for (int ibe = 0; ibe < Th.nbe; ibe++) {
        const B &K(Th.be(ibe));
        label = K.lab;
        if (!bigEndian) FreeFEM::SwapBytes((char *)&label, sizeof(int), 1);
        fwrite(&label, sizeof(int), 1, fp);
      }
    }
  } else {
    for (int it = 0; it < Th.nt; it++) {
      const T &K(Th.elements[it]);
      fprintf(fp, "%d\n", K.lab);
    }

    if (surface) {
      for (int ibe = 0; ibe < Th.nbe; ibe++) {
        const B &K(Th.be(ibe));
        fprintf(fp, "%d\n", K.lab);
      }
    }
  }
  fprintf(fp, "\n");

  int size_list = 0;
  list< int >::const_iterator ilist;

  for (ilist = list_label_Elem.begin( ); ilist != list_label_Elem.end( ); ilist++) {
    size_list++;
  }

  fprintf(fp, "LOOKUP_TABLE FreeFempp_table %d\n", size_list);
  {
    list< int >::const_iterator ilist;

    for (ilist = list_label_Elem.begin( ); ilist != list_label_Elem.end( ); ilist++) {
      if (binary) {
        int tab[4];
        tab[0] = (int)ColorTable[abs(*ilist) % NbColorTable][0] * 255;
        tab[1] = (int)ColorTable[abs(*ilist) % NbColorTable][1] * 255;
        tab[2] = (int)ColorTable[abs(*ilist) % NbColorTable][2] * 255;
        tab[3] = 255;

        for (int itab = 0; itab < 4; itab++) {
          char newvalue[sizeof(int)];
          int bid0 = sprintf(newvalue, "%s", (char *)&tab[itab]);
          fwrite(&newvalue, sizeof(unsigned char), 1, fp);
        }
      } else {
        float tab[4];
        tab[0] = ColorTable[abs(*ilist) % NbColorTable][0];
        tab[1] = ColorTable[abs(*ilist) % NbColorTable][1];
        tab[2] = ColorTable[abs(*ilist) % NbColorTable][2];
        tab[3] = 1.0;

        fprintf(fp, "%.8f %.8f %.8f %.8f\n", tab[0], tab[1], tab[2], tab[3]);
      }
    }
  }
  fprintf(fp, "\n");
}

template< class MMesh >
AnyType VTK_WriteMeshT_Op< MMesh >::operator( )(Stack stack) const {
  string *pffname = GetAny< string * >((*filename)(stack));
  MMesh *pTh = GetAny< MMesh * >((*eTh)(stack));

  ffassert(pTh);
  MMesh &Th = *pTh;
  string *dataname;
  int nbofsol = l.size( );
  KN< int > order(nbofsol);
  char *nameofuser[nbofsol];

  for (int ii = 0; ii < nbofsol; ii++) order[ii] = 0;
  if (nargs[0]) {
    dataname = GetAny< string * >((*nargs[0])(stack));
  }
  bool surface(arg(1, stack, true));
  if (nargs[2]) {
    order = GetAny< KN_< long > >((*nargs[2])(stack));
  }
  bool floatmesh(arg(3, stack, false));
  bool floatsol(arg(4, stack, false));
  bool binary(arg(5, stack, true));
  bool swap(arg(6, stack, isBigEndian( )));
  int datasize = floatmesh ? sizeof(float) : sizeof(double);
  int datasizeSol = floatsol ? sizeof(float) : sizeof(double);

  int iii = 0;
  if (nargs[0]) {
    // char *data = newcopy(dataname->c_str());
    char *data = newcopy(dataname);
    char *name = strtok(data, " \t\n");

    nameofuser[iii] = newcopy(name);
    if (verbosity > 5)
      cout << "   iovtk writeMesh: value of iii  =" << iii << "  \"" << nameofuser[iii] << "\"\n";
    iii++;

    while ((name = strtok(NULL, " \t\n\0"))) {
      if (iii >= nbofsol) {
        if (verbosity > 5)
          cout << " iovtk writeMesh: The number of data name is too large " << endl;
        break;
      }
      nameofuser[iii] = newcopy(name);
      if (verbosity > 5)
        cout << "   iovtk writeMesh: value of iii  =" << iii << "  \"" << nameofuser[iii] << "\"\n";
      iii++;
    }

    if (iii < nbofsol)
      if (verbosity)
        cout << "   iovtk writeMesh:  The number of data name is too small, we give default name "
             << endl;
  }

  if (iii < nbofsol)
    for (int iiii = iii; iiii < nbofsol; iiii++)
      nameofuser[iiii] = newcopy(l[iiii].name.c_str( ));    // dataff;

  // determination of number of order 0 et 1.
  int Norder0 = 0;

  for (int ii = 0; ii < nbofsol; ii++)
    if (order[ii] == 0) Norder0++;

  // lecture du nom des variables

  FILE *fp = fopen((*pffname).c_str( ), "wb");
  // FILE *fp = fopen( newcopy(pffname), "wb");
  if (!fp) {
    cerr << "Unable to open file " << (*pffname).c_str( ) << endl;
    ExecError("error in reading vtk file");
  }

  // type of VTK FILE
  char *pch = newcopy(pffname);
  int VTK_FILE = 0;
  int ls = 0;
  int lll = strlen(pch);
  if (!strcmp(pch + lll - (ls = 4), ".vtk"))
    VTK_FILE = 1;
  else if (!strcmp(pch + lll - (ls = 4), ".vtu"))
    VTK_FILE = 2;

  if (verbosity) cout << " " << pffname << " VTK_FILE " << VTK_FILE << endl;

  if (VTK_FILE == 1) {
    // CAS VTK
    VTK_WRITE_MESHT< MMesh >(*pffname, fp, Th, binary, datasize, surface, swap);

    // WRITE SOLUTIONS
    if (datasizeSol == sizeof(float)) {
      if (Norder0 > 0) {
        fprintf(fp, "FIELD FieldData %d\n", Norder0);

        for (int ii = 0; ii < nbofsol; ii++) {
          if (order[ii] == 0) {
            int nsol;
            if (surface) {
              nsol = Th.nt + Th.nbe;
            } else {
              nsol = Th.nt;
            }

            fprintf(fp, "%s %ld %d float\n", nameofuser[ii], l[ii].nbfloat, nsol);
            if (verbosity > 5)
              cout << "   iovtk writeMesh: name of data(" << ii << ")=" << nameofuser[ii] << endl;
            // changement ecriture solution
            l[ii].writesolutionP0_float(fp, Th, stack, surface, binary, swap);
          }
        }
      }

      if (Norder0 < nbofsol) {
        fprintf(fp, "POINT_DATA %d\n", Th.nv);
        fprintf(fp, "FIELD FieldData %d\n", nbofsol - Norder0);

        for (int ii = 0; ii < nbofsol; ii++) {
          if (order[ii] == 1) {
            // fprintf(fp,"%s %d %d float\n",l[ii].name.c_str(),l[ii].nbfloat,Th.nv);
            fprintf(fp, "%s %ld %d float\n", nameofuser[ii], l[ii].nbfloat, Th.nv);
            if (verbosity > 5) {
              cout << "iovtk writeMeshS: name of data(" << ii << ")=" << nameofuser[ii] << " "
                   << l[ii].name << endl;
            }

            // changement ecriture solution
            l[ii].writesolutionP1_float(fp, Th, stack, binary, swap);
          }
        }
      }
    }

    if (datasizeSol == sizeof(double)) {
      if (Norder0 > 0) {
        fprintf(fp, "FIELD FieldData %d\n", Norder0);

        for (int ii = 0; ii < nbofsol; ii++) {
          if (order[ii] == 0) {
            int nsol;
            if (surface) {
              nsol = Th.nt + Th.nbe;
            } else {
              nsol = Th.nt;
            }

            fprintf(fp, "%s %ld %d double\n", nameofuser[ii], l[ii].nbfloat, nsol);
            if (verbosity > 5)
              cout << "   iovtk writeMeshS:name of data(" << ii << ")=" << nameofuser[ii] << endl;

            // changement ecriture solution
            l[ii].writesolutionP0_double(fp, Th, stack, surface, binary, swap);
          }
        }
      }
      if (Norder0 < nbofsol) {
        fprintf(fp, "POINT_DATA %d\n", Th.nv);
        fprintf(fp, "FIELD FieldData %d\n", nbofsol - Norder0);

        for (int ii = 0; ii < nbofsol; ii++) {
          if (order[ii] == 1) {
            fprintf(fp, "%s %ld %d double\n", nameofuser[ii], l[ii].nbfloat, Th.nv);
            if (verbosity > 5) {
              cout << "iovtk writeMeshS: name of data(" << ii << ")=" << nameofuser[ii] << endl;
            }

            // changement ecriture solution
            l[ii].writesolutionP1_double(fp, Th, stack, binary, swap);
          }
        }
      }
    }
  } else if (VTK_FILE == 2) {

    VTU_WRITE_MESHT< MMesh >(fp, Th, binary, datasize, surface);

    // Solution Order
    // order 1
    if (Norder0 != nbofsol) {
      BEGINTYPE_VTU(fp, "PointData");

      for (int ii = 0; ii < nbofsol; ii++) {
        if (order[ii] == 1) {
          if (datasize == sizeof(float)) {
            VTU_DATA_ARRAY(fp, "Float32", nameofuser[ii], l[ii].nbfloat, binary);
            l[ii].writesolutionP1_float(fp, Th, stack, binary, swap, 1);
          } else if (datasize == sizeof(double)) {
            VTU_DATA_ARRAY(fp, "Float64", nameofuser[ii], l[ii].nbfloat, binary);
            l[ii].writesolutionP1_double(fp, Th, stack, binary, swap, 1);
          }

          ENDTYPE_VTU(fp, "DataArray");
        }
      }

      ENDTYPE_VTU(fp, "PointData");
    }

    // order 0
    if (Norder0 > 0) {
      BEGINTYPE_VTU(fp, "CellData");

      for (int ii = 0; ii < nbofsol; ii++) {
        if (order[ii] == 0) {
          if (datasize == sizeof(float)) {
            VTU_DATA_ARRAY(fp, "Float32", nameofuser[ii], l[ii].nbfloat, binary);
            l[ii].writesolutionP0_float(fp, Th, stack, surface, binary, swap, 1);
          } else if (datasize == sizeof(double)) {
            VTU_DATA_ARRAY(fp, "Float64", nameofuser[ii], l[ii].nbfloat, binary);
            l[ii].writesolutionP0_double(fp, Th, stack, surface, binary, swap, 1);
          }

          ENDTYPE_VTU(fp, "DataArray");
        }
      }

      ENDTYPE_VTU(fp, "CellData");
    }

    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
  } else {
    cout << " iovtk extension file is not correct (" << VTK_FILE << " != 1 or 2 ) " << endl;
    ExecError(" iovtk : extension file");
  }

  // close file fp
  fclose(fp);

  for (int iiii = 0; iiii < nbofsol; iiii++) {
    delete[] nameofuser[iiii];
  }

  delete[] pch;
  return (MMesh *)NULL;
}

//= =============================================
// FIN ECRITURE DE FICHIER .vtk (3D SURFACE)
//= =============================================

template< class MMesh >
class VTK_LoadMeshT_Op : public E_F0mps {
 public:
  Expression filename;
  static const int n_name_param = 9;    //
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  int arg(int i, Stack stack, int a) const {
    return nargs[i] ? GetAny< int >((*nargs[i])(stack)) : a;
  }
  bool arg(int i, Stack stack, bool a) const {
    return nargs[i] ? GetAny< bool >((*nargs[i])(stack)) : a;
  }
  double arg(int i, Stack stack, double a) const {
    return nargs[i] ? GetAny< double >((*nargs[i])(stack)) : a;
  }
  KN<KN<double> >* arg(int i, Stack stack,KN<KN<double> >*p) const {
    return nargs[i] ? GetAny< KN<KN<double> >* >((*nargs[i])(stack)) : p;
  }
 
 public:
  VTK_LoadMeshT_Op(const basicAC_F0 &args, Expression ffname) : filename(ffname) {
    if (verbosity) {
      cout << "Load mesh given by VTK " << endl;
    }

    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator( )(Stack stack) const;
};

template<>
basicAC_F0::name_and_type VTK_LoadMeshT_Op< MeshS >::name_param[] = {
  {"reftri", &typeid(long)},         
  {"swap", &typeid(bool)},
  {"refedge", &typeid(long)},        
  {"namelabel", &typeid(string)},
  {"cleanmesh", &typeid(bool)},      
  {"removeduplicate", &typeid(bool)},
  {"precisvertice", &typeid(double)},
  {"ridgeangledetection", &typeid(double)},
  {"fields", &typeid( KN<KN<double> >*)}
};

template<>
basicAC_F0::name_and_type VTK_LoadMeshT_Op< MeshL >::name_param[] = {
  {"refedge", &typeid(long)},        
  {"swap", &typeid(bool)},
  {"refbdpoint", &typeid(long)},     
  {"namelabel", &typeid(string)},
  {"cleanmesh", &typeid(bool)},      
  {"removeduplicate", &typeid(bool)},
  {"precisvertice", &typeid(double)},
  {"ridgeangledetection", &typeid(double)},
  {"fields", &typeid( KN<KN<double> >*)}
};

template< class MMesh >
class VTK_LoadMeshT : public OneOperator {
 public:
  typedef const MMesh *ppmesh;
  VTK_LoadMeshT( ) : OneOperator(atype< ppmesh >( ), atype< string * >( )) {}

  E_F0 *code(const basicAC_F0 &args) const {
    return new VTK_LoadMeshT_Op< MMesh >(args, t[0]->CastTo(args[0]));
  }
};

template< class MMesh >
MMesh *VTK_LoadT(const string &filename, bool bigEndian, bool cleanmesh, bool removeduplicate,
                 double precisvertice, double ridgeangledetection, KN<KN<double> >* pfields) {
  // swap = bigEndian or not bigEndian
  // variable freefem++
  typedef typename MMesh::Element T;
  typedef typename MMesh::BorderElement B;
  typedef typename MMesh::Vertex V;

  int nv, nt = 0, nbe = 0;
  char *res;
  // Reading Mesh in vtk formats

  FILE *fp = fopen(filename.c_str( ), "rb");
  if (!fp) {
    cerr << "Unable to open file " << filename.c_str( ) << endl;
    exit(1);
  }

  char buffer[256], buffer2[256], buffer3[256];;

  res = fgets(buffer, sizeof(buffer), fp);    // version line
  res = fgets(buffer, sizeof(buffer), fp);    // title

  fscanf(fp, "%s", buffer);    // ASCII or BINARY
  bool binary = false;
  if (!strncmp(buffer, "BINARY", 6)) binary = true;

  if (fscanf(fp, "%s %s", buffer, buffer2) != 2) {
    cerr << "error in reading vtk files" << filename << endl;
    ExecError("error in reading vtk file");
  }

  if (strcmp(buffer, "DATASET") || strcmp(buffer2, "UNSTRUCTURED_GRID")) {
    cout << "VTK reader can only read unstructured datasets" << endl;
    ExecError("error in reading vtk file");
    exit(1);
  }

  // read mesh vertices
  if (fscanf(fp, "%s %d %s\n", buffer, &nv, buffer2) != 3) {
    cout << "error in reading vtk files" << endl;
    ExecError("error in reading vtk file");
    exit(1);
  }

  if (strcmp(buffer, "POINTS") || !nv) {
    cerr << "No points in dataset" << endl;
    ExecError("error in reading vtk file: No points in dataset");
    exit(1);
  }

  int datasize;
  if (!strncmp(buffer2, "double", 6))
    datasize = sizeof(double);
  else if (!strncmp(buffer2, "float", 5))
    datasize = sizeof(float);
  else {
    cout << "VTK reader only accepts float or double datasets" << endl;
    ExecError("error in reading vtk file: VTK reader only accepts float or double datasets");
    exit(1);
  }

  if (verbosity > 3)
    cout << "Reading points" << nv << ", buffer2 " << buffer2 << ", binary " << binary << " "
         << datasize << " " << sizeof(float) << endl;

  V *vff = new V[nv];

  for (int i = 0; i < nv; i++) {
    if (verbosity > 9) cout << " i=" << i << endl;
    double xyz[3];
    if (binary) {
      if (datasize == sizeof(float)) {
        float f[3];
        if (fread(f, sizeof(float), 3, fp) != 3) {
          ExecError("error in reading vtk file: VTK reader only accepts float or double datasets");
          ExecError("error in reading vtk file");
        }

        if (!bigEndian) FreeFEM::SwapBytes((char *)f, sizeof(float), 3);

        for (int j = 0; j < 3; j++) xyz[j] = f[j];

      } else if (fread(xyz, sizeof(double), 3, fp) != 3) {
        cout << "error in reading vtk files" << endl;
        ExecError("error in reading vtk file");
      }

      if (!bigEndian) FreeFEM::SwapBytes((char *)xyz, sizeof(double), 3);

    } else {
      if (verbosity > 10) cout << datasize << " " << sizeof(float) << endl;
      if (datasize == sizeof(float)) {
        if (fscanf(fp, "%lf %lf %lf", &xyz[0], &xyz[1], &xyz[2]) != 3) {
          cout << "error in reading vtk files (float)" << endl;
          ExecError("error in reading vtk file");
        }
      } else if (fscanf(fp, "%lf %lf %lf", &xyz[0], &xyz[1], &xyz[2]) != 3) {
        cout << "error in reading vtk files" << endl;
        ExecError("error in reading vtk file");
      }
    }

    vff[i].x = xyz[0];
    vff[i].y = xyz[1];
    vff[i].z = xyz[2];
    vff[i].lab = 1;

    if (verbosity > 9) {
      printf("xyz = %f %f %f\n", xyz[0], xyz[1], xyz[2]);
    }
  }

  // read mesh elements
  int numElements, numElements2, totalNumInt;
  if (fscanf(fp, "%s %d %d\n", buffer, &numElements, &totalNumInt) != 3) {
    cout << "error in reading vtk files" << endl;
    ExecError("error in reading vtk file");
  }

  if (verbosity > 3) printf("reading parameter %s %d %d\n", buffer, numElements, totalNumInt);

  if (strcmp(buffer, "CELLS") || !numElements) {
    cout << "No cells in dataset" << endl;
    ExecError("error in reading vtk file");
  }

  if (verbosity > 3) cout << "Reading cells" << numElements << endl;

  int *IntCells = new int[totalNumInt - numElements];
  int *firstCell = new int[numElements + 1];
  int *TypeCells = new int[numElements];
  int numIntCells = 0;

  for (unsigned int i = 0; i < numElements; i++) {
    int numVerts, n[100];
    if (verbosity > 9) cout << "i=" << i << " " << numElements << endl;

    for (int ii = 0; ii < 100; ii++) {
      n[ii] = -1;
    }

    if (binary) {
      if (fread(&numVerts, sizeof(int), 1, fp) != 1) {
        cout << "error in reading VTK files " << endl;
        ExecError("error in reading vtk file");
      }

      if (!bigEndian) FreeFEM::SwapBytes((char *)&numVerts, sizeof(int), 1);

      if ((int)fread(n, sizeof(int), numVerts, fp) != numVerts) {
        cout << "error in reading VTK files " << endl;
        ExecError("error in reading vtk file");
      }

      if (!bigEndian) FreeFEM::SwapBytes((char *)n, sizeof(int), numVerts);
    } else {
      if (fscanf(fp, "%d", &numVerts) != 1) {
        cout << "error in reading VTK files " << endl;
        ExecError("error in reading vtk file");
      }
      if (verbosity > 9) cout << "numVerts " << numVerts << endl;

      for (int j = 0; j < numVerts; j++) {
        if (fscanf(fp, "%d", &n[j]) != 1) {
          cout << "error in reading VTK files " << endl;
          ExecError("error in reading vtk file");
        }
        if (verbosity > 9) cout << "  n[j]" << n[j] << endl;
      }
    }

    firstCell[i] = numIntCells;

    for (int j = 0; j < numVerts; j++) {
      if (n[j] >= 0 && n[j] < nv) {
        IntCells[numIntCells] = n[j];
        numIntCells++;
      } else {
        cout << "Bad vertex index" << endl;
        ExecError("error in reading vtk file");
      }
    }
  }

  firstCell[numElements] = totalNumInt - numElements;

  if (fscanf(fp, "%s %d\n", buffer, &numElements2) != 2) {
    cout << " Error in reading CELL_TYPES ARGUMENT " << endl;
    ExecError("error in reading vtk file");
  }

  if (strcmp(buffer, "CELL_TYPES") || numElements2 != (int)numElements) {
    cout << "No or invalid number of cells types" << endl;
    ExecError("error in reading vtk file");
  }

  if (verbosity > 3) printf("reading parameter %s %d\n", buffer, numElements2);

  for (unsigned int i = 0; i < numElements; i++) {
    int type;
    if (binary) {
      if (fread(&type, sizeof(int), 1, fp) != 1) {
        cout << "bug in readings cell types" << endl;
        ExecError("error in reading vtk file");
      }

      if (!bigEndian) {
        FreeFEM::SwapBytes((char *)&type, sizeof(int), 1);
      }
    } else {
      if (fscanf(fp, "%d", &type) != 1) {
        cout << "bug in readings cell types" << endl;
        ExecError("error in reading vtk file");
      }
    }

    TypeCells[i] = type;

    switch (type) {
      case 1:    // Vertex
        if (std::is_same< MMesh, MeshL >::value) nbe++;
        break;

      case 3:    // Edge/line
        if (std::is_same< MMesh, MeshS >::value)
          nbe++;
        else if (std::is_same< MMesh, MeshL >::value)
          nt++;
        break;
      case 5:    // Triangle
        if (std::is_same< MMesh, MeshS >::value)
          nt++;
        else if (std::is_same< MMesh, MeshL >::value) {
          cout << "We are loading a three dimensional CURVE mesh. Three is no triangle." << endl;
          ExecError("error in reading vtk file");
        }
        break;
      case 10:    // Tetrah�dre
        if (std::is_same< MMesh, MeshS >::value)
          cout << "We are loading a three dimensional SURFACE mesh. Three is no tetrahedron."
               << endl;
        else if (std::is_same< MMesh, MeshL >::value)
          cout << "We are loading a three dimensional CURVE mesh. Three is no tetrahedron." << endl;
        ExecError("error in reading vtk file");
        break;
      default:
        cout << "Error :: This type of cell is not considered in Freefem++ " << type << endl;
        ExecError("error in reading vtk file");
        break;
    }
  }

  int nbp=0,nbf=0, err=0;
    if (fscanf(fp, "%s %d", buffer, &nbp) != 2)
        {   cout << "error in reading vtk files pfields" << endl;
            err++;}
    int startdatapoint=0;
    if(err==0)
    {
        int nf=-1;
        /*
         CELL_DATA 209726
         Scalars  Label int 1
         LOOKUP_TABLE FreeFempp_table
         ....
         LOOKUP_TABLE FreeFempp_table 7
         4*7 value
         */
        
        if (strcmp(buffer, "CELL_DATA"))  { //  read region number if exist
            if (strcmp(buffer, "POINT_DATA"))  {
                cout << "VTK reader can only read CELL_DATA or POINT_DATA datasets:  not " << buffer<< " " << nbp<<  endl;
                err=1;
            }
            else startdatapoint=1;
        }
        else {
            if ((!err) &&(fscanf(fp, "%s %s %s %d\n", buffer, buffer2,buffer3,&nbf) != 4)) {
                cout << "error in reading vtk files FIELD FieldData" << endl;
                err++;
            }}
        
        if( strcmp(buffer3, "int") !=0)// not integer
            err++;
        if ((!err) &&(fscanf(fp, "%s %s\n", buffer, buffer2) != 2))
            err++;
        // read nbf
        cout << " err= " << err << " read nbp "<< nbp << endl;
        if(err==0)
            for( nf=0 ; nf < nbp; nf++)
            {
                int ii[1];
                if (binary)
                { if (fread(ii, sizeof(int), 1, fp) != 1) err++;}
                else
                {      if (fscanf(fp, "%d", ii) != 1) err++;}
                if(err) break;
            }
        if(err) cout << " err reading CELL_DATA  at " << nf << endl;
        if ((!err) &&(fscanf(fp, "%s %s %d\n", buffer, buffer2,&nbf) != 3 ) ) err++;
        nf =-1;
        if(err==0)
            for( nf=0 ; nf < nbf; nf++)
            {
                float f[4];
                char cc[4];
                if (binary)
                { if (fread(cc, sizeof(char), 4, fp) != 4) err++;}
                else
                {      if (fscanf(fp, "%f %f %f %fa",f+0,f+1,f+2,f+3) != 4) err++;}
                if(err) break;
            }
        if(err&& nf>=0) cout << " err LOOKUP_TABLE FreeFempp_table at " << nf << " " << err << endl;
        
        startdatapoint=0;
    }
           
    
  if(pfields && err==0) {
    if(verbosity>1)   cout << " try  reading POINT_DATA  " << startdatapoint << endl;
    /*
    POINT_DATA 32436
    FIELD FieldData 2
    Velocity 3 32436 float
      
    */
    nbp =0;//  no POINT_DATA
    nbf=0;
    if (startdatapoint==0)
      if (fscanf(fp, "%s", buffer) != 1) {
        cout << "error in reading vtk files pfields" << endl;
        err++;
      }
        
      if (strcmp(buffer, "POINT_DATA")==0) {
        if (fscanf(fp, "%d", &nbp) != 1) err++;
        if (fscanf(fp, "%s", buffer) != 1) err++;
      }
      if(err == 0 && strcmp(buffer, "FIELD")!=0) {
        cout << "VTK reader can only read FIELD/POINT_DATA datasets:  not " << buffer<<  endl;
        err++;
      }
        
      if ((!err) &&(fscanf(fp, "%s %d", buffer2,&nbf) != 2)) {
        cout << "error in reading vtk files FIELD FieldData" << endl;
        err++;
      }
        
      if( err) nbf=0;
      else pfields->resize(nbf);
      for(int nf=0 ; nf < nbf; nf++) {
                
        int m,nv;
        // read mesh vertices
        if (fscanf(fp, "%s %d %d %s\n", buffer,&m, &nv, buffer2) != 4) {
          cout << "error in reading vtk files " << endl;
          err++;
          break;
        }
        int n = m*nv;
        if(verbosity) cout << " reading "<< buffer << " "<< m << " " << nv << " "<< buffer2 << endl;
        int datasize;
        if (!strncmp(buffer2, "double", 6))
          datasize = sizeof(double);
        else if (!strncmp(buffer2, "float", 5))
           datasize = sizeof(float);
        else {
          cout << "VTK reader only accepts float or double datasets" << endl;
          err++;
          break;
        }
        //  read data ..
        if(err) break;
        (*pfields)[nf].resize(n);
        double* pv=&(*pfields)[nf][0];
        for(int i=0; i<n; ++i) {
          if (binary) {
            double f[1];
            if (datasize == sizeof(float)){
              if (fread(f, sizeof(float), 1, fp) != 1) {
                cout << "error in reading  vtk fields float at " << nf << " / " << i << endl;
                err++;
                break;
              }
            }
            else {
              if (fread(f, sizeof(double), 1, fp) != 1) {
                cout << "error in reading  vtk fields double at " << nf << " / " << i << endl;
                err++;
                break;
              }
            }
            if (!bigEndian)
              FreeFEM::SwapBytes((char *)f, sizeof(float), 1);
            *pv++=*f;
        
          }
          else
            if (fscanf(fp, "%lf", pv++) != 1) {
              cout << "error in reading vtk files ascii fields" << nf << " / " << i <<endl;
              err++;
            }
        if(err) break;
        }
        if(err) break;
      }
    if( err ) (*pfields).resize(0);
  } // end if(pfields)
  fclose(fp);

        
  T *tff;
  if (nt > 0)
    tff = new T[nt];
  else
    ExecError("error in reading vtk file: Not element");
  T *ttff = tff;
  B *bff;

  if (nbe > 0)
    bff = new B[nbe];
  else
    if (!std::is_same< MMesh, Mesh3 >::value) bff=NULL;
    else ExecError("error in reading vtk file: Not border element");
  B *bbff = bff;

  for (unsigned int i = 0; i < numElements; i++) {
    int type = TypeCells[i];
    int ivb[B::nv], ivt[T::nv];
    int label = 1;

    switch (type) {
      case 1:    // Vertex
        assert((firstCell[i + 1] - firstCell[i]) == 1);
        for (int j = firstCell[i]; j < firstCell[i + 1]; j++) ivb[j - firstCell[i]] = IntCells[j];
        if (std::is_same< MMesh, MeshL >::value) (bbff++)->set(vff, ivb, label);
        break;
      case 3:    // Edge
        assert((firstCell[i + 1] - firstCell[i]) == 2);
        if (std::is_same< MMesh, MeshS >::value) {
          for (int j = firstCell[i]; j < firstCell[i + 1]; j++) ivb[j - firstCell[i]] = IntCells[j];
          (bbff++)->set(vff, ivb, label);
        } else if (std::is_same< MMesh, MeshL >::value) {
          for (int j = firstCell[i]; j < firstCell[i + 1]; j++) ivt[j - firstCell[i]] = IntCells[j];
          (ttff++)->set(vff, ivt, label);
        }
        break;

      case 5:    // Triangle
        assert((firstCell[i + 1] - firstCell[i]) == 3);
        for (int j = firstCell[i]; j < firstCell[i + 1]; j++) ivt[j - firstCell[i]] = IntCells[j];
        if (std::is_same< MMesh, MeshS >::value) (ttff++)->set(vff, ivt, label);
        break;
      default:
        break;
    }
  }

  delete[] IntCells;
  delete[] firstCell;
  delete[] TypeCells;

  MMesh *pTh = new MMesh(nv, nt, nbe, vff, tff, bff, cleanmesh || (nbe==0), removeduplicate, (nbe==0), precisvertice);
  return pTh;
}

template< class MMesh >
AnyType VTK_LoadMeshT_Op< MMesh >::operator( )(Stack stack) const {
  string *pffname = GetAny< string * >((*filename)(stack));

  bool swap(arg(1, stack, false));
  string *DataLabel;
  if (nargs[3]) {
    DataLabel = GetAny< string * >((*nargs[3])(stack));
  }

  bool cleanmesh(arg(4, stack, false));
  bool removeduplicate(arg(5, stack, false));
  double precisvertice(arg(6, stack, 1e-6));
  double ridgeangledetection(arg(7, stack, 8.*atan(1.)/9.));
  KN<KN<double> > * pfields=0;
  pfields=arg(8, stack,pfields);
  
  MMesh *Th = VTK_LoadT< MMesh >(*pffname, swap, cleanmesh, removeduplicate, precisvertice,ridgeangledetection, pfields);
  // TODO Axel
  // MMesh *Th = VTK_LoadT< MMesh >(*pffname, swap, cleanmesh || (nbe==0), removeduplicate, (nbe==0), precisvertice,ridgeangledetection, pfields);

  // A faire fonction pour changer le label

  Add2StackOfPtr2FreeRC(stack, Th);

  return Th;
}

//= ======================
// FIN 3D SURFACE Fichier .vtk
//= ======================

void saveMatlab(const string &file, const Mesh &Th) {
  // ErrorInsaveMesh e;
  {
    ofstream pf(file.c_str( ));
    ffassert(pf);
    typedef Mesh::Element Element;

    for (int k = 0; k < Th.nt; ++k) {
      Element &K = Th[k];
      pf << "x = [ ";

      for (size_t n = 0; n < 3; n++) {
        pf << std::setprecision(5) << setw(18) << K[n].x << " ";
      }

      pf << std::setprecision(5) << setw(18) << K[0].x << " ]; ";
      pf << "y = [ ";

      for (size_t n = 0; n < 3; n++) {
        pf << std::setprecision(5) << setw(18) << K[n].y << " ";
      }

      pf << std::setprecision(5) << setw(18) << K[0].y << " ]; ";
      pf << "line(x,y);" << endl;
    }

    pf.close( );
  }
}

void saveTecplot(const string &file, const Mesh &Th) {
  string shape;
  ofstream pf(file.c_str( ));
  size_t n, m;

  pf << "TITLE = \" \"\n";
  pf << "VARIABLES = \"X\", \"Y\"";
  if (Th.dim == 3) {
    pf << ", \"Z\"";
  }

  pf << endl;
  if (Th.dim == 2) {
    m = 3;
    shape = "TRIANGLE";
  }
  /*      else if (el->getShape()==LINE) {
   *      m = 2;
   *      shape = "LINESEG";
   *      }
   */
  else if (Th.dim == 3) {
    m = 4;
    shape = "TETRAHEDRON";
  }

  pf << "ZONE N=" << Th.nv << ", E=" << Th.nt << ", F=FEPOINT, ET=" << shape << endl;

  for (int i = 0; i < Th.nv; i++) {
    pf << std::setprecision(5) << setw(18) << (R2)Th(i) << " \n";
  }

  for (int k = 0; k < Th.nt; ++k) {
    for (n = 0; n < m; n++) {
      pf << Th(k, n) + 1 << "  ";
    }

    pf << endl;
  }

  pf.close( );
}

/*  class Init1 { public:
 * Init1();
 * };
 *
 * $1 */

#ifndef COMMON_HPDDM_PARALLEL_IO
static void Load_Init( ) {    // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
  typedef Mesh *pmesh;
  typedef Mesh3 *pmesh3;
  typedef MeshS *pmeshS;
  typedef MeshL *pmeshL;

  // if (verbosity)
  if (verbosity && (mpirank == 0)) {
    cout << " load: iovtk " << endl;
  }

  if (!Global.Find("savevtk").NotNull( )) {
    Global.Add("savevtk", "(", new OneOperatorCode< VTK_WriteMesh_Op >);
    Global.Add("savevtk", "(", new OneOperatorCode< VTK_WriteMesh3_Op >);
  }
  Global.Add("savevtk", "(", new OneOperatorCode< VTK_WriteMeshT_Op< MeshS > >);
  Global.Add("savevtk", "(", new OneOperatorCode< VTK_WriteMeshT_Op< MeshL > >);
  Global.Add("vtkload", "(", new VTK_LoadMesh);
  Global.Add("vtkload3", "(", new VTK_LoadMesh3);
  Global.Add("vtkloadS", "(", new VTK_LoadMeshT< MeshS >);
  Global.Add("vtkloadL", "(", new VTK_LoadMeshT< MeshL >);
}

LOADFUNC(Load_Init)
#endif
