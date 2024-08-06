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
// SUMMARY : Gmsh - FreeFem++ link
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
//
// Vincent HUBER - vincent.huber@cemosis.fr   October 2014
// manage verbosity levels
// Add SaveGMH "d  FH
//
#include "ff++.hpp"

using namespace Fem2D;

// Table of number of vertex for an element type of gmsh
static const int nvElemGmsh[30] = {2, 3, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
// we considerer only edges, triangles and tetrahedrons in Freefem++
// 15 :: Vertex Corner
// 1  :: Edge/line
// 2  :: triangles
// 4  :: tetrahedrons
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

class GMSH_LoadMesh_Op : public E_F0mps {
 public:
  Expression filename;
  static const int n_name_param = 2;    //
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];

 public:
  GMSH_LoadMesh_Op(const basicAC_F0 &args, Expression ffname) : filename(ffname) {
    if (verbosity > 1) {
      cout << "Load mesh given by GMSH " << endl;
    }

    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type GMSH_LoadMesh_Op::name_param[] = {{"reftri", &typeid(long)},
                                                            {"renum", &typeid(long)}};

class GMSH_LoadMesh : public OneOperator {
 public:
  GMSH_LoadMesh( ) : OneOperator(atype< pmesh >( ), atype< string * >( )) {}

  E_F0 *code(const basicAC_F0 &args) const {
    return new GMSH_LoadMesh_Op(args, t[0]->CastTo(args[0]));
  }
};

Mesh *GMSH_Load(const string &filename) {
  // variable freefem++
  int nv, nt = 0, nbe = 0, ret;
  Mesh::Vertex *vff;

  map< int, int > mapnumv;

  // loading mesh and reading mesh in gmsh are in the file GModelIO_Mesh.cpp (directory Geo)
  char str[256] = "ZZZ", *res;
  double version = 2.0;
  bool binary = false, swap = false, postpro = false;
  FILE *fp = fopen(filename.c_str( ), "rb");
  if (!fp) {
    cerr << "Unable to open file " << filename.c_str( ) << endl;
    exit(1);
  }

  while (!feof(fp)) {
    res = fgets(str, sizeof(str), fp);
    if (verbosity > 10 && res == nullptr) cout << "fgets error or EOF" << endl;
    if (str[0] == '$') {
      if (!strncmp(&str[1], "MeshFormat", 10)) {
        if (!fgets(str, sizeof(str), fp)) {
          exit(1);
        }

        int format, size;
        if (sscanf(str, "%lf %d %d", &version, &format, &size) != 3) {
          exit(1);
        }

        if (verbosity > 1) {
          cout << "Mesh Format is " << format << endl;
        }

        if (format) {
          binary = true;
          if (verbosity > 2) {
            cout << "Mesh is in binary format" << endl;
          }

          int one;
          if (fread(&one, sizeof(int), 1, fp) != 1) {
            exit(1);
          }

          if (one != 1) {
            swap = true;
            if (verbosity > 2) {
              cout << "Swapping bytes from binary file" << endl;
            }
          }
        }
      } else if (!strncmp(&str[1], "PhysicalNames", 13)) {
        if (verbosity > 1) {
          cout << " PhysicalNames is not considered in freefem++ " << endl;
        }
      } else if (!strncmp(&str[1], "NO", 2) || !strncmp(&str[1], "Nodes", 5) ||
                 !strncmp(&str[1], "ParametricNodes", 15)) {
        const bool parametric = !strncmp(&str[1], "ParametricNodes", 15);
        if (parametric == true) {
          cerr << " ParametricNodes is not considered yet in freefem++" << endl;
          exit(1);
        }

        if (!fgets(str, sizeof(str), fp)) {
          exit(1);
        }

        if (sscanf(str, "%d", &nv) != 1) {
          exit(1);
        }

        if (verbosity > 1) {
          printf("%d vertices\n", nv);
        }

        // local variables freefem++
        vff = new Mesh::Vertex[nv];

        for (int i = 0; i < nv; i++) {
          int num;
          double xyz[3] = {0}, uv[2];

          // if (!parametric){
          if (!binary) {
            if (fscanf(fp, "%d %lf %lf %lf", &num, &xyz[0], &xyz[1], &xyz[2]) != 4) {
              exit(1);
            }
          } else {
            if (fread(&num, sizeof(int), 1, fp) != 1) {
              exit(1);
            }

            if (swap) {
              SwapBytes((char *)&num, sizeof(int), 1);
            }

            if (fread(xyz, sizeof(double), 3, fp) != 3) {
              exit(1);
            }

            if (swap) {
              SwapBytes((char *)xyz, sizeof(double), 3);
            }
          }

          assert(abs(xyz[2]) < 1.e-7);
          vff[i].x = xyz[0];
          vff[i].y = xyz[1];
          vff[i].lab = 1;
          mapnumv[num] = i;
        }
      } else if (!strncmp(&str[1], "ELM", 3) || !strncmp(&str[1], "Elements", 8)) {
        if (!fgets(str, sizeof(str), fp)) {
          exit(1);
        }

        int numElements;
        sscanf(str, "%d", &numElements);

        if (verbosity > 2) {
          cout << "Loading elements\n" << endl;
        }

        if (!binary) {
          for (int i = 0; i < numElements; i++) {
            int num, type, physical = 0, elementary = 0, partition = 0, numVertices;
            if (version <= 1.0) {
              ret = fscanf(fp, "%d %d %d %d %d", &num, &type, &physical, &elementary, &numVertices);
              if (ret == EOF) cerr << "fscanf error" << endl;
            } else {
              int numTags;
              ret = fscanf(fp, "%d %d %d", &num, &type, &numTags);
              if (ret == EOF) cerr << "fscanf error" << endl;
              for (int j = 0; j < numTags; j++) {
                int tag;
                ret = fscanf(fp, "%d", &tag);
                if (ret == EOF) cerr << "fscanf error" << endl;
                if (j == 0) {
                  physical = tag;
                } else if (j == 1) {
                  elementary = tag;
                } else if (j == 2) {
                  partition = tag;
                }

                // ignore any other tags for now
              }
	      if (verbosity > 99) cout << type << endl;
              ffassert(type >= 1 && type <= 31);
              if ((numVertices = nvElemGmsh[type - 1]) == 0) {
                cerr << "Element of type " << type << " is not considered in Freefem++" << endl;
                exit(1);
              }
            }

            if (type == 1) {
              nbe++;
            }

            if (type == 2) {
              nt++;
            }

            if (type == 4) {
              if (verbosity > 1) cerr << "We are loading a two dimensional mesh " << endl;
              exit(1);
            }

            int indices[60];

            for (int j = 0; j < numVertices; j++) {
              ret = fscanf(fp, "%d", &indices[j]);
              if (ret == EOF) cerr << "fscanf error" << endl;
            }
          }
        } else {
          int numElementsPartial = 0;

          while (numElementsPartial < numElements) {
            int header[3] = {};
            if (fread(header, sizeof(int), 3, fp) != 3) {
              exit(1);
            }

            if (swap) {
              SwapBytes((char *)header, sizeof(int), 3);
            }

            int type = header[0];
            int numElms = header[1];
            int numTags = header[2];
            int numVertices;
            assert(type >= 1 && type <= 31);
            if ((numVertices = nvElemGmsh[type - 1]) == 0) {
              cerr << "Element of type " << type << " is not considered in FreeFem++" << endl;
              exit(1);
            }

            unsigned int n = 1 + numTags + numVertices;
            int *data = new int[n];

            for (int i = 0; i < numElms; i++) {
              if (fread(data, sizeof(int), n, fp) != n) {
                exit(1);
              }

              if (swap) {
                SwapBytes((char *)data, sizeof(int), n);
              }

              int physical = (numTags > 0) ? data[4 - numTags] : 0;
              int elementary = (numTags > 1) ? data[4 - numTags + 1] : 0;
              int partition = (numTags > 2) ? data[4 - numTags + 2] : 0;
              int *indices = &data[numTags + 1];

              if (type == 1) {
                nbe++;
              }

              if (type == 2) {
                nt++;
              }

              if (type == 4) {
                if (verbosity > 1) cerr << "We are loading a two dimensional mesh " << endl;
                exit(1);
              }
            }

            delete[] data;
            numElementsPartial += numElms;
          }
        }

        break;
      }
    }
  }

  fclose(fp);

  Mesh::Triangle *tff = new Mesh::Triangle[nt];
  Mesh::Triangle *ttff = tff;
  Mesh::BorderElement *bff = new Mesh::BorderElement[nbe];
  Mesh::BorderElement *bbff = bff;

  // second reading
  fp = fopen(filename.c_str( ), "rb");
  if (fp == nullptr) cerr << "fopen error" << endl;
  while (!feof(fp)) {
    res = fgets(str, sizeof(str), fp);
    if (verbosity > 10 && res == nullptr) cout << "fgets error or EOF" << endl;
    if (str[0] == '$') {
      if (!strncmp(&str[1], "ELM", 3) || !strncmp(&str[1], "Elements", 8)) {
        if (!fgets(str, sizeof(str), fp)) {
          exit(1);
        }

        int numElements;
        sscanf(str, "%d", &numElements);

        if (!binary) {
          int ie = 0;
          int it = 0;

          for (int i = 0; i < numElements; i++) {
            int num, type, physical = 0, elementary = 0, partition = 0, numVertices;
            if (version <= 1.0) {
              ret = fscanf(fp, "%d %d %d %d %d", &num, &type, &physical, &elementary, &numVertices);
              if (ret == EOF) cerr << "fscanf error" << endl;
            } else {
              int numTags;
              ret = fscanf(fp, "%d %d %d", &num, &type, &numTags);
              if (ret == EOF) cerr << "fscanf error" << endl;

              for (int j = 0; j < numTags; j++) {
                int tag;
                ret = fscanf(fp, "%d", &tag);
                if (ret == EOF) cerr << "fscanf error" << endl;
                if (j == 0) {
                  physical = tag;
                } else if (j == 1) {
                  elementary = tag;
                } else if (j == 2) {
                  partition = tag;
                }

                // ignore any other tags for now
              }

              assert(type >= 1 && type <= 31);
              if ((numVertices = nvElemGmsh[type - 1]) == 0) {
                cerr << "Element of type " << type << " is not considered in FreeFem++" << endl;
                exit(1);
              }
            }

            int indices[60];

            for (int j = 0; j < numVertices; j++) {
              ret = fscanf(fp, "%d", &indices[j]);
              if (ret == EOF) cerr << "fscanf error" << endl;
            }

            if (type == 1) {
              int iv0, iv1;
              iv0 = mapnumv[indices[0]];
              iv1 = mapnumv[indices[1]];
              if (verbosity > 2) {
                cout << "Elem " << ie + 1 << " " << iv0 + 1 << " " << iv1 + 1 << endl;
              }

              (bbff++)->set(vff, iv0, iv1, physical);
              ie++;
            }

            if (type == 2) {
              int iv0, iv1, iv2;
              iv0 = mapnumv[indices[0]];
              iv1 = mapnumv[indices[1]];
              iv2 = mapnumv[indices[2]];
              if (verbosity > 2) {
                cout << "Triangles " << it + 1 << " " << iv0 + 1 << " " << iv1 + 1 << " " << iv2 + 1
                     << endl;
              }

              (ttff++)->set(vff, iv0, iv1, iv2, physical);
              if (verbosity > 2) {
                cout << "mes=" << tff[it].area << endl;
              }

              it++;
            }
          }

          assert(it == nt);
          assert(ie == nbe);
        } else {
          int ie = 0;
          int it = 0;
          int numElementsPartial = 0;

          while (numElementsPartial < numElements) {
            int header[3];
            if (fread(header, sizeof(int), 3, fp) != 3) {
              exit(1);
            }

            if (swap) {
              SwapBytes((char *)header, sizeof(int), 3);
            }

            int type = header[0];
            int numElms = header[1];
            int numTags = header[2];
            int numVertices;
            assert(type >= 1 && type <= 31);
            if ((numVertices = nvElemGmsh[type - 1]) == 0) {
              cerr << "Element of type " << type << " is not considered in FreeFem++" << endl;
              exit(1);
            }

            unsigned int n = 1 + numTags + numVertices;
            int *data = new int[n];

            for (int i = 0; i < numElms; i++) {
              if (fread(data, sizeof(int), n, fp) != n) {
                exit(1);
              }

              if (swap) {
                SwapBytes((char *)data, sizeof(int), n);
              }

              int physical = (numTags > 0) ? data[4 - numTags] : 0;
              int elementary = (numTags > 1) ? data[4 - numTags + 1] : 0;
              int *indices = &data[numTags + 1];

              if (type == 1) {
                int iv0, iv1;
                iv0 = mapnumv[indices[0]];
                iv1 = mapnumv[indices[1]];
                (bbff++)->set(vff, iv0, iv1, physical);
                ie++;
              }

              if (type == 2) {
                double mes = -1;
                int iv0, iv1, iv2;
                iv0 = mapnumv[indices[0]];
                iv1 = mapnumv[indices[1]];
                iv2 = mapnumv[indices[2]];
                (ttff++)->set(vff, iv0, iv1, iv2, physical, mes);

                it++;
              }
            }

            delete[] data;
            numElementsPartial += numElms;
          }

          assert(it == nt);
          assert(ie == nbe);
        }
      } else if (!strncmp(&str[1], "NodeData", 8)) {
        if (verbosity > 1) {
          cout << " NodeData is not considered in freefem++ " << endl;
        }
      } else if (!strncmp(&str[1], "ElementData", 11) || !strncmp(&str[1], "ElementNodeData", 15)) {
        if (verbosity > 1) {
          cout << " ElementData/ElementNodeData is not considered in freefem++ " << endl;
        }
      }
    }
  }

  fclose(fp);

  Mesh *pTh = new Mesh(nv, nt, nbe, vff, tff, bff);
  R2 Pn, Px;
  pTh->BoundingBox(Pn, Px);
  if (!pTh->quadtree) {
    pTh->quadtree = new Fem2D::FQuadTree(pTh, Pn, Px, pTh->nv);
  }
  return pTh;
}

AnyType GMSH_LoadMesh_Op::operator( )(Stack stack) const {
  string *pffname = GetAny< string * >((*filename)(stack));
  int renumsurf = 0;

  if (nargs[1]) {
    renumsurf = GetAny< long >((*nargs[1])(stack));
  }

  assert(renumsurf <= 1 && renumsurf >= 0);

  Mesh *Th = GMSH_Load(*pffname);

  Add2StackOfPtr2FreeRC(stack, Th);

  return Th;
}

// Load three dimensionnal mesh
// Mesh3

class GMSH_LoadMesh3_Op : public E_F0mps {
 public:
  Expression filename;
  static const int n_name_param = 5;    //
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }
  bool arg(int i, Stack stack, bool a) const {
    return nargs[i] ? GetAny< bool >((*nargs[i])(stack)) : a;
  }
  double arg(int i, Stack stack, double a) const {
    return nargs[i] ? GetAny< double >((*nargs[i])(stack)) : a;
  }

 public:
  GMSH_LoadMesh3_Op(const basicAC_F0 &args, Expression ffname) : filename(ffname) {
    if (verbosity > 1) {
      cout << "Load mesh given by GMSH " << endl;
    }

    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type GMSH_LoadMesh3_Op::name_param[] = {{"reftet", &typeid(long)},
                                                             {"renum", &typeid(long)},
                                                             {"cleanmesh", &typeid(bool)},
                                                             {"removeduplicate", &typeid(bool)},
                                                             {"precisvertice", &typeid(double)}

};

class GMSH_LoadMesh3 : public OneOperator {
 public:
  GMSH_LoadMesh3( ) : OneOperator(atype< pmesh3 >( ), atype< string * >( )) {}

  E_F0 *code(const basicAC_F0 &args) const {
    return new GMSH_LoadMesh3_Op(args, t[0]->CastTo(args[0]));
  }
};

Mesh3 *GMSH_Load3(const string &filename, bool cleanmesh, bool removeduplicate,
                  double precisvertice) {
  // variable freefem++
  int nv, nt = 0, nbe = 0, ret;
  Vertex3 *vff;

  map< int, int > mapnumv;

  // loading mesh and reading mesh in gmsh are in the file GModelIO_Mesh.cpp (directory Geo)
  char str[256] = "ZZZ", *res;
  double version = 2.0;
  bool binary = false, swap = false, postpro = false;
  FILE *fp = fopen(filename.c_str( ), "rb");

  if (!fp) {
    cerr << "Unable to open file " << filename.c_str( ) << endl;
    exit(1);
  }

  while (!feof(fp)) {
    res = fgets(str, sizeof(str), fp);
    if (verbosity > 10 && res == nullptr) cout << "fgets error or EOF" << endl;
    if (str[0] == '$') {
      if (!strncmp(&str[1], "MeshFormat", 10)) {
        if (!fgets(str, sizeof(str), fp)) {
          exit(1);
        }

        int format, size;
        if (sscanf(str, "%lf %d %d", &version, &format, &size) != 3) {
          exit(1);
        }

        if (format) {
          binary = true;
          if (verbosity > 1) {
            cout << "Mesh is in binary format" << endl;
          }

          int one;
          if (fread(&one, sizeof(int), 1, fp) != 1) {
            exit(1);
          }

          if (one != 1) {
            swap = true;
            if (verbosity > 1) {
              cout << "Swapping bytes from binary file" << endl;
            }
          }
        }
      } else if (!strncmp(&str[1], "PhysicalNames", 13)) {
        if (verbosity > 1) {
          cout << " PhysicalNames is not considered in freefem++ " << endl;
        }
      } else if (!strncmp(&str[1], "NO", 2) || !strncmp(&str[1], "Nodes", 5) ||
                 !strncmp(&str[1], "ParametricNodes", 15)) {
        const bool parametric = !strncmp(&str[1], "ParametricNodes", 15);
        if (parametric == true) {
          cerr << " ParametricNodes is not considered yet in FreeFem++" << endl;
          exit(1);
        }

        if (!fgets(str, sizeof(str), fp)) {
          exit(1);
        }

        if (sscanf(str, "%d", &nv) != 1) {
          exit(1);
        }

        if (verbosity > 1) {
          printf("%d vertices\n", nv);
        }

        // local variables freefem++
        vff = new Vertex3[nv];
        // map<int,int> mapnumv;

        for (int i = 0; i < nv; i++) {
          int num;
          double xyz[3], uv[2];

          // if (!parametric){
          if (!binary) {
            if (fscanf(fp, "%d %lf %lf %lf", &num, &xyz[0], &xyz[1], &xyz[2]) != 4) {
              exit(1);
            }
          } else {
            if (fread(&num, sizeof(int), 1, fp) != 1) {
              exit(1);
            }

            if (swap) {
              SwapBytes((char *)&num, sizeof(int), 1);
            }

            if (fread(xyz, sizeof(double), 3, fp) != 3) {
              exit(1);
            }

            if (swap) {
              SwapBytes((char *)xyz, sizeof(double), 3);
            }
          }

          vff[i].x = xyz[0];
          vff[i].y = xyz[1];
          vff[i].z = xyz[2];
          vff[i].lab = 1;
          mapnumv[num] = i;
        }
      } else if (!strncmp(&str[1], "ELM", 3) || !strncmp(&str[1], "Elements", 8)) {
        if (!fgets(str, sizeof(str), fp)) {
          exit(1);
        }

        int numElements;
        sscanf(str, "%d", &numElements);

        if (!binary) {
          for (int i = 0; i < numElements; i++) {
            int num, type, physical = 0, elementary = 0, partition = 0, numVertices;
            if (version <= 1.0) {
              ret = fscanf(fp, "%d %d %d %d %d", &num, &type, &physical, &elementary, &numVertices);
              if (ret == EOF) cerr << "fscanf error" << endl;
            } else {
              int numTags;
              ret = fscanf(fp, "%d %d %d", &num, &type, &numTags);
              if (ret == EOF) cerr << "fscanf error" << endl;

              for (int j = 0; j < numTags; j++) {
                int tag;
                ret = fscanf(fp, "%d", &tag);
                if (ret == EOF) cerr << "fscanf error" << endl;
                if (j == 0) {
                  physical = tag;
                } else if (j == 1) {
                  elementary = tag;
                } else if (j == 2) {
                  partition = tag;
                }

                // ignore any other tags for now
              }

              assert(type >= 1 && type <= 31);
              if ((numVertices = nvElemGmsh[type - 1]) == 0) {
                cerr << "Element of type " << type << " is not considered in FreeFem++" << endl;
                exit(1);
              }
            }

            if (type == 1) {
              if (i == 0) {
                if (verbosity > 0) {
                  cout << "edges in 3D mesh are not considered yet in freefem++, skeep data"
                       << endl;
                }
              }
            }

            if (type == 2) {
              nbe++;
            }

            if (type == 4) {
              nt++;
            }

            int indices[60];

            for (int j = 0; j < numVertices; j++) {
              ret = fscanf(fp, "%d", &indices[j]);
              if (ret == EOF) cerr << "fscanf error" << endl;
            }
          }
        } else {
          int numElementsPartial = 0;

          while (numElementsPartial < numElements) {
            int header[3];
            if (fread(header, sizeof(int), 3, fp) != 3) {
              exit(1);
            }

            if (swap) {
              SwapBytes((char *)header, sizeof(int), 3);
            }

            int type = header[0];
            int numElms = header[1];
            int numTags = header[2];
            int numVertices;
            assert(type >= 1 && type <= 31);
            if ((numVertices = nvElemGmsh[type - 1]) == 0) {
              cerr << "Element of type " << type << " is not considered in FreeFem++" << endl;
              exit(1);
            }

            unsigned int n = 1 + numTags + numVertices;
            int *data = new int[n];

            for (int i = 0; i < numElms; i++) {
              if (fread(data, sizeof(int), n, fp) != n) {
                exit(1);
              }

              if (swap) {
                SwapBytes((char *)data, sizeof(int), n);
              }

              int num = data[0];
              int physical = (numTags > 0) ? data[4 - numTags] : 0;
              int elementary = (numTags > 1) ? data[4 - numTags + 1] : 0;
              int partition = (numTags > 2) ? data[4 - numTags + 2] : 0;
              int *indices = &data[numTags + 1];

              if (type == 1 && i == 0) {
                if (verbosity>0) cout << " Edges in 3D mesh are not used in FreeFem++; skip data" << endl;
                // exit(1);
              }

              if (type == 2) {
                nbe++;
              }

              if (type == 4) {
                nt++;
              }
            }

            delete[] data;
            numElementsPartial += numElms;
          }
        }

        break;
      }
    }
  }

  fclose(fp);

  if (verbosity > 1) {
    cout << "closing file " << nt << " " << nbe << endl;
  }

  Tet *tff = new Tet[nt];
  Tet *ttff = tff;
  Triangle3 *bff = new Triangle3[nbe];
  Triangle3 *bbff = bff;

  // second reading
  fp = fopen(filename.c_str( ), "rb");

  while (!feof(fp)) {
    res = fgets(str, sizeof(str), fp);
    if (verbosity > 10 && res == nullptr) cout << "fgets error or EOF" << endl;
    if (str[0] == '$') {
      if (!strncmp(&str[1], "ELM", 3) || !strncmp(&str[1], "Elements", 8)) {
        if (!fgets(str, sizeof(str), fp)) {
          exit(1);
        }

        int numElements;
        sscanf(str, "%d", &numElements);

        if (verbosity > 0) {
          printf("%d tetrahedrons\n", nt);
          printf("%d triangles\n", nbe);
          printf("%d numElements\n", numElements);
        }

        if (!binary) {
          int ie = 0;
          int it = 0;

          for (int i = 0; i < numElements; i++) {
            int num, type, physical = 0, elementary = 0, partition = 0, numVertices;
            if (version <= 1.0) {
              ret = fscanf(fp, "%d %d %d %d %d", &num, &type, &physical, &elementary, &numVertices);
              if (ret == EOF) cerr << "fscanf error" << endl;
            } else {
              int numTags;
              ret = fscanf(fp, "%d %d %d", &num, &type, &numTags);
              if (ret == EOF) cerr << "fscanf error" << endl;

              for (int j = 0; j < numTags; j++) {
                int tag;
                ret = fscanf(fp, "%d", &tag);
                if (ret == EOF) cerr << "fscanf error" << endl;

                if (j == 0) {
                  physical = tag;
                } else if (j == 1) {
                  elementary = tag;
                } else if (j == 2) {
                  partition = tag;
                }

                // ignore any other tags for now
              }

              assert(type >= 1 && type <= 31);
              if ((numVertices = nvElemGmsh[type - 1]) == 0) {
                cerr << "Element of type " << type << " is not considered in FreeFem++" << endl;
                exit(1);
              }
            }

            int indices[60];

            for (int j = 0; j < numVertices; j++) {
              ret = fscanf(fp, "%d", &indices[j]);
              if (ret == EOF) cerr << "fscanf error" << endl;
            }

            if (type == 2) {
              int ivff[3];

              for (int ii = 0; ii < numVertices; ii++) {
                ivff[ii] = mapnumv[indices[ii]];
                assert(ivff[ii] >= 0 && ivff[ii] < nv);
              }

              (bbff++)->set(vff, ivff, physical);
              ie++;
            }

            if (type == 4) {
              int ivff[4];

              for (int ii = 0; ii < numVertices; ii++) {
                ivff[ii] = mapnumv[indices[ii]];
                assert(ivff[ii] >= 0 && ivff[ii] < nv);
              }

              (ttff++)->set(vff, ivff, physical);
              it++;
            }
          }

          assert(it == nt);
          assert(ie == nbe);
        } else {
          int ie = 0;
          int it = 0;
          int numElementsPartial = 0;

          while (numElementsPartial < numElements) {
            int header[3] = {};
            if (fread(header, sizeof(int), 3, fp) != 3) {
              exit(1);
            }

            if (swap) {
              SwapBytes((char *)header, sizeof(int), 3);
            }

            int type = header[0];
            int numElms = header[1];
            int numTags = header[2];
            int numVertices;
            assert(type >= 1 && type <= 31);
            if ((numVertices = nvElemGmsh[type - 1]) == 0) {
              cerr << "Element of type " << type << " is not considered in FreeFem++" << endl;
              exit(1);
            }

            unsigned int n = 1 + numTags + numVertices;
            int *data = new int[n];

            for (int i = 0; i < numElms; i++) {
              if (fread(data, sizeof(int), n, fp) != n) {
                exit(1);
              }

              if (swap) {
                SwapBytes((char *)data, sizeof(int), n);
              }

              int num = data[0];
              int physical = (numTags > 0) ? data[4 - numTags] : 0;
              int elementary = (numTags > 1) ? data[4 - numTags + 1] : 0;
              int partition = (numTags > 2) ? data[4 - numTags + 2] : 0;
              int *indices = &data[numTags + 1];

              if (type == 2) {
                int ivff[3];

                for (int ii = 0; ii < numVertices; ii++) {
                  ivff[ii] = mapnumv[indices[ii]];
                }

                (bbff++)->set(vff, ivff, physical);
                ie++;
              }

              if (type == 4) {
                int ivff[4];

                for (int ii = 0; ii < numVertices; ii++) {
                  ivff[ii] = mapnumv[indices[ii]];
                }

                (ttff++)->set(vff, ivff, physical);
                it++;
              }
            }

            delete[] data;
            numElementsPartial += numElms;
          }

          assert(it == nt);
          assert(ie == nbe);
        }
      } else if (!strncmp(&str[1], "NodeData", 8)) {
        if (verbosity) {
          cout << " NodeData is not considered in FreeFem++ " << endl;
        }
      } else if (!strncmp(&str[1], "ElementData", 11) || !strncmp(&str[1], "ElementNodeData", 15)) {
        if (verbosity) {
          cout << " ElementData/ElementNodeData is not considered in FreeFem++ " << endl;
        }
      }
    }
  }

  fclose(fp);

  if (nt == 0) {
    cout << " Return type false. Your mesh is MeshS type; use gmshloadS function." << endl;
    ffassert(0);
    
  } else {
    Mesh3 *Th3 = new Mesh3(nv, nt, nbe, vff, tff, bff, cleanmesh|| (nbe==0), removeduplicate,(nbe==0), precisvertice);
    return Th3;
  }
}

AnyType GMSH_LoadMesh3_Op::operator( )(Stack stack) const {
  string *pffname = GetAny< string * >((*filename)(stack));
  int renumsurf = 0;

  if (nargs[1]) renumsurf = GetAny< long >((*nargs[1])(stack));

  bool cleanmesh(arg(2, stack, false));
  bool removeduplicate(arg(3, stack, false));
  double precisvertice(arg(4, stack, 1e-6));

  assert(renumsurf <= 1 && renumsurf >= 0);

  Mesh3 *Th3_t = GMSH_Load3(*pffname, cleanmesh, removeduplicate, precisvertice);

  // Th3_t->getTypeMesh3() = 1;
  Th3_t->BuildGTree( );
  Add2StackOfPtr2FreeRC(stack, Th3_t);

  return Th3_t;
}


template <class MMesh>
class GMSH_LoadMeshT_Op : public E_F0mps {
 public:
  Expression filename;
  static const int n_name_param = 6;    //
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }
  bool arg(int i, Stack stack, bool a) const {
    return nargs[i] ? GetAny< bool >((*nargs[i])(stack)) : a;
  }
  double arg(int i, Stack stack, double a) const {
    return nargs[i] ? GetAny< double >((*nargs[i])(stack)) : a;
  }

 public:
  GMSH_LoadMeshT_Op(const basicAC_F0 &args, Expression ffname) : filename(ffname) {
    if (verbosity > 1) {
      cout << "Load mesh given by GMSH." << endl;
    }

    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator( )(Stack stack) const;
};

template < >
basicAC_F0::name_and_type GMSH_LoadMeshT_Op<MeshS>::name_param[] = {
	{"reftri", &typeid(long)},
    {"renum", &typeid(long)},
	{"cleanmesh", &typeid(bool)},
	{"removeduplicate", &typeid(bool)},
	{"precisvertice", &typeid(double)},
	{"ridgeangledetection", &typeid(double) }};

template < >
basicAC_F0::name_and_type GMSH_LoadMeshT_Op<MeshL>::name_param[] = {
	{"refedge", &typeid(long)},
	{"renum", &typeid(long)},
	{"cleanmesh", &typeid(bool)},
	{"removeduplicate", &typeid(bool)},
	{"precisvertice", &typeid(double)},
	{"ridgeangledetection", &typeid(double) }};



template<class MMesh>
class GMSH_LoadMeshT : public OneOperator {
 public:
  typedef const MMesh *ppmesh;
  GMSH_LoadMeshT( ) : OneOperator(atype< ppmesh >( ), atype< string * >( )) {}

  E_F0 *code(const basicAC_F0 &args) const {
    return new GMSH_LoadMeshT_Op<MMesh>(args, t[0]->CastTo(args[0]));
  }
};



template<class MMesh>
MMesh *GMSH_LoadT(const string &filename, bool cleanmesh, bool removeduplicate,
                  double precisvertice, double ridgeangledetection) {
  
  
  typedef typename MMesh::Element T;
  typedef typename MMesh::BorderElement B;
  typedef typename MMesh::Vertex V;
  typedef typename MMesh::Element::RdHat TRdHat;
  typedef typename MMesh::BorderElement::RdHat BRdHat;
					  
					  
  int nv, nt = 0, nbe = 0;
  V *vff;

  map< int, int > mapnumv;

  // loading mesh and reading mesh in gmsh are in the file GModelIO_Mesh.cpp (directory Geo)
  char str[256] = "ZZZ";
  double version = 2.0;
  bool binary = false, swap = false, postpro = false;
  FILE *fp = fopen(filename.c_str( ), "rb");
  if (!fp) {
    cerr << "Unable to open file " << filename.c_str( ) << endl;
    exit(1);
  }

  while (!feof(fp)) {
    fgets(str, sizeof(str), fp);
    if (str[0] == '$') {
      if (!strncmp(&str[1], "MeshFormat", 10)) {
        if (!fgets(str, sizeof(str), fp)) {
          exit(1);
        }

        int format, size;
        if (sscanf(str, "%lf %d %d", &version, &format, &size) != 3) {
          exit(1);
        }

        if (format) {
          binary = true;
          if (verbosity > 1) {
            cout << "Mesh is in binary format" << endl;
          }

          int one;
          if (fread(&one, sizeof(int), 1, fp) != 1) {
            exit(1);
          }

          if (one != 1) {
            swap = true;
            if (verbosity > 1) {
              cout << "Swapping bytes from binary file" << endl;
            }
          }
        }
      } else if (!strncmp(&str[1], "PhysicalNames", 13)) {
        if (verbosity > 1) {
          cout << " PhysicalNames is not considered in freefem++ " << endl;
        }
      } else if (!strncmp(&str[1], "NO", 2) || !strncmp(&str[1], "Nodes", 5) ||
                 !strncmp(&str[1], "ParametricNodes", 15)) {
        const bool parametric = !strncmp(&str[1], "ParametricNodes", 15);
        if (parametric == true) {
          cerr << " ParametricNodes is not considered yet in FreeFem++" << endl;
          exit(1);
        }

        if (!fgets(str, sizeof(str), fp)) {
          exit(1);
        }

        if (sscanf(str, "%d", &nv) != 1) {
          exit(1);
        }

        if (verbosity > 1) {
          printf("%d vertices\n", nv);
        }

        // local variables freefem++
        vff = new V[nv];
        // map<int,int> mapnumv;

        int minVertex = nv + 1, maxVertex = -1;

        for (int i = 0; i < nv; i++) {
          int num;
          double xyz[3], uv[2];

          // if (!parametric){
          if (!binary) {
            if (fscanf(fp, "%d %lf %lf %lf", &num, &xyz[0], &xyz[1], &xyz[2]) != 4) {
              exit(1);
            }
          } else {
            if (fread(&num, sizeof(int), 1, fp) != 1) {
              exit(1);
            }

            if (swap) {
              SwapBytes((char *)&num, sizeof(int), 1);
            }

            if (fread(xyz, sizeof(double), 3, fp) != 3) {
              exit(1);
            }

            if (swap) {
              SwapBytes((char *)xyz, sizeof(double), 3);
            }
          }

          vff[i].x = xyz[0];
          vff[i].y = xyz[1];
          vff[i].z = xyz[2];
          vff[i].lab = 1;
          mapnumv[num] = i;
        }
      } else if (!strncmp(&str[1], "ELM", 3) || !strncmp(&str[1], "Elements", 8)) {
        if (!fgets(str, sizeof(str), fp)) {
          exit(1);
        }

        int numElements;
        sscanf(str, "%d", &numElements);

        if (!binary) {
          for (int i = 0; i < numElements; i++) {
            int num, type, physical = 0, elementary = 0, partition = 0, numVertices;
            if (version <= 1.0) {
              fscanf(fp, "%d %d %d %d %d", &num, &type, &physical, &elementary, &numVertices);
            } else {
              int numTags;
              fscanf(fp, "%d %d %d", &num, &type, &numTags);

              for (int j = 0; j < numTags; j++) {
                int tag;
                fscanf(fp, "%d", &tag);
                if (j == 0) {
                  physical = tag;
                } else if (j == 1) {
                  elementary = tag;
                } else if (j == 2) {
                  partition = tag;
                }

                // ignore any other tags for now
              }

              assert(type >= 1 && type <= 31);
              if ((numVertices = nvElemGmsh[type - 1]) == 0) {
                cerr << "Element of type " << type << ", nv elem =" << nvElemGmsh[type - 1] << " is not considered in FreeFem++" << endl;
                exit(1);
              }
            }

            if(type == 1 && is_same< MMesh, MeshS >::value) nbe++;
			else if(type == 1 && is_same< MMesh, MeshL >::value) nt++;
            
            if (type == 2 && is_same< MMesh, MeshS >::value) nt++;

            if (type == 4) {
              cerr << "We are loading a three dimensional SURFACE mesh." << endl;
              exit(1);
            }

            int indices[60];

            for (int j = 0; j < numVertices; j++) {
              fscanf(fp, "%d", &indices[j]);
            }
          }
        } else {
          int numElementsPartial = 0;

          while (numElementsPartial < numElements) {
            int header[3];
            if (fread(header, sizeof(int), 3, fp) != 3) {
              exit(1);
            }

            if (swap) {
              SwapBytes((char *)header, sizeof(int), 3);
            }

            int type = header[0];
            int numElms = header[1];
            int numTags = header[2];
            int numVertices;
            assert(type >= 1 && type <= 31);
            if ((numVertices = nvElemGmsh[type - 1]) == 0) {
              cerr << "Element of type " << type << " is not considered in FreeFem++" << endl;
              exit(1);
            }

            unsigned int n = 1 + numTags + numVertices;
            int *data = new int[n];

            for (int i = 0; i < numElms; i++) {
              if (fread(data, sizeof(int), n, fp) != n) {
                exit(1);
              }

              if (swap) {
                SwapBytes((char *)data, sizeof(int), n);
              }

              int num = data[0];
              int physical = (numTags > 0) ? data[4 - numTags] : 0;
              int elementary = (numTags > 1) ? data[4 - numTags + 1] : 0;
              int partition = (numTags > 2) ? data[4 - numTags + 2] : 0;
              int *indices = &data[numTags + 1];

              if(type == 1 && is_same< MMesh, MeshS >::value) nbe++;
  			  else if(type == 1 && is_same< MMesh, MeshL >::value) nt++;
            
              if (type == 2 && is_same< MMesh, MeshS >::value) nt++;

              if (type == 4) {
                cerr << "We are loading a three dimensional SURFACE mesh " << endl;
                exit(1);
              }
            }

            delete[] data;
            numElementsPartial += numElms;
          }
        }

        break;
      }
    }
  }

  fclose(fp);

  if (verbosity > 1) {
    cout << "closing file " << nt << " " << nbe << endl;
  }

  T *tff = new T[nt];
  T *ttff = tff;
  B *bff = new B[nbe];
  B *bbff = bff;

  // second reading
  fp = fopen(filename.c_str( ), "rb");

  while (!feof(fp)) {
    fgets(str, sizeof(str), fp);
    if (str[0] == '$') {
      if (!strncmp(&str[1], "ELM", 3) || !strncmp(&str[1], "Elements", 8)) {
        if (!fgets(str, sizeof(str), fp)) {
          exit(1);
        }

        int numElements;
        sscanf(str, "%d", &numElements);

        if (verbosity > 0) {
          printf("%d elements\n", nt);
          printf("%d border elements\n", nbe);
          printf("%d numElements\n", numElements);
        }

        if (!binary) {
          int ie = 0;
          int it = 0;

          for (int i = 0; i < numElements; i++) {
            int num, type, physical = 0, elementary = 0, partition = 0, numVertices;
            if (version <= 1.0) {
              fscanf(fp, "%d %d %d %d %d", &num, &type, &physical, &elementary, &numVertices);
            } else {
              int numTags;
              fscanf(fp, "%d %d %d", &num, &type, &numTags);

              for (int j = 0; j < numTags; j++) {
                int tag;
                fscanf(fp, "%d", &tag);

                if (j == 0) {
                  physical = tag;
                } else if (j == 1) {
                  elementary = tag;
                } else if (j == 2) {
                  partition = tag;
                }

                // ignore any other tags for now
              }

              assert(type >= 1 && type <= 31);
              if ((numVertices = nvElemGmsh[type - 1]) == 0) {
                cerr << "Element of type " << type << " is not considered in FreeFem++" << endl;
                exit(1);
              }
            }

            int indices[60];

            for (int j = 0; j < numVertices; j++) {
              fscanf(fp, "%d", &indices[j]);
            }

            if (type == 1) {
              int iv[2];
              for (int i = 0; i < 2; i++) iv[i] = mapnumv[indices[i]];
              if (verbosity > 2) {
                cout << "Elem " << ie + 1 << " " << iv[0] + 1 << " " << iv[1] + 1 << endl;
              }

              if(is_same< MMesh, MeshS >::value) {
			  	(bbff++)->set(vff, iv, physical);
              	ie++;
			  }
              else if(is_same< MMesh, MeshL >::value) {
			  	(ttff++)->set(vff, iv, physical);
              	it++;
			  }
			
            }

            if (type == 2 && is_same< MMesh, MeshS >::value) {
              int iv[3];
              for (int i = 0; i < 3; i++) iv[i] = mapnumv[indices[i]];
              if (verbosity > 2) {
                cout << "Triangles " << it + 1 << " " << iv[0] + 1 << " " << iv[1] + 1 << " "
                     << iv[2] + 1 << endl;
              }

              (ttff++)->set(vff, iv, physical);
              if (verbosity > 2) {
                cout << "mes=" << tff[it].mesure( ) << endl;
              }

              if (tff[it].mesure( ) < 1e-8) {
                cerr << "bug : measure < 1e-8 !" << endl;
                exit(1);
              }

              it++;
            }
          }

          assert(it == nt);
          assert(ie == nbe);
        } else {
          int ie = 0;
          int it = 0;
          int numElementsPartial = 0;

          while (numElementsPartial < numElements) {
            int header[3];
            if (fread(header, sizeof(int), 3, fp) != 3) {
              exit(1);
            }

            if (swap) {
              SwapBytes((char *)header, sizeof(int), 3);
            }

            int type = header[0];
            int numElms = header[1];
            int numTags = header[2];
            int numVertices;
            assert(type >= 1 && type <= 31);
            if ((numVertices = nvElemGmsh[type - 1]) == 0) {
              cerr << "Element of type " << type << " is not considered in Freefem++" << endl;
              exit(1);
            }

            unsigned int n = 1 + numTags + numVertices;
            int *data = new int[n];

            for (int i = 0; i < numElms; i++) {
              if (fread(data, sizeof(int), n, fp) != n) {
                exit(1);
              }

              if (swap) 
                SwapBytes((char *)data, sizeof(int), n);
              

              int num = data[0];
              int physical = (numTags > 0) ? data[4 - numTags] : 0;
              int elementary = (numTags > 1) ? data[4 - numTags + 1] : 0;
              int partition = (numTags > 2) ? data[4 - numTags + 2] : 0;
              int *indices = &data[numTags + 1];

              if (type == 1) {
                int iv[2];
                for (int i = 0; i < 2; i++) iv[i] = mapnumv[indices[i]];
                if(is_same< MMesh, MeshS >::value) {
					(bbff++)->set(vff, iv, physical);
                    ie++;
				}
				else if(is_same< MMesh, MeshL >::value) {
					double mes = -1;
	                (ttff++)->set(vff, iv, physical, mes);
					it++;	
				}
              }

              if (type == 2) {
                double mes = -1;
                int iv[3];
                for (int i = 0; i < 3; i++) iv[i] = mapnumv[indices[i]];
				if(is_same< MMesh, MeshS >::value) {
                	(ttff++)->set(vff, iv, physical, mes);
					it++;
				}
			}
            }

            delete[] data;
            numElementsPartial += numElms;
          }

          assert(it == nt);
          assert(ie == nbe);
        }
      } else if (!strncmp(&str[1], "NodeData", 8)) {
        if (verbosity) {
          cout << " NodeData is not considered in FreeFem++ " << endl;
        }
      } else if (!strncmp(&str[1], "ElementData", 11) || !strncmp(&str[1], "ElementNodeData", 15)) {
        if (verbosity) {
          cout << " ElementData/ElementNodeData is not considered in FreeFem++ " << endl;
        }
      }
    }
  }

  fclose(fp);

  MMesh *Th = new MMesh(nv, nt, nbe, vff, tff, bff, cleanmesh, removeduplicate, precisvertice, ridgeangledetection);
  return Th;
}

template<class MMesh>
AnyType GMSH_LoadMeshT_Op<MMesh>::operator( )(Stack stack) const {
  string *pffname = GetAny< string * >((*filename)(stack));
  int renumsurf = 0;

  if (nargs[1]) renumsurf = GetAny< long >((*nargs[1])(stack));
  bool cleanmesh(arg(2, stack, false));
  bool removeduplicate(arg(3, stack, false));
  double precisvertice(arg(4, stack, 1e-6));
  double ridgeangledetection(arg(5, stack, 8.*atan(1.)/9.));

  assert(renumsurf <= 1 && renumsurf >= 0);

  MMesh *Th_t = GMSH_LoadT<MMesh>(*pffname, cleanmesh, removeduplicate, precisvertice, ridgeangledetection);

  Th_t->BuildGTree( );
  Add2StackOfPtr2FreeRC(stack, Th_t);

  return Th_t;
}

/*  class Init1 { public:
 * Init1();
 * };
 *
 * $1 */

bool SaveGMSH(pmesh3 pTh, string *filewoext) {
  /*  Code from Lo Sala <salalo80@gmail.com>
   * May 3th 2016 ...
   * /This function save a
   * func int SaveGmsh(mesh3 & msh, string namewoextension)
   * {
   * string nameoffile=namewoextension+".msh";
   * //mesh3 msh=amsh;//To jump over a "bug"
   * ofstream f1(nameoffile);
   * //f1.scientific;
   * f1.precision(6);
   * int nbvertices=msh.nv;
   * f1<<"$MeshFormat"<<endl;
   * f1<<"2.2 0 8"<<endl;
   * f1<<"$EndMeshFormat"<<endl;
   * f1<<"$Nodes"<<endl;
   * f1<<nbvertices<<endl;
   * for(int i=0;i<nbvertices;++i)
   * {
   * f1<<(i+1)<<" "<<msh(i).x<<" "<<msh(i).y<<" "<<msh(i).z<<endl;
   * }
   * f1<<"$EndNodes"<<endl;
   * f1<<"$Elements"<<endl;
   * f1<<msh.nt+msh.nbe<<endl;
   * for(int i=0;i<msh.nbe;++i)
   * {
   * //2 is a triangle
   * f1<<(i+1)<<" 2 ";
   * //two tags: the label
   * f1<<"2 "<<msh.be(i).label<<" "<<msh.be(i).label<<" ";
   * //list of nodes
   * f1<<msh.be(i)[0]+1<<" "<<msh.be(i)[1]+1<<" "<<msh.be(i)[2]+1<<endl;
   * }
   *
   * for(int i=0;i<msh.nt;++i)
   * {
   * //4 is a tethrahedron
   * f1<<(msh.nbe+i+1)<<" 4 ";
   * //two tags: the label
   * f1<<"2 "<<msh[i].label<<" "<<msh[i].label<<" ";
   * //list of nodes
   * f1<<msh[i][0]+1<<" "<<msh[i][1]+1<<" "<<msh[i][2]+1<<" "<<msh[i][3]+1<<endl;
   * }
   *
   * f1<<"$EndElements"<<endl;
   *
   * return 0;
   * }
   *
   */
  string file = *filewoext + ".msh";
  ofstream f1(file.c_str( ));

  if (!f1) {
    ffassert(f1);
    return 1;
  }

  f1.precision(15);
  const Mesh3 &msh = *pTh;
  long nbvertices = msh.nv;
  f1 << "$MeshFormat" << endl;
  f1 << "2.2 0 8" << endl;
  f1 << "$EndMeshFormat" << endl;
  f1 << "$Nodes" << endl;
  f1 << nbvertices << endl;

  for (int i = 0; i < nbvertices; ++i) {
    f1 << (i + 1) << " " << msh(i).x << " " << msh(i).y << " " << msh(i).z << endl;
  }

  f1 << "$EndNodes" << endl;
  f1 << "$Elements" << endl;
  f1 << msh.nt + msh.nbe << endl;

  for (int i = 0; i < msh.nbe; ++i) {
    // 2 is a triangle
    f1 << (i + 1) << " 2 ";
    // two tags: the label
    f1 << "2 " << msh.be(i).lab << " " << msh.be(i).lab << " ";
    // list of nodes
    f1 << msh(msh.be(i)[0]) + 1 << " " << msh(msh.be(i)[1]) + 1 << " " << msh(msh.be(i)[2]) + 1
       << endl;
  }

  for (int i = 0; i < msh.nt; ++i) {
    // 4 is a tethrahedron
    f1 << (msh.nbe + i + 1) << " 4 ";
    // two tags: the label
    f1 << "2 " << msh[i].lab << " " << msh[i].lab << " ";
    // list of nodes
    f1 << msh(msh[i][0]) + 1 << " " << msh(msh[i][1]) + 1 << " " << msh(msh[i][2]) + 1 << " "
       << msh(msh[i][3]) + 1 << endl;
  }

  f1 << "$EndElements" << endl;
  return 0;    // OK ..
}

bool SaveGMSH(pmeshS pTh, string *filewoext) {

  string file = *filewoext + ".msh";
  ofstream f1(file.c_str( ));

  if (!f1) {
      cout << " Error Opening file " << file << endl;
      ExecError("Error Opening file");
    return 1;
  }

  f1.precision(15);
  const MeshS &msh = *pTh;
  long nbvertices = msh.nv;
  f1 << "$MeshFormat" << endl;
  f1 << "2.2 0 8" << endl;
  f1 << "$EndMeshFormat" << endl;
  f1 << "$Nodes" << endl;
  f1 << nbvertices << endl;

  for (int i = 0; i < nbvertices; ++i) {
    f1 << (i + 1) << " " << msh(i).x << " " << msh(i).y << " " << msh(i).z << endl;
  }

  f1 << "$EndNodes" << endl;
  f1 << "$Elements" << endl;
  f1 << msh.nt + msh.nbe << endl;

  for (int i = 0; i < msh.nbe; ++i) {
    // 1 is an edge
    f1 << (i + 1) << " 1 ";
    // two tags: the label
    f1 << "1 " << msh.be(i).lab << " ";    //<< msh.be(i).lab << " ";
    // list of nodes
    f1 << msh(msh.be(i)[0]) + 1 << " " << msh(msh.be(i)[1]) + 1 << endl;
  }

  for (int i = 0; i < msh.nt; ++i) {
    // 2 is a triangle
    f1 << (msh.nbe + i + 1) << " 2 ";
    // two tags: the label
    f1 << "2 " << msh[i].lab << " " << msh[i].lab << " ";
    // list of nodes
    f1 << msh(msh[i][0]) + 1 << " " << msh(msh[i][1]) + 1 << " " << msh(msh[i][2]) + 1 << endl;
  }

  f1 << "$EndElements" << endl;
  return 0;    // OK ..
}

bool SaveGMSH(pmeshL pTh, string *filewoext) {

  string file = *filewoext + ".msh";
  ofstream f1(file.c_str( ));

  if (!f1) {
    ffassert(f1);
    return 1;
  }

  f1.precision(15);
  const MeshL &msh = *pTh;
  long nbvertices = msh.nv;
  f1 << "$MeshFormat" << endl;
  f1 << "2.2 0 8" << endl;
  f1 << "$EndMeshFormat" << endl;
  f1 << "$Nodes" << endl;
  f1 << nbvertices << endl;

  for (int i = 0; i < nbvertices; ++i) {
    f1 << (i + 1) << " " << msh(i).x << " " << msh(i).y << " " << msh(i).z << endl;
  }

  f1 << "$EndNodes" << endl;
  f1 << "$Elements" << endl;
  f1 << msh.nt << endl;


  for (int i = 0; i < msh.nt; ++i) {
    // 1 is an edge
    f1 << (i + 1) << " 1 ";
    // two tags: the label
    f1 << "1 " << msh[i].lab << " ";
    // list of nodes
    f1 << msh(msh[i][0]) + 1 << " " << msh(msh[i][1]) + 1  << endl;
  }

  f1 << "$EndElements" << endl;
  return 0;    // OK ..
}





static void Load_Init( ) {    // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
  // if (verbosity)
  if (verbosity > 1 && (mpirank == 0)) {
    cout << " load: gmsh " << endl;
  }

  Global.Add("gmshload3", "(", new GMSH_LoadMesh3);
  Global.Add("gmshloadS", "(", new GMSH_LoadMeshT<MeshS>);
  Global.Add("gmshloadL", "(", new GMSH_LoadMeshT<MeshL>);
  Global.Add("gmshload", "(", new GMSH_LoadMesh);
  Global.Add("savegmsh", "(", new OneOperator2< bool, pmesh3, string * >(SaveGMSH));
  Global.Add("savegmsh", "(", new OneOperator2< bool, pmeshS, string * >(SaveGMSH));
  Global.Add("savegmsh", "(", new OneOperator2< bool, pmeshL, string * >(SaveGMSH));
  
}

LOADFUNC(Load_Init)
