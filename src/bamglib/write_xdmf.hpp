#ifndef write_xdmf_h_
#define write_xdmf_h_

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
using std::string;
using std::vector;
typedef long Int4;

class WriteXdmf {
 private:
  const char* WXffname;
  ofstream xdmf_file;
  char* Elemtype;
  Int4 nbofelem;
  Int4 nbofvertex;
  Int4 nbvperelem;
  char* xdmf_filename;
  int dimension;

 public:
  WriteXdmf(const char* ffname, Int4 nbelem, Int4 nbvertex);
  virtual ~WriteXdmf( );
  void WriteXdmfMeshFile2D( );
  // void WriteXdmfMeshFile3D();
  void WriteXdmfSolFile2DInit( );
  void WriteXdmfSolFile2DAddField(string* fieldname, int data_type, int result_order,
                                  int trans_dim);
  void WriteXdmfSolFile2DFinalize( );
  void WriteXdmfSolFile3DInit( );
  void WriteXdmfSolFile3DAddField(string* fieldname, int data_type, int result_order,
                                  int trans_dim);
  void WriteXdmfSolFile3DFinalize( );
};

#endif
