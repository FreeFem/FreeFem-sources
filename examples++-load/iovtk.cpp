// ORIG-DATE:   September 2009
// -*- Mode : c++ -*%
//
// SUMMARY  : READ/WRITE MESH AND SOLUTION IN FORMAT VTK    
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
// AUTHOR   : Jacques Morice
// E-MAIL   : jacques.morice@ann.jussieu.fr
//

/* 
 This file is part of Freefem++
 
 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.
 
 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

 Thank to the ARN ()  FF2A3 grant
 ref:ANR-07-CIS7-002-01 
 */

//  FH   July 2009
//   comment all
//       Th3_t->BuildBound();
//       Th3_t->BuildAdj();
//       Th3_t->Buildbnormalv();  
//       Th3_t->BuildjElementConteningVertex();
//   is now in the constructor of Mesh3 to be consistante. 
//  FH  nev 2009
//  correct  gestion of nameofuser    variable

#include  <iostream>
#include  <cfloat>
#include <cmath>
#include <complex>
using namespace std;

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
//#include "LayerMesh.hpp"
//#include "TransfoMesh_v2.hpp"
#include "msh3.hpp"
//#include "GQuadTree.hpp"
//#include "lex.hpp"
#include <set>
#include <vector>
#include <list>
#include <fstream>

using namespace Fem2D;

char * newcopy(const char * s)
{
  char *r(new char  [strlen(s)+1]);
  strcpy(r, s);
  return r;
}

// Tables of element type of VTK considered in Freefem++
static const int nvElemVTK[25] = {  1, 0, 2, 0, 3, 
				    0, 0, 0, 0, 4, 
				    0, 0, 0, 0, 0, 
				    0, 0, 0, 0, 0, 
				    0, 0, 0, 0, 0 };   

// we considerer only vertex, edges, triangles and tetrahedrons in Freefem++
//  1  :: Vertex Corner
//  3  :: Edge/line
//  5  :: triangles
// 10  :: tetrahedrons

enum FFppCells{VTK_EDGE=3, VTK_TRI=5, VTK_TET=10};

static const int NbColorTable=30;
// Table of colors for Labels of elements :: RGB
static const float ColorTable[30][3] = {  
  {1.0,0.0,0.0},  /* red    */
  {1.0,1.0,0.0},  /* yellow */
  {1.0,0.0,1.0},  /* ????   */ 
  {0.0,1.0,0.0},  /* green  */ 
  {0.0,1.0,1.0},  /* cyan   */ 
  {0.0,0.0,1.0},  /* blue   */ 
  {1.0,  0.5, 0.0},  /* orange */
  {0.5,  0.0, 1.0},  /* violet */ 
  {1.0,  0.0, 0.5},  /* ???  */
  {0.5,  1.0, 0.0},  /* ???  */
  {0.0,  1.0, 0.5},  /* ???  */
  {0.0,  0.5, 1.0},  /* ???  */
  {1.0,  0.5, 0.5},  /*  ???    */
  {1.0,  1.0, 0.5},  /*  ???    */
  {1.0,  0.5, 1.0},  /*  ???    */ 
  {0.5,  1.0, 0.5},  /*  ???    */ 
  {0.5,  1.0, 1.0},  /*  ???   */ 
  {0.5,  0.5, 1.0},  /*  ???   */ 
  {0.4,  0.0, 0.0},  /* dark blue    */
  {0.4,  0.4, 0.0},  /* dark yellow */
  {0.4,  0.0, 0.4},  /* dark ????   */ 
  {0.0,  0.4, 0.0},  /* dark green  */ 
  {0.0,  0.4, 0.4},  /* dark cyan   */ 
  {0.0,  0.0, 0.4},  /* dark blue   */
  {0.8,  0.0, 0.0},  /*  ???    */
  {0.8,  0.8, 0.0},  /*  ???    */
  {0.8,  0.0, 0.8},  /*  ???    */ 
  {0.0,  0.8, 0.0},  /*  ???    */ 
  {0.0,  0.8, 0.8},  /*  ???    */ 
  {0.0,  0.0, 0.8},  /*  ???    */
};  // a voir F.Hecht

void SwapBytes(char *array, int size, int n)
  {
    char *x = new char[size];
    for(int i = 0; i < n; i++) {
      char *a = &array[i * size];
      memcpy(x, a, size);
      for(int c = 0; c < size; c++)
        a[size - 1 - c] = x[c];
    }
    delete [] x;
  }

// two dimensional case

// LOAD fichier.vtk
class VTK_LoadMesh_Op : public E_F0mps 
{
public:
  Expression filename;
  static const int n_name_param = 4; // 
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
public:
  VTK_LoadMesh_Op(const basicAC_F0 &  args,Expression ffname) 
    : filename(ffname)
  {
    if(verbosity) cout << "Load mesh given by VTK " << endl;
    args.SetNameParam(n_name_param,name_param,nargs);
  } 
  
  AnyType operator()(Stack stack)  const ;
};
basicAC_F0::name_and_type VTK_LoadMesh_Op::name_param[]= {
  {  "reft", &typeid(long)},
  {  "swap", &typeid(bool)},
  {  "refe", &typeid(long)},
  {  "namelabel",&typeid(string)}
};


class  VTK_LoadMesh : public OneOperator { public:  
    VTK_LoadMesh() : OneOperator(atype<pmesh>(),atype<string *>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
    return  new VTK_LoadMesh_Op( args,t[0]->CastTo(args[0]) ); 
  }
};

Mesh * VTK_Load(const string & filename, bool bigEndian)
// swap = bigEndian or not bigEndian
{
  // variable freefem++
  int nv, nt=0, nbe=0;
  Mesh::Vertex   *vff;
  map<int,int> mapnumv;

  // Reading Mesh in vtk formats 

  FILE *fp = fopen(filename.c_str(), "rb");
  if(!fp){
    cerr << "Unable to open file " << filename.c_str() << endl;
    exit(1);
  }

  char buffer[256], buffer2[256];
    
  fgets(buffer, sizeof(buffer), fp); // version line
  fgets(buffer, sizeof(buffer), fp); // title
  
  fscanf(fp, "%s", buffer); // ASCII or BINARY
  bool binary = false;
  if(!strncmp(buffer, "BINARY",6)) binary = true;
    
  if(fscanf(fp, "%s %s", buffer, buffer2) != 2) { 
    cerr << "error in reading vtk files" << filename << endl;
    ExecError("load vtk files");
  }
  if(strcmp(buffer, "DATASET") || strcmp(buffer2, "UNSTRUCTURED_GRID")){
    cout << "VTK reader can only read unstructured datasets" << endl;
    ExecError("load  vtk files");
    exit(1);
  }
    
  // read mesh vertices
  if(fscanf(fp, "%s %d %s\n", buffer, &nv, buffer2) != 3)
    { cout << "error in reading vtk files" << endl;
      ExecError("load  vtk files reah vertices");
      exit(1); }
  if(strcmp(buffer, "POINTS") || !nv){
    cerr << "No points in dataset" << endl;
      ExecError("load  vtk files: No points in dataset");
    exit(1);
  }
  int datasize;
  if(!strncmp(buffer2, "double",6))
    datasize = sizeof(double);
  else if(!strncmp(buffer2, "float",5))
    datasize = sizeof(float);
  else {
    cout <<"VTK reader only accepts float or double datasets" << endl;
    ExecError("load  vtk files VTK reader only accepts float or double datasets");
    exit(1);
  }
  if(verbosity>1)  
    cout << "   vtkio: Reading %d points" << nv << endl;
  vff = new Mesh::Vertex[nv];
    
  for(int i = 0 ; i < nv; i++){
    double xyz[3];
    if(binary){
      if(datasize == sizeof(float)){
	float f[3];
	if(fread(f, sizeof(float), 3, fp) != 3){ 
	  ExecError("load  vtk files VTK reader only accepts float or double datasets");
	  ExecError("error in reading vtk file"); 
	} 
	if(!bigEndian) SwapBytes((char*)f, sizeof(float), 3);
	for(int j = 0; j < 3; j++) xyz[j] = f[j];
      }
      else{
	if(fread(xyz, sizeof(double), 3, fp) != 3) { cout << "error in reading vtk files" << endl; ExecError("error in reading vtk file"); }
	if(!bigEndian) SwapBytes((char*)xyz, sizeof(double), 3);
      }
    }
    else{
      if(fscanf(fp, "%lf %lf %lf", &xyz[0], &xyz[1], &xyz[2]) != 3){ cout << "error in reading vtk files" << endl; ExecError("error in reading vtk file"); }
    }    
    vff[i].x = xyz[0];
    vff[i].y = xyz[1];
    if( abs(xyz[2]) > 1.e-7){
      cout << "we are plotted a two dimensional mesh: z coordinate must be 0" << endl;
      ExecError("error in reading vtk file,we are plotted a two dimensional mesh: z coordinate must be 0");
    }    
    vff[i].lab = 1;
  }    

  // read mesh elements
  int numElements, numElements2, totalNumInt;
  if(fscanf(fp, "%s %d %d\n", buffer, &numElements, &totalNumInt) != 3){ cout << "error in reading vtk files" << endl; ExecError("error in reading vtk file"); }
  if(strcmp(buffer, "CELLS") || !numElements){
    cout << "No cells in dataset" << endl;
    ExecError("error in reading vtk file");
  }
  cout << "Reading cells" << numElements << endl;
  
  int *IntCells   = new int[totalNumInt-numElements];
  int *firstCell  = new int[numElements+1];
  int *TypeCells   = new int[numElements];
  int numIntCells = 0;
  
  for(unsigned int i = 0; i < numElements; i++){
    int numVerts, n[100];
    for(int ii = 0; ii < 100; ii++) n[ii]=-1;
    if(binary){
      if( fread(&numVerts, sizeof(int), 1, fp) != 1 ) {
	cout << "error in reading VTK files " << endl;
	ExecError("error in reading vtk file");
	}
      if( !bigEndian) SwapBytes((char*)&numVerts, sizeof(int), 1);
      if((int)fread(n, sizeof(int), numVerts, fp) != numVerts){
	cout << "error in reading VTK files " << endl;
	ExecError("error in reading vtk file");
      }
      if(!bigEndian) SwapBytes((char*)n, sizeof(int), numVerts);
    }
    else{
      if(fscanf(fp, "%d", &numVerts) != 1){
	cout << "error in reading VTK files " << endl;
	  ExecError("error in reading vtk file");
      }
      for(int j = 0; j < numVerts; j++){
	if(fscanf(fp, "%d", &n[j]) != 1){
	  cout << "error in reading VTK files " << endl;
	  ExecError("error in reading vtk file");
	}
      }
    }
    firstCell[i] = numIntCells;
    for(int j = 0; j < numVerts; j++){
      if(n[j] >= 0 && n[j] < nv){
	IntCells[numIntCells] = n[j];
	numIntCells++;
      }	  
      else{
	cout << "Bad vertex index" << endl;
	ExecError("error in reading vtk file");
	}
    }
  }
  firstCell[numElements] = totalNumInt-numElements;  
  
  if(fscanf(fp, "%s %d\n", buffer, &numElements2) != 2){
    cout << " Error in reading CELL_TYPES ARGUMENT " << endl;
    ExecError("error in reading vtk file");
  }
  if(strcmp(buffer, "CELL_TYPES") || numElements2 != (int)numElements){
    cout <<"No or invalid number of cells types" << endl;
    ExecError("error in reading vtk file");
  }
  
  
  // 2D
  
  for(unsigned int i = 0; i < numElements; i++){
    int type;
    if(binary){
	if(fread(&type, sizeof(int), 1, fp) != 1){
	  cout <<"bug in readings cell types" << endl;
	  ExecError("error in reading vtk file");
	}
	if(!bigEndian) SwapBytes((char*)&type, sizeof(int), 1);
    }
    else{
      if(fscanf(fp, "%d", &type) != 1){
	cout <<"bug in readings cell types" << endl;
	ExecError("error in reading vtk file");
      }
    }
    TypeCells[i] = type;
    switch(type){
    case 1:  // Vertex
      cout << "this type of cell is not taking account in Freefem++ " << endl;
      break;
    case 3:  // Edge/line
      nbe++; // 2D	
      break;
    case 5:  // Triangle
      nt++; // 2D
      break;   
    case 10: // Tetrahèdre
      cout << "We are loading a three dimensional mesh. Three is no tetrahedron." << endl;
      ExecError("error in reading vtk file");
      break;  
    default: 
      cout << "Error :: This type of cell is not considered in Freefem++"<< endl;
      ExecError("error in reading vtk file");
      break;
      }
  }  
  /*
  // 3D 

    for(unsigned int i = 0; i < numElements; i++){
      int type;
      if(binary){
	if(fread(&type, sizeof(int), 1, fp) != 1){
	  cout <<"bug in readings cell types" << endl;
	  ExecError("error in reading vtk file");
	}
	if(!bigEndian) SwapBytes((char*)&type, sizeof(int), 1);
      }
      else{
	if(fscanf(fp, "%d", &type) != 1){
	  cout <<"bug in readings cell types" << endl;
	  ExecError("error in reading vtk file");
	}
      }
      TypeCells[i] = type;
      switch(type){
      case 1:  // Vertex
	cout << "this type of cell is not taking account in Freefem++ " << endl;
	break;
      case 3:  // Edge/line
	cout << "this type of cell is not taking account in Freefem++ for a two dimensional mesh" << endl; // 3D
	break;
      case 5:  // Triangle
	nbe++; // 3D
	break;   
      case 10: // Tetrahèdre
	nt++;
	break;  
      default: 
	cout << "Error :: This type of cell is not considered in Freefem++"<< endl;
	ExecError("error in reading vtk file");
	break;
	
      }
    }  
    */
    fclose(fp);

    // 2D Versions
    Mesh::Triangle *tff; 
    if(nt>0) tff = new Mesh::Triangle[nt];
    Mesh::Triangle *ttff = tff;
    
    Mesh::BorderElement *bff;
    if(nbe>0) bff= new Mesh::BorderElement[nbe];
    Mesh::BorderElement *bbff = bff;

    for(unsigned int i = 0; i < numElements; i++){
      int type=TypeCells[i];
      int iv[3];
      int label=1;
      switch(type){
      case 1:  // Vertex
	cout << "this type of cell is not taking account in Freefem++ " << endl;
	break;
      case 3:  // Edge/line
	assert( (firstCell[i+1]-firstCell[i]) == 2 );
	for(int j=firstCell[i]; j<firstCell[i+1]; j++){
    	  iv[j-firstCell[i]] = IntCells[j];
	}
	(bbff++)->set(vff, iv[0], iv[1], label);
	break;
      case 5:  // Triangle
	assert( (firstCell[i+1]-firstCell[i]) == 3 );
	for(int j=firstCell[i]; j<firstCell[i+1]; j++){
    	  iv[j-firstCell[i]] = IntCells[j];
	}
	(ttff++)->set(vff, iv[0], iv[1], iv[2], label);
	break;   
      default: 
	break;	
      } 
    }
    /*
    // 3D versions
    Tet *tff  = new Tet[nt];
    Tet *ttff = tff;

    Triangle3 *bff  = new Triangle3[nbe];
    Triangle3 *bbff = bff;

    for(unsigned int i = 0; i < numElements; i++){
      int type=TypeCells[i];
      int ivb[3],ivt[4];
      int label=1;
      switch(type){
      case 1:  // Vertex
	cout << "this type of cell is not taking account in Freefem++ " << endl;
	break;
      case 3:  // Edge/line
	break;
      case 5:  // Triangle
	assert( (firstCell[i+1]-firstCell[i]) == 3 );
	for(int j=firstCell[i]; j<firstCell[i+1]; j++){
    	  ivb[j-firstCell[i]] = IntCells[j];
	}
	(bbff++)->set(vff, ivb, label);
	break;   
      case 10: // Tetrahèdre
	assert( (firstCell[i+1]-firstCell[i]) == 4 );
	for(int j=firstCell[i]; j<firstCell[i+1]; j++){
    	  ivt[j-firstCell[i]] = IntCells[j];
	}
	(ttff++)->set(vff, ivt, label);
	break;  
      default: 
	break;	
      } 
    }
    */
    delete [] IntCells;   
    delete [] firstCell; 
    delete [] TypeCells;
    
    Mesh *pTh = new Mesh(nv,nt,nbe,vff,tff,bff);
    return pTh;
}

AnyType VTK_LoadMesh_Op::operator()(Stack stack)  const 
{
 
  string * pffname= GetAny<string *>((*filename)(stack));
  bool swap = false;
  int reftri = 1;
  int refedges = 1;
  if( nargs[0] )  reftri = GetAny< int >((*nargs[0])(stack));
  if( nargs[1] )  swap = GetAny< bool >((*nargs[1])(stack));
  if( nargs[2] )  refedges = GetAny< int >((*nargs[2])(stack));

  string *DataLabel;
  if( nargs[3] ) DataLabel = GetAny<string *>((*nargs[3])(stack));

  Mesh * Th = VTK_Load( *pffname , swap); 
 
  // A faire fonction pour changer le label

  Add2StackOfPtr2FreeRC(stack,Th);
		
  return Th;
}

//==============================================
// ECRITURE DE FICHIER .vtk
//==============================================

class VTK_WriteMesh_Op : public E_F0mps 
{
public:
  typedef long  Result;
  Expression eTh;
  Expression filename;
  
  struct Expression2 {
    string name;
    long what; // 1 scalar, 2 vector, 3 symtensor
    long nbfloat; // 1 scalar, 2 vector (3D), 3 symtensor(3D)
    Expression e[3];
    Expression2() {e[0]=0; e[1]=0; e[2]=0;  what=0; nbfloat=0; };
    Expression &operator[](int i){return e[i];}
    double eval(int i,Stack stack) const  { 
    if (e[i]) {
      return GetAny< double >( (*e[i])(stack) );
    }
    else 
      return 0;
    }
  };
  vector<Expression2> l;
  static const int n_name_param = 7;  
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  long    arg(int i,Stack stack,long a ) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
  
public:
  VTK_WriteMesh_Op(const basicAC_F0 &  args) : l( args.size()-2 )
  {
    int nbofsol;
    int ddim=2;
    int stsize=3;

    int sca=0,vec=0,ten=0;
    
    string scas("scalaire"); 
    string vecs("vector"); 
    string tens("tensor"); 

    if(verbosity) cout << "Write Mesh and Solutions in VTK Formats" << endl;
    args.SetNameParam(n_name_param,name_param,nargs);

    if (BCastTo<string *>(args[0])) filename = CastTo<string *>(args[0]);
    if (BCastTo<pmesh>(args[1])) eTh= CastTo<pmesh>(args[1]);
  
    nbofsol = l.size();
    for (size_t i=2;i<args.size();i++){
      size_t jj=i-2;
      
      if ( BCastTo<double>( args[i] ))
	{
	  l[jj].what=1;
	  l[jj].nbfloat=1;
	  l[jj][0]=to<double>( args[i] );
	  
	  char number[16];
	  sprintf(number,"%li",jj+1);
	  l[jj].name=scas;
	  l[jj].name+=number;
	  sca++;	  
	}
      else if ( args[i].left()==atype<E_Array>() )
	{
	  const E_Array * a0  = dynamic_cast<const E_Array *>( args[i].LeftValue() );
	  //cout << "taille" << a0->size() << endl;
	  if (a0->size() != ddim && a0->size() != stsize) 
	    CompileError("savesol in 2D: vector solution is 2 composant, tensor solution is 3 composant");
	  
	  if( a0->size() == ddim){
	    // vector solution
	    l[jj].what=2;
	    l[jj].nbfloat=ddim;
	   
	    for(int j=0; j<ddim; j++){
	      l[jj][j] = to<double>( (*a0)[j]);
	    }

	    char number[16];
	    sprintf(number,"%li",jj+1);
	    l[jj].name=vecs;
	    l[jj].name+=number;
	    vec++;	  
	  }
	  else if( a0->size() == stsize){
	    // symmetric tensor solution
	    l[jj].what=3;
	    l[jj].nbfloat=stsize;
	    for(int j=0; j<stsize; j++){
	      l[jj][j] = to<double>( (*a0)[j]);
	    }
	    char number[16];
	    sprintf(number,"%li",jj+1);
	    l[jj].name=tens;
	    l[jj].name+=number;
	    ten++;	
	  }
	  
	}
      else {
	cout << " arg " << i << " " << args[i].left() << endl;
	CompileError("save solution in 2D in format VTK: Sorry no way to save this kind of data");
      }
    }
    
  }
  static ArrayOfaType  typeargs() { return  ArrayOfaType( atype<string *>(), atype<pmesh>(), true); }// all type
  static  E_F0 * f(const basicAC_F0 & args) { return new VTK_WriteMesh_Op(args);}   
  AnyType operator()(Stack stack)  const ;
};
basicAC_F0::name_and_type VTK_WriteMesh_Op::name_param[]= {
  {  "dataname", &typeid(string*)},
  {  "withsurfacemesh", &typeid(bool)},
  {  "order", &typeid(KN_<long>)},
  // A rajouter dans le 3D
  {  "floatmesh", &typeid(bool)},
  {  "floatsol", &typeid(bool)},
  {  "bin", &typeid(bool)},
  {  "swap", &typeid(bool)}
};

void VTK_WRITE_MESH( const string &filename, FILE *fp, const Mesh &Th, bool binary, int datasize, bool surface, bool bigEndian){
     
  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "%s, Created by Freefem++ \n", filename.c_str());
  if(binary)
    fprintf(fp, "BINARY\n");
  else
    fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
  // get all the entities in the model

  // write mesh vertices
  
  if(datasize == sizeof(float)) 
    {      
      fprintf(fp, "POINTS %d float\n", Th.nv);
      for(unsigned int i = 0; i < Th.nv; i++){
	const  Mesh::Vertex & P = Th.vertices[i];
	float f[3];
	f[0]=P.x;
	f[1]=P.y;
	f[2]=0;  // P.z; 3D case
	if(binary){
	  if(!bigEndian) SwapBytes((char*)&f, sizeof(float), 3);
	  fwrite(&f, sizeof(float), 3, fp);
	}
      else{
	fprintf(fp,"%.8g %.8g %.8g\n",P.x,P.y,0.);
      }
    }  
  }
  else if(datasize == sizeof(double)){
    fprintf(fp, "POINTS %d double\n", Th.nv);
    for(unsigned int i = 0; i < Th.nv; i++){
      const  Mesh::Vertex & P = Th.vertices[i];
      double f[3];
      f[0]=P.x;
      f[1]=P.y;
      f[2]=0;  // P.z; 3D case
      if(binary){
	if(!bigEndian) SwapBytes((char*)&f, sizeof(float), 3);
	fwrite(&f, sizeof(float), 3, fp);
      }
      else{
	fprintf(fp,"%.15lg %.15lg %.15lg\n",f[0],f[1],f[2]);
      }
    }  
  }
  fprintf(fp,"\n");
  if(verbosity > 1) printf("writting vertices is finish\n");
  if(verbosity > 1) printf("writting elements now\n");
  
  //================
  //    CELL 
  //================
  // loop over all elements we need to save and count vertices
  int numElements, totalNumInt;
  if(surface){
    numElements = Th.nt+Th.neb; 
    totalNumInt = Th.nt*3+Th.neb*2+numElements;
  }
  else{
    numElements = Th.nt; 
    totalNumInt = Th.nt*3+numElements;
  }
  if(verbosity > 1) printf("writting cells \n");
  // print vertex indices in ascii or binary
  fprintf(fp, "CELLS %d %d\n", numElements, totalNumInt);

  if(binary){
    int IntType=3;
    if(verbosity > 1) printf("writting triangle elements \n");
    for(int it=0; it< Th.nt; it++){
     
      const Mesh::Triangle & K( Th.t(it) );
      int iv[IntType+1];
      
      iv[0] = IntType;
      for(int ii=0; ii<IntType; ii++){
	iv[ii+1] = Th.operator()(K[ii]);
      }
      
      if(!bigEndian) SwapBytes((char*)&iv, sizeof(int), IntType+1);
      fwrite(&iv, sizeof(int), IntType+1, fp);
      
    }
    if(surface){
      if(verbosity > 1) printf("writting edge elements \n");
      IntType=2;
      for(int ibe=0; ibe<Th.neb; ibe++){
	 const Mesh::BorderElement &K( Th.be(ibe) );
      
	 int iv[IntType+1];    
	 iv[0] = IntType;
	 for(int ii=0; ii<IntType; ii++){
	   iv[ii+1] = Th.operator()(K[ii]);
	 }
	 
	 if(!bigEndian) SwapBytes((char*)&iv, sizeof(int), IntType+1);
	 fwrite(&iv, sizeof(int), IntType+1, fp);
      }
    }
  }
  else{
    int IntType=3;
    if(verbosity > 1) printf("writting triangle elements \n");
    for(int it=0; it< Th.nt; it++){
      const Mesh::Triangle &K( Th.t(it) );

      int iv[IntType+1];      
      iv[0] = IntType;
      for(int ii=0; ii<IntType; ii++){
	iv[ii+1] = Th.operator()(K[ii]);
      }
      fprintf(fp,"%d %d %d %d\n", iv[0],iv[1],iv[2],iv[3]);
    }
    if(surface){
      if(verbosity > 1) printf("writting edge elements \n");
      IntType=2;
      for(int ibe=0; ibe<Th.neb; ibe++){
	const Mesh::BorderElement &K( Th.be(ibe) );
	
	int iv[IntType+1];    
	iv[0] = IntType;
	for(int ii=0; ii<IntType; ii++){
	  iv[ii+1] = Th.operator()(K[ii]);
	}
	
	fprintf(fp,"%d %d %d\n", iv[0],iv[1],iv[2]);
      }
    }
  }
  fprintf(fp, "\n");


  // CELL_TYPE
  // print element types in ascii or binary  
  
  fprintf(fp, "CELL_TYPES %d\n", numElements); 
  if(binary){
    int type;    
    type = VTK_TRI;  
    for(int it=0; it< Th.nt; it++){
     
      if(!bigEndian) SwapBytes((char*)&type, sizeof(int), 1);
      fwrite(&type, sizeof(int), 1, fp);
    }
    if(surface){
      type=VTK_EDGE;
      for(int ibe=0; ibe<Th.neb; ibe++){
	 
	if(!bigEndian) SwapBytes((char*)&type, sizeof(int), 1);
	fwrite(&type, sizeof(int), 1, fp);
      }
    }
  }
  else{
    int type;
    type= VTK_TRI;
    for(int it=0; it< Th.nt; it++){      
      fprintf(fp,"%d ",type);
    }
    if(surface){
      type=VTK_EDGE;
      for(int ibe=0; ibe<Th.neb; ibe++){
	fprintf(fp,"%d%c",type,(ibe%10==9)? '\n' : ' ');
      }
    }
  } 

  fprintf(fp, "\n");

  //=================================
  //  WRITE SOLUTION IN FORMAT VTK
  //     LABEL OF ELEMENTS
  //=================================

  list<int> list_label_Elem;
  //list<int> list_label_Border_Elem;
  {
    
    for(int it=0; it< Th.nt; it++){
      const Mesh::Triangle &K( Th.t(it) );
      list<int>::const_iterator ilist;
      int labOk=0;
      for( ilist=list_label_Elem.begin(); ilist!=list_label_Elem.end(); ilist++){
	  if( *ilist == K.lab ){ labOk = 1;   break; } 
	}
	if( labOk == 0){
	  list_label_Elem.push_back(K.lab);
	}
    }
    
    if(surface){
      for(int ibe=0; ibe<Th.neb; ibe++){
	const Mesh::BorderElement &K( Th.be(ibe) );
	list<int>::const_iterator ilist;
	int labOk=0;
	for( ilist=list_label_Elem.begin(); ilist!=list_label_Elem.end(); ilist++){
	  if( *ilist == K.lab ){ labOk = 1;   break; } 
	}
	if( labOk == 0){
	  list_label_Elem.push_back(K.lab);
	}
      
      }
    }
  
  }
  list_label_Elem.sort();

  //=================================
  //=================================

  fprintf(fp, "CELL_DATA %d\n", numElements);
  int cell_fd=1;
  int cell_lab=1;
  fprintf(fp, "Scalars  Label int %d\n", cell_fd);
  fprintf(fp, "LOOKUP_TABLE FreeFempp_table\n"); 
  // Determination des labels
  if(binary){
    int label;    
    for(int it=0; it< Th.nt; it++){
      const Mesh::Triangle &K( Th.t(it) );
      label =K.lab;
      if(!bigEndian) SwapBytes((char*)&label, sizeof(int), 1);
      fwrite(&label, sizeof(int), 1, fp);
    }
    if(surface){
      for(int ibe=0; ibe<Th.neb; ibe++){
	const Mesh::BorderElement &K( Th.be(ibe) );
	label =K.lab;
	if(!bigEndian) SwapBytes((char*)&label, sizeof(int), 1);
	fwrite(&label, sizeof(int), 1, fp);
      }
    }
  }
  else{
    int label;
    for(int it=0; it< Th.nt; it++){ 
      const Mesh::Triangle &K( Th.t(it) );
      label =K.lab;
      fprintf(fp,"%d\n",label);
    }
    if(surface){
      for(int ibe=0; ibe<Th.neb; ibe++){
	const Mesh::BorderElement &K( Th.be(ibe) );
	label =K.lab;
	fprintf(fp,"%d\n",label);
      }
    }
  } 
  fprintf(fp,"\n");
  int size_list=0;
  list<int>::const_iterator ilist;
  for( ilist=list_label_Elem.begin(); ilist!=list_label_Elem.end(); ilist++) size_list++;

  fprintf(fp, "LOOKUP_TABLE FreeFempp_table %d\n",size_list);
  { list<int>::const_iterator ilist;
    for( ilist=list_label_Elem.begin(); ilist!=list_label_Elem.end(); ilist++){
      float tab[4];
      tab[0] = ColorTable[abs(*ilist)%NbColorTable][0];
      tab[1] = ColorTable[abs(*ilist)%NbColorTable][1];
      tab[2] = ColorTable[abs(*ilist)%NbColorTable][2];
      tab[3] = 1.0;
      if(binary){
	if(!bigEndian) SwapBytes((char*)&tab, sizeof(float), 4);
	fwrite(&tab, sizeof(float), 4, fp);
      }
      else
	fprintf(fp,"%.8f %.8f %.8f %.8f\n",tab[0],tab[1],tab[2],tab[3]);
    }
  }
}

AnyType VTK_WriteMesh_Op::operator()(Stack stack)  const 
{
 
  string * pffname= GetAny<string *>((*filename)(stack));
  Mesh * pTh= GetAny<Mesh *>((*eTh)(stack));
  ffassert(pTh);
  Mesh &Th=*pTh;
  bool swap = false;
  bool bigEndian=false;
  bool binary = false;
  bool surface = true;
  bool floatmesh = true;
  bool floatsol = true;
  int datasize = sizeof(float);
  int datasizeSol = sizeof(float);
  string *dataname;
  int nbofsol = l.size();
  KN<int> order(nbofsol);
  
  char *nameofuser[nbofsol];

  for(int ii=0; ii<nbofsol; ii++) order[ii] = 0;

  if( nargs[0] ) dataname  = GetAny< string* >( (*nargs[0])(stack) );  
  if( nargs[1] ) surface   = GetAny< bool >( (*nargs[1])(stack) );
  if( nargs[2] ) order     = GetAny< KN_<long> >( (*nargs[2])(stack) );
  if( nargs[3] ) floatmesh = GetAny< bool >( (*nargs[3])(stack) );
  if( nargs[4] ) floatsol  = GetAny< bool  >( (*nargs[4])(stack) );
  if( nargs[5] ) binary    = GetAny< bool  >( (*nargs[5])(stack) );
  if( nargs[6] ) bigEndian = GetAny< bool  >( (*nargs[6])(stack) );

  swap = bigEndian;

  if( !floatmesh ) datasize = sizeof(double);
  if( !floatsol ) datasizeSol= sizeof(double);

  int iii=0;
  if( nargs[0])
    {
      char *data = newcopy(dataname->c_str());
      char * name = strtok(data," \t\n");
      
      nameofuser[iii] = newcopy(name);
      if(verbosity>5)
	cout << "   iovtk : value of iii  =" << iii << "  \""<<  nameofuser[iii] <<  "\"\n";
      iii++;
      {
	
	while( (name = strtok(NULL," \t\n\0")) ){
	  
	  if( iii >= nbofsol ){
	    if(verbosity>5)
	      cout << "   iovtk : The number of data name istoo large " << endl;
	    break;
	  }
	  nameofuser[iii] = newcopy(name);
	  if(verbosity>5)
	    cout << "   iovtk : value of iii  =" << iii << "  \""<<  nameofuser[iii] <<  "\"\n";
	  iii++;
	}
	if( iii < nbofsol){	
	  if(verbosity>6)
	    cout << "   iovtk:  The number of data name is too small, we give default name " << endl;
	}
	delete [] data;
      }
    }
  if( iii < nbofsol ){
    for( int iiii=iii; iiii<nbofsol; iiii++){
      //      char *dataff = new char[l[iii].name.size()+1];
      //strcpy(dataff, l[iii].name.c_str());
      nameofuser[iiii] = newcopy(l[iii].name.c_str());//dataff;
      
    }
  }

  
  FILE *fp = fopen( (*pffname).c_str(), "wb");
  if(!fp){
    cerr << "Unable to open file " << (*pffname).c_str() << endl;
    ExecError("error in reading vtk file");
  }
  
  VTK_WRITE_MESH( *pffname, fp, Th, binary, datasize, surface, swap);    		

  
  
  // determination of number of order 0 et 1.
  int Norder0=0;
  for(int ii=0; ii< nbofsol; ii++)
    if(order[ii] == 0) Norder0++;


  if( datasizeSol == sizeof(float) ){
    
    if( Norder0 >0){
      fprintf(fp, "FIELD FieldData %d\n", Norder0);
      for(int ii=0; ii< nbofsol; ii++){
	if(order[ii] == 0){
	  int nsol;
	  if(surface){
	    nsol = Th.nt+Th.neb;
	  }
	  else{
	    nsol = Th.nt;
	  }
	  	  
	  //fprintf(fp,"%s %d %d float\n",l[ii].name.c_str(),l[ii].nbfloat,nsol); 
	  fprintf(fp,"%s %ld %d float\n",nameofuser[ii],l[ii].nbfloat,nsol); 
	  if(verbosity>5)
	  cout << "name of data("<< ii <<")=" << nameofuser[ii] << " " << l[ii].name << endl;
	 
//	  if(l[ii].what == 1) fprintf(fp,"%s %d %d float\n",l[ii].name,ii,l[ii].nbfloat,nsol); 
//        if(l[ii].what == 2) fprintf(fp,"%s %d %d float\n",l[ii].name,ii,l[ii].nbfloat,nsol); 
//        if(l[ii].what == 3) fprintf(fp,"%s %d %d float\n",l[ii].name,ii,l[ii].nbfloat,nsol); 
// 	  }
// 	  else{
// 	    if(l[ii].what == 1) fprintf(fp,"scalaire%d %d %d float\n",ii,l[ii].nbfloat,nsol); 
// 	    if(l[ii].what == 2) fprintf(fp,"vector%d %d %d float\n",ii,l[ii].nbfloat,nsol); 
// 	    if(l[ii].what == 3) fprintf(fp,"Tensor%d %d %d float\n",ii,l[ii].nbfloat,nsol); 
// 	  }

	  MeshPoint *mp3(MeshPointStack(stack)); 
	  R2 Cdg_hat = R2(1./3.,1./3.);  
	
	  for (int it=0;it<Th.nt;it++){
	    const Mesh::Triangle  & K(Th.t(it));
	    mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);

	    for(int j=0;j<l[ii].nbfloat;j++){
	      float value = l[ii].eval(j,stack);
	      if(binary){
		if(!bigEndian) SwapBytes((char*)&value, sizeof(float), 1);
		fwrite(&value, sizeof(float), 1, fp);
	      }
	      else{
		fprintf(fp,"%lf ",value); 
	      }
	    } 
	  }
	  if( surface ){
	    for (int ibe=0;ibe<Th.neb;ibe++){
	      // determination du triangle contenant cette edge
	      int ie;
	      int it = Th.BoundaryElement( ibe, ie); 
	      const Mesh::Triangle  & K(Th.t(it));
	      mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);

	      for(int j=0;j<l[ii].nbfloat;j++){
		float value = l[ii].eval(j,stack);
		if(binary){
		  if(!bigEndian) SwapBytes((char*)&value, sizeof(float), 1);
		  fwrite(&value, sizeof(float), 1, fp);
		}
		else{
		  fprintf(fp,"%lf ",value); 
		}
	      } 

	    }
	  }

	  fprintf(fp,"\n");
	}
      }
    }
    if( Norder0 < nbofsol ){
      fprintf(fp, "POINT_DATA %d\n", Th.nv);
      fprintf(fp, "FIELD FieldData %d\n", nbofsol-Norder0);
      for(int ii=0; ii< nbofsol; ii++){
	if(order[ii] == 1){
	
	  //fprintf(fp,"%s %d %d float\n",l[ii].name.c_str(),l[ii].nbfloat,Th.nv); 
	  fprintf(fp,"%s %ld %d float\n",nameofuser[ii],l[ii].nbfloat,Th.nv); 
	  if(verbosity>5)
	  cout << "name of data("<< ii <<")=" << nameofuser[ii] << " " << l[ii].name << endl;
	  
	  /*
	  if(l[ii].what == 1){	  
	    //fprintf(fp, "Scalars  pressure%d float %d\n",ii,l[ii].what);
	    //fprintf(fp, "LOOKUP_TABLE default\n");
	    fprintf(fp,"pressure%d %d %d float\n",ii,l[ii].nbfloat,Th.nv);
	  }
	  if(l[ii].what == 2){
	    //fprintf(fp, "FIELD FieldData 1\n" );
	    fprintf(fp,"Vitesse%d %d %d float\n",ii,l[ii].nbfloat,Th.nv); 
	  }
	  if(l[ii].what == 3){
	    //fprintf(fp, "FIELD FieldData 1\n" );
	    fprintf(fp,"Tensor%d %d %d float\n",ii,l[ii].nbfloat,Th.nv); 
	  }
	  */
	  MeshPoint *mp3(MeshPointStack(stack)); 
	  KN<double> valsol(Th.nv*l[ii].nbfloat);
	  KN<int> takemesh(Th.nv);
	  takemesh =0;
	  valsol   =0.;
	  for(int it=0;  it<Th.nt; it++){
	    for(int iv=0; iv<3; iv++){
	      int i=Th(it,iv);
	      mp3->setP(&Th,it,iv);
 
	      for(int j=0;j<l[ii].nbfloat;j++){
		valsol[ i*l[ii].nbfloat+j ] = valsol[ i*l[ii].nbfloat+j ] + l[ii].eval(j,stack);	     
	      }
	      takemesh[i] = takemesh[i]+1;    
	    }
	  }
	
	  for (int iv=0;iv<Th.nv;iv++)
	    {
	    for(int j=0;j<l[ii].nbfloat;j++)
	      {
		valsol[iv*l[ii].nbfloat+j] = valsol[ iv*l[ii].nbfloat+j ]/takemesh[iv];
		float value = valsol[iv*l[ii].nbfloat+j];
		if(binary)
		  {
		    if(!bigEndian) SwapBytes((char*)&value, sizeof(float), 1);
		    fwrite(&value, sizeof(float), 1, fp);
		  }
		else
		  {
		  fprintf(fp,"%.8lf ",value); 
		  }
	      }	 
	    if(!binary)
	      fprintf(fp,"\n");  
	    }
	  if(!binary)	  
	    fprintf(fp,"\n");
	}
      }
    }    
  }

  // datasizeSol == sizeof(double)


  if( datasizeSol == sizeof(double) ){
    
    if( Norder0 >0){
      fprintf(fp, "FIELD FieldData %d\n", Norder0);
      for(int ii=0; ii< nbofsol; ii++){
	if(order[ii] == 0){
	  int nsol;
	  if(surface){
	    nsol = Th.nt+Th.neb;
	  }
	  else{
	    nsol = Th.nt;
	  }
	  
	  fprintf(fp,"%s %ld %d float\n",nameofuser[ii],l[ii].nbfloat,nsol); 
	  if(verbosity>5)
	  cout << "name of data("<< ii <<")=" << nameofuser[ii]  << endl;
	  /*
	  if(l[ii].what == 1) fprintf(fp,"pressure%d %d %d double\n",ii,l[ii].nbfloat,nsol); 
	  if(l[ii].what == 2) fprintf(fp,"Vitesse%d %d %d double\n",ii,l[ii].nbfloat,nsol); 
	  if(l[ii].what == 3) fprintf(fp,"Tensor%d %d %d double\n",ii,l[ii].nbfloat,nsol); 
	  */
	  MeshPoint *mp3(MeshPointStack(stack)); 
	R2 Cdg_hat = R2(1./3.,1./3.);  
	
	for (int it=0;it<Th.nt;it++){
	  const Mesh::Triangle  & K(Th.t(it));
	  mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);

	  for(int j=0;j<l[ii].nbfloat;j++)
	    {
	      double value = l[ii].eval(j,stack);
	      if(binary){
		if(!bigEndian) SwapBytes((char*)&value, sizeof(double), 1);
		fwrite(&value, sizeof(double), 1, fp);
	      }
	      else{
		fprintf(fp,"%.8f\n",value); 
	      }
	    } 
	}
	if( surface ){
	  for (int ibe=0;ibe<Th.neb;ibe++){
	    // determination du triangle contenant cette edge
	    int ie;
	    int it = Th.BoundaryElement( ibe, ie); 
	    const Mesh::Triangle  & K(Th.t(it));
	    mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);

	    for(int j=0;j<l[ii].nbfloat;j++){
	      double value = l[ii].eval(j,stack);
	      if(binary){
		if(!bigEndian) SwapBytes((char*)&value, sizeof(double), 1);
		fwrite(&value, sizeof(double), 1, fp);
	      }
	      else{
		fprintf(fp,"%.8f\n",value); 
	      }
	    } 

	  }
	}

	fprintf(fp,"\n");
	}
      }
    }
    if( Norder0 < nbofsol ){
      fprintf(fp, "POINT_DATA %d\n", Th.nv);
      fprintf(fp, "FIELD FieldData %d\n", nbofsol-Norder0);
      for(int ii=0; ii< nbofsol; ii++){
	if(order[ii] == 1){

	  
	  fprintf(fp,"%s %ld %d float\n",nameofuser[ii],l[ii].nbfloat,Th.nv); 
	  if(verbosity>5)
	  cout << "name of data("<< ii <<")=" << nameofuser[ii]  << endl;
	  /*
	    if(nargs[0]){
	    if(l[ii].what == 1) fprintf(fp,"%s %d %d double\n",l[ii].name,ii,l[ii].nbfloat,Th.nv); 
	    if(l[ii].what == 2) fprintf(fp,"%s %d %d double\n",l[ii].name,ii,l[ii].nbfloat,Th.nv); 
	    if(l[ii].what == 3) fprintf(fp,"%s %d %d double\n",l[ii].name,ii,l[ii].nbfloat,Th.nv); 
	    }
	    else{
	    if(l[ii].what == 1){	  
	    fprintf(fp,"scalaire%d %d %d double \n",ii,l[ii].nbfloat,Th.nv);
	    }
	    if(l[ii].what == 2){
	    fprintf(fp,"vector%d %d %d double \n",ii,l[ii].nbfloat,Th.nv); 
	    }
	    if(l[ii].what == 3){   
	    fprintf(fp,"tensor%d %d %d double \n",ii,l[ii].nbfloat,Th.nv); 
	    }
	    }
	  */
	  MeshPoint *mp3(MeshPointStack(stack)); 
	  KN<double> valsol(Th.nv*l[ii].nbfloat);
	  KN<int> takemesh(Th.nv);
	  takemesh =0;
	  valsol   =0.;
	  for(int it=0;  it<Th.nt; it++){
	    for(int iv=0; iv<3; iv++){
	      int i=Th(it,iv);
	      mp3->setP(&Th,it,iv);
	      
	      for(int j=0;j<l[ii].nbfloat;j++){
		valsol[ i*l[ii].nbfloat+j ] = valsol[ i*l[ii].nbfloat+j ] + l[ii].eval(j,stack);	     
	      }
	    takemesh[i] = takemesh[i]+1;    
	    }
	  }
	  
	  for (int iv=0;iv<Th.nv;iv++){
	    for(int j=0;j<l[ii].nbfloat;j++){
	      valsol[iv*l[ii].nbfloat+j] = valsol[ iv*l[ii].nbfloat+j ]/takemesh[iv];
	      double value = valsol[iv*l[ii].nbfloat+j];
	      if(binary){
		if(!bigEndian) SwapBytes((char*)&value, sizeof(double), 1);
		fwrite(&value, sizeof(double), 1, fp);
	      }
	      else{
		fprintf(fp,"%.8f\n",value); 
	      }
	    }	   
	  }
	fprintf(fp,"\n");
	}
      }
    }
  }

  fclose(fp);


  for( int iiii=0; iiii<nbofsol; iiii++){
    delete [] nameofuser[iiii];
  }
}



//==============================================
// FIN ECRITURE DE FICHIER .vtk (2D)
//==============================================

//=======================
// FIN 2D Fichier .vtk
//=======================

// three dimensional case
// LOAD fichier.vtk
class VTK_LoadMesh3_Op : public E_F0mps 
{
public:
  Expression filename;
  static const int n_name_param = 4; // 
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
public:
  VTK_LoadMesh3_Op(const basicAC_F0 &  args,Expression ffname) 
    : filename(ffname)
  {
    if(verbosity) cout << "Load mesh given by VTK " << endl;
    args.SetNameParam(n_name_param,name_param,nargs);   
  } 
  
  AnyType operator()(Stack stack)  const ;
};
basicAC_F0::name_and_type VTK_LoadMesh3_Op::name_param[]= {
  {  "reftet", &typeid(long)},
  {  "swap", &typeid(bool)},
  {  "refface", &typeid(long)},
  {  "namelabel",&typeid(string)}
};


class  VTK_LoadMesh3 : public OneOperator { public:  
    VTK_LoadMesh3() : OneOperator(atype<pmesh3>(),atype<string *>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
    return  new VTK_LoadMesh3_Op( args,t[0]->CastTo(args[0]) ); 
  }
};

Mesh3 * VTK_Load3(const string & filename, bool bigEndian)
// swap = bigEndian or not bigEndian
{
  // variable freefem++
  int nv, nt=0, nbe=0;
 
  // Reading Mesh in vtk formats 
  FILE *fp = fopen(filename.c_str(), "rb");
  if(!fp){
    cerr << "Unable to open file " << filename.c_str() << endl;
    ExecError("error in reading vtk file");
  }

  char buffer[256], buffer2[256];
    
  fgets(buffer, sizeof(buffer), fp); // version line
  fgets(buffer, sizeof(buffer), fp); // title
  
  fscanf(fp, "%s", buffer); // ASCII or BINARY
  bool binary = false;
  if( !strcmp(buffer, "BINARY") ) binary = true;
    
  if(fscanf(fp, "%s %s", buffer, buffer2) != 2){ cout << "error in reading vtk files" << endl; ExecError("error in reading vtk file"); }
  if(strcmp(buffer, "DATASET") || strcmp(buffer2, "UNSTRUCTURED_GRID")){
    cout << "VTK reader can only read unstructured datasets" << endl;
    ExecError("error in reading vtk file");
  }
    
  // read mesh vertices
  if(fscanf(fp, "%s %d %s\n", buffer, &nv, buffer2) != 3){ cout << "error in reading vtk files" << endl; ExecError("error in reading vtk file"); } 
  if(strcmp(buffer, "POINTS") || !nv){
    cerr << "No points in dataset" << endl;
    ExecError("error in reading vtk file");
  }
  int datasize;
  if( !strncmp(buffer2, "double",6))
    datasize = sizeof(double);
  else if( !strncmp(buffer2, "float",5))
    datasize = sizeof(float);
  else{
    cout << "VTK reader only accepts float or double datasets" << endl;
    ExecError("error in reading vtk file");
  }
    
  cout << "Reading %d points" << nv << " buffer2" <<  buffer2 << "binary" << binary << " " << datasize << " "<< sizeof(float) << endl;
  Vertex3   *vff = new Vertex3[nv];
    
  for(int i = 0 ; i < nv; i++){
    cout << " i=" << i << endl;
    double xyz[3];
    if(binary){
      if(datasize == sizeof(float)){
	float f[3];
	if(fread(f, sizeof(float), 3, fp) != 3){ cout << "error in reading vtk files" << endl; ExecError("error in reading vtk file"); }
	if(!bigEndian) SwapBytes((char*)f, sizeof(float), 3);
	for(int j = 0; j < 3; j++) xyz[j] = f[j];
      }
      else{
	if(fread(xyz, sizeof(double), 3, fp) != 3){ cout << "error in reading vtk files" << endl; ExecError("error in reading vtk file"); }
	if(!bigEndian) SwapBytes((char*)xyz, sizeof(double), 3);
      }
    }
    else{
      cout << datasize << " "<< sizeof(float) << endl;
      if(datasize == sizeof(float)){ 
	if(fscanf(fp, "%lf %lf %lf", &xyz[0], &xyz[1], &xyz[2]) != 3){ cout << "error in reading vtk files (float)" << endl; ExecError("error in reading vtk file"); }
      }
      else{
	if(fscanf(fp, "%lf %lf %lf", &xyz[0], &xyz[1], &xyz[2]) != 3){ cout << "error in reading vtk files (double)" << endl; ExecError("error in reading vtk file"); } 
      }
    }
    vff[i].x = xyz[0];
    vff[i].y = xyz[1];
    vff[i].z = xyz[2];
    vff[i].lab = 1;
    
    printf( "xyz = %f %f %f\n", xyz[0],xyz[1], xyz[2]);
  }    

  // read mesh elements
  int numElements, numElements2, totalNumInt;
  if(fscanf(fp, "%s %d %d\n", buffer, &numElements, &totalNumInt) != 3){ cout << "error in reading vtk files" << endl; ExecError("error in reading vtk file"); }
  printf("reading parameter %s %d %d\n", buffer, numElements, totalNumInt);
  if(strncmp(buffer, "CELLS",5) || !numElements){
    cout << "No cells in dataset" << endl;
    ExecError("error in reading vtk file");
  }
  cout << "Reading cells" << numElements << endl;
  
  int *IntCells   = new int[totalNumInt-numElements];
  int *firstCell  = new int[numElements+1];
  int *TypeCells   = new int[numElements];
  int numIntCells = 0;
  
  for(unsigned int i = 0; i < numElements; i++){
    int numVerts, n[100];
    cout << "i=" << i << " " << numElements << endl;
    for(int ii = 0; ii < 100; ii++) n[ii]=-1;
    if(binary){
      if( fread(&numVerts, sizeof(int), 1, fp) != 1 ) {
	cout << "error in reading VTK files " << endl;
	ExecError("error in reading vtk file");
      }
      if( !bigEndian) SwapBytes((char*)&numVerts, sizeof(int), 1);
      if((int)fread(n, sizeof(int), numVerts, fp) != numVerts){
	cout << "error in reading VTK files " << endl;
	ExecError("error in reading vtk file");
      }
      if(!bigEndian) SwapBytes((char*)n, sizeof(int), numVerts);
    }
    else{
      if(fscanf(fp, "%d", &numVerts) != 1){
	cout << "error in reading VTK files " << endl;
	ExecError("error in reading vtk file");
      }
      cout << "numVerts" << numVerts << endl;
      for(int j = 0; j < numVerts; j++){
	if(fscanf(fp, "%d", &n[j]) != 1){
	  cout << "error in reading VTK files " << endl;
	  ExecError("error in reading vtk file");
	}
	cout << "n[j]" << n[j] << endl;
      }
    }
    firstCell[i] = numIntCells;
    for(int j = 0; j < numVerts; j++){
      if(n[j] >= 0 && n[j] < nv){
	IntCells[numIntCells] = n[j];
	numIntCells++;
      }	  
      else{
	cout << "Bad vertex index" << endl;
	ExecError("error in reading vtk file");
      }
    }
  }
  firstCell[numElements] = totalNumInt-numElements;  
  
  if(fscanf(fp, "%s %d\n", buffer, &numElements2) != 2){
    cout << " Error in reading CELL_TYPES ARGUMENT " << endl;
    ExecError("error in reading vtk file");
  }
  if(strcmp(buffer, "CELL_TYPES") || numElements2 != (int)numElements){
    cout <<"No or invalid number of cells types" << endl;
    ExecError("error in reading vtk file");
  }
  
  printf( "reading parameter %s %d\n", buffer, numElements2);
  
  // 3D 
  
  for(unsigned int i = 0; i < numElements; i++){
    int type;
    if(binary){
      if(fread(&type, sizeof(int), 1, fp) != 1){
	cout <<"bug in readings cell types" << endl;
	ExecError("error in reading vtk file");
      }
      if(!bigEndian) SwapBytes((char*)&type, sizeof(int), 1);
    }
    else{
      if(fscanf(fp, "%d", &type) != 1){
	cout <<"bug in readings cell types" << endl;
	ExecError("error in reading vtk file");
      }
    }
    TypeCells[i] = type;
      switch(type){
      case 1:  // Vertex
	cout << "this type of cell is not taking account in Freefem++ " << endl;
	break;
      case 3:  // Edge/line
	cout << "this type of cell is not taking account in Freefem++ for a two dimensional mesh" << endl; // 3D
	break;
      case 5:  // Triangle
	nbe++; // 3D
	break;   
      case 10: // Tetrahèdre
	nt++;
	break;  
      default: 
	cout << "Error :: This type of cell is not considered in Freefem++"<< endl;
	ExecError("error in reading vtk file");
	break;
	
      }
    }  

    fclose(fp);

    
    // 3D versions
    
    Tet *tff;
    if( nt >0) tff = new Tet[nt];
    Tet *ttff = tff;

    Triangle3 *bff  = new Triangle3[nbe];
    Triangle3 *bbff = bff;

    for(unsigned int i = 0; i < numElements; i++){
      int type=TypeCells[i];
      int ivb[3],ivt[4];
      int label=1;
      switch(type){
      case 1:  // Vertex
	cout << "this type of cell is not taking account in Freefem++ " << endl;
	break;
      case 3:  // Edge/line
	break;
      case 5:  // Triangle
	cout << i << " " << firstCell[i+1] << " " << firstCell[i] << endl;
	assert( (firstCell[i+1]-firstCell[i]) == 3 );
	for(int j=firstCell[i]; j<firstCell[i+1]; j++){
    	  ivb[j-firstCell[i]] = IntCells[j];
	}
	(bbff++)->set(vff, ivb, label);
	break;   
      case 10: // Tetrahèdre
	assert( (firstCell[i+1]-firstCell[i]) == 4 );
	for(int j=firstCell[i]; j<firstCell[i+1]; j++){
    	  ivt[j-firstCell[i]] = IntCells[j];
	}
	(ttff++)->set(vff, ivt, label);
	break;  
      default: 
	break;	
      } 
    }
    
    delete [] IntCells;   
    delete [] firstCell; 
    delete [] TypeCells;

    if(nt>0){
      Mesh3 *pTh = new Mesh3(nv,nt,nbe,vff,tff,bff);
      return pTh;
    }
    else{
      Mesh3 *pTh = new Mesh3(nv,nbe,vff,bff);
      return pTh;
    }
    
}

AnyType VTK_LoadMesh3_Op::operator()(Stack stack)  const 
{
 
  string * pffname= GetAny<string *>((*filename)(stack));
  bool swap = false;
  int reftetra=1; 
  int reftri = 1;
  
  if( nargs[0] )  reftetra = GetAny<int >((*nargs[0])(stack));
  if( nargs[1] )  swap = GetAny< bool >((*nargs[1])(stack));
  if( nargs[2] )  reftri = GetAny<int >((*nargs[2])(stack));

  string *DataLabel;
  if( nargs[3] ) DataLabel = GetAny<string *>((*nargs[3])(stack));

  Mesh3 * Th = VTK_Load3( *pffname, swap); 

  // A faire fonction pour changer le label

  Add2StackOfPtr2FreeRC(stack,Th);
		
  return Th;
}

//==============================================
// ECRITURE DE FICHIER .vtk (3D)
//==============================================

class VTK_WriteMesh3_Op : public E_F0mps 
{
public:
  typedef long  Result;
  Expression eTh;
  Expression filename;
  
  struct Expression2 {
    string name;
    long what; // 1 scalar, 2 vector, 3 symtensor
    long nbfloat; // 1 scalar, 2 vector (3D), 3 symtensor(3D)
    Expression e[6];
    Expression2() {e[0]=0; e[1]=0; e[2]=0; e[3]=0; e[4]=0; e[5]=0; what=0; nbfloat=0;};
    Expression &operator[](int i){return e[i];}
    double eval(int i,Stack stack) const  { 
    if (e[i]) {
      return GetAny< double >( (*e[i])(stack) );
    }
    else 
      return 0;
    }
  };
  vector<Expression2> l;
  static const int n_name_param = 7;  
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  long    arg(int i,Stack stack,long a ) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
  
public:
  VTK_WriteMesh3_Op(const basicAC_F0 &  args) : l( args.size()-2 )
  {
    int nbofsol;
    int ddim=3;
    int stsize=6;

     int sca=0,vec=0,ten=0;
    
    string scas("scalaire"); 
    string vecs("vector"); 
    string tens("tensor"); 

    if(verbosity) cout << "Write Mesh and Solutions in VTK Formats" << endl;
    args.SetNameParam(n_name_param,name_param,nargs);

    if (BCastTo<string *>(args[0])) filename = CastTo<string *>(args[0]);
    if (BCastTo<pmesh3>(args[1])) eTh= CastTo<pmesh3>(args[1]);
    
    nbofsol = l.size();
    for (size_t i=2;i<args.size();i++){
      size_t jj=i-2;
      
      if ( BCastTo<double>( args[i] ))
	{
	  l[jj].what=1;
	  l[jj].nbfloat=1;
	  l[jj][0]=to<double>( args[i] );
	  
	  char number[16];
	  sprintf(number,"%li",jj+1);
	  l[jj].name=scas;
	  l[jj].name+=number;
	  sca++;

	}
      else if ( args[i].left()==atype<E_Array>() )
	{
	  const E_Array * a0  = dynamic_cast<const E_Array *>( args[i].LeftValue() );
	  //cout << "taille" << a0->size() << endl;
	  if (a0->size() != ddim && a0->size() != stsize) 
	    CompileError("savesol in 3D: vector solution is 3 composant, tensor solution is 6 composant");
	  
	  if( a0->size() == ddim){
	    // vector solution
	    l[jj].what=2;
	    l[jj].nbfloat=ddim;
	    for(int j=0; j<ddim; j++){
	      l[jj][j] = to<double>( (*a0)[j]);
	    }
	     char number[16];
	    sprintf(number,"%li",jj+1);
	    l[jj].name=vecs;
	    l[jj].name+=number;
	    vec++;	  
	    
	  }
	  else if( a0->size() == stsize){
	    // symmetric tensor solution
	    l[jj].what=3;
	    l[jj].nbfloat=stsize;
	    for(int j=0; j<stsize; j++){
	      l[jj][j] = to<double>( (*a0)[j]);
	    }
	    char number[16];
	    sprintf(number,"%li",jj+1);
	    l[jj].name=tens;
	    l[jj].name+=number;
	    ten++;
	  }
	  
	}
      else {
	cout << " arg " << i << " " << args[i].left() << endl;
	CompileError("savesol in 2D: Sorry no way to save this kind of data");
      }
    }
  }
  static ArrayOfaType  typeargs() { return  ArrayOfaType( atype<string *>(), atype<pmesh3>(), true); }// all type
  static  E_F0 * f(const basicAC_F0 & args) { return new VTK_WriteMesh3_Op(args);}   
  AnyType operator()(Stack stack)  const ;
};
basicAC_F0::name_and_type VTK_WriteMesh3_Op::name_param[]= {
  {  "dataname", &typeid(string*)},
  {  "withsurfacemesh", &typeid(bool)},
  {  "order", &typeid(KN_<long>)},
  // A rajouter dans le 3D
  {  "floatmesh", &typeid(bool)},
  {  "floatsol", &typeid(bool)},
  {  "bin", &typeid(bool)},
  {  "swap", &typeid(bool)}
};
  
void VTK_WRITE_MESH3( const string &filename, FILE *fp, const Mesh3 &Th, bool binary, int datasize, bool surface, bool bigEndian){
     
  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "%s, Created by Freefem++ \n", filename.c_str());
  if(binary)
    fprintf(fp, "BINARY\n");
  else
    fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
  // get all the entities in the model

  // write mesh vertices
  fprintf(fp, "POINTS %d double\n", Th.nv);
  
  if(datasize == sizeof(float)) {      
    for(unsigned int i = 0; i < Th.nv; i++){
      const  Vertex3 & P = Th.vertices[i];
      float f[3];
      f[0]=P.x;
      f[1]=P.y;
      f[2]=P.z; 
      if(binary){
	if(!bigEndian) SwapBytes((char*)&f, sizeof(float), 3);
	fwrite(&f, sizeof(float), 3, fp);
      }
      else{
	fprintf(fp,"%.8f %.8f %.8f\n",f[0],f[1],f[2]);
      }
    }  
  }
  else if(datasize == sizeof(double)){
    for(unsigned int i = 0; i < Th.nv; i++){
      const  Vertex3 & P = Th.vertices[i];
      double f[3];
      f[0]=P.x;
      f[1]=P.y;
      f[2]=P.z;    // 3D case
      if(binary){
	if(!bigEndian) SwapBytes((char*)&f, sizeof(float), 3);
	fwrite(&f, sizeof(float), 3, fp);
      }
      else{
	fprintf(fp,"%lf %lf %lf\n",f[0],f[1],f[2]);
      }
    }  
  }
  fprintf(fp,"\n");
  if(verbosity > 1) printf("writting vertices is finish\n");
  if(verbosity > 1) printf("writting elements now\n");
  
  //================
  //    CELL 
  //================
  // loop over all elements we need to save and count vertices
  int numElements, totalNumInt;
  if(surface){
    numElements = Th.nt+Th.nbe; 
    totalNumInt = Th.nt*4+Th.nbe*3+numElements;
  }
  else{
    numElements = Th.nt; 
    totalNumInt = Th.nt*4+numElements;
  }
  if(verbosity > 1) printf("writting cells \n");
  // print vertex indices in ascii or binary
  fprintf(fp, "CELLS %d %d\n", numElements, totalNumInt);

  if(binary){
    int IntType=4;
    if(verbosity > 1) printf("writting tetrahedron elements \n");
    for(int it=0; it< Th.nt; it++){    
      const Tet & K( Th.elements[it] );
      int iv[IntType+1];
      
      iv[0] = IntType;
      for(int ii=0; ii<IntType; ii++){
	iv[ii+1] = Th.operator()(K[ii]);
      }
      
      if(!bigEndian) SwapBytes((char*)&iv, sizeof(int), IntType+1);
      fwrite(&iv, sizeof(int), IntType+1, fp);
      
    }
    if(surface){
      if(verbosity > 1) printf("writting triangle elements \n");
      IntType=3;
      for(int ibe=0; ibe<Th.nbe; ibe++){
	 const Triangle3 &K( Th.be(ibe) );
      
	 int iv[IntType+1];    
	 iv[0] = IntType;
	 for(int ii=0; ii<IntType; ii++){
	   iv[ii+1] = Th.operator()(K[ii]);
	 }
	 
	 if(!bigEndian) SwapBytes((char*)&iv, sizeof(int), IntType+1);
	 fwrite(&iv, sizeof(int), IntType+1, fp);
      }
    }
  }
  else{
    int IntType=4;
    if(verbosity > 1) printf("writting tetrahedron elements \n");
    for(int it=0; it< Th.nt; it++){
      const Tet &K( Th.elements[it] );

      int iv[IntType+1];      
      iv[0] = IntType;
      for(int ii=0; ii<IntType; ii++){
	iv[ii+1] = Th.operator()(K[ii]);
      }
      fprintf(fp,"%d %d %d %d %d\n", iv[0],iv[1],iv[2],iv[3],iv[4]);
    }
    if(surface){
      if(verbosity > 1) printf("writting triangle elements \n");
      IntType=3;
      for(int ibe=0; ibe<Th.nbe; ibe++){
	const Triangle3 &K( Th.be(ibe) );
	
	int iv[IntType+1];    
	iv[0] = IntType;
	for(int ii=0; ii<IntType; ii++){
	  iv[ii+1] = Th.operator()(K[ii]);
	}
	
	fprintf(fp,"%d %d %d %d\n", iv[0],iv[1],iv[2],iv[3]);
      }
    }
  }
  fprintf(fp, "\n");


  // CELL_TYPE
  // print element types in ascii or binary  
  
  fprintf(fp, "CELL_TYPES %d\n", numElements); 
  if(binary){
    int type;    
    type = VTK_TET;  
    for(int it=0; it< Th.nt; it++){
     
      if(!bigEndian) SwapBytes((char*)&type, sizeof(int), 1);
      fwrite(&type, sizeof(int), 1, fp);
    }
    if(surface){
      type=VTK_TRI;
      for(int ibe=0; ibe<Th.nbe; ibe++){
	 
	if(!bigEndian) SwapBytes((char*)&type, sizeof(int), 1);
	fwrite(&type, sizeof(int), 1, fp);
      }
    }
  }
  else{
    int type;
    type= VTK_TET;
    for(int it=0; it< Th.nt; it++){      
      fprintf(fp,"%d ",type);
    }
    if(surface){
      type=VTK_TRI;
      for(int ibe=0; ibe<Th.nbe; ibe++){
	fprintf(fp,"%d ",type);
      }
    }
  } 

  fprintf(fp, "\n");

  //=================================
  //  WRITE SOLUTION IN FORMAT VTK
  //     LABEL OF ELEMENTS
  //=================================

  list<int> list_label_Elem;
  //list<int> list_label_Border_Elem;
  {
    
    for(int it=0; it< Th.nt; it++){
      const Tet &K( Th.elements[it] );
      list<int>::const_iterator ilist;
      int labOk=0;
      for( ilist=list_label_Elem.begin(); ilist!=list_label_Elem.end(); ilist++){
	  if( *ilist == K.lab ){ labOk = 1;   break; } 
	}
	if( labOk == 0){
	  list_label_Elem.push_back(K.lab);
	}
    }
    
    if(surface){
      for(int ibe=0; ibe<Th.nbe; ibe++){
	const Triangle3 &K( Th.be(ibe) );
	list<int>::const_iterator ilist;
	int labOk=0;
	for( ilist=list_label_Elem.begin(); ilist!=list_label_Elem.end(); ilist++){
	  if( *ilist == K.lab ){ labOk = 1;   break; } 
	}
	if( labOk == 0){
	  list_label_Elem.push_back(K.lab);
	}
      
      }
    }
  
  }
  list_label_Elem.sort();

  //=================================
  //=================================

  fprintf(fp, "CELL_DATA %d\n", numElements);
  int cell_fd=1;
  int cell_lab=1;
  fprintf(fp, "Scalars  Label int%d\n", cell_fd);
  fprintf(fp, "LOOKUP_TABLE FreeFempp_table\n"); 
  // Determination des labels
  if(binary){
    int label;    
    for(int it=0; it< Th.nt; it++){
      const Tet &K( Th.elements[it] );
      label =K.lab;
      if(!bigEndian) SwapBytes((char*)&label, sizeof(int), 1);
      fwrite(&label, sizeof(int), 1, fp);
    }
    if(surface){
      for(int ibe=0; ibe<Th.nbe; ibe++){
	const Triangle3 &K( Th.be(ibe) );
	label =K.lab;
	if(!bigEndian) SwapBytes((char*)&label, sizeof(int), 1);
	fwrite(&label, sizeof(int), 1, fp);
      }
    }
  }
  else{
    int label;
    for(int it=0; it< Th.nt; it++){ 
      const Tet &K( Th.elements[it] );
      label =K.lab;
      fprintf(fp,"%d\n",label);
    }
    if(surface){
      for(int ibe=0; ibe<Th.nbe; ibe++){
	const Triangle3 &K( Th.be(ibe) );
	label =K.lab;
	fprintf(fp,"%d\n",label);
      }
    }
  } 

  fprintf(fp,"\n");

  int size_list=0;
  list<int>::const_iterator ilist;
  for( ilist=list_label_Elem.begin(); ilist!=list_label_Elem.end(); ilist++) size_list++;

  fprintf(fp, "LOOKUP_TABLE FreeFempp_table %d\n",size_list);
  { list<int>::const_iterator ilist;
    for( ilist=list_label_Elem.begin(); ilist!=list_label_Elem.end(); ilist++){
      float tab[4];
      tab[0] = ColorTable[abs(*ilist)%NbColorTable][0];
      tab[1] = ColorTable[abs(*ilist)%NbColorTable][1];
      tab[2] = ColorTable[abs(*ilist)%NbColorTable][2];
      tab[3] = 1.0;
    
      if(binary){
	if(!bigEndian) SwapBytes((char*)&tab, sizeof(float), 4);
	fwrite(&tab, sizeof(float), 4, fp);
      }
      else
	fprintf(fp,"%.8f %.8f %.8f %.8f\n",tab[0],tab[1],tab[2],tab[3]);
    }
  }
  fprintf(fp,"\n");
}

AnyType VTK_WriteMesh3_Op::operator()(Stack stack)  const 
{
 
  string * pffname= GetAny<string *>((*filename)(stack));
  Mesh3 * pTh= GetAny<Mesh3 *>((*eTh)(stack));
  ffassert(pTh);
  Mesh3 &Th=*pTh;
  bool swap = false;
  bool bigEndian = false;
  bool binary = false;
  bool surface = true;
  bool floatmesh = true;
  bool floatsol = true;
  int datasize = sizeof(float);
  int datasizeSol = sizeof(float);
  string *dataname;
  int nbofsol = l.size();
  KN<int> order(nbofsol);
  
  char *nameofuser[nbofsol];

  for(int ii=0; ii<nbofsol; ii++) order[ii] = 0;

  if( nargs[0] ) dataname  = GetAny< string* >( (*nargs[0])(stack) );  
  if( nargs[1] ) surface   = GetAny< bool >( (*nargs[1])(stack) );
  if( nargs[2] ) order     = GetAny< KN_<long> >( (*nargs[2])(stack) );
  if( nargs[3] ) floatmesh = GetAny< bool >( (*nargs[3])(stack) );
  if( nargs[4] ) floatsol  = GetAny< bool  >( (*nargs[4])(stack) );
  if( nargs[5] ) binary    = GetAny< bool  >( (*nargs[5])(stack) );
  if( nargs[6] ) bigEndian = GetAny< bool  >( (*nargs[6])(stack) );

  swap = bigEndian;

  if( !floatmesh ) datasize = sizeof(double);
  if( !floatsol ) datasizeSol= sizeof(double);


  int iii=0;
  if( nargs[0]){
    char *data = newcopy(dataname->c_str());
    if(verbosity>5)
      cout << "   iovtk writeMesh3: names  \""<< data <<"\"" <<  endl;
    char * name =strtok(data," \n\0\t");
    nameofuser[iii] = newcopy(name);
    if(verbosity>5)
      cout << "   iovtk writeMesh3:value of iii=" << iii << " " << nameofuser[iii] <<endl;
    iii++;
    {
      
      while(( name= strtok(NULL," \n\0\t")) ){
	if( iii >= nbofsol )
	  {
	    if(verbosity)
	      cout << "   iovtk writeMesh3: The number of data name is too large " << endl;
	    break;
	  }
	nameofuser[iii] = newcopy(name);
	if(verbosity>5)
	  cout << "   iovtk writeMesh3:value of iii=" << iii << " " << nameofuser[iii] <<endl;
	iii++;
      }
      if( iii < nbofsol){	
	if(verbosity)
	cout << "   iovtk writeMesh3: The number of data name is too small, we give default name " << endl;
      }
      
    }
  }
  if( iii < nbofsol ){
    for( int iiii=iii; iiii<nbofsol; iiii++)
      nameofuser[iiii] = newcopy(l[iii].name.c_str());
    
  }


  // lecture du nom des variables

  FILE *fp = fopen( (*pffname).c_str(), "wb");
  if(!fp){
    cerr << "Unable to open file " << (*pffname).c_str() << endl;
    ExecError("error in reading vtk file")  ;
  }
  
  VTK_WRITE_MESH3( *pffname, fp, Th, binary, datasize, surface, swap); 
   	
  // determination of number of order 0 et 1.
  int Norder0=0;
  for(int ii=0; ii< nbofsol; ii++)
    if(order[ii] == 0) Norder0++;


  if( datasizeSol == sizeof(float) ){
    
    if( Norder0 >0){
      fprintf(fp, "FIELD FieldData %d\n", Norder0);
      for(int ii=0; ii< nbofsol; ii++){
	if(order[ii] == 0){
	  int nsol;
	  if(surface){
	    nsol = Th.nt+Th.nbe;
	  }
	  else{
	    nsol = Th.nt;
	  }
	  	  
	  fprintf(fp,"%s %ld %d float\n",nameofuser[ii],l[ii].nbfloat,nsol); 
	  if(verbosity>5)	    
	  cout << "   iovtk writeMesh3: name of data("<< ii <<")=" << nameofuser[ii]  << endl;
	
	  MeshPoint *mp3(MeshPointStack(stack)); 
	  R3 Cdg_hat = R3(1./4.,1./4.,1./4.);  
	
	  for (int it=0;it<Th.nt;it++)
	    {
	      const Tet  & K(Th.t(it));
	      mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);
	      
	      for(int j=0;j<l[ii].nbfloat;j++)
		{
		  float value = l[ii].eval(j,stack);
		  if(binary){
		    if(!bigEndian) SwapBytes((char*)&value, sizeof(float), 1);
		    fwrite(&value, sizeof(float), 1, fp);
		  }
		  else{
		    fprintf(fp,"%lf ",value); 
		  }
		}
	      if(!binary)  fprintf(fp,"\n");
	    }
	  if( surface ){
	    for (int ibe=0;ibe<Th.nbe;ibe++)
	      {
		// determination du triangle contenant cette edge
		int ie;
		int it = Th.BoundaryElement( ibe, ie); 
		const Tet  & K(Th.t(it));
		mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);
		
		for(int j=0;j<l[ii].nbfloat;j++){
		  float value = l[ii].eval(j,stack);
		  if(binary){
		    if(!bigEndian) SwapBytes((char*)&value, sizeof(float), 1);
		    fwrite(&value, sizeof(float), 1, fp);
		  }
		  else{
		    fprintf(fp,"%.8lf ",value); 
		  }
		} 
		if(binary)  fprintf(fp,"\n");
		
	      }
	  }
	  
	  fprintf(fp,"\n");
	}
      }
    }
    if( Norder0 < nbofsol )
      {
	fprintf(fp, "POINT_DATA %d\n", Th.nv);
	fprintf(fp, "FIELD FieldData %d\n", nbofsol-Norder0);
	for(int ii=0; ii< nbofsol; ii++){
	  if(order[ii] == 1)
	    {
	      
	      fprintf(fp,"%s %ld %d float\n",nameofuser[ii],l[ii].nbfloat,Th.nv); 
	      if(verbosity>5)
		cout << "   iovtk writeMesh3:name of data("<< ii <<")=" << nameofuser[ii]  << endl;
	    
	      
	      
	      MeshPoint *mp3(MeshPointStack(stack)); 
	      KN<double> valsol(Th.nv*l[ii].nbfloat);
	      KN<int> takemesh(Th.nv);
	      takemesh =0;
	      valsol   =0.;
	      for(int it=0;  it<Th.nt; it++)
		{
		  for(int iv=0; iv<4; iv++)
		    {
		      int i=Th(it,iv);
		      mp3->setP(&Th,it,iv);
		      
		      for(int j=0;j<l[ii].nbfloat;j++)
			{
			  valsol[ i*l[ii].nbfloat+j ] = valsol[ i*l[ii].nbfloat+j ] + l[ii].eval(j,stack);	     
			}
		      takemesh[i] = takemesh[i]+1;    
		    }
		}
	      
	      for (int iv=0;iv<Th.nv;iv++)
		{
		  for(int j=0;j<l[ii].nbfloat;j++)
		    {
		      valsol[iv*l[ii].nbfloat+j] = valsol[ iv*l[ii].nbfloat+j ]/takemesh[iv];
		      float value = valsol[iv*l[ii].nbfloat+j];
		      if(binary)
			{
			  if(!bigEndian) SwapBytes((char*)&value, sizeof(float), 1);
			  fwrite(&value, sizeof(float), 1, fp);
			}
		      else
			{
			  fprintf(fp,"%lf\n",value); 
			}
		    }	   
		}
	      if(!binary)
		fprintf(fp,"\n");
	    }
	}
      }    
  }
  
  // datasizeSol == sizeof(double)
  
  
  if( datasizeSol == sizeof(double) ){
    
    if( Norder0 >0){
      fprintf(fp, "FIELD FieldData %d\n", Norder0);
      for(int ii=0; ii< nbofsol; ii++){
	if(order[ii] == 0){
	  int nsol;
	  if(surface){
	    nsol = Th.nt+Th.nbe;
	  }
	  else{
	    nsol = Th.nt;
	  }
	  
	  fprintf(fp,"%s %ld %d float\n",nameofuser[ii],l[ii].nbfloat,nsol); 
	  if(verbosity>5)
	  cout << "   iovtk writeMesh3:name of data("<< ii <<")=" << nameofuser[ii]  << endl;
	    
	  MeshPoint *mp3(MeshPointStack(stack)); 
	  R3 Cdg_hat = R3(1./4.,1./4.,1./4.);  
	
	  for (int it=0;it<Th.nt;it++){
	    const Tet & K(Th.t(it));
	    mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);
	    
	    for(int j=0;j<l[ii].nbfloat;j++){
	      double value = l[ii].eval(j,stack);
	      if(binary){
		if(!bigEndian) SwapBytes((char*)&value, sizeof(double), 1);
		fwrite(&value, sizeof(double), 1, fp);
	      }
	      else{
		fprintf(fp,"%.8f ",value); 
	      }
	    } 
	  }
	  if( surface ){
	    for (int ibe=0;ibe<Th.nbe;ibe++){
	      // determination du triangle contenant cette edge
	      int ie;
	      int it = Th.BoundaryElement( ibe, ie); 
	      const Tet  & K(Th.t(it));
	      mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);
	      
	      for(int j=0;j<l[ii].nbfloat;j++){
		double value = l[ii].eval(j,stack);
		if(binary){
		  if(!bigEndian) SwapBytes((char*)&value, sizeof(double), 1);
		  fwrite(&value, sizeof(double), 1, fp);
		}
		else{
		  fprintf(fp,"%.8f ",value); 
		}
	      } 
	      
	    }
	  }
	  
	  fprintf(fp,"\n");
	}
      }
    }
    if( Norder0 < nbofsol ){
      fprintf(fp, "POINT_DATA %d\n", Th.nv);
      fprintf(fp, "FIELD FieldData %d\n", nbofsol-Norder0);
      for(int ii=0; ii< nbofsol; ii++){
	if(order[ii] == 1){

	  fprintf(fp,"%s %ld %d float\n",nameofuser[ii],l[ii].nbfloat,Th.nv); 
	  if(verbosity>5)
	  cout << "   iovtk writeMesh3:name of data("<< ii <<")=" << nameofuser[ii]  << endl;
	 
	  MeshPoint *mp3(MeshPointStack(stack)); 
	  KN<double> valsol(Th.nv*l[ii].nbfloat);
	  KN<int> takemesh(Th.nv);
	  takemesh =0;
	  valsol   =0.;
	  for(int it=0;  it<Th.nt; it++){
	    for(int iv=0; iv<4; iv++){
	      int i=Th(it,iv);
	      mp3->setP(&Th,it,iv);
	      
	      for(int j=0;j<l[ii].nbfloat;j++){
		valsol[ i*l[ii].nbfloat+j ] = valsol[ i*l[ii].nbfloat+j ] + l[ii].eval(j,stack);	     
	      }
	      takemesh[i] = takemesh[i]+1;    
	    }
	  }
	  
	  for (int iv=0;iv<Th.nv;iv++){
	    for(int j=0;j<l[ii].nbfloat;j++){
	      valsol[iv*l[ii].nbfloat+j] = valsol[ iv*l[ii].nbfloat+j ]/takemesh[iv];
	      double value = valsol[iv*l[ii].nbfloat+j];
	      if(binary){
		if(!bigEndian) SwapBytes((char*)&value, sizeof(double), 1);
		fwrite(&value, sizeof(double), 1, fp);
	      }
	      else{
		fprintf(fp,"%.8f ",value); 
	      }
	    }	   
	  }
	  fprintf(fp,"\n");
	}
      }
    }
  }

  fclose(fp);

  for( int iiii=0; iiii<nbofsol; iiii++){
    delete [] nameofuser[iiii];
  }


}

//==============================================
// FIN ECRITURE DE FICHIER .vtk (3D)
//==============================================

//=======================
// FIN 3D Fichier .vtk
//=======================
void saveMatlab(const string &file, const Mesh &Th)
{
  //   ErrorInsaveMesh e;
   {
     ofstream pf(file.c_str());
     ffassert(pf);
     typedef  Mesh::Element Element;
     for(int k=0; k<Th.nt;++k)
       {
	 Element &K = Th[k];
	 pf << "x = [ ";
	 for (size_t n=0; n<3; n++)
	   pf << std::setprecision(5) << setw(18) << K[n].x << " ";
	 pf << std::setprecision(5) << setw(18) << K[0].x << " ]; ";
	 pf << "y = [ ";
	 for (size_t n=0; n<3; n++)
	   pf << std::setprecision(5) << setw(18) << K[n].y << " ";
	 pf << std::setprecision(5) << setw(18) << K[0].y  << " ]; ";
	   pf << "line(x,y);" << endl;
       }
     pf.close();
   }
}

void saveTecplot(const string &file, const Mesh &Th)
{
  string shape;
  ofstream pf(file.c_str());
  size_t n, m;
  
  pf << "TITLE = \" \"\n";
  pf << "VARIABLES = \"X\", \"Y\"";
  if (Th.dim==3)
    pf << ", \"Z\"";
  pf << endl;
  if (Th.dim==2) {
    m = 3;
    shape = "TRIANGLE";
  }
  /*      else if (el->getShape()==LINE) {
	  m = 2;
	  shape = "LINESEG";
	  }
  */
  
  else if (Th.dim==3)
	 {
	   m = 4;
	   shape = "TETRAHEDRON";
	 }
  
  pf << "ZONE N=" << Th.nv  << ", E=" << Th.nt << ", F=FEPOINT, ET=" << shape << endl;
  for (int i=0;i<Th.nv;i++)
    pf << std::setprecision(5) << setw(18) << (R2) Th(i)  << " \n" ;
  
  
  for (int k=0;k<Th.nt;++k)
    {
      for (n=0; n<m; n++)
	pf << Th(k,n)+1  << "  ";
      pf << endl;
    }
  
  pf.close();
}


class Init1 { public:
  Init1();
};

static Init1 init1;  //  une variable globale qui serat construite  au chargement dynamique 

Init1::Init1(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 
  
  typedef Mesh *pmesh;
  //typedef Mesh2 *pmesh2;
  typedef Mesh3 *pmesh3;
  //if (verbosity)
  if(verbosity) cout << " load: iovtk " << endl;
  Global.Add("savevtk","(",new OneOperatorCode<VTK_WriteMesh_Op>);
  Global.Add("savevtk","(",new OneOperatorCode<VTK_WriteMesh3_Op>);
  Global.Add("vtkload","(",new VTK_LoadMesh3);
  Global.Add("vtkload","(",new VTK_LoadMesh);
  
}
