// ORIG-DATE:   September 2009
// -*- Mode : c++ -*%
//
// SUMMARY  : interface avec le logiciel gmsh    
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
//   
#include "ff++.hpp"

using namespace Fem2D;

// Table of number of vertex for an element type of gmsh 
static const int nvElemGmsh[30] = { 2, 3, 0, 4, 0, 
				    0, 0, 0, 0, 0, 
				    0, 0, 0, 0, 1, 
				    0, 0, 0, 0, 0, 
				    0, 0, 0, 0, 0, 
				    0, 0, 0, 0, 0 };   
// we considerer only edges, triangles and tetrahedrons in Freefem++
// 15 :: Vertex Corner
// 1  :: Edge/line
// 2  :: triangles
// 4  :: tetrahedrons
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


class GMSH_LoadMesh_Op : public E_F0mps 
{
public:
  Expression filename;
  static const int n_name_param = 2; // 
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
public:
  GMSH_LoadMesh_Op(const basicAC_F0 &  args,Expression ffname) 
    : filename(ffname)
  {
    if(verbosity) cout << "Load mesh given by GMSH " << endl;
    args.SetNameParam(n_name_param,name_param,nargs);
  } 
  
  AnyType operator()(Stack stack)  const ;
};
basicAC_F0::name_and_type GMSH_LoadMesh_Op::name_param[]= {
  {  "reftri", &typeid(long)},
  {  "renum", &typeid(long)}
};


class  GMSH_LoadMesh : public OneOperator { public:  
    GMSH_LoadMesh() : OneOperator(atype<pmesh>(),atype<string *>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
    return  new GMSH_LoadMesh_Op( args,t[0]->CastTo(args[0]) ); 
  }
};

Mesh * GMSH_Load(const string & filename)
  {
    // variable freefem++
    int nv, nt=0, nbe=0;
    Mesh::Vertex   *vff;
    map<int,int> mapnumv;

    // loading mesh and reading mesh in gmsh are in the file GModelIO_Mesh.cpp (directory Geo) 
    char str[256]  = "ZZZ";
    double version = 2.0;
    bool binary = false, swap = false, postpro = false;
    
    FILE *fp = fopen(filename.c_str(), "rb");
    if(!fp){
      cerr << "Unable to open file " << filename.c_str() << endl;
      exit(1);
    }
    
    while( !feof(fp) ){
      fgets(str, sizeof(str), fp);  
      if( str[0] == '$'){
	if(!strncmp(&str[1], "MeshFormat", 10)) {
	  if(!fgets(str, sizeof(str), fp)) exit(1);
	  int format, size;
	  if(sscanf(str, "%lf %d %d", &version, &format, &size) != 3) exit(1);
	  cout << "Mesh Format is " <<  format << endl;
	  if(format){
	    binary = true;
	    cout << "Mesh is in binary format" << endl;
	    int one;
	    if(fread(&one, sizeof(int), 1, fp) != 1) exit(1);
	    if(one != 1){
	      swap = true;
	      cout << "Swapping bytes from binary file" << endl;
	    }
	  }
	}
	else if(!strncmp(&str[1], "PhysicalNames", 13)) {	  
	  cout << " PhysicalNames is not considered in freefem++ " << endl;
	}
	
	else if(!strncmp(&str[1], "NO", 2) || !strncmp(&str[1], "Nodes", 5) ||
		!strncmp(&str[1], "ParametricNodes", 15)) {
	  
	  const bool parametric = !strncmp(&str[1], "ParametricNodes", 15);
	  if( parametric == true ){
	    cerr << " ParametricNodes is not considered yet in freefem++" << endl;
	    exit(1);
	  }
	  if(!fgets(str, sizeof(str), fp)) exit(1);
	  
	  if(sscanf(str, "%d", &nv) != 1) exit(1);
	  printf("%d vertices\n", nv);
	  
	  // local variables freefem++
	  vff = new Mesh::Vertex[nv]; 

	  int minVertex = nv + 1, maxVertex = -1;
	  for(int i = 0; i < nv; i++) {
	    int num;	
	    double xyz[3], uv[2];
	    
	    //if (!parametric){		
	    if(!binary){
	      if (fscanf(fp, "%d %lf %lf %lf", &num, &xyz[0], &xyz[1], &xyz[2]) != 4)
		exit(1);	  
	    }
	    else{
	      if(fread(&num, sizeof(int), 1, fp) != 1) exit(1);
	      if(swap) SwapBytes((char*)&num, sizeof(int), 1);
	      if(fread(xyz, sizeof(double), 3, fp) != 3) exit(1);
	      if(swap) SwapBytes((char*)xyz, sizeof(double), 3);
	    }
	    
	    assert( abs(xyz[2]) < 1.e-7);
	    vff[i].x = xyz[0];
	    vff[i].y = xyz[1];
	    vff[i].lab = 1;
	    mapnumv[num] = i;
	  }
	}
	else if(!strncmp(&str[1], "ELM", 3) || !strncmp(&str[1], "Elements", 8)) {

	  if(!fgets(str, sizeof(str), fp)) exit(1);
	  int numElements;
	  sscanf(str, "%d", &numElements);
	   
	  if(verbosity) cout << "Loading elements\n" << endl; 
	  if(!binary){
	    for(int i = 0; i < numElements; i++) {
	      int num, type, physical = 0, elementary = 0, partition = 0, numVertices;
	      if(version <= 1.0){
		fscanf(fp, "%d %d %d %d %d", &num, &type, &physical, &elementary, &numVertices);		
	      }
	      else{
		int numTags;
		fscanf(fp, "%d %d %d", &num, &type, &numTags);
		for(int j = 0; j < numTags; j++){
		  int tag;
		  fscanf(fp, "%d", &tag);       
		  if(j == 0)      physical = tag;
		  else if(j == 1) elementary = tag;
		  else if(j == 2) partition = tag;
		  // ignore any other tags for now
		}
		assert(type>=1 && type <=31);
		if( (numVertices = nvElemGmsh[type-1]) == 0){
		  cerr << "Element of type " << type  << " is not considered in Freefem++" << endl;
		  exit(1);
		}
	      }
	      
	      if( type == 1 ) nbe++;
	      if( type == 2 ) nt++;
	      if( type == 4 ){
		cout << "We are loading a two dimensionnal mesh " << endl;
		exit(1);
	      }
	      
	      int indices[60];
	      for(int j = 0; j < numVertices; j++) fscanf(fp, "%d", &indices[j]);
	      
	    }
	  }
	  else{
	    int numElementsPartial = 0;
	    while(numElementsPartial < numElements){
	      int header[3];
	      if( fread(header, sizeof(int), 3, fp) != 3 ) exit(1);
	      if(swap) SwapBytes((char*)header, sizeof(int), 3);
	      int type = header[0];
	      int numElms = header[1];
	      int numTags = header[2];
	      int numVertices; 
	      assert(type>=1 && type <=31);
	      if( (numVertices = nvElemGmsh[type-1]) == 0){
		cerr << "Element of type " << type  << " is not considered in Freefem++" << endl;
		exit(1);
	      }
	      unsigned int n = 1 + numTags + numVertices;
	      int *data = new int[n];

	      for(int i = 0; i < numElms; i++) {
		if(fread(data, sizeof(int), n, fp) != n) exit(1);
		if(swap) SwapBytes((char*)data, sizeof(int), n);
		int num = data[0];
		int physical = (numTags > 0) ? data[4 - numTags] : 0;
		int elementary = (numTags > 1) ? data[4 - numTags + 1] : 0;
		int partition = (numTags > 2) ? data[4 - numTags + 2] : 0;
		int *indices = &data[numTags + 1];
		
		if( type == 1 )  nbe++;
		if( type == 2 )  nt++;
		if( type == 4 ){
		  cout << "We are loading a two dimensionnal mesh " << endl;
		  exit(1);
		}
		
	      }
	      delete [] data;
	      numElementsPartial += numElms;
	      
	    }
	  }
	  break;
	}
      }
    }
    fclose(fp);
    
    Mesh::Triangle *tff  = new Mesh::Triangle[nt];
    Mesh::Triangle *ttff = tff;

    Mesh::BorderElement *bff  = new Mesh::BorderElement[nbe];
    Mesh::BorderElement *bbff = bff;

    // second reading
    fp = fopen(filename.c_str(), "rb");
    
    while( !feof(fp) ){
      fgets(str, sizeof(str), fp);  
      if( str[0] == '$'){
	
	if(!strncmp(&str[1], "ELM", 3) || !strncmp(&str[1], "Elements", 8)) {

	  if(!fgets(str, sizeof(str), fp)) exit(1);
	  int numElements;
	  sscanf(str, "%d", &numElements);
	   
	 
	  if(!binary){

	    int ie=0; 
	    int it=0;

	    for(int i = 0; i < numElements; i++) {
	      int num, type, physical = 0, elementary = 0, partition = 0, numVertices;
	      if(version <= 1.0){
		fscanf(fp, "%d %d %d %d %d", &num, &type, &physical, &elementary, &numVertices);		
	      }
	      else{
		int numTags;
		fscanf(fp, "%d %d %d", &num, &type, &numTags);
		for(int j = 0; j < numTags; j++){
		  int tag;
		  fscanf(fp, "%d", &tag);       
		  if(j == 0)      physical = tag;
		  else if(j == 1) elementary = tag;
		  else if(j == 2) partition = tag;
		  // ignore any other tags for now
		}
		assert(type>=1 && type <=31);
		if( (numVertices = nvElemGmsh[type-1]) == 0){
		  cerr << "Element of type " << type  << " is not considered in Freefem++" << endl;
		  exit(1);
		}
	      }

	      int indices[60];
	      for(int j = 0; j < numVertices; j++) fscanf(fp, "%d", &indices[j]);

	      if( type == 1 ){
		int iv0,iv1;
		iv0 = mapnumv[ indices[0] ];	
		iv1 = mapnumv[ indices[1] ];
		cout << "Elem " << ie+1 << " " << iv0+1 << " " << iv1+1 << endl;
		(bbff++)->set(vff, iv0, iv1, physical);
		ie++;
	      }
	      if( type == 2 ){
		int iv0,iv1,iv2;
		iv0 = mapnumv[ indices[0] ];	
		iv1 = mapnumv[ indices[1] ];
		iv2 = mapnumv[ indices[2] ];	
		cout << "Triangles " << it+1 << " " << iv0+1 << " " << iv1+1 << " " << iv2+1 << endl;
		
		(ttff++)->set(vff, iv0, iv1, iv2, physical);
		cout << "mes=" << tff[it].area << endl;
		if( tff[it].area < 1e-8 ){
		  cout << "bug" << endl;
		  exit(1);
		}
		it++;
	      }
	    }
	    assert(it==nt);
	    assert(ie==nbe);
	  }
	  else{
	    
	    int ie=0; 
	    int it=0;

	    int numElementsPartial = 0;
	    while(numElementsPartial < numElements){
	      int header[3];
	      if( fread(header, sizeof(int), 3, fp) != 3 ) exit(1);
	      if(swap) SwapBytes((char*)header, sizeof(int), 3);
	      int type = header[0];
	      int numElms = header[1];
	      int numTags = header[2];
	      int numVertices;
	      assert(type>=1 && type <=31);
	      if( (numVertices = nvElemGmsh[type-1]) == 0){
		cerr << "Element of type " << type  << " is not considered in Freefem++" << endl;
		exit(1);
	      }
	      unsigned int n = 1 + numTags + numVertices;
	      int *data = new int[n];
	      for(int i = 0; i < numElms; i++) {
		if(fread(data, sizeof(int), n, fp) != n) exit(1);
		if(swap) SwapBytes((char*)data, sizeof(int), n);
		int num = data[0];
		int physical = (numTags > 0) ? data[4 - numTags] : 0;
		int elementary = (numTags > 1) ? data[4 - numTags + 1] : 0;
		int partition = (numTags > 2) ? data[4 - numTags + 2] : 0;
		int *indices = &data[numTags + 1];
		
		if( type == 1 ){
		  int iv0,iv1;
		  iv0 = mapnumv[ indices[0] ];	
		  iv1 = mapnumv[ indices[1] ];
		  (bbff++)->set(vff, iv0, iv1, physical);
		  ie++;
		}
		if( type == 2 ){
		  double mes=-1;
		  int iv0,iv1,iv2;
		  iv0 = mapnumv[ indices[0] ];	
		  iv1 = mapnumv[ indices[1] ];
		  iv2 = mapnumv[ indices[2] ];		
		  (ttff++)->set(vff, iv0, iv1, iv2, physical,mes);
		  
		  it++;
		}
		
	      }
	      delete [] data;
	      numElementsPartial += numElms;
	  
	      
	    }
	    assert(it==nt);
	    assert(ie==nbe);
	  }
	}
        else if(!strncmp(&str[1], "NodeData", 8)) {
	  cout << " NodeData is not considered in freefem++ " << endl;
	}
	else if(!strncmp(&str[1], "ElementData", 11) ||
		!strncmp(&str[1], "ElementNodeData", 15)){

	  cout << " ElementData/ElementNodeData is not considered in freefem++ " << endl;
	}
      }
    }
    fclose(fp);

    Mesh *pTh = new Mesh(nv,nt,nbe,vff,tff,bff);
    return pTh;

  }



AnyType GMSH_LoadMesh_Op::operator()(Stack stack)  const 
{
 
  string * pffname= GetAny<string *>((*filename)(stack));
  int renumsurf = 0; 
  if( nargs[1] )  renumsurf = GetAny< int >( (*nargs[1])(stack) );
  assert( renumsurf <=1 && renumsurf >= 0);

  Mesh * Th = GMSH_Load( *pffname ); 
   
  Add2StackOfPtr2FreeRC(stack,Th);
		
  return Th;
}

// Load three dimensionnal mesh

class GMSH_LoadMesh3_Op : public E_F0mps 
{
public:
  Expression filename;
  static const int n_name_param = 2; // 
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
public:
  GMSH_LoadMesh3_Op(const basicAC_F0 &  args,Expression ffname) 
    : filename(ffname)
  {
    if(verbosity) cout << "Load mesh given by GMSH " << endl;
    args.SetNameParam(n_name_param,name_param,nargs);
  } 
  
  AnyType operator()(Stack stack)  const ;
};
basicAC_F0::name_and_type GMSH_LoadMesh3_Op::name_param[]= {
  {  "reftet", &typeid(long)},
  {  "renum", &typeid(long)}
};


class  GMSH_LoadMesh3 : public OneOperator { public:  
    GMSH_LoadMesh3() : OneOperator(atype<pmesh3>(),atype<string *>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
    return  new GMSH_LoadMesh3_Op( args,t[0]->CastTo(args[0]) ); 
  }
};

Mesh3 * GMSH_Load3(const string & filename)
  {
    // variable freefem++
    int nv, nt=0, nbe=0;
    Vertex3   *vff;
    map<int,int> mapnumv; 
	  
    // loading mesh and reading mesh in gmsh are in the file GModelIO_Mesh.cpp (directory Geo) 
    char str[256]  = "ZZZ";
    double version = 2.0;
    bool binary = false, swap = false, postpro = false;
    
    FILE *fp = fopen(filename.c_str(), "rb");
    if(!fp){
      cerr << "Unable to open file " << filename.c_str() << endl;
      exit(1);
    }
    
    while( !feof(fp) ){
      fgets(str, sizeof(str), fp);  
      if( str[0] == '$'){
	if(!strncmp(&str[1], "MeshFormat", 10)) {
	  if(!fgets(str, sizeof(str), fp)) exit(1);
	  int format, size;
	  if(sscanf(str, "%lf %d %d", &version, &format, &size) != 3) exit(1);
	  if(format){
	    binary = true;
	    cout << "Mesh is in binary format" << endl;
	    int one;
	    if(fread(&one, sizeof(int), 1, fp) != 1) exit(1);
	    if(one != 1){
	      swap = true;
	      cout << "Swapping bytes from binary file" << endl;
	    }
	  }
	}
	else if(!strncmp(&str[1], "PhysicalNames", 13)) {
	  cout << " PhysicalNames is not considered in freefem++ " << endl;
	}
	
	else if(!strncmp(&str[1], "NO", 2) || !strncmp(&str[1], "Nodes", 5) ||
		!strncmp(&str[1], "ParametricNodes", 15)) {
	  
	  const bool parametric = !strncmp(&str[1], "ParametricNodes", 15);
	  if( parametric == true ){
	    cerr << " ParametricNodes is not considered yet in freefem++" << endl;
	    exit(1);
	  }
	  if(!fgets(str, sizeof(str), fp)) exit(1);
	  
	  if(sscanf(str, "%d", &nv) != 1) exit(1);
	  printf("%d vertices\n", nv);
	  
	  // local variables freefem++
	  vff = new Vertex3[nv]; 
	  //map<int,int> mapnumv; 
	  
	  int minVertex = nv + 1, maxVertex = -1;
	  for(int i = 0; i < nv; i++) {
	    int num;	
	    double xyz[3], uv[2];
	    
	    //if (!parametric){		
	    if(!binary){
	      if (fscanf(fp, "%d %lf %lf %lf", &num, &xyz[0], &xyz[1], &xyz[2]) != 4)
		exit(1);	  
	    }
	    else{
	      if(fread(&num, sizeof(int), 1, fp) != 1) exit(1);
	      if(swap) SwapBytes((char*)&num, sizeof(int), 1);
	      if(fread(xyz, sizeof(double), 3, fp) != 3) exit(1);
	      if(swap) SwapBytes((char*)xyz, sizeof(double), 3);
	    }
	    
	    vff[i].x = xyz[0];
	    vff[i].y = xyz[1];
	    vff[i].z = xyz[2];
	    vff[i].lab = 1;
	    mapnumv[num] = i;
	  }
	}
	else if(!strncmp(&str[1], "ELM", 3) || !strncmp(&str[1], "Elements", 8)) {

	  if(!fgets(str, sizeof(str), fp)) exit(1);
	  int numElements;
	  sscanf(str, "%d", &numElements);
	   
	  if(!binary){
	    for(int i = 0; i < numElements; i++) 
	      {
		int num, type, physical = 0, elementary = 0, partition = 0, numVertices;
		if(version <= 1.0)
		  {
		    fscanf(fp, "%d %d %d %d %d", &num, &type, &physical, &elementary, &numVertices);		
		  }
		else{
		  int numTags;
		  fscanf(fp, "%d %d %d", &num, &type, &numTags);
		  for(int j = 0; j < numTags; j++){
		    int tag;
		    fscanf(fp, "%d", &tag);       
		    if(j == 0)      physical = tag;
		    else if(j == 1) elementary = tag;
		    else if(j == 2) partition = tag;
		    // ignore any other tags for now
		  }
		  assert(type>=1 && type <=31);
		  if( (numVertices = nvElemGmsh[type-1]) == 0){		  
		    cerr << "Element of type " << type  << " is not considered in Freefem++" << endl;
		    exit(1);
		  }
		}
		
		if( type == 1 ){
		  if(i==0)
		    cout << "edges in 3D mesh are not considered yet in freefem++, skeep data" << endl;	
		}
		if( type == 2 ) nbe++;
		if( type == 4 ) nt++;
		
		int indices[60];
		for(int j = 0; j < numVertices; j++) fscanf(fp, "%d", &indices[j]);
		
	      }
	  }
	  else
	    {
	      int numElementsPartial = 0;
	      while(numElementsPartial < numElements){
		int header[3];
		if( fread(header, sizeof(int), 3, fp) != 3 ) exit(1);
		if(swap) SwapBytes((char*)header, sizeof(int), 3);
		int type = header[0];
		int numElms = header[1];
		int numTags = header[2];
		int numVertices; 
		assert(type>=1 && type <=31);
		if( (numVertices = nvElemGmsh[type-1]) == 0){
		  cout << "Element of type " << type  << " is not considered in Freefem++" << endl;
		  exit(1);
	      }
		unsigned int n = 1 + numTags + numVertices;
		int *data = new int[n];
		
		for(int i = 0; i < numElms; i++) {
		  if(fread(data, sizeof(int), n, fp) != n) exit(1);
		  if(swap) SwapBytes((char*)data, sizeof(int), n);
		  int num = data[0];
		  int physical = (numTags > 0) ? data[4 - numTags] : 0;
		  int elementary = (numTags > 1) ? data[4 - numTags + 1] : 0;
		  int partition = (numTags > 2) ? data[4 - numTags + 2] : 0;
		  int *indices = &data[numTags + 1];
		  
		  if( type == 1 && i==0 ){
		    cout << "edges in 3D mesh are not used in freefem++,skeep data" << endl;
		    //exit(1);
		  }
		  if( type == 2 ) nbe++;
		  if( type == 4 ) nt++;
		}
		delete [] data;
		numElementsPartial += numElms;
		
	    }
	  }
	  break;
	}
      }
    }
    fclose(fp);
    
    cout << "closing file " << nt << " " << nbe << endl;

    Tet *tff  = new Tet[nt];
    Tet *ttff = tff;

    Triangle3 *bff  = new Triangle3[nbe];
    Triangle3 *bbff = bff;

    // second reading
    fp = fopen(filename.c_str(), "rb");
    
    while( !feof(fp) ){
      fgets(str, sizeof(str), fp);  
      if( str[0] == '$'){
	
	if(!strncmp(&str[1], "ELM", 3) || !strncmp(&str[1], "Elements", 8)) {

	  if(!fgets(str, sizeof(str), fp)) exit(1);
	  int numElements;
	  sscanf(str, "%d", &numElements);

	  printf("%d tetrahedrons\n", nt);
	  printf("%d triangles\n", nbe);
	  printf("%d numElements\n", numElements);
	  if(!binary){

	    int ie=0; 
	    int it=0;

	    for(int i = 0; i < numElements; i++) {
	      int num, type, physical = 0, elementary = 0, partition = 0, numVertices;
	      if(version <= 1.0){
		fscanf(fp, "%d %d %d %d %d", &num, &type, &physical, &elementary, &numVertices);		
	      }
	      else{
		int numTags;
		fscanf(fp, "%d %d %d", &num, &type, &numTags);
		for(int j = 0; j < numTags; j++){
		  int tag;
		  fscanf(fp, "%d", &tag);  
	 
		  if(j == 0)      physical = tag;
		  else if(j == 1) elementary = tag;
		  else if(j == 2) partition = tag;
		  // ignore any other tags for now
		}
		assert(type>=1 && type <=31);
		if( (numVertices = nvElemGmsh[type-1]) == 0){
		  cerr << "Element of type " << type  << " is not considered in Freefem++" << endl;
		  exit(1);
		}
	      }

	      int indices[60];
	      for(int j = 0; j < numVertices; j++){
		fscanf(fp, "%d", &indices[j]);
	      }
	      if( type == 2 ){
		int ivff[3];
		for(int ii=0; ii < numVertices; ii++){
		  ivff[ii] = mapnumv[ indices[ii] ];
		  assert( ivff[ii]>=0 && ivff[ii]< nv );
		}
		(bbff++)->set(vff,ivff,physical);
		ie++;
	      }
	      if( type == 4 ){
		int ivff[4];
		for(int ii=0; ii < numVertices; ii++){
		  ivff[ii] = mapnumv[ indices[ii] ];
		  assert( ivff[ii]>=0 && ivff[ii]< nv );
		}
		(ttff++)->set(vff,ivff,physical);
		it++;
	      }
	    }
	    assert( it==nt );
	    assert( ie==nbe );
	  }
	  else{
	    
	    int ie=0; 
	    int it=0;

	    int numElementsPartial = 0;
	    while(numElementsPartial < numElements){
	      int header[3];
	      if( fread(header, sizeof(int), 3, fp) != 3 ) exit(1);
	      if(swap) SwapBytes((char*)header, sizeof(int), 3);
	      int type = header[0];
	      int numElms = header[1];
	      int numTags = header[2];
	      int numVertices;
	      assert(type>=1 && type <=31);
	      if( (numVertices = nvElemGmsh[type-1]) == 0){
		cerr << "Element of type " << type  << " is not considered in Freefem++" << endl;
		exit(1);
	      }
	      unsigned int n = 1 + numTags + numVertices;
	      int *data = new int[n];
	      for(int i = 0; i < numElms; i++) {
		if(fread(data, sizeof(int), n, fp) != n) exit(1);
		if(swap) SwapBytes((char*)data, sizeof(int), n);
		int num = data[0];
		int physical = (numTags > 0) ? data[4 - numTags] : 0;
		int elementary = (numTags > 1) ? data[4 - numTags + 1] : 0;
		int partition = (numTags > 2) ? data[4 - numTags + 2] : 0;
		int *indices = &data[numTags + 1];

		if( type == 2 ){
		  int ivff[3];
		  for(int ii=0; ii < numVertices; ii++) ivff[ii] = mapnumv[ indices[ii] ];		
		  (bbff++)->set(vff,ivff,physical);
		  ie++;
		}
		if( type == 4 ){
		  int ivff[4];
		  for(int ii=0; ii < numVertices; ii++) ivff[ii] = mapnumv[ indices[ii] ];		
		  (ttff++)->set(vff,ivff,physical);
		  it++;
		}
		
	      }
	      delete [] data;
	      numElementsPartial += numElms;
	     
	    }
	    assert(it==nt);
	    assert(ie==nbe);
	  }
	}
        else if(!strncmp(&str[1], "NodeData", 8)) {
	  cout << " NodeData is not considered in freefem++ " << endl;
	}
	else if(!strncmp(&str[1], "ElementData", 11) ||
		!strncmp(&str[1], "ElementNodeData", 15)){

	  cout << " ElementData/ElementNodeData is not considered in freefem++ " << endl;
	}
      }
    }
    fclose(fp);

    if(nt==0){
      Mesh3 *Th3 = new Mesh3(nv,nbe,vff,bff);
      return Th3;
    }
    else{
      Mesh3 *Th3 = new Mesh3(nv,nt,nbe,vff,tff,bff);
      return Th3;  
    }
    

  }



AnyType GMSH_LoadMesh3_Op::operator()(Stack stack)  const 
{
 
  string * pffname= GetAny<string *>((*filename)(stack));
  int renumsurf = 0; 
  if( nargs[1] )  renumsurf = GetAny< int >( (*nargs[1])(stack) );
  assert( renumsurf <=1 && renumsurf >= 0);

  Mesh3 * Th3_t = GMSH_Load3( *pffname); 
 
  Th3_t->BuildGTree();  
  Add2StackOfPtr2FreeRC(stack,Th3_t);
		
  return Th3_t;
}

class Init1 { public:
  Init1();
};

static Init1 init1;  //  une variable globale qui serat construite  au chargement dynamique 

Init1::Init1(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 
  
  //if (verbosity)
  if(verbosity) cout << " load: gmsh " << endl;
  Global.Add("gmshload3","(",new GMSH_LoadMesh3);
  Global.Add("gmshload","(",new GMSH_LoadMesh);
  if(verbosity) cout << " load: gmsh  " << endl;
}
