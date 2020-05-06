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
// SUMMARY : READ/WRITE MESH IN FORMAT PLY - Polygon File Format / Stanford Triangle Format
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS :
// E-MAIL  :

#include "ff++.hpp"

using namespace Fem2D;

// =====================
// Test for BigEndian
// =====================
bool isBigEndian( ) {
    unsigned int x = 1;
    char *ptr = (char*)&x;
    if (ptr[0] == 1) {
        if (verbosity>1) cout << "machine is little endian" << endl;
        return false;
    }
    else {
        if (verbosity>1)  cout << "machine is big endian" << endl;
        return true;
    }
}

namespace FreeFEM {
void SwapBytes(char *array, int size, int n) {
    char *x = new char[size];
    
    for (int i = 0; i < n; i++) {
        char *a = &array[i * size];
        memcpy(x, a, size);
        
        for (int c = 0; c < size; c++)
            a[size - 1 - c] = x[c];
    }
    delete[] x;
}
}



//= =============================================
// LOAD DE FICHIER .ply for meshS
//= =============================================

template <class MMesh>
class PLY_LoadMeshT_Op : public E_F0mps {
public:
    Expression filename;
    static const int n_name_param = 4;    //
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
    
public:
    PLY_LoadMeshT_Op(const basicAC_F0 &args, Expression ffname) : filename(ffname) {
        if (verbosity) {
            cout << "Load mesh given by PLY " << endl;
        }
        
        args.SetNameParam(n_name_param, name_param, nargs);
    }
    
    AnyType operator( )(Stack stack) const;
};

template<>
basicAC_F0::name_and_type PLY_LoadMeshT_Op< Mesh3 >::name_param[] = {
    {"swap", &typeid(bool)}, {"cleanmesh", &typeid(bool)},
    {"removeduplicate", &typeid(bool)}, {"precisvertice", &typeid(double)}};
template<>
basicAC_F0::name_and_type PLY_LoadMeshT_Op< MeshS >::name_param[] = {
    {"swap", &typeid(bool)}, {"cleanmesh", &typeid(bool)},
    {"removeduplicate", &typeid(bool)}, {"precisvertice", &typeid(double)}};
template<>
basicAC_F0::name_and_type PLY_LoadMeshT_Op< MeshL >::name_param[] = {
    {"swap", &typeid(bool)}, {"cleanmesh", &typeid(bool)},
    {"removeduplicate", &typeid(bool)}, {"precisvertice", &typeid(double)}};



template< class MMesh >
class PLY_LoadMeshT : public OneOperator {
public:
    typedef const MMesh *ppmesh;
    PLY_LoadMeshT( ) : OneOperator(atype< ppmesh >( ), atype< string * >( )) {}
    
    E_F0 *code(const basicAC_F0 &args) const {
        return new PLY_LoadMeshT_Op< MMesh >(args, t[0]->CastTo(args[0]));
    }
};

template< class MMesh >
MMesh *PLY_LoadT(const string &filename, bool bigEndian, bool cleanmesh, bool removeduplicate,
                 double precisvertice) {
    
    typedef typename MMesh::Element T;
    typedef typename MMesh::BorderElement B;
    typedef typename MMesh::Vertex V;

    int nv, nt = 0, nbe = 0;
    int nerr = 0;
    char *res;
    char buffer[256], buffer2[256], buffer3[256];
    
    FILE *fp = fopen(filename.c_str( ), "rb");
    if (!fp) {
        cerr << "Unable to open file " << filename.c_str( ) << endl;
        ExecError("Unable to open ply file");
    }
    fscanf(fp, "%s", buffer); // file identifiant
    if (strcmp(buffer, "ply")){
        cout << "Invalid ply file id" << endl;
        ExecError("Invalid ply file id");
    }
    
    // variable to check info given
    bool havelabVert = false, havelabElem = false, havelabBdElem = false; //label
    bool havecolorVert = false, havecolorElem = false, havecolorBdElem = false; // color in RGB
    bool haveElem = false, haveBdElem = false, indexVertFace = false;  // vertice index given by a uchar
    bool haveNormal = false, haveMaterialIndex = false;// normals, material property, amount of transparency
    bool haveTransparencyVert = false, haveTransparencyElem = false, haveTransparencyBdElem = false;
    bool end_header = false;
    int datasize=0; // properties of coordinates
    
    
    double version;
    if (fscanf(fp, "%s %s %lf", buffer, buffer2, &version) != 3) {
        cout << "error in reading header ply file" << endl;
        ExecError("error in reading header ply file");
    }
    
    if (version > 1.){
        cout << "only PLY format version 1.0 supported" << endl;
        ExecError("only PLY format version 1.0 supported");
    }
    
    // data formats: ascii / binary_little_endian / binary_big_endian
    // swap = bigEndian or not bigEndian
    bool binary = false, swap=false;
    
    if (!strcmp(buffer2, "binary_little_endian")) {
        binary = true;
        swap = bigEndian ? 1 : 0;
        
    }
    else if (!strcmp(buffer2, "binary_big_endian")) {
        binary = true;
        swap = bigEndian ? 0 : 1;
    }
    
    // read header
    while (!end_header) {
        
        fscanf(fp, "%s", buffer);
        if (!strcmp(buffer, "end_header")) end_header = true;
        
        // skip comment: information about the data
        // element description
        if (!strcmp(buffer, "element")) {
            fscanf(fp, "%s", buffer);
            if (!strcmp(buffer, "vertex"))
                fscanf(fp, "%d", &nv);
            else if (!strcmp(buffer, "face")) {
                if(!haveElem) {
                    fscanf(fp, "%d", &nt);
                    haveElem=true;}
                else {
                    fscanf(fp, "%d", &nbe);
                    haveBdElem=true;}
                
            }
            else if (!strcmp(buffer, "edge")) {
                if (haveBdElem) {
                    cout << "error in format border elements ply file" << endl;
                    ExecError("error in format border elements ply file");}
                else {
                    fscanf(fp, "%d", &nbe);
                    haveBdElem=true;}
            }
        }
            // properties about element section
            if (!strcmp(buffer, "property")) {
                fscanf(fp, "%s", buffer);
                if (!strcmp(buffer, "float") || !strcmp(buffer, "float16")) {
                    fscanf(fp, "%s", buffer2);
                    if (!strcmp(buffer2, "x") || !strcmp(buffer, "y") || !strcmp(buffer, "z") ) {
                        if(!datasize) datasize = sizeof(float);
                        else if(datasize!=sizeof(float) ) {
                            cout << "error in format vertex ply file" << endl;
                            ExecError("error in reading vertex  ply file");
                        }
                    }
                    else if (!strcmp(buffer2, "nx") || !strcmp(buffer, "ny") || !strcmp(buffer, "nz") )
                        haveNormal = true;
                }
                else if (!strcmp(buffer, "double") || !strcmp(buffer, "float32")) {
                    fscanf(fp, "%s", buffer2);
                    if (!strcmp(buffer2, "x") || !strcmp(buffer, "y") || !strcmp(buffer, "z") ){
                        if(!datasize) datasize = sizeof(double);
                        else if(datasize!=sizeof(double) ) {
                            cout << "error in format vertex ply file" << endl;
                            ExecError("error in reading vertex  ply file");
                        }
                    }
                    else if (!strcmp(buffer2, "nx") || !strcmp(buffer, "ny") || !strcmp(buffer, "nz") )
                        haveNormal = true;
                }
                // description about the list provides by element
                else if (!strcmp(buffer, "list")) {
                    fscanf(fp, "%s %s %s", buffer,buffer2, buffer3);
                    if (!strcmp(buffer3, "vertex_indices")) indexVertFace=true; //havelabel=true;  //number of vertices for each face
                }
                else if (!strcmp(buffer, "int")) {
                    fscanf(fp, "%s", buffer);
                    if (!strcmp(buffer, "flags")) {
                        if (nv && !nt && !nbe)
                            havelabVert=true;
                        else if( nv && nt && !nbe)
                            havelabElem=true;
                        else if( nv && nt && nbe)
                            havelabBdElem=true;
                    }
                }
                
                // color given in section element
                else if (!strcmp(buffer, "uint8") || !strcmp(buffer, "int") || !strcmp(buffer, "uchar") ) {
                    fscanf(fp, "%s", buffer);
                    if (!strcmp(buffer, "red") || !strcmp(buffer, "green") || !strcmp(buffer, "blue")) {
                        if (nv && !nt && !nbe)
                            havecolorVert=true;
                        else if( nv && nt && !nbe)
                            havecolorElem=true;
                        else if( nv && nt && nbe)
                            havecolorBdElem=true;
                    }
                    else if (!strcmp(buffer, "alpha") ) {
                        if (nv && !nt && !nbe)
                            haveTransparencyVert=true;
                        else if( nv && nt && !nbe)
                            haveTransparencyElem=true;
                        else if( nv && nt && nbe)
                            haveTransparencyBdElem=true;
                    }
                }
            }
        }
        
        if (verbosity >9) {
            cout << " ***** info about ply mesh ***** " << endl;
            cout << " havelabVert= " << havelabVert <<  " havelabElem= " << havelabElem <<  endl;
            cout << " havelabBdElem= " << havelabBdElem << endl;
            cout << " havecolorVert= " << havecolorVert <<  " havecolorElem= " << havecolorElem <<  " havecolorBdElem= " << havecolorBdElem << endl;
            cout << " haveNormal= " << haveNormal << " haveElem= " << haveElem << " haveBdElem= " << haveBdElem << endl;
            cout << " haveMaterialIndex= " << haveMaterialIndex << endl;
        }
        
        V *vff = new V[nv];
        T *tff = new T[nt];
        T *ttff = tff;
        B *bff = new B[nbe];
        B *bbff = bff;
        
        // read vertex
        if (verbosity > 3)
            cout << "Reading " << nv << " points, " << nt << " elements, " << nbe  << " border element " << endl;
        
        int numVerts=0;
        
        if (binary) {
            int hack;
            fscanf(fp, "%d", &hack);  // hack
        }
        
        for (int i = 0; i < nv; i++) {
            if (verbosity > 9)
                cout << " i=" << i << endl;
            double xyz[3], n[3];
            unsigned int color[3];
            if (binary) {
                if (datasize == sizeof(float)) {
                    float f[3], nn[3];
                    
                    if (fread(f, sizeof(float), 3, fp) != 3) {
                        cout << "error in reading PLY file" << endl;
                        ExecError("error in reading PLY file");
                    }
                    if (swap)
                        FreeFEM::SwapBytes((char *)f, sizeof(float), 3);
                    
                    if (verbosity>5)
                        printf("-- xyz %d = %lf %lf %lf\n", i , f[0], f[1], f[2] );
                    for (int j = 0; j < 3; j++)
                        xyz[j] = f[j];
                    
                    // don't need
                    if (haveNormal) {
                        if (fread(nn, sizeof(float), 3, fp) != 3) {
                            cout << "error in reading PLY file" << endl;
                            ExecError("error in reading PLY file");
                        }
                    }
                    
                }
                else {
                    if (fread(xyz, sizeof(double), 3, fp) != 3) {
                        cout << "error in reading PLY file" << endl;
                        ExecError("error in reading PLY file");
                    }
                    if (swap)
                        FreeFEM::SwapBytes((char *)xyz, sizeof(double), 3);
                    
                    // don't need
                    if (haveNormal) {
                        if (fread(n, sizeof(float), 3, fp) != 3) {
                            cout << "error in reading PLY file" << endl;
                            ExecError("error in reading PLY file");
                        }
                    }
                }
                
                // don't need
                if (havecolorVert) {
                    int color2[3];
                    if (fread(color2, sizeof(char), 3, fp) != 3) {
                        cout << "error in reading PLY file" << endl;
                        ExecError("error in reading PLY file");
                    }
                }
                // don't need
                if (haveTransparencyVert) {
                    char transparencyVert;
                    if (fread(&transparencyVert, sizeof(char), 1, fp) != 1) {
                        cout << "error in reading PLY file" << endl;
                        ExecError("error in reading PLY file");
                    }
                }
            }
            else {
                if (fscanf(fp, "%lf %lf %lf", &xyz[0], &xyz[1], &xyz[2]) != 3) {
                    cout << "error in reading PLY file (float)" << endl;
                    ExecError("error in reading PLY file (float)");
                }
                if (haveNormal) {
                    if (fscanf(fp, "%lf %lf %lf\n", &n[0], &n[1], &n[2]) != 3) {
                        cout << "error in reading PLY file (normal)" << endl;
                        ExecError("error in reading PLY file");
                    }
                }
                if (havecolorVert) {
                    if (fscanf(fp, "%u %u %u\n", &color[0], &color[1], &color[2]) != 3) {
                        cout << "error in reading ply file (color vert)" << endl;
                        ExecError("error in reading ply file (color vert");
                    }
                }
                if (haveTransparencyVert) {
                    int transparencyVert;
                    if (fscanf(fp, "%d\n", &transparencyVert) != 1) {
                        cout << "error in reading ply file (alpha)" << endl;
                        ExecError("error in reading ply file (alpha)");
                    }
                }
            }
            int lab=1;
            vff[i].x = xyz[0];
            vff[i].y = xyz[1];
            vff[i].z = xyz[2];
            vff[i].lab = lab;
            if (verbosity > 9)
                printf("xyz %d = %f %f %f\n", i , xyz[0], xyz[1], xyz[2]);
        }
        
        int label[nt], ivt[T::nt], ivb[B::nv], color[3];
        
        for (int it = 0; it < nt; ++it) {
            label[it]=0;
            if (verbosity > 9)
                cout << "it=" << it << " " << nt << endl;
            
            if (binary) {
                if (indexVertFace && (fread(&numVerts, sizeof(char), 1, fp) != 1)) {
                    cout << "error in reading PLY file (numVerts element)" << endl;
                    ExecError("error in reading PLY file (numVerts element)");
                }
                if (swap)
                    FreeFEM::SwapBytes((char *)&numVerts, sizeof(int), 1);
                if (verbosity > 9)
                    cout << "numVerts= " << numVerts << endl;
                
                if (fread(ivt, sizeof(int), numVerts, fp) != numVerts) {
                    cout << "error in reading PLY file (element)" << endl;
                    ExecError("error in reading PLY file (element)");
                }
                if(verbosity>5 && T::nv==4)
                    cout << "test " << ivt[0] << " " <<ivt[1] << " "<<ivt[2] << " "<<ivt[3] << endl;
                else if(verbosity>5 && T::nv==3)
                cout << "test " << ivt[0] << " " <<ivt[1] << " "<<ivt[2] << endl;
                
                if (swap)
                    FreeFEM::SwapBytes((char *)ivt, sizeof(int), numVerts);
                
                if (havelabElem) {
                    
                    if( fread(&label[it], sizeof(int), 1, fp) != 1) {
                        cout << "error in reading PLY file (lab element) " << endl;
                        ExecError("error in reading PLY file (lab element)");
                    }
                    if (swap)
                        FreeFEM::SwapBytes((char *)&label[it], sizeof(int), 1);
                }
                // don't need
                if (havecolorElem) {
                    char color2[3];
                    if( fread(color2, sizeof(char), 3, fp) != 3) {
                        cout << "error in reading PLY file (color element)" << endl;
                        ExecError("error in reading PLY file (color element)");
                    }
                }
                // don't need
                if (haveTransparencyVert) {
                    char transparencyElem;
                    if (fread(&transparencyElem, sizeof(char), 1, fp) != 1) {
                        cout << "error in reading PLY file (transparency element)" << endl;
                        ExecError("error in reading PLY file (transparency element)");
                    }
                }
            }
            
            else {
                if (indexVertFace && (fscanf(fp, "%d", &numVerts) != 1) ) {
                    
                    cout << "error in reading PLY file (numVerts element)" << endl;
                    ExecError("error in reading PLY file (numVerts element)");
                }
                if(verbosity>5)
                    cout << " numVerts test " << numVerts << endl;
                ffassert(numVerts==T::nv);   /////// pb
                
                for (int j = 0; j < numVerts; ++j) {
                    if (fscanf(fp, "%d", &ivt[j]) != 1) {
                        cout << "error in reading PLY file (element)" << endl;
                        ExecError("error in reading PLY file (element)'");
                    }
                    if (verbosity > 9)
                        cout << "  ivt[j]" << ivt[j] << endl;
                }
                if (havelabElem && (fscanf(fp, "%d", &label[it]) != 1) ){
                    cout << "error in reading PLY file (lab element)" << endl;
                    ExecError("error in reading PLY file (lab element)");
                }
                
                if (havelabElem && verbosity > 9)
                    cout << "  label[it]" << label[it] << endl;
                if (havecolorElem)
                    fscanf(fp, "%d %d %d", &color[0], &color[1], &color[2]);
                
                if (haveTransparencyElem) {
                    int transparencyElem;
                    fscanf(fp, "%d", &transparencyElem);
                }
            }
            
            (ttff++)->set(vff, ivt, label[it]);
        }
        
        if(B::nv==2) {
            for (int ibe = 0; ibe < nbe; ++ibe) {
                
                if(binary) {
                    if (fread(ivb, sizeof(int), 2, fp) != 2) {
                        cout << "error in reading PLY file (edge)" << endl;
                        ExecError("error in reading PLY file (edge)");
                    }
                    if (swap)
                        FreeFEM::SwapBytes((char *)ivb, sizeof(int), 2);
                    // don't need
                    if (havecolorBdElem && (fread(color, sizeof(int), 3, fp) != 3) ) {
                        cout << "error in reading PLY file (edge)" << endl;
                        ExecError("error in reading ply file (edge)");
                    }
                    // don't need
                    if (haveTransparencyVert) {
                        int transparencyBdElem;
                        if (fread(&transparencyBdElem, sizeof(int), 1, fp) != 1) {
                            cout << "error in reading PLY file" << endl;
                            ExecError("error in reading PLY file");
                        }
                    }
                }
                else {
                    if (fscanf(fp, "%d %d", &ivb[0], &ivb[1]) != 2) {
                        cout << "error in reading edge PLY files " << endl;
                        ExecError("error in reading edge PLY file");
                    }
                    if (verbosity>5)
                        cout << "test ivb " <<ivb[0] << " " <<ivb[1]  << endl;
                    
                    if (havecolorBdElem)
                        fscanf(fp, "%d %d %d", &color[0], &color[1], &color[2]);
                    
                    if (haveTransparencyBdElem ) {
                        int transparencyBdElem;
                        fscanf(fp, "%d", &transparencyBdElem);
                    }
                }
                int lab = 1;
                (bbff++)->set(vff, ivb, lab);
            }
        }
        else if(B::nv==3) {
            
            for (int ibe = 0; ibe < nbe; ++ibe) {
                label[ibe]=0;
                if (verbosity > 9)
                    cout << "ibe=" << ibe << " " << nbe << endl;
                
                if (binary) {
                    if (indexVertFace && (fread(&numVerts, sizeof(char), 1, fp) != 1)) {
                        cout << "error in reading PLY file (numVerts element)" << endl;
                        ExecError("error in reading PLY file (numVerts element)");
                    }
                    
                    if (swap)
                        FreeFEM::SwapBytes((char *)&numVerts, sizeof(int), 1);
                    if (verbosity > 9)
                        cout << "numVerts= " << numVerts << endl;
                    
                    if (fread(ivb, sizeof(int), numVerts, fp) != numVerts) {
                        cout << "error in reading PLY file (element)" << endl;
                        ExecError("error in reading PLY file (element)");
                    }
                    if(verbosity>5)
                        cout << "test " << ivb[0] << " " <<ivb[1] << " "<<ivb[2] << endl;
                    
                    if (swap)
                        FreeFEM::SwapBytes((char *)ivt, sizeof(int), numVerts);
                    
                    if (havelabElem) {
                        
                        if( fread(&label[ibe], sizeof(int), 1, fp) != 1) {
                            cout << "error in reading PLY file (lab element) " << endl;
                            ExecError("error in reading PLY file (lab element)");
                        }
                        if (swap)
                            FreeFEM::SwapBytes((char *)&label[ibe], sizeof(int), 1);
                    }
                    // don't need
                    if (havecolorElem) {
                        char color2[3];
                        if( fread(color2, sizeof(char), 3, fp) != 3) {
                            cout << "error in reading PLY file (color element)" << endl;
                            ExecError("error in reading PLY file (color element)");
                        }
                    }
                    // don't need
                    if (haveTransparencyVert) {
                        char transparencyElem;
                        if (fread(&transparencyElem, sizeof(char), 1, fp) != 1) {
                            cout << "error in reading PLY file (transparency element)" << endl;
                            ExecError("error in reading PLY file (transparency element)");
                        }
                    }
                }
                
                
                
                else {
                    if (indexVertFace && (fscanf(fp, "%d", &numVerts) != 1) ) {
                        
                        cout << "error in reading PLY file (numVerts element)" << endl;
                        ExecError("error in reading PLY file (numVerts element)");
                    }
                    if(verbosity>5)
                        cout << " numVerts test " << numVerts << endl;
                    ffassert(numVerts==B::nv);   /////// pb
                    
                    for (int j = 0; j < numVerts; ++j) {
                        if (fscanf(fp, "%d", &ivb[j]) != 1) {
                            cout << "error in reading PLY file (element)" << endl;
                            ExecError("error in reading PLY file (element)'");
                        }
                        if (verbosity > 9)
                            cout << "  ivb[j]" << ivb[j] << endl;    ////
                    }
                    if (havelabElem && (fscanf(fp, "%d", &label[ibe]) != 1) ){
                        cout << "error in reading PLY file (lab element)" << endl;
                        ExecError("error in reading PLY file (lab element)");
                    }
                    
                    if (havelabElem && verbosity > 9)
                        cout << "  label[it]" << label[ibe] << endl;
                    if (havecolorElem)
                        fscanf(fp, "%d %d %d", &color[0], &color[1], &color[2]);
                    
                    if (haveTransparencyElem) {
                        int transparencyElem;
                        fscanf(fp, "%d", &transparencyElem);
                    }
                }
                (bbff++)->set(vff, ivb, label[ibe]);
            }
    }
        fclose(fp);
        MMesh *pTh = new MMesh(nv, nt, nbe, vff, tff, bff, cleanmesh, removeduplicate, precisvertice);
        return pTh;
    }
    
    
    template <class MMesh>
    AnyType PLY_LoadMeshT_Op< MMesh >::operator( )(Stack stack) const {
        string *pffname = GetAny< string * >((*filename)(stack));
        bool swap(arg(0, stack, false));
        bool cleanmesh(arg(1, stack, false));
        bool removeduplicate(arg(2, stack, false));
        double precisvertice(arg(3, stack, 1e-6));
        
        MMesh *Th = PLY_LoadT<MMesh>(*pffname, swap, cleanmesh, removeduplicate, precisvertice);
        Add2StackOfPtr2FreeRC(stack, Th);
        return Th;
    }
    
    template< class MMesh >
    class PLY_WriteMeshT_Op : public E_F0mps {
    public:
        typedef long Result;
        Expression eTh;
        typedef const MMesh *ppmesh;
        Expression filename;
        static const int n_name_param = 2;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        
        bool arg(int i, Stack stack, bool a) const {
            return nargs[i] ? GetAny< bool >((*nargs[i])(stack)) : a;
        }
        
    public:
        PLY_WriteMeshT_Op(const basicAC_F0 &args) {
            
            if ((std::is_same< MMesh, Mesh3 >::value) && verbosity > 2)
                cout << "Write Mesh3 in PLY Format" << endl;
            else if ((std::is_same< MMesh, MeshS >::value) && verbosity > 2)
                cout << "Write MeshS in PLY Format" << endl;
            else if ((std::is_same< MMesh, MeshL >::value) && verbosity > 2)
                cout << "Write MeshL in PLY Format" << endl;
            args.SetNameParam(n_name_param, name_param, nargs);
            if (BCastTo< string * >(args[0]))
                filename = CastTo< string * >(args[0]);
            if (BCastTo< ppmesh >(args[1]))
                eTh = CastTo< ppmesh >(args[1]);
        }
        
        static ArrayOfaType typeargs( ) {
            return ArrayOfaType(atype< string * >( ), atype< ppmesh >( ), true);
        }
        static E_F0 *f(const basicAC_F0 &args) { return new PLY_WriteMeshT_Op<MMesh>(args); }
        
        AnyType operator( )(Stack stack) const;
    };
    
    template<>
    basicAC_F0::name_and_type PLY_WriteMeshT_Op<Mesh3>::name_param[] = {
        {"floatmesh", &typeid(bool)}, {"bin", &typeid(bool)}};
    template<>
    basicAC_F0::name_and_type PLY_WriteMeshT_Op<MeshS>::name_param[] = {
        {"floatmesh", &typeid(bool)}, {"bin", &typeid(bool)}};
    template<>
    basicAC_F0::name_and_type PLY_WriteMeshT_Op<MeshL>::name_param[] = {
        {"floatmesh", &typeid(bool)}, {"bin", &typeid(bool)}};
    
    
    template< class MMesh >
    void PLY_WRITE_MESHT(const string &filename, FILE *fp, const MMesh &Th, bool binary, int datasize,
                         bool bigEndian) {
        
        typedef typename MMesh::Element T;
        typedef typename MMesh::BorderElement B;
        typedef typename MMesh::Vertex V;
        
        bool swap = bigEndian ? 1 : 0;
        
        fprintf(fp, "ply\n");
        if (binary && !bigEndian)
            fprintf(fp, "format binary_little_endian 1.0\n");
        else if (binary && bigEndian)
            fprintf(fp, "format binary_big_endian 1.0\n");
        else
            fprintf(fp, "format ascii 1.0\n");
        
        fprintf(fp, "comment generated by FreeFEM, %s\n", filename.c_str( ));
        fprintf(fp, "element vertex %d\n", Th.nv);
        if (datasize == sizeof(float)) {
            fprintf(fp, "property float x\n");
            fprintf(fp, "property float y\n");
            fprintf(fp, "property float z\n");
        }
        else if (datasize == sizeof(double)) {
            fprintf(fp, "property double x\n");
            fprintf(fp, "property double y\n");
            fprintf(fp, "property double z\n");
        }
        fprintf(fp, "element face %d\n", Th.nt);
        fprintf(fp, "property list uchar int vertex_indices\n");
        fprintf(fp, "property int flags\n");
        
        
        if(B::nv == 3) {
            fprintf(fp, "element face %d\n", Th.nbe);
            fprintf(fp, "property list uchar int vertex_indices\n");
            fprintf(fp, "property int flags\n");
        }
        
        else if(B::nv == 2) {
            fprintf(fp, "element edge %d\n", Th.nbe);
            fprintf(fp, "property int vertex1\n");
            fprintf(fp, "property int vertex2\n");
        }
        
        fprintf(fp, "end_header\n");
        
        // write mesh vertices
        if (verbosity > 1)
            printf("writing vertex \n");
        if (datasize == sizeof(float)) {
            
            for (unsigned int i = 0; i < Th.nv; i++) {
                const Vertex3 &P = Th.vertices[i];
                float f[3];
                f[0] = P.x;
                f[1] = P.y;
                f[2] = P.z;
                if (binary) {
                    if (swap)
                        FreeFEM::SwapBytes((char *)&f, sizeof(float), 3);
                    fwrite(&f, sizeof(float), 3, fp);
                } else
                    fprintf(fp, "%.8g %.8g %.8g\n", P.x, P.y, P.z);
            }
        }
        else if (datasize == sizeof(double)) {
            for (unsigned int i = 0; i < Th.nv; i++) {
                const Vertex3 &P = Th.vertices[i];
                double f[3];
                f[0] = P.x;
                f[1] = P.y;
                f[2] = P.z;
                if (binary) {
                    if (swap)
                        FreeFEM::SwapBytes((char *)&f, sizeof(double), 3);
                    fwrite((unsigned char *)&f, sizeof(double), 3, fp);
                } else
                    fprintf(fp, "%.15lg %.15lg %.15lg\n", f[0], f[1], f[2]);
            }
        }
        // write  element
        int numVerts=-1;
        if(T::nv == 4 ) numVerts = 4;
        else if(T::nv == 3 ) numVerts = 3;
        else if(T::nv == 2 ) numVerts = 2;
        
        if (verbosity > 1) printf("writing elements \n");
        if (binary) {
            for (int it = 0; it < Th.nt; it++) {
                const T &K(Th.t(it));
                if (swap)
                    FreeFEM::SwapBytes((char *)&numVerts, sizeof(char), 1);
                fwrite(&numVerts, sizeof(char), 1, fp);
                int iv[numVerts + 1];
                for (int ii = 0; ii < numVerts; ii++)
                    iv[ii] = Th.operator( )(K[ii]);
                iv[numVerts]=K.lab;
                if (swap)
                    FreeFEM::SwapBytes((char *)&iv, sizeof(int), numVerts + 1);
                fwrite(&iv, sizeof(int), numVerts + 1, fp);
            }
        }
        else {
            for (int it = 0; it < Th.nt; it++) {
                const T &K(Th.t(it));
                int iv[numVerts + 1];
                iv[0] = numVerts;
                for (int ii = 0; ii < numVerts; ii++)
                    iv[ii + 1] = Th.operator( )(K[ii]);
                if(T::nv == 4 ) fprintf(fp, "%d %d %d %d %d %d\n", iv[0], iv[1], iv[2], iv[3], iv[4], K.lab);
                if(T::nv == 3 ) fprintf(fp, "%d %d %d %d %d\n", iv[0], iv[1], iv[2], iv[3], K.lab);
            }
        }
        
        // border element
        // triangle3
        if(B::nv==3) {
            
            numVerts = 3;
            if (verbosity > 1) printf("writing border elements \n");
            if (binary) {
                for (int it = 0; it < Th.nbe; it++) {
                    const B &K(Th.be(it));
                    if (swap)
                        FreeFEM::SwapBytes((char *)&numVerts, sizeof(char), 1);
                    fwrite(&numVerts, sizeof(char), 1, fp);
                    int iv[numVerts + 1];
                    for (int ii = 0; ii < numVerts; ii++)
                        iv[ii] = Th.operator( )(K[ii]);
                    iv[numVerts]=K.lab;
                    if (swap)
                        FreeFEM::SwapBytes((char *)&iv, sizeof(int), numVerts + 1);
                    fwrite(&iv, sizeof(int), numVerts + 1, fp);
                }
            }
            else {
                for (int it = 0; it < Th.nbe; it++) {
                    const B &K(Th.be(it));
                    int iv[numVerts + 1];
                    iv[0] = numVerts;
                    for (int ii = 0; ii < numVerts; ii++)
                        iv[ii + 1] = Th.operator( )(K[ii]);
                    fprintf(fp, "%d %d %d %d %d\n", iv[0], iv[1], iv[2], iv[3], K.lab);
                }
            }
        }
       // BoundaryEdgeS
        if(B::nv==2) {
            // write edge bd element
            if (verbosity > 1)
                printf("writing edge elements \n");
            numVerts = 2;
            if (binary) {
                for (int ibe = 0; ibe < Th.nbe; ibe++) {
                    const B &K(Th.be(ibe));
                    int iv[numVerts];
                    for (int ii = 0; ii < numVerts; ii++)
                        iv[ii] = Th.operator( )(K[ii]);
                    if (swap)
                        FreeFEM::SwapBytes((char *)&iv, sizeof(int), numVerts + 2);
                    fwrite(&iv, sizeof(int), numVerts, fp);
                }
            }
            else {
                for (int ibe = 0; ibe < Th.nbe; ibe++) {
                    const B &K(Th.be(ibe));
                    int iv[numVerts];
                    for (int ii = 0; ii < numVerts; ii++)
                        iv[ii] = Th.operator( )(K[ii]);
                    fprintf(fp, "%d %d\n", iv[0], iv[1]);
                }
            }
            if (verbosity > 1)
                printf("end writing ply file\n");
            
        }
        
    }
    
    
    template< class MMesh >
    AnyType PLY_WriteMeshT_Op< MMesh >::operator( )(Stack stack) const {
        string *pffname = GetAny< string * >((*filename)(stack));
        MMesh *pTh = GetAny< MMesh * >((*eTh)(stack));
        ffassert(pTh);
        MMesh &Th = *pTh;
        bool bigEndian = isBigEndian( );
        bool swap = bigEndian;
        //string *dataname;
        bool floatmesh(arg(0, stack, false));
        bool binary(arg(1, stack, true));
        int datasize = floatmesh ? sizeof(float) : sizeof(double);
        FILE *fp = fopen((*pffname).c_str( ), "wb");
        if (!fp) {
            cerr << "Unable to open file " << (*pffname).c_str( ) << endl;
            ExecError("error in reading vtk file");
        }
        PLY_WRITE_MESHT<MMesh>(*pffname, fp, Th, binary, datasize, swap);
        return (MMesh *)NULL;
    }
    
#ifndef COMMON_HPDDM_PARALLEL_IO
    static void Load_Init( ) {
        
        
        // if (verbosity)
        if (verbosity && (mpirank == 0))
            cout << " load: ioply " << endl;
        
        Global.Add("saveply", "(", new OneOperatorCode< PLY_WriteMeshT_Op <Mesh3> >);
        Global.Add("saveply", "(", new OneOperatorCode< PLY_WriteMeshT_Op <MeshS> >);
        Global.Add("saveply", "(", new OneOperatorCode< PLY_WriteMeshT_Op <MeshL> >);
        Global.Add("plyload3", "(", new PLY_LoadMeshT <Mesh3>);
        Global.Add("plyloadS", "(", new PLY_LoadMeshT <MeshS>);
        Global.Add("plyloadL", "(", new PLY_LoadMeshT <MeshL>);
            
    }
    
    LOADFUNC(Load_Init)
#endif
