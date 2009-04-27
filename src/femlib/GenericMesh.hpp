// ORIG-DATE:     Dec 2007
// -*- Mode : c++ -*-
//
// SUMMARY  :  Model of generic mesh 1d,2d,3d    
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curi, Paris,  FRANCE 
// AUTHOR   : Frederic Hecht
// E-MAIL   : frederic.hecht@ann.jussieu.fr
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
#ifndef GENERICMESH_HPP_
#define GENERICMESH_HPP_

// la regle de programmation 3  
extern long verbosity;

#include "cassert" 
#include "assertion.hpp" 
#include <cstdlib>
#include <utility>

//#include <algorithm>
//#include <Functional>

#include "RefCounter.hpp"


using namespace ::std;


#include "Serialize.hpp"

#include "GQuadTree.hpp"
// definition R
namespace Fem2D  {
#include "R3.hpp"
#include "Label.hpp"
#include "HashTable.hpp"

  const double UnSetMesure=-1e+200;


  // struct R {}; 
template<int d> struct typeRd {typedef R0 Rd;};
template<> struct typeRd<1> {typedef R1 Rd;};
template<> struct typeRd<2> {typedef R2 Rd;};
template<> struct typeRd<3> {typedef R3 Rd;};

const int NbTypeItemElement = 4;

const int TypeVertex =0;
const int TypeEdge   =1;
const int TypeFace   =2;
const int TypeVolume =3;
    //  add FH  ... april 2009  for peroidic Boundary Condition.
    //  gestion of the permutation 1,2,3
    // here just user order     
    //  NumPerm  : number of the permutation  
    //   p  permutation      n= NumPerm(p)
    //   p1 permutation inv  n1 = NumPerm(p1)
    //   p1[p[i]]=i 
    //   =>  n1   number of the perm /  p[p1[i]]  increase <=> i
    //  SetNumPerm: set the permutation form number 
    
    template<int d> inline int NumPerm(int *p) {ffassert(0);}  
    template<int d> inline int NumPerm1(int *p) {ffassert(0);}  // num perm inverse 
    template<> inline int NumPerm<1>(int *p) { return 0;}
    template<> inline int NumPerm1<1>(int *p) { return 0;}
    template<> inline int NumPerm<2>(int *p) { return p[0] > p[1] ;}
    template<> inline int NumPerm1<2>(int *p) { return p[0] > p[1] ;}
    
    template<> inline int NumPerm1<3>(int *p) { 
	// signe + depart*2
	int k=0,i0=0,i1=1,i2=2,j[3];
	if(p[i0]> p[i1]) swap(i0,i1),k +=1;
	if(p[i1]> p[i2]) swap(i1,i2),k +=1;
	if(p[i0]> p[i1]) swap(i0,i1),k +=1;
	assert(p[i0] < p[i1] && p[i1] < p[i2]);
	// j is inv of i1,i2,i3
	j[i0]=0;j[i1]=1;j[i2]=2;
	return (k%2)+i0*2; // signe + depart*2
    }
    
    template<> inline int NumPerm<3>(int *p) { 
	// signe + depart*2
	int k=0,i0=0,i1=1,i2=2,j[3];
	if(p[i0]> p[i1]) swap(i0,i1),k +=1;
	if(p[i1]> p[i2]) swap(i1,i2),k +=1;
	if(p[i0]> p[i1]) swap(i0,i1),k +=1;
	assert(p[i0] < p[i1] && p[i1] < p[i2]);
	// j is inv of i1,i2,i3
	j[i0]=0;j[i1]=1;j[i2]=2; 
	return (k%2)+ ((j[0]+3)%3)*2; // signe + depart*2
    }    
    // build de permutation     
    template<int d> inline int SetNumPerm(int n,int *p) { ffassert(0);return 0; }// a error}    
    template<int d> inline int SetNumPerm1(int n,int *p) { ffassert(0);return 0; }// a error}    
    
    template<> inline int SetNumPerm<1>(int n,int *p) { p[0]=0;} // a error}    
    template<> inline int SetNumPerm<2>(int n,int *p) { p[0]=n;p[1]=1-n;} // a error}    
    
    // build perm inverse
    template<> inline int SetNumPerm1<1>(int n,int *p) { p[0]=0;} // a error}    
    template<> inline int SetNumPerm1<2>(int n,int *p) { p[0]=n;p[1]=1-n;} // a error} 
    
    template<> inline  int SetNumPerm1<3>(int n,int *p) { 
	int i=n/2, j= n%2 ? 2:1;
	p[i]=0;p[(i+j)%3]=1;p[(i+j+j)%3]=2;
	assert( n == NumPerm1<3>(p));
    }    
    
    template<> inline   int SetNumPerm<3>(int n,int *p) { 
	int i=n/2, j= n%2 ? 2:1;
	p[0]=i;p[1]=(i+j)%3;p[2]=(i+j+j)%3;
	assert( n == NumPerm<3>(p));
    } 
    // ---  end add periodic 
    
class DataFENodeDF {
    int  * nbref; // pointer  on common ref counter 
public:
  int ndfon[4];
  const int NbOfElements;
  const  int NbOfNodes;
  const  int NbOfDF;  
  const int * const NodesOfElement;
  const int * const FirstDfOfNodeData;
  const int * const FirstNodeOfElement; //  0 
  const int MaxNbNodePerElement;
  const int MaxNbDFPerElement;

  int ndfonVertex()const {return ndfon[0];}
  int ndfonEdge()const {return ndfon[1];}
  int ndfonFace()const {return ndfon[2];}
  int ndfonTet()const {return ndfon[3];}

  DataFENodeDF(const DataFENodeDF & m)
    :
    nbref( m.nbref ) ,
    NbOfElements(m.NbOfElements),
    NbOfNodes(m.NbOfNodes),
    NbOfDF(m.NbOfDF),
    NodesOfElement(m.NodesOfElement),
    FirstDfOfNodeData(m.FirstDfOfNodeData),
    FirstNodeOfElement(m.FirstNodeOfElement),
    MaxNbNodePerElement(m.MaxNbNodePerElement),
    MaxNbDFPerElement(m.MaxNbDFPerElement) 
    {
      for(int i=0;i<NbTypeItemElement;++i)
	    ndfon[i]=m.ndfon[i];
	(*nbref)++; // add one to the ref counter 
    }
  DataFENodeDF( 
	       int andfon[NbTypeItemElement],
	       int aNbOfElements,
	       int aNbOfNodes,
	       int aNbOfDF,
	       const int * aNodesOfElement,
	       const int * aFirstDfOfNodeData,
	       int aMaxNbNodePerElement,
	       int aMaxNbDFPerElement)
    :
    nbref( new int(0)),// new ref counter 
    NbOfElements(aNbOfElements),
    NbOfNodes(aNbOfNodes),
    NbOfDF(aNbOfDF),
    NodesOfElement(aNodesOfElement),
    FirstDfOfNodeData(aFirstDfOfNodeData),
    FirstNodeOfElement(0),
    MaxNbNodePerElement(aMaxNbNodePerElement),
    MaxNbDFPerElement(aMaxNbDFPerElement) 
  {
    for(int i=0;i<NbTypeItemElement;++i)
      ndfon[i]=andfon[i];
  }
 ~DataFENodeDF()  
  {
    if ((*nbref)==0) // remove if nbref ==0
       {
	 delete nbref;
         delete [] NodesOfElement;
         delete [] FirstDfOfNodeData;
         delete [] FirstNodeOfElement;
       }
	else  (*nbref)--;
  }
private:
	void operator=(const DataFENodeDF &) ;    
	
};

template<typename Rn>
class GenericVertex : public Rn,public Label
{
    
  template<typename T,typename B,typename V>
        friend class GenericMesh;
  friend inline ostream& operator <<(ostream& f, const GenericVertex & v )
  { f << (const Rn &) v << ' ' << (const Label &) v   ; return f; }
  friend inline istream& operator >> (istream& f,  GenericVertex & v )
        { f >> (Rn &) v >> (Label &) v ; return f; }
    
  Rn *normal; // pointeur sur la normal exterieur pour filtre des points de departs
   
public:
  typedef Rn Rd;
  static const int d=Rd::d;  
  GenericVertex() : Rd(),Label(),normal(0) {};
  GenericVertex(const Rd & P,int r=0): Rd(P),Label(r),normal(0){}
    
  void SetNormal(Rd *&n,const Rd & N)
  { if (normal) { 
      Rd NN=*normal+N; 
      *normal= NN/NN.norme(); }
    else *(normal=n++)=N;}
    
  Rd Ne() const {return normal ? *normal: Rd();}
  bool ninside(const Rd & P) const
  {
    return normal? (Rd(*this,P),*normal)<=0: true;
  }
    
private: // pas de copie pour ne pas prendre l'adresse
  GenericVertex(const GenericVertex &);
  void operator=(const GenericVertex &); 
  
};
  inline  R1 ExtNormal( GenericVertex<R1> *const v[2],int const f[1])  {   return f[0]==0 ? R1(-1):R1(1);  }
  inline  R2 ExtNormal( GenericVertex<R2> *const v[3],int const f[2])  {   return R2(*v[f[1]],*v[f[0]]).perp();  }
  inline  R3 ExtNormal( GenericVertex<R3> *const v[4],int const f[3])  {   return R3(*v[f[1]],*v[f[0]])^R3(*v[f[2]],*v[f[0]]) ;  }
    

template<typename Data>  
class GenericElement: public Label {
public:
  typedef typename Data::V Vertex;
  typedef typename Data::V::Rd Rd;
  typedef typename Data::RdHat RdHat;// for parametrization 
  typedef typename Data::RdHatBord RdHatBord;// for parametrization 

  typedef typename Rd::R R;

  static const int nv=Data::NbOfVertices;  // nb  of vertices 
  static const int ne=Data::NbOfEdges;  // nb  of edges
  static const int nf=Data::NbOfFaces;  // nb of faces
  static const int nt=Data::NT;  // nb of tets
  static const int nitem=nv+ne+nf+nt;
  static const int nva=Data::NbOfVertexOnHyperFace;  
  static const int nea=Data::NbOfAdjElem;  
  static const int d=Rd::d;  
  static const int (* const nvedge)[2] ;//
  static const int (* const nvface)[3] ;//
  static const int (* const onWhatBorder)[nitem] ;//
    
  static const int (* const nvadj)[nva] ;//  
  static const int nitemdim[4]; //  nv,ne,nf,nt 
  // variable prive 
private:
  Vertex *vertices[nv]; // an array of 3 pointer to vertex
  R mes; 
public:
  GenericElement() {}
  const Vertex & operator[](int i) const  {
    ASSERTION(i>=0 && i <nv);
    return *vertices[i];} // to see triangle as a array of vertex

  Vertex & operator[](int i)  {
    ASSERTION(i>=0 && i <nv);
    return *vertices[i];} // to see triangle as a array of vertex

  const Vertex& at(int i) const   { return *vertices[i];} // to see triangle as a array of vert

  Vertex& at(int i)    { return *vertices[i];} // to see triangle as a array of vert

  void set(Vertex * v0,int * iv,int r,double mss=UnSetMesure) 
  { 
    for(int i=0;i<nv;++i)	
      vertices[i]=v0+iv[i];
    mes=(mss!=UnSetMesure) ? mss : Data::mesure(vertices);
    lab=r;
    ASSERTION(mss==UnSetMesure && mes>0);
  }

  
  istream & Read1(istream & f,Vertex * v0,int n)
  {
    int iv[nv],ir,err=0;
    for (int i=0;i<nv;++i) 
      {
	f >> iv[i];
	iv[i]--;
	if (  ! (iv[i]>=0 && iv[i]<n)) err++;
      }
    f >> ir;
    if(err || ! f.good() )
      {
	cerr << " Erreur GenericElement::Read1 " << nv <<  " " << n << "  : " ;
	for (int j=0;j<nv;++j) 
	  cerr << iv[j] <<  " ";
	cerr << " , " << ir <<  endl;
	abort();
      }

    set(v0,iv,ir); 
    return f;
  }
     	
  Rd Edge(int i) const {ASSERTION(i>=0 && i <ne);
    return Rd(at(nvedge[i][0]),at(nvedge[i][1]));}// opposite edge vertex i

  Rd N(int i) const  { return ExtNormal(vertices,nvadj[i]);}
  Rd PBord(int i,RdHatBord P) const   { return Data::PBord(nvadj[i],P);}  

  Rd operator()(const RdHat & Phat) const {
    Rd r= (1.-Phat.sum())*(*(Rd*) vertices[0]);    
    for (int i=1;i<nv;++i)
      r+=  Phat[i-1]*(*(Rd*) vertices[i]);
    return r;
  }


  int faceOrientation(int i) const 
  {// def the permutatution of orient the face
    int fo =0;
    Vertex * f[3]={&at(nvface[i][0]), &at(nvface[i][1]), &at(nvface[i][2])}; 
    if(f[0]>f[1]) fo+=1,Exchange(f[0],f[1]); 
    if(f[1]>f[2]) { fo+=2,Exchange(f[1],f[2]); 
    if(f[0]>f[1]) fo+=4,Exchange(f[0],f[1]); }
    return fo;
  }

  bool   EdgeOrientation(int i) const 
    { return &at(nvedge[i][0]) < &at(nvedge[i][1]);}
    
  R lenEdge(int i) const {ASSERTION(i>=0 && i <3);
    Rd E=Edge(i);return sqrt((E,E));}

  R  mesure() const {return mes;}


    
  static  int NbNodes(int c)  // use the bit i of c to say if node in objet of dim  i existe
  { int c0=(c&1)!=0, c1=(c&2)!=0, c2=(c&4)!=0, c3=(c&8)!=0;  
    return nv*c0 +ne*c1 +nf*c2 + nt*c3 ;}

    static  int NbNodes(const int c[4])  // use the bit i of c to say if node in objet of dim  i existe
  { int c0=(c[0])!=0, c1=(c[1])!=0, c2=(c[2])!=0, c3=(c[3])!=0;  
    return nv*c0 +ne*c1 +nf*c2 + nt*c3 ;}

  void Renum(Vertex   *v0, int * r)  
  { 
    for (int i=0;i<nv;i++) 
      vertices[i]=v0+r[vertices[i]-v0];
  }

  void Change(Vertex   *vold, Vertex *vnew) 
  { 
    for (int i=0;i<nv;i++) 
      vertices[i]=vnew+vertices[i]-vold;
  }

  //Rd n(int i) const //  unit exterior normal
  //  {Rd E=Edge(i);return Rd(E.y,-E.x)/Norme2(E);} 
   
    
private:
  // pas de copie
  GenericElement(const  GenericElement &);
  GenericElement &operator = (const  GenericElement &);
};

 
    template<int N> inline void PermI2J(const void **I,const void **J,int *S)
    {
	ffassert(0);
    }
    template<> inline void PermI2J<1>(const void **I,const void **J,int *S)
    {
	S[0]=0;
    }
    template<> inline void PermI2J<2>(const void **I,const void **J,int *S)
    {
	if(I[0]==J[0])
	  { assert(I[1]==J[1]);
	    S[0]=0;S[1]=1;}
	else 
	  { assert(I[1]==J[0]&&I[0]==J[1]);
	    S[0]=1;S[1]=0;}	
    }
    template<> inline void PermI2J<3>(const void **I,const void **J,int *S)
    {
	if(I[0]==J[0]) S[0]=0;
	else if(I[0]==J[1]) S[0]=1;
	else {S[0]=2; assert(I[0]==J[2]) ;}
	if(I[1]==J[0]) S[1]=0;
	else if(I[1]==J[1]) S[1]=1;
	else {S[1]=2; assert(I[1]==J[2]) ; }
	S[2]=3-S[0]-S[1];
	assert(I[2]==J[3-S[0]-S[1]]);
    }
    
template<typename T,typename B,typename V>
class GenericMesh : public RefCounter
{ 
public:
  typedef GenericMesh GMesh;
  typedef T Element;   
  typedef typename V::Rd Rd;
  typedef typename Rd::R R;
  typedef V  Vertex;
  typedef B BorderElement;
  typedef  EF23::GTree<V> GTree;
  typedef typename Element::RdHat RdHat;// for parametrization                                        

  int nt,nv,nbe;
  R mes,mesb;
  //private:
  V *vertices;
  T *elements;
  B *borderelements;
  Rd  * bnormalv; //  boundary vertex normal 
  Rd Pmin,Pmax; // // the bound  of the domain  see BuildBound 
  static const int nea=T::nea; //  numbering of adj (4 in Tet,  3 in Tria, 2 in seg) 
  static const int nva=T::nva; //  numbering of vertex in Adj element 
  static int kfind,kthrough; //  number of search and number of throught element. 
  int *TheAdjacencesLink; // to store the adj link  k*nea+i -> k'*nea+i' 
  int *BoundaryElementHeadLink; //
  int *ElementConteningVertex;  
  GTree *gtree;
public:
  const T & operator[](int i) const {return elements[CheckT(i)];}
  const V& operator()(int i) const {return vertices[CheckV(i)];}
  const B& be(int i) const {return borderelements[CheckBE(i)];}
  
  T & t(int i)  {return elements[CheckT(i)];}
  V & v(int i)  {return vertices[CheckV(i)];}
  B & be(int i) {return borderelements[CheckBE(i)];}
  
  GenericMesh()
    : nt(0),nv(0),nbe(0),  mes(0.),mesb(0.) ,
      vertices(0),elements(0),borderelements(0),bnormalv(0),
      TheAdjacencesLink(0),BoundaryElementHeadLink(0),
      ElementConteningVertex(0), gtree(0)
  {} 
  
  void set(int mv,int mt,int mbe) 
  {
    assert(nt==0 && nv==0 && nbe ==0);
    nt=mt;
    nv=mv;
    nbe=mbe;
    vertices=new V[nv];
    elements= new T[nt];
    borderelements = new B[nbe]; 
    assert( nt >=0 && elements);
    assert( nv >0 && vertices);
    
  }
 
 
  int operator()(const T & t) const {return CheckT(&t - elements);}
  int operator()(const T * t) const {return CheckT(t - elements);}
  int operator()(const V & v) const {return CheckV(&v - vertices);}
  int operator()(const V  * v) const{return CheckV(v - vertices);}
  int operator()(const B & k) const {return CheckBE(&k - borderelements);}
  int operator()(const B  * k) const{return CheckBE(k - borderelements);}
  int operator()(int it,int j) const {return operator()(elements[it][j]);}// Nu vertex j of triangle it
  int be(int it,int j) const {return operator()(borderelements[it][j]);}// Nu vertex j of triangle it
  
  int CheckV(int i) const { ASSERTION(i>=0 && i < nv); return i;} 
  int CheckT(int i) const { ASSERTION(i>=0 && i < nt); return i;}
  int CheckBE(int i) const { ASSERTION(i>=0 && i < nbe); return i;}
  
   
  int Contening(const Vertex * v) const{ return ElementConteningVertex[ v  - vertices];} 
  void BuildAdj();
  void BuildSurfaceAdj();  // Add J. Morice function that give the TheAdjacencesSurfaceLink
  void Buildbnormalv();
  void BuildBound();
  void BuildjElementConteningVertex();
  void BuildGTree() {if(gtree==0)  gtree=new GTree(vertices,Pmin,Pmax,nv);}    
  DataFENodeDF BuildDFNumbering(int dfon[NbTypeItemElement],int nbequibe=0,int *equibe=0) const ;
    DataFENodeDF BuildDFNumbering(int ndfv,int ndfe,int ndff,int ndft,int nbequibe=0,int *equibe=0) const 
  { int dfon[NbTypeItemElement]={ndfv,ndfe,ndff,ndft};
    return  BuildDFNumbering(dfon,nbequibe,equibe);
  }
  
  int ElementAdj(int k,int &j) const  {
    int p=TheAdjacencesLink[nea*k+j];
    j=p%nea;
    return p>=0 ? p/nea: -1;}

  int ElementAdj(int k,int &j,Rd& PHat) const  
    {
    //   return the kk the number of adj element k  to hyperface j (opposite to vertex j)
    //   out j: is the new hyperface number in element kk.
    // and 
    // in : Pt is the point  on hyperface j on element k on ref element K hat.
    //  remark   lb[j]==0  at enter
    // you get the new 	point Pt (in  on hyperface j on element kk
    //  and lb[j] ==0 at return (j have change).
    int p=TheAdjacencesLink[nea*k+j];
    if(p>=0) 
      {
	
	  R lb[Rd::d+1];//{1.-PHat.sum(),PHat}; 
	  R lbb[Rd::d+1];//{1.-PHat.sum(),PHat}; 
	  PHat.toBary(lb); // R1 R2 R3 
	  if(Abs(lb[j])>1e-10)
	   assert(Abs(lb[j])<1e-10);
	int sigma[T::nva];	  
	const void * nvkj[T::nva], *nvkkjj[T::nva];
	int jj=p%nea;
	int kk=p/nea;

	Element & K(elements[CheckT(k)]);
	Element & KK(elements[CheckT(kk)]);
	Rd Pin=K(PHat);  
	for (int l=0;l<T::nva;++l)
	    nvkj[l] =&K[T::nvadj[j][l]];
	for (int l=0;l<T::nva;++l)
	    nvkkjj[l] = &KK[T::nvadj[jj][l]];
	//  il faut permute ll.
	PermI2J<nva>(nvkj,nvkkjj,sigma);
	for (int l=0;l<T::nva;++l)
	    lbb[T::nvadj[jj][l]]=lb[T::nvadj[j][sigma[l]]];
	lbb[jj]=0;
#ifdef DEBUG	  
	Rd PH=PHat;  
#endif
	PHat=Rd(lbb+1);
#ifdef DEBUG	  
	Rd Pout=KK(PHat);
	if( (Pin-Pout).norme2() > 1e-10 )
	    {
		for (int l=0;l<=T::nva;++l)
		    cout << lbb[l] <<" < -- " << lb[l] << endl;
		for (int l=0;l<T::nva;++l)
		    cout <<l << " :    o=  " << nvkkjj[l]  << "   i= " << nvkj[l] << " " <<  sigma[l] 
		         << " -- " << &KK[T::nvadj[jj][l]]  << " == " << &K[T::nvadj[j][sigma[l]]] 
		    << " -- " << &K[T::nvadj[j][l]]  << " == " << &KK[T::nvadj[jj][sigma[l]]] 
		    << " -- " << lbb[T::nvadj[jj][l]] << " == " << lb[T::nvadj[j][sigma[l]]]
		    << " ++ " << T::nvadj[jj][l] << " <-- " << T::nvadj[j][sigma[l]] 
		    << endl;
		cout << "Adj:  j= " << j << " ," << Pin << " != " << Pout << " , " << PH << " -> " << PHat << "  jj = " << jj <<  endl;
		assert(0);
	    }
#endif	  
	j=jj;
	return kk;
      }
    return -1;//  on border 
  }
    
  int GetAllElementAdj(int it,int *tabk) const
  { //  get the tab of all adj element (max ne)
    //  and return the size of the tab 
    int i=0;
    for(int j=0;j<nea;++j)
      {
	tabk[i]=TheAdjacencesLink[3*it+j]/3;
	if(tabk[i] >=0 && tabk[i]!=it) i++; 
      }
    return i;
  }

  int BoundaryElement(int be,int & ItemInK) const {
    int i= BoundaryElementHeadLink[be]; 
    ItemInK = i%nea; 
    return i/nea;}
  
  // Add J. Morice 
  template<int N,int M>
  SortArray<int,N> itemadjs(const int (* const  nu )[N],int k,int i, int *sens) 
  {
    int nv[N];
    B & K(borderelements[CheckBE(k)]);
    ASSERTION(i>=0 && i <M);
    for (int j=0;j<N;++j){
      nv[j] = operator()(K[nu[i][j]]);
    }
    if(nv[0] > nv[1] )
      *sens = 1;
    else
      *sens =-1;
    return SortArray<int,N>(nv);
  }

  SortArray<int,B::nva> items(int k,int i,int *sens) 
  {
    return itemadjs<B::nva,B::nv>(B::nvadj,k,i,sens);
  }

  
  template<int N,int M>
  SortArray<int,N> iteme(const int (* const  nu )[N],int k,int i) 
  {
    int nv[N];
    Element & K(elements[CheckT(k)]);
    ASSERTION(i>=0 && i <M);
    for (int j=0;j<N;++j){
      nv[j] = operator()(K[nu[i][j]]);
    }

    return SortArray<int,N>(nv);
  }

  SortArray<int,B::nv> itemadj(int k,int i) 
  {
    return iteme<B::nv,T::nea>(T::nvadj,k,i);
  }
  
  SortArray<int,B::nv> itembe(int k) 
  {
    int nv[B::nv];
    B & K(borderelements[CheckBE(k)]);
    
    for (int j=0;j<B::nv;++j){
      nv[j] = operator()(K[j]);
    }

    return SortArray<int,B::nv>(nv);
  }

  //  const Element * Find(const Rd & P) const ;
  const Element * Find(Rd P, RdHat & Phat,bool & outside,const Element * tstart=0) const  
  {return EF23::Find<GMesh>(*this,this->gtree,P,Phat,outside,tstart);}
  
  R mesure(){ return mes;}
  R bordermesure(){ return mesb;}
  virtual ~GenericMesh() { 
    cout << "~GenericMesh\n";
   
    delete [] ElementConteningVertex;
    delete [] TheAdjacencesLink;
    delete [] BoundaryElementHeadLink;
    delete [] borderelements;
    if(nt>0) delete [] elements;
    delete [] vertices;
    delete [] bnormalv;
    if(gtree) delete gtree;
    
  }

  Serialize serialize() const;
  
private:
  GenericMesh(const GenericMesh &); // pas de construction par copie
   void operator=(const GenericMesh &);// pas affectation par copy 
};

template<typename T,typename B,typename V>
void GenericMesh<T,B,V>::BuildjElementConteningVertex()
{
  const int nkv= T::nv;
  if(!ElementConteningVertex) ElementConteningVertex = new int[nv];

    for(int i=0;i<nv;++i)
	ElementConteningVertex[i]=-1; 
    
    for (int k=0;k<nt;++k)
	for (int i=0;i<nkv;++i)
	    ElementConteningVertex[operator()(elements[k][i])]=k ;
    int kerr=0;
    for(int i=0;i<nv;++i)
	if (ElementConteningVertex[i]<0) 
	    kerr++; 
    assert(kerr==0);

} 
template<typename T,typename B,typename V>
void GenericMesh<T,B,V>::BuildAdj()
{
  const int nva   = T::nva;
  const int nea   = T::nea;
  assert(TheAdjacencesLink==0);
  TheAdjacencesLink = new int[nea*nt];
  BoundaryElementHeadLink = new int[nbe];
  HashTable<SortArray<int,nva>,int> h(nea*nt,nv);
  int nk=0,nba=0;
  int err=0;

  cout << "nva=// nea=" << nva << " " << nea << " "<< nbe << endl;
  for (int k=0;k<nt;++k)
    for (int i=0;i<nea;++i)
      {
        SortArray<int,nva> a(itemadj(k,i));
	//cout << " ### "   << " item(k,i)= " << itemadj(k,i) << " a= " << a << " k " << k << " i " << i << endl;
	typename HashTable<SortArray<int,nva>,int>::iterator p= h.find(a);
	if(!p) 
	  { 
	    h.add(a,nk);
	    TheAdjacencesLink[nk]=-1;
	    nba++;
	  }
	else 
	  {	  
	    ASSERTION(p->v>=0);
	    TheAdjacencesLink[nk]=p->v;
	    TheAdjacencesLink[p->v]=nk;
	    p->v=-1-nk;
	    nba--;
	  }
	++nk;
      }
    
  for (int k=0;k<nbe;++k)
     {
	SortArray<int,nva> a(itembe(k));

	typename HashTable<SortArray<int,nva>,int>::iterator p= h.find(a);
	//cout << k << " ### "   << " item(k,i)= " << itembe(k) << " a= " << a << endl;
	if(!p) { err++;
	if(err==1) cerr << "Err  Border element not in mesh \n";
	if (err<10)  cerr << " \t " << k << " " << a << endl;
	}
	 else
	   {
	     BoundaryElementHeadLink[k] = p->v <0 ? -p->v-1 : p->v;
	     #ifndef NDEBUG
	     int t=BoundaryElementHeadLink[k]/nea;
	     int e=BoundaryElementHeadLink[k]%nea;
	     //cout << k << " ### "   << a << " = " << itemadj(t,e) << " t " << t << " e " << e << endl;
	     assert(itemadj(t,e)==a);
	     #endif
	   }
     }

    
  assert(err==0);
  int na= h.n;
  cout << "  -- Nb adj  = "<< na << " on border " << nba << " nea = " << nea << " nva = " << nva ;
  if(nea==2)
    cout << " Const d'Euler: " << nt - na + nv << endl;
  else
    cout << endl;
}

template<typename T,typename B,typename V>
void GenericMesh<T,B,V>::BuildSurfaceAdj()
{
 
  // assert(TheSurfaceAdjacencesLink==0); plus tard
  int *TheSurfaceAdjacencesLink = new int[B::nea*nbe];
  HashTable<SortArray<int,B::nva>,int> h(B::nea*nbe,nv);
  int nk=0;
  int err=0;
  int sens;

  cout << "nea/nva" << B::nea << " "  << B::nva << endl;
  for (int k=0;k<nbe;++k)
    for (int i=0;i<B::nea;++i)
      {
        SortArray<int,B::nva> a(items(k,i,&sens));

	typename HashTable<SortArray<int,B::nva>,int>::iterator p= h.find(a);
	if(!p) 
	  { 
	    h.add(a,nk);
	    TheSurfaceAdjacencesLink[nk]=sens;
	  }
	else 
	  {	  
	    ASSERTION(p->v>=0);
	    if( sens*TheSurfaceAdjacencesLink[p->v] != -1){
	 
	      B & K(borderelements[CheckBE(k)]);
	      int firstVertex  =  operator()(K[B::nvadj[i][0]])+1;
	      int secondVertex =  operator()(K[B::nvadj[i][1]])+1;
	      cout << " The edges defined by vertex is " << firstVertex << "-" << secondVertex << " is oriented twicd in the same direction "<< endl;
	      err++;
	    }
	    TheSurfaceAdjacencesLink[nk]=p->v;
	    TheSurfaceAdjacencesLink[p->v]=nk;    
	  }
	if( err > 10 ) 
	  exit(1); 
	nk++;
      }
    
  assert(err==0);
  delete [ ] TheSurfaceAdjacencesLink; 
  cout << "number of adjacents edges " << nk << endl; 
}



template<typename T,typename B,typename V>
DataFENodeDF GenericMesh<T,B,V>::BuildDFNumbering(int ndfon[NbTypeItemElement],int nbequibe,int *equibe) const
{
/*
   nbequibe nb of  borderelement with equi boundary condition
  for i =0, 2*nbequibe, i+= 2)
    
     be0= equibe[i]/8 <=>  be1=equibe[i+1] /8
     equibe[i]%8    given the permuation  p0 compare to sort array. 
     equibe[i+1]%8 given the permuation   p1 compare to sort array.
    the  numbering of perumation 
     SetNumPerm<nkf>(p0,equibe[i+1]%8);
     SetNumPerm<nkf>(p1,equibe[i+1]%8);
    
    so a level of point with always have:
      be0[p0[j]] <==>  be1[p1[j]]
 
 
 */
  const GenericMesh & Th(*this);
  int nnodeK = T::NbNodes(ndfon);
  int *p = 0, *pp=0; 
  unsigned int tinfty=-1;
    
  const int nkv= T::nv;
  const int nkf= T::nf;
  const int nke= T::ne;
  const int nkt= T::nt;
  const int nbev= B::nv;
  const int nbef= B::nf;
  const int nbee= B::ne;
  const int nk[]={nkv,nke,nkf,nkt};
  int MaxNbNodePerElement=0;
  int MaxNbDFPerElement=0;
  int nbNodes=0;
  int NbOfDF=0;
  int n=0;
  int minndf=100000000;
  int maxndf=0;
  int nbnzero=0;
  for (int dd=0;dd<NbTypeItemElement;++dd)
    if(ndfon[dd])
      {
	
        nbnzero++;
	minndf=Min(minndf,ndfon[dd]);
	maxndf=Max(maxndf,ndfon[dd]);
	MaxNbDFPerElement   += ndfon[dd]*nk[dd];
	MaxNbNodePerElement += nk[dd];
      }
  bool constndfpernode = minndf == maxndf;
  bool nodearevertices = ( nbnzero ==1  && ndfon[0]) && nbequibe==0 ;
 
  assert(maxndf>0);
  const int nkeys=4+6+4+1;
  assert(nnodeK<= nkeys);
  
  if(nodearevertices)
    {
      nbNodes=nv;
      NbOfDF=nbNodes*ndfon[0];
    }
  else
    { 
      p =  new int[nnodeK*nt];
      typedef SortArray<unsigned int,2> Key;
	
      Key keys[nkeys*2];
      int keysdim[nkeys*2];
      
      int of = Th.nv+1;
      int ndim[NbTypeItemElement]={0,0,0,0};
      NbOfDF=0;
      {
	HashTable<Key,int> h(nnodeK*nt,of+nt);
	int  nbmaxeq = 1+nnodeK*nbequibe;
	int  nbhasheq = nbequibe ? of+nt : 1;
	HashTable<Key,Key> equi(nbmaxeq,nbhasheq); 
	  //  constuction of item translation for 
	  if(verbosity>2)
	  cout << " nb equi be :  " << nbequibe << endl;
	  for(int ieq=0,keq=0;keq<nbequibe;++keq)
	    {
		int p1[nbev],p2[nbev];
		int v1[nbev],v2[nbev];
		int be1=equibe[ieq]/8,pe1=equibe[ieq++]%8;
		int be2=equibe[ieq]/8,pe2=equibe[ieq++]%8;
		int itemb1,k1=BoundaryElement(be1,itemb1);
		int itemb2,k2=BoundaryElement(be2,itemb2);
		const B& b1(Th.be(be1));
	        const B& b2(Th.be(be2));
		SetNumPerm<nbev>(pe1,p1);
		SetNumPerm<nbev>(pe2,p2);
		int m=0;
		for(int i=0;i<nbev;++i)
		  {
		      v1[i]=Th(b1[p1[i]]);
		      v2[i]=Th(b2[p2[i]]);
		  }
		
		int dimb= B::RdHat::d;
		if(ndfon[dimb])
		  { //  simple border element
		  }
                if( ndfon[0] )
		    for(int i=0;i<nbev;++i)
		  {
		    keysdim[m]=0,keys[m++]=Key(v1[i],tinfty);
		    keysdim[m]=0,keys[m++]=Key(v2[i],tinfty);
		      if(verbosity >100) cout << be1<< " " << be2 << "  " <<  v1[i] << " <--> " << v2[i] 
			                 << " /  " << Th(v1[i]) << " <=> " << Th(v2[i]) << endl;
		  }
		if( ndfon[1] )
		  for(int i=0;i<nbee;++i)
		  {
		      keysdim[m]=1,keys[m++]=Key(v1[B::nvedge[i][0]],v1[B::nvedge[i][1]]);
		      keysdim[m]=1,keys[m++]=Key(v2[B::nvedge[i][0]],v2[B::nvedge[i][1]]);
		  }
		if( ndfon[2] && dimb ==2)
		  {
		      assert(nbef==1 && nkf != 1);
		      int ii;
		      keysdim[m]=2,keys[m++]=Key(k1+of,ElementAdj(k1,ii=itemb1) +of);
		      keysdim[m]=2,keys[m++]=Key(k2+of,ElementAdj(k2,ii=itemb2) +of);		     
		  }
		        
		for(int j=0;j<m;)
		  {
		      int i0=j++,i1=j++; 
		      if(keys[i1]<keys[i0]) swap( keys[i0],keys[i1]);
		      typename HashTable<Key,Key>::iterator pe = equi.add(keys[i0],keys[i1]);
		     // if(pe) assert(pe->k == keys[i0]);
		  }
		
	    }
	  
	  //  to find the final equivalent  key ... 
	  //  in change of chaibe of equi key
	  for (int it=0,change=1;change;it++)
	    { 
	    change=0;
	    assert(it<10);
	    for (typename HashTable<Key,Key>::iterator qe,pe=equi.begin() ; pe != equi.end(); ++pe)
	      { 
		  
		  assert( pe->k < pe->v); 
		  qe=equi.find(pe->v);
		  if(qe) 
		    { change++;
		     assert( qe->k < qe->v);
		      pe->v = qe->v;
		  }
	    
		  
	      }
	      cout << " iteration in final equivalent " << it << " nb change " << change << endl;
	    }
	    
	
	  //  construction of nodes numbering 
	  
	  // ------------
	  
	for(int k=0;k<nt;++k)
	  {
	    const T & K(Th[k]);
	    int m=0;
	    if( ndfon[0] )//  node on vertex
	      for(int i=0;i<nkv;++i)
		keysdim[m]=0,keys[m++]=Key(Th(K[i]),tinfty);
	    if(  ndfon[1] )//  node on Edge
	      for(int i=0;i<nke;++i)
		keysdim[m]=1,keys[m++]=Key(Th(K[T::nvedge[i][0]]),Th(K[T::nvedge[i][1]]));
	    if(  ndfon[2])//  node on Face
	      if (nkf==1) 
		keysdim[m]=2,keys[m++]=Key(k+of,tinfty);
	      else
		for(int ii,i=0;i<nkf;++i)
		  keysdim[m]=2,keys[m++]=Key(k+of,ElementAdj(k,ii=i) +of);
	    if(  ndfon[3] )//  node on Tet
	      if(nkt==1)
		keysdim[m]=3,keys[m++]=Key(k+of,tinfty);
	    
	    if(k<0)
	    {
	      for(int i=0;i<nke;++i)
		cout << " e= " << T::nvedge[i][0] << " " << T::nvedge[i][1] << endl;
	      cout << ndfon[0] << " " << ndfon[1] << " " << ndfon[2] << " " << ndfon[3] << ": "
		   <<  " m = "<< m << "  " <<nnodeK
		   << " " << T::nv
		   << " " << T::ne 
		   << " " << T::nf 
		   << " " << T::nt  
		   << endl;
	    }
	    assert(m==nnodeK);
	    for(int i=0;i<m;i++)
	      {
		Key ki=keys[i],kio=ki; 
		typename HashTable<Key,Key>::iterator pe= equi.find(ki);
		  if(pe) {   ki= pe->v;  }
		typename HashTable<Key,int>::iterator pk= h.find(ki);
		if(!pk)
		 {
		  pk = h.add(ki,nbNodes++);
		  NbOfDF += ndfon[keysdim[i]];
		 }
		if(verbosity>100 && pe ) cout << kio << " -> " << pe->k  << " :: " <<  pk->v << endl;   
		p[n++] = pk->v;
		ndim[keysdim[i]]++;
			  
	      }
	  }
      }
      cout << " nb of Nodes " << nbNodes << endl;
      cout << " nb of DoF   " << NbOfDF << endl ;
      if( ! constndfpernode) 
	{ 
	  pp=new int[nbNodes+1];
	  int kk=0,nn=0;
	  for(int k=0;k<nt;++k)
	    for(int i=0;i<nnodeK;i++)
		pp[p[nn++]]=ndfon[keysdim[i]];;
	  for(int n=0;n<nbNodes;++n)	
	    {  
	      int ndfn=pp[n];
	      pp[n]=kk;
	      kk+=ndfn;
	    }
	  pp[nbNodes]=NbOfDF;//  add last 
	  assert(kk==NbOfDF);
	}
    }
  return DataFENodeDF(ndfon,nt,nbNodes,NbOfDF,p,pp,MaxNbNodePerElement,MaxNbDFPerElement);
}
template<typename T,typename B,typename V>
void GenericMesh<T,B,V>::BuildBound()
{
    if(vertices && (nv>0))
    {
	Pmin=vertices[0];
	Pmax=vertices[0];
	 for(int i=1;i<nv;++i)
	 {
	     Pmin=Minc(Pmin,vertices[i]);
	     Pmax=Maxc(Pmax,vertices[i]);

	 }
    }
	
}

template<typename T,typename B,typename V>
void GenericMesh<T,B,V>::Buildbnormalv()
{
    const int nkv= T::nv;
   // const int nkf= T::nf;
   // const int nke= T::ne;
   // const int nkt= T::nt;
    
    if (bnormalv) 
      {assert(0);return;}
    int nb=0;
    for (int k=0;k<nt;k++)
	for (int i=0;i<nkv;i++)
	{  
	    int ii(i),kk;
	    kk=ElementAdj(k,ii);
	    if (kk<0 || kk==k) nb++;
	}
    if(verbosity>2)
	cout << " number of real boundary  " << nb << endl;
    bnormalv= new Rd[nb];
    Rd *n=bnormalv;
    for (int k=0;k<nt;k++)
	for (int i=0;i<nea;i++)
	{  
	    int ii,kk=ElementAdj(k,ii=i);
	    if (kk<0 || kk==k) {
		Element & K(elements[k]);
		Rd N;//=K.n(i);
		for(int j=0;j<nva;++j)
		{
		    K[Element::nvadj[i][j]].SetNormal(n,N);
		}
		    
	    }
	}
    // cout << nb << " == " << n-bnormalv << endl;
    assert(n - bnormalv <= nb );
}

}
#endif
