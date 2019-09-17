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
extern long searchMethod; //pichon
#include <map>  // Add J. Morice

#include "cassert"
#include "assertion.hpp"
#include <cstdlib>
#include <utility>
#include <limits>

#include "RefCounter.hpp"


using namespace ::std;


#include "Serialize.hpp"

#include "GQuadTree.hpp"
// definition R
namespace Fem2D  {
#include "R3.hpp"
#include "Label.hpp"
#include "HashTable.hpp"

inline int  randwalk(int ncas)
{
        const long long a = 314125421, m=  777777;
        static long xn = 19999999%m;
        if(ncas <=0) xn=19999999%m;
        long long xxn = xn;
        xn = a*xxn%m +1;
        xn %= m;
        return  (xn*ncas)/m;
}

  const double UnSetMesure=-1e+200;

inline int maxdfon(const int *dfon){ return max(max(dfon[0],dfon[1]),max(dfon[2],dfon[3]));}
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

    template<int d> inline int NumPerm(int *) {ffassert(0);}
    template<int d> inline int NumPerm1(int *) {ffassert(0);}  // num perm inverse
    template<> inline int NumPerm<1>(int *) { return 0;}
    template<> inline int NumPerm1<1>(int *) { return 0;}
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
    template<int d> inline void SetNumPerm(int n,int *p) { ffassert(0); }// a error}
    template<int d> inline void SetNumPerm1(int n,int *p) { ffassert(0); }// a error}

    template<> inline void SetNumPerm<1>(int ,int *p) { p[0]=0;} // a error}
    template<> inline void SetNumPerm<2>(int n,int *p) { p[0]=n;p[1]=1-n;} // a error}

    // build perm inverse
    template<> inline void SetNumPerm1<1>(int ,int *p) { p[0]=0;} // a error}
    template<> inline void SetNumPerm1<2>(int n,int *p) { p[0]=n;p[1]=1-n;} // a error}

    template<> inline  void SetNumPerm1<3>(int n,int *p) {
	int i=n/2, j= n%2 ? 2:1;
	p[i]=0;p[(i+j)%3]=1;p[(i+j+j)%3]=2;
	assert( n == NumPerm1<3>(p));
    }

    template<> inline   void SetNumPerm<3>(int n,int *p) {
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
  const int MaxNbDFPerNode;
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
    MaxNbDFPerElement(m.MaxNbDFPerElement) ,
    MaxNbDFPerNode(maxdfon(m.ndfon))
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
    MaxNbDFPerElement(aMaxNbDFPerElement) ,
    MaxNbDFPerNode(maxdfon(andfon))
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
  template<int d> inline  R1 ExtNormal( GenericVertex<R1> *const v[2],int const f[1])  {  static_assert ( d== 1 ,"dim =1"); return f[0]==0 ? R1(-1):R1(1);  }
  template<int d> inline  R2 ExtNormal( GenericVertex<R2> *const v[3],int const f[2])  {    static_assert ( d== 2 ,"dim=2");  return R2(*v[f[1]],*v[f[0]]).perp();  }
    // correct signe N in 3d mai 2009 (FH)
  template<int d> inline  R3 ExtNormal( GenericVertex<R3> *const v[4],int const f[3])  {  static_assert ( d== 3 ,"dim=3"); return R3(*v[f[0]],*v[f[2]])^R3(*v[f[0]],*v[f[1]]) ;  }
  template<> // pour axel exterior  Normal of surface ..
    inline  R3 ExtNormal<2>( GenericVertex<R3> *const v[3],int const f[2])  {   return R3(*v[f[0]],*v[f[1]])   ^(R3(*v[0],*v[1])^R3(*v[0],*v[2])) ;
       // R3(*v[0],*v[2])^R3(*v[0],*v[1])  ^ R3(*v[f[0]],*v[f[1]]) ;

    }// module 2 aire*l


    // Clever the orientation in case of only 1 vertex april 2019 (Hard to find !!!!! FH and PHT)
    template<int NN,typename V> struct  SwapOrient { static void SwapO(V  **w){swap(w[0],w[1]);}} ;
    template<typename V> struct  SwapOrient<1,V> { static void SwapO(V **){}} ;


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

  GenericElement & set(Vertex * v0,int * iv,int r,double mss=UnSetMesure)
  {
    for(int i=0;i<nv;++i)
      vertices[i]=v0+iv[i];
    mes=(mss!=UnSetMesure) ? mss : Data::mesure(vertices);
    lab=r;
    ASSERTION(mss==UnSetMesure && mes>0);
    return *this;
  }

    void changeOrientation() { SwapOrient<nv,Vertex>::SwapO(vertices);}

  istream & Read1(istream & f,Vertex * v0,int n)
  {
    int iv[nv]={},ir,err=0;
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

  Rd N(int i) const  { return ExtNormal<RdHat::d> (vertices,nvadj[i]);}
  RdHat PBord(int i,RdHatBord P) const   { return Data::PBord(nvadj[i],P);} // Correction FH  mars 2019 For Axel

  Rd operator()(const RdHat & Phat) const {
    Rd r= (1.-Phat.sum())*(*(Rd*) vertices[0]);
    for (int i=1;i<nv;++i)
      r+=  Phat[i-1]*(*(Rd*) vertices[i]);
    return r;
  }

 int faceOrient(int i) const
    {// def the permutatution of orient the face
	int fo =1;
	const Vertex * f[3]={&at(nvface[i][0]), &at(nvface[i][1]), &at(nvface[i][2])};
	if(f[0]>f[1]) fo = -fo,Exchange(f[0],f[1]);
	if(f[1]>f[2]) { fo = -fo,Exchange(f[1],f[2]);
	if(f[0]>f[1]) fo = -fo,Exchange(f[0],f[1]); }
	return fo;
    }

  int facePermutation(int i) const
  {// def the permutatution of orient the face
    int fo =0;
    const Vertex * f[3]={&at(nvface[i][0]), &at(nvface[i][1]), &at(nvface[i][2])};
    if(f[0]>f[1]) fo+=1,Exchange(f[0],f[1]);
    if(f[1]>f[2]) { fo+=2,Exchange(f[1],f[2]);
    if(f[0]>f[1]) fo+=4,Exchange(f[0],f[1]); }
    return fo;
  }

  int   EdgeOrientation(int i) const
    {  return 2*(&at(nvedge[i][0]) < &at(nvedge[i][1]))-1;}// return -1 or 1 FH: Change jan 2018

  R lenEdge(int i) const {ASSERTION(i>=0 && i <ne);
    Rd E=Edge(i);return sqrt((E,E));}
  R lenEdgesmax() const
    {
        R lx2 = 0;
        for (int i=0; i< ne; ++i)
            lx2 = max(lx2,lenEdge2(i));
        return sqrt(lx2);

    }
  R lenEdge2(int i) const {ASSERTION(i>=0 && i <ne);
        Rd E=Edge(i);return ((E,E));}

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

  int nt,nv,nbe,nadjnomanifold;
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
    int nbElmts() const {return nt;}
    int nbBrdElmts() const {return nbe;}
    int nbVertices() const {return nv;}
    const T & operator[](int i) const {return elements[CheckT(i)];}
  const V& operator()(int i) const {return vertices[CheckV(i)];}
  const B& be(int i) const {return borderelements[CheckBE(i)];}
  void  BoundingBox(Rd &pmin,Rd &pmax) const {pmin=Pmin;pmax=Pmax;}
  T & t(int i)  {return elements[CheckT(i)];}
  V & v(int i)  {return vertices[CheckV(i)];}
  B & be(int i) {return borderelements[CheckBE(i)];}
  const  T & t(int i)  const {return elements[CheckT(i)];}
  const  V & v(int i)  const {return vertices[CheckV(i)];}


  GenericMesh()
    : nt(0),nv(0),nbe(0),  mes(0.),mesb(0.) ,
      vertices(0),elements(0),borderelements(0),bnormalv(0),
      TheAdjacencesLink(0),BoundaryElementHeadLink(0),
      ElementConteningVertex(0), gtree(0)
  {
      
  }

  GenericMesh(const  Serialize &serialized) ;

  void set(int mv,int mt,int mbe)
  {
    assert(nt==0 && nv==0 && nbe ==0);
    nt=mt;
    nv=mv;
    nbe=mbe;
    vertices=new V[nv];
    if(nt) elements= new T[nt];
    if(nbe>0) borderelements = new B[nbe];
    assert( nt >=0 && elements);
    assert( nv >0 && vertices);

  }


  int operator()(const T & tt) const {return CheckT(&tt - elements);}
  int operator()(const T * tt) const {return CheckT(tt - elements);}
  int operator()(const V & vv) const {return CheckV(&vv - vertices);}
  int operator()(const V  * vv) const{return CheckV(vv - vertices);}
  int operator()(const B & k) const {return CheckBE(&k - borderelements);}
  int operator()(const B  * k) const{return CheckBE(k - borderelements);}
  int operator()(int it,int j) const {return operator()(elements[it][j]);}// Nu vertex j of triangle it
  int be(int it,int j) const {return operator()(borderelements[it][j]);}// Nu vertex j of triangle it

  int CheckV(int i) const { ASSERTION(i>=0 && i < nv); return i;}
  int CheckT(int i) const { ASSERTION(i>=0 && i < nt); return i;}
  int CheckBE(int i) const { ASSERTION(i>=0 && i < nbe); return i;}


  int Contening(const Vertex * vv) const{ return ElementConteningVertex[ vv  - vertices];}
  void BuildAdj();
  void BuildBoundaryElementAdj();  // Add J. Morice function that give the TheAdjacencesSurfaceLink :: Version avec un manifold
  void BuildBoundaryElementAdj(const int &nbsurf, int* firstDefSurface, int* labelDefSurface, int* senslabelDefSurface); // version avec plusieurs vari�t�s
  
  //template<>
  void SameVertex(const double precis_mesh, Vertex *v, Element *t, int nv, int nt, int *ind_nv_t, int *Numero_Som, int &new_nv, bool removeduplicate);
  void VertexInElement( Vertex *vertice, Element *list, int &nv, int *(&ind_nv), int nt, int *ind_nt, int *(&old2new) );
  void clean_mesh(const double precis_mesh, int &nv, int &nt, int &nbe, V *(&v), T *(&t), B *(&b), bool removeduplicate, bool rebuildboundary, int orientation);
 
  void Buildbnormalv();
  void BuildBound();
  void BuildjElementConteningVertex();
  void BuildGTree() {if(gtree==0)  gtree=new GTree(vertices,Pmin,Pmax,nv);}
  DataFENodeDF BuildDFNumbering(int dfon[NbTypeItemElement],int nbequibe=0,int *equibe=0) const ;
    DataFENodeDF BuildDFNumbering(int ndfv,int ndfe,int ndff,int ndft,int nbequibe=0,int *equibe=0) const
  { int dfon[NbTypeItemElement]={ndfv,ndfe,ndff,ndft};
    return  BuildDFNumbering(dfon,nbequibe,equibe);
  }
  int nElementonB(int k,int j) const // add v4 F.H 12/2018 (not sure for internal boundary !!!)
    { int kk= TheAdjacencesLink[nea*k+j]; return (kk>=0) && (kk%nea  != k) ? 2 : 1;}

  int ElementAdj(int k,int &j) const  {
    int p=TheAdjacencesLink[nea*k+j];
    if(p>=0) j=p%nea;
    return p>=0 ? p/nea: -1-j;}// modif FH. to change the code of copule k,kadj on border element..
    //  correct bug of 23/05/2013 : 1  dof on RT0 3d...



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

  int BoundaryElement(int bbe,int & ItemInK) const {
    int i= BoundaryElementHeadLink[bbe];
    ItemInK = i%nea;
    return i/nea;}

  // Add J. Morice
  template<int N,int M>
  SortArray<int,N> itemadjs(const int (* const  nu )[N],int k,int i, int *sens) const
  {
    int nnv[N];
    const B & K(borderelements[CheckBE(k)]);
    ASSERTION(i>=0 && i <M);
    for (int j=0;j<N;++j){
      nnv[j] = operator()(K[nu[i][j]]);
    }
     return SortArray<int,N>(nnv,sens);
  }

  SortArray<int,B::nva> items(int k,int i,int *sens) const
  {
    return itemadjs<B::nva,B::nv>(B::nvadj,k,i,sens);
  }


  template<int N,int M>
  SortArray<int,N> iteme(const int (* const  nu )[N],int k,int i,int *psens=0)  const
  {
    int nnv[N];
    const Element & K(elements[CheckT(k)]);
    ASSERTION(i>=0 && i <M);
    for (int j=0;j<N;++j){
      nnv[j] = operator()(K[nu[i][j]]);
    }

    return SortArray<int,N>(nnv,psens);
  }

  SortArray<int,B::nv> itemadj(int k,int i,int *psens=0) const
  {
    return iteme<B::nv,T::nea>(T::nvadj,k,i,psens);
  }

  SortArray<int,B::nv> itembe(int k,int *psens=0) const
  {
    int nnv[B::nv];
    const B & K(borderelements[CheckBE(k)]);

    for (int j=0;j<B::nv;++j){
      nnv[j] = operator()(K[j]);
    }

    return SortArray<int,B::nv>(nnv,psens);
  }

  //  const Element * Find(const Rd & P) const ;
  const Element * Find(Rd P, RdHat & Phat,bool & outside,const Element * tstart=0) const
  {return EF23::Find<GMesh>(*this,this->gtree,P,Phat,outside,tstart);}

  R mesure(){ return mes;}
  R bordermesure(){ return mesb;}
  virtual ~GenericMesh() {

    delete [] ElementConteningVertex;
    delete [] TheAdjacencesLink;
    delete [] BoundaryElementHeadLink;
    if(nt>0) delete [] elements;
    if(nbe>0) delete [] borderelements;
    delete [] vertices;
    delete [] bnormalv;
    if(gtree) delete gtree;
      ElementConteningVertex=0;
      TheAdjacencesLink=0;
      BoundaryElementHeadLink=0;
      borderelements=0;
      elements=0;
      vertices=0;
      bnormalv=0;
      gtree=0;
      nt=(0);
      nv=(0);
      nbe=(0);
      mes=(0.);
      mesb=(0.);
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
    int lerr[10]={};
  if(!ElementConteningVertex) ElementConteningVertex = new int[nv];

    for(int i=0;i<nv;++i)
	ElementConteningVertex[i]=-1;

    for (int k=0;k<nt;++k)
	for (int i=0;i<nkv;++i)
	    ElementConteningVertex[operator()(elements[k][i])]=k ;
    int kerr=0;
    for(int i=0;i<nv;++i)
	if (ElementConteningVertex[i]<0)
	    if(kerr<10) lerr[kerr++]=i;
    if(kerr)
      {
	cerr << " Fatal error: some vertex are not at least in one element  \n       :   " ;
	for(int i=0;i<kerr;++i)
	    cerr << " " << lerr[i];
	cerr << endl;
      }
    ffassert(kerr==0);//  Sure Error.

}
    template<typename T,typename B,typename V>
    void GenericMesh<T,B,V>::BuildAdj()
    {
        int verb = verbosity ;
        if(TheAdjacencesLink!=0) return ;//  already build ...
        TheAdjacencesLink = new int[nea*nt];
        BoundaryElementHeadLink = new int[nbe];
        HashTable<SortArray<int,nva>,int> h(nea*nt,nv);
        int nk=0,nba=0;
        int err=0;
        if(verbosity>5)
            cout << "   -- BuildAdj:nva= " << nva << " " << nea << " "<< nbe << endl;
        nadjnomanifold=0;
        for (int k=0;k<nt;++k)
            for (int i=0;i<nea;++i) {
                int sens;
                SortArray<int,nva> a(itemadj(k,i,&sens));//  warning the face of tet given interieon normal FH.
                if(verbosity>99) cout <<nk << " T "<< k << "### "   << " item(k,i)= " << itemadj(k,i) << " a= " << a << " k " << k << " i " << i << endl;
                typename HashTable<SortArray<int,nva>,int>::iterator p= h.find(a);
                if(!p) {
                    h.add(a,nk);
                    TheAdjacencesLink[nk]=-1;
                    nba++;
                }
                else {
                    if(p->v<0 ) {// no manifold TO DO
                        // clean adj
                        int nk1=-1-p->v;
                        int nk2= TheAdjacencesLink[nk1];// next
                        if(nk2>=0) { // firt time remove existing link ...
                            nadjnomanifold++;
                            TheAdjacencesLink[nk1]=nk  ; // inserting between nk1 and nk2
                            TheAdjacencesLink[nk]=nk2 ;
                            //  on no manifold border .
                            if( verb>99 ) cout << " Border manifold " << k << " "<< i << " ::  " << itemadj(k,i)<< " ::  " << nk1/nea << " " << nk2/nea
                            << " :: " <<TheAdjacencesLink[nk1]/nea << '.' << TheAdjacencesLink[nk]/nea << '.' << TheAdjacencesLink[nk2]/nea << endl;
                            //nba--;
                            // no manifold TO DO if false
                        }
                    }
                    else {
                        // cout << " test p->v "<< p->v << endl;
                        TheAdjacencesLink[nk]=p->v;
                        TheAdjacencesLink[p->v]=nk;
                        p->v=-1-nk;
                        nba--;
                    }
                   
                    
                }
                ++nk;
            }
        int nbordnomanifold=0;
        if(nadjnomanifold)
        {
            
            // remove all adj of no manifold border
            int k3= nt*nea;
            for (int p=0; p< k3; ++p )
            {
                int pp=TheAdjacencesLink[p];
                if(pp>0 && TheAdjacencesLink[pp]>=0 && TheAdjacencesLink[pp] != p )
                    { // border of no manifold
                        if( verb>19 ) cout << "  -- " <<  p/nea << " " ;
                        ++nbordnomanifold;
                       // remove link ...  put -2 in all list ...
                        TheAdjacencesLink[p]=-2;
                        while (pp>=0 && TheAdjacencesLink[pp]>=0 )
                        {
                            if( verb>19 ) cout <<  pp/nea << " " ;
                            int ppp=pp;
                            pp=TheAdjacencesLink[pp];
                            TheAdjacencesLink[ppp]=-2;// break the list ...
                        }
                        if( verb>19 ) cout << " . " << endl;
                    }
            }
            
        }
        if(verbosity&& nadjnomanifold) cerr << "  --- Warning manifold obj nb:" << nbordnomanifold<< " adj "<< nadjnomanifold << " of  dim =" << T::RdHat::d << endl;
        int kerr=0,kerrf=0,nbei=0,fwarn=0;
        map<pair<int,int>,pair<int,int> > mapfs;
        int uncorrect =0, nbchangeorient=0;
        for(int step=0; step<2; ++step)
        {
            for (int ke=0;ke<nbe;++ke)
            {
                int sens,s=0,ss=0;
                SortArray<int,nva> a(itembe(ke,&sens));

                typename HashTable<SortArray<int,nva>,int>::iterator p= h.find(a);

                if(verbosity>99) cout << "B " << ke << " ### "   << " item(k,i)= " << itembe(ke) << " a= " << a << endl;

                if(!p) { err++;
                    if(err==1) cerr << "Err  Border element not in mesh \n";
                    if (err<10)  cerr << " \t " << ke << " " << a << endl;
                }
                else
                {
                    int nk = p->v <0 ? -p->v-1 : p->v;
                    int nkk= TheAdjacencesLink[nk];
                    if( nkk>=0)
                    {
                        nbei ++;
                        // choise le bon .. too get the correct normal
                        int k= nk/nea, e=nk%nea;
                        int kk= nkk/nea, ee=nkk%nea;
                        itemadj(k,e,&s); if(verbosity>15) cout << " item(k,e)= " << itemadj(k,e) <<  " k " << k << " e " << e << " s " << s << endl;
                        itemadj(kk,ee,&ss); if(verbosity>15) cout << " item(kk,ee)= " << itemadj(kk,ee) << " kk " << kk << " ee " << ee << " ss " << ss << endl;
                       //assert(s && ss && s== -ss);
                         if (!(s && ss && s== -ss)) {
                            cerr << " Bad orientation: The adj border element  defined by [ " << itemadj(k,e) << " ]  is oriented in the same direction in element "
                        << k << " and in the element " << kk << " ****** bug in mesh construction? orientation parameter? "<< endl;
                        ffassert(0);
                        }
                        if( sens == s) {swap(nk,nkk);swap(k,kk);} //  autre cote
                        //  verif sens normal
                        int regk= elements[k].lab;
                        int regkk= elements[kk].lab;

                        if(regk != regkk)
                        { //  verif orientation  interior normal
                            int sr =1;
                            if(step==0)
                            {
                                if(regk>regkk)
                                    mapfs[make_pair(regkk,regk)].second++;// increme le + grand
                                else
                                    mapfs[make_pair(regk,regkk)].first++;// increme le + grand
                            }
                            else { // correct the sens of face
                                if(regk>regkk) sr=-1,swap(regk,regkk);
                                // change orientation de du plus petit nombre
                                map<pair<int,int>,pair<int,int> >::iterator p= mapfs.find(make_pair(regk,regkk));
                                int nfk = p->second.first;
                                int nfkk =p->second.second;
                                //  event regkk
                                bool change=0;
                                if( (nfk < nfkk) && sr == 1 ) change=1;
                                if( (nfk > nfkk) && sr == -1 ) change=1;
                                if( change)
                                {
                                    nbchangeorient++;
                                    borderelements[ke].changeOrientation();
                                    if(verbosity>2) cout << " changeOrientation face:"<< ke << endl;
                                    swap(nk,nkk);
                                    swap(k,kk);
                                }
                            }

                        }
                        
                    }
                    else
                    {//  verif if the face is in correct sens ...
                      int sk,k= nk/nea, e=nk%nea; //   same oreintatio ???
                      itemadj(k,e,&sk);
                       if( sk != sens &&  (step == 0) )
                       {
                           fwarn++;
                           if( verbosity>4 && fwarn < 10) cout << "  --  warning true  boundary element "<<ke << " is no in correct orientation " <<endl;
                           borderelements[ke].changeOrientation();
                           nbchangeorient++;
                       }
                    }

            


                BoundaryElementHeadLink[ke] = nk;
                int sk;
                int tt=BoundaryElementHeadLink[ke]/nea;
                int ee=BoundaryElementHeadLink[ke]%nea;
                itemadj(tt,ee,&sk);
                if(!(itemadj(tt,ee)==a)) {
                    if(kerrf< 10) cout << " err face  not in element "<< ke << " =="
                        << tt << " " << ee << endl;
                    kerrf++;
                }
            }
        }

        for(map<pair<int,int>,pair<int,int> >::iterator p=mapfs.begin(); p != mapfs.end(); ++p)
        {
            if( p->second.first && p->second.second)
            {
                if(verbosity>2 && step==0)  cout << " error in orientation of internal face beetwen region "
                    << p->first.first << " , " << p->first.second << " to no zero value "
                    << p->second.first << "  " << p->second.second << endl;
                uncorrect++;
            }
        }
        if(uncorrect==0) break;
    }
    if( nbchangeorient && verbosity>2) cout << " Warning change orientation of " << nbchangeorient << " faces \n";
    if( fwarn && verbosity>2) cout << " Warning error in boundary oriention  " << fwarn  << " faces \n";
    if( kerr || kerrf ) {
        cout << " Erreur in boundary orientation  bug in mesh or bug in ff++ "  << kerr  << " / " <<nbei  << "\n\n";
        cout << "  or Erreur in face    "  << kerrf  << " / " <<nbei  << "\n\n";

        }
        ffassert(kerr==0 && kerrf==0);
        ffassert(err==0);
        int na= h.n;
        if(verbosity>1)
        {
            cout << "  -- BuildAdj: nb Elememt " << nt << " nb vertices " << nv << endl;
            cout << "             : nb adj  = "<< na << " on border " << nba << " nea = " << nea << " nva = " << nva
                 << " nb no manifold border " << nadjnomanifold << endl;
            if(nea==2)
                cout << " Const d'Euler: " << nt - na + nv << endl;
            else
                cout << endl;
        }
}


template<typename T,typename B,typename V>
void GenericMesh<T,B,V>::BuildBoundaryElementAdj()
{
    // Return in TheBorderElementAjacencesLink
    //  if exist a link :: sign(nk_link)*(nk_link+1)
    //  else            :: sign(nk)*(nk)

    int *TheBoundaryElementAdjacencesLink = new int[B::nea*nbe];
    HashTable<SortArray<int,B::nva>,int> h(B::nea*nbe,nv);
    int nk=0;
    int err=0,errm=0;
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
                TheBoundaryElementAdjacencesLink[nk] = sens*(nk+1)   ;  // sens;
            }
            else
            {
                ASSERTION(p->v>=0);
                if( sens*TheBoundaryElementAdjacencesLink[p->v] > 0 ){

                    B & K(borderelements[CheckBE(k)]);
                    // Bug before is here nea : nb of border element for adj : , nva  = nb of vertex of border element for adl

                    cout << " The adj border element  defined by [ " ;
                    for(int ia=0; ia<nva; ++ia)
                        cout <<  operator()(K[B::nvadj[i][ia]])+1<< " ";
                    cout << " ]  is oriented in the same direction in element " << k+1 <<
                    " and in element "<<  1+(p->v/B::nea) << endl;
                    err++;

                }
                if( abs(TheBoundaryElementAdjacencesLink[p->v]) != 1+p->v ){

                    B & K(borderelements[CheckBE(k)]);
                    if(errm<10)
                    {
                        cout << " The adj border element defined by vertex [ ";
                        for(int ia=0; ia<nva; ++ia)
                            cout <<  operator()(K[B::nvadj[i][ia]])+1<< " ";
                        cout << " ] is belong to the three border elements ::"
                        << 1+(p->v)/B::nea <<", "<< k+1 <<" and "<< 1+(abs(TheBoundaryElementAdjacencesLink[p->v])-1)/B::nea << endl;
                    }
                    errm++;
                }
                if( errm)  cout << " The border "<<B::RdHat::d <<  "d  " << " is none  manifold" << endl;
                if( err)  cout << " The border "<<B::RdHat::d <<  "d  " << " is badly oriented !!!!  " << endl;
                ffassert(err==0 && errm==0);
                TheBoundaryElementAdjacencesLink[nk]= TheBoundaryElementAdjacencesLink[p->v];
                TheBoundaryElementAdjacencesLink[p->v]= sens*(nk+1);

            }
            if( err > 10 )
                exit(1);
            nk++;
        }

    assert(err==0);
    delete [ ] TheBoundaryElementAdjacencesLink;
    if(verbosity) cout << "number of adjacents edges " << nk << endl;
}


template<typename T,typename B,typename V>
void GenericMesh<T,B,V>::BuildBoundaryElementAdj(const int &nbsurf, int* firstDefSurface, int* labelDefSurface, int* senslabelDefSurface)
{

   // Return in TheBoundaryElementAdjacences
  //  if exist a link :: sign(nk_link)*(nk_link+1)
  //  else            :: sign(nk)*(nk)


  for(int isurf=0; isurf < nbsurf; isurf++){

    //######################################
    // Trop operations if ===> a changer

    int nbe_surf=0; // number in the surface

    for(int k=0; k<nbe; k++){
      B & K(borderelements[CheckBE(k)]);
      int label = K.lab;
      for(int iii=firstDefSurface[isurf]; iii< firstDefSurface[isurf+1];iii++)
	if( label == labelDefSurface[iii] ) nbe_surf++;
    }

    int facek=0;
    int *surf_be=new int[nbe_surf];
    int *orientation_surf_be=new int[nbe_surf];
    for(int k=0; k<nbe; k++){
      B & K(borderelements[CheckBE(k)]);
      int label = K.lab;
      for(int iii=firstDefSurface[isurf]; iii< firstDefSurface[isurf+1];iii++)
	if( label == labelDefSurface[iii] ) {
	  surf_be[facek]=k ;
	  orientation_surf_be[facek]=senslabelDefSurface[iii];
	  facek++;
	}
    }

    //######################################

    int *TheBoundaryElementAdjacencesLink = new int[B::nea*nbe_surf];
    HashTable<SortArray<int,B::nva>,int> h(B::nea*nbe_surf,nv);
    int nk=0;
    int err=0;
    int sens;
      if(verbosity>4) cout << "BuildBoundaryElementAdj: nea/nva " << B::nea << " "  << B::nva << endl;
    for (int k=0;k<nbe_surf;++k)
      for (int i=0;i<B::nea;++i)
	{
	  SortArray<int,B::nva> a(items( surf_be[k],i,&sens));
	  sens=sens*orientation_surf_be[k];
	  typename HashTable<SortArray<int,B::nva>,int>::iterator p= h.find(a);
	  if(!p)
	    {
	      h.add(a,nk);
	      TheBoundaryElementAdjacencesLink[nk]=sens*(nk+1);
	      // nk est un nombre locale depend de la surfaces choisie
	      // element du bord est donn�e par ::  surf_be[nk/3];
	      // arrete corespondante locale de l'element :: nk%3;
	    }
	  else
	    {

	      ASSERTION(p->v>=0);
	      if( sens*TheBoundaryElementAdjacencesLink[p->v] > 0 ){

		B & K(borderelements[CheckBE(surf_be[k])]);
		int firstVertex  =  operator()(K[B::nvadj[0][i]])+1;
		int secondVertex =  operator()(K[B::nvadj[1][i]])+1;
		cout << " The edges, defined by vertex is " << firstVertex << "-" << secondVertex <<  ", is oriented in the same direction in element " << surf_be[k]+1 <<
		  " and in element "<<  1+surf_be[(p->v/B::nea)] << endl;

		err++;
	      }

	      if( abs(TheBoundaryElementAdjacencesLink[p->v]) != 1+p->v ){

		B & K(borderelements[CheckBE(k)]);
		int firstVertex  =  operator()(K[B::nvadj[0][i]])+1;
		int secondVertex =  operator()(K[B::nvadj[1][i]])+1;
		cout << " The edges defined by vertex is " << firstVertex << "-" << secondVertex << "belong to the three border elements ::"
		     << 1+surf_be[(p->v)/B::nea] <<", "<< surf_be[k]+1 <<" and  "<< 1+surf_be[(abs(TheBoundaryElementAdjacencesLink[p->v])-1)/B::nea] << endl;
		cout << " The "<< isurf+1 << " Surface contains these edges is not a manifold" << endl;
		err++;
		assert(err==0);
	      }


	      TheBoundaryElementAdjacencesLink[nk]   = TheBoundaryElementAdjacencesLink[p->v];;
	      TheBoundaryElementAdjacencesLink[p->v] = sens*(nk+1);

	    }
	  if( err > 10 )
	    exit(1);
	  nk++;
	}

    assert(err==0);
    delete [] surf_be;
    delete [] orientation_surf_be;
    delete [ ] TheBoundaryElementAdjacencesLink;
    if(verbosity) cout << "number of adjacents edges " << nk << endl;
  }
}

    
template<typename T,typename B,typename V>
void GenericMesh<T,B,V>::VertexInElement(V *vertice, T *list, int &nv, int *(&ind_nv), int nt, int *ind_nt, int *(&old2new) ) {
    
    // new map and list vertices after clean multiple elements

    int map[nv];
    int indv[nv];
  
    int takev[nv], takenewv[nv] ;
    for (int i=0;i<nv;i++) {
        indv[i]=-1;
        map[i]=i; // identidy by default
        takev[i]=-1;
        takenewv[i]=-1;
    }
    int ik,oldik, np=0;
    // loop on new elements to check used vertice
    for (int it=0;it<nt;it++) {
       int itri=ind_nt[it];
        const T &K(list[itri]);
        // the vertices in new triangles (vertice in new num
        for (int j=0;j<T::nv;j++){
            ik=old2new[&(K[j]) - vertice]; // new indi v
            oldik=&(K[j]) - vertice;  // origin ind v
          
            if ( takenewv[ik]==-1){
                map[oldik]=np;
                indv[np]=ind_nv[ik] ;
                takenewv[ik]=np;
                np++;
            }
           map[oldik]=takenewv[ik];
    
        }
    }
    if(verbosity>5)  cout << " real used vertice:" << np << endl;
 
    for(int i=0;i<nv;i++) {
        ind_nv[i]=indv[i];
        old2new[i]=map[i];
    }
    nv=np;
        
}
    
 
    
// output int *ind_nv_t, int *old2new, int &new_nv
template<typename T,typename B,typename V>
void GenericMesh<T,B,V>::SameVertex(const double precis_mesh, V *vertice, T *element, int nv, int nt, int *ind_nv_t, int *old2new, int &new_nv, bool removmult) {

    if (verbosity > 2)
    cout << "clean mesh: remove multiple vertices, elements, border elements and check border elements " << endl;
    double precispt=0;
    // remove the multiple vertices
    // 1/ compute bmin, bmax, hmin with the criteria precis_mesh
    // 2/ build a gtree to extract the simple vertices (not multiple)

    double hmin = 1e10;
    Rd bmin, bmax;
    if (verbosity > 2) {cout << " BuilBound " << endl;}
    //Mesh::Vertex ?
    bmin.x = vertice[0].x;
    bmin.y = vertice[0].y;
    bmin.z = vertice[0].z;
    bmax.x = bmin.x;
    bmax.y = bmin.y;
    bmax.z = bmin.z;

    if (verbosity > 1) {cout << " determination of bmin and bmax" << endl;}

    for (int i = 1; i < nv; i++) {
        bmin.x = min(bmin.x, vertice[i].x);
        bmin.y = min(bmin.y, vertice[i].y);
        bmin.z = min(bmin.z, vertice[i].z);
        
        bmax.x = max(bmax.x, vertice[i].x);
        bmax.y = max(bmax.y, vertice[i].y);
        bmax.z = max(bmax.z, vertice[i].z);
    }

    double longmini_box;
    longmini_box = pow(bmax.x - bmin.x, 2) + pow(bmax.y - bmin.y, 2) + pow(bmax.z - bmin.z, 2);
    longmini_box = sqrt(longmini_box);

    if (verbosity > 1) {
        cout << " bmin := " << bmin.x << " " << bmin.y << " " << bmin.z << endl;
        cout << " bmax := " << bmax.x << " " << bmax.y << " " << bmax.z << endl;
        cout << " box volume :=" << longmini_box << endl;
    }

    if (precis_mesh < 0)
        precispt = longmini_box * 1e-7;
    else
        precispt = precis_mesh;

    // determination de hmin
    for (int i = 0; i < nt ; i++) {
        const T &K(element[i]);
        double longedge=0.;
        int iv[T::nea];
        
        for (int j = 0; j < T::nea ; j++)
            iv[j] = operator()(K[j]);
        for (int j = 0; j < T::nea ; j++)
            for (int k = j + 1; k < T::nea; k++) {
                int &i1 = iv[j];
                int &i2 = iv[k];
                longedge = pow(vertice[i1].x - vertice[i2].x, 2)
                + pow(vertice[i1].y - vertice[i2].y, 2)
                + pow(vertice[i1].z - vertice[i2].z, 2);
                longedge = sqrt(longedge);
                if (longedge > precispt) hmin = min(hmin, longedge);
            }
        
        }

    
        if (verbosity > 5) {
            cout << "    longmini_box" << longmini_box << endl;
            cout << "    hmin =" << hmin << endl;
        }

        assert(hmin < longmini_box);
        if (verbosity > 5)
            cout << "    Norme2(bmin-bmax)=" << Norme2(bmin - bmax) << endl;
    
        // assertion pour la taille de l octree
        assert(hmin > Norme2(bmin - bmax) / 1e9);
    
    
        double hseuil = hmin / 10.;
        if (verbosity > 3)
            cout << "    hseuil=" << hseuil << endl;

        V *vv = new V[nv];
        GTree *gtree = new GTree(vv, bmin, bmax, 0);
    
        if (verbosity > 2) {
            cout << "  -- taille de la boite " << endl;
            cout << "\t" << bmin.x << " " << bmin.y << " " << bmin.z << endl;
            cout << "\t" << bmax.x << " " << bmax.y << " " << bmax.z << endl;
        }
    
        // creation of octree
        for (int i = 0; i < nv; i++) {
            const V &vi(vertice[i]);
            V *pvi = gtree->ToClose(vi, hseuil);
            if (!pvi) {
                vv[new_nv].x = vi.x;
                vv[new_nv].y = vi.y;
                vv[new_nv].z = vi.z;
                vv[new_nv].lab = vi.lab;
                ind_nv_t[new_nv] = i;
                old2new[i] = new_nv;
                gtree->Add(vv[new_nv]);
                new_nv++;
            }
            // treatment of multiple vertice made in sameElement
            else
                old2new[i] = pvi - vv;
                
        }
        delete gtree;
        delete [] vv;
    
}
    
    
template<typename T,typename B,typename V>
void GenericMesh<T,B,V>::clean_mesh(double precis_mesh, int &nv, int &nt, int &nbe, V *(&v), T *(&t), B *(&b), bool removeduplicate, bool rebuildboundary, int orientation) {

  
    // array for the index of new vertices, element, borderelement in the old numbering
    int new_nv=0, new_nt=0, new_nbe=0;
   // double mes = 0, mesb = 0;
    int *ind_nv=new int[nv];
    int *ind_nt=new int[nt];
    int *ind_nbe=new int[nbe];
    // mapping old to new vertices index
    int *old2new=new int[nv];

    // clean multiples vertice
    SameVertex(precis_mesh, v, t, nv, nt, ind_nv, old2new, new_nv, removeduplicate);
   
    if(!removeduplicate)
        nv=new_nv;
    
    // clean multiple elements and border elements
    SameElement(v, t, nt, ind_nt, old2new, new_nt, removeduplicate);
    SameElement(v, b, nbe, ind_nbe, old2new, new_nbe, removeduplicate);
    
    // if removeduplicate, must recheck the vertice -> could have some vertex are not at least in one element
    
    if(removeduplicate)
        VertexInElement(v, t, nv, ind_nv, new_nt, ind_nt, old2new);
  
    V *vv;
    vv=new V[nv];
    
    for(int i=0;i<nv;i++) {
        int iv=ind_nv[i];
        vv[i].x=v[iv].x;
        vv[i].y=v[iv].y;
        vv[i].z=v[iv].z;
        vv[i].lab=v[iv].lab;
    }
  
  
    double mes=0., mesb=0.;
    nt=new_nt;
    T *tt=new T[nt];
    for(int i=0;i<nt;i++) {
        int &it=ind_nt[i];
        const T &K(t[it]);
        mes=K.mesure();
        int iv[T::nv];
        for (int j = 0; j < T::nea ; j++) {
            iv[j] =  old2new[ &(K[j]) - v];
            assert(iv[j] >= 0 && iv[j] < nv);
        }
        if (orientation<0)
            swap(iv[1], iv[2]);
        
        tt[i].set(vv, iv, K.lab);
        mes+=tt[i].mesure();
    }
    
    nbe=new_nbe;
    B *bb=new B[nbe];
    for(int i=0;i<nbe;i++) {
        int &ie=ind_nbe[i];
        const B &K(b[ie]);
        int iv[B::nv];
        for (int j = 0; j < B::nea ; j++) {
            iv[j] =  old2new[ &(K[j]) - v];
            assert(iv[j] >= 0 && iv[j] < nv);
        }
        if (orientation<0)
            swap(iv[(B::nv)-2], iv[(B::nv)-1]);
        bb[i].set(vv, iv, K.lab);
        mesb+=bb[i].mesure();
    }
 
  
    if(verbosity>2)
        cout << "after clean mesh, nv = " << nv << " nt = " << nt << " nbe = " << nbe << endl;
   
    if (mes<0) {
        cerr << " Error of mesh orientation , current orientation = " << orientation << endl;
        cerr << " mesure element mesh = " << mes << endl;
        cerr << " mesure border mesh = " << mesb << endl;
    }
  
    delete [] v;
    delete [] t;
    delete [] b;
    v=vv;
    t=tt;
    b=bb;
    delete []ind_nv;
    delete []ind_nt;
    delete []ind_nbe;
    delete []old2new;

}


    
    
    
    
    

// template<typename T,typename B,typename V>
// void GenericMesh<T,B,V>::BuildBoundaryElementAdj_V2(const int &nbsurf,int *firstDefSurface, int *labelDefSurface, int *senslabelDefSurface)
// {
// //   assert(firstDefSurface.N() == nbsurf+1);
// //   assert(labelDefSurface.N() == firstDefSurface[nbsurf]);
// //   assert(senslabelDefSurface.N() == firstDefSurface[nbsurf]);

//   // determination des labels des surfaces
//   map<int, int> maplabel;
//   int numero_label=0;
//   for(int ii=0; ii< firstDefSurface[nbsurf]; ii++){
//     map<int,int>::const_iterator imap=maplabel.find( abs(labelDefSurface[ii]) );
//     //cout << "K.lab= " << K.lab << endl;
//     if(imap == maplabel.end()){
//       maplabel[ abs(labelDefSurface[ii]) ] = numero_label;
//       numero_label = numero_label+1;
//     }
//   }

//   int *nbe_label=new int[numero_label];
//   for(int ii=0; ii< numero_label; ii++) nbe_label[ii] = 0;
//   for(int k=0; k<nbe; k++){
//     B & K(borderelements[CheckBE(k)]);
//     map<int,int>::const_iterator imap=maplabel.find( K.lab );

// //  if(imap == maplabel.end()){
// //       printf("The label %d given for Definition of different surface is not in the border element mesh\n",K.lab);
// //       exit(1);
// //     }
// //     else{
//     nbe_label[(*imap).second]++;
//     //    }
//   }

//   int all_nbe_label=0;
//   for(int k=0; k<numero_label; k++){
//     all_nbe_label=all_nbe_label+nbe_label[k];
//   }
//   /*
//     if(all_nbe_label != nbe){
//     cerr << "some element in the border element are not references in the Surface description" << endl;
//     exit(1);
//     }
//     assert(all_nbe_label == nbe);  // autrement cela veut dire que certain element du bord n'ont pas �t� mis dans le descriptif
//   */
//   int *organisation_be_label;
//   organisation_be_label = new int[all_nbe_label];

//   int *count_nbe_label =new int[numero_label];
//   int *debut_nbe_label =new int[numero_label+1];

//   for(int ii=0; ii< numero_label; ii++)
//     count_nbe_label[ii] =0;

//   debut_nbe_label[0]=0;
//   for(int ii=1; ii< numero_label; ii++)
//     debut_nbe_label[ii] = debut_nbe_label[ii-1]+nbe_label[ii-1];
//   debut_nbe_label[numero_label] = all_nbe_label;

//   for(int k=0; k<nbe; k++){
//     B & K(borderelements[CheckBE(k)]);
//     map<int,int>::const_iterator imap=maplabel.find( K.lab );
//     assert(imap != maplabel.end());
//     organisation_be_label[ debut_nbe_label[(*imap).second] + count_nbe_label[(*imap).second] ] = k ;
//     count_nbe_label[(*imap).second ]++;
//   }

//   for(int ii=0; ii< numero_label; ii++)
//     assert( count_nbe_label[ii] == nbe_label[ii] );

//   delete [] count_nbe_label;

//   for(int isurf=0; isurf < nbsurf; isurf++){

//     int nbe_surf=0; // number in the surface
//     for( int iii=firstDefSurface[isurf]; iii< firstDefSurface[isurf+1];iii++ ){
//       map<int,int>::const_iterator imap=maplabel.find( abs(labelDefSurface[iii]) );
//       nbe_surf=nbe_surf+nbe_label[ (*imap).second ];
//     }

//     // assert(TheBoundaryElementAdjacencesLink==0); plus tard
//     int *TheBoundaryElementAdjacencesLink = new int[B::nea*nbe_surf];
//     HashTable<SortArray<int,B::nva>,int> h(B::nea*nbe_surf,nv);
//     int nk=0;
//     int err=0;
//     int sens;

//     int count_sbe;
//     int *surf_be = new int[nbe_surf];

//     count_sbe=0;
//     for( int iii=firstDefSurface[isurf]; iii< firstDefSurface[isurf+1];iii++ ){
//       map<int,int>::const_iterator imap=maplabel.find( abs(labelDefSurface[iii]) );

//       for( int jjj= debut_nbe_label[(*imap).second]; jjj < debut_nbe_label[(*imap).second+1]; jjj++ ){
// 	int k=organisation_be_label[jjj];
// 	surf_be[count_sbe] = k;
// 	count_sbe++;

// 	for (int i=0;i<B::nea;++i)
// 	  {
// 	    SortArray<int,B::nva> a(items( k,i,&sens));
// 	    sens=sens*senslabelDefSurface[iii];
// 	    typename HashTable<SortArray<int,B::nva>,int>::iterator p= h.find(a);
// 	    if(!p)
// 	      {
// 		h.add(a,nk);
// 		TheBoundaryElementAdjacencesLink[nk] = sens*(nk+1);
//  	      }
// 	    else
// 	      {

// 		ASSERTION(p->v>=0);
// 		if( sens*TheBoundaryElementAdjacencesLink[p->v] > 0){

// 		  B & K(borderelements[CheckBE(k)]);
// 		  int firstVertex  =  operator()(K[B::nvadj[i][0]])+1;
// 		  int secondVertex =  operator()(K[B::nvadj[i][1]])+1;
// 		  cout << " The edges, defined by vertex is " << firstVertex << "-" << secondVertex << ", is oriented in the same direction in element " << k+1 <<
// 		    " and in element "<<  1+surf_be[(p->v/B::nea)] << endl;
// 		  err++;
// 		}

// 		if( abs(TheBoundaryElementAdjacencesLink[p->v]) != 1+p->v ){

// 		  B & K(borderelements[CheckBE(k)]);
// 		  int firstVertex  =  operator()(K[B::nvadj[i][0]])+1;
// 		  int secondVertex =  operator()(K[B::nvadj[i][1]])+1;
// 		  cout << " The edges defined by vertex is " << firstVertex << "-" << secondVertex << "belong to the three border elements ::"
// 		       << 1+surf_be[(p->v)/B::nea] <<", "<< surf_be[k]+1 <<" and  "<< 1+surf_be[(abs(TheBoundaryElementAdjacencesLink[p->v])-1)/B::nea] << endl;
// 		  cout << " The "<< isurf+1 << " Surface contains these edges is not a manifold" << endl;
// 		  err++;
// 		  assert(err==0);
// 		}

// 		TheBoundaryElementAdjacencesLink[nk]=TheBoundaryElementAdjacencesLink[p->v];
// 		TheBoundaryElementAdjacencesLink[p->v]=sens*(nk+1);
// 	      }

// 	    if( err > 10 )
// 	      exit(1);
// 	    nk++;
// 	  }
//       }
//     }

//     assert(err==0);
//     delete [ ] TheBoundaryElementAdjacencesLink;
//     delete [ ] surf_be;
//     if(verbosity) cout << "number of adjacents edges " << nk << endl;
//   }

//   delete [] organisation_be_label;
//   delete [] debut_nbe_label;
//   delete [] nbe_label;
// }


template<typename T,typename B,typename V>
DataFENodeDF GenericMesh<T,B,V>::BuildDFNumbering(int ndfon[NbTypeItemElement],int nbequibe,int *equibe) const
{
/*  Numbering in 1d, 2d ok if no manifold object
 Warning:   wrong in 3d if no manifold object but improble.
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
    unsigned int tinfty=std::numeric_limits<unsigned int>::max()   ;

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

      int of = Th.nv+10;// Modif FH 28/05/2013
      int ndim[NbTypeItemElement]={0,0,0,0};
      NbOfDF=0;
      {
	HashTable<Key,int> h(nnodeK*nt,of+nt);
	int  nbmaxeq = 1+nnodeK*nbequibe;
	int  nbhasheq = nbequibe ? of+nt : 1;
	HashTable<Key,Key> equi(nbmaxeq,nbhasheq);
	  //  constuction of item translation for
	  if(verbosity>9)
	  cout << " nb equi be :  " << nbequibe << endl;
	  for(int ieq=0,keq=0;keq<nbequibe;++keq)
	    {
		int p1[nbev],p2[nbev];
		int v1[nbev]={},v2[nbev]={};
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
                      if(keys[i0]< keys[i1] ) // not equal ... Add nov. 2014
		        equi.add(keys[i0],keys[i1]);
		  }

	    }

	  //  to find the final equivalent  key ...
	  //  in change of chaibe of equi key
	  for (int it=0,change=1;change;it++)
	    {
	    change=0;
	    ffassert(it++<100);
	    for (typename HashTable<Key,Key>::iterator qe,pe=equi.begin() ; pe != equi.end(); ++pe)
	      {
                  if( verbosity>9999) cout << pe->k << " " << pe->v << endl;
		  ffassert( pe->k < pe->v);
		  qe=equi.find(pe->v);
		  if(qe)
                  {
                      if( verbosity>9999) cout << pe->k << " " << pe->v << " <=> " << qe->k <<endl;

                     change++;
		     ffassert( qe->k < qe->v);
		      pe->v = qe->v;
		  }


	      }
	     if(verbosity>5)
	      cout << "     -- BuildDF: iteration in final equivalent " << it << " nb change " << change << endl;
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
	      {
	      if (nkf==1)
		keysdim[m]=2,keys[m++]=Key(k+of,tinfty);
	      else
		for(int ii,i=0;i<nkf;++i)
		  keysdim[m]=2,keys[m++]=Key(k+of,ElementAdj(k,ii=i) +of);
	      }
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
       if(verbosity)
	 {   cout << "  -- Build Nodes/DF on mesh :   n.v. " << nv <<  ", n. elmt. " << nt << ", n b. elmt. " <<nbe << endl;
	     cout << "     nb of Nodes " << nbNodes << "    nb of DoF   " << NbOfDF <<"  DFon="
	          << ndfon[0]<<ndfon[1]<<ndfon[2]<<ndfon[3]<<endl ;
	 }
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
    mes=0.;
    mesb=0.;

    for (int i=0;i<nt;i++)
	mes += this->elements[i].mesure();

    for (int i=0;i<nbe;i++)
	mesb += this->be(i).mesure();

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
  if(verbosity>3)
    cout << "  -- GMesh" << V::d << " , n V: " << nv << " , n Elm: " << nt << " , n B Elm: " << nbe << "mes " << mes << " " << mesb
	 << " , bb: (" << Pmin << ") , (" << Pmax << ")\n";
}

template<typename T,typename B,typename V>
void GenericMesh<T,B,V>::Buildbnormalv()
{
    const int nkv= T::nv;

    if (bnormalv)
      {return;}
    int nb=0;
    for (int k=0;k<nt;k++)
	for (int i=0;i<nkv;i++)
	{
	    int ii(i),kk;
	    kk=ElementAdj(k,ii);
	    if (kk<0 || kk==k) nb++;
	}
    if(verbosity>4)
	cout << " number of real boundary element " << nb << endl;
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
    assert(n - bnormalv <= nb );
}

static const char * GenericMesh_magicmesh="GenericMesh v0";
template<typename T,typename B,typename V>
Serialize GenericMesh<T,B,V>::serialize() const
{
    const int nve = T::nv;
    const int nvbe = B::nv;
    const int d = Rd::d;
    long long  l=0;
    l += sizeof(long long);
    l += 6*sizeof(int);
    l += nt*(sizeof(int)*(nve + 1));
    l += nv*( sizeof(int) + sizeof(double)*d);
    l += nbe*(sizeof(int)*(nvbe+1));
    if(verbosity>3)
       cout << "Serialize gmesh " << l << " " << nve << " " << nvbe << endl;
    Serialize  serialized(l,GenericMesh_magicmesh);
    // cout << l << magicmesh << endl;
    size_t pp=0;
    serialized.put(pp, l);
    serialized.put( pp,d);
    serialized.put( pp,nve);
    serialized.put( pp,nvbe);
    serialized.put( pp,nt);
    serialized.put( pp,nv);
    serialized.put( pp,nbe);
    if (verbosity>9)
	cout << " GenericMesh Serialized : " << l << " "  << nt << " " << nv << " " << nbe << endl;
    for (int i=0;i<nv;i++)
      {
	  for(int j=0;j<d;++j)
	  serialized.put(pp,vertices[i][j]);
	  serialized.put(pp,vertices[i].lab);
      }
    for (int i=0;i<nt;i++)
      {

	  const Element & K(elements[i]);
	  for(int j=0;j<nve;++j)
	      serialized.put(pp,(int) operator()(K[j]));
	  serialized.put(pp, K.lab);
      }
    for (int i=0;i<nbe;i++)
      {
	  const BorderElement & K(borderelements[i]);
	  for(int j=0;j<nvbe;++j)
	      serialized.put(pp,(int) operator()(K[j]));
	  serialized.put(pp, K.lab);
      }
    assert(pp==serialized.size());
    return serialized;
}

/*    GenericMesh()
    : nt(0),nv(0),nbe(0),  mes(0.),mesb(0.) ,
    vertices(0),elements(0),borderelements(0),bnormalv(0),
    TheAdjacencesLink(0),BoundaryElementHeadLink(0),
    ElementConteningVertex(0), gtree(0)
    {}
 */
    template<typename T,typename B,typename V>
    GenericMesh<T,B,V>::GenericMesh(const  Serialize &serialized)
    : nt(0),nv(0),nbe(0),  mes(0.),mesb(0.) ,
    vertices(0),elements(0),borderelements(0),bnormalv(0),
    TheAdjacencesLink(0),BoundaryElementHeadLink(0),
    ElementConteningVertex(0), gtree(0)
    {
	const int nve = T::nv;
	const int nvbe = B::nv;
	const int d = Rd::d;
	int dd,nnve,nnvbe,nnt,nnv,nnbe;
	long long  l=0;
	size_t pp=0;
	serialized.get(pp, l);
	serialized.get( pp,dd);
	serialized.get( pp,nnve);
	serialized.get( pp,nnvbe);
	serialized.get( pp,nnt);
	serialized.get( pp,nnv);
	serialized.get( pp,nnbe);
	ffassert(d==dd && nve == nnve && nvbe == nnvbe);
	set(nnv,nnt,nnbe);
	for (int i=0;i<nv;i++)
	  {
	      for(int j=0;j<d;++j)
		  serialized.get(pp,vertices[i][j]);
	      serialized.get(pp,vertices[i].lab);
	  }
	mes=0.;
	for (int i=0;i<nt;i++)
	  {
	      int ii[nve],lab;
	      for(int j=0;j<nve;++j)
		  serialized.get(pp,ii[j]);
	      serialized.get(pp,lab);
	     mes += elements[i].set(vertices,ii,lab).mesure();

	  }
	mesb=0;
	for (int i=0;i<nbe;i++)
	  {
	      int ii[nvbe],lab;
	      for(int j=0;j<nvbe;++j)
		  serialized.get(pp,ii[j]);
	      serialized.get(pp, lab);
	      mesb += borderelements[i].set(vertices,ii,lab).mesure();
	  }
       assert(pp==serialized.size());
    }

    
    template<typename TypeGenericElement, typename TypeVert >
    void SameElement( TypeVert *vertice, TypeGenericElement *list, int nelt, int *(&ind_nv_t), int *old2new, int &new_nelt, bool removeduplicate) {

        new_nelt=0;
        static const int nk=TypeGenericElement::nv;
        int iv[nk];
        
        HashTable<SortArray<int,nk>,int> h(nk*nelt,nelt);
        
        int *originmulti=new int[nelt];
        int *allmulti=new int[nelt];
        
        int originTypeGenericElement=0;
        int multiTypeGenericElement=0;
        int indice[nelt];
        
        for (int i=0;i<nelt;i++) {
            allmulti[i]=-1;
            originmulti[i]=-1;
            indice[i]=-1;
        }
        
        for (int i=0;i<nelt;i++) {
            
            const TypeGenericElement &K(list[i]);
            for (int j=0;j<nk;j++)
            iv[j] = old2new[ &(K[j]) - vertice ];    // vertice of element in new numbering
            int sens;
            SortArray<int,nk> a(iv,&sens);
            typename HashTable<SortArray<int,nk>,int>::iterator p= h.find(a);
            // check multiple element
            // 1/ keep the original
            if (!commonValue(a)) {   // if iv[0] != iv[1] (!= iv[2] != iv[3])
                if(!p) {
                    h.add(a,new_nelt);
                    indice[new_nelt]=i;
                    new_nelt++;
                }
                // or 2/ rm all multiples or keep the orinal
                else if(removeduplicate) {
                    originmulti[i]=p->v;  // the double elt the current
                    multiTypeGenericElement++;
               
                    if(originmulti[p->v]==-1) {   // the origin elt
                        originmulti[p->v]=p->v;
                        originTypeGenericElement++;
                    }
                }
            }
        }
        // rebuild the index list if remove all multiples
        if(removeduplicate) {
         int cmp=0;
            for(int i=0;i<nelt;i++) {
                if(originmulti[i]==-1) {
                    ind_nv_t[cmp]=i;
                    cmp++;
                }
            }
           assert((nelt-originTypeGenericElement-multiTypeGenericElement)==cmp);
            new_nelt=cmp;
            if (verbosity>2)
                cout << "no duplicate elements: "<< cmp << " and remove " << multiTypeGenericElement << " multiples elements, corresponding to " << originTypeGenericElement << " original elements " << endl;
            
        }
        else {
            for(int i=0;i<nelt;i++)
               ind_nv_t[i]=indice[i];
            if (verbosity>2)
                cout << " Warning, the mesh could contain multiple same elements, keep a single copy in mesh...option removemulti in the operator mesh is avalaible" << endl;
        }
 
        delete [] originmulti;
        delete [] allmulti;

    }
    
}
#endif
