#ifndef FEM2_H_
#define FEM2_H_
#include <string> 
#include "RefCounter.hpp"
#include "Serialize.hpp"
// some usefull function 


template<class K> class KN_;

const  double Pi =  3.14159265358979323846264338328;

template<class T> inline T Min (const T &a,const T &b){return a < b ? a : b;}
template<class T> inline T Max (const T &a,const T & b){return a > b ? a : b;}
template<class T> inline T Abs (const T &a){return a <0 ? -a : a;}

template<class T> inline void Exchange (T& a,T& b) {T c=a;a=b;b=c;}
template<class T> inline T Max (const T &a,const T & b,const T & c){return Max(Max(a,b),c);}
template<class T> inline T Min (const T &a,const T & b,const T & c){return Min(Min(a,b),c);}


// definition R
typedef double R;
namespace Fem2D {

// The class R2
class R2 {
  friend ostream& operator <<(ostream& f, const R2 & P )
  { f << P.x << ' ' << P.y   ; return f; }
  friend istream& operator >>(istream& f,  R2 & P)
  { f >>  P.x >>  P.y  ; return f; }

public:  
  typedef R Rcross; 
  R x,y;
  R2 () :x(0),y(0) {};
  R2 (R a,R b):x(a),y(b)  {}
  R2 (R2 A,R2 B) :x(B.x-A.x),y(B.y-A.y)  {}
  R2   operator+(R2 P)const   {return R2(x+P.x,y+P.y);}
  R2   operator+=(R2 P)  {x += P.x;y += P.y;return *this;}
  R2   operator-(R2 P)const   {return R2(x-P.x,y-P.y);}
  R2   operator-=(R2 P) {x -= P.x;y -= P.y;return *this;}
  R2   operator-()const  {return R2(-x,-y);}
  R2   operator+()const  {return *this;}
  R   operator,(R2 P)const  {return  x*P.x+y*P.y;} // produit scalaire
  Rcross   operator^(R2 P)const {return  x*P.y-y*P.x;} // produit mixte
  R2   operator*(R c)const {return R2(x*c,y*c);}
  R2   operator/(R c)const {return R2(x/c,y/c);}
  R2    operator*=(R c)  {x *= c;y *= c;return *this;}
  R2    operator/=(R c)  {x /= c;y /= c;return *this;}
  R  &  operator[](int i){ return (&x)[i];}
  R2   perp() {return R2(-y,x);} // the perpendicular
  friend R2 operator*(R c,R2 P) {return P*c;}
  friend R2 operator/(R c,R2 P) {return P/c;}
};

// The class R3
class R3: public R2 {
  friend ostream& operator <<(ostream& f, const R3 & P )
  { f << P.x << ' ' << P.y << ' ' << P.z   ; return f; }
  friend istream& operator >>(istream& f,  R3 & P)
  { f >>  P.x >>  P.y >>  P.z  ; return f; }

public:  
  typedef R3 Rcross; 
  R z;
 
  R3 () :z(0) {};
  R3 (R a,R b,R c):R2(a,b),z(c)  {}
  R3 (R2 P2):R2(P2),z(0)  {}
  R3 (const R3 & A,const R3 & B) :R2(B.x-A.x,B.y-A.y),z(B.z-A.z)  {}
  R3 & operator=(R2 P2) {x=P2.x;y=P2.y;z=0;return *this;}
  R3   operator+(R3 P)const   {return R3(x+P.x,y+P.y,z+P.z);}
  R3 & operator+=(R3 P)  {x += P.x;y += P.y;z += P.z;return *this;}
  R3   operator-(R3 P)const   {return R3(x-P.x,y-P.y,z-P.z);}
  R3 & operator-=(R3 P) {x -= P.x;y -= P.y;z -= P.z;return *this;}
  R3   operator-()const  {return R3(-x,-y,-z);}
  R3   operator+()const  {return *this;}
  R    operator,(R3 P)const  {return  x*P.x+y*P.y+z*P.z;} // produit scalaire
  Rcross  operator^(R3 P)const {return R3(y*P.z-z*P.y ,P.x*z-x*P.z, x*P.y-y*P.x);} // produit mixte
  R3   operator*(R c)const {return R3(x*c,y*c,z*c);}
  R3   operator/(R c)const {return R3(x/c,y/c,z/c);}
  R  &  operator[](int i){ return (&x)[i];}
  friend R3 operator*(R c,R3 P) {return P*c;}
  friend R3 operator/(R c,R3 P) {return P/c;}
};

inline void MoveTo(R2 P) { rmoveto((float) P.x,(float)P.y);}
inline void LineTo(R2 P) { rlineto((float)P.x,(float)P.y);}

inline R Area2(const R2 A,const R2 B,const R2 C){return (B-A)^(C-A);} 
inline R Norme2_2(const R2 & A){ return (A,A);}
inline R Norme2(const R2 & A){ return sqrt((A,A));}
inline R Norme2_2(const R3 & A){ return (A,A);}
inline R Norme2(const R3 & A){ return sqrt((A,A));}
inline R Norme_infty(const R2 & A){return Max(Abs(A.x),Abs(A.y));}
inline R Norme_infty(const R3 & A){return Max(Abs(A.x),Abs(A.y),Abs(A.z));}
inline R Theta(R2 P){ return atan2(P.y,P.x);}

// def de numerotation dans un triangles direct sens (trigo)
// the edge is oposite of the vertex
//  [3] is a edge
 const short VerticesOfTriangularEdge[3][2] = {{1,2},{2,0},{0,1}};
 //  [3] is a vertices 
 const short EdgesVertexTriangle[3][2] = {{1,2},{2,0},{0,1}};
 const short OppositeVertex[3] = {0,1,2};
 const short OppositeEdge[3] =  {0,1,2};
 const short NextEdge[3] = {1,2,0};
 const short PreviousEdge[3] = {2,0,1};
 const short NextVertex[3] = {1,2,0};
 const short PreviousVertex[3] = {2,0,1};
 //  array to same is some onwhat is on edge for boundary condition
 // is on a edge i  onwhat is in [0,7[ 
 // 2  on edge 1 and   
 const int onWhatIsEdge[3][7] = { {0,1,3, 2,0,0, 0}, // edge 0 
                                   {3,0,1, 0,2,0, 0},
                                   {1,3,0, 0,0,2, 0}};
                                  
const   R2 TriangleHat[3]= { R2(0,0),R2(1,0),R2(0,1) } ; 
                                  
// --------
class Label {  // reference number for the physics
public: 
  int lab;
  Label(int r=0):lab(r){}
  int operator!() const{return !lab;} 
  bool operator<(const Label & r) const {return lab < r.lab;} 
  bool operator==(const Label & r) const {return lab == r.lab;} 
  bool operator>(const Label & r) const { return lab > r.lab;} 
};
inline ostream& operator <<(ostream& f,const Label & r  )
  { f <<  r.lab ; return f; }
inline istream& operator >>(istream& f, Label & r  )
  { f >>  r.lab ; return f; }
  

class Vertex : public R2,public Label {
 friend class Mesh;
   R2 *normal; // pointeur sur la normal exterieur pour filtre les
   // point de depart 
public:
  Vertex() : R2(),Label(),normal(0){};
  Vertex(R2 P,int r=0): R2(P),Label(r),normal(0){}
  bool ninside(const R2 & P) const { return normal? (R2(*this,P),*normal)<=0: true;}
  void SetNormal(R2 *&n,const R2 & N)
     { if (normal) { 
         R2 NN=*normal+N; 
         *normal= NN/Norme2(NN); }
       else *(normal=n++)=N;}
  R2 Ne() const {return normal ? *normal: R2(0,0);}
//  void operator=(const Vertex & v) { x=v.x;y=v.y;lab=v.lab;normal=0;}
};
inline ostream& operator <<(ostream& f, const Vertex & v )
  { f << (R2) v << ' ' << (Label) v   ; return f; }
inline istream& operator >> (istream& f,  Vertex & v )
  { f >> (R2&) v >> (Label&) v ; return f; }



class Triangle: public Label {
  Vertex *vertices[3]; // an array of 3 pointer to vertex
public:
  R area;
  Triangle(){};              // constructor empty for array
  Vertex & operator[](int i) const// to see triangle as a array of vertex
       {return *vertices[i];} 
       
  Vertex *& operator()(int i) // to see triangle as a array of vertex
       {return vertices[i];} 
  
  Triangle(Vertex * v0,int i0,int i1,int i2,int r,R a=0.0): Label(r) 
     { R2 A=*(vertices[0]=v0+i0);
       R2 B=*(vertices[1]=v0+i1);
       R2 C=*(vertices[2]=v0+i2); 
       area = a ==0 ? (( B-A)^(C-A))*0.5 : a;
       throwassert(area>0);}

  void set(Vertex * v0,int i0,int i1,int i2,int r,R a=0.0)
     { lab=r; 
       R2 A=*(vertices[0]=v0+i0);
       R2 B=*(vertices[1]=v0+i1);
       R2 C=*(vertices[2]=v0+i2); 
       area = a ==0 ? (( B-A)^(C-A))*0.5 : a;
       throwassert(area>0);}
       
  R2 Edge(int i) const // opposite edge vertex i
     {return (R2) *vertices[(i+2)%3]-(R2) *vertices[(i+1)%3];}
     
  Vertex & Edge(int j,int i) const // Vertex j of edge i
     {throwassert(j==0 || j==1 );return  *vertices[(i+j+1)%3];}
  
  R2 n(int i) const //  unit exterior normal
     {R2 E=Edge(i);return R2(E.y,-E.x)/Norme2(E);} 
  
  R2 H(int i)  const  // heigth ($\nabla \lambda_i$ 
     {R2 E=Edge(i);return R2(-E.y,E.x)/(2*area);} 
     
  R lenEdge(int i) const {R2 E=Edge(i);return sqrt((E,E));}
  R h() const { return Max(lenEdge(0),lenEdge(1),lenEdge(2));}
  
  void Renum(Vertex   *v0, long * r)  { 
    for (int i=0;i<3;i++) 
      vertices[i]=v0+r[vertices[i]-v0];}
      
  Vertex & VerticeOfEdge(int i,int j) const  // vertex j of edge i 
    {return  *vertices[(i+1+j)%3];}  // vertex j of edge i 

    R EdgeOrientation(int i) const { // return +1 or -1 
     R Orient[2]={-1.,1.};
    return  Orient[vertices[ (i+1)%3] < vertices[ (i+2)%3] ] ;}
    
  bool intersect(R2 P,R2 Q) const 
   { 
     const R2 &A(*vertices[0]);
     const R2 &B(*vertices[1]);
     const R2 &C(*vertices[2]);
     R2 mn(Min(A.x,B.x,C.x),Min(A.y,B.y,C.y)),
        mx(Max(A.x,B.x,C.x),Max(A.y,B.y,C.y));
     assert(P.x < Q.x && P.y < Q.y ); 
     return (mx.x >= P.x) && (mn.x <= Q.x) &&  (mx.y >= P.y) && (mn.y <= Q.y)  ;

   }
     
//  const Vertex & VerticeOfEdge(int i,int j) const {return  *vertices[(i+1+j)%3];}      
  void Draw(double shink=1) const;
  void Fill(int color) const;
  void Draw(int edge,double shink=1) const;
  void SetVertex(int j,Vertex *v){vertices[j]=v;}
  R2 operator()  (const R2 & P) const{ // local to Global in triangle 
      return    (const R2 &) *vertices[0] * (1-P.x-P.y) 
             +  (const R2 &) *vertices[1] * (P.x) 
             +  (const R2 &) *vertices[2] * (P.y) ;}
private:
  Triangle(const Triangle &);  //  no copy of triangle
  void operator=(const Triangle &);             

};

class BoundaryEdge: public Label {
public:
  Vertex *vertices[2];
  BoundaryEdge(Vertex * v0,int i0,int i1,int r): Label(r) 
  { vertices[0]=v0+i0; vertices[1]=v0+i1; }
  void set(Vertex * v0,int i0,int i1,int r)
  { lab=r,vertices[0]=v0+i0; vertices[1]=v0+i1; }
  bool in(const Vertex * pv) const {return pv == vertices[0] || pv == vertices[1];}
  BoundaryEdge(){}; // constructor empty for array 
  void Draw() const;
  Vertex & operator[](int i) const {return *vertices[i];}
  R length() const { return Norme2(R2(*vertices[0],*vertices[1]));}
  void Renum(Vertex   *v0, long * r) { 
    for (int i=0;i<2;i++) 
      vertices[i]=v0+r[vertices[i]-v0];}
};

typedef BoundaryEdge Edge;
typedef Triangle Tet;  // just to play

class Mortar { 
  public:
  friend class Mesh;
  friend class ConstructDataFElement;
   Mesh * Th;
   int nleft,nright;
   int *left,*right;
   Mortar(): Th(0),left(0),right(0),nleft(0),nright(0){}
   void Draw() const;
   public:
     int NbLeft() const{return nleft;} 
     int NbRight() const{return nright;} 
     int TLeft(int i) const { return left[i]/3;}
     //   
     int NbT() const {return nleft+nright;} 
     int T_e(int i,int & e)  const { // return the triangle number + the edge number 
        throwassert(i>=0 && i < nleft+nright);
        int k= (i<nleft ? left[i]: right[i-nleft]);
        e=k%3;
        return k/3;}
     int ELeft(int i) const { return left[i]%3;}
     int TRight(int i) const { return right[i]/3;}
     int ERight(int i) const { return right[i]%3;}
     
     //  warning  i is in [0, nleft]
     Vertex & VLeft(int i) const ;
     Vertex & VRight(int i) const;
     
     
};

 
class FQuadTree;

class Mesh: public RefCounter { public:
   static const char magicmesh[8]  ;
  int dim; 
  int nt,nv,neb,ne,ntet;
  R area;
  R volume;
  static int kthrough,kfind;
  FQuadTree *quadtree; 
  Vertex *vertices;
  Triangle *triangles;
  BoundaryEdge  *bedges;
  Edge  *edges;  // edge element 
  Tet * tet; // 
  
    
  int NbMortars,NbMortarsPaper;
  Mortar *mortars; //  list of mortar
  int    *datamortars;
  R2  * bnormalv; //  boundary vertex normal 
  //Triangle * adj;
  Triangle & operator[](int i) const {throwassert(i>=0 && i<nt);return triangles[i];}
 // const Triangle & operator[](int i) const {return triangles[i];}
  Vertex & operator()(int i) const {throwassert(i>=0 && i<nv);return vertices[i];}
  Mesh(const char * filename) {read(filename);} // read on a file
  Mesh(const string s) {read(s.c_str());}
  Mesh( const Serialize & ) ;
  Mesh(int nbv,R2 * P);
  
  Serialize serialize() const;
  Mesh(int nbv,int nbt,int nbeb,Vertex *v,Triangle *t,BoundaryEdge  *b);  
  Mesh(const Mesh & Thold,int *split,bool WithMortar=true,int label=1);
  ~Mesh();
  int number(const Triangle & t) const {return &t - triangles;}
  int number(const Triangle * t) const {return t  - triangles;}
  int number(const Vertex & v)   const {return &v - vertices;}
  int number(const Vertex * v)   const {return v  - vertices;}
  int operator()(const Triangle & t) const {return &t - triangles;}
  int operator()(const Triangle * t) const {return t  - triangles;}
  int operator()(const Vertex & v)   const {return &v - vertices;}
  int operator()(const Vertex * v)   const {return v  - vertices;}
  int operator()(int it,int j) const {return number(triangles[it][j]);}
        // Nu vertex j of triangle it
  void BoundingBox(R2 & Pmin,R2 &Pmax) const;
  void InitDraw() const ;
  void Draw(int init=2,bool fill=false) const;
  void DrawBoundary() const; 
  int Contening(const Vertex * v) const{ return TriangleConteningVertex[ v  - vertices];}
  int renum();
  int gibbsv (long* ptvoi,long* vois,long* lvois,long* w,long* v);
  int TriangleAdj(int it,int &j) const 
      {int i=TheAdjacencesLink[3*it+j];j=i%3;return i/3;}
  int nTonEdge(int it,int e) const { int k=3*it+e;return k==TheAdjacencesLink[k] ? 1 : 2;}
      
  void VerticesNumberOfEdge(const Triangle & T,int j,int & j0,int & j1) const 
      {j0 =  number(T[(j+1)%3]),j1=  number(T[(j+ 2)%3]);}
  bool SensOfEdge(const Triangle & T,int j) const 
      { return  number(T[(j+1)%3]) <number(T[(j+ 2)%3]);}
      
  int GetAllTriangleAdj(int it,int *tabk) const
   { //  get the tab of all adj triangle to a traingle (max 3)
     //  and return the size of the tab 
      int i=0;
      tabk[i]=TheAdjacencesLink[3*it+0]/3;
      if(tabk[i] >=0 && tabk[i]!=it) i++; 
      tabk[i]=TheAdjacencesLink[3*it+1]/3;
      if(tabk[i] >=0 && tabk[i]!=it) i++; 
      tabk[i]=TheAdjacencesLink[3*it+2]/3;
      if(tabk[i] >=0 && tabk[i]!=it) i++; 
      return i;
   }
      
  int BoundaryTriangle(int be,int & edgeInT) const {
     int i= BoundaryEdgeHeadLink[be]; edgeInT = i%3; 
     return i/3;}
     
  Triangle * Find(const R2 & P) const ;
  const Triangle * Find(R2 P, R2 & Phat,bool & outside,const Triangle * tstart=0) const  ;
  
  BoundaryEdge * TheBoundaryEdge(int i,int j)  const
   {  
    int p2;
    for (int p=BoundaryAdjacencesHead[i];p>=0;p=BoundaryAdjacencesLink[p])
     { 
     if ( bedges[p2=p/2].in(vertices+j) )   return bedges+p2;
     }
     
    return 0;}
 void destroy() {RefCounter::destroy();}
 void MakeQuadTree() ;
private:
  void read(const char * filename); 
  void BuildBoundaryAdjacences();
  void ConsAdjacence();
  void Buildbnormalv(); 

  // to construct the adj triangle               
  int *TheAdjacencesLink;
  int *BoundaryEdgeHeadLink;
  int *BoundaryAdjacencesHead;
  int *BoundaryAdjacencesLink; 
  int *TriangleConteningVertex;               
};


//  2 routines to compute the caracteristic
int WalkInTriangle(const Mesh & Th,int it, double *lambda,
                   const  KN_<R> & U,const  KN_<R> & V, R & dt);
int  WalkInTriangle(const Mesh & Th,int it, double *lambda,
                     R u, R v, R & dt); 
                   
int Walk(const Mesh & Th,int & it, R *l,
         const KN_<R>  & U,const KN_<R>  & V, R dt) ;
 
void DrawMark(R2 P,R k=0.02);

         
template<class T>
void  HeapSort(T *c,long n)
{
  long l,j,r,i;
  T crit;
  c--; // on decale de 1 pour que le tableau commence a 1
  if( n <= 1) return;
  l = n/2 + 1;
  r = n;
  while (1) { // label 2
    if(l <= 1 ) { // label 20
      crit = c[r];
      c[r--] = c[1];
    if ( r == 1 ) { c[1]=crit; return;}
    } else  crit = c[--l]; 
    j=l;
    while (1) {// label 4
      i=j;
      j=2*j;
      if  (j>r) {c[i]=crit;break;} // L8 -> G2
      if ((j<r) && (c[j] < c[j+1])) j++; // L5
      if (crit < c[j]) c[i]=c[j]; // L6+1 G4
      else {c[i]=crit;break;} //L8 -> G2
    }
  }
}

template<class T,class T1,class T2>
void  HeapSort(T *c,T1 *c1,T2 *c2,long n)
{
  long l,j,r,i;
  T crit;  T1 crit1;  T2 crit2;
  c--;c1--;c2--; // on decale de 1 pour que le tableau commence a 1
  if( n <= 1) return;
  l = n/2 + 1;
  r = n;
  while (1) { // label 2
    if(l <= 1 ) { // label 20
      crit = c[r];     crit1 = c1[r];  crit2 = c2[r];
      c2[r] = c2[1];   c1[r] = c1[1];  c[r--] = c[1];
    if ( r == 1 ) { c2[1] = crit2,   c1[1] = crit1; c[1]=crit; return;}
    } else  {crit = c[--l];crit1=c1[l];crit2=c2[l]; }
    j=l;
    while (1) {// label 4
      i=j;
      j=2*j;
      if  (j>r) {c[i]=crit;c1[i]=crit1;c2[i]=crit2;break;} // L8 -> G2
      if ((j<r) && (c[j] < c[j+1])) j++; // L5
      if (crit < c[j]) {c[i]=c[j];c1[i]=c1[j];c2[i]=c2[j];} // L6+1 G4
      else {c[i]=crit;c1[i]=crit1;c2[i]=crit2;break;} //L8 -> G2
    }
  }
}
          
 R2 SubTriangle(const int N,const int n,const int l);
 
      //  warning  i is in [0, nleft]
 inline     Vertex & Mortar::VLeft(int i) const 
  { throwassert(i>=0 && i <= nleft);
    return i< nleft ? (*Th)[TLeft( i  )][VerticesOfTriangularEdge[ELeft(i)][0]]
                    : (*Th)[TLeft( i-1)][VerticesOfTriangularEdge[ELeft(i-1)][1]];}
 inline Vertex & Mortar::VRight(int i) const 
  { throwassert(i>=0 && i <= nright);
    return i< nright ? (*Th)[TRight( i  )][VerticesOfTriangularEdge[ERight(i)][1]]
                     : (*Th)[TRight( i-1)][VerticesOfTriangularEdge[ERight(i-1)][0]];}

 void DrawCommentaire(const char *cm,float x=0.1,float y=0.97) ;
 
inline R2 minmax(const R2 & a,const R2 & b) 
  {return R2(Min(a.x,b.x),Max(a.y,b.y));}


}



#include "FQuadTree.hpp"

namespace Fem2D {

inline void Mesh::MakeQuadTree() 
  {
    if (!quadtree) {
        R2 Pn,Px;
        BoundingBox(Pn,Px);
        quadtree=new FQuadTree(this,Pn,Px,nv);
    } 
 }

}

#endif
