// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY:  Bamg: Bidimensional Anisotrope Mesh Generator
// RELEASE: 0 
// USAGE  : You may copy freely these files and use it for    
//          teaching or research. These or part of these may   
//          not be sold or used for a commercial purpose with- 
//          out our consent : fax (33) 1 39 63 55 14       
//
// AUTHOR:   F. Hecht,    
// ORG    :  INRIA
// E-MAIL :   Frederic.Hecht@Inria.fr   
//
// ORIG-DATE:     Dec 97

namespace Fem2D {
const int MaxDeep = 30;
typedef  long  IntQuad;
const IntQuad MaxISize = ( 1L << MaxDeep);


class Mesh;
//typename Vertex;

class FQuadTree {
 public:
  class I2 { public:
     static bool INTER_SEG1d(long a,long b,long x,long y) { return (((y) > (a)) && ((x) <(b)));}
     long x,y;
     I2() {}
     I2(long i,long j): x(i),y(j) {}
     I2(const I2 &pp,long k,long l): x(pp.x+(( k&1) ? l : 0)),y(pp.y+(( k&2) ? l : 0)) {}
     void Add(long k,long l) { x+= (( k&1) ? l : 0) ; y+= (( k&2) ? l : 0);}
     I2(const I2 &A,const I2 &B) : x(B.x-A.x),y(B.y-A.y) {}

     long Case(long l) const  { return ( ( y & l) ? (( x & l) ? 3 : 2 ) :( ( x & l)? 1 : 0 ));}
     long norm() const { return Max(abs(x),abs(y));}
     bool less(I2 h) const  { return abs(x) <h.x && abs(y) <h.y;}
     bool interseg(I2 pp,long hb,long h) const { 
        return INTER_SEG1d(x,x+hb,pp.x-h,pp.x+h) && INTER_SEG1d(y,y+hb,pp.y-h,pp.y+h);}
     bool interseg(const I2 &pp,long hb,const I2 &h) const { 
        return INTER_SEG1d(x,x+hb,pp.x-h.x,pp.x+h.x) && INTER_SEG1d(y,y+hb,pp.y-h.y,pp.y+h.y);}
  };

  class QuadTreeBox { 
  public:

    long n; // if n < 4 => Vertex else =>  QuadTreeBox;
    union {
      QuadTreeBox *b[4];
      Vertex * v[4];
    };
    

  }; // end class QuadTreeBox  /////////////////

  class StorageQuadTreeBox {
  public:
    QuadTreeBox *b,*bc,*be;
    long len;
    StorageQuadTreeBox *n; // next StorageQuadTreeBox
    StorageQuadTreeBox(long ,StorageQuadTreeBox * =0);
    ~StorageQuadTreeBox()
    { //cout <<  "~StorageQuadTreeBox " << this << " n " << n << " b " << b << endl;
      if(n) delete n;
      delete [] b;
    }
    long  SizeOf() const {
      return len*sizeof(QuadTreeBox)+sizeof(StorageQuadTreeBox)+ (n?n->SizeOf():0);
    }
  }; // end class  StorageQuadTreeBox 
  
  StorageQuadTreeBox * sb;
  
  
  long  lenStorageQuadTreeBox;

public:
  QuadTreeBox * root;
  Mesh *th;
  
  R2 cMin,cMax; //  box of QuadTree
  R coef; //	
  long XtoI(R x) { return  (long) ((Max(Min(x,cMax.x),cMin.x)-cMin.x)*coef);}
  long YtoJ(R y) { return  (long) ((Max(Min(y,cMax.y),cMin.y)-cMin.y)*coef);}
  R ItoX(long i){ return double(i)*coef+cMin.x ;}
  R ItoY(long j){ return double(j)*coef+cMin.y ;}
  I2  R2ToI2(const R2 &P) { return I2(XtoI(P.x),YtoJ(P.y));}
  I2  R2ToI2(const R2 *P) { return I2(XtoI(P->x),YtoJ(P->y));}
  long NbQuadTreeBoxSearch,NbVerticesSearch;
  long NbQuadTreeBox,NbVertices;
  
  Vertex * NearestVertex(const R2 & P) { return NearestVertex(XtoI(P.x),YtoJ(P.y));}
  Vertex * NearestVertexWithNormal(const R2 & P);
  // { return NearestVertexWithNormal(XtoI(P.x),YtoJ(P.y));}
  
  Vertex * NearestVertex(long i,long j);
  
//  Vertex *  NearestVertexWithNormal(long i,long j); // new version
  Vertex *  ToClose(const R2 & ,R ,long,long);
  Vertex *  ToClose(const R2 & P,R delta){
      long hx = (long) (coef*delta);
      return ToClose(P,delta,hx,hx);}
  long SizeOf() const {return sizeof(FQuadTree)+sb->SizeOf();}

#ifdef DRAWING
  void  IMoveTo(long i,long j)
{  rmoveto(float(i)/coef+cMin.x ,float(j)/coef+cMin.y  );}

void  ILineTo(long i,long j)
 {rlineto(float(i)/coef+cMin.x ,float(j)/coef+cMin.y  );}
    void Draw();
  void PlotQuad(I2 pp,long hb);
  void PlotX(I2 pp,long hb);

#endif

  void  Add( Vertex & w);

  QuadTreeBox* NewQuadTreeBox()
  {
    ///cout << "NewQuadTreeBox " << sb << " " << sb->bc << " " 
    //<< sb->be << " " <<lenStorageQuadTreeBox <<endl;
    if(! (sb->bc<sb->be)) 
	sb=new StorageQuadTreeBox(lenStorageQuadTreeBox,sb);

    throwassert(sb && (sb->bc->n == 0));
    NbQuadTreeBox++;
    return sb->bc++;
  }
  ~FQuadTree();
  FQuadTree(Mesh * t,R2 Pmin,R2 Pmax,long nbv=-1);
  FQuadTree();
  friend ostream& operator <<(ostream& f, const  FQuadTree & qt);


	
};
//#undef IJ


} // name space 

