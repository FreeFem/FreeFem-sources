
//typedef double R;
namespace  Fem2D {
class MeshPointBase { public:
 R3 P;
  R3 PHat;
  const Mesh * Th;  
  const Triangle * T;
  long region, t,v,f,e,gsens; // triangle,vertex, face or edge
  long  label;  
  R3 N; //  if on boundary 
  bool outside;
  void set(const R2 &P2,const R2 & P_Hat,const  baseFElement & K,int ll,const R2 &NN,int iedge)
   { 
     P.x=P2.x;
     P.y=P2.y;
     P.z=0;
     PHat.x=P_Hat.x;
     PHat.y=P_Hat.y;
     PHat.z;
     T=&K.T;
     Th=&K.Vh.Th; 
     region = T->lab;
     label = ll;
     v=f=-1; 
     e=iedge;
     t=(*Th)(T); 
     throwassert( Abs( (NN,NN) -1.0) < 1e-5 );
     N.x=NN.x;   
     N.y=NN.y;   
     N.z=0;   
   }
  void set(const Mesh & aTh,const R2 &P2,const R2 & P_Hat,const  Triangle & aK,int ll,const R2 &NN,int iedge)
   { 
     P.x=P2.x;
     P.y=P2.y;
     P.z=0;
     PHat.x=P_Hat.x;
     PHat.y=P_Hat.y;
     PHat.z;
     T=&aK;
     Th=&aTh; 
     region = T->lab;
     label = ll;
     v=f=-1;
     t=(*Th)(T); 
     e=iedge; 
     throwassert( Abs( (NN,NN) -1.0) < 1e-5 );
     N.x=NN.x;   
     N.y=NN.y;   
     N.z=0;   
   }
   
  void set(const R2 &P2,const R2 & P_Hat,const  baseFElement & K,int ll)
   { 
     P.x=P2.x;
     P.y=P2.y;
     P.z=0;
     PHat.x=P_Hat.x;
     PHat.y=P_Hat.y;
     PHat.z;
     T=&K.T;
     Th=&K.Vh.Th; 
     region = T->lab;
     label = ll;
     t=(*Th)(T);
     v=f=e=-1;  
     N.x=0;   
     N.y=0;   
     N.z=0;   
   }
   
     void set(const R2 &P2,const R2 & P_Hat,const  baseFElement & K)
   { 
     P.x=P2.x;
     P.y=P2.y;
     P.z=0;
     PHat.x=P_Hat.x;
     PHat.y=P_Hat.y;
     PHat.z;
     T=&K.T;
     Th=&K.Vh.Th; 
     region = T->lab;
     v=f=e=-1;  
     int l0(0),l1(0),l2(0),ke(-1),kk(0);
     if ( P_Hat.x<1.e-6) l1=1,ke=1,kk++;
     if ( P_Hat.y<1.e-6) l2=1,ke=2,kk++;
     if ( P_Hat.y+P_Hat.x>0.999999) l0=1,ke=0,kk++;     
     if (kk==0) label=0;
     else if (kk==2)
      { 
      if (!l0) label=(*T)[0].lab;
      else if (!l1) label=(*T)[1].lab;
      else if (!l2) label=(*T)[2].lab;
      } 
     else  {
       int i1,i2;
       Th->VerticesNumberOfEdge(K.T,ke,i1,i2);
       const BoundaryEdge * be=Th->TheBoundaryEdge(i1,i2);
       label= be ? be->lab : 0;
//       cout << "lab =" <<  label << " " << ke << " " <<  kk << " " << P_Hat 
//            << ": " <<  K.number << " , " << (R2) P << endl;
      }
   
     t=(*Th)(T);
     N.x=0;   
     N.y=0;   
     N.z=0;   
   }

  void set(const  Mesh &aTh, const R2 &P2,const R2 & P_Hat,const  Triangle & aK,const int ll,bool coutside=false)
   { 
     P.x=P2.x;
     P.y=P2.y;
     P.z=0;
     PHat.x=P_Hat.x;
     PHat.y=P_Hat.y;
     PHat.z;
     T=&aK;
     Th=&aTh; 
     region = T->lab;
     label = ll;
     t=v=f=e=-1;  
     N.x=0;   
     N.y=0;   
     N.z=0;   
     outside=coutside;
   }
   
  void setP(Mesh * pTh,int tt,int ss)
   { 
     T=&(*pTh)[tt];
     Vertex & V=(*T)[ss];
     (R2 &) P= V ;
     P.z=0;
     (R2 &) PHat = TriangleHat[ss];
     PHat.z=0;
     Th=pTh; 
     region = T->lab;
     label = V.lab;
     t=v=f=e=0;
     v=ss;
   }
   
  void change(const R2 & PH,Triangle & tt,int ll)
   { 
     T= &tt;
     (R2 &) PHat = PH;
     (R2 &) P = (*T)(PH);
     region = T->lab;
     label = ll;
     t=v=f=e=0;
     
   }
   void unset() 
   {
     P.x=-1e30;
     P.y=-1e30;
     P.z=-1e30;
     T=0;
     Th=0;
      label =0;
     
   }
   bool isUnset() const { return P.x == -1e30;} // BofBof   
   void set(R x=0.0,R y=0.0,R z=0.0) 
   {
     P.x=x;
     P.y=y;
     P.z=z;
     T=0;
     Th=0;
      label =0;
      t=f=e=v=-1; 
     
     
   }  
};
class MeshPoint : public MeshPointBase { public:
  MeshPointBase other;
  void unset() {  MeshPointBase::unset(); other.unset();}
  void set(const R2 &P2,const R2 & P_Hat,const  baseFElement & K,int ll,const R2 &NN,int iedge) {
     MeshPointBase::set(P2,P_Hat,K,ll,NN,iedge);
     other.unset();}   
  void set(const Mesh & aTh,const R2 &P2,const R2 & P_Hat,const  Triangle & aK,int ll,const R2 &NN,int iedge) {
    MeshPointBase::set(aTh,P2,P_Hat,aK,ll,NN,iedge);
    other.unset();}
  void set(const R2 &P2,const R2 & P_Hat,const  baseFElement & K) {
     MeshPointBase::set(P2,P_Hat,K);
     other.unset();
    }
  void set(const  Mesh &aTh, const R2 &P2,const R2 & P_Hat,const  Triangle & aK,
           const int ll,bool coutside=false) {
      MeshPointBase::set(aTh,P2,P_Hat,aK,ll,coutside);
      other.unset();
           }
    
  void set(R x=0.0,R y=0.0,R z=0.0) {  
     MeshPointBase::set(x,y,z);
     other.unset();} 
  void change(const R2 & PH,Triangle & tt,int ll) {
     MeshPointBase::change(PH,tt,ll);
     other.unset(); }
   void setP(Mesh * pTh,int tt,int ss) { 
      MeshPointBase::setP(pTh,tt,ss); 
      other.unset(); } 
   bool operator==(const MeshPoint & mp) const {
      return T == mp.T &&  P.x == mp.P.x && P.y == mp.P.y 
          && P.z == mp.P.z ;}
  bool  SetAdj() {
     assert(Th && T && t >=0 && e>=0);
     int ieo=e,to=t,ie=e;
     
     t=Th->TriangleAdj(t,ie);
     e=ie;
     if ( t == to && t >= 0 ) return false;
     int io0=VerticesOfTriangularEdge[ieo][0];
     int io1=VerticesOfTriangularEdge[ieo][1];
     int i0=VerticesOfTriangularEdge[e][0];
     int i1=VerticesOfTriangularEdge[e][1];
     T= &(*Th)[t];
     R l[3];
     l[1]=PHat.x;
     l[2]=PHat.y;
     l[0]=1-PHat.x-PHat.y;
     R le=l[io0]; 
     l[i1]=le;
     l[i0]=1-le;
     l[3-i1-i0]=0;
     PHat.x=l[1];
     PHat.y=l[2];
     gsens = -gsens;     
     return true;         
   }
  };

ostream & operator << ( ostream &,const  MeshPoint & )  ;    
inline static MeshPoint* MeshPointStack(Stack s) {void * p= static_cast<void **>(s)[0];throwassert(p); return static_cast<MeshPoint*>( p);}
inline static void MeshPointStack(Stack s,MeshPoint* mp) {*static_cast<MeshPoint**>(static_cast<void *>(s)) = mp;}

}
