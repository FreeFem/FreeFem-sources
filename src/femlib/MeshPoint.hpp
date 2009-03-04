
//typedef double R;
namespace  Fem2D {
class MeshPointBase { public:
  R3 P;
  R3 PHat;
  union {
  const Mesh * Th;  
  const Mesh3 * Th3;  
  };
  union{
  const Triangle * T;
  const Tet * T3;
  };
  long region, t,v,f,e,gsens; // triangle,vertex, face or edge
  long  label;  
  R3 N; //  if on boundary 
  bool outside;
  int VF; 
  int d;
  void set(const R2 &P2,const R2 & P_Hat,const  baseFElement & K,int ll,const R2 &NN,int iedge)
   { 
     P.x=P2.x;
     P.y=P2.y;
     P.z=0;
     PHat.x=P_Hat.x;
     PHat.y=P_Hat.y;
     PHat.z=0.;
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
     VF=0;
     d=2;
   }
  void set(const Mesh & aTh,const R2 &P2,const R2 & P_Hat,const  Triangle & aK,int ll,const R2 &NN,int iedge,int VFF=0)
   { 
     P.x=P2.x;
     P.y=P2.y;
     P.z=0;
     PHat.x=P_Hat.x;
     PHat.y=P_Hat.y;
     PHat.z=0.;
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
     VF=VFF;
     d=2;
   }
   
  void set(const R2 &P2,const R2 & P_Hat,const  baseFElement & K,int ll)
   { 
     P.x=P2.x;
     P.y=P2.y;
     P.z=0;
     PHat.x=P_Hat.x;
     PHat.y=P_Hat.y;
     PHat.z=0.;
     T=&K.T;
     Th=&K.Vh.Th; 
     region = T->lab;
     label = ll;
     t=(*Th)(T);
     v=f=e=-1;  
     N.x=0;   
     N.y=0;   
     N.z=0;   
     VF=0;  
     d=2;
   }
   
     void set(const R2 &P2,const R2 & P_Hat,const  baseFElement & K)
   { 
     P.x=P2.x;
     P.y=P2.y;
     P.z=0;
     PHat.x=P_Hat.x;
     PHat.y=P_Hat.y;
     PHat.z=0.;
     T=&K.T;
     Th=&K.Vh.Th; 
     region = T->lab;
     v=f=e=-1;  
     N.x=0;   
     N.y=0;   
     N.z=0;   
     VF=0;  
     int ll[3],kk(0);
     if ( P_Hat.x<1.e-6) ll[kk++]=1;
     if ( P_Hat.y<1.e-6) ll[kk++]=2;
     if ( P_Hat.y+P_Hat.x>0.999999) ll[kk++]=0;    
     if (kk==0) label=0;
     else if (kk==2)
      { 
	 v = 3-ll[0]-ll[1];// 3 = 0+1+2 sommet oppose
        label=(*T)[v].lab;
      } 
     else  {
       e = ll[0]; 
       int i1,i2;
       Th->VerticesNumberOfEdge(K.T,e,i1,i2);
       const BoundaryEdge * be=Th->TheBoundaryEdge(i1,i2);
       label= be ? be->lab : 0;
      // R2 E(K.T.Edge(ke));
      // (R2 &) N = E.perp()/Norme2(E);
       
   
      //cout << "lab =" <<  label << " " << e << " " <<  kk << " " << P_Hat 
       //    << ": " <<  K.number << " , " << (R2) P << " " << N << endl;
      }
   
     t=(*Th)(T);
     d=2;
   }

  void set(const  Mesh &aTh, const R2 &P2,const R2 & P_Hat,const  Triangle & aK,const int ll,bool coutside=false)
   { 
     P.x=P2.x;
     P.y=P2.y;
     P.z=0;
     PHat.x=P_Hat.x;
     PHat.y=P_Hat.y;
     PHat.z=0.;
     T=&aK;
     Th=&aTh; 
     region = T->lab;
     label = ll;
     t=v=f=e=-1;  
     N.x=0;   
     N.y=0;   
     N.z=0;   
     outside=coutside;
     VF=0;  
     d=2;
   }
   
  void setP(const Mesh * pTh,int tt,int ss)
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
     VF=0;
     d=2;  
   }
   
  void change(const R2 & PH,const Triangle & tt,int ll)
   { 
     T= &tt;
     (R2 &) PHat = PH;
     (R2 &) P = (*T)(PH);
     region = T->lab;
     label = ll;
     t=v=f=e=0;
     VF=0;
     d=2;  
     
   }
    void change(const R3 & PH,const Tet & tt,int ll)
    { 
	T3= &tt;
	(R3 &) PHat = PH;
	(R3 &) P = (*T3)(PH);
	region = T3->lab;
	label = ll;
	t=v=f=e=0;
	VF=0;
	d=2;  
	
    }
    
   void unset() 
   {
     P.x=-1e30;
     P.y=-1e30;
     P.z=-1e30;
     T=0;
     Th=0;
      label =0;
     VF=0;  
      d=0;
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
      VF=0;  
      d=0;
   }



  // ------- 3D
  void set(const R3 &P2,const R3 & P_Hat,const  baseFElement3 & K,int ll,const R3 &NN,int iface)
   { 
     P=P2;
     PHat=P_Hat;
     T3=&K.T;
     Th3=&K.Vh.Th; 
     region = T3->lab;
     label = ll;
     e=v=-1; 
     f=iface;
     t=(*Th3)(T3); 
     assert( Abs( (NN,NN) -1.0) < 1e-5 );
     N=NN;   
     VF=0;
     d=3;
   }
  void set(const Mesh3 & aTh,const R3 &P2,const R3 & P_Hat,const  Tet & aK,int ll,const R3 &NN,int iface,int VFF=0)
   { 
     P=P2;
     PHat=P_Hat;
     T3=&aK;
     Th3=&aTh; 
     region = T3->lab;
     label = ll;
     v=f=-1;
     t=(*Th3)(T3); 
     v=e=-1;
     f=iface; 
     assert( Abs( (NN,NN) -1.0) < 1e-5 );
     N=NN;   
     VF=VFF;
     d=3;
   }
   
  void set(const R3 &P2,const R3 & P_Hat,const  baseFElement3 & K,int ll)
   { 
     P=P2;
     PHat=P_Hat;
     T3=&K.T;
     Th3=&K.Vh.Th; 
     region = T3->lab;
     label = ll;
     t=(*Th)(T);
     v=f=e=-1;  
     N.x=0;   
     N.y=0;   
     N.z=0;   
     VF=0;  
     d=3;
   }
   
     void set(const R3 &P2,const R3 & P_Hat,const  baseFElement3 & K)
   { 
     P=P2;
     PHat=P_Hat;
     T3=&K.T;
     Th3=&K.Vh.Th; 
     region = T3->lab;
     v=f=e=-1;  
     N.x=0;   
     N.y=0;   
     N.z=0;   
     VF=0;  
     int ll[4],kk(0);
     if ( P_Hat.x<1.e-6) ll[kk++]=1;
     if ( P_Hat.y<1.e-6) ll[kk++]=2;
     if ( P_Hat.z<1.e-6) ll[kk++]=3;
     if ( P_Hat.x+P_Hat.y+P_Hat.z>0.999999) ll[kk++]=0;    
     if (kk==0) label=0;
     else if (kk==3)
      { 
	 v = 6-ll[0]-ll[1]-ll[2];// 3 = 0+1+2 sommet oppose
	 label=(*T)[v].lab;
      } 
     else  {
       //  on edge
       //ffassert(0); // a faire 
       /*
       e = ll[0]; 
       int i1,i2,I3;
       /*
       Th3->VerticesNumberOfEdge(K.T,e,i1,i2);
       const BoundaryEdge * be=Th3->TheBoundaryEdge(i1,i2);
       label= be ? be->lab : 0;
       */
       label=-1;// to say 
      }
   
     t=(*Th3)(T3);
     d=3;
   }

  void set(const  Mesh3 &aTh, const R3 &P2,const R3 & P_Hat,const  Tet & aK,const int ll,bool coutside=false)
   { 
     P=P2;
     PHat=P_Hat;
     T3=&aK;
     Th3=&aTh; 
     region = T3->lab;
     label = ll;
     t=v=f=e=-1;  
     N.x=0;   
     N.y=0;   
     N.z=0;   
    // area=mes;
     outside=coutside;
     VF=0;  
     d=3;
   }
   
  void setP(const Mesh3 * pTh,int tt,int ss)
   { 
     T3=&(*pTh)[tt];
     const Mesh3::Vertex & V=(*T3)[ss];
     P= V ;
     PHat = TetHat[ss];
     Th3=pTh; 
     region = T3->lab;
     label = V.lab;
     t=v=f=e=-1;
     v=ss;
     VF=0;
     d=3;  
   }
   

  // --------  
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
  void set(const Mesh & aTh,const R2 &P2,const R2 & P_Hat,const  Triangle & aK,int ll,const R2 &NN,int iedge,int VFF) {
    MeshPointBase::set(aTh,P2,P_Hat,aK,ll,NN,iedge,VFF);
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

  // 3d
  void set(const R3 &P2,const R3 & P_Hat,const  baseFElement3 & K,int ll,const R3 &NN,int iedge) {
     MeshPointBase::set(P2,P_Hat,K,ll,NN,iedge);
     other.unset();}   
  void set(const Mesh3 & aTh,const R3 &P2,const R3 & P_Hat,const  Tet & aK,int ll,const R3 &NN,int iedge) {
    MeshPointBase::set(aTh,P2,P_Hat,aK,ll,NN,iedge);
    other.unset();}
  void set(const Mesh3 & aTh,const R3 &P2,const R3 & P_Hat,const  Tet & aK,int ll,const R3 &NN,int iedge,int VFF) {
    MeshPointBase::set(aTh,P2,P_Hat,aK,ll,NN,iedge,VFF);
    other.unset();}
  void set(const R3 &P2,const R3 & P_Hat,const  baseFElement3 & K) {
     MeshPointBase::set(P2,P_Hat,K);
     other.unset();
    }
  void set(const  Mesh3 &aTh, const R3 &P2,const R3 & P_Hat,const  Tet & aK,
           const int ll,bool coutside=false) {
    MeshPointBase::set(aTh,P2,P_Hat,aK,ll,coutside);
    other.unset();
  }
  // fin 3d
  void set(R x=0.0,R y=0.0,R z=0.0) {  
     MeshPointBase::set(x,y,z);
     other.unset();} 
  void change(const R2 & PH,const Triangle & tt,int ll) {
     MeshPointBase::change(PH,tt,ll);
     other.unset(); }
    void change(const R3 & PH,const Tet & tt,int ll) {
	MeshPointBase::change(PH,tt,ll);
    other.unset(); }
  void setP(const Mesh * pTh,int tt,int ss) { 
      MeshPointBase::setP(pTh,tt,ss); 
      other.unset(); } 
  void setP(const Mesh3 * pTh,int tt,int ss) {  // 3d
      MeshPointBase::setP(pTh,tt,ss); 
      other.unset(); } 

   bool operator==(const MeshPoint & mp) const {
      return T == mp.T &&  P.x == mp.P.x && P.y == mp.P.y 
          && P.z == mp.P.z ;}
  bool  SetAdj() {
     if (!(Th && T && t >=0 && e>=0)) return false;//  modif 
     if(VF==0)
     {
     int ieo=e,to=t,ie=e;
     
     t=Th->ElementAdj(t,ie);
     e=ie;
     if ( t == to && t >= 0 ) return false;
     int io0=VerticesOfTriangularEdge[ieo][0];
     //     int io1=VerticesOfTriangularEdge[ieo][1];
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
     } 
     else
     { //  
       VF = 1 + VF%2; 
       ffassert(0); // a faire 
     }  
     return true;         
   }
  };

ostream & operator << ( ostream &,const  MeshPoint & )  ;    
inline static MeshPoint* MeshPointStack(Stack s) {void * p= static_cast<void **>(s)[MeshPointStackOffset];throwassert(p); return static_cast<MeshPoint*>( p);}
inline static void MeshPointStack(Stack s,MeshPoint* mp) {*static_cast<MeshPoint**>(static_cast<void *>(s)) = mp;}
#ifdef NEWFFSTACK
inline static MeshPoint* MeshPointStack(void * s) {void * p= static_cast<void **>(s)[MeshPointStackOffset];throwassert(p); return static_cast<MeshPoint*>( p);}
inline static void MeshPointStack( void *s,MeshPoint* mp) {*static_cast<MeshPoint**>(static_cast<void *>(s)) = mp;}

#endif

}
