#include <cmath>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <fstream>
//#include <strstream.h>
//using namespace std;  //introduces namespace std
#include "RNM.hpp"
#include "rgraph.hpp"
#include "fem.hpp"
#include "FQuadTree.hpp"

inline int NbOfSubTriangle(int k)
{  
   assert(k>0);
   return  k*k;
}


inline int NbOfSubInternalVertices(int k)
{ 
   assert(k>0);
   int  r= (k+2)*(k+1)/2;
   assert(r>=0);
   return r;
}


//  for the  mortar elements
 Mesh::Mesh(const Mesh & Th,int * split,bool WithMortar)
 { //  routine complique 
   //  count the number of elements
   
    area=Th.area;
    BoundaryAdjacencesHead=0;
    BoundaryAdjacencesLink=0;
    BoundaryEdgeHeadLink=0;
    quadtree =0;
    NbMortars=0;
    mortars=0;
    TheAdjacencesLink =0;
    quadtree =0;
    bedges=0;
    neb =0;
   R2 Pmin,Pmax;
   Th.BoundingBox(Pmin,Pmax);
   nt=0;
   for (int i=0;i<Th.nt;i++)
     nt += NbOfSubTriangle(split[i]);
   triangles = new Triangle[nt];
   assert(triangles);
   //  computation of thee numbers of vertices
   //  on decoupe tous betement
   // et on recolle le points ensuite 
   // avec le quadtree qui sont sur les aretes ou les sommets 
   //  calcul du magorant du nombre de sommets 
   int nvmax=Th.nv;
//  attention il faut supprime les aretes frontiere  interne    
   for (int ieb=0,jj;ieb<Th.neb;ieb++)
       neb += split[Th.BoundaryTriangle(ieb,jj)];
   int nebmax=neb;
   
   for (int i=0;i<Th.nt;i++)
     nvmax += NbOfSubInternalVertices(split[i]) -3;

  //  compute the minimal Hsize of the new mesh  
  R hm = Th[0].h();
//  cout << " hm " << hm << endl;
  for (int it=0;it<Th.nt;it++)
    { 
      assert(split[it]>0 && split[it]<=64);
//      cout << " it = " <<it << " h " <<  Th[it].h() << " " << split[it] << " hm = " << hm <<  endl;
      hm=Min(hm,Th[it].h()/(R) split[it]);
    }
   R seuil=hm/4.0;
   cout << " seuil = " <<  seuil << " hmin = " << hm <<  endl; 
   assert(seuil>1e-15);
   vertices = new Vertex[nvmax];
   assert( vertices && Th.nv <= nvmax);
   
   //  copy of the old vertices 
   for (int i=0;i<Th.nv;i++)
     vertices[i] = Th(i);
     
   nv = Th.nv;
   
   quadtree = new FQuadTree(this,Pmin,Pmax,nv); //  put all the old vertices in the quadtree 
   bedges = new BoundaryEdge[nebmax];
   assert(bedges);
//  generation of the boundary edges
   neb =0;
   for (int ieb=0,jj;ieb<Th.neb;ieb++)
       {  
          int kold = Th.BoundaryTriangle(ieb,jj);
          int n=split[kold];
          
          int n1=n+1;
          BoundaryEdge & be(Th.bedges[ieb]);
          Vertex *pva= vertices + (&be[0]-Th.vertices);
          Vertex *pvb= vertices + (&be[1]-Th.vertices);
          Vertex *pv0=pva;
          R2 A(*pva),B(*pvb);
          R la = 1,lb=0, delta=1.0/n;
          
          for (int j=1;j<n;j++) 
           { 
             la-=delta;
             lb+=delta;
             assert(nv<nvmax);
             Vertex *pv1= vertices + nv;
             
             (R2 &) *pv1 = A*la + B*lb;
             (Reference &) *pv1 = (Reference) be;
             quadtree->Add(*pv1);
             nv++;
             assert(neb<nebmax);
             bedges[neb].vertices[0]=pv0;
             bedges[neb].vertices[1]=pv1;
             (Reference &) bedges[neb]= be;
             neb++;
             pv0=pv1;
           }
         assert(neb<nebmax);
         bedges[neb].vertices[0]=pv0;
         bedges[neb].vertices[1]=pvb;
         (Reference &) bedges[neb]= be;
         neb++;
       }    
    assert(neb==nebmax);
    
//   create the new vertices and the new triangle 
   int kerr=0,kt=0;
   for (int K=0;K<Th.nt;K++)
    { // cout << K << endl;
      Triangle &T(Th[K]);
      R2 A(T[0]);
      R2 B(T[1]);
      R2 C(T[2]);
      long N=split[K];
      long N2=N*N;
      int vt[3];
      for (int n=0;n<N2;n++,kt++) //  loop on all sub triangle
        {
          //(Reference&) triangles[kt] = (Reference&) T; 
          for(int j=0;j<3;j++) //  Loop on the 3 vertices
           { R2 PTj=SubTriangle(N,n,j);
             R la=1-PTj.x-PTj.y,lb=PTj.x,lc=PTj.y;
             R2 Pj= A*la+B*lb+C*lc; 
             Vertex *pV;
             pV=quadtree->ToClose(Pj,seuil);
             if(!pV)
             { // new vertex
               //  cout << "    -- " << nv << "  New Vertices " << Pj << " n=" << n << " j=" << j << " " << PTj << endl;
                 assert(nv<nvmax);
                 vertices[nv]=Pj;
                 (Reference&) vertices[nv]=0; //  Internal vertices 
                 pV = vertices + nv;
                 quadtree->Add(*pV);
                 nv++;                  
               }  //  end of new vertex
               vt[j]=number(pV); 
           //  triangles[kt].SetVertex(j,pV);
           } // for(int j=0;j<3;j++)
          // cout << kt << " " << n << endl;
           R2 A=vertices[vt[0]];
           R2 B=vertices[vt[1]];
           R2 C=vertices[vt[2]]; 
           R a = (( B-A)^(C-A))*0.5;
           if (a>0) 
             triangles[kt] = Triangle(vertices,vt[0],vt[1],vt[2],T.ref);
           else 
             triangles[kt] = Triangle(vertices,vt[0],vt[2],vt[1],T.ref);
           
        }   // end loop on all sub triangle
   
     } //  end 
    
   cout << " -- Nb of vertices       " << nv << endl; 
   cout << " -- Nb of triangle       " << nt << endl; 
   cout << " -- Nb of boundary edges " << neb << endl; 
//  
   if (!WithMortar)   
    { // 
      long nba = neb;
     // 
     long nbsd = 1; // bofbof 
      // faux just pour un test 
      long  *sd;
      nbsd=1;
      sd=new long[2*nbsd];
      sd[0]=-1;
      sd[1]=Th[0].ref;         
       
      long nbs=nv;
      long nbsmax=nv;
      long            err = 0, nbsold = nbs;       
      long           *c = 0;
      long           *tri = 0;
      long           *nu = 0;
      long           *reft = 0;
      float          *cr = 0;
      float          *h = 0;
      long nbtmax = 2 * nbsmax;
      long * arete  = new long[2*nba]; 
      nu = new long[6*nbtmax];
      c = new long[2*nbsmax];
      tri = new long[(4 * nbsmax + 2 * nbsd)];
      reft = new long[nbtmax];
      cr = new float[(2 * nbsmax + 2)];
      h = new float[nbsmax];
      for(int i=0,k=0; i<nv; i++)
       { 
          cr[k++]  =vertices[i].x;
          cr[k++]=vertices[i].y;
          h[i]=1; 
       }
      for (int i=0,k=0 ;i<neb;i++)
       {
          arete[k++] =number(bedges[i][0])+1;
          arete[k++] =number(bedges[i][1])+1;
       }
       
       
extern int 
mshptg_ (float *cr, float *h, long *c, long *nu, long *nbs, long nbsmx, long *tri,
	 long *arete, long nba, long *sd,
	 long nbsd, long *reft, long *nbt, float coef, float puis, long *err);
      
      long nbt=0;
      mshptg_ (cr, h, c, nu, &nbs, nbs, tri, arete, nba, (long *) sd, nbsd, reft, &nbt, .25, .75, &err);
      assert(err==0 && nbt !=0);
      delete [] triangles;
      nt = nbt;
      cout << " Nb Triangles = " << nbt << endl;
      triangles = new Triangle[nt];
      for(int i=0,k=0;i<nt;i++,k+=3)
         triangles[i]=Triangle(vertices,nu[k]-1,nu[k+1]-1,nu[k+2]-1,reft[i]);
      delete [] arete;
      delete [] nu;
      delete [] c;
      delete [] tri;
      delete [] reft;
      delete [] cr;
      delete [] h;
      delete [] sd;
      
    }
   ConsAdjacence();
}
