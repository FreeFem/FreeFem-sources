// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY: Bamg: Bidimensional Anisotrope Mesh Generator
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

// #define TRACETRIANGLE 3
//#pragma dont_inline on
//#pragma global_optimizer off
//#pragma inline_depth(0)

#undef NDEBUG
extern long verbosity ;
//#define strcasecmp strcmp
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include <iomanip>
using namespace std;

#include "Meshio.h"
#include "Mesh2.h"
#include "QuadTree.h"
#include "SetOfE4.h"

#include "rgraph.hpp"
#include "fem.hpp"
#include "AFunction.hpp"
#include "BamgFreeFem.hpp"
#include "RNM.hpp"
#include "FESpace.hpp"
#include "MeshPoint.hpp"


Fem2D::Mesh *bamg2msh( bamg::Triangles* tTh,bool renumbering)
{ 
  using namespace bamg;
  bamg::Triangles & th (*tTh);
  tTh->ReNumberingTheTriangleBySubDomain(!renumbering);//  just compress 
  //tTh->NbRef++;
  Int4  i,j,k=0;
  int nv  =  tTh->nbv;
  int nt  =   tTh->nbt - tTh->NbOutT;
  int neb =   tTh->nbe;
  
  int nbcrakev = 0;
  tTh->ReMakeTriangleContainingTheVertex();
  Fem2D::Triangle * t =  new Fem2D::Triangle[nt]  ;
  Fem2D::BoundaryEdge * b_e = new Fem2D::BoundaryEdge[neb];
  
  Fem2D::Vertex vbase;        
  Fem2D::Vertex *vb(&vbase);
  if (verbosity>5)
    cout << "  -- Before cracking mesh:  Nb Triangles = " << nt << " Nb of Vertices " << nv << endl;
  for ( i=0;i<nt;i++) // unset all triangles 
    for (j=0;j<3;j++)
      t[i](j)=0;
  //  nv=0;    
  for (int iv=0;iv<th.nbv;iv++) // vertex 
    {
      // cout << iv << " : " ;
      const Vertex & v(th[iv]); 
      int kk=0; // nb cracked
      int kc=0; 
      int kkk =0; // nb triangle  with same number 
      Triangle * tbegin = v.t;
      Fem2D::Vertex * vv = vb+iv;
      int i  = v.vint;       
      throwassert(tbegin && (i >= 0 ) && (i <3));
      // turn around the vertex v
      TriangleAdjacent ta(tbegin,EdgesVertexTriangle[i][0]);// previous edge
      int k=0;
      do {
        int kv = VerticesOfTriangularEdge[ta][1];
        k++; 
        Triangle * tt (ta);
        throwassert( &v == & (*  tt)[kv] );            
        if ( ta.Cracked() ) 
          {   // cout << " || "    ;                    
            if ( kk == 0) tbegin=ta,kkk=0;  //  begin by a cracked edge  => restart                
            if (  kkk ) { kc =1;vv = vb +  nv++;  kkk = 0; } // new vertex if use 
            kk++;
            // number of cracked edge view                 
          }
        if ( tt->link ) { // if good triangles store the value 
          int it = th.Number(tt);
          throwassert(it < nt);
          int iiv=vv-vb;
          t[it](kv) = vv;
          /*
          cout << it << " " << kv << " "<< iiv  << endl;
          if (&th(it)[kv] != &th[iiv])
             cout << it << " " << kv << " "<< iiv << " != " << th.Number(th(it)[kv]) << endl ;
          */
          kkk++;
        } else if (kk) { // crack + boundary 
          if (  kkk ) { kc =1;vv = vb +  nv++;  kkk = 0; } // new vertex if use 
        }
        
        ta = Next(ta).Adj(); 
      } while ( (tbegin != ta)); 
      throwassert(k);
      if (kc)  nbcrakev++;
    }
  Fem2D::Vertex * v = new Fem2D::Vertex[nv];
  //  set the vertices --
  for (i=0;i<nt;i++)
    { 
      for (j=0;j<3;j++)
        {
          throwassert( t[i](j) );             
          int k = t[i](j) - vb;
          t[i](j) = v+ k;
          throwassert(k>=0 && k < nv);
          Vertex & thv(th(i)[j]);
          v[k].x    =  thv.r.x;
          v[k].y    =  thv.r.y;
          v[k].lab  =  thv.ref();        
        }  
    }
  // warning in cracked edges 
  // construction of the edges --
  
  if (nbcrakev && verbosity>2)
    cout << " -- Nb of craked vertices = " << nbcrakev << " Nb of created vertices " << nv - th.nbv << endl;
  
  
  for (i=0;i<tTh->nbe;i++)
    {
      int i0=tTh->Number(tTh->edges[i][0]),i1=tTh->Number(tTh->edges[i][1]);
      throwassert(i0>=0 && i0 <nv);
      throwassert(i1>=0 && i1 <nv);      
      b_e[i]=Fem2D::BoundaryEdge(v,i0,i1,tTh->edges[i].ref);
    }      
  Int4 *reft = new Int4[tTh->nbt];
  Int4 nbref = tTh->ConsRefTriangle(reft);
  for( i=0,k=0;i<tTh->nbt;i++)
    if(tTh->triangles[i].link)
      { 
        
        Fem2D::R2 A(t[k][0]),B(t[k][1]),C(t[k][2]);
        t[k].area = (( B-A)^(C-A))*0.5 ;
        t[k].lab = tTh->subdomains[reft[i]].ref;  // a faire
        throwassert(k == i);
        k++;
      }
  delete [] reft;
  throwassert ( nt == k);
  tTh->ReMakeTriangleContainingTheVertex();
  
  if (verbosity)
    cout << "  --  mesh:  Nb of Triangles = "  << setw(6) <<  nt << ", Nb of Vertices " << nv << endl;
  
  {  
    Fem2D::Mesh *m = new Fem2D::Mesh(nv,nt,neb,v,t,b_e);
    if (renumbering) m->renum();
    m->MakeQuadTree();
    return m;
  }
}

Fem2D::Mesh *bamg2msh(const bamg::Geometry &Gh)
{ 
  // ------------------
  int nv=  Gh.nbv;
  int neb=Gh.nbe;
  Fem2D::Triangle * t = 0 ;
  Fem2D::BoundaryEdge * b_e = new Fem2D::BoundaryEdge[neb];
  Fem2D::Vertex *v = new Fem2D::Vertex[nv]  ;
  for (int i=0;i<nv;i++)
    {
      const bamg::GeometricalVertex & vg( Gh[i]);
      v[i].x=vg.r.x;
      v[i].y=vg.r.y;
      v[i].lab=vg.ref();
      
    }     
  for (int ie=0;ie<neb;ie++)
    {
      const bamg::GeometricalEdge & e= Gh(ie);
      int i0=Gh.Number(e[0]),i1=Gh.Number(e[0]);
      b_e[ie]= Fem2D::BoundaryEdge(v,i0,i1,e.ref);
    }   
  
  {  
    Fem2D::Mesh *m = new Fem2D::Mesh(nv,0,neb,v,t,b_e);
    m->MakeQuadTree();
    return m;
  }
  // ------------------
}





bamg::Triangles * msh2bamg(const Fem2D::Mesh & Th,double cutoffradian)
  
{
  using namespace bamg;
  Triangles *Tn=new Triangles(Th.nv);
  
  Tn->nbv = Th.nv;
  Tn->nbt = Th.nt;
  Tn->nbe = Th.neb;
  Tn->name= new char[strlen("msh2bamg")+1];
  strcpy(Tn->name,"msh2bamg");
  //  Tn->triangles = new Triangle [Tn->nbtx];
  throwassert(Tn->triangles);
  //  Tn->vertices = new Vertex [Tn->nbvx];
  //  Tn->ordre = new (Vertex* [Tn->nbvx]);
  Tn->edges = new Edge [Th.neb];
 
  Int4 i;
  for (i = 0; i < Th.nv; i++)
    {
      Tn->vertices[i].r.x = Th(i).x;
      Tn->vertices[i].r.y = Th(i).y;
      Tn->vertices[i].ReferenceNumber = Th(i).lab;
    }
  
  //  Int4 i1 [nbt],i2 [nbt],i3 [nbt];
  for (i = 0; i < Th.nt; i++)
    {
      int i1 = Th(Th[i][0]);
      int i2 = Th(Th[i][1]);
      int i3 = Th(Th[i][2]);
      Tn->triangles[i]= Triangle( Tn,i1 ,i2 ,i3 );
      Tn->triangles[i].color = Th[i].lab;
    }
  for (i = 0; i < Th.neb; i++)
    {
      Tn->edges[i].v[0] = Tn->vertices + Th(Th.bedges[i][0]);
      Tn->edges[i].v[1] = Tn->vertices + Th(Th.bedges[i][1]);
      Tn->edges[i].ref = Th.bedges[i].lab;
    }
  //  Real8 cutoffradian = -1;
  Tn->ConsGeometry(cutoffradian);
  Tn->Gh.AfterRead();    
  Tn->SetIntCoor();
  Tn->FillHoleInMesh();
  return Tn;
}


bamg::Triangles * msh2bamg(const Fem2D::Mesh & Th,double cutoffradian, 
                           int  nbdfv, int * ndfv,int  nbdfe, int * ndfe)
{
  using namespace bamg;
  Triangles *Tn=new Triangles(Th.nv);
  assert(nbdfe==0 && nbdfv==0); // a faire pour les maillages periodique 
  
  Tn->nbv = Th.nv;
  Tn->nbt = Th.nt;
  Tn->nbe = Th.neb;
  Tn->name= new char[strlen("msh2bamg")+1];
  strcpy(Tn->name,"msh2bamg");
  //  Tn->triangles = new Triangle [Tn->nbtx];
  throwassert(Tn->triangles);
  //  Tn->vertices = new Vertex [Tn->nbvx];
  //  Tn->ordre = new (Vertex* [Tn->nbvx]);
  Tn->edges = new Edge [Th.neb];
  
  Int4 i;
  for (i = 0; i < Th.nv; i++)
    {
      Tn->vertices[i].r.x = Th(i).x;
      Tn->vertices[i].r.y = Th(i).y;
      Tn->vertices[i].ReferenceNumber = Th(i).lab;
    }
  
  //  Int4 i1 [nbt],i2 [nbt],i3 [nbt];
  for (i = 0; i < Th.nt; i++)
    {
      int i1 = Th(Th[i][0]);
      int i2 = Th(Th[i][1]);
      int i3 = Th(Th[i][2]);
      Tn->triangles[i]= Triangle( Tn,i1 ,i2 ,i3 );
      Tn->triangles[i].color = Th[i].lab;
    }
  for (i = 0; i < Th.neb; i++)
    {
      Tn->edges[i].v[0] = Tn->vertices + Th(Th.bedges[i][0]);
      Tn->edges[i].v[1] = Tn->vertices + Th(Th.bedges[i][1]);
      Tn->edges[i].ref = Th.bedges[i].lab;
    }
  //  Real8 cutoffradian = -1;
  Tn->ConsGeometry(cutoffradian);
  Tn->Gh.AfterRead();    
  Tn->SetIntCoor();
  Tn->FillHoleInMesh();
  return Tn;
}



Fem2D::Mesh *  BuildMesh(Stack stack, E_BorderN const * const & b,bool justboundary) 
{
  using namespace bamg;
  using bamg::Abs;
  using bamg::Max;
  using bamg::Min;
  using bamg::Pi;
  Fem2D::MeshPoint & mp (*Fem2D::MeshPointStack(stack)), mps = mp;
  
  int nbvx=0,nbe=0,nbsd=0;
  for (E_BorderN const * k=b;k;k=k->next)
    {long n=  Max(1L,Abs(k->Nbseg(stack))); ;
    nbvx += n+1;
    nbe += n;
    nbsd++;
    }
  
  Geometry * Gh =  new Geometry;
  
  if(verbosity>2)
    cout <<"\t\t"  << "  Begin: ConstGeometry from Border"  << endl;
  //throwassert(empty());
  const char * filename = "FREEFEM.gh";
  Gh->name=new char [strlen(filename)+1];
  strcpy(Gh->name,filename);
  Real8 Hmin = HUGE_VAL;// the infinie value 
  Int4 hvertices =0;
  Int4 i,nn,n;
  Int4 dim=0;
  Gh->MaximalAngleOfCorner =30.00*Pi/180.0;
  Gh->nbv = 0;
  Gh->nbvx = nbvx;
  
  Gh->nbe = nbe;
  Gh->edges = new GeometricalEdge[Gh->nbe];
  bamg::Vertex *vertices =  new GeometricalVertex[Gh->nbvx];
  double lmin= HUGE_VAL;
  //  generation des points et des lignes 
  i=0;
  for (E_BorderN const * k=b;k;k=k->next)
    {
      assert(k->b->xfrom); // a faire 
      double & t = *  k->var(stack);
      double a(k->from(stack)),b(k->to(stack));
      n=Max(Abs(k->Nbseg(stack)),1L);     
      t=a;
      double delta = (b-a)/n;
      for ( nn=0;nn<=n;nn++,i++, t += delta)
        {
          if (nn==n) t=b; // to remove roundoff error 
          mp.label = k->label();
          k->code(stack); // compute x,y, label
          // cout << " ----- " << i << " " << mp.P.x << " " << mp.P.y << endl;
          vertices[i].r.x=mp.P.x;
          vertices[i].r.y=mp.P.y;
          vertices[i].ReferenceNumber=  mp.label;
          vertices[i].color = i;
          if (nn>0) {
            lmin=min(lmin,Norme2_2( vertices[i].r-vertices[i-1].r));
          }     
        }
    }
  lmin = sqrt(lmin);
  double eps = (lmin)/16.; 
  int nbvprev = i;
  long nbv=0;
  Gh->pmin =  vertices[0].r;
  Gh->pmax =  vertices[0].r;
  // recherche des extrema des vertices pmin,pmax
  for (i=0;i<nbvprev;i++) 
    {
      Gh->pmin.x = Min(Gh->pmin.x,vertices[i].r.x);
      Gh->pmin.y = Min(Gh->pmin.y,vertices[i].r.y);
      Gh->pmax.x = Max(Gh->pmax.x,vertices[i].r.x);
      Gh->pmax.y = Max(Gh->pmax.y,vertices[i].r.y);
    }
  Gh->coefIcoor= (MaxICoor)/(Max(Gh->pmax.x-Gh->pmin.x,Gh->pmax.y-Gh->pmin.y));
  Icoor1 epsI = (Icoor1) (Gh->coefIcoor*eps);
  throwassert(Gh->coefIcoor >0);
  if (verbosity>2) 
    {
      cout <<"\t\t"  << "     Geom: min="<< Gh->pmin << "max ="<< Gh->pmax 
           << " hmin = " << Gh->MinimalHmin() <<  endl;
    }
  nbv = 0;        
  {  // find common point 
    QuadTree quadtree;
    Metric Id(1.);
    for ( i=0;i<nbvprev;i++) 
      { 
        vertices[i].i = Gh->toI2(vertices[i].r);
        vertices[i].m = Id;
        Vertex *v= quadtree.ToClose(vertices[i],eps,epsI,epsI) ;
        // quadtree.NearestVertex(vertices[i].i.x,vertices[i].i.y); 
        if( v && Norme1(v->r - vertices[i]) < eps )
          { vertices[i].color=v->color; }
        else  {quadtree.Add(vertices[i]);
        vertices[i].color = nbv++;}
      }
    /*   
         if (nbvprev-nbv==0)
         {
         reffecran();
         bamg::R2 O((Gh->pmin+Gh->pmax)/2),D(Gh->pmax-Gh->pmin);
         cadreortho(O.x,O.y,Max(D.x,D.y)*1.1);
         xGrafCoef = Gh->coefIcoor;
         yGrafCoef = Gh->coefIcoor;
         xGrafOffSet = Gh->pmin.x;
         yGrafOffSet = Gh->pmin.y;  
         quadtree.Draw();
         for (int i=0;i<nbvprev;i++)
         {
         rmoveto(vertices[i].r.x,vertices[i].r.y);
         
         char buf[100];
         sprintf(buf,"%d",i);
         plotstring(buf);
         }
         rattente(1);
         }
    */   
  } // to delete quadtree
  cout << " Nb of common points " << nbvprev-nbv <<endl;
  
  Gh->nbvx = nbv;
  Gh->nbv = nbv;
  
  Gh->vertices = new GeometricalVertex[nbv];
  throwassert(Gh->nbvx >= Gh->nbv);
  Gh->nbiv = Gh->nbv;
  Int4 k=0;
  const Direction NoDirOfSearch;
  //  compression of points    
  int kkk;     
  for ( i=0,kkk=0;kkk<nbvprev;kkk++) 
    {
      if (vertices[kkk].color == i) // if new points 
        {
          Gh->vertices[i].r.x = vertices[kkk].r.x ;
          Gh->vertices[i].r.y = vertices[kkk].r.y;
          //Gh->vertices[i].link = Gh->vertices + i;
          throwassert(Gh->vertices[i].IsThe());
          Gh->vertices[i].ReferenceNumber = vertices[kkk].ReferenceNumber  ;
          Gh->vertices[i].DirOfSearch = NoDirOfSearch;
          Gh->vertices[i].color =0;
          Gh->vertices[i].Set();
          //  vertices[i].SetCorner();
          //  vertices[i].SetRequired();
          i++;
        }
    }
  throwassert(i==nbv);
  R2 zero2(0,0);
  if(verbosity>5) 
    cout <<"\t\t"  << "     Record Edges: Nb of Edge " << Gh->nbe <<endl;
  throwassert(Gh->edges);
  throwassert (Gh->nbv >0); 
  Real4 *len =0;
  if (!hvertices) 
    {
      len = new Real4[Gh->nbv];
      for(i=0;i<Gh->nbv;i++)
        len[i]=0;
    }
  int nnn=0;
  i=0;
  for (E_BorderN const * k=b;k;k=k->next)
    {
      double & t = *  k->var(stack);
      double a(k->from(stack)),b(k->to(stack));
      n=Max(Abs(k->Nbseg(stack)),1L);     
      double delta = (b-a)/n;
      t=a+delta/2; 
      for ( nn=0;nn<n;nn++,i++, t += delta)
        {
          
          mp.label = k->label();
          k->code(stack); 
          Int4 i1 =  vertices[nnn].color, i2 =  vertices[++nnn].color;
          throwassert(i1 >= 0 && i1 < nbv);
          throwassert(i2 >= 0 && i2 < nbv);
          Gh->edges[i].ref = mp.label;
          Gh->edges[i].v[0]=  Gh->vertices + i1;
          Gh->edges[i].v[1]=  Gh->vertices + i2;
          R2 x12 = Gh->vertices[i2].r-Gh->vertices[i1].r;
          Real8 l12=Norme2(x12);
          Gh->edges[i].tg[0]=zero2;
          Gh->edges[i].tg[1]=zero2;
          Gh->edges[i].SensAdj[0] = Gh->edges[i].SensAdj[1] = -1;
          Gh->edges[i].Adj[0] = Gh->edges[i].Adj[1] = 0;
          Gh->edges[i].flag = 0;
          Gh->edges[i].link=0;
          if (!hvertices) 
            {
              Gh->vertices[i1].color++;
              Gh->vertices[i2].color++;
              len[i1] += l12;
              len[i2] += l12;
            }
          
          Hmin = Min(Hmin,l12);
        }
      nnn++; 
    }
  
  delete [] vertices; vertices=0;
  
  throwassert(nnn==nbvprev);
  throwassert(i==Gh->nbe);
  // definition  the default of the given mesh size 
  if (!hvertices) 
    {
      for (i=0;i<Gh->nbv;i++) 
        if (Gh->vertices[i].color > 0) 
          Gh->vertices[i].m=  Metric(len[i] /(Real4) Gh->vertices[i].color);
        else 
          Gh->vertices[i].m=  Metric(Hmin);
      delete [] len;
      
      if(verbosity>3) 
        cout <<"\t\t"  << "     Geom Hmin " << Hmin << endl;
    }
  
  Gh->NbSubDomains=nbsd;
  if (Gh->NbSubDomains>0) 
    {
      Gh->subdomains = new GeometricalSubDomain[  Gh->NbSubDomains];
      Int4 i1=0;
      i=0;
      for (E_BorderN const * k=b;k;k=k->next,i++)
        { long Nbseg =k->Nbseg(stack);
        long n=  Max(1L,Abs(Nbseg)); 
        Gh->subdomains[i].sens = Nbseg >0 ? 1 : -1;
        Gh->subdomains[i].edge=Gh->edges + i1;
        Gh->subdomains[i].ref = i;
        i1 += n;
        }
    }
  Gh->NbEquiEdges=0;
  Gh->NbCrackedEdges=0;
  Fem2D::Mesh * m=0;
  if (justboundary)
    m=bamg2msh(*Gh);
  else 
    {
      Gh->AfterRead();   
      Triangles *Th = new Triangles( (Gh->nbv*Gh->nbv)/9 +1000 ,*Gh);
      m=bamg2msh(Th,true);      
      delete Th;
    }
  delete Gh;
  /* deja fait  dans bamg2msh
     Fem2D::R2 Pn,Px;
     m->BoundingBox(Pn,Px);
     m->quadtree=new Fem2D::FQuadTree(m,Pn,Px,m->nv);
  ---------- */
  mp=mps;   
  m->decrement();
  return m;   
}

void E_BorderN::BoundingBox(Stack stack,double  &xmin,double & xmax, double & ymin,double & ymax) const
{  
  Fem2D::MeshPoint & mp (*Fem2D::MeshPointStack(stack)), mps = mp;
  for (E_BorderN const * k=this;k;k=k->next)
    {
      assert(k->b->xfrom); // a faire 
      double & t = *  k->var(stack);
      double a(k->from(stack)),b(k->to(stack));
      long n=Max(Abs(k->Nbseg(stack)),1L);     
      t=a;
      double delta = (b-a)/n;
      for (int  nn=0;nn<=n;nn++, t += delta)
        {
          if (nn==n) t=b; // to remove roundoff error 
          mp.label = k->label();
          k->code(stack); // compute x,y, label
          xmin=Min(xmin,mp.P.x);
          xmax=Max(xmax,mp.P.x);
          ymin=Min(ymin,mp.P.y);
          ymax=Max(ymax,mp.P.y);
        }
    }
  mp=mps; 
}  
void E_BorderN::Plot(Stack stack) const
{  
  using Fem2D::R2;
  
  Fem2D::MeshPoint & mp (*Fem2D::MeshPointStack(stack)), mps = mp;
  float x0,x1,y0,y1;
  getcadre(x0,x1,y0,y1);
  float h= (x1-x0)*0.01;
  
  for (E_BorderN const * k=this;k;k=k->next)
    {
      assert(k->b->xfrom); // a faire 
      double & t = *  k->var(stack);
      double a(k->from(stack)),b(k->to(stack));
      long n=Max(Abs(k->Nbseg(stack)),1L);     
      t=a;
      double delta = (b-a)/n;
      R2 P,Po;
      for (int  nn=0;nn<=n;nn++, t += delta)
        {
          if (nn==n) t=b; // to remove roundoff error 
          mp.label = k->label();
          mp.P.z=0;
          k->code(stack); // compute x,y, label
          P=mp.P;
          couleur(2+mp.label);
          if(nn!=0) { LineTo(P);
          R2 uv(Po,P);
          double l = Max(sqrt((uv,uv)),1e-20);
          
          R2 dd = uv*(-h/l);
          R2 dn = dd.perp()*0.5;
          
          LineTo(P+dd+dn);
          MoveTo(P+dd-dn);
          LineTo(P);}
          else {
            DrawMark(mp.P,0.01);
            MoveTo(mp.P);
            
          }
          
          
          //  cout << k->label()<< " " << nn << ", x,y = " << mp.P.x  << " , " << mp.P.y << endl;
          Po=P;
        }
      DrawMark(mp.P,0.01);
      MoveTo(mp.P);
    }
  
  mp=mps; 
}

Fem2D::Mesh *  BuildMeshBorder(Stack stack, E_BorderN const * const & b) 
{
  return BuildMesh(stack,b,true);
}
Fem2D::Mesh *  BuildMesh(Stack stack, E_BorderN const * const & b) 
{
  return BuildMesh(stack,b,false);
}

Fem2D::Mesh *  ReadTriangulate( string  * const & s) {
  using namespace Fem2D;
  KN<R2> xy;
  char c;
  int nv;
  for(int step=0;step<2;step++)
    {
      nv=0;
      ifstream f(s->c_str());
      if(!f) {cerr <<" Error openning file " << *s << endl;
      ExecError("Openning file ");}
      while (f.good())
        {
          R2 P;
          f >> P ;
          if (!f.good()) break;
          if (step) xy[nv]=P;
          nv++;
          while (f.get(c) &&  (c!='\n' && c!='\r' ) ) 0; // eat until control (new line
        } 
      if (!step && nv ) xy.init(nv); // alloc the array        
    }
  if(verbosity)
    cout << " we read  " << nv << " coordinates  xy "<< endl;
  
  Mesh * m=new Mesh(nv,xy); 
  m->MakeQuadTree();
  m->decrement();
  delete s;
  return m;
  
}
Fem2D::Mesh *  ReadMeshbamg( string * const & s) {
  using bamg::Triangles;
  Triangles * bTh= new Triangles(s->c_str());
  // bTh->inquire();
  Fem2D::Mesh * m=bamg2msh(bTh,false);// no renum
  delete bTh;
  /*  deja fait 
      Fem2D::R2 Pn,Px;
      m->BoundingBox(Pn,Px);
      m->quadtree=new Fem2D::FQuadTree(m,Pn,Px,m->nv); 
  */
  delete s;
  m->decrement();
  return m;
}
