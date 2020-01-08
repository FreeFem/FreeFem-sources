// ORIG-DATE:     Dec 97
// -*- Mode : c++ -*-
//
// SUMMARY  :
// USAGE    :
// ORG      :
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
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
 */

extern long verbosity ;
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <fstream>
#include "RNM.hpp"

using namespace std;

#include "Meshio.h"
#include "Mesh2.h"
#include "QuadTree.h"
#include "SetOfE4.h"

#include "rgraph.hpp"
#include "fem.hpp"
#include "AFunction.hpp"
#include "BamgFreeFem.hpp"
#include "FESpace.hpp"
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "MeshPoint.hpp"
#include "PlotStream.hpp"
#include <set>
const Fem2D::Mesh *bamg2msh( bamg::Triangles* tTh,bool renumbering)
{
  using namespace bamg;
  bamg::Triangles & th (*tTh);
  tTh->ReNumberingTheTriangleBySubDomain(!renumbering);//  just compress
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
  for (int iv=0;iv<th.nbv;iv++) // vertex
    {
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
          {
            if ( kk == 0) tbegin=ta,kkk=0;  //  begin by a cracked edge  => restart
            if (  kkk ) { kc =1;vv = vb +  nv++;  kkk = 0; } // new vertex if use
            kk++;
            // number of cracked edge view
          }
        if ( tt->link ) { // if good triangles store the value
          int it = th.Number(tt);
          throwassert(it < nt);
          t[it](kv) = vv;
          kkk++;
        } else if (kk) { // crack + boundary
          if (  kkk ) { kc =1;vv = vb +  nv++;  kkk = 0; } // new vertex if use
        }

        ta = Next(ta).Adj();
      } while ( (tbegin != ta));
      throwassert(k);
      if (kc)  nbcrakev++;
    }
  double badvalue=12345e100;
  Fem2D::Vertex * v = new Fem2D::Vertex[nv];
    for(int i=0; i<nv; ++i)
        v[i].x    =  badvalue;
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
    int kerr=0;

    for(int i=0; i<nv; ++i)
        if(v[i].x    ==  badvalue) kerr++;
    if(kerr)
    {
        cerr << " bamg2msh: Error: missing "<<kerr << " vertices ???? " << endl;
        cerr << " I  May be: Some giving point are outside the domain ???? " << endl;
        ExecError("Buildmesh : bamg2msh missing vertices");
    }
  // warning in cracked edges
  // construction of the edges --

  if (nbcrakev && verbosity>2)
    cout << "  -- Nb of craked vertices = " << nbcrakev << " Nb of created vertices " << nv - th.nbv << endl;


  for (i=0;i<tTh->nbe;i++)
    {
      int i0=tTh->Number(tTh->edges[i][0]),i1=tTh->Number(tTh->edges[i][1]);
      throwassert(i0>=0 && i0 <nv);
      throwassert(i1>=0 && i1 <nv);
      b_e[i]=Fem2D::BoundaryEdge(v,i0,i1,tTh->edges[i].ref);
    }
  Int4 *reft = new Int4[tTh->nbt];
  //Int4 nbref =
  long nerr=0;
  tTh->ConsRefTriangle(reft);
  for( i=0,k=0;i<tTh->nbt;i++)
    if(tTh->triangles[i].link)
      {

        Fem2D::R2 A(t[k][0]),B(t[k][1]),C(t[k][2]);
        t[k].area = (( B-A)^(C-A))*0.5 ;
          nerr += t[k].area <=0;
        t[k].lab = tTh->subdomains[reft[i]].ref;  // a faire
        throwassert(k == i);
        k++;
      }
  delete [] reft;
  throwassert ( nt == k);
  tTh->ReMakeTriangleContainingTheVertex();

  if (verbosity)
    cout << "  --  mesh:  Nb of Triangles = "  << setw(6) <<  nt << ", Nb of Vertices " << nv << endl;
  if(nerr)
  {
      cerr << "Fatal error: Number of Negative triangles "<< nerr << endl;
      delete [] v;
      delete [] t;
      delete [] b_e;
      ExecError("bamg2mesh error negative triangle ");
  }
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


 bamg::Triangles * msh2bamg(const Fem2D::Mesh & Th,double cutoffradian,long * reqedgeslab,int nreqedgeslab)

{
  using namespace bamg;
  Triangles *Tn=new Triangles(Th.nv);
  Tn->nbv = Th.nv;
  Tn->nbt = Th.nt;
  Tn->nbe = Th.neb;
  Tn->name= new char[strlen("msh2bamg")+1];
  strcpy(Tn->name,"msh2bamg");
  throwassert(Tn->triangles);
  Tn->edges = new Edge [Th.neb];

  Int4 i;
  Metric Mid(1.);
  for (i = 0; i < Th.nv; i++)
    {
      Tn->vertices[i].r.x = Th(i).x;
      Tn->vertices[i].r.y = Th(i).y;
	Tn->vertices[i].m=Mid;
      Tn->vertices[i].ReferenceNumber = Th(i).lab;
    }

  for (i = 0; i < Th.nt; i++)
    {
      int i1 = Th(Th[i][0]);
      int i2 = Th(Th[i][1]);
      int i3 = Th(Th[i][2]);
      Tn->triangles[i]= Triangle( Tn,i1 ,i2 ,i3 );
      Tn->triangles[i].color = Th[i].lab;
    }
    //  Real8 cutoffradian = -1;
    // add code   un change boundary part ...  frev 2009 JYU FH
    set<int> labreq;
    if(nreqedgeslab && verbosity) cout << " label of required edges " ;
    for (int i=0; i <nreqedgeslab;++i)
      {
	  if(verbosity)
	      cout << " " << reqedgeslab[i];
	  labreq.insert(reqedgeslab[i]);
      }
    bamg::GeometricalEdge paszero;  // add JYU    fevr 2009   for  required edge ....
    if(nreqedgeslab && verbosity) cout << endl;
    int k=0;
    for (i = 0; i < Th.neb; i++)
    {
      Tn->edges[i].v[0] = Tn->vertices + Th(Th.bedges[i][0]);
      Tn->edges[i].v[1] = Tn->vertices + Th(Th.bedges[i][1]);
      Tn->edges[i].ref = Th.bedges[i].lab;
      Tn->edges[i].on = 0;
      if( labreq.find( Tn->edges[i].ref) != labreq.end())
	{
	    k++;
	    Tn->edges[i].on = &paszero;
	}

    }
  if(verbosity)cout << "  number of required edges : "<< k << endl;


  Tn->ConsGeometry(cutoffradian);
  Tn->Gh.AfterRead();
  Tn->SetIntCoor();
  Tn->FillHoleInMesh();
  return Tn;
}


 bamg::Triangles * msh2bamg(const Fem2D::Mesh & Th,double cutoffradian,
                           int  nbdfv, int * ndfv,int  nbdfe, int * ndfe,
			   long * reqedgeslab,int nreqedgeslab)
{
  using namespace bamg;
  Triangles *Tn=new Triangles(Th.nv);
  KN<int> equiedges(Th.neb);
  for(int i=0;i<Th.neb;i++)
    equiedges[i]=2*i;
  if(nbdfe !=0 )
  {
  KN<int>  kk(Th.neb),kn(Th.neb);
   kk=0;
    for(int i=0;i<Th.neb;i++)
      {
         int df=ndfe[i];
         kk[df]++;
         if(kk[df]==1) kn[df]=i;
         else {
           int k=kn[df],sens=0;
           int di0=ndfv[Th(Th.bedges[i][0])];
           int di1=ndfv[Th(Th.bedges[i][1])];
           int dk0=ndfv[Th(Th.bedges[k][0])];
           int dk1=ndfv[Th(Th.bedges[k][1])];
           if ((di0==dk0) &&(di1==dk1) ) sens=0;
           else if ((di1==dk0) &&(di0==dk1) ) sens=1;
           else  {
              cout << "Error in periodic mesh " << di0 << " " << di1 << " <=> " << dk0 << " " << dk1 << endl;
              ExecError("bug periodic mesh in ??? ");
           }
           equiedges[i]=2*k+sens;

         }
      }

  }; // a faire pour les maillages periodique

  Tn->nbv = Th.nv;
  Tn->nbt = Th.nt;
  Tn->nbe = Th.neb;
  Tn->name= new char[strlen("msh2bamg")+1];
  strcpy(Tn->name,"msh2bamg");
  throwassert(Tn->triangles);
  Tn->edges = new Edge [Th.neb];

  Int4 i;
    Metric Mid(1.);
  for (i = 0; i < Th.nv; i++)
    {
      Tn->vertices[i].r.x = Th(i).x;
      Tn->vertices[i].r.y = Th(i).y;
      Tn->vertices[i].ReferenceNumber = Th(i).lab;
      Tn->vertices[i].m=Mid;
    }

  for (i = 0; i < Th.nt; i++)
    {
      int i1 = Th(Th[i][0]);
      int i2 = Th(Th[i][1]);
      int i3 = Th(Th[i][2]);
      Tn->triangles[i]= Triangle( Tn,i1 ,i2 ,i3 );
      Tn->triangles[i].color = Th[i].lab;
    }

    // add code   un change boundary part ...  frev 2009 JYU FH
    set<int> labreq;
    if(nreqedgeslab && verbosity) cout << " label of required edges " ;
    for (int i=0; i <nreqedgeslab;++i)
      {
	  if(verbosity)
	      cout << " " << reqedgeslab[i];
	  labreq.insert(reqedgeslab[i]);
      }
    bamg::GeometricalEdge paszero;  // add JYU    fevr 2009   for  required edge ....
    if(nreqedgeslab && verbosity) cout << endl;
    int k=0;

  for (i = 0; i < Th.neb; i++)
    {
      Tn->edges[i].v[0] = Tn->vertices + Th(Th.bedges[i][0]);
      Tn->edges[i].v[1] = Tn->vertices + Th(Th.bedges[i][1]);
      Tn->edges[i].ref = Th.bedges[i].lab;
      Tn->edges[i].on = 0;
      if( labreq.find( Tn->edges[i].ref) != labreq.end())
	  {
	      k++;
	      Tn->edges[i].on = &paszero;
	  }
    }
  Tn->ConsGeometry(cutoffradian,equiedges);
  Tn->Gh.AfterRead();
  Tn->SetIntCoor();
  Tn->FillHoleInMesh();
  return Tn;
}

extern double  genrand_res53(void) ;
double NormalDistrib(double sigma)
{
    if( sigma <=0.) return 0.;
    double rand = genrand_res53()  ;
    const double TWOPI = 3.14159265358979323846264338328*2.;
    return  sigma*sqrt(-2.0*log(rand))*cos(TWOPI*rand);
}    // (stack,b,true,0,true); bool Requiredboundary=true, KNM<double> *pintern=0, double alea=0)
const Fem2D::Mesh *  BuildMesh(Stack stack, E_BorderN const * const & b,bool justboundary,int nbvmax,bool Requiredboundary,KNM<double> *pintern,double alea)
{
    if(alea) Requiredboundary=1;
    int nbvinter=0;
    if( pintern)
    {
        nbvinter=pintern->N();
        if((pintern->M() != 2 ) && ( pintern->M()!=3))
        {
            cout << " point m = " <<pintern->M()<<endl;
            ExecError("Errror: BuildMesh number of column of internal point (point=)  must be 2 or 3!");
        }
    }
    int brefintp= -2000000000;
  using namespace bamg;
  using bamg::Abs;
  using bamg::Max;
  using bamg::Min;
  using bamg::Pi;
  Fem2D::MeshPoint & mp (*Fem2D::MeshPointStack(stack)), mps = mp;

  int nbvx=nbvinter,nbe=0,nbsd=0;
  for (E_BorderN const * k=b;k;k=k->next)
  {
      int nbd = k->NbBorder(stack);
      for(int index=0; index<nbd; ++index )
    {long n=  Max(1L,Abs(k->Nbseg(stack,index))); ;
    nbvx += n+1;
    nbe += n;
    nbsd++;
    }
  }
  Geometry * Gh =  new Geometry;

  if(verbosity>2)
    cout <<"\t\t"  << "  Begin: ConstGeometry from nb Border  "  << nbsd <<endl;
  const char * filename = "FREEFEM.gh";
  Gh->name=new char [strlen(filename)+1];
  strcpy(Gh->name,filename);
  Real8 Hmin = HUGE_VAL;// the infinie value
  Int4 hvertices =0;
  Int4 i,nn,n;
  Gh->MaximalAngleOfCorner =30.00*Pi/180.0;
  Gh->nbv = 0;
  Gh->nbvx = nbvx;

  Gh->nbe = nbe;
  Gh->edges = new GeometricalEdge[Gh->nbe];
  bamg::Vertex *vertices =  new Vertex[Gh->nbvx];// correction 2009/07/03
  double lmin= HUGE_VAL;
  //  generation des points et des lignes
  i=0;
  for (E_BorderN const * k=b;k;k=k->next)
    {
    int nbd = k->NbBorder(stack);
    for(int index=0; index<nbd; ++index )
    {
      assert(k->b->xfrom); // a faire
      double & t = *  k->var(stack),tt;
      double a(k->from(stack)),b(k->to(stack));
      long * indx = k->index(stack);
      if(indx) *indx = index;
      else ffassert(index==0);
      n=Max(Abs(k->Nbseg(stack,index)),1L);
      tt=t=a;
      double delta = (b-a)/n;
      for ( nn=0;nn<=n;nn++,i++, tt += delta)
        {
          t = tt;
          if (nn==n) t=b; // to remove roundoff error
            if( nn >0 && nn < n) { t += NormalDistrib(alea); } // Add F. Hecht Juin 2018 for J-M Sac Epee:  jean-marc.sac-epee@univ-lorraine.fr
          mp.label = k->label();
          k->code(stack); // compute x,y, label
          vertices[i].r.x=mp.P.x;
          vertices[i].r.y=mp.P.y;
          vertices[i].ReferenceNumber=  mp.label;
          vertices[i].color = i;
          if (nn>0) {
            lmin=min(lmin,Norme2_2( vertices[i].r-vertices[i-1].r));
          }
        }
    }
}
 // add interna point
    if(pintern)
    {
        for(int k=0; k<nbvinter; ++k)
        {
            ffassert(i <= nbvx);
            vertices[i].r.x=(*pintern)((long) k,0L);
            vertices[i].r.y=(*pintern)((long) k,1L);
            vertices[i].ReferenceNumber=  brefintp + k ;// code internal point ..
            vertices[i].color = i;
            i++;
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

  double diameter=Max(Gh->pmax.x-Gh->pmin.x,Gh->pmax.y-Gh->pmin.y);
  Gh->coefIcoor= (MaxICoor)/diameter;
  Icoor1 epsI = (Icoor1) (Gh->coefIcoor*eps);
  ffassert(Gh->coefIcoor >0);

 if(lmin<diameter*1e-7) {
    ExecError(" Error points  border points to close < diameter*1e-7 ");}

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
  } // to delete quadtree
  if (verbosity>1)
  cout << " Nb of common points " << nbvprev-nbv <<endl;

  Gh->nbvx = nbv;
  Gh->nbv = nbv;

  Gh->vertices = new GeometricalVertex[nbv];
  throwassert(Gh->nbvx >= Gh->nbv);
  Gh->nbiv = Gh->nbv;
  const Direction NoDirOfSearch;
  //  compression of points
  int kkk;
  for ( i=0,kkk=0;kkk<nbvprev;kkk++)
    {
      if (vertices[kkk].color == i) // if new points
        {
          Gh->vertices[i].r.x = vertices[kkk].r.x ;
          Gh->vertices[i].r.y = vertices[kkk].r.y;
          throwassert(Gh->vertices[i].IsThe());
          Gh->vertices[i].ReferenceNumber = vertices[kkk].ReferenceNumber  ;
          Gh->vertices[i].DirOfSearch = NoDirOfSearch;
          Gh->vertices[i].color =0;
          Gh->vertices[i].Set();
          if(Requiredboundary)
           Gh->vertices[i].SetRequired();

          if(Gh->vertices[i].ReferenceNumber < 0)
            Gh->vertices[i].SetRequired();

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

  {    int nbd = k->NbBorder(stack);
      for(int index=0; index<nbd; ++index )
      {
      double & t = *  k->var(stack);
      double a(k->from(stack)),b(k->to(stack));
      n=Max(Abs(k->Nbseg(stack,index)),1L);
      long * indx = (k->index(stack));
      if(indx) *indx = index;
      else ffassert(index==0);

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
	  if(Requiredboundary)
	  Gh->edges[i].SetRequired();

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
      }}

  delete [] vertices; vertices=0;

  throwassert(nnn+nbvinter==nbvprev);
  throwassert(i==Gh->nbe);
  // definition  the default of the given mesh size
  if (!hvertices)
    {
        bool hvint = pintern ? pintern->M() ==3 : 0;
      for (i=0;i<Gh->nbv;i++)
      {
        if(hvint &&Gh->vertices[i].ReferenceNumber <brefintp +nbvinter)
        {
            long k =Gh->vertices[i].ReferenceNumber-brefintp;
            Gh->vertices[i].m=Metric( (*pintern)(k ,2L));
            Gh->vertices[i].ReferenceNumber = -1; //++ bof bof FH ..

        }
        else if (Gh->vertices[i].color > 0)
          Gh->vertices[i].m=  Metric(len[i] /(Real4) Gh->vertices[i].color);
        else
          Gh->vertices[i].m=  Metric(Hmin);
      }
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
        for (E_BorderN const * k=b;k;k=k->next)
        {
            int nbd = k->NbBorder(stack);
            for(int index=0; index<nbd; ++index,i++)
            {
                long Nbseg =k->Nbseg(stack,index);
                long n=  Max(1L,Abs(Nbseg));
                Gh->subdomains[i].sens = Nbseg >0 ? 1 : -1;
                Gh->subdomains[i].edge=Gh->edges + i1;
                Gh->subdomains[i].ref = i;
                i1 += n;
            }}
    }
  Gh->NbEquiEdges=0;
  Gh->NbCrackedEdges=0;
  const Fem2D::Mesh * m=0;
  if (justboundary)
    m=bamg2msh(*Gh);
  else
  {
      Gh->AfterRead();
      int nbtx= nbvmax ? nbvmax :  (Gh->nbv*Gh->nbv)/9 +1000;
      if(verbosity> 99) cout << " ** Gh = " << endl << *Gh << endl << " *** " <<endl; ;
      Triangles *Th = 0;
      try {
	  Th =new Triangles( nbtx ,*Gh);
          if(alea) //  Add F. Hecht Juin 2018 for J-M Sac Epee:  jean-marc.sac-epee@univ-lorraine.fr
          {
              Th->SetVertexFieldOn();
              for( int i=0;i<Th->nbv;++i)
              {
                  VertexOnGeom *on=0;
                  if( !(Th->vertices[i].on) ) // we are non on geometry
                  {
                      // move a little the  points
                      Th->vertices[i].r.x += NormalDistrib(alea);
                      Th->vertices[i].r.y += NormalDistrib(alea);

                  }
              }

          }
	  if(0)
	    {


	      Th->SetVertexFieldOn();
	      for( int i=0;i<Th->nbv;++i)
		{
		  VertexOnGeom *on=0;
		  if( (on =Th->vertices[i].on) ) // we are on geometrie
		    {
		      if(on->abscisse <0) {
			  bamg::GeometricalVertex * gv= on->gv;
		      }
		      else {// erreur car un point est sur un arete en non un sommet
			  bamg::GeometricalEdge * ge= on->ge;
		      }
		    }
		}
	    }
          m=bamg2msh(Th,true);

      }
      catch(...)
      {
	  Gh->NbRef=0;
	  delete Gh;
          if(m) delete m;
          if(Th) delete Th;// clean memory ???
	  cout << " catch Err bamg "  << endl;
	  throw ;
     }
      delete Th;
  }

  delete Gh;
  /* deja fait  dans bamg2msh
     Fem2D::R2 Pn,Px;
     m->BoundingBox(Pn,Px);
     m->quadtree=new Fem2D::FQuadTree(m,Pn,Px,m->nv);
  ---------- */
  mp=mps;
  return m;
}

void E_BorderN::BoundingBox(Stack stack,double  &xmin,double & xmax, double & ymin,double & ymax, double & zmin,double & zmax) const
{
  Fem2D::MeshPoint & mp (*Fem2D::MeshPointStack(stack)), mps = mp;
    
  for (E_BorderN const * k=this;k;k=k->next)
    {
        int nbd = k->NbBorder(stack);
        for(int index=0; index<nbd; ++index )
        {
      assert(k->b->xfrom); // a faire
      double & t = *  k->var(stack);
      double a(k->from(stack)),b(k->to(stack));
      long * indx = (k->index(stack));
      if(indx) *indx = index;
      else ffassert(index==0);

      long n=Max(Abs(k->Nbseg(stack,index)),1L);
      t=a;
      double delta = (b-a)/n;
      for (int  nn=0;nn<=n;nn++, t += delta)
        {
          if (nn==n) t=b; // to remove roundoff error
          mp.P=Fem2D::R3();// reset , x,y, z to 0 . Jan 2020 FH.
          mp.label = k->label();
          k->code(stack); // compute x,y, label
          xmin=Min(xmin,mp.P.x);
          xmax=Max(xmax,mp.P.x);
          ymin=Min(ymin,mp.P.y);
          ymax=Max(ymax,mp.P.y);
          zmin=Min(zmin,mp.P.z);
          zmax=Max(zmax,mp.P.z);
        }
        }}
  mp=mps;
}


void E_BorderN::Plot(Stack stack) const
{
  using Fem2D::R2;

  Fem2D::MeshPoint & mp (*Fem2D::MeshPointStack(stack)), mps = mp;
  float x0,x1,y0,y1;
  getcadre(x0,x1,y0,y1);
  float h= (x1-x0)*0.01;
  int nbd=0;
  for (E_BorderN const * k=this;k;k=k->next)
    {
      int nbdr = k->NbBorder(stack);
      for(int index=0; index<nbdr; ++index )
     {
      nbd++;
      assert(k->b->xfrom); // a faire
      double & t = *  k->var(stack);
      double a(k->from(stack)),b(k->to(stack));
      long n=Max(Abs(k->Nbseg(stack,index)),1L);
      long * indx = (k->index(stack));
      if(indx) *indx = index;
      else ffassert(index==0);

      t=a;
      double delta = (b-a)/n;
      R2 P,Po;
      for (int  nn=0;nn<=n;nn++, t += delta)
        {
          if (nn==n) t=b; // to remove roundoff error
          mp.label = k->label();
          mp.P.z=0;
          k->code(stack); // compute x,y, label
          P=mp.P.p2();
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
            DrawMark(mp.P.p2(),0.01);
            MoveTo(mp.P.p2());

          }
          Po=P;
        }
      DrawMark(mp.P.p2(),0.01);
      MoveTo(mp.P.p2());
    }
}
  if(verbosity>9) cout << "  -- Plot size : " << nbd << " Border \n";
  mp=mps;
}
void E_BorderN::SavePlot(Stack stack,PlotStream & plot ) const
{
    using Fem2D::R2;

    Fem2D::MeshPoint & mp (*Fem2D::MeshPointStack(stack)), mps = mp;

    long nbd1=0;// nb of sub border
    for (E_BorderN const * k=this;k;k=k->next)
    {
     int nbdr = k->NbBorder(stack);
     for(int index=0; index<nbdr; ++index )
	nbd1++;
    }
   plot << nbd1;
   int nbd=0;
    for (E_BorderN const * k=this;k;k=k->next)
      {
          int nbdr = k->NbBorder(stack);
          for(int index=0; index<nbdr; ++index )
          {
	  nbd++;
	  assert(k->b->xfrom); // a faire
	  double & t = *  k->var(stack);
	  double a(k->from(stack)),b(k->to(stack));
	  long n=Max(Abs(k->Nbseg(stack,index)),1L);
          long * indx = (k->index(stack));
          if(indx) *indx = index;
          else ffassert(index==0);

	  t=a;
	  double delta = (b-a)/n;
	  R2 P,Po;
	  plot<< (long) n;
	  for (int  nn=0;nn<=n;nn++, t += delta)
	    {
		if (nn==n) t=b; // to remove roundoff error
		mp.label = k->label();
		mp.P.z=0;
		k->code(stack); // compute x,y, label
		P=mp.P.p2();
		plot << (long) mp.label <<P.x << P.y;
	    }

          }}
    assert(nbd==nbd1);
    if(verbosity>9) cout << "  -- Plot size : " << nbd << " Border \n";
    mp=mps;
}

const Fem2D::Mesh *  BuildMeshBorder(Stack stack, E_BorderN const * const & b)
{
  return BuildMesh(stack,b,true,0,true);
}
const Fem2D::Mesh *  BuildMesh(Stack stack, E_BorderN const * const & b,bool Requiredboundary)
{
  return BuildMesh(stack,b,false,0,Requiredboundary);
}

const Fem2D::Mesh *  ReadTriangulate( string  * const & s) {
  using namespace Fem2D;
  KN<R2> xy;
  char c;
  int nv=0;
  for(int step=0;step<2;step++)
    {
      nv=0;
      ifstream f(s->c_str());
      if(!f) {cerr <<" Error opening file " << *s << endl;
      ExecError("Openning file ");}
      while (f.good())
        {
          R2 P;
          f >> P ;
          if (!f.good()) break;
          if (step) xy[nv]=P;
          nv++;
          while (f.get(c) &&  (c!='\n' && c!='\r' ) ) (void) 0; // eat until control (new line
        }
      if (!step && nv ) xy.init(nv); // alloc the array
    }
  if(verbosity)
    cout << " we read  " << nv << " coordinates  xy "<< endl;

  Mesh * m=new Mesh(nv,xy);
  m->MakeQuadTree();
  return m;

}
const Fem2D::Mesh *  Triangulate( const  KN_<double> & xx,const  KN_<double> & yy)
{
    using namespace Fem2D;
    ffassert(xx.N()==yy.N());
    int nv=xx.N();
    KN<R2> xy(nv);
    for(int i=0;i<nv;++i)
       xy[i]= R2(xx[i],yy[i]);
    Mesh * m=new Mesh(nv,xy);
    m->MakeQuadTree();
    return m;

}
const Fem2D::Mesh *  ReadMeshbamg( string * const & s) {
  using bamg::Triangles;
  Triangles * bTh= new Triangles(s->c_str());
 const Fem2D::Mesh * m=bamg2msh(bTh,false);// no renum
  delete bTh;
  return m;
}

const Fem2D::Mesh *  buildmeshbamg( string * const & s, int nbvxin) {

  using bamg::Triangles;
  using bamg::Geometry;
  Geometry Gh(s->c_str());
  int nbvx = nbvxin ? nbvxin : ((Gh.nbv*Gh.nbv)/9 +1000);
  Triangles * bTh=  new Triangles(nbvx,Gh);
  const Fem2D::Mesh * m=bamg2msh(bTh,false);// no renum
  delete bTh;
  return m;
}
