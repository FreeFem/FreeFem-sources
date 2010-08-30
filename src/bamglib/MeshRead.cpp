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
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "Meshio.h"
#include "Mesh2.h"
#include "QuadTree.h"
#include "SetOfE4.h"
#ifdef __MWERKS__
#pragma optimization_level 2
//#pragma inline_depth 1
#endif

#ifdef DRAWING1
  extern bool withrgraphique ;
#endif

namespace bamg {

static const  Direction NoDirOfSearch=Direction();

void Triangles::Read(MeshIstream & f_in,int Version,Real8 cutoffradian)
{
  Real8 hmin = HUGE_VAL;// the infinie value 
  //  Real8 MaximalAngleOfCorner = 10*Pi/180;// 
  Int4 i;
  Int4 dim=0;
  Int4 hvertices =0;
  Int4 ifgeom=0;
  Metric M1(1);
  if (verbosity>1)
    cout << "  -- ReadMesh " << f_in.CurrentFile  << " Version = " << Version << endl;
  int field=0;
  int showfield=0;
  while (f_in.cm()) 
    {  
      field=0;
      char  fieldname[256];
      if(f_in.eof()) break;
      f_in.cm() >> fieldname ;;
      if(f_in.eof()) break;
      f_in.err() ;
      if(verbosity>9)
	cout <<  "    " << fieldname << endl;
      if (!strcmp(fieldname,"MeshVersionFormatted") )
         f_in >>   Version ;
      else if(!strcmp(fieldname,"End"))
         break;
      else if (!strcmp(fieldname,"Dimension"))
         {
           f_in >>   dim ;
           assert(dim ==2);
         }
      else if  (!strcmp(fieldname,"Geometry"))
	{ 
	  char * fgeom ;
	  f_in >> fgeom;
	  //	  if (cutoffradian>=0) => bug if change edit the file geometry 
	  //  Gh.MaximalAngleOfCorner = cutoffradian;
	  if (strlen(fgeom))
	      Gh.ReadGeometry(fgeom);
	  else 
	    { 
	      // include geometry 
	      f_in.cm();
	      Gh.ReadGeometry(f_in,fgeom);
	    }

	  Gh.AfterRead();
	  ifgeom=1;
	  delete [] fgeom;
	}
      else if (!strcmp(fieldname,"Identifier"))
	{
	  if (identity) delete [] identity;
	  f_in >> identity;
	}
      else if (!strcmp(fieldname,"hVertices"))
         { hvertices =1;
	   Real4 h;
           for (i=0;i< nbv;i++) {
            f_in >>  h ;
	    vertices[i].m = Metric(h);}
         }
      else if (!strcmp(fieldname,"Vertices"))
         { 
           assert(dim ==2);
           f_in >>   nbv ;
	   if(verbosity>3)
	     cout << "   Nb of Vertices = " << nbv << endl;
	   nbvx=nbv;
	   vertices=new Vertex[nbvx];
	   assert(vertices);
	   ordre=new (Vertex* [nbvx]);
	   assert(ordre);
	   
	   nbiv = nbv;
           for (i=0;i<nbv;i++) {
             f_in >> 
	       	 	 vertices[i].r.x >>  
	           vertices[i].r.y >> 
             vertices[i].ReferenceNumber  ;
             vertices[i].DirOfSearch =NoDirOfSearch;
	           vertices[i].m=M1;
	     vertices[i].color =0;}
	   nbtx =  2*nbv-2; // for filling The Holes and quadrilaterals 
	   triangles =new Triangle[nbtx];
	   assert(triangles);
	   nbt =0;
         }
      else if (!strcmp(fieldname,"Triangles"))
	{ 
	  if (dim !=2)
	    cerr << "ReadMesh: Dimension <> 2" <<endl , f_in.ShowIoErr(0);
	  if(! vertices || !triangles  || !nbv )
	    cerr << "ReadMesh:Triangles before Vertices" <<endl,
	      f_in.ShowIoErr(0);
	  int NbOfTria;
	  f_in >>  NbOfTria ;
	  if(verbosity>3)
	    cout << "   NbOfTria = " << NbOfTria << endl;
	  if (nbt+NbOfTria >= nbtx)
	    cerr << "ReadMesh: We must have 2*NbOfQuad + NbOfTria  = "
		 << nbt + NbOfTria<<"  < 2*nbv-2 ="  << nbtx << endl,
	      f_in.ShowIoErr(0); 
	  //	  begintria = nbt;
	  for (i=0;i<NbOfTria;i++)  
	    {
	      Int4 i1,i2,i3,iref;
	      Triangle & t = triangles[nbt++];
	      f_in >>  i1 >>  i2 >>   i3 >>   iref ;
	      t = Triangle(this,i1-1,i2-1,i3-1);
	      t.color=iref;
	    }
	  // endtria=nbt;
	}
      else if (!strcmp(fieldname,"Quadrilaterals"))
        { 

	  if (dim !=2)
	    cerr << "ReadMesh: Dimension <> 2" <<endl , f_in.ShowIoErr(0);
	  if(! vertices || !triangles  || !nbv )
	    cerr << "ReadMesh:Quadrilaterals before Vertices" <<endl,
	      f_in.ShowIoErr(0);
	  f_in >>   NbOfQuad ;
	  if(verbosity>3)
	    cout << "   NbOfQuad= " << NbOfQuad << endl;
	  if (nbt+2*NbOfQuad >= nbtx)
	    cerr << "ReadMesh: We must have 2*NbOfQuad + NbOfTria  = "
		 << nbt + 2*NbOfQuad <<"  < 2*nbv-2 ="  << nbtx << endl,
	      f_in.ShowIoErr(0);
	  //  beginquad=nbt;
	  for (i=0;i<NbOfQuad;i++)  
	    { 
	      Int4 i1,i2,i3,i4,iref;
	      Triangle & t1 = triangles[nbt++];
	      Triangle & t2 = triangles[nbt++];
	      f_in >>  i1 >>  i2 >>   i3 >> i4 >>    iref ;
	      t1 = Triangle(this,i1-1,i2-1,i3-1);
	      t1.color=iref;
	      t2 = Triangle(this,i3-1,i4-1,i1-1);
	      t2.color=iref;
	      t1.SetHidden(OppositeEdge[1]); // two time  because the adj was not created 
	      t2.SetHidden(OppositeEdge[1]); 
	    }
	  // endquad=nbt;
	}
      else if (!strcmp(fieldname,"VertexOnGeometricVertex"))
	{ 
           f_in  >> NbVerticesOnGeomVertex ;
 	   if(verbosity>5)
	     cout << "   NbVerticesOnGeomVertex = "   << NbVerticesOnGeomVertex << endl
		  << " Gh.vertices " << Gh.vertices << endl;
	   if( NbVerticesOnGeomVertex) 
	     {
	       VerticesOnGeomVertex= new  VertexOnGeom[NbVerticesOnGeomVertex] ;
	       if(verbosity>7)
		 cout << "   VerticesOnGeomVertex = " << VerticesOnGeomVertex << endl
		      << "   Gh.vertices " << Gh.vertices << endl;
	       assert(VerticesOnGeomVertex);
	       for (Int4 i0=0;i0<NbVerticesOnGeomVertex;i0++)  
		 { 
		   Int4  i1,i2;
		   //VertexOnGeom & v =VerticesOnGeomVertex[i0];
		   f_in >>  i1 >> i2 ;
		   VerticesOnGeomVertex[i0]=VertexOnGeom(vertices[i1-1],Gh.vertices[i2-1]);
		 }
	     }
	 }
      else if (!strcmp(fieldname,"VertexOnGeometricEdge"))
         { 
           f_in >>  NbVerticesOnGeomEdge ;
 	   if(verbosity>3)
	     cout << "   NbVerticesOnGeomEdge = " << NbVerticesOnGeomEdge << endl;
	   if(NbVerticesOnGeomEdge) 
	     {
	       VerticesOnGeomEdge= new  VertexOnGeom[NbVerticesOnGeomEdge] ;
	       assert(VerticesOnGeomEdge);
	       for (Int4 i0=0;i0<NbVerticesOnGeomEdge;i0++)  
		 { 
		   Int4  i1,i2;
		   Real8 s;
		   //VertexOnGeom & v =VerticesOnGeomVertex[i0];
		   f_in >>  i1 >> i2  >> s;
		   VerticesOnGeomEdge[i0]=VertexOnGeom(vertices[i1-1],Gh.edges[i2-1],s);
		 }
	     }
         }
      else if (!strcmp(fieldname,"Edges"))
	{ 
          Int4 i,j, i1,i2;
          f_in >>  nbe ;
          edges = new Edge[nbe];
	  if(verbosity>5)
	    cout << "     Record Edges: Nb of Edge " << nbe << " edges " <<  edges << endl;
          assert(edges);
	  Real4 *len =0;
          if (!hvertices) 
	    {
	      len = new Real4[nbv];
	      for(i=0;i<nbv;i++)
		len[i]=0;
	    }
          for (i=0;i<nbe;i++) 
	    {
	      f_in  >> i1  >> i2  >> edges[i].ref ;
                 
	      assert(i1 >0 && i2 >0);
	      assert(i1 <= nbv && i2 <= nbv);
	      i1--;
	      i2--;
	      edges[i].v[0]= vertices +i1;
	      edges[i].v[1]= vertices +i2;
	      edges[i].adj[0]=0;
	      edges[i].adj[1]=0;

	      R2 x12 = vertices[i2].r-vertices[i1].r;
	      Real8 l12=sqrt( (x12,x12));        
	      if (!hvertices) {
		vertices[i1].color++;
		vertices[i2].color++;
		len[i1]+= l12;
		len[i2] += l12;}
	      hmin = Min(hmin,l12);
	    }
	  // definition  the default of the given mesh size 
          if (!hvertices) 
	    {
            for (i=0;i<nbv;i++) 
	      if (vertices[i].color > 0) 
		vertices[i].m=  Metric(len[i] /(Real4) vertices[i].color);
	      else 
		vertices[i].m=  Metric(hmin);
	      delete [] len;
	    }
	  if(verbosity>5)
	    cout << "     hmin " << hmin << endl;
	  // construction of edges[].adj 
	    for (i=0;i<nbv;i++) 
	      vertices[i].color = (vertices[i].color ==2) ? -1 : -2;
	    for (i=0;i<nbe;i++)
	      for (j=0;j<2;j++) 
		{ 
		  Vertex *v=edges[i].v[j];
		  Int4 i0=v->color,j0;
		  if(i0==-1)
		    v->color=i*2+j;
		  else if (i0>=0) {// i and i0 edge are adjacent by the vertex v
		    j0 =  i0%2;
		    i0 =  i0/2;
		    assert( v ==  edges[i0 ].v[j0]);
		    edges[i ].adj[ j ] =edges +i0;
		    edges[i0].adj[ j0] =edges +i ;
		    assert(edges[i0].v[j0] == v);
		    //	    if(verbosity>8)
		    //  cout << " edges adj " << i0 << " "<< j0 << " <-->  "  << i << " " << j << endl;
		      v->color = -3;}
		}

	}
	
/*   ne peut pas marche car il n'y a pas de flag dans les aretes du maillages
       else if (!strcmp(fieldname,"RequiredEdges"))
        { 
          int i,j,n;
          f_in  >> n ;

          for (i=0;i<n;i++) {     
            f_in  >>  j ;
            assert( j <= nbe );
            assert( j > 0 );
            j--;
            edges[j].SetRequired();  }
      }
*/	
      else if (!strcmp(fieldname,"EdgeOnGeometricEdge"))
	{ 
	  assert(edges);
	  int i1,i2,i,j;
	  f_in >> i2 ;
	  if(verbosity>3)
	    cout << "     Record EdgeOnGeometricEdge: Nb " << i2 <<endl;
	  for (i1=0;i1<i2;i1++)
	    {
	      f_in  >> i >> j ;
	      if(!(i>0 && j >0 && i <= nbe && j <= Gh.nbe))
		{
		  cerr << "      Record EdgeOnGeometricEdge i=" << i << " j = " << j;
		  cerr << " nbe = " << nbe << " Gh.nbe = " <<  Gh.nbe << endl;
		  cerr << " We must have : (i>0 && j >0 && i <= nbe && j <= Gh.nbe) ";
		  cerr << " Fatal error in file " << name << " line " << f_in.LineNumber << endl;
		  MeshError(1);
		}
		
		
	      edges[i-1].on = Gh.edges + j-1;
	    }
	  //	  cout << "end EdgeOnGeometricEdge" << endl;
	}
      else if (!strcmp(fieldname,"SubDomain") || !strcmp(fieldname,"SubDomainFromMesh") )
	{ 
	  
	  f_in >> NbSubDomains ;
	  subdomains = new SubDomain [ NbSubDomains ];
	    for (i=0;i< NbSubDomains;i++)
	      { Int4 i3,head,sens;
	      f_in  >>  i3 >>	head >> sens  >> subdomains[i].ref ;
		assert (i3==3);
		head --;
		assert(head < nbt && head >=0);
		subdomains[i].head = triangles+head;		
	      }
	}
      else
	{ // unkown field
	  field = ++showfield;
	  if(showfield==1) // just to show one time 
	    if (verbosity>5)
	      cout << "     Warning we skip the field " << fieldname << " at line " << f_in.LineNumber << endl;
	}
      showfield=field; // just to show one time 
    }
    
    if (ifgeom==0)
      {
	if (verbosity)
	  cout  << " ## Warning: Mesh without geometry we construct a geometry (theta =" 
		<< cutoffradian*180/Pi << " degres )" << endl;
	ConsGeometry(cutoffradian);	
	Gh.AfterRead();
      }	  
}




void Triangles::Read_am_fmt(MeshIstream & f_in)
{
  Metric M1(1);

  if (verbosity>1)
    cout << "  -- ReadMesh .am_fmt file " <<  f_in.CurrentFile  << endl;
   	
     Int4 i;
     f_in.cm() >> nbv >> nbt ;
     if (verbosity>3)
       cout << "    nbv = " << nbv  << " nbt = " << nbt << endl;
     f_in.eol() ;// 
     nbvx = nbv;
     nbtx =  2*nbv-2; // for filling The Holes and quadrilaterals 
     triangles =new Triangle[nbtx];
     assert(triangles);
     vertices=new Vertex[nbvx];
     ordre=new (Vertex* [nbvx]);
     
     for (     i=0;i<nbt;i++)
       {
	 Int4 i1,i2,i3;
	 f_in >>  i1 >>  i2 >>   i3 ;
	 triangles[i]  = Triangle(this,i1-1,i2-1,i3-1);	 
       }
      f_in.eol() ;// 
      
      for ( i=0;i<nbv;i++)
	f_in >> vertices[i].r.x >>   vertices[i].r.y,
	  vertices[i].m = M1,vertices[i].DirOfSearch =NoDirOfSearch;

      f_in.eol() ;// 
      
      for ( i=0;i<nbt;i++)  
	f_in >> triangles[i].color;
      f_in.eol() ;// 
      
      for ( i=0;i<nbv;i++)  
	     f_in >> vertices[i].ReferenceNumber;
      
      
}

////////////////////////

void  Triangles::Read_am(MeshIstream &ff)
{
  if (verbosity>1)
    cout << "  -- ReadMesh .am_fmt file " <<  ff.CurrentFile  << endl;
    Metric M1(1);	
  
  IFortranUnFormattedFile f_in(ff);

  Int4  l=f_in.Record();
  assert(l==2*sizeof(Int4));
  f_in >> nbv >> nbt ;
  l=f_in.Record();
  assert((size_t) l==nbt*sizeof(long)*4 + nbv*(2*sizeof(float)+sizeof(long)));
  if (verbosity>3)
    cout << "    nbv = " << nbv  << " nbt = " << nbt << endl;
  
  nbvx = nbv;
  nbtx =  2*nbv-2; // for filling The Holes and quadrilaterals 
  triangles =new Triangle[nbtx];
  assert(triangles);
  vertices=new Vertex[nbvx];
  ordre=new (Vertex* [nbvx]);
  

  Int4 i;
  for (     i=0;i<nbt;i++) {
    long i1,i2,i3;
    f_in >>  i1 >>  i2 >>   i3 ;
    triangles[i]  = Triangle(this,i1-1,i2-1,i3-1); }
  
  for ( i=0;i<nbv;i++) {
    float x,y;
    f_in >> x >> y;
    vertices[i].r.x =x;
    vertices[i].r.y=y;
    vertices[i].m=M1;}
  
  for ( i=0;i<nbt;i++) {
    long i;
    f_in >> i;
    triangles[i].color=i;}
  
  for ( i=0;i<nbv;i++) {
    long i;
    f_in >> i;
    vertices[i].ReferenceNumber=i;}
}

//////////////////////////////////

void  Triangles::Read_nopo(MeshIstream & ff)
{

 if (verbosity>1)
    cout << "  -- ReadMesh .nopo file " <<  ff.CurrentFile  << endl;
 IFortranUnFormattedFile f_in(ff);
 
 
 Int4  l,i,j;
 l=f_in.Record(); 
 l=f_in.Record();
 f_in >> i;
 assert(i==32);
 Int4 niveau,netat,ntacm;

 char titre[80+1],  date[2*4+1], nomcre[6*4+1], typesd[5];
 f_in.read4(titre,20);
 f_in.read4(date,2);
 f_in.read4(nomcre,6);
 f_in.read4(typesd,1);


   f_in >> niveau>>netat>>ntacm;
 if (strcmp("NOPO",typesd))
   {
     cout << " where in record  " << f_in.where() << " " << strcmp("NOPO",typesd) << endl;
     cerr << " not a  nopo file but `" << typesd <<"`"<< " len = " << strlen(typesd) << endl;
     cerr << (int) typesd[0] << (int) typesd[1] << (int) typesd[2] << (int) typesd[3] << (int) typesd[4] << endl;
     cout << " nomcre :" << nomcre << endl;
     cout << " date   :" << date << endl;
     cout << " titre  :" << titre<< endl;
     MeshError(112);
   }
 if(verbosity>2)
   cout << "    nb de tableau associe : " << ntacm << " niveau =" << niveau << endl;
 
 for (i=0;i<ntacm;i++)
   f_in.Record();

 f_in.Record();
 f_in >> l;
 assert(l == 27);
 Int4 nop2[27];
 for (i=0;i<27;i++)
   f_in >> nop2[i];
 Int4 ndim = nop2[0];
 Int4 ncopnp = nop2[3];
 Int4 ne = nop2[4];
 Int4 ntria = nop2[7];
 Int4 nquad = nop2[8];
 Int4 np = nop2[21];
 // Int4 nef = nop2[13];
 Metric M1(1);
 if(verbosity>2) 
   cout << "    ndim = " << ndim << " ncopnp= " << ncopnp << " ne = " << ne 
	<< "    ntri = " << ntria << " nquad = " << nquad << " np = " << np << endl;
 nbv = np;
 nbt = 2*nquad + ntria;
 if (ne != nquad+ntria || ndim !=2 || ncopnp != 1 )
   {
     cerr << " not only tria & quad in nopo mesh on dim != 2 ou ncopnp != 1 " << endl;
     MeshError(113);
   }
 if( nop2[24]>=0)  f_in.Record();
  NbOfQuad = nquad;
  nbvx = nbv;
  nbtx =  2*nbv-2; // for filling The Holes and quadrilaterals 
  triangles =new Triangle[nbtx];
  assert(triangles);
  vertices=new Vertex[nbvx];
  ordre=new (Vertex* [nbvx]);


 f_in >> l;

  if(verbosity>9)
    cout << " Read cnop4 nb of float  " << l << endl;
  
  assert(l==2*np);
  for (i=0;i<np;i++)
    {  float x,y;
    f_in >>  x>> y;
    vertices[i].r.x=x;
    vertices[i].r.y=y;
    vertices[i].m=M1;
    vertices[i].DirOfSearch =NoDirOfSearch;

    }
  f_in.Record();
  // lecture de nop5 bonjour les degats
  f_in >> l;
  if(verbosity>9)
    cout << " Read nop5  nb of int4 " << l << endl;
 Int4 k=0;
 Int4 nbe4 =  3*ntria + 4*nquad;
 // cout << " nbv = " << nbv << " nbe4 " << nbe4 << endl;
 SetOfEdges4 * edge4= new SetOfEdges4(nbe4,nbv); 
 Int4 * refe = new Int4[nbe4];
 Int4 kr =0;
 for (i=0;i<ne;i++)
   {
     // Int4 ng[4]={0,0,0,0};
     Int4 np[4],rv[4],re[4];
     Int4 ncge,nmae,ndsde,npo;
     f_in >> ncge >> nmae >> ndsde >> npo ;
     //cout << " element " << i << " " << ncge << " "  
     // << nmae <<" " <<  npo << endl;
     if (ncge != 3 && ncge != 4)
       {
	 cerr << " read nopo type element[" << i << "] =" 
	      << ncge << " not 3 or 4 " << endl;
	 MeshError(115);
       }
     if (npo != 3 && npo != 4)
       {
	 cerr << " read nopo element[" << i << "] npo = "  
	      << npo << " not 3 or 4 " << endl;
	 MeshError(115);
       }
     
     for( j=0;j<npo;j++)
       {f_in >>np[j];np[j]--;}
     
     if (ncopnp !=1) 
       {
	 f_in >> npo;
	 if (npo != 3 || npo != 4)
	   {
	     cerr << " read nopo type element[" << i << "]= "  
		  << ncge << " not 3 or 4 " << endl;
	     MeshError(115);
	   }
	 
	 for(j=0;j<npo;j++)
	   {f_in >>np[j];np[j]--;}
	 
       }
     if (nmae>0) 
       {
	 Int4  ining; // no ref 
	 
	 f_in>>ining;
	 if (ining==1)
	   MeshError(116);
	 if (ining==2)
	   for (j=0;j<npo;j++)
	     f_in >>re[j];
	 for (j=0;j<npo;j++)
	   f_in >>rv[j];
	 
	 
	 // set the ref on vertex and the shift np of -1 to start at 0
	 for (j=0;j<npo;j++)
	   vertices[np[j]].ReferenceNumber = rv[j];
	 
	 if (ining==2)
	   for (j=0;j<npo;j++)
	     if (re[j])
	       {
		 kr++;
		 Int4 i0 = np[j];
		 Int4 i1 = np[(j+1)%npo];
		 // cout << kr << " ref  edge " << i0 << " " << i1 << " " << re[j] << endl;
		 refe[edge4->addtrie(i0,i1)]=re[j];
	       }
       }
     
     if (npo==3) 
       { // triangles 
	 triangles[k]  = Triangle(this,np[0],np[1],np[2]); 
	 triangles[k].color = ndsde;
	 k++;
       }
     else if (npo==4)
       { // quad 
	 Triangle & t1 = triangles[k++];
	 Triangle & t2 = triangles[k++];
	 t1 = Triangle(this,np[0],np[1],np[2]);
	 t2 = Triangle(this,np[2],np[3],np[0]);
	 t1.SetHidden(OppositeEdge[1]); // two time  because the adj was not created 
	 t2.SetHidden(OppositeEdge[1]); 
	 t1.color = ndsde;
	 t2.color = ndsde;
       }
     else
       {
	 cerr << " read nopo type element =" << npo << " not 3 or 4 " << endl;
	 MeshError(114);
       }
     
     
   }
 // cout << k << " == " << nbt << endl;
 assert(k==nbt);

 nbe = edge4->nb();
 if (nbe)
   {
     if (verbosity>7)
     cout << " Nb of ref edges = " << nbe << endl;
     if (edges)
       delete [] edges;
     edges = new Edge[nbe];
     for (i=0;i<nbe;i++)
       {
	 edges[i].v[0] = vertices + edge4->i(i);
	 edges[i].v[1] = vertices + edge4->j(i);
	 edges[i].ref = refe[i];
	 //	 cout << i << " " <<  edge4->i(i) << " " <<  edge4->j(i) << endl;
       }
      if (verbosity>7)
	cout << " Number of reference edge in the  mesh = " << nbe << endl;
   }
 delete [] refe;
 delete edge4;
}
  void  Triangles::Read_ftq(MeshIstream & f_in)
{
  //  
  if (verbosity>1)
    cout << "  -- ReadMesh .ftq file " <<  f_in.CurrentFile  << endl;
  
  Int4 i,ne,nt,nq;
  f_in.cm() >> nbv >> ne >> nt >> nq ;
  if (verbosity>3)
    cout << "    nbv = " << nbv  << " nbtra = " << nt << " nbquad = " << nq << endl;
  nbt = nt+2*nq;
  

  nbvx = nbv;
  nbtx =  2*nbv-2; // for filling The Holes and quadrilaterals 
  triangles =new Triangle[nbtx];
  assert(triangles);
  vertices=new Vertex[nbvx];
  ordre=new (Vertex* [nbvx]);
  Int4 k=0;
  
  for ( i=0;i<ne;i++) 
    {
      long ii,i1,i2,i3,i4,ref;
      f_in >>  ii;
	if (ii==3) 
	  { // triangles 
	    f_in >> i1>>  i2 >>   i3 >> ref ;
	    triangles[k]  = Triangle(this,i1-1,i2-1,i3-1); 
	    triangles[k++].color = ref;
	  }
	else if (ii==4)
	  { // quad 
	    f_in >> i1>>  i2 >>   i3 >> i4 >> ref ;
	    Triangle & t1 = triangles[k++];
	    Triangle & t2 = triangles[k++];
	    t1 = Triangle(this,i1-1,i2-1,i3-1);
	    t1.color=ref;
	    t2 = Triangle(this,i3-1,i4-1,i1-1);
	    t2.color=ref;
	    t1.SetHidden(OppositeEdge[1]); // two time  because the adj was not created 
	    t2.SetHidden(OppositeEdge[1]); 	  
	  }
	else
	  {
	    cout << " read ftq type element =" << ii << " not 3 or 4 " << endl;
	    MeshError(111);
	  }
    }
  assert(k==nbt);
  Metric M1(1);
  for ( i=0;i<nbv;i++)
    {
      f_in >> vertices[i].r.x >>   vertices[i].r.y >> vertices[i].ReferenceNumber;
       vertices[i].DirOfSearch =NoDirOfSearch;
       vertices[i].m = M1;

    }
}

///////////////////////////////////////////////

void  Triangles::Read_msh(MeshIstream &f_in)
{
    Metric M1(1.);
  if (verbosity>1)
    cout << "  -- ReadMesh .msh file " <<  f_in.CurrentFile  << endl;
   	
     Int4 i;
     f_in.cm() >> nbv >> nbt ;
     while (f_in.in.peek()==' ')
         f_in.in.get();
     if(isdigit(f_in.in.peek())) 
       f_in >> nbe;
     if (verbosity>3)
       cout << "    nbv = " << nbv  << " nbt = " << nbt << " nbe = " << nbe << endl;
     nbvx = nbv;
     nbtx =  2*nbv-2; // for filling The Holes and quadrilaterals 
     triangles =new Triangle[nbtx];
     assert(triangles);
     vertices=new Vertex[nbvx];
     ordre=new (Vertex* [nbvx]);
      edges = new Edge[nbe];
     for ( i=0;i<nbv;i++)
	{
	 f_in >> vertices[i].r.x >>   vertices[i].r.y >> vertices[i].ReferenceNumber;
	    vertices[i].on=0;
	    vertices[i].m=M1;
         //if(vertices[i].ReferenceNumber>NbRef)	NbRef=vertices[i].ReferenceNumber;  
  	}
     for (     i=0;i<nbt;i++)
       {
	 Int4 i1,i2,i3,r;
	 f_in >>  i1 >>  i2 >>   i3 >> r;
	 triangles[i]  = Triangle(this,i1-1,i2-1,i3-1);
	 triangles[i].color = r;	 
       }
     for (i=0;i<nbe;i++)
       {
	 Int4 i1,i2,r;
	 f_in >>  i1 >>  i2 >> r;
	      edges[i].v[0]= vertices +i1-1;
	      edges[i].v[1]= vertices +i2-1;
	      edges[i].adj[0]=0;
	      edges[i].adj[1]=0;
	      edges[i].ref = r;
	      edges[i].on=0;
       }
      
}

//////////////////////////////////////////////////

void  Triangles::Read_amdba(MeshIstream &f_in )
{
  Metric M1(1);
  if (verbosity>1)
    cout << "  -- ReadMesh .amdba file " <<  f_in.CurrentFile  << endl;
   	
     Int4 i;
     f_in.cm() >> nbv >> nbt ;
     //    if (verbosity>3)
       cout << "    nbv = " << nbv  << " nbt = " << nbt << endl;
     f_in.eol() ;// 
     nbvx = nbv;
     nbtx =  2*nbv-2; // for filling The Holes and quadrilaterals 
     triangles =new Triangle[nbtx];
     assert(triangles);
     vertices=new Vertex[nbvx];
     ordre=new (Vertex* [nbvx]);
     Int4 j;
      for ( i=0;i<nbv;i++)
	{
	  f_in >> j ;
	  assert( j >0 && j <= nbv);
	  j--;
	  f_in >> vertices[j].r.x >>   vertices[j].r.y >> vertices[j].ReferenceNumber;
	   vertices[j].m=M1;
	   vertices[j].DirOfSearch=NoDirOfSearch;
	}
     
      for (     i=0;i<nbt;i++)
       {
	 Int4 i1,i2,i3,ref;
	   f_in >> j ;
	   assert( j >0 && j <= nbt);
	   j--;
	   f_in >> i1 >>  i2 >>   i3 >> ref;
	 triangles[j]  = Triangle(this,i1-1,i2-1,i3-1);
	 triangles[j].color =ref;
       }
      f_in.eol() ;// 
      
  // cerr << " a faire " << endl;
  //MeshError(888);
}


Triangles::Triangles(const char * filename,Real8 cutoffradian) 
: Gh(*(new Geometry())), BTh(*this)
{ 
#ifdef DRAWING1
   if (!withrgraphique) {initgraphique();withrgraphique=true;}   
#endif

  //  Int4 beginquad=0,begintria=0;
  // Int4 endquad=0;endtria=0;
  //int type_file=0;

  int lll = strlen(filename);
  int  am_fmt = !strcmp(filename + lll - 7,".am_fmt");
  int  amdba = !strcmp(filename + lll - 6,".amdba");
  int  am = !strcmp(filename + lll - 3,".am");
  int  nopo = !strcmp(filename + lll - 5,".nopo");
  int  msh = !strcmp(filename + lll - 4,".msh");
  int  ftq = !strcmp(filename + lll - 4,".ftq");

 // cout << " Lecture type  :" << filename + lll - 7 <<":" <<am_fmt<<  endl;

  char * cname = new char[lll+1];
  strcpy(cname,filename);
  Int4 inbvx =0;
  PreInit(inbvx,cname);
  OnDisk = 1;
  //  allocGeometry = &Gh; // after Preinit ; 

  MeshIstream f_in (filename);
  
  if (f_in.IsString("MeshVersionFormatted"))
    {
      int version ;
      f_in >> version ;
      Read(f_in,version,cutoffradian);
    }
  else {     
    if (am_fmt) Read_am_fmt(f_in);
    else if (am) Read_am(f_in);
    else if (amdba) Read_amdba(f_in);
    else if (msh) Read_msh(f_in);
    else if (nopo) Read_nopo(f_in);
    else if (ftq) Read_ftq(f_in);
    else 
      { 
	cerr << " Unkown type mesh " << filename << endl;
	MeshError(2);
      }
      ConsGeometry(cutoffradian);
      Gh.AfterRead();    
  }

  SetIntCoor();
  FillHoleInMesh();
  // Make the quad ---
   
}

void Geometry::ReadGeometry(const char * filename)
{
  OnDisk = 1;
  if(verbosity>1)
      cout << "  -- ReadGeometry " << filename << endl;
  MeshIstream f_in (filename);
  ReadGeometry(f_in,filename);
}



void Geometry::ReadGeometry(MeshIstream & f_in,const char * filename)  
{
  if(verbosity>1)
    cout << "  -- ReadGeometry " << filename << endl;
  assert(empty());
  nbiv=nbv=nbvx=0;
  nbe=nbt=nbtx=0;
  NbOfCurves=0;
 // BeginOfCurve=0;
  name=new char [strlen(filename)+1];
  strcpy(name,filename);
  Real8 Hmin = HUGE_VAL;// the infinie value 
//  Real8 MaximalAngleOfCorner = 10*Pi/180; ; 
  Int4 hvertices =0;
  Int4 i;
  Int4 Version,dim=0;
  int field=0;
  int showfield=0;
  int NbErr=0;

  while (f_in.cm()) 
    { 
      field=0;
      // warning ruse for on allocate fiedname at each time 
      char fieldname[256];
      f_in.cm() >> fieldname ;
      f_in.err();
      if(f_in.eof()) break;
//      cout <<  fieldname <<  " line " << LineNumber  <<endl;
      if (!strcmp(fieldname,"MeshVersionFormatted") )
        f_in  >> Version ;
      else if (!strcmp(fieldname,"End"))
	break;
      else if (!strcmp(fieldname,"end"))
	break;
      else if (!strcmp(fieldname,"Dimension"))
        {
          f_in   >>  dim ;
	  if(verbosity>5) 
	    cout << "     Geom Record Dimension dim = " << dim << endl;        
          assert(dim ==2);
         }
       else if (!strcmp(fieldname,"hVertices"))
         { 
	   if (nbv <=0) {
	     cerr<<"Error: the field Vertex is not found before hVertices " << filename<<endl;
	     NbErr++;}       
	   if(verbosity>5) 
	    cout << "     Geom Record hVertices nbv=" << nbv <<  endl;
	   hvertices =1;
           for (i=0;i< nbv;i++) 
	     {
	       Real4 h;
	       f_in  >>  h ; 
	       vertices[i].m = Metric(h);
	     }
	 }
       else if (!strcmp(fieldname,"MetricVertices"))
         { hvertices =1;
	   if (nbv <=0) {
	     cerr<<"Error: the field Vertex is not found before MetricVertices " << filename<<endl;
	     NbErr++;}       
           if(verbosity>5) 
	     cout << "     Geom Record MetricVertices nbv =" << nbv <<  endl;
           for (i=0;i< nbv;i++) 
	     {
	       Real4 a11,a21,a22;
	       f_in  >>  a11 >> a21 >> a22  ; 
	       vertices[i].m = Metric(a11,a21,a22);
	     }
	 }
       else if (!strcmp(fieldname,"h1h2VpVertices"))
         { hvertices =1;
	   if (nbv <=0) {
	     cerr<<"Error: the field Vertex is not found before h1h2VpVertices " << filename<<endl;
	     NbErr++;}       
           if(verbosity>5) 
	     cout << "     Geom Record h1h2VpVertices nbv=" << nbv << endl;

           for (i=0;i< nbv;i++) 
	     {
	       Real4 h1,h2,v1,v2;
	       f_in  >> h1 >> h2 >>v1 >>v2 ; 
	       vertices[i].m = Metric(MatVVP2x2(1/(h1*h1),1/(h2*h2),D2(v1,v2)));
	     }
	 }
      else if (!strcmp(fieldname,"Vertices"))
        { 
          assert(dim ==2);
          f_in   >>  nbv ;
	  //          if(LineError) break;
          nbvx = nbv;
          
          vertices = new GeometricalVertex[nbvx];
	  if(verbosity>5) 
	    cout << "     Geom Record Vertices nbv = " << nbv << "vertices = " << vertices<<endl;
          assert(nbvx >= nbv);
          nbiv = nbv;
          for (i=0;i<nbv;i++) {
            f_in  >> vertices[i].r.x  ;
            // if(LineError) break;
            f_in  >> vertices[i].r.y ;
	    // if(LineError) break;
            f_in >>   vertices[i].ReferenceNumber   ;
            vertices[i].DirOfSearch=NoDirOfSearch;
	    //            vertices[i].m.h = 0;
            vertices[i].color =0;
            vertices[i].Set();}
	  // if(LineError) break;
	    pmin =  vertices[0].r;
	    pmax =  vertices[0].r;
	    // recherche des extrema des vertices pmin,pmax
	    for (i=0;i<nbv;i++) {
	      pmin.x = Min(pmin.x,vertices[i].r.x);
	      pmin.y = Min(pmin.y,vertices[i].r.y);
	      pmax.x = Max(pmax.x,vertices[i].r.x);
	      pmax.y = Max(pmax.y,vertices[i].r.y);
	    }
	    
	      R2 DD05 = (pmax-pmin)*0.05;
	      pmin -=  DD05;
	      pmax +=  DD05;
	    
	    coefIcoor= (MaxICoor)/(Max(pmax.x-pmin.x,pmax.y-pmin.y));
	    assert(coefIcoor >0);
	    if (verbosity>5) {
	      cout << "     Geom: min="<< pmin << "max ="<< pmax << " hmin = " << MinimalHmin() <<  endl;}
        }
      else if(!strcmp(fieldname,"MaximalAngleOfCorner")||!strcmp(fieldname,"AngleOfCornerBound"))
        {         
          f_in >> MaximalAngleOfCorner;
              
	  if(verbosity>5) 
	    cout << "     Geom Record MaximalAngleOfCorner " << MaximalAngleOfCorner <<endl;
          MaximalAngleOfCorner *= Pi/180;
        }
      else if (!strcmp(fieldname,"Edges"))
        {
	  if (nbv <=0) {
	    cerr<<"Error: the field edges is not found before MetricVertices " << filename<<endl;
	    NbErr++;}   
	  else 
	    {
	      int i1,i2;
	      R2 zero2(0,0);
	      f_in   >>  nbe ;
	      
	      edges = new GeometricalEdge[nbe];
	      if(verbosity>5) 
		cout << "     Record Edges: Nb of Edge " << nbe <<endl;
	      assert(edges);
	      assert (nbv >0); 
	      Real4 *len =0;
	      if (!hvertices) 
		{
		  len = new Real4[nbv];
		  for(i=0;i<nbv;i++)
		    len[i]=0;
		}
	      
	      for (i=0;i<nbe;i++) 
		{
		  f_in  >> i1   >> i2 >>  edges[i].ref  ;
		  
		  i1--;i2--; // for C index
		  edges[i].v[0]=  vertices + i1;
		  edges[i].v[1]=  vertices + i2;
		  R2 x12 = vertices[i2].r-vertices[i1].r;
		  Real8 l12=sqrt((x12,x12));
		  edges[i].tg[0]=zero2;
		  edges[i].tg[1]=zero2;
		  edges[i].SensAdj[0] = edges[i].SensAdj[1] = -1;
		  edges[i].Adj[0] = edges[i].Adj[1] = 0;
		  edges[i].flag = 0;
		  if (!hvertices) 
		    {
		      vertices[i1].color++;
		      vertices[i2].color++;
		      len[i1] += l12;
		      len[i2] += l12;
		    }
		  
		  Hmin = Min(Hmin,l12);
		}
	      // definition  the default of the given mesh size 
	      if (!hvertices) 
		{
		  for (i=0;i<nbv;i++) 
		    if (vertices[i].color > 0) 
		      vertices[i].m=  Metric(len[i] /(Real4) vertices[i].color);
		    else 
		      vertices[i].m=  Metric(Hmin);
		  delete [] len;
		  
		  if(verbosity>3) 
		    cout << "     Geom Hmin " << Hmin << endl;
		}
	      
	    }
	}
      else if (!strcmp(fieldname,"EdgesTangence") ||!strcmp(fieldname,"TangentAtEdges")  )
        { 
          int n,i,j,k;
          R2 tg;
          f_in  >> n ;
          
	  if(verbosity>5) 
	    cout << "     Record TangentAtEdges: Nb of Edge " << n <<endl;
          
          for (k=0;k<n;k++)
            {
	      f_in  >>  i  >> j ;
	      f_in >> tg.x  >> tg.y ;
	      assert( i <= nbe );
	      assert( i > 0 );
	      assert ( j == 1 || j==2 );
	      i--;j--;// for C index
	      edges[i].tg[j] = tg;
            }
        }
      else if (!strcmp(fieldname,"Corners"))
        { 
          int i,j,n;
          f_in  >> n ;
	  if(verbosity>5) 
	    cout << "     Record Corner: Nb of Corner " << n <<endl;
          
          for (i=0;i<n;i++) {     
            f_in  >>  j ;
            assert( j <= nbv );
            assert( j > 0 );
            j--;
            vertices[j].SetCorner();
            vertices[j].SetRequired();  }
        }
      else if (!strcmp(fieldname,"RequiredVertices"))
        { 
          int i,j,n;
          f_in  >> n ;

          for (i=0;i<n;i++) {     
            f_in  >>  j ;
            assert( j <= nbv );
            assert( j > 0 );
            j--;
            vertices[j].SetRequired();  }
      }
      else if (!strcmp(fieldname,"RequiredEdges"))
        { 
          int i,j,n;
          f_in  >> n ;

          for (i=0;i<n;i++) {     
            f_in  >>  j ;
            assert( j <= nbe );
            assert( j > 0 );
            j--;
            edges[j].SetRequired();  }
      }
    else if (!strcmp(fieldname,"SubDomain") || !strcmp(fieldname,"SubDomainFromGeom"))
      { 
	f_in   >>  NbSubDomains ;
	if (NbSubDomains>0) 
	  {
	    subdomains = new GeometricalSubDomain[  NbSubDomains];
	    Int4 i0,i1,i2,i3;
	    for (i=0;i<NbSubDomains;i++) 
	      {
		f_in  >> i0  >>i1 
		      >> i2  >>i3 ; 
		
		assert(i0 == 2);
		assert(i1<=nbe && i1>0);
		subdomains[i].edge=edges + (i1-1);
		subdomains[i].sens = (int) i2;
		subdomains[i].ref = i3;
	      }
	  }
      }
      else
	{ // unkown field
	  field = ++showfield;
	  if(showfield==1) // just to show one time 
	    if (verbosity>3)
	      cout << "    Warning we skip the field " << fieldname << " at line " << f_in.LineNumber << endl;
	}
      showfield=field; // just to show one time 
    } // while !eof()
  // generation  de la geometrie 
  // 1 construction des aretes 
  // construire des aretes en chaque sommets 
  
  if (nbv <=0) {
    cerr<<"Error: the field Vertex is not found in " << filename<<endl;
    NbErr++;}
  if(nbe <=0) {
    cerr <<"Error: the field Edges is not found in "<< filename<<endl
      ;NbErr++;}
  if(NbErr) MeshError(1);

 
}


}  // end of namespace bamg 
