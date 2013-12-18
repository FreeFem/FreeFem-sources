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
#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include "Meshio.h"
#include "Mesh2.h"
#include "QuadTree.h"
#include "SetOfE4.h"


namespace bamg {

void Triangles::Write(const char * filename,const TypeFileMesh typein )
{
  TypeFileMesh type = typein;
  const char * gsuffix=".gmsh";
  int ls=0;
  int lll = strlen(filename);
  if (type==AutoMesh)
    {
      type = BDMesh;
      if      (!strcmp(filename + lll - (ls=7),".am_fmt")) type = am_fmtMesh;
      else if (!strcmp(filename + lll - (ls=6),".amdba"))  type = amdbaMesh;
      else if (!strcmp(filename + lll - (ls=3),".am"))     type = amMesh;
      else if (!strcmp(filename + lll - (ls=5),".nopo"))   type = NOPOMesh;
      else if (!strcmp(filename + lll - (ls=4),".msh"))    type = mshMesh;
      else if (!strcmp(filename + lll - (ls=4),".ftq"))    type = ftqMesh;
      else if (!strcmp(filename + lll - (ls=7),".AM_FMT")) type = am_fmtMesh;
      else if (!strcmp(filename + lll - (ls=6),".AMDBA"))  type = amdbaMesh;
      else if (!strcmp(filename + lll - (ls=3),".AM"))     type = amMesh;
      else if (!strcmp(filename + lll - (ls=5),".NOPO"))   type = NOPOMesh;
      else if (!strcmp(filename + lll - (ls=4),".MSH"))    type = mshMesh;
      else if (!strcmp(filename + lll - (ls=4),".FTQ"))    type = ftqMesh;
      else ls=0;
    } 
  if (verbosity>1)
    {
      cout << "  -- Writing the file " << filename << " of type " ;
     switch (type) 
       {
       case BDMesh     :  cout << " BD Mesh "  ; break;
       case NOPOMesh   :  cout << " NOPO "     ; break;
       case amMesh     :  cout << " am "       ; break;
       case am_fmtMesh :  cout << " am_fmt "   ; break;
       case amdbaMesh  :  cout << " amdba "    ; break;
       case ftqMesh    :  cout << " ftq "      ; break;
       case mshMesh    :  cout << " msh "      ; break;
	default: 
	  cerr << endl 
	       <<  " Unkown type mesh file " << (int) type << " for Writing " << filename <<endl;
	  MeshError(1);
       }     
     Int4 NbOfTria =  nbt-2*NbOfQuad-NbOutT ;
     if (NbOfTria)      cout << " NbOfTria = " << NbOfTria;
     if (NbOfQuad)      cout << " NbOfQuad = " << NbOfQuad;
     if (nbe)   	cout << " NbOfRefEdge = " << nbe ;
     cout    << endl;

    }
  ofstream f(filename /*,ios::trunc*/);
 f.precision(12);

   if (f)
     switch (type) 
       {
       case BDMesh     : 
          {
             if ( ! Gh.OnDisk)
              {
                 delete [] Gh.name;
                 Gh.name = new char[lll+1+strlen(gsuffix)];
                 strcpy(Gh.name,filename);
                 if (Gh.name[lll-ls-1]=='.') strcpy(Gh.name+lll-ls, gsuffix+1);
                 else strcpy(Gh.name+lll-ls,gsuffix);
                cout << " write geo in " << Gh.name << endl;
                  ofstream f(Gh.name) ;
                  f << Gh ;
                  Gh.OnDisk=true;
              }
	     f << *this     ;
	     break;
           }
       case NOPOMesh   :  Write_nopo(f)  ; break;
       case amMesh     :  Write_am(f)    ; break;
       case am_fmtMesh :  Write_am_fmt(f); break;
       case amdbaMesh  :  Write_amdba(f) ; break;
       case ftqMesh    :  Write_ftq(f)   ; break;
       case mshMesh    :  Write_msh(f)   ; break;
	default: 
	  cerr << " Unkown type mesh file " << (int) type << " for Writing " << filename <<endl;
	  MeshError(1);
       }
   else
     {
       cerr << " Error openning file " << filename << endl;
       MeshError(1);
     }
   if(verbosity>5)
   cout << "end write" << endl;
      
}
void Triangles::Write_nop5(OFortranUnFormattedFile * f,
			   Int4 &lnop5,Int4 &nef,Int4 &lgpdn,Int4 ndsr) const 
{
  ndsr =0;
  Int4 * reft = new Int4[nbt];
  //Int4 nbInT = ;
  ConsRefTriangle(reft);
  Int4 no5l[20];

  Int4 i5 = 0;
  Int4 i,j,k=0,l5;
  //  Int4 ining=0;
  Int4 imax,imin;

  lgpdn = 0;
  nef=0;
  // construction of a liste linked  of edge
  Edge ** head  = new Edge *[nbv];
  Edge  ** link = new Edge * [nbe];
  for (i=0;i<nbv;i++)
    head[i]=0; // empty liste
  
  for (i=0;i<nbe;i++)
    { 
      j = Min(Number(edges[i][0]),Number(edges[i][1]));
      link[i]=head[j];
      head[j]=edges +i;
    }
  for ( i=0;i<nbt;i++)
    {
      no5l[0]    = 0;
      Int4 kining=0;
      Int4 ining=0;
      Int4 nmae =0;
      Int4 np=0;
      l5 = 2;
      Triangle & t = triangles[i];
      Triangle * ta; // 
      Vertex *v0,*v1,*v2,*v3;
      if (reft[i]<0) continue;
      ta = t.Quadrangle(v0,v1,v2,v3);
      if (!ta)
	{ // a triangles
	  no5l[l5++] = Max(subdomains[reft[i]].ref,(Int4) 1);
	  np = 3;
	  no5l[l5++] = np;
	  no5l[0]    = np;
	  no5l[l5++] = Number(triangles[i][0]) +1;
	  no5l[l5++] = Number(triangles[i][1]) +1;
	  no5l[l5++] = Number(triangles[i][2]) +1;
	  imax = Max3(no5l[l5-1],no5l[l5-2],no5l[l5-3]);
	  imin = Min3(no5l[l5-1],no5l[l5-2],no5l[l5-3]);
	  lgpdn = Max(lgpdn,imax-imin);
          kining=l5++;
	  // ref of 3 edges 
	  for (j=0;j<3;j++)
	    {
	      no5l[l5] = 0;
	      int i0 = (int) j;
	      int i1 = (i0+1) %3;
	      Int4 j1= no5l[4+i0];
	      Int4 j2= no5l[4+i1];
	      Int4 ji = Min(j1,j2)-1;
	      Int4 ja = j1+j2-2;
	      Edge * e=head[ji];
	      while (e)
		if(Number((*e)[0])+Number((*e)[1]) == ja) 
		   {
		     no5l[l5] = e->ref;
		     break;
		   }
		else
		  e = link[Number(e)];
	      l5++;
	    }
	  if ( no5l[l5-1] || no5l[l5-2] || no5l[l5-3]  )
	    ining=2;
	  else 
	    l5 -= 3;
	  
	  no5l[l5++] = triangles[i][0].ref();
	  no5l[l5++] = triangles[i][1].ref();
	  no5l[l5++] = triangles[i][2].ref();
	  if (ining ||  no5l[l5-1] || no5l[l5-2] || no5l[l5-3]  )
            ining= ining ? ining : 3;
	  else
	    l5 -= 3,ining=0;

	  
	}
      else if ( &t<ta)
	{ 
	  k++;
	  no5l[l5++] = Max(subdomains[reft[i]].ref,(Int4) 1);
	  np =4;
	  no5l[l5++] = np;
	  no5l[0]    = np;
	  
	  no5l[l5++] = Number(v0) +1;
	  no5l[l5++] = Number(v1) +1;
	  no5l[l5++] = Number(v2) +1;
	  no5l[l5++] = Number(v3) +1;
	  
	  imax = Max(Max(no5l[l5-1],no5l[l5-2]),Max(no5l[l5-3],no5l[l5-4]));
	  imin = Min(Min(no5l[l5-1],no5l[l5-2]),Min(no5l[l5-3],no5l[l5-4]));
	  lgpdn = Max(lgpdn,imax-imin);


          kining=l5++;
	  // ref of the 4 edges 
	  // ref of 3 edges 
	  for (j=0;j<4;j++)
	    {
	      no5l[l5] = 0;
	      int i0 = (int) j;
	      int i1 = (i0+1) %4;
	      Int4 j1= no5l[4+i0];
	      Int4 j2= no5l[4+i1];
	      Int4 ji = Min(j1,j2)-1;
	      Int4 ja = j1+j2-2;
	      Edge *e=head[ji];
	      while (e)
		if(Number((*e)[0])+Number((*e)[1]) == ja) 
		   {
		     no5l[l5] = e->ref;
		     break;
		   }
		else
		  e = link[Number(e)];
	      l5++;
	    }
	  if ( no5l[l5-1] || no5l[l5-2] || no5l[l5-3] || no5l[l5-4] )
	    ining=2;
	  else 
	    l5 -= 4;

	  no5l[l5++] = v0->ref();
	  no5l[l5++] = v1->ref();
	  no5l[l5++] = v2->ref();
	  no5l[l5++] = v3->ref();
	  if (ining || no5l[l5-1] || no5l[l5-2] || no5l[l5-3] || no5l[l5-4] )
            ining= ining ? ining : 3;
          else
	    l5 -= 4;

	}
      else l5=0;

     if (l5)
       {
	 if (ining)
	   {
	     nef ++;
	     nmae = l5-kining;
	     no5l[kining] = ining;
	   }
	 else l5--;
	 no5l[1]=nmae;
	 // for all ref  
	 for (j=kining+1;j<l5;j++)
	   {
	     no5l[j] = Abs(no5l[j]);
	     ndsr = Max(ndsr,no5l[j]);
	   }
	 
	 if (f && i < 10 && verbosity > 10)
	   { 
	     cout << " e[ " << i << " @" << i5 << "]=";
	     for (j=0;j<l5;j++)
	       cout << " " << no5l[j]; 
	     cout << endl;
	   }
	 
	 if (f)
	   for (j=0;j<l5;j++)
	     *f << no5l[j]; 
	 i5 += l5;
       }
    }
  if(verbosity>10)
  cout << "   fin write nopo 5 i5=" << i5 << " " << i5*4 << endl;
  lnop5=i5; 
  lgpdn++; // add 1
  delete [] reft;
  delete [] head;
  delete [] link;
  
}

void Triangles::Write_nopo(ostream &ff) const

{
 Int4  nef=0;
 Int4 lpgdn=0;
 Int4 ndsr=0;
 Int4 i;
 Int4 ndsd=1;
 Int4 lnop5=0;
 
 OFortranUnFormattedFile f(ff);
 
 for (i=0;i<NbSubDomains ;i++)
   ndsd=Max(ndsd,subdomains[i].ref);
 
 // to compute the lnop5,nef,lpgdn,ndsr parameter 
 Write_nop5(0,lnop5,nef,lpgdn,ndsr);
 
 f.Record();
 
 f <<  Int4(13)<<Int4(6)<<Int4(32)<<Int4(0)<<Int4(27)<<Int4(0) ;
 f << Int4(nbv+nbv) ;
 f << lnop5;
 f << Int4(1 )<<Int4(1)<<Int4(1 )<<Int4(1)<<Int4(2)<<Int4(1);

 f.Record(33*sizeof(Int4)); 

 f << Int4(32) ;
 
 //char *c=identity;
 time_t timer =time(0);
 char buf[10];
 strftime(buf ,10,"%y/%m/%d",localtime(&timer));
 f.write4(identity,20);
 f.write4(buf,2);
 f.write4("created with BAMG",6);
 f.write4("NOPO",1);

 
 f << Int4(0) << Int4(1) << Int4(0) ;
 f.Record();
 Int4 nbquad= NbOfQuad;
 Int4 nbtria= nbt-NbOutT - 2*NbOfQuad;

 cout << " lnop5      = " << lnop5 << endl;
 cout << " nbquad     = " << nbquad << endl;
 cout << " nbtrai     = " << nbtria << endl;
 cout << " lpgdn      = " << lpgdn << endl;
 cout << " nef        = " << nef  << endl;
 cout << " np         = " << nbv  << endl;
 cout << " ndsr       = " << ndsr << endl;
 f << Int4(27)  
   << Int4(2)  << ndsr     << ndsd    << Int4(1) << nbtria+nbquad
   << Int4(0)  << Int4(0)  << nbtria  << nbquad  << Int4(0)
   << Int4(0)  << Int4(0)  << Int4(0) << nef     << Int4(nbv)
   << Int4(0)  << Int4(0)  << Int4(0) << Int4(0) << Int4(0)
   << Int4(0)  << nbv      << Int4(2) << lpgdn   << Int4(0)
   << lnop5    << Int4(1) ;
  f.Record();
  f << (Int4) 2*nbv;
  for (i=0;i<nbv;i++)
    f << (float) vertices[i].r.x <<  (float) vertices[i].r.y;
  f.Record();
  f << lnop5;
  Write_nop5(&f,lnop5,nef,lpgdn,ndsr);
  // cout << "fin write nopo" << endl;
}

void Triangles::Write_am_fmt(ostream &f) const 
{
  Int4 i,j;
  assert(this && nbt);
  Int4 * reft = new Int4[nbt];
  Int4 nbInT =    ConsRefTriangle(reft);
  f.precision(12);
  f << nbv << " " << nbInT << endl;
  for (i=0;i<nbt;i++)
    if(reft[i]>=0)
      {
	f << Number(triangles[i][0]) +1 << " " ;
	f << Number(triangles[i][1]) +1 << " " ;
	f << Number(triangles[i][2]) +1 << " " ;
	f << endl;
      }
  for (i=0;i<nbv;i++)
      f << vertices[i].r.x << " " << vertices[i].r.y << endl;
   for (j=i=0;i<nbt;i++) 
     if (reft[i]>=0)
       f << subdomains[reft[i]].ref  << (j++%10 == 9 ?  '\n' : ' ');
   f << endl;
   for (i=0;i<nbv;i++)
     f << vertices[i].ref()  << (i%10 == 9 ?  '\n' : ' ');
   f << endl;
   delete [] reft;


}

void Triangles::Write_am(ostream &ff) const 
{
  OFortranUnFormattedFile f(ff);  
  Int4 i,j;
  assert(this && nbt);
  Int4 * reft = new Int4[nbt];
  Int4 nbInT =    ConsRefTriangle(reft);
  f.Record();
  f << nbv << nbInT ;
  f.Record();
  for (i=0;i<nbt;i++)
    if(reft[i]>=0)
      {
	f << Number(triangles[i][0]) +1 ;
	f << Number(triangles[i][1]) +1 ;
	f << Number(triangles[i][2]) +1 ;
      }
  for (i=0;i<nbv;i++)
    {
      float x= vertices[i].r.x;
      float y= vertices[i].r.y;
      f << x << y ;
    }
  for (j=i=0;i<nbt;i++) 
    if (reft[i]>=0)
      f << subdomains[reft[i]].ref;
  for (i=0;i<nbv;i++)
    f << vertices[i].ref() ;
  delete [] reft;
}

void Triangles::Write_ftq(ostream &f) const 
{

  Int4 i;
  assert(this && nbt);
  Int4 * reft = new Int4[nbt];
  Int4 nbInT =    ConsRefTriangle(reft);
  f.precision(12);
  Int4 nele = nbInT-NbOfQuad;
  Int4 ntri =  nbInT-2*NbOfQuad;
  Int4 nqua =  NbOfQuad;

  f << nbv << " " << nele << " " << ntri <<  " " << nqua << endl;
  Int4 k=0;
  for( i=0;i<nbt;i++)
    { 
      Triangle & t = triangles[i];
      Triangle * ta; // 
      Vertex *v0,*v1,*v2,*v3;
      if (reft[i]<0) continue;
      ta = t.Quadrangle(v0,v1,v2,v3);
      if (!ta)
	{ // a triangles
	  f << "3 " 
	    << Number(triangles[i][0]) +1 << " " 
	    << Number(triangles[i][1]) +1 << " " 
	    << Number(triangles[i][2]) +1 << " " 
	    << subdomains[reft[i]].ref << endl;
	  k++;
	}
      if ( &t<ta)
	{ 
	  k++;
	  f << "4 " << Number(v0)+1 << " " << Number(v1)+1  << " "  
	    << Number(v2)+1 << " "  << Number(v3)+1 << " "  
	    << subdomains[reft[i]].ref << endl;
	}
    }
  assert(k == nele);
  
  for (i=0;i<nbv;i++)
    f << vertices[i].r.x << " " << vertices[i].r.y 
      << " " <<  vertices[i].ref() << endl;
  delete [] reft;
  
  
}
void Triangles::Write_msh(ostream &f) const 
{
  Int4 i;
  assert(this && nbt);
  Int4 * reft = new Int4[nbt];
  Int4 nbInT =    ConsRefTriangle(reft);
  f.precision(12);
  f << nbv << " " << nbInT << " " << nbe <<  endl;

  for (i=0;i<nbv;i++)
    f << vertices[i].r.x << " " << vertices[i].r.y << " " 
      << vertices[i].ref() <<   endl;

  for (i=0;i<nbt;i++)
    if(reft[i]>=0)
      f << Number(triangles[i][0]) +1 << " " 
	<< Number(triangles[i][1]) +1 << " " 
	<< Number(triangles[i][2]) +1 << " " 
	<< subdomains[reft[i]].ref << endl;
  

  for (i=0;i<nbe;i++)
    f << Number(edges[i][0]) +1 << " "  << Number(edges[i][1]) +1 
      << " " << edges[i].ref << endl;
      
   delete [] reft;

}

void Triangles::Write_amdba(ostream &f) const 
{
  assert(this && nbt);

  Int4 i,j;
  Int4 * reft = new Int4[nbt];
  Int4 nbInT =    ConsRefTriangle(reft);
  f << nbv << " " << nbInT << endl;
  cout.precision(12);
  for (i=0;i<nbv;i++)
    f << i+1 << " " 
      << vertices[i].r.x 
      << " " << vertices[i].r.y 
      << " " << vertices[i].ref() << endl;
  j=1;
  for (i=0;i<nbt;i++)
    if(reft[i]>=0)
	f << j++ << " " 
	  << Number(triangles[i][0]) +1 << " " 
	  << Number(triangles[i][1]) +1 << " " 
	  << Number(triangles[i][2]) +1 << " " 
	  << subdomains[reft[i]].ref  << endl ;
	f << endl;
   delete [] reft;


}

void Triangles::Write(const char * filename)
{
  ofstream f(filename);
  if (f)
    {
       if (name) delete name;
       name = new char[strlen(filename)+1];
       strcpy(name,filename);
       OnDisk =1;
       f << *this;
    }
}
void Triangles::WriteElements(ostream& f,Int4 * reft ,Int4 nbInT) const
   { 
     const Triangles & Th= *this;
     // do triangle and quad 
     if(verbosity>9) 
       cout  << " In Triangles::WriteElements " << endl
	     << "   Nb of In triangles " << nbInT-Th.NbOfQuad*2 << endl
	     << "   Nb of Quadrilaterals " <<  Th.NbOfQuad << endl
	     << "   Nb of in+out+quad  triangles " << Th.nbt << " " << nbInT << endl;
	 
     Int4 k=nbInT-Th.NbOfQuad*2;
     Int4 num =0;
     if (k>0) {
       f << "\nTriangles\n"<< k << endl;
       for(Int4 i=0;i<Th.nbt;i++)
	 { 
	   Triangle & t = Th.triangles[i];
	   if (reft[i]>=0 && !( t.Hidden(0) || t.Hidden(1) || t.Hidden(2) ))
	    { k--;
	      f << Th.Number(t[0])+1 << " " << Th.Number(t[1])+1 
		   << " "  << Th.Number(t[2])+1  << " " << Th.subdomains[reft[i]].ref << endl;
		   reft[i] = ++num;
		 }
	 }
     } 
     if (Th.NbOfQuad>0) {
       f << "\nQuadrilaterals\n"<<Th.NbOfQuad << endl;
       k = Th.NbOfQuad;
       for(Int4 i=0;i<Th.nbt;i++)
	 { 
	   Triangle & t = Th.triangles[i];
	   Triangle * ta; // 
	   Vertex *v0,*v1,*v2,*v3;
	   if (reft[i]<0) continue;
	   if ((ta=t.Quadrangle(v0,v1,v2,v3)) !=0 && &t<ta)
	      { 
		k--;
		f << Th.Number(v0)+1 << " " << Th.Number(v1)+1  << " "  
		  << Th.Number(v2)+1 << " "  << Th.Number(v3)+1 << " "  
		  << Th.subdomains[reft[i]].ref << endl;
		  reft[i] = ++num;
		  reft[Number(ta)] = num;
	      }
	 }
       assert(k==0);
     }
     // warning reft is now the element number 
   }

ostream& operator <<(ostream& f, const   Triangles & Th) 
 {
  //  Th.FindSubDomain();
   // warning just on say the class is on the disk
  //  ((Triangles *) &Th)->OnDisk = 1;

   Int4 * reft = new Int4[Th.nbt];
   Int4 nbInT =    Th.ConsRefTriangle(reft);
   {
     f << "MeshVersionFormatted 0" <<endl;
     f << "\nDimension\n"  << 2 << endl;
     f << "\nIdentifier\n" ;
     WriteStr(f,Th.identity);
     f << "\n\nGeometry\n" ;
     if( Th.Gh.OnDisk)
       WriteStr(f,Th.Gh.name),     f <<endl;
     else
       { // empty file name -> geom in same file
	 f << "\"\"" << endl << endl;
	 f << "# BEGIN of the include geometry file because geometry is not on the disk"
	   << Th.Gh << endl;
	 f << "End" << endl 
	   << "# END of the include geometrie file because geometry is not on the disk"
	   << endl ;
       }
   }
   { 
     f.precision(12);
     f << "\nVertices\n" << Th.nbv <<endl;
     for (Int4 i=0;i<Th.nbv;i++)
       {
	 Vertex & v =  Th.vertices[i];
	 f << v.r.x << " " << v.r.y << " " << v.ref() << endl;
       }
   }
  Int4 ie; 
   {
     f << "\nEdges\n"<< Th.nbe << endl;
     for(ie=0;ie<Th.nbe;ie++)
       { 
	 Edge & e = Th.edges[ie];
	 f << Th.Number(e[0])+1 << " " << Th.Number(e[1])+1;
	f << " " << e.ref <<endl;
       }
     if(Th.NbCrackedEdges)
       {
	 f << "\nCrackedEdges\n"<< Th.NbCrackedEdges << endl;
	 for( ie=0;ie<Th.NbCrackedEdges;ie++)
	   { 
	     Edge & e1 = *Th.CrackedEdges[ie].a.edge;
	     Edge & e2 = *Th.CrackedEdges[ie].b.edge;
	     f << Th.Number(e1)+1 << " " << Th.Number(e2)+1 <<endl;;
	   }
       }
   }

   Th.WriteElements(f,reft,nbInT);
   {
     f << "\nSubDomainFromMesh\n" << Th.NbSubDomains<< endl ;
     for (Int4 i=0;i<Th.NbSubDomains;i++)
       f << 3 << " " << reft[Th.Number(Th.subdomains[i].head)] << " " << 1 << " " 
	 <<  Th.subdomains[i].ref << endl;
     
   }
   if (Th.Gh.NbSubDomains)
     {
        f << "\nSubDomainFromGeom\n" << Th.Gh.NbSubDomains << endl ;
       for (Int4 i=0;i<Th.NbSubDomains;i++)
	 {  
	 f << 2 << " " << Th.Number(Th.subdomains[i].edge)+1 << " " 
	   <<  Th.subdomains[i].sens  << " " <<  Th.Gh.subdomains[i].ref << endl;
	 } 
   }
   {
     f << "\nVertexOnGeometricVertex\n"<<  Th.NbVerticesOnGeomVertex << endl;
     for (Int4 i0=0;i0<Th.NbVerticesOnGeomVertex;i0++)
       {
	 VertexOnGeom & v =Th.VerticesOnGeomVertex[i0];
	 assert(v.OnGeomVertex()) ;
	 f << " " << Th.Number(( Vertex *)v)+1  
	   << " " << Th.Gh.Number(( GeometricalVertex * )v)+1 
	   << endl;
       }
   }
   { 
     f << "\nVertexOnGeometricEdge\n"<<  Th.NbVerticesOnGeomEdge << endl;
     for (Int4 i0=0;i0<Th.NbVerticesOnGeomEdge;i0++)
       {
	 const VertexOnGeom & v =Th.VerticesOnGeomEdge[i0];
	 assert(v.OnGeomEdge()) ;   
	 f << " " << Th.Number((Vertex * )v)+1  ;
	 f << " " << Th.Gh.Number((const  GeometricalEdge * )v)+1  ;
	 f << " " << (Real8 ) v << endl;
       }
   }
   {
     Int4 i0,k=0;

     for (i0=0;i0<Th.nbe;i0++)
       if ( Th.edges[i0].on ) k++;
     
     f << "\nEdgeOnGeometricEdge\n"<< k << endl;
      for (i0=0;i0<Th.nbe;i0++)
	if ( Th.edges[i0].on ) 
	  f << (i0+1) << " "  << (1+Th.Gh.Number(Th.edges[i0].on)) <<  endl;
      if (Th.NbCrackedEdges)
	{
	  f << "\nCrackedEdges\n"<< Th.NbCrackedEdges << endl;	  
	  for(i0=0;i0< Th.NbCrackedEdges; i0++) 
	    {
	      f << Th.Number(Th.CrackedEdges[i0].a.edge) << " " ;
	      f  << Th.Number(Th.CrackedEdges[i0].b.edge) << endl;
	    }
	}
   }  
   if (&Th.BTh != &Th && Th.BTh.OnDisk && Th.BTh.name) 
     {
       int *mark=new int[Th.nbv];
       Int4 i;
       for (i=0;i<Th.nbv;i++)
	 mark[i]=-1;
       f << "\nMeshSupportOfVertices\n" <<endl;
       WriteStr(f,Th.BTh.name);
       f <<endl;
       f << "\nIdentityOfMeshSupport" << endl;
       WriteStr(f,Th.BTh.identity);
       f<<endl;

       f << "\nVertexOnSupportVertex" << endl;
       f<< Th.NbVertexOnBThVertex << endl;
       for(i=0;i<Th.NbVertexOnBThVertex;i++) {
	 const VertexOnVertex & vov = Th.VertexOnBThVertex[i];
	 Int4 iv = Th.Number(vov.v);
	 mark[iv] =0;
	 f << iv+1<< " " << Th.BTh.Number(vov.bv)+1 << endl;}

       f << "\nVertexOnSupportEdge" << endl;
       f << Th.NbVertexOnBThEdge << endl;
       for(i=0;i<Th.NbVertexOnBThEdge;i++) {
	 const VertexOnEdge & voe = Th.VertexOnBThEdge[i];
	 Int4 iv = Th.Number(voe.v);
	 //	 assert(mark[iv] == -1]);
	 mark[iv] = 1;
	 f << iv+1 << " " << Th.BTh.Number(voe.be)+1 << " " << voe.abcisse <<  endl;}
       
       f << "\nVertexOnSupportTriangle" << endl;   
       Int4 k = Th.nbv -  Th.NbVertexOnBThEdge - Th.NbVertexOnBThVertex;
       f << k << endl;
       //       Int4 kkk=0;
       CurrentTh=&Th.BTh;
       for (i=0;i<Th.nbv;i++) 
	 if (mark[i] == -1) {
	   k--;
	   Icoor2 dete[3];
	   I2 I = Th.BTh.toI2(Th.vertices[i].r);
	   Triangle * tb = Th.BTh.FindTriangleContening(I,dete);
	   if (tb->link) // a true triangle
	     {
	       Real8 aa= (Real8) dete[1]/ tb->det, bb= (Real8) dete[2] / tb->det;
	       f << i+1 << " " << Th.BTh.Number(tb)+1 << " " << aa << " " << bb << endl ;
	     }
	   else 
	     {
	       double aa,bb,det[3];
	       TriangleAdjacent ta=CloseBoundaryEdgeV2(I,tb,aa,bb);
	       int k = ta;
	       det[VerticesOfTriangularEdge[k][1]] =aa;
	       det[VerticesOfTriangularEdge[k][0]] = bb;
	       det[OppositeVertex[k]] = 1- aa -bb;
	       Triangle * tb = ta;
	       f << i+1 << Th.BTh.Number(tb)+1 << " " << det[1] << " " << det[2] <<endl;
	     }
	 }
       assert(!k);
       delete [] mark;
	 

     }
   f << "\nEnd" << endl;
   //  Th.ConsLinkTriangle();
   delete [] reft;
   return f;
   
}



void Geometry::Write(const char * filename)
{
  ofstream f(filename);
  if (f)
    {
      if(verbosity>1)
	cout << "  -- write geometry in file " << filename << endl;
       if (name) delete name;
       name = new char[strlen(filename)+1];
       strcpy(name,filename);
       OnDisk =1;
       f << *this;
    }
}

ostream& operator <<(ostream& f, const   Geometry & Gh) 
{
   Int4  NbCorner=0;
   {
     f << "MeshVersionFormatted 0" <<endl;
     f << "\nDimension\n"  << 2 << endl;
//     f << "\nIdentifier\n" ;
//     WriteStr(f,Gh.identity);
//     f <<endl;
   }
   int nbreqv=0;
   { 
     
     f.precision(12);
     f << "\nVertices\n" << Gh.nbv <<endl;
     for (Int4 i=0;i<Gh.nbv;i++)
       {
	 GeometricalVertex & v =  Gh.vertices[i];
	 if (v.Required()) nbreqv++;
	 f << v.r.x << " " << v.r.y << " " << v.ref() << endl;
	 if (v.Corner()) NbCorner++;
       }
   }
   
   int nbcracked=0;

   {
     int nbreq=0;
     f << "\nEdges\n"<< Gh.nbe << endl;
     for(Int4 ie=0;ie<Gh.nbe;ie++)
       { 
	 
	 GeometricalEdge & e = Gh.edges[ie];
	 if (e.Required()) nbreq++;
	 if (e.Cracked()) { 
	   Int4 ie1 = Gh.Number(e.link);
	   if (ie <= ie1)  ++nbcracked;}
	 f << Gh.Number(e[0])+1 << " " << Gh.Number(e[1])+1;
	 f << " " << e.ref <<endl;
       }
     
     if (nbcracked)
       {
	 f << "\nCrackedEdges\n"<< nbcracked<< endl;
	 for(Int4 ie=0;ie<Gh.nbe;ie++)
	   {
	     GeometricalEdge & e = Gh.edges[ie];
	     if (e.Cracked()) { 
	       Int4  ie1 = Gh.Number(e.link);
	       if (ie <= ie1)  f << ie+1 << " " << ie1+1<< endl;
	     }
	   }
       }
     if(nbreq)
       {
	 f << "\nRequiredEdges\n"<< nbreq<< endl;
         for(Int4 ie=0;ie<Gh.nbe;ie++)
           {
             GeometricalEdge & e = Gh.edges[ie];
             if (e.Required()) 
	       f << ie+1 << endl;
	   }
       }
     
     
     
   }

    f << "\nAngleOfCornerBound\n" 
     << Gh.MaximalAngleOfCorner*180/Pi << endl;
    if (NbCorner) 
      {
	f << "\nCorners\n" << NbCorner << endl;
	for (Int4 i=0,j=0;i<Gh.nbv;i++)
	  {
	    GeometricalVertex & v =  Gh.vertices[i];
	    if (v.Corner()) 
	      j++,f << Gh.Number(v)+1 << (j % 5 ? ' ' : '\n');
	  }
        
      
      }

    if(nbreqv)
      {
	f << "\nRequiredVertices\n"<< nbreqv<< endl;
	for (Int4 j=0,i=0;i<Gh.nbv;i++)
	  {
	    GeometricalVertex & v =  Gh.vertices[i];
	    
	    if (v.Required()) 
	      j++,f << i+1 << (j % 5 ? ' ' : '\n');
	  }
	f << endl;
      }
    
    { 
       Int4 i;
       f << "\nSubDomainFromGeom\n" ;
       f << Gh.NbSubDomains<< endl;
       for (i=0;i<Gh.NbSubDomains;i++) 
         f << "2 " << Gh.Number(Gh.subdomains[i].edge)+1 << " " << Gh.subdomains[i].sens 
           << " " << Gh.subdomains[i].ref << endl;        
     }
     {
       Int4 n=0,i;

       for(i=0;i< Gh.nbe;i++)
	 {
	   if(Gh.edges[i].TgA() && Gh.edges[i][0].Corner() ) 
	     n++;
	   if(Gh.edges[i].TgB() && Gh.edges[i][1].Corner() ) 
	     n++;
	 }
       if (n) {
	 f << "TangentAtEdges " << n << endl;
	 for(i=0;i< Gh.nbe;i++)
	   {
	     if (Gh.edges[i].TgA() && Gh.edges[i][0].Corner() ) 
	       f << i+1 << " 1 " << Gh.edges[i].tg[0].x 
		 << " " << Gh.edges[i].tg[0].y << endl;
	     if (Gh.edges[i].TgB() && Gh.edges[i][1].Corner() ) 
	       f << i+1 << " 2 " << Gh.edges[i].tg[1].x 
		 << " " << Gh.edges[i].tg[1].y << endl;
	   }
	 
       }}
     //  f << " Not Yet Implemented" << endl;
     
     return f;
}


} // end of namespace bamg 
