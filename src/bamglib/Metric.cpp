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
#include "Meshio.h"
#include "Mesh2.h"


namespace bamg {


inline Real8 det3x3(Real8 A[3] ,Real8 B[3],Real8 C[3])
{ return    A[0] * ( B[1]*C[2]-B[2]*C[1])
          - A[1] * ( B[0]*C[2]-B[2]*C[0])
          + A[2] * ( B[0]*C[1]-B[1]*C[0]);
}

SaveMetricInterpole  LastMetricInterpole;

void ReductionSimultanee( MetricAnIso M1,  MetricAnIso M2,double & l1,double & l2, D2xD2 & V) 
{
  double a11=M1.a11,a21=M1.a21,a22=M1.a22;
  double b11=M2.a11,b21=M2.a21,b22=M2.a22;
  //  M1 v = l M2 v
  // (M1 - l M2) v =0
  // det (M1 - l M2) =0
  // det (M1 - l M2) = a l^2 + b l + c;
  // = (a11 - l * b11) * (a22 - l * b22) - (a21 - l * b21 ) ^2
  //  const double eps = 1.e-5;
  const double /*c11 = a11*a11,*/ c21= a21*a21;
  const double /*d11 = b11*b11,*/ d21= b21*b21;
  const double a=b11*b22 - d21;
  const double b=-a11*b22-a22*b11+2*a21*b21;
  const double c=-c21+a11*a22;
  const double bb = b*b,ac= a*c;
  const double delta = bb - 4 * ac;
  //  const double kk=c11+c22+c21+d11+d21+d22;
  // modif F Hecht feb 1998 
  // cerr.precision(14);
  //cerr  <<  bb << " " << ac << " " <<  bb <<  " " <<a << endl;
  // cerr << a11 << " " << a21 << " " << a22 << endl;
  //cerr << b11 << " " << b21 << " " << b22 << endl;
  if (bb + Abs(ac) < 1.0e-20 || (delta< 1.0E-4 * bb ) )
   {
   // racine double;
     // cerr << "double " << endl ;
    if (Abs(a) < 1.e-30 )
     l1 = l2 = 0;
    else 
     l1=l2=-b/(2*a); 
    V= D2xD2(1,0,0,1);
   }
  else {
    // cerr << "  -- " << a << endl ;
    const double delta2 = sqrt(delta);
    l1= (-b - delta2)/(2*a);
    l2= (-b + delta2)/(2*a);
    // M1 v = l M2 v
    //  ( (M1 - I M2) x,y)  = (x,(M1 - I M2) y) \forall y
    // so Ker((M1 - I M2)) = Im((M1 - I M2))^\perp
      double v0 = a11-l1*b11, v1 = a21-l1*b21,v2 = a22 - l1*b22;
      double s0 = v0*v0 + v1*v1, s1 = v1*v1 +v2*v2;
      double vp1x,vp1y,vp2x,vp2y;

      if(s1 < s0)
	s0=sqrt(s0),vp1x=v1/s0,vp1y=-v0/s0;
      else
	s1=sqrt(s1),vp1x=v2/s1,vp1y=-v1/s1;

      v0 = a11-l2*b11, v1 = a21-l2*b21,v2 = a22 - l2*b22;
      s0 = v0*v0 + v1*v1, s1 = v1*v1 +v2*v2;
      if(s1 < s0)
	s0=sqrt(s0),vp2x=v1/s0,vp2y=-v0/s0;
      else
	s1=sqrt(s1),vp2x=v2/s1,vp2y=-v1/s1;
#ifdef DEBUG
      assert(Abs(vp1y)+Abs(vp2y)>0);
#endif
      V=D2xD2(vp1x,vp2x,vp1y,vp2y);
  }
  return;

}

MetricAnIso Intersection(const MetricAnIso M1,const MetricAnIso M2) ;
MetricAnIso Intersection(const MetricAnIso M1,const MetricAnIso M2) 
{
      D2xD2 M;
      double l1,l2;
      ReductionSimultanee(M1,M2,l1,l2,M);
      R2 v0(M.x.x,M.y.x);
      R2 v1(M.x.y,M.y.y);
      D2xD2 M_1(M.inv());
      D2xD2 D(Max(M1(v0,v0),M2(v0,v0)),0,0,Max(M1(v1,v1),M2(v1,v1)));
      D2xD2 Mi(M_1.t()*D*M_1);
      return MetricAnIso(Mi.x.x,0.5*(Mi.x.y+Mi.y.x),Mi.y.y);
}

MetricAnIso::MetricAnIso(const Real8  a[3],const  MetricAnIso m0,
	   const  MetricAnIso m1,const  MetricAnIso m2 )
{
  MetricAnIso mab(a[0]*m0.a11 + a[1]*m1.a11 + a[2]*m2.a11,
		  a[0]*m0.a21 + a[1]*m1.a21 + a[2]*m2.a21,
		  a[0]*m0.a22 + a[1]*m1.a22 + a[2]*m2.a22);
  
  MatVVP2x2 vab(mab);
 
  R2 v1(vab.v.x,vab.v.y);
  R2 v2(-v1.y,v1.x);
  
  Real8 h1 = a[0] / m0(v1) + a[1] / m1(v1) + a[2] / m2(v1);
  Real8 h2 = a[0] / m0(v2) + a[1] / m1(v2) + a[2] / m2(v2);

  vab.lambda1 =  1 / (h1*h1);
  vab.lambda2 =  1 / (h2*h2);
  *this = vab;
}

 MetricAnIso::MetricAnIso( Real8  a,const  MetricAnIso ma,
	                   Real8  b,const  MetricAnIso mb)
{ 
  MetricAnIso mab(a*ma.a11+b*mb.a11,a*ma.a21+b*mb.a21,a*ma.a22+b*mb.a22);
  MatVVP2x2 vab(mab);
  
  R2 v1(vab.v.x,vab.v.y);
  R2 v2(-v1.y,v1.x);
  

  Real8 h1 = a / ma(v1) + b / mb(v1);
  Real8 h2 = a / ma(v2) + b / mb(v2);
  vab.lambda1 =  1 / (h1*h1);
  vab.lambda2 =  1 / (h2*h2);
  *this = vab;
}



 MatVVP2x2::MatVVP2x2(const MetricAnIso M) 
{
  double a11=M.a11,a21=M.a21,a22=M.a22;
  const double eps = 1.e-5;
  double c11 = a11*a11, c22 = a22*a22, c21= a21*a21;
  double b=-a11-a22,c=-c21+a11*a22;
  double   delta = b*b - 4 * c ;
  double n2=(c11+c22+c21);
  if ( n2 < 1e-30) 
     lambda1=lambda2=0,v.x=1,v.y=0;
  else if (delta < eps*n2)
    { 
      lambda1=lambda2=-b/2, v.x=1,v.y=0;
    }
  else 
    {  //    ---  construction  de 2 vecteur dans (Im ( A - D(i) Id) ortogonal 
      delta = sqrt(delta);
      lambda1 = (-b-delta)/2.0,lambda2 = (-b+delta)/2.0;
      double v0 = a11-lambda1, v1 = a21,v2 = a22 - lambda1;
      double s0 = v0*v0 + v1*v1, s1 = v1*v1 +v2*v2;
 
      if(s1 < s0)
	s0=sqrt(s0),v.x=v1/s0,v.y=-v0/s0;
      else
	s1=sqrt(s1),v.x=v2/s1,v.y=-v1/s1;
    };
}


 int MetricAnIso::IntersectWith(const MetricAnIso M2) 
{
  //cerr << " - " << *this << M2 <<  endl;
      int r=0;
      MetricAnIso & M1 = *this;
      D2xD2 M;
      double l1,l2;
      
      ReductionSimultanee(*this,M2,l1,l2,M);
      // cerr << M << endl;
      R2 v1(M.x.x,M.y.x);
      R2 v2(M.x.y,M.y.y);
      double l11=M1(v1,v1);
      double l12=M1(v2,v2);
      double l21=M2(v1,v1);
      double l22=M2(v2,v2);
      if ( l11 < l21 )  r=1,l11=l21;
      if ( l12 < l22 )  r=1,l12=l22; 
      // cerr << r << endl;
      if (r) { // change
        D2xD2 M_1(M.inv());
        D2xD2 D(l11,0,0,l12); 
        D2xD2 Mi(M_1.t()*D*M_1);
        a11=Mi.x.x;
        a21=0.5*(Mi.x.y+Mi.y.x);
        a22=Mi.y.y; }
      return r;
}
void Triangles::IntersectGeomMetric(const Real8 err=1,const int iso=0)

{
  if(verbosity>1)
    cout << "  -- IntersectGeomMetric geometric err=" << err << (iso ? " iso " : " aniso "  ) << endl;
  Real8 ss[2]={0.00001,0.99999};
  Real8 errC = 2*sqrt(2*err);
  Real8 hmax = Gh.MaximalHmax();
  Real8 hmin = Gh.MinimalHmin();
  Real8 maxaniso = 1e6;
  assert(hmax>0);
  SetVertexFieldOn();
  if (errC > 1) errC = 1;
  for (Int4  i=0;i<nbe;i++)
   for (int j=0;j<2;j++)
    {
      
      Vertex V;
      VertexOnGeom GV;
      // cerr << Number(edges[i]) << " " << ss[j] << endl;
      Gh.ProjectOnCurve(edges[i],ss[j],V,GV);
	{
	  GeometricalEdge * eg = GV;
	  Real8 s = GV;
	  R2 tg;
	  //	   cerr << i << " " << j << " " << Number(V) << " on = " 
	  //	<< Gh.Number(eg) << " at s = " << s << " " << endl;
	  Real8  R1= eg->R1tg(s,tg);
	  // cerr << " R = " << 1/Max(R1,1e-20) << tg << " on x " 
	  //    << V.r << errC/ Max(R1,1e-20) <<  " hold=" <<V.m(tg) << " "  << endl;
	  Real8 ht = hmax;
          if (R1>1.0e-20) 
	    {  // err relative to the length of the edge
	      ht = Min(Max(errC/R1,hmin),hmax);
	    }
	  Real8 hn = iso? ht : Min(hmax,ht*maxaniso);
	  //cerr << ht << " " << hn << "m=" << edges[i][j].m <<  endl;
	  assert(ht>0 && hn>0);
	  MatVVP2x2 Vp(1/(ht*ht),1/(hn*hn),tg);
	  //cerr << " : " ;
	  Metric MVp(Vp);
	  // cerr << " : "  << MVp  << endl;
	  edges[i][j].m.IntersectWith(MVp);
	  //cerr << " . " << endl;
	}

    }
  // the problem is for the vertex on vertex 
     
}
/*
void  Triangles::BoundAnisotropy(Real8 anisomax)
{
  if (verbosity > 1) 
    cout << "  -- BoundAnisotropy by  " << anisomax << endl; 
  Real8 h1=1.e30,h2=1e-30,rx=0;
  Real8 coef = 1./(anisomax*anisomax);
  Real8 hn1=1.e30,hn2=1e-30,rnx =1.e-30;  
  for (Int4 i=0;i<nbv;i++)
    {

      MatVVP2x2 Vp(vertices[i]);
      
      h1=Min(h1,Vp.lmin());
      h2=Max(h2,Vp.lmax());
      rx = Max(rx,Vp.Aniso2());
      
      Vp.BoundAniso2(coef);
      
      hn1=Min(hn1,Vp.lmin());
      hn2=Max(hn2,Vp.lmax());
      rnx = Max(rnx,Vp.Aniso2());

      
      vertices[i].m = Vp;

    }

  if (verbosity>2)
    {
      cout << "     input :  Hmin = " << sqrt(1/h2)  << " Hmax = " << sqrt(1/h1) 
	   << " factor of anisotropy max  = " << sqrt(rx) << endl;
      cout << "     output:  Hmin = " << sqrt(1/hn2) << " Hmax = " << sqrt(1/hn1) 
	   << " factor of anisotropy max  = " << sqrt(rnx) << endl;
    }
}
*/
void  Triangles::BoundAnisotropy(Real8 anisomax,Real8 hminaniso)
{
  double lminaniso = 1/ (Max(hminaniso*hminaniso,1e-100));
  if (verbosity > 1) 
    cout << "  -- BoundAnisotropy by  " << anisomax << endl; 
  Real8 h1=1.e30,h2=1e-30,rx=0;
  Real8 coef = 1./(anisomax*anisomax);
  Real8 hn1=1.e30,hn2=1e-30,rnx =1.e-30;  
  for (Int4 i=0;i<nbv;i++)
    {

      MatVVP2x2 Vp(vertices[i]);
      double lmax=Vp.lmax();
      h1=Min(h1,Vp.lmin());
      h2=Max(h2,Vp.lmax());
      rx = Max(rx,Vp.Aniso2());

      Vp *= Min(lminaniso,lmax)/lmax;
      
      Vp.BoundAniso2(coef);
      
      hn1=Min(hn1,Vp.lmin());
      hn2=Max(hn2,Vp.lmax());
      rnx = Max(rnx,Vp.Aniso2());

      
      vertices[i].m = Vp;

    }

  if (verbosity>2)
    {
      cout << "     input :  Hmin = " << sqrt(1/h2)  << " Hmax = " << sqrt(1/h1) 
	   << " factor of anisotropy max  = " << sqrt(rx) << endl;
      cout << "     output:  Hmin = " << sqrt(1/hn2) << " Hmax = " << sqrt(1/hn1) 
	   << " factor of anisotropy max  = " << sqrt(rnx) << endl;
    }
}
void Triangles::IntersectConsMetric(const double * s,const Int4 nbsol,const int * typsols,
				    const  Real8 hmin1,const Real8 hmax1,const Real8 coef,
				    const Real8 anisomax ,const Real8 CutOff,const int NbJacobi,
				    const int DoNormalisation,const double power,const int choice)
{ //  the array of solution s is store    
  // sol0,sol1,...,soln    on vertex 0
  //  sol0,sol1,...,soln   on vertex 1
  //  etc.
  //  choise = 0 =>  H is computed with green formule
  //   otherwise  => H is computed from P2 on 4T 
  const int dim = 2;
  
  int sizeoftype[] = { 1, dim ,dim * (dim+1) / 2, dim * dim } ; 

  // computation of the nb of field 
  Int4 ntmp = 0;
  if (typsols)
    {
      for (Int4 i=0;i<nbsol;i++)
	     ntmp += sizeoftype[typsols[i]];
    }
  else
    ntmp = nbsol;

  // n is the total number of fields

  const Int4 n = ntmp;

  Int4 i,k,iA,iB,iC,iv;
  R2 O(0,0);
  int RelativeMetric = CutOff>1e-30;
  Real8 hmin = Max(hmin1,MinimalHmin());
  Real8 hmax = Min(hmax1,MaximalHmax());
  Real8 coef2 = 1/(coef*coef);

  if(verbosity>1) 
    {
      cout << "  -- Construction of Metric: Nb of field. " << n << " nbt = " 
	   << nbt << " nbv= " << nbv 
	   << " coef = " << coef << endl
	   << "     hmin = " << hmin << " hmax=" << hmax 
	   << " anisomax = " << anisomax <<  " Nb Jacobi " << NbJacobi << " Power = " << power ;
      if (RelativeMetric)
	cout << " RelativeErr with CutOff= "  <<  CutOff << endl;
      else
	cout << " Absolute Err" <<endl;
    }
  double *ss=(double*)s;//, *ssiii = ss;

  double sA,sB,sC;

  Real8 *detT = new Real8[nbt];
  Real8 *Mmass= new Real8[nbv];
  Real8 *Mmassxx= new Real8[nbv];
  Real8 *dxdx= new Real8[nbv];
  Real8 *dxdy= new Real8[nbv];
  Real8 *dydy= new Real8[nbv];
  Real8 *workT= new Real8[nbt];
  Real8 *workV= new Real8[nbv];
  int *OnBoundary = new int[nbv];
  for (iv=0;iv<nbv;iv++)
    {
      Mmass[iv]=0;
      OnBoundary[iv]=0;
      Mmassxx[iv]=0;
    }

  for (i=0;i<nbt;i++) 
    if(triangles[i].link) // the real triangles 
      {
	const Triangle &t=triangles[i];
	// coor of 3 vertices 
	R2 A=t[0];
	R2 B=t[1];
	R2 C=t[2];


	// number of the 3 vertices
	iA = Number(t[0]);
	iB = Number(t[1]);
	iC = Number(t[2]);
	
	Real8 dett = bamg::Area2(A,B,C);
	detT[i]=dett;
	dett /= 6;

	// construction of on boundary 
	int nbb =0;
	for(int j=0;j<3;j++)
          {
	    Triangle *ta=t.Adj(j);
	    if ( ! ta || !ta->link) // no adj triangle => edge on boundary
	      OnBoundary[Number(t[VerticesOfTriangularEdge[j][0]])]=1,
		OnBoundary[Number(t[VerticesOfTriangularEdge[j][1]])]=1,
		nbb++;
	  }
	
	workT[i] = nbb;
	Mmass[iA] += dett;
	Mmass[iB] += dett;
	Mmass[iC] += dett;
	
	if((nbb==0)|| !choice)
	  {
	    Mmassxx[iA] += dett;
	    Mmassxx[iB] += dett;
	    Mmassxx[iC] += dett;
	  }
      }
  else
    workT[i]=-1;

//  for (Int4 kcount=0;kcount<n;kcount++,ss++)
    for (Int4 nusol=0;nusol<nbsol;nusol++)
    { //for all Solution  

      Real8 smin=ss[0],smax=ss[0];
      
      Real8 h1=1.e30,h2=1e-30,rx=0;
      Real8 coef = 1./(anisomax*anisomax);
      Real8 hn1=1.e30,hn2=1e-30,rnx =1.e-30;  
      int nbfield = typsols? sizeoftype[typsols[nusol]] : 1; 
      if (nbfield == 1) 
       for ( iv=0,k=0; iv<nbv; iv++,k+=n )
				{
				  dxdx[iv]=dxdy[iv]=dydy[iv]=0;
				  smin=Min(smin,ss[k]);
				  smax=Max(smax,ss[k]);
				 }
			  else
			   {
         //  cas vectoriel 
          for ( iv=0,k=0; iv<nbv; iv++,k+=n )
          {	
           double v=0;		     
			     for (int i=0;i<nbfield;i++) 
			         v += ss[k+i]*ss[k+i];
			     v = sqrt(v);
				   smin=Min(smin,v);
				   smax=Max(smax,v);
			    }
			   }
      Real8 sdelta = smax-smin;
      Real8 absmax=Max(Abs(smin),Abs(smax));
      Real8 cnorm = DoNormalisation ? coef2/sdelta : coef2;
      
      if(verbosity>2) 
	     cout << "    Solution " << nusol <<  " Min = " << smin << " Max = " 
	       << smax << " Delta =" << sdelta << " cnorm = " << cnorm <<  " Nb of fields =" << nbfield << endl;

      
      if ( sdelta < 1.0e-10*Max(absmax,1e-20) && (nbfield ==1)) 
				{
				  if (verbosity>2)
				    cout << "      Solution " << nusol << " is constant. We skip. " 
					 << " Min = " << smin << " Max = " << smax << endl;
				continue;
				}
				
	 double *sf  = ss; 
	 for (Int4 nufield=0;nufield<nbfield;nufield++,ss++) 
	   {
	     for ( iv=0,k=0; iv<nbv; iv++,k+=n )
		       dxdx[iv]=dxdy[iv]=dydy[iv]=0;
       for (i=0;i<nbt;i++) 
	      if(triangles[i].link)
	  {// for real all triangles 
	    // coor of 3 vertices 
	    R2 A=triangles[i][0];
	    R2 B=triangles[i][1];
	    R2 C=triangles[i][2];
	    
	    
	    // warning the normal is internal and the 
	    //   size is the length of the edge
	    R2 nAB = Orthogonal(B-A);
	    R2 nBC = Orthogonal(C-B);
	    R2 nCA = Orthogonal(A-C);
	    // remark :  nAB + nBC + nCA == 0 

	    // number of the 3 vertices
	    iA = Number(triangles[i][0]);
	    iB = Number(triangles[i][1]);
	    iC = Number(triangles[i][2]);
	    
	    // for the test of  boundary edge
	    // the 3 adj triangles 
	    Triangle *tBC = triangles[i].TriangleAdj(OppositeEdge[0]);
	    Triangle *tCA = triangles[i].TriangleAdj(OppositeEdge[1]);
	    Triangle *tAB = triangles[i].TriangleAdj(OppositeEdge[2]);

	    // value of the P1 fonction on 3 vertices 
	    sA = ss[iA*n];
	    sB = ss[iB*n];
	    sC = ss[iC*n];

	    R2 Grads = (nAB * sC + nBC * sA + nCA * sB ) /detT[i] ;
	    if(choice) 
	      {
		int nbb = 0;
		Real8 dd = detT[i];
		Real8 lla,llb,llc,llf;
		Real8  taa[3][3],bb[3];
		// construction of the trans of lin system
		for (int j=0;j<3;j++)
		  {
		    int ie = OppositeEdge[j];
		    TriangleAdjacent ta = triangles[i].Adj(ie);
		    Triangle *tt = ta;
		    if (tt && tt->link)
		      {
			Vertex &v = *ta.OppositeVertex();
			R2 V = v;
			Int4 iV = Number(v);
			Real8 lA  = bamg::Area2(V,B,C)/dd;
			Real8 lB  = bamg::Area2(A,V,C)/dd;
			Real8 lC  = bamg::Area2(A,B,V)/dd;
			taa[0][j] =  lB*lC;
			taa[1][j] =  lC*lA;
			taa[2][j] =  lA*lB;
			//Real8 xx = V.x-V.y;
			//Real8 yy = V.x + V.y;
			//cout << " iv " << ss[iV*n] << " == " << (8*xx*xx+yy*yy)
			//     << " l = " << lA << " " << lB << " " << lC 
			//     << " = " << lA+lB+lC << " " <<  V << " == " << A*lA+B*lB+C*lC << endl;
			
			lla = lA,llb=lB,llc=lC,llf=ss[iV*n] ;

			bb[j]     =  ss[iV*n] - ( sA*lA + sB*lB + sC*lC ) ;
		      }
		    else
		      {
			nbb++;
			taa[0][j]=0;
			taa[1][j]=0;
			taa[2][j]=0;
			taa[j][j]=1;
			bb[j]=0;
		      }
		  }

		// resolution of 3x3 lineaire system transpose
		Real8 det33 =  det3x3(taa[0],taa[1],taa[2]);		
		Real8 cBC   =  det3x3(bb,taa[1],taa[2]);
		Real8 cCA   =  det3x3(taa[0],bb,taa[2]);
		Real8 cAB   =  det3x3(taa[0],taa[1],bb);
		
		assert(det33);
		//	det33=1;
		// verif
		//	cout << " " << (taa[0][0]*cBC +  taa[1][0]*cCA + taa[2][0] * cAB)/det33 << " == " << bb[0] ;
		//	cout << " " << (taa[0][1]*cBC +  taa[1][1]*cCA + taa[2][1] * cAB)/det33 << " == " << bb[1];
		//	cout << " " << (taa[0][2]*cBC +  taa[1][2]*cCA + taa[2][2] * cAB)/det33 << " == " << bb[2] 
		//	     << "  -- " ;
		//cout << lla*sA + llb*sB+llc*sC+ (lla*llb* cAB +  llb*llc* cBC + llc*lla*cCA)/det33 
		//   << " == " << llf <<  endl;
		// computation of the gradient in the element 
		
		// H( li*lj) = grad li grad lj + grad lj grad lj
		// grad li = njk  / detT ; with i j k ={A,B,C)
		Real8 Hxx = cAB * ( nBC.x*nCA.x) +  cBC * ( nCA.x*nAB.x) + cCA * (nAB.x*nBC.x);
		Real8 Hyy = cAB * ( nBC.y*nCA.y) +  cBC * ( nCA.y*nAB.y) + cCA * (nAB.y*nBC.y);
		Real8 Hxy = cAB * ( nBC.y*nCA.x) +  cBC * ( nCA.y*nAB.x) + cCA * (nAB.y*nBC.x) 
		          + cAB * ( nBC.x*nCA.y) +  cBC * ( nCA.x*nAB.y) + cCA * (nAB.x*nBC.y);
		Real8 coef = 1.0/(3*dd*det33);
		Real8 coef2 = 2*coef;
		//	cout << " H = " << Hxx << " " << Hyy << " " <<  Hxy/2 << " coef2 = " << coef2 << endl;
		Hxx *= coef2;
		Hyy *= coef2;
		Hxy *= coef2;
		//cout << i  << " H = " << 3*Hxx/dd << " " << 3*Hyy/dd << " " <<  3*Hxy/(dd*2) << " nbb = " << nbb << endl;
		if(nbb==0)
		  {
		    dxdx[iA] += Hxx;
		    dydy[iA] += Hyy;
		    dxdy[iA] += Hxy;
		    
		    dxdx[iB] += Hxx;
		    dydy[iB] += Hyy;
		    dxdy[iB] += Hxy;
		    
		    dxdx[iC] += Hxx;
		    dydy[iC] += Hyy;
		    dxdy[iC] += Hxy;
		  }
		
	      }
	    else
	      {
		
		// if edge on boundary no contribution  => normal = 0
		if ( ! tBC || ! tBC->link ) nBC = O;
		if ( ! tCA || ! tCA->link ) nCA = O;
		if ( ! tAB || ! tAB->link ) nAB = O;
	    
		// remark we forgot a 1/2 because
		//       $\\int_{edge} w_i = 1/2 $ if $i$ is in edge 
		//                          0  if not
		// if we don't take the  boundary 
		// dxdx[iA] += ( nCA.x + nAB.x ) *Grads.x;
		
		dxdx[iA] += ( nCA.x + nAB.x ) *Grads.x;
		dxdx[iB] += ( nAB.x + nBC.x ) *Grads.x;
		dxdx[iC] += ( nBC.x + nCA.x ) *Grads.x;
		
		// warning optimization (1) the divide by 2 is done on the metrix construction
		dxdy[iA] += (( nCA.y + nAB.y ) *Grads.x + ( nCA.x + nAB.x ) *Grads.y) ;
		dxdy[iB] += (( nAB.y + nBC.y ) *Grads.x + ( nAB.x + nBC.x ) *Grads.y) ;
		dxdy[iC] += (( nBC.y + nCA.y ) *Grads.x + ( nBC.x + nCA.x ) *Grads.y) ; 
		
		dydy[iA] += ( nCA.y + nAB.y ) *Grads.y;
		dydy[iB] += ( nAB.y + nBC.y ) *Grads.y;
		dydy[iC] += ( nBC.y + nCA.y ) *Grads.y;
	      }
	    
	  } // for real all triangles 
     Int4 kk=0;
      for ( iv=0,k=0 ; iv<nbv; iv++,k+=n )
	if(Mmassxx[iv]>0) 
	  {
	    dxdx[iv] /= 2*Mmassxx[iv];
	    // warning optimization (1) on term dxdy[iv]*ci/2 
	    dxdy[iv] /= 4*Mmassxx[iv];
	    dydy[iv] /= 2*Mmassxx[iv];
	    // Compute the matrix with abs(eigen value)
	    Metric M(dxdx[iv], dxdy[iv], dydy[iv]);
	    MatVVP2x2 Vp(M);
	    //cout <<iv <<  "  M  = " <<  M <<  " aniso= " << Vp.Aniso() ;
	    Vp.Abs();
	    M = Vp;
	      dxdx[iv] = M.a11;
	      dxdy[iv] = M.a21;
	      dydy[iv] = M.a22;
	      //  cout << " (abs)  iv M  = " <<  M <<  " aniso= " << Vp.Aniso() <<endl;
	  }
	else kk++;
      
      
      // correction of second derivate
      // by a laplacien

      Real8 *d2[3] = { dxdx, dxdy, dydy};
      Real8 *dd;
      for (int xy = 0;xy<3;xy++)
	{
	  dd = d2[xy];
      // do leat 2 iteration for boundary problem
	  for (int ijacobi=0;ijacobi<Max(NbJacobi,2);ijacobi++)
	    {
	      for (i=0;i<nbt;i++) 
		if(triangles[i].link) // the real triangles 
		  {
		    // number of the 3 vertices
		    iA = Number(triangles[i][0]);
		    iB = Number(triangles[i][1]);
		    iC = Number(triangles[i][2]);
		    Real8 cc=3;
		    if(ijacobi==0)
		      cc = Max((Real8) ((Mmassxx[iA]>0)+(Mmassxx[iB]>0)+(Mmassxx[iC]>0)),1.);
		    workT[i] = (dd[iA]+dd[iB]+dd[iC])/cc;
		  }
	      for (iv=0;iv<nbv;iv++)
		workV[iv]=0;

	      for (i=0;i<nbt;i++) 
		if(triangles[i].link) // the real triangles 
		  {
		    // number of the 3 vertices
		    iA = Number(triangles[i][0]);
		    iB = Number(triangles[i][1]);
		    iC = Number(triangles[i][2]);
		    Real8 cc =  workT[i]*detT[i];
		    workV[iA] += cc;
		    workV[iB] += cc;
		    workV[iC] += cc;
		  }

	      for (iv=0;iv<nbv;iv++)
		if( ijacobi<NbJacobi || OnBoundary[iv])
		  dd[iv] = workV[iv]/(Mmass[iv]*6);
	      

	    }

	  
	}

      // constuction  of the metrix from the Hessian dxdx. dxdy,dydy

      Real8 rCutOff=CutOff*absmax;// relative cut off 

      for ( iv=0,k=0 ; iv<nbv; iv++,k+=n )
	{ // for all vertices 
	  //{
	  //Metric M(dxdx[iv], dxdy[iv], dydy[iv]);
	  // MatVVP2x2 Vp(M);	  
	  //cout << " iv M="<<  M << "  Vp = " << Vp << " aniso  " << Vp.Aniso() << endl;
	  //}
	  MetricIso Miso;
// new code to compute ci ---	  
	  Real8 ci ;
	  if (RelativeMetric)
	    { //   compute the norm of the solution
	       double xx =0,*sfk=sf+k; 
	       for (int ifield=0;ifield<nbfield;ifield++,sfk++)
	          xx += *sfk* *sfk;	       
	       xx=sqrt(xx);
	       ci = coef2/Max(xx,rCutOff);
	    }
	  else ci = cnorm;
	  
 // old 
//	  Real8 ci = RelativeMetric ? coef2/(Max(Abs(ss[k]),rCutOff)) : cnorm ;
 //   modif F Hecht 101099
	  Metric Miv(dxdx[iv]*ci, dxdy[iv]*ci,  dydy[iv]*ci);
	  MatVVP2x2 Vp(Miv);

	  Vp.Abs();
	 if(power!=1.0) 
	      Vp.pow(power);
	  


	  h1=Min(h1,Vp.lmin());
	  h2=Max(h2,Vp.lmax());

	  Vp.Maxh(hmin);
	  Vp.Minh(hmax);

	  rx = Max(rx,Vp.Aniso2());

	  Vp.BoundAniso2(coef);

	  hn1=Min(hn1,Vp.lmin());
	  hn2=Max(hn2,Vp.lmax());
	  rnx = Max(rnx,Vp.Aniso2());

	  Metric MVp(Vp);
	  vertices[iv].m.IntersectWith(MVp);
	}// for all vertices 
      if (verbosity>2)
	{ 
	  cout << "              Field " << nufield << " of solution " << nusol  << endl;
	  cout << "              before bounding :  Hmin = " << sqrt(1/h2) << " Hmax = " 
	       << sqrt(1/h1)  << " factor of anisotropy max  = " << sqrt(rx) << endl;
	  cout << "              after  bounding :  Hmin = " << sqrt(1/hn2) << " Hmax = " 
	       << sqrt(1/hn1)  << " factor of anisotropy max  = " << sqrt(rnx) << endl;
	}
	 } //  end of for all field
    }// end for all solution 

  delete [] detT;
  delete [] Mmass;
  delete [] dxdx;
  delete [] dxdy;
  delete [] dydy;
  delete []  workT;
  delete [] workV;
  delete [] Mmassxx;
  delete []  OnBoundary;
 
}


void Triangles::ReadMetric(const char * fmetrix,const Real8 hmin1=1.0e-30,const Real8 hmax1=1.0e30,const Real8 coef=1)
{
  Real8 hmin = Max(hmin1,MinimalHmin());
  Real8 hmax = Min(hmax1,MaximalHmax());
  MeshIstream f_metrix(fmetrix);
  Int4 k,j;
  f_metrix >>  k >> j ;
  if(verbosity>1)
    cout << " metrix: open " << fmetrix 
	 << ", le coef = " << coef
	 << ", hmin = " << hmin 
	 << ", hmax = " << hmax 
	 << (  (j == 1)? " Iso " : " AnIso " )<< endl;
  
  if (k != nbv || !(j == 1 || j == 3)) 
    {
      cerr << " Error Pb metrix " << k << " <> " 
	   <<  nbv << " or  1 or 3 <> " << j << endl;
      MeshError(1002);
    }
  
  cout << " j = " << j << endl;
  //  Int4 nberr = 0;
  for (Int4 iv=0;iv<nbv;iv++)
    {
      Real8 h;
      if (j == 1) 
	{
	f_metrix >>  h ;
	vertices[iv].m=Metric(Max(hmin,Min(hmax, h*coef)));
	}
      else if (j==3) 
	{
	  Real8 a,b,c;	     
	  f_metrix >>  a >> b >> c  ;
	  MetricAnIso M(a,b,c);
	  MatVVP2x2 Vp(M/coef);
	  
	  Vp.Maxh(hmin);
	  Vp.Minh(hmax);
	  vertices[iv].m = Vp;
	  
	}
    }
 
}

void Triangles::WriteMetric(ostream & f,int iso)
{
  if (iso)
    {
      f <<  nbv <<" " << 1 << endl ;
      for (Int4 iv=0;iv<nbv;iv++)
	{
	  MatVVP2x2 V=vertices[iv].m;
	  f <<  V.hmin()  << endl;
	}
    }
else
  {
    f <<  nbv <<" " << 3 << endl ;
    for (Int4 iv=0;iv<nbv;iv++)
      f <<  vertices[iv].m.a11 << " " 
	<<  vertices[iv].m.a21 << " " 
	<<  vertices[iv].m.a22 << endl;
  }
}
void  Triangles::MaxSubDivision(Real8 maxsubdiv)
{
const  Real8 maxsubdiv2 = maxsubdiv*maxsubdiv;
#ifdef DRAWING2
  inquire();
#endif	    
  if(verbosity>1)
    cout << "  -- Limit the subdivision of a edges in the new mesh by " << maxsubdiv <<   endl  ;
  // for all the edges 
  // if the len of the edge is to long 
  Int4 it,nbchange=0;    
  Real8 lmax=0;
  for (it=0;it<nbt;it++)
    {
      Triangle &t=triangles[it];
      for (int j=0;j<3;j++)
	{
	  Triangle &tt = *t.TriangleAdj(j);
	  if ( ! &tt ||  it < Number(tt) && ( tt.link || t.link)) 
	    {
		Vertex &v0 = t[VerticesOfTriangularEdge[j][0]];
		Vertex &v1 = t[VerticesOfTriangularEdge[j][1]];
		R2 AB= (R2) v1-(R2) v0;
		Metric M = v0;
		Real8 l = M(AB,AB);
		lmax = Max(lmax,l);
		if(l> maxsubdiv2)
		  { R2 AC = M.Orthogonal(AB);// the ortogonal vector of AB in M
		    Real8 lc = M(AC,AC);
		    D2xD2 Rt(AB,AC);// Rt.x = AB , Rt.y = AC;
		    D2xD2 Rt1(Rt.inv());
		    D2xD2 D(maxsubdiv2,0,0,lc);
		    D2xD2 MM = Rt1*D*Rt1.t();
#ifdef DRAWING1
		    v0.m.Draw(v0);
#endif	    
		    v0.m =  M = MetricAnIso(MM.x.x,MM.y.x,MM.y.y);
#ifdef DRAWING1
		    v0.m.Draw(v0);
#endif	    
		    //		    cout << " M(AB,AB) = " << M(AB,AB) << " == " << maxsubdiv 
		    //	 << " M(AC,AC) = " << M(AC,AC) << " == " << lc << endl; 
		    nbchange++;
		  }
		M = v1;
		l = M(AB,AB);
		lmax = Max(lmax,l);
		if(l> maxsubdiv2)
		  { R2 AC = M.Orthogonal(AB);// the ortogonal vector of AB in M
		    Real8 lc = M(AC,AC);
		    D2xD2 Rt(AB,AC);// Rt.x = AB , Rt.y = AC;
		    D2xD2 Rt1(Rt.inv());
		    D2xD2 D(maxsubdiv2,0,0,lc);
		    D2xD2  MM = Rt1*D*Rt1.t();
#ifdef DRAWING1
		    v1.m.Draw(v1);
#endif	    
		    v1.m =  M = MetricAnIso(MM.x.x,MM.y.x,MM.y.y);
#ifdef DRAWING1
		    v1.m.Draw(v1);
		    inquire();
#endif	    
		    // cout << " M(AB,AB) = " << M(AB,AB) << " == " << maxsubdiv 
		    //	 << " M(AC,AC) = " << M(AC,AC) << " == " << lc << endl; 
		    nbchange++;
		  }
		
		
	    }
	}
    }
  if(verbosity>3)
  cout << "    Nb of metric change = " << nbchange 
       << " Max  of the subdivision of a edges before change  = " << sqrt(lmax) << endl;
#ifdef DRAWING2
  inquire();
#endif	    

}

void Triangles::SmoothMetric(Real8 raisonmax) 
{ 
  if(raisonmax<1.1) return;
  if(verbosity > 1)
     cout << "  -- Triangles::SmoothMetric raisonmax = " << raisonmax << " " <<nbv <<endl;
  ReMakeTriangleContainingTheVertex();
  Int4 i,j,kch,kk,ip;
  Int4 *first_np_or_next_t0 = new Int4[nbv];
  Int4 *first_np_or_next_t1 = new Int4[nbv];
  Int4 Head0 =0,Head1=-1;
  Real8 logseuil= log(raisonmax);

  for(i=0;i<nbv-1;i++)
    first_np_or_next_t0[i]=i+1; 
  first_np_or_next_t0[nbv-1]=-1;// end;
  for(i=0;i<nbv;i++)
    first_np_or_next_t1[i]=-1;
  kk=0;
  while (Head0>=0&& kk++<100) {
    kch=0;
    for (i=Head0;i>=0;i=first_np_or_next_t0[ip=i],first_np_or_next_t0[ip]=-1)
      {  //  pour tous les triangles autour du sommet s
	// 	cout << kk << " i = " << i << " " << ip << endl;
	register Triangle * t= vertices[i].t;
	assert(t);
	Vertex & vi = vertices[i];
	TriangleAdjacent ta(t,EdgesVertexTriangle[vertices[i].vint][0]);
	Vertex *pvj0 = ta.EdgeVertex(0);
	while (1) {
	  //	  cout << i << " " <<  Number(ta.EdgeVertex(0)) << " "
	  //      << Number(ta.EdgeVertex(1)) << "  ---> " ;
	  ta=Previous(Adj(ta));
	  // cout <<  Number(ta.EdgeVertex(0)) << " " << Number(ta.EdgeVertex(1)) << endl;
	  assert(vertices+i == ta.EdgeVertex(1));
	  Vertex & vj = *(ta.EdgeVertex(0));
	  if ( &vj ) {
	    j= &vj-vertices;
	    assert(j>=0 && j < nbv);
	    R2 Aij = (R2) vj - (R2) vi;
	    Real8 ll =  Norme2(Aij);
	    if (0) {  
	      Real8 hi = ll/vi.m(Aij);
	      Real8 hj = ll/vj.m(Aij);
	      if(hi < hj)
		{
		  Real8 dh=(hj-hi)/ll;
		  //cout << " dh = " << dh << endl;
		  if (dh>logseuil) {
		    vj.m.IntersectWith(vi.m/(1 +logseuil*ll/hi));
		    if(first_np_or_next_t1[j]<0)
		      kch++,first_np_or_next_t1[j]=Head1,Head1=j;
		  }
		}
	    } 
	    else
	      {
		Real8 li = vi.m(Aij);
		//Real8 lj = vj.m(Aij);
		//		if ( i == 2 || j == 2)
		//  cout << " inter " << i << " " << j << " " << ((1 +logseuil*li)) <<  endl;
	      	if( vj.m.IntersectWith(vi.m/(1 +logseuil*li)) )
		  //if( vj.m.IntersectWith(vi.m*(lj/li/(1 +logseuil*lj))) )
		  if(first_np_or_next_t1[j]<0) // if the metrix change 
		    kch++,first_np_or_next_t1[j]=Head1,Head1=j;
	      }
	  }
	  if  ( &vj ==  pvj0 ) break;
	}
      }
    Head0 = Head1;
    Head1 = -1;
    Exchange(first_np_or_next_t0,first_np_or_next_t1);
    if(verbosity>5)
    cout << "     Iteration " << kk << " Nb de  vertices with change  " << kch << endl;
  }
  if(verbosity>2 && verbosity < 5) 
    cout << "    Nb of Loop " << kch << endl;
  delete [] first_np_or_next_t0;
  delete [] first_np_or_next_t1;
}

void Geometry::ReadMetric(const char * fmetrix,Real8 hmin=1.0e-30,Real8 hmax=1.0e30,Real8 coef=1)
{
  hmin = Max(hmin,MinimalHmin());
  MeshIstream f_metrix(fmetrix);
  Int4 k,j;
  f_metrix >>  k >> j ;
  if(verbosity>1)
    cout << "  -- ReadMetric  " << fmetrix 
	 << ",  coef = " << coef
	 << ", hmin = " << hmin 
	 << ", hmax = " << hmax 
	 << (  (j == 1)? " Iso " : " AnIso " ) << endl;
  
  if (k != nbv ||  !(j == 1 || j == 3)) {
    cerr << " Error Pb metrix " << k << " <> " 
	 <<  nbv << " or  1 or 3  <> " << j << endl;
    MeshError(1003);}
	 
  
  //  Int4 nberr = 0;
  for (Int4 iv=0;iv<nbv;iv++)
    {
    Real8 h;
    if (j == 1) 
      {
      f_metrix >>  h ;
      vertices[iv].m=Metric(Max(hmin,Min(hmax, h*coef)));
      }
    else if (j==3) 
      {
	Real8 a,b,c;	     
	f_metrix >>  a >> b >> c  ;
	MetricAnIso M(a,b,c);
	MatVVP2x2 Vp(M/coef);
	Vp.Maxh(hmin);
      Vp.Minh(hmax);
      vertices[iv].m = Vp;
      }
    }
  
}


Real8 LengthInterpole(const MetricAnIso Ma,const  MetricAnIso Mb, R2 AB)
{
  Real8 k=1./2.;
  int level=0;
  static int kkk=0;
  static  Metric Ms1[32],Ms2[32];
  static Real8 lMs1[32],lMs2[32];
  static double K[32];
  Real8 l=0,sss=0;
  Ms1[level]=Ma;
  Ms2[level]=Mb;
  Real8 sa =  Ma(AB);
  Real8 sb =  Mb(AB);
  lMs1[level]=sa;
  lMs2[level]=sb;
  K[level]=k;
  level++;
  int i=0;
  Real8 * L= LastMetricInterpole.L, *S = LastMetricInterpole.S;
  Real8  sstop = 0.1; // Max(0.6,(sa+sb)/5000);
  while (level) {
    level--;
    Metric M1=Ms1[level];
    Metric M2=Ms2[level];
    k=K[level];
    Real8 s1=  lMs1[level];
    Real8 s2=  lMs2[level];

    Real8 s= (s1+s2)*k;
//    if (level >20  && i < 2030-level)
//    cout << "                  level " << level << " " << i << " " << s << " " << k <<endl;
    if( s > sstop   && level < 30 && i < 500-level ) {
      Metric Mi(0.5,M1,0.5,M2);
      Real8 si = Mi(AB);
      if( Abs((s1+s2)-(si+si)) > s1*0.001) 
	{
	  k=k/2;
	  // we begin by the end to walk in the correct sens from a to b
	       // due to the stack 
	       Ms1[level]=Mi;
	  Ms2[level]=M2;
	  lMs1[level]=si;
	  lMs2[level]=s2;
	  K[level]=k;
	  level++;
	  Ms1[level]=M1;
	  Ms2[level]=Mi;
	  lMs1[level]=s1;
	  lMs2[level]=si;
	  K[level]=k;
	  level++;
	}
      else
	L[i]= l += s,S[i]=sss+=k,i++;
    }
    else 
      L[i]= l += s,S[i]=sss+=k,i++;//cout << i << " l = " << l << " sss = " << sss << endl;
  }
  // warning for optimisation S is in [0:0.5] not in [0:1]
  assert(i<512);
  LastMetricInterpole.lab=l;
  LastMetricInterpole.opt=i;
  if (i>200 && kkk++<10)
     cout << "Warning LengthInterpole: ( i = " << i << " l = " << l << " sss " << sss << " ) " << sstop <<endl;
  return l;
}

Real8 abscisseInterpole(const MetricAnIso Ma,const  MetricAnIso Mb, R2 AB,Real8 s,int optim)
{ 
  if(!optim)  LengthInterpole(Ma,Mb,AB);
  Real8 l  = s* LastMetricInterpole.lab,r;
  int j=LastMetricInterpole.opt-1,i=0,k;
  
  Real8 * L= LastMetricInterpole.L, *S = LastMetricInterpole.S;
  // warning for optimisation S is the abcisse in [0:0.5]
  // and L is le lenght 
  if(l<=L[0])
    r=2*S[0]*l/L[0];
  else if (l>=L[j])
    r=1;
  else 
    {
      while (j-i>1)
	{
	  k= (i+j)/2;
	  if(l<=L[k])
	    j=k;// l<=L[j] 
	  else
	    i=k; //  L[i]<l
	};
      //   cout  << i << " " << j  <<" " << L[i] << " " << L[j] << " " << S[i] << " " <<  S[j]  << " l=" << l << endl;
      if (i==j)
	r = 2*S[i];
      else
	r =  2*(S[i]*(L[j]-l)+ S[j]*(l-L[i]))/(L[j]-L[i]);
    }
  assert(r<=1 && r>=0);
  return r ;
    
}


#ifdef DRAWING 

void MetricAnIso::Draw(R2 c) const 
{ 
  float x= c.x,y= c.y;
  if (InPtScreen(x,y)) {
  R2 X(cos(0.0),sin(0.0));
  X = X / operator()(X);
  rmoveto(x+X.x,y+X.y);
  
  for (int i=1;i<=100;i++)
   { double t= 2*Pi*i/100.0;
     R2 X(cos(t),sin(t));
     X = X / Max(operator()(X),1.0e-5);
     rlineto(x+X.x,y+X.y); } 
  }
}

void MetricIso::Draw(R2 c) const 
{ 
  float x= c.x,y= c.y;
  if (InPtScreen(x,y)) {
  rmoveto(x+h,y);
  for (int i=1;i<=40;i++)
   { double t= Pi*i/20.0;
     rlineto(x+h*cos(t),y+h*sin(t)); } 
  }
}


#endif
}   // end of namespace bamg 
