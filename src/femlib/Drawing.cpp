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

#include <cmath>
#include "error.hpp"
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "rgraph.hpp"


#include "RNM.hpp"
#include "fem.hpp"
#include "FESpacen.hpp" 
#include "FESpace.hpp" 

namespace Fem2D {
void NewSetColorTable(int nb,float *colors=0,int nbcolors=0,bool hsv=true);

void NewSetColorTable(int nb,float *colors,int nbcolors,bool hsv)
{
  if(colors && nbcolors)
    SetColorTable1(nb,hsv,nbcolors,colors);
  else
    SetColorTable(nb);
}
int dichotomie(const RN_ &viso,R v) 
{
  int i=0,j=viso.N(),k;
  if  (v <viso[0] || v >viso[j-1]) 
    return -1;  
  while (i<j-1)    
   if ( viso[k=(i+j)/2]> v) j=k;
   else i=k;
  return i;
}

void plot(long i)
{
  char buf[24];
  snprintf(buf,24,"%ld",i);
  plotstring(buf);
}
void plot(int i)
{
  char buf[24];
  snprintf(buf,24,"%d",i);
  plotstring(buf);
}
void plot(double i)
{
  char buf[24];
  snprintf(buf,24,"%g",i);
  plotstring(buf);
}

void DrawIsoT(const R2 Pt[3],const R ff[3],const RN_ & Viso);
void DrawIsoTfill(const R2 Pt[3],const R ff[3],const RN_ & Viso,double rapz);

void FillRect(float x0,float y0, float x1, float y1)
 {
     float r[8];
     r[0]=x0;r[1]=y0;
     r[2]=x1;r[3]=y0;
     r[4]=x1;r[5]=y1;
     r[6]=x0;r[7]=y1;
     fillpoly(4,r);
 }

void PlotValue(const RN_ & Viso,int  k0,const char * cmm)
{
   float xmin,xmax,ymin,ymax;
 //  cout << "PlotValue " << Viso << endl;
   // int ix,iy;
  // GetSizeScreen(ix,iy);   
   getcadre(xmin,xmax,ymin,ymax);
   float dx=(xmax-xmin);
   float dy=(ymax-ymin);
   //  10 points  
  // int kk = Max(30,iy/10);
   float h=GetHeigthFont();
   float x0=xmin+dx*0.8;
   float y= ymin+dy*0.95;
  // cout << x0 << " " << y << " " << h <<  endl;
   y -= k0* h*1.4;
   couleur(0);
   FillRect(x0-h*0.5,y-h*(1.4*Viso.N()+0.3),x0+h*9,y+h*1.5);
   couleur(1);
   rmoveto(x0+1.2*h,y);
   plotstring(cmm);
   y -=  h*1.4;
   for (int i=0;i<Viso.N();i++)
    {
       couleur(i+4);
       FillRect(x0,y,x0+h,y+h);
       rmoveto(x0+1.2*h,y);
       y -=  h*1.4;
       plot(Viso[i]);                 
    }
}

void DrawCommentaire(const char * cm,float x,float y) 
{
   float xmin,xmax,ymin,ymax;
   getcadre(xmin,xmax,ymin,ymax);
   float dx=(xmax-xmin);
   float dy=(ymax-ymin);
   
   rmoveto(xmin+dx*x,ymin+dy*y);
   plotstring(cm);   
}

void Mesh::DrawBoundary() const
{  
    for(int i=0;i<neb;i++)
       bedges[i].Draw();
}
void Mesh::InitDraw() const 
   {
     R2 Pmin,Pmax;
     BoundingBox(Pmin,Pmax);
     R2 O((Pmin+Pmax)/2);
     R r =  (Max(Pmax.x-Pmin.x,Pmax.y-Pmin.y)*0.55);
     showgraphic();
     cadreortho((float)O.x,(float)(O.y+r*0.05),(float) r);
    }
inline  R Square(R x){return x*x;}

void  SetDefaultIsoValue(const RN_& U,RN_ & Viso)
 {
   R umn = U.min();
   R umx = U.max();
   int N = Viso.N();
   R d = (umx-umn)/(N+1);
   R x = umn+d/2;
   for (int i = 0;i < N;i++)
     {Viso[i]=x;x +=d; }
   NewSetColorTable(N+4) ;  
 }
void  SetDefaultIsoValue(const RN_& u,const RN_& v,RN_ & Viso)
 {
   R s = (Square(u[0])+Square(v[0]));
   R umn =s ;
   R umx = s;    
   int n(v.N());
    for (int i=1;i<n;i++)
     {   s = (Square(u[i])+Square(v[i]));
        umn = Min(umn,s),umx=Max(umn,s);
      }  
   umx=sqrt(umx);
   umn=sqrt(umn);
   
   int N = Viso.N();
   R d = (umx-umn)/(N-1);
   R x = umn;
   for (int i = 0;i < N;i++)
     {Viso[i]=x;x +=d; }
   Viso[N-1]= umx; 
   NewSetColorTable(N+4) ;  
 }
 
 void Mesh::Draw(int init,bool fill) const
 {
   
   if (init>0)  InitDraw();
   
   SetColorTable(16) ;  
   for (int i=0;i<nt;i++)
    {
      if (fill) 
       {
         triangles[i].Fill(2+Abs(triangles[i].lab));         
         couleur(triangles[i].lab?1:0);
         triangles[i].Draw();
       }
      else 
       {
       couleur(1+Abs(triangles[i].lab));
       triangles[i].Draw();}
    }
     

   couleur(1);
   penthickness(2); 
   for (int i=0;i<NbMortars;i++)
    {     
      mortars[i].Draw();
    }
   penthickness(3);
   double l= 1e100;
   for (int i=0;i<neb;i++)
     {
       couleur(1+Abs( bedges[i].lab));
       bedges[i].Draw();
       l=Min(l,bedges[i].length());
     }
   for (int i=0;i<nv;i++)
     { 
       R2 P=vertices[i], N=vertices[i].Ne();
       if (Norme2_2(N) ) 
         { MoveTo(P); LineTo(P+N*l/5);  }
     }
   penthickness(1); 
   if(init%2)
   for (int i=0;i<nv;i++)
    {
      couleur(1+Abs(vertices[i].lab));
      MoveTo(vertices[i]);
      plot(i);
    }
   couleur(1);
//   rattente(0);
 }
    // to remove link of blas with graphic tools dec. 2018 FH
    double dot(const RN_& U,const RN_& V)
    {
        double s=0;
        for(int i=0; i<U.N(); ++i)
            s += U[i]*V[i];
        return s;
    }
void FElement::Draw(const RN_& U,const RN_ & Viso,int composante) const
{   
  int nsb = nbsubdivision();
  int nsb2 = NbOfSubTriangle(nsb);
  int nbdf=NbDoF();
  RN fk(nbdf);
  for (int i=0;i<nbdf;i++) // get the local value
    fk[i] = U[operator()(i)];
  RNMK fb(nbdf,N,3); //  the value for basic fonction
  RN ff(3);
  R2 Pt,P[3],A(T[0]),B(T[1]),C(T[2]);
  bool whatd[last_operatortype];
  initwhatd(whatd,0);
  
  for (int k=0;k<nsb2;k++)
   {   
   //  ff=0.0;
     for (int j=0;j<3;j++)
      { 
     //   cout << nbdf << endl;
        Pt=SubTriangle(nsb,k,j); //  current point 
        BF(whatd,Pt,fb);
        ff[j] = dot(fb('.',composante,0),fk);
        P[j] = A*(1-Pt.x-Pt.y)+B*Pt.x + C*Pt.y;
       }
     DrawIsoT(P,ff,Viso);  
   } 
} 
    
void FElement::Drawfill(const RN_& U,const RN_ & Viso,int composante, double rapz) const
{   
  int nsb = nbsubdivision();
  int nsb2 =  NbOfSubTriangle(nsb);
  int nbdf=NbDoF();
  RN fk(nbdf);
  for (int i=0;i<nbdf;i++) // get the local value
    fk[i] = U[operator()(i)];
  RNMK fb(nbdf,N,3); //  the value for basic fonction
  RN ff(3);
  R2 Pt,P[3],A(T[0]),B(T[1]),C(T[2]);
  bool whatd[last_operatortype];
  initwhatd(whatd,0);
  
  for (int k=0;k<nsb2;k++)
   {   
   //  ff=0.0;
     for (int j=0;j<3;j++)
      { 
     //   cout << nbdf << endl;
        Pt=SubTriangle(nsb,k,j); //  current point 
        BF(whatd,Pt,fb);
        ff[j] = dot(fb('.',composante,0),fk);
        P[j] = A*(1-Pt.x-Pt.y)+B*Pt.x + C*Pt.y;
       }
     DrawIsoTfill(P,ff,Viso,rapz);  
   } 
} 

R2 FElement::MinMax(const RN_& U,const RN_& V,int i0,int i1) const
{
  R2 minmax(1e100,-1e100);
  int nsb = nbsubdivision();
  int nsb2 =  NbOfSubTriangle(nsb);
  int nbdf=NbDoF();
  RN fk(nbdf);
  RN gk(nbdf);
  RNMK fb(nbdf,N,3); //  the value for basic fonction
  RN f0(3);
  RN f1(3);
  bool whatd[last_operatortype];
  initwhatd(whatd,0);  
  R2 Pt,A(T[0]),B(T[1]),C(T[2]);
  for (int i=0;i<nbdf;i++) // get the local value
    fk[i] = U[operator()(i)];
  for (int i=0;i<nbdf;i++) // get the local value
    gk[i] = V[operator()(i)];
  for (int k=0;k<nsb2;k++)
   {   
   //  ff=0.0;
     for (int j=0;j<3;j++)
      { 
     //   cout << nbdf << endl;
        R2 Pt=SubTriangle(nsb,k,j); //  current point 
        BF(whatd,Pt,fb);
        f0[j] = dot(fb('.',i0,0),fk);
        f1[j] = dot(fb('.',i1,0),gk);
        R2 uv(f0[j],f1[j]);
        R uv2 = (uv,uv);
        minmax.x=Min(minmax.x,uv2);
        minmax.y=Max(minmax.y,uv2);
        
       }
    }
   return minmax;
}

R2 FElement::MinMax(const RN_& U,int i0) const
{
  R2 minmax(1e100,-1e100);
  int nsb = nbsubdivision();
  int nsb2 =  NbOfSubTriangle(nsb);
  int nbdf=NbDoF();
  RN fk(nbdf);
  RNMK fb(nbdf,N,3); //  the value for basic fonction
  RN f0(3);
  RN f1(3);
  R2 Pt,A(T[0]),B(T[1]),C(T[2]);
   bool whatd[last_operatortype];
  initwhatd(whatd,0);
 
  for (int i=0;i<nbdf;i++) // get the local value
    fk[i] = U[operator()(i)];
  for (int k=0;k<nsb2;k++)
   {   
   //  ff=0.0;
     for (int j=0;j<3;j++)
      { 
     //   cout << nbdf << endl;
        R2 Pt=SubTriangle(nsb,k,j); //  current point 
        BF(whatd,Pt,fb);
        f0[j] = dot(fb('.',i0,0),fk);
        minmax.x=Min(minmax.x,f0[j]);
        minmax.y=Max(minmax.y,f0[j]);        
       }
    }
   return minmax;
}
void FElement::Draw(const RN_& U,const RN_& V,const RN_ & Viso,R coef,int i0,int i1,double ArrowSize) const
{
  bool whatd[last_operatortype];
  initwhatd(whatd,0);
   
  float xmin,xmax,ymin,ymax;
  getcadre(xmin,xmax,ymin,ymax);
  float d= Max(ymax-ymin,xmax-xmin);
  R kk = d*0.005;
  if(ArrowSize>0) kk=ArrowSize/100.;
  R cc = d*0.05;
  int nsb = nbsubdivision();
  int nsb2 = NbOfSubTriangle(nsb);
  int nbdf=NbDoF();
  RN fk(nbdf);
  RN gk(nbdf);
  for (int i=0;i<nbdf;i++) // get the local value
    fk[i] = U[operator()(i)];
  for (int i=0;i<nbdf;i++) // get the local value
    gk[i] = V[operator()(i)];
    
  RNMK fb(nbdf,N,3); //  the value for basic fonction
  RN f0(3);
  RN f1(3);
  R2 Pt,P[3],A(T[0]),B(T[1]),C(T[2]);
  for (int k=0;k<nsb2;k++)
   {   
   //  ff=0.0;
     for (int j=0;j<3;j++)
      { 
     //   cout << nbdf << endl;
        Pt=SubTriangle(nsb,k,j); //  current point 
        BF(whatd,Pt,fb);
        f0[j] = dot(fb('.',i0,0),fk);
        f1[j] = dot(fb('.',i1,0),gk);
     //   if(number<2) 
     //   cout << number << " " << fk << " " << fb('.',0,0) << " :  " << f0[j] << " " <<  f1[j] << endl;
        P[j] = A*(1-Pt.x-Pt.y)+B*Pt.x + C*Pt.y;
        R2 uv(f0[j],f1[j]);
        R  l = Max(sqrt((uv,uv)),1e-30) ;
        
         {
             int col = 4+dichotomie(Viso,l);
             if(verbosity>99) cout << " color vect " << l << " " << col << endl; 
          couleur(col);
          uv = coef*uv;
          l *= coef;
          R2 dd = uv*(-kk/l);// modif F.H size of arraow 08/14 FH.l
          
          R2 dn = dd.perp()*0.5;
          if (l*10000.< kk) continue;
          if (l < kk) 
	        uv = uv*(kk/l);
          else if (l> cc)
	        uv = uv*(cc/l);	   
          MoveTo(P[j]);       
          LineTo(P[j]+uv);
          if (l>kk) {
            LineTo(P[j]+uv+dd+dn);
            MoveTo(P[j]+uv+dd-dn);
            LineTo(P[j]+uv);}
          }
      }
   } 
} 

void DrawIsoT(const R2 Pt[3],const R ff[3],const RN_ & Viso)
{
  R2 PQ[5];
  int NbIso = Viso.N();
  for(int l=0;l< NbIso;l++)  /*    loop on the level curves */
    {
      R xf = Viso[l];
      int im=0;
      for(int i=0;i<3;i++) // for the  3 edges 
	{
	  int j = (i+1)%3;
	  R fi=(ff[i]);
	  R fj=(ff[j]);
	  
	  if(((fi<=xf)&&(fj>=xf))||((fi>=xf)&&(fj<=xf)))
	    {
	      if (Abs(fi-fj)<=0.1e-10) 	/* one side must be drawn */
		{
		  couleur(l+4);
		  MoveTo(Pt[i]);
		  LineTo(Pt[j]);
		}
	      else
		{
		  R  xlam=(fi-xf)/(fi-fj);
		  PQ[im++]   = Pt[i] * (1.F-xlam)  +  Pt[j]* xlam;
		}
	    }
	}
      
      if (im>=2) /*    draw one segment */
	{
	  couleur(l+4);
	  MoveTo(PQ[0]);
	  LineTo(PQ[1]);
	}
    }
  
} 
void DrawIsoTfill(const R2 Pt[3],const R ff[3],const RN_ & Viso,double rapz)
{
  R2 PQ[10];
  int NbIso = Viso.N();
  R eps= (Viso[NbIso-1]-Viso[0])*1e-6;
  for(int l=1;l< NbIso;l++)  //   loop on the level curves 
    {
      R xfb = Viso[l-1];
      R xfh = Viso[l];
      int im=0;
      for(int i=0;i<3;i++) // for the  3 edges 
	{
          int j=(i+1)%3;
          R fi=(ff[i]);
          R fj=(ff[j]);
          R xxfb =  xfb;
          R xxfh =  xfh;
	  if (fj<fi ) Exchange(xxfb,xxfh);
          R xf  = xxfb;
	  if(((fi<=xf)&&(fj>=xf))||((fi>=xf)&&(fj<=xf)))
	    {
	      if (Abs(fi-fj)>=0.1e-20)
		{
		  R  xlam=(fi-xf)/(fi-fj);
		  PQ[im++]   = Pt[i] * (1.F-xlam)  +  Pt[j]* xlam;
		}
	    }
          xf = xxfh;	  
	  if(((fi<=xf)&&(fj>=xf))||((fi>=xf)&&(fj<=xf)))
	    {
	      if (Abs(fi-fj)>=0.1e-20)
		{
		  R  xlam=(fi-xf)/(fi-fj);
		  PQ[im++]   = Pt[i] * (1.F-xlam)  +  Pt[j]* xlam;
		}
	    }
	   if (  xfb-eps <=fj  && fj <= xfh+eps) 
	     PQ[im++] = Pt[j];
	}
       if (im>2) 
         {
           float f[20];
           int k=0;
           for (int i=0;i<im;i++)
            f[k++]=PQ[i].x,f[k++]=PQ[i].y;
           f[k++]=f[0];
           f[k++]=f[1];
           
           couleur(3+l);
           fillpoly(im,f);
         }
     }
} 

template<class R2>
void TBoundaryEdge<R2>::Draw() const 
{
  couleur(2+(lab%6));
  MoveTo(*vertices[0]);
  LineTo(*vertices[1]);
}

    template<class R>   
    void FElement::SaveDraw(const KN_<R> & U,int composante,R* Usave) const 
    {
	int nsb = nbsubdivision();
	int nsbv = NbOfSubInternalVertices(nsb);
	int nbdf=NbDoF();
	KN<R> fk(nbdf);
	for (int i=0;i<nbdf;i++) // get the local value
	    fk[i] = U[operator()(i)];
	R2 Pt;
	bool whatd[last_operatortype];
	initwhatd(whatd,0);
	
	RNMK fb(nbdf,N,3); //  the value for basic fonction
	for (int k=0;k<nsbv;k++)
	  {
	      Pt=SubInternalVertex(nsb,k); 
	      BF(whatd,Pt,fb);
	    Usave[k] = R();// (fb('.',composante,0),fk);
	    for( int i=0; i<nbdf;++i)
		Usave[k]   +=  fb(i,composante,0)*fk[i];
	    
	  }
	
    }
     template<class R> 
    void FElement::SaveDraw(const KN_<R> & U,const KN_<R> & V,int iU,int iV,R * Usave) const 
    {
	int nsb = nbsubdivision();
	int nsbv = NbOfSubInternalVertices(nsb);
	int nbdf=NbDoF();
	KN<R> fk(nbdf);
	KN<R> gk(nbdf);
	for (int i=0;i<nbdf;i++) // get the local valu
	  {
	      fk[i] = U[operator()(i)];
	      gk[i] = V[operator()(i)];
	  }
	R2 Pt;
	bool whatd[last_operatortype];
	initwhatd(whatd,0);
	
	RNMK fb(nbdf,N,3); //  the value for basic fonction
	for (int k=0,k2=0;k<nsbv;k++)
	  {
	      Pt=SubInternalVertex(nsb,k); 
	      BF(whatd,Pt,fb);
	    Usave[k2]=R();
	    Usave[k2+1]=R();
	    for( int i=0; i<nbdf;++i)
	      {
		Usave[k2]   +=  fb(i,iU,0)*fk[i];
		Usave[k2+1] +=  fb(i,iV,0)*gk[i];
	      }
	    k2+=2;
	  }
	
    } 
    
template<class R>    
KN<R>  FESpace::newSaveDraw(const KN_<R> & U,int composante,int & lg,int & nsb) const 
    {
        nsb = TFE[0]->nbsubdivision; 
	int nsbv = NbOfSubInternalVertices(nsb);
	lg = nsbv*Th.nt;
	if(verbosity>99)
	cout << "           ++ newSaveDraw what: nt " << Th.nt << " " << nsbv << " " << lg << endl;
	KN<R> v(lg);
	ffassert(v);
        for (int k=0,i=0;k<Th.nt;k++)
	  {
	    (*this)[k].SaveDraw( U,composante,&v[i]);	
	      i+=nsbv;
	  }
	return KN<R>(true,v);// to remove the copy.
    }
    
    
template<class R>      
KN<R>  FESpace::newSaveDraw(const KN_<R> & U,const KN_<R> & V,int iU,int iV,int & lg,int & nsb) const 
    {
	nsb = TFE[0]->nbsubdivision;
	int nsbv = NbOfSubInternalVertices(nsb)*2;
	lg = nsbv*Th.nt;
	
	KN<R> v(lg);
	
        for (int k=0,i=0;k<Th.nt;k++)
	  {
	      (*this)[k].SaveDraw( U,V,iU,iV,&v[i]);	
	      i+=nsbv;
	  }
	return  KN<R>(true,v);// to remove the copy.
    }
    
typedef complex<double> Complex;   
template     
    KN<double>  FESpace::newSaveDraw<double>(const KN_<double> & U,const KN_<double> & V,int iU,int iV,int & lg,int & nsb) const ;
template    
    KN<double>  FESpace::newSaveDraw<double>(const KN_<double> & U,int composante,int & lg,int & nsb) const ;
    template     
    KN<Complex>  FESpace::newSaveDraw<Complex>(const KN_<Complex> & U,const KN_<Complex> & V,int iU,int iV,int & lg,int & nsb) const ;
    template    
    KN<Complex>  FESpace::newSaveDraw<Complex>(const KN_<Complex> & U,int composante,int & lg,int & nsb) const ;
    
   
void  FESpace::Draw(const RN_& U,const RN_ & Viso,int j,float *colors,int nbcolors,bool hsv,bool drawborder) const 
{
  showgraphic();
  NewSetColorTable(Viso.N()+4,colors,nbcolors,hsv);
  for (int k=0;k<Th.nt;k++) 
    (*this)[k].Draw( U,Viso,j);
  NewSetColorTable(2+6,colors,nbcolors,hsv);
  if(drawborder) Th.DrawBoundary();
  NewSetColorTable(Viso.N()+4,colors,nbcolors,hsv);

}

void  FESpace::Drawfill(const RN_& U,const RN_ & Viso,int j,double rapz,float *colors,int nbcolors,bool hsv,bool drawborder) const 
{
  showgraphic();
  NewSetColorTable(Viso.N()+4,colors,nbcolors,hsv);
  for (int k=0;k<Th.nt;k++) 
    (*this)[k].Drawfill( U,Viso,j,rapz);
  NewSetColorTable(2+6,colors,nbcolors,hsv);
  if(drawborder) Th.DrawBoundary();
  NewSetColorTable(Viso.N()+4,colors,nbcolors,hsv);
}



R2 FESpace::MinMax(const KN_<R>& U,const KN_<R>& V,int j0,int j1,bool  bb) const 
{
  R2 Pminmax(1e100,-1e100);
    for (int k=0;k<Th.nt;k++)
     { 
        const Triangle & K(Th[k]);
        R2 A(K[0]),B(K[1]),C(K[2]);
        if (bb  || InRecScreen(Min(A.x,B.x,C.x),Min(A.y,B.y,C.y),Max(A.x,B.x,C.x),Max(A.y,B.y,C.y))) 
         Pminmax=minmax(Pminmax,(*this)[k].MinMax(U,V,j0,j1));
     }
   return Pminmax;
}
R2 FESpace::MinMax(const KN_<R>& U,int j0,bool bb) const 
{
  R2 Pminmax(1e100,-1e100);
    for (int k=0;k<Th.nt;k++) 
      {
        const Triangle & K(Th[k]);
        R2 A(K[0]),B(K[1]),C(K[2]);
        if (bb  || InRecScreen(Min(A.x,B.x,C.x),Min(A.y,B.y,C.y),Max(A.x,B.x,C.x),Max(A.y,B.y,C.y))) 
          Pminmax=minmax(Pminmax,(*this)[k].MinMax(U,j0));
      }
   return Pminmax;
}
void  FESpace::Draw(const KN_<R>& U,const RN_ & Viso, R coef,int j0,int j1,float *colors,int nbcolors,bool hsv,bool drawborder,double ArrowSize) const
{ 
  showgraphic();
  NewSetColorTable(Viso.N()+5,colors,nbcolors,hsv);
  for (int k=0;k<Th.nt;k++) 
    (*this)[k].Draw( U,U,Viso,coef,j0,j1,ArrowSize);
  NewSetColorTable(2+6,colors,nbcolors,hsv);
  if(drawborder) Th.DrawBoundary();
  NewSetColorTable(Viso.N()+5,colors,nbcolors,hsv);

}

void  FESpace::Draw(const KN_<R>& U,const KN_<R>& V,const RN_ & Viso, R coef,int iu,int iv,float *colors,int nbcolors,bool hsv,bool drawborder,double ArrowSize) const
{ 
  showgraphic();
  NewSetColorTable(Viso.N()+5,colors,nbcolors,hsv);
  for (int k=0;k<Th.nt;k++) 
    (*this)[k].Draw( U,V,Viso,coef,iu,iv,ArrowSize);
  NewSetColorTable(2+6,colors,nbcolors,hsv);
  if(drawborder) Th.DrawBoundary();
   NewSetColorTable(Viso.N()+5,colors,nbcolors,hsv);

}
   

template<class R2>
void TTriangle<R2>::Draw(double skrink) const
{
  const TTriangle<R2> & K(*this);
  R2 A(K[0]),B(K[1]),C(K[2]),G((A+B+C)/3.);
  A = G + (A-G)*skrink;
  B = G + (B-G)*skrink;
  C = G + (C-G)*skrink;
  MoveTo(A);
  LineTo(B);
  LineTo(C);
  LineTo(A);
}

template<class R2>
void TTriangle<R2>::Draw(int edge,double skrink) const 
{
  const TTriangle<R2> & K(*this);
  R2 A(K[0]),B(K[1]),C(K[2]),G((A+B+C)/3.);
  MoveTo(G+(*vertices[(edge+1)%3]-G)*skrink);
  LineTo(G+(*vertices[(edge+2)%3]-G)*skrink);
}


template  void TTriangle<R2>::Draw(double skrink) const;
template  void TTriangle<R2>::Draw(int edge,double skrink) const;

template<class R2>
void TTriangle<R2>::Fill(int color) const 
{
  const TTriangle<R2> & K(*this);
  R2 A(K[0]),B(K[1]),C(K[2]);
  float p[]={(float)A.x,(float)A.y,(float)B.x,(float)B.y,(float)C.x,(float)C.y};
  int c=LaCouleur(); 
  couleur(color); 
  fillpoly(3,p);
  couleur(c); 
}

void DrawMark(R2 P,R k)
 {
   float x0,x1,y0,y1;
   getcadre(x0,x1,y0,y1);
   float h= (x1-x0)*(float)k;
   rmoveto((float)P.x+h,(float)P.y-h);
   rlineto((float)P.x+h,(float)P.y+h);
   rlineto((float)P.x-h,(float)P.y+h);
   rlineto((float)P.x-h,(float)P.y-h);
   rlineto((float)P.x+h,(float)P.y-h);
 }
 
 Triangle * Mesh::Find(const R2 & P) const
 {
    // brute force
    
    for (int i=0;i<nt;i++)
     {
      kthrough++;
      const Triangle & K(triangles[i]);
      R2 A(K[0]),B(K[1]),C(K[2]);
      R a=Area2(P,B,C);
      R b=Area2(A,P,C);
      R c=Area2(A,B,P);
      R s=a+b+c;
      R eps=s*1e-6;
      if (a>-eps && b >-eps && c >-eps) return triangles + i;
      
     }
   return 0; // outside 
 }
    
  

template<class R2>
 void TMortar<R2>::Draw() const {
   throwassert(Th);
   for (int i=0;i<nleft;i++)
     (*Th)[left[i]/3].Draw(left[i]%3,0.8);     
   for (int i=0;i<nright;i++)
     (*Th)[right[i]/3].Draw(right[i]%3,0.8);
     
 }
}
