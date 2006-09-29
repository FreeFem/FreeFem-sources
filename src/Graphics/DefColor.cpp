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
using namespace std;

void DefColor(float & r, float & g, float & b,
              int k,int nb, bool hsv,bool grey,int nbcolors,float *colors);


 void hsvToRgb (float h, float s, float v, float & r, float & g, float & b)
{
  int i;
  float aa, bb, cc, f;

  if (s == 0) /* Grayscale */
    r = g = b = v;
  else {
    h = h - floor(h);
    if (h == 1.0) h = 0;
    h *= 6.0;
    i =  int(h);
    f = h - i;
    aa = v * (1 - s);
    bb = v * (1 - (s * f));
    cc = v * (1 - (s * (1 - f)));
    switch (i) {
    case 0: r = v;  g = cc; b = aa; break;
    case 1: r = bb; g = v;  b = aa; break;
    case 2: r = aa; g = v;  b = cc; break;
    case 3: r = aa; g = bb; b = v;  break;
    case 4: r = cc; g = aa; b = v;  break;
    case 5: r = v;  g = aa; b = bb; break;
    }
  }
}
//  def des couleurs de la tables 
void DefColor(float & r, float & g, float & b,
              int k,int nb, bool hsv,bool grey,int nbcolors,float *colors)
{
  if(k<=0) {  r=g=b=1.;} //  white
  else if (k==1)  { r=g=b=0.; } // black
  else if (k >= nb)   {  r=g=b=0.;} // black
  else if (grey) { float gg = 0.1+0.9*float(k-2)/(nb-3); r=g=b=gg;} 
  else if (nbcolors<=1) {  
     float h=float(k-2)/(nb-2),s=1.,v=1.;
     hsvToRgb(h,s,v,r,g,b); 
     return;}     
  else   { //  interpolation dans la table hsv    
      int i= (k-2); 
      int j0= i*(nbcolors-1) / (nb-2);
      int j1=j0+1;
      int i0=  j0*(nb-2)/(nbcolors-1);
      int i1=  j1*(nb-2)/(nbcolors-1);
      int j03=j0*3,j13=j1*3;
      float a=float(i1-i)/(i1-i0),a1=1-a;
      if (hsv)
       {
        float h = colors[j03+0]*a + colors[j13+0]*a1;
        float s = colors[j03+1]*a + colors[j13+1]*a1;
        float v = colors[j03+2]*a + colors[j13+2]*a1;
        hsvToRgb(h,s,v,r,g,b); }
      else 
      {
       r = colors[j03+0]*a + colors[j13+0]*a1;
       g = colors[j03+1]*a + colors[j13+1]*a1;
       b = colors[j03+2]*a + colors[j13+2]*a1;
      }
    }     
 
}

