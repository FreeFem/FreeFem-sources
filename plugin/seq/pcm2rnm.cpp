/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : Tools to read ppm file
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-universite.fr

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep: pcm.cpp
/* clang-format on */

/*  use in freefem++ edp
 * see :
 * real[int,int] ff1("tt.pmm"); // read  image and set to an array.
 * real[int]  ff(ff1.nx*ff1.ny);
 * ff=ff1;
 */
// tools to read ppm file
/*  use in freefem++ edp file:
 * -----------------------------
 * complex[int,int] cc(1,1);
 * readpcm("tt.pcm",cc); // read  the flow image and set to un complex matrix array.
 * or
 * real[int,int] u(1,1),v(1,1);
 * readpcm("tt.pcm",u,v); // read the flow image and set to 2  real  matrix array.
 */
#include "pcm.hpp"
#include <iostream>
#include <cfloat>

using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
using namespace std;

#include "RNM.hpp"
#include <cmath>

long read1(const long &, const long &) { return 1; }

KNM< Complex > *read_pcm(string *filename, KNM< Complex > *p) {
  PCM pcm(filename->c_str( ));

  p->resize(pcm.width, pcm.height);
  pcm_complex *pc = pcm.image;

  for (int j = 0; j < pcm.height; ++j) {
    for (int i = 0; i < pcm.width; ++i, pc++) {
      (*p)(i, j) = Complex(pc->r, pc->i);
    }
  }

  return p;
}

long read_pcm(string *const &filename, KNM< double > *const &u, KNM< double > *const &v) {
  PCM pcm(filename->c_str( ));

  cout << " pcm  " << filename->c_str( ) << " : " << pcm.width << " x " << pcm.height << endl;
  u->resize(pcm.width, pcm.height);
  v->resize(pcm.width, pcm.height);
  pcm_complex *pc;
  float x1 = -1e+30, x2 = -1e+30;

  for (int j = 0; j < pcm.height; ++j) {
    for (int i = 0; i < pcm.width; ++i) {
      pc = pcm.Get(i, j);
      if (pc) {
        (*u)(i, j) = pc->r;
        (*v)(i, j) = pc->i;

        x1 = max(x1, pc->r);
        x2 = max(x2, pc->i);
        if (i < 0 && j < 0) {
          cout << i << " " << j << " " << pc->r << " " << pc->i << endl;
        }
      }
    }
  }

  cout << " max uv : " << x1 << " " << x2 << endl;
  return (long)pcm.width * pcm.height;
}

template<class T> inline T max (const T &a,const T & b,const T & c){return max(max(a,b),c);}
template<class T> inline T min (const T &a,const T & b,const T & c){return min(min(a,b),c);}

void rgb2hsv( double  r,double g, double b,double & h,double &s, double &v )
{

    double vmax = max(r, g, b);
    double vmin = min(r, g, b);
    double delta =vmax - vmin;
    v = vmax;
    

    if (vmax == 0.0) {
        s = 0;
        h = 0;
    }
    else if (delta < 1e-4) {
        s = 0;
        h = 0;
    }
    else {
        s = delta / vmax;

        if (vmax == r) {
            h =  ((g - b) / delta) + 0.;
        }
        else if (vmax == g) {
            h =  ((b - r) / delta) + 2.;
        }
        else {
            h =  ((r - g) / delta) + 4.;
        }
    }

    if (h < 0) h += 6.;
    h /= 6. ;

}
/*
void hsv2rgb(const unsigned char &src_h, const unsigned char &src_s, const unsigned char &src_v, unsigned char &dst_r, unsigned char &dst_g, unsigned char &dst_b)
{
    double h = src_h *   2.0f; // 0-360
    double s = src_s / 255.0f; // 0.0-1.0
    double v = src_v / 255.0f; // 0.0-1.0

    double r, g, b; // 0.0-1.0

    int   hi = (int)(h / 60.0f) % 6;
    double f  = (h / 60.0f) - hi;
    double p  = v * (1.0f - s);
    double q  = v * (1.0f - s * f);
    double t  = v * (1.0f - s * (1.0f - f));

    switch(hi) {
        case 0: r = v, g = t, b = p; break;
        case 1: r = q, g = v, b = p; break;
        case 2: r = p, g = v, b = t; break;
        case 3: r = p, g = q, b = v; break;
        case 4: r = t, g = p, b = v; break;
        case 5: r = v, g = p, b = q; break;
    }

    dst_r = (unsigned char)(r * 255); // dst_r : 0-255
    dst_g = (unsigned char)(g * 255); // dst_r : 0-255
    dst_b = (unsigned char)(b * 255); // dst_r : 0-255
}
 */
static void Load_Init( ) {
  cout << " load: init pcm2rmn  " << endl;

  Global.Add("readpcm", "(",
             new OneOperator2< KNM< Complex > *, string *, KNM< Complex > * >(&read_pcm),
             new OneOperator3_< long, string *, KNM< double > *, KNM< double > * >(&read_pcm));
}

LOADFUNC(Load_Init)
