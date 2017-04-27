// ORIG-DATE:     March 2015
// -*- Mode : c++ -*-
//
// SUMMARY  : liaison medit freefem++ : adaptmesh in 3d
// USAGE    : LGPL
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE
// AUTHOR   : Jacques Morice
// E-MAIL   : jacques.morice@ann.jussieu.fr
//
//   for automatic  compilation with ff-c++
//ff-c++-LIBRARY-dep:  mshmet libMesh
//ff-c++-cpp-dep:
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
 
 Thank to the ARN ()  FF2A3 grant
 ref:ANR-07-CIS7-002-01
 */
//ff-c++-LIBRARY-dep:   libMesh
//compile: ff-c++ -auto aniso.cpp
#include "ff++.hpp"
#include "msh3.hpp"
#include "eigenv.h"
//#define ADAPTLIBRARY
/*
 int eigenv(int symmat,double *mat,double lambda[3],double v[3][3]);
 int eigen2(double *mm,double *lambda,double vp[2][2]);

 */
static const int wrapperMetric[6]={0,1,3,2,4,5};

void met3dwrap(double m[6],double mw[6])
{
  for(int i=0;i<6;++i)
      mw[i]=m[wrapperMetric[i]];
}
void unmet3dwrap(double mw[6],double m[6])
{
    for(int i=0;i<6;++i)
        m[wrapperMetric[i]]=mw[i];
}


int BoundAniso2d(double m[3],double cmin)
{
    double vp[2][2];
    double l[2];
    int nv=eigen2(m,l,vp);
    if(nv)
    {
    double vpmx=max(l[0],l[1]);
    double vpmn = vpmx/cmin;
    l[0]=max(l[0],vpmn);
    l[1]=max(l[1],vpmn);
    m[0] = l[0]*vp[0][0]*vp[0][0] + l[1]*vp[1][0]*vp[1][0];
    m[1] = l[0]*vp[0][0]*vp[0][1] + l[1]*vp[1][0]*vp[1][1];
    m[2] = l[0]*vp[0][1]*vp[0][1] + l[1]*vp[1][1]*vp[1][1];
    }
    return nv ? 0: 1;
}
int BoundAniso3d(double mff[6],double cmin)
{
    double m[6];
    double vp[3][3];
    double l[3];
    met3dwrap(mff,m);
    int nv=eigenv(1,m,l,vp);
    if(nv)
    {
    double vpmx=max(max(l[0],l[1]),l[2]);
    double vpmn = vpmx/cmin;
    l[0]=max(l[0],vpmn);
    l[1]=max(l[1],vpmn);
    l[2]=max(l[2],vpmn);
   
    for (int k=0,i=0; i<3; i++)
        for (int j=i; j<3; j++)
            m[k++] = l[0]*vp[0][i]*vp[0][j] + l[1]*vp[1][i]*vp[1][j] \
            + l[2]*vp[2][i]*vp[2][j];
    }
    unmet3dwrap(m,mff);
    return nv ? 0: 1;
}
long Boundaniso(long const &k,KN<double> * const &pm,double const &animax)
{
    KN<double> & m(*pm);
    long ns = m.N()/k;
    ffassert( ns*k== m.N());
    double lmin = sqrt(animax);
    ffassert( k==3 || k==6);
    int err=0;
    if( k==3)
    {// <[m11,m12,m22]
        for(int i=0; i<ns;++i)
            err+=BoundAniso2d(&m[i*k],animax);
            
    }
    else
    {
     // 3d [m11,m21,m22,m31,m32,m33
        for(int i=0; i<ns;++i)
            err+=BoundAniso3d(&m[i*k],animax);

    }
    ffassert(err==0);
    return ns;
}

static void Load_Init(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
    
    //if (verbosity)
    if(verbosity) cout << " load: aniso  " << endl;
    
    Global.Add("boundaniso","(",new OneOperator3_<long,long,KN<double>*,double>(Boundaniso));
    
}


#define  WITH_NO_INIT
#include "msh3.hpp"
LOADFUNC(Load_Init)
