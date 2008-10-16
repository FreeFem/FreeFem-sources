// ********** DO NOT REMOVE THIS BANNER **********
// ORIG-DATE:     Jan 2008
// -*- Mode : c++ -*-
//
// SUMMARY  :Generic Pk Lagrange finite element class
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curi, Paris,  FRANCE 
// AUTHOR   : Frederic Hecht
// E-MAIL   : frederic.hecht@ann.jussieu.fr
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

 Thank to the ARN   FF2A3 grant
 ref:ANR-07-CIS7-002-01 
 */

#include "FESpacen.hpp"

namespace Fem2D {

template<class Rd,class E>
static void SetPtPk(Rd *Pt,const int *dfon,int nn)
{  // P0 P1 et P2 
    const int d= E::Rd::d;
    int k=0;
    
    if(dfon[0])
      {
	for(int i=0;i<=d;++i)
	  Pt[k++]=Rd();
	
	for(int i=0;i<d;++i)
	  Pt[i+1][i]=1.;
      }

    if(dfon[1]&& d !=1)
	for(int i=0;i<E::ne;++i)
	    Pt[k++] = (Pt[E::nvedge[i][0]]+Pt[E::nvedge[i][1]])*0.5;
    
    if(dfon[d]==1) 
	Pt[k++]=Rd::diag(1./(d+1));
    if(nn != k)
      { 
	cout << nn << " == " << k << " d = "<< d << " " << dfon[0]<< dfon[1]<<dfon[2]<<dfon[3]<<" "<< E::ne << endl;
	assert(nn==k); 
      }  
    if(verbosity>9)
      cout << " Pk = " << KN_<Rd>(Pt,nn)<<"\n";
    
}

//  a class of Lagrange Pk finite element 
template<class MMesh>
class TypeOfFE_Lagrange: public  GTypeOfFE<MMesh>
{
  //typedef typename  MMesh Mesh;
public:
  typedef   MMesh Mesh;
  typedef typename  Mesh::Element Element;
  typedef typename  Element::Rd Rd;
  typedef typename  Element::RdHat RdHat;
  static const int d=Rd::d;
  struct A4 {
    int dfon[4];
    
    A4(int k) {
      if(k==0) 
	{// P0
	  dfon[0]=dfon[1]=dfon[2]=dfon[3]=0;
	  dfon[d]=1;
	}
      else
	{
	  dfon[0]=1;
	  dfon[1]=max(k-1,0);
	  dfon[2]=d>1?max(k-2,0):0;
	  dfon[3]=d>2?max(k-3,0):0;}
    if(verbosity>9)      
      cout << "A4 "<<   k<< " "   <<dfon[0]<< dfon[1]<<dfon[2]<<dfon[3]<<endl;
    }
    operator const  int  * () const {return dfon;}
  };
  
  RdHat *Pt;
  TypeOfFE_Lagrange(int k):
    GTypeOfFE<Mesh>(A4(k),1,Max(k,1),k<=2,k==0)
  {
    int n=this->NbDoF;
    if(verbosity>9)    
    cout << "\n +++ P"<<k<<" : ndof : "<< n <<endl;
    SetPtPk<RdHat,Element> (this->PtInterpolation,this->ndfOn(),this->NbDoF);
    if(verbosity>9)    cout << this->PtInterpolation<< endl;
    for (int i=0;i<n;i++) 
      {
	this->pInterpolation[i]=i;
	this->cInterpolation[i]=0;
	this->dofInterpolation[i]=i;
	this->coefInterpolation[i]=1.;
      }
  }
  ~TypeOfFE_Lagrange(){ } //cout << "TypeOfFE_Lagrange"<< this->NbDoF<<endl;}
private:
  TypeOfFE_Lagrange( const TypeOfFE_Lagrange &) ;
  void operator=( const TypeOfFE_Lagrange &) ;
};


 }
