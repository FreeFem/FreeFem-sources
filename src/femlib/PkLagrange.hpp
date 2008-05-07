#include "FESpacen.hpp"

namespace Fem2D {

template<class Rd,class E>
static void SetPtPk(Rd *Pt,const int *dfon,int nn)
{  // P0 P1 et P2 
    const int d= E::Rd::d;
    int k=0;
    
    if(dfon[0])
	for(int i=0;i<=d;++i)
	    Pt[k++]=Rd();
    
    for(int i=0;i<d;++i)
	Pt[i+1][i]=1.;
    
    if(dfon[1])
	for(int i=0;i<E::ne;++i)
	    Pt[k++] = (Pt[E::nvedge[i][0]]+Pt[E::nvedge[i][1]])*0.5;
    
    if(dfon[d]==1) 
	Pt[k++]=Rd::diag(1./(d+1));
    
    assert(nn==k);   
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
      
      cout << "A4 "<<   k<< " "   <<dfon[0]<< dfon[1]<<dfon[2]<<dfon[3]<<endl;
    }
    operator const  int  * () const {return dfon;}
  };
  
  RdHat *Pt;
  TypeOfFE_Lagrange(int k):
    GTypeOfFE<Mesh>(A4(k),1,k,k<=2)
  {
    int n=this->NbDoF;
    
    cout << "\n +++ P"<<k<<" : ndof : "<< n <<endl;
    SetPtPk<RdHat,Element> (this->PtInterpolation,this->ndfOn(),this->NbDoF);
    cout << this->PtInterpolation<< endl;
    for (int i=0;i<n;i++) 
      {
	this->pInterpolation[i]=i;
	this->cInterpolation[i]=0;
	this->dofInterpolation[i]=i;
	this->coefInterpolation[i]=1.;
      }
  }
  
private:
  TypeOfFE_Lagrange( const TypeOfFE_Lagrange &) ;
  void operator=( const TypeOfFE_Lagrange &) ;
};


 }
