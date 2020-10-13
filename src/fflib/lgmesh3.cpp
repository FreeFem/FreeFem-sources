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

#include "ff++.hpp"
#include "array_resize.hpp"
#include "AFunction_ext.hpp"
#include "PlotStream.hpp"

using Fem2D::Mesh;
using Fem2D::MeshPoint;

extern bool NoWait; 

typedef Mesh const * pmesh;
typedef Mesh3 const * pmesh3;
typedef MeshS const * pmeshS;
typedef MeshL const * pmeshL;

map<pair<int,int>,int>::iterator closeto(map<pair<int,int>,int> & m, pair<int,int> & k)
{
    map<pair<int,int>,int>::iterator it= m.end();
    for (int i=-1;i<2;++i)
	for (int j=-1;j<2;++j)
	  {
	      pair<int,int>  kk(k.first+i,k.second+j);// bug corrigie april 2011 FH. tanks . J. Morice
	      it=  m.find(kk);
	      if(it != m.end()) return it;
	  }
    return it;
}

 inline void Perm3_I2J(const int *I,const int*J,int *S)
{
    if(I[0]==J[0]) S[0]=0;
    else if(I[0]==J[1]) S[0]=1;
    else {S[0]=2; assert(I[0]==J[2]) ;}
    if(I[1]==J[0]) S[1]=0;
    else if(I[1]==J[1]) S[1]=1;
    else {S[1]=2; assert(I[1]==J[2]) ; }
    S[2]=3-S[0]-S[1];
    assert(I[2]==J[3-S[0]-S[1]]);
}
template<class Mesh>
bool BuildPeriodic( 
		   int nbcperiodic ,
		   Expression *periodic,
		   const Mesh &Th,Stack stack,
		   KN<int> & ndfe) 
{ 
    typedef typename Mesh::BorderElement BE;
    /*
     build numbering of vertex form 0 to nbdfv-1
     and build numbering  of  edge form 0 to nbdfe-1
     we removing common vextex or common edge    
     --  we suppose one df by vertex 
     nbdfv number of df on vertex 
     ndfv[i]  given the numero of the df of the vertex 
     -- we suppose 1 df
     */ 
    //////////////////////////typedef Mesh3::BorderElement BE;
    //typedef Smallvect<int,2> int2;    
    if (nbcperiodic ) {
	
	//    KN<int> ndfv(Th.nv);
	//   KN<int> ndfe(Th.nbe);
	//ffassert(ndfv.N()==Th.nv);
	//ffassert(ndfe.N()==Th.nbe);
        
	MeshPoint *mp=MeshPointStack(stack),smp=*mp;   
	int n= nbcperiodic;
	if (verbosity >2)
	    cout << " Nb of pair of periodic conditions (3d) : = " << n <<  endl;
	int * link1=0;
	int * link2=0;
	KN<int*> plk1(n),plk2(n);
	KN<int> nlk1(n),nlk2(n);
	KN<int> lab1(n),lab2(n);
#ifndef  HUGE_VAL      
	const double infty= numeric_limits<double>::infinity();
#else
	const double infty= HUGE_VAL;
#endif       
	int nblink1, nblink2;
	int *plink1 , *plink2;
        for (int step=0;step<2;step++)
	  {
	      nblink1=0,     nblink2=0;
	      plink1=link1,  plink2=link2;
	      for (int ip=0, k=0;ip<n;ip++,k+=6)
		{
		    int label1=GetAny<long>((*periodic[k+0])(stack));
		    int label2=GetAny<long>((*periodic[k+3])(stack));
		    lab1[ip]=label1;
		    lab2[ip]=label2;
		    
		    int l1=nblink1;
		    int l2=nblink2;
		    plk1[ip]= plink1;
		    plk2[ip]= plink2;
		    for (int ke=0;ke<Th.nbe;ke++)
		      {
			  if (Th.be(ke).lab==label1)
			    {
				if (plink1) *plink1++=ke;
				nblink1++;
			    }
			  else if (Th.be(ke).lab==label2)
			    {
				if (plink2) *plink2++=ke;
				nblink2++;
			    }
		      }
		    nlk1[ip]= nblink1-l1;
		    nlk2[ip]= nblink2-l2;              
		}
	      if(step) break; // no reallocl 
	      if (verbosity >3)
		  cout << "  Periodic = " << nblink1 << " " << nblink2 << " step=" << step << endl;
	      if(nblink1 != nblink2)
		    ExecError("Periodic 3d:  the both number of face is not the same ");
		
	      
	      ndfe.resize(nblink1*2);
	      link1 = new int[nblink1];
	      link2 = new int[nblink2];
	  }
        if ( nblink1 >0) 
	  {
	      int indfe=0;
	      for (int ip=0, k=0;ip<n;ip++,k+=6)
		{
		    map<pair<int,int>,int> m;
		    const int kkx1=1,kky1=2,kkx2=4,kky2=5;
		    int label1=lab1[ip],label2=lab2[ip];
		    int n1=nlk1[ip],n2=nlk2[ip];
		    int *pke1=plk1[ip], *pke2=plk2[ip];
		    //int oip=pke1-link1;
		    typedef HashTable<SortArray<int,3>,int> HTable;
		    typedef HTable::iterator HTiterator;
		    HashTable<SortArray<int,3>,int> table1(n1,Th.nv); //  Table of face lab 1
		    R2 P[3];
		    R2 Pmin(infty,infty),Pmax(-infty,-infty);
		    double hmn=infty;
		    int iface[3];
		    if (verbosity >1)
			cout << "  --Update: periodic  couple label1= " << label1 
			     << ", n faces= " << n1 << "; "
			<< ", label2= " << label2<<  ", n faces= " << n2 <<endl; 
		    if (n1 != n2) 
		      ExecError("periodic 3D BC:  the number of set of faces is not the same ");
		    //  compute the hmn size to find common point 
		    for (int i1=0;i1<n1;i1++)
		      {
			  const BE & e =Th.be(pke1[i1]);
			  
			  assert(e.lab==label1) ;
			  {   
			      for(int ee=0;ee<3;++ee)
				{
				    iface[ee]=Th(e[ee]);
				    mp->set(e[ee].x,e[ee].y,e[ee].z);
				    P[ee].x=GetAny<double>((*periodic[k+kkx1])(stack));
				    P[ee].y=GetAny<double>((*periodic[k+kky1])(stack));
				}
			      
			      HTiterator hte=table1.add(SortArray<int,3>(iface),pke1[i1]);
			      ffassert(hte-table1.begin() == i1);
			      for(int ee=0,eo=2;ee<3;eo=ee++)
				{
				    double l = (P[ee]-P[eo]).norme2();
				    Pmin=Minc(Pmin,P[ee]);
				    Pmax=Maxc(Pmax,P[ee]);
				    hmn=Min(hmn,l);
				}
			      
			      
			  }                                
		      }
		    hmn=sqrt(hmn);
		    ffassert(hmn>1.0e-20);
		    double coef = 8/hmn;
		    double x0 = Pmin.x;
		    double y0 = Pmin.y;
		    if (verbosity > 2)
			cout << "  --Update: periodic " << Pmin << " " << Pmax << " " << " h=" << hmn 
			<< " ,  coef = "<< coef << " / " << (Pmax-Pmin).norme2()*coef*coef << endl;
		    ffassert(!n1 || (coef>1e-10 && (Pmax-Pmin).norme2()*coef*coef < 1.e14 )); // correct  FH mars 2013 
		    
		    //  map construction ----
		    for (int i1=0;i1<n1;i1++)
		      {
			  int ie=pke1[i1];
			  const BE & e =Th.be(ie);
			  assert (e.lab==label1);
			  for (int ne=0;ne<3;ne++)
			    {
				int kv=Th(e[ne]);
				
				// cout << ne << " " << kv << " " << " " << e[ne] << " ";
				mp->set(e[ne].x,e[ne].y,e[ne].z);
				double xx=GetAny<double>((*periodic[k+kkx1])(stack));
				double yy=GetAny<double>((*periodic[k+kky1])(stack));
				pair<int,int> ij((int) ((xx-x0)*coef),(int) ((yy-y0)*coef));				    
				map<pair<int,int>,int>::iterator im=closeto(m,ij);
				if (im==m.end())
				  {
				      if (verbosity >50)
					  cout << kv << " " << xx << " " << yy << " ->   " << ij << " :: " << ie << endl;
				      im=m.insert(pair<pair<int,int>,int>(ij,kv)).first;
				  }
				else {
				    if(im->second != kv)
					cout << kv << " ==  " << im->second << " " << xx << " " << yy << " ->   " << ij << " == " << im->first << endl;
				    ffassert( im->second == kv);
				}
				
			    }                                
		      }
		    //  find ...  face of list 2 in list 1 .... 
		    int err=0;
		    for (int i2=0;i2<n2;i2++)
		      {
			  int ie2=pke2[i2];
			  const BE & e =Th.be(ie2);
			  assert (e.lab==label2);
			{
			    if (verbosity >50)
				cout << ie2 << " : " <<Th(e[0]) << " " << Th(e[1]) << " " << Th(e[2]) << ":: ";
			    R2 P[3];
			    pair<int,int> I[3];
			    map<pair<int,int>,int>::iterator im;
			    int i2to1[3];
			    
			    for(int ee=0;ee<3;++ee)
			      {
				  mp->set(e[ee].x,e[ee].y,e[ee].z);
				  
				  P[ee].x=GetAny<double>((*periodic[k+kkx2])(stack));
				  P[ee].y=GetAny<double>((*periodic[k+kky2])(stack));
				  I[ee].first = int((P[ee].x-x0)*coef);
				  I[ee].second= int((P[ee].y-y0)*coef);
				  im=closeto(m,I[ee]);
				  
				  if(im == m.end() )
				    {
					cout << " vertex : " << Th(e[ee]) << " " <<e[ee]<< "  Abscisses: s = "<< P[ee]   << "  " <<  I[ee] << endl; 
					ExecError("periodic: Sorry one vertex of face is losted "); 
				    }
				  i2to1[ee] = im->second;
				  
			      }
			    
			    //int ie1=-1;
			    SortArray<int,3> sf(i2to1);
			    HTiterator ht=table1.find(sf);
			    if( ! ht ) {
				err++;
				cerr << " missing face " << ie2 << " " << sf <<  endl;
			    }
			    else
			      {
				
				  int ie1 = ht->v; 
				  assert(ie1>=0  );
				  const BE & eo =Th.be(ie1);
				  int fo[3]={Th(eo[0]),Th(eo[1]),Th(eo[2])};
				  
				  int np1= NumPerm1<3>(fo); //  number of 
				  int np2= NumPerm1<3>(i2to1);
				  ndfe[indfe++]=ie1*8+np1; 
				  ndfe[indfe++]=ie2*8+np2; 
				  int p1[3],p2[3];
				  SetNumPerm<3>(np1,p1);
				  SetNumPerm<3>(np2,p2);
				  if(verbosity>50)
				      cout <<"  " << ie1 << " ==  " << ie2  << ":  " <<  fo[p1[0]] << " " << fo[p1[1]] << " " << fo[p1[2]] << " == " 
				       << i2to1[p2[0]] << " " << i2to1[p2[1]] << " " << i2to1[p2[2]]  << "  e ="
				       << Th(e[p2[0]]) << " " << Th(e[p2[1]]) << " " << Th(e[p2[2]])  << " "
 				      
				           << " nu= " << np1 << " , " << np2  << endl;
				 
			      }
			    
			}
		      }
		    if(err ) {
			cerr << " 3d periodic FE: number of  missing periodic faces " << err << endl; 
			ExecError(" 3d periodic condition missing common face ");
		    }
		    
		    
		}
	      *mp = smp;
	      ffassert(indfe==ndfe.size());

	      /*  rm  FH :
	      for (int i=0;i<Th.nbe;i++)
		  ndfe[i]=i;// circular link
	      for (int i=0;i<Th.nv;i++)
		  ndfv[i]=i;// circular link
	      for (int i=0;i<nblink1;i++)
		{
		    int ie1=-link1[i]/8;
		    int ie2=link2[i];
		    int np=-link1[i]%8;
		    int p[3]={np%3,np/3,3-np%3-np/3};
		    if (verbosity >50)
			cout << " face " << ie1 << " <==> " << ie2 << endl;
		    ffassert(ie1!=ie2);
		    if(!InCircularList(ndfe,ie1,ie2))   // merge of two list 
			Exchange(ndfe[ie1],ndfe[ie2]);
		    int kv1[3],kv2[3];
		    // kv1 == kv2 (p) 
		    const BE & e1=Th.be(ie1), &e2=Th.be(ie2);
		    for (int ee=0;ee<3;ee++)
		      {
			  kv1[ee] =Th(e1[ee]);
			  kv2[ee] =Th(e2[ee]);		           
		      }
		    for (int ee=0;ee<3;ee++)
		      {
			  int iv1=kv1[ee];
			  int iv2=kv2[p[ee]];
			  //  l'orientation de la face 2 / face 1  est  dans p.
			  
			  if (!InCircularList(ndfv,iv1,iv2)) {  // merge of two list 
			      Exchange(ndfv[iv2],ndfv[iv1]);
			      
			      if (verbosity >50)
				{ 
				    int ii=iv1,l=1;
				    while ( (ii=ndfv[ii]) != iv1 && l++<10) (void) 0;
				    if( (verbosity >50) || ( l > 2 && verbosity >40)) 
				     cout << l << "  vertex " << iv1 <<  "<==> " << iv2 << " list : " << iv1 << " " << Th(iv1) << " <=> " << Th(iv2);
				    int i=iv1,k=0;
				    while ( (i=ndfv[i]) != iv1 && k++<10)
					cout << ", "<< i ; 
				    cout << endl;
				    
				}}                  
		      }
		    
		} 
	      // generation de numero de dlt
	      
	      nbdfv = numeroteclink(ndfv) ; 
	      nbdfe = numeroteclink(ndfe) ; 
	      if (verbosity>2) 
		  cout << "  -- nb df on vertices " << nbdfv << endl;
	      */
	      delete [] link1;
	      delete [] link2;
	      return true; //new FESpace(**ppTh,*tef,nbdfv,ndfv,nbdfe,ndfe);
	  }
        else {
	      delete [] link1;
	      delete [] link2;
        }
	
    }
    return false;   
}


bool  v_fes3::buildperiodic(Stack stack, KN<int> & ndfe) { 
    return BuildPeriodic(nbcperiodic,periodic,**ppTh,stack,ndfe);
    
}

bool  v_fesS::buildperiodic(Stack stack, KN<int> & ndfe) {
    return BuildPeriodic(nbcperiodic,periodic,**ppTh,stack,ndfe);
    
}

bool  v_fesL::buildperiodic(Stack stack, KN<int> & ndfe) {
    return BuildPeriodic(nbcperiodic,periodic,**ppTh,stack,ndfe);
    
}

template<class Mesh> 
class GlgVertex {
public:
  typedef double R;
  typedef typename Mesh::Rd Rd;
  typedef typename Mesh::Vertex Vertex;
  CountPointer<const Mesh> pTh;
  const Vertex *v;
  void Check() const {   if (!v || !pTh) { ExecError("Too bad! Unset Vertex!"); } }
  void init() { v=0;pTh.init();}
  GlgVertex(const Mesh * Th,long kk): pTh(Th),v( &(*pTh)(kk)) {}
  GlgVertex(const Mesh * Th,const Vertex * kk): pTh(Th),v(kk) {}
  operator int() const { Check(); return (* pTh)(v);} 
  operator Rd*(){ Check(); return v;} 
  R x() const {Check() ; return v->X();}
  R y() const {Check() ; return v->Y();}
  R z() const {Check() ; return v->Z();}
  long lab() const {Check() ; return v->lab;}
  void destroy()  {pTh.destroy();}
};

template<class Mesh>
class GlgElement { public:
typedef typename Mesh::Element Element;

    struct Adj {// 
	const Mesh *pTh;
	const Element *k;
	Adj(const GlgElement<Mesh> & pp) : pTh(pp.pTh),k(pp.k) {}
	GlgElement<Mesh> adj(long & e) const  {
	    int ee,ko;
	    ffassert(pTh && k && e >=0 && e < Element::nv );
	    long kk=pTh->ElementAdj(ko=(*pTh)(k),ee=e);
	    if(ee>=0) e=ee;//  ok adj exist
	    else  return  GlgElement<Mesh>(pTh,ko);// return same .. 
	    
	return  GlgElement<Mesh>(pTh,kk);}
    }; 
    
    CountPointer<const Mesh> pTh;
  const Element *k;
  
  GlgElement():  k(0) {}
  void  Check() const  {   if (!k || !pTh) { ExecError("Unset Triangle,Sorry!"); } }
  void init() { k=0;pTh.init();}
  void destroy() {pTh.destroy();}
  GlgElement(const Mesh * Th,long kk): pTh(Th),k( &(*pTh)[kk]) {}
  GlgElement(const Mesh * Th,Element * kk): pTh(Th),k(kk) {}
  operator int() const { Check(); return (* pTh)(k);} 
  GlgVertex<Mesh> operator [](const long & i) const { Check(); return GlgVertex<Mesh>(pTh,&(*k)[i]);}   
  long lab() const {Check() ; return k ? k->lab : LONG_MIN;}
  double mes() const {Check() ; return k->mesure() ;}
  long n() const { return k ? Element::nv: 0 ;}
    
  bool operator==(const GlgElement & l) const { return pTh==l.pTh && k == l.k;}
  bool operator!=(const GlgElement & l) const { return pTh!=l.pTh || k != l.k;}
  bool operator<(const GlgElement & l) const { return pTh==l.pTh && k <l.k;}
  bool operator<=(const GlgElement & l) const { return pTh==l.pTh && k <=l.k;}
    
};

template<class Mesh>
class GlgBoundaryElement { public:
    
    typedef typename Mesh::Element Element;
    typedef typename Mesh::BorderElement BorderElement;
    
    struct BE {
	const Mesh * p;
	BE(const Mesh *pp) : p(pp) {}
	BE(const Mesh **pp) : p(*pp) {}
	operator const Mesh * () const {return p;}
    };
    
    CountPointer<const Mesh> pTh;
    const BorderElement *k;
    
    GlgBoundaryElement():  k(0) {}
    void  Check() const  {   if (!k || !pTh) { ExecError("Unset BoundaryEdge,Sorry!"); } }
    void init() { k=0;pTh.init();}
    void destroy() {pTh.destroy();}
    GlgBoundaryElement(const Mesh * Th,long kk): pTh(Th),k( &(*pTh).be(kk)) {}
    GlgBoundaryElement(const Mesh * Th,BoundaryEdge * kk): pTh(Th),k(kk) {}
    GlgBoundaryElement(const BE & be,long kk): pTh(be.p),k( &(*pTh).be(kk)) {}
    GlgBoundaryElement(const BE & be,BoundaryEdge * kk): pTh(be.p),k(kk) {}
    operator int() const { Check(); return (* pTh)(k);} 
    GlgVertex<Mesh> operator [](const long & i) const { Check(); return GlgVertex<Mesh>(pTh,&(*k)[i]);}   
    long lab() const {Check() ; return k ? k->lab : 0;}
    double length() const {Check() ; return k->length()  ;}
    long n() const { return k ? BorderElement::nv : 0 ;}
    GlgElement<Mesh> element() const {Check() ;int ee; return GlgElement<Mesh>(pTh,(*pTh).BoundaryElement((*pTh)(k),ee));}
    long nuBoundaryElement() const {Check() ;int ee;  (*pTh).BoundaryElement((*pTh)(k),ee);return ee;}
    
};


// Tools for 3D volume mesh
GlgBoundaryElement<Mesh3> get_element(GlgBoundaryElement<Mesh3>::BE const & a, long const & n){  return GlgBoundaryElement<Mesh3>(a,n);}
GlgVertex<Mesh3> get_element(GlgBoundaryElement<Mesh3> const & a, long const & n){  return a[n];}

GlgElement<Mesh3> get_adj(GlgElement<Mesh3>::Adj const & a, long  * const & n){return  a.adj(*n);}

GlgElement<Mesh3> get_element(pmesh3 const & a, long const & n) {  return GlgElement<Mesh3>(a,n);}
GlgElement<Mesh3> get_element(pmesh3 *const & a, long const & n) {  return GlgElement<Mesh3>(*a,n);}

GlgVertex<Mesh3> get_vertex(pmesh3 const & a, long const & n){ return GlgVertex<Mesh3>(a,n);}
GlgVertex<Mesh3> get_vertex(pmesh3 *const & a, long const & n){ return GlgVertex<Mesh3>(*a,n);}
GlgVertex<Mesh3> get_element(GlgElement<Mesh3> const & a, long const & n) {  return a[n];}

GlgElement<Mesh3> getElement(GlgBoundaryElement<Mesh3> const & a)
{    return a.element();}

long NuElement(GlgBoundaryElement<Mesh3> const & a)
{    return a.nuBoundaryElement(); }

R getx(GlgVertex<Mesh3> const & a){  return a.x();}
R gety(GlgVertex<Mesh3> const & a){  return a.y();}
R getz(GlgVertex<Mesh3> const & a){  return a.z();}
long  getlab(GlgVertex<Mesh3> const & a){  return a.lab();}
long getlab(GlgElement<Mesh3> const & a){  return a.lab();}
long getlab(GlgBoundaryElement<Mesh3> const & a){  return a.lab();}
R getmes(GlgElement<Mesh3> const & a){  return a.mes();}

double pmesh_mes(pmesh3 * p) { ffassert(p) ;  return *p ? (**p).mes : 0.0;}
double pmesh_mesb(pmesh3 * p) { ffassert(p) ;  return *p ? (**p).mesb : 0.0;}
long pmesh_nt(pmesh3 * p) { ffassert(p) ;  return *p ? (**p).nt : 0;}
long pmesh_nv(pmesh3 * p) { ffassert(p) ;  return *p ? (**p).nv : 0;}
long pmesh_nbe(pmesh3 * p) { ffassert(p) ;  return *p ? (**p).nbe : 0;}

double pmesh_hmax(pmesh3 * p)
{ if(p && *p) {
    double hmax2 =0;
    const Mesh3 & Th = **p;
    for(int k=0; k< Th.nt; ++k)
        for(int e=0; e<6; ++e)
            hmax2=max(hmax2,Th[k].Edge(e).norme2());
    return sqrt(hmax2); }
  else return 0.0; }

double pmesh_hmin(pmesh3 * p)
{ if(p && *p) {
    double hmin2 =1e100;
    const Mesh3 & Th = **p;
    for(int k=0; k< Th.nt; ++k)
        for(int e=0; e<6; ++e)
            hmin2=min(hmin2,Th[k].Edge(e).norme2());
    return sqrt(hmin2); }
  else return 0.0; }


pmeshS pmesh3_gamma(Stack stack, pmesh3 * const & p)
{ throwassert(p && *p) ;
  const Mesh3 & Th = **p;
  const MeshS *ThS = Th.meshS;
  if(ThS==NULL) cout << "The meshS member is empty!, use before <mesh3>=buildBdMesh(<mesh3>)" << endl;
  return (ThS);
}

pmeshL pmeshS_gamma(Stack stack, pmeshS * const & p)
{ throwassert(p && *p) ;
    const MeshS & Th = **p;
    const MeshL *ThL = Th.meshL;
    if(ThL==NULL) cout << "The meshL member is empty!, use before <meshL>=buildBdMesh(<meshL>)" << endl;
    return (ThL);
}


long pmesh_nadjnomanifold(pmesh3 * p) { ffassert(p) ;  return *p ? ((**p).meshS)->nadjnomanifold : 0;}

long pmesh_nadjnomanifold(pmeshS * p) { ffassert(p) ;  return *p ? (**p).nadjnomanifold : 0;}


// Tools for 3D surface mesh
GlgBoundaryElement<MeshS> get_element(GlgBoundaryElement<MeshS>::BE const & a, long const & n){  return GlgBoundaryElement<MeshS>(a,n);}
GlgVertex<MeshS> get_element(GlgBoundaryElement<MeshS> const & a, long const & n){  return a[n];}

GlgElement<MeshS> get_adj(GlgElement<MeshS>::Adj const & a, long  * const & n){return  a.adj(*n);}

GlgElement<MeshS> get_element(pmeshS const & a, long const & n) {  return GlgElement<MeshS>(a,n);}
GlgElement<MeshS> get_element(pmeshS *const & a, long const & n) {  return GlgElement<MeshS>(*a,n);}

GlgVertex<MeshS> get_vertex(pmeshS const & a, long const & n){ return GlgVertex<MeshS>(a,n);}
GlgVertex<MeshS> get_vertex(pmeshS *const & a, long const & n){ return GlgVertex<MeshS>(*a,n);}
GlgVertex<MeshS> get_element(GlgElement<MeshS> const & a, long const & n) {  return a[n];}

GlgElement<MeshS> getElement(GlgBoundaryElement<MeshS> const & a)
{    return a.element();}

long NuElement(GlgBoundaryElement<MeshS> const & a)
{    return a.nuBoundaryElement(); }

R getx(GlgVertex<MeshS> const & a){  return a.x();}
R gety(GlgVertex<MeshS> const & a){  return a.y();}
R getz(GlgVertex<MeshS> const & a){  return a.z();}
long  getlab(GlgVertex<MeshS> const & a){  return a.lab();}
long getlab(GlgElement<MeshS> const & a){  return a.lab();}
long getlab(GlgBoundaryElement<MeshS> const & a){  return a.lab();}
R getmes(GlgElement<MeshS> const & a){  return a.mes();}

double pmesh_mes(pmeshS * p) { ffassert(p) ;  return *p ? (**p).mes : 0.0;}
double pmesh_mesb(pmeshS * p) { ffassert(p) ;  return *p ? (**p).mesb : 0.0;}
long pmesh_nt(pmeshS * p) { ffassert(p) ;  return *p ? (**p).nt : 0;}
long pmesh_nv(pmeshS * p) { ffassert(p) ;  return *p ? (**p).nv : 0;}
long pmesh_nbe(pmeshS * p) { ffassert(p) ;  return *p ? (**p).nbe : 0;}

double pmesh_hmax(pmeshS * p)
{ ffassert(p && *p) ;
    double hmax2 =0;
    const MeshS & Th = **p;
    for(int k=0; k< Th.nt; ++k)
        for(int e=0; e<3; ++e)
            hmax2=max(hmax2,Th[k].Edge(e).norme2());
    return sqrt(hmax2);}

double pmesh_hmin(pmeshS * p)
{ throwassert(p && *p) ;
    double hmin2 =1e100;
    const MeshS & Th = **p;
    for(int k=0; k< Th.nt; ++k)
        for(int e=0; e<3; ++e)
            hmin2=min(hmin2,Th[k].Edge(e).norme2());
    return sqrt(hmin2);}

template<class MMesh>
class Op3_MeshDmp : public quad_function< const MMesh* *, R, R, R, MeshPoint * > {
 public:
  class Op : public E_F0mps {
   public:
    typedef const MMesh *ppmesh;
    typedef typename MMesh::RdHat RdHat;
    typedef typename MMesh::Element T;
    Expression a, b, c, d;
    Op(Expression aa, Expression bb, Expression cc, Expression dd) : a(aa), b(bb), c(cc), d(dd) {}
    AnyType operator( )(Stack s) const {
      R xx(GetAny< R >((*b)(s)));
      R yy(GetAny< R >((*c)(s)));
      R zz(GetAny< R >((*d)(s)));
      ppmesh *ppTh(GetAny< ppmesh * >((*a)(s)));
      if (!ppTh || !*ppTh) ExecError("Op3_Mesh3mp unset mesh ??");
      ppmesh pTh(*ppTh);
      MeshPoint *mp = new MeshPoint( );
      mp->set(xx, yy, zz);
      RdHat PHat;
      bool outside;
      const T *K = pTh->Find(mp->P, PHat, outside);
      mp->set(*pTh, mp->P, PHat, *K, 0, outside);
      return mp;
    }
  };
};

GlgBoundaryElement<MeshL> get_element(GlgBoundaryElement<MeshL>::BE const & a, long const & n){  return GlgBoundaryElement<MeshL>(a,n);}
GlgVertex<MeshL> get_element(GlgBoundaryElement<MeshL> const & a, long const & n){  return a[n];}

GlgElement<MeshL> get_adj(GlgElement<MeshL>::Adj const & a, long  * const & n){return  a.adj(*n);}

GlgElement<MeshL> get_element(pmeshL const & a, long const & n) {  return GlgElement<MeshL>(a,n);}
GlgElement<MeshL> get_element(pmeshL *const & a, long const & n) {  return GlgElement<MeshL>(*a,n);}

GlgVertex<MeshL> get_vertex(pmeshL const & a, long const & n){ return GlgVertex<MeshL>(a,n);}
GlgVertex<MeshL> get_vertex(pmeshL *const & a, long const & n){ return GlgVertex<MeshL>(*a,n);}
GlgVertex<MeshL> get_element(GlgElement<MeshL> const & a, long const & n) {  return a[n];}

GlgElement<MeshL> getElement(GlgBoundaryElement<MeshL> const & a)
{    return a.element();}

long NuElement(GlgBoundaryElement<MeshL> const & a)
{    return a.nuBoundaryElement(); }

R getx(GlgVertex<MeshL> const & a){  return a.x();}
R gety(GlgVertex<MeshL> const & a){  return a.y();}
R getz(GlgVertex<MeshL> const & a){  return a.z();}
long  getlab(GlgVertex<MeshL> const & a){  return a.lab();}
long getlab(GlgElement<MeshL> const & a){  return a.lab();}
long getlab(GlgBoundaryElement<MeshL> const & a){  return a.lab();}
R getmes(GlgElement<MeshL> const & a){  return a.mes();}

double pmesh_mes(pmeshL * p) { ffassert(p) ;  return *p ? (**p).mes : 0.0;}
double pmesh_mesb(pmeshL * p) { ffassert(p) ;  return *p ? (**p).mesb : 0.0;}
long pmesh_nt(pmeshL * p) { ffassert(p) ;  return *p ? (**p).nt : 0;}
long pmesh_nv(pmeshL * p) { ffassert(p) ;  return *p ? (**p).nv : 0;}
long pmesh_nbe(pmeshL * p) { ffassert(p) ;  return *p ? (**p).nbe : 0;}

double pmesh_hmax(pmeshL * p)
{ ffassert(p && *p) ;
    double hmax2 =0;
    const MeshL & Th = **p;
    for(int k=0; k< Th.nt; ++k)
            hmax2=max(hmax2,Th[k].Edge(0).norme2());
    return sqrt(hmax2);}

double pmesh_hmin(pmeshL * p)
{ throwassert(p && *p) ;
    double hmin2 =1e100;
    const MeshL & Th = **p;
    for(int k=0; k< Th.nt; ++k)
            hmin2=min(hmin2,Th[k].Edge(0).norme2());
    return sqrt(hmin2);}










// 3D volume
pf3rbase* get_element(pf3rbasearray *const & a, long const & n)
{
    return (**a)[n];
}

pf3r get_element(pf3rarray const & a, long const & n)
{  //cout << " ************ " << n << " " << a.second << endl;
    return pf3r( *(*a.first)[n],a.second);
}

//  complex case 
pf3cbase* get_element(pf3cbasearray *const & a, long const & n)
{
    return (**a)[n];
}
pf3c get_element(pf3carray const & a, long const & n)
{
    return pf3c( *(*a.first)[n],a.second);
}
//  end complex case

// 3D surface

pfSrbase* get_element(pfSrbasearray *const & a, long const & n)
{
    return (**a)[n];
}

pfSr get_element(pfSrarray const & a, long const & n)
{  //cout << " ************ " << n << " " << a.second << endl;
    return pfSr( *(*a.first)[n],a.second);
}

//  complex case
pfScbase* get_element(pfScbasearray *const & a, long const & n)
{
    return (**a)[n];
}
pfSc get_element(pfScarray const & a, long const & n)
{
    return pfSc( *(*a.first)[n],a.second);
}

// 3D line

pfLrbase* get_element(pfLrbasearray *const & a, long const & n)
{
    return (**a)[n];
}

pfLr get_element(pfLrarray const & a, long const & n)
{  //cout << " ************ " << n << " " << a.second << endl;
    return pfLr( *(*a.first)[n],a.second);
}

//  complex case
pfLcbase* get_element(pfLcbasearray *const & a, long const & n)
{
    return (**a)[n];
}
pfLc get_element(pfLcarray const & a, long const & n)
{
    return pfLc( *(*a.first)[n],a.second);
}

class MoveMesh3 :  public E_F0mps { public:  
 
   typedef pmesh  Result;
   Expression getmesh;
   Expression U,V;
   int nbsol;    
    vector<Expression> sol;
   
    MoveMesh3(const basicAC_F0 & args) :nbsol(args.size()-2),sol(args.size()-2)
    {   
  using  Fem2D::Triangle;
  using  Fem2D::Vertex;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
 // using  Fem2D::R;
  using  Fem2D::MeshPointStack;
      args.SetNameParam();
      getmesh=to<pmesh>(args[0]); 
      const E_Array * a = dynamic_cast<const E_Array *>(args[1].LeftValue());
      
      ffassert(a);
      if (a->size() !=2) CompileError("movemesh(Th,[u,v],...) need 2 componate in array ",atype<pmesh>());
      U=to<double>( (*a)[0]);
      V=to<double>( (*a)[1]);
      
      for (int i=2;i<args.size();i++)
        sol[i-2]=to<double>(args[i]);      
    }   
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmesh>(),atype<E_Array>(),true);}
    static  E_F0 * f(const basicAC_F0 & args){ return new MoveMesh3(args);} 
    AnyType operator()(Stack s) const ;
  operator aType () const { return atype<Result>();} 

};

// 3D volume readMesh
class ReadMesh3 :  public E_F0 { public:  
    
  Expression filename; 
  typedef pmesh3  Result;
  ReadMesh3(const basicAC_F0 & args) 
  {   
    args.SetNameParam(); 
    filename=to<string*>(args[0]);   
  }   
  static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<string*>());}
  static  E_F0 * f(const basicAC_F0 & args){ return new ReadMesh3(args);} 
  AnyType operator()(Stack stack) const;
};


AnyType ReadMesh3::operator()(Stack stack) const 
{
  using  Fem2D::MeshPointStack;
 
  string * fn =  GetAny<string*>((*filename)(stack));
  if(verbosity > 2)
      cout << "ReadMesh3 " << *fn << endl;
  Mesh3 *Thh = new Mesh3(*fn);

  Thh->BuildGTree();
    if (Thh->meshS) Thh->meshS->BuildGTree();
  Add2StackOfPtr2FreeRC(stack,Thh);
  return SetAny<pmesh3>(Thh);;
  
}

// 3D surface readmesh

class ReadMeshS :  public E_F0 { public:
    
    Expression filename;
    typedef pmeshS  Result;
    ReadMeshS(const basicAC_F0 & args)
    {
        args.SetNameParam();
        filename=to<string*>(args[0]);
    }
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<string*>());}
    static  E_F0 * f(const basicAC_F0 & args){ return new ReadMeshS(args);}
    AnyType operator()(Stack stack) const;
};


AnyType ReadMeshS::operator()(Stack stack) const
{
    using  Fem2D::MeshPointStack;
    
    string * fn =  GetAny<string*>((*filename)(stack));
    if(verbosity > 2)
        cout << "ReadMeshS " << *fn << endl;
    MeshS *Th = new MeshS(*fn);
    Th->BuildGTree();
    Add2StackOfPtr2FreeRC(stack,Th);
    return SetAny<pmeshS>(Th);
    
}



class ReadMeshTEST_Op:  public E_F0mps
{
    public:
    
        Expression filename;
    //typedef pmeshS  Result;
        static const int n_name_param = 3;
        static basicAC_F0::name_and_type name_param [];
        Expression nargs[n_name_param];
        bool arg (int i, Stack stack, bool a) const {return nargs[i] ? GetAny<bool>((*nargs[i])(stack)) : a;}
    
    public:
        ReadMeshTEST_Op (const basicAC_F0 &args,Expression ffname)
        : filename(ffname) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }
        AnyType operator()(Stack stack) const;
    };


    basicAC_F0::name_and_type ReadMeshTEST_Op::name_param [] = {
        {"cleanmesh", &typeid(bool)},
        {"removeduplicate", &typeid(bool)},
        {"rebuildboundary", &typeid(bool)}
    };



class ReadMeshTEST: public OneOperator {
public:
    ReadMeshTEST (): OneOperator(atype<pmeshS>(), atype<string *>()) {}
    
    E_F0*code (const basicAC_F0 &args) const {
        return new ReadMeshTEST_Op(args,t[0]->CastTo(args[0]));
    }
};

AnyType ReadMeshTEST_Op::operator () (Stack stack)  const {

    string * fn =  GetAny<string*>((*filename)(stack));
    bool cleanmesh(arg(0, stack, false));
    bool removeduplicate(arg(1, stack, false));
    bool rebuildboundary(arg(2, stack, false));
    
    MeshS *Th = new MeshS(*fn, cleanmesh, removeduplicate, rebuildboundary);//, int orientation=1, double precis_mesh=1e-7););
    Th->BuildGTree();
    Add2StackOfPtr2FreeRC(stack,Th);
    return Th;
}



/*
AnyType ReadMeshTEST_Op::operator () (Stack stack)  const {
    using  Fem2D::MeshPointStack;
 
    bool cleanmesh(arg(0, stack, false));
    bool removeduplicate(arg(1, stack, false));
    bool rebuildboundary(arg(2, stack, false));
    //string * fn =  GetAny<string*>((*filename)(stack));
    if(verbosity > 2)
        cout << "ReadMeshS " << *fn << endl;
    MeshS *Th = new MeshS(*fn,cleanmesh,removeduplicate,rebuildboundary);
    Th->BuildGTree();
    Add2StackOfPtr2FreeRC(stack,Th);
    return SetAny<pmeshS>(Th);

}
 
*/


class ReadMeshL :  public E_F0 { public:
    
    Expression filename;
    typedef pmeshL  Result;
    ReadMeshL(const basicAC_F0 & args)
    {
        args.SetNameParam();
        filename=to<string*>(args[0]);
    }
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<string*>());}
    static  E_F0 * f(const basicAC_F0 & args){ return new ReadMeshL(args);}
    AnyType operator()(Stack stack) const;
};


AnyType ReadMeshL::operator()(Stack stack) const
{
    using  Fem2D::MeshPointStack;
    
    string * fn =  GetAny<string*>((*filename)(stack));
    if(verbosity > 2)
        cout << "ReadMeshL " << *fn << endl;
    MeshL *Th = new MeshL(*fn,0);  // param 0-> initialize just the meshS
    Th->BuildGTree();
    Add2StackOfPtr2FreeRC(stack,Th);
    return SetAny<pmeshL>(Th);
    
}









// 3D Volume SaveMesh
class SaveMesh3 :  public E_F0 { public:  
 
   typedef pmesh3  Result;
   Expression getmesh;
   Expression filename; 
   Expression xx,yy,zz;  
   SaveMesh3(const basicAC_F0 & args) 
    {   
      xx=0;
      yy=0;
      zz=0;
      args.SetNameParam();
      getmesh=to<pmesh3>(args[0]); 
      filename=to<string*>(args[1]);
      // what's it mean in 3d?
      /*if (args.size() >2)
        {
          const E_Array * a = dynamic_cast<const E_Array *>(args[2].LeftValue());
          if (!a) CompileError("savemesh(Th,\"filename\",[u,v,w],...");
          int k=a->size() ;
         // cout << k << endl;
          if ( k!=2 && k !=3) CompileError("savemesh(Th,\"filename\",[u,v,w]) need 2 or 3  componate in array ",atype<pmesh>());
          xx=to<double>( (*a)[0]);
          yy=to<double>( (*a)[1]);
          if(k==3)
           zz=to<double>( (*a)[2]);
         }*/
      
   }   
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmesh3>(),atype<string*>(),true);}
    static  E_F0 * f(const basicAC_F0 & args){ return new SaveMesh3(args);} 
    AnyType operator()(Stack s) const ;
  
};


AnyType SaveMesh3::operator()(Stack stack) const 
{
  using  Fem2D::MeshPointStack;
  
  
   pmesh3 Thh = GetAny<pmesh3>((*getmesh)(stack));
   string * fn =  GetAny<string*>((*filename)(stack));
   
   if (verbosity > 2)
       cout << "SaveMesh3 " << *fn << " " << Thh << endl;
   
    int ret=Thh->Save(*fn);
    
    if( ret!=0) {ExecError("PB Write error !");}
   return SetAny<pmesh3>(Thh);

}

// 3D surface SaveMesh
class SaveMeshS :  public E_F0 { public:
    
    typedef pmeshS Result;
    Expression getmesh;
    Expression filename;
    //Expression xx,yy,zz;
    SaveMeshS(const basicAC_F0 & args)
    {
        //xx=0;
        //yy=0;
        //zz=0;
        args.SetNameParam();
        getmesh=to<pmeshS>(args[0]);
        filename=to<string*>(args[1]);
        // what's it mean in 3d?
        /* if (args.size() >2)
        {
            const E_Array * a = dynamic_cast<const E_Array *>(args[2].LeftValue());
            if (!a) CompileError("savemesh(Th,\"filename\",[u,v,w],...");
            int k=a->size() ;
            // cout << k << endl;
            if ( k!=2 && k !=3) CompileError("savemesh(Th,\"filename\",[u,v,w]) need 2 or 3  componate in array ",atype<pmeshS>());
            xx=to<double>( (*a)[0]);
            yy=to<double>( (*a)[1]);
            if(k==3)
                zz=to<double>( (*a)[2]);
        }*/
        
    }
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmeshS>(),atype<string*>(),true);}
    static  E_F0 * f(const basicAC_F0 & args){ return new SaveMeshS(args);}
    AnyType operator()(Stack s) const ;
    
};

AnyType SaveMeshS::operator()(Stack stack) const
{
    using  Fem2D::MeshPointStack;
    
    
    pmeshS Thh = GetAny<pmeshS>((*getmesh)(stack));
    string * fn =  GetAny<string*>((*filename)(stack));
    
    if (verbosity > 2)
        cout << "SaveMeshS " << *fn << " " << Thh << endl;
    int ret=Thh->Save(*fn);
    if( ret!=0) {ExecError("PB Write error !");}
    return SetAny<pmeshS>(Thh);
    
}


class SaveMeshL :  public E_F0 { public:
    
    typedef pmeshL Result;
    Expression getmesh;
    Expression filename;
    //Expression xx,yy,zz;
    SaveMeshL(const basicAC_F0 & args)
    {
        args.SetNameParam();
        getmesh=to<pmeshL>(args[0]);
        filename=to<string*>(args[1]);
        
    }
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmeshL>(),atype<string*>(),true);}
    static  E_F0 * f(const basicAC_F0 & args){ return new SaveMeshL(args);}
    AnyType operator()(Stack s) const ;
    
};


// 3D line SaveMesh
AnyType SaveMeshL::operator()(Stack stack) const
{
    using  Fem2D::MeshPointStack;
    
    
    pmeshL Thh = GetAny<pmeshL>((*getmesh)(stack));
    string * fn =  GetAny<string*>((*filename)(stack));
    
    if (verbosity > 2)
        cout << "SaveMeshL " << *fn << " " << Thh << endl;
    int ret=Thh->Save(*fn);
    if( ret!=0) {ExecError("PB Write error !");}
    return SetAny<pmeshL>(Thh);
    
}


class SaveSurfaceMesh3 :  public E_F0 { public:  
 
    typedef pmesh  Result;
  Expression getmesh;
  Expression filename; 
  Expression filename1;  
  int pointsfaces;
  SaveSurfaceMesh3(const basicAC_F0 & args) 
  {   
    args.SetNameParam();
    getmesh=to<pmesh3>(args[0]); 
    filename=to<string*>(args[1]); 
    pointsfaces=0;
    if (args.size() >2) 
      {
	if(BCastTo<string *>(args[2])){ 
	  pointsfaces=1;
	  filename1=CastTo<string*>(args[2]);
	}  
	else{
	  CompileError("savesurfmesh(Th,filename.points,filename.faces)");
	}
      }
    
  }   
  static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmesh3>(),atype<string*>(),true);}
  static  E_F0 * f(const basicAC_F0 & args){ return new SaveSurfaceMesh3(args);} 
  AnyType operator()(Stack s) const ;
  
};


AnyType SaveSurfaceMesh3::operator()(Stack stack) const 
{
  using  Fem2D::MeshPointStack;
  
  
  pmesh3 Thh = GetAny<pmesh3>((*getmesh)(stack));
  string * fn =  GetAny<string*>((*filename)(stack));
  if( pointsfaces==0 ){
    cout << "SaveSurfaceMesh3 " << *fn << " " << Thh << endl;
    Thh->SaveSurface(*fn);
  }
  else{
    string * fn1 =  GetAny<string*>((*filename1)(stack));
    cout << "SaveSurfaceMesh3 " << *fn << " " << *fn1 << " " << Thh << endl;
    Thh->SaveSurface(*fn,*fn1); 
  }
  return SetAny<pmesh3>(Thh);
  
}


// version 3d of buildmeshborder
const MeshL* BuildMeshCurve3(Stack stack, E_BorderN const * const & b)   //  ,bool justboundary,int nbvmax=0,bool Requiredboundary       ,KNM<double> *pintern,double alea)
{
    int nbvinter=0;
    double precis_mesh=1.e-7;
 
    int brefintp= -2000000000;

    MeshPoint *mp(MeshPointStack(stack));
    mp->set();
    MeshPoint mps = *mp;
    
    int Gnbv=0,Gnbt=0,nbsd=0;    // vertice, edges, nb subdomains
    for (E_BorderN const *k=b;k;k=k->next) {
        int nbd = k->NbBorder(stack);
        for(int index=0;index<nbd;++index ) {
            long n=Max(1L,Abs(k->Nbseg(stack,index)));
            Gnbv+=n+1;
            Gnbt+=n;
            nbsd++;
        }
    }
    
    MeshL *Th =  new MeshL();
    
    if(verbosity>2)
        cout <<" Begin: MeshL from nb Border  "  << nbsd <<endl;

    Vertex3 *vertices =  new Vertex3[Gnbv];
    
    double lmin=0.;

    //  generation des points et des lignes
    long i=0,n=0;
    for (E_BorderN const * k=b;k;k=k->next) {
        int nbd = k->NbBorder(stack);
        for(int index=0; index<nbd; ++index ) {
            assert(k->b->xfrom);
            double & t = *  k->var(stack),tt;
            double a(k->from(stack)),b(k->to(stack));
            long * indx = k->index(stack);
            if(indx) *indx = index;
            else ffassert(index==0);
            n=Max(Abs(k->Nbseg(stack,index)),1L);
            tt=t=a;
            double delta = (b-a)/n;
            for ( int nn=0;nn<=n;++nn,++i, tt += delta) {
                t = tt;
                if (nn==n) t=b; // to remove roundoff error
                mp->label = k->label();
                k->code(stack); // compute x,y,z label
                vertices[i].x=mp->P.x;
                vertices[i].y=mp->P.y;
                vertices[i].z=mp->P.z;
                vertices[i].lab= mp->label;
               
                //vertices[i].color = i;
                if (nn>0) {
                    lmin=min(lmin,Norme2(vertices[i]-vertices[i-1]));
                }
            }
        }
    }
    lmin = sqrt(lmin);
    double eps = (lmin)/16.;
    int nbvprev = i;
    //long nbv=0;
    R3 bmin, bmax;
    bmin.x = vertices[0].x;
    bmin.y = vertices[0].y;
    bmin.z = vertices[0].z;
    bmax.x = bmin.x;
    bmax.y = bmin.y;
    bmax.z = bmin.z;

    // recherche des extrema des vertices pmin,pmax

    for (int ii = 1; ii < nbvprev; ii++) {
        bmin.x = min(bmin.x, vertices[ii].x);
        bmin.y = min(bmin.y, vertices[ii].y);
        bmin.z = min(bmin.z, vertices[ii].z);
        
        bmax.x = max(bmax.x, vertices[ii].x);
        bmax.y = max(bmax.y, vertices[ii].y);
        bmax.z = max(bmax.z, vertices[ii].z);
    }

    
    double longmini_box = pow(bmax.x - bmin.x, 2) + pow(bmax.y - bmin.y, 2) + pow(bmax.z - bmin.z, 2);
    longmini_box = sqrt(longmini_box);
    
    if (verbosity > 1) {
        cout << " bmin := " << bmin.x << " " << bmin.y << " " << bmin.z << endl;
        cout << " bmax := " << bmax.x << " " << bmax.y << " " << bmax.z << endl;
        cout << " box volume :=" << longmini_box << endl;
    }
    
    if (precis_mesh < 0)
        precis_mesh = longmini_box * 1e-7;
    else
        precis_mesh = precis_mesh;
    
    // determination de hmin
    double hmin =longmini_box;
   
    if (verbosity > 5)
        cout << "    Norme2(bmin-bmax)=" << Norme2(bmin - bmax) << endl;
    
    // assertion pour la taille de l octree
    assert(hmin > Norme2(bmin - bmax) / 1e9);
    
    
    double hseuil = hmin *1e-7;
    if (verbosity > 3)
        cout << "    hseuil=" << hseuil << endl;
    
    Vertex3 *vv = new Vertex3[nbvprev];
    EF23::GTree<Vertex3> *gtree = new EF23::GTree<Vertex3>(vv, bmin, bmax, 0);
    
    if (verbosity > 2) {
        cout << "  -- taille de la boite " << endl;
        cout << "\t" << bmin.x << " " << bmin.y << " " << bmin.z << endl;
        cout << "\t" << bmax.x << " " << bmax.y << " " << bmax.z << endl;
    }
    
    // creation of octree
    int new_nv=0;
    int *old2new=new int[nbvprev];
    for(int i = 0; i < nbvprev; ++i)
        old2new[i]=i;
    
    for (int ii = 0; ii < nbvprev; ++ii) {
        const Vertex3 &vi(vertices[ii]);
        Vertex3 *pvi = gtree->ToClose(vi, hseuil);
        if (!pvi) {
            vv[new_nv].x = vi.x;
            vv[new_nv].y = vi.y;
            vv[new_nv].z = vi.z;
            vv[new_nv].lab = vi.lab;
            old2new[ii] = new_nv;
            gtree->Add(vv[new_nv]);
            new_nv++;
        }
        else
            old2new[ii] = pvi - vv;
    }
    
   Th->nv = new_nv;
   Th->vertices = new Vertex3[new_nv];
   for(int ii=0;ii<new_nv;++ii) {
       Th->vertices[ii].x = vv[ii].x;
       Th->vertices[ii].y = vv[ii].y;
       Th->vertices[ii].z = vv[ii].z;
       Th->vertices[ii].lab = vv[ii].lab;
   }

   int nnn=0;
   i=0;
    
   Th->elements = new EdgeL[Gnbt];
    
   for (E_BorderN const * k=b;k;k=k->next) {
       int nbd = k->NbBorder(stack);
       for(int index=0; index<nbd; ++index ) {
           double & t = *  k->var(stack);
           double a(k->from(stack)),b(k->to(stack));
           n=Max(Abs(k->Nbseg(stack,index)),1L);
           long * indx = (k->index(stack));
           if(indx)
               *indx = index;
           else
               ffassert(index==0);
           double delta = (b-a)/n;
           t=a+delta/2;
           for (int nn=0;nn<n;nn++,i++, t += delta) {
               
               mp->label = k->label();
               k->code(stack);
               int iv[2];
               iv[0]=old2new[nnn];
               iv[1]=old2new[++nnn];
               Th->elements[i].set(Th->vertices,iv,mp->label);
           }
           nnn++;
       }
  }
  Th->nt=i;
  Th->nbe=0;  //Th->borderelements = new BoundaryPointL[0];
  
  delete [] vv;
  delete []vertices;
  delete gtree;
  delete [] old2new;
    
  Th->BuildBound();
  Th->BuildAdj();
  // Th->Buildbnormalv(); // no normal with 3d curve
  Th->BuildjElementConteningVertex();
  Th->BuildBdElem();
    
  Add2StackOfPtr2FreeRC(stack, Th);
  *mp=mps;
   
  return Th;
}


AnyType MoveMesh3::operator()(Stack stack) const 
{
  ffassert(0);
    /*
  using  Fem2D::Triangle;
  using  Fem2D::Vertex;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
 // using  Fem2D::R;
  using  Fem2D::MeshPointStack;
   MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
   Mesh * Thh = GetAny<pmesh>((*getmesh)(stack));
   ffassert(Thh);
   long nbv=Thh->nv;
   long nbt=Thh->nt;
   KN<double> u(nbv),v(nbv);
   double infini=DBL_MAX;
   u=infini;
   for (int it=0;it<nbt;it++)
    for (int iv=0;iv<3;iv++)
    {
      int i=(*Thh)(it,iv);
      if ( u[i]==infini) { // if nuset the set 
        mp->setP(Thh,it,iv);
        u[i]=GetAny<double>((*U)(stack));
        v[i]=GetAny<double>((*V)(stack));
      }
    }
    
   Mesh * pth= MoveTheMesh(*Thh,u,v);
   if (pth)
     for (size_t i=0;i<sol.size();i++)
       { //  ale 
          pair<FEbase<double>,int> * s = GetAny<pair<FEbase<double>,int>*>( (*sol[i])(stack));
          ffassert(s->first.Vh);
          ffassert( &s->first.Vh->Th == Thh); // same old mesh
          ffassert(0); // a faire ????
       }
   *mp=mps;
    pth->decrement();   
    return SetAny<pmesh>(pth);
    */
}

//3D volume
inline pmesh3 *  initMesh(pmesh3 * const & p, string * const & s) {
  Mesh3 * m;
  if(verbosity > 2)
      cout << " initMesh " << *s << endl;
    *p= m =new Mesh3(*s);
  m->BuildGTree();
 //  delete s;  modif mars 2006 auto del ptr 
  return p;
 }

//3D surface
inline pmeshS *  initMesh(pmeshS * const & p, string * const & s) {
    MeshS * m;
    if(verbosity > 2)
        cout << " initMesh " << *s << endl;
    *p= m =new MeshS(*s);
    m->BuildGTree();
    //  delete s;  modif mars 2006 auto del ptr
    return p;
}

//3D line
inline pmeshL *  initMesh(pmeshL * const & p, string * const & s) {
    MeshL * m;
    if(verbosity > 2)
        cout << " initMesh " << *s << endl;
    *p= m =new MeshL(*s);
    m->BuildGTree();
    //  delete s;  modif mars 2006 auto del ptr
    return p;
}

/*inline pmeshS *  initMesh(pmeshS * const & p, string * const & s) {
Mesh3 * m;
    MeshS * mS;
    if(verbosity > 2)
        cout << " initMesh " << *s << endl;
    
    m =new Mesh3(*s);
    *p=m->getMeshS();
    m->getMeshS()->BuildGTree();
    //  delete s;  modif mars 2006 auto del ptr
    return p;
}*/


/*
class CheckMoveMesh :  public E_F0mps { public:  
 
   typedef double  Result;
   Expression getmesh;
   Expression U,V;
   int nbsol;    
    vector<Expression> sol;
   
    CheckMoveMesh(const basicAC_F0 & args) :nbsol(args.size()-2),sol(args.size()-2)
    {   
  using  Fem2D::Triangle;
  using  Fem2D::Vertex;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
 // using  Fem2D::R;
  using  Fem2D::MeshPointStack;
      args.SetNameParam();
      getmesh=to<pmesh>(args[0]); 
      const E_Array * a = dynamic_cast<const E_Array *>(args[1].LeftValue());
      
      ffassert(a);
      if (a->size() !=2) CompileError("CheckMoveMesh(Th,[u,v]) need 2 componate in array ",atype<pmesh>());
      U=to<double>( (*a)[0]);
      V=to<double>( (*a)[1]);
      
      for (int i=2;i<args.size();i++)
        sol[i-2]=to<double>(args[i]);      
    }   
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmesh>(),atype<E_Array>(),false);}
    static  E_F0 * f(const basicAC_F0 & args){ return new CheckMoveMesh(args);} 
    AnyType operator()(Stack s) const ;
    operator aType () const { return atype<double>();}         
  
};
AnyType CheckMoveMesh::operator()(Stack stack) const 
{
 
  using  Fem2D::Triangle;
  using  Fem2D::Vertex;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
 // using  Fem2D::R;
  using  Fem2D::MeshPointStack;
   MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
   Mesh * Thh = GetAny<pmesh>((*getmesh)(stack));
   Mesh & Th(*Thh);
   ffassert(Thh);
   long nbv=Thh->nv;
   long nbt=Thh->nt;
   KN<double> u(nbv),v(nbv);
   double infini=DBL_MAX;
   u=infini;
   for (int it=0;it<nbt;it++)
    for (int iv=0;iv<3;iv++)
    {
      int i=(*Thh)(it,iv);
      if ( u[i]==infini) { // if nuset the set 
        mp->setP(Thh,it,iv);
        u[i]=GetAny<double>((*U)(stack));
        v[i]=GetAny<double>((*V)(stack));
      }
    }
     double minarea=DBL_MAX;
    for (int t=0;t<Th.nt;t++)
     {
      int i0=Th(t,0),i1=Th(t,1),i2=Th(t,2);
      minarea=Min(minarea,Area2(R2(u[i0],v[i0]), R2(u[i1],v[i1]),R2(u[i2],v[i2])));
     }
    *mp=mps;
    return SetAny<double>(minarea/2.);

}
*/

template<class R, class v_fes>
AnyType set_fe3 (Stack s,Expression ppfe, Expression e)
{ 
  typedef typename  v_fes::pfes pfes;
  typedef typename  v_fes::FESpace FESpace;
  typedef typename  FESpace::Mesh Mesh;
  typedef typename  FESpace::FElement FElement;
  typedef typename  Mesh::Element Element;
  typedef typename  Mesh::Vertex Vertex;  
  typedef typename  Mesh::RdHat RdHat;  
  typedef typename  Mesh::Rd Rd;  
  long kkff = Mesh::kfind,  kkth = Mesh::kthrough;
  StackOfPtr2Free * sptr = WhereStackOfPtr2Free(s);
  
  
  MeshPoint *mps=MeshPointStack(s),mp=*mps;  
  pair<FEbase<R,v_fes> *,int>  pp=GetAny<pair<FEbase<R,v_fes> *,int> >((*ppfe)(s));
  FEbase<R,v_fes> & fe(*pp.first);
  const  FESpace & Vh(*fe.newVh());
  KN<R> gg(Vh.MaximalNbOfDF()); 
  const  Mesh & Th(Vh.Th);
  //   R F[100]; // buffer 
  TabFuncArg tabexp(s,Vh.N);
  tabexp[0]=e;
  
  if(Vh.N!=1)
    {  cerr << " Try to set a  vectorial  FE function  (nb  componant=" <<  Vh.N << ") with one scalar " << endl;
       ExecError(" Error interploation (set)  FE function (vectorial) with a scalar");
    }
  KN<R> * y=new  KN<R>(Vh.NbOfDF);
  KN<R> & yy(*y);
  // interpoler
  int npPh = Vh.maxNbPtforInterpolation;
  KNM<R>   Vp(npPh,1);
  KN<R>  Vdf(Vh.MaxNbDFPerElement);
  const E_F0 & ff(* (const  E_F0 *) e ) ;
 //   cout << "Vh.isFEMesh() :" <<Vh.isFEMesh() << endl;
  if (Vh.isFEMesh() )
    {
      
      ffassert(Vh.NbOfDF == Th.nv && Vh.N == 1 );
      for (int iv=0;iv<Th.nv;iv++)
	{
	  const Vertex & v(Th(iv));
	  int ik=Th.Contening(&v);
	  const Element & K(Th[ik]);
	  int il=-1;
	  for(int k=0;k<Element::nv;++k)
	    if  ( &K[k] == &v) il=k;
	  assert(il>=0);
	  mps->setP(&Th,ik,il);
	  yy[iv] = GetAny<R>( ff(s) );
	  sptr->clean(); // modif FH mars 2006  clean Ptr
	}
      
    }
  else     
    {
      InterpolationMatrix<RdHat> ipmat(Vh);    
      for (int t=0;t<Th.nt;t++)
	{
	  FElement K(Vh[t]);
	  int nbdf=K.NbDoF() ;	  
	  ipmat.set(K);
	  for (int p=0;p<ipmat.np;p++)
	    { 
	      const RdHat & PtHat(ipmat.P[p]);
	      mps->set(K.T(PtHat),PtHat,K);
	      Vp[p]=GetAny<R>( ff(s) );
	    }
	  K.Pi_h(Vp,Vdf,ipmat);
	 //   cout << "Vp::: " << Vp << " " << Vdf << endl;
	  for (int df=0;df<nbdf;df++)         
	    (*y)[K(df)] =  Vdf[df] ;
	  sptr->clean(); // modif FH mars 2006  clean Ptr	
	}
    }
  *mps=mp;
  fe=y;
  kkff = Mesh::kfind - kkff;
  kkth = Mesh::kthrough -kkth;
  
  if(verbosity>1)
    ShowBound(*y,cout) 
      << " " << kkth << "/" << kkff << " =  " << double(kkth)/Max<double>(1.,kkff) << endl;
  return SetAny<FEbase<R,v_fes>*>(&fe); 
}

















template<class K,class v_fes>
E_set_fev3<K,v_fes>::E_set_fev3(const E_Array * a,Expression pp) 
  :aa(*a),ppfe(pp),optimize(true),
   where_in_stack_opt(),optiexp0(),optiexpK() 
   { 
     aa.map(to<K>) ;
     bool kdump=false;
     if(optimize)
       { // new code Optimized  -------
	 int n=aa.size();
	 deque<pair<Expression,int> > ll;
	 MapOfE_F0 m;
	 where_in_stack_opt.resize(n);
	 size_t top = currentblock->OffSet(0), topbb=top; // FH. bofbof ??? 
	 for (int i=0; i<n; i++)
	   {
	     Expression ee= aa[i].LeftValue();
	     if (kdump)
	       cout << "Optimize OneOperatorMakePtrFE:  type exp: " << typeid(*ee).name() << " "<<endl;
	     where_in_stack_opt[i]=ee->Optimize(ll, m, top);
	     if (kdump)
	       cout  << "\n\t\t"<< i  << ": " << where_in_stack_opt[i] << endl;
	   }
	 
	 currentblock->OffSet(top-topbb);
	 //  
	 int k=ll.size(),k0=0,k1=0;
	 for (int i=0;i<k;i++)
	   if (ll[i].first->MeshIndependent()) k0++;
	 deque<pair<Expression,int> > l0(k0),l1(k-k0);
	 k0=0,k1=0;
	 for (int i=0;i<k;i++)
	   if (ll[i].first->MeshIndependent()) 
	     {
	       if (kdump)
		 cout << " mi " << ll[i].second << " " << *(ll[i].first) << endl;
	       l0[k0++]=ll[i];
	     }
	   else 
	     {
	       if (kdump)
		 cout << " md " << ll[i].second << " " << *(ll[i].first) << endl;
	       l1[k1++]=ll[i];
	     }
	 if (k0)      
	   optiexp0 = new E_F0_Optimize(l0,m,0);  // constant part
	 if (k1) 
	   optiexpK = new E_F0_Optimize(l1,m,0);  // none constant part
	 
       }
     
     
   }


template<class K,class v_fes>   
AnyType E_set_fev3<K,v_fes>::operator()(Stack s)  const
{  
  StackOfPtr2Free * sptr = WhereStackOfPtr2Free(s);     
  MeshPoint *mps=MeshPointStack(s), mp=*mps;   
  FEbase<K,v_fes> ** pp=GetAny< FEbase<K,v_fes> **>((*ppfe)(s));
  FEbase<K,v_fes> & fe(**pp);
  const  FESpace & Vh(*fe.newVh());
 // KN<K> gg(Vh.MaximalNbOfDF()); 
  
  const  Mesh & Th(Vh.Th); 
  const int dim=Vh.N;
  K ** copt=0;
  if (optimize)   copt= new K *[dim];
  if(copt) {
    assert((size_t) dim== where_in_stack_opt.size());
    for (int i=0;i<dim;i++)
      {
        int offset=where_in_stack_opt[i];
        assert(offset>10);
        copt[i]= static_cast<K *>(static_cast<void *>((char*)s+offset));
        *(copt[i])=0;
      }
    if (optiexp0) (*optiexp0)(s); // init 
  }
  
  ffassert(dim<100);
  //   R F[100]; // buffer 
  
  TabFuncArg tabexp(s,Vh.N);
  //   const E_Array * aa = dynamic_cast<const E_Array *>(e);
  ffassert( aa.size() == Vh.N);
  for (int i=0;i<dim;i++)
    tabexp[i]=aa[i]; 
  
  KN<K> * y=new  KN<K>(Vh.NbOfDF);
  KN<K> & yy(*y);
  int npPh = Vh.maxNbPtforInterpolation;
  KNM<K>   Vp(npPh,dim);
  KN<K>  Vdf(Vh.MaxNbDFPerElement);
  
  if (Vh.isFEMesh() )
    {
      // crrect bug 29/08/2011 (thanks to rychet@fzu.cz)
      // remove wrong bulid of KHat (memory out of bound)
      ffassert(Vh.NbOfDF == Th.nv && dim == 1 );
      for (int iv=0;iv<Th.nv;iv++)
	{
	  const E_F0 & ff(* (const  E_F0 *) aa[0]  ) ;
	  const Vertex & v(Th(iv));
	  int ik=Th.Contening(&v);
	  const Element & Kv(Th[ik]);
          int il=-1;
          for(int k=0;k<Element::nv;++k)
            if  ( &Kv[k] == &v) il=k;
          assert(il>=0);
          mps->set(Th,v,RdHat::KHat[il],Kv,v.lab);
	  if (copt) {
	    if (optiexpK) (*optiexpK)(s); 
	    yy[iv] =  *(copt[0]);
	  }
	  else 
	    yy[iv] = GetAny<K>( ff(s) );
	  sptr->clean(); // modif FH mars 2006  clean Ptr
       }
    }
  else
    {
       InterpolationMatrix<RdHat> ipmat(Vh);    
       
       for (int t=0;t<Th.nt;t++)
	 {
	   FElement Kt(Vh[t]);
	   int nbdf=Kt.NbDoF();
	   
	   //gg=K();
	   
	   for (int p=0;p<ipmat.np;p++)
	     {
	       const RdHat & PtHat(ipmat.P[p]);
	       mps->set(Kt.T(PtHat),PtHat,Kt);
	       
	       if (copt) { // optimize  version 
		 if (optiexpK) (*optiexpK)(s);
		 for (int j=0;j<dim;j++)
		   Vp(p,j) = *(copt[j]);}
	       else  // old version 
		 for (int j=0;j<dim;j++)
		   if (tabexp[j]) 
		     Vp(p,j)=GetAny<K>( (*tabexp[j])(s) );
		   else Vp(p,j)=0;
	       
	     }
             if(verbosity>10000) cout << " E_set_fev3: " << ipmat << endl;
	   Kt.Pi_h(Vp,Vdf,ipmat);  
	   //  cout << "Vp --- "<< Vp  << "  Vdf;; " << Vdf[0] << endl;
	   for (int df=0;df<nbdf;df++)         
	     yy[Kt(df)] =  Vdf[df] ;
	   
	   sptr->clean(); // modif FH mars 2006  clean Ptr          
	 } 
    }
  fe=y;
  if (copt) delete [] copt;
  *MeshPointStack(s) = mp;
  if(verbosity>1)
    ShowBound(*y,cout) << endl ;
  return Nothing;
}

// 3D volume
template<class K>
inline FEbase<K,v_fes> * MakePtrFE3_(pfes3 * const &  a){ 
  FEbase<K,v_fes3> * p=new FEbase<K,v_fes3>(a);
  return p ;}
  
template<class K>
inline FEbase<K,v_fes3> ** MakePtrFE3_2(FEbase<K,v_fes3> * * const &  p,pfes3 * const &  a){ 
  *p=new FEbase<K,v_fes3>(a);
  return p ;}

template<class K>  
inline FEbaseArray<K,v_fes3> ** MakePtrFE3_3(FEbaseArray<K,v_fes3> * * const &  p,pfes3 * const &  a,const long & N){ 
  *p=new FEbaseArray<K,v_fes3>(a,N);
  return p ;}


// 3D surface
template<class K>
inline FEbase<K,v_fesS> * MakePtrFE3_(pfesS * const &  a){
    FEbase<K,v_fesS> * p=new FEbase<K,v_fesS>(a);
    return p ;}

template<class K>
inline FEbase<K,v_fesS> ** MakePtrFE3_2(FEbase<K,v_fesS> * * const &  p,pfesS * const &  a){
    *p=new FEbase<K,v_fesS>(a);
    return p ;}

template<class K>
inline FEbaseArray<K,v_fesS> ** MakePtrFE3_3(FEbaseArray<K,v_fesS> * * const &  p,pfesS * const &  a,const long & N){
    *p=new FEbaseArray<K,v_fesS>(a,N);
    return p ;}


// 3D curve
template<class K>
inline FEbase<K,v_fesL> * MakePtrFE3_(pfesL * const &  a){
    FEbase<K,v_fesL> * p=new FEbase<K,v_fesL>(a);
    return p ;}

template<class K>
inline FEbase<K,v_fesL> ** MakePtrFE3_2(FEbase<K,v_fesL> * * const &  p,pfesL * const &  a){
    *p=new FEbase<K,v_fesL>(a);
    return p ;}

template<class K>
inline FEbaseArray<K,v_fesL> ** MakePtrFE3_3(FEbaseArray<K,v_fesL> * * const &  p,pfesL * const &  a,const long & N){
    *p=new FEbaseArray<K,v_fesL>(a,N);
    return p ;}



template<class K,class v_fes>
class  OneOperatorMakePtrFE3 : public OneOperator 
{
public:
  typedef typename  v_fes::pfes pfes;
  typedef typename  v_fes::FESpace FESpace;
  typedef typename  FESpace::Mesh Mesh;
  typedef typename  FESpace::FElement FElement;
  typedef typename  Mesh::Element Element;
  typedef typename  Mesh::Vertex Vertex;  
  typedef typename  Mesh::RdHat RdHat;  
  typedef typename  Mesh::Rd Rd;  

  // il faut Optimize 
  // typedef double K;
  typedef  FEbase<K,v_fes> ** R;
  typedef pfes* B;
  class CODE : public E_F0mps  
  {
  public:
    Expression fer,fes;
    E_set_fev3<K,v_fes> * e_set_fev3;
    const E_Array * v;
    CODE(const basicAC_F0 & args) 
      : 
      fer(to<R>(args[0])),
      fes(to<B>(args[1])),
      e_set_fev3(0) 
    {
      if (BCastTo<K>(args[2]) )
	v = new E_Array(basicAC_F0_wa(to<K>(args[2]))); 
      else 
	v = dynamic_cast<const E_Array *>( args[2].LeftValue() );
      if (!v) {
	cout << "Error: type of arg :" << *args[2].left()  << " in " << typeid(K).name() << " case " << endl;
	ErrorCompile(" We wait  a double/complex expression or a array expression",1);
      }
      e_set_fev3=  new   E_set_fev3<K,v_fes>(v,fer);
      
    }
    
    AnyType operator()(Stack stack)  const {
      R  p = GetAny<R>( (*fer)(stack));
      B  a = GetAny<B>( (*fes)(stack)); 
      *p=new FEbase<K,v_fes>(a);
      (*e_set_fev3)(stack); 
      return SetAny<R>(p);
    }
    operator aType () const { return atype<R>();}         
    
    
  };
  
  E_F0 * code(const basicAC_F0 & args) const 
  { return  new CODE(args);}
  OneOperatorMakePtrFE3(aType tt):  // tt= aType<double>() or aType<E_Array>()  
    OneOperator(map_type[typeid(R).name()],map_type[typeid(R).name()],map_type[typeid(B).name()],tt)
  {}
};


/*
template<class K,class v_fes>    
KN<K> * pf3r2vect( pair<FEbase<K,v_fes> *,int> p)
 {  
    KN<K> * x=p.first->x();
    if ( !x) {  // defined 
      FESpace3 * Vh= p.first->newVh();     
      throwassert( Vh);
      *p.first = x = new KN<K>(Vh->NbOfDF);
      *x=K(); 
    }
    return x;}
*/
// 3D volume
template<class K>        
long pf3r_nbdf(pair<FEbase<K,v_fes3> *,int> p)
 {  
   if (!p.first->Vh) p.first->Vh= p.first->newVh();
   throwassert( !!p.first->Vh);
   return p.first->Vh->NbOfDF;
 }
template<class K>
pmesh3 pf3r_Th(pair<FEbase<K,v_fes3> *,int> p)
{
    if (!p.first->Vh) p.first->Vh= p.first->newVh();
    throwassert( !!p.first->Vh);
    return & p.first->Vh->Th;
}

// 3D surface
template<class K>
long pfSr_nbdf(pair<FEbase<K,v_fesS> *,int> p)
{
    if (!p.first->Vh) p.first->Vh= p.first->newVh();
    throwassert( !!p.first->Vh);
    return p.first->Vh->NbOfDF;
}
template<class K>
pmeshS pfSr_Th(pair<FEbase<K,v_fesS> *,int> p)
{
    if (!p.first->Vh) p.first->Vh= p.first->newVh();
    throwassert( !!p.first->Vh);
    return & p.first->Vh->Th;
}

// 3D curve
template<class K>
long pfLr_nbdf(pair<FEbase<K,v_fesL> *,int> p)
{
    if (!p.first->Vh) p.first->Vh= p.first->newVh();
    throwassert( !!p.first->Vh);
    return p.first->Vh->NbOfDF;
}
template<class K>
pmeshL pfLr_Th(pair<FEbase<K,v_fesL> *,int> p)
{
    if (!p.first->Vh) p.first->Vh= p.first->newVh();
    throwassert( !!p.first->Vh);
    return & p.first->Vh->Th;
}

//3D volume
long pVh3_ndof(pfes3 * p)
 { throwassert(p && *p);
   FESpace3 *fes=**p; ;  return fes->NbOfDF ;}
long pVh3_nt(pfes3 * p)
 { throwassert(p && *p);
   FESpace3 *fes=**p; ;  return fes->NbOfElements ;}
long pVh3_ndofK(pfes3 * p)
 { throwassert(p && *p);
   FESpace3 *fes=**p;   return (*fes)[0].NbDoF() ;}
pmesh3 pVh3_Th(pfes3 * p)
{ throwassert(p && *p);
    FESpace3 *fes=**p; ;  return &fes->Th ;}

//3D surface
long pVhS_ndof(pfesS * p)
{ throwassert(p && *p);
    FESpaceS *fes=**p; ;  return fes->NbOfDF ;}
long pVhS_nt(pfesS * p)
{ throwassert(p && *p);
    FESpaceS *fes=**p; ;  return fes->NbOfElements ;}
long pVhS_ndofK(pfesS * p)
{ throwassert(p && *p);
    FESpaceS *fes=**p;   return (*fes)[0].NbDoF() ;}
pmeshS pVhS_Th(pfesS * p)
{ throwassert(p && *p);
    FESpaceS *fes=**p; ;  return &fes->Th ;}

//3D curve
long pVhL_ndof(pfesL * p)
{ throwassert(p && *p);
    FESpaceL *fes=**p; ;  return fes->NbOfDF ;}
long pVhL_nt(pfesL * p)
{ throwassert(p && *p);
    FESpaceL *fes=**p; ;  return fes->NbOfElements ;}
long pVhL_ndofK(pfesL * p)
{ throwassert(p && *p);
    FESpaceL *fes=**p;   return (*fes)[0].NbDoF() ;}
pmeshL pVhL_Th(pfesL * p)
{ throwassert(p && *p);
    FESpaceL *fes=**p; ;  return &fes->Th ;}



template<class R,int dd,class v_fes>
AnyType pf3r2R(Stack s,const AnyType &a)
{
  typedef typename  v_fes::pfes pfes;
  typedef typename  v_fes::FESpace FESpace;
  typedef typename  FESpace::Mesh Mesh;
  typedef typename  FESpace::FElement FElement;
  typedef typename  Mesh::Element Element;
  typedef typename  Mesh::Vertex Vertex;  
  typedef typename  Mesh::RdHat RdHat;  
  typedef typename  Mesh::Rd Rd;  

  
  pair< FEbase<R,v_fes> *  ,int> ppfe=GetAny<pair< FEbase<R,v_fes> *,int> >(a);
  FEbase<R,v_fes> & fe( *ppfe.first);
  int componante=ppfe.second;
  //  cout << "+++++++++++ " << componante <<endl;
  if ( !fe.x()) {
    if ( !fe.x()){
      // CompileError(" Sorry unset fem array ");
      return   SetAny<R>(0.0);
    }
  }
  
  const FESpace & Vh(*fe.Vh);
  const Mesh & Th(Vh.Th);
  MeshPoint & mp = *MeshPointStack(s);
  const Element *K;
  RdHat PHat;
  bool outside=false;
  bool qnu=true;
  if ( mp.Th3 && mp.Th3->elements == Th.elements && mp.T)
   {
     qnu=false;
     K=mp.T3;
     PHat=mp.PHat;
   }
  else if ( mp.other.Th3 
            && (mp.other.Th3->elements ==  Th.elements)
            && (mp.other.P.x == mp.P.x) && (mp.other.P.y == mp.P.y) && (mp.other.P.z == mp.P.z)   )
    {
      K=mp.other.T3;
      PHat=mp.other.PHat;
      outside = mp.other.outside;
    } 
  else {
    if (mp.isUnset()) ExecError("Try to get unset x,y, ...");
    K=Th.Find(mp.P,PHat,outside);
    mp.other.set(Th,mp.P,PHat,*K,0,outside);
  }
  if(verbosity>1000)
    {
    if(outside)
	cout << "  ---  " << qnu << "  " << mp.P << " out=" << mp.outside <<  " out=" << outside << " K= " << K << " " << PHat << endl;
    else 
	cout << "  ---  " << qnu << " P=  " << mp.P << " out=" << mp.outside <<  " out=" << outside << " K(PHat) == P =  " <<  (*K)(PHat) << " PHat = " << PHat << endl;
    }
  const FElement KK(Vh[Th(K)]);
  if (outside && KK.tfe->discontinue)
  {   if(verbosity>=10000) cout << " outside f() = 0. " << KK.tfe->discontinue << "\n";
    return   SetAny<R>(0.0);
  }
#ifndef NDEBUG
  if (!outside) 
    {
      if ( Norme2_2( (*K)(PHat) - mp.P ) > 1e-12 )
      {
        cout << "bug ??  " << Norme2_2( (*K)(PHat) - mp.P ) << " , " << mp.P << " , " << (*K)(PHat) << endl;
          verbosity=100000;
        K=Th.Find(mp.P,PHat,outside);
        ffassert(0 );
      }
    }
    
#endif
/*  int nbdf=KK.NbDoF();
  
  int N= KK.N;
  KN_<R> U(*fe.x());
  KNMK<R> fb(nbdf,N,3); //  the value for basic fonction
  KN<R> fk(nbdf);
  for (int i=0;i<nbdf;i++) // get the local value
    fk[i] = U[KK(i)];
    //  get value of basic function
  KK.BF(PHat,fb);  
  
  R r = (fb('.',componante,dd),fk);  
*/
  
 const R rr = KK(PHat,*fe.x(),componante,dd);
//  cout << " " << rr << endl;
//  R2 B(mp.P);
/*   if ( r < 0.08001 &&  Norme2_2(mp.P) > 0.05  && Norme2_2(mp.P) < 0.4*0.4  ) 
     {
     int vv=verbosity;
      cout << " f()  triangle  " << Th(K) << " " << mp.P << " " << PHat << " =  " << r << " " <<outside ;
      if (outside) {  verbosity = 200;
         K=Th.Find(mp.P,PHat,outside);
          cout << Th(K) << " " << outside << endl;}
       cout << endl; verbosity=vv;
     } */
//  if ( qnu )   
 //  cout << " f()  triangle       " << Th(K) << " " << mp.P << " " << PHat << " =  " << r <<  endl;
 if(verbosity>=10000)
      cout << componante<< " "<< dd << " f()  Tet:  " << Th(K) << " " << mp.P << " " << PHat << " =  " << rr <<  endl; 
  return SetAny<R>(rr);
}


template<class R,int dd,class v_fes>
AnyType pfSr2R(Stack s,const AnyType &a)
{
  typedef typename  v_fes::pfes pfes;
  typedef typename  v_fes::FESpace FESpace;
  typedef typename  FESpace::Mesh Mesh;
  typedef typename  FESpace::FElement FElement;
  typedef typename  Mesh::Element Element;
  typedef typename  Mesh::Vertex Vertex;
  typedef typename  Mesh::RdHat RdHat;
  typedef typename  Mesh::Rd Rd;

  pair< FEbase<R,v_fes> *  ,int> ppfe=GetAny<pair< FEbase<R,v_fes> *,int> >(a);
  FEbase<R,v_fes> & fe( *ppfe.first);
  int componante=ppfe.second;
  if ( !fe.x()) return   SetAny<R>(0.0);
  const FESpace & Vh(*fe.Vh);
  const Mesh & Th(Vh.Th);
  MeshPoint & mp = *MeshPointStack(s);
  const Element *K;
  RdHat PHat;
  bool outside=false;
  bool qnu=true;
  if (mp.ThS && mp.ThS->elements == Th.elements && mp.TS)
   {
     qnu=false;
     K=mp.TS;
     PHat.x = mp.PHat.x;
     PHat.y = mp.PHat.y;
   }
  else if (mp.other.ThS
            && (mp.other.ThS->elements ==  Th.elements)
            && (mp.other.P.x == mp.P.x) && (mp.other.P.y == mp.P.y) && (mp.other.P.z == mp.P.z)   )
    {
      K=mp.other.TS;
      PHat.x = mp.other.PHat.x;
      PHat.y = mp.other.PHat.y;
      outside = mp.other.outside;
    }
  else {
    if (mp.isUnset()) ExecError("Try to get unset x,y, ...");
    K=Th.Find(mp.P,PHat,outside);
    mp.other.set(Th,mp.P,PHat,*K,0,outside);
  }
  if(verbosity>1000)
    {
    if(outside)
	cout << "  ---  " << qnu << "  " << mp.P << " out=" << mp.outside <<  " out=" << outside << " K= " << K << " " << PHat << endl;
    else
	cout << "  ---  " << qnu << " P=  " << mp.P << " out=" << mp.outside <<  " out=" << outside << " K(PHat) == P =  " <<  (*K)(PHat) << " PHat = " << PHat << endl;
    }
  const FElement KK(Vh[Th(K)]);
  if (outside && KK.tfe->discontinue)
  {   if(verbosity>=10000) cout << " outside f() = 0. " << KK.tfe->discontinue << "\n";
    return   SetAny<R>(0.0);
  }
#ifndef NDEBUG
  if (!outside) 
    {
      if ( Norme2_2( (*K)(PHat) - mp.P ) > K->mesure()/1000 )
      {
        cout << "bug ??  " << Norme2_2( (*K)(PHat) - mp.P ) << " " << mp.P << " " << (*K)(PHat) << endl;
       // ffassert(0);
      }
    }
#endif

 const R rr = KK(PHat,*fe.x(),componante,dd);

 if(verbosity>=10000)
      cout << componante<< " "<< dd << " f()  Tet:  " << Th(K) << " " << mp.P << " " << PHat << " =  " << rr <<  endl;
  return SetAny<R>(rr);
}

template<class R,int dd,class v_fes>
AnyType pfLr2R(Stack s,const AnyType &a)
{
    typedef typename  v_fes::pfes pfes;
    typedef typename  v_fes::FESpace FESpace;
    typedef typename  FESpace::Mesh Mesh;
    typedef typename  FESpace::FElement FElement;
    typedef typename  Mesh::Element Element;
    typedef typename  Mesh::Vertex Vertex;
    typedef typename  Mesh::RdHat RdHat;
    typedef typename  Mesh::Rd Rd;
    
    pair< FEbase<R,v_fes> *  ,int> ppfe=GetAny<pair< FEbase<R,v_fes> *,int> >(a);
    FEbase<R,v_fes> & fe( *ppfe.first);
    int componante=ppfe.second;
    if ( !fe.x()) {
        if ( !fe.x()){
            return   SetAny<R>(0.0);
        }
    }
    
    const FESpace & Vh(*fe.Vh);
    const Mesh & Th(Vh.Th);
    MeshPoint & mp = *MeshPointStack(s);
    const Element *K;
    R1 PHat;
    bool outside=false;
    bool qnu=true;
    if (mp.ThL && mp.ThL->elements == Th.elements && mp.TL)
    {
        qnu=false;
        K=mp.TL;
        PHat.x = mp.PHat.x;
    }
    else if (mp.other.ThL
             && (mp.other.ThL->elements ==  Th.elements)
             && (mp.other.P.x == mp.P.x) && (mp.other.P.y == mp.P.y) && (mp.other.P.z == mp.P.z)   )
    {
        K=mp.other.TL;
        PHat.x = mp.other.PHat.x;
        outside = mp.other.outside;
    }
    else {
        if (mp.isUnset()) ExecError("Try to get unset x,y, ...");
        K=Th.Find(mp.P,PHat,outside);
        mp.other.set(Th,mp.P,PHat,*K,0,outside);
    }
    if(verbosity>1000)
    {
        if(outside)
            cout << "  ---  " << qnu << "  " << mp.P << " out=" << mp.outside <<  " out=" << outside << " K= " << K << " " << PHat << endl;
        else
            cout << "  ---  " << qnu << " P=  " << mp.P << " out=" << mp.outside <<  " out=" << outside << " K(PHat) == P =  " <<  (*K)(PHat) << " PHat = " << PHat << endl;
    }
    const FElement KK(Vh[Th(K)]);
    if (outside && KK.tfe->discontinue)
    {   if(verbosity>=10000) cout << " outside f() = 0. " << KK.tfe->discontinue << "\n";
        return   SetAny<R>(0.0);
    }
#ifndef NDEBUG
    if (!outside)
    {
        if ( Norme2_2( (*K)(PHat) - mp.P ) > 1e-12 )
        {
            cout << "bug ??  " << Norme2_2( (*K)(PHat) - mp.P ) << " " << mp.P << " " << (*K)(PHat) << endl;
        ffassert(0);
        }
    }
#endif
    
    const R rr = KK(PHat,*fe.x(),componante,dd);
    
    if(verbosity>=10000)
        cout << componante<< " "<< dd << " f()  Tet:  " << Th(K) << " " << mp.P << " " << PHat << " =  " << rr <<  endl;
    return SetAny<R>(rr);
}

template<class K,class v_fes>
class Op4_pf32K : public quad_function<pair<FEbase<K,v_fes> *,int>,R,R,R,K> { public:
    
    
    class Op : public E_F0mps { public:
	Expression a,b,c,d;
      Op(Expression aa,Expression bb,Expression cc,Expression dd) 
	: a(aa),b(bb),c(cc),d(dd) {}       
      AnyType operator()(Stack s)  const 
      { 
	
	R xx(GetAny<R>((*b)(s)));
	R yy(GetAny<R>((*c)(s)));
	R zz(GetAny<R>((*d)(s)));
	MeshPoint & mp = *MeshPointStack(s),mps=mp;
	mp.set(xx,yy,zz);
	AnyType ret = pf3r2R<K,0,v_fes>(s,(*a)(s));
	mp=mps;
	return  ret;}
      
    };
};

template<class K,class v_fes>
class Op4_pfS2K : public quad_function<pair<FEbase<K,v_fes> *,int>,R,R,R,K> { public:
    
    
    class Op : public E_F0mps { public:
    Expression a,b,c,d;
      Op(Expression aa,Expression bb,Expression cc,Expression dd)
    : a(aa),b(bb),c(cc),d(dd) {}
      AnyType operator()(Stack s)  const
      {
    
    R xx(GetAny<R>((*b)(s)));
    R yy(GetAny<R>((*c)(s)));
    R zz(GetAny<R>((*d)(s)));
    MeshPoint & mp = *MeshPointStack(s),mps=mp;
    mp.set(xx,yy,zz);
    AnyType ret = pfSr2R<K,0,v_fes>(s,(*a)(s));
    mp=mps;
    return  ret;}
      
    };
};


template<class K,class v_fes>
class Op4_pfL2K : public quad_function<pair<FEbase<K,v_fes> *,int>,R,R,R,K> { public:
    
    
    class Op : public E_F0mps { public:
    Expression a,b,c,d;
      Op(Expression aa,Expression bb,Expression cc,Expression dd)
    : a(aa),b(bb),c(cc),d(dd) {}
      AnyType operator()(Stack s)  const
      {
    
    R xx(GetAny<R>((*b)(s)));
    R yy(GetAny<R>((*c)(s)));
    R zz(GetAny<R>((*d)(s)));
    MeshPoint & mp = *MeshPointStack(s),mps=mp;
    mp.set(xx,yy,zz);
    AnyType ret = pfLr2R<K,0,v_fes>(s,(*a)(s));
    mp=mps;
    return  ret;}
      
    };
};


template<class K,class v_fes>    
KN<K> * pf3r2vect( pair<FEbase<K,v_fes> *,int> p)
{  
  //  cout << "  pf3r2vect " << endl;
  typedef typename  v_fes::FESpace FESpace;
  KN<K> * x=p.first->x();
  if ( !x) {  // defined 
    FESpace * Vh= p.first->newVh();     
    throwassert( Vh);
    *p.first = x = new KN<K>(Vh->NbOfDF);
    *x=K(); 
  }
  return x;}

// add 19 jan 2009 FH

//3D volume
class pVh3_ndf : public ternary_function<pfes3 *,long,long,long> { public:
    
    
    class Op : public E_F0mps { public:
	Expression a,b,c;
	Op(Expression aa,Expression bb,Expression cc) : a(aa),b(bb),c(cc) {}       
	AnyType operator()(Stack s)  const 
        { 
	    pfes3 * p(GetAny<pfes3 *>((*a)(s)));
	    long  k(GetAny<long>((*b)(s)));
	    long  i(GetAny<long>((*c)(s)));
	    throwassert(p && *p);
	    FESpace3 *fes=**p;
	    throwassert(fes && k >=0 && k < fes->NbOfElements );
	    FESpace3::FElement K=(*fes)[k];
	    throwassert(i>=0 && i <K.NbDoF() );
	    long ret(K(i));
	    return  ret;
	}
	
    };
};

//3D surface
class pVhS_ndf : public ternary_function<pfesS *,long,long,long> { public:
    
    
    class Op : public E_F0mps { public:
        Expression a,b,c;
        Op(Expression aa,Expression bb,Expression cc) : a(aa),b(bb),c(cc) {}
        AnyType operator()(Stack s)  const
        {
            pfesS * p(GetAny<pfesS *>((*a)(s)));
            long  k(GetAny<long>((*b)(s)));
            long  i(GetAny<long>((*c)(s)));
            throwassert(p && *p);
            FESpaceS *fes=**p;
            throwassert(fes && k >=0 && k < fes->NbOfElements );
            FESpaceS::FElement K=(*fes)[k];
            throwassert(i>=0 && i <K.NbDoF() );
            long ret(K(i));
            return  ret;
        }
        
    };
};

//3D line
class pVhL_ndf : public ternary_function<pfesL *,long,long,long> { public:
    
    
    class Op : public E_F0mps { public:
        Expression a,b,c;
        Op(Expression aa,Expression bb,Expression cc) : a(aa),b(bb),c(cc) {}
        AnyType operator()(Stack s)  const
        {
            pfesL * p(GetAny<pfesL *>((*a)(s)));
            long  k(GetAny<long>((*b)(s)));
            long  i(GetAny<long>((*c)(s)));
            throwassert(p && *p);
            FESpaceL *fes=**p;
            throwassert(fes && k >=0 && k < fes->NbOfElements );
            FESpaceL::FElement K=(*fes)[k];
            throwassert(i>=0 && i <K.NbDoF() );
            long ret(K(i));
            return  ret;
        }
        
    };
};


/*  no trivial ....   FH jan 2009..  */
/*class Op4_Mesh32mp : public quad_function<pmesh3*,R,R,R,MeshPoint *> { public:
    class Op : public E_F0mps { public:
	Expression a,b,c,d;
	Op(Expression aa,Expression bb,Expression cc,Expression dd) : a(aa),b(bb),c(cc),d(dd) {}       
	AnyType operator()(Stack s)  const 
        { 
	    R xx(GetAny<R>((*b)(s)));
	    R yy(GetAny<R>((*c)(s)));
	    R zz(GetAny<R>((*d)(s)));
	    pmesh3 *ppTh(GetAny<pmesh3*>((*a)(s)));
	    if( !ppTh || !*ppTh) ExecError("Op3_Mesh32mp unset mesh ??");
	    pmesh3 pTh(*ppTh);
	    MeshPoint * mp = new MeshPoint();
	    mp->set(xx,yy,zz);
	    R3 PHat;
	    bool outside;
	    const Tet * K=pTh->Find(mp->P,PHat,outside);
	    int n=(*pTh)(K);
	    if(verbosity>200)
	    cout << " n Tet = " << n << " " << K << " " <<  mp->P <<  endl;
	    mp->set(*pTh,mp->P,PHat,*K,0,outside);
	return mp;}
	
    };
};*/


// FH

// add Feb. 2010 FH 
template<class A> inline AnyType DestroyKN(Stack,const AnyType &x){
    KN<A> * a=GetAny<KN<A>*>(x);
    for (int i=0;i<a->N(); i++)
	(*a)[i]->destroy();
    a->destroy(); 
    return  Nothing;
}
template<class RR,class A,class B>  
RR * get_elementp_(const A & a,const B & b){ 
    if( b<0 || a->N() <= b) 
      { cerr << " Out of bound  0 <=" << b << " < "  << a->N() << " array type = " << typeid(A).name() << endl;
	  ExecError("Out of bound in operator []");}
    return  &((*a)[b]);}

template<class R>  R * set_initinit( R* const & a,const long & n){ 
    SHOWVERB( cout << " set_init " << typeid(R).name() << " " << n << endl);
    a->init(n);
    for (int i=0;i<n;i++)
	(*a)[i]=0;
    return a;}

//3d Volume
void init_mesh3_array()
{
    Dcl_Type<KN<pmesh3> *>(0,::DestroyKN<pmesh3> );
    atype<KN<pmesh3>* >()->Add("[","",new OneOperator2_<pmesh3*,KN<pmesh3>*,long >(get_elementp_<pmesh3,KN<pmesh3>*,long>));
    TheOperators->Add("<-", 
		      new OneOperator2_<KN<pmesh3> *,KN<pmesh3> *,long>(&set_initinit));
    map_type_of_map[make_pair(atype<long>(),atype<pmesh3>())]=atype<KN<pmesh3>*>(); // vector
    
    Dcl_Type< Resize<KN<pmesh3> > > ();
    Add<KN<pmesh3>*>("resize",".",new OneOperator1< Resize<KN<pmesh3> >,KN<pmesh3>*>(to_Resize));
    Add<Resize<KN<pmesh3> > >("(","",new OneOperator2_<KN<pmesh3> *,Resize<KN<pmesh3> > , long   >(resizeandclean1));
    
    
}

//3D surface
void init_meshS_array()
{
    Dcl_Type<KN<pmeshS> *>(0,::DestroyKN<pmeshS> );
    atype<KN<pmeshS>* >()->Add("[","",new OneOperator2_<pmeshS*,KN<pmeshS>*,long >(get_elementp_<pmeshS,KN<pmeshS>*,long>));
    TheOperators->Add("<-",
                      new OneOperator2_<KN<pmeshS> *,KN<pmeshS> *,long>(&set_initinit));
    map_type_of_map[make_pair(atype<long>(),atype<pmeshS>())]=atype<KN<pmeshS>*>(); // vector
    
    Dcl_Type< Resize<KN<pmeshS> > > ();
    Add<KN<pmeshS>*>("resize",".",new OneOperator1< Resize<KN<pmeshS> >,KN<pmeshS>*>(to_Resize));
    Add<Resize<KN<pmeshS> > >("(","",new OneOperator2_<KN<pmeshS> *,Resize<KN<pmeshS> > , long   >(resizeandclean1));
  
    
}


//3D Line
void init_meshL_array()
{
    Dcl_Type<KN<pmeshL> *>(0,::DestroyKN<pmeshL> );
    atype<KN<pmeshL>* >()->Add("[","",new OneOperator2_<pmeshL*,KN<pmeshL>*,long >(get_elementp_<pmeshL,KN<pmeshL>*,long>));
    TheOperators->Add("<-",
                      new OneOperator2_<KN<pmeshL> *,KN<pmeshL> *,long>(&set_initinit));
    map_type_of_map[make_pair(atype<long>(),atype<pmeshL>())]=atype<KN<pmeshL>*>(); // vector
    
    Dcl_Type< Resize<KN<pmeshL> > > ();
    Add<KN<pmeshL>*>("resize",".",new OneOperator1< Resize<KN<pmeshL> >,KN<pmeshL>*>(to_Resize));
    Add<Resize<KN<pmeshL> > >("(","",new OneOperator2_<KN<pmeshL> *,Resize<KN<pmeshL> > , long   >(resizeandclean1));
    
    
}


template<class K,class v_fes>
bool pfer_refresh3(pair<FEbase<K,v_fes> *,int> p)
{
    typedef typename v_fes::FESpace FESpace;
    int n = 0;
    if( p.first->x() ) n = p.first->x()->N();
    const FESpace * Vho = p.first->Vh;
    p.first->Vh= p.first->newVh();
    if ( n && n !=  p.first->Vh->NbOfDF)
        * p.first  = new KN<K>(p.first->Vh->NbOfDF,K());
    return Vho != (const FESpace * ) p.first->Vh;// FE change !!!
}

template<class Mesh>
KN_<long> listoflabel(Stack s, Mesh const * const & pTh)
{
    set<long> sl;
    for(int e=0; e< pTh->nbBrdElmts();++e)
        sl.insert(pTh->be(e).lab);
    
    long n = sl.size(),ii=0;
    long * p = Add2StackOfPtr2FreeA<long>(s,new long[n]); //   mark to be delete ..
    KN_<long> A(p,n);
    for(auto  i=sl.begin();i!=sl.end();++i)
        A[ii++]= *i;
    ffassert(ii== n);
  if(verbosity>99)   cout << " -- listoflabels "<< A << endl;
    return A;
};
template<class Mesh>
KN_<long> listofregion(Stack s, Mesh  const * const  & pTh)
{
    const Mesh & Th=*pTh;
    set<long> sl;
    for(int t=0; t< Th.nt;++t)
        sl.insert(Th[t].lab);
    
    long n = sl.size(),ii=0;
    long * p = Add2StackOfPtr2FreeA<long>(s,new long[n]); //   mark to be delete ..
    KN_<long> A(p,n);
    for(auto  i=sl.begin();i!=sl.end();++i)
        A[ii++]= *i;
     if(verbosity>99)   cout << " -- listofregions "<< A << endl;
    return A;
};
/*
// to see write u(x,y,z) 31 jan 2020 FH...
// ZZZ
template< class K,class v_fes >
class Op4_pfeK : public quad_function< pair< FEbase< K, v_fesL > *, int >, R, R,R, K > {
public:
    class Op : public E_F0mps {
    public:
        Expression a, b, c,d;
        Op(Expression aa, Expression bb, Expression cc,Expression dd)
        : a(aa), b(bb), c(cc), d(dd) {
        }
        AnyType operator( )(Stack s) const {
            R xx(GetAny< R >((*b)(s)));
            R yy(GetAny< R >((*c)(s)));
            R zz(GetAny< R >((*d)(s)));
            MeshPoint &mp = *MeshPointStack(s), mps = mp;
            mp.set(xx, yy, zz);
            AnyType ret = std::is_same<v_fes, v_fesS>::value ?
              pfSr2R< K, 0,v_fesS >(s, (*a)(s)) :pfLr2R< K, 0,v_fesL >(s, (*a)(s))  ;
            mp = mps;
            return ret;
        }
    };
};
*/

void init_lgmesh3() {
  if(verbosity&&(mpirank==0))  cout <<"lg_mesh3 ";

    //   Global.Add("buildmesh","(",new OneOperatorCode<classBuildMesh3>);
    // Global.Add("buildmesh","(",new OneOperatorCode<BuildMeshFile3>);

  //3D volume
  atype<pmesh3>()->AddCast( new E_F1_funcT<pmesh3,pmesh3*>(UnRef<pmesh3 >));
  atype<pfes3 >()->AddCast(  new E_F1_funcT<pfes3,pfes3*>(UnRef<pfes3>));
  atype<pf3rbase>()->AddCast(  new E_F1_funcT<pf3rbase,pf3rbase>(UnRef<pf3rbase>));
  atype<pf3cbase>()->AddCast(  new E_F1_funcT<pf3cbase,pf3cbase>(UnRef<pf3cbase>));
  // 3D surface
  atype<pmeshS>()->AddCast( new E_F1_funcT<pmeshS,pmeshS*>(UnRef<pmeshS >));
  atype<pfesS >()->AddCast(  new E_F1_funcT<pfesS,pfesS*>(UnRef<pfesS>));
  atype<pfSrbase>()->AddCast(  new E_F1_funcT<pfSrbase,pfSrbase>(UnRef<pfSrbase>));
  atype<pfScbase>()->AddCast(  new E_F1_funcT<pfScbase,pfScbase>(UnRef<pfScbase>));
  // 3D curve
  atype<pmeshL>()->AddCast( new E_F1_funcT<pmeshL,pmeshL*>(UnRef<pmeshL>));
  atype<pfesL >()->AddCast(  new E_F1_funcT<pfesL,pfesL*>(UnRef<pfesL>));
  atype<pfLrbase>()->AddCast(  new E_F1_funcT<pfLrbase,pfLrbase>(UnRef<pfLrbase>));
  atype<pfLcbase>()->AddCast(  new E_F1_funcT<pfLcbase,pfLcbase>(UnRef<pfLcbase>));
  //3D volume
  Add<pf3r>("[]",".",new OneOperator1<KN<double> *,pf3r>(pf3r2vect<R,v_fes3>));
  Add<pf3c>("[]",".",new OneOperator1<KN<Complex> *,pf3c>(pf3r2vect<Complex,v_fes3>));
  Add<pf3r>("(","",new OneQuadOperator<Op4_pf32K<R,v_fes3>,Op4_pf32K<R,v_fes3>::Op> );
  Add<pf3c>("(","",new OneQuadOperator<Op4_pf32K<Complex,v_fes3>,Op4_pf32K<Complex,v_fes3>::Op> );
  //3D surface
  Add<pfSr>("[]",".",new OneOperator1<KN<double> *,pfSr>(pf3r2vect<R,v_fesS>));
  Add<pfSc>("[]",".",new OneOperator1<KN<Complex> *,pfSc>(pf3r2vect<Complex,v_fesS>));
  Add<pfSr>("(","",new OneQuadOperator<Op4_pfS2K<R,v_fesS>,Op4_pfS2K<R,v_fesS>::Op> );
  Add<pfSc>("(","",new OneQuadOperator<Op4_pfS2K<Complex,v_fesS>,Op4_pfS2K<Complex,v_fesS>::Op> );
    
  //3D curve
  Add<pfLr>("[]",".",new OneOperator1<KN<double> *,pfLr>(pf3r2vect<R,v_fesL>));
  Add<pfLc>("[]",".",new OneOperator1<KN<Complex> *,pfLc>(pf3r2vect<Complex,v_fesL>));
  Add<pfLr>("(","",new OneQuadOperator<Op4_pfL2K<R,v_fesL>,Op4_pfL2K<R,v_fesL>::Op> );
  Add<pfLc>("(","",new OneQuadOperator<Op4_pfL2K<Complex,v_fesL>,Op4_pfL2K<Complex,v_fesL>::Op> );
    
  Add<double>("(","",new OneQuadOperator<Op4_K2R<R>,Op4_K2R<R>::Op> );
  // Add<long>("(","",new OneTernaryOperator<Op3_K2R<long>,Op3_K2R<long>::Op> ); // FH stupide
  Add<Complex>("(","",new OneQuadOperator<Op4_K2R<Complex>,Op4_K2R<Complex>::Op> );

 // Th(x,y,z)  first case ... remove  old case
 //  Add<pmesh3 *>("(","",new OneQuadOperator<Op4_Mesh32mp,Op4_Mesh32mp::Op> );
 // Add<pmesh3 *>("(","",new OneQuadOperator<Op4_Mesh32mp,Op4_Mesh32mp::Op> );
    
    

  //3D volume
  TheOperators->Add("<-", new OneOperator2_<pmesh3*,pmesh3*,string* >(&initMesh));
  // use for :   mesh Th = readmesh ( ...);
  TheOperators->Add("<-", new OneOperator2_<pmesh3*,pmesh3*,pmesh3 >(&set_copy_incr));
  TheOperators->Add("=", new OneOperator2<pmesh3*,pmesh3*,pmesh3 >(&set_eqdestroy_incr));
  //3D surface
  TheOperators->Add("<-", new OneOperator2_<pmeshS*,pmeshS*,string* >(&initMesh));
  // use for :   mesh Th = readmesh ( ...);
  TheOperators->Add("<-", new OneOperator2_<pmeshS*,pmeshS*,pmeshS >(&set_copy_incr));
  TheOperators->Add("=", new OneOperator2<pmeshS*,pmeshS*,pmeshS >(&set_eqdestroy_incr));
  //3D line
  TheOperators->Add("<-", new OneOperator2_<pmeshL*,pmeshL*,string* >(&initMesh));
  // use for :   mesh Th = readmesh ( ...);
  TheOperators->Add("<-", new OneOperator2_<pmeshL*,pmeshL*,pmeshL >(&set_copy_incr));
  TheOperators->Add("=", new OneOperator2<pmeshL*,pmeshL*,pmeshL >(&set_eqdestroy_incr));

  Global.Add("readmesh3","(",new OneOperatorCode<ReadMesh3>);
  Global.Add("readmeshS","(",new OneOperatorCode<ReadMeshS>);
  Global.Add("readmeshL","(",new OneOperatorCode<ReadMeshL>);
  Global.Add("savemesh","(",new OneOperatorCode<SaveMesh3>);
  Global.Add("savemesh","(",new OneOperatorCode<SaveMeshS>);
  Global.Add("savemesh","(",new OneOperatorCode<SaveMeshL>);
  Global.Add("savesurfacemesh","(",new OneOperatorCode<SaveSurfaceMesh3>);
    
  Global.Add("buildmeshL","(",new OneOperator1s_<pmeshL,const E_BorderN *>(BuildMeshCurve3));

    
  // 3D volume
  Dcl_Type<GlgVertex<Mesh3> >();
  Dcl_Type<GlgElement<Mesh3> >( );
  Dcl_Type<GlgElement<Mesh3>::Adj>( );
    
  Dcl_Type<GlgBoundaryElement<Mesh3>::BE >( );
  Dcl_Type<GlgBoundaryElement<Mesh3> >( );
    
  atype<long>()->AddCast(
		  new E_F1_funcT<long,GlgVertex<Mesh3> >(Cast<long,GlgVertex<Mesh3> >),
		  new E_F1_funcT<long,GlgElement<Mesh3> >(Cast<long,GlgElement<Mesh3> >),
		  new E_F1_funcT<long,GlgBoundaryElement<Mesh3> >(Cast<long,GlgBoundaryElement<Mesh3> >)

  );

 Add<pmesh3>("[","",new OneOperator2_<GlgElement<Mesh3>,pmesh3,long>(get_element));
 Add<pmesh3*>("[","",new OneOperator2_<GlgElement<Mesh3>,pmesh3*,long>(get_element));
 Add<pmesh3>("(","",new OneOperator2_<GlgVertex<Mesh3>,pmesh3,long>(get_vertex));
 Add<pmesh3*>("(","",new OneOperator2_<GlgVertex<Mesh3>,pmesh3*,long>(get_vertex));
//   second case of Th(x,y,z) ???? FH  attention deja fait dans Op4_Mesh32mp now remove
 Add<pmesh3* >("(", "", new OneQuadOperator< Op3_MeshDmp<Mesh3>, Op3_MeshDmp<Mesh3>::Op >);
    
 Add<pmesh3*>("be",".",new OneOperator1_<GlgBoundaryElement<Mesh3>::BE,pmesh3*>(Build));
 Add<GlgElement<Mesh3> >("adj",".",new OneOperator1_<GlgElement<Mesh3>::Adj,GlgElement<Mesh3> >(Build));
 Add<GlgBoundaryElement<Mesh3>::BE>("(","",new OneOperator2_<GlgBoundaryElement<Mesh3>,GlgBoundaryElement<Mesh3>::BE,long>(get_element));
 Add<GlgElement<Mesh3>::Adj>("(","",new OneOperator2_<GlgElement<Mesh3>,GlgElement<Mesh3>::Adj,long*>(get_adj));
 TheOperators->Add("==", new OneBinaryOperator<Op2_eq<GlgElement<Mesh3>,GlgElement<Mesh3> > >);
 TheOperators->Add("!=", new OneBinaryOperator<Op2_ne<GlgElement<Mesh3>,GlgElement<Mesh3> > >);
 TheOperators->Add("<", new OneBinaryOperator<Op2_lt<GlgElement<Mesh3>,GlgElement<Mesh3> > >);
 TheOperators->Add("<=", new OneBinaryOperator<Op2_le<GlgElement<Mesh3>,GlgElement<Mesh3> > >);

 Add<GlgBoundaryElement<Mesh3> >("label",".",new OneOperator1_<long,GlgBoundaryElement<Mesh3> >(getlab));
 Add<GlgBoundaryElement<Mesh3> >("Element",".",new OneOperator1_<GlgElement<Mesh3> ,GlgBoundaryElement<Mesh3> >(getElement));
 Add<GlgBoundaryElement<Mesh3> >("whoinElement",".",new OneOperator1_<long,GlgBoundaryElement<Mesh3> >(NuElement));
   

 Add<GlgElement<Mesh3> >("[","",new OneOperator2_<GlgVertex<Mesh3> ,GlgElement<Mesh3> ,long>(get_element));
 Add<GlgBoundaryElement<Mesh3> >("[","",new OneOperator2_<GlgVertex<Mesh3>,GlgBoundaryElement<Mesh3>,long>(get_element));
 
 Add<GlgVertex<Mesh3> >("x",".",new OneOperator1_<R,GlgVertex<Mesh3> >(getx));
 Add<GlgVertex<Mesh3> >("y",".",new OneOperator1_<R,GlgVertex<Mesh3> >(gety));
 Add<GlgVertex<Mesh3> >("z",".",new OneOperator1_<R,GlgVertex<Mesh3> >(getz));
 Add<GlgVertex<Mesh3> >("label",".",new OneOperator1_<long,GlgVertex<Mesh3> >(getlab));
 Add<GlgElement<Mesh3> >("label",".",new OneOperator1_<long,GlgElement<Mesh3> >(getlab));
 Add<GlgElement<Mesh3> >("region",".",new OneOperator1_<long,GlgElement<Mesh3> >(getlab));
 Add<GlgElement<Mesh3> >("mesure",".",new OneOperator1_<double,GlgElement<Mesh3> >(getmes));
 Add<GlgElement<Mesh3> >("measure",".",new OneOperator1_<double,GlgElement<Mesh3> >(getmes));
 Add<pmesh3*>("mesure",".",new OneOperator1<double,pmesh3*>(pmesh_mes));
 Add<pmesh3*>("measure",".",new OneOperator1<double,pmesh3*>(pmesh_mes));
 Add<pmesh3*>("bordermesure",".",new OneOperator1<double,pmesh3*>(pmesh_mesb));
 Add<pmesh3*>("bordermeasure",".",new OneOperator1<double,pmesh3*>(pmesh_mesb));
 Add<pmesh3*>("nt",".",new OneOperator1<long,pmesh3*>(pmesh_nt));
 Add<pmesh3*>("nv",".",new OneOperator1<long,pmesh3*>(pmesh_nv));
 Add<pmesh3*>("nbe",".",new OneOperator1<long,pmesh3*>(pmesh_nbe));
 Add<pmesh3*>("hmax",".",new OneOperator1<double,pmesh3*>(pmesh_hmax));
 Add<pmesh3*>("hmin",".",new OneOperator1<double,pmesh3*>(pmesh_hmin));

 Add<pmesh3*>("Gamma",".",new OneOperator1s_<pmeshS,pmesh3*>(pmesh3_gamma));
 Add<pmesh3*>("nbnomanifold",".",new OneOperator1<long,pmesh3*>(pmesh_nadjnomanifold));
 
 Add<pmeshS*>("Gamma",".",new OneOperator1s_<pmeshL,pmeshS*>(pmeshS_gamma));
    
 //3D surface
 Dcl_Type<GlgVertex<MeshS> >();
 Dcl_Type<GlgElement<MeshS> >( );
 Dcl_Type<GlgElement<MeshS>::Adj>( );
    
 Dcl_Type<GlgBoundaryElement<MeshS>::BE >( );
 Dcl_Type<GlgBoundaryElement<MeshS> >( );
    
 atype<long>()->AddCast(
                        new E_F1_funcT<long,GlgVertex<MeshS> >(Cast<long,GlgVertex<MeshS> >),
                        new E_F1_funcT<long,GlgElement<MeshS> >(Cast<long,GlgElement<MeshS> >),
                        new E_F1_funcT<long,GlgBoundaryElement<MeshS> >(Cast<long,GlgBoundaryElement<MeshS> >)
                       );
    
 Add<pmeshS>("[","",new OneOperator2_<GlgElement<MeshS>,pmeshS,long>(get_element));
 Add<pmeshS*>("[","",new OneOperator2_<GlgElement<MeshS>,pmeshS*,long>(get_element));
 Add<pmeshS>("(","",new OneOperator2_<GlgVertex<MeshS>,pmeshS,long>(get_vertex));
 Add<pmeshS*>("(","",new OneOperator2_<GlgVertex<MeshS>,pmeshS*,long>(get_vertex));
 Add<pmeshS* >("(", "", new OneQuadOperator< Op3_MeshDmp<MeshS>, Op3_MeshDmp<MeshS>::Op >);
    
 Add<pmeshS*>("be",".",new OneOperator1_<GlgBoundaryElement<MeshS>::BE,pmeshS*>(Build));
 Add<GlgElement<MeshS> >("adj",".",new OneOperator1_<GlgElement<MeshS>::Adj,GlgElement<MeshS> >(Build));
 Add<GlgBoundaryElement<MeshS>::BE>("(","",new OneOperator2_<GlgBoundaryElement<MeshS>,GlgBoundaryElement<MeshS>::BE,long>(get_element));
 Add<GlgElement<MeshS>::Adj>("(","",new OneOperator2_<GlgElement<MeshS>,GlgElement<MeshS>::Adj,long*>(get_adj));
 TheOperators->Add("==", new OneBinaryOperator<Op2_eq<GlgElement<MeshS>,GlgElement<MeshS> > >);
 TheOperators->Add("!=", new OneBinaryOperator<Op2_ne<GlgElement<MeshS>,GlgElement<MeshS> > >);
 TheOperators->Add("<", new OneBinaryOperator<Op2_lt<GlgElement<MeshS>,GlgElement<MeshS> > >);
 TheOperators->Add("<=", new OneBinaryOperator<Op2_le<GlgElement<MeshS>,GlgElement<MeshS> > >);
 
 Add<GlgBoundaryElement<MeshS> >("label",".",new OneOperator1_<long,GlgBoundaryElement<MeshS> >(getlab));
 Add<GlgBoundaryElement<MeshS> >("Element",".",new OneOperator1_<GlgElement<MeshS> ,GlgBoundaryElement<MeshS> >(getElement));
 Add<GlgBoundaryElement<MeshS> >("whoinElement",".",new OneOperator1_<long,GlgBoundaryElement<MeshS> >(NuElement));
  
 Add<GlgElement<MeshS> >("[","",new OneOperator2_<GlgVertex<MeshS> ,GlgElement<MeshS> ,long>(get_element));
 Add<GlgBoundaryElement<MeshS> >("[","",new OneOperator2_<GlgVertex<MeshS>,GlgBoundaryElement<MeshS>,long>(get_element));
  
 Add<GlgVertex<MeshS> >("x",".",new OneOperator1_<R,GlgVertex<MeshS> >(getx));
 Add<GlgVertex<MeshS> >("y",".",new OneOperator1_<R,GlgVertex<MeshS> >(gety));
 Add<GlgVertex<MeshS> >("z",".",new OneOperator1_<R,GlgVertex<MeshS> >(getz));
 Add<GlgVertex<MeshS> >("label",".",new OneOperator1_<long,GlgVertex<MeshS> >(getlab));
 Add<GlgElement<MeshS> >("label",".",new OneOperator1_<long,GlgElement<MeshS> >(getlab));
 Add<GlgElement<MeshS> >("region",".",new OneOperator1_<long,GlgElement<MeshS> >(getlab));
 Add<GlgElement<MeshS> >("mesure",".",new OneOperator1_<double,GlgElement<MeshS> >(getmes));
 Add<GlgElement<MeshS> >("measure",".",new OneOperator1_<double,GlgElement<MeshS> >(getmes));
 Add<pmeshS*>("mesure",".",new OneOperator1<double,pmeshS*>(pmesh_mes));
 Add<pmeshS*>("measure",".",new OneOperator1<double,pmeshS*>(pmesh_mes));
 Add<pmeshS*>("bordermesure",".",new OneOperator1<double,pmeshS*>(pmesh_mesb));
 Add<pmeshS*>("bordermeasure",".",new OneOperator1<double,pmeshS*>(pmesh_mesb));
 Add<pmeshS*>("nt",".",new OneOperator1<long,pmeshS*>(pmesh_nt));
 Add<pmeshS*>("nv",".",new OneOperator1<long,pmeshS*>(pmesh_nv));
 Add<pmeshS*>("nbe",".",new OneOperator1<long,pmeshS*>(pmesh_nbe));
 Add<pmeshS*>("hmax",".",new OneOperator1<double,pmeshS*>(pmesh_hmax));
 Add<pmeshS*>("hmin",".",new OneOperator1<double,pmeshS*>(pmesh_hmin));
 Add<pmeshS*>("nbnomanifold",".",new OneOperator1<long,pmeshS*>(pmesh_nadjnomanifold));
  
    
 //3D line
 Dcl_Type<GlgVertex<MeshL> >();
 Dcl_Type<GlgElement<MeshL> >( );
 Dcl_Type<GlgElement<MeshL>::Adj>( );
    
 Dcl_Type<GlgBoundaryElement<MeshL>::BE >( );
 Dcl_Type<GlgBoundaryElement<MeshL> >( );
    
 atype<long>()->AddCast(
                        new E_F1_funcT<long,GlgVertex<MeshL> >(Cast<long,GlgVertex<MeshL> >),
                        new E_F1_funcT<long,GlgElement<MeshL> >(Cast<long,GlgElement<MeshL> >),
                        new E_F1_funcT<long,GlgBoundaryElement<MeshL> >(Cast<long,GlgBoundaryElement<MeshL> >)
                        );
    
 Add<pmeshL>("[","",new OneOperator2_<GlgElement<MeshL>,pmeshL,long>(get_element));
 Add<pmeshL*>("[","",new OneOperator2_<GlgElement<MeshL>,pmeshL*,long>(get_element));
 Add<pmeshL>("(","",new OneOperator2_<GlgVertex<MeshL>,pmeshL,long>(get_vertex));
 Add<pmeshL*>("(","",new OneOperator2_<GlgVertex<MeshL>,pmeshL*,long>(get_vertex));
 Add<pmeshL* >("(", "", new OneQuadOperator< Op3_MeshDmp<MeshL>, Op3_MeshDmp<MeshL>::Op >);
    
 Add<pmeshL*>("be",".",new OneOperator1_<GlgBoundaryElement<MeshL>::BE,pmeshL*>(Build));
 Add<GlgElement<MeshL> >("adj",".",new OneOperator1_<GlgElement<MeshL>::Adj,GlgElement<MeshL> >(Build));
 Add<GlgBoundaryElement<MeshL>::BE>("(","",new OneOperator2_<GlgBoundaryElement<MeshL>,GlgBoundaryElement<MeshL>::BE,long>(get_element));
 Add<GlgElement<MeshL>::Adj>("(","",new OneOperator2_<GlgElement<MeshL>,GlgElement<MeshL>::Adj,long*>(get_adj));
 
 TheOperators->Add("==", new OneBinaryOperator<Op2_eq<GlgElement<MeshL>,GlgElement<MeshL> > >);
 TheOperators->Add("!=", new OneBinaryOperator<Op2_ne<GlgElement<MeshL>,GlgElement<MeshL> > >);
 TheOperators->Add("<", new OneBinaryOperator<Op2_lt<GlgElement<MeshL>,GlgElement<MeshL> > >);
 TheOperators->Add("<=", new OneBinaryOperator<Op2_le<GlgElement<MeshL>,GlgElement<MeshL> > >);
    
 Add<GlgBoundaryElement<MeshL> >("label",".",new OneOperator1_<long,GlgBoundaryElement<MeshL> >(getlab));
 Add<GlgBoundaryElement<MeshL> >("Element",".",new OneOperator1_<GlgElement<MeshL> ,GlgBoundaryElement<MeshL> >(getElement));
 Add<GlgBoundaryElement<MeshL> >("whoinElement",".",new OneOperator1_<long,GlgBoundaryElement<MeshL> >(NuElement));
 Add<GlgElement<MeshL> >("[","",new OneOperator2_<GlgVertex<MeshL> ,GlgElement<MeshL> ,long>(get_element));
 Add<GlgBoundaryElement<MeshL> >("[","",new OneOperator2_<GlgVertex<MeshL>,GlgBoundaryElement<MeshL>,long>(get_element));
 
 Add<GlgVertex<MeshL> >("x",".",new OneOperator1_<R,GlgVertex<MeshL> >(getx));
 Add<GlgVertex<MeshL> >("y",".",new OneOperator1_<R,GlgVertex<MeshL> >(gety));
 Add<GlgVertex<MeshL> >("z",".",new OneOperator1_<R,GlgVertex<MeshL> >(getz));
 Add<GlgVertex<MeshL> >("label",".",new OneOperator1_<long,GlgVertex<MeshL> >(getlab));
 Add<GlgElement<MeshL> >("label",".",new OneOperator1_<long,GlgElement<MeshL> >(getlab));
 Add<GlgElement<MeshL> >("region",".",new OneOperator1_<long,GlgElement<MeshL> >(getlab));
 Add<GlgElement<MeshL> >("mesure",".",new OneOperator1_<double,GlgElement<MeshL> >(getmes));
 Add<GlgElement<MeshL> >("measure",".",new OneOperator1_<double,GlgElement<MeshL> >(getmes));
 
 Add<pmeshL*>("mesure",".",new OneOperator1<double,pmeshL*>(pmesh_mes));
 Add<pmeshL*>("measure",".",new OneOperator1<double,pmeshL*>(pmesh_mes));
 Add<pmeshL*>("bordermesure",".",new OneOperator1<double,pmeshL*>(pmesh_mesb));
 Add<pmeshL*>("bordermeasure",".",new OneOperator1<double,pmeshL*>(pmesh_mesb));
 Add<pmeshL*>("nt",".",new OneOperator1<long,pmeshL*>(pmesh_nt));
 Add<pmeshL*>("nv",".",new OneOperator1<long,pmeshL*>(pmesh_nv));
 Add<pmeshL*>("nbe",".",new OneOperator1<long,pmeshL*>(pmesh_nbe));
 Add<pmeshL*>("hmax",".",new OneOperator1<double,pmeshL*>(pmesh_hmax));
 Add<pmeshL*>("hmin",".",new OneOperator1<double,pmeshL*>(pmesh_hmin));
 //Add<pmeshL*>("nbnomanifold",".",new OneOperator1<long,pmeshL*>(pmesh_nadjnomanifold));

    
    
    
    
    
    
 // 3D volume
 TheOperators->Add("<-",
       new OneOperator2_<pf3rbase*,pf3rbase*,pfes3* >(MakePtrFE3_2),
       new OneOperator3_<pf3rbasearray*,pf3rbasearray*,pfes3*,long >(MakePtrFE3_3),  
       new OneOperator2_<pf3cbase*,pf3cbase*,pfes3* >(MakePtrFE3_2),
       new OneOperator3_<pf3cbasearray*,pf3cbasearray*,pfes3*,long >(MakePtrFE3_3) //,
       //  new OneOperator2_<pmesharray*,pmesharray*,long >(MakePtr)
       );
 TheOperators->Add("<-",
	   new OneOperatorMakePtrFE3<double,v_fes3>(atype<double>()),  //  scalar case
	   new OneOperatorMakePtrFE3<double,v_fes3>(atype<E_Array>()),  //  vect case
	   new OneOperatorMakePtrFE3<Complex,v_fes3>(atype<Complex>()),  //  scalar complex  case
	   new OneOperatorMakePtrFE3<Complex,v_fes3>(atype<E_Array>())  //  vect complex case
       );
 TheOperators->Add("<-",
       new OneOperator2_<pfes3*,pfes3*,pfes3>(&set_copy_incr));
 TheOperators->Add("=",
	   new OneOperator2<pfes3*,pfes3*,pfes3>(&set_eqdestroy_incr)
	   );
 TheOperators->Add("=",
       new OneOperator2_<pf3r,pf3r,double,E_F_StackF0F0opt2<double> >(set_fe3<double,v_fes3>) ,      // modif template
       new OneOperator2_<pf3c,pf3c,Complex,E_F_StackF0F0opt2<Complex> >(set_fe3<Complex,v_fes3>)     // modif template
       ) ;     

// 3D surface
TheOperators->Add("<-",
      new OneOperator2_<pfSrbase*,pfSrbase*,pfesS* >(MakePtrFE3_2),
      new OneOperator3_<pfSrbasearray*,pfSrbasearray*,pfesS*,long >(MakePtrFE3_3),
      new OneOperator2_<pfScbase*,pfScbase*,pfesS* >(MakePtrFE3_2),
      new OneOperator3_<pfScbasearray*,pfScbasearray*,pfesS*,long >(MakePtrFE3_3) //,
      //  new OneOperator2_<pmesharray*,pmesharray*,long >(MakePtr)
      );
TheOperators->Add("<-",
      new OneOperatorMakePtrFE3<double,v_fesS>(atype<double>()),  //  scalar case
      new OneOperatorMakePtrFE3<double,v_fesS>(atype<E_Array>()),  //  vect case
      new OneOperatorMakePtrFE3<Complex,v_fesS>(atype<Complex>()),  //  scalar complex  case
      new OneOperatorMakePtrFE3<Complex,v_fesS>(atype<E_Array>())  //  vect complex case
      );
TheOperators->Add("<-",
      new OneOperator2_<pfesS*,pfesS*,pfesS>(&set_copy_incr));
TheOperators->Add("=",
      new OneOperator2<pfesS*,pfesS*,pfesS>(&set_eqdestroy_incr)
      );
TheOperators->Add("=",
      new OneOperator2_<pfSr,pfSr,double,E_F_StackF0F0opt2<double> >(set_fe3<double,v_fesS>) ,   // modif/ use template
      new OneOperator2_<pfSc,pfSc,Complex,E_F_StackF0F0opt2<Complex> >(set_fe3<Complex,v_fesS>) // modif/ use template
     ) ;

 // 3D curve
 TheOperators->Add("<-",
                   new OneOperator2_<pfLrbase*,pfLrbase*,pfesL* >(MakePtrFE3_2),
                   new OneOperator3_<pfLrbasearray*,pfLrbasearray*,pfesL*,long >(MakePtrFE3_3),
                   new OneOperator2_<pfLcbase*,pfLcbase*,pfesL* >(MakePtrFE3_2),
                   new OneOperator3_<pfLcbasearray*,pfLcbasearray*,pfesL*,long >(MakePtrFE3_3) //,
                   //  new OneOperator2_<pmesharray*,pmesharray*,long >(MakePtr)
                   );
 TheOperators->Add("<-",
                   new OneOperatorMakePtrFE3<double,v_fesL>(atype<double>()),  //  scalar case
                   new OneOperatorMakePtrFE3<double,v_fesL>(atype<E_Array>()),  //  vect case
                   new OneOperatorMakePtrFE3<Complex,v_fesL>(atype<Complex>()),  //  scalar complex  case
                   new OneOperatorMakePtrFE3<Complex,v_fesL>(atype<E_Array>())  //  vect complex case
                      );
 TheOperators->Add("<-",
                   new OneOperator2_<pfesL*,pfesL*,pfesL>(&set_copy_incr));
 TheOperators->Add("=",
                   new OneOperator2<pfesL*,pfesL*,pfesL>(&set_eqdestroy_incr)
                   );
 TheOperators->Add("=",
                   new OneOperator2_<pfLr,pfLr,double,E_F_StackF0F0opt2<double> >(set_fe3<double,v_fesL>) ,   // modif/ use template
                    new OneOperator2_<pfLc,pfLc,Complex,E_F_StackF0F0opt2<Complex> >(set_fe3<Complex,v_fesL>) // modif/ use template
                   ) ;
  // to write u(x,y,z) when u si FE function on Surface or Line .. FH. jan 2020. 
  /*  Add< pfSr >("(", "", new OneQuadOperator< Op4_pfeK< R, v_fesS>, Op4_pfeK< R,v_fesS >::Op >);
    Add< pfSc >("(", "", new OneQuadOperator< Op4_pfeK< Complex,v_fesS >, Op4_pfeK< Complex,v_fesS >::Op >);
    Add< pfLr >("(", "", new OneQuadOperator< Op4_pfeK< R,v_fesL >, Op4_pfeK< R,v_fesL >::Op >);
    Add< pfLc >("(", "", new OneQuadOperator< Op4_pfeK< Complex, v_fesL>, Op4_pfeK< Complex,v_fesL >::Op >);
*/
    
 map_type[typeid(double).name()]->AddCast(
   new E_F1_funcT<double,pf3r>(pf3r2R<R,0,v_fes3>)
   );
   
 map_type[typeid(Complex).name()]->AddCast(
   new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,0,v_fes3>)
 );

 map_type[typeid(double).name()]->AddCast(
   new E_F1_funcT<double,pfSr>(pfSr2R<R,0,v_fesS>)
   );

 map_type[typeid(Complex).name()]->AddCast(
   new E_F1_funcT<Complex,pfSc>(pfSr2R<Complex,0,v_fesS>)
 );
    
 map_type[typeid(double).name()]->AddCast(
   new E_F1_funcT<double,pfLr>(pfLr2R<R,0,v_fesL>)
 );
    
 map_type[typeid(Complex).name()]->AddCast(
   new E_F1_funcT<Complex,pfLc>(pfLr2R<Complex,0,v_fesL>)
 );

 Global.Add("dz","(",new OneOperatorCode<CODE_Diff<Ftest,op_dz> >);
 Global.Add("dxz","(",new OneOperatorCode<CODE_Diff<Ftest,op_dxz> >);
 Global.Add("dyz","(",new OneOperatorCode<CODE_Diff<Ftest,op_dyz> >);
 Global.Add("dzx","(",new OneOperatorCode<CODE_Diff<Ftest,op_dzx> >);
 Global.Add("dzx","(",new OneOperatorCode<CODE_Diff<Ftest,op_dzy> >);
 Global.Add("dzz","(",new OneOperatorCode<CODE_Diff<Ftest,op_dzz> >);


 Global.Add("dz","(",new OneOperatorCode<CODE_Diff<Finconnue,op_dz> >);
 Global.Add("dxz","(",new OneOperatorCode<CODE_Diff<Finconnue,op_dxz> >);
 Global.Add("dyz","(",new OneOperatorCode<CODE_Diff<Finconnue,op_dyz> >);
 Global.Add("dzx","(",new OneOperatorCode<CODE_Diff<Finconnue,op_dzx> >);
 Global.Add("dzx","(",new OneOperatorCode<CODE_Diff<Finconnue,op_dzy> >);
 Global.Add("dzz","(",new OneOperatorCode<CODE_Diff<Finconnue,op_dzz> >);


 
// bof
 // volume 3D real
 Global.Add("dx","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dx,v_fes3>));
 Global.Add("dy","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dy,v_fes3>));
 Global.Add("dz","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dz,v_fes3>));
 Global.Add("dxx","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dxx,v_fes3>));
 Global.Add("dyy","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dyy,v_fes3>));
 Global.Add("dxy","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dxy,v_fes3>));
 Global.Add("dyx","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dyx,v_fes3>));
 Global.Add("dzx","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dzx,v_fes3>));
 Global.Add("dzy","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dzy,v_fes3>));
 Global.Add("dzz","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dzz,v_fes3>));
 Global.Add("dxz","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dxz,v_fes3>));
 Global.Add("dyz","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dyz,v_fes3>));
 // volume 3D complex
 Global.Add("dx","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dx,v_fes3>));
 Global.Add("dy","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dy,v_fes3>));
 Global.Add("dz","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dz,v_fes3>));
 Global.Add("dxx","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dxx,v_fes3>));
 Global.Add("dyy","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dyy,v_fes3>));
 Global.Add("dxy","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dxy,v_fes3>));
 Global.Add("dyx","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dyx,v_fes3>));
 Global.Add("dzx","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dzx,v_fes3>));
 Global.Add("dzy","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dzy,v_fes3>));
 Global.Add("dzz","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dzz,v_fes3>));
 Global.Add("dxz","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dxz,v_fes3>));
 Global.Add("dyz","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dyz,v_fes3>));
 // surface 3D real
 Global.Add("dx","(",new E_F1_funcT<double,pfSr>(pfSr2R<R,op_dx,v_fesS>));
 Global.Add("dy","(",new E_F1_funcT<double,pfSr>(pfSr2R<R,op_dy,v_fesS>));
 Global.Add("dz","(",new E_F1_funcT<double,pfSr>(pfSr2R<R,op_dz,v_fesS>));
 Global.Add("dxx","(",new E_F1_funcT<double,pfSr>(pfSr2R<R,op_dxx,v_fesS>));
 Global.Add("dyy","(",new E_F1_funcT<double,pfSr>(pfSr2R<R,op_dyy,v_fesS>));
 Global.Add("dxy","(",new E_F1_funcT<double,pfSr>(pfSr2R<R,op_dxy,v_fesS>));
 Global.Add("dyx","(",new E_F1_funcT<double,pfSr>(pfSr2R<R,op_dyx,v_fesS>));
 Global.Add("dzx","(",new E_F1_funcT<double,pfSr>(pfSr2R<R,op_dzx,v_fesS>));
 Global.Add("dzy","(",new E_F1_funcT<double,pfSr>(pfSr2R<R,op_dzy,v_fesS>));
 Global.Add("dzz","(",new E_F1_funcT<double,pfSr>(pfSr2R<R,op_dzz,v_fesS>));
 Global.Add("dxz","(",new E_F1_funcT<double,pfSr>(pfSr2R<R,op_dxz,v_fesS>));
 Global.Add("dyz","(",new E_F1_funcT<double,pfSr>(pfSr2R<R,op_dyz,v_fesS>));
 // surface 3D complex
 Global.Add("dx","(",new E_F1_funcT<Complex,pfSc>(pfSr2R<Complex,op_dx,v_fesS>));
 Global.Add("dy","(",new E_F1_funcT<Complex,pfSc>(pfSr2R<Complex,op_dy,v_fesS>));
 Global.Add("dz","(",new E_F1_funcT<Complex,pfSc>(pfSr2R<Complex,op_dz,v_fesS>));
 Global.Add("dxx","(",new E_F1_funcT<Complex,pfSc>(pfSr2R<Complex,op_dxx,v_fesS>));
 Global.Add("dyy","(",new E_F1_funcT<Complex,pfSc>(pfSr2R<Complex,op_dyy,v_fesS>));
 Global.Add("dxy","(",new E_F1_funcT<Complex,pfSc>(pfSr2R<Complex,op_dxy,v_fesS>));
 Global.Add("dyx","(",new E_F1_funcT<Complex,pfSc>(pfSr2R<Complex,op_dyx,v_fesS>));
 Global.Add("dzx","(",new E_F1_funcT<Complex,pfSc>(pfSr2R<Complex,op_dzx,v_fesS>));
 Global.Add("dzy","(",new E_F1_funcT<Complex,pfSc>(pfSr2R<Complex,op_dzy,v_fesS>));
 Global.Add("dzz","(",new E_F1_funcT<Complex,pfSc>(pfSr2R<Complex,op_dzz,v_fesS>));
 Global.Add("dxz","(",new E_F1_funcT<Complex,pfSc>(pfSr2R<Complex,op_dxz,v_fesS>));
 Global.Add("dyz","(",new E_F1_funcT<Complex,pfSc>(pfSr2R<Complex,op_dyz,v_fesS>));
 // curve 3D real
 Global.Add("dx","(",new E_F1_funcT<double,pfLr>(pfLr2R<R,op_dx,v_fesL>));
 Global.Add("dy","(",new E_F1_funcT<double,pfLr>(pfLr2R<R,op_dy,v_fesL>));
 Global.Add("dz","(",new E_F1_funcT<double,pfLr>(pfLr2R<R,op_dz,v_fesL>));
 Global.Add("dxx","(",new E_F1_funcT<double,pfLr>(pfLr2R<R,op_dxx,v_fesL>));
 Global.Add("dyy","(",new E_F1_funcT<double,pfLr>(pfLr2R<R,op_dyy,v_fesL>));
 Global.Add("dxy","(",new E_F1_funcT<double,pfLr>(pfLr2R<R,op_dxy,v_fesL>));
 Global.Add("dyx","(",new E_F1_funcT<double,pfLr>(pfLr2R<R,op_dyx,v_fesL>));
 Global.Add("dzx","(",new E_F1_funcT<double,pfLr>(pfLr2R<R,op_dzx,v_fesL>));
 Global.Add("dzy","(",new E_F1_funcT<double,pfLr>(pfLr2R<R,op_dzy,v_fesL>));
 Global.Add("dzz","(",new E_F1_funcT<double,pfLr>(pfLr2R<R,op_dzz,v_fesL>));
 Global.Add("dxz","(",new E_F1_funcT<double,pfLr>(pfLr2R<R,op_dxz,v_fesL>));
 Global.Add("dyz","(",new E_F1_funcT<double,pfLr>(pfLr2R<R,op_dyz,v_fesL>));
 // curve 3D complex
 Global.Add("dx","(",new E_F1_funcT<Complex,pfLc>(pfLr2R<Complex,op_dx,v_fesL>));
 Global.Add("dy","(",new E_F1_funcT<Complex,pfLc>(pfLr2R<Complex,op_dy,v_fesL>));
 Global.Add("dz","(",new E_F1_funcT<Complex,pfLc>(pfLr2R<Complex,op_dz,v_fesL>));
 Global.Add("dxx","(",new E_F1_funcT<Complex,pfLc>(pfLr2R<Complex,op_dxx,v_fesL>));
 Global.Add("dyy","(",new E_F1_funcT<Complex,pfLc>(pfLr2R<Complex,op_dyy,v_fesL>));
 Global.Add("dxy","(",new E_F1_funcT<Complex,pfLc>(pfLr2R<Complex,op_dxy,v_fesL>));
 Global.Add("dyx","(",new E_F1_funcT<Complex,pfLc>(pfLr2R<Complex,op_dyx,v_fesL>));
 Global.Add("dzx","(",new E_F1_funcT<Complex,pfLc>(pfLr2R<Complex,op_dzx,v_fesL>));
 Global.Add("dzy","(",new E_F1_funcT<Complex,pfLc>(pfLr2R<Complex,op_dzy,v_fesL>));
 Global.Add("dzz","(",new E_F1_funcT<Complex,pfLc>(pfLr2R<Complex,op_dzz,v_fesL>));
 Global.Add("dxz","(",new E_F1_funcT<Complex,pfLc>(pfLr2R<Complex,op_dxz,v_fesL>));
 Global.Add("dyz","(",new E_F1_funcT<Complex,pfLc>(pfLr2R<Complex,op_dyz,v_fesL>));
 
    
    
    
    
    
    
    
    
    
    
 // 3d volume
 Global.Add("int3d","(",new OneOperatorCode<CDomainOfIntegration3d>);
 Global.Add("int2d","(",new OneOperatorCode<CDomainOfIntegrationBorder3d>);
 Global.Add("intallfaces","(",new OneOperatorCode<CDomainOfIntegrationAllFaces>);

 // 3d surface
 Global.Add("int2d","(",new OneOperatorCode<CDomainOfIntegrationS>);
 Global.Add("int1d","(",new OneOperatorCode<CDomainOfIntegrationBorderS>);
 Global.Add("intalledges","(",new OneOperatorCode<CDomainOfIntegrationAllEdgesS>);
    
 // 3d curve
 Global.Add("int1d","(",new OneOperatorCode<CDomainOfIntegrationL>);
 Global.Add("int0d","(",new OneOperatorCode<CDomainOfIntegrationBorderL>); // point int0d??
 /*decommente par J. Morice 14/01/09*/

 Add<pf3r>("refresh",".",new OneOperator1<bool,pf3r>(pfer_refresh3<R,v_fes3>));
 Add<pf3c>("refresh",".",new OneOperator1<bool,pf3c>(pfer_refresh3<Complex,v_fes3>));

 // 3d volume FE
 Add<pf3r>("n",".",new OneOperator1<long,pf3r>(pf3r_nbdf<R>));
 Add<pf3c>("n",".",new OneOperator1<long,pf3c>(pf3r_nbdf<Complex>));
    
 Add<pf3r>("Th",".",new OneOperator1<pmesh3 ,pf3r>(pf3r_Th<R>));
 Add<pf3c>("Th",".",new OneOperator1<pmesh3,pf3c>(pf3r_Th<Complex>));
    
 Add<pfes3*>("ndof",".",new OneOperator1<long,pfes3*>(pVh3_ndof));
 Add<pfes3*>("nt",".",new OneOperator1<long,pfes3*>(pVh3_nt));
 Add<pfes3*>("ndofK",".",new OneOperator1<long,pfes3*>(pVh3_ndofK));
 Add<pfes3*>("Th",".",new OneOperator1<pmesh3,pfes3*>(pVh3_Th));
 
 // 3d surface FE
 Add<pfSr>("n",".",new OneOperator1<long,pfSr>(pfSr_nbdf<R>));
 Add<pfSc>("n",".",new OneOperator1<long,pfSc>(pfSr_nbdf<Complex>));
 Add<pfSr>("Th",".",new OneOperator1<pmeshS ,pfSr>(pfSr_Th<R>));
 Add<pfSc>("Th",".",new OneOperator1<pmeshS,pfSc>(pfSr_Th<Complex>));
 Add<pfesS*>("ndof",".",new OneOperator1<long,pfesS*>(pVhS_ndof));
 Add<pfesS*>("nt",".",new OneOperator1<long,pfesS*>(pVhS_nt));
 Add<pfesS*>("ndofK",".",new OneOperator1<long,pfesS*>(pVhS_ndofK));
// Add<pfesS*>("Th",".",new OneOperator1<pmesh3,pfesS*>(pVhS_Th));
 
// 3d curve FE
Add<pfLr>("n",".",new OneOperator1<long,pfLr>(pfLr_nbdf<R>));
Add<pfLc>("n",".",new OneOperator1<long,pfLc>(pfLr_nbdf<Complex>));
Add<pfLr>("Th",".",new OneOperator1<pmeshL ,pfLr>(pfLr_Th<R>));
Add<pfLc>("Th",".",new OneOperator1<pmeshL,pfLc>(pfLr_Th<Complex>));
Add<pfesL*>("ndof",".",new OneOperator1<long,pfesL*>(pVhL_ndof));
Add<pfesL*>("nt",".",new OneOperator1<long,pfesL*>(pVhL_nt));
Add<pfesL*>("ndofK",".",new OneOperator1<long,pfesL*>(pVhL_ndofK));
    
    
 //Add<pf3rbasearray*>("[","",new OneOperator2_<pf3rbase*,pf3rbasearray*,long>(get_element));
 //Add<pf3rarray>("[","",new OneOperator2_<pf3r,pf3rarray,long>(get_element));
 //Add<pf3carray>("[","",new OneOperator2_<pf3c,pf3carray,long>(get_element));
 
 // 3d volume
 Add<pf3cbasearray*>("[","",new OneOperator2_<pf3cbase*,pf3cbasearray*,long>(get_element));  // use ???? FH sep. 2009
 Add<pf3rbasearray*>("[","",new OneOperator2_<pf3rbase*,pf3rbasearray*,long>(get_element));  //  use ???? FH sep. 2009
 Add<pf3rarray>("[","",new OneOperator2_FE_get_elmnt<double,v_fes3>());// new version FH sep 2009
 Add<pf3carray>("[","",new OneOperator2_FE_get_elmnt<Complex,v_fes3>());
 Add<pfes3*>("(","", new OneTernaryOperator<pVh3_ndf,pVh3_ndf::Op>  );
 // 3d surface
 Add<pfScbasearray*>("[","",new OneOperator2_<pfScbase*,pfScbasearray*,long>(get_element));  // use ???? FH sep. 2009
 Add<pfSrbasearray*>("[","",new OneOperator2_<pfSrbase*,pfSrbasearray*,long>(get_element));  //  use ???? FH sep. 2009
 Add<pfSrarray>("[","",new OneOperator2_FE_get_elmnt<double,v_fesS>());// new version FH sep 2009
 Add<pfScarray>("[","",new OneOperator2_FE_get_elmnt<Complex,v_fesS>());
 Add<pfesS*>("(","", new OneTernaryOperator<pVhS_ndf,pVhS_ndf::Op>  );
 // 3d surface
 Add<pfLcbasearray*>("[","",new OneOperator2_<pfLcbase*,pfLcbasearray*,long>(get_element));  // use ???? FH sep. 2009
 Add<pfLrbasearray*>("[","",new OneOperator2_<pfLrbase*,pfLrbasearray*,long>(get_element));  //  use ???? FH sep. 2009
 Add<pfLrarray>("[","",new OneOperator2_FE_get_elmnt<double,v_fesL>());// new version FH sep 2009
 Add<pfLcarray>("[","",new OneOperator2_FE_get_elmnt<Complex,v_fesL>());
 Add<pfesL*>("(","", new OneTernaryOperator<pVhL_ndf,pVhL_ndf::Op>  );
    
 
 init_mesh3_array(); //3D vomlume
 init_meshS_array();  //3D surface
 init_meshL_array();  //3D line

 // Add jan 2019 F.H ..to get a sorted the array of label and region of a mesh.
 Global.Add("labels","(",new OneOperator1s_<KN_<long>,pmeshL>(listoflabel));
 Global.Add("labels","(",new OneOperator1s_<KN_<long>,pmeshS>(listoflabel));
 Global.Add("labels","(",new OneOperator1s_<KN_<long>,pmesh3>(listoflabel));
 Global.Add("labels","(",new OneOperator1s_<KN_<long>,pmesh>(listoflabel));
 Global.Add("regions","(",new OneOperator1s_<KN_<long>,pmeshL>(listofregion));
 Global.Add("regions","(",new OneOperator1s_<KN_<long>,pmeshS>(listofregion));
 Global.Add("regions","(",new OneOperator1s_<KN_<long>,pmesh3>(listofregion));
 Global.Add("regions","(",new OneOperator1s_<KN_<long>,pmesh>(listofregion));
  
    
}

//#include "InitFunct.hpp"
//static addingInitFunct TheaddingInitFunct(-10,init_lgmesh3);
template E_set_fev3<double,v_fes3>::E_set_fev3(const E_Array * a,Expression pp) ;
template E_set_fev3<Complex,v_fes3>::E_set_fev3(const E_Array * a,Expression pp) ;
template E_set_fev3<double,v_fesS>::E_set_fev3(const E_Array * a,Expression pp) ;
template E_set_fev3<Complex,v_fesS>::E_set_fev3(const E_Array * a,Expression pp) ;
template E_set_fev3<double,v_fesL>::E_set_fev3(const E_Array * a,Expression pp) ;
template E_set_fev3<Complex,v_fesL>::E_set_fev3(const E_Array * a,Expression pp) ;
