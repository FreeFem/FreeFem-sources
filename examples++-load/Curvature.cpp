// ORIG-DATE: Novembre 2016
//  Calcul de la coubure d'un courbe discrete
// Compilation: ff-c++ -auto courbure.cpp

// re parametrisation
// -*- Mode : c++ -*-
//
// SUMMARY  :
// USAGE    : LGPL
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE
// AUTHOR   :  F. Hecht Jan. 2016
// E-MAIL   : frederic.hecht@upmc.fr
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

/*
 
 Usage dans freefem++:
 mesh Th;
 int[int] labs=[1,2,3];
 int lab =1;
 fespace Vh(Th,P1);
 Vh courbure;
 courbure[]= courbure(Th,labs); // courbure des bord de labs de Th
 courbure[]= courbure(Th,lab); // courbure du bord de lab de Th
 int np = 100;
 real[int,int] SL(3,np);
 // .. set SL here
// SL(0,.) === x
// SL(1,.) === y
// SL(3,.) === absisse curviligne
 real lg=reparametrage(SL); // calcul de absisse curviligne de SL 
 
 // Border defini avec SL..
 border bSL(t=0,1)
 {
 P= Curve(SL,t);
 label=2;
 }
 // rediscretisation par equirepation de absisse curviligne de SL  en np point
 SL=equiparametre(SL,np);
 
 real a11=3,a12=2,a22=1;
 cout << vp1(a11,a12,a22); // premier vp de [[a11,a12][a12,a22]]
 */



#ifndef WITH_NO_INIT
#include "ff++.hpp"
#include "AFunction_ext.hpp"
#include "eigenv.h"

#endif


using namespace std;

using namespace  Fem2D;



// fonction determinant les points d'intersection
static int debug =0;


R3  * courbe(Stack stack,const KNM_<double> &b,const  long &li0,const  long & li1,const double & ss, long *const&   pi)
{
    assert(b.N() >=3);
    int i0=li0,i1=li1,im;
    if(i0<0) i0=0;
    if(i1<0) i1=b.M()-1;
    double lg=b(2,i1);
    R3 Q;
    ffassert(lg>0 && b(2,0)==0.);
    double s = ss*lg;
    int k=0,k1=i1;
    while(i0 < i1-1)
    {
        ffassert(k++ < k1);
        im = (i0+i1)/2;
        if(s <b(2,im)  )
        { i1=im;
        }
        else if(s>b(2,im)  )
        { i0= im;
        }
        else {  Q=R3(b(0,im),b(1,im),0);  i0=i1=im;break;}
    }
    if(i0<i1)
    {
        ffassert(b(2,i0) <= s );
        ffassert(b(2,i1) >= s );
        R2 A(b(0,i0),b(1,i0));
        R2 B(b(0,i1),b(1,i1));
        double l1=(b(2,i1)-s);
        double l0=s-b(2,i0);
        Q= (l1*A + l0*B)/(l1+l0);
    }
    if(pi) *pi=i0;
    R3 *pQ = Add2StackOfPtr2Free(stack,new R3(Q));
    // MeshPoint &mp= *MeshPointStack(stack); // the struct to get x,y, normal , value
    //mp.P.x=Q.x; // get the current x value
    //mp.P.y=Q.y; // get the current y value
    return pQ;
}
R3   * courbe(Stack stack,const KNM_<double> &b,const double & ss)
{
    return courbe(stack,b,-1,-1,ss,0);
}
/*
   Remarque sur la courbure 2d
   a = angle tangent / ox
   R rayon de coubure 
   s abs curviline
   c courbure = 1/R
   ds = R d a
   c = da/ds 
   => dans le cas d'une ligne brise
   c = [a] \delta (dirac) en ds 
   
   Courbure 3d  c = 1/R1 + 1/R2 = cr+ ca
 avec :
  r == x
  z == y
  donc en axi  
   cr = courbure 2d classsique
   soit N=(Nr,Nz) normal au point (r,z) unitaire
 
   alors Ra = r/Nr; car la composante r de   (r/Nr  N)  =  r; 
   donc ca = Nr/r ;
 
 
 
 */
KN<double>* courbure(Stack stack,pmesh const & pTh,KN<long> *const & lab,bool axi)
{
    const double pi= M_PI;
    const double twopi= 2*M_PI;
    const Mesh &Th=*pTh;
    KN<double> *pc=new KN<double>(Th.nv);
    KN<double> & c=*pc;
    c=0.;
    map<long,int> mlab;
    for (int i=0; i< lab->N(); ++i)
    {
        long l=(*lab)[i];
        if(verbosity>9)
        cout << i << " lab "<< l << endl;
        mlab[l]=i;
    }
    KN<int> cn(Th.nv);
    KN<double> le(Th.nv);
    KN<double> nr(Th.nv);
    cn=0;
    le=0.;
    nr=0.;
    for(int ee=0; ee< Th.neb; ++ee)
    {
        BoundaryEdge e = Th.be(ee);
        if( mlab.find(e.lab) != mlab.end() )
        {
            int ie,i =Th.BoundaryElement(ee,ie);
            const Triangle & K(Th[i]);
            R2 E=K.Edge(ie);
            double lE=E.norme();
            int iv[]= {Th(K[VerticesOfTriangularEdge[ie][0]]),
                Th(K[VerticesOfTriangularEdge[ie][1]])};
            for (int j=0; j<2;++j)
            {
                int k = iv[j];
                double aa = atan2(E.y,E.x);
                if(j==1) aa = -aa;
                c[k] += aa;
                cn[k]+= 1+(j==1);
                le[k]+=lE;
                if( axi) nr[k]+= E.y/lE;
                 // cout << k << " " << lE << " "<< le[k] << endl;
            }
        }
    }
    double epsTh = sqrt(Th.area/Th.nt)*1e-6;// bofbof
    for(int i=0;i<Th.nv;++i)
        if( cn[i]==3)
        {
            
            //  cout << i << " courbure " << c[i] << " " << le[i] << " " << c[i]/le[i]/2. <<endl;
            if( c[i]> pi) c[i]-=twopi;
            else if (c[i]< -pi ) c[i]+=twopi;
            
            if(axi)
            {
                // Bof Bof  moyen des normal:
                //  c[i] angle entre n1 et n2 :
                //   ( n1 + n2 )/ || n1+n2|| et on a :|| n1+n2|| = 1+cos(abs(c[i]))
                double Nr = nr[i]/(1.+cos(abs(c[i])));// ok
                double r=Th(i).x;
                if(verbosity>9999)
                cout << " R1 " << r/Nr << " R2 " << 1/(c[i]/(le[i]/2.)) << " da=" << c[i] << " le " << le[i]/2 << endl;

                c[i] = r*c[i]/(le[i]/2.)+ Nr ;
            }
            else c[i] /= (le[i]/2.);
        }
        else if ((cn[i]>0) && (abs(Th(i).x) < epsTh && axi))
        {
            // extermite r=0
           // cout << "axi 0" << nr[i] <<  endl;
            double ci=c[i];
            if(cn[i]==2) ci=pi+ci;
            ci *=2;
            double rm=le[i]*sqrt(1-nr[i]*nr[i])/2;
            if(verbosity>999)
              cout <<  Th(i).y << " R1 "<< ci/le[i] << " ci " << ci << " le " << le[i] << " cn " << cn[i] << " " << c[i] << " " << nr[i] << " " << asin(nr[i]) <<endl;
            c[i] = ci*2*rm;//2*nr[i];
        }
        else
            c[i]=0;
    //cout << c << endl;
    return Add2StackOfPtr2FreeRC(stack,pc);
}


KN<double>* courbure(Stack stack,pmesh const & pTh,KN<long> *const & lab)
{
    return courbure(stack,pTh,lab,0);
}
KN<double>* courbure(Stack stack,pmesh const & pTh,const long & lab)
{
    KN<long> ll(1); ll=lab;
    return courbure(stack,pTh,&ll,0);
}
KN<double>* courbureaxi(Stack stack,pmesh const & pTh,KN<long> *const & lab)
{
    return courbure(stack,pTh,lab,1);
}
KN<double>* courbureaxi(Stack stack,pmesh const & pTh,const long & lab)
{
    KN<long> ll(1); ll=lab;
    return courbure(stack,pTh,&ll,1);
}

double  reparametrage(Stack stack,const KNM_<double> &bb)
{
    KNM_<double> b=bb;
    ffassert(b.N()>=3);
    R2 P(b(0,0),b(1,0));
    double s =0;
    b(2,0)=  s;
    for(int i=1;i<b.M();++i)
    {
        R2 Q(b(0,i),b(1,i));
        s += R2(P,Q).norme() ;
        b(2,i)=  s;
        P = Q;
    }
    return s;
}

KNM<double>  * equiparametre(Stack stack,const KNM_<double> &bb,const long & n)
{
    double lg=reparametrage(stack,bb);
    KNM_<double> b=bb;
    KNM<double> *pc =new KNM<double>(3,n);
    
    KNM<double> & c=*pc;
    int m =b.M();
    int n1=n-1,m1=m-1;
    double delta = 1./n1;
    ffassert(b.N()==3);
    R2 P(b(0,0),b(1,0));
    double s =0;
    c(':',0)=b(':',0);
    c(':',n1) = b(':',m1);
    for(int i=1;i<n1;++i)
    {
        double s=i*delta;
        R3 P= *courbe(stack,bb,s);
        c(0,i)=P.x;
        c(1,i)=P.y;
        c(2,i)=s*lg;
        if(debug) cout << i << " " << P << " " << s << endl;
    }
    return  Add2StackOfPtr2FreeRC(stack,pc);
}

double ExtractBorder(Stack stack,pmesh const & pTh,KN_<long> const & lab, KNM<double> *const &bb)
{
    const Mesh &Th=*pTh;
    map<long,int> mlab;
    for (int i=0; i< lab.N(); ++i)
    {
        long l=lab[i];
        if(verbosity>9)
            cout << i << " lab "<< l << endl;
        mlab[l]=i;
    }
    KN<long> nx(Th.nv),nee(Th.neb*2),mark(Th.nv);
    nx=-1;
    mark=-1;
    int nel=0;
    for(int ee=0,k=0; ee< Th.neb; ++ee)
    {
        BoundaryEdge e = Th.be(ee);
        if( mlab.find(e.lab) != mlab.end() )
        {
            int ie,i =Th.BoundaryElement(ee,ie);
            const Triangle & K(Th[i]);
            R2 E=K.Edge(ie);
            int iv[]= {Th(K[VerticesOfTriangularEdge[ie][0]]),
                Th(K[VerticesOfTriangularEdge[ie][1]])};
            nx[iv[0]]=nel;
            mark[iv[1]]=nel;
            nee[nel++]=iv[1];
            nee[nel++]=iv[0];
            if(verbosity>99) cout << " "<< nel/2 << " : "<<iv[1] << " " << iv[0] <<endl;
        }
    }
    if(verbosity>9)
        cout << " n edge  "<< nel/2   << endl;
    // recherech depart
    int bg =-1,nbg=0;
    if(nel==0) return 0;
    for( int k=0; k<nel;++k)
    {
        int j= nee[k];
        if( nx[j] >=0 && mark[j]<0)
        { nbg++;
            bg = j;
        }
        
        
    }
    int np = nel/2+1;
    if(nbg ==0) {// bord ferme
        bg= nee[0];
        nee[0]=-1; // on ouvre 
    }
    else
    { if(( nbg != 1)|| (verbosity > 4)) cout << " error (no connexe boundary be carefull with internal boundary (pb of sens) ) : nb start = " << nbg << endl;
      ffassert( nbg==1); // un depart pas plus
    }
    bb->resize(3,np);
    KNM<double> &b(*bb);
    int i=0,iv=bg;
    while (iv>=0 && i< np)
    {
        
        b(0,i)=Th(iv).x;
        b(1,i)=Th(iv).y;
        b(2,i)=0; // compute after
        i++;
	if( nx[iv] < 0) break;
        iv =  nee[nx[iv]];
       
    }
    if(nbg==0 ) {
        // on ferme
        ffassert(i+1==np);
        iv=bg;
        b(0,i)=Th(iv).x;
        b(1,i)=Th(iv).y;
        b(2,i)=0; // compute after
        i++;
    }
    ffassert(i==np); // bonne longueur => sinon bug plus d'une composant connexe ????
   return reparametrage(stack,b);
}
double ExtractBorder(Stack stack,pmesh const & pTh,long const & lab, KNM<double> *const &bb)
{
    KN<long> tab(1);
    tab=lab;
    return ExtractBorder(stack,pTh,tab,bb);
}


#define  EPSD           1.e-15
#define  EPSD2          1.e-10
#define  EPS6           5.e-06
#define  EPS            1.e-06
#define  EPSX2          2.e-06
#define  MAXTOU         50
/* check if numbers are equal */
#define egal(x,y)   ( \
(  ((x) == 0.0f) ? (fabs(y) < EPS) : \
( ((y) == 0.0f) ? (fabs(x) < EPS) : \
(fabs((x)-(y)) / (fabs(x) + fabs(y)) < EPSX2) )  ) )

/* eigen value + vector extraction */
int eigen2(double *mm,double *lambda,double vp[2][2]) {
    double   m[3],dd,a1,xn,ddeltb,rr1,rr2,ux,uy;
    
    /* init */
    ux = 1.0;
    uy = 0.0;
    
    /* normalize */
    memcpy(m,mm,3*sizeof(double));
    xn = fabs(m[0]);
    if ( fabs(m[1]) > xn )  xn = fabs(m[1]);
    if ( fabs(m[2]) > xn )  xn = fabs(m[2]);
    if ( xn < EPSD2 ) {
        lambda[0] = lambda[1] = 0.0;
        vp[0][0] = 1.0;
        vp[0][1] = 0.0;
        vp[1][0] = 0.0;
        vp[1][1] = 1.0;
        return(1);
    }
    xn = 1.0 / xn;
    m[0] *= xn;
    m[1] *= xn;
    m[2] *= xn;
    
    if ( egal(m[1],0.0) ) {
        rr1 = m[0];
        rr2 = m[2];
        goto vect;
    }
    
    /* eigenvalues of jacobian */
    a1	 = -(m[0] + m[2]);
    ddeltb = a1*a1 - 4.0 * (m[0]*m[2] - m[1]*m[1]);
    
    if ( ddeltb < 0.0 ) {
        fprintf(stderr,"  Delta: %f\n",ddeltb);
        ddeltb = 0.0;
    }
    ddeltb = sqrt(ddeltb);
    
    if ( fabs(a1) < EPS ) {
        rr1 = 0.5 * sqrt(ddeltb);
        rr2 = -rr1;
    }
    else if ( a1 < 0.0 ) {
        rr1 = 0.5 * (-a1 + ddeltb);
        rr2 = (-m[1]*m[1] + m[0]*m[2]) / rr1;
    }
    else if ( a1 > 0.0 ) {
        rr1 = 0.5 * (-a1 - ddeltb);
        rr2 = (-m[1]*m[1] + m[0]*m[2]) / rr1;
    }
    else {
        rr1 = 0.5 * ddeltb;
        rr2 = -rr1;
    }
    
vect:
    xn = 1.0 / xn;
    lambda[0] = rr1 * xn;
    lambda[1] = rr2 * xn;
    
    /* eigenvectors */
    a1 = m[0] - rr1;
    if ( fabs(a1)+fabs(m[1]) < EPS ) {
        if (fabs(lambda[1]) < fabs(lambda[0]) ) {
            ux = 1.0;
            uy = 0.0;
        }
        else {
            ux = 0.0;
            uy = 1.0;
        }
    }
    else if ( fabs(a1) < fabs(m[1]) ) {
        ux = 1.0;
        uy = -a1 / m[1];
    }
    else if ( fabs(a1) > fabs(m[1]) ) {
        ux = -m[1] / a1;
        uy = 1.0;
    }
    else if ( fabs(lambda[1]) > fabs(lambda[0]) ) {
        ux = 0.0;
        uy = 1.0;
    }
    else {
        ux = 1.0;
        uy = 0.0;
    }
    
    dd = sqrt(ux*ux + uy*uy);
    dd = 1.0 / dd;
    if ( fabs(lambda[0]) > fabs(lambda[1]) ) {
        vp[0][0] =  ux * dd;
        vp[0][1] =  uy * dd;
    }
    else {
        vp[0][0] =  uy * dd;
        vp[0][1] = -ux * dd;
    }
    
    /* orthogonal vector */
    vp[1][0] = -vp[0][1];
    vp[1][1] =  vp[0][0];
    
    return(1);
}


double vp1(const double &a11,const double &a12,const double &a22)
{
    double vp[2][2];
    double l[2];
    double m[3]={a11,a12,a22};
    int nv=eigen2(m,l,vp);
    return l[0];
}
double Tresca(const double &a11,const double &a12,const double &a22)
{
    double vp[2][2];
    double l[2];
    double m[3]={a11,a12,a22};
    int nv=eigen2(m,l,vp);
    return max(fabs(l[0]-l[1]),max(fabs(l[0]),fabs(l[1])));
}
double Tresca(const double &arr,const double &arz,const double &azz,const double &att)
{
    double vp[2][2];
    double l[3];
    double m[3]={arr,arz,azz};
    int nv=eigen2(m,l,vp);
    l[2]=att;
    
    return max(fabs(l[0]-l[1]),max(fabs(l[0]-l[2]),fabs(l[1]-l[2])));
}
double VonMises(const double &a11,const double &a12,const double &a22)
{
    double vp[2][2];
    double l[3];
    double m[3]={a11,a12,a22};
    int nv=eigen2(m,l,vp);
    l[2]=0;
    double s1 = l[1]-l[0];
    double s2 = l[0]-l[2];
    double s3 = l[1]-l[2];
    return sqrt((s1*s1+s2*s2+s3*s3)/2.);
}
double VonMises(const double &a11,const double &a12,const double &a22,const double &att)
{
    double vp[2][2];
    double l[3];
    double m[3]={a11,a12,a22};
    int nv=eigen2(m,l,vp);
    l[2]=att;
    double s1 = l[1]-l[0];
    double s2 = l[0]-l[1];
    double s3 = l[1]-l[2];
    return sqrt((s1*s1+s2*s2+s3*s3)/2.);
}

static void finit()
{
    
    Global.Add("extractborder","(",new OneOperator3s_<double,pmesh,KN_<long>,KNM<double>*>(ExtractBorder));
    Global.Add("extractborder","(",new OneOperator3s_<double,pmesh,long,KNM<double>*>(ExtractBorder));
    
    Global.Add("curvature","(",new OneOperator2s_<KN<double>*,pmesh ,KN<long> *>(courbure));
    Global.Add("curvature","(",new OneOperator2s_<KN<double>*,pmesh ,long>(courbure));
    Global.Add("raxicurvature","(",new OneOperator2s_<KN<double>*,pmesh ,KN<long> *>(courbureaxi));
    Global.Add("raxicurvature","(",new OneOperator2s_<KN<double>*,pmesh ,long>(courbureaxi));
    Global.Add("curves","(",new OneOperator2s_<R3*,KNM_<double>,double>(courbe));
    Global.Add("setcurveabcisse","(",new OneOperator1s_<double,KNM_<double> >(reparametrage));
    Global.Add("equiparameter","(",new OneOperator2s_<KNM<double>* ,KNM_<double>, long  >(equiparametre));
   // Global.Add("vp1","(",new OneOperator3_<double, double,double,double >(vp1));
    Global.Add("Tresca","(",new OneOperator3_<double, double,double,double >(Tresca));
    Global.Add("VonMises","(",new OneOperator3_<double, double,double,double >(VonMises));
    Global.Add("Tresca","(",new OneOperator4_<double, double,double,double,double >(Tresca));
    Global.Add("VonMises","(",new OneOperator4_<double, double,double,double,double >(VonMises));
    
    
}

LOADFUNC(finit);  //  une variable globale qui serat construite  au chargement dynamique
