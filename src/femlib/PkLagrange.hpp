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
// SUMMARY : Generic Pk Lagrange finite element class (Jan. 2008)
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-universite.fr

/*
 Thank to the ARN FF2A3 grant
 ref: ANR-07-CIS7-002-01
 */

#ifndef _PK_LAGRANGE_HPP_
#define _PK_LAGRANGE_HPP_
#include "splitsimplex.hpp"

#include "FESpacen.hpp"
#include <vector>
namespace Fem2D {


template<class Rd,class E>
static void SetPtPk(Rd *Pt, const int *dfon, int nn) {
    // P0, P1, P1b & P2
    typedef typename E::RdHat RdHat;
    const int dHat = E::RdHat::d;
    //const RdHat *KHat=RdHat::KHat;
    // sorry on some arch this is unset so rebuild ... FH 
    static RdHat KHat[dHat+1]; // Bug if no static on ubuntu ??? FH
    for(int i=0; i<=dHat;++i)
    KHat[i+1][i]=1. ;
    int k = 0;
    if (dfon[0]) {
        
        for (int i = 0; i <= dHat; ++i)
         Pt[k++] = KHat[i];
        
    }
    
    if (dfon[1] && dHat !=1)
        for(int i = 0; i < E::ne; ++i)
           Pt[k++] = (KHat[E::nvedge[i][0]] + KHat[E::nvedge[i][1]])*0.5;
    
    if (dfon[dHat] == 1)
        Pt[k++] = Rd::diag(1./(dHat+1));
    
    if (nn != k) {
        cout << " " << nn << " != " << k << " - dHat = " << dHat << " " << dfon[0] << dfon[1] << dfon[2] << dfon[3] << " " << E::ne << endl;
        ffassert(nn == k);
    }
    if (verbosity > 9)
        cout << " Pk = " << KN_<Rd>(Pt, nn) << "\n";
}
template<class Rd,class E>
static void SetPtPkDc(int kk,Rd *Pt, const int *dfon, int nn)
{
    typedef typename E::RdHat RdHat;
    const int dHat = E::RdHat::d;
    static RdHat KHat[dHat+1]; // Bug if no static on ubunti ??? FH
    // sorry on some arch this is unset so rebuild ... FH
    for(int i=0; i<=dHat;++i)
    KHat[i+1][i]=1. ;

    int k=0;
    if (kk>0) {
        int n= dfon[dHat];
        double ckk=double(kk);
        for(int dof=0; dof<n;++dof)
         {
           RdHat P;
           invNumSimplex(dof,P);
             P /= ckk;
           Pt[k++] = P;
           if (verbosity > 9)
             cout << k << " " << Pt[k-1] << endl;
         }
    }
    
    else  if (kk==-1 ) // edge
        for(int i = 0; i < E::ne; ++i)
           Pt[k++] = (KHat[E::nvedge[i][0]] + KHat[E::nvedge[i][1]])*0.5;
    else if (kk==-2)  // face
         for(int i = 0; i < E::nf; ++i)
               Pt[k++] = (KHat[E::nvface[i][0]] + KHat[E::nvface[i][1]]+ KHat[E::nvface[i][2]])/3.;
    else if (kk==-3) // vertices
        for(int i = 0; i < E::nv; ++i)
              Pt[k++] = KHat[i];

    if (nn != k)
     {
        cout << " " << nn << " != " << k << " - dHat = " << dHat << " " << dfon[0] << dfon[1] << dfon[2] << dfon[3] << " " << E::ne << endl;
         cout << " BIG BUG FH !!!!! "<< endl;
        ffassert(nn == k);
    }
    if (verbosity > 9)
        cout << " Pk = " << KN_<Rd>(Pt, nn) << "\n";
}
inline void invNumSimplex (int n,R1 &P){ P.x=n;}// Identity ..
inline void invNumSimplex (int n,R2 &P){int i1,i2; invNumSimplex2(n,i1,i2);P.x=i1;P.y=i2;}
inline void invNumSimplex (int n,R3 &P){int i1,i2,i3; invNumSimplex3(n,i1,i2,i3);P.x=i1;P.y=i2;P.z=i3;}

//  a class of Lagrange Pk finite element
template<class MMesh>
class TypeOfFE_Lagrange: public GTypeOfFE<MMesh> {
public:
    typedef MMesh Mesh;
    typedef typename Mesh::Element Element;
    typedef typename Element::RdHat RdHat;
    typedef typename Element::Rd Rd;

    static const int dHat=RdHat::d;
    struct A4 {
        int dfon[4];
        A4(int k,int ttdc=0) {
            dfon[0] = dfon[1] = dfon[2] = dfon[3] = 0;

            if (k == 0) {// P0
               
                dfon[dHat] = 1;
            }
            else  if( ttdc==0)
            {
                if (k == -1) { // P1b
                    dfon[0] = 1;
                    dfon[dHat] = 1;
                }
                else if (k == -2 && dHat==2) { // P2b
                    dfon[0] = 1;
                    dfon[1] = 1;
                    dfon[dHat] = 1;
                }

                else {
                    dfon[0] = 1;
                    dfon[1] = max(k - 1, 0);
                    dfon[2] = dHat > 1 ? max(k - 2, 0) : 0;
                    dfon[3] = dHat > 2 ? max(k - 3, 0) : 0;
                }
            }
            else if(ttdc==1)
            {
              if (k == -1 )  // P0edge
                     dfon[dHat] = Element::ne;
               else if (k == -2 )  // P0face
                     dfon[dHat] = Element::nf;
               else if (k == -3 )  // P0vertex
                     dfon[dHat] = Element::nv;
               else if(k>0) {
                
                int fd =1, nd=1;
                for(int i=0; i< dHat; ++i)
                {
                    fd *= i+1;// dHat!
                    nd *= (k+i+1); // (k+1)  .. (k+dHat)
                }
                //  ndof =  (k+1)  .. (k+dHat) / dHat! :
                //  d=3,  dim :P1 = 4, dim P2 10, dim P3 : 20
                dfon[dHat] =  nd/fd;
                if(verbosity>9)
                    cout << " TypeOfFE_Lagrange dim=" << dHat << " P_"<< k << " " << ttdc  << " = " <<  dfon[dHat] << " !!" << fd << " " << nd <<  endl;
                ffassert( nd % fd ==0);
               }
                
            }
            else  if(ttdc==2) //  Continious
               {
                   if (k == -1 )  // P0edge
                          dfon[1] = 1;
                    else if (k == -2 )  // P0face
                           dfon[2] = 1;
                    else if (k == -3 )  // P0vertex
                           dfon[0] = 1;
                    else ffassert(0);

               }
            else ffassert(0);
               
             if (verbosity > 9) // debile verbosity is set after (just for debug FH.)
                cout << " A4 d=" << dHat << " k= " << k <<  " " << ttdc << ":: " << dfon[0] << dfon[1] << dfon[2] << dfon[3] << endl;
        }
        operator const int * () const {return dfon;}
    };
    
    //RdHat *Pt;
    const RdHat G;
    R cshrink;
    R cshrink1;

    RdHat Shrink(const RdHat& P) const { return (P-G)*cshrink+G;}
    RdHat Shrink1(const RdHat& P)const { return (P-G)*cshrink1+G;}
    
    TypeOfFE_Lagrange(int k,int ttdc=0,R ccshrink=0):
    GTypeOfFE<Mesh>(A4(k,ttdc), 1, k == -1 ? -1 : Max(k, 1), k <= 2, (k == 0)|| ttdc ),
    cshrink(1.-ccshrink),cshrink1(1./cshrink),G(RdHat::diag(1./(RdHat::d+1))) {
        int n = this->NbDoF;
        if (verbosity > 9)
            cout << "\n +++ P" << k << ": ndof = " << n << " " << ttdc << " schrk: " << ccshrink << endl;
        ffassert(cshrink>0 && cshrink<=1.);
        if(ttdc>0)  // for discontinuous FE.
        {
            SetPtPkDc<RdHat, Element>(k,this->PtInterpolation, this->ndfOn(), this->NbDoF);
            for (int i = 0; i < n; i++)
              this->PtInterpolation[i]=Shrink(this->PtInterpolation[i]);
        }
        else
            SetPtPk<RdHat, Element>(this->PtInterpolation, this->ndfOn(), this->NbDoF);
        if (verbosity > 9)
            cout << this->PtInterpolation << endl;
        for (int i = 0; i < n; i++) {
            this->pInterpolation[i] = i;
            this->cInterpolation[i] = 0.;
            this->dofInterpolation[i] = i;
            this->coefInterpolation[i] = 1.;
        }
    }
    ~TypeOfFE_Lagrange(){}
    
private:
    TypeOfFE_Lagrange(const TypeOfFE_Lagrange &);
    void operator = (const TypeOfFE_Lagrange &);
};
/*
template<class Mesh>
const typename Mesh::Element::RdHat TypeOfFE_Lagrange<Mesh>::G=Mesh::Element::RdHat::diag(1./(Mesh::Element::RdHat::d+1)) ;
*/
template<class Mesh>
class TypeOfFE_LagrangeDC : public  TypeOfFE_Lagrange<Mesh> {
public:
    typedef typename Mesh::Element Element;
    typedef  Element E;
    typedef typename Element::RdHat RdHat;
    typedef typename Element::Rd Rd;
    typedef GFElement<Mesh> FElement;
    constexpr static int  dHat =  RdHat::d;
    constexpr static int  d =  Rd::d;
    const int k;
    
    vector<double> ml,mc;
    vector<int> mi;
    TypeOfFE_LagrangeDC(int kk,double ccskrink=0):
    TypeOfFE_Lagrange<Mesh>(kk,1,ccskrink),k(kk),
    ml(this->NbDoF*k),mc(this->NbDoF*k),mi(this->NbDoF*k)
    {
        // we have  k monome pedof
        //  where l[d+1]  is the barycentric coor..
        //  m(l)_i = (ml[i]*l[mi[i]]+mc[i])
        // fb_j(l)  = Prod_{i in [k*j,k*(j+1)[}  m(l)_i
        const int n = this->NbDoF;
        
        int km=0;
        for(int dof=0; dof<n;++dof)
        {
            RdHat P;
            int ijk[dHat+1];
          
            invNumSimplex(dof,P);
            long l, s=0;
            for(int j=0; j<dHat;++j)
            {
                l = lround(P[j]);
                s+=l;
                ijk[j+1]=l;
            }
            ijk[0]=k-s;
            // construction of the monome
            int m=0;
            for(int li=0;li<=dHat;++li)
            {
                int ii=ijk[m++];
                for(int i=0; i<ii; ++i)
                {
                    //  (k*l[li]-i)/(ii-i);;
                    mi[km] = li;
                    ml[km] = double(k)/(ii-i);
                    mc[km] = -double(i)/(ii-i);
                    km++;
                }
            }
            
            
        }
        if(verbosity>9)
        for(int dof=0,m=0; dof<n;++dof)
        {
            cout.setf(ios::showpos);
            cout << dof << " : " ;
            for(int l=0; l<k;++l,++m)
            cout << "("<< ml[m]<< " l_" << mi[m] << " +"<<  mc[m] <<")" ;
            cout << endl;
            cout.unsetf(ios::showpos);
        }
        ffassert(km==k*this->NbDoF);
        
        
    }
    
    void FB(const What_d whatd,const Mesh & Th,const Element & K,const RdHat &PHat1, RNMK_ & val) const
    {
        RdHat PHat= TypeOfFE_Lagrange<Mesh>::Shrink1(PHat1);
        double l[dHat+1];
        PHat.toBary(l);
        // l = 0 en i/k , 1 en ii/k  => monome : ( k*l - i )/( ii-i) ok..
        RN_ f0(val('.',0,op_id));
        if (whatd & Fop_D0)
            for( int dof = 0,m=0; dof < this->NbDoF; ++dof)
            {
            double f=1.;
            for(int i=0; i<k;++i,++m)
            f *= (ml[m]*l[mi[m]]+mc[m]);
            f0[dof]= f;
            }
        else if(whatd & (Fop_D1|Fop_D2))
          {
            bool d2= Fop_D2 && k>1;//  calcul de
            const unsigned int fop[3]={Fop_dx,Fop_dy,Fop_dz};
            const  int op[3]={op_dx,op_dy,op_dz};
              const  int dop[9]={op_dxx,op_dxy,op_dxz, op_dyx,op_dyy,op_dyz, op_dzx,op_dzy,op_dzz};

            Rd Dl[dHat+1];
            Rd DDl[dHat+1][d];
            K.Gradlambda(Dl);
            KN<Rd> df(this->NbDoF),ddf(this->NbDoF*dHat);
            
            for( int dof = 0,m=0; dof < this->NbDoF; ++dof)
                {
                double f=1.;
                Rd Df;
                Rd DDf[d] ;
                for(int i=0; i<k;++i,++m)
                    { // f = f*b // df = df*b+db*f; ddf  = ddf*b + 2*df*db + f * ddb
                     int im=mi[m];
                     double b =(ml[m]*l[im]+mc[m]);
                     Rd Db=this->cshrink1*ml[m]*Dl[im];
                     if(d2)
                        for(int l=0; l<d;++l)
                           DDf[l] = b*DDf[l]+ Db[l]*Df + Db*Df[l] ;
                     Df = b*Df+ (f*Db);
                     f *= b;
                    }
                  df[dof] = Df;
                  if(d2)
                  {
                      for(int l=0; l<d;++l)
                         ddf[dof*3+l] = DDf[l];
                  }
                }
            // copy data d  D
            for (int dd=0; dd< Rd::d;++dd)
                {
                    if (whatd & fop[dd])
                    {
                        RN_ dfdd(val('.',0,op[dd]));
                        for(int i=0;i<this->NbDoF;++i)
                          dfdd[i]= df[i][dd];
                    }
                }
            // copy data  DD
            if(d2)
                {
                   for (int id=0; id< Rd::d;++id)
                    for (int jd=0; jd<= id;++jd)
                    {
                        int op = dop[id*3+jd];
                        const unsigned int fopij= 1<< op;
                        if (whatd & fopij)
                        {
                          RN_ dfdd(val('.',0,op));
                          for(int i=0;i<this->NbDoF;++i)
                           dfdd[i]= ddf[i*3+id][jd];
                        }
                    }
                }
            
         }
    }
    
} ;

template<class Mesh>
class TypeOfFE_ConstDC : public  TypeOfFE_Lagrange<Mesh> {
public:
    typedef typename Mesh::Element Element;
    typedef  Element E;
    typedef typename Element::RdHat RdHat;
    typedef typename Element::Rd Rd;
    typedef GFElement<Mesh> FElement;
    constexpr static int  dHat =  RdHat::d;
    constexpr static int  d =  Rd::d;
    const int k;
    //  dc == 2 : some continuite, dc ==1 : just on element
    TypeOfFE_ConstDC(int kk,int dc):
    TypeOfFE_Lagrange<Mesh>(-kk,dc,(dc==1) ? 0.001:0)  ,k(kk)
    {
        // k = 1 const close to edge
        // k = 2 const close to  face
        // k = 3 const close to  vertex

    }
    
    void FB(const What_d whatd,const Mesh & Th,const Element & K,const RdHat &PHat, RNMK_ & val) const
    {
        double l[d+1];
        PHat.toBary(l);
        // l = 0 en i/k , 1 en ii/k  => monome : ( k*l - i )/( ii-i) ok..
        RN_ f0(val('.',0,op_id));
        if (whatd & Fop_D0)
            if(k==1) // edge
            {   // close edge ie = arg min_e l[e2]+l[e1]
                int ie = 0;
                double lm = 2;
                for(int i=0;i<Element::ne; ++i)
                {
                    double li = (l[E::nvedge[i][0]] + l[E::nvedge[i][1]]);
                    if( li < lm)  { ie= i;lm =li;}
                }
                f0[ie]=1;
            }
            else if(k==2) // face
            {
                // close face ie = arg min_i l[i]
                int ie = 0;
                double lm = 2;
                for(int i=0;i<Element::nv; ++i)
                {
                    double li = l[i];
                    if(li < lm){ ie= i;lm =li;}
                }
                f0[ie]=1;
                
            }
            else if(k==3) // vertex
            {
                // close vertex =  // close vertex ie = arg max_i l[i]
                int ie = 0;
                double lm = -2;
                for(int i=0;i<Element::nv; ++i)
                {
                    double li = l[i];
                    if( li > lm){ ie= i;lm =li;}
                }
                f0[ie]=1;
                
            }
            else
                ffassert(0);
    }
                            
    
} ;
    
} // namespace Fem2D


#endif //_PK_LAGRANGE_HPP_
