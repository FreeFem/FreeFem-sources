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
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

#include "ff++.hpp"
#include "AddNewFE.h"
#include <iostream>
// Attention probleme de numerotation des inconnues
// -------------------------------------------------
// dans freefem, il y a un noeud par objets  sommet, arete, element.
// et donc la numerotation des dl dans l'element depend
// de l'orientation des aretes
//
/// ---------------------------------------------------------------
namespace Fem2D {

  class TypeOfFE_P2pnc_3d : public GTypeOfFE< Mesh3 > {
   public:
   typedef Mesh3 Mesh;
   typedef Mesh3::Element Element;
   typedef GFElement< Mesh3 > FElement;
      static int dfon[];
      static const int d = Mesh::Rd::d;
      static const GQuadratureFormular< R2 > &QFf;    // quadrature formula on a face
      static const GQuadratureFormular< R3 > &QFt;    // quadrature formula on a face
      TypeOfFE_P2pnc_3d( );                      // constructor
      void FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K, const RdHat &PHat,
              RNMK_ &val) const;
      void set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M, int ocoef, int odf,
               int *nump) const;
    };
int TypeOfFE_P2pnc_3d::dfon[] = {0, 0, 3, 1};    // 2 dofs on each edge, 2 dofs on each face


const GQuadratureFormular< R2 > &TypeOfFE_P2pnc_3d::QFf= QuadratureFormular_T_5 ;
const GQuadratureFormular< R3 > &TypeOfFE_P2pnc_3d::QFt= QuadratureFormular_Tet_5 ;

TypeOfFE_P2pnc_3d::TypeOfFE_P2pnc_3d( )
: GTypeOfFE< Mesh3 >(TypeOfFE_P2pnc_3d::dfon, 1, 3,
                     4 * QFf.n * 3 +  QFt.n,// N coef InterPP
                     4 * QFf.n +  QFt.n,// N Pt Inter
                     false, true) {
   static  R3 Pt[] = {R3(0., 0., 0.), R3(1., 0., 0.), R3(0., 1., 0.),
               R3(0., 0., 1.)};    // 4 ref tetrahedron vertices
   static const int  nvfo[4][3]  ={{1,2,3}, {0,2,3},{0,1,3},{ 0,1,2}};
    if(verbosity > 0 && mpirank == 0) cout << "TypeOfFE_P2pnc_3d QFf exact:"<< QFf.exact << ", QFt exact " <<QFt.exact<< endl;
    int ipt=0;
    int doft = 4*3; // last dof
    
    // construction de point d'interpolation sur les face
    for (int f = 0; f < Element::nf; ++f) {
        for (int q = 0; q < QFf.n; ++q, ++ipt) {
            double x = QFf[q].x;
            double y = QFf[q].y;
            this->PtInterpolation[ipt] = Pt[nvfo[f][0]] * (1. - x - y) + Pt[nvfo[f][1]] * x + Pt[nvfo[f][2]] * y;
        }}
    
    for (int q = 0; q < QFt.n; ++q, ++ipt)
      this->PtInterpolation[ipt] = QFt[q];

    ffassert(ipt == this->NbPtforInterpolation); // verif
    if(verbosity>99 && mpirank == 0)
    for( int i=0; i< ipt;++i)
        cout << i << " P/i " <<  this->PtInterpolation[i] << endl;
    {
        int i =0,ipt=0; //
        
        
        //  coef d'interpolation ... sur les face ...
        //  int li*f
        for (int f = 0; f < Element::nf; ++f)
            for (int qq = 0; qq < QFf.n; qq++,++ipt)
            {
              double ll[4]; // dans Khat
              this->PtInterpolation[ipt].toBary(ll);

              for (int kf = 0; kf < 3; ++kf,i++)
              {
                  int kfK = nvfo[f][kf];// vertex dans K ..
                  int dof = 3*f+kf;
                  if(qq==0)
                  if (verbosity > 0 && mpirank == 0) cout << " dof " << dof << " " << f << " "<< kfK << endl;
                  {
                      this->pInterpolation[i] = ipt;          // pk in (13.1)
                      this->cInterpolation[i] = 0;          // jk in (13.1)
                      this->dofInterpolation[i] = dof;        // ik in (13.1)
                      this->coefInterpolation[i] = QFf[qq].a*ll[kfK];
                      
                  }
              }
        }
        
        for (int q = 0; q < QFt.n; ++q, ++i,ipt++)
        {
            this->pInterpolation[i] = ipt;          // pk in (13.1)
            this->cInterpolation[i] = 0;          // jk in (13.1)
            this->dofInterpolation[i] = doft;        // ik in (13.1)
            this->coefInterpolation[i] = QFt[q].a;      // alfak: we will fill them with 'set' (below)
        }

        ffassert(i==this->NbcoefforInterpolation);
        
    }
}

void Setp3(int *p,int n) {
// n =  signe + depart*2
    int i=n/2, j= n%2, k = 3-i;
    if( i==0){ if(j==1) swap(p[1],p[2]); }
    else {
        swap(p[0],p[i]);
        int k =3 -i;// autre
        if(j==0) swap(p[i],p[k]);
        }
}
void TypeOfFE_P2pnc_3d::set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M,
                            int ocoef, int odf, int *nump) const {
    static int count =0;count++;
    int k = Th(K);
    int verb = (mpirank == 0 && (verbosity>99 || (verbosity>9 && count < 4 )));//||  ocoef+odf || nump  ;
    int n = this->NbDoF;
    int np=M.np;
    int ncoef=M.ncoef;
    //int *p = M.p+odf;// correction FH mai 2020 ...
    //  update the coef to be ajdacent compatible
    //
    int i =ocoef; //
    
 //   cout << " set " << Th(K) << endl;
    
    //  coef d'interpolation ... sur les face ...
    //  int li*f
    int ip=0;
    for (int f = 0; f < Element::nf; ++f)
    {
        int fp = K.facePermutation(f);// orientation de la face ..
        int p3[3]={0,1,2};
        SetNumPerm<3>(fp,p3);
        if(verb ) cout << " set: " << f << " " << p3[0] << " "<< p3[1] << " " << p3[2]
           << "  ::: " << Th(k,Element::nvface[f][p3[0]]) << " " <<  Th(k,Element::nvface[f][p3[1]]) << " " <<  Th(k,Element::nvface[f][p3[2]]) << " " <<  endl;
        for (int qq = 0; qq < QFf.n; qq++,ip++)
        {
            
            int ipt = nump ? nump[ip] : ip;
            ffassert(ipt<=np);
            //cout << ip << " -> "<< ipt << "  f= " << f <<" " << qq << " i "<< i << " " << this->pInterpolation[i]
            //<< "ipt: " << this->PtInterpolation[ipt] << " ::ip  "<< this->PtInterpolation[ip] <<endl;
            //ffassert( ipt == this->pInterpolation[i]) ; ce ne marche pas !!!!!! mais le point sont les  memem
            double ll[4]; // dans Khat
            M.P[ipt].toBary(ll);//  point sur la face f ????
            if(verb) cout << " P " << M.P[ipt] << " " << i  << endl;
            for (int kf = 0; kf < 3; ++kf,i++)
            {
                if(verb) cout << i <<  " " <<f<< kf  <<  " " <<  ipt << " " << ll[f] << " " << M.p[i] << " // " << ocoef << " " << odf << " " << nump << endl;
               // ffassert(ipt == this->pInterpolation[i]);
                ffassert(abs(ll[f])<1e-10);

                int kfK = Element::nvface[f][p3[kf]];// vertex dans K ..
                M.coef[i] = QFf[qq].a*ll[kfK];
             }
        }
    }
    
    for (int q = 0; q < QFt.n; ++q, ++i)
    {
        M.coef[i] = QFt[q].a;      // alfak: we will fill them with 'set' (below)
    }
   
    ffassert(i<= ncoef);

}
void TypeOfFE_P2pnc_3d::FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K,
                           const RdHat &PHat, RNMK_ &val) const {
    
    assert(val.N( ) >= 13);    // 13 degrees of freedom
    assert(val.M( ) == 1);     // 1 components
/*
   Dof  Numbering  3 dof / face   the dof a associated to a vertex face
  and the vertex a numbering in increase way.
 
 we have to numbering    î on ref element with
 
   So we need of a permutation p to insure the compatibilit between adjacent element
   
 */
    //  generated  with file Element_P2pnc_3d.edp
    // Warning  p(^i) =i
   // CC(i,j) = dof(j)(mo(i)); //
   // phi_k =sum_i C1(k,i) mo_i
    double C1[13][13] = {
      {9, -15, -15, 3, -60, -60, 60, -30, 30, 30, 0, 180, 180},
       {0, 0, 18, 0, 0, 30, -30, 0, 0, -30, 0, -180, 0},
       {0, 18, 0, 0, 30, 0, -30, 0, -30, 0, 0, 0, -180},
  
      {9, -9, -3, -3, -90, 0, 30, 0, 30, 0, 0, 0, 180},
       {0, 0, 18, 0, 0, -30, 30, 0, 0, -30, -180, 0, 0},
       {-6, 18, -12, 6, 60, 30, -90, 0, -60, 30, 180, 0, -180},
  
      {-3, -3, 27, 9, 0, 30, 30, 0, 0, -90, -180, -180, 0},
       {-3, 9, -9, -3, 0, 30, -90, 0, 0, 30, 180, 0, 0},
       {9, -3, -9, -3, 0, -90, 30, 0, 0, 30, 0, 180, 0},
  
      {0, 12, 12, 6, 60, 60, -60, -30, -30, -30, 0, -180, -180},
       {3, -9, 9, 3, -60, -30, 90, 0, 0, -30, -180, 0, 180},
       {0, 0, -18, 0, 0, -30, -30, 0, 0, 30, 180, 180, 0},
  
      {-5, -5, -5, -5, 20, 20, 20, 20, 20, 20, 0, 0, 0}};


    int pp[13]={0,1,2,3,4,5,6,7,8,9,10,11,12}; // Permutaion de dof
    int p[13]; // Permutaion de dof
    p[12]=pp[12];
    static int ccount =0; ccount++;
    int verb = (ccount < 2 && verbosity>9 ) || ( verbosity>99);
    for (int ff = 0,k=0; ff < Element::nf; ff++, k+=3 ) {
      // orientation de la face a envert
        int fp = K.facePermutation(ff);
        int p3[3]={0,1,2};
        SetNumPerm<3>(fp,p3);
        if( verb)
        {
            int i3[3]={Th(K[Element::nvface[ff][p3[0]]]),Th(K[Element::nvface[ff][p3[1]]]),Th(K[Element::nvface[ff][p3[2]]])};
            cout << "k " << k <<" " << ff <<  " / "<< p3[0] << " "<< p3[1] << " "<< p3[2]
            << " / " << i3[0] << " " << i3[1] << " " << i3[2]
            //   << " / " << Element::nvface[ff][p3[0]]<< " " <<Element::nvface[ff][p3[1]]<< " " <<Element::nvface[ff][p3[2]]<< " "
            //    << " / " << Element::nvface[ff][0]<< " " <<Element::nvface[ff][1]<< " " <<Element::nvface[ff][2]<< " "
            <<  endl;
        }
        p[k+0] = pp[k+p3[0]];;
        p[k+1] = pp[k+p3[1]];
        p[k+2] = pp[k+p3[2]];
    }
    if( verb )
    {  cout << "\n K = " << Th(K) << " p : " ;
        for(int i=0; i<13;++i)
            cout << p[i] << ",";
        cout << endl;
    }
    double l[4],mo[13];
    PHat.toBary(l);
    {// construction des monone
        int k=0;
        for(int i=0;i<4;i++)
            mo[k++]=l[i];
        for(int i=0;i<4;i++)
            for(int j=0;j<i;j++)
                mo[k++] = l[i]*l[j];
        mo[k++] = l[1]*l[2]*l[2];
        mo[k++] = l[0]*l[2]*l[2];
        mo[k++] = l[0]*l[1]*l[1];
    }
    RN_ f0(val('.', 0, op_id));
    if (whatd & Fop_D0) {
        double v[13];
        // v = c1*mo;
        for(int i=0; i<13;++i)
            v[i]=0.;
        for(int i=0; i<13;++i)
            for(int j=0; j<13;++j)
                v[i]+= C1[p[i]][j]*mo[j];
       // do pernumation
        for(int i=0;i<13;++i)
            f0[i]= v[i];
    }
    if (whatd & (Fop_D1 | Fop_D2)) {
        const unsigned int iDx =whatd & Fop_dx , iDy =whatd & Fop_dy , iDz =whatd & Fop_dz; 
        R3 Dl[4];
        K.Gradlambda(Dl);
        R3 Dmo[13],Dv[13];
   //     static int ccc=0;
   //     if(ccc++ < 100) cout << " IDxyz =" << iDx<< " " <<iDy << " " << iDz << " " << whatd  << " " <<Dl[0] <<  endl;
        int k=0;
        for(int i=0;i<4;i++)
            Dmo[k++]=Dl[i];
        for(int i=0;i<4;i++)
            for(int j=0;j<i;j++)
                Dmo[k++] = Dl[i]*l[j]+ l[i]*Dl[j];
        Dmo[k++] = Dl[1]*(l[2]*l[2]) + (2.*l[1]*l[2])*Dl[2];
        Dmo[k++] = Dl[0]*(l[2]*l[2]) + (2.*l[0]*l[2])*Dl[2];
        Dmo[k++] = Dl[0]*(l[1]*l[1]) + (2.*l[0]*l[1])*Dl[1];
        ffassert(k==13);
        for(int i=0; i<13;++i)
            Dv[i]=R3();
        for(int i=0; i<13;++i)
            for(int j=0; j<13;++j)
                Dv[i]+= C1[p[i]][j]*Dmo[j];
        RN_ f0x(val('.', 0, op_dx));
        RN_ f0y(val('.', 0, op_dy));
        RN_ f0z(val('.', 0, op_dz));
        for(int i=0; i<13;++i)
        {
            if(iDx) f0x[i] = Dv[i].x;
            if(iDy) f0y[i] = Dv[i].y;
            if(iDz) f0z[i] = Dv[i].z;
        }
        if (whatd & ( Fop_D2))
             ffassert(0); // do do ..
    }

}
 static TypeOfFE_P2pnc_3d P2pnc_3d;
 GTypeOfFE< Mesh3 > &Elm_P2pnc_3d(P2pnc_3d);

 static AddNewFE3 TFE_P2pnc_3d("P2pnc3d", &Elm_P2pnc_3d);
  
}    // namespace Fem2D

// --- fin --
