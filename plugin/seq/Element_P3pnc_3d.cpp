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

  class TypeOfFE_P3pnc_3d : public GTypeOfFE< Mesh3 > {
   public:
   typedef Mesh3 Mesh;
   typedef Mesh3::Element Element;
   typedef GFElement< Mesh3 > FElement;
      static int dfon[];
      static const int d = Mesh::Rd::d;
      static const GQuadratureFormular< R2 > &QFf;    // quadrature formula on a face
      static const GQuadratureFormular< R3 > &QFt;    // quadrature formula on a face
      TypeOfFE_P3pnc_3d( );                      // constructor
      void FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K, const RdHat &PHat,
              RNMK_ &val) const;
      void set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M, int ocoef, int odf,
               int *nump) const;
    };
int TypeOfFE_P3pnc_3d::dfon[] = {0, 0, 6, 4}; //0 dof on node, 0 on edge, 6 on face, 4 on element//

const GQuadratureFormular< R2 > &TypeOfFE_P3pnc_3d::QFf= QuadratureFormular_T_7 ;
const GQuadratureFormular< R3 > &TypeOfFE_P3pnc_3d::QFt= QuadratureFormular_Tet_5 ;

TypeOfFE_P3pnc_3d::TypeOfFE_P3pnc_3d( )
: GTypeOfFE< Mesh3 >(TypeOfFE_P3pnc_3d::dfon, 1, 3,
                     4 * QFf.n * 6 + 4 * QFt.n,// N coef InterPP
                     4 * QFf.n + 1 * QFt.n,// N Pt Inter
                     false, true) {
   static  R3 Pt[] = {R3(0., 0., 0.), R3(1., 0., 0.), R3(0., 1., 0.),
               R3(0., 0., 1.)};    // 4 ref tetrahedron vertices
   static const int  nvfo[4][3]  ={{1,2,3}, {0,2,3},{0,1,3},{ 0,1,2}};
    cout << "TypeOfFE_P3pnc_3d QFf exact:"<< QFf.exact << ", QFt exact " <<QFt.exact<< endl;
    int ipt=0;
    int doft0 = 4*6; // nb dof on face

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
    if(verbosity>99)
    for( int i=0; i< ipt;++i)
        cout << i << " P/i " <<  this->PtInterpolation[i] << endl;
    {
        int i =0,ipt=0; //


        //  coef d'interpolation ... sur les face ...
        //  int li*f
        for (int f = 0; f < Element::nf; ++f)   //loop on the 4 faces
            for (int qq = 0; qq < QFf.n; qq++,++ipt)
            {
              double ll[4]; // dans Khat
              this->PtInterpolation[ipt].toBary(ll);

              for (int kf = 0; kf < 3; ++kf,i++)
              {
                  int kfK = nvfo[f][kf];// vertex dans K ..
                  int dof = 6*f+kf;
                  if(qq==0)
                  cout << " dof " << dof << " " << f << " "<< kfK << endl;
                  {
                      this->pInterpolation[i] = ipt;          // pk in (13.1)
                      this->cInterpolation[i] = 0;          // jk in (13.1)
                      this->dofInterpolation[i] = dof;        // ik in (13.1)
                      this->coefInterpolation[i] = QFf[qq].a*ll[kfK];
                  }
              }

              for (int kf = 0; kf < 3; ++kf,i++)
              { //edge next  of vertex .. kf. pas terrible pour la numerotation !!!
                  int kfK = nvfo[f][(kf+2)%3];// vertex dans K opose a kf
                  int kfK1 = nvfo[f][(kf+1)%3];// vertex dans K ..
                  int dof = 6*f+kf+3;
                  if(qq==0)
                  cout << " dof " << dof << " " << f << " "<< kfK << endl;
                  {
                      this->pInterpolation[i] = ipt;          // pk in (13.1)
                      this->cInterpolation[i] = 0;          // jk in (13.1)
                      this->dofInterpolation[i] = dof;        // ik in (13.1)
                      this->coefInterpolation[i] = QFf[qq].a*ll[kfK]*ll[kfK1];
                  }
              }

          } // end loop on face
        cout << " i =" << i << " " <<QFf.n*6*4 <<  endl;
        for (int q = 0; q < QFt.n; ++q,ipt++) //loop on quadrature point on element
        {  double ll[4]; // dans Khat
            this->PtInterpolation[ipt].toBary(ll);
          for (int kt = 0; kt < 4; ++kt,i++) //loop on the 4 dof in the element
          {
            this->pInterpolation[i] = ipt;          // pk in (13.1)
            this->cInterpolation[i] = 0;          // jk in (13.1)
            this->dofInterpolation[i] = doft0+kt;        // ik in (13.1)
            this->coefInterpolation[i] = QFt[q].a*ll[kt];      // alfak: we will fill them with 'set' (below)
          }
        }
        cout << i << "== " << this->NbcoefforInterpolation << " " << QFt.n*4 + QFf.n*6*4 << endl;;
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
void TypeOfFE_P3pnc_3d::set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M,
                            int ocoef, int odf, int *nump) const {
    static int count =0;count++;
    int k = Th(K);
    int verb = verbosity>99 || (verbosity>9 && count < 4 ) ;//||  ocoef+odf || nump  ;
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
            //func p2 on face
             for (int kf = 0; kf < 3; ++kf,i++){
                 int kfK = Element::nvface[f][p3[(kf+2)%3]];// vertex dans K ..
                 int kfK1 = Element::nvface[f][p3[(kf+1)%3]];// vertex dans K ..
                 M.coef[i] = QFf[qq].a*ll[kfK]*ll[kfK1];
            }

        } //end loop en quadrature
    } //end loop on the four faces

    double ll[4]; // dans Khat

    for (int q = 0; q < QFt.n; ++q, ++ip)
    {   int ipt=ip;
        double ll[4]; // dans Khat
        M.P[ipt].toBary(ll);//
        M.coef[i++] = QFt[q].a*ll[0];
        M.coef[i++] = QFt[q].a*ll[1];
        M.coef[i++] = QFt[q].a*ll[2];
        M.coef[i++] = QFt[q].a*ll[3];        // alfak: we will fill them with 'set' (below)
    }

  //  ffassert(i==  ncoef+ocoef);

}
void TypeOfFE_P3pnc_3d::FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K,
                           const RdHat &PHat, RNMK_ &val) const {

    assert(val.N( ) >= 28);    // 28 degrees of freedom
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
    
    double C1[28][28] = {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -150, -120, 30, -120, 60, 30, 960, -120, -120, 60, 0, 0, 0, -1680, 0, 0, 0, 0},
         {0, 0, 0, 0, -30, 120, 150, 120, -960, -60, 0, 0, 0, 120, 300, -960, 0, 0, 150, 0, 0, 0, 1680, 0, 0, 0, 0, 1680},
         {0, -150, 960, 60, 0, -120, -120, 0, 30, 0, 0, -120, -120, 0, 60, 0, 0, 30, 0, 0, 0, 0, 0, 0, -1680, 0, 0, 0},
    
        {-96, 882, -3168, -156, 522, -576, -18, -648, 2862, 204, -2628, 504, 432, -216, -936, 2952, 2592, -378, -738, 84, 0, 0, -5040, -3360, 5040, 3360, 0, -5040},
         {144, 2562, -3408, -156, -708, 504, 672, 672, -2058, -156, 2562, -336, 342, 504, 264, -2058, -3408, 342, 672, -156, -3360, 0, 3360, 5040, 5040, -3360, 0, 3360},
         {-96, -2628, 2592, 84, 522, -216, -738, -648, 2952, 204, 882, 504, -378, -576, -936, 2862, -3168, 432, -18, -156, 3360, 0, -5040, 5040, -3360, 0, 0, -5040},
    
        {60, 960, -150, 0, -120, -120, 0, 30, 0, 0, -120, -120, 0, 60, 0, 0, 30, 0, 0, 0, -1680, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 30, -120, -150, -120, 960, 60, 0, 0, 0, 60, -120, -120, 0, 0, 30, 0, 0, 0, 0, 0, 0, 0, 0, -1680},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 150, 120, -30, 300, 120, 150, -960, 120, -960, -60, 0, 0, 0, 1680, 0, 0, 1680, 0},
    
        {84, 2592, -2628, -96, -378, 504, 882, 432, -3168, -156, -738, -216, 522, -936, -576, -18, 2952, -648, 2862, 204, -3360, 0, 0, -5040, 3360, 0, -5040, 5040},
         {-156, -3168, 882, -96, 432, 504, -2628, -378, 2592, 84, -18, -576, 522, -936, -216, -738, 2862, -648, 2952, 204, 5040, 3360, 0, -5040, 0, 0, -5040, -3360},
         {-156, -3408, 2562, 144, 342, -336, 2562, 342, -3408, -156, 672, 504, -708, 264, 504, 672, -2058, 672, -2058, -156, 5040, -3360, 0, 3360, -3360, 0, 3360, 5040},
    
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 30, 60, 30, -120, -120, -150, -120, -120, 960, 60, 0, 0, 0, 0, 0, 0, -1680, 0},
         {0, 30, -120, 60, 0, -120, 960, 0, -150, 0, 0, 60, -120, 0, -120, 0, 0, 30, 0, 0, 0, -1680, 0, 0, 0, 0, 0, 0},
         {-60, -960, 150, 0, 120, 120, 0, -30, 0, 0, -960, 300, 0, 120, 0, 0, 150, 0, 0, 0, 1680, 0, 0, 0, 0, 1680, 0, 0},
    
        {204, 2862, -18, -156, -648, -576, -3168, 522, 882, -96, 2952, -936, 432, -216, 504, -2628, -738, -378, 2592, 84, -5040, 5040, 3360, 0, 0, -5040, -3360, 0},
         {204, 2952, -738, 84, -648, -216, 2592, 522, -2628, -96, 2862, -936, -378, -576, 504, 882, -18, 432, -3168, -156, -5040, -3360, 0, 0, 0, -5040, 5040, 3360},
         {-156, -2058, 672, -156, 672, 504, -3408, -708, 2562, 144, -2058, 264, 342, 504, -336, 2562, 672, 342, -3408, -156, 3360, 5040, -3360, 0, 0, 3360, 5040, -3360},
    
        {60, -120, 30, 0, -120, 60, 0, 30, 0, 0, 960, -120, 0, -120, 0, 0, -150, 0, 0, 0, 0, 0, 0, 0, 0, -1680, 0, 0},
         {0, 150, -960, -60, 0, 300, -960, 0, 150, 0, 0, 120, 120, 0, 120, 0, 0, -30, 0, 0, 0, 1680, 0, 0, 1680, 0, 0, 0},
         {0, 0, 0, 0, 30, 60, 30, -120, -120, 60, 0, 0, 0, -120, -120, 960, 0, 0, -150, 0, 0, 0, -1680, 0, 0, 0, 0, 0},
    
        {84, -738, 2952, 204, -378, -936, 2862, 432, -18, -156, 2592, -216, -648, 504, -576, -3168, -2628, 522, 882, -96, 0, -5040, 5040, 3360, -5040, -3360, 0, 0},
         {-156, 672, -2058, -156, 342, 264, -2058, 342, 672, -156, -3408, 504, 672, -336, 504, -3408, 2562, -708, 2562, 144, 0, 3360, 5040, -3360, 3360, 5040, -3360, 0},
         {-156, -18, 2862, 204, 432, -936, 2952, -378, -738, 84, -3168, -576, -648, 504, -216, 2592, 882, 522, -2628, -96, 0, -5040, -3360, 0, -5040, 5040, 3360, 0},
    
        {-48, 316, -904, -8, 236, 192, -44, -344, 596, 72, 316, 352, 36, 192, -448, 596, -904, 36, -44, -8, 0, 0, -1120, 1120, 1120, 0, 0, -1120},
         {-8, -904, 316, -48, 36, 352, 316, 36, -904, -8, -44, 192, 236, -448, 192, -44, 596, -344, 596, 72, 1120, 0, 0, -1120, 0, 0, -1120, 1120},
         {72, 596, -44, -8, -344, 192, -904, 236, 316, -48, 596, -448, 36, 192, 352, 316, -44, 36, -904, -8, -1120, 1120, 0, 0, 0, -1120, 1120, 0},
    
        {-8, -44, 596, 72, 36, -448, 596, 36, -44, -8, -904, 192, -344, 352, 192, -904, 316, 236, 316, -48, 0, -1120, 1120, 0, -1120, 1120, 0, 0}};



    int pp[28]={0,1,2,3,4,5,6,7,8,9,10,11,12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27}; // Permutaion de dof
    int p[28]; // Permutaion de dof
    for(int k=24;k<28;++k) // internal dof !!!
      p[k]=pp[k];
    static int ccount =0; ccount++;
    int verb = (ccount < 2 && verbosity>9 ) || ( verbosity>99);
    for (int ff = 0,k=0; ff < Element::nf; ff++, k+=6 ) {// face dof ..
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
        p[k+3] = pp[k+3+p3[0]];;
        p[k+4] = pp[k+3+p3[1]];
        p[k+5] = pp[k+3+p3[2]];
    }
    if( verb )
    {  cout << "\n K = " << Th(K) << " p : " ;
        for(int i=0; i<28;++i)
            cout << p[i] << ",";
        cout << endl;
    }
    // verif ... FH
    KN<int> p1(28);
    p1 = -1;
    for( int i=0; i< 28;++i)
        p1[p[i]]=1;
    ffassert(p1.min()>=0);
    double l[4],mo[28];
    PHat.toBary(l);
    {// construction des monone
      int k=0;
      for(int i=0;i<4;i++){
        for(int j=0;j<i+1;j++){
          for(int s=0;s<j+1;s++){
            mo[k++]=l[i]*l[j]*l[s];}}}  //hierarchical basis of P3, all monomes of degree 3
        mo[k++] = l[0]*l[0]*l[0]*l[1];
        mo[k++] = l[1]*l[1]*l[1]*l[2];
        mo[k++] = l[2]*l[2]*l[2]*l[3];
        mo[k++] = l[3]*l[3]*l[3]*l[0];
        mo[k++] = l[1]*l[1]*l[1]*l[0];
        mo[k++] = l[0]*l[0]*l[0]*l[3];
        mo[k++] = l[3]*l[3]*l[3]*l[2];
        mo[k++] = l[2]*l[2]*l[2]*l[1];
    }
    RN_ f0(val('.', 0, op_id));
    if (whatd & Fop_D0) {
        double v[28];
        // v = c1*mo;
        for(int i=0; i<28;++i)
            v[i]=0.;
        for(int i=0; i<28;++i)
            for(int j=0; j<28;++j)
                v[i]+= C1[p[i]][j]*mo[j];
       // do pernumation
        for(int i=0;i<28;++i)
            f0[i]= v[i];
    }
    if (whatd & (Fop_D1 | Fop_D2)) {
        const unsigned int iDx =whatd & Fop_dx , iDy =whatd & Fop_dy , iDz =whatd & Fop_dz;
        R3 Dl[4];
        K.Gradlambda(Dl);
        R3 Dmo[28],Dv[28];
   //     static int ccc=0;
   //     if(ccc++ < 100) cout << " IDxyz =" << iDx<< " " <<iDy << " " << iDz << " " << whatd  << " " <<Dl[0] <<  endl;
        int k=0;
        for(int i=0;i<4;i++){
          for(int j=0;j<i+1;j++){
           for(int s=0;s<j+1;s++){
             Dmo[k++]=Dl[i]*l[j]*l[s]+l[i]*Dl[j]*l[s]+l[i]*l[j]*Dl[s];}}}


        Dmo[k++] = 3*Dl[0]*l[0]*l[0]*l[1] + l[0]*l[0]*l[0]*Dl[1];
        Dmo[k++] = 3*Dl[1]*l[1]*l[1]*l[2] + l[1]*l[1]*l[1]*Dl[2];
        Dmo[k++] = 3*Dl[2]*l[2]*l[2]*l[3] + l[2]*l[2]*l[2]*Dl[3];
        Dmo[k++] = 3*Dl[3]*l[3]*l[3]*l[0] + l[3]*l[3]*l[3]*Dl[0];
        Dmo[k++] = 3*Dl[1]*l[1]*l[1]*l[0] + l[1]*l[1]*l[1]*Dl[0];
        Dmo[k++] = 3*Dl[0]*l[0]*l[0]*l[3] + l[0]*l[0]*l[0]*Dl[3];
        Dmo[k++] = 3*Dl[3]*l[3]*l[3]*l[2] + l[3]*l[3]*l[3]*Dl[2];
        Dmo[k++] = 3*Dl[2]*l[2]*l[2]*l[1] + l[2]*l[2]*l[2]*Dl[1];

        ffassert(k==28);
        for(int i=0; i<28;++i)
            Dv[i]=R3();
        for(int i=0; i<28;++i)
            for(int j=0; j<28;++j)
                Dv[i]+= C1[p[i]][j]*Dmo[j];
        RN_ f0x(val('.', 0, op_dx));
        RN_ f0y(val('.', 0, op_dy));
        RN_ f0z(val('.', 0, op_dz));
        for(int i=0; i<28;++i)
        {
            if(iDx) f0x[i] = Dv[i].x;
            if(iDy) f0y[i] = Dv[i].y;
            if(iDz) f0z[i] = Dv[i].z;
        }
        if (whatd & ( Fop_D2))
             ffassert(0); // do do ..
    }

}
 static TypeOfFE_P3pnc_3d P3pnc_3d;
 GTypeOfFE< Mesh3 > &Elm_P3pnc_3d(P3pnc_3d);

 static AddNewFE3 TFE_P3pnc_3d("P3pnc3d", &Elm_P3pnc_3d);

}    // namespace Fem2D

// --- fin --
