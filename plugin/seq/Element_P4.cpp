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
// Attention probleme de numerotation des inconnues
// -------------------------------------------------
// dans freefem, il y a un noeud par objets  sommet, arete, element.
// et donc la numerotation des dl dans l'element depend
// de l'orientation des aretes
//
/// ---------------------------------------------------------------
namespace Fem2D {
  // ------ P4  Hierarchical (just remove P1 node of the P2 finite element)  --------
  class TypeOfFE_P4Lagrange : public TypeOfFE {
   public:
    static const int k = 4;
    static const int ndf = (k + 2) * (k + 1) / 2;
    static int Data[];
    static double Pi_h_coef[];
    static const int nn[15][4];
    static const int aa[15][4];
    static const int ff[15];
    static const int il[15];
    static const int jl[15];
    static const int kl[15];

    TypeOfFE_P4Lagrange( ) : TypeOfFE(3 + 3 * 3 + 3, 1, Data, 4, 1, 15 + 6, 15, 0) {
      static const R2 Pt[15] = {R2(0 / 4., 0 / 4.), R2(4 / 4., 0 / 4.), R2(0 / 4., 4 / 4.),
                                R2(3 / 4., 1 / 4.), R2(2 / 4., 2 / 4.), R2(1 / 4., 3 / 4.),
                                R2(0 / 4., 3 / 4.), R2(0 / 4., 2 / 4.), R2(0 / 4., 1 / 4.),
                                R2(1 / 4., 0 / 4.), R2(2 / 4., 0 / 4.), R2(3 / 4., 0 / 4.),
                                R2(1 / 4., 2 / 4.), R2(2 / 4., 1 / 4.), R2(1 / 4., 1 / 4.)};

      // 3,4,5, 6,7,8, 9,10,11,
      int other[15] = {0, 1, 2, 5, 4, 3, 8, 7, 6, 11, 10, 9, 12, 13, 14};
      int kk = 0;

      for (int i = 0; i < NbDoF; i++) {
        pij_alpha[kk++] = IPJ(i, i, 0);
        if (other[i] != i) {
          pij_alpha[kk++] = IPJ(i, other[i], 0);
        }

        P_Pi_h[i] = Pt[i];
      }

      assert(P_Pi_h.N( ) == NbDoF);
      assert(pij_alpha.N( ) == kk);
    }

    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const RdHat &PHat,
            RNMK_ &val) const;
    void Pi_h_alpha(const baseFElement &K, KN_< double > &v) const {
      for (int i = 0; i < 15 + 6; ++i) {
        v[i] = 1;
      }

      int e0 = K.EdgeOrientation(0);
      int e1 = K.EdgeOrientation(1);
      int e2 = K.EdgeOrientation(2);
      int ooo[6] = {e0, e0, e1, e1, e2, e2};
      /*   3,4
       *   5,
       *   6,7
       *   8,9,
       *   10,
       *   11,12,
       *   13,14,
       *   15
       *   16,17
       */
      int iii[6] = {3, 6, 8, 11, 13, 16};
      int jjj[6] = {};

      for (int i = 0; i < 6; ++i) {
        jjj[i] = iii[i] + 1;    // si orient = -1
      }

      for (int i = 0; i < 6; ++i) {
        if (ooo[i] == 1) {
          v[jjj[i]] = 0;
        } else {
          v[iii[i]] = 0;
        }
      }
    }
  };

  // on what     nu df on node node of df
  int TypeOfFE_P4Lagrange::Data[] = {
    0, 1, 2, 3, 3, 3, 4, 4, 4, 5, 5,  5,  6,  6,  6,    // the support number  of the node of the df
    0, 0, 0, 0, 1, 2, 0, 1, 2, 0, 1,  2,  0,  1,  2,    // the number of the df on  the node
    0, 1, 2, 3, 3, 3, 4, 4, 4, 5, 5,  5,  6,  6,  6,    // the node of the df
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0,  0,    // the df come from which FE (generaly 0)
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,    // which are de df on sub FE
    0,    // for each compontant $j=0,N-1$ it give the sub FE associated
    0, 15};
  void TypeOfFE_P4Lagrange::FB(const bool *whatd, const Mesh &, const Triangle &K,
                               const RdHat &PHat, RNMK_ &val) const {
    R2 A(K[0]), B(K[1]), C(K[2]);
    R l0 = 1. - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
    R L[3] = {l0 * k, l1 * k, l2 * k};

    throwassert(val.N( ) >= 10);
    throwassert(val.M( ) == 1);
    // Attention il faut renumeroter les fonction de bases
    // car dans freefem++, il y a un node par sommet, arete or element
    // et la numerotation naturelle  mais 2 noud pas arete
    // donc p est la perumation
    // echange de numerotation si les arete sont dans le mauvais sens
    int p[15] = {};

    for (int i = 0; i < 15; ++i) {
      p[i] = i;
    }

    if (K.EdgeOrientation(0) < 0) {
      Exchange(p[3], p[5]);    // 3,4
    }

    if (K.EdgeOrientation(1) < 0) {
      Exchange(p[6], p[8]);    // 5,6
    }

    if (K.EdgeOrientation(2) < 0) {
      Exchange(p[9], p[11]);    // 7,8
    }

    val = 0;
    /*
     * //  les fonction de base du Pk Lagrange sont
     */

    if (whatd[op_id]) {
      RN_ f0(val('.', 0, op_id));

      for (int df = 0; df < ndf; df++) {
        int pdf = p[df];
        R f = 1. / ff[df];

        for (int i = 0; i < k; ++i) {
          f *= L[nn[df][i]] - aa[df][i];
        }

        f0[pdf] = f;
      }
    }

    if (whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] || whatd[op_dxy]) {
      R2 D[] = {K.H(0) * k, K.H(1) * k, K.H(2) * k};
      if (whatd[op_dx] || whatd[op_dy]) {
        for (int df = 0; df < ndf; df++) {
          int pdf = p[df];
          R fx = 0., fy = 0., f = 1. / ff[df];

          for (int i = 0; i < k; ++i) {
            int n = nn[df][i];
            R Ln = L[n] - aa[df][i];
            fx = fx * Ln + f * D[n].x;
            fy = fy * Ln + f * D[n].y;
            f = f * Ln;
          }

          if (whatd[op_dx]) {
            val(pdf, 0, op_dx) = fx;
          }

          if (whatd[op_dy]) {
            val(pdf, 0, op_dy) = fy;
          }
        }
      }

      if (whatd[op_dyy] || whatd[op_dxy] || whatd[op_dxx]) {
        for (int df = 0; df < ndf; df++) {
          int pdf = p[df];
          R fx = 0., fy = 0., f = 1. / ff[df];
          R fxx = 0., fyy = 0., fxy = 0.;

          for (int i = 0; i < k; ++i) {
            int n = nn[df][i];
            R Ln = L[n] - aa[df][i];
            fxx = fxx * Ln + 2. * fx * D[n].x;
            fyy = fyy * Ln + 2. * fy * D[n].y;
            fxy = fxy * Ln + fx * D[n].y + fy * D[n].x;
            fx = fx * Ln + f * D[n].x;
            fy = fy * Ln + f * D[n].y;
            f = f * Ln;
          }

          if (whatd[op_dxx]) {
            val(pdf, 0, op_dxx) = fxx;
          }

          if (whatd[op_dyy]) {
            val(pdf, 0, op_dyy) = fyy;
          }

          if (whatd[op_dxy]) {
            val(pdf, 0, op_dxy) = fxy;
          }
        }
      }
    }
  }


#include "Element_PkL.hpp" // for Pk_L
static TypeOfFE_Pk_L P4_L(4);
GTypeOfFE< MeshL > &Elm_P4_L(P4_L);


// Author: F. Hecht , P-H Tournier, Jet Hoe Tang jethoe.tang@googlemail.com
// Jan 2017
// in tets
class TypeOfFE_P4_3d : public GTypeOfFE< Mesh3 > {
public:
    typedef Mesh3 Mesh;
    typedef Mesh3::Element Element;
    typedef GFElement< Mesh3 > FElement;
    static const int kp = 4;    // P4
    static const int ndof = (kp + 3) * (kp + 2) * (kp + 1) / 6;
    static int dfon[];
    static int nl[ndof][kp];
    static int cl[ndof][kp];
    static int cp[ndof];
    static int pp[ndof][4];
    static const int d = Mesh::Rd::d;
    
    TypeOfFE_P4_3d( );    // constructor
    void FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K, const RdHat &PHat,
            RNMK_ &val) const;
    void set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M, int ocoef, int odf,
             int *nump) const;
    void OrientDoFOfFace(int odf,const Element &K,int f,int *pp,int inv) const
    {
        const int doff = odf+22+3*f; // first dof of face f
        int np=K.facePermutation(f);
        int i=np/2, j= np%2 ? 2:1;
        int i0=i,i1=(i+j)%3,i2=(i+j+j)%3;
        if( inv)
        {
            pp[doff+i0] = doff+0;
            pp[doff+i1] = doff+1;
            pp[doff+i2] = doff+2;
        }
        else
        {
            pp[doff+0] = doff+i0;
            pp[doff+1] = doff+i1;
            pp[doff+2] = doff+i2;

        }
/*        static int count = 0;
        if( count++ < 8)
        cout << " OrientDoFOfFace " << np << " " << inv << " f " << " " << f << " ::::" << i0 << i1 << i2 << " : " << pp[doff+0] << " " << pp[doff+1] << " " << pp[doff+2] << endl;*/
    }
 
};

int  TypeOfFE_P4_3d::pp[35][4] = {
  {4, 0, 0, 0},
   {0, 4, 0, 0},
   {0, 0, 4, 0},
   {0, 0, 0, 4},
    
   {1, 3, 0, 0},
   {2, 2, 0, 0},
   {3, 1, 0, 0},
    
   {1, 0, 3, 0},
   {2, 0, 2, 0},
   {3, 0, 1, 0},
    
   {1, 0, 0, 3},
   {2, 0, 0, 2},
   {3, 0, 0, 1},
    
   {0, 1, 3, 0},
   {0, 2, 2, 0},
   {0, 3, 1, 0},
    
   {0, 1, 0, 3},
   {0, 2, 0, 2},
   {0, 3, 0, 1},
    
   {0, 0, 1, 3},
   {0, 0, 2, 2},
   {0, 0, 3, 1},
    
   {0, 2, 1, 1},
   {0, 1, 2, 1},
   {0, 1, 1, 2},
    
   {1, 0, 1, 2},
   {1, 0, 2, 1},
   {2, 0, 1, 1},
    
   {2, 1, 0, 1},
   {1, 2, 0, 1},
  {1, 1, 0, 2},
    
   {1, 1, 2, 0},
   {1, 2, 1, 0},
   {2, 1, 1, 0},
    
   {1, 1, 1, 1}};



int  TypeOfFE_P4_3d::nl[35][4] = {
  {0, 0, 0, 0},
   {1, 1, 1, 1},
   {2, 2, 2, 2},
   {3, 3, 3, 3},
   {0, 1, 1, 1},
   {0, 0, 1, 1},
   {0, 0, 0, 1},
   {0, 2, 2, 2},
   {0, 0, 2, 2},
   {0, 0, 0, 2},

  {0, 3, 3, 3},
   {0, 0, 3, 3},
   {0, 0, 0, 3},
   {1, 2, 2, 2},
   {1, 1, 2, 2},
   {1, 1, 1, 2},
   {1, 3, 3, 3},
   {1, 1, 3, 3},
   {1, 1, 1, 3},
   {2, 3, 3, 3},

   {2, 2, 3, 3},
   {2, 2, 2, 3},
   {1, 1, 2, 3},
   {1, 2, 2, 3},
   {1, 2, 3, 3},
   {0, 2, 3, 3},
   {0, 2, 2, 3},
   {0, 0, 2, 3},
   {0, 0, 1, 3},
   {0, 1, 1, 3},

   {0, 1, 3, 3},
   {0, 1, 2, 2},
   {0, 1, 1, 2},
   {0, 0, 1, 2},
   {0, 1, 2, 3}};



int  TypeOfFE_P4_3d::cl[35][4] = {
  {0, 1, 2, 3},
   {0, 1, 2, 3},
   {0, 1, 2, 3},
   {0, 1, 2, 3},
   {0, 0, 1, 2},
   {0, 1, 0, 1},
   {0, 1, 2, 0},
   {0, 0, 1, 2},
   {0, 1, 0, 1},
   {0, 1, 2, 0},

  {0, 0, 1, 2},
   {0, 1, 0, 1},
   {0, 1, 2, 0},
   {0, 0, 1, 2},
   {0, 1, 0, 1},
   {0, 1, 2, 0},
   {0, 0, 1, 2},
   {0, 1, 0, 1},
   {0, 1, 2, 0},
   {0, 0, 1, 2},

  {0, 1, 0, 1},
   {0, 1, 2, 0},
   {0, 1, 0, 0},
   {0, 0, 1, 0},
   {0, 0, 0, 1},
   {0, 0, 0, 1},
   {0, 0, 1, 0},
   {0, 1, 0, 0},
   {0, 1, 0, 0},
   {0, 0, 1, 0},

  {0, 0, 0, 1},
   {0, 0, 0, 1},
   {0, 0, 1, 0},
   {0, 1, 0, 0},
   {0, 0, 0, 0}};



int  TypeOfFE_P4_3d::cp[35] = {
24, 24, 24, 24, 6, 4, 6, 6, 4, 6, 6, 4, 6, 6, 4, 6, 6, 4, 6, 6, 4, 6, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1};
int TypeOfFE_P4_3d::dfon[] = {1, 3, 3, 1};    // 2 dofs on each edge, 2 dofs on each face

TypeOfFE_P4_3d::TypeOfFE_P4_3d( ) : GTypeOfFE< Mesh >(TypeOfFE_P4_3d::dfon, 1, 3, false, false)
{
    typedef Element E;
    int n = this->NbDoF;
    bool dd = verbosity > 5;
    if (dd) {
        cout << "\n +++ P4  : ndof : " << n << " " << this->PtInterpolation.N( ) << endl;
    }
    
    R3 *Pt = this->PtInterpolation;
    // construction of interpolation ppoint
    
    {
        double cc = 1. / 4.;
        
        for (int i = 0; i < ndof; ++i) {
            Pt[i] = R3::KHat[0] * cc * pp[i][0] + R3::KHat[1] * cc * pp[i][1] +
            R3::KHat[2] * cc * pp[i][2] + R3::KHat[3] * cc * pp[i][3];
        }
        
        if (dd) {
            cout << this->PtInterpolation << endl;
        }
    }
    
    for (int i = 0; i < n; i++) {
        this->pInterpolation[i] = i;
        this->cInterpolation[i] = 0;
        this->dofInterpolation[i] = i;
        this->coefInterpolation[i] = 1.;
    }
}

void TypeOfFE_P4_3d::set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M,
                         int ocoef, int odf, int *nump) const {
    //  faux nump don la numerotation des p -> local
    // ne marche que le cas scalaire ??? FH...
    /*
     for (int k=0;k<M.ncoef;++k)
     vdf[M.dofe[k]] += M.coef[k]*vpt(M.p[k],M.comp[k]);
     */
    int Np=M.p.N();
    int n = this->NbDoF;
    int *p = M.p;
    
    //  for (int i = 0; i < n; ++i) {
    //    M.p[i] = i;
    //  }
    
    if (verbosity > 9) {
        cout << " P4  set:" << odf << " : ";
    }
    
    int dof = 4+odf;
    
    for (int e = 0; e < 6; ++e) {
        int oe = K.EdgeOrientation(e);
        int i1=dof;
        int i2=dof+2;
        //cout << e <<"  "<< i1 << " " << i2 << "  , " << oe <<"  p: " <<  p[i1] << " " <<  p[i2] <<" " << Np << endl;
        ffassert( i1>=0 && i2 >=0);
        ffassert( i1<Np && i2 <Np );
        
        if ((oe < 0) && (p[i1] < p[i2])) {
            swap(p[i1], p[i2]);
        } else if ((oe > 0) && (p[i1] > p[i2])) {
            swap(p[i1], p[i2]);
        }
        dof += 3;
    }
    for (int f = 0; f < 4; ++f)
        OrientDoFOfFace(odf,K,f,p,1);

    if (verbosity > 99) {
        cout << " " << M.p << endl;  ;
    }
    
}

void TypeOfFE_P4_3d::FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K,
                        const RdHat &PHat, RNMK_ &val) const {
    assert(val.N( ) >= 35);    // 23 degrees of freedom
    assert(val.M( ) == 1);     // 3 components
    // int n = this->NbDoF;
    // -------------
    // perm: the permutation for which the 4 tetrahedron vertices are listed with increasing GLOBAL
    // number (i.e. perm[0] is the local number of the vertex with the smallest global number, ...
    // perm[3] is the local number of the vertex with the biggest global number.)
    // -------------
    R ld[4];
    PHat.toBary(ld);
    ld[0] *= 4.;
    ld[1] *= 4.;
    ld[2] *= 4.;
    ld[3] *= 4.;
    
    int p[35];
    for(int i=0; i<35;++i)
        p[i]=i;
    
    {
        int dof = 4;
        
        for (int e = 0; e < 6; ++e) {
            int oe = K.EdgeOrientation(e);
            if (oe < 0) {
                swap(p[dof], p[dof + 2]);
            }
            
            dof += 3;
        }
        for (int f = 0; f < 4; ++f)
          OrientDoFOfFace(0,K,f,p,0);

    }
    //    static int ddd = 100;
    //    ddd++;
    val = 0.;
    RN_ f0(val('.', 0, op_id));
    //  if (ddd < 20) {
    //     cout << ld[0] << " " << ld[1] << " " << ld[2] << " " << ld[3] << " ::";
    //  }
    
    if (whatd & Fop_D0) {
        for (int i = 0; i < 35; ++i) {
            R fi = 1. / cp[i];
            
            for (int l = 0; l < 4; ++l) {
                fi *= ld[nl[i][l]] - cl[i][l];
            }
            
            //      if (ddd < 20) {
            //       cout << " " << fi;
            //    }
            
            f0[p[i]] = fi;
        }
        
        ///      if (ddd < 20) {
        ///       cout << endl;
        ///     }
    }
    
    if (whatd & (Fop_D1 | Fop_D2)) {
        R3 Dld[4], Df[35];
        K.Gradlambda(Dld);
        Dld[0] *= 4.;
        Dld[1] *= 4.;
        Dld[2] *= 4.;
        Dld[3] *= 4.;
        
        for (int i = 0; i < 35; ++i) {
            R fi = 1. / cp[i];
            R3 &dfi = Df[p[i]];
            
            for (int l = 0; l < 4; ++l) {
                double ci = ld[nl[i][l]] - cl[i][l];
                dfi *= ci;
                dfi += fi * Dld[nl[i][l]];
                fi *= ci;
            }
            
            RN_ f0x(val('.', 0, op_dx));
            RN_ f0y(val('.', 0, op_dy));
            RN_ f0z(val('.', 0, op_dz));
            if (whatd & Fop_dx) {
                for (int i = 0; i < 35; ++i) {
                    f0x[i] = Df[i].x;
                }
            }
            
            if (whatd & Fop_dy) {
                for (int i = 0; i < 35; ++i) {
                    f0y[i] = Df[i].y;
                }
            }
            
            if (whatd & Fop_dz) {
                for (int i = 0; i < 35; ++i) {
                    f0z[i] = Df[i].z;
                }
            }
            
            ffassert(!(whatd & Fop_D2));    // no D2 to do !!!
        }
    }
}

class TypeOfFE_P4_S : public GTypeOfFE< MeshS > {
public:
    typedef MeshS Mesh;
    typedef MeshS::Element Element;
    typedef GFElement< MeshS > FElement;
    typedef Mesh::RdHat RdHat;
    typedef Mesh::Rd Rd;
    
    static const int kp = 4;    // P4
    static const int ndof =  (kp + 2) * (kp + 1) / 2;//  15
    constexpr static int dfon[]= {1, 3, 3, 0};
    
    
    
    static const int d = Rd::d;
    static const int nn[15][4];
    static const int aa[15][4];
    static const int ff[15];
    static const int il[15];
    static const int jl[15];
    static const int kl[15];

    
    constexpr static const int dHat = RdHat::d;
    
    TypeOfFE_P4_S( );    // constructor
    void FB(const What_d whatd, const MeshS &Th, const MeshS::Element &K, const RdHat &PHat,
            RNMK_ &val) const;
    void set(const MeshS &Th, const Element &K, InterpolationMatrix< RdHat > &M, int ocoef, int odf,
             int *nump) const;
};
TypeOfFE_P4_S::TypeOfFE_P4_S( ) : GTypeOfFE< MeshS >(TypeOfFE_P4_3d::dfon, 1, 3, false, false)
{
 
    
    typedef Element E;
    int n = this->NbDoF;
    bool dd = verbosity > 5;
    if (dd) {
        cout << "\n +++ P4  : ndof : " << n << " " << this->PtInterpolation.N( ) << endl;
    }
    
    RdHat *Pt = this->PtInterpolation;
    // construction of interpolation ppoint
    
    {
        double cc = 1. / 4.;
        
        for (int i = 0; i < ndof; ++i)
        Pt[i] = RdHat::KHat[0] * cc * il[i]+ RdHat::KHat[1] * cc * jl[i] + RdHat::KHat[2] * cc * kl[i] ;
        
        
        if (dd) {
            cout << this->PtInterpolation << endl;
        }
    }
    
    for (int i = 0; i < n; i++) {
        this->pInterpolation[i] = i;
        this->cInterpolation[i] = 0;
        this->dofInterpolation[i] = i;
        this->coefInterpolation[i] = 1.;
    }
}

void TypeOfFE_P4_S::set(const MeshS &Th, const TypeOfFE_P4_S::Element &K, InterpolationMatrix< TypeOfFE_P4_S::RdHat > &M,
                        int ocoef, int odf, int *nump) const {
    //  faux nump don la numerotation des p -> local
    // ne marche que le cas scalaire ??? FH...
    /*
     for (int k=0;k<M.ncoef;++k)
     vdf[M.dofe[k]] += M.coef[k]*vpt(M.p[k],M.comp[k]);
     */
    int Np=M.p.N();
    int n = this->NbDoF;
    int *p = M.p;
    
    if (verbosity > 9) {
        cout << " P4_S  set:" << odf << " : ";
    }
    
    int dof = 3+odf;
    
    for (int e = 0; e < 3; ++e) {
        int oe = K.EdgeOrientation(e);
        int i1=dof;
        int i2=dof+2;
        //cout << e <<"  "<< i1 << " " << i2 << "  , " << oe <<"  p: " <<  p[i1] << " " <<  p[i2] <<" " << Np << endl;
        ffassert( i1>=0 && i2 >=0);
        ffassert( i1<Np && i2 <Np );
        
        if ((oe < 0) && (p[i1] < p[i2])) {
            swap(p[i1], p[i2]);
        } else if ((oe > 0) && (p[i1] > p[i2])) {
            swap(p[i1], p[i2]);
        }
        
        dof += 3;
    }
    if (verbosity > 99) {
        cout << " " << M.p << endl;  ;
    }
    
}

void TypeOfFE_P4_S::FB(const What_d whatd, const TypeOfFE_P4_S::Mesh &Th, const TypeOfFE_P4_S::Element &K,
                       const TypeOfFE_P4_S::RdHat &PHat, RNMK_ &val) const {
    
      
    
    assert(val.N( ) >= 15);    // 10 degrees of freedom
    assert(val.M( ) == 1);     // 1 components
    // int n = this->NbDoF;
    // -------------
    // perm: the permutation for which the 4 tetrahedron vertices are listed with increasing GLOBAL
    // number (i.e. perm[0] is the local number of the vertex with the smallest global number, ...
    // perm[3] is the local number of the vertex with the biggest global number.)
    // -------------
    const int ndf=15;
    R L[3];
    PHat.toBary(L);
    L[0] *= 4.;
    L[1] *= 4.;
    L[2] *= 4.;
    
    
    throwassert(val.N( ) >= 15);
    throwassert(val.M( ) == 1);
    // Attention il faut renumeroter les fonction de bases
    // car dans freefem++, il y a un node par sommet, arete or element
    // et la numerotation naturelle  mais 2 noud pas arete
    // donc p est la perumation
    // echange de numerotation si les arete sont dans le mauvais sens
    int p[15] = {};
    
    for (int i = 0; i < 15; ++i) {
        p[i] = i;
    }
    
    if (K.EdgeOrientation(0) < 0) {
        Exchange(p[3], p[5]);    // 3,5
    }
    
    if (K.EdgeOrientation(1) < 0) {
        Exchange(p[6], p[8]);    // 6,8
    }
    
    if (K.EdgeOrientation(2) < 0) {
        Exchange(p[9], p[11]);    // 9,11
    }
    
    val = 0;
    
    if (whatd & Fop_D0) {
        RN_ f0(val('.', 0, op_id));
        
        for (int df = 0; df < ndf; df++) {
            int pdf = p[df];
            R f = 1. / ff[df];
            
            for (int i = 0; i < kp; ++i) {
                f *= L[nn[df][i]] - aa[df][i];
            }
            
            f0[pdf] = f;
        }
    }
    
    if(whatd & (Fop_D1|Fop_D2))
    {
        
        R3 D[3] ;
        K.Gradlambda(D);
        D[0]*= kp;
        D[1]*= kp;
        D[2]*= kp;
        if (whatd & (Fop_D1|Fop_D2)) {
            for (int df = 0; df < ndf; df++) {
                int pdf = p[df];
                R fx = 0., fy = 0.,fz = 0., f = 1. / ff[df];
                
                for (int i = 0; i < kp; ++i) {
                    int n = nn[df][i];
                    R Ln = L[n] - aa[df][i];
                    fx = fx * Ln + f * D[n].x;
                    fy = fy * Ln + f * D[n].y;
                    fz = fz * Ln + f * D[n].z;
                    f = f * Ln;
                }
                
                if (whatd & Fop_dx) {
                    val(pdf, 0, op_dx) = fx;
                }
                
                if (whatd & Fop_dy) {
                    val(pdf, 0, op_dy) = fy;
                }
                if (whatd & Fop_dz) {
                    val(pdf, 0, op_dz) = fz;
                }
            }
        }
        
        if (whatd &Fop_D2) {
            for (int df = 0; df < ndf; df++) {
                int pdf = p[df];
                R fx = 0., fy = 0.,fz=0., f = 1. / ff[df];
                R fxx = 0., fyy = 0., fzz=0.,  fxy = 0., fxz = 0., fyz = 0. ;
                
                for (int i = 0; i < kp; ++i) {
                    int n = nn[df][i];
                    R Ln = L[n] - aa[df][i];
                    fxx = fxx * Ln + 2. * fx * D[n].x;
                    fyy = fyy * Ln + 2. * fy * D[n].y;
                    fzz = fzz * Ln + 2. * fz * D[n].z;
                    fxy = fxy * Ln + fx * D[n].y + fy * D[n].x;
                    fxz = fxz * Ln + fx * D[n].z + fz * D[n].x;
                    fyz = fyz * Ln + fy * D[n].z + fz * D[n].y;
                    fx = fx * Ln + f * D[n].x;
                    fy = fy * Ln + f * D[n].y;
                    fz = fz * Ln + f * D[n].z;
                    f = f * Ln;
                }
                
                if (whatd & Fop_dxx) {
                    val(pdf, 0, op_dxx) = fxx;
                }
                
                if (whatd & Fop_dyy) {
                    val(pdf, 0, op_dyy) = fyy;
                }
                if (whatd & Fop_dzz) {
                    val(pdf, 0, op_dzz) = fzz;
                }
                
                if (whatd & Fop_dxy) {
                    val(pdf, 0, op_dxy) = fxy;
                }
                if (whatd & Fop_dxz) {
                    val(pdf, 0, op_dxz) = fxz;
                }
                if (whatd & Fop_dyz) {
                    val(pdf, 0, op_dyz) = fyz;
                }
            }
        }
    }
}
#include "Element_P4.hpp"


// link with FreeFem++
static TypeOfFE_P4Lagrange P4LagrangeP4;
// a static variable to add the finite element to freefem++
//  static AddNewFE P4Lagrange("P4", &P4LagrangeP4);
static TypeOfFE_P4_3d P4_3d;
static TypeOfFE_P4_S P4_S;

GTypeOfFE< Mesh3 > &Elm_P4_3d(P4_3d);
GTypeOfFE< MeshS > &Elm_P4_S(P4_S);



static void init( ) {
    if(verbosity && mpirank ==0 )
    cout << " load : P4 ";
  AddNewFE("P4", &P4LagrangeP4);
  static ListOfTFE FE_P4("P4", &P4LagrangeP4); // to add P4 in list of Common FE

  AddNewFE3("P43d", &Elm_P4_3d,"P4");
  AddNewFES("P4S", &Elm_P4_S,"P4");
  AddNewFEL("P4L", &Elm_P4_L,"P4");
}

}    // namespace Fem2D
LOADFUNC(Fem2D::init);
// --- fin --
