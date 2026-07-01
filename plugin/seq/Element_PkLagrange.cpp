#include "ff++.hpp"
#include "AddNewFE.h"
#include "polynomial_PKLagrange.hpp"
#include "utils_PKLagrange.hpp"
#include "Element_PkLagrange3d.hpp"

//This file is strongly inspired by the files elementsP3 and p4.
using namespace std;
// Attention probleme de numerotation des inconnues
// -------------------------------------------------
// dans freefem, il y a un noeud par objets  sommet, arete, element.
// et donc la numerotation des dl dans l'element depend
// de l'orientation des aretes
//
/// ---------------------------------------------------------------
namespace Fem2D {

////////////////////////////////////////////////////////////////////////
///////////////////////////////  2D  ///////////////////////////////////
////////////////////////////////////////////////////////////////////////

// ------ PK  Hierarchical  --------
class TypeOfFE_PkLagrange : public TypeOfFE {

  private:
    // static int * pdata;
    // FUnction to Fill Data table
    /*int *PrepareData(int PK) {
        int ndof = (PK + 2) * (PK + 1) / 2;
        //static vector<int> data(5 * ndof + 3, 0);
        // Fill data
        int *Data= FillDataLagrange(PK);
        //return data.data();
        return Data;
    }*/

  public:
    int k; // The degree
    int ndf; // Number of DOfs
    vector<double> Pi_h_coef;
    // Polynomial parameters(Basis fcts):
    //nn for the barycentric function, aa for the shifts, and ff for the denominator
    vector<vector<long>> nn; //[ndf][k];
    vector<vector<long>> aa; //[ndf][k];
    vector<long> ff;         //[ndf];
    // Points coordinates
    vector<long> il; //[ndf];
    vector<long> jl; //[ndf];
    vector<long> kl; //[ndf];
    // Constuctor
    TypeOfFE_PkLagrange(int PK)
        : TypeOfFE((PK + 1) * (PK + 2) / 2, 1, FillDataLagrange(PK), PK, PK,
                   (PK + 1) * (PK + 2) / 2 + ((PK % 2) ? (PK - 1) * 3 : (PK - 2) * 3), (PK + 1) * (PK + 2) / 2, 0) {
        k = PK;
        ndf = (k + 2) * (k + 1) / 2;
        nn.resize(ndf);
        aa.resize(ndf);
        for (long i = 0; i < ndf; i++) {
            aa[i].resize(k, -1);
            nn[i].resize(k, -1);
        }
        ff.resize(ndf, -1);
        il.resize(ndf, -1);
        jl.resize(ndf, -1);
        kl.resize(ndf, -1);
        // Fill Data Lagrange
        // FillDataLagrange( k, Data, Pi_h_coef );
        Pi_h_coef.resize(ndf, 1.);

        // Fill ff,il,jl,kl,nn,aa
        BasisFctPK(k, nn, aa, ff, il, jl, kl);
        vector<R2> Pt = PtConstruction(k);
        // 
        vector<int> other(ndf, 0); //[ndf] = {-1, -1, -1, 4, 3, 6, 5, 8, 7, -1}; in P3//
        FillOther(other, k, ndf);

        int kk = 0;

        for (int i = 0; i < NbDoF; i++) {
            pij_alpha[kk++] = IPJ(i, i, 0);
            if (other[i] != i) {
                pij_alpha[kk++] = IPJ(i, other[i], 0);
            }

            P_Pi_h[i] = Pt[i];
        }
        assert(P_Pi_h.N() == NbDoF);
        assert(pij_alpha.N() == kk);
    }

    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const RdHat &PHat, RNMK_ &val) const;
    void Pi_h_alpha(const baseFElement &K, KN_<double> &v) const {

        // Generic code
        int nbrPerm = (k % 2) ? (k - 1) * 3 : (k - 2) * 3;         // Number of permutation
        vector<vector<int>> indicesIJ(2, vector<int>(nbrPerm, 0)); // IJ indices
        vector<int> ooo(nbrPerm, 0);
        int KK = (k % 2) ? (k - 1) : (k - 2);
        for (int i = 0; i < nbrPerm; i++) {
            if (k % 2 == 1)
                indicesIJ[0][i] = (3 + 2 * i);
            else {
                int arete = i / (k - 2);
                int idxlocal = i % (k - 2);
                if (k - 2 > 2)
                    idxlocal /= ((k - 2) / 2);
                indicesIJ[0][i] = (3 + 2 * i + arete + idxlocal);
            }
            indicesIJ[1][i] = (1 + indicesIJ[0][i]);
            ooo[i] = K.EdgeOrientation(i / KK);
        }
        for (int i = 0; i < ndf + nbrPerm; ++i) {
            v[i] = 1;
        }
        for (int i = 0; i < nbrPerm; ++i) {
            if (ooo[i] == 1) {
                v[indicesIJ[1][i]] = 0;
            } else {
                v[indicesIJ[0][i]] = 0;
            }
        }
    }
};

void TypeOfFE_PkLagrange::FB(const bool *whatd, const Mesh &, const Triangle &K, const RdHat &PHat, RNMK_ &val) const {
    R2 A(K[0]), B(K[1]), C(K[2]);
    R l0 = 1 - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
    R L[3] = {l0 * k, l1 * k, l2 * k};

    /*throwassert(val.N( ) >= 10);
    throwassert(val.M( ) == 1);*/
    vector<int> p(ndf, 0);

    for (int i = 0; i < ndf; ++i) {
        p[i] = i;
    }

    // Get index exchange
    vector<pair<int, int>> ExId = ExchangeidxVector(k);

    int nbrPermAxis = (k % 2) ? (k - 1) : (k - 2); // Number of permutation per axis
    for (int i = 0; i < 3; i++) {
        if (K.EdgeOrientation(i) < 0) {
            for (int j = i * nbrPermAxis / 2; j < (i + 1) * nbrPermAxis / 2; j++) {
                Exchange(p[ExId[j].first], p[ExId[j].second]);
            }
        }
    }
    val = 0;

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

////////////////////////////////////////////////////////////////////////
///////////////////////////////  3D  ///////////////////////////////////
////////////////////////////////////////////////////////////////////////

class TypeOfFE_PKLagrange_3d : public GTypeOfFE< Mesh3 > {
   public:
    typedef Mesh3 Mesh;
    typedef Mesh3::Element Element;
    typedef GFElement< Mesh3 > FElement;
    int kp;                        // PK
    int ndof;                      //= (kp + 3) * (kp + 2) * (kp + 1) / 6;
    int dfon[4];                   // Topological info
    vector< vector< int > > nl;    //[20][3] for P3; // returns the lambdas used for polynomial
    vector< vector< int > > cl;    //[20][3]for P3;  //  polynomial form = prod 1/cp (nl-cl)
    vector< int > cp;              //[20]for P3; //  polynomial factors
    vector< vector< int > > pp;    //[20][4] for P3;  //bary. coordinates
    static const int d = Mesh::Rd::d;
    TypeOfFE_PKLagrange_3d(int PK);    // constructor
    void FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K, const RdHat &PHat, RNMK_ &val) const;
    void set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M, int ocoef, int odf, int *nump) const;
    void OrientDoFOfFace(int odf, const Element &K, int f, int *pp, int inv) const;

   private:
    static int *Preparedfon(int PK) {
      // Fct needed in the parent class constructor
      static int locdfon[4];
      // Fill data
      locdfon[0] = 1;                                     // 1 dof per vertice
      locdfon[1] = PK - 1;                                // dof per edge
      locdfon[2] = (PK - 1) * (PK - 2) / 2;               // dof on the surface
      locdfon[3] = (PK - 1) * (PK - 2) * (PK - 3) / 6;    // dof in the volume
      return locdfon;
    }
};

TypeOfFE_PKLagrange_3d::TypeOfFE_PKLagrange_3d(int PK) : GTypeOfFE< Mesh >(Preparedfon(PK), 1, PK, false, false) {
    kp = PK;                                      // order
    ndof = (PK + 1) * (PK + 2) * (PK + 3) / 6;    // # dofs
    nl.resize(ndof, vector< int >(PK));           // nl cl and cp desribe the shape fcts
    cl.resize(ndof, vector< int >(PK));
    cp.resize(ndof);

    pp.resize(ndof, vector< int >(4));    // barycoordinates

    // Filling dfon explicitly
    dfon[0] = 1;                                     // 1 dof per vertice
    dfon[1] = PK - 1;                                // dof per edge
    dfon[2] = (PK - 1) * (PK - 2) / 2;               // dof on the surface
    dfon[3] = (PK - 1) * (PK - 2) * (PK - 3) / 6;    // dof in the volume*/

    // Fill nl cl cp pp
    BasisFctPK3D(PK, nl, cl, cp, pp);

    // Fill pp
    // pp = generate_barycoordinates(PK);

    typedef Element E;
    int n = this->NbDoF;
    bool dd = verbosity > 5;
    if (dd) {
      cout << "\n +++ PK  : ndof : " << n << " " << this->PtInterpolation.N( ) << endl;
    }

    R3 *Pt = this->PtInterpolation;
    // construction of interpolation ppoint

    {
      double cc = 1. / PK;

      for (int i = 0; i < ndof; ++i) {
        Pt[i] = R3::KHat[0] * cc * pp[i][0] + R3::KHat[1] * cc * pp[i][1] + R3::KHat[2] * cc * pp[i][2] + R3::KHat[3] * cc * pp[i][3];
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

// Adapted only to P4 de REF
void TypeOfFE_PKLagrange_3d::OrientDoFOfFace(int odf, const Element &K, int f, int *pp, int inv) const {
    const int doff = odf + 4 * dfon[0] + 6 * dfon[1] + f * dfon[2];    // first dof of face f 4 vertices + 6 edges+ x* number face

    int np = K.facePermutation(f);    // Returns first vertice of the face and the orientation
    int i = np / 2, j = np % (2) ? 2 : 1;
    int i0 = i, i1 = (i + j) % 3, i2 = (i + j + j) % 3;
    // position of the internel vertices of the faces
    int spos[3] = {0, kp - 3, dfon[2] - 1};
    if (inv) {
      std::vector< int > permut(3, 0);    //={i0,i1,i2};
      permut[i0] = 0;
      permut[i1] = 1;
      permut[i2] = 2;

      // cout<<"I's "<<i0<<i1<<i2<<endl;
      //  The internal dofs of the face
      //  can be seen as a (kp-3) 2d element
      //  Dofs are numberer diagonally
      //  The loops handles the dofs to be permuted
      //  using basy coordinates and the face orientation

      for (int i = 0; i <= kp - 3; i++) {
        for (int j = 0; j <= kp - 3 - i; j++) {
          int k = (kp - 3 - j - i);
          std::vector< int > barcoord = {i, j, k};
          int oldnumdof = numerotationLocIntern(kp, barcoord[0], barcoord[1], barcoord[2]);
          int newnumdof = numerotationPermutation(kp, barcoord, permut);
          // cout<<" doff perm: "<<oldnumdof <<" | "<<newnumdof<<endl;
          pp[doff + newnumdof] = doff + oldnumdof;
        }
      }
    }
    // This was used for the old version of the P4 (We don't need it)
    /*else {
        pp[doff + spos[0]] = doff + spos[i0];
        pp[doff + spos[1]] = doff + spos[i1];
        pp[doff + spos[2]] = doff + spos[i2];
        if(dfon[2]>3){
        pp[doff + apos[0][1]] = doff + apos[i0][i1];
        pp[doff + apos[0][2]] = doff + apos[i0][i2];
        pp[doff + apos[1][2]] = doff + apos[i1][i2];
        }
    }*/
}

void TypeOfFE_PKLagrange_3d::set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M, int ocoef, int odf, int *nump) const {
    //  faux nump don la numerotation des p -> local
    // ne marche que le cas scalaire ??? FH...
    /*
     for (int k=0;k<M.ncoef;++k)
     vdf[M.dofe[k]] += M.coef[k]*vpt(M.p[k],M.comp[k]);
     */
    int Np = M.p.N( );
    int n = this->NbDoF;
    int *p = M.p;

    for (int i = 0; i < n; ++i) {
      M.p[i] = i;
    }

    if (verbosity > 9) {
      cout << " PK  set:" << odf << " : ";
    }

    // int dof = 4+odf;
    int dof = odf;

    int nbrPermAxis = (kp % 2) ? (kp - 1) : (kp - 2);             // Number of permutation per edge
    vector< pair< int, int > > ExId = ExchangeidxVector3D(kp);    // indices to exchange

    for (int e = 0; e < 6; ++e) {
      int oe = K.EdgeOrientation(e);
      for (int j = e * nbrPermAxis / 2; j < (e + 1) * nbrPermAxis / 2; j++) {
        int i1 = dof + ExId[j].first;
        int i2 = dof + ExId[j].second;
        // cout << e <<"  "<< i1 << " " << i2 << "  , " << oe <<"  p: " <<  p[i1] << " " <<  p[i2] <<" " << Np << endl;
        ffassert(i1 >= 0 && i2 >= 0);
        ffassert(i1 < Np && i2 < Np);

        if ((oe < 0) && (p[i1] < p[i2])) {
          swap(p[i1], p[i2]);
        } else if ((oe > 0) && (p[i1] > p[i2])) {
          swap(p[i1], p[i2]);
        }
      }
    }

    if (dfon[2] > 1) {
      for (int f = 0; f < 4; ++f) OrientDoFOfFace(odf, K, f, p, 1);
    }
    if (verbosity > 99) {
      cout << " " << M.p << endl;
      ;
    }
}

void TypeOfFE_PKLagrange_3d::FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K, const RdHat &PHat, RNMK_ &val) const {
    assert(val.N( ) >= ndof);    //  degrees of freedom
    assert(val.M( ) == 1);       // 3 components
    // int n = this->NbDoF;
    // -------------
    // perm: the permutation for which the 4 tetrahedron vertices are listed with increasing GLOBAL
    // number (i.e. perm[0] is the local number of the vertex with the smallest global number, ...
    // perm[3] is the local number of the vertex with the biggest global number.)
    // -------------
    R ld[4];
    PHat.toBary(ld);
    ld[0] *= R(kp);
    ld[1] *= R(kp);
    ld[2] *= R(kp);
    ld[3] *= R(kp);

    vector< int > p(ndof, 0);
    for (int cpt = 0; cpt < p.size( ); cpt++) p[cpt] = cpt;
    {
      // int dof = 4;
      int dof = 0;

      vector< pair< int, int > > ExId = ExchangeidxVector3D(kp);
      int nbrPermAxis = (kp % 2) ? (kp - 1) : (kp - 2);    // Number of permutation per edge
      for (int e = 0; e < 6; ++e) {
        int oe = K.EdgeOrientation(e);
        if (oe < 0) {
          for (int j = e * nbrPermAxis / 2; j < (e + 1) * nbrPermAxis / 2; j++) swap(p[ExId[j].first], p[ExId[j].second]);
        }

        dof += kp - 1;
      }

      if (dfon[2] > 1) {
        for (int f = 0; f < 4; ++f) OrientDoFOfFace(0, K, f, p.data( ), 1);
      }
    }
    //    static int ddd = 100;
    //    ddd++;
    val = 0.;
    RN_ f0(val('.', 0, op_id));
    //  if (ddd < 20) {
    //     cout << ld[0] << " " << ld[1] << " " << ld[2] << " " << ld[3] << " ::";
    //  }

    if (whatd & Fop_D0) {
      for (int i = 0; i < ndof; ++i) {
        R fi = 1. / cp[i];

        for (int l = 0; l < kp; ++l) {
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
      R3 Dld[4], Df[ndof];
      K.Gradlambda(Dld);
      Dld[0] *= R(kp);
      Dld[1] *= R(kp);
      Dld[2] *= R(kp);
      Dld[3] *= R(kp);

      for (int i = 0; i < ndof; ++i) {
        R fi = 1. / cp[i];
        R3 &dfi = Df[p[i]];

        for (int l = 0; l < kp; ++l) {
          double ci = ld[nl[i][l]] - cl[i][l];
          dfi *= ci;
          dfi += fi * Dld[nl[i][l]];
          fi *= ci;
        }

        RN_ f0x(val('.', 0, op_dx));
        RN_ f0y(val('.', 0, op_dy));
        RN_ f0z(val('.', 0, op_dz));
        if (whatd & Fop_dx) {
          for (int i = 0; i < ndof; ++i) {
            f0x[i] = Df[i].x;
          }
        }

        if (whatd & Fop_dy) {
          for (int i = 0; i < ndof; ++i) {
            f0y[i] = Df[i].y;
          }
        }

        if (whatd & Fop_dz) {
          for (int i = 0; i < ndof; ++i) {
            f0z[i] = Df[i].z;
          }
        }

        ffassert(!(whatd & Fop_D2));    // no D2 to do !!!
      }
    }
}

static TypeOfFE_PkLagrange PKLagrange5(5);
static TypeOfFE_PkLagrange PKLagrange6(6);
static TypeOfFE_PkLagrange PKLagrange7(7);
static TypeOfFE_PkLagrange PKLagrange8(8);
static TypeOfFE_PkLagrange PKLagrange9(9);
static TypeOfFE_PkLagrange PKLagrange10(10);
static TypeOfFE_PkLagrange PKLagrange11(11);

#include "Element_PkL.hpp"
  // if you need to add a newFE (higher order $ORDER for instance)
  // you can simply do it by following the steps:
  /*
    1) Declare static TypeOfFE_PKLagrange_3d PKLagrangeP$ORDER($ORDER);
    2) GTypeOfFE< Mesh3 > &Elm_PKL$ORDER_3d(PKLagrangeP$ORDER);
    3) then,  you can add your newFE using:
        AddNewFE3("P$ORDER3d", &Elm_PKL$ORDER_3d); in the init fct
    4) finally recompile the file using ff-c++
  */

  /*
  static TypeOfFE_PKLagrange_3d PKLagrangeP3(3);
  GTypeOfFE< Mesh3 > &Elm_PKL3_3d(PKLagrangeP3);

  static TypeOfFE_PKLagrange_3d PKLagrangeP4(4);
  GTypeOfFE< Mesh3 > &Elm_PKL4_3d(PKLagrangeP4);
  */

  static TypeOfFE_PKLagrange_3d PKLagrangeP5(5);
  GTypeOfFE< Mesh3 > &Elm_PKL5_3d(PKLagrangeP5);

  static TypeOfFE_PKLagrange_3d PKLagrangeP6(6);
  GTypeOfFE< Mesh3 > &Elm_PKL6_3d(PKLagrangeP6);

  static TypeOfFE_PKLagrange_3d PKLagrangeP7(7);
  GTypeOfFE< Mesh3 > &Elm_PKL7_3d(PKLagrangeP7);

  static TypeOfFE_PKLagrange_3d PKLagrangeP8(8);
  GTypeOfFE< Mesh3 > &Elm_PKL8_3d(PKLagrangeP8);

  static TypeOfFE_PKLagrange_3d PKLagrangeP9(9);
  GTypeOfFE< Mesh3 > &Elm_PKL9_3d(PKLagrangeP9);

  static TypeOfFE_PKLagrange_3d PKLagrangeP10(10);
  GTypeOfFE< Mesh3 > &Elm_PKL10_3d(PKLagrangeP10);

  static TypeOfFE_PKLagrange_3d PKLagrangeP11(11);
  GTypeOfFE< Mesh3 > &Elm_PKL11_3d(PKLagrangeP11);

static void init() {

    AddNewFE("P5", &PKLagrange5);
    static ListOfTFE FE_PK5("P5", &PKLagrange5); // to add PKLagrange in list of Common FE
    AddNewFE("P6", &PKLagrange6);
    static ListOfTFE FE_PK6("P6", &PKLagrange6);
    AddNewFE("P7", &PKLagrange7);
    static ListOfTFE FE_PK7("P7", &PKLagrange7);
    AddNewFE("P8", &PKLagrange8);
    static ListOfTFE FE_PK8("P8", &PKLagrange8);
    AddNewFE("P9", &PKLagrange9);
    static ListOfTFE FE_PK9("P9", &PKLagrange9);
    AddNewFE("P10", &PKLagrange10);
    static ListOfTFE FE_PK10("P10", &PKLagrange10);
    AddNewFE("P11", &PKLagrange11);
    static ListOfTFE FE_PK11("P11", &PKLagrange11);

    /*
    AddNewFE3("P33d",  &Elm_PKL3_3d,  "P3");
    AddNewFE3("P43d",  &Elm_PKL4_3d,  "P4");
    */
    AddNewFE3("P53d",  &Elm_PKL5_3d,  "P5");
    AddNewFE3("P63d",  &Elm_PKL6_3d,  "P6");
    AddNewFE3("P73d",  &Elm_PKL7_3d,  "P7");
    AddNewFE3("P83d",  &Elm_PKL8_3d,  "P8");
    AddNewFE3("P93d",  &Elm_PKL9_3d,  "P9");
    AddNewFE3("P103d", &Elm_PKL10_3d, "P10");
    AddNewFE3("P113d", &Elm_PKL11_3d, "P11");
}

} // namespace Fem2D
LOADFUNC(Fem2D::init);
