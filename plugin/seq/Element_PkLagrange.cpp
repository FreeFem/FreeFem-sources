#include "polynomial_PKLagrange.hpp"
#include "utils_PKLagrange.hpp"
#include "AddNewFE.h"

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
static TypeOfFE_PkLagrange PKLagrange5(5);
static TypeOfFE_PkLagrange PKLagrange6(6);
static TypeOfFE_PkLagrange PKLagrange7(7);
static TypeOfFE_PkLagrange PKLagrange8(8);

static void init() {

    AddNewFE("PKLagrange5", &PKLagrange5);
    static ListOfTFE FE_PK5("PKLagrange5", &PKLagrange5); // to add PKLagrange in list of Common FE
    AddNewFE("PKLagrange6", &PKLagrange6);
    static ListOfTFE FE_PK6("PKLagrange6", &PKLagrange6); // to add PKLagrange in list of Common FE
    AddNewFE("PKLagrange7", &PKLagrange7);
    static ListOfTFE FE_PK7("PKLagrange7", &PKLagrange7); // to add PKLagrange in list of Common FE
    AddNewFE("PKLagrange8", &PKLagrange8);
    static ListOfTFE FE_PK8("PKLagrange8", &PKLagrange8); // to add PKLagrange in list of Common FE*/
}

} // namespace Fem2D
LOADFUNC(Fem2D::init);
