#include "ff++.hpp"
#include "AFunction_ext.hpp"
#include "compositeFESpace.hpp"
#ifdef PARALLELE
#include <mpi.h>
#endif

extern long mpisize, mpirank;
static const int TAG_XFER_HDR = 3000;
static const int TAG_XFER_BODY = 3001;
static const int TAG_XFER_DOF = 3002;

static const double COVER_EMPTY_TOL = 1e-12;
static const double COVER_ONE_TOL = 1e-9;
static double transferTolN = 0.0;
static double transferTolT = 0.0;

template<class Mesh>
double localMaxEdgeLength(const Mesh& Th) {
    double hmax = 0;
    for (int k = 0; k < Th.nt; ++k) {
        for (int i = 0; i < Mesh::Element::nv; ++i) {
            for (int j = i+1; j < Mesh::Element::nv; ++j) {
                hmax = std::max(hmax, (Th[k][i] - Th[k][j]).norme());
            }
        }
    }
    return hmax;
}

struct BBox { double pmin[3], pmax[3]; bool empty;};

template<class Mesh>
static BBox rawBoxOf(const Mesh& Th) {
    typename Mesh::Rd pmin, pmax;
    Th.BoundingBox(pmin, pmax);
    BBox b = {};
    b.empty = (Th.nt == 0);
    for (int d=0; d<3; ++d){
        b.pmin[d] = (d < Mesh::Rd::d) ? pmin[d] : 0.0;
        b.pmax[d] = (d<Mesh::Rd::d) ? pmax[d] : 0.0;
    }
    return b;
}

static void inflate(BBox& b, double h) {
    if (b.empty) return;
    for (int d=0; d<3; ++d) {
        b.pmin[d] -= h;
        b.pmax[d] += h;
    }
}

static bool bboxOverlap(const BBox& a, const BBox& b) {
    if (a.empty || b.empty) return false;
    for (int d = 0; d < 3; ++d) {
        if (a.pmax[d] < b.pmin[d] || b.pmax[d] < a.pmin[d]) return false;
    }
    return true;
}

template<class Rd>
static bool pointInBox(const BBox& b, const Rd& p) {
    if (b.empty) return false;
    for (int d = 0; d<Rd::d; ++d) {
        if (p[d] < b.pmin[d] || p[d] > b.pmax[d]) return false;
    }
    return true;
}

static bool boxContains(const BBox& outer, const BBox& inner) {
    if (outer.empty) return false;
    if (inner.empty) return true;
    
    for (int d = 0; d < 3; ++d) {
        if (inner.pmin[d] < outer.pmin[d] || inner.pmax[d] > outer.pmax[d]) return false;
    }
    return true;
}

struct CoverStats {
    long nRows = 0;
    long nEmpty = 0;
    long nPartial = 0;
    double vmin = 1.0;
};

static CoverStats coverageOf(const KN<double>& cover) {
    CoverStats s; s.nRows = cover.n;
    if (cover.n == 0) return s;
    s.vmin = cover[0];
    for (int i = 0; i < cover.n; ++i) {
        const double c = cover[i];
        s.vmin = std::min(s.vmin, c);
        if (std::abs(c) <= COVER_EMPTY_TOL) s.nEmpty++;
        else if (std::abs(c-1.0) > COVER_ONE_TOL) s.nPartial++;
    }
    return s;
}

static bool coversAll(const CoverStats& s) { return s.nEmpty == 0 && s.nPartial == 0; }

static KN<double> rowSums(const MatriceMorse<double>* M, int nrow) {
    KN<double> r(nrow, 0.0);
    for (size_t k = 0; k < M->nnz; ++k) r[M->i[k]] += M->aij[k];
    return r;
}

// Localisation in simplex of dimension dHat
template<int dHat> struct LocateByDim;

template<> struct LocateByDim<3> {
    template<class E>
    static bool refCoords(const E& K, const R3& x, R3& xhat, double epsT, double, double& dist2) {
        const R3 &A = K[0], &B = K[1], &C = K[2], &D = K[3];
        const double detK = 6.0*K.mesure();
        double l[4];
        l[1] = det(A, x, C, D) / detK;
        l[2] = det(A, B, x, D) / detK;
        l[3] = det(A, B, C, x) / detK;
        l[0] = 1.0 - l[1] - l[2] - l[3];
        dist2 = 0.0;
        for (int i = 0; i < 4; ++i) if (l[i] < -epsT) return false;
        xhat = R3(l[1], l[2], l[3]);
        return true;
    }
};

template<> struct LocateByDim<2> {
    template<class E>
    static bool refCoords(const E& K, const R3& x, R2& xhat, double epsT, double epsN, double& dist2) {
        const R3 &A = K[0], &B = K[1], &C = K[2];
        const R3 AB(A, B), AC(A, C), AP(A, x);
        const R3 N = AB^AC;
        const double N2 = (N, N);
        if (!(N2>0)) return false; // triangle dégénéré
        const double pn = (AP, N);
        if (pn*pn > epsN*epsN*N2) return false;
        dist2 = pn*pn/N2;
        double l[3];
        l[1] = det(AP, AC, N) / N2;
        l[2] = det(AB, AP, N) / N2;
        l[0] = 1.0 - l[1] - l[2];
        for (int i = 0; i < 3; ++i) if (l[i] < -epsT) return false;
        xhat = R2(l[1], l[2]);
        return true;
    }
};

template<> struct LocateByDim<1> {          // EdgeL dans R^3
    template<class E>
    static bool refCoords(const E& K, const R3& x, R1& xhat,
                          double epsT, double epsN, double& dist2) {
        const R3 &A = K[0], &B = K[1];
        const R3 AB(A, B), AP(A, x);
        const double ab2 = (AB, AB);
        if (!(ab2 > 0)) return false; //segment dégénéré
        const double l1  = (AP, AB) / ab2;
        if (l1 < -epsT || l1 > 1.0 + epsT) return false;
        const R3 Pj = A + l1*AB;
        dist2 = R3(Pj, x).norme2();                   // distance a la DROITE
        if (dist2 > epsN*epsN) return false;
        xhat = R1(l1);
        return true;
    }
};

template<class Mesh>
struct ElementLocator {
    typedef typename Mesh::Element Element;
    typedef typename Mesh::Rd Rd;
    typedef typename Element::RdHat RdHat;
    static bool refCoords(const Element& K, const Rd& x, RdHat& xhat, double epsT, double epsN, double& dist2)
    { return LocateByDim<RdHat::d>::refCoords(K, x, xhat, epsT, epsN, dist2); }
};

template<class Mesh>
struct FragLocator {
    double org[3], invh[3];
    int    nc[3];
    double epsT = 0.0;
    double epsN = 0.0;
    double hCell = 0.0;
    std::vector<int> cellStart;         // nCells+1, CSR
    std::vector<int> fragOf, elemOf;
    const std::vector<Mesh*>* frags = nullptr;

    void cellRange(const typename Mesh::Element& K, int* a, int* b) const {
        double lo[3], hi[3];
        for (int d = 0; d < 3; ++d) { lo[d] = 1e300; hi[d] = -1e300; }
        for (int i = 0; i < Mesh::Element::nv; ++i)
            for (int d = 0; d < 3; ++d) {
                const double c = K[i][d];
                lo[d] = std::min(lo[d], c); hi[d] = std::max(hi[d], c);
            }
        for (int d = 0; d < 3; ++d) { lo[d] -= epsN; hi[d] += epsN; }
        for (int d = 0; d < 3; ++d) {
            a[d] = (int)((lo[d] - org[d]) * invh[d]);
            b[d] = (int)((hi[d] - org[d]) * invh[d]);
            if (a[d] < 0) a[d] = 0;
            if (b[d] > nc[d]-1) b[d] = nc[d]-1;
            if (b[d] < a[d]) b[d] = a[d];
        }
    }


    void build(const std::vector<Mesh*>& F, double epsNRequested = 0.0, double epsTRequested = 0.0, double epsNFracOfH = 0.0);
    bool locate(const typename Mesh::Rd& x, int& j, int& k,
                typename Mesh::Element::RdHat& xhat, double* d2out = nullptr) const;
};

template<class Mesh>
void FragLocator<Mesh>::build(const std::vector<Mesh*>& F, double epsNRequested, double epsTRequested, double epsNFracOfH) {
    frags = &F;
    cellStart.clear(); fragOf.clear(); elemOf.clear();

    // 1. bbox de l'union + comptage
    double lo[3] = { 1e300, 1e300, 1e300 }, hi[3] = { -1e300, -1e300, -1e300 };
    long Ntot = 0;
    for (size_t j = 0; j < F.size(); ++j) {
        if (!F[j]) continue;
        Ntot += F[j]->nt;
        for (int iv = 0; iv < F[j]->nv; ++iv)
            for (int d = 0; d < 3; ++d) {
                const double c = (*F[j])(iv)[d];
                lo[d] = std::min(lo[d], c); hi[d] = std::max(hi[d], c);
            }
    }
    if (Ntot == 0) { nc[0]=nc[1]=nc[2]=0; cellStart.assign(1,0); return; }

    double ext[3]; int nActive = 0; double prodExt = 1.0, L2 = 0.0;
    for (int d = 0; d < 3; ++d) {
        ext[d] = std::max(0.0, hi[d] - lo[d]);
        L2 += ext[d]*ext[d];
        if (ext[d] > 0) { ++nActive; prodExt *= ext[d]; }
    }
    epsN = std::max(epsNRequested, 1e-12*std::sqrt(L2));
    epsT = std::max(epsTRequested, 1e-10);

    double mesTot = 0.0;
    for (size_t j = 0; j < F.size(); ++j) {
        if (!F[j]) continue;
        for (int k = 0; k < F[j]->nt; ++k) mesTot += std::abs((*F[j])[k].mesure());
    }

    const int dH = Mesh::Element::RdHat::d;
    double h = (mesTot > 0) ? std::pow(mesTot/(double)Ntot, 1.0/dH) : 0.0;

    const double nCellsMax = 4.0*(double)Ntot;
    if (nActive > 0) {
        const double hMin = std::pow(prodExt/nCellsMax, 1.0/nActive);
        if (!(h > hMin)) h = hMin;
    }
    if (!(h>0)) h = 1.0;
    hCell = h;
    if (epsNFracOfH > 0) epsN = std::max(epsN, epsNFracOfH*h);

    for (int d = 0; d < 3; ++d) {
        org[d] = lo[d];
        if (ext[d] > 0) {
            const double t = ext[d]/h + 1.0;
            nc[d]   = (t > 1e7) ? 10000000 : std::max(1, (int)t);
            invh[d] = nc[d]/ext[d];
        } else { nc[d] = 1; invh[d] = 0.0; }      // cf. 2.5
    }
    const long nCells = (long)nc[0]*nc[1]*nc[2];

    // 3. deux passes CSR ; l'ordre (j croissant, k croissant) fait le determinisme
    cellStart.assign(nCells + 1, 0);
    for (int pass = 0; pass < 2; ++pass) {
        for (size_t j = 0; j < F.size(); ++j) {
            if (!F[j]) continue;
            for (int k = 0; k < F[j]->nt; ++k) {
                int a[3], b[3];
                cellRange((*F[j])[k], a, b);          // bbox de l'element -> plage de cellules
                for (int z = a[2]; z <= b[2]; ++z)
                for (int y = a[1]; y <= b[1]; ++y)
                for (int x = a[0]; x <= b[0]; ++x) {
                    const long c = (long)x + nc[0]*((long)y + nc[1]*(long)z);
                    if (pass == 0) cellStart[c+1]++;
                    else { const int q = cellStart[c]++; fragOf[q] = (int)j; elemOf[q] = k; }
                }
            }
        }
        if (pass == 0) {
            for (long c = 0; c < nCells; ++c) cellStart[c+1] += cellStart[c];
            fragOf.resize(cellStart[nCells]); elemOf.resize(cellStart[nCells]);
        }
    }
    // la passe 2 a decale cellStart : le remettre en place
    for (long c = nCells; c > 0; --c) cellStart[c] = cellStart[c-1];
    cellStart[0] = 0;

    if (verbosity > 5) {
        const long nCellsTot = (long)nc[0]*nc[1]*nc[2];
        long nOcc = 0, mx = 0;
        for (long c = 0; c < nCellsTot; ++c) {
            const long e = cellStart[c+1] - cellStart[c];
            if (e > 0) { ++nOcc; mx = std::max(mx, e); }
        }
        cout << " -- FragLocator: dHat " << dH << " nt " << Ntot << " h " << hCell
            << " nc " << nc[0] << "x" << nc[1] << "x" << nc[2]
            << " cellules " << nCellsTot << " occupees " << nOcc
            << " balayage moyen " << (double)fragOf.size()/std::max(1L,nOcc)
            << " max " << mx << " epsN " << epsN << endl;
    }


}

template<class Mesh>
bool FragLocator<Mesh>::locate(const typename Mesh::Rd& x, int& j, int& k,
                               typename Mesh::Element::RdHat& xhat, double* d2out) const
{
    typedef typename Mesh::Element::RdHat RdHat;
    if (cellStart.size() <= 1) return false;
    int c[3];
    for (int d = 0; d < 3; ++d) {
        c[d] = (invh[d] > 0) ? (int)((x[d] - org[d]) * invh[d]) : 0;
        if (c[d] < 0) c[d] = 0;
        else if (c[d] >= nc[d]) c[d] = nc[d] - 1;
    }
    const long cc = (long)c[0] + nc[0]*((long)c[1] + nc[1]*(long)c[2]);

    int    jb = -1, kb = -1;
    double d2b = 1e300;                 // sentinelle : un NaN ne peut jamais gagner
    RdHat  xb;
    for (int q = cellStart[cc]; q < cellStart[cc+1]; ++q) {
        const int jj = fragOf[q], kk = elemOf[q];
        RdHat  xh;
        double d2;
        if (!ElementLocator<Mesh>::refCoords((*(*frags)[jj])[kk], x, xh, epsT, epsN, d2))
            continue;
        if (d2 < d2b) { jb = jj; kb = kk; d2b = d2; xb = xh; }
        if (d2 <= 0.0) break;           // optimal : aucun candidat ne peut faire mieux
    }
    if (jb < 0) return false;
    j = jb; k = kb; xhat = xb;
    if (d2out) *d2out = d2b;
    return true;
}


template<class Mesh>
static bool allDstVerticesInside(const GFESpace<Mesh>& dstVh, const FragLocator<Mesh>& L) {
    typename Mesh::Element::RdHat xh;
    int j = -1, k = -1;
    const Mesh& Th = dstVh.Th;
    for (int iv = 0; iv < Th.nv; ++iv){
        if (!L.locate(Th(iv), j, k , xh)) return false;
    }
    return true;
}

template<class Mesh1, class Mesh2>
void computeOverlapRankPairs(pcommworld comm, const DistributedMesh<Mesh1>& Dsrc, const DistributedMesh<Mesh2>& Ddst, KN<int>& sendToRanks, KN<int>& recvFromRanks, std::vector<BBox>& allSrc, std::vector<BBox>& allDst) {
    ffassert(Dsrc.comm == Ddst.comm);
    const int nproc = mpisize > 0 ? (int)mpisize : 1;
    std::vector<BBox> both(2*nproc);
    const double h = std::max(localMaxEdgeLength(*Dsrc.LocalMesh), localMaxEdgeLength(*Ddst.LocalMesh));
    BBox mySrc = rawBoxOf(*Dsrc.LocalMesh); inflate(mySrc, h);
    BBox myDst = rawBoxOf(*Ddst.LocalMesh); inflate(myDst, h);
    BBox packed[2] = {mySrc, myDst};
    #ifdef PARALLELE
        MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
        MPI_Allgather(packed, 2*sizeof(BBox), MPI_BYTE, both.data(), 2*sizeof(BBox), MPI_BYTE, cw);
    #else
        both[0] = packed[0]; both[1] = packed[1];
    #endif

    allSrc.assign(nproc, BBox()); allDst.assign(nproc, BBox());
    for (int r = 0; r < nproc; ++r) {
        allSrc[r] = both[2*r]; allDst[r] = both[2*r+1];
    }

    std::vector<int> snd, rcv;
    for (int r = 0; r < nproc; ++r) {
        if (bboxOverlap(mySrc, allDst[r])) snd.push_back(r);
        if (bboxOverlap(allSrc[r], myDst)) rcv.push_back(r);
    }
    sendToRanks = KN<int>(snd.size());
    for (int i=0; i<(int)snd.size(); ++i) sendToRanks[i] = snd[i];
    recvFromRanks = KN<int>(rcv.size());
    for (int i=0; i<(int)rcv.size(); ++i) recvFromRanks[i] = rcv[i];
}

template<class Mesh>
static Mesh* buildFragmentFor(const Mesh& Loc, const BBox& target, KN<int>& n2o, const KN<int>* usable = nullptr) {
  const int nv = Mesh::Element::nv;
  KN<int> mask(Loc.nt, 0);
  int nkept = 0;

  for (int k = 0; k < Loc.nt; ++k) {
    typename Mesh::Rd g;
    for (int i = 0; i < nv; ++i) g += Loc[k][i];
    g /= nv;
    if (pointInBox(target, g) && (!usable || (*usable)[k] >= 0)) { mask[k] = 1; ++nkept; }
  }
  if (nkept == 0) { n2o = KN<int>(0); return nullptr; }

  Mesh* frag = truncmesh(Loc, 1, (int*)mask, false, 9999, 1e-7, 1L, true, false);
  ffassert(frag && frag->nt > 0);

  //frag->BuildGTree();
  n2o = n2oFromSplit(Loc, *frag, mask);
  return frag;
}

template<class Mesh, class R>
static KN<R> fragmentDofVector(const GFESpace<Mesh>& srcVh, const Mesh& frag, const KN<int>& n2o, FragInj* injOut = nullptr)
{
    ffassert(srcVh.TFE.N() == 1);            // espaces composites hors perimetre
    GFESpace<Mesh> fragVh(frag, *srcVh.TFE[0]);   
    KN<long> inj = restrictDOFPartial(fragVh, srcVh, n2o);  
    const int nd = fragVh.NbOfDF;
    std::vector<int> vpos;
    std::vector<long> vsrc;
    if (injOut) { vpos.reserve(nd); vsrc.reserve(nd); }

    KN<R> sub(nd);
    for (int d = 0; d < nd; ++d) { 
        if (inj[d] < 0) sub[d] = R();
        else { 
            sub[d] = R(1.0); // masque 0/1 
            if (injOut) { vpos.push_back(d); vsrc.push_back(inj[d]); }
        }
    }
    
    if (injOut) {
        injOut->nd = nd;
        const long np = (long)vpos.size();
        injOut->pos.resize(np);
        injOut->src.resize(np);
        for (long k = 0; k < np; ++k) {
            injOut->pos[k] = vpos[k];
            injOut->src[k] = vsrc[k];
        }
    }    
    return sub;
}

#ifdef PARALLELE
template<class Mesh, class R, class OnFragment>
static void exchangeFragments(pcommworld comm,
                              const KN<int>& sendToRanks, const KN<int>& recvFromRanks,
                              const std::vector<Mesh*>& fragOut,      
                              const std::vector<KN<R> >& subOut,   
                              OnFragment&& onArrive, std::vector<Mesh*>* keepFrags = nullptr, std::vector<KN<R>>* keepDof = nullptr) // onArrive(int j, Mesh&, const KN<R>&)
{
    MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
    const int nS = sendToRanks.n, nR = recvFromRanks.n;
    ffassert(!keepFrags || keepDof);

    if (keepFrags) {
        keepFrags->assign(nR, nullptr);
        keepDof->assign(nR, KN<R>());
    }

    std::vector<long long> hdrIn(2*nR, 0), hdrOut(2*nS, 0);
    std::vector<Serialize*> ser(nS, nullptr);

    std::vector<MPI_Request> rqH(nR + nS);
    for (int j = 0; j < nR; ++j)
        MPI_Irecv(&hdrIn[2*j], 2, MPI_LONG_LONG, recvFromRanks[j], TAG_XFER_HDR, cw, &rqH[j]);

    for (int i = 0; i < nS; ++i) {
        if (fragOut[i]) {
            ser[i] = new Serialize(fragOut[i]->serialize());   // copie refcountee, pas cher
            hdrOut[2*i]   = (long long)ser[i]->size();
            hdrOut[2*i+1] = (long long)subOut[i].n;
            ffassert(hdrOut[2*i] < (long long)INT_MAX);        // un seul Isend, count en int
            ffassert(hdrOut[2*i+1] * (long long)sizeof(R) < (long long)INT_MAX);
        }
        MPI_Isend(&hdrOut[2*i], 2, MPI_LONG_LONG, sendToRanks[i], TAG_XFER_HDR, cw, &rqH[nR+i]);
    }
    MPI_Waitall(nR + nS, rqH.data(), MPI_STATUSES_IGNORE);

    std::vector<Serialize*> bufIn(nR, nullptr);
    std::vector<KN<R> >     dofIn(nR);          // construits par defaut : unset,
                                                // donc resize() alloue bien
    std::vector<int>        remaining(nR, 0);

    std::vector<MPI_Request> rqR;   rqR.reserve(2*nR);
    std::vector<int>         owner; owner.reserve(2*nR);
    std::vector<MPI_Request> rqS;   rqS.reserve(2*nS);

    // Tous les Irecv AVANT les Isend 
    for (int j = 0; j < nR; ++j) {
        if (hdrIn[2*j] == 0) continue;
        bufIn[j] = new Serialize((size_t)hdrIn[2*j], Fem2D::GenericMesh_magicmesh);
        dofIn[j].resize((long)hdrIn[2*j+1]);
        remaining[j] = 2;
        rqR.resize(rqR.size() + 2);
        owner.push_back(j);
        owner.push_back(j);
        MPI_Irecv((char*)*bufIn[j], (int)hdrIn[2*j], MPI_BYTE,
                  recvFromRanks[j], TAG_XFER_BODY, cw, &rqR[rqR.size()-2]);
        MPI_Irecv((R*)dofIn[j], (int)(hdrIn[2*j+1]*sizeof(R)), MPI_BYTE,
                  recvFromRanks[j], TAG_XFER_DOF, cw, &rqR[rqR.size()-1]);
    }
    for (int i = 0; i < nS; ++i) {
        if (!ser[i]) continue;
        rqS.resize(rqS.size() + 2);
        MPI_Isend((char*)*ser[i], (int)hdrOut[2*i], MPI_BYTE,
                  sendToRanks[i], TAG_XFER_BODY, cw, &rqS[rqS.size()-2]);
        MPI_Isend((R*)subOut[i], (int)(subOut[i].n*sizeof(R)), MPI_BYTE,
                  sendToRanks[i], TAG_XFER_DOF, cw, &rqS[rqS.size()-1]);
    }

    auto cleanup = [&]() {
        if (!rqR.empty()) MPI_Waitall((int)rqR.size(), rqR.data(), MPI_STATUSES_IGNORE);
        if(!rqS.empty()) MPI_Waitall((int)rqS.size(), rqS.data(), MPI_STATUSES_IGNORE);
        for (int i = 0; i < nS; ++i) { delete ser[i]; ser[i] = nullptr; }
        for (int j = 0; j < nR; ++j) { delete bufIn[j]; bufIn[j] = nullptr; }
    };

    // --- boucle d'arrivee : on interpole pendant que le reste arrive --------
    try {
        size_t done = 0;
        while (done < rqR.size()) {
            int t = MPI_UNDEFINED;
            MPI_Waitany((int)rqR.size(), rqR.data(), &t, MPI_STATUS_IGNORE);
            if (t == MPI_UNDEFINED) break;          // plus aucune requete active
            ++done;
            const int j = owner[t];
            if (--remaining[j] > 0) continue;       // l'autre moitie n'est pas arrivee

            Mesh* frag = new Mesh(*bufIn[j]);
            if (!keepFrags) frag->BuildGTree();
            delete bufIn[j];
            bufIn[j] = nullptr;

            if (!keepFrags) onArrive(j, *frag, dofIn[j]);

            if (keepFrags) {
                (*keepFrags)[j] = frag;
                (*keepDof)[j].resize(dofIn[j].n);
                (*keepDof)[j] = dofIn[j];
            }
            else { frag->destroy(); dofIn[j].resize(0); }
        }
    }
    catch(...) {
        cleanup();
        throw;
    }
    cleanup();
}
#else
template<class Mesh, class R, class OnFragment>
static void exchangeFragments(pcommworld,
                              const KN<int>& sendToRanks, const KN<int>& recvFromRanks,
                              const std::vector<Mesh*>& fragOut,
                              const std::vector<KN<R> >& subOut,
                              OnFragment&& onArrive, std::vector<Mesh*>* keepFrags = nullptr, std::vector<KN<R>>* keepDof = nullptr)
{
    ffassert(!keepFrags || keepDof);
    if (keepFrags) {
        keepFrags->assign(recvFromRanks.n, nullptr);
        keepDof->assign(recvFromRanks.n, KN<R>());
    }
    if (recvFromRanks.n && sendToRanks.n && fragOut[0]) {
        Serialize s = fragOut[0]->serialize();
        Mesh* frag = new Mesh(s);
        if (!keepFrags) frag->BuildGTree();
        if (!keepFrags) onArrive(0, *frag, subOut[0]);
        ffassert(!keepFrags || keepDof);
        if (keepFrags) {
            (*keepFrags)[0] = frag;
            (*keepDof)[0] = subOut[0];
        }
        else {
            frag->destroy();
        }
    }
}
#endif

#ifdef PARALLELE
template<class R>
static void exchangeDofOnly(pcommworld comm,
                            const KN<int>& sendToRanks, const KN<int>& recvFromRanks,
                            const std::vector<KN<R> >& subOut,
                            const std::vector<int>& ndSend,
                            const std::vector<int>& ndFrag,
                            std::vector<KN<R> >& dofIn)
{
    MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
    const int nS = sendToRanks.n, nR = recvFromRanks.n;
    ffassert((int)ndSend.size() == nS && (int)ndFrag.size() == nR);
    dofIn.assign(nR, KN<R>());

    std::vector<MPI_Request> rq;
    rq.reserve(nR + nS);                       // AUCUNE reallocation ensuite

    for (int j = 0; j < nR; ++j) {             // tous les Irecv AVANT les Isend
        if (ndFrag[j] == 0) continue;
        dofIn[j].resize(ndFrag[j]);
        rq.resize(rq.size() + 1);
        MPI_Irecv((R*)dofIn[j], (int)(ndFrag[j]*sizeof(R)), MPI_BYTE,
                  recvFromRanks[j], TAG_XFER_DOF, cw, &rq[rq.size()-1]);
    }
    for (int i = 0; i < nS; ++i) {
        if (ndSend[i] == 0) continue;
        ffassert(subOut[i].n >= ndSend[i]);
        rq.resize(rq.size() + 1);
        MPI_Isend((R*)subOut[i], (int)(ndSend[i]*(long)sizeof(R)), MPI_BYTE,
                  sendToRanks[i], TAG_XFER_DOF, cw, &rq[rq.size()-1]);
    }
    if (!rq.empty()) MPI_Waitall((int)rq.size(), rq.data(), MPI_STATUSES_IGNORE);
}
#else
template<class R>
static void exchangeDofOnly(pcommworld, const KN<int>&, const KN<int>& recvFromRanks,
                            const std::vector<KN<R> >& subOut,
                            const std::vector<int>& ndSend,
                            const std::vector<int>& ndFrag,
                            std::vector<KN<R> >& dofIn)
{
    dofIn.assign(recvFromRanks.n, KN<R>());
    if (recvFromRanks.n && !subOut.empty() && ndFrag[0] > 0 && ndSend[0] > 0) {
        dofIn[0].resize(ndFrag[0]);
        dofIn[0] = subOut[0];
    }
}
#endif


template<class Mesh1, class Mesh2, class OnFragment>
static void collectFragments(const DistributedMesh<Mesh1>& Dsrc,
                             const GFESpace<Mesh1>& srcVh,
                             const DistributedMesh<Mesh2>& Ddst,
                             KN<int>& sendToRanks, KN<int>& recvFromRanks,
                             OnFragment&& onArrive, KN<long>* sentCounts = nullptr, KN<long>* recvCounts = nullptr, std::vector<FragInj>* injOut = nullptr, std::vector<Mesh1*>* keepFrags = nullptr, std::vector<KN<double>>* keepDof = nullptr)
{
    std::vector<BBox> allSrc, allDst;
    computeOverlapRankPairs(Dsrc.comm, Dsrc, Ddst, sendToRanks, recvFromRanks, allSrc, allDst);
    if (injOut) injOut->resize(sendToRanks.n);
    if (recvCounts) { recvCounts->resize(recvFromRanks.n); *recvCounts = 0L; }
    const Mesh1* srcGeom = Dsrc.CoverMesh ? Dsrc.CoverMesh : Dsrc.LocalMesh;
    KN<int> cover2local(srcGeom->nt, -1);
    if (Dsrc.CoverMesh){
        ffassert(Dsrc.localToCoverElement.n <= srcGeom->nt);
        for (int k = 0; k < Dsrc.localToCoverElement.n; ++k){
            ffassert(Dsrc.localToCoverElement[k] >= 0 && Dsrc.localToCoverElement[k] < srcGeom->nt);
            cover2local[Dsrc.localToCoverElement[k]] = k;
        }
    }
    else{
        for (int k = 0; k < srcGeom->nt; ++k) cover2local[k] = k;
    }

    KN<int> usable(cover2local);                 // copie : cover2local reste requis plus bas
    const bool haveOwn = (Dsrc.CoverMesh && Dsrc.coverPartition.n == srcGeom->nt);
    if (haveOwn) {
        const int me = (int)mpirank;
        for (int k = 0; k < srcGeom->nt; ++k)
            if (Dsrc.coverPartition[k] != me) usable[k] = -1;
    }

    if (sentCounts) sentCounts->resize(sendToRanks.n);

    std::vector<Mesh1*> fragOut(sendToRanks.n, nullptr);
    std::vector<KN<double> > subOut(sendToRanks.n);
    for (int i = 0; i < sendToRanks.n; ++i) {
        KN<int> n2oCover;
        fragOut[i] = buildFragmentFor(*srcGeom, allDst[sendToRanks[i]], n2oCover, &usable);
        if (sentCounts) (*sentCounts)[i] = fragOut[i] ? fragOut[i]->nt : 0;
        if (fragOut[i]) {
            KN<int> n2oLocal(n2oCover.n);
            for (int kc = 0; kc < n2oCover.n; ++kc){
                n2oLocal[kc] = cover2local[n2oCover[kc]];
            }
            subOut[i] = fragmentDofVector<Mesh1, double>(srcVh, *fragOut[i], n2oLocal, injOut ? &(*injOut)[i] : nullptr);
        }
    }

    exchangeFragments(Dsrc.comm, sendToRanks, recvFromRanks, fragOut, subOut, onArrive, keepFrags, keepDof);

    for (int i = 0; i < sendToRanks.n; ++i)
        if (fragOut[i]) fragOut[i]->destroy();
}

static KN<int> dataInterpolate(int N) {
    KN<int> data(4 + N);
    data[0] = 0;
    data[1] = op_id;
    data[2] = 1;
    data[3] = 0;
    for (int c = 0; c < N; ++c) data[4 + c] = c; 
    return data;
}

struct SearchMethodGuard {
    long saved;
    SearchMethodGuard() : saved(searchMethod) { searchMethod = 0; }
    ~SearchMethodGuard() { searchMethod = saved; }
};

template<class R>
static void applyPlainOp(const FragOp& F, const KN<R>& u, KN<R>& dstU){
    for (long k = 0; k < F.ii.n; ++k)
        dstU[F.ii[k]] += F.aij[k]*u[F.jj[k]];
}

static void reportCoverage(const CoverStats& s, pcommworld comm, const char* where) {
    long glo[2] = { s.nEmpty, s.nPartial };
    double gmin = s.vmin;
    int rank = 0;
    #ifdef PARALLELE
    {
        long loc[2] = { s.nEmpty, s.nPartial };
        double vmin = s.vmin;
        MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
        MPI_Allreduce(loc,  glo,  2, MPI_LONG,   MPI_SUM, cw);
        MPI_Allreduce(&vmin,&gmin,1, MPI_DOUBLE, MPI_MIN, cw);
        MPI_Comm_rank(cw, &rank);
    }
    #endif
    if (rank != 0) return;
    static int warned = 0;
    if (glo[0] > 0 && warned < 3) {
        cerr << "Warning: " << where << " : " << glo[0]
             << " destination dof(s) are not covered by any source fragment"
                " and are set to zero (min cover = " << gmin << ")." << endl;
        if (++warned == 3) cerr << "Warning: further coverage warnings suppressed." << endl;
    }
    if (glo[1] > 0 && verbosity > 0)
        cout << " -- transferPlan: " << glo[1]
             << " destination dof(s) only partially covered (min cover = " << gmin << ")" << endl;
}

template<class Mesh>
static void buildFragmentOperators(const GFESpace<Mesh>& dstVh, const GFESpace<Mesh>& srcVh,
                                   const std::vector<Mesh*>& frags,
                                   const std::vector<KN<double> >& dofIn,
                                   const int* data,
                                   std::vector<FragOp*>& Mout, std::vector<int>& ndFrag,
                                   KN<double>& cover, double epsNRequested = 0.0, double epsTRequested = 0.0, const FragLocator<Mesh>* preBuilt = nullptr)
{
    typedef typename Mesh::Element          Element;
    typedef typename Element::RdHat         RdHat;
    typedef typename GFESpace<Mesh>::FElement FElement;

    const int nR = (int)frags.size();
    Mout.assign(nR, nullptr);
    ndFrag.assign(nR, 0);
    cover.resize(dstVh.NbOfDF); cover = 0.0;

    ffassert(srcVh.TFE.N() == 1 && dstVh.TFE.N() == 1);
    const GTypeOfFE<Mesh>* tfeSrc = srcVh.TFE[0];

    FragLocator<Mesh> Lown;
    if (!preBuilt) Lown.build(frags, epsNRequested, epsTRequested);
    const FragLocator<Mesh>& L = preBuilt ? *preBuilt : Lown;

    std::vector<GFESpace<Mesh>*> fragVh(nR, nullptr);
    std::vector<char> owned(nR, 0);
    for (int j = 0; j < nR; ++j)
        if (frags[j]) {
            if (frags[j] == &srcVh.Th) {
                fragVh[j] = const_cast<GFESpace<Mesh>*>(&srcVh);
                owned[j]  = 0;
            } else {
                fragVh[j] = new GFESpace<Mesh>(*frags[j], *tfeSrc);
                owned[j]  = 1;
            }
            ndFrag[j] = fragVh[j]->NbOfDF;
        }

    std::vector<std::vector<int> >    cooI(nR), cooJ(nR);
    std::vector<std::vector<double> > cooV(nR);

    // parametres d'interpolation : cf. lgmat.cpp:891-896
    int op = data[1];
    const int* iU2V = data + 4;
    op = (op == 3) ? op_dz : op;
    const What_d whatd = (What_d)(1 << op);
    const double eps = 1.0e-10;

    const int nbdfVK = tfeSrc->NbDoF;
    const int NVh    = tfeSrc->N;
    const int sfb1   = NVh * last_operatortype * nbdfVK;

    InterpolationMatrix<RdHat> ipmat(dstVh);
    const int nbp = ipmat.np;

    KN<double> kv(sfb1 * nbp);
    double* v = kv;
    std::vector<int>   jf(nbp, -1), kf(nbp, -1);
    std::vector<char>  found(nbp, 0);
    std::vector<RdHat> xh(nbp);

    KN<bool> fait(dstVh.NbOfDF); fait = false;

    for (int it = 0; it < dstVh.Th.nt; ++it) {
        FElement KU = dstVh[it];

        // si tous les DDL de l'element sont deja traites, aucune localisation
        bool todo = false;
        for (int df = 0; df < KU.NbDoF(); ++df) if (!fait[KU(df)]) { todo = true; break; }
        if (!todo) continue;

        ipmat.set(KU);
        const Element& TU = dstVh.Th[it];

        for (int p = 0; p < nbp; ++p) {
            int j = -1, k = -1;
            found[p] = L.locate(TU(ipmat.P[p]), j, k, xh[p]) ? 1 : 0;
            jf[p] = j; kf[p] = k;
            if (found[p]) {
                KNMK_<double> fb(v + p*sfb1, nbdfVK, NVh, last_operatortype);
                tfeSrc->FB(whatd, *frags[j], (*frags[j])[k], xh[p], fb);
            }
        }

        for (int i = 0; i < ipmat.ncoef; ++i) {
            const int dfu = KU(ipmat.dofe[i]);
            if (fait[dfu]) continue;
            const int p = ipmat.p[i];
            if (!found[p]) continue;                 // remplace le test intV[p] (inside)
            const int jU = ipmat.comp[i];
            const int jV = iU2V ? iU2V[jU] : jU;
            if (jV < 0 || jV >= NVh) continue;
            const double aipj = ipmat.coef[i];
            const int j = jf[p], k = kf[p];
            FElement KV = (*fragVh[j])[k];
            KNMK_<double> fb(v + p*sfb1, nbdfVK, NVh, last_operatortype);
            KN_<double> fbj(fb('.', jV, op));
            for (int idfv = 0; idfv < nbdfVK; ++idfv)
                if (std::abs(fbj[idfv]) > eps) {
                    const double c = fbj[idfv] * aipj;
                    if (std::abs(c) > eps) {
                        cooI[j].push_back(dfu);
                        cooJ[j].push_back(KV(idfv));
                        cooV[j].push_back(c);
                    }
                }
        }

        // marquage INCONDITIONNEL de tous les DDL de l'element : cf.
        // lgmat.cpp:967-971. Ne marquer que les DDL resolus resoudrait
        // davantage de lignes et ferait bouger les comptes de couverture.
        for (int df = 0; df < KU.NbDoF(); ++df) fait[KU(df)] = true;
    }

    // compactage en FragOp + cover = somme_j M_j * chi_j
    for (int j = 0; j < nR; ++j) {
        if (cooI[j].empty()) continue;
        FragOp* F = new FragOp;
        F->nrow = dstVh.NbOfDF;
        F->ncol = ndFrag[j];
        const long nz = (long)cooI[j].size();
        F->ii.resize(nz); F->jj.resize(nz); F->aij.resize(nz);
        for (long q = 0; q < nz; ++q) {
            F->ii[q] = cooI[j][q]; F->jj[q] = cooJ[j][q]; F->aij[q] = cooV[j][q];
        }
        Mout[j] = F;
        // dofIn[j] = masque 0/1, un par DDl du fragment
        if (dofIn[j].n >= (long)ndFrag[j])
            for (long q = 0; q < nz; ++q)
                cover[F->ii[q]] += F->aij[q] * dofIn[j][F->jj[q]];
    }

    for (int j = 0; j < nR; ++j) if (owned[j]) delete fragVh[j];
}

template<class Mesh>
static bool detectTolerances(const GFESpace<Mesh>& dstVh, const std::vector<Mesh*>& frags,
                             pcommworld comm, double& tolN, double& tolT)
{
    double dMax = 0.0, rel = 0.0, hRef = 0.0;
    long   nOK = 0;
    {
        FragLocator<Mesh> L;
        L.build(frags, 0.0, 0.0, 0.25);      // epsN = h/4, epsT au plancher
        hRef = L.hCell;
        if (hRef > 0) {
            const Mesh& Th = dstVh.Th;
            typename Mesh::Element::RdHat xh;
            int j = -1, k = -1;
            for (int iv = 0; iv < Th.nv; ++iv) {
                double d2 = 0.0;
                if (L.locate(Th(iv), j, k, xh, &d2)) { ++nOK; dMax = std::max(dMax, std::sqrt(d2)); }
            }
            rel = dMax/hRef;
        }
    }

    double loc[3] = { dMax, rel, hRef }, glo[3] = { loc[0], loc[1], loc[2] };
    long   lc = nOK, gc = nOK;
    #ifdef PARALLELE
    if (mpisize > 1) {
        MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
        MPI_Allreduce(loc, glo, 3, MPI_DOUBLE, MPI_MAX, cw);
        MPI_Allreduce(&lc, &gc, 1, MPI_LONG, MPI_SUM, cw);
    }
    #endif
    hRef = glo[2];                            // hRef GLOBAL : la decision doit etre
    if (!(hRef > 0)) return false;    

    if (gc == 0) {                            // aucun echantillon : indecidable
        if (verbosity > 0 && mpirank == 0)
            cerr << "Warning: transferPlan: tolerances non mesurables (aucun point"
                    " destination interieur a la source) ; tolN/tolT inchangees." << endl;
        return false;
    }
    const double dg = glo[0]; rel = glo[1];
    if (dg <= 1e-10*hRef) return true;        // conforme : les planchers de build() suffisent

    tolN = std::max(tolN, std::pow(2.0, std::ceil(std::log2(4.0*dg))));
    tolT = std::max(tolT, std::pow(2.0, std::ceil(std::log2(32.0*rel*rel))));
    if (verbosity > 0 && mpirank == 0)
        cout << " -- transferPlan: ecart geometrique " << dg << " (h " << hRef
             << ", relatif " << rel << ", " << gc << " echantillons)"
             << " -> tolN " << tolN << " tolT " << tolT << endl;
    return true;
}

static FragOp* morseToFragOp(const MatriceMorse<double>& A) {
    ffassert(A.fortran == 0 && A.half == 0);
    FragOp* F = new FragOp;
    F->nrow = A.n;
    F->ncol = A.m;
    const long nz = (long)A.nnz;
    F->ii.resize(nz); F->jj.resize(nz); F->aij.resize(nz);
    for (long k = 0; k < nz; ++k) {
        F->ii[k] = A.i[k]; F->jj[k] = A.j[k]; F->aij[k] = A.aij[k];
    }
    return F;
}

template<class Mesh>
static TransferPlan<Mesh>* buildTransferPlan(const DistributedMesh<Mesh>& Dsrc, const GFESpace<Mesh>& srcVh, const DistributedMesh<Mesh>& Ddst, const GFESpace<Mesh>& dstVh)
{

    ffassert(srcVh.N == dstVh.N);
    ffassert(srcVh.TFE.N() == 1 && dstVh.TFE.N() == 1);

    TransferPlan<Mesh>* P = new TransferPlan<Mesh>();
    P->Dsrc = &Dsrc; Dsrc.increment();
    P->Ddst = &Ddst; Ddst.increment();
    P->srcVh = &srcVh; P->srcTh = &srcVh.Th;
    P->dstVh = &dstVh; P->dstTh = &dstVh.Th;
    P->nSrcDof = srcVh.NbOfDF;
    P->nDstDof = dstVh.NbOfDF;
    P->comm = Dsrc.comm;

    KN<int> data = dataInterpolate(dstVh.N);

    if (&Dsrc == &Ddst && &srcVh.Th == &dstVh.Th) {
        MatriceMorse<double>* A = buildInterpolationMatrixT(dstVh, srcVh, (int*)data);
        P->M.assign(1, morseToFragOp(*A));
        delete A;
        P->ndFrag.assign(1, srcVh.NbOfDF);
        P->path = XFER_LOCAL;
        return P;
    }

    SearchMethodGuard sg;
    double tolN = std::max(0.0, transferTolN);
    double tolT = std::max(0.0, transferTolT);
    #ifdef PARALLELE
    if (Dsrc.comm || mpisize > 1) {
        MPI_Comm cw = Dsrc.comm ? *(MPI_Comm*)Dsrc.comm : MPI_COMM_WORLD;
        double t[2] = {tolN, tolT}, g[2];
        MPI_Allreduce(&t, g, 2, MPI_DOUBLE, MPI_MAX, cw);
        tolN = g[0];
        tolT = g[1];
    }
    #endif

    const bool autoTol = (tolN <= 0.0 && tolT <= 0.0);
    if (autoTol) {
        std::vector<Mesh*> oneT(1, const_cast<Mesh*>(&srcVh.Th));
        detectTolerances(dstVh, oneT, Dsrc.comm, tolN, tolT);
    }

    const int nproc = mpisize > 0 ? (int)mpisize : 1;
    if (nproc == 1) {
        P->chi = interpolatePoU(Dsrc, srcVh);
        std::vector<Mesh*> one1(1, const_cast<Mesh*>(&srcVh.Th));
        std::vector<KN<double>> mask1(1);
        mask1[0].resize(srcVh.NbOfDF); mask1[0] = P->chi;
        buildFragmentOperators(dstVh, srcVh, one1, mask1, data, P->M, P->ndFrag, P->cover, tolN, tolT);
        if (P->M.empty())  P->M.assign(1, (FragOp*)nullptr);
        if (!P->M[0]) {
            FragOp* F0 = new FragOp;
            F0->nrow = P->nDstDof;
            F0->ncol = srcVh.NbOfDF;
            P->M[0] = F0;
        }
        reportCoverage(coverageOf(P->cover), nullptr, "interpolateD");
        P->path = XFER_SINGLE_RANK;
        return P;        
    }

    const BBox bs = rawBoxOf(srcVh.Th), bd = rawBoxOf(dstVh.Th);
    std::vector<Mesh*>       one(1, const_cast<Mesh*>(&srcVh.Th));
    std::vector<KN<double> > mask(1);
    FragLocator<Mesh>        Lv;
    bool haveLv  = false;
    int  localCand = 0;

    const bool codim = ((int)Mesh::Element::RdHat::d < (int)Mesh::Rd::d);
    if (boxContains(bs, bd)) {
        Lv.build(one, tolN, tolT); haveLv = true;
        bool cand;
        if (!codim && dstVh.TFE[0]->ndfonVertex == 0) cand = true;   // volume, aucun DDL sommet
        else cand = allDstVerticesInside(dstVh, Lv);
        localCand = cand ? 1 : 0;
    }


    int globalCand = localCand;
    #ifdef PARALLELE
    {
        MPI_Comm cw = Dsrc.comm ? *(MPI_Comm*)Dsrc.comm : MPI_COMM_WORLD;
        MPI_Allreduce(&localCand, &globalCand, 1, MPI_INT, MPI_MIN, cw);
    }
    #endif

    std::vector<FragOp*> Mtry;
    std::vector<int>     ndTry;
    KN<double>           coverTry;
    int globalOK = 0;
    if (globalCand) {
        mask[0].resize(srcVh.NbOfDF); mask[0] = 1.0;
        buildFragmentOperators(dstVh, srcVh, one, mask, data,
                               Mtry, ndTry, coverTry,
                               tolN, tolT, haveLv ? &Lv : nullptr);
        const int localOK = coversAll(coverageOf(coverTry)) ? 1 : 0;
        globalOK = localOK;
        #ifdef PARALLELE
        {
            MPI_Comm cw = Dsrc.comm ? *(MPI_Comm*)Dsrc.comm : MPI_COMM_WORLD;
            MPI_Allreduce(&localOK, &globalOK, 1, MPI_INT, MPI_MIN, cw);
        }
        #endif
    }


    if (globalOK) {
        P->M = Mtry;
        P->ndFrag = ndTry;
        P->path = XFER_LOCAL;
        P->chi.resize(0);
        P->cover.resize(0);
        return P;
    }
    for (size_t j = 0; j < Mtry.size(); ++j) delete Mtry[j];

    std::vector<Mesh*>        keepFrags;
    std::vector<KN<double> >  keepDof;

    collectFragments(Dsrc, srcVh, Ddst, P->sendToRanks, P->recvFromRanks,
                     [](int, const Mesh&, const KN<double>&) {}, 
                     nullptr, nullptr, &P->inj, &keepFrags, &keepDof);


    if (autoTol) detectTolerances(dstVh, keepFrags, Dsrc.comm, tolN, tolT);
    try {
        buildFragmentOperators(dstVh, srcVh, keepFrags, keepDof, data,
                               P->M, P->ndFrag, P->cover,
                               tolN, tolT);
    }
    catch (...) {
        for (size_t j = 0; j < keepFrags.size(); ++j)
            if (keepFrags[j]) keepFrags[j]->destroy();
        throw;
    }
    for (size_t j = 0; j < keepFrags.size(); ++j)
        if (keepFrags[j]) keepFrags[j]->destroy();

    ffassert((int)P->M.size() == P->recvFromRanks.n);
    ffassert((int)P->ndFrag.size() == P->recvFromRanks.n);

    reportCoverage(coverageOf(P->cover), P->comm, "transferPlan");
    P->path = XFER_GENERAL;
    return P;
}

template<class Mesh, class R>
static void applyTransferPlan(const TransferPlan<Mesh>& P, const KN<R>& srcU, KN<R>& dstU)
{
    ffassert(srcU.n == P.nSrcDof && dstU.n == P.nDstDof);
    dstU = R();

    if (P.path == XFER_LOCAL) {                     // ni chi ni renormalisation
        applyPlainOp(*P.M[0], srcU, dstU);
        return;
    }

    if (P.path == XFER_SINGLE_RANK) {
        KN<R> scaledU(srcU.n);
        for (int d = 0; d < srcU.n; ++d) scaledU[d] = srcU[d] * P.chi[d];
        applyPlainOp(*P.M[0], scaledU, dstU); 
    } else {
        std::vector<KN<R> > subOut(P.sendToRanks.n);
        std::vector<int> ndSend(P.sendToRanks.n, 0);
        ffassert((int)P.inj.size() == P.sendToRanks.n);
        for (int i = 0; i < P.sendToRanks.n; ++i) {
            const FragInj& J = P.inj[i];
            ndSend[i] = J.nd;
            subOut[i].resize(J.nd);  if (J.nd > 0) subOut[i] = R(); 
            for (long k = 0; k < J.pos.n; ++k)
                subOut[i][J.pos[k]] = srcU[J.src[k]];
        }
        std::vector<KN<R> > dofIn;
        exchangeDofOnly(P.comm, P.sendToRanks, P.recvFromRanks, subOut, ndSend, P.ndFrag, dofIn);
        for (int j = 0; j < P.recvFromRanks.n; ++j)
            if (P.M[j]) applyPlainOp(*P.M[j], dofIn[j], dstU);
    }

    for (int i = 0; i < dstU.n; ++i)
        if (std::abs(P.cover[i]) > COVER_EMPTY_TOL) dstU[i] /= P.cover[i];
}

template<class Mesh, class T>
static void exchangeRawOnPlan(const TransferPlan<Mesh>& P, const KN<T>& srcValues, T fillHoles, std::vector<KN<T> >& perFragment) {
    ffassert(srcValues.n == P.nSrcDof);
    ffassert(P.path == XFER_GENERAL);

    std::vector<KN<T>> subOut(P.sendToRanks.n);
    std::vector<int> ndSend(P.sendToRanks.n,0);
    ffassert((int)P.inj.size() == P.sendToRanks.n);
    for (int i = 0; i < P.sendToRanks.n; ++i) {
        const FragInj& J = P.inj[i];
        ndSend[i] = J.nd;
        subOut[i].resize(J.nd);
        if (J.nd > 0) subOut[i] = fillHoles;
        for (long k = 0; k < J.pos.n; ++k)
            subOut[i][J.pos[k]] = srcValues[J.src[k]];
    }
    exchangeDofOnly(P.comm, P.sendToRanks, P.recvFromRanks, subOut, ndSend, P.ndFrag, perFragment);
}

namespace {
    struct XTrip {int i; long g; double v; };
    inline bool xtripLess(const XTrip& a, const XTrip& b) {
        return a.i != b.i ? a.i < b.i : a.g < b.g;
    }
}

template<class Mesh>
static MatriceMorse<double>* assembleTransferMatrix(const TransferPlan<Mesh>& P, const KN<long>& srcNumbering, KN<long>& colGlobal) {
    ffassert(srcNumbering.n == P.nSrcDof);
    const bool weighted = (P.path != XFER_LOCAL);

    std::vector<XTrip> T;

    if (P.path == XFER_GENERAL) {
        std::vector<KN<long> > globFrag;
        exchangeRawOnPlan(P, srcNumbering, -1L, globFrag);

        for (int j = 0; j < P.recvFromRanks.n; ++j) {
            if (!P.M[j]) continue;
            const FragOp& F = *P.M[j];
            for (long k = 0; k < F.ii.n; ++k) {
                const int c = F.jj[k];
                const long g = globFrag[j][c];
                if (g < 0) continue;          // trou : c'est le masque, poids nul
                XTrip t; t.i = F.ii[k]; t.g = g; t.v = F.aij[k];   // poids 1
                T.push_back(t);
            }
        }
    }
    
    else {
        ffassert(P.M.size() == 1 && P.M[0]);
        const FragOp& F = *P.M[0];
        for (long k = 0; k < F.ii.n; ++k) {
            const int c = F.jj[k];
            const double w = weighted ? P.chi[c] : 1.0;
            if (w == 0.0) continue;
            XTrip t; t.i = F.ii[k]; t.g = srcNumbering[c]; t.v = F.aij[k]*w;
            T.push_back(t);
        }
    }

    std::sort(T.begin(), T.end(), xtripLess);

    size_t nw = 0;
    for (size_t r = 0; r < T.size(); ) { // avancement par s en fin de corps
        size_t s = r; double acc = 0.0;
        while (s < T.size() && T[s].i == T[r].i && T[s].g == T[r].g) { acc += T[s].v; ++s; }
        T[nw].i = T[r].i; T[nw].g = T[r].g; T[nw].v = acc; ++nw;
        r = s;
    }
    T.resize(nw);

    if (weighted) {
        size_t keep = 0;
        for (size_t k = 0; k < T.size(); ++k) {
            const double cv = P.cover[T[k].i];
            if (std::abs(cv) <= COVER_EMPTY_TOL) continue;
            T[keep] = T[k]; T[keep].v /= cv; ++keep;
        }
        T.resize(keep);
    }

    std::map<long, int> pos;
    for (int d = 0; d < P.nSrcDof; ++d) pos[srcNumbering[d]] = d;

    std::vector<long> extra;
    for (size_t k = 0; k < T.size(); ++k) {
        if (pos.find(T[k].g) == pos.end()) {
            pos[T[k].g] = P.nSrcDof + (int)extra.size();
            extra.push_back(T[k].g);
        }
    }

    const long m = (long)P.nSrcDof + (long)extra.size();
    colGlobal.resize(m);
    for (int d = 0; d < P.nSrcDof; ++d)          colGlobal[d] = srcNumbering[d];
    for (size_t e = 0; e < extra.size(); ++e)    colGlobal[(long)P.nSrcDof + (long)e] = extra[e];

    MatriceMorse<double>* A = new MatriceMorse<double>(P.nDstDof, (int)m, 0, 0);
    for (size_t k = 0; k < T.size(); ++k)
        (*A)(T[k].i, pos[T[k].g]) += T[k].v;
    return A;
}


template<class Mesh, class R>
static int interpolateDistributed(const DistributedMesh<Mesh>& Dsrc,
                                   const GFESpace<Mesh>& srcVh, const KN<R>& srcU,
                                   const DistributedMesh<Mesh>& Ddst,
                                   const GFESpace<Mesh>& dstVh, KN<R>& dstU)
{
    ffassert(srcVh.N == dstVh.N);
    ffassert(srcVh.TFE.N() == 1 && dstVh.TFE.N() == 1);
    ffassert(srcU.n == srcVh.NbOfDF);
    ffassert(dstU.n == dstVh.NbOfDF);

    TransferPlan<Mesh>* P = buildTransferPlan(Dsrc, srcVh, Ddst, dstVh);
    applyTransferPlan(*P, srcU, dstU);
    const int path = P->path;
    P->destroy();
    return path;
}


template<class Mesh>
long transferFragmentCounts(const DistributedMesh<Mesh>** const & Dsrc,
                            const DistributedMesh<Mesh>** const & Ddst,
                            KN<long>* const & nEnvoyes, KN<long>* const & nRecus)
{
    throwassert(Dsrc && *Dsrc && Ddst && *Ddst && nEnvoyes && nRecus);
    const DistributedMesh<Mesh>& A = **Dsrc;
    const DistributedMesh<Mesh>& B = **Ddst;

    GFESpace<Mesh> srcVh(*A.LocalMesh, DataFE<Mesh>::P1);

    KN<int> snd, rcv;
    collectFragments(A, srcVh, B, snd, rcv, [&](int j, const Mesh& frag, const KN<double>&) { (*nRecus)[j] = frag.nt; }, nEnvoyes, nRecus);
    return 0L;
}


template<class Mesh>
long overlapRankPairs(const DistributedMesh<Mesh>** const & Dsrc, const DistributedMesh<Mesh>** const & Ddst, KN<long>* const & sendToRanks, KN<long>* const & recvFromRanks) {
    throwassert(Dsrc && *Dsrc && Ddst && *Ddst);
    KN<int> snd, rcv;
    std::vector<BBox> a, b;
    computeOverlapRankPairs((**Dsrc).comm, **Dsrc, **Ddst, snd, rcv, a, b);
    sendToRanks->resize(snd.n);
    for (int i = 0; i < snd.n; ++i) (*sendToRanks)[i] = snd[i]; 
    recvFromRanks->resize(rcv.n);
    for (int i = 0; i < rcv.n; ++i) (*recvFromRanks)[i] = rcv[i];
    return 0L;
}

template<class Mesh>
long transferP1(const DistributedMesh<Mesh>** const & Dsrc, KN<double>* const & uSrc,
                const DistributedMesh<Mesh>** const & Ddst, KN<double>* const & uDst)
{
    throwassert(Dsrc && *Dsrc && Ddst && *Ddst && uSrc && uDst);
    GFESpace<Mesh> srcVh(*(**Dsrc).LocalMesh, DataFE<Mesh>::P1);
    GFESpace<Mesh> dstVh(*(**Ddst).LocalMesh, DataFE<Mesh>::P1);
    ffassert(uSrc->n == srcVh.NbOfDF);
    uDst->resize(dstVh.NbOfDF);            // resize, pas operator= : cf. le bug du jalon 2
    int xfer = interpolateDistributed(**Dsrc, srcVh, *uSrc, **Ddst, dstVh, *uDst);
    return (long)xfer;
}

template<class Mesh, class R>
static int interpolateFE(v_dfes<Mesh>* pSrc, const GFESpace<Mesh>& srcVh, const KN<R>& uSrc, v_dfes<Mesh>* pDst, const GFESpace<Mesh>& dstVh, KN<R>& uDst){
    ffassert(pSrc && pDst && pSrc->DTh && pDst->DTh);
    ffassert(pSrc->DTh->comm == pDst->DTh->comm);
    ffassert(uSrc.n == srcVh.NbOfDF && uDst.n == dstVh.NbOfDF);
    return interpolateDistributed(*pSrc->DTh, srcVh, uSrc, *pDst->DTh, dstVh, uDst);
}

template<class Mesh, class R>
long interpolateD(v_dfes<Mesh>** const & ppSrc, KN<R>* const & uSrc,
                  v_dfes<Mesh>** const & ppDst, KN<R>* const & uDst)
{
    throwassert(ppSrc && *ppSrc && ppDst && *ppDst && uSrc && uDst);
    v_dfes<Mesh>* pSrc = *ppSrc;
    v_dfes<Mesh>* pDst = *ppDst;

    const GFESpace<Mesh>& srcVh = **pSrc;           // operator FESpace*() : lgfem.hpp:427
    const GFESpace<Mesh>& dstVh = **pDst;

    int xfer = interpolateFE(pSrc, srcVh, *uSrc, pDst, dstVh, *uDst);
    return (long)xfer;
}

template<class Mesh>
const TransferPlan<Mesh>* makeTransferPlan(v_dfes<Mesh>** const& ppSrc,
                                           v_dfes<Mesh>** const& ppDst)
{
    throwassert(ppSrc && *ppSrc && ppDst && *ppDst);
    v_dfes<Mesh>* pSrc = *ppSrc;
    v_dfes<Mesh>* pDst = *ppDst;
    throwassert(pSrc->DTh && pDst->DTh);
    ffassert(pSrc->DTh->comm == pDst->DTh->comm);
    const GFESpace<Mesh>& srcVh = **pSrc;      // operator FESpace*() : lgfem.hpp:428
    const GFESpace<Mesh>& dstVh = **pDst;
    return buildTransferPlan(*pSrc->DTh, srcVh, *pDst->DTh, dstVh);
}

template<class Mesh, class R>
long interpolateDPlan(const TransferPlan<Mesh>** const& ppP,
                      KN<R>* const& uSrc, KN<R>* const& uDst)
{
    throwassert(ppP && *ppP && uSrc && uDst);
    const TransferPlan<Mesh>& P = **ppP;
    if (uSrc->n != P.nSrcDof || uDst->n != P.nDstDof)
        ExecError("interpolateD(plan, u[], v[]) : outdated plan — "
                  "the FE spaces have changed size since transferPlan()");
    applyTransferPlan(P, *uSrc, *uDst);
    return (long)P.path;
}

template<class Mesh>
long probeRawTransfer(const TransferPlan<Mesh>** const& ppP,
                      KN<long>* const& srcNumbering,
                      KN<long>* const& out)
{
    throwassert(ppP && *ppP && srcNumbering && out);
    const TransferPlan<Mesh>& P = **ppP;
    out->resize(5); *out = 0L; (*out)[3] = LONG_MAX; (*out)[4] = -1L;
    if (P.path != XFER_GENERAL) return (long)P.path;

    std::vector<KN<long> > globFrag;
    exchangeRawOnPlan(P, *srcNumbering, -1L, globFrag);

    for (int j = 0; j < P.recvFromRanks.n; ++j)
        for (long c = 0; c < globFrag[j].n; ++c) {
            const long g = globFrag[j][c];
            (*out)[0]++;
            if (g < 0) { (*out)[1]++; continue; }
            (*out)[3] = std::min((*out)[3], g);
            (*out)[4] = std::max((*out)[4], g);
        }
    return (long)P.path;
}

template<class Mesh>
long transferCoverage(const TransferPlan<Mesh>** const& ppP, KN<double>* const& out) {
    throwassert(ppP && *ppP && out);
    const TransferPlan<Mesh>& P = **ppP;
    out->resize(4);
    if (P.path == XFER_LOCAL) {
        (*out)[0] = 0; (*out)[1] = 0; (*out)[2] = 1.0; (*out)[3] = P.nDstDof;
        return (long)P.path;
    }
    const CoverStats s = coverageOf(P.cover);
    (*out)[0] = (double)s.nEmpty;
    (*out)[1] = (double)s.nPartial;
    (*out)[2] = s.vmin;
    (*out)[3] = (double)s.nRows;
    return (long)P.path;
}

template<class Mesh>
long interpolateMat(const TransferPlan<Mesh>** const & ppP, KN<long>* const& srcNumbering, Matrice_Creuse<double>* const& out, KN<double>* const& colNumbering) {
    throwassert(ppP && *ppP && srcNumbering && out && colNumbering);
    const TransferPlan<Mesh>& P = **ppP;
    if (srcNumbering->n != P.nSrcDof)
        ExecError("interpolateMat : source numbering not compatible with the plan");

    KN<long> colGlobal;
    MatriceMorse<double>* A = assembleTransferMatrix(P, *srcNumbering, colGlobal);
    out->A.master(A);

    colNumbering->resize(colGlobal.n);
    for (long k = 0; k <colGlobal.n; ++k) (*colNumbering)[k] = (double)colGlobal[k];
    return (long)P.path;
}

template<class Mesh>
static bool sameDistributedSpace(const v_dfes<Mesh>* a, const GFESpace<Mesh>& Va, const v_dfes<Mesh>* b, const GFESpace<Mesh>& Vb) {
    if (a==b) return true;
    return (a->DTh && a->DTh == b->DTh && a->nbcperiodic == 0 && b->nbcperiodic == 0 && Va.TFE.N() == 1 && Vb.TFE.N() == 1 && Va.TFE[0] == Vb.TFE[0]);
}

template<class Mesh, class R>
static void assignFEDistributed(FEbase<R, v_dfes<Mesh>>* dst, FEbase<R, v_dfes<Mesh>>* src){
    if (dst == src) return;

    v_dfes<Mesh>* pS = *src->pVh; v_dfes<Mesh>* pD = *dst->pVh;
    ffassert(pS && pS->DTh);
    ffassert(pD && pD->DTh);

    if (pS->N != pD->N) ExecError("u = v: incompatible number of components");

    KN<R>* xs = src->x();
    if (!xs) {
        const GFESpace<Mesh>* V = src->newVh();
        *src = xs = new KN<R>(V->NbOfDF);
        *xs = R();
    }
    const GFESpace<Mesh>& VhS = (src->Vh) ? *src->Vh : *src->newVh();
    const GFESpace<Mesh>& VhD = *dst->newVh();

    if (xs->N() != VhS.NbOfDF) ExecError("u = v: outdated source FE function");

    if (sameDistributedSpace(pS, VhS, pD, VhD)) {
        *dst = new KN<R>(*xs);
        return;
    }

    KN<R>* y = new KN<R>(VhD.NbOfDF);
    try { interpolateFE(pS, VhS, *xs, pD, VhD, *y); }
    catch (...) { delete y; throw; }
    *dst = y;
    return;
}

template<class Mesh, class R>
static std::pair<FEbase<R, v_dfes<Mesh>>*, int>
setFEDistributed(const std::pair<FEbase<R, v_dfes<Mesh>>*, int>& dst, const std::pair<FEbase<R, v_dfes<Mesh>>*, int>& src) {
    if ((*src.first->pVh)->N != 1 || (*dst.first->pVh)->N != 1) ExecError("b1 = a1 on a vector FE function: use [b1,...]=[a1,...]");
    assignFEDistributed<Mesh, R>(dst.first, src.first);
    return dst;
}

template<class Mesh, class R>
static FEbase<R, v_dfes<Mesh>>** initFEDistributed(FEbase<R, v_dfes<Mesh>>** const& p, v_dfes<Mesh>** const& a, const std::pair<FEbase<R, v_dfes<Mesh>>*, int>& src) {
    *p = new FEbase<R, v_dfes<Mesh>>(a);
    setFEDistributed<Mesh, R>(std::make_pair(*p, 0), src);
    return p;
}

template<class Mesh, class R>
class E_assignFEDistributed : public E_F0mps {
    typedef FEbase<R, v_dfes<Mesh>> FE;
    Expression dst, src;
public:
    E_assignFEDistributed(Expression d, Expression s) : dst(d), src(s) {}

    AnyType operator()(Stack s) const { 
        FE* d = *GetAny<FE**>((*dst)(s));
        FE* r = *GetAny<FE**>((*src)(s));
        assignFEDistributed<Mesh, R>(d, r); 
        return Nothing;
    }
    operator aType() const { return atype<void>(); }
};

template<class Mesh, class R>
Expression newAssignFEDistributed(Expression dst, Expression src) {
    return new E_assignFEDistributed<Mesh, R>(dst, src);
}

template Expression newAssignFEDistributed<Mesh3, double >(Expression, Expression);
template Expression newAssignFEDistributed<Mesh3, Complex>(Expression, Expression);
template Expression newAssignFEDistributed<MeshS, double >(Expression, Expression);
template Expression newAssignFEDistributed<MeshS, Complex>(Expression, Expression);
template Expression newAssignFEDistributed<MeshL, double >(Expression, Expression);
template Expression newAssignFEDistributed<MeshL, Complex>(Expression, Expression);



static void registerTransferGlobals() {
        static bool done = false;
        if (done) return;
        done = true;
        Global.New("transferTolN", CPValue<double>(transferTolN));
        Global.New("transferTolT", CPValue<double>(transferTolT));
    }

template<class Mesh>
void registerTransferInterpolateOps() {
    registerTransferGlobals();
    typedef const DistributedMesh<Mesh>** DMP;
    Global.Add("overlapRankPairs", "(",
        new OneOperator4_<long, DMP, DMP, KN<long>*, KN<long>*>(overlapRankPairs<Mesh>));
    Global.Add("transferFragmentCounts", "(",
        new OneOperator4_<long, DMP, DMP, KN<long>*, KN<long>*>(transferFragmentCounts<Mesh>));
    Global.Add("transferP1", "(",
        new OneOperator4_<long, DMP, KN<double>*, DMP, KN<double>*>(transferP1<Mesh>));
    Global.Add("interpolateD", "(",
        new OneOperator4_<long, v_dfes<Mesh>**, KN<double>*, v_dfes<Mesh>**, KN<double>*>(
            interpolateD<Mesh,double>));
    Global.Add("interpolateD", "(",
        new OneOperator4_<long, v_dfes<Mesh>**, KN<Complex>*, v_dfes<Mesh>**, KN<Complex>*>(
            interpolateD<Mesh,Complex>));

    typedef const TransferPlan<Mesh>*  TP;
    typedef const TransferPlan<Mesh>** TPP;

    TheOperators->Add("<-", new OneOperator2_<TP*, TP*, TP>(&set_copy_incr));
    TheOperators->Add("=",  new OneOperator2 <TP*, TP*, TP>(&set_eqdestroy_incr));

    Global.Add("transferPlan", "(",
        new OneOperator2_<TP, v_dfes<Mesh>**, v_dfes<Mesh>**,
                          E_F_F0F0_Add2RC<TP, v_dfes<Mesh>**, v_dfes<Mesh>**> >(
            makeTransferPlan<Mesh>));

    Global.Add("interpolateD", "(",
        new OneOperator3_<long, TPP, KN<double>*,  KN<double>* >(interpolateDPlan<Mesh,double>));
    Global.Add("interpolateD", "(",
        new OneOperator3_<long, TPP, KN<Complex>*, KN<Complex>*>(interpolateDPlan<Mesh,Complex>));

    Global.Add("probeRawTransfer", "(",
        new OneOperator3_<long, TPP, KN<long>*, KN<long>*>(probeRawTransfer<Mesh>));

    Global.Add("transferCoverage", "(",
        new OneOperator2_<long, TPP, KN<double>*>(transferCoverage<Mesh>));


    Global.Add("interpolateMat", "(",
        new OneOperator4_<long, TPP, KN<long>*, Matrice_Creuse<double>*, KN<double>*>(
            interpolateMat<Mesh>));

    typedef FEbase<double, v_dfes<Mesh>>* pdfRbase; typedef std::pair<pdfRbase, int> pdfR;
    typedef FEbase<Complex, v_dfes<Mesh>>* pdfCbase; typedef std::pair<pdfCbase, int> pdfC;

    TheOperators->Add("=",
        new OneOperator2_<pdfR, pdfR, pdfR>(setFEDistributed<Mesh, double>),
        new OneOperator2_<pdfC, pdfC, pdfC>(setFEDistributed<Mesh, Complex>));
    TheOperators->Add("<-",
        new OneOperator3_<pdfRbase*, pdfRbase*, v_dfes<Mesh>**, pdfR>(initFEDistributed<Mesh, double>),
        new OneOperator3_<pdfCbase*, pdfCbase*, v_dfes<Mesh>**, pdfC>(initFEDistributed<Mesh, Complex>));
}

template void registerTransferInterpolateOps<Mesh3>();
template void registerTransferInterpolateOps<MeshS>();
template void registerTransferInterpolateOps<MeshL>();
