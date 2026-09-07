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

static bool coversAllRows(const MatriceMorse<double>* M, int nrow) {
    KN<int> cnt(nrow, 0);
    KN<double> row(nrow, 0.0);
    for (size_t k = 0; k < M->nnz; ++k) { cnt[M->i[k]]++; row[M->i[k]] += M->aij[k]; }
    for (int i = 0; i < nrow; ++i)
        if (cnt[i] == 0 || std::abs(row[i] - 1.0) > 1e-12) return false;

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
static Mesh* buildFragmentFor(const Mesh& Loc, const BBox& target, KN<int>& n2o) {
  const int nv = Mesh::Element::nv;
  KN<int> mask(Loc.nt, 0);
  int nkept = 0;

  for (int k = 0; k < Loc.nt; ++k) {
    typename Mesh::Rd g;
    for (int i = 0; i < nv; ++i) g += Loc[k][i];
    g /= nv;
    if (pointInBox(target, g)) { mask[k] = 1; ++nkept; }
  }
  if (nkept == 0) { n2o = KN<int>(0); return nullptr; }

  Mesh* frag = truncmesh(Loc, 1, (int*)mask, false, 9999, 1e-7, 1L, true, false);
  ffassert(frag && frag->nt > 0);

  frag->BuildGTree();
  n2o = n2oFromSplit(Loc, *frag, mask);
  return frag;
}

template<class Mesh, class R>
static KN<R> fragmentDofVector(const GFESpace<Mesh>& srcVh, const KN<R>& scaledU,
                               const KN<double>& chi, const Mesh& frag, const KN<int>& n2o)
{
    ffassert(srcVh.TFE.N() == 1);            // espaces composites hors perimetre
    GFESpace<Mesh> fragVh(frag, *srcVh.TFE[0]);   
    KN<long> inj = restrictDOFPartial(fragVh, srcVh, n2o);  
    const int nd = fragVh.NbOfDF; 
    KN<R> sub(2*nd);
    for (int d = 0; d < nd; ++d) { 
        if (inj[d] < 0) { sub[d] = R(); sub[nd + d] = R(); }
        else { sub[d] = scaledU[inj[d]]; sub[nd + d] = R(chi[inj[d]]); }
    }
    return sub;
}

template<class Mesh, class R>
struct RecvFragment {
    Mesh*  mesh = nullptr;   
    KN<R>  dof;             
};

#ifdef PARALLELE
template<class Mesh, class R>
static void exchangeFragments(pcommworld comm,
                              const KN<int>& sendToRanks, const KN<int>& recvFromRanks,
                              const std::vector<Mesh*>& fragOut,      // aligne sur sendToRanks
                              const std::vector<KN<R> >& subOut,      // idem, .n == 0 si vide
                              std::vector<RecvFragment<Mesh,R> >& in) // aligne sur recvFromRanks
{
    MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
    const int nS = sendToRanks.n, nR = recvFromRanks.n;
    in.assign(nR, RecvFragment<Mesh,R>());

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
        }
        MPI_Isend(&hdrOut[2*i], 2, MPI_LONG_LONG, sendToRanks[i], TAG_XFER_HDR, cw, &rqH[nR+i]);
    }
    MPI_Waitall(nR + nS, rqH.data(), MPI_STATUSES_IGNORE);

    std::vector<Serialize*> bufIn(nR, nullptr);
    std::vector<MPI_Request> rq;
    rq.reserve(2*(nR + nS));

    for (int j = 0; j < nR; ++j) {
        if (hdrIn[2*j] == 0) continue;
        bufIn[j] = new Serialize((size_t)hdrIn[2*j], Fem2D::GenericMesh_magicmesh);
        in[j].dof.resize((long)hdrIn[2*j+1]);
        rq.resize(rq.size() + 2);
        MPI_Irecv((char*)*bufIn[j], (int)hdrIn[2*j], MPI_BYTE,
                  recvFromRanks[j], TAG_XFER_BODY, cw, &rq[rq.size()-2]);
        MPI_Irecv((R*)in[j].dof, (int)(hdrIn[2*j+1]*sizeof(R)), MPI_BYTE,
                  recvFromRanks[j], TAG_XFER_DOF, cw, &rq[rq.size()-1]);
    }
    for (int i = 0; i < nS; ++i) {
        if (!ser[i]) continue;
        rq.resize(rq.size() + 2);
        MPI_Isend((char*)*ser[i], (int)hdrOut[2*i], MPI_BYTE,
                  sendToRanks[i], TAG_XFER_BODY, cw, &rq[rq.size()-2]);
        MPI_Isend((R*)subOut[i], (int)(subOut[i].n*sizeof(R)), MPI_BYTE,
                  sendToRanks[i], TAG_XFER_DOF, cw, &rq[rq.size()-1]);
    }
    if (!rq.empty()) MPI_Waitall((int)rq.size(), rq.data(), MPI_STATUSES_IGNORE);

    // ---- reconstruction ----------------------------------------------
    for (int j = 0; j < nR; ++j) {
        if (!bufIn[j]) continue;
        in[j].mesh = new Mesh(*bufIn[j]);
        in[j].mesh->BuildGTree();
        delete bufIn[j];
    }
    for (int i = 0; i < nS; ++i) delete ser[i];
}
#else
template<class Mesh, class R>
static void exchangeFragments(pcommworld,
                              const KN<int>& sendToRanks, const KN<int>& recvFromRanks,
                              const std::vector<Mesh*>& fragOut,
                              const std::vector<KN<R> >& subOut,
                              std::vector<RecvFragment<Mesh,R> >& in)
{
    in.assign(recvFromRanks.n, RecvFragment<Mesh,R>());
    if (recvFromRanks.n == 0 || sendToRanks.n == 0) return;
    ffassert(recvFromRanks.n == 1 && sendToRanks.n == 1);
    if (!fragOut[0]) return;

    Serialize s = fragOut[0]->serialize();
    in[0].mesh = new Mesh(s);
    in[0].mesh->BuildGTree();
    in[0].dof.resize(subOut[0].n);
    in[0].dof = subOut[0];
}
#endif

template<class Mesh1, class Mesh2, class R>
static void collectFragments(const DistributedMesh<Mesh1>& Dsrc,
                             const GFESpace<Mesh1>& srcVh, const KN<R>& srcU,
                             const DistributedMesh<Mesh2>& Ddst,
                             KN<int>& sendToRanks, KN<int>& recvFromRanks,
                             std::vector<RecvFragment<Mesh1,R> >& in, KN<long>* sentCounts = nullptr)
{
    std::vector<BBox> allSrc, allDst;
    computeOverlapRankPairs(Dsrc.comm, Dsrc, Ddst, sendToRanks, recvFromRanks, allSrc, allDst);
    const Mesh1* srcGeom = Dsrc.CoverMesh ? Dsrc.CoverMesh : Dsrc.LocalMesh;
    KN<int> cover2local(srcGeom->nt, -1);
    if (Dsrc.CoverMesh){
        for (int k = 0; k < Dsrc.localToCoverElement.n; ++k)
            cover2local[Dsrc.localToCoverElement[k]] = k;
    }
    else{
        for (int k = 0; k < srcGeom->nt; ++k) cover2local[k] = k;
    }

    if (sentCounts) sentCounts->resize(sendToRanks.n);
    KN<double> chi = interpolatePoU(Dsrc, srcVh);   
    ffassert(chi.n == srcVh.NbOfDF && srcU.n == srcVh.NbOfDF);
    KN<R> scaledU(srcU.n);
    for (int d = 0; d < srcU.n; ++d) scaledU[d] = srcU[d] * chi[d];

    std::vector<Mesh1*> fragOut(sendToRanks.n, nullptr);
    std::vector<KN<R> > subOut(sendToRanks.n);
    for (int i = 0; i < sendToRanks.n; ++i) {
        KN<int> n2oCover;
        fragOut[i] = buildFragmentFor(*srcGeom, allDst[sendToRanks[i]], n2oCover);
        if (sentCounts) (*sentCounts)[i] = fragOut[i] ? fragOut[i]->nt : 0;
        if (fragOut[i]) {
            KN<int> n2oLocal(n2oCover.n);
            for (int kc = 0; kc < n2oCover.n; ++kc){
                n2oLocal[kc] = cover2local[n2oCover[kc]];
            }
            subOut[i] = fragmentDofVector(srcVh, scaledU, chi, *fragOut[i], n2oLocal);
        }
    }

    exchangeFragments(Dsrc.comm, sendToRanks, recvFromRanks, fragOut, subOut, in);

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

 template<class Mesh, class R>
 static void accumulateFragment(const GFESpace<Mesh>& dstVh, const GFESpace<Mesh>& fragVh, const KN<R>& payload, const int* data, KN<R>& dstU, KN<R>& cover) {
    const int nd = fragVh.NbOfDF;
    ffassert(payload.n == 2*nd);

    MatriceMorse<double>* M =
        buildInterpolationMatrixT<GFESpace<Mesh>, GFESpace<Mesh> >(dstVh, fragVh, (int*)data);

    for (size_t k = 0; k < M->nnz; ++k) {
        const int ii = M->i[k], cc = M->j[k];
        dstU[ii] += M->aij[k] * payload[cc];
        cover[ii] += M->aij[k] * payload[nd+cc];
    }
    delete M;
 }

 enum TransferPath { XFER_GENERAL = 0, XFER_SINGLE_RANK = 1, XFER_LOCAL = 2 };

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

    const int N = dstVh.N;
    KN<int> data = dataInterpolate(N);
    dstU = R();
    KN<R> cover(dstU.n, R());
    const long saveSearch = searchMethod;
    searchMethod = 0;

    const int nproc = mpisize > 0 ? (int)mpisize : 1;
    if (nproc == 1) {
        KN<double> chi = interpolatePoU(Dsrc, srcVh);
        KN<R> payload(2*srcVh.NbOfDF);
        for (int d = 0; d < srcVh.NbOfDF; ++d){
            payload[d] = srcU[d]*chi[d];
            payload[srcVh.NbOfDF + d] = R(chi[d]);
        }
        accumulateFragment(dstVh, srcVh, payload, data, dstU, cover);
        searchMethod = saveSearch;
        for (int i = 0; i < dstU.n; ++i){
            if (std::abs(cover[i]) > 1e-14) dstU[i] /= cover[i];
        }
        return XFER_SINGLE_RANK;
    }

    MatriceMorse<double>* Mloc = nullptr;
    int localOK = 0;
    {
        const BBox bs = rawBoxOf(srcVh.Th), bd = rawBoxOf(dstVh.Th);
        if (boxContains(bs, bd)) {
            Mloc = buildInterpolationMatrixT<GFESpace<Mesh>, GFESpace<Mesh> >(dstVh, srcVh, (int*)data);
            localOK = coversAllRows(Mloc, dstU.n) ? 1 : 0;
        }
    }

    int globalOK = localOK;
    #ifdef PARALLELE
    {
        MPI_Comm cw = Dsrc.comm ? *(MPI_Comm*)Dsrc.comm : MPI_COMM_WORLD;
        MPI_Allreduce(&localOK, &globalOK, 1, MPI_INT, MPI_MIN, cw);
    }
    #endif

    if (globalOK) {
        for (size_t k = 0; k < Mloc->nnz; ++k)
            dstU[Mloc->i[k]] += Mloc->aij[k]*srcU[Mloc->j[k]];
        delete Mloc;
        searchMethod = saveSearch;
        return XFER_LOCAL;
    }
    delete Mloc;

    KN<int> snd, rcv;
    std::vector<RecvFragment<Mesh,R> > in;
    collectFragments(Dsrc, srcVh, srcU, Ddst, snd, rcv, in);

    for (size_t j = 0; j < in.size(); ++j) {
        if (!in[j].mesh) continue;
        GFESpace<Mesh> fragVh(*in[j].mesh, *srcVh.TFE[0]);
        accumulateFragment(dstVh, fragVh, in[j].dof, data, dstU, cover);
    }
    for (int i = 0; i < dstU.n; ++i)
        if (std::abs(cover[i]) > 1e-14) dstU[i] /= cover[i];

    searchMethod = saveSearch;

    for (size_t j = 0; j < in.size(); ++j)
        if (in[j].mesh) in[j].mesh->destroy();

    return XFER_GENERAL;
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
    KN<double> u(srcVh.NbOfDF, 1.0);

    KN<int> snd, rcv;
    std::vector<RecvFragment<Mesh,double> > in;
    collectFragments(A, srcVh, u, B, snd, rcv, in, nEnvoyes);

    nRecus->resize((long)in.size());
    for (int j = 0; j < (int)in.size(); ++j) {
        (*nRecus)[j] = in[j].mesh ? in[j].mesh->nt : 0;
        if (in[j].mesh) in[j].mesh->destroy();
    }
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
long interpolateD(v_dfes<Mesh>** const & ppSrc, KN<R>* const & uSrc,
                  v_dfes<Mesh>** const & ppDst, KN<R>* const & uDst)
{
    throwassert(ppSrc && *ppSrc && ppDst && *ppDst && uSrc && uDst);
    v_dfes<Mesh>* pSrc = *ppSrc;
    v_dfes<Mesh>* pDst = *ppDst;
    throwassert(pSrc->DTh && pDst->DTh);
    ffassert(pSrc->DTh->comm == pDst->DTh->comm);   // meme limite que detectDistributionMode

    const GFESpace<Mesh>& srcVh = **pSrc;           // operator FESpace*() : lgfem.hpp:427
    const GFESpace<Mesh>& dstVh = **pDst;

    ffassert(uSrc->n == srcVh.NbOfDF);
    ffassert(uDst->n == dstVh.NbOfDF);

    int xfer = interpolateDistributed(*pSrc->DTh, srcVh, *uSrc, *pDst->DTh, dstVh, *uDst);
    return (long)xfer;
}

template<class Mesh>
void registerTransferInterpolateOps() {
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


}

template void registerTransferInterpolateOps<Mesh3>();
template void registerTransferInterpolateOps<MeshS>();
template void registerTransferInterpolateOps<MeshL>();
