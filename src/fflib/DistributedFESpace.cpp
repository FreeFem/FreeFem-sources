#include "ff++.hpp"

// Utility functions for building distributed dof data
namespace {
template<class Mesh>
KN<long> buildInterfaceDofList(const DistributedMesh<Mesh>& D, int j, const GFESpace<Mesh>&Vhi){
    KN<int> n2oS;
    Mesh* S = buildInterfaceSubmesh(D, j, n2oS);
    if (!S) return KN<long>(0);

    const int me = D.rank;
    const int nb = D.neighborRanks[j];

    KN<int> coverToLocal(D.CoverMesh->nt, -1);
    for (int kl = 0; kl < D.localToCoverElement.n; ++kl){
        coverToLocal[D.localToCoverElement[kl]] = kl;
    }

    KN<long> out;
    {
        GFESpace<Mesh> Vs(*S, *Vhi.TFE[0]);
        KN<int> byI(Vs.NbOfDF, 0), byJ(Vs.NbOfDF, 0);
        for (int kS = 0; kS < Vs.NbOfElements; ++kS){
            const int p = D.coverPartition[n2oS[kS]];
            GFElement<Mesh> KS(&Vs, kS);
            for (int df = 0; df < KS.NbDoF(); ++df){
                const int d = KS(df);
                if (p == me) byI[d] = 1;
                else if (p == nb) byJ[d] = 1;
            }
        }
        KN<int> n2oLocal(Vs.NbOfElements);
        for (int kS = 0; kS < Vs.NbOfElements; ++kS){
            n2oLocal[kS] = coverToLocal[n2oS[kS]];
        }
        KN<long> mapS2L = restrictDOFPartial(Vs, Vhi, n2oLocal);

        int n = 0;
        for (int d = 0; d < Vs.NbOfDF; ++d) if (byI[d] && byJ[d]) ++n;
        out.resize(n);
        int c = 0;
        for (int d = 0; d < Vs.NbOfDF; ++d){
            if (byI[d] && byJ[d]) { ffassert(mapS2L[d] >= 0); out[c++] = mapS2L[d]; }
        }
    }
    S->destroy();
    return out;
}

template<class Mesh>
KN<KN<long>> buildDofIntersection(const DistributedMesh<Mesh>& D, const GFESpace<Mesh>& Vhi){
    const int nN = D.neighborRanks.n;
    KN<KN<long>> dofI(1 + nN);
    dofI[0].resize(nN);

    for (int j = 0; j< nN; ++j) dofI[0][j] = D.neighborRanks[j];

    for (int j = 0; j < nN; ++j){
        if (D.overlap == 0) {
            dofI[1+j] = buildInterfaceDofList(D, j, Vhi);
            continue;
        }
        KN<int> n2o;
        Mesh* WhMesh = buildIntersectionSubmesh(D, j, n2o);
        {
            GFESpace<Mesh> Wh(*WhMesh, *Vhi.TFE[0]);
            dofI[1+j] = restrictDOF(Wh, Vhi, n2o);
        }
        WhMesh->destroy();
    }
    return dofI;
}

static KN<KN<long>> purgeEmptyIntersections(const KN<KN<long>>& raw) {
    const int nN = raw.n - 1;
    int nK = 0;
    for (int j = 0; j < nN; ++j) if (raw[1+j].n > 0) ++nK;
    KN<KN<long>> dofI(1+nK);
    dofI[0].resize(nK);
    int k = 0;
    for (int j = 0; j < nN; ++j){
        if (raw[1+j].n > 0) { 
            dofI[0][k] = raw[0][j];
            dofI[1+k] = raw[1+j];
            ++k;
        }
    }
    return dofI;
}

template<class Mesh>
KN<long> buildLoc2globDof(const DistributedMesh<Mesh>& D, const GFESpace<Mesh>& Vhi){
    ffassert(D.TrueGlobalMesh && D.localToGlobalElement.n == D.LocalMesh->nt);
    GFESpace<Mesh> Vhglob(*D.TrueGlobalMesh, *Vhi.TFE[0]);
    return restrictDOF(Vhi, Vhglob, D.localToGlobalElement);
}
}

// buildDistributedDofData() needed for the constructor
template<class Mesh>
void DGFESpace<Mesh>::buildDistributedDofData(){
    const DistributedMesh<Mesh>& D = DTh();
    const bool collective = (checkDfespace != 0) && (D.size > 1);

    KN<KN<long>> raw = buildDofIntersection(D, *this);

    if (D.overlap == 0){
        Ddof.resize(this->NbOfDF);
        Ddof = 1.0;
        for (int j = 0; j < raw.n -1; ++j){
            for (long p = 0; p < raw[1+j].n; ++p) Ddof[raw[1+j][p]] += 1.0;
        }
        for (int d = 0; d < this->NbOfDF; ++d) Ddof[d] = 1.0 / Ddof[d];
    }
    else {
        Ddof = interpolatePoU(D, *this);
    }

    if (collective){ 
        const int bad = checkIntersectionSymmetry(D.comm, raw);
        if (bad >= 0){
            cerr << "fespace distribute: ddl intersection asymmetric on the rank" << raw[0][bad] << " (local " << raw[1+bad].n << ")" << endl;
        }
        if (agreeOnStatus(bad >= 0 ? 1 : 0, D.comm)){
            distributedError(D.comm, "fespace distribute: incoherent intersection map");
        }
    }
    dofIntersectionDof = purgeEmptyIntersections(raw);

    pouResidual = (D.size <= 1) ? 0.0 : -1.0;
    if (collective){
        pouResidual = checkPartitionOfUnity(D.comm, dofIntersectionDof, Ddof, this->NbOfDF);
        if (pouResidual > 1.0e-8) {
            if (D.rank == 0){
                cout << "[fespace distribute] Problem with the partition of unity: max|sum -1| = " << pouResidual << endl;
            }
        }
    }
}

template<class Mesh>
void DGFESpace<Mesh>::ensureLoc2glob() const {
    if (loc2globBuilt) return;
    const DistributedMesh<Mesh>& D = DTh();
    if (!D.TrueGlobalMesh){
        distributedError(D.comm, "Udh.loc2glob : unavailable: no global mesh in scatter mode or distribute(..., keepGlobal=false)");
    }
    loc2globDof = buildLoc2globDof(D, *this);
    loc2globBuilt = true;
}

template<class Mesh>
void DGFESpace<Mesh>::ensureNumbering() const {
    if (numberingBuilt) return;
    const DistributedMesh<Mesh>& D = DTh();
    numberingDof = distributedDofNumbering(D.comm, dofIntersectionDof, Ddof, this->NbOfDF, globalNdof);
    numberingBuilt = true;
}

// constructor & destructor
template<class Mesh>
DGFESpace<Mesh>::DGFESpace(const DistributedMesh<Mesh>& D, const GTypeOfFE<Mesh>& tfe, int nbequibe, int* equibe)
    : GFESpace<Mesh>(*D.LocalMesh, tfe, nbequibe, equibe), dmesh(D) {
    buildDistributedDofData();
}

template<class Mesh>
DGFESpace<Mesh>::~DGFESpace() {
    if (exchangeHandle) g_distributedBackend->destroyExchangeHandle(exchangeHandle);
}

// dist() and factory functions to make the connection to lgfem.hpp
template<class Mesh>
const DGFESpace<Mesh>& dist(const GFESpace<Mesh>& Vh) {
    return static_cast<const DGFESpace<Mesh>&>(Vh);
}

template<class Mesh>
GFESpace<Mesh>* newDistributedFESpace(const DistributedMesh<Mesh>& DTh, const GTypeOfFE<Mesh>& tfe, int nbequibe, int* equibe) {
    return new DGFESpace<Mesh>(DTh, tfe, nbequibe, equibe);
}

// Instanciations
template class DGFESpace<Mesh3>;
template const DGFESpace<Mesh3>& dist(const GFESpace<Mesh3>& Vh);
template GFESpace<Mesh3>* newDistributedFESpace(const DistributedMesh<Mesh3>& DTh, const GTypeOfFE<Mesh3>& tfe, int nbequibe, int* equibe);

template class DGFESpace<MeshS>;
template const DGFESpace<MeshS>& dist(const GFESpace<MeshS>& Vh);
template GFESpace<MeshS>* newDistributedFESpace(const DistributedMesh<MeshS>& DTh, const GTypeOfFE<MeshS>& tfe, int nbequibe, int* equibe);

template class DGFESpace<MeshL>;
template const DGFESpace<MeshL>& dist(const GFESpace<MeshL>& Vh);
template GFESpace<MeshL>* newDistributedFESpace(const DistributedMesh<MeshL>& DTh, const GTypeOfFE<MeshL>& tfe, int nbequibe, int* equibe);

