#include "ff++.hpp"
#include "metis.h"

#ifdef FF_WITH_SCOTCH
#include "PartitionScotch.hpp"
#endif

#ifdef PARALLELE
#include <mpi.h>
#endif

#if defined(PARALLELE) && defined(FF_WITH_PARMETIS)
#include "PartitionParmetis.hpp"
#endif


extern long mpisize, mpirank;
static const size_t CHUNK = size_t(512) << 20;
static const int TAG_SCATTER_HDR = 2000;
static const int TAG_SCATTER_BUF = 2001;
static const int TAG_SCATTER_PART = 2002;
static const int TAG_POU_CHECK = 1003;
static const int TAG_SYM_CHECK = 1004; 
static const int TAG_DOF_POU = 1001;
static const int TAG_DOF_NUM = 1002;

template<class Mesh>
int computeGlobalPartition(const Mesh &Th, KN<int> &part, const std::string &method, pcommworld comm, bool broadcast){
    const int nbt = Th.nt, nbv = Th.nv;
    const int nve = Mesh::Element::nv;

    // Séquentiel
    if (mpisize <= 1){
        part = 0;
        return DIST_OK;
    }

    int status = checkPartitionMethod(method);
    if (!status && mpisize >= nbt) status = DIST_NT_TOO_SMALL;
    if (status) return status;
    
    if (method == "parmetis") {
#if defined(PARALLELE) && defined(FF_WITH_PARMETIS)
      MPI_Comm cwp = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
      if (ffParmetisPartFaceDual(Th, (int)mpisize, (int*)part, cwp) != 0)
        status = DIST_PART_FAILED;
      return status;                 // tous les rangs ont part : pas de broadcast
#endif
  }

    if (mpirank == 0) {
      if (method == "metis") {
        idx_t nt = nbt, nv = nbv;
        KN<idx_t> eptr(nt + 1), elmnts(nve * nt), npart(nv), epart(nt, 0);
        for (idx_t k = 0, i = 0; k < nt; ++k) {
          eptr[k] = i;
          for (idx_t j = 0; j < nve; j++)
            elmnts[i++] = Th(k, j);
          eptr[k + 1] = i;
        }
        idx_t nparts = mpisize, edgecut, ncommon = 1;
        METIS_PartMeshDual(&nt, &nv, eptr, (idx_t *)elmnts, 0, 0, &ncommon, &nparts, 0, 0,
                          &edgecut, (idx_t *)epart, (idx_t *)npart);
        for (int k = 0; k < nbt; ++k) part[k] = (int)epart[k];
      }
      else {  // "scotch"
#ifdef FF_WITH_SCOTCH
        KN<SCOTCH_Num> epart(nbt);
        SCOTCH_randomReset();
        if (ffScotchPartFaceDual(Th, (SCOTCH_Num)mpisize, (SCOTCH_Num *)epart) != 0){
          status = DIST_PART_FAILED;
        }
        else {
          for (int k = 0; k < nbt; ++k) part[k] = (int)epart[k];
        }
#endif
      }
    }

  // --- Broadcast depuis le rang 0 (uniquement en build MPI) ---
#ifdef PARALLELE
  if (broadcast){
    MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
    MPI_Bcast((int*)part, nbt, MPI_INT, 0, cw);
  }
#endif
  return status;
}

const char* distributeStatusMessage(int status) {
  switch (status) {
    case DIST_PART_SIZE: return "distribute: partition[] requires Th.nt values";
    case DIST_PART_RANGE : return "distribute: partition[] has a value outside [0, mpisize[";
    case DIST_PART_EMPTY : return "distribute: at least one rank gets no element: reduce mpisize or refine the mesh";
    case DIST_NT_TOO_SMALL : return "distribute: mpisize must be < Th.nt";
    case DIST_PART_FAILED : return "distribute: partitioning failure";
    case DIST_METHOD_BAD : return "distribute: unknown partmethod (expected \"metis\", \"scotch\" or \"parmetis\")";
    case DIST_SCOTCH_NA : return "distribute: scotch partitioner is not available in this build";
    case DIST_PARMETIS_NA : return "distribute: parmetis partitioner is not available in this build";
    case DIST_PARMETIS_SCAT : return "distribute: parmetis is not compatible with scatter mode";
    case DIST_ARG_METHOD : return "distribute: partmethod differs between ranks";
    case DIST_ARG_OVERLAP : return "distribute: overlap differs between ranks";
    case DIST_ARG_KEEPGLOBAL : return "distribute: keepGlobal differs between ranks";
    case DIST_ARG_SCATTER : return "distribute: scatter argument does not agree with detected mode";
  }
  return "distribute: unknown status";
}

int checkPartitionMethod(const std::string &method) {
  if (method != "metis" && method != "scotch" && method != "parmetis") return DIST_METHOD_BAD;
#ifndef FF_WITH_SCOTCH
  if (method == "scotch") return DIST_SCOTCH_NA;
#endif
#if !defined(PARALLELE) || !defined(FF_WITH_PARMETIS)
  if (method == "parmetis") return DIST_PARMETIS_NA;
#endif
  return DIST_OK;
}

int checkPartitionUsable(const KN<int>& part, int nt, bool allowEmpty){
  if (part.n != nt) return DIST_PART_SIZE;
  if (mpisize <= 1) return DIST_OK;
  KN<int> count((int)mpisize, 0);
  for (int k = 0; k < nt; ++k){
    const int p = part[k];
    if (p < 0 || p >= (int)mpisize) return DIST_PART_RANGE;
    count[p]++;
  }
  if (!allowEmpty){
    for (int r = 0; r < (int)mpisize; ++r) if (count[r] == 0) return DIST_PART_EMPTY;
  }
  return DIST_OK;
}

#ifdef PARALLELE
template<class T, class Vec>
static void exchangeOnIntersections(MPI_Comm cw, const KN<int>& nbr, const KN<KN<long>>& dofI, const Vec& v, MPI_Datatype dt, int tag, std::vector<std::vector<T>>& rcv)
{
  const int nN = nbr.n;
  std::vector<std::vector<T>> snd(nN);
  rcv.assign(nN, std::vector<T>());
  std::vector<MPI_Request> rq(2*nN);
  for (int j = 0; j < nN; ++j){
    const KN<long>& L = dofI[1+j];
    snd[j].resize(L.n);
    rcv[j].resize(L.n);
    for (long p = 0; p < L.n; ++p) snd[j][p] = v[L[p]];
    MPI_Irecv(rcv[j].data(), (int)L.n, dt, nbr[j], tag, cw, &rq[j]);
    MPI_Isend(snd[j].data(), (int)L.n, dt, nbr[j], tag, cw, &rq[nN+j]);
  }
  MPI_Waitall(2*nN, rq.data(), MPI_STATUSES_IGNORE);
}

static double pouResidualLocal(const KN<KN<long>>& dofI, const KN<double>& Ddof, int nLocDof, const std::vector<std::vector<double>>& rcvD){
  std::vector<double> sum(nLocDof);
  const int nN = dofI.n - 1;
  for (int d = 0; d < nLocDof; ++d) sum[d] = Ddof[d];
  for (int j = 0; j < nN; ++j){
    const KN<long>& L = dofI[1+j];
    for (long p = 0; p < L.n; ++p) sum[L[p]] += rcvD[j][p];
  }
  double e = 0.0;
  for (int d = 0; d < nLocDof; ++d) e = std::max(e, fabs(sum[d] - 1.0));
  return e;
}

int detectDistributionMode(int localNt, pcommworld comm)
{
  MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
  int size;
  MPI_Comm_size(cw, &size);

  if (size != (int)mpisize){
    ExecError("distribute: sub-communicator not yet supported hence comm must cover all ranks");
  }

  int has = (localNt > 0) ? 1 : 0;
  KN<int> hasMesh(size);
  MPI_Allgather(&has, 1, MPI_INT, (int*)hasMesh, 1, MPI_INT, cw);

  int nHas = hasMesh.sum();
  if (nHas == size){
    return DM_REPLICATED;
  }
  else if (nHas == 1 && hasMesh[0] == 1){
    return DM_SCATTER;
  }
  else if (nHas == 0){
    ExecError("distribute: no rank has a mesh");
  }
  else if (hasMesh[0] == 0){
    ExecError("distribute: to use scatter mode, the global mesh must be built on rank 0");
  }
  else{
    ExecError("distribute: mixed mode not supported: either all ranks have the global mesh, or only rank 0");
  }
  return DM_REPLICATED;
}

int agreeOnStatus(int local, pcommworld comm){
  MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, cw);
  return global;
}

static bool hasAttached(const Mesh3& Th) { return Th.meshS != nullptr; }
static bool hasAttached(const MeshS&) { return false; }
static bool hasAttached(const MeshL&) { return false; }
static void rebuildAttached(Mesh3* Th, bool f) {if (f) Th->BuildMeshS();}
static void rebuildAttached(MeshS*, bool) {}
static void rebuildAttached(MeshL*, bool) {}

template<class Mesh>
void sendMesh(const Mesh& Th, int dest, pcommworld comm){
  MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
  Serialize buf = Th.serialize();
  long long hdr[2] = { (long long)buf.size(), hasAttached(Th) ? 1LL : 0LL};

  MPI_Send(hdr, 2, MPI_LONG_LONG, dest, TAG_SCATTER_HDR, cw);

  char* raw = (char*)buf;
  size_t total = (size_t)hdr[0], off = 0;
  while (off < total) {
    int n = (int)std::min(CHUNK, total - off);
    MPI_Send(raw+off, n, MPI_BYTE, dest, TAG_SCATTER_BUF, cw);
    off +=n;
  }
}

template<class Mesh>
Mesh* recvMesh(int src, pcommworld comm){
  MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
  MPI_Status st;
  long long hdr[2];
  MPI_Recv(hdr, 2, MPI_LONG_LONG, src, TAG_SCATTER_HDR, cw, &st);
  ffassert(hdr[0]>0);

  Serialize buf((size_t)hdr[0], Fem2D::GenericMesh_magicmesh);
  char* raw = (char*)buf;
  size_t total = (size_t)hdr[0], off = 0;

  while (off < total){
    int n = (int)std::min(CHUNK, total-off);
    MPI_Recv(raw+off, n, MPI_BYTE, src, TAG_SCATTER_BUF, cw, &st);
    int got;
    MPI_Get_count(&st, MPI_BYTE, &got);
    off +=got;
  }

  Mesh* Th = new Mesh(buf);
  rebuildAttached(Th, hdr[1] != 0);
  return Th;
}

void sendPartition(const KN<int>& part, int dest, pcommworld comm){
  MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
  MPI_Send((int*)part, part.n, MPI_INT, dest, TAG_SCATTER_PART, cw);
}

KN<int> recvPartition(int n, int src, pcommworld comm){
  MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
  KN<int> partition(n);
  MPI_Status st;
  MPI_Recv((int*)partition, n, MPI_INT, src, TAG_SCATTER_PART, cw, &st);
  return partition;
}
#else
int detectDistributionMode(int, pcommworld){
  return DM_REPLICATED;
}

int agreeOnStatus(int local, pcommworld){
  return local;
}

template<class Mesh>
void sendMesh(const Mesh&, int, pcommworld){
  ExecError("sendMesh is used only in a parallel setting.");
}
template<class Mesh>
Mesh* recvMesh(int, pcommworld){
  ExecError("recvMesh is used only in a parallel setting.");
  return nullptr;
}
void sendPartition(const KN<int>&, int, pcommworld){
  ExecError("sendPartition is only used in a parallel setting");
}
KN<int> recvPartition(int, int, pcommworld){
  ExecError("recvPartition is only used in a parallel setting");
  return KN<int>(0);
}
#endif

static KN<long> trivialNumbering(int nLocDof, long& globalNdof){
  KN<long> num(nLocDof);
  for (int d = 0; d < nLocDof; ++d) num[d] = d;
  globalNdof = nLocDof;
  return num;
}

#ifndef PARALLELE
KN<long> distributedDofNumbering(pcommworld, const KN<KN<long> > &, const KN<double> &, int nLocDof, long &globalNdof)
{
  return trivialNumbering(nLocDof, globalNdof);
}

double checkPartitionOfUnity(pcommworld, const KN<KN<long> >&, const KN<double>&, int){ return 0.0; }

int checkIntersectionSymmetry(pcommworld, const KN<KN<long>>&) { return -1; }

int checkPartitionConsistency(pcommworld, const KN<int> &) { return 0; }

#else
KN<long> distributedDofNumbering(pcommworld comm, const KN<KN<long>>& dofI, const KN<double>& Ddof, int nLocDof, long& globalNdof)
{
  if (mpisize <=1) return trivialNumbering(nLocDof, globalNdof);

  MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
  const int nN = dofI.n - 1;
  ffassert(dofI[0].n == nN);
  KN<int> nbr(nN);
  for (int j = 0; j < nN; ++j) nbr[j] = (int)dofI[0][j];
  int rank, size;
  MPI_Comm_rank(cw, &rank);
  MPI_Comm_size(cw, &size);

  std::vector<std::vector<double>> rcvD;
  exchangeOnIntersections<double>(cw, nbr, dofI, Ddof, MPI_DOUBLE, TAG_DOF_POU, rcvD);

  const double EPS = 1.0e-12; // HPDDM_EPS
  std::vector<char> owned(nLocDof, 1);

  for (int j = 0; j < nN; ++j){
    const KN<long>& L = dofI[1+j];
    for (long p = 0; p <L.n; ++p){
      const long d = L[p];
      const double v = Ddof[d] - rcvD[j][p];
      const bool win = (v > EPS) || (fabs(v) < EPS && rank > nbr[j]);
      if (!win) owned[d] = 0;
    }
  }

  long nOwned = 0;
  for (int d = 0; d < nLocDof; ++d) if (owned[d]) ++nOwned;

  std::vector<long> counts(size);
  MPI_Allgather(&nOwned, 1, MPI_LONG, counts.data(), 1, MPI_LONG, cw);

  long start = 0;
  for (int r = 0; r < rank; ++r) start += counts[r];
  globalNdof = 0;
  for (int r = 0; r < size; ++r) globalNdof += counts[r];

  KN<long> num(nLocDof, -1L);
  long c = start;
  for (int d = 0; d < nLocDof; ++d) if (owned[d]) num[d] = c++;

  std::vector<std::vector<long>> rcvN;
  exchangeOnIntersections<long>(cw, nbr, dofI, num, MPI_LONG, TAG_DOF_NUM, rcvN);

  for (int j = 0; j < nN; ++j){
    const KN<long>& L = dofI[1+j];
    for (long p = 0; p < L.n; ++p) {
      const long d = L[p];
      if (num[d] < 0 && rcvN[j][p] >= 0) num[d] = rcvN[j][p];
    }
  }
  if (nLocDof > 0) ffassert(num.min() >= 0);
  return num;
}

double checkPartitionOfUnity(pcommworld comm, const KN<KN<long>>& dofI, const KN<double>& Ddof, int nLocDof){
  if (mpisize <= 1) return 0.0;

  MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
  const int nN = dofI.n-1;
  ffassert(dofI[0].n == nN);
  ffassert(Ddof.n == nLocDof);
  KN<int> nbr(nN);
  for (int j = 0; j < nN; ++j) nbr[j] = (int)dofI[0][j];

  std::vector<std::vector<double>> rcvD;
  exchangeOnIntersections<double>(cw, nbr, dofI, Ddof, MPI_DOUBLE, TAG_POU_CHECK, rcvD);

  const double e = pouResidualLocal(dofI, Ddof, nLocDof, rcvD);

  double eg = 0.0;
  MPI_Allreduce(&e, &eg, 1, MPI_DOUBLE, MPI_MAX, cw);
  return eg;
}

// Return -1 if all sizes match, else the local index of the faulty neighbour
// Application before purge, to have all Irecv matched
int checkIntersectionSymmetry(pcommworld comm, const KN<KN<long>>& dofI) {
  if (mpisize <= 1) return -1;

  MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
  const int nN = dofI.n -1;
  ffassert(dofI[0].n == nN);

  std::vector<MPI_Request> rq(2*nN);
  std::vector<long> mine(nN), his(nN);
  for (int j = 0; j < nN; ++j) mine[j] = dofI[1+j].n;
  for (int j = 0; j < nN; ++j){
    const int nb = (int)dofI[0][j];
    MPI_Irecv(&his[j], 1, MPI_LONG, nb, TAG_SYM_CHECK, cw, &rq[j]);
    MPI_Isend(&mine[j], 1, MPI_LONG, nb, TAG_SYM_CHECK, cw, &rq[nN+j]);
  }
  MPI_Waitall(2*nN, rq.data(), MPI_STATUSES_IGNORE);

  for (int j = 0; j < nN; ++j) if (mine[j] != his[j]) return j;
  return -1;
}

int checkPartitionConsistency(pcommworld comm, const KN<int>& part){
  if (mpisize <= 1) return 0;

  MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(cw, &rank);

  int bad = (agreeOnStatus(part.n, comm) != part.n) ? 1 : 0;
  if (agreeOnStatus(bad, comm)) return 1;

  KN<int> ref(part.n);
  if (rank == 0) ref = part;
  MPI_Bcast((int*)ref, part.n, MPI_INT, 0, cw);
  for (int k = 0; k < part.n; ++k){
    if (ref[k] != part[k]) {bad = 1; break; }
  }

  return agreeOnStatus(bad, comm);
}

#endif

// Instanciations explicites (miroir de registerDistributedMeshOps<MeshS/Mesh3/MeshL>)
template int computeGlobalPartition<MeshS>(const MeshS&, KN<int>&, const std::string&, pcommworld, bool);
template int computeGlobalPartition<Mesh3>(const Mesh3&, KN<int>&, const std::string&, pcommworld, bool);
template int computeGlobalPartition<MeshL>(const MeshL&, KN<int>&, const std::string&, pcommworld, bool);

template void sendMesh<MeshS>(const MeshS&, int, pcommworld);
template MeshS* recvMesh<MeshS>(int, pcommworld);
template void sendMesh<Mesh3>(const Mesh3&, int, pcommworld);
template Mesh3* recvMesh<Mesh3>(int, pcommworld);
template void sendMesh<MeshL>(const MeshL&, int, pcommworld);
template MeshL* recvMesh<MeshL>(int, pcommworld);