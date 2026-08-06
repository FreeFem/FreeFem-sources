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

template<class Mesh>
void computeGlobalPartition(const Mesh &Th, KN<int> &part, const std::string &method, pcommworld comm){
    const int nbt = Th.nt, nbv = Th.nv;
    const int nve = Mesh::Element::nv;

    // Séquentiel
    if (mpisize <= 1){
        part = 0;
        return;
    }

    if (method != "metis" && method != "scotch" && method != "parmetis")
      ExecError("distribute: partmethod inconnu (attendu \"metis\", \"scotch\" ou \"parmetis\")");
    if (mpisize >= nbt)
      ExecError("distribute: mpisize doit etre strictement inferieur a Th.nt");

      if (method == "parmetis") {
#if defined(PARALLELE) && defined(FF_WITH_PARMETIS)
        MPI_Comm cwp = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
        if (ffParmetisPartFaceDual(Th, (int)mpisize, (int*)part, cwp) != 0)
          ExecError("distribute: echec du partitionnement ParMETIS");
        return;                 // tous les rangs ont part : pas de broadcast
#else
        ExecError("distribute: partmethod=\"parmetis\" indisponible ; reconfigurer "
                "FreeFEM avec ParMETIS (cf. PARMETIS_CORE dans la sortie de configure)");
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
        if (ffScotchPartFaceDual(Th, (SCOTCH_Num)mpisize, (SCOTCH_Num *)epart) != 0)
          ExecError("distribute: echec du partitionnement Scotch");
        for (int k = 0; k < nbt; ++k) part[k] = (int)epart[k];
        #else
        ExecError("distribute: partmethod=\"scotch\" indisponible ; "
                  "reconfigurer FreeFEM avec Scotch (cf. SCOTCH_CORE dans la sortie de configure)");
        #endif
      }
    }

  // --- Broadcast depuis le rang 0 (uniquement en build MPI) ---
#ifdef PARALLELE
  MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
  MPI_Bcast((int*)part, nbt, MPI_INT, 0, cw);
#endif
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
#endif

static KN<long> trivialNumbering(int nLocDof, long& globalNdof){
  KN<long> num(nLocDof);
  for (int d = 0; d < nLocDof; ++d) num[d] = d;
  globalNdof = nLocDof;
  return num;
}

#ifndef PARALLELE
KN<long> distributedDofNumbering(pcommworld, const KN<int>&, const KN<KN<long> > &, const KN<double> &, int nLocDof, long &globalNdof)
{
  return trivialNumbering(nLocDof, globalNdof);
}
#else
KN<long> distributedDofNumbering(pcommworld comm, const KN<int>& neighborRanks, const KN<KN<long>>& dofI, const KN<double>& Ddof, int nLocDof, long& globalNdof)
{
  if (mpisize <=1) return trivialNumbering(nLocDof, globalNdof);

  MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
  const int nN = neighborRanks.n;
  ffassert(dofI.n == 1 + nN);
  int rank, size;
  MPI_Comm_rank(cw, &rank);
  MPI_Comm_size(cw, &size);

  std::vector<std::vector<double>> rcvD;
  exchangeOnIntersections<double>(cw, neighborRanks, dofI, Ddof, MPI_DOUBLE, 1001, rcvD);

  const double EPS = 1.0e-12; // HPDDM_EPS
  std::vector<char> owned(nLocDof, 1);

  for (int j = 0; j < nN; ++j){
    const KN<long>& L = dofI[1+j];
    for (long p = 0; p <L.n; ++p){
      const long d = L[p];
      const double v = Ddof[d] - rcvD[j][p];
      const bool win = (v > EPS) || (fabs(v) < EPS && rank > neighborRanks[j]);
      if (!win) owned[d] = 0;
    }
  }

  if (verbosity > 4) {
    std::vector<double> sum(nLocDof);
    for (int d = 0; d < nLocDof; ++d) sum[d] = Ddof[d];
    for (int j = 0; j < nN; ++j){
      const KN<long>& L = dofI[1+j];
      for (long p = 0; p < L.n; ++p) sum[L[p]] += rcvD[j][p];
    }
    double e = 0.0;
    for (int d = 0; d < nLocDof; ++d) e = std::max(e, fabs(sum[d] - 1.0));
    double eg = 0.0;
    MPI_Reduce(&e, &eg, 1, MPI_DOUBLE, MPI_MAX, 0, cw);
    if (rank == 0) cout << "--distributedDofNumbering: max|sum PoU - 1| = " << eg << endl;
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
  exchangeOnIntersections<long>(cw, neighborRanks, dofI, num, MPI_LONG, 1002, rcvN);

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
#endif

// Instanciations explicites (miroir de registerDistributedMeshOps<MeshS/Mesh3/MeshL>)
template void computeGlobalPartition<MeshS>(const MeshS&, KN<int>&, const std::string&, pcommworld);
template void computeGlobalPartition<Mesh3>(const Mesh3&, KN<int>&, const std::string&, pcommworld);
template void computeGlobalPartition<MeshL>(const MeshL&, KN<int>&, const std::string&, pcommworld);