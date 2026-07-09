#include "ff++.hpp"
#include "metis.h"

#ifdef FF_WITH_SCOTCHMETIS
extern "C" {
    int  SCOTCH_METIS_PartMeshDual(idx_t*, idx_t*, idx_t*, idx_t*, idx_t*, idx_t*, idx_t*,
                                 idx_t*, real_t*, idx_t*, idx_t*, idx_t*, idx_t*);
    void SCOTCH_randomReset(void);
}
#endif

#ifdef PARALLELE
#include <mpi.h>
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

    // Partition computed on rank 0
    if (mpirank == 0) {
    idx_t nt = nbt, nv = nbv;
    KN<idx_t> eptr(nt + 1), elmnts(nve * nt), npart(nv), epart(nt, 0);
    for (idx_t k = 0, i = 0; k < nt; ++k) {
      eptr[k] = i;
      for (idx_t j = 0; j < nve; j++)
        elmnts[i++] = Th(k, j);
      eptr[k + 1] = i;
    }
    idx_t nparts = mpisize, edgecut, ncommon = 1;

    if (method == "metis") {
      METIS_PartMeshDual(&nt, &nv, eptr, (idx_t *)elmnts, 0, 0, &ncommon, &nparts, 0, 0,
                         &edgecut, (idx_t *)epart, (idx_t *)npart);
    }
    else if (method == "scotch") {
      #ifdef FF_WITH_SCOTCHMETIS
        SCOTCH_randomReset();
        SCOTCH_METIS_PartMeshDual(&nt, &nv, eptr, (idx_t *)elmnts, 0, 0, &ncommon, &nparts, 0, 0,
                         &edgecut, (idx_t *)epart, (idx_t *)npart);
      #else
        ExecError("distribute: partitioner=\"scotch\" indisponible");
      #endif
    }
    else {
      ExecError("distribute: partitioner inconnu (attendu \"metis\" ou \"scotch\")");
    }
    for (int k = 0; k < nbt; ++k) part[k] = (int)epart[k];
  }

  // --- Broadcast depuis le rang 0 (uniquement en build MPI) ---
#ifdef PARALLELE
  MPI_Comm cw = comm ? *(MPI_Comm*)comm : MPI_COMM_WORLD;
  MPI_Bcast((int*)part, nbt, MPI_INT, 0, cw);
#endif
}

// Instanciations explicites (miroir de registerDistributedMeshOps<MeshS/Mesh3/MeshL>)
template void computeGlobalPartition<MeshS>(const MeshS&, KN<int>&, const std::string&, pcommworld);
template void computeGlobalPartition<Mesh3>(const Mesh3&, KN<int>&, const std::string&, pcommworld);
template void computeGlobalPartition<MeshL>(const MeshL&, KN<int>&, const std::string&, pcommworld);