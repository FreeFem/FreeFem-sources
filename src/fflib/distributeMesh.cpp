#include "ff++.hpp"
#include "metis.h"

#ifdef FF_WITH_SCOTCH
#include "PartitionScotch.hpp"
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

    if (method != "metis" && method != "scotch")
      ExecError("distribute: partitioner inconnu (attendu \"metis\" ou \"scotch\")");
    if (mpisize >= nbt)
      ExecError("distribute: mpisize doit etre strictement inferieur a Th.nt");

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
        ExecError("distribute: partitioner=\"scotch\" indisponible ; "
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

// Instanciations explicites (miroir de registerDistributedMeshOps<MeshS/Mesh3/MeshL>)
template void computeGlobalPartition<MeshS>(const MeshS&, KN<int>&, const std::string&, pcommworld);
template void computeGlobalPartition<Mesh3>(const Mesh3&, KN<int>&, const std::string&, pcommworld);
template void computeGlobalPartition<MeshL>(const MeshL&, KN<int>&, const std::string&, pcommworld);