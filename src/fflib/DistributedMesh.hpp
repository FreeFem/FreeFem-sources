#include <vector>
#include "RNM.hpp"
#include <string>
#ifndef FF_PCOMMWORLD_TYPEDEF
#define FF_PCOMMWORLD_TYPEDEF
typedef void* pcommworld;
#endif

enum DistributionMode {DM_SCATTER = 0, DM_REPLICATED = 1 };

template <class Mesh>
class DistributedMesh : public RefCounter {
public: 
  Mesh * LocalMesh;
  Mesh* BorderMesh;
  const Mesh* CoverMesh = nullptr; // referent for truncated LocalMesh
  const Mesh* TrueGlobalMesh = nullptr; // true global mesh; null in scatter mode

  int overlap;
  int interfaceLabel;

  KN<int> neighborRanks;
  KN<double> partitionOfUnity;
  KN<int> localToCoverElement; // n2o : LocalMesh -> CoverMesh
  KN<int> localToGlobalElement; // n2o (LocalMesh -> Th)
  KN<int> coverPartition; // parititon restreinte à CoverMesh

  pcommworld comm = nullptr;
  
  DistributedMesh() : LocalMesh(nullptr), BorderMesh(nullptr), overlap(0), interfaceLabel(10),
  neighborRanks(), partitionOfUnity(), coverPartition() {}
//  DistributedMesh(Mesh * locmesh) : LocalMesh(locmesh) {}
  void setReferenceMeshes(const Mesh* cover, const Mesh* global){
    if (cover) cover->increment();
    if (global) global->increment();
    if (CoverMesh) CoverMesh->destroy();
    if (TrueGlobalMesh) TrueGlobalMesh->destroy();
    CoverMesh = cover;
    TrueGlobalMesh = global;
  }
  ~DistributedMesh() {
    if (LocalMesh) LocalMesh->destroy();
    if (BorderMesh) BorderMesh->destroy();
    if (CoverMesh) CoverMesh->destroy();
    if (TrueGlobalMesh) TrueGlobalMesh->destroy();
  }

private:
  DistributedMesh(const DistributedMesh &); // pas de construction par copie
   void operator=(const DistributedMesh &);// pas affectation par copy
};

// Helpers géométriques (liens géométrie (msh3.cpp) - FESpace (lgmat.cpp))
template<class Mesh>
KN<int> keptElements(const Mesh& Th, const KN<int>& coverPartition, long sizeoverlaps, int rank);

// Sous maillage d'intersections
template<class Mesh>
Mesh* buildIntersectionSubmesh(const DistributedMesh<Mesh>& D, int j, KN<int>& n2o);

template<class Mesh>
Mesh* buildInterfaceSubmesh(const DistributedMesh<Mesh>& D, int j, KN<int>& n2o);

// Partition globale
template<class Mesh>
int computeGlobalPartition(const Mesh& Th, KN<int>& part, const std::string& method, pcommworld comm = nullptr, bool broadcast = true, int nWorkers = 0);

int detectDistributionMode(int localNt, pcommworld comm);

int agreeOnStatus(int local, pcommworld comm);
// Collective error codes
enum DistributeStatus {
  DIST_OK = 0,
  DIST_PART_SIZE = 1, // partition[] does not have Th.nt values
  DIST_PART_RANGE = 2, // partition[k] outside [0, mpisize[
  DIST_PART_EMPTY = 3, // At least one rank without element
  DIST_NT_TOO_SMALL = 4, // mpisize >= Th.nt
  DIST_PART_FAILED = 5, // partitioner return error
  DIST_METHOD_BAD = 6, // unknown partmethod
  DIST_SCOTCH_NA = 7, // build without scotch
  DIST_PARMETIS_NA = 8, // build without parmetis
  DIST_PARMETIS_SCAT = 9, // parmetis incompatible with scatter mode
  DIST_ARG_METHOD = 10, // partmethod different between ranks
  DIST_ARG_OVERLAP = 11, // overlap differs between ranks
  DIST_ARG_KEEPGLOBAL = 12, // keepGlobal differs between ranks
  DIST_ARG_SCATTER = 13, // scatter does not agree with detected mode
  DIST_ARG_PARTWORKERS = 14 // partworkers differs between ranks
};
const char* distributeStatusMessage(int status);
// 0 if method known and available in build
int checkPartitionMethod(const std::string& method);
// 0 if part is usable
int checkPartitionUsable(const KN<int>& part, int nt, bool allowEmpty = false);

template<class Mesh>
void sendMesh(const Mesh& Th, int dest, pcommworld comm);

template<class Mesh>
Mesh* recvMesh(int src, pcommworld comm);

void sendPartition(const KN<int>& part, int dest, pcommworld comm);

KN<int> recvPartition(int n, int src, pcommworld comm);

KN<long> distributedDofNumbering(pcommworld comm, const KN<KN<long>>& dofI, const KN<double>& Ddof, int nLocDof, long& globalNdof);

double checkPartitionOfUnity(pcommworld comm, const KN<KN<long>>& dofI, const KN<double>& Ddof, int nLocDof);

int checkIntersectionSymmetry(pcommworld comm, const KN<KN<long>>& dofI);

// 0 if user partition equal on all ranks, 1 else
int checkPartitionConsistency(pcommworld comm, const KN<int>& part);

template<class R>
R distributedReduce(pcommworld comm, R local);
