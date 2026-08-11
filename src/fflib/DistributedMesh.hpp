#include <vector>
#include "RNM.hpp"
#include <string>
#ifndef FF_PCOMMWORLD_TYPEDEF
#define FF_PCOMMWORLD_TYPEDEF
typedef void* pcommworld;
#endif

template <class Mesh>
class DistributedMesh : public RefCounter {
public: 
  Mesh * LocalMesh;
  Mesh* BorderMesh;
  const Mesh* GlobalMesh = nullptr;

  int overlap;
  int interfaceLabel;

  KN<int> neighborRanks;
  KN<double> partitionOfUnity;
  KN<int> localToGlobalElement; // n2o (LocalMesh -> Th)
  KN<int> globalPartition;

  pcommworld comm = nullptr;
  
  DistributedMesh() : LocalMesh(nullptr), BorderMesh(nullptr), overlap(0), interfaceLabel(10),
  neighborRanks(), partitionOfUnity(), globalPartition() {}
//  DistributedMesh(Mesh * locmesh) : LocalMesh(locmesh) {}
  ~DistributedMesh() {
    if (LocalMesh) LocalMesh->destroy();
    if (BorderMesh) BorderMesh->destroy();
  }

private:
  DistributedMesh(const DistributedMesh &); // pas de construction par copie
   void operator=(const DistributedMesh &);// pas affectation par copy
};

// Helpers géométriques (liens géométrie (msh3.cpp) - FESpace (lgmat.cpp))
template<class Mesh>
KN<int> keptGlobalElements(const Mesh& Th, const KN<int>& globalPartition, long sizeoverlaps, int rank);

// Sous maillage d'intersections
template<class Mesh>
Mesh* buildIntersectionSubmesh(const DistributedMesh<Mesh>& D, int j, KN<int>& n2o);

// Partition globale
template<class Mesh>
void computeGlobalPartition(const Mesh& Th, KN<int>& part, const std::string& method, pcommworld comm = nullptr);

template<class Mesh>
void sendMesh(const Mesh& Th, int dest, pcommworld comm);

template<class Mesh>
Mesh* recvMesh(int src, pcommworld comm);

KN<long> distributedDofNumbering(pcommworld comm, const KN<KN<long>>& dofI, const KN<double>& Ddof, int nLocDof, long& globalNdof);