#include <vector>
#include "RNM.hpp"

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