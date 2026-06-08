#include <vector>
#include "RNM.hpp"

template <class Mesh>
class DistributedMesh : public RefCounter {
public: 
  Mesh * LocalMesh;
  Mesh* BorderMesh;

  int overlap;
  int interfaceLabel;

  KN<int> neighborRanks;
  KN<double> partitionOfUnity;
  KN<int> localToGLobalVertex;
  std::vector<KN<double>> intersectionSupport;
  KN<int> globalPartition;

  DistributedMesh() : LocalMesh(nullptr), BorderMesh(nullptr), overlap(0), interfaceLabel(10),
  neighborRanks(0), partitionOfUnity(0), localToGLobalVertex(0), globalPartition(0) {}
//  DistributedMesh(Mesh * locmesh) : LocalMesh(locmesh) {}
  ~DistributedMesh() {
    if (LocalMesh) LocalMesh->destroy();
    if (BorderMesh) BorderMesh->destroy();
  }

private:
  DistributedMesh(const DistributedMesh &); // pas de construction par copie
   void operator=(const DistributedMesh &);// pas affectation par copy
};