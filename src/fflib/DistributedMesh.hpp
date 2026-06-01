template <class Mesh>
class DistributedMesh : public RefCounter {
public: 
  Mesh * LocalMesh;
  DistributedMesh() : LocalMesh(nullptr) {}
  DistributedMesh(Mesh * locmesh) : LocalMesh(locmesh) {}
  ~DistributedMesh() {
    if (LocalMesh) LocalMesh->destroy();
  }

private:
  DistributedMesh(const DistributedMesh &); // pas de construction par copie
   void operator=(const DistributedMesh &);// pas affectation par copy
};