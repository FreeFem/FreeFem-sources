#pragma once
#include "DistributedMesh.hpp"

template<class Mesh>
class DGFESpace : public GFESpace<Mesh> {
public:
    CountPointer<const DistributedMesh<Mesh>> dmesh;
    // dof data (moved from v_dfes in lgfem.hpp)
    KN<KN<long>> dofIntersectionDof;
    KN<double> Ddof;
    double pouResidual = -1.0;
    mutable KN<long> loc2globDof, numberingDof;
    mutable long globalNdof = 0;
    mutable bool numberingBuilt = false, loc2globBuilt = false;
    mutable void* exchangeHandle = nullptr;

    DGFESpace(const DistributedMesh<Mesh>& D, const GTypeOfFE<Mesh>& tfe, int nbequibe = 0, int* equibe = nullptr);
    ~DGFESpace();
    const DistributedMesh<Mesh>& DTh() const { return *dmesh; }
    void ensureLoc2glob() const;
    void ensureNumbering() const;
    // Copy forbidden
    DGFESpace(const DGFESpace&) = delete;
    void operator=(const DGFESpace&) = delete;

private:
    void buildDistributedDofData();
};