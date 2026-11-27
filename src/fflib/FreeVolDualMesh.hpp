/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

//////////////////////////////////////////////////////////////////////////////
// FreeVol over a dual mesh: the base class FreeVol_DualModel of the models with P1 unknowns (one value by vertex
// = by control volume of the dual mesh), with the routines common to these models (the slopes, their limiter,
// the reconstructed values on the interfaces, the time step with a domain decomposition)
// the method is given by the type of the dual mesh, so by the dimension of the mesh:
//   1D (meshdualL), 3D (meshdual3): the median dual mesh (class dualmeshN of fvdualmesh.hpp)
//   2D (meshdual): the dual mesh of F. Bouchut (class dualmesh of fvdualmesh.hpp, the centers of the inscribed
//       circles), copied in the format of dualmeshN (FromBouchut), so that the routines and the models are the same
//       in 1D, 2D, 3D, as FreeVol_Geom for the primal mesh; the 2D method keeps its slopes (Grad: the slopes in the
//       triangles with the centers of gravity of the control volumes, then their mean, meanslope) and its time step
//       on the dual edges (DtEdge)
// a model: class FreeVol_X : public FreeVol_DualModel with the 2 constructors FreeVol_DualModel(dm) (dm a dualmesh
// in 2D, a dualmeshN in 1D, 3D); its creator: d == 2 ? *(const dualmesh *)dual : *FreeVol_DualModel::DualN(dual, d)
//////////////////////////////////////////////////////////////////////////////

#ifndef FREEVOLDUALMESH_HPP
#define FREEVOLDUALMESH_HPP

#include "FreeVol.hpp"
#include "fvdualmesh.hpp"   // the dual meshes: class dualmesh (2D) and class dualmeshN (1D, 3D)

class FreeVol_DualModel : public FreeVol_model
{
public:
    const dualmesh *dmF;    // 2D: the dual mesh of F. Bouchut, given to fvscheme (0 in 1D, 3D)
    const dualmeshN *dmN;   // 1D, 3D, surface: the median dual mesh, given to fvscheme (0 in 2D)
    dualmeshN dF;           // 2D: the dual mesh of F. Bouchut in the format of dualmeshN (FromBouchut)
    const dualmeshN *pD;    // the dual mesh in the format of dualmeshN: &dF (2D) or dmN (1D, 3D, surface)
    mutable KN<double> cvn; // surface (meshS): [3i+c] the unit normal of the control volume i (the mean of the
                            // normals of the triangles around the vertex i, weights their areas; periodic: of the
                            // triangles around all the vertices of its class, the normal of the union)
    mutable KN<double> cvdef;// surface: [3i+c] the sum on the interfaces of the control volume i of area*N, N in
                            // the plane of i as in the fluxes (0 on a plane: not 0 on a curved surface; periodic:
                            // without the boundary pieces of the periodic borders, inside the union)

    FreeVol_DualModel(const dualmesh &dm) : FreeVol_model(), dmF(&dm), dmN(0), pD(&dF) { FromBouchut(dm, dF); }
    FreeVol_DualModel(const dualmeshN &dm) : FreeVol_model(), dmF(0), dmN(&dm), pD(&dm) { if (dm.dim == 4) SurfNormals(); }
    const dualmeshN &DM() const { return *pD; }
    bool dual2D() const { return dmF != 0; }
    bool surf() const { return pD->dim == 4; }// surface (meshS)

    // periodic (fvscheme FV(..., periodic=...), see FreeVol_model::SetPeriodic): the control volumes of the vertices
    // of a periodic class (geom.vper, the vertices glued by SetPeriodic on the primal mesh) are one control volume:
    // the boundary pieces on the periodic borders (geom.plab) are inside it (no flux: also right for the non
    // conservative fluxes), U, the slopes and Up are the means on the class weighted by the volumes (Up: the flux
    // balance of the union / its volume); the models which do it say PeriodicOK() (advectiondm, igddm, swtopodm, eulerdm)
    // built at the first use (the model is built again after a change of mesh)
    mutable bool pinit = false;
    mutable KN<long> pcls;      // [i] the class of the vertex i (index in pbeg), -1 if not periodic
    mutable KN<long> pbeg, pmem;// the classes: the vertices pmem[pbeg[c] .. pbeg[c+1]-1]
    mutable KN<char> pskip;     // [m-npe] the boundary piece m is on a periodic border
    mutable KN<double> pgrav;   // [3i+c] the center of gravity of the control volume of i (the reconstruction, the
                                // limiter); periodic: of the union of the control volumes of the class, seen from i
    // the coordinates of the vertex i of the primal mesh
    R3 PVertex(long i) const {
        if (geom.fedim == 2) { const R2 &P = (*(const Mesh *)geom.pTh)(i); return R3(P.x, P.y, 0.); }
        if (geom.fedim == 3) return (*(const Mesh3 *)geom.pTh)(i);
        if (geom.fedim == 1) return (*(const MeshL *)geom.pTh)(i);
        return (*(const MeshS *)geom.pTh)(i);
    }
    bool PeriodicOK() const { return false; }    // a dual model must use the functions below (see advectiondm)
    void PeriodicInit() const {
        if (pinit) return;
        pinit = true;
        const dualmeshN &D = DM();
        const long nv = D.nv;
        pcls.resize(nv); pcls = -1L;
        pbeg.resize(1); pbeg[0] = 0; pmem.resize(0); pskip.resize(0);
        pgrav.resize(3 * nv); pgrav = D.cvgrav;
        if (!geom.vper.N() || geom.nv != nv) return;
        // the classes with 2 vertices or more
        KN<long> cnt(nv); cnt = 0L;
        for (long i = 0; i < nv; ++i) cnt[geom.vper[i]]++;
        long nc = 0, nm = 0;
        KN<long> cof(nv); cof = -1L;
        for (long i = 0; i < nv; ++i) if (cnt[i] > 1) { cof[i] = nc++; nm += cnt[i]; }
        pbeg.resize(nc + 1); pmem.resize(nm);
        pbeg = 0L;
        for (long i = 0; i < nv; ++i) { long c = cof[geom.vper[i]]; if (c >= 0) pbeg[c + 1]++; }
        for (long c = 0; c < nc; ++c) pbeg[c + 1] += pbeg[c];
        KN<long> fill(nc); fill = 0L;
        for (long i = 0; i < nv; ++i) { long c = cof[geom.vper[i]]; if (c >= 0) { pmem[pbeg[c] + fill[c]++] = i; pcls[i] = c; } }
        // the center of gravity of the union of the control volumes of a class, seen from each vertex i of the class
        // (the vertices of a class are the translations of the same point: the control volume k translated by P_i-P_k)
        for (long c = 0; c < nc; ++c) {
            R3 gm; double sv = 0.;
            for (long p = pbeg[c]; p < pbeg[c + 1]; ++p) {
                long k = pmem[p];
                gm += (R3(D.cvgrav[3 * k], D.cvgrav[3 * k + 1], D.cvgrav[3 * k + 2]) - PVertex(k)) * D.cvvol[k];
                sv += D.cvvol[k];
            }
            gm /= sv;
            for (long p = pbeg[c]; p < pbeg[c + 1]; ++p) {
                long i = pmem[p];
                R3 gi = gm + PVertex(i);
                for (int l = 0; l < 3; ++l) pgrav[3 * i + l] = gi[l];
            }
        }
        // the boundary pieces on the periodic borders
        pskip.resize(D.nde - D.npe); pskip = (char)0;
        for (long m = D.npe; m < D.nde; ++m)
            for (long l = 0; l < geom.plab.N(); ++l)
                if (D.dbelab[m - D.npe] == geom.plab[l] && pcls[D.cvLR[2 * m]] >= 0) { pskip[m - D.npe] = 1; break; }// glued
        if (surf() && nc) SurfNormals();// surface: the normals and cvdef of the unions of the control volumes
        if (verbosity > 1 && !mpirank) cout << " fvscheme periodic (dual mesh): " << nc << " classes of control volumes" << endl;
    }
    // the boundary piece m is inside a periodic control volume: no flux
    bool PSkip(long m) const {
        PeriodicInit();
        return pskip.N() && m >= DM().npe && pskip[m - DM().npe];
    }
    // the mean of w (nc components by vertex: w[nc*i+c]) on each periodic class, weighted by the volumes
    void PeriodicMean(KN_<double> w, int nc) const {
        PeriodicInit();
        const dualmeshN &D = DM();
        for (long c = 0; c + 1 < pbeg.N(); ++c)
            for (int k = 0; k < nc; ++k) {
                double sv = 0., sw = 0.;
                for (long p = pbeg[c]; p < pbeg[c + 1]; ++p) { long i = pmem[p]; sv += D.cvvol[i]; sw += D.cvvol[i] * w[nc * i + k]; }
                for (long p = pbeg[c]; p < pbeg[c + 1]; ++p) w[nc * pmem[p] + k] = sw / sv;
            }
    }
    // U at the beginning of the flux computation (the same value on the vertices of a class), Up at the end (the
    // flux balance of the union of the control volumes of a class / its volume)
    void PeriodicU(KN_d *U) const { PeriodicInit(); if (pbeg.N() > 1) for (int j = 0; j < Neq; ++j) PeriodicMean(U[j], 1); }
    void PeriodicUp(KN_d *Up) const { PeriodicUp(Up, Neq); }
    void PeriodicUp(KN_d *Up, int n) const { PeriodicInit(); if (pbeg.N() > 1) for (int j = 0; j < n; ++j) PeriodicMean(Up[j], 1); }
    bool PeriodicVertex(long i) const { PeriodicInit(); return pcls[i] >= 0; }

    // the dual mesh given to fvscheme in 1D (a dualmeshL*), 3D (a dualmesh3*) or on a surface (a dualmeshS*), d the
    // mesh of the unknowns (1 meshL, 3 mesh3, 4 meshS)
    static const dualmeshN *DualN(const void *dual, int d) {
        if (d == 1) return (const dualmeshL *)dual;
        if (d == 4) return (const dualmeshS *)dual;
        return (const dualmesh3 *)dual;
    }

    // surface (meshS): the momentum of a control volume is in its tangent plane (the plane of the normal cvn); on
    // an interface the state of the right control volume is turned into the plane of the left one (the rotation
    // which turns its normal into the one of the left), the flux computed in the plane of the left control volume,
    // and the momentum flux turned back into the plane of the right one (as the cells of FreeVol_swtopo.cpp)
    // (built again by PeriodicInit with the periodic classes: the vertices of a class have the normal of the union,
    // the boundary pieces of the periodic borders are not in cvdef, as in the fluxes)
    void SurfNormals() const {
        const dualmeshN &D = DM();
        const MeshS &Th = *(const MeshS *)D.primal;
        const bool per = pinit && pbeg.N() > 1;// periodic classes (PeriodicInit)
        cvn.resize(3*D.nv); cvn = 0.;
        for (long k=0; k<Th.nt; ++k) {
            R3 P0=Th(Th(k,0)), P1=Th(Th(k,1)), P2=Th(Th(k,2)), Nk=(P1-P0)^(P2-P0);// 2 |K| n_K
            for (int l=0; l<3; ++l) for (int c=0; c<3; ++c) cvn[3*Th(k,l)+c] += Nk[c];
        }
        if (per)// the sum on the class (the triangles around all its vertices: the union, translated)
            for (long cl=0; cl+1<pbeg.N(); ++cl)
                for (int c=0; c<3; ++c) {
                    double sn = 0.;
                    for (long p=pbeg[cl]; p<pbeg[cl+1]; ++p) sn += cvn[3*pmem[p]+c];
                    for (long p=pbeg[cl]; p<pbeg[cl+1]; ++p) cvn[3*pmem[p]+c] = sn;
                }
        for (long i=0; i<D.nv; ++i) {
            double n = sqrt(cvn[3*i]*cvn[3*i]+cvn[3*i+1]*cvn[3*i+1]+cvn[3*i+2]*cvn[3*i+2]);
            for (int c=0; c<3; ++c) cvn[3*i+c] /= n;
        }
        // the sum of area*N on the interfaces of each control volume, with the normals of the fluxes (SurfLR): in
        // the plane of the left control volume, turned into the plane of the right one
        cvdef.resize(3*D.nv); cvdef = 0.;
        double Ul[4] = {0.,0.,0.,0.}, Ur[4] = {0.,0.,0.,0.};
        for (long m=0; m<D.nde; ++m) {
            if (per && pskip.N() && m>=D.npe && pskip[m-D.npe]) continue;// inside the union (PSkip)
            long ivL=D.cvLR[2*m], ivR=D.cvLR[2*m+1];
            R3 N(D.normal[3*m],D.normal[3*m+1],D.normal[3*m+2]);
            SurfLR(ivL, m<D.npe ? ivR : ivL, N, Ul, Ur);
            double n[3] = {N.x,N.y,N.z};
            for (int c=0; c<3; ++c) cvdef[3*ivL+c] += D.area[m]*n[c];
            if (m<D.npe) {
                RotCV(ivL,ivR,n);
                for (int c=0; c<3; ++c) cvdef[3*ivR+c] -= D.area[m]*n[c];
            }
        }
    }
    // surface: the pressure P_i of each control volume times the sum of area*N on its interfaces removed from the
    // momentum balance: a state at rest (constant pressure, the lake at rest) stays at rest, as on a plane (the
    // tangent part of the integral of p N on the border of the control volume is the integral of the surface
    // gradient of p: this sum is an error of the discretization of the curved surface by the control volumes)
    void SurfDefect(KN_d *Up, const KN<double> &P) const {
        const dualmeshN &D = DM();
        for (long i=0; i<D.nv; ++i)
            for (int c=0; c<3; ++c) Up[1+c][i] -= P[i]*cvdef[3*i+c]*(1./D.cvvol[i]);
    }
    // the vector q (3 components) in the plane of the control volume i (its component on the normal removed)
    void ProjCV(long i, double *q) const {
        double qn = q[0]*cvn[3*i]+q[1]*cvn[3*i+1]+q[2]*cvn[3*i+2];
        for (int c=0; c<3; ++c) q[c] -= qn*cvn[3*i+c];
    }
    // the rotation which turns the normal of the control volume i into the one of j, applied to q (Rodrigues)
    void RotCV(long i, long j, double *q) const {
        R3 a(cvn[3*i],cvn[3*i+1],cvn[3*i+2]), b(cvn[3*j],cvn[3*j+1],cvn[3*j+2]), v(q[0],q[1],q[2]);
        R3 k = a^b;// sin(angle) times the axis
        double c = (a,b), s2 = (k,k);
        if (s2 < 1.e-28) return;// the same plane
        R3 w = v*c + (k^v) + k*((k,v)*(1.-c)/s2);
        for (int l=0; l<3; ++l) q[l] = w[l];
    }
    // the momentum U[1..3] of each control volume in its plane
    void SurfU(KN_d *U) const {
        for (long i=0; i<DM().nv; ++i) {
            double q[3] = {U[1][i],U[2][i],U[3][i]};
            ProjCV(i,q);
            for (int c=0; c<3; ++c) U[1+c][i] = q[c];
        }
    }
    // the states of an interface in the plane of its left control volume ivL: Ul (momentum Ul[1..3]) in the plane of
    // ivL, Ur in the plane of ivR then turned into the plane of ivL (a boundary piece: ivR = ivL, Ur = uBC); the
    // normal N in the plane of ivL (unit)
    void SurfLR(long ivL, long ivR, R3 &N, double *Ul, double *Ur) const {
        ProjCV(ivL,Ul+1);
        ProjCV(ivR,Ur+1);
        if (ivR != ivL) RotCV(ivR,ivL,Ur+1);
        double n[3] = {N.x,N.y,N.z};
        ProjCV(ivL,n);
        N = R3(n[0],n[1],n[2]);
        N /= N.norme();
    }
    // the dual mesh of the model: the one given to fvscheme FV(Thstar, [V1, ...], ...), it must be still built (the
    // model built again after a change of mesh, Thstar destroyed), or without the dual mesh (dual = 0): the dual
    // mesh built in the script on the mesh of V1
    // (meshdual Thstar = builddualmesh(Th) in 2D, meshdualL in 1D, meshdual3 in 3D), found in the registry of
    // fvdualmesh (the dual meshes and their primal mesh); it must be the only one of this mesh, and of its size (a
    // check only: a dual mesh keeps its primal mesh, so a new mesh cannot have the address of a mesh with a dual)
    static const void *FindDual(const void *given, const vector< FreeVol_FE > &V, const char *eq, const char *ex) {
        if (given) {
            if (!DualIsBuilt(given)) {
                cout << "fvscheme " << eq << ": the dual mesh given to fvscheme is destroyed (end of the block of meshdual Thstar = builddualmesh(Th)?), it must exist while the fvscheme is used" << endl;
                ExecError("FreeVol: the dual mesh of the fvscheme is destroyed");
            }
            return given;
        }
        int d = V.size() ? V[0].dim : 0;
        const void *primal = V.size() ? V[0].mesh() : 0;
        long nv = 0, nt = 0;
        if (d == 1) { nv = ((const MeshL *)primal)->nv; nt = ((const MeshL *)primal)->nt; }
        else if (d == 4) { nv = ((const MeshS *)primal)->nv; nt = ((const MeshS *)primal)->nt; }
        else if (d == 2) { nv = ((const Mesh *)primal)->nv; nt = ((const Mesh *)primal)->nt; }
        else if (d == 3) { nv = ((const Mesh3 *)primal)->nv; nt = ((const Mesh3 *)primal)->nt; }
        std::vector< const void * > duals;
        if (primal && d >= 1 && d <= 4) DualsOfPrimal(primal, d, duals);
        const void *dual = 0;
        int n = 0;
        for (size_t i = 0; i < duals.size(); ++i) {
            bool ok;
            if (d == 2) {
                const dualmesh *F = (const dualmesh *)duals[i];
                ok = F->cvarea && F->edgenumbering && F->cvarea->N() == nv && F->edgenumbering->N() == 3*nt;
            }
            else ok = ((const dualmeshN *)duals[i])->nv == nv && ((const dualmeshN *)duals[i])->nt == nt;
            if (ok) { dual = duals[i]; n++; }
        }
        if (n == 0) {
            cout << "fvscheme " << eq << ": no dual mesh of the mesh of " << ex << ", build it before fvscheme, ex. in 2D: meshdual Thstar = builddualmesh(Th); fvscheme FV(" << ex << ", eq=\"" << eq << "\", ...) (meshdualL in 1D, meshdual3 in 3D)" << endl;
            ExecError("FreeVol: no dual mesh of the mesh of the unknowns (builddualmesh)");
        }
        if (n > 1) {
            cout << "fvscheme " << eq << ": " << n << " dual meshes of the mesh of " << ex << ", give the one of the scheme: fvscheme FV(Thstar, " << ex << ", eq=\"" << eq << "\", ...)" << endl;
            ExecError("FreeVol: several dual meshes of the mesh of the unknowns, give it to fvscheme");
        }
        return dual;
    }
    // the dual mesh of the script (fvscheme) must be the dual mesh of the mesh of the model: in init and at each
    // fluxupdate (the dual mesh of the script built again on another mesh, or destroyed: end of its block)
    void CheckDual(const char *what) const {
        if (!DualIsBuilt(dmF ? (const void *)dmF : (const void *)dmN)) {
            cout << what << ": the dual mesh of the fvscheme is destroyed (end of the block of meshdual Thstar = builddualmesh(Th)?), it must exist while the fvscheme is used" << endl;
            ExecError("FreeVol: the dual mesh of the fvscheme is destroyed");
        }
        const void *primal = dmF ? (const void *)dmF->primal : dmN->primal;
        if (primal != geom.pTh || DM().dim != geom.fedim) {
            cout << what << ": the dual mesh is not the dual mesh of the mesh of the model (meshdual / meshdualL / meshdual3 / meshdualS Thstar = builddualmesh(Th), after a change of mesh too; resetfvscheme)" << endl;
            ExecError("FreeVol: the dual mesh is not the dual mesh of the mesh of the model");
        }
    }

    // 2D: the dual mesh of F. Bouchut in the format of dualmeshN: the interfaces m < npe (dual edges) and the
    // boundary pieces npe <= m < nde (2 by boundary edge: the uBC at the piece of the vertex ivL, the other vertex
    // ivR of the boundary edge, with the weights 0.75, 0.25), the vectors with 3 components (z = 0)
    static void FromBouchut(const dualmesh &F, dualmeshN &D) {
        const Mesh &Th = *F.primal;
        const KN<long> &nbedtab = *(F.nbed);
        const long npe = nbedtab[0], nde = nbedtab[1], nv = Th.nv;
        const KN<long> &cvLR = *(F.cvLR);
        const KN<double> &normals = *((F.extras)->dmnormal), &lengths = *((F.extras)->dmlengthedge);
        const KN<double> &midedge = *((F.extras)->dmmidedge), &cvgrav = *(F.cvgrav), &cvarea = *(F.cvarea);
        const KN<double> &Nbdryvert = *((F.extras)->dmNbdryvert);
        const KN<long> &dbelab = *((F.normalized)->dbelab);
        D.destroy(); D.dim = 2; D.primal = &Th; D.nv = nv; D.nt = Th.nt; D.npe = npe; D.nde = nde;
        D.cvLR.resize(2*nde); D.normal.resize(3*nde); D.area.resize(nde); D.mid.resize(3*nde);
        D.dbelab.resize(nde-npe); D.bvert.resize(3*(nde-npe)); D.bweight.resize(3*(nde-npe));
        D.cvvol.resize(nv); D.cvgrav.resize(3*nv); D.Nbdryvert.resize(3*nv);
        for (long m = 0; m < nde; ++m) {
            D.cvLR[2*m] = cvLR[2*m]; D.cvLR[2*m+1] = cvLR[2*m+1];
            D.normal[3*m] = normals[2*m]; D.normal[3*m+1] = normals[2*m+1]; D.normal[3*m+2] = 0.;
            D.mid[3*m] = midedge[2*m]; D.mid[3*m+1] = midedge[2*m+1]; D.mid[3*m+2] = 0.;
            D.area[m] = lengths[m];
        }
        for (long m = npe; m < nde; ++m) {
            long b = 3*(m-npe), ivL = cvLR[2*m];
            long ivR = ((m-npe)%2 == 0) ? cvLR[2*(m+1)] : cvLR[2*(m-1)];// the other vertex of the boundary edge
            D.dbelab[m-npe] = dbelab[m-npe];
            D.bvert[b] = ivL; D.bvert[b+1] = ivR; D.bvert[b+2] = ivL;
            D.bweight[b] = 0.75; D.bweight[b+1] = 0.25; D.bweight[b+2] = 0.;
        }
        for (long i = 0; i < nv; ++i) {
            D.cvvol[i] = cvarea[i];
            D.cvgrav[3*i] = cvgrav[2*i]; D.cvgrav[3*i+1] = cvgrav[2*i+1]; D.cvgrav[3*i+2] = 0.;
        }
        // the mean unit outward normal at the vertices of the physical boundary: the one of F. Bouchut, given at the
        // vertex of the first piece of each boundary edge
        D.Nbdryvert = 0.;
        for (long ib = 0; ib < (nde-npe)/2; ++ib) {
            long i = cvLR[2*(npe+2*ib)];
            D.Nbdryvert[3*i] = Nbdryvert[2*ib]; D.Nbdryvert[3*i+1] = Nbdryvert[2*ib+1];
        }
        D.buildlists();
    }

    // the centered slopes g[3*i+c] of the P1 function u in the control volumes: 1D, 3D: the mean of the gradients
    // in the elements around the vertex (matgrad); 2D: the method of F. Bouchut (the slopes in the triangles with
    // the values at the centers of gravity of the control volumes, mat_dxygrav, then their mean, meanslope)
    void Grad(const KN<double> &u, KN<double> &g) const {
        const long nv = DM().nv;
        g.resize(3*nv);
        if (dmF) {
            KN<double> dxyt(2*dmF->primal->nt), g2(2*nv);
            PeriodicInit();
            if (pbeg.N() > 1) {// periodic: the slopes in the triangles as mat_dxygrav, with the centers of gravity of
                               // the unions of the control volumes of the classes (pgrav)
                const Mesh &Th = *dmF->primal;
                for (long k = 0; k < Th.nt; ++k) {
                    long v0 = Th(k,0), v1 = Th(k,1), v2 = Th(k,2);
                    double vec1x = pgrav[3*v1]-pgrav[3*v0], vec1y = pgrav[3*v1+1]-pgrav[3*v0+1];
                    double vec2x = pgrav[3*v2]-pgrav[3*v0], vec2y = pgrav[3*v2+1]-pgrav[3*v0+1];
                    double vec1sq = vec1x*vec1x+vec1y*vec1y, vec2sq = vec2x*vec2x+vec2y*vec2y, vec1vec2 = vec1x*vec2x+vec1y*vec2y;
                    double den = vec1sq*vec2sq-vec1vec2*vec1vec2;
                    double com1 = vec1sq/den, com2 = vec2sq/den, com12 = vec1vec2/den;
                    double a1 = com2*u[v1]-com12*u[v2]-(com2-com12)*u[v0], a2 = -com12*u[v1]+com1*u[v2]+(com12-com1)*u[v0];
                    dxyt[2*k] = a1*vec1x+a2*vec2x;
                    dxyt[2*k+1] = a1*vec1y+a2*vec2y;
                }
            }
            else dxyt = (*((dmF->matrices)->mat_dxygrav))*u;// slopes in triangles
            meanslope(dxyt,g2);// centered slopes in control volumes
            for (long i = 0; i < nv; ++i) { g[3*i] = g2[2*i]; g[3*i+1] = g2[2*i+1]; g[3*i+2] = 0.; }
        }
        else g = (*(dmN->matgrad))*u;
        PeriodicInit();
        if (pbeg.N() > 1) PeriodicMean(g, 3);// periodic: the slope of the union of the control volumes
    }
    // the increment of the slope g in the control volume i from its center of gravity to the centroid of the
    // interface m (g . GE)
    double IncGE(const KN<double> &g, long i, long m) const {
        const dualmeshN &D = DM();
        PeriodicInit();
        double inc = 0.;
        for (int c = 0; c < 3; ++c) inc += g[3*i+c]*(D.mid[3*m+c]-pgrav[3*i+c]);// periodic: the center of the union
        return inc;
    }
    // the values of the control volume i reconstructed at the centroid of the interface m (gU[j]: the slopes of
    // the unknown j, used if ord > 1)
    void GetU(long m,long i,const KN_d *U,const KN<double> *gU,double *u) const {
        for(int j=0;j<Neq;++j){
            double uij = U[j][i];
            if(ord>1) uij += IncGE(gU[j],i,m);
            u[j] = uij;
        }
    }
    // the momentum reflected by a wall of unit normal N (q = h*u, d components): q - 2 (q.N) N; 2D as the scheme of
    // F. Bouchut (in the frame of the normal and the tangent)
    void Reflect(int d,const R3 &Nn,double h,const double *q,double *qR) const {
        if (d == 2) {
            double u = q[0]/max(h,1.e-13), v = q[1]/max(h,1.e-13);
            double un = u * Nn.x + v * Nn.y;// normal velocity
            double ut = -u * Nn.y + v * Nn.x;// tangential velocity (unchanged)
            qR[0] = h * (-un * Nn.x + ut * (-Nn.y));
            qR[1] = h * (-un * Nn.y + ut * ( Nn.x));
            return;
        }
        double un = 0.;// normal velocity
        for(int c=0; c<d; ++c) un += q[c]/max(h,1.e-13)*Nn[c];
        for(int c=0; c<d; ++c) qR[c] = h*(q[c]/max(h,1.e-13) - 2.*un*Nn[c]);
    }

    // the limitation of the slopes Gw[3*i+c] in the control volumes
    // ilimiter = 0: no limiter, 1: maximum principle type, 2: minmod, 3: Sweby, 4: Van Albada type, 12=1+2, 13=1+3, 14=1+4
    double SlopeLimiterVF(const KN<double> &w, const KN<double> &Gw, KN<double> &alpha, const long &ilimiter) const
    {
        const dualmeshN &D = DM();
        // the min and the max of w on the neighbours (periodic: on the neighbours of all the control volumes of the class)
        KN<double> wmin(D.nv), wmax(D.nv);
        for(long i=0;i<D.nv;++i){
            double wi=w[i], wmini = wi, wmaxi = wi;
            for(long p=D.ifbeg[i];p<D.ifbeg[i+1];++p){//loop on interfaces around control volume i
                long m=D.iflist[p];
                if (m<D.npe) {
                    long ie=(D.iforient[p]==1)? D.cvLR[2*m+1] : D.cvLR[2*m];//opposite control volume
                    wmini = min(wmini,w[ie]);
                    wmaxi = max(wmaxi,w[ie]);
                }
            }
            wmin[i] = wmini; wmax[i] = wmaxi;
        }
        PeriodicInit();
        for (long c = 0; c + 1 < pbeg.N(); ++c) {
            double a = 1e300, b = -1e300;
            for (long p = pbeg[c]; p < pbeg[c + 1]; ++p) { a = min(a, wmin[pmem[p]]); b = max(b, wmax[pmem[p]]); }
            for (long p = pbeg[c]; p < pbeg[c + 1]; ++p) { wmin[pmem[p]] = a; wmax[pmem[p]] = b; }
        }
        for(long i=0;i<D.nv;++i){//loop on control volumes
            R3 G(pgrav[3*i],pgrav[3*i+1],pgrav[3*i+2]);//gravity center of control volume i (periodic: of the union)
            R3 Gwi(Gw[3*i],Gw[3*i+1],Gw[3*i+2]);//gradient in control volume i
            double wi=w[i], wmini = wmin[i], wmaxi = wmax[i];
            double a = 1.;
            if (ilimiter!=0){
                for(long p=D.ifbeg[i];p<D.ifbeg[i+1];++p){
                    long m=D.iflist[p];
                    if (PSkip(m)) continue;// periodic: a boundary piece inside the union of the control volumes
                    R3 E(D.mid[3*m],D.mid[3*m+1],D.mid[3*m+2]);//centroid of the interface
                    R3 GE=E-G;
                    double incwe = (Gwi,GE);//increment w(E)-w(G)
                    if ((ilimiter==1)||(ilimiter>10)){
                        double ae=1., we = wi + incwe;
                        if (we > wi) ae = .9*(wmaxi-wi)/(we-wi);
                        else if(we < wi ) ae = .9*(wi-wmini)/(wi-we);
                        a = min(a,ae);
                    }
                    if ((ilimiter>1)&&(m<D.npe)){// interface not on the boundary
                        long ie=(D.iforient[p]==1)? D.cvLR[2*m+1] : D.cvLR[2*m];//opposite control volume
                        R3 Gwie(Gw[3*ie],Gw[3*ie+1],Gw[3*ie+2]);//gradient in the opposite control volume
                        double incwee = (Gwie,GE);//increment computed with the neighboring gradient
                        double ae=limiterfunc(incwe,incwee,ilimiter%10)*incwe/(incwe*incwe+1.e-20);
                        a = min(a,ae);
                    }
                }
            }
            alpha[i]=a;
        }//end of loop on control volumes
        for (long c = 0; c + 1 < pbeg.N(); ++c) {// periodic: the same limitation in the control volumes of a class
            double a = 1.;
            for (long p = pbeg[c]; p < pbeg[c + 1]; ++p) a = min(a, alpha[pmem[p]]);
            for (long p = pbeg[c]; p < pbeg[c + 1]; ++p) alpha[pmem[p]] = a;
        }
        return 0.;
    }//end of SlopeLimiterVF

    // the slopes at the vertices of the physical boundary: 0 (ibdryslope = 0) or tangent to the boundary (1)
    void BdrySlope(KN<double> &Gw, const long &ibdryslope) const {
        const dualmeshN &D = DM();
        for(long i=0;i<D.nv;++i){
            R3 Nb(D.Nbdryvert[3*i],D.Nbdryvert[3*i+1],D.Nbdryvert[3*i+2]);//mean unit normal, or 0 if not on physical boundary
            if ((Nb,Nb)<=1.e-30) continue;
            if (PeriodicVertex(i)) continue;// periodic: not a physical boundary
            if (ibdryslope==0) for(int c=0;c<3;++c) Gw[3*i+c]=0.;
            else if (ibdryslope==1) {
                double GwN=Gw[3*i]*Nb.x+Gw[3*i+1]*Nb.y+Gw[3*i+2]*Nb.z;
                for(int c=0;c<3;++c) Gw[3*i+c]-=GwN*Nb[c];// imposing Gw.N=0 on boundary
            }
        }
    }

    // the control volumes not in the time step with a domain decomposition (DdmCutCV2: 2D, DdmCutCVN: 1D, 3D)
    void DdmCutCV(KN<long> &cut) const {
        if (dmF) DdmCutCV2(cut); else DdmCutCVN(cut);
    }
    // the control volumes cut by the border of the local mesh of a domain decomposition (a boundary dual edge
    // with the label -111111): their area and their interfaces are not complete, they are not in the time step
    // (overlap vertices, set again by the exchange; each vertex is complete on one process at least)
    void DdmCutCV2(KN<long> &cut) const {
        const dualmeshN &D = DM();
        const KN<long> &cvLR = D.cvLR;
        const long npe = D.npe, nde = D.nde;
        cut.resize(D.nv); cut=0L;
        for(long m=npe;m<nde;++m) if(D.dbelab[m-npe]==-111111) cut[cvLR[2*m]]=1;
        if(ord>1) {// at the order 2, their neighbours (2): the slopes in the triangles use the centers of gravity of
                   // the control volumes, wrong for a cut one, so the slope of a neighbour (meanslope) is wrong
            KN<long> c1(cut);
            for(long m=0;m<npe;++m) { long l=cvLR[2*m], r=cvLR[2*m+1]; if(c1[l]==1 && !cut[r]) cut[r]=2; if(c1[r]==1 && !cut[l]) cut[l]=2; }
        }
    }
    // the time step of the dual edge m (the area of its control volumes / its length / speed), not for the control
    // volumes cut by the border of the local mesh (DdmCutCV): dt unchanged if they are all cut; at the order 2, dt
    // unchanged if one of them is cut or a neighbour of a cut one (its slope is wrong, so the speed of the edge)
    double DtEdge(double dt,long m,long npe,long ivL,long ivR,const KN<long> &cut,const KN<double> &cvarea,double le,double speed) const {
        bool cL = cut[ivL], cR = (m<npe) && cut[ivR];
        if (cL && (cR || m>=npe)) return dt;
        if (ord>1 && (cL || cR)) return dt;
        double am = (m>=npe) ? cvarea[ivL] : (cL ? cvarea[ivR] : (cR ? cvarea[ivL] : min(cvarea[ivL],cvarea[ivR])));
        return min(dt,am/le/(speed+1.e-20));
    }
    // the control volumes cut by the border of the local mesh of a domain decomposition (a boundary piece with
    // the label -111111): their volume and their interfaces are not complete, they are not in the time step
    // (overlap vertices, set again by the exchange; each vertex is complete on one process at least), at the
    // order 2 also the control volumes at a distance 1 or 2 of them
    void DdmCutCVN(KN<long> &cut) const {
        const dualmeshN &D = DM();
        cut.resize(D.nv); cut=0L;
        for(long m=D.npe;m<D.nde;++m) if(D.dbelab[m-D.npe]==-111111) cut[D.cvLR[2*m]]=1;
        if(ord>1) {// at the order 2, also their neighbours and the neighbours of these ones: the slope of a cut control
                   // volume is not complete, so the speed of its interfaces, and with the limiters 2, 3, 4 (they use the
                   // slopes of the neighbours) the slopes of its neighbours (overlap of 3 layers of triangles)
            for(int layer=0; layer<2; ++layer) {
                KN<long> c1(cut);
                for(long m=0;m<D.npe;++m) { long l=D.cvLR[2*m], r=D.cvLR[2*m+1]; if(c1[l]) cut[r]=1; if(c1[r]) cut[l]=1; }
            }
        }
    }

// ============================ 2D: the dual mesh of F. Bouchut (dmF) ============================
    // computes mean slopes over control volumes from slopes in triangles around
    double meanslope(const KN<double> &aa,KN<double> &cc) const {
        const dualmesh &Thstar = *dmF;
        const Mesh *pTh = (const Mesh *)dmF->primal;
     //aa: input=gradient in triangles
     //cc: result=gradient in dualmesh control volumes
     ffassert((aa.N()==2*(pTh->nt))&&(cc.N()==2*(pTh->nv)));
     const KN<long> &nbedtab=*(Thstar.nbed);// (size 3)
     const long maxinterf=nbedtab[2];
     const MatriceMorse<double> &mattri2cvprod=*((Thstar.matrices)->mat_tri2cvprod);
     long Thnt=pTh->nt;
     long Thnv=pTh->nv;
     KN<double> dxft(Thnt),dyft(Thnt);
     for (long k=0;k<Thnt;k++){
      dxft(k)=aa(2*k);
      dyft(k)=aa(2*k+1);
     }
     KN<double> vecd(Thnv*maxinterf);
     KN<double> dxres(Thnv);
     KN<double> dyres(Thnv);
     KN<long> nneit(Thnv);
     KN<double> weight(Thnv*maxinterf);
     nneit=*((Thstar.extras)->dmnbneit);
     weight=*((Thstar.extras)->dmvolcvprod);
     const long normalizee=1;
     vecd=mattri2cvprod*dxft;sumwtabprod(vecd,nneit,weight,normalizee,dxres);
     vecd=mattri2cvprod*dyft;sumwtabprod(vecd,nneit,weight,normalizee,dyres);
     for (long i=0;i<Thnv;i++){
      cc(2*i)=dxres(i);
      cc(2*i+1)=dyres(i);
     }
     return 0.;
    }//end of meanslope
};

#endif
