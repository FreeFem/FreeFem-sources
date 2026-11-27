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

#include "ff++.hpp"
#include "AFunction_ext.hpp"
#include "AFunction.hpp"
#include "lgmesh3.hpp"    // the FE functions on mesh3 (pf3r) and meshL (pfLr)
#include <array>
#include <map>
#include <string>
#include <vector>

struct KN_d : public KN_<double>  {
    typedef double R;
    KN_d() : KN_(0,0,0){} // make array
    KN_d(R *u,long nn,long s): KN_<double>(u,n,s) {}
    KN_d & operator=(KN_<double> s){ set(s);return *this;}//  copy the pointer, ...
    KN_d & operator=(double d) {KN_<double>::operator=(d); return *this;}
};

// a FE function of FreeFEM (a component), in 1D (meshL), 2D (mesh) or 3D (mesh3):
// the velocity of a model or the unknowns of fluxupdate
struct FreeVol_FE {
    int dim;    // 1 meshL, 2 mesh, 3 mesh3, 4 meshS (surface in 3D) (0: unset)
    pfer f2;
    pf3r f3;
    pfLr fL;
    pfSr fS;
    FreeVol_FE() : dim(0), f2(0,-1), f3(0,-1), fL(0,-1), fS(0,-1) {}
    KN<double> *x() const;           // the values of the FE function (all the components), 0 if unset
    int comp() const;                // the component of this function
    int N() const;                   // number of components of the fespace
    long NbOfDF() const;             // number of dof of the fespace (all the components)
    long NbOfNodes() const;          // number of nodes
    int MaxNbDFPerElement() const;
    long dof(long k, int i) const;   // the node i of the element k
    const void *mesh() const;        // the mesh of the fespace (Mesh, Mesh3, MeshL or MeshS)
    const RefCounter *meshrc() const;// the same, with its reference counter
    KN_<double> values() const;      // the values of this component
};

// the geometry of a mesh for the finite volume schemes, the same in 1D, 2D and 3D:
// cells (segments, triangles, tetrahedra) and their faces (dim+1 faces by cell, face e opposite
// to the vertex e), kf = k*nf+e
// surface (meshS: triangles in 3D): dim = 3 (the vectors, the gradients have 3 components), tdim = 2, the faces
// are the edges, N = the unit outward conormal (in the plane of the triangle), Nc = the unit normal of the triangle
struct FreeVol_Geom {
    int dim, nf, nvf;      // dimension, faces by cell (dim+1), vertices by face (dim)
    int tdim;              // topological dimension: dim, 2 for a surface
    int fedim;             // the mesh: 1 meshL, 2 mesh, 3 mesh3, 4 meshS (as FreeVol_FE::dim)
    bool surf;             // surface (meshS)
    KN<R3> Nc;             // [k]  surface: unit normal of the triangle (empty otherwise)
    long nt, nv;           // number of cells, of vertices
    const void *pTh;       // the mesh (Mesh, Mesh3, MeshL or MeshS)
    KN<double> meas;       // [k]  measure of the cell
    KN<R3> G;              // [k]  centroid of the cell
    KN<long> adj;          // [kf] the neighbour cell, -1 on the boundary
    KN<int> adjf;          // [kf] the face in the neighbour cell
    KN<double> fmeas;      // [kf] measure of the face (1 in 1D)
    KN<R3> N;              // [kf] unit outward normal
    KN<R3> F;              // [kf] centroid of the face
    KN<long> flab;         // [kf] label of a boundary face (0 elsewhere)
    KN<long> fvert;        // [kf*nvf+l] the vertices of the face
    long nbdry;            // number of faces on the true boundary
    long nfaces;           // number of faces (interior + boundary)
    // periodic (fvscheme FV(..., periodic=[[l1,f1],[l2,f2], ...])): the faces of the labels l1 and l2 glued as
    // interior faces (adj, adjf: the cell and the face on the other side, flab = 0)
    long nper;             // number of glued faces (0 without periodic)
    KN<R3> pshift;         // [kf] glued face: the translation from the face of the other side to kf (empty otherwise)
    KN<long> vper;         // [i] the representative of the periodic class of the vertex i (empty without periodic)
    KN<long> plab;         // the labels of the periodic borders (empty without periodic)
    void build(const Mesh &Th);
    void build(const Mesh3 &Th);
    void build(const MeshL &Th);
    void build(const MeshS &Th);
    void build(const FreeVol_FE &fe);   // build on the mesh of a FE function
};

// forward declaration (registry of the models)
class FreeVol_model;
// the dual meshes of the plugin fvdualmesh (fvdualmesh.hpp), given to fvscheme FV(Thstar, [V1,V2], ...)
class dualmesh; class dualmeshN; class dualmeshL; class dualmesh3; class dualmeshS;

class FreeVol_model  {
public:
    static const long Neqmax = 100;
    long Neq;// number of equations in the system
    FreeVol_Geom geom;// the geometry of the mesh (1D, 2D or 3D)
    long numflux,ord,limiter;
    KN<double> parameters;
    KN<long> intparameters;
    KN<long> labelBC;
    MatriceMorse<R> *matgrad; // the gradients: [dim*k+c] row of the component c of the gradient in k
    string *equation;// name of the equation

    // the model: physical flux and eigenvalues on the face kf = k*nf+e, unit normal N
    virtual void Flux(long kf,const R3 &N,const double *U,R *F) const =0;
    virtual void lambda(long kf,const R3 &N,const double *U,double *l) const =0;
    virtual double SlopeLimiterVF(const KN<double> &w, const KN<double> &Gw, KN<double> &alpha, const long &ilimiter) const ;
    virtual void BdrySlope(KN<double> &Gw, const long &ibdryslope) const ;
    double absvp(long kf,const R3 &N,const double *U) const ;
    double lambdamin(long kf,const R3 &N,const double *U) const ;
    double lambdamax(long kf,const R3 &N,const double *U) const ;
    void GetU(long kf,const KN_d *U,const KN<double> *gU,double *u) const ;
    void GetUBC(long kf,const KN_d *U,double *u) const ;
    // the limited slopes, the same for all the models and all the dimensions
    // (virtual: a model can compute its own slopes, ex. igd on rho and u)
    virtual void LimitedSlopes(KN_d *U, KN<double> *gU) const;   // gU[j][dim*k+c]
    // the fluxes of the cells and the time step, defined by each model
    virtual double computedtFlux(KN_d *U,KN_d *Up, KN_d *uBC) const =0;

    FreeVol_model() : Neq(0),numflux(0),ord(0),limiter(0),labelBC(0),matgrad(0),equation(0),parameters(0),intparameters(0) {}
    // the mesh is the mesh of the FE function fe (the velocity of the model)
    void set(long n,long ordd,long numfluxa,const FreeVol_FE &fe, KN<long> labelBCC, const char *equationn, long limiterr, KN<double> parameterss, KN<long> intparameterss) ;
    // periodic= accepted by the model (the primal models: the glued faces are interior faces; a dual model: it
    // must skip the boundary pieces on the periodic borders and take the mean on the periodic classes)
    virtual bool PeriodicOK() const { return true; }
    // periodic boundary conditions: glue the faces of the pairs of labels (after set, by fvscheme)
    // periodic: the expressions of periodic= (as for a fespace), [label1, f1, label2, f2] by pair in 2D
    void SetPeriodic(int nbcperiodic, Expression *periodic, Stack s);
    virtual ~FreeVol_model() ;
};

// the parameters given to fvscheme (or resetfvscheme): to build the model again on a new mesh
struct FreeVol_Params {
    string eq;
    vector< FreeVol_FE > V;    // the FE functions given to fvscheme (the velocity of the advection)
    long order, numflux, limiter;
    KN<long> labelBC;
    KN<double> parameters;
    KN<long> intparameters;
    bool rebuild;    // a new mesh (adaptmesh, ...): build the model again (true) or stop (false)
    const void *dual;          // the dual mesh given to fvscheme (meshdualL, meshdual, meshdual3), 0 if not given
    int nbcperiodic;           // periodic=: the number of pairs of labels (0 without periodic)
    Expression *periodic;      // their expressions (of the instruction fvscheme, not freed here)
};

class FreeVol { // encapsulation of FreeVol_model ptr
public:
    // made by FreeFEM with init() (InitP), not with the constructor: the members are set in init()
    FreeVol_model *FVM;
    FreeVol_Params *par;          // the parameters of the definition
    const RefCounter *pThref;     // the mesh of the model, with a reference: not freed while used, so a
                                  // new mesh cannot have the same address (sure detection of a change)
    FreeVol() : FVM(0), par(0), pThref(0) {}
    ~FreeVol() {destroy();}
    void destroy(){
     if (FVM) {delete FVM; FVM=0;}
     if (par) {delete par; par=0;}
     if (pThref) {pThref->destroy(); pThref=0;}
    }
    void init() { FVM=0; par=0; pThref=0; if(verbosity>99)  cout << " FreeVol init " << this << endl; }
    operator FreeVol_model * () {return FVM;}
};

// Registry of the models: each model (separate plugin FreeVol_<eq>, ex. FreeVol_advection) registers in its Load_Init
// a function which creates it, under the name given by eq= in the script:
//   fvscheme FV([V1,V2], eq="advection", ...)
// V: the FE functions given to fvscheme (meshL, mesh or mesh3: the dimension of the problem)
// the creator returns a new model; numflux = -99 if numflux= is not given (default of the model)
// dual: the dual mesh given to fvscheme, fvscheme FV(Thstar, [V1,V2], eq=..., ...), 0 if not given: a dualmeshL*
// (meshdualL) in 1D, a dualmesh* (meshdual) in 2D, a dualmesh3* (meshdual3) in 3D (plugin fvdualmesh)
typedef FreeVol_model *(*FreeVol_Creator)(const vector< FreeVol_FE > &V, long order,
                                          long numflux, KN<long> labelBC, long limiter,
                                          KN<double> parameters, KN<long> intparameters,
                                          const void *dual);
void FreeVol_RegisterModel(const string &eq, FreeVol_Creator creator);
FreeVol *FreeVol_CreateModel(FreeVol *a, const string &eq, const vector< FreeVol_FE > &V,
                             long order, long numflux, KN<long> labelBC, long limiter,
                             KN<double> parameters, KN<long> intparameters, bool rebuild = true,
                             const void *dual = 0, int nbcperiodic = 0, Expression *periodic = 0, Stack s = 0);
// the mesh of the unknowns is not the mesh of the model (adaptmesh, ...): build the model again
// with the same parameters, or stop if rebuild=false
void FreeVol_MeshChanged(FreeVol *a, const void *pThnew, Stack s);

double sumwtabprod(const KN<double> &aa,const KN<long> &bb,const KN<double> &ww,const long &normalize,KN<double> &cc);
double limiterfunc(double amin,double amax,const long &ilimitertype);
double numfluxgd(const double &gamma,const double* statele,const double* stateri,const double* normalp,double &absvp,double* Fluxesvect);
