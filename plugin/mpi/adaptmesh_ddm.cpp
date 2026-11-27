// adaptmesh_ddm: parallel 2d adaptmesh (MPI) by domain decomposition
//
// Base: Adaptation::operator() of src/fflib/lgmesh.cpp (adaptmesh, bamg).
//
// Principle:
//   pass 1. the mesh (the same on all the processes) is cut into mpisize
//      subdomains by recursive bisection of the barycentres. Each process
//      computes the metric as adaptmesh on its subdomain extended by
//      20 + nbjacobi layers (same metric as in sequential at its own vertices),
//      normalized with the global min/max of the solutions, then adapts its
//      subdomain with bamg, the interface edges being required;
//      the pieces are gathered (MPI_Allgatherv) and glued by the interface
//      vertices; the metric at the vertices of the input mesh is kept.
//   pass 2. a band of nblayer layers around the interfaces of pass 1
//      is adapted again (bisection rotated by pi/4), with the metric of pass 1
//      interpolated in the input mesh (as bamg with its background mesh).
//   pass 3. the same around the cuts of pass 2, each connected component
//      being given entirely to one process: the whole mesh has been adapted at
//      least once without a frozen edge inside (no star at the crossings).
//   With one process, the result is exactly that of adaptmesh.
//
// Usage (FreeFem++-mpi):
//   load "adaptmesh_ddm"
//   Th = adaptmeshddm(Th, u, [vx,vy], err=0.01, nbpass=3, nblayer=6, comm=mpiCommWorld);
//   All the parameters of adaptmesh (metric=, periodic=, nomeshgeneration=,
//   splitin2, nbvx, ...); with one process, result identical to adaptmesh.
//   periodic= in parallel: periodic boundaries frozen during the passes, then
//   a last pass adapts again a band around them on a single process.
//   splitin2 in parallel: splitting of the global mesh at the end.
//   nbvx: limit of the total number of vertices, shared between the processes
//   according to the number of vertices predicted by the metric.
//
// Distributed mode (as parmmg3d with nodeCommunicators), with owned=:
//   Th = LOCAL mesh with overlap (MatCreate), u defined on Th,
//   owned = P0 array (PartitionCreate); returns the adapted local mesh WITHOUT
//   overlap (interfaces: label -111111), rebalanced between the processes:
//     fespace Ph(Th, P0); Ph part; PartitionCreate(Th, part[], P0);
//     Th = adaptmeshddm(Th, u, err=0.01, owned=part[], comm=mpiCommWorld);
//     DmeshReconstruct(Th); MatCreate(Th, A, P1);
//   periodic= in this mode: periodic boundaries frozen during the passes, then a band around them gathered on
//   process 0 and adapted again with the periodicity; this band stays on process 0 (the 2 sides of a periodic
//   boundary in the same local mesh, as macro ThPeriodicity of macro_ddm.idp), the rest is cut between the
//   other processes (freevol_advection2d_ddm_adaptddm_periodic.edp)
//
// Compilation: see build.sh

#include <mpi.h>
#include "ff++.hpp"
#include "AFunction_ext.hpp"
#include "Mesh2.h"
#include "BamgFreeFem.hpp"
#include <vector>
#include <unordered_map>
#include <map>
#include <set>
#include <functional>
#include <climits>
#include <algorithm>

using Fem2D::Mesh;
using Fem2D::MeshPoint;
typedef Mesh const *pmesh;
typedef void *pcommworld;

class AdaptationDDM : public E_F0mps {
 public:
  typedef pmesh Result;

  static basicAC_F0::name_and_type name_param[];
  static const int n_name_param = 32;

  int nbsol;
  Expression nargs[n_name_param];
  Expression getmesh;
  int typesol[1000];
  vector<Expression> sol;
  Expression em11, em22, em12;    // metric=[m11,m12,m22]
  int nbcperiodic;    // periodic=
  Expression *periodic;

  double arg(int i, Stack stack, double a) const { return nargs[i] ? GetAny<double>((*nargs[i])(stack)) : a; }
  long arg(int i, Stack stack, long a) const { return nargs[i] ? GetAny<long>((*nargs[i])(stack)) : a; }
  long *arg(int i, Stack stack, long *a) const { return nargs[i] ? GetAny<long *>((*nargs[i])(stack)) : a; }
  bool arg(int i, Stack stack, bool a) const { return nargs[i] ? GetAny<bool>((*nargs[i])(stack)) : a; }

  AdaptationDDM(const basicAC_F0 &args) : nbsol(args.size() - 1), sol(args.size() - 1) {
    args.SetNameParam(n_name_param, name_param, nargs);
    getmesh = to<pmesh>(args[0]);
    int ksol = 0;
    ffassert(nbsol < (int)(sizeof(typesol) / sizeof(int)));
    for (int i = 1; i < nbsol + 1; i++)
      if (args[i].left() == atype<E_Array>()) {
        const E_Array *a = dynamic_cast<const E_Array *>(args[i].LeftValue());
        ffassert(a);
        ksol += a->size();
      } else
        ksol++;
    sol.resize(ksol);
    ksol = 0;
    for (int i = 1; i < nbsol + 1; i++)
      if (args[i].left() == atype<E_Array>()) {
        const E_Array *a = dynamic_cast<const E_Array *>(args[i].LeftValue());
        ffassert(a);
        int N = a->size();
        typesol[i - 1] = N - 1;
        if (N <= 4) {
          for (int j = 0; j < N; j++) sol[ksol++] = to<double>((*a)[j]);
        } else
          lgerror(" adaptmeshddm : vector with more than 4 components");
      } else {
        typesol[i - 1] = 0;
        sol[ksol++] = to<double>(args[i]);
      }
    // metric= and periodic=: as Adaptation (lgmesh.cpp)
    em11 = em22 = em12 = 0;
    const E_Array *expmetrix = dynamic_cast<const E_Array *>(nargs[24]);
    if (expmetrix) {
      if (expmetrix->nbitem() != 3) ExecError("\nSorry we wait an array with 3 componants in: metrix=[m11,m12,m22]");
      em11 = (*expmetrix)[0];
      em12 = (*expmetrix)[1];
      em22 = (*expmetrix)[2];
      for (int j = 0; j < 3; j++)
        if ((*expmetrix)[j].left() != atype<KN<double> *>())
          CompileError("Sorry the componants in metrix=[m11,m12,m22] must be vectors");
    }
    nbcperiodic = 0;
    periodic = 0;
    GetPeriodic(2, nargs[25], nbcperiodic, periodic);
  }

  static ArrayOfaType typeargs() { return ArrayOfaType(atype<pmesh>(), true); }
  static E_F0 *f(const basicAC_F0 &args) { return new AdaptationDDM(args); }
  AnyType operator()(Stack s) const;
  AnyType Distributed(Stack s) const;    // owned=: distributed mode
  operator aType() const { return atype<pmesh>(); }
};

basicAC_F0::name_and_type AdaptationDDM::name_param[AdaptationDDM::n_name_param] = {
  {"hmin", &typeid(double)},
  {"hmax", &typeid(double)},
  {"err", &typeid(double)},
  {"errg", &typeid(double)},
  {"nbvx", &typeid(long)},    // 4
  {"nbsmooth", &typeid(long)},
  {"nbjacoby", &typeid(long)},
  {"ratio", &typeid(double)},
  {"omega", &typeid(double)},
  {"iso", &typeid(bool)},    // 9
  {"abserror", &typeid(bool)},
  {"cutoff", &typeid(double)},
  {"verbosity", &typeid(long)},
  {"inquire", &typeid(bool)},
  {"splitpbedge", &typeid(bool)},    // 14
  {"maxsubdiv", &typeid(double)},
  {"anisomax", &typeid(double)},
  {"rescaling", &typeid(bool)},
  {"keepbackvertices", &typeid(bool)},
  {"IsMetric", &typeid(bool)},    // 19
  {"power", &typeid(double)},    // 20
  {"thetamax", &typeid(double)},
  {"splitin2", &typeid(bool)},
  {"nomeshgeneration", &typeid(bool)},
  {"metric", &typeid(E_Array)},    // 24
  {"periodic", &typeid(E_Array)},    // 25
  {"requirededges", &typeid(KN_<long>)},    // 26
  {"warning", &typeid(long *)},    // 27
  {"comm", &typeid(pcommworld)},    // 28
  {"nbpass", &typeid(long)},    // 29
  {"nblayer", &typeid(long)},    // 30
  {"owned", &typeid(KN<double> *)}    // 31: distributed mode
};

// Partition of the triangles into np parts by recursive bisection of the
// barycentres, in the frame rotated by the angle theta.
static void RCB(const vector<double> &u, const vector<double> &v, int *idx, int n, int np, int p0,
                vector<int> &part) {
  if (np <= 1) {
    for (int i = 0; i < n; i++) part[idx[i]] = p0;
    return;
  }
  double umin = 1e300, umax = -1e300, vmin = 1e300, vmax = -1e300;
  for (int i = 0; i < n; i++) {
    umin = min(umin, u[idx[i]]), umax = max(umax, u[idx[i]]);
    vmin = min(vmin, v[idx[i]]), vmax = max(vmax, v[idx[i]]);
  }
  const vector<double> &w = (umax - umin >= vmax - vmin) ? u : v;
  int np1 = np / 2;
  int n1 = (int)((long)n * np1 / np);
  std::nth_element(idx, idx + n1, idx + n, [&w](int a, int b) { return w[a] < w[b] || (w[a] == w[b] && a < b); });
  RCB(u, v, idx, n1, np1, p0, part);
  RCB(u, v, idx + n1, n - n1, np - np1, p0 + np1, part);
}

// The triangles outside act (if act is given) have part = -1 (frozen).
static void Partition(const Mesh &Th, int np, double theta, vector<int> &part, const vector<char> *act = 0) {
  int nt = Th.nt;
  vector<double> u(nt), v(nt);
  vector<int> idx;
  double c = cos(theta), s = sin(theta);
  part.assign(nt, -1);
  for (int k = 0; k < nt; k++)
    if (!act || (*act)[k]) {
      R2 G = (Th[k][0] + Th[k][1] + Th[k][2]) / 3.;
      u[k] = c * G.x + s * G.y;
      v[k] = -s * G.x + c * G.y;
      idx.push_back(k);
    }
  RCB(u, v, idx.data(), idx.size(), np, 0, part);
}

// Connected components (through the vertices) of the triangles of act, given
// entirely to the processes, the largest first to the least loaded.
static void PartitionComponents(const Mesh &Th, int np, const vector<char> &act, vector<int> &part) {
  int nt = Th.nt;
  part.assign(nt, -1);
  vector<vector<int> > vt(Th.nv);    // active triangles of each vertex
  for (int k = 0; k < nt; k++)
    if (act[k])
      for (int j = 0; j < 3; j++) vt[Th(k, j)].push_back(k);
  vector<int> comp(nt, -1);
  vector<vector<int> > comps;
  for (int k0 = 0; k0 < nt; k0++)
    if (act[k0] && comp[k0] < 0) {
      int c = comps.size();
      comps.push_back(vector<int>(1, k0));
      comp[k0] = c;
      for (size_t q = 0; q < comps[c].size(); q++) {
        int k = comps[c][q];
        for (int j = 0; j < 3; j++)
          for (int kk : vt[Th(k, j)])
            if (comp[kk] < 0) comp[kk] = c, comps[c].push_back(kk);
      }
    }
  vector<int> order(comps.size());
  for (size_t c = 0; c < comps.size(); c++) order[c] = c;
  std::stable_sort(order.begin(), order.end(),
                   [&comps](int a, int b) { return comps[a].size() > comps[b].size(); });
  vector<long> load(np, 0);
  for (int c : order) {
    int r = std::min_element(load.begin(), load.end()) - load.begin();
    load[r] += comps[c].size();
    for (int k : comps[c]) part[k] = r;
  }
}

// Search of the vertices by coordinates, up to a tolerance (grid).
struct PointHash {
  double tol;
  std::unordered_map<long long, vector<int> > cells;
  vector<R2> pts;
  PointHash(double t) : tol(t) {}
  long long key(long long i, long long j) const { return i * 2147483647LL + j; }
  int find(const R2 &P) const {
    long long i = (long long)floor(P.x / tol), j = (long long)floor(P.y / tol);
    for (long long di = -1; di <= 1; di++)
      for (long long dj = -1; dj <= 1; dj++) {
        auto it = cells.find(key(i + di, j + dj));
        if (it != cells.end())
          for (int k : it->second)
            if (fabs(pts[k].x - P.x) <= tol && fabs(pts[k].y - P.y) <= tol) return k;
      }
    return -1;
  }
  int add(const R2 &P) {
    int k = pts.size();
    pts.push_back(P);
    cells[key((long long)floor(P.x / tol), (long long)floor(P.y / tol))].push_back(k);
    return k;
  }
};

// msh2bamg of src/femlib/BamgFreeFem.cpp (v4.15), plus: the ends of the
// interface edges (label labi) are corners of the geometry. Otherwise bamg
// computes a tangent through the point where the interface meets the boundary of
// the domain and bends the boundary (vertices added off the original boundary).
// With ndfv/ndfe (BuildPeriodic2), the periodic edges are declared
// equivalent as in the periodic variant of msh2bamg.
// Thg, l2g: global mesh and global numbers of the vertices of Th. Where an
// interface meets the boundary of the domain, the tangent of the boundary is the one
// bamg would compute on the global mesh (curved boundaries, thetamax).
static bamg::Triangles *msh2bamgDDM(const Mesh &Th, double cutoffradian, long *reqedgeslab, int nreqedgeslab,
                                    long labi, const Mesh *Thg, const vector<int> *l2g, int nbdfe = 0,
                                    int *ndfv = 0, int *ndfe = 0) {
  using namespace bamg;
  KN<int> equiedges(Th.neb);
  for (int i = 0; i < Th.neb; i++) equiedges[i] = 2 * i;
  if (nbdfe != 0) {
    KN<int> kk(Th.neb), kn(Th.neb);
    kk = 0;
    for (int i = 0; i < Th.neb; i++) {
      int df = ndfe[i];
      kk[df]++;
      if (kk[df] == 1)
        kn[df] = i;
      else {
        int k = kn[df], sens = 0;
        int di0 = ndfv[Th(Th.bedges[i][0])], di1 = ndfv[Th(Th.bedges[i][1])];
        int dk0 = ndfv[Th(Th.bedges[k][0])], dk1 = ndfv[Th(Th.bedges[k][1])];
        if ((di0 == dk0) && (di1 == dk1))
          sens = 0;
        else if ((di1 == dk0) && (di0 == dk1))
          sens = 1;
        else
          ExecError("adaptmeshddm : bug periodic mesh");
        equiedges[i] = 2 * k + sens;
      }
    }
  }
  Triangles *Tn = new Triangles(Th.nv);
  Tn->nbv = Th.nv;
  Tn->nbt = Th.nt;
  Tn->nbe = Th.neb;
  Tn->name = new char[strlen("msh2bamg") + 1];
  strcpy(Tn->name, "msh2bamg");
  throwassert(Tn->triangles);
  Tn->edges = new bamg::Edge[Th.neb];

  Int4 i;
  Metric Mid(1.);
  for (i = 0; i < Th.nv; i++) {
    Tn->vertices[i].r.x = Th(i).x;
    Tn->vertices[i].r.y = Th(i).y;
    Tn->vertices[i].m = Mid;
    Tn->vertices[i].ReferenceNumber = Th(i).lab;
  }
  for (i = 0; i < Th.nt; i++) {
    int i1 = Th(Th[i][0]);
    int i2 = Th(Th[i][1]);
    int i3 = Th(Th[i][2]);
    Tn->triangles[i] = bamg::Triangle(Tn, i1, i2, i3);
    Tn->triangles[i].color = Th[i].lab;
  }
  set<int> labreq;
  for (int i = 0; i < nreqedgeslab; ++i) labreq.insert(reqedgeslab[i]);
  static bamg::GeometricalEdge paszero;    // required edge
  for (i = 0; i < Th.neb; i++) {
    Tn->edges[i].v[0] = Tn->vertices + Th(Th.bedges[i][0]);
    Tn->edges[i].v[1] = Tn->vertices + Th(Th.bedges[i][1]);
    Tn->edges[i].ref = Th.bedges[i].lab;
    Tn->edges[i].on = 0;
    if (labreq.find(Tn->edges[i].ref) != labreq.end()) Tn->edges[i].on = &paszero;
  }
  if (nbdfe != 0)
    Tn->ConsGeometry(cutoffradian, (int *)equiedges);
  else
    Tn->ConsGeometry(cutoffradian);
  Geometry &Gh(Tn->Gh);
  std::set<const GeometricalVertex *> vifc;    // ends of the interface edges
  for (i = 0; i < Gh.nbe; i++)
    if (Gh.edges[i].ref == labi) {
      Gh.edges[i].v[0]->SetCorner(), Gh.edges[i].v[1]->SetCorner();
      vifc.insert(Gh.edges[i].v[0]), vifc.insert(Gh.edges[i].v[1]);
    }
  if (Thg && !vifc.empty()) {
    // global number of a boundary vertex (by its coordinates)
    std::map<std::pair<double, double>, int> num;
    for (i = 0; i < Th.neb; i++)
      for (int j = 0; j < 2; j++) {
        const Fem2D::Vertex &v = Th.bedges[i][j];
        num[std::make_pair(v.x, v.y)] = (*l2g)[Th(v)];
      }
    // global boundary edges of each vertex
    std::map<int, vector<int> > gbe;
    for (int ib = 0; ib < Thg->neb; ib++)
      for (int j = 0; j < 2; j++) gbe[(*Thg)(Thg->bedges[ib][j])].push_back(ib);
    for (i = 0; i < Gh.nbe; i++)
      if (Gh.edges[i].ref != labi)
        for (int jj = 0; jj < 2; jj++) {
          GeometricalVertex *V = Gh.edges[i].v[jj], *W = Gh.edges[i].v[1 - jj];
          if (!vifc.count(V)) continue;
          auto iv = num.find(std::make_pair(V->r.x, V->r.y)), iw = num.find(std::make_pair(W->r.x, W->r.y));
          if (iv == num.end() || iw == num.end()) continue;
          auto ig = gbe.find(iv->second);
          if (ig == gbe.end() || ig->second.size() != 2) continue;    // global corner (order != 2)
          int U = -1;
          for (int ib : ig->second) {
            int a = (*Thg)(Thg->bedges[ib][0]), b = (*Thg)(Thg->bedges[ib][1]);
            int o = a == iv->second ? b : a;
            if (o != iw->second) U = o;
          }
          if (U < 0) continue;
          bamg::R2 PU((*Thg)(U).x, (*Thg)(U).y), PV = V->r, PW = W->r;
          bamg::R2 A = PV - PU, B = PW - PV;
          double c = (A, B) / (Norme2(A) * Norme2(B));
          if (acos(std::max(-1., std::min(1., c))) >= Gh.MaximalAngleOfCorner) continue;    // global corner
          Gh.edges[i].tg[jj] = PW - PU;    // tangent of the global mesh
        }
  }
  Tn->Gh.AfterRead();
  Tn->SetIntCoor();
  Tn->FillHoleInMesh();
  return Tn;
}

// Submesh of the triangles of in (vertices in the global order); the edges
// towards the outside of in have the label labi.
static Mesh *SubMesh(const Mesh &Th, const vector<char> &in, long labi, vector<int> &l2g) {
  vector<int> g2l(Th.nv, -1);
  int nt = 0;
  l2g.clear();
  for (int k = 0; k < Th.nt; k++)
    if (in[k]) {
      nt++;
      for (int j = 0; j < 3; j++) g2l[Th(k, j)] = 1;
    }
  for (int iv = 0; iv < Th.nv; iv++)
    if (g2l[iv] > 0) g2l[iv] = l2g.size(), l2g.push_back(iv);
  vector<int> be, belab;
  for (int i = 0; i < Th.neb; i++) {
    int e, k = Th.BoundaryElement(i, e);
    int ee = e, kk = Th.ElementAdj(k, ee);
    if (kk < 0) kk = k;
    if (in[k] && in[kk])
      be.push_back(g2l[Th(Th.bedges[i][0])]), be.push_back(g2l[Th(Th.bedges[i][1])]), belab.push_back(Th.bedges[i].lab);
  }
  for (int k = 0; k < Th.nt; k++)
    if (in[k])
      for (int e = 0; e < 3; e++) {
        int ee = e, kk = Th.ElementAdj(k, ee);
        if (kk < 0 || kk == k || in[kk]) continue;
        be.push_back(g2l[Th(k, (e + 1) % 3)]), be.push_back(g2l[Th(k, (e + 2) % 3)]), belab.push_back(labi);
      }
  int nv = l2g.size(), nbe = belab.size();
  Fem2D::Vertex *v = new Fem2D::Vertex[nv];
  Fem2D::Triangle *t = new Fem2D::Triangle[nt];
  Fem2D::BoundaryEdge *b = new Fem2D::BoundaryEdge[nbe];
  for (int i = 0; i < nv; i++) v[i].x = Th(l2g[i]).x, v[i].y = Th(l2g[i]).y, v[i].lab = Th(l2g[i]).lab;
  for (int k = 0, kl = 0; k < Th.nt; k++)
    if (in[k]) t[kl++].set(v, g2l[Th(k, 0)], g2l[Th(k, 1)], g2l[Th(k, 2)], Th[k].lab);
  for (int i = 0; i < nbe; i++) b[i].set(v, be[2 * i], be[2 * i + 1], belab[i]);
  return new Mesh(nv, nt, nbe, v, t, b);
}

// Adapted submesh of a process, as arrays.
struct Piece {
  vector<double> xy;    // 2*nv
  vector<int> vlab, iface;    // nv: label, interface vertex (0/1)
  vector<int> tri, tlab;    // 3*nt, nt
  vector<int> bed, blab;    // 2*nbe, nbe (edges to keep)
  vector<int> bidx;    // pass 1: vertices of the input mesh and
  vector<double> bmet;    // their metric (a11, a21, a22)
};

// Error on one process -> error on all (otherwise the others wait).
static void CheckErrors(MPI_Comm comm, int mpisize, int lerr, const string &msg) {
  int gerr = lerr;
  if (mpisize > 1) MPI_Allreduce(&lerr, &gerr, 1, MPI_INT, MPI_MAX, comm);
  if (gerr) {
    if (lerr) cerr << " adaptmeshddm : " << msg << endl;
    ExecError("adaptmeshddm : error in the adaptation of a subdomain");
  }
}

// Gathers arrays of different sizes on all the processes.
template <class T>
static void AllGather(MPI_Comm comm, int mpisize, const vector<T> &loc, vector<T> &all, vector<int> &cnt,
                      vector<int> &dsp, MPI_Datatype type) {
  int n = loc.size();
  cnt.assign(mpisize, n);
  if (mpisize > 1) MPI_Allgather(&n, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
  dsp.assign(mpisize + 1, 0);
  for (int r = 0; r < mpisize; r++) dsp[r + 1] = dsp[r] + cnt[r];
  all.resize(dsp[mpisize]);
  if (mpisize > 1)
    MPI_Allgatherv((void *)loc.data(), n, type, all.data(), cnt.data(), dsp.data(), type, comm);
  else
    std::copy(loc.begin(), loc.end(), all.begin());
}

AnyType AdaptationDDM::operator()(Stack stack) const {
  if (nargs[31]) return Distributed(stack);
  using namespace bamg;
  using bamg::Abs;
  using bamg::Max;
  using bamg::Min;

  Real8 err = arg(2, stack, 0.01);
  Real8 errg = arg(3, stack, Min(0.01, err));
  long nbsx = Max(100L, arg(4, stack, 9000L));
  long nbsmooth = arg(5, stack, 3L);
  long nbjacobi = arg(6, stack, 0L);
  const Real8 raison = arg(7, stack, 1.8);
  const Real8 omega = arg(8, stack, 1.0);
  bool iso = arg(9, stack, false);
  bool AbsError = arg(10, stack, true);
  Real8 CutOff = arg(11, stack, 1.0e-6);
  verbosity = arg(12, stack, (long)verbosity);
  bool SplitEdgeWith2Boundary = arg(14, stack, true);
  double maxsubdiv = Max(Min(arg(15, stack, 10.0), 10.0), 0.1);
  double anisomax = Max((double)arg(16, stack, 1.0e6), 1.0);
  bool rescaling = arg(17, stack, true);
  bool KeepBackVertices = arg(18, stack, true);
  int givenmetric = arg(19, stack, false);
  double powerM = arg(20, stack, 1.0);
  double cutoffradian = arg(21, stack, -1.0) * bamg::Pi / 180.;
  bool split = arg(22, stack, false);
  bool nomeshgeneration = arg(23, stack, false);
  long lwarning = 0;
  long *pwarning = arg(27, stack, &lwarning);
  long &warning = *pwarning;
  long nbpass = Max(1L, arg(29, stack, 3L));
  long nblayer = Max(0L, arg(30, stack, 6L));
  // inquire= (13): bamg drawing, no effect without DRAWING, as adaptmesh

  pcommworld pcomm = nargs[28] ? GetAny<pcommworld>((*nargs[28])(stack)) : nullptr;
  int mpiinit = 0;
  MPI_Initialized(&mpiinit);
  MPI_Comm comm = pcomm ? *(MPI_Comm *)pcomm : MPI_COMM_WORLD;
  int mpirank = 0, mpisize = 1;
  if (mpiinit) {
    MPI_Comm_rank(comm, &mpirank);
    MPI_Comm_size(comm, &mpisize);
  }

  KN<long> reqedges0;
  KN<long> reqedges(nargs[26] ? GetAny<KN_<long> >((*nargs[26])(stack)) : (KN_<long>)reqedges0);

  const Mesh *Thh = GetAny<pmesh>((*getmesh)(stack));
  ffassert(Thh);

  // label of the interface edges: largest boundary label + 1
  long labi = 0;
  for (int i = 0; i < Thh->neb; i++) labi = max(labi, (long)Thh->bedges[i].lab);
  for (int i = 0; i < reqedges.N(); i++) labi = max(labi, reqedges[i]);
  labi++;
  KN<long> reqlab(reqedges.N() + 1);
  for (int i = 0; i < reqedges.N(); i++) reqlab[i] = reqedges[i];
  reqlab[reqedges.N()] = labi;

  // splitin2: as adaptmesh (SplitElement of bamg) on one process, unless
  // the mesh has an internal boundary: SplitElement then crashes (adaptmesh
  // too, bamg bug); the global mesh is split at the end, as in parallel
  bool internalbe = false;
  for (int i = 0; i < Thh->neb && !internalbe; i++) {
    int e, k = Thh->BoundaryElement(i, e);
    int ee = e, kk = Thh->ElementAdj(k, ee);
    internalbe = kk >= 0 && kk != k;
  }
  const bool splitbamg = split && mpisize == 1 && !internalbe;

  // metric=[m11,m12,m22]: intersected with the computed metric, then replaced
  // by the final metric at the vertices of the input mesh (as adaptmesh)
  bool mtx = em11 && em22 && em12;
  KN<double> *mm11 = 0, *mm12 = 0, *mm22 = 0;
  if (mtx) {
    mm11 = GetAny<KN<double> *>((*em11)(stack));
    mm22 = GetAny<KN<double> *>((*em22)(stack));
    mm12 = GetAny<KN<double> *>((*em12)(stack));
    if (mm11->N() != Thh->nv || mm22->N() != Thh->nv || mm12->N() != Thh->nv)
      ExecError("The size of3  metrics array must be equal to nb of vertex");
  }
  // periodic=: on one process, as adaptmesh. On several, the periodic
  // boundaries are frozen during the passes, then a last pass
  // adapts again a band around them on a single process, with the periodicity.
  vector<long> plabs;
  for (int i = 0; i < nbcperiodic; i++)
    plabs.push_back(GetAny<long>((*periodic[4 * i])(stack))), plabs.push_back(GetAny<long>((*periodic[4 * i + 2])(stack)));
  const bool perextra = nbcperiodic && mpisize > 1;
  KN<long> reqlabP(reqlab.N() + plabs.size());    // + frozen periodic labels
  for (int i = 0; i < reqlab.N(); i++) reqlabP[i] = reqlab[i];
  for (int i = 0; i < (int)plabs.size(); i++) reqlabP[reqlab.N() + i] = plabs[i];
  // metric of the input mesh to gather (pass 1)
  const bool needmet0 = nbpass > 1 || mtx || nomeshgeneration || perextra;

  // default hmin, hmax: those of adaptmesh on the global mesh
  // (bounding box extended by 5 % on each side, as in bamg)
  Fem2D::R2 Pmin, Pmax;
  Thh->BoundingBox(Pmin, Pmax);
  Real8 extent = 1.1 * Max(Pmax.x - Pmin.x, Pmax.y - Pmin.y);
  Real8 hmax = 0.3 * extent;
  Real8 hmin = 2.0 * extent / MaxICoor;
  hmin = Max(hmin, arg(0, stack, hmin));
  hmax = Min(hmax, arg(1, stack, hmax));
  if (iso) anisomax = 1;

  int ksol = 0;
  for (int i = 0; i < nbsol; i++) ksol += typesol[i] + 1;

  if (verbosity > 2 && mpirank == 0) {
    cout << endl << " \t\t ## adaptmeshddm : nbsol= " << nbsol << ", nbsx = " << nbsx << ", err = " << err;
    cout << ", hmin = " << hmin << ", hmax = " << hmax << endl;
    cout << " \t\t    ratio  = " << raison << ", nbsmooth = " << nbsmooth << ", omega = " << omega
         << ", iso = " << iso << endl;
    cout << " \t\t    AbsError =" << AbsError << ", CutOff = " << CutOff << ", nbjacobi = " << nbjacobi << endl;
    cout << " \t\t    maxsubdiv = " << maxsubdiv << " splitpbedge = " << SplitEdgeWith2Boundary << endl;
    cout << " \t\t    anisomax = " << anisomax << ", rescaling = " << rescaling << ", power = " << powerM
         << ", KeepBackvertices = " << KeepBackVertices << " IsMetric = " << givenmetric << endl;
    cout << " \t\t    mpisize = " << mpisize << ", nbpass = " << nbpass << ", nblayer = " << nblayer
         << ", label interface = " << labi << endl
         << endl;
  }

  MeshPoint *mp(Fem2D::MeshPointStack(stack)), mps = *mp;
  const long verb0 = verbosity;
  const Mesh *cur = Thh;
  long lwarn = 0;
  vector<char> ifacev;    // vertices of cur on the interfaces of the previous pass
  // metric computed in pass 1 at the vertices of the input mesh: the next
  // passes interpolate it (as bamg with its background mesh), instead of
  // computing again a Hessian on an intermediate mesh already finer
  vector<double> gmet0;
  if (nbpass > 1 || perextra) const_cast<Mesh *>(Thh)->MakeQuadTree();
  const int npasses = nbpass + (perextra ? 1 : 0);

  for (int pass = 0; pass < npasses; ++pass) {
    const Mesh &Th(*cur);
    const bool ppass = perextra && pass == nbpass;    // pass of the periodic boundaries
    const KN<long> &rq = (perextra && !ppass) ? reqlabP : reqlab;
    double tm[6];    // timings of the phases (verbosity > 1)
    tm[0] = MPI_Wtime();
    vector<int> part;
    // vertices touching an adapted triangle (to glue the frozen triangles)
    vector<char> vact(Th.nv, 1);
    int nact = Th.nt;
    if (pass == 0)
      Partition(Th, mpisize, 0., part);
    else if (ppass) {
      // band of nblayer layers around the periodic boundaries, on process 0
      vector<char> act(Th.nt, 0);
      vact.assign(Th.nv, 0);
      for (int i = 0; i < Th.neb; i++)
        if (std::find(plabs.begin(), plabs.end(), (long)Th.bedges[i].lab) != plabs.end())
          vact[Th(Th.bedges[i][0])] = vact[Th(Th.bedges[i][1])] = 1;
      for (int l = 0; l <= nblayer; l++) {
        for (int k = 0; k < Th.nt; k++)
          if (vact[Th(k, 0)] || vact[Th(k, 1)] || vact[Th(k, 2)]) act[k] = 1;
        for (int k = 0; k < Th.nt; k++)
          if (act[k]) vact[Th(k, 0)] = vact[Th(k, 1)] = vact[Th(k, 2)] = 1;
      }
      nact = 0;
      part.assign(Th.nt, -1);
      for (int k = 0; k < Th.nt; k++)
        if (act[k]) part[k] = 0, nact++;
      if (nact == 0) break;
    } else {
      // next passes: only a band of nblayer layers of triangles
      // around the previous interfaces is adapted, the rest is frozen
      vector<char> act(Th.nt, 0);
      vact = ifacev;
      for (int l = 0; l <= nblayer; l++) {
        for (int k = 0; k < Th.nt; k++)
          if (vact[Th(k, 0)] || vact[Th(k, 1)] || vact[Th(k, 2)]) act[k] = 1;
        for (int k = 0; k < Th.nt; k++)
          if (act[k]) vact[Th(k, 0)] = vact[Th(k, 1)] = vact[Th(k, 2)] = 1;
      }
      nact = 0;
      for (int k = 0; k < Th.nt; k++) nact += act[k];
      if (nact == 0) {    // no interface (a single process)
        if (perextra) {
          pass = nbpass - 1;    // -> pass of the periodic boundaries
          continue;
        }
        break;
      }
      if (pass == 1)
        Partition(Th, mpisize, bamg::Pi / 4, part, &act);
      else
        PartitionComponents(Th, mpisize, act, part);
    }

    // ---- local submesh. The vertices keep the global order and the boundary
    // edges the order of Th.bedges (bamg depends on the numbering: on one
    // process we get exactly adaptmesh).
    vector<int> g2l(Th.nv, -1), l2g;
    int ntl = 0;
    for (int k = 0; k < Th.nt; k++)
      if (part[k] == mpirank) {
        ntl++;
        for (int j = 0; j < 3; j++) g2l[Th(k, j)] = 1;
      }
    for (int iv = 0; iv < Th.nv; iv++)
      if (g2l[iv] > 0) g2l[iv] = l2g.size(), l2g.push_back(iv);
    // boundary edges of the submesh: boundary of the domain and internal edges of
    // the user, then interfaces (label labi). For an interface edge
    // we keep the original label (-1 if none) and whether this process keeps it.
    typedef std::pair<double, double> XY;
    typedef std::pair<XY, XY> Seg;
    std::map<Seg, std::pair<int, int> > ifaceorig;
    vector<int> be, belab;
    for (int i = 0; i < Th.neb; i++) {
      int e, k = Th.BoundaryElement(i, e);
      int ee = e, kk = Th.ElementAdj(k, ee);
      if (kk < 0) kk = k;
      if (part[k] == mpirank && part[kk] == mpirank)
        be.push_back(g2l[Th(Th.bedges[i][0])]), be.push_back(g2l[Th(Th.bedges[i][1])]),
          belab.push_back(Th.bedges[i].lab);
    }
    for (int k = 0; k < Th.nt; k++)
      if (part[k] == mpirank)
        for (int e = 0; e < 3; e++) {
          int ee = e, kk = Th.ElementAdj(k, ee);
          if (kk < 0 || kk == k || part[kk] == mpirank) continue;
          int i0 = Th(k, (e + 1) % 3), i1 = Th(k, (e + 2) % 3);
          Fem2D::BoundaryEdge *pbe = Th.TheBoundaryEdge(i0, i1);
          XY A(Th(i0).x, Th(i0).y), B(Th(i1).x, Th(i1).y);
          // second: bit 0 = this process keeps the edge, bit 1 = the other side is adapted
          ifaceorig[Seg(min(A, B), max(A, B))] =
            make_pair(pbe ? (int)pbe->lab : -1, (int)(part[kk] < 0 || mpirank < part[kk]) | (part[kk] >= 0 ? 2 : 0));
          be.push_back(g2l[i0]), be.push_back(g2l[i1]), belab.push_back(labi);
        }

    tm[1] = MPI_Wtime();
    // ---- phase 1: bamg conversion and values of the solutions at the vertices
    Mesh *loc = 0, *gloc = 0;
    Triangles *oTh = 0, *gTh = 0;
    double *lessol = 0;
    int nvl = l2g.size();
    // pass 1: the metric is computed on the subdomain extended by nghost
    // layers (gTh, vertices gl2g): the own vertices then have the same
    // neighbourhood as in sequential, hence the same Hessian as adaptmesh; the
    // 20 layers also cover the gradation (SmoothMetric), which goes further
    const int nghost = 20 + nbjacobi;
    vector<int> gl2g;
    // per solution: min, max (of the norm for a vector), |.|max on the
    // own vertices (-> global values) and on the extended domain (s*g)
    vector<double> smin(nbsol, 1e300), smax(nbsol, -1e300), sabs(nbsol, 0.);
    vector<double> sming(nbsol, 1e300), smaxg(nbsol, -1e300), sabsg(nbsol, 0.);
    int lerr = 0;
    string msg;
    if (ntl) try {
        int nbel = belab.size();
        Fem2D::Vertex *vl = new Fem2D::Vertex[nvl];
        Fem2D::Triangle *tl = new Fem2D::Triangle[ntl];
        Fem2D::BoundaryEdge *bl = new Fem2D::BoundaryEdge[nbel];
        for (int i = 0; i < nvl; i++) {
          vl[i].x = Th(l2g[i]).x;
          vl[i].y = Th(l2g[i]).y;
          vl[i].lab = Th(l2g[i]).lab;
        }
        for (int k = 0, kl = 0; k < Th.nt; k++)
          if (part[k] == mpirank) tl[kl++].set(vl, g2l[Th(k, 0)], g2l[Th(k, 1)], g2l[Th(k, 2)], Th[k].lab);
        for (int i = 0; i < nbel; i++) bl[i].set(vl, be[2 * i], be[2 * i + 1], belab[i]);
        if (mpirank) verbosity = 0;    // only process 0 talks
        loc = new Mesh(nvl, ntl, nbel, vl, tl, bl);
        if (nbcperiodic && (mpisize == 1 || ppass)) {    // equivalent periodic edges
          KN<int> ndfv(loc->nv), ndfe(loc->neb);
          for (int i = 0; i < ndfv.N(); ++i) ndfv[i] = i;
          int nbdfv = 0, nbdfe = 0;
          BuildPeriodic2(nbcperiodic, periodic, *loc, stack, nbdfv, ndfv, nbdfe, ndfe);
          *mp = mps;
          oTh = msh2bamgDDM(*loc, cutoffradian, (long *)rq, rq.N(), labi, &Th, &l2g, nbdfe, (int *)ndfv, (int *)ndfe);
        } else
          oTh = msh2bamgDDM(*loc, cutoffradian, (long *)rq, rq.N(), labi, &Th, &l2g);
        oTh->ReMakeTriangleContainingTheVertex();
        ffassert(oTh->nbv == nvl);

        if (pass == 0) {
          vector<char> in(Th.nt, 0), vin(Th.nv, 0);
          for (int k = 0; k < Th.nt; k++) in[k] = part[k] == mpirank;
          for (int l = 0; l < nghost; l++) {
            for (int k = 0; k < Th.nt; k++)
              if (in[k]) vin[Th(k, 0)] = vin[Th(k, 1)] = vin[Th(k, 2)] = 1;
            for (int k = 0; k < Th.nt; k++)
              if (vin[Th(k, 0)] || vin[Th(k, 1)] || vin[Th(k, 2)]) in[k] = 1;
          }
          int ntg = 0;
          for (int k = 0; k < Th.nt; k++) ntg += in[k];
          if (ntg == ntl)    // a single process: no layer added
            gTh = oTh, gl2g = l2g;
          else {
            gloc = SubMesh(Th, in, labi, gl2g);
            gTh = msh2bamgDDM(*gloc, cutoffradian, (long *)rq, rq.N(), labi, &Th, &gl2g);
            gTh->ReMakeTriangleContainingTheVertex();
            ffassert(gTh->nbv == (int)gl2g.size());
          }
          // values of the solutions at the vertices of gTh (numbers of msh2bamg);
          // evaluation on the current global mesh, as adaptmesh
          int nvg = gl2g.size();
          vector<int> gg2l(Th.nv, -1);
          for (int i = 0; i < nvg; i++) gg2l[gl2g[i]] = i;
          lessol = new double[nvg * ksol];
          vector<int> done(nvg, 0);
          for (int k = 0; k < Th.nt; k++)
            if (in[k])
              for (int jt = 0; jt < 3; jt++) {
                int iv = gg2l[Th(k, jt)];
                if (!done[iv]) {
                  done[iv] = 1;
                  mp->setP(&Th, k, jt);
                  double *ss = lessol + ksol * iv;
                  for (int j = 0; j < ksol; j++) *ss++ = GetAny<double>((*sol[j])(stack));
                }
              }
          *mp = mps;
          for (int i = 0, k0 = 0; i < nbsol; k0 += typesol[i] + 1, i++)
            for (int iv = 0; iv < nvg; iv++) {
              double *s = lessol + ksol * iv + k0, v = s[0];
              if (typesol[i]) {
                v = 0;
                for (int j = 0; j <= typesol[i]; j++) v += s[j] * s[j];
                v = sqrt(v);
              }
              sming[i] = Min(sming[i], v), smaxg[i] = Max(smaxg[i], v), sabsg[i] = Max(sabsg[i], Abs(v));
              if (g2l[gl2g[iv]] >= 0)
                smin[i] = Min(smin[i], v), smax[i] = Max(smax[i], v), sabs[i] = Max(sabs[i], Abs(v));
            }
        }
      } catch (std::exception &e) {
        lerr = 1, msg = e.what();
      } catch (...) {
        lerr = 1, msg = "unknown error";
      }
    verbosity = verb0;
    *mp = mps;
    CheckErrors(comm, mpisize, lerr, msg);

    // global min/max of the solutions (normalization of the metric)
    vector<double> gmin(smin), gmax(smax), gabs(sabs);
    if (mpisize > 1 && nbsol) {
      MPI_Allreduce(smin.data(), gmin.data(), nbsol, MPI_DOUBLE, MPI_MIN, comm);
      MPI_Allreduce(smax.data(), gmax.data(), nbsol, MPI_DOUBLE, MPI_MAX, comm);
      MPI_Allreduce(sabs.data(), gabs.data(), nbsol, MPI_DOUBLE, MPI_MAX, comm);
    }

    tm[2] = MPI_Wtime();
    // ---- phase 2: metric, adapted mesh of the subdomain
    Piece mine;
    double lmass = 0;
    if (ntl) try {
        if (mpirank) verbosity = 0;
        Triangles &BTh(*oTh);
        if (pass > 0) {    // metric of pass 1 interpolated in the input mesh
          for (int iv = 0; iv < BTh.nbv; iv++) {
            Fem2D::R2 Pv(BTh[iv].r.x, BTh[iv].r.y), Phat;
            bool outside;
            const Fem2D::Triangle *K = Thh->Find(Pv, Phat, outside);
            int k = Thh->number(K);
            Real8 l[3] = {1 - Phat.x - Phat.y, Phat.x, Phat.y};
            MetricAnIso M[3];
            for (int j = 0; j < 3; j++) {
              const double *m = &gmet0[3 * (*Thh)(k, j)];
              M[j] = MetricAnIso(m[0], m[1], m[2]);
            }
            BTh[iv].m = MetricAnIso(l, M[0], M[1], M[2]);
          }
          BTh.SmoothMetric(raison);    // gradation across the old interfaces
        } else {
          Triangles &GTh(*gTh);    // extended domain
          int nvg = gl2g.size();
          Metric Mhmax(hmax);
          for (int iv = 0; iv < GTh.nbv; iv++) GTh[iv].m = Mhmax;
          if (mtx)
            for (int iv = 0; iv < GTh.nbv; iv++) {
              int g = gl2g[iv];
              GTh[iv].m.IntersectWith(MetricAnIso((*mm11)[g], (*mm12)[g], (*mm22)[g]));
            }
          if (givenmetric) {
            if (ksol == 1) {
              for (int iv = 0; iv < GTh.nbv; iv++) GTh[iv].m.IntersectWith(Metric(lessol[iv]));
            } else if (ksol == 3) {
              for (int iv = 0, k = 0; iv < GTh.nbv; iv++, k += 3) {
                Metric MM(lessol[k], lessol[k + 1], lessol[k + 2]);
                MatVVP2x2 vp(MM);
                vp.Abs();
                GTh[iv].m.IntersectWith(vp);
              }
            } else
              ExecError("adaptmeshddm : IsMetric and ksol != 1 or 3");
          } else {
            // IntersectConsMetric of bamg, one solution at a time.
            // In absolute error with rescaling, bamg divides the Hessian by
            // the local range dl = max - min of the solution on the subdomain:
            // we pass coef = sqrt(err * dg / dl) to get the division by
            // the global range dg, as adaptmesh (the solution is not changed).
            // In relative error, cutoff is relative to the global max.
            for (int i = 0, k0 = 0; i < nbsol; k0 += typesol[i] + 1, i++) {
              int nf = typesol[i] + 1;
              double dl = smaxg[i] - sming[i], dg = gmax[i] - gmin[i];
              double coef = sqrt(err), cut = 0;
              if (AbsError) {
                if (dl <= 1e-10 * Max(sabsg[i], 1e-20)) continue;    // locally constant
                if (rescaling) coef = sqrt(err * dg / dl);
              } else {
                if (sabsg[i] <= 0) continue;
                cut = CutOff * gabs[i] / sabsg[i];
              }
              vector<double> s1(nvg * nf);
              for (int iv = 0; iv < nvg; iv++)
                for (int j = 0; j < nf; j++) s1[iv * nf + j] = lessol[iv * ksol + k0 + j];
              GTh.IntersectConsMetric(s1.data(), 1, typesol + i, hmin, hmax, coef, anisomax, cut, nbjacobi,
                                      rescaling, powerM, 0);
            }
          }
          // all the rest of adaptmesh on the extended domain, then copy
          GTh.IntersectGeomMetric(errg, iso);
          GTh.SmoothMetric(raison);
          GTh.MaxSubDivision(maxsubdiv);
          GTh.BoundAnisotropy(anisomax);
          if (gTh != oTh) {    // metric of the own vertices
            vector<int> gg2l(Th.nv, -1);
            for (int i = 0; i < nvg; i++) gg2l[gl2g[i]] = i;
            for (int iv = 0; iv < BTh.nbv; iv++) BTh[iv].m = GTh[gg2l[l2g[iv]]].m;
          }
          if (needmet0)
            for (int iv = 0; iv < BTh.nbv; iv++) {
              const MetricAnIso &m = BTh[iv].m;
              mine.bidx.push_back(l2g[iv]);
              mine.bmet.push_back(m.a11), mine.bmet.push_back(m.a21), mine.bmet.push_back(m.a22);
            }
        }
        // number of vertices predicted by the metric: integral of sqrt(det M)
        for (int k = 0; k < loc->nt; k++) {
          double d = 0;
          for (int j = 0; j < 3; j++) {
            const MetricAnIso &m = BTh[(*loc)(k, j)].m;
            d += sqrt(Max(m.a11 * m.a22 - m.a21 * m.a21, 0.)) / 3;
          }
          lmass += (*loc)[k].area * d;
        }
      } catch (std::exception &e) {
        lerr = 1, msg = e.what();
      } catch (...) {
        lerr = 1, msg = "unknown error";
      }
    verbosity = verb0;
    CheckErrors(comm, mpisize, lerr, msg);

    // metric of the input mesh (pass 1): metric= on output, next passes
    if (pass == 0 && needmet0) {
      vector<int> abidx, ci0, di0;
      vector<double> abmet;
      AllGather(comm, mpisize, mine.bidx, abidx, ci0, di0, MPI_INT);
      AllGather(comm, mpisize, mine.bmet, abmet, ci0, di0, MPI_DOUBLE);
      gmet0.assign(3 * Thh->nv, 0.);
      vector<char> seen(Thh->nv, 0);
      for (size_t i = 0; i < abidx.size(); i++)
        if (!seen[abidx[i]]) {    // interface vertex: metric of the first process
          seen[abidx[i]] = 1;
          for (int j = 0; j < 3; j++) gmet0[3 * abidx[i] + j] = abmet[3 * i + j];
        }
      if (mtx)
        for (int i = 0; i < Thh->nv; i++)
          (*mm11)[i] = gmet0[3 * i], (*mm12)[i] = gmet0[3 * i + 1], (*mm22)[i] = gmet0[3 * i + 2];
      if (nomeshgeneration) {    // metric only, as adaptmesh
        delete[] lessol;
        if (gTh != oTh) delete gTh;
        delete oTh;
        delete loc;
        delete gloc;
        return SetAny<pmesh>(Thh);
      }
    }

    // nbvx limits the total number of vertices: the budget (nbvx - frozen vertices)
    // is shared in proportion to the number of vertices predicted by the metric, with
    // a floor per process = vertices of the background submesh + 10 (bamg
    // can keep them all, keepbackvertices, and overflows if nbvx is too
    // small); the budget of the processes at the floor is removed from the others
    long nbsxl = nbsx;
    if (mpisize > 1) {
      double lmf[2] = {lmass, ntl ? nvl + 10. : 0.};
      vector<double> amf(2 * mpisize);
      MPI_Allgather(lmf, 2, MPI_DOUBLE, amf.data(), 2, MPI_DOUBLE, comm);
      long nfix = 0;
      for (int iv = 0; iv < Th.nv; iv++) nfix += !vact[iv];
      vector<double> cap(mpisize, 0.);
      vector<char> fixe(mpisize, 0);
      for (int r = 0; r < mpisize; r++) fixe[r] = amf[2 * r + 1] == 0;    // without triangle
      for (bool change = true; change;) {
        change = false;
        double budget = nbsx - nfix, m = 0;
        for (int r = 0; r < mpisize; r++)
          if (fixe[r])
            budget -= amf[2 * r + 1];
          else
            m += amf[2 * r];
        for (int r = 0; r < mpisize; r++)
          if (!fixe[r]) {
            cap[r] = m > 0 ? budget * amf[2 * r] / m : 0;
            if (cap[r] < amf[2 * r + 1]) fixe[r] = 1, change = true;
          }
      }
      nbsxl = fixe[mpirank] ? (long)amf[2 * mpirank + 1] : Max((long)llround(cap[mpirank]), (long)nvl + 10);
    }

    if (ntl) try {
        if (mpirank) verbosity = 0;
        Triangles &BTh(*oTh);
        Triangles *nTh = new Triangles(nbsxl, BTh, KeepBackVertices);    // Adaptation
        if (splitbamg) nTh->SplitElement(1);    // as adaptmesh (otherwise: at the end)
        if (SplitEdgeWith2Boundary) nTh->SplitInternalEdgeWithBorderVertices();
        if (verbosity > 3) nTh->ShowHistogram();
        if (nbsmooth) nTh->SmoothingVertex(nbsmooth, omega);
        if (verbosity > 2 && nbsmooth) nTh->ShowHistogram();
        if (verbosity > 0) nTh->ShowRegulaty();
        lwarn += nTh->warning;
        const Mesh *a = bamg2msh(nTh, true);
        delete nTh;

        // flattening of the adapted submesh
        int nv = a->nv;
        mine.xy.resize(2 * nv);
        mine.vlab.resize(nv);
        mine.iface.assign(nv, 0);
        for (int i = 0; i < nv; i++) {
          mine.xy[2 * i] = (*a)(i).x;
          mine.xy[2 * i + 1] = (*a)(i).y;
          mine.vlab[i] = (*a)(i).lab;
        }
        for (int k = 0; k < a->nt; k++) {
          for (int j = 0; j < 3; j++) mine.tri.push_back((*a)(k, j));
          mine.tlab.push_back((*a)[k].lab);
        }
        int nmiss = 0;
        for (int i = 0; i < a->neb; i++) {
          int i0 = (*a)(a->bedges[i][0]), i1 = (*a)(a->bedges[i][1]), lab = a->bedges[i].lab;
          if (lab == labi) {
            XY A((*a)(i0).x, (*a)(i0).y), B((*a)(i1).x, (*a)(i1).y);
            auto it = ifaceorig.find(Seg(min(A, B), max(A, B)));
            if (it == ifaceorig.end()) {
              nmiss++;
              continue;
            }
            // iface: 1 = interface with the frozen zone, 2 = between two adapted zones
            int f = (it->second.second & 2) ? 2 : 1;
            mine.iface[i0] = max(mine.iface[i0], f), mine.iface[i1] = max(mine.iface[i1], f);
            if (it->second.first < 0 || !(it->second.second & 1)) continue;
            lab = it->second.first;    // internal edge of the user on the interface
          }
          mine.bed.push_back(i0), mine.bed.push_back(i1), mine.blab.push_back(lab);
        }
        delete a;
        if (nmiss) ExecError("adaptmeshddm : an interface edge was modified by bamg");
      } catch (std::exception &e) {
        lerr = 1, msg = e.what();
      } catch (...) {
        lerr = 1, msg = "unknown error";
      }
    verbosity = verb0;
    delete[] lessol;
    if (gTh != oTh) delete gTh;
    delete oTh;
    delete loc;
    delete gloc;
    CheckErrors(comm, mpisize, lerr, msg);

    tm[3] = MPI_Wtime();
    // ---- phase 3: gathering and gluing by the interface vertices
    vector<int> nums(3), anums;
    nums[0] = mine.vlab.size(), nums[1] = mine.tlab.size(), nums[2] = mine.blab.size();
    vector<int> cnt, dsp, ci, di;
    AllGather(comm, mpisize, nums, anums, cnt, dsp, MPI_INT);
    vector<double> axy;
    vector<int> avlab, aiface, atri, atlab, abed, ablab;
    AllGather(comm, mpisize, mine.xy, axy, cnt, dsp, MPI_DOUBLE);
    AllGather(comm, mpisize, mine.vlab, avlab, ci, di, MPI_INT);
    AllGather(comm, mpisize, mine.iface, aiface, ci, di, MPI_INT);
    AllGather(comm, mpisize, mine.tri, atri, ci, di, MPI_INT);
    AllGather(comm, mpisize, mine.tlab, atlab, ci, di, MPI_INT);
    AllGather(comm, mpisize, mine.bed, abed, ci, di, MPI_INT);
    AllGather(comm, mpisize, mine.blab, ablab, ci, di, MPI_INT);

    tm[4] = MPI_Wtime();
    PointHash ph(1e-10 * extent);
    vector<Fem2D::R2> P;
    vector<int> vlab;
    vector<int> num;    // global number of the vertex i of the process r
    int ov = 0, ot = 0, ob = 0;
    vector<int> tri, tlab, bed, blab;
    int nmerge = 0;
    // frozen triangles: taken from cur by each process, without communication
    if (nact < Th.nt) {
      vector<int> fnum(Th.nv, -1);
      for (int k = 0; k < Th.nt; k++)
        if (part[k] < 0) {
          for (int j = 0; j < 3; j++) {
            int iv = Th(k, j);
            if (fnum[iv] < 0) {
              Fem2D::R2 Q(Th(iv).x, Th(iv).y);
              int g = P.size();
              if (vact[iv])
                ph.add(Q);    // vertex on the boundary of the adapted zone
              else
                ph.pts.push_back(Q);
              P.push_back(Q), vlab.push_back(Th(iv).lab);
              fnum[iv] = g;
            }
            tri.push_back(fnum[iv]);
          }
          tlab.push_back(Th[k].lab);
        }
      for (int i = 0; i < Th.neb; i++) {
        int e, k = Th.BoundaryElement(i, e);
        int ee = e, kk = Th.ElementAdj(k, ee);
        if (kk < 0) kk = k;
        if (part[k] < 0 && part[kk] < 0)
          bed.push_back(fnum[Th(Th.bedges[i][0])]), bed.push_back(fnum[Th(Th.bedges[i][1])]),
            blab.push_back(Th.bedges[i].lab);
      }
    }
    vector<char> nifacev;
    for (int r = 0; r < mpisize; r++) {
      int nv = anums[3 * r], nt = anums[3 * r + 1], nbe = anums[3 * r + 2];
      num.resize(nv);
      for (int i = 0; i < nv; i++) {
        Fem2D::R2 Q(axy[2 * (ov + i)], axy[2 * (ov + i) + 1]);
        int lab = avlab[ov + i], g = -1;
        if (aiface[ov + i]) {
          g = ph.find(Q);
          if (g >= 0) {
            nmerge++;
            if (vlab[g] == labi) vlab[g] = lab;
          } else {
            g = ph.add(Q);
            ffassert(g == (int)P.size());
          }
        } else {
          g = P.size();
          ph.pts.push_back(Q);    // keeps the numbering of ph = that of P
        }
        if (g == (int)P.size()) P.push_back(Q), vlab.push_back(lab);
        num[i] = g;
        if (aiface[ov + i] == 2) {    // the next pass adapts again around the cuts
          if ((int)nifacev.size() <= g) nifacev.resize(g + 1, 0);
          nifacev[g] = 1;
        }
      }
      for (int k = 0; k < nt; k++) {
        for (int j = 0; j < 3; j++) tri.push_back(num[atri[3 * (ot + k) + j]]);
        tlab.push_back(atlab[ot + k]);
      }
      for (int i = 0; i < nbe; i++) {
        bed.push_back(num[abed[2 * (ob + i)]]), bed.push_back(num[abed[2 * (ob + i) + 1]]);
        blab.push_back(ablab[ob + i]);
      }
      ov += nv, ot += nt, ob += nbe;
    }
    int nv = P.size(), nt = tlab.size(), nbe = blab.size();
    Fem2D::Vertex *v = new Fem2D::Vertex[nv];
    Fem2D::Triangle *t = new Fem2D::Triangle[nt];
    Fem2D::BoundaryEdge *b = new Fem2D::BoundaryEdge[nbe];
    for (int i = 0; i < nv; i++) {
      v[i].x = P[i].x;
      v[i].y = P[i].y;
      v[i].lab = vlab[i] == labi ? 0 : vlab[i];
    }
    for (int k = 0; k < nt; k++) t[k].set(v, tri[3 * k], tri[3 * k + 1], tri[3 * k + 2], tlab[k]);
    for (int i = 0; i < nbe; i++) b[i].set(v, bed[2 * i], bed[2 * i + 1], blab[i]);
    if (mpirank) verbosity = 0;
    Mesh *g = new Mesh(nv, nt, nbe, v, t, b);
    verbosity = verb0;
    tm[5] = MPI_Wtime();
    nifacev.resize(nv, 0);
    ifacev.swap(nifacev);
    if (verbosity > 1 && mpirank == 0)
      cout << "  -- adaptmeshddm pass " << pass << " : " << nact << " / " << Th.nt << " triangles adapted -> "
           << nt << " triangles, " << nv << " vertices, " << nmerge << " interface vertices merged" << endl
           << "     time: partition " << tm[1] - tm[0] << " s, metric " << tm[2] - tm[1] << " s, mesh "
           << tm[3] - tm[2] << " s, gather " << tm[4] - tm[3] << " s, merge " << tm[5] - tm[4] << " s" << endl;
    if (cur != Thh) delete cur;
    cur = g;
  }

  Mesh *res = const_cast<Mesh *>(cur);
  if (split && !splitbamg) {    // split of each triangle into 4, on the global mesh
    KN<int> s(res->nt);
    s = 2;
    Mesh *r2 = new Mesh(*res, (int *)s, false, 0);
    if (res != Thh) delete res;
    res = r2;
  }
  if (res == Thh) return SetAny<pmesh>(Thh);
  double tr = MPI_Wtime();
  if (mpisize > 1) res->renum();
  res->MakeQuadTree();
  if (verbosity > 1 && mpirank == 0) cout << "     time: renum + quadtree " << MPI_Wtime() - tr << " s" << endl;
  long gwarn = lwarn;
  if (mpisize > 1) MPI_Allreduce(&lwarn, &gwarn, 1, MPI_LONG, MPI_SUM, comm);
  warning = gwarn;
  Add2StackOfPtr2FreeRC(stack, res);
  return SetAny<pmesh>(res);
}

// =====================================================================
// Distributed mode (as parmmg3d with nodeCommunicators): owned= given.
// Input: LOCAL mesh with overlap (MatCreate / DmeshCreate) and
//   owned = P0 array (nonzero = triangle owned by this process, PartitionCreate).
// Output: adapted local mesh WITHOUT overlap, interface edges with
//   label -111111, to pass to DmeshReconstruct. There is never a global mesh.
//   pass 1: each process remeshes its own triangles (interfaces
//             frozen), metric computed on the whole local mesh (overlap
//             included);
//   pass 2: band around the interfaces, cut by bisection rotated by pi/4,
//             triangles sent to the process which remeshes it (MPI_Alltoallv), the
//             metric is carried with the vertices;
//   pass 3: band around the cuts of pass 2, connected components (even
//             spanning several processes) given entirely to one process.
// =====================================================================

// mesh piece of a process (+ metric at the vertices)
struct DPiece {
  vector<double> xy, met;    // 2 and 3 per vertex
  vector<int> vlab, tri, tlab;
  vector<int> bed, blab;    // edges of the user (boundary of the domain, internal)
  int nv() const { return (int)vlab.size(); }
  int nt() const { return (int)tlab.size(); }
  int addv(double x, double y, int lab, const double *m) {
    xy.push_back(x), xy.push_back(y), vlab.push_back(lab);
    met.push_back(m[0]), met.push_back(m[1]), met.push_back(m[2]);
    return nv() - 1;
  }
};
typedef std::pair<int, int> EKey;
static inline EKey Ek(int a, int b) { return a < b ? EKey(a, b) : EKey(b, a); }
typedef std::pair<double, double> DXY;
typedef std::pair<DXY, DXY> DSeg;
static inline DSeg Sg(double x0, double y0, double x1, double y1) {
  DXY A(x0, y0), B(x1, y1);
  return A < B ? DSeg(A, B) : DSeg(B, A);
}

// FreeFEM Mesh of a piece. Boundary edges: those of the user, except the
// internal edges (label in intlabs) on the boundary of the piece, which get labb (frozen)
// and are recorded in keep to recover their label; the other boundary edges
// get labb.
static Mesh *PieceMesh(const DPiece &P, long labb, const std::set<long> &intlabs, std::map<DSeg, int> *keep) {
  std::map<EKey, int> cnt;
  std::map<EKey, int> first;    // edge -> its first vertex in the counterclockwise order of its triangle
  for (int k = 0; k < P.nt(); k++)
    for (int e = 0; e < 3; e++) {
      int a = P.tri[3 * k + (e + 1) % 3], b = P.tri[3 * k + (e + 2) % 3];
      cnt[Ek(a, b)]++;
      first[Ek(a, b)] = a;
    }
  std::set<EKey> ub;
  vector<int> be, belab;
  for (int i = 0; i < (int)P.blab.size(); i++) {
    int a = P.bed[2 * i], b = P.bed[2 * i + 1];
    EKey k = Ek(a, b);
    if (!ub.insert(k).second) continue;
    // boundary of the piece: the domain on the left (the orientation is lost by Migrate, which stores the edges by
    // sorted vertex numbers; the equivalent periodic edges of bamg need consistent orientations along a side)
    if (cnt[k] == 1 && first[k] != a) std::swap(a, b);
    int lab = P.blab[i];
    if (cnt[k] == 1 && intlabs.count(lab)) {
      if (keep) (*keep)[Sg(P.xy[2 * a], P.xy[2 * a + 1], P.xy[2 * b], P.xy[2 * b + 1])] = lab;
      lab = labb;
    }
    be.push_back(a), be.push_back(b), belab.push_back(lab);
  }
  for (int k = 0; k < P.nt(); k++)
    for (int e = 0; e < 3; e++) {
      int a = P.tri[3 * k + (e + 1) % 3], b = P.tri[3 * k + (e + 2) % 3];
      if (cnt[Ek(a, b)] == 1 && !ub.count(Ek(a, b))) be.push_back(a), be.push_back(b), belab.push_back(labb);
    }
  int nv = P.nv(), nt = P.nt(), nbe = belab.size();
  Fem2D::Vertex *v = new Fem2D::Vertex[nv];
  Fem2D::Triangle *t = new Fem2D::Triangle[nt];
  Fem2D::BoundaryEdge *b = new Fem2D::BoundaryEdge[nbe];
  for (int i = 0; i < nv; i++) v[i].x = P.xy[2 * i], v[i].y = P.xy[2 * i + 1], v[i].lab = P.vlab[i];
  for (int k = 0; k < nt; k++) t[k].set(v, P.tri[3 * k], P.tri[3 * k + 1], P.tri[3 * k + 2], P.tlab[k]);
  for (int i = 0; i < nbe; i++) b[i].set(v, be[2 * i], be[2 * i + 1], belab[i]);
  return new Mesh(nv, nt, nbe, v, t, b);
}

struct GenParams {
  bool kbv, splitpbedge;
  long nbsmooth;
  double omega;
};

// Adapted mesh of a piece: background loc / oTh (same vertex numbers, metric
// already in oTh). Output: piece with the metric interpolated at the new vertices
// in the background; onLabi = vertices on a frozen edge (label labi).
static void GenerateDPiece(Mesh *loc, bamg::Triangles *oTh, long cap, const GenParams &g, long labi,
                           const std::map<DSeg, int> &keep, DPiece &out, vector<char> &onLabi, long &lwarn) {
  using namespace bamg;
  vector<MetricAnIso> m0(oTh->nbv);
  for (int i = 0; i < oTh->nbv; i++) m0[i] = (*oTh)[i].m;
  Triangles *nTh = new Triangles(cap, *oTh, g.kbv);
  if (g.splitpbedge) nTh->SplitInternalEdgeWithBorderVertices();
  if (g.nbsmooth) nTh->SmoothingVertex(g.nbsmooth, g.omega);
  lwarn += nTh->warning;
  const Mesh *a = bamg2msh(nTh, true);
  delete nTh;
  loc->MakeQuadTree();
  out = DPiece();
  for (int i = 0; i < a->nv; i++) {
    Fem2D::R2 P((*a)(i).x, (*a)(i).y), Phat;
    bool outside;
    const Fem2D::Triangle *K = loc->Find(P, Phat, outside);
    int k = loc->number(K);
    Real8 l[3] = {1 - Phat.x - Phat.y, Phat.x, Phat.y};
    MetricAnIso M = MetricAnIso(l, m0[(*loc)(k, 0)], m0[(*loc)(k, 1)], m0[(*loc)(k, 2)]);
    double m[3] = {M.a11, M.a21, M.a22};
    out.addv(P.x, P.y, (*a)(i).lab, m);
  }
  for (int k = 0; k < a->nt; k++) {
    for (int j = 0; j < 3; j++) out.tri.push_back((*a)(k, j));
    out.tlab.push_back((*a)[k].lab);
  }
  onLabi.assign(a->nv, 0);
  for (int i = 0; i < a->neb; i++) {
    int i0 = (*a)(a->bedges[i][0]), i1 = (*a)(a->bedges[i][1]), lab = a->bedges[i].lab;
    if (lab == labi) {
      onLabi[i0] = onLabi[i1] = 1;
      auto it = keep.find(Sg((*a)(i0).x, (*a)(i0).y, (*a)(i1).x, (*a)(i1).y));
      if (it == keep.end()) continue;
      lab = it->second;    // internal edge of the user, frozen
    }
    out.bed.push_back(i0), out.bed.push_back(i1), out.blab.push_back(lab);
  }
  delete a;
  // labels of the vertices: bamg leaves some of them uninitialized (random
  // values, adaptmesh too); here the label of a boundary edge, 0 otherwise
  out.vlab.assign(out.nv(), 0);
  for (int i = 0; i < (int)out.blab.size(); i++)
    for (int j = 0; j < 2; j++)
      if (!out.vlab[out.bed[2 * i + j]]) out.vlab[out.bed[2 * i + j]] = out.blab[i];
}

// Share of nbvx of a process: budget shared in proportion to lmass, floor
// flo (0 = process without triangle), see the global mode.
static long ShareBudget(MPI_Comm comm, int size, int rank, double lmass, double flo, double budget) {
  double lmf[2] = {lmass, flo};
  vector<double> amf(2 * size);
  MPI_Allgather(lmf, 2, MPI_DOUBLE, amf.data(), 2, MPI_DOUBLE, comm);
  vector<double> cap(size, 0.);
  vector<char> fixe(size, 0);
  for (int r = 0; r < size; r++) fixe[r] = amf[2 * r + 1] == 0;
  for (bool change = true; change;) {
    change = false;
    double b = budget, m = 0;
    for (int r = 0; r < size; r++)
      if (fixe[r])
        b -= amf[2 * r + 1];
      else
        m += amf[2 * r];
    for (int r = 0; r < size; r++)
      if (!fixe[r]) {
        cap[r] = m > 0 ? b * amf[2 * r] / m : 0;
        if (cap[r] < amf[2 * r + 1]) fixe[r] = 1, change = true;
      }
  }
  return fixe[rank] ? (long)amf[2 * rank + 1] : std::max((long)llround(cap[rank]), (long)flo);
}

// number of vertices predicted by the metric of a piece (integral of sqrt(det M))
static double PieceMass(const DPiece &P) {
  double s = 0;
  for (int k = 0; k < P.nt(); k++) {
    const double *A = &P.xy[2 * P.tri[3 * k]], *B = &P.xy[2 * P.tri[3 * k + 1]], *C = &P.xy[2 * P.tri[3 * k + 2]];
    double area = 0.5 * fabs((B[0] - A[0]) * (C[1] - A[1]) - (B[1] - A[1]) * (C[0] - A[0])), d = 0;
    for (int j = 0; j < 3; j++) {
      const double *m = &P.met[3 * P.tri[3 * k + j]];
      d += sqrt(std::max(m[0] * m[2] - m[1] * m[1], 0.)) / 3;
    }
    s += area * d;
  }
  return s;
}

// Sends the triangles of P to tgt (>= 0); F = remaining triangles (tgt < 0),
// B = received triangles, glued. seedB = vertices of B on a boundary edge which
// is not against the frozen part (cut between two pieces): seeds of the
// next pass.
static void Migrate(const DPiece &P, const vector<int> &tgt, MPI_Comm comm, int size, double tol, DPiece &F,
                    DPiece &B, std::set<DXY> &seedB) {
  std::map<EKey, std::pair<int, int> > et;    // edge -> 2 triangles
  for (int k = 0; k < P.nt(); k++)
    for (int e = 0; e < 3; e++) {
      auto &p = et.insert(make_pair(Ek(P.tri[3 * k + (e + 1) % 3], P.tri[3 * k + (e + 2) % 3]), make_pair(-1, -1))).first->second;
      (p.first < 0 ? p.first : p.second) = k;
    }
  std::map<EKey, int> ulab;
  for (int i = 0; i < (int)P.blab.size(); i++) ulab[Ek(P.bed[2 * i], P.bed[2 * i + 1])] = P.blab[i];
  // packets: 15 reals and 10 integers per triangle
  vector<vector<double> > sd(size);
  vector<vector<int> > si(size);
  for (int k = 0; k < P.nt(); k++) {
    int t = tgt[k];
    if (t < 0) continue;
    for (int j = 0; j < 3; j++) {
      int v = P.tri[3 * k + j];
      sd[t].push_back(P.xy[2 * v]), sd[t].push_back(P.xy[2 * v + 1]);
      for (int c = 0; c < 3; c++) sd[t].push_back(P.met[3 * v + c]);
      si[t].push_back(P.vlab[v]);
    }
    si[t].push_back(P.tlab[k]);
    for (int e = 0; e < 3; e++) {
      EKey key = Ek(P.tri[3 * k + (e + 1) % 3], P.tri[3 * k + (e + 2) % 3]);
      auto u = ulab.find(key);
      si[t].push_back(u == ulab.end() ? INT_MIN : u->second);
      const auto &p = et[key];
      int nb = p.first == k ? p.second : p.first, flag = 0;
      if (nb >= 0) flag = tgt[nb] < 0 ? 1 : (tgt[nb] != t ? 2 : 3);
      si[t].push_back(flag);
    }
  }
  vector<int> scd(size), sci(size), rcd(size), rci(size), sdd(size + 1, 0), sdi(size + 1, 0), rdd(size + 1, 0),
    rdi(size + 1, 0);
  vector<double> bd;
  vector<int> bi;
  for (int r = 0; r < size; r++) {
    scd[r] = sd[r].size(), sci[r] = si[r].size();
    bd.insert(bd.end(), sd[r].begin(), sd[r].end()), bi.insert(bi.end(), si[r].begin(), si[r].end());
    sdd[r + 1] = sdd[r] + scd[r], sdi[r + 1] = sdi[r] + sci[r];
  }
  MPI_Alltoall(scd.data(), 1, MPI_INT, rcd.data(), 1, MPI_INT, comm);
  MPI_Alltoall(sci.data(), 1, MPI_INT, rci.data(), 1, MPI_INT, comm);
  for (int r = 0; r < size; r++) rdd[r + 1] = rdd[r] + rcd[r], rdi[r + 1] = rdi[r] + rci[r];
  vector<double> rd(rdd[size]);
  vector<int> ri(rdi[size]);
  MPI_Alltoallv(bd.data(), scd.data(), sdd.data(), MPI_DOUBLE, rd.data(), rcd.data(), rdd.data(), MPI_DOUBLE, comm);
  MPI_Alltoallv(bi.data(), sci.data(), sdi.data(), MPI_INT, ri.data(), rci.data(), rdi.data(), MPI_INT, comm);

  // frozen part
  F = DPiece();
  vector<int> fn(P.nv(), -1);
  for (int k = 0; k < P.nt(); k++)
    if (tgt[k] < 0) {
      for (int j = 0; j < 3; j++) {
        int v = P.tri[3 * k + j];
        if (fn[v] < 0) fn[v] = F.addv(P.xy[2 * v], P.xy[2 * v + 1], P.vlab[v], &P.met[3 * v]);
        F.tri.push_back(fn[v]);
      }
      F.tlab.push_back(P.tlab[k]);
    }
  for (int i = 0; i < (int)P.blab.size(); i++) {
    const auto &p = et[Ek(P.bed[2 * i], P.bed[2 * i + 1])];
    bool fz = (p.first >= 0 && tgt[p.first] < 0) || (p.second >= 0 && tgt[p.second] < 0);
    if (fz) F.bed.push_back(fn[P.bed[2 * i]]), F.bed.push_back(fn[P.bed[2 * i + 1]]), F.blab.push_back(P.blab[i]);
  }

  // received piece
  B = DPiece();
  PointHash ph(tol);
  std::map<EKey, std::pair<int, int> > einfo;    // edge -> (user label, flag), counter
  std::map<EKey, int> ecnt;
  int ntr = ri.size() / 10;
  for (int q = 0; q < ntr; q++) {
    const double *d = &rd[15 * q];
    const int *it = &ri[10 * q];
    int vv[3];
    for (int j = 0; j < 3; j++) {
      R2 Q(d[5 * j], d[5 * j + 1]);
      int g = ph.find(Q);
      if (g < 0) {
        g = ph.add(Q);
        B.addv(Q.x, Q.y, it[j], d + 5 * j + 2);
      }
      vv[j] = g;
    }
    for (int j = 0; j < 3; j++) B.tri.push_back(vv[j]);
    B.tlab.push_back(it[3]);
    for (int e = 0; e < 3; e++) {
      EKey key = Ek(vv[(e + 1) % 3], vv[(e + 2) % 3]);
      int lab = it[4 + 2 * e], flag = it[5 + 2 * e];
      if (ecnt[key]++ == 0)
        einfo[key] = make_pair(lab, flag);
      else if (lab != INT_MIN)
        einfo[key].first = lab;
    }
  }
  for (auto &e : ecnt) {
    const auto &inf = einfo[e.first];
    int a = e.first.first, b = e.first.second;
    if (inf.first != INT_MIN) B.bed.push_back(a), B.bed.push_back(b), B.blab.push_back(inf.first);
    if (e.second == 1 && inf.first == INT_MIN && inf.second != 1) {    // cut: seeds
      seedB.insert(DXY(B.xy[2 * a], B.xy[2 * a + 1]));
      seedB.insert(DXY(B.xy[2 * b], B.xy[2 * b + 1]));
    }
  }
}

// Gluing of two pieces of the same process (common vertices by coordinates).
static DPiece MergePieces(const DPiece &A, const DPiece &B, double tol) {
  DPiece P = A;
  PointHash ph(tol);
  for (int i = 0; i < A.nv(); i++) ph.add(R2(A.xy[2 * i], A.xy[2 * i + 1]));
  vector<int> nb(B.nv());
  for (int i = 0; i < B.nv(); i++) {
    R2 Q(B.xy[2 * i], B.xy[2 * i + 1]);
    int g = ph.find(Q);
    if (g < 0) g = ph.add(Q), P.addv(Q.x, Q.y, B.vlab[i], &B.met[3 * i]);
    nb[i] = g;
  }
  for (int k = 0; k < B.nt(); k++) {
    for (int j = 0; j < 3; j++) P.tri.push_back(nb[B.tri[3 * k + j]]);
    P.tlab.push_back(B.tlab[k]);
  }
  std::set<EKey> seen;
  for (int i = 0; i < (int)A.blab.size(); i++) seen.insert(Ek(A.bed[2 * i], A.bed[2 * i + 1]));
  for (int i = 0; i < (int)B.blab.size(); i++) {
    int a = nb[B.bed[2 * i]], b = nb[B.bed[2 * i + 1]];
    if (seen.insert(Ek(a, b)).second) P.bed.push_back(a), P.bed.push_back(b), P.blab.push_back(B.blab[i]);
  }
  return P;
}

// Connected components (through the vertices) of the triangles act, joined between
// processes by the common boundary vertices, given entirely (the largest
// first, to the least loaded process). tgt = -1 outside act.
static void DistributedComponents(const DPiece &P, const vector<char> &act, MPI_Comm comm, int size, int rank,
                                  double tol, vector<int> &tgt) {
  int nt = P.nt();
  vector<vector<int> > vt(P.nv());
  for (int k = 0; k < nt; k++)
    if (act[k])
      for (int j = 0; j < 3; j++) vt[P.tri[3 * k + j]].push_back(k);
  vector<int> comp(nt, -1);
  int nc = 0;
  vector<long> csize;
  for (int k0 = 0; k0 < nt; k0++)
    if (act[k0] && comp[k0] < 0) {
      vector<int> st(1, k0);
      comp[k0] = nc;
      long n = 0;
      while (!st.empty()) {
        int k = st.back();
        st.pop_back();
        n++;
        for (int j = 0; j < 3; j++)
          for (int kk : vt[P.tri[3 * k + j]])
            if (comp[kk] < 0) comp[kk] = nc, st.push_back(kk);
      }
      csize.push_back(n);
      nc++;
    }
  // boundary vertices of the piece touching a component
  std::map<EKey, int> cnt;
  for (int k = 0; k < nt; k++)
    for (int e = 0; e < 3; e++) cnt[Ek(P.tri[3 * k + (e + 1) % 3], P.tri[3 * k + (e + 2) % 3])]++;
  vector<char> bv(P.nv(), 0);
  for (auto &c : cnt)
    if (c.second == 1) bv[c.first.first] = bv[c.first.second] = 1;
  vector<double> pts;
  vector<int> pc;
  for (int v = 0; v < P.nv(); v++)
    if (bv[v] && !vt[v].empty()) pts.push_back(P.xy[2 * v]), pts.push_back(P.xy[2 * v + 1]), pc.push_back(comp[vt[v][0]]);
  vector<int> ncs(size);
  MPI_Allgather(&nc, 1, MPI_INT, ncs.data(), 1, MPI_INT, comm);
  vector<int> off(size + 1, 0);
  for (int r = 0; r < size; r++) off[r + 1] = off[r] + ncs[r];
  for (int &c : pc) c += off[rank];
  vector<double> apts;
  vector<int> apc, cnt1, dsp1, cnt2, dsp2;
  AllGather(comm, size, pts, apts, cnt1, dsp1, MPI_DOUBLE);
  AllGather(comm, size, pc, apc, cnt2, dsp2, MPI_INT);
  vector<long> acs;
  {
    vector<long> lcs(csize), tmp;
    vector<int> c3, d3;
    AllGather(comm, size, lcs, acs, c3, d3, MPI_LONG);
  }
  int NC = off[size];
  vector<int> uf(NC);
  for (int i = 0; i < NC; i++) uf[i] = i;
  std::function<int(int)> root = [&](int i) { return uf[i] == i ? i : uf[i] = root(uf[i]); };
  PointHash ph(tol);
  vector<int> pcomp;
  for (int i = 0; i < (int)apc.size(); i++) {
    R2 Q(apts[2 * i], apts[2 * i + 1]);
    int g = ph.find(Q);
    if (g < 0)
      ph.add(Q), pcomp.push_back(apc[i]);
    else {
      int a = root(pcomp[g]), b = root(apc[i]);
      if (a != b) uf[std::max(a, b)] = std::min(a, b);
    }
  }
  vector<long> rs(NC, 0);
  for (int i = 0; i < NC; i++) rs[root(i)] += acs[i];
  vector<int> order;
  for (int i = 0; i < NC; i++)
    if (root(i) == i) order.push_back(i);
  std::stable_sort(order.begin(), order.end(), [&rs](int a, int b) { return rs[a] > rs[b]; });
  vector<long> load(size, 0);
  vector<int> assign(NC, 0);
  for (int c : order) {
    int r = std::min_element(load.begin(), load.end()) - load.begin();
    load[r] += rs[c];
    assign[c] = r;
  }
  tgt.assign(nt, -1);
  for (int k = 0; k < nt; k++)
    if (act[k]) tgt[k] = assign[root(comp[k] + off[rank])];
}

// Band of nlayer layers of triangles (through the vertices) around the boundary edges of P whose label is in labs
// (the periodic boundaries). The band is the same as on the global mesh: after each layer, the active vertices are
// shared between the processes by their coordinates (AllGather + PointHash). Without this exchange, a triangle
// touching a periodic boundary only by a vertex, on a process which has no periodic edge there, is left out of the
// band: the band is pinched at this vertex, its inner (frozen) boundary meets the periodic boundary on one side
// only, the 2 sides are cut into different curves and the periodic remeshing of bamg crashes (GeomToTriangles1).
static void BandAround(const DPiece &P, const vector<long> &labs, long nlayer, vector<char> &act, MPI_Comm comm,
                       int mpisize, double tol) {
  int nt = P.nt(), nv = P.nv();
  act.assign(nt, 0);
  vector<char> vact(nv, 0);
  auto share = [&]() {    // vact <- vact on all the processes (vertices with the same coordinates)
    if (mpisize <= 1) return;
    vector<double> loc, all;
    vector<int> cnt, dsp;
    for (int i = 0; i < nv; i++)
      if (vact[i]) loc.push_back(P.xy[2 * i]), loc.push_back(P.xy[2 * i + 1]);
    AllGather(comm, mpisize, loc, all, cnt, dsp, MPI_DOUBLE);
    PointHash ph(tol);
    for (size_t k = 0; k + 1 < all.size(); k += 2) ph.add(R2(all[k], all[k + 1]));
    for (int i = 0; i < nv; i++)
      if (!vact[i] && ph.find(R2(P.xy[2 * i], P.xy[2 * i + 1])) >= 0) vact[i] = 1;
  };
  for (int i = 0; i < (int)P.blab.size(); i++)
    if (std::find(labs.begin(), labs.end(), (long)P.blab[i]) != labs.end()) vact[P.bed[2 * i]] = vact[P.bed[2 * i + 1]] = 1;
  share();
  for (int l = 0; l <= nlayer; l++) {
    for (int k = 0; k < nt; k++)
      if (vact[P.tri[3 * k]] || vact[P.tri[3 * k + 1]] || vact[P.tri[3 * k + 2]]) act[k] = 1;
    for (int k = 0; k < nt; k++)
      if (act[k]) vact[P.tri[3 * k]] = vact[P.tri[3 * k + 1]] = vact[P.tri[3 * k + 2]] = 1;
    if (l < nlayer) share();
  }
}

// Distributed recursive bisection of the barycentres: tgt = process of each
// triangle. The medians are found by bisection on counts
// (MPI_Allreduce), without gathering the triangles.
// fixed0 (periodic boundaries): these triangles go to process 0, the others are cut between the processes 1 to size-1.
static void DistributedRCB(const DPiece &P, MPI_Comm comm, int size, vector<int> &tgt, const vector<char> *fixed0 = 0) {
  int nt = P.nt();
  vector<double> g(2 * nt);
  for (int k = 0; k < nt; k++)
    for (int c = 0; c < 2; c++)
      g[2 * k + c] = (P.xy[2 * P.tri[3 * k] + c] + P.xy[2 * P.tri[3 * k + 1] + c] + P.xy[2 * P.tri[3 * k + 2] + c]) / 3;
  const int p0 = (fixed0 && size > 1) ? 1 : 0;
  vector<int> lo(nt, p0), np(nt, size - p0);
  if (p0)
    for (int k = 0; k < nt; k++)
      if ((*fixed0)[k]) lo[k] = 0, np[k] = 1;
  for (int level = 0; (1 << level) < size - p0; level++) {
    // group = first process lo; bounding box and number of triangles
    vector<double> bmin(2 * size, 1e300), bmax(2 * size, -1e300), gmin(2 * size), gmax(2 * size);
    vector<long> n(size, 0), gn(size);
    for (int k = 0; k < nt; k++)
      if (np[k] > 1) {
        n[lo[k]]++;
        for (int c = 0; c < 2; c++)
          bmin[2 * lo[k] + c] = std::min(bmin[2 * lo[k] + c], g[2 * k + c]),
                           bmax[2 * lo[k] + c] = std::max(bmax[2 * lo[k] + c], g[2 * k + c]);
      }
    MPI_Allreduce(bmin.data(), gmin.data(), 2 * size, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(bmax.data(), gmax.data(), 2 * size, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(n.data(), gn.data(), size, MPI_LONG, MPI_SUM, comm);
    vector<int> axis(size, 0), gnp(size, 0);
    for (int k = 0; k < nt; k++) gnp[lo[k]] = np[k];
    {
      vector<int> t(gnp);
      MPI_Allreduce(t.data(), gnp.data(), size, MPI_INT, MPI_MAX, comm);
    }
    vector<double> a(size), b(size), cut(size);
    vector<long> want(size, 0);
    for (int r = 0; r < size; r++)
      if (gn[r] > 0 && gnp[r] > 1) {
        axis[r] = (gmax[2 * r] - gmin[2 * r] >= gmax[2 * r + 1] - gmin[2 * r + 1]) ? 0 : 1;
        a[r] = gmin[2 * r + axis[r]], b[r] = gmax[2 * r + axis[r]];
        want[r] = (long)((double)gn[r] * (gnp[r] / 2) / gnp[r]);
      }
    for (int it = 0; it < 60; it++) {    // bisection: number of barycentres < cut
      for (int r = 0; r < size; r++) cut[r] = 0.5 * (a[r] + b[r]);
      vector<long> c(size, 0), gc(size);
      for (int k = 0; k < nt; k++)
        if (np[k] > 1 && g[2 * k + axis[lo[k]]] < cut[lo[k]]) c[lo[k]]++;
      MPI_Allreduce(c.data(), gc.data(), size, MPI_LONG, MPI_SUM, comm);
      for (int r = 0; r < size; r++)
        if (gc[r] < want[r])
          a[r] = cut[r];
        else
          b[r] = cut[r];
    }
    for (int k = 0; k < nt; k++)
      if (np[k] > 1) {
        int h = np[k] / 2;
        if (g[2 * k + axis[lo[k]]] < b[lo[k]])
          np[k] = h;
        else
          lo[k] += h, np[k] -= h;
      }
  }
  tgt = lo;
}

AnyType AdaptationDDM::Distributed(Stack stack) const {
  using namespace bamg;
  using bamg::Abs;
  using bamg::Max;
  using bamg::Min;

  Real8 err = arg(2, stack, 0.01);
  Real8 errg = arg(3, stack, Min(0.01, err));
  long nbsx = Max(100L, arg(4, stack, 9000L));
  long nbsmooth = arg(5, stack, 3L);
  long nbjacobi = arg(6, stack, 0L);
  const Real8 raison = arg(7, stack, 1.8);
  const Real8 omega = arg(8, stack, 1.0);
  bool iso = arg(9, stack, false);
  bool AbsError = arg(10, stack, true);
  Real8 CutOff = arg(11, stack, 1.0e-6);
  verbosity = arg(12, stack, (long)verbosity);
  bool SplitEdgeWith2Boundary = arg(14, stack, true);
  double maxsubdiv = Max(Min(arg(15, stack, 10.0), 10.0), 0.1);
  double anisomax = Max((double)arg(16, stack, 1.0e6), 1.0);
  bool rescaling = arg(17, stack, true);
  bool KeepBackVertices = arg(18, stack, true);
  int givenmetric = arg(19, stack, false);
  double powerM = arg(20, stack, 1.0);
  double cutoffradian = arg(21, stack, -1.0) * bamg::Pi / 180.;
  bool split = arg(22, stack, false);
  bool nomeshgeneration = arg(23, stack, false);
  long lwarning = 0;
  long *pwarning = arg(27, stack, &lwarning);
  long &warning = *pwarning;
  long nbpass = Max(1L, arg(29, stack, 3L));
  long nblayer = Max(0L, arg(30, stack, 6L));

  pcommworld pcomm = nargs[28] ? GetAny<pcommworld>((*nargs[28])(stack)) : nullptr;
  int mpiinit = 0;
  MPI_Initialized(&mpiinit);
  if (!mpiinit) ExecError("adaptmeshddm : owned= requires MPI (FreeFem++-mpi)");
  MPI_Comm comm = pcomm ? *(MPI_Comm *)pcomm : MPI_COMM_WORLD;
  int mpirank = 0, mpisize = 1;
  MPI_Comm_rank(comm, &mpirank);
  MPI_Comm_size(comm, &mpisize);

  KN<long> reqedges0;
  KN<long> reqedges(nargs[26] ? GetAny<KN_<long> >((*nargs[26])(stack)) : (KN_<long>)reqedges0);
  const Mesh *Thh = GetAny<pmesh>((*getmesh)(stack));
  ffassert(Thh);
  const Mesh &Th(*Thh);
  KN<double> *pown = GetAny<KN<double> *>((*nargs[31])(stack));
  if (!pown || pown->N() != Th.nt) ExecError("adaptmeshddm : owned= must have one element per triangle (P0)");
  vector<char> own(Th.nt);
  for (int k = 0; k < Th.nt; k++) own[k] = fabs((*pown)[k]) > 1e-2;

  const int fake = -111111;    // label of the artificial boundaries of macro_ddm
  // label of the frozen edges (largest label + 1, on all the processes)
  long labi = 0;
  for (int i = 0; i < Th.neb; i++) labi = max(labi, (long)Th.bedges[i].lab);
  for (int i = 0; i < reqedges.N(); i++) labi = max(labi, reqedges[i]);
  {
    long g = labi;
    MPI_Allreduce(&labi, &g, 1, MPI_LONG, MPI_MAX, comm);
    labi = g + 1;
  }
  KN<long> reqlab(reqedges.N() + 1);
  for (int i = 0; i < reqedges.N(); i++) reqlab[i] = reqedges[i];
  reqlab[reqedges.N()] = labi;
  // periodic=: on one process, as adaptmesh. On several, the periodic boundaries are frozen during the passes
  // (reqlabP), then a band around them is gathered on process 0 and adapted again with the periodicity, and it
  // stays on process 0 at the rebalancing: the 2 sides of a periodic boundary are in the same local mesh (as
  // macro ThPeriodicity of macro_ddm.idp, needed by fvscheme(..., periodic=...)); the other processes have no
  // periodic boundary
  vector<long> plabs;
  for (int i = 0; i < nbcperiodic; i++)
    plabs.push_back(GetAny<long>((*periodic[4 * i])(stack))), plabs.push_back(GetAny<long>((*periodic[4 * i + 2])(stack)));
  const bool perextra = nbcperiodic && mpisize > 1;
  KN<long> reqlabP(reqlab.N() + plabs.size());    // + frozen periodic labels
  for (int i = 0; i < reqlab.N(); i++) reqlabP[i] = reqlab[i];
  for (int i = 0; i < (int)plabs.size(); i++) reqlabP[reqlab.N() + i] = plabs[i];
  const KN<long> &rqP = perextra ? reqlabP : reqlab;    // required labels of the passes
  // labels of the internal edges of the user (on all the processes)
  std::set<long> intlabs;
  {
    vector<int> li, all, c, d;
    for (int i = 0; i < Th.neb; i++) {
      int e, k = Th.BoundaryElement(i, e);
      int ee = e, kk = Th.ElementAdj(k, ee);
      if (kk >= 0 && kk != k && Th.bedges[i].lab != fake) li.push_back(Th.bedges[i].lab);
    }
    AllGather(comm, mpisize, li, all, c, d, MPI_INT);
    intlabs.insert(all.begin(), all.end());
  }
  // global bounding box: default hmin, hmax and tolerance of the gluing
  Fem2D::R2 Pmin, Pmax;
  Th.BoundingBox(Pmin, Pmax);
  double bb[4] = {-Pmin.x, -Pmin.y, Pmax.x, Pmax.y}, gbb[4];
  MPI_Allreduce(bb, gbb, 4, MPI_DOUBLE, MPI_MAX, comm);
  Real8 extent = 1.1 * Max(gbb[2] + gbb[0], gbb[3] + gbb[1]);
  Real8 hmax = 0.3 * extent;
  Real8 hmin = 2.0 * extent / MaxICoor;
  hmin = Max(hmin, arg(0, stack, hmin));
  hmax = Min(hmax, arg(1, stack, hmax));
  if (iso) anisomax = 1;
  const double tol = 1e-10 * extent;
  int ksol = 0;
  for (int i = 0; i < nbsol; i++) ksol += typesol[i] + 1;
  GenParams gp = {KeepBackVertices, SplitEdgeWith2Boundary, nbsmooth, omega};

  bool mtx = em11 && em22 && em12;
  KN<double> *mm11 = 0, *mm12 = 0, *mm22 = 0;
  if (mtx) {
    mm11 = GetAny<KN<double> *>((*em11)(stack));
    mm22 = GetAny<KN<double> *>((*em22)(stack));
    mm12 = GetAny<KN<double> *>((*em12)(stack));
    if (mm11->N() != Th.nv || mm22->N() != Th.nv || mm12->N() != Th.nv)
      ExecError("The size of3  metrics array must be equal to nb of vertex");
  }

  MeshPoint *mp(Fem2D::MeshPointStack(stack)), mps = *mp;
  const long verb0 = verbosity;
  long lwarn = 0;
  int lerr = 0;
  string msg;
  double t0 = MPI_Wtime();

  // ---- pass 1: own triangles, metric on the whole local mesh
  DPiece P;
  vector<char> seeds;
  vector<int> l2g, gl2g;
  Mesh *loc = 0, *gloc = 0;
  Triangles *oTh = 0, *gTh = 0;
  double *lessol = 0;
  vector<double> smin(nbsol, 1e300), smax(nbsol, -1e300), sabs(nbsol, 0.);
  vector<double> sming(nbsol, 1e300), smaxg(nbsol, -1e300), sabsg(nbsol, 0.);
  std::map<DSeg, int> keep;    // internal edges of the user on an interface
  int ntl = 0;
  for (int k = 0; k < Th.nt; k++) ntl += own[k];
  try {
    if (mpirank) verbosity = 0;
    // whole local mesh (overlap included): artificial boundaries -> labi
    {
      vector<char> all(Th.nt, 1);
      gloc = SubMesh(Th, all, labi, gl2g);
      for (int i = 0; i < gloc->neb; i++)
        if (gloc->bedges[i].lab == fake) gloc->bedges[i].lab = labi;
      gTh = msh2bamgDDM(*gloc, cutoffradian, (long *)reqlab, reqlab.N(), labi, &Th, &gl2g);
      gTh->ReMakeTriangleContainingTheVertex();
    }
    // own triangles: edges towards the other processes -> labi
    if (ntl) {
      loc = SubMesh(Th, own, labi, l2g);
      for (int k = 0; k < Th.nt; k++)
        if (own[k])
          for (int e = 0; e < 3; e++) {
            int ee = e, kk = Th.ElementAdj(k, ee);
            if (kk < 0 || kk == k || own[kk]) continue;
            int i0 = Th(k, (e + 1) % 3), i1 = Th(k, (e + 2) % 3);
            Fem2D::BoundaryEdge *pbe = Th.TheBoundaryEdge(i0, i1);
            if (pbe && pbe->lab != fake) keep[Sg(Th(i0).x, Th(i0).y, Th(i1).x, Th(i1).y)] = pbe->lab;
          }
      for (int i = 0; i < loc->neb; i++)
        if (loc->bedges[i].lab == fake) loc->bedges[i].lab = labi;
      if (nbcperiodic && !perextra) {    // one process: the equivalent periodic edges, as adaptmesh
        KN<int> ndfv(loc->nv), ndfe(loc->neb);
        for (int i = 0; i < ndfv.N(); ++i) ndfv[i] = i;
        int nbdfv = 0, nbdfe = 0;
        BuildPeriodic2(nbcperiodic, periodic, *loc, stack, nbdfv, ndfv, nbdfe, ndfe);
        *mp = mps;
        oTh = msh2bamgDDM(*loc, cutoffradian, (long *)reqlab, reqlab.N(), labi, &Th, &l2g, nbdfe, (int *)ndfv, (int *)ndfe);
      } else
        oTh = msh2bamgDDM(*loc, cutoffradian, (long *)rqP, rqP.N(), labi, &Th, &l2g);
      oTh->ReMakeTriangleContainingTheVertex();
    }
    // solutions at the vertices of the local mesh (numbers of Th = numbers of gTh)
    int nvg = gl2g.size();
    lessol = new double[nvg * ksol];
    vector<char> done(nvg, 0), vown(Th.nv, 0);
    for (int k = 0; k < Th.nt; k++)
      if (own[k]) vown[Th(k, 0)] = vown[Th(k, 1)] = vown[Th(k, 2)] = 1;
    for (int k = 0; k < Th.nt; k++)
      for (int jt = 0; jt < 3; jt++) {
        int iv = Th(k, jt);
        if (!done[iv]) {
          done[iv] = 1;
          mp->setP(&Th, k, jt);
          double *ss = lessol + ksol * iv;
          for (int j = 0; j < ksol; j++) *ss++ = GetAny<double>((*sol[j])(stack));
        }
      }
    *mp = mps;
    for (int i = 0, k0 = 0; i < nbsol; k0 += typesol[i] + 1, i++)
      for (int iv = 0; iv < nvg; iv++) {
        double *s = lessol + ksol * iv + k0, v = s[0];
        if (typesol[i]) {
          v = 0;
          for (int j = 0; j <= typesol[i]; j++) v += s[j] * s[j];
          v = sqrt(v);
        }
        sming[i] = Min(sming[i], v), smaxg[i] = Max(smaxg[i], v), sabsg[i] = Max(sabsg[i], Abs(v));
        if (vown[iv]) smin[i] = Min(smin[i], v), smax[i] = Max(smax[i], v), sabs[i] = Max(sabs[i], Abs(v));
      }
  } catch (std::exception &e) {
    lerr = 1, msg = e.what();
  } catch (...) {
    lerr = 1, msg = "unknown error";
  }
  verbosity = verb0;
  *mp = mps;
  CheckErrors(comm, mpisize, lerr, msg);
  vector<double> gmin(smin), gmax(smax), gabs(sabs);
  if (nbsol) {
    MPI_Allreduce(smin.data(), gmin.data(), nbsol, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(smax.data(), gmax.data(), nbsol, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(sabs.data(), gabs.data(), nbsol, MPI_DOUBLE, MPI_MAX, comm);
  }

  // metric: as adaptmesh, on the whole local mesh
  double lmass = 0;
  try {
    if (mpirank) verbosity = 0;
    Triangles &GTh(*gTh);
    int nvg = gl2g.size();
    Metric Mhmax(hmax);
    for (int iv = 0; iv < GTh.nbv; iv++) GTh[iv].m = Mhmax;
    if (mtx)
      for (int iv = 0; iv < GTh.nbv; iv++) {
        int g = gl2g[iv];
        GTh[iv].m.IntersectWith(MetricAnIso((*mm11)[g], (*mm12)[g], (*mm22)[g]));
      }
    if (givenmetric) {
      if (ksol == 1) {
        for (int iv = 0; iv < GTh.nbv; iv++) GTh[iv].m.IntersectWith(Metric(lessol[iv]));
      } else if (ksol == 3) {
        for (int iv = 0, k = 0; iv < GTh.nbv; iv++, k += 3) {
          Metric MM(lessol[k], lessol[k + 1], lessol[k + 2]);
          MatVVP2x2 vp(MM);
          vp.Abs();
          GTh[iv].m.IntersectWith(vp);
        }
      } else
        ExecError("adaptmeshddm : IsMetric and ksol != 1 or 3");
    } else
      for (int i = 0, k0 = 0; i < nbsol; k0 += typesol[i] + 1, i++) {
        int nf = typesol[i] + 1;
        double dl = smaxg[i] - sming[i], dg = gmax[i] - gmin[i];
        double coef = sqrt(err), cut = 0;
        if (AbsError) {
          if (dl <= 1e-10 * Max(sabsg[i], 1e-20)) continue;
          if (rescaling) coef = sqrt(err * dg / dl);
        } else {
          if (sabsg[i] <= 0) continue;
          cut = CutOff * gabs[i] / sabsg[i];
        }
        vector<double> s1(nvg * nf);
        for (int iv = 0; iv < nvg; iv++)
          for (int j = 0; j < nf; j++) s1[iv * nf + j] = lessol[iv * ksol + k0 + j];
        GTh.IntersectConsMetric(s1.data(), 1, typesol + i, hmin, hmax, coef, anisomax, cut, nbjacobi, rescaling,
                                powerM, 0);
      }
    GTh.IntersectGeomMetric(errg, iso);
    GTh.SmoothMetric(raison);
    GTh.MaxSubDivision(maxsubdiv);
    GTh.BoundAnisotropy(anisomax);
    if (mtx)    // metric= on output: at the vertices of the local mesh
      for (int iv = 0; iv < GTh.nbv; iv++) {
        int g = gl2g[iv];
        (*mm11)[g] = GTh[iv].m.a11, (*mm12)[g] = GTh[iv].m.a21, (*mm22)[g] = GTh[iv].m.a22;
      }
    if (ntl && !nomeshgeneration) {
      Triangles &BTh(*oTh);
      for (int iv = 0; iv < BTh.nbv; iv++) BTh[iv].m = GTh[l2g[iv]].m;    // gl2g = identity
      for (int k = 0; k < loc->nt; k++) {
        double d = 0;
        for (int j = 0; j < 3; j++) {
          const MetricAnIso &m = BTh[(*loc)(k, j)].m;
          d += sqrt(Max(m.a11 * m.a22 - m.a21 * m.a21, 0.)) / 3;
        }
        lmass += (*loc)[k].area * d;
      }
    }
  } catch (std::exception &e) {
    lerr = 1, msg = e.what();
  } catch (...) {
    lerr = 1, msg = "unknown error";
  }
  verbosity = verb0;
  delete[] lessol;
  delete gTh;
  delete gloc;
  CheckErrors(comm, mpisize, lerr, msg);
  if (nomeshgeneration) {
    delete oTh;
    delete loc;
    return SetAny<pmesh>(Thh);
  }

  long cap = ShareBudget(comm, mpisize, mpirank, lmass, ntl ? l2g.size() + 10. : 0., nbsx);
  try {
    if (mpirank) verbosity = 0;
    if (ntl) GenerateDPiece(loc, oTh, cap, gp, labi, keep, P, seeds, lwarn);
  } catch (std::exception &e) {
    lerr = 1, msg = e.what();
  } catch (...) {
    lerr = 1, msg = "unknown error";
  }
  verbosity = verb0;
  delete oTh;
  delete loc;
  CheckErrors(comm, mpisize, lerr, msg);
  if (verbosity > 1) {
    long n = P.nt(), gn;
    MPI_Allreduce(&n, &gn, 1, MPI_LONG, MPI_SUM, comm);
    if (mpirank == 0)
      cout << "  -- adaptmeshddm (distributed) pass 0 : " << gn << " triangles, " << MPI_Wtime() - t0 << " s" << endl;
  }

  // remeshing of a received piece B (metric carried with the vertices) -> Bn; per: the periodic boundaries of B
  // adapted with their equivalent edges (BuildPeriodic2), otherwise frozen (rqP)
  auto RemeshPiece = [&](const DPiece &B, long capb, bool per, DPiece &Bn) {
    vector<char> onl;
    try {
      if (mpirank) verbosity = 0;
      if (B.nt()) {
        std::map<DSeg, int> keepB;
        Mesh *bl = PieceMesh(B, labi, intlabs, &keepB);
        Triangles *bTh = 0;
        if (per) {
          KN<int> ndfv(bl->nv), ndfe(bl->neb);
          for (int i = 0; i < ndfv.N(); ++i) ndfv[i] = i;
          int nbdfv = 0, nbdfe = 0;
          BuildPeriodic2(nbcperiodic, periodic, *bl, stack, nbdfv, ndfv, nbdfe, ndfe);
          *mp = mps;
          bTh = msh2bamgDDM(*bl, cutoffradian, (long *)reqlab, reqlab.N(), labi, 0, 0, nbdfe, (int *)ndfv, (int *)ndfe);
        } else
          bTh = msh2bamgDDM(*bl, cutoffradian, (long *)rqP, rqP.N(), labi, 0, 0);
        bTh->ReMakeTriangleContainingTheVertex();
        for (int iv = 0; iv < bTh->nbv; iv++)
          (*bTh)[iv].m = MetricAnIso(B.met[3 * iv], B.met[3 * iv + 1], B.met[3 * iv + 2]);
        bTh->SmoothMetric(raison);    // gradation across the old interfaces
        GenerateDPiece(bl, bTh, capb, gp, labi, keepB, Bn, onl, lwarn);
        delete bTh;
        delete bl;
      }
    } catch (std::exception &e) {
      lerr = 1, msg = e.what();
    } catch (...) {
      lerr = 1, msg = "unknown error";
    }
    verbosity = verb0;
    *mp = mps;
    CheckErrors(comm, mpisize, lerr, msg);
  };

  // ---- passes 2 and 3: bands around the interfaces, sent and remeshed
  for (int pass = 1; pass < nbpass; ++pass) {
    t0 = MPI_Wtime();
    int nt = P.nt();
    vector<char> act(nt, 0), vact(seeds);
    vact.resize(P.nv(), 0);
    for (int l = 0; l <= nblayer; l++) {
      for (int k = 0; k < nt; k++)
        if (vact[P.tri[3 * k]] || vact[P.tri[3 * k + 1]] || vact[P.tri[3 * k + 2]]) act[k] = 1;
      for (int k = 0; k < nt; k++)
        if (act[k]) vact[P.tri[3 * k]] = vact[P.tri[3 * k + 1]] = vact[P.tri[3 * k + 2]] = 1;
    }
    long nact = 0, gnact;
    for (int k = 0; k < nt; k++) nact += act[k];
    MPI_Allreduce(&nact, &gnact, 1, MPI_LONG, MPI_SUM, comm);
    if (gnact == 0) break;
    vector<int> tgt(nt, -1);
    if (pass == 1) {    // bisection rotated by pi/4 of the barycentres of the band
      double c = cos(bamg::Pi / 4), s = sin(bamg::Pi / 4);
      vector<double> uv, auv;
      for (int k = 0; k < nt; k++)
        if (act[k]) {
          double gx = 0, gy = 0;
          for (int j = 0; j < 3; j++) gx += P.xy[2 * P.tri[3 * k + j]] / 3, gy += P.xy[2 * P.tri[3 * k + j] + 1] / 3;
          uv.push_back(c * gx + s * gy), uv.push_back(-s * gx + c * gy);
        }
      vector<int> cn, ds;
      AllGather(comm, mpisize, uv, auv, cn, ds, MPI_DOUBLE);
      int N = auv.size() / 2;
      vector<double> u(N), v(N);
      vector<int> idx(N), part(N);
      for (int i = 0; i < N; i++) u[i] = auv[2 * i], v[i] = auv[2 * i + 1], idx[i] = i;
      RCB(u, v, idx.data(), N, mpisize, 0, part);
      for (int k = 0, q = ds[mpirank] / 2; k < nt; k++)
        if (act[k]) tgt[k] = part[q++];
    } else
      DistributedComponents(P, act, comm, mpisize, mpirank, tol, tgt);

    DPiece F, B, Bn;
    std::set<DXY> seedB;
    Migrate(P, tgt, comm, mpisize, tol, F, B, seedB);

    // budget: nbvx minus the vertices of the frozen part (all processes)
    long nf = F.nv(), gnf;
    MPI_Allreduce(&nf, &gnf, 1, MPI_LONG, MPI_SUM, comm);
    long capb = ShareBudget(comm, mpisize, mpirank, PieceMass(B), B.nt() ? B.nv() + 10. : 0., nbsx - gnf);
    RemeshPiece(B, capb, false, Bn);
    P = MergePieces(F, Bn, tol);
    seeds.assign(P.nv(), 0);
    for (int i = 0; i < P.nv(); i++) seeds[i] = seedB.count(DXY(P.xy[2 * i], P.xy[2 * i + 1])) > 0;
    if (verbosity > 1) {
      long n[2] = {P.nt(), nact}, gn[2];
      MPI_Allreduce(n, gn, 2, MPI_LONG, MPI_SUM, comm);
      if (mpirank == 0)
        cout << "  -- adaptmeshddm (distributed) pass " << pass << " : " << gn[1] << " triangles adapted again -> " << gn[0]
             << " triangles, " << MPI_Wtime() - t0 << " s" << endl;
    }
  }

  // ---- pass of the periodic boundaries (several processes): a band of nblayer layers around them gathered on
  // process 0 and adapted again with the equivalent periodic edges (the 2 sides discretized in the same way)
  const long nbandp = Max(nblayer, 3L);    // >= the overlap + 1: the overlap of the other processes stays off them
  if (perextra) {
    t0 = MPI_Wtime();
    vector<char> act;
    BandAround(P, plabs, nbandp, act, comm, mpisize, tol);
    long nact = 0, gnact;
    for (int k = 0; k < P.nt(); k++) nact += act[k];
    MPI_Allreduce(&nact, &gnact, 1, MPI_LONG, MPI_SUM, comm);
    if (gnact) {
      vector<int> tgt(P.nt(), -1);
      for (int k = 0; k < P.nt(); k++)
        if (act[k]) tgt[k] = 0;
      DPiece F, B, Bn;
      std::set<DXY> seedB;
      Migrate(P, tgt, comm, mpisize, tol, F, B, seedB);
      long nf = F.nv(), gnf;
      MPI_Allreduce(&nf, &gnf, 1, MPI_LONG, MPI_SUM, comm);
      long capb = ShareBudget(comm, mpisize, mpirank, PieceMass(B), B.nt() ? B.nv() + 10. : 0., nbsx - gnf);
      RemeshPiece(B, capb, true, Bn);
      P = MergePieces(F, Bn, tol);
      if (verbosity > 1) {
        long n = P.nt(), gn;
        MPI_Allreduce(&n, &gn, 1, MPI_LONG, MPI_SUM, comm);
        if (mpirank == 0)
          cout << "  -- adaptmeshddm (distributed) periodic pass : " << gnact << " triangles adapted again -> " << gn
               << " triangles, " << MPI_Wtime() - t0 << " s" << endl;
      }
    }
  }

  // ---- final rebalancing (as parmmg): distributed bisection, sending; periodic: the band around the periodic
  // boundaries on process 0, the rest cut between the other processes
  {
    t0 = MPI_Wtime();
    vector<int> tgt;
    vector<char> band;
    if (perextra) BandAround(P, plabs, nbandp, band, comm, mpisize, tol);
    DistributedRCB(P, comm, mpisize, tgt, perextra ? &band : 0);
    DPiece F, B;
    std::set<DXY> sd;
    Migrate(P, tgt, comm, mpisize, tol, F, B, sd);
    P = B;
    if (verbosity > 1) {
      int n = P.nt(), nmin, nmax;
      MPI_Allreduce(&n, &nmin, 1, MPI_INT, MPI_MIN, comm);
      MPI_Allreduce(&n, &nmax, 1, MPI_INT, MPI_MAX, comm);
      if (mpirank == 0)
        cout << "  -- adaptmeshddm (distributed) rebalancing: " << nmin << " to " << nmax << " triangles per process, "
             << MPI_Wtime() - t0 << " s" << endl;
    }
  }

  // ---- local mesh without overlap (interfaces: label -111111)
  for (int i = 0; i < P.nv(); i++)
    if (P.vlab[i] == labi) P.vlab[i] = 0;
  {
    int vide = P.nt() == 0, gvide;
    MPI_Allreduce(&vide, &gvide, 1, MPI_INT, MPI_MAX, comm);
    if (gvide) ExecError("adaptmeshddm : a process has no triangle (mesh too small for mpisize)");
  }
  std::set<long> none;
  Mesh *res = PieceMesh(P, fake, none, 0);
  if (split) {
    KN<int> s(res->nt);
    s = 2;
    Mesh *r2 = new Mesh(*res, (int *)s, false, 0);
    delete res;
    res = r2;
  }
  res->MakeQuadTree();
  long gwarn = lwarn;
  MPI_Allreduce(&lwarn, &gwarn, 1, MPI_LONG, MPI_SUM, comm);
  warning = gwarn;
  Add2StackOfPtr2FreeRC(stack, res);
  return SetAny<pmesh>(res);
}

// ============================================================================================================
// remapP0ddm(Tho, uo[], Do[], Th, u[]): the conservative transfer of P0 (finite volume) unknowns from a distributed
// mesh to an other one (domain decomposition, MPI_COMM_WORLD), the parallel version of remapP0 (msh3.cpp)
//   Tho: the old local mesh (with its overlap), uo: P0 on Tho (also on the overlap: after the exchange), Do: P0 on
//   Tho, 1 on the cells owned by this process (the partition of unity A.D of macro_ddm for P0), 0 elsewhere;
//   Th: the new local mesh, u: P0 on Th (all its cells, also its overlap)
//   uo is reconstructed affine in each owned old cell (least squares on the cells which share a vertex with it,
//   exact for an affine function, limited: Barth-Jespersen, no new extremum); each process sends its owned old
//   cells (vertices, values, slopes) to the processes whose new local mesh meets them (bounding boxes,
//   MPI_Alltoallv), and u on each new cell is the exact mean of these reconstructions (the exact intersections of
//   the triangles): the mass on the owned new cells is conserved exactly (each old cell is owned by one process)
//   a vector FE function [P0,...,P0] of m components (uo[m*k+c]): the m components together, the same limiter for
//   the components of a cell (a constant affine relation between them stays constant: the lake at rest H + zb)
// ============================================================================================================
namespace RemapDDM {
using Fem2D::R2;
static inline double Det(const R2 &a, const R2 &b, const R2 &c) { return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x); }
// the part of the convex polygon P (n points, counterclockwise) on the left of the line a -> b
static int Clip(const R2 *P, int n, const R2 &a, const R2 &b, R2 *Q) {
  int m = 0;
  for (int i = 0; i < n; ++i) {
    const R2 &p = P[i], &q = P[(i + 1)%n];
    double dp = Det(a, b, p), dq = Det(a, b, q);
    if (dp >= 0) Q[m++] = p;
    if ((dp > 0 && dq < 0) || (dp < 0 && dq > 0)) Q[m++] = p + (dp/(dp - dq))*(q - p);
  }
  return m;
}
// the area and the centroid C of the intersection of two triangles (counterclockwise)
static double TriInter(const R2 *A, const R2 *B, R2 &C) {
  R2 P[12], Q[12];
  int n = 3;
  for (int i = 0; i < 3; ++i) P[i] = A[i];
  for (int j = 0; j < 3 && n > 2; ++j) {
    n = Clip(P, n, B[j], B[(j + 1)%3], Q);
    for (int i = 0; i < n; ++i) P[i] = Q[i];
  }
  double s = 0;
  R2 c(0., 0.);
  for (int i = 1; i + 1 < n; ++i) {
    double d = Det(P[0], P[i], P[i + 1]);
    s += d;
    c += (P[0] + P[i] + P[i + 1])*(d/3.);
  }
  C = s > 0 ? c/s : A[0];
  return n > 2 ? max(0., s/2) : 0.;
}
}// namespace RemapDDM

static double RemapP0DDM(const pmesh &pTho, KN_< double > const &uo, KN_< double > const &Do, const pmesh &pTh, KN< double > *const &pu) {
  using namespace RemapDDM;
  if (!pTho || !pTh) ExecError("remapP0ddm: a mesh is not defined");
  const Mesh &Tho = *pTho, &Th = *pTh;
  int mpirank = 0, mpisize = 1;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &mpirank);
  MPI_Comm_size(comm, &mpisize);
  const int nto = Tho.nt, nvo = Tho.nv;
  int m = nto ? uo.N( )/nto : 1;
  {// the number of components: the same on all the processes (a process without old cell)
    int mg = 0, ml = nto ? m : 0;
    MPI_Allreduce(&ml, &mg, 1, MPI_INT, MPI_MAX, comm);
    m = max(mg, 1);
  }
  int lerr = (nto && uo.N( ) != m*nto) || Do.N( ) != nto, gerr = 0;
  MPI_Allreduce(&lerr, &gerr, 1, MPI_INT, MPI_MAX, comm);
  if (gerr) ExecError("remapP0ddm(Tho, uo[], Do[], Th, u[]): uo must be P0 on Tho (or [P0,...,P0]), Do P0 on Tho");
  // the limited slopes of the owned old cells (the cells which share a vertex: the overlap gives them)
  KN< int > vbeg(nvo + 1), vtri(3*nto);
  vbeg = 0;
  for (int k = 0; k < nto; ++k) for (int l = 0; l < 3; ++l) vbeg[Tho(k, l) + 1]++;
  for (int i = 0; i < nvo; ++i) vbeg[i + 1] += vbeg[i];
  {
    KN< int > fill(nvo);
    fill = 0;
    for (int k = 0; k < nto; ++k) for (int l = 0; l < 3; ++l) { int i = Tho(k, l); vtri[vbeg[i] + fill[i]++] = k; }
  }
  KN< R2 > G(nto);
  for (int k = 0; k < nto; ++k) G[k] = (R2(Tho[k][0]) + R2(Tho[k][1]) + R2(Tho[k][2]))/3.;
  // the records sent: the 3 vertices (counterclockwise), then for each component its value and its slope
  const int rs = 6 + 3*m;
  vector< double > own;// the records of the owned old cells
  vector< double > obox;// their bounding boxes
  {
    KN< int > mark(nto);
    mark = -1;
    vector< int > nb;
    vector< double > r1(m), r2(m), umin(m), umax(m), gx(m), gy(m);
    for (int k = 0; k < nto; ++k) {
      if (Do[k] < 0.5) continue;
      nb.clear( );
      mark[k] = k;
      for (int l = 0; l < 3; ++l) {
        int i = Tho(k, l);
        for (int p = vbeg[i]; p < vbeg[i + 1]; ++p) if (mark[vtri[p]] != k) { mark[vtri[p]] = k; nb.push_back(vtri[p]); }
      }
      double m11 = 0, m12 = 0, m22 = 0;
      for (int c = 0; c < m; ++c) { r1[c] = r2[c] = 0; umin[c] = umax[c] = uo[m*k + c]; gx[c] = gy[c] = 0; }
      for (size_t q = 0; q < nb.size( ); ++q) {
        int j = nb[q];
        R2 d = G[j] - G[k];
        double w = 1./(d.x*d.x + d.y*d.y);
        m11 += w*d.x*d.x; m12 += w*d.x*d.y; m22 += w*d.y*d.y;
        for (int c = 0; c < m; ++c) {
          double du = uo[m*j + c] - uo[m*k + c];
          r1[c] += w*d.x*du; r2[c] += w*d.y*du;
          umin[c] = min(umin[c], uo[m*j + c]); umax[c] = max(umax[c], uo[m*j + c]);
        }
      }
      double dt = m11*m22 - m12*m12;
      if (nb.size( ) >= 2 && dt > 1e-12*(m11 + m22)*(m11 + m22)) {
        double a = 1.;// the limiter, the same for all the components
        for (int c = 0; c < m; ++c) {
          gx[c] = (m22*r1[c] - m12*r2[c])/dt; gy[c] = (m11*r2[c] - m12*r1[c])/dt;
          for (int l = 0; l < 3; ++l) {
            R2 P(Tho[k][l]);
            double dv = gx[c]*(P.x - G[k].x) + gy[c]*(P.y - G[k].y), uk = uo[m*k + c];
            if (dv > 1e-300) a = min(a, (umax[c] - uk)/dv);
            else if (dv < -1e-300) a = min(a, (umin[c] - uk)/dv);
          }
        }
        a = max(0., a);
        for (int c = 0; c < m; ++c) { gx[c] *= a; gy[c] *= a; }
      }
      R2 P[3] = {Tho[k][0], Tho[k][1], Tho[k][2]};
      if (Det(P[0], P[1], P[2]) < 0) swap(P[1], P[2]);
      for (int l = 0; l < 3; ++l) { own.push_back(P[l].x); own.push_back(P[l].y); }
      for (int c = 0; c < m; ++c) { own.push_back(uo[m*k + c]); own.push_back(gx[c]); own.push_back(gy[c]); }
      obox.push_back(min(P[0].x, min(P[1].x, P[2].x))); obox.push_back(max(P[0].x, max(P[1].x, P[2].x)));
      obox.push_back(min(P[0].y, min(P[1].y, P[2].y))); obox.push_back(max(P[0].y, max(P[1].y, P[2].y)));
    }
  }
  // the bounding boxes of the new local meshes (all the processes)
  double bb[4] = {1e300, -1e300, 1e300, -1e300};
  for (int i = 0; i < Th.nv; ++i) {
    bb[0] = min(bb[0], Th(i).x); bb[1] = max(bb[1], Th(i).x);
    bb[2] = min(bb[2], Th(i).y); bb[3] = max(bb[3], Th(i).y);
  }
  double eps = 1e-9*max(1., max(bb[1] - bb[0], bb[3] - bb[2]));
  if (Th.nv) { bb[0] -= eps; bb[1] += eps; bb[2] -= eps; bb[3] += eps; }
  vector< double > abb(4*mpisize);
  MPI_Allgather(bb, 4, MPI_DOUBLE, abb.data( ), 4, MPI_DOUBLE, comm);
  // the records sent to each process: the owned old cells which meet its box
  const int no = obox.size( )/4;
  vector< int > scnt(mpisize, 0), sdsp(mpisize, 0), rcnt(mpisize, 0), rdsp(mpisize, 0);
  vector< double > sbuf;
  for (int r = 0; r < mpisize; ++r) {
    const double *B = &abb[4*r];
    sdsp[r] = sbuf.size( );
    for (int q = 0; q < no; ++q) {
      const double *b = &obox[4*q];
      if (b[1] < B[0] || b[0] > B[1] || b[3] < B[2] || b[2] > B[3]) continue;
      sbuf.insert(sbuf.end( ), own.begin( ) + q*rs, own.begin( ) + (q + 1)*rs);
    }
    scnt[r] = sbuf.size( ) - sdsp[r];
  }
  MPI_Alltoall(scnt.data( ), 1, MPI_INT, rcnt.data( ), 1, MPI_INT, comm);
  long rtot = 0;
  for (int r = 0; r < mpisize; ++r) { rdsp[r] = rtot; rtot += rcnt[r]; }
  vector< double > rbuf(max(rtot, 1L));
  MPI_Alltoallv(sbuf.data( ), scnt.data( ), sdsp.data( ), MPI_DOUBLE, rbuf.data( ), rcnt.data( ), rdsp.data( ), MPI_DOUBLE, comm);
  const int nr = rtot/rs;
  // a grid of buckets on the box of the new local mesh: the received cells in each bucket
  const int ng = max(1, int(sqrt(double(max(nr, 1)))/2) + 1);
  const double gx0 = bb[0], gy0 = bb[2], hx = max(bb[1] - bb[0], 1e-300)/ng, hy = max(bb[3] - bb[2], 1e-300)/ng;
  auto cellx = [&](double x) { return min(ng - 1, max(0, int((x - gx0)/hx))); };
  auto celly = [&](double y) { return min(ng - 1, max(0, int((y - gy0)/hy))); };
  vector< vector< int > > bucket(ng*ng);
  vector< R2 > rG(nr);
  for (int q = 0; q < nr; ++q) {
    const double *p = &rbuf[q*rs];
    rG[q] = R2((p[0] + p[2] + p[4])/3., (p[1] + p[3] + p[5])/3.);
    int i0 = cellx(min(p[0], min(p[2], p[4]))), i1 = cellx(max(p[0], max(p[2], p[4])));
    int j0 = celly(min(p[1], min(p[3], p[5]))), j1 = celly(max(p[1], max(p[3], p[5])));
    for (int i = i0; i <= i1; ++i) for (int j = j0; j <= j1; ++j) bucket[i + ng*j].push_back(q);
  }
  // the exact means on the new cells
  KN< double > &u = *pu;
  u.resize(m*Th.nt);
  vector< int > seen(nr, -1);
  vector< double > su(m);
  for (int kn = 0; kn < Th.nt; ++kn) {
    R2 A[3] = {Th[kn][0], Th[kn][1], Th[kn][2]};
    if (Det(A[0], A[1], A[2]) < 0) swap(A[1], A[2]);
    int i0 = cellx(min(A[0].x, min(A[1].x, A[2].x))), i1 = cellx(max(A[0].x, max(A[1].x, A[2].x)));
    int j0 = celly(min(A[0].y, min(A[1].y, A[2].y))), j1 = celly(max(A[0].y, max(A[1].y, A[2].y)));
    double sa = 0;
    for (int c = 0; c < m; ++c) su[c] = 0;
    for (int i = i0; i <= i1; ++i)
      for (int j = j0; j <= j1; ++j)
        for (int q : bucket[i + ng*j]) {
          if (seen[q] == kn) continue;
          seen[q] = kn;
          const double *p = &rbuf[q*rs];
          R2 B[3] = {R2(p[0], p[1]), R2(p[2], p[3]), R2(p[4], p[5])}, C;
          double a = TriInter(A, B, C);
          if (a <= 0) continue;
          sa += a;
          for (int c = 0; c < m; ++c) {
            const double *v = p + 6 + 3*c;
            su[c] += a*(v[0] + v[1]*(C.x - rG[q].x) + v[2]*(C.y - rG[q].y));
          }
        }
    if (sa > 0) { for (int c = 0; c < m; ++c) u[m*kn + c] = su[c]/sa; continue; }
    // not covered by the old cells (a new border): the nearest received cell
    R2 Gk = (A[0] + A[1] + A[2])/3.;
    int qb = -1;
    double db = 1e300;
    for (int q = 0; q < nr; ++q) { double d = (rG[q].x - Gk.x)*(rG[q].x - Gk.x) + (rG[q].y - Gk.y)*(rG[q].y - Gk.y); if (d < db) { db = d; qb = q; } }
    for (int c = 0; c < m; ++c) u[m*kn + c] = qb >= 0 ? rbuf[qb*rs + 6 + 3*c] : 0.;
  }
  return 0.;
}

static void Load_Init() {
  Global.Add("adaptmeshddm", "(", new OneOperatorCode<AdaptationDDM>);
  Global.Add("remapP0ddm", "(", new OneOperator5_< double, pmesh, KN_< double >, KN_< double >, pmesh, KN< double > * >(RemapP0DDM));
  if (verbosity > 1 && !mpirank) cout << " load: adaptmesh_ddm (adaptmeshddm)" << endl;
}

LOADFUNC(Load_Init)
