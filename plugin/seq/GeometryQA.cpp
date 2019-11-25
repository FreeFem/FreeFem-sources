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
// AUTHORS : Jean-Marie Mirebeau
// E-MAIL  : jean-marie.mirebeau@math.u-psud.fr

// *INDENT-OFF* //
// ff-c++-LIBRARY-dep:
// ff-c++-cpp-dep:
// *INDENT-ON* //
#define _USE_MATH_DEFINES
#include <iostream>
#include <cfloat>
#include <cmath>
using namespace std;

#include "ff++.hpp"
using namespace Fem2D;

// #include <boost/operators.hpp>
namespace mir {
#define _FLAGGED_BOUNDARY_
#include "Geometry.hpp"

  ostream_math coutMath = cout << Mathematica;

  template<>
  const string R2::name = "R2";
  template<>
  const string Z2::name = "Z2";
  template<>
  const R2 R2::NABiDim =
    R2(DBL_MAX,
       DBL_MAX);    // NAN could be a better alternative, but equality tests become problematic,
  template<>
  const Z2 Z2::NABiDim = Z2(INT_MAX, INT_MAX);    // and NAN does not exists for integers.
  template<>
  const string R3::name = "R3";
  template<>
  const string Z3::name = "Z3";

#include "ExampleMetrics.h"

  Edge *Edge::refine(Tab< Edge > &EdgeAllocator, Tab< Vertex > &VertexAllocator,
                     const Metric2 &metric, refinement_priority priority) {
    Edge *const wf = which_first(priority);

    if (wf != this) {
      wf->refine(EdgeAllocator, VertexAllocator, metric, priority);
    }

    if (sister != NULL) {
      Edge *const swf = sister->which_first(priority);
      if (swf != sister) {
        swf->refine(EdgeAllocator, VertexAllocator, metric, priority);
      }
    }

    Vertex const *w = next->v;
    Vertex *t = VertexAllocator.next( );
    const int newGen =
      1 + max(max(u->getGen( ), v->getGen( )),
              max(next->v->getGen( ), sister == NULL ? -1 : sister->next->v->getGen( )));
    *t = Vertex((*u + *v) / 2, newGen, metric);
    Edge *const wt = EdgeAllocator.next( ), *const tw = EdgeAllocator.next( ),
                *const ut = EdgeAllocator.next( );
    *wt = Edge(w, t, tw, this);
    *tw = Edge(t, w, wt, prev( ));

#ifdef _FLAGGED_BOUNDARY_
    *ut = Edge(u, t, NULL, tw, onBoundary( ));
#else
    *ut = Edge(u, t, NULL, tw);
#endif

    u = t;
    prev( )->next = ut;
    next->next = wt;

    if (sister == NULL) {
      return ut;
    }

    Vertex const *x = sister->next->v;
    Edge *xt = EdgeAllocator.next( ), *tx = EdgeAllocator.next( ), *vt = EdgeAllocator.next( );
    *xt = Edge(x, t, tx, sister);
    *tx = Edge(t, x, xt, sister->prev( ));

#ifdef _FLAGGED_BOUNDARY_
    *vt = Edge(v, t, this, tx, onBoundary( ));
#else
    *vt = Edge(v, t, this, tx);
#endif

    sister->u = t;
    sister->prev( )->next = vt;
    sister->next->next = xt;

    ut->sister = sister;
    sister->sister = ut;
    sister = vt;

    return ut;
  }

  bool Edge::hRefine3(double h, Tab< Edge > &EdgeAllocator, Tab< Vertex > &VertexAllocator,
                      const Metric2 &metric, refinement_priority priority) {
    Edge *const wf = which_first(priority);

    if (this != wf) {
      return wf->hRefine3(h, EdgeAllocator, VertexAllocator, metric, priority);
    }

    const double maxLen = max(max(vec( ).norm( ), next->vec( ).norm( )), prev( )->vec( ).norm( ));
    R2 const *w = next->v;
    double minSize = sqrt(metric(*w).invNorm( ));

    if (metric.lip == 0) {
      if (h * minSize < maxLen) {
        refine(EdgeAllocator, VertexAllocator, metric, priority);
        return true;
      }

      return false;
    }

    // If the metric is not constant, it is sampled until resolution allows to make a sensible
    // choice. Smallest eigenvalue used only.
    for (int pow = 1; h * (minSize - metric.lip * maxLen / (2 * pow)) < maxLen / 2;
         pow *= 2) {    // !!! condition should be checked
      for (int i = 0; i <= pow; i++) {
        for (int j = 0; i + j <= pow; j++) {
          if (i % 2 == 0 && j % 2 == 0) {
            continue;
          }

          minSize =
            min(minSize, sqrt(metric((*u * i + *v * j + *w * (pow - i - j)) / pow).invNorm( )));
          if (h * minSize < maxLen) {
            refine(EdgeAllocator, VertexAllocator, metric, priority);
            return true;
          }
        }
      }
    }

    return false;
  }

  Edge *Edge::hRefine2(double h, Tab< Edge > &EdgeAllocator, Tab< Vertex > &VertexAllocator,
                       const Metric2 &metric, safe_vector< Edge * > *recursive, bool exaggerate) {
    // in the case recursive=!NULL of a recursive split, the safe_vector recursive contains all
    // newly created edges which are oriented like *this.
    const R2 Vec = vec( );
    sym2 m = metric(*u);

    if (exaggerate) {
      m = m.exaggerate( );
    }

    double minSize = 1 / m.norm(Vec);

    if (metric.lip == 0) {
      if (h * minSize < 1) {
        Edge *const e = refine(EdgeAllocator, VertexAllocator, metric, selected_edge_first);
        if (recursive) {
          hRefine2(h, EdgeAllocator, VertexAllocator, metric, recursive, exaggerate);
          e->hRefine2(h, EdgeAllocator, VertexAllocator, metric, recursive, exaggerate);
          recursive->push_back(e);
        }

        return e;
      }

      return NULL;
    }

    // If the metric is not constant, it is sampled until resolution allows to make a sensible
    // choice (based on lip constant). Smallest eigenvalue used only.
    for (int pow = 1; h * (minSize - metric.lip / (2 * pow)) < 0.5;
         pow *= 2) {    // !!! condition should be checked
      for (int i = 0; i <= pow; i++) {
        if (i % 2 == 0) {
          continue;
        }

        m = metric((*u * i + *v * (pow - i)) / pow);
        if (exaggerate) {
          m = m.exaggerate( );
        }

        minSize = min(minSize, 1 / m.norm(Vec));
        if (h * minSize < 1) {
          Edge *const e = refine(EdgeAllocator, VertexAllocator, metric, selected_edge_first);
          if (recursive) {
            hRefine2(h, EdgeAllocator, VertexAllocator, metric, recursive, exaggerate);
            e->hRefine2(h, EdgeAllocator, VertexAllocator, metric, recursive, exaggerate);
            recursive->push_back(e);
          }

          return e;
        }
      }
    }

    return NULL;
  }

  Edge *Edge::which_first(refinement_priority priority) {
    if (priority == selected_edge_first) {
      return this;
    } else if (priority == newest_vertex_first) {
      const int ugen = u->getGen( ), vgen = v->getGen( ), wgen = next->v->getGen( );
      return ugen > vgen ? (ugen > wgen ? next : this) : (vgen > wgen ? prev( ) : this);
    } else {    // euclidean_longest_edge_first
      const double this_len = vec( ).norm( ), next_len = next->vec( ).norm( ),
                   prev_len = prev( )->vec( ).norm( );
      return next_len > prev_len ? (next_len > this_len ? next : this)
                                 : (prev_len > this_len ? prev( ) : this);
    }
  }

  bool Edge::cut(Vertex const *start, Vertex const *end, Tab< Edge > &EdgeAllocator,
                 Tab< Vertex > &VertexAllocator, const Metric2 &metric, vector< Edge * > &aligned) {
    // aligned collects created edges which are aligned with the cut
    if (start == v) {
      return next->cut(start, end, EdgeAllocator, VertexAllocator, metric,
                       aligned);    // looking for the starting edge, in the wrong direction
    }

    if (start == u) {    // looking for the starting edge, in the good direction
      if (end == v) {
        return false;
      }

      const R2 se = *end - *start;
      Edge *currentEdge = this;
      double currentDet;
      double prevDet = -det(vec( ), se);

      do {
        currentDet = -prevDet;
        prevDet = det(currentEdge->prev( )->vec( ), se);
        if (currentDet > 0 && prevDet > 0) {
          return currentEdge->cut(start, end, NULL, EdgeAllocator, VertexAllocator, metric,
                                  aligned);
        }

        currentEdge = currentEdge->prev( )->sister;
        if (currentEdge == this) {
          return false;
        }
      } while (currentEdge != NULL);

      currentEdge = this;
      currentDet = det(vec( ), se);

      while (currentEdge->sister != NULL) {
        currentEdge = currentEdge->sister->next;
        if (currentEdge == this) {
          return false;
        }

        prevDet = -currentDet;
        currentDet = det(currentEdge->vec( ), se);
        if (currentDet > 0 && prevDet > 0) {
          return currentEdge->cut(start, end, NULL, EdgeAllocator, VertexAllocator, metric,
                                  aligned);
        }
      }
    }

    return false;    // not supposed to happen
  }

  bool Edge::cut(Vertex const *start, Vertex const *end, Edge *oldSister,
                 Tab< Edge > &EdgeAllocator, Tab< Vertex > &VertexAllocator, const Metric2 &metric,
                 vector< Edge * > &aligned) {
    Vertex *t = next->intersect(start, end, VertexAllocator, metric);

    if (oldSister == NULL) {
      if (t == NULL) {
        return false;
      }

      Edge *const tw = next;
      Edge *const wu = tw->next;
      Edge *const vt = EdgeAllocator.next( );
      Edge *const ut = EdgeAllocator.next( ), *const tu = EdgeAllocator.next( );
      tw->u = t;
      wu->next = ut;
      *ut = Edge(u, t, tu, next);
      aligned.push_back(ut);
      *tu = Edge(t, u, ut, this);
      *vt = Edge(v, t, next->sister, tu);

      vt->sister->sister = vt;
      next = vt;
      return vt->sister->cut(start, end, tw, EdgeAllocator, VertexAllocator, metric, aligned);
    }

    Vertex const *const w = next->v;
    Vertex const *const s = sister->v;
    if (t != NULL) {
      Edge *const tw = next;
      Edge *const wu = next->next;
      Edge *const ws = EdgeAllocator.next( );
      Edge *const sw = EdgeAllocator.next( );
      Edge *const st = EdgeAllocator.next( );
      Edge *const ts = EdgeAllocator.next( );
      Edge *const us = EdgeAllocator.next( );
      Edge *const vt = EdgeAllocator.next( );

      tw->u = t;
      tw->next = ws;
      wu->next = us;
      *ws = Edge(w, s, sw, st);
      *sw = Edge(s, w, ws, wu);
      *st = Edge(s, t, ts, tw);
      aligned.push_back(st);
      *ts = Edge(t, s, st, this);
      *us = Edge(u, s, oldSister, sw);
      oldSister->sister = us;
      *vt = Edge(v, t, tw->sister, ts);

      vt->sister->sister = vt;
      next = vt;
      u = s;

      return vt->sister->cut(start, end, tw, EdgeAllocator, VertexAllocator, metric, aligned);
    }

    t = next->next->intersect(start, end, VertexAllocator, metric);
    if (t != NULL) {
      Edge *const vw = next;
      Edge *const tu = next->next;
      Edge *const sw = EdgeAllocator.next( );
      Edge *const ws = EdgeAllocator.next( );
      Edge *const ts = EdgeAllocator.next( );
      Edge *const st = EdgeAllocator.next( );
      Edge *const us = EdgeAllocator.next( );
      Edge *const wt = EdgeAllocator.next( );

      vw->next = ws;
      tu->u = t;
      tu->next = us;
      *sw = Edge(s, w, ws, wt);
      *ws = Edge(w, s, sw, this);
      *ts = Edge(t, s, st, sw);
      *st = Edge(s, t, ts, tu);
      aligned.push_back(st);
      *us = Edge(u, s, oldSister, st);
      oldSister->sister = us;
      *wt = Edge(w, t, tu->sister, ts);

      wt->sister->sister = wt;
      u = s;

      return wt->sister->cut(start, end, tu, EdgeAllocator, VertexAllocator, metric, aligned);
    }

    if (w == end) {    // last bisected edge
      Edge *const vw = next;
      Edge *const wu = next->next;
      Edge *const us = EdgeAllocator.next( );
      Edge *const sw = EdgeAllocator.next( );
      Edge *ws = EdgeAllocator.next( );

      vw->next = ws;
      wu->next = us;
      *us = Edge(u, s, oldSister, sw);
      oldSister->sister = us;
      *sw = Edge(s, w, ws, wu);
      aligned.push_back(sw);
      *ws = Edge(w, s, sw, this);

      u = s;

      return true;
    }

    return false;
  }

  Vertex *Edge::intersect(Vertex const *start, Vertex const *end, Tab< Vertex > &VertexAllocator,
                          const Metric2 &metric) {
    if (start == end || start == u || start == v || end == u || end == v || u == v) {
      return NULL;
    }

    const R2 es = *start - *end, uv = vec( );
    const R2 diff = (*end + *start) - (*v + *u);
    if (det(uv, es) == 0) {
      return NULL;
    }

    const R2 coef = diff.lin_solve(uv, es);
    if (coef.x <= -1 || coef.x >= 1 || coef.y <= -1 || coef.y >= 1 || coef == R2::NABiDim) {
      return NULL;
    }

    *VertexAllocator.next( ) = Vertex(*u * (1 - coef.x) / 2 + *v * (1 + coef.x) / 2,
                                      1 + max(u->getGen( ), v->getGen( )), metric);
    return &VertexAllocator[VertexAllocator.max_accessed_pos];
  }

  ostream_math operator<<(ostream_math f, Edge *e) {
    if (e != NULL) {
      return f << *e;
    }

    return f;
  }

  Triangulation::Triangulation(int N, const Metric2 &Metric) : metric(Metric) {
    for (int i = 0; i <= N; i++) {
      for (int j = 0; j <= N; j++) {
        vertices[i + (N + 1) * j] =
          Vertex(R2(i / double(N), j / double(N)), abs(N - i - j), metric);
      }
    }

    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        const int index = i + N * j;
#ifdef _FLAGGED_BOUNDARY_
        edges[6 * index + 0] =
          Edge(&vertices[i + (N + 1) * j], &vertices[i + 1 + (N + 1) * j],
               j > 0 ? &edges[6 * (index - N) + 3] : NULL, &edges[6 * index + 1], j == 0 ? 1 : 0);
        edges[6 * index + 1] =
          Edge(&vertices[i + 1 + (N + 1) * j], &vertices[i + (N + 1) * (j + 1)],
               &edges[6 * index + 4], &edges[6 * index + 2]);
        edges[6 * index + 2] =
          Edge(&vertices[i + (N + 1) * (j + 1)], &vertices[i + (N + 1) * j],
               i > 0 ? &edges[6 * (index - 1) + 5] : NULL, &edges[6 * index + 0], i == 0 ? 4 : 0);
        edges[6 * index + 3] = Edge(
          &vertices[i + 1 + (N + 1) * (j + 1)], &vertices[i + (N + 1) * (j + 1)],
          j + 1 < N ? &edges[6 * (index + N)] : NULL, &edges[6 * index + 4], j + 1 == N ? 3 : 0);
        edges[6 * index + 4] =
          Edge(&vertices[i + (N + 1) * (j + 1)], &vertices[i + 1 + (N + 1) * j],
               &edges[6 * index + 1], &edges[6 * index + 5]);
        edges[6 * index + 5] =
          Edge(&vertices[i + 1 + (N + 1) * j], &vertices[i + 1 + (N + 1) * (j + 1)],
               i + 1 < N ? &edges[6 * (index + 1) + 2] : NULL, &edges[6 * index + 3],
               i + 1 == N ? 2 : 0);
#else
        edges[6 * index + 0] =
          Edge(&vertices[i + (N + 1) * j], &vertices[i + 1 + (N + 1) * j],
               j > 0 ? &edges[6 * (index - N) + 3] : NULL, &edges[6 * index + 1]);
        edges[6 * index + 1] =
          Edge(&vertices[i + 1 + (N + 1) * j], &vertices[i + (N + 1) * (j + 1)],
               &edges[6 * index + 4], &edges[6 * index + 2]);
        edges[6 * index + 2] =
          Edge(&vertices[i + (N + 1) * (j + 1)], &vertices[i + (N + 1) * j],
               i > 0 ? &edges[6 * (index - 1) + 5] : NULL, &edges[6 * index + 0]);
        edges[6 * index + 3] =
          Edge(&vertices[i + 1 + (N + 1) * (j + 1)], &vertices[i + (N + 1) * (j + 1)],
               j + 1 < N ? &edges[6 * (index + N)] : NULL, &edges[6 * index + 4]);
        edges[6 * index + 4] =
          Edge(&vertices[i + (N + 1) * (j + 1)], &vertices[i + 1 + (N + 1) * j],
               &edges[6 * index + 1], &edges[6 * index + 5]);
        edges[6 * index + 5] =
          Edge(&vertices[i + 1 + (N + 1) * j], &vertices[i + 1 + (N + 1) * (j + 1)],
               i + 1 < N ? &edges[6 * (index + 1) + 2] : NULL, &edges[6 * index + 3]);
#endif
      }
    }

    movie_init( );
  }

  void Triangulation::Delaunay_ordered(const vector< bool > &toExclude) {
    vector< double > gains;
    gains.resize(ne_oriented( ));
    SortedList< RZ > toFlip;

    for (int i = 0; i < ne_oriented( ); i++) {    // first sort and first computation of the angles
      if (!edges[i].isRepresentative( )) {
        continue;
      }

      const double gain = toExclude[i] ? 0 : edges[i].flipGain( );
      gains[i] = gain;
      if (gain > 0) {
        toFlip.insert(RZ(gain, i));
      }
    }

    while (toFlip.Card( ) > 0) {
      const RZ current = toFlip.pop( );
      Edge *affected[4];
      if (!edges[current.number].flip(affected)) {
        continue;
      }

      movie_frame( );

      for (int i = 0; i < 4; i++) {
        const int index = edges.index(affected[i]);
        if (gains[index] > 0) {
          toFlip.remove(RZ(gains[index], index));
        }

        const double gain = toExclude[index] ? 0 : edges[index].flipGain( );
        gains[index] = gain;
        if (gain > 0) {
          toFlip.insert(RZ(gain, index));
        }
      }
    }
  }

  void Triangulation::hRefineQA(double h, unsigned int flag, Edge::refinement_priority priority) {
    if (h <= 0) {
      return;
    }

    if (!(flag & hRQA_noIsoRef)) {
      hRefine(h, priority);
    }

    const bool exportIntermediateData = flag & hRQA_exportIntermediateData;
    const bool finalRefine = flag & hRQA_finalRefine;

    // const double ratio=0.5; //over_refinement of the extremal half edges. Replaced with an
    // exageration of the metric.

    if (exportIntermediateData) {
      vertices.export_content("vertices_iso.txt");
      Tab< Z2 > connectivity;
      Connectivity(connectivity);
      connectivity.export_content("connectivity_iso.txt");
    }

    // collecting the extremal edges associated to each vertex.
    const int nv_iso = nv( );
    const int ne_iso = ne_oriented( );
    cout << "Triangulation::hRefineQA : Intermediate isotropic triangulation contains " << nt( )
         << " triangles" << endl;

    double *minDet = new double[nv_iso]( );    // set to zero
    double *maxDet = new double[nv_iso]( );
    Edge **minEdge = new Edge *[nv_iso]( );
    Edge **maxEdge = new Edge *[nv_iso]( );
    R2 *eigenVec = new R2[nv_iso];
    cout << " nv = " << nv( ) << endl;

    for (int i = 0; i < nv_iso; i++) {
      eigenVec[i] = vertices[i].getm( ).eigensys( );
    }

    for (int i = 0; i < ne_iso; i++) {
      Edge *e = &edges[i];
      if (!e->isRepresentative( )) {
        continue;
      }

      R2 vec = e->vec( );
      vec /= vec.norm( );

      for (int j = 0; j <= 1; j++) {
        const int index = vertices.index(j == 0 ? e->getu( ) : e->getv( ));
        const double currentDet = det(vec, eigenVec[index]) * (1 - 2 * j);
        if (currentDet < minDet[index]) {
          minDet[index] = currentDet;
          minEdge[index] = e;
        } else if (currentDet > maxDet[index]) {
          maxDet[index] = currentDet;
          maxEdge[index] = e;
        }
      }
    }

    delete[] minDet;
    delete[] maxDet;
    delete[] eigenVec;

    // refining these edges
    Tab< Edge > halfEdges;
    Tab< int > endVertex;
    safe_vector< Edge * >
      toExcludePtr;    // these edges should be excluded from the flipping process

    for (int i = 0; i < ne_iso; i++) {
      Edge *e = &edges[i];
      if (!e->isRepresentative( )) {
        continue;
      }

      const int indexu = vertices.index(e->getu( )), indexv = vertices.index(e->getv( ));
      bool extru = indexu < nv_iso && ((e == minEdge[indexu]) || (e == maxEdge[indexu]));
      bool extrv = indexv < nv_iso && ((e == minEdge[indexv]) || (e == maxEdge[indexv]));

      if (!extru && !extrv) {
        continue;
      }

      Edge *const f = e->hRefine2(h, edges, vertices, metric, NULL);    // non recursive split.

      if (f == NULL) {
        continue;
      }

      const int indexw =
        vertices
          .max_accessed_pos;    // index of the mid point between u and v that was just created

      if (exportIntermediateData) {
        if (extru) {
          *halfEdges.next( ) = *f;
        }

        if (extrv) {
          *halfEdges.next( ) = *e;
        }
      }

      if (finalRefine) {
        if (extru && !extrv) {
          endVertex[indexw] =
            indexv + 1;    // "+1" required to distinguish from values initialized to zero
        }

        if (extrv && !extru) {
          endVertex[indexw] = indexu + 1;
        }
      }

      if (extru) {    // recursive split, exaggerate //*ratio
        f->hRefine2(h, edges, vertices, metric, &toExcludePtr, true);
        toExcludePtr.push_back(f);
      }

      if (extrv) {
        e->hRefine2(h, edges, vertices, metric, &toExcludePtr, true);
        toExcludePtr.push_back(e);
      }

      // Note : in the theoretical algorithm, the edge is split in 3, and the midpart is refined if
      // extru || extrv
      movie_frame( );
    }

    delete[] minEdge;
    delete[] maxEdge;
    if (exportIntermediateData) {
      halfEdges.export_content("halfEdges.txt");
      edges.export_content("edgesBeforeDelaunay.txt");
    }

    cout << "Triangulation::hRefineQA : Intermediate anisotropic triangulation contains " << nt( )
         << " triangles." << endl;

    vector< bool > toExclude;
    toExclude.resize(ne_oriented( ));

    for (vector< Edge * >::const_iterator e = toExcludePtr.begin( ); e != toExcludePtr.end( );
         ++e) {
      if (!(*e)->onBoundary( )) {
        toExclude[edges.index(*e)] = true;
        toExclude[edges.index((*e)->getSister( ))] = true;
      }    // boundary edges are fixed anyway
    }

    toExcludePtr.clear( );

    if (movie_name.size( ) > 0) {
      cout << "Beginning main flip. Movie frame : " << movie_frame_number;
    }

    Delaunay_ordered(toExclude);
    if (movie_name.size( ) > 0) {
      cout << "; Finishing : " << movie_frame_number << endl;
    }

    if (finalRefine) {
      for (int i = 0; i < ne_oriented( ); ++i) {
        Edge *e = &edges[i];
        const int indexw = vertices.index(e->getu( ));
        if (endVertex[indexw] > 0) {
          e->cut(e->getu( ), &vertices[endVertex[indexw] - 1], edges, vertices, metric,
                 toExcludePtr);
          endVertex[indexw] = 0;
          movie_frame( );
        }
      }

      cout << "Triangulation::hRefineQA : Triangulation contains " << nt( )
           << " triangles after optional refinement (which eliminates the remaining large angles.)"
           << endl;
    }

    // return;
    // now taking care of the boundary
    for (int i = 0; i < ne_oriented( ); ++i) {
      if (edges[i].onBoundary( )) {
        if (edges[i].hRefine2(h, edges, vertices, metric, &toExcludePtr, true)) {
          movie_frame( );
        }
      }
    }

    toExclude.resize(ne_oriented( ));

    for (vector< Edge * >::const_iterator e = toExcludePtr.begin( ); e != toExcludePtr.end( );
         ++e) {
      if (!(*e)->onBoundary( )) {
        toExclude[edges.index(*e)] = true;
        toExclude[edges.index((*e)->getSister( ))] = true;
      }    // boundary edges are fixed anyway
    }

    Delaunay_ordered(toExclude);
    cout << "Triangulation::hRefineQA : Final triangulation contains " << nt( )
         << " triangles after refinement of the boundary." << endl;
  }

  Fem2D::Mesh *Triangulation::export_to_Mesh( ) const {
    typedef Fem2D::Triangle FFT;
    typedef Fem2D::Vertex FFV;
    using Fem2D::R2;
    typedef Fem2D::BoundaryEdge FFBE;
    using Fem2D::Mesh;
    using namespace Fem2D;
    // using  Fem2D::R;

    vector< bool > onBoundary;
    onBoundary.resize(nv( ));
    int boundaryEdges = 0;

    for (int i = 0; i < ne_oriented( ); i++) {
      const Edge &e = edges[i];
      if (!e.onBoundary( ) || !e.isRepresentative( )) {
        continue;
      }

      onBoundary[vertices.index(e.getu( ))] = true;
      onBoundary[vertices.index(e.getv( ))] = true;
      boundaryEdges++;
    }

    int nbv = nv( );            // nombre de sommet
    int nbt = nt( );            // nombre de triangles
    int neb = boundaryEdges;    // nombre d'aretes fontiere
    // allocation des nouveaux items du maillage
    FFV *v = new FFV[nbv + nbt];
    FFT *t = new FFT[nbt * 3];
    FFBE *b = new FFBE[neb];
    // generation des nouveaus sommets
    FFV *vv = v;

    for (int i = 0; i < nbv; i++) {
      vv->x = vertices[i].x;
      vv->y = vertices[i].y;
      vv->lab = onBoundary[i];
      vv++;
    }

    // generation des triangles
    FFT *tt = t;
    int nberr = 0;
    Vertex const *triangle[3];

    for (int i = 0; i < ne_oriented( ); i++) {
      if (edges[i].isRepresentative3(triangle)) {
        int i0 = vertices.index(triangle[0]);
        int i1 = vertices.index(triangle[1]);
        int i2 = vertices.index(triangle[2]);
        (*tt++).set(v, i0, i1, i2, 0);
      }
    }

    // les arete frontieres qui n'ont pas change
    FFBE *bb = b;

    // edges
    for (int i = 0; i < ne_oriented( ); i++) {
      const Edge &e = edges[i];
      if (!e.onBoundary( ) || !e.isRepresentative( )) {
        continue;
      }

      // data_out << 1+vertices.index(e.getu())
      int i1 = vertices.index(e.getu( ));
      int i2 = vertices.index(e.getv( ));
      int lab = e.onBoundary( );
      (*bb++).set(v, i1, i2, lab);
    }

    Mesh *m = new Mesh(nbv, nbt, neb, v, t, b);

    return m;
  }

  void Triangulation::export_to_FreeFem(const char *filename) const {
    ofstream data_out;

    data_out.open(filename);

    vector< bool > onBoundary;
    onBoundary.resize(nv( ));
    int boundaryEdges = 0;

    for (int i = 0; i < ne_oriented( ); i++) {
      const Edge &e = edges[i];
      if (!e.onBoundary( ) || !e.isRepresentative( )) {
        continue;
      }

      onBoundary[vertices.index(e.getu( ))] = true;
      onBoundary[vertices.index(e.getv( ))] = true;
      boundaryEdges++;
    }

    data_out << nv( ) << " " << nt( ) << " " << boundaryEdges << endl;

    // vertices
    for (int i = 0; i < nv( ); i++) {
      data_out << vertices[i] << " " << onBoundary[i] << endl;
    }

    // triangles
    Vertex const *triangle[3];

    for (int i = 0; i < ne_oriented( ); i++) {
      if (edges[i].isRepresentative3(triangle)) {
        data_out << 1 + vertices.index(triangle[0]) << " " << 1 + vertices.index(triangle[1]) << " "
                 << 1 + vertices.index(triangle[2]) << " " << 0 << endl;
      }
    }

    cout << "Exporting edges" << endl;

    // edges
    for (int i = 0; i < ne_oriented( ); i++) {
      const Edge &e = edges[i];
      if (!e.onBoundary( ) || !e.isRepresentative( )) {
        continue;
      }

#ifdef _FLAGGED_BOUNDARY_
      data_out << 1 + vertices.index(e.getu( )) << " " << 1 + vertices.index(e.getv( )) << " "
               << e.onBoundary( ) << endl;
#else
      data_out << 1 + vertices.index(e.getu( )) << " " << 1 + vertices.index(e.getv( )) << " " << 1
               << endl;
#endif
    }

    data_out.close( );
  }

  string Triangulation::movie_frame_name( ) const {
    ostringstream oss;

    oss << movie_name << "_";
    if (movie_frame_number < 10) {
      oss << 0;
    }

    if (movie_frame_number < 100) {
      oss << 0;
    }

    if (movie_frame_number < 1000) {
      oss << 0;
    }

    oss << movie_frame_number++ << ".txt";
    return oss.str( );
  }

#ifdef FF___HPP_
  Triangulation::Triangulation(const Fem2D::Mesh &Th, const Metric2 &Metric) : metric(Metric) {
    for (int i = 0; i < Th.nv; ++i) {
      const Fem2D::R2 Point = Th(i);
      vertices[i] = Vertex(R2(Point.x, Point.y), 0, metric);
    }

    cout << "Hello ???" << endl;

    // I collect boundary edge labels, since I do not understand how Th.BorderElementAdj(ui,vi)
    // works.
    std::map< pair< int, int >, int > BorderElementLabels;

    for (int i = 0; i < Th.nbBrdElmts( ); ++i) {
      const int ui = Th(Th.be(i)[0]), vi = Th(Th.be(i)[1]);
      const int label = Th.be(i).lab;
      BorderElementLabels[pair< int, int >(ui, vi)] = label;
      BorderElementLabels[pair< int, int >(vi, ui)] = label;
    }

    for (int i = 0; i < Th.nt; ++i) {
      for (int j = 0; j < 3; ++j) {
        int ip;
        int jp = j;
        ip = Th.ElementAdj(i, jp);
#ifdef _FLAGGED_BOUNDARY_
        Edge *const sister = (i == ip) ? NULL : &edges[3 * ip + jp];
        const int ui = Th(i, (j + 1) % 3), vi = Th(i, (j + 2) % 3);
        const std::map< pair< int, int >, int >::iterator it =
          BorderElementLabels.find(pair< int, int >(ui, vi));
        const int label = (it == BorderElementLabels.end( )) ? 0 : it->second;
        edges[3 * i + j] =
          Edge(&vertices[ui], &vertices[vi], sister, &edges[3 * i + (j + 1) % 3], label);
#else
        edges[3 * i + j] = Edge(&vertices[Th(i, (j + 1) % 3)], &vertices[Th(i, (j + 2) % 3)],
                                i == ip ? NULL : &edges[3 * ip + jp], &edges[3 * i + (j + 1) % 3]);
#endif
      }
    }
  }

#endif
}    // end namespace mir
