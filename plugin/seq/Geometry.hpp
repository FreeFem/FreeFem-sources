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
// SUMMARY : DelaunayFlip
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Jean-Marie Mirebeau
// E-MAIL  : jean-marie.mirebeau@math.u-psud.fr

#ifndef GEOMETRY
#define GEOMETRY

#include <iostream>
#include "RZ.h"
#include "SortedList.h"

using namespace std;

/**************** Vertex ******************/
class Vertex : public R2 {    // heritage is perhaps not totally justified here, but...
  sym2 m;
  int gen;
  // friend int main(int argc, const char * argv[]); //for Debug purposes

 public:
  Vertex( ) : R2( ), m( ){};
  Vertex(R2 u, int Gen, sym2 metric = sym2( )) : R2(u), gen(Gen), m(metric) {}

  Vertex(R2 u, int Gen, const Metric2 &metric) : R2(u), gen(Gen), m(metric(u)) {}

  const sym2 &getm( ) const { return m; }

  int getGen( ) const { return gen; }

  sym3 homogeneousDistance( ) const {
    const R2 a(*this);
    const R2 ma = m * a;
    const double ama = a.scal(ma);
    return sym3(m.xx, m.yy, ama, m.xy, -ma.y, -ma.x);
  }

  double dist2(const R2 &P) const { return m.norm2(*this - P); }
};

/*inline ostream& operator<<(ostream &f, const Vertex &u){f<<R2(u)<<u.getGen(); return f;}
 * inline ostream_math operator<<(ostream_math f, const Vertex &u){
 *  if(f.format==Mathematica) f<<"{"<<R2(u)<<","<<u.getGen()<<"}"; else f.os<<u; return f;}
 */
inline ostream &operator<<(ostream &f, const Vertex &u) { return f << R2(u); }

inline ostream_math operator<<(ostream_math f, const Vertex &u) {
  return f << "{" << R2(u) << "," << u.getGen( ) << "," << u.getm( ) << "}";
}

/****************** Edge ****************/
class Edge {    // edges are oriented. t is on the left.
                // bool refinable;
  Vertex const *u;
  Vertex const *
    v;    // Note : routines tend to leave v constant when possible (i.e. except construct and flip)
  Edge *next, *sister;    // next edge on triangle and reversed edge
  Edge *prev( ) const { return next->next; }

  bool cut(Vertex const *start, Vertex const *end, Edge *oldSister, Tab< Edge > &EdgeAllocator,
           Tab< Vertex > &VertexAllocator, const Metric2 &metric, vector< Edge * > &aligned);
  // friend int main(int argc, const char * argv[]); //for Debug purposes
#ifdef _FLAGGED_BOUNDARY_
  int onBoundary_;
#endif

 public:
  Edge( ) : u(NULL), v(NULL), next(NULL), sister(NULL){};

#ifdef _FLAGGED_BOUNDARY_
  Edge(Vertex const *U, Vertex const *V, Edge *S, Edge *N, int OnBoundary_ = 0)
    : u(U), v(V), sister(S), next(N), onBoundary_(OnBoundary_){};
  int onBoundary( ) const {
    return onBoundary_;
  }    // {if(sister!=NULL) return 0; if(onBoundary_>0) return onBoundary_; return 10;}

#else
  Edge(Vertex const *U, Vertex const *V, Edge *S, Edge *N) : u(U), v(V), next(N), sister(S){};
  bool onBoundary( ) const { return sister == NULL; }

#endif

  inline bool flipable( ) const;
  inline double flipGain( ) const;    // gain brought by a flip in terms of square cosine
  inline bool flip( );
  inline bool flip(Edge *affected[4]);
  inline bool flip_resolve( );
  Vertex const *getu( ) const { return u; }

  Vertex const *getv( ) const { return v; }

  R2 vec( ) const { return *v - *u; }    // R2(v->x - u->x, v->y - u->y)

  bool isRepresentative( ) const {
    return (sister == NULL) || (u->x < v->x) || (u->x == v->x && u->y < v->y);
  }

  Edge *representative( ) { return isRepresentative( ) ? this : sister; }

  bool isRepresentative3( ) const {
    R2 v = vec( );
    return v < next->vec( ) && v < prev( )->vec( );
  }

  bool isRepresentative3(Vertex const *triangle[3]) const {
    triangle[0] = u;
    triangle[1] = v;
    triangle[2] = next->v;
    return isRepresentative3( );
  }

  enum refinement_priority {
    selected_edge_first,
    newest_vertex_first,
    euclidean_longest_edge_first
  };
  Edge *which_first(refinement_priority priority);
  // refine splits the edge and returns a pointer to the other half of the splitted edge //bool
  // newestVertexFirst=true
  Edge *refine(Tab< Edge > &EdgeAllocator, Tab< Vertex > &VertexAllocator, const Metric2 &metric,
               refinement_priority priority);
  // splits associated triangle if larger than size h/sqrt(smallEigenVal) (metric taken at the worst
  // point on the triangle). Non recursive.
  bool hRefine3(double h, Tab< Edge > &EdgeAllocator, Tab< Vertex > &VertexAllocator,
                const Metric2 &metric, refinement_priority priority);
  // splits the edge if its length in the metric (taken at the worst point on the edge) is larger
  // than h.
  Edge *hRefine2(double h, Tab< Edge > &EdgeAllocator, Tab< Vertex > &VertexAllocator,
                 const Metric2 &metric, safe_vector< Edge * > *recursive = NULL,
                 bool exaggerate = false);
  bool check( ) const;

  Edge *getNext( ) const { return next; }

  Edge *getSister( ) const { return sister; }

  bool cut(Vertex const *start, Vertex const *end, Tab< Edge > &EdgeAllocator,
           Tab< Vertex > &VertexAllocator, const Metric2 &metric, vector< Edge * > &aligned);
  Vertex *intersect(Vertex const *start, Vertex const *end, Tab< Vertex > &VertexAllocator,
                    const Metric2 &metric);
};

inline ostream &operator<<(ostream &f, const Edge &e) {
  f << R2(*(e.getu( ))) << " " << R2(*(e.getv( )));
  return f;
}

inline ostream_math operator<<(ostream_math f, const Edge &e) {
  if (f.format == Mathematica) {
    f << "{" << R2(*(e.getu( ))) << "," << R2(*(e.getv( ))) << "}";
  } else {
    f.os << e;
  }

  return f;
}

inline bool Edge::check( ) const {
  if (u == NULL || v == NULL) {
    cout << "Edge::check : Invalid extremities";
  } else if (u == v) {
    cout << "Edge::check : identical extremities";
  } else if (next == NULL || next->next == NULL) {
    cout << "Edge::check : Missing edge connections";
  } else if (next->next->next != this) {
    cout << "Edge::check : not a triangle";
  } else if (next->u != v) {
    cout << "Edge::check : invalid next edge (next->u!=v)";
  } else if (sister != NULL && sister->u != v) {
    cout << "Edge::check : invalid sister edge";
  }
  // else if(sister!=NULL && refinable && !sister->refinable) cout<<"Edge::check : flag refinable
  // inconsistent with sister edge"<<endl;
  else if (isRepresentative3( ) && det(vec( ), next->vec( )) < 0) {
    cout << "Edge::check : trigonometric order not respected";
  }

#ifdef _FLAGGED_BOUNDARY_
  else if (sister == NULL && onBoundary_ == 0) {
    cout << "Edge::check : Interior edge without sister !" << endl;
  }
#endif
  else {
    return true;
  }
  coutMath << " " << *this << *next << *next->next << endl;
  return false;
}

inline bool Edge::flipable( ) const {
  return !onBoundary( ) && det(sister->prev( )->vec( ), next->vec( )) > 0 &&
         det(prev( )->vec( ), sister->next->vec( )) > 0;
}

inline double Edge::flipGain( ) const {    // edge is assumed to be flipable, hence angles are <pi.
                                           // Result positive <=> flip increases minimum angle.
  if (!flipable( )) {
    return 0;
  }

  Vertex const *s = next->v;
  Vertex const *t = sister->next->v;
  const R2 uv = *v - *u, st = *t - *s, vs = *s - *v, su = *u - *s, ut = *t - *u, tv = *v - *t;
  const sym2 &mu = u->getm( ), &mv = v->getm( ), ms = s->getm( ), mt = t->getm( );
  return min(min(min(-ms.cos(st, vs), ms.cos(st, su)), min(mt.cos(st, ut), -mt.cos(st, tv))),
             min(-mu.cos(su, ut), -mv.cos(tv, vs))) -
         min(min(min(-mu.cos(uv, su), mu.cos(uv, ut)), min(mv.cos(uv, tv), -mv.cos(uv, vs))),
             min(-ms.cos(vs, su), -mt.cos(ut, tv)));
}

inline bool Edge::flip( ) {
  if (sister == NULL) {
    return false;    // triangle inversion is unchecked
  }

  Edge *e = this, *f = sister;
  Edge *ep = prev( ), *fp = sister->prev( );
  Vertex const *u = ep->u;
  Vertex const *v = fp->u;

  e->u = u;
  e->v = v;
  f->u = v;
  f->v = u;

  e->next->next = e;
  f->next->next = f;
  ep->next = f->next;
  fp->next = e->next;
  e->next = fp;
  f->next = ep;
  return true;
}

inline bool Edge::flip(Edge *affected[4]) {    // representatives of affected edges
  if (flip( )) {
    affected[0] = next->representative( );
    affected[1] = prev( )->representative( );
    affected[2] = sister->next->representative( );
    affected[3] = sister->prev( )->representative( );
    return true;
  }

  return false;
}

inline bool Edge::flip_resolve( ) {
  if (flipGain( ) <= 0) {
    return false;
  }

  Edge *affected[4];
  flip(affected);

  for (int i = 0; i < 4; i++) {
    affected[i]->flip_resolve( );
  }

  return true;
}

/****************************************************/

class Triangulation {
  Tab< Vertex > vertices;
  Tab< Edge > edges;
  friend int main(int argc, const char *argv[]);    // for Debug purposes

 public:
  const Metric2 &metric;
  const Tab< Vertex > &getVertices( ) { return vertices; }

  const Tab< Edge > &getEdges( ) { return edges; }

  int nv( ) const { return vertices.card( ); }

  int ne_oriented( ) const { return edges.card( ); }

  int nt( ) const { return ne_oriented( ) / 3; }

  Triangulation(int N, const Metric2 &Metric);    // basic NxN square triangulation
  // Triangulation(const Tab<Vertex> &Vertices, const Tab<Edge> &Edges, sym2 (*Metric)(const
  // R2&)=NULL, double Lip=5): vertices(Vertices), edges(Edges), metric(Metric), lip(Lip)
  // {if(!check()) cout << "Invalid triangulation !" << endl; movie_init();} //pointer mismatch
#ifdef FF___HPP_
  Triangulation(const Fem2D::Mesh &Th, const Metric2 &Metric);
#endif

  void Delaunay_ordered(const vector< bool > &toExclude);
  void Delaunay_ordered( ) {
    vector< bool > toExclude;
    toExclude.resize(ne_oriented( ));
    Delaunay_ordered(toExclude);
  }

  void Delaunay_unordered( ) {
    for (int i = 0; i < ne_oriented( ); ++i) {
      edges[i].flip_resolve( );
    }
  }

  int Connectivity(Tab< Z2 > &connectivity) const {
    int counter = 0;

    for (int i = 0; i < ne_oriented( ); ++i) {
      if (edges[i].isRepresentative( )) {
        connectivity[counter++] =
          Z2(vertices.index(edges[i].getu( )), vertices.index(edges[i].getv( )));
      }
    }

    return nt( );    // number of triangles
  }

  // Refines triangles isotropically, based on the small eigenvalue of the metric.
  void hRefine(double h = 1,
               Edge::refinement_priority priority = Edge::euclidean_longest_edge_first) {
    if (h <= 0) {
      return;
    }

    for (int i = 0; i < ne_oriented( ); ++i) {
      if (edges[i].hRefine3(h, edges, vertices, metric, priority)) {
        movie_frame( );
      }
    }
  }

  void hRefineQA(double h = 1, unsigned int flag = 0,
                 Edge::refinement_priority priority = Edge::euclidean_longest_edge_first);
  enum hRefineQA_opt { hRQA_finalRefine = 1, hRQA_exportIntermediateData = 2, hRQA_noIsoRef = 4 };

  // Affects only the position, not the metric
  void moveMesh(R2 (*vec)(const R2 &), double amplification = 1) {
    for (int i = 0; i < nv( ); i++) {
      vertices[i] += vec(vertices[i]) * amplification;
    }
  }

  bool check( ) const {
    bool passed = true;

    for (int i = 0; i < edges.card( ); i++) {
      passed = passed && edges[i].check( );
    }

    return passed;
  }

  void export_to_FreeFem(const char *filename) const;
  Fem2D::Mesh *export_to_Mesh( ) const;
  void export_to_Mathematica(const char *filename) const {
    ofstream data_out;
    data_out.open(filename);
    data_out << Mathematica << edges;
    data_out.close( );
  }

  void export_to_Mathematica_Metric(const char *filename) const {
    ofstream data_out;
    data_out.open(filename);
    data_out << Mathematica << vertices;
    data_out.close( );
  }

  void export_to(Format_Math format, const char *filename) const {
    if (format == Mathematica) {
      export_to_Mathematica(filename);
    } else {
      export_to_FreeFem(filename);
    }
  }

  // movie functionality
  string movie_name;
  Format_Math movie_format;
  mutable int movie_frame_number;
  void movie_frame( ) const {
    if (movie_name.size( ) == 0) {
      return;
    }

    export_to(movie_format, movie_frame_name( ).c_str( ));
  }

 private:
  void movie_init( ) {
    movie_name = "";
    movie_frame_number = 0;
    movie_format = Mathematica;
  }

  string movie_frame_name( ) const;
};

#endif
