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

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

// fvdualmesh.hpp: the dual meshes for the finite volumes with the unknowns at the vertices (P1, FreeVol)
//   class dualmesh  (type meshdual, 2D): the dual mesh of F. Bouchut, the control volumes built with the centers of
//                   the inscribed circles of the triangles (the default of builddualmesh(Th))
//   class dualmeshN (types meshdualL 1D, meshdual3 3D): the median dual mesh, the middles of the edges and the
//                   centroids of the elements (at the end of this file)

#ifndef FVDUALMESH_HPP
#define FVDUALMESH_HPP

#include <ff++.hpp>
#include <map>
#include <utility>
#include <vector>
#include "AFunction.hpp"
#include "AFunction_ext.hpp"    // Extension of "AFunction.hpp" to deal with more than 3 parameters function
using namespace Fem2D;

// normalized part of dualmesh
class dmnormalized {
public :
  // attributes (number limited to 10)
  KN<long>* ndvcv; //(array of size 3) containing  [0] ndv: number of dual vertices, [1] ncv: number of dual control volumes, [2] ndbe: number of dual boundary elements
  KN<double>* dvxy; //array (size 2*(Th.nt+2*nbdrypieces)) of x,y coordinates of dual vertices : incenters of each triangle of primal mesh, and boundary old+new vertices (x,y components written adjacently, the list is by couple)
  KN<long>* ninterf; //array (size Th.nv) of number of interfaces around each vertex=control volume
  KNM<long>* dvcv; //array (size (sizeneip,Th.nv)), for each control volume and each of its interface, dual vertex number of origin in anticlockwise order. One has dvcv(j,i)=dualpiece01((1-orientation(j,i))/2,interf(j,i))
  KNM<long>* dvbe; //array (size (2,ndbe)), dual vertices for dual boundary elements
  KN<long>* dbelab; //array (size ndbe), dual boundary element labels
  KNM<long>* interf; //array (size (sizeneip,Th.nv)), for each control volume, list of the interfaces around it, in anticlockwise order
  KNM<long>* orientation; //array (size (sizeneip,Th.nv)), for each control volume and each of its interface, orientation = +1 if the control volume is Left with respect to interface, -1 if it is Right.
  // functions
  // constructors
  dmnormalized() : ndvcv(0),dvxy(0),ninterf(0),dvcv(0)
                   ,dvbe(0),dbelab(0),interf(0),orientation(0) {
   if ((!mpirank)&&(verbosity>1)) cout << " dmnormalized::init " << this << endl;}
  // destructors
  ~dmnormalized(){destroy();}
  void destroy(){
   if ((!mpirank)&&(verbosity>1)) cout << " dmnormalized::destroy " << this << endl;
   if (ndvcv) {delete ndvcv; ndvcv=0;}
   if (dvxy) {delete dvxy; dvxy=0;}
   if (ninterf) {delete ninterf; ninterf=0;}
   if (dvcv) {delete dvcv; dvcv=0;}
   if (dvbe) {delete dvbe; dvbe=0;}
   if (dbelab) {delete dbelab; dbelab=0;}
   if (interf) {delete interf; interf=0;}
   if (orientation) {delete orientation; orientation=0;}
  }
};// end of class dmnormalized

// extra quantities associated to dualmesh
class dmextras {
public :
  // attributes (number limited to 10)
  KN<double>* dmcutedge;//array (size 2*nedge) of x,y components of point of cut on edge (written adjacently, the list is by couple)
  KN<double>* dmnormal;//unitary normal to interfaces of control volumes, oriented from Left vertex to Right vertex
  KN<double>* dmlengthedge;//length of interface
  KN<double>* dmareaminLR;//smallest area between L and R control volumes (used for CFL condition)
  KN<long>*   dmnbneit;//
  KN<double>* dmvolcvprod;//
  KN<double>* dmNbdryvert;//mean unitary normal to primal boundary vertices, or 0 if not on physical boundary
  KN<double>* dmmidedge;//array (size 2*ndualedge) of x,y components of middle of dual edge (written adjacently, the list is by couple)
  // functions
  // constructors
  dmextras() : dmcutedge(0),dmnormal(0),dmlengthedge(0),dmareaminLR(0)
               ,dmnbneit(0),dmvolcvprod(0),dmNbdryvert(0),dmmidedge(0) {
   if ((!mpirank)&&(verbosity>1)) cout << " dmextras::init " << this << endl;}
  // destructors
  ~dmextras(){destroy();}
  void destroy(){
   if ((!mpirank)&&(verbosity>1)) cout << " dmextras::destroy " << this << endl;
   if (dmcutedge) {delete dmcutedge; dmcutedge=0;}
   if (dmnormal) {delete dmnormal; dmnormal=0;}
   if (dmlengthedge) {delete dmlengthedge; dmlengthedge=0;}
   if (dmareaminLR) {delete dmareaminLR; dmareaminLR=0;}
   if (dmnbneit) {delete dmnbneit; dmnbneit=0;}
   if (dmvolcvprod) {delete dmvolcvprod; dmvolcvprod=0;}
   if (dmNbdryvert) {delete dmNbdryvert; dmNbdryvert=0;}
   if (dmmidedge) {delete dmmidedge; dmmidedge=0;}
  }
};// end of class dmextras

// matrices associated to dualmesh
class dmmatrices {
public :
  // attributes (number limited to 10)
  MatriceMorse<double>* matvert2grav; //matrix to change vertex values to gravity centers values, making conversion from FE representation to FV representation
  MatriceMorse<double>* mat_piece2cvL; //Left connectivity matrix between pieces of boundary of control volumes (=interfaces) and control volumes
  MatriceMorse<double>* mat_piece2cvR; //Right connectivity matrix between pieces of boundary of control volumes (=interfaces) and control volumes
  MatriceMorse<double>* mat_dxygrav; //matrix for computing slopes in triangles in terms of values at gravity centers of control volumes
                                     // this means dx(), dy() but with respect to gravity centers instead of vertices of Th
  MatriceMorse<double>* mat_piece2cvdxyL; //Left,Right matrices for computing values at middle of interfaces (used for second-order reconstruction)
  MatriceMorse<double>* mat_piece2cvdxyR; //defined by vector joining gravity center of control volume to the middle of the interface
  MatriceMorse<double>* mat_piece2cvprodL; //matrices for gathering the interfaces left/right values at vertices to perform limitation (used for second-order reconstruction)
  MatriceMorse<double>* mat_piece2cvprodR;
  MatriceMorse<double>* mat_tri2cvprod; //matrix for gathering the triangle slopes at vertices to perform slope limitation (used for second-order reconstruction)
  // functions
  // constructors
  dmmatrices() : matvert2grav(0),mat_piece2cvL(0),mat_piece2cvR(0),mat_dxygrav(0)
                 ,mat_piece2cvdxyL(0),mat_piece2cvdxyR(0),mat_piece2cvprodL(0),mat_piece2cvprodR(0),mat_tri2cvprod(0) {
   if ((!mpirank)&&(verbosity>1)) cout << " dmmatrices::init " << this << endl;}
  // destructors
  ~dmmatrices(){destroy();}
  void destroy(){
   if ((!mpirank)&&(verbosity>1)) cout << " dmmatrices::destroy " << this << endl;
   if (matvert2grav) {delete matvert2grav; matvert2grav=0;}
   if (mat_piece2cvL) {delete mat_piece2cvL; mat_piece2cvL=0;}
   if (mat_piece2cvR) {delete mat_piece2cvR; mat_piece2cvR=0;}
   if (mat_dxygrav) {delete mat_dxygrav; mat_dxygrav=0;}
   if (mat_piece2cvdxyL) {delete mat_piece2cvdxyL; mat_piece2cvdxyL=0;}
   if (mat_piece2cvdxyR) {delete mat_piece2cvdxyR; mat_piece2cvdxyR=0;}
   if (mat_piece2cvprodL) {delete mat_piece2cvprodL; mat_piece2cvprodL=0;}
   if (mat_piece2cvprodR) {delete mat_piece2cvprodR; mat_piece2cvprodR=0;}
   if (mat_tri2cvprod) {delete mat_tri2cvprod; mat_tri2cvprod=0;}
  }
};// end of class dmmatrices
//--------------------------------------------------
// definition of class dualmesh
//------------------------------------------------
// the mesh of the control volumes of a dual mesh (fespace Wstar(Thstar, P0)) is no more a dual mesh (for the plot:
// fvdualmesh.cpp)
void DualSubForget(const void *sub);
// the dual meshes and their primal mesh (fvdualmesh.cpp), from build to destroy: fvscheme FV([V1, ...], eq=...)
// without the dual mesh finds the dual mesh built in the script on the mesh of V1 (FreeVol)
void DualOfPrimalAdd(const void *primal, const void *dual, int dim);
void DualOfPrimalForget(const void *dual);
void DualsOfPrimal(const void *primal, int dim, std::vector< const void * > &duals);// the dual meshes of primal
bool DualIsBuilt(const void *dual);// the dual mesh is built and not destroyed (fvscheme: the dual mesh still used)
// a dual mesh keeps its primal mesh with a reference (increment in build, destroy in destroy): the primal mesh is
// not freed while its dual mesh is used, so a new mesh cannot have its address
class dualmesh {
public :
  // attributes (number limited to 10)
  const Mesh * primal; // pointer to primal mesh
  dmnormalized* normalized;// normalized part of dualmesh
  KN<long>* nbed; // (array of size 3) containing  [0] nedge: number of edges of the primal mesh, ndualedge: number of edges of the dual mesh. [1] ndualedge=nedge+2*nbdrypieces, with nbdrypieces the number of "true boundary" elements of the primal mesh. [2] maximum number sizeneip of interfaces around a control volume
  KN<long>* edgenumbering; // (size 3*Th.nt) index of primal edges, in terms of triangles and their 3 edges (each triangle has 3 adjacent values in array)
  KN<long>* cvLR; // (size 2*ndualedge) array of indices of Left and Right control volumes, in terms of dual edge index (Left and Right are written adjacently, the list is by couple)
  KNM<long>* dualpiece01; // (size (2,ndualedge)) origin and end of piece of boundary of control volume, in terms of dual edge index (value is dual vertex index, origin and end indexed by the first int (0 or 1))
  KN<double>* cvarea; // (size Th.nv) area of control volumes
  KN<double>* cvgrav; // (size 2*Th.nv) coordinates of centers of gravity of control volumes (written adjacently, the list is by couple)
  dmextras* extras; // extra quantities of dualmesh
  dmmatrices* matrices; // matrices of dualmesh
  const Mesh *sub; // the control volumes as a mesh of FreeFEM (each one a fan of triangles, label = number + 1):
                   // fespace Wstar(Thstar, P0), built when it is used (fvdualmesh.cpp), 0 before
  // functions
  // constructors
  dualmesh() : primal(0),normalized(0),nbed(0),edgenumbering(0),cvLR(0)
               ,dualpiece01(0),cvarea(0),cvgrav(0),extras(0),matrices(0),sub(0) {
   if ((!mpirank)&&(verbosity>1)) cout << " dualmesh::init " << this << endl;}
  void init(){
   if ((!mpirank)&&(verbosity>1)) cout << " dualmesh::init " << this << endl;
   primal=0;normalized=0;nbed=0;edgenumbering=0;cvLR=0;
   dualpiece01=0;cvarea=0;cvgrav=0;extras=0;matrices=0;sub=0;
  }
  void build(const Mesh *const &pTh);
  // destructors
  ~dualmesh(){destroy();}
  void destroy(){
   if ((!mpirank)&&(verbosity>1)) cout << " dualmesh::destroy " << this << endl;
   if (normalized) {delete normalized; normalized=0;}
   if (nbed) {delete nbed; nbed=0;}
   if (edgenumbering) {delete edgenumbering; edgenumbering=0;}
   if (cvLR) {delete cvLR; cvLR=0;}
   if (dualpiece01) {delete dualpiece01; dualpiece01=0;}
   if (cvarea) {delete cvarea; cvarea=0;}
   if (cvgrav) {delete cvgrav; cvgrav=0;}
   if (extras) {delete extras; extras=0;}
   if (matrices) {delete matrices; matrices=0;}
   if (sub) {DualSubForget(sub); sub->destroy(); sub=0;}// destroyed when it is no more used by a fespace
   DualOfPrimalForget(this);
   if (primal) {primal->destroy(); primal=0;}// the reference of build
  }
};// end of class dualmesh

//---------------------------------
// build attributes of dualmesh
//---------------------------------
void dualmesh::build(const Mesh *const &pTh){
 if ((!mpirank)&&(verbosity>1)) cout << " building dualmesh " << endl;
 ffassert(pTh);
 pTh->increment();// the primal mesh kept while the dual mesh is used
 if (primal) primal->destroy();// built again without destroy (builddualmesh(Thstar, Th))
 primal=pTh;
 DualOfPrimalAdd(pTh,this,2);
 const Mesh &Th=*pTh;
 if (!(nbed)) nbed=new KN<long>;
 KN<long> &nbedtab=*(nbed);
 long nedge,ndualedge;
 KN<double> vxy;
 if (!(edgenumbering)) edgenumbering=new KN<long>;
 KN<long> &edgesmesh=*(edgenumbering);
 if (!(cvLR)) cvLR=new KN<long>;
 KN<long> &ctLR=*(cvLR);
 if (!(dualpiece01)) dualpiece01=new KNM<long>;
 KNM<long> &cvbdry=*(dualpiece01);
 if (!(cvarea)) cvarea=new KN<double>;
 KN<double> &volcv=*(cvarea);
 if (!(cvgrav)) cvgrav=new KN<double>;
 KN<double> &gravcv=*(cvgrav);
 KN<long> bepiece;// (size ndualedge) index of boundary piece as primal mesh boundary element (valid if piece is on true boundary)

 nedge=0;//number of edges of the primal mesh, the true value is computed below
 vxy.resize(2*Th.nt);//array of x,y coordinates of dual vertices : incenters of each triangle of primal mesh, and boundary old+new vertices (x,y components written adjacently, the list is by couple). The size will be increased afterwards
 edgesmesh.resize(3*Th.nt);//index of primal edges, in terms of triangles and their 3 edges (each triangle has 3 adjacent values in array)
 ctLR.resize(2*(3*Th.nt+2*Th.neb));//array of indices of Left and Right control volumes, in terms of dual edge index (Left and Right are written adjacently, the list is by couple). The size will be reduced afterwards
 ctLR=-1;// initialized to -1, this value means not (yet) defined
 cvbdry.resize(2,3*Th.nt+2*Th.neb);//origin and end of piece of boundary of control volume, in terms of dual edge index (value is dual vertex index, origin and end indexed by the first int (0 or 1)). The size will be reduced afterwards
 volcv.resize(Th.nv);// area of control volumes
 volcv=0.;
 gravcv.resize(2*Th.nv);// coordinates of centers of gravity of control volumes (written adjacently, the list is by couple)
 gravcv=0.;
 for (long k=0;k<Th.nt;++k) {//begin loop on triangles
  double v0x=Th(Th(k,0)).x;
  double v0y=Th(Th(k,0)).y;
  double v1x=Th(Th(k,1)).x;
  double v1y=Th(Th(k,1)).y;
  double v2x=Th(Th(k,2)).x;
  double v2y=Th(Th(k,2)).y;
  // incenter of the triangle
  double length0,length1,length2,w0x,w0y,w1x,w1y;
  length0=sqrt(pow(v2x-v1x,2)+pow(v2y-v1y,2));
  length1=sqrt(pow(v0x-v2x,2)+pow(v0y-v2y,2));
  length2=sqrt(pow(v1x-v0x,2)+pow(v1y-v0y,2));
  w0x=(v1x-v0x)/length2+(v2x-v0x)/length1;
  w0y=(v1y-v0y)/length2+(v2y-v0y)/length1;
  w1x=(v2x-v1x)/length0+(v0x-v1x)/length2;
  w1y=(v2y-v1y)/length0+(v0y-v1y)/length2;
  double lambda0=((v1x-v0x)*w1y-(v1y-v0y)*w1x)/(w0x*w1y-w0y*w1x);
  vxy[2*k]=v0x+lambda0*w0x;
  vxy[2*k+1]=v0y+lambda0*w0y;//end of def of incenter
  // primal edges (cut by line between two centers)
  int ee;
  for (int j=0;j<3;j++) {
   ee=j;
   long kee=Th.ElementAdj(k,ee);
   if (int(kee)>=k) {// other side of edge not yet known
    edgesmesh(3*k+j)=nedge;
    ctLR(2*nedge)=Th(k,(j+1)%3);// associate neighbouring vertex (->left control volume)
//    trofdof(nedge,1)=k;// associate triangle k
//    trofdof(nedge,2)=(j+1)%3;// associate vertex j
    if (int(kee)==k){// edge on boundary
     ctLR(2*nedge+1)=Th(k,(j+2)%3);// associate neighbouring vertex (->right control volume)
//     trofdof(nedge,4)=k;// associate triangle k
//     trofdof(nedge,5)=(j+2)%3;// associate vertex j
     cvbdry(0L,nedge)=-2;// temporary value!
     cvbdry(1L,nedge)=k;
     double ck1x=vxy[2*k];
     double ck1y=vxy[2*k+1];
     double ck2x=(Th(Th(k,(j+1)%3)).x+Th(Th(k,(j+2)%3)).x)/2.;
     double ck2y=(Th(Th(k,(j+1)%3)).y+Th(Th(k,(j+2)%3)).y)/2.;
     double ck3x=Th(Th(k,(j+1)%3)).x;
     double ck3y=Th(Th(k,(j+1)%3)).y;
     double areatri=abs((ck2x-ck1x)*(ck3y-ck1y)-(ck2y-ck1y)*(ck3x-ck1x))/2.;
     volcv(Th(k,(j+1)%3))+=areatri;
     gravcv(2*Th(k,(j+1)%3))+=areatri*(ck1x+ck2x+ck3x)/3.;
     gravcv(2*Th(k,(j+1)%3)+1)+=areatri*(ck1y+ck2y+ck3y)/3.;
     ck3x=Th(Th(k,(j+2)%3)).x;
     ck3y=Th(Th(k,(j+2)%3)).y;
     areatri=abs((ck2x-ck1x)*(ck3y-ck1y)-(ck2y-ck1y)*(ck3x-ck1x))/2.;
     volcv(Th(k,(j+2)%3))+=areatri;
     gravcv(2*Th(k,(j+2)%3))+=areatri*(ck1x+ck2x+ck3x)/3.;
     gravcv(2*Th(k,(j+2)%3)+1)+=areatri*(ck1y+ck2y+ck3y)/3.;
    }
    nedge++;
   }
   else {// the edge has already been visited from the other side
    ee=j;
    long k0=int(Th.ElementAdj(k,ee));
    double centerkx0=vxy[2*k];
    double centerky0=vxy[2*k+1];
    double centeradjx0=vxy[2*k0];
    double centeradjy0=vxy[2*k0+1];
    double scalinters0=((Th(Th(k,(j+1)%3)).x-centerkx0)*(Th(Th(k,(j+2)%3)).y-Th(Th(k,(j+1)%3)).y)-(Th(Th(k,(j+1)%3)).y-centerky0)*(Th(Th(k,(j+2)%3)).x-Th(Th(k,(j+1)%3)).x))/((centeradjx0-centerkx0)*(Th(Th(k,(j+2)%3)).y-Th(Th(k,(j+1)%3)).y)-(centeradjy0-centerky0)*(Th(Th(k,(j+2)%3)).x-Th(Th(k,(j+1)%3)).x));
    long noki=edgesmesh(3*k0+ee);
    edgesmesh(3*k+j)=noki;
    double cutedges2noki=centerkx0+scalinters0*(centeradjx0-centerkx0);
    double cutedges2nokip1=centerky0+scalinters0*(centeradjy0-centerky0);
    ctLR(2*noki+1)=Th(k,(j+1)%3);// associate neighbouring vertex (->right control volume)
//    trofdof(noki,4)=k;// associate triangle k
//    trofdof(noki,5)=(j+1)%3;// associate vertex j
    cvbdry(0L,noki)=k;
    cvbdry(1L,noki)=k0;
    double ck1x=centerkx0;
    double ck1y=centerky0;
    double ck2x=cutedges2noki;
    double ck2y=cutedges2nokip1;
    double ck3x=Th(Th(k,(j+1)%3)).x;
    double ck3y=Th(Th(k,(j+1)%3)).y;
    double areatri=abs((ck2x-ck1x)*(ck3y-ck1y)-(ck2y-ck1y)*(ck3x-ck1x))/2.;
    volcv(Th(k,(j+1)%3))+=areatri;
    gravcv(2*Th(k,(j+1)%3))+=areatri*(ck1x+ck2x+ck3x)/3.;
    gravcv(2*Th(k,(j+1)%3)+1)+=areatri*(ck1y+ck2y+ck3y)/3.;
    ck3x=Th(Th(k,(j+2)%3)).x;
    ck3y=Th(Th(k,(j+2)%3)).y;
    areatri=abs((ck2x-ck1x)*(ck3y-ck1y)-(ck2y-ck1y)*(ck3x-ck1x))/2.;
    volcv(Th(k,(j+2)%3))+=areatri;
    gravcv(2*Th(k,(j+2)%3))+=areatri*(ck1x+ck2x+ck3x)/3.;
    gravcv(2*Th(k,(j+2)%3)+1)+=areatri*(ck1y+ck2y+ck3y)/3.;
    ck1x=centeradjx0;
    ck1y=centeradjy0;
    areatri=abs((ck2x-ck1x)*(ck3y-ck1y)-(ck2y-ck1y)*(ck3x-ck1x))/2.;
    volcv(Th(k,(j+2)%3))+=areatri;
    gravcv(2*Th(k,(j+2)%3))+=areatri*(ck1x+ck2x+ck3x)/3.;
    gravcv(2*Th(k,(j+2)%3)+1)+=areatri*(ck1y+ck2y+ck3y)/3.;
    ck3x=Th(Th(k,(j+1)%3)).x;
    ck3y=Th(Th(k,(j+1)%3)).y;
    areatri=abs((ck2x-ck1x)*(ck3y-ck1y)-(ck2y-ck1y)*(ck3x-ck1x))/2.;
    volcv(Th(k,(j+1)%3))+=areatri;
    gravcv(2*Th(k,(j+1)%3))+=areatri*(ck1x+ck2x+ck3x)/3.;
    gravcv(2*Th(k,(j+1)%3)+1)+=areatri*(ck1y+ck2y+ck3y)/3.;
   }
  }//end of loop on j
 }//end of loop on triangles

 // centers of gravity of control volumes
 for (long i=0;i<Th.nv;i++) {
  gravcv(2*i)=gravcv(2*i)/volcv(i);
  gravcv(2*i+1)=gravcv(2*i+1)/volcv(i);
 }

 //boundary pieces
 bepiece.resize(nedge+2*Th.neb);//index of boundary piece as primal mesh boundary element (valid if piece is on true boundary). The size will be reduced afterwards
 bepiece=-1;//-1 means undefined. This will be the value for indices 0:nedge-1
 long nbdrypieces=0;// initialize the number of boundary pieces
 vxy.resize(2*(Th.nt+2*Th.neb));
 static const long maxcrossdepth=3L;//max number of crossing boundary components that will be sorted
 KNM<long> primalbdryvert2dualbdryvert(Th.nv,maxcrossdepth);//dual boundary vertex in terms of primal boundary vertex. The second index enables to manage with the case of crossing boundary
 primalbdryvert2dualbdryvert=-1L;
 KN<long> crossdepth(Th.nv);crossdepth=0L;// for each boundary primal vertex, number of crossing boundary components
 for (long ib=0;ib<Th.neb;ib++){// loop on boundary elements
  int iib;
  long kib=Th.BoundaryElement(ib, iib);//neighbouring triangle kib, edge in triangle iib
  int ee=iib;
  long kee=Th.ElementAdj(kib,ee);
  if (int(kee)==kib) {// if edge is really on boundary
//   ffassert (Th(kib,(iib+1)%3)==Th(Th.bedges[ib][0]));// problem of bug of freefem on the orientation of boundary elements for square()
//   ffassert (Th(kib,(iib+2)%3)==Th(Th.bedges[ib][1]));// problem of bug of freefem on the orientation of boundary elements for square()
   ctLR(2*(nedge+2*nbdrypieces))=Th(kib,(iib+1)%3);// associate neighbouring vertex (->left control volume)
//   trofdof(nedge+2*nbdrypieces,1)=kib;// associate triangle
//   trofdof(nedge+2*nbdrypieces,2)=(iib+1)%3;// associate vertex
   cvbdry(0L,nedge+2*nbdrypieces)=Th.nt+2*nbdrypieces;// gives the dual vertex number
   cvbdry(1L,nedge+2*nbdrypieces)=Th.nt+2*nbdrypieces+1;// gives the dual vertex number
   bepiece(nedge+2*nbdrypieces)=ib;// associate index of boundary element
   vxy(2*(Th.nt+2*nbdrypieces))=Th(Th(kib,(iib+1)%3)).x;
   vxy(2*(Th.nt+2*nbdrypieces)+1)=Th(Th(kib,(iib+1)%3)).y;
   ctLR(2*(nedge+2*nbdrypieces+1))=Th(kib,(iib+2)%3);// associate neighbouring vertex (->left control volume)
//   trofdof(nedge+2*nbdrypieces+1,1)=kib;// associate triangle
//   trofdof(nedge+2*nbdrypieces+1,2)=(iib+2)%3;// associate vertex
   cvbdry(0L,nedge+2*nbdrypieces+1)=Th.nt+2*nbdrypieces+1;// gives the dual vertex number
   cvbdry(1L,nedge+2*nbdrypieces+1)=Th(kib,(iib+2)%3);//temporary value : primal vertex
   bepiece(nedge+2*nbdrypieces+1)=ib;// associate index of boundary element
   vxy(2*(Th.nt+2*nbdrypieces+1))=(Th(Th(kib,(iib+1)%3)).x+Th(Th(kib,(iib+2)%3)).x)/2.;
   vxy(2*(Th.nt+2*nbdrypieces+1)+1)=(Th(Th(kib,(iib+1)%3)).y+Th(Th(kib,(iib+2)%3)).y)/2.;
   long m=edgesmesh(3*kib+iib);//index as internal edge
   cvbdry(0L,m)=Th.nt+2*nbdrypieces+1;//dual vertex index on middle of dual boundary edge
   long vert=Th(kib,(iib+1)%3);//primal vertex
   if (crossdepth(vert)<maxcrossdepth) {
    primalbdryvert2dualbdryvert(vert,crossdepth(vert))=Th.nt+2*nbdrypieces;
    crossdepth(vert)++;
   } else {//more than maxcrossdepth components for a crossing boundary
    crossdepth(vert)++;
    //cout << "builddualmesh warning: crossing boundary with multiple components? mpirank=" << mpirank << " vertex " << vert << " x=" << Th(vert).x << " y=" << Th(vert).y << endl;
   }
   nbdrypieces++;
  }
 }//end of loop on boundary elements
 // nbdrypieces can be less than Th.neb if there are internal boundaries
 ndualedge=nedge+2*nbdrypieces;//ndualedge and nbdrypieces have now their final value
 ctLR.resize(2*ndualedge);
 cvbdry.resize(2,ndualedge);
 bepiece.resize(ndualedge);
 nbedtab.resize(3);
 nbedtab[0]=nedge;
 nbedtab[1]=ndualedge;
 vxy.resize(2*(Th.nt+2*nbdrypieces));
 // completing the array cvbdry for the endpoints of dual boundary elements
 for (long ibp=0;ibp<nbdrypieces;ibp++) {// loop on real boundary elements
  long vert=cvbdry(1L,nedge+2*ibp+1);//primal vertex at end of boundary element
  if (crossdepth(vert)==1L) {
   cvbdry(1L,nedge+2*ibp+1)=primalbdryvert2dualbdryvert(vert,0L);//sets dual vertex number in the generic case of no crossing boundary
  } else {// crossing boundary
   long crossdp=min(crossdepth(vert),maxcrossdepth);
   double thetamin=2.*Pi;
   long lcompmin=0;
   long dvert=cvbdry(1L,nedge+2*ibp);
   double vectx=vxy(2*dvert)-Th(vert).x;
   double vecty=vxy(2*dvert+1)-Th(vert).y;
   double rvect=sqrt(vectx*vectx+vecty*vecty);
   vectx=vectx/rvect;
   vecty=vecty/rvect;
   for (long lcomp=0;lcomp<crossdp;lcomp++) {
    long dvert1=primalbdryvert2dualbdryvert(vert,lcomp);//dual vertex number
    double vect1x=vxy(2*(dvert1+1))-vxy(2*dvert1);
    double vect1y=vxy(2*(dvert1+1)+1)-vxy(2*dvert1+1);
    double rvect1=sqrt(vect1x*vect1x+vect1y*vect1y);
    vect1x=vect1x/rvect1;
    vect1y=vect1y/rvect1;
    double thet1=atan2(vect1y*vectx-vect1x*vecty,vect1x*vectx+vect1y*vecty);
    if (thet1<0.) thet1=thet1+2.*Pi;
    if (thet1<thetamin) {
     thetamin=thet1;
     lcompmin=lcomp;
    }
   }
   double dvertres=primalbdryvert2dualbdryvert(vert,lcompmin);
   cvbdry(1L,nedge+2*ibp+1)=dvertres;
   cout << "builddualmesh: crossing boundary mpirank " << mpirank << " vertex " << vert << " x=" << Th(vert).x << " y=" << Th(vert).y << " dual boundary element " << ibp;
   cout << ", dual vertex " << dvertres << ", x=" << vxy(2*dvert) << ", y=" << vxy(2*dvert+1) << endl;//coordinates of the middle of dual boundary element ibp
  }
 }
 KNM<long> dbe(2,2*nbdrypieces);//dual vertices for dual boundary elements
 for (long idbe=0;idbe<2*nbdrypieces;idbe++) {// loop on dual boundary elements
  dbe(0L,idbe)=cvbdry(0L,nedge+idbe);
  dbe(1L,idbe)=cvbdry(1L,nedge+idbe);
  //dbe(0L,2*ibp)=Th.nt+2*ibp;
  //dbe(1L,2*ibp)=Th.nt+2*ibp+1;
  //dbe(0L,2*ibp+1)=Th.nt+2*ibp+1;
  //dbe(1L,2*ibp+1)=cvbdry(1L,nedge+2*ibp+1);
 }

 KN<long> nbneipiece(Th.nv);// number of interfaces around each vertex=control volume
 nbneipiece=0L;
 long sizeneip=6;//max number of neighbouring interfaces (will be enlarged if necessary)
 KNM<long> listneip(Th.nv,sizeneip);//for each control volume, list of the interfaces around it
 listneip=-1L;
 KNM<long> orientvertinterf(Th.nv,sizeneip);//for each control volume and each of its interface, orientation = +1 if the control volume is Left with respect to interface, -1 if it is Right, 3 means indefined.
 orientvertinterf=3*1L;
 for (long m=0;m<ndualedge;m++) {// loop on left interfaces
  long vertm=ctLR(2*m);//trofdof(m,0);
  if (nbneipiece(vertm)+1>sizeneip) {
   sizeneip++;
   listneip.resize(Th.nv,sizeneip);
   orientvertinterf.resize(Th.nv,sizeneip);
  }
  listneip(vertm,nbneipiece(vertm))=m;
  orientvertinterf(vertm,nbneipiece(vertm))=1;
  nbneipiece(vertm)++;
 }
 for (long m=0;m<nedge;m++) {// loop on right interfaces. No right vertex on boundary
  long vertm=ctLR(2*m+1);//trofdof(m,3);
  if (nbneipiece(vertm)+1>sizeneip) {
   sizeneip++;
   listneip.resize(Th.nv,sizeneip);
   orientvertinterf.resize(Th.nv,sizeneip);
  }
  listneip(vertm,nbneipiece(vertm))=m;
  orientvertinterf(vertm,nbneipiece(vertm))=-1;
  nbneipiece(vertm)++;
 }
 nbedtab[2]=sizeneip;
 // reordering the interfaces so that one follows each other anticlockwise
 KN<long> listvertint(Th.nt+2*nbdrypieces);//working array indexed by dual vertices
 KN<long> listvertint01(Th.nt+2*nbdrypieces);//idem
 KNM<long> listneipord(sizeneip,Th.nv);//for each control volume, list of the interfaces around it, in anticlockwise order
 listneipord=-1L;
 KNM<long> orientvertinterford(sizeneip,Th.nv);//for each control volume and each of its interface, orientation = +1 if the control volume is Left with respect to interface, -1 if it is Right, 3 means indefined.
 orientvertinterford=3*1L;
 KNM<long> dualvertord(sizeneip,Th.nv);//for each control volume and each of its interface, dual vertex number of origin in anticlockwise order
 dualvertord=-1L;
 listvertint=-0L;// initialization is not mandatory, it avoids seg error in case of strange mesh
 listvertint01=-0L;// initialization is not mandatory, it avoids seg error in case of strange mesh
 for (long i=0;i<Th.nv;++i) {// loop on control volumes
  //listvertint=-0L;// initialization is not mandatory, it takes time, but avoids seg error in case of strange mesh
  //listvertint01=-0L;// initialization is not mandatory, it takes time, but avoids seg error in case of strange mesh
  if (crossdepth(i)<=maxcrossdepth) {//sorting the edges of the control volume i
   long ibegin=-1;
   for (long j=0;j<nbneipiece(i);++j) {// loop on interfaces of control volume i
    long m=listneip(i,j);//interface index
    long orim=orientvertinterf(i,j);//orientation +-1
    for (long i01=0;i01<2;i01++){// loop on origin and end of interface
     if ((cvbdry(i01,m)>=Th.nt)&&((cvbdry(i01,m)-Th.nt)%2==0)) ibegin=cvbdry(i01,m);// choose as first if it is a primal vertex on boundary
     if (i01==((1-orim)/2)) {//find the correct interface and origin/end to follow up
      listvertint(cvbdry(i01,m))=m;//associate interface index to dual vertex index
      listvertint01(cvbdry(i01,m))=i01;//associate origin or end
     }
    }
   }//end of loop on interfaces
   if (ibegin==-1) ibegin=cvbdry(0L,listneip(i,0L));//choose any dual vertex index to begin
   long dvertind=ibegin;//index of current dual vertex for loop
   for (long j=0;j<nbneipiece(i);++j) {// loop on dual vertices of control volume i
    long m=listvertint(dvertind);
    long i01=listvertint01(dvertind);
    listneipord(j,i)=m;
    orientvertinterford(j,i)=1-2*i01;
    dualvertord(j,i)=dvertind;
    dvertind=cvbdry(1L-i01,m);
   }//end of loop on dual vertices of control volume i
   //check that we have made the full loop around control volume
   if (dvertind!=ibegin) {cout << "problem in building dualmesh, is there an intricate crossing boundary? " << endl; cout << "mpirank="<< mpirank << " vertex " << i << " x=" << Th(i).x << " y=" << Th(i).y; ffassert(0);}
  } else {// keep non ordered data
   for (long j=0;j<nbneipiece(i);++j) {// loop on dual vertices of control volume i
    listneipord(j,i)=listneip(i,j);//keep unordered, but take care of interverted indices
    orientvertinterford(j,i)=orientvertinterf(i,j);//keep unordered, but take care of interverted indices
    dualvertord(j,i)=cvbdry((1-orientvertinterf(i,j))/2,listneip(i,j));
   }
   cout << "builddualmesh: mpirank " << mpirank << " warning non ordered control volume " << i << " x=" << Th(i).x << " y=" << Th(i).y << endl;
  }
 }//end of loop on control volumes

//////////////////////////////////////////////////////////
// normalized part
//////////////////////////////////////////////////////
 if (!(normalized)) normalized=new dmnormalized;
 dmnormalized &dmnor=*(normalized);
 if (!(dmnor.ndvcv)) dmnor.ndvcv=new KN<long>;
 KN<long> &nvcv=*(dmnor.ndvcv);
 nvcv.resize(3);
 nvcv[0]=Th.nt+ndualedge-nedge;//number of dual vertices
 nvcv[1]=Th.nv;//number of dual control volumes
 nvcv[2]=ndualedge-nedge;//number of dual boundary elements
 if (!(dmnor.dvxy)) dmnor.dvxy=new KN<double>;
 KN<double> &nvxy=*(dmnor.dvxy);
 nvxy.resize(vxy.N());
 nvxy=vxy;
 if (!(dmnor.ninterf)) dmnor.ninterf=new KN<long>;
 KN<long> &ninterfa=*(dmnor.ninterf);
 ninterfa.resize(Th.nv);
 ninterfa=nbneipiece;
 if (!(dmnor.dvcv)) dmnor.dvcv=new KNM<long>;
 KNM<long> &vcv=*(dmnor.dvcv);
 vcv.resize(sizeneip,Th.nv);
 vcv=dualvertord;
 if (!(dmnor.dvbe)) dmnor.dvbe=new KNM<long>;
 KNM<long> &vbe=*(dmnor.dvbe);
 vbe.resize(2,2*nbdrypieces);
 vbe=dbe;
 if (!(dmnor.dbelab)) dmnor.dbelab=new KN<long>;
 KN<long> &belab=*(dmnor.dbelab);
 belab.resize(2*nbdrypieces);
 for (long m=nedge;m<ndualedge;++m) {
  belab(m-nedge)=Th.bedges[bepiece(m)].lab;
 }
 if (!(dmnor.interf)) dmnor.interf=new KNM<long>;
 KNM<long> &interfa=*(dmnor.interf);
 interfa.resize(sizeneip,Th.nv);
 interfa=listneipord;
 if (!(dmnor.orientation)) dmnor.orientation=new KNM<long>;
 KNM<long> &orient=*(dmnor.orientation);
 orient.resize(sizeneip,Th.nv);
 orient=orientvertinterford;

//////////////////////////////////////////////////////////
// extra quantities
//////////////////////////////////////////////////////
 if (!(extras)) extras=new dmextras;
 dmextras &dmext=*(extras);

// array (size 2*nedge) of x,y components of point of cut on edge (written adjacently, the list is by couple)
 if ((dmext.dmcutedge)) delete dmext.dmcutedge;
 dmext.dmcutedge=new KN<double>(2*nedge);
 KN<double> &cutedges=*(dmext.dmcutedge);
 for (long m=0;m<nedge;++m){
  double centerkx0=vxy(2*cvbdry(0L,m));
  double centerky0=vxy(2*cvbdry(0L,m)+1);
  double centeradjx0=vxy(2*cvbdry(1L,m));
  double centeradjy0=vxy(2*cvbdry(1L,m)+1);
  long vertL=ctLR(2*m),vertR=ctLR(2*m+1);
  double scalinters0=((Th(vertL).x-centerkx0)*(Th(vertR).y-Th(vertL).y)-(Th(vertL).y-centerky0)*(Th(vertR).x-Th(vertL).x))/((centeradjx0-centerkx0)*(Th(vertR).y-Th(vertL).y)-(centeradjy0-centerky0)*(Th(vertR).x-Th(vertL).x));
  cutedges(2*m)=centerkx0+scalinters0*(centeradjx0-centerkx0);
  cutedges(2*m+1)=centerky0+scalinters0*(centeradjy0-centerky0);
 }
 //end of definition of cutedges

// unitary normal to interfaces of control volumes, oriented from Left vertex to Right vertex
 if ((dmext.dmnormal)) delete dmext.dmnormal;
 dmext.dmnormal=new KN<double>(2*ndualedge);
 KN<double> &normals=*(dmext.dmnormal);
 if ((dmext.dmlengthedge)) delete dmext.dmlengthedge;
 dmext.dmlengthedge=new KN<double>(ndualedge);
 KN<double> &lengthedge=*(dmext.dmlengthedge);
 for (long m=0;m<ndualedge;m++) {// components of the tangent vector on the boundary of control volume
  double tangx=vxy(2*cvbdry(1L,m))-vxy(2*cvbdry(0L,m));
  double tangy=vxy(2*cvbdry(1L,m)+1)-vxy(2*cvbdry(0L,m)+1);
  double normtang=sqrt(tangx*tangx+tangy*tangy);
  lengthedge(m)=normtang;
  normals(2*m)=tangy/normtang;
  normals(2*m+1)=-tangx/normtang;
 }
 //end of definition of normals and lengthedge

// smallest area between L and R control volumes (used for CFL condition)
 if ((dmext.dmareaminLR)) delete dmext.dmareaminLR;
 dmext.dmareaminLR=new KN<double>(ndualedge);
 KN<double> &areaminLR=*(dmext.dmareaminLR);
 for (long m=0;m<ndualedge;m++){// loop on interfaces
  areaminLR(m)=volcv(ctLR(2*m));
  if (ctLR(2*m+1)!=-1){// if interface is internal
   areaminLR(m)=min(areaminLR(m),volcv(ctLR(2*m+1)));
  }
 }
 //end of definition of areaminLR

//////////////////////////////////////////////////////////
// matrices
//////////////////////////////////////////////////////
 if (!(matrices)) matrices=new dmmatrices;
 dmmatrices &dmmat=*(matrices);

//matrix to change vertex values to gravity centers values, making conversion from FE representation to FV representation
 if ((dmmat.matvert2grav)) delete dmmat.matvert2grav;
 dmmat.matvert2grav=new MatriceMorse<double>(Th.nv,Th.nv,0,0);
 MatriceMorse<double> &Amatvert2grav=*(dmmat.matvert2grav);
 for (long i=0;i<Th.nv;i++) {
  R2 P(gravcv(2*i),gravcv(2*i+1));
  R2 PHat;
  bool outside;
//  MeshPoint mp;
//  mp.set(gravcv(2*i),gravcv(2*i+1),0.);
//  const Triangle*K = Th.Find(mp.P.p2(),PHat,outside);
  const Triangle*K=Th.Find(P,PHat,outside);
  long k=Th(*K);// long k=Th(gravcvx(i),gravcvy(i)).nuTriangle;
  double xp0=Th(Th(k,0)).x;
  double yp0=Th(Th(k,0)).y;
  double xp1=Th(Th(k,1)).x;
  double yp1=Th(Th(k,1)).y;
  double xp2=Th(Th(k,2)).x;
  double yp2=Th(Th(k,2)).y;
  double v0x=xp1-xp2;
  double v0y=yp1-yp2;
  double v1x=xp0-xp2;
  double v1y=yp0-yp2;
  double v2x=xp0-xp1;
  double v2y=yp0-yp1;
  double cg0=((gravcv(2*i)-xp1)*v0y-(gravcv(2*i+1)-yp1)*v0x)/((xp0-xp1)*v0y-(yp0-yp1)*v0x);
  double cg1=((gravcv(2*i)-xp2)*v1y-(gravcv(2*i+1)-yp2)*v1x)/((xp1-xp2)*v1y-(yp1-yp2)*v1x);
  double cg2=((gravcv(2*i)-xp0)*v2y-(gravcv(2*i+1)-yp0)*v2x)/((xp2-xp0)*v2y-(yp2-yp0)*v2x);
//  double errcg=abs(cg1-PHat.x)+abs(cg2-PHat.y);
//  if (outside) {cout << "outside, k " << k << endl; cout << "PHat " << PHat.x << " " << PHat.y << endl;}
//  if (errcg>1.e-3) {
//  cout << "errcg " << abs(errcg) << endl;
//  cout << "x,y,k " << gravcv(2*i) << " " << gravcv(2*i+1) << " " << k << endl;
//  cout << "cg0,cg1,cg2 " << cg0 << " " << cg1 << " " << cg2 << endl;
//  }
  Amatvert2grav(i,Th(k,0))=cg0;//1.-PHat.x-PHat.y;//  barycentric coordinates can be negative if the gravity center P of the control volume falls out of the mesh.
  Amatvert2grav(i,Th(k,1))=cg1;//PHat.x;//
  Amatvert2grav(i,Th(k,2))=cg2;//PHat.y;//
 }
 //end of definition of Amatvert2grav

//connectivity matrices between pieces of boundary of control volumes (=interfaces) and control volumes
 if ((dmmat.mat_piece2cvL)) delete dmmat.mat_piece2cvL;
 dmmat.mat_piece2cvL=new MatriceMorse<double>(Th.nv,ndualedge,0,0);
 MatriceMorse<double> &Amat_piece2cvL=*(dmmat.mat_piece2cvL);
 for (long m=0;m<ndualedge;m++){// loop on interfaces of control volumes
  Amat_piece2cvL(ctLR(2*m),m)=1.;
 }
 //end of definition of Amat_piece2cvL
 if ((dmmat.mat_piece2cvR)) delete dmmat.mat_piece2cvR;
 dmmat.mat_piece2cvR=new MatriceMorse<double>(Th.nv,ndualedge,0,0);
 MatriceMorse<double> &Amat_piece2cvR=*(dmmat.mat_piece2cvR);
 for (long m=0;m<nedge;m++){// loop on INTERNAL interfaces of control volumes
   Amat_piece2cvR(ctLR(2*m+1),m)=1.;
 }
 //end of definition of Amat_piece2cvR

// matrix for computing slopes in triangles in terms of values at gravity centers of control volumes
// this means dx(), dy() but with respect to gravity centers instead of vertices of Th
 if ((dmmat.mat_dxygrav)) delete dmmat.mat_dxygrav;
 dmmat.mat_dxygrav=new MatriceMorse<double>(2*Th.nt,Th.nv,0,0);
 MatriceMorse<double> &Amat_dxygrav=*(dmmat.mat_dxygrav);
 for (long k=0;k<Th.nt;k++) {
  long vert0=Th(k,0);
  long vert1=Th(k,1);
  long vert2=Th(k,2);
  double grav0x=gravcv(2*vert0),grav0y=gravcv(2*vert0+1);
  double grav1x=gravcv(2*vert1),grav1y=gravcv(2*vert1+1);
  double grav2x=gravcv(2*vert2),grav2y=gravcv(2*vert2+1);
  double vec1x=grav1x-grav0x,vec1y=grav1y-grav0y;
  double vec2x=grav2x-grav0x,vec2y=grav2y-grav0y;
  double vec1sq=vec1x*vec1x+vec1y*vec1y;
  double vec2sq=vec2x*vec2x+vec2y*vec2y;
  double vec1vec2=vec1x*vec2x+vec1y*vec2y;
  double com1=vec1sq/(vec1sq*vec2sq-vec1vec2*vec1vec2);
  double com2=vec2sq/(vec1sq*vec2sq-vec1vec2*vec1vec2);
  double com12=vec1vec2/(vec1sq*vec2sq-vec1vec2*vec1vec2);
  Amat_dxygrav(2*k,vert0)=(-com2+com12)*vec1x+(com12-com1)*vec2x;
  Amat_dxygrav(2*k,vert1)=com2*vec1x-com12*vec2x;
  Amat_dxygrav(2*k,vert2)=-com12*vec1x+com1*vec2x;
  Amat_dxygrav(2*k+1,vert0)=(-com2+com12)*vec1y+(com12-com1)*vec2y;
  Amat_dxygrav(2*k+1,vert1)=com2*vec1y-com12*vec2y;
  Amat_dxygrav(2*k+1,vert2)=-com12*vec1y+com1*vec2y;
 }
 //end of definition of Amat_dxygrav

// matrices for computing values at middle of interfaces (used for second-order reconstruction)
// defined by vector joining gravity center of control volume to the middle of the interface
 if ((dmmat.mat_piece2cvdxyL)) delete dmmat.mat_piece2cvdxyL;
 dmmat.mat_piece2cvdxyL=new MatriceMorse<double>(ndualedge,2*Th.nv,0,0);
 MatriceMorse<double> &Amat_piece2cvdxyL=*(dmmat.mat_piece2cvdxyL);
 if ((dmext.dmmidedge)) delete dmext.dmmidedge;
 dmext.dmmidedge=new KN<double>(2*ndualedge);
 KN<double> &midedge=*(dmext.dmmidedge);
 for (long m=0;m<ndualedge;m++){// loop on interfaces of control volumes
  Amat_piece2cvdxyL(m,2*ctLR(2*m))=(vxy(2*cvbdry(0L,m))+vxy(2*cvbdry(1L,m)))/2.-gravcv(2*ctLR(2*m));
  Amat_piece2cvdxyL(m,2*ctLR(2*m)+1)=(vxy(2*cvbdry(0L,m)+1)+vxy(2*cvbdry(1L,m)+1))/2.-gravcv(2*ctLR(2*m)+1);
  midedge(2*m)=(vxy(2*cvbdry(0L,m))+vxy(2*cvbdry(1L,m)))/2.;
  midedge(2*m+1)=(vxy(2*cvbdry(0L,m)+1)+vxy(2*cvbdry(1L,m)+1))/2.;
 }
 //end of definition of Amat_piece2cvdxyL
 if ((dmmat.mat_piece2cvdxyR)) delete dmmat.mat_piece2cvdxyR;
 dmmat.mat_piece2cvdxyR=new MatriceMorse<double>(ndualedge,2*Th.nv,0,0);
 MatriceMorse<double> &Amat_piece2cvdxyR=*(dmmat.mat_piece2cvdxyR);
 for (long m=0;m<nedge;m++){// loop on INTERNAL interfaces of control volumes
  Amat_piece2cvdxyR(m,2*ctLR(2*m+1))=(vxy(2*cvbdry(0L,m))+vxy(2*cvbdry(1L,m)))/2.-gravcv(2*ctLR(2*m+1));
  Amat_piece2cvdxyR(m,2*ctLR(2*m+1)+1)=(vxy(2*cvbdry(0L,m)+1)+vxy(2*cvbdry(1L,m)+1))/2.-gravcv(2*ctLR(2*m+1)+1);
 }
 //end of definition of Amat_piece2cvdxyR

// matrices for gathering the interfaces left/right values at vertices to perform limitation (used for second-order reconstruction)
 if ((dmmat.mat_piece2cvprodL)) delete dmmat.mat_piece2cvprodL;
 dmmat.mat_piece2cvprodL=new MatriceMorse<double>(sizeneip*Th.nv,ndualedge,0,0);
 MatriceMorse<double> &Amat_piece2cvprodL=*(dmmat.mat_piece2cvprodL);
 {
  KN<long> nbneipiece(Th.nv);// number of interfaces around each vertex=control volume
  nbneipiece=0L;
  for (long m=0;m<ndualedge;m++) {// loop on left interfaces
   long vertm=ctLR(2*m);
   Amat_piece2cvprodL(vertm*sizeneip+nbneipiece(vertm),m)=1.;
   nbneipiece(vertm)++;
  }
  for (long m=0;m<nedge;m++) {// loop on right interfaces. No right vertex on boundary
   long vertm=ctLR(2*m+1);
   nbneipiece(vertm)++;
  }
 }
 //end of definition of Amat_piece2cvprodL
 if ((dmmat.mat_piece2cvprodR)) delete dmmat.mat_piece2cvprodR;
 dmmat.mat_piece2cvprodR=new MatriceMorse<double>(sizeneip*Th.nv,ndualedge,0,0);
 MatriceMorse<double> &Amat_piece2cvprodR=*(dmmat.mat_piece2cvprodR);
 {
  KN<long> nbneipiece(Th.nv);// number of interfaces around each vertex=control volume
  nbneipiece=0L;
  for (long m=0;m<ndualedge;m++) {// loop on left interfaces
   long vertm=ctLR(2*m);
   nbneipiece(vertm)++;
  }
  for (long m=0;m<nedge;m++) {// loop on right interfaces. No right vertex on boundary
   long vertm=ctLR(2*m+1);
   Amat_piece2cvprodR(vertm*sizeneip+nbneipiece(vertm),m)=1.;
   nbneipiece(vertm)++;
  }
 }
 //end of definition of Amat_piece2cvprodR

// matrix for gathering the triangle slopes at vertices to perform slope limitation (used for second-order reconstruction)
 if ((dmmat.mat_tri2cvprod)) delete dmmat.mat_tri2cvprod;
 dmmat.mat_tri2cvprod=new MatriceMorse<double>(sizeneip*Th.nv,Th.nt,0,0);
 MatriceMorse<double> &Amat_tri2cvprod=*(dmmat.mat_tri2cvprod);
 if ((dmext.dmnbneit)) delete dmext.dmnbneit;
 dmext.dmnbneit=new KN<long>(Th.nv);
 KN<long> &nbneit=*(dmext.dmnbneit);
 if ((dmext.dmvolcvprod)) delete dmext.dmvolcvprod;
 dmext.dmvolcvprod=new KN<double>(sizeneip*Th.nv);
 KN<double> &volcvprod=*(dmext.dmvolcvprod);
 nbneit=0L;
 KN<long> tricvprod(sizeneip*Th.nv);// indices of triangle and its vertex
 for (long k=0;k<Th.nt;k++) {
  for (int j=0;j<3;j++) {
   Amat_tri2cvprod(Th(k,j)*sizeneip+nbneit(Th(k,j)),k)=1.;
   tricvprod(Th(k,j)*sizeneip+nbneit(Th(k,j)))=3*k+j;
   nbneit(Th(k,j))++;
  }
 }
 for (long i=0;i<Th.nv;i++) {
  for (long j=0;j<nbneit(i);j++) {
   long k=tricvprod(i*sizeneip+j)/3;
   long jv=tricvprod(i*sizeneip+j)%3;
   long e1=edgesmesh(3*k+(jv+1)%3);
   long e2=edgesmesh(3*k+(jv+2)%3);
   double ck1x=cutedges(2*e1);
   double ck1y=cutedges(2*e1+1);
   double ck2x=Th(Th(k,jv)).x;
   double ck2y=Th(Th(k,jv)).y;
   double ck3x=vxy(2*k);
   double ck3y=vxy(2*k+1);
   double areatri=abs((ck2x-ck1x)*(ck3y-ck1y)-(ck2y-ck1y)*(ck3x-ck1x))/2.;
   volcvprod(i*sizeneip+j)=areatri;
   ck1x=cutedges(2*e2);
   ck1y=cutedges(2*e2+1);
   areatri=abs((ck2x-ck1x)*(ck3y-ck1y)-(ck2y-ck1y)*(ck3x-ck1x))/2.;
   volcvprod(i*sizeneip+j)+=areatri;
  }
 }
 //end of definition of Amat_tri2cvprod and nbneit and volcvprod

// mean unitary normal to primal boundary vertices, or 0 if not on physical boundary
 if ((dmext.dmNbdryvert)) delete dmext.dmNbdryvert;
 dmext.dmNbdryvert=new KN<double>;
 KN<double> &Nbdryvert=*(dmext.dmNbdryvert);
 const long ndbe=2*nbdrypieces;
 KN<double> Nsbdryvert(4*ndbe/2);//0, 1 or 2 boundary normals, for each boundary vertex
 Nsbdryvert=0.;
 for (long ib=0;ib<ndbe/2;ib++) {// loop on couples of boundary pieces, or equivalently primal vertices on boundary, or primal true boundary element
  long m=nedge+2*ib;// index of first piece, the second is m+1
  long labelonbdry=belab(m-nedge);//label on boundary
  if (!(labelonbdry==-111111)) {//if on physical boundary
   long vert1=ctLR(2*m);// primal vertex of origin of primal boundary element
   long vert2=ctLR(2*(m+1));// primal vertex of end of primal boundary element
   double ex=Th(vert2).x-Th(vert1).x;//coordinates of edge vector
   double ey=Th(vert2).y-Th(vert1).y;
   double Nx=ey,Ny=-ex;//external normal (not unitary)
   if (Nsbdryvert(4*ib)*Nsbdryvert(4*ib)+Nsbdryvert(4*ib+1)*Nsbdryvert(4*ib+1)<1.e-30){
    Nsbdryvert(4*ib)=Nx;Nsbdryvert(4*ib+1)=Ny;
   } else if (Nsbdryvert(4*ib+2)*Nsbdryvert(4*ib+2)+Nsbdryvert(4*ib+3)*Nsbdryvert(4*ib+3)<1.e-30){
    Nsbdryvert(4*ib+2)=Nx;Nsbdryvert(4*ib+3)=Ny;
   }
   long ibnext=(vbe(1L,2*ib+1)-Th.nt)/2;
   if (Nsbdryvert(4*ibnext)*Nsbdryvert(4*ibnext)+Nsbdryvert(4*ibnext+1)*Nsbdryvert(4*ibnext+1)<1.e-30){
    Nsbdryvert(4*ibnext)=Nx;Nsbdryvert(4*ibnext+1)=Ny;
   } else if (Nsbdryvert(4*ibnext+2)*Nsbdryvert(4*ibnext+2)+Nsbdryvert(4*ibnext+3)*Nsbdryvert(4*ibnext+3)<1.e-30){
    Nsbdryvert(4*ibnext+2)=Nx;Nsbdryvert(4*ibnext+3)=Ny;
   }
  }
 }
 Nbdryvert.resize(2*ndbe/2);
 for (long ib=0;ib<ndbe/2;ib++) {//loop on primal vertices on boundary
  if (Nsbdryvert(4*ib+2)*Nsbdryvert(4*ib+2)+Nsbdryvert(4*ib+3)*Nsbdryvert(4*ib+3)>1.e-30){
   double ww=sqrt((Nsbdryvert(4*ib)+Nsbdryvert(4*ib+2))*(Nsbdryvert(4*ib)+Nsbdryvert(4*ib+2))+(Nsbdryvert(4*ib+1)+Nsbdryvert(4*ib+3))*(Nsbdryvert(4*ib+1)+Nsbdryvert(4*ib+3)));
   Nbdryvert(2*ib)=(Nsbdryvert(4*ib)+Nsbdryvert(4*ib+2))/ww;
   Nbdryvert(2*ib+1)=(Nsbdryvert(4*ib+1)+Nsbdryvert(4*ib+3))/ww;
  } else {
   double ww=sqrt(Nsbdryvert(4*ib)*Nsbdryvert(4*ib)+Nsbdryvert(4*ib+1)*Nsbdryvert(4*ib+1));
   if (ww>1.e-30){
    Nbdryvert(2*ib)=Nsbdryvert(4*ib)/ww;
    Nbdryvert(2*ib+1)=Nsbdryvert(4*ib+1)/ww;
   } else{
    Nbdryvert(2*ib)=0.;
    Nbdryvert(2*ib+1)=0.;
   }
  }
 }
 //end of definition of Nbdryvert


}// end of dualmesh::build

// ============================================================================================================
// the median dual mesh (class dualmeshN) of a meshL (1D) or of a mesh3 (3D), for the finite volumes with the unknowns
// at the vertices (P1); in 2D the dual mesh is the one above (class dualmesh)
//   1D: the control volume of a vertex goes from the middle of the segment on its left to the middle of the
//       segment on its right; one interface by segment (at its middle), one boundary piece by boundary point
//   3D: the median dual: in a tetrahedron K, the interface of the edge [a,b] is made of the two triangles
//       (middle of [a,b], centroid of a face of K containing [a,b], centroid of K); the interface of an edge is
//       the union of these triangles in the tetrahedra around the edge; a control volume gets |K|/4 from each
//       tetrahedron K around its vertex; a boundary triangle is cut in 3 pieces (one by vertex, area/3), each
//       piece is the quadrilateral (vertex, middles of the 2 edges, centroid of the triangle)
//   surface (meshS, triangles in R^3): the median dual: in a triangle K, the interface of the edge [a,b] is the
//       segment (middle of [a,b], centroid of K), in the plane of K; the interface of an edge is the union of
//       these segments in the triangles around the edge (one interface by edge, its normal the sum of l*n, n the
//       unit conormal in the plane of each triangle); a control volume gets |K|/3 from each triangle K around its
//       vertex; a boundary edge is cut in 2 pieces (one by vertex: the half edge, outward conormal in the plane
//       of its triangle, as the dual mesh of F. Bouchut in 2D)
//   the interfaces m < npe are interior (one by primal edge), npe <= m < nde are boundary pieces
// all the vectors have 3 components ([3*i+c]), also in 1D (the meshL can be in 3D)

class dualmeshN {
public:
  int dim;                  // the mesh, as FreeVol_FE::dim: 1 meshL, 3 mesh3, 4 meshS (2: the dual mesh of
                            // F. Bouchut copied in this format by FreeVol), 0: not built
  const void *primal;       // the primal mesh (MeshL, Mesh3 or MeshS)
  long nv, nt;              // vertices (= control volumes), elements of the primal mesh
  long npe, nde;            // interior interfaces (= primal edges), all the interfaces (+ boundary pieces)
  KN<long> cvLR;            // [2m],[2m+1]: left and right control volumes of the interface m (boundary piece: the vertex twice)
  KN<double> normal;        // [3m+c]: unit normal of the interface m, from left to right (boundary piece: outward)
  KN<double> area;          // [m]: measure of the interface m (1 in 1D)
  KN<double> mid;           // [3m+c]: centroid of the interface m
  KN<long> dbelab;          // [m-npe]: label of the boundary piece m
  KN<long> bvert;           // [3(m-npe)+l]: the vertices of the boundary element of the piece m (1D: the vertex)
  KN<double> bweight;       // [3(m-npe)+l]: the P1 weights of these vertices at the centroid of the piece m
  KN<double> cvvol;         // [i]: measure of the control volume i
  KN<double> cvgrav;        // [3i+c]: center of gravity of the control volume i
  KN<long> ifbeg, iflist;   // interfaces around the control volume i: iflist[ifbeg[i] .. ifbeg[i+1]-1]
  KN<long> iforient;        // [same index]: +1 if i is the left control volume of the interface, -1 if right
  KN<double> Nbdryvert;     // [3i+c]: mean unit outward normal of the physical boundary at the vertex i (0 elsewhere)
  MatriceMorse<double> *matgrad;// [3i+c] row: the gradient of a P1 function in the control volume i (the mean of
                            // its gradients in the elements around i, with the weights |K|/(dim+1))

  const MeshL *subL;        // the control volumes as a mesh of FreeFEM (label = number + 1): 1D the half segments,
  const Mesh3 *sub3;        // 3D 24 tetrahedra by tetrahedron; fespace Wstar(Thstar, P0), built when it is used
  const MeshS *subS;        // surface: 6 triangles by triangle (vertex, middle of an edge, centroid)
  const RefCounter *primalref;// the primal mesh with the reference of build (0: not built by build)

  dualmeshN() : dim(0), primal(0), nv(0), nt(0), npe(0), nde(0), matgrad(0), subL(0), sub3(0), subS(0), primalref(0) {}
  // a variable of the script (meshdualL, meshdual3) is made by FreeFEM with init() and freed with destroy(), not
  // with the constructor and the destructor: init() sets the arrays (empty), destroy() frees them
  void init() { dim=0; primal=0; nv=nt=npe=nde=0; matgrad=0; subL=0; sub3=0; subS=0; primalref=0;
    cvLR.init(); normal.init(); area.init(); mid.init(); dbelab.init(); bvert.init(); bweight.init();
    cvvol.init(); cvgrav.init(); ifbeg.init(); iflist.init(); iforient.init(); Nbdryvert.init(); }
  template<class T> static void ResetKN(KN<T> &a) { a.destroy(); a.init(); }// delete [], then empty
  void destroy() {
    ResetKN(cvLR); ResetKN(normal); ResetKN(area); ResetKN(mid); ResetKN(dbelab); ResetKN(bvert); ResetKN(bweight);
    ResetKN(cvvol); ResetKN(cvgrav); ResetKN(ifbeg); ResetKN(iflist); ResetKN(iforient); ResetKN(Nbdryvert);
    if (matgrad) { delete matgrad; matgrad=0; }
    if (subL) { subL->destroy(); subL=0; }// destroyed when it is no more used by a fespace
    if (sub3) { DualSubForget(sub3); sub3->destroy(); sub3=0; }
    if (subS) { DualSubForget(subS); subS->destroy(); subS=0; }
    DualOfPrimalForget(this);
    if (primalref) { primalref->destroy(); primalref=0; }
    dim=0; primal=0;
  }
  ~dualmeshN() { destroy(); }

  // the true boundary elements of the mesh, without an element on the other side (not the internal boundary
  // elements: the interfaces between 2 regions, the seams of a meshS made of pieces)
  template<class M> static void TrueBoundary(const M &Th, KN<long> &tb) {
    long n=0;
    tb.resize(Th.nbe);
    for (long ib=0; ib<Th.nbe; ++ib) {
      int ie; long k=Th.BoundaryElement(ib,ie);
      int j=ie;
      if (Th.ElementAdj(k,j)<0) tb[n++]=ib;
    }
    tb.resize(n);
  }

  // the interfaces around each control volume (from cvLR)
  void buildlists() {
    KN<long> cnt(nv); cnt=0L;
    for (long m=0; m<nde; ++m) { cnt[cvLR[2*m]]++; if (m<npe) cnt[cvLR[2*m+1]]++; }
    ifbeg.resize(nv+1); ifbeg[0]=0;
    for (long i=0; i<nv; ++i) ifbeg[i+1]=ifbeg[i]+cnt[i];
    iflist.resize(ifbeg[nv]); iforient.resize(ifbeg[nv]);
    cnt=0L;
    for (long m=0; m<nde; ++m) {
      long l=cvLR[2*m]; iflist[ifbeg[l]+cnt[l]]=m; iforient[ifbeg[l]+cnt[l]]=1; cnt[l]++;
      if (m<npe) { long r=cvLR[2*m+1]; iflist[ifbeg[r]+cnt[r]]=m; iforient[ifbeg[r]+cnt[r]]=-1; cnt[r]++; }
    }
  }
  // the mean unit outward normal at the vertices of the physical boundary (label != -111111: not an interface
  // of the domain decomposition)
  void buildNbdryvert() {
    Nbdryvert.resize(3*nv); Nbdryvert=0.;
    KN<double> sa(nv); sa=0.;// the boundary area around each vertex
    for (long m=npe; m<nde; ++m) {
      if (dbelab[m-npe]==-111111) continue;
      long i=cvLR[2*m];
      for (int c=0; c<3; ++c) Nbdryvert[3*i+c]+=area[m]*normal[3*m+c];
      sa[i]+=area[m];
    }
    for (long i=0; i<nv; ++i) {
      double n=sqrt(Nbdryvert[3*i]*Nbdryvert[3*i]+Nbdryvert[3*i+1]*Nbdryvert[3*i+1]+Nbdryvert[3*i+2]*Nbdryvert[3*i+2]);
      // a singular vertex (2 parts touching at a vertex or an edge) where the outward normals cancel: no mean
      // normal, the slope is not projected (as at an interior vertex)
      if (n>1.e-8*sa[i]) for (int c=0; c<3; ++c) Nbdryvert[3*i+c]/=n;
      else for (int c=0; c<3; ++c) Nbdryvert[3*i+c]=0.;
    }
  }

  // 1D: meshL
  void build(const MeshL &Th) {
    destroy(); dim=1; primal=&Th; nv=Th.nv; nt=Th.nt;
    DualOfPrimalAdd(&Th,this,1);
    Th.increment(); primalref=&Th;
    KN<long> tb; TrueBoundary(Th,tb);
    const long nb=tb.N();
    npe=nt; nde=npe+nb;
    cvLR.resize(2*nde); normal.resize(3*nde); area.resize(nde); mid.resize(3*nde);
    dbelab.resize(nde-npe); bvert.resize(3*(nde-npe)); bweight.resize(3*(nde-npe));
    cvvol.resize(nv); cvgrav.resize(3*nv); cvvol=0.; cvgrav=0.;
    for (long k=0; k<nt; ++k) {// one interface by segment: at its middle, from its vertex 0 to its vertex 1
      long a=Th(k,0), b=Th(k,1);
      R3 A=Th(a), B=Th(b), t=B-A;
      double len=t.norme(); t/=len;
      cvLR[2*k]=a; cvLR[2*k+1]=b; area[k]=1.;
      R3 M=(A+B)*0.5;
      for (int c=0; c<3; ++c) { normal[3*k+c]=t[c]; mid[3*k+c]=M[c]; }
      R3 Ga=A*0.75+B*0.25, Gb=A*0.25+B*0.75;// the centers of the two halves of the segment
      cvvol[a]+=.5*len; cvvol[b]+=.5*len;
      for (int c=0; c<3; ++c) { cvgrav[3*a+c]+=.5*len*Ga[c]; cvgrav[3*b+c]+=.5*len*Gb[c]; }
    }
    for (long q=0; q<nb; ++q) {// one boundary piece by boundary point: the point, outward normal
      long ib=tb[q];
      int ie; long k=Th.BoundaryElement(ib,ie);
      long i=Th(Th.be(ib)[0]);
      long o=(Th(k,0)==i) ? Th(k,1) : Th(k,0);// the other vertex of the segment
      R3 P=Th(i), n=P-R3(Th(o)); n/=n.norme();
      long m=npe+q;
      cvLR[2*m]=cvLR[2*m+1]=i; area[m]=1.;
      for (int c=0; c<3; ++c) { normal[3*m+c]=n[c]; mid[3*m+c]=P[c]; }
      dbelab[q]=Th.be(ib).lab;
      for (int l=0; l<3; ++l) { bvert[3*q+l]=i; bweight[3*q+l]= l==0 ? 1. : 0.; }
    }
    for (long i=0; i<nv; ++i) for (int c=0; c<3; ++c) cvgrav[3*i+c]/=cvvol[i];
    // the gradient in the control volumes: the mean of the gradients in the segments around, weights |K|/2
    matgrad=new MatriceMorse<double>(3*nv,nv,0,0);
    MatriceMorse<double> &A=*matgrad;
    for (long k=0; k<nt; ++k) {
      long v[2]={Th(k,0),Th(k,1)};
      R3 t=R3(Th(v[1]))-R3(Th(v[0]));
      double len=t.norme();
      R3 g[2]={t*(-1./(len*len)),t*(1./(len*len))};// the gradients of the P1 basis functions
      for (int l=0; l<2; ++l) { long i=v[l]; double w=.5*len/cvvol[i];
        for (int j=0; j<2; ++j) for (int c=0; c<3; ++c) if (g[j][c]!=0.) A(3*i+c,v[j])+=w*g[j][c]; }
    }
    buildlists(); buildNbdryvert();
  }

  // 3D: mesh3, the median dual
  void build(const Mesh3 &Th) {
    destroy(); dim=3; primal=&Th; nv=Th.nv; nt=Th.nt;
    DualOfPrimalAdd(&Th,this,3);
    Th.increment(); primalref=&Th;
    static const int ed[6][4]={{0,1,2,3},{0,2,1,3},{0,3,1,2},{1,2,0,3},{1,3,0,2},{2,3,0,1}};// edge (i,j), other vertices p,q
    std::map< std::pair<long,long>, long > edges;// the primal edges (min,max) -> interface number
    for (long k=0; k<nt; ++k) for (int e=0; e<6; ++e) {
      long a=Th(k,ed[e][0]), b=Th(k,ed[e][1]);
      std::pair<long,long> key(min(a,b),max(a,b));
      if (edges.find(key)==edges.end()) { long m=edges.size(); edges[key]=m; }
    }
    KN<long> tb; TrueBoundary(Th,tb);
    const long nb=tb.N();
    npe=edges.size(); nde=npe+3*nb;
    cvLR.resize(2*nde); normal.resize(3*nde); area.resize(nde); mid.resize(3*nde);
    dbelab.resize(nde-npe); bvert.resize(3*(nde-npe)); bweight.resize(3*(nde-npe));
    cvvol.resize(nv); cvgrav.resize(3*nv); cvvol=0.; cvgrav=0.;
    for (std::map< std::pair<long,long>, long >::iterator it=edges.begin(); it!=edges.end(); ++it) {
      cvLR[2*it->second]=it->first.first; cvLR[2*it->second+1]=it->first.second;// left: the smaller vertex number
    }
    KN<double> S(3*npe), C(3*npe), sa(npe); S=0.; C=0.; sa=0.;// area vectors, centroids*areas, areas
    matgrad=new MatriceMorse<double>(3*nv,nv,0,0);
    KN<double> gw(4*3*nt);// the gradients of the P1 basis functions in each tetrahedron
    for (long k=0; k<nt; ++k) {
      long v[4]; R3 P[4], G(0.,0.,0.);
      for (int l=0; l<4; ++l) { v[l]=Th(k,l); P[l]=Th(v[l]); G+=P[l]*0.25; }
      // the interfaces of the 6 edges in K
      for (int e=0; e<6; ++e) {
        int i=ed[e][0], j=ed[e][1], p=ed[e][2], q=ed[e][3];
        long a=v[i], b=v[j];
        long m=edges[std::pair<long,long>(min(a,b),max(a,b))];
        R3 LR = R3(Th(cvLR[2*m+1]))-R3(Th(cvLR[2*m]));// from left to right
        R3 M=(P[i]+P[j])*0.5, F[2]={(P[i]+P[j]+P[p])/3.,(P[i]+P[j]+P[q])/3.};
        for (int f=0; f<2; ++f) {
          R3 s=((F[f]-M)^(G-M))*0.5;
          if ((s,LR)<0.) s=-s;
          double as=s.norme();
          R3 c=(M+F[f]+G)/3.;
          for (int cc=0; cc<3; ++cc) { S[3*m+cc]+=s[cc]; C[3*m+cc]+=as*c[cc]; }
          sa[m]+=as;
        }
      }
      // the part of K in the control volume of each vertex: |K|/4, its centroid by the 6 tetrahedra
      // (G, vertex, middle of an edge, centroid of a face) on the 3 faces of K containing the vertex
      for (int i=0; i<4; ++i) {
        for (int f=0; f<3; ++f) {// the faces (i,j,l) of K containing i
          int o[3], n=0; for (int l=0; l<4; ++l) if (l!=i) o[n++]=l;
          int j=o[f], l=o[(f+1)%3];
          R3 F=(P[i]+P[j]+P[l])/3., Mij=(P[i]+P[j])*0.5, Mil=(P[i]+P[l])*0.5;
          R3 T[2][3]={{P[i],Mij,F},{P[i],F,Mil}};
          for (int s=0; s<2; ++s) {
            double vs=fabs(((T[s][0]-G),((T[s][1]-G)^(T[s][2]-G))))/6.;
            R3 cs=(G+T[s][0]+T[s][1]+T[s][2])*0.25;
            cvvol[v[i]]+=vs;
            for (int c=0; c<3; ++c) cvgrav[3*v[i]+c]+=vs*cs[c];
          }
        }
      }
      // the gradients of the P1 basis functions: the rows of the inverse of [P1-P0,P2-P0,P3-P0]
      R3 e1=P[1]-P[0], e2=P[2]-P[0], e3=P[3]-P[0];
      double det=(e1,(e2^e3));
      R3 g1=(e2^e3)/det, g2=(e3^e1)/det, g3=(e1^e2)/det, g0=-(g1+g2+g3);
      R3 g[4]={g0,g1,g2,g3};
      for (int l=0; l<4; ++l) for (int c=0; c<3; ++c) gw[12*k+3*l+c]=g[l][c];
    }
    for (long m=0; m<npe; ++m) {
      double n=sqrt(S[3*m]*S[3*m]+S[3*m+1]*S[3*m+1]+S[3*m+2]*S[3*m+2]);
      area[m]=n;
      for (int c=0; c<3; ++c) { normal[3*m+c]=S[3*m+c]/n; mid[3*m+c]=C[3*m+c]/sa[m]; }
    }
    // the boundary triangles: 3 pieces, outward normal (away from the opposite vertex of the tetrahedron)
    for (long q=0; q<nb; ++q) {
      long ib=tb[q];
      int ie; long k=Th.BoundaryElement(ib,ie);
      long w[3]; R3 Q[3];
      for (int l=0; l<3; ++l) { w[l]=Th(Th.be(ib)[l]); Q[l]=Th(w[l]); }
      R3 n=(Q[1]-Q[0])^(Q[2]-Q[0]);
      R3 Popp=Th(Th(k,ie));
      if ((n,Q[0]-Popp)<0.) n=-n;
      double A=.5*n.norme(); n/=n.norme();
      R3 F=(Q[0]+Q[1]+Q[2])/3.;
      for (int l=0; l<3; ++l) {
        long m=npe+3*q+l, b=3*q+l;
        int l1=(l+1)%3, l2=(l+2)%3;
        cvLR[2*m]=cvLR[2*m+1]=w[l]; area[m]=A/3.;
        // the centroid of the quadrilateral (Q[l], middle of [l,l1], F, middle of [l,l2]): weights 11/18, 7/36, 7/36
        R3 c=Q[l]*(11./18.)+Q[l1]*(7./36.)+Q[l2]*(7./36.);
        for (int cc=0; cc<3; ++cc) { normal[3*m+cc]=n[cc]; mid[3*m+cc]=c[cc]; }
        dbelab[m-npe]=Th.be(ib).lab;
        bvert[3*(m-npe)]=w[l]; bvert[3*(m-npe)+1]=w[l1]; bvert[3*(m-npe)+2]=w[l2];
        bweight[3*(m-npe)]=11./18.; bweight[3*(m-npe)+1]=7./36.; bweight[3*(m-npe)+2]=7./36.;
      }
    }
    for (long i=0; i<nv; ++i) for (int c=0; c<3; ++c) cvgrav[3*i+c]/=cvvol[i];
    // the gradient in the control volumes: the mean of the gradients in the tetrahedra around, weights |K|/4
    MatriceMorse<double> &A=*matgrad;
    for (long k=0; k<nt; ++k) {
      double vk=Th[k].mesure();
      for (int l=0; l<4; ++l) { long i=Th(k,l); double wl=.25*vk/cvvol[i];
        for (int j=0; j<4; ++j) for (int c=0; c<3; ++c) A(3*i+c,Th(k,j))+=wl*gw[12*k+3*j+c]; }
    }
    buildlists(); buildNbdryvert();
  }

  // surface: meshS (triangles in R^3), the median dual
  void build(const MeshS &Th) {
    destroy(); dim=4; primal=&Th; nv=Th.nv; nt=Th.nt;
    DualOfPrimalAdd(&Th,this,4);
    Th.increment(); primalref=&Th;
    static const int ed[3][2]={{1,2},{2,0},{0,1}};// the edge e of a triangle (opposite to the vertex e)
    std::map< std::pair<long,long>, long > edges;// the primal edges (min,max) -> interface number
    for (long k=0; k<nt; ++k) for (int e=0; e<3; ++e) {
      long a=Th(k,ed[e][0]), b=Th(k,ed[e][1]);
      std::pair<long,long> key(min(a,b),max(a,b));
      if (edges.find(key)==edges.end()) { long m=edges.size(); edges[key]=m; }
    }
    KN<long> tb; TrueBoundary(Th,tb);
    const long nb=tb.N();
    npe=edges.size(); nde=npe+2*nb;
    cvLR.resize(2*nde); normal.resize(3*nde); area.resize(nde); mid.resize(3*nde);
    dbelab.resize(nde-npe); bvert.resize(3*(nde-npe)); bweight.resize(3*(nde-npe));
    cvvol.resize(nv); cvgrav.resize(3*nv); cvvol=0.; cvgrav=0.;
    for (std::map< std::pair<long,long>, long >::iterator it=edges.begin(); it!=edges.end(); ++it) {
      cvLR[2*it->second]=it->first.first; cvLR[2*it->second+1]=it->first.second;// left: the smaller vertex number
    }
    KN<double> S(3*npe), C(3*npe), sl(npe); S=0.; C=0.; sl=0.;// sum of l*n, of l*centers, of l
    KN<double> gw(3*3*nt);// the gradients of the P1 basis functions in each triangle (in its plane)
    for (long k=0; k<nt; ++k) {
      long v[3]; R3 P[3];
      for (int l=0; l<3; ++l) { v[l]=Th(k,l); P[l]=Th(v[l]); }
      R3 G=(P[0]+P[1]+P[2])/3., Nk=(P[1]-P[0])^(P[2]-P[0]);
      double ak=.5*Nk.norme(); Nk/=Nk.norme();// the area and the unit normal of the triangle
      for (int l=0; l<3; ++l) {// the part of the vertex l: the quadrilateral (P_l, middles of its 2 edges, G)
        int l1=(l+1)%3, l2=(l+2)%3;
        R3 c=(P[l]*2.+(P[l]+P[l1])*.5+(P[l]+P[l2])*.5+G*2.)/6.;// the mean of the centroids of its 2 triangles
        cvvol[v[l]]+=ak/3.;
        for (int cc=0; cc<3; ++cc) cvgrav[3*v[l]+cc]+=ak/3.*c[cc];
        R3 g=(Nk^(P[l2]-P[l1]))/(2.*ak);
        for (int cc=0; cc<3; ++cc) gw[9*k+3*l+cc]=g[cc];
      }
      for (int e=0; e<3; ++e) {// the segment (middle of the edge, G) of the interface of the edge e
        long a=Th(k,ed[e][0]), b=Th(k,ed[e][1]);
        long i=min(a,b), j=max(a,b), m=edges[std::pair<long,long>(i,j)];
        R3 Pi=Th(i), Pj=Th(j), M=(Pi+Pj)*.5, t=G-M;
        double l=t.norme();
        R3 n=t^Nk; n/=n.norme();
        if ((n,Pj-Pi)<0.) n=-n;// from the left vertex i to the right one j
        R3 cs=(M+G)*.5;
        for (int cc=0; cc<3; ++cc) { S[3*m+cc]+=l*n[cc]; C[3*m+cc]+=l*cs[cc]; }
        sl[m]+=l;
      }
    }
    for (long m=0; m<npe; ++m) {
      R3 s(S[3*m],S[3*m+1],S[3*m+2]);
      area[m]=s.norme();
      for (int c=0; c<3; ++c) { normal[3*m+c]=s[c]/area[m]; mid[3*m+c]=C[3*m+c]/sl[m]; }
    }
    for (long q=0; q<nb; ++q) {// 2 boundary pieces by boundary edge: the half edges, outward conormal
      long ib=tb[q];
      int ie; long k=Th.BoundaryElement(ib,ie);
      long w[2]={Th(Th.be(ib)[0]),Th(Th.be(ib)[1])};
      R3 Q[2]={Th(w[0]),Th(w[1])}, Popp=Th(Th(k,ie));
      R3 Nk=(R3(Th(Th(k,1)))-R3(Th(Th(k,0))))^(R3(Th(Th(k,2)))-R3(Th(Th(k,0))));
      R3 t=Q[1]-Q[0], n=t^Nk;
      double len=t.norme(); n/=n.norme();
      if ((n,Q[0]-Popp)<0.) n=-n;
      for (int l=0; l<2; ++l) {
        long m=npe+2*q+l, b=3*(m-npe);
        cvLR[2*m]=cvLR[2*m+1]=w[l]; area[m]=.5*len;
        R3 c=Q[l]*.75+Q[1-l]*.25;// the middle of the half edge
        for (int cc=0; cc<3; ++cc) { normal[3*m+cc]=n[cc]; mid[3*m+cc]=c[cc]; }
        dbelab[m-npe]=Th.be(ib).lab;
        bvert[b]=w[l]; bvert[b+1]=w[1-l]; bvert[b+2]=w[l];
        bweight[b]=.75; bweight[b+1]=.25; bweight[b+2]=0.;
      }
    }
    for (long i=0; i<nv; ++i) for (int c=0; c<3; ++c) cvgrav[3*i+c]/=cvvol[i];
    // the gradient in the control volumes: the mean of the gradients in the triangles around, weights |K|/3
    matgrad=new MatriceMorse<double>(3*nv,nv,0,0);
    MatriceMorse<double> &A=*matgrad;
    for (long k=0; k<nt; ++k) {
      double ak=Th[k].mesure();
      for (int l=0; l<3; ++l) { long i=Th(k,l); double wl=ak/3./cvvol[i];
        for (int j=0; j<3; ++j) for (int c=0; c<3; ++c) A(3*i+c,Th(k,j))+=wl*gw[9*k+3*j+c]; }
    }
    buildlists(); buildNbdryvert();
  }
};

// the types of the plugin fvdualmesh: meshdualL (dual mesh of a meshL), meshdual3 (dual mesh of a mesh3),
// meshdualS (dual mesh of a meshS)
class dualmeshL : public dualmeshN {};
class dualmesh3 : public dualmeshN {};
class dualmeshS : public dualmeshN {};

#endif
