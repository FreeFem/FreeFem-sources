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

// fvdualmesh.cpp Build the dual mesh of a given mesh,
// and associated quantities needed for finite volume use
// Francois Bouchut, 2025

#include "fvdualmesh.hpp"
#include "lgmesh3.hpp"     // meshL, mesh3
namespace Fem2D { using std::max; }// the max of 2 int (hidden by the max of FESpace.hpp in PkLagrange.hpp)
#include "PkLagrange.hpp"  // TypeOfFE_Lagrange (the P0 of the dual mesh in 1D and 3D)


//-------------------------------------------------------
// transfert functions for .edp code
//-------------------------------------------------------

// the mesh of the control volumes (fespace Wstar(Thstar, P0), int2d(Thstar), ...): built when it is used; when the
// dual mesh is built again (Thstar = builddualmesh(Th) after adaptmesh), a fespace on the old one keeps the address of
// the pointer: it is built again at once, and the fespace follows the new dual mesh (as after Th = adaptmesh(Th, ...))
static pmesh *Dual_Mesh_sub(dualmesh *const &pdm);
static pmeshL *DualL_sub(dualmeshL *const &pdm);
static pmesh3 *Dual3_sub(dualmesh3 *const &pdm);
static pmeshS *DualS_sub(dualmeshS *const &pdm);
static void DualSubAgain(dualmesh *pdm) { Dual_Mesh_sub(pdm); }
static void DualSubAgain(dualmeshL *pdm) { DualL_sub(pdm); }
static void DualSubAgain(dualmesh3 *pdm) { Dual3_sub(pdm); }
static void DualSubAgain(dualmeshS *pdm) { DualS_sub(pdm); }
static bool DualHasSub(const dualmesh *pdm) { return pdm->sub != 0; }
static bool DualHasSub(const dualmeshN *pdm) { return pdm->subL || pdm->sub3 || pdm->subS; }

double Build_Dual_Mesh(dualmesh* const &pdm, const Mesh *const &pTh){
 ffassert(pTh);
 bool sub = DualHasSub(pdm);
 pdm->build(pTh);
 if (sub) DualSubAgain(pdm);
 return 0.;
}

// builddualmesh(Th): the dual mesh of Th as a value, for meshdual Thstar = builddualmesh(Th); (operator <-) and
// Thstar = builddualmesh(Th); (operator =, the previous dual mesh destroyed); only the primal mesh is kept here,
// the dual mesh is built in the meshdual variable
struct DualMeshOf { const Mesh *pTh; };
DualMeshOf Dual_Mesh_of(const Mesh *const &pTh){
 DualMeshOf d;
 d.pTh = pTh;
 return d;
}
dualmesh *Dual_Mesh_init(dualmesh *const &pdm, const DualMeshOf &d){
 ffassert(d.pTh);
 pdm->init();
 pdm->build(d.pTh);
 return pdm;
}
dualmesh *Dual_Mesh_set(dualmesh *const &pdm, const DualMeshOf &d){
 ffassert(d.pTh);
 bool sub = DualHasSub(pdm);
 pdm->destroy();
 pdm->init();
 pdm->build(d.pTh);
 if (sub) DualSubAgain(pdm);
 return pdm;
}

//-------------------------------------------------------
// the median dual mesh in 1D (meshdualL, of a meshL), in 3D (meshdual3, of a mesh3) and on a surface (meshdualS, of a
// meshS): class dualmeshN of fvdualmesh.hpp
//-------------------------------------------------------
// the classes dualmeshL, dualmesh3 and dualmeshS (types meshdualL, meshdual3, meshdualS): fvdualmesh.hpp
template<class D, class M> double Build_DualN(D* const &pdm, const M* const &pTh){
 ffassert(pTh);
 bool sub = DualHasSub(pdm);
 pdm->build(*pTh);
 if (sub) DualSubAgain(pdm);
 return 0.;
}
// builddualmesh(Th) as a value: meshdualL Thstar = builddualmesh(ThL); (operator <-), Thstar = builddualmesh(ThL); (operator =)
template<class M> struct DualMeshOfN { const M *pTh; };
template<class M> DualMeshOfN<M> Dual_MeshN_of(const M* const &pTh){
 DualMeshOfN<M> d;
 d.pTh = pTh;
 return d;
}
template<class D, class M> D *Dual_MeshN_init(D* const &pdm, const DualMeshOfN<M> &d){
 ffassert(d.pTh);
 pdm->init();
 pdm->build(*d.pTh);
 return pdm;
}
template<class D, class M> D *Dual_MeshN_set(D* const &pdm, const DualMeshOfN<M> &d){
 ffassert(d.pTh);
 bool sub = DualHasSub(pdm);
 pdm->destroy();
 pdm->init();
 pdm->build(*d.pTh);
 if (sub) DualSubAgain(pdm);
 return pdm;
}
template<class D> long DualN_ncv(D* const pdm){ if (!pdm->dim) ExecError("meshdualL/meshdual3/meshdualS: the dual mesh is not built (builddualmesh)"); return pdm->nv; }
template<class D> long DualN_npe(D* const pdm){ if (!pdm->dim) ExecError("meshdualL/meshdual3/meshdualS: the dual mesh is not built (builddualmesh)"); return pdm->npe; }
template<class D> long DualN_nde(D* const pdm){ if (!pdm->dim) ExecError("meshdualL/meshdual3/meshdualS: the dual mesh is not built (builddualmesh)"); return pdm->nde; }
template<class D> long DualN_ndbe(D* const pdm){ if (!pdm->dim) ExecError("meshdualL/meshdual3/meshdualS: the dual mesh is not built (builddualmesh)"); return pdm->nde-pdm->npe; }
template<class D> double DualN_cvarea(KN<double>* const &pvolcv, D* const &pdm){
 if (!pdm->dim) ExecError("dualmeshcvarea: the dual mesh is not built (builddualmesh)");
 pvolcv->resize(pdm->nv);
 *pvolcv=pdm->cvvol;
 return 0.;
}

long Dual_Mesh_npe(dualmesh* const pdm){//number of primal edges
 ffassert(pdm->primal);
 const KN<long> &nbedtab=*(pdm->nbed);
 return nbedtab[0];
}

long Dual_Mesh_nde(dualmesh* const pdm){//number of dual edges
 ffassert(pdm->primal);
 const KN<long> &nbedtab=*(pdm->nbed);
 return nbedtab[1];
}

long Dual_Mesh_ndv(dualmesh* const pdm){//number of dual vertices
 ffassert(pdm->primal);
// const KN<long> &nbedtab=*(pdm->nbed);
// return (*(pdm->primal)).nt+nbedtab[1]-nbedtab[0];
 const KN<long> &nvcv=*((*(pdm->normalized)).ndvcv);
 return nvcv[0];
}

long Dual_Mesh_ncv(dualmesh* const pdm){//number of dual control volumes
 ffassert(pdm->primal);
 const KN<long> &nvcv=*((*(pdm->normalized)).ndvcv);
 return nvcv[1];// or (*(pdm->primal)).nv;
}

long Dual_Mesh_ndbe(dualmesh* const pdm){//number of dual boundary elements
 ffassert(pdm->primal);
// const KN<long> &nbedtab=*(pdm->nbed);
// return nbedtab[1]-nbedtab[0];
 const KN<long> &nvcv=*((*(pdm->normalized)).ndvcv);
 return nvcv[2];
}

long Dual_Mesh_maxinterf(dualmesh* const pdm){//maximal number of neighbours around control volumes
 ffassert(pdm->primal);
 const KN<long> &nbedtab=*(pdm->nbed);
 return nbedtab[2];
}

double Dual_Mesh_dvxy(KN<double>* const &pdvxy, dualmesh* const &pdm){
 KN<double> &vxy=*pdvxy;
 long sizetab=(*((*(pdm->normalized)).dvxy)).N();
 vxy.resize(sizetab);
 vxy=*((*(pdm->normalized)).dvxy);
 return 0.;
}

double Dual_Mesh_ninterf(KN<long>* const &pninterf, dualmesh* const &pdm){
 KN<long> &ninterfa=*pninterf;
 long sizetab=(*((*(pdm->normalized)).ninterf)).N();
 ninterfa.resize(sizetab);
 ninterfa=*((*(pdm->normalized)).ninterf);
 return 0.;
}

double Dual_Mesh_dvcv(KNM<long>* const &pdvcv, dualmesh* const &pdm){
 KNM<long> &vcv=*pdvcv;
 long sizetabn=(*((*(pdm->normalized)).dvcv)).N();
 long sizetabm=(*((*(pdm->normalized)).dvcv)).M();
 vcv.resize(sizetabn,sizetabm);
 vcv=*((*(pdm->normalized)).dvcv);
 return 0.;
}

double Dual_Mesh_dvbe(KNM<long>* const &pdvbe, dualmesh* const &pdm){
 KNM<long> &vbe=*pdvbe;
 long sizetabn=(*((*(pdm->normalized)).dvbe)).N();
 long sizetabm=(*((*(pdm->normalized)).dvbe)).M();
 vbe.resize(sizetabn,sizetabm);
 vbe=*((*(pdm->normalized)).dvbe);
 return 0.;
}

double Dual_Mesh_dbelab(KN<long>* const &pdbelab, dualmesh* const &pdm){
 KN<long> &belab=*pdbelab;
 long sizetab=(*((*(pdm->normalized)).dbelab)).N();
 belab.resize(sizetab);
 belab=*((*(pdm->normalized)).dbelab);
 return 0.;
}

double Dual_Mesh_interf(KNM<long>* const &pinterf, dualmesh* const &pdm){
 KNM<long> &interfa=*pinterf;
 long sizetabn=(*((*(pdm->normalized)).interf)).N();
 long sizetabm=(*((*(pdm->normalized)).interf)).M();
 interfa.resize(sizetabn,sizetabm);
 interfa=*((*(pdm->normalized)).interf);
 return 0.;
}

double Dual_Mesh_orientation(KNM<long>* const &porientation, dualmesh* const &pdm){
 KNM<long> &orient=*porientation;
 long sizetabn=(*((*(pdm->normalized)).orientation)).N();
 long sizetabm=(*((*(pdm->normalized)).orientation)).M();
 orient.resize(sizetabn,sizetabm);
 orient=*((*(pdm->normalized)).orientation);
 return 0.;
}

double Dual_Mesh_edgenumbering(KN<long>* const &pedgenumbering, dualmesh* const &pdm){
 KN<long> &edgesmesh=*pedgenumbering;
 long sizetab=(*(pdm->edgenumbering)).N();
 edgesmesh.resize(sizetab);
 edgesmesh=*(pdm->edgenumbering);
 return 0.;
}

double Dual_Mesh_cvLR(KN<long>* const &pcvLR, dualmesh* const &pdm){
 KN<long> &ctLR=*pcvLR;
 long sizetab=(*(pdm->cvLR)).N();
 ctLR.resize(sizetab);
 ctLR=*(pdm->cvLR);
 return 0.;
}

double Dual_Mesh_dualpiece01(KNM<long>* const &pcvbdry, dualmesh* const &pdm){
 KNM<long> &cvbdry=*pcvbdry;
 long sizetabn=(*(pdm->dualpiece01)).N();
 long sizetabm=(*(pdm->dualpiece01)).M();
 cvbdry.resize(sizetabn,sizetabm);
 cvbdry=*(pdm->dualpiece01);
 return 0.;
}

double Dual_Mesh_cvarea(KN<double>* const &pvolcv, dualmesh* const &pdm){
 KN<double> &volcv=*pvolcv;
 long sizetab=(*(pdm->cvarea)).N();
 volcv.resize(sizetab);
 volcv=*(pdm->cvarea);
 return 0.;
}

double Dual_Mesh_cvgrav(KN<double>* const &pgravcv, dualmesh* const &pdm){
 KN<double> &gravcv=*pgravcv;
 long sizetab=(*(pdm->cvgrav)).N();
 gravcv.resize(sizetab);
 gravcv=*(pdm->cvgrav);
 return 0.;
}

// array (size 2*nedge) of x,y components of point of cut on edge (written adjacently, the list is by couple)
double Dual_Mesh_cutedge(KN<double>* const &pcutedges, dualmesh* const &pdm){
 ffassert(pdm->primal);
 const long nedge=(*(pdm->nbed))[0];
 KN<double> &cutedges=*pcutedges;
 cutedges.resize(2*nedge);
 cutedges=*((pdm->extras)->dmcutedge);
 return 0.;
}//end of Dual_Mesh_cutedge

// computation of unitary normal to interfaces of control volumes, oriented from Left vertex to Right vertex
double Dual_Mesh_normals(KN<double>* const &pnormals, KN<double>* const &plengthedge, dualmesh* const &pdm){
 ffassert(pdm->primal);
 const KN<long> &nbedtab=*(pdm->nbed);
 const long ndualedge=nbedtab[1];
 KN<double> &normals=*pnormals;// unitary normal at boundary of control volume, oriented from Left vertex to Right vertex
 normals.resize(2*ndualedge);
 KN<double> &lengthedge=*plengthedge;// length of interface of control volume
 lengthedge.resize(ndualedge);
 normals=*((pdm->extras)->dmnormal);
 lengthedge=*((pdm->extras)->dmlengthedge);
 return 0.;
}// end of Dual_Mesh_normals

// smallest area between L and R control volumes (used for CFL condition)
double Dual_Mesh_areaminLR(KN<double>* const &pareaminLR, dualmesh* const &pdm){
 ffassert(pdm->primal);
 const KN<long> &nbedtab=*(pdm->nbed);
 const long ndualedge=nbedtab[1];
 KN<double> &areaminLR=*pareaminLR;
 areaminLR.resize(ndualedge);
 areaminLR=*((pdm->extras)->dmareaminLR);
 return 0.;
}// end of Dual_Mesh_areaminLR

// write dualmesh to a file. Only sequential use .vtk format for paraview with "Surface With Edges" type of plot
double Dual_Mesh_save(dualmesh* const &pdm, string* const &pfilename){
 const string &filename=*pfilename;
 ffassert(pdm->primal);
 const KN<long> &nvcv=*((*(pdm->normalized)).ndvcv);
 const KN<double> &vxy=*((*(pdm->normalized)).dvxy);
 const KN<long> &ninterfa=*((*(pdm->normalized)).ninterf);
 const KNM<long> &vcv=*((*(pdm->normalized)).dvcv);
 const KNM<long> &vbe=*((*(pdm->normalized)).dvbe);
 const KN<long> &belab=*((*(pdm->normalized)).dbelab);
 long ndv=nvcv[0],ncv=nvcv[1],ndbe=nvcv[2];
 {
  const char*pch = filename.c_str();
  int lll = strlen(pch);
  if (strcmp(pch + lll - 4, ".vtk")) {
   ExecError(" savedualmesh: only .vtk file format is supported");
   return 1;
  }
  ofstream f1(pch);
  if (!f1) {
   cout << " Error Opening file " << filename << endl;
   ExecError(" Error Opening file");
   return 1;
  }
  if (verbosity>1) cout << " writing dualmesh to " << filename << endl;
  f1.precision(15);
  f1 << "# vtk DataFile Version 3.0" << endl;
  f1 << "2D polygonal mesh with boundaries and labels" << endl;
  f1 << "ASCII" << endl;
  f1 << "DATASET UNSTRUCTURED_GRID" << endl;
  f1 << endl;

  f1 << "POINTS " << ndv << " float" << endl;// number of dual vertices
  for (long i = 0; i < ndv; ++i) {// loop on dual vertices
    f1 << vxy(2*i) << " " << vxy(2*i+1) << " " << "0." << endl;// coordinates of dual vertex
  }
  f1 << endl;

  long mysum=0;
  for (long i=0; i < ncv; ++i) {
   mysum+=ninterfa(i)+1;
  }
  mysum+=3*ndbe;
  f1 << "CELLS " << ncv+ndbe << " " << mysum << endl;// number of control volumes + dual boundary elements, and number of data
  for (long i=0; i < ncv; ++i) {// loop on control volumes
   f1 << ninterfa(i) << " ";// number of dual vertices around control volume
   for (long j=0; j<ninterfa(i); ++j) {// loop on dual vertices around control volume
    f1 << " " << vcv(j,i);// index of dual vertex
   }
   f1 << endl;
  }
  for (long idbe=0;idbe<ndbe;idbe++) {// loop on dual boundary elements
   f1 << 2 << "  " << vbe(0L,idbe) << " " << vbe(1L,idbe) << endl;// indices of the 2 dual vertices origin and end of dual boundary element
  }
  f1 << endl;

  f1 << "CELL_TYPES " << ncv+ndbe << endl;
  for (long i=0; i < ncv; ++i) {// loop on control volumes
   f1 << "7" << endl;// polygon
  }
  for (long idbe=0;idbe<ndbe;idbe++) {// loop on dual boundary elements
   f1 << "3" << endl;// line
  }
  f1 << endl;

  f1 << "CELL_DATA " << ncv+ndbe << endl;
  f1 << "SCALARS bord_label int" << endl;
  f1 << "LOOKUP_TABLE default" << endl;
  for (long i=0; i < ncv; ++i) {// loop on control volumes
   f1 << "0" << endl;// region set to 0
  }
  for (long idbe=0;idbe<ndbe;idbe++) {// loop on dual boundary elements
   f1 << belab[idbe] << endl;// boundary label
  }
 };//end of writing to file
 return 0.;
}// end of Dual_Mesh_save

//matrix to change vertex values to gravity centers values, making conversion from FE representation to FV representation
newpMatrice_Creuse<double> Mat_vert2grav(Stack stack,dualmesh* const &pdm) {
 ffassert(pdm->primal);
 const long nv=(pdm->primal)->nv;
 MatriceMorse<double> &A = *new MatriceMorse<double>(nv,nv,0,0);
 A=*((pdm->matrices)->matvert2grav);
 return newpMatrice_Creuse<double>(stack,&A);
}// end of Mat_vert2grav

//connectivity matrices between pieces of boundary of control volumes (=interfaces) and control volumes
newpMatrice_Creuse<double> Mat_piece2cvL(Stack stack,dualmesh* const &pdm) {
 ffassert(pdm->primal);
 const long nv=(pdm->primal)->nv;
 const long ndualedge=(*(pdm->nbed))[1];
 MatriceMorse<double> &A = *new MatriceMorse<double>(nv,ndualedge,0,0);
 A=*((pdm->matrices)->mat_piece2cvL);
 return newpMatrice_Creuse<double>(stack,&A);
}// end of Mat_piece2cvL
newpMatrice_Creuse<double> Mat_piece2cvR(Stack stack,dualmesh* const &pdm) {
 ffassert(pdm->primal);
 const long nv=(pdm->primal)->nv;
 const long ndualedge=(*(pdm->nbed))[1];
 MatriceMorse<double> &A = *new MatriceMorse<double>(nv,ndualedge,0,0);
 A=*((pdm->matrices)->mat_piece2cvR);
 return newpMatrice_Creuse<double>(stack,&A);
}// end of Mat_piece2cvR

// matrix for computing slopes in triangles in terms of values at gravity centers of control volumes
// this means dx(), dy() but with respect to gravity centers instead of vertices of Th
newpMatrice_Creuse<double> Mat_dxygrav(Stack stack,dualmesh* const &pdm) {
 ffassert(pdm->primal);
 const long nv=(pdm->primal)->nv;
 const long nt=(pdm->primal)->nt;
 MatriceMorse<double> &A = *new MatriceMorse<double>(2*nt,nv,0,0);
 A=*((pdm->matrices)->mat_dxygrav);
 return newpMatrice_Creuse<double>(stack,&A);
}// end of Mat_dxygrav

// matrices for computing values at middle of interfaces (used for second-order reconstruction)
// defined by vector joining gravity center of control volume to the middle of the interface
newpMatrice_Creuse<double> Mat_piece2cvdxyL(Stack stack,dualmesh* const &pdm) {
 ffassert(pdm->primal);
 const long nv=(pdm->primal)->nv;
 const long ndualedge=(*(pdm->nbed))[1];
 MatriceMorse<double> &A = *new MatriceMorse<double>(ndualedge,2*nv,0,0);
 A=*((pdm->matrices)->mat_piece2cvdxyL);
 return newpMatrice_Creuse<double>(stack,&A);
}// end of Mat_piece2cvdxyL
newpMatrice_Creuse<double> Mat_piece2cvdxyR(Stack stack,dualmesh* const &pdm) {
 ffassert(pdm->primal);
 const long nv=(pdm->primal)->nv;
 const long ndualedge=(*(pdm->nbed))[1];
 MatriceMorse<double> &A = *new MatriceMorse<double>(ndualedge,2*nv,0,0);
 A=*((pdm->matrices)->mat_piece2cvdxyR);
 return newpMatrice_Creuse<double>(stack,&A);
}// end of Mat_piece2cvdxyR

// matrices for gathering the interfaces left/right values at vertices to perform limitation (used for second-order reconstruction)
newpMatrice_Creuse<double> Mat_piece2cvprodL(Stack stack,dualmesh* const &pdm) {
 ffassert(pdm->primal);
 const long nv=(pdm->primal)->nv;
 const long ndualedge=(*(pdm->nbed))[1];
 const long sizeneip=(*(pdm->nbed))[2];
 MatriceMorse<double> &A = *new MatriceMorse<double>(sizeneip*nv,ndualedge,0,0);
 A=*((pdm->matrices)->mat_piece2cvprodL);
 return newpMatrice_Creuse<double>(stack,&A);
}// end of Mat_piece2cvprodL
newpMatrice_Creuse<double> Mat_piece2cvprodR(Stack stack,dualmesh* const &pdm) {
 ffassert(pdm->primal);
 const long nv=(pdm->primal)->nv;
 const long ndualedge=(*(pdm->nbed))[1];
 const long sizeneip=(*(pdm->nbed))[2];
 MatriceMorse<double> &A = *new MatriceMorse<double>(sizeneip*nv,ndualedge,0,0);
 A=*((pdm->matrices)->mat_piece2cvprodR);
 return newpMatrice_Creuse<double>(stack,&A);
}// end of Mat_piece2cvprodR

// matrix for gathering the triangle slopes at vertices to perform slope limitation (used for second-order reconstruction)
newpMatrice_Creuse<double> Mat_tri2cvprod(Stack stack, KN<long>* const &pnbneit, KN<double>* const &pvolcvprod, dualmesh* const &pdm) {
 ffassert(pdm->primal);
 const long nv=(pdm->primal)->nv;
 const long nt=(pdm->primal)->nt;
 const long sizeneip=(*(pdm->nbed))[2];
 KN<long> &nbneit=*pnbneit;// number of triangles around each vertex=control volume
 nbneit.resize(nv);
 nbneit=*((pdm->extras)->dmnbneit);
 KN<double> &volcvprod=*pvolcvprod;// part of area of triangle which lies in control volume
 volcvprod.resize(nv*sizeneip);
 volcvprod=*((pdm->extras)->dmvolcvprod);
 MatriceMorse<double> &A = *new MatriceMorse<double>(nv*sizeneip,nt,0,0);
 A=*((pdm->matrices)->mat_tri2cvprod);
 return newpMatrice_Creuse<double>(stack,&A);
}// end of Mat_tri2cvprod

// computation of mean unitary normal to primal boundary vertices, or 0 if not on physical boundary
double Dual_Mesh_Nbdryvert(KN<double>* const &pNbdryvert, dualmesh* const &pdm){
 ffassert(pdm->primal);
 const long nedge=(*(pdm->nbed))[0];
 const long ndualedge=(*(pdm->nbed))[1];
 const long ndbe=ndualedge-nedge;
 KN<double> &Nbdryvert=*pNbdryvert;
 const KN<double> &dm_Nbdryvert=*((pdm->extras)->dmNbdryvert);
 Nbdryvert.resize(dm_Nbdryvert.N());
 Nbdryvert=dm_Nbdryvert;
 return 0.;
}// end of Dual_Mesh_Nbdryvert

//array (size 2*ndualedge) of x,y components of middle of dual edge (written adjacently, the list is by couple)
double Dual_Mesh_Midedge(KN<double>* const &pmidedge, dualmesh* const &pdm){
 ffassert(pdm->primal);
 KN<double> &midedge=*pmidedge;
 const KN<double> &dm_midedge=*((pdm->extras)->dmmidedge);
 midedge.resize(dm_midedge.N());
 midedge=dm_midedge;
 return 0.;
}// end of Dual_Mesh_Midedge

//-------------------------------------------------------
// plot of a dual mesh with ffglut: plot(Thstar) and plot(Thstar, u[]) (one value by control volume)
// sent to ffglut as a new item of the plot: what = 60 (the mesh) or 61 (with the values)
//-------------------------------------------------------
#include "PlotStream.hpp"
extern FILE *ThePlotStream;
extern bool NoWait;

class PlotDualMesh : public E_F0mps {
 public:
  typedef long Result;
  // the items of the plot: a meshdual (e) with its values (eu, or 0), or a mesh of FreeFEM (e)
  struct Item {
    bool dual;
    Expression e, eu;
  };
  vector< Item > items;
  // same names and numbers as the options of plot (see Plot::name_param in lgfem.cpp)
  static const int n_name_param = 21;
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  PlotDualMesh(const basicAC_F0 &args) {
    args.SetNameParam(n_name_param, name_param, nargs);
    if (nargs[8] || nargs[10] || nargs[12] || nargs[15] || nargs[18] || nargs[20])
      CompileError("plot of a meshdual: options bb, nbarrow, varrow, hsv, add, ech not available");
    for (int i = 0; i < args.size( ); ++i) {
      Item it = {false, 0, 0};
      if (BCastTo< dualmesh * >(args[i])) {
        it.dual = true;
        it.e = to< dualmesh * >(args[i]);
        if (i + 1 < args.size( ) && BCastTo< KN_< double > >(args[i + 1]))
          it.eu = to< KN_< double > >(args[++i]);
      } else if (BCastTo< pmesh >(args[i]))
        it.e = to< pmesh >(args[i]);
      else
        CompileError("plot with a meshdual: only meshdual (with its values) and mesh");
      items.push_back(it);
    }
  }
  AnyType operator( )(Stack s) const;
};

// plot(Thstar), plot(Thstar, u[]), plot(Th, Thstar), plot(Thstar, Th), plot(Th, Thstar, u[]), plot(Thstar, u[], Th)
template< int N >
class PlotDualMeshN : public PlotDualMesh {
 public:
  PlotDualMeshN(const basicAC_F0 &args) : PlotDualMesh(args) {}
  static ArrayOfaType typeargs( ) {
    aType D = atype< dualmesh * >( ), U = atype< KN_< double > >( ), M = atype< pmesh >( );
    switch (N) {
      case 1: return ArrayOfaType(D);
      case 2: return ArrayOfaType(D, U);
      case 3: return ArrayOfaType(M, D);
      case 4: return ArrayOfaType(D, M);
      case 5: return ArrayOfaType(M, D, U);
      default: return ArrayOfaType(D, U, M);
    }
  }
  static E_F0 *f(const basicAC_F0 &args) { return new PlotDualMeshN(args); }
};

basicAC_F0::name_and_type PlotDualMesh::name_param[] = {
  {"coef", &typeid(double)},       {"cmm", &typeid(string *)},      {"ps", &typeid(string *)},
  {"wait", &typeid(bool)},         {"fill", &typeid(bool)},         {"value", &typeid(bool)},
  {"clean", &typeid(bool)},        {"aspectratio", &typeid(bool)},  {"bb", &typeid(E_Array)},
  {"nbiso", &typeid(long)},        {"nbarrow", &typeid(long)},      {"viso", &typeid(KN_< double >)},
  {"varrow", &typeid(KN_< double >)}, {"bw", &typeid(bool)},        {"grey", &typeid(bool)},
  {"hsv", &typeid(KN_< double >)}, {"boundary", &typeid(bool)},     {"dim", &typeid(long)},
  {"add", &typeid(bool)},          {"prev", &typeid(bool)},         {"ech", &typeid(double)}};

// the data of the item of a dual mesh for ffglut (what = 60 or 61, see OnePlotMeshD in ffglut.cpp):
// dual vertices, centers (primal vertices), polygons, border, values
static void SendDualMesh(PlotStream &theplot, dualmesh *pdm, const KN< double > &u) {
  const Mesh &Th = *pdm->primal;
  const dmnormalized &dn = *pdm->normalized;
  const KN< long > &nvcv = *dn.ndvcv;
  long ndv = nvcv[0], ncv = nvcv[1], ndbe = nvcv[2];
  const KN< long > &ninterf = *dn.ninterf;
  const KNM< long > &dvcv = *dn.dvcv;
  const KNM< long > &dvbe = *dn.dvbe;
  const KN< long > &dbelab = *dn.dbelab;
  KN< double > xy(*dn.dvxy), cxy(2 * ncv), nv(ncv), be(2 * ndbe), belab(ndbe);
  long ntot = 0;
  for (long i = 0; i < ncv; ++i) ntot += ninterf[i];
  KN< double > cv(ntot);
  for (long i = 0, k = 0; i < ncv; ++i) {
    cxy[2 * i] = Th(i).x;
    cxy[2 * i + 1] = Th(i).y;
    nv[i] = ninterf[i];
    for (long j = 0; j < ninterf[i]; ++j) cv[k++] = dvcv(j, i);
  }
  for (long i = 0; i < ndbe; ++i) {
    be[2 * i] = dvbe(0L, i);
    be[2 * i + 1] = dvbe(1L, i);
    belab[i] = dbelab[i];
  }
  ffassert(xy.N( ) == 2 * ndv);
  theplot << (long)(u.N( ) ? 61 : 60);
  theplot << xy << cxy << nv << cv << be << belab << u;
}
static void SendDualPlotArgs(PlotStream &theplot, const Expression *nargs, Stack s);

AnyType PlotDualMesh::operator( )(Stack s) const {
  // evaluate the items: the meshdual, their values, the meshes (numbered from 1, as in lgfem.cpp)
  size_t n = items.size( );
  vector< dualmesh * > pdm(n, (dualmesh *)0);
  vector< KN< double > > u(n);
  vector< long > imesh(n, 0);
  map< const Mesh *, long > mapth;
  for (size_t i = 0; i < n; ++i)
    if (items[i].dual) {
      pdm[i] = GetAny< dualmesh * >((*items[i].e)(s));
      if (!pdm[i] || !pdm[i]->primal) ExecError("plot: the meshdual is not built (builddualmesh)");
      if (items[i].eu) {
        KN_< double > uu = GetAny< KN_< double > >((*items[i].eu)(s));
        long ncv = (*pdm[i]->normalized->ndvcv)[1];
        if (uu.N( ) != ncv)
          ExecError("plot(meshdual, u): u must have one value by control volume (Th.nv)");
        u[i].resize(ncv);
        u[i] = uu;
      }
    } else {
      const Mesh *pTh = GetAny< pmesh >((*items[i].e)(s));
      if (!pTh) ExecError("plot: the mesh is not defined");
      if (mapth.find(pTh) == mapth.end( )) {
        long k = (long)mapth.size( ) + 1;
        mapth[pTh] = k;
      }
      imesh[i] = mapth[pTh];
    }
  if (!ThePlotStream) return 0L;

  PlotStream theplot(ThePlotStream);
  SendDualPlotArgs(theplot, nargs, s);
  // the meshes of FreeFEM (the dual meshes are sent in their item)
  theplot.SendMeshes( );
  theplot << (long)mapth.size( );
  for (map< const Mesh *, long >::const_iterator i = mapth.begin( ); i != mapth.end( ); ++i)
    theplot << i->second << *i->first;
  theplot.SendPlots( );
  theplot << (long)n;
  // first the dual meshes in colors: the meshes are drawn over them (plot(Th, Thstar, u[]))
  for (int pass = 0; pass < 2; ++pass)
    for (size_t i = 0; i < n; ++i)
      if ((pdm[i] && u[i].N( )) != (pass == 0))
        continue;
      else if (pdm[i])
        SendDualMesh(theplot, pdm[i], u[i]);
      else
        theplot << 0L << imesh[i];    // what = 0: a mesh, its number
  theplot.SendEndPlot( );
  return 0L;
}

// the options of the plot (same numbers as Plot::name_param in lgfem.cpp), common to meshdual and meshdual3
static void SendDualPlotArgs(PlotStream &theplot, const Expression *nargs, Stack s) {
  theplot.SendNewPlot( );
  if (nargs[0]) (theplot << 0L) <= GetAny< double >((*nargs[0])(s));
  if (nargs[1]) (theplot << 1L) <= GetAny< string * >((*nargs[1])(s));
  if (nargs[2]) (theplot << 2L) <= GetAny< string * >((*nargs[2])(s));
  (theplot << 3L) <= (bool)(!NoWait && nargs[3] && GetAny< bool >((*nargs[3])(s)));
  for (long k = 4; k <= 7; ++k)
    if (nargs[k]) (theplot << k) <= GetAny< bool >((*nargs[k])(s));
  if (nargs[9]) (theplot << 9L) <= GetAny< long >((*nargs[9])(s));
  if (nargs[11]) (theplot << 11L) <= GetAny< KN_< double > >((*nargs[11])(s));
  if (nargs[13]) (theplot << 13L) <= GetAny< bool >((*nargs[13])(s));
  if (nargs[14]) (theplot << 14L) <= GetAny< bool >((*nargs[14])(s));
  if (nargs[16]) (theplot << 16L) <= GetAny< bool >((*nargs[16])(s));
  if (nargs[17]) (theplot << 17L) <= GetAny< long >((*nargs[17])(s));
  if (nargs[19]) (theplot << 19L) <= GetAny< bool >((*nargs[19])(s));
  theplot.SendEndArgPlot( );
}

//-------------------------------------------------------
// plot of a 3d dual mesh (meshdual3, the median dual of a mesh3, class dualmeshN) with ffglut:
// plot(Thstar), plot(Thstar, u[]) (one value by control volume = by vertex of the mesh3), with the mesh3:
// plot(Th3, Thstar), plot(Thstar, Th3), plot(Th3, Thstar, u[]), plot(Thstar, u[], Th3)
// sent to ffglut as a new item: what = 62 (the mesh) or 63 (with the values), the number of the primal
// mesh3 (sent with the meshes of the plot) and the values; ffglut builds the control volumes from the mesh3
//-------------------------------------------------------
class PlotDualMesh3 : public E_F0mps {
 public:
  typedef long Result;
  struct Item {
    bool dual;
    Expression e, eu;
  };
  vector< Item > items;
  Expression nargs[PlotDualMesh::n_name_param];
  PlotDualMesh3(const basicAC_F0 &args) {
    args.SetNameParam(PlotDualMesh::n_name_param, PlotDualMesh::name_param, nargs);
    if (nargs[8] || nargs[10] || nargs[12] || nargs[15] || nargs[18] || nargs[20])
      CompileError("plot of a meshdual3: options bb, nbarrow, varrow, hsv, add, ech not available");
    for (int i = 0; i < args.size( ); ++i) {
      Item it = {false, 0, 0};
      if (BCastTo< dualmesh3 * >(args[i])) {
        it.dual = true;
        it.e = to< dualmesh3 * >(args[i]);
        if (i + 1 < args.size( ) && BCastTo< KN_< double > >(args[i + 1]))
          it.eu = to< KN_< double > >(args[++i]);
      } else if (BCastTo< pmesh3 >(args[i]))
        it.e = to< pmesh3 >(args[i]);
      else
        CompileError("plot with a meshdual3: only meshdual3 (with its values) and mesh3");
      items.push_back(it);
    }
  }
  AnyType operator( )(Stack s) const;
};

template< int N >
class PlotDualMesh3N : public PlotDualMesh3 {
 public:
  PlotDualMesh3N(const basicAC_F0 &args) : PlotDualMesh3(args) {}
  static ArrayOfaType typeargs( ) {
    aType D = atype< dualmesh3 * >( ), U = atype< KN_< double > >( ), M = atype< pmesh3 >( );
    switch (N) {
      case 1: return ArrayOfaType(D);
      case 2: return ArrayOfaType(D, U);
      case 3: return ArrayOfaType(M, D);
      case 4: return ArrayOfaType(D, M);
      case 5: return ArrayOfaType(M, D, U);
      default: return ArrayOfaType(D, U, M);
    }
  }
  static E_F0 *f(const basicAC_F0 &args) { return new PlotDualMesh3N(args); }
};

AnyType PlotDualMesh3::operator( )(Stack s) const {
  // the mesh3 of the items (the primal meshes of the meshdual3 and the mesh3 of the plot), numbered from 1
  size_t n = items.size( );
  vector< long > imesh(n, 0), what(n, 5);
  vector< KN< double > > u(n);
  map< const Mesh3 *, long > mapth;
  for (size_t i = 0; i < n; ++i) {
    const Mesh3 *pTh = 0;
    if (items[i].dual) {
      dualmesh3 *pdm = GetAny< dualmesh3 * >((*items[i].e)(s));
      if (!pdm || pdm->dim != 3 || !pdm->primal) ExecError("plot: the meshdual3 is not built (builddualmesh)");
      pTh = (const Mesh3 *)pdm->primal;
      what[i] = 62;
      if (items[i].eu) {
        KN_< double > uu = GetAny< KN_< double > >((*items[i].eu)(s));
        if (uu.N( ) != pdm->nv)
          ExecError("plot(meshdual3, u): u must have one value by control volume (Th.nv)");
        u[i].resize(pdm->nv);
        u[i] = uu;
        what[i] = 63;
      }
    } else {
      pTh = GetAny< pmesh3 >((*items[i].e)(s));
      if (!pTh) ExecError("plot: the mesh3 is not defined");
    }
    if (mapth.find(pTh) == mapth.end( )) {
      long k = (long)mapth.size( ) + 1;
      mapth[pTh] = k;
    }
    imesh[i] = mapth[pTh];
  }
  if (!ThePlotStream) return 0L;

  PlotStream theplot(ThePlotStream);
  SendDualPlotArgs(theplot, nargs, s);
  theplot.SendMeshes( );
  theplot << 0L;
  theplot.SendMeshes3( );
  theplot << (long)mapth.size( );
  for (map< const Mesh3 *, long >::const_iterator i = mapth.begin( ); i != mapth.end( ); ++i)
    theplot << i->second << *i->first;
  theplot.SendPlots( );
  theplot << (long)n;
  for (size_t i = 0; i < n; ++i) {
    theplot << what[i] << imesh[i];    // what = 5: a mesh3, 62, 63: a meshdual3 (its primal mesh3)
    if (what[i] == 63) theplot << u[i];
  }
  theplot.SendEndPlot( );
  return 0L;
}

//-------------------------------------------------------
// plot of a surface dual mesh (meshdualS, the median dual of a meshS, class dualmeshN) with ffglut:
// plot(Thstar), plot(Thstar, u[]) (one value by control volume = by vertex of the meshS), with the meshS:
// plot(ThS, Thstar), plot(Thstar, ThS), plot(ThS, Thstar, u[]), plot(Thstar, u[], ThS)
// sent to ffglut as a new item: what = 64 (the mesh) or 65 (with the values), the number of the primal
// meshS (sent with the meshes of the plot) and the values; ffglut builds the control volumes from the meshS
//-------------------------------------------------------
class PlotDualMeshS : public E_F0mps {
 public:
  typedef long Result;
  struct Item {
    bool dual;
    Expression e, eu;
  };
  vector< Item > items;
  Expression nargs[PlotDualMesh::n_name_param];
  PlotDualMeshS(const basicAC_F0 &args) {
    args.SetNameParam(PlotDualMesh::n_name_param, PlotDualMesh::name_param, nargs);
    if (nargs[8] || nargs[10] || nargs[12] || nargs[15] || nargs[18] || nargs[20])
      CompileError("plot of a meshdualS: options bb, nbarrow, varrow, hsv, add, ech not available");
    for (int i = 0; i < args.size( ); ++i) {
      Item it = {false, 0, 0};
      if (BCastTo< dualmeshS * >(args[i])) {
        it.dual = true;
        it.e = to< dualmeshS * >(args[i]);
        if (i + 1 < args.size( ) && BCastTo< KN_< double > >(args[i + 1]))
          it.eu = to< KN_< double > >(args[++i]);
      } else if (BCastTo< pmeshS >(args[i]))
        it.e = to< pmeshS >(args[i]);
      else
        CompileError("plot with a meshdualS: only meshdualS (with its values) and meshS");
      items.push_back(it);
    }
  }
  AnyType operator( )(Stack s) const;
};

template< int N >
class PlotDualMeshSN : public PlotDualMeshS {
 public:
  PlotDualMeshSN(const basicAC_F0 &args) : PlotDualMeshS(args) {}
  static ArrayOfaType typeargs( ) {
    aType D = atype< dualmeshS * >( ), U = atype< KN_< double > >( ), M = atype< pmeshS >( );
    switch (N) {
      case 1: return ArrayOfaType(D);
      case 2: return ArrayOfaType(D, U);
      case 3: return ArrayOfaType(M, D);
      case 4: return ArrayOfaType(D, M);
      case 5: return ArrayOfaType(M, D, U);
      default: return ArrayOfaType(D, U, M);
    }
  }
  static E_F0 *f(const basicAC_F0 &args) { return new PlotDualMeshSN(args); }
};

AnyType PlotDualMeshS::operator( )(Stack s) const {
  // the meshS of the items (the primal meshes of the meshdualS and the meshS of the plot), numbered from 1
  size_t n = items.size( );
  vector< long > imesh(n, 0), what(n, 50);
  vector< KN< double > > u(n);
  map< const MeshS *, long > mapth;
  for (size_t i = 0; i < n; ++i) {
    const MeshS *pTh = 0;
    if (items[i].dual) {
      dualmeshS *pdm = GetAny< dualmeshS * >((*items[i].e)(s));
      if (!pdm || pdm->dim != 4 || !pdm->primal) ExecError("plot: the meshdualS is not built (builddualmesh)");
      pTh = (const MeshS *)pdm->primal;
      what[i] = 64;
      if (items[i].eu) {
        KN_< double > uu = GetAny< KN_< double > >((*items[i].eu)(s));
        if (uu.N( ) != pdm->nv)
          ExecError("plot(meshdualS, u): u must have one value by control volume (Th.nv)");
        u[i].resize(pdm->nv);
        u[i] = uu;
        what[i] = 65;
      }
    } else {
      pTh = GetAny< pmeshS >((*items[i].e)(s));
      if (!pTh) ExecError("plot: the meshS is not defined");
    }
    if (mapth.find(pTh) == mapth.end( )) {
      long k = (long)mapth.size( ) + 1;
      mapth[pTh] = k;
    }
    imesh[i] = mapth[pTh];
  }
  if (!ThePlotStream) return 0L;

  PlotStream theplot(ThePlotStream);
  SendDualPlotArgs(theplot, nargs, s);
  theplot.SendMeshes( );
  theplot << 0L;
  theplot.SendMeshesS( );
  theplot << (long)mapth.size( );
  for (map< const MeshS *, long >::const_iterator i = mapth.begin( ); i != mapth.end( ); ++i)
    theplot << i->second << *i->first;
  theplot.SendPlots( );
  theplot << (long)n;
  for (size_t i = 0; i < n; ++i) {
    theplot << what[i] << imesh[i];    // what = 50: a meshS, 64, 65: a meshdualS (its primal meshS)
    if (what[i] == 65) theplot << u[i];
  }
  theplot.SendEndPlot( );
  return 0L;
}

//-------------------------------------------------------
// dual mesh for medit: each control volume as a fan of triangles (vertex 0 = the center of the
// control volume, the primal vertex; label = number of the control volume + 1), with the option
// "|dual" ffmedit draws only the edges of the control volumes (see fvdualmesh.idp, meditD)
//-------------------------------------------------------
// the fan of triangles of the control volumes; with primal = true, each interface crossed by an
// edge of the primal mesh is cut at the crossing point, so the primal edge is made of two edges
// of the fan (center - crossing point - center), kept as internal edges of label lprimal
struct DualFan {
  vector< R2 > v;          // the vertices
  vector< long > t;        // the triangles (3 vertices: center, a, b)
  vector< long > tcv;      // the control volume of each triangle
  vector< long > be;       // the edges (2 vertices) and their label
  vector< long > belab;
  long lprimal;            // the label of the primal edges
  DualFan(dualmesh *pdm, bool primal);
};

// the label of the primal edges in the fan: greatest label of the border + 1
static long DualPrimalLabel(dualmesh *pdm) {
  const KN< long > &dbelab = *pdm->normalized->dbelab;
  long l = 0;
  for (long i = 0; i < dbelab.N( ); ++i) l = max(l, abs(dbelab[i]));
  return l + 1;
}

DualFan::DualFan(dualmesh *pdm, bool primal) {
  if (!pdm || !pdm->primal) ExecError("dualmeshfan: the meshdual is not built (builddualmesh)");
  const Mesh &Th = *pdm->primal;
  const dmnormalized &dn = *pdm->normalized;
  const KN< long > &nvcv = *dn.ndvcv;
  long ndv = nvcv[0], ncv = nvcv[1], ndbe = nvcv[2];
  const KN< double > &xy = *dn.dvxy;
  const KN< long > &ninterf = *dn.ninterf;
  const KNM< long > &dvcv = *dn.dvcv;
  const KNM< long > &dvbe = *dn.dvbe;
  const KN< long > &dbelab = *dn.dbelab;
  R2 Pn, Px;
  Th.BoundingBox(Pn, Px);
  double eps = 1e-10 * (Px - Pn).norme( );
  lprimal = DualPrimalLabel(pdm);
  for (long k = 0; k < ndv; ++k) v.push_back(R2(xy[2 * k], xy[2 * k + 1]));
  // the center: the vertex of the polygon at the primal vertex (boundary), or a new vertex
  vector< long > center(ncv);
  for (long i = 0; i < ncv; ++i) {
    center[i] = -1;
    for (long j = 0; j < ninterf[i] && center[i] < 0; ++j) {
      long a = dvcv(j, i);
      if ((v[a] - Th(i)).norme( ) < eps) center[i] = a;
    }
    if (center[i] < 0) {
      center[i] = v.size( );
      v.push_back(Th(i));
    }
  }
  // the control volume on the other side of the edge (a, b) of a polygon: the owner of (b, a)
  map< pair< long, long >, long > owner;
  for (long i = 0; i < ncv; ++i)
    for (long j = 0; j < ninterf[i]; ++j)
      owner[make_pair(dvcv(j, i), dvcv((j + 1) % ninterf[i], i))] = i;
  // the crossing point of the primal edge (center i, center i2) and of the interface (a, b)
  map< pair< long, long >, long > cut;
  long nbent = 0;    // crossing point outside the primal edge (very obtuse triangles)
  if (primal)
    for (long i = 0; i < ncv; ++i)
      for (long j = 0; j < ninterf[i]; ++j) {
        long a = dvcv(j, i), b = dvcv((j + 1) % ninterf[i], i);
        map< pair< long, long >, long >::const_iterator o = owner.find(make_pair(b, a));
        if (o == owner.end( ) || o->second < i) continue;    // the border, or already done
        long i2 = o->second;
        R2 C1 = v[center[i]], C2 = v[center[i2]], A = v[a], B = v[b];
        R2 AB = B - A, CC = C2 - C1;
        double det = AB.x * CC.y - AB.y * CC.x;
        if (fabs(det) < 1e-30) continue;
        // A + s AB = C1 + r CC
        double s = ((C1.x - A.x) * CC.y - (C1.y - A.y) * CC.x) / det;
        if (s < 1e-6 || s > 1 - 1e-6) continue;    // at a vertex of the polygon (primal border)
        double r = (AB.y * (C1.x - A.x) - AB.x * (C1.y - A.y)) / det;    // = C1 + r CC
        if (r < -1e-9 || r > 1 + 1e-9) nbent++;
        long m = v.size( );
        v.push_back(A + AB * s);
        cut[make_pair(min(a, b), max(a, b))] = m;
        be.push_back(center[i]);
        be.push_back(m);
        belab.push_back(lprimal);
        be.push_back(m);
        be.push_back(center[i2]);
        belab.push_back(lprimal);
      }
  if (nbent && verbosity > 1)
    cout << " dualmeshfan: " << nbent << " primal edges drawn with a bend (obtuse triangles)" << endl;
  // the triangles (center, a, b), cut in two at the crossing point
  for (long i = 0; i < ncv; ++i) {
    long n = ninterf[i], c = center[i];
    for (long j = 0; j < n; ++j) {
      long a = dvcv(j, i), b = dvcv((j + 1) % n, i);
      // the edges of the polygon at the center (boundary); the center can be several vertices
      // of the polygon (vertex where several parts of the border meet: one copy by part)
      if ((v[a] - v[c]).norme( ) < eps || (v[b] - v[c]).norme( ) < eps) continue;
      map< pair< long, long >, long >::const_iterator q = cut.find(make_pair(min(a, b), max(a, b)));
      long tri[2][3] = {{c, a, b}, {-1, -1, -1}};
      if (q != cut.end( )) {
        tri[0][2] = q->second;
        tri[1][0] = c, tri[1][1] = q->second, tri[1][2] = b;
      }
      for (int h = 0; h < 2; ++h)
        if (tri[h][0] >= 0) {
          for (int l = 0; l < 3; ++l) t.push_back(tri[h][l]);
          tcv.push_back(i);
        }
    }
  }
  for (long e = 0; e < ndbe; ++e) {
    be.push_back(dvbe(0L, e));
    be.push_back(dvbe(1L, e));
    belab.push_back(dbelab[e]);
  }
}

// dualmeshfan(Thstar) or dualmeshfan(Thstar, 1): with the edges of the primal mesh (label
// dualmeshprimallabel(Thstar)); vertex 0 of the triangles = the center, label = control volume + 1
static Mesh *BuildFanMesh(dualmesh *pdm, bool primal) {
  DualFan F(pdm, primal);
  long nv = F.v.size( ), nt = F.tcv.size( ), nbe = F.belab.size( );
  Mesh::Vertex *v = new Mesh::Vertex[nv];
  Mesh::Triangle *t = new Mesh::Triangle[nt];
  Mesh::BorderElement *b = new Mesh::BorderElement[nbe];
  for (long k = 0; k < nv; ++k) {
    v[k].x = F.v[k].x;
    v[k].y = F.v[k].y;
    v[k].lab = 0;
  }
  for (long k = 0; k < nt; ++k) t[k].set(v, F.t[3 * k], F.t[3 * k + 1], F.t[3 * k + 2], F.tcv[k] + 1);
  for (long e = 0; e < nbe; ++e) b[e].set(v, F.be[2 * e], F.be[2 * e + 1], F.belab[e]);
  Mesh *pTh = new Mesh(nv, nt, nbe, v, t, b);
  R2 Pn, Px;
  pTh->BoundingBox(Pn, Px);
  if (!pTh->quadtree) pTh->quadtree = new Fem2D::FQuadTree(pTh, Pn, Px, pTh->nv);
  return pTh;
}
pmesh Dual_Mesh_fan(Stack stack, dualmesh *const &pdm, long const &primal) {
  return Add2StackOfPtr2FreeRC(stack, BuildFanMesh(pdm, primal != 0));
}
pmesh Dual_Mesh_fan1(Stack stack, dualmesh *const &pdm) { return Dual_Mesh_fan(stack, pdm, 0L); }

// the value of each triangle of the fan (P0 on dualmeshfan(Thstar, primal)): the value of its control volume
double Dual_Mesh_fanvalues(KN< double > *const &pw, dualmesh *const &pdm, KN_< double > const &u,
                           long const &primal) {
  DualFan F(pdm, primal != 0);
  if (u.N( ) != (*pdm->normalized->ndvcv)[1])
    ExecError("dualmeshfanvalues: u must have one value by control volume");
  KN< double > &w = *pw;
  w.resize(F.tcv.size( ));
  for (long k = 0; k < (long)F.tcv.size( ); ++k) w[k] = u[F.tcv[k]];
  return 0.;
}
double Dual_Mesh_fanvalues1(KN< double > *const &pw, dualmesh *const &pdm, KN_< double > const &u) {
  return Dual_Mesh_fanvalues(pw, pdm, u, 0L);
}
long Dual_Mesh_primallabel(dualmesh *const &pdm) { return DualPrimalLabel(pdm); }

//-------------------------------------------------------
// the P0 finite element on the dual mesh: fespace Wstar(Thstar, P0) (Thstar a meshdual, meshdualL or meshdual3)
//-------------------------------------------------------
// the control volumes as a mesh of FreeFEM, the label of an element = the number of its control volume + 1 (so a
// function P0 on this mesh is constant on a control volume if it has the same value on its elements: dualtoP0):
//   2D (meshdual, the dual mesh of F. Bouchut): each control volume is a fan of triangles (dualmeshfan)
//   1D (meshdualL, median): one segment by control volume, between the middles of the 2 segments around its vertex
//      (from the boundary point to the middle of its segment at a boundary point): the vertices of this mesh are the
//      middles of the segments and the boundary points, the element i is the control volume of the vertex i, so
//      Wstar.ndof = Th.nv and w[] = u[]; a vertex with more than 2 segments (a branching meshL) keeps one half
//      segment by segment (elements after the first Th.nv ones)
//   3D (meshdual3, median): each tetrahedron cut in 24 tetrahedra (vertex, middle of an edge, centroid of a face,
//      centroid of the tetrahedron), the boundary triangles in 6 triangles
//   surface (meshdualS, median): each triangle cut in 6 triangles
// built when it is used (the cast of the dual mesh to a mesh, for fespace, int2d, savevtk, ...) and kept in the
// dual mesh (destroyed with it, when no fespace uses it)
static MeshL *BuildSubMeshL(const dualmeshN &D) {
  const MeshL &Th = *(const MeshL *)D.primal;
  // the segments around each vertex: the first ones (sa, its local vertex la) and the second ones (sb, lb)
  KN< long > deg(Th.nv), sa(Th.nv), sb(Th.nv), la(Th.nv), lb(Th.nv);
  deg = 0L; sa = -1L; sb = -1L;
  for (long k = 0; k < Th.nt; ++k)
    for (int l = 0; l < 2; ++l) {
      long i = Th(k, l);
      if (deg[i] == 0) { sa[i] = k; la[i] = l; } else if (deg[i] == 1) { sb[i] = k; lb[i] = l; }
      deg[i]++;
    }
  // the vertices: the middles of the segments (0 .. nt-1), then the vertices of Th with a single segment (the
  // boundary points) or with more than 2 segments (branching)
  KN< long > iv(Th.nv); iv = -1L;
  long nv = Th.nt;
  for (long i = 0; i < Th.nv; ++i) if (deg[i] != 2) iv[i] = nv++;
  long nt = 0;
  for (long i = 0; i < Th.nv; ++i) nt += deg[i] <= 2 ? 1 : deg[i];
  long nbe = Th.nbe;
  MeshL::Vertex *v = new MeshL::Vertex[nv];
  MeshL::Element *t = new MeshL::Element[nt];
  MeshL::BorderElement *b = new MeshL::BorderElement[nbe];
  for (long k = 0; k < Th.nt; ++k) {
    R3 M = (R3(Th(Th(k, 0))) + R3(Th(Th(k, 1))))*0.5;
    v[k].x = M.x; v[k].y = M.y; v[k].z = M.z; v[k].lab = 0;
  }
  for (long i = 0; i < Th.nv; ++i)
    if (iv[i] >= 0) { long m = iv[i]; v[m].x = Th(i).x; v[m].y = Th(i).y; v[m].z = Th(i).z; v[m].lab = Th(i).lab; }
  // the element i: the control volume of the vertex i, oriented as the segments of Th
  long ne = Th.nv;
  for (long i = 0; i < Th.nv; ++i) {
    int e[2];
    if (deg[i] == 2) {// between the middles of its 2 segments: the segment where i is the end first
      if (la[i] == 1) { e[0] = (int)sa[i]; e[1] = (int)sb[i]; } else { e[0] = (int)sb[i]; e[1] = (int)sa[i]; }
      t[i].set(v, e, i + 1);
    } else if (deg[i] == 1) {// a boundary point: from the point to the middle of its segment
      if (la[i] == 0) { e[0] = (int)iv[i]; e[1] = (int)sa[i]; } else { e[0] = (int)sa[i]; e[1] = (int)iv[i]; }
      t[i].set(v, e, i + 1);
    } else {// branching (or isolated) vertex: one half segment by segment around it
      long first = 1;
      for (long k = 0; k < Th.nt; ++k)
        for (int l = 0; l < 2; ++l)
          if (Th(k, l) == i) {
            e[0] = (int)iv[i]; e[1] = (int)k;
            t[first ? i : ne++].set(v, e, i + 1);
            first = 0;
          }
      if (first) ExecError("builddualmesh (meshL): a vertex without segment");
    }
  }
  for (long ib = 0; ib < nbe; ++ib) {
    long i = Th(Th.be(ib)[0]);
    if (iv[i] < 0) ExecError("builddualmesh (meshL): a boundary point inside the curve");
    int ib1[1] = {(int)iv[i]};
    b[ib].set(v, ib1, Th.be(ib).lab);
  }
  MeshL *pTh = new MeshL(nv, nt, nbe, v, t, b);
  pTh->BuildGTree( );
  return pTh;
}
static Mesh3 *BuildSubMesh3(const dualmeshN &D) {
  const Mesh3 &Th = *(const Mesh3 *)D.primal;
  // the vertices: the vertices of Th, the middles of the edges, the centroids of the faces, of the tetrahedra
  map< pair< int, int >, int > edges;
  map< vector< int >, int > faces;
  vector< R3 > P;
  for (int i = 0; i < Th.nv; ++i) P.push_back(R3(Th(i)));
  for (int k = 0; k < Th.nt; ++k)
    for (int i = 0; i < 4; ++i)
      for (int j = i + 1; j < 4; ++j) {
        int a = Th(k, i), c = Th(k, j);
        pair< int, int > e(min(a, c), max(a, c));
        if (edges.find(e) == edges.end( )) { edges[e] = P.size( ); P.push_back((P[a] + P[c])*0.5); }
      }
  for (int k = 0; k < Th.nt; ++k)
    for (int f = 0; f < 4; ++f) {
      vector< int > q;
      for (int l = 0; l < 4; ++l) if (l != f) q.push_back(Th(k, l));
      vector< int > key(q);
      sort(key.begin( ), key.end( ));
      if (faces.find(key) == faces.end( )) { faces[key] = P.size( ); P.push_back((P[q[0]] + P[q[1]] + P[q[2]])/3.); }
    }
  int ng = P.size( );
  for (int k = 0; k < Th.nt; ++k) P.push_back((P[Th(k, 0)] + P[Th(k, 1)] + P[Th(k, 2)] + P[Th(k, 3)])*0.25);
  long nv = P.size( ), nt = 24*Th.nt, nbe = 6*Th.nbe;
  Mesh3::Vertex *v = new Mesh3::Vertex[nv];
  Mesh3::Element *t = new Mesh3::Element[nt];
  Mesh3::BorderElement *b = new Mesh3::BorderElement[nbe];
  for (long i = 0; i < nv; ++i) { v[i].x = P[i].x; v[i].y = P[i].y; v[i].z = P[i].z; v[i].lab = i < Th.nv ? Th(i).lab : 0; }
  long it = 0;
  for (int k = 0; k < Th.nt; ++k) {
    int g = ng + k;
    for (int l = 0; l < 4; ++l) {// the part of the vertex l: 6 tetrahedra on the 3 faces of K containing it
      int a = Th(k, l);
      for (int f = 0; f < 4; ++f) {
        if (f == l) continue;// the face f (opposite to the vertex f) contains a
        vector< int > key;
        for (int m = 0; m < 4; ++m) if (m != f) key.push_back(Th(k, m));
        sort(key.begin( ), key.end( ));
        int F = faces[key];
        for (int m = 0; m < 4; ++m) {
          if (m == l || m == f) continue;// the edge (a, Th(k,m)) of the face f
          int c = Th(k, m), M = edges[pair< int, int >(min(a, c), max(a, c))];
          int iv[4] = {a, M, F, g};
          if (det(P[iv[0]], P[iv[1]], P[iv[2]], P[iv[3]]) < 0) swap(iv[1], iv[2]);
          t[it++].set(v, iv, a + 1);
        }
      }
    }
  }
  long ib = 0;
  for (int e = 0; e < Th.nbe; ++e) {// the boundary triangle: 6 triangles, the orientation of the triangle
    const Mesh3::BorderElement &B = Th.be(e);
    int q[3] = {(int)Th(B[0]), (int)Th(B[1]), (int)Th(B[2])};
    vector< int > key(q, q + 3);
    sort(key.begin( ), key.end( ));
    int F = faces[key];
    R3 n0 = (P[q[1]] - P[q[0]]) ^ (P[q[2]] - P[q[0]]);
    for (int l = 0; l < 3; ++l)
      for (int h = 1; h <= 2; ++h) {
        int a = q[l], c = q[(l + h)%3], M = edges[pair< int, int >(min(a, c), max(a, c))];
        int iv[3] = {a, M, F};
        if ((((P[M] - P[a]) ^ (P[F] - P[a])), n0) < 0) swap(iv[1], iv[2]);
        b[ib++].set(v, iv, B.lab);
      }
  }
  Mesh3 *pTh = new Mesh3(nv, nt, nbe, v, t, b);
  pTh->BuildGTree( );
  return pTh;
}
//   surface (meshdualS, median): each triangle cut in 6 triangles (vertex, middle of an edge, centroid), the
//      boundary edges in 2 halves
static MeshS *BuildSubMeshS(const dualmeshN &D) {
  const MeshS &Th = *(const MeshS *)D.primal;
  // the vertices: the vertices of Th, the middles of the edges, the centroids of the triangles
  map< pair< int, int >, int > edges;
  vector< R3 > P;
  for (int i = 0; i < Th.nv; ++i) P.push_back(R3(Th(i)));
  for (int k = 0; k < Th.nt; ++k)
    for (int i = 0; i < 3; ++i) {
      int a = Th(k, i), c = Th(k, (i + 1) % 3);
      pair< int, int > e(min(a, c), max(a, c));
      if (edges.find(e) == edges.end( )) { edges[e] = P.size( ); P.push_back((P[a] + P[c])*0.5); }
    }
  int ng = P.size( );
  for (int k = 0; k < Th.nt; ++k) P.push_back((P[Th(k, 0)] + P[Th(k, 1)] + P[Th(k, 2)])/3.);
  long nv = P.size( ), nt = 6*Th.nt, nbe = 2*Th.nbe;
  MeshS::Vertex *v = new MeshS::Vertex[nv];
  MeshS::Element *t = new MeshS::Element[nt];
  MeshS::BorderElement *b = new MeshS::BorderElement[nbe];
  for (long i = 0; i < nv; ++i) { v[i].x = P[i].x; v[i].y = P[i].y; v[i].z = P[i].z; v[i].lab = i < Th.nv ? Th(i).lab : 0; }
  for (int k = 0; k < Th.nt; ++k) {// the part of the vertex a in K: (a, middle of [a,b], G), (a, G, middle of [a,c])
    int G = ng + k;
    for (int l = 0; l < 3; ++l) {
      int a = Th(k, l), c1 = Th(k, (l + 1) % 3), c2 = Th(k, (l + 2) % 3);
      int M1 = edges[pair< int, int >(min(a, c1), max(a, c1))], M2 = edges[pair< int, int >(min(a, c2), max(a, c2))];
      int t0[3] = {a, M1, G}, t1[3] = {a, G, M2};// the vertex of the control volume first (vertex 0, see P0dual)
      t[6*k+2*l].set(v, t0, a + 1);
      t[6*k+2*l+1].set(v, t1, a + 1);
    }
  }
  for (int ib = 0; ib < Th.nbe; ++ib) {// the 2 halves of each boundary edge
    int a = Th(Th.be(ib)[0]), c = Th(Th.be(ib)[1]);
    int M = edges[pair< int, int >(min(a, c), max(a, c))];
    int e0[2] = {a, M}, e1[2] = {M, c};
    b[2*ib].set(v, e0, Th.be(ib).lab);
    b[2*ib+1].set(v, e1, Th.be(ib).lab);
  }
  MeshS *pTh = new MeshS(nv, nt, nbe, v, t, b);
  pTh->BuildGTree( );
  return pTh;
}
// the P0 finite element of the dual mesh: fespace Wstar(Thstar, P0) (P0 replaced by P0dual in lgfem.cpp when the
// mesh is a dual mesh): the P0 of FreeFEM, interpolated at the vertex 0 of the element (the vertex of its control
// volume) and not at its centroid, so the function is constant on each control volume (Wstar w = x+y: the value at
// the vertex of the control volume, as the unknowns of the finite volumes)
class TypeOfFE_P0dual : public TypeOfFE {// 2D
 public:
  static int Data[];
  static double Pi_h_coef[];
  TypeOfFE_P0dual( ) : TypeOfFE(0, 0, 1, 1, Data, 1, 1, 1, 1, Pi_h_coef) {
    pij_alpha[0] = IPJ(0, 0, 0);
    P_Pi_h[0] = R2(0., 0.);// the vertex 0
  }
  void FB(const bool *whatd, const Mesh &, const Triangle &, const R2 &, RNMK_ &val) const {
    val = 0;
    if (whatd[op_id]) val(0, 0, 0) = 1;
  }
};
int TypeOfFE_P0dual::Data[] = {6, 0, 0, 0, 0, 0, 0, 1};
double TypeOfFE_P0dual::Pi_h_coef[] = {1.0};
// the interpolation point of P0dual: 3D, surface: the vertex 0 of the element (the vertex of its control volume);
// 1D: the centre of the element (the element is the whole control volume, its vertices are the middles of the
// segments: the centre is the vertex of the control volume on a uniform mesh)
template< class MMesh > typename MMesh::Element::RdHat P0dualPoint( ) { return typename MMesh::Element::RdHat( ); }
template<> R1 P0dualPoint< MeshL >( ) { return R1(0.5); }
template< class MMesh >
class TypeOfFE_P0dualN : public TypeOfFE_Lagrange< MMesh > {// 1D (meshL), 3D (mesh3), surface (meshS)
 public:
  typedef typename MMesh::Element Element;
  typedef typename Element::RdHat RdHat;
  typedef GFElement< MMesh > FElement;
  TypeOfFE_P0dualN( ) : TypeOfFE_Lagrange< MMesh >(0) { this->PtInterpolation[0] = P0dualPoint< MMesh >( ); }
  void FB(const What_d whatd, const MMesh &, const Element &, const RdHat &, RNMK_ &val) const {
    val = 0;
    RN_ f0(val('.', 0, op_id));
    if (whatd & Fop_D0) f0[0] = 1;
  }
  R operator( )(const FElement &K, const RdHat &, const KN_< R > &u, int, int op) const { return op == 0 ? u(K(0)) : 0.; }
};
static TypeOfFE_P0dual The_P0dual2;
static TypeOfFE_P0dualN< Mesh3 > The_P0dual3;
static TypeOfFE_P0dualN< MeshL > The_P0dualL;
static TypeOfFE_P0dualN< MeshS > The_P0dualS;
const TypeOfFE *P0dualFE2 = &The_P0dual2;
const GTypeOfFE< Mesh3 > *P0dualFE3 = &The_P0dual3;
const GTypeOfFE< MeshL > *P0dualFEL = &The_P0dualL;
const GTypeOfFE< MeshS > *P0dualFES = &The_P0dualS;

// the casts of the dual mesh to its mesh of the control volumes (fespace Wstar(Thstar, P0), int2d(Thstar), ...)
// the meshes of the control volumes of the 2D dual meshes and their dual mesh: plot(w) with w in Wstar(Thstar, P0)
// is the plot of the dual mesh (the polygons of the control volumes), not the one of a function on the fans
static map< const void *, dualmesh * > DualSubs;
static map< const void *, dualmesh3 * > DualSubs3;// the same in 3D
static map< const void *, dualmeshS * > DualSubsS;// the same on a surface
void DualSubForget(const void *sub) { DualSubs.erase(sub); DualSubs3.erase(sub); DualSubsS.erase(sub); }
// the dual meshes and their primal mesh: dual -> (primal, dimension), from build to destroy; allocated with the
// first dual mesh, freed with the last one (not a static map: a dual mesh can be destroyed at the end of the
// program, after the static variables)
typedef map< const void *, pair< const void *, int > > DualPrimalMap;
static DualPrimalMap *DualPrimals = 0;
void DualOfPrimalAdd(const void *primal, const void *dual, int dim) {
  if (!DualPrimals) DualPrimals = new DualPrimalMap;
  (*DualPrimals)[dual] = make_pair(primal, dim);
}
void DualOfPrimalForget(const void *dual) {
  if (!DualPrimals) return;
  DualPrimals->erase(dual);
  if (DualPrimals->empty()) { delete DualPrimals; DualPrimals = 0; }
}
bool DualIsBuilt(const void *dual) { return DualPrimals && DualPrimals->count(dual); }
void DualsOfPrimal(const void *primal, int dim, std::vector< const void * > &duals) {
  duals.clear();
  if (!DualPrimals) return;
  for (DualPrimalMap::const_iterator it = DualPrimals->begin(); it != DualPrimals->end(); ++it)
    if (it->second.first == primal && it->second.second == dim) duals.push_back(it->first);
}
static pmesh *Dual_Mesh_sub(dualmesh *const &pdm) {
  if (!pdm || !pdm->primal) ExecError("meshdual: the dual mesh is not built (builddualmesh)");
  if (!pdm->sub) { pdm->sub = BuildFanMesh(pdm, false); DualSubs[pdm->sub] = pdm; }
  return &pdm->sub;
}
static pmesh Dual_Mesh_subv(dualmesh *const &pdm) { return *Dual_Mesh_sub(pdm); }
static pmeshL *DualL_sub(dualmeshL *const &pdm) {
  if (!pdm || !pdm->dim) ExecError("meshdualL: the dual mesh is not built (builddualmesh)");
  if (!pdm->subL) pdm->subL = BuildSubMeshL(*pdm);
  return &pdm->subL;
}
static pmeshL DualL_subv(dualmeshL *const &pdm) { return *DualL_sub(pdm); }
static pmesh3 *Dual3_sub(dualmesh3 *const &pdm) {
  if (!pdm || !pdm->dim) ExecError("meshdual3: the dual mesh is not built (builddualmesh)");
  if (!pdm->sub3) { pdm->sub3 = BuildSubMesh3(*pdm); DualSubs3[pdm->sub3] = pdm; }
  return &pdm->sub3;
}
static pmesh3 Dual3_subv(dualmesh3 *const &pdm) { return *Dual3_sub(pdm); }
static pmeshS *DualS_sub(dualmeshS *const &pdm) {
  if (!pdm || !pdm->dim) ExecError("meshdualS: the dual mesh is not built (builddualmesh)");
  if (!pdm->subS) { pdm->subS = BuildSubMeshS(*pdm); DualSubsS[pdm->subS] = pdm; }
  return &pdm->subS;
}
static pmeshS DualS_subv(dualmeshS *const &pdm) { return *DualS_sub(pdm); }
// the functions of the casts (signature of the casts of FreeFEM, as UnRef in lgfem.cpp)
template< class R, class D, R (*f)(D *const &) >
static AnyType DualCast(Stack, const AnyType &a) { return SetAny< R >(f(GetAny< D * >(a))); }
// the measure of an element
static double ElemMes(const Mesh::Element &K) { return K.area; }
template< class E > static double ElemMes(const E &K) { return K.mesure( ); }
// dualtoP0(w[], Thstar, u[]): u one value by control volume (the unknowns of the finite volumes, P1 of Th) -> w the
// function P0 of fespace Wstar(Thstar, P0) (the value of each element = the value of its control volume)
// dualfromP0(u[], Thstar, w[]): the mean of w on each control volume
template< class M >
static double DualToP0(KN< double > *const &pw, const M *pTh, long ncv, KN_< double > const &u) {
  if (u.N( ) != ncv) ExecError("dualtoP0: u must have one value by control volume (Th.nv)");
  KN< double > &w = *pw;
  w.resize(pTh->nt);
  for (long k = 0; k < pTh->nt; ++k) w[k] = u[(*pTh)[k].lab - 1];
  return 0.;
}
template< class M >
static double DualFromP0(KN< double > *const &pu, const M *pTh, long ncv, KN_< double > const &w) {
  if (w.N( ) != pTh->nt) ExecError("dualfromP0: w must be a function of fespace Wstar(Thstar, P0)");
  KN< double > &u = *pu, mes(ncv);
  u.resize(ncv); u = 0.; mes = 0.;
  for (long k = 0; k < pTh->nt; ++k) {
    long i = (*pTh)[k].lab - 1;
    double mk = ElemMes((*pTh)[k]);
    u[i] += mk*w[k]; mes[i] += mk;
  }
  for (long i = 0; i < ncv; ++i) if (mes[i] > 0.) u[i] /= mes[i];
  return 0.;
}
static double Dual_toP0(KN< double > *const &pw, dualmesh *const &pdm, KN_< double > const &u) {
  return DualToP0(pw, *Dual_Mesh_sub(pdm), pdm->primal->nv, u);
}
static double Dual_fromP0(KN< double > *const &pu, dualmesh *const &pdm, KN_< double > const &w) {
  return DualFromP0(pu, *Dual_Mesh_sub(pdm), pdm->primal->nv, w);
}
static double DualL_toP0(KN< double > *const &pw, dualmeshL *const &pdm, KN_< double > const &u) {
  return DualToP0(pw, *DualL_sub(pdm), pdm->nv, u);
}
static double DualL_fromP0(KN< double > *const &pu, dualmeshL *const &pdm, KN_< double > const &w) {
  return DualFromP0(pu, *DualL_sub(pdm), pdm->nv, w);
}
static double Dual3_toP0(KN< double > *const &pw, dualmesh3 *const &pdm, KN_< double > const &u) {
  return DualToP0(pw, *Dual3_sub(pdm), pdm->nv, u);
}
static double Dual3_fromP0(KN< double > *const &pu, dualmesh3 *const &pdm, KN_< double > const &w) {
  return DualFromP0(pu, *Dual3_sub(pdm), pdm->nv, w);
}
static double DualS_toP0(KN< double > *const &pw, dualmeshS *const &pdm, KN_< double > const &u) {
  return DualToP0(pw, *DualS_sub(pdm), pdm->nv, u);
}
static double DualS_fromP0(KN< double > *const &pu, dualmeshS *const &pdm, KN_< double > const &w) {
  return DualFromP0(pu, *DualS_sub(pdm), pdm->nv, w);
}

// plot(w) with w in fespace Wstar(Thstar, P0) (2D, called by Plot in lgfem.cpp): sent as the item of the dual mesh
// with the value of each control volume (what = 61: the polygons, the key m of ffglut draws their edges, as
// plot(Thstar, u[])); false if th is not the mesh of the control volumes of a dual mesh or w is not P0 on it
bool IsDualSubP0(const Mesh *th, long ndof) { return th && DualSubs.count(th) && ndof == th->nt; }
bool DualPlotP0(PlotStream &theplot, const Mesh *th, const KN_< double > &w) {
  if (!IsDualSubP0(th, w.N( ))) return false;
  dualmesh *pdm = DualSubs[th];
  KN< double > u(pdm->primal->nv);
  u = 0.;
  for (long k = 0; k < th->nt; ++k) u[(*th)[k].lab - 1] = w[k];
  SendDualMesh(theplot, pdm, u);
  return true;
}
// plot(w) with w in fespace Wstar(Thstar, P0), 3D (Thstar a meshdual3, called by Plot in lgfem.cpp): the item of the
// dual mesh, what = 63 with the number of its primal mesh3 (sent with the meshes of the plot in place of the mesh of
// the control volumes) and the value of each control volume (the key m of ffglut draws the edges of the control
// volumes, F1 the cutting plane, as plot(Thstar, u[])); DualSub3Primal: the primal mesh3, 0 if th is not the mesh
// of the control volumes of a meshdual3 or w is not P0 on it
const Mesh3 *DualSub3Primal(const Mesh3 *th, long ndof) {
  if (!th || !DualSubs3.count(th) || ndof != th->nt) return 0;
  return (const Mesh3 *)DualSubs3[th]->primal;
}
bool DualPlotP0_3(PlotStream &theplot, const Mesh3 *th, const KN_< double > &w, long imesh) {
  if (!DualSub3Primal(th, w.N( ))) return false;
  KN< double > u(DualSubs3[th]->nv);
  u = 0.;
  for (long k = 0; k < th->nt; ++k) u[(*th)[k].lab - 1] = w[k];
  theplot << 63L << imesh << u;
  return true;
}
// the same on a surface (Thstar a meshdualS): what = 65 with the number of its primal meshS
const MeshS *DualSubSPrimal(const MeshS *th, long ndof) {
  if (!th || !DualSubsS.count(th) || ndof != th->nt) return 0;
  return (const MeshS *)DualSubsS[th]->primal;
}
bool DualPlotP0_S(PlotStream &theplot, const MeshS *th, const KN_< double > &w, long imesh) {
  if (!DualSubSPrimal(th, w.N( ))) return false;
  KN< double > u(DualSubsS[th]->nv);
  u = 0.;
  for (long k = 0; k < th->nt; ++k) u[(*th)[k].lab - 1] = w[k];
  theplot << 65L << imesh << u;
  return true;
}

void Load_Init_fvdualmesh () {
    if (verbosity && !mpirank) cout << "fvdualmesh ";// a kernel module: its name in the line " Load: lg_fem ..." (lg.ypp)
//    zzzfff->Add("R2", atype< R2 * >( ));// see src/fflib/lgfem.cpp line 6316
    Dcl_Type< dualmesh * >(InitP< dualmesh >, Destroy< dualmesh >);
    zzzfff->Add("dualmesh",atype<dualmesh*>());
    zzzfff->Add("meshdual",atype<dualmesh*>());// the name of the type, as meshL, meshS, mesh3
    // priority 1: before the plot of FreeFEM, which accepts all the types
    Global.Add("plot", "(", new OneOperatorCode< PlotDualMeshN< 1 >, 1 >);// plot(Thstar)
    Global.Add("plot", "(", new OneOperatorCode< PlotDualMeshN< 2 >, 1 >);// plot(Thstar, u[])
    Global.Add("plot", "(", new OneOperatorCode< PlotDualMeshN< 3 >, 1 >);// plot(Th, Thstar)
    Global.Add("plot", "(", new OneOperatorCode< PlotDualMeshN< 4 >, 1 >);// plot(Thstar, Th)
    Global.Add("plot", "(", new OneOperatorCode< PlotDualMeshN< 5 >, 1 >);// plot(Th, Thstar, u[])
    Global.Add("plot", "(", new OneOperatorCode< PlotDualMeshN< 6 >, 1 >);// plot(Thstar, u[], Th)
    Global.Add("dualmeshfan", "(", new OneOperator1s_< pmesh, dualmesh * const >(Dual_Mesh_fan1));
    Global.Add("dualmeshfan", "(", new OneOperator2s_< pmesh, dualmesh * const, long >(Dual_Mesh_fan));
    Global.Add("dualmeshfanvalues", "(", new OneOperator3_< double, KN< double > * const, dualmesh * const, KN_< double > >(Dual_Mesh_fanvalues1));
    Global.Add("dualmeshfanvalues", "(", new OneOperator4_< double, KN< double > * const, dualmesh * const, KN_< double >, long >(Dual_Mesh_fanvalues));
    Global.Add("dualmeshprimallabel", "(", new OneOperator1_< long, dualmesh * const >(Dual_Mesh_primallabel));
    Global.Add("builddualmesh", "(", new OneOperator2_<double, dualmesh * const, const Mesh * const>  (Build_Dual_Mesh));
    Dcl_Type< DualMeshOf >();
    Global.Add("builddualmesh", "(", new OneOperator1_<DualMeshOf, const Mesh * const>  (Dual_Mesh_of));
    TheOperators->Add("<-", new OneOperator2_<dualmesh *, dualmesh *, DualMeshOf>  (Dual_Mesh_init));// meshdual Thstar = builddualmesh(Th);
    TheOperators->Add("=", new OneOperator2_<dualmesh *, dualmesh *, DualMeshOf>  (Dual_Mesh_set));// Thstar = builddualmesh(Th);
    // 1D (meshdualL, of a meshL) and 3D (meshdual3, of a mesh3): the median dual mesh
    Dcl_Type< dualmeshL * >(InitP< dualmeshL >, Destroy< dualmeshL >);
    zzzfff->Add("meshdualL",atype<dualmeshL*>());
    Dcl_Type< dualmesh3 * >(InitP< dualmesh3 >, Destroy< dualmesh3 >);
    zzzfff->Add("meshdual3",atype<dualmesh3*>());
    Dcl_Type< dualmeshS * >(InitP< dualmeshS >, Destroy< dualmeshS >);// surface (meshS): the median dual mesh
    zzzfff->Add("meshdualS",atype<dualmeshS*>());
    Dcl_Type< DualMeshOfN<MeshL> >();
    Dcl_Type< DualMeshOfN<Mesh3> >();
    Dcl_Type< DualMeshOfN<MeshS> >();
    Global.Add("builddualmesh", "(", new OneOperator2_<double, dualmeshL * const, const MeshL * const>  (Build_DualN<dualmeshL,MeshL>));
    Global.Add("builddualmesh", "(", new OneOperator2_<double, dualmesh3 * const, const Mesh3 * const>  (Build_DualN<dualmesh3,Mesh3>));
    Global.Add("builddualmesh", "(", new OneOperator1_<DualMeshOfN<MeshL>, const MeshL * const>  (Dual_MeshN_of<MeshL>));
    Global.Add("builddualmesh", "(", new OneOperator1_<DualMeshOfN<Mesh3>, const Mesh3 * const>  (Dual_MeshN_of<Mesh3>));
    Global.Add("builddualmesh", "(", new OneOperator2_<double, dualmeshS * const, const MeshS * const>  (Build_DualN<dualmeshS,MeshS>));
    Global.Add("builddualmesh", "(", new OneOperator1_<DualMeshOfN<MeshS>, const MeshS * const>  (Dual_MeshN_of<MeshS>));
    TheOperators->Add("<-", new OneOperator2_<dualmeshL *, dualmeshL *, DualMeshOfN<MeshL> >  (Dual_MeshN_init<dualmeshL,MeshL>));// meshdualL Thstar = builddualmesh(ThL);
    TheOperators->Add("=", new OneOperator2_<dualmeshL *, dualmeshL *, DualMeshOfN<MeshL> >  (Dual_MeshN_set<dualmeshL,MeshL>));
    TheOperators->Add("<-", new OneOperator2_<dualmesh3 *, dualmesh3 *, DualMeshOfN<Mesh3> >  (Dual_MeshN_init<dualmesh3,Mesh3>));// meshdual3 Thstar = builddualmesh(Th3);
    TheOperators->Add("=", new OneOperator2_<dualmesh3 *, dualmesh3 *, DualMeshOfN<Mesh3> >  (Dual_MeshN_set<dualmesh3,Mesh3>));
    TheOperators->Add("<-", new OneOperator2_<dualmeshS *, dualmeshS *, DualMeshOfN<MeshS> >  (Dual_MeshN_init<dualmeshS,MeshS>));// meshdualS Thstar = builddualmesh(ThS);
    TheOperators->Add("=", new OneOperator2_<dualmeshS *, dualmeshS *, DualMeshOfN<MeshS> >  (Dual_MeshN_set<dualmeshS,MeshS>));
    Add< dualmeshL * >("ncv", ".", new OneOperator1<long, dualmeshL * const>  (DualN_ncv<dualmeshL>));
    Add< dualmeshL * >("npe", ".", new OneOperator1<long, dualmeshL * const>  (DualN_npe<dualmeshL>));
    Add< dualmeshL * >("nde", ".", new OneOperator1<long, dualmeshL * const>  (DualN_nde<dualmeshL>));
    Add< dualmeshL * >("ndbe", ".", new OneOperator1<long, dualmeshL * const>  (DualN_ndbe<dualmeshL>));
    Add< dualmesh3 * >("ncv", ".", new OneOperator1<long, dualmesh3 * const>  (DualN_ncv<dualmesh3>));
    Add< dualmesh3 * >("npe", ".", new OneOperator1<long, dualmesh3 * const>  (DualN_npe<dualmesh3>));
    Add< dualmesh3 * >("nde", ".", new OneOperator1<long, dualmesh3 * const>  (DualN_nde<dualmesh3>));
    Add< dualmesh3 * >("ndbe", ".", new OneOperator1<long, dualmesh3 * const>  (DualN_ndbe<dualmesh3>));
    Add< dualmeshS * >("ncv", ".", new OneOperator1<long, dualmeshS * const>  (DualN_ncv<dualmeshS>));
    Add< dualmeshS * >("npe", ".", new OneOperator1<long, dualmeshS * const>  (DualN_npe<dualmeshS>));
    Add< dualmeshS * >("nde", ".", new OneOperator1<long, dualmeshS * const>  (DualN_nde<dualmeshS>));
    Add< dualmeshS * >("ndbe", ".", new OneOperator1<long, dualmeshS * const>  (DualN_ndbe<dualmeshS>));
    Global.Add("dualmeshcvarea", "(", new OneOperator2_<double, KN<double> * const, dualmeshL * const>  (DualN_cvarea<dualmeshL>));
    Global.Add("dualmeshcvarea", "(", new OneOperator2_<double, KN<double> * const, dualmesh3 * const>  (DualN_cvarea<dualmesh3>));
    Global.Add("dualmeshcvarea", "(", new OneOperator2_<double, KN<double> * const, dualmeshS * const>  (DualN_cvarea<dualmeshS>));
    // the P0 finite element on the dual mesh: the casts of the dual mesh to the mesh of its control volumes
    // (fespace Wstar(Thstar, P0), int2d(Thstar)(w), savevtk(..., Thstar, w), ...), dualtoP0, dualfromP0
    atype< pmesh * >( )->AddCast(new E_F1_funcT< pmesh *, dualmesh * >(DualCast< pmesh *, dualmesh, Dual_Mesh_sub >));
    atype< pmesh >( )->AddCast(new E_F1_funcT< pmesh, dualmesh * >(DualCast< pmesh, dualmesh, Dual_Mesh_subv >));
    atype< pmeshL * >( )->AddCast(new E_F1_funcT< pmeshL *, dualmeshL * >(DualCast< pmeshL *, dualmeshL, DualL_sub >));
    atype< pmeshL >( )->AddCast(new E_F1_funcT< pmeshL, dualmeshL * >(DualCast< pmeshL, dualmeshL, DualL_subv >));
    atype< pmesh3 * >( )->AddCast(new E_F1_funcT< pmesh3 *, dualmesh3 * >(DualCast< pmesh3 *, dualmesh3, Dual3_sub >));
    atype< pmesh3 >( )->AddCast(new E_F1_funcT< pmesh3, dualmesh3 * >(DualCast< pmesh3, dualmesh3, Dual3_subv >));
    atype< pmeshS * >( )->AddCast(new E_F1_funcT< pmeshS *, dualmeshS * >(DualCast< pmeshS *, dualmeshS, DualS_sub >));
    atype< pmeshS >( )->AddCast(new E_F1_funcT< pmeshS, dualmeshS * >(DualCast< pmeshS, dualmeshS, DualS_subv >));
    Global.Add("dualtoP0", "(", new OneOperator3_< double, KN< double > * const, dualmesh * const, KN_< double > >(Dual_toP0));
    Global.Add("dualtoP0", "(", new OneOperator3_< double, KN< double > * const, dualmeshL * const, KN_< double > >(DualL_toP0));
    Global.Add("dualtoP0", "(", new OneOperator3_< double, KN< double > * const, dualmesh3 * const, KN_< double > >(Dual3_toP0));
    Global.Add("dualfromP0", "(", new OneOperator3_< double, KN< double > * const, dualmesh * const, KN_< double > >(Dual_fromP0));
    Global.Add("dualfromP0", "(", new OneOperator3_< double, KN< double > * const, dualmeshL * const, KN_< double > >(DualL_fromP0));
    Global.Add("dualfromP0", "(", new OneOperator3_< double, KN< double > * const, dualmesh3 * const, KN_< double > >(Dual3_fromP0));
    Global.Add("dualtoP0", "(", new OneOperator3_< double, KN< double > * const, dualmeshS * const, KN_< double > >(DualS_toP0));
    Global.Add("dualfromP0", "(", new OneOperator3_< double, KN< double > * const, dualmeshS * const, KN_< double > >(DualS_fromP0));
    Global.Add("plot", "(", new OneOperatorCode< PlotDualMesh3N< 1 >, 1 >);// plot(Thstar)
    Global.Add("plot", "(", new OneOperatorCode< PlotDualMesh3N< 2 >, 1 >);// plot(Thstar, u[])
    Global.Add("plot", "(", new OneOperatorCode< PlotDualMesh3N< 3 >, 1 >);// plot(Th3, Thstar)
    Global.Add("plot", "(", new OneOperatorCode< PlotDualMesh3N< 4 >, 1 >);// plot(Thstar, Th3)
    Global.Add("plot", "(", new OneOperatorCode< PlotDualMesh3N< 5 >, 1 >);// plot(Th3, Thstar, u[])
    Global.Add("plot", "(", new OneOperatorCode< PlotDualMesh3N< 6 >, 1 >);// plot(Thstar, u[], Th3)
    Global.Add("plot", "(", new OneOperatorCode< PlotDualMeshSN< 1 >, 1 >);// plot(Thstar), Thstar a meshdualS
    Global.Add("plot", "(", new OneOperatorCode< PlotDualMeshSN< 2 >, 1 >);// plot(Thstar, u[])
    Global.Add("plot", "(", new OneOperatorCode< PlotDualMeshSN< 3 >, 1 >);// plot(ThS, Thstar)
    Global.Add("plot", "(", new OneOperatorCode< PlotDualMeshSN< 4 >, 1 >);// plot(Thstar, ThS)
    Global.Add("plot", "(", new OneOperatorCode< PlotDualMeshSN< 5 >, 1 >);// plot(ThS, Thstar, u[])
    Global.Add("plot", "(", new OneOperatorCode< PlotDualMeshSN< 6 >, 1 >);// plot(Thstar, u[], ThS)
    Add< dualmesh * >("npe", ".", new OneOperator1<long, dualmesh * const>  (Dual_Mesh_npe));
    Add< dualmesh * >("nde", ".", new OneOperator1<long, dualmesh * const>  (Dual_Mesh_nde));
    Add< dualmesh * >("ndv", ".", new OneOperator1<long, dualmesh * const>  (Dual_Mesh_ndv));
    Add< dualmesh * >("ncv", ".", new OneOperator1<long, dualmesh * const>  (Dual_Mesh_ncv));
    Add< dualmesh * >("ndbe", ".", new OneOperator1<long, dualmesh * const>  (Dual_Mesh_ndbe));
    Add< dualmesh * >("maxinterf", ".", new OneOperator1<long, dualmesh * const>  (Dual_Mesh_maxinterf));
    Global.Add("dualmeshdvxy", "(", new OneOperator2_<double, KN<double> * const, dualmesh * const>  (Dual_Mesh_dvxy));
    Global.Add("dualmeshninterf", "(", new OneOperator2_<double, KN<long> * const, dualmesh * const>  (Dual_Mesh_ninterf));
    Global.Add("dualmeshdvcv", "(", new OneOperator2_<double, KNM<long> * const, dualmesh * const>  (Dual_Mesh_dvcv));
    Global.Add("dualmeshdvbe", "(", new OneOperator2_<double, KNM<long> * const, dualmesh * const>  (Dual_Mesh_dvbe));
    Global.Add("dualmeshdbelab", "(", new OneOperator2_<double, KN<long> * const, dualmesh * const>  (Dual_Mesh_dbelab));
    Global.Add("dualmeshinterf", "(", new OneOperator2_<double, KNM<long> * const, dualmesh * const>  (Dual_Mesh_interf));
    Global.Add("dualmeshorientation", "(", new OneOperator2_<double, KNM<long> * const, dualmesh * const>  (Dual_Mesh_orientation));
    Global.Add("dualmeshedgenumbering", "(", new OneOperator2_<double, KN<long> * const, dualmesh * const>  (Dual_Mesh_edgenumbering));
    Global.Add("dualmeshcvLR", "(", new OneOperator2_<double, KN<long> * const, dualmesh * const>  (Dual_Mesh_cvLR));
    Global.Add("dualmeshdualpiece01", "(", new OneOperator2_<double, KNM<long> * const, dualmesh * const>  (Dual_Mesh_dualpiece01));
    Global.Add("dualmeshcvarea", "(", new OneOperator2_<double, KN<double> * const, dualmesh * const>  (Dual_Mesh_cvarea));
    Global.Add("dualmeshcvgrav", "(", new OneOperator2_<double, KN<double> * const, dualmesh * const>  (Dual_Mesh_cvgrav));
    Global.Add("dualmeshcutedge", "(", new OneOperator2_<double, KN<double> * const, dualmesh * const>  (Dual_Mesh_cutedge));
    Global.Add("dualmeshnormals", "(", new OneOperator3_<double, KN<double> * const, KN<double> * const, dualmesh * const>  (Dual_Mesh_normals));
    Global.Add("dualmeshareaminLR", "(", new OneOperator2_<double, KN<double> * const, dualmesh * const>  (Dual_Mesh_areaminLR));
    Global.Add("savedualmesh", "(", new OneOperator2_<double, dualmesh * const, string * const>  (Dual_Mesh_save));
    Global.Add("dualmeshmatvert2grav", "(", new OneOperator1s_<newpMatrice_Creuse<double>, dualmesh * const>  (Mat_vert2grav));
    Global.Add("dualmeshmatpiece2cvL", "(", new OneOperator1s_<newpMatrice_Creuse<double>, dualmesh * const>  (Mat_piece2cvL));
    Global.Add("dualmeshmatpiece2cvR", "(", new OneOperator1s_<newpMatrice_Creuse<double>, dualmesh * const>  (Mat_piece2cvR));
    Global.Add("dualmeshmatdxygrav", "(", new OneOperator1s_<newpMatrice_Creuse<double>, dualmesh * const>  (Mat_dxygrav));
    Global.Add("dualmeshmatpiece2cvdxyL", "(", new OneOperator1s_<newpMatrice_Creuse<double>, dualmesh * const>  (Mat_piece2cvdxyL));
    Global.Add("dualmeshmatpiece2cvdxyR", "(", new OneOperator1s_<newpMatrice_Creuse<double>, dualmesh * const>  (Mat_piece2cvdxyR));
    Global.Add("dualmeshmatpiece2cvprodL", "(", new OneOperator1s_<newpMatrice_Creuse<double>, dualmesh * const>  (Mat_piece2cvprodL));
    Global.Add("dualmeshmatpiece2cvprodR", "(", new OneOperator1s_<newpMatrice_Creuse<double>, dualmesh * const>  (Mat_piece2cvprodR));
    Global.Add("dualmeshmattri2cvprod", "(", new OneOperator3s_<newpMatrice_Creuse<double>, KN<long> * const, KN<double> * const, dualmesh * const>  (Mat_tri2cvprod));
    Global.Add("dualmeshNbdryvert", "(", new OneOperator2_<double, KN<double> * const, dualmesh * const>  (Dual_Mesh_Nbdryvert));
    Global.Add("dualmeshmidedge", "(", new OneOperator2_<double, KN<double> * const, dualmesh * const>  (Dual_Mesh_Midedge));
}
