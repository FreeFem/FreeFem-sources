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
// Under Linux (with FreeFEM++), compile into a plugin with
// ff-c++ -shared -fPIC FreeVol_advectiondm.cpp -I$FREEFEMPP_DIR/include -L$FREEFEMPP_DIR/lib FreeVol.so
// or simply ff-c++ -shared -fPIC FreeVol_advectiondm.cpp FreeVol.so
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...
/* ---------------------------------------------------------

   finite volume method to solve the advection equation
   U_t + div([V1*U,V2*U]) = 0, U scalar,
   where V1,V2 depend on space position.
   second-order in space by slope reconstruction
   numflux = 0 : HLL, numflux = 1 : Rusanov, numflux = 2 : upwind
   unknown U is P1 (one value by vertex), the control volumes are the cells of the dual mesh of the mesh of the
   velocity, built in the script and given to fvscheme; the method is given by the type of the dual mesh:
   1D (meshL):  meshdualL Thstar = builddualmesh(Th), the median dual mesh, fvscheme FV(Thstar, [V1], eq="advectiondm", ...),
                the velocity in P1
   2D (mesh):   meshdual Thstar = builddualmesh(Th), the dual mesh of F. Bouchut, FV(Thstar, [V1,V2], ...), the
                velocity in P0, P1, or a fespace with one dof by edge (P0edge or P1nc)
   3D (mesh3):  meshdual3 Thstar = builddualmesh(Th), the median dual mesh, FV(Thstar, [V1,V2,V3], ...), the velocity in P1
   surface (meshS): meshdualS Thstar = builddualmesh(Th), the median dual mesh, FV(Thstar, [V1,V2,V3], ...), the velocity in P1 (tangent to the surface)

   --------------------------------------------------------*/

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include "FreeVolDualMesh.hpp"  // FreeVol.hpp, the dual meshes (fvdualmesh.hpp), the base class FreeVol_DualModel

// advection by a given vector field over the dual mesh (P1 unknown), 1D (meshL), 2D (mesh), 3D (mesh3)
// the velocity: 2D (the dual mesh of F. Bouchut) in P0, P1 (the exact mean on the interfaces), or a fespace with
// one dof by edge (P0edge, P1nc); 1D, 3D (the median dual mesh) in P1
// the slopes in the control volumes, limited, the values reconstructed at the centroids of the interfaces
class FreeVol_AdvDM : public FreeVol_DualModel
{
public:
    int nV;             // the number of components of the velocity: the dimension
    mutable KN_d V[3];  // the components of the velocity (their values: read again at each fluxupdate)
    vector< FreeVol_FE > Vf;  // the FE functions of the velocity, given to fvscheme
    long nnV;           // the number of dof by component of the velocity
    int vtype;          // the fespace of the velocity: 0 P0, 1 P1, 2 one dof by edge (P0edge, P1nc; 2D only)
    // 2D (the dual mesh of F. Bouchut)
    KN<long> dofV[2];   // vtype 2: dofV[c][3*k+e] the dof of V[c] on the edge e of the triangle k
    KN<long> k3eind; // one of the two indices k3e of the dual edge considered
    // P1 velocity on an interior interface: the exact mean of the P1 velocity on the interface (the segment between
    // the 2 centers is cut where it crosses the edge, the velocity is linear on each part): exact balance of a
    // constant state (|C_i| div V for a linear velocity)
    KN<long> ipk;       // P1: [2m+p] the triangle of the part p of the interface m (-1: no part)
    KN<double> ipl;     // P1: [2m+p] the length of the part
    KN<double> ipx;     // P1: [4m+2p+c] the middle of the part
    mutable KN<double> VIx, VIy;// P1: the mean of the velocity on each interior interface

    FreeVol_AdvDM(const dualmesh &dm, const vector< FreeVol_FE > &Vfe,int oord,int snumflux,KN<long> labelBCC,int limiterr,KN<double> parameterss,KN<long> intparameterss)
    : FreeVol_DualModel(dm), nV(0), nnV(0), vtype(1)
    { init(Vfe,oord,snumflux,labelBCC,limiterr,parameterss,intparameterss); }
    FreeVol_AdvDM(const dualmeshN &dm, const vector< FreeVol_FE > &Vfe,int oord,int snumflux,KN<long> labelBCC,int limiterr,KN<double> parameterss,KN<long> intparameterss)
    : FreeVol_DualModel(dm), nV(0), nnV(0), vtype(1)
    { init(Vfe,oord,snumflux,labelBCC,limiterr,parameterss,intparameterss); }

    void init(const vector< FreeVol_FE > &Vfe,int oord,int snumflux,KN<long> labelBCC,int limiterr,KN<double> parameterss,KN<long> intparameterss){
        int fe = Vfe[0].dim;// 1 meshL, 2 mesh, 3 mesh3, 4 meshS
        int d = (fe == 4) ? 3 : fe;// surface (meshS): 3 components of the velocity, tangent to the surface
        nV = Vfe.size();
        if (fe < 1 || fe > 4 || nV != d)
            ExecError("fvscheme advectiondm: the velocity has 1 component on a meshL (1D), 2 on a mesh (2D), 3 on a mesh3 (3D) and on a meshS (surface)");
        for (int c = 0; c < nV; ++c) {
            if (!Vfe[c].x()) { cout << " Sorry V" << c+1 << " is unset "<< endl; ffassert(0);}
            if (Vfe[c].dim != fe || Vfe[c].mesh() != Vfe[0].mesh()) ExecError("fvscheme advectiondm: the components of the velocity must be on the same mesh");
        }
        int neq = 1;//defines the number of equations (unknowns) of the system
        set(neq,oord,snumflux,Vfe[0],labelBCC,"advectiondm",limiterr,parameterss,intparameterss);
        CheckDual("fvscheme advectiondm");
        Vf = Vfe;
        if (!dual2D()) {// 1D, 3D, surface: P1
            nnV = geom.nv;
            for (int c = 0; c < nV; ++c) {
                const FreeVol_FE &f = Vfe[c];
                if (f.NbOfDF() / f.N() != geom.nv || f.MaxNbDFPerElement() / f.N() != (fe == 4 ? 3 : d+1))
                    ExecError("fvscheme advectiondm: in 1D, 3D and on a surface, the velocity must be P1");
                V[c].set( f.values() );
            }
            return;
        }
        // 2D: the dual mesh of F. Bouchut
        const FreeVol_Geom &g = geom;
        const Mesh &Th = *(const Mesh *)g.pTh;
        const dualmesh &Thstar = *dmF;
        for(int c=0; c<2; ++c){
            const FreeVol_FE &f = Vfe[c];
            long nn  = f.NbOfDF() / f.N();                 // nb dofs per component
            int ldf = f.MaxNbDFPerElement() / f.N();       // 1 : P0, 3 : P1 or one dof by edge
            if((verbosity>4)&&(!mpirank)) cout << "comp=" << f.comp() << " N=" << f.N() << " nn=" << nn << " ldf=" << ldf << endl;
            int t;
            if (ldf == 1 && nn == Th.nt) t = 0;                     // P0
            else if (ldf == 3 && nn == Th.nv) t = 1;                // P1
            else if (ldf == 3 && nn == g.nfaces) {                  // one dof by edge: P0edge, P1nc
                t = 2;
                dofV[c].resize(3*Th.nt);
                for(long k=0; k<Th.nt; ++k)
                    for(int e=0; e<3; ++e) dofV[c][3*k+e] = f.dof(k,e);
            }
            else {cout << "fespace must be built on P0 or P1 or a fespace with one dof by edge (P0edge, P1nc)" << endl;ffassert(0);}
            if (c == 0) vtype = t;
            else if (t != vtype) {cout << " error: the components of V must be of the same type" << endl;ffassert(0);}
            V[c].set( f.values() );
            if (c == 0) nnV = nn;
        }
        Vf = Vfe;

        const long npe=DM().npe, nde=DM().nde;
        const KN<long> &edgenumbering=*(Thstar.edgenumbering);// (size 3*Th.nt) index of primal edges, in terms of triangles and their 3 edges (each triangle has 3 adjacent values in array)
        const KNM<long> &dualpiece01=*(Thstar.dualpiece01);// (size (2,ndualedge)) origin and end of piece of boundary of control volume, in terms of dual edge index (value is dual vertex index, origin and end indexed by the first int (0 or 1))
        k3eind.resize(nde);
        for (long m=0;m<npe;m++){//loop in internal dual edges
         long k=dualpiece01(1L,m);//triangle
         for (long e=0;e<3;++e){
          if (edgenumbering(3*k+e)==m) k3eind(m)=3*k+e;
         }
        }
        long countbe=0;
        for (long ib=0;ib<Th.neb;++ib){//loop on the boundary elements of primal mesh
            int iib;
            long kib = Th.BoundaryElement(ib,iib);
            int ee=iib;
            long kee=Th.ElementAdj(kib,ee);
            if(int(kee)==kib) {//really on boundary
             k3eind(npe+2*countbe)=3*kib+iib;
             k3eind(npe+2*countbe+1)=3*kib+iib;
             countbe++;
            }
        }
        // P1 velocity: the parts of the interior interfaces in the triangles
        if (vtype == 1) {
            const KN<double> &dvxy=*((Thstar.normalized)->dvxy);// the dual vertices: the centers of the triangles, then the points of the boundary
            const KN<long> &cvLR=*(Thstar.cvLR);
            ipk.resize(2*npe); ipl.resize(2*npe); ipx.resize(4*npe); VIx.resize(npe); VIy.resize(npe);
            ipk=-1L; ipl=0.; ipx=0.;
            for (long m=0; m<npe; ++m) {
                long a=dualpiece01(0L,m), b=dualpiece01(1L,m);
                R2 A(dvxy[2*a],dvxy[2*a+1]), B(dvxy[2*b],dvxy[2*b+1]);
                long ka = a < Th.nt ? a : -1, kb = b < Th.nt ? b : -1;
                R2 P[2][2]; long kp[2]; int np=0;
                if (ka >= 0 && kb >= 0) {// the crossing point C of [A,B] with the edge (cvLR) of the 2 triangles
                    R2 E0=Th(cvLR[2*m]), E1=Th(cvLR[2*m+1]), d=B-A, e=E1-E0, w=E0-A;
                    double det = -d.x*e.y + e.x*d.y;
                    double sc = fabs(det) > 1e-30 ? (-w.x*e.y + e.x*w.y)/det : 0.5;
                    sc = max(0., min(1., sc));
                    R2 C = A + d*sc;
                    kp[0]=ka; P[0][0]=A; P[0][1]=C; kp[1]=kb; P[1][0]=C; P[1][1]=B; np=2;
                }
                else { kp[0] = ka >= 0 ? ka : kb; P[0][0]=A; P[0][1]=B; np=1; }// the interface in one triangle
                for (int p=0; p<np; ++p) {
                    R2 M=(P[p][0]+P[p][1])*0.5, D=P[p][1]-P[p][0];
                    ipk[2*m+p]=kp[p]; ipl[2*m+p]=sqrt(D.x*D.x+D.y*D.y);
                    ipx[4*m+2*p]=M.x; ipx[4*m+2*p+1]=M.y;
                }
            }
        }
    }//end of init
    // P1 velocity: its value at the point P of the triangle k (barycentric coordinates)
    R2 VP1(long k, double px, double py) const {
        const Mesh &Th = *dmF->primal;
        R2 A=Th[k][0], B=Th[k][1], C=Th[k][2];
        double det=(B.x-A.x)*(C.y-A.y)-(B.y-A.y)*(C.x-A.x);
        double l1=((px-A.x)*(C.y-A.y)-(py-A.y)*(C.x-A.x))/det, l2=((B.x-A.x)*(py-A.y)-(B.y-A.y)*(px-A.x))/det, l0=1.-l1-l2;
        long i0=Th(k,0), i1=Th(k,1), i2=Th(k,2);
        return R2(l0*V[0][i0]+l1*V[0][i1]+l2*V[0][i2], l0*V[1][i0]+l1*V[1][i1]+l2*V[1][i2]);
    }
    // P1 velocity: its mean on each interior interface (the velocity can change at each fluxupdate)
    void UpdateVI() const {
        for (long m=0; m<DM().npe; ++m) {
            R2 S(0.,0.); double L=0.;
            for (int p=0; p<2; ++p) if (ipk[2*m+p] >= 0) {
                S += VP1(ipk[2*m+p], ipx[4*m+2*p], ipx[4*m+2*p+1]) * ipl[2*m+p];
                L += ipl[2*m+p];
            }
            VIx[m] = S.x/L; VIy[m] = S.y/L;
        }
    }

    // the values of the velocity: read again before each computation, the velocity can be changed after
    // fvscheme (V1 = ... in the script makes a new array of values: the old one is freed)
    void UpdateV() const {
        for(int c=0; c<nV; ++c){
            if (!Vf[c].x() || Vf[c].NbOfDF() / Vf[c].N() != nnV)
                ExecError("fluxupdate advectiondm: the velocity is not on its fespace of the definition of fvscheme (resetfvscheme)");
            V[c].set( Vf[c].values() );
        }
    }
    double VN2(long m2LR,const R3 &N) const {//2D: normal velocity at each side of edge
        const Mesh &Th = *dmF->primal;
        const long npe=DM().npe;
        const KN<long> &cvLR=*(dmF->cvLR);//(size 2*nde) indices of Left and Right control volumes, in terms of dual edge index (Left and Right are written adjacently, the list is by couple)
        const KNM<long> &dualpiece01=*(dmF->dualpiece01);// (size (2,ndualedge)) origin and end of piece of boundary of control volume, in terms of dual edge index (value is dual vertex index, origin and end indexed by the first int (0 or 1))
        long m=m2LR/2;//index of dual edge
        R2 Vm;//velocity on the edge
        if (vtype==0) {             // P0
          long k0,k1;//neighboring triangles
          if (m<npe) {
           k0=dualpiece01(0L,m);
           k1=dualpiece01(1L,m);
           if (k0>Th.nt) k0=k1;
          } else {
           k1=k3eind[m]/3;
           k0=k1;
          }
          Vm=R2( .5*(V[0][k0]+V[0][k1]), .5*(V[1][k0]+V[1][k1]) );
        } else if (vtype==1) {      // P1
          long iv=cvLR(m2LR);//index of control volume
          if (m<npe) {              // the exact mean of the velocity on the interface
           Vm=R2( VIx[m], VIy[m] );
          } else {
           long ivR=((m-npe)%2==0)? cvLR(2*(m+1)) : cvLR(2*(m-1));//index of the other primal vertex on the boundary
           Vm=R2( 0.75*V[0][iv]+0.25*V[0][ivR], 0.75*V[1][iv]+0.25*V[1][ivR] );
          }
        } else {                    // one dof by edge
          long k3e=k3eind[m];
          Vm=R2( V[0][dofV[0][k3e]], V[1][dofV[1][k3e]] );
        }
        return Vm.x*N.x+Vm.y*N.y;
    }
    double VNN(long m2LR,const R3 &N) const {//1D, 3D: normal velocity at each side of the interface
        const dualmeshN &D = DM();
        long m=m2LR/2;//index of the interface
        R3 Vm(0.,0.,0.);
        if (m<D.npe) {// interior: the velocity at the vertex of the side
            long iv=D.cvLR[m2LR];
            for(int c=0; c<nV; ++c) Vm[c]=V[c][iv];
        } else {// boundary piece: the velocity consistent with the interior interfaces (mean of the 2 vertices =
                // value at the middle of the edge): 1D the vertex, 3D 1/2 the vertex + 1/4 the middles of its 2 edges
                // of the boundary triangle (weights 3/4, 1/8, 1/8; not the centroid of the piece, 11/18, 7/36, 7/36):
                // then the balance of a constant state is exactly |C_i| div V for a P1 velocity, also at the boundary
            long b=3*(m-D.npe);
            static const double w3[3]={.75,.125,.125};
            for(int c=0; c<nV; ++c)
                if (nV==1) Vm[c]=V[c][D.bvert[b]];
                else for(int l=0; l<3; ++l) Vm[c]+=w3[l]*V[c][D.bvert[b+l]];
        }
        return (Vm,N);
    }
    double VN(long m2LR,const R3 &N) const { return dual2D() ? VN2(m2LR,N) : VNN(m2LR,N); }
    void Flux(long m2LR,const R3 &N,const double *U,double *F) const {
        double vn= VN(m2LR,N);
        F[0]= vn*U[0];
    }
    void lambda(long m2LR,const R3 &N,const  double *U,double *l) const {
        double vn= VN(m2LR,N);
        l[0]= vn;
    }
    // numerical flux of an interior face: numflux = 0 HLL, 1 Rusanov, 2 upwind
    void NFlux(long m2LRl,long m2LRr,const R3 &Nn,const double *Ul,const double *Ur,double &speed,double *phi) const {
        if(numflux != 2) {// 0 HLL, 1 Rusanov
            double Fr[Neqmax],Fl[Neqmax];
            if(numflux == 0){// HLL
                double sl = min(lambdamin(m2LRl,Nn,Ul),lambdamin(m2LRr,Nn,Ur));
                double sr = max(lambdamax(m2LRl,Nn,Ul),lambdamax(m2LRr,Nn,Ur));
                speed=max(abs(sl),abs(sr));
                Flux(m2LRl,Nn,Ul,Fl);
                Flux(m2LRr,Nn,Ur,Fr);
                for(int j=0; j<Neq;++j)
                    phi[j] = (sl>=0) ? Fl[j] : ((sr<=0) ? Fr[j] : (sr*Fl[j] - sl*Fr[j] + sl*sr*(Ur[j]-Ul[j]))/(sr-sl+1.e-10));
            }
            else if(numflux == 1){// Rusanov
                double cRus = max(absvp(m2LRl,Nn,Ul),absvp(m2LRr,Nn,Ur));
                speed=cRus;
                Flux(m2LRl,Nn,Ul,Fl);
                Flux(m2LRr,Nn,Ur,Fr);
                for(int j=0; j<Neq;++j)
                    phi[j] = 0.5*((Fr[j]+Fl[j]) - cRus*(Ur[j]-Ul[j]));
            }
            else ExecError("FreeVol advectiondm: numflux not defined");
            return;
        }
        double lL[Neqmax],lR[Neqmax];
        lambda(m2LRl,Nn,Ul,lL);
        lambda(m2LRr,Nn,Ur,lR);
        double vit=.5*(lL[0]+lR[0]);
        speed=abs(vit);
        if (vit>0)
            phi[0]=vit*Ul[0];
        else
            phi[0]=vit*Ur[0];
    }
    // numerical flux of a boundary face: labelBC = 0 Dirichlet (uBC), 1 Neumann (the default)
    void NFlux_BC(long kf,int labelement,const R3 &Nn,const double *UL,double &speed,double *phi, const double *UBCL) const {
        int choice_BC = 1;// 0 Dirichlet, 1 Neumann(default)
        for(int i=0;i<labelBC.n/2;i++)
            if(labelement==labelBC[2*i]) choice_BC = labelBC[2*i+1];
        if(choice_BC == 0) NFlux(kf,kf,Nn,UL,UBCL,speed,phi);
        else if(choice_BC == 1) NFlux(kf,kf,Nn,UL,UL,speed,phi);
        else {cout<<"choice_BC = "<<choice_BC<<" is not defined for advectiondm"<<endl;ExecError("FreeVol advectiondm: labelBC not defined");}
    }

    // the scheme on the dual mesh: the flux balance of each control volume and the time step (CFL = 1)
    bool PeriodicOK() const { return true; }// periodic= (PSkip, PeriodicU, PeriodicUp, FreeVolDualMesh.hpp)
    double computedtFlux(KN_d *U,KN_d *Up, KN_d *uBC) const {
        CheckDual("fluxupdate advectiondm");// first: UpdateVI uses the dual mesh
        PeriodicU(U);// periodic: the same values in the control volumes of a class
        UpdateV();
        if (dual2D() && vtype==1) UpdateVI();
        const dualmeshN &D = DM();
        const long nv=D.nv, npe=D.npe, nde=D.nde;
        double maxspeed;
        Up[0]=0.;
        //left and right states and fluxes for an interface
        double Ul[Neqmax],Ur[Neqmax],phiel[Neqmax],UBCl[Neqmax];
        vector< KN<double> > gU(1);// the slopes in the control volumes
        if(ord>1){
            KN<double> u(nv), alpha_u(nv);
            u=U[0];
            Grad(u,gU[0]);// centered slopes in the control volumes
            BdrySlope(gU[0],1L);
            SlopeLimiterVF(u,gU[0],alpha_u,limiter);
            for(long i=0;i<nv;++i)
                for(int c=0;c<3;++c) gU[0][3*i+c] *= alpha_u[i];
        }
        KN<long> cut;// the control volumes cut by the border of the local mesh (domain decomposition): not in dt
        DdmCutCV(cut);
        KN<double> sumf(dual2D() ? 0 : nv);// 1D, 3D: sum on the interfaces of the control volume of area*speed
        if (!dual2D()) sumf=0.;
        double dt=1.e4;
        for (long m=0;m<nde;++m){//loop on the interfaces of the dual mesh
            if (PSkip(m)) continue;// periodic: a boundary piece inside a control volume (no flux)
            long ivL=D.cvLR[2*m],ivR=D.cvLR[2*m+1];
            R3 N(D.normal[3*m],D.normal[3*m+1],D.normal[3*m+2]);
            GetU(m,ivL,U,&gU[0],Ul);
            if (m<npe) {//if interface is internal
                GetU(m,ivR,U,&gU[0],Ur);
                NFlux(2*m,2*m+1,N,Ul,Ur,maxspeed,phiel);// numerical flux
            } else {//boundary piece: uBC (P1) at its centroid
                long b=3*(m-npe);
                for(int j=0;j<Neq;++j){
                    UBCl[j]=0.;
                    for(int l=0;l<3;++l) UBCl[j]+=D.bweight[b+l]*uBC[j][D.bvert[b+l]];
                }
                NFlux_BC(2*m,D.dbelab[m-npe],N,Ul,maxspeed,phiel,UBCl);
            }
            double le=D.area[m];
            if (dual2D()) dt=DtEdge(dt,m,npe,ivL,ivR,cut,D.cvvol,le,maxspeed);// 2D: by dual edge
            else { sumf[ivL]+=le*maxspeed; if(m<npe) sumf[ivR]+=le*maxspeed; }
            // distribution of the fluxes to the control volumes
            for(int j=0; j<Neq;++j){
                phiel[j] *= le;
                Up[j][ivL] += phiel[j]*(1./D.cvvol[ivL]);
                if(m<npe) {//if interface is internal
                    Up[j][ivR] -= phiel[j]*(1./D.cvvol[ivR]);
                }
            }
        }
        if (!dual2D())// 1D, 3D: min on the control volumes of volume/(sum on its interfaces of area*speed)
            for(long i=0;i<nv;++i) if(!cut[i]) dt=min(dt,D.cvvol[i]/(sumf[i]+1.e-20));
        PeriodicUp(Up);// periodic: the flux balance of the union of the control volumes of a class
        return dt;
    }//end of computedtFlux

    ~FreeVol_AdvDM() {if(verbosity>99)  cout << " delete FreeVol_AdvDM "<< this << endl;}
};

// creation of the model for fvscheme FV([V1,V2], eq="advectiondm", ...) (see the registry in FreeVol.cpp)
static FreeVol_model *FreeVol_AdvDM_Create(const vector< FreeVol_FE > &V, long order, long numflux,
                                           KN<long> labelBC, long limiter, KN<double> parameters,
                                           KN<long> intparameters,
                                           const void *dual) {
    if (numflux == -99) numflux = 2;    // upwind by default
    if ((verbosity) && (!mpirank)) {
        if ((parameters.N() > 1) || ((parameters.N() == 1) && parameters[0] != 0.)) {    // parameters is not yet used
            cout << "parameters = [ ";
            for (int i = 0; i < parameters.N() - 1; i++) cout << parameters[i] << " , ";
            cout << parameters[parameters.N() - 1] << " ]" << endl;
        }
        if ((intparameters.N() > 1) || ((intparameters.N() == 1) && intparameters[0] != 0)) {    // intparameters is not yet used
            cout << "intparameters = [ ";
            for (int i = 0; i < intparameters.N() - 1; i++) cout << intparameters[i] << " , ";
            cout << intparameters[intparameters.N() - 1] << " ]" << endl;
        }
        for (int i = 0; i < labelBC.n / 2; i++) {
            if (labelBC[2 * i + 1] == 0) cout << "on boundary " << labelBC[2 * i] << " labelBC=0 : Dirichlet BC" << endl;
            else if (labelBC[2 * i + 1] == 1) cout << "on boundary " << labelBC[2 * i] << " labelBC=1 : Neumann BC" << endl;
            else cout << "on boundary " << labelBC[2 * i] << " labelBC=" << labelBC[2 * i + 1] << " is not defined" << endl;
        }
        if (numflux == 0) cout << "numflux = 0 : HLL" << endl;
        else if (numflux == 1) cout << "numflux = 1 : Rusanov" << endl;
        else if (numflux == 2) cout << "numflux = 2 : Upwind" << endl;
        if (limiter == 0) cout << "limiter = 0 : no limiter is used" << endl;
        else if (limiter == 1) cout << "limiter = 1 : maximum principle type" << endl;
        else if (limiter == 2) cout << "limiter = 2 : minmod" << endl;
        else if (limiter == 3) cout << "limiter = 3 : Sweby" << endl;
        else if (limiter == 4) cout << "limiter = 4 : Van Albada type" << endl;
        else cout << "limiter = " << limiter << endl;
    }
    // a numerical flux not defined for the advection: stop (whatever the verbosity)
    if (numflux < 0 || numflux > 2) {
        cout << "numflux = " << numflux << " is not defined for advectiondm (0 HLL, 1 Rusanov, 2 Upwind)" << endl;
        ExecError("FreeVol advectiondm: numflux not defined");
    }
    // the dual mesh built in the script: a dualmeshL* (1D), a dualmesh* (2D), a dualmesh3* (3D), given to fvscheme
    // FV(Thstar, [V1,V2], ...), or found from the mesh of [V1,V2]: fvscheme FV([V1,V2], ...)
    dual = FreeVol_DualModel::FindDual(dual, V, "advectiondm", "[V1,V2]");
    // the dual mesh given to fvscheme: 2D the dual mesh of F. Bouchut (meshdual), 1D, 3D the median dual mesh
    // (meshdualL, meshdual3)
    int d = V.size() ? V[0].dim : 2;
    if (d == 2)
        return new FreeVol_AdvDM(*(const dualmesh *)dual, V, order, numflux, labelBC, limiter, parameters, intparameters);
    return new FreeVol_AdvDM(*FreeVol_DualModel::DualN(dual, d), V, order, numflux, labelBC, limiter, parameters, intparameters);
}

static void Load_Init() {
    if (verbosity && !mpirank) cout << " load: FreeVol_advectiondm" << endl;
    FreeVol_RegisterModel("advectiondm", FreeVol_AdvDM_Create);
}
LOADFUNC(Load_Init)
