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
// ff-c++ -shared -fPIC FreeVol_swtopo.cpp -I$FREEFEMPP_DIR/include -L$FREEFEMPP_DIR/lib FreeVol.so
// or simply ff-c++ -shared -fPIC FreeVol_swtopo.cpp FreeVol.so
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...
/* ---------------------------------------------------------

   finite volume method to solve the shallow water system with topography
   U_t + div(F(U)) + gravity*H*(0,grad zb) = 0, U=[H,H*u], H>=0, u the velocity (d components)
   F(U) = ( H*u , H*u x u + gravity/2*H^2 Id ), zb is a given space dependent topography, zb_t=0
   in 1D, 2D or 3D, the dimension d is given by the mesh of the unknowns:
     1D: meshL (segments on the x axis), fvscheme FV([H,Hu,zb], eq="swtopo", ...)
     2D: mesh (triangles),               fvscheme FV([H,Hu,Hv,zb], eq="swtopo", ...)
     3D: mesh3 (tetrahedra),             fvscheme FV([H,Hu,Hv,Hw,zb], eq="swtopo", ...)
     (in 3D, it is the barotropic gas dynamics with p = gravity/2*H^2 and a potential force)
     surface: meshS (triangles in 3D),   fvscheme FV([H,Hu,Hv,Hw,zb], eq="swtopo", ...), the momentum tangent to
       the surface (gravity normal to the surface, zb along the normal): the momentum of each cell is in the plane
       of its triangle; on an edge, the state of the neighbour is turned about the edge into the plane of the cell,
       and the momentum flux is turned back into the plane of the neighbour (each triangle is flat: the lake at
       rest is exact as in 2D)
   gravity=parameters[0]  --> parameters=[gravity]
   second-order in space by slope reconstruction (on H, H+zb and u)
   hydrostatic reconstruction is used to handle topography
   numflux = 0 : HLL, numflux = 1 : Rusanov, numflux = 2 : Suliciu/HLLC
   labelBC = 0 : Dirichlet, 1 : Neumann, 2 : Wall, 3 : mix Dirichlet and wall
   unknown (U,zb) is P0 (one value by cell)

   --------------------------------------------------------*/

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include "FreeVol.hpp"

// shallow water system with topography, in 1D, 2D or 3D
// the unknowns: U[0] = H, U[1..d] = H*u (d components), U[d+1] = zb
class FreeVol_SWTOPO: public FreeVol_model
{
public:
    int d;     // dimension, d components of the velocity
    int iz;    // the index of zb in U: d+1

    FreeVol_SWTOPO(const vector< FreeVol_FE > &Vfe,int oord,int snumflux,KN<long> labelBCC,int limiterr,KN<double> parameterss,KN<long> intparameterss)
    : FreeVol_model(), d(0), iz(0)
    { init(Vfe,oord,snumflux,labelBCC,limiterr,parameterss,intparameterss);}

    void init(const vector< FreeVol_FE > &Vfe,int oord,int snumflux,KN<long> labelBCC,int limiterr,KN<double> parameterss,KN<long> intparameterss){
        int nU = Vfe.size();
        d = nU ? Vfe[0].dim : 0;
        if (d == 4) d = 3;// surface (meshS): 3 components of the momentum, tangent to the surface
        if (d < 1 || nU != d+2)
            ExecError("fvscheme swtopo: the unknowns are [H,Hu,zb] in 1D (meshL), [H,Hu,Hv,zb] in 2D (mesh), [H,Hu,Hv,Hw,zb] in 3D (mesh3) and on a surface (meshS)");
        for (int c = 0; c < nU; ++c)
            if (!Vfe[c].x()) { cout << " Sorry unknown " << c+1 << " is unset "<< endl; ffassert(0);}
        const FreeVol_FE &f = Vfe[0];
        if (f.MaxNbDFPerElement() / f.N() != 1) ExecError("fvscheme swtopo: the unknowns must be P0");
        iz = d+1;
        int neq = d+2;//defines the number of equations (unknowns) of the system
        set(neq,oord,snumflux,f,labelBCC,"swtopo",limiterr,parameterss,intparameterss);
        if (f.NbOfNodes() != geom.nt) ExecError("fvscheme swtopo: the unknowns must be P0");
    }//end of init

    // Flux F(U).N with included topo
    void Flux(long kf,const R3 &N,const double *U,double *F) const {
        double h  = U[0];
        for (int j = 0; j < Neq; ++j) F[j] = 0.;// no flux for the topography
        if(h <= 1.e-13) return;
        double un = 0.;
        for (int c = 0; c < d; ++c) un += U[1+c]/h*N[c];
        F[0] = h*un;
        for (int c = 0; c < d; ++c) F[1+c] = U[1+c]*un + 0.5*parameters[0]*h*h*N[c];
    }
    // Local eigenvalues: un-c, un (d-1 times), un+c, 0 (topography)
    void lambda(long kf,const R3 &N,const double *U,double *l) const {
        double h  = U[0];
        for (int j = 0; j < Neq; ++j) l[j] = 0.;
        if(h <= 1.e-13) return;
        double un = 0.;
        for (int c = 0; c < d; ++c) un += U[1+c]/h*N[c];
        double cs = sqrt(parameters[0]*h);
        l[0] = un - cs;
        for (int c = 1; c < d; ++c) l[c] = un;
        l[d] = un + cs;
        l[iz] = 0.; // topographie
    }
    // the limited slopes: on H, H+zb and u, then limitation for H positive and for the velocity
    void SWSlopes(KN_d *U, KN<double> *gU) const {
        const FreeVol_Geom &g = geom;
        long nt = g.nt;
        KN<double> h(nt),hpzb(nt),vel[3],alpha(nt);
        KN<double> grad_h(d*nt),grad_hpzb(d*nt),grad_vel[3];
        for(int c=0; c<d; ++c) { vel[c].resize(nt); grad_vel[c].resize(d*nt); }
        for(long k=0; k< nt; ++k){
            h[k] = U[0][k];
            for(int c=0; c<d; ++c) vel[c][k] = U[1+c][k]/max(U[0][k],1.e-13);
        }
        hpzb=h+U[iz];
        MatriceMorse<R> &MATGRAD = *(matgrad);
        grad_h = MATGRAD * U[0];
        grad_hpzb = MATGRAD * hpzb;
        for(int c=0; c<d; ++c) grad_vel[c] = MATGRAD * vel[c];
        BdrySlope(grad_h,1L);
        BdrySlope(grad_hpzb,1L);
        for(int c=0; c<d; ++c) BdrySlope(grad_vel[c],1L);
        SlopeLimiterVF(h,grad_h,alpha,limiter);
        for(long k=0; k< d*nt; ++k) grad_h[k] *= alpha[k/d];
        SlopeLimiterVF(hpzb,grad_hpzb,alpha,limiter);
        for(long k=0; k< d*nt; ++k) grad_hpzb[k] *= alpha[k/d];
        for(int c=0; c<d; ++c) {
            SlopeLimiterVF(vel[c],grad_vel[c],alpha,limiter);
            for(long k=0; k< d*nt; ++k) grad_vel[c][k] *= alpha[k/d];
        }
        for(long k=0; k< nt; ++k){
            if (U[0][k]<1.e-13) {// dry cell
                for(int i=0; i<d; ++i) {
                    grad_h[d*k+i] = 0.;
                    for(int c=0; c<d; ++c) grad_vel[c][d*k+i] = 0.;
                }
            }
        }
        for(int j=0; j<Neq; ++j) gU[j].resize(d*nt);
        for(long k=0; k< nt; ++k){
            double hk = U[0][k];
            // limiter for H positive
            double alphah = 1.;
            for (int e=0;e<g.nf;++e) {
                R3 GE = g.F[k*g.nf+e]-g.G[k];// centroid of k to the centroid of face e
                double he = hk;
                for(int i=0; i<d; ++i) he += grad_h[d*k+i]*GE[i];
                alphah = min(alphah,1./max((1.-he/max(hk,1.e-13))/0.9,1.));
            }
            for(int i=0; i<d; ++i) grad_h[d*k+i] *= alphah;
            // limitation for "positivity of internal energy"
            double maxvellim = 0.;
            for (int e=0;e<g.nf;++e) {
                R3 GE = g.F[k*g.nf+e]-g.G[k];
                double he = hk;
                for(int i=0; i<d; ++i) he += grad_h[d*k+i]*GE[i];
                double inc2 = 0.;
                for(int c=0; c<d; ++c) {
                    double incue = 0.;
                    for(int i=0; i<d; ++i) incue += grad_vel[c][d*k+i]*GE[i];
                    incue *= hk;
                    incue /= max(he,1.e-13);
                    inc2 += incue*incue;
                }
                maxvellim = max(maxvellim,sqrt(.5*inc2/max(.5*parameters[0]*he,1.e-30)));
            }
            double alphavel = 1./max(maxvellim/.5,1.);
            //define the final slopes
            for(int i=0; i<d; ++i) {
                long ki = d*k+i;
                for(int c=0; c<d; ++c) grad_vel[c][ki] *= alphavel;
                gU[0][ki] = grad_h[ki];
                for(int c=0; c<d; ++c) gU[1+c][ki] = grad_h[ki]*vel[c][k]+grad_vel[c][ki]*hk;
                gU[iz][ki] = grad_hpzb[ki]-grad_h[ki];
            }
        }
    }
    // hydrostatic reconstruction: the state at the face with the topography of the other side
    void HydroRec(const double *U, double zbo, double *Us) const {
        double hs = max(U[0] - max(zbo - U[iz],0.),0.);
        Us[0] = hs;
        for(int c=0; c<d; ++c) Us[1+c] = U[0]>0. ? U[1+c]/U[0]*hs : 0.;
        Us[iz] = U[iz];
    }
    // surface: the momentum q (3 components) in the plane of the cell k (its component on the normal removed)
    void ProjPlane(long k, double *q) const {
        const R3 &n = geom.Nc[k];
        double qn = q[0]*n.x + q[1]*n.y + q[2]*n.z;
        for(int c=0; c<3; ++c) q[c] -= qn*n[c];
    }
    // surface: the rotation about the edge of the face kf (of k, kfo the same face in its neighbour ko) which turns
    // the plane of ko into the plane of k (dir = 1) or the plane of k into the plane of ko (dir = -1): the frame
    // (N of k, edge, N of k x edge) of k corresponds to the frame (-N of ko, edge, -N of ko x edge) of ko
    void RotEdge(long kf, long kfo, double *q, int dir) const {
        const FreeVol_Geom &g = geom;
        const MeshS &Th = *(const MeshS *)g.pTh;
        R3 A = Th(g.fvert[2*kf]), B = Th(g.fvert[2*kf+1]), t = B - A;// the direction of the edge
        t /= t.norme();
        R3 a1 = g.N[kf], a2 = -g.N[kfo], b1 = a1^t, b2 = a2^t;
        R3 v(q[0],q[1],q[2]), w;
        if(dir > 0) w = a1*(v,a2) + t*(v,t) + b1*(v,b2);
        else w = a2*(v,a1) + t*(v,t) + b2*(v,b1);
        for(int c=0; c<3; ++c) q[c] = w[c];
    }
    double computedtFlux(KN_d *U,KN_d *Up, KN_d *uBC) const {
        const FreeVol_Geom &g = geom;
        double gr = parameters[0];
        double dt=1.e4,maxspeed;
        KN<double> sumf(g.dim==2 ? 0 : g.nt);// 1D, 3D: sum on the faces of the cell of fmeas*speed
        // domain decomposition: not in the time step the speeds which are not the ones of the whole mesh (overlap
        // cells, set again by the exchange; each face and cell is complete on one process at least: overlap of 1 cell
        // at the order 1, 2 at the order 2):
        //   cut: the cells with a face on the border of the local mesh (label -111111), the speed of this face is not
        //        the one of the true neighbour, at the order 2 their slopes are not complete
        //   2D (dt by face): the faces -111111, at the order 2 the faces of a cut cell
        //   1D, 3D (dt by cell, sum on its faces): the cut cells, at the order 2 also their neighbours (dtcut)
        KN<long> cut(g.nt), dtcut(g.nt); cut=0L;
        for(long kf=0; kf<g.nt*g.nf; ++kf) if(g.adj[kf]<0 && g.flab[kf]==-111111) cut[kf/g.nf]=1;
        dtcut=cut;
        if(ord>1)
            for(long kf=0; kf<g.nt*g.nf; ++kf) if(g.adj[kf]>=0 && cut[g.adj[kf]]) dtcut[kf/g.nf]=1;
        if(g.dim!=2) sumf=0.;
        for(int j=0;j<Neq;++j) Up[j]=0.;
        //left and right states and fluxes for an interface
        double Ul[Neqmax],Ur[Neqmax],Uls[Neqmax],Urs[Neqmax],phie[Neqmax],phiel[Neqmax],phier[Neqmax],UBCl[Neqmax];
        for(long k=0; k< g.nt; ++k){
            if(U[0][k]<1.e-13){
                for(int c=0; c<d; ++c) U[1+c][k] = 0.;
                U[0][k] = 0.;
            }
            if(g.surf){// surface: the momentum in the plane of the triangle
                double q[3]={U[1][k],U[2][k],U[3][k]};
                ProjPlane(k,q);
                for(int c=0; c<3; ++c) U[1+c][k] = q[c];
            }
        }
        vector< KN<double> > gU(Neq);
        if(ord>1) SWSlopes(U,&gU[0]);
        const KN<double> *pgU = ord>1 ? &gU[0] : 0;
        for(long k=0; k<g.nt; ++k){
            for(int e=0;e<g.nf;++e){
                long kf=k*g.nf+e;
                long ko=g.adj[kf];
                if(ko<0 || ko<k){// once per face ...
                    bool intern = ko >= 0;
                    const R3 &N=g.N[kf];
                    GetU(kf,U,pgU,Ul);
                    double areamin=g.meas[k];
                    if(intern){
                        long kfo=ko*g.nf+g.adjf[kf];
                        GetU(kfo,U,pgU,Ur);
                        if(g.surf){// surface: each momentum in the plane of its cell, the one of ko turned into the plane of k
                            ProjPlane(k,Ul+1);
                            ProjPlane(ko,Ur+1);
                            RotEdge(kf,kfo,Ur+1,1);
                        }
                        HydroRec(Ul,Ur[iz],Uls);
                        HydroRec(Ur,Ul[iz],Urs);
                        NFlux(kf,kfo,N,Uls,Urs,maxspeed,phie);// numerical flux
                        areamin=min(areamin,g.meas[ko]);
                        phiel[0] = phier[0] = phie[0];
                        for(int c=0; c<d; ++c) {
                            phiel[1+c] = phie[1+c] + .5*gr*(Ul[0]*Ul[0] - Uls[0]*Uls[0])*N[c];
                            phier[1+c] = phie[1+c] + .5*gr*(Ur[0]*Ur[0] - Urs[0]*Urs[0])*N[c];
                        }
                        phiel[iz] = phier[iz] = 0.;
                        if(ord>1){// the centered source term in the cell for the reconstruction
                            for(int c=0; c<d; ++c) {
                                phiel[1+c] += gr*(Ul[iz] - U[iz][k]) * (U[0][k] + Ul[0])/2.*N[c];
                                phier[1+c] += gr*(Ur[iz] - U[iz][ko]) * (U[0][ko] + Ur[0])/2.*N[c];
                            }
                        }
                        if(g.surf) RotEdge(kf,kfo,phier+1,-1);// surface: the momentum flux of ko in its plane
                    } else { // boundary face
                        if(g.surf) ProjPlane(k,Ul+1);
                        GetUBC(kf,uBC,UBCl);//here we consider that uBC[iz]=Ul[iz]
                        NFlux_BC(kf,g.flab[kf],N,Ul,maxspeed,phiel,UBCl);
                        phiel[iz] = 0.;
                        if(ord>1)// the same centered source term on a boundary face (lake at rest near the boundary)
                            for(int c=0; c<d; ++c)
                                phiel[1+c] += gr*(Ul[iz] - U[iz][k]) * (U[0][k] + Ul[0])/2.*N[c];
                    }
                    double le=g.fmeas[kf];
                    if(g.dim==2) {
                        bool okdt = (ord>1) ? (!cut[k] && !(intern && cut[ko])) : (intern || g.flab[kf]!=-111111);
                        if(okdt) dt=min(dt,areamin/le/(maxspeed+1.e-20));
                    }
                    else { sumf[k] += le*maxspeed; if(intern) sumf[ko] += le*maxspeed; }
                    // distribution of the fluxes to the cells
                    for(int j=0; j<Neq;++j){
                        phiel[j] *= le;
                        Up[j][k] += phiel[j]*(1./g.meas[k]);
                        if(intern) {
                            phier[j] *= le;
                            Up[j][ko] -= phier[j]*(1./g.meas[ko]);
                        }
                    }
                }
            }
        }
        if(g.dim!=2)
            for(long k=0; k<g.nt; ++k) if(!dtcut[k]) dt=min(dt,g.meas[k]/(sumf[k]+1.e-20));
        return dt;
    }//end of computedtFlux
    // numerical flux of an interior face: numflux = 0 HLL, 1 Rusanov, 2 Suliciu/HLLC
    void NFlux(long kfl,long kfr,const R3 &Nn,const double *Ul,const double *Ur,double &speed,double *phi) const {
        if(numflux != 2) {// 0 HLL, 1 Rusanov
            double Fr[Neqmax],Fl[Neqmax];
            if(numflux == 0){// HLL
                double sl = min(lambdamin(kfl,Nn,Ul),lambdamin(kfr,Nn,Ur));
                double sr = max(lambdamax(kfl,Nn,Ul),lambdamax(kfr,Nn,Ur));
                speed=max(abs(sl),abs(sr));
                Flux(kfl,Nn,Ul,Fl);
                Flux(kfr,Nn,Ur,Fr);
                for(int j=0; j<Neq;++j)
                    phi[j] = (sl>=0) ? Fl[j] : ((sr<=0) ? Fr[j] : (sr*Fl[j] - sl*Fr[j] + sl*sr*(Ur[j]-Ul[j]))/(sr-sl+1.e-10));
            }
            else if(numflux == 1){// Rusanov
                double cRus = max(absvp(kfl,Nn,Ul),absvp(kfr,Nn,Ur));
                speed=cRus;
                Flux(kfl,Nn,Ul,Fl);
                Flux(kfr,Nn,Ur,Fr);
                for(int j=0; j<Neq;++j)
                    phi[j] = 0.5*((Fr[j]+Fl[j]) - cRus*(Ur[j]-Ul[j]));
            }
            else ExecError("FreeVol swtopo: numflux not defined");
            return;
        }
        double Np[2],stateL[4],stateR[4],F[4];
        double hl = max(Ul[0],1.e-13), hr = max(Ur[0],1.e-13);
        stateL[0]=Ul[0]; stateL[3]=0.5*parameters[0]*Ul[0];
        stateR[0]=Ur[0]; stateR[3]=0.5*parameters[0]*Ur[0];
        if(d<=2){// the 2D flux (1D: the velocity on x, tangential velocity 0)
            Np[0]=Nn.x;Np[1]=d==2 ? Nn.y : 0.;
            stateL[1]=Ul[1]/hl; stateL[2]=d==2 ? Ul[2]/hl : 0.;
            stateR[1]=Ur[1]/hr; stateR[2]=d==2 ? Ur[2]/hr : 0.;
            numfluxgd(2.,stateL,stateR,Np,speed,F);
            for(int c=0; c<=d; ++c) phi[c]=F[c];
        } else {// 3D: in the frame (N,T1,T2), the 2D flux with each tangential velocity (passive)
            R3 T1 = fabs(Nn.x) < 0.9 ? R3(1.,0.,0.)^Nn : R3(0.,1.,0.)^Nn;
            T1 /= T1.norme();
            R3 T2 = Nn^T1, Tt[2]={T1,T2};
            R3 ul(Ul[1]/hl,Ul[2]/hl,Ul[3]/hl), ur(Ur[1]/hr,Ur[2]/hr,Ur[3]/hr), Fm=Nn*0.;
            Np[0]=1.;Np[1]=0.;
            stateL[1]=(ul,Nn); stateR[1]=(ur,Nn);
            for(int t=0; t<2; ++t){
                stateL[2]=(ul,Tt[t]); stateR[2]=(ur,Tt[t]);
                numfluxgd(2.,stateL,stateR,Np,speed,F);
                Fm += Tt[t]*F[2];
            }
            Fm += Nn*F[1];
            phi[0]=F[0];
            for(int c=0; c<3; ++c) phi[1+c]=Fm[c];
        }
        phi[iz]=0.;
        if (speed>1.e7) {
            cout << "blowup of speed " << speed << endl;
            cout << "left state " ; for(int j=0;j<Neq;++j) cout << Ul[j] << " "; cout << endl;
            cout << "right state "; for(int j=0;j<Neq;++j) cout << Ur[j] << " "; cout << endl;
            ExecError("FreeVol swtopo: blowup of speed");
        }
    }
    void NFlux_BC(long kf,int labelement,const R3 &Nn,const double *UL,double &speed,double *phi, const double *UBCL) const {
        int choice_BC = 1;// 0 Dirichlet, 1 Neumann(default), 2 Wall, 3 mix Dirichlet and wall
        for(int i=0;i<labelBC.n/2;i++){
            if(labelement==labelBC[2*i]){
                choice_BC = labelBC[2*i+1];
            }
        }
        // Inner state UL = [h, hu, (hv, hw,) zb]
        double h  = UL[0];
        double un = 0.;// normal velocity
        for(int c=0; c<d; ++c) un += UL[1+c]/max(h,1.e-13)*Nn[c];
        // Outer state (ghost)
        double UR[Neqmax];
        if(choice_BC == 0){
            for (int j = 0; j < Neq; ++j) UR[j]=UBCL[j];
        } else if(choice_BC == 1){
            for (int j = 0; j < Neq; ++j) UR[j]=UL[j];
        } else if(choice_BC == 2 || choice_BC == 3){
            // build the "mirror" state (wall: un reversed, tangential velocity unchanged)
            UR[0] = choice_BC == 2 ? h : UBCL[0];
            for(int c=0; c<d; ++c) {
                double hvR = h*(UL[1+c]/max(h,1.e-13) - 2.*un*Nn[c]);
                UR[1+c] = (choice_BC == 2 || UBCL[1+c]==0.) ? hvR : UBCL[1+c];
            }
            UR[iz] = UL[iz];
        } else {cout<<"choice_BC = "<<choice_BC<<" is not defined for swtopo"<<endl;ExecError("FreeVol swtopo: labelBC not defined");}
        NFlux(kf,kf,Nn,UL,UR,speed,phi);
    }
    ~FreeVol_SWTOPO() { if(verbosity>99) cout << " delete FreeVol_SWTOPO " << this << endl; }
};

// creation of the model for fvscheme FV([H,Hu,Hv,zb], eq="swtopo", ...) (see the registry in FreeVol.cpp)
static FreeVol_model *FreeVol_SWTOPO_Create(const vector< FreeVol_FE > &V, long order, long numflux,
                                            KN<long> labelBC, long limiter, KN<double> parameters,
                                            KN<long> intparameters,
                                            const void *dual) {
    if (dual) ExecError("fvscheme eq=\"swtopo\": this model does not take a dual mesh (fvscheme FV([...], eq=\"swtopo\", ...))");
    if (numflux == -99) numflux = 2;    // Suliciu/HLLC by default
    // parameters[0] is the gravity constant
    if (parameters.N() < 1 || parameters[0] <= 0.) {
        cout << "gravity parameter is not valid, it must be given in parameters=" << endl;
        ExecError("FreeVol swtopo: parameters=[gravity] expected");
    }
    if ((verbosity) && (!mpirank)) {
        cout << "gravity : " << parameters[0] << endl;
        if ((intparameters.N() > 1) || ((intparameters.N() == 1) && intparameters[0] != 0)) {    // intparameters is not yet used
            cout << "intparameters = [ ";
            for (int i = 0; i < intparameters.N() - 1; i++) cout << intparameters[i] << " , ";
            cout << intparameters[intparameters.N() - 1] << " ]" << endl;
        }
        for (int i = 0; i < labelBC.n / 2; i++) {
            if (labelBC[2 * i + 1] == 0) cout << "on boundary " << labelBC[2 * i] << " labelBC=0 : Dirichlet BC" << endl;
            else if (labelBC[2 * i + 1] == 1) cout << "on boundary " << labelBC[2 * i] << " labelBC=1 : Neumann BC" << endl;
            else if (labelBC[2 * i + 1] == 2) cout << "on boundary " << labelBC[2 * i] << " labelBC=2 : Wall BC" << endl;
            else if (labelBC[2 * i + 1] == 3) cout << "on boundary " << labelBC[2 * i] << " labelBC=3 : Dirichlet Wall BC" << endl;
            else cout << "on boundary " << labelBC[2 * i] << " labelBC=" << labelBC[2 * i + 1] << " is not defined" << endl;
        }
        if (numflux == 0) cout << "numflux = 0 : HLL" << endl;
        else if (numflux == 1) cout << "numflux = 1 : Rusanov" << endl;
        else if (numflux == 2) cout << "numflux = 2 : Suliciu/HLLC" << endl;
        if (limiter == 0) cout << "limiter = 0 : no limiter is used" << endl;
        else if (limiter == 1) cout << "limiter = 1 : maximum principle type" << endl;
        else if (limiter == 2) cout << "limiter = 2 : minmod" << endl;
        else if (limiter == 3) cout << "limiter = 3 : Sweby" << endl;
        else if (limiter == 4) cout << "limiter = 4 : Van Albada type" << endl;
        else cout << "limiter = " << limiter << endl;
    }
    // a numerical flux not defined for swtopo: stop (whatever the verbosity)
    if (numflux < 0 || numflux > 2) {
        cout << "numflux = " << numflux << " is not defined for swtopo (0 HLL, 1 Rusanov, 2 Suliciu/HLLC)" << endl;
        ExecError("FreeVol swtopo: numflux not defined");
    }
    return new FreeVol_SWTOPO(V, order, numflux, labelBC, limiter, parameters, intparameters);
}

static void Load_Init() {
    if (verbosity && !mpirank) cout << " load: FreeVol_swtopo" << endl;
    FreeVol_RegisterModel("swtopo", FreeVol_SWTOPO_Create);
}
LOADFUNC(Load_Init)
