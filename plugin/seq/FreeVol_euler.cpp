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
// ff-c++ -shared -fPIC FreeVol_euler.cpp -I$FREEFEMPP_DIR/include -L$FREEFEMPP_DIR/lib FreeVol.so
// or simply ff-c++ -shared -fPIC FreeVol_euler.cpp FreeVol.so
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...
/* ---------------------------------------------------------

   finite volume method to solve the compressible Euler system (full gas dynamics with internal energy)
   U_t + div(F(U)) = 0, U=[rho,rho*u,E], E = rho*|u|^2/2+rho*e, rho>=0, e>=0, u the velocity (d components)
   F(U) = ( rho*u , rho*u x u + p Id , (E+p) u ), p=(gamma-1)*rho*e
   in 1D, 2D or 3D, the dimension d is given by the mesh of the unknowns:
     1D: meshL (segments on the x axis), fvscheme FV([rho,rhou,E], eq="euler", ...)
     2D: mesh (triangles),               fvscheme FV([rho,rhou,rhov,E], eq="euler", ...)
     3D: mesh3 (tetrahedra),             fvscheme FV([rho,rhou,rhov,rhow,E], eq="euler", ...)
   gamma=parameters[0]  --> parameters=[gamma]
   second-order in space by slope reconstruction (on rho, u and rho*e, with limitations for rho and rho*e
   positive and for the velocity: positivity of the internal energy), the generic scheme of FreeVol_model for
   the fluxes and the time step
   numflux = 0 : HLL, numflux = 1 : Rusanov, numflux = 2 : Suliciu/HLLC (3D: numfluxgd in the frame (N,T1,T2)
   with the tangential velocity T1, the velocity T2 passive: upwind by the mass flux)
   labelBC = 0 : Dirichlet, 1 : Neumann, 2 : Wall, 3 : Dirichlet and wall (a velocity component 0 in the
   Dirichlet data: the one of the wall)
   unknown U is P0 (one value by cell)

   --------------------------------------------------------*/

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include "FreeVol.hpp"

// system of compressible Euler equations, in 1D, 2D or 3D
// the unknowns: U[0] = rho, U[1..d] = rho*u (d components), U[d+1] = E
class FreeVol_Euler: public FreeVol_model
{
public:
    int d;     // dimension, d components of the velocity
    int iE;    // index of the energy: d+1

    FreeVol_Euler(const vector< FreeVol_FE > &Vfe,int oord,int snumflux,KN<long> labelBCC,int limiterr,KN<double> parameterss,KN<long> intparameterss)
    : FreeVol_model(), d(0), iE(0)
    { init(Vfe,oord,snumflux,labelBCC,limiterr,parameterss,intparameterss);}

    void init(const vector< FreeVol_FE > &Vfe,int oord,int snumflux,KN<long> labelBCC,int limiterr,KN<double> parameterss,KN<long> intparameterss){
        int nU = Vfe.size();
        d = nU ? Vfe[0].dim : 0;
        if (d < 1 || nU != d+2)
            ExecError("fvscheme euler: the unknowns are [rho,rhou,E] in 1D (meshL), [rho,rhou,rhov,E] in 2D (mesh), [rho,rhou,rhov,rhow,E] in 3D (mesh3)");
        for (int c = 0; c < nU; ++c)
            if (!Vfe[c].x()) { cout << " Sorry unknown " << c+1 << " is unset "<< endl; ffassert(0);}
        const FreeVol_FE &f = Vfe[0];
        if (f.MaxNbDFPerElement() / f.N() != 1) ExecError("fvscheme euler: the unknowns must be P0");
        iE = d+1;
        int neq = d+2;//defines the number of equations (unknowns) of the system
        set(neq,oord,snumflux,f,labelBCC,"euler",limiterr,parameterss,intparameterss);
        if (f.NbOfNodes() != geom.nt) ExecError("fvscheme euler: the unknowns must be P0");
    }//end of init

    // the kinetic energy rho*|u|^2/2
    double kinetic(const double *U) const {
        double q2 = 0.;
        for (int c = 0; c < d; ++c) q2 += U[1+c]*U[1+c];
        return .5*q2/max(U[0],1.e-13);
    }
    // E >= rho*|u|^2/2 + 1e-12 rho (a positive internal energy), the state V modified (a copy of U)
    void FixEnergy(const double *U, double *V) const {
        for (int j = 0; j < Neq; ++j) V[j] = U[j];
        V[iE] = max(V[iE], kinetic(V)+1.e-12*max(V[0],1.e-13));
    }
    // Flux F(U).N
    void Flux(long kf,const R3 &N,const double *U,double *F) const {
        double h = U[0];
        for (int j = 0; j < Neq; ++j) F[j] = 0.;
        if (h <= 1.e-13) return;
        double un = 0.;
        for (int c = 0; c < d; ++c) un += U[1+c]/h*N[c];
        double e = max((U[iE]-kinetic(U))/h,1.e-12);
        double p = (parameters[0]-1.)*h*e;
        F[0] = h*un;
        for (int c = 0; c < d; ++c) F[1+c] = U[1+c]*un + p*N[c];
        F[iE] = (U[iE]+p)*un;
    }
    // Local eigenvalues: un-c, un (d times), un+c
    void lambda(long kf,const R3 &N,const double *U,double *l) const {
        double h = U[0];
        for (int j = 0; j < Neq; ++j) l[j] = 0.;
        if (h <= 1.e-13) return;
        double un = 0.;
        for (int c = 0; c < d; ++c) un += U[1+c]/h*N[c];
        double e = max((U[iE]-kinetic(U))/h,1.e-12);
        double cs = sqrt(parameters[0]*(parameters[0]-1.)*e);
        l[0] = un - cs;
        for (int c = 1; c <= d; ++c) l[c] = un;
        l[iE] = un + cs;
    }
    // the limited slopes: on rho, u and rho*e, then limitations for rho and rho*e positive and for the velocity
    // (positivity of the internal energy); used by computedtFlux
    void LimitedSlopes(KN_d *U, KN<double> *gU) const {
        const FreeVol_Geom &g = geom;
        long nt = g.nt;
        KN<double> h(nt),hei(nt),vel[3],alpha(nt);
        KN<double> grad_h(d*nt),grad_hei(d*nt),grad_vel[3];
        for(int c=0; c<d; ++c) { vel[c].resize(nt); grad_vel[c].resize(d*nt); }
        for(long k=0; k< nt; ++k){
            double Uk[Neqmax];
            for(int j=0; j<Neq; ++j) Uk[j] = U[j][k];
            h[k] = U[0][k];
            for(int c=0; c<d; ++c) vel[c][k] = U[1+c][k]/max(U[0][k],1.e-13);
            hei[k] = U[iE][k]-kinetic(Uk);
        }
        MatriceMorse<R> &MATGRAD = *(matgrad);
        grad_h = MATGRAD * h;
        grad_hei = MATGRAD * hei;
        for(int c=0; c<d; ++c) grad_vel[c] = MATGRAD * vel[c];
        BdrySlope(grad_h,1L);
        BdrySlope(grad_hei,1L);
        for(int c=0; c<d; ++c) BdrySlope(grad_vel[c],1L);
        SlopeLimiterVF(h,grad_h,alpha,limiter);
        for(long k=0; k< d*nt; ++k) grad_h[k] *= alpha[k/d];
        SlopeLimiterVF(hei,grad_hei,alpha,limiter);
        for(long k=0; k< d*nt; ++k) grad_hei[k] *= alpha[k/d];
        for(int c=0; c<d; ++c) {
            SlopeLimiterVF(vel[c],grad_vel[c],alpha,limiter);
            for(long k=0; k< d*nt; ++k) grad_vel[c][k] *= alpha[k/d];
        }
        for(long k=0; k< nt; ++k){
            if (U[0][k]<1.e-13) {// vacuum
                for(int i=0; i<d; ++i) {
                    grad_h[d*k+i] = 0.;
                    grad_hei[d*k+i] = 0.;
                    for(int c=0; c<d; ++c) grad_vel[c][d*k+i] = 0.;
                }
            }
        }
        for(int j=0; j<Neq; ++j) gU[j].resize(d*nt);
        for(long k=0; k< nt; ++k){
            double hk = U[0][k];
            // limiter for rho positive
            double alphah = 1.;
            for (int e=0;e<g.nf;++e) {
                R3 GE = g.F[k*g.nf+e]-g.G[k];// centroid of k to the centroid of face e
                double he = hk;
                for(int i=0; i<d; ++i) he += grad_h[d*k+i]*GE[i];
                alphah = min(alphah,1./max((1.-he/max(hk,1.e-13))/0.9,1.));
            }
            for(int i=0; i<d; ++i) grad_h[d*k+i] *= alphah;
            // limiter for rho*e positive
            double alphahei = 1.;
            for (int e=0;e<g.nf;++e) {
                R3 GE = g.F[k*g.nf+e]-g.G[k];
                double heie = hei[k];
                for(int i=0; i<d; ++i) heie += grad_hei[d*k+i]*GE[i];
                alphahei = min(alphahei,1./max((1.-heie/max(hei[k],1.e-30))/0.9,1.));
            }
            for(int i=0; i<d; ++i) grad_hei[d*k+i] *= alphahei;
            // limitation on the velocity for the positivity of the internal energy
            double maxvellim = 0.;
            for (int e=0;e<g.nf;++e) {
                R3 GE = g.F[k*g.nf+e]-g.G[k];
                double he = hk, heie = hei[k];
                for(int i=0; i<d; ++i) { he += grad_h[d*k+i]*GE[i]; heie += grad_hei[d*k+i]*GE[i]; }
                double inc2 = 0.;
                for(int c=0; c<d; ++c) {
                    double incue = 0.;
                    for(int i=0; i<d; ++i) incue += grad_vel[c][d*k+i]*GE[i];
                    incue *= hk;
                    incue /= max(he,1.e-13);
                    inc2 += incue*incue;
                }
                maxvellim = max(maxvellim,sqrt(.5*inc2/max(heie/max(he,1.e-13),1.e-30)));
            }
            double alphavel = 1./max(maxvellim/.5,1.);
            //define the final slopes
            double v2 = 0.;
            for(int c=0; c<d; ++c) v2 += vel[c][k]*vel[c][k];
            for(int i=0; i<d; ++i) {
                long ki = d*k+i;
                for(int c=0; c<d; ++c) grad_vel[c][ki] *= alphavel;
                gU[0][ki] = grad_h[ki];
                double gE = .5*v2*grad_h[ki] + grad_hei[ki];
                for(int c=0; c<d; ++c) {
                    gU[1+c][ki] = grad_h[ki]*vel[c][k]+grad_vel[c][ki]*hk;
                    gE += hk*vel[c][k]*grad_vel[c][ki];
                }
                gU[iE][ki] = gE;
            }
        }
    }
    // vacuum: rho < 1e-13 -> U = 0; else E >= rho*|u|^2/2 + 1e-12 rho; then the generic scheme
    void Vacuum(KN_d *U) const {
        for(long k=0; k< geom.nt; ++k){
            if(U[0][k]<1.e-13) { for(int j=0; j<Neq; ++j) U[j][k] = 0.; }
            else {
                double Uk[Neqmax];
                for(int j=0; j<Neq; ++j) Uk[j] = U[j][k];
                U[iE][k] = max(U[iE][k], kinetic(Uk)+1.e-12*U[0][k]);
            }
        }
    }
    // the fluxes and the time step (CFL = 1), the speed of each face given by its numerical flux:
    // 2D: min on the edges of area/(length*speed) (enough in 2D: for a velocity constant in a triangle,
    // the greatest flux of an edge is the sum of the outgoing fluxes);
    // 1D, 3D: min on the cells of measure/sum on the faces of (measure of the face*speed): a
    // tetrahedron can have 2 outgoing faces, each one with half of the outgoing flux
    double computedtFlux(KN_d *U,KN_d *Up, KN_d *uBC) const {
        Vacuum(U);
        const FreeVol_Geom &g = geom;
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
        for(long j=0; j<Neq; ++j) Up[j]=0.;
        //left and right states and fluxes for an interface
        double Ul[Neqmax],Ur[Neqmax],phie[Neqmax],UBCl[Neqmax];
        vector< KN<double> > gU(Neq);
        if(ord>1) LimitedSlopes(U,&gU[0]);
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
                        NFlux(kf,kfo,N,Ul,Ur,maxspeed,phie);// numerical flux
                        areamin = min(areamin,g.meas[ko]);
                    } else { // boundary face
                        GetUBC(kf,uBC,UBCl);
                        NFlux_BC(kf,g.flab[kf],N,Ul,maxspeed,phie,UBCl);
                    }
                    double le=g.fmeas[kf];
                    if(g.dim==2) {
                        bool okdt = (ord>1) ? (!cut[k] && !(intern && cut[ko])) : (intern || g.flab[kf]!=-111111);
                        if(okdt) dt=min(dt,areamin/le/(maxspeed+1.e-20));
                    }
                    else { sumf[k] += le*maxspeed; if(intern) sumf[ko] += le*maxspeed; }
                    // distribution of the fluxes to the cells
                    for(long j=0; j<Neq; ++j){
                        phie[j] *= le;
                        Up[j][k] += phie[j]*(1./g.meas[k]);
                        if(intern) Up[j][ko] -= phie[j]*(1./g.meas[ko]);
                    }
                }
            }
        }
        if(g.dim!=2)
            for(long k=0; k<g.nt; ++k) if(!dtcut[k]) dt=min(dt,g.meas[k]/(sumf[k]+1.e-20));
        return dt;
    }
    // numerical flux of an interior face, the reconstructed states with E >= rho*|u|^2/2 + 1e-12 rho:
    // numflux = 0 HLL, 1 Rusanov, 2 Suliciu/HLLC
    void NFlux(long kfl,long kfr,const R3 &Nn,const double *Ul0,const double *Ur0,double &speed,double *phi) const {
        double Ul[Neqmax], Ur[Neqmax];
        FixEnergy(Ul0,Ul);
        FixEnergy(Ur0,Ur);
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
            else ExecError("FreeVol euler: numflux not defined");
            return;
        }
        double gamma = parameters[0];
        // numfluxgd: the 2D flux (mass, 2 momentum, energy) with the states [rho,u,v,e]
        double Np[2],stateL[4],stateR[4],F[4];
        double hl = max(Ul[0],1.e-13), hr = max(Ur[0],1.e-13);
        stateL[0]=Ul[0]; stateL[3]=(Ul[iE]-kinetic(Ul))/hl;// the internal energy (all the velocity components)
        stateR[0]=Ur[0]; stateR[3]=(Ur[iE]-kinetic(Ur))/hr;
        if(d<=2){// 1D: the velocity on x, tangential velocity 0
            Np[0]=Nn.x; Np[1]=d==2 ? Nn.y : 0.;
            stateL[1]=Ul[1]/hl; stateL[2]=d==2 ? Ul[2]/hl : 0.;
            stateR[1]=Ur[1]/hr; stateR[2]=d==2 ? Ur[2]/hr : 0.;
            numfluxgd(gamma,stateL,stateR,Np,speed,F);
            phi[0]=F[0];
            for(int c=0; c<d; ++c) phi[1+c]=F[1+c];
            phi[iE]=F[3];
        } else {// 3D: in the frame (N,T1,T2), the 2D flux with the tangential velocity T1; in numfluxgd the
                // tangential velocity is passive (flux of mass * the upwind value, the side given by the sign of
                // the mass flux): the velocity T2 the same, its momentum flux F0*t2 and its kinetic energy F0*t2^2/2
            R3 T1 = fabs(Nn.x) < 0.9 ? R3(1.,0.,0.)^Nn : R3(0.,1.,0.)^Nn;
            T1 /= T1.norme();
            R3 T2 = Nn^T1;
            R3 ul(Ul[1]/hl,Ul[2]/hl,Ul[3]/hl), ur(Ur[1]/hr,Ur[2]/hr,Ur[3]/hr);
            Np[0]=1.; Np[1]=0.;
            stateL[1]=(ul,Nn); stateL[2]=(ul,T1);
            stateR[1]=(ur,Nn); stateR[2]=(ur,T1);
            numfluxgd(gamma,stateL,stateR,Np,speed,F);
            double t2 = F[0] >= 0. ? (ul,T2) : (ur,T2);// upwind
            R3 Fm = Nn*F[1] + T1*F[2] + T2*(F[0]*t2);
            phi[0]=F[0];
            for(int c=0; c<3; ++c) phi[1+c]=Fm[c];
            phi[iE]=F[3] + F[0]*.5*t2*t2;
        }
        if (speed>1.e7) {
            cout << "blowup of speed " << speed << endl;
            cout << "left state " ; for(int j=0;j<Neq;++j) cout << Ul[j] << " "; cout << endl;
            cout << "right state "; for(int j=0;j<Neq;++j) cout << Ur[j] << " "; cout << endl;
            ExecError("FreeVol euler: blowup of speed");
        }
    }
    void NFlux_BC(long kf,int labelement,const R3 &Nn,const double *UL0,double &speed,double *phi, const double *UBCL) const {
        int choice_BC = 1;// 0 Dirichlet, 1 Neumann(default), 2 Wall, 3 mix Dirichlet and wall
        for(int i=0;i<labelBC.n/2;i++){
            if(labelement==labelBC[2*i]){
                choice_BC = labelBC[2*i+1];
            }
        }
        double UL[Neqmax];
        FixEnergy(UL0,UL);
        double h  = UL[0];
        double un = 0.;// normal velocity
        for(int c=0; c<d; ++c) un += UL[1+c]/max(h,1.e-13)*Nn[c];
        // Outer state (ghost)
        double UR[Neqmax];
        if(choice_BC == 0){
            for (int j = 0; j < Neq; ++j) UR[j]=UBCL[j];
        } else if(choice_BC == 1){
            for (int j = 0; j < Neq; ++j) UR[j]=UL[j];
        } else if(choice_BC == 2){
            // build the "mirror" state (wall: un reversed, tangential velocity unchanged)
            UR[0] = h;
            for(int c=0; c<d; ++c) UR[1+c] = h*(UL[1+c]/max(h,1.e-13) - 2.*un*Nn[c]);
            UR[iE] = UL[iE];
        } else if(choice_BC == 3){
            // Dirichlet, a velocity component 0 in the data: the one of the wall
            UR[0] = UBCL[0];
            for(int c=0; c<d; ++c) UR[1+c] = (UBCL[1+c]==0.) ? h*(UL[1+c]/max(h,1.e-13) - 2.*un*Nn[c]) : UBCL[1+c];
            UR[iE] = UBCL[iE];
        } else {cout<<"choice_BC = "<<choice_BC<<" is not defined for euler"<<endl;ExecError("FreeVol euler: labelBC not defined");}
        NFlux(kf,kf,Nn,UL,UR,speed,phi);
    }
    ~FreeVol_Euler() {if(verbosity>99)  cout << " delete FreeVol_Euler "<< this << endl;}
};

// creation of the model for fvscheme FV([rho,rhou,rhov,E], eq="euler", ...) (see the registry in FreeVol.cpp)
static FreeVol_model *FreeVol_Euler_Create(const vector< FreeVol_FE > &V, long order, long numflux,
                                           KN<long> labelBC, long limiter, KN<double> parameters,
                                           KN<long> intparameters,
                                           const void *dual) {
    if (dual) ExecError("fvscheme eq=\"euler\": this model does not take a dual mesh (fvscheme FV([...], eq=\"euler\", ...))");
    if (numflux == -99) numflux = 2;    // Suliciu/HLLC by default
    // parameters[0] is the exponent gamma of the gas law p=(gamma-1)*rho*e
    if (parameters.N() < 1 || parameters[0] <= 1.) {
        cout << "parameter gamma > 1 must be given in parameters=[gamma]" << endl;
        ExecError("FreeVol euler: parameters=[gamma] expected");
    }
    if ((verbosity) && (!mpirank)) {
        cout << "gamma : " << parameters[0] << endl;
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
    // a numerical flux not defined for euler: stop (whatever the verbosity)
    if (numflux < 0 || numflux > 2) {
        cout << "numflux = " << numflux << " is not defined for euler (0 HLL, 1 Rusanov, 2 Suliciu/HLLC)" << endl;
        ExecError("FreeVol euler: numflux not defined");
    }
    return new FreeVol_Euler(V, order, numflux, labelBC, limiter, parameters, intparameters);
}

static void Load_Init() {
    if (verbosity && !mpirank) cout << " load: FreeVol_euler" << endl;
    FreeVol_RegisterModel("euler", FreeVol_Euler_Create);
}
LOADFUNC(Load_Init)
