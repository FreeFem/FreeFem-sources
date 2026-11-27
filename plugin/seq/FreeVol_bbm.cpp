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
// ff-c++ -shared -fPIC FreeVol_bbm.cpp -I$FREEFEMPP_DIR/include -L$FREEFEMPP_DIR/lib FreeVol.so
// or simply ff-c++ -shared -fPIC FreeVol_bbm.cpp FreeVol.so
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...
/* ---------------------------------------------------------

   the hyperbolic part of the BBM-BBM (Boussinesq) system, for a finite element / finite volume method:
     A_BBM dE/dt + div([F(E),G(E)]) = 0,  A_BBM = diag(I - b Delta, I - d Delta, I - d Delta),
     E = [eta,u,v], F(E) = ((1+eta) u, eta + |u|^2/2, 0), G(E) = ((1+eta) v, 0, eta + |u|^2/2)
   this model gives R = div([F(E),G(E)]) by the finite volume method (fluxupdate): the script solves then
   A_BBM W = -R by finite elements and updates E (see freevol_bbm1d.edp, freevol_bbm2d.edp)
   in 1D, 2D or 3D (d components of the velocity u), the dimension is given by the mesh of the unknowns:
     1D: meshL, fvscheme FV([eta,u], eq="bbm", ...)    2D: mesh, fvscheme FV([eta,u,v], eq="bbm", ...)
     3D: mesh3, fvscheme FV([eta,u,v,w], eq="bbm", ...)
   the eigenvalues of the jacobian of F(E).N: u.N - c, 0 (d-1 times), u.N + c, with c = sqrt(1+eta), 1+eta > 0
   second-order in space by slope reconstruction (the generic scheme of FreeVol_model)
   numflux = 0 : HLL (default), numflux = 1 : Rusanov
   labelBC = 0 : Dirichlet, 1 : Neumann, 2 : Wall (u.N = 0)
   unknown E is P0 (one value by cell)

   --------------------------------------------------------*/

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include "FreeVol.hpp"

// the hyperbolic part of the BBM-BBM system, in 1D, 2D or 3D
// the unknowns: U[0] = eta, U[1..d] = u (d components)
class FreeVol_BBM: public FreeVol_model
{
public:
    int d;     // dimension, d components of the velocity

    FreeVol_BBM(const vector< FreeVol_FE > &Vfe,int oord,int snumflux,KN<long> labelBCC,int limiterr,KN<double> parameterss,KN<long> intparameterss)
    : FreeVol_model(), d(0)
    { init(Vfe,oord,snumflux,labelBCC,limiterr,parameterss,intparameterss);}

    void init(const vector< FreeVol_FE > &Vfe,int oord,int snumflux,KN<long> labelBCC,int limiterr,KN<double> parameterss,KN<long> intparameterss){
        int nU = Vfe.size();
        d = nU ? Vfe[0].dim : 0;
        if (d < 1 || nU != d+1)
            ExecError("fvscheme bbm: the unknowns are [eta,u] in 1D (meshL), [eta,u,v] in 2D (mesh), [eta,u,v,w] in 3D (mesh3)");
        for (int c = 0; c < nU; ++c)
            if (!Vfe[c].x()) { cout << " Sorry unknown " << c+1 << " is unset "<< endl; ffassert(0);}
        const FreeVol_FE &f = Vfe[0];
        if (f.MaxNbDFPerElement() / f.N() != 1) ExecError("fvscheme bbm: the unknowns must be P0");
        int neq = d+1;//defines the number of equations (unknowns) of the system
        set(neq,oord,snumflux,f,labelBCC,"bbm",limiterr,parameterss,intparameterss);
        if (f.NbOfNodes() != geom.nt) ExecError("fvscheme bbm: the unknowns must be P0");
    }//end of init

    // Flux F(E).N = ( (1+eta) u.N , (eta + |u|^2/2) N )
    void Flux(long kf,const R3 &N,const double *U,double *F) const {
        double un = 0., u2 = 0.;
        for (int c = 0; c < d; ++c) { un += U[1+c]*N[c]; u2 += U[1+c]*U[1+c]; }
        F[0] = (1.+U[0])*un;
        for (int c = 0; c < d; ++c) F[1+c] = (U[0] + 0.5*u2)*N[c];
    }
    // the eigenvalues: u.N - c, 0 (d-1 times), u.N + c, c = sqrt(1+eta)
    void lambda(long kf,const R3 &N,const double *U,double *l) const {
        double un = 0.;
        for (int c = 0; c < d; ++c) un += U[1+c]*N[c];
        if (1.+U[0] <= 0.) ExecError("FreeVol bbm: 1+eta <= 0 (the system is not hyperbolic)");
        double cs = sqrt(1.+U[0]);
        for (int j = 0; j < Neq; ++j) l[j] = 0.;
        l[0] = un - cs;
        l[d] = un + cs;
    }
    // the fluxes and the time step (CFL = 1), the speed of each face given by its numerical flux:
    // 2D: min on the edges of area/(length*speed) (enough in 2D: for a velocity constant in a triangle,
    // the greatest flux of an edge is the sum of the outgoing fluxes);
    // 1D, 3D: min on the cells of measure/sum on the faces of (measure of the face*speed): a
    // tetrahedron can have 2 outgoing faces, each one with half of the outgoing flux
    double computedtFlux(KN_d *U,KN_d *Up, KN_d *uBC) const {
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
    // numerical flux of an interior face: numflux = 0 HLL, 1 Rusanov
    void NFlux(long kfl,long kfr,const R3 &Nn,const double *Ul,const double *Ur,double &speed,double *phi) const {
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
        else ExecError("FreeVol bbm: numflux not defined");
    }
    void NFlux_BC(long kf,int labelement,const R3 &Nn,const double *UL,double &speed,double *phi, const double *UBCL) const {
        int choice_BC = 1;// 0 Dirichlet, 1 Neumann(default), 2 Wall
        for(int i=0;i<labelBC.n/2;i++)
            if(labelement==labelBC[2*i]) choice_BC = labelBC[2*i+1];
        double un = 0.;
        for(int c=0; c<d; ++c) un += UL[1+c]*Nn[c];
        double UR[Neqmax];// the state outside
        if(choice_BC == 0){
            for (int j = 0; j < Neq; ++j) UR[j]=UBCL[j];
        } else if(choice_BC == 1){
            for (int j = 0; j < Neq; ++j) UR[j]=UL[j];
        } else if(choice_BC == 2){// wall: the normal velocity reversed
            UR[0] = UL[0];
            for(int c=0; c<d; ++c) UR[1+c] = UL[1+c] - 2.*un*Nn[c];
        } else {cout<<"choice_BC = "<<choice_BC<<" is not defined for bbm"<<endl;ExecError("FreeVol bbm: labelBC not defined");}
        NFlux(kf,kf,Nn,UL,UR,speed,phi);
    }
    ~FreeVol_BBM() {if(verbosity>99)  cout << " delete FreeVol_BBM "<< this << endl;}
};

// creation of the model for fvscheme FV([eta,u,v], eq="bbm", ...) (see the registry in FreeVol.cpp)
static FreeVol_model *FreeVol_BBM_Create(const vector< FreeVol_FE > &V, long order, long numflux,
                                         KN<long> labelBC, long limiter, KN<double> parameters,
                                         KN<long> intparameters,
                                         const void *dual) {
    if (dual) ExecError("fvscheme eq=\"bbm\": this model does not take a dual mesh (fvscheme FV([...], eq=\"bbm\", ...))");
    if (numflux == -99) numflux = 0;    // HLL by default
    if ((verbosity) && (!mpirank)) {
        for (int i = 0; i < labelBC.n / 2; i++) {
            if (labelBC[2 * i + 1] == 0) cout << "on boundary " << labelBC[2 * i] << " labelBC=0 : Dirichlet BC" << endl;
            else if (labelBC[2 * i + 1] == 1) cout << "on boundary " << labelBC[2 * i] << " labelBC=1 : Neumann BC" << endl;
            else if (labelBC[2 * i + 1] == 2) cout << "on boundary " << labelBC[2 * i] << " labelBC=2 : Wall BC" << endl;
            else cout << "on boundary " << labelBC[2 * i] << " labelBC=" << labelBC[2 * i + 1] << " is not defined" << endl;
        }
        if (numflux == 0) cout << "numflux = 0 : HLL" << endl;
        else if (numflux == 1) cout << "numflux = 1 : Rusanov" << endl;
        if (limiter == 0) cout << "limiter = 0 : no limiter is used" << endl;
        else if (limiter == 1) cout << "limiter = 1 : maximum principle type" << endl;
        else cout << "limiter = " << limiter << endl;
    }
    if (numflux < 0 || numflux > 1) {
        cout << "numflux = " << numflux << " is not defined for bbm (0 HLL, 1 Rusanov)" << endl;
        ExecError("FreeVol bbm: numflux not defined");
    }
    return new FreeVol_BBM(V, order, numflux, labelBC, limiter, parameters, intparameters);
}

static void Load_Init() {
    if (verbosity && !mpirank) cout << " load: FreeVol_bbm" << endl;
    FreeVol_RegisterModel("bbm", FreeVol_BBM_Create);
}
LOADFUNC(Load_Init)
