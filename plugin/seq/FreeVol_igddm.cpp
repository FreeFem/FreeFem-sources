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
// ff-c++ -shared -fPIC FreeVol_igddm.cpp -I$FREEFEMPP_DIR/include -L$FREEFEMPP_DIR/lib FreeVol.so
// or simply ff-c++ -shared -fPIC FreeVol_igddm.cpp FreeVol.so
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...
/* ---------------------------------------------------------

   finite volume method to solve the isentropic gas dynamics system
   U_t + div([F(U),G(U)]) = 0, U=[rho,rho*u,rho*v], rho>=0,
          (          rho*u            )         (          rho*v            )
   F(U) = ( rho*u^2 + kappa*rho^gamma ), G(U) = (         rho*u*v           )
          (         rho*u*v           )         ( rho*v^2 + kappa*rho^gamma )
   gamma=parameters[0], kappa=parameters[1]  --> parameters=[gamma,kappa]
   second-order in space by slope reconstruction
   numflux = 0 : HLL, numflux = 1 : Rusanov, numflux = 2 : HLLC
   unknown U is P1 (one value by vertex), the control volumes are the cells of the dual mesh of the mesh of the
   unknowns, built in the script and given to fvscheme; the method is given by the type of the dual mesh:
   1D (meshL):  meshdualL Thstar = builddualmesh(Th), the median dual mesh, fvscheme FV(Thstar, [rho,rhou], eq="igddm", ...)
   2D (mesh):   meshdual Thstar = builddualmesh(Th), the dual mesh of F. Bouchut, FV(Thstar, [rho,rhou,rhov], ...)
   3D (mesh3):  meshdual3 Thstar = builddualmesh(Th), the median dual mesh, FV(Thstar, [rho,rhou,rhov,rhow], ...)
   surface (meshS): meshdualS Thstar = builddualmesh(Th), the median dual mesh, FV(Thstar, [rho,rhou,rhov,rhow], ...), the momentum tangent
     (to the control volume: the state of the neighbour turned into its plane on each interface; a state at
     rest stays at rest: SurfDefect)
   labelBC = 0 : Dirichlet, 1 : Neumann, 2 : Wall

   --------------------------------------------------------*/

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include "FreeVolDualMesh.hpp"  // FreeVol.hpp, the dual meshes (fvdualmesh.hpp), the base class FreeVol_DualModel

// system of isentropic gas dynamics over the dual mesh (P1 unknowns), 1D (meshL), 2D (mesh), 3D (mesh3)
// the unknowns: U[0] = rho, U[1..d] = rho*u (d components); the fluxes of FreeVol_igd.cpp (3D in the frame of the
// normal; 2D: the operations of the scheme of F. Bouchut), the slopes on rho and u, limited, then limited for rho
// positive and for the velocity ("positivity of internal energy")
class FreeVol_IGDDM: public FreeVol_DualModel
{
public:
    int d;     // dimension, d components of the velocity
    FreeVol_IGDDM(const dualmesh &dm, const vector< FreeVol_FE > &Vfe,int oord,int snumflux,KN<long> labelBCC,int limiterr,KN<double> parameterss,KN<long> intparameterss)
    : FreeVol_DualModel(dm), d(0)
    { init(Vfe,oord,snumflux,labelBCC,limiterr,parameterss,intparameterss); }
    FreeVol_IGDDM(const dualmeshN &dm, const vector< FreeVol_FE > &Vfe,int oord,int snumflux,KN<long> labelBCC,int limiterr,KN<double> parameterss,KN<long> intparameterss)
    : FreeVol_DualModel(dm), d(0)
    { init(Vfe,oord,snumflux,labelBCC,limiterr,parameterss,intparameterss); }

    void init(const vector< FreeVol_FE > &Vfe,int oord,int snumflux,KN<long> labelBCC,int limiterr,KN<double> parameterss,KN<long> intparameterss){
        int nU = Vfe.size();
        int fe = nU ? Vfe[0].dim : 0;// 1 meshL, 2 mesh, 3 mesh3, 4 meshS
        d = (fe == 4) ? 3 : fe;// surface (meshS): 3 components of the velocity, tangent to the surface
        if (fe < 1 || fe > 4 || nU != d+1)
            ExecError("fvscheme igddm: the unknowns are [rho,rhou] in 1D (meshL), [rho,rhou,rhov] in 2D (mesh), [rho,rhou,rhov,rhow] in 3D (mesh3) and [rho,rhou,rhov,rhow] on a surface (meshS)");
        for (int c = 0; c < nU; ++c) {
            if (!Vfe[c].x()) { cout << " Sorry unknown " << c+1 << " is unset "<< endl; ffassert(0);}
            if (Vfe[c].dim != fe || Vfe[c].mesh() != Vfe[0].mesh()) ExecError("fvscheme igddm: the unknowns must be on the same mesh");
        }
        const FreeVol_FE &f = Vfe[0];
        int neq = d+1;//defines the number of equations (unknowns) of the system
        set(neq,oord,snumflux,f,labelBCC,"igddm",limiterr,parameterss,intparameterss);
        if (f.NbOfDF() / f.N() != geom.nv || f.MaxNbDFPerElement() / f.N() != (fe == 4 ? 3 : d+1))
            ExecError("fvscheme igddm: the unknowns must be P1");
        CheckDual("fvscheme igddm");
    }//end of init

    double pressure(double rho) const { return parameters[1]*pow(rho,parameters[0]); }
    double energy(double rho) const { return parameters[1]/(parameters[0]-1.)*pow(rho,parameters[0]-1.); }// p/((gamma-1) rho)

    // Flux F(U).N: F[0] = q.N, F[1+c] = sum on c2 of (q_c u_c2 + p delta) N_c2
    void Flux(long m2LR,const R3 &N,const double *U,double *F) const {
        double h  = U[0];
        for (int j = 0; j < Neq; ++j) F[j] = 0.;
        if (h <= 1.e-13) return;
        double vel[3];
        for (int c = 0; c < d; ++c) vel[c] = U[1+c]/h;
        double p = pressure(h);
        for (int c = 0; c < d; ++c) F[0] += U[1+c]*N[c];
        for (int c = 0; c < d; ++c)
            for (int c2 = 0; c2 < d; ++c2) F[1+c] += (U[1+c]*vel[c2] + (c2 == c ? p : 0.))*N[c2];
    }
    // Local eigenvalues: un-c, un (d-1 times), un+c
    void lambda(long m2LR,const R3 &N,const double *U,double *l) const {
        double h  = U[0];
        for (int j = 0; j < Neq; ++j) l[j] = 0.;
        if (h <= 1.e-13) return;
        double un = 0.;
        for (int c = 0; c < d; ++c) un += U[1+c]/h*N[c];
        double cs = sqrt(parameters[0]*parameters[1]*pow(h,parameters[0]-1.));
        l[0] = un - cs;
        for (int c = 1; c < d; ++c) l[c] = un;
        l[d] = un + cs;
    }
    // the limited slopes gU[j][3*i+c] in the control volumes: on rho and u, then limitation for rho positive and
    // for the velocity ("positivity of internal energy")
    void LimitedSlopes(KN_d *U, KN<double> *gU) const {
        const dualmeshN &D = DM();
        long nv = D.nv;
        KN<double> h(nv),vel[3],alpha(nv);
        KN<double> grad_h,grad_vel[3];
        for(int c=0; c<d; ++c) vel[c].resize(nv);
        for(long i=0; i< nv; ++i){
            h[i] = U[0][i];
            for(int c=0; c<d; ++c) vel[c][i] = U[1+c][i]/max(U[0][i],1.e-13);
        }
        Grad(h,grad_h);// centered slopes in the control volumes
        for(int c=0; c<d; ++c) Grad(vel[c],grad_vel[c]);
        BdrySlope(grad_h,1L);
        for(int c=0; c<d; ++c) BdrySlope(grad_vel[c],1L);
        SlopeLimiterVF(h,grad_h,alpha,limiter);
        for(long k=0; k< 3*nv; ++k) grad_h[k] *= alpha[k/3];
        for(int c=0; c<d; ++c) {
            SlopeLimiterVF(vel[c],grad_vel[c],alpha,limiter);
            for(long k=0; k< 3*nv; ++k) grad_vel[c][k] *= alpha[k/3];
        }
        for(int j=0; j<Neq; ++j) { gU[j].resize(3*nv); gU[j]=0.; }
        for(long i=0; i< nv; ++i){
            double hi = U[0][i];
            if (hi<1.e-13) continue;// vacuum: slopes 0
            // limitation for rho positive
            double alphah = 1.;
            for(long p=D.ifbeg[i];p<D.ifbeg[i+1];++p){//loop on interfaces around control volume i
                double he = hi + IncGE(grad_h,i,D.iflist[p]);
                alphah = min(alphah,1./max((1.-he/max(hi,1.e-13))/0.9,1.));
            }
            for(int c=0; c<3; ++c) grad_h[3*i+c] *= alphah;
            // limitation for "positivity of internal energy"
            double maxvellim = 0.;
            for(long p=D.ifbeg[i];p<D.ifbeg[i+1];++p){
                long m=D.iflist[p];
                double he = hi + IncGE(grad_h,i,m);
                double inc2 = 0.;
                for(int c=0; c<d; ++c) {
                    double incue = IncGE(grad_vel[c],i,m);
                    incue *= hi;
                    incue /= max(he,1.e-13);
                    inc2 += incue*incue;
                }
                maxvellim = max(maxvellim,sqrt(.5*inc2/max(energy(he),1.e-30)));
            }
            double alphavel = 1./max(maxvellim/.5,1.);
            //define the final slopes
            for(int l=0; l<3; ++l) {
                long il = 3*i+l;
                for(int c=0; c<d; ++c) grad_vel[c][il] *= alphavel;
                gU[0][il] = grad_h[il];
                for(int c=0; c<d; ++c) gU[1+c][il] = grad_h[il]*vel[c][i]+grad_vel[c][il]*hi;
            }
        }
    }

    // the scheme on the dual mesh: the flux balance of each control volume and the time step (CFL = 1)
    bool PeriodicOK() const { return true; }// periodic= (PSkip, PeriodicU, PeriodicUp, FreeVolDualMesh.hpp)
    double computedtFlux(KN_d *U,KN_d *Up, KN_d *uBC) const {
        CheckDual("fluxupdate igddm");
        PeriodicU(U);// periodic: the same values in the control volumes of a class
        const dualmeshN &D = DM();
        const long nv=D.nv, npe=D.npe, nde=D.nde;
        for(long i=0; i<nv; ++i)// vacuum: rho < 1e-13 -> rho = rho*u = 0
            if(U[0][i]<1.e-13) { for(int j=0; j<Neq; ++j) U[j][i] = 0.; }
        if (surf()) SurfU(U);// surface: the momentum of each control volume in its plane
        double maxspeed;
        for(int j=0; j<Neq; ++j) Up[j]=0.;
        //left and right states and fluxes for an interface
        double Ul[Neqmax],Ur[Neqmax],phiel[Neqmax],UBCl[Neqmax];
        vector< KN<double> > gU(Neq);// the slopes in the control volumes
        if(ord>1) LimitedSlopes(U,&gU[0]);
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
                if (surf()) SurfLR(ivL,ivR,N,Ul,Ur);// surface: the states and N in the plane of ivL
                NFlux(2*m,2*m+1,N,Ul,Ur,maxspeed,phiel);// numerical flux
            } else {//boundary piece: uBC (P1) at its centroid
                long b=3*(m-npe);
                for(int j=0;j<Neq;++j){
                    UBCl[j]=0.;
                    for(int l=0;l<3;++l) UBCl[j]+=D.bweight[b+l]*uBC[j][D.bvert[b+l]];
                }
                if (surf()) SurfLR(ivL,ivL,N,Ul,UBCl);// surface: the states and N in the plane of ivL
                NFlux_BC(2*m,D.dbelab[m-npe],N,Ul,maxspeed,phiel,UBCl);
            }
            double le=D.area[m];
            if (dual2D()) dt=DtEdge(dt,m,npe,ivL,ivR,cut,D.cvvol,le,maxspeed);// 2D: by dual edge
            else { sumf[ivL]+=le*maxspeed; if(m<npe) sumf[ivR]+=le*maxspeed; }
            // distribution of the fluxes to the control volumes
            double phier[Neqmax];// the flux for ivR: phiel, on a surface its momentum turned into the plane of ivR
            for(int j=0; j<Neq;++j) phier[j] = phiel[j];
            if (surf() && m<npe) RotCV(ivL,ivR,phier+1);
            for(int j=0; j<Neq;++j){
                phiel[j] *= le;
                Up[j][ivL] += phiel[j]*(1./D.cvvol[ivL]);
                if(m<npe) {//if interface is internal
                    phier[j] *= le;
                    Up[j][ivR] -= phier[j]*(1./D.cvvol[ivR]);
                }
            }
        }
        if (!dual2D())// 1D, 3D: min on the control volumes of volume/(sum on its interfaces of area*speed)
            for(long i=0;i<nv;++i) if(!cut[i]) dt=min(dt,D.cvvol[i]/(sumf[i]+1.e-20));
        if (surf()) {// surface: a gas at rest stays at rest
            KN<double> P(nv);
            for(long i=0;i<nv;++i) P[i] = U[0][i]>1.e-13 ? pressure(U[0][i]) : 0.;
            SurfDefect(Up,P);
        }
        PeriodicUp(Up);// periodic: the flux balance of the union of the control volumes of a class
        return dt;
    }//end of computedtFlux
    // numerical flux of an interior face: numflux = 0 HLL, 1 Rusanov, 2 Suliciu/HLLC (as FreeVol_igd.cpp)
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
            else ExecError("FreeVol igddm: numflux not defined");
            return;
        }
        double gamma = parameters[0];
        // numfluxgd returns 4 values (mass, 2 momentum, energy): in F, not in phi (Neq = d+1 values)
        double Np[2],stateL[4],stateR[4],F[4];
        double hl = max(Ul[0],1.e-13), hr = max(Ur[0],1.e-13);
        stateL[0]=Ul[0]; stateL[3]=energy(Ul[0]);
        stateR[0]=Ur[0]; stateR[3]=energy(Ur[0]);
        if(d<=2){// the 2D flux (1D: the velocity on x, tangential velocity 0)
            Np[0]=Nn.x;Np[1]=d==2 ? Nn.y : 0.;
            stateL[1]=Ul[1]/hl; stateL[2]=d==2 ? Ul[2]/hl : 0.;
            stateR[1]=Ur[1]/hr; stateR[2]=d==2 ? Ur[2]/hr : 0.;
            numfluxgd(gamma,stateL,stateR,Np,speed,F);
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
                numfluxgd(gamma,stateL,stateR,Np,speed,F);
                Fm += Tt[t]*F[2];
            }
            Fm += Nn*F[1];
            phi[0]=F[0];
            for(int c=0; c<3; ++c) phi[1+c]=Fm[c];
        }
        if (speed>1.e7) {
            cout << "blowup of speed " << speed << endl;
            cout << "left state " ; for(int j=0;j<Neq;++j) cout << Ul[j] << " "; cout << endl;
            cout << "right state "; for(int j=0;j<Neq;++j) cout << Ur[j] << " "; cout << endl;
            ExecError("FreeVol igddm: blowup of speed");
        }
    }
    // numerical flux of a boundary face: labelBC = 0 Dirichlet (uBC), 1 Neumann (the default), 2 wall
    void NFlux_BC(long m2LR,int labelement,const R3 &Nn,const double *UL,double &speed,double *phi, const double *UBCL) const {
        int choice_BC = 1;// 0 Dirichlet, 1 Neumann(default), 2 Wall
        for(int i=0;i<labelBC.n/2;i++)
            if(labelement==labelBC[2*i]) choice_BC = labelBC[2*i+1];
        double h  = UL[0];
        // Outer state (ghost)
        double UR[Neqmax];
        if(choice_BC == 0){
            for (int j = 0; j < Neq; ++j) UR[j]=UBCL[j];
        } else if(choice_BC == 1){
            for (int j = 0; j < Neq; ++j) UR[j]=UL[j];
        } else if(choice_BC == 2){
            // build the "mirror" state (wall: un reversed, tangential velocity unchanged)
            UR[0] = h;
            Reflect(d,Nn,h,UL+1,UR+1);
        } else {cout<<"choice_BC = "<<choice_BC<<" is not defined for igddm"<<endl;ExecError("FreeVol igddm: labelBC not defined");}
        NFlux(m2LR,m2LR,Nn,UL,UR,speed,phi);
    }
    ~FreeVol_IGDDM() {if(verbosity>99)  cout << " delete FreeVol_IGDDM "<< this << endl;}
};

// creation of the model for fvscheme FV(Thstar, [rho,rhou,rhov], eq="igddm", ...) (see the registry in FreeVol.cpp)
static FreeVol_model *FreeVol_IGDDM_Create(const vector< FreeVol_FE > &V, long order, long numflux,
                                           KN<long> labelBC, long limiter, KN<double> parameters,
                                           KN<long> intparameters,
                                           const void *dual) {
    // the dual mesh built in the script: a dualmeshL* (1D), a dualmesh* (2D), a dualmesh3* (3D), given to fvscheme
    // FV(Thstar, [rho,rhou,rhov], ...), or found from the mesh of [rho,rhou,rhov]: fvscheme FV([rho,rhou,rhov], ...)
    dual = FreeVol_DualModel::FindDual(dual, V, "igddm", "[rho,rhou,rhov]");
    if (numflux == -99) numflux = 2;    // Suliciu/HLLC by default
    // parameters[0] is the exponent gamma of the gas law p=kappa*rho^gamma, parameters[1] is the constant kappa
    if (parameters.N() < 2 || parameters[0] <= 1. || parameters[1] <= 0.) {
        cout << "parameters gamma > 1 and kappa > 0 must be given in parameters=[gamma,kappa]" << endl;
        ExecError("FreeVol igddm: parameters=[gamma,kappa] expected");
    }
    if ((verbosity) && (!mpirank)) {
        cout << "gamma : " << parameters[0] << ", kappa : " << parameters[1] << endl;
        if ((intparameters.N() > 1) || ((intparameters.N() == 1) && intparameters[0] != 0)) {    // intparameters is not yet used
            cout << "intparameters = [ ";
            for (int i = 0; i < intparameters.N() - 1; i++) cout << intparameters[i] << " , ";
            cout << intparameters[intparameters.N() - 1] << " ]" << endl;
        }
        for (int i = 0; i < labelBC.n / 2; i++) {
            if (labelBC[2 * i + 1] == 0) cout << "on boundary " << labelBC[2 * i] << " labelBC=0 : Dirichlet BC" << endl;
            else if (labelBC[2 * i + 1] == 1) cout << "on boundary " << labelBC[2 * i] << " labelBC=1 : Neumann BC" << endl;
            else if (labelBC[2 * i + 1] == 2) cout << "on boundary " << labelBC[2 * i] << " labelBC=2 : Wall BC" << endl;
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
    // a numerical flux not defined for igddm: stop (whatever the verbosity)
    if (numflux < 0 || numflux > 2) {
        cout << "numflux = " << numflux << " is not defined for igddm (0 HLL, 1 Rusanov, 2 Suliciu/HLLC)" << endl;
        ExecError("FreeVol igddm: numflux not defined");
    }
    // the dual mesh given to fvscheme: 2D the dual mesh of F. Bouchut (meshdual), 1D, 3D the median dual mesh
    // (meshdualL, meshdual3)
    int d = V.size() ? V[0].dim : 2;
    if (d == 2)
        return new FreeVol_IGDDM(*(const dualmesh *)dual, V, order, numflux, labelBC, limiter, parameters, intparameters);
    return new FreeVol_IGDDM(*FreeVol_DualModel::DualN(dual, d), V, order, numflux, labelBC, limiter, parameters, intparameters);
}

static void Load_Init() {
    if (verbosity && !mpirank) cout << " load: FreeVol_igddm" << endl;
    FreeVol_RegisterModel("igddm", FreeVol_IGDDM_Create);
}
LOADFUNC(Load_Init)
