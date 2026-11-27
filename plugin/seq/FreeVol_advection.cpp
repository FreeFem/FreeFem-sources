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
// ff-c++ -shared -fPIC FreeVol_advection.cpp -I$FREEFEMPP_DIR/include -L$FREEFEMPP_DIR/lib FreeVol.so
// or simply ff-c++ -shared -fPIC FreeVol_advection.cpp FreeVol.so
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...
/* ---------------------------------------------------------

   finite volume method to solve the advection equation
   U_t + div(V U) = 0, U scalar, V given, in 1D, 2D or 3D:
     1D: meshL (segments), fvscheme FV([V1], ...) (or [V1,V2] / [V1,V2,V3] for a curve)
     2D: mesh (triangles), fvscheme FV([V1,V2], ...)
     3D: mesh3 (tetrahedra), fvscheme FV([V1,V2,V3], ...)
     surface: meshS (triangles in 3D), fvscheme FV([V1,V2,V3], ...), V tangent to the surface: on an edge, the
       normal velocity with the mean of the conormals of the 2 triangles (the same flux for both: conservative)
   the dimension is given by the mesh of V1, V2, ...; each component in P0, P1, or a fespace
   with one dof by face (P0edge or P1nc in 2D); second-order in space by slope reconstruction
   numflux = 0 : HLL, numflux = 1 : Rusanov, numflux = 2 : upwind
   unknown U is P0 (one value by cell)

   --------------------------------------------------------*/

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include "FreeVol.hpp"

// advection by a given vector field, in 1D, 2D or 3D
class FreeVol_Adv : public FreeVol_model
{
public:
    int nV;             // number of components of the velocity (1 to 3)
    mutable KN_d V[3];  // the components of the velocity (their values: read again at each fluxupdate)
    vector< FreeVol_FE > Vf;  // the FE functions of the velocity, given to fvscheme
    long nnV;           // the number of dof by component of the velocity
    int vtype;          // the fespace of the velocity: 0 P0, 1 P1, 2 one dof by face (P0edge, P1nc)
    KN<long> dofV[3];   // vtype 2: dofV[c][kf] the dof of V[c] on the face kf = k*nf+e

    FreeVol_Adv(const vector< FreeVol_FE > &Vfe, int oord,int snumflux,KN<long> labelBCC,int limiterr,KN<double> parameterss,KN<long> intparameterss)
    : FreeVol_model(), nV(0), vtype(0), nnV(0)
    { init(Vfe,oord,snumflux,labelBCC,limiterr,parameterss,intparameterss);}

    void init(const vector< FreeVol_FE > &Vfe, int oord,int snumflux,KN<long> labelBCC,int limiterr,KN<double> parameterss,KN<long> intparameterss){
        nV = Vfe.size();
        int dim = nV ? Vfe[0].dim : 0;
        if (nV < 1 || nV > 3 || (dim == 2 && nV != 2) || (dim >= 3 && nV != 3))
            ExecError("fvscheme advection: the velocity has 1 component in 1D (meshL, or 2, 3 for a curve), 2 in 2D (mesh), 3 in 3D (mesh3) and on a surface (meshS)");
        for (int c = 0; c < nV; ++c)
            if (!Vfe[c].x()) { cout << " Sorry V" << c+1 << " is unset "<< endl; ffassert(0);}
        for (int c = 1; c < nV; ++c)
            if (Vfe[c].mesh() != Vfe[0].mesh()) ExecError("fvscheme advection: the components of the velocity must be on the same mesh");
        int neq = 1;//defines the number of equations (unknowns) of the system
        set(neq,oord,snumflux,Vfe[0],labelBCC,"advection",limiterr,parameterss,intparameterss);
        const FreeVol_Geom &g = geom;
        for(int c=0; c<nV; ++c){
            const FreeVol_FE &f = Vfe[c];
            long nn  = f.NbOfDF() / f.N();                 // nb dofs per component
            int ldf = f.MaxNbDFPerElement() / f.N();       // 1 : P0, dim+1 : P1 or one dof by face
            if((verbosity>4)&&(!mpirank)) cout << "comp=" << f.comp() << " N=" << f.N() << " nn=" << nn << " ldf=" << ldf << endl;
            int t;
            if (ldf == 1 && nn == g.nt) t = 0;                      // P0
            else if (ldf == g.nf && nn == g.nv) t = 1;              // P1 (in 1D, also a dof by face = vertex)
            else if (ldf == g.nf && nn == g.nfaces) {               // one dof by face: P0edge, P1nc
                t = 2;
                dofV[c].resize(g.nt*g.nf);
                for(long k=0; k<g.nt; ++k)
                    for(int e=0; e<g.nf; ++e) dofV[c][k*g.nf+e] = f.dof(k,e);
            }
            else {cout << "fespace must be built on P0 or P1 or a fespace with one dof by face (P0edge, P1nc)" << endl;ffassert(0);}
            if (c == 0) vtype = t;
            else if (t != vtype) {cout << " error: the components of V must be of the same type" << endl;ffassert(0);}
            V[c].set( f.values() );
            if (c == 0) nnV = nn;
        }
        Vf = Vfe;
    }//end of init

    // the values of the velocity: read again before each computation, the velocity can be changed after
    // fvscheme (V1 = ... in the script makes a new array of values: the old one is freed)
    void UpdateV() const {
        for(int c=0; c<nV; ++c){
            if (!Vf[c].x() || Vf[c].NbOfDF() / Vf[c].N() != nnV)
                ExecError("fluxupdate advection: the velocity is not on its fespace of the definition of fvscheme (resetfvscheme)");
            V[c].set( Vf[c].values() );
        }
    }
    // the fluxes and the time step (CFL = 1), the speed of each face given by its numerical flux:
    // 2D: min on the edges of area/(length*speed) (enough in 2D: for a velocity constant in a triangle,
    // the greatest flux of an edge is the sum of the outgoing fluxes);
    // 1D, 3D: min on the cells of measure/sum on the faces of (measure of the face*speed): a
    // tetrahedron can have 2 outgoing faces, each one with half of the outgoing flux
    double computedtFlux(KN_d *U,KN_d *Up, KN_d *uBC) const {
        UpdateV();
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
                        if(g.surf){// surface: the mean of the conormals of the 2 triangles (the same flux for both)
                            R3 Ns = g.N[kf]-g.N[kfo];
                            Ns /= Ns.norme();
                            NFlux(kf,kfo,Ns,Ul,Ur,maxspeed,phie);
                        }
                        else NFlux(kf,kfo,N,Ul,Ur,maxspeed,phie);// numerical flux
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

    double VN(long kf,const R3 &N) const {//normal velocity on the face kf, at each side
        const FreeVol_Geom &g = geom;
        long k = kf/g.nf;
        R3 Vm(0.,0.,0.);//velocity on the face kf
        for(int c=0; c<nV; ++c){
            if (vtype==0) {             // P0
                Vm[c] = V[c][k];
            } else if (vtype==1) {      // P1: mean of the values at the vertices of the face
                double s = 0.;
                for(int l=0; l<g.nvf; ++l) s += V[c][g.fvert[kf*g.nvf+l]];
                Vm[c] = g.nvf==2 ? 0.5*s : s/g.nvf;
            } else {                    // one dof by face
                Vm[c] = V[c][dofV[c][kf]];
            }
        }
        return (Vm,N);
    }
    void Flux(long kf,const R3 &N,const double *U,double *F) const {
        double vn= VN(kf,N);
        F[0]= vn*U[0];
    }
    void lambda(long kf,const R3 &N,const  double *U,double *l) const {
        double vn= VN(kf,N);
        l[0]= vn;
    }
    // numerical flux of an interior face: numflux = 0 HLL, 1 Rusanov, 2 upwind
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
            else ExecError("FreeVol advection: numflux not defined");
            return;
        }
        double lL[Neqmax],lR[Neqmax];
        lambda(kfl,Nn,Ul,lL);
        lambda(kfr,Nn,Ur,lR);
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
        else {cout<<"choice_BC = "<<choice_BC<<" is not defined for advection"<<endl;ExecError("FreeVol advection: labelBC not defined");}
    }
    ~FreeVol_Adv() {if(verbosity>99)  cout << " delete FreeVol_Adv "<< this << endl;}
};

// creation of the model for fvscheme FV([V1,V2], eq="advection", ...) (see the registry in FreeVol.cpp)
static FreeVol_model *FreeVol_Adv_Create(const vector< FreeVol_FE > &V, long order, long numflux,
                                         KN<long> labelBC, long limiter, KN<double> parameters,
                                         KN<long> intparameters,
                                         const void *dual) {
    if (dual) ExecError("fvscheme eq=\"advection\": this model does not take a dual mesh (fvscheme FV([...], eq=\"advection\", ...))");
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
        cout << "numflux = " << numflux << " is not defined for advection (0 HLL, 1 Rusanov, 2 Upwind)" << endl;
        ExecError("FreeVol advection: numflux not defined");
    }
    return new FreeVol_Adv(V, order, numflux, labelBC, limiter, parameters, intparameters);
}

static void Load_Init() {
    if (verbosity && !mpirank) cout << " load: FreeVol_advection" << endl;
    FreeVol_RegisterModel("advection", FreeVol_Adv_Create);
}
LOADFUNC(Load_Init)
