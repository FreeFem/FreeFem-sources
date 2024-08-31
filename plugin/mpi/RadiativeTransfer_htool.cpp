//ff-c++-LIBRARY-dep: mpi pthread htool
#include <htool/htool.hpp>
#include <ff++.hpp>
#include <AFunction_ext.hpp>

using namespace htool;

double dic = 0.03; // make dic=1 if kappa is constant to speed up

class KappaGrid {
public:
    KNMK<float> * tab;

    int nx, ny, nz;
    double xend, xstart;
    double yend, ystart;
    double zend, zstart;
    double dx, dy, dz;

    void init() {tab = 0;}
    void destroy() {delete tab;}
};

KappaGrid *init_KappaGrid(Stack stack, KappaGrid *const &a, KN<double> *const& bounds, double const& h, Expression kappa0) {
    KN<double>& bndsg = *bounds;

    int ixu = (int)ceil((bndsg[1]-bndsg[0])/h);
    int iyu = (int)ceil((bndsg[3]-bndsg[2])/h);
    int izu = (int)ceil((bndsg[5]-bndsg[4])/h);

    a->nx = ixu+1;
    a->ny = iyu+1;
    a->nz = izu+1;

    int sz = a->nx * a->ny * a->nz;
    a->xstart = bndsg[0];
    a->xend = bndsg[0] + ixu*h;
    a->dx = (a->xend - a->xstart)/(a->nx-1);
    a->ystart = bndsg[2];
    a->yend = bndsg[2] + iyu*h;
    a->dy = (a->yend - a->ystart)/(a->ny-1);
    a->zstart = bndsg[4];
    a->zend = bndsg[4] + izu*h;
    a->dz = (a->zend - a->zstart)/(a->nz-1);

    a->tab = new KNMK< float >(a->nx, a->ny, a->nz);

    MeshPoint* mp(Fem2D::MeshPointStack(stack));

    int ix, iy, iz;
    double xs, ys, zs;
    for (iy = 0; iy < a->ny; iy++) {
        ys = a->ystart + a->dy * iy;
        for (ix = 0; ix < a->nx; ix++) {
            xs = a->xstart + a->dx * ix;
            for (iz = 0; iz < a->nz; iz++) {
                zs = a->zstart + a->dz * iz;
                mp->set(xs, ys, zs);
                (*a->tab)(ix, iy, iz) = GetAny<double>( (*kappa0)(stack) );
            }
        }
    }

    return a;
}

double KappaGrid_eval(KappaGrid *const &a, const double &xi, const double &yi, const double &zi) {
    int ix = (a->nx-1) * (xi - a->xstart + a->dx/2) / (a->xend - a->xstart);
    int iy = (a->ny-1) * (yi - a->ystart + a->dy/2) / (a->yend - a->ystart);
    int iz = (a->nz-1) * (zi - a->zstart + a->dz/2) / (a->zend - a->zstart);
    ix = max(0, min(ix, a->nx - 1));
    iy = max(0, min(iy, a->ny - 1));
    iz = max(0, min(iz, a->nz - 1));

    return (*a->tab)(ix, iy, iz);
}

//////////////////////////// used for the volume integral ////////////////////////////
class Generator_Volume: public VirtualGenerator<double>{
public:
    const Mesh3& Th;
    KN<int> heade,nexte;
    KappaGrid* kappa0;

    typedef SortArray<int,2> SortArray2;
    mutable HashTable<SortArray2,int> edges;
    const int e2[10][2] = {{ 0,0},{ 0,1},{ 0,2},{ 0,3}, {1,1}, {1,2}, {1,3}, {2,2},{2,3}, {3,3}};

    Generator_Volume(pmesh3 pth3, KappaGrid* k):
    VirtualGenerator<double>(), Th(*pth3), kappa0(k),
        edges(pth3->nt*3+pth3->nv,pth3->nv), heade(), nexte(pth3->nt*10) {
        // pour i,j -> liste de tet de sommet i,j
        // i,j -> liste des tet contenant i et j  ??
        // k, e ->  i,j -> numero d'arete
        // inverse tableau i,j -> e ->
        int ne=0;
        for(int k=0; k<Th.nt; ++k)
        for(int e=0; e<10; ++e) {
            int i = Th(k,e2[e][0]);
            int j = Th(k,e2[e][1]);
            SortArray2 ij(i,j);
            auto nKV = edges.find(ij);
            if( nKV==0) nKV = edges.add(ij, ne++); // new
        }
        if (mpirank == 0 && verbosity) cout << " nb edge + sommets = "<< ne << endl;

        heade.resize(ne);
        heade = -1;
        // build change
        for(int k=0; k<Th.nt; ++k)
        for(int e=0; e<10; ++e) {
            int i = Th(k,e2[e][0]);
            int j = Th(k,e2[e][1]);
            SortArray2 ij(i,j);
            auto nKV = edges.find(ij);
            ffassert(nKV);
            int ie = nKV->v;
            int ke = k*10+e; // ke -> e
            int h = heade[ie];
            nexte[ke] = h;
            heade[ie] = ke;
        }
    }

    double get_coef(const int& i, const int& j) const {
        using Fem2D::R3;
        const double fourpi = Pi*4.;
        // pour les tet K de sommet i;j  ik;jk num de K
        // choisir le formule en fonct de le dist de i,j

        R3 I = Th(i); // vertex i
        R3 J = Th(j); // vertex j
        R3 IJ = R3(I,J);
        double lIJ2 = IJ.norme2();
        double lIJ = sqrt(lIJ2);

        int exact = 5; // choose the quadrature formula
        SortArray2 jj(j,j);
        SortArray2 ij(i,j);
        auto nKV = edges.find(jj);
        auto nKVij = edges.find(ij);
        if (nKVij == 0) exact = 2; // i and j not in same element -> use lower order quadrature
        const GQuadratureFormular<R3> * pQF = QF_Simplex<R3> (exact), QF = *pQF;
        int ie = nKV->v;

        double kappa_ij = 0; // mean value of kappa on the (i,j) line segment

        if (i == j) {
            kappa_ij = KappaGrid_eval(kappa0, I[0], I[1], I[2]);
        }
        else {
            int cpt = 0;
            // quadrature over segment IJ with step size dic to compute kappa_ij
            for(double ic = 0; ic<lIJ; ic+=dic, cpt++) {
                R3 aux = I + ic/lIJ*IJ;
                kappa_ij += KappaGrid_eval(kappa0, aux[0], aux[1], aux[2]);
            }
            kappa_ij /= cpt;
        }

        double a_ij = 0;

        // loop over all elements in the support of basis function j to compute the integral
        for (int k10=heade[ie]; k10>=0; k10 = nexte[k10]) {
            int k = k10/10;
            int e = k10%10;
            int jk = e2[e][0]; // local index of j in element k

            double l3[4];
            const Mesh3::Element & K = Th[k];
            const double ck = K.mesure()/fourpi;
            // quadrature loop
            for(int p=0; p < QF.n; ++p) {
                double cc = QF[p].a*ck;
                R3 Qp = K(QF[p]); // the quadrature point
                QF[p].toBary(l3); // barycentric coordinates of the quadrature point
                R3 IQ = R3(I,Qp);
                double lIQ2 = IQ.norme2();
                double lIQ = sqrt(lIQ2);
                double wj = l3[jk];
                double kappa_j = KappaGrid_eval(kappa0, Qp[0], Qp[1], Qp[2]);
                a_ij += cc * wj * kappa_j * exp(-kappa_ij*lIQ)/(lIQ2);
            }
        }

        return a_ij;
    }

    void copy_submatrix(int m, int n, const int *const rows, const int *const cols, double *ptr) const {
        std::fill_n(ptr,m*n,0);
        for (int i=0; i<m; i++)
        for (int j=0; j<n; j++) {
            ptr[i+m*j] = get_coef(rows[i], cols[j]);
        }
    }
    ~Generator_Volume(){}
};

VirtualGenerator<double>** init_Generator_Volume(VirtualGenerator<double>** const &  ppf, 
                                                 pmesh3 const & pTh, KappaGrid * const & kappa0) {
    *ppf = dynamic_cast<VirtualGenerator<double>*> (new Generator_Volume(pTh,kappa0));
    return ppf;
}

//////////////////////////// used for the boundary integral ////////////////////////////
class Generator_Boundary: public VirtualGenerator<double>{
public:
    const Mesh3 & Th3;
    const MeshS & ThS;
    KN_<double> seeface;
    KappaGrid * kappa0;
    KN<int> headv, nextv;

    Generator_Boundary(pmesh3 pth3, pmeshS pthS, KN_<double> see, KappaGrid* k):
    VirtualGenerator<double>(), Th3(*pth3), ThS(*pthS), seeface(see),
        kappa0(k), headv(pthS->nv,-1), nextv(pthS->nt*3) {

        for(int k=0; k<ThS.nt; ++k)
        for(int v=0; v<3; ++v) {
            int j = ThS(k,v);
            int kv = k*3+v;
            int h = headv[j];
            nextv[kv] = h;
            headv[j] = kv;
        }
    }

    double get_coef(const int& i, const int& j) const {
        using Fem2D::R3;
        const double fourpi = Pi*4.;

        R3 I = Th3(i);
        R3 J = ThS(j);
        R3 IJ = R3(I,J);
        //if ((IJ,N) <= 0) return 0.;
        double lIJ2 = IJ.norme2();
        double lIJ = sqrt(lIJ2);

        int exact = 25; // choose the quadrature formula
        const GQuadratureFormular<Fem2D::R2> * pQF= QF_Simplex<Fem2D::R2>(exact), QF = *pQF;

        double kappa_ij = 0; // mean value of kappa on the (i,j) line segment

        int cpt = 0;
        // quadrature over segment IJ with step size dic to compute kappa_ij
        for(double ic = 0; ic<lIJ; ic+=dic, cpt++) {
            R3 aux = I + ic/lIJ*IJ;
            kappa_ij += KappaGrid_eval(kappa0, aux[0], aux[1], aux[2]);
        }
        kappa_ij /= cpt;

        double a_ij = 0;

        // loop over all elements in the support of basis function j to compute the integral
        for (int k3=headv[j]; k3>=0; k3=nextv[k3]) {
            int k = k3/3;
            int v = k3%3;
            int jjk = ThS(k,v);  // local index of j in element k

            double l3[3];
            const MeshS::Element & K = ThS[k];
            const double ck = K.mesure()/fourpi;
            // quadrature loop
            for(int p=0; p < QF.n; ++p) {
                double cc = QF[p].a*ck;
                R3 Qp = K(QF[p]); // the quadrature point
                QF[p].toBary(l3); // barycentric coordinates of the quadrature point
                R3 IQ = R3(I,Qp);
                double lIQ2 = IQ.norme2();
                double lIQ = sqrt(lIQ2);
                double wj = l3[v];
                a_ij += cc * wj * seeface[j] * pow(min((IQ,-1*K.NormalTUnitaire()),0.)/lIQ2, 2)
                        * exp(-kappa_ij*lIQ);
            }
        }

        return a_ij;
    }

    void copy_submatrix(int m, int n, const int *const rows, const int *const cols, double *ptr) const {
        std::fill_n(ptr,m*n,0);
        for (int i=0; i<m; i++)
        for (int j=0; j<n; j++) {
            ptr[i+m*j] = get_coef(rows[i], cols[j]);
        }
    }
    ~Generator_Boundary(){}
};

VirtualGenerator<double> **init_Generator_Boundary(VirtualGenerator<double>** const &  ppf,
    pmesh3 const & pth3, pmeshS const & pthS, KN<double>* const & seeface,
    KappaGrid * const & kappa0) {
    *ppf = dynamic_cast<VirtualGenerator<double>*> (new Generator_Boundary(pth3,pthS,*seeface,kappa0));
    return ppf;
}

template<class R, class A0, class A1, class A2, class A3, class E=E_F0> // extend (4th arg.)
class E_F_F0F0F0es_ : public E { 
    public:
    typedef R (*func)(Stack, const A0 &, const A1 &, const A2 &, Expression); // extend (4th arg.)
    func f;
    Expression a0,a1,a2,a3; // extend
    E_F_F0F0F0es_(func ff,
        Expression aa0,
        Expression aa1,
        Expression aa2,
        Expression aa3)   // extend
    : f(ff), a0(aa0), a1(aa1), a2(aa2), a3(aa3) {} // extend (4th arg.)
    AnyType operator()(Stack s) const {
        return SetAny<R>( f(s, GetAny<A0>((*a0)(s)),
            GetAny<A1>((*a1)(s)),
            GetAny<A2>((*a2)(s)),
            a3 ) ); // extend (4th arg.)
    }

    virtual size_t nbitem() const {return a3->nbitem();} // modif
    bool MeshIndependent() const {
        return E::MeshIndependent() && a0->MeshIndependent() && a1->MeshIndependent() && a2->MeshIndependent()
            && a3->MeshIndependent();
    } 
};

template<class R, class A=R, class B=A, class C=B, class D=C , class CODE=E_F_F0F0F0es_<R,A,B,C,D,E_F0> >
class OneOperator3es_ : public OneOperator { // 3->4
    aType r; // return type
    typedef typename CODE::func func;
    func f;
    public:
    E_F0 * code(const basicAC_F0 & args) const {
        if (args.named_parameter && !args.named_parameter->empty())
            CompileError("They are used Named parameter");

        return new CODE(f,
            t[0]->CastTo(args[0]),
            t[1]->CastTo(args[1]),
            t[2]->CastTo(args[2]),
            t[3]->CastTo(args[3])); // extend
    }
    OneOperator3es_(func ff) : // 3->4
    OneOperator(map_type[typeid(R).name()],
        map_type[typeid(A).name()],
        map_type[typeid(B).name()],
        map_type[typeid(C).name()],
        map_type[typeid(D).name()]), // extend
    f(ff) {}
};

double Bnu(double nu, double T) {
    return nu*nu*nu/(exp(nu/T)-1);
}

static void Init_RT() {
    Global.Add("Bnu","(",new OneOperator2<double>(Bnu));

    Dcl_Type< KappaGrid * >(InitP< KappaGrid >, Destroy< KappaGrid >);
    zzzfff->Add("KappaGrid", atype< KappaGrid * >( ));

    TheOperators->Add("<-", new OneOperator3_<VirtualGenerator<double> **, VirtualGenerator<double> **,
    pmesh3, KappaGrid *>(init_Generator_Volume));
    TheOperators->Add("<-", new OneOperator5_<VirtualGenerator<double> **, VirtualGenerator<double> **,
    pmesh3, pmeshS, KN<double>*, KappaGrid *>(init_Generator_Boundary));

    TheOperators->Add(
        "<-", new OneOperator3es_<KappaGrid *, KappaGrid *, KN<double>*, double, double>(init_KappaGrid));
    atype< KappaGrid * >()->Add(
        "(", "", new OneOperator4_<double, KappaGrid *, double, double, double>(KappaGrid_eval));
}

LOADFUNC(Init_RT)
