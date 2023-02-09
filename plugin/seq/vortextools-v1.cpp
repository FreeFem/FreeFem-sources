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
// AUTHORS : V. KALT, G.SADAKA, I. DANAILA, F. HECHT
// PAPER : Identification of vortices in quantum fluids: finite element algorithms and programs
// E-MAIL  : ...

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

// to compile: ff-c++ vortextools.cpp
// WARNING: do not compile under windows

#include "ff++.hpp"
#include "AFunction_ext.hpp"
#include <array>

R3 get(KNM<double> &b,int i){
    return R3(b(i,0),b(i,1),b(i,2));
}

// Builds a Bspline with 4 control points evaluates it at a given value
// Inputs:
//    p0: first control point
//    p1: second control point
//    p2: third control point
//    p3: fourth control point
//    t: value at which the Bspline is evaluated
// Returns:
//    r: the value of the Bspline at parameter t
R3  BSp(R3 &p0,R3 &p1,R3 &p2, R3 &p3,double t){
    double t2=t*t, t3=t2*t, tm13=pow((1.-t),3);
    R3 r = tm13/6. * p0 + (3.*t3-6.*t2+4.)/6. * p1 + (-3.*t3+3.*t2+3.*t+1.)/6. * p2 + t3/6. * p3;
    return r;
}

// Evaluates the derivative of the Bspline at a given value
// Inputs:
//    p0: first control point
//    p1: second control point
//    p2: third control point
//    p3: fourth control point
//    t: value at which the Bspline is evaluated
// Returns:
//    r: the value of the Bspline derivative at parameter t
R3  dBSp(R3 &p0,R3 &p1,R3 &p2, R3 &p3,double t){
    double t2=t*t, tm12=pow((1.-t),2);
    R3 r = - tm12/2. * p0 + (3.*t2-4.*t)/2. * p1 + (-3.*t2+2.*t+1.)/2. * p2 + t2/2. * p3;
    return r;
}

// Evaluates the second derivative of the Bspline at a given value
// Inputs:
//    p0: first control point
//    p1: second control point
//    p2: third control point
//    p3: fourth control point
//    t: value at which the Bspline is evaluated
// Returns:
//    r: the value of the Bspline second derivative at parameter t
R3  ddBSp(R3 &p0,R3 &p1,R3 &p2, R3 &p3,double t){
    R3 r = (1.-t) * p0 + (3.*t-2.) * p1 + (-3.*t+1.) * p2 + t * p3;
    return r;
}

// Evaluates the third derivative of the Bspline at a given value
// Inputs:
//    p0: first control point
//    p1: second control point
//    p2: third control point
//    p3: fourth control point
//    t: value at which the Bspline is evaluated
// Returns:
//    r: the value of the Bspline third derivative at parameter t
R3  dddBSp(R3 &p0,R3 &p1,R3 &p2, R3 &p3,double t){
    R3 r = - p0 + 3. * p1 - 3. * p2 + p3;
    return r;
}

// Computes the curvature kappa and the torsion tau of the Bspline at parameter value t
// Inputs:
//    p0: first control point
//    p1: second control point
//    p2: third control point
//    p3: fourth control point
//    t: value at which the Bspline is evaluated
// Returns:
//    A R2 vector containing the curvature and the torsion of the Bspline at parameter t
R2 kappatau(R3 &P0,R3 &P1,R3 &P2,R3 &P3,double t)
{
    R3 Sp, Spp, Sppp, SpxSpp;
    Sp = dBSp(P0,P1,P2,P3,t);
    Spp = ddBSp(P0,P1,P2,P3,t);
    Sppp = dddBSp(P0,P1,P2,P3,t);
    SpxSpp = Sp^Spp;
    double norm2=SpxSpp.norme2(),norm=sqrt(norm2);
    return R2(norm/pow(Sp.norme(),3),(SpxSpp,Sppp)/norm2);
}

// Smooth a Curve through a 5-point moving average, apply niter smoothing iterations
// Inputs:
//    pc: list of points representing the curve
//    niter: number of smoothing iterations
// Returns:
//    0
long smoothCurve(KNM<double> * const &pc,long const & niter)
{
    double omega=.8,omega1=1.-omega;
    KNM<double> &p=*pc; // pour retirer le const ..
    long n= p.N()-1;
    long m = p.M();
    ffassert(m==3);
    
    KNM<double> q(p); // copy p;
    bool loop=false;
    if(p(0,0)==p(n,0L) && p(0,1)==p(n,1L) && p(0,2)==p(n,2L)) loop = true;// Check if the curve is closed
    
    for(int iter=0; iter<niter; ++iter){
        for(int j = 0; j<m;++j){
            for(int i=0;i<=n;i++){
                if(loop){
                    if(i<n){ // If the curve is closed, apply moving average on all points
                        int ipp1 = (i+n-2)%n;
                        int ip1 = (i+n-1)%n;
                        int imm1 = (i+n+2)%n;
                        int im1 = (i+n+1)%n;
                        q(i,j) = 1./5. * (p(ipp1,j) + p(ip1,j) + p(i,j) + p(im1,j) + p(imm1,j));
                    } else {
                        q(i,j) = q(0,j);
                    }
                } else { // else do not change end points and use a 3-point average for i=1 and i=N-1
                    if(i==0 | i==n){
                        q(i,j) = p(i,j);
                    } else if(i==1 | i==(n-1)) {
                        q(i,j) = 1./3. * (p(i-1,j) + p(i,j) + p(i+1,j));
                    } else {
                        q(i,j) = 1./5. * (p(i-2,j) + p(i-1,j) + p(i,j) + p(i+1,j) + p(i+2,j));
                    }
                }
            }
        }
        p *= omega1;
        q *= omega;
        p += q;
    }
    return 0L;
}

// Interpolate a curve using a bspline
// Input: 
//    pc: a list of points defining the curve
//    np: the number of points added between existing points
// Returns:
//    pk: the new list of points
KNM_<double> BSp(Stack stack,KNM_<double> const & pc, long const & np)
{
    KNM_<double> p = pc;
    long n = p.N()-1;
    long m = p.M();
    long nbp = n;
    long nbp2;
    long ntvec = max(2L,np);
    double dt = 1./ ntvec;
    ffassert(m==3);
    
    bool loop=false;// Check if the curve is closed
    if(p(0,0)==p(n,0L) && p(0,1)==p(n,1L) && p(0,2)==p(n,2L)) loop = true;
    
    // Endpoints need to be repeated mutliple times to ensure they do not change when building the Bspline
    KNM<double> b(1,m);
    if(loop){
        nbp2 = nbp+3;
        b.resize(nbp2,3);
        b(0,':') = p(nbp-2,':');
        b(1,':') = p(nbp-1,':');
        for(int i=0;i<=nbp;i++){
            b(i+2,':') = p(i,':');
        }
    } else {
        nbp2 = nbp+5;
        b.resize(nbp2,3);
        b(0,':') = p(0,':');
        b(1,':') = p(0,':');
        for(int i=0;i<nbp2-2;i++){
            b(i+2,':') = p(i,':');
        }
        b(nbp2-2,':') = p(nbp,':');
        b(nbp2-1,':') = p(nbp,':');
    }
    int newnbp = nbp2-3;
    int taille = newnbp*(ntvec-1)+1;
    int npk=6;
    KNM_<double> pk(Add2StackOfPtr2FreeA(stack,new double[taille*npk]), taille, npk);// pour les PB de memoire ...
    int ik=0;
    R2 kt;
    R3 bsp, bsp_old;
    double s_tot = 0.;
    for(int i=0;i<newnbp;i++)
    {
        R3 P0(get(b,i)),P1(get(b,i+1)),P2(get(b,i+2)),P3(get(b,i+3));//Existing points
        for(int ni=0;ni<ntvec-1;ni++)// add new points
        {
            double t = ni*dt;
            bsp = BSp(P0,P1,P2,P3,t);
            if(i != 0) s_tot = s_tot + (bsp - bsp_old).norme();
            bsp_old = bsp;
            pk(ik,0) = bsp.x;
            pk(ik,1) = bsp.y;
            pk(ik,2) = bsp.z;
            pk(ik,3) = s_tot;
            kt = kappatau(P0,P1,P2,P3,t); // Compute curvature and torsion
            pk(ik,4) = kt.x;
            pk(ik,5) = kt.y;
            ik++;
        }
        if(i==newnbp-1){
            bsp = BSp(P0,P1,P2,P3,1.);
            s_tot = s_tot + (bsp - bsp_old).norme();
            pk(ik,0) = bsp.x;
            pk(ik,1) = bsp.y;
            pk(ik,2) = bsp.z;
            pk(ik,3) = s_tot;
            kt = kappatau(P0,P1,P2,P3,1.);
            pk(ik,4) = kt.x;
            pk(ik,5) = kt.y;
            ik++;
        }
    }
    if(!loop){
        int ind=5*(ntvec-1);
        for(int ni=0;ni<=ind;ni++){
            pk(ni,4)=pk(ind+1,4);
            pk(ni,5)=pk(ind+1,5);
            pk(ik-1-ni,4)=pk(ik-1-ind-1,4);
            pk(ik-1-ni,5)=pk(ik-1-ind-1,5);
        }
    }
    return pk;
}

// Computes the position of a vortex point in a triangle
// Inputs: 
//    u: Array of wavefunction values on triangle vertices
// Returns:
//    A R2 vector containing the coordinates of the zero
R2 zero(Complex u[3])
{
    typedef Complex C;
    C O(0.,0.);
    auto det=[](C &a, C & b, C &c) { C ab = b-a, ac = c-a; return ab.real()*ac.imag() - ab.imag() * ac.real();};
    double d= det(u[0],u[1],u[2]);
    double x= det(u[0],O,u[2]);
    double y= det(u[0],u[1],O);
    return R2(x/d,y/d);
}

// Computes the position of a vortex point in a triangle
// Inputs: 
//    u0: wavefunction value on the first vertex of the triangle
//    u1: wavefunction value on the second vertex of the triangle
//    u2: wavefunction value on the third vertex of the triangle
// Returns:
//    A R3 vector containing the coordinates of the zero
R3 * zero3(Stack stack,Complex const & u0,Complex const & u1,Complex const & u2)
{
    Complex u[]={u0,u1,u2};
    R3 Q;
    Q=zero(u);
    cout << " P= "<< Q << endl;
    R3 *pQ = Add2StackOfPtr2Free(stack, new R3(Q));
    return pQ;
}

// Check if point P is inside a triangle
// Inputs: 
//    P: Coordinates of the point
//    eps: tolerance value for the search of the vortex points
// Returns:
//    b: A boolean indicating if the point is inside the triangle
bool in( R2 P,double eps)
{
    return  (P.x > -eps) && (P.y > -eps) && (1.-P.x-P.y > -eps) ;
}

// Check if a vortex point exists in a triangle
// Inputs: 
//    u: Array of wavefunction values on triangle vertices
// Returns:
//    b: A boolean indicating the presence of a vortex point
bool in(Complex u[3])
{
    double xmin = min(min(u[0].real(),u[1].real()),u[2].real());
    double xmax = max(max(u[0].real(),u[1].real()),u[2].real());
    double ymin = min(min(u[0].imag(),u[1].imag()),u[2].imag());
    double ymax = max(max(u[0].imag(),u[1].imag()),u[2].imag());
    bool b = (xmin < 0.) && (xmax > 0.) &&  (ymin < 0.) && (ymax > 0.);
    return b;
}

// Check if a vortex point exists in a triangle
// Search vortices in a 2D complex field:
// Inputs: 
//    u: Array of wavefunction values on triangle vertices
//    eps: tolerance value for the search of the vortex points
// Outputs:
//    P: Vortex coordinates
// Returns:
//    b: A boolean indicating if the vortex point is inside the triangle
bool in(Complex u[],R2 & P,double eps)
{ // bof bof on value of eps !!!! 1e-15
    // possible 2 value one on u and one barycentrique ???
    // because not the same range ...
    double xmin = min(min(u[0].real(),u[1].real()),u[2].real());
    double xmax = max(max(u[0].real(),u[1].real()),u[2].real());
    double ymin = min(min(u[0].imag(),u[1].imag()),u[2].imag());
    double ymax = max(max(u[0].imag(),u[1].imag()),u[2].imag());
    bool b = (xmin < eps) && (xmax > -eps) &&  (ymin < eps) && (ymax > -eps);
    if( b )//Compute the vortex position and check if it is inside the triangle
    {
        P = zero(u);
        b = in(P,eps);// check if the point are inside or on border of the triangle
    }
    else P=R2(-1.,-1.); // outside
    return b;
}

// Search vortices in a 3D complex field:
// Inputs: 
//    fu: mesh and complex P1 field representing the wavefunction
//    eps: tolerance value for the search of the vortex points
// Outputs:
//    fuc: P0 function: =1 in tetrahedrons containing a vortex point, =0 otherwise
// Returns:
//    nkk: the number of tetrahedrons crossed by vortex lines
long uZero(pf3c const & fu, pf3r const & fuc, double const &eps)
{
    typedef Mesh3::Element Element;
    typedef   v_fes3::FESpace FESpace;
    typedef Complex K;
    pf3cbase bu=fu.first;
    pf3rbase buc=fuc.first;
    
    ffassert(fu.second==0);
    ffassert(fuc.second==0);
    
    KN<K> * pu=bu->x();
    KN<double> * puc=buc->x();
    
    FESpace *pUh = fu.first->Vh ;
    FESpace *pUch = fuc.first->newVh( );
    ffassert(pUh && pUch);
    ffassert(&pUh->Th ==  &pUch->Th );
    if(pu ==0 || pu->N() != pUh->NbOfDF) {
        ffassert(0); //  u undef ...
    }
    if(puc ==0 || puc->N() != pUch->NbOfDF) {
        if(!mpirank && verbosity>0) cout << "  FE create or recreate " << puc <<  endl;
        if(puc) delete [] puc;
        *fuc.first = puc = new KN< double >(pUch->NbOfDF);
        *puc = double( );
    }
    const Mesh3 & Th=pUh->Th;
    KN<K> & u=*pu;
    KN<double> & uc=*puc;
    FESpace & Uh = *pUh;
    int nkk=0;
    const double twopi = 2.*Pi;
    double charge;
    uc=0.;
    for(int k=0; k<Th.nt;++k)
    {
        const Element & K = Th[k];
        int i0 = Uh(k,0);
        int i1 = Uh(k,1);
        int i2 = Uh(k,2);
        int i3 = Uh(k,3);
        
        double xmn = min(min(min(u[i0].real(),u[i1].real()),u[i2].real()),u[i3].real());
        double xmx = max(max(max(u[i0].real(),u[i1].real()),u[i2].real()),u[i3].real());
        double ymn = min(min(min(u[i0].imag(),u[i1].imag()),u[i2].imag()),u[i3].imag());
        double ymx = max(max(max(u[i0].imag(),u[i1].imag()),u[i2].imag()),u[i3].imag());
        
        bool b = (xmn < eps) && (xmx > -eps) &&  (ymn < eps) && (ymx > -eps);
        if(b)
        {
            int ucn=0;// number of faces of charge = 1
            for(int i=0; i<4;++i)
            {
                int i0 = Th(K[Element::nvface[i][0]]);
                int i1 = Th(K[Element::nvface[i][1]]);
                int i2 = Th(K[Element::nvface[i][2]]);
                Complex u2[3]={u[i0],u[i1],u[i2]};
                R2 P2;
                if(in(u2,P2,eps)){
                    Complex icharge = (log(u2[1]/u2[0]) + log(u2[2]/u2[1]) + log(u2[0]/u2[2]))/twopi;
                    charge = imag(icharge);
                    if(isnan(real(icharge)) || isnan(imag(icharge))){// We divide by zero
                        charge=1.;
                    }
                    if(abs(charge)>.98) ucn++;// compute for each tetraedra the number of face of abs(charge) == 1
                }
            }
            uc[k]=ucn;
            nkk++;
        }
    }
    return nkk;
}

// Search vortices in a 2D complex field:
// Inputs: 
//    pTh, pu: mesh and complex P1 field representing the wavefunction
// Outputs:
//    ppoints: array storing vortex coordinates
//    pucharge: P0 function: =1 in triangles containing a vortex point, =0 otherwise
//    pdmin: minimum distance between the vortices
// Returns:
//    nbc: number of vortices
long uZero2D(const Mesh * const & pTh,KNM<double>*const &ppoints,KN<Complex>*const &pu, KN<double>*const &pucharge, double* const & pdmin)
{
    typedef Mesh::Element Element;
    typedef Element::Vertex Vertex;
    KN<Complex> &u = *pu;
    KN<double> &ucharge = *pucharge;
    const Mesh &Th = *pTh;
    KNM<double> &pts=*ppoints;
    KNM<double> points(Th.nt,2);
    double &dmin = *pdmin;
    ffassert(u.N()==Th.nv);
    ffassert(ucharge.N()==Th.nt);
    const double twopi = 2.*Pi;
    Complex u0,u1,u2;// value of u over the vertex of the tetrahedron
    double charge;
    double l0,l1,l2;
    int nbc=0;
    for (int k=0; k<Th.nt; k++){
        const Element & K = Th[k];
        int i0 = Th(k,0);
        int i1 = Th(k,1);
        int i2 = Th(k,2);
        Complex u2[3]={u[i0],u[i1],u[i2]};
        R2 P2;
        if(in(u2,P2,1.e-15)){
            Complex icharge = (log(u2[1]/u2[0]) + log(u2[2]/u2[1]) + log(u2[0]/u2[2]))/twopi;
            charge = imag(icharge);
            if(isnan(real(icharge)) || isnan(imag(icharge))){// We divide by zero
                charge=1.;
            }
            if(abs(charge)>.98){
                l1=P2.x; l2=P2.y; l0 = 1-l1-l2;
                points(nbc,0) = Th(i0).x * l0 + Th(i1).x * l1 +  Th(i2).x * l2;
                points(nbc,1) = Th(i0).y * l0 + Th(i1).y * l1 +  Th(i2).y * l2;
                ucharge[k] = 1.;
                nbc++;
            }
        }
        // remove duplicated points
        for(int i=0;i<nbc-1;i++) if(points(i,0) == points(nbc-1,0) && points(i,1) == points(nbc-1,1)) nbc--;
    }
    points.resize(nbc,2);
    pts.resize(nbc,2);
    pts=points;
    KN<Vertex> Pf(nbc);
    for (int i=0; i< nbc;++i)
    {   Pf[i].x =points(i,0);
        Pf[i].y =points(i,1);
    }
    R2 Pmin,Pmax;
    Th.BoundingBox(Pmin,Pmax);
    FQuadTree *gtree=new FQuadTree(Pf,Pmin,Pmax,0);
    double mindist=1e100;
    gtree->Add(Pf[0]);
    for (int i=1; i<nbc; i++){
        const Vertex * pvi=gtree->TrueNearestVertex(Pf[i]);
        R2 d(Pf[i],*pvi);
        mindist = min(mindist,d.norme2());
        gtree->Add(Pf[i]);
    }
    dmin = sqrt(mindist);
    delete gtree;
    return (long)nbc;
}

// Build a graph of vortex points from a complex wavefunction
// Inputs:
//    fu: the wavefunction
//    eps: tolerance value for the search of the vortex points
// Outputs:
//    ppoints: list of vortex points
//    pbe: list of line begining and end
//    ploop: list of boolean indicating if the lines are closed or not
// Returns:
//    b: the number of vortex lines
long ZeroLines(pf3c const & fu,double const & eps, KNM<double>*const &ppoints,KN<long>*const &pbe,KN<long>* const & ploop)
{
    
    typedef   v_fes3::FESpace FESpace;
    typedef Complex K;
    pf3cbase bu=fu.first;
    
    ffassert(fu.second==0);
    KN<K> * pu=bu->x();
    FESpace *pUh = fu.first->Vh ;
    const Mesh3 & Th=pUh->Th;
    KN<K> & u=*pu;
    R3 D=Th.Pmax-Th.Pmin;
    R3 M =(Th.Pmax+Th.Pmin)*.5;
    double Dmesh = max(max(D.x,D.y),D.z);
    double lb = Dmesh*0.6;
    R3 Bmin = M - R3(lb,lb,lb), Bmax= M + R3(lb,lb,lb);
    double hseuil= Dmesh*1e-12;// ????
    typedef Mesh3::Element Element;
    typedef Element::Vertex Vertex;
    
    KNM<double> &pts=*ppoints;
    KN<long> &be=*pbe;
    KN<long> &loop = *ploop;
    
    long nt =Th.nt;
    long nv = Th.nv;
    int nbp2=0;
    for (int k=0; k<Th.nt; k++)
    {
        const Element & K = Th[k];
        int fi[4],kf[4],nfi=0;
        R3 PF[4];
        for(int i=0; i< 4;++i)
        {
            int i0 =Th(K[Element::nvface[i][0]]);
            int i1 =Th(K[Element::nvface[i][1]]);
            int i2 =Th(K[Element::nvface[i][2]]);
            Complex u2[]={u[i0],u[i1],u[i2]};
            R2 P2;
            if (in(u2,P2, eps))
            {nbp2++;}
        }
    }
    int nbpx = nbp2;
    KN<Vertex> Pf(nbpx);
    int nbp=0;
    EF23::GTree<Vertex> *gtree=new EF23::GTree<Vertex>(Pf,Bmin,Bmax,0);
    hseuil = 2./gtree->coef; //  seuil minimal dans gtree
    if(verbosity>9) cout << " hseuil minimal "<< hseuil << endl;
    long ihseuil=gtree->coef*hseuil;
    ffassert(ihseuil);
    int nbarc =0;
    typedef std::array<int,2> Arc;
    vector<Arc> arc;
    double hseuil2=hseuil*hseuil;
    double charge;
    const double twopi = 2.*Pi;
    for (int k=0; k<Th.nt; k++)
    {
        const Element & K = Th[k];
        int fi[4],kf[4],ip[4],nfi=0;
        R3 PF[4];
        /*
         // IN TEST
        int ku40 = 0, ku4[4];
        int j0 =Th(K[0]);
        int j1 =Th(K[1]);
        int j2 =Th(K[2]);
        int j3 =Th(K[3]);
        Complex u4[]={u[j0],u[j1],u[j2],u[j3]};
        
        if(abs(u4[0])==0) ku4[ku40++] = 0;
        if(abs(u4[1])==0) ku4[ku40++] = 1;
        if(abs(u4[2])==0) ku4[ku40++] = 2;
        if(abs(u4[3])==0) ku4[ku40++] = 3;
        if(ku40>2) cout<<" we have 2 zero face on the tetraedre "<<k<<" "<<ku40<<endl;
        */
        for(int i=0; i< 4;++i)
        {
            int i0 =Th(K[Element::nvface[i][0]]);
            int i1 =Th(K[Element::nvface[i][1]]);
            int i2 =Th(K[Element::nvface[i][2]]);
            Complex u2[]={u[i0],u[i1],u[i2]};
            /*
             // IN TEST
            int ku0 = 0, ku[3];
            if(abs(u2[0])==0) ku[ku0++] = 0;
            if(abs(u2[1])==0) ku[ku0++] = 1;
            if(abs(u2[2])==0) ku[ku0++] = 2;
            if(ku0>1) cout<<" we have a zero face "<<k<<" "<<i<<" "<<ku0<<endl;
            */
            R2 P2;
            if (in(u2,P2, eps))
            {
                Complex icharge = (log(u2[1]/u2[0]) + log(u2[2]/u2[1]) + log(u2[0]/u2[2]))/twopi;
                charge = imag(icharge);
                if(isnan(real(icharge)) || isnan(imag(icharge))){// We divide by zero
                    charge=1.;
                }
                if(abs(charge)>.98){
                    R3 P=K(K.PBord(i,P2));//
                    Vertex * pvi=gtree->ToClose(P,hseuil,true);
                    // verif brute force
                    if(!pvi && nbp)
                    {
                        pvi=gtree->NearestVertex(P,true);
                        int j = pvi - Pf;
                        double l2 = R3(*pvi,P).norme2();
                        if(l2 > hseuil2) pvi=0;
                        else    if( l2 < hseuil2) {
                            if(verbosity>9) cout << " bug in ToClose ??? " << k << " "<< i  << " == " << j << " :  " << P << " j " << Pf[j] << " / " << l2 <<endl;
                        }
                        
                    }
                    if(0) //  force brute ???
                        for(int j=0; j< nbp; ++j)
                    {
                        double l2 = R3(P,Pf[j]).norme2();
                        if( l2 < hseuil2) {
                            if(verbosity>9) cout << " bug " << k << " "<< i  << " == " << j << " :  " << P << " j " << Pf[j] << " / " << l2 <<endl;
                        }
                        
                    }
                    if(! pvi)
                    {
                        ffassert(nbp<nbpx);
                        Pf[nbp].x = P.x;
                        Pf[nbp].y = P.y;
                        Pf[nbp].z = P.z;
                        Pf[nbp].lab = 4*k+i;
                        gtree->Add(Pf[nbp]);
                        pvi = &Pf[nbp];
                        nbp++;
                    }
                    fi[nfi]=i;
                    ip[nfi]= pvi- (Vertex *) Pf;
                    nfi++;
                }
            }
            // compress ip
            sort(ip,ip+nfi);
            int nfi2=unique (ip, ip+nfi)-ip;
            if( nfi2 != nfi && verbosity>9) cout << " nfi "<< nfi << " " << nfi2 << endl;
            ffassert( nfi2 <=2);
            if(nfi2==2)
            {
                int i0= ip[0],i1=ip[1];
                Arc A={i0,i1};
                arc.push_back(A);
            }
        }
        // in ArcF(i0,i1) i \in [i0,i1], if i%5==4 => i/5 numéro thétraèdre else (i%5) numéro de face dans le tétraèdre i/5
    }
    delete gtree;
    sort(arc.begin(),arc.end());
    int nbua=unique(arc.begin(),arc.end())-arc.begin();
    if(verbosity>9) cout << " nbua "<< nbua << " " << arc.size() << " nbp " << nbp << " / " << nbp2 << endl;
    arc.resize(nbua);
    vector<vector<int>> av(nbp);
    if(verbosity>99) cout << " nb arc = "<< arc.size() << endl;
    for(int i=0;i<arc.size(); ++i)
    {
        if(verbosity>99) cout << " arc "<< i <<" " <<  arc[i][0] <<" " << arc[i][1] << endl;
        for(int j=0;j<2; ++j)
        {
            av[arc[i][j]].push_back(2*i+j);
        }
    }
    // recherche des banche du graphe
    vector<int> b; // debut de branche et fin de branche
    KN<long> next(nbp);
    next=-1L;
    int nca =0;
    KN<int> ca(arc.size(),0);
    auto which = [&arc] (int a,int v) {ffassert( arc[a][0] ==v ||  arc[a][1] ==v); return arc[a][0] ==v ? 0 : 1 ;};
    auto nexta = [&arc,&av] (int i,int a) {
        int aa=-1;
        if(av[i].size() ==2)
        {
            if(verbosity>99) cout << i << " " << a << " : " << av[i][0]/2 << " "<< av[i][1]/2 << endl;
            int k = av[i][1]/2 !=  a ;
            
            ffassert( av[i][1-k]/2 == a);
            aa = av[i][k];
        }
        return aa;};
    auto branch = [&ca , &arc, & next, &which, &nexta ,  &nca, & b, &av](int a,int i) {
        // depart a, sommet ai de a
        if( ca[a]) return ;
        int ia=which(a,i);
        int lg=1;
        // depart de branch ..
        b.push_back(2*a+ia);
        ca[a]=1;
        int v = arc[a][1-ia];
        if(verbosity>9) cout << " branch "<< a << " " << ia << " s= "<< i << " o= " << av[i].size() << " -> "<< v << " " ;
        while(1) {
            int aa = nexta(v,a);
            if(aa <0) { if(verbosity>9) cout << " fin branch s= " <<v << " o=  " << av[v].size() << " lg = " << lg << endl;
                break;
            }
            if( ca[aa/2]) { if(verbosity>9 ) cout << " fin loop "<< aa/2 << " s= " << v <<  " o=" << av[v].size()   << " lg = " << lg << endl;
                break;
            }
            a =aa/2;
            ia = aa%2;
            int vv = arc[a][1-ia];
            next[v] = vv;
            if(verbosity>99) cout<< " -> " << vv << endl;
            v = vv;
            ffassert(nca < arc.size())  ;
            nca++;
            lg++;
            ca[a]=1;
        }
        
        return;
    };
    
    for( int i=0; i<nbp; ++i)
    if( av[i].size()>0 && av[i].size()!=2 ) {
        // un depart possible
        if(verbosity>9) cout << " depart a "<< av[i].size() << " a= " <<av[i][0]/2<< " " << endl;
        for(int j=0; j<av[i].size();++j)
        branch(av[i][j]/2,i);
    }
    int nbline = b.size();
    loop.resize(nbline);
    for(int i=0; i< nbline; ++i)
    loop(i)=0;
    //  recheche des autre depart (les boucle)
    for( int aa=0; aa<arc.size(); ++aa)
    {
        if(ca[aa]==0) // no colorie , un depart possible dans le 2 sens
        {
            if(verbosity>9) cout << " depart boucle " << aa << arc[aa][0] << endl;
            branch(aa,arc[aa][0]);
        }
    }
    loop.resize(b.size());
    for(int i=nbline; i< b.size(); ++i)
    loop(i)=1;
    
    KNM<double> points(nbp+b.size(),3);
    
    int nbe=0;
    int nbc=0;
    {
        if(verbosity>99) cout << " parcours de branch .." << endl;
        for(int i=0; i< b.size(); ++i)
        {
            int a = b[i]/2, j =b[i]%2, s0=arc[a][j] , s= arc[a][1-j],ss;
            if(verbosity>99) cout << " branch "<< i << " s0 = "<< s0 << " " << av[s0].size() << ": "  << endl;
            nbc++;
            be.resize(nbc*2);
            be(2*(nbc-1))=nbe;
            R3 p0=Pf[s0];
            points(nbe,0)=p0.x; points(nbe,1)=p0.y; points(nbe,2)=p0.z; nbe++;
            while ((ss=next[s])>=0)
            {
                R3 pss=Pf[ss];
                points(nbe,0)=pss.x; points(nbe,1)=pss.y; points(nbe,2)=pss.z; nbe++;
                if(verbosity>99) cout << s << " -> ";
                s=ss;
            }
            be(2*(nbc-1)+1)=nbe-1;
            if(verbosity>99) cout  << " :  " << av[s].size() << "-1" << endl;
        }
    }
    points.resize(nbe,3);
    pts.resize(nbe,3);
    pts=points;
    return b.size();
}

// Computes the curvature of a curve
// Inputs:
//    pTH: meshL representing the curve
// Outputs:
//    pC: list storing the curvature at points in pTH
// Returns:
//    0
long curvatureL(pmeshL const &pTh, KN<double>*const &pc)
{
    typedef MeshL::Element Element;
    KN<double> &c = *pc;
    const MeshL &Th = *pTh;
    int nt=Th.nt,nv=Th.nv;
    ffassert(c.N()==nv);
    
    KN< int > cn(nv);
    KN< double > le(nv);
    c = 0.;
    cn = 0;
    le = 0.;
    
    for(int k=0;k<nt;k++){
        const Element & E(Th[k]);
        double lE = E.mesure();
        R3 EV(E[0],E[1]);// vector 0->1
        int ii=1, kk = Th.ElementAdj(k,ii);
        if(kk<0) continue;
        const Element & EAdj(Th[kk]);
        R3 EVAdj(EAdj[0],EAdj[1]);// vector 0->1
        double lAdj = EAdj.mesure();
        double scalprod = (EV,EVAdj);
        double cosa = scalprod/lE/lAdj;
        double aa = acos(cosa);
        c[Th(k,1)] = 2.* aa / (lE + lAdj);
    }
    return 0L;
}

// Computes the arc-length of a curve
// Inputs:
//    x: list of x coordinates for points in the curve
//    y: list of y coordinates for points in the curve
//    z: list of z coordinates for points in the curve
// Outputs:
//    ss: Arc-length at points in xyz
// Returns:
//    l: the length of the curve
double abscisses(KN_<double> const &  x,KN_<double> const &  y,KN_<double> const &  z,KN_<double> const&  ss)
{
    assert( x.N()==ss.N());
    assert( y.N()==ss.N());
    assert( z.N()==ss.N());
    KN_<double> s=ss;
    double l=0;
    s[0]=l;
    R3 P(x[0],y[0],z[0]);
    for(int i=1; i<ss.N();++i)
    {
        R3 Q(x[i],y[i],z[i]);
        l += R3(P,Q).norme();
        s[i]=l;
        P=Q;
    }
    return l;
}

// Use linear interpolation to send data from irregular discretization to regular discretization
// Inputs:
//    so: abscissa values that form an irregular discretization
//    xo: Values on the irregular discretization
// Outputs:
//    xn: Interpolated values on regular discretization
// Returns:
//    l: the length of the curve
double interpol(KN_<double> const &  so,KN_<double> const &  xo,KN_<double> const &  xn)
{
    int N = xn.N();
    int M = so.N();
    double l = so[M-1];
    double dl = l/(N-1.),si=0, si1;
    ffassert(so.N()==xo.N());
    int i0 = 0;
    for(int i=0; i<N;++i)
    {
        si = i*dl;
        // find i0  such that  [so[i],so[i+1] [
        while (i0+2<M)
        { // cout << i0<<" " << si << " " << so[i0+1] << " " << (si < so[i0+1] ) <<endl;
            if(si < so[i0+1]) break;
            else ++i0;
         }
        
        double si0 = so[i0];
        double si1 = so[i0+1];
        ffassert(si0 <= si && si <= si1);
        double l = (si-si0)/(si1-si0);
        xn[i] = xo[i0]*(1-l)+xo[i0+1]*(l);
    }
    return l;
}

static void inittt( ) {
    Global.Add("uZero2D", "(",new OneOperator5_<long,const Mesh *,KNM<double> *,KN<Complex> *,KN<double> *,double* >(uZero2D));
    Global.Add("uZero", "(",new OneOperator3_<long,pf3c,pf3r,double>(uZero));
    Global.Add("ZeroLines", "(",new OneOperator5_<long,pf3c,double,KNM<double> *,KN<long>*,KN<long>*> (ZeroLines) );
    Global.Add("BSp", "(",new OneOperator2s_<KNM_<double>,KNM_<double>,long >(BSp));
    Global.Add("curvatureL", "(", new OneOperator2_< long, pmeshL, KN<double> * >(curvatureL));
    Global.Add("smoothCurve", "(",new OneOperator2_<long,KNM<double>*,long  >(smoothCurve));
    Global.Add("zero3", "(",new OneOperator3s_<R3*,Complex,Complex,Complex>(zero3));
    Global.Add("interpol", "(",new OneOperator3_<double ,KN_<double> > (interpol) );
    Global.Add("abscisses", "(",new OneOperator4_<double ,KN_<double> > (abscisses) );
}

LOADFUNC(inittt);
