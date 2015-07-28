//  add some mixte finite element IN 3D ...
// RT1// F. Hecht  DEC  2014
// ------------------------------------------------------------
//   Pk = P1^3  + P1h (x,y,z)  ( dim Pk = 4*3 + 3) = 15
//   3 dof / fac + 
//   test 
/*

--  edp script associed: 
 LaplaceRT1.edp
 lame-TD-NSS.edp
 test-ElementMixte.edp
 
 */

//ff-c++-LIBRARY-dep:

#include "ff++.hpp"
#include "AddNewFE.h"

/*
 Per vedere se ci sono errori di compilazione e per poter usare elementi in un file edp con load:
 scrivo in terminale ff-c++ P012_3d_Modif.cpp
 (cio' crea P012_3d_Modif.o e P012_3d_Modif.dylib)
 */

namespace Fem2D {
    
    //------------------- New: --------------------//
    
    // Nedelec order 1 (degree 2) 3D
    // TypeOfFE_Edge1_3d derived from GTypeOfFE<> (which is defined in FESpacen.hpp)
    
    class TypeOfFE_Edge1_3d : public GTypeOfFE<Mesh3>
    {
    public:
        typedef Mesh3  Mesh;
        typedef Mesh3::Element  Element;
        typedef GFElement<Mesh3>  FElement;
        
        static int dfon[];
        static const int d=Mesh::Rd::d;
        static const GQuadratureFormular<R1> QFe; // quadrature formula on an edge
        static const GQuadratureFormular<R2> QFf; // quadrature formule on a face
        int edgeface[4][3]; // not used?!
        TypeOfFE_Edge1_3d(); // constructor
        void FB(const What_d whatd,const Mesh & Th,const Mesh3::Element & K,const Rd & P,RNMK_ & val) const;
        void set(const Mesh & Th,const Element & K,InterpolationMatrix<RdHat> & M,int ocoef,int odf,int *nump) const;
    };
    
    int TypeOfFE_Edge1_3d:: dfon[] = {0,2,2,0}; // 2 dofs on each edge, 2 dofs on each face
    
    // Quadrature formula on an edge, exact for degree 3 (ok: int_e (deg2*t *lambda))
    const GQuadratureFormular<R1> TypeOfFE_Edge1_3d:: QFe(-1+2*2,2,GaussLegendre(2),true);
    // argomenti: exact, num pti integrazione, pti integrazione, clean (~GQuadratureFormular() {if(clean) delete [] p;})
    // GaussLegendre definita in QuadratureFormular.cpp (nel commento dice che è su [0,1])
    
    // Quadrature formule on a face, exact for degree 2 (ok: int_f (deg2*t)), internal quadrature points
    static GQuadraturePoint<R2> P_QuadratureFormular_T_2_intp[3] = {
        GQuadraturePoint<R2>(1./3.,R2(1./6.,4./6.)) ,
        GQuadraturePoint<R2>(1./3.,R2(4./6.,1./6.)) ,
        GQuadraturePoint<R2>(1./3.,R2(1./6.,1./6.)) };
    // GQuadratureFormular<R2> const QuadratureFormular_T_2_int(2,3,P_QuadratureFormular_T_2_intp);
    const GQuadratureFormular<R2> TypeOfFE_Edge1_3d:: QFf(2,3,P_QuadratureFormular_T_2_intp);
    
    // In Mesh3dn.cpp:
    // static const int  nvedgeTet[6][2] = { {0,1},{0,2},{0,3},{1,2},{1,3},{2,3} };
    // static const int  nvfaceTet[4][3] = { {3,2,1},{0,2,3},{3,1,0},{0,1,2} };
    // In GenericMesh.hpp:
    // Vertex& at(int i) {return *vertices[i];}
    // penso che vertices venga riempito nell'ordine dato dal file di mesh per ogni tetraedro
    // In GenericMesh.hpp:
    // Rd Edge(int i) const {ASSERTION(i>=0 && i <ne);
    // return Rd(at(nvedge[i][0]),at(nvedge[i][1]));}
    // In GenericMesh.hpp:
    // bool   EdgeOrientation(int i) const
    // { return &at(nvedge[i][0]) < &at(nvedge[i][1]);}
    
    // Constructor
    // (usando parte di quello madre)
    // dfon, d, nsub ('per grafica'),
    // kPi = NbcoefforInterpolation e' num di alfa di (13.1) (chapter 13 ff++doc)
    //     = 3(=numComp)*QFe.n(num quad pts edge)*2*ne(2 fncts per edge) +
    //       3(=numComp)*QFf.n(num quad pts face)*2*nf(2 fncts per face)
    // npPi = NbPtforInterpolation = ne*QFe.n+nf*QFf.n
    // invariantinterpolationMatrix = false (i.e. it depends on the tetrahedron), discon=true
    TypeOfFE_Edge1_3d:: TypeOfFE_Edge1_3d(): GTypeOfFE<Mesh3>(TypeOfFE_Edge1_3d::dfon,d,1, Element::ne*2*3*QFe.n+Element::nf*2*3*QFf.n, Element::ne*QFe.n+Element::nf*QFf.n, false,true)
    {
        assert(QFe.n); // a cosa serve?
        assert(QFf.n);
        R3 Pt[] = {R3(0.,0.,0.),R3(1.,0.,0.),R3(0.,1.,0.),R3(0.,0.,1.)}; // 4 ref tetrahedron vertices
        
        {
            // We build the interpolation pts on the edges of the reference tetrahedron:
            int i;
            i=0;
            for(int e=0; e<Element::ne; ++e)
                for(int q=0; q<QFe.n; ++q,++i)
                {
                    double x=QFe[q].x;
                    this->PtInterpolation[i] = Pt[Element::nvedge[e][0]]*(1.-x)+Pt[Element::nvedge[e][1]]*(x);
                    // rispetto all'originale ho scambiato x e 1-x ! ?? (forse non e' importante l'ordine perche' ho gruppo di pti di quadratura simmetrico ?)
                }
            // We build the interpolation pts on the faces of the reference tetrahedron:
            // (the index i mustn't be reinitialised!)
            for(int f=0; f<Element::nf; ++f)
                for(int q=0; q<QFf.n; ++q,++i)
                {
                    double x=QFf[q].x;
                    double y=QFf[q].y;
                    this->PtInterpolation[i] = Pt[Element::nvface[f][0]]*(1.-x-y) + Pt[Element::nvface[f][1]]*x + Pt[Element::nvface[f][2]]*y;
                    // penso che non sia importante l'ordine perche' i tre punti di quadratura sono simmetrici
                }
        }
        {
            // We build the indices in (13.1) : edge basis fncts
            int i=0, p=0; // i is the k in (13.1) (chapter 13 ff++doc)
            int e; // we will need e also below, in the part referred to faces
            for(e=0; e<(Element::ne)*2; e++) // loop on the 12 edge basis fcts
            {
                if (e%2==1) {p = p-QFe.n;} // if I consider an 'even' basis fct, the quad pts are the ones of the previous basis fnct (they correspond to the same edge)
                for(int q=0; q<QFe.n; ++q,++p) // loop on the 2 edge quadrature pts
                    for (int c=0; c<3; c++,i++) // loop on the 3 components
                    {
                        this->pInterpolation[i]=p; // pk in (13.1)
                        this->cInterpolation[i]=c; // jk in (13.1)
                        this->dofInterpolation[i]=e; // ik in (13.1)
                        this->coefInterpolation[i]=0.; // alfak: we will fill them with 'set' (below) because they depend on the tetrahedron
                    }
            }
            // We build the indices in (13.1) : face basis fncts
            // (the indices i and p mustn't be reinitialised)
            for(int f=0; f<(Element::nf)*2; f++) // loop on the 8 face basis fcts
            {
                if (f%2==1) {p = p-QFf.n;} // if I consider an 'even' basis fct, the quad pts are the ones of the previous basis fnct (they correspond to the same face)
                for(int q=0; q<QFf.n; ++q,++p) // loop on the 3 face quadrature pts
                    for (int c=0; c<3; c++,i++) // loop on the 3 components
                    {
                        this->pInterpolation[i]=p; // pk in (13.1)
                        this->cInterpolation[i]=c; // jk in (13.1)
                        this->dofInterpolation[i]=e+f; // ik in (13.1)
                        this->coefInterpolation[i]=0.; // alphak: we will fill them with 'set' (below) because they depend on the tetrahedron
                    }
            }
        }
        // cout <<  " ++ TypeOfFE_Edge1_3d():"<< this->PtInterpolation << endl; // a cosa serve?
    }
    
    // For the coefficients of interpolation alphak in (13.1) (same 'for loops' as above) (vedi foglio)
    void TypeOfFE_Edge1_3d:: set(const Mesh & Th,const Element & K,InterpolationMatrix<RdHat> & M,int ocoef,int odf,int *nump) const
    {
        int i=ocoef, p=0; // Edge0_3d: parto da ocoef a riempire M perche' 'potrei avere un elemento finito dentro un altro elemento finito' (non capito) ??qui??
        
        // Edge basis fncts:
        int e; // we will need e also below, in the part referred to faces
        for(e=0; e<(Element::ne)*2; e++)
        {
            int ee=e/2;
            R3 E=K.Edge(ee); // the edge local number is given by the integer division between e and 2
            int eo = K.EdgeOrientation(ee);
            if(!eo) E=-E;
            if (e%2==1) {p = p-QFe.n;} // if I consider an 'even' basis fct, the quad pts are the ones of the previous basis fnct (they correspond to the same edge)
            for(int q=0; q<QFe.n; ++q,++p)
            {
                double ll=QFe[q].x; // value of lambda_0 or lambda_1 // why not (1-QFe[q].x) ??
                if( (e%2+eo) == 1 ) ll = 1-ll; // if exactly one between e%2 and eo is equal to 1 (so the sum is equal to 1), take the other lambda
                // i.e. if I'm considering the 2nd dof of the edge (or) if the edge is badly oriented, take the other lambda (but not if both)
                for(int c=0; c<3; c++,i++)
                {
                    M.coef[i] = E[c]*QFe[q].a*ll;
                    // QFe[q].a e' il peso su [0,1] del pto di integrazione q
                    // QFe[q].x e' la x su [0,1] del pto di integrazione q,
                    // quindi la sua prima lambda e' 1-QFe[q].x e la sua seconda lambda e' QFe[q].x
                }
            }
        }
        // Face basis fncts:
        // (the indices i and p mustn't be reinitialised)
        for(int f=0; f<(Element::nf)*2; f++) // loop on the 8 face basis fcts
        {
            int ff=f/2; // the face number
            int iff=f%2; // 1st or 2nd dof of the face
            const Element::Vertex * fV[3] = {& K.at(Element::nvface[ff][0]), & K.at(Element::nvface[ff][1]), & K.at(Element::nvface[ff][2])};
            // We 'order' the 3 vertices of the face according to their global numbers:
            // i0 will be the local number in the FACE of the vertex with the smallest global number
            // i1 will be the local number in the face of the vertex with the second smallest global number
            // i2 will be the local number in the face of the vertex with the largest global number
            int i0=0, i1=1, i2=2;
            if(fV[i0]>fV[i1]) Exchange(i0,i1);
            if(fV[i1]>fV[i2]) { Exchange(i1,i2);
                if(fV[i0]>fV[i1]) Exchange(i0,i1); }
            // now local numbers in the tetrahedron:
            i0 = Element::nvface[ff][i0], i1 = Element::nvface[ff][i1], i2 = Element::nvface[ff][i2];
            int ie0=i0, ie1 = iff==0? i1 : i2; // edge for the face dof (its endpoints local numbers)
            R3 E(K[ie0],K[ie1]);
            if (iff) {p = p-QFf.n;} // if I consider an 'even' basis fct, the quad pts are the ones of the previous basis fnct (they correspond to the same face)
            for(int q=0; q<QFf.n; ++q,++p) // loop on the 3 face quadrature pts
                for (int c=0; c<3; c++,i++) // loop on the 3 components
                {
                    M.coef[i] = E[c]*QFf[q].a;
                    // con questa formula di quadratura QFf i pesi sono tutti uguali,
                    // se fossero diversi bisognerebbe vedere se i punti di quadratura si incollano bene su una faccia in comune a due tetraedri ??
                }
        }
        ffassert(i==M.ncoef && M.np == p);
    }
    
    // In Mesh3dn.hpp:
    //void Gradlambda(R3 * GradL) const
    //{
    //    R3 V1(at(0),at(1));
    //    R3 V2(at(0),at(2));
    //    R3 V3(at(0),at(3));
    //    R det1=1./(6.*mesure());
    //    GradL[1]= (V2^V3)*det1;
    //    GradL[2]= (V3^V1)*det1;
    //    GradL[3]= (V1^V2)*det1;
    //    GradL[0]=-GradL[1]-GradL[2]-GradL[3];
    //}
    
    /*
     FB: We are on a tetrahedron K
     (the local numbering of its 4 vertices is given by how they are listed where K is described in the mesh file).
     
     AT FIRST we BUILD the 'pre-basis' functions OMEGA in the order given the GLOBAL numbering of the tetrahedron vertices,
     that is the 1st examined edge is from the vertex with the 1st smallest GLOBAL number
     to the one the 2nd smallest GLOBAL number, and so on (see nvedege),
     and the 1st examined face is the one opposite the vertex with the smallest GLOBAL number
     (and 2 edges of the face are chosen and oriented looking at its 3 global numbers), and so on.
     In this way:
     - a dof which is common to 2 adjacent tetrahedra is the same in the two tetraedra
     (since the orientation of each edge is chosen using the global numbers),
     - we can use the coefficients giving a dual basis that were calculated for the reference tetrahedron
     since the structure of orientation of the edges is the same (up to a rotation) as the one used in the reference tetrahedron.
     ([On the contrary, in the previous not working version we built the omegas in the order given the LOCAL numbering,
     changing the orientation of each edge using the GLOBAL numbers, but in that way the structure of orientation of the edges
     was not the same as the one used in the reference tetrahedron, and actually we couldn't use the coefficients calculated.
     Indeed, it didn't work for the tetrahedra whose vertices were not listed with increasing global number e.g. Th2.mesh,
     but it worked for tetrahedra whose vertices were listed with increasing global number e.g. Th1 and Th2b,
     since in this second case the orientation of each edge was actually left unchanged.])
     
     BUT THEN, when we build the DUAL basis functions PHI, we use the permutation p20 to go back to the FreeFem numbering of the dofs
     (which follows the LOCAL numbering).
     For instance the 1st examined edge can be the 3rd edge looking at the local numbering.
     */
    
    // val contains the values of the basis functions and of their derivatives at the point of K corresponding to the point P of the reference tetrahedron, by components
    void TypeOfFE_Edge1_3d:: FB(const What_d whatd,const Mesh & Th,const Mesh3::Element & K,const Rd &P,RNMK_ & val) const
    {
        assert(val.N()>=20); // 20 degrees of freedom, why >= ??
        assert(val.M()==3); // 3 components
        // -------------
        // perm: the permutation for which the 4 tetrahedron vertices are listed with increasing GLOBAL number
        // (i.e. perm[0] is the local number of the vertex with the smallest global number, ...
        //       perm[3] is the local number of the vertex with the biggest global number.)
        const Element::Vertex * tV[4] = {& K.at(0), & K.at(1), & K.at(2), & K.at(3)};
        int k0=0, k1=1, k2=2, k3=3;
        if(tV[k0]>tV[k1]) Exchange(k0,k1);
        if(tV[k1]>tV[k2]) Exchange(k1,k2);
        if(tV[k2]>tV[k3]) Exchange(k2,k3);
        if(tV[k0]>tV[k1]) Exchange(k0,k1);
        if(tV[k1]>tV[k2]) Exchange(k1,k2);
        if(tV[k0]>tV[k1]) Exchange(k0,k1);
        int perm[4] = {k0,k1,k2,k3};
        // -------------
        // We build mynvface to be used here instead of the FreeFem nvface,
        // (in order to exploit the results of perm and write a better code),
        // in mynvface in all the triplets the numbers are increasing
        static const int  mynvface[4][3] = { {1,2,3},{0,2,3},{0,1,3},{0,1,2} };
        // -------------
        // If [a,b] is the i-th edge (where a,b are its vertices local numbers), edgesMap[(a+1)*(b+1)] = i
        int edgesMap[13] = {-1,-1,0,1,2,-1,3,-1,4,-1,-1,-1,5}; // a map<int,int> would be more slow
        // -------------
        // the 4 barycentric coordinates for the reference tetrahedron evaluated at the point P
        // (they have the same value at the real tetrahedron's point corresponding to the reference tetrahedron's point P)
        R l[] = {1.-P.sum(),P.x,P.y,P.z};
        R3 D[4];
        K.Gradlambda(D); // (riempie un array di 4 R3)
        val = 0;
        // -----
        int p20[20]; // the permutation from the dofs numbering of my tetrahedron (numbered using the GLOBAL vertex numbers) to the dofs numbering of FreeFem !!!!!
        for(int i=0; i<6; ++i) // edges
        {
            // see below
            int ii0 = Element::nvedge[i][0], ii1 = Element::nvedge[i][1];
            int i0 = perm[ii0]; int i1 = perm[ii1];
            int iEdge = edgesMap[(i0+1)*(i1+1)]; // i of the edge [i0,i1]
            p20[i*2] = iEdge*2;
            p20[i*2+1] = iEdge*2+1;
        }
        for(int j=0; j<4; ++j) // faces
        {
            // using perm the 1st examined face is the one opposite the vertex with the smallest GLOBAL number, and so on
            // (see below)
            int jFace = perm[j];
            p20[12+j*2] = 12+jFace*2;
            p20[12+j*2+1] = 12+jFace*2+1;
        }
        // -----
        
        if (whatd & Fop_D0) // Fop_D0 definito in FESpacen.hpp ma non capito bene il <<
        {
            R3 X = K(P); //usato nella riga commentata più sotto
            // First, the functions omega (they don't constitute a dual basis! only a basis)
            R3 omega[20];
            // 12 edge functions:
            for(int i=0; i<6; ++i)
            {
                int ii0 = Element::nvedge[i][0], ii1 = Element::nvedge[i][1];
                int i0 = perm[ii0]; int i1 = perm[ii1];
                // ! :
                // using perm, [i0,i1] is already from the smallest global number to the greatest global number,
                // since nvedge always gives indices which read perm from left to right
                //cout << i << " oriented edge (" << K.at(i0) << "->" << K.at(i1) << ")" << endl;
                omega[i*2] = l[i0]*(l[i0]*D[i1]-l[i1]*D[i0]);
                omega[i*2+1] = l[i1]*(l[i0]*D[i1]-l[i1]*D[i0]);
            }
            // 8 face functions:
            for(int j=0; j<4; ++j)
            {
                // In (my)nvface the face opposite the vertex k is (my)nvface[k][],
                // using perm the 1st examined face is the one opposite the vertex with the smallest GLOBAL number, and so on,
                // and, since in mynvface the numbers always increase, [i0,i1,i2] are already ordered with increasing global number
                int ii0 = mynvface[j][0];  int ii1 = mynvface[j][1];  int ii2 = mynvface[j][2];
                int i0 = perm[ii0];  int i1 = perm[ii1];  int i2 = perm[ii2];
                omega[12+j*2] = l[i2]*(l[i0]*D[i1]-l[i1]*D[i0]);
                omega[12+j*2+1] = l[i1]*(l[i0]*D[i2]-l[i2]*D[i0]);
            }
            
            // Now, the functions phi that really constitute a dual basis
            R3 phi[20];
            phi[p20[0]] = +4*omega[0]-2*omega[1]-4*omega[16]+2*omega[17]-4*omega[18]+2*omega[19];
            phi[p20[1]] = -2*omega[0]+4*omega[1]-2*omega[16]-2*omega[17]-2*omega[18]-2*omega[19];
            phi[p20[2]] = +4*omega[2]-2*omega[3]-4*omega[14]+2*omega[15]+2*omega[18]-4*omega[19];
            phi[p20[3]] = -2*omega[2]+4*omega[3]-2*omega[14]-2*omega[15]-2*omega[18]-2*omega[19];
            phi[p20[4]] = +4*omega[4]-2*omega[5]+2*omega[14]-4*omega[15]+2*omega[16]-4*omega[17];
            phi[p20[5]] = -2*omega[4]+4*omega[5]-2*omega[14]-2*omega[15]-2*omega[16]-2*omega[17];
            phi[p20[6]] = +4*omega[6]-2*omega[7]-4*omega[12]+2*omega[13]+2*omega[18]-4*omega[19];
            phi[p20[7]] = -2*omega[6]+4*omega[7]-2*omega[12]-2*omega[13]+4*omega[18]-2*omega[19];
            phi[p20[8]] = +4*omega[8]-2*omega[9]+2*omega[12]-4*omega[13]+2*omega[16]-4*omega[17];
            phi[p20[9]] = -2*omega[8]+4*omega[9]-2*omega[12]-2*omega[13]+4*omega[16]-2*omega[17];
            phi[p20[10]] = +4*omega[10]-2*omega[11]+2*omega[12]-4*omega[13]+2*omega[14]-4*omega[15];
            phi[p20[11]] = -2*omega[10]+4*omega[11]+4*omega[12]-2*omega[13]+4*omega[14]-2*omega[15];
            phi[p20[12]] = +8*omega[12]-4*omega[13];
            phi[p20[13]] = -4*omega[12]+8*omega[13];
            phi[p20[14]] = +8*omega[14]-4*omega[15];
            phi[p20[15]] = -4*omega[14]+8*omega[15];
            phi[p20[16]] = +8*omega[16]-4*omega[17];
            phi[p20[17]] = -4*omega[16]+8*omega[17];
            phi[p20[18]] = +8*omega[18]-4*omega[19];
            phi[p20[19]] = -4*omega[18]+8*omega[19];
            
            for(int k=0; k<20; ++k)
            {
                val(k,0,op_id) = phi[k].x;
                val(k,1,op_id) = phi[k].y;
                val(k,2,op_id) = phi[k].z;
            }
        }
        
        if (whatd & Fop_D1) // Derivatives wrt x,y,z
        {
            R3 omegadx[20];
            R3 omegady[20];
            R3 omegadz[20];
            // 12 edge functions:
            for(int i=0; i<6; ++i)
            {
                int ii0 = Element::nvedge[i][0], ii1 = Element::nvedge[i][1];
                int i0 = perm[ii0]; int i1 = perm[ii1];
                // using perm, [i0,i1] is already from the smallest global number to the greatest global number,
                // since nvedge always gives indices which read perm from left to right
                if (whatd & Fop_dx)
                {
                    omegadx[i*2] = D[i0].x*(l[i0]*D[i1]-l[i1]*D[i0]) + l[i0]*(D[i0].x*D[i1]-D[i1].x*D[i0]);
                    omegadx[i*2+1] = D[i1].x*(l[i0]*D[i1]-l[i1]*D[i0]) + l[i1]*(D[i0].x*D[i1]-D[i1].x*D[i0]);
                }
                if (whatd & Fop_dy)
                {
                    omegady[i*2] = D[i0].y*(l[i0]*D[i1]-l[i1]*D[i0]) + l[i0]*(D[i0].y*D[i1]-D[i1].y*D[i0]);
                    omegady[i*2+1] = D[i1].y*(l[i0]*D[i1]-l[i1]*D[i0]) + l[i1]*(D[i0].y*D[i1]-D[i1].y*D[i0]);
                }
                if (whatd & Fop_dz)
                {
                    omegadz[i*2] = D[i0].z*(l[i0]*D[i1]-l[i1]*D[i0]) + l[i0]*(D[i0].z*D[i1]-D[i1].z*D[i0]);
                    omegadz[i*2+1] = D[i1].z*(l[i0]*D[i1]-l[i1]*D[i0]) + l[i1]*(D[i0].z*D[i1]-D[i1].z*D[i0]);
                }
            }
            // 8 face functions:
            for(int j=0; j<4; ++j)
            {
                // In (my)nvface the face opposite the vertex k is (my)nvface[k][],
                // using perm the 1st examined face is the one opposite the vertex with the smallest GLOBAL number, and so on,
                // and, since in mynvface the numbers always increase, [i0,i1,i2] are already ordered with increasing global number
                int ii0 = mynvface[j][0];  int ii1 = mynvface[j][1];  int ii2 = mynvface[j][2];
                int i0 = perm[ii0];  int i1 = perm[ii1];  int i2 = perm[ii2];
                if (whatd & Fop_dx)
                {
                    omegadx[12+j*2] = D[i2].x*(l[i0]*D[i1]-l[i1]*D[i0]) + l[i2]*(D[i0].x*D[i1]-D[i1].x*D[i0]);
                    omegadx[12+j*2+1] = D[i1].x*(l[i0]*D[i2]-l[i2]*D[i0]) + l[i1]*(D[i0].x*D[i2]-D[i2].x*D[i0]);
                }
                if (whatd & Fop_dy)
                {
                    omegady[12+j*2] = D[i2].y*(l[i0]*D[i1]-l[i1]*D[i0]) + l[i2]*(D[i0].y*D[i1]-D[i1].y*D[i0]);
                    omegady[12+j*2+1] = D[i1].y*(l[i0]*D[i2]-l[i2]*D[i0]) + l[i1]*(D[i0].y*D[i2]-D[i2].y*D[i0]);
                }
                if (whatd & Fop_dz)
                {
                    omegadz[12+j*2] = D[i2].z*(l[i0]*D[i1]-l[i1]*D[i0]) + l[i2]*(D[i0].z*D[i1]-D[i1].z*D[i0]);
                    omegadz[12+j*2+1] = D[i1].z*(l[i0]*D[i2]-l[i2]*D[i0]) + l[i1]*(D[i0].z*D[i2]-D[i2].z*D[i0]);
                }
            }
            // Now, the functions phi that really constitute a dual basis
            R3 phidx[20];
            if (whatd & Fop_dx)
            {
                phidx[p20[0]] = +4*omegadx[0]-2*omegadx[1]-4*omegadx[16]+2*omegadx[17]-4*omegadx[18]+2*omegadx[19];
                phidx[p20[1]] = -2*omegadx[0]+4*omegadx[1]-2*omegadx[16]-2*omegadx[17]-2*omegadx[18]-2*omegadx[19];
                phidx[p20[2]] = +4*omegadx[2]-2*omegadx[3]-4*omegadx[14]+2*omegadx[15]+2*omegadx[18]-4*omegadx[19];
                phidx[p20[3]] = -2*omegadx[2]+4*omegadx[3]-2*omegadx[14]-2*omegadx[15]-2*omegadx[18]-2*omegadx[19];
                phidx[p20[4]] = +4*omegadx[4]-2*omegadx[5]+2*omegadx[14]-4*omegadx[15]+2*omegadx[16]-4*omegadx[17];
                phidx[p20[5]] = -2*omegadx[4]+4*omegadx[5]-2*omegadx[14]-2*omegadx[15]-2*omegadx[16]-2*omegadx[17];
                phidx[p20[6]] = +4*omegadx[6]-2*omegadx[7]-4*omegadx[12]+2*omegadx[13]+2*omegadx[18]-4*omegadx[19];
                phidx[p20[7]] = -2*omegadx[6]+4*omegadx[7]-2*omegadx[12]-2*omegadx[13]+4*omegadx[18]-2*omegadx[19];
                phidx[p20[8]] = +4*omegadx[8]-2*omegadx[9]+2*omegadx[12]-4*omegadx[13]+2*omegadx[16]-4*omegadx[17];
                phidx[p20[9]] = -2*omegadx[8]+4*omegadx[9]-2*omegadx[12]-2*omegadx[13]+4*omegadx[16]-2*omegadx[17];
                phidx[p20[10]] = +4*omegadx[10]-2*omegadx[11]+2*omegadx[12]-4*omegadx[13]+2*omegadx[14]-4*omegadx[15];
                phidx[p20[11]] = -2*omegadx[10]+4*omegadx[11]+4*omegadx[12]-2*omegadx[13]+4*omegadx[14]-2*omegadx[15];
                phidx[p20[12]] = +8*omegadx[12]-4*omegadx[13];
                phidx[p20[13]] = -4*omegadx[12]+8*omegadx[13];
                phidx[p20[14]] = +8*omegadx[14]-4*omegadx[15];
                phidx[p20[15]] = -4*omegadx[14]+8*omegadx[15];
                phidx[p20[16]] = +8*omegadx[16]-4*omegadx[17];
                phidx[p20[17]] = -4*omegadx[16]+8*omegadx[17];
                phidx[p20[18]] = +8*omegadx[18]-4*omegadx[19];
                phidx[p20[19]] = -4*omegadx[18]+8*omegadx[19];
                
                for(int k=0; k<20; ++k)
                {
                    val(k,0,op_dx) = phidx[k].x;
                    val(k,1,op_dx) = phidx[k].y;
                    val(k,2,op_dx) = phidx[k].z;
                }
            }
            R3 phidy[20];
            if (whatd & Fop_dy)
            {
                phidy[p20[0]] = +4*omegady[0]-2*omegady[1]-4*omegady[16]+2*omegady[17]-4*omegady[18]+2*omegady[19];
                phidy[p20[1]] = -2*omegady[0]+4*omegady[1]-2*omegady[16]-2*omegady[17]-2*omegady[18]-2*omegady[19];
                phidy[p20[2]] = +4*omegady[2]-2*omegady[3]-4*omegady[14]+2*omegady[15]+2*omegady[18]-4*omegady[19];
                phidy[p20[3]] = -2*omegady[2]+4*omegady[3]-2*omegady[14]-2*omegady[15]-2*omegady[18]-2*omegady[19];
                phidy[p20[4]] = +4*omegady[4]-2*omegady[5]+2*omegady[14]-4*omegady[15]+2*omegady[16]-4*omegady[17];
                phidy[p20[5]] = -2*omegady[4]+4*omegady[5]-2*omegady[14]-2*omegady[15]-2*omegady[16]-2*omegady[17];
                phidy[p20[6]] = +4*omegady[6]-2*omegady[7]-4*omegady[12]+2*omegady[13]+2*omegady[18]-4*omegady[19];
                phidy[p20[7]] = -2*omegady[6]+4*omegady[7]-2*omegady[12]-2*omegady[13]+4*omegady[18]-2*omegady[19];
                phidy[p20[8]] = +4*omegady[8]-2*omegady[9]+2*omegady[12]-4*omegady[13]+2*omegady[16]-4*omegady[17];
                phidy[p20[9]] = -2*omegady[8]+4*omegady[9]-2*omegady[12]-2*omegady[13]+4*omegady[16]-2*omegady[17];
                phidy[p20[10]] = +4*omegady[10]-2*omegady[11]+2*omegady[12]-4*omegady[13]+2*omegady[14]-4*omegady[15];
                phidy[p20[11]] = -2*omegady[10]+4*omegady[11]+4*omegady[12]-2*omegady[13]+4*omegady[14]-2*omegady[15];
                phidy[p20[12]] = +8*omegady[12]-4*omegady[13];
                phidy[p20[13]] = -4*omegady[12]+8*omegady[13];
                phidy[p20[14]] = +8*omegady[14]-4*omegady[15];
                phidy[p20[15]] = -4*omegady[14]+8*omegady[15];
                phidy[p20[16]] = +8*omegady[16]-4*omegady[17];
                phidy[p20[17]] = -4*omegady[16]+8*omegady[17];
                phidy[p20[18]] = +8*omegady[18]-4*omegady[19];
                phidy[p20[19]] = -4*omegady[18]+8*omegady[19];
                
                for(int k=0; k<20; ++k)
                {
                    val(k,0,op_dy) = phidy[k].x;
                    val(k,1,op_dy) = phidy[k].y;
                    val(k,2,op_dy) = phidy[k].z;
                }
            }
            R3 phidz[20];
            if (whatd & Fop_dz)
            {
                phidz[p20[0]] = +4*omegadz[0]-2*omegadz[1]-4*omegadz[16]+2*omegadz[17]-4*omegadz[18]+2*omegadz[19];
                phidz[p20[1]] = -2*omegadz[0]+4*omegadz[1]-2*omegadz[16]-2*omegadz[17]-2*omegadz[18]-2*omegadz[19];
                phidz[p20[2]] = +4*omegadz[2]-2*omegadz[3]-4*omegadz[14]+2*omegadz[15]+2*omegadz[18]-4*omegadz[19];
                phidz[p20[3]] = -2*omegadz[2]+4*omegadz[3]-2*omegadz[14]-2*omegadz[15]-2*omegadz[18]-2*omegadz[19];
                phidz[p20[4]] = +4*omegadz[4]-2*omegadz[5]+2*omegadz[14]-4*omegadz[15]+2*omegadz[16]-4*omegadz[17];
                phidz[p20[5]] = -2*omegadz[4]+4*omegadz[5]-2*omegadz[14]-2*omegadz[15]-2*omegadz[16]-2*omegadz[17];
                phidz[p20[6]] = +4*omegadz[6]-2*omegadz[7]-4*omegadz[12]+2*omegadz[13]+2*omegadz[18]-4*omegadz[19];
                phidz[p20[7]] = -2*omegadz[6]+4*omegadz[7]-2*omegadz[12]-2*omegadz[13]+4*omegadz[18]-2*omegadz[19];
                phidz[p20[8]] = +4*omegadz[8]-2*omegadz[9]+2*omegadz[12]-4*omegadz[13]+2*omegadz[16]-4*omegadz[17];
                phidz[p20[9]] = -2*omegadz[8]+4*omegadz[9]-2*omegadz[12]-2*omegadz[13]+4*omegadz[16]-2*omegadz[17];
                phidz[p20[10]] = +4*omegadz[10]-2*omegadz[11]+2*omegadz[12]-4*omegadz[13]+2*omegadz[14]-4*omegadz[15];
                phidz[p20[11]] = -2*omegadz[10]+4*omegadz[11]+4*omegadz[12]-2*omegadz[13]+4*omegadz[14]-2*omegadz[15];
                phidz[p20[12]] = +8*omegadz[12]-4*omegadz[13];
                phidz[p20[13]] = -4*omegadz[12]+8*omegadz[13];
                phidz[p20[14]] = +8*omegadz[14]-4*omegadz[15];
                phidz[p20[15]] = -4*omegadz[14]+8*omegadz[15];
                phidz[p20[16]] = +8*omegadz[16]-4*omegadz[17];
                phidz[p20[17]] = -4*omegadz[16]+8*omegadz[17];
                phidz[p20[18]] = +8*omegadz[18]-4*omegadz[19];
                phidz[p20[19]] = -4*omegadz[18]+8*omegadz[19];
                
                for(int k=0; k<20; ++k)
                {
                    val(k,0,op_dz) = phidz[k].x;
                    val(k,1,op_dz) = phidz[k].y;
                    val(k,2,op_dz) = phidz[k].z;
                }
            }
        }
        
    }
    
    static TypeOfFE_Edge1_3d  Edge1_3d; // TypeOfFE_Edge1_3d e' il nome della classe appena definita
    GTypeOfFE<Mesh3> & GEdge13d(Edge1_3d); // GTypeOfFE<Mesh3> e' la classe madre
    static AddNewFE3 TypeOfFE_Edge1_3d("Edge13d",&GEdge13d); // Edge13d sara' il nome usato dall'utente
    
    //    // Come in BernardiRaugel.cpp (MA se scrivo in terminale ff-c++ P012_3d_Modif.cpp, ERRORE compilazione)
    //    //  ----   cooking to add the finite elemet to freefem table --------
    //    // a static variable to def the finite element
    //    static TypeOfFE_Edge1_3d Edge1_3d;
    //    //  now adding   FE in FreeFEm++  table
    //    static AddNewFE Edge13d("Edge13d",&Edge1_3d);
    //    // --- end cooking
    
    //    // Come a pag 345 docu (MA se scrivo in terminale ff-c++ P012_3d_Modif.cpp, ERRORI compilazione)
    //    static TypeOfFE_Edge1_3d Edge1_3d;
    //    static AddNewFE("Edge13d", Edge1_3d);
    //    static AddNewFE("Edge13d",&GEdge13d);
    
    
    class TypeOfFE_RT1_3d : public GTypeOfFE<Mesh3>
    {
    public:
        typedef Mesh3  Mesh;
        typedef Mesh3::Element  Element;
        typedef GFElement<Mesh3>  FElement;
        
        static int dfon[];
        static const int d=Mesh::Rd::d;
        static const GQuadratureFormular<R3> & QFt; // quadrature formula on an tet
        static const GQuadratureFormular<R2> QFf; // quadrature formule on a face
        TypeOfFE_RT1_3d(); // constructor
        void FB(const What_d whatd,const Mesh & Th,const Mesh3::Element & K,const Rd & P,RNMK_ & val) const;
       // void set(const Mesh & Th,const Element & K,InterpolationMatrix<RdHat> & M,int ocoef,int odf,int *nump) const;
    };
    int TypeOfFE_RT1_3d:: dfon[] = {0,0,3,3}; // 2 dofs on each edge, 2 dofs on each face
    
    // Quadrature formula on an edge, exact for degree 3 (ok: int_e (deg2*t *lambda))
    const GQuadratureFormular<R3> & TypeOfFE_RT1_3d:: QFt( QuadratureFormular_Tet_1);
    
     const GQuadratureFormular<R2> TypeOfFE_RT1_3d:: QFf(2,3,P_QuadratureFormular_T_2_intp);
    
     TypeOfFE_RT1_3d:: TypeOfFE_RT1_3d(): GTypeOfFE<Mesh3>(TypeOfFE_Edge1_3d::dfon,d,1,Element::nf*QFf.n+QFt.n, Element::nf*QFf.n, false,true)
    {
        ffassert(0);
    }
    void TypeOfFE_RT1_3d:: FB(const What_d whatd,const Mesh & Th,const Mesh3::Element & K,const Rd &P,RNMK_ & val) const
    {
        
        ffassert(0);
    }
    
} // chiude namespace Fem2D {



// --- fin -- 


