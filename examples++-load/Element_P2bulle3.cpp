
/*
 P2 bulle 3 bulle per face + on e full bull
 
 */

//ff-c++-LIBRARY-dep:

#include "ff++.hpp"
#include "AddNewFE.h"
#include <iostream>



namespace Fem2D {
    
    // Author: F. Hecht , P-H Tournier, Jet Hoe Tang jethoe.tang@googlemail.com
    //   Jan 2017
    // in tets
    class TypeOfFE_P2_bulle3_3d : public GTypeOfFE<Mesh3>
    {
    public:
        typedef Mesh3  Mesh;
        typedef Mesh3::Element  Element;
        typedef GFElement<Mesh3>  FElement;
        
        static int dfon[];
        static const int d=Mesh::Rd::d;
        static const GQuadratureFormular<R1> QFe; // quadrature formula on an edge
        static const GQuadratureFormular<R2> QFf; // quadrature formula on a face
        TypeOfFE_P2_bulle3_3d(); // constructor
        void FB(const What_d whatd,const Mesh & Th,const Mesh3::Element & K,const Rd & P,RNMK_ & val) const;
        void set(const Mesh & Th,const Element & K,InterpolationMatrix<RdHat> & M,int ocoef,int odf,int *nump) const;
    };
    
    int TypeOfFE_P2_bulle3_3d:: dfon[] = {1,1,3,1}; // 2 dofs on each edge, 2 dofs on each face
    
    TypeOfFE_P2_bulle3_3d::TypeOfFE_P2_bulle3_3d():  GTypeOfFE<Mesh>(TypeOfFE_P2_bulle3_3d:: dfon,1,3,false,false)
    {
        typedef Element E;
        const  R alpha = 0.1885804846964451;// ; (7.-sqrt(13.))/18.;
       const   double a1= alpha, a0=(1-2*a1);
        int n=this->NbDoF;
        const int d= E::Rd::d;
        bool dd = verbosity>9;
        if(dd)
            cout << "\n +++ P2 3bulle : ndof : "<< n <<endl;
        R3 *Pt= this->PtInterpolation;
        // construction of interpolation ppoint

        {int k=0;
            for(int i=0;i<=d;++i)
                Pt[k++]=Rd();
            for(int i=0;i<d;++i)
                Pt[i+1][i]=1.;
            for(int i=0;i<E::ne;++i)
                Pt[k++] = (Pt[E::nvedge[i][0]]+Pt[E::nvedge[i][1]])*0.5;
            for(int i=0;i<E::nf;++i)
                for(int j=0; j<3;++j)
                {
                   
                    int j1=(j+1)%3,j2=(j+2)%3;
                    Pt[k++] = Pt[E::nvface[i][j]]*a0 +Pt[E::nvface[i][j1]]*a1 +Pt[E::nvface[i][j2]]*a1;
                    if(dd) cout << k-1<< " "  << Pt << endl;
                }
            Pt[k++]=Rd(0.25,0.25,0.25);
            if(dd)
            {
            cout << endl;
            cout << "k = " << k <<endl;
            cout << "alpha = " << alpha << endl;
            cout << "  interpolation nodes:" << endl;
            for (int i = 0 ; i < k; i ++ )            
            {                
                cout << "  "<<"["<< Pt[i].x<<" "<< Pt[i].y << " " << Pt[i].z << "]"<< endl;
            }
            cout << endl;
            }
            ffassert(k==n);
            cout << this->PtInterpolation<< endl;
        }
       
        for (int i=0;i<n;i++)
        {
            this->pInterpolation[i]=i;
            this->cInterpolation[i]=0;
            this->dofInterpolation[i]=i;
            this->coefInterpolation[i]=1.;
        }
    }
    void TypeOfFE_P2_bulle3_3d:: set(const Mesh & Th,const Element & K,InterpolationMatrix<RdHat> & M,int ocoef,int odf,int *nump) const
    {
        int n=this->NbDoF;
        int *p =M.p;
        for(int i=0; i< n;++i)
            M.p[i]=i;
        int k=10;
       if(verbosity>9) cout << " P2 3 bulle set:";
        for(int ff=0; ff<Element::nf; ff++,k+=3)
        {
            // oriantation de la face  a endroit
            int fp=K.facePermutation(ff);
             if(verbosity>9)  cout << fp << " ";
           // fp=0; // No perm
            if(fp&1) Exchange(p[k],p[k+1]);
            if(fp&2) Exchange(p[k+1],p[k+2]);
            if(fp&4) Exchange(p[k],p[k+1]);
        }
         if(verbosity>9)  cout << endl;
        /*
         for (int k=0;k<M.ncoef;++k)
         vdf[M.dofe[k]] += M.coef[k]*vpt(M.p[k],M.comp[k]);
         donc on a:
         vdf[i]= pt[p[i]]
         */
    }
    
    void TypeOfFE_P2_bulle3_3d:: FB(const What_d whatd,const Mesh & Th,const Mesh3::Element & K,const Rd &P,RNMK_ & val) const
    {
        assert(val.N()>=10+3*4+1); // 23 degrees of freedom
        assert(val.M()==1); // 3 components
        // -------------
        // perm: the permutation for which the 4 tetrahedron vertices are listed with increasing GLOBAL number
        // (i.e. perm[0] is the local number of the vertex with the smallest global number, ...
        //       perm[3] is the local number of the vertex with the biggest global number.)
        const Element::Vertex * tV[4] = {& K.at(0), & K.at(1), & K.at(2), & K.at(3)};
        static const int  nvf[4][3] ={{3,2,1}, {0,2,3},{ 3,1,0},{ 0,1,2}} ;
        static const int  nve[6][2] = { {0,1},{0,2},{0,3},{1,2},{1,3},{2,3} };
        // -------------
        // -------------
        // the 4 barycentric coordinates for the reference tetrahedron evaluated at the point P
        // (they have the same value at the real tetrahedron's point corresponding to the reference tetrahedron's point P)
      /*  R l[] = {1.-P.sum(),P.x,P.y,P.z};
        R lp[] = {l[1]*l[2]*l[3] , l[0]*l[2]*l[3] ,l[0]*l[1]*l[3] ,l[0]*l[1]*l[2]  };
        R lb = l[0]*l[1]*l[2]*l[3];
        R alpha = (7.-sqrt(13.))/18.;
        R alphap = 1-alpha;
        R aapap = alpha*alphap*alphap;*/
        R w = 1.-P.sum();
        R alpha = 0.1885804846964451; // alpha = (7-sqrt(13))/18;
        R c = 1./(alpha*alpha*(2.0*alpha-1.0)*(3.0*alpha-1.0));
        R l[] = {
            w*(2*w-1),
            P.x*(2*P.x-1),
            P.y*(2*P.y-1),
            P.z*(2*P.z-1),
            4*w*P.x,//
            4*w*P.y,
            4*w*P.z,
            4*P.x*P.y,
            4*P.x*P.z,
            4*P.y*P.z,
            
            (P.x*P.y*P.z*(P.z-alpha))*c,// 10 , v3
            (P.x*P.y*P.z*(P.y-alpha))*c,  // v2
            (P.x*P.y*P.z*(P.x-alpha))*c,   // v1
            
            (w*P.y*P.z*(  w-alpha))*c,//13,  v0
            (w*P.y*P.z*(P.y-alpha))*c,// 2
            (w*P.y*P.z*(P.z-alpha))*c,// 3
            
            (w*P.x*P.z*(P.z-alpha))*c,//16,   v3
            (w*P.x*P.z*(P.x-alpha))*c,// v1
            (w*P.x*P.z*(  w-alpha))*c,// v0
            
            (w*P.x*P.y*(  w-alpha))*c,//19,   0
            (w*P.x*P.y*(P.x-alpha))*c,// 1
            (w*P.x*P.y*(P.y-alpha))*c,// 2
            256*w*P.x*P.y*P.z
        };
        // 1 PB numerotation des arete
      
        int k=10;
        int p[23]  = {   0,1,2,3,4,5,6,7,8,9, 10,11,12, 13,14,15, 16,17,18, 19,20,21, 22};
        for(int ff=0; ff<Element::nf; ff++,k+=3)
        {
            // orientation de la face a envert
            int fp=K.facePermutation(ff);
            //fp=0; // No perm
            if(fp&4) Exchange(p[k],p[k+1]);
            if(fp&2) Exchange(p[k+1],p[k+2]);
            if(fp&1) Exchange(p[k],p[k+1]);
       }
        
        assert(val.N() >=E::nv+E::ne);
        assert(val.M()==1 );
        
        val=0;
        RN_ f0(val('.',0,op_id));
        //
        R beta1 = -0.1530178854880987;
        R beta2 =  0.1174552862797528;
        R beta3 = -0.4698211451190111;
        R beta4 = -0.1422503968333846;
        R beta5 =  0.1004882060140154;
        R beta6 = -0.0341147839626256;
        R beta7 = -0.0997720100233590;
        if (whatd & Fop_D0)
        {
            // corner
            f0[p[0]]=  l[0] + beta1*(l[13]+l[18]+l[19]) + beta2*(l[14]+l[15]+ l[17]+l[16] +l[20]+l[21]) + beta5*l[22];
            f0[p[1]]=  l[1] + beta1*(l[12]+l[17]+l[20]) + beta2*(l[10]+l[11]+ l[16]+l[18] +l[19]+l[21]) + beta5*l[22];
            f0[p[2]]=  l[2] + beta1*(l[11]+l[14]+l[21]) + beta2*(l[10]+l[12]+ l[13]+l[15] +l[19]+l[20]) + beta5*l[22];
            f0[p[3]]=  l[3] + beta1*(l[10]+l[15]+l[16]) + beta2*(l[11]+l[12]+ l[14]+l[13]+ l[17]+l[18]) + beta5*l[22];
            // edge
            f0[p[4]]=  l[4] + beta4*(l[16]+l[21] )+ beta3*(l[17]+l[18] +l[20]+l[19] )  + beta6*l[22];// 01 faces 2, 3
            f0[p[5]]=  l[5] + beta4*(l[15]+l[20] )+ beta3*(l[13]+l[14] +l[21]+l[19] )  + beta6*l[22];// 02 f 1 3
            f0[p[6]]=  l[6] + beta4*(l[14]+l[17] )+ beta3*(l[13]+l[15] +l[16]+l[18] )  + beta6*l[22];//  03 f 12
            f0[p[7]]=  l[7] + beta4*(l[10]+l[19] )+ beta3*(l[11]+l[12] +l[20]+l[21] )  + beta6*l[22];;//  12 f 0 3
            f0[p[8]]=  l[8] + beta4*(l[11]+l[18] )+ beta3*(l[10]+l[12] +l[16]+l[17] )  + beta6*l[22];//  13 // f 02
            f0[p[9]]=  l[9] + beta4*(l[12]+l[13] )+ beta3*(l[10]+l[11] +l[14]+l[15] )  + beta6*l[22];// 23  // 01
            // face
            f0[p[10]]=  l[10] + beta7*l[22];
            f0[p[11]]=  l[11] + beta7*l[22];
            f0[p[12]]=  l[12] + beta7*l[22];
            
            f0[p[13]]=  l[13] + beta7*l[22];
            f0[p[14]]=  l[14] + beta7*l[22];
            f0[p[15]]=  l[15] + beta7*l[22];
            
            f0[p[16]]=  l[16] + beta7*l[22];
            f0[p[17]]=  l[17] + beta7*l[22];
            f0[p[18]]=  l[18] + beta7*l[22];
            
            f0[p[19]]=  l[19] + beta7*l[22];
            f0[p[20]]=  l[20] + beta7*l[22];
            f0[p[21]]=  l[21] + beta7*l[22];
            // center
            f0[p[22]]=  l[22];
        }
        
//         cout << "BEGIN OF CPP  * * * * * * * * * * * * * * * * * * *  " << endl;
//         cout << "- f0:" << endl << f0  << endl;
//         cout << "- basis function:" << endl;
//         int ii;
//         for(ii = 0 ; ii < 23 ; ii++)
//         {
//             cout << l[ii] << ",";
//         }
//         cout << endl;
//         cout << "- permutation:" << endl;
//         for(ii = 0 ; ii < 23; ii++ )
//         {
//             cout << p[ii] << ",";
//         }
//         cout << endl;
//         cout << "- point:" << endl << "[" << P.x << "," << P.y << "," << P.z << "]"     << endl;
//         cout << "END OF CPP * * * * * * * * * * * * * * * * * * * * * " << endl;
    }
    
    
    
    
    static TypeOfFE_P2_bulle3_3d  P2Bulle3_3d;
    GTypeOfFE<Mesh3> & Elm_P2Bulle3_3d(P2Bulle3_3d);
    
    static AddNewFE3  TFE_P2Bulle3_3d("P2b33d",&Elm_P2Bulle3_3d);
    
    
} // closes namespace Fem2D {



// --- fin -- 


