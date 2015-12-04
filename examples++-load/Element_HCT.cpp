//  Sur l'implementation des elements finis de Hsieh-Clough-Tocher complet et reduit
// http://opac.inria.fr/search~S12*frf?/aHastie%2C+Trevor+J./ahastie+trevor+j/-53%2C-1%2C0%2CB/frameset&FF=ahassan+kamal&5%2C%2C5
// F. Hecht and Mme. Hanen Narje : ferchichihanen@gmail.com
// -----------------------  related files:
//  to check  and validate  : testFEHCT.edp
//  to get a real example   :  bilapHCT.edp
// ------------------------------------------------------------

#include "ff++.hpp"
#include "AddNewFE.h"

// #include "problem.hpp"
namespace  Fem2D {
    
    // ------
    class TypeOfFE_HCT : public  TypeOfFE { public:
        static int Data[];
        // double Pi_h_coef[];
        
        TypeOfFE_HCT(): TypeOfFE(12 ,
                                 3,// hack   u, u_x, u_y for interpolation
                                 Data,
                                 2,
                                 1,
                                 9+6, // nb coef to build interpolation
                                 6, // np point to build interpolation
                                 0)
        {
            const double gauss1=(1.-sqrt(1./3.))/2;
            const double gauss2=1.-gauss1;
            
            const R2 Pt[] = { R2(0,0),R2(1,0),R2(0,1),R2(0.5,0.5), R2(0,0.5),R2(0.5,0)};
            // for the 3 vertices 3 coef => 9 coef ..
            int kk=0;
            for (int p=0;p<3;p++)
            {
                P_Pi_h[p]=Pt[p];
                pij_alpha[kk]= IPJ(kk,p,0);kk++;// VALUE
                pij_alpha[kk]= IPJ(kk,p,1);kk++;//DX
                pij_alpha[kk]= IPJ(kk,p,2);kk++;//DY
            }
            // for
            int p=3;
            for (int e=0;e<3;++e)
            { // point d'integration sur l'arete e
                
                P_Pi_h[p]= Pt[p];
                //	  cout <<"\n" <<  p << " --  " << P_Pi_h[p] << " ::  " << A << " " << B << endl;
                pij_alpha[kk++]= IPJ(9+e,p,1); // coef =  ne_x * sge
                pij_alpha[kk++]= IPJ(9+e,p,2); // coef =  ne_y * sge
                p++;
            }
            assert(P_Pi_h.N()==p);
            assert(pij_alpha.N()==kk);
            
        }
        void FB(const bool * whatd, const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
        void Pi_h_alpha(const baseFElement & K,KN_<double> & v) const;
        //     R operator()(const FElement & K,const  R2 & PHat,const KN_<R> & u,int componante,int op) const;
        
    } ;
    /*
     R TypeOfFE_HCT::operator()(const FElement & K,const  R2 & PHat,const KN_<R> & u,int componante,int op) const
     {
     R v[10000],vf[1000];
     assert(N*3*NbDoF<=10000 && NbDoF <1000 );
     KNMK_<R> fb(v,NbDoF,N,op+1); //  the value for basic fonction
     KN_<R> fk(vf,NbDoF);
     for (int i=0;i<NbDoF;i++) // get the local value
     fk[i] = u[K(i)];
     //  get value of basic function
     bool whatd[last_operatortype];
     for (int i=0;i<last_operatortype;i++)
     whatd[i]=false;
     whatd[op]=true;
     FB(whatd,K.Vh.Th,K.T,PHat,fb);
     cout << " N ===" <<N << " " <<NbDoF << " " <<op+1 << endl;
     cout << " fk = "<< fk << endl;
     cout << " bf = " <<componante << " " << op << "==="  << fb('.',componante,op)<< endl;
     R r = (fb('.',componante,op),fk);
     return r;
     }
     */
    //                     on what     nu df on node node of df
    int TypeOfFE_HCT::Data[]={
        0,0,0, 1,1,1,  2,2,2,  3,4,5, //  support on what
        0,1,2, 0,1,2,  0,1,2,  0,0,0, //df on node
        0,0,0, 1,1,1,  2,2,2,  3,4,5, // th node
        0,0,0, 0,0,0,  0,0,0,  0,0,0, // df of prevoius FE
        0,1,2,  3,4,5, 6,7,8, 9,10, 11, // ++ which df on prevoiu
        0,0,0,  // ???
        0,0,0,
        12,12,12
    };
    
    void TypeOfFE_HCT::Pi_h_alpha(const baseFElement & K,KN_<double> & v) const
    {
        const Triangle & T(K.T);
        int k=0;
        // coef pour les 3 sommets  fois le 3 composantes
        for (int i=0;i<9;i++)
            v[k++]=1;
        //   integration sur les aretes
        for (int i=0;i<3;i++)
        {
            
            R2 N(T.Edge(i).perp());
            N  *= T.EdgeOrientation(i)/N.norme();
            v[k++]= N.x;
            v[k++]= N.y;
        }
       // cout << " v =" << v << endl;
        ffassert(v.N()==k);
    }
    void set2zero(double *p,int k)
    {
        for(int i=0;i<k;++i)
            p[i]=0.;
    }
    
#define P3(a,b,c) a*b*c
    //#define P3aaa(a) a*a*a
    //#define P3abb(a,b) a*b*b
    //#define P3aaax(a) 3*a*a*a##x
    //#define P3aaaxy(a) 6*a*a##x*a##y
#define P3abcx(x,a,b,c) a##x*b*c+a*b##x*c+a*b*c##x
#define P3abcxy(x,y,a,b,c) a##x*b##y*c+a##y*b##x*c+a##y*b*c##x  +  a##x*b*c##y+a*b##x*c##y+a*b##y*c##x
#define P3abcxyz(x,y,z,a,b,c) a##x*b##y*c##z+a##y*b##x*c##z+a##y*b##z*c##x  +  a##x*b##z*c##y+a##z*b##x*c##y+a##z*b##y*c##x
    //#define P3abbx(x,a,b) a##x*b*b+2*a*b*b##x
    //#define P3abbxy(x,y,a,b) 2*a##x*b##y*b + 2*a##y*b*b##x + 2*a*b##y*b##x
    
#define P3X(a,b,c) P3abcx(x,a,b,c)
#define P3Y(a,b,c) P3abcx(y,a,b,c)
#define P3XY(a,b,c) P3abcxy(x,y,a,b,c)
#define P3XX(a,b,c) P3abcxy(x,x,a,b,c)
#define P3YY(a,b,c) P3abcxy(y,y,a,b,c)
#define P3XXX(a,b,c) P3abcxyz(x,x,x,a,b,c)
#define P3YYY(a,b,c) P3abcxyz(y,y,y,a,b,c)
#define P3XXY(a,b,c) P3abcxyz(x,x,y,a,b,c)
#define P3XYY(a,b,c) P3abcxyz(x,y,y,a,b,c)
    
#define LL10(P3) {P3(li,li,li),   P3(li1,li1,li1),   P3(li2,li2,li2), \
P3(li,li,li2),  P3(li,li,li1), \
P3(li1,li1,li), P3(li1,li1,li2), \
P3(li2,li2,li1),  P3(li2,li2,li), \
P3(li,li1,li2)} \

    void TypeOfFE_HCT::FB(const bool * whatd,const Mesh & ,const Triangle & K,const R2 & P,RNMK_ & val) const
    {
        typedef double R;
        double area= K.area;
        int Nop=val.K();
       // cout << "Nop = " << Nop << " " << whatd[op_id] <<  whatd[op_dx] <<  whatd[op_dy] << whatd[op_dxx] <<   whatd[op_dxy] << whatd[op_dyy]<< endl;
        R2 A(K[0]), B(K[1]),C(K[2]);
        R l[3] = {1-P.x-P.y,P.x,P.y};
        
        R2 Dl[3]=  {K.H(0), K.H(1), K.H(2) };
        
        
        R2 E[3]={ K.Edge(0),K.Edge(1),K.Edge(2)};
        R lg2[3] = { E[0].norme2(), E[1].norme2(),E[2].norme2()};
        R lg[3] = { sqrt(lg2[0]),sqrt(lg2[1]) ,sqrt(lg2[2])  };
        
        R eta[3]= { (lg2[2]-lg2[1])/lg2[0], (lg2[0]-lg2[2])/lg2[1], (lg2[1]-lg2[0])/lg2[2] };
        //    double l2E[3]={  (E[0],E[0]),  (E[1],E[1]),  (E[2],E[2]) };
        // double lE[3]={  sqrt(l2E[0]), sqrt(l2E[1]), sqrt(l2E[2]) };
        double sgE[3]={ K.EdgeOrientation(0), K.EdgeOrientation(1), K.EdgeOrientation(2)};
        // $ w_{3+i} = ccc[i] * ( li-2*li*li) $
        //   donc  $  D(w_i) =  ccc[i]  (1-2*li) Dl[i] $
        //  we must have $$ int_{e_i} dn(w_{3+j) ) =  \delta_{ij} $
        // $int_e_i dn(w_{3+i} )  = ccc[i] (Dl[i],Ne[i]) = 1 $
        //  $ ccc[i] = 1/  (Dl[i],Ne[i]) $
        R2 Ne[3]= {
            E[0].perp() *sgE[0],
            E[1].perp() *sgE[1],
            E[2].perp() *sgE[2]
        };
        val=0;
        
        throwassert( val.N()>=6);
        throwassert(val.M()==3);

        int i0 = 0;
        if( l[1] <l[i0]) i0=1;
        if( l[2] <l[i0]) i0=2;
        int i1 = (i0+1)%3, i2 = (i0+2)%3;
        // cout << " NMK= "<<val.N() << " " << val.M() << " " << val.K() << endl;
        // cout << " Ki :  <<"<< l[0] << " " << l[1] << " " << l[2] << " i0 = Ki == " <<  i0 << endl;
        double etai = eta[i0],etai1 = eta[i1],etai2 = eta[i2];
        double  li = l[i0], li1=l[i1], li2 =l[i2];
        double  lix = Dl[i0].x, li1x=Dl[i1].x, li2x = Dl[i2].x;
        double  liy = Dl[i0].y, li1y=Dl[i1].y, li2y = Dl[i2].y;
        //  i0,i1,i2,
        int p12[12]={3*i0,  3*i1, 3*i2,
            3*i0+2,3*i0+1,
            3*i1+2,3*i1+1,
            3*i2+2,3*i2+1,
            9+i0, 9+i1, 9+i2
            
        }; // renumerotation DL .. ff-> paper
        // int q12[12]; for(int i=0; i<12; ++i) q12[p12[i]]=i; // num  paper -> ff invers ...
        
        /*
         double ll[10]={ P3(li,li,li),   P3(li1,li1,li1),   P3(li2,li2,li2),
         P3(li,li,li2),  P3(li,li,li1),   P3(li1,li1,li),
         P3(li1,li1,li2), P3(li2,li2,li1),  P3(li2,li2,li),
         P3(li,li1,li2) };
         */
        // paper DOF fig 1.5.1 corresponding array Ai:
        // ( P(q_ij), Dp(q_ij-1 - q_ij), Dp(q_ij+1 - q_ij), (j=0,1,2) (9 DOF)
        //  Dp(b_ij)(h_ij)  ou bij = middle of edge ij (j=0,1,2) (3 DOF)
        //  Warning
        double ccc[] = { 1./(Dl[0],Ne[0]), 1./(Dl[1],Ne[1]), 1./(Dl[2],Ne[2]) };
        double c12 =1./12.;
        double Ai[12][10] = {
            {(-0.5)*(etai1-etai2)  ,  0,  0, (1.5)*(3+etai1),     (1.5)*(3-etai2),    0,     0,     0,     0,     0 },
            {(0.5)*(1-2*etai-etai2),  1,  0, (-1.5)*(1-etai),     (1.5)*(etai+etai2), 3,     3,     0,     0,     3*(1-etai)},
            {(0.5)*(1+2*etai+etai1),  0,  1, (-1.5)*(etai+etai1), (-1.5)*(1+etai),    0,     0,     3,     3,     3*(1+etai) },
            {(-c12)*(1+etai1)      ,  0,  0, (0.25)*(7+etai1),    (-0.5),             0,     0,     0,     0,     0 },
            {(-c12)*(1-etai2)      ,  0,  0, (-0.5),              (0.25)*(7-etai2),   0,     0,     0,     0,     0 },
            {(-c12)*(7+etai2)      ,  0,  0, (0.5),               (0.25)*(5+etai2),   1,     0,     0,     0,     -1},
            {(1./6.)*(4-etai)      ,  0,  0, (-0.25)*(3-etai),    (-0.25)*(5-etai),   0,     1,     0,     0,     (0.5)*(3-etai) },
            {(1./6.)*(4+etai)      ,  0,  0, (-0.25)*(5+etai),    (-0.25)*(3+etai),   0,     0,     1,     0,     (0.5)*(3+etai) },
            {(-c12)*(7-etai1)      ,  0,  0, (0.25)*(5-etai1),    (0.5),              0,     0,     0,     1,    -1},
            {(4./3.)               ,  0,  0, -2,                   -2,                0,     0,     0,     0,     4 },
            {(-2./3.)              ,  0,  0,  2,                   0,                 0,     0,     0,     0,     0 },
            {(-2./3.)              ,  0,  0,  0,                   2,                 0,     0,     0,     0,     0 }};
        //
        const int nnzdd= 6*3;// nb coef..
        double add[]={1., 1.,0.,0.,1.   ,1., 1.,0.,0.,1.,  1., 1.,0.,0.,1., 1.,1.,1.};
        int idd[]={0,1,1, 2,2, 3,4,4,5,5, 6,7,7,8,8,   9,10, 11};
        int jdd[]={0,1,2, 1,2, 3,4,5,4,5, 6,7,8,7,8,   9,10, 11};
        
        
        
        for(int i=0,kb=0,kc=15; i<3; ++i,kb += 5,++kc)
        {
            int ii=i*5+1;
            int ip=(i+1)%3;
            int is=(i+2)%3;
            R2 Es=  E[is];
            R2 Ep=  -E[ip];
            double dd[4]={Es.x,Es.y, Ep.x,Ep.y};//
            
            add[ii]=dd[0];
            add[ii+1]=dd[2];
            add[ii+2]=dd[1];
            add[ii+3]=dd[3];
            
            add[kc]= (sgE[i]*2*area/lg[i]);// hauteur
        }
        //       cout << endl;
        //       for (int k=0; k<nnzdd;++k)
        //           cout << idd[k] << " " << jdd[k] << " " << add[k] << endl;
        //       cout << endl;
        
        double i012[3]={1,0,0}, i112[3]={12,12,11};
        
        double AAA[12][10];
        set2zero(&AAA[0][0],120);
        double AA[12][10];
        set2zero(&AA[0][0],120);
        
        for (int jj=0;jj<10;++jj)
        {
            for (int i=0;i<12;++i)
            {
                AAA[p12[i]][jj]  += Ai[i][jj] ;
                AA[p12[i]][jj]  += Ai[i][jj] ;
            }
        }
        
        if(1)
        {
            set2zero(&AA[0][0],120);
            for(int k=0; k< nnzdd;++k)
            {
                int i= idd[k];
                int j= jdd[k];
                double dij=add[k];
                for(int jj=0; jj< 10; ++jj)
                    AA[i][jj] += dij*AAA[j][jj];
            }
        }
        
        
        
    
        
        
        if (whatd[op_id] || whatd[op_dx] ||  whatd[op_dy])
        {
          //  cout << "id dx dy"<< endl;

            double ll[10]= LL10(P3);
            double llx[10]= LL10(P3X);
            double lly[10]= LL10(P3Y);
            
            RN_ f(val('.',0,op_id));
            RN_ fx(val('.',1,op_id));
            RN_ fy(val('.',2,op_id));
             
            for(int i=0;i<12;++i)
                for(int j=0;j<10;++j)
                {
                    f[i] += AA[i][j]*ll[j];
                    fx[i] +=AA[i][j]*llx[j];
                    fy[i] +=AA[i][j]*lly[j];
                }
            if( whatd[op_dx]  ) val('.',0,op_dx)=fx;
            if( whatd[op_dy]  ) val('.',0,op_dy)=fy;
            
            
            //cout << "ll=" <<  KN_<double>(ll,10) << endl;
            
          //  cout << "f=" <<  f << endl;
        }
        
        if (whatd[op_dx] || whatd[op_dxx] ||  whatd[op_dxy])
        {
           // cout << "dx dxx dxy"<< endl;
            double ll[10]= LL10(P3X);
            double llx[10]= LL10(P3XX);
            double lly[10]= LL10(P3XY);
            
            RN_ f (val('.',0,op_dx));
            RN_ fx(val('.',1,op_dx));
            RN_ fy(val('.',2,op_dx));
            f=0.;
            for(int i=0;i<12;++i)
                for(int j=0;j<10;++j)
                {
                    f [i] += AA[i][j]*ll[j];
                    fx[i] +=AA[i][j]*llx[j];
                    fy[i] +=AA[i][j]*lly[j];
                }
            if( whatd[op_dxx]  ) val('.',0,op_dxx)=fx;
            if( whatd[op_dxy] )  val('.',0,op_dxy)=fy;
            
        }
        if (whatd[op_dy] || whatd[op_dyy] )
        {
            
          //  cout << "dy dyy "<< endl;
            double ll[10]= LL10(P3Y);
            double llx[10]= LL10(P3XY);
            double lly[10]= LL10(P3YY);
            
            RN_ f (val('.',0,op_dy));
            RN_ fx(val('.',1,op_dy));
            RN_ fy(val('.',2,op_dy));
            f=0.;
            for(int i=0;i<12;++i)
                for(int j=0;j<10;++j)
                {
                    f [i] +=AA[i][j]*ll[j];
                    fx[i] +=AA[i][j]*llx[j];
                    fy[i] +=AA[i][j]*lly[j];
                }
            if( whatd[op_dyy] ) val('.',0,op_dyy)=fy;
         }
      //  double NAN=nan("");
        if(Nop>op_dxx)
        {val('.',1,op_dxx)= NAN;
            val('.',2,op_dxx)= NAN;
        }
        if(Nop>op_dyy)
        {
            val('.',1,op_dyy)= NAN;
            val('.',2,op_dyy)= NAN;
        }
        if(Nop>op_dxy)
        {
            val('.',1,op_dxy)= NAN;
            val('.',2,op_dxy)= NAN;
        }
        
        
    }
    
    
    // a static variable to add the finite element to freefem++
    static TypeOfFE_HCT  Lagrange_HCT;
    static AddNewFE  FE_HCT("HCT",&Lagrange_HCT);
} // FEM2d namespace
// --- fin --


