#include "ff++.hpp"
#include "AddNewFE.h"
// Fini element on Quadrature Formula ...
//  to opimize computation ...
/// ---------------------------------------------------------------
namespace  Fem2D {

    
    // -------------------
    class TypeOfFE_QF2d : public  TypeOfFE { public:
        static int *DataQF2d(int n)
        {
            int *d = new int[n*5+3];
            int *p=d;
            for(int i=0;i<n;++i)
                *p++=6;
            for(int i=0;i<n;++i)
                *p++=i;
            for(int i=0;i<n;++i)
                *p++=0;
            for(int i=0;i<n;++i)
                *p++=0;
            for(int i=0;i<n;++i)
                *p++=i;
            *p++=0;
            *p++=0;
            *p++=n;
            return d;
        }
     //   static double *cfill(int n) {  double *d = new double[n]]; for(int i=0; i<n; ++i) d[i]1.; return d;}
        typedef GQuadratureFormular<R2> QF;
        int *Data;
        
            int m;
            KN<int> w;
            
        TypeOfFE_QF2d(const QF *qf): TypeOfFE(0,0,qf->n,1,DataQF2d(qf->n),1,1,qf->n,qf->n,new double[qf->n] ),m(qf->n <15 ?5:7), w(m*m)
	{
            int debug=verbosity>99;
	    for (int i=0;i<NbDoF;i++) {
		pij_alpha[i]= IPJ(i,i,0);
		P_Pi_h[i]=(*qf)[i];
                coef_Pi_h_alpha[i]=1.;
	    }
            double dx = 1./m, dy=1./m;
           if(debug)  cout << " n " << NbDoF << endl;
            for( int i=0; i< m;++i)
              for( int j=0; j< m;++j)
              {
                  R2 P(dx*i+dx/2,dy*j+dy/2);
                  int k= i*m+j;
                  double dd = 100000;
                  
                  for (int l=0; l<NbDoF;++l)
                  {
                      R2 PQ(P,P_Pi_h[l]);
                      double dPQ = (PQ,PQ);
                      if( dd > dPQ) { dd=dPQ; w[k]= l;}
                  }
                 if(debug) cout << k << " " << i << " " << j << " :  " << w[k] << " ( " << P_Pi_h[w[k]] << " )" << dd << " || " << P <<  endl;
              }
            // verif
            int err=0;
             for (int l=0; l<NbDoF;++l)
                    {
                        int k = ijP(P_Pi_h[l]);
                        if(k !=l)
                        {
                            err++;
                            cout << " Erreur search TypeOfFE_QF2d loose point " << l << " / " << NbDoF << " m = " << m << endl;
                        }
                    }
            if(err)
            {
                ErrorExec("TypeOfFE_QF2d: increase m in TypeOfFE_QF2d ",0);
                ffassert(0);
            }
            

	}
        int ijP(const R2 &P) const {
            int i=  P.x*m ;
            int j= P.y*m ;
            i= min(i,m-1);
            j= min(j,m-1);
            return  w(i*m+j);
        }
	void FB(const bool * whatd,const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const
            {
                int k=ijP(P);
                val=0;
                if (whatd[op_id]) {
                 assert(val.K()>op_id);
                 val(k,0,op_id) =  1;
                }
            }
            
       /*     R operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const
            {
                
            }
        */
   
} ;


    class TypeOfFE_QF3d : public GTypeOfFE<Mesh3>
    {
    public:
        typedef Mesh3  Mesh;
        typedef Mesh3::Element  Element;
        typedef GFElement<Mesh3>  FElement;
        int np,*dfon,m;
        KN<int> w;
        static const int d=Mesh::Rd::d;
        typedef GQuadratureFormular<R1> QFe; // quadrature formula on an edge
        typedef GQuadratureFormular<R2> QFf; // quadrature formule on a face
        typedef GQuadratureFormular<R3> QFk; // quadrature formule on a element
        TypeOfFE_QF3d(const TypeOfFE_QF3d::QFk &qf  ); // constructor
        void FB(const What_d whatd,const Mesh & Th,const Mesh3::Element & K,const Rd & P,RNMK_ & val) const;
    //    void set(const Mesh & Th,const Element & K,InterpolationMatrix<RdHat> & M,int ocoef,int odf,int *nump) const;
        static int * pdfon(int n);
        ~TypeOfFE_QF3d(){ delete [] dfon;}
        
        int ijP(const R3 &P) const {
            int i=  P.x*m ;
            int j= P.y*m ;
            int k= P.z*m ;
            i= min(i,m-1);
            j= min(j,m-1);
            k= min(k,m-1);
            return  w(i*m*m+j*m+k);
        }

    };
    
    int * TypeOfFE_QF3d::pdfon(int n) { int *p=new int[4]; p[0]=p[1]=p[2]=0; p[3]=n; return p;}
    
    TypeOfFE_QF3d::TypeOfFE_QF3d(const TypeOfFE_QF3d::QFk &qf):
      GTypeOfFE<Mesh3>(dfon=pdfon(qf.n),1,1, qf.n , qf.n , true,true),
      np(qf.n),m(5),w(m*m*m)
    {
        for (int i=0;i<np;++i) {
            this->PtInterpolation[i]=qf[i];
            this->pInterpolation[i]=i;
            this->cInterpolation[i]=0;
            this->dofInterpolation[i]=i;
            this->coefInterpolation[i]=1;
        }
        //
        int debug = verbosity>99;
        if(debug)  cout << " n " << NbDoF << endl;
        double dd= 1./m;
        for( int i=0; i< m;++i)
            for( int j=0; j< m;++j)
                for( int k=0; k< m;++k)
            {
                R3 P(dd*i+dd/2,dd*j+dd/2,dd*k+dd/2);
                int kk= i*m*m+j*m +k;
                double dd = 1000000;
                
                for (int l=0; l<NbDoF;++l)
                {
                    R3 PQ(P,this->PtInterpolation[l]);
                    double dPQ = (PQ,PQ);
                    if( dd > dPQ) { dd=dPQ; w[kk]= l;}
                }
                if(debug) cout << kk << " " << i << " " << j << " " << k << " :  " << w[k] << " ( " << this->PtInterpolation[w[kk]] << " )" << dd << " || " << P <<  endl;
            }
        // verif
        int err=0;
        for (int l=0; l<NbDoF;++l)
        {
            int kk = ijP(this->PtInterpolation[l]);
            if(kk !=l)
            {
                err++;
                cout << " Erreur search TypeOfFE_QF3d loose point " << l <<" NbDoF=" << NbDoF <<  " m = " << m << endl;
            }
        }
        if(err)
        {
            
            ErrorExec("TypeOfFE_QF2d: increase m in TypeOfFE_QF3d ",0);
            ffassert(0);
        }
        
        
    }
    void  TypeOfFE_QF3d::FB(const What_d whatd,const Mesh & Th,const Mesh3::Element & K,const Rd &P, RNMK_ & val) const
    {
        assert(val.N() >=4);
        assert(val.M()==1 );
        int kqf=ijP(P);
        // wi = signe * (x - qi)/ (volume*d)
        val=0;
        if (whatd & Fop_D0)
             val(kqf,0,op_id) = 1 ;

    }

    
 

// link with FreeFem++ 
static TypeOfFE_QF2d TypeOfFE_QF2d1(&QuadratureFormular_T_1);
static TypeOfFE_QF2d TypeOfFE_QF2d2(&QuadratureFormular_T_2);
static TypeOfFE_QF2d TypeOfFE_QF2d5(&QuadratureFormular_T_5);
static TypeOfFE_QF2d TypeOfFE_QF2d7(&QuadratureFormular_T_7);
static TypeOfFE_QF2d TypeOfFE_QF2d9(&QuadratureFormular_T_9);

    static TypeOfFE_QF3d TypeOfFE_QF3d1(QuadratureFormular_Tet_1);
    static TypeOfFE_QF3d TypeOfFE_QF3d2(QuadratureFormular_Tet_2);
    static TypeOfFE_QF3d TypeOfFE_QF3d5(QuadratureFormular_Tet_5);
    
static AddNewFE  AddNewFE_QF2d1("FEQF1",&TypeOfFE_QF2d1);
static AddNewFE  AddNewFE_QF2d2("FEQF2",&TypeOfFE_QF2d2);
static AddNewFE  AddNewFE_QF2d5("FEQF5",&TypeOfFE_QF2d5);
static AddNewFE  AddNewFE_QF2d7("FEQF7",&TypeOfFE_QF2d7);
static AddNewFE  AddNewFE_QF2d9("FEQF9",&TypeOfFE_QF2d9);
static AddNewFE  AddNewFE_QF2ddef("FEQF",&TypeOfFE_QF2d5);
    
static AddNewFE3  AddNewFE_QF3d1("FEQF13d",&TypeOfFE_QF3d1);
static AddNewFE3  AddNewFE_QF3d2("FEQF23d",&TypeOfFE_QF3d2);
static AddNewFE3  AddNewFE_QF3d5("FEQF53d",&TypeOfFE_QF3d5);
static AddNewFE3  AddNewFE_QF3ddef("FEQF3d",&TypeOfFE_QF3d5);
    
} // FEM2d namespace 


// --- fin -- 

