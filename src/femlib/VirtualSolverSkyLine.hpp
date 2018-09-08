#include <iostream>
#include <cmath>
#include "HashMatrix.hpp"
#include <vector>
#include "VirtualSolver.hpp"

#include <complex>
#include "SkyLineSolver.hpp"
#include <queue>
#include <vector>
#include <utility>
//#include <climits>

template <class Z>
Z order_CutHill_McKee(Z n, Z *Ap, Z* Ai,Z*p)
{// CutHill-McKee
    std::vector<Z> mark(n,-1),r(n);
    Z nr =0;
    unsigned long long pft=0,pft0=0;
    std::multimap<Z,Z> q,qn;
    bool debug =0; 
    unsigned long long  pfso=0 ,pfs=0,pfss=0,pfs0=0;
    for(Z id=0; id<n; ++id)
        if(mark[id]==-1)
        { // start point
            Z nri=nr;
            Z s=id;
            
            for(int step=0; step < 5; ++step)
            {
                if(step) s = r[nr-1]; // last point ..
                 for (int i= nri; i< nr;++i)
                    mark[r[i]]=-1; // clean
                nr=nri; //
                mark[s]=nr;
                r[nr++] = s;
                Z sc = Ap[s+1]-Ap[s];
                q.clear();
                q.insert(std::make_pair(-sc,s));
                if(debug)  cout << " s= "<< s << " "<< nr << endl;
                int lvl=0;
                while (q.size()>0)
                {
                    if(debug)  cout << lvl++ <<" " << nr <<  " : ";
                    
                    qn.clear();
                    for(auto iq=q.begin(); iq != q.end(); ++iq )
                    {
                        int s = iq->second;
                        for(int k=  Ap[s]; k<Ap[s+1];++k)
                        {
                            Z ss= Ai[k];
                            if(mark[ss]<0)
                            {
                                Z ssc = Ap[ss+1]-Ap[ss];
                                qn.insert(std::make_pair(-ssc,ss));
                                mark[ss]=nr;
                                r[nr++] = ss;
                                if(debug)   cout << ss <<" ";
                            }
                            
                        }
                    }
                    if(debug)   cout << endl;
                    swap(q,qn);
                }
                
                //  compute SkyLine
                
                pfs =0;
                pfs0= 0;
                for(int i= nri; i<nr; ++i)
                {
                    Z s = r[i];
                    Z imn = s,i0mn=i;
                    Z imx = s,i0mx=i;
                    for(int k=  Ap[s]; k<Ap[s+1];++k)
                    {
                        Z ss0=Ai[k];
                        Z ss= mark[ss0];
                        imn = min(imn,ss);
                        imx = max(imx,ss);
                        i0mn = min(i0mn,ss0);
                        i0mx = max(i0mx,ss0);
                    }
                    pfs += s-imn;
                    pfs0 += i-i0mn;
                }
                if(debug)  cout << "         ** " << step << " " << pfso << " " << pfs << " " << pfs0 << endl;
                if(!pfso) pfso= pfs0;
                if(pfs <pfso )
                {
                    for(int i= nri; i<nr; ++i)
                      p[r[i]]=mark[r[i]];
                    cout << step<< "            pfs  " << pfs << " " << pfso << endl;
                    pfss=pfs;
                }
                else if(step==0)
                {
                  if(debug)  cout << step<< "       ID pfs  " << pfs << " " << pfs0 << endl;
                    for(int i= nri; i<nr; ++i)
                      p[r[i]]=r[i];
                }
                if((pfso <pfs) && (step > 5) ) break;
                pfso=min(pfs,pfso);
            }
          
            
          pft0 += pfs0;
          pft += pfss;
          if(verbosity>1) cout << "order_CutHill_McKee: pft = " << pft << " pft id = " << pft0 << endl;
    }
    if(verbosity>9)
    {
        for(int i=0; i< n; ++i)
        {  if(i%10==0) cout << endl;
            cout << p[i] << " ";
        }
        cout << endl<< endl;
    }
    return (Z) pft;
}

    template<class Z=int,class K=double>
    class VirtualSolverSkyLine: public VirtualSolver<Z,K> {
    public:
        typedef HashMatrix<Z,K>  HMat;
        typedef SkyLineMatrix<Z,K>  SLMat;
        int typesolver;
        HMat *A;
        SLMat *SL;
        int  cs,cn;
        Z *p; // permuation old -> new
        K *v;
        double tol_pivot;
        
        mutable int status;
        
        VirtualSolverSkyLine(HMat  *AA,int ttypesolver,double tol_pivott=1e-15)
        :typesolver(ttypesolver),A(AA),SL(0),cs(0),cn(0),p(0),v(0),
        tol_pivot(tol_pivott<0 ?  1e-15 :tol_pivott)  {}
        
        
        void SetState(){
            if( A->re_do_numerics ) cn++;
            if( A->re_do_symbolic) cs++;
            this->CheckState(A->n,cs,cn);
        }
        void dosolver(K *x,K*b,int N,int trans) {
            
            for(int k=0,oo=0; k<N;++k, oo+= A->n)
            {
                K * xx = x+oo,  *bb =  b+oo;
                for(Z i=0; i<A->n; ++i)
                    v[i] = bb[p[i]];
                SL->solve(v,trans);
                for(Z i=0; i<A->n; ++i)
                    xx[p[i]] = v[i];
                
            }
        }
        void fac_init(){
            
            if(p) delete p;
            p = new Z[A->n];
            if(v) delete v;
            v = new K[A->n];
        }
        void fac_symbolic(){
            // build permuation
            
            if( A->n !=A->m) MATERROR(1,"VirtualSolverSkyLine non squre mat");
            Z *Ap,*Ai;
            K * Ax;
            
            A->CSC(Ap,Ai,Ax);
            
            order_CutHill_McKee(A->n,Ap,Ai,p);
            //  rebuild perumation
            
        }
        void fac_numeric(){
            if(SL) delete SL;
            SL = new SLMat(A,p,typesolver);
            if(verbosity>1) std::cout << " size of Skyline mat ="<< SL->size() << " nz :" <<SL->pL[A->n] << " nzz " << A->nnz << endl;
            if(verbosity>99)
            {
                 cout << *SL << endl;
                cout << *A << endl;
            }
            SL->factorize(tol_pivot);
            
        }
        ~VirtualSolverSkyLine(){
            if(p) delete [] p;
            if(v) delete [] v;
            if(SL) delete SL;
        }
        
    };
    
    
    
    
