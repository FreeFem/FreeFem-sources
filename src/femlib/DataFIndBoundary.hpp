#ifndef _DataFindBoundary_hpp__
#define _DataFindBoundary_hpp__

#include <RNM.hpp>
#include <fstream>
// for M_PI 
#ifdef __STRICT_ANSI__
#undef __STRICT_ANSI__
#endif

#include <cmath>

template<typename Mesh>
struct GenericDataFindBoundary
{
    typedef typename Mesh::Vertex Vertex;
    
    typedef typename Mesh::Element Element;
    typedef typename Mesh::BorderElement BorderElement;
    typedef  typename Mesh::Rd Rd;
    typedef  typename Mesh::RdHat RdHat;
    //using   EF23::GTree;

    static const int d = Rd::d;
    static const int dHat = RdHat::d;
    static const bool bborder = d == dHat; //  build border .??
    const Mesh *pTh;
    EF23::GTree<Vertex> *tree;
    KN<Vertex> P;// the barycentre of element
    KN<double> delta;// de dist of recheche
    KN<Vertex*> lp;// buffer of recheche ..
    long debug ;
    GenericDataFindBoundary(Mesh const * pTh,int ddebug=0);
    ~GenericDataFindBoundary() ;//{delete tree;}
    int Find(Rd P,double *l,int & outside) const ;
    void gnuplot(const string & fn);
};

template<typename Mesh>
GenericDataFindBoundary<Mesh>::~GenericDataFindBoundary()
{
    delete tree;
}
template<typename Mesh>
void GenericDataFindBoundary<Mesh>::gnuplot(const string & fn)
{ // for debugging ..
    if(dHat < 3)
    {
    ofstream gp(fn.c_str());
    ffassert(gp);
    //
    const Mesh &Th = *pTh;
    for(int be=0; be<Th.nbe; ++be)
    {
        const BorderElement &B=Th.be(be);
        int e,k = Th.BoundaryElement(be,e);
        {
            int ee=e, kk=  Th.ElementAdj(k,ee);
            if ( kk >=0 || k != kk)
            {
                for(int j=0; j< BorderElement::nv;++j)
                 gp  << (Rd) B[j] << endl;
               gp  << "\n\n";
            }
        }
    }
    if( dHat==2)
    {
    for(int i=0; i<P.N(); ++i)
    {
        int N=100;
        double dt = M_PI*2./N, r = delta[i];
        for(int j=0;j<=N; ++j)
        {
            double x = P[i].x+r*cos(dt*j);
            double y=  P[i].y+r*sin(dt*j);
            gp << x << " " << y << endl;
        }
        gp << "\n\n";
    }}
    }
}
template<typename Mesh>
int GenericDataFindBoundary<Mesh>::Find(typename Mesh::Rd PP,double *l,int & outside) const
{  // FH: outside : 0 inside, 1 out close, 2, out fare, , -1 inside
    typedef double R;
    int nu=-1,ne=-1;
    R dnu= 1e200;
    R dl[dHat+1];
    outside = 0;
    Vertex *p =tree->TrueNearestVertex(PP);
    int i = p-P;
    long lvp=tree->ListNearestVertex(lp,lp.N(), delta[i],P[i]);
    for(int j=0; j<lvp; ++j)
    {
        int k = lp[j]->lab/3;
        int e = lp[j]->lab%3;
        if(debug) cout << "    -- k = "<< k << " " << e << " " << j << endl;
        
        const Element & K=(*pTh)[k];
        int nl[3],n=0;
        Rd & A(K[0]), & B(K[1]), & C(K[2]);
        R l[3]={0,0,0};
        R area2= K.area*2;
        R eps =  -area2*1e-6;
        l[0] = Area2(PP,B,C);
        l[1] = Area2(A,PP,C);
        l[2] = area2-l[0]-l[1];
        if (l[0] < eps) nl[n++]=0;
        if (l[1] < eps) nl[n++]=1;
        if (l[2] < eps) nl[n++]=2;
        if( n == 0) {
            l[0] /=area2;
            l[1] /=area2;
            l[2] /=area2;
            if(debug) cout << "   -- in nu "<< nu << " , " << dnu << " :  " << l[1] << " " << l[2] << endl;
            
            return k;
        }
        if(nu<0) nu=k;
        { // calcul dist
            R dn[3];
            int ee[3];
            R de[3];
            for(int j=0; j< n;++j)
            {
                int jj= nl[j], j0=(jj+1)%3, j1=(jj+2)%3;
                
                Rd AB=R2(K[j0],K[j1]),  AP( K[j0],PP), BP(K[j1],PP);
                R la=  (AB,AP);
                R lb= -(AB,BP);
                R lab2 =AB.norme2();
                ee[j]=jj;
                if( la <=0) de[j]=0, dn[j]= AP.norme();
                else if( lb <=0) de[j]=1, dn[j]= BP.norme();
                else de[j]=la/lab2,dn[j]=-l[jj]/sqrt(lab2); //
                
                if(debug) {
                    R dl[3];
                    dl[jj]=0;
                    dl[j0] = 1-de[j];
                    dl[j1] = de[j];
                    Rd PQ(PP,K(Rd(dl[1],dl[2])));
                    R lp = PQ.norme();
                    cout << " \t\t  " << jj<< " " << de[j] <<",  " << dn[j] << " : " <<-l[jj]/(AB.norme()) << " " << AP.norme() << " " << BP.norme() << " : " << lp  << " ??? \n" ;
                }
            }
            int j=0;
            if( n==2 && dn[1]< dn[0]) j=1;
            if( dnu > dn[j] ) {
                nu = k;
                ne= e;
                int jj= nl[j], j0=(jj+1)%3, j1=(jj+2)%3;
                dnu=dn[j];
                dl[jj]=0;
                dl[j0] = 1-de[j];
                dl[j1] = de[j];
            }
            if(debug) {
                Rd Ph=Rd(dl),PQ(PP,(*pTh)[nu](Ph));
                cout<< "     " <<dnu << " (" << nu << " " << k << ") " << dn[j] << " : " << j <<" " << nl[j]
                <<  " n " << n << " " << de[j] << " |" << PQ.norme() << endl;
            }
        }
        
    }
    if (l[ne] > 0)  outside = -1 ; // restart  go to inside ...
    else     outside = (dnu<= delta[i] )? 1: 2;// fare point
    l[0]=dl[0];
    l[1]=dl[1];
    l[2]=dl[2];
    if(debug)   cout << "  -- out nu "<< nu << " "<< ne << " , " << dnu <<" d_i " << delta[i] << " :  "
        << l[1] << " " << l[2] << " "<< outside<<  endl;
    return nu;
}
template<typename Mesh>
int  TrueBorder(const Mesh &Th,typename Mesh::Vertex *P,double *delta)
{
    typedef typename Mesh::Vertex Vertex;
    
    typedef typename Mesh::Element Element;
    typedef typename Mesh::BorderElement BorderElement;
    typedef  typename Mesh::Rd Rd;
    typedef  typename Mesh::BorderElement::RdHat RdHat;
    static const int d = Rd::d;
    static const int dHat = RdHat::d;
  
    int nv =0;
    RdHat GHat(RdHat::diag(1./(dHat+1)));

for(int be=0; be<Th.nbe; ++be)
{
    const BorderElement &E=Th.be(be);
    int e,k = Th.BoundaryElement(be,e);
    {
        int ee=e, kk=  Th.ElementAdj(k,ee);
        if ( kk >=0 || k != kk)
        {
            E(GHat);
            Rd G(E(GHat));
            double l = 0;// 1.5 to be sure .. FH
            for(int i=0; i< BorderElement::nv ;++i)
                l = max(l, Rd(G,E[i]).norme2()) ;
                delta[nv]=l;
                P[nv].lab= Element::ne*k+e;//  element and edge
                (Rd &) P[nv++]=G;
                
                
                }
    }
}
    return nv;
}

template<typename Mesh>
GenericDataFindBoundary<Mesh>::GenericDataFindBoundary(Mesh const * _pTh,int ddebug)
: pTh(_pTh),tree(0), P(bborder ? pTh->nbe: pTh->nt),delta(P.N()),lp(0),debug(ddebug)
{
    const int nvE= Element::nv;
    const int nvB = BorderElement::nv;
    const int nvK = bborder ? nvB : nvE;
    const Mesh &Th = *pTh;
    // extract true Border if d != dHat
    // othesize keep all mesh
    RdHat GHat(RdHat::diag(1./(dHat+1)));

    int nv =0;
    //  warning in case of meshL ,  bord is points  => code bborder stupide..
    if(bborder)
       nv =  TrueBorder(Th,P,delta);
    else
    { //
        for(int k=0; k<Th.nt; ++k)
        {
            const Element& K= Th[k];
            {
                    Rd G(K(GHat));
                    double l = 0;// 1.5 to be sur .. FH
                    for(int i=0; i< Element::nv ;++i)
                        l = max(l, Rd(G,K[i]).norme2()) ;
                    delta[nv]=l;
                    P[nv].lab= k;//  element and edge
                    (Rd &) P[nv++]=G;
                    
                    
                }
            }
        }
    
    //P.resize(nv); no resize because no copy of vertices ...
    delta.resize(nv);
    lp.resize(nv);
    if(debug>7)  gnuplot("dfb0.gp");
    Vertex * P0= &P[0];
    KN<double> d0(delta);
    delta= 0.;
    Rd Pn, Px;
    Th.BoundingBox(Pn, Px);
    double col=0;
    tree=new EF23::GTree<Vertex> (&P[0], Pn, Px,nv);// build quadtree
    
    for(int i=0;i<nv; ++i)
    {
        if(debug>9)   cout << i << " " << d0[i] << endl;
        int lvp=tree->ListNearestVertex(lp,nv, d0[i]*3,P[i]);
        for(int j=0,k; j<lvp; ++j)
        {
            k= lp[j]-P0;
            double dij = Rd(P[i],*lp[j]).norme();
            delta[k]=max(delta[k],d0[i]+dij);
            delta[i]=max(delta[i],d0[i]+dij);

            if(debug>9) cout << k << " "<< delta[k] << ", ";
        }
        if(debug>9) cout << endl;
    }
    if(debug>9)
        for(int i=0;i<nv; ++i)
            cout  << i << " " << d0[i] << " " <<delta[i] << endl;
    if(debug>5)      gnuplot("dfb1.gp");
}
            // Bof Bof pas sur du tout.
/*
template<typename Mesh>
void BuildDataFindBoundary<Mesh>() const
{
    static int count =0;
    if( dfb ==0) {
        dfb=new DataFindBoundary(this);//,count++==0?9:0);
        dfb->debug=0;
    }
    
    }
 */
#endif

