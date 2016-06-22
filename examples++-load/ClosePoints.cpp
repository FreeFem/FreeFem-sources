//ff-c++-LIBRARY-dep:  [flann]

#include <ff++.hpp>
#include <AFunction_ext.hpp>
//#include "HashTable.hpp"
using namespace Fem2D;

static bool debug = false;

class R2close {
public:
    double *data;// minimal points+ delta
    typedef double  * Point;
    int n,nx,offset;
    Point *P;
    const double EPSILON;
    double x0,y0,x1,y1, coef;// boundin box
    R2close() :data(0), n(0), nx(1000000), P(new Point[nx]),EPSILON(1e-6),offset(0)  {InitialiserListe();}
    R2close(double *dd,int mx,double eps=1e-6,int offsett=1) :
    data(dd), n(0), nx(mx), P(new Point[nx]),EPSILON(eps),offset(offsett)
    {InitialiserListe();}
    Point  operator[](int i) const { return P[i];}
    
    int N,m;  // le nombre de code = ncase()%n
    int *head; // pointeur tableau de dimension  m    contenant les tete de liste pour chaque code
    int *next; // pointeur tableau de dimension nx donnant le point suivant de mÍme code
    static const int NotaPoint =-1;
    void InitialiserListe()
    {
        if(data)
        {
            x0=data[0];
            y0=data[1];
            x1=data[2];
            y1=data[3];
            
        }
        else
        {
            x0=0;
            y0=1;
            x1=0;
            y1=1;
        }
        
        coef  = 1./max(x1-x0,y1-y0);
        if(verbosity>10)
            cout << "     bounding box ClosePoints  Pmin=[" << x0 << ", "<< y0 << "], Pmax=[ " << x1 << " " << y1 << "] " << "eps= " <<  EPSILON <<endl;
        
        N=max(sqrt(nx),10.);
        m=max(nx/10,100);
        next = new int[nx];
        head = new int[m];
        for(int i=0; i< m;++i)
            head[i]=NotaPoint;
    }
    int AddSimple(double *p)
    {
        double x = p[0], y=p[offset];
        assert(n<nx);
        P[n]=p;
        
        int k=ncase(x,y)%m;
        assert(k>=0);
        next[n]=head[k];
        head[k]=n;
        if(debug)   cout << "  AddSimple " << n << " <- " << k << " / " << x << " " << y << " / "<< offset << endl;
        return n++;
    }
private:
    Point * Exist(double *p) const
    {
        double x = p[0], y=p[offset];
        
        for(int i=0;i<n;++i)
            if(Equivalent(x,y,P[i]))
                return P+i;
        return 0;
    }
public:
    bool Equivalent(double x0,double y0,const Point Q) const
    {
        return  ((x0-Q[0])*(x0-Q[0])+(y0-Q[offset])*(y0-Q[offset])) < EPSILON*EPSILON;
    }
    
    Point * Exist(double x, double y,int k) const
    {
        for(int i=head[k%m] ; i!=NotaPoint; i=next[i] )
        {   //cout << i << " ";
            if(Equivalent(x,y,P[i]))
                return P+i;
        }
        return 0;
    }
    
    int Add(double *p)
    {
        Point *q=Exist(p);
        if( q) return  q-P; // numero du point
        else return AddSimple(p);
    }
 
    int  FindAll(double x,double y,int * p)
    {
        int nbp=0;
        if(debug)   cout << " Find " << x << " "<< y  << " " << EPSILON << " " << ncase(x,y) << ": " ;
        // double x = p[0], y=p[offset];
        double eps=EPSILON/2;
        Point *q=0;
        int kk[9],k=0,nc;
        for(int i=-1;i<2;++i)
            for(int j=-1;j<2;++j)
            {
                nc = ncase(x+eps*i,y+eps*j);
                //cout <<x+eps*i << " " << y+eps*j << " " << nc << " . ";
                if(nc>=0)
                    for(int i=0;i<k;++i) // remove double cas e..
                        if(kk[i]==nc)
                        {
                            nc=-1;
                            break;
                        }
                if(nc>=0) kk[k++]=nc;
            }
        if(k>4) {
            cout << "   ClosePoints: FindAll Bug ??? : " << k << " : " ;
            for(int i=0; i<k;++i)
                cout << " " << kk[i];
            cout << endl;
        }
        assert(k <=4);
        for(int i=0;i<k;++i)
        {
            int k=kk[i];
            for(int i=head[k%m] ; i!=NotaPoint; i=next[i] )
            {
                if(Equivalent(x,y,P[i]))
                    p[nbp++]=i;
            }
         }
        return nbp;
    }
    
    Point * Find(double x,double y)
    {
      if(debug)   cout << " Find " << x << " "<< y  << " " << EPSILON << " " << ncase(x,y) << ": " ;
       // double x = p[0], y=p[offset];
        double eps=EPSILON/2;
        Point *q=0;
        int kk[9],k=0,nc;
        for(int i=-1;i<2;++i)
            for(int j=-1;j<2;++j)
            {
                nc = ncase(x+eps*i,y+eps*j);
                //cout <<x+eps*i << " " << y+eps*j << " " << nc << " . ";
                if(nc>=0)
                for(int i=0;i<k;++i) // remove double cas e..
                    if(kk[i]==nc)
                    {
                        nc=-1;
                        break;
                    }
                if(nc>=0) kk[k++]=nc;
            }
        if(k>4) {
            cout << "   ClosePoints: Bug ??? : " << k << " : " ;
            for(int i=0; i<k;++i)
                cout << " " << kk[i];
            cout << endl;
        }
        assert(k <=4);
        for(int i=0;i<k;++i)
        {
            q=Exist(x,y,kk[i]) ;
           if(debug)  cout<< " " << kk[i]  ;

            if(q) break;
        }
     if(debug)    cout << " q= "<< q << endl;
        if( q)
            return  q; // numero du point
        else return 0;
    }
    Point * Find(double *p)
    {
        return Find(*p,*(p+offset));
    }
    int AddOpt(double *p)
    {
        Point *q=Find(p);
        if( q)
            return  q-P; // numero du point
        else
            return AddSimple(p);
    }
    
    long ncase(double x,double y)
    {
        if(x <x0 || x >=x1 || y <y0 || y >=y1 ) return -1; // dehors
        else return  long((x-x0)/EPSILON/2) +  long((y-y0)/EPSILON/2)*N;// indice de la case contenant (x,y).
    }
    ~R2close()
    {
        delete []P;
        delete [] head;
        delete [] next;
    }
    // pas de copie de  cette structure
private:
    R2close(const R2close&);
    R2close&operator=(const R2close&);
};

double dist2(int n,double *p,double *q)
{
    double s=0;
    for(int i=0; i<n;++n)
        s += (q[i]-p[i])*(q[i]-p[i]);
    return s;
}
long Hcode(int n,double eps,double *p,double *p0);

KN<long>* CloseTo2(Stack stack,double const & eps,KNM_<double> const &  P,KNM_<double>  const &  Q)
{
    long Pn0=P.N();
    long Pm0=P.M();
    long Qn0=Q.N();
    long Qm0=Q.M();
    double * p0=&(P(0,0));
    int po10 = P.step*P.shapei.step;
    int po01 = P.step*P.shapej.step;
    double * q0=&(Q(0,0));
    int qo10 =Q.step*Q.shapei.step;
    int qo01 =Q.step*Q.shapej.step;

    double x0= P(0,':').min();
    double y0= P(1,':').min();
    double x1= P(0,':').max();
    double y1= P(1,':').max();
    
    // add cc
    double dd= max(x1-x0,y1-y0)*0.01;
    if (dd==0) dd = max(abs(x0),abs(y0))*1e-8;
    if (dd==0) dd=1e-8;
    double data[]={x0-dd,y0-dd,x1+dd,y1+dd};

    R2close S(data,Pm0,eps,po10);
    for (int i=0; i<Pm0;++i)
    {
        if(verbosity>19 )
            cout << i << " :: " << P(0,i) << " " << P(1,i) << endl;
        int j= S.AddSimple(&P(0,i));
    }
    KN<long>* pr= new KN<long>(Qm0);
    for (int i=0; i<Qm0;++i)
    {
        R2close::Point * p= S.Find(Q(0,i),Q(1,i));
        if(p)
          (*pr)[i]= p-S.P;
        else
            (*pr)[i]=-1;
    }

    return Add2StackOfPtr2FreeRC(stack,pr);

}
KN<long>* CloseTo2t(Stack stack,double const & eps,KNM_<double> const &  P,KNM_<double>  const &  Q)
{
    return CloseTo2(stack,eps,P.t(),Q.t());
}

KN<long>* CloseTo(Stack stack,double const & eps,KNM_<double> const &  P,KNM<double> * const &  q,bool tq,bool inv=0)
{
    
    long n0=P.N();
    long m0=P.M();
    if(verbosity>2 )
        cout << " -ClosePoints Size array;   n0 " << n0 << " m0 " << m0  << endl;
    // if(t) std::swap(n0,m0);
    ffassert(n0==2);
    //   bool tq=t; // bofbof ...
    // KNM_<double>  Po =*p;
    // KNM_<double>  Pt =p->t();
    //KNM_<double> P(t ? Pt: Po);
    KNM<double> & Qo=*q;
    double * p0=&(P(0,0));
    
    int offset10 = P.step*P.shapei.step;
    int offset01 = P.step*P.shapej.step;
    if(verbosity>10 )
    {
        cout << "     offset of 0 1 :  "<< offset01  << endl;
        cout << "     offset of 1 0 :  "<< offset10 << endl;
    }
    MeshPoint &mp= *MeshPointStack(stack); // the struct to get x,y, normal , value
    
    double x0= P(0,':').min();
    double y0= P(1,':').min();
    double x1= P(0,':').max();
    double y1= P(1,':').max();
    
    // add cc
    double dd= max(x1-x0,y1-y0)*0.01;
    double data[]={x0-dd,y0-dd,x1+dd,y1+dd};
    KN<long>* pr=0;
    if(inv) pr = new KN<long>(m0);
    
    R2close S(data,m0,eps,offset10);
    for (int i=0; i<m0;++i)
    {
        if(verbosity>19 )
            cout << i << " :: " << P(0,i) << " " << P(1,i) << endl;
        int j= S.AddOpt(&P(0,i));
        if(pr) (*pr)[i]=j;
    }
    
    if( pr==0)
        pr=new KN<long>(S.n);
    if(q)
    {
        if(tq)
            Qo.resize(S.n,n0);
        else
            Qo.resize(n0,S.n);
        KNM_<double> Q(tq ? Qo.t() : Qo);
        for(int i=0;i<S.n;++i)
        {
            double * pi=S[i];
            if(!inv) (*pr)[i]= (pi-p0)/offset01;
            for(int j=0; j<n0;++j)
                Q(j,i)=pi[j*offset10];
        }
    }
    else if(!inv)
    {
        for(int i=0;i<S.n;++i)
        {
            double * pi=S[i];
            (*pr)[i]= (pi-p0)/offset01;
        }
    }
    if(verbosity>2 ) cout << "  - ClosePoint: nb of common points " << m0-S.n;
    //   double *p0=&(*p)(0,0);
    
    
    return Add2StackOfPtr2FreeRC(stack,pr);
}

template<bool inv>
KN<long>* CloseTo(Stack stack,double const & eps,KNM<double> * const &  p,KNM<double> * const &  q)
{
    KNM_<double> P(*p);
    return CloseTo(stack,eps,P,q,false,inv);
    
}

template<bool inv>
KN<long>* CloseTo(Stack stack,double const & eps,Transpose<KNM<double>  *>  const &  p,KNM<double> * const &  q)
{
    KNM_<double> P(p.t->t());
    return CloseTo(stack,eps,P,q,false,inv);
    
}

template<bool inv>
KN<long>* CloseTo(Stack stack,double const & eps,Transpose<KNM<double>  *>  const &  p,Transpose<KNM<double>  *>  const &  q)
{
    KNM_<double> P(p.t->t());
    return CloseTo(stack,eps,P,q,true,inv);
    
}

template<bool inv>
KN<long>* CloseTo(Stack stack,double const & eps,KNM<double> * const &  p,Transpose<KNM<double>  *> const &  q)
{
    KNM_<double> P(*p);
    return CloseTo(stack,eps,P,q.t,true,inv);
    
}

template<bool inv>
KN<long>* CloseTo(Stack stack,double const & eps,KNM<double> * const &  p)
{
    KNM_<double> P(*p);
    return CloseTo(stack,eps,P,0,false,inv);
    
}

template<bool inv>
KN<long>* CloseTo(Stack stack,double const & eps,KNM_<double>  const &  p)
{
    KNM_<double> P(p);
    if(verbosity>5)
        cout << " CloseTo KNM_ "<< P.N() << " " << P.M() << endl;
    return CloseTo(stack,eps,P,0,false,inv);
    
}

template<bool inv>
KN<long>* CloseTo(Stack stack,double const & eps,Transpose<KNM<double>  *>  const &  p)
{
    KNM_<double> P(p.t->t());
    return CloseTo(stack,eps,P,0,false,inv);
}

template<bool inv>
KN<long>* CloseTo(Stack stack,double const & eps,pmesh  const &  pTh,KNM<double> * const &  pq)
{
    ffassert(pTh && pq);
    const  Mesh & Th=*pTh;
    const KNM<double> Q=*pq;
    int np=Q.N();
    assert(Q.M()>=2);
    KN<long>* pr=new KN<long>(inv? Th.nv :np);
    KN<int> b(Th.nv);
    b=0;
    for (int i=0;i<Th.nv;i++)
        if(Th(i).lab!=0) b[i]=1;
    for (int i=0;i<Th.neb;i++)
    {
        b[ Th(Th.bedges[i][0])]=1;
        b[ Th(Th.bedges[i][1])]=1;
    }
    
    
    *pr=-1L;
    R2 Pmin,Pmax;
    Th.BoundingBox(Pmin,Pmax);
    int nv=b.sum();
    if(verbosity>9) cout << " Th.nv " <<Th.nv << " " << nv << "/ " << Pmin << " "<< Pmax <<endl;
    FQuadTree * quadtree = new FQuadTree(pTh,Pmin,Pmax,nv); //  put all the old vertices in the quadtree
    //  copy of the old vertices
    for (int i=0;i<Th.nv;i++)
    {
        if(b[i])
        {  cout << i << " " << Th.vertices[i] << endl;
            quadtree->Add(Th.vertices[i]);
        }
        
    }
    cout << " After quadterr" << endl;
    for (int j=0;j<np; ++j)
    {
        R2 P(Q(j,0),Q(j,1));
        
        Vertex * pV=quadtree->ToClose(P,eps);
        if(pV)
        {
            pV=quadtree->NearestVertex(P);
            long k= Th(pV);
            // cout << j << " " << k << " "  << P << " " << pV <<endl;
            if(inv) (*pr)[k]=j;
            else (*pr)[j]=k;
            
        }
    }
    delete quadtree;
    
    return Add2StackOfPtr2FreeRC(stack,pr);
}
bool InterAB_Disq(R2 A,R2 B,R2 C,double r)
{
    double ACB = det(A,C,B);
    if (ACB <0) return false;
    R2 AB(A,B);
    R2 AC(A,C),CB(C,B);
    double LAB= AB.norme();
    double H = ACB/4/LAB;
    if( H > r) return false;
    double AC2 = AC.norme2(),r2=r*r;
    if(AC2 < r2) return true;
    double CB2 = CB.norme2();
    if(CB2 < r2) return true;
    double ABAC=(AB,AC),ABCB=(AB,CB);
    return   ABAC >0 && ABCB >0;
}

void Add(KN<long>& I,int i)
{
    int n= I.N();
    int m =-I[n-1];
    if( m <=0) {
        m=n;
        n = m*2;
        I.resize(n);
        I[n-1]=-m-1;
    }
    else --m;
   if(debug)   cout << " add "<< m <<" "<< i <<" "<< n << endl;
    I[m]=i;
    if(m<n-1) --I[n-1];
}
void Clean(KN<long> &I)
{
    int n= I.N();
    int m =-I[n-1];
    if( m >=0) I.resize(m-1);
}
long Voisinage( KNM_<double> const &  P ,KNM_<double> const &  Q, double const &  eps,KN<KN<long> >  * const & IJ)
{
    debug = (verbosity>999);
    int np=P.N();
    int nq=Q.N();
    int mp=P.M();
    int mq=Q.M();

    double *p= &P(0,0);
    double *q= &Q(0,0);
    int offset01= P.step*P.shapej.step;;
    int offset10= P.step*P.shapei.step;;
    double * q0=&(Q(0,0));
    int qoffset01 = Q.step*Q.shapej.step;;
    int qoffset10 = Q.step*Q.shapei.step;;
    ffassert( mp ==2);
    ffassert( mq ==2);
    KN<int> lp(np);
    
    IJ->resize(nq);
    for (int i=0;i<nq ;i++)
    {
        (*IJ)[i].resize(2);
        (*IJ)[i][0]=-1;
        (*IJ)[i][1]=-1;
    }
    if(verbosity>99)  cout << " offset01 " <<offset01 << " " << offset10 << " p" << p << " " << np << " "<< P.M()  <<  endl;
    // store - the size in last value of IJ[i][j]
    double data[4];
    data[0] = data[2] = p[0];
    data[1] = data[3] = p[offset01];
    for(int i=0,k=0; i<np; ++i, k+= offset10)
    {
        data[0] = min(data[0], p[k]);
        data[2] = max(data[2], p[k]);
        data[1] = min(data[1], p[k+offset01]);
        data[3] = max(data[3], p[k+offset01]);
        
    }
    double eps2 = eps+eps;
    data[0] -=  eps2;
    data[2] +=  eps2;
    data[1] -=  eps2;
    data[3] +=  eps2;
    
    R2close SP(data,np,eps,offset01);
    for (int i=0; i<np;++i)
        SP.AddSimple(&P(i,0));
    
    for (int j=0; j<nq;++j)
    {
        int nlp =SP.FindAll(Q(j,0),Q(j,1),lp);
        for(int k=0; k<nlp;++k)
        {  //
            int i = lp[k];
            if(verbosity>99)  cout << " Add to i=" << i << " -> j " << j << endl;
            Add((*IJ)[i],j);
        }
    }
    for (int j=0; j<nq;++j)
        Clean((*IJ)[j]);
    debug = 0;
    return 0;
}
#ifdef WITH_flann
#include <flann/flann.hpp>
long ff_flann_search( KNM_<double> const &  P ,KNM_<double> const &  Q, double const &  eps,KN<KN<long> >  * const & pIJ)
{
    KN<KN<long> > &IJ=*pIJ;
    int np=P.M();
    int nq=Q.M();
    
    int mp=P.N();
    int mq=Q.N();
    int nn=nq;
    double *p= &P(0,0);
    double *q= &Q(0,0);
    int offset01= P.step*P.shapej.step;;
    int offset10= P.step*P.shapei.step;;
    double * q0=&(Q(0,0));
    int qoffset01 = Q.step*Q.shapej.step;;
    int qoffset10 = Q.step*Q.shapei.step;;
    cout << np << " " << nq << " po 01, 10 " << offset01 << " " << offset10 << " " << &P(0,1) - p
                            << " qo 01 , 01 =: "  << qoffset01 << " "<< qoffset10 << endl;
    ffassert( mp == mq && offset10==1 && qoffset10==1);
    //ffassert( mq ==2);
    
    
    IJ.resize(nq);
    cout << np << " " << mq << endl;
    cout << nq << " " << mq << endl;
    flann::Matrix<double> dataset(p,np,nq);
    flann::Matrix<double> query(q,nq,mq);
    std::vector< std::vector<int> > indices;
    std::vector<std::vector<double> > dists;
    flann::SearchParams params(128);
    params.eps= eps;
    flann::Index<flann::L2<double> > index(dataset, flann::KDTreeIndexParams(4));
    index.buildIndex();
    index.radiusSearch(query, indices, dists, eps , params);
    for (int j=0; j<indices.size() ;++j)
    {
         IJ[j].resize(indices[j].size());
        for(int k=0;k<indices[j].size(); ++k)
           IJ[j][k]=indices[j][k];
    }
   for (int j=0; j<indices.size() ;++j)
    {
        int k = j*mq;
        cout << j << " [ " << q[k] << " " << q[k+1] <<"] : ";
        for(int k=0;k<indices[j].size(); ++k)
            cout <<  indices[j][k] << " " << dists[j][k] << ", ";
        cout << endl;
    }
    
    return 0;
}


#endif
static void init()
{
    #ifdef WITH_flann
    Global.Add("radiusSearch","(",new OneOperator4_<long, KNM_<double> , KNM_<double>  ,double,KN<KN<long> > *   >(ff_flann_search));
#endif
    Global.Add("Voisinage","(",new OneOperator4_<long, KNM_<double> , KNM_<double>  ,double,KN<KN<long> > *   >(Voisinage));
    Global.Add("ClosePoints2","(",new OneOperator3s_<KN<long>*,double, KNM_<double> , KNM_<double>   >(CloseTo2));
  //s  Global.Add("ClosePoints2t","(",new OneOperator3s_<KN<long>*,double, KNM_<double> , KNM_<double>   >(CloseTo2t));
    
    //   numbering ..
   // Global.Add("ClosePoints","(",new OneOperator2s_<KN<long>*,double, Transpose<KNM<double>  *>  >(CloseTo<false>));
   // Global.Add("ClosePoints","(",new OneOperator2s_<KN<long>*,double, KNM<double>   *  >(CloseTo<false>,1));
    Global.Add("ClosePoints","(",new OneOperator2s_<KN<long>*,double, KNM_<double>     >(CloseTo<false>,0));
    //Global.Add("ClosePoints","(",new OneOperator3s_<KN<long>*,double, Transpose<KNM<double>  *> ,KNM<double> * >(CloseTo<false>));
   // Global.Add("ClosePoints","(",new OneOperator3s_<KN<long>*,double, Transpose<KNM<double>  *> ,Transpose<KNM<double>  *> >(CloseTo<false>));
   // Global.Add("ClosePoints","(",new OneOperator3s_<KN<long>*,double,KNM<double> * ,KNM<double> * >(CloseTo<false>));
   // Global.Add("ClosePoints","(",new OneOperator3s_<KN<long>*,double,KNM<double> * ,Transpose<KNM<double>  *> >(CloseTo<false>));
    Global.Add("ClosePoints","(",new OneOperator3s_<KN<long>*,double,pmesh,KNM<double> * >(CloseTo<false>));
    // inv  numbering ..
   // Global.Add("ClosePoints1","(",new OneOperator2s_<KN<long>*,double, Transpose<KNM<double>  *>  >(CloseTo<true>));
   // Global.Add("ClosePoints1","(",new OneOperator2s_<KN<long>*,double, KNM<double>   *  >(CloseTo<true>,1));
    Global.Add("ClosePoints1","(",new OneOperator2s_<KN<long>*,double, KNM_<double>     >(CloseTo<true>,0));
  //  Global.Add("ClosePoints1","(",new OneOperator3s_<KN<long>*,double, Transpose<KNM<double>  *> ,KNM<double> * >(CloseTo<true>));
  //  Global.Add("ClosePoints1","(",new OneOperator3s_<KN<long>*,double, Transpose<KNM<double>  *> ,Transpose<KNM<double>  *> >(CloseTo<true>));
  //  Global.Add("ClosePoints1","(",new OneOperator3s_<KN<long>*,double,KNM<double> * ,KNM<double> * >(CloseTo<true>));
  //  Global.Add("ClosePoints1","(",new OneOperator3s_<KN<long>*,double,KNM<double> * ,Transpose<KNM<double>  *> >(CloseTo<true>));
    Global.Add("ClosePoints1","(",new OneOperator3s_<KN<long>*,double,pmesh,KNM<double> * >(CloseTo<true>));
    
}

LOADFUNC(init);
