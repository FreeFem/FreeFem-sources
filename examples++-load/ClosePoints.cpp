#include <ff++.hpp>
#include <AFunction_ext.hpp>
//#include "HashTable.hpp"
using namespace Fem2D;



class R2close {
public:
    double *data;// minimal points+ delta
    typedef double  * Point;
    int n,nx,offset;
    Point *P;
    const double EPSILON;
    double x0,y0,x1,y1, coef;// boundin box
    R2close() :data(0), n(0), nx(1000000), P(new Point[nx]),EPSILON(1e-6)  {InitialiserListe();}
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
        
        coef  = 1./max(x1-x0,y1-y1);
        cout << " bb  [" << x0 << ", "<< y0 << "], [ " << x1 << " " << y1 << "], " << coef << endl;
        
        N=max(sqrt(nx),10.);
        m=max(nx/10,10);
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
        return n++;
    }
    Point * Exist(double *p) const
    {
        double x = p[0], y=p[offset];

        for(int i=0;i<n;++i)
            if(Equivalent(x,y,P[i]))
                return P+i;
        return 0;
    }
    
    bool Equivalent(double x0,double y0,const Point Q) const
    {
        return  max(abs(x0-Q[0]),abs(y0-Q[offset])) < EPSILON;
    }
    
    Point * Exist(double x, double y,int k) const
    {
        for(int i=head[k%m] ; i!=NotaPoint; i=next[i] )
            if(Equivalent(x,y,P[i]))
                return P+i;
        return 0;
    }
    
    int Add(double *p)
    {
        Point *q=Exist(p);
        if( q) return  q-P; // numero du point
        else return AddSimple(p);
    }
    
    int AddOpt(double *p)
    {
        double x = p[0], y=p[offset];
        double eps=EPSILON/2;
        Point *q=0;
        int kk[9],k=0,nc;
        for(int i=-1;i<2;++i)
            for(int j=-1;j<2;++j)
            {
                nc = ncase(x+eps*i,y+eps*j);
                for(int i=0;i<k;++i) // remove double cas e..
                    if(kk[i]==nc)
                    {
                        nc=-1;
                        break;
                    }
                if(nc>=0) kk[k++]=nc;
                
            }
        assert(k <=4);
        for(int i=0;i<k;++i)
        {
            q=Exist(x,y,kk[i]) ;
            if(q) break;
        }
        if( q)
            return  q-P; // numero du point
        else
            return AddSimple(p);
    }
    
    long ncase(double x,double y)
    {
        if(x <x0 || x >=x1 || y <y0 || y >=y1 ) return -1; // dehors
        else return  long((x-x0)*coef*N)*N +  long((y-y0)*coef*N);// indice de la case contenant (x,y).
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

KN<long>* CloseTo(Stack stack,double const & eps,bool t,KNM<double> * const &  p,KNM<double> * const &  q,bool tq)
{
    
  long n0=p->N();
  long m0=p->M();
  cout << " n0 " << n0 << " m0 " << m0 << " t = " << t << endl;
  if(t) std::swap(n0,m0);
  ffassert(n0==2);
 //   bool tq=t; // bofbof ...
  KNM_<double>  Po =*p;
  KNM_<double>  Pt =p->t();
  KNM_<double> P(t ? Pt: Po);
  KNM<double> & Qo=*q;
   double * p0=&(P(0,0));
   int offset10 =( &(P(1,0)) - p0);
   int offset01 =( &(P(0,1)) - p0);
    cout << " offset of 0 1 :  "<< offset01  << endl;
    cout << " offset of 1 0 :  "<< offset10 << endl;
  MeshPoint &mp= *MeshPointStack(stack); // the struct to get x,y, normal , value
    
    double x0= P(0,':').min();
    double y0= P(1,':').min();
    double x1= P(0,':').max();
    double y1= P(1,':').max();
    
    // add cc
    double dd= max(x1-x0,y1-y0)*0.01;
    double data[]={x0-dd,y0-dd,x1+dd,y1+dd};
    
    
  R2close S(data,m0,eps,offset10);
  for (int i=0; i<m0;++i)
    {
        cout << i << " :: " << P(0,i) << " " << P(1,i) << endl;
        S.AddOpt(&P(0,i));
    }
    
    KN<long>* pr=new KN<long>(S.n);
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
       (*pr)[i]= (pi-p0)/offset01;
       for(int j=0; j<n0;++j)
           Q(j,i)=pi[j*offset10];
   }
    }
    else
    {
        for(int i=0;i<S.n;++i)
        {
            double * pi=S[i];
            (*pr)[i]= (pi-p0)/offset01;
        }
    }
    
 //   double *p0=&(*p)(0,0);
    
    
  return Add2StackOfPtr2FreeRC(stack,pr);
}

KN<long>* CloseTo(Stack stack,double const & eps,KNM<double> * const &  p,KNM<double> * const &  q)
{
    return CloseTo(stack,eps,false,p,q,false);
    
}
KN<long>* CloseTo(Stack stack,double const & eps,Transpose<KNM<double>  *>  const &  p,KNM<double> * const &  q)
{
    return CloseTo(stack,eps,true,p.t,q,false);
    
}
KN<long>* CloseTo(Stack stack,double const & eps,Transpose<KNM<double>  *>  const &  p,Transpose<KNM<double>  *>  const &  q)
{
    return CloseTo(stack,eps,true,p.t,q,true);
    
}

KN<long>* CloseTo(Stack stack,double const & eps,KNM<double> * const &  p,Transpose<KNM<double>  *> const &  q)
{
    return CloseTo(stack,eps,false,p,q.t,true);
    
}

KN<long>* CloseTo(Stack stack,double const & eps,KNM<double> * const &  p)
{
    return CloseTo(stack,eps,false,p,0,false);
    
}

KN<long>* CloseTo(Stack stack,double const & eps,Transpose<KNM<double>  *>  const &  p)
{
    return CloseTo(stack,eps,true,p.t,0,false);
}


KN<long>* CloseTo(Stack stack,double const & eps,pmesh  const &  pTh,KNM<double> * const &  pq)
{
  ffassert(pTh && pq);
 const  Mesh & Th=*pTh;
  const KNM<double> Q=*pq;
  int np=Q.N();
  assert(Q.M()>=2);
  KN<long>* pr=new KN<long>(np);
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
    cout << " Th.nv " <<Th.nv << " " << nv << "/ " << Pmin << " "<< Pmax <<endl;
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
            cout << j << " " << k << " "  << P << " " << pV <<endl;

          (*pr)[j]=k;
        }
    }
    delete quadtree;
    
  return Add2StackOfPtr2FreeRC(stack,pr);
}


void init(){
    Global.Add("ClosePoints","(",new OneOperator2s_<KN<long>*,double, Transpose<KNM<double>  *>  >(CloseTo));
    Global.Add("ClosePoints","(",new OneOperator2s_<KN<long>*,double, KNM<double>   *  >(CloseTo));
   
   Global.Add("ClosePoints","(",new OneOperator3s_<KN<long>*,double, Transpose<KNM<double>  *> ,KNM<double> * >(CloseTo));
    Global.Add("ClosePoints","(",new OneOperator3s_<KN<long>*,double, Transpose<KNM<double>  *> ,Transpose<KNM<double>  *> >(CloseTo));
   Global.Add("ClosePoints","(",new OneOperator3s_<KN<long>*,double,KNM<double> * ,KNM<double> * >(CloseTo));
    Global.Add("ClosePoints","(",new OneOperator3s_<KN<long>*,double,KNM<double> * ,Transpose<KNM<double>  *> >(CloseTo));
  Global.Add("ClosePoints","(",new OneOperator3s_<KN<long>*,double,pmesh,KNM<double> * >(CloseTo));
}

LOADFUNC(init);
