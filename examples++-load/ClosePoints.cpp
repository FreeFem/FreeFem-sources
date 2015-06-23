#include <ff++.hpp>
#include <AFunction_ext.hpp>
using namespace Fem2D;


double dist2(int n,double *p,double *q)
{
    double s=0;
    for(int i=0; i<n;++n)
        s += (q[i]-p[i])*(q[i]-p[i]);
    return s;
}
long Hcode(int n,double eps,double *p,double *p0);

KN<long>* CloseTo(Stack stack,double const & eps,KNM<double> * const &  p,KNM<double> * const &  q)
{
    
  long n=q->N();
  KN<long>* pr=new KN<long>(q->N());
  MeshPoint &mp= *MeshPointStack(stack); // the struct to get x,y, normal , value
  
  //HashTable<long,int> table(n);
    double *p0=&(*p)(0,0);
    
    
  return Add2StackOfPtr2FreeRC(stack,pr);
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
  //  Global.Add("ClosePoints","(",new OneOperator3s_<KN<long>*,double,KNM<double> * ,KNM<double> * >(CloseTo));
  Global.Add("ClosePoints","(",new OneOperator3s_<KN<long>*,double,pmesh,KNM<double> * >(CloseTo));
}

LOADFUNC(init);
