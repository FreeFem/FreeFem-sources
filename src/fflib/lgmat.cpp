#ifdef __MWERKS__
#pragma optimization_level 0
//#pragma inline_depth(0)
#endif

#include  <cmath>
#include  <iostream>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
#include "rgraph.hpp"
#include <cstdio>
#include "MatriceCreuse_tpl.hpp"

//#include "fem3.hpp"
#include "MeshPoint.hpp"
#include <complex>
#include "Operator.hpp" 

#include <set>
#include <vector>

#include "lex.hpp"
#include "lgfem.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"
#include "CGNL.hpp"
namespace bamg { class Triangles; }
namespace Fem2D { void DrawIsoT(const R2 Pt[3],const R ff[3],const RN_ & Viso); }

#include "BamgFreeFem.hpp"



// Debut FH Houston -------- avril 2004 ---------
//  class for operator on sparce matrix 
//   ------------------------------------
//  pour le addition de matrice ----matrice creuse in french)
//  MatriceCreuse    real class for matrix sparce
//  Matrice_Creuse   class for matrix sparce +  poiteur on FE space def 
//         to recompute matrice in case of mesh change
//  list<pair<R,MatriceCreuse<R> *> > * is liste of 
//  \sum_i a_i*A_i  where a_i is a scalare and A_i is a Sparce matrix
//

template<class R> 
list<pair<R,MatriceCreuse<R> *> > * to(Matrice_Creuse<R> * M)
{
  list<pair<R,MatriceCreuse<R> *> >  * l=new list<pair<R,MatriceCreuse<R> *> >;
   l ->push_back(make_pair<R,MatriceCreuse<R> *>(1,M->A));
  return  l;
}

template<class R> 
struct Op2_ListCM: public binary_function<R,Matrice_Creuse<R> *,list<pair<R,MatriceCreuse<R> *> > *> 
 { 
   typedef pair<R,MatriceCreuse<R>*> P;
   typedef list<P> L;
   typedef L * RR;
   typedef R  AA;
   typedef Matrice_Creuse<R> * BB;
   
 static  RR f(const   AA & a,const   BB & b)  
  {
    RR r=  new list<P> ;
    r ->push_back(make_pair<R,MatriceCreuse<R> *>(a,b->A));
    return r;}
};

template<class R> 
struct Op2_ListMC: public binary_function<Matrice_Creuse<R> *,R,list<pair<R,MatriceCreuse<R> *> > *> 
 { 
   typedef pair<R,MatriceCreuse<R>*> P;
   typedef list<P> L;
   typedef L * RR;
   typedef R  AA;
   typedef Matrice_Creuse<R> * BB;
   
 static  RR f(const   BB & b,const   AA & a)  
  {
    RR r=  new list<P> ;
    r ->push_back(make_pair<R,MatriceCreuse<R> *>(a,b->A));
    return r;}
};


template<class R> 
struct Op2_ListCMCMadd: public binary_function<list<pair<R,MatriceCreuse<R> *> > *,
                                               list<pair<R,MatriceCreuse<R> *> > *,
                                               list<pair<R,MatriceCreuse<R> *> > *  >
{  //  ... + ...
   typedef pair<R,MatriceCreuse<R>*> P;
   typedef list<P> L;
   typedef L * RR;

  static   RR f(const RR & a,const RR & b)  
  { 
    a->insert(a->end(),b->begin(),b->end());
    delete b;
    return a;
  }    
    
};

template<class R> 
struct Op2_ListMCMadd: public binary_function<Matrice_Creuse<R> *,
                                              list<pair<R,MatriceCreuse<R> *> > *,                                               
                                               list<pair<R,MatriceCreuse<R> *> > *  >
{  //  M + ....
   typedef pair<R,MatriceCreuse<R>*> P;
   typedef list<P> L;
   typedef L * RR;
   typedef Matrice_Creuse<R> * MM;

  static   RR f(const MM & a,const RR & b)  
  { 
    
    b->push_front(make_pair<R,MatriceCreuse<R> *>(R(1.),a->A));
    return b;
  }
    
    
};

template<class R> 
struct Op2_ListCMMadd: public binary_function< list<pair<R,MatriceCreuse<R> *> > *,                                                                                              
                                               Matrice_Creuse<R> * ,
                                               list<pair<R,MatriceCreuse<R> *> > *>
{  //   .... + M
   typedef pair<R,MatriceCreuse<R>*> P;
   typedef list<P> L;
   typedef L * RR;
   typedef Matrice_Creuse<R> * MM;

  static   RR f(const RR & a,const MM & b)  
  { 
    
    a->push_back(make_pair<R,MatriceCreuse<R> *>(R(1.),b->A));
    return a;
  }
    
    
};

template<class R> 
struct Op2_ListMMadd: public binary_function< Matrice_Creuse<R> *,
                                              Matrice_Creuse<R> * ,
                                              list<pair<R,MatriceCreuse<R> *> > *>
{  //  M + M
   typedef pair<R,MatriceCreuse<R>*> P;
   typedef list<P> L;
   typedef L * RR;
   typedef Matrice_Creuse<R> * MM;

  static   RR f(const MM & a,const MM & b)  
  { 
    L * l=to(a);
    l->push_back(make_pair<R,MatriceCreuse<R> *>(R(1.),b->A));
    return l;
  }
    
    
};

// Fin Add FH Houston --------


void buildInterpolationMatrix(MatriceMorse<R> * m,const FESpace & Uh,const FESpace & Vh,void *data)
{  //  Uh = Vh 

  int op=op_id; //  value of the function
  bool transpose=false;
  bool inside=false;
  if (data)
   {
     int * idata=static_cast<int*>(data);
     transpose=idata[0];
     op=idata[1];
     inside=idata[2];
     ffassert(op>=0 && op < 4);
   }
  if(verbosity>2) 
    {
      cout << " -- buildInterpolationMatrix   transpose =" << transpose << endl
           << "              value, dx , dy          op = " << op << endl
           << "                            just  inside = " << inside << endl;
    }
  using namespace Fem2D;
  int n=Uh.NbOfDF;
  int mm=Vh.NbOfDF;
  if(transpose) Exchange(n,mm);
  m->symetrique = false;
  m->dummy=false;
  m->a=0;
  m->lg=0;
  m->cl=0;
  m->nbcoef=0;
  m->n=n;
  m->m=mm;
  int n1=n+1;
  const  Mesh & ThU =Uh.Th; // line 
  const  Mesh & ThV =Vh.Th; // colunm
  bool samemesh =  &Uh.Th == &Vh.Th;  // same Mesh
  int thecolor =0;
  int nbn_u = Uh.NbOfNodes;
  int nbn_v = Vh.NbOfNodes;
  
  int nbcoef =0;
  
  KN<int> color(ThV.nt);
  KN<int> mark(n);
  mark=0;
  
  int *lg = new int [n1];
  int * cl = 0;
  double *a=0;
  
  color=thecolor++;
  FElement Uh0 = Uh[0];
  FElement Vh0 = Vh[0];
  
  FElement::aIPJ ipjU(Uh0.Pi_h_ipj()); 
  FElement::aR2  PtHatU(Uh0.Pi_h_R2()); 
  
  FElement::aIPJ ipjV(Vh0.Pi_h_ipj()); 
  FElement::aR2  PtHatV(Vh0.Pi_h_R2()); 
  
  int nbdfVK= Vh0.NbDoF();
  int nbdfUK= Uh0.NbDoF();
  int NVh= Vh0.N;
  int NUh= Uh0.N;
  
  ffassert(NVh==NUh); 
  
  
  int nbp= PtHatU.N(); // 
  int nbc= ipjU.N(); // 
  KN<R2> PV(nbp);   //  the PtHat in ThV mesh
  KN<int> itV(nbp); // the Triangle number
  KN<bool> intV(nbp); // ouside or not 
  KN<R>   AipjU(ipjU.N());
  
  KNM<R> aaa(nbp,nbdfVK); 


   const R eps = 1.0e-10;
   const int sfb1=Vh0.N*last_operatortype*Vh0.NbDoF();
   KN<R> kv(sfb1*nbp);
   R * v = kv;
    KN<int> ik(nbp); // the Triangle number
//   KNMK_<R> fb(v+ik[i],VhV0.NbDoF,VhV0.N,1); //  the value for basic fonction

   bool whatd[last_operatortype];
   for (int i=0;i<last_operatortype;i++) 
     whatd[i]=false;
   whatd[op]=true; // the value of function
   KN<bool> fait(Uh.NbOfDF);
   fait=false;
  map< pair<int,int> , double > sij; 
  
  for (int step=0;step<2;step++) 
   {
      
	  for (int it=0;it<ThU.nt;it++)
	    {
	      thecolor++; //  change the current color
	      const Triangle & TU(ThU[it]);
	      FElement KU(Uh[it]);
	      KU.Pi_h(AipjU);
	      int nbkU = 0;
	      if (samemesh)
	        {
	          nbkU = 1; 
	          PV = PtHatU;
	          itV = it;
	        }
	      else 
	       {  
	          const Triangle *ts=0;
	          for (int i=0;i<nbp;i++)
	            {
	              bool outside;
	              ts=ThV.Find(TU(PtHatU[i]),PV[i],outside,ts); 
	              itV[i]= ThV(ts);
	              intV[i]=outside && inside; //  ouside and inside flag 
	            }
	       }
	       
	      
	      for (int p=0;p<nbp;p++)
	         { 

	              KNMK_<R> fb(v+p*sfb1,nbdfVK,NVh,last_operatortype); // valeur de fonction de base de Vh 
	              // ou:   fb(idf,j,0) valeur de la j composante de la fonction idf 
	              Vh0.tfe->FB(whatd,ThV,ThV[itV[p]],PV[p],fb);  
	          }

	      for (int i=0;i<ipjU.N();i++) 
	          { // pour tous le terme 
	           const FElement::IPJ &ipj_i(ipjU[i]);
	           assert(ipj_i.j==0); // car  Vh.N=0
	           int dfu = KU(ipj_i.i); // le numero de df global 
	           if(fait[dfu]) continue;
	           int j = ipj_i.j; // la composante
	           int p=ipj_i.p;  //  le points
	           if (intV[p]) continue; //  ouside and inside flag => next 
	           R aipj = AipjU[i];
	           FElement KV(Vh[itV[p]]);
	           
	            KNMK_<R> fb(v+p*sfb1,nbdfVK,NVh,last_operatortype); 
	            KN_<R> fbj(fb('.',j,op)); 
	            
	           for (int idfv=0;idfv<nbdfVK;idfv++) 
	             if (Abs(fbj[idfv])>eps) 
	              {
	                int dfv=KV(idfv);
	                int i=dfu, j=dfv;
	                if(transpose) Exchange(i,j);
	                // le term dfu,dfv existe dans la matrice
	                R c= fbj[idfv]*aipj;
	                //  cout << " Mat inter " << i << " , "<< j << " = " << c << " " <<step << " " << it << " " <<  endl; 
	                if(Abs(c)>eps)
	                 //
	                   sij[make_pair(i,j)] = c;
	                /*   
	                 if(step==0)
	                   sij.insert(make_pair(i,j));
	                 else	                  
	                   (*m)(i,j)=c;
                        */
	              }
	              
	                      
	          }
	          
	      for (int df=0;df<KU.NbDoF();df++) 
	          {  
	           int dfu = KU(df); // le numero de df global 
	           fait[dfu]=true;
	           }
	       
	       
	    }
	    if (step==0)
	     {
	         nbcoef = sij.size();
	         cl = new int[nbcoef];
	         a = new double[nbcoef];
	         int k=0;
	         for(int i=0;i<n1;i++)
	           lg[i]=0;
	          
	         for (map<pair<int,int>, double >::iterator kk=sij.begin();kk!=sij.end();++kk)
	          { 
	            int i= kk->first.first;
	            int j= kk->first.second;
	           // cout << " Mat inter " << i << " , "<< j  << endl;
	            cl[k]=j;
	            a[k]= kk->second;	            
	            lg[i+1]=++k;
	          }
	         assert(k==nbcoef);
	         //  on bouche les ligne vide   lg[i]=0; 
	         //  lg est un tableau croissant  =>
	         for(int i=1;i<n1;i++)
	              lg[i]=max(lg[i-1],lg[i]) ;
	         m->lg=lg;
             m->cl=cl;
             m->a=a;
             m->nbcoef=nbcoef;
             fait=false;
	     }
    }
    sij.clear();
  //assert(0); // a faire to do
}

AnyType SetMatrixInterpolation(Stack stack,Expression emat,Expression einter)
{
  using namespace Fem2D;
  
  Matrice_Creuse<R> * sparce_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
  const MatrixInterpolation::Op * mi(dynamic_cast<const MatrixInterpolation::Op *>(einter));
  ffassert(einter);

  pfes * pUh = GetAny< pfes * >((* mi->a)(stack));
  pfes * pVh = GetAny<  pfes * >((* mi->b)(stack));
  FESpace * Uh = **pUh;
  FESpace * Vh = **pVh;
  ffassert(Vh);
  ffassert(Uh);
  int data[ MatrixInterpolation::Op::n_name_param];
  data[0]=mi->arg(0,stack,false); // transpose not
  data[1]=mi->arg(1,stack,(long) op_id); ; // get just value
  data[2]=mi->arg(2,stack,false); ; // get just value
  
  sparce_mat->pUh=pUh;
  sparce_mat->pVh=pVh;
  sparce_mat->typemat=TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
  sparce_mat->A.master(new MatriceMorse<R>(*Uh,*Vh,buildInterpolationMatrix,data));

  return sparce_mat;
}

template<class RA,class RB,class RAB>
AnyType ProdMat(Stack stack,Expression emat,Expression prodmat)
{
  using namespace Fem2D;
  
  Matrice_Creuse<RAB> * sparce_mat =GetAny<Matrice_Creuse<RA>* >((*emat)(stack));
  const Matrix_Prod<RA,RB>  AB = GetAny<Matrix_Prod<RA,RB> >((*prodmat)(stack));
  sparce_mat->pUh=AB.A->pUh;
  sparce_mat->pVh=AB.B->pVh;
  MatriceMorse<RA> *mA= AB.A->A->toMatriceMorse(AB.ta);
  MatriceMorse<RB> *mB= AB.B->A->toMatriceMorse(AB.tb);
  if( !mA && ! mB) ExecError(" Sorry error: in MatProd,  pb trans in MorseMat");
  if( mA->m != mB->n) {
    cerr << " -- Error dim ProdMat A*B : tA =" << AB.ta << " = tB " << AB.tb << endl;
    cerr << " --MatProd " << mA->n<< " "<< mA->m << " x " << mB->n<< " "<< mB->m <<  endl;
    ExecError(" Wrong mat dim in MatProd");
  }
  MatriceMorse<RAB> *mAB=new MatriceMorse<RA>();
  mA->prod(*mB,*mAB);
  
  sparce_mat->typemat=(mA->n == mB->m) ? TypeSolveMat(TypeSolveMat::GMRES) : TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
  sparce_mat->A.master(mAB);
  delete mA;
  delete mB;
  return sparce_mat;
}

template<class R>
AnyType CombMat(Stack stack,Expression emat,Expression combMat)
{
  using namespace Fem2D;
  
  Matrice_Creuse<R> * sparce_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
  list<pair<R,MatriceCreuse<R> *> > *  lcB = GetAny<list<pair<R,MatriceCreuse<R> *> >*>((*combMat)(stack));
  sparce_mat->pUh=0;
  sparce_mat->pVh=0; 
   MatriceCreuse<R> * AA=BuildCombMat<R>(*lcB);
  sparce_mat->A.master(AA);
  sparce_mat->typemat=(AA->n == AA->m) ? TypeSolveMat(TypeSolveMat::GMRES) : TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
  delete lcB;
  return sparce_mat;
}

template<class R>
AnyType DiagMat(Stack stack,Expression emat,Expression edia)
{
  using namespace Fem2D;
  KN<R> * diag=GetAny<KN<R>* >((*edia)(stack));
  Matrice_Creuse<R> * sparce_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
  sparce_mat->pUh=0;
  sparce_mat->pVh=0;
  sparce_mat->typemat=TypeSolveMat(TypeSolveMat::GC); //  none square matrice (morse)
  sparce_mat->A.master(new MatriceMorse<R>(diag->N(),*diag));
  return sparce_mat;
}



template<class Rin,class Rout>
 struct  ChangeMatriceMorse {
 static  MatriceMorse<Rout> *f(MatriceMorse<Rin> *mr)
 {
    MatriceMorse<Rout>*  mrr=new MatriceMorse<Rout>(*mr);
    delete mr;
    return mrr;
 }
 };

template<class R> 
 struct  ChangeMatriceMorse<R,R> {
 static MatriceMorse<R>* f(MatriceMorse<R>* mr)
 {
   return mr;
 }
 };
template<class R,class RR>
AnyType CopyMat_tt(Stack stack,Expression emat,Expression eA,bool transp)
{
  using namespace Fem2D;
  Matrice_Creuse<R> * Mat;
 
  if(transp)
   {
    Matrice_Creuse_Transpose<R>  tMat=GetAny<Matrice_Creuse_Transpose<R> >((*eA)(stack));
    Mat=tMat; 
   }
  else   Mat =GetAny<Matrice_Creuse<R>*>((*eA)(stack));
  MatriceMorse<R> * mr=Mat->A->toMatriceMorse(transp,false);
  MatriceMorse<RR> * mrr = ChangeMatriceMorse<R,RR>::f(mr);
  
  Matrice_Creuse<RR> * sparce_mat =GetAny<Matrice_Creuse<RR>* >((*emat)(stack));
  sparce_mat->pUh=Mat->pUh;
  sparce_mat->pVh=Mat->pUh;;
  sparce_mat->typemat=TypeSolveMat(TypeSolveMat::GC); //  none square matrice (morse)
  sparce_mat->A.master(mrr);
  return sparce_mat;
}

template<class R,class RR>
AnyType CopyTrans(Stack stack,Expression emat,Expression eA)
{
 return CopyMat_tt<R,RR>(stack,emat,eA,true);
}
template<class R,class RR>
AnyType CopyMat(Stack stack,Expression emat,Expression eA)
{
 return CopyMat_tt<R,RR>(stack,emat,eA,false);
}

template<class R>
AnyType MatFull2Sparce(Stack stack,Expression emat,Expression eA)
{
  KNM<R> * A=GetAny<KNM<R>* >((*eA)(stack));
  Matrice_Creuse<R> * sparce_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
  sparce_mat->pUh=0;
  sparce_mat->pVh=0;
  sparce_mat->typemat=TypeSolveMat(TypeSolveMat::GMRES); //  none square matrice (morse)
  sparce_mat->A.master(new MatriceMorse<R>(*A,0.0));
  
 return sparce_mat;
}

template<class RR,class AA=RR,class BB=AA> 
struct Op2_pair: public binary_function<AA,BB,RR> { 
  static RR f(const AA & a,const BB & b)  
  { return RR( a, b);} 
};


template<class R>
long get_mat_n(Matrice_Creuse<R> * p)
 { ffassert(p ) ;  return p->A ?p->A->n: 0  ;}

template<class R>
long get_mat_m(Matrice_Creuse<R> * p)
 { ffassert(p ) ;  return p->A ?p->A->m: 0  ;}



template<class R>
long get_diag(Matrice_Creuse<R> * p, KN<R> * x)
 { ffassert(p && x ) ;  return p->A ?p->A->getdiag(*x): 0  ;}
template<class R>
long set_diag(Matrice_Creuse<R> * p, KN<R> * x)
 { ffassert(p && x ) ;  return p->A ?p->A->setdiag(*x): 0  ;}

 
template<class R>  
R * get_elementp2mc(Matrice_Creuse<R> * const  & ac,const long & b,const long & c){ 
   MatriceCreuse<R> * a= ac ? ac->A:0 ;
  if(  !a || a->n <= b || c<0 || a->m <= c  ) 
   { cerr << " Out of bound  0 <=" << b << " < "  << a->n << ",  0 <= " << c << " < "  << a->m
           << " array type = " << typeid(ac).name() << endl;
     cerr << ac << " " << a << endl;
     ExecError("Out of bound in operator Matrice_Creuse<R> (,)");}
    return  &((*a)(b,c));}


template<class RR,class AA=RR,class BB=AA> 
struct Op2_mulAv: public binary_function<AA,BB,RR> { 
  static RR f(const AA & a,const BB & b)  
  { return (*a->A * *b );} 
};

template<class RR,class AA=RR,class BB=AA> 
struct Op2_mulvirtAv: public binary_function<AA,BB,RR> { 
  static RR f(const AA & a,const BB & b)  
  { return RR( (*a).A, *b );} 
};
 
 
template<class K>
class OneBinaryOperatorA_inv : public OneOperator { public:  
  OneBinaryOperatorA_inv() : OneOperator(atype<Matrice_Creuse_inv<K> >(),atype<Matrice_Creuse<K> *>(),atype<long>()) {}
    E_F0 * code(const basicAC_F0 & args) const 
     { Expression p=args[1];
       if ( ! p->EvaluableWithOutStack() ) 
        { 
          bool bb=p->EvaluableWithOutStack();
          cout << bb << " " <<  * p <<  endl;
          CompileError(" A^p, The p must be a constant == -1, sorry");}
       long pv = GetAny<long>((*p)(0));
        if (pv !=-1)   
         { char buf[100];
           sprintf(buf," A^%d, The pow must be  == -1, sorry",pv);
           CompileError(buf);}     
       return  new E_F_F0<Matrice_Creuse_inv<K>,Matrice_Creuse<K> *>(Build<Matrice_Creuse_inv<K>,Matrice_Creuse<K> *>,t[0]->CastTo(args[0])); 
    }
};

template<class K>
class Psor :  public E_F0 { public:  
 
   typedef double  Result;
   Expression mat;
   Expression xx,gmn,gmx,oomega;  
   Psor(const basicAC_F0 & args) 
    {   
      args.SetNameParam();
      mat=to<Matrice_Creuse<K> *>(args[0]); 
      gmn=to<KN<K>*>(args[1]); 
      gmx=to<KN<K>*>(args[2]); 
      xx=to<KN<K>*>(args[3]); 
      oomega=to<double>(args[4]); 
      
   }   
    static ArrayOfaType  typeargs() { 
      return  ArrayOfaType( atype<double>(),
                            atype<Matrice_Creuse<K> *>(),
                            atype<KN<K>*>(),
                            atype<KN<K>*>(),
                            atype<KN<K>*>(),
                            atype<double>(),false);}
                            
    static  E_F0 * f(const basicAC_F0 & args){ return new Psor(args);} 
    
    AnyType operator()(Stack s) const {
      Matrice_Creuse<K>* A= GetAny<Matrice_Creuse<K>* >( (*mat)(s) );
      KN<K>* gmin = GetAny<KN<K>* >( (*gmn)(s) );
      KN<K>* gmax = GetAny<KN<K>* >( (*gmx)(s) );
      KN<K>* x = GetAny<KN<K>* >( (*xx)(s) );
      double omega = GetAny<double>((*oomega)(s));
      return A->A->psor(*gmin,*gmax,*x,omega);
    }
  
};
template <class R>
 struct TheDiagMat {
  Matrice_Creuse<R> * A; 
  TheDiagMat(Matrice_Creuse<R> * AA) :A(AA) {ffassert(A);}
  void   get_mat_daig( KN_<R> & x) { ffassert(A && A->A && x.N() == A->A->n  && A->A->n == A->A->m );
     A->A->getdiag(x);}
  void   set_mat_daig(const  KN_<R> & x) { ffassert(A && A->A && x.N() == A->A->n  && A->A->n == A->A->m );
     A->A->setdiag(x);}
 };
 
template<class R>
TheDiagMat<R> thediag(Matrice_Creuse<R> * p)
 {  return  TheDiagMat<R>(p);}
 
template<class R>
TheDiagMat<R> set_mat_daig(TheDiagMat<R> dm,KN<R> * x)
{
  dm.set_mat_daig(*x);
  return dm;
}
template<class R>
KN<R> * get_mat_daig(KN<R> * x,TheDiagMat<R> dm)
{
  dm.get_mat_daig(*x);
  return x;
}
 
template <class R>
void AddSparceMat()
{
 Dcl_Type<TheDiagMat<R> >();

TheOperators->Add("*", 
        new OneBinaryOperator<Op2_mulvirtAv<typename VirtualMatrice<R>::plusAx,Matrice_Creuse<R>*,KN<R>* > >,
        new OneBinaryOperator<Op2_mulvirtAv<typename VirtualMatrice<R>::plusAtx,Matrice_Creuse_Transpose<R>,KN<R>* > >,
        new OneBinaryOperator<Op2_mulvirtAv<typename VirtualMatrice<R>::solveAxeqb,Matrice_Creuse_inv<R>,KN<R>* > >     
        );
TheOperators->Add("^", new OneBinaryOperatorA_inv<R>());
  
// matrix new code   FH (Houston 2004)        
 TheOperators->Add("=",
//       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation::Op*,E_F_StackF0F0>(SetMatrixInterpolation),
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const Matrix_Prod<R,R>,E_F_StackF0F0>(ProdMat<R,R,R>),
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,KN<R> *,E_F_StackF0F0>(DiagMat<R>),       
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,Matrice_Creuse_Transpose<R>,E_F_StackF0F0>(CopyTrans<R,R>), 
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,Matrice_Creuse<R>*,E_F_StackF0F0>(CopyMat<R,R>) ,
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,KNM<R>*,E_F_StackF0F0>(MatFull2Sparce<R>) ,
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,list<pair<R,MatriceCreuse<R> *> > *,E_F_StackF0F0>(CombMat<R>) 
       );
       
 TheOperators->Add("<-",
//       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation::Op*,E_F_StackF0F0>(SetMatrixInterpolation),
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const Matrix_Prod<R,R>,E_F_StackF0F0>(ProdMat<R,R,R>),
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,KN<R> *,E_F_StackF0F0>(DiagMat<R>)  ,
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,Matrice_Creuse_Transpose<R>,E_F_StackF0F0>(CopyTrans<R,R>), 
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,Matrice_Creuse<R>*,E_F_StackF0F0>(CopyMat<R,R>) ,
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,KNM<R>*,E_F_StackF0F0>(MatFull2Sparce<R>) ,
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,list<pair<R,MatriceCreuse<R> *> > *,E_F_StackF0F0>(CombMat<R>) 
       
       );
TheOperators->Add("*", 
        new OneBinaryOperator<Op2_pair<Matrix_Prod<R,R>,Matrice_Creuse<R>*,Matrice_Creuse<R>*> >,
        new OneBinaryOperator<Op2_pair<Matrix_Prod<R,R>,Matrice_Creuse_Transpose<R>,Matrice_Creuse<R>* > >, 
        new OneBinaryOperator<Op2_pair<Matrix_Prod<R,R>,Matrice_Creuse_Transpose<R>,Matrice_Creuse_Transpose<R> > >,
        new OneBinaryOperator<Op2_pair<Matrix_Prod<R,R>,Matrice_Creuse<R>*,Matrice_Creuse_Transpose<R> > > ,
        new OneBinaryOperator<Op2_ListCM<R> >  , 
        new OneBinaryOperator<Op2_ListMC<R> >   
        );
TheOperators->Add("+", 
        new OneBinaryOperator<Op2_ListCMCMadd<R> >,
        new OneBinaryOperator<Op2_ListCMMadd<R> >,
        new OneBinaryOperator<Op2_ListMCMadd<R> >,
        new OneBinaryOperator<Op2_ListMMadd<R> >
       
       ); 
        
 Add<Matrice_Creuse<R> *>("n",".",new OneOperator1<long,Matrice_Creuse<R> *>(get_mat_n<R>) );
 Add<Matrice_Creuse<R> *>("m",".",new OneOperator1<long,Matrice_Creuse<R> *>(get_mat_m<R>) );
 
 Add<Matrice_Creuse<R> *>("diag",".",new OneOperator1<TheDiagMat<R> ,Matrice_Creuse<R> *>(thediag<R>) );
// Add<Matrice_Creuse<R> *>("setdiag",".",new OneOperator2<long,Matrice_Creuse<R> *,KN<R> *>(set_diag<R>) );
 TheOperators->Add("=", new OneOperator2<KN<R>*,KN<R>*,TheDiagMat<R> >(get_mat_daig<R>) );
 TheOperators->Add("=", new OneOperator2<TheDiagMat<R>,TheDiagMat<R>,KN<R>*>(set_mat_daig<R>) );
 
 Global.Add("set","(",new SetMatrix<R>);
 //Global.Add("psor","(",new  OneOperatorCode<Psor<R> > );
 
 atype<Matrice_Creuse<R> * >()->Add("(","",new OneOperator3_<R*,Matrice_Creuse<R> *,long,long >(get_elementp2mc<R>));



      
//  --- end  
}

//extern Map_type_of_map map_type_of_map ; //  to store te type 
//extern Map_type_of_map map_pair_of_type ; //  to store te type 
extern int lineno(); 
class  PrintErrorCompile : public OneOperator {
    public: 
    const char * cmm;
    E_F0 * code(const basicAC_F0 & ) const 
     { ErrorCompile(cmm,lineno());
      return 0;} 
    PrintErrorCompile(const char * cc): OneOperator(map_type[typeid(R).name()]),cmm(cc){}

};

class PrintErrorCompileIM :  public E_F0info { public:  
 typedef double  Result;
 static E_F0 *   f(const basicAC_F0 & args)  
    {   
     lgerror("\n\n *** change interplotematrix in interpole.\n  *** Bad name in previous version,\n *** sorry FH.\n\n");
     return 0;  }   
    static ArrayOfaType  typeargs() {return  ArrayOfaType(true);}
    operator aType () const { return atype<double>();} 

};
void  init_lgmat() 

{
   map_type_of_map[make_pair(atype<Matrice_Creuse<double>* >(),atype<double*>())]=atype<Matrice_Creuse<double> *>();
   map_type_of_map[make_pair(atype<Matrice_Creuse<double>* >(),atype<Complex*>())]=atype<Matrice_Creuse<Complex> *>();
    AddSparceMat<double>();
    AddSparceMat<Complex>();
 
 Add<const MatrixInterpolation::Op *>("<-","(", new MatrixInterpolation);
 Global.Add("interpolate","(",new MatrixInterpolation);
 Global.Add("interplotematrix","(",new  OneOperatorCode<PrintErrorCompileIM>);
       

 // pour compatibiliter 

  TheOperators->Add("=",
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation::Op*,E_F_StackF0F0>(SetMatrixInterpolation));
       
 TheOperators->Add("<-",
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation::Op*,E_F_StackF0F0>(SetMatrixInterpolation));
 // construction of complex matrix form a double matrix
 TheOperators->Add("=", new OneOperator2_<Matrice_Creuse<Complex>*,Matrice_Creuse<Complex>*,Matrice_Creuse<double>*,E_F_StackF0F0>(CopyMat<R,Complex>)
                 );
                    
 TheOperators->Add("<-", new OneOperator2_<Matrice_Creuse<Complex>*,Matrice_Creuse<Complex>*,Matrice_Creuse<double>*,E_F_StackF0F0>(CopyMat<R,Complex>)
                 );
}

