#include  <iostream>
using namespace std;

#include "rgraph.hpp"
#include "error.hpp"
#include "AFunction.hpp"

//#include "lex.hpp"
#include "MatriceCreuse_tpl.hpp"
#include "MeshPoint.hpp"
#include "lgfem.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"
#include <set>

basicAC_F0::name_and_type  CDomainOfIntegration::name_param[]= {
     "qft", &typeid(const Fem2D::QuadratureFormular *),
     "qfe", &typeid(const Fem2D::QuadratureFormular1d *),
     "qforder",&typeid(long),
     "qfnbpT",&typeid(long),
     "qfnbpE",&typeid(long),
     "optimize",&typeid(bool),
     "binside",&typeid(double)
};


basicAC_F0::name_and_type  Problem::name_param[]= {
     "init", &typeid(bool),
     "solver", &typeid(TypeSolveMat*),
     "eps", &typeid(double)  ,
     "precon",&typeid(Polymorphic*), 
     "dimKrylov",&typeid(long),
     "bmat",&typeid(Matrice_Creuse<R>* ),
     "tgv",&typeid(double ),
     "strategy",&typeid(long )
     
};
basicAC_F0::name_and_type  OpCall_FormBilinear<C_args>::name_param[]= {
     "init", &typeid(bool),
     "solver", &typeid(TypeSolveMat*),
     "eps", &typeid(double)  ,
     "precon",&typeid(Polymorphic*), 
     "dimKrylov",&typeid(long),
     "bmat",&typeid(Matrice_Creuse<R>* ),
     "tgv",&typeid(double ),
     "factorize",&typeid(bool),
     "strategy",&typeid(long )
          
};
basicAC_F0::name_and_type  OpCall_FormBilinear<FormBilinear>::name_param[]= {
     "init", &typeid(bool),
     "solver", &typeid(TypeSolveMat*),
     "eps", &typeid(double)  ,
     "precon",&typeid(Polymorphic*), 
     "dimKrylov",&typeid(long),
     "bmat",&typeid(Matrice_Creuse<R>* ),
     "tgv",&typeid(double ),
     "factorize",&typeid(bool),
     "strategy",&typeid(long )     
};
basicAC_F0::name_and_type  OpCall_FormLinear<FormLinear>::name_param[]= {
     "tgv",&typeid(double )
};
basicAC_F0::name_and_type  OpCall_FormLinear2<FormBilinear>::name_param[]= {
     "tgv",&typeid(double )
};
basicAC_F0::name_and_type  OpCall_FormLinear<C_args>::name_param[]= {
     "tgv",&typeid(double )
};
basicAC_F0::name_and_type  OpCall_FormLinear2<C_args>::name_param[]= {
     "tgv",&typeid(double )
};

namespace Fem2D {

void Check(const Opera &Op,int N,int  M)
 {
   int err=0;
   for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++)
     {  // attention la fonction test donne la ligne 
       //  et la fonction test est en second      
       BilinearOperator::K ll(*l);
       pair<int,int> jj(ll.first.first),ii(ll.first.second);
       if (ii.first <0 || ii.first >= M) err++;
       if (jj.first <0 || jj.first >= N) err++;
       
     }
   if (err) {
     cout << "Check Bilinear Operator" << N << " " << M << endl;
     for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++)
       {  // attention la fonction test donne la ligne 
         //  et la fonction test est en second      
         BilinearOperator::K ll(*l);
         pair<int,int> jj(ll.first.first),ii(ll.first.second);
         cout << " +  " << jj.first << " " << jj.second << "*" << ii.first << " " << ii.second << endl;
          }
     ExecError("Check BilinearOperator N M");
   }
 }
 void Check(const  BC_set<R> * bc,int N)
 {
   int err=0;
   int kk=bc->bc.size();
   for (int k=0;k<kk;k++)
     {
       pair<int,Expression> xx=bc->bc[k];
       if (xx.first >= N) { 
         err++;
         cerr << " Sorry : just " << N << " componant in FE space \n"
              << "   and Boundary condition refere to " << xx.first << "+1 componant " << endl;
       }
     }
   if (err) 
     ExecError("Incompatibility beetwen  boundary condition  and FE space");
 }
 
  void Check(const  Ftest * fl,int N)
 {
   assert(fl);
   int err=0;
   Ftest::const_iterator kk= fl->v.end(),k;
   int ii=0;
   for (k=fl->v.begin();k<kk;k++)
     { 
       ii++;
       int j=k->first.first;
       if (  j >= N) { 
         err++;
         cerr << " Sorry : just " << N << " componant in FE space \n"
              << " and linear var form  refere to " << j << "+1 componant (part " << ii << ")" << endl;
       }
     }
   if (err) 
     ExecError("Incompatibility beetwen linear varf  and FE space");
 }
 //---------------------------------------------------------------------------------------
   void  Element_OpVF(MatriceElementairePleine<R> & mat,
                      const FElement & Ku,const FElement & KKu,
                      const FElement & Kv,const FElement & KKv,
   R * p,int ie,int iie, int label,void *bstack)
  {
   
   pair<void *,double *> * bs=static_cast<pair<void *,double *> *>(bstack);   
   void * stack= bs->first;
   double binside = *bs->second; // truc FH pour fluide de grad2 (decentrage bizard)
   assert(mat.onFace); //   Finite Volume or discontinous Galerkine 
    assert(ie>=0 && ie < 3); //  int on edge 
    MeshPoint mp= *MeshPointStack(stack);
    double ** copt = Stack_Ptr<double*>(stack,ElemMatPtrOffset);

  bool same = &Ku == & Kv;
  assert(same); 
  const Triangle & T  = Ku.T;
  int nTonEdge =  &Ku == &KKu ? 1 : 2; 
  double cmean = 1./nTonEdge;

  throwassert(&T == &Kv.T);  
  const QuadratureFormular & FI = mat.FIT;
  const QuadratureFormular1d & FIb = mat.FIE;
  long npi;
  R *a=mat.a;
  R *pa=a;
  long i,j;
  long n= mat.n,m=mat.m,nx=n*m;
  assert(nx<=mat.lga);
  long N= Kv.N;
  long M= Ku.N;
  
  long mu=Ku.NbDoF();
  long mmu=KKu.NbDoF();
  long nv=Kv.NbDoF();
  long nnv=Kv.NbDoF();
  assert(mu==mmu && nv == nnv) ; 
  
  
  
  const Opera &Op(*mat.bilinearform);
  bool classoptm = copt && Op.optiexpK;
//  if (Ku.number<1 && verbosity/100 && verbosity % 10 == 2) 
  if (Ku.number<1 && ( verbosity > 1 ) )
     cout << "Element_OpVF P: copt = " << copt << " " << classoptm << " binside (For FH) =" << binside <<endl;
  
    
  KN<bool> Dop(last_operatortype); //  sinon ca plate bizarre 
  Op.DiffOp(Dop);  
  int lastop=1+Dop.last(binder1st<equal_to<bool> >(equal_to<bool>(),true));
 //assert(lastop<=3);
  int lffv = nv*N*last_operatortype;  
  int lffu = mu*M*last_operatortype;
  int loffset =  same ? 0 :  (nv+nnv)*N*last_operatortype;
  
  RNMK_ fv(p,nv,N,lastop); //  the value for basic fonction in K
  RNMK_ ffv(p + lffv ,nnv,N,lastop); //  the value for basic fonction in KK
  RNMK_ fu(  (R*) fv   + loffset  ,mu,M,lastop); //  the value for basic fonction
  RNMK_ ffu( (R*) fu  + lffu  ,mmu,M,lastop); //  the value for basic fonction
  
  R2 E=T.Edge(ie);
  R le = sqrt((E,E));
  R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
           PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]),
           PC(TriangleHat[OppositeVertex[ie]]);
  // warning the to edge are in opposite sens         
  R2 PPA(TriangleHat[VerticesOfTriangularEdge[iie][1]]),
     PPB(TriangleHat[VerticesOfTriangularEdge[iie][0]]),
     PPC(TriangleHat[OppositeVertex[ie]]);
  R2 Normal(E.perp()/-le); 
  for (npi=0;npi<FIb.n;npi++) // loop on the integration point
      {
        pa =a;
        QuadratureFormular1d::Point pi( FIb[npi]);        
        R coef = le*pi.a;
        R sa=pi.x,sb=1-sa;
        R2 Pt(PA*sa+PB*sb ); //
        R2 PPt(PPA*sa+PPB*sb ); //  
        if (binside) {
            Pt   = (1-binside)*Pt + binside*PC; 
            PPt  = (1-binside)*PPt + binside*PPC; }
        Ku.BF(Dop,Pt,fu);
        KKu.BF(Dop,PPt,ffu);
        if (!same) { Kv.BF(Dop,Pt,fv); KKv.BF(Dop,PPt,ffv); }     
        // int label=-999999; // a passer en argument 
        MeshPointStack(stack)->set(T(Pt),Pt,Kv,label, Normal,ie);
        if (classoptm) (*Op.optiexpK)(stack); // call optim version 
        
        
        for ( i=0;  i<n;   i++ )  
            { 
              int ik= mat.nik[i];
              int ikk=mat.nikk[i]; 
               
              RNM_ wi(fv(Max(ik,0),'.','.'));   
              RNM_ wwi(ffv(Max(ikk,0),'.','.'));   
                  
              for ( j=0;  j<m;   j++,pa++ ) 
                { 
                  int jk= mat.njk[j];
                  int jkk=mat.njkk[j];
                    
                  RNM_ wj(fu(Max(jk,0),'.','.'));
                  RNM_ wwj(ffu(Max(jkk,0),'.','.'));
                  
                  int il=0;
                  for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                      {       
                        BilinearOperator::K ll(*l);
                        pair<int,int> jj(ll.first.first),ii(ll.first.second);
                        int iis = ii.second, jjs=jj.second;
                        
                        int iicase  = iis / last_operatortype;
                        int jjcase  = jjs / last_operatortype;
                        
                         iis %= last_operatortype;
                         jjs %= last_operatortype;
                        R w_i=0,w_j=0,ww_i=0,ww_j=0;
                        
                        if(ik>=0) w_i =   wi(ii.first,iis ); 
                        if(jk>=0) w_j =   wj(jj.first,jjs );
                        
                        if( iicase>0 && ikk>=0) ww_i =  wwi(ii.first,iis ); 
                        if( jjcase>0 && jkk>=0) ww_j =  wwj(jj.first,jjs );
                       
                        
                        if       (iicase==Code_Jump) w_i = ww_i-w_i; // jump
                        else  if (iicase==Code_Mean) { 
                            
                            w_i = cmean*  (w_i + ww_i );} // average
                        else  if (iicase==Code_OtherSide) w_i = ww_i;  // valeur de autre cote
                        
                        if      (jjcase==Code_Jump) w_j = ww_j-w_j; // jump
                        else if (jjcase==Code_Mean) w_j = cmean*  (w_j +ww_j ); // average
                        else if (jjcase==Code_OtherSide) w_j = ww_j;  //  valeur de l'autre cote    
                        
                       // R ccc = GetAny<R>(ll.second.eval(stack));
                       
                        R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                       if ( copt && Kv.number <1)
                        {
                         double cc  =  GetAny<double>(ll.second.eval(stack));
                         if ( ccc != cc) { 
                          cerr << cc << " != " << ccc << " => ";
                         cerr << "Sorry error in Optimization (b) add:  int2d(Th,optimize=0)(...)" << endl;
                         ExecError("In Optimized version "); }
                 }
                         *pa += coef * ccc * w_i*w_j;
                      }
                }
            } 
         // else pa += m;
      }
  
  
     pa=a;
     if (verbosity > 55 && (Ku.number <=0 || KKu.number <=0 )) { 
      cout <<endl  << " edge between " << Ku.number << " , " <<  KKu.number   << " =  "<<  T[0] << ", " << T[1] << ", " << T[2] << " " << nx << endl;
      cout << " K u, uu =  " << Ku.number << " " << KKu.number << " " <<  " K v, vv =  " << Kv.number << " " << KKv.number << " " <<endl; 
      for (int i=0;i<n;i++)
      {
      cout << setw(2) << i << setw(4) << mat.ni[i] <<  setw(4) << mat.nik[i] << setw(4) << mat.nikk[i]  <<  " :";
      for (int j=0;j<m;j++)
      cout << setw(5)  << (*pa++) << " ";
      cout << endl;
      } } 
  
  *MeshPointStack(stack) = mp;
  }  

 //--------------------------------------------------------------------------------------
 
 
 
 
 
void  Element_Op(MatriceElementairePleine<R> & mat,const FElement & Ku,const FElement & Kv,R * p,int ie,int label,void *stack)
  {
    MeshPoint mp= *MeshPointStack(stack);
    double ** copt = Stack_Ptr<double*>(stack,ElemMatPtrOffset);

  bool same = &Ku == & Kv;
  const Triangle & T  = Ku.T;
  throwassert(&T == &Kv.T);  
  const QuadratureFormular & FI = mat.FIT;
  const QuadratureFormular1d & FIb = mat.FIE;
  long npi;
  R *a=mat.a;
  R *pa=a;
  long i,j;
  long n= mat.n,m=mat.m,nx=n*m;
  long N= Kv.N;
  long M= Ku.N;
  
  
  
  
  
  const Opera &Op(*mat.bilinearform);
  bool classoptm = copt && Op.optiexpK;
  if (Ku.number<1 && verbosity/100 && verbosity % 10 == 2) 
     cout << "Element_Op P: copt = " << copt << " " << classoptm << endl;
    assert(Op.MaxOp() <last_operatortype);
//
  
    
  KN<bool> Dop(last_operatortype);
  Op.DiffOp(Dop);  
  int lastop=1+Dop.last(binder1st<equal_to<bool> >(equal_to<bool>(),true));
 //assert(lastop<=3);

  RNMK_ fv(p,n,N,lastop); //  the value for basic fonction
  RNMK_ fu(p+ (same ?0:n*N*lastop) ,m,M,lastop); //  the value for basic fonction
  
  for (i=0;i< nx;i++) 
    *pa++ = 0.; 
  if (ie<0)    
    for (npi=0;npi<FI.n;npi++) // loop on the integration point
      {
        QuadraturePoint pi(FI[npi]);
        R coef = T.area*pi.a;
        R2 Pt(pi);
        pa =a;
        Ku.BF(Dop,Pt,fu);
        MeshPointStack(stack)->set(T(Pt),Pt,Kv);
        if (classoptm) (*Op.optiexpK)(stack); // call optim version 
        
        if (!same) Kv.BF(Dop,Pt,fv);      
        for ( i=0;  i<n;   i++ )  
          { 
            
            // attention la fonction test donne la ligne 
            //  et la fonction test est en second      
            
            RNM_ wi(fv(i,'.','.'));         
            for ( j=0;  j<m;   j++,pa++ ) 
              { 
                RNM_ wj(fu(j,'.','.'));
                int il=0;
                for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                  {  // attention la fonction test donne la ligne 
                    //  et la fonction test est en second      
                    BilinearOperator::K ll(*l);
                    pair<int,int> jj(ll.first.first),ii(ll.first.second);
                    R w_i =  wi(ii.first,ii.second); 
                    R w_j =  wj(jj.first,jj.second);
                    R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                if ( copt && Kv.number <1)
                 {
                     double cc  =  GetAny<double>(ll.second.eval(stack));
                     //cout << *(copt[il]) << " == " <<  cc << endl;
                     if ( ccc != cc) { 
                        cerr << cc << " != " << ccc << " => ";
                       cerr << "Sorry error in Optimization (a) add:  int2d(Th,optimize=0)(...)" << endl;
                       ExecError("In Optimized version "); }
                 }
                    
                    
                    *pa += coef * ccc * w_i*w_j;
                  }
              }
          }
      }
  else // int on edge ie 
    for (npi=0;npi<FIb.n;npi++) // loop on the integration point
      {
        pa =a;
        QuadratureFormular1d::Point pi( FIb[npi]);
        R2 E=T.Edge(ie);
        R le = sqrt((E,E));
        R coef = le*pi.a;
        R sa=pi.x,sb=1-sa;
        R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
          PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]);
        R2 Pt(PA*sa+PB*sb ); //  
        Ku.BF(Dop,Pt,fu);
        if (!same) Kv.BF(Dop,Pt,fv);      
        // int label=-999999; // a passer en argument 
        MeshPointStack(stack)->set(T(Pt),Pt,Kv,label,R2(E.y,-E.x)/le,ie);
        if (classoptm) (*Op.optiexpK)(stack); // call optim version 
        
        
        for ( i=0;  i<n;   i++ )  
         // if (onWhatIsEdge[ie][Kv.DFOnWhat(i)]) // juste the df on edge bofbof generaly wrong FH dec 2003
            { 
              RNM_ wi(fv(i,'.','.'));       
              for ( j=0;  j<m;   j++,pa++ ) 
                { 
                  RNM_ wj(fu(j,'.','.'));
                  int il=0;
                  for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                   // if (onWhatIsEdge[ie][Kv.DFOnWhat(j)]) // juste the df on edge bofbof generaly wrong FH dec 2003
                      {       
                        BilinearOperator::K ll(*l);
                        pair<int,int> jj(ll.first.first),ii(ll.first.second);
                        
                        R w_i =  wi(ii.first,ii.second); 
                        R w_j =  wj(jj.first,jj.second);
                       // R ccc = GetAny<R>(ll.second.eval(stack));
                       
                    R ccc = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                if ( copt && Kv.number <1)
                 {
                     double cc  =  GetAny<double>(ll.second.eval(stack));
                     if ( ccc != cc) { 
                        cerr << cc << " != " << ccc << " => ";
                       cerr << "Sorry error in Optimization (b) add:  int2d(Th,optimize=0)(...)" << endl;
                       ExecError("In Optimized version "); }
                 }
                         *pa += coef * ccc * w_i*w_j;
                      }
                }
            } 
         // else pa += m;  FH dec 2003
      }
  
  
  /*  pa=a;
      if (Ku.Vh.Th(T) >=0 ) {
      cout <<endl  << " Triangle " << Ku.Vh.Th(T) << " =  "<<  T[0] << ", " << T[1] << ", " << T[2] << " " << nx << endl;
      for (int i=0;i<n;i++)
      {
      cout << setw(2) << i << setw(4) << mat.ni[i] << " :";
      for (int j=0;j<m;j++)
      cout << setw(5)  << (*pa++) << " ";
      cout << endl;
      } } 
  */ 
  *MeshPointStack(stack) = mp;
  }  
  
  
  
  void  Element_Op(MatriceElementaireSymetrique<R> & mat,const FElement & Ku,R * p,int ie,int label, void * stack)
  {
    MeshPoint mp= *MeshPointStack(stack);
    double ** copt = Stack_Ptr<double*>(stack,ElemMatPtrOffset);
    const Triangle & T  = Ku.T;
    //  const QuadratureFormular & FI = QuadratureFormular_T_2;
    //  const QuadratureFormular1d & FIb = QF_GaussLegendre2;
    const QuadratureFormular & FI = mat.FIT;
    const QuadratureFormular1d & FIb = mat.FIE;
    long npi;
    R *a=mat.a;
    R *pa=a;
    long i,j;
    long n= mat.n,m=mat.m,nx=n*(m+1)/2;
    long N= Ku.N;
    long M=N;
    bool show = Ku.Vh.Th(T)==0;
    char * xxx[] ={" u"," v"," p"," q"," r"};
    char * xxxx[] ={" u'"," v'"," p'"," q'"," r'"};
    char * yyy[] ={" ","_x ","_y "};
    
    
    throwassert(mat.bilinearform);
    
    const Opera &Op(*mat.bilinearform);
    bool classoptm = copt && Op.optiexpK;
    assert(  (copt !=0) ||  (Op.where_in_stack_opt.size() !=0) );
  if (Ku.number<1  && verbosity/100 && verbosity % 10 == 2 ) 
     cout << "Element_Op S: copt = " << copt << " " << classoptm << endl;
    assert(Op.MaxOp() <last_operatortype);


    KN<bool> Dop(last_operatortype);
    Op.DiffOp(Dop);  
    int lastop=1+Dop.last(binder1st<equal_to<bool> >(equal_to<bool>(),true));
   // assert(lastop<=3);

    RNMK_ fu(p,n,N,lastop); //  the value for basic fonction
    
    pa =a;
    for (i=0;i< nx;i++) 
      *pa++ = 0.; 
    
    if (ie<0)   
      for (npi=0;npi<FI.n;npi++) // loop on the integration point
        {
          QuadraturePoint pi(FI[npi]);
          R coef = T.area*pi.a;
          R2 Pt(pi);
          pa =a;
          Ku.BF(Dop,Pt,fu);
          MeshPointStack(stack)->set(T(pi),pi,Ku);
          if (classoptm) (*Op.optiexpK)(stack); // call optim version 
          for ( i=0;  i<n;   i++ )  
            { 
              RNM_ wi(fu(i,'.','.'));
              //    if (Ku.Vh.Th(T) < 1) cout << i <<" " <<Pt<< "wi =" << wi ;
              for ( j=0;  j<=i;  j++,pa++ ) // 
                {
                  
                  RNM_ wj(fu(j,'.','.'));
                  //   if (Ku.Vh.Th(T) < 1) cout << j <<" " <<Pt<< "wj =" << wj ;
                  int il=0;
                  for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                    {       
                      const  BilinearOperator::K & ll(*l);
                      pair<int,int> ii(ll.first.first),jj(ll.first.second);
                      R w_i =  wi(ii.first,ii.second);
                      R w_j =  wj(jj.first,jj.second);
                      
                      R c = copt ? *(copt[il]): GetAny<double>(ll.second.eval(stack));
                if ( copt && Ku.number <1)
                 {
                     double cc  =  GetAny<double>(ll.second.eval(stack));
                    // cout << *(copt[il]) << " == " <<  cc << endl;
                     if ( c != cc) { 
                       cerr << c << " != " << cc << " => ";
                       cerr << "Sorry error in Optimization (c) add:  int2d(Th,optimize=0)(...)" << endl;
                       ExecError("In Optimized version "); }
                 }
                      
                      *pa += coef * c * w_i*w_j;
                      /*    if (Ku.Vh.Th(T) < 1 && npi < 1 && i < 1 && j < 1 ) 
                            cout <<" + " << c << " (" <<coef << " " << w_i << " " << w_j << " " << jj.first << " " << jj.second << ") " 
                            << xxx[ii.first]  << yyy[ii.second] << " "
                            << xxxx[jj.first] << yyy[jj.second] << " " ;*/
                    }
                }
              
            }
          
        } 
    else // int on edge ie 
      for (npi=0;npi<FIb.n;npi++) // loop on the integration point
        {
          
          pa =a;
          QuadratureFormular1d::Point pi( FIb[npi]);
          R2 E=T.Edge(ie);
          R le = sqrt((E,E));
          R coef = le*pi.a;
          R sa=pi.x,sb=1-sa;
          R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
            PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]);
          R2 Pt(PA*sa+PB*sb ); //  
          Ku.BF(Dop,Pt,fu);
          // int label=-999999; // a passer en argument 
          MeshPointStack(stack)->set(T(Pt),Pt,Ku,label,R2(E.y,-E.x)/le,ie);
          if (classoptm) (*Op.optiexpK)(stack); // call optim version 
          
          
          for ( i=0;  i<n;   i++ )  
           // if ( onWhatIsEdge[ie][Ku.DFOnWhat(i)]) // generaly wrong FH dec 2003
              { 
                RNM_ wi(fu(i,'.','.'));     
                for ( j=0;  j<=i;   j++,pa++ ) 
                  { 
                    RNM_ wj(fu(j,'.','.'));
                    int il=0;
                    for (BilinearOperator::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                     // if (onWhatIsEdge[ie][Ku.DFOnWhat(j)]) // generaly wrong FH dec 2003
                        {       
                          BilinearOperator::K ll(*l);
                          pair<int,int> ii(ll.first.first),jj(ll.first.second);
                          R w_i =  wi(ii.first,ii.second); 
                          R w_j =  wj(jj.first,jj.second);
                         // R ccc = GetAny<R>(ll.second.eval(stack));
                      R ccc = copt ? *(copt[il]): GetAny<double>(ll.second.eval(stack));
                if ( copt && Ku.number <1)
                 {
                     double cc  =  GetAny<double>(ll.second.eval(stack));
                     if ( ccc != cc) { 
                       cerr << ccc << " != " << cc << ", xy = "<< T(Pt) << " => ";
                       cerr << "Sorry error in Optimization (d)  add:  int2d(Th,optimize=0)(...)" << endl;
                       ExecError("In Optimized version "); }
                 }
                          
                          *pa += coef * ccc * w_i*w_j;
                        }
                  }
              } //else pa+= i+1;
        }
    
    /*   
         pa=a;
         if (Ku.Vh.Th(T) <=0 ) {
         cout <<endl  << " Triangle " << Ku.Vh.Th(T) << " =  "<<  T[0] << ", " << T[1] << ", " << T[2] << " " << nx << endl;
         for (int i=0;i<n;i++)
         {
         cout << setw(2) << i << setw(4) << mat.ni[i] << " :";
         for (int j=0;j<=i;j++)
         cout << setw(5)  << (*pa++) << " ";
         cout << endl;
         } } 
         pa=a;
         for (int i=0;i<n;i++)
         cout << mat.ni[i] << " " ;
         for (int i=0;i<n;i++)
         for (int j=0;j<n;j++,pa++)
         if ( mat.ni[i]==150 && mat.nj[j] == 150)
         cout << "a_150,150 = "<< *pa ;
         cout << endl; 
    */
    
    *MeshPointStack(stack) = mp;
    
  }    
  
  
  
  // #pragma optimization_level 0
  void  Element_rhs(const FElement & Kv,const LOperaD &Op,R * p,void * stack,RN_ & B,
                    const QuadratureFormular & FI = QuadratureFormular_T_2)
  {
    MeshPoint mp=*MeshPointStack(stack) ;
    double ** copt = Stack_Ptr<double*>(stack,ElemMatPtrOffset);
    const Triangle & T  = Kv.T;
    //  const QuadratureFormular & FI = QuadratureFormular_T_2;
    //  const QuadratureFormular & FI = QuadratureFormular_T_2;
    long npi;
    long i,n=Kv.NbDoF(),N=Kv.N;
    
    //  bool show = Kv.Vh.Th(T)==0;
    char * xxx[] ={" u"," v,"," p"," q"," r"};
    char * xxxx[] ={" u'"," v',"," p'"," q'"," r'"};
    char * yyy[] ={" ","_x ","_y "};
    
    bool classoptm = copt && Op.optiexpK;
    assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
    if (Kv.number<1  && verbosity/100 && verbosity % 10 == 2) 
     cout << "Element_rhs S0: copt = " << copt << " " << classoptm << endl;


    KN<bool> Dop(last_operatortype);
    Op.DiffOp(Dop);  
    int lastop=1+Dop.last(binder1st<equal_to<bool> >(equal_to<bool>(),true));
    assert(Op.MaxOp() <last_operatortype);
    
  //  assert(lastop<=3);

    
    RNMK_ fu(p,n,N,lastop); //  the value for basic fonction
    
    for (npi=0;npi<FI.n;npi++) // loop on the integration point
      {
        QuadraturePoint pi(FI[npi]);
        R coef = T.area*pi.a;
        R2 Pt(pi);
        Kv.BF(Dop,Pt,fu);
        MeshPointStack(stack)->set(T(Pt),Pt,Kv);
        if (classoptm) (*Op.optiexpK)(stack); // call optim version 
        for ( i=0;  i<n;   i++ )  
          { 
            RNM_ wi(fu(i,'.','.'));
            int il=0;
            for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
              {       
                LOperaD::K ll(*l);
                pair<int,int> ii(ll.first);
                R w_i =  wi(ii.first,ii.second);
                //copt=0;
                R c = copt ? *(copt[il]) : GetAny<double>(ll.second.eval(stack)); //GetAny<double>(ll.second.eval(stack));
                if ( copt && Kv.number <1)
                 {
                     double cc  =  GetAny<double>(ll.second.eval(stack));
                     if ( c != cc) { 
                       cerr << c << " != " << cc << " => ";
                       cerr << "Sorry error in Optimization add:  int2d(Th,optimize=0)(...)" << endl;
                       ExecError("In Optimized version "); }
                 }
                //if (Kv.number<5) cout << il<< " " << i << "  c== " <<  c << endl;
                R a = coef * c * w_i;
                B[Kv(i)] += a;
              }
          }
        
        
      }  
    *MeshPointStack(stack) = mp;
    
    
  }  
  
  // #pragma optimization_level 0
  void  Element_rhs(const  Mesh & ThI,const Triangle & KI,
                    const FESpace & Vh,const LOperaD &Op,R * p,void * stack,RN_ & B,
                    const QuadratureFormular & FI = QuadratureFormular_T_2)
  {
    MeshPoint mp=*MeshPointStack(stack) ;
    double ** copt = Stack_Ptr<double*>(stack,ElemMatPtrOffset);
//    int maxd = Op.MaxOp();
//    assert(maxd<last_operatortype);
    const Triangle * Kp=0;

    bool classoptm = copt && Op.optiexpK;
    assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
    if (ThI(KI)<1 && verbosity/100 && verbosity % 10 == 2)

     cout << "Element_rhs 3: copt = " << copt << " " << classoptm << endl;

    KN<bool> Dop(last_operatortype);
    Op.DiffOp(Dop);  
    int lastop=1+Dop.last(binder1st<equal_to<bool> >(equal_to<bool>(),true));
    assert(Op.MaxOp() <last_operatortype);

   // assert(lastop<=3);
      
    for (long npi=0;npi<FI.n;npi++) // loop on the integration point
      {
        QuadraturePoint pi(FI[npi]);
        R2 PI(KI(pi));      
        R coef = KI.area*pi.a;
        MeshPointStack(stack)->set(ThI,PI,pi,KI,KI.lab);
        if (classoptm) (*Op.optiexpK)(stack); // call optim version 
        bool outside;
        R2 Pt;
        const Triangle & K  = *Vh.Th.Find(PI,Pt,outside,Kp);
        if ( ! outside) 
          {
            const  FElement  Kv= Vh[K];
            long i,n=Kv.NbDoF(),N=Kv.N;
            RNMK_ fu(p,n,N,lastop); //  the value for basic fonction
            Kv.BF(Dop,Pt,fu);
            
            for ( i=0;  i<n;   i++ )  
              { 
                RNM_ wi(fu(i,'.','.'));
                int il=0;
                for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                  {       
                    LOperaD::K ll(*l);
                    pair<int,int> ii(ll.first);
                    
                    R w_i =  wi(ii.first,ii.second);
                    
                    R c = copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));;//GetAny<double>(ll.second.eval(stack));
                if ( copt && ThI(KI) <1)
                 {
                     double cc  =  GetAny<double>(ll.second.eval(stack));
                     if ( c != cc) { 
                       cerr << c << " != " << cc << " => ";
                       cerr << "Sorry error in Optimization add:  int2d(Th,optimize=0)(...)" << endl;
                       ExecError("In Optimized version "); }
                 }
                    
                    R a = coef * c * w_i;
                    B[Kv(i)] += a;
                  }
              }
          }
        Kp = & K; 
      }  
    *MeshPointStack(stack) = mp;
    
    
  }  
  
  
  void  Element_rhs(const FElement & Kv,int ie,int label,const LOperaD &Op,R * p,void * stack,RN_ & B,
                    const QuadratureFormular1d & FI = QF_GaussLegendre2,bool alledges=false)
  {
    MeshPoint mp=*MeshPointStack(stack) ;
    double ** copt = Stack_Ptr<double*>(stack,ElemMatPtrOffset);
    const Triangle & T  = Kv.T;
    // const QuadratureFormular1d & FI = QF_GaussLegendre2;
    long npi;
    long i,n=Kv.NbDoF(),N=Kv.N;
    
    //  bool show = Kv.Vh.Th(T)==0;
    char * xxx[] ={" u"," v,"," p"," q"," r"};
    char * xxxx[] ={" u'"," v',"," p'"," q'"," r'"};
    char * yyy[] ={" ","_x ","_y "};

    bool classoptm = copt && Op.optiexpK;
    assert(  (copt !=0) ==  (Op.where_in_stack_opt.size() !=0) );
    if (Kv.number<1 && verbosity/100 && verbosity % 10 == 2) 
     cout << "Element_rhs S: copt = " << copt << " " << classoptm << endl;
    KN<bool> Dop(last_operatortype);
    Op.DiffOp(Dop);  
    int lastop=1+Dop.last(binder1st<equal_to<bool> >(equal_to<bool>(),true));
    assert(Op.MaxOp() <last_operatortype);
   // assert(lastop<=3);

    RNMK_ fu(p,n,N,lastop); //  the value for basic fonction
    
    for (npi=0;npi<FI.n;npi++) // loop on the integration point
      {
        QuadratureFormular1d::Point pi( FI[npi]);
        R2 E=T.Edge(ie);
        R le = sqrt((E,E));
        R coef = le*pi.a;
        R sa=pi.x,sb=1-sa;
        R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
          PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]);
        R2 Pt(PA*sa+PB*sb ); //  
        Kv.BF(Dop,Pt,fu);
        MeshPointStack(stack)->set(T(Pt),Pt,Kv,label,R2(E.y,-E.x)/le,ie);
        if (classoptm) (*Op.optiexpK)(stack); // call optim version         
        
        for ( i=0;  i<n;   i++ )  
         // if (alledges || onWhatIsEdge[ie][Kv.DFOnWhat(i)]) // bofbof faux si il y a des derives ..
            { 
              RNM_ wi(fu(i,'.','.'));
              int il=0;
              for (LOperaD::const_iterator l=Op.v.begin();l!=Op.v.end();l++,il++)
                {       
                  LOperaD::K ll(*l);
                  pair<int,int> ii(ll.first);
                  R w_i =  wi(ii.first,ii.second);
                  R c =copt ? *(copt[il]) : GetAny<R>(ll.second.eval(stack));
                if ( copt && Kv.number<1 <1)
                 {
                     double cc  =  GetAny<double>(ll.second.eval(stack));
                     if ( c != cc) { 
                       cerr << "Sorry orrer in Optimization add:  int2d(Th,optimize=0)(...)" << endl;
                       ExecError("In Optimized version "); }
                 }
                  
                  
                  //= GetAny<double>(ll.second.eval(stack));
                 
                  B[Kv(i)] += coef * c * w_i;
                }
            }
        
        
      }  
    *MeshPointStack(stack) = mp;
    
  } 
  
  
  
  void AssembleBilinearForm(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                            MatriceCreuse<R>  & A, const  FormBilinear * b  )
    
  {
    const CDomainOfIntegration & di= *b->di;
    const Mesh * pThdi = GetAny<pmesh>( (* di.Th)(stack));
    if ( pThdi != &Th) { 
      ExecError("No way to compute bilinear form with integrale of on mesh \n"
                "  test  or unkwon function  defined on an other mesh! sorry to hard.   ");
    }
    SHOWVERB(cout << " FormBilinear " << endl);
    MatriceElementaireSymetrique<R> *mates =0;
    MatriceElementairePleine<R> *matep =0;
    const bool useopt=di.UseOpt(stack);    
    double binside=di.binside(stack);

    const vector<Expression>  & what(di.what);             
    CDomainOfIntegration::typeofkind  kind = di.kind;
    set<int> setoflab;
    bool all=true; 
    const QuadratureFormular1d & FIE = di.FIE(stack);
    const QuadratureFormular & FIT = di.FIT(stack);
    bool VF=b->VF();  // finite Volume or discontinous Galerkin
    if (verbosity>2) cout << "  -- discontinous Galerkin  =" << VF << " size of Mat =" << A.size()<< " Bytes\n";
    if (verbosity>3) 
      if (CDomainOfIntegration::int1d==kind) cout << "  -- boundary int border  " ;
      else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all edges, "   ;
      else cout << "  --  int  in  " ;
    for (int i=0;i<what.size();i++)
      {long  lab  = GetAny<long>( (*what[i])(stack));
      setoflab.insert(lab);
      if ( verbosity>3) cout << lab << " ";
      all=false;
      }
     if (verbosity>3) cout <<" Optimized = "<< useopt << ", ";
  const E_F0 & optiexp0=*b->b->optiexp0;
  const E_F0 & optiexpK=*b->b->optiexpK;
  int n_where_in_stack_opt=b->b->where_in_stack_opt.size();
  double** where_in_stack =0;
  if (n_where_in_stack_opt && useopt)
    where_in_stack = new double * [n_where_in_stack_opt];
  if (where_in_stack)
   {
    assert(b->b->v.size()==n_where_in_stack_opt);
    for (int i=0;i<n_where_in_stack_opt;i++)
    {
      int offset=b->b->where_in_stack_opt[i];
      assert(offset>10);
      where_in_stack[i]= static_cast<double *>(static_cast<void *>((char*)stack+offset));
      *(where_in_stack[i])=0;
     }
    if(&optiexp0) optiexp0(stack);
    if (verbosity/100 && verbosity % 10 == 2)
     {
    int il=0;
    for (BilinearOperator::const_iterator l=b->b->v.begin();l!=b->b->v.end();l++,il++)
       cout << il << " coef (" << l->first << ") = " << *(where_in_stack[il]) << " offset=" << b->b->where_in_stack_opt[il] <<endl;
    }
    }
    Stack_Ptr<double*>(stack,ElemMatPtrOffset) =where_in_stack;
    void *paramate=stack;
    pair<void *,double *> parammatElement_OpVF;  
    parammatElement_OpVF.first = stack;
    parammatElement_OpVF.second= & binside;
    
    if (verbosity >3) 
      if (all) cout << " all " << endl ;
      else cout << endl;
    if(VF) {
      if(&Uh != &Vh || sym)
       ExecError("To Day in bilinear form with discontinous Galerkin:   \n"
                "  test or unkwon function must be  defined on the same FEspace, \n"
                "  and the matrice is not symetric. \n" 
                " To do other case in a future (F. Hecht) dec. 2003 ");
      
      matep= new MatriceElementairePleine<R>(Uh,VF,FIT,FIE);
      matep->faceelement = Element_OpVF;   
      paramate= &parammatElement_OpVF;            
    }
    else if (sym) {
      mates= new MatriceElementaireSymetrique<R>(Uh,FIT,FIE);
      mates->element = Element_Op;               
    }
    else {
      matep= new MatriceElementairePleine<R>(Uh,Vh,FIT,FIE);
      matep->element = Element_Op;               
    }
    MatriceElementaire<R> & mate(*( sym? (MatriceElementaire<R> *)mates : (MatriceElementaire<R> *) matep));
    
    
    mate.bilinearform=b->b;
    
    Check(*mate.bilinearform,mate.Uh.N,mate.Vh.N);
    
    if (di.kind == CDomainOfIntegration::int1d )
      {
        for( int e=0;e<Th.neb;e++)
          {
            if (all || setoflab.find(Th.bedges[e].lab) != setoflab.end())   
              {                  
                int ie,i =Th.BoundaryTriangle(e,ie);
                A += mate(i,ie,Th.bedges[e].lab,stack);  
              }
          }
      }
    else if (di.kind == CDomainOfIntegration::intalledges)
      {
        for (int i=0;i< Th.nt; i++) 
          {
            if ( all || setoflab.find(Th[i].lab) != setoflab.end())
             for (int ie=0;ie<3;ie++)   
              A += mate(i,ie,Th[i].lab,paramate);   
          }
         
      }      
    else if (di.kind == CDomainOfIntegration::int2d )
      {
        for (int i=0;i< Th.nt; i++) 
          {
            if ( all || setoflab.find(Th[i].lab) != setoflab.end())  
              A += mate(i,-1,Th[i].lab,stack);   
            // AA += mate;
          }
      } 
    else 
      InternalError(" kind of CDomainOfIntegration unkown");
      
    if (where_in_stack) delete [] where_in_stack;
    delete &mate;
  }
  
  bool AssembleVarForm(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                       MatriceCreuse<R>  * A,KN<R> * B,const list<C_F0> &largs)
  { // return true if BC 
    bool ret=false; 
    typedef DotStar_KN_<R> DotStar;
    typedef DotSlash_KN_<R> DotSlash;
    list<C_F0>::const_iterator ii,ib=largs.begin(),
      ie=largs.end();
    aType tFB( atype<const  FormBilinear *>() );                       
    aType tMat( atype<Matrice_Creuse<R>*>() );                       
    aType tFL( atype<const  FormLinear *>() );                       
    aType tTab( atype<KN<R> *>() );                       
    aType tMatX( atype<VirtualMatrice<R>::plusAx >() );  
    aType tMatTX( atype<VirtualMatrice<R>::plusAtx >() );  
    aType tDotStar(atype<DotStar >() );
    aType tBC( atype<const  BC_set<R>  *>()) ;                    
    for (ii=ib;ii != ie;ii++)
      {
        Expression e=ii->LeftValue();
        aType r = ii->left();
        if (r==tFB) 
          { if (A)
            AssembleBilinearForm( stack,Th,Uh,Vh,sym,*A,dynamic_cast<const  FormBilinear *>(e));
          
          }
        else if (r==tMat)
          {
            if (A)
              InternalError(" Add sparce matrice; to do, sorry");
          }
        else if (r==tFL)
          {
            if (B)
              AssembleLinearForm( stack,Th, Vh, B,dynamic_cast<const  FormLinear *>(e)) ;
          }
        else if (r==tTab)
          {
            if ( B) 
              *B += *GetAny<KN<R> *>( (*e)(stack) );
          }
        else if (r==tDotStar)
          {
            if ( B) 
              {
                DotStar ab=GetAny<DotStar>( (*e)(stack) );
                *B += ab;
              }
          }
        else if (r==tMatX)
          {
            if ( B) 
              { 
                *B += GetAny<VirtualMatrice<R>::plusAx >( (*e)(stack) )  ;
              }
          }
        else if (r==tMatTX)
          {
            if ( B) 
              { 
                *B += GetAny<VirtualMatrice<R>::plusAtx >( (*e)(stack) )  ;
              }
          }
        else if (r== tBC) 
          ret=true;
        else 
          { 
            cerr << "Problem:operator() unkwon type " << * r <<  endl;
            throw(ErrorExec("Problem:operator() unkwon type",1));
          }
      }
    return ret;
  }                            
  
  void AssembleBC(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                  MatriceCreuse<R>  * A,KN<R> * B,KN<R> * X, const list<C_F0> &largs , R tgv  )
  {
    list<C_F0>::const_iterator ii,ib=largs.begin(),
      ie=largs.end();
    aType tBC( atype<const  BC_set<R>  *>()) ;                    
    for (ii=ib;ii != ie;ii++)
      {
        Expression e=ii->LeftValue();
        aType r = ii->left();
        if (r==tBC)
          AssembleBC(stack,Th,Uh,Vh,sym,A,B,X, dynamic_cast<const  BC_set<R> *>(e),tgv);
      }
    
  }
  void AssembleBC(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                  MatriceCreuse<R>  * A,KN<R> * B,KN<R> * X, const  BC_set<R> * bc, R tgv  )
    
  {
    MeshPoint *mps= MeshPointStack(stack),mp=*mps;
    
    int ktbc=0, nbon =0;
    bool Aii = A && A->n == A->m;
    int Nbcomp=Vh.N;
    Check(bc,Nbcomp);
    assert(Vh.N == Uh.N);
    TabFuncArg tabexp(stack,Vh.N);
    RN buf(Vh.MaximalNbOfDF()*3*Vh.N),gg(buf);
    if ( B && B->N() != Vh.NbOfDF) ExecError("AssembleBC size rhs and nb of DF of Vh");
    if(verbosity>99) cout << " Problem : BC_set " ;
    nbon =bc->on.size();
    set<long> on;
    for (int i=0;i<nbon;i++)
      {
        long  lab  = GetAny<long>( (*bc->on[i])(stack));
        if(verbosity>99) cout << lab << " " ;
        on.insert(lab);
      }
    if(verbosity>99) 
      cout << endl;
    int kk=bc->bc.size();
    
    const int dim=Vh.N;
    FElement::aIPJ ipj(Vh[0].Pi_h_ipj()); 
    FElement::aR2  PtHat(Vh[0].Pi_h_R2()); 
    
    KN<int> PtonB(PtHat.N());
    
    KN<R>   Aipj(ipj.N());
    KNM<R>  Vp(dim,PtHat.N());
    
    
    for (int ib=0;ib<Th.neb;ib++)
      {
        int ie;
        int it = Th.BoundaryTriangle(ib,ie);
        int r =Th.bedges[ib].lab;
        if (on.find(r) != on.end() ) 
          {
            const FElement K(Uh[it]);
            R2 E=K.T.Edge(ie);
            R le = sqrt((E,E));
            
            ktbc++;
            if(verbosity>99)   cout << "BC " << it << " " << ie << " lab=" << r <<  ":\t"
                                    << K.T[VerticesOfTriangularEdge[ie][0]] << "; " 
                                    << K.T[VerticesOfTriangularEdge[ie][1]] << " E=" << K.T.Edge(ie) << endl;
            
            for (int k=0;k<kk;k++)
              {
                gg=0;
                pair<int,Expression> xx=bc->bc[k];
                tabexp=0; 
                int comp = xx.first;
                tabexp[comp]=xx.second;
                // while  (comp+1 <Nbcomp && which_uh[comp+1] == which_uh[comp])   
                while  (comp+1 <Nbcomp && Uh.dim_which_sub_fem[comp+1] == Uh.dim_which_sub_fem[comp])   
                  {  // the right 
                    k++; // NEXT COMP 
                    comp++;
                    if (k<kk && (comp == bc->bc[k].first) )
                      tabexp[comp]=bc->bc[k].second;
                    else
                      CompileError("In Boundary condition the vector FESpace , we must have:"
                                   " all componant, in the right order");
                    
                  }
#ifdef OLDPih                    
                K.Pi_h(gg,F_Pi_h,buf,&tabexp);
                
#else
         K.Pi_h(Aipj);
         PtonB = 0;
        for (int i=0;i<Aipj.N();i++)
           PtonB[ipj[i].p] += onWhatIsEdge[ie][K.DFOnWhat(ipj[i].i)] ;
                       
         for (int p=0;p<PtHat.N();p++)
          if (PtonB[p]) // in on boundary 
          { 
            mps->set(K.T(PtHat[p]),PtHat[p],K,r,R2(E.y,-E.x)/le,ie); // la normal bofbof ?
            KN_<double> Vpp(Vp('.',p));
            for (int j=0;j<dim;j++)
             if (tabexp[j]) 
               Vpp[j]=GetAny<double>( (*tabexp[j])(stack) );
              else Vpp[j]=0;
           }
           
         for (int i=0;i<Aipj.N();i++)
          { 
           const FElement::IPJ &ipj_i(ipj[i]);
           gg[ipj_i.i] += Aipj[i]*Vp(ipj_i.j,ipj_i.p);           
          }
#endif
                int nbdf = K.NbDoF();
                for (int df=0;df<nbdf;df++)
                  // if (K.FromFE(df)==which_uh[xx.first] && onWhatIsEdge[ie][K.DFOnWhat(df)] ) 
                  if (K.FromFE(df)==Uh.dim_which_sub_fem[xx.first] && onWhatIsEdge[ie][K.DFOnWhat(df)] ) 
                    {
                      // cout << k << " df=" << df <<  " g= " << gg[df] <<" " << gg(FromTo(0,2)) << endl;
                      int ddf=K(df);
                      // AA(ddf,ddf) =tgv;
                      if (Aii)  A->diag(ddf)=tgv;
                      if (B) (*B)[ddf]=tgv*gg[df]; 
                      if (X) (*X)[ddf]=gg[df];
                    }
              }
          }
      }
    if (! ktbc  && nbon && verbosity ) 
      {
        cout << " Warning: -- Your set of boundary condition is incompatible with the mesh label." << endl;
      }
    *mps =mp;            
  }
  
  void AssembleLinearForm(Stack stack,const Mesh & Th,const FESpace & Vh,KN<double> * B,const  FormLinear * l )
  {
    Check(l->l,Vh.N);
    if ( B && B->N() != Vh.NbOfDF) ExecError("AssembleLinearForm size rhs and nb of DF of Vh");
   
    KN<double> buf(Vh.MaximalNbOfDF()*3*Vh.N);
    
    //            const  FormLinear * l=dynamic_cast<const  FormLinear *>(e);
    const CDomainOfIntegration & di= *l->di;
    const Mesh & ThI = * GetAny<pmesh>( (* di.Th)(stack));
    bool sameMesh = &ThI == &Th;
    
    SHOWVERB(cout << " FormLinear " << endl);
    const vector<Expression>  & what(di.what);
    
    CDomainOfIntegration::typeofkind  kind = di.kind;
    const QuadratureFormular1d & FIE = di.FIE(stack);
    const QuadratureFormular & FIT = di.FIT(stack);
    const bool useopt=di.UseOpt(stack);    
    
  //  cout << "AssembleLinearForm " << l->l->v.size() << endl; 
    set<int> setoflab;
    bool all=true; 
    if ( verbosity>3) 
      if (kind==CDomainOfIntegration::int1d) cout << "  -- boundary int border " ;
      else if (kind==CDomainOfIntegration::intalledges) cout << "  -- boundary int all edges " ;
      else cout << "  -- boundary int  " ;
    for (int i=0;i<what.size();i++)
      {long  lab  = GetAny<long>( (*what[i])(stack));
      setoflab.insert(lab);
      if ( verbosity>3) cout << lab << " ";
      all=false;
      }
     if (verbosity>3) cout << " Optimized = "<< useopt << ", ";
      
  const E_F0 & optiexp0=*l->l->optiexp0;
  const E_F0 & optiexpK=*l->l->optiexpK;
  int n_where_in_stack_opt=l->l->where_in_stack_opt.size();
  double** where_in_stack =0;
  if (n_where_in_stack_opt && useopt)
    where_in_stack = new double * [n_where_in_stack_opt];
  if (where_in_stack)
   {
    assert(l->l->v.size()==n_where_in_stack_opt);
    for (int i=0;i<n_where_in_stack_opt;i++)
    {
      int offset=l->l->where_in_stack_opt[i];
      assert(offset>10);
      where_in_stack[i]= static_cast<double *>(static_cast<void *>((char*)stack+offset));
      *(where_in_stack[i])=0;
     }
     if(&optiexp0) optiexp0(stack);
    
    if( (verbosity/100) && verbosity % 10 == 2)
    {
        int il=0;

        for (LinearOperatorD::const_iterator ll=l->l->v.begin();ll!=l->l->v.end();ll++,il++)
         cout << il << " coef (" << ll->first << ") = " << *(where_in_stack[il]) << " offset=" << l->l->where_in_stack_opt[il] <<endl;
    
        for (int i=0;i<n_where_in_stack_opt;i++)
         cout << "const coef " << i << " = " << *(where_in_stack[i]) << endl;
       }
    }
    Stack_Ptr<double*>(stack,ElemMatPtrOffset) =where_in_stack;
    
    
    if (verbosity >3) 
      if (all) cout << " all " << endl ;
      else cout << endl;
    
    if (kind==CDomainOfIntegration::int1d)
      {
        
        for( int e=0;e<ThI.neb;e++)
          {
            if (all || setoflab.find(ThI.bedges[e].lab) != setoflab.end())   
              {                  
                int ie,i =Th.BoundaryTriangle(e,ie);
                if ( sameMesh) 
                  Element_rhs(Vh[i],ie,Th.bedges[e].lab,*l->l,buf,stack,*B,FIE,false); 
                else 
                  InternalError("To Do") ;
              }
          }
      }
    else if (kind==CDomainOfIntegration::intalledges)
     {
      for (int i=0;i< ThI.nt; i++) 
        if (all || setoflab.find(ThI[i].lab) != setoflab.end()) 
         for (int ie=0;ie<3;ie++)
            if ( sameMesh) 
                Element_rhs(Vh[i],ie,Th[i].lab,*l->l,buf,stack,*B,FIE,true); 
             else 
                InternalError("To Do") ;
     }
    else {
      
      for (int i=0;i< ThI.nt; i++) 
        if (all || setoflab.find(ThI[i].lab) != setoflab.end()) 
          if ( sameMesh ) 
            Element_rhs(Vh[i],*l->l,buf,stack,*B,FIT); 
          else 
            Element_rhs(ThI,ThI[i],Vh,*l->l,buf,stack,*B,FIT);
    }  
    
    if (n_where_in_stack_opt) delete [] where_in_stack;
             
  }
  
  
  
}

bool isVF(const list<C_F0> & largs)  // true => VF type of Matrix   
{
  list<C_F0>::const_iterator ii,ib=largs.begin(),
    ie=largs.end();
    
  bool VVF =false;   
  for (ii=ib;ii != ie;ii++)
    {
      Expression e=ii->LeftValue();
      aType r = ii->left();
      if (r==atype<const  FormBilinear *>()) 
        {
          const  FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
          bool vvf  = bb->VF();
          if( vvf &&  bb->di->kind != CDomainOfIntegration::intalledges )
            CompileError("Sorry, no  jump or moy in bilinear form no of type intalledges  ");
           VVF = vvf || VVF;
          }
    }
  return VVF;
} 
AnyType Problem::operator()(Stack stack) const
{  
  using namespace Fem2D;
  MeshPoint *mps= MeshPointStack(stack),mp=*mps;
  long NbSpace = 50; 
  long itmax=0; 
  double eps=1e-6;
//  bool VF=false;
//  VF=isVF(op->largs);
 // assert(!VF); 
  R tgv = 1e30;
// type de matrice par default
#ifdef UMFPACKxxxx         
     TypeSolveMat tmat(TypeSolveMat::UMFpack); 
#else            
    TypeSolveMat tmat(TypeSolveMat::LU);
#endif    
     
   TypeSolveMat    *typemat=&tmat;
  bool initmat=true;
  int umfpackstrategy=0; 
  if (nargs[0]) initmat= ! GetAny<bool>((*nargs[0])(stack));
  if (nargs[1]) typemat= GetAny<TypeSolveMat *>((*nargs[1])(stack));
  if (nargs[2]) eps= GetAny<double>((*nargs[2])(stack));
  // 3 precon 
  if (nargs[4]) NbSpace= GetAny<long>((*nargs[4])(stack));
  if (nargs[6]) tgv= GetAny<double>((*nargs[6])(stack));
  if (nargs[7]) umfpackstrategy = GetAny<long>((*nargs[7])(stack));
  
  bool sym = typemat->sym;
  
  list<C_F0>::const_iterator ii,ib=op->largs.begin(),
    ie=op->largs.end();
  int Nbcomp2=var.size(),Nbcomp=Nbcomp2/2; // nb de composante 
  throwassert(Nbcomp2==2*Nbcomp);
  Data *data= dataptr(stack);
  //   data->init();
  KN<int>  which_comp(Nbcomp2),which_uh(Nbcomp2);
  
  TabFuncArg tabexp(stack,Nbcomp);
  
  vector< pfer > u_hh(Nbcomp2); 
  for (int i=0;i<var.size();i++)
    u_hh[i] = GetAny< pfer  >( (*(var[i]))(stack));
  for (int i=0;i<var.size();i++)
    u_hh[i].first->newVh(); 
  //   compression pour les cas vectoriel
  int kkk=0;
  for (int i=0;i<Nbcomp2;i++)
    {
      if ( u_hh[i].second==0) 
        kkk++;
      else {          
        throwassert(u_hh[i].second==(u_hh[i-1].second+1));}
      which_uh[i]=kkk-1;  
      which_comp[i]=u_hh[i].second;  
    }
  
  vector< pferbase > u_h(kkk); 
  kkk= 0;
  for (int i=0;i<Nbcomp2;i++)
    if ( u_hh[i].second==0) u_h[kkk++]=u_hh[i].first;
  const int  Nb2 = kkk, Nb=Nb2/2; // nb of FESpace 
  throwassert(Nb2==2*Nb);
  
  const FESpace ** LL = new  const FESpace *[var.size()];
  for (int i=0;i<Nb2;i++)
    LL[i]= (*(u_h[i])).newVh();
  SHOWVERB(cout << "Problem  " << Nb << endl);
  
  //   const de  
  
  //  const FESpace * Uhh , *Vhh;
  const Mesh * pTh= &LL[0]->Th;
  for (int i=0;i<Nb2;i++)
    if ( &LL[0]->Th != pTh)
      ExecError("all the finites elements spaces must be defined on the same mesh in solve");
  if ( pTh != data->pTh ) 
    {
       initmat = true;
       data->pTh=pTh;
       if (Nb==1) 
         { //  cas scalaire
           data->Uh=LL[0];
           data->Vh=LL[1]; }
       else 
         { //  cas vectoriel 
           bool same=true;
           for (int i=0;i<Nb;i++)
             if ( LL[i] != LL[Nb+i] )
               {
                 same = false;
                 break;
               }
           if(!same)
             InternalError("Methode de Galerkine (ˆ faire)");
           else
             {
               
               bool unique=true;
               for (int i=1;i<Nb;i++)
                 if ( LL[0] != LL[i]) 
                   {
                     unique = false;
                     break;
                   }
               if (unique) 
                 data->Uh.master( new FESpace(*LL[0],Nb));
               else 
                 data->Uh.master(new FESpace(LL,Nb));
               data->Vh=data->Uh;
             }
           
         }
    }                     
  const FESpace & Uh(*data->Uh);
  const FESpace & Vh(*data->Vh);
  throwassert(Nbcomp==Uh.N && Nbcomp==Vh.N); 
  RN *B=new RN(Vh.NbOfDF);
  RN *X=B;
  
  *B=0;
  
  bool initx = typemat->t==TypeSolveMat::GC;
  
  const  Mesh & Th(Uh.Th);
  
  if (initx) 
    {
      X=new KN<R>(B->N());
      const FEbase<R> & u_h0 = *(u_h[0]);
      const FESpace  * u_Vh = &*u_h0.Vh ;
      
      if ( u_Vh==0  || &(*(u_h[0])).Vh->Th != &Th )
        {
          *X=0;
          if(verbosity>1)
            cout << "   -- Change of Mesh " << (u_Vh ? & (*(u_h[0])).Vh->Th: 0 ) 
                 << "  " << &Th <<  endl;
        }
      else
        {
          if (Nb==1) 
            *X= * u_h[0]->x();
          else { // dispatch the solution  
            const FElement ** sK= new const FElement * [Nb];  
            RN ** sol= new RN * [Nb];
            for (int i=0;i<Nb;i++) {
              
              sol[i] = (*(u_h[i])).x() ;
            }
            
            for (int it=0;it<Th.nt;it++)
              {
                const FElement K(Uh[it]);
                const int nbdf=K.NbDoF();
                for (int i=0;i<Nb;i++)
                  sK[i]= new FElement( (*LL[i])[it]) ;
                for (int df=0;df< nbdf;df++)
                  {  int kfe=K.FromFE(df);
                  int kdf=K.FromDF(df);
                  const FElement & SK(*sK[kfe]);
                  (*X)[K(df)]= (*sol[kfe])[SK(kdf)] ;
                  }
                for (int i=0;i<Nb;i++)
                  delete sK[i];        
              }
            delete [] sol;
            delete [] sK; 
          }}
    }
  if(verbosity>2) cout << "   Problem(): initmat " << initmat << " VF (discontinuous Galerkin) = " << VF << endl;
  if (initmat) 
   {
    if (typemat->profile) 
      data->A.master(new MatriceProfile<R>(Vh,VF));
    else 
      {
        if ( &Uh == & Vh )
          data->A.master(new MatriceMorse<R>(Vh,sym,VF));
        else 
          data->A.master(new MatriceMorse<R>(Vh,Uh,VF));
      }
      MatriceCreuse<R>  & AA(data->A);
     if(verbosity<1) cout <<  "   -- size of Matrix " << AA.size()<< " Bytes" << " skyline =" <<typemat->profile << endl;
    }
  MatriceCreuse<R>  & A(data->A);
  if  (AssembleVarForm( stack,Th,Uh,Vh,sym, initmat ? &A:0 , B, op->largs)) 
    { 
      *B = - *B; 
      AssembleBC     ( stack,Th,Uh,Vh,sym, initmat ? &A:0 , B, initx ? X:0,  op->largs, tgv );
    }
  else 
    *B = - *B;
  
  if (initmat)
    if (typemat->profile)
      {
        MatriceProfile<R> & AA(dynamic_cast<MatriceProfile<R> &>(A));
        throwassert(&AA);
        switch (typemat->t) {
        case TypeSolveMat::LU       : AA.LU(); break;
        case TypeSolveMat::CROUT    : AA.crout(); break;
        case TypeSolveMat::CHOLESKY : AA.cholesky(); break;
        default:
          cerr << " type resolution " << typemat->t << endl;
          CompileError("type resolution profile inconnue"); break;       
        }
      }
    else 
      {
        if(verbosity>5) cout << " Matrice morse GC Precond diag" << endl;
        MatriceMorse<R> & AA(dynamic_cast<MatriceMorse<R> &>(A));
        throwassert(&AA);
        switch (typemat->t) {
        case    TypeSolveMat::GC:   
          if (precon)
            AA.SetSolverMaster(new SolveGCPrecon<R>(AA,precon,stack,eps));
          else 
            AA.SetSolverMaster(new SolveGCDiag<R>(AA,eps));
          break; 
        case TypeSolveMat::GMRES :
          if (precon)
            AA.SetSolverMaster(new SolveGMRESPrecon<R>(AA,precon,stack,NbSpace,itmax,eps));
          else 
            AA.SetSolverMaster(new SolveGMRESDiag<R>(AA,NbSpace,itmax,eps));
         break;
#ifdef UMFPACK         
        case TypeSolveMat::UMFpack :
            AA.SetSolverMaster(new SolveUMFPack<R>(AA,umfpackstrategy,tgv,eps));
         break;
           
#endif         
        default:
          cerr << " type resolution " << typemat->t << endl;
          CompileError("type resolution inconnue"); break;       
        }
        
      }
  if(verbosity>3) cout << "   B  min " << B->min() << " ,  max = " << B->max() << endl;
  if (verbosity>99)
   {
    cout << " X= " << *X << endl;
    cout << " B= " << *B << endl;
    }
    
  A.Solve(*X,*B);
  if (verbosity>99)
   {
    cout << " X= " << *X << endl;
   }
  
  
  // dispatch the solution 
  if (Nb==1)  {
    *(u_h[0])=X;
    if (X != B ) delete B;  }
  else {
    const FElement ** sK= new const FElement * [Nb];
    
    RN ** sol= new RN * [Nb];
    for (int i=0;i<Nb;i++) {
      sol[i]= new RN( LL[i]->NbOfDF) ;
      *(u_h[i]) = sol[i];
    }
    
    for (int it=0;it<Th.nt;it++)
      {
        const FElement K(Uh[it]);
        const int nbdf=K.NbDoF();
        for (int i=0;i<Nb;i++)
          sK[i]= new FElement( (*LL[i])[it]) ;
        for (int df=0;df< nbdf;df++)
          {  int kfe=K.FromFE(df);
          int kdf=K.FromDF(df);
          const FElement & SK(*sK[kfe]);
          (*sol[kfe])[SK(kdf)] = (*X)[K(df)];
          }
        for (int i=0;i<Nb;i++)
          delete sK[i];
        
      }
    
     delete [] sK;
     delete [] sol;
     if (X != B && X ) delete X; 
     delete B; 
  }
  
  if (verbosity) 
    {cout << " -- Solve : " ; 
    for (int i=0;i<Nb;i++) 
      cout  << "          min " << (u_h[i])->x()->min() << "  max " << (u_h[i])->x()->max() << endl ;
    }
    delete [] LL;
    *mps=mp;
    return SetAny<const Problem *>(this);
}





void GetBilinearParam(const ListOfId &l,basicAC_F0::name_and_type *name_param,int n_name_param,
                      Expression *nargs,int & N,int & M,  vector<Expression> & var )
{
  
  for (int i=0;i<n_name_param;i++)
    nargs[i]=0;
  int nb=l.size(),n=0,nbarray=0;
  ListOfId * array[2];
  for (int i=0;i<nb;i++)
    if (l[i].r == 0 && l[i].re  == 0 && l[i].array == 0) 
      n++;
       else if (l[i].array) array[Min(nbarray++,1)] = l[i].array;
    else 
      {  
        bool ok=false;
        for (int j=0;j<n_name_param;j++)
          if (!strcmp(l[i].id,name_param[j].name))
            {
              ok = !nargs[j];
              nargs[j]= map_type[name_param[j].type->name()]->CastTo(C_F0(l[i].e,l[i].re));
              break;
            }
        if (!ok) 
          {
            cerr << " Error name argument " << l[i].id << " the kown arg : "; 
            for (int k=0;k<n_name_param;k++)
              cerr << name_param[k].name << " ";
            cerr << endl;
                 CompileError("Unkown name argument or two times same name argument ");
          }
      }
  
  if (nbarray)
    { // new version ok
      if(nbarray!=2) 
        CompileError(" Must have 2 array, one for unknow functions, one for test functions");
      N = array[0]->size();
      M = array[1]->size();
      var.resize(N+M);
      for (int k=0,j=0;k<2;k++)
        for  (int i=0;i<array[k]->size();i++)
          { 
            const UnId & idi((*array[k])[i]);
            if (idi.r == 0 && idi.re  == 0 && idi.array==0 ) 
              { C_F0 c=::Find( idi.id);
              var[j++]=CastTo<pfer>(c); }         
            else 
              CompileError(" Just Variable in array parameter ");   
          }   
    }
  else
    { // old version
      assert(n%2==0);
      N=n/2;
      M=N; 
      var.resize(N+M);
      for  (int i=0,j=0;i<l.size();i++)
        if (l[i].r == 0 && l[i].re  == 0 && l[i].array==0 ) 
          {
            C_F0 c=::Find(l[i].id);
            var[j++]=CastTo<pfer>(c);
          } 
      
    }
  
}

 


bool FieldOfForm(const C_args *op)  // true => complex problem 
{
  list<C_F0>::const_iterator ii,ib=op->largs.begin(),
    ie=op->largs.end();
  bool complextype =false;   
  for (ii=ib;ii != ie;ii++)
    {
      Expression e=ii->LeftValue();
      aType r = ii->left();
      if (r==atype<const  FormBilinear *>()) 
        {
          const  FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
          if (! bb->b->mappable(BCastToR)) 
            complextype=true;
          }
      else if (r==atype<const  FormLinear *>())
        {
          const  FormLinear * ll=dynamic_cast<const  FormLinear *>(e);
          if (! ll->l->mappable(BCastToR)) 
            complextype=true;
        }
    }
  for (ii=ib;ii != ie;ii++)
    {
      Expression e=ii->LeftValue();
      aType r = ii->left();
      if (r==atype<const  FormBilinear *>()) 
        {
          const  FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
          Foperator * b=const_cast<  Foperator *>(bb->b);
          
          if (complextype)  b->mapping(&CCastToC);
          else b->mapping(&CCastToR) ;                
        }
      else if (r==atype<const  FormLinear *>())
        {
          const  FormLinear * ll=dynamic_cast<const  FormLinear *>(e);
          Ftest * l= const_cast<Ftest *>(ll->l);
          if (complextype)  l->mapping(&CCastToC) ;
          else l->mapping(&CCastToR) ;                
        }
    } 
  return complextype;
}  


Problem::Problem(const C_args * ca,const ListOfId &l,size_t & top) :
  op(ca),var(l.size()),offset(align8(top)),VF(false) 
{
  SHOWVERB(cout << "Problem : -----------------------------" << top << endl);
  top = offset + sizeof(Data);
  GetBilinearParam(l,name_param,n_name_param,nargs, Nitem,Mitem,var);
  
        
  precon = 0; //  a changer 
  if ( nargs[3])
    {
      const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(nargs[3]);
      assert(op);
      precon = op->Find("(",ArrayOfaType(atype<KN<double>* >(),false));
   }

  VF=isVF(op->largs);   
 // cout << " Problem ) VF = " << VF << endl;
  complextype =  FieldOfForm(op);
  if( complextype) 
    CompileError("Problem just Real problem no complex do day");
  
}

Expression IsFebaseArray(Expression f)
{
  assert(f);
  size_t N=f->nbitem();
  const E_Array * vvi(dynamic_cast<const E_Array *>(f));
  if ( ! vvi) return 0;
  const E_Array & vi(*vvi);
  Expression febase=0;
  for (int i=0;i<N;i++)
    { 
      assert(vi[i].left() == atype<pfer>() );
      const E_FEcomp<R> * comp=dynamic_cast<const E_FEcomp<R> *>( vi[i].LeftValue()) ;   
      if (!(comp && comp->comp == i  && comp->N == N)) return 0; 
      if (!febase) febase = comp->a0;
      else if(comp->a0 != febase) return 0;        
    }   
  return febase;  
}

Call_FormBilinear::Call_FormBilinear(Expression * na,Expression  BB,Expression fi, Expression fj)
  : nargs(na),largs(),N(fi->nbitem()),M(fj->nbitem()), 
     euh(fi), evh(fj)
{
  assert(nargs );
  const C_args * LLL=dynamic_cast<const C_args *>(BB);
  if (!LLL) 
    CompileError("Sorry the variationnal form (varf)  is not a the variationnal form (type const C_args *)");
  largs=LLL->largs;
}

Call_FormLinear::Call_FormLinear(Expression *na,Expression  LL, Expression ft)
  : nargs(na),largs(),N(ft->nbitem()),
  ppfes(ft)//IsFebaseArray(ft)) 
{
  const C_args * LLL=dynamic_cast<const C_args *>(LL);
  if ( !LLL) CompileError("The parameter of a LinearForm must be a array of all componate of FE function");
  largs=LLL->largs;
}
bool C_args::IsLinearOperator() const {
  int n=largs.size();
  aType tRn =atype<KN<R>* >();
  for (const_iterator i=largs.begin(); i != largs.end();i++) 
    {  C_F0  c= *i;Expression e=c; aType r=c.left();
    if (     ( r != atype<const  FormLinear *>() )
             &&  ( r != atype<const  BC_set<R> *>() )
             &&  ( r != atype<VirtualMatrice<R>::plusAx >() )
             &&  ( r != atype<VirtualMatrice<R>::plusAtx >() )
             &&  ( r != tRn) 
             ) return false;
    }
  return true;}
bool C_args::IsBilinearOperator() const {
      int n=largs.size();
      aType tRn =atype<Matrice_Creuse<R>* >();
      for (const_iterator i=largs.begin(); i != largs.end();i++) 
        {  C_F0  c= *i;Expression e=c; aType r=c.left();
        if (     ( r!= atype<const  FormBilinear *>() )
                 &&  ( r != atype<const  BC_set<R> *>() )
                 &&  ( r != tRn) 
                 ) return false;
        }
      return true;}


void SetArgsFormLinear(const ListOfId *lid,int ordre)
{  
  //  the local parameter are 
  //  ordre ==2 => bilinear form  unknown (newU_) and test function (newV_)
  //  ordre ==1 =>   linear form just  test function (newV_)
  // ---------------------
  throwassert(ordre >0 && ordre <=2 && (lid || lid->size()>0 ) );
  const ListOfId & l(*lid);
  int nb=l.size();
  int n=0;
  C_F0 type,init;
  int nbarray=0;
  ListOfId * array[2];
  aType uh=atype<const finconnue*>(),vh=atype<const ftest*>();
  
  for (int i=0;i<nb;i++)        
    if (l[i].r == 0 &&  l[i].re == 0 && l[i].id  ) n++; 
         else if (l[i].array) 
           array[Min(nbarray++,2)] = l[i].array;  
  if (nbarray && n==0)
    {  // 
      
      if(nbarray!=ordre) 
        { cerr << " form " << ordre << " == " << nbarray << " Nb of Array "<<endl;
        CompileError(" Must have 1 or 2 array, one for unknow functions, one for test functions");
        }
      for (int k=0,j=0;k<ordre;k++)
        for  (int i=0,iend=array[k]->size();i<iend;i++)
          { 
            const UnId & idi((*array[k])[i].id);
            if (idi.r == 0 && idi.re  == 0 && idi.array==0 ) 
              {
                if (k==ordre-2)  //  unknow function just in case of bilinear form  
                  currentblock->NewID(uh,idi.id,C_F0(newU_(i),uh));
                else   //  test function 
                  currentblock->NewID(vh,idi.id,C_F0(newV_(i),vh));
              }
            else 
              CompileError(" Just Variable in array parameter "); 
          }  
    }
  else if (nbarray==0)
    {    // a supprimer  to remove   in case of bilinear    
      
      SHOWVERB(cout << "SetArgs:: form  set parameter " << endl);
      throwassert( ordre==1 || n%2==0);
      int nn=ordre==1 ? 0 : n/2; // ordre == 1 => no unknown function just test function
      
      for (int i=0,j=0;i<nb;i++)
        if (l[i].r == 0 && l[i].re  == 0 && l[i].array==0) 
          {          
            SHOWVERB(cout <<"  " <<  l[i].id  << " " << (j<nn) << endl); 
            if (j<nn)
              currentblock->NewID(uh,l[i].id,C_F0(newU_(j%nn),uh));
            else 
              currentblock->NewID(vh,l[i].id,C_F0(newV_(j%nn),vh));
            j++;
          }
    }
  else
    {
       CompileError(" Sorry you mixte formulation with and without array "); 
    }
}

const Fem2D::QuadratureFormular & CDomainOfIntegration::FIT(Stack stack) const 
{
  using namespace Fem2D;
  if (nargs[0]) return  *GetAny<const Fem2D::QuadratureFormular *>((*nargs[0])(stack));
  int exact = 5;
  if (nargs[2]) exact=  GetAny<long>((*nargs[2])(stack))-1;
  if( QuadratureFormular_T_1.exact >= exact ) return QuadratureFormular_T_1;
  if( QuadratureFormular_T_2.exact >= exact ) return QuadratureFormular_T_2;
  if( QuadratureFormular_T_5.exact >= exact ) return QuadratureFormular_T_5;
  cerr << " Ordre of the Integration Formular ordre " << exact+1 << " exact = " << exact << endl;
  ExecError(" We find  no Integration Formular on Triangle for this  order to hight");
  return QuadratureFormular_T_1;
}       
const Fem2D::QuadratureFormular1d & CDomainOfIntegration::FIE(Stack stack) const 
{
  using namespace Fem2D;
  if (nargs[1]) return  *GetAny<const Fem2D::QuadratureFormular1d *>((*nargs[1])(stack));
  int exact = 5;
  if (nargs[2]) exact=  GetAny<long>((*nargs[2])(stack))-1;
  if( 1 >= exact ) return QF_GaussLegendre1;
  if( 3 >= exact ) return QF_GaussLegendre2;
  if( 5 >= exact ) return QF_GaussLegendre3;
  cerr << " Ordre of the Integration Formular on Edge, order = " << exact+1 << " exact = " << exact << endl;
  ExecError(" We find  no Integration Formular on Edge  for this  order to hight");
  return QF_GaussLegendre1;
}

void buildInterpolationMatrix(MatriceMorse<R> * m,const FESpace & Uh,const FESpace & Vh)
{  //  Uh = Vh 
  using namespace Fem2D;
  m->symetrique = false;
  m->dummy=false;
  m->a=0;
  m->lg=0;
  m->cl=0;
  m->nbcoef=0;
  const  Mesh & ThU =Uh.Th; // line 
  const  Mesh & ThV =Vh.Th; // colunm
  bool samemesh =  &Uh.Th == &Vh.Th;  // same Mesh
  int thecolor =0;
  int nbn_u = Uh.NbOfNodes;
  int nbn_v = Vh.NbOfNodes;
  int n= m->n;
  int nbcoef =0;
  
  KN<int> color(ThV.nt);
  KN<int> mark(n);
  mark=0;
  
  int *lg = new int [n+1];
  int * cl = 0;
  double *a=0;
  /*
  for (int i=0;i<Uh.NbOfDF+1;i++)
    lg[i]=-1;
  KN<int> nl(n);
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
  
  assert(NVh==NUh); 
  
  
  int nbp= PtHatU.N(); // 
  int nbc= ipjU.N(); // 
  KN<R2> PV(nbp);   //  the PtHat in ThV mesh
  KN<int> itV(nbp); // the Triangle number
  
  KNM<R> aaa(nbp,nbdfVK); 

//   R v[1000],vf[100];
//   assert(N*3*NbDoF<=1000 && NbDoF <100 );
//   KNMK_<R> fb(v,NbDoF,N,op+1); //  the value for basic fonction
//   KN_<R> fk(vf,NbDoF);
  // for (int i=0;i<NbDoF;i++) // get the local value
   // fk[i] = u[K(i)];
    //  get value of basic function
 //  bool whatd[last_operatortype];
 //  for (int i=0;i<last_operatortype;i++) 
 //    whatd[i]=false;
 //  whatd[op]=true;
 //  FB(whatd,K.Vh.Th,K.T,PHat,fb);  
 //  R r = (fb('.',componante,op),fk);  
 //  return r;

   const R eps = 1.0e-10;
   const int sfb1=Vh0.N*1*Vh0.NbDoF();
   KN<R> kv(sfb1*nbp);
   R * v = kv;
    KN<int> ik(nbp); // the Triangle number
//   KNMK_<R> fb(v+ik[i],VhV0.NbDoF,VhV0.N,1); //  the value for basic fonction

   bool whatd[last_operatortype];
   for (int i=0;i<last_operatortype;i++) 
     whatd[i]=false;
   whatd[op_id]=true; // the value of function
  
  for (int step=0;step<2;step++) 
   {
      
	  for (int it=0;it<ThU.nt;it++)
	    {
	      thecolor++; //  change the current color
	      const Triangle & TU(ThU[it]);
	      FElement KU(Uh[it]);
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
	            }
	       }
	       
	      
	      for (int i=0;i<nbp;i++)
	         { 
	              KNMK_<R> fb(v+i*sfb1,nbdfVK,NV,1); // valeur de fonction de base de Vh 
	              // ou:   fb(idf,j,0) valeur de la j composante de la fonction idf 
	              Vh0.FB(whatd,ThV,ThV[itV[i]],PV[i],fb);  
	          }

	      for (int i=0;i<ipj.N();i++) 
	          { // pour tous le terme 
	           const FElement::IPJ &ipj_i(ipjU[i]);
	           assert(ipj_i.j==0); // car  Vh.N=0
	           int dfu = KU(ipj_i.i); // le numero de df global 
	           int j = ipj_i.j; // la composante
	           int p=ipj_i.p;  //  le points
	           
	            KNMK_<R> fb(v+i*sfb1,nbdfVK,NV,1); 
	            KN_<R> fbj(fb('.',j,1)); 
	            
	           for (int idfv=0;idfv<nbdfVK;idfv++) 
	             if (Abs(fbj[idfv])>eps) 
	              {
	                int dfv=KV(idfv);
	                // le term dfu,dfv existe dans la matrice
	                R c= fbj[idfv]; 
	              }
	                      
	          }
	    
	      //   for (int df=0;df<nbdf;df++)         
	       //     (*y)[K(df)] =  gg[df] ;
	       
	       
	    }
	    if (step==0)
	     {
	       assert(0); // do allocation 
	     }
    }
    */
  m->lg=lg;
  m->cl=cl;
  m->a=a;
  m->nbcoef=nbcoef;
  assert(0); // a faire to do
}

AnyType SetMatrixInterpolation(Stack stack,Expression emat,Expression einter)
{
  using namespace Fem2D;
  
  Matrice_Creuse<R> * sparce_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack));
  const MatrixInterpolation::Op * mi(dynamic_cast<const MatrixInterpolation::Op *>(einter));
  assert(einter);

  pfes * pUh = GetAny< pfes * >((* mi->a)(stack));
  pfes * pVh = GetAny<  pfes * >((* mi->b)(stack));
  FESpace * Uh = **pUh;
  FESpace * Vh = **pVh;
  assert(Vh);
  assert(Uh);
  sparce_mat->pUh=pUh;
  sparce_mat->pVh=pVh;
  sparce_mat->typemat=TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
  sparce_mat->A.master(new MatriceMorse<R>(*Uh,*Vh,buildInterpolationMatrix));

  return sparce_mat;
}
