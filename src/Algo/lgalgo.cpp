#include  <iostream>
#include  <cfloat>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
#include "rgraph.hpp"
#include "RNM.hpp"
#include "MatriceCreuse_tpl.hpp"
#include "MeshPoint.hpp"
#include "lgfem.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"
#include "NRJ.hpp"
#include "RosenBrock.hpp"
#include "LineSearch.hpp"
#include "CubicLS.hpp"
#include  "BrentLS.hpp"
#include "Optima.hpp"
#include "BFGS.hpp"
//#include  "CG.hpp"
#include  "NewtonRaphson.hpp"

//template<class R>
extern Block *currentblock;

typedef double R;

class OptimAlgo : public OneOperator 
{ public:
  typedef KN<R> Kn;
  typedef KN_<R> Kn_;
  typedef R REAL;
  typedef Param<REAL> PARAM;
  typedef KN<REAL> VECT;
  typedef KNM<REAL> MAT;
  typedef VirtualMatrice<REAL> VMAT;

   const int cas;


 

  class E_LCG: public E_F0mps { public:
   const int cas;
   static basicAC_F0::name_and_type name_param[] ;
   static const int n_name_param =3;
   Expression nargs[n_name_param];
   Expression X;
   C_F0 inittheparam,theparam,closetheparam; 
   Expression J,dJ,hJ;
   long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
   R arg(int i,Stack stack,R a) const{ return nargs[i] ? GetAny<R>( (*nargs[i])(stack) ): a;}

class lgNRJ : public NRJ<PARAM,VECT,VMAT,REAL> {
private:
   Stack stack;
   Expression J,dJ,hJ,theparame;

protected:
  void setparam( const Param<R>& x )
    {
     KN<double> *p=GetAny<KN<double> *>( (*theparame)(stack) );
     assert( p->N() == x.N());
     *p =x;
    }
public:
  lgNRJ(Stack s,int n,Expression t,Expression JJ,Expression dJJ,Expression hJJ) 
   : NRJ<PARAM,VECT,VMAT,REAL>(n), stack(s),theparame(t),J(JJ),dJ(dJJ),hJ(hJJ)  
    {  }
  
  REAL Val(const Param<R>& x) {
     setparam(x);
     assert(J);
     return GetAny<R>( (*J)(stack));}
     
  VECT* Gradient(const Param<R>& x) {
     setparam(x);
     
     return dJ ? GetAny<Kn*>( (*dJ)(stack)) :0;}
     
  VMAT* Hessian(const Param<R>& x) {
     setparam(x);
     if (!hJ) return 0;
      Matrice_Creuse<R> * M=  GetAny<Matrice_Creuse<R> *>( (*hJ)(stack));
      assert(M && M->A );
      return (VirtualMatrice<R>*) M->A;}
  
};  
  E_LCG(const basicAC_F0 & args,int cc) :
    cas(cc)
   {
      int nbj= args.size()-1;
      currentblock = new Block(currentblock); // make a new block to 
      X = to<Kn*>(args[nbj]);
      C_F0 X_n(args[nbj],"n");
      //  the expression to init the theparam of all 
      inittheparam = currentblock->NewVar<LocalVariable>("the parameter",atype<KN<R> *>(),X_n);
      theparam = currentblock->Find("the parameter"); //  the expression for the parameter
      args.SetNameParam(n_name_param,name_param,nargs);
      const  Polymorphic * opJ=0;
      const  Polymorphic * opdJ=0;
      const  Polymorphic * ophJ=0;
      if (nbj>0)
      {  opJ=  dynamic_cast<const  Polymorphic *>(args[0].LeftValue());
         assert(opJ); }
      if (nbj>1)
      {   opdJ=  dynamic_cast<const  Polymorphic *>(args[1].LeftValue());
         assert(opdJ); }
      if (nbj>2)   
      {  ophJ=  dynamic_cast<const  Polymorphic *>(args[2].LeftValue());
         assert(ophJ); }
      J=dJ=hJ=0;
    
      J= to<R>(C_F0(opJ,"(",theparam));
      if(opdJ)
        dJ= to<Kn*>(C_F0(opdJ,"(",theparam));
      if(ophJ)
       hJ= to< Matrice_Creuse<R> *>(C_F0(ophJ,"(",theparam));
         
      closetheparam=currentblock->close(currentblock);   // the cleanning block expression 
/*         
      if (nargs[2]) 
      {  const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(nargs[2]);
         assert(op); 
         C = op->Find("(",ArrayOfaType(atype<Kn* >(),false)); }
       else  C =0; */
   }
     
     virtual AnyType operator()(Stack stack)  const {
     
    typedef LineSearch<PARAM,VECT,MAT,VMAT,REAL> LS;
    typedef CubicLineSearch<LS> CLS;
    
     R tol=arg(0,stack,1E-6);
     int itMax=arg(1,stack,100L);
     int itMaxLine=arg(2,stack,100L);
     bool verbose=verbosity>3;
     
     Kn &x = *GetAny<Kn *>((*X)(stack));
     const int n=x.N();
     Kn * para = GetAny<KN<double>*>( inittheparam.eval(stack) ) ; // do allocation 
     
     Param<R> param(x);
     lgNRJ nrj1(stack,n,theparam,J,dJ,hJ);
     KN<double> * delta =0;
     if (!dJ) delta = new KN<double>(n); 
     CLS ls( & nrj1,itMaxLine,delta);
     REAL initialNrj = nrj1.getVal(param);
     Optima<LS> *opt=0;
     if(cas==1) 
        opt=new BFGS<LS>((CLS*)&ls,itMax,tol,verbose);
      else if(cas==2)
         opt=new Newt<LS>((CLS*)&ls,itMax,tol,verbose);
      else
        ErrorExec("lgalgo: Case no available Sorry internal bug",cas);
        


      param = opt->optimizer(param);
      
      REAL finalNrj = nrj1.getVal(param);
      
      
      if(verbosity)
        cout <<endl<<"*** RESULTS SUMMARY ***"<<endl;

     if(verbosity>1) {
     cout <<"  The number of call to NRJ : "<< nrj1.Appel_Val() << endl;
     cout <<"  The number of call to gradient : "<< nrj1.Appel_Grad() << endl;
     cout <<"  The number of call to hessian : "<< nrj1.Appel_Hess() << endl; }
     if(verbosity>2)
      {
     cout <<"Normalized residue : ";
     affiche(opt->allResidue());
      }
     if(verbosity) {
     cout <<"  Initial NRJ value : " << initialNrj << endl;
     cout <<"  Final NRJ value : " << finalNrj << endl;}
     delete opt;
     
     /// cout<<"The final optimized parameters : "<<param;
     x=param;
     if( delta) delete delta; 
     closetheparam.eval(stack); // clean memory 
     
     return SetAny<long>(0);
     }  
    operator aType () const { return atype<long>();}         
     
  };
  
  E_F0 * code(const basicAC_F0 & args) const {
    return new E_LCG(args,cas);}
  OptimAlgo(int c) :   OneOperator(atype<long>(),
                             atype<Polymorphic*>(),
                             atype<KN<R> *>()),cas(c){}
  OptimAlgo(int c,int cc) :   OneOperator(atype<long>(),
                             atype<Polymorphic*>(),
                             atype<Polymorphic*>(),
                             atype<KN<R> *>()),cas(c){}
  OptimAlgo(int c,int cc,int ccc) :   OneOperator(atype<long>(),
                             atype<Polymorphic*>(),
                             atype<Polymorphic*>(),
                             atype<Polymorphic*>(),
                             atype<KN<R> *>()),cas(c){}
   
};
//template<class R>
basicAC_F0::name_and_type  OptimAlgo::E_LCG::name_param[]= {
     "eps", &typeid(double)  ,
     "nbiter",&typeid(long) ,
     "nbiterline",&typeid(long)
};

void init_algo()
{
 Global.Add("BFGS","(",new OptimAlgo(1,1));  //  j + dJ
 Global.Add("Newtow","(",new OptimAlgo(2,2,2));  //  j + dJ
}
