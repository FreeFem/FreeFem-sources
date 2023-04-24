// -*- Mode : c++ -*-
//
// SUMMARY  :
// USAGE    :
// ORG      :
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
//

/*

 This file is part of Freefem++

 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef PROBLEM_HPP_
#define PROBLEM_HPP_

extern Block *currentblock;

template<class K> class Matrice_Creuse;
//template<class K>  using MatriceMap=map<pair<int,int>,K>;
template<class K>  using MatriceMap=HashMatrix<int,K>;

//template<class K> class MatriceCreuse;
namespace  Fem2D {
  template<class K> class SolveGCPrecon;
  template<class K> class SolveGMRESPrecon;
  template<class K> class SolveGMRESDiag;
//  int IsoLineK(double *f,R2 *Q,double eps);

}

#include "P1IsoValue.hpp"


template<class K> class SolveGCDiag;
class Plot;
class v_fes; // [[file:lgfem.hpp::v_fes]]

// real [[file:lgfem.hpp::FEbase]] using [[file:lgfem.hpp::v_fes]]
typedef FEbase<double,v_fes> * pferbase ;		// <<pferbase>>
typedef FEbaseArray<double,v_fes> * pferbasearray ;	// <<pferbasearray>>
typedef pair<pferbase,int> pfer ;			// <<pfer>>
typedef pair<pferbasearray,int> pferarray ;		// <<pferarray>>

// complex [[file:lgfem.hpp::FEbase]] using [[file:lgfem.hpp::v_fes]]
typedef FEbase<Complex,v_fes> * pfecbase ;		// <<pfecbase>>
typedef FEbaseArray<Complex,v_fes> * pfecbasearray ;	// <<pfecbasearray>>
typedef pair<pfecbase,int> pfec ;			// <<pfec>>
typedef pair<pfecbasearray,int> pfecarray ;		// <<pfecarray>>


//typedef pair<pmesh *,int> pmesharray ;

typedef   LinearComb<MGauche,C_F0> Finconnue;
typedef   LinearComb<MDroit,C_F0> Ftest;
typedef   LinearComb<pair<MGauche,MDroit>,C_F0> Foperator;

inline int intOp(const MGauche &i) {return i.second;}
inline int intOp(const MDroit &i) {return i.second;}
inline int intOp(pair<MGauche,MDroit> & p) {return Max(intOp(p.first),intOp(p.second));}

inline void SetOp(KN_<bool> & d,const MGauche &i)
          { d[i.second% last_operatortype]=true;}
inline void SetOp(KN_<bool> & d,const MDroit &i)
          { d[(int) i.second % last_operatortype]=true;}
inline void SetOp(KN_<bool> & d,const pair<MGauche,MDroit> & p)
          {SetOp(d,p.first);SetOp(d,p.second);}

inline unsigned int GetDiffOp(const MGauche &i, int& lastop)
   {int op=(i.second% last_operatortype);
     lastop=max(lastop,op) ;
     return 1<<op;}
inline unsigned int GetDiffOp(const MDroit &i, int& lastop)
   {int op=(i.second% last_operatortype);
     lastop=max(lastop,op) ;
     return 1<<op;}
inline unsigned int GetDiffOp(const  pair<MGauche,MDroit> &p, int& lastop)
{ return GetDiffOp(p.first,lastop)|GetDiffOp(p.second,lastop);}

typedef  const Finconnue  finconnue;
typedef  const Ftest ftest;
typedef  const Foperator foperator;

Expression IsFebaseArray(Expression f);

void SetArgsFormLinear(const ListOfId *lid,int ordre);
/*
inline ostream & operator<<(ostream & f,const  TypeSolveMat & tm)
{
  switch(tm.t) {
   case TypeSolveMat::NONESQUARE:  f << "No Square (Sparse Morse)"; break;
   case TypeSolveMat::LU:  f << "LU (Skyline)"; break;
   case TypeSolveMat::CROUT:  f << "CROUT (Skyline)"; break;
   case TypeSolveMat::CHOLESKY:  f << "CHOLESKY (Skyline)"; break;
   case TypeSolveMat::GC:  f << "CG (Sparse Morse)"; break;
   case TypeSolveMat::GMRES:  f << "GMRES (Sparse Morse)"; break;
   case TypeSolveMat::SparseSolver:  f << "SparseSolver (Sparse Morse)"; break;
   default: f << "Unknown  bug???";
   }
  return f;
}

*/
class C_args: public E_F0mps  {public:
  typedef const C_args *  Result;
  list<C_F0> largs;
  typedef list<C_F0> ::const_iterator const_iterator ;
  // il faut expendre
  C_args() :largs(){}
  C_args(C_F0 c) : largs() { if(!c.Zero() )largs.push_back(c);}
  C_args(  const basicAC_F0 & args) :largs(){
    int n=args.size();
    for (int i=0;i< n;i++)
      {
       if(args[i].Zero()) ; //  skip zero term ...
       else  if (args[i].left() == atype<const C_args *>())
          {
            const C_args * a = dynamic_cast<const C_args *>(args[i].LeftValue());
            if (a == NULL) printf("dynamic_cast error\n");
            for (list<C_F0>::const_iterator i=a->largs.begin();i!=a->largs.end();i++)
              if( ! i->Zero()) // skip Zero term
              largs.push_back(*i);
          }
        else
          largs.push_back(args[i]);
      };}
  static ArrayOfaType  typeargs() { return ArrayOfaType(true);}
  AnyType operator()(Stack ) const  { return SetAny<const C_args *>(this);}
  operator aType () const { return atype<const C_args *>();}

  static  E_F0 * f(const basicAC_F0 & args) { return new C_args(args);}
  bool Zero()  const { return largs.empty();} // BIG WARNING April and wrong functon FH v 3.60 .......
  bool IsLinearOperator() const;
  bool IsBilinearOperator() const;
  bool IsBemBilinearOperator() const;
  bool IsMixedBilinearOperator() const;
};

class C_args_minus: public C_args  {public:
  C_args_minus(  const basicAC_F0 & args) {
    int n=args.size();
    ffassert(n==2);
    if (args[0].left() == atype<const C_args *>())
      {
        const C_args * a = dynamic_cast<const C_args *>(args[0].LeftValue());
        ffassert(a);
        for (list<C_F0>::const_iterator i=a->largs.begin();i!=a->largs.end();i++)
          largs.push_back(*i);
      }
    else
      largs.push_back(args[0]);

    largs.push_back(C_F0(TheOperators,"-",args[1]));
     }

  static ArrayOfaType  typeargs() { return ArrayOfaType(atype<const C_args *>(),true);}
  static  E_F0 * f(const basicAC_F0 & args) { return new C_args_minus(args);}
};

bool isVF(const list<C_F0> & largs);

template<typename F>
class Minus_Form: public  E_F0mps    {public:
  typedef const F * Result;
  static ArrayOfaType  typeargs() { return ArrayOfaType(atype<const F *>());}
  static  E_F0 * f(const basicAC_F0 & args) {
    int n=args.size();
    ffassert(n==1);
    aType tF=atype<Result>();
    ffassert(args[0].left() == tF);
    Result f = dynamic_cast<Result>(args[0].LeftValue());
    ffassert(f);
    // F mf = -*f;
    F * rf=new F(-*f);
    return  rf;
  }
    operator aType () const { return atype<Result>();}


};

//template<class RR=double>
class BC_set : public E_F0mps { public:
   bool  complextype;
  typedef const BC_set* Result;
  vector<Expression> on;
  vector<int> onis;

  vector<pair<int,Expression> > bc; //  n� de l'inconnue+ valeur
  BC_set(  const basicAC_F0 & args)
    :on(args.size()),onis(args.size())
    {
    int n = args.size();
    ffassert(args.named_parameter);
    AC_F0::const_iterator ii=args.named_parameter->begin();
    AC_F0::const_iterator ie=args.named_parameter->end();
    bc.resize(args.named_parameter->size());
    complextype=false;
    for (int kk=0;ii!=ie;kk++,ii++)
     {
      if( ! BCastTo<double>(ii->second))
       complextype = true;
     }
    ii=args.named_parameter->begin();
    for (int kk=0;ii!=ie;kk++,ii++)
      { //
        C_F0 x=Find(ii->first);
        if (x.left() != atype<const finconnue *>())
          CompileError("We expected an unknown  u=... of the problem");
        const finconnue * uu = dynamic_cast<const finconnue *>(x.LeftValue());
        ffassert(uu);
        const MGauche *ui=uu->simple();
        ffassert(ui && ui->second == op_id);
        //if(verbosity>9)
         cout <<"bc=" << kk << " on : " << ii->first << " n " << ui->first <<   " = ? " << endl;
        if (complextype)
        bc[kk]= make_pair(ui->first,CastTo<Complex>(ii->second));
        else
        bc[kk]= make_pair(ui->first,CastTo<double>(ii->second));
        //ii->second;
        cout << " ii->second type=" << ii->second.left() << endl;
      }
      cout<< "args.size()=" << args.size() << endl;
    //  sort bc / num de composante
      cout << "atype< double *>" << atype<double *>() << " %%%% " << endl;
        std::sort(bc.begin(),bc.end());
        //if(verbosity>9)
        for (vector<pair<int,Expression> >::iterator i=bc.begin(); i !=bc.end();++i)
            cout <<"bc:  on " <<  i->first << " " << i->second << endl;

    for (int i=0;i<n;i++)
      if( ! BCastTo<KN_<long> >(args[i]))
         {
	   on[i]=CastTo<long>(args[i]);
	   onis[i]=0;
         }
	 else
	 {
	   on[i]=CastTo<KN_<long> >(args[i]);
	   onis[i]=1;
	 }
  }

  // operator by copy : used for FESpace composite to renum block index of finconnue et ftest.
  // Remark: it used also in other part of Freefem
  BC_set(  const BC_set & bcold ): complextype(bcold.complextype), on(bcold.on), 
                                  onis(bcold.onis), bc(bcold.bc){} 
  int nbtrue(bool *ok) const 
     { 
        int k=0; 
        for (size_t i=0;i<bc.size();i++) 
          if (ok[i]) k++;
         return k;
     }
  BC_set(  const BC_set & bcold, bool * ok ): complextype(bcold.complextype), on(bcold.on), 
                                  onis(bcold.onis){
                                    bc.resize( bcold.nbtrue(ok) );
                                    int kk=0;
                                    for (size_t i=0;i<bcold.bc.size();i++){
                                      if(ok[i]){
                                          bc[kk].first  = bcold.bc[i].first;
                                          bc[kk].second = bcold.bc[i].second;
                                          kk++;
                                      }
                                    }
                                    ffassert( kk== bc.size() );
                                  } 

    template<class K>
     void CastToK()
    {
      aType rr =  complextype ? atype<Complex>() : atype<double>();
      if (rr == atype<Complex>()) complextype= true;
      if(verbosity > 10) cout << " CastToK => " << complextype <<endl;
     for ( vector<pair<int,Expression> >::iterator k=bc.begin();k!=bc.end();k++)
	k->second=CastTo<K>(C_F0(k->second,rr)) ;
    }
/* De
  // ajout modif FH  mai 2007 XXXXXXXXXXXXX....
 void mappingC(C_F0 (*f)(const C_F0 &)) {

      for ( vector<pair<int,Expression> >::iterator k=bc.begin();k!=bc.end();k++)
	  k->second=CastTo<Complex>(C_F0(k->second,rr)) ;}
  // fin ajout
*/
  static ArrayOfaType  typeargs() { return ArrayOfaType(/*atype<long>(),*/true);} //  change frev 2011 FH...
  AnyType operator()(Stack ) const  { return SetAny<Result>(this);}
  operator aType () const { return atype<Result>();}

  static  E_F0 * f(const basicAC_F0 & args) { return new BC_set(args);}
  //    void init(Stack stack) const {}

};

class CDomainOfIntegration: public E_F0mps {
public:
  static const int n_name_param =12;
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  enum typeofkind  { int2d=0, int1d=1, intalledges=2,intallVFedges=3, int3d = 4, intallfaces= 5,intallVFfaces=6,int0d=7,intall0d=8 } ; //3d
  typeofkind  kind; //  0
  int d; // 3d
  int dHat;//
//  bool isMeshS,isMeshL;
  typedef const CDomainOfIntegration* Result;
  Expression Th;
  Expression mapt[3],mapu[3];
  vector<Expression> what;
  vector<int> whatis; // 0 -> long , 1 -> array ???
  CDomainOfIntegration( const basicAC_F0 & args,typeofkind b=int2d,int ddim=2,int ddHat=0) // 3d
    :kind(b),d(ddim),dHat(ddHat==0 ? d : ddHat), Th(0), what(args.size()-1),whatis(args.size()-1)

  {
      
  //  isMeshS=surface;
 //   isMeshL=curve;
    mapt[0]=mapt[1]=mapt[2]=0; // no map of intergration points for test function
    mapu[0]=mapu[1]=mapu[2]=0; // no map of intergration points for unknows function
    args.SetNameParam(n_name_param,name_param,nargs);
    if(d==2) // 3d
      Th=CastTo<pmesh>(args[0]);
    else if(d==3 && dHat==3){
        Th=CastTo<pmesh3>(args[0]);}
    else if(d==3 && dHat==2){
        Th=CastTo<pmeshS>(args[0]);}
    else if(d==3 && dHat==1){
        Th=CastTo<pmeshL>(args[0]);}
    else ffassert(0); // a faire
    int n=args.size();

    for (int i=1;i<n;i++)
      if(!BCastTo<KN_<long> >(args[i]) )
	 {
	   whatis[i-1]=0;
	   what[i-1]=CastTo<long>(args[i]);
	 }
      else
	 {
	   whatis[i-1]=1;
	   what[i-1]=CastTo<KN_<long> >(args[i]);
	 }
      const E_Array *pmapt = dynamic_cast<E_Array *>(nargs[10]);
      const E_Array *pmapu = dynamic_cast<E_Array *>(nargs[11]);
      if(pmapt )
      {
          if( pmapt->size() != d ) ErrorCompile("mapt bad arry size ",1);
          for(int i=0; i<d; ++i)
              mapt[i]=CastTo<double >((*pmapt)[i]);
      }
      if(pmapu )
      {
          if( pmapu->size() != d ) ErrorCompile("mapu bad arry size ",1);
          for(int i=0; i<d; ++i)
              mapu[i]=CastTo<double >((*pmapu)[i]);

      }
    // cout << " CDomainOfIntegration " << this << endl;
  }
  static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmesh>(), true);} // all type
  AnyType operator()(Stack ) const  { return SetAny<const CDomainOfIntegration *>(this);}
  operator aType () const { return atype<const CDomainOfIntegration *>();}

  static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args);}
  const Fem2D::QuadratureFormular & FIT(Stack) const ;
  const Fem2D::QuadratureFormular1d & FIE(Stack) const ;
  const Fem2D::GQuadratureFormular<R3> & FIV(Stack) const ;  // 3d
  long UseOpt(Stack s) const  {  return nargs[5] ? GetAny<long>( (*(nargs[5]))(s) )  : 1;}
  double  binside(Stack s) const { return nargs[6] ? GetAny<double>( (*(nargs[6]))(s) )  : 0;} // truc pour FH
  bool intmortar(Stack s) const { return nargs[7] ? GetAny<bool>( (*(nargs[7])) (s) )  : 1;} // truc  pour
  double levelset(Stack s) const { return nargs[9] ? GetAny<double>( (*(nargs[9]))(s) )  : 0;}
  bool  islevelset() const { return nargs[9] != 0; }
  bool withmap() const {return mapu[0] || mapt[0]; }
};


class CDomainOfIntegrationBorder: public CDomainOfIntegration {
public:
  CDomainOfIntegrationBorder( const basicAC_F0 & args) :CDomainOfIntegration(args,int1d) {}
  static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args,int1d);}
};

class CDomainOfIntegrationAllEdges: public CDomainOfIntegration {
public:
  CDomainOfIntegrationAllEdges( const basicAC_F0 & args) :CDomainOfIntegration(args,intalledges) {}
  static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args,intalledges);}
};

class CDomainOfIntegrationVFEdges: public CDomainOfIntegration {
public:
  CDomainOfIntegrationVFEdges( const basicAC_F0 & args) :CDomainOfIntegration(args,intallVFedges) {}
  static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args,intallVFedges);}
};

// 3D Volume
class CDomainOfIntegration3d: public CDomainOfIntegration {
public:
  CDomainOfIntegration3d( const basicAC_F0 & args) :CDomainOfIntegration(args,int3d,3) {}
  static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args,int3d,3);}
  static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmesh3>(), true);} // all type
};

class CDomainOfIntegrationBorder3d: public CDomainOfIntegration {
public:
  CDomainOfIntegrationBorder3d( const basicAC_F0 & args) :CDomainOfIntegration(args,int2d,3) {}
  static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args,int2d,3);}
  static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmesh3>(), true);} // all type
};

class CDomainOfIntegrationAllFaces: public CDomainOfIntegration {
public:
  CDomainOfIntegrationAllFaces( const basicAC_F0 & args) :CDomainOfIntegration(args,intallfaces,3) {}
  static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args,intallfaces,3);}
  static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmesh3>(), true);} // all type
};
//  MeshL 
class CDomainOfIntegrationAll0d: public CDomainOfIntegration {
public:
    CDomainOfIntegrationAll0d( const basicAC_F0 & args) :CDomainOfIntegration(args,intall0d,3) {}
  static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args,intall0d,3,1);}
  static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmeshL>(), true);} // all type
};

// 3D surface

class CDomainOfIntegrationS: public CDomainOfIntegration {
public:
    CDomainOfIntegrationS( const basicAC_F0 & args) :CDomainOfIntegration(args,int2d,3,2) {}
    static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args,int2d,3,2);}
    static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmeshS>(), true);} // all type
};

class CDomainOfIntegrationBorderS: public CDomainOfIntegration {
public:
    CDomainOfIntegrationBorderS( const basicAC_F0 & args) :CDomainOfIntegration(args,int1d,3,2) {}
    static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args,int1d,3,2);}
    static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmeshS>(), true);} // all type
};

class CDomainOfIntegrationAllEdgesS: public CDomainOfIntegration {
public:
    CDomainOfIntegrationAllEdgesS( const basicAC_F0 & args) :CDomainOfIntegration(args,intalledges,3,2) {}
    static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args,intalledges,3,2);}
    static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmeshS>(), true);} // all type
};


// 3D curve

class CDomainOfIntegrationL: public CDomainOfIntegration {
public:
    CDomainOfIntegrationL( const basicAC_F0 & args) :CDomainOfIntegration(args,int1d,3,1) {}
    static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args,int1d,3,1);}
    static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmeshL>(), true);} // all type
};

class CDomainOfIntegrationBorderL: public CDomainOfIntegration {
public:
    CDomainOfIntegrationBorderL( const basicAC_F0 & args) :CDomainOfIntegration(args,int0d,3,1) {}
    static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args,int0d,3,1);}
    static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmeshL>(), true);} // all type
};

// hack build template
template<class T> struct CadnaType{
   typedef T  Scalaire;
 };
#ifdef  HAVE_CADNA
#include <cadnafree.h>
// specialisation
template<> struct CadnaType<complex<double> >{
   typedef complex<double_st> Scalaire;
 };
template<> struct CadnaType<double> {
   typedef double_st Scalaire;
 };
inline double_st conj(const double_st &x){ return x;};
inline complex<double_st> conj(const complex<double_st> &x){ return complex<double_st>(x.real(),-x.imag());}

 inline double norm(complex<double_st> x){return x.real()*x.real()+x.imag()*x.imag();}
 inline double norm(double_st x){return x*x;}
inline int cestac(const complex<double_st> & z)
{return min(cestac(z.real()),cestac(z.imag()));}
#endif

/*
// A construire
class CompositeOperator{
vector<BlockOperator> block_op; 
};
class BlockOperator{
  long offsetI;
  long offsetJ;
  vector<VirtualMatrix *> pBlockMatrix;
};
*/

class Problem  : public Polymorphic {
  //  typedef double R;
  static basicAC_F0::name_and_type name_param[] ;
  static const int n_name_param =3+NB_NAME_PARM_MAT; // modi FH oct 2005 add tol_pivot 02/ 2007 add nbiter
  int Nitem,Mitem;
  const int dim;
public:
 template<class FESpace>
  struct Data {
    typedef typename FESpace::Mesh Mesh;
    const Mesh * pTh;
    CountPointer<const FESpace> Uh;
    CountPointer<const FESpace> Vh;
    CountPointer<MatriceCreuse<double> > AR;
    CountPointer<MatriceCreuse<Complex> > AC;
    typedef  CadnaType<double>::Scalaire double_st;
    typedef  CadnaType<complex<double> >::Scalaire cmplx_st;
    MatriceCreuse<double_st>  * AcadnaR;
    MatriceCreuse<cmplx_st>  * AcadnaC;

    void init()  {pTh=0; AcadnaR=0;AcadnaC=0; Uh.init(),Vh.init();AR.init();AC.init();}
    void destroy() {
      pTh=0;
      Uh.destroy();
      Vh.destroy();
      AR.destroy();
      AC.destroy();

      if(AcadnaR) AcadnaR->destroy();
      if(AcadnaC) AcadnaC->destroy();
      }
  } ;
  
  struct DataComposite{
    /*
    const int Nfe; // number of FESpace in Vh
    const int Mfe; // number of FESpace in Uh
    vector< void * > *pFESpaceUh; // pointer to FESpace 
    vector< void * > *pFESpaceVh; // pointer to FESpace 

    vector< int > *typeUh; // int of the type of Vh
    vector< int > *typeVh; // int of the type of Vh

    vector< void * > *pThU; 
    vector< void * > *pThV; 

    vector<long> *sizeUh;
    vector<long> *sizeVh;

    vector<long> *offsetUh;
    vector<long> *offsetVh; 
    */

    vector< void * > *pThU; 
    vector< void * > *pThV; 

    CountPointer<MatriceCreuse<double> > ARglobal;
    CountPointer<MatriceCreuse<Complex> > ACglobal;
    void init()  {ARglobal.init();ACglobal.init();}
    void destroy() {
      if(ARglobal){ ARglobal->SetSolver(); ARglobal.destroy(); }
      if(ACglobal){ ACglobal->SetSolver(); ACglobal.destroy(); }
    }
  };
  
  const OneOperator *precon;

  C_args *op; // the list of all operator
  KNM<list<C_F0>> block_largs;
  // list<C_F0> largs;  
  mutable vector<Expression> var; // list des var pour les solutions et test
  mutable vector<int> type_var;   // list des type des FESpace for each var pour les solutions et test (used only different type of FESpace)
  bool complextype,VF;

  //   Expression noinit,type,epsilon;
  Expression nargs[n_name_param];

  const size_t  offset;
  Problem(const C_args * ca,const ListOfId &l,size_t & top) ;
  static ArrayOfaType  typeargs() { return ArrayOfaType(true);}// all type

  Data<FESpace> * dataptr (Stack stack) const {return   (Data<FESpace> *) (void *) (((char *) stack)+offset);}
  Data<FESpace3> * dataptr3 (Stack stack) const {return   (Data<FESpace3> *) (void *) (((char *) stack)+offset);}
  Data<FESpaceS> * dataptrS (Stack stack) const {return   (Data<FESpaceS> *) (void *) (((char *) stack)+offset);}
  Data<FESpaceL> * dataptrL (Stack stack) const {return   (Data<FESpaceL> *) (void *) (((char *) stack)+offset);}

  DataComposite *dataptrCompo(Stack stack) const {return  (DataComposite *) (void *) (((char *) stack)+offset);}
    
  void init(Stack stack) const  {
      //cout << " init  " << (char *) dataptr(stack) - (char*) stack  << " " << offset <<  endl;
      if(dim==2)
      dataptr(stack)->init();
      else if(dim==3)
      dataptr3(stack)->init();
      else if(dim==4)
      dataptrS(stack)->init();
      else if(dim==5)
      dataptrL(stack)->init();
      else if(dim==6) dataptrCompo(stack)->init(); // allocation pour que les offset soit correct
  }
  void destroy(Stack stack)  const  {
      if(dim==2) dataptr(stack)->destroy();
      else if(dim==3) dataptr3(stack)->destroy();
      else if(dim==4) dataptrS(stack)->destroy();
      else if(dim==5) dataptrL(stack)->destroy();
      else if(dim==6) dataptrCompo(stack)->destroy(); // allocation pour que les offset soit correct
  }

  template<class R,class FESpace,class v_fes>
  AnyType eval(Stack stack,Data<FESpace> * data,CountPointer<MatriceCreuse<R> > & dataA,
  MatriceCreuse< typename CadnaType<R>::Scalaire >  * & dataCadna) const;
  template<class R>    // TODO if coupling FE wit problem
  AnyType evalComposite(Stack stack,CountPointer<MatriceCreuse<R> > & dataA) const;
  AnyType operator()(Stack stack) const; // move in Problem.cpp  dec. 2018 FH

  bool Empty() const {return false;}
  size_t nbitem() const { return Nitem;}
};


class Solve : public Problem { public:
  // just a problem with implicit solve
  Solve(const C_args * ca,const ListOfId &l,size_t & top)
    : Problem(new C_args(*ca),l,top) {}
};

class FormBilinear : public E_F0mps { public:
  typedef const FormBilinear* Result;
  typedef const CDomainOfIntegration * A;
  typedef const foperator  *  B;
  A  di;
  Foperator * b;
  FormBilinear(const basicAC_F0 & args) {
   di= dynamic_cast<A>(CastTo<A>(args[0]));
   B  bb= dynamic_cast<B>(CastTo<B>(args[1]));
  // b = bb->Optimize(currentblock); // FH1004
  b=new Foperator(*bb); // FH1004  no optimisation here because we don't the type of the bilinear form here.
  //  the opimisation is done after in FieldOfForm routine
  // to find if the form is real or complex

  // delete bb; il ne faut pas detruire .. car bb peut etre dans une variable
    ffassert(di && b); }

  static ArrayOfaType  typeargs() { return ArrayOfaType(atype<A>(),atype<B>());}// all type
  AnyType operator()(Stack ) const { return SetAny<Result>(this);}
   operator aType () const { return atype<Result>();}

  static  E_F0 * f(const basicAC_F0 & args) { return new FormBilinear(args);}
  FormBilinear(A a,Expression bb) : di(a),b(new Foperator(*dynamic_cast<B>(bb))/*->Optimize(currentblock) FH1004 */)
  {ffassert(b);}

  // Added J. Morice to recreate largs
  // FormBilinear(A a, Foperator &bb) : di(a),b(&bb)/*->Optimize(currentblock) FH1004 */
  // {ffassert(b);}
  // Fin Added J. Morice 

  FormBilinear operator-() const { return  FormBilinear(di,C_F0(TheOperators,"-",C_F0(b,atype<B>())));}
  bool VF() const { return MaxOp(b) >= last_operatortype;}
  int dim() const {return di->d;}
  FormBilinear(const FormBilinear & fb) : di(fb.di),b(new Foperator(*fb.b) ) {}
  //  void init(Stack stack) const {}
};


//template<class v_fes>
class FormLinear : public E_F0mps { public:
  typedef const FormLinear* Result;
  typedef const CDomainOfIntegration * A;
  typedef const ftest  *  B;
  A  di;
  Ftest * l;

  FormLinear(const basicAC_F0 & args) {
    di= dynamic_cast<A>(CastTo<A>(args[0]));
    assert(di);
    Expression a1=CastTo<B>(args[1]);
    assert(a1);
   // cout << "  ---FormLinear: "<< a1 << "  " << typeid(*a1).name() << *a1 <<endl;
    B ll= dynamic_cast<B>(a1);
    assert(ll);
    l = new Ftest(*ll); // FH1004 ->Optimize(currentblock);  same as bilinear
    // delete ll; // il ne faut pas detruire car ll peut etre dans une variable
    assert(l);
    ffassert(di && l);
  }
  bool VF() const { return MaxOp(l) >= last_operatortype;}

  static ArrayOfaType  typeargs() { return ArrayOfaType(atype<A>(),atype<B>());}// all type
  AnyType operator()(Stack ) const  { return SetAny<Result>(this);}
   operator aType () const { return atype<Result>();}
  int dim() const {return di->d;}
  static  E_F0 * f(const basicAC_F0 & args) { return new FormLinear(args);}
  FormLinear(A a,Expression bb) : di(a),l(new Ftest(*dynamic_cast<B>(bb))/*->Optimize(currentblock) FH1004 */) {ffassert(l);}
  FormLinear operator-() const { return  FormLinear(di,C_F0(TheOperators,"-",C_F0(l,atype<B>())));}
  //  void init(Stack stack) const {}
  FormLinear(const FormLinear & fb) : di(fb.di),l(new Ftest(*fb.l) ) {}

};

template<class VFES>
class Call_FormLinear: public E_F0mps
{
public:
  list<C_F0> largs;
  Expression *nargs;
  typedef list<C_F0>::const_iterator const_iterator;
  const int N;
  Expression ppfes;

  Call_FormLinear(Expression * na,Expression  LL, Expression ft) ;
  AnyType operator()(Stack stack) const
  { InternalError(" bug: no eval of Call_FormLinear ");}
   operator aType () const { return atype<void>();}

};

template<class VFES1,class VFES2>
class Call_FormBilinear: public E_F0mps
{
public:
  Expression *nargs;
  list<C_F0> largs;
  typedef list<C_F0>::const_iterator const_iterator;

  const int N,M;
  Expression euh,evh;
  Call_FormBilinear(Expression * na,Expression  LL, Expression fi,Expression fj) ;
  AnyType operator()(Stack stack) const
  { InternalError(" bug: no eval of Call_FormBilinear ");}
  
operator aType () const { return atype<void>();}
};

template<class VFES1,class VFES2>
class Call_CompositeFormBilinear: public E_F0mps
{
public:
  Expression *nargs;
  //list<C_F0> varflargs;
  //list<C_F0> largs;
  KNM<list<C_F0>> block_largs;

  typedef list<C_F0>::const_iterator const_iterator;

  const int N,M;
  Expression euh,evh;
  Call_CompositeFormBilinear(Expression * na,Expression  LL, Expression fi,Expression fj) ;
  AnyType operator()(Stack stack) const
  { InternalError(" bug: no eval of Call_CompositeFormBilinear ");}
  
operator aType () const { return atype<void>();}

};


struct OpCall_FormLinear_np {
  static basicAC_F0::name_and_type name_param[] ;
  static const int n_name_param =1;

};

struct OpCall_FormBilinear_np {
  static basicAC_F0::name_and_type name_param[] ;
    static const int n_name_param =1+NB_NAME_PARM_MAT+NB_NAME_PARM_HMAT; // 9-> 11 FH 31/10/2005  11->12 nbiter 02/2007  // 12->22 MUMPS+ Autre Solveur 02/08  // 34->40 param bem solver
};

template<class T,class v_fes>
struct OpCall_FormLinear
  : public OneOperator,
    public OpCall_FormLinear_np
{
  typedef v_fes *pfes;
  E_F0 * code(const basicAC_F0 & args) const
  {
    Expression * nargs = new Expression[n_name_param];
    args.SetNameParam(n_name_param,name_param,nargs);

    return  new Call_FormLinear<v_fes>(nargs,to<const C_args*>(args[0]),to<pfes*>(args[1]));}
  OpCall_FormLinear() :
    OneOperator(atype<const Call_FormLinear<v_fes>*>(),atype<const T*>(),atype<pfes*>()) {}
};


template<class T,class v_fes>
struct OpCall_FormLinear2
  : public OneOperator,
    public OpCall_FormLinear_np
{
  typedef v_fes *pfes;
  E_F0 * code(const basicAC_F0 & args) const
  {
    Expression * nargs = new Expression[this->n_name_param];
    args.SetNameParam(this->n_name_param,this->name_param,nargs);

    Expression p=args[1];
    if ( ! p->EvaluableWithOutStack() )
      { CompileError("  a(long,Vh) , The long  must be a constant, and  = 0, sorry");}
    long pv = GetAny<long>((*p)(NullStack));
    if ( pv )
      { CompileError("  a(long,Vh) , The long  must be a constant == 0, sorry");}
    return  new Call_FormLinear<v_fes>(nargs,to<const C_args*>(args[0]),to<pfes*>(args[2]));}
  OpCall_FormLinear2() :
    OneOperator(atype<const Call_FormLinear<v_fes>*>(),atype<const T*>(),atype<long>(),atype<pfes*>()) {}
};

template<class T,class v_fes1, class v_fes2>
struct OpCall_FormBilinear
  : public OneOperator ,
  OpCall_FormBilinear_np
{
  typedef v_fes1 *pfes1; typedef v_fes2 *pfes2;

  E_F0 * code(const basicAC_F0 & args) const
  { Expression * nargs = new Expression[n_name_param];
    args.SetNameParam(n_name_param,name_param,nargs);
    // cout << " OpCall_FormBilinear " << *args[0].left() << " " << args[0].LeftValue() << endl;
    return  new Call_FormBilinear<v_fes1, v_fes2>(nargs,to<const C_args*>(args[0]),to<pfes1*>(args[1]),to<pfes2*>(args[2]));}
  OpCall_FormBilinear() :
    OneOperator(atype<const Call_FormBilinear<v_fes1,v_fes2>*>(),atype<const T *>(),atype<pfes1*>(),atype<pfes2*>()) {}
};

template<class T,class v_fes1, class v_fes2>
struct OpCall_CompositeFormBilinear
  : public OneOperator ,
  OpCall_FormBilinear_np
{
  typedef v_fes1 *pfes1; typedef v_fes2 *pfes2;
  
  E_F0 * code(const basicAC_F0 & args) const
  { Expression * nargs = new Expression[n_name_param];
    args.SetNameParam(n_name_param,name_param,nargs);

    return  new Call_CompositeFormBilinear<v_fes1, v_fes2>(nargs,to<const C_args*>(args[0]),to<pfes1*>(args[1]),to<pfes2*>(args[2]));

  }
  
  OpCall_CompositeFormBilinear() :
    OneOperator(atype<const Call_CompositeFormBilinear<v_fes1,v_fes2>*>(),atype<const T *>(),atype<pfes1*>(),atype<pfes2*>()) {}
};



bool FieldOfForm( list<C_F0> & largs ,bool complextype);
template<class A>  struct IsComplexType { static const bool value=false;};
template<>  struct IsComplexType<Complex> { static const bool value=true;};

template<class R,class MMesh,class v_fes>  //  to make   x=linearform(x)
struct OpArraytoLinearForm
  : public OneOperator
{
  typedef typename Call_FormLinear<v_fes>::const_iterator const_iterator;
  const bool isptr;
  const bool init;
  const bool zero;
  class Op : public E_F0mps
  {
  public:
    Call_FormLinear<v_fes> *l;
    Expression x;
    const  bool isptr;
    const  bool init;
    const  bool zero;
    AnyType operator()(Stack s)  const ;
    Op(Expression xx,Expression  ll,bool isptrr,bool initt,bool zzero)
      : l(new Call_FormLinear<v_fes>(*dynamic_cast<const Call_FormLinear<v_fes> *>(ll))),
        x(xx),
        isptr(isptrr),init(initt),zero(zzero)
        {assert(l);

	bool iscmplx=FieldOfForm(l->largs,IsComplexType<R>::value);
	//cout<< "FieldOfForm:iscmplx " << iscmplx << " " << IsComplexType<R>::value << " " <<( (iscmplx) == IsComplexType<R>::value) << endl;
	ffassert( (iscmplx) == IsComplexType<R>::value);
}
     operator aType () const { return atype<KN<R> *>();}

  };
  E_F0 * code(const basicAC_F0 & args) const
  { if(isptr) return new Op(to<KN<R> *>(args[0]),args[1],isptr,init,zero);
    else      return new Op(to<KN_<R> >(args[0]),args[1],isptr,init,zero);}


  //  OpArraytoLinearForm(const basicForEachType * tt) :
  //    OneOperator(atype<KN_<R> >(),tt,atype<const Call_FormLinear*>()),init(false),isptr(false) {}

  OpArraytoLinearForm(const basicForEachType * tt,bool isptrr, bool initt,bool zzero=1) :
    OneOperator(atype<KN_<R> >(),tt,atype<const Call_FormLinear<v_fes>*>()),
    isptr(isptrr), init(initt),zero(zzero) {}

};


template<class R,class MMesh,class v_fes1,class v_fes2>  //  to make   A=linearform(x)
struct OpMatrixtoBilinearForm
  : public OneOperator
{
  typedef typename Call_FormBilinear<v_fes1,v_fes2>::const_iterator const_iterator;
  int init;
  class Op : public E_F0mps {
  public:
    Call_FormBilinear<v_fes1,v_fes2> *b;
    Expression a;
    int init;
    AnyType operator()(Stack s)  const ;

    Op(Expression aa,Expression  bb,int initt)
      : b(new Call_FormBilinear<v_fes1,v_fes2>(* dynamic_cast<const Call_FormBilinear<v_fes1,v_fes2> *>(bb))),a(aa),init(initt)
  { assert(b && b->nargs);
    bool iscmplx=FieldOfForm(b->largs,IsComplexType<R>::value)  ;
     // cout<< "FieldOfForm:iscmplx " << iscmplx << " " << IsComplexType<R>::value << " " << ((iscmplx) == IsComplexType<R>::value) << endl;
    ffassert( (iscmplx) == IsComplexType<R>::value);
}
  operator aType () const { return atype<Matrice_Creuse<R>  *>();}

  };
  E_F0 * code(const basicAC_F0 & args) const
  { return  new Op(to<Matrice_Creuse<R>*>(args[0]),args[1],init);}
  OpMatrixtoBilinearForm(int initt=0) :
    OneOperator(atype<Matrice_Creuse<R>*>(),atype<Matrice_Creuse<R>*>(),atype<const Call_FormBilinear<v_fes1,v_fes2>*>()),
    init(initt)
 {}
};

template<class R>  //  to make   x=linearform(x)
struct OpArraytoLinearFormVG
  : public OneOperator
{
  typedef typename Call_FormLinear<vect_generic_v_fes>::const_iterator const_iterator;
  const bool isptr;
  const bool init;
  const bool zero;
  class Op : public E_F0mps
  {
  public:
    Call_FormLinear<vect_generic_v_fes> *l;
    Expression x;
    const  bool isptr;
    const  bool init;
    const  bool zero;

    AnyType operator()(Stack s)  const ;
    Op(Expression xx,Expression  ll,bool isptrr,bool initt,bool zzero)
      : l(new Call_FormLinear<vect_generic_v_fes>(*dynamic_cast<const Call_FormLinear<vect_generic_v_fes> *>(ll))),
        x(xx),
        isptr(isptrr),init(initt),zero(zzero)
        {
          assert(l);

        	bool iscmplx=FieldOfForm(l->largs,IsComplexType<R>::value);
	        //cout<< "FieldOfForm:iscmplx " << iscmplx << " " << IsComplexType<R>::value << " " <<( (iscmplx) == IsComplexType<R>::value) << endl;
	        ffassert( (iscmplx) == IsComplexType<R>::value);

        }
    operator aType () const { return atype<KN<R> *>();}

  };
  E_F0 * code(const basicAC_F0 & args) const
  { if(isptr) return new Op(to<KN<R> *>(args[0]),args[1],isptr,init,zero);
    else      return new Op(to<KN_<R> >(args[0]),args[1],isptr,init,zero);}


  //  OpArraytoLinearForm(const basicForEachType * tt) :
  //    OneOperator(atype<KN_<R> >(),tt,atype<const Call_FormLinear*>()),init(false),isptr(false) {}

  OpArraytoLinearFormVG(const basicForEachType * tt,bool isptrr, bool initt,bool zzero=1) :
    OneOperator(atype<KN_<R> >(),tt,atype<const Call_FormLinear<vect_generic_v_fes>*>()),
    isptr(isptrr), init(initt),zero(zzero) {}

};

/*
template<class R>  //  to make   A=linearform(x)
struct OpMatrixtoBilinearFormVG
  : public OneOperator
{
  typedef typename Call_FormBilinear<vect_generic_v_fes,vect_generic_v_fes>::const_iterator const_iterator;
  int init;
  
  class Op : public E_F0mps {
    public:
      Call_FormBilinear<vect_generic_v_fes,vect_generic_v_fes> *b;
      Expression a;
      int init;
      //AnyType operator()(Stack s)  const;
      
      Op(Expression aa,Expression  bb,int initt)
        : b(new Call_FormBilinear<vect_generic_v_fes,vect_generic_v_fes>(* dynamic_cast<const Call_FormBilinear<vect_generic_v_fes,vect_generic_v_fes> *>(bb))),a(aa),init(initt)
    { 
      assert(b && b->nargs);
      bool iscmplx=FieldOfForm(b->largs,IsComplexType<R>::value)  ;
        // cout<< "FieldOfForm:iscmplx " << iscmplx << " " << IsComplexType<R>::value << " " << ((iscmplx) == IsComplexType<R>::value) << endl;
      ffassert( (iscmplx) == IsComplexType<R>::value);
    }
    operator aType () const { return atype<Matrice_Creuse<R>  *>();}

    AnyType operator()(Stack s)  const;
  };

  E_F0 * code(const basicAC_F0 & args) const
  { return  new Op(to<Matrice_Creuse<R>*>(args[0]),args[1],init); }
  OpMatrixtoBilinearFormVG(int initt=0) :
    OneOperator(atype<Matrice_Creuse<R>*>(),atype<Matrice_Creuse<R>*>(),atype<const Call_FormBilinear<vect_generic_v_fes,vect_generic_v_fes>*>()),
    init(initt){};

};
*/
template<class R>
class IntFunction  : public E_F0mps { public:
  typedef R Result;
  typedef const CDomainOfIntegration * A;
  typedef R  B;
  A  di;
  Expression fonc;
  IntFunction(const basicAC_F0 & args) {
    di= dynamic_cast<A>(CastTo<A>(args[0]));
    fonc= CastTo<B>(args[1]);
    ffassert(di && fonc); }

  static ArrayOfaType  typeargs() { return ArrayOfaType(atype<A>(),atype<B>());}// all type
  AnyType operator()(Stack ) const;
  static  E_F0 * f(const basicAC_F0 & args) { return new IntFunction(args);}
  //  IntFunction(A a,Expression bb) : di(a),fonc(bb) {}
   operator aType () const { return atype<Result>();}

};


extern Block *currentblock;

class TypeFormOperator: public ForEachType<const C_args*> {
public:
  TypeFormOperator() : ForEachType<const C_args*>(0,0) {}
  void SetArgs(const ListOfId *lid) const {
    SetArgsFormLinear(lid,2);    }

  Type_Expr SetParam(const C_F0 & c,const ListOfId *l,size_t & top) const
  { return Type_Expr(this,CastTo(c));}

  inline  C_F0 Initialization(const Type_Expr & e) const {return C_F0();}

};

class TypeFormBilinear: public ForEachType<const FormBilinear*> {
public:
  TypeFormBilinear() : ForEachType<const FormBilinear*>(0,0) {}
  void SetArgs(const ListOfId *lid) const {
    SetArgsFormLinear(lid,2);
  }

  Type_Expr SetParam(const C_F0 & c,const ListOfId *l,size_t & top) const
  { return Type_Expr(this,CastTo(c));}


  C_F0 Initialization(const Type_Expr & e) const
  {
    // cout << "Initialization " << *e.first << endl;
    return C_F0(); }  // nothing to initialize
  Type_Expr construct(const Type_Expr & e) const
  {
    //cout << "construct " << *e.first << endl;
    return e; }

};

template<bool exec_init,class Problem>
class TypeSolve : public ForEachType<const Problem*>   {
public:
  TypeSolve() : ForEachType<const Problem*>(0,0) {}

  void SetArgs(const ListOfId *lid) const {
    SetArgsFormLinear(lid,2);


  }
  Type_Expr SetParam(const C_F0 & c,const ListOfId *l,size_t & top) const
  {   if (c.left() != atype<const C_args*>())
    CompileError(" Problem  a(...) = invalid type ",c.left());
  const C_args * ca = dynamic_cast<const C_args *>(c.LeftValue());
  Problem * pb=new Problem(ca,*l,top);
  SHOWVERB(cout << "solve:SetParam " << ca << " pb=" << pb << endl);

  return Type_Expr(this,pb);

  }

  class SolveInit: public  E_F0 { public:
    const Problem * a;
    AnyType operator()(Stack s)  const {
      a->init(s);
      return exec_init ? (*a)(s) : Nothing ;
    }
    SolveInit(const Type_Expr &  te) : a(dynamic_cast<const Problem *>(te.second))
    {  SHOWVERB(cout << "SolveInit " << te.second << endl);
    ffassert(a);}
  };

  class SolveDel: public  E_F0 { public:
    const Problem * a;
    SolveDel(const C_F0 & c) : a(dynamic_cast<const Problem *>(c.LeftValue()))
    {
      SHOWVERB(cout << "SolveDel " << c.left()  << endl);
      ffassert(a);}

    AnyType operator()(Stack s)  const {
      a->destroy(s);
      return Nothing;
    }};

  Expression Destroy(const C_F0 & c) const
  { return new SolveDel(c);}

  bool ExistDestroy() const {return true;}

  C_F0 Initialization(const Type_Expr & e) const
  {  return C_F0( new SolveInit(e) ,atype<void>()); }
};




class TypeFormLinear: public ForEachType<const FormLinear*> {
public:
  TypeFormLinear() : ForEachType<const FormLinear*>(0,0) {}

  void SetArgs(const ListOfId *lid) const {
    SetArgsFormLinear(lid,1);  }

  Type_Expr SetParam(const C_F0 & c,const ListOfId *l,size_t & top) const
  { return Type_Expr(this,CastTo(c));}
  //  {  return Type_Expr(c.left(),c.LeftValue());  } //

  C_F0 Initialization(const Type_Expr & e) const
  {  return C_F0(); }  // nothing to initialize

};



template<class K> class Matrice_Creuse
{
public:
  UniqueffId Uh,Vh; // pour la reconstruction
    // MatriceCreuse<K> == VirtualMatrix<int,K>
 typedef  VirtualMatrix<int,K> VMat;
  typedef  HashMatrix<int,K> HMat;

 CountPointer<MatriceCreuse<K> > A;
    int typemat; //
     static const int TS_SYM=1,TS_DEF_POS=2,TS_PARA=4;
  //TypeSolveMat typemat;
    size_t count;
  void init() {
      count=0;
    A.init();Uh.init();Vh.init();
      typemat=0 ; }//
  Matrice_Creuse() { init();}
  void destroy() {// Correct Oct 2015 FH (avant test a 'envert) !!!!
    if(verbosity>99999)
        cerr << " ## DEL MC " << this <<" " << count <<" " << A <<  endl;
    if(count--==0)
    {
        VMat *pvm=pMC();
        if(pvm) pvm->SetSolver();
        pvm=0; // del solver before the del of the real matrix...
       A.destroy();
    }
//else count--;
    //    Uh.destroy();
    //Vh.destroy();
  }
  void copysolver(Matrice_Creuse *a)
    {
        VMat *pvm=pMC(), *pvam=a->pMC();
        if( pvm)
        {
            //typename VMat::VSolver s= pvam ?   pvam->CloneSolver() : 0;
            //  to hard  to CloneSolver
            pvm->SetSolver();//So  Solver ..
            pvm=0;
        }
    }
  Matrice_Creuse( MatriceCreuse<K> * aa)//,const pfes  *ppUh,const pfes  *ppVh)
    :A(aa){}//,pUh(ppUh),pVh(ppVh),Uh(*ppUh),Vh(*ppVh) {}
  Matrice_Creuse( MatriceCreuse<K> * aa,const UniqueffId *pUh,const UniqueffId *pVh)//,const pfes  *ppUh,const pfes  *ppVh)
    :A(aa),Uh(*pUh),Vh(*pVh) {}//,pUh(ppUh),pVh(ppVh),Uh(*ppUh),Vh(*ppVh) {}
  long N() const {return  A ? A->n : 0;}
  long M() const { return A ? A->m : 0;}
  void resize(int n,int m) {
      if(A) A->resize(n,m);
      else {//  matrice vide a cree
          HashMatrix<int,K> *phm= new HashMatrix<int,K>(n,m,0,0);
          MatriceCreuse<K> *pmc(phm);
          A.master(pmc);
      }
      
  }
  void increment(){ count++;}
    VMat *pMC()  {return A ? ( MatriceCreuse<K> *)A:0; }
    HMat *pHM()  {return dynamic_cast<HashMatrix<int,K> *>(pMC());}
};

template<class K> class newpMatrice_Creuse
{
public:
    MatriceCreuse<K> *pmc;
    newpMatrice_Creuse(Stack s,HashMatrix<int,K> *pvm) :pmc(pvm)
    {

        if(verbosity>99999)  cerr << " newpMatrice_Creuse Add2StackOfPtr2FreeRC "<< pmc  << endl;
        Add2StackOfPtr2FreeRC(s,pmc);
       // Add2StackOfPtr2Free(s,pmc);
    }
    Matrice_Creuse<K> * set(Matrice_Creuse<K> *pmcc,int init)   {
        if(init) pmcc->init() ;
        pmc->increment() ;
        pmcc->A.master(pmc);
      //  pmcc->A.cswap(pmc);

        if(verbosity>99999)   cerr << "newpMatrice_Creuse  set " << pmcc << " " << pmcc->count <<" " << pmcc->A
        << " to " << pmc  << " init: "<< init << endl;
       // pmc->dump(cerr) << endl;
         pmc=0;
        return  pmcc;
    }
    Matrice_Creuse<K> * add(Matrice_Creuse<K> *pmcc,double cc=1)   {
        //  pmcc->A.cswap(pmc);
        HashMatrix<int,K> *pA= pmcc->pHM();
        HashMatrix<int,K> *pC= dynamic_cast<HashMatrix<int,K> *>(pmc);
        
        if(pA==0) { cerr<< " error A  += B (A) empty matrix" << endl; }
        if(pC==0) { cerr<< " error += B  (B)  empty matrix" << endl; }
        ffassert(pA && pC);
        HashMatrix<int,K> A=*pA;
        pA->Add(pC,cc);
      // to do.. XXXX  July 2017 FH.
        if(verbosity>99999)   cerr << "newpMatrice_Creuse  add " << pmcc << " " << pmcc->count <<" " << pmcc->A
            << " to " << pmc  << endl;
        
        pmc=0;//  pcm is delete after instriction
        return  pmcc;
    }
  //  ~newpMatrice_Creuse() { if(pmc) delete pmc;pmc=0; }
};

template<class K> class Matrice_Creuse_Transpose;

 template<class KA,class KB>   class Matrix_Prod { public:
  Matrice_Creuse<KA> *A;
  Matrice_Creuse<KB> *B;
  bool ta,tb;
  Matrix_Prod(Matrice_Creuse<KA> *AA,Matrice_Creuse<KB> *BB) : A(AA),B(BB),ta(false),tb(false) {assert(AA && BB);}
  Matrix_Prod(Matrice_Creuse_Transpose<KA> AA,Matrice_Creuse<KB> *BB)           : A(AA),B(BB),ta(true),tb(false) {assert(AA && BB);}
  Matrix_Prod(Matrice_Creuse<KA> *AA,Matrice_Creuse_Transpose<KB> BB)           : A(AA),B(BB),ta(false),tb(true) {assert(AA && BB);}
  Matrix_Prod(Matrice_Creuse_Transpose<KA> AA,Matrice_Creuse_Transpose<KB> BB) : A(AA),B(BB),ta(true),tb(true) {assert(AA && BB);}
 };

template<class K>  ostream & operator << (ostream & f,const Matrice_Creuse<K> & A)
{ if ( !A.A) f << " unset sparse matrix " << endl;
  else A.A->dump(f);  ;
 return f;  }

template<class K>  istream & operator >> (istream & f,Matrice_Creuse<K> & A)
{
    int wm=WhichMatrix(f);
    if ( wm>0 )
    {
        MatriceMorse<K> *HA =0;
  	A.A.master(HA=new MatriceMorse<K>(f,wm));
        A.typemat=HA->sym() ;//(A.A->n == A.A->m) ? TypeSolveMat(TypeSolveMat::GMRES) : TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
    }
    else {
	cerr << " unknown type of matrix " << wm <<endl;
	ExecError("Error reading the matrix ");
	A.A =0; }
    return f;  }

template<class K> class Matrice_Creuse_Transpose  { public:
  Matrice_Creuse<K> * A;

  Matrice_Creuse_Transpose(Matrice_Creuse<K> * AA) : A(AA) {assert(A);}
  operator MatriceCreuse<K> & () const {return *A->A;}
  operator Matrice_Creuse<K> * () const {return A;}
};

template<class K> class Matrice_Creuse_inv  { public:
  Matrice_Creuse<K> * A;
  Matrice_Creuse_inv(Matrice_Creuse<K> * AA) : A(AA) {assert(A);}
  operator MatriceCreuse<K> & () const {return *A->A;}
  operator Matrice_Creuse<K> * () const {return A;}
};

template<class K> class Matrice_Creuse_inv_trans  { public:// add aug 2018 FH.
    Matrice_Creuse<K> * A;
    Matrice_Creuse_inv_trans(Matrice_Creuse<K> * AA) : A(AA) {assert(A);}
    Matrice_Creuse_inv_trans(Matrice_Creuse_Transpose<K> * AA) : A(AA) {assert(A);}
    Matrice_Creuse_inv_trans(const Matrice_Creuse_Transpose<K> & AA) : A(AA.A) {assert(A);}
    operator MatriceCreuse<K> & () const {return *A->A;}
    operator Matrice_Creuse<K> * () const {return A;}
};
/*
//==== PAS TESTER A ENLEVER MORICE  ====//
class qOperateurMatrice{
  public:
  // classe operateur quelquonque pour la construction d'une matrice
  void * mat;
  // mat : pointer on a HMatrix or a Matrice_Creuse
  void * get_mat(){return mat;};
};

// structure for composite Matrix or block  matrix
template<class R>
class BlockCompositeMatrice {
  int nbb; // number of blocks
  
  void ** data; // point to list of array of matrix
  // Matrix can be a HMatrix or a Matrice_Creuse
  
  KN<int> offsetI;
  KN<int> offsetJ;

  BlockCompositeMatrice(int n_nbb, KN<int> n_offsetI, KN<int> n_offsetJ, void ** m_data ) :
    nbb(n_nbb), offsetI(n_offsetI), offsetJ(n_offsetJ), data(m_data)
  {
    // check type of the composite matrix
  }

  ~BlockCompositeMatrice(){ 
    // on ne peut pas dealllouer les matrices car dans certain cas on en aurra besoin.
    data=nullptr;   
  }

  Matrice_Creuse<R> * getGlobalMatrix(){
    // construct the globel matrix correspond to the block matrix

    Matrice_Creuse<R> *A = new Matrice_Creuse<R>() ;
    A->resize( offsetJ.sum(), offsetI.sum() ); // test function (Vh) are the line and inconnu function (Uh) are the column

    A->init();                           
    cout << " A.N=" <<  A->N() << endl;
    cout << " A.M=" <<  A->M() << endl;

    // loop over the block
    for( int i=0; i<nbb; i++){
        if( data[i] == atype < Matrice_Creuse<R>* > ){
          A->pHM()->Add( data[i]->pHM(), R(1), false, offsetJ[j], offsetI[i] );
        }
#ifndef FFLANG
#ifdef PARALLELE
        if( data[i] == atype < HMatrixVirt<R> ** Hmat>  ){
          // creation de la matrice dense 
          KNM<R>* M= HMatrixVirtToDense< KNM<R>, R >( data[i] );

          HashMatrix<int,R> *phm= new HashMatrix<int,R>(*M);
          MatriceCreuse<R> *pmc(phm);

          Matrice_Creuse<R> BBB;
          BBB.A=0;
          BBB.A.master(pmc);
          A->pHM()->Add( BBB.pHM(), R(1), false, offsetJ[j], offsetI[i] );
        }
        else{
#endif // PARALLELE
#endif // FFLANG
        else{
          cerr << "Global Matrix Construction :: error of the type of one block of the matrix" << endl;
          ffassert(0);
        }
      
    }
    return A;
  }

};
*/


namespace Fem2D {

  inline void F_Pi_h(R* v, const R2 & P,const baseFElement & K,int i,const R2 & Phat,void * arg)
  {
    TabFuncArg &tabe(*(TabFuncArg*)arg);
    //MeshPoint & mp = *MeshPointStack(tabe.s);
    MeshPointStack(tabe.s)->set(P,Phat,K);
    tabe.eval(v);
    // if (Norme2_2(P-mp.P) > 1e-10)
    //  cout << " bug?? F_Pi_h " << endl;

  }

  inline void FoX_1_Pi_h(R* v, const R2 & P,const baseFElement & K,int i,const R2 & Phat,void * arg)
  {
    TabFuncArg &tabe(*(TabFuncArg*)arg);
    MeshPointStack(tabe.s)->set(P,Phat,K);
    R2 X=tabe.eval_X();
    MeshPointStack(tabe.s)->set(X.x,X.y);
    tabe.eval_2(v);
  }

//general templates for 2d and 3d volume // for Surf version ...
template<class R,typename MC,class MMesh,class FESpace1,class FESpace2>  bool AssembleVarForm(Stack stack,const MMesh & Th,
                                                                   const FESpace1 & Uh,const FESpace2 & Vh,bool sym,
                                                                   MC  * A,KN_<R> * B,const list<C_F0> &largs );

template<class R,class MMesh,class FESpace1,class FESpace2>   void AssembleBC(Stack stack,const MMesh & Th,
                                                  const FESpace1 & Uh,const FESpace2 & Vh,bool sym,
                                                  MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X,
                                                  const list<C_F0> &largs , double tgv  );
template<class R,class MMesh,class FESpace1,class FESpace2>  void AssembleBC(Stack stack,const MMesh & Th,const FESpace1 & Uh,
                                                  const FESpace2 & Vh,bool sym, MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X,
                                                  const list<C_F0> &largs , double tgv  );


// 2d case
template<class R>   void AssembleLinearForm(Stack stack,const Mesh & Th,const FESpace & Vh,KN_<R> * B,const  FormLinear * const l);
template<class R>   void  Element_rhs(const FElement & Kv,int ie,int label,const LOperaD &Op,double * p,void * stack,KN_<R> & B,bool all,int optim);
template<class R>   void  Element_rhs(const FElement & Kv,const LOperaD &Op,double * p,void * stack,KN_<R> & B,int optim);
template<class R>   void  Element_Op(MatriceElementairePleine<R,FESpace> & mat,const FElement & Ku,const FElement & Kv,double * p,
                                     int ie,int label, void *stack,R2 *B);
template<class R>   void  Element_Op(MatriceElementaireSymetrique<R,FESpace> & mat,const FElement & Ku,double * p,int ie,int label,
                                     void * stack,R2 *B);
template<class R>   void AssembleBC(Stack stack,const Mesh & Th3,const FESpace & Uh3,const FESpace & Vh3,bool sym,
                                    MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X, const  BC_set * bc , double tgv   );


// 3d volume case
template<class R>   void AssembleLinearForm(Stack stack,const Mesh3 & Th,const FESpace3 & Vh,KN_<R> * B,const  FormLinear * const l);
template<class R>   void  Element_rhs(const FElement3 & Kv,int ie,int label,const LOperaD &Op,double * p,void * stack,KN_<R> & B,bool all,int optim);
template<class R>   void  Element_rhs(const FElement3 & Kv,const LOperaD &Op,double * p,void * stack,KN_<R> & B,int optim);
template<class R>   void  Element_Op(MatriceElementairePleine<R,FESpace3> & mat,const FElement3 & Ku,const FElement3 & Kv,double * p,
                                     int ie,int label, void *stack,R3 *B);
template<class R>   void  Element_Op(MatriceElementaireSymetrique<R,FESpace3> & mat,const FElement3 & Ku,double * p,int ie,int label,
                                     void * stack,R3 *B);

// 3d surface case
template<class R>   void AssembleLinearForm(Stack stack,const MeshS & Th,const FESpaceS & Vh,KN_<R> * B,const  FormLinear * const l);
template<class R>   void  Element_rhs(const FElementS & Kv,int ie,int label,const LOperaD &Op,double * p,void * stack,KN_<R> & B,bool all,int optim);
template<class R>   void  Element_rhs(const FElementS & Kv,const LOperaD &Op,double * p,void * stack,KN_<R> & B,int optim);
template<class R>   void  Element_Op(MatriceElementairePleine<R,FESpaceS> & mat,const FElementS & Ku,const FElementS & Kv,double * p,
                                     int ie,int label, void *stack,R3 *B);
template<class R>   void  Element_Op(MatriceElementaireSymetrique<R,FESpaceS> & mat,const FElementS & Ku,double * p,int ie,int label,
                                    void * stack,R3 *B);

// 3d curve case
template<class R>   void AssembleLinearForm(Stack stack,const MeshL & Th,const FESpaceL & Vh,KN_<R> * B,const  FormLinear * const l);
template<class R>   void  Element_rhs(const FElementL & Kv,int ie,int label,const LOperaD &Op,double * p,void * stack,KN_<R> & B,bool all,int optim);
template<class R>   void  Element_rhs(const FElementL & Kv,const LOperaD &Op,double * p,void * stack,KN_<R> & B,int optim);
template<class R>   void  Element_Op(MatriceElementairePleine<R,FESpaceL> & mat,const FElementL & Ku,const FElementL & Kv,double * p,
                                     int ie,int label, void *stack,R3 *B);
template<class R>   void  Element_Op(MatriceElementaireSymetrique<R,FESpaceL> & mat,const FElementL & Ku,double * p,int ie,int label,
                                     void * stack,R3 *B);
    
}


template<class R,class MMesh,class v_fes>
AnyType OpArraytoLinearForm<R,MMesh,v_fes>::Op::operator()(Stack stack)  const
{
  typedef v_fes *pfes;
  typedef typename  v_fes::FESpace FESpaceT;
  typedef typename  FESpaceT::FElement FElementT;
  /*typedef typename  MMesh::Element ElementT;
  typedef typename  MMesh::Vertex VertexT;
  typedef typename  MMesh::RdHat RdHatT;
  typedef typename  MMesh::Rd RdT;*/

  pfes  &  pp= *GetAny<pfes * >((*l->ppfes)(stack));
  FESpaceT * pVh = *pp ;
  FESpaceT & Vh = *pVh ;
  double tgv= ff_tgv;
  if (l->nargs[0]) tgv= GetAny<double>((*l->nargs[0])(stack));
  long NbOfDF =  pVh ? Vh.NbOfDF: 0;
  KN<R> *px=0;
  if(isptr)
    {
     px = GetAny<KN<R> * >((*x)(stack) );
     if(init )
       px->init(NbOfDF);

     if(px->N() != NbOfDF) //add Dec 2009
     {
         if(!zero ) ExecError("Error in OpArraytoLinearForm   += not correct size:  n != NbOfDF !");
	 px->resize(NbOfDF);
     }
     ffassert(px->N() == NbOfDF);
   }
  KN_<R>  xx( px ? *(KN_<R> *) px : GetAny<KN_<R> >((*x)(stack) ));
  if(zero && NbOfDF )
   xx=R();

  if (  pVh && AssembleVarForm<R,MatriceCreuse<R>,MMesh,FESpaceT,FESpaceT>(stack,Vh.Th,Vh,Vh,false,0,&xx,l->largs) )
    AssembleBC<R,MMesh,FESpaceT,FESpaceT>(stack,Vh.Th,Vh,Vh,false,0,&xx,0,l->largs,tgv);
  return SetAny<KN_<R> >(xx);
}

/**
       *  @brief  Builds a new largs  whit each element are included in one block
       *  @param  largs list of argument of the initial Forms 
       *  @param  NpVh  number of FESpace in Vh 
       *  @param  indexBlockVh give the index of the block for a given component of FESpace Vh  
       * 
       */

// 
inline KN< list<C_F0> > creationLinearFormCompositeFESpace( const list<C_F0> & largs, const int& NpVh, const KN<int> & indexBlockVh
                                                           ,const KN<int> &localIndexInTheBlockVh )
{ 
  // Vh is the FESpace
  // ================================================

  list<C_F0> newlargs; // creation de la nouvelle list largs
  KN< list<C_F0> > block_largs( (long)NpVh ); 

  list<C_F0>::const_iterator b_ii,b_ib=largs.begin(),b_ie=largs.end(); 
  int count=0;
  for (b_ii=b_ib;b_ii != b_ie;b_ii++){
    Expression e=b_ii->LeftValue();
    aType r = b_ii->left();
    // Case FormLinear
    if (r==atype<const  FormLinear *>() ){
      const FormLinear * ll=dynamic_cast<const  FormLinear *>(e);
      LOperaD * Op = const_cast<LOperaD *>(ll->l);
      if (Op == NULL) {
        if(mpirank == 0) cout << "dynamic_cast error" << endl; 
        ffassert(0);
      }

      size_t Opsize= Op->v.size();
      
      // creation the vector of the indexBlock for each element OpChange
      std::vector< int > indexBlock(Opsize);

      for(size_t jj=0; jj<Opsize; jj++){
        indexBlock[jj] = indexBlockVh[ Op->v[jj].first.first ];
      }

      LOperaD * OpBloc= new LOperaD();

      // put the term inside OpChange in the good block
      for(int jbloc=0; jbloc<NpVh; jbloc++){
        int countOP=0;
        //LOperaD * OpBloc= new LOperaD();
        for(size_t jj=0; jj<Opsize; jj++){
          if( indexBlock[jj] == jbloc ){
              if (countOP == 0){
                ffassert( OpBloc->v.empty() );
              } 
              OpBloc->add(Op->v[jj].first,Op->v[jj].second); // Add the LinearOperator to bloc (ibloc,jbloc).
              countOP += 1;
          }
        }
        if( countOP > 0 ){   
          // <<  countOP << " voila titi " << "OpBloc->v.size()= " << OpBloc->v.size() << endl; 
          ffassert( OpBloc->v.size() > 0);     
          for(size_t jj=0; jj<OpBloc->v.size(); jj++){
            OpBloc->v[jj].first.first = localIndexInTheBlockVh( OpBloc->v[jj].first.first );
          }

          newlargs.push_back( C_F0( new FormLinear( (ll->di), OpBloc ), r ) );
          block_largs(jbloc).push_back( C_F0( new FormLinear( (ll->di), OpBloc ), r ) );

          // cout << "OpBloc->v.size()=" << OpBloc->v.size() << ", OpBloc->v.empty()" << OpBloc->v.empty() << endl;
          // cout << "clear v" << endl;
          OpBloc->v.clear(); 
          // cout << "OpBloc->v.size()=" << OpBloc->v.size() << ", OpBloc->v.empty()" << OpBloc->v.empty() << endl;
          ffassert( ( OpBloc->v.empty() == true)  ); // check if OpBloc is empty after clear();
        }
      }
      delete OpBloc;
    }      
    // case BC_set 
    else if(r == atype<const  BC_set  *>()){
      const BC_set * oldbc=dynamic_cast< const BC_set *>(e); // on ne peut pas utiliser " const BC_set * " ou autrement erreur ce ompilation:  Morice
      
      BC_set * bc = new BC_set(*oldbc);
    
      int sizebc=bc->bc.size();
     
      std::vector< int > indexBlock(sizebc);
      // calculate the index of the componenent where the bloc
      for (int k=0; k<sizebc; k++)
      {
        pair<int,Expression> &xx=bc->bc[k];
        indexBlock[k] = indexBlockVh[xx.first]; 
        xx.first = localIndexInTheBlockVh(xx.first);
      }
      
      for (int k=0; k<sizebc; k++)
      {
        cout << "bc["<<k<<"]="<< bc->bc[k].first << endl;
      }
      
      // Add the bc too the correct block
      bool addBC = false;
      bool *okBC =new bool[NpVh];
      for(int ibloc=0; ibloc<NpVh; ibloc++){
        addBC = false;
        // construction of okBC for ibloc 
        for (int k=0; k<sizebc; k++){
          if(indexBlock[k] == ibloc){
            okBC[k] = true;
            addBC = true;
          }
          else{
            okBC[k] = false;
          }
        }
        if(addBC){
          cout << "addBC to block: " << ibloc << endl; 
          for(int jj=0; jj<NpVh; jj++ ){ cout << "okBC["<<jj<<"]=: " << okBC[jj] << endl;}
        }
        // add the BC_set correspond to the ibloc of the composite FESpace 
        if(addBC) newlargs.push_back( C_F0( new BC_set(*bc,okBC), r) ); 
        if(addBC) block_largs(ibloc).push_back( C_F0( new BC_set(*bc,okBC), r) );
      }
      delete [] okBC;
    }
  }
  //return newlargs;
  return block_largs;
}

template<class R>
AnyType OpArraytoLinearFormVG<R>::Op::operator()(Stack stack)  const
{
  // typedef typename  v_fes::FESpace FESpaceT;
  // get tgv value 
  double tgv= ff_tgv;
  if (l->nargs[0]) tgv= GetAny<double>((*l->nargs[0])(stack));

  // get the composite FESpace
  pvectgenericfes  *  pCompoVh= GetAny<pvectgenericfes * >((*l->ppfes)(stack));
  ffassert( *pCompoVh );
  (*pCompoVh)->update();

// get info the composite FESpace
  int  NpVh = (*pCompoVh)->N;
  KN<int> VhNbItem = (*pCompoVh)->vectOfNbitem();
  KN<int> VhNbOfDf = (*pCompoVh)->vectOfNbOfDF();
  // get x from the stack 
  long totalNbOfDF =  (*pCompoVh) ? (*pCompoVh)->totalNbOfDF(): 0;
  KN<R> *px=0;
  if(isptr)
    {
     px = GetAny<KN<R> * >((*x)(stack) );
     if(init )
       px->init(totalNbOfDF);

     if(px->N() != totalNbOfDF) //add Dec 2009
     {
        if(!zero ) ExecError("Error in OpArraytoLinearForm   += not correct size:  n != NbOfDF !");
	      px->resize(totalNbOfDF);
     }
     ffassert(px->N() == totalNbOfDF);
   }
  // construction of the Array of Composite FESpace
  KN_<R>  xx( px ? *(KN_<R> *) px : GetAny<KN_<R> >((*x)(stack) ));
  if(zero && totalNbOfDF )
   xx=R();

  // get info the composite FESpace
  /*
  int  NpVh = (*pCompoVh)->N;
  KN<int> VhNbItem = (*pCompoVh)->vectOfNbitem();
  KN<int> VhNbOfDf = (*pCompoVh)->vectOfNbOfDF();
  */
  int VhtotalNbItem=0;
  for(int i=0; i< NpVh; i++) VhtotalNbItem += VhNbItem[i];


  KN<int> beginBlockVh( NpVh );
  KN<int> indexBlockVh( VhtotalNbItem );
  KN<int> localIndexInTheBlockVh( VhtotalNbItem );

  int current_index=0;
  for(int i=0; i<NpVh; i++){
    beginBlockVh[i] = current_index;
    for(int j=0; j< VhNbItem[i]; j++){
      indexBlockVh[current_index] = i;
      localIndexInTheBlockVh[current_index] = j;
      current_index++;
    }

  }
  ffassert( current_index == VhtotalNbItem );


  // creation of the LinearForm corresponding to each block
  KN< list<C_F0> > block_largs = creationLinearFormCompositeFESpace( l->largs, NpVh, indexBlockVh, localIndexInTheBlockVh);
  
  // ===  loop over the block ===//
  int offsetVh = 0;
  for( int j=0; j<NpVh; j++){
    // get the information of the block
    const list<C_F0> & b_largs=block_largs(j); 
    int M_block = VhNbOfDf[j];

    if( b_largs.size() > 0 ){
      cout << "construction of the block "<< j <<" of the varf" << endl;
      //KN<R> * xxblock = new KN<R>(M_block);
      //KN_<R> &xxblock2 = *(KN_<R> *) xxblock;

      KN<R> xxblock2(M_block);
      xxblock2=R();   // initiallise the block to zero ??? (get previous value of xxblock2).

      if( (*pCompoVh)->typeFE[j] == 2 ){
        // 2d Mesh 
        FESpace * pVh = (FESpace*) (*pCompoVh)->vect[j]->getpVh();;
        FESpace & Vh = *pVh;

        if (  pVh && AssembleVarForm<R,MatriceCreuse<R>,Mesh,FESpace,FESpace>(stack,Vh.Th,Vh,Vh,false,0,&xxblock2,b_largs) )
          AssembleBC<R,Mesh,FESpace,FESpace>(stack,Vh.Th,Vh,Vh,false,0,&xxblock2,0,b_largs,tgv);

      }else if( (*pCompoVh)->typeFE[j] == 3 ){
        // 3d Mesh
        FESpace3 * pVh = (FESpace3*) (*pCompoVh)->vect[j]->getpVh();;
        FESpace3 & Vh = *pVh;

        if (  pVh && AssembleVarForm<R,MatriceCreuse<R>,Mesh3,FESpace3,FESpace3>(stack,Vh.Th,Vh,Vh,false,0,&xxblock2,b_largs) )
          AssembleBC<R,Mesh3,FESpace3,FESpace3>(stack,Vh.Th,Vh,Vh,false,0,&xxblock2,0,b_largs,tgv);

      }else if( (*pCompoVh)->typeFE[j] == 4 ){
        // 3d Surface Mesh 
        FESpaceS * pVh = (FESpaceS*) (*pCompoVh)->vect[j]->getpVh();;
        FESpaceS & Vh = *pVh;

        if (  pVh && AssembleVarForm<R,MatriceCreuse<R>,MeshS,FESpaceS,FESpaceS>(stack,Vh.Th,Vh,Vh,false,0,&xxblock2,b_largs) )
          AssembleBC<R,MeshS,FESpaceS,FESpaceS>(stack,Vh.Th,Vh,Vh,false,0,&xxblock2,0,b_largs,tgv);

      }else if( (*pCompoVh)->typeFE[j] == 5 ){
        // 3d Curve Mesh
        FESpaceL * pVh = (FESpaceL*) (*pCompoVh)->vect[j]->getpVh();;
        FESpaceL & Vh = *pVh;

        if (  pVh && AssembleVarForm<R,MatriceCreuse<R>,MeshL,FESpaceL,FESpaceL>(stack,Vh.Th,Vh,Vh,false,0,&xxblock2,b_largs) )
          AssembleBC<R,MeshL,FESpaceL,FESpaceL>(stack,Vh.Th,Vh,Vh,false,0,&xxblock2,0,b_largs,tgv);

      }
      else{
        cerr << "Error in the definition in the type of FESpace."  << endl;
        ffassert(0);
      }
      for(int indexJ=0; indexJ<M_block; indexJ++){
        xx[indexJ+offsetVh] += xxblock2[indexJ]; 
      }
      //delete xxblock;
    }
    offsetVh += M_block;
  }

  return SetAny<KN_<R> >(xx);
}



template<typename KK,typename vv_fes,typename CC>
struct FF_L_args {
    typedef  KK K;
    typedef vv_fes v_fes;
    typedef v_fes *pfes;
    typedef typename  v_fes::FESpace FESpaceT;
    typedef typename  FESpaceT::Mesh MeshT;
    typedef FEbase<K,v_fes> ** R;
    typedef R  A;
    typedef pfes* B;
    typedef const CC *  C;
    typedef  CC *  MC;
    static bool Check(MC l) {return IsComplexType<K>::value==FieldOfForm(l->largs,IsComplexType<K>::value);}
    static MC  Clone(Expression ll){C l = dynamic_cast<C>(ll);ffassert(l);return new CC(*l); }
    static void f(KN<K>   *x,MC l,A pp,Stack stack)
    {
        ffassert(l);
        double tgv=ff_tgv;
        if (l->nargs[0]) tgv= GetAny<double>((*l->nargs[0])(stack));
        FESpaceT * pVh= (*pp)->newVh();
        KN_<K>  xx=*x;
        if (  pVh && AssembleVarForm<K,MatriceCreuse<K>,MeshT,FESpaceT,FESpaceT>(stack,pVh->Th,*pVh,*pVh,false,0,&xx,l->largs) )
            AssembleBC<K,MeshT,FESpaceT,FESpaceT>(stack,pVh->Th,*pVh,*pVh,false,0,&xx,0,l->largs,tgv);
    }
};
template<class R=double>
struct CGMatVirtPreco : CGMatVirt<int,R>
{
    int n;
    MatriceMorse<R> *A;
    CGMatVirtPreco(Stack stack,const OneOperator* pprecon,MatriceMorse<R> *HA);
     R * addmatmul(R *x,R *Ax)
    {

    }
};

template<class R,class MMESH, class FESpace1, class FESpace2>
void creationBlockOfMatrixToBilinearForm( const FESpace1 * PUh, const FESpace2 * PVh, const int &sym, const double &tgv, 
                             const list<C_F0> & largs, Stack stack, Matrice_Creuse<R> &A);

template<class R,class MMesh,class v_fes1,class v_fes2>
AnyType OpMatrixtoBilinearForm<R,MMesh,v_fes1,v_fes2>::Op::operator()(Stack stack)  const
{
  typedef typename  MMesh::Element Element;
  typedef typename  MMesh::Vertex Vertex;
  typedef typename  MMesh::RdHat RdHat;
  typedef typename  MMesh::Rd Rd;
    
  typedef typename  v_fes1::pfes pfes1;
  typedef typename  v_fes1::FESpace FESpace1;
  typedef typename  FESpace1::Mesh Mesh1;
  typedef typename  FESpace1::FElement FElement1;
    
  typedef typename  v_fes2::pfes pfes2;
  typedef typename  v_fes2::FESpace FESpace2;
  typedef typename  FESpace2::Mesh Mesh2;
  typedef typename  FESpace2::FElement FElement2;
  assert(b && b->nargs);// *GetAny<pfes * >
  pfes1  * pUh= GetAny<pfes1 *>((*b->euh)(stack));
  pfes2  * pVh= GetAny<pfes2 *>((*b->evh)(stack));
  const FESpace1 * PUh =  (FESpace1*) **pUh ;
  const FESpace2 * PVh =  (FESpace2*) **pVh ;
  bool A_is_square= (void*)PUh == (void*)PVh || (PUh->NbOfDF) == (PVh->NbOfDF) ;

  bool VF=isVF(b->largs);
  Data_Sparse_Solver ds;
  ds.factorize=0;
  ds.initmat=true;
  int np = OpCall_FormBilinear_np::n_name_param - NB_NAME_PARM_HMAT;
  SetEnd_Data_Sparse_Solver<R>(stack,ds, b->nargs,np);

  if (! A_is_square )
   {
     if(verbosity>3) cout << " -- the solver  is un set  on rectangular matrix  " << endl;
   }
   WhereStackOfPtr2Free(stack)=new StackOfPtr2Free(stack);// FH aout 2007

  Matrice_Creuse<R> & A( * GetAny<Matrice_Creuse<R>*>((*a)(stack)));
  // 
  if(init) A.init(); //
  if( ! PUh || ! PVh) return SetAny<Matrice_Creuse<R>  *>(&A);

  creationBlockOfMatrixToBilinearForm<R,MMesh,FESpace1,FESpace2>( PUh, PVh, ds.sym, ds.tgv, b->largs, stack, A);
  /*
  const FESpace1 & Uh =  *PUh ;
  const FESpace2 & Vh =  *PVh ;
  const MMesh* pTh = (is_same< Mesh1, Mesh2 >::value) ? (MMesh*)&PUh->Th : 0;
  const MMesh &Th= *pTh ;    // integration Th
  bool same=isSameMesh(b->largs,&Uh.Th,&Vh.Th,stack);
  if ( same)
   {
     if ( A.Uh != Uh  || A.Vh != Vh )
      { // reconstruct all the matrix
        A.A=0; // to delete  old  matrix ADD FH 16112005
        A.Uh=Uh;
        A.Vh=Vh;
        if (ds.sym ){
          A.A.master( new  MatriceMorse<R>(Vh.NbOfDF,Vh.NbOfDF,ds.sym) );
	        ffassert( (void*)&Uh == (void*)&Vh);
        }
	      else
	        A.A.master( new  MatriceMorse<R>(Vh.NbOfDF,Uh.NbOfDF,Vh.NbOfDF*2,0) ); // lines corresponding to test functions
          // reset the solver ...
      }
     *A.A=R(); // reset value of the matrix

     if ( AssembleVarForm<R,MatriceCreuse<R>,MMesh,FESpace1,FESpace2 >( stack,Th,Uh,Vh,ds.sym>0,A.A,0,b->largs) )
       AssembleBC<R,MMesh,FESpace1,FESpace2>( stack,Th,Uh,Vh,ds.sym>0,A.A,0,0,b->largs,ds.tgv);
   }
  else
   { // add FH 17 06 2005  int on different meshes.
#ifdef V3__CODE
     MatriceMap<R>   AAA;
     MatriceMorse<R> *pMA =   new  MatriceMorse<R>(Vh.NbOfDF,Uh.NbOfDF,AAA.size(),ds.sym>0);
       bool bc=AssembleVarForm<R,MatriceMap<R>,MMesh,FESpace1,FESpace2>( stack,Th,Uh,Vh,ds.sym>0,&AAA,0,b->largs);
     pMA->addMap(1.,AAA);
#else
       MatriceMorse<R> *pMA =   new  MatriceMorse<R>(Vh.NbOfDF,Uh.NbOfDF,0,ds.sym);
       MatriceMap<R>  &  AAA = *pMA;
       bool bc=AssembleVarForm<R,MatriceMap<R>,MMesh,FESpace1,FESpace2>( stack,Th,Uh,Vh,ds.sym>0,&AAA,0,b->largs);

#endif
       A.A.master(pMA ) ;

       if (bc)
           AssembleBC<R>( stack,Th,Uh,Vh,ds.sym>0,A.A,0,0,b->largs,ds.tgv);

   }
   A.pHM()->half = ds.sym;
   if (A_is_square)
        SetSolver(stack,VF,*A.A,ds);
   */ 
  return SetAny<Matrice_Creuse<R>  *>(&A);

}
/*
list<C_F0>  creationLargsForCompositeFESpace( const list<C_F0> & largs, const int &NpUh, const int &NpVh, 
                                            const KN<int> &indexBlockUh, const KN<int> &indexBlockVh );

KNM< list<C_F0> > computeBlockLargs( const list<C_F0> & largs, const int &NpUh, const int &NpVh, 
                                    const KN<int> &indexBlockUh, const KN<int> &indexBlockVh );

void changeComponentFormCompositeFESpace( const KN<int> &localIndexInTheBlockUh, const KN<int> &localIndexInTheBlockVh, 
        KNM< list<C_F0> > & block_largs );
*/
template<class R,class MMesh, class FESpace1, class FESpace2>
void creationBlockOfMatrixToBilinearForm( const FESpace1 * PUh, const FESpace2 * PVh, const int &sym, const double &tgv, 
                             const list<C_F0> & largs, Stack stack, Matrice_Creuse<R> &A){
  typedef typename  FESpace1::Mesh Mesh1;
  typedef typename  FESpace2::Mesh Mesh2;

  // this lines must be defined outside this function
  // if(init) A.init();
  // if( ! PUh || ! PVh) return SetAny<Matrice_Creuse<R>  *>(&A);
  const FESpace1 & Uh =  *PUh ;
  const FESpace2 & Vh =  *PVh ;
  const MMesh* pTh = (is_same< Mesh1, Mesh2 >::value) ? (MMesh*)&PUh->Th : 0;
  const MMesh &Th= *pTh ;    // integration Th
  bool same=isSameMesh( largs, &Uh.Th, &Vh.Th, stack);
  if ( same)
   {
     if ( A.Uh != Uh  || A.Vh != Vh )
      { // reconstruct all the matrix
        A.A=0; // to delete  old  matrix ADD FH 16112005
        A.Uh=Uh;
        A.Vh=Vh;
        if (sym ){
          A.A.master( new  MatriceMorse<R>(Vh.NbOfDF,Vh.NbOfDF,sym) );
	        ffassert( (void*)&Uh == (void*)&Vh);
        }
	      else
	        A.A.master( new  MatriceMorse<R>(Vh.NbOfDF,Uh.NbOfDF,Vh.NbOfDF*2,0) ); // lines corresponding to test functions
          // reset the solver ...
      }
     *A.A=R(); // reset value of the matrix
     if ( AssembleVarForm<R,MatriceCreuse<R>,MMesh,FESpace1,FESpace2 >( stack,Th,Uh,Vh,sym>0,A.A,0,largs) )
       AssembleBC<R,MMesh,FESpace1,FESpace2>( stack,Th,Uh,Vh,sym>0,A.A,0,0,largs,tgv);
   }
  else 
   { // add FH 17 06 2005  int on different meshes.
#ifdef V3__CODE
     MatriceMap<R>   AAA;
     MatriceMorse<R> *pMA =   new  MatriceMorse<R>(Vh.NbOfDF,Uh.NbOfDF,AAA.size(),sym>0);
       bool bc=AssembleVarForm<R,MatriceMap<R>,MMesh,FESpace1,FESpace2>( stack,Th,Uh,Vh,sym>0,&AAA,0,largs);
     pMA->addMap(1.,AAA);
#else
       MatriceMorse<R> *pMA =   new  MatriceMorse<R>(Vh.NbOfDF,Uh.NbOfDF,0,sym);
       MatriceMap<R>  &  AAA = *pMA;
       bool bc=AssembleVarForm<R,MatriceMap<R>,MMesh,FESpace1,FESpace2>( stack,Th,Uh,Vh,sym>0,&AAA,0,largs);

#endif
       A.A.master(pMA ) ;

       if (bc)
           AssembleBC<R>( stack,Th,Uh,Vh,sym>0,A.A,0,0,largs,tgv);

   }
}
#include "compositeFESpace.hpp"
/*
template<class R,class MMesh, class FESpace1, class FESpace2>
void varfToCompositeBlockLinearSystem(bool initmat, bool initx, const FESpace1 * PUh, const FESpace2 * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<R> *B, KN_<R> *X, MatriceCreuse<R> &A);


// Mesh - Mesh
template void varfToCompositeBlockLinearSystem< double, Mesh, FESpace, FESpace>
                              (bool initmat, bool initx, const FESpace * PUh, const FESpace * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<double> *B, KN_<double> *X, MatriceCreuse<double> &A);

template void varfToCompositeBlockLinearSystem< Complex, Mesh, FESpace, FESpace>
                              (bool initmat, bool initx, const FESpace * PUh, const FESpace * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<Complex> *B, KN_<Complex> *X, MatriceCreuse<Complex> &A);
// MeshL - MeshL
template void varfToCompositeBlockLinearSystem< double, MeshL, FESpaceL, FESpaceL>
                              (bool initmat, bool initx, const FESpaceL * PUh, const FESpaceL * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<double> *B, KN_<double> *X, MatriceCreuse<double> &A);

template void varfToCompositeBlockLinearSystem< Complex, MeshL, FESpaceL, FESpaceL>
                              (bool initmat, bool initx, const FESpaceL * PUh, const FESpaceL * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<Complex> *B, KN_<Complex> *X, MatriceCreuse<Complex> &A);

// Mesh - MeshL
template void varfToCompositeBlockLinearSystem< double, MeshL, FESpace, FESpaceL>
                              (bool initmat, bool initx, const FESpace * PUh, const FESpaceL * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<double> *B, KN_<double> *X, MatriceCreuse<double> &A);

template void varfToCompositeBlockLinearSystem< Complex, MeshL, FESpace, FESpaceL>
                              (bool initmat, bool initx, const FESpace * PUh, const FESpaceL * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<Complex> *B, KN_<Complex> *X, MatriceCreuse<Complex> &A);

// MeshL - Mesh
template void varfToCompositeBlockLinearSystem< double, MeshL, FESpaceL, FESpace>
                              (bool initmat, bool initx, const FESpaceL * PUh, const FESpace * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<double> *B, KN_<double> *X, MatriceCreuse<double> &A);

template void varfToCompositeBlockLinearSystem< Complex, MeshL, FESpaceL, FESpace>
                              (bool initmat, bool initx, const FESpaceL * PUh, const FESpace * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<Complex> *B, KN_<Complex> *X, MatriceCreuse<Complex> &A);
                              */
/*
// MeshL - MeshL
template void varfToCompositeBlockLinearSystem<class complex,class MMesh, class FESpace1, class FESpace2>(bool initmat, bool initx, const FESpace1 * PUh, const FESpace2 * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<R> *B, KN_<R> *X, MatriceCreuse<R> &A);
*/


/*
template<class R,class MMesh, class FESpace1, class FESpace2>
void varfToCompositeBlockLinearSystem(bool initmat, bool initx, const FESpace1 * PUh, const FESpace2 * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<R> *B, KN_<R> *X, MatriceCreuse<R> &A)
                              {
  typedef typename  FESpace1::Mesh Mesh1;
  typedef typename  FESpace2::Mesh Mesh2;

  // this lines must be defined outside this function
  // if(init) A.init();
  // if( ! PUh || ! PVh) return SetAny<Matrice_Creuse<R>  *>(&A);
  const FESpace1 & Uh =  *PUh ;
  const FESpace2 & Vh =  *PVh ;
  const MMesh* pTh = (is_same< Mesh1, Mesh2 >::value) ? (MMesh*)&PUh->Th : 0;
  const MMesh &Th= *pTh ;    // integration Th
  bool same=isSameMesh( largs, &Uh.Th, &Vh.Th, stack);

  // PAC(e) :  on retourne le type MatriceCreuse Ici et non Matrice_Creuse<R> dans creationBlockOfMatrixToBilinearForm.
  // Attention pour la generalisation
  if(same){
    if  (AssembleVarForm<R,MatriceCreuse<R>,MMesh,FESpace1,FESpace2 >( stack,Th,Uh,Vh,sym, initmat ? &A:0 , B, largs))
    {
      if( B ){
        *B = - *B;
        // hach FH
        for (int i=0, n= B->N(); i< n; i++)
        if( abs((*B)[i]) < 1.e-60 ) (*B)[i]=0;
      }
      // AssembleBC<R,MMesh,FESpace1,FESpace2> ( stack,Th,Uh,Vh,sym, initmat ? &A:0 , B, initx ? X:0,  largs, tgv );   // TODO with problem
      AssembleBC<R> ( stack,Th,Uh,Vh,sym, initmat ? &A:0 , B, initx ? X:0,  largs, tgv );   // TODO with problem
    }
    else{
      if( B ) *B = - *B;
    }
  }else{
#ifdef V3__CODE
    MatriceMap<R>   AAA;
    cout << "V3__CODE=" << AAA.size() << endl;
    ffassert(0); // code a faire
    MatriceMorse<R> *pMA =   new  MatriceMorse<R>(Vh.NbOfDF,Uh.NbOfDF,AAA.size(),sym>0);
    bool bc=AssembleVarForm<R,MatriceMap<R>,MMesh,FESpace1,FESpace2>( stack,Th,Uh,Vh,sym>0,initmat ? &AAA:0,B,largs);
    pMA->addMap(1.,AAA);
#else
    MatriceMorse<R> *pMA =  dynamic_cast< HashMatrix<int,R>*>(&A);// new  MatriceMorse<R>(Vh.NbOfDF,Uh.NbOfDF,0,sym);
    MatriceMap<R>  &  AAA = *pMA;
    bool bc=AssembleVarForm<R,MatriceMap<R>,MMesh,FESpace1,FESpace2>( stack,Th,Uh,Vh,sym>0,initmat ? &AAA:0,B,largs);
#endif
    //cout << "AAA == matrice map=" << AAA << endl;
    //(*pMA) = AAA;
    //cout << "  *pMA == matrice map=" <<  *pMA << endl;
    //cout << "     A == matrice map=" <<  *dynamic_cast< HashMatrix<int,R>*>(&A) << endl;
    if (bc){
      if( B ){
        *B = - *B;
        // hach FH
        for (int i=0, n= B->N(); i< n; i++)
        if( abs((*B)[i]) < 1.e-60 ) (*B)[i]=0;
      }
      AssembleBC<R> ( stack,Th,Uh,Vh,sym, initmat ? &A:0 , B, initx ? X:0,  largs, tgv );   // TODO with problem
    }
    else{
      if( B ) *B = - *B;
    }
  }
}
*/

//bool SetGMRES();
//bool SetCG();
#ifdef HAVE_LIBUMFPACK
//bool SetUMFPACK();
#endif

namespace FreeFempp {

template<class R>
class TypeVarForm { public:
    aType tFB;
    aType tMat;
    aType tMat3;
    aType tFL;
  //aType tFL3;
    aType tTab;
    aType tMatX;
    aType tMatTX;
    aType tDotStar;
    aType tBC ;
  // aType tBC3 ;
TypeVarForm() :
  tFB( atype<const  FormBilinear *>() ),
 // tBemKFB( atype<const BemKFormBilinear *>() ),
  //tFB3( atype<const  FormBilinear<v_fes3> *>() ),
  tMat( atype<Matrice_Creuse<R>*>() ),
  //  tMat3( atype<Matrice_Creuse<R,v_fes3>*>() ),
  tFL( atype<const  FormLinear *>() ),
  //tFL3( atype<const  FormLinear<v_fes3> *>() ),
  tTab( atype<KN<R> *>() ),
  tMatX( atype<typename RNM_VirtualMatrix<R>::plusAx >() ),
  tMatTX( atype<typename RNM_VirtualMatrix<R>::plusAtx >() ),
  tDotStar(atype< DotStar_KN_<R> >() ),
  tBC( atype<const  BC_set  *>())
  {  }


 static TypeVarForm *Global;
};

}
#endif
