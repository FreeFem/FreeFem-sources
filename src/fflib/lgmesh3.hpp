#ifndef LGMESH3_HPP
#define LGMESH3_HPP
// 3d real (2d equivalent at [[file:problem.hpp::pferbase]])
typedef FEbase<double,v_fes3> * pf3rbase ;
typedef FEbaseArray<double,v_fes3> * pf3rbasearray ;
typedef pair<pf3rbase,int> pf3r ;
typedef pair<pf3rbasearray,int> pf3rarray ;

// 3d complex (2d equivalent at [[file:problem.hpp::pfecbase]])
typedef FEbase<Complex,v_fes3> * pf3cbase ;
typedef FEbaseArray<Complex,v_fes3> * pf3cbasearray ;
typedef pair<pf3cbase,int> pf3c ;
typedef pair<pf3cbasearray,int> pf3carray ;
// fin

// Surf real (2d equivalent at [[file:problem.hpp::pferbase]])
typedef FEbase<double,v_fesS> * pfSrbase ;
typedef FEbaseArray<double,v_fesS> * pfSrbasearray ;
typedef pair<pfSrbase,int> pfSr ;
typedef pair<pfSrbasearray,int> pfSrarray ;

// Surf complex (2d equivalent at [[file:problem.hpp::pfecbase]])
typedef FEbase<Complex,v_fesS> * pfScbase ;
typedef FEbaseArray<Complex,v_fesS> * pfScbasearray ;
typedef pair<pfScbase,int> pfSc ;
typedef pair<pfScbasearray,int> pfScarray ;

// Curve real (2d equivalent at [[file:problem.hpp::pferbase]])
typedef FEbase<double,v_fesL> * pfLrbase ;
typedef FEbaseArray<double,v_fesL> * pfLrbasearray ;
typedef pair<pfLrbase,int> pfLr ;
typedef pair<pfLrbasearray,int> pfLrarray ;

// Curve complex (2d equivalent at [[file:problem.hpp::pfecbase]])
typedef FEbase<Complex,v_fesL> * pfLcbase ;
typedef FEbaseArray<Complex,v_fesL> * pfLcbasearray ;
typedef pair<pfLcbase,int> pfLc ;
typedef pair<pfLcbasearray,int> pfLcarray ;       

// fin

// 3d distributed volume
typedef FEbase<double, v_dfes3>* pdf3rbase;   typedef FEbaseArray<double, v_dfes3>* pdf3rbasearray;
typedef pair<pdf3rbase,int> pdf3r;            typedef pair<pdf3rbasearray,int> pdf3rarray;
typedef FEbase<Complex, v_dfes3>* pdf3cbase;  typedef FEbaseArray<Complex, v_dfes3>* pdf3cbasearray;
typedef pair<pdf3cbase,int> pdf3c;            typedef pair<pdf3cbasearray,int> pdf3carray;
// 3d distributed surface
typedef FEbase<double, v_dfesS>* pdfSrbase;   typedef FEbaseArray<double, v_dfesS>* pdfSrbasearray;
typedef pair<pdfSrbase,int> pdfSr;            typedef pair<pdfSrbasearray,int> pdfSrarray;
typedef FEbase<Complex, v_dfesS>* pdfScbase;  typedef FEbaseArray<Complex, v_dfesS>* pdfScbasearray;
typedef pair<pdfScbase,int> pdfSc;            typedef pair<pdfScbasearray,int> pdfScarray;
// 3d distributed curve
typedef FEbase<double, v_dfesL>* pdfLrbase;   typedef FEbaseArray<double, v_dfesL>* pdfLrbasearray;
typedef pair<pdfLrbase,int> pdfLr;            typedef pair<pdfLrbasearray,int> pdfLrarray;
typedef FEbase<Complex, v_dfesL>* pdfLcbase;  typedef FEbaseArray<Complex, v_dfesL>* pdfLcbasearray;
typedef pair<pdfLcbase,int> pdfLc;            typedef pair<pdfLcbasearray,int> pdfLcarray;
// fin



bool isSameMesh(const list<C_F0> & largs,const void * Thu,const void * Thv,Stack stack) ; // true => VF type of Matrix   
  //bool isSameMesh(const list<C_F0> & largs,const Mesh * Thu,const Mesh * Thv,Stack stack)  ;


inline C_F0 CCastToR(const C_F0 & f){ return C_F0(atype<double>()->CastTo(f),atype<double>());}
inline bool BCastToR(const C_F0 & f){ return atype<double>()->CastingFrom(f.left());}



inline C_F0 CCastToC(const C_F0 & f){ return C_F0(atype<Complex>()->CastTo(f),atype<Complex>());}
inline bool BCastToC(const C_F0 & f){ return atype<Complex>()->CastingFrom(f.left());}

template<class Result>
inline Expression CastTo(const C_F0 & f) { return atype<Result>()->CastTo(f);}

// <<BCastTo>>
template<class Result>
inline bool BCastTo(const C_F0 & f) { return atype<Result>()->CastingFrom(f.left());}

inline void Check(bool  v,const char * mess)
{
  if (!v) { cerr << " Error " << mess ;
  throw(ErrorExec(mess,1));
  }
}           

 
#endif
