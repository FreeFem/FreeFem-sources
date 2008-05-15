#ifndef LGMESH3_HPP
#define LGMESH3_HPP
// 3d
#define AFAIRE(cmm) { cerr << "A Faire/ To Do  "<< cmm << " file " << __FILE__<< " line "<< __LINE__ << endl; InternalError(cmm);}
typedef FEbase<double,v_fes3> * pf3rbase ;
typedef FEbaseArray<double,v_fes3> * pf3rbasearray ;
typedef pair<pf3rbase,int> pf3r ;
typedef pair<pf3rbasearray,int> pf3rarray ;

typedef FEbase<Complex,v_fes3> * pf3cbase ;
typedef FEbaseArray<Complex,v_fes3> * pf3cbasearray ;
typedef pair<pf3cbase,int> pf3c ;
typedef pair<pf3cbasearray,int> pf3carray ;
// fin
bool isSameMesh(const list<C_F0> & largs,const void * Thu,const void * Thv,Stack stack) ; // true => VF type of Matrix   
  //bool isSameMesh(const list<C_F0> & largs,const Mesh * Thu,const Mesh * Thv,Stack stack)  ;


inline C_F0 CCastToR(const C_F0 & f){ return C_F0(atype<double>()->CastTo(f),atype<double>());}
inline bool BCastToR(const C_F0 & f){ return atype<double>()->CastingFrom(f.left());}



inline C_F0 CCastToC(const C_F0 & f){ return C_F0(atype<Complex>()->CastTo(f),atype<Complex>());}
inline bool BCastToC(const C_F0 & f){ return atype<Complex>()->CastingFrom(f.left());}

template<class Result>
inline Expression CastTo(const C_F0 & f) { return atype<Result>()->CastTo(f);}

template<class Result>
inline bool BCastTo(const C_F0 & f) { return atype<Result>()->CastingFrom(f.left());}

inline void Check(bool  v,const char * mess)
{
  if (!v) { cerr << " Error " << mess ;
  throw(ErrorExec(mess,1));
  }
}           

 
#endif
