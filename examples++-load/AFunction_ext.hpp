// In order to use functions with the stack and 2 parameters , and more  new classes (OneOperator2s_, OneOperator3s_, etc.) must
// In order to use functions with more than 3 parameters, new classes (OneOperator4_, OneOperator5_, etc.) must
// be defined. See example code in include/AFunction.hpp
// Two classes must be defined (here we show an example for a function accepting 4 arguments):
//
//  class OneOperator4_ 
//  class E_F_F0F0F0F0_
//
// Note: in file includeAFunction.hpp, the class "OneOperator" (around line 400) mut be modified.
// ************************************************
// Add F. Hecht  oct 2009
// ****  2 paramters with the stack
//  class OneOperator2s_
//  class E_F_F0F0s_                                                                                                                                                      

template<class R,class A0,class A1, class E=E_F0>   // extend (4th arg.)
class E_F_F0F0s_ :public  E { public:                               // extend 
    typedef  R (*func)(Stack s,const  A0 &,const  A1 & ) ; // extend (statck +2th arg.)
  func f;
  Expression a0,a1;          // extend
  E_F_F0F0s_(func ff,
	     Expression aa0,
	     Expression aa1)
    : f(ff), a0(aa0), a1(aa1) {}  // extend (2th arg.)
  AnyType operator()(Stack s)  const 
  {return SetAny<R>( f( s,
			GetAny<A0>((*a0)(s)),
			GetAny<A1>((*a1)(s))  ) );}   // extend (2th arg.)
  virtual size_t nbitem() const {return a1->nbitem(); } // modif ??? 
  bool MeshIndependent() const 
  {return a0->MeshIndependent() && a1->MeshIndependent() ;} // extend (2th arg.)
  
};

template<class R,class A=R,class B=A,class C=B, class D=C ,class CODE=E_F_F0F0s_<R,A,B,E_F0> >    // extend (4th arg.)
class  OneOperator2s_ : public OneOperator {     // 
  aType r; //  return type 
  typedef typename  CODE::func  func;
  func f;
public: 
  E_F0 * code(const basicAC_F0 & args) const 
  { return  new CODE(f,
		     t[0]->CastTo(args[0]),
		     t[1]->CastTo(args[1]));}     // extend
  OneOperator2s_(func  ff):                        // 3->4
    OneOperator(map_type[typeid(R).name()],
		map_type[typeid(A).name()],
		map_type[typeid(B).name()]),      // extens
    f(ff){}
};



// ****  2 paramters with the stack
//  class OneOperator2s_
//  class E_F_F0F0s_                                                                                                                                                      

template<class R,class A0,class A1,class A2, class E=E_F0>   // extend (4th arg.)
class E_F_F0F0F0s_ :public  E { public:                               // extend 
    typedef  R (*func)(Stack s,const  A0 &,const  A1 &,const A2 & ) ; // extend (statck +2th arg.)
  func f;
  Expression a0,a1,a2;          // extend
  E_F_F0F0F0s_(func ff,
	     Expression aa0,
	     Expression aa1,
	     Expression aa2)
    : f(ff), a0(aa0), a1(aa1), a2(aa2) {}  // extend (2th arg.)
  AnyType operator()(Stack s)  const 
  {return SetAny<R>( f( s,
			GetAny<A0>((*a0)(s)),
			GetAny<A0>((*a1)(s)),
			GetAny<A1>((*a2)(s))  ) );}   // extend (3th arg.)
  virtual size_t nbitem() const {return a2->nbitem(); } // modif ??? 
  bool MeshIndependent() const 
  {return a0->MeshIndependent() && a1->MeshIndependent() && a2->MeshIndependent() ;} // extend (2th arg.)
  
};

template<class R,class A=R,class B=A,class C=B, class D=C ,class CODE=E_F_F0F0F0s_<R,A,B,C,E_F0> >    // extend (3th arg.)
class  OneOperator3s_ : public OneOperator {     // 
  aType r; //  return type 
  typedef typename  CODE::func  func;
  func f;
public: 
  E_F0 * code(const basicAC_F0 & args) const 
  { return  new CODE(f,
		     t[0]->CastTo(args[0]),
		     t[0]->CastTo(args[1]),
		     t[1]->CastTo(args[2]));}     // extend
  OneOperator3s_(func  ff):                        // 2->
    OneOperator(map_type[typeid(R).name()],
		map_type[typeid(A).name()],
		map_type[typeid(B).name()],
		map_type[typeid(C).name()]),      // extend
    f(ff){}
};
// ***********************************************

// ***********************************************
// **** 4 parameters
// ***********************
template<class R,class A0,class A1,class A2, class A3, class E=E_F0>   // extend (4th arg.)
class E_F_F0F0F0F0_ :public  E { public:                               // extend 
   typedef  R (*func)(const  A0 &,const  A1 & , const A2 &, const A3 & ) ; // extend (4th arg.)
  func f;
  Expression a0,a1,a2,a3;          // extend
  E_F_F0F0F0F0_(func ff,
		Expression aa0,
		Expression aa1,
		Expression aa2,
		Expression aa3)   // extend 
    : f(ff), a0(aa0), a1(aa1), a2(aa2), a3(aa3) {}  // extend (4th arg.)
  AnyType operator()(Stack s)  const 
    {return SetAny<R>( f( GetAny<A0>((*a0)(s)),
			  GetAny<A1>((*a1)(s)),
			  GetAny<A2>((*a2)(s)),
			  GetAny<A3>((*a3)(s))  ) );}   // extend (4th arg.)
  virtual size_t nbitem() const {return a3->nbitem(); } // modif
      bool MeshIndependent() const 
      {return a0->MeshIndependent() && a1->MeshIndependent()&& a2->MeshIndependent()&& a3->MeshIndependent();} // extend (4th arg.)

};

template<class R,class A=R,class B=A,class C=B, class D=C ,class CODE=E_F_F0F0F0F0_<R,A,B,C,D,E_F0> >    // extend (4th arg.)
class  OneOperator4_ : public OneOperator {     // 3->4
  aType r; //  return type 
    typedef typename  CODE::func  func;
  func f;
public: 
  E_F0 * code(const basicAC_F0 & args) const 
  { return  new CODE(f,
		     t[0]->CastTo(args[0]),
		     t[1]->CastTo(args[1]),
		     t[2]->CastTo(args[2]),
		     t[3]->CastTo(args[3]));}     // extend
  OneOperator4_(func  ff):                        // 3->4
    OneOperator(map_type[typeid(R).name()],
		map_type[typeid(A).name()],
		map_type[typeid(B).name()],
		map_type[typeid(C).name()],
		map_type[typeid(D).name()]),      // extens
    f(ff){}
};

// ***********************************************
// **** 5 parameters
// ***********************
//
//  NOTE: add the following line in AFunction.hpp
//
//    OneOperator(aType rr,aType  a,aType  b,aType c,aType d,aType e) 
//      : r(rr),ArrayOfaType(a,b,c,d,e,false),next(0),pref(0) {throwassert(rr && a && b && c && d);} 

template<class R,class A0,class A1,class A2, class A3, class A4, class E=E_F0>   // extend AX
class E_F_F0F0F0F0F0_ :public  E { public:                               // extend 
   typedef  R (*func)(const  A0 &,const  A1 & , const A2 &, const A3 &, const A4 & ) ; // extend AX
  func f;
  Expression a0,a1,a2,a3,a4;          // extend aX
  E_F_F0F0F0F0F0_(func ff,            // extend F0
		  Expression aa0,
		  Expression aa1,
		  Expression aa2,
		  Expression aa3,
		  Expression aa4)   // extend 
    : f(ff), a0(aa0), a1(aa1), a2(aa2), a3(aa3), a4(aa4) {}  // extend aX
  AnyType operator()(Stack s)  const 
    {return SetAny<R>( f( GetAny<A0>((*a0)(s)),
			  GetAny<A1>((*a1)(s)),
			  GetAny<A2>((*a2)(s)),
			  GetAny<A3>((*a3)(s)),
			  GetAny<A4>((*a4)(s))  ) );}   // extend aX
    virtual size_t nbitem() const {return a4->nbitem(); } 
      bool MeshIndependent() const 
      {return a0->MeshIndependent() && a1->MeshIndependent()&& a2->MeshIndependent()
	 && a3->MeshIndependent()&& a4->MeshIndependent();} // extend aX

};

template<class R,class A=R,class B=A,class C=B, class D=C ,class E=D ,class CODE=E_F_F0F0F0F0F0_<R,A,B,C,D,E,E_F0> >    // extend  
class  OneOperator5_ : public OneOperator {     // 3->4
  aType r; //  return type 
    typedef typename  CODE::func  func;
  func f;
public: 
  E_F0 * code(const basicAC_F0 & args) const 
  { return  new CODE(f,
		     t[0]->CastTo(args[0]),
		     t[1]->CastTo(args[1]),
		     t[2]->CastTo(args[2]),
		     t[3]->CastTo(args[3]),
		     t[4]->CastTo(args[4]));}     // extend
  OneOperator5_(func  ff):                        // 3->4
    OneOperator(map_type[typeid(R).name()],
		map_type[typeid(A).name()],
		map_type[typeid(B).name()],
		map_type[typeid(C).name()],
		map_type[typeid(D).name()],
		map_type[typeid(E).name()]),      // extend
    f(ff){}
};

// ***********************************************
// **** 6 parameters
// ***********************
//
//  NOTE: add the following line in AFunction.hpp
//    OneOperator(aType rr,aType  a,aType  b,aType c,aType d,aType e,aType f) 
//       : r(rr),ArrayOfaType(a,b,c,d,e,f,false),next(0),pref(0) {throwassert(rr && a && b && c && d && f);} 

template<class R,class A0,class A1,class A2, class A3, class A4, class A5, class E=E_F0>   // extend AX
class E_F_F0F0F0F0F0F0_ :public  E { public:                               // extend 
   typedef  R (*func)(const  A0 &,const  A1 & , const A2 &, const A3 &, const A4 &, const A5 & ) ; // extend AX
  func f;
  Expression a0,a1,a2,a3,a4,a5;          // extend aX
  E_F_F0F0F0F0F0F0_(func ff,            // extend F0
		    Expression aa0,
		    Expression aa1,
		    Expression aa2,
		    Expression aa3,
		    Expression aa4,
		    Expression aa5)   // extend 
    : f(ff), a0(aa0), a1(aa1), a2(aa2), a3(aa3), a4(aa4), a5(aa5) {}  // extend aX
  AnyType operator()(Stack s)  const 
    {return SetAny<R>( f( GetAny<A0>((*a0)(s)),
			  GetAny<A1>((*a1)(s)),
			  GetAny<A2>((*a2)(s)),
			  GetAny<A3>((*a3)(s)),
			  GetAny<A4>((*a4)(s)),
			  GetAny<A5>((*a5)(s)) ) );}   // extend aX
    virtual size_t nbitem() const {return a5->nbitem(); } 
      bool MeshIndependent() const 
      {return a0->MeshIndependent() && a1->MeshIndependent()&& a2->MeshIndependent()
	 && a3->MeshIndependent()&& a4->MeshIndependent()&& a5->MeshIndependent();} // extend aX

};

template<class R,class A=R,class B=A,class C=B, class D=C ,class E=D ,class F=E ,class CODE=E_F_F0F0F0F0F0F0_<R,A,B,C,D,E,F,E_F0> >    // extend  
class  OneOperator6_ : public OneOperator {     // 3->4
  aType r; //  return type 
    typedef typename  CODE::func  func;
  func f;
public: 
  E_F0 * code(const basicAC_F0 & args) const 
  { return  new CODE(f,
		     t[0]->CastTo(args[0]),
		     t[1]->CastTo(args[1]),
		     t[2]->CastTo(args[2]),
		     t[3]->CastTo(args[3]),
		     t[4]->CastTo(args[4]),
		     t[5]->CastTo(args[5]));}     // extend
  OneOperator6_(func  ff):                        // 3->4
    OneOperator(map_type[typeid(R).name()],
		map_type[typeid(A).name()],
		map_type[typeid(B).name()],
		map_type[typeid(C).name()],
		map_type[typeid(D).name()],
		map_type[typeid(E).name()],
		map_type[typeid(F).name()]),      // extend
    f(ff){}
};


// ***********************************************
// **** 7 parameters
// ***********************
//
//  NOTE: add the following line in AFunction.hpp
//    OneOperator(aType rr,aType  a,aType  b,aType c,aType d,aType e,aType f) 
//       : r(rr),ArrayOfaType(a,b,c,d,e,f,false),next(0),pref(0) {throwassert(rr && a && b && c && d && f);} 

template<class R,class A0,class A1,class A2, class A3, class A4, class A5, class A6, class E=E_F0>   // extend AX
class E_F_F0F0F0F0F0F0F0_ :public  E { public:                               // extend 
   typedef  R (*func)(const  A0 &,const  A1 & , const A2 &, const A3 &, const A4 &, const A5 &, const A6 & ) ; // extend AX
  func f;
  Expression a0,a1,a2,a3,a4,a5,a6;          // extend aX
  E_F_F0F0F0F0F0F0F0_(func ff,            // extend F0
		      Expression aa0,
		      Expression aa1,
		      Expression aa2,
		      Expression aa3,
		      Expression aa4,
		      Expression aa5,
		      Expression aa6)   // extend 
    : f(ff), a0(aa0), a1(aa1), a2(aa2), a3(aa3), a4(aa4), a5(aa5), a6(aa6) {}  // extend aX
  AnyType operator()(Stack s)  const 
    {return SetAny<R>( f( GetAny<A0>((*a0)(s)),
			  GetAny<A1>((*a1)(s)),
			  GetAny<A2>((*a2)(s)),
			  GetAny<A3>((*a3)(s)),
			  GetAny<A4>((*a4)(s)),
			  GetAny<A5>((*a5)(s)),
			  GetAny<A6>((*a6)(s)) ) );}   // extend aX
  virtual size_t nbitem() const {return a6->nbitem(); } // modif
      bool MeshIndependent() const 
      {return a0->MeshIndependent() && a1->MeshIndependent()&& a2->MeshIndependent()
	 && a3->MeshIndependent()&& a4->MeshIndependent()&& a5->MeshIndependent()&& a6->MeshIndependent();} // extend aX

};

template<class R,class A=R,class B=A,class C=B, class D=C ,class E=D ,class F=E ,class G=F ,class CODE=E_F_F0F0F0F0F0F0F0_<R,A,B,C,D,E,F,G,E_F0> >    // extend  
class  OneOperator7_ : public OneOperator {     // 3->4
  aType r; //  return type 
    typedef typename  CODE::func  func;
  func f;
public: 
  E_F0 * code(const basicAC_F0 & args) const 
  { return  new CODE(f,
		     t[0]->CastTo(args[0]),
		     t[1]->CastTo(args[1]),
		     t[2]->CastTo(args[2]),
		     t[3]->CastTo(args[3]),
		     t[4]->CastTo(args[4]),
		     t[5]->CastTo(args[5]),
		     t[6]->CastTo(args[6]));}     // extend
  OneOperator7_(func  ff):                        // 3->4
    OneOperator(map_type[typeid(R).name()],
		map_type[typeid(A).name()],
		map_type[typeid(B).name()],
		map_type[typeid(C).name()],
		map_type[typeid(D).name()],
		map_type[typeid(E).name()],
		map_type[typeid(F).name()],
		map_type[typeid(G).name()]),      // extend
    f(ff){}
};



// ***********************************************
// **** 8 parameters
// ***********************
//

template<class R,class A0,class A1,class A2, class A3, class A4, class A5, class A6, class A7, class E=E_F0>   // extend AX
class E_F_F0F0F0F0F0F0F0F0_ :public  E { public:                               // extend 
   typedef  R (*func)(const  A0 &,const  A1 & , const A2 &, const A3 &, const A4 &, const A5 &, const A6 &, const A7 & ) ; // extend AX
  func f;
  Expression a0,a1,a2,a3,a4,a5,a6,a7;          // extend aX
  E_F_F0F0F0F0F0F0F0F0_(func ff,            // extend F0
			Expression aa0,
			Expression aa1,
			Expression aa2,
			Expression aa3,
			Expression aa4,
			Expression aa5,
			Expression aa6,
			Expression aa7)   // extend 
    : f(ff), a0(aa0), a1(aa1), a2(aa2), a3(aa3), a4(aa4), a5(aa5), a6(aa6), a7(aa7) {}  // extend aX
  AnyType operator()(Stack s)  const 
    {return SetAny<R>( f( GetAny<A0>((*a0)(s)),
			  GetAny<A1>((*a1)(s)),
			  GetAny<A2>((*a2)(s)),
			  GetAny<A3>((*a3)(s)),
			  GetAny<A4>((*a4)(s)),
			  GetAny<A5>((*a5)(s)),
			  GetAny<A6>((*a6)(s)),
			  GetAny<A7>((*a7)(s)) ) );}   // extend aX
  virtual size_t nbitem() const {return a7->nbitem(); } // modif
      bool MeshIndependent() const 
      {return a0->MeshIndependent() && a1->MeshIndependent()&& a2->MeshIndependent()
	 && a3->MeshIndependent()&& a4->MeshIndependent()&& a5->MeshIndependent()&& a6->MeshIndependent()&& a7->MeshIndependent();} // extend aX

};

template<class R,class A=R,class B=A,class C=B, class D=C ,class E=D ,class F=E ,class G=F ,class H=G , class CODE=E_F_F0F0F0F0F0F0F0F0_<R,A,B,C,D,E,F,G,H,E_F0> >    // extend  
class  OneOperator8_ : public OneOperator {     // 3->4
  aType r; //  return type 
    typedef typename  CODE::func  func;
  func f;
public: 
  E_F0 * code(const basicAC_F0 & args) const 
  { return  new CODE(f,
		     t[0]->CastTo(args[0]),
		     t[1]->CastTo(args[1]),
		     t[2]->CastTo(args[2]),
		     t[3]->CastTo(args[3]),
		     t[4]->CastTo(args[4]),
		     t[5]->CastTo(args[5]),
		     t[6]->CastTo(args[6]),
		     t[7]->CastTo(args[7]));}     // extend
  OneOperator8_(func  ff):                        // 3->4
    OneOperator(map_type[typeid(R).name()],
		map_type[typeid(A).name()],
		map_type[typeid(B).name()],
		map_type[typeid(C).name()],
		map_type[typeid(D).name()],
		map_type[typeid(E).name()],
		map_type[typeid(F).name()],
		map_type[typeid(G).name()],
		map_type[typeid(H).name()]),      // extend
    f(ff){}
};


// ***********************************************
// **** 9 parameters
// ***********************
//

/*
template<class R,class A0,class A1,class A2, class A3, class A4, class A5, class A6, class A7, class A8, class E=E_F0>   // extend AX
class E_F_F0F0F0F0F0F0F0F0F0_ :public  E { public:                               // extend 
   typedef  R (*func)(const  A0 &,const  A1 & , const A2 &, const A3 &, const A4 &, const A5 &, const A6 &, const A7 &, const A8 & ) ; // extend AX
  func f;
  Expression a0,a1,a2,a3,a4,a5,a6,a7,a8;          // extend aX
  E_F_F0F0F0F0F0F0F0F0F0_(func ff,            // extend F0
			  Expression aa0,
			  Expression aa1,
			  Expression aa2,
			  Expression aa3,
			  Expression aa4,
			  Expression aa5,
			  Expression aa6,
			  Expression aa7,
			  Expression aa8)   // extend 
    : f(ff), a0(aa0), a1(aa1), a2(aa2), a3(aa3), a4(aa4), a5(aa5), a6(aa6), a7(aa7), a8(aa8) {}  // extend aX
  AnyType operator()(Stack s)  const 
    {return SetAny<R>( f( GetAny<A0>((*a0)(s)),
			  GetAny<A1>((*a1)(s)),
			  GetAny<A2>((*a2)(s)),
			  GetAny<A3>((*a3)(s)),
			  GetAny<A4>((*a4)(s)),
			  GetAny<A5>((*a5)(s)),
			  GetAny<A6>((*a6)(s)),
			  GetAny<A7>((*a7)(s)),
			  GetAny<A8>((*a8)(s)) ) );}   // extend aX
  virtual size_t nbitem() const {return a8->nbitem(); } // modif
      bool MeshIndependent() const 
      {return a0->MeshIndependent() && a1->MeshIndependent()&& a2->MeshIndependent()
	 && a3->MeshIndependent()&& a4->MeshIndependent()&& a5->MeshIndependent()&& a6->MeshIndependent()&& a7->MeshIndependent()&& a8->MeshIndependent();} // extend aX

};

template<class R,class A=R,class B=A,class C=B, class D=C ,class E=D ,class F=E ,class G=F ,class H=G ,class I=H , class CODE=E_F_F0F0F0F0F0F0F0F0F0_<R,A,B,C,D,E,F,G,H,I,E_F0> >    // extend  
class  OneOperator9_ : public OneOperator {     // 3->4
  aType r; //  return type 
    typedef typename  CODE::func  func;
  func f;
public: 
  E_F0 * code(const basicAC_F0 & args) const 
  { return  new CODE(f,
		     t[0]->CastTo(args[0]),
		     t[1]->CastTo(args[1]),
		     t[2]->CastTo(args[2]),
		     t[3]->CastTo(args[3]),
		     t[4]->CastTo(args[4]),
		     t[5]->CastTo(args[5]),
		     t[6]->CastTo(args[6]),
		     t[7]->CastTo(args[7]),
		     t[8]->CastTo(args[8]));}     // extend
  OneOperator9_(func  ff):                        // 3->4
    OneOperator(map_type[typeid(R).name()],
		map_type[typeid(A).name()],
		map_type[typeid(B).name()],
		map_type[typeid(C).name()],
		map_type[typeid(D).name()],
		map_type[typeid(E).name()],
		map_type[typeid(F).name()],
		map_type[typeid(G).name()],
		map_type[typeid(H).name()],
		map_type[typeid(I).name()]),      // extend
    f(ff){}
};
*/

/*

//
// ***********************************************
// **** 10 parameters
// ***********************
//

template<class R,class A0,class A1,class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9, class E=E_F0>   // extend AX
class E_F_F0F0F0F0F0F0F0F0F0F0_ :public  E { public:                               // extend 
   typedef  R (*func)(const  A0 &,const  A1 & , const A2 &, const A3 &, const A4 &, const A5 &, const A6 &, const A7 &, const A8 &, const A9 & ) ; // extend AX
  func f;
  Expression a0,a1,a2,a3,a4,a5,a6,a7,a8,a9;          // extend aX
  E_F_F0F0F0F0F0F0F0F0F0F0_(func ff,            // extend F0
			    Expression aa0,
			    Expression aa1,
			    Expression aa2,
			    Expression aa3,
			    Expression aa4,
			    Expression aa5,
			    Expression aa6,
			    Expression aa7,
			    Expression aa8,
			    Expression aa9)   // extend 
    : f(ff), a0(aa0), a1(aa1), a2(aa2), a3(aa3), a4(aa4), a5(aa5), a6(aa6), a7(aa7), a8(aa8), a9(aa9) {}  // extend aX
  AnyType operator()(Stack s)  const 
    {return SetAny<R>( f( GetAny<A0>((*a0)(s)),
			  GetAny<A1>((*a1)(s)),
			  GetAny<A2>((*a2)(s)),
			  GetAny<A3>((*a3)(s)),
			  GetAny<A4>((*a4)(s)),
			  GetAny<A5>((*a5)(s)),
			  GetAny<A6>((*a6)(s)),
			  GetAny<A7>((*a7)(s)),
			  GetAny<A8>((*a8)(s)),
			  GetAny<A9>((*a9)(s)) ) );}   // extend aX
  virtual size_t nbitem() const {return a9->nbitem(); } // modif 
      bool MeshIndependent() const 
      {return a0->MeshIndependent() && a1->MeshIndependent()&& a2->MeshIndependent()
	 && a3->MeshIndependent()&& a4->MeshIndependent()&& a5->MeshIndependent()&& a6->MeshIndependent()
	 && a7->MeshIndependent()&& a8->MeshIndependent()&& a9->MeshIndependent();} // extend aX

};

template<class R,class A=R,class B=A,class C=B, class D=C ,class E=D ,class F=E ,class G=F ,class H=G ,class I=H , class J=I , class CODE=E_F_F0F0F0F0F0F0F0F0F0F0_<R,A,B,C,D,E,F,G,H,I,J,E_F0> >    // extend  
class  OneOperator10_ : public OneOperator {     // 3->4
  aType r; //  return type 
    typedef typename  CODE::func  func;
  func f;
public: 
  E_F0 * code(const basicAC_F0 & args) const 
  { return  new CODE(f,
		     t[0]->CastTo(args[0]),
		     t[1]->CastTo(args[1]),
		     t[2]->CastTo(args[2]),
		     t[3]->CastTo(args[3]),
		     t[4]->CastTo(args[4]),
		     t[5]->CastTo(args[5]),
		     t[6]->CastTo(args[6]),
		     t[7]->CastTo(args[7]),
		     t[8]->CastTo(args[8]),
		     t[9]->CastTo(args[9]));}     // extend
  OneOperator10_(func  ff):                        // 3->4
    OneOperator(map_type[typeid(R).name()],
		map_type[typeid(A).name()],
		map_type[typeid(B).name()],
		map_type[typeid(C).name()],
		map_type[typeid(D).name()],
		map_type[typeid(E).name()],
		map_type[typeid(F).name()],
		map_type[typeid(G).name()],
		map_type[typeid(H).name()],
		map_type[typeid(I).name()],
		map_type[typeid(J).name()]),      // extend
    f(ff){}
};

*/
